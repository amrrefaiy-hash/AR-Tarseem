# ==========================================================================
# AR-Tarseem — QGIS Plugin
# Version 1.1.0 — November 2026
# ==========================================================================
# Copyright (C) 2026 Ahmed M. Raslan, PhD & Amr R. Refaiy, PhD
#                    Water Engineering Experts
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <https://www.gnu.org/licenses/>.
# A copy is distributed with this plugin in the file LICENSE.
# ==========================================================================
"""
================================================================================
AR-TARSEEM  ─  QGIS HYDROLOGICAL DELINEATION PLUGIN CORE v20
================================================================================
Changes vs. v17:
    • Reorganised hydrograph outputs. The flat ``Hydrographs/``
      directory is replaced by a structured ``Modeling_Outputs/``
      tree with one folder per element type, each split into raw
      CSVs and interactive HTML visualisations:

          <OUTPUT_DIR>/Modeling_Outputs/
              ├── Catchments/   ├── CSVs/   W{pid}_{rp}yr.csv
              │                 └── HTMLs/  W{pid}_{rp}yr.html
              ├── SubBasins/    ├── CSVs/   W{pid}-Sub{sid}_{rp}yr.csv
              │                 └── HTMLs/  W{pid}-Sub{sid}_{rp}yr.html
              ├── Reaches/      ├── CSVs/   Reach_W{pid}-R{sid}_{rp}yr.csv
              │                 └── HTMLs/  Reach_W{pid}-R{sid}_{rp}yr.html
              └── Junctions/    ├── CSVs/   W{pid}-J{sid}_{rp}yr.csv
                                └── HTMLs/  W{pid}-J{sid}_{rp}yr.html

      Each HTML is a self-contained interactive visualisation built
      with Plotly (loaded from CDN). They show the hydrograph time
      series with annotated peak Q, peak time, and total runoff
      volume. The four element types each get a tailored layout:

          Catchments / SubBasins   Three stacked panels:
              • Rainfall hyetograph (P / Loss / Excess, mm/min)
              • Cumulative depth   (P / Loss / Excess, mm)
              • Discharge          (Q, m³/s)
          Reaches                  Inflow vs. routed outflow overlay
                                   with peak attenuation %.
          Junctions                Combined-flow line plot with peak
                                   annotation.

      All HTMLs are responsive (zoom, pan, hover, save-as-PNG via
      Plotly's modebar). They fall back gracefully if the Plotly
      CDN is unreachable: the page still renders the title and
      summary card, with a notice that the chart could not load.

Earlier versions:
    v17 – Sub-basin hydrographs, reach routing (Muskingum/LAG),
          Junctions point layer with combined hydrographs at every
          sub-basin outlet, topological accumulation.
    v16 – Automatic QML styling: CONFIG['STYLES_DIR'] applies
          built-in QML files to every output layer.
    v15 – Areal Reduction Factor (ARF) applied to design rainfall via
          CONFIG['ARF_METHOD']: 'NONE' (default), 'TRRL_UK' or
          'RIYADH'; ARF computed from full upstream catchment area at
          each pour point, inherited by all sub-basins downstream.
          New ARF column on Catchments / LFP_PL / SubBasins / SubLFPs
          when active; P_<rp>_yrs values stored already-reduced.
    v14 – Holistic Kirpich Tc/Tlag (TC_KIRPH / TLAG_KIRPH); CONFIG
          TC_METHOD selects which Tc drives the engine; CONFIG
          PEAK_RATE_FACTOR (Gamma SCS UH); per-catchment design-storm
          hydrograph CSVs and QP_<rp>_yrs peak-flow attributes.
    v13 – Optional design-rainfall / SCS-CN runoff / volume columns
          (P_<n>_yrs / R_<n>_yrs / V_<n>yrs_1000m3) on Catchments and
          LFP_PL via CONFIG['RAINFALL_DEPTHS_MM'].
    v12 – Hierarchical naming: W<pid>, W<pid>-Sub<sid>, LFP<pid>-W<pid>, ...
    v11 – BASIN_SLOPE_FROM toggle; Catchments+LFP_PL shared schema.
    v10 – NRCS lag/Tc on LFP_PL & SubLFPs (when CN raster provided).
    v9  – Conditional CN attributes; CN_clipped raster output.
    v8  – Optional Curve Number raster auto-aligned to DEM grid.
    v7  – Reaches mirrors LFP_PL minus TLAG.
    v6  – Optional sub-catchment subdivision (HEC-HMS style).
    v5  – Catchments gained mean basin slope (Horn 1981).
    v4  – LFP_PL gained Kirpich Tc / segment-slope attributes.
    v3  – LFP terminates at the catchment outlet.
    v2  – Pour points auto-reprojected to DEM CRS; 'Snapped' verification layer.
    v1  – Initial pure-NumPy / GDAL implementation (no pip required).

Uses ONLY libraries that ship with QGIS: numpy, osgeo (GDAL), QGIS API.
================================================================================
"""

import os
import re
import heapq
import math
from collections import deque

import numpy as np
from osgeo import gdal, ogr, osr

from qgis.core import (
    QgsProject, QgsRasterLayer, QgsVectorLayer,
    QgsField, QgsFeature, QgsGeometry, QgsPointXY,
    QgsCoordinateReferenceSystem, QgsCoordinateTransform,
)
# v1.0.1 — Qt 5 / Qt 6 compatibility. ``QVariant.Type`` no longer exists
# under PyQt6 (QGIS 4), where QgsField takes a ``QMetaType.Type`` instead.
# ``qt_compat`` resolves the right constant for whichever Qt is running.
try:
    from .qt_compat import (FT_DOUBLE as _FT_DOUBLE, FT_INT as _FT_INT,
                            FT_LONGLONG as _FT_LONGLONG,
                            FT_STRING as _FT_STRING)
except ImportError:  # executed outside the plugin package
    from qt_compat import (FT_DOUBLE as _FT_DOUBLE, FT_INT as _FT_INT,
                           FT_LONGLONG as _FT_LONGLONG,
                           FT_STRING as _FT_STRING)

# Reference numerical tables (SCS storm distributions at 1-min
# resolution). Package-relative import when running as the plugin;
# plain import when this file is executed from the QGIS console with
# scs_data.py sitting next to it; graceful fallback otherwise.
try:
    from .scs_data import SCS_CUM_1MIN as _SCS_CUM_TABLES
except ImportError:
    try:
        from scs_data import SCS_CUM_1MIN as _SCS_CUM_TABLES
    except ImportError:
        _SCS_CUM_TABLES = None   # coarse Type II fallback
try:
    from .scs_data import DIM_UH_TR as _DIM_UH_TR, DIM_UH_QR as _DIM_UH_QR
except ImportError:
    try:
        from scs_data import DIM_UH_TR as _DIM_UH_TR, DIM_UH_QR as _DIM_UH_QR
    except ImportError:
        _DIM_UH_TR = _DIM_UH_QR = None   # Gamma UH fallback

# HEC-HMS .basin writer (v19). Same dual import pattern: package-
# relative when running as the plugin, plain when the file is run
# from the QGIS console with hec_hms_export.py next to it.
try:
    from .hec_hms_export import write_basin_file as _write_hms_basin
except ImportError:
    try:
        from hec_hms_export import write_basin_file as _write_hms_basin
    except ImportError:
        _write_hms_basin = None


# =============================================================================
# 1. CONFIGURATION  ── edit only this block
# =============================================================================
CONFIG = {
    'DEM_LAYER'          : 'FathomDEM',
    'POINTS_LAYER'       : 'PourPoints',
    'POINT_ID_FIELD'     : None,

    # ---- Curve Number (SCS) raster -------------------------------------
    # Name of a CN raster layer already loaded in the QGIS project.
    # Set to None (or '') to skip CN zonal statistics entirely.
    # When set:
    #   • The raster is reprojected on-the-fly to match the DEM grid
    #     exactly (same CRS, extent, and cell size — nearest-neighbour
    #     resampling, since CN is integer-valued).
    #   • An area-weighted mean CN is added to Catchments AND SubBasins
    #     under the field name CN_AVG, plus CN_COV_P (% of basin cells
    #     covered by valid CN values).
    'CN_LAYER'           : None,
    # v1.0 — alternative to CN_LAYER: a single constant Curve Number
    # applied to the whole grid (1–100). When set, it takes precedence
    # over CN_LAYER. Use it when no land-use/soil raster is available.
    'CN_VALUE'           : None,

    # The tool creates an 'ar_tarseem_output' subfolder inside this
    # directory, so C:/temp/ -> C:/temp/ar_tarseem_output/.
    'OUTPUT_DIR'         : r'C:/temp/',

    'DEPRESSION_EPSILON' : 1e-6,
    'STREAM_THRESHOLD'   : 1000,
    'SNAP_RADIUS_CELLS'  : 25,

    # ---- Basin slope source ---------------------------------------------
    # Which DEM to use when computing per-cell terrain slope (used by
    # zonal mean basin slope BAS_SLP / BAS_SLP_P, and by the NRCS lag
    # equation):
    #   'original' : original (unfilled) DEM. Reflects real topography.
    #                Standard choice for the SCS/NRCS empirical formulas.
    #                Default — recommended for most cases.
    #   'filled'   : pit-filled DEM used for flow direction. Self-
    #                consistent with flow paths but flattens depressions
    #                to ε-graded plateaus, which deflates the basin slope.
    'BASIN_SLOPE_FROM'   : 'original',

    # ---- Design rainfall depths for SCS-CN runoff & volume --------------
    # Provide rainfall depths in MILLIMETRES for any subset of the eight
    # standard return periods (2, 5, 10, 20, 25, 50, 100, 200 years).
    #
    #   • Each entry adds three columns to BOTH Catchments and LFP_PL:
    #         P_<n>_yrs        the rainfall depth you provided (mm)
    #         R_<n>_yrs        excess rainfall (= runoff depth, mm)
    #                          via SCS-CN method:
    #                            S  = 25.4 * (1000/CN - 10)        [mm]
    #                            Ia = 0.2 * S                       [mm]
    #                            R  = (P-0.2S)^2 / (P+0.8S)  if P > 0.2S
    #                                 0                       otherwise
    #         V_<n>yrs_1000m3  runoff volume = R * AREA_KM2          [10^3 m³]
    #
    #   • Entries left out (or set to None) are SKIPPED — no columns added.
    #   • R and V are only computed when CN_LAYER is also provided.
    #     When CN is absent but rainfall depths are present, only the
    #     P_<n>_yrs columns are emitted (rainfall depths alone).
    #   • Order of insertion is fixed: 2, 5, 10, 20, 25, 50, 100, 200.
    'RAINFALL_DEPTHS_MM' : {
        # 2   : None,
        # 5   : None,
        # 10  : None,
        # 20  : None,
        # 25  : None,
        # 50  : None,
        # 100 : None,
        # 200 : None,
    },

    # ---- Time of concentration method (v14) -----------------------------
    # Three Tc/Tlag estimates are ALWAYS exposed as attributes on every
    # catchment and longest-flow-path feature (when the inputs allow):
    #
    #   TC_KIRPH / TLAG_KIRPH   Kirpich 1940 — holistic
    #                           Tc = 0.0195 · L^0.77 · S^-0.385
    #                           with L  = full longest flow path length
    #                           and  S  = (ELE_US - ELE_DS) / L
    #
    #   TC_KIRP  / TLAG_KIRP    Kirpich 1940 — segmented (10-85 method)
    #                           Tc = sum of Kirpich Tc on three sub-
    #                           segments along the LFP:
    #                               US → P85 (0-15% from start)
    #                               P85 → P10 (15-90%)
    #                               P10 → DS  (90-100%)
    #
    #   TC_NRCS  / TLAG_NRCS    NRCS (SCS) lag equation — requires CN.
    #                           Tlag = L^0.8 · (2540-22.86·CN)^0.7
    #                                 / (1410 · CN^0.7 · Y^0.5)  · 60
    #                           Only emitted when CN_LAYER is set.
    #
    # TC_METHOD selects which set is consumed by the v14 peak-flow and
    # hydrograph engine. The other two remain visible as attributes.
    # All three are full-fledged outputs — no method is preferred over
    # the others; only the default for the engine is set here.
    # Choices:
    #     'KIRPICH_HOLISTIC'   default. Uses TC_KIRPH / TLAG_KIRPH.
    #     'KIRPICH_10_85'      uses the 10-85 segmented Kirpich form.
    #     'NRCS'               uses NRCS lag (needs CN). If a catchment
    #                          has no valid CN_AVG, that catchment
    #                          falls back to KIRPICH_HOLISTIC silently
    #                          (with a console warning).
    'TC_METHOD'          : 'KIRPICH_HOLISTIC',

    # ---- SCS Unit Hydrograph peak rate factor (PRF) (v14) ---------------
    # The SCS dimensionless unit hydrograph is parameterised by a peak
    # rate factor relating peak height to volume / time-to-peak. The
    # standard SCS curve has PRF = 484 (English) ≡ 0.208 metric. Lower
    # values describe flatter, attenuated responses; higher values
    # describe sharp-peaking, steep catchments.
    #
    # Recommended values (NRCS / SCS guidance):
    #     100-200    swamp / coastal flat
    #     300        moderately flat / mild relief
    #     484        standard mixed-relief watersheds  (DEFAULT)
    #     600        steep / mountainous / urban
    #
    # The tool implements arbitrary PRF via a Gamma-function UH whose
    # shape parameter m is solved numerically so that the dimensionless
    # qp · Tp / V matches the requested PRF. Volume is rigorously
    # conserved across PRF choices — only the hydrograph SHAPE changes.
    #
    # Valid range enforced: 100 ≤ PEAK_RATE_FACTOR ≤ 600. Out-of-range
    # values are clamped with a console warning.
    'PEAK_RATE_FACTOR'   : 484,

    # ---- Hydrograph computation (v14) -----------------------------------
    # When all of the following are true, the tool computes a full
    # design-storm hydrograph for every catchment for every supplied
    # return period:
    #   - HYDROGRAPHS_ENABLED is True
    #   - CN_LAYER is set (so SCS-CN losses can be evaluated)
    #   - RAINFALL_DEPTHS_MM contains at least one positive entry
    #
    # Pipeline:
    #   1. Build a 24-h SCS Type II rainfall hyetograph at 1-min
    #      resolution, scaled to the supplied P_<rp>_yrs depth.
    #   2. Apply SCS-CN cumulative method using CN_AVG to derive the
    #      excess hyetograph.
    #   3. Build an SCS Unit Hydrograph (Gamma form) with
    #          Tp   = TLAG selected by TC_METHOD
    #          PRF  = PEAK_RATE_FACTOR
    #          A    = AREA_KM2 of the catchment
    #   4. Convolve excess hyetograph with the UH to obtain Q(t) at
    #      1-min resolution.
    #
    # Outputs:
    #   - Peak Q per return period appears in the attribute table as
    #         QP_<n>_yrs              [m³/s]
    #     (alongside the existing P_<n>_yrs / R_<n>_yrs / V_<n>yrs_1000m3
    #     columns), on Catchments AND LFP_PL.
    #   - Full hydrograph CSVs are written to
    #         <OUTPUT_DIR>/Hydrographs/W{pid}_{rp}yr.csv
    #     with columns:
    #         time_min , P_mm , Loss_mm , Excess_mm , Q_m3s
    #
    # In v14 sub-basins do NOT receive peak/hydrograph outputs — those
    # require reach-routing convolution which will be added later.
    'HYDROGRAPHS_ENABLED'  : True,
    'STORM_DURATION_HR'    : 24,
    'HYD_TIME_STEP_MIN'    : 1,

    # v1.0 — hydrograph window length (CSV / HTML completeness).
    # Every hydrograph is sized to
    #     storm duration + HYD_TAIL_FACTOR · Tlag · PRF-stretch (+ reach Tc)
    # capped at storm + HYD_MAX_TAIL_HR hours. Increase these values if
    # the console warns that a hydrograph has not fully receded to zero
    # by the end of the output window (both are editable in the dialog).
    'HYD_TAIL_FACTOR'      : 20.0,
    'HYD_MAX_TAIL_HR'      : 96.0,
    # v19: all four standard SCS/NRCS 24-h distributions,
    # at 1-minute tabular resolution:
    #   'SCS_TYPE_I' | 'SCS_TYPE_IA' | 'SCS_TYPE_II' | 'SCS_TYPE_III'
    #   'CUSTOM'  — user-supplied hyetograph read from
    #               CUSTOM_HYETOGRAPH_CSV (v19.1). Two
    #               columns: time (minutes), depth (mm); a header row
    #               is skipped automatically; depths may be either
    #               CUMULATIVE or INCREMENTAL (auto-detected: a
    #               non-decreasing column is treated as cumulative).
    #               The curve is normalised to a mass curve and
    #               rescaled to each return period's total depth P,
    #               before use.
    #               Storm duration is taken from the CSV's last time
    #               (STORM_DURATION_HR is ignored for CUSTOM).
    'STORM_DISTRIBUTION'   : 'SCS_TYPE_II',
    'CUSTOM_HYETOGRAPH_CSV': '',

    # ---- Initial abstraction ratio lambda = Ia/S (v19) -------
    # Classic SCS-CN uses Ia = 0.2*S; newer NRCS guidance (Hawkins et
    # al.) recommends 0.05 for many watersheds. Range 0.01-0.5,
    # clamped with a console warning. Generalised runoff form:
    #     R = (P - lam*S)^2 / (P - lam*S + S)    for P > lam*S
    # which reduces to (P-0.2S)^2/(P+0.8S) at lam = 0.2.
    'SCS_LAMBDA'           : 0.2,
    # v1.0 — lower bound on the time of concentration (minutes) used by
    # the peak-flow / hydrograph engine. 0 = no minimum. The floor is
    # applied to Tc; the matching Tlag = 0.6 x Tc is used with it, so
    # small or steep basins cannot produce an unrealistically sharp
    # unit hydrograph. All three Tc estimates are still written to the
    # attribute tables unchanged.
    'TC_MIN_MINUTES'       : 0.0,

    # ---- Areal Reduction Factor (ARF) (v15) -----------------------------
    # Watershed-area-based reduction applied to all P_<rp>_yrs design
    # depths. Point rainfall (from a gauge) overestimates the area-mean
    # rainfall over a large catchment, since storms are not spatially
    # uniform. The ARF discounts the supplied point depth to obtain
    # the depth that should drive runoff calculations.
    #
    # Choices:
    #   'NONE'      ARF = 1.0 — no reduction. Design depths used as-is.
    #               (Default — preserves v14 behaviour exactly.)
    #
    #   'TRRL_UK'   Transport & Road Research Laboratory, UK
    #               (Riyadh Hydraulic Design Manual ref. {16}):
    #                   ARF = 1 - 0.044 · A^0.275       (A in km²)
    #               Empirical, derived for African watersheds; no
    #               threshold area, falls smoothly from ~0.96 at
    #               1 km² to ~0.70 at 1000 km².
    #
    #   'RIYADH'    Riyadh Development High Commission
    #               (Riyadh Hydraulic Design Manual eq. 2-4):
    #                   ARF = 1.0                        if A ≤ 133 km²
    #                   ARF = 1.4167 - 0.0852 · ln(A)    if A > 133 km²
    #               Calibrated for Riyadh-region watersheds; conservative
    #               (no reduction below 133 km², ~0.83 at 1000 km²,
    #               ~0.70 at 4500 km²).
    #
    # ARF is evaluated from the FULL upstream catchment area at the
    # pour point (AREA_KM2 of the catchment polygon). Every sub-basin
    # downstream of that pour point INHERITS the same ARF — sub-basins
    # do NOT recompute their own ARF from local area, because the
    # design storm is defined over the parent catchment.
    #
    # When ARF_METHOD ≠ 'NONE':
    #   - ARF appears as a new attribute on Catchments / LFP_PL /
    #     SubBasins / SubLFPs.
    #   - The value stored in P_<rp>_yrs IS the reduced depth
    #     (= supplied depth × ARF). All downstream R / V / QP and
    #     hydrograph CSVs are computed from this reduced depth, so
    #     the attribute table remains internally self-consistent
    #     (R = SCS-CN(P_table, CN), V = R · A, etc.).
    #
    # Note: TRRL_UK and RIYADH formulas are duration-independent (they
    # implicitly target the design-storm duration of 24 h that this
    # tool uses). Duration-aware ARF curves (e.g. US Weather Bureau's
    # 30-min / 1-h / 3-h / 6-h / 24-h family) are not implemented in
    # v15; the supplied storm is always treated as the 24-h reference.
    #   'FIXED'  — apply the user-defined constant ARF_FIXED_VALUE
    #              (0 < ARF <= 1) to every catchment regardless of
    #              area (v20).
    'ARF_METHOD'         : 'NONE',
    'ARF_FIXED_VALUE'    : 1.0,

    # ---- Sub-catchment subdivision (v6) ---------------------------------
    # Target sub-basin area in km² (rough; actual sizes vary with the
    # locations of stream junctions). 0 disables subdivision and the tool
    # behaves exactly as v5 (one catchment + one LFP per pour point).
    'SUB_BASIN_AREA_KM2' : 0.0,

    # ---- Reservoirs / dams (v23) -----------------------------------------
    # Optional point layer of dams / reservoirs. Every dam point is
    # snapped to the nearest main flow path (highest flow accumulation
    # stream cell within SNAP_RADIUS_CELLS — the same snapping rule the
    # pour points use) and ENFORCED as a junction during sub-catchment
    # subdivision, so a sub-basin boundary is always created exactly at
    # the dam. A dam layer therefore also FORCES subdivision on, even
    # when SUB_BASIN_AREA_KM2 == 0 (in that case the only junctions are
    # the outlet + the dams).
    #
    #   DAMS_LAYER      — project layer name of the dam points ('' = off)
    #   DAM_NAME_FIELD  — field holding the dam name; the junction at
    #                     the dam is renamed  J{pid}-Sub{sid}(DamName)
    #                     and the name is carried into the HEC-HMS
    #                     .basin junction element as well.
    #   DAM_ATTR_FIELDS — list of dam fields to copy onto the target
    #                     junction feature (each lands in a 'DAM_'-
    #                     prefixed column of the Junctions layer).
    'DAMS_LAYER'         : '',
    'DAM_NAME_FIELD'     : '',
    'DAM_ATTR_FIELDS'    : [],

    # v24 — "Special Delineation Points" extensions:
    #   SP_POINTS_ONLY — when True, sub-catchment division is done at
    #       the special points ONLY (the area-based subdivision from
    #       SUB_BASIN_AREA_KM2 is ignored; junctions = outlet + the
    #       snapped points). When False, the special points are
    #       COMBINED with the Delineation-tab options (Stream
    #       threshold, Snap radius, Sub-basin target area).
    #   SP_RESERVOIRS — {point name: reservoir definition}. Any special
    #       point can be declared a Reservoir; its combined inflow is
    #       then routed through modified-Puls storage routing (same
    #       engine as AR-Awdiya) before continuing downstream, and the
    #       Junctions layer gains inflow/outflow/attenuation/storage/
    #       stage outputs. Definition dict (AR-Awdiya-compatible):
    #         mode       : 'structures' | 'table'
    #         ev_text    : "elev_m, volume_10³m³" rows (structures)
    #         structures : list of
    #             {'type':'spillway','crest':m,'C':1.7,'width':m} and/or
    #             {'type':'pipe','invert':m,'n':1,'radius':m,'Cd':0.6}
    #         table_text : "storage_10³m³, outflow_m³/s" rows (table)
    #         init_mode  : 'inflow'    -> initial outflow equals inflow
    #                      'elevation' -> pool starts at init_elev
    #         init_elev  : initial water-surface elevation (m)
    # v1.0 — AR-Awdiya project export (.json). The exported file opens
    # in AR-Awdiya with "Load project" and is ready to Run: the whole
    # delineated network, its hydrology and the reservoir definitions
    # are carried over.
    #   AWDIYA_FREQ_MODE  'all'    -> frequency analysis over every
    #                                 return period (AR-Awdiya applies
    #                                 each depth to all catchments)
    #                     'single' -> one storm, each catchment keeps
    #                                 its OWN depth (Thiessen variation
    #                                 preserved)
    #   AWDIYA_PRIMARY_RP  return period whose depths fill each
    #                      catchment's P (None = largest available)
    'EXPORT_AWDIYA'      : False,
    'AWDIYA_PROJECT_NAME': '',
    'AWDIYA_FREQ_MODE'   : 'all',
    'AWDIYA_PRIMARY_RP'  : None,

    'SP_POINTS_ONLY'     : False,
    'SP_RESERVOIRS'      : {},

    # v1.0 — Thiessen (Voronoi) polygons of the rainfall stations, the
    # very areas used to weight the station depths. Produced only when a
    # rainfall-stations layer is supplied.
    'OUT_THIESSEN'       : 'Rainfall_Thiessen_Polygons',

    'OUT_STREAMS'        : 'Streams',
    'OUT_CATCHMENTS'     : 'Catchments',
    'OUT_LFP_LINES'      : 'LFP_PL',
    'OUT_LFP_POINTS'     : 'LFP_PTs',
    'OUT_SNAPPED'        : 'Snapped',     # NEW — verification layer
    # Sub-catchment outputs (only populated when SUB_BASIN_AREA_KM2 > 0)
    'OUT_SUBBASINS'      : 'SubBasins',
    'OUT_SUBLFPS'        : 'SubLFPs',
    'OUT_REACHES'        : 'Reaches',
    'OUT_JUNCTIONS'      : 'Junctions',   # v17 — sub-basin outlets

    # ---- Reach routing of hydrographs (v17) -----------------------------
    # When sub-basins are enabled, the hydrograph from each sub-basin
    # outlet (junction) must propagate downstream through one or more
    # reaches before reaching the catchment outlet. Each reach attenuates
    # and delays the inflow.
    #
    # REACH_ROUTING_METHOD chooses the routing approach:
    #   'MUSKINGUM' (default) — classical Muskingum with two parameters:
    #                K [hours] = travel time on the reach, set per
    #                            reach to its computed Tc;
    #                X [-]     = MUSKINGUM_X (default 0.2). 0.0 gives
    #                            full attenuation (linear reservoir),
    #                            0.5 gives pure translation (no
    #                            attenuation). Typical range 0.1-0.3
    #                            for natural rivers.
    #   'LAG'                 — pure translation: shift the inflow
    #                            hydrograph by Tc minutes. No
    #                            attenuation (peak preserved).
    #   'NONE'                — no routing at all; outflow == inflow.
    #                            (Useful for debugging or when reaches
    #                            are very short and routing effects
    #                            are negligible.)
    #   'MUSKINGUM_CUNGE' (v19): variable-parameter
    #                Muskingum-Cunge on a trapezoidal channel; K and X
    #                recomputed each step from Manning hydraulics, the
    #                reach subdivided into N = L/(c*dt) segments (<=40).
    #                Reach length & slope come from the delineated reach
    #                geometry; section shape from the MC_* keys below.
    'REACH_ROUTING_METHOD' : 'MUSKINGUM',
    # v1.0 — lower bound on the reach travel time (minutes) that drives
    # Muskingum K / the Lag shift. 0 = no minimum. Short headwater
    # reaches can otherwise produce travel times of a fraction of a
    # minute, which routes as a near-instantaneous translation.
    'REACH_TC_MIN_MINUTES' : 0.0,

    'MUSKINGUM_X'          : 0.2,
    # v1.0 — Muskingum subreach count applied to EVERY reach:
    #   0  = Automatic (default): per-reach N = ceil(2·K·X/Δt), the
    #        stability-based value (also written to the HEC-HMS
    #        .basin 'Subreaches:' line).
    #   N>0 = fixed count for all reaches. Physically: more subreaches
    #        → less attenuation and shorter travel time per segment
    #        (the cascade approaches pure translation); fewer
    #        subreaches → more attenuation and longer travel time per
    #        segment (approaches a single linear reservoir).
    'MUSK_SUBREACHES'      : 0,

    # ---- Muskingum-Cunge channel geometry (v19) --------------
    'MC_MANNING_N'         : 0.035,   # Manning roughness n
    'MC_BOTTOM_WIDTH_M'    : 5.0,     # bottom width b (m)
    'MC_SIDE_SLOPE'        : 2.0,     # side slope z (H:V)

    # REACH_TC_METHOD selects the reach Tc that drives K in the
    # Muskingum routing (or the lag in 'LAG' routing). All three
    # variants are computed and stored as attributes on every reach;
    # this only chooses which one is consumed:
    #   'KIRPICH_HOLISTIC' (default) — TC_KIRPH on the reach
    #   'KIRPICH_10_85'              — TC_TOTAL  on the reach
    'REACH_TC_METHOD'      : 'KIRPICH_HOLISTIC',

    # ---- QGIS layer styling (v16) ---------------------------------------
    # Directory containing QML style files. When set to a non-empty
    # path, the tool looks up QML files there by name and applies them
    # to every output layer it creates (memory layers added to the
    # project AND raster layers written to disk and re-loaded).
    #
    # Required QML file names (case sensitive; missing files are
    # warned about but never abort the run):
    #     Catchments.qml                            → Catchments
    #     Longest_Flow_Paths.qml                    → LFP_PL
    #     Streams.qml                               → Streams
    #     Pour_Points.qml                           → Snapped
    #     Subbasins.qml                             → SubBasins
    #     Longest_Flow_Paths_-_Subbasins.qml        → SubLFPs
    #     Reaches.qml                               → Reaches
    #     Rainfall_Thiessen_Polygons.qml            → Rainfall_Thiessen_Polygons
    #     DEM.qml                                   → DEM_filled raster
    #
    # When STYLES_DIR is empty (default), no styling is applied and
    # layers receive QGIS's default appearance (v15 behaviour).
    'STYLES_DIR'         : '',

    # ---- HEC-HMS interchange (v19) ---------------------------
    # When EXPORT_HEC_HMS is True the delineated network of every
    # catchment (Subbasins with Area/CN/SCS lag, Reaches with the
    # selected routing method, Junctions, Sink at the pour point) is
    # written as <OUTPUT_DIR>/HEC_HMS/<name>.basin in the exact element
    # format HEC-HMS / WMS reads. Un-subdivided catchments export as a
    # single Subbasin -> Sink pair. HEC_HMS_VERSION: '4.7'..'4.13'
    # (4.11+ writes 'Unregulated Output Required'; older versions the
    # 'Enable Quality Routing' header line, mirroring HEC-HMS itself).
    'EXPORT_HEC_HMS'     : False,
    'HEC_HMS_VERSION'    : '4.12',
    # v1.0 — base name for the exported HEC-HMS files. '' keeps the
    # default 'AR_Tarseem'. The name is used for the .basin, the
    # companion .pdata, the .met files and the basin/project name
    # written inside them, mirroring the AR-Awdiya project-name field.
    'HEC_HMS_NAME'       : '',
    # v22: also write one HEC-HMS .met file per return period
    # (Hypothetical Storm, per-subbasin Storm Depth overrides).
    'EXPORT_HEC_HMS_MET' : False,

    # ---- Spatially-variable rainfall via Thiessen polygons (v22) --------
    # Optional point layer of rainfall stations. When set, rainfall
    # depths per return period are read from the layer's attribute
    # fields (RAIN_STATIONS_FIELDS maps return period -> field name)
    # and every catchment / sub-basin receives its own area-weighted
    # depth from a Thiessen (nearest-station / Voronoi) assignment of
    # the basin's DEM cells. The uniform RAINFALL_DEPTHS_MM table is
    # ignored while this mode is active (its return-period keys are
    # replaced by the mapping's keys).
    'RAIN_STATIONS_LAYER'  : '',
    'RAIN_STATIONS_FIELDS' : {},     # e.g. {100: 'P100', 50: 'P50'}
}


# =============================================================================
# 2. D8 CONSTANTS
# =============================================================================
D8_CODES = [
    ( 64, -1,  0, 1.0),
    (128, -1,  1, math.sqrt(2)),
    (  1,  0,  1, 1.0),
    (  2,  1,  1, math.sqrt(2)),
    (  4,  1,  0, 1.0),
    (  8,  1, -1, math.sqrt(2)),
    ( 16,  0, -1, 1.0),
    ( 32, -1, -1, math.sqrt(2)),
]
CODE_TO_OFFSET = {c: (dr, dc) for c, dr, dc, _ in D8_CODES}
INV_REQUIRED_CODE = {
    (-1,  0):   4, (-1,  1):   8, ( 0,  1):  16, ( 1,  1):  32,
    ( 1,  0):  64, ( 1, -1): 128, ( 0, -1):   1, (-1, -1):   2,
}


# =============================================================================
# 3. RASTER I/O HELPERS
# =============================================================================
def read_dem(path):
    ds = gdal.Open(path, gdal.GA_ReadOnly)
    if ds is None:
        raise IOError(f"Cannot open {path}")
    band   = ds.GetRasterBand(1)
    arr    = band.ReadAsArray().astype(np.float64)
    nodata = band.GetNoDataValue()
    gt     = ds.GetGeoTransform()
    proj   = ds.GetProjection()
    ds = None
    return arr, gt, proj, nodata


def write_raster(arr, gt, proj, out_path,
                 dtype=gdal.GDT_Float32, nodata=None):
    rows, cols = arr.shape
    drv = gdal.GetDriverByName('GTiff')
    ds  = drv.Create(out_path, cols, rows, 1, dtype,
                     options=['COMPRESS=DEFLATE', 'TILED=YES'])
    ds.SetGeoTransform(gt)
    ds.SetProjection(proj)
    band = ds.GetRasterBand(1)
    band.WriteArray(arr)
    if nodata is not None:
        band.SetNoDataValue(float(nodata))
    band.FlushCache()
    ds = None
    return out_path


def cell_to_xy(row, col, gt, centre=True):
    off = 0.5 if centre else 0.0
    x = gt[0] + (col + off) * gt[1] + (row + off) * gt[2]
    y = gt[3] + (col + off) * gt[4] + (row + off) * gt[5]
    return x, y


def xy_to_cell(x, y, gt):
    col = int((x - gt[0]) / gt[1])
    row = int((y - gt[3]) / gt[5])
    return row, col


# =============================================================================
# 4. HYDROLOGY PRIMITIVES
# =============================================================================
def priority_flood_fill(dem, nodata=None, epsilon=1e-6):
    rows, cols = dem.shape
    filled = dem.copy().astype(np.float64)
    closed = np.zeros((rows, cols), dtype=bool)
    if nodata is not None:
        closed[dem == nodata] = True
        filled[dem == nodata] = np.nan

    pq = []
    counter = 0
    for r in range(rows):
        for c in (0, cols - 1):
            if not closed[r, c]:
                heapq.heappush(pq, (filled[r, c], counter, r, c))
                counter += 1
                closed[r, c] = True
    for c in range(1, cols - 1):
        for r in (0, rows - 1):
            if not closed[r, c]:
                heapq.heappush(pq, (filled[r, c], counter, r, c))
                counter += 1
                closed[r, c] = True

    while pq:
        elev, _, r, c = heapq.heappop(pq)
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                nr, nc = r + dr, c + dc
                if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
                    continue
                if closed[nr, nc]:
                    continue
                new_elev = max(filled[nr, nc], elev + epsilon)
                filled[nr, nc] = new_elev
                heapq.heappush(pq, (new_elev, counter, nr, nc))
                counter += 1
                closed[nr, nc] = True
    return filled


def flow_direction_d8(dem, nodata_mask=None):
    rows, cols = dem.shape
    fdir       = np.zeros((rows, cols), dtype=np.uint8)
    max_slope  = np.full((rows, cols), -np.inf)
    for code, dr, dc, dist in D8_CODES:
        r0 = max(0, -dr); r1 = rows - max(0, dr)
        c0 = max(0, -dc); c1 = cols - max(0, dc)
        src = dem[r0:r1, c0:c1]
        tgt = dem[r0 + dr : r1 + dr, c0 + dc : c1 + dc]
        slope = (src - tgt) / dist
        cur = max_slope[r0:r1, c0:c1]
        better = slope > cur
        max_slope[r0:r1, c0:c1] = np.where(better, slope, cur)
        fdir[r0:r1, c0:c1] = np.where(better, code,
                                      fdir[r0:r1, c0:c1]).astype(np.uint8)
    if nodata_mask is not None:
        fdir[nodata_mask] = 0
    return fdir


def flow_accumulation(fdir):
    rows, cols = fdir.shape
    indeg = np.zeros((rows, cols), dtype=np.int32)
    for code, dr, dc, _ in D8_CODES:
        rs, cs = np.where(fdir == code)
        trs, tcs = rs + dr, cs + dc
        valid = (trs >= 0) & (trs < rows) & (tcs >= 0) & (tcs < cols)
        np.add.at(indeg, (trs[valid], tcs[valid]), 1)

    facc  = np.ones((rows, cols), dtype=np.int64)
    queue = deque()
    rs0, cs0 = np.where(indeg == 0)
    for r, c in zip(rs0.tolist(), cs0.tolist()):
        queue.append((r, c))

    while queue:
        r, c = queue.popleft()
        d = int(fdir[r, c])
        if d not in CODE_TO_OFFSET:
            continue
        dr, dc = CODE_TO_OFFSET[d]
        nr, nc = r + dr, c + dc
        if 0 <= nr < rows and 0 <= nc < cols:
            facc[nr, nc] += facc[r, c]
            indeg[nr, nc] -= 1
            if indeg[nr, nc] == 0:
                queue.append((nr, nc))
    return facc


def delineate_catchment(fdir, outlet_row, outlet_col):
    rows, cols = fdir.shape
    catch = np.zeros((rows, cols), dtype=bool)
    catch[outlet_row, outlet_col] = True
    q = deque([(outlet_row, outlet_col)])
    while q:
        r, c = q.popleft()
        for (dr, dc), need_code in INV_REQUIRED_CODE.items():
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                if not catch[nr, nc] and fdir[nr, nc] == need_code:
                    catch[nr, nc] = True
                    q.append((nr, nc))
    return catch


def flow_distance(fdir, outlet_row, outlet_col, dx, dy):
    rows, cols = fdir.shape
    dist = np.full((rows, cols), np.nan)
    dist[outlet_row, outlet_col] = 0.0
    q = deque([(outlet_row, outlet_col)])
    diag = math.sqrt(dx * dx + dy * dy)
    while q:
        r, c = q.popleft()
        cur = dist[r, c]
        for (dr, dc), need_code in INV_REQUIRED_CODE.items():
            nr, nc = r + dr, c + dc
            if 0 <= nr < rows and 0 <= nc < cols:
                if np.isnan(dist[nr, nc]) and fdir[nr, nc] == need_code:
                    if dr != 0 and dc != 0:
                        seg = diag
                    elif dr != 0:
                        seg = dy
                    else:
                        seg = dx
                    dist[nr, nc] = cur + seg
                    q.append((nr, nc))
    return dist


def trace_flow_path(fdir, row, col, max_steps=200000,
                    stop_row=None, stop_col=None):
    """
    Follow the D8 chain from (row, col) downstream.
    Stops when:
        - the path reaches (stop_row, stop_col) — the basin outlet, OR
        - the cell has no valid flow direction (sink), OR
        - the path leaves the grid, OR
        - a loop is detected (safety net), OR
        - max_steps is reached.
    Passing the outlet via stop_row/stop_col ensures the LFP terminates at
    the catchment outlet and does NOT continue downstream past it.
    """
    rows, cols = fdir.shape
    path = [(row, col)]
    seen = {(row, col)}
    if stop_row is not None and row == stop_row and col == stop_col:
        return path
    for _ in range(max_steps):
        d = int(fdir[row, col])
        if d not in CODE_TO_OFFSET:
            break
        dr, dc = CODE_TO_OFFSET[d]
        nr, nc = row + dr, col + dc
        if nr < 0 or nr >= rows or nc < 0 or nc >= cols:
            break
        if (nr, nc) in seen:
            break
        row, col = nr, nc
        path.append((row, col))
        seen.add((row, col))
        if stop_row is not None and row == stop_row and col == stop_col:
            break
    return path


def snap_to_stream(streams, facc, row, col, radius):
    rows, cols = streams.shape
    r0 = max(0, row - radius); r1 = min(rows, row + radius + 1)
    c0 = max(0, col - radius); c1 = min(cols, col + radius + 1)
    win_s = streams[r0:r1, c0:c1]
    if win_s.sum() == 0:
        return row, col, False
    win_f = facc[r0:r1, c0:c1]
    masked = np.where(win_s > 0, win_f, -1)
    idx = np.unravel_index(int(np.argmax(masked)), masked.shape)
    return r0 + idx[0], c0 + idx[1], True


def horn_slope(dem, dx, dy, nodata_mask=None):
    """
    Per-cell terrain slope using Horn's algorithm (1981) — the standard
    used by ArcGIS 'Slope' and 'gdaldem slope'.

    For a 3×3 window
        a b c
        d e f
        g h i
    centred on cell e:
        dz/dx = ((c + 2f + i) - (a + 2d + g)) / (8 * dx)
        dz/dy = ((g + 2h + i) - (a + 2b + c)) / (8 * dy)
        slope_mm = sqrt((dz/dx)^2 + (dz/dy)^2)        # tan(angle), m/m
    Edge cells (1-pixel border) are returned as NaN. Nodata cells are NaN.
    """
    rows, cols = dem.shape
    slope = np.full((rows, cols), np.nan, dtype=np.float64)
    if rows < 3 or cols < 3:
        return slope

    a = dem[0:-2, 0:-2]; b = dem[0:-2, 1:-1]; c = dem[0:-2, 2:  ]
    d = dem[1:-1, 0:-2];                      f = dem[1:-1, 2:  ]
    g = dem[2:  , 0:-2]; h = dem[2:  , 1:-1]; i = dem[2:  , 2:  ]

    dzdx = ((c + 2.0 * f + i) - (a + 2.0 * d + g)) / (8.0 * dx)
    dzdy = ((g + 2.0 * h + i) - (a + 2.0 * b + c)) / (8.0 * dy)
    slope[1:-1, 1:-1] = np.sqrt(dzdx * dzdx + dzdy * dzdy)

    if nodata_mask is not None:
        slope[nodata_mask] = np.nan
        # also blank out cells whose 3×3 window touches a nodata cell
        nd = nodata_mask
        nd_neighbour = (
            nd[0:-2, 0:-2] | nd[0:-2, 1:-1] | nd[0:-2, 2:  ] |
            nd[1:-1, 0:-2] | nd[1:-1, 1:-1] | nd[1:-1, 2:  ] |
            nd[2:  , 0:-2] | nd[2:  , 1:-1] | nd[2:  , 2:  ]
        )
        sub = slope[1:-1, 1:-1]
        sub[nd_neighbour] = np.nan
        slope[1:-1, 1:-1] = sub
    return slope


def basin_mean_slope(slope_grid, catch_mask):
    """
    Mean basin slope (m/m) over all valid (non-NaN) catchment cells.
    Returns 0.0 if no valid cells.
    """
    vals = slope_grid[catch_mask]
    vals = vals[~np.isnan(vals)]
    if vals.size == 0:
        return 0.0
    return float(np.mean(vals))


def reproject_to_match_dem(src_path, dem_gt, dem_proj, dem_shape,
                           dst_path, resample='near'):
    """
    Reproject and resample a raster so it has IDENTICAL grid alignment
    to the DEM: same CRS, same origin, same cell size, same shape.
    Uses GDAL Warp with nearest-neighbour by default (correct choice
    for categorical / integer-valued rasters like CN). Writes a GeoTIFF
    at dst_path and returns the path.
    """
    rows, cols = dem_shape
    x0    = dem_gt[0]
    y0    = dem_gt[3]
    px    = dem_gt[1]
    py    = dem_gt[5]                          # negative for north-up
    x_max = x0 + cols * px
    y_min = y0 + rows * py
    out_bounds = (min(x0, x_max), min(y0, y_min),
                  max(x0, x_max), max(y0, y_min))

    warp_opts = gdal.WarpOptions(
        format='GTiff',
        dstSRS=dem_proj,
        outputBounds=out_bounds,
        width=cols, height=rows,
        resampleAlg=resample,
        creationOptions=['COMPRESS=DEFLATE', 'TILED=YES'],
        multithread=True,
    )
    out_ds = gdal.Warp(dst_path, src_path, options=warp_opts)
    if out_ds is None:
        raise RuntimeError(f"gdal.Warp failed for {src_path}")
    out_ds = None
    return dst_path


def basin_mean_cn(cn_grid, cn_nodata, catch_mask):
    """
    Area-weighted arithmetic mean CN over the catchment.
    Excludes cells where CN equals the source nodata value (or is NaN).
    Returns (cn_mean, coverage_pct):
        cn_mean        average CN, or 0.0 if no valid cells
        coverage_pct   percentage of catchment cells with a valid CN
    """
    n_total = int(catch_mask.sum())
    if n_total == 0:
        return 0.0, 0.0
    vals = cn_grid[catch_mask].astype(np.float64)
    if cn_nodata is not None:
        valid = vals != cn_nodata
    else:
        valid = np.ones_like(vals, dtype=bool)
    valid &= ~np.isnan(vals)
    if not valid.any():
        return 0.0, 0.0
    cn_mean      = float(vals[valid].mean())
    coverage_pct = 100.0 * float(valid.sum()) / float(n_total)
    return cn_mean, coverage_pct


# =============================================================================
# 5. STREAM POLYLINES
# =============================================================================
def streams_to_polylines(fdir, streams, facc, gt, min_length_cells=2):
    rows, cols = streams.shape
    visited = np.zeros_like(streams, dtype=bool)
    rs, cs = np.where(streams.astype(bool))
    if len(rs) == 0:
        return []
    order = np.argsort(facc[rs, cs])
    rs, cs = rs[order], cs[order]

    polylines = []
    for r0, c0 in zip(rs.tolist(), cs.tolist()):
        if visited[r0, c0]:
            continue
        r, c = r0, c0
        line = []
        while True:
            line.append((r, c))
            visited[r, c] = True
            d = int(fdir[r, c])
            if d not in CODE_TO_OFFSET:
                break
            dr, dc = CODE_TO_OFFSET[d]
            nr, nc = r + dr, c + dc
            if not (0 <= nr < rows and 0 <= nc < cols):
                break
            if not streams[nr, nc]:
                break
            if visited[nr, nc]:
                line.append((nr, nc))
                break
            r, c = nr, nc
        if len(line) >= min_length_cells:
            polylines.append([cell_to_xy(rr, cc, gt) for rr, cc in line])
    return polylines


# =============================================================================
# 6. POLYGONIZE A BOOLEAN MASK
# =============================================================================
def mask_to_qgs_geometry(mask, gt, proj):
    rows, cols = mask.shape
    src_ds = gdal.GetDriverByName('MEM').Create('', cols, rows, 1, gdal.GDT_Byte)
    src_ds.SetGeoTransform(gt)
    src_ds.SetProjection(proj)
    src_band = src_ds.GetRasterBand(1)
    src_band.WriteArray(mask.astype(np.uint8))
    src_band.SetNoDataValue(0)

    drv  = ogr.GetDriverByName('Memory')
    vds  = drv.CreateDataSource('out')
    srs  = osr.SpatialReference()
    if proj:
        srs.ImportFromWkt(proj)
    layer = vds.CreateLayer('m', srs, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn('val', ogr.OFTInteger))
    gdal.Polygonize(src_band, src_band, layer, 0, [])

    geoms = []
    layer.ResetReading()
    for feat in layer:
        wkt = feat.GetGeometryRef().ExportToWkt()
        geoms.append(QgsGeometry.fromWkt(wkt))
    src_ds = None
    vds    = None
    if not geoms:
        return None
    if len(geoms) == 1:
        return geoms[0]
    out = geoms[0]
    for g in geoms[1:]:
        out = out.combine(g)
    return out


# =============================================================================
# 7. LFP CHARACTERISTICS
# =============================================================================
def sample_dem(xy, dem_layer):
    if xy is None:
        return None
    val, ok = dem_layer.dataProvider().sample(xy, 1)
    return val if ok else None


def kirpich_tc(length_m, slope_mm):
    """
    Kirpich (1940) time of concentration — SI form.
        Tc[min] = 0.0195 * L[m]^0.77 * S^(-0.385)
    L in metres, S as a dimensionless decimal (m/m).
    Returns 0.0 (with caller-side warning) if length<=0 or slope<=0,
    because the formula is undefined for non-positive slope.
    """
    if length_m <= 0 or slope_mm <= 0:
        return 0.0
    return 0.0195 * (length_m ** 0.77) * (slope_mm ** -0.385)


def nrcs_lag(length_m, cn, basin_slope_pct):
    """
    NRCS (SCS) lag-time equation — metric form, returns minutes.
        Tlag[min] = L^0.8 * (2540 - 22.86*CN)^0.7
                    / (1410 * CN^0.7 * Y^0.5)  *  60
    where:
        L  = longest flow path length in METRES
        CN = curve number (dimensionless, valid range ~30-100)
        Y  = mean basin slope in PERCENT (e.g. 4.5, not 0.045)

    NRCS guidance: most reliable for watersheds < 10 km², slopes ≤ 5%,
    and CN between 50-95. Outside these bounds the formula is still
    evaluated but results should be treated with caution.

    Returns 0.0 when any input is non-positive, when CN is non-physical
    (CN > 100), or when (2540 - 22.86*CN) <= 0.
    """
    if length_m <= 0 or cn <= 0 or basin_slope_pct <= 0:
        return 0.0
    if cn > 100:                          # non-physical CN
        return 0.0
    s_term = 2540.0 - 22.86 * cn          # potential maximum retention term
    if s_term <= 0:
        return 0.0
    numerator   = (length_m ** 0.8) * (s_term ** 0.7)
    denominator = 1410.0 * (cn ** 0.7) * (basin_slope_pct ** 0.5)
    if denominator <= 0:
        return 0.0
    return (numerator / denominator) * 60.0


def scs_cn_runoff(p_mm, cn, lam=0.2):
    """
    SCS-CN excess rainfall (direct runoff depth) — metric form, returns mm.

        S  = 25.4 * (1000/CN - 10)                        [mm]
        Ia = 0.2 * S                                       [mm]
        R  = (P - 0.2*S)^2 / (P + 0.8*S)   if P > 0.2*S    [mm]
             0                              otherwise

    Below the initial-abstraction threshold (P <= Ia) no runoff is
    generated. Above it, the curve number method assumes the ratio
    of actual retention to potential maximum retention equals the
    ratio of actual runoff to potential maximum runoff.

    Parameters
    ----------
    p_mm : float
        Rainfall depth in millimetres.
    cn : float
        Runoff curve number (0 < CN ≤ 100). Non-physical values are
        rejected (returns 0).

    Returns
    -------
    R : float
        Excess rainfall in millimetres (≥ 0).
    """
    if p_mm is None or cn is None:
        return 0.0
    if p_mm <= 0 or cn <= 0 or cn > 100:
        return 0.0
    lam  = min(0.5, max(0.01, float(lam) if lam is not None else 0.2))
    s_mm = 25.4 * (1000.0 / cn - 10.0)
    ia   = lam * s_mm
    if p_mm <= ia:
        return 0.0
    # Generalised SCS-CN: R = (P-Ia)^2/(P-Ia+S); identical to the
    # classic (P-0.2S)^2/(P+0.8S) at lam = 0.2.
    return ((p_mm - ia) ** 2) / (p_mm - ia + s_mm)


# =============================================================================
# 7a. AREAL REDUCTION FACTOR (ARF) (v15)
# =============================================================================
def arf_trrl_uk(area_km2):
    """
    Transport and Road Research Laboratory, UK areal reduction factor.

        ARF = 1 - 0.044 · A^0.275          [A in km²]

    Empirical formula derived for African watersheds (TRRL Road Note
    31; Riyadh Hydraulic Design Manual reference {16}). Result is
    clamped to [0, 1] for numerical safety on extreme inputs.

    Sample values:
        A = 1     km² → 0.956
        A = 100   km² → 0.844
        A = 1000  km² → 0.703
        A = 10000 km² → 0.446
    """
    if area_km2 is None or area_km2 <= 0:
        return 1.0
    arf = 1.0 - 0.044 * (float(area_km2) ** 0.275)
    return max(0.0, min(1.0, arf))


def arf_riyadh(area_km2):
    """
    Riyadh ARF (Saudi Arabia) per the Riyadh Development High
    Commission technical specification (Riyadh Hydraulic Design
    Manual equation 2-4):

        ARF = 1                                if A ≤ 133 km²
            = 1.4167 - 0.0852 · ln(A)          if A > 133 km²

    The 133 km² no-reduction threshold and the linear-in-ln(A)
    recession were fitted by the Commission's experts for Riyadh-
    region watersheds, with deliberately conservative behaviour
    (no reduction at all on small catchments).

    Sample values matching the source document:
        A = 100   km² → 1.000  (below threshold)
        A = 1000  km² → 0.828  (document-cited example)
        A = 2500  km² → 0.750  (document-cited point)
        A = 4500  km² → 0.700
    """
    if area_km2 is None or area_km2 <= 0:
        return 1.0
    A = float(area_km2)
    if A <= 133.0:
        return 1.0
    arf = 1.4167 - 0.0852 * math.log(A)
    return max(0.0, min(1.0, arf))


def compute_arf(area_km2, method, fixed_value=1.0):
    """
    Dispatcher returning the ARF for a given catchment area and method.

    Parameters
    ----------
    area_km2 : float
        Full upstream catchment area at the pour point, in km².
    method : str
        One of {'NONE', 'TRRL_UK', 'RIYADH', 'FIXED'},
        case-insensitive. 'FIXED' returns ``fixed_value`` clamped
        to (0, 1] for every catchment, ignoring the area.
        Anything else (including None) is treated as 'NONE' and
        returns 1.0 silently.

    Returns
    -------
    arf : float
        Areal reduction factor, in [0, 1].
    """
    if method is None:
        return 1.0
    m = str(method).upper().strip()
    if m == 'NONE':
        return 1.0
    if m == 'TRRL_UK':
        return arf_trrl_uk(area_km2)
    if m == 'RIYADH':
        return arf_riyadh(area_km2)
    if m == 'FIXED':
        try:
            v = float(fixed_value)
        except (TypeError, ValueError):
            v = 1.0
        if not (0.0 < v <= 1.0):
            print(f"   ⚠ ARF_FIXED_VALUE={fixed_value!r} outside "
                  f"(0, 1] — using 1.0.")
            v = 1.0
        return v
    return 1.0


# =============================================================================
# 7a-bis. QGIS LAYER STYLING (v16)
# =============================================================================
# When CONFIG['STYLES_DIR'] points to a directory holding QML files,
# the tool applies them to its outputs. The mappings below are fixed
# and reflect the supplied QML filenames; users do NOT need to edit
# them — only set STYLES_DIR.

# Vector layers: maps the CONFIG key holding the layer name to the
# QML file basename. None means "no QML supplied; leave unstyled".
QML_STYLE_MAP_VECTOR = {
    'OUT_CATCHMENTS' : 'Catchments.qml',
    'OUT_LFP_LINES'  : 'Longest_Flow_Paths.qml',
    'OUT_LFP_POINTS' : None,
    'OUT_STREAMS'    : 'Streams.qml',
    'OUT_SNAPPED'    : 'Pour_Points.qml',
    'OUT_SUBBASINS'  : 'Subbasins.qml',
    'OUT_SUBLFPS'    : 'Longest_Flow_Paths_-_Subbasins.qml',
    'OUT_REACHES'    : 'Reaches.qml',
    'OUT_JUNCTIONS'  : 'Junctions.qml',    # v17.1 — supplied QML
    'OUT_THIESSEN'   : 'Rainfall_Thiessen_Polygons.qml',   # v1.0
}

# Raster layers: maps the layer NAME (as added to the project) to the
# QML basename. The DEM-style is applied to DEM_filled (a continuous-
# elevation raster, identical visual category to the user's DEM).
QML_STYLE_MAP_RASTER = {
    'DEM_filled' : 'DEM.qml',
    # 'FlowAcc'    : None,    # no QML supplied
    # 'CN_clipped' : None,    # no QML supplied
}


def apply_qml_style(layer, styles_dir, qml_basename, verbose=True):
    """Apply a QML style file to a QgsMapLayer.

    Soft-fails (prints a warning, returns False) when:
      - styles_dir is empty
      - qml_basename is None
      - the file does not exist on disk
      - QGIS reports an error parsing / applying the style

    On success, calls layer.triggerRepaint() so the change is
    visible immediately when the layer is already in the project.
    Never raises — failures are logged and processing continues.
    """
    if not styles_dir:
        return False
    if not qml_basename:
        return False
    qml_path = os.path.join(styles_dir, qml_basename)
    if not os.path.isfile(qml_path):
        if verbose:
            print(f"   ⚠ style file not found: {qml_path}")
        return False
    try:
        msg, ok = layer.loadNamedStyle(qml_path)
    except Exception as e:
        if verbose:
            print(f"   ⚠ style load raised for '{layer.name()}': {e}")
        return False
    if not ok:
        if verbose:
            print(f"   ⚠ style load failed for '{layer.name()}': {msg}")
        return False
    try:
        layer.triggerRepaint()
    # B110: styling is cosmetic — a missing or invalid QML must not fail the run
    except Exception:  # nosec B110
        pass
    if verbose:
        print(f"   ✓ styled '{layer.name()}' with {qml_basename}")
    return True


def style_vector_layer(layer, cfg_key, styles_dir, verbose=True):
    """Convenience wrapper: look up the QML for a CONFIG output key and apply.
    cfg_key is a key in QML_STYLE_MAP_VECTOR (e.g. 'OUT_CATCHMENTS').
    """
    qml = QML_STYLE_MAP_VECTOR.get(cfg_key)
    return apply_qml_style(layer, styles_dir, qml, verbose=verbose)


def style_raster_layer(layer, layer_name_key, styles_dir, verbose=True):
    """Convenience wrapper for raster layers."""
    qml = QML_STYLE_MAP_RASTER.get(layer_name_key)
    return apply_qml_style(layer, styles_dir, qml, verbose=verbose)


# =============================================================================
# 7b. HOLISTIC KIRPICH + SCS DESIGN-STORM HYDROGRAPH (v14)
# =============================================================================
def kirpich_tc_holistic(length_m, ele_us, ele_ds):
    """
    Holistic Kirpich (1940) Tc — treats the longest flow path as a
    single segment. Uses the full path length and the end-to-end
    average slope:
        S       = (ELE_US - ELE_DS) / L     [m/m]
        Tc[min] = 0.0195 * L[m]^0.77 * S^-0.385
    Returns 0.0 if the path is degenerate (L <= 0 or non-positive
    drop).
    """
    if length_m <= 0:
        return 0.0
    drop = float(ele_us) - float(ele_ds)
    if drop <= 0:
        return 0.0
    slope = drop / length_m
    return 0.0195 * (length_m ** 0.77) * (slope ** -0.385)


# ---- SCS Type II 24-hour rainfall distribution (cumulative fraction) -----
# Standard NRCS NEH-630 / TR-55 tabulation. Values read at 0.5-h
# resolution where the storm is gentle and at finer resolution near the
# peak burst (hours 9.5-13.5). Cumulative depth at t = SCS_TYPE_II_HOURS[i]
# divided by total storm depth equals SCS_TYPE_II_CUM[i].
SCS_TYPE_II_HOURS = (
    0.0,  0.5,  1.0,  1.5,  2.0,  2.5,  3.0,  3.5,  4.0,  4.5,
    5.0,  5.5,  6.0,  6.5,  7.0,  7.5,  8.0,  8.5,  9.0,  9.5,
    9.75, 10.0, 10.5, 10.75, 11.0, 11.5, 11.75, 12.0,
   12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 16.5,
   17.0, 17.5, 18.0, 18.5, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0,
)
SCS_TYPE_II_CUM = (
    0.000, 0.005, 0.011, 0.017, 0.022, 0.029, 0.035, 0.041, 0.048, 0.056,
    0.064, 0.072, 0.080, 0.090, 0.098, 0.110, 0.120, 0.133, 0.147, 0.163,
    0.172, 0.181, 0.204, 0.218, 0.235, 0.283, 0.357, 0.663,
    0.735, 0.772, 0.799, 0.820, 0.835, 0.850, 0.860, 0.880, 0.892,
    0.902, 0.913, 0.922, 0.930, 0.938, 0.952, 0.964, 0.976, 0.988, 1.000,
)


def load_custom_hyetograph(csv_path):
    """Load a user hyetograph CSV -> dict(T_min, F, duration_hr) or None.

    v19.1 (custom rainfall). Format: two comma/semicolon/
    tab-separated columns, time (minutes) and depth (mm). A non-numeric
    first row (header) is skipped. Depth may be cumulative or
    incremental — auto-detected: a non-decreasing column is treated as
    cumulative, anything else is cumulated. The result is normalised to
    a dimensionless mass curve F in [0, 1]; the total depth applied at
    run time is each return period's P (mm), so a frequency-run depth
    rescales the custom curve.
    """
    if not csv_path:
        return None
    try:
        rows = []
        with open(csv_path, 'r', encoding='utf-8-sig') as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                for sep in (',', ';', '\t'):
                    if sep in line:
                        parts = line.split(sep)
                        break
                else:
                    parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    t = float(parts[0].strip())
                    d = float(parts[1].strip())
                except ValueError:
                    continue                    # header or junk row
                rows.append((t, d))
        if len(rows) < 2:
            print(f"   \u26a0 custom hyetograph '{csv_path}': fewer than "
                  f"2 numeric rows — ignored.")
            return None
        rows.sort(key=lambda r: r[0])
        t_arr = np.asarray([r[0] for r in rows], dtype=float)
        d_arr = np.asarray([r[1] for r in rows], dtype=float)
        if np.any(d_arr < 0):
            print(f"   \u26a0 custom hyetograph: negative depths — "
                  f"clipped to 0.")
            d_arr = np.clip(d_arr, 0.0, None)
        # auto-detect cumulative vs incremental
        if np.all(np.diff(d_arr) >= -1e-9):
            cum = d_arr.copy()
            kind = 'cumulative'
        else:
            cum = np.cumsum(d_arr)
            kind = 'incremental'
        # anchor at t=0 with zero depth
        if t_arr[0] > 0:
            t_arr = np.insert(t_arr, 0, 0.0)
            cum   = np.insert(cum,   0, 0.0)
        total = float(cum[-1])
        if total <= 0:
            print(f"   \u26a0 custom hyetograph: zero total depth — "
                  f"ignored.")
            return None
        F = np.clip(cum / total, 0.0, 1.0)
        F = np.maximum.accumulate(F)             # enforce monotone
        dur_hr = float(t_arr[-1]) / 60.0
        print(f"   Custom hyetograph: {len(t_arr)} points, "
              f"{kind} depths, duration {dur_hr:.2f} h, "
              f"native total {total:.1f} mm (rescaled to each "
              f"return-period P)")
        return dict(T_min=t_arr, F=F, duration_hr=dur_hr)
    except OSError as e:
        print(f"   \u26a0 custom hyetograph '{csv_path}' unreadable: {e}")
        return None


def design_storm_cum_fraction(t_minutes, distribution='SCS_TYPE_II',
                              duration_hr=24.0, custom_curve=None):
    """
    Cumulative rainfall fraction at given times (in minutes) for the
    selected design-storm distribution.

    Parameters
    ----------
    t_minutes : 1-D numpy array
        Times in minutes from storm start.
    distribution : str
        Currently only 'SCS_TYPE_II' is supported.
    duration_hr : float
        Storm duration in hours (default 24). The dimensionless curve
        is stretched to this duration.

    Returns
    -------
    F : 1-D numpy array
        Cumulative fraction in [0, 1] at the requested times.
    """
    distribution = (distribution or 'SCS_TYPE_II').upper().strip()
    if distribution == 'CUSTOM':
        if custom_curve is not None:
            t_m = np.asarray(t_minutes, dtype=float)
            F = np.interp(t_m, custom_curve['T_min'], custom_curve['F'],
                          left=0.0, right=1.0)
            return np.clip(F, 0.0, 1.0)
        print("   \u26a0 STORM_DISTRIBUTION='CUSTOM' but no usable "
              "hyetograph — using SCS_TYPE_II.")
        distribution = 'SCS_TYPE_II'
    valid = ('SCS_TYPE_I', 'SCS_TYPE_IA', 'SCS_TYPE_II', 'SCS_TYPE_III')
    if distribution not in valid:
        print(f"   \u26a0 STORM_DISTRIBUTION='{distribution}' not "
              f"recognised \u2014 using SCS_TYPE_II.")
        distribution = 'SCS_TYPE_II'
    if _SCS_CUM_TABLES is not None and distribution in _SCS_CUM_TABLES:
        # Exact reference tables: cumulative fraction at 1-min
        # resolution over the 24-h storm (1441 points).
        table_h = np.arange(1441, dtype=float) / 60.0     # hours
        table_F = np.asarray(_SCS_CUM_TABLES[distribution], dtype=float)
    else:
        if distribution != 'SCS_TYPE_II':
            print("   \u26a0 scs_data module unavailable \u2014 "
                  "falling back to the coarse Type II table.")
        table_h = np.asarray(SCS_TYPE_II_HOURS, dtype=float)
        table_F = np.asarray(SCS_TYPE_II_CUM,   dtype=float)
    # rescale dimensionless time-axis from 24 h to the requested duration
    table_h_scaled = table_h * (duration_hr / 24.0)
    t_h = np.asarray(t_minutes, dtype=float) / 60.0
    F = np.interp(t_h, table_h_scaled, table_F, left=0.0, right=1.0)
    return np.clip(F, 0.0, 1.0)


def _shape_const_for_m(m):
    """qp · Tp / V for a Gamma UH with shape parameter m. Proportional
    to the SCS peak rate factor (PRF). Strictly increasing in m for
    m > 0."""
    if m <= 0:
        return 0.0
    return (m ** (m + 1.0)) / (math.exp(m) * math.gamma(m + 1.0))


# Cached: shape constant of the standard SCS UH (m = 3.7, PRF = 484)
def _np_trapz(y, x=None, dx=1.0):
    """Trapezoidal integration that works on both NumPy 1.x (trapz)
    and NumPy 2.x (trapezoid); QGIS ships either depending on release."""
    f = getattr(np, 'trapz', None) or getattr(np, 'trapezoid')
    return f(y, x) if x is not None else f(y, dx=dx)


_SC_STD  = _shape_const_for_m(3.7)
_PRF_STD = 484.0


def gamma_shape_for_prf(prf, m_lo=1e-3, m_hi=20.0):
    """
    Solve for the Gamma-UH shape parameter m that produces the supplied
    SCS peak rate factor (PRF). PRF=484 → m≈3.7 (standard SCS UH).
    Bisection over a wide bracket (no scipy dependency).
    """
    target = (float(prf) / _PRF_STD) * _SC_STD
    lo, hi = m_lo, m_hi
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if _shape_const_for_m(mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def gamma_uh_ordinates(prf, area_km2, tp_min, dt_min=1.0,
                       tail_factor=20.0):
    """
    Build a Gamma-function SCS Unit Hydrograph at dt_min resolution.
    Returns ordinates in m³/s per millimetre of unit excess.

    The ordinate at relative time τ = t / Tp is
        q(τ)/qp = τ^m · exp(m · (1 − τ))
    The peak qp (m³/s per mm of excess) is set so that the discrete
    integral ∫q dt equals the unit volume = 1000 · A km² · 1 mm = 1000·A m³.
    This guarantees rigorous mass conservation regardless of PRF.

    Parameters
    ----------
    prf : float
        Peak rate factor (100-600 typical; 484 = standard SCS).
    area_km2 : float
        Catchment area in km².
    tp_min : float
        Time to peak in minutes.
    dt_min : float
        Time step of the resulting UH (default 1 minute).
    tail_factor : float
        UH duration in multiples of Tp; the tail is scaled inversely
        with PRF so low-PRF UHs (very flat) get a longer tail.

    Returns
    -------
    uh : 1-D numpy array
        UH ordinates in m³/s per mm.
    t  : 1-D numpy array
        Times in minutes.
    """
    if area_km2 <= 0 or tp_min <= 0:
        return np.zeros(1, dtype=float), np.zeros(1, dtype=float)

    m = gamma_shape_for_prf(prf)
    # Stretch the tail proportional to (484 / PRF) so flat UHs have a
    # long enough recession; cap at a generous upper bound.
    stretch    = max(1.0, _PRF_STD / max(float(prf), 1.0))
    n_steps    = int(round(tail_factor * tp_min * stretch / dt_min)) + 1
    n_steps    = max(n_steps, int(round(60 * 24 / dt_min)))   # ≥ 24 h
    t          = np.arange(n_steps, dtype=float) * dt_min
    tau        = t / tp_min
    with np.errstate(invalid='ignore'):
        ratio = np.where(t > 0,
                         (tau ** m) * np.exp(m * (1.0 - tau)),
                         0.0)
    # Calibrate so that mass is exactly 1000·A m³ for a 1-mm pulse.
    vol_target_m3 = 1000.0 * float(area_km2)            # 1 mm × A km² = 1000 A m³
    integral      = ratio.sum() * (dt_min * 60.0)       # dt is in minutes ⇒ s
    if integral <= 0:
        return np.zeros_like(t), t
    qp_per_mm = vol_target_m3 / integral
    return qp_per_mm * ratio, t



def hms_uh_ordinates(prf, area_km2, tlag_min, dt_min=1.0):
    """SCS Unit Hydrograph matching HEC-HMS output.

    v21 (HEC-HMS parity): instead of the analytic Gamma shape, use the
    301-point dimensionless curve (t/Tp = 0..6, q/qp) obtained by
    deconvolving HEC-HMS's own PRECIP-EXCESS/FLOW results, together
    with the HEC-HMS time-to-peak convention

        Tp = dt/2 + Tlag

    and the standard metric peak  Up = 0.208 * (PRF/484) * A / Tp[hr]
    (m3/s per mm of excess). Ordinates are then rescaled so the
    discrete UH volume equals exactly 1000*A m3 per mm of excess —
    the same unit-volume normalisation HEC-HMS applies — preserving
    rigorous mass conservation.

    Falls back to the Gamma builder when the dimensionless tables are
    unavailable (console use without scs_data.py) or when the inputs
    are degenerate.
    """
    if (_DIM_UH_TR is None or _DIM_UH_QR is None
            or area_km2 <= 0 or tlag_min <= 0):
        return gamma_uh_ordinates(prf, area_km2, tlag_min, dt_min=dt_min)

    tp_min = 0.5 * float(dt_min) + float(tlag_min)     # HEC-HMS convention
    tp_hr  = tp_min / 60.0

    tr = np.asarray(_DIM_UH_TR, dtype=float)
    qr = np.asarray(_DIM_UH_QR, dtype=float)

    # ---- v1.0: make the PRF change the SHAPE, not just a scale -------
    # The peak rate factor is defined by
    #       Up = 0.208 * (PRF/484) * A / Tp
    # with the unit hydrograph always enclosing unit volume. Scaling
    # the ordinates by PRF/484 and then normalising the volume (as the
    # earlier code did) cancels exactly, leaving PRF with no effect at
    # all. The peak can only change if the TIME BASE changes with it,
    # so the recession limb is stretched by a factor s while the
    # rising limb and Tp = dt/2 + Tlag stay put (PRF changes how long
    # the hydrograph takes to recede, not when it peaks):
    #
    #       q(x) = q_r(x)                     for x <= 1
    #       q(x) = q_r(1 + (x - 1)/s)         for x >  1     (x = t/Tp)
    #
    # Unit volume then gives  peak ∝ 1 / (A_rise + s·A_rec), so
    # requiring peak/peak_484 = PRF/484 fixes s in closed form. At
    # PRF = 484 this yields s = 1 and the curve is bit-for-bit the
    # HEC-HMS-matching one, preserving the verified parity.
    _prf = max(1.0, float(prf))
    a_rise = (float(_np_trapz(qr[tr <= 1.0], tr[tr <= 1.0]))
              if np.count_nonzero(tr <= 1.0) > 1 else 0.0)
    _m = tr >= 1.0
    a_rec = (float(_np_trapz(qr[_m], tr[_m]))
             if np.count_nonzero(_m) > 1 else 0.0)
    if a_rec > 0:
        s = ((_PRF_STD / _prf) * (a_rise + a_rec) - a_rise) / a_rec
        s = min(20.0, max(0.05, s))
    else:
        s = 1.0

    t_end   = (1.0 + s * (float(tr[-1]) - 1.0)) * tp_min
    n_steps = int(np.ceil(t_end / dt_min)) + 1
    t  = np.arange(n_steps, dtype=float) * dt_min
    x  = t / tp_min
    x_eff = np.where(x <= 1.0, x, 1.0 + (x - 1.0) / s)
    up = 0.208 * (_prf / _PRF_STD) * float(area_km2) / tp_hr
    uh = up * np.interp(x_eff, tr, qr, left=0.0, right=0.0)

    # exact unit-volume normalisation: 1 mm over A km2 = 1000*A m3
    vol = float(np.sum(uh)) * float(dt_min) * 60.0
    if vol > 0:
        uh *= (1000.0 * float(area_km2)) / vol
    return uh, t

def compute_design_hydrograph(p_total_mm, cn, area_km2, tp_min, prf,
                              storm_duration_hr=24.0,
                              dt_min=1.0,
                              storm_distribution='SCS_TYPE_II',
                              total_minutes=None,
                              scs_lambda=0.2,
                              custom_hyetograph=None,
                              tail_factor=20.0,
                              max_tail_hr=96.0):
    """
    Compute a design-storm hydrograph for a single catchment / return
    period via:
        SCS Type II rainfall  →  SCS-CN losses  →  Gamma SCS UH

    Parameters
    ----------
    p_total_mm : float
        Total storm depth in mm.
    cn : float
        Curve number (50-95 typical).
    area_km2 : float
        Catchment area in km².
    tp_min : float
        Time to peak (= Tlag) in minutes for the SCS UH. Selected
        upstream by TC_METHOD.
    prf : float
        SCS peak rate factor (100-600).
    storm_duration_hr : float
        Storm duration (default 24 hours).
    dt_min : float
        Time step in minutes (default 1).
    storm_distribution : str
        Only 'SCS_TYPE_II' supported in v14.
    total_minutes : int or None
        Desired total length of the time series. If None, defaults to
        storm_duration_hr · 60 + UH tail length.

    Returns
    -------
    out : dict with arrays
        time_min   times since storm start, length N
        P_mm       incremental precipitation,        length N
        Loss_mm    incremental loss (P − Excess),    length N
        Excess_mm  incremental excess,               length N
        Q_m3s      total flow,                       length N
        peak_Q     scalar, m³/s
        peak_t     scalar, minutes from storm start
        R_total    scalar, total runoff depth in mm
    """
    if p_total_mm is None or p_total_mm <= 0:
        return None
    if cn is None or cn <= 0 or cn > 100:
        return None
    if tp_min is None or tp_min <= 0:
        return None
    if area_km2 <= 0:
        return None

    # 1. Rainfall hyetograph (incremental, mm per dt).
    # v19.1: a CUSTOM hyetograph fixes the storm duration from its own
    # time axis (STORM_DURATION_HR is ignored in that case).
    if (str(storm_distribution or '').upper().strip() == 'CUSTOM'
            and custom_hyetograph is not None):
        storm_duration_hr = float(custom_hyetograph['duration_hr'])
    storm_min = int(round(storm_duration_hr * 60.0))
    # Storm + sufficient post-storm window to capture the recession.
    if total_minutes is None:
        total_minutes = int(round(storm_min
                                  + float(tail_factor) * tp_min *
                                  max(1.0, _PRF_STD / max(float(prf), 1.0))))
        # Cap to keep runtimes sane (user-editable, v1.0)
        total_minutes = min(total_minutes,
                            int(round(storm_min
                                      + 60.0 * float(max_tail_hr))))
    n_steps_storm = int(round(storm_min / dt_min)) + 1
    t_storm   = np.arange(n_steps_storm, dtype=float) * dt_min
    F_storm   = design_storm_cum_fraction(t_storm,
                                          distribution=storm_distribution,
                                          duration_hr=storm_duration_hr,
                                          custom_curve=custom_hyetograph)
    cum_P     = F_storm * float(p_total_mm)         # mm, length n_steps_storm
    incr_P    = np.diff(cum_P, prepend=0.0)          # mm/step

    # 2. SCS-CN losses applied to the cumulative rainfall.
    lam   = min(0.5, max(0.01, float(scs_lambda)
                         if scs_lambda is not None else 0.2))
    s_mm  = 25.4 * (1000.0 / cn - 10.0)
    ia_mm = lam * s_mm
    cum_R = np.where(cum_P > ia_mm,
                     (cum_P - ia_mm) ** 2
                     / np.maximum(cum_P - ia_mm + s_mm, 1e-9),
                     0.0)
    incr_R = np.diff(cum_R, prepend=0.0)
    incr_R = np.clip(incr_R, 0.0, None)
    incr_loss = incr_P - incr_R                     # mm/step

    # 3. SCS Unit Hydrograph. v21: HEC-HMS-matching dimensionless
    #    curve with Tp = dt/2 + Tlag (the tp_min argument carries the
    #    lag, as at every call site); Gamma fallback without scs_data.
    uh, _t_uh = hms_uh_ordinates(prf, area_km2, tp_min, dt_min=dt_min)

    # 4. Convolve excess hyetograph with the UH.
    Q = np.convolve(incr_R, uh, mode='full')

    # 5. Pad / trim every column to a common length.
    n_total = int(round(total_minutes / dt_min)) + 1
    n_total = max(n_total, len(Q))
    def _pad(a, n):
        if len(a) >= n:
            return a[:n]
        out = np.zeros(n, dtype=float)
        out[:len(a)] = a
        return out
    incr_P    = _pad(incr_P,    n_total)
    incr_R    = _pad(incr_R,    n_total)
    incr_loss = _pad(incr_loss, n_total)
    Q         = _pad(Q,         n_total)
    time_min  = np.arange(n_total, dtype=float) * dt_min

    return dict(
        time_min  = time_min,
        P_mm      = incr_P,
        Loss_mm   = incr_loss,
        Excess_mm = incr_R,
        Q_m3s     = Q,
        peak_Q    = float(Q.max()),
        peak_t    = float(time_min[int(Q.argmax())]),
        R_total   = float(cum_R[-1]) if len(cum_R) else 0.0,
    )


def select_tlag_for_method(method, ch, tlag_nrcs, min_tc_min=0.0):
    """
    Pick the (Tlag, Tc, label) triple driven into the peak-flow /
    hydrograph engine according to CONFIG['TC_METHOD'].

    Falls back to KIRPICH_HOLISTIC (with the returned label suffixed
    by ' [fallback]') when NRCS is requested but no usable CN-based
    Tlag is available.
    """
    method = (method or 'KIRPICH_HOLISTIC').upper().strip()
    holistic_tlag = ch.get('TLAG_KIRPH', 0.0)
    holistic_tc   = ch.get('TC_KIRPH',   0.0)
    if method == 'KIRPICH_HOLISTIC':
        out = (holistic_tlag, holistic_tc, 'KIRPICH_HOLISTIC')
    elif method == 'KIRPICH_10_85':
        out = (ch.get('TLAG', 0.0), ch.get('TC_TOTAL', 0.0),
               'KIRPICH_10_85')
    elif method == 'NRCS':
        if tlag_nrcs and tlag_nrcs > 0:
            tc_n = tlag_nrcs / 0.6 if tlag_nrcs > 0 else 0.0
            out = (float(tlag_nrcs), float(tc_n), 'NRCS')
        else:
            out = (holistic_tlag, holistic_tc,
                   'KIRPICH_HOLISTIC [fallback]')
    else:
        out = (holistic_tlag, holistic_tc,
               f'KIRPICH_HOLISTIC [unknown: {method}]')
    return _apply_tc_floor(out, min_tc_min)


def _apply_tc_floor(triple, min_tc_min):
    """v1.0 — raise (Tlag, Tc) to the user minimum Tc when set.

    The floor is expressed on Tc (the familiar quantity); the matching
    Tlag = 0.6 x Tc is used with it, keeping the pair consistent. The
    label is suffixed so the attribute tables and the log show clearly
    where the minimum was binding."""
    tlag, tc, label = triple
    try:
        floor = float(min_tc_min or 0.0)
    except (TypeError, ValueError):
        floor = 0.0
    if floor > 0 and float(tc or 0.0) < floor:
        return 0.6 * floor, floor, f'{label} [Tc floored to {floor:g} min]'
    return tlag, tc, label


def write_hydrograph_csv(path, hyd):
    """Write a single hydrograph dict to CSV with header line."""
    import csv as _csv
    n = len(hyd['time_min'])
    with open(path, 'w', newline='') as fh:
        w = _csv.writer(fh)
        w.writerow(['time_min', 'P_mm', 'Loss_mm', 'Excess_mm', 'Q_m3s'])
        for i in range(n):
            w.writerow([
                f"{hyd['time_min'][i]:g}",
                f"{hyd['P_mm'][i]:.5f}",
                f"{hyd['Loss_mm'][i]:.5f}",
                f"{hyd['Excess_mm'][i]:.5f}",
                f"{hyd['Q_m3s'][i]:.5f}",
            ])


def write_q_only_csv(path, time_min, q_m3s, qin_m3s=None):
    """Write a routed-reach or junction-combined hydrograph to CSV.

    Sub-basin local hydrographs use write_hydrograph_csv() because they
    have full P/Loss/Excess columns. Reaches and junctions only have
    Q values (and optionally an upstream inflow Q for reaches), so a
    leaner schema is used here:
        time_min, Q_in_m3s [optional], Q_out_m3s
    """
    import csv as _csv
    n = len(time_min)
    has_in = (qin_m3s is not None)
    with open(path, 'w', newline='') as fh:
        w = _csv.writer(fh)
        if has_in:
            w.writerow(['time_min', 'Q_in_m3s', 'Q_out_m3s'])
        else:
            w.writerow(['time_min', 'Q_m3s'])
        for i in range(n):
            row = [f"{time_min[i]:g}"]
            if has_in:
                row.append(f"{float(qin_m3s[i]):.5f}")
            row.append(f"{float(q_m3s[i]):.5f}")
            w.writerow(row)


def write_reservoir_csv(path, time_min, qin, qout, storage, stage=None):
    """v24 — AR-Awdiya-style reservoir routing table:
        time_min, Inflow_m3s, Outflow_m3s, Storage_m3 [, Stage_m]
    (instantaneous pool storage, not cumulative)."""
    import csv as _csv
    has_stage = stage is not None
    with open(path, 'w', newline='') as fh:
        w = _csv.writer(fh)
        head = ['time_min', 'Inflow_m3s', 'Outflow_m3s', 'Storage_m3']
        if has_stage:
            head.append('Stage_m')
        w.writerow(head)
        for i in range(len(time_min)):
            row = [f"{time_min[i]:g}",
                   f"{float(qin[i]):.5f}",
                   f"{float(qout[i]):.5f}",
                   f"{float(storage[i]):.2f}"]
            if has_stage:
                row.append(f"{float(stage[i]):.3f}")
            w.writerow(row)


# =============================================================================
# 7b-bis. MODELING OUTPUTS DIRECTORY + HTML VISUALISATIONS (v18)
# =============================================================================
# v18 reorganises hydrograph outputs into one folder per element type
# (Catchments / SubBasins / Reaches / Junctions), each with a CSVs/
# subfolder for raw data and an HTMLs/ subfolder for interactive
# Plotly-based visualisations. The HTMLs are self-contained: they
# pull plotly.min.js from the CDN at first load and then work offline.

# Element-type → directory-name mapping. Used everywhere the tool
# decides where to write a CSV or HTML.
MODELING_DIR_MAP = {
    'CATCHMENT' : 'Catchments',
    'SUBBASIN'  : 'SubBasins',
    'REACH'     : 'Reaches',
    'JUNCTION'  : 'Junctions',
}


def make_modeling_outputs_dirs(out_dir):
    """Create the v18 ``Modeling_Outputs/`` tree under ``out_dir``.

    Returns a dict keyed by element-type code with ``csv`` and
    ``html`` subdirectory paths:

        {
          'CATCHMENT': {'csv': '.../Catchments/CSVs',
                        'html': '.../Catchments/HTMLs'},
          'SUBBASIN' : {'csv': ..., 'html': ...},
          ...
        }

    Idempotent — safe to call repeatedly. Returns None if ``out_dir``
    is None.
    """
    if not out_dir:
        return None
    root = os.path.join(out_dir, 'Modeling_Outputs')
    paths = {}
    for code, name in MODELING_DIR_MAP.items():
        d_csv  = os.path.join(root, name, 'CSVs')
        d_html = os.path.join(root, name, 'HTMLs')
        os.makedirs(d_csv,  exist_ok=True)
        os.makedirs(d_html, exist_ok=True)
        paths[code] = {'csv': d_csv, 'html': d_html, 'root': os.path.join(root, name)}
    paths['ROOT'] = root
    return paths


# ---- HTML template -------------------------------------------------------
# A minimal, robust HTML wrapper. Plotly is loaded from CDN; if it
# fails to load (e.g. the user opens the HTML offline before the CDN
# was ever cached by their browser), the page still shows the title
# and summary card, with a friendly notice instead of a blank chart.
_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>{title}</title>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<style>
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; padding: 24px;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI",
                 Roboto, Helvetica, Arial, sans-serif;
    color: #1f2933; background: #f5f7fa;
  }}
  h1 {{
    margin: 0 0 4px 0; font-size: 22px; font-weight: 600;
    color: #102a43;
  }}
  .subtitle {{
    color: #486581; font-size: 14px; margin-bottom: 20px;
  }}
  .summary {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
    gap: 12px; margin-bottom: 24px;
  }}
  .card {{
    background: white; border: 1px solid #d9e2ec;
    border-radius: 8px; padding: 14px 16px;
    box-shadow: 0 1px 2px rgba(0,0,0,0.04);
  }}
  .card .label {{
    text-transform: uppercase; font-size: 11px;
    letter-spacing: 0.06em; color: #627d98; font-weight: 600;
  }}
  .card .value {{
    font-size: 20px; font-weight: 600; color: #102a43;
    margin-top: 4px;
  }}
  .card .unit {{
    font-size: 12px; color: #829ab1; font-weight: 400;
    margin-left: 4px;
  }}
  #chart {{
    background: white; border: 1px solid #d9e2ec;
    border-radius: 8px; padding: 12px; min-height: 520px;
  }}
  #fallback {{
    display: none; padding: 24px; text-align: center;
    color: #829ab1; font-size: 14px;
  }}
  footer {{
    margin-top: 16px; color: #829ab1; font-size: 11px;
    text-align: right;
  }}
</style>
</head>
<body>
  <h1>{title}</h1>
  <div class="subtitle">{subtitle}</div>
  <div class="summary">{summary_cards}</div>
  <div id="chart"></div>
  <div id="fallback">
    Plotly could not be loaded from the CDN.
    The CSV next to this file contains the raw data.
  </div>
  <footer>Generated by AR-Tarseem v1.0.0</footer>

<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"
        onerror="document.getElementById('chart').style.display='none';
                 document.getElementById('fallback').style.display='block';">
</script>
<script>
  var DATA = {data_json};
  var LAYOUT = {layout_json};
  var CONFIG = {{ responsive: true, displaylogo: false,
                 modeBarButtonsToRemove: ['lasso2d', 'select2d'] }};
  if (typeof Plotly !== 'undefined') {{
    Plotly.newPlot('chart', DATA, LAYOUT, CONFIG);
  }} else {{
    document.getElementById('chart').style.display = 'none';
    document.getElementById('fallback').style.display = 'block';
  }}
</script>
</body>
</html>
"""


def _summary_card(label, value, unit=''):
    """Tiny HTML snippet for a stat card."""
    unit_html = (f'<span class="unit">{unit}</span>' if unit else '')
    return (f'<div class="card"><div class="label">{label}</div>'
            f'<div class="value">{value}{unit_html}</div></div>')


def _hours_axis(time_min):
    """Convert minutes-from-start to hours for nicer x-axis ticks."""
    return [t / 60.0 for t in time_min]


def write_hydrograph_html(path, kind, title, subtitle, hyd_data,
                          extra_summary=None):
    """Write a self-contained interactive HTML for one element.

    Parameters
    ----------
    path : str
        Destination .html file.
    kind : {'CATCHMENT', 'SUBBASIN', 'REACH', 'JUNCTION'}
        Which template to use.
    title : str
        Bold heading shown on the page.
    subtitle : str
        Smaller line under the title (e.g. "100-yr design storm").
    hyd_data : dict
        Time-series data, structure depends on ``kind``:
            CATCHMENT / SUBBASIN:
                time_min, P_mm, Loss_mm, Excess_mm, Q_m3s
            REACH:
                time_min, Q_in_m3s, Q_out_m3s
            JUNCTION:
                time_min, Q_m3s
    extra_summary : list of (label, value, unit) tuples, optional
        Additional cards to show above the chart (peak Q is added
        automatically).
    """
    import json

    # Build summary cards
    cards = []
    time_min = list(hyd_data.get('time_min', []))
    if not time_min:
        cards.append(_summary_card('No data', '—'))
    else:
        # Always show peak Q
        if kind == 'REACH':
            q_in  = hyd_data.get('Q_in_m3s', [])
            q_out = hyd_data.get('Q_out_m3s', [])
            qp_in  = max(q_in)  if q_in  else 0.0
            qp_out = max(q_out) if q_out else 0.0
            atten  = (100.0 * (1 - qp_out / qp_in)) if qp_in > 0 else 0.0
            cards.append(_summary_card('Peak inflow',  f'{qp_in:.2f}',  'm³/s'))
            cards.append(_summary_card('Peak outflow', f'{qp_out:.2f}', 'm³/s'))
            cards.append(_summary_card('Attenuation',  f'{atten:.1f}',  '%'))
        else:
            q = hyd_data.get('Q_m3s', [])
            qp = max(q) if q else 0.0
            ip = q.index(qp) if (q and qp > 0) else 0
            cards.append(_summary_card('Peak Q', f'{qp:.2f}', 'm³/s'))
            cards.append(_summary_card('Time to peak',
                                       f'{time_min[ip]/60.0:.2f}', 'h'))
            if 'P_mm' in hyd_data and hyd_data['P_mm']:
                p_total = sum(hyd_data['P_mm'])
                cards.append(_summary_card('Total rainfall',
                                           f'{p_total:.1f}', 'mm'))
            if 'Excess_mm' in hyd_data and hyd_data['Excess_mm']:
                r_total = sum(hyd_data['Excess_mm'])
                cards.append(_summary_card('Runoff depth',
                                           f'{r_total:.1f}', 'mm'))
            # Volume = ∫ Q dt = sum(Q) × dt(s)
            if q:
                # infer dt
                dt_min = (time_min[1] - time_min[0]) if len(time_min) > 1 else 1.0
                vol = sum(q) * dt_min * 60.0
                cards.append(_summary_card('Total volume',
                                           f'{vol/1000.0:.1f}', '×10³ m³'))
    if extra_summary:
        for lab, val, unit in extra_summary:
            cards.append(_summary_card(lab, val, unit))
    summary_html = ''.join(cards)

    # Build chart traces + layout per kind
    t_h = _hours_axis(time_min)

    def safelist(key):
        return list(hyd_data.get(key, []))

    if kind in ('CATCHMENT', 'SUBBASIN'):
        # Three stacked panels:
        #   Row 1: P / Loss / Excess (mm/min)
        #   Row 2: Cumulative depth  (mm)
        #   Row 3: Q (m³/s)
        P    = safelist('P_mm')
        Loss = safelist('Loss_mm')
        Exc  = safelist('Excess_mm')
        Q    = safelist('Q_m3s')

        # cumulative
        def cumsum(xs):
            s = 0.0; out = []
            for x in xs:
                s += float(x); out.append(s)
            return out
        Pcum    = cumsum(P)
        Losscum = cumsum(Loss)
        Excum   = cumsum(Exc)

        data = [
            {'type': 'scatter', 'x': t_h, 'y': P,
             'mode': 'lines', 'name': 'Rainfall',
             'line': {'color': '#1f78b4', 'width': 1.0},
             'fill': 'tozeroy', 'fillcolor': 'rgba(31,120,180,0.18)',
             'xaxis': 'x', 'yaxis': 'y'},
            {'type': 'scatter', 'x': t_h, 'y': Loss,
             'mode': 'lines', 'name': 'Loss',
             'line': {'color': '#b2182b', 'width': 1.0, 'dash': 'dot'},
             'xaxis': 'x', 'yaxis': 'y'},
            {'type': 'scatter', 'x': t_h, 'y': Exc,
             'mode': 'lines', 'name': 'Excess',
             'line': {'color': '#4daf4a', 'width': 1.4},
             'xaxis': 'x', 'yaxis': 'y'},

            {'type': 'scatter', 'x': t_h, 'y': Pcum,
             'mode': 'lines', 'name': 'Cumulative P',
             'line': {'color': '#1f78b4', 'width': 1.6},
             'xaxis': 'x2', 'yaxis': 'y2', 'showlegend': False},
            {'type': 'scatter', 'x': t_h, 'y': Losscum,
             'mode': 'lines', 'name': 'Cumulative Loss',
             'line': {'color': '#b2182b', 'width': 1.6, 'dash': 'dot'},
             'xaxis': 'x2', 'yaxis': 'y2', 'showlegend': False},
            {'type': 'scatter', 'x': t_h, 'y': Excum,
             'mode': 'lines', 'name': 'Cumulative Excess',
             'line': {'color': '#4daf4a', 'width': 1.6},
             'xaxis': 'x2', 'yaxis': 'y2', 'showlegend': False},

            {'type': 'scatter', 'x': t_h, 'y': Q,
             'mode': 'lines', 'name': 'Discharge',
             'line': {'color': '#222', 'width': 2.0},
             'fill': 'tozeroy', 'fillcolor': 'rgba(34,34,34,0.10)',
             'xaxis': 'x3', 'yaxis': 'y3'},
        ]
        layout = {
            'grid': {'rows': 3, 'columns': 1, 'pattern': 'independent',
                     'roworder': 'top to bottom'},
            'xaxis':  {'title': '', 'showticklabels': False, 'domain': [0, 1]},
            'xaxis2': {'title': '', 'showticklabels': False, 'domain': [0, 1]},
            'xaxis3': {'title': 'Time (hours from storm start)',
                       'domain': [0, 1]},
            'yaxis':  {'title': 'Rate (mm/min)',  'rangemode': 'tozero'},
            'yaxis2': {'title': 'Depth (mm)',     'rangemode': 'tozero'},
            'yaxis3': {'title': 'Discharge (m³/s)', 'rangemode': 'tozero'},
            'height': 760, 'margin': {'l': 70, 'r': 30, 't': 30, 'b': 60},
            'hovermode': 'x unified',
            'legend': {'orientation': 'h', 'y': 1.06, 'x': 0.5,
                       'xanchor': 'center'},
            'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
        }

    elif kind == 'REACH':
        # Inflow vs routed outflow overlay
        q_in  = safelist('Q_in_m3s')
        q_out = safelist('Q_out_m3s')
        data = [
            {'type': 'scatter', 'x': t_h, 'y': q_in,
             'mode': 'lines', 'name': 'Inflow (upstream)',
             'line': {'color': '#1f78b4', 'width': 2.0, 'dash': 'dot'}},
            {'type': 'scatter', 'x': t_h, 'y': q_out,
             'mode': 'lines', 'name': 'Routed outflow (downstream)',
             'line': {'color': '#e41a1c', 'width': 2.4},
             'fill': 'tonexty', 'fillcolor': 'rgba(228,26,28,0.10)'},
        ]
        layout = {
            'xaxis': {'title': 'Time (hours from storm start)'},
            'yaxis': {'title': 'Discharge (m³/s)', 'rangemode': 'tozero'},
            'height': 540, 'margin': {'l': 70, 'r': 30, 't': 30, 'b': 60},
            'hovermode': 'x unified',
            'legend': {'orientation': 'h', 'y': 1.08, 'x': 0.5,
                       'xanchor': 'center'},
            'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
        }

    else:   # JUNCTION
        Q = safelist('Q_m3s')
        data = [
            {'type': 'scatter', 'x': t_h, 'y': Q,
             'mode': 'lines', 'name': 'Combined discharge',
             'line': {'color': '#33a02c', 'width': 2.4},
             'fill': 'tozeroy', 'fillcolor': 'rgba(51,160,44,0.12)'},
        ]
        layout = {
            'xaxis': {'title': 'Time (hours from storm start)'},
            'yaxis': {'title': 'Discharge (m³/s)', 'rangemode': 'tozero'},
            'height': 540, 'margin': {'l': 70, 'r': 30, 't': 30, 'b': 60},
            'hovermode': 'x unified',
            'showlegend': False,
            'plot_bgcolor': 'white', 'paper_bgcolor': 'white',
        }

    html = _HTML_TEMPLATE.format(
        title         = title,
        subtitle      = subtitle,
        summary_cards = summary_html,
        data_json     = json.dumps(data),
        layout_json   = json.dumps(layout),
    )
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(html)


# =============================================================================
# 7c. REACH ROUTING + TOPOLOGICAL JUNCTION COMBINATION (v17)
# =============================================================================
def lag_route(inflow, lag_min, dt_min=1.0):
    """
    Pure-translation routing: shift the inflow hydrograph forward by
    ``lag_min`` minutes. No attenuation.

    v1.0 (HEC-HMS parity fix) — the shift is applied EXACTLY, with
    linear interpolation between ordinates when the lag is not a whole
    number of time steps:

        O(t) = I(t − lag)        (I linearly interpolated)

    which is precisely what HEC-HMS "Specified Lag" routing does.
    Previously the lag was rounded to the nearest whole step
    (±dt/2 per reach); across a network of nested reaches those
    rounding errors accumulate differently along each tributary path,
    mutually misaligning the superposed peaks at downstream junctions
    and systematically lowering the combined peak by a fraction of a
    percent relative to HEC-HMS. Exact fractional translation removes
    that artifact; integer lags reproduce the old behaviour
    bit-for-bit.

    Before water arrives, the outflow holds the inflow's initial value
    (hydrographs start at zero, so in practice the leading samples are
    zero); the tail that would arrive after the simulation window ends
    is dropped, as before.

    Parameters
    ----------
    inflow : 1-D numpy array
        Inflow hydrograph at the upstream end of the reach (m³/s).
    lag_min : float
        Travel time on the reach in minutes.
    dt_min : float
        Time step of the inflow series (default 1 minute).

    Returns
    -------
    outflow : 1-D numpy array, same length as inflow.
    """
    inflow = np.asarray(inflow, dtype=float)
    if lag_min is None or lag_min <= 0:
        return inflow.copy()
    n = len(inflow)
    lag = float(lag_min)
    dt  = float(dt_min)
    if lag / dt >= n:
        return np.full(n, inflow[0], dtype=float)
    t = np.arange(n, dtype=float) * dt
    # exact translation: sample the inflow at (t − lag), linear
    # interpolation between ordinates, clamped to I[0] before arrival
    return np.interp(t - lag, t, inflow,
                     left=float(inflow[0]), right=float(inflow[-1]))


def muskingum_subreaches(K_min, X, dt_min=1.0):
    """v1.0 — number of Muskingum subreaches N so that every subreach
    satisfies the classical accuracy/positivity bound

        2·(K/N)·X  ≤  Δt          (⇒ C0 ≥ 0)

    which, together with X ≤ 0.5, also guarantees Δt ≤ 2·(K/N)·(1−X)
    (⇒ C2 ≥ 0). This mirrors how HEC-HMS treats long reaches: the
    reach is solved as a cascade of subreaches instead of one step,
    which removes the negative-coefficient numerical dispersion that
    otherwise lowers routed peaks when Δt ≪ 2·K·X."""
    if K_min is None or K_min <= 0 or X is None or X <= 0:
        return 1
    return max(1, int(math.ceil(2.0 * float(K_min) * float(X)
                                / max(float(dt_min), 1e-9))))


def muskingum_route(inflow, K_min, X, dt_min=1.0, n_sub=None):
    """
    Classical Muskingum reach routing, solved as a subreach cascade.

    Storage:        S = K · [X · I + (1-X) · O]
    Continuity:     dS/dt = I - O
    →  O[t+1] = C0·I[t+1] + C1·I[t] + C2·O[t]
    where the standard discretized coefficients are:
        denom = K - K·X + Δt/2
        C0    = (-K·X + Δt/2) / denom
        C1    = ( K·X + Δt/2) / denom
        C2    = ( K - K·X - Δt/2) / denom
    and C0 + C1 + C2 = 1 exactly (mass-balance preserving).

    v1.0 — the reach is subdivided into N = ceil(2·K·X / Δt) subreaches
    of K/N each (see muskingum_subreaches), so every subreach satisfies
    2·K'·X ≤ Δt ≤ 2·K'·(1−X) and all three coefficients are
    non-negative. The single-step scheme previously used at
    Δt ≪ 2·K·X produced a negative C0 whose numerical dispersion
    systematically lowered routed peaks by a fraction of a percent;
    the cascade removes that artifact and makes the routed hydrograph
    insensitive to the computational time step. The same N is written
    to the HEC-HMS .basin export ('Subreaches:'), keeping the
    in-plugin results and the HEC-HMS re-computation aligned.

    Parameters
    ----------
    inflow : 1-D numpy array
        Reach inflow time series (m³/s) at uniform Δt.
    K_min : float
        Muskingum K for the WHOLE reach, in MINUTES (set per reach to
        that reach's Tc).
    X : float
        Muskingum weighting in [0, 0.5]. 0.0 = linear reservoir
        (full attenuation), 0.5 = pure translation (no attenuation).
        Typical natural-reach value: 0.2.
    dt_min : float
        Time step in minutes (default 1).

    Returns
    -------
    outflow : 1-D numpy array, same length as inflow.
    """
    inflow = np.asarray(inflow, dtype=float)
    n = len(inflow)
    if K_min is None or K_min <= 0:
        return inflow.copy()
    K_tot = float(K_min)
    X  = float(X)
    dt = float(dt_min)
    # v1.0 — n_sub override: a positive value fixes the subreach count
    # (dialog option "Muskingum subreaches"); None/0 keeps the
    # automatic stability-based estimate.
    n_sub = (int(n_sub) if n_sub else
             muskingum_subreaches(K_tot, X, dt))
    K = K_tot / n_sub
    denom = K - K * X + 0.5 * dt
    if denom <= 0:
        # numerically degenerate (extremely small K with very large X);
        # fall back to lag routing
        return lag_route(inflow, K_min, dt_min)
    C0 = (-K * X + 0.5 * dt) / denom
    C1 = ( K * X + 0.5 * dt) / denom
    C2 = ( K - K * X - 0.5 * dt) / denom
    # ---- fast exact cascade via impulse response ---------------------
    # A single subreach is the linear filter
    #     out[t] = C0·I[t] + C1·I[t-1] + C2·out[t-1]
    # whose unit-impulse response is
    #     h[0] = C0,  h[1] = C1 + C2·C0,  h[t≥2] = C2·h[t-1]
    # The N-subreach cascade response is the N-fold convolution
    # h_N = h * h * … * h, and (with the steady-state initial
    # condition out(0) = inflow(0), exact because C0+C1+C2 = 1)
    #     out = I0 + (inflow − I0) ⊛ h_N
    # Everything runs in C-speed numpy convolutions, so the cascade
    # costs no more than the old single-step loop.
    h1 = [C0, C1 + C2 * C0]
    _t = h1[-1]
    _c2 = abs(C2)
    while abs(_t) > 1e-15 and len(h1) < n:
        _t *= C2
        h1.append(_t)
        if _c2 < 1e-12:
            break
    h1 = np.asarray(h1, dtype=float)
    hN = h1
    for _s in range(n_sub - 1):
        hN = np.convolve(hN, h1)
        if len(hN) > n:
            hN = hN[:n]
        # prune the negligible tail to keep the convolutions small
        _mx = np.max(np.abs(hN))
        if _mx > 0:
            _nz = np.nonzero(np.abs(hN) > 1e-15 * _mx)[0]
            if len(_nz):
                hN = hN[:int(_nz[-1]) + 1]
    i0 = float(inflow[0])
    out = i0 + np.convolve(inflow - i0, hN)[:n]
    return out


def _trapezoid_manning_q(y, b, z, n, s0):
    """Manning discharge (m3/s) for a trapezoidal section of depth y."""
    if y <= 0:
        return 0.0
    a = (b + z * y) * y
    p = b + 2.0 * y * math.sqrt(1.0 + z * z)
    if p <= 0 or a <= 0:
        return 0.0
    r = a / p
    return a / n * (r ** (2.0 / 3.0)) * math.sqrt(s0)


def _normal_depth(q, b, z, n, s0):
    """Normal depth (m) by bisection of the Manning relation.

    Solver: expand the upper bracket until it
    exceeds the target discharge (cap 80 m), then 45 bisection steps.
    """
    lo, hi = 1e-4, 0.5
    while _trapezoid_manning_q(hi, b, z, n, s0) < q and hi < 80.0:
        hi *= 2.0
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        if _trapezoid_manning_q(mid, b, z, n, s0) < q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def muskingum_cunge_route(inflow, dt_min, length_m, slope,
                          manning_n=0.035, bottom_width_m=5.0,
                          side_slope=2.0):
    """Variable-parameter Muskingum-Cunge routing .

    The reach is split into N = round(L / (c0*dt)) sub-reaches
    (1 <= N <= 40) using a reference celerity at 2/3 of the inflow
    peak. Within each sub-reach, at every time step, K = dx/c and
    X = 0.5*(1 - Q/(T*S0*c*dx)) are recomputed from the local
    three-point reference discharge, X clipped to [0, 0.5], and the
    standard Muskingum coefficients applied.
    """
    inflow = np.asarray(inflow, dtype=float)
    L  = max(1.0,   float(length_m or 0.0))
    s0 = max(1e-5,  float(slope or 0.0))
    n  = max(0.008, float(manning_n))
    b  = max(0.05,  float(bottom_width_m))
    z  = max(0.0,   float(side_slope))
    dt_hr  = float(dt_min) / 60.0
    dt_sec = dt_hr * 3600.0

    qp = float(inflow.max()) if inflow.size else 0.0
    if qp <= 0:
        return inflow.copy()
    q_ref0 = max(0.05, 0.667 * qp)
    y0 = _normal_depth(q_ref0, b, z, n, s0)
    a0 = (b + z * y0) * y0
    c0 = max(0.05, 5.0 / 3.0 * q_ref0 / max(a0, 1e-9))
    n_seg = int(min(40, max(1, round(L / (c0 * dt_sec)))))
    dx = L / n_seg

    cur = inflow
    clipped = 0
    for _seg in range(n_seg):
        out = np.zeros_like(cur)
        out[0] = cur[0]
        for i in range(1, len(cur)):
            q_r = max(0.001, (cur[i] + cur[i - 1] + out[i - 1]) / 3.0)
            y = _normal_depth(q_r, b, z, n, s0)
            a = (b + z * y) * y
            tw = b + 2.0 * z * y
            c = max(0.05, 5.0 / 3.0 * q_r / max(a, 1e-9))
            K = dx / c / 3600.0                      # hours
            X = 0.5 * (1.0 - q_r / (tw * s0 * c * dx))
            if X < 0.0:
                X = 0.0
                clipped += 1
            elif X > 0.5:
                X = 0.5
            den = 2.0 * K * (1.0 - X) + dt_hr
            c0_ = (dt_hr - 2.0 * K * X) / den
            c1_ = (dt_hr + 2.0 * K * X) / den
            c2_ = (2.0 * K * (1.0 - X) - dt_hr) / den
            out[i] = max(0.0, c0_ * cur[i] + c1_ * cur[i - 1]
                         + c2_ * out[i - 1])
        cur = out
    if clipped > 0.3 * n_seg * inflow.size:
        print("     \u26a0 Muskingum-Cunge X heavily clipped \u2014 check "
              "channel geometry / slope or reduce HYD_TIME_STEP_MIN")
    return cur


def route_reach(inflow, tc_min, method='MUSKINGUM',
                muskingum_x=0.2, dt_min=1.0,
                reach_length_m=0.0, reach_slope=0.0,
                mc_manning_n=0.035, mc_bottom_width_m=5.0,
                mc_side_slope=2.0, muskingum_n_sub=None):
    """Dispatcher: route an inflow hydrograph through one reach.

    method ∈ {'MUSKINGUM', 'LAG', 'NONE'} — anything unrecognised falls
    back to 'NONE' (outflow == inflow) with a console warning.
    """
    inflow = np.asarray(inflow, dtype=float)
    m = (method or 'NONE').upper().strip()
    if m == 'NONE':
        return inflow.copy()
    if m == 'LAG':
        return lag_route(inflow, tc_min, dt_min=dt_min)
    if m == 'MUSKINGUM':
        return muskingum_route(inflow, tc_min, muskingum_x, dt_min=dt_min,
                               n_sub=muskingum_n_sub)
    if m == 'MUSKINGUM_CUNGE':
        return muskingum_cunge_route(
            inflow, dt_min, reach_length_m, reach_slope,
            manning_n=mc_manning_n,
            bottom_width_m=mc_bottom_width_m,
            side_slope=mc_side_slope)
    print(f"   ⚠ REACH_ROUTING_METHOD={method!r} not recognised — "
          f"using NONE (no attenuation)")
    return inflow.copy()


def select_reach_tc(rch, method='KIRPICH_HOLISTIC', min_tc_min=0.0):
    """Pick the reach Tc value (minutes) consumed by route_reach().

    Falls back to TC_KIRPH if the requested method's value is missing
    or non-positive.
    """
    method = (method or 'KIRPICH_HOLISTIC').upper().strip()
    tc_kirph = float(rch.get('TC_KIRPH') or 0.0)
    tc_total = float(rch.get('TC_TOTAL') or 0.0)
    if method == 'KIRPICH_10_85' and tc_total > 0:
        out = (tc_total, 'KIRPICH_10_85')
    elif method == 'KIRPICH_HOLISTIC' and tc_kirph > 0:
        out = (tc_kirph, 'KIRPICH_HOLISTIC')
    elif tc_kirph > 0:
        out = (tc_kirph, f'KIRPICH_HOLISTIC [fallback from {method}]')
    elif tc_total > 0:
        out = (tc_total, 'KIRPICH_10_85 [fallback]')
    else:
        out = (0.0, 'NONE [no Tc available]')
    # v1.0 — minimum reach travel time (drives Muskingum K / the Lag
    # shift). Applied last so it also lifts a missing/zero Tc.
    try:
        floor = float(min_tc_min or 0.0)
    except (TypeError, ValueError):
        floor = 0.0
    if floor > 0 and float(out[0] or 0.0) < floor:
        return floor, f'{out[1]} [floored to {floor:g} min]'
    return out


def topological_order_subbasins(sid_of, ds_sid_of):
    """
    Kahn's algorithm: produce sub-basin SUB_IDs in upstream→downstream
    order so each sub-basin is processed only AFTER all of its
    upstream contributors are finished.

    Parameters
    ----------
    sid_of : dict
        junction (r,c) tuple → SUB_ID. Already keyed in the catchment.
    ds_sid_of : dict
        SUB_ID → SUB_ID of the immediately downstream sub-basin (or 0
        if this is the outlet sub-basin = sid 1).

    Returns
    -------
    order : list of int
        SUB_IDs in valid topological order. The outlet (sid=1) is
        always the LAST element.
    """
    all_sids = sorted(set(sid_of.values()))
    indeg = {s: 0 for s in all_sids}
    for s in all_sids:
        ds = ds_sid_of.get(s, 0)
        if ds in indeg:
            indeg[ds] += 1
    # start with junctions that have no upstream contributors
    queue = [s for s in all_sids if indeg[s] == 0]
    order = []
    while queue:
        # pop smallest sid first for deterministic order
        queue.sort()
        s = queue.pop(0)
        order.append(s)
        ds = ds_sid_of.get(s, 0)
        if ds in indeg:
            indeg[ds] -= 1
            if indeg[ds] == 0:
                queue.append(ds)
    if len(order) != len(all_sids):
        # cycle — shouldn't happen for a true DAG drainage network
        missing = [s for s in all_sids if s not in order]
        print(f"   ⚠ topological sort incomplete; missing SUB_IDs={missing}")
        order.extend(missing)
    return order


# Canonical ordering of return periods that ever appear in the schema.
# This list is the single source of truth for column ordering; CONFIG
# entries that are not in this list are ignored with a console warning.
_DESIGN_RETURN_PERIODS = (2, 5, 10, 20, 25, 50, 100, 200)


def parse_design_rainfall(cfg_dict):
    """
    Read CONFIG['RAINFALL_DEPTHS_MM'] and return the ordered list of
    (return_period, P_mm) pairs for which the user supplied a positive
    numeric value. Anything missing/None/non-positive is silently skipped.
    Unknown return periods produce a console warning.
    """
    if not isinstance(cfg_dict, dict) or not cfg_dict:
        return []
    # v22: ANY positive integer return period is accepted (the classic
    # 2..200 tuple is only the dialog's convenience default). Keys are
    # coerced to int; depths must be positive numbers.
    out = []
    for rp_key, v in cfg_dict.items():
        try:
            rp = int(round(float(rp_key)))
        except (TypeError, ValueError):
            print(f"⚠ RAINFALL_DEPTHS_MM: ignoring non-numeric "
                  f"return-period key {rp_key!r}.")
            continue
        if rp <= 0:
            print(f"⚠ RAINFALL_DEPTHS_MM: ignoring non-positive "
                  f"return period {rp_key!r}.")
            continue
        if v is None:
            continue
        try:
            v = float(v)
        except (TypeError, ValueError):
            continue
        if v > 0:
            out.append((rp, v))
    out.sort(key=lambda x: x[0])
    return out


def design_rainfall_field_blocks(active_periods, with_cn, with_qpeak=False):
    """
    Build the QgsField definitions for the design-rainfall block, given
    the user-supplied (rp, P_mm) pairs and whether CN is available.

    Returns: list of QgsField, in the order:
        [P_<n>_yrs,
         R_<n>_yrs (if CN),
         V_<n>yrs_1000m3 (if CN),
         QP_<n>_yrs (if CN and with_qpeak)]
        per period.

    QP_<n>_yrs is added only when both CN is available AND
    with_qpeak=True (i.e. v14's HYDROGRAPHS_ENABLED). Sub-basin layers
    pass with_qpeak=False so they only carry the P/R/V columns.
    """
    fields = []
    for rp, _ in active_periods:
        fields.append(
            QgsField(f'P_{rp}_yrs', _FT_DOUBLE, 'double', 10, 2)
        )
        if with_cn:
            fields.append(
                QgsField(f'R_{rp}_yrs', _FT_DOUBLE, 'double', 10, 2)
            )
            fields.append(
                QgsField(f'V_{rp}yrs_1000m3',
                         _FT_DOUBLE, 'double', 14, 3)
            )
            if with_qpeak:
                fields.append(
                    QgsField(f'QP_{rp}_yrs',
                             _FT_DOUBLE, 'double', 12, 3)
                )
    return fields


def design_rainfall_attr_block(active_periods, cn_avg, area_km2, with_cn,
                               qpeaks=None, scs_lambda=0.2):
    """
    Build the matching attribute values for the design-rainfall block.
    Uses scs_cn_runoff for R and Q*A for V; volume is in 10^3 m³.

    Parameters
    ----------
    active_periods : list of (rp, P_mm)
    cn_avg : float or None
    area_km2 : float
    with_cn : bool
    qpeaks : dict {rp: float} or None
        If provided AND with_cn AND any rp present here, the matching
        QP_<rp>_yrs column is appended after V (peak Q in m³/s).
        Missing periods get None (NULL in the table).

    Returns a list of values in the same order as
    design_rainfall_field_blocks(... , with_qpeak=qpeaks is not None).
    """
    out = []
    with_qpeak = qpeaks is not None
    for rp, p_mm in active_periods:
        out.append(float(p_mm))
        if with_cn:
            r_mm = scs_cn_runoff(p_mm, cn_avg, lam=scs_lambda)
            v_1000m3 = r_mm * float(area_km2)   # mm·km² == 10^3 m³
            out.append(float(r_mm))
            out.append(float(v_1000m3))
            if with_qpeak:
                qp = qpeaks.get(rp) if isinstance(qpeaks, dict) else None
                out.append(float(qp) if qp is not None else None)
    return out


def compute_lfp_characteristics(line_geom, dem_layer):
    """
    Computes longest-flow-path metrics, three sub-segment slopes, and the
    Kirpich time of concentration for each sub-segment plus totals.

    Sub-segments along the LFP from upstream (US) to outlet (DS):
        Segment 1: US (0% from start)  → P85 (15% from start)  length = 0.15·L
        Segment 2: P85 (15%)           → P10 (90%)             length = 0.75·L
        Segment 3: P10 (90%)           → DS  (100%)            length = 0.10·L
    """
    if line_geom is None or line_geom.isNull():
        return None
    if line_geom.isMultipart():
        parts = line_geom.asMultiPolyline()
        if not parts:
            return None
        line_geom = QgsGeometry.fromPolylineXY(parts[0])
    L = line_geom.length()
    if L == 0:
        return None

    verts = line_geom.asPolyline()
    pt_us = QgsPointXY(verts[0])
    pt_ds = QgsPointXY(verts[-1])

    d10_outlet = L * 0.10
    d85_outlet = L * 0.85
    d10_start  = L - d10_outlet     # 90 % from start (P10)
    d85_start  = L - d85_outlet     # 15 % from start (P85)
    g10 = line_geom.interpolate(d10_start)
    g85 = line_geom.interpolate(d85_start)
    pt_10 = QgsPointXY(g10.asPoint()) if g10 and not g10.isNull() else None
    pt_85 = QgsPointXY(g85.asPoint()) if g85 and not g85.isNull() else None

    ele_us = sample_dem(pt_us, dem_layer) or 0.0
    ele_ds = sample_dem(pt_ds, dem_layer) or 0.0
    ele_10 = (sample_dem(pt_10, dem_layer) or 0.0) if pt_10 else 0.0
    ele_85 = (sample_dem(pt_85, dem_layer) or 0.0) if pt_85 else 0.0

    # --- segment lengths along the LFP ----------------------------------
    L_us_85 = 0.15 * L
    L_85_10 = 0.75 * L
    L_10_ds = 0.10 * L

    # --- segment slopes (m/m).  Drop is upstream-elev minus downstream-elev,
    #     so a positive slope means the path actually descends in that segment.
    slp_us_85 = (ele_us - ele_85) / L_us_85 if L_us_85 > 0 else 0.0
    slp_85_10 = (ele_85 - ele_10) / L_85_10 if L_85_10 > 0 else 0.0
    slp_10_ds = (ele_10 - ele_ds) / L_10_ds if L_10_ds > 0 else 0.0

    # --- Kirpich Tc per segment (minutes), plus totals ------------------
    tc_us_85 = kirpich_tc(L_us_85, slp_us_85)
    tc_85_10 = kirpich_tc(L_85_10, slp_85_10)
    tc_10_ds = kirpich_tc(L_10_ds, slp_10_ds)
    tc_total = tc_us_85 + tc_85_10 + tc_10_ds
    t_lag    = 0.6 * tc_total

    # --- Holistic Kirpich Tc (v14) — single segment using full L and the
    #     end-to-end average slope (ELE_US - ELE_DS) / L ------------------
    tc_kirph   = kirpich_tc_holistic(L, ele_us, ele_ds)
    tlag_kirph = 0.6 * tc_kirph

    # legacy fields kept for backwards compatibility (== SLP_85_10)
    slope    = slp_85_10
    slope_p  = slope * 100

    return dict(
        # geometry & elevations
        L=L, ELE_US=ele_us, ELE_DS=ele_ds, ELE_10=ele_10, ELE_85=ele_85,
        PT_US=pt_us, PT_DS=pt_ds, PT_10=pt_10, PT_85=pt_85,
        D10=d10_start, D85=d85_start,
        # legacy 10-85 slope
        SLOPE=slope, SLOPE_PCT=slope_p,
        # NEW: segment slopes (m/m)
        SLP_US_85=slp_us_85, SLP_85_10=slp_85_10, SLP_10_DS=slp_10_ds,
        # NEW: Kirpich Tc (minutes)
        TC_US_85=tc_us_85, TC_85_10=tc_85_10, TC_10_DS=tc_10_ds,
        TC_TOTAL=tc_total, TLAG=t_lag,
        # v14: holistic Kirpich Tc / Tlag (minutes)
        TC_KIRPH=tc_kirph, TLAG_KIRPH=tlag_kirph,
    )


# =============================================================================
# 8. SUB-BASIN SUBDIVISION
# =============================================================================
# All functions below operate on the parent catchment AFTER it has been
# delineated. They never need to re-run pit-fill, flow direction, or flow
# accumulation — those are reused as-is.
#
# Vocabulary used throughout this section:
#   junction         a (row, col) cell that is treated as an internal pour
#                    point. The set of junctions always includes the basin
#                    outlet, plus every confluence inside the basin streams,
#                    plus any cells inserted to keep sub-basins close to
#                    the user-specified target area.
#   sub-basin        the LOCAL contributing area of one junction = its full
#                    upstream catchment minus the catchments of any
#                    junctions further upstream. Every cell in the parent
#                    catchment belongs to exactly one sub-basin.
#   reach            the channel polyline from one junction down to its
#                    immediately-downstream junction (one per non-outlet
#                    sub-basin).
# =============================================================================

def find_confluences_in_basin(in_basin_streams, fdir):
    """Return set of (r,c) cells that are stream cells inside the basin
    AND have ≥2 upstream neighbours that are also in-basin stream cells.
    Those are the natural sub-basin boundaries."""
    rows, cols = fdir.shape
    confluences = set()
    rs, cs = np.where(in_basin_streams)
    for r, c in zip(rs.tolist(), cs.tolist()):
        n_up = 0
        for (dr, dc), need_code in INV_REQUIRED_CODE.items():
            nr, nc = r + dr, c + dc
            if (0 <= nr < rows and 0 <= nc < cols
                    and in_basin_streams[nr, nc]
                    and fdir[nr, nc] == need_code):
                n_up += 1
                if n_up >= 2:
                    confluences.add((r, c))
                    break
    return confluences


def place_junctions(catchment_mask, streams, fdir, facc,
                    outlet_rc, target_cells):
    """
    Build the set of junctions inside a catchment.

    Algorithm:
      1. Start with {outlet} ∪ {every confluence inside the basin streams}.
      2. BFS from the outlet upward along the stream tree. While walking
         each stem from one junction j upstream, track the local area
         (= facc[j] - facc[current_cell]). If it crosses target_cells
         before another junction is reached, place a new junction at the
         current cell and continue from there.

    Returns: list of (row, col) tuples (order is irrelevant for downstream
             uses; topology is figured out separately).
    """
    in_basin_streams = (streams.astype(bool) & catchment_mask)
    rows, cols = fdir.shape

    junctions = set(find_confluences_in_basin(in_basin_streams, fdir))
    junctions.add(tuple(outlet_rc))

    # If subdivision is disabled (target_cells <= 0), return only the outlet.
    if target_cells <= 0:
        return [tuple(outlet_rc)]

    queue = deque([tuple(outlet_rc)])
    processed = {tuple(outlet_rc)}

    while queue:
        j = queue.popleft()
        # Find each immediately-upstream stream neighbour of j → start of a stem
        for (dr, dc), need_code in INV_REQUIRED_CODE.items():
            nr, nc = j[0] + dr, j[1] + dc
            if not (0 <= nr < rows and 0 <= nc < cols):
                continue
            if not in_basin_streams[nr, nc] or fdir[nr, nc] != need_code:
                continue

            # ----------------------------------------------------------------
            # Walk this single-stem branch upstream from (nr, nc) toward the
            # divide. The reference for "local area accumulated since the
            # last junction" is the facc value of the first cell ABOVE the
            # current downstream junction on THIS branch — NOT facc[j],
            # because at a confluence j receives flow from multiple branches
            # and facc[j] mixes them all.
            # local_area = facc[branch_start_cell] - facc[current_cell]
            # is then valid as long as we stay on a single stem (which we do
            # because we stop at every confluence).
            # ----------------------------------------------------------------
            current = (nr, nc)
            ref_facc = int(facc[current[0], current[1]])  # branch base

            while True:
                if current in junctions:
                    if current not in processed:
                        processed.add(current)
                        queue.append(current)
                    break

                local_area = ref_facc - int(facc[current[0], current[1]])
                if local_area >= target_cells:
                    junctions.add(current)
                    processed.add(current)
                    queue.append(current)
                    ref_facc = int(facc[current[0], current[1]])

                # find unique upstream stream cell
                ups = []
                for (dr2, dc2), code2 in INV_REQUIRED_CODE.items():
                    nr2, nc2 = current[0] + dr2, current[1] + dc2
                    if (0 <= nr2 < rows and 0 <= nc2 < cols
                            and in_basin_streams[nr2, nc2]
                            and fdir[nr2, nc2] == code2):
                        ups.append((nr2, nc2))
                if not ups:                       # headwater
                    break
                if len(ups) >= 2:                 # confluence (already in set)
                    break
                current = ups[0]

    return list(junctions)


def trace_to_next_junction(start_rc, junctions_set, fdir, gt,
                           max_steps=1000000):
    """
    Trace flow direction downstream from start_rc until we hit a cell that
    is in junctions_set (excluding start_rc itself). Returns:
        (downstream_junction_rc or None,
         list of (row,col) cells on the path including both endpoints,
         list of (x,y) coords for the same cells)
    Returns (None, [start_rc], [start_xy]) if no junction is found
    downstream (i.e. start_rc IS the outlet of the parent catchment).
    """
    rows, cols = fdir.shape
    r, c = start_rc
    path_rc = [start_rc]
    seen = {start_rc}
    for _ in range(max_steps):
        d = int(fdir[r, c])
        if d not in CODE_TO_OFFSET:
            break
        dr, dc = CODE_TO_OFFSET[d]
        nr, nc = r + dr, c + dc
        if not (0 <= nr < rows and 0 <= nc < cols):
            break
        if (nr, nc) in seen:
            break
        r, c = nr, nc
        path_rc.append((r, c))
        seen.add((r, c))
        if (r, c) in junctions_set and (r, c) != start_rc:
            ds = (r, c)
            coords = [cell_to_xy(rr, cc, gt) for rr, cc in path_rc]
            return ds, path_rc, coords
    coords = [cell_to_xy(rr, cc, gt) for rr, cc in path_rc]
    return None, path_rc, coords


def compute_subbasin_id_grid(junctions_list, junctions_set,
                             fdir, parent_catch):
    """
    Assign each cell in parent_catch the SUB_ID of the FIRST junction
    encountered when walking downstream along the flow direction grid.

    Returns: int32 array, same shape as fdir, with values:
        0           cells not in parent_catch
        1, 2, 3, …  the SUB_ID of the owner junction (= index into
                    junctions_list + 1).

    SUB_ID assignment is intentionally arbitrary at this stage. Callers
    re-number afterwards so SUB_ID 1 is always the outlet sub-basin.
    """
    rows, cols = fdir.shape
    j_to_sid = {j: i + 1 for i, j in enumerate(junctions_list)}

    sid_grid = np.zeros((rows, cols), dtype=np.int32)
    # seed: each junction owns its own cell
    for j, sid in j_to_sid.items():
        sid_grid[j[0], j[1]] = sid

    rs, cs = np.where(parent_catch)
    for r0, c0 in zip(rs.tolist(), cs.tolist()):
        if sid_grid[r0, c0] != 0:
            continue
        path = [(r0, c0)]
        r, c = r0, c0
        sid = 0
        while True:
            d = int(fdir[r, c])
            if d not in CODE_TO_OFFSET:
                break
            dr, dc = CODE_TO_OFFSET[d]
            nr, nc = r + dr, c + dc
            if not (0 <= nr < rows and 0 <= nc < cols):
                break
            if sid_grid[nr, nc] != 0:
                sid = int(sid_grid[nr, nc])
                break
            path.append((nr, nc))
            r, c = nr, nc
        for pr, pc in path:
            sid_grid[pr, pc] = sid
    return sid_grid


def build_subbasin_topology(junctions_list, junctions_set, fdir,
                            outlet_rc, gt):
    """
    For every junction, find its immediately-downstream junction and the
    geometry of the channel reach connecting them.

    Returns three dicts keyed by junction (row, col):
        ds_junction[j]   — (row, col) of immediately-downstream junction,
                            or None if j IS the outlet.
        reach_coords[j]  — list of (x,y) along the flow path from j down
                            to ds_junction[j] (inclusive on both ends).
                            Empty list for the outlet.
        reach_cells[j]   — same as above but as (row, col) tuples.
    """
    ds_junction  = {}
    reach_coords = {}
    reach_cells  = {}
    for j in junctions_list:
        if j == tuple(outlet_rc):
            ds_junction[j]  = None
            reach_coords[j] = []
            reach_cells[j]  = []
            continue
        ds, cells_rc, coords = trace_to_next_junction(
            j, junctions_set, fdir, gt)
        ds_junction[j]  = ds
        reach_cells[j]  = cells_rc
        reach_coords[j] = coords
    return ds_junction, reach_cells, reach_coords


def renumber_junctions_outlet_first(junctions_list, ds_junction,
                                    outlet_rc):
    """
    Build the SUB_ID mapping so that SUB_ID 1 is always the outlet
    sub-basin and successive IDs walk upstream in BFS order. Also
    invert ds_junction to give DS_SUB_ID for every junction.

    Returns:
        sid_of_junction  : dict (row,col) → SUB_ID (int >= 1)
        ds_sid_of_sid    : dict SUB_ID  → DS_SUB_ID (0 for outlet)
    """
    outlet = tuple(outlet_rc)
    # who flows directly into j ?
    upstream_of = {j: [] for j in junctions_list}
    for j, ds in ds_junction.items():
        if ds is not None:
            upstream_of[ds].append(j)

    sid_of_junction = {}
    order = []
    queue = deque([outlet])
    sid = 1
    while queue:
        j = queue.popleft()
        if j in sid_of_junction:
            continue
        sid_of_junction[j] = sid
        order.append(j)
        sid += 1
        for up in upstream_of.get(j, []):
            queue.append(up)
    # cover any junctions not connected (shouldn't happen, but safe)
    for j in junctions_list:
        if j not in sid_of_junction:
            sid_of_junction[j] = sid
            order.append(j)
            sid += 1

    ds_sid_of_sid = {}
    for j, this_sid in sid_of_junction.items():
        ds = ds_junction.get(j)
        ds_sid_of_sid[this_sid] = sid_of_junction[ds] if ds is not None else 0
    return sid_of_junction, ds_sid_of_sid


def reach_attributes(reach_cells, dem_layer, gt):
    """
    Length (m), upstream/downstream elevation (m), slope (m/m and %).
    `reach_cells` is a list of (row, col) from upstream junction down to
    downstream junction. Length is computed along the flow path with
    proper diagonal segment weighting.
    """
    if not reach_cells or len(reach_cells) < 2:
        return None
    # length along the polyline (uses cell-centre coordinates)
    pts = [cell_to_xy(r, c, gt) for r, c in reach_cells]
    L = 0.0
    for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:]):
        L += math.hypot(x2 - x1, y2 - y1)

    pt_us = QgsPointXY(*pts[0])
    pt_ds = QgsPointXY(*pts[-1])
    ele_us = sample_dem(pt_us, dem_layer) or 0.0
    ele_ds = sample_dem(pt_ds, dem_layer) or 0.0
    slope = ((ele_us - ele_ds) / L) if L > 0 else 0.0
    return dict(L=L, ELE_US=ele_us, ELE_DS=ele_ds,
                SLOPE=slope, SLOPE_PCT=slope * 100.0)


def lfp_within_mask(sub_mask, outlet_rc, fdir, dx, dy, gt, dem_layer):
    """
    Longest flow path WITHIN sub_mask, ending at outlet_rc. Returns the
    QgsGeometry line plus the dict from compute_lfp_characteristics, or
    (None, None) if no valid path is found.
    """
    fdist = flow_distance(fdir, outlet_rc[0], outlet_rc[1], dx, dy)
    fdist_in = np.where(sub_mask, fdist, np.nan)
    if np.all(np.isnan(fdist_in)):
        return None, None
    max_idx = int(np.nanargmax(fdist_in))
    src_r, src_c = np.unravel_index(max_idx, fdist_in.shape)
    path_cells = trace_flow_path(fdir, src_r, src_c,
                                 stop_row=outlet_rc[0], stop_col=outlet_rc[1])
    if len(path_cells) < 2:
        return None, None
    coords = [cell_to_xy(r, c, gt) for r, c in path_cells]
    line_geom = QgsGeometry.fromPolylineXY(
        [QgsPointXY(x, y) for x, y in coords])
    ch = compute_lfp_characteristics(line_geom, dem_layer)
    return line_geom, ch


# =============================================================================
# 9. MAIN WORKFLOW
# =============================================================================
def thiessen_index_grid(stations, gt, shape, max_cells=4_000_000):
    """v1.0 — nearest-station index grid (1-based), the raster form of
    the Voronoi diagram. Pure NumPy, so it carries no GDAL/QGIS
    dependency and can be verified directly against the weights that
    thiessen_basin_depths derives from the same rule.

    Returns (idx, geotransform, step); the grid is coarsened by `step`
    when it would otherwise exceed max_cells."""
    if not stations or gt is None:
        return None, None, 1
    nrows, ncols = int(shape[0]), int(shape[1])
    if nrows < 2 or ncols < 2:
        return None, None, 1
    step = 1
    while (nrows // step) * (ncols // step) > max_cells:
        step += 1
    r_idx = np.arange(0, nrows, step, dtype=np.float64)
    c_idx = np.arange(0, ncols, step, dtype=np.float64)
    gt_o = (gt[0], gt[1] * step, gt[2] * step,
            gt[3], gt[4], gt[5] * step)
    sx = np.asarray([s['x'] for s in stations], dtype=np.float64)
    sy = np.asarray([s['y'] for s in stations], dtype=np.float64)
    idx = np.zeros((r_idx.size, c_idx.size), dtype=np.int32)
    block = max(1, int(2_000_000 / max(1, c_idx.size * len(stations))))
    for i0 in range(0, r_idx.size, block):
        rr = r_idx[i0:i0 + block]
        x = gt[0] + (c_idx[None, :] + 0.5) * gt[1] \
            + (rr[:, None] + 0.5) * gt[2]
        y = gt[3] + (c_idx[None, :] + 0.5) * gt[4] \
            + (rr[:, None] + 0.5) * gt[5]
        d2 = ((x[None, ...] - sx[:, None, None]) ** 2
              + (y[None, ...] - sy[:, None, None]) ** 2)
        idx[i0:i0 + block, :] = np.argmin(d2, axis=0).astype(np.int32) + 1
    return idx, gt_o, step


def build_thiessen_layer(stations, gt, shape, crs_authid,
                        layer_name='Thiessen', max_cells=4_000_000):
    """v1.0 — Voronoi (Thiessen) polygons of the rainfall stations.

    Every cell of the DEM grid is assigned to its nearest station and
    the resulting index grid is polygonised, so the polygons are
    exactly the areas that weight the station depths (the same
    nearest-station rule thiessen_basin_depths applies). On very large
    grids the diagram is built on a coarsened grid — which changes only
    the smoothness of the drawn boundary, never the weights, since
    those are always computed at full resolution.

    Returns a memory QgsVectorLayer, or None when it cannot be built.
    """
    if not stations or gt is None:
        return None
    try:
        idx, gt_o, step = thiessen_index_grid(stations, gt, shape,
                                              max_cells=max_cells)
        if idx is None:
            return None
        # ---- polygonise the index grid ------------------------------
        mem = gdal.GetDriverByName('MEM').Create(
            '', int(idx.shape[1]), int(idx.shape[0]), 1, gdal.GDT_Int32)
        mem.SetGeoTransform(gt_o)
        srs = osr.SpatialReference()
        try:
            srs.SetFromUserInput(str(crs_authid))
            mem.SetProjection(srs.ExportToWkt())
        except Exception:
            srs = None
        mem.GetRasterBand(1).WriteArray(idx)
        ogr_ds = ogr.GetDriverByName('Memory').CreateDataSource('thiessen')
        ogr_lyr = ogr_ds.CreateLayer('t', srs=srs, geom_type=ogr.wkbPolygon)
        ogr_lyr.CreateField(ogr.FieldDefn('STA_IDX', ogr.OFTInteger))
        gdal.Polygonize(mem.GetRasterBand(1), None, ogr_lyr, 0)

        # ---- assemble the QGIS layer --------------------------------
        rps = sorted(int(r) for r in (stations[0].get('depths') or {}).keys())
        vl = QgsVectorLayer(f"Polygon?crs={crs_authid}", layer_name, "memory")
        flds = [QgsField('STA_IDX', _FT_INT),
                QgsField('STATION', _FT_STRING, 'string', 64),
                QgsField('AREA_KM2', _FT_DOUBLE, 'double', 15, 4)]
        for rp in rps:
            flds.append(QgsField(f'P_{rp}_yrs',
                                 _FT_DOUBLE, 'double', 10, 2))
        vl.dataProvider().addAttributes(flds)
        vl.updateFields()

        feats = []
        for f in ogr_lyr:
            g = f.GetGeometryRef()
            if g is None:
                continue
            i = int(f.GetField('STA_IDX')) - 1
            if i < 0 or i >= len(stations):
                continue
            qg = QgsGeometry.fromWkt(g.ExportToWkt())
            if qg is None or qg.isEmpty():
                continue
            st = stations[i]
            attrs = [i + 1, str(st.get('name') or f'Station {i + 1}'),
                     float(qg.area()) / 1.0e6]
            for rp in rps:
                v = (st.get('depths') or {}).get(rp)
                attrs.append(float(v) if v is not None else None)
            qf = QgsFeature()
            qf.setGeometry(qg)
            qf.setAttributes(attrs)
            feats.append(qf)
        if not feats:
            return None
        vl.dataProvider().addFeatures(feats)
        vl.updateExtents()
        if step > 1:
            print(f"     (diagram drawn on every {step}th cell for size; "
                  f"the weights themselves use the full grid)")
        return vl
    except Exception as e:
        print(f"     \u26a0 Thiessen polygon layer not created: {e}")
        return None


def thiessen_basin_depths(mask, gt, stations, chunk=2_000_000):
    """Area-weighted rainfall depths for one basin via Thiessen polygons.

    v22 (spatially-variable rainfall). Every cell of the basin ``mask``
    is assigned to its nearest rainfall station (squared Euclidean
    distance in map units) — an exact raster-space evaluation of the
    Thiessen / Voronoi diagram, so weights are the fraction of basin
    cells falling in each station's polygon. Depths are then the
    weight-averaged station values per return period.

    Parameters
    ----------
    mask : 2-D bool array          basin cells on the DEM grid
    gt   : GDAL geotransform       for cell-centre coordinates
    stations : list of dict(x=..., y=..., name=..., depths={rp: mm})
    chunk : int                    cells per distance batch (memory cap)

    Returns
    -------
    (depths, weights) : ({rp: mm}, {station_name: fraction}) or (None, None)
    """
    if not stations:
        return None, None
    rows, cols = np.where(mask)
    n = rows.size
    if n == 0:
        return None, None
    sx = np.asarray([s['x'] for s in stations], dtype=np.float64)
    sy = np.asarray([s['y'] for s in stations], dtype=np.float64)
    counts = np.zeros(len(stations), dtype=np.int64)
    for i0 in range(0, n, chunk):
        r = rows[i0:i0 + chunk].astype(np.float64)
        c = cols[i0:i0 + chunk].astype(np.float64)
        x = gt[0] + (c + 0.5) * gt[1] + (r + 0.5) * gt[2]
        y = gt[3] + (c + 0.5) * gt[4] + (r + 0.5) * gt[5]
        d2 = ((x[None, :] - sx[:, None]) ** 2
              + (y[None, :] - sy[:, None]) ** 2)
        nearest = np.argmin(d2, axis=0)
        counts += np.bincount(nearest, minlength=len(stations))
    w = counts.astype(np.float64) / float(n)
    weights = {stations[i]['name']: float(w[i])
               for i in range(len(stations)) if w[i] > 0}
    # all stations share the same return-period keys (validated on load)
    rps = list(stations[0]['depths'].keys())
    depths = {}
    for rp in rps:
        vals = np.asarray([float(s['depths'][rp]) for s in stations])
        depths[rp] = float(np.sum(w * vals))
    return depths, weights


def mask_centroid_xy(mask, gt):
    """Centroid of a raster mask in map coordinates (v22.1).

    Mean of the cell-centre coordinates of all True cells — the exact
    discrete centroid of the basin on the DEM grid. Used to place
    HEC-HMS Subbasin icons at each sub-basin's centroid rather than on
    its outlet junction.
    """
    rows, cols = np.where(mask)
    if rows.size == 0:
        return None
    r = rows.mean()
    c = cols.mean()
    x = gt[0] + (c + 0.5) * gt[1] + (r + 0.5) * gt[2]
    y = gt[3] + (c + 0.5) * gt[4] + (r + 0.5) * gt[5]
    return (float(x), float(y))


def load_rain_stations(layer_name, fields_map, dem_crs, project):
    """Load the rainfall-stations point layer for Thiessen mode (v22).

    Returns (stations, rps) where stations is a list of
    dict(x, y, name, depths={rp: mm}) in the DEM CRS, or (None, None)
    when the mode is inactive or the layer/fields are unusable. Any
    station missing a positive numeric value in a mapped field is
    skipped with a warning.
    """
    if not layer_name or not fields_map:
        return None, None
    try:
        lyr = get_layer(layer_name)
    except Exception as e:
        print(f"⚠ Rain-stations layer '{layer_name}' not found ({e}) "
              f"— Thiessen mode disabled.")
        return None, None
    # normalise mapping: rp(int) -> field(str)
    fmap = {}
    for k, fname in fields_map.items():
        try:
            rp = int(round(float(k)))
        except (TypeError, ValueError):
            continue
        if rp > 0 and fname:
            fmap[rp] = str(fname)
    if not fmap:
        print("⚠ RAIN_STATIONS_FIELDS empty/invalid — Thiessen "
              "mode disabled.")
        return None, None
    missing = [f for f in fmap.values()
               if lyr.fields().indexOf(f) < 0]
    if missing:
        print(f"⚠ Rain-stations layer lacks field(s) {missing} — "
              f"Thiessen mode disabled.")
        return None, None
    xform = None
    if lyr.crs() != dem_crs:
        xform = QgsCoordinateTransform(lyr.crs(), dem_crs, project)
    stations = []
    for i, f in enumerate(lyr.getFeatures()):
        g = f.geometry()
        if g is None or g.isEmpty():
            continue
        pt = g.asPoint() if not g.isMultipart() else g.asMultiPoint()[0]
        if xform is not None:
            pt = xform.transform(pt)
        depths = {}
        ok = True
        for rp, fname in fmap.items():
            v = f[fname]
            try:
                v = float(v)
            except (TypeError, ValueError):
                ok = False
                break
            if not (v > 0) or not np.isfinite(v):
                ok = False
                break
            depths[rp] = v
        name = None
        for cand in ('name', 'Name', 'NAME', 'station', 'Station', 'ID',
                     'id', 'Id'):
            if lyr.fields().indexOf(cand) >= 0 and f[cand] is not None:
                name = str(f[cand])
                break
        if name is None:
            name = f"Station_{i + 1}"
        if ok:
            stations.append(dict(x=float(pt.x()), y=float(pt.y()),
                                 name=name, depths=depths))
        else:
            print(f"   ⚠ station '{name}' skipped — missing/"
                  f"invalid depth in one of {list(fmap.values())}.")
    if not stations:
        print("⚠ No usable rainfall stations — Thiessen mode "
              "disabled.")
        return None, None
    return stations, sorted(fmap.keys())


# =============================================================================
# v24 — Reservoir storage routing (modified Puls / storage-indication)
# -----------------------------------------------------------------------------
# Faithful port of the AR-Awdiya reservoir engine, so a dam defined in
# the "Special Delineation Points" tab routes IDENTICALLY here and in
# AR-Awdiya after .basin export. Two definition modes:
#   'structures' — elevation–volume curve + outflow structures
#                  (broad-crested spillway weir and/or pipe orifice);
#                  stage–discharge is summed across structures and
#                  converted to a storage–discharge curve via the EV
#                  curve (80 stage samples + every breakpoint).
#   'table'      — direct storage–discharge table.
# Routing solves  (2S₂/Δt + O₂) = (I₁ + I₂) + (2S₁/Δt − O₁).
# =============================================================================
G_ACCEL = 9.81


def _lerp_curve(x, xs, ys):
    """Clamped piecewise-linear interpolation (xs ascending)."""
    if not xs:
        return 0.0
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(1, len(xs)):
        if xs[i] >= x:
            span = xs[i] - xs[i - 1]
            t = (x - xs[i - 1]) / span if span else 0.0
            return ys[i - 1] + t * (ys[i] - ys[i - 1])
    return ys[-1]


def parse_storage_table(text):
    """Rows 'storage(10³ m³), outflow(m³/s)' → (S[m³], O[m³/s]) ascending,
    deduplicated, anchored at (0, 0)."""
    rows = []
    import re as _re
    for line in _re.split(r'[\n;]+', text or ''):
        p = [s for s in _re.split(r'[\s,\t]+', line.strip()) if s]
        if len(p) >= 2:
            try:
                s_, o_ = float(p[0]), float(p[1])
            except ValueError:
                continue
            if np.isfinite(s_) and np.isfinite(o_) and s_ >= 0 and o_ >= 0:
                rows.append((s_ * 1000.0, o_))         # 10³ m³ → m³
    rows.sort(key=lambda r: r[0])
    S, O = [], []  # noqa: E741
    last = -1.0
    for s_, o_ in rows:
        if s_ == last:
            continue
        S.append(s_); O.append(o_); last = s_
    if S and S[0] > 0:
        S.insert(0, 0.0); O.insert(0, 0.0)
    return S, O


def parse_ev(text):
    """Rows 'elevation(m), volume(10³ m³)' → (E[m] ascending, V[m³])."""
    rows = []
    import re as _re
    for line in _re.split(r'[\n;]+', text or ''):
        p = [s for s in _re.split(r'[\s,\t]+', line.strip()) if s]
        if len(p) >= 2:
            try:
                e_, v_ = float(p[0]), float(p[1])
            except ValueError:
                continue
            if np.isfinite(e_) and np.isfinite(v_) and v_ >= 0:
                rows.append((e_, v_ * 1000.0))
    rows.sort(key=lambda r: r[0])
    E, V = [], []
    last = None
    for e_, v_ in rows:
        if e_ == last:
            continue
        E.append(e_); V.append(v_); last = e_
    return E, V


def orifice_geometry(st):
    """v1.0 — (area m2, centre elevation m) of an orifice outlet.

    Circular  : A = pi*r^2,  centre = invert + r
    Rectangular: A = S*H,    centre = invert + H/2
    The centre elevation is the head datum for the orifice equation and
    is exactly what HEC-HMS asks for as 'Centerline Elevation'."""
    inv = float(st.get('invert', 0.0))
    if str(st.get('shape', 'circle')).lower().startswith('rect'):
        s = max(1e-6, float(st.get('span', 0.0)))
        h = max(1e-6, float(st.get('height', 0.0)))
        return s * h, inv + 0.5 * h
    r = max(1e-4, float(st.get('radius', 0.0)))
    return math.pi * r * r, inv + r


def structure_q(st, stage):
    """Discharge (m³/s) from one outflow structure at water-surface
    elevation `stage` (m). Spillway: Q = C·L·H^1.5 (broad-crested
    weir). Pipe: HEC-HMS orifice Q = n·Cd·A·√(2g·Hc), Hc measured
    above the orifice CENTRE (invert + radius)."""
    if not st:
        return 0.0
    typ = st.get('type')
    if typ == 'spillway':
        H = stage - float(st.get('crest', 0.0))
        if H <= 0:
            return 0.0
        return (max(0.0, float(st.get('C', 1.7)))
                * max(0.0, float(st.get('width', 0.0)))
                * H ** 1.5)
    if typ == 'pipe':
        n = max(1, int(round(float(st.get('n', 1) or 1))))
        cd = max(0.05, float(st.get('Cd', 0.6) or 0.6))
        a_pipe, centre = orifice_geometry(st)
        hc = stage - centre
        if hc <= 0:
            return 0.0
        return n * cd * a_pipe * math.sqrt(2.0 * G_ACCEL * hc)
    return 0.0


def reservoir_sd(res):
    """Build the storage–discharge routing curve for one reservoir
    definition. Returns dict(S, O, warn, info, Saxis, Eaxis,
    has_stage, ev_rows) or None when `res` is falsy. S/O with fewer
    than 2 points means the definition is unusable (see warn)."""
    if not res:
        return None
    if res.get('mode') == 'structures':
        E, V = parse_ev(res.get('ev_text'))
        if len(E) < 2:
            return dict(S=[], O=[], warn='elevation–volume curve needs '
                        '≥2 rows', info='reservoir (no EV curve)',
                        Saxis=None, Eaxis=None, has_stage=False,
                        ev_rows=None)
        sts = [s for s in (res.get('structures') or [])
               if s and s.get('type')]
        if not sts:
            return dict(S=[], O=[], warn='no outflow structures defined',
                        info='reservoir (no structures)',
                        Saxis=None, Eaxis=None, has_stage=False,
                        ev_rows=list(zip(E, V)))
        emin, emax = E[0], E[-1]
        stage_set = {emin + (emax - emin) * i / 80.0 for i in range(81)}
        stage_set.update(E)
        for s in sts:
            if s['type'] == 'spillway':
                stage_set.add(float(s.get('crest', 0.0)))
            elif s['type'] == 'pipe':
                stage_set.add(orifice_geometry(s)[1])
        stages = sorted(e for e in stage_set
                        if emin - 1e-9 <= e <= emax + 1e-9)
        v0 = _lerp_curve(emin, E, V)
        S, O = [], []  # noqa: E741
        last_s = -1.0
        for e in stages:
            s_ = max(0.0, _lerp_curve(e, E, V) - v0)
            o_ = sum(structure_q(st, e) for st in sts)
            if O and abs(s_ - last_s) < 1e-6:
                O[-1] = max(O[-1], o_)
                continue
            S.append(s_); O.append(o_); last_s = s_
        Saxis, Eaxis = [], []
        last_se = -1.0
        for e in stages:
            s_ = max(0.0, _lerp_curve(e, E, V) - v0)
            if Saxis and abs(s_ - last_se) < 1e-6:
                continue
            Saxis.append(s_); Eaxis.append(e); last_se = s_
        if S and S[0] > 0:
            S.insert(0, 0.0); O.insert(0, 0.0)
        cap = ' + '.join(
            (f"weir C={s.get('C')}, L={s.get('width')} m, "
             f"crest {s.get('crest')} m") if s['type'] == 'spillway'
            else (f"{s.get('n', 1)}×orifice "
                  + (f"{s.get('span')}×{s.get('height')} m rect"
                     if str(s.get('shape', 'circle')).lower()
                     .startswith('rect')
                     else f"r={s.get('radius')} m")
                  + f", Cd={s.get('Cd', 0.6)}, invert {s.get('invert')} m")
            for s in sts)
        return dict(S=S, O=O, warn=None,
                    info=f'reservoir (outflow structures: {cap})',
                    Saxis=Saxis, Eaxis=Eaxis, has_stage=True,
                    ev_rows=list(zip(E, V)))
    S, O = parse_storage_table(res.get('table_text'))  # noqa: E741
    return dict(S=S, O=O,
                warn=('storage–discharge table needs ≥2 rows'
                      if len(S) < 2 else None),
                info=f'reservoir (storage–discharge table, '
                     f'{len(S)}-pt)',
                Saxis=None, Eaxis=None, has_stage=False, ev_rows=None)


def storage_at_elevation(sd, elev):
    """v1.0 — pool storage (m3, above the curve datum) at water-surface
    elevation `elev`, or None when the definition carries no elevation
    information (storage-discharge table mode). Values outside the
    curve are clamped to its ends."""
    if not sd or not sd.get('has_stage'):
        return None
    E, S = sd.get('Eaxis'), sd.get('Saxis')
    if not E or not S or len(E) < 2:
        return None
    return max(0.0, _lerp_curve(float(elev), E, S))


# NOTE on the `# noqa: E741` markers in the three functions below:
# `I` and `O` are the standard symbols of the modified-Puls
# (storage-indication) formulation, and `O` is also the published key of
# the storage-discharge dict that hec_hms_export, awdiya_export and the
# dialog all read. Renaming them would obscure the hydrology and break
# that contract, so the style rule is suppressed rather than obeyed.
def route_reservoir(inflow, sd, dt_sec, init_storage=None):
    """Modified-Puls (storage-indication) routing of `inflow` through
    the storage–discharge curve `sd` (from reservoir_sd). Returns
    dict(O, storage, stage, max_storage, max_stage, inflow_peak,
    warn, changed). When the curve is unusable the inflow is passed
    through unrouted (changed=False)."""
    I = np.asarray(inflow, dtype=float)  # noqa: E741
    if not sd or not sd.get('S') or len(sd['S']) < 2:
        return dict(O=I, storage=None, stage=None, max_storage=None,
                    max_stage=None, inflow_peak=None,
                    warn=(sd or {}).get('warn')
                    or 'reservoir not defined — passed through unrouted',
                    changed=False,
                    info=(sd or {}).get('info') or 'reservoir (undefined)')
    S = list(sd['S']); O = list(sd['O']); n = len(S)  # noqa: E741
    A = [2.0 * S[i] / dt_sec + O[i] for i in range(n)]   # 2S/Δt + O ↑

    def o_from_a(a):
        if a <= A[0]:
            return max(0.0, O[0] * a / A[0]) if A[0] > 0 else 0.0
        if a >= A[-1]:
            t = (a - A[-2]) / ((A[-1] - A[-2]) or 1.0)
            return O[-2] + t * (O[-1] - O[-2])
        for i in range(1, n):
            if A[i] >= a:
                t = (a - A[i - 1]) / ((A[i] - A[i - 1]) or 1.0)
                return O[i - 1] + t * (O[i] - O[i - 1])
        return O[-1]

    m = len(I)
    out = np.zeros(m); stor = np.zeros(m)
    # v1.0 — initial condition. Default: the pool sits at the curve
    # datum (empty flood storage above it). With an initial water
    # level the run starts from the storage that level implies, and
    # the matching outflow is read off the storage-discharge curve —
    # a partly full reservoir has less room to attenuate the flood.
    if init_storage is not None and float(init_storage) > 0.0:
        s0 = float(init_storage)
        o0 = _lerp_curve(s0, S, O)
        a1 = 2.0 * s0 / dt_sec + o0
        o1 = o0
    else:
        a1, o1 = A[0], O[0]
    out[0] = o1
    stor[0] = max(0.0, (a1 - o1) * dt_sec / 2.0)
    ip = I[0] if m else 0.0
    max_stor = stor[0] if m else 0.0
    for k in range(1, m):
        if I[k] > ip:
            ip = I[k]
        a2 = (I[k - 1] + I[k]) + (a1 - 2.0 * o1)
        if a2 < 0:
            a2 = 0.0
        o2 = max(0.0, o_from_a(a2))
        out[k] = o2
        stor[k] = max(0.0, (a2 - o2) * dt_sec / 2.0)
        if stor[k] > max_stor:
            max_stor = stor[k]
        a1, o1 = a2, o2
    stage = max_stage = None
    if sd.get('has_stage') and sd.get('Saxis') and len(sd['Saxis']) >= 2:
        stage = np.array([_lerp_curve(s_, sd['Saxis'], sd['Eaxis'])
                          for s_ in stor])
        max_stage = float(stage.max()) if m else None
    max_o = O[-1]
    capped = bool(np.any(out > max_o * 1.001))
    return dict(O=out, storage=stor, stage=stage,
                max_storage=float(max_stor),
                max_stage=max_stage,
                inflow_peak=float(round(ip, 3)),
                warn=('reservoir outflow extrapolated beyond the '
                      'defined curve — extend the elevation–volume '
                      'curve / table' if capped else None),
                changed=True, info=sd.get('info'))


def load_dams(layer_name, name_field, attr_fields, dem_crs, project):
    """Load the reservoirs/dams point layer (v23).

    Returns (dams, out_attr_fields) where:
        dams            — list of dict(x, y, name, attrs={field: value})
                          with coordinates already in the DEM CRS, or
                          None when the feature is inactive/unusable.
        out_attr_fields — list of (out_column_name, QgsField) pairs for
                          the dam attribute columns to append to the
                          Junctions layer schema. Column names are the
                          source field names prefixed with 'DAM_' (the
                          name field itself is excluded — it already
                          lands in the dedicated DAM_NAME column).
    """
    if not layer_name:
        return None, []
    try:
        lyr = get_layer(layer_name)
    except Exception as e:
        print(f"⚠ Dams layer '{layer_name}' not found ({e}) — "
              f"Reservoirs feature disabled.")
        return None, []
    if not name_field or lyr.fields().indexOf(name_field) < 0:
        print(f"⚠ Dam name field '{name_field}' not found in layer "
              f"'{layer_name}' — Reservoirs feature disabled.")
        return None, []

    # ---- resolve which attribute fields to copy ------------------------
    wanted = [str(f) for f in (attr_fields or [])
              if f and str(f) != str(name_field)]
    out_attr_fields = []
    seen_out = set()
    for fname in wanted:
        idx = lyr.fields().indexOf(fname)
        if idx < 0:
            print(f"   ⚠ dam field '{fname}' not found — skipped.")
            continue
        src = lyr.fields().at(idx)
        out_name = f"DAM_{fname}"
        if out_name in seen_out:
            continue
        seen_out.add(out_name)
        fld = QgsField(src)
        fld.setName(out_name)
        out_attr_fields.append((fname, out_name, fld))

    xform = None
    if lyr.crs() != dem_crs:
        xform = QgsCoordinateTransform(lyr.crs(), dem_crs, project)

    dams = []
    for i, f in enumerate(lyr.getFeatures()):
        g = f.geometry()
        if g is None or g.isEmpty():
            continue
        pt = g.asPoint() if not g.isMultipart() else g.asMultiPoint()[0]
        if xform is not None:
            try:
                pt = xform.transform(pt)
            except Exception as e:
                print(f"   ⚠ dam #{i + 1} reprojection failed ({e}) — "
                      f"skipped.")
                continue
        nm = f[name_field]
        nm = (str(nm).strip() if nm is not None else '') or f"Dam_{i + 1}"
        attrs = {}
        for fname, out_name, _fld in out_attr_fields:
            v = f[fname]
            # QVariant NULL → None for clean attribute writing
            attrs[out_name] = None if (v is None or str(v) == 'NULL') else v
        dams.append(dict(x=float(pt.x()), y=float(pt.y()),
                         name=nm, attrs=attrs))
    if not dams:
        print(f"⚠ Dams layer '{layer_name}' has no usable point "
              f"features — Reservoirs feature disabled.")
        return None, []
    return dams, out_attr_fields


def get_layer(name):
    layers = QgsProject.instance().mapLayersByName(name)
    if not layers:
        raise ValueError(f"Layer '{name}' not in project")
    return layers[0]


def run(cfg=None, progress_cb=None):
    """Main entry point.

    Parameters
    ----------
    progress_cb : callable or None
        v1.0 — optional ``progress_cb(fraction, message)`` invoked at
        each major stage (fraction in [0, 1]). Used by the dialog to
        drive a progress bar and keep the UI responsive between
        stages. Exceptions raised by the callback are swallowed.
    cfg : dict or None
        Configuration dict with the same keys as the module-level
        ``CONFIG`` template at the top of this file. When ``None``
        (the default) the module-level ``CONFIG`` is used as-is, so
        the script remains runnable from the QGIS Python console.
    """
    if cfg is None:
        cfg = CONFIG
    project    = QgsProject.instance()
    dem_layer  = get_layer(cfg['DEM_LAYER'])
    pts_layer  = get_layer(cfg['POINTS_LAYER'])
    dem_crs    = dem_layer.crs()
    pts_crs    = pts_layer.crs()
    crs_authid = dem_crs.authid()
    # v21: all outputs are written into a dedicated subfolder named
    # 'ar_tarseem_output' inside the user-chosen directory, so runs
    # never scatter files across an existing folder. If the user
    # already pointed at a folder with that exact name, use it as-is
    # (avoids ar_tarseem_output/ar_tarseem_output nesting on re-runs).
    out_dir    = cfg['OUTPUT_DIR']
    if os.path.basename(os.path.normpath(out_dir)) != 'ar_tarseem_output':
        out_dir = os.path.join(out_dir, 'ar_tarseem_output')
    os.makedirs(out_dir, exist_ok=True)
    print(f"[CFG] Output folder: {out_dir}")

    # ---- v16: QGIS layer styling --------------------------------------
    # Resolved EARLY (before any output layer is created) because the
    # very first layers added to the project — DEM_filled, FlowAcc and
    # Streams — happen in steps [4] and [5] below, well before the
    # main per-pour-point loop. The on-screen [CFG] announcement is
    # emitted later, alongside the other CFG echoes, for grouping.
    styles_dir = str(cfg.get('STYLES_DIR', '') or '').strip()
    styles_dir = os.path.expanduser(styles_dir) if styles_dir else ''
    if styles_dir and not os.path.isdir(styles_dir):
        print(f"⚠ STYLES_DIR '{styles_dir}' not found or not a "
              f"directory — proceeding without styling.")
        styles_dir = ''

    # ------------------------------------------------------------------ CRS
    if dem_crs.isGeographic():
        raise ValueError(
            f"DEM CRS is geographic ({dem_crs.authid()}). "
            "Reproject the DEM to a metric CRS (e.g. UTM) before running."
        )

    print(f"[CRS] DEM         : {dem_crs.authid()} ({dem_crs.description()})")
    print(f"[CRS] Pour points : {pts_crs.authid()} ({pts_crs.description()})")

    if dem_crs.authid() != pts_crs.authid():
        print(f"[CRS] ⚠ Mismatch detected — pour points will be reprojected.")
        xform = QgsCoordinateTransform(pts_crs, dem_crs, project)
    else:
        print(f"[CRS] OK — same CRS, no reprojection needed.")
        xform = None

    # --------------------------------------------------------------- raster
    print("\n[0] Reading DEM ...")
    dem, gt, proj, nodata = read_dem(dem_layer.source())
    dx, dy = abs(gt[1]), abs(gt[5])
    print(f"    shape={dem.shape}  cell={dx:.2f} × {dy:.2f}")
    # extents in DEM CRS, for diagnostic
    xmin = gt[0]; ymax = gt[3]
    xmax = xmin + dem.shape[1] * gt[1]
    ymin = ymax + dem.shape[0] * gt[5]
    print(f"    extent  X:[{xmin:.1f}, {xmax:.1f}]   Y:[{ymin:.1f}, {ymax:.1f}]")

    # ---- v1.0: progress reporting ------------------------------------
    def _prog(frac, msg):
        if progress_cb is None:
            return
        try:
            progress_cb(max(0.0, min(1.0, float(frac))), str(msg))
        # B110: progress reporting must never abort the computation
        except Exception:  # nosec B110
            pass

    _prog(0.02, "Reading DEM & inputs")

    _prog(0.05, "Pit filling (priority flood)")
    print("[1] Pit-fill / depression breach (priority flood) ...")
    nd_mask = (dem == nodata) if nodata is not None else None
    dem_filled = priority_flood_fill(dem, nodata=nodata,
                                     epsilon=cfg['DEPRESSION_EPSILON'])
    if nd_mask is not None:
        dem_filled[nd_mask] = nodata

    _prog(0.14, "D8 flow direction")
    print("[2] D8 flow direction ...")
    fdir = flow_direction_d8(dem_filled, nodata_mask=nd_mask)

    _prog(0.22, "Flow accumulation")
    print("[3] Flow accumulation ...")
    facc = flow_accumulation(fdir)

    _prog(0.30, "Stream extraction")
    print("[4] Stream extraction ...")
    streams = (facc > cfg['STREAM_THRESHOLD']).astype(np.uint8)
    print(f"    stream cells: {int(streams.sum()):,}")

    print("[4b] Per-cell terrain slope (Horn 1981) ...")
    slope_source = str(cfg.get('BASIN_SLOPE_FROM', 'original')).lower()
    if slope_source not in ('original', 'filled'):
        print(f"    ⚠ BASIN_SLOPE_FROM='{slope_source}' not recognised; "
              f"falling back to 'original'.")
        slope_source = 'original'
    if slope_source == 'filled':
        print("    Using FILLED DEM (depressions raised to ε-graded plateaus).")
        slope_input = dem_filled
    else:
        print("    Using ORIGINAL DEM (depressions preserved at real elevations).")
        slope_input = dem
    slope_grid = horn_slope(slope_input, dx, dy, nodata_mask=nd_mask)
    valid = slope_grid[~np.isnan(slope_grid)]
    if valid.size:
        print(f"    slope range  : {valid.min():.4f} – {valid.max():.4f} m/m")
        print(f"    slope mean   : {valid.mean():.4f} m/m  ({valid.mean()*100:.2f} %)")

    # ---- v23: Reservoirs / dams — load + snap to main flow path ---------
    # Each dam point is snapped to the highest-flow-accumulation stream
    # cell within SNAP_RADIUS_CELLS (identical rule to pour-point
    # snapping), i.e. to the nearest MAIN flow path. The snapped cell is
    # later enforced as a junction during sub-catchment subdivision.
    dams_cfg_layer = str(cfg.get('DAMS_LAYER', '') or '').strip()
    dam_cells      = {}      # (row, col) -> dam dict (merged if shared)
    dam_attr_cols  = []      # [(src_field, out_name, QgsField), ...]
    res_any        = False   # v24: any point with a USABLE reservoir def
    if dams_cfg_layer:
        print(f"\n[4d] Special delineation points layer "
              f"'{dams_cfg_layer}' ...")
        _dams, dam_attr_cols = load_dams(
            dams_cfg_layer,
            str(cfg.get('DAM_NAME_FIELD', '') or '').strip(),
            cfg.get('DAM_ATTR_FIELDS') or [],
            dem_crs, project)
        sp_res_defs = cfg.get('SP_RESERVOIRS') or {}
        for d in (_dams or []):
            # v24 — attach the reservoir definition (if the user marked
            # this point as a reservoir) and prebuild its routing curve.
            _rdef = sp_res_defs.get(d['name'])
            d['res'] = _rdef
            d['sd_curve'] = None
            if _rdef:
                _sd = reservoir_sd(_rdef)
                if _sd and len(_sd.get('S') or []) >= 2:
                    d['sd_curve'] = _sd
                    res_any = True
                    print(f"     '{d['name']}': RESERVOIR — {_sd['info']}")
                    # v1.0 — initial condition
                    d['init_storage'] = None
                    if str(_rdef.get('init_mode', 'inflow')) == 'elevation':
                        _ie = _rdef.get('init_elev')
                        _s0 = (storage_at_elevation(_sd, _ie)
                               if _ie is not None else None)
                        if _s0 is None:
                            print(f"       ⚠ initial elevation needs an "
                                  f"elevation–volume curve; this "
                                  f"reservoir is defined by a storage–"
                                  f"discharge table, so the run starts "
                                  f"with outflow = inflow instead.")
                        else:
                            d['init_storage'] = _s0
                            _emin = _sd['Eaxis'][0]
                            _emax = _sd['Eaxis'][-1]
                            if not (_emin - 1e-9 <= float(_ie)
                                    <= _emax + 1e-9):
                                print(f"       ⚠ initial elevation "
                                      f"{float(_ie):g} m is outside the "
                                      f"elevation–volume curve "
                                      f"({_emin:g}–{_emax:g} m) — "
                                      f"clamped.")
                            print(f"       initial elevation "
                                  f"{float(_ie):g} m → starting storage "
                                  f"{_s0 / 1000.0:.1f}×10³ m³")
                else:
                    print(f"     ⚠ '{d['name']}': reservoir definition "
                          f"unusable ({(_sd or {}).get('warn')}) — the "
                          f"point is still enforced as a junction but "
                          f"flow passes through UNROUTED.")
            r0d, c0d = xy_to_cell(d['x'], d['y'], gt)
            if not (0 <= r0d < dem.shape[0] and 0 <= c0d < dem.shape[1]):
                print(f"     ⚠ dam '{d['name']}' outside DEM extent — "
                      f"skipped.")
                continue
            srd, scd, snapped_d = snap_to_stream(
                streams, facc, r0d, c0d, cfg['SNAP_RADIUS_CELLS'])
            if not snapped_d:
                print(f"     ⚠ dam '{d['name']}': no main flow path "
                      f"within {cfg['SNAP_RADIUS_CELLS']} cells — "
                      f"skipped (no junction enforced).")
                continue
            sxd, syd = cell_to_xy(srd, scd, gt)
            d['snap_xy'] = (sxd, syd)
            d['rc']      = (srd, scd)
            print(f"     dam '{d['name']}' snapped → cell "
                  f"({srd},{scd})  ({sxd:.1f}, {syd:.1f})  "
                  f"facc={int(facc[srd, scd])}")
            if (srd, scd) in dam_cells:
                # two points collapse onto the same stream cell — merge
                # names; the first point's attributes win, but a usable
                # reservoir definition survives from either point.
                prev = dam_cells[(srd, scd)]
                prev['name'] = f"{prev['name']}+{d['name']}"
                if prev.get('sd_curve') is None and d.get('sd_curve'):
                    prev['res'] = d['res']
                    prev['sd_curve'] = d['sd_curve']
                print(f"     ⚠ point '{d['name']}' shares a snapped cell "
                      f"with '{prev['name']}' — names merged, first "
                      f"point's attributes kept.")
            else:
                dam_cells[(srd, scd)] = d
        if dam_cells:
            print(f"     {len(dam_cells)} dam junction(s) will be "
                  f"enforced during subdivision.")

    # ---- optional Curve Number raster -----------------------------------
    cn_layer_name = cfg.get('CN_LAYER')
    cn_grid       = None
    cn_nodata     = None
    # v1.0 — constant Curve Number takes precedence over the raster.
    _cn_value = None
    try:
        if cfg.get('CN_VALUE') not in (None, ''):
            _cn_value = float(cfg.get('CN_VALUE'))
    except (TypeError, ValueError):
        _cn_value = None
    if _cn_value is not None and _cn_value > 0:
        if not (1.0 <= _cn_value <= 100.0):
            _cn_clamped = min(max(_cn_value, 1.0), 100.0)
            print(f"\n[4c] ⚠ CN_VALUE={_cn_value} outside 1–100 — "
                  f"clamped to {_cn_clamped}.")
            _cn_value = _cn_clamped
        print(f"\n[4c] Constant Curve Number CN = {_cn_value:g} applied "
              f"to the whole grid (no CN raster needed).")
        cn_grid   = np.full(dem.shape, _cn_value, dtype=np.float32)
        cn_nodata = None
        use_cn    = True
    else:
        use_cn = bool(cn_layer_name)
    if use_cn and cn_grid is None:
        print(f"\n[4c] Curve Number raster '{cn_layer_name}' ...")
        cn_layers = QgsProject.instance().mapLayersByName(cn_layer_name)
        if not cn_layers:
            print(f"    ⚠ CN layer '{cn_layer_name}' not found in project — "
                  f"skipping CN zonal statistics.")
            use_cn = False
        else:
            cn_layer = cn_layers[0]
            cn_src   = cn_layer.source()
            cn_crs   = cn_layer.crs().authid()
            print(f"    source: {cn_src}")
            print(f"    CRS: {cn_crs}  (DEM CRS: {dem_crs.authid()})")

            need_reproject = (cn_crs != dem_crs.authid()
                              or cn_layer.rasterUnitsPerPixelX() != dx
                              or cn_layer.rasterUnitsPerPixelY() != dy)
            if need_reproject:
                print("    reprojecting + resampling CN to match DEM grid "
                      "(nearest-neighbour) ...")
                cn_aligned = os.path.join(out_dir, 'CN_aligned.tif')
                reproject_to_match_dem(cn_src, gt, proj, dem.shape,
                                       cn_aligned, resample='near')
                cn_path = cn_aligned
            else:
                print("    CN raster already aligned with DEM grid.")
                cn_path = cn_src

            cn_arr, cn_gt, _, cn_nodata = read_dem(cn_path)
            if cn_arr.shape != dem.shape:
                # Pad / truncate as a last resort (shouldn't happen after warp)
                print(f"    ⚠ CN shape {cn_arr.shape} ≠ DEM shape {dem.shape}; "
                      f"forcing reprojection.")
                cn_aligned = os.path.join(out_dir, 'CN_aligned.tif')
                reproject_to_match_dem(cn_src, gt, proj, dem.shape,
                                       cn_aligned, resample='near')
                cn_arr, cn_gt, _, cn_nodata = read_dem(cn_aligned)

            cn_grid = cn_arr
            cn_valid = cn_grid.copy().astype(np.float64)
            if cn_nodata is not None:
                cn_valid[cn_grid == cn_nodata] = np.nan
            v = cn_valid[~np.isnan(cn_valid)]
            if v.size:
                print(f"    CN range : {v.min():.1f} – {v.max():.1f}  "
                      f"(mean {v.mean():.2f}, {v.size:,} valid cells)")
            else:
                print("    ⚠ CN raster has no valid cells inside DEM — "
                      "skipping zonal statistics.")
                use_cn = False
                cn_grid = None
    else:
        print("\n[4c] No Curve Number source (raster or constant) — "
              "skipping Curve Number step.")

    p_filled = write_raster(dem_filled, gt, proj,
                            os.path.join(out_dir, 'DEM_filled.tif'),
                            dtype=gdal.GDT_Float32, nodata=nodata)
    p_fdir   = write_raster(fdir.astype(np.int32), gt, proj,
                            os.path.join(out_dir, 'FlowDir.tif'),
                            dtype=gdal.GDT_Int32, nodata=0)
    p_facc   = write_raster(facc.astype(np.float32), gt, proj,
                            os.path.join(out_dir, 'FlowAcc.tif'),
                            dtype=gdal.GDT_Float32)
    for n, p in (('DEM_filled', p_filled), ('FlowAcc', p_facc)):
        rl = QgsRasterLayer(p, n)
        if rl.isValid():
            project.addMapLayer(rl)
            # v16 — apply QML style if one is registered for this raster
            style_raster_layer(rl, n, styles_dir, verbose=True)

    _prog(0.36, "Vectorising stream network")
    print("[5] Vectorising stream network ...")
    stream_lines = streams_to_polylines(fdir, streams, facc, gt)
    streams_out = QgsVectorLayer(f"LineString?crs={crs_authid}",
                                 cfg['OUT_STREAMS'], "memory")
    streams_out.dataProvider().addAttributes([QgsField('FID', _FT_INT)])
    streams_out.updateFields()
    sf = []
    for i, coords in enumerate(stream_lines):
        f = QgsFeature()
        f.setGeometry(QgsGeometry.fromPolylineXY(
            [QgsPointXY(x, y) for x, y in coords]))
        f.setAttributes([i])
        sf.append(f)
    streams_out.dataProvider().addFeatures(sf)
    streams_out.updateExtents()
    project.addMapLayer(streams_out)
    # v16 — apply QML style if one is registered
    style_vector_layer(streams_out, 'OUT_STREAMS', styles_dir, verbose=True)
    print(f"    polylines: {len(sf)}")

    # --------------------------------------------------- output layers ----
    # As of v11, Catchments and LFP_PL share an identical attribute schema
    # (the union of catchment metrics + LFP metrics) so either layer alone
    # is self-contained. Field order is locked so the same attribute list
    # populates either feature type.

    # ------- shared schema for Catchments AND LFP_PL ----------------------
    catch_lfp_fields = [
        # Identity / catchment metrics
        QgsField('PID',       _FT_INT),
        QgsField('NAME',      _FT_STRING, 'string', 32),
        QgsField('AREA_KM2',  _FT_DOUBLE, 'double', 15, 4),
        QgsField('BAS_SLP',   _FT_DOUBLE, 'double', 12, 6),
        QgsField('BAS_SLP_P', _FT_DOUBLE, 'double', 10, 3),
        # LFP geometry & elevations
        QgsField('L',         _FT_DOUBLE, 'double', 15, 2),
        QgsField('ELE_US',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('ELE_DS',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('ELE_10',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('ELE_85',    _FT_DOUBLE, 'double', 10, 2),
        # Legacy 10-85 slope
        QgsField('SLOPE',     _FT_DOUBLE, 'double', 10, 6),
        QgsField('SLOPE_P',   _FT_DOUBLE, 'double', 10, 2),
        # Segment slopes (m/m)
        QgsField('SLP_US_85', _FT_DOUBLE, 'double', 10, 6),
        QgsField('SLP_85_10', _FT_DOUBLE, 'double', 10, 6),
        QgsField('SLP_10_DS', _FT_DOUBLE, 'double', 10, 6),
        # Kirpich Tc per segment + total + Kirpich Tlag (minutes)
        QgsField('TC_US_85_K', _FT_DOUBLE, 'double', 12, 3),
        QgsField('TC_85_10_K', _FT_DOUBLE, 'double', 12, 3),
        QgsField('TC_10_DS_K', _FT_DOUBLE, 'double', 12, 3),
        QgsField('TC_KIRP',    _FT_DOUBLE, 'double', 12, 3),
        QgsField('TLAG_KIRP',  _FT_DOUBLE, 'double', 12, 3),
        # v14: Holistic Kirpich Tc / Tlag (whole-path, single segment)
        QgsField('TC_KIRPH',   _FT_DOUBLE, 'double', 12, 3),
        QgsField('TLAG_KIRPH', _FT_DOUBLE, 'double', 12, 3),
        # v1.0 — SCS unit-hydrograph peak rate factor actually used
        QgsField('PRF',        _FT_INT),
    ]
    if use_cn:
        catch_lfp_fields += [
            QgsField('CN_AVG',    _FT_DOUBLE, 'double', 8, 2),
            QgsField('CN_COV_P',  _FT_DOUBLE, 'double', 7, 2),
            QgsField('TLAG_NRCS', _FT_DOUBLE, 'double', 12, 3),
            QgsField('TC_NRCS',   _FT_DOUBLE, 'double', 12, 3),
        ]

    # Design-rainfall block (P / R / V per return period, in fixed order).
    # Only periods the user supplied are emitted. R and V require CN.
    active_periods = parse_design_rainfall(cfg.get('RAINFALL_DEPTHS_MM'))

    # ---- v22: spatially-variable rainfall (Thiessen stations) ----------
    thiessen_out = None          # v1.0 Voronoi polygons layer
    rain_stations, station_rps = load_rain_stations(
        str(cfg.get('RAIN_STATIONS_LAYER', '') or '').strip(),
        cfg.get('RAIN_STATIONS_FIELDS') or {},
        dem_crs, project)
    if rain_stations:
        # Return periods come from the field mapping; the uniform table
        # is superseded. The 'depth' put here is the plain station mean
        # — a placeholder only; every catchment / sub-basin gets its
        # own Thiessen-weighted depth downstream.
        active_periods = [
            (rp, float(np.mean([s['depths'][rp] for s in rain_stations])))
            for rp in station_rps]
        print(f"[CFG] Spatially-variable rainfall: "
              f"{len(rain_stations)} station(s), return periods "
              f"{station_rps} (Thiessen polygons on the DEM grid; "
              f"uniform RAINFALL_DEPTHS_MM table superseded)")
        for s in rain_stations:
            print(f"       {s['name']}: "
                  + ", ".join(f"{rp}yr={s['depths'][rp]:g}mm"
                              for rp in station_rps))
        # v1.0 — the Voronoi/Thiessen polygons actually used for the
        # weighting, as a map layer with each station's depths.
        thiessen_out = build_thiessen_layer(
            rain_stations, gt, dem.shape, crs_authid,
            layer_name=str(cfg.get('OUT_THIESSEN',
                                   'Rainfall_Thiessen_Polygons')
                           or 'Rainfall_Thiessen_Polygons'))
        if thiessen_out is not None:
            print(f"     Thiessen polygon layer built "
                  f"({thiessen_out.featureCount()} polygon(s))")

    # ---- v14: hydrograph engine settings -------------------------------
    # The peak-flow column QP_<rp>_yrs is added to Catchments and LFP_PL
    # only when:  HYDROGRAPHS_ENABLED  AND  use_cn  AND  active_periods.
    # Sub-basins do NOT receive QP_<rp>_yrs columns in v14 (reach
    # routing required — deferred to a later version).
    hyd_enabled_cfg = bool(cfg.get('HYDROGRAPHS_ENABLED', True))
    do_hydrographs  = bool(hyd_enabled_cfg and use_cn and active_periods)

    # Validate / clamp peak rate factor.
    prf_cfg = float(cfg.get('PEAK_RATE_FACTOR', 484.0))
    if prf_cfg < 100 or prf_cfg > 600:
        prf_clamped = max(100.0, min(600.0, prf_cfg))
        print(f"⚠ PEAK_RATE_FACTOR={prf_cfg} outside [100, 600] — "
              f"clamped to {prf_clamped}.")
        prf_cfg = prf_clamped

    tc_method_cfg = str(cfg.get('TC_METHOD', 'KIRPICH_HOLISTIC')).upper()
    if tc_method_cfg not in ('KIRPICH_HOLISTIC', 'KIRPICH_10_85', 'NRCS'):
        print(f"⚠ TC_METHOD='{cfg.get('TC_METHOD')}' not recognised — "
              f"defaulting to KIRPICH_HOLISTIC.")
        tc_method_cfg = 'KIRPICH_HOLISTIC'

    # ---- v15: ARF method ----------------------------------------------
    arf_method_cfg = str(cfg.get('ARF_METHOD', 'NONE')).upper().strip()
    if arf_method_cfg not in ('NONE', 'TRRL_UK', 'RIYADH', 'FIXED'):
        print(f"⚠ ARF_METHOD='{cfg.get('ARF_METHOD')}' not recognised — "
              f"defaulting to NONE (no ARF reduction).")
        arf_method_cfg = 'NONE'
    arf_fixed_cfg = float(cfg.get('ARF_FIXED_VALUE', 1.0) or 1.0)
    if arf_method_cfg == 'FIXED' and not (0.0 < arf_fixed_cfg <= 1.0):
        print(f"⚠ ARF_FIXED_VALUE={arf_fixed_cfg} outside (0, 1] — "
              f"using 1.0 (no reduction).")
        arf_fixed_cfg = 1.0
    arf_active = (arf_method_cfg != 'NONE')
    if arf_method_cfg == 'FIXED':
        print(f"   ARF: fixed user value {arf_fixed_cfg:g} applied to "
              f"all catchments")

    # ---- v17: reach routing & reach Tc method --------------------------
    reach_routing_cfg = str(
        cfg.get('REACH_ROUTING_METHOD', 'MUSKINGUM')).upper().strip()
    if reach_routing_cfg not in ('MUSKINGUM', 'MUSKINGUM_CUNGE',
                                 'LAG', 'NONE'):
        print(f"⚠ REACH_ROUTING_METHOD='{cfg.get('REACH_ROUTING_METHOD')}'"
              f" not recognised — defaulting to MUSKINGUM.")
        reach_routing_cfg = 'MUSKINGUM'

    # v1.0 — user minimums (0 = disabled)
    try:
        tc_floor_cfg = max(0.0, float(cfg.get('TC_MIN_MINUTES', 0.0) or 0.0))
    except (TypeError, ValueError):
        tc_floor_cfg = 0.0
    try:
        reach_tc_floor_cfg = max(0.0, float(
            cfg.get('REACH_TC_MIN_MINUTES', 0.0) or 0.0))
    except (TypeError, ValueError):
        reach_tc_floor_cfg = 0.0
    if tc_floor_cfg > 0:
        print(f"[CFG] minimum Tc = {tc_floor_cfg:g} min "
              f"(Tlag = 0.6 x Tc applied with it)")
    if reach_tc_floor_cfg > 0:
        print(f"[CFG] minimum reach travel time = "
              f"{reach_tc_floor_cfg:g} min")

    musk_x_cfg = float(cfg.get('MUSKINGUM_X', 0.2))
    # v1.0 — 0 = automatic per-reach subreaches; N>0 = fixed for all
    try:
        musk_nsub_cfg = max(0, int(cfg.get('MUSK_SUBREACHES', 0) or 0))
    except (TypeError, ValueError):
        musk_nsub_cfg = 0

    # ---- v19: lambda = Ia/S ratio ---------------------------
    scs_lambda_cfg = float(cfg.get('SCS_LAMBDA', 0.2) or 0.2)
    if not (0.01 <= scs_lambda_cfg <= 0.5):
        _cl = max(0.01, min(0.5, scs_lambda_cfg))
        print(f"\u26a0 SCS_LAMBDA={scs_lambda_cfg} outside [0.01, 0.5] "
              f"\u2014 clamped to {_cl}.")
        scs_lambda_cfg = _cl
    if scs_lambda_cfg != 0.2:
        print(f"   SCS-CN initial abstraction: Ia = "
              f"{scs_lambda_cfg:g}\u00b7S (non-classic lambda)")

    # ---- v19: Muskingum-Cunge channel geometry --------------
    mc_n_cfg = max(0.008, float(cfg.get('MC_MANNING_N', 0.035) or 0.035))
    mc_b_cfg = max(0.05,  float(cfg.get('MC_BOTTOM_WIDTH_M', 5.0) or 5.0))
    mc_z_cfg = max(0.0,   float(cfg.get('MC_SIDE_SLOPE', 2.0) or 2.0))
    if reach_routing_cfg == 'MUSKINGUM_CUNGE':
        print(f"   Muskingum-Cunge channel: n={mc_n_cfg:g}, "
              f"b={mc_b_cfg:g} m, z={mc_z_cfg:g} (H:V)")

    # ---- v19.1: custom hyetograph ----------------------------
    storm_dist_cfg = str(cfg.get('STORM_DISTRIBUTION',
                                 'SCS_TYPE_II')).upper().strip()
    custom_hyeto_cfg = None
    if storm_dist_cfg == 'CUSTOM':
        custom_hyeto_cfg = load_custom_hyetograph(
            str(cfg.get('CUSTOM_HYETOGRAPH_CSV', '') or '').strip())
        if custom_hyeto_cfg is None:
            print("\u26a0 CUSTOM storm selected but the hyetograph CSV "
                  "could not be used — falling back to SCS_TYPE_II.")

    # ---- v19: HEC-HMS export --------------------------------
    hms_export_cfg  = bool(cfg.get('EXPORT_HEC_HMS', False))
    # v1.0 — the network payload feeds BOTH the HEC-HMS files and the
    # AR-Awdiya project, so collect it when either export is requested.
    awdiya_export_cfg = bool(cfg.get('EXPORT_AWDIYA', False))
    net_collect_cfg = hms_export_cfg or awdiya_export_cfg
    hms_version_cfg = str(cfg.get('HEC_HMS_VERSION', '4.12') or '4.12')
    # v1.0 — sanitize the user-supplied HEC-HMS output name
    _hms_nm = str(cfg.get('HEC_HMS_NAME', '') or '').strip()
    hms_name_cfg = (re.sub(r'[^A-Za-z0-9._-]+', '_', _hms_nm).strip('_')
                    or 'AR_Tarseem')
    hms_catchments  = []   # per-catchment export payloads
    if not (0.0 <= musk_x_cfg <= 0.5):
        clamped = max(0.0, min(0.5, musk_x_cfg))
        print(f"⚠ MUSKINGUM_X={musk_x_cfg} outside [0, 0.5] — "
              f"clamped to {clamped}.")
        musk_x_cfg = clamped

    reach_tc_cfg = str(
        cfg.get('REACH_TC_METHOD', 'KIRPICH_HOLISTIC')).upper().strip()
    if reach_tc_cfg not in ('KIRPICH_HOLISTIC', 'KIRPICH_10_85'):
        print(f"⚠ REACH_TC_METHOD='{cfg.get('REACH_TC_METHOD')}' not "
              f"recognised — defaulting to KIRPICH_HOLISTIC.")
        reach_tc_cfg = 'KIRPICH_HOLISTIC'

    # ---- v6: sub-basin subdivision toggle -----------------------------
    # Resolved here (early) because the v17 [CFG] echoes below check
    # do_subdivide. The actual subdivision target_cells (which depends
    # on dx, dy from the DEM) is computed later in run(), once the
    # DEM has been loaded.
    target_area_km2 = float(cfg.get('SUB_BASIN_AREA_KM2', 0.0))
    do_subdivide    = target_area_km2 > 0
    # v23 — a dams layer FORCES subdivision so that a sub-catchment
    # boundary is always produced at every dam, even when the area-
    # based subdivision is off (SUB_BASIN_AREA_KM2 == 0). In that case
    # the junction set is simply {outlet} ∪ {dams in the catchment}.
    if dam_cells and not do_subdivide:
        do_subdivide = True
        print("\n[CFG] Reservoirs present — sub-catchment subdivision "
              "FORCED ON (junctions at dams + outlet only, since "
              "SUB_BASIN_AREA_KM2 = 0).")

    # ---- v16: QGIS layer styling — announcement only ------------------
    # styles_dir was already resolved at the top of run() (near out_dir);
    # here we simply log what will happen, grouped with the other CFG
    # echoes for readability.
    if styles_dir:
        print(f"\n[CFG] STYLES_DIR={styles_dir}  → QML styles will "
              f"be applied to outputs")
    else:
        print(f"\n[CFG] STYLES_DIR is empty — outputs will use QGIS "
              f"default styles")

    if active_periods:
        rp_summary = ", ".join(f"{rp}yr={p:g}mm" for rp, p in active_periods)
        if use_cn:
            if do_hydrographs:
                print(f"\n[CFG] design rainfall: {rp_summary}"
                      f"  → P/R/V/QP columns + hydrograph CSVs")
                print(f"[CFG] TC_METHOD={tc_method_cfg}  "
                      f"PEAK_RATE_FACTOR={prf_cfg:g}")
                # v17 routing settings (only meaningful if subdivision is on)
                if do_subdivide:
                    print(f"[CFG] REACH_ROUTING={reach_routing_cfg}"
                          f"  REACH_TC_METHOD={reach_tc_cfg}"
                          + (f"  X={musk_x_cfg}"
                             if reach_routing_cfg == 'MUSKINGUM' else ""))
            else:
                print(f"\n[CFG] design rainfall: {rp_summary}"
                      f"  → P/R/V columns (hydrographs disabled)")
        else:
            print(f"\n[CFG] design rainfall: {rp_summary}  → P columns only "
                  f"(no CN source set, so R and V cannot be computed)")
        # v15 — Areal reduction factor column (added BEFORE the rainfall
        # block so it appears immediately to the left of P_<rp>_yrs).
        if arf_active:
            print(f"[CFG] ARF_METHOD={arf_method_cfg}  "
                  f"→ P_<rp>_yrs values will be reduced; ARF column added")
            catch_lfp_fields.append(
                QgsField('ARF', _FT_DOUBLE, 'double', 8, 4))
        catch_lfp_fields += design_rainfall_field_blocks(
            active_periods, use_cn, with_qpeak=do_hydrographs)
    elif arf_active:
        # ARF requested but no rainfall to reduce — warn and continue.
        print(f"\n[CFG] ARF_METHOD={arf_method_cfg} requested but no "
              f"RAINFALL_DEPTHS_MM supplied; ARF will not be applied.")

    # Hydrograph output directory (only if hydrographs are enabled).
    # v18: structured Modeling_Outputs/ tree with one folder per element
    # type, each split into CSVs/ and HTMLs/.
    hyd_out_dir = None    # legacy alias retained: == modeling_dirs['ROOT']
    modeling_dirs = None
    if do_hydrographs:
        modeling_dirs = make_modeling_outputs_dirs(out_dir)
        hyd_out_dir = modeling_dirs['ROOT']

    catch_out = QgsVectorLayer(f"MultiPolygon?crs={crs_authid}",
                               cfg['OUT_CATCHMENTS'], "memory")
    catch_out.dataProvider().addAttributes(catch_lfp_fields)
    catch_out.updateFields()

    lfp_out = QgsVectorLayer(f"LineString?crs={crs_authid}",
                             cfg['OUT_LFP_LINES'], "memory")
    lfp_out.dataProvider().addAttributes(catch_lfp_fields)
    lfp_out.updateFields()

    pts_out = QgsVectorLayer(f"Point?crs={crs_authid}",
                             cfg['OUT_LFP_POINTS'], "memory")
    pts_out.dataProvider().addAttributes([
        QgsField('PID',        _FT_INT),
        QgsField('NAME',       _FT_STRING, 'string', 32),
        QgsField('POINT_TYPE', _FT_STRING),
        QgsField('ELEVATION',  _FT_DOUBLE, 'double', 10, 2),
        QgsField('DIST_FS',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('TOTAL_L',    _FT_DOUBLE, 'double', 15, 2),
        QgsField('SLOPE_P',    _FT_DOUBLE, 'double', 10, 2),
    ])
    pts_out.updateFields()

    snap_out = QgsVectorLayer(f"Point?crs={crs_authid}",
                              cfg['OUT_SNAPPED'], "memory")
    snap_out.dataProvider().addAttributes([
        QgsField('PID',     _FT_INT),
        QgsField('NAME',    _FT_STRING, 'string', 32),
        QgsField('SNAPPED', _FT_INT),
        QgsField('FACC',    _FT_LONGLONG),
    ])
    snap_out.updateFields()

    # ---- subdivision: target cell count (depends on DEM resolution) ----
    # do_subdivide and target_area_km2 were resolved early (alongside
    # other CFG flags). target_cells needs dx/dy from the DEM, so it
    # waits until here.
    # v24 — "special points ONLY": suppress the area-based subdivision
    # so the catchment is split exclusively at the outlet + the special
    # points (Stream threshold and Snap radius from the Delineation tab
    # still govern the network and the snapping).
    sp_points_only = bool(cfg.get('SP_POINTS_ONLY', False)) \
        and bool(dam_cells)
    target_cells = (int(round(target_area_km2 * 1e6 / (dx * dy)))
                    if (do_subdivide and not sp_points_only) else 0)
    if sp_points_only:
        print(f"\n[CFG] Special-points-ONLY subdivision — sub-basin "
              f"target area ignored; junctions = outlet + "
              f"{len(dam_cells)} special point(s).")
    elif do_subdivide:
        print(f"\n[CFG] sub-basin subdivision ENABLED — target = "
              f"{target_area_km2:.3f} km² ({target_cells:,} cells)")
    else:
        print(f"\n[CFG] sub-basin subdivision disabled "
              f"(SUB_BASIN_AREA_KM2 = 0).")

    # ------- shared schema for SubBasins AND SubLFPs ---------------------
    sub_sublfp_fields = [
        QgsField('PID',       _FT_INT),
        QgsField('SUB_ID',    _FT_INT),
        QgsField('DS_SUB_ID', _FT_INT),    # 0 if drains to pour-point
        QgsField('IS_OUTLET', _FT_INT),
        QgsField('NAME',      _FT_STRING, 'string', 32),
        QgsField('DS_NAME',   _FT_STRING, 'string', 32),
        # Sub-basin metrics
        QgsField('AREA_KM2',  _FT_DOUBLE, 'double', 15, 4),
        QgsField('BAS_SLP',   _FT_DOUBLE, 'double', 12, 6),
        QgsField('BAS_SLP_P', _FT_DOUBLE, 'double', 10, 3),
        # Sub-LFP geometry & elevations
        QgsField('L',         _FT_DOUBLE, 'double', 15, 2),
        QgsField('ELE_US',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('ELE_DS',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('ELE_10',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('ELE_85',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('SLOPE',     _FT_DOUBLE, 'double', 10, 6),
        QgsField('SLOPE_P',   _FT_DOUBLE, 'double', 10, 2),
        QgsField('SLP_US_85', _FT_DOUBLE, 'double', 10, 6),
        QgsField('SLP_85_10', _FT_DOUBLE, 'double', 10, 6),
        QgsField('SLP_10_DS', _FT_DOUBLE, 'double', 10, 6),
        QgsField('TC_US_85_K', _FT_DOUBLE, 'double', 12, 3),
        QgsField('TC_85_10_K', _FT_DOUBLE, 'double', 12, 3),
        QgsField('TC_10_DS_K', _FT_DOUBLE, 'double', 12, 3),
        QgsField('TC_KIRP',    _FT_DOUBLE, 'double', 12, 3),
        QgsField('TLAG_KIRP',  _FT_DOUBLE, 'double', 12, 3),
        # v14: Holistic Kirpich Tc / Tlag (whole-path, single segment)
        QgsField('TC_KIRPH',   _FT_DOUBLE, 'double', 12, 3),
        QgsField('TLAG_KIRPH', _FT_DOUBLE, 'double', 12, 3),
        # v1.0 — SCS unit-hydrograph peak rate factor actually used
        QgsField('PRF',        _FT_INT),
    ]
    if use_cn:
        sub_sublfp_fields += [
            QgsField('CN_AVG',    _FT_DOUBLE, 'double', 8, 2),
            QgsField('CN_COV_P',  _FT_DOUBLE, 'double', 7, 2),
            QgsField('TLAG_NRCS', _FT_DOUBLE, 'double', 12, 3),
            QgsField('TC_NRCS',   _FT_DOUBLE, 'double', 12, 3),
        ]
    # v15 — sub-basins inherit the parent catchment's ARF.
    if arf_active:
        sub_sublfp_fields.append(
            QgsField('ARF', _FT_DOUBLE, 'double', 8, 4))
    # v17 — design rainfall block on sub-basins, identical to catchments.
    # Each sub-basin gets P/R/V (and QP when hydrographs are enabled),
    # computed using its OWN local CN_AVG and area, but the parent's
    # ARF-reduced design depth.
    if active_periods:
        sub_sublfp_fields += design_rainfall_field_blocks(
            active_periods, use_cn, with_qpeak=do_hydrographs)

    sub_out = QgsVectorLayer(f"MultiPolygon?crs={crs_authid}",
                             cfg['OUT_SUBBASINS'], "memory")
    sub_out.dataProvider().addAttributes(sub_sublfp_fields)
    sub_out.updateFields()

    sublfp_out = QgsVectorLayer(f"LineString?crs={crs_authid}",
                                cfg['OUT_SUBLFPS'], "memory")
    sublfp_out.dataProvider().addAttributes(sub_sublfp_fields)
    sublfp_out.updateFields()

    reach_out = QgsVectorLayer(f"LineString?crs={crs_authid}",
                               cfg['OUT_REACHES'], "memory")
    reach_out.dataProvider().addAttributes([
        QgsField('PID',       _FT_INT),
        QgsField('SUB_ID',    _FT_INT),    # the upstream sub-basin
        QgsField('DS_SUB_ID', _FT_INT),    # the downstream sub-basin
        QgsField('NAME',      _FT_STRING, 'string', 32),
        QgsField('DS_NAME',   _FT_STRING, 'string', 32),
        # --- mirror LFP_PL schema (Tc and segment slopes) ---
        # Note: the 10/85 segmentation and Kirpich Tc are computed on the
        # reach treated as its own flow path. This matches HEC-HMS reach-
        # element timing parameters. TLAG_KIRP / TLAG_KIRPH are
        # intentionally NOT included here (reach-routing lag is not the
        # same concept as basin lag).
        QgsField('L',         _FT_DOUBLE, 'double', 15, 2),
        QgsField('ELE_US',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('ELE_DS',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('ELE_10',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('ELE_85',    _FT_DOUBLE, 'double', 10, 2),
        QgsField('SLOPE',     _FT_DOUBLE, 'double', 10, 6),
        QgsField('SLOPE_P',   _FT_DOUBLE, 'double', 10, 2),
        QgsField('SLP_US_85', _FT_DOUBLE, 'double', 10, 6),
        QgsField('SLP_85_10', _FT_DOUBLE, 'double', 10, 6),
        QgsField('SLP_10_DS', _FT_DOUBLE, 'double', 10, 6),
        QgsField('TC_US_85_K', _FT_DOUBLE, 'double', 12, 3),
        QgsField('TC_85_10_K', _FT_DOUBLE, 'double', 12, 3),
        QgsField('TC_10_DS_K', _FT_DOUBLE, 'double', 12, 3),
        QgsField('TC_KIRP',    _FT_DOUBLE, 'double', 12, 3),
        # v14: holistic Kirpich Tc on the reach
        QgsField('TC_KIRPH',   _FT_DOUBLE, 'double', 12, 3),
    ])
    # v17 — routing-method-tagged Tc (which one drives K) + per-RP routed
    # peak flow & % attenuation (only when hydrographs are enabled).
    reach_route_fields = [
        QgsField('RT_METHOD', _FT_STRING, 'string', 16),
        QgsField('RT_TC_MIN', _FT_DOUBLE, 'double', 12, 3),
        # v1.0 — Muskingum parameters actually used on this reach
        # (NULL when the routing method is not Muskingum):
        #   MUSK_K_HR — K in hours (= reach Tc), as exported to HMS
        #   SUBREACH  — subreach count of the routing cascade (fixed
        #               value from the dialog, or the automatic
        #               stability-based estimate)
        QgsField('MUSK_K_HR', _FT_DOUBLE, 'double', 12, 4),
        QgsField('SUBREACH',  _FT_INT),
    ]
    if active_periods and do_hydrographs:
        for rp, _ in active_periods:
            reach_route_fields.append(
                QgsField(f'QP_INF_{rp}_yrs',
                         _FT_DOUBLE, 'double', 12, 3))
            reach_route_fields.append(
                QgsField(f'QP_RTD_{rp}_yrs',
                         _FT_DOUBLE, 'double', 12, 3))
            reach_route_fields.append(
                QgsField(f'ATTEN_{rp}_pct',
                         _FT_DOUBLE, 'double', 8, 2))
    reach_out.dataProvider().addAttributes(reach_route_fields)
    reach_out.updateFields()

    # ---- v17: Junctions point layer ------------------------------------
    # One junction per sub-basin outlet. Same point is also the start of
    # that sub-basin's outflow reach AND the end of any upstream reach
    # flowing in. Naming: W{pid}-J{sid}.
    junc_out = QgsVectorLayer(f"Point?crs={crs_authid}",
                              cfg['OUT_JUNCTIONS'], "memory")
    junc_fields = [
        QgsField('PID',         _FT_INT),
        QgsField('SUB_ID',      _FT_INT),
        # v23 — widened to 80 chars: junction names now carry the dam
        # name suffix, e.g. "J1-Sub9(DamName)".
        QgsField('NAME',        _FT_STRING, 'string', 80),
        QgsField('DS_NAME',     _FT_STRING, 'string', 80),
        QgsField('IS_OUTLET',   _FT_INT),    # 1 for catchment outlet junction
        # v20.1: a single area column — the TRUE cumulative catchment
        # area upstream of the junction (this sub-basin's local area
        # plus every sub-basin that eventually drains to it, followed
        # through the ds_sid topology). No local/tributary split.
        QgsField('A_UPS_KM2',   _FT_DOUBLE, 'double', 15, 4),
    ]
    if arf_active:
        junc_fields.append(
            QgsField('ARF',     _FT_DOUBLE, 'double', 8, 4))
    if active_periods and do_hydrographs:
        # Combined peak Q at the junction = local + Σ(routed inflows),
        # plus (v20) the hydrograph volume at the junction, in 10³ m³
        # (same convention as V_<n>yrs_1000m3 on Catchments/SubBasins).
        for rp, _ in active_periods:
            junc_fields.append(
                QgsField(f'QP_{rp}_yrs',
                         _FT_DOUBLE, 'double', 12, 3))
            junc_fields.append(
                QgsField(f'V_{rp}yrs_1000m3',
                         _FT_DOUBLE, 'double', 15, 3))
    # v23 — reservoirs: dedicated dam-name column + one 'DAM_'-prefixed
    # copy of every user-selected dam attribute field. NULL on
    # junctions that are not dam locations.
    if dam_cells:
        junc_fields.append(
            QgsField('DAM_NAME', _FT_STRING, 'string', 64))
        for _src, _out, _fld in dam_attr_cols:
            junc_fields.append(QgsField(_fld))
    # v24 — reservoir routing outputs (only when at least one special
    # point has a usable reservoir definition). For a reservoir
    # junction, QP_<rp>_yrs / V_<rp>yrs_1000m3 above hold the routed
    # OUTFLOW peak/volume (what continues downstream); the columns
    # below add the AR-Awdiya-style detail. NULL on non-reservoir
    # junctions.
    #   RES_INFO            — routing-curve description
    #   QIN_<rp>_yrs        — peak reservoir inflow            [m³/s]
    #   ATT_<rp>_pct        — peak attenuation 100·(1−Qout/Qin)  [%]
    #   SMAX_<rp>_1000m3    — peak pool storage             [10³ m³]
    #   WLMX_<rp>_m         — peak water-surface elevation       [m]
    #                         (NULL for storage–discharge-table mode)
    if res_any:
        junc_fields.append(
            QgsField('RES_INFO', _FT_STRING, 'string', 160))
        if active_periods:
            for rp, _ in active_periods:
                junc_fields.append(
                    QgsField(f'QIN_{rp}_yrs',
                             _FT_DOUBLE, 'double', 12, 3))
                junc_fields.append(
                    QgsField(f'ATT_{rp}_pct',
                             _FT_DOUBLE, 'double', 8, 1))
                junc_fields.append(
                    QgsField(f'SMAX_{rp}_1000m3',
                             _FT_DOUBLE, 'double', 15, 3))
                junc_fields.append(
                    QgsField(f'WLMX_{rp}_m',
                             _FT_DOUBLE, 'double', 10, 2))
    junc_out.dataProvider().addAttributes(junc_fields)
    junc_out.updateFields()

    cf, lf, pf, snapf = [], [], [], []
    sub_feats, sublfp_feats, reach_feats, junc_feats = [], [], [], []

    # v23 — dam snap-verification points (PID 0 marks dam entries, so
    # they are distinguishable from pour-point snaps in the same layer).
    for _rc, _dam in dam_cells.items():
        _sx, _sy = _dam['snap_xy']
        _snf = QgsFeature()
        _snf.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(_sx, _sy)))
        _snf.setAttributes([0, f"Dam:{_dam['name']}", 1,
                            int(facc[_rc[0], _rc[1]])])
        snapf.append(_snf)

    # Union of all catchment masks — used at the end to write a CN raster
    # clipped to the area covered by the user's basins (only if CN is enabled).
    catch_union_mask = (np.zeros(dem.shape, dtype=bool)
                        if use_cn else None)
    id_field = cfg['POINT_ID_FIELD']

    _prog(0.40, "Processing pour points")
    print("\n[6] Processing pour points ...")
    _n_pts = max(1, int(pts_layer.featureCount() or 1))
    for i, feat in enumerate(pts_layer.getFeatures()):
        _prog(0.40 + 0.45 * (i / _n_pts),
              f"Catchment {i + 1} of {_n_pts} — delineation & "
              f"hydrographs")
        pid = feat[id_field] if id_field else (i + 1)
        gxy = feat.geometry().asPoint()
        x_native, y_native = float(gxy.x()), float(gxy.y())

        # ---- transform to DEM CRS if needed ---------------------------
        if xform is not None:
            try:
                p_xf = xform.transform(QgsPointXY(x_native, y_native))
                x0, y0 = float(p_xf.x()), float(p_xf.y())
                print(f"\n  -- Point {pid}  native ({x_native:.1f}, {y_native:.1f})"
                      f"  →  DEM-CRS ({x0:.1f}, {y0:.1f})")
            except Exception as e:
                print(f"\n  -- Point {pid}  reprojection FAILED: {e}")
                continue
        else:
            x0, y0 = x_native, y_native
            print(f"\n  -- Point {pid}  ({x0:.1f}, {y0:.1f})")

        r0, c0 = xy_to_cell(x0, y0, gt)
        if not (0 <= r0 < dem.shape[0] and 0 <= c0 < dem.shape[1]):
            print(f"     point outside DEM extent — skipping")
            continue

        sr, sc, snapped = snap_to_stream(streams, facc, r0, c0,
                                         cfg['SNAP_RADIUS_CELLS'])
        sx, sy = cell_to_xy(sr, sc, gt)
        if snapped:
            print(f"     snapped → cell ({sr},{sc})  ({sx:.1f}, {sy:.1f})  "
                  f"facc={int(facc[sr, sc])}")
        else:
            print(f"     ⚠ no stream within {cfg['SNAP_RADIUS_CELLS']} cells "
                  f"— using raw cell")

        # ---- Hierarchical names for this watershed -----------------------
        #   catch_name  : watershed identifier               e.g. "W1"
        #   lfp_name    : longest flow path identifier       e.g. "LFP1-W1"
        # The LFP number always equals the watershed (PID).
        catch_name = f"W{pid}"
        lfp_name   = f"LFP{pid}-{catch_name}"

        # save snap point for verification (tagged with catchment name)
        snf = QgsFeature()
        snf.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(sx, sy)))
        snf.setAttributes([int(pid), catch_name, int(snapped),
                           int(facc[sr, sc])])
        snapf.append(snf)

        # catchment
        catch = delineate_catchment(fdir, sr, sc)
        n_cells = int(catch.sum())
        if n_cells == 0:
            print(f"     empty catchment, skipping")
            continue
        area_km2 = n_cells * dx * dy / 1e6

        # mean basin slope from the pre-computed Horn slope grid
        bas_slp   = basin_mean_slope(slope_grid, catch)
        bas_slp_p = bas_slp * 100.0

        # CN zonal mean over the whole catchment (only if CN_LAYER set)
        if use_cn and cn_grid is not None:
            cn_avg, cn_cov = basin_mean_cn(cn_grid, cn_nodata, catch)
            cn_msg = f"  CN_avg={cn_avg:.2f} (cov={cn_cov:.1f}%)"
        else:
            cn_avg, cn_cov = None, None
            cn_msg = ""

        print(f"     catchment: {n_cells:,} cells  ({area_km2:.3f} km²)  "
              f"basin slope = {bas_slp:.4f} m/m  ({bas_slp_p:.2f} %){cn_msg}")

        # accumulate the catchment into the global union (for CN clipping)
        if catch_union_mask is not None:
            catch_union_mask |= catch

        # ---- LFP --------------------------------------------------------
        fdist = flow_distance(fdir, sr, sc, dx, dy)
        fdist_in = np.where(catch, fdist, np.nan)
        if np.all(np.isnan(fdist_in)):
            print(f"     could not compute LFP")
            continue
        max_idx = int(np.nanargmax(fdist_in))
        src_r, src_c = np.unravel_index(max_idx, fdist_in.shape)
        path_cells = trace_flow_path(fdir, src_r, src_c,
                                     stop_row=sr, stop_col=sc)
        lfp_coords = [cell_to_xy(rr, cc, gt) for rr, cc in path_cells]

        line_geom = QgsGeometry.fromPolylineXY(
            [QgsPointXY(x, y) for x, y in lfp_coords])
        ch = compute_lfp_characteristics(line_geom, dem_layer)
        if ch is None:
            continue
        print(f"     LFP L={ch['L']:.1f}m  US={ch['ELE_US']:.1f}  "
              f"DS={ch['ELE_DS']:.1f}  S(10-85)={ch['SLOPE_PCT']:.2f}%")
        print(f"     Kirpich: Tc(US-85)={ch['TC_US_85']:.2f}  "
              f"Tc(85-10)={ch['TC_85_10']:.2f}  "
              f"Tc(10-DS)={ch['TC_10_DS']:.2f}  "
              f"Tc_total={ch['TC_TOTAL']:.2f} min  "
              f"Tlag={ch['TLAG']:.2f} min")
        # holistic Kirpich (v14)
        print(f"     Kirpich (holistic): Tc={ch['TC_KIRPH']:.2f} min  "
              f"Tlag={ch['TLAG_KIRPH']:.2f} min")
        # warn on non-positive segment slopes (DEM noise → Tc=0 for that seg)
        for label, sval in (('SLP_US_85', ch['SLP_US_85']),
                            ('SLP_85_10', ch['SLP_85_10']),
                            ('SLP_10_DS', ch['SLP_10_DS'])):
            if sval <= 0:
                print(f"     ⚠ {label}={sval:.6f} (≤0) — Tc set to 0 for that segment")

        # ---- NRCS lag (only when CN is enabled) ------------------------
        tlag_nrcs = tc_nrcs = None
        if use_cn:
            tlag_nrcs = nrcs_lag(ch['L'], cn_avg, bas_slp_p)
            tc_nrcs   = tlag_nrcs / 0.6 if tlag_nrcs > 0 else 0.0
            print(f"     NRCS: CN={cn_avg:.2f}  Y={bas_slp_p:.2f}%  "
                  f"Tlag={tlag_nrcs:.2f} min  Tc={tc_nrcs:.2f} min")
            warn_msgs = []
            if cn_avg < 50 or cn_avg > 95:
                warn_msgs.append(f"CN={cn_avg:.1f} outside 50-95")
            if bas_slp_p > 5.0:
                warn_msgs.append(f"slope={bas_slp_p:.2f}% > 5%")
            if area_km2 > 10.0:
                warn_msgs.append(f"area={area_km2:.2f} km² > 10")
            if warn_msgs:
                print(f"     ⚠ NRCS validity: " + ", ".join(warn_msgs))

        # ---- v15: Areal Reduction Factor for this catchment -------------
        # ARF depends only on the FULL upstream catchment area at the
        # pour point. The same value will be inherited by every sub-
        # basin downstream of this pour point (sub-basins share the
        # design storm of their parent catchment).
        arf_value = compute_arf(area_km2, arf_method_cfg,
                                fixed_value=arf_fixed_cfg)
        if arf_active:
            print(f"     ARF ({arf_method_cfg}): A={area_km2:.2f} km² "
                  f"→ ARF={arf_value:.4f}")
            if arf_value < 1.0:
                print(f"       all P_<rp>_yrs values will be reduced by "
                      f"this factor before any downstream calculation")

        # Pre-compute the ARF-reduced rainfall pairs once. These drive
        # both the design_rainfall_attr_block call AND the hydrograph
        # engine, so every column in the resulting row is internally
        # self-consistent (P_table → R via SCS-CN → V → Qpeak).
        # v22: in Thiessen mode every catchment gets its own area-
        # weighted depth from the station Voronoi assignment; the ARF
        # then applies on top, exactly as for uniform rainfall.
        catch_thiessen = None
        if rain_stations:
            catch_thiessen, catch_weights = thiessen_basin_depths(
                catch, gt, rain_stations)
        if catch_thiessen:
            active_periods_eff = [
                (rp, catch_thiessen.get(rp, p_mm) * arf_value)
                for rp, p_mm in active_periods]
            w_str = ", ".join(f"{k}:{v * 100:.1f}%"
                              for k, v in sorted(catch_weights.items(),
                                                 key=lambda x: -x[1]))
            print(f"     Thiessen weights: {w_str}")
            print("     P (Thiessen): "
                  + ", ".join(f"{rp}yr={p / max(arf_value, 1e-9):.1f}mm"
                              for rp, p in active_periods_eff))
        else:
            active_periods_eff = (
                [(rp, p_mm * arf_value) for rp, p_mm in active_periods]
                if active_periods else []
            )

        # ---- v14: design-storm hydrographs per return period -----------
        # Computed only when:  HYDROGRAPHS_ENABLED, CN available, and at
        # least one rainfall return period was supplied. The selected
        # Tlag (TC_METHOD) and PRF (PEAK_RATE_FACTOR) drive the SCS UH.
        # In v15 the rainfall depths fed in are ALREADY ARF-reduced.
        qpeaks_by_rp = {} if do_hydrographs else None
        if do_hydrographs:
            tlag_sel, tc_sel, label_sel = select_tlag_for_method(
                tc_method_cfg, ch, tlag_nrcs, tc_floor_cfg)
            if tlag_sel <= 0:
                print(f"     ⚠ selected Tlag is non-positive "
                      f"(method={label_sel}); hydrographs skipped for {catch_name}.")
            else:
                print(f"     Hydrograph engine: Tlag={tlag_sel:.2f} min "
                      f"(method={label_sel}), PRF={prf_cfg:g}")
                for rp, p_eff in active_periods_eff:
                    hyd = compute_design_hydrograph(
                        p_total_mm = p_eff,
                        cn         = cn_avg,
                        area_km2   = area_km2,
                        tp_min     = tlag_sel,
                        prf        = prf_cfg,
                        storm_duration_hr  = float(cfg.get('STORM_DURATION_HR', 24.0)),
                        dt_min             = float(cfg.get('HYD_TIME_STEP_MIN', 1.0)),
                        storm_distribution = str(cfg.get('STORM_DISTRIBUTION',
                                                            'SCS_TYPE_II')),
                        scs_lambda         = scs_lambda_cfg,
                        custom_hyetograph  = custom_hyeto_cfg,
                        tail_factor        = float(cfg.get('HYD_TAIL_FACTOR', 20.0)),
                        max_tail_hr        = float(cfg.get('HYD_MAX_TAIL_HR', 96.0)),
                    )
                    if hyd is None:
                        qpeaks_by_rp[rp] = None
                        continue
                    qpeaks_by_rp[rp] = hyd['peak_Q']
                    # v1.0 — completeness check (recession to ~zero)
                    _qe = float(hyd['Q_m3s'][-1])
                    if hyd['peak_Q'] > 0 and _qe > 0.005 * hyd['peak_Q']:
                        print(f"       ⚠ {catch_name} {rp}-yr: "
                              f"hydrograph NOT fully receded at the "
                              f"end of the output window (Q_end = "
                              f"{_qe:.3f} m³/s = "
                              f"{100.0*_qe/hyd['peak_Q']:.1f}% of "
                              f"peak). Increase 'Hydrograph tail "
                              f"length' / 'Max tail' in the dialog.")
                    csv_name = f"{catch_name}_{rp}yr.csv"
                    csv_path = os.path.join(
                        modeling_dirs['CATCHMENT']['csv'], csv_name)
                    write_hydrograph_csv(csv_path, hyd)
                    # v18 — interactive HTML next to the CSV
                    html_path = os.path.join(
                        modeling_dirs['CATCHMENT']['html'],
                        csv_name[:-4] + '.html')
                    write_hydrograph_html(
                        html_path, kind='CATCHMENT',
                        title=f"{catch_name} — Catchment hydrograph "
                              f"({rp}-yr storm)",
                        subtitle=f"P={p_eff:.1f} mm  •  CN={cn_avg:.1f}  •  "
                                 f"A={area_km2:.2f} km²  •  "
                                 f"Tlag={tlag_sel:.0f} min "
                                 f"({label_sel})",
                        hyd_data={
                            'time_min':  list(hyd['time_min']),
                            'P_mm':      list(hyd['P_mm']),
                            'Loss_mm':   list(hyd['Loss_mm']),
                            'Excess_mm': list(hyd['Excess_mm']),
                            'Q_m3s':     list(hyd['Q_m3s']),
                        },
                    )
                    print(f"       {rp}-yr  P={p_eff:.1f}mm "
                          f"{'(ARF-reduced)' if arf_active and arf_value < 1.0 else ''}"
                          f"  R={hyd['R_total']:.2f}mm  "
                          f"Qp={hyd['peak_Q']:.2f} m³/s @ "
                          f"t={hyd['peak_t']:.0f} min  →  CSV+HTML")

        # ---- BUILD COMBINED ATTRIBUTE LIST (Catchments + LFP_PL share) ---
        # Order MUST match catch_lfp_fields exactly.
        # NAME (index 1) is set per-layer below, since the polygon and line
        # use different naming conventions ("W{pid}" vs "LFP{pid}-W{pid}").
        combined_attrs = [
            int(pid), None,        # NAME placeholder; overwritten per layer
            float(area_km2), float(bas_slp), float(bas_slp_p),
            ch['L'],
            ch['ELE_US'], ch['ELE_DS'], ch['ELE_10'], ch['ELE_85'],
            ch['SLOPE'], ch['SLOPE_PCT'],
            ch['SLP_US_85'], ch['SLP_85_10'], ch['SLP_10_DS'],
            ch['TC_US_85'], ch['TC_85_10'], ch['TC_10_DS'],
            ch['TC_TOTAL'], ch['TLAG'],
            # v14 — holistic Kirpich
            float(ch['TC_KIRPH']), float(ch['TLAG_KIRPH']),
            int(round(prf_cfg)),        # v1.0 — adopted PRF
        ]
        if use_cn:
            combined_attrs += [
                float(cn_avg), float(cn_cov),
                float(tlag_nrcs), float(tc_nrcs),
            ]
        # v15 — ARF column (sits immediately before the rainfall block)
        if arf_active:
            combined_attrs.append(float(arf_value))

        # Design-rainfall block (P / R / V [/ QP]) — order must match the schema
        # In v15 we pass the ARF-reduced pairs so P_<rp>_yrs in the
        # attribute table reflects the depth that was actually used.
        if active_periods:
            combined_attrs += design_rainfall_attr_block(
                active_periods_eff,
                cn_avg if use_cn else None,
                area_km2,
                use_cn,
                qpeaks=qpeaks_by_rp,
                scs_lambda=scs_lambda_cfg)
            if use_cn:
                # brief console echo (uses the ARF-reduced depths)
                rs = ", ".join(
                    f"{rp}yr P={p:.1f}mm R={scs_cn_runoff(p, cn_avg):.2f}mm "
                    f"V={scs_cn_runoff(p, cn_avg)*area_km2:.2f}×10³m³"
                    for rp, p in active_periods_eff)
                print(f"     SCS-CN: {rs}")

        # ---- emit Catchments feature ------------------------------------
        poly_geom = mask_to_qgs_geometry(catch, gt, proj)
        if poly_geom is not None:
            attrs_for_catch = list(combined_attrs)
            attrs_for_catch[1] = catch_name
            f = QgsFeature()
            f.setGeometry(poly_geom)
            f.setAttributes(attrs_for_catch)
            cf.append(f)

        # ---- emit LFP_PL feature (same attributes, line geometry, LFP name)
        attrs_for_lfp = list(combined_attrs)
        attrs_for_lfp[1] = lfp_name
        fl = QgsFeature()
        fl.setGeometry(line_geom)
        fl.setAttributes(attrs_for_lfp)
        lf.append(fl)

        # ---- emit LFP_PTs (4 points: START, END, P10, P85) --------------
        for ptype, pxy, pele, pdist in [
            ('START', ch['PT_US'], ch['ELE_US'], 0.0),
            ('END',   ch['PT_DS'], ch['ELE_DS'], ch['L']),
            ('P10',   ch['PT_10'], ch['ELE_10'], ch['D10']),
            ('P85',   ch['PT_85'], ch['ELE_85'], ch['D85']),
        ]:
            if pxy is None:
                continue
            pp = QgsFeature()
            pp.setGeometry(QgsGeometry.fromPointXY(pxy))
            pp.setAttributes([int(pid), lfp_name, ptype, float(pele),
                              float(pdist), ch['L'], ch['SLOPE_PCT']])
            pf.append(pp)

        # ============================================================
        #   SUB-BASIN SUBDIVISION (only if SUB_BASIN_AREA_KM2 > 0)
        # ============================================================
        if not do_subdivide:
            continue

        outlet_rc = (sr, sc)
        junctions = place_junctions(catch, streams, fdir, facc,
                                    outlet_rc, target_cells)
        junctions_set = set(junctions)

        # ---- v23: enforce a junction at every dam inside this ----------
        # catchment. The dam cell is already snapped onto the main flow
        # path (highest-facc stream cell within the Snap radius), so it
        # is guaranteed to sit on the in-basin stream network whenever
        # the catchment mask covers it.
        dam_of_cell = {}          # (r,c) -> dam dict, THIS catchment only
        n_dam_junc  = 0
        for _rc, _dam in dam_cells.items():
            if not catch[_rc[0], _rc[1]]:
                continue          # dam belongs to another catchment
            dam_of_cell[_rc] = _dam
            if _rc not in junctions_set:
                junctions.append(_rc)
                junctions_set.add(_rc)
                n_dam_junc += 1
                print(f"     dam '{_dam['name']}' → junction ENFORCED "
                      f"at cell {_rc}")
            else:
                print(f"     dam '{_dam['name']}' coincides with an "
                      f"existing junction at cell {_rc} — junction "
                      f"tagged with the dam.")

        n_junc = len(junctions)
        print(f"     subdivision: {n_junc} junction(s) "
              f"({n_junc - 1 - n_dam_junc} confluences/inserts + "
              f"{n_dam_junc} dam(s) + 1 outlet)")

        # topology: who flows into whom + reach geometries
        ds_junction, reach_cells, reach_coords = build_subbasin_topology(
            junctions, junctions_set, fdir, outlet_rc, gt)

        # number sub-basins so SUB_ID 1 == outlet sub-basin
        sid_of, ds_sid_of = renumber_junctions_outlet_first(
            junctions, ds_junction, outlet_rc)

        # v23 — dam lookup by SUB_ID (a dam junction IS the outlet of
        # the sub-basin it creates, so sid_of[dam_cell] identifies the
        # junction feature that must carry the dam name + attributes).
        dam_of_sid = {sid_of[_rc]: _dam
                      for _rc, _dam in dam_of_cell.items()
                      if _rc in sid_of}

        # which cells belong to which sub-basin
        sid_grid = compute_subbasin_id_grid(
            junctions, junctions_set, fdir, catch)
        # remap arbitrary internal IDs → outlet-first IDs
        old_to_new = np.zeros(n_junc + 2, dtype=np.int32)
        for idx, j in enumerate(junctions):
            old_id = idx + 1
            old_to_new[old_id] = sid_of[j]
        sid_grid_renumbered = old_to_new[sid_grid]

        # ============================================================
        # v17 — Three-pass sub-basin / reach / junction processing.
        #
        # Pass 0 (v17.1): compute a UNIFORM length budget for every
        #         hydrograph in this catchment, based on the largest
        #         sub-basin Tlag plus the slowest reach Tc. All
        #         sub-basin hydrographs and routed reach outflows in
        #         this catchment will share that length. This is what
        #         makes element-wise summation in pass 3 safe.
        #
        # Pass 1: per-sub-basin geometry, local CN, local hydrograph.
        #         Stored per SUB_ID in `sub_data[sid]`. Schema attrs
        #         not yet appended (P/R/V/QP block needs the local
        #         hydrograph result that we compute here).
        #
        # Pass 2: reach geometry per non-outlet sub-basin, plus the
        #         routing-method metadata (Tc + method label). The
        #         actual routing is deferred to pass 3 so we can use
        #         topological order.
        #
        # Pass 3: walk SUB_IDs in upstream→downstream topological
        #         order. For each junction:
        #           Q_combined(t) = Q_local(t) + Σ Q_routed(t) from
        #                           all upstream sub-basins draining
        #                           directly to this junction.
        #         If this is NOT the outlet sub-basin, route the
        #         combined Q through this sub-basin's reach to obtain
        #         Q_routed for delivery to the downstream junction.
        #         Emit SubBasins, SubLFPs, Reaches and Junctions
        #         features once everything is known.
        # ============================================================

        # ----- Pass 0 (v17.1): compute uniform hydrograph length -----
        # Without this, compute_design_hydrograph picks total_minutes
        # individually per sub-basin (storm + 20·Tp·stretch), so two
        # sub-basins with different Tlags end up with arrays of
        # different shape and Q_local + Q_routed broadcast-fails.
        # We size every hydrograph in this catchment to the worst-case
        # tail seen by ANY sub-basin or reach.
        catch_n_total = None
        if do_hydrographs and use_cn and active_periods:
            storm_min = int(round(
                float(cfg.get('STORM_DURATION_HR', 24.0)) * 60.0))
            dt_min_cfg = float(cfg.get('HYD_TIME_STEP_MIN', 1.0))
            stretch = max(1.0, _PRF_STD / max(float(prf_cfg), 1.0))

            # Pre-compute every sub-basin's selected Tlag and every
            # reach's Tc so we can pick the maximum tail.
            tlag_max = 0.0
            for j in junctions:
                this_sid = sid_of[j]
                sub_mask_pre = (sid_grid_renumbered == this_sid)
                if int(sub_mask_pre.sum()) == 0:
                    continue
                # Use the sub-basin LFP characteristics for the lag
                _line_pre, ch_pre = lfp_within_mask(
                    sub_mask_pre, j, fdir, dx, dy, gt, dem_layer)
                if ch_pre is None:
                    continue
                if use_cn and cn_grid is not None:
                    cn_pre, _ = basin_mean_cn(cn_grid, cn_nodata,
                                              sub_mask_pre)
                else:
                    cn_pre = None
                if cn_pre is None:
                    continue
                bslp_pre_p = basin_mean_slope(slope_grid,
                                              sub_mask_pre) * 100.0
                tlag_nrcs_pre = (nrcs_lag(ch_pre['L'], cn_pre, bslp_pre_p)
                                 if use_cn else None)
                tlag_pre, _, _ = select_tlag_for_method(
                    tc_method_cfg, ch_pre, tlag_nrcs_pre, tc_floor_cfg)
                if tlag_pre and tlag_pre > tlag_max:
                    tlag_max = tlag_pre

            # Reach Tc contribution to tail: the slowest reach adds
            # roughly its own travel time to the recession length.
            tc_max = 0.0
            for j in junctions:
                this_sid = sid_of[j]
                if this_sid == 1:
                    continue
                coords_pre = reach_coords.get(j, [])
                if len(coords_pre) < 2:
                    continue
                rg_pre = QgsGeometry.fromPolylineXY(
                    [QgsPointXY(x, y) for x, y in coords_pre])
                rch_pre = compute_lfp_characteristics(rg_pre, dem_layer)
                if rch_pre is None:
                    continue
                tc_min_pre, _ = select_reach_tc(
                    rch_pre, reach_tc_cfg, reach_tc_floor_cfg)
                if tc_min_pre and tc_min_pre > tc_max:
                    tc_max = tc_min_pre

            # Tail = storm + factor·Tp·stretch + reach lag, capped at
            # storm + HYD_MAX_TAIL_HR (both editable, v1.0)
            _tailf = float(cfg.get('HYD_TAIL_FACTOR', 20.0))
            _tailcap = float(cfg.get('HYD_MAX_TAIL_HR', 96.0))
            tlag_for_tail = max(tlag_max, 1.0)   # avoid zero
            total_minutes_budget = int(round(
                storm_min + _tailf * tlag_for_tail * stretch + tc_max))
            total_minutes_budget = min(total_minutes_budget,
                                       int(round(storm_min
                                                 + 60.0 * _tailcap)))
            catch_n_total = int(round(
                total_minutes_budget / dt_min_cfg)) + 1
            print(f"     v17 length budget: tlag_max={tlag_max:.1f}  "
                  f"tc_max={tc_max:.1f}  → {catch_n_total} steps "
                  f"({catch_n_total*dt_min_cfg/60:.1f} h)")

        # Helper: align a 1-D numpy array to the catchment-uniform
        # length (pad with zeros, trim if longer).
        def _align_len(arr, n):
            arr = np.asarray(arr, dtype=float)
            if len(arr) == n:
                return arr
            if len(arr) > n:
                return arr[:n]
            out = np.zeros(n, dtype=float)
            out[:len(arr)] = arr
            return out

        # ----- Pass 1: per-sub-basin geometry + local hydrograph -----
        sub_data = {}      # sid -> dict(name, area_km2, mask, geom,
                           #             cn, line, ch, lfp_name, ds_*,
                           #             local_hyds={rp: hyd dict},
                           #             local_qpeaks={rp: float},
                           #             tlag_used)
        for j in junctions:
            this_sid = sid_of[j]
            this_ds  = ds_sid_of[this_sid]
            sub_mask = (sid_grid_renumbered == this_sid)
            n_sub_cells = int(sub_mask.sum())
            if n_sub_cells == 0:
                continue
            sub_area = n_sub_cells * dx * dy / 1e6
            sub_bslp = basin_mean_slope(slope_grid, sub_mask)
            sub_bslp_p = sub_bslp * 100.0

            sub_name    = f"W{pid}-Sub{this_sid}"
            sublfp_name = f"LFP{pid}-{sub_name}"
            if this_ds > 0:
                ds_sub_name    = f"W{pid}-Sub{this_ds}"
                ds_sublfp_name = f"LFP{pid}-{ds_sub_name}"
            else:
                ds_sub_name    = catch_name      # outlet sub → catchment
                ds_sublfp_name = lfp_name        # outlet sub-LFP → parent LFP

            if use_cn and cn_grid is not None:
                sub_cn, sub_cn_cov = basin_mean_cn(cn_grid, cn_nodata, sub_mask)
            else:
                sub_cn, sub_cn_cov = None, None

            sub_line, sub_ch = lfp_within_mask(
                sub_mask, j, fdir, dx, dy, gt, dem_layer)
            if sub_line is None or sub_ch is None:
                continue

            sub_tlag_nrcs = sub_tc_nrcs = None
            if use_cn:
                sub_tlag_nrcs = nrcs_lag(sub_ch['L'], sub_cn, sub_bslp_p)
                sub_tc_nrcs   = (sub_tlag_nrcs / 0.6
                                 if sub_tlag_nrcs > 0 else 0.0)

            # ---- v17: per-sub-basin LOCAL hydrograph -----------------
            # Driver: parent catchment's ARF-reduced design depth (each
            # sub-basin shares the same storm), but with sub-basin CN,
            # sub-basin area, and sub-basin Tlag (per TC_METHOD).
            # v22: per-sub-basin Thiessen depth — this is the point of
            # the spatially-variable mode: each sub-basin's storms use
            # ITS OWN area-weighted depth, not the catchment mean.
            sub_periods_eff = active_periods_eff
            if rain_stations:
                _sub_th, _sub_w = thiessen_basin_depths(
                    sub_mask, gt, rain_stations)
                if _sub_th:
                    sub_periods_eff = [
                        (rp, _sub_th.get(rp,
                                         p / max(arf_value, 1e-9))
                         * arf_value)
                        for rp, p in active_periods_eff]
            local_hyds   = {}
            local_qpeaks = {}
            sub_tlag_sel = sub_tc_label = None
            if do_hydrographs and use_cn and (sub_cn is not None):
                sub_tlag_sel, _, sub_tc_label = select_tlag_for_method(
                    tc_method_cfg, sub_ch, sub_tlag_nrcs, tc_floor_cfg)
                if sub_tlag_sel and sub_tlag_sel > 0:
                    for rp, p_eff in sub_periods_eff:
                        h = compute_design_hydrograph(
                            p_total_mm = p_eff,
                            cn         = sub_cn,
                            area_km2   = sub_area,
                            tp_min     = sub_tlag_sel,
                            prf        = prf_cfg,
                            storm_duration_hr = float(
                                cfg.get('STORM_DURATION_HR', 24.0)),
                            dt_min            = float(
                                cfg.get('HYD_TIME_STEP_MIN', 1.0)),
                            storm_distribution = str(
                                cfg.get('STORM_DISTRIBUTION',
                                           'SCS_TYPE_II')),
                            scs_lambda         = scs_lambda_cfg,
                            custom_hyetograph  = custom_hyeto_cfg,
                            # v17.1 — pin the length to the catchment-
                            # wide budget so all sub-basin Q arrays in
                            # this catchment can be summed elementwise.
                            total_minutes = (catch_n_total - 1
                                             if catch_n_total else None),
                        )
                        if h is None:
                            local_qpeaks[rp] = None
                            continue
                        local_hyds[rp]   = h
                        local_qpeaks[rp] = h['peak_Q']
                        # write per-sub-basin local hydrograph CSV + HTML (v18)
                        csv_name = f"{sub_name}_{rp}yr.csv"
                        write_hydrograph_csv(
                            os.path.join(modeling_dirs['SUBBASIN']['csv'],
                                         csv_name), h)
                        html_path = os.path.join(
                            modeling_dirs['SUBBASIN']['html'],
                            csv_name[:-4] + '.html')
                        write_hydrograph_html(
                            html_path, kind='SUBBASIN',
                            title=f"{sub_name} — Sub-basin local hydrograph "
                                  f"({rp}-yr storm)",
                            subtitle=f"P={p_eff:.1f} mm  •  "
                                     f"CN={sub_cn:.1f}  •  "
                                     f"A={sub_area:.2f} km²  •  "
                                     f"Tlag={sub_tlag_sel:.0f} min "
                                     f"({sub_tc_label})",
                            hyd_data={
                                'time_min':  list(h['time_min']),
                                'P_mm':      list(h['P_mm']),
                                'Loss_mm':   list(h['Loss_mm']),
                                'Excess_mm': list(h['Excess_mm']),
                                'Q_m3s':     list(h['Q_m3s']),
                            },
                        )

            sub_data[this_sid] = dict(
                sid          = this_sid,
                ds_sid       = this_ds,
                is_outlet    = (1 if this_sid == 1 else 0),
                name         = sub_name,
                ds_name      = ds_sub_name,
                lfp_name     = sublfp_name,
                ds_lfp_name  = ds_sublfp_name,
                area_km2     = sub_area,
                bslp         = sub_bslp,
                bslp_p       = sub_bslp_p,
                mask         = sub_mask,
                line         = sub_line,
                ch           = sub_ch,
                cn           = sub_cn,
                cn_cov       = sub_cn_cov,
                tlag_nrcs    = sub_tlag_nrcs,
                tc_nrcs      = sub_tc_nrcs,
                tlag_used    = sub_tlag_sel,
                tc_label     = sub_tc_label,
                local_hyds   = local_hyds,
                local_qpeaks = local_qpeaks,
                periods_eff  = sub_periods_eff,
                junction_xy  = cell_to_xy(j[0], j[1], gt),
                dam          = dam_of_sid.get(this_sid),   # v23
            )

        # ----- Pass 2: reach geometry + Tc selection -----------------
        # reach_data[sid] only exists for non-outlet sub-basins.
        reach_data = {}
        for j in junctions:
            this_sid = sid_of[j]
            if this_sid == 1:
                continue                           # outlet has no reach
            coords = reach_coords.get(j, [])
            if len(coords) < 2:
                continue
            rg = QgsGeometry.fromPolylineXY(
                [QgsPointXY(x, y) for x, y in coords])
            rch = compute_lfp_characteristics(rg, dem_layer)
            if rch is None:
                continue
            this_ds_sid = ds_sid_of[this_sid]
            tc_min, tc_label = select_reach_tc(
                rch, reach_tc_cfg, reach_tc_floor_cfg)
            # v1.0 — reach naming W{pid}-R{sid}; the downstream
            # reference is the JUNCTION the reach discharges into
            # (dam-name suffix included when that junction is a
            # reservoir), matching the HEC-HMS topology.
            r_name    = f"W{pid}-R{this_sid}"
            _dsd = (sub_data.get(this_ds_sid) or {}).get('dam')
            r_ds_name = ((f"W{pid}-J{this_ds_sid}"
                          + (f"({_dsd['name']})" if _dsd else ""))
                         if this_ds_sid > 0 else catch_name)
            reach_data[this_sid] = dict(
                sid       = this_sid,
                ds_sid    = this_ds_sid,
                geom      = rg,
                rch       = rch,
                tc_min    = tc_min,
                tc_label  = tc_label,
                name      = r_name,
                ds_name   = r_ds_name,
                # outflow hydrographs filled in pass 3
                routed    = {},   # rp -> routed Q array
                inflows   = {},   # rp -> inflow Q array (combined at junction)
            )

        # ----- Pass 3: topological combine + reach routing -----------
        # combined_at_junction[sid][rp] = Q array at the junction at
        # the OUTLET of sub-basin `sid` for return period `rp`.
        # That is: sum of (local runoff of sid) + (routed inflows from
        # upstream sub-basins draining directly to sid).
        combined_at_junction = {s: {} for s in sub_data}
        do_routing = (do_hydrographs and use_cn and active_periods)
        if do_routing:
            order = topological_order_subbasins(sid_of, ds_sid_of)
            # Length is fixed catchment-wide by Pass 0; fall back to
            # whatever the first available local hydrograph reports.
            n_total = catch_n_total
            if n_total is None:
                for sd in sub_data.values():
                    for rp_, h in sd['local_hyds'].items():
                        n_total = len(h['Q_m3s'])
                        break
                    if n_total is not None:
                        break

            for sid in order:
                sd = sub_data.get(sid)
                if sd is None:
                    continue
                for rp, p_eff in active_periods_eff:
                    h = sd['local_hyds'].get(rp)
                    if h is None or n_total is None:
                        continue
                    Q_local = _align_len(h['Q_m3s'], n_total)
                    # accumulate routed inflows from upstream sub-basins
                    Q_in_sum = np.zeros(n_total, dtype=float)
                    for sd_up in sub_data.values():
                        if sd_up['ds_sid'] != sid:
                            continue
                        rd_up = reach_data.get(sd_up['sid'])
                        if rd_up is None:
                            continue
                        Q_routed = rd_up['routed'].get(rp)
                        if Q_routed is None:
                            continue
                        Q_in_sum += _align_len(Q_routed, n_total)
                    Q_comb = Q_local + Q_in_sum
                    combined_at_junction[sid][rp] = Q_comb
                    # ---- v24: reservoir storage routing ----------
                    # If this junction is a special point with a
                    # usable reservoir definition, the combined
                    # inflow is routed through the modified-Puls
                    # engine; the OUTFLOW (not the inflow) is what
                    # continues down the reach.
                    Q_down = Q_comb
                    _dam_p3 = sd.get('dam')
                    if _dam_p3 and _dam_p3.get('sd_curve'):
                        _dt_s3 = float(cfg.get('HYD_TIME_STEP_MIN',
                                               1.0)) * 60.0
                        _rr = route_reservoir(
                            Q_comb, _dam_p3['sd_curve'], _dt_s3,
                            init_storage=_dam_p3.get('init_storage'))
                        sd.setdefault('res_results', {})[rp] = _rr
                        if _rr['changed']:
                            Q_down = _rr['O']
                            _qpi = float(np.max(Q_comb))
                            _qpo = float(np.max(Q_down))
                            _att = (100.0 * (1.0 - _qpo / _qpi)
                                    if _qpi > 0 else 0.0)
                            print(f"       Sub{sid} "
                                  f"[{_dam_p3['name']}] {rp}-yr: "
                                  f"reservoir Qin={_qpi:.2f} → "
                                  f"Qout={_qpo:.2f} m³/s "
                                  f"(attenuation {_att:.0f}%, "
                                  f"Smax="
                                  f"{_rr['max_storage']/1000.0:.1f}"
                                  f"×10³ m³"
                                  + (f", WLmax={_rr['max_stage']:.2f} m"
                                     if _rr['max_stage'] is not None
                                     else "") + ")")
                        if _rr.get('warn'):
                            print(f"       ⚠ Sub{sid} "
                                  f"[{_dam_p3['name']}] {rp}-yr: "
                                  f"{_rr['warn']}")
                    # if this sub-basin has a downstream reach, route
                    # Q_down through it to populate rd['routed'][rp]
                    rd = reach_data.get(sid)
                    if rd is not None:
                        rd['inflows'][rp] = Q_down
                        Q_out = route_reach(
                            Q_down,
                            rd['tc_min'],
                            method=reach_routing_cfg,
                            muskingum_x=musk_x_cfg,
                            muskingum_n_sub=(musk_nsub_cfg or None),
                            dt_min=float(cfg.get('HYD_TIME_STEP_MIN', 1.0)),
                            # v19: geometry for M-Cunge —
                            # L and slope from the delineated reach,
                            # trapezoid shape from CONFIG.
                            reach_length_m=float(rd['rch'].get('L') or 0.0),
                            reach_slope=max(
                                float(rd['rch'].get('SLOPE') or 0.0), 1e-5),
                            mc_manning_n=mc_n_cfg,
                            mc_bottom_width_m=mc_b_cfg,
                            mc_side_slope=mc_z_cfg,
                        )
                        rd['routed'][rp] = _align_len(Q_out, n_total)

        # ----- Console summary of routing for this catchment ---------
        if do_routing:
            print(f"     v17: {len(sub_data)} sub-basin(s) processed; "
                  f"{len(reach_data)} reach(es) routed via "
                  f"{reach_routing_cfg}"
                  + (f" (X={musk_x_cfg})"
                     if reach_routing_cfg == 'MUSKINGUM' else ""))


        # ============================================================
        # v19: collect this catchment for HEC-HMS export
        # ============================================================
        if net_collect_cfg:
            try:
                _tlag_c = locals().get('tlag_sel') or 0.0
                if not _tlag_c and ch:
                    _tlag_c = float(ch.get('TLAG_KIRPH') or
                                    ch.get('TLAG') or 0.0)
                _subs_x, _reaches_x = [], []
                for _sid, _sd in (sub_data or {}).items():
                    _stl = float(_sd.get('tlag_used') or 0.0)
                    if not _stl:
                        _sch = _sd.get('ch') or {}
                        _stl = float(_sch.get('TLAG_KIRPH') or
                                     _sch.get('TLAG') or 0.0)
                    _sub_cxy = None
                    if _sd.get('mask') is not None:
                        _sub_cxy = mask_centroid_xy(_sd['mask'], gt)
                    _subs_x.append(dict(
                        sid=_sid, ds_sid=_sd['ds_sid'],
                        name=_sd['name'],
                        centroid_xy=_sub_cxy,
                        depths=dict(_sd.get('periods_eff')
                                    or active_periods_eff or []),
                        # v1.0 — ARF actually applied to this
                        # sub-basin. 'depths' stay EFFECTIVE (already
                        # reduced) so the HEC-HMS export is untouched;
                        # the AR-Awdiya export divides them back out
                        # and exports the ARF separately.
                        arf=float(arf_value or 1.0),
                        area_km2=float(_sd.get('area_km2') or 0.0),
                        cn=float(_sd.get('cn') or cn_avg or 0.0),
                        lag_min=_stl,
                        junction_xy=_sd['junction_xy'],
                        # v23 — reservoir/dam name (None when the
                        # junction is not a dam location); appended to
                        # the .basin junction element name.
                        dam_name=((_sd.get('dam') or {}).get('name')
                                  or None),
                        # v24 — reservoir definition for the .basin
                        # Reservoir element + .pdata Elevation-Storage
                        # table (None → plain junction).
                        reservoir=(lambda _dd: (dict(
                            mode=(_dd['res'] or {}).get('mode', 'table'),
                            ev_rows=_dd['sd_curve'].get('ev_rows'),
                            sd_rows=list(zip(_dd['sd_curve']['S'],
                                             _dd['sd_curve']['O'])),
                            structures=((_dd['res'] or {})
                                        .get('structures') or []),
                            init_mode=str((_dd['res'] or {}).get(
                                'init_mode', 'inflow')),
                            init_elev=(_dd['res'] or {}).get('init_elev'),
                        ) if (_dd and _dd.get('sd_curve')) else None))(
                            _sd.get('dam')),
                    ))
                for _sid, _rd in (reach_data or {}).items():
                    _rch_x = _rd.get('rch') or {}
                    _routing = dict(method=reach_routing_cfg)
                    if reach_routing_cfg == 'LAG':
                        _routing['lag_min'] = float(_rd.get('tc_min') or 0.0)
                    elif reach_routing_cfg == 'MUSKINGUM':
                        # plugin routes with K = reach Tc; export the same
                        _routing['K_hr'] = float(_rd.get('tc_min')
                                                 or 0.0) / 60.0
                        _routing['X'] = musk_x_cfg
                        # v1.0 — same subreach count the in-plugin
                        # cascade used, so HEC-HMS reproduces the
                        # plugin's routed hydrograph
                        _routing['subreaches'] = (
                            musk_nsub_cfg or muskingum_subreaches(
                                float(_rd.get('tc_min') or 0.0),
                                musk_x_cfg,
                                float(cfg.get('HYD_TIME_STEP_MIN', 1.0))))
                    elif reach_routing_cfg == 'MUSKINGUM_CUNGE':
                        _routing.update(
                            L_m=float(_rch_x.get('L') or 0.0),
                            S0=max(float(_rch_x.get('SLOPE') or 0.0), 1e-5),
                            b_m=mc_b_cfg, z=mc_z_cfg, n=mc_n_cfg)
                    _mid = None
                    _g = _rd.get('geom')
                    if _g is not None and not _g.isEmpty():
                        _c = _g.interpolate(_g.length() * 0.5)
                        if _c and not _c.isNull():
                            _p = _c.asPoint()
                            _mid = (float(_p.x()), float(_p.y()))
                    _reaches_x.append(dict(
                        sid=_sid, ds_sid=_rd['ds_sid'], name=_rd['name'],
                        routing=_routing, mid_xy=_mid))
                hms_catchments.append(dict(
                    pid=pid, name=catch_name,
                    centroid_xy=mask_centroid_xy(catch, gt),
                    depths=dict(active_periods_eff or []),
                    arf=float(arf_value or 1.0),      # v1.0
                    outlet_xy=(float(sx), float(sy)),
                    area_km2=float(area_km2 or 0.0),
                    cn=float(cn_avg or 0.0) if use_cn else 0.0,
                    lag_min=float(_tlag_c or 0.0),
                    subs=_subs_x, reaches=_reaches_x,
                ))
            except Exception as _e_hms:
                print(f"     ⚠ HEC-HMS collection failed for "
                      f"{catch_name}: {_e_hms}")

        # ============================================================
        # Build & emit features (Pass 3b)
        # ============================================================

        # ---- emit one feature per sub-basin (shared schema with SubLFPs) -
        for sid, sd in sub_data.items():
            this_ds   = sd['ds_sid']
            sub_ch    = sd['ch']
            sub_area  = sd['area_km2']
            sub_bslp  = sd['bslp']
            sub_bslp_p = sd['bslp_p']
            sub_cn    = sd['cn']
            sub_cn_cov= sd['cn_cov']
            sub_tlag_nrcs = sd['tlag_nrcs']
            sub_tc_nrcs   = sd['tc_nrcs']
            sub_combined_attrs = [
                int(pid), int(sid), int(this_ds), int(sd['is_outlet']),
                None, None,                 # NAME / DS_NAME placeholders
                float(sub_area), float(sub_bslp), float(sub_bslp_p),
                sub_ch['L'],
                sub_ch['ELE_US'], sub_ch['ELE_DS'],
                sub_ch['ELE_10'], sub_ch['ELE_85'],
                sub_ch['SLOPE'], sub_ch['SLOPE_PCT'],
                sub_ch['SLP_US_85'], sub_ch['SLP_85_10'], sub_ch['SLP_10_DS'],
                sub_ch['TC_US_85'], sub_ch['TC_85_10'], sub_ch['TC_10_DS'],
                sub_ch['TC_TOTAL'], sub_ch['TLAG'],
                # v14 — holistic Kirpich
                float(sub_ch['TC_KIRPH']), float(sub_ch['TLAG_KIRPH']),
                int(round(prf_cfg)),        # v1.0 — adopted PRF
            ]
            if use_cn:
                sub_combined_attrs += [
                    float(sub_cn) if sub_cn is not None else None,
                    float(sub_cn_cov) if sub_cn_cov is not None else None,
                    float(sub_tlag_nrcs) if sub_tlag_nrcs is not None else None,
                    float(sub_tc_nrcs)   if sub_tc_nrcs   is not None else None,
                ]
            if arf_active:
                sub_combined_attrs.append(float(arf_value))
            # v17 — design rainfall block on sub-basin
            if active_periods:
                sub_combined_attrs += design_rainfall_attr_block(
                    sd.get('periods_eff') or active_periods_eff,
                    sub_cn if use_cn else None,
                    sub_area,
                    use_cn,
                    qpeaks=(sd['local_qpeaks'] if do_hydrographs else None),
                    scs_lambda=scs_lambda_cfg)

            # SubBasins feature (polygon, sub-basin naming)
            sub_geom = mask_to_qgs_geometry(sd['mask'], gt, proj)
            if sub_geom is not None:
                attrs_for_sub = list(sub_combined_attrs)
                attrs_for_sub[4] = sd['name']
                attrs_for_sub[5] = sd['ds_name']
                sf = QgsFeature()
                sf.setGeometry(sub_geom)
                sf.setAttributes(attrs_for_sub)
                sub_feats.append(sf)

            # SubLFPs feature (line, sub-LFP naming)
            attrs_for_sublfp = list(sub_combined_attrs)
            attrs_for_sublfp[4] = sd['lfp_name']
            attrs_for_sublfp[5] = sd['ds_lfp_name']
            sl = QgsFeature()
            sl.setGeometry(sd['line'])
            sl.setAttributes(attrs_for_sublfp)
            sublfp_feats.append(sl)

        # ---- emit one reach per non-outlet sub-basin -----------------
        for sid, rd in reach_data.items():
            rch = rd['rch']
            r_attrs = [
                int(pid), int(sid), int(rd['ds_sid']),
                rd['name'], rd['ds_name'],
                rch['L'],
                rch['ELE_US'], rch['ELE_DS'], rch['ELE_10'], rch['ELE_85'],
                rch['SLOPE'], rch['SLOPE_PCT'],
                rch['SLP_US_85'], rch['SLP_85_10'], rch['SLP_10_DS'],
                rch['TC_US_85'], rch['TC_85_10'], rch['TC_10_DS'],
                rch['TC_TOTAL'],
                # v14 — holistic Kirpich on the reach
                float(rch['TC_KIRPH']),
            ]
            # v17 — routing-method tag and per-rp routed peaks + attenuation
            r_attrs += [reach_routing_cfg, float(rd['tc_min'])]
            # v1.0 — Muskingum K (hours) + subreach count used
            if reach_routing_cfg == 'MUSKINGUM':
                _k_hr = float(rd['tc_min']) / 60.0
                _n_sr = (musk_nsub_cfg or muskingum_subreaches(
                    float(rd['tc_min']), musk_x_cfg,
                    float(cfg.get('HYD_TIME_STEP_MIN', 1.0))))
                r_attrs += [_k_hr, int(_n_sr)]
            else:
                r_attrs += [None, None]
            if active_periods and do_hydrographs:
                for rp, _ in active_periods:
                    Q_in  = rd['inflows'].get(rp)
                    Q_out = rd['routed'].get(rp)
                    if Q_in is None or Q_out is None:
                        r_attrs += [None, None, None]
                    else:
                        qp_in  = float(np.max(Q_in))
                        qp_out = float(np.max(Q_out))
                        atten  = (100.0 * (1.0 - qp_out / qp_in)
                                  if qp_in > 0 else 0.0)
                        r_attrs += [qp_in, qp_out, atten]
                        # write routed reach hydrograph CSV + HTML (v18)
                        csv_name = f"Reach_{rd['name']}_{rp}yr.csv"
                        n_q = len(Q_out)
                        t_q = np.arange(n_q, dtype=float) * float(
                            cfg.get('HYD_TIME_STEP_MIN', 1.0))
                        write_q_only_csv(
                            os.path.join(modeling_dirs['REACH']['csv'],
                                         csv_name),
                            t_q, Q_out, qin_m3s=Q_in)
                        html_path = os.path.join(
                            modeling_dirs['REACH']['html'],
                            csv_name[:-4] + '.html')
                        write_hydrograph_html(
                            html_path, kind='REACH',
                            title=f"Reach {rd['name']} → {rd['ds_name']} "
                                  f"({rp}-yr storm)",
                            subtitle=f"Routing: {reach_routing_cfg}"
                                     + (f" (X={musk_x_cfg})"
                                        if reach_routing_cfg == 'MUSKINGUM'
                                        else "")
                                     + f"  •  Tc={rd['tc_min']:.1f} min "
                                     f"({rd['tc_label']})  •  L={rch['L']:.0f} m",
                            hyd_data={
                                'time_min':  list(t_q),
                                'Q_in_m3s':  list(Q_in),
                                'Q_out_m3s': list(Q_out),
                            },
                        )
            rf = QgsFeature()
            rf.setGeometry(rd['geom'])
            rf.setAttributes(r_attrs)
            reach_feats.append(rf)

        # ---- emit one junction feature per sub-basin outlet ----------
        # v20.1: true cumulative upstream area per junction. Follow
        # the ds_sid topology transitively (memoised DFS) so chained
        # sub-basins accumulate correctly, not just direct inflowers.
        _ups_children = {s: [] for s in sub_data}
        for _s, _d in sub_data.items():
            _ds = _d['ds_sid']
            if _ds in _ups_children:
                _ups_children[_ds].append(_s)
        _cum_area_memo = {}
        def _cum_area(s):
            if s in _cum_area_memo:
                return _cum_area_memo[s]
            a = float(sub_data[s]['area_km2'] or 0.0)
            for u in _ups_children[s]:
                a += _cum_area(u)
            _cum_area_memo[s] = a
            return a
        # v23 — junction display name, with the dam-name suffix when the
        # junction was enforced by a reservoir. v1.0 naming:
        # "W{pid}-J{sid}", e.g. "W1-J9(DamName)".
        def _junc_display_name(s):
            base = f"W{pid}-J{s}"
            _d = (sub_data.get(s) or {}).get('dam')
            return f"{base}({_d['name']})" if _d else base

        for sid, sd in sub_data.items():
            jx, jy = sd['junction_xy']
            j_name = _junc_display_name(sid)
            if sid == 1:
                # outlet junction → drains to the catchment
                j_ds_name = catch_name
            else:
                j_ds_name = _junc_display_name(sd['ds_sid'])
            # v20.1: single area attribute — total catchment area
            # upstream of this junction (transitive, memoised).
            a_ups = _cum_area(sid)
            j_attrs = [
                int(pid), int(sid), j_name, j_ds_name,
                int(sd['is_outlet']),
                float(a_ups),
            ]
            if arf_active:
                j_attrs.append(float(arf_value))
            if active_periods and do_hydrographs:
                for rp, _ in active_periods:
                    Q_comb = combined_at_junction.get(sid, {}).get(rp)
                    if Q_comb is None:
                        j_attrs.append(None)
                        j_attrs.append(None)   # V_<rp>yrs_1000m3
                    else:
                        # v24 — reservoir junction: QP/V report the
                        # routed OUTFLOW (what continues downstream);
                        # Q_comb remains the inflow.
                        _rr = (sd.get('res_results') or {}).get(rp)
                        _is_res = bool(_rr and _rr.get('changed'))
                        Q_show = _rr['O'] if _is_res else Q_comb
                        qp_comb = float(np.max(Q_show))
                        # v1.0 — completeness check: the output window
                        # should extend until the flow recedes to
                        # (practically) zero.
                        if qp_comb > 0 and \
                                float(Q_show[-1]) > 0.005 * qp_comb:
                            print(f"     ⚠ {j_name} {rp}-yr: hydrograph "
                                  f"NOT fully receded at the end of the "
                                  f"output window (Q_end = "
                                  f"{float(Q_show[-1]):.3f} m³/s = "
                                  f"{100.0*float(Q_show[-1])/qp_comb:.1f}% "
                                  f"of peak). Increase 'Hydrograph tail "
                                  f"length' / 'Max tail' in the dialog "
                                  f"to capture the full recession.")
                        j_attrs.append(qp_comb)
                        # v20: junction hydrograph volume, trapezoidal
                        # integration of Q over the run window,
                        # expressed in 10³ m³ for consistency with the
                        # V_<n>yrs_1000m3 catchment/sub-basin columns.
                        _dt_s = float(cfg.get('HYD_TIME_STEP_MIN',
                                              1.0)) * 60.0
                        _v_m3 = float(np.trapz(Q_show, dx=_dt_s)) \
                            if hasattr(np, 'trapz') \
                            else float(np.trapezoid(Q_show, dx=_dt_s))
                        j_attrs.append(_v_m3 / 1000.0)
                        # write junction hydrograph CSV + HTML (v18);
                        # reservoir junctions get the AR-Awdiya-style
                        # routing table + inflow/outflow chart (v24).
                        csv_name = f"{j_name}_{rp}yr.csv"
                        n_q = len(Q_show)
                        t_q = np.arange(n_q, dtype=float) * float(
                            cfg.get('HYD_TIME_STEP_MIN', 1.0))
                        html_path = os.path.join(
                            modeling_dirs['JUNCTION']['html'],
                            csv_name[:-4] + '.html')
                        if _is_res:
                            write_reservoir_csv(
                                os.path.join(
                                    modeling_dirs['JUNCTION']['csv'],
                                    csv_name),
                                t_q, Q_comb, Q_show, _rr['storage'],
                                _rr['stage'])
                            _qpi = float(np.max(Q_comb))
                            _att = (100.0 * (1.0 - qp_comb / _qpi)
                                    if _qpi > 0 else 0.0)
                            _extra = [
                                ('Peak storage',
                                 f"{_rr['max_storage']/1000.0:.1f}",
                                 '×10³ m³')]
                            if _rr['max_stage'] is not None:
                                _extra.append(
                                    ('Peak stage',
                                     f"{_rr['max_stage']:.2f}", 'm'))
                            write_hydrograph_html(
                                html_path, kind='REACH',
                                title=f"Reservoir {j_name} → "
                                      f"{j_ds_name} ({rp}-yr storm)",
                                subtitle=f"Storage routing (modified "
                                         f"Puls)  •  "
                                         f"{(_rr.get('info') or '')}"
                                         f"  •  Upstream area "
                                         f"A_ups={a_ups:.2f} km²",
                                hyd_data={
                                    'time_min':  list(t_q),
                                    'Q_in_m3s':  list(Q_comb),
                                    'Q_out_m3s': list(Q_show),
                                },
                                extra_summary=_extra,
                            )
                        else:
                            write_q_only_csv(
                                os.path.join(
                                    modeling_dirs['JUNCTION']['csv'],
                                    csv_name),
                                t_q, Q_show)
                            write_hydrograph_html(
                                html_path, kind='JUNCTION',
                                title=f"Junction {j_name} → "
                                      f"{j_ds_name} ({rp}-yr storm)",
                                subtitle=f"Combined flow at sub-basin "
                                         f"outlet  •  Upstream area "
                                         f"A_ups={a_ups:.2f} km²",
                                hyd_data={
                                    'time_min': list(t_q),
                                    'Q_m3s':    list(Q_show),
                                },
                            )
            # v23 — dam characteristics on the target junction. Order
            # matches the schema tail: DAM_NAME, then one column per
            # selected dam field (NULL on non-dam junctions).
            if dam_cells:
                _d = sd.get('dam')
                j_attrs.append(_d['name'] if _d else None)
                for _src, _out, _fld in dam_attr_cols:
                    j_attrs.append(_d['attrs'].get(_out)
                                   if _d else None)
            # v24 — reservoir routing outputs (order mirrors the
            # schema tail: RES_INFO then per-rp QIN/ATT/SMAX/WLMX).
            if res_any:
                _d = sd.get('dam')
                _sdc = (_d or {}).get('sd_curve')
                j_attrs.append(_sdc['info'] if _sdc else None)
                if active_periods:
                    for rp, _ in active_periods:
                        _rr = (sd.get('res_results') or {}).get(rp)
                        if _rr and _rr.get('changed'):
                            _qi = float(np.max(
                                combined_at_junction[sid][rp]))
                            _qo = float(np.max(_rr['O']))
                            j_attrs.append(_qi)
                            j_attrs.append(
                                100.0 * (1.0 - _qo / _qi)
                                if _qi > 0 else 0.0)
                            j_attrs.append(
                                _rr['max_storage'] / 1000.0
                                if _rr['max_storage'] is not None
                                else None)
                            j_attrs.append(_rr['max_stage'])
                        else:
                            j_attrs += [None, None, None, None]
            jf = QgsFeature()
            jf.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(jx, jy)))
            jf.setAttributes(j_attrs)
            junc_feats.append(jf)

        n_subs    = len(sub_data)
        n_reaches = len(reach_data)
        n_juncs   = len(sub_data)
        print(f"     → {n_subs} sub-basin(s), {n_reaches} reach(es), "
              f"{n_juncs} junction(s) emitted")

    catch_out.dataProvider().addFeatures(cf)
    lfp_out.dataProvider().addFeatures(lf)
    pts_out.dataProvider().addFeatures(pf)
    snap_out.dataProvider().addFeatures(snapf)
    sub_out.dataProvider().addFeatures(sub_feats)
    sublfp_out.dataProvider().addFeatures(sublfp_feats)
    reach_out.dataProvider().addFeatures(reach_feats)
    junc_out.dataProvider().addFeatures(junc_feats)

    # v16 — pair each layer with its CONFIG key for QML lookup. Order of
    # pairs is preserved for the addMapLayer step, matching v15 behaviour.
    _prog(0.90, "Adding layers to the project")
    base_layers = [
        (catch_out, 'OUT_CATCHMENTS'),
        (lfp_out,   'OUT_LFP_LINES'),
        (pts_out,   'OUT_LFP_POINTS'),
        (snap_out,  'OUT_SNAPPED'),
    ]
    if thiessen_out is not None:                    # v1.0
        base_layers.append((thiessen_out, 'OUT_THIESSEN'))
    sub_layers = [
        (sub_out,    'OUT_SUBBASINS'),
        (sublfp_out, 'OUT_SUBLFPS'),
        (reach_out,  'OUT_REACHES'),
        (junc_out,   'OUT_JUNCTIONS'),    # v17 — junctions point layer
    ] if do_subdivide else []
    for lyr, cfg_key in base_layers + sub_layers:
        lyr.updateExtents()
        project.addMapLayer(lyr)
        # apply the registered QML style (silent no-op when STYLES_DIR
        # is empty or no QML is registered for this layer)
        style_vector_layer(lyr, cfg_key, styles_dir, verbose=True)

    # ---- emit the catchment-clipped CN raster (only if CN is enabled) ---
    cn_clipped_path = None
    if use_cn and cn_grid is not None and catch_union_mask is not None:
        n_in = int(catch_union_mask.sum())
        if n_in > 0:
            _prog(0.88, "Writing output rasters")
            print(f"\n[7] Writing CN raster clipped to catchment union "
                  f"({n_in:,} cells) ...")
            # Choose a clip nodata. Use the source nodata if it exists,
            # otherwise -9999 (kept distinct from typical CN range 30-100).
            clip_nodata = cn_nodata if cn_nodata is not None else -9999
            cn_clip = np.where(catch_union_mask, cn_grid,
                               clip_nodata).astype(np.float32)
            # Where the source CN itself was nodata, preserve it.
            if cn_nodata is not None:
                cn_clip[cn_grid == cn_nodata] = clip_nodata
            cn_clipped_path = os.path.join(out_dir, 'CN_clipped.tif')
            write_raster(cn_clip, gt, proj, cn_clipped_path,
                         dtype=gdal.GDT_Float32, nodata=clip_nodata)
            cn_rl = QgsRasterLayer(cn_clipped_path, 'CN_clipped')
            if cn_rl.isValid():
                project.addMapLayer(cn_rl)
                print(f"    -> {cn_clipped_path} added to project as "
                      f"'CN_clipped'.")
            else:
                print(f"    ⚠ Clipped CN raster wrote OK but failed to "
                      f"register as a QGIS layer.")

    # ---- v19: write the HEC-HMS .basin file -----------------
    hms_basin_path = None
    if hms_export_cfg:
        if _write_hms_basin is None:
            print("⚠ EXPORT_HEC_HMS requested but the "
                  "hec_hms_export module is unavailable — skipped.")
        elif not hms_catchments:
            print("⚠ EXPORT_HEC_HMS: no catchments collected — "
                  "nothing to write.")
        else:
            try:
                hms_basin_path = _write_hms_basin(
                    os.path.join(out_dir, 'HEC_HMS',
                                 f'{hms_name_cfg}.basin'),
                    hms_name_cfg,
                    hms_catchments,
                    version=hms_version_cfg,
                    prf=prf_cfg,          # v1.0 — SCS unitgraph type
                )
                try:
                    from .hec_hms_export import hms_unitgraph_type
                except ImportError:
                    from hec_hms_export import hms_unitgraph_type
                _ug_kw, _ug_prf = hms_unitgraph_type(prf_cfg)
                print(f"    Unitgraph Type: {_ug_kw} "
                      f"(peak rate factor {_ug_prf:g})")
                if abs(_ug_prf - prf_cfg) > 0.5:
                    print(f"    \u26a0 PEAK_RATE_FACTOR={prf_cfg:g} is not "
                          f"one of the peak rate factors HEC-HMS offers "
                          f"(284, 484 and 100\u2013600 in steps of 50), "
                          f"so the export uses the nearest, "
                          f"{_ug_prf:g} ({_ug_kw}). AR-Tarseem's own "
                          f"hydrographs still use {prf_cfg:g} \u2014 pick "
                          f"a supported value if the two must agree.")
                n_sub_x = sum(max(1, len(c['subs']))
                              for c in hms_catchments)
                n_rch_x = sum(len(c['reaches']) for c in hms_catchments)
                _prog(0.93, "HEC-HMS export")
                print(f"\n[8] HEC-HMS export (v{hms_version_cfg}): "
                      f"{len(hms_catchments)} basin(s), "
                      f"{n_sub_x} subbasin(s), {n_rch_x} reach(es)")
                print(f"    -> {hms_basin_path}")
                print("    Open in HEC-HMS via File → Import → "
                      "Basin Model, or in WMS.")
                # ---- v22: companion .met files (one per return period)
                if bool(cfg.get('EXPORT_HEC_HMS_MET', False)):
                    try:
                        from .hec_hms_export import write_met_files
                    except ImportError:
                        from hec_hms_export import write_met_files
                    # storm duration: a CUSTOM hyetograph fixes it from
                    # its own time axis, otherwise STORM_DURATION_HR
                    _met_dur_hr = float(
                        cfg.get('STORM_DURATION_HR', 24.0))
                    if (str(cfg.get('STORM_DISTRIBUTION', '')).upper()
                            == 'CUSTOM'
                            and custom_hyeto_cfg is not None):
                        _met_dur_hr = float(
                            custom_hyeto_cfg['duration_hr'])
                    _met_paths = write_met_files(
                        os.path.join(out_dir, 'HEC_HMS'),
                        hms_name_cfg,
                        hms_catchments,
                        storm_distribution=str(
                            cfg.get('STORM_DISTRIBUTION', 'SCS_TYPE_II')),
                        version=hms_version_cfg,
                        storm_duration_hr=_met_dur_hr)
                    if _met_paths:
                        print(f"    + {len(_met_paths)} .met file(s) "
                              f"(one per return period, per-subbasin "
                              f"storm depths)")

                print(f"    For one-to-one comparison with HEC-HMS, "
                      f"set the HMS control interval equal to "
                      f"HYD_TIME_STEP_MIN "
                      f"({cfg.get('HYD_TIME_STEP_MIN', 1.0):g} min): "
                      f"peaks of fast basins are sensitive to the "
                      f"computation interval.")
            except Exception as _e_hms:
                print(f"⚠ HEC-HMS export failed: {_e_hms}")

    # ---- v1.0: AR-Awdiya project (.json) ----------------------------
    # Independent of the HEC-HMS export: it reuses the same collected
    # network payload and writes a project AR-Awdiya can open and run.
    if bool(cfg.get('EXPORT_AWDIYA', False)) and hms_catchments:
        try:
            try:
                from .awdiya_export import write_awdiya_project
            except ImportError:
                from awdiya_export import write_awdiya_project
            _aw_dur_hr = float(cfg.get('STORM_DURATION_HR', 24.0))
            if (str(cfg.get('STORM_DISTRIBUTION', '')).upper() == 'CUSTOM'
                    and custom_hyeto_cfg is not None):
                _aw_dur_hr = float(custom_hyeto_cfg['duration_hr'])
            _aw_rp = cfg.get('AWDIYA_PRIMARY_RP')
            _aw_path = write_awdiya_project(
                os.path.join(out_dir, 'AR_Awdiya',
                             'AR_Tarseem_project.json'),
                hms_catchments,
                project_name=(str(cfg.get('AWDIYA_PROJECT_NAME')
                                  or '').strip() or 'AR-Tarseem model'),
                dt_min=float(cfg.get('HYD_TIME_STEP_MIN', 1.0)),
                storm_distribution=str(
                    cfg.get('STORM_DISTRIBUTION', 'SCS_TYPE_II')),
                storm_duration_hr=_aw_dur_hr,
                prf=prf_cfg,        # v1.0 fix: was cfg['SCS_UH_PRF'] (never set)
                scs_lambda=float(cfg.get('SCS_LAMBDA', 0.2)),
                freq_mode=str(cfg.get('AWDIYA_FREQ_MODE', 'all')),
                primary_rp=(int(_aw_rp) if _aw_rp else None),
                hms_version=str(cfg.get('HEC_HMS_VERSION', '4.11')))
            _n_c = sum(len(_c.get('subs') or []) or 1
                       for _c in hms_catchments)
            print(f"\n[8b] AR-Awdiya project: {_n_c} catchment(s), "
                  f"{len(hms_catchments)} outlet(s)")
            print(f"    -> {_aw_path}")
            print('    Open it in AR-Awdiya with "Load project" — '
                  'the model is ready to Run.')
        except Exception as _e_aw:
            print(f"    \u26a0 AR-Awdiya export failed: {_e_aw}")

    print(f"\n{'='*60}")
    print(f"✅ DONE  catchments={len(cf)}  LFPs={len(lf)}  "
          f"LFP points={len(pf)}  streams={streams_out.featureCount()}  "
          f"snapped={len(snapf)}")
    if do_subdivide:
        print(f"   sub-basins={len(sub_feats)}  "
              f"sub-LFPs={len(sublfp_feats)}  "
              f"reaches={len(reach_feats)}  "
              f"junctions={len(junc_feats)}")
    if cn_clipped_path:
        print(f"   CN_clipped: {cn_clipped_path}")
    if do_hydrographs and hyd_out_dir:
        # v18 — count CSVs and HTMLs across the structured tree
        n_csv = n_html = 0
        for root_, _dirs, files_ in os.walk(hyd_out_dir):
            for f in files_:
                if f.endswith('.csv'):
                    n_csv += 1
                elif f.endswith('.html'):
                    n_html += 1
        print(f"   Modeling_Outputs/: {n_csv} CSV(s) + {n_html} "
              f"HTML(s) in {hyd_out_dir}")
        print(f"     TC_METHOD={tc_method_cfg}  PRF={prf_cfg:g}")
        if do_subdivide:
            print(f"     REACH_ROUTING={reach_routing_cfg}"
                  + (f" (X={musk_x_cfg})"
                     if reach_routing_cfg == 'MUSKINGUM' else "")
                  + f"  REACH_TC_METHOD={reach_tc_cfg}")
    if arf_active:
        print(f"   ARF_METHOD={arf_method_cfg} applied to design rainfall")
    if styles_dir:
        print(f"   Styling: QML files from {styles_dir}")
    if hms_basin_path:
        print(f"   HEC-HMS: {hms_basin_path}")
    print(f"{'='*60}")
    _prog(1.0, "Done")


# =============================================================================
# When executed directly we still call run(); when imported as a module
# (e.g. by the plugin) nothing fires.
# =============================================================================
if __name__ == '__main__':
    run()
