# AR-Tarseem — QGIS Plugin (v1.0.0)

**Copyright (C) 2026 Ahmed M. Raslan, PhD & Amr R. Refaiy, PhD —
Water Engineering Experts**

Licensed under the **GNU General Public License v3.0 or later**
(GPL-3.0-or-later). See [LICENSE](LICENSE) for the full text.

AR-Tarseem is a DEM-based watershed, sub-basin, longest-flow-path,
SCS-CN/UH design-hydrograph, reservoir-routing, reach-routing and
HEC-HMS interchange tool for QGIS. Pure NumPy + GDAL — no external
Python packages beyond what ships with QGIS.

## Capabilities

1. **Delineation.** Streams (flow-accumulation threshold), one
   catchment per pour point, longest flow paths with segment slopes
   and three time-of-concentration estimates (Kirpich holistic,
   Kirpich 10-85, NRCS lag), optional HEC-HMS-style subdivision into
   sub-basins with sub-LFPs, reaches and junctions.

2. **Special delineation points & reservoirs.** Enforce junctions at
   dam or structure locations (snapped to the main flow path), divide
   sub-catchments at the points only or combined with the area-based
   subdivision, and define any point as a reservoir with modified-Puls
   storage routing — elevation–volume curve plus spillway/orifice
   structures, or a direct storage–discharge table. Outputs include
   inflow, outflow, attenuation, peak storage and peak stage.

3. **Design storms.** All four SCS/NRCS 24-h distributions (Type I,
   IA, II, III) at 1-minute tabular resolution, or a **custom
   hyetograph** CSV: two columns — time (minutes), depth (mm); header
   optional; comma/semicolon/tab/space separated; cumulative vs
   incremental depths auto-detected; the curve is normalised and
   rescaled to each return period's total P, and its last time defines
   the storm duration.

4. **Losses.** SCS-CN with Curve Number from a raster (auto-aligned to
   the DEM) or a single constant value, and adjustable
   initial-abstraction ratio λ = Ia/S (0.01–0.5; 0.2 classic, 0.05 per
   newer NRCS guidance), using the generalised runoff equation
   R = (P−λS)²/(P−λS+S).

5. **Areal Reduction Factor.** TRRL-UK or Riyadh equations computed
   per catchment from its full upstream area — or a **direct fixed
   value** (0 < ARF ≤ 1) applied to all catchments.

6. **Hydrographs.** SCS Unit Hydrograph (Gamma form, PRF 100–600) per
   catchment / sub-basin / junction / reach; CSVs plus interactive
   HTML visualisations in a `Modeling_Outputs/` tree, with an editable
   output window so every hydrograph recedes to zero.

7. **Rainfall distribution.** Custom return periods and
   spatially-variable rainfall through Thiessen station weighting.

8. **Reach routing.** Lag, Muskingum (K from reach Tc), or
   variable-parameter **Muskingum-Cunge** on a trapezoidal channel
   (Manning n, bottom width, side slope; length & slope from the
   delineated reach; K and X recomputed every time step; the reach
   subdivided into N = L/(cΔt) segments, ≤ 40).

9. **HEC-HMS export.** Writes `<OUTPUT_DIR>/HEC_HMS/AR_Tarseem.basin`:
   Subbasins (Area, SCS CN, SCS lag), Reaches, Junctions, true
   Reservoir elements and Sinks at real map coordinates, for HEC-HMS
   4.7–4.13 (version-dependent headers, verified against a genuine
   HEC-HMS 4.12 basin file). Plus per-return-period `.met` files and a
   `.pdata` carrying the reservoir Elevation-Storage tables. Open via
   *File → Import → Basin Model*, or in WMS.

10. **AR-Awdiya export.** Writes the whole model as a ready-to-run
    AR-Awdiya project (`.json`) — catchments, junctions, reach
    routing, reservoirs and design storms included.

## Install

**From the QGIS Plugin Manager:** *Plugins → Manage and Install
Plugins → All*, search for "AR-Tarseem", and click Install.

**From ZIP:** *Plugins → Manage and Install Plugins → Install from
ZIP* → select the plugin zip.

## Requirements

QGIS 3.22 or later. No external Python packages — the plugin uses only
`numpy`, `osgeo` (GDAL) and the QGIS API, all of which ship with QGIS.

## Files

- `hydro_core.py` — delineation + hydrology engine; `run(cfg)` entry
  point, `CONFIG` retained for QGIS-console use.
- `scs_data.py` — SCS distribution tables (4 × 1441 points) and the
  dimensionless UH reference curve.
- `hec_hms_export.py` — the `.basin` / `.met` / `.pdata` writer.
- `awdiya_export.py` — the AR-Awdiya project (`.json`) writer.
- `hydro_dialog.py` — tabbed dialog (Inputs / Delineation / Tc &
  losses / Rainfall & hydrographs / Routing & ARF / HEC-HMS /
  AR-Awdiya / Styling / Output names / About).
- `hydro_task.py` — background execution via `QgsTask`.
- `hydro_plugin.py`, `__init__.py`, `metadata.txt`, `icon.svg`,
  `icon.png`.

## Console use

`hydro_core.py` still runs standalone from the QGIS Python console:
edit the `CONFIG` block and execute the file (put `scs_data.py` and
`hec_hms_export.py` next to it to keep all features; without them it
falls back to the coarse Type II table and skips HEC-HMS export).

## Output location

All results are written into `<chosen folder>/ar_tarseem_output/`
(vector/raster outputs, `Modeling_Outputs/`, `HEC_HMS/`). Re-running
with the same folder reuses it without nesting.

## Matching HEC-HMS results

The hydrograph engine reproduces HEC-HMS numerics: SCS-CN losses and
the HEC-HMS dimensionless unit hydrograph (Tp = Δt/2 + Tlag, unit
volume normalised). To compare one-to-one, set the HEC-HMS control
interval equal to `HYD_TIME_STEP_MIN` — peak discharges of fast basins
are sensitive to the computation interval (a 1-min run reads several
percent higher than a 10-min run on short-lag basins).

## Progress & logging

The engine runs as a QGIS background task, so QGIS stays responsive
during a run, with a live progress bar reporting every stage and
catchment. Detailed messages appear in *View → Panels → Log Messages →
"AR-Tarseem"*.

## Issues

Please report bugs and feature requests on the
[issue tracker](https://github.com/amrrefaiy-hash/AR-Tarseem/issues).

## License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.
