# AR-Tarseem — QGIS Plugin (v1.1.0)

**Copyright (C) 2026 Ahmed M. Raslan, PhD & Amr R. Refaiy, PhD -
Water Engineering Experts.**

*Free software: you may redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free
Software Foundation, either version 2 of the License or (at your
option) any later version. It is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY, without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file
[`LICENSE`](LICENSE) for the full text.*

AR-Tarseem is a DEM-based watershed, sub-basin, longest-flow-path,
SCS-CN/UH design-hydrograph, reach-routing and HEC-HMS interchange
tool for QGIS. Pure NumPy + GDAL — no external Python packages.

## Capabilities

1. **Delineation.** Streams (flow-accumulation threshold), one
   catchment per pour point, longest flow paths with segment slopes
   and three time-of-concentration estimates (Kirpich holistic,
   Kirpich 10-85, NRCS lag), optional HEC-HMS-style subdivision into
   sub-basins with sub-LFPs, reaches and junctions.

2. **Design storms.** All four SCS/NRCS 24-h distributions (Type I,
   IA, II, III) at 1-minute tabular resolution, or a **custom
   hyetograph** CSV: two columns — time (minutes), depth (mm); header
   optional; comma/semicolon/tab/space separated; cumulative vs
   incremental depths auto-detected; the curve is normalised and
   rescaled to each return period's total P, and its last time defines
   the storm duration.

3. **Losses.** SCS-CN with adjustable initial-abstraction ratio
   λ = Ia/S (0.01–0.5; 0.2 classic, 0.05 per newer NRCS guidance),
   using the generalised runoff equation R = (P−λS)²/(P−λS+S).

4. **Areal Reduction Factor.** TRRL_UK or Riyadh equations computed
   per catchment from its full upstream area — or a **direct fixed
   value** (0 < ARF ≤ 1) entered by the user and applied to all
   catchments.

5. **Hydrographs.** SCS Unit Hydrograph (Gamma form, PRF 100–600) per
   catchment / sub-basin / junction / reach; CSVs plus interactive
   Plotly HTML visualisations in a `Modeling_Outputs/` tree.

6. **Reach routing.** Lag, Muskingum (K from reach Tc), or
   variable-parameter **Muskingum-Cunge** on a trapezoidal channel
   (Manning n, bottom width, side slope; length & slope from the
   delineated reach; K and X recomputed every time step; the reach
   subdivided into N = L/(cΔt) segments, ≤ 40).

7. **HEC-HMS export.** `EXPORT_HEC_HMS = True` writes
   `<OUTPUT_DIR>/HEC_HMS/AR_Tarseem.basin`: Subbasins (Area, SCS CN,
   SCS lag), Reaches (matching your routing settings), Junctions and
   Sinks at real map coordinates, for HEC-HMS 4.7–4.13 (correct
   version-dependent headers, verified against a genuine HEC-HMS 4.12
   basin file). Open via *File → Import → Basin Model*, or in WMS.
   Reservoirs are modelled downstream in HEC-HMS after export.

## Compatibility

One package, both QGIS generations:

| QGIS | Qt / bindings | Status |
|------|---------------|--------|
| 3.22 → 3.4x | Qt 5 / PyQt5 | Supported; the run executes as a background `QgsTask` and QGIS stays responsive. |
| 4.0 → 4.2 LTR and newer | Qt 6 / PyQt6 | Supported; the run executes on the main thread (see *Run modes* below). |

`qt_compat.py` resolves every Qt enumeration in both its scoped
(PyQt6: `Qt.CheckState.Checked`) and unscoped (PyQt5: `Qt.Checked`)
spelling, and builds attribute-table field types from `QMetaType.Type`
on Qt 6 and `QVariant.Type` on Qt 5. No separate QGIS 4 build exists or
is needed.

### Run modes

Under Qt 6, Qt enforces object thread affinity strictly enough that
handing it layers built on a worker thread can abort the process — a
hard crash with no Python traceback. Because the engine (`hydro_core`)
is deliberately never modified, its layers cannot be built on the GUI
thread instead, so on Qt 6 the whole run executes on the main thread:
QGIS is busy until it finishes, but it cannot crash. The progress bar
is live in both modes.

To override the automatic choice (for testing only):

```
AR_TARSEEM_BACKGROUND=1   # force the background task
AR_TARSEEM_BACKGROUND=0   # force foreground execution
```

or set the QSettings key `AR-Tarseem/background_run` to `always` /
`never`.

## Install

QGIS → *Plugins → Manage and Install Plugins → Install from ZIP* →
select the plugin zip — or install it from the official QGIS plugin
repository, where AR-Tarseem is published under the GPL.

## Files

- `hydro_core.py` — delineation + hydrology engine; `run(cfg)` entry
  point, `CONFIG` retained for QGIS-console use.
- `scs_data.py` — SCS distribution tables (4 × 1441 points) and the
  dimensionless UH reference curve.
- `hec_hms_export.py` — the .basin writer.
- `hydro_dialog.py` — tabbed dialog (Inputs / Delineation / Tc &
  losses / Rainfall & hydrographs / Routing & ARF / HEC-HMS /
  Styling / Output names / About).
- `hydro_plugin.py`, `__init__.py`, `metadata.txt`, `icon.svg`.

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

The hydrograph engine reproduces HEC-HMS numerics: SCS-CN losses,
the HEC-HMS dimensionless unit hydrograph (Tp = Δt/2 + Tlag, unit
volume normalised). To compare one-to-one, set the HEC-HMS control
interval equal to `HYD_TIME_STEP_MIN` — peak discharges of fast
basins are sensitive to the computation interval (a 1-min run reads
several percent higher than a 10-min run on short-lag basins).

## Progress & logging

On QGIS 3 the run is a background `QgsTask`: QGIS stays responsive and
the native task bar reports progress. On QGIS 4 the run is synchronous
and the QGIS window is busy until it finishes. Either way a live
progress bar reports each stage, and progress messages appear in
*View → Panels → Log Messages → "AR-Tarseem"*.
