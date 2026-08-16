# ==========================================================================
# AR-Tarseem - QGIS Plugin
# Version 1.0.0 - July 2026
# ==========================================================================
# Copyright (C) 2026 Ahmed M. Raslan, PhD & Amr R. Refaiy, PhD - Water Engineering Experts
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-3.0-or-later
# ==========================================================================
"""HEC-HMS .basin export (v19).

Copyright (c) 2026 Ahmed M. Raslan, PhD & Amr R. Refaiy, PhD - Water Engineering Experts

Writes the delineated network as a HEC-HMS ``.basin`` file in the exact
element format HEC-HMS / WMS reads. The writer was cross-checked against a genuine HEC-HMS 4.12 basin
file:

* header block with version-dependent lines — 4.11+ writes
  ``Unregulated Output Required``; earlier versions write
  ``Enable Quality Routing``;
* ``Subbasin`` elements with SCS loss (Curve Number) and SCS transform
  (Lag in minutes, standard unitgraph);
* ``Reach`` elements for routed links — ``Lag`` / ``Muskingum`` /
  ``Muskingum-Cunge`` (trapezoid channel) matching the plugin's
  REACH_ROUTING_METHOD;
* ``Junction`` elements at sub-basin outlets;
* one ``Sink`` per catchment at the pour point;
* the trailing schematic/layer blocks HEC-HMS expects.

Canvas coordinates are the real map coordinates from the delineation
(same convention as HEC-HMS files produced from georeferenced GIS
data), so the network lays out geographically when opened.
"""

from __future__ import annotations

import datetime
import os
import math


_MONTHS = ('January', 'February', 'March', 'April', 'May', 'June',
           'July', 'August', 'September', 'October', 'November',
           'December')


def _stamp():
    now = datetime.datetime.now()
    dstr = f"{now.day} {_MONTHS[now.month - 1]} {now.year}"
    tstr = now.strftime("%H:%M:%S")
    return (f"     Last Modified Date: {dstr}\n"
            f"     Last Modified Time: {tstr}\n")


def _sanitize(name):
    return str(name).strip().replace(' ', '_') or 'Element'


def _version_number(version):
    parts = str(version).split('.')
    try:
        return int(parts[0]) * 100 + int(parts[1])
    except (ValueError, IndexError):
        return 412


def _header(basin_name, version):
    vnum = _version_number(version)
    if vnum >= 411:
        extra = ("     Unregulated Output Required: No\n"
                 "\n"
                 "     Enable Sediment Routing: No\n")
    else:
        extra = ("\n"
                 "     Enable Sediment Routing: No\n"
                 "\n"
                 "     Enable Quality Routing: No\n")
    return (f"Basin: {basin_name}\n"
            f"{_stamp()}"
            f"     Version: {version}\n"
            f"     Filepath Separator: \\\n"
            f"     Unit System: Metric\n"
            f"     Missing Flow To Zero: No\n"
            f"     Enable Flow Ratio: No\n"
            f"     Compute Local Flow At Junctions: No\n"
            f"{extra}"
            f"End:\n\n")


# HEC-HMS offers a fixed menu of SCS peak rate factors. Anything else
# has to be snapped to the nearest one on export.
HMS_PRF_OPTIONS = (
    (484.0, 'STANDARD'),          # Standard (PRF 484)
    (284.0, 'DELMARVA'),          # Delmarva (PRF 284)
    (100.0, 'PRF100'), (150.0, 'PRF150'), (200.0, 'PRF200'),
    (250.0, 'PRF250'), (300.0, 'PRF300'), (350.0, 'PRF350'),
    (400.0, 'PRF400'), (450.0, 'PRF450'), (500.0, 'PRF500'),
    (550.0, 'PRF550'), (600.0, 'PRF600'),
)


def hms_unitgraph_type(prf):
    """v1.0 — (keyword, prf_used) for the HEC-HMS 'Unitgraph Type' line.

    Returns the exact keyword when the peak rate factor is one HEC-HMS
    supports, otherwise the nearest available one, so the exported file
    is always valid. The caller reports any substitution."""
    try:
        p = float(prf)
    except (TypeError, ValueError):
        return 'STANDARD', 484.0
    for v, kw in HMS_PRF_OPTIONS:
        if abs(p - v) < 0.5:
            return kw, v
    v, kw = min(HMS_PRF_OPTIONS, key=lambda o: abs(o[0] - p))
    return kw, v


def _subbasin(name, x, y, area_km2, cn, lag_min, downstream, prf=484.0):
    return (f"Subbasin: {name}\n"
            f"{_stamp()}"
            f"     Canvas X: {x:.4f}\n"
            f"     Canvas Y: {y:.4f}\n"
            f"     Area: {area_km2:.4f}\n"
            f"     Downstream: {downstream}\n"
            f"\n"
            f"     Discretization: None\n"
            f"\n"
            f"     Canopy: None\n"
            f"     Allow Simultaneous Precip Et: No\n"
            f"     Plant Uptake Method: None\n"
            f"\n"
            f"     Surface: None\n"
            f"\n"
            f"     LossRate: SCS\n"
            f"     Percent Impervious Area: 0.0\n"
            f"     Curve Number: {cn:.1f}\n"
            f"\n"
            f"     Transform: SCS\n"
            f"     Lag: {lag_min:.2f}\n"
            # v1.0 — the SCS peak rate factor, in HEC-HMS's own
            # vocabulary: STANDARD (484), DELMARVA (284) or PRF<n>.
            f"     Unitgraph Type: {hms_unitgraph_type(prf)[0]}\n"
            f"\n"
            f"     Baseflow: None\n"
            f"End:\n\n")


def _reach(name, x, y, downstream, routing):
    """routing: dict(method=..., plus method-specific keys)."""
    body = (f"Reach: {name}\n"
            f"{_stamp()}"
            f"     Canvas X: {x:.4f}\n"
            f"     Canvas Y: {y:.4f}\n"
            f"     Downstream: {downstream}\n")
    m = (routing.get('method') or 'NONE').upper()
    if m == 'LAG':
        body += ("     Route: Lag\n"
                 "     Lag Method: Specified Lag\n"
                 f"     Lag: {float(routing.get('lag_min') or 0.0):.2f}\n")
    elif m == 'MUSKINGUM':
        body += ("     Route: Muskingum\n"
                 f"     Muskingum K: {float(routing.get('K_hr') or 0.0):.4f}\n"
                 f"     Muskingum X: {float(routing.get('X') or 0.2):.3f}\n"
                 f"     Subreaches: "
                 f"{max(1, int(routing.get('subreaches') or 1))}\n")
    elif m == 'MUSKINGUM_CUNGE':
        body += ("     Route: Muskingum-Cunge\n"
                 "     Channel: Trapezoid\n"
                 f"     Length: {float(routing.get('L_m') or 0.0):.2f}\n"
                 f"     Energy Slope: {float(routing.get('S0') or 0.0):.5f}\n"
                 f"     Bottom Width: {float(routing.get('b_m') or 0.0):.2f}\n"
                 f"     Side Slope: {float(routing.get('z') or 0.0):.2f}\n"
                 f"     Mannings n: {float(routing.get('n') or 0.0):.4f}\n")
    # 'NONE' → no Route lines: HEC-HMS treats it as no-routing reach; we
    # skip emitting the element entirely upstream of this function.
    body += "End:\n\n"
    return body


def _junction(name, x, y, downstream):
    return (f"Junction: {name}\n"
            f"{_stamp()}"
            f"     Canvas X: {x:.4f}\n"
            f"     Canvas Y: {y:.4f}\n"
            f"     Downstream: {downstream}\n"
            f"End:\n\n")


def _res_table_name(index):
    """v1.0 — Elevation-Storage table name.

    HEC-HMS names paired-data tables plainly ('Table 1', 'Table 2', …)
    and reuses that name inside the DSS pathname, where brackets and
    other punctuation are not valid. Naming the table after the
    element (e.g. 'W1-J10(aaa)_ES') therefore produced a table HMS
    could not resolve, leaving the reservoir without a routing curve.
    """
    return f"Table {int(index)}"


def _reservoir(name, x, y, downstream, res, table_index=1):
    """v24 — Reservoir element, AR-Awdiya-compatible: Route: Controlled
    Outflow with an Elevation-Storage table reference (table lives in
    the companion .pdata) plus inline Spillway / Conduit blocks when
    the reservoir is defined by outflow structures."""
    tbl = _res_table_name(table_index)
    _init_elev = None
    if str((res or {}).get('init_mode', 'inflow')) == 'elevation':
        try:
            _init_elev = float((res or {}).get('init_elev'))
        except (TypeError, ValueError):
            _init_elev = None
    body = (f"Reservoir: {name}\n"
            f"{_stamp()}"
            f"     Canvas X: {x:.4f}\n"
            f"     Canvas Y: {y:.4f}\n"
            f"     Downstream: {downstream}\n"
            f"\n"
            f"     Route: Controlled Outflow\n"
            f"     Routing Curve: Elevation-Storage\n"
            # v1.0 — initial condition: either a starting water level
            # ('Initial Elevation') or the default assumption that the
            # first outflow equals the first inflow. HEC-HMS expects
            # exactly one of the two lines, in this position.
            + (f"     Initial Elevation: {float(_init_elev):g}\n"
               if _init_elev is not None
               else "     Initial Outflow Equals Inflow: Yes\n")
            + f"     Elevation-Storage Table: {tbl}\n"
            f"     Adaptive Control: On\n"
            f"     Main Tailwater Condition: None\n"
            f"     Auxiliary Tailwater Condition: None\n")
    if (res or {}).get('mode') == 'structures':
        # conduits first, then spillways — the order HEC-HMS itself
        # writes, so exported files diff cleanly against GUI-made ones
        sts = [s for s in (res.get('structures') or []) if s and s.get('type')]
        for st in sorted(sts, key=lambda s: 0 if s['type'] == 'pipe' else 1):
            if st['type'] == 'pipe':
                # Area and CENTRE elevation follow the outlet shape:
                #   circular    A = pi*r^2 , centre = invert + r
                #   rectangular A = S*H    , centre = invert + H/2
                if str(st.get('shape', 'circle')).lower().startswith('rect'):
                    s_w = float(st.get('span', 0.0))
                    s_h = float(st.get('height', 0.0))
                    area = s_w * s_h
                    centre = float(st.get('invert', 0.0)) + 0.5 * s_h
                else:
                    r = float(st.get('radius', 0.0))
                    area = math.pi * r * r
                    centre = float(st.get('invert', 0.0)) + r
                body += (f"\n     Conduit: Orifice\n"
                         f"     Conduit Outlet: Main\n"
                         f"     Orifice Coefficient: "
                         f"{float(st.get('Cd', 0.6) or 0.6)}\n"
                         f"     Orifice Area: {area:.4f}\n"
                         f"     Centerline Elevation: {centre:.2f}\n"
                         f"     Number Barrels: "
                         f"{max(1, int(round(float(st.get('n', 1) or 1))))}\n"
                         f"     End Conduit:\n")
            elif st['type'] == 'spillway':
                body += (f"\n     Spillway: Broad-Crested Spillway\n"
                         f"     Spillway Outlet: Main\n"
                         f"     Spillway Crest Length: "
                         f"{float(st.get('width', 0.0))}\n"
                         f"     Spillway Crest Elevation: "
                         f"{float(st.get('crest', 0.0))}\n"
                         f"     Spillway Coefficient: "
                         f"{float(st.get('C', 1.7))}\n"
                         f"     End Spillway:\n")
    body += ("\n     Evaporation Method: Zero Evaporation\n"
             "     End Evaporation:\n"
             "End:\n\n")
    return body


def _sink(name, x, y):
    return (f"Sink: {name}\n"
            f"{_stamp()}"
            f"     Canvas X: {x:.4f}\n"
            f"     Canvas Y: {y:.4f}\n"
            f"End:\n\n")


_TAIL = ("Basin Layer Properties:\n"
         "     Element Layer:\n"
         "          Name: Icons\n"
         "          Layer shown: Yes\n"
         "     End Layer:\n"
         "End:\n"
         "\n"
         "Basin Spatial Properties:\n"
         "End:\n"
         "\n"
         "Basin Schematic Properties:\n"
         "     Extent Method: Manual\n"
         "     Buffer: 0\n"
         "     Draw Icons: Yes\n"
         "     Draw Icon Labels: Name\n"
         "     Draw Map Objects: No\n"
         "     Draw Gridlines: No\n"
         "End:\n")


def write_basin_file(path, basin_name, catchments, version='4.12',
                     prf=484.0):
    """Write one .basin file covering all delineated catchments.

    Parameters
    ----------
    path : str
        Output .basin path.
    basin_name : str
        HEC-HMS basin-model name.
    catchments : list of dict
        One entry per pour point, as collected by hydro_core.run():
          pid, name, outlet_xy=(x, y), area_km2, cn, lag_min,
          subs = list of dict(sid, name, ds_sid, area_km2, cn, lag_min,
                              junction_xy, dam_name=None)
                 — ``dam_name`` (v23), when set, is appended to the
                 junction element name: ``W{pid}-J{sid}(DamName)``.
          reaches = list of dict(sid, ds_sid, name, routing=dict,
                                 mid_xy)
        When ``subs`` is empty the catchment exports as a single
        Subbasin → Sink pair.
    version : str
        HEC-HMS version string ('4.7' … '4.13').
    """
    s = _header(_sanitize(basin_name), version)
    subs_txt, reach_txt, junc_txt, sink_txt = [], [], [], []
    _pdata_tables = []          # v24: (element_name, reservoir dict)

    for c in catchments:
        pid = c['pid']
        ox, oy = c['outlet_xy']
        sink_name = _sanitize(f"Outlet_W{pid}")

        if not c.get('subs'):
            # Lumped: single subbasin straight into the sink. Icon at
            # the catchment centroid (v22.1), sink at the pour point.
            cx, cy = c.get('centroid_xy') or (ox, oy)
            subs_txt.append(_subbasin(
                _sanitize(c['name']), cx, cy,
                float(c['area_km2'] or 0.0),
                float(c['cn'] or 0.0),
                float(c['lag_min'] or 0.0),
                sink_name, prf))
            sink_txt.append(_sink(sink_name, ox, oy))
            continue

        # Subdivided: Subbasin -> its junction; junction -> reach ->
        # downstream junction; outlet junction (sid 1) -> sink.
        sub_by_sid = {s_['sid']: s_ for s_ in c['subs']}
        reach_by_sid = {r_['sid']: r_ for r_ in (c.get('reaches') or [])}

        def junc_name(sid):
            # v23 — junctions enforced by a reservoir carry the dam
            # name, matching the QGIS Junctions layer:
            # "W1-J9(DamName)". Used consistently for the junction
            # element itself AND as every Downstream reference to it.
            # v1.0 naming: W{pid}-J{sid}, e.g. "W1-J9(DamName)"
            base = f"W{pid}-J{sid}"
            dn = (sub_by_sid.get(sid) or {}).get('dam_name')
            return _sanitize(f"{base}({dn})" if dn else base)

        for s_ in c['subs']:
            sid = s_['sid']
            jx, jy = s_['junction_xy']
            # v22.1: Subbasin icon at the sub-basin CENTROID; only the
            # junction element stays on the outlet point.
            scx, scy = s_.get('centroid_xy') or (jx, jy)
            subs_txt.append(_subbasin(
                _sanitize(s_['name']), scx, scy,
                float(s_['area_km2'] or 0.0),
                float(s_['cn'] or 0.0),
                float(s_['lag_min'] or 0.0),
                junc_name(sid), prf))

            # junction downstream target
            if sid == 1:
                ds_target = sink_name
            else:
                r_ = reach_by_sid.get(sid)
                ds_sid = s_['ds_sid']
                ds_junc = (junc_name(ds_sid) if ds_sid in sub_by_sid
                           else sink_name)
                if r_ is not None and \
                        (r_['routing'].get('method') or 'NONE') != 'NONE':
                    # v1.0 naming: W{pid}-R{sid}
                    rname = _sanitize(f"W{pid}-R{sid}")
                    mx, my = r_.get('mid_xy') or s_['junction_xy']
                    reach_txt.append(_reach(rname, mx, my, ds_junc,
                                            r_['routing']))
                    ds_target = rname
                else:
                    ds_target = ds_junc
            _res = s_.get('reservoir')
            if _res:
                _tbl_i = len(_pdata_tables) + 1
                junc_txt.append(_reservoir(junc_name(sid), jx, jy,
                                           ds_target, _res, _tbl_i))
                _pdata_tables.append((junc_name(sid), _res, _tbl_i))
            else:
                junc_txt.append(_junction(junc_name(sid), jx, jy,
                                          ds_target))

        sink_txt.append(_sink(sink_name, ox, oy))

    s += "".join(subs_txt) + "".join(reach_txt) \
         + "".join(junc_txt) + "".join(sink_txt) + _TAIL

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\r\n") as f:
        f.write(s)

    # v24 — companion .pdata (Paired Data Manager) with each
    # reservoir's Elevation-Storage table, the format HEC-HMS
    # actually reads (and AR-Awdiya imports/exports). Written next
    # to the .basin with the same stem.
    if _pdata_tables:
        pdata_path = os.path.splitext(path)[0] + '.pdata'
        write_pdata_file(pdata_path, _sanitize(basin_name),
                         _pdata_tables, version)
        print(f"    -> reservoir Elevation-Storage table(s) written "
              f"to {pdata_path}")
        _csvs = write_paired_data_csvs(os.path.dirname(path),
                                       _pdata_tables)
        if _csvs:
            print(f"    -> {len(_csvs) - 1} elevation-storage CSV(s) + "
                  f"PAIRED_DATA_README.txt (HEC-HMS keeps these values "
                  f"in the project DSS)")
    return path


def write_pdata_file(path, project_name, tables, version='4.12'):
    """v1.0 — Paired Data Manager file with one Elevation-Storage table
    per reservoir.

    Written to match a file produced by HEC-HMS itself: plain table
    names ('Table 1', 'Table 2', …) and metadata only. HEC-HMS keeps
    paired-data VALUES in the project DSS file, not in the .pdata, so
    the elevation-storage pairs are exported alongside as CSV (see
    write_paired_data_csvs) for pasting/importing into the Paired Data
    Manager.

    tables : list of (element_name, reservoir_dict, table_index)
    """
    now = datetime.datetime.now()
    dstr = now.strftime('%d %B %Y')
    tstr = now.strftime('%H:%M:%S')
    s = (f"Paired Data Manager: {project_name}\n"
         f"     Version: {version}\n"
         f"     Filepath Separator: \\\n"
         f"End:\n\n")
    for _el_name, _res, idx in tables:
        tbl = _res_table_name(idx)
        s += (f"Table: {tbl}\n"
              f"     Table Type: Elevation-Storage\n"
              f"     Last Modified Date: {dstr}\n"
              f"     Last Modified Time: {tstr}\n"
              f"     X-Units: M\n"
              f"     Y-Units: THOU M3\n"
              f"     Use External DSS File: NO\n"
              f"     DSS File: {project_name}.dss\n"
              f"     Pathname: //{tbl}/ELEVATION-STORAGE///TABLE/\n"
              f"     Interpolation: Linear Interpolation\n"
              f"End:\n\n")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\r\n") as f:
        f.write(s)
    return path


def write_paired_data_csvs(out_dir, tables):
    """v1.0 — one CSV per reservoir holding the Elevation-Storage pairs
    that belong in the HEC-HMS Paired Data Manager, plus a short
    README explaining how to get them in (HMS stores these values in
    the project DSS, so they cannot travel inside the .pdata).

    Returns the list of files written."""
    os.makedirs(out_dir, exist_ok=True)
    written = []
    for el_name, res, idx in tables:
        tbl = _res_table_name(idx)
        rows = []
        if (res or {}).get('ev_rows'):
            rows = [(float(e), float(v) / 1000.0)
                    for e, v in res['ev_rows']]
        else:
            rows = [(float(i), float(s_) / 1000.0)
                    for i, (s_, _o) in enumerate(res.get('sd_rows') or [])]
        p = os.path.join(out_dir,
                         f"{tbl.replace(' ', '_')}_"
                         f"{_sanitize(el_name)}.csv")
        with open(p, 'w', encoding='utf-8', newline='') as fh:
            fh.write("Elevation (M),Storage (THOU M3)\n")
            for x_, y_ in rows:
                fh.write(f"{x_:.2f},{y_:.4f}\n")
        written.append(p)
    if written:
        readme = os.path.join(out_dir, 'PAIRED_DATA_README.txt')
        with open(readme, 'w', encoding='utf-8') as fh:
            fh.write(
                "HEC-HMS reservoir Elevation-Storage tables\n"
                "=========================================\n\n"
                "The .basin file defines each reservoir and points it at\n"
                "a paired-data table (Table 1, Table 2, ...). The .pdata\n"
                "file creates those tables with the correct name, units\n"
                "(M / THOU M3) and interpolation.\n\n"
                "HEC-HMS stores the VALUES of a paired-data table in the\n"
                "project DSS file rather than in the .pdata, so the\n"
                "elevation-storage pairs are supplied here as CSV, one\n"
                "file per reservoir:\n\n")
            for p in written:
                fh.write(f"    {os.path.basename(p)}\n")
            fh.write(
                "\nTo load them: open the project in HEC-HMS, go to\n"
                "Components > Paired Data Manager > Elevation-Storage,\n"
                "select the table, and paste the two columns into the\n"
                "table editor (or import the pairs from the CSV).\n"
                "Everything else - the reservoir, its outlet structures\n"
                "and the network - is already complete.\n")
        written.append(readme)
    return written


# Hypothetical-Storm 'Storm Type' names exactly as HEC-HMS 4.12
# *stores* them in the .met file (cross-checked against a genuine
# HMS 4.12 export): 'SCS Type II' etc. — Roman numerals, not the
# arabic digits shown in some GUI labels.
_STORM_TYPE_MAP = {
    'SCS_TYPE_I':   'SCS Type I',
    'SCS_TYPE_IA':  'SCS Type IA',
    'SCS_TYPE_II':  'SCS Type II',
    'SCS_TYPE_III': 'SCS Type III',
}


def write_met_files(out_dir, basin_name, catchments,
                    storm_distribution='SCS_TYPE_II', version='4.12',
                    storm_duration_hr=24.0):
    """Write one HEC-HMS .met per return period (v22.1).

    Emits the *file* keywords HEC-HMS 4.12 actually stores (verified
    against a genuine HMS 4.12 export), which differ from the GUI
    labels:

    * ``Precipitation Method: Point Depth`` + ``Storm Type: SCS Type
      II`` (Roman numerals) instead of a GUI-style ``Method:`` line;
    * ``Storm Duration: <hr>`` — required by HMS, taken from
      ``storm_duration_hr`` (wire to STORM_DURATION_HR);
    * ``Depth-Area Reduction Method: No Reduction`` (any ARF is
      already baked into the exported depths);
    * ``Uniform Depth: No`` — HMS's stored flag for the GUI's
      'Spatial Distribution: Variable By Subbasin'. Written always,
      with every ``Subbasin`` block carrying its own ``Depth:`` line;
      this is valid whether or not the depths happen to be equal, so
      uniform (non-Thiessen) exports import identically.

    ``catchments`` is the same list handed to :func:`write_basin_file`,
    where every catchment and subbasin dict now carries
    ``depths={rp: mm}``. Returns the list of written paths.
    """
    storm_type = _STORM_TYPE_MAP.get(
        str(storm_distribution or '').upper().strip(), 'SCS Type II')

    # gather (name, area, depths) for every exported subbasin
    subs = []
    rps = set()
    for c in catchments:
        if c.get('subs'):
            for s in c['subs']:
                d = dict(s.get('depths') or c.get('depths') or {})
                subs.append((_sanitize(s['name']),
                             float(s.get('area_km2') or 0.0), d))
                rps.update(d.keys())
        else:
            d = dict(c.get('depths') or {})
            subs.append((_sanitize(c['name']),
                         float(c.get('area_km2') or 0.0), d))
            rps.update(d.keys())
    if not subs or not rps:
        return []

    try:
        dur_hr = float(storm_duration_hr)
    except (TypeError, ValueError):
        dur_hr = 24.0
    # HMS writes whole hours without a decimal point
    dur_txt = (f"{dur_hr:g}")

    paths = []
    for rp in sorted(rps):
        met_name = f"Met_{rp}yr"
        # Meteorology block — full method list exactly as HMS 4.12
        # writes it (incl. 'Shortwave'/'Longwave', one word).
        s = (f"Meteorology: {met_name}\n"
             f"{_stamp()}"
             f"     Version: {version}\n"
             f"     Unit System: Metric\n"
             f"     Set Missing Data to Default: No\n"
             f"     Precipitation Method: Hypothetical Storm\n"
             f"     Air Temperature Method: None\n"
             f"     Atmospheric Pressure Method: None\n"
             f"     Dew Point Method: None\n"
             f"     Wind Speed Method: None\n"
             f"     Shortwave Radiation Method: None\n"
             f"     Longwave Radiation Method: None\n"
             f"     Snowmelt Method: None\n"
             f"     Evapotranspiration Method: No Evapotranspiration\n"
             f"     Use Basin Model: {_sanitize(basin_name)}\n"
             f"End:\n\n")
        # Precip parameters block — stored keywords, not GUI labels.
        s += (f"Precip Method Parameters: Hypothetical Storm\n"
              f"{_stamp()}"
              f"     Precipitation Method: Point Depth\n"
              f"     Storm Type: {storm_type}\n"
              f"     Storm Duration: {dur_txt}\n"
              f"     Depth-Area Reduction Method: No Reduction\n"
              f"     Uniform Depth: No\n"
              f"End:\n\n")
        # One block per subbasin, each with its own Depth (mm).
        for name, _a, d in subs:
            s += (f"Subbasin: {name}\n"
                  f"{_stamp()}")
            if rp in d:
                s += f"     Depth: {float(d[rp]):.2f}\n"
            s += "End:\n\n"
        path = os.path.join(out_dir, f"{met_name}.met")
        os.makedirs(out_dir, exist_ok=True)
        with open(path, "w", encoding="utf-8", newline="\r\n") as fh:
            fh.write(s)
        paths.append(path)
    return paths
