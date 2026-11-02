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
"""AR-Awdiya project export (.json).

Writes the delineated network as a complete AR-Awdiya project file that
opens with "Load project" and is ready to Run — no manual re-entry.

Element mapping
---------------
    AR-Tarseem                     AR-Awdiya
    ---------------------------    -------------------------------------
    sub-basin  W{p}-Sub{s}     ->  catchment  (cats[])   drains to its
                                   own junction, direct connection
    junction   W{p}-J{s}       ->  junction   (juncs[])  carries the
                                   REACH routing of that sub-basin on
                                   the link to the downstream junction
    reservoir  W{p}-J{s}(Dam)  ->  junction with kind='reservoir' plus
                                   the elevation-volume curve / outflow
                                   structures or storage-discharge table
    catchment outlet W{p}      ->  outlet     (outs[])

Every catchment (pour point) becomes its own outlet, so several pour
points export as several independent networks in one project.

Hydrology carried over
----------------------
    area (km2), Curve Number, lag/Tc (min), peak rate factor, lambda,
    storm distribution + duration, per-sub-basin design depth as the
    UNREDUCED point depth together with the areal reduction factor in
    the catchment's ARF field (AR-Awdiya then applies P x ARF itself,
    reproducing the effective depth without double counting), reach
    routing method and its parameters (Lag / Muskingum K,X /
    Muskingum-Cunge geometry), and the reservoir definitions.

Return periods
--------------
    'single' : one storm; each catchment keeps its OWN depth, so the
               Thiessen spatial variation of rainfall is preserved.
    'all'    : AR-Awdiya frequency analysis; the depth list is written
               to the frequency table and AR-Awdiya then applies each
               depth uniformly to all catchments (its own convention).
"""

from __future__ import annotations

import datetime
import json
import math
import os

PROJECT_APP = 'SCS_Watershed_Model'
PROJECT_VERSION = 1

# canvas geometry used by AR-Awdiya's network view
_CANVAS_X0, _CANVAS_X1 = 90.0, 860.0
_ROW_TOP, _ROW_STEP, _ROW_MAX = 60.0, 78.0, 470.0

_STORM_MAP = {
    'SCS_TYPE_I': 'typeI', 'SCS_TYPE_IA': 'typeIA',
    'SCS_TYPE_II': 'typeII', 'SCS_TYPE_III': 'typeIII',
    'CUSTOM': 'custom',
}


def _reach(routing, tc_min):
    """AR-Tarseem routing dict -> AR-Awdiya reach object."""
    r = dict(method='none', lagMin=30, K=2, X=0.2,
             L=2000, S0=0.005, n=0.035, b=5, z=2)
    if not routing:
        return r
    m = str(routing.get('method', '')).upper()
    if m == 'LAG':
        r['method'] = 'lag'
        r['lagMin'] = round(float(routing.get('lag_min')
                                  or tc_min or 0.0), 4)
    elif m == 'MUSKINGUM':
        r['method'] = 'musk'
        r['K'] = round(float(routing.get('K_hr')       # hours, as AR-Awdiya
                             or (tc_min or 0.0) / 60.0), 5)
        r['X'] = float(routing.get('X') or 0.2)
    elif m in ('MUSKINGUM_CUNGE', 'MC'):
        r['method'] = 'mc'
        r['L'] = round(float(routing.get('length_m') or 2000.0), 2)
        r['S0'] = float(routing.get('slope') or 0.005) or 0.005
        r['n'] = float(routing.get('manning_n') or 0.035)
        r['b'] = float(routing.get('bottom_width_m') or 5.0)
        r['z'] = float(routing.get('side_slope') or 2.0)
    return r


# rectangular outlets converted for AR-Awdiya, reported in warning.txt
_RECT_NOTES = []
# reservoirs whose initial water level cannot travel to AR-Awdiya
_INIT_NOTES = []


def _reservoir_fields(res, owner_name=''):
    """Reservoir definition -> AR-Awdiya junction fields."""
    if not res:
        return {}
    # AR-Awdiya always starts a reservoir empty (at the curve datum);
    # it has no initial-water-level option, so record the loss.
    if str(res.get('init_mode', 'inflow')) == 'elevation' \
            and res.get('init_elev') is not None:
        _INIT_NOTES.append((owner_name, float(res['init_elev'])))
    mode = str(res.get('mode') or 'structures')
    out = {}
    if mode == 'structures' and res.get('ev_rows'):
        # elevation (m), volume (10^3 m3) — AR-Awdiya's row format
        out['resMode'] = 'structures'
        out['evText'] = '\n'.join(
            f"{float(e):g}, {float(v) / 1000.0:g}"
            for e, v in res['ev_rows'])
        sts = []
        for st in (res.get('structures') or []):
            if not st or not st.get('type'):
                continue
            st = dict(st, _owner=owner_name)
            if st['type'] == 'spillway':
                sts.append(dict(type='spillway',
                                crest=float(st.get('crest', 0.0)),
                                C=float(st.get('C', 1.7)),
                                width=float(st.get('width', 5.0))))
            else:
                # AR-Awdiya models orifice outlets as CIRCULAR only.
                # A rectangular outlet is therefore exported with an
                # equivalent radius of H/2, which preserves the centre
                # elevation (invert + H/2) that governs the head; the
                # flow area differs (pi*(H/2)^2 instead of S*H), so the
                # substitution is recorded in warning.txt.
                if str(st.get('shape', 'circle')).lower().startswith(
                        'rect'):
                    _h = float(st.get('height', 0.0))
                    _s = float(st.get('span', 0.0))
                    _r = 0.5 * _h
                    _RECT_NOTES.append(
                        (str(st.get('_owner') or ''), _s, _h, _r))
                else:
                    _r = float(st.get('radius', 0.45))
                sts.append(dict(type='pipe',
                                invert=float(st.get('invert', 0.0)),
                                n=int(st.get('n', 1) or 1),
                                radius=_r,
                                Cd=float(st.get('Cd', 0.6)),
                                Cw=float(st.get('Cw', 1.5))))
        out['structures'] = sts
    else:
        # storage (10^3 m3), outflow (m3/s)
        out['resMode'] = 'storage'
        out['storageText'] = '\n'.join(
            f"{float(s) / 1000.0:g}, {float(o):g}"
            for s, o in (res.get('sd_rows') or []))
    return out


def _point_depth(obj, rp, fallback=None):
    """Un-reduced point depth for a return period: the payload carries
    ARF-reduced depths (as HEC-HMS needs them), so divide the factor
    back out and let AR-Awdiya re-apply it through the ARF field."""
    d = (obj.get('depths') or {}).get(rp)
    if d is None and fallback is not None:
        d = (fallback.get('depths') or {}).get(rp)
        obj = fallback
    if d is None:
        return None
    a = float(obj.get('arf') or 1.0)
    return float(d) / (a if a > 1e-9 else 1.0)


def _arf_text(obj):
    """ARF as AR-Awdiya expects it: a plain multiplier string, or ''
    when no areal reduction was applied."""
    a = float(obj.get('arf') or 1.0)
    return '' if abs(a - 1.0) < 1e-9 else f"{a:.4f}"


def _layout(catchments):
    """Assign canvas coordinates: real easting sets the column order,
    topological depth (hops to the outlet) sets the row, so the network
    reads top-to-bottom exactly the way AR-Awdiya draws it."""
    xs = []
    for c in catchments:
        for s in c.get('subs') or []:
            jxy = s.get('junction_xy') or (0.0, 0.0)
            cxy = s.get('centroid_xy') or jxy
            xs.extend([float(jxy[0]), float(cxy[0])])
        xs.append(float((c.get('outlet_xy') or (0.0, 0.0))[0]))
    if not xs:
        xs = [0.0, 1.0]
    x_min, x_max = min(xs), max(xs)
    span = (x_max - x_min) or 1.0

    def sx(x):
        f = (float(x) - x_min) / span
        return round(_CANVAS_X0 + f * (_CANVAS_X1 - _CANVAS_X0), 1)

    def sy(depth):
        return round(min(_ROW_TOP + depth * _ROW_STEP, _ROW_MAX), 1)

    return sx, sy


def _depth_of(subs_by_sid, sid, cache):
    """Hops from sub-basin `sid` down to the outlet sub-basin."""
    if sid in cache:
        return cache[sid]
    seen, d, cur = set(), 0, sid
    while cur and cur in subs_by_sid and cur not in seen:
        seen.add(cur)
        nxt = int(subs_by_sid[cur].get('ds_sid') or 0)
        if nxt <= 0:
            break
        d += 1
        cur = nxt
    cache[sid] = d
    return d


def build_awdiya_project(catchments,
                         project_name='AR-Tarseem model',
                         dt_min=1.0,
                         total_hr=None,
                         storm_distribution='SCS_TYPE_II',
                         storm_duration_hr=24.0,
                         prf=484.0,
                         scs_lambda=0.2,
                         freq_mode='all',
                         primary_rp=None,
                         hyeto_text='',
                         hms_version='4.11'):
    """Return the AR-Awdiya project dict (ready for json.dump)."""
    del _RECT_NOTES[:]          # per-export collectors
    del _INIT_NOTES[:]
    cats, juncs, outs = [], [], []
    uid = juid = ouid = 0
    sxf, syf = _layout(catchments)

    # return periods present (sorted); the 'primary' one fills each
    # catchment's own P value
    rps = []
    for c in catchments:
        for s in (c.get('subs') or []):
            rps.extend(int(r) for r in (s.get('depths') or {}).keys())
        rps.extend(int(r) for r in (c.get('depths') or {}).keys())
    rps = sorted(set(rps))
    if primary_rp is None:
        primary_rp = rps[-1] if rps else None

    rain_mode = _STORM_MAP.get(str(storm_distribution).upper(), 'typeII')
    dur_min = int(round(float(storm_duration_hr) * 60.0))

    for c in catchments:
        pid = c.get('pid')
        subs = list(c.get('subs') or [])
        reach_by_sid = {int(r.get('sid')): r
                        for r in (c.get('reaches') or [])
                        if r.get('sid') is not None}

        # ---- outlet -------------------------------------------------
        ouid += 1
        o_key = f'O{ouid}'
        oxy = c.get('outlet_xy') or (0.0, 0.0)
        outs.append(dict(id=ouid,
                         name=f"Outlet {c.get('name') or pid}",
                         nx=sxf(oxy[0]), ny=_ROW_MAX + 40.0))

        if not subs:
            # no subdivision: the whole catchment is one AR-Awdiya
            # catchment discharging straight into its outlet
            uid += 1
            dep = _point_depth(c, primary_rp)
            cats.append(dict(
                id=uid, name=str(c.get('name') or f'W{pid}'),
                A=round(float(c.get('area_km2') or 0.0), 6),
                CN=round(float(c.get('cn') or 0.0), 2),
                P=round(float(dep), 3) if dep is not None else 0.0,
                type='Lag', val=round(float(c.get('lag_min') or 0.0), 4),
                rainMode=rain_mode, hyetText=hyeto_text,
                hyetCumulative=False, durMin=dur_min,
                prf=float(prf), lambda_=float(scs_lambda),
                arf=_arf_text(c), ds=o_key,
                reach=_reach(None, None),
                nx=sxf((c.get('centroid_xy') or oxy)[0]), ny=_ROW_TOP))
            # 'lambda' is a python keyword — fix the key afterwards
            cats[-1]['lambda'] = cats[-1].pop('lambda_')
            continue

        subs_by_sid = {int(s['sid']): s for s in subs}
        depth_cache = {}
        j_key_of = {}

        # ---- one junction per sub-basin -----------------------------
        for s in subs:
            sid = int(s['sid'])
            juid += 1
            j_key_of[sid] = f'J{juid}'
            jxy = s.get('junction_xy') or (0.0, 0.0)
            dam = s.get('dam_name')
            j = dict(id=juid,
                     name=(f"W{pid}-J{sid}" + (f"({dam})" if dam else "")),
                     kind='reservoir' if s.get('reservoir') else 'junction',
                     ds=o_key,                       # fixed in second pass
                     reach=_reach(reach_by_sid.get(sid, {}).get('routing'),
                                  s.get('lag_min')),
                     nx=sxf(jxy[0]),
                     ny=syf(_depth_of(subs_by_sid, sid, depth_cache) + 1))
            j.update(_reservoir_fields(s.get('reservoir'),
                                       j['name']))
            juncs.append(j)

        # link each junction to the downstream junction (or the outlet)
        for s, j in zip(subs, juncs[-len(subs):]):
            ds_sid = int(s.get('ds_sid') or 0)
            j['ds'] = j_key_of.get(ds_sid, o_key) if ds_sid > 0 else o_key

        # ---- one catchment per sub-basin ----------------------------
        for s in subs:
            sid = int(s['sid'])
            uid += 1
            dep = _point_depth(s, primary_rp, fallback=c)
            cxy = s.get('centroid_xy') or s.get('junction_xy') or (0.0, 0.0)
            cat = dict(
                id=uid, name=str(s.get('name') or f'W{pid}-Sub{sid}'),
                A=round(float(s.get('area_km2') or 0.0), 6),
                CN=round(float(s.get('cn') or 0.0), 2),
                P=round(float(dep), 3) if dep is not None else 0.0,
                type='Lag', val=round(float(s.get('lag_min') or 0.0), 4),
                rainMode=rain_mode, hyetText=hyeto_text,
                hyetCumulative=False, durMin=dur_min,
                prf=float(prf), arf=_arf_text(s),
                ds=j_key_of[sid],          # sub-basin -> its own junction
                reach=_reach(None, None),  # direct connection
                nx=sxf(cxy[0]),
                ny=syf(_depth_of(subs_by_sid, sid, depth_cache)))
            cat['lambda'] = float(scs_lambda)
            cats.append(cat)

    # ---- global settings --------------------------------------------
    if total_hr is None:
        total_hr = 30
        try:
            _tl = max([float(s.get('lag_min') or 0.0)
                       for c in catchments for s in (c.get('subs') or [])]
                      or [60.0])
            total_hr = int(max(30, round(float(storm_duration_hr)
                                         + 20.0 * _tl / 60.0)))
        # B110: optional depth lookup; a bad row just omits it
        except Exception:  # nosec B110
            pass

    freq_rows = '\n'.join(
        f"{rp}-yr, {round(float(_uniform_depth(catchments, rp)), 3)}"
        for rp in rps) if rps else ''

    return {
        'app': PROJECT_APP,
        'version': PROJECT_VERSION,
        'saved': datetime.datetime.now().isoformat(timespec='seconds'),
        'source': 'AR-Tarseem v1.1.0 (QGIS)',
        'counters': {'uid': uid, 'juid': juid, 'ouid': ouid},
        'cats': cats, 'juncs': juncs, 'outs': outs,
        'settings': {
            'dtMin': str(dt_min),
            'totalHr': str(total_hr),
            'globalRain': 'per',
            'freqOn': bool(freq_mode == 'all' and rps),
            'freqRows': freq_rows,
            'basinUnit': 'KM2',
            'hmsVersion': str(hms_version),
            'projectName': str(project_name),
        },
    }


def _uniform_depth(catchments, rp):
    """Area-weighted mean POINT depth for a return period (used for
    the frequency table; AR-Awdiya applies it to all catchments and
    multiplies each by that catchment's own ARF)."""
    num = den = 0.0
    for c in catchments:
        for s in (c.get('subs') or []):
            a = float(s.get('area_km2') or 0.0)
            d = _point_depth(s, rp, fallback=c)
            if d is None:
                continue
            num += a * float(d)
            den += a
        if not (c.get('subs') or []):
            a = float(c.get('area_km2') or 0.0)
            d = _point_depth(c, rp)
            if d is not None:
                num += a * float(d)
                den += a
    return (num / den) if den else 0.0


def write_awdiya_project(path, catchments, **kwargs):
    """Build and write the AR-Awdiya .json project. Returns the path.

    When a reservoir uses a RECTANGULAR low-level outlet, AR-Awdiya —
    which supports circular orifices only — receives an equivalent
    circular outlet of radius H/2 and a warning.txt is written next to
    the project describing exactly what was substituted."""
    data = build_awdiya_project(catchments, **kwargs)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(data, fh, indent=1)
    if _RECT_NOTES or _INIT_NOTES:
        wpath = os.path.join(os.path.dirname(path), 'warning.txt')
        parts = ["AR-Awdiya export - WARNINGS\n",
                 "===========================\n\n"]
        n = 0
        if _INIT_NOTES:
            n += 1
            parts.append(
                f"{n}) INITIAL RESERVOIR LEVEL NOT TRANSFERRED\n"
                "-----------------------------------------\n\n"
                "AR-Awdiya always starts a reservoir EMPTY, at the\n"
                "datum of its elevation-volume curve; it has no\n"
                "initial-water-level setting. The reservoir(s) below\n"
                "start from a specified level in AR-Tarseem:\n\n")
            for nm, el in _INIT_NOTES:
                parts.append(f"  * {nm or '(reservoir)'}: initial "
                             f"elevation {el:g} m\n")
            parts.append(
                "\nConsequence: in AR-Awdiya the pool begins with MORE\n"
                "empty flood storage than intended, so the attenuation\n"
                "is over-estimated and the peak outflow under-estimated\n"
                "compared with AR-Tarseem and the HEC-HMS export, both\n"
                "of which start from the level you specified.\n\n"
                "Work-around: raise the datum of the elevation-volume\n"
                "curve in AR-Awdiya to the initial level (drop the rows\n"
                "below it and re-zero the volumes); that reproduces the\n"
                "same starting storage.\n\n")
        if _RECT_NOTES:
            n += 1
            parts.append(
                f"{n}) RECTANGULAR OUTLET(S) CONVERTED TO CIRCULAR\n"
                "---------------------------------------------\n\n"
                "AR-Awdiya defines a reservoir low-level outlet as a\n"
                "CIRCULAR orifice only. The reservoir(s) below use a\n"
                "RECTANGULAR outlet in AR-Tarseem, so an equivalent\n"
                "circular orifice of radius R = H / 2 was exported:\n\n")
            for nm, s_, h_, r_ in _RECT_NOTES:
                parts.append(
                    f"  * {nm or '(reservoir)'}\n"
                    f"      AR-Tarseem : rectangular, span S = {s_:g} m,"
                    f" height H = {h_:g} m (area {s_ * h_:g} m2)\n"
                    f"      AR-Awdiya  : circular, radius R = {r_:g} m "
                    f"(area {math.pi * r_ * r_:.4f} m2)\n\n")
            parts.append(
                "Why R = H / 2: it preserves the CENTRE elevation of\n"
                "the opening (invert + H/2), which sets the head on the\n"
                "orifice and therefore the timing of the outflow. The\n"
                "flow AREA is not preserved: a circle of radius H/2 has\n"
                "area pi*(H/2)^2 instead of S*H, so AR-Awdiya computes\n"
                "a different outlet discharge unless S = pi*H/4.\n\n"
                "If the low-level outlet matters for your results, set\n"
                "the radius in AR-Awdiya to sqrt(S*H/pi) instead, which\n"
                "preserves the AREA rather than the centre elevation.\n\n")
        parts.append(
            "The HEC-HMS export is NOT affected by any of the above:\n"
            "it receives the true outlet geometry (area S*H, centre\n"
            "invert + H/2) and the initial reservoir elevation.\n")
        with open(wpath, 'w', encoding='utf-8') as fh:
            fh.write("".join(parts))
        _what = " and ".join(
            ([f"{len(_INIT_NOTES)} initial reservoir level(s)"]
             if _INIT_NOTES else [])
            + ([f"{len(_RECT_NOTES)} rectangular outlet(s)"]
               if _RECT_NOTES else []))
        print(f"    \u26a0 AR-Awdiya cannot represent {_what} exactly "
              f"- see {wpath}")
    return path
