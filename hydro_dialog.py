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
"""AR-Tarseem dialog.

Copyright (c) 2026 Ahmed M. Raslan, PhD & Amr R. Refaiy, PhD - Water Engineering Experts

Replaces the ``CONFIG`` dict with a tabbed Qt dialog using QGIS-native
widgets. v19+ capabilities: initial-abstraction
ratio λ = Ia/S, all four SCS 24-h storm distributions (Type I / IA /
II / III at 1-min resolution), Muskingum-Cunge reach routing with
trapezoidal channel geometry, and HEC-HMS .basin export.

Execution model: synchronous on the GUI thread (same as the no-install
script); prints are redirected to QgsMessageLog tab "AR-Tarseem".
"""

from __future__ import annotations

import contextlib
import io
import os
import sys
import traceback

from qgis.PyQt.QtCore import Qt
from qgis.PyQt.QtGui import QCursor
from qgis.PyQt.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPlainTextEdit,
    QProgressDialog,
    QPushButton,
    QSpinBox,
    QStackedWidget,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from qgis.core import (
    Qgis,
    QgsApplication,
    QgsMapLayerProxyModel,
    QgsMessageLog,
    QgsProject,
    QgsRasterLayer,
    QgsVectorLayer,
)
from qgis.gui import (
    QgsFieldComboBox,
    QgsFileWidget,
    QgsMapLayerComboBox,
)

from .hydro_task import ARTarseemRunTask
from .hydro_core import (
    CONFIG as DEFAULT_CONFIG,
    _DESIGN_RETURN_PERIODS,
    reservoir_sd,
    run as hydro_run,
)


LOG_TAG = "AR-Tarseem"


# =============================================================================
# stdout redirector — sends print() from hydro_core to QgsMessageLog
# =============================================================================
class _LogStream(io.TextIOBase):
    """Line-buffered stream pushing complete lines to QgsMessageLog."""

    def __init__(self, level: int = Qgis.Info):
        super().__init__()
        self._level = level
        self._buf: list[str] = []

    def writable(self) -> bool:
        return True

    def write(self, s: str) -> int:
        if not s:
            return 0
        self._buf.append(s)
        if "\n" in s:
            joined = "".join(self._buf)
            lines = joined.split("\n")
            self._buf = [lines[-1]] if lines[-1] else []
            for line in lines[:-1]:
                if line:
                    QgsMessageLog.logMessage(line, LOG_TAG, self._level)
        return len(s)

    def flush(self) -> None:
        if self._buf:
            tail = "".join(self._buf)
            self._buf = []
            if tail.strip():
                QgsMessageLog.logMessage(tail, LOG_TAG, self._level)


@contextlib.contextmanager
def _redirect_prints_to_log():
    out_stream = _LogStream(Qgis.Info)
    err_stream = _LogStream(Qgis.Warning)
    old_out, old_err = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = out_stream, err_stream
    try:
        yield
    finally:
        out_stream.flush()
        err_stream.flush()
        sys.stdout, sys.stderr = old_out, old_err


# =============================================================================
# Dialog
# =============================================================================
class HydroDelineationDialog(QDialog):
    """Tabbed dialog mirroring the v19 ``CONFIG`` block."""

    PLUGIN_VERSION = "1.0.0"
    AUTHOR = "Ahmed M. Raslan, PhD & Amr R. Refaiy, PhD - Water Engineering Experts"

    def __init__(self, iface, parent=None):
        super().__init__(parent)
        self.iface = iface
        self._task = None            # v1.0 — background QgsTask handle
        self.setWindowTitle(f"AR-Tarseem v{self.PLUGIN_VERSION}")
        self.setMinimumWidth(700)
        self.setMinimumHeight(560)

        root = QVBoxLayout(self)
        self.tabs = QTabWidget(self)
        root.addWidget(self.tabs)

        self._build_tab_inputs()
        self._build_tab_delineation()
        self._build_tab_reservoirs()
        self._build_tab_tc_loss()
        self._build_tab_rainfall_hydrographs()
        self._build_tab_routing_arf()
        self._build_tab_awdiya()
        self._build_tab_hms()
        self._build_tab_styling()
        self._build_tab_outputs()
        self._build_tab_about()

        btn_row = QDialogButtonBox(self)
        btn_row.addButton("Run", QDialogButtonBox.AcceptRole)
        btn_row.addButton("Close", QDialogButtonBox.RejectRole)
        btn_row.accepted.connect(self.on_run)
        btn_row.rejected.connect(self.reject)
        root.addWidget(btn_row)

        self._wire_dependencies()

    # =========================================================================
    # Tab 1 — Inputs
    # =========================================================================
    def _build_tab_inputs(self):
        w = QWidget(self)
        f = QFormLayout(w)

        self.cb_dem = QgsMapLayerComboBox()
        self.cb_dem.setFilters(QgsMapLayerProxyModel.RasterLayer)
        self.cb_dem.setAllowEmptyLayer(False)

        self.cb_points = QgsMapLayerComboBox()
        self.cb_points.setFilters(QgsMapLayerProxyModel.PointLayer)

        self.cb_point_id = QgsFieldComboBox()
        self.cb_point_id.setAllowEmptyFieldName(True)

        self.cb_cn = QgsMapLayerComboBox()
        self.cb_cn.setFilters(QgsMapLayerProxyModel.RasterLayer)
        self.cb_cn.setAllowEmptyLayer(True)
        self._preselect_by_name(self.cb_cn, DEFAULT_CONFIG.get("CN_LAYER"))

        # v1.0 — Curve Number can be a raster OR a single constant value
        self.cb_cn_source = QComboBox()
        self.cb_cn_source.addItem("Raster layer", "raster")
        self.cb_cn_source.addItem("Constant value", "value")
        self.sb_cn_value = QDoubleSpinBox()
        self.sb_cn_value.setRange(1.0, 100.0)
        self.sb_cn_value.setDecimals(1)
        self.sb_cn_value.setSingleStep(1.0)
        self.sb_cn_value.setValue(85.0)

        def _cn_source_changed():
            is_val = (self.cb_cn_source.currentData() == "value")
            self.cb_cn.setEnabled(not is_val)
            self.sb_cn_value.setEnabled(is_val)
        self.cb_cn_source.currentIndexChanged.connect(_cn_source_changed)
        _cn_source_changed()

        self.fw_outdir = QgsFileWidget()
        self.fw_outdir.setStorageMode(QgsFileWidget.GetDirectory)
        self.fw_outdir.setFilePath(str(DEFAULT_CONFIG.get("OUTPUT_DIR", "")))

        self._preselect_by_name(self.cb_dem, DEFAULT_CONFIG.get("DEM_LAYER"))
        self._preselect_by_name(self.cb_points, DEFAULT_CONFIG.get("POINTS_LAYER"))

        f.addRow("DEM raster", self.cb_dem)
        f.addRow("Pour-points layer", self.cb_points)
        f.addRow("Pour-point ID field (optional)", self.cb_point_id)
        f.addRow("Curve Number source", self.cb_cn_source)
        f.addRow("Curve Number raster (optional)", self.cb_cn)
        f.addRow("Constant Curve Number value", self.sb_cn_value)
        f.addRow("Output directory", self.fw_outdir)

        hint = QLabel(
            "<i>The DEM must be in a metric (projected) CRS. Pour points "
            "in a different CRS will be reprojected on the fly. "
            "A CN raster is aligned to the DEM grid automatically; alternatively "
            "choose 'Constant value' as the Curve Number source to apply a "
            "single CN to the whole area.</i>"
        )
        hint.setWordWrap(True)
        f.addRow(hint)

        self.tabs.addTab(w, "Inputs")

    # =========================================================================
    # Tab 2 — Delineation
    # =========================================================================
    def _build_tab_delineation(self):
        w = QWidget(self)
        f = QFormLayout(w)

        self.sb_eps = QDoubleSpinBox()
        self.sb_eps.setDecimals(8)
        self.sb_eps.setRange(1e-9, 1.0)
        self.sb_eps.setSingleStep(1e-6)
        self.sb_eps.setValue(float(DEFAULT_CONFIG.get("DEPRESSION_EPSILON", 1e-6)))

        self.sb_thresh = QSpinBox()
        self.sb_thresh.setRange(1, 10_000_000)
        self.sb_thresh.setValue(int(DEFAULT_CONFIG.get("STREAM_THRESHOLD", 1000)))
        self.sb_thresh.setSuffix(" cells")

        self.sb_snap = QSpinBox()
        self.sb_snap.setRange(0, 1000)
        self.sb_snap.setValue(int(DEFAULT_CONFIG.get("SNAP_RADIUS_CELLS", 25)))
        self.sb_snap.setSuffix(" cells")

        self.cb_slope_from = QComboBox()
        self.cb_slope_from.addItem("original (real topography)", "original")
        self.cb_slope_from.addItem("filled (consistent with flow paths)", "filled")
        idx = self.cb_slope_from.findData(
            str(DEFAULT_CONFIG.get("BASIN_SLOPE_FROM", "original")).lower())
        if idx >= 0:
            self.cb_slope_from.setCurrentIndex(idx)

        self.sb_subarea = QDoubleSpinBox()
        self.sb_subarea.setRange(0.0, 100_000.0)
        self.sb_subarea.setDecimals(3)
        self.sb_subarea.setSingleStep(0.1)
        self.sb_subarea.setValue(float(DEFAULT_CONFIG.get("SUB_BASIN_AREA_KM2", 0.0)))
        self.sb_subarea.setSuffix(" km²")
        self.sb_subarea.setSpecialValueText("0  (disabled — one catchment per pour point)")

        f.addRow("Depression epsilon", self.sb_eps)
        f.addRow("Stream threshold", self.sb_thresh)
        f.addRow("Snap radius", self.sb_snap)
        f.addRow("Basin slope source", self.cb_slope_from)
        f.addRow("Sub-basin target area", self.sb_subarea)

        hint = QLabel(
            "<i>Stream threshold is in cells of contributing area (flow "
            "accumulation). Lower → denser stream network. "
            "Sub-basin area &gt; 0 enables HEC-HMS-style subdivision and "
            "produces SubBasins, SubLFPs, Reaches and Junctions layers.</i>")
        hint.setWordWrap(True)
        f.addRow(hint)

        self.tabs.addTab(w, "Delineation")

    # =========================================================================
    # Tab 2b — Special Delineation Points (v23: enforced junctions; v24:
    #          points-only mode + per-point reservoir storage routing)
    # =========================================================================
    def _build_tab_reservoirs(self):
        w = QWidget(self)
        v = QVBoxLayout(w)

        # ---- point layer + name field + attribute copy --------------
        gb_src = QGroupBox("Special points source")
        f = QFormLayout(gb_src)

        self.cb_dams = QgsMapLayerComboBox()
        self.cb_dams.setFilters(QgsMapLayerProxyModel.PointLayer)
        self.cb_dams.setAllowEmptyLayer(True)

        self.cb_dam_name = QgsFieldComboBox()
        self.cb_dam_name.setAllowEmptyFieldName(True)

        self.lst_dam_fields = QListWidget()
        self.lst_dam_fields.setSelectionMode(QListWidget.NoSelection)
        self.lst_dam_fields.setMaximumHeight(110)

        self._preselect_by_name(self.cb_dams,
                                DEFAULT_CONFIG.get("DAMS_LAYER"))

        f.addRow("Special points layer (optional)", self.cb_dams)
        f.addRow("Point name field", self.cb_dam_name)
        f.addRow("Fields to copy to junction", self.lst_dam_fields)

        # ---- subdivision mode ---------------------------------------
        self.chk_sp_only = QCheckBox(
            "Delineate sub-catchments at these points ONLY "
            "(ignore the Sub-basin target area)")
        self.chk_sp_only.setToolTip(
            "Ticked: junctions = outlet + special points; the "
            "area-based subdivision from the Delineation tab is "
            "ignored (Stream threshold and Snap radius still apply).\n"
            "Unticked: the special points are COMBINED with the "
            "Delineation-tab options — area-based junctions AND "
            "point-enforced junctions.")
        f.addRow(self.chk_sp_only)
        v.addWidget(gb_src)

        # ---- per-point reservoir definition (v24) --------------------
        self.gb_res = QGroupBox(
            "Reservoir definition per point — modified-Puls storage "
            "routing")
        rf = QFormLayout(self.gb_res)

        self.cb_sp_point = QComboBox()
        rf.addRow("Special point", self.cb_sp_point)

        self.chk_res = QCheckBox(
            "Define this point as a Reservoir (route its inflow "
            "through storage)")
        rf.addRow(self.chk_res)

        self.cb_res_mode = QComboBox()
        self.cb_res_mode.addItem("Elevation–volume + outflow structures",
                                 "structures")
        self.cb_res_mode.addItem("Storage–discharge table", "table")
        rf.addRow("Definition mode", self.cb_res_mode)

        self.stk_res = QStackedWidget()

        # -- page 0: EV curve + structures --
        pg_st = QWidget()
        stf = QFormLayout(pg_st)
        stf.setContentsMargins(0, 0, 0, 0)
        self.txt_ev = QPlainTextEdit()
        self.txt_ev.setPlaceholderText(
            "elevation (m), volume (10\u00b3 m\u00b3) — one row per "
            "line, e.g.\n410.0, 0\n412.5, 850\n415.0, 1850")
        self.txt_ev.setMaximumHeight(90)
        stf.addRow("Elevation–volume curve", self.txt_ev)

        self.gb_spill = QGroupBox("Spillway (broad-crested weir: "
                                  "Q = C·L·H^1.5)")
        self.gb_spill.setCheckable(True)
        self.gb_spill.setChecked(True)
        spf = QFormLayout(self.gb_spill)
        self.sb_sp_crest = QDoubleSpinBox()
        self.sb_sp_crest.setRange(-1000.0, 10000.0)
        self.sb_sp_crest.setDecimals(2)
        self.sb_sp_C = QDoubleSpinBox()
        self.sb_sp_C.setRange(0.1, 5.0); self.sb_sp_C.setDecimals(2)
        self.sb_sp_C.setValue(1.7)
        self.sb_sp_L = QDoubleSpinBox()
        self.sb_sp_L.setRange(0.1, 10000.0); self.sb_sp_L.setDecimals(2)
        self.sb_sp_L.setValue(5.0)
        spf.addRow("Crest elevation (m)", self.sb_sp_crest)
        spf.addRow("Coefficient C", self.sb_sp_C)
        spf.addRow("Crest length L (m)", self.sb_sp_L)
        stf.addRow(self.gb_spill)

        self.gb_pipe = QGroupBox("Low-level outlet (orifice: "
                                 "Q = n·Cd·A·\u221a(2gH))")
        self.gb_pipe.setCheckable(True)
        self.gb_pipe.setChecked(False)
        ppf = QFormLayout(self.gb_pipe)
        self.sb_pp_inv = QDoubleSpinBox()
        self.sb_pp_inv.setRange(-1000.0, 10000.0)
        self.sb_pp_inv.setDecimals(2)
        self.sb_pp_n = QSpinBox()
        self.sb_pp_n.setRange(1, 50); self.sb_pp_n.setValue(1)
        self.sb_pp_r = QDoubleSpinBox()
        self.sb_pp_r.setRange(0.01, 20.0); self.sb_pp_r.setDecimals(3)
        self.sb_pp_r.setValue(0.45)
        self.sb_pp_cd = QDoubleSpinBox()
        self.sb_pp_cd.setRange(0.05, 1.0); self.sb_pp_cd.setDecimals(2)
        self.sb_pp_cd.setValue(0.6)
        # v1.0 — outlet shape: circular (radius R) or rectangular
        # (span S x height H). HEC-HMS receives area and CENTRE
        # elevation accordingly: circle pi*R^2 / invert+R,
        # rectangle S*H / invert+H/2.
        self.cb_pp_shape = QComboBox()
        self.cb_pp_shape.addItem("Circular (radius R)", "circle")
        self.cb_pp_shape.addItem("Rectangular (span S \u00d7 height H)",
                                 "rect")
        self.sb_pp_span = QDoubleSpinBox()
        self.sb_pp_span.setRange(0.01, 100.0)
        self.sb_pp_span.setDecimals(3)
        self.sb_pp_span.setValue(1.0)
        self.sb_pp_h = QDoubleSpinBox()
        self.sb_pp_h.setRange(0.01, 100.0)
        self.sb_pp_h.setDecimals(3)
        self.sb_pp_h.setValue(1.0)

        ppf.addRow("Invert elevation (m)", self.sb_pp_inv)
        ppf.addRow("Number of barrels n", self.sb_pp_n)
        ppf.addRow("Outlet shape", self.cb_pp_shape)
        self._row_pp_r = ("Radius r (m)", self.sb_pp_r)
        ppf.addRow(*self._row_pp_r)
        ppf.addRow("Span S (m)", self.sb_pp_span)
        ppf.addRow("Height H (m)", self.sb_pp_h)
        ppf.addRow("Discharge coeff. Cd", self.sb_pp_cd)

        def _pp_shape_changed():
            rect = (self.cb_pp_shape.currentData() == "rect")
            for w in (self.sb_pp_r,):
                w.setEnabled(not rect)
            for w in (self.sb_pp_span, self.sb_pp_h):
                w.setEnabled(rect)
            lbl_r = ppf.labelForField(self.sb_pp_r)
            if lbl_r:
                lbl_r.setEnabled(not rect)
            for w in (self.sb_pp_span, self.sb_pp_h):
                lb = ppf.labelForField(w)
                if lb:
                    lb.setEnabled(rect)
        self.cb_pp_shape.currentIndexChanged.connect(
            lambda *_a: _pp_shape_changed())
        _pp_shape_changed()
        stf.addRow(self.gb_pipe)

        # v1.0 — initial condition of the pool
        self.cb_res_init = QComboBox()
        self.cb_res_init.addItem(
            "Outflow equals inflow (pool at curve datum)", "inflow")
        self.cb_res_init.addItem("Initial elevation (m)", "elevation")
        self.sb_res_init_el = QDoubleSpinBox()
        self.sb_res_init_el.setRange(-1000.0, 10000.0)
        self.sb_res_init_el.setDecimals(2)
        self.sb_res_init_el.setEnabled(False)
        self.sb_res_init_el.setToolTip(
            "Water-surface elevation at the start of the storm. The "
            "storage it implies is read from the elevation\u2013volume "
            "curve, so the flood arrives on a partly full pool with "
            "less room to attenuate.\nExported to HEC-HMS as "
            "'Initial Elevation'. Requires an elevation\u2013volume "
            "curve (not available with a storage\u2013discharge "
            "table).")
        self.cb_res_init.currentIndexChanged.connect(
            lambda *_a: self.sb_res_init_el.setEnabled(
                self.cb_res_init.currentData() == "elevation"))
        stf.addRow("Initial condition", self.cb_res_init)
        stf.addRow("Initial elevation (m)", self.sb_res_init_el)

        self.stk_res.addWidget(pg_st)

        # -- page 1: storage–discharge table --
        pg_tb = QWidget()
        tbf = QFormLayout(pg_tb)
        tbf.setContentsMargins(0, 0, 0, 0)
        self.txt_sd = QPlainTextEdit()
        self.txt_sd.setPlaceholderText(
            "storage (10\u00b3 m\u00b3), outflow (m\u00b3/s) — one row "
            "per line, e.g.\n0, 0\n500, 2.5\n1850, 45.0")
        self.txt_sd.setMaximumHeight(140)
        tbf.addRow("Storage–discharge table", self.txt_sd)
        self.stk_res.addWidget(pg_tb)

        rf.addRow(self.stk_res)
        v.addWidget(self.gb_res)

        hint = QLabel(
            "<i>Each special point is snapped to the nearest <b>main "
            "flow path</b> (highest flow accumulation within the "
            "Delineation-tab <b>Snap radius</b>) and enforced as a "
            "junction — a sub-basin boundary is always created there. "
            "The junction is renamed <b>J1-Sub9(PointName)</b> and the "
            "ticked fields are copied onto it as <b>DAM_*</b> columns. "
            "Points defined as <b>Reservoirs</b> route their combined "
            "inflow through modified-Puls (storage-indication) routing: "
            "the Junctions layer gains "
            "QIN/ATT/SMAX/WLMX columns, QP/V report the routed "
            "outflow, the junction CSV becomes a full routing table "
            "(inflow · outflow · storage · stage), and the HEC-HMS "
            "export writes a true <b>Reservoir</b> element with its "
            "Elevation-Storage table in a companion .pdata. Points "
            "with no stream within the snap radius are skipped with a "
            "warning.</i>")
        hint.setWordWrap(True)
        v.addWidget(hint)
        v.addStretch(1)

        # ---- reservoir-definition state store ------------------------
        self._sp_res = {}            # point name -> definition dict
        self._sp_current = None      # name shown in the editor

        self.cb_res_mode.currentIndexChanged.connect(
            self.stk_res.setCurrentIndex)
        self.chk_res.toggled.connect(self._res_editor_enable)
        self.cb_sp_point.currentIndexChanged.connect(
            self._on_sp_point_changed)
        self._res_editor_enable(False)

        self.tabs.addTab(w, "Special Delineation Points")

    # ---- Special Delineation Points helpers (v24) --------------------
    def _res_editor_enable(self, on):
        self.cb_res_mode.setEnabled(on)
        self.stk_res.setEnabled(on)

    def _res_editor_load(self, name):
        """Populate the reservoir editor with the stored definition of
        `name` (or a blank default)."""
        d = self._sp_res.get(name) or {}
        self.chk_res.setChecked(bool(d.get('enabled')))
        mode = d.get('mode', 'structures')
        self.cb_res_mode.setCurrentIndex(1 if mode == 'table' else 0)
        self.txt_ev.setPlainText(d.get('ev_text', ''))
        self.txt_sd.setPlainText(d.get('table_text', ''))
        k = self.cb_res_init.findData(d.get('init_mode', 'inflow'))
        self.cb_res_init.setCurrentIndex(k if k >= 0 else 0)
        self.sb_res_init_el.setValue(float(d.get('init_elev', 0.0) or 0.0))
        self.sb_res_init_el.setEnabled(
            self.cb_res_init.currentData() == "elevation")
        sp = pp = None
        for st in (d.get('structures') or []):
            if st.get('type') == 'spillway' and sp is None:
                sp = st
            elif st.get('type') == 'pipe' and pp is None:
                pp = st
        self.gb_spill.setChecked(sp is not None or not d)
        if sp:
            self.sb_sp_crest.setValue(float(sp.get('crest', 0.0)))
            self.sb_sp_C.setValue(float(sp.get('C', 1.7)))
            self.sb_sp_L.setValue(float(sp.get('width', 5.0)))
        self.gb_pipe.setChecked(pp is not None)
        if pp:
            self.sb_pp_inv.setValue(float(pp.get('invert', 0.0)))
            self.sb_pp_n.setValue(int(pp.get('n', 1)))
            k = self.cb_pp_shape.findData(pp.get('shape', 'circle'))
            self.cb_pp_shape.setCurrentIndex(k if k >= 0 else 0)
            self.sb_pp_r.setValue(float(pp.get('radius', 0.45)))
            self.sb_pp_span.setValue(float(pp.get('span', 1.0)))
            self.sb_pp_h.setValue(float(pp.get('height', 1.0)))
            self.sb_pp_cd.setValue(float(pp.get('Cd', 0.6)))

    def _res_editor_save(self, name):
        """Persist the editor state for `name` into self._sp_res."""
        if not name:
            return
        structures = []
        if self.gb_spill.isChecked():
            structures.append(dict(
                type='spillway',
                crest=float(self.sb_sp_crest.value()),
                C=float(self.sb_sp_C.value()),
                width=float(self.sb_sp_L.value())))
        if self.gb_pipe.isChecked():
            structures.append(dict(
                type='pipe',
                invert=float(self.sb_pp_inv.value()),
                n=int(self.sb_pp_n.value()),
                shape=self.cb_pp_shape.currentData(),
                radius=float(self.sb_pp_r.value()),
                span=float(self.sb_pp_span.value()),
                height=float(self.sb_pp_h.value()),
                Cd=float(self.sb_pp_cd.value())))
        self._sp_res[name] = dict(
            enabled=self.chk_res.isChecked(),
            mode=self.cb_res_mode.currentData(),
            ev_text=self.txt_ev.toPlainText(),
            structures=structures,
            table_text=self.txt_sd.toPlainText(),
            init_mode=self.cb_res_init.currentData(),
            init_elev=float(self.sb_res_init_el.value()))

    def _on_sp_point_changed(self):
        if self._sp_current:
            self._res_editor_save(self._sp_current)
        name = self.cb_sp_point.currentText() or None
        self._sp_current = name
        if name:
            self._res_editor_load(name)

    def _refresh_sp_points(self):
        """Repopulate the special-point selector from the layer + the
        chosen name field."""
        if self._sp_current:
            self._res_editor_save(self._sp_current)
        self.cb_sp_point.blockSignals(True)
        self.cb_sp_point.clear()
        layer = self.cb_dams.currentLayer()
        fld = self.cb_dam_name.currentField()
        names = []
        if isinstance(layer, QgsVectorLayer) and fld \
                and layer.fields().indexOf(fld) >= 0:
            seen = set()
            for i, ft in enumerate(layer.getFeatures()):
                nm = ft[fld]
                nm = (str(nm).strip() if nm is not None else '') \
                    or f"Dam_{i + 1}"
                if nm not in seen:
                    seen.add(nm)
                    names.append(nm)
        self.cb_sp_point.addItems(names)
        self.cb_sp_point.blockSignals(False)
        self.gb_res.setEnabled(bool(names))
        self._sp_current = None
        if names:
            self._on_sp_point_changed()

    # =========================================================================
    # Tab 3 — Time of concentration & losses
    # =========================================================================
    def _build_tab_tc_loss(self):
        w = QWidget(self)
        f = QFormLayout(w)

        self.cb_tc = QComboBox()
        self.cb_tc.addItem("Kirpich — holistic (default)", "KIRPICH_HOLISTIC")
        self.cb_tc.addItem("Kirpich — 10-85 segmented",   "KIRPICH_10_85")
        self.cb_tc.addItem("NRCS lag (requires a CN source)", "NRCS")
        idx = self.cb_tc.findData(
            str(DEFAULT_CONFIG.get("TC_METHOD", "KIRPICH_HOLISTIC")).upper())
        if idx >= 0:
            self.cb_tc.setCurrentIndex(idx)

        # v1.0 — the peak rate factors HEC-HMS itself offers, so the
        # plugin and the exported .basin always describe the same unit
        # hydrograph. 'Custom' still allows any value for AR-Tarseem's
        # own computation; the export then snaps to the nearest
        # HEC-HMS option and says so in the log.
        self.cb_prf = QComboBox()
        self.cb_prf.addItem("Standard (PRF 484)", 484)
        self.cb_prf.addItem("Delmarva (PRF 284)", 284)
        for _p in (100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600):
            self.cb_prf.addItem(f"Peak Rate Factor {_p}", _p)
        _k = self.cb_prf.findData(
            int(DEFAULT_CONFIG.get("PEAK_RATE_FACTOR", 484)))
        self.cb_prf.setCurrentIndex(_k if _k >= 0 else 0)
        self.cb_prf.setToolTip(
            "The peak rate factors HEC-HMS provides. Only these values "
            "are offered, so the plugin and the exported .basin always "
            "describe the same unit hydrograph (Unitgraph Type: "
            "STANDARD, DELMARVA or PRF<n>).\n"
            "Lower values give a flatter, longer unit hydrograph; "
            "higher values a sharper, shorter one. The adopted value "
            "is written to the PRF field of the Catchments and "
            "SubBasins layers.")

        # v19 — initial abstraction ratio λ = Ia/S
        self.sb_lambda = QDoubleSpinBox()
        self.sb_lambda.setRange(0.01, 0.5)
        self.sb_lambda.setDecimals(2)
        self.sb_lambda.setSingleStep(0.01)
        self.sb_lambda.setValue(float(DEFAULT_CONFIG.get("SCS_LAMBDA", 0.2)))

        f.addRow("Tc / Tlag method (drives Q-peak)", self.cb_tc)
        f.addRow("SCS UH peak rate factor (PRF)", self.cb_prf)
        f.addRow("λ = Ia/S  (initial abstraction ratio)", self.sb_lambda)

        # v1.0 — lower bound on the Tc that feeds the hydrograph engine
        self.sb_tc_min = QDoubleSpinBox()
        self.sb_tc_min.setRange(0.0, 100000.0)
        self.sb_tc_min.setDecimals(2)
        self.sb_tc_min.setSingleStep(5.0)
        self.sb_tc_min.setSuffix(" min")
        self.sb_tc_min.setSpecialValueText("No minimum")
        self.sb_tc_min.setValue(
            float(DEFAULT_CONFIG.get("TC_MIN_MINUTES", 0.0) or 0.0))
        self.sb_tc_min.setToolTip(
            "Lower bound on the time of concentration used by the "
            "peak-flow / hydrograph engine (0 = no minimum).\n"
            "The floor is applied to Tc and the matching "
            "Tlag = 0.6 \u00d7 Tc is used with it. The three computed "
            "Tc estimates are still written to the attribute tables "
            "unchanged.")
        f.addRow("Minimum Tc  (Tlag = 0.6 × Tc)", self.sb_tc_min)

        hint = QLabel(
            "<i>All three Tc estimates are always written to attributes; "
            "the selection only decides which one feeds the peak-flow / "
            "hydrograph engine. PRF: 100-200 = swampy/flat, 484 = "
            "standard, 600 = steep/urban.<br><br>"
            "λ = Ia/S: 0.2 is the "
            "classic SCS assumption; 0.05 is recommended by newer NRCS "
            "guidance for many watersheds. Runoff uses the generalised "
            "form R = (P−λS)² / (P−λS+S), which reduces to the familiar "
            "(P−0.2S)²/(P+0.8S) at λ = 0.2. Applies to attribute-table "
            "R/V columns AND the hydrograph engine.</i>")
        hint.setWordWrap(True)
        f.addRow(hint)

        self.tabs.addTab(w, "Tc && losses")

    # =========================================================================
    # Tab 4 — Rainfall & hydrographs
    # =========================================================================
    def _build_tab_rainfall_hydrographs(self):
        w = QWidget(self)
        v = QVBoxLayout(w)

        gb_rain = QGroupBox("Design rainfall (return periods fully editable)")
        rain_v = QVBoxLayout(gb_rain)

        # v22: BOTH columns are editable and rows can be added/removed, so
        # any custom return period (e.g. 3, 15, 500 yrs) can be defined.
        # The third column holds the station attribute field used in the
        # spatially-variable (Thiessen) mode.
        self.tbl_rainfall = QTableWidget(len(_DESIGN_RETURN_PERIODS), 3)
        self.tbl_rainfall.setHorizontalHeaderLabels(
            ["Return period (yrs)", "Depth (mm)", "Station field (Thiessen)"])
        self.tbl_rainfall.verticalHeader().setVisible(False)
        existing = DEFAULT_CONFIG.get("RAINFALL_DEPTHS_MM") or {}
        for row, rp in enumerate(_DESIGN_RETURN_PERIODS):
            self.tbl_rainfall.setItem(row, 0, QTableWidgetItem(str(rp)))
            depth = existing.get(rp)
            txt = "" if depth in (None, 0, 0.0) else str(depth)
            self.tbl_rainfall.setItem(row, 1, QTableWidgetItem(txt))
            self.tbl_rainfall.setItem(row, 2, QTableWidgetItem(""))
        self.tbl_rainfall.resizeColumnsToContents()
        rain_v.addWidget(self.tbl_rainfall)

        btn_row = QHBoxLayout()
        self.btn_rp_add = QPushButton("+ Add return period")
        self.btn_rp_del = QPushButton("\u2212 Remove selected")
        btn_row.addWidget(self.btn_rp_add)
        btn_row.addWidget(self.btn_rp_del)
        btn_row.addStretch(1)
        rain_v.addLayout(btn_row)

        def _add_rp_row():
            r = self.tbl_rainfall.rowCount()
            self.tbl_rainfall.insertRow(r)
            for c_ in range(3):
                self.tbl_rainfall.setItem(r, c_, QTableWidgetItem(""))
            self.tbl_rainfall.setCurrentCell(r, 0)

        def _del_rp_rows():
            rows = sorted({i.row() for i in
                           self.tbl_rainfall.selectedIndexes()},
                          reverse=True)
            for r in rows:
                self.tbl_rainfall.removeRow(r)

        self.btn_rp_add.clicked.connect(_add_rp_row)
        self.btn_rp_del.clicked.connect(_del_rp_rows)

        rain_v.addWidget(QLabel(
            "<i>Any positive integer return period is accepted. Leave a "
            "depth blank (or 0) to skip that row. R and V columns require "
            "the CN source (raster or constant value).</i>"))
        v.addWidget(gb_rain)

        # ---- v22: spatially-variable rainfall (Thiessen) -----------------
        self.gb_thiessen = QGroupBox(
            "Spatially-variable rainfall \u2014 Thiessen polygons of "
            "rainfall stations")
        self.gb_thiessen.setCheckable(True)
        self.gb_thiessen.setChecked(
            bool(str(DEFAULT_CONFIG.get("RAIN_STATIONS_LAYER", "")).strip()))
        th_f = QFormLayout(self.gb_thiessen)

        self.cb_stations = QgsMapLayerComboBox()
        self.cb_stations.setFilters(QgsMapLayerProxyModel.PointLayer)
        self.cb_stations.setAllowEmptyLayer(True)
        self._preselect_by_name(self.cb_stations,
                                DEFAULT_CONFIG.get("RAIN_STATIONS_LAYER"))
        th_f.addRow("Rainfall stations layer", self.cb_stations)
        th_f.addRow(QLabel(
            "<i>Point layer with one depth field per return "
            "period. Thiessen (Voronoi) weighting gives every "
            "catchment and sub-basin its own areal depth, "
            "overriding the design-rainfall table. The "
            "polygons used for that weighting are added to "
            "the project as the <b>Thiessen</b> layer, each "
            "carrying its station name and depths.</i>"))
        v.addWidget(self.gb_thiessen)

        gb_hyd = QGroupBox("Design-storm hydrograph engine")
        hyd_f = QFormLayout(gb_hyd)

        self.chk_hyd = QCheckBox("Compute hydrograph CSVs + interactive HTMLs")
        self.chk_hyd.setChecked(bool(DEFAULT_CONFIG.get("HYDROGRAPHS_ENABLED", True)))
        hyd_f.addRow(self.chk_hyd)

        self.sb_storm_hr = QDoubleSpinBox()
        self.sb_storm_hr.setRange(1.0, 240.0)
        self.sb_storm_hr.setDecimals(1)
        self.sb_storm_hr.setSingleStep(1.0)
        self.sb_storm_hr.setSuffix(" h")
        self.sb_storm_hr.setValue(float(DEFAULT_CONFIG.get("STORM_DURATION_HR", 24.0)))
        hyd_f.addRow("Storm duration", self.sb_storm_hr)

        self.sb_dt_min = QDoubleSpinBox()
        self.sb_dt_min.setRange(0.1, 60.0)
        self.sb_dt_min.setDecimals(2)
        self.sb_dt_min.setSingleStep(0.1)
        self.sb_dt_min.setSuffix(" min")
        self.sb_dt_min.setValue(float(DEFAULT_CONFIG.get("HYD_TIME_STEP_MIN", 1.0)))
        hyd_f.addRow("Time step", self.sb_dt_min)

        # v1.0 — editable output window so the CSV/HTML hydrographs are
        # complete (flow recedes to ~zero before the window ends)
        self.sb_tail_factor = QDoubleSpinBox()
        self.sb_tail_factor.setRange(5.0, 500.0)
        self.sb_tail_factor.setDecimals(0)
        self.sb_tail_factor.setSingleStep(5.0)
        self.sb_tail_factor.setSuffix(" × Tlag")
        self.sb_tail_factor.setValue(
            float(DEFAULT_CONFIG.get("HYD_TAIL_FACTOR", 20.0)))
        self.sb_tail_factor.setToolTip(
            "Post-storm window length as a multiple of the lag time "
            "(PRF-adjusted). Increase if the console warns that a "
            "hydrograph has not fully receded to zero.")
        hyd_f.addRow("Hydrograph tail length", self.sb_tail_factor)

        self.sb_tail_cap = QDoubleSpinBox()
        self.sb_tail_cap.setRange(1.0, 2000.0)
        self.sb_tail_cap.setDecimals(0)
        self.sb_tail_cap.setSingleStep(24.0)
        self.sb_tail_cap.setSuffix(" h")
        self.sb_tail_cap.setValue(
            float(DEFAULT_CONFIG.get("HYD_MAX_TAIL_HR", 96.0)))
        self.sb_tail_cap.setToolTip(
            "Upper cap on the post-storm window. Raise it together "
            "with the tail length for very slow basins so the output "
            "hydrograph is complete (ends at zero flow).")
        hyd_f.addRow("Max tail after storm", self.sb_tail_cap)

        # v19 — all four SCS/NRCS 24-h distributions
        self.cb_dist = QComboBox()
        self.cb_dist.addItem("SCS Type I (24-hour) — Pacific maritime", "SCS_TYPE_I")
        self.cb_dist.addItem("SCS Type IA (24-hour) — low-intensity coastal", "SCS_TYPE_IA")
        self.cb_dist.addItem("SCS Type II (24-hour) — standard/default", "SCS_TYPE_II")
        self.cb_dist.addItem("SCS Type III (24-hour) — Gulf/Atlantic tropical", "SCS_TYPE_III")
        self.cb_dist.addItem("Custom hyetograph (CSV)", "CUSTOM")
        idx = self.cb_dist.findData(
            str(DEFAULT_CONFIG.get("STORM_DISTRIBUTION", "SCS_TYPE_II")).upper())
        if idx >= 0:
            self.cb_dist.setCurrentIndex(idx)
        hyd_f.addRow("Storm distribution", self.cb_dist)

        # v19.1 — custom hyetograph CSV
        self.fw_custom_csv = QgsFileWidget()
        self.fw_custom_csv.setStorageMode(QgsFileWidget.GetFile)
        self.fw_custom_csv.setFilter("CSV files (*.csv *.txt);;All files (*)")
        self.fw_custom_csv.setFilePath(
            str(DEFAULT_CONFIG.get("CUSTOM_HYETOGRAPH_CSV", "")))
        hyd_f.addRow("Custom hyetograph CSV", self.fw_custom_csv)

        v.addWidget(gb_hyd)
        v.addWidget(QLabel(
            "<i>The four SCS distributions use exact 1-minute "
            "cumulative tables. <b>Custom "
            "hyetograph</b>: CSV with two columns — time (min), depth "
            "(mm); cumulative or incremental depths are auto-detected; "
            "the curve is normalised and rescaled to each return "
            "period's P, and its last time defines the storm duration "
            "(the duration spinbox is ignored). Outputs land in "
            "<b>&lt;output dir&gt;/Modeling_Outputs/</b>.</i>"))

        self.tabs.addTab(w, "Rainfall && hydrographs")

    # =========================================================================
    # Tab 5 — Reach routing & ARF
    # =========================================================================
    def _build_tab_routing_arf(self):
        w = QWidget(self)
        v = QVBoxLayout(w)

        gb_route = QGroupBox("Reach routing of hydrographs (sub-basins only)")
        route_f = QFormLayout(gb_route)

        self.cb_route = QComboBox()
        self.cb_route.addItem("Muskingum (K from reach Tc, X below)", "MUSKINGUM")
        self.cb_route.addItem("Muskingum-Cunge (channel geometry)", "MUSKINGUM_CUNGE")
        self.cb_route.addItem("Lag (pure translation)", "LAG")
        self.cb_route.addItem("None (outflow == inflow)", "NONE")
        idx = self.cb_route.findData(
            str(DEFAULT_CONFIG.get("REACH_ROUTING_METHOD", "MUSKINGUM")).upper())
        if idx >= 0:
            self.cb_route.setCurrentIndex(idx)

        self.sb_musk_x = QDoubleSpinBox()
        self.sb_musk_x.setRange(0.0, 0.5)
        self.sb_musk_x.setDecimals(3)
        self.sb_musk_x.setSingleStep(0.05)
        self.sb_musk_x.setValue(float(DEFAULT_CONFIG.get("MUSKINGUM_X", 0.2)))

        # v1.0 — Muskingum subreach count: Automatic (per-reach,
        # stability-based, as before) or a fixed value for all reaches
        self.sb_musk_nsub = QSpinBox()
        self.sb_musk_nsub.setRange(0, 5000)
        self.sb_musk_nsub.setSpecialValueText("Automatic (per reach)")
        self.sb_musk_nsub.setValue(
            int(DEFAULT_CONFIG.get("MUSK_SUBREACHES", 0) or 0))
        self.sb_musk_nsub.setToolTip(
            "Number of Muskingum subreaches applied to EVERY reach "
            "(also written to the HEC-HMS .basin 'Subreaches:' line).\n"
            "Automatic = per-reach N = ceil(2·K·X/Δt), the "
            "stability-based estimate.\n"
            "Physically: increasing the subreaches reduces the "
            "attenuation and travel time per segment (the cascade "
            "approaches pure translation); decreasing them increases "
            "the attenuation and travel time per segment (approaches "
            "a single linear reservoir).")

        self.cb_reach_tc = QComboBox()
        self.cb_reach_tc.addItem("Kirpich — holistic (default)", "KIRPICH_HOLISTIC")
        self.cb_reach_tc.addItem("Kirpich — 10-85 segmented",   "KIRPICH_10_85")
        idx = self.cb_reach_tc.findData(
            str(DEFAULT_CONFIG.get("REACH_TC_METHOD", "KIRPICH_HOLISTIC")).upper())
        if idx >= 0:
            self.cb_reach_tc.setCurrentIndex(idx)

        route_f.addRow("Routing method", self.cb_route)
        route_f.addRow("Muskingum X", self.sb_musk_x)
        route_f.addRow("Muskingum subreaches", self.sb_musk_nsub)

        # v1.0 — lower bound on the reach travel time (drives K / lag)
        self.sb_reach_tc_min = QDoubleSpinBox()
        self.sb_reach_tc_min.setRange(0.0, 100000.0)
        self.sb_reach_tc_min.setDecimals(2)
        self.sb_reach_tc_min.setSingleStep(1.0)
        self.sb_reach_tc_min.setSuffix(" min")
        self.sb_reach_tc_min.setSpecialValueText("No minimum")
        self.sb_reach_tc_min.setValue(
            float(DEFAULT_CONFIG.get("REACH_TC_MIN_MINUTES", 0.0) or 0.0))
        self.sb_reach_tc_min.setToolTip(
            "Lower bound on the reach travel time that drives "
            "Muskingum K (K = reach Tc) and the Lag shift "
            "(0 = no minimum).\nVery short headwater reaches "
            "otherwise route almost instantaneously.")
        route_f.addRow("Minimum reach travel time", self.sb_reach_tc_min)
        route_f.addRow("Reach Tc source (drives K / lag)", self.cb_reach_tc)
        v.addWidget(gb_route)

        # v19 — Muskingum-Cunge trapezoidal channel geometry
        self.gb_mc = QGroupBox("Muskingum-Cunge channel geometry (trapezoid)")
        mc_f = QFormLayout(self.gb_mc)

        self.sb_mc_n = QDoubleSpinBox()
        self.sb_mc_n.setRange(0.008, 0.2)
        self.sb_mc_n.setDecimals(3)
        self.sb_mc_n.setSingleStep(0.005)
        self.sb_mc_n.setValue(float(DEFAULT_CONFIG.get("MC_MANNING_N", 0.035)))

        self.sb_mc_b = QDoubleSpinBox()
        self.sb_mc_b.setRange(0.05, 500.0)
        self.sb_mc_b.setDecimals(2)
        self.sb_mc_b.setSingleStep(0.5)
        self.sb_mc_b.setSuffix(" m")
        self.sb_mc_b.setValue(float(DEFAULT_CONFIG.get("MC_BOTTOM_WIDTH_M", 5.0)))

        self.sb_mc_z = QDoubleSpinBox()
        self.sb_mc_z.setRange(0.0, 10.0)
        self.sb_mc_z.setDecimals(2)
        self.sb_mc_z.setSingleStep(0.5)
        self.sb_mc_z.setValue(float(DEFAULT_CONFIG.get("MC_SIDE_SLOPE", 2.0)))

        mc_f.addRow("Manning n", self.sb_mc_n)
        mc_f.addRow("Bottom width b", self.sb_mc_b)
        mc_f.addRow("Side slope z (H:V)", self.sb_mc_z)
        mc_f.addRow(QLabel(
            "<i>Variable-parameter M-C: K and X "
            "recomputed each step from Manning hydraulics; the reach is "
            "subdivided into N = L/(c·Δt) segments (≤ 40). Reach length "
            "and slope come from the delineated reach geometry; this "
            "section shape applies to all routed reaches.</i>"))
        v.addWidget(self.gb_mc)

        gb_arf = QGroupBox("Areal Reduction Factor (ARF)")
        arf_f = QFormLayout(gb_arf)

        self.cb_arf = QComboBox()
        self.cb_arf.addItem("None — use point depths as-is (default)", "NONE")
        self.cb_arf.addItem("TRRL_UK: 1 - 0.044·A^0.275", "TRRL_UK")
        self.cb_arf.addItem("Riyadh: ARF=1 if A≤133 km², else 1.4167-0.0852·ln(A)", "RIYADH")
        self.cb_arf.addItem("Fixed — user-defined value (all catchments)", "FIXED")
        idx = self.cb_arf.findData(
            str(DEFAULT_CONFIG.get("ARF_METHOD", "NONE")).upper())
        if idx >= 0:
            self.cb_arf.setCurrentIndex(idx)
        arf_f.addRow("ARF method", self.cb_arf)

        # v20 — direct user-defined ARF value
        self.sb_arf_fixed = QDoubleSpinBox()
        self.sb_arf_fixed.setRange(0.01, 1.0)
        self.sb_arf_fixed.setDecimals(3)
        self.sb_arf_fixed.setSingleStep(0.01)
        self.sb_arf_fixed.setValue(
            float(DEFAULT_CONFIG.get("ARF_FIXED_VALUE", 1.0)))
        arf_f.addRow("Fixed ARF value", self.sb_arf_fixed)
        arf_f.addRow(QLabel(
            "<i>With <b>Fixed</b>, this exact value multiplies the "
            "design rainfall of every catchment and its sub-basins, "
            "regardless of area. With the equation methods the value "
            "is computed per catchment from its full upstream "
            "area.</i>"))
        v.addWidget(gb_arf)
        v.addStretch(1)

        self.tabs.addTab(w, "Routing && ARF")

    # =========================================================================
    # Tab 6 — HEC-HMS export (v19)
    # =========================================================================
    # =========================================================================
    # Tab — AR-Awdiya (v1.0): export the whole model as a ready-to-run
    #        AR-Awdiya project (.json)
    # =========================================================================
    def _build_tab_awdiya(self):
        w = QWidget(self)
        f = QFormLayout(w)

        self.chk_awdiya = QCheckBox(
            "Export the model as an AR-Awdiya project (.json)")
        self.chk_awdiya.setChecked(
            bool(DEFAULT_CONFIG.get("EXPORT_AWDIYA", False)))
        f.addRow(self.chk_awdiya)

        self.ed_awdiya_name = QLineEdit()
        self.ed_awdiya_name.setPlaceholderText(
            "Project name shown in AR-Awdiya (optional)")
        f.addRow("Project name", self.ed_awdiya_name)

        self.cb_awdiya_freq = QComboBox()
        self.cb_awdiya_freq.addItem(
            "All return periods (frequency analysis)", "all")
        self.cb_awdiya_freq.addItem(
            "Single storm — keep each catchment's own depth", "single")
        self.cb_awdiya_freq.setToolTip(
            "All return periods: the depth table is exported and "
            "AR-Awdiya runs every return period, applying each depth "
            "to all catchments.\n"
            "Single storm: one event in which every catchment keeps its "
            "OWN depth, so Thiessen rainfall variation and the areal "
            "reduction stay exactly as computed here.")
        f.addRow("Storm scenarios", self.cb_awdiya_freq)

        self.cb_awdiya_rp = QComboBox()
        self.cb_awdiya_rp.addItem("Largest available", 0)
        f.addRow("Per-catchment depth from", self.cb_awdiya_rp)

        f.addRow(QLabel(
            "<i>Writes <b>&lt;output dir&gt;/AR_Awdiya/"
            "AR_Tarseem_project.json</b>. Open it in AR-Awdiya with "
            "<b>Load project</b> — the network is ready to Run: every "
            "sub-basin becomes a catchment (area, Curve Number, lag, "
            "peak rate factor, \u03bb, storm type and depth), every "
            "junction carries its reach routing (Lag / Muskingum K,X / "
            "Muskingum-Cunge geometry), reservoirs keep their "
            "elevation\u2013volume curves and outflow structures, and "
            "each pour point becomes an outlet. Storm depths are exported "
            "as the <b>unreduced point depths</b> with the areal "
            "reduction factor placed in each catchment's <b>ARF</b> "
            "field, so AR-Awdiya reproduces the same effective depth "
            "while both values stay visible and editable (the HEC-HMS "
            "export is unaffected). In AR-Awdiya you can "
            "then edit the network interactively and generate Word and "
            "PDF report documents of all outputs.</i>"))

        self.tabs.addTab(w, "AR-Awdiya")

    def _build_tab_hms(self):
        w = QWidget(self)
        f = QFormLayout(w)

        self.chk_hms = QCheckBox("Export the delineated network as a "
                                 "HEC-HMS .basin file")
        self.chk_hms.setChecked(bool(DEFAULT_CONFIG.get("EXPORT_HEC_HMS", False)))
        f.addRow(self.chk_hms)

        self.ed_hms_name = QLineEdit()
        self.ed_hms_name.setPlaceholderText("AR_Tarseem")
        self.ed_hms_name.setToolTip(
            "Base name for the exported HEC-HMS files: "
            "<name>.basin, <name>.pdata and the .met files, and the "
            "basin/project name written inside them. Leave empty to "
            "use 'AR_Tarseem'.")
        f.addRow("Output file name", self.ed_hms_name)

        self.cb_hms_ver = QComboBox()
        for v_ in ("4.7", "4.8", "4.9", "4.10", "4.11", "4.12", "4.13"):
            self.cb_hms_ver.addItem(f"HEC-HMS v{v_}", v_)
        idx = self.cb_hms_ver.findData(
            str(DEFAULT_CONFIG.get("HEC_HMS_VERSION", "4.12")))
        if idx >= 0:
            self.cb_hms_ver.setCurrentIndex(idx)
        f.addRow("Target HEC-HMS version", self.cb_hms_ver)

        self.chk_hms_met = QCheckBox(
            "Also export .met files (one per return period, "
            "per-subbasin storm depths)")
        self.chk_hms_met.setChecked(
            bool(DEFAULT_CONFIG.get("EXPORT_HEC_HMS_MET", False)))
        f.addRow(self.chk_hms_met)

        f.addRow(QLabel(
            "<i>Writes "
            "<b>&lt;output dir&gt;/HEC_HMS/AR_Tarseem.basin</b> "
            "containing:<br>"
            "&nbsp;&nbsp;• one <b>Subbasin</b> per sub-basin (Area, SCS "
            "Curve Number, SCS transform with Lag in minutes),<br>"
            "&nbsp;&nbsp;• one <b>Reach</b> per routed link — Lag / "
            "Muskingum / Muskingum-Cunge, matching the routing tab,<br>"
            "&nbsp;&nbsp;• <b>Junctions</b> at sub-basin outlets and one "
            "<b>Sink</b> per pour point,<br>"
            "&nbsp;&nbsp;• real map coordinates as canvas positions, so "
            "the network lays out geographically.<br><br>"
            "Catchments without subdivision export as a single Subbasin "
            "→ Sink pair. Version matters: 4.11+ uses the newer header "
            "block ('Unregulated Output Required'); older versions get "
            "'Enable Quality Routing' — mirroring HEC-HMS itself. Open "
            "via <b>File → Import → Basin Model</b> in HEC-HMS, or in "
            "WMS.</i>"))

        self.tabs.addTab(w, "HEC-HMS")

    # =========================================================================
    # Tab 7 — Styling
    # =========================================================================
    def _build_tab_styling(self):
        w = QWidget(self)
        f = QFormLayout(w)

        self.fw_styles_dir = QgsFileWidget()
        self.fw_styles_dir.setStorageMode(QgsFileWidget.GetDirectory)
        self.fw_styles_dir.setFilePath(str(DEFAULT_CONFIG.get("STYLES_DIR", "")))
        f.addRow("QML styles directory", self.fw_styles_dir)

        f.addRow(QLabel(
            "<i>If set, QML files named Catchments.qml, "
            "Longest_Flow_Paths.qml, Streams.qml, Pour_Points.qml, "
            "Subbasins.qml, Longest_Flow_Paths_-_Subbasins.qml, "
            "Reaches.qml, Junctions.qml and DEM.qml are applied to the "
            "corresponding output layers. Missing files are warned about "
            "but never abort the run.</i>"))

        self.tabs.addTab(w, "Styling")

    # =========================================================================
    # Tab 8 — Output names
    # =========================================================================
    def _build_tab_outputs(self):
        w = QWidget(self)
        f = QFormLayout(w)

        self.le_out_streams = QLineEdit(DEFAULT_CONFIG.get("OUT_STREAMS", "Streams"))
        self.le_out_catch   = QLineEdit(DEFAULT_CONFIG.get("OUT_CATCHMENTS", "Catchments"))
        self.le_out_lfp_pl  = QLineEdit(DEFAULT_CONFIG.get("OUT_LFP_LINES", "LFP_PL"))
        self.le_out_lfp_pt  = QLineEdit(DEFAULT_CONFIG.get("OUT_LFP_POINTS", "LFP_PTs"))
        self.le_out_snap    = QLineEdit(DEFAULT_CONFIG.get("OUT_SNAPPED", "Snapped"))
        self.le_out_thies   = QLineEdit(
            DEFAULT_CONFIG.get("OUT_THIESSEN",
                               "Rainfall_Thiessen_Polygons"))
        self.le_out_subs    = QLineEdit(DEFAULT_CONFIG.get("OUT_SUBBASINS", "SubBasins"))
        self.le_out_sublfp  = QLineEdit(DEFAULT_CONFIG.get("OUT_SUBLFPS", "SubLFPs"))
        self.le_out_reach   = QLineEdit(DEFAULT_CONFIG.get("OUT_REACHES", "Reaches"))
        self.le_out_junc    = QLineEdit(DEFAULT_CONFIG.get("OUT_JUNCTIONS", "Junctions"))

        f.addRow("Streams layer name",    self.le_out_streams)
        f.addRow("Catchments layer name", self.le_out_catch)
        f.addRow("LFP polyline name",     self.le_out_lfp_pl)
        f.addRow("LFP points name",       self.le_out_lfp_pt)
        f.addRow("Snapped points name",   self.le_out_snap)
        f.addRow("Thiessen polygons name", self.le_out_thies)
        f.addRow("SubBasins name",        self.le_out_subs)
        f.addRow("SubLFPs name",          self.le_out_sublfp)
        f.addRow("Reaches name",          self.le_out_reach)
        f.addRow("Junctions name",        self.le_out_junc)

        hint = QLabel("<i>Advanced — only change on name collisions.</i>")
        hint.setWordWrap(True)
        f.addRow(hint)

        self.tabs.addTab(w, "Output names")

    # =========================================================================
    # Tab 9 — About
    # =========================================================================
    def _build_tab_about(self):
        w = QWidget(self)
        v = QVBoxLayout(w)
        v.setContentsMargins(20, 20, 20, 20)

        title = QLabel(f"<h2>AR-Tarseem v{self.PLUGIN_VERSION}</h2>")
        title.setTextFormat(Qt.RichText)
        v.addWidget(title)

        body = QLabel(
            f"<p><b>Copyright (C) 2026 {self.AUTHOR}</b><br>"
            "Licensed under the GNU General Public License v3.0 "
            "or later (GPL-3.0-or-later). This program comes with "
            "ABSOLUTELY NO WARRANTY; see the LICENSE file for "
            "details.</p>"
            "<p>AR-Tarseem is a complete DEM-to-hydrograph watershed "
            "modelling tool for QGIS, built purely on NumPy and GDAL.</p>"
            "<p><b>Capabilities</b></p>"
            "<ul>"
            "<li><b>Watershed delineation</b> — pit filling, flow "
            "direction/accumulation, streams, one catchment per pour "
            "point, longest flow paths, optional sub-basin subdivision "
            "with reaches and junctions.</li>"
            "<li><b>Special delineation points &amp; reservoirs</b> — "
            "enforce junctions at dam locations (snapped to the main "
            "flow path), divide sub-catchments at the points only or "
            "combined with the area-based subdivision, and define any "
            "point as a reservoir with modified-Puls storage routing "
            "(elevation\u2013volume curve + spillway/orifice "
            "structures, or a storage\u2013discharge table) — with "
            "inflow, outflow, attenuation, storage and stage "
            "outputs.</li>"
            "<li><b>Terrain analytics</b> — basin slope, channel segment "
            "slopes, and three time-of-concentration estimates (Kirpich "
            "holistic, Kirpich 10-85, NRCS lag).</li>"
            "<li><b>Design storms</b> — SCS Types I / IA / II / III at "
            "1-minute resolution, or a custom hyetograph CSV; any "
            "custom return periods; areal reduction by TRRL-UK / "
            "Riyadh equations or a fixed user value.</li>"
            "<li><b>Spatially-variable rainfall</b> — rainfall-"
            "station point layers with per-return-period depth "
            "fields; Thiessen/Voronoi weighting gives every "
            "catchment and sub-basin its own depth.</li>"
            "<li><b>Rainfall-runoff</b> — SCS-CN losses with adjustable "
            "\u03bb = Ia/S and the Curve Number taken from a raster or "
            "a single constant value; SCS Unit Hydrograph with "
            "HEC-HMS-matching dimensionless curve and adjustable peak "
            "rate factor.</li>"
            "<li><b>Routing</b> — Lag, Muskingum, and variable-parameter "
            "Muskingum-Cunge on a trapezoidal channel.</li>"
            "<li><b>Results</b> — attribute tables with peaks and "
            "volumes per return period; hydrograph CSVs and interactive "
            "HTML charts for catchments, sub-basins, reaches and "
            "junctions, with an editable output window (tail length "
            "and cap) so every hydrograph is complete and recedes to "
            "zero flow; automatic QML styling.</li>"
            "<li><b>Interoperability</b> — HEC-HMS .basin export of the "
            "whole delineated network (versions 4.7\u20134.13) "
            "including true Reservoir elements, plus .met files per "
            "return period carrying each subbasin's own storm depth "
            "and a .pdata file with the reservoir Elevation-Storage "
            "tables.</li>"
            "</ul>"
            "<p><b>Contact</b><br>"
            "Amr R. Refaiy, PhD — "
            "<a href='mailto:amrgaballa93@gmail.com'>amrgaballa93@gmail.com</a><br>"
            "Ahmed M. Raslan, PhD — "
            "<a href='mailto:raslan.ahmed2033@gmail.com'>raslan.ahmed2033@gmail.com</a></p>"
            "<p><b>Recommended companion — AR-Awdiya</b><br>"
            "For tighter integration, download <b>AR-Awdiya</b>, the "
            "hydrologic network modelling companion of AR-Tarseem. "
            "Open the exported .basin / .pdata files there to edit and "
            "re-run the network interactively, and to generate "
            "polished <b>Word and PDF documents</b> of all the outputs "
            "(hydrographs, reservoir routing tables, summaries) ready "
            "for engineering deliverables.</p>"
            "<p><i>Runs execute as a background task \u2014 QGIS stays responsive and a live progress bar reports each stage. Progress log: View \u2192 Panels \u2192 Log Messages \u2192 \"AR-Tarseem\". Outputs are written to &lt;chosen folder&gt;/ar_tarseem_output/.</i></p>")
        body.setWordWrap(True)
        body.setTextFormat(Qt.RichText)
        body.setOpenExternalLinks(True)   # clickable mailto links
        v.addWidget(body)
        v.addStretch(1)

        self.tabs.addTab(w, "About")

    # =========================================================================
    # Wiring
    # =========================================================================
    def _wire_dependencies(self):
        def _refresh_field_combo():
            layer = self.cb_points.currentLayer()
            if isinstance(layer, QgsVectorLayer):
                self.cb_point_id.setLayer(layer)
            else:
                self.cb_point_id.setLayer(None)

        self.cb_points.layerChanged.connect(_refresh_field_combo)
        _refresh_field_combo()

        # v23 — Reservoirs tab: dam layer drives the name-field combo
        # and the checkable attribute-field list (all ticked by default).
        def _refresh_dam_widgets():
            layer = self.cb_dams.currentLayer()
            self.lst_dam_fields.clear()
            if isinstance(layer, QgsVectorLayer):
                self.cb_dam_name.setLayer(layer)
                for fld in layer.fields():
                    it = QListWidgetItem(fld.name())
                    it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
                    it.setCheckState(Qt.Checked)
                    self.lst_dam_fields.addItem(it)
            else:
                self.cb_dam_name.setLayer(None)
            self._refresh_sp_points()   # v24

        # v1.0 — AR-Awdiya: keep the return-period selector in sync
        # with the design-rainfall table (rows can be added/removed)
        def _refresh_awdiya_rps():
            cur = self.cb_awdiya_rp.currentData()
            self.cb_awdiya_rp.blockSignals(True)
            self.cb_awdiya_rp.clear()
            self.cb_awdiya_rp.addItem("Largest available", 0)
            for row in range(self.tbl_rainfall.rowCount()):
                it = self.tbl_rainfall.item(row, 0)
                txt = (it.text().strip() if it else "")
                try:
                    rp = int(float(txt))
                except (TypeError, ValueError):
                    continue
                self.cb_awdiya_rp.addItem(f"{rp}-year storm", rp)
            k = self.cb_awdiya_rp.findData(cur)
            self.cb_awdiya_rp.setCurrentIndex(k if k >= 0 else 0)
            self.cb_awdiya_rp.blockSignals(False)

        self._refresh_awdiya_rps = _refresh_awdiya_rps
        self.tbl_rainfall.itemChanged.connect(
            lambda *_a: _refresh_awdiya_rps())
        _refresh_awdiya_rps()

        self.cb_dams.layerChanged.connect(_refresh_dam_widgets)
        # v24 — a different name field yields different point names
        self.cb_dam_name.fieldChanged.connect(
            lambda _f: self._refresh_sp_points())
        _refresh_dam_widgets()

        # fixed-ARF spinbox only active when FIXED is selected
        def _arf_toggle():
            self.sb_arf_fixed.setEnabled(
                self.cb_arf.currentData() == "FIXED")
        self.cb_arf.currentIndexChanged.connect(_arf_toggle)
        _arf_toggle()

        # custom-CSV picker only active when CUSTOM is selected; the
        # storm-duration spinbox is ignored for CUSTOM (duration comes
        # from the CSV), so grey it out too
        def _dist_toggle():
            is_custom = self.cb_dist.currentData() == "CUSTOM"
            self.fw_custom_csv.setEnabled(is_custom)
            self.sb_storm_hr.setEnabled(not is_custom)
        self.cb_dist.currentIndexChanged.connect(_dist_toggle)
        _dist_toggle()

        # grey the M-C geometry group unless Muskingum-Cunge is selected
        def _mc_toggle():
            self.gb_mc.setEnabled(
                self.cb_route.currentData() == "MUSKINGUM_CUNGE")
        self.cb_route.currentIndexChanged.connect(_mc_toggle)
        _mc_toggle()

    @staticmethod
    def _preselect_by_name(combo: QgsMapLayerComboBox, name):
        # v1.0 — optional layer combos must start EMPTY. QgsMapLayerComboBox
        # otherwise selects the first matching layer in the project, which
        # silently feeds an unintended layer into the run (a stray raster
        # read as Curve Numbers, a point layer read as rainfall stations,
        # or any point layer treated as special delineation points).
        if not name:
            try:
                combo.setLayer(None)
                if combo.allowEmptyLayer() and combo.currentIndex() != 0:
                    combo.setCurrentIndex(0)
            except Exception:
                pass
            return
        for layer in QgsProject.instance().mapLayers().values():
            if layer.name() == name:
                combo.setLayer(layer)
                return

    # =========================================================================
    # Build cfg dict from widgets
    # =========================================================================
    def _collect_cfg(self):
        dem_layer = self.cb_dem.currentLayer()
        pts_layer = self.cb_points.currentLayer()

        if dem_layer is None:
            QMessageBox.warning(self, "Missing input", "Pick a DEM raster.")
            return None
        if pts_layer is None:
            QMessageBox.warning(self, "Missing input",
                                "Pick a pour-points layer.")
            return None
        if not isinstance(dem_layer, QgsRasterLayer):
            QMessageBox.warning(self, "Wrong layer type",
                                "The DEM must be a raster layer.")
            return None
        out_dir = self.fw_outdir.filePath().strip()
        if not out_dir:
            QMessageBox.warning(self, "Missing input",
                                "Pick an output directory.")
            return None

        thiessen_on = self.gb_thiessen.isChecked()
        stations_layer = (self.cb_stations.currentLayer()
                          if thiessen_on else None)
        if thiessen_on and stations_layer is None:
            QMessageBox.warning(
                self, "Missing input",
                "Spatially-variable rainfall is on \u2014 pick the "
                "rainfall-stations point layer (or untick the group).")
            return None

        rainfall = {}
        station_fields = {}
        seen_rps = set()
        for row in range(self.tbl_rainfall.rowCount()):
            def _cell(c_):
                it = self.tbl_rainfall.item(row, c_)
                return (it.text() if it else "").strip()
            rp_txt = _cell(0)
            if not rp_txt:
                continue
            try:
                rp = int(round(float(rp_txt.replace(",", "."))))
                if rp <= 0:
                    raise ValueError
            except ValueError:
                QMessageBox.warning(
                    self, "Bad return period",
                    f"Row {row + 1}: '{rp_txt}' is not a positive "
                    f"integer return period.")
                return None
            if rp in seen_rps:
                QMessageBox.warning(
                    self, "Duplicate return period",
                    f"Return period {rp} appears more than once.")
                return None

            if thiessen_on:
                fld = _cell(2)
                if not fld:
                    continue          # row not mapped -> skipped
                if stations_layer.fields().indexOf(fld) < 0:
                    QMessageBox.warning(
                        self, "Unknown station field",
                        f"Row {row + 1} ({rp}-yr): field '{fld}' does "
                        f"not exist in layer "
                        f"'{stations_layer.name()}'.")
                    return None
                station_fields[rp] = fld
                rainfall[rp] = 1.0    # placeholder; superseded by Thiessen
                seen_rps.add(rp)
            else:
                d_txt = _cell(1).replace(",", ".")
                if not d_txt:
                    continue
                try:
                    v = float(d_txt)
                except ValueError:
                    QMessageBox.warning(
                        self, "Bad rainfall depth",
                        f"Row {row + 1} ({rp}-yr): '{d_txt}' is not a "
                        f"number.")
                    return None
                if v > 0:
                    rainfall[rp] = v
                    seen_rps.add(rp)

        if thiessen_on and not station_fields:
            QMessageBox.warning(
                self, "No mapped return periods",
                "Spatially-variable rainfall is on but no row has a "
                "station field name in the third column.")
            return None

        custom_csv = self.fw_custom_csv.filePath().strip()
        if self.cb_dist.currentData() == "CUSTOM":
            if not custom_csv:
                QMessageBox.warning(
                    self, "Missing input",
                    "Custom hyetograph selected — pick the CSV file "
                    "(two columns: time in minutes, depth in mm).")
                return None
            if not os.path.isfile(custom_csv):
                QMessageBox.warning(
                    self, "File not found",
                    f"Custom hyetograph CSV does not exist:\n{custom_csv}")
                return None

        # v1.0 — Curve Number as raster OR constant value
        cn_is_value = (self.cb_cn_source.currentData() == "value")
        cn_layer = self.cb_cn.currentLayer() if not cn_is_value else None
        cn_name = cn_layer.name() if cn_layer is not None else None
        cn_value = (float(self.sb_cn_value.value())
                    if cn_is_value else None)
        point_id_field = self.cb_point_id.currentField() or None

        # ---- v23/v24: Special Delineation Points tab -------------------
        dams_layer = self.cb_dams.currentLayer()
        dam_name_field = ""
        dam_attr_fields = []
        sp_points_only = False
        sp_reservoirs = {}
        if dams_layer is not None:
            dam_name_field = (self.cb_dam_name.currentField() or "").strip()
            if not dam_name_field:
                QMessageBox.warning(
                    self, "Missing input",
                    "A Special Delineation Points layer is selected — "
                    "pick the field holding the point name in the "
                    "Special Delineation Points tab (or clear the "
                    "layer).")
                return None
            for i in range(self.lst_dam_fields.count()):
                it = self.lst_dam_fields.item(i)
                if (it.checkState() == Qt.Checked
                        and it.text() != dam_name_field):
                    dam_attr_fields.append(it.text())
            sp_points_only = self.chk_sp_only.isChecked()
            # persist the definition currently in the editor, then
            # keep only points actually flagged as reservoirs
            if self._sp_current:
                self._res_editor_save(self._sp_current)
            for nm, d in self._sp_res.items():
                if not d.get('enabled'):
                    continue
                sp_reservoirs[nm] = dict(
                    mode=d.get('mode', 'structures'),
                    ev_text=d.get('ev_text', ''),
                    structures=d.get('structures') or [],
                    table_text=d.get('table_text', ''),
                    init_mode=d.get('init_mode', 'inflow'),
                    init_elev=d.get('init_elev'))
            # early sanity check so the user learns about an unusable
            # reservoir before the (long) run, not after
            for nm, rdef in sp_reservoirs.items():
                _sd = reservoir_sd(rdef)
                if not _sd or len(_sd.get('S') or []) < 2:
                    QMessageBox.warning(
                        self, "Incomplete reservoir definition",
                        f"Reservoir '{nm}': "
                        f"{(_sd or {}).get('warn') or 'undefined'}.\n\n"
                        "Fix it in the Special Delineation Points tab, "
                        "or untick 'Define this point as a Reservoir'. "
                        "(If you run anyway, the point is still "
                        "enforced as a junction but flow passes "
                        "through unrouted.)")
                    return None

        cfg = {
            # ---- inputs ----
            "DEM_LAYER":          dem_layer.name(),
            "POINTS_LAYER":       pts_layer.name(),
            "POINT_ID_FIELD":     point_id_field,
            "CN_LAYER":           cn_name,
            "CN_VALUE":           cn_value,
            "OUTPUT_DIR":         out_dir,

            # ---- delineation ----
            "DEPRESSION_EPSILON": float(self.sb_eps.value()),
            "STREAM_THRESHOLD":   int(self.sb_thresh.value()),
            "SNAP_RADIUS_CELLS":  int(self.sb_snap.value()),
            "BASIN_SLOPE_FROM":   self.cb_slope_from.currentData(),
            "SUB_BASIN_AREA_KM2": float(self.sb_subarea.value()),

            # ---- special delineation points (v23/v24) ----
            "DAMS_LAYER":         (dams_layer.name()
                                   if dams_layer is not None else ""),
            "DAM_NAME_FIELD":     dam_name_field,
            "DAM_ATTR_FIELDS":    dam_attr_fields,
            "SP_POINTS_ONLY":     sp_points_only,
            "SP_RESERVOIRS":      sp_reservoirs,

            # ---- rainfall + hydrographs ----
            "RAINFALL_DEPTHS_MM":  rainfall,
            "TC_METHOD":           self.cb_tc.currentData(),
            "PEAK_RATE_FACTOR":    int(self.cb_prf.currentData()),
            "TC_MIN_MINUTES":      float(self.sb_tc_min.value()),
            "SCS_LAMBDA":          float(self.sb_lambda.value()),
            "HYDROGRAPHS_ENABLED": bool(self.chk_hyd.isChecked()),
            "STORM_DURATION_HR":   float(self.sb_storm_hr.value()),
            "HYD_TIME_STEP_MIN":   float(self.sb_dt_min.value()),
            "HYD_TAIL_FACTOR":     float(self.sb_tail_factor.value()),
            "HYD_MAX_TAIL_HR":     float(self.sb_tail_cap.value()),
            "STORM_DISTRIBUTION":  self.cb_dist.currentData(),
            "CUSTOM_HYETOGRAPH_CSV": custom_csv,
            "RAIN_STATIONS_LAYER":  (stations_layer.name()
                                     if thiessen_on else ""),
            "RAIN_STATIONS_FIELDS": station_fields,

            # ---- ARF & reach routing ----
            "ARF_METHOD":           self.cb_arf.currentData(),
            "ARF_FIXED_VALUE":      float(self.sb_arf_fixed.value()),
            "REACH_ROUTING_METHOD": self.cb_route.currentData(),
            "REACH_TC_MIN_MINUTES": float(
                self.sb_reach_tc_min.value()),
            "MUSKINGUM_X":          float(self.sb_musk_x.value()),
            "MUSK_SUBREACHES":     int(self.sb_musk_nsub.value()),
            "REACH_TC_METHOD":      self.cb_reach_tc.currentData(),
            "MC_MANNING_N":         float(self.sb_mc_n.value()),
            "MC_BOTTOM_WIDTH_M":    float(self.sb_mc_b.value()),
            "MC_SIDE_SLOPE":        float(self.sb_mc_z.value()),

            # ---- HEC-HMS export (v19) ----
            # ---- AR-Awdiya project export (v1.0) ----
            "EXPORT_AWDIYA":       self.chk_awdiya.isChecked(),
            "AWDIYA_PROJECT_NAME": self.ed_awdiya_name.text().strip(),
            "AWDIYA_FREQ_MODE":    self.cb_awdiya_freq.currentData(),
            "AWDIYA_PRIMARY_RP":   (int(self.cb_awdiya_rp.currentData())
                                    or None),
            "EXPORT_HEC_HMS":       bool(self.chk_hms.isChecked()),
            "EXPORT_HEC_HMS_MET":   bool(self.chk_hms.isChecked()
                                         and self.chk_hms_met.isChecked()),
            "HEC_HMS_NAME":        self.ed_hms_name.text().strip(),
            "HEC_HMS_VERSION":      self.cb_hms_ver.currentData(),

            # ---- styling ----
            "STYLES_DIR":           self.fw_styles_dir.filePath().strip(),

            # ---- output names ----
            "OUT_STREAMS":     self.le_out_streams.text().strip() or "Streams",
            "OUT_CATCHMENTS":  self.le_out_catch.text().strip()   or "Catchments",
            "OUT_LFP_LINES":   self.le_out_lfp_pl.text().strip()  or "LFP_PL",
            "OUT_LFP_POINTS":  self.le_out_lfp_pt.text().strip()  or "LFP_PTs",
            "OUT_SNAPPED":     self.le_out_snap.text().strip()    or "Snapped",
            "OUT_THIESSEN":        self.le_out_thies.text().strip(),
            "OUT_SUBBASINS":   self.le_out_subs.text().strip()    or "SubBasins",
            "OUT_SUBLFPS":     self.le_out_sublfp.text().strip()  or "SubLFPs",
            "OUT_REACHES":     self.le_out_reach.text().strip()   or "Reaches",
            "OUT_JUNCTIONS":   self.le_out_junc.text().strip()    or "Junctions",
        }
        return cfg

    # =========================================================================
    # Run
    def reject(self):
        """v1.0 — Esc / Cancel while a run is active is ignored (see
        closeEvent)."""
        if self._task is not None:
            self.iface.messageBar().pushMessage(
                "AR-Tarseem",
                "A run is still in progress — the dialog will close "
                "itself when the run finishes.",
                level=Qgis.Info, duration=4)
            return
        super().reject()

    def closeEvent(self, event):
        """v1.0 — while a background run is active, keep the dialog
        alive: the completion callback needs it, and closing mid-run
        would orphan the progress reporting."""
        if self._task is not None:
            self.iface.messageBar().pushMessage(
                "AR-Tarseem",
                "A run is still in progress — the dialog will close "
                "itself when the run finishes.",
                level=Qgis.Info, duration=4)
            event.ignore()
            return
        super().closeEvent(event)

    # =========================================================================
    def on_run(self):
        cfg = self._collect_cfg()
        if cfg is None:
            return

        try:
            os.makedirs(cfg["OUTPUT_DIR"], exist_ok=True)
        except OSError as e:
            QMessageBox.critical(self, "Output directory",
                                 f"Can't create '{cfg['OUTPUT_DIR']}':\n{e}")
            return

        QgsMessageLog.logMessage("─── Run starting ───", LOG_TAG, Qgis.Info)
        self.setEnabled(False)
        self.iface.messageBar().pushMessage(
            "AR-Tarseem",
            "Running in the background — QGIS stays responsive. "
            "Progress: task bar / Log Messages → \"AR-Tarseem\".",
            level=Qgis.Info, duration=5)

        # v1.0 — the engine runs as a background QgsTask: QGIS remains
        # fully responsive, progress is reported live (native task
        # manager + this dialog), and the produced layers are
        # registered on the GUI thread when the task finishes.
        # hydro_core itself is untouched — see hydro_task.py for the
        # deferred-registration mechanism. The progress dialog is
        # non-modal so the map canvas stays usable during the run.
        # v1.0 — parent the progress dialog to the QGIS MAIN WINDOW
        # (not to this settings dialog, which is hidden during the
        # run) and keep it strictly non-modal so the canvas, layer
        # panel and attribute tables remain fully interactive.
        prog = QProgressDialog("Preparing…", None, 0, 100,
                               self.iface.mainWindow())
        prog.setWindowTitle("AR-Tarseem — running")
        prog.setWindowModality(Qt.NonModal)
        prog.setCancelButton(None)          # run is not interruptible
        prog.setMinimumDuration(0)
        prog.setMinimumWidth(420)
        prog.setAutoClose(False)
        prog.setAutoReset(False)
        prog.setValue(0)
        prog.show()

        def _on_done(ok, err_text):
            try:
                prog.close()
                prog.deleteLater()
            except Exception:
                pass
            self._task = None
            try:
                self.setEnabled(True)
            except Exception:
                return                      # dialog was closed meanwhile
            if ok:
                msg = "Done. Layers added; outputs in Modeling_Outputs/."
                if cfg.get("EXPORT_HEC_HMS"):
                    msg = ("Done. Layers added; HEC-HMS .basin in "
                           "<output>/HEC_HMS/.")
                self.iface.messageBar().pushMessage(
                    "AR-Tarseem", msg, level=Qgis.Success, duration=6)
                self.accept()
            else:
                self.show()            # bring the settings back
                self.raise_()
                self.activateWindow()
                if err_text.strip():
                    QgsMessageLog.logMessage(err_text, LOG_TAG,
                                             Qgis.Critical)
                last = (err_text.strip().splitlines()[-1]
                        if err_text.strip() else "Unknown error.")
                QMessageBox.critical(
                    self, "AR-Tarseem — error",
                    f"The run failed:\n\n{last}\n\n"
                    "Full traceback in Log Messages → AR-Tarseem.")

        self._task = ARTarseemRunTask(cfg, _redirect_prints_to_log,
                                      _on_done)
        self._task.progressChanged.connect(
            lambda p: prog.setValue(int(round(p))))
        self._task.bridge.message.connect(prog.setLabelText)
        QgsApplication.taskManager().addTask(self._task)
        # v1.0 — hide the settings dialog while the engine runs; even
        # disabled it would sit in front of the map. Only the small
        # non-modal progress window remains, and QGIS is fully usable.
        self.hide()
