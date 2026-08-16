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
"""Plugin class: registers a menu/toolbar action and opens the dialog.

Copyright (c) 2026 Ahmed M. Raslan, PhD & Amr R. Refaiy, PhD - Water Engineering Experts
"""

from __future__ import annotations

import os

from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction


class HydroDelineationPlugin:
    """Thin wrapper around QGIS plugin lifecycle hooks."""

    MENU_NAME = "&Hydrology"
    ACTION_TEXT = "AR-Tarseem…"

    def __init__(self, iface):
        self.iface = iface
        self.action: QAction | None = None
        self._dialog = None

    def initGui(self):
        icon_path = os.path.join(os.path.dirname(__file__), "icon.png")
        self.action = QAction(QIcon(icon_path), self.ACTION_TEXT,
                              self.iface.mainWindow())
        self.action.setObjectName("ARTarseemAction")
        self.action.triggered.connect(self.show_dialog)
        self.iface.addPluginToMenu(self.MENU_NAME, self.action)
        self.iface.addToolBarIcon(self.action)

    def unload(self):
        if self.action is not None:
            self.iface.removePluginMenu(self.MENU_NAME, self.action)
            self.iface.removeToolBarIcon(self.action)
            self.action.deleteLater()
            self.action = None
        if self._dialog is not None:
            try:
                self._dialog.close()
            except Exception:
                pass
            self._dialog = None

    def show_dialog(self):
        from .hydro_dialog import HydroDelineationDialog
        self._dialog = HydroDelineationDialog(self.iface,
                                              self.iface.mainWindow())
        # v1.0 — non-modal: exec_() would make the dialog
        # application-modal, blocking all interaction with the QGIS
        # main window even while the engine runs as a background task.
        self._dialog.show()
        self._dialog.raise_()
        self._dialog.activateWindow()
