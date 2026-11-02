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
"""AR-Tarseem — QGIS plugin entry point.

Copyright (c) 2026 Ahmed M. Raslan, PhD & Amr R. Refaiy, PhD - Water Engineering Experts

QGIS calls ``classFactory(iface)`` once at startup. Everything else
(dialog, core algorithm) is loaded lazily inside the plugin class.
"""


def classFactory(iface):  # noqa: N802 — QGIS-mandated name
    from .hydro_plugin import HydroDelineationPlugin
    return HydroDelineationPlugin(iface)
