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
"""Qt 5 / Qt 6 compatibility layer — one plugin, both QGIS generations.

QGIS 4.x ("Norrköping", March 2026 onward) is built on Qt 6 and binds
through **PyQt6**, which differs from PyQt5 in two ways that break almost
every QGIS 3 plugin:

1. **Scoped enumerations.**  PyQt5 exposed enum members directly on the
   owning class (``Qt.Checked``, ``QgsFileWidget.GetDirectory``).  PyQt6
   exposes them only on the nested enumeration type
   (``Qt.CheckState.Checked``, ``QgsFileWidget.StorageMode.GetDirectory``),
   and the unscoped spelling raises ``AttributeError``.

2. **No int → enum coercion.**  PyQt5 silently accepted a plain integer
   wherever an enum or a flags object was expected; PyQt6 raises
   ``TypeError``.

3. **``QVariant.Type`` is gone.**  Field types for ``QgsField`` must be
   given as ``QMetaType.Type`` members under Qt 6.

Rather than fork the code base, every affected symbol is resolved **once,
at import time**, through the helpers below.  Each resolver tries the
Qt 6 (scoped) spelling first and falls back to the Qt 5 (unscoped) one, so
the same module runs unmodified on QGIS 3.22 → 3.4x and on QGIS 4.x.

Nothing in this module imports anything that does not exist in both
generations, so importing it can never be the thing that fails.
"""

from __future__ import annotations

from qgis.PyQt import QtCore
from qgis.PyQt.QtCore import Qt
from qgis.PyQt.QtWidgets import QAbstractItemView, QDialogButtonBox

# --------------------------------------------------------------------------
# Which Qt are we on?
# --------------------------------------------------------------------------
try:
    QT_VERSION_INT = int(QtCore.QT_VERSION)
except Exception:                                    # pragma: no cover
    try:
        QT_VERSION_INT = int(QtCore.QT_VERSION_STR.split('.')[0]) << 16
    except Exception:
        QT_VERSION_INT = 0x050000

QT_MAJOR = (QT_VERSION_INT >> 16) & 0xFF
IS_QT6 = QT_MAJOR >= 6
IS_QT5 = not IS_QT6

QT_VERSION_STR = getattr(QtCore, 'QT_VERSION_STR', '%d.x' % QT_MAJOR)


# --------------------------------------------------------------------------
# Generic enum resolution
# --------------------------------------------------------------------------
def enum_member(owner, nested, name, default=None, required=True):
    """Return ``owner.nested.name`` (Qt 6) or ``owner.name`` (Qt 5).

    Parameters
    ----------
    owner : type
        The class that owns the enumeration (``Qt``, ``QgsFileWidget`` …).
    nested : str
        Name of the nested enumeration type as Qt 6 / PyQt6 spells it
        (``'CheckState'``, ``'StorageMode'`` …).  May be ``None`` to skip
        straight to the scan.
    name : str
        The member name, which is identical in both generations.
    required : bool
        When False, return ``default`` instead of raising if nothing
        matches — used for symbols that genuinely do not exist in one of
        the two generations.

    The last resort is a scan of the owner's nested types, which covers
    classes whose enumeration was *renamed* between Qt 5 and Qt 6.
    """
    # 1. Qt 6 / scoped spelling.
    if nested:
        scope = getattr(owner, nested, None)
        if scope is not None:
            member = getattr(scope, name, None)
            if member is not None:
                return member

    # 2. Qt 5 / unscoped spelling (also hit on PyQt6 builds that keep the
    #    old aliases for backward compatibility).
    member = getattr(owner, name, None)
    if member is not None:
        return member

    # 3. Scan the owner's nested enumeration types.
    for attr in dir(owner):
        if attr.startswith('_'):
            continue
        try:
            scope = getattr(owner, attr)
        # B112: probing attributes of a C++ binding; unreadable ones are skipped
        except Exception:  # nosec B112
            continue
        if isinstance(scope, type):
            member = getattr(scope, name, None)
            if member is not None:
                return member

    if required:
        raise AttributeError(
            "AR-Tarseem: cannot resolve %s.%s on Qt %s — please report this "
            "with your QGIS version." % (getattr(owner, '__name__', owner),
                                         name, QT_VERSION_STR))
    return default


# --------------------------------------------------------------------------
# Field types:  QVariant.Type (Qt 5)  ->  QMetaType.Type (Qt 6)
# --------------------------------------------------------------------------
# ``QgsField`` has taken a ``QMetaType::Type`` since QGIS 3.38 and takes
# ONLY that from QGIS 4.  On Qt 5 we keep ``QVariant.Type`` so the plugin
# still installs on QGIS 3.22, which predates the QMetaType overload.
_QMETATYPE_NAMES = {
    'Int':      'Int',
    'LongLong': 'LongLong',
    'Double':   'Double',
    'String':   'QString',
    'Bool':     'Bool',
    'Date':     'QDate',
    'DateTime': 'QDateTime',
}


def _build_field_types():
    if IS_QT6:
        from qgis.PyQt.QtCore import QMetaType
        scope = getattr(QMetaType, 'Type', QMetaType)
        return {k: getattr(scope, v) for k, v in _QMETATYPE_NAMES.items()}
    from qgis.PyQt.QtCore import QVariant
    return {k: getattr(QVariant, k) for k in _QMETATYPE_NAMES}


_FIELD_TYPES = _build_field_types()

#: Field-type constants to hand to ``QgsField(...)`` on either generation.
FT_INT = _FIELD_TYPES['Int']
FT_LONGLONG = _FIELD_TYPES['LongLong']
FT_DOUBLE = _FIELD_TYPES['Double']
FT_STRING = _FIELD_TYPES['String']
FT_BOOL = _FIELD_TYPES['Bool']


# --------------------------------------------------------------------------
# Qt core enums used by the dialog
# --------------------------------------------------------------------------
QT_CHECKED = enum_member(Qt, 'CheckState', 'Checked')
QT_UNCHECKED = enum_member(Qt, 'CheckState', 'Unchecked')
QT_ITEM_IS_USER_CHECKABLE = enum_member(Qt, 'ItemFlag', 'ItemIsUserCheckable')
QT_RICH_TEXT = enum_member(Qt, 'TextFormat', 'RichText')
QT_NON_MODAL = enum_member(Qt, 'WindowModality', 'NonModal')
QT_WAIT_CURSOR = enum_member(Qt, 'CursorShape', 'WaitCursor')
QT_ALIGN_TOP = enum_member(Qt, 'AlignmentFlag', 'AlignTop')

#: ``QListWidget.NoSelection`` lives on ``QAbstractItemView`` in Qt 6.
NO_SELECTION = enum_member(QAbstractItemView, 'SelectionMode', 'NoSelection')

ACCEPT_ROLE = enum_member(QDialogButtonBox, 'ButtonRole', 'AcceptRole')
REJECT_ROLE = enum_member(QDialogButtonBox, 'ButtonRole', 'RejectRole')


# --------------------------------------------------------------------------
# QGIS enums
# --------------------------------------------------------------------------
from qgis.core import Qgis, QgsMapLayerProxyModel        # noqa: E402

MSG_INFO = enum_member(Qgis, 'MessageLevel', 'Info')
MSG_WARNING = enum_member(Qgis, 'MessageLevel', 'Warning')
MSG_CRITICAL = enum_member(Qgis, 'MessageLevel', 'Critical')
MSG_SUCCESS = enum_member(Qgis, 'MessageLevel', 'Success')

# QGIS 4 moved the layer-filter members onto a nested ``Filter`` type;
# QGIS 3 spelled them straight on the class.  Both are covered.
FILTER_RASTER = enum_member(QgsMapLayerProxyModel, 'Filter', 'RasterLayer')
FILTER_POINT = enum_member(QgsMapLayerProxyModel, 'Filter', 'PointLayer')

# ``qgis.gui`` is absent from headless / qgis_process contexts, where the
# engine is still importable and useful; these two constants are only
# needed by the dialog, so a missing GUI library must not break the import.
try:
    from qgis.gui import QgsFileWidget                    # noqa: E402

    GET_DIRECTORY = enum_member(QgsFileWidget, 'StorageMode', 'GetDirectory')
    GET_FILE = enum_member(QgsFileWidget, 'StorageMode', 'GetFile')
except ImportError:                                       # pragma: no cover
    GET_DIRECTORY = None
    GET_FILE = None


# --------------------------------------------------------------------------
# Misc helpers
# --------------------------------------------------------------------------
def exec_dialog(dlg):
    """``exec()`` on Qt 6, ``exec_()`` on Qt 5 (both exist on recent PyQt5,
    but ``exec_`` was dropped in PyQt6)."""
    runner = getattr(dlg, 'exec', None) or getattr(dlg, 'exec_')
    return runner()


__all__ = [
    'IS_QT5', 'IS_QT6', 'QT_MAJOR', 'QT_VERSION_STR', 'enum_member',
    'FT_INT', 'FT_LONGLONG', 'FT_DOUBLE', 'FT_STRING', 'FT_BOOL',
    'QT_CHECKED', 'QT_UNCHECKED', 'QT_ITEM_IS_USER_CHECKABLE',
    'QT_RICH_TEXT', 'QT_NON_MODAL', 'QT_WAIT_CURSOR', 'QT_ALIGN_TOP',
    'NO_SELECTION', 'ACCEPT_ROLE', 'REJECT_ROLE',
    'MSG_INFO', 'MSG_WARNING', 'MSG_CRITICAL', 'MSG_SUCCESS',
    'FILTER_RASTER', 'FILTER_POINT', 'GET_DIRECTORY', 'GET_FILE',
    'exec_dialog',
]
