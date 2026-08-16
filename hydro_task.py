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
"""Background execution of the AR-Tarseem engine via QgsTask.

Design constraint honoured here: ``hydro_core`` is NOT modified in any
way. The engine keeps calling ``QgsProject.instance().addMapLayer(...)``
exactly as it always did; this module makes that safe from a worker
thread by injecting a thin project proxy (module-attribute patch, not a
code change) that:

  1. moves each produced layer's thread affinity to the main thread
     (``moveToThread`` — legal here because it is invoked from the
     thread that currently owns the object, i.e. the worker), and
  2. DEFERS the actual project registration: layers are collected
     during the run and registered on the MAIN thread inside
     ``QgsTask.finished()``, which QGIS guarantees to execute on the
     GUI thread.

Everything the engine does in between (feature writing, QML styling of
the not-yet-registered layers, GDAL I/O, NumPy work) is thread-safe or
touches only worker-owned / not-yet-rendered objects.

Benefits delivered: QGIS stays fully responsive during a run, progress
is reported through the native task manager (status-bar spinner + our
progress dialog), and a failure in the engine is contained and reported
instead of freezing the interface.
"""

from __future__ import annotations

import sys
import traceback

from qgis.core import (QgsApplication, QgsCoordinateTransform,
                       QgsProject, QgsTask)
from qgis.PyQt.QtCore import QObject, pyqtSignal

from . import hydro_core


class _DeferredProjectProxy:
    """Stands in for the live QgsProject during a background run.

    Read access (``mapLayersByName``, ``transformContext``, CRS
    queries, ...) passes straight through to the real project — those
    are the only project reads the engine performs, and they are safe
    from the worker. ``addMapLayer`` is intercepted: the layer is
    re-affined to the main thread and queued for registration after
    the task finishes."""

    def __init__(self, real_project):
        self._real = real_project
        self.pending_layers = []

    def __getattr__(self, name):
        return getattr(self._real, name)

    def addMapLayer(self, layer, *args, **kwargs):
        try:
            layer.moveToThread(QgsApplication.instance().thread())
        except Exception:
            pass
        self.pending_layers.append(layer)
        return layer


class _ShimProjectClass:
    """Replaces the ``QgsProject`` *name* inside hydro_core's namespace
    for the duration of one background run (restored afterwards), so
    ``QgsProject.instance()`` yields the deferred proxy. hydro_core's
    source is untouched."""

    _proxy = None

    @classmethod
    def instance(cls):
        return cls._proxy


class _XformShim:
    """Stands in for the ``QgsCoordinateTransform`` *name* inside
    hydro_core during a background run. The engine constructs
    transforms as ``QgsCoordinateTransform(src_crs, dst_crs, project)``
    — but ``project`` is our deferred proxy during the run, which the
    sip bindings reject ("argument has unexpected type"). This shim
    swaps the proxy for the REAL project in any argument position and
    delegates everything else (including class attributes / enums) to
    the genuine class."""

    def __init__(self, proxy):
        self._proxy = proxy

    def __call__(self, *args, **kwargs):
        args = tuple(self._proxy._real if a is self._proxy else a
                     for a in args)
        kwargs = {k: (self._proxy._real if v is self._proxy else v)
                  for k, v in kwargs.items()}
        return QgsCoordinateTransform(*args, **kwargs)

    def __getattr__(self, name):
        return getattr(QgsCoordinateTransform, name)


class TaskBridge(QObject):
    """Cross-thread signal bridge: worker emits, GUI thread receives
    (Qt auto-queues signals across threads)."""
    message = pyqtSignal(str)


class ARTarseemRunTask(QgsTask):
    """QgsTask wrapper around hydro_core.run().

    Parameters
    ----------
    cfg : dict
        The configuration collected by the dialog.
    log_ctx_factory : callable
        Zero-argument callable returning the print-redirection context
        manager (the dialog's ``_redirect_prints_to_log``); passed in
        to avoid a circular import.
    on_done : callable(ok: bool, err_text: str)
        Invoked on the MAIN thread from ``finished()`` after the
        deferred layers have been registered.
    """

    def __init__(self, cfg, log_ctx_factory, on_done):
        # flags = empty Flags(): no Cancel button — the engine is not
        # interruptible mid-run (by design: hydro_core is untouched);
        # offering a cancel that silently keeps computing would
        # mislead users. NOTE: PyQGIS requires a QgsTask.Flags object
        # here — a plain int raises TypeError on QGIS 3.3x.
        try:
            _flags = QgsTask.Flags()
        except Exception:                     # very old/new sip variants
            _flags = QgsTask.Flag(0)
        super().__init__("AR-Tarseem — watershed modelling run", _flags)
        self._cfg = cfg
        self._log_ctx_factory = log_ctx_factory
        self._on_done = on_done
        self._err_text = ""
        self._proxy = None
        self.bridge = TaskBridge()

    # ---- worker thread ------------------------------------------------
    def run(self):
        saved_qgsproject = hydro_core.QgsProject
        saved_xform = hydro_core.QgsCoordinateTransform
        # v1.0 — GUI smoothness under the GIL: the engine's pure-
        # Python loops hold the interpreter lock, and the GUI thread
        # only gets scheduled every switch interval (default 5 ms).
        # Shortening it during the run keeps the canvas and panels
        # fluid at a negligible throughput cost; restored afterwards.
        saved_switch = sys.getswitchinterval()
        try:
            sys.setswitchinterval(0.001)
        except Exception:
            pass
        self._proxy = _DeferredProjectProxy(QgsProject.instance())
        _ShimProjectClass._proxy = self._proxy
        hydro_core.QgsProject = _ShimProjectClass
        hydro_core.QgsCoordinateTransform = _XformShim(self._proxy)
        try:
            def _cb(frac, msg):
                # QgsTask.setProgress is safe from the worker; the
                # message goes through a queued cross-thread signal.
                self.setProgress(max(0.0, min(100.0, 100.0 * frac)))
                try:
                    self.bridge.message.emit(str(msg))
                except Exception:
                    pass

            with self._log_ctx_factory():
                hydro_core.run(self._cfg, progress_cb=_cb)
            return True
        except Exception as exc:
            self._err_text = "".join(traceback.format_exception(
                type(exc), exc, exc.__traceback__))
            return False
        finally:
            hydro_core.QgsProject = saved_qgsproject
            hydro_core.QgsCoordinateTransform = saved_xform
            _ShimProjectClass._proxy = None
            try:
                sys.setswitchinterval(saved_switch)
            except Exception:
                pass

    # ---- main (GUI) thread --------------------------------------------
    def finished(self, ok):
        n_added = 0
        try:
            if ok and self._proxy is not None:
                real = QgsProject.instance()
                for lyr in self._proxy.pending_layers:
                    try:
                        real.addMapLayer(lyr)
                        n_added += 1
                    except Exception as e:      # keep going; report once
                        self._err_text += (
                            f"\n⚠ could not register layer "
                            f"{getattr(lyr, 'name', lambda: '?')()}: {e}")
        finally:
            if self._proxy is not None:
                self._proxy.pending_layers = []
            try:
                self._on_done(bool(ok), self._err_text)
            except Exception:
                pass

    @property
    def error_text(self):
        return self._err_text
