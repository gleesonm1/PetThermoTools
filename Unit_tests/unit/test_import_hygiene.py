"""
Tier A regression test for the `petthermotools.Path` submodule shadowing bug.

`petthermotools/__init__.py` does `from petthermotools.Path import *` for
every submodule, with no `__all__` anywhere to scope what gets re-exported.
Path.py itself does `from pathlib import Path`, so that star-import
overwrites the `Path` attribute Python's import machinery normally sets on
the `petthermotools` package after importing the `petthermotools.Path`
submodule. The net effect: `import petthermotools.Path as x` and
`from petthermotools import Path` both silently resolve to `pathlib.Path`
instead of the calculation-path submodule.

This was flagged as a suspected-but-untested risk during a code audit; it is
confirmed here. Fixing it (Part(i)#6 of the cleanup plan) is scoped as
*non-breaking*: nothing below should be "fixed" by renaming the submodule or
changing what `ptt.Path` means for existing users, only by not letting an
internal star-import clobber it -- so this test pins the bug rather than
prescribing the fix.
"""
import importlib
import pathlib

import petthermotools as ptt


def test_petthermotools_dot_path_attribute_is_shadowed_by_pathlib():
    """Characterizes the current (buggy) state: ptt.Path is pathlib.Path,
    not the petthermotools.Path submodule. If this ever starts failing, the
    shadowing has been fixed -- update/remove this test and the workaround
    in test_multi_path_supervisor.py at the same time."""
    assert ptt.Path is pathlib.Path


def test_petthermotools_path_submodule_is_reachable_via_sys_modules():
    """The real submodule is never lost -- it is registered in
    sys.modules under its dotted name regardless of the __init__.py
    star-import clobbering the package attribute, so this is the reliable
    way to reach it today."""
    real_path_module = importlib.import_module("petthermotools.Path")
    assert hasattr(real_path_module, "multi_path")
    assert real_path_module is not pathlib.Path
