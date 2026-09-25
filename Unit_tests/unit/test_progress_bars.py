"""
Tier A regression test: progress bars must work outside Jupyter.

Found by the real-alphaMELTS golden test on a clean Windows install:
`findLiq_multi` crashed with "ImportError: IProgress not found. Please update
jupyter and ipywidgets". Five modules did `from tqdm.notebook import tqdm,
trange` (Liq, Melting, PhaseDiagrams, Path, Saturation); tqdm.notebook needs
ipywidgets, which setup.py does not declare, so any environment without it
(a clean `pip install .`, CI, a terminal or script user) crashed in every code
path that actually draws a bar (findLiq_multi, findCO2_multi, the legacy
findSatPressure). Machines that happen to have Jupyter installed never saw it.

FIXED (Phase 1): the modules now use `from tqdm.auto import tqdm, trange`, which
gives the notebook widget bar in Jupyter (with ipywidgets) and a text bar
everywhere else. This is the regression test.

The absence of ipywidgets is simulated by patching tqdm.notebook's widget
lookup, so the outcome is the same on every machine, including ones that do
have ipywidgets.
"""
import gc
import importlib

import pytest

_MODULES = ["Liq", "Melting", "PhaseDiagrams", "Path", "Saturation"]


@pytest.mark.parametrize("module", _MODULES)
def test_progress_bar_works_without_ipywidgets(module, monkeypatch):
    tqdm_notebook = importlib.import_module("tqdm.notebook")
    monkeypatch.setattr(tqdm_notebook, "IProgress", None)  # what tqdm sees without ipywidgets
    # The bar that fails to construct is left half-built, and tqdm's __del__ then
    # calls close(), which raises AttributeError ('disp'). Left to the garbage
    # collector that runs at an arbitrary later moment -- e.g. inside pytest's
    # own ast.parse when it formats another test's failure, where it crashed the
    # whole session with `SystemError: AST constructor recursion depth mismatch`
    # on Python 3.11.0. So neutralise close() and collect while still patched.
    monkeypatch.setattr(tqdm_notebook.tqdm, "close", lambda self, *args, **kwargs: None)

    mod = importlib.import_module(f"petthermotools.{module}")
    failure = None
    try:
        for _ in mod.tqdm(range(2)):
            pass
    except ImportError as exc:
        failure = f"petthermotools.{module}.tqdm raised {type(exc).__name__}: {exc}"
    gc.collect()
    if failure:
        # Raised outside the except block: the ImportError (and, through its
        # traceback, the half-built bar) is not kept alive by the report.
        raise AssertionError(failure)
