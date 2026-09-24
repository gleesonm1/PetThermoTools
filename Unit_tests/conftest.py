"""
Shared pytest configuration for the PetThermoTools test suite.

Test tiers
----------
- ``Unit_tests/unit/``        : Tier A. Pure Python, no alphaMELTS or Julia
  required. Must run on every OS/Python version, including CI.
- ``Unit_tests/integration/alphamelts/``    : Tier B1. Needs a working alphaMELTS
  (``meltsdynamic``) install. Skipped automatically if not present.
- ``Unit_tests/integration/magemin/``  : Tier B2. Needs a working Julia +
  MAGEMinCalc install (``~/.petthermotools_julia_env``). Skipped
  automatically if not present.

Per team decision, tiers B1/B2 are local-only for now (not run in CI).
"""
import importlib.util
import sys
from pathlib import Path

import pytest

# Allow running the test suite against the working tree without installing
# the package first (`pip install -e .`), by putting src/ on sys.path.
_SRC = Path(__file__).resolve().parent.parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

_JULIA_ENV_DIR = Path.home() / ".petthermotools_julia_env"


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "melts: requires a working alphaMELTS (meltsdynamic) install"
    )
    config.addinivalue_line(
        "markers", "magemin: requires a working Julia + MAGEMinCalc install"
    )


def _melts_importable() -> bool:
    # NOTE: deliberately a find_spec() check, not a real `import meltsdynamic`.
    # This function runs at collection time for every test session, and
    # actually importing/initializing the C engine here would be an
    # unwanted side effect on a check that is only supposed to answer
    # "is it available".
    try:
        return importlib.util.find_spec("meltsdynamic") is not None
    except (ImportError, ValueError):
        return False


def _magemin_available() -> bool:
    # NOTE: deliberately a find_spec() check, not `import juliacall`.
    # Importing juliacall embeds Julia in this process, which is exactly the
    # state Path.py:174-175 treats as "MAGEMin already initiated" and
    # refuses to run MELTS after (on Windows). A collection-time
    # availability probe must not trigger that side effect itself, or every
    # MELTS test in the same session would start failing on Windows.
    try:
        if importlib.util.find_spec("juliacall") is None:
            return False
    except (ImportError, ValueError):
        return False
    # find_spec() only proves the Python-side package is installed; PTT's
    # own MAGEMin backend also needs its dedicated Julia project env to
    # exist (see GenFuncs._ensure_julia_ready).
    return _JULIA_ENV_DIR.exists()


@pytest.fixture(scope="session")
def melts_available() -> bool:
    """True if alphaMELTS (meltsdynamic) can be imported in this environment."""
    return _melts_importable()


@pytest.fixture(scope="session")
def magemin_available() -> bool:
    """True if Julia + MAGEMinCalc are set up in this environment."""
    return _magemin_available()


def pytest_collection_modifyitems(config, items):
    """Auto-skip integration tests when their backend isn't available, so
    `pytest` run with no args does the right thing everywhere (CI included)
    without every test file needing its own importorskip boilerplate."""
    melts_ok = _melts_importable()
    magemin_ok = _magemin_available()

    skip_melts = pytest.mark.skip(reason="alphaMELTS (meltsdynamic) not available")
    skip_magemin = pytest.mark.skip(reason="Julia + MAGEMinCalc not available")

    for item in items:
        if "melts" in item.keywords and not melts_ok:
            item.add_marker(skip_melts)
        if "magemin" in item.keywords and not magemin_ok:
            item.add_marker(skip_magemin)
