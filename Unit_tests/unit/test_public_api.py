"""
Tier A: the public API freeze.

A guard rail for MAINTAINERS, not a restriction on users: it lives in Unit_tests/
(not shipped with the package), changes nothing at run time, and only fails when a
change to the code would break something users call. Adding functions, adding
optional parameters at the end, changing behaviour or defaults, and touching
underscore names are all free.

Two tiers (helpers in _api_scan.py):

1. Documented usage. Every MAINTAINED notebook is scanned and must only use
   package names that still exist, with arguments the functions still accept.
   "Maintained" is what ReadTheDocs publishes (the toctrees in docs/index.rst) plus
   Workshops/, docs/teaching_materials/ and docs/PaperFigures/; notebooks in
   docs/conf.py `exclude_patterns` are retired and ignored. ReadTheDocs never runs
   notebooks (nbsphinx_execute = 'never'), so this is the only automatic check that
   the published examples still match the code.

2. The snapshot (Unit_tests/api/public_api.json): every public name today and the
   parameters of every function and class defined in the package. Removing or
   renaming a name or parameter, reordering leading positional parameters, making an
   optional parameter required, or adding a new required parameter fails. Names that
   only leak in from other libraries (`ptt.np`, ...) are recorded but not enforced;
   names from optional packages (alphaMELTS, ...) are never required.

A deliberate change to the public surface is made by re-running
Unit_tests/api/capture_public_api.py and committing the diff (see Unit_tests/README.md).
"""
import inspect
import warnings

import pytest
import _api_scan
from _api_scan import (
    OPTIONAL_PACKAGES, REPO, SNAPSHOT, SNAPSHOT_FORMAT, CallUse, check_call, compare_surface,
    encode_param, is_excluded_from_docs, load_snapshot, maintained_notebooks,
    public_surface, scan_notebook, toctree_notebooks,
)

warnings.simplefilter("ignore")
import petthermotools as ptt  # noqa: E402

NOTEBOOKS = maintained_notebooks()


def _relative(path):
    return path.relative_to(REPO.resolve()).as_posix()


# --------------------------------------------------------------------------
# 1. Documented usage
# --------------------------------------------------------------------------

@pytest.mark.parametrize("path", NOTEBOOKS, ids=_relative)
def test_notebook_only_uses_the_current_api(path):
    usage = scan_notebook(path)
    problems = []
    for name, cells in sorted(usage.names.items()):
        if not hasattr(ptt, name):
            problems.append(f"cell {cells[0]}: petthermotools has no attribute {name!r}")
    for call in usage.calls:
        if hasattr(ptt, call.name):
            problems += check_call(getattr(ptt, call.name), call)
    assert not problems, f"{_relative(path)} uses the package in a way that no longer works:\n  " + "\n  ".join(problems)


def test_the_scan_really_covers_the_maintained_notebooks():
    """Guards against the scan silently checking nothing (a broken toctree parser,
    a moved folder, notebooks that stopped parsing)."""
    assert len(toctree_notebooks()) >= 15, "the toctree parser found almost no published notebooks"
    assert len([p for p in NOTEBOOKS if "/Workshops/" in p.as_posix()]) >= 8
    assert len([p for p in NOTEBOOKS if "/teaching_materials/" in p.as_posix()]) >= 3
    assert len([p for p in NOTEBOOKS if "/PaperFigures/" in p.as_posix()]) >= 4

    parsed = skipped = calls = 0
    names = set()
    for path in NOTEBOOKS:
        usage = scan_notebook(path)
        parsed += usage.cells_parsed
        skipped += len(usage.cells_skipped)
        calls += len(usage.calls)
        names |= set(usage.names)
    assert parsed >= 500, f"only {parsed} code cells could be parsed"
    assert skipped <= 10, f"{skipped} cells could not be parsed (exotic syntax?); the scan is missing them"
    assert len(names) >= 20 and calls >= 100, f"scan found {len(names)} names and {calls} calls: too few"


def test_no_maintained_notebook_is_excluded_from_the_docs_build():
    """A maintained notebook listed in conf.py exclude_patterns would silently stop being published."""
    excluded = [_relative(p) for p in NOTEBOOKS if (REPO / "docs").resolve() in p.parents and is_excluded_from_docs(p)]
    assert not excluded, f"maintained notebooks excluded from the docs build: {excluded}"


def _fn(source):
    namespace = {}
    exec(source, namespace)
    return namespace["f"]


def _call(n_positional=0, keywords=(), star=False, double_star=False):
    return CallUse("f", 0, n_positional, tuple(keywords), star, double_star)


def test_check_call_rules_on_synthetic_functions():
    f = _fn("def f(a, b=1, *, c=2): pass")
    assert check_call(f, _call(n_positional=2, keywords=["c"])) == []
    assert check_call(f, _call(keywords=["a", "b", "c"])) == []
    assert any("no longer accepts" in p for p in check_call(f, _call(keywords=["zzz"])))
    assert any("positional arguments" in p for p in check_call(f, _call(n_positional=3)))
    # a **kwargs / *args catch-all accepts anything
    assert check_call(_fn("def f(a, **kw): pass"), _call(keywords=["anything"])) == []
    assert check_call(_fn("def f(*args): pass"), _call(n_positional=9)) == []
    # a call that unpacks its arguments cannot be judged on positional count
    assert check_call(f, _call(n_positional=5, star=True)) == []
    # omitting required arguments is not judged (a notebook call is not always complete)
    assert check_call(_fn("def f(a, b): pass"), _call()) == []


# --------------------------------------------------------------------------
# 2. The snapshot
# --------------------------------------------------------------------------

def _entry(source):
    f = _fn(source)
    return {"kind": "function", "package": "petthermotools", "origin": "own",
            "parameters": [encode_param(p) for p in inspect.signature(f).parameters.values()]}


@pytest.mark.parametrize("before, after, expected", [
    ("def f(a, b=1): pass", "def f(a, b=1, c=2): pass", None),                 # optional parameter added at the end
    ("def f(a, b=1): pass", "def f(a, b=1, **kw): pass", None),                # catch-all added
    ("def f(a, b=1): pass", "def f(a): pass", "'b' was removed or renamed"),
    ("def f(a, b=1): pass", "def f(a, c=1): pass", "'b' was removed or renamed"),
    ("def f(a, b=1): pass", "def f(a, b): pass", "used to be optional"),
    ("def f(a, b=1): pass", "def f(a, b=1, c=None, d=3, *, e): pass", "new required parameter 'e'"),
    ("def f(a, b=1): pass", "def f(b=1, a=None): pass", "order of the leading positional"),
    ("def f(a, b=1): pass", "def f(a, *, b=1): pass", "can no longer be passed both"),
    ("def f(a, *, b=1): pass", "def f(a, b=1): pass", None),                    # keyword-only relaxed to normal: fine
    ("def f(a, *, b=1): pass", "def f(a): pass", "'b' was removed or renamed"),
    ("def f(*args): pass", "def f(a): pass", "no longer accepts *args"),
])
def test_snapshot_compatibility_rules(before, after, expected):
    problems = compare_surface({"f": _entry(before)}, {"f": _entry(after)})
    if expected is None:
        assert problems == []
    else:
        assert any(expected in p for p in problems), problems


def test_snapshot_rules_for_names_and_kinds():
    entry = {"kind": "function", "package": "petthermotools", "origin": "own", "parameters": ["a"]}
    assert any("no longer available" in p for p in compare_surface({"f": entry}, {}))
    now_module = {"f": {"kind": "module", "package": "petthermotools", "origin": "own"}}
    assert any("is now a module" in p for p in compare_surface({"f": entry}, now_module))
    # function <-> class is the same to a caller
    as_class = {"f": {"kind": "class", "package": "petthermotools", "origin": "own", "parameters": ["a"]}}
    assert compare_surface({"f": entry}, as_class) == []
    # a name we define that is silently replaced by another library's object is reported
    replaced = {"f": {"kind": "function", "package": "numpy", "origin": "external"}}
    assert any("now comes from 'numpy'" in p for p in compare_surface({"f": entry}, replaced))
    # names from optional packages, and environment-dependent names, are never required
    assert compare_surface({"MELTSdynamic": {"kind": "function", "package": "meltsdynamic", "origin": "external"}}, {}) == []
    assert compare_surface({"pipes": {"kind": "other", "package": "builtins", "origin": "own", "type": "NoneType"}},
                           {"pipes": {"kind": "function", "package": "wurlitzer", "origin": "external"}}) == []
    assert "meltsdynamic" in OPTIONAL_PACKAGES


def test_leaked_names_are_informational_unless_enforced(monkeypatch):
    leaked = {"np": {"kind": "module", "package": "numpy", "origin": "external"}}
    monkeypatch.setattr(_api_scan, "ENFORCE_EXTERNAL_NAMES", False)
    assert compare_surface(leaked, {}) == []
    monkeypatch.setattr(_api_scan, "ENFORCE_EXTERNAL_NAMES", True)
    assert any("no longer available" in p for p in compare_surface(leaked, {}))


def test_the_snapshot_file_is_well_formed():
    snapshot = load_snapshot()
    assert snapshot["format"] == SNAPSHOT_FORMAT
    names = snapshot["names"]
    assert len(names) >= 100, "the snapshot looks truncated"
    for name, entry in names.items():
        assert entry["kind"] in ("function", "class", "module", "other"), name
        assert entry["origin"] in ("own", "external"), name
    # the calculation entry points users rely on are recorded with their parameters
    for key in ("isobaric_crystallisation", "polybaric_crystallisation_path", "equilibrate_multi",
                "phaseDiagram_calc", "saturation_pressure", "harker"):
        assert key in names and names[key].get("parameters"), f"{key} is not recorded with parameters"


def test_the_public_api_still_honours_the_snapshot():
    problems = compare_surface(load_snapshot()["names"], public_surface(ptt))
    assert not problems, (
        f"{len(problems)} breaking change(s) to the public API:\n  " + "\n  ".join(problems[:40])
        + f"\nIf this is deliberate (after a deprecation period), run Unit_tests/api/capture_public_api.py "
        f"and commit the change to {SNAPSHOT.relative_to(REPO)}."
    )
