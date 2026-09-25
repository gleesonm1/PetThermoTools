"""
Helpers for the public-API freeze (see test_public_api.py and
Unit_tests/api/capture_public_api.py). Not a test module.

Two things are protected:

1. Documented usage: which package names, and which arguments, the MAINTAINED
   notebooks actually use (read straight from the notebooks).
2. The snapshot: every name that is public today (`dir(petthermotools)` without a
   leading underscore) and, for functions and classes defined in the package,
   their parameters.

"Maintained" is not a hand-kept list: it is what ReadTheDocs publishes (the
notebooks reachable from the toctrees in docs/index.rst), plus everything in
Workshops/, docs/teaching_materials/ and docs/PaperFigures/. Notebooks listed in
docs/conf.py `exclude_patterns` are retired and not checked.
"""
import ast
import fnmatch
import inspect
import json
import re
import warnings
from pathlib import Path
from typing import NamedTuple

REPO = Path(__file__).resolve().parents[2]
DOCS = REPO / "docs"
SNAPSHOT = REPO / "Unit_tests" / "api" / "public_api.json"
SNAPSHOT_FORMAT = 1

# Names that only exist (or only have their real type) when an optional package is
# installed. They are recorded in the snapshot but never required, so a minimal
# environment (CI) does not fail on them.
OPTIONAL_PACKAGES = {"meltsdynamic", "wurlitzer", "ngibbs", "juliacall", "juliapkg",
                     "torch", "IPython", "ipywidgets", "pyMelt"}

# Names that always exist but whose TYPE depends on the environment: `pipes` (Path.py)
# is a wurlitzer function inside a Jupyter kernel and None everywhere else. Only their
# existence is required.
ENVIRONMENT_DEPENDENT_NAMES = {"pipes"}

# Names that merely leak into the namespace from other libraries (`ptt.np`, `ptt.Path`,
# `ptt.Process`, ... via star-imports) are recorded for information but not enforced:
# they are not PetThermoTools functions, and requiring them would lock the clutter in.
# A leaked name that a maintained notebook uses is still protected by the
# documented-usage tier. Set to True to enforce them as well.
ENFORCE_EXTERNAL_NAMES = False


# ---------------------------------------------------------------------------
# Which notebooks are maintained
# ---------------------------------------------------------------------------

def _toctree_entries(rst):
    """(entry, glob) pairs from every `.. toctree::` block in an .rst file."""
    lines = rst.read_text(encoding="utf-8", errors="replace").split("\n")
    out, i = [], 0
    while i < len(lines):
        m = re.match(r"^(\s*)\.\. toctree::", lines[i])
        if not m:
            i += 1
            continue
        base, glob = len(m.group(1)), False
        i += 1
        while i < len(lines):
            line = lines[i]
            if line.strip() == "":
                i += 1
                continue
            if len(line) - len(line.lstrip()) <= base:
                break
            text = line.strip()
            if text.startswith(":"):
                glob = glob or text.startswith(":glob:")
            elif not text.startswith(".."):  # `.. entry` is a commented-out entry
                titled = re.match(r".*<(.+)>$", text)
                out.append((titled.group(1).strip() if titled else text, glob))
            i += 1
    return out


def _resolve(rst, entry, glob):
    base = DOCS if entry.startswith("/") else rst.parent
    entry = entry.lstrip("/")
    if glob or any(c in entry for c in "*?["):
        hits = {p for pattern in (entry, entry + "*") for p in base.glob(pattern)}
        return sorted(p for p in hits if p.suffix in (".rst", ".ipynb"))
    for ext in ("", ".rst", ".ipynb"):
        p = (base / (entry + ext)).resolve()
        if p.is_file() and p.suffix in (".rst", ".ipynb"):
            return [p]
    return []


def toctree_notebooks():
    """Notebooks reachable from docs/index.rst through the toctrees: what ReadTheDocs links."""
    seen, notebooks, queue = set(), set(), [DOCS / "index.rst"]
    while queue:
        rst = queue.pop()
        if rst in seen or not rst.is_file():
            continue
        seen.add(rst)
        for entry, glob in _toctree_entries(rst):
            for p in _resolve(rst, entry, glob):
                if p.suffix == ".ipynb":
                    notebooks.add(p.resolve())
                else:
                    queue.append(p)
    return notebooks


def maintained_notebooks():
    """The notebooks whose use of the package is protected, sorted."""
    notebooks = set(toctree_notebooks())
    for folder in ("Workshops", "docs/teaching_materials", "docs/PaperFigures"):
        notebooks |= {p.resolve() for p in (REPO / folder).rglob("*.ipynb")
                      if ".ipynb_checkpoints" not in p.parts}
    return sorted(notebooks)


def conf_exclude_patterns():
    """`exclude_patterns` from docs/conf.py, read without executing the config."""
    tree = ast.parse((DOCS / "conf.py").read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "exclude_patterns" for t in node.targets):
            return ast.literal_eval(node.value)
    return []


def is_excluded_from_docs(path):
    rel = path.resolve().relative_to(DOCS.resolve()).as_posix()
    return any(fnmatch.fnmatch(rel, pattern) for pattern in conf_exclude_patterns())


# ---------------------------------------------------------------------------
# What a notebook uses
# ---------------------------------------------------------------------------

class CallUse(NamedTuple):
    name: str
    cell: int
    n_positional: int
    keywords: tuple
    has_star: bool      # f(*args)
    has_double_star: bool  # f(**kwargs)


class Usage:
    def __init__(self):
        self.names = {}          # package attribute -> [cells that use it]
        self.calls = []          # CallUse
        self.cells_parsed = 0
        self.cells_skipped = []  # cells that could not be parsed (exotic syntax)


def _sanitize(source):
    """Replace IPython magics (%, !, trailing ?) with `pass` so the cell parses."""
    out = []
    for line in source.split("\n"):
        stripped = line.lstrip()
        is_magic = stripped.startswith(("%", "!")) or (stripped.rstrip().endswith("?") and not stripped.startswith("#"))
        out.append(line[:len(line) - len(stripped)] + "pass" if is_magic else line)
    return "\n".join(out)


def scan_notebook(path):
    """Which package names, and which call arguments, a notebook uses.

    The alias (`import petthermotools as ptt`) is collected over the whole notebook
    first, because the import is normally in an earlier cell than its uses."""
    usage = Usage()
    data = json.loads(path.read_text(encoding="utf-8"))
    trees = []
    for index, cell in enumerate(data.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        source = cell.get("source", "")
        source = "".join(source) if isinstance(source, list) else source
        try:
            with warnings.catch_warnings():
                # notebook code often has strings like '\degree'; that warning is the notebook's, not ours
                warnings.simplefilter("ignore")
                trees.append((index, ast.parse(_sanitize(source))))
        except SyntaxError:
            usage.cells_skipped.append(index)
            continue
        usage.cells_parsed += 1

    aliases, direct = set(), {}
    for _, tree in trees:
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    if a.name == "petthermotools":
                        aliases.add(a.asname or "petthermotools")
            elif isinstance(node, ast.ImportFrom) and node.module == "petthermotools":
                for a in node.names:
                    if a.name != "*":
                        direct[a.asname or a.name] = a.name

    for index, tree in trees:
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id in aliases:
                usage.names.setdefault(node.attr, []).append(index)
            if isinstance(node, ast.ImportFrom) and node.module == "petthermotools":
                for a in node.names:
                    if a.name != "*":
                        usage.names.setdefault(a.name, []).append(index)
            if isinstance(node, ast.Call):
                func, name = node.func, None
                if isinstance(func, ast.Attribute) and isinstance(func.value, ast.Name) and func.value.id in aliases:
                    name = func.attr
                elif isinstance(func, ast.Name) and func.id in direct:
                    name = direct[func.id]
                if name:
                    usage.calls.append(CallUse(
                        name=name, cell=index,
                        n_positional=sum(1 for a in node.args if not isinstance(a, ast.Starred)),
                        keywords=tuple(k.arg for k in node.keywords if k.arg),
                        has_star=any(isinstance(a, ast.Starred) for a in node.args),
                        has_double_star=any(k.arg is None for k in node.keywords)))
    return usage


def check_call(obj, call):
    """Problems (strings) if `call` could not be bound to `obj`'s current signature.
    Required-argument omissions are NOT checked: a notebook call is not always complete."""
    try:
        params = inspect.signature(obj).parameters
    except (TypeError, ValueError):
        return []
    P = inspect.Parameter
    var_positional = any(p.kind is P.VAR_POSITIONAL for p in params.values())
    var_keyword = any(p.kind is P.VAR_KEYWORD for p in params.values())
    positional = [p for p in params.values() if p.kind in (P.POSITIONAL_ONLY, P.POSITIONAL_OR_KEYWORD)]
    problems = []
    if not call.has_star and not var_positional and call.n_positional > len(positional):
        problems.append(f"cell {call.cell}: {call.name}() is called with {call.n_positional} positional arguments "
                        f"but now accepts at most {len(positional)}")
    if not var_keyword:
        for keyword in call.keywords:
            p = params.get(keyword)
            if p is None or p.kind in (P.POSITIONAL_ONLY, P.VAR_POSITIONAL, P.VAR_KEYWORD):
                problems.append(f"cell {call.cell}: {call.name}() is called with {keyword}=... "
                                f"but it no longer accepts that argument")
    return problems


# ---------------------------------------------------------------------------
# The snapshot of the public surface
# ---------------------------------------------------------------------------

def encode_param(p):
    """One short string per parameter: 'x' required, 'x=' has a default,
    'kw:x' keyword-only, 'pos:x' positional-only, '*args', '**kwargs'."""
    if p.kind is p.VAR_POSITIONAL:
        return "*" + p.name
    if p.kind is p.VAR_KEYWORD:
        return "**" + p.name
    prefix = {p.POSITIONAL_ONLY: "pos:", p.KEYWORD_ONLY: "kw:"}.get(p.kind, "")
    return prefix + p.name + ("=" if p.default is not p.empty else "")


def decode_param(text):
    """-> (name, kind, has_default) with kind in pos, pos_or_kw, kw, var_pos, var_kw."""
    if text.startswith("**"):
        return text[2:], "var_kw", False
    if text.startswith("*"):
        return text[1:], "var_pos", False
    kind = "pos_or_kw"
    if text.startswith("pos:"):
        kind, text = "pos", text[4:]
    elif text.startswith("kw:"):
        kind, text = "kw", text[3:]
    has_default = text.endswith("=")
    return text.rstrip("="), kind, has_default


def _kind_of(obj):
    if inspect.ismodule(obj):
        return "module"
    if inspect.isclass(obj):
        return "class"
    return "function" if callable(obj) else "other"


def _package_of(obj, kind):
    name = obj.__name__ if kind == "module" else (getattr(obj, "__module__", None) or type(obj).__module__)
    return (name or "builtins").split(".")[0]


def public_surface(ptt):
    """{name: entry} for every public-looking name of the package, in this environment."""
    surface = {}
    for name in sorted(n for n in dir(ptt) if not n.startswith("_")):
        obj = getattr(ptt, name)
        kind = _kind_of(obj)
        package = _package_of(obj, kind)
        entry = {"kind": kind, "package": package,
                 "origin": "own" if package in ("petthermotools", "builtins") else "external"}
        if kind == "other":
            entry["type"] = type(obj).__name__
        if kind in ("function", "class") and package == "petthermotools":
            try:
                entry["parameters"] = [encode_param(p) for p in inspect.signature(obj).parameters.values()]
            except (TypeError, ValueError):
                pass
        surface[name] = entry
    return surface


def write_snapshot(surface, version, path=SNAPSHOT):
    lines = ["{", f'"format": {SNAPSHOT_FORMAT},', f'"petthermotools_version": {json.dumps(version)},', '"names": {']
    items = list(surface.items())
    for i, (name, entry) in enumerate(items):
        lines.append(f"{json.dumps(name)}: {json.dumps(entry)}{',' if i < len(items) - 1 else ''}")
    lines += ["}", "}"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_snapshot(path=SNAPSHOT):
    return json.loads(path.read_text(encoding="utf-8"))


def _compare_parameters(name, recorded, current):
    problems = []
    rec = [decode_param(t) for t in recorded]
    cur = [decode_param(t) for t in current]
    cur_by_name = {n: (k, d) for n, k, d in cur}

    for pname, kind, has_default in rec:
        if kind in ("var_pos", "var_kw"):
            if not any(k == kind for _, k, _ in cur):
                problems.append(f"{name}(): no longer accepts {'*' if kind == 'var_pos' else '**'}{pname}")
            continue
        if pname not in cur_by_name:
            problems.append(f"{name}(): parameter '{pname}' was removed or renamed")
            continue
        now_kind, now_default = cur_by_name[pname]
        if kind == "pos_or_kw" and now_kind != "pos_or_kw":
            problems.append(f"{name}(): parameter '{pname}' can no longer be passed both by position and by keyword")
        elif kind == "pos" and now_kind not in ("pos", "pos_or_kw"):
            problems.append(f"{name}(): positional-only parameter '{pname}' can no longer be passed by position")
        elif kind == "kw" and now_kind not in ("kw", "pos_or_kw"):
            problems.append(f"{name}(): keyword-only parameter '{pname}' can no longer be passed by keyword")
        if has_default and not now_default:
            problems.append(f"{name}(): parameter '{pname}' used to be optional and is now required")

    rec_pos = [n for n, k, _ in rec if k in ("pos", "pos_or_kw")]
    cur_pos = [n for n, k, _ in cur if k in ("pos", "pos_or_kw")]
    if cur_pos[:len(rec_pos)] != rec_pos and all(n in cur_by_name for n in rec_pos):
        problems.append(f"{name}(): the order of the leading positional parameters changed "
                        f"({', '.join(rec_pos[:6])}{', ...' if len(rec_pos) > 6 else ''}); callers passing them by position break")

    recorded_names = {n for n, _, _ in rec}
    for pname, kind, has_default in cur:
        if kind not in ("var_pos", "var_kw") and not has_default and pname not in recorded_names:
            problems.append(f"{name}(): new required parameter '{pname}' (existing callers do not pass it)")
    return problems


def compare_surface(recorded_names, current_surface):
    """Problems (strings) where the current public surface no longer honours the snapshot.
    Additions are fine: new names, new optional parameters at the end."""
    problems = []
    for name, entry in recorded_names.items():
        optional = entry.get("package") in OPTIONAL_PACKAGES
        if entry.get("origin") == "external" and not ENFORCE_EXTERNAL_NAMES:
            continue
        if name not in current_surface:
            if not optional:
                problems.append(f"{name}: no longer available as petthermotools.{name}")
            continue
        if optional or name in ENVIRONMENT_DEPENDENT_NAMES:
            continue
        now = current_surface[name]
        callable_kinds = ("function", "class")
        same_kind = entry["kind"] == now["kind"] or (entry["kind"] in callable_kinds and now["kind"] in callable_kinds)
        if not same_kind:
            problems.append(f"{name}: was a {entry['kind']}, is now a {now['kind']}")
            continue
        if entry["origin"] == "own" and now["origin"] == "external":
            # e.g. replaced by a function from another library: its parameters can no longer be checked
            problems.append(f"{name}: used to be defined in petthermotools, now comes from '{now['package']}'")
            continue
        if entry.get("parameters") is not None and now.get("parameters") is not None:
            problems += _compare_parameters(name, entry["parameters"], now["parameters"])
    return problems
