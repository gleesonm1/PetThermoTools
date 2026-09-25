"""
(Re)generate Unit_tests/api/public_api.json, the snapshot of the public API that
Unit_tests/unit/test_public_api.py protects.

    python Unit_tests/api/capture_public_api.py                 # from the COMMITTED code (HEAD)
    python Unit_tests/api/capture_public_api.py --working-tree  # from your working tree (rarely right)

Run it ONLY for a deliberate change to the public surface (a name or parameter is
being added, or, after a deprecation period, removed) and commit the diff: the
diff is the review of exactly what changed for users. Additions are always safe;
a removal or rename here is a breaking change.

Why HEAD by default: CI and every other user see the committed code. A snapshot
taken from a working tree with uncommitted work records names that do not exist
there (this happened: an uncommitted `trace_engine` import put five names into
the snapshot and CI failed on every run). So the package is imported from a
clean checkout of HEAD in a temporary git worktree; uncommitted changes under
src/ are NOT included. Commit the code change first, then capture.

Names that come from optional packages (alphaMELTS, wurlitzer, ngibbs, ...) are
marked by their package, so a minimal environment (CI) does not require them.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
UNIT = REPO / "Unit_tests" / "unit"
sys.path.insert(0, str(UNIT))

from _api_scan import SNAPSHOT, write_snapshot  # noqa: E402

# Runs in a separate interpreter so the package is imported from the chosen source tree.
_PROBE = r"""
import json, sys, warnings
warnings.simplefilter("ignore")
sys.path.insert(0, sys.argv[1])   # <source tree>/src
sys.path.insert(0, sys.argv[2])   # Unit_tests/unit (for _api_scan)
import petthermotools as ptt
from _api_scan import public_surface
print("RESULT=" + json.dumps({"version": ptt.__version__, "surface": public_surface(ptt)}))
"""


def _surface_from(source_root):
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONIOENCODING": "utf-8"}
    run = subprocess.run(
        [sys.executable, "-c", _PROBE, str(source_root / "src"), str(UNIT)],
        capture_output=True, text=True, encoding="utf-8", errors="replace", env=env, cwd=source_root,
    )
    lines = [ln for ln in run.stdout.splitlines() if ln.startswith("RESULT=")]
    if run.returncode != 0 or not lines:
        sys.exit(f"could not import the package from {source_root}:\n{run.stderr[-1500:]}")
    result = json.loads(lines[-1][len("RESULT="):])
    return result["surface"], result["version"]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--working-tree", action="store_true",
                        help="import the package from the working tree instead of the committed code")
    args = parser.parse_args()

    if args.working_tree:
        surface, version = _surface_from(REPO)
        source = "the WORKING TREE (uncommitted changes included)"
    else:
        head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "--", "src"],
                               capture_output=True, text=True, check=True).stdout.strip()
        with tempfile.TemporaryDirectory(prefix="ptt_api_head_") as tmp:
            checkout = Path(tmp) / "checkout"
            subprocess.run(["git", "-C", str(REPO), "worktree", "add", "--detach", "--quiet", str(checkout), "HEAD"],
                           check=True)
            try:
                surface, version = _surface_from(checkout)
            finally:
                subprocess.run(["git", "-C", str(REPO), "worktree", "remove", "--force", str(checkout)],
                               capture_output=True)
        source = f"the committed code at {head}"
        if dirty:
            print(f"note: uncommitted changes under src/ are NOT included:\n" +
                  "\n".join("   " + ln for ln in dirty.splitlines()[:8]))

    write_snapshot(surface, version)
    kinds = {}
    for entry in surface.values():
        kinds[entry["kind"]] = kinds.get(entry["kind"], 0) + 1
    print(f"wrote {SNAPSHOT.relative_to(REPO)} from {source}: {len(surface)} public names {kinds}")


if __name__ == "__main__":
    main()
