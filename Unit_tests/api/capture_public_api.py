"""
(Re)generate Unit_tests/api/public_api.json, the snapshot of the public API that
Unit_tests/unit/test_public_api.py protects.

    python Unit_tests/api/capture_public_api.py

Run it ONLY for a deliberate change to the public surface (a name or parameter is
being added, or, after a deprecation period, removed) and commit the diff: the
diff is the review of exactly what changed for users. Additions are always safe;
a removal or rename here is a breaking change.

The snapshot is written from whatever environment runs this, and names that come
from optional packages (alphaMELTS, wurlitzer, ngibbs, ...) are marked by their
package so a minimal environment (CI) does not require them.
"""
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "Unit_tests" / "unit"))
warnings.simplefilter("ignore")

import petthermotools as ptt  # noqa: E402

from _api_scan import SNAPSHOT, public_surface, write_snapshot  # noqa: E402

if __name__ == "__main__":
    surface = public_surface(ptt)
    write_snapshot(surface, ptt.__version__)
    kinds = {}
    for entry in surface.values():
        kinds[entry["kind"]] = kinds.get(entry["kind"], 0) + 1
    print(f"wrote {SNAPSHOT.relative_to(REPO)}: {len(surface)} public names {kinds}")
