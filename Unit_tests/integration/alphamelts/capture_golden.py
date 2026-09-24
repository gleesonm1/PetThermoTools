"""
(Re)generate the golden reference outputs used by test_golden_melts.py.

    python Unit_tests/integration/alphamelts/capture_golden.py            # all cases
    python Unit_tests/integration/alphamelts/capture_golden.py polybaric_v102 ...

Run this ONLY when you have decided that a change in results is correct (for
example, a deliberate constant change such as the Fe2O3 molar mass): the
golden files are the definition of "unchanged", so regenerating them and
committing the diff is how a deliberate numerical change is recorded.

Each case is written to Unit_tests/golden/alphamelts/<case>/ as numbered CSV files
(text, so they are independent of the pandas/numpy versions that wrote them)
plus `_frames.json` (file -> result key) and `_meta.json` (what produced it).
Runs from a scratch directory because MELTS drops files in the CWD.
"""
import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import warnings
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
GOLDEN = REPO / "Unit_tests" / "golden" / "alphamelts"

sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(HERE))
warnings.simplefilter("ignore")

import numpy  # noqa: E402
import pandas  # noqa: E402

import _cases  # noqa: E402


def _git_commit():
    try:
        out = subprocess.run(
            ["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() or None
    except Exception:  # noqa: BLE001
        return None


def main(names):
    import petthermotools

    unknown = [n for n in names if n not in _cases.CASES]
    if unknown:
        sys.exit(f"unknown case(s): {unknown}; available: {list(_cases.CASES)}")

    os.chdir(tempfile.mkdtemp(prefix="ptt_golden_"))  # MELTS writes files into the CWD
    meta = {
        "captured_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "pandas": pandas.__version__,
        "numpy": numpy.__version__,
        "petthermotools": getattr(petthermotools, "__version__", None),
        "git_commit": _git_commit(),
    }

    for name in names:
        with _cases.quiet_fds():
            result = _cases.run_case(name)
        frames = _cases.flatten_frames(result)

        case_dir = GOLDEN / name
        shutil.rmtree(case_dir, ignore_errors=True)
        case_dir.mkdir(parents=True)
        index = {}
        for i, (key, df) in enumerate(sorted(frames.items())):
            filename = f"{i:03d}.csv"
            df.to_csv(case_dir / filename, float_format="%.12g")
            index[filename] = key
        (case_dir / "_frames.json").write_text(json.dumps(index, indent=1))
        (case_dir / "_meta.json").write_text(json.dumps(meta, indent=1))
        size_kb = sum(p.stat().st_size for p in case_dir.iterdir()) / 1024
        print(f"captured {name:<34} {len(frames):>3} frames  {size_kb:7.0f} KB", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cases", nargs="*", help="case names (default: all)")
    args = parser.parse_args()
    main(args.cases or list(_cases.CASES))
