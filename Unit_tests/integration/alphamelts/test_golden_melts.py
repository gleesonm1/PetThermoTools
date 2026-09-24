"""
Tier B1 golden-output tests: real alphaMELTS calculations must reproduce the
recorded results.

These are the safety net for the cleanup and speed work. Every refactor of
Path.py / MELTS.py / GenFuncs.py / the wrappers must leave these unchanged;
a deliberate numerical change is recorded by re-running capture_golden.py and
committing the diff.

Needs a working alphaMELTS install (`meltsdynamic`): auto-skipped otherwise
(see Unit_tests/conftest.py); local-only by decision, not run in CI. Each case
takes a few seconds and runs its calculation in worker processes exactly as a
user's call would, from a scratch directory (MELTS drops files in the CWD).

Tolerance: two runs on one machine are bit-for-bit identical (checked for both
the single-run and the batch path). Different platforms (macOS vs Windows) may
differ in the last digits, so comparison uses rtol=1e-6 by default; set
PTT_GOLDEN_RTOL to tighten or loosen it. The set of phases/tables produced
must match exactly.
"""
import json
import os
from pathlib import Path

import pandas as pd
import pytest

from _cases import CASES, flatten_frames, quiet_fds, run_case

pytestmark = pytest.mark.melts

GOLDEN = Path(__file__).resolve().parents[2] / "golden" / "alphamelts"
RTOL = float(os.environ.get("PTT_GOLDEN_RTOL", "1e-6"))
ATOL = 1e-8


def _load_golden(name):
    case_dir = GOLDEN / name
    index_file = case_dir / "_frames.json"
    if not index_file.exists():
        pytest.fail(
            f"no golden data for {name!r} in {case_dir}. Capture it with:\n"
            f"    python Unit_tests/integration/alphamelts/capture_golden.py {name}"
        )
    index = json.loads(index_file.read_text())
    return {key: pd.read_csv(case_dir / filename, index_col=0) for filename, key in index.items()}


@pytest.mark.parametrize("name", list(CASES))
def test_matches_golden(name, tmp_path, monkeypatch):
    expected = _load_golden(name)

    monkeypatch.chdir(tmp_path)  # MELTS writes *_tbl.txt / .inp files into the CWD
    # The engine writes ~420 KB of progress text per run to stderr from several
    # workers at once. Left alone, pytest shows it interleaved character by
    # character in the failure report and buries the assertion (seen on Windows).
    # quiet_fds() silences it, and still shows a worker's traceback if the run raises.
    with quiet_fds():
        actual = flatten_frames(run_case(name))

    assert set(actual) == set(expected), (
        f"tables differ. only in result: {sorted(set(actual) - set(expected))}; "
        f"only in golden: {sorted(set(expected) - set(actual))}"
    )

    for key in sorted(expected):
        got = actual[key].copy()
        got.columns = got.columns.astype(str)  # CSV round-trips column labels as text
        try:
            pd.testing.assert_frame_equal(
                got, expected[key],
                check_dtype=False, check_index_type=False, check_column_type=False,
                check_names=False, rtol=RTOL, atol=ATOL,
            )
        except AssertionError as exc:
            raise AssertionError(f"{name}: table {key!r} differs from golden\n{exc}") from None
