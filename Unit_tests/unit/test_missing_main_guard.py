"""
Tier A regression test: a user script with no `if __name__ == "__main__":`
guard must fail fast, not hang.

Found in the wild on Windows: a diagnostic script without the guard made
`multi_path` respawn workers for ~15 minutes until it was killed by hand. It is
a very common way to write a script, and it affects macOS too (and Linux under
`spawn`/`forkserver`, the latter being the default from Python 3.14).

Mechanism: every worker is a fresh interpreter that re-imports the main script.
Unguarded, it re-runs the calculation while still importing; Python refuses
("An attempt has been made to start a new process before the current process
has finished its bootstrapping phase") and the worker dies before reporting
anything. The parent's loop only advances on messages from workers, so it
respawns forever: the same root cause as
test_multi_path_does_not_hang_when_worker_never_reports, with a much more
common trigger. Meanwhile the user sees the same traceback repeat endlessly.

The script (_no_guard_repro.py) uses a fake healthy worker, so no MELTS is
needed, and selects `spawn` explicitly so Linux behaves like Windows/macOS.
The guarded twin proves the hang is caused by the missing guard and not by the
harness.
"""
from pathlib import Path

import pytest
from _proc_helpers import run_script_capped

_SCRIPT = Path(__file__).resolve().parent / "_no_guard_repro.py"

# Seconds the script may keep running once it has started. A round (spawn,
# re-import, fail) takes a few seconds, so this covers several respawns.
_CAP_S = 20


def _run(tmp_path, *extra):
    run = run_script_capped([_SCRIPT, *extra], tmp_path, cap_s=_CAP_S)
    if not run.started:
        pytest.fail(f"harness problem: the script never got as far as its calculation.\n{run.err[-2000:]}")
    return run


def test_script_with_main_guard_completes(tmp_path):
    run = _run(tmp_path, "guarded")
    assert run.finished, f"guarded script still running {_CAP_S}s after it started"
    assert run.returncode == 0, f"guarded script failed:\n{run.err[-2000:]}"
    assert "RETURNED=1" in run.out


@pytest.mark.xfail(
    raises=AssertionError,  # only a genuine hang counts as "expected"
    reason=(
        "Known bug: without an `if __name__ == '__main__':` guard every worker "
        "dies while re-importing the script, and multi_path (which only "
        "advances on worker messages) respawns workers forever. See the "
        "cleanup plan, Part (ii)#1: a supervisor that reports crashed workers "
        "after bounded retries, plus a clear message naming the missing guard."
    ),
    strict=False,
)
def test_script_without_main_guard_fails_fast(tmp_path):
    run = _run(tmp_path)
    assert run.finished, (
        f"the script was still running {_CAP_S}s after it started: workers die "
        f"while bootstrapping ({run.err.count('bootstrapping phase')} "
        "'bootstrapping phase' errors so far) and are respawned indefinitely"
    )
