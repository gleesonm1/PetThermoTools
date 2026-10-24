"""
Tier A regression test for Mechanism 1 of the "stray python processes left
in Task Manager/Activity Monitor after a failed calculation" symptom.

Originally written against multi_path_v1's batch drain loop (Path.py, now
multi_path_v1), which only called `p.join()` on the *timeout* exit path, never
on the natural-exit path taken every time a worker finishes or crashes before
the round's timeout -- the common case. On POSIX this leaves the already-
exited child as a zombie: an entry in the process table (visible in
`ps`/Activity Monitor, ~0% CPU, harmless individually) that lingers until
something calls join()/wait() on it. On Windows the equivalent is an unclosed
process handle.

multi_path is v2 now (mp_utils.run_tasks-based), with its own, different
natural-exit path: `_retire()`, called from the `_DONE` handler once a worker
reports `retiring=True`. The same bug was found there during this test's
v1->v2 update (2026-10-21): `_retire()` only called `.join()` when
`process.is_alive()` was still true, which is false by the time `_retire()`
runs for the common case (a worker that already exited on its own) -- fixed
in mp_utils.py alongside this test.

The actual check runs in a separate interpreter (_process_reaping_repro.py),
bounded by killing the whole process tree if it doesn't return in time --
this test used to call multi_path directly inside the pytest process with
nothing to catch a hang beyond CI's own 60-minute job timeout, and kept
hanging there in CI (ubuntu/3.10, never reproduced locally) even after the
dead-monkeypatch fix made the rest of this file fast. Whatever the exact
remaining cause, this bounds the damage to one failed test instead of a
stuck hour-long job. See _process_reaping_repro.py's own docstring for why
it spies on mp_utils._CTX.Process rather than Path.Process/Path.path_multi
(same dead-monkeypatch class of bug already fixed in test_multi_path_
supervisor.py and test_timeout_handling.py).
"""
from pathlib import Path

from _proc_helpers import run_script_capped

_REPRO = Path(__file__).resolve().parent / "_process_reaping_repro.py"

# A round (spawn, run, report) should take at most a few seconds for these
# fast, non-hanging fake workers; generous headroom for a busy CI machine.
_CAP_S = 60


def test_every_worker_process_is_joined_after_natural_exit(tmp_path):
    run = run_script_capped([_REPRO], tmp_path, cap_s=_CAP_S)
    if not run.started:
        raise AssertionError(
            f"harness problem: the script never got as far as calling multi_path.\n"
            f"stdout:\n{run.out}\nstderr:\n{run.err[-2000:]}"
        )
    assert run.finished, (
        f"the script was still running {_CAP_S}s after it started -- "
        f"multi_path hung instead of returning.\nstdout:\n{run.out}\nstderr:\n{run.err[-2000:]}"
    )
    assert run.returncode == 0, f"script raised:\n{run.err[-2000:]}"

    result = next((line.split("=", 1)[1] for line in run.out.splitlines() if line.startswith("RESULT=")), None)
    assert result is not None, f"no RESULT= line in output:\n{run.out}"

    if result == "NO_WORKERS":
        raise AssertionError("no worker processes were created -- test setup is broken")
    if result.startswith("STILL_ALIVE"):
        raise AssertionError(f"worker(s) still running after multi_path returned ({result})")
    if result.startswith("UNJOINED"):
        n, total = result.split(":", 1)[1].split("/")
        raise AssertionError(
            f"{n}/{total} worker process(es) exited without ever having .join() called on "
            "them -- left as zombies / unclosed handles."
        )
    assert result == "OK", f"unexpected RESULT={result}"
