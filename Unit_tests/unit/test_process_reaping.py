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

multi_path is v2 now (mp_utils.run_tasks-based). Its supervisor calls
`is_alive()` on every worker it retires, and CPython's `is_alive()` reaps an
exited child (waitpid on POSIX) -- so v2 does not have the v1 bug. This test
checks the outcome (no worker left unreaped after multi_path returns), not
whether `.join()` specifically was called.

The actual check runs in a separate interpreter (_process_reaping_repro.py),
bounded by killing the whole process tree if it doesn't return in time. An
earlier version of the repro swapped a Process subclass into
mp_utils._CTX.Process and hung on ubuntu/Python 3.10 CI only; see the repro's
docstring for the replacement approach and the faulthandler stack dumps it
now prints if it ever hangs again.
"""
from pathlib import Path

from _proc_helpers import run_script_capped

_REPRO = Path(__file__).resolve().parent / "_process_reaping_repro.py"

# A round (spawn, run, report) should take at most a few seconds for these
# fast, non-hanging fake workers; generous headroom for a busy CI machine.
_CAP_S = 60


def test_every_worker_process_is_reaped_after_natural_exit(tmp_path):
    run = run_script_capped([_REPRO, tmp_path / "worker_pids.txt"], tmp_path, cap_s=_CAP_S)
    if not run.started:
        raise AssertionError(
            f"harness problem: the script never got as far as calling multi_path.\n"
            f"stdout:\n{run.out}\nstderr:\n{run.err[-2000:]}"
        )
    # stderr holds faulthandler stack dumps of every still-running process
    # (see _process_reaping_repro.py), so show plenty of it on a hang.
    assert run.finished, (
        f"the script was still running {_CAP_S}s after it started -- "
        f"multi_path hung instead of returning.\nstdout:\n{run.out}\nstderr:\n{run.err[-12000:]}"
    )
    assert run.returncode == 0, f"script raised:\n{run.err[-2000:]}"

    result = next((line.split("=", 1)[1] for line in run.out.splitlines() if line.startswith("RESULT=")), None)
    assert result is not None, f"no RESULT= line in output:\n{run.out}"

    if result == "NO_WORKERS":
        raise AssertionError("no worker processes were created -- test setup is broken")
    if result.startswith("STILL_ALIVE"):
        raise AssertionError(f"worker(s) still running after multi_path returned ({result})")
    if result.startswith("UNJOINED"):
        n = result.split(":", 1)[1]
        raise AssertionError(
            f"{n} worker process(es) exited without ever having .join() called on "
            "them -- left as zombies / unclosed handles."
        )
    assert result == "OK", f"unexpected RESULT={result}"
