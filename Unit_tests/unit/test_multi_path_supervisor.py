"""
Tier A characterization tests for the multi_path batch supervisor in Path.py.

These replace the real MELTS/MAGEMin worker (`path_multi`) with a fake so
they can exercise the *parent*-process orchestration logic (queue draining,
timeout handling, round retries) without alphaMELTS or Julia installed.
`Process(target=path_multi, ...)` resolves `path_multi` as a module-global
name at call time, so patching `petthermotools.Path.path_multi` redirects
every subsequently-spawned worker to the fake.

The motivating failure mode is the one reported by the maintainer: the
alphaMELTS C library cannot be re-initialized in a process after a failed
calculation. In the current code that means a fresh worker process can die
*before ever writing anything to its queue* (its `MELTSdynamic(...)`
construction at Path.py:792-802 has no try/except around it).

Each scenario runs in a *separate interpreter* (see _multi_path_repro.py) and
the whole process tree is killed if it exceeds its time bound. An earlier
version ran multi_path on a daemon thread inside the pytest process; when
the hang was reproduced, that thread could not be killed and kept
respawning workers into later tests (and, once the patch was undone, ran
the real path_multi in them).
"""
import subprocess
import sys
from pathlib import Path

import pytest
from _proc_helpers import kill_process_tree, utf8_env

_REPRO = Path(__file__).resolve().parent / "_multi_path_repro.py"


def _run_repro(mode, timeout_s):
    """Run `_multi_path_repro.py <mode>` in a fresh interpreter.

    Returns (finished, returncode, stderr). If it has not exited within
    timeout_s, the whole process tree (parent + every spawned worker) is
    killed and finished=False is returned.
    """
    proc = subprocess.Popen(
        [sys.executable, str(_REPRO), mode],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
        errors="replace",
        env=utf8_env(),
    )
    try:
        _, err = proc.communicate(timeout=timeout_s)
        return True, proc.returncode, err
    except subprocess.TimeoutExpired:
        kill_process_tree(proc.pid)
        _, err = proc.communicate(timeout=30)
        return False, None, err


@pytest.mark.xfail(
    raises=AssertionError,  # only the intended failure counts as "expected"
    reason=(
        "Known bug: Path.py's batch loop (Path.py:572-666) only advances "
        "index_out from messages placed on the worker's queue. A worker "
        "whose engine construction fails before its first q.put -- e.g. the "
        "alphaMELTS C library refusing to reinitialize after a prior crash "
        "-- is silently respawned forever. See the cleanup plan, Part(ii)#1 "
        "(persistent worker pool with bounded retries)."
    ),
    strict=False,
)
def test_multi_path_does_not_hang_when_worker_never_reports():
    finished, _, _ = _run_repro("never_reports", timeout_s=30)

    assert finished, (
        "multi_path did not return within 30s when every worker fails before "
        "its first queue message -- this reproduces the alphaMELTS "
        "reinitialization failure as an infinite respawn loop."
    )


def test_multi_path_completes_when_worker_reports_promptly():
    """Sanity check on the harness itself: a well-behaved fake worker must
    NOT also appear to hang, so the xfail above is pinned on the "never
    reports" behaviour and not on some unrelated harness bug."""
    finished, returncode, err = _run_repro("completes", timeout_s=120)

    assert finished, "multi_path hung even though every worker reported promptly"
    assert returncode == 0, f"multi_path raised in the repro script:\n{err}"
