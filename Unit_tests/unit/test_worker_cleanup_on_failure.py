"""
Tier A regression test: worker processes must not outlive a failed parent.

Path.py's wait loops (Path.py:416-446 single run, 602-633 batches) have no
try/finally. The only `except` catches `Empty`, so a KeyboardInterrupt
(stopping a notebook cell) or any other exception in the parent propagates
straight out without terminating or joining the workers it started. Those
workers are non-daemon Process objects whose parent -- the notebook kernel --
stays alive, so they keep running. Worse, the timeout is enforced *by the
parent's loop*, so once that loop is gone nothing ever stops a hung worker:
it stays in Task Manager / Activity Monitor, holding memory and CPU, until
the kernel is restarted.

Each scenario runs in its own interpreter (_orphan_repro.py); the test kills
any survivors it finds so the test run itself never leaves processes behind.

Harness/setup problems use pytest.fail() (not AssertionError) so the
xfail(raises=AssertionError) below can only ever be satisfied by the real
bug, never by a broken harness.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest
from _proc_helpers import kill_process_tree, utf8_env

_REPRO = Path(__file__).resolve().parent / "_orphan_repro.py"

_EXPECTED_EXCEPTION = {"interrupt": "KeyboardInterrupt", "raises": "RuntimeError"}


def _run_orphan_repro(mode, timeout_s=90):
    """Returns (caught, alive_pids) reported by the repro script.

    Output goes to temp files, not pipes: the surviving workers inherit the
    script's stdout, so with a pipe `communicate()` would wait for EOF until
    those very workers died -- exactly the condition being tested for.
    """
    utf8 = {"encoding": "utf-8", "errors": "replace"}  # see _proc_helpers.utf8_env
    with tempfile.TemporaryFile("w+", **utf8) as out_f, tempfile.TemporaryFile("w+", **utf8) as err_f:
        proc = subprocess.Popen(
            [sys.executable, str(_REPRO), mode], stdout=out_f, stderr=err_f, env=utf8_env()
        )
        alive = []
        try:
            try:
                proc.wait(timeout=timeout_s)
            except subprocess.TimeoutExpired:
                pytest.fail(f"_orphan_repro.py {mode} did not finish within {timeout_s}s")

            out_f.seek(0)
            err_f.seek(0)
            out, err = out_f.read(), err_f.read()

            lines = dict(line.split("=", 1) for line in out.splitlines() if "=" in line)
            if "CAUGHT" not in lines or "ALIVE" not in lines:
                pytest.fail(f"_orphan_repro.py {mode} produced no report.\nstdout:\n{out}\nstderr:\n{err}")

            alive = [int(pid) for pid in lines["ALIVE"].split(",") if pid]
            return lines["CAUGHT"], alive
        finally:
            # Whatever the outcome, don't leave the workers running.
            for pid in alive:
                kill_process_tree(pid)
            kill_process_tree(proc.pid)


@pytest.mark.parametrize("mode", ["interrupt", "raises"])
@pytest.mark.xfail(
    raises=AssertionError,  # only the intended failure counts as "expected"
    reason=(
        "Known bug: multi_path's wait loops have no try/finally, so if the "
        "parent is interrupted or raises mid-run, the workers it started are "
        "never terminated or joined -- and with the parent's loop gone, "
        "nothing enforces their timeout. See the cleanup plan, Part(ii)#1 "
        "(supervisor that reaps every worker in a finally block)."
    ),
    strict=False,
)
def test_workers_are_cleaned_up_when_parent_fails_mid_run(mode):
    caught, alive = _run_orphan_repro(mode)

    if caught != _EXPECTED_EXCEPTION[mode]:
        pytest.fail(
            f"harness problem: expected the parent to fail with "
            f"{_EXPECTED_EXCEPTION[mode]} but it reported {caught!r}"
        )

    assert not alive, (
        f"{len(alive)} worker process(es) were still running after the parent "
        f"failed with {caught} (pids: {alive}) -- unsupervised, with no timeout "
        "enforcement."
    )
