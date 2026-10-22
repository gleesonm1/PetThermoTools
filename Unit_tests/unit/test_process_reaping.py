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
in mp_utils.py alongside this test. This test now spies on the actual class
mp_utils.run_tasks uses to spawn workers (`mp_utils._CTX.Process`), not
`Path.Process`/`Path.path_multi`, which multi_path v2 never touches (see
test_multi_path_supervisor.py and test_timeout_handling.py for the same
dead-monkeypatch class of bug, fixed the same way).
"""
import importlib

import numpy as np

# See test_import_hygiene.py: `import petthermotools.Path` is shadowed by
# pathlib.Path via the package's own `from .Path import *`.
ptt_path = importlib.import_module("petthermotools.Path")
ptt_mp_utils = importlib.import_module("petthermotools.mp_utils")

COMP = {
    'SiO2_Liq': 52.0, 'TiO2_Liq': 2.0, 'Al2O3_Liq': 13.0, 'FeOt_Liq': 9.0,
    'MgO_Liq': 9.0, 'CaO_Liq': 10.0, 'Na2O_Liq': 2.0, 'K2O_Liq': 0.4,
    'P2O5_Liq': 0.2, 'Fe3Fet_Liq': 0.2, 'H2O_Liq': 2.0,
}


# Under spawn, Process.start() pickles the WHOLE Process instance (not just
# target/args/kwargs) to hand to the child, so the tracking subclass has to
# be a real, picklable module-level class -- a class defined inside the test
# function has no qualified name reduction.dump() can resolve and fails with
# "Can't get local object ...". It subclasses mp_utils._CTX.Process as it
# stands at import time (the real spawn-context Process class, e.g.
# SpawnProcess) rather than the platform-default `multiprocessing.Process`
# -- those are not always the same class (Linux's bare `multiprocessing.
# Process` defaults to ForkProcess), and swapping in a fork-flavoured class
# under a spawn-context Queue/pipe setup would be a real, separate bug in
# the test itself, not the thing under test.
class _JoinTrackingProcess(ptt_mp_utils._CTX.Process):
    """A drop-in replacement that records every instance created and
    whether .join() was ever called on it. Only the parent-side wrapper is
    subclassed -- the target/args/kwargs sent to the child across the spawn
    boundary are unaffected, so this does not change what actually runs in
    the worker."""

    created = []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.join_called = False
        _JoinTrackingProcess.created.append(self)

    def join(self, *args, **kwargs):
        self.join_called = True
        return super().join(*args, **kwargs)


def _fake_worker_completes_instantly(q, task, **kwargs):
    """Fake _path_task_v2: reports one step for the row and exits normally
    -- i.e. it always takes the natural-exit branch, never the timeout
    branch, isolating Mechanism 1 from the terminate()/timeout behaviour
    covered separately in test_timeout_handling.py."""
    q.put([task.row_id, 0, {"Conditions": {}, "_dummy_keys": []}])


def test_every_worker_process_is_joined_after_natural_exit():
    _JoinTrackingProcess.created.clear()
    orig_cls = ptt_mp_utils._CTX.Process
    orig_task_fn = ptt_path._path_task_v2
    ptt_mp_utils._CTX.Process = _JoinTrackingProcess
    ptt_path._path_task_v2 = _fake_worker_completes_instantly
    try:
        ptt_path.multi_path(
            Model="MELTSv1.0.2",
            comp=COMP,
            T_C=1200.0,
            P_bar=np.array([1000.0, 2000.0, 3000.0, 4000.0, 5000.0]),
            cores=2,
            timeout=5,
            multi_processing=True,
            Print_suppress=True,
        )
    finally:
        ptt_mp_utils._CTX.Process = orig_cls
        ptt_path._path_task_v2 = orig_task_fn

    assert _JoinTrackingProcess.created, "no worker processes were created -- test setup is broken"

    for proc in _JoinTrackingProcess.created:
        assert not proc.is_alive(), "worker should have exited on its own well within the timeout"

    unjoined = [p.pid for p in _JoinTrackingProcess.created if not p.join_called]
    assert not unjoined, (
        f"{len(unjoined)}/{len(_JoinTrackingProcess.created)} worker process(es) "
        f"exited without ever having .join() called on them (pids: {unjoined}) "
        "-- left as zombies / unclosed handles."
    )
