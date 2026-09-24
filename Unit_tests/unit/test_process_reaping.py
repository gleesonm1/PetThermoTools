"""
Tier A regression test for Mechanism 1 of the "stray python processes left
in Task Manager/Activity Monitor after a failed calculation" symptom.

Path.py's batch drain loop (Path.py:602-633) only calls `p.join()` on the
*timeout* exit path:

    except Empty:
        if not p.is_alive():
            break                          # <-- no p.join() here
        if time.time() - start > timeout:
            p.terminate()
            p.join()                       # <-- only reached on timeout
            ...

The natural-exit branch -- taken every time a worker finishes or crashes
*before* the round's timeout, which is the common case -- never joins the
process. On POSIX this leaves the already-exited child as a zombie: an
entry in the process table (visible in `ps`/Activity Monitor, ~0% CPU,
harmless individually) that lingers until something calls join()/wait() on
it. On Windows the equivalent is an unclosed process handle. Over a long
session with many batches/rounds this accumulates -- a plausible explanation
for stray python processes that "don't go away" after failed calculations.

This test spies on multiprocessing.Process to record whether `.join()` was
actually called on every worker multi_path creates, using a fake worker (see
test_multi_path_supervisor.py) so no alphaMELTS/Julia install is needed.
"""
import importlib
import multiprocessing

import numpy as np
import pytest

# See test_import_hygiene.py: `import petthermotools.Path` is shadowed by
# pathlib.Path via the package's own `from .Path import *`.
ptt_path = importlib.import_module("petthermotools.Path")

COMP = {
    'SiO2_Liq': 52.0, 'TiO2_Liq': 2.0, 'Al2O3_Liq': 13.0, 'FeOt_Liq': 9.0,
    'MgO_Liq': 9.0, 'CaO_Liq': 10.0, 'Na2O_Liq': 2.0, 'K2O_Liq': 0.4,
    'P2O5_Liq': 0.2, 'Fe3Fet_Liq': 0.2, 'H2O_Liq': 2.0,
}


class _JoinTrackingProcess(multiprocessing.Process):
    """A drop-in replacement for multiprocessing.Process that records every
    instance created and whether .join() was ever called on it. Only the
    parent-side wrapper is subclassed -- the target/args/kwargs sent to the
    child across the spawn boundary are unaffected, so this does not change
    what actually runs in the worker."""

    created = []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.join_called = False
        _JoinTrackingProcess.created.append(self)

    def join(self, *args, **kwargs):
        self.join_called = True
        return super().join(*args, **kwargs)


def _fake_worker_completes_instantly(q, index, **kwargs):
    """A well-behaved worker that finishes well within any reasonable
    timeout -- i.e. it always takes the natural-exit branch, never the
    timeout branch, isolating Mechanism 1 from the terminate()/timeout
    behaviour covered separately in test_multi_path_supervisor.py."""
    for i in list(index)[:3]:
        q.put([i, -1, "Start"])
        q.put([i, 0, {"Conditions": {}, "_dummy_keys": []}])


@pytest.mark.xfail(
    raises=AssertionError,  # only the intended failure counts as "expected"
    reason=(
        "Known bug (Mechanism 1): Path.py's batch drain loop only calls "
        "p.join() on the timeout exit path (Path.py:622-623). The natural-"
        "exit path (Path.py:616-617, `if not p.is_alive(): break`) never "
        "joins the process, leaving an unreaped zombie/handle behind for "
        "every worker that finishes or crashes before the round timeout. "
        "See the cleanup plan, Part(ii)#1 (persistent worker pool that "
        "reaps every task regardless of how it exited)."
    ),
    strict=False,
)
def test_every_worker_process_is_joined_after_natural_exit():
    _JoinTrackingProcess.created.clear()
    orig_process = ptt_path.Process
    orig_worker = ptt_path.path_multi
    ptt_path.Process = _JoinTrackingProcess
    ptt_path.path_multi = _fake_worker_completes_instantly
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
        ptt_path.Process = orig_process
        ptt_path.path_multi = orig_worker

    assert _JoinTrackingProcess.created, "no worker processes were created -- test setup is broken"

    for proc in _JoinTrackingProcess.created:
        assert not proc.is_alive(), "worker should have exited on its own well within the timeout"

    unjoined = [p.pid for p in _JoinTrackingProcess.created if not p.join_called]
    assert not unjoined, (
        f"{len(unjoined)}/{len(_JoinTrackingProcess.created)} worker process(es) "
        f"exited without ever having .join() called on them (pids: {unjoined}) "
        "-- left as zombies / unclosed handles."
    )
