"""
Helper script for test_worker_cleanup_on_failure.py -- NOT a test module.

Simulates what a notebook user experiences when a cell running multi_path is
stopped (KeyboardInterrupt) or the parent raises partway through: the
exception is caught by the caller (the "kernel" lives on), and we then report
which worker processes are still running.

Usage:  python _orphan_repro.py interrupt | raises

Prints two lines for the test to parse:
    CAUGHT=<exception class name, or None>
    ALIVE=<comma-separated PIDs of workers still running, possibly empty>

Exits with os._exit so Python's exit-time join of non-daemon children can't
hang the script; the test kills whatever ALIVE reports.
"""
import ctypes
import importlib
import multiprocessing
import os
import sys
import threading
import time
import _thread
from pathlib import Path

import numpy as np

# Test the working tree, not a stale installed copy (mirrors conftest.py).
_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

COMP = {
    'SiO2_Liq': 52.0, 'TiO2_Liq': 2.0, 'Al2O3_Liq': 13.0, 'FeOt_Liq': 9.0,
    'MgO_Liq': 9.0, 'CaO_Liq': 10.0, 'Na2O_Liq': 2.0, 'K2O_Liq': 0.4,
    'P2O5_Liq': 0.2, 'Fe3Fet_Liq': 0.2, 'H2O_Liq': 2.0,
}


def _worker_hangs(q, index, **kwargs):
    """Reports 'Start' then hangs, like a calculation stuck in the engine.
    The sleep is finite so a worker that escapes cleanup can't live forever."""
    for i in list(index)[:3]:
        q.put([i, -1, "Start"])
    time.sleep(120)


def _fail_the_parent_soon(mode, main_thread_id):
    # Wait until workers exist and the parent is inside its wait loop.
    deadline = time.time() + 60
    while time.time() < deadline:
        if multiprocessing.active_children():
            break
        time.sleep(0.05)
    time.sleep(3)
    if mode == "interrupt":
        _thread.interrupt_main()  # what stopping a notebook cell does
    else:
        ctypes.pythonapi.PyThreadState_SetAsyncExc(
            ctypes.c_ulong(main_thread_id), ctypes.py_object(RuntimeError)
        )


def main(mode):
    ptt_path = importlib.import_module("petthermotools.Path")
    ptt_path.path_multi = _worker_hangs

    threading.Thread(
        target=_fail_the_parent_soon,
        args=(mode, threading.main_thread().ident),
        daemon=True,
    ).start()

    caught = None
    try:
        ptt_path.multi_path(
            Model="MELTSv1.0.2",
            comp=COMP,
            T_C=1200.0,
            P_bar=np.array([1000.0, 2000.0, 3000.0, 4000.0, 5000.0]),
            cores=2,
            # Large, so multi_path's own timeout can't be what stops the workers.
            timeout=600,
            multi_processing=True,
            Print_suppress=True,
        )
    except BaseException as exc:  # a notebook catches this and carries on
        caught = type(exc).__name__

    time.sleep(0.5)  # let any cleanup the parent did take effect
    alive = [p.pid for p in multiprocessing.active_children() if p.is_alive()]
    print(f"CAUGHT={caught}")
    print("ALIVE=" + ",".join(str(pid) for pid in alive), flush=True)
    os._exit(0)


if __name__ == "__main__":
    main(sys.argv[1])
