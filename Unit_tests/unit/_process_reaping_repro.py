"""
Helper script for test_process_reaping.py -- NOT a test module (the leading
underscore keeps pytest from collecting it).

Runs multi_path with a fast fake worker in its own interpreter (so a hang is
bounded by killing the process tree from outside), then checks that every
worker process it started has been reaped.

How "reaped" is detected, without replacing the Process class:
multiprocessing keeps every started child in `multiprocessing.process._children`
until the parent collects its exit status -- `.join()`, or `is_alive()`/a later
`start()` finding it has exited (each calls waitpid on POSIX). Any worker still
listed after multi_path returns is unreaped (a zombie on POSIX). This checks
the outcome, not which call did the reaping: run_tasks calls is_alive() on
every worker it retires, which already reaps an exited one, so its extra
join() is not what prevents zombies. An earlier version swapped a Process subclass into
mp_utils._CTX.Process and imported petthermotools at module level (so every
spawned child re-imported it during bootstrap); that combination hung on
ubuntu/Python 3.10 CI only, while the identical multi_path path without it
(_multi_path_repro.py "completes") passes there.

faulthandler: this module's top level also runs in every spawned child, so a
process still alive after _DUMP_AFTER_S prints every thread's stack to stderr
-- if this ever hangs again, the test failure shows where.

Prints:
    STARTED=1                     just before multi_path is called
    RESULT=OK                     every worker exited and was reaped
    RESULT=NO_WORKERS             the fake worker never ran
    RESULT=STILL_ALIVE:<n>        n worker(s) still alive after multi_path returned
    RESULT=UNJOINED:<n>           n exited worker(s) were never reaped

Usage:  python _process_reaping_repro.py <pid-log-file>
"""
import faulthandler
import importlib
import multiprocessing
import os
import sys
from pathlib import Path

import numpy as np

_DUMP_AFTER_S = 40   # < test_process_reaping._CAP_S, so the dump lands before the kill
faulthandler.dump_traceback_later(_DUMP_AFTER_S, exit=False)

# Test the working tree, not a stale installed copy (mirrors conftest.py).
_SRC = Path(__file__).resolve().parents[2] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

COMP = {
    'SiO2_Liq': 52.0, 'TiO2_Liq': 2.0, 'Al2O3_Liq': 13.0, 'FeOt_Liq': 9.0,
    'MgO_Liq': 9.0, 'CaO_Liq': 10.0, 'Na2O_Liq': 2.0, 'K2O_Liq': 0.4,
    'P2O5_Liq': 0.2, 'Fe3Fet_Liq': 0.2, 'H2O_Liq': 2.0,
}


def _fake_worker_completes_instantly(q, task, **kwargs):
    """Fake _path_task_v2: records that it ran, reports one step, exits normally
    -- always the natural-exit branch, never the timeout branch (that is
    covered by test_timeout_handling.py)."""
    with open(os.environ["PTT_REAPING_PID_LOG"], "a") as f:
        f.write(f"{os.getpid()}\n")
    q.put([task.row_id, 0, {"Conditions": {}, "_dummy_keys": []}])


def main(pid_log):
    os.environ["PTT_REAPING_PID_LOG"] = pid_log   # inherited by spawned workers
    ptt_path = importlib.import_module("petthermotools.Path")
    ptt_path._path_task_v2 = _fake_worker_completes_instantly

    print("STARTED=1", flush=True)
    ptt_path.multi_path(
        Model="MELTSv1.0.2",
        comp=COMP,
        T_C=1200.0,
        P_bar=np.array([1000.0, 2000.0]),  # one row per core, see _multi_path_repro.py
        cores=2,
        timeout=30,
        multi_processing=True,
        Print_suppress=True,
    )

    # Snapshot before anything that polls children (is_alive() reaps a zombie).
    leftover = list(multiprocessing.process._children)
    ran = Path(pid_log).read_text().split() if Path(pid_log).exists() else []

    if not ran:
        result = "NO_WORKERS"
    elif leftover:
        alive = sum(1 for p in leftover if p.is_alive())
        result = f"STILL_ALIVE:{alive}" if alive else f"UNJOINED:{len(leftover)}"
    else:
        result = "OK"

    print(f"RESULT={result}", flush=True)
    os._exit(0)


if __name__ == "__main__":
    main(sys.argv[1])
