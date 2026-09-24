"""
Wall-clock and memory benchmark for MELTS crystallisation runs through
petthermotools.multi_path -- the baseline the speed work (cleanup plan, Part
(ii)) is measured against.

    python benchmarks/bench_melts.py                          # defaults
    python benchmarks/bench_melts.py --cores 1 2 4 8 --runs 1 8 24 --repeats 3
    python benchmarks/bench_melts.py --label after-worker-pool

Each scenario is N isobaric-crystallisation runs (1300 -> 1100 C in 10 C steps
at N pressures between 0.5 and 5 kbar; N=1 is a single run at 1 kbar) on
`cores` worker processes. Recorded per scenario: wall time (perf_counter around
multi_path), peak number of processes in the tree, and peak summed resident
memory of the tree (summed per-process RSS, so shared pages are counted more
than once: use it to compare before/after, not as an absolute figure).

One untimed warm-up call runs first (cold start is visibly slower). Results are
printed and written to benchmarks/results/<host>_<utc>_<label>.json together
with the environment (platform, python, pandas, numpy, git commit), so runs on
different machines/versions can be compared later. For meaningful numbers close
other heavy applications; the timings include process start-up on purpose,
because that is a cost users pay.

Runs from a scratch directory (MELTS drops files in the CWD) with the engine's
output (hundreds of KB per run, on stderr) sent to /dev/null. Not to a shared
capture file: with many workers writing to one file the writes serialise and
the benchmark then measures its own capture (seen on Windows: 8 workers took
9.4s instead of 4.7s).
"""
import argparse
import atexit
import json
import os
import platform
import shutil
import socket
import statistics
import subprocess
import sys
import tempfile
import threading
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "Unit_tests" / "integration" / "alphamelts"))
warnings.simplefilter("ignore")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import psutil  # noqa: E402

import _cases  # noqa: E402


class _TreeSampler(threading.Thread):
    """Samples this process and all descendants: peak count and peak summed RSS."""

    def __init__(self, interval=0.1):
        super().__init__(daemon=True)
        self.interval = interval
        self.stop_event = threading.Event()
        self.peak_procs = 0
        self.peak_rss_mb = 0.0

    def run(self):
        me = psutil.Process()
        while not self.stop_event.is_set():
            try:
                procs = [me] + me.children(recursive=True)
                rss = 0
                for p in procs:
                    try:
                        rss += p.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                self.peak_procs = max(self.peak_procs, len(procs))
                self.peak_rss_mb = max(self.peak_rss_mb, rss / 1e6)
            except psutil.NoSuchProcess:
                pass
            self.stop_event.wait(self.interval)


def run_scenario(n_runs, cores):
    import petthermotools as ptt

    kwargs = {
        'Model': "MELTSv1.0.2", 'bulk': _cases.COMP, 'H2O_init': 0.5, 'Fe3Fet_init': 0.15,
        'T_start_C': 1300, 'T_end_C': 1100, 'dt_C': 10, 'find_liquidus': True,
        'P_bar': np.linspace(500.0, 5000.0, n_runs) if n_runs > 1 else 1000.0,
        'cores': cores, 'timeout': 300, 'Print_suppress': True,
    }
    sampler = _TreeSampler()
    sampler.start()
    t0 = time.perf_counter()
    # stderr to /dev/null, not to a shared capture file: see quiet_fds().
    with _cases.quiet_fds(keep_stderr=False):
        result = ptt.multi_path(**kwargs)
    wall = time.perf_counter() - t0
    sampler.stop_event.set()
    sampler.join(timeout=2)
    completed = len([k for k in result if str(k).startswith("Run ")]) if n_runs > 1 else int(bool(result))
    return {"wall_s": wall, "peak_procs": sampler.peak_procs,
            "peak_rss_mb": sampler.peak_rss_mb, "completed_runs": completed}


def _git_commit():
    try:
        out = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() or None
    except Exception:  # noqa: BLE001
        return None


def main(args):
    import petthermotools as ptt
    from petthermotools import core_config

    workdir = tempfile.mkdtemp(prefix="ptt_bench_")
    os.chdir(workdir)  # MELTS writes files into the CWD
    # Leave the scratch dir first (Windows cannot delete a directory that is the CWD).
    atexit.register(lambda: (os.chdir(REPO), shutil.rmtree(workdir, ignore_errors=True)))

    env = {
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "host": socket.gethostname(), "platform": platform.platform(),
        "processor": platform.processor(), "cpu_logical": os.cpu_count(),
        "cpu_physical": psutil.cpu_count(logical=False),
        "python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__,
        "petthermotools": getattr(ptt, "__version__", None), "git_commit": _git_commit(),
        "package_default_workers": getattr(core_config, "MAX_WORKERS", None),
    }
    print("environment:", json.dumps(env))

    print("warm-up (untimed) ...", flush=True)
    run_scenario(1, 1)

    scenarios = []
    for n_runs in args.runs:
        for cores in ([1] if n_runs == 1 else args.cores):
            walls, peaks = [], []
            for _ in range(args.repeats):
                r = run_scenario(n_runs, cores)
                if r["completed_runs"] != n_runs:
                    print(f"WARNING: {n_runs} runs on {cores} cores completed only {r['completed_runs']}")
                walls.append(r["wall_s"])
                peaks.append(r)
            row = {
                "runs": n_runs, "cores": cores, "repeats": args.repeats,
                "wall_s_all": [round(w, 3) for w in walls],
                "wall_s_median": round(statistics.median(walls), 3),
                "wall_s_min": round(min(walls), 3),
                "s_per_run_median": round(statistics.median(walls) / n_runs, 3),
                "peak_procs": max(p["peak_procs"] for p in peaks),
                "peak_rss_mb": round(max(p["peak_rss_mb"] for p in peaks)),
            }
            scenarios.append(row)
            print(f"runs={n_runs:<3} cores={cores:<2} median={row['wall_s_median']:>7.2f}s "
                  f"min={row['wall_s_min']:>7.2f}s  {row['s_per_run_median']:>6.2f} s/run  "
                  f"peak procs={row['peak_procs']:<3} peak RSS={row['peak_rss_mb']} MB", flush=True)

    base = {s["runs"]: s["wall_s_median"] for s in scenarios if s["cores"] == 1}
    print("\nspeed-up vs 1 core (median wall time):")
    for s in scenarios:
        if s["runs"] in base and s["runs"] > 1:
            print(f"  runs={s['runs']:<3} cores={s['cores']:<2} x{base[s['runs']] / s['wall_s_median']:.2f}")

    out_dir = REPO / "benchmarks" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    label = f"_{args.label}" if args.label else ""
    out_path = out_dir / f"{env['host']}_{stamp}{label}.json"
    out_path.write_text(json.dumps({"environment": env, "scenarios": scenarios}, indent=1))
    print(f"\nwrote {out_path.relative_to(REPO)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cores", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument("--runs", type=int, nargs="+", default=[1, 8])
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--label", default="", help="tag for the results file, e.g. 'baseline'")
    main(parser.parse_args())
