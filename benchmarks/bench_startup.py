"""
Worker start-up cost, with and without contention -- no MELTS calculation
involved, so it takes about a minute and can be run on any machine.

    python benchmarks/bench_startup.py
    python benchmarks/bench_startup.py --concurrency 1 2 4 8 16 --repeats 5

Why: every worker the package spawns (macOS/Windows use the "spawn" start
method) is a fresh interpreter that re-imports petthermotools, which pulls in
torch (via ngibbs), scipy, matplotlib, pandas and IPython/ipywidgets even
though a MELTS worker needs none of them. In the MELTS benchmark, wall time
for small batches is dominated by this, and it gets WORSE with more workers
because they all import at once (competing for disk, CPU and, on Windows,
antivirus scanning of the same DLLs).

For each concurrency level N this launches N fresh interpreters simultaneously
and times how long until all N have finished importing:
    package : `import petthermotools`                (what a worker pays today)
    slim    : numpy + pandas (+ meltsdynamic if importable)   (roughly what a
              MELTS worker actually needs)
The gap between the two is the saving available from a lighter worker
(e.g. a lazy package __init__ or a minimal worker module); the growth of each
row with N is the contention cost.

Results are printed and written to benchmarks/results/<host>_<utc>_startup.json.
"""
import argparse
import json
import os
import platform
import socket
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

_CODE = {
    "package": "import petthermotools",
    "slim": "import numpy, pandas\ntry:\n    import meltsdynamic\nexcept ImportError:\n    pass",
}


def launch_and_wait(code, n):
    env = {**os.environ, "PYTHONIOENCODING": "utf-8"}
    env["PYTHONPATH"] = os.pathsep.join([str(REPO / "src"), env.get("PYTHONPATH", "")]).rstrip(os.pathsep)
    t0 = time.perf_counter()
    procs = [
        subprocess.Popen([sys.executable, "-W", "ignore", "-c", code], env=env,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for _ in range(n)
    ]
    codes = [p.wait() for p in procs]
    elapsed = time.perf_counter() - t0
    if any(codes):
        raise RuntimeError(f"a child exited with failure codes {codes}: {code!r}")
    return elapsed


def main(args):
    env = {
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "host": socket.gethostname(), "platform": platform.platform(),
        "cpu_logical": os.cpu_count(), "python": platform.python_version(),
    }
    print("environment:", json.dumps(env))

    for code in _CODE.values():  # warm the OS file cache so runs compare fairly
        launch_and_wait(code, 1)

    rows = []
    for target, code in _CODE.items():
        for n in args.concurrency:
            times = [launch_and_wait(code, n) for _ in range(args.repeats)]
            rows.append({"target": target, "concurrent": n, "median_s": round(statistics.median(times), 3),
                         "min_s": round(min(times), 3), "all_s": [round(t, 3) for t in times]})
    base = {r["target"]: r["median_s"] for r in rows if r["concurrent"] == 1}
    print(f"\n{'target':<8} {'N':>3} {'median':>8} {'min':>8}   x slower than N=1")
    for r in rows:
        ratio = r["median_s"] / base[r["target"]] if r["target"] in base else float("nan")
        print(f"{r['target']:<8} {r['concurrent']:>3} {r['median_s']:>7.2f}s {r['min_s']:>7.2f}s   x{ratio:.2f}")

    out_dir = REPO / "benchmarks" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{env['host']}_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}_startup.json"
    path.write_text(json.dumps({"environment": env, "rows": rows}, indent=1))
    print(f"\nwrote {path.relative_to(REPO)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--concurrency", type=int, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--repeats", type=int, default=3)
    main(parser.parse_args())
