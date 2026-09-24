# Benchmarks

Timing of the multiprocessing path, used to check that speed changes actually help. Run from the repo root; both scripts
work from a scratch directory and clean it up.

| Script | Needs | What it measures |
|---|---|---|
| `bench_melts.py` | alphaMELTS | wall time, peak process count and summed memory for 1 run and for batches of N runs on C worker processes |
| `bench_startup.py` | nothing (no MELTS) | time for N fresh interpreters started at once to import `petthermotools` (what every spawned worker pays) versus a slim import (numpy + pandas [+ meltsdynamic]) |

```
python benchmarks/bench_melts.py --label before-change
python benchmarks/bench_melts.py --cores 1 2 4 8 --runs 1 8 24 --repeats 3
python benchmarks/bench_startup.py --concurrency 1 2 4 8 16 32
```

Raw output goes to `benchmarks/results/` (git-ignored: it contains host names). Curated, anonymised copies of the
reference measurements are in `baselines/`; copy a result there (replace the `host` field) to keep it.

## Baselines (in `baselines/`)

Measured before any change to the package's calculation code (HEAD `f23a508` plus the new tests). Engine stderr was
sent to the null device (see below). Medians of 2 repeats after an untimed warm-up.

| | macOS arm64, 8 cores, py3.11.0, pandas 3.0.2 (ngibbs/torch installed) | Windows 11, Ryzen Threadripper PRO 5955WX 16C/32T, py3.13.9, pandas 3.0.6 (no torch) |
|---|---|---|
| 1 run | 3.3 s | 2.3 s |
| 8 runs, 1 core | 10.8 s | 9.0 s |
| 8 runs, 2 cores | 8.4 s | 6.3 s |
| 8 runs, 4 cores | 5.5 s | 4.3 s |
| 8 runs, 8 cores | 9.4 s | 4.9 s |
| import `petthermotools`, 1 at a time / 8 at once | 2.9 s / 6.5 s | 1.3 s / 1.7 s (3.9 s at 32) |
| slim import, 1 at a time / 8 at once | 0.37 s / 0.73 s | 0.47 s / 0.65 s (1.5 s at 32) |

Reading: for these small runs almost all the wall time is starting worker processes and importing the package, not
MELTS. More workers make it worse once the start-ups compete for CPU (and, on the Mac, memory: 3.1 GB at 8 workers).
The size of the start-up cost depends on what is installed (torch adds about 1.5 s per worker).

## Pitfalls

- **Silence the engine's stderr with the null device when timing.** The MELTS C library writes about 420 KB per run to
  file descriptor 2. Capturing it to one shared file serialises the workers' writes (on Windows it made 8 workers take
  9.4 s instead of 4.7 s). In a Jupyter kernel the same output goes to the kernel's own stderr pipe and cost about 40 %
  of an 8-run call (22 s -> 13 s on Windows). `bench_melts.py` uses `quiet_fds(keep_stderr=False)`.
- Close other heavy applications; Windows Defender scanning and a busy machine both change the numbers. The first import
  after an install or reboot can be several times slower than a warm one.
- Compare like with like: same machine, same installed packages (torch in particular), same `--cores`.
