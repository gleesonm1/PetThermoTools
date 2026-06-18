"""
Benchmark: nGibbs ForwardMB vs. PetThermoTools property calculators
(serial calc_phase_props_MELTS and parallel calc_phase_props_MELTS_parallel)
over batch sizes from 128 to 32768.

Grid layout for each batch size N = 2^k:
  n_P = 2^(k//2),  n_T = N // n_P   (balanced P-T grid)

Results are printed as a table and saved to benchmark_results_ngibbs_props.csv
alongside this script.

Adjust NGIBBS_ROOT below if your nGibbs checkout is in a different location.
"""

from pathlib import Path
import sys
import warnings
import multiprocessing

# ---- Path setup ---------------------------------------------------------- #
SCRIPT_DIR = Path(__file__).parent
NGIBBS_ROOT = Path("C:/Git_Repositories/nGibbs")          # adjust if needed

PTT_SRC = Path("C:/Git_Repositories/PetThermoTools/src")

for p in [str(PTT_SRC), str(NGIBBS_ROOT / "src"), str(NGIBBS_ROOT / "src" / "ngibbs")]:
    if p not in sys.path:
        sys.path.insert(0, p)


if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

sys.path.append("C:\\Users\\dashf\\Downloads\\alphaMELTSEnsemble-master\\BJT2025\\PTT_Test\\PetThermoTools-0.2.39\\PetThermoTools-0.2.39\\MELTS\\") # On BJT local machine


# ---- Guard required for ProcessPoolExecutor on Windows ------------------- #
if __name__ == "__main__":
    multiprocessing.freeze_support()

    # ---- Imports --------------------------------------------------------- #
    import numpy as np
    import pandas as pd
    from time import perf_counter
    from tqdm import tqdm

    import torch
    import petthermotools as ptt
    from petthermotools.MELTS import calc_phase_props_MELTS, calc_phase_props_MELTS_parallel

    from ngibbs import MELTS102EmulatorCPU
    from ngibbs.utils.math_utils import grid_sample
    from petthermotools.core_config import MAX_WORKERS

    # ---- Composition ----------------------------------------------------- #
    MORB = {
        "SiO2":  48.68,
        "TiO2":   1.01,
        "Al2O3": 17.64,
        "Cr2O3":  0.03,
        "FeO":    7.59,
        "Fe2O3":  0.89,
        "MgO":    9.10,
        "CaO":   12.45,
        "Na2O":   2.65,
        "K2O":    0.03,
        "P2O5":   0.08,
        "H2O":    0.20,
    }

    comp_array   = np.array([[v for v in MORB.values()]], dtype=np.float32)
    input_headers = ["Pressure(System_main)", "Temperature(System_main)"] + list(MORB.keys())

    # ---- Pre-initialise serial MELTS (exclude library load from timings) - #
    print("Initialising serial MELTS...", end=" ", flush=True)
    from meltsdynamic import MELTSdynamic
    melts_serial = MELTSdynamic(1)
    print("done.\n")

    print(f"Parallel workers: {MAX_WORKERS}\n")

    # ---- Batch sizes ----------------------------------------------------- #
    # 2^7 = 128  …  2^15 = 32 768
    BATCH_SIZES = [2**k for k in range(7, 16)]

    def grid_dims(N: int):
        k = int(round(np.log2(N)))
        n_P = 2 ** (k // 2)
        n_T = N // n_P
        return n_P, n_T

    # ---- Run benchmark --------------------------------------------------- #
    records = []

    for N in tqdm(BATCH_SIZES, desc="batch sizes", unit="batch"):
        n_P, n_T = grid_dims(N)

        nGibbs_input = grid_sample(
            [[2000.0, 10000.0, n_P], [1400.0, 900.0, n_T]],
            comp_array,
        )

        # -- nGibbs timing (ForwardMB only, no table splitting) -------------
        tqdm.write(f"\n  N={N:>6d} ({n_P}P x {n_T}T)")
        tqdm.write(f"    nGibbs...", end=" ")
        t0 = perf_counter()
        nGibbs_output = MELTS102EmulatorCPU.ForwardMB(
            nGibbs_input, headers=input_headers, outputs=["ptt_out"]
        )
        t_ngibbs = perf_counter() - t0
        tqdm.write(f"{t_ngibbs:.2f}s")

        ptt_full = nGibbs_output["ptt_out"]

        # -- Serial: single call on the full undivided table ----------------
        tqdm.write(f"    serial...", end=" ")
        t0 = perf_counter()
        calc_phase_props_MELTS(ptt_full, melts=melts_serial)
        t_serial = perf_counter() - t0
        tqdm.write(f"{t_serial:.2f}s")

        if N == BATCH_SIZES[0]:
            ptt_full['Conditions'].to_csv(SCRIPT_DIR / "serial_N128_Conditions.csv", index=False)
            if 'liquid1_prop' in ptt_full:
                ptt_full['liquid1_prop'].to_csv(SCRIPT_DIR / "serial_N128_liquid1_prop.csv", index=False)

        # Re-run nGibbs to get a fresh table for parallel (serial modified it)
        nGibbs_output2 = MELTS102EmulatorCPU.ForwardMB(
            nGibbs_input, headers=input_headers, outputs=["ptt_out"]
        )
        ptt_full2 = nGibbs_output2["ptt_out"]

        # -- Parallel: single call, one pool, rows distributed across workers
        tqdm.write(f"    parallel...", end=" ")
        t0 = perf_counter()
        calc_phase_props_MELTS_parallel(ptt_full2)
        t_parallel = perf_counter() - t0
        tqdm.write(f"{t_parallel:.2f}s")

        if N == BATCH_SIZES[0]:
            ptt_full2['Conditions'].to_csv(SCRIPT_DIR / "parallel_N128_Conditions.csv", index=False)
            if 'liquid1_prop' in ptt_full2:
                ptt_full2['liquid1_prop'].to_csv(SCRIPT_DIR / "parallel_N128_liquid1_prop.csv", index=False)

        speedup = t_serial / t_parallel if t_parallel > 0 else float("inf")
        ratio_serial = t_serial / t_ngibbs if t_ngibbs > 0 else float("inf")

        tqdm.write(
            f"    speedup {speedup:.1f}x  serial/nGibbs {ratio_serial:.1f}x"
        )

        records.append(dict(
            N=N, n_P=n_P, n_T=n_T,
            t_ngibbs_s=t_ngibbs,
            t_serial_s=t_serial,
            t_parallel_s=t_parallel,
            speedup=speedup,
        ))

    # ---- Summary table --------------------------------------------------- #
    header = (
        f"{'N':>8}  {'n_P':>5}  {'n_T':>6}  "
        f"{'nGibbs (s)':>12}  {'serial (s)':>12}  {'parallel (s)':>14}  "
        f"{'speedup':>9}  {'serial/nGibbs':>15}"
    )
    print("\n" + header)
    print("-" * len(header))
    for r in records:
        ratio = r['t_serial_s'] / r['t_ngibbs_s'] if r['t_ngibbs_s'] > 0 else float('inf')
        print(
            f"{r['N']:>8d}  {r['n_P']:>5d}  {r['n_T']:>6d}  "
            f"{r['t_ngibbs_s']:>12.4f}  {r['t_serial_s']:>12.4f}  {r['t_parallel_s']:>14.4f}  "
            f"{r['speedup']:>8.1f}x  {ratio:>14.1f}x"
        )

    # ---- Save results ---------------------------------------------------- #
    out_path = SCRIPT_DIR / "benchmark_results_ngibbs_props.csv"
    pd.DataFrame(records).to_csv(out_path, index=False)
    print(f"\nResults saved to {out_path}")
