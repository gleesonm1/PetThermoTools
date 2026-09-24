"""
Case registry shared by the MELTS golden-output tests, the golden capture
script and the benchmark harness. Not a test module (leading underscore).

Each case is (public function name, kwargs). The cases are small (a few
seconds each) and cover the calculation families the cleanup plan touches:
isobaric / fractional / polybaric / decompression paths, a batch run (which
exercises the multiprocessing supervisor), and single-point equilibrate and
liquidus.

Running a case calls the real alphaMELTS engine, so it needs `meltsdynamic`.
The engine writes hundreds of KB of progress text to stderr (its C library
writes directly to file descriptor 2) and MELTS drops *_tbl.txt / .inp files in
the working directory, so callers should use `quiet_fds()` and run from a
scratch dir.
"""
import contextlib
import os
import sys
import tempfile

import numpy as np
import pandas as pd

COMP = {
    'SiO2_Liq': 52.72, 'TiO2_Liq': 2.08, 'Al2O3_Liq': 13.26, 'FeOt_Liq': 9.21,
    'MgO_Liq': 9.32, 'CaO_Liq': 10.26, 'Na2O_Liq': 2.24, 'K2O_Liq': 0.43,
    'P2O5_Liq': 0.22, 'Fe3Fet_Liq': 0.15, 'H2O_Liq': 0.5,
}

_PATH = {'Model': "MELTSv1.0.2", 'bulk': COMP, 'H2O_init': 0.5, 'Fe3Fet_init': 0.15, 'timeout': 120}

CASES = {
    "isobaric_v102_1kbar": ("isobaric_crystallisation", {
        **_PATH, 'T_start_C': 1300, 'T_end_C': 1100, 'dt_C': 10, 'P_bar': 1000, 'find_liquidus': True}),
    "isobaric_fractional_v102_2kbar": ("isobaric_crystallisation", {
        **_PATH, 'T_start_C': 1300, 'T_end_C': 1100, 'dt_C': 10, 'P_bar': 2000, 'find_liquidus': True,
        'Frac_solid': True}),
    "isobaric_pmelts_10kbar": ("isobaric_crystallisation", {
        **_PATH, 'Model': "pMELTS", 'T_start_C': 1450, 'T_end_C': 1250, 'dt_C': 10, 'P_bar': 10000,
        'find_liquidus': True}),
    "isobaric_v120_wet_buffered": ("isobaric_crystallisation", {
        **_PATH, 'Model': "MELTSv1.2.0", 'H2O_init': 3.0, 'T_start_C': 1250, 'T_end_C': 1050, 'dt_C': 10,
        'P_bar': 3000, 'find_liquidus': True, 'fO2_buffer': "FMQ", 'fO2_offset': 0.0}),
    "polybaric_v102": ("polybaric_crystallisation_path", {
        **_PATH, 'T_start_C': 1300, 'T_end_C': 1100, 'dt_C': 10, 'P_start_bar': 5000, 'P_end_bar': 1000,
        'dp_bar': 200, 'find_liquidus': True}),
    "isochoric_v102": ("isochoric_crystallisation", {
        **_PATH, 'T_start_C': 1300, 'T_end_C': 1100, 'dt_C': 10, 'P_bar': 1000, 'find_liquidus': True}),
    "isothermal_decompression_v102": ("isothermal_decompression", {
        **_PATH, 'T_C': 1250, 'P_start_bar': 5000, 'P_end_bar': 500, 'dp_bar': 500}),
    "isentropic_decompression_v102": ("isentropic_decompression", {
        **_PATH, 'T_C': 1300, 'P_start_bar': 5000, 'P_end_bar': 500, 'dp_bar': 500}),
    # Batch: several runs across worker processes (exercises path_multi and the
    # supervisor). `cores` is left to the package default, so this also shows
    # the results do not depend on the core count.
    "batch_isobaric_v102_6P": ("isobaric_crystallisation", {
        **_PATH, 'T_start_C': 1300, 'T_end_C': 1100, 'dt_C': 10,
        'P_bar': np.array([500., 1000., 2000., 3000., 4000., 5000.]), 'find_liquidus': True}),
    "equilibrate_v102": ("equilibrate_multi", {
        'Model': "MELTSv1.0.2", 'bulk': COMP, 'T_C': 1150.0, 'P_bar': 1000.0,
        'H2O_init': 0.5, 'Fe3Fet_init': 0.15, 'timeout': 60}),
    "findliq_v102": ("findLiq_multi", {
        'Model': "MELTSv1.0.2", 'bulk': COMP, 'T_initial_C': 1300.0, 'P_bar': 1000.0,
        'H2O_Liq': 0.5, 'Fe3Fet_Liq': 0.15}),
}


@contextlib.contextmanager
def quiet_fds(keep_stderr=True):
    """Silence the MELTS engine's progress text (hundreds of KB per run).

    The C library writes it to file descriptor 2 (stderr) -- measured: ~420 KB
    on stderr and nothing on stdout for one small run -- and worker processes
    inherit the descriptors, so redirecting sys.stdout/sys.stderr alone would
    not silence it. stdout goes to /dev/null.

    keep_stderr=True (default; for correctness runs): stderr is captured to a
    temp file so that, if the block raises, the lines that are NOT engine
    progress (e.g. a worker's traceback) are shown instead of being lost.

    keep_stderr=False (use when TIMING): stderr also goes to /dev/null.
    Capturing to a file makes every worker write its ~420 KB of small
    unbuffered writes to the SAME file, which serialises them: measured on
    Windows it made an 8-run/8-worker batch take 9.4s instead of 4.7s."""
    sys.stdout.flush()
    sys.stderr.flush()
    saved = {fd: os.dup(fd) for fd in (1, 2)}
    devnull = os.open(os.devnull, os.O_WRONLY)
    log = tempfile.TemporaryFile("w+b") if keep_stderr else None

    def restore():
        sys.stdout.flush()
        sys.stderr.flush()
        for fd, dup in saved.items():
            os.dup2(dup, fd)

    try:
        os.dup2(devnull, 1)
        os.dup2(log.fileno() if log is not None else devnull, 2)
        try:
            yield
        except BaseException:
            restore()
            if log is not None:
                log.seek(0)
                keep = [ln for ln in log.read().decode("utf-8", "replace").splitlines()
                        if ln.strip() and not ln.startswith(("...", "Processed line"))]
                if keep:
                    print("--- captured stderr (engine progress text removed) ---", file=sys.stderr)
                    print("\n".join(keep[-30:]), file=sys.stderr)
            raise
    finally:
        restore()
        for dup in saved.values():
            os.close(dup)
        os.close(devnull)
        if log is not None:
            log.close()


def run_case(name):
    """Run one case against the real engine and return the package's result."""
    import petthermotools as ptt

    fn_name, kwargs = CASES[name]
    return getattr(ptt, fn_name)(**kwargs)


def flatten_frames(result, prefix=""):
    """Reduce a result (DataFrame, Series, or arbitrarily nested dict of them)
    to {"a/b": DataFrame}. Non-tabular entries (e.g. the echoed 'Input' dict
    of batch runs) are ignored."""
    frames = {}
    if isinstance(result, pd.DataFrame):
        frames[prefix or "result"] = result
    elif isinstance(result, pd.Series):
        frames[prefix or "result"] = result.to_frame()
    elif isinstance(result, dict):
        for key, value in result.items():
            frames.update(flatten_frames(value, f"{prefix}/{key}" if prefix else str(key)))
    return frames
