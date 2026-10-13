# core_config.py
import platform
import subprocess
import psutil
import logging
import multiprocessing
import warnings
import numbers


def get_performance_core_count():
    system = platform.system()
    arch = platform.machine()
    
    if system == "Darwin":
        if arch == "arm64":
            try:
                cmd = ["sysctl", "-n", "hw.perflevel0.physicalcpu"]
                return int(subprocess.check_output(cmd).decode().strip())
            except:
                return psutil.cpu_count(logical=False) or 1
        return psutil.cpu_count(logical=False) or 1

    lt = psutil.cpu_count(logical=True) or 1
    lf = psutil.cpu_count(logical=False) or 1
    if lt > lf:
        p_cores = lt - lf
        return p_cores if p_cores > 0 else lf
    return lf

# Set the constant here
_p_cores = max(1, get_performance_core_count())

logger = logging.getLogger(__name__)

if _p_cores > 32:
    logger.warning("Greater than 32 performance cores identified. This may indicate you're using a shared environment and we've accidentally identified the number of cores on the entire cluster! Defaulting to 8.")
    MAX_WORKERS = 8
else:
    MAX_WORKERS = max(1, _p_cores)

def resolve_cores(cores):
    '''
    Resolves the `cores` argument accepted by every parallel calculation
    function in PetThermoTools (multi_path, equilibrate_multi, findLiq_multi,
    findCO2_multi, AdiabaticDecompressionMelting, findSatPressure/_multi,
    mineral_cosaturation/find_mineral_cosaturation, phaseDiagram_calc,
    saturation_pressure). Centralised here so all of them accept the same
    values and give the same warning, rather than each re-implementing it.

    Parameters
    ----------
    cores : None, int, or {"performance", "all"}
        None
            The package default -- MAX_WORKERS (performance cores only;
            see get_performance_core_count's own docstring). Leaves
            efficiency cores free for other work on the machine.
        int
            Used exactly as given (unchanged, existing behaviour) -- must
            be >= 1.
        "performance"
            Explicit alias for the default (MAX_WORKERS). Same value,
            spelled out for clarity in a script, no warning.
        "all"
            Every logical core (performance + efficiency), via
            multiprocessing.cpu_count(). Can speed up a batch, but leaves
            no spare cores for anything else running on the machine at the
            same time, and may reduce battery life / raise fan noise on a
            laptop. For calculations that still wait for the slowest
            worker in a group before starting the next one (most of these
            functions, pre-dating the persistent-worker-pool rewrite of
            findLiq_multi -- see that function's own docstring and
            MAINTENANCE_NOTES.md), contention from something else using an
            efficiency core at the same time can occasionally make the
            batch SLOWER than the default, not faster. Prints a warning
            every time it's used, not just the first.

    Returns
    -------
    int
        The resolved worker/thread count to actually use.
    '''
    if cores is None:
        return MAX_WORKERS

    if isinstance(cores, str):
        key = cores.strip().lower()
        if key == "performance":
            return MAX_WORKERS
        if key == "all":
            all_cores = multiprocessing.cpu_count()
            warnings.warn(
                f"cores='all' requested -- using all {all_cores} logical cores "
                f"(performance + efficiency) instead of the default, "
                f"performance-core-only {MAX_WORKERS}. This can speed up this "
                f"calculation, but leaves no spare cores for anything else "
                f"running on this machine at the same time, and may reduce "
                f"battery life or increase heat/fan noise on a laptop. If "
                f"something else ends up competing for a core while this "
                f"batch is running, it can even come out SLOWER than the "
                f"default, not faster, for calculations that wait for the "
                f"slowest worker each round.",
                stacklevel=3,
            )
            return all_cores
        raise ValueError(
            f"cores={cores!r} is not a recognised value. Use an integer, "
            f"'performance' ({MAX_WORKERS} on this machine), or 'all' "
            f"({multiprocessing.cpu_count()} on this machine)."
        )

    # numbers.Integral (not plain `int`) so numpy integer scalars pass too --
    # e.g. `cores = some_array[i]` is np.int64, not a Python int, but is a
    # perfectly good core count. bool is technically an int/Integral
    # subclass in Python, so it needs its own explicit exclusion either way.
    if isinstance(cores, bool) or not isinstance(cores, numbers.Integral):
        raise ValueError(
            f"cores must be an int, 'performance', or 'all'; got {type(cores).__name__}."
        )
    if cores < 1:
        raise ValueError(f"cores must be >= 1, got {cores}.")
    return int(cores)  # normalise numpy/other Integral types to plain int