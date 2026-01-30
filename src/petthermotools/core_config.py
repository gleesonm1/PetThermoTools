# core_config.py
import platform
import subprocess
import psutil
import logging


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
    logger.warning("Greater than 32 performance cores identified. This may indicate you're using a shared environment and we've accidentally Defaulting to 8.")
    MAX_WORKERS = 8
else:
    MAX_WORKERS = max(1, _p_cores)