"""
Shared helpers for tests that run multi_path scenarios in a separate
interpreter. Not a test module (leading underscore keeps pytest from
collecting it); importable from sibling test files because pytest puts this
directory on sys.path.
"""


import os
import subprocess
import sys
import time


class ScriptRun:
    """Outcome of run_script_capped."""

    def __init__(self):
        self.started = False    # the script printed its start marker
        self.finished = False   # it exited within cap_s of starting
        self.returncode = None
        self.out = ""
        self.err = ""


def run_script_capped(script_args, workdir, cap_s, start_marker="STARTED=1", start_timeout_s=180):
    """Run `python <script_args>` in its own interpreter and bound how long it
    may run. The clock starts when the script prints `start_marker` (after its
    imports), so a slow import on a busy machine cannot look like a hang. If
    it has not exited cap_s after that, the whole process tree is killed and
    `finished` stays False. Output goes to files (not pipes: surviving
    children would keep a pipe open) and is read back as UTF-8."""
    out_path, err_path = workdir / "run.out", workdir / "run.err"
    run = ScriptRun()
    with open(out_path, "w") as out, open(err_path, "w") as err:
        proc = subprocess.Popen(
            [sys.executable, "-W", "ignore", *map(str, script_args)],
            stdout=out, stderr=err, env=utf8_env(),
        )
    t_launch = time.time()
    started_at = None
    try:
        while proc.poll() is None:
            now = time.time()
            if started_at is None and start_marker in read_utf8(out_path):
                started_at = now
            if started_at is not None and now - started_at > cap_s:
                break
            if now - t_launch > start_timeout_s:
                break
            time.sleep(0.2)
        run.finished = proc.poll() is not None
        run.returncode = proc.poll()
        run.started = started_at is not None or start_marker in read_utf8(out_path)
    finally:
        kill_process_tree(proc.pid)  # no-op if it already exited
        run.out, run.err = read_utf8(out_path), read_utf8(err_path)
    return run


def utf8_env():
    """Environment for child interpreters whose output the harness reads back.

    Without this, a child's stdout encoding follows the machine (cp1252 on a
    default Windows shell, UTF-8 if PYTHONIOENCODING/PYTHONUTF8 happens to be
    set), so a single non-cp1252 character in a child's output either crashes
    the child or makes the parent's read fail -- depending on the environment
    the tests are run from. Pin both ends to UTF-8.
    """
    return {**os.environ, "PYTHONIOENCODING": "utf-8"}


def read_utf8(path):
    """Read a file a child wrote (see utf8_env). Never raises on bad bytes."""
    return path.read_text(encoding="utf-8", errors="replace")


def kill_process_tree(pid):
    """Kill `pid` and every descendant. Safe to call on a PID that is
    already gone."""
    import psutil  # a declared install_requires dependency of petthermotools

    try:
        root = psutil.Process(pid)
        root.suspend()  # stop it respawning workers while we enumerate them
        victims = root.children(recursive=True) + [root]
    except psutil.NoSuchProcess:
        return
    for proc in victims:
        try:
            proc.kill()
        except psutil.NoSuchProcess:
            pass
    psutil.wait_procs(victims, timeout=10)
