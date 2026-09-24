"""
Tier A failure-injection tests: what multi_path does when a worker hangs past
the timeout.

Each scenario (see _timeout_repro.py) runs in its own interpreter with a fake
worker that hangs at a different point, and the scenarios run concurrently
from one module-scoped fixture to keep the cost down. Each scenario's clock
starts only once its interpreter has finished importing, so slow imports on a
busy CI machine cannot masquerade as a hang.

What happens today (macOS; each scenario is one test below):
  - workers that hang AFTER reporting are terminated and reaped (good), but the
    call then returns an empty dict (nothing reported a result) or a result
    that looks complete (some steps arrived) with no error, no warning and no
    status: the only signal is a print();
  - a worker that hangs BEFORE reporting anything is never noticed: the loop
    only advances on queue messages, so it respawns forever;
  - on macOS/Linux a worker that ignores SIGTERM blocks the parent forever:
    Process.terminate() is SIGTERM and is followed by an unbounded join().

Harness problems use pytest.fail() (not AssertionError) so the
xfail(raises=AssertionError) marks can only ever be satisfied by the real bug.
"""
import subprocess
import sys
import time
from pathlib import Path

import pytest
from _proc_helpers import kill_process_tree, read_utf8, utf8_env

_CHILD = Path(__file__).resolve().parent / "_timeout_repro.py"

# Seconds a scenario may run inside multi_path before it counts as "did not
# return". Several multiples of a full round (2 x multi_path's timeout of 4s).
_RETURN_CAP_S = 25
# Whole-fixture safety net (covers interpreter start-up under heavy load).
_OVERALL_CAP_S = 240

_ALL = [
    "batch_report_then_hang",
    "batch_step_then_hang",
    "batch_hang_before_report",
    "batch_ignores_sigterm",
    "single_hang_before_report",
    "single_step_then_hang",
]
_SCENARIOS = [m for m in _ALL if not (m == "batch_ignores_sigterm" and sys.platform == "win32")]


class _Run:
    def __init__(self, name, workdir):
        self.name = name
        self.path = workdir / f"{name}.out"
        with open(self.path, "w") as out, open(workdir / f"{name}.err", "w") as err:
            self.proc = subprocess.Popen(
                [sys.executable, str(_CHILD), name], stdout=out, stderr=err, env=utf8_env()
            )
        self.started_at = None
        self.returned = False   # the child finished and printed its report
        self.hung = False       # still inside multi_path after _RETURN_CAP_S
        self.report = {}

    def text(self):
        try:
            return read_utf8(self.path)
        except OSError:
            return ""

    def parse(self):
        self.report = dict(
            line.split("=", 1) for line in self.text().splitlines() if "=" in line
        )

    # convenience accessors (valid once .returned)
    @property
    def raised(self):
        return self.report.get("OUTCOME") == "RAISED"

    @property
    def exc(self):
        return self.report.get("EXC", "")

    @property
    def empty(self):
        return self.report.get("EMPTY") == "1"

    @property
    def alive(self):
        return [int(p) for p in self.report.get("ALIVE", "").split(",") if p]

    @property
    def signals(self):
        return [s for s in self.report.get("SIGNALS", "").split(" || ") if s]

    def describe(self):
        return (
            f"{self.report.get('OUTCOME')} {self.exc} empty={self.empty} "
            f"signals={self.signals}"
        )


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    workdir = tmp_path_factory.mktemp("timeout_scenarios")
    all_runs = {name: _Run(name, workdir) for name in _SCENARIOS}
    killed = []
    try:
        t_launch = time.time()
        while True:
            pending = [r for r in all_runs.values() if not (r.returned or r.hung)]
            if not pending:
                break
            now = time.time()
            for r in pending:
                if r.proc.poll() is not None:
                    r.parse()
                    if "OUTCOME" not in r.report:
                        pytest.fail(
                            f"harness problem: {r.name} exited without a report.\n"
                            f"stdout:\n{r.text()}\nstderr:\n{read_utf8(workdir / (r.name + '.err'))}"
                        )
                    r.returned = True
                    continue
                if r.started_at is None and "STARTED=1" in r.text():
                    r.started_at = now
                if r.started_at is not None and now - r.started_at > _RETURN_CAP_S:
                    r.hung = True
                    kill_process_tree(r.proc.pid)
                    killed.append(r.name)
            if now - t_launch > _OVERALL_CAP_S:
                pytest.fail(
                    "harness problem: scenarios never got as far as calling "
                    f"multi_path within {_OVERALL_CAP_S}s: "
                    f"{[r.name for r in pending if r.started_at is None]}"
                )
            time.sleep(0.2)
        yield all_runs
    finally:
        # Whatever the outcome, leave nothing running.
        for r in all_runs.values():
            for pid in r.alive:
                kill_process_tree(pid)
            kill_process_tree(r.proc.pid)


def _returned(runs, name):
    """The scenario ended with multi_path returning or raising (not hanging)."""
    r = runs[name]
    if r.hung:
        pytest.fail(f"harness/precondition problem: {name} did not return, so it cannot be inspected")
    return r


def _is_signalled(r):
    """A timeout that kills runs must be detectable through a channel a script
    can see: a deliberate exception or a warning. print() does not count.
    Exactly which of these (or a per-run status in the result) is the right
    design is left open; this only asserts it is not silent."""
    return (r.raised and r.exc in ("TimeoutError", "RuntimeError", "ValueError")) or bool(r.signals)


# --------------------------------------------------------------------------
# Works today: a timeout does kill and reap the workers
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name", [
    "batch_report_then_hang",
    "batch_step_then_hang",
    "single_hang_before_report",
    "single_step_then_hang",
])
def test_timed_out_workers_are_terminated_and_reaped(runs, name):
    r = _returned(runs, name)
    assert not r.alive, f"{len(r.alive)} worker(s) still running after the timeout (pids: {r.alive})"


# --------------------------------------------------------------------------
# Proposed behaviour (xfail today)
# --------------------------------------------------------------------------

_XFAIL = pytest.mark.xfail(
    raises=AssertionError,  # only a genuine mismatch counts as "expected"
    strict=False,
)


@_XFAIL
@pytest.mark.parametrize("name", ["batch_report_then_hang", "single_hang_before_report"])
def test_timeout_that_loses_every_run_is_not_silent(runs, name):
    """Today: returns an empty result; the only trace is a print()."""
    r = _returned(runs, name)
    assert _is_signalled(r), (
        f"every run timed out and nothing was returned, but multi_path gave no "
        f"error or warning ({r.describe()})"
    )


@_XFAIL
@pytest.mark.parametrize("name", ["batch_step_then_hang", "single_step_then_hang"])
def test_timeout_that_truncates_runs_is_signalled(runs, name):
    """Today: returns the partial results as if the runs had completed."""
    r = _returned(runs, name)
    assert _is_signalled(r), (
        f"runs were killed mid-calculation but multi_path returned their partial "
        f"results without any error or warning ({r.describe()})"
    )


@_XFAIL
def test_hang_before_first_report_does_not_loop_forever(runs):
    """Same failure as a worker that crashes before reporting (see
    test_multi_path_supervisor.py): the loop only advances on queue messages,
    so a worker that never sends one is timed out and respawned forever."""
    r = runs["batch_hang_before_report"]
    assert r.returned, (
        f"multi_path was still running {_RETURN_CAP_S}s after it started (a "
        "worker that hangs before reporting is respawned indefinitely)"
    )


@_XFAIL
@pytest.mark.skipif(sys.platform == "win32", reason="Windows terminate() is TerminateProcess, which cannot be ignored")
def test_worker_ignoring_sigterm_does_not_block_forever(runs):
    """Process.terminate() sends SIGTERM and is followed by an unbounded
    join(); a worker stuck where SIGTERM is not acted on blocks the parent
    indefinitely. A fix escalates to kill() after a grace period."""
    r = runs["batch_ignores_sigterm"]
    assert r.returned, (
        f"multi_path was still blocked {_RETURN_CAP_S}s after it started: "
        "terminate() was ignored and join() has no timeout"
    )
