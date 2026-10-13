"""
mp_utils.py -- a shared multiprocessing supervisor for PetThermoTools's
"one task = one freshly-initialised engine process" calculations.

STATUS: standalone sketch, not wired into anything yet. Nothing in the rest
of the package imports this module. It is the concrete form of the Part (ii)
#1 "persistent worker pool" item from the cleanup plan, written so its
structure can be read and exercised in isolation (see the toy example in
``if __name__ == "__main__":`` at the bottom, and the illustrative
``# multi_path, sketched against run_tasks`` section) before any of
multi_path / equilibrate_multi / findCO2_multi / findLiq_multi /
mineral_cosaturation / AdiabaticDecompressionMelting / phaseDiagram_calc are
changed to use it.

Design constraint this preserves: the alphaMELTS C library cannot be
reinitialised in-process after a failure, so a worker that fails or times
out always exits -- it never attempts a second task -- and the supervisor
replaces it with a freshly spawned process. This module only changes how
that spawn/queue/timeout/retry bookkeeping is done, not the one-process-per-
risky-unit-of-work architecture itself.
"""

import time
import multiprocessing as mp
from dataclasses import dataclass
from enum import Enum
from queue import Empty
from typing import Any, Callable, Dict, Iterable, List, Optional


# Explicit "spawn" context everywhere in this module, so macOS / Windows /
# Linux behave identically regardless of the platform default (Linux forks
# by default today; Python 3.14 moves its default to forkserver).
_CTX = mp.get_context("spawn")

# Tags used on the single shared results queue to distinguish supervisor
# bookkeeping messages from a task's own progress reports. task_fn never
# sees or produces these tags itself -- see _ProgressProxy below.
_STARTED = "__mp_utils_started__"
_DONE = "__mp_utils_done__"
_PROGRESS = "__mp_utils_progress__"

# Placed on the task queue to tell an idle worker "no more work, exit
# cleanly" -- identity-checked, so it can never collide with a real index.
_SHUTDOWN = object()


class TaskStatus(str, Enum):
    OK = "ok"
    TIMEOUT = "timeout"    # ran past the absolute per-task `timeout`
    STALLED = "stalled"    # ran past `stall_timeout` with no activity reported
    CRASHED = "crashed"
    GAVE_UP = "gave_up"    # retries exhausted (started life as TIMEOUT/STALLED/CRASHED)


@dataclass
class TaskResult:
    index: Any
    status: TaskStatus
    error: Optional[str] = None  # set for anything other than OK
    attempts: int = 1


@dataclass
class _Tracked:
    """Supervisor-side bookkeeping for the task a given worker is running."""
    index: Any
    process: Any
    started: float
    last_activity: float  # updated on _STARTED and on every _PROGRESS item


class _ProgressProxy:
    """
    What `task_fn` actually gets instead of the raw results_queue: something
    that looks like a queue (``.put(item)``) but transparently tags every
    item with this task's pid/index before it goes on the wire, so the
    supervisor can always attribute a progress message to the right task
    for stall-detection -- regardless of whatever format task_fn itself
    chooses to put() (path_MELTS's `[run_index, step_i, data]`, this
    module's own toy `("progress", i, "halfway")`, or anything else).
    `on_progress` callbacks still see exactly the original, un-wrapped item.
    """
    __slots__ = ("_q", "_pid", "_index")

    def __init__(self, q, pid, index):
        self._q = q
        self._pid = pid
        self._index = index

    def put(self, item):
        self._q.put((_PROGRESS, self._pid, self._index, item))


def _worker_loop(task_fn, task_queue, results_queue, task_kwargs, max_tasks):
    """
    Runs inside a spawned child process for the lifetime of that process.

    Pulls one task index at a time from `task_queue` and calls
    ``task_fn(q, index, **task_kwargs)``, where `q` is a `_ProgressProxy`
    task_fn can call ``q.put(item)`` on exactly like a real
    multiprocessing.Queue -- the same calling convention Path.path() and
    the per-index body of Path.path_multi() already use, so an existing
    single-run worker body is most of the way to being a valid task_fn (it
    would lose its own ``for i in index:`` loop and its 3-index cap, since
    looping across tasks is now this function's job, not the task_fn's).

    Exit conditions, all reported to the supervisor via a `retiring=True`
    flag on the `_DONE` message so it can clean this process up THIS
    iteration rather than waiting for a later liveness-sweep pass to
    notice it has died (see run_tasks's _DONE handling for why that lag
    matters: a worker left in `workers` after it has already exited is a
    real, observed resource leak, not a hypothetical one):
    - cleanly, after pulling the `_SHUTDOWN` sentinel (supervisor is
      winding down and has no more work for this worker) -- no `_DONE` is
      sent for this case, there is no task in flight to report on.
    - immediately, if `task_fn` raises for the task it is currently
      running. This is deliberate: "new process = clean C library" is
      preserved by never attempting a second task in a process that just
      failed one. The supervisor sees the CRASHED report, retires this
      process, and spawns a fresh replacement if there is still work left.
    - after `max_tasks` *successful* tasks (if `max_tasks` is not None).
      This is the important one: `task_fn` reporting success only means it
      didn't raise -- it does NOT mean the engine state behind it is
      trustworthy. path_MELTS, for one real example, wraps its own
      liquidus search in a bare `except:` and returns *normally* on
      failure (MELTS.py: `except: print('Liquidus calculation failed.');
      return Results`), so a task_fn built on it could report OK on a
      process whose MELTS engine just had an internal failure. Capping
      tasks-per-worker bounds how many *subsequent* tasks can silently run
      on that same, possibly-bad state before a fresh process is forced --
      matching the safety margin today's Path.path_multi already relies on
      via its own hard-coded `index[:3]` cap. Raising or disabling this
      cap should only follow an audit of task_fn's own call chain
      confirming every real failure actually raises instead of being
      swallowed -- not assumed.
    """
    completed = 0
    pid = _CTX.current_process().pid
    while True:
        index = task_queue.get()
        if index is _SHUTDOWN:
            return

        results_queue.put((_STARTED, index, pid))
        proxy = _ProgressProxy(results_queue, pid, index)
        try:
            task_fn(proxy, index, **task_kwargs)
            completed += 1
            retiring = max_tasks is not None and completed >= max_tasks
            results_queue.put((_DONE, index, TaskStatus.OK, None, retiring))
            if retiring:
                return
        except Exception as exc:  # noqa: BLE001 -- intentionally broad: see
            # module docstring. Reporting *what* failed is the point --
            # today's `except: idx.append(i); break` in path_multi reports
            # nothing back to the parent at all.
            results_queue.put((_DONE, index, TaskStatus.CRASHED, str(exc), True))
            return


def run_tasks(
    task_fn: Callable,
    indices: Iterable[Any],
    *,
    task_kwargs: Optional[Dict[str, Any]] = None,
    cores: int = 1,
    timeout: float = 240.0,
    stall_timeout: Optional[float] = None,
    max_retries: int = 0,
    kill_grace: float = 5.0,
    max_tasks_per_worker: Optional[int] = 3,
    on_progress: Optional[Callable[[Any], None]] = None,
) -> Dict[Any, TaskResult]:
    """
    Run ``task_fn(q, index, **task_kwargs)`` once per item in `indices`,
    across up to `cores` persistent worker processes, and return one
    TaskResult per index.

    This is the one implementation meant to replace the near-identical
    "spawn `cores` fresh processes per round, wait for the whole round,
    repeat" loop currently duplicated in multi_path, equilibrate_multi,
    findCO2_multi, findLiq_multi, mineral_cosaturation,
    AdiabaticDecompressionMelting and phaseDiagram_calc. Workers are
    started once and each pulls the next index off a shared queue as soon
    as it is free, so one slow task no longer idles the other `cores - 1`
    workers for the rest of a round.

    Parameters
    ----------
    task_fn : callable
        ``task_fn(q, index, **task_kwargs)``. May push its own progress
        items via ``q.put(item)`` while it runs (exactly as path_MELTS
        already does, once per path step) -- see `stall_timeout` and
        `on_progress`. Must raise on failure rather than swallowing its
        own exceptions, so this function can see and report it.
    indices : iterable
        One entry per unit of work (e.g. a row index into a composition
        DataFrame). Order is not preserved in the no-failure case because
        workers race for the next item; `TaskResult.index` says which
        input each result belongs to.
    task_kwargs : dict, optional
        Keyword arguments passed to every call of `task_fn`, in addition
        to `q` and `index`.
    cores : int, default 1
        Upper bound on concurrent workers. Clamped to ``len(indices)`` --
        a 2-item batch never starts more than 2 workers, regardless of
        machine size.
    timeout : float, default 240.0
        Absolute per-task budget, in seconds, measured from when a worker
        actually starts that task (not from when the round began, and not
        derived from any other task's size) -- see module docstring's
        rationale versus today's per-round, smallest-group-sized timeout.
        This is a backstop for "technically still progressing, just slow
        overall" -- `stall_timeout` is almost always what you want for
        catching a genuinely stuck calculation sooner.
    stall_timeout : float or None, default None
        Kill and retry a task if no progress item has arrived for this
        many seconds -- measured from the *last* activity (a fresh task
        starting counts as activity too), not from the task's start. This
        is what actually catches "stuck inside one C call that never
        returns" quickly: path_MELTS reports progress once per path step
        (see MELTS.py's per-step `q.put([run_index, step_i, data])`), so a
        hang inside any single step's calcEquilibriumState() call is
        detectable as soon as `stall_timeout` seconds pass with no new
        step reported -- it does not need to wait for the much larger
        absolute `timeout` to elapse. Left as None (disabled) by default
        because it only works if task_fn actually reports progress
        periodically; a task_fn that only reports once at the very end
        would look "stalled" for its entire runtime.
    max_retries : int, default 0
        How many additional attempts a CRASHED/TIMEOUT/STALLED task gets
        before its TaskResult is finalised as GAVE_UP. Defaults to 0:
        retrying a TIMEOUT or STALLED task (as opposed to a CRASHED one)
        is frequently just paying the same wait again for the same stuck
        calculation, especially before `stall_timeout` is tuned to catch
        it quickly. Raise this once `stall_timeout` is set low enough that
        a retry is actually cheap, or if task_fn's failures are genuinely
        transient (e.g. sensitive to a randomised initial guess).
    kill_grace : float, default 5.0
        Seconds to wait after `terminate()` (SIGTERM) before escalating to
        `kill()` (SIGKILL) and then bounding the final `join()` by the same
        amount -- replaces today's unbounded `p.join()` after `terminate()`,
        which can hang forever if SIGTERM is ignored.
    max_tasks_per_worker : int or None, default 3
        A worker voluntarily exits (and is replaced) after this many
        SUCCESSFUL tasks, not just after a CRASHED/TIMEOUT one. task_fn
        reporting success only means it didn't raise -- it does not mean
        the engine state behind it is trustworthy (see _worker_loop's
        docstring: path_MELTS swallows its own liquidus-search failures
        and returns normally). This bounds how many further tasks can run
        on a silently-bad process before a fresh one is forced, matching
        the safety margin today's Path.path_multi gets from its hard-coded
        `index[:3]` cap. Set to None only once task_fn's whole call chain
        is confirmed to always raise on a genuine engine-level failure.
    on_progress : callable, optional
        Called with each item a task_fn pushes via `q.put(item)` (e.g. a
        single path_MELTS step), unwrapped back to exactly what task_fn
        sent. Lets a caller keep today's "partial results survive a
        timeout" behaviour instead of only getting the final per-task
        outcome.

    Returns
    -------
    dict[index -> TaskResult]
        Exactly one entry per input index, always -- including GAVE_UP
        entries for anything that never succeeded, so a truncated/failed
        run is visible in the return value instead of silently absent
        (contrast with today's multi_path, where del Results[run][-1]
        followed by `if Results[run]:` can make a killed-early run vanish
        from the output with no trace).
    """
    task_kwargs = dict(task_kwargs or {})
    indices = list(indices)
    if not indices:
        return {}
    cores = max(1, min(cores, len(indices)))

    task_queue = _CTX.Queue()
    results_queue = _CTX.Queue()

    attempts: Dict[Any, int] = {i: 1 for i in indices}  # 1 = the initial attempt
    outstanding = set(indices)          # indices without a finalised TaskResult yet
    results: Dict[Any, TaskResult] = {}
    tracked_by_pid: Dict[int, _Tracked] = {}
    workers: List[Any] = []

    def _spawn_worker():
        p = _CTX.Process(
            target=_worker_loop,
            args=(task_fn, task_queue, results_queue, task_kwargs, max_tasks_per_worker),
            daemon=True,
        )
        p.start()
        workers.append(p)
        return p

    def _terminate(process):
        process.terminate()
        process.join(kill_grace)
        if process.is_alive():
            process.kill()
            process.join(kill_grace)

    def _retire(pid):
        """Stop tracking `pid`'s task, drop its process from `workers`, and
        make sure it's actually dead. Does NOT decide whether to spawn a
        replacement -- callers do that, since one is only needed while
        there's still outstanding work."""
        t = tracked_by_pid.pop(pid, None)
        if t is None:
            return
        if t.process in workers:
            workers.remove(t.process)
        # t.process can be None: see the _STARTED handler below for when
        # and why (a worker that dies extremely fast -- before its own
        # _STARTED message is processed -- can have _Tracked created with
        # no process reference). Nothing to terminate in that case; the
        # worker is already known dead, which is exactly why this
        # situation can arise. Defensive guard only -- the _STARTED
        # handler is the actual fix, this just stops it from ever being
        # able to crash the supervisor again if some other path also
        # produces a None reference.
        if t.process is not None and t.process.is_alive():
            _terminate(t.process)

    def _finalise_or_retry(index, status: TaskStatus, error: Optional[str]):
        if index not in outstanding:
            return  # already finalised; ignore a stale/racing report
        if status is not TaskStatus.OK and attempts[index] < 1 + max_retries:
            attempts[index] += 1
            task_queue.put(index)  # goes back in for any idle worker to pick up
            return
        if status is not TaskStatus.OK:
            status = TaskStatus.GAVE_UP
        results[index] = TaskResult(index, status, error, attempts=attempts[index])
        outstanding.discard(index)

    for _ in range(cores):
        _spawn_worker()
    for index in indices:
        task_queue.put(index)

    try:
        while outstanding:
            # 1. Liveness sweep: replace any worker that has died, whether
            #    or not it had a task in flight. A worker can die before it
            #    ever reports _STARTED for anything -- e.g. it fails during
            #    its own process bootstrap (an unpicklable task_fn/argument,
            #    a missing import, an OS-level spawn failure) -- and
            #    without this sweep the supervisor has no way to notice:
            #    there is no timeout to expire for a task that was never
            #    recorded as started, so outstanding work would simply wait
            #    forever for a worker that no longer exists. (First found
            #    by this module's own __main__ smoke test: a task_fn
            #    defined inside `if __name__ == "__main__":` isn't a
            #    picklable module attribute under spawn, so every initial
            #    worker died at this exact bootstrap step.)
            for w in list(workers):
                if not w.is_alive():
                    pid = w.pid
                    workers.remove(w)
                    t = tracked_by_pid.pop(pid, None)
                    if t is not None:
                        _finalise_or_retry(t.index, TaskStatus.CRASHED,
                                            "worker process died unexpectedly")
                    if outstanding:
                        _spawn_worker()

            # 2. Per-task timeout/stall check, every poll regardless of
            #    which worker reports next -- a single shared results_queue
            #    (as opposed to one queue per worker) is what makes
            #    draining every worker's progress "at once" rather than one
            #    worker at a time, today's actual bottleneck, essentially
            #    free. `stall_timeout` (measured since the LAST activity)
            #    is what actually catches a hung calcEquilibriumState()
            #    call quickly; `timeout` (measured since the task started)
            #    is just the backstop for "slow but still moving".
            now = time.time()
            for pid, t in list(tracked_by_pid.items()):
                if stall_timeout is not None and now - t.last_activity > stall_timeout:
                    index = t.index
                    _retire(pid)
                    if outstanding:
                        _spawn_worker()
                    _finalise_or_retry(index, TaskStatus.STALLED,
                                        f"no activity for over {stall_timeout}s")
                elif now - t.started > timeout:
                    index = t.index
                    _retire(pid)
                    if outstanding:
                        _spawn_worker()
                    _finalise_or_retry(index, TaskStatus.TIMEOUT, "timed out")

            # 3. Drain EVERYTHING currently sitting in the queue, not just one
            #    message, before looping back to the liveness sweep above.
            #    This closes (shrinks to effectively zero, not just "one
            #    message at a time") a real race: a worker that finishes and
            #    exits sends its success (_PROGRESS then _DONE) and then
            #    terminates -- if only one message were drained per
            #    iteration, a BURST of near-simultaneous completions (e.g.
            #    several identical-speed workers finishing together) could
            #    leave most of their messages undrained while the NEXT
            #    iteration's liveness sweep already sees those workers as
            #    "not alive" and wrongly reports them CRASHED, discarding a
            #    result that actually succeeded (confirmed live via direct
            #    instrumentation: "LIVENESS SWEEP: pid=... not alive,
            #    finalising index=... as CRASHED" firing immediately after a
            #    real, successful completion). Draining the whole backlog
            #    every iteration means the liveness sweep essentially never
            #    gets to run against a worker whose own messages are still
            #    queued. The first `get` still blocks briefly (so this loop
            #    doesn't busy-spin when there is genuinely nothing to do);
            #    every subsequent one in the same iteration is non-blocking.
            got_any = False
            while True:
                try:
                    item = results_queue.get(timeout=0 if got_any else 0.1)
                except Empty:
                    break
                got_any = True

                tag = item[0]
                if tag == _STARTED:
                    _, index, pid = item
                    proc = next((w for w in workers if w.pid == pid), None)
                    if proc is None:
                        # The worker that sent this has ALREADY died and been
                        # removed from `workers` by an earlier liveness-sweep
                        # pass this same loop -- its own _STARTED message can
                        # arrive late if it crashes extremely fast (e.g. right
                        # after _worker_loop's first task_queue.get()).
                        # Finalise this task now instead of recording a
                        # _Tracked with no process reference (which used to
                        # crash _retire() with AttributeError: 'NoneType'
                        # object has no attribute 'is_alive' the first time
                        # this task was next looked at -- confirmed live, not
                        # hypothetical) and instead of silently sitting until
                        # `timeout` lapses for a process that's already long
                        # gone.
                        _finalise_or_retry(index, TaskStatus.CRASHED,
                                            "worker died before its own start could be registered")
                    else:
                        now = time.time()
                        tracked_by_pid[pid] = _Tracked(index, proc, now, now)
                elif tag == _PROGRESS:
                    _, pid, index, inner_item = item
                    if pid in tracked_by_pid:
                        tracked_by_pid[pid].last_activity = time.time()
                    if on_progress is not None:
                        on_progress(inner_item)
                elif tag == _DONE:
                    _, index, status, error, retiring = item
                    pid = next((p for p, t in tracked_by_pid.items() if t.index == index), None)
                    if retiring:
                        # The worker has already exited on its own -- either
                        # it CRASHED, or it hit max_tasks_per_worker after
                        # this success. Clean it up THIS iteration rather
                        # than waiting for the next liveness-sweep pass: if
                        # this was also the last outstanding task, `while
                        # outstanding:` can exit before that sweep gets
                        # another turn, which is exactly how a worker got
                        # left in `workers` forever -- confirmed via `ps aux`
                        # after a real run (see MAINTENANCE_NOTES.md), not a
                        # hypothetical. Leaving a dead Process object in
                        # `workers` is also a latent PID-reuse hazard: a
                        # later _STARTED lookup (`next(w for w in workers if
                        # w.pid == pid)`) could match that stale,
                        # already-dead entry instead of the new process the
                        # OS just gave the same pid to.
                        if pid is not None:
                            _retire(pid)
                            if outstanding:
                                _spawn_worker()
                    else:
                        # Still alive and will pull its next task on its own
                        # -- only stop tracking this *task*, don't touch the
                        # worker.
                        tracked_by_pid.pop(pid, None)
                    _finalise_or_retry(index, status, error)

        return results
    finally:
        # Guaranteed cleanup: whatever happens above -- including the
        # caller hitting Ctrl-C, or an exception raised out of this loop --
        # every worker this call started gets a shutdown signal and is
        # terminated/joined. Nothing is left running with its timeout no
        # longer enforced, which is what happens today if the parent dies
        # mid-round.
        for _ in workers:
            task_queue.put(_SHUTDOWN)
        for w in list(workers):
            if w.is_alive():
                _terminate(w)
            else:
                w.join()


# ---------------------------------------------------------------------------
# multi_path, sketched against run_tasks (illustrative only -- Path.py is
# NOT changed, and this block is not imported or executed by anything).
# ---------------------------------------------------------------------------
#
# Today, Path.multi_path's batch branch spawns `cores` fresh processes per
# round running Path.path_multi (which itself loops over up to 3 indices
# before exiting), repeating rounds until every index has reported. With
# run_tasks, the per-index body of that loop becomes the task_fn -- closer
# to today's single-run Path.path() than to path_multi, since looping across
# many indices is now run_tasks's job:
#
#     def _path_task(q, i, *, Model, comp, T_path_C, P_path_bar, ...):
#         # exactly the body already inside path_multi's `for i in index:`,
#         # minus the addNodeAfter() reuse (new process per task => always a
#         # fresh MELTSdynamic(), so the "does addNodeAfter leak state after
#         # a failure" question from the plan doesn't need answering here).
#         Results, tr = path_MELTS(Model=Model, comp=comp, T_path_C=T_path_C[i],
#                                   P_path_bar=P_path_bar[i], ..., q=q, index=i)
#         if tr is False:
#             raise RuntimeError("path_MELTS reported trial failure")
#
#     results = run_tasks(
#         _path_task, range(len(T_path_C)),
#         task_kwargs=dict(Model=Model, comp=comp, T_path_C=T_path_C, P_path_bar=P_path_bar, ...),
#         cores=cores, timeout=timeout_main, max_retries=2,
#         on_progress=lambda item: Results.setdefault(item[0], {}).__setitem__(item[1], item[2]),
#     )
#     # results[i].status is TaskStatus.OK / TIMEOUT / CRASHED / GAVE_UP for
#     # every i in range(len(T_path_C)) -- including the ones that never
#     # succeeded, instead of today's silent omission from the output.


def _toy_task(q, i, *, fail_indices=(), hang_indices=(), steady_indices=()):
    """
    A deliberately trivial task_fn for the smoke test below. Must be a
    module-level function (not defined inside `if __name__ == "__main__":`,
    where the first version of this test accidentally put it) -- under the
    `spawn` context every worker re-imports this module and looks up
    `task_fn` by its qualified name, so it has to already exist as a normal
    module attribute. This is the same constraint today's Path.path_multi /
    Path.path already satisfy by being ordinary top-level functions in
    Path.py.
    """
    import random
    if i in steady_indices:
        # Several genuine progress reports, each gap well under
        # stall_timeout -- proves steady (if individually slow) progress
        # is never mistaken for a stall, even once the *cumulative* time
        # exceeds stall_timeout.
        for step in range(5):
            time.sleep(0.15)
            q.put(("progress", i, step))
        return
    if i in hang_indices:
        q.put(("progress", i, "started"))  # one report, then total silence
        time.sleep(5)  # >> stall_timeout, but < the overall `timeout` below
    if i in fail_indices:
        raise ValueError(f"synthetic failure for task {i}")
    q.put(("progress", i, "halfway"))
    time.sleep(random.uniform(0.01, 0.05))


def _quick_task(q, i, *, fail_indices=()):
    """Module-level (see _toy_task's docstring for why) task_fn for the
    high-churn zombie-leak regression test."""
    if i in fail_indices:
        raise ValueError(f"synthetic failure for task {i}")


def _flaky_task(q, i, *, fail_frac=0.3, hang_frac=0.2):
    """
    Module-level task_fn for the short-timeout/high-churn race regression
    test: a deterministic (seeded by `i`) mix of near-instant successes,
    near-instant crashes, and hangs that will always exceed the short
    timeout used below. Combined with a short `timeout` and fast
    `kill_grace`, this maximises worker spawn/kill/respawn churn -- the
    condition that produced AttributeError: 'NoneType' object has no
    attribute 'is_alive' live (a worker dying fast enough that the
    liveness sweep removes it from `workers` before its own _STARTED
    message is processed; see run_tasks's _STARTED handling).
    """
    import random
    import time
    r = random.Random(i)
    x = r.random()
    if x < hang_frac:
        time.sleep(30)
    elif x < hang_frac + fail_frac:
        raise ValueError(f"synthetic failure for task {i}")


if __name__ == "__main__":
    # Toy, MELTS-free smoke test of the pool-management logic itself:
    # plain successes, a task that always raises (max_retries=0 -> GAVE_UP
    # on the first failure, no wasted retry), a task that goes silent after
    # one report (caught by stall_timeout, much sooner than the 10s overall
    # `timeout`), and a task that makes slow-but-steady progress (never
    # mistaken for a stall). Run directly: `python -m petthermotools.mp_utils`.
    progress_seen = []
    start = time.time()
    out = run_tasks(
        _toy_task,
        range(9),
        task_kwargs=dict(fail_indices={3}, hang_indices={5}, steady_indices={8}),
        cores=4,
        timeout=10.0,
        stall_timeout=0.5,
        max_retries=0,
        on_progress=progress_seen.append,
    )
    elapsed = time.time() - start

    for i in range(9):
        print(i, out[i].status.value, out[i].attempts, out[i].error)
    print(f"elapsed: {elapsed:.1f}s")

    assert out[3].status is TaskStatus.GAVE_UP and out[3].attempts == 1
    assert "synthetic failure" in out[3].error

    assert out[5].status is TaskStatus.GAVE_UP and out[5].attempts == 1
    assert "no activity" in out[5].error  # killed by stall_timeout, not the 10s timeout

    for i in (0, 1, 2, 4, 6, 7, 8):
        assert out[i].status is TaskStatus.OK, (i, out[i])

    assert len(progress_seen) >= 6 + 5  # >=1 per quick success, 5 for the steady task

    # The whole run must finish in well under the 10s overall `timeout` --
    # proves task 5 was actually killed by stall_timeout (~0.5s after its
    # one report), not by waiting out the much larger absolute budget.
    assert elapsed < 5.0, f"took {elapsed:.1f}s -- stall_timeout does not seem to be firing"

    print("mp_utils smoke test: OK")

    # --- High-churn test for the "voluntary retirement leaks a worker"
    # bug: with max_tasks_per_worker=1, EVERY successful task retires its
    # worker, which used to be cleaned up only by the liveness sweep's
    # *next* pass -- a pass that may never come if that was also the last
    # outstanding task. Found for real, not hypothesised: a live
    # findLiq_multi run (named findLiq_multi_v2 at the time; later promoted
    # to replace the old findLiq_multi, see MAINTENANCE_NOTES.md 2026-10-10)
    # left several worker processes running (and
    # burning ~80% CPU) for 30-50+ minutes after the parent script had
    # already returned (see MAINTENANCE_NOTES.md). cores=8 with 60 tasks
    # forces dozens of spawn/retire cycles, maximising the chance of
    # catching a leak if the fix regresses.
    import multiprocessing as _mp

    fail_set = set(range(0, 60, 7))  # a scattered mix of crashes and successes
    out2 = run_tasks(
        _quick_task,
        range(60),
        task_kwargs=dict(fail_indices=fail_set),
        cores=8,
        timeout=5.0,
        max_retries=0,
        max_tasks_per_worker=1,
    )
    leaked = _mp.active_children()
    print(f"\nhigh-churn test: {len(out2)} tasks, {len(fail_set)} designed to fail, "
          f"{len(leaked)} worker process(es) still alive after run_tasks returned.")
    assert len(leaked) == 0, f"LEAK: {leaked}"
    for i in fail_set:
        assert out2[i].status is TaskStatus.GAVE_UP, (i, out2[i])
    for i in set(range(60)) - fail_set:
        assert out2[i].status is TaskStatus.OK, (i, out2[i])
    print("high-churn zombie-leak test: OK")

    # --- Short-timeout/high-churn race regression test: see _flaky_task's
    # docstring. Repeated several times since the race that produced
    # AttributeError: 'NoneType' object has no attribute 'is_alive' (found
    # live via a real findLiq_multi run against a real dataset, 10-30s
    # timeouts) is timing-dependent, not guaranteed to reproduce every
    # single run -- the point is many independent chances, none of which
    # should ever crash the supervisor itself (individual task failures,
    # i.e. GAVE_UP, are expected and fine).
    for trial in range(5):
        out3 = run_tasks(
            _flaky_task,
            range(80),
            task_kwargs=dict(fail_frac=0.3, hang_frac=0.2),
            cores=8,
            timeout=2.0,
            max_retries=1,
            max_tasks_per_worker=5,
            kill_grace=1.0,
        )
        assert len(out3) == 80
    leaked3 = _mp.active_children()
    assert len(leaked3) == 0, f"LEAK: {leaked3}"
    print("short-timeout/high-churn race regression test (5 trials, 400 tasks total): OK")
