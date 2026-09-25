# Tests

## Layout

| Path | Tier | Needs | Runs in CI |
|---|---|---|---|
| `unit/` | A | nothing beyond the package's own dependencies | yes (ubuntu / macos / windows x Python 3.10-3.13, plus one pandas<3 job) |
| `integration/alphamelts/` | B1 | a working alphaMELTS (`meltsdynamic`) | no: skipped automatically when alphaMELTS is not found |
| `integration/magemin/` | B2 (not written yet) | Julia + MAGEMinCalc | no |
| `golden/alphamelts/` | data | recorded results for the B1 tests | - |
| `api/` | data | snapshot of the public API (`public_api.json`) used by `unit/test_public_api.py`, and the script that regenerates it | - |

`conftest.py` puts `src/` on `sys.path` (tests run against the working tree, not an installed copy), registers the
`melts` / `magemin` markers and skips those tests when the backend is missing. Availability is checked with
`importlib.util.find_spec`, never by importing: importing `juliacall` would make `Path.py` refuse to run MELTS on
Windows (it treats that as "MAGEMin already initiated").

Older files kept as they were: `test_MELTS.py`, `test_ngibbs_geometry.py`, `SCSS_unit_Tests.ipynb`, and the
vendored `MELTS/` folder (git-ignored).

## Running

```
python -m pytest Unit_tests/unit                 # Tier A, no alphaMELTS needed (about 2 min)
python -m pytest Unit_tests/integration          # real MELTS golden tests (about 45 s)
python -m pytest Unit_tests/unit Unit_tests/integration -q
```

## An xfail here means "documented bug", not "flaky"

Most of `unit/` is a set of failure-injection tests that pin bugs found in the multiprocessing supervisor
(`Path.py`) and in input handling. Each is `xfail(raises=AssertionError, strict=False)`:
- it fails today, for the reason written in its `reason=`;
- `raises=AssertionError` means any other exception (a harness problem, an import error) is a real failure, not an
  expected one; harness problems use `pytest.fail()` for the same reason;
- when the bug is fixed the test starts to pass (XPASS): remove the marker and it becomes a hard assertion.

| File | Bug pinned |
|---|---|
| `test_multi_path_supervisor.py` | a worker that dies before its first queue message is respawned forever |
| `test_missing_main_guard.py` | same failure, triggered by a user script without `if __name__ == "__main__":` |
| `test_timeout_handling.py` | runs killed by a timeout are lost or returned as complete with no error or warning; a worker that hangs before reporting loops forever; a worker that ignores SIGTERM blocks `join()` (POSIX only) |
| `test_process_reaping.py` | workers that exit on their own are never `join()`ed |
| `test_worker_cleanup_on_failure.py` | if the parent is interrupted or raises, live workers are left running and their timeout is no longer enforced |
| `test_input_validation.py` | unknown or mistyped `Model`, and bad compositions, are silently routed, accepted, or fail with obscure errors (the file also pins routing and existing clear errors that must not regress) |
| `test_progress_bars.py` | FIXED in Phase 1 (now an ordinary passing regression test): `tqdm.notebook` needed ipywidgets, which `setup.py` does not declare, so `findLiq_multi` etc. crashed on a clean install outside Jupyter; the modules now use `tqdm.auto` |
| `test_import_hygiene.py` | `petthermotools.Path` is shadowed by `pathlib.Path` (star-imports); pins current behaviour |

## Source hygiene guard (`unit/test_no_invalid_escape_sequences.py`)

Fails if any package source file contains an invalid escape sequence (for example `'$\degree$'` without an `r` prefix).
These warn on every import (DeprecationWarning up to Python 3.11, SyntaxWarning from 3.12) and will become errors. Write
such strings as raw strings, `r'...'`; if the string also needs a real escape such as `\n`, double only the offending
backslash. All 14 existing ones were fixed in Phase 1, each proven to leave the string's value unchanged.

## Public API freeze (`unit/test_public_api.py`)

A guard rail for maintainers, not a restriction on users: nothing here ships with the package or runs at run time, and
it fails only when a change would break something users call.

- **Documented usage.** Every *maintained* notebook is scanned; each `ptt.<name>` it uses must still exist and each call
  must still fit the function's signature (keywords and positional count). Maintained = the notebooks ReadTheDocs links
  from the toctrees in `docs/index.rst`, plus everything in `Workshops/`, `docs/teaching_materials/` and
  `docs/PaperFigures/`. Notebooks listed in `exclude_patterns` in `docs/conf.py` are retired and ignored. ReadTheDocs
  never runs notebooks (`nbsphinx_execute = 'never'`), so this is the only automatic check that the published examples
  still match the code. A new docs notebook is covered as soon as it is in a toctree.
- **Snapshot** (`api/public_api.json`): every public name and, for functions and classes defined in the package, their
  parameters. It fails on: a removed or renamed name or parameter, reordered leading positional parameters, a parameter
  that used to be optional and is now required, or a new required parameter. Free: new names, new optional parameters,
  changes to behaviour or default values, and underscore names. Names that only leak in from other libraries
  (`ptt.np`, ...) are recorded but not enforced (`ENFORCE_EXTERNAL_NAMES` in `unit/_api_scan.py`); names from optional
  packages (alphaMELTS, ...) are never required, so a minimal environment such as CI passes.
- **Changing the public API on purpose:** run `python Unit_tests/api/capture_public_api.py` and commit the diff; the diff
  is the review of what changed for users. Deprecate first (keep the old name and warn) before removing anything.
  The script imports the package from the **committed** code (a temporary checkout of `HEAD`), never from your working
  tree, so uncommitted work cannot leak into the snapshot (an uncommitted `trace_engine` import once put five names in
  it and CI failed): commit the code change first, then capture. `--working-tree` overrides this.

## Golden tests (`integration/alphamelts/`)

Eleven small real-MELTS calculations (isobaric, fractional, pMELTS, v1.2.0 with fO2 buffer, polybaric, isochoric,
isothermal/isentropic decompression, a 6-pressure batch, `equilibrate`, `findLiq`) are compared with recorded results.
They are the safety net for refactors: a change that alters results fails here. Two runs of a case on one machine are
bit-for-bit identical, for the single-run and the batch path; the macOS-recorded data also reproduces on Windows at
rtol 1e-6 under pandas 2 and 3.

- Cases live in `_cases.py` (used by the tests, the capture script and the benchmarks).
- `capture_golden.py` rewrites the reference data. Run it only when a change in results is deliberate (for example a
  constant is corrected) and commit the diff. Never run it casually or from a second machine.
- Tolerance is rtol 1e-6; set `PTT_GOLDEN_RTOL` to change it. The set of tables (phases) must match exactly.
  For a refactor that must not change any result, run locally with `PTT_GOLDEN_RTOL=1e-11`: the CSVs store 12
  significant digits, so that is the practical exact gate (`0` can never pass, because the stored values are rounded).
- Deliberate result changes are re-recorded with `capture_golden.py` in the same PR (example: correcting the Fe2O3
  molar mass to 159.69 everywhere moved primary results by a median 5e-6, 99th percentile 8e-4, with the same phases).
- MELTS drops `*_tbl.txt` / `.inp` files in the current directory: the tests run in a temporary directory. This is
  also why `docs/Examples/**` has stray `*_tbl.txt` files after a notebook is run there. The six file names the
  engine writes (`Bulk_comp_tbl.txt`, `Liquid_comp_tbl.txt`, `Phase_main_tbl.txt`, `Solid_comp_tbl.txt`,
  `System_main_tbl.txt`, `liquid-model-batch.inp`) are git-ignored and no longer tracked, so running a notebook does
  not change `git status`.
- The engine writes about 420 KB per run to **stderr** (file descriptor 2). `_cases.quiet_fds()` silences it;
  use `keep_stderr=False` when timing (capturing it to one shared file serialises the workers).

## Other machines (especially Windows)

- A virtual environment does not see alphaMELTS unless the folder containing `meltsdynamic.py` is on the path (the
  `my_MELTS_path.pth` written by `install_alphaMELTS()` lives in the base environment). Set `PYTHONPATH` to that folder.
- Run from a copy outside cloud-synced folders (Google Drive): it avoids `__pycache__`/`.pytest_cache` and DLL reads in
  the synced folder, and Drive can briefly present files as 0-byte placeholders.
- In a venv, `sys.executable` on Windows is a launcher: the process tree is launcher -> interpreter -> workers, and the
  workers run under the base interpreter. Match processes by parent PID, not by the venv path.
- `test_worker_ignoring_sigterm_does_not_block_forever` is skipped on Windows (`terminate()` cannot be ignored there).
- Expected on Windows without ipywidgets: the `findliq_v102` golden test fails with the ImportError above until the
  progress-bar bug is fixed.

## Writing multiprocessing or failure tests here (lessons learned)

- Run each hang scenario in its **own interpreter** and kill the whole process tree on a time cap
  (`_proc_helpers.run_script_capped`, `kill_process_tree`). Running it on a thread inside pytest leaks a respawn loop
  into every later test.
- Start a scenario's clock after its imports (a start marker), so slow imports on a busy machine cannot look like a hang.
- Redirect a child's output to **files, not pipes** (surviving workers inherit the pipe and block `communicate()`), and
  read it back as UTF-8 with `utf8_env()` / `read_utf8()`: the default encoding differs by machine.
- Helper scripts need an `if __name__ == "__main__":` guard (spawn re-imports them), except the one that tests its absence.
- Do not let a test leave half-built objects with failing `__del__`, or a live exception traceback, behind. A leaked
  broken tqdm bar crashed a whole pytest session with `SystemError: AST constructor recursion depth mismatch` on Python
  3.11.0 (INTERNALERROR), at the next failure report.
- Do not name a folder `melts/` anywhere in the repo: `.gitignore` has `MELTS/`, and git's ignore matching is
  case-insensitive on macOS and Windows (not Linux), so the folder would be silently ignored. Hence `alphamelts/`.

See also `benchmarks/README.md`.
