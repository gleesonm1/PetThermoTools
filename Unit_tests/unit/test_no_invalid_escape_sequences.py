"""
Tier A: the package source must not contain invalid escape sequences.

A string like 'T ($\\degree$C)' written without the doubled backslash, or without an
r prefix, is an "invalid escape sequence". Python keeps the backslash today, so it
works, but it warns on every import (DeprecationWarning up to 3.11, SyntaxWarning
from 3.12) and is slated to become an error. It showed up as a wall of warnings on
Windows and in CI. All 14 occurrences were fixed in Phase 1 (each fix proven to leave
the string's value unchanged); this keeps them from coming back.

Fix a new one by writing the string as a raw string, r'...', for LaTeX labels,
regular expressions and Windows paths. If the string also needs a real escape such as
\\n, double just the offending backslash instead.
"""
import warnings
from pathlib import Path

SRC = Path(__file__).resolve().parents[2] / "src" / "petthermotools"


def invalid_escapes(source, name="<source>"):
    """(line, message) for every invalid escape sequence Python reports when compiling `source`."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        compile(source, name, "exec")
    return [(w.lineno, str(w.message)) for w in caught if "invalid escape sequence" in str(w.message)]


def test_the_package_source_has_no_invalid_escape_sequences():
    problems = []
    files = sorted(SRC.rglob("*.py"))
    assert len(files) >= 15, "found almost no package source files: the path is wrong"
    for path in files:
        for line, message in invalid_escapes(path.read_text(encoding="utf-8"), str(path)):
            problems.append(f"{path.relative_to(SRC).as_posix()}:{line}: {message}")
    assert not problems, (
        f"{len(problems)} invalid escape sequence(s); write the string as r'...' (or double the backslash):\n  "
        + "\n  ".join(problems)
    )


def test_the_checker_itself_detects_a_bad_escape():
    """Guards against the test above passing only because the warning is not being recorded."""
    assert invalid_escapes("label = 'T ($\\degree$C)'\n")            # '\d' is invalid
    assert not invalid_escapes("label = r'T ($\\degree$C)'\n")       # raw string: fine
    assert not invalid_escapes("text = 'line one\\nline two'\n")     # a real escape: fine
