"""
Tier A test for the harness itself: reading back a child's output must work
whatever the machine's default encoding is.

Found on Windows: a child's stdout follows the console/locale encoding (cp1252
on a default shell, UTF-8 if PYTHONIOENCODING/PYTHONUTF8 happens to be set),
and the harness read files back with the locale default. A character
cp1252 cannot encode -- the mineral-chemistry ones (Fe³⁺, ρ, ...) are exactly
what this package's output could one day contain -- then crashed the child, or
made the parent's read raise UnicodeDecodeError (a ValueError, not the OSError
the harness caught).
"""
import subprocess
import sys

from _proc_helpers import read_utf8, utf8_env

# Characters cp1252 cannot represent (plus one it can, the em dash used in
# the package's own "Timeout — terminating process" message).
_TEXT = "Fe³⁺ ρ — °C"
# The child's code is pure ASCII (escapes, not literal characters), so how
# the command line itself is encoded cannot affect the result.
_CODE = r"print('Fe³⁺ ρ — °C', flush=True)"


def test_non_cp1252_child_output_round_trips(tmp_path):
    out_path = tmp_path / "child.out"
    with open(out_path, "w") as out:
        proc = subprocess.run(
            [sys.executable, "-c", _CODE],
            stdout=out, stderr=subprocess.PIPE, env=utf8_env(), timeout=60,
        )
    assert proc.returncode == 0, proc.stderr.decode("utf-8", "replace")
    assert read_utf8(out_path).strip() == _TEXT
