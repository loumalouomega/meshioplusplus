"""Classify why a notebook failed to execute, from the log nbconvert wrote.

Used by ``.github/workflows/notebooks.yml`` (roadmap 3.11.4). A notebook that
needs a library or a build option the CI wheel does not have (CutCells,
OpenCASCADE, Netgen, MMG, an optional Python package) is *skipped*: it says so in
the job summary but does not fail the job. Anything else is a *fail*: an example
that no longer runs against the library.

    python tools/classify_notebook_error.py LOG     # prints "skip|fail<TAB>message"
"""

from __future__ import annotations

import re
import sys

_ANSI = re.compile(r"\x1b\[[0-9;]*m")
_ERROR_LINE = re.compile(r"^(?:[A-Za-z_][\w.]*)(?:Error|Exception|Exit): ?.*$")
_MISSING = re.compile(
    r"this build has no|needs the OpenCASCADE|is not installed|"
    r"ModuleNotFoundError|ImportError|needs [A-Za-z+ ]+, which is not|"
    r"-DMESHIOPLUSPLUS_WITH_|"
    # an external tool the operation drives (QuadriFlow, Instant Meshes, ...)
    r"no [\w ]+ is installed|install \S+ .* on PATH",
    re.I,
)


def last_error(text: str) -> str:
    """The last ``SomeError: message`` line of an nbconvert log, or ``""``."""
    found = ""
    runner = ""  # nbclient's own errors (a timeout, a dead kernel): the fallback
    for line in _ANSI.sub("", text).splitlines():
        line = line.strip()
        # nbconvert prefixes log lines with a timestamp in CI; keep what follows it.
        line = re.sub(r"^\d{4}-\d\d-\d\dT[\d:.]+Z ", "", line)
        if not _ERROR_LINE.match(line):
            continue
        if "nbclient.exceptions" in line:
            runner = line
        else:
            found = line
    return found or runner


def classify(text: str) -> tuple[str, str]:
    message = last_error(text) or "no error line found in the log"
    return ("skip" if _MISSING.search(message) else "fail"), message


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    with open(argv[0], encoding="utf-8", errors="replace") as f:
        kind, message = classify(f.read())
    print(f"{kind}\t{message[:300]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
