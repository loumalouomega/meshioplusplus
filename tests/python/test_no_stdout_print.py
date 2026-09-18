"""Guards roadmap §1's "reader fallback prints to stdout" fix from regressing.

``_helpers.py``'s ambiguous-extension fallback loop used to ``print(e)`` a
declined reader's exception straight to stdout, corrupting any CLI subcommand
that writes JSON there. The fix (routing declines through ``logging`` instead)
leaves nothing to test the *absence* of except by re-introducing the bug, so
this is a static guard: no bare, unredirected ``print(`` call may exist under
the library package outside the two surfaces that are allowed to talk to
stdout on purpose (the CLI and the MCP HTTP dashboard).

A ``rich`` ``Console(stderr=True).print(...)`` call is a *method* call (it
follows a ``.``), not the builtin, and is unaffected by this check; likewise
``print(..., file=sys.stderr)`` is an explicit stderr write, not a stdout one.
"""

import pathlib
import re

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
PKG = REPO / "src" / "python" / "meshioplusplus"

# print( not preceded by a '.' (a rich Console.print(...) method call).
_BARE_PRINT = re.compile(r"(?<!\.)\bprint\(")

_ALLOWED_DIRS = ("_cli", "mcp")


def _iter_python_files():
    if not PKG.is_dir():
        return
    for path in PKG.rglob("*.py"):
        rel = path.relative_to(PKG)
        if rel.parts and rel.parts[0] in _ALLOWED_DIRS:
            continue
        yield path


def test_no_bare_stdout_print_outside_cli_and_mcp():
    violations = []
    for path in _iter_python_files():
        text = path.read_text()
        for lineno, line in enumerate(text.splitlines(), start=1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            for m in _BARE_PRINT.finditer(line):
                # Explicit stderr redirection is fine (physicsnemo/train.py's
                # error-reporting path uses this).
                tail = line[m.end() :]
                if "file=sys.stderr" in tail or "file=stderr" in tail:
                    continue
                violations.append(f"{path.relative_to(REPO)}:{lineno}: {stripped}")

    assert not violations, (
        "bare print() writes to stdout, which corrupts any CLI subcommand "
        "emitting JSON there; use `logging`/`warnings` (see _helpers.py's "
        "ambiguous-extension fallback for the pattern) or an explicit "
        "file=sys.stderr instead:\n" + "\n".join(violations)
    )
