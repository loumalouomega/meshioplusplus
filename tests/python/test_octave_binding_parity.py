"""Keeps the Octave/MATLAB binding in step with the R binding.

``bindings/octave/meshioplusplus`` ports the R binding entry point for entry
point: every ``R_mio_x`` .Call target has an ``O_mio_x`` MEX command with the
same argument count, and every exported ``mio_x`` R function has an
``inst/+mio/x.m`` wrapper. These static checks fail a change that extends the
R binding (or the gateway) without the other side, and a ``.m`` file that uses
Octave-only syntax MATLAB would reject. They need neither Octave nor R.
"""

import pathlib
import re

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
R_PKG = REPO / "bindings" / "r" / "meshioplusplus"
OCT_PKG = REPO / "bindings" / "octave" / "meshioplusplus"
MIO_DIR = OCT_PKG / "inst" / "+mio"

# Commands the gateway implements itself (src/mio_mex.c), not ports of R.
_BUILTIN_COMMANDS = {"__commands__", "__release__", "__is_open__"}


def _r_exports():
    text = (R_PKG / "NAMESPACE").read_text()
    return sorted(set(re.findall(r"^export\((mio_[A-Za-z0-9_]+)\)", text, re.M)))


def _r_calldefs():
    text = (R_PKG / "src" / "init.c").read_text()
    return {
        name[2:]: int(n)
        for name, n in re.findall(r"CALLDEF\((R_mio_\w+),\s*(\d+)\)", text)
    }


def _octave_entries():
    entries = {}
    for path in sorted((OCT_PKG / "src").glob("*.c")):
        for name, n in re.findall(
            r"MIO_O_ENTRY\((mio_\w+),\s*(\d+)\)", path.read_text()
        ):
            assert name not in entries, f"{name} registered twice"
            entries[name] = int(n)
    return entries


def _m_files():
    return sorted(OCT_PKG.glob("inst/**/*.m")) + sorted(OCT_PKG.glob("tests/*.m"))


def test_every_r_export_has_an_octave_wrapper():
    missing = [
        e for e in _r_exports() if not (MIO_DIR / f"{e[len('mio_'):]}.m").exists()
    ]
    assert not missing, f"R exports without an inst/+mio/<name>.m wrapper: {missing}"


def test_every_r_entry_point_has_a_mex_command_with_the_same_arity():
    r_defs = _r_calldefs()
    oct_entries = _octave_entries()
    missing = sorted(set(r_defs) - set(oct_entries))
    assert not missing, f"R .Call entry points without a MEX command: {missing}"
    wrong = {
        k: (r_defs[k], oct_entries[k]) for k in r_defs if oct_entries[k] != r_defs[k]
    }
    assert not wrong, f"argument counts differ (R, Octave): {wrong}"


def test_every_gateway_command_used_by_a_wrapper_exists():
    known = set(_octave_entries()) | _BUILTIN_COMMANDS
    unknown = {}
    for path in _m_files():
        for cmd in re.findall(r"meshioplusplus_mex\(\s*'([^']+)'", path.read_text()):
            if cmd not in known:
                unknown.setdefault(cmd, []).append(path.name)
    assert not unknown, f"wrappers call unregistered MEX commands: {unknown}"


def test_no_case_insensitive_name_collisions():
    # A class file (MeshHandle.m) and a function file must not differ only
    # by case: macOS and Windows file systems would merge them.
    seen = {}
    for path in sorted(MIO_DIR.glob("*.m")):
        key = path.stem.lower()
        assert key not in seen, f"{path.name} collides with {seen[key]}"
        seen[key] = path.name


# Octave accepts these; MATLAB does not.
_OCTAVE_ONLY = [
    (re.compile(r"^\s*#"), "'#' comment"),
    (re.compile(r"!="), "'!=' (use '~=')"),
    (
        re.compile(r"\bend(function|if|for|while|switch|_try_catch|_unwind_protect)\b"),
        "Octave end keyword",
    ),
    (re.compile(r"\+\+|--(?=\s*;)|[+\-*/]="), "increment/compound assignment"),
    (re.compile(r"(^|[^\w.])printf\("), "printf (use fprintf)"),
    (re.compile(r"(^|[^\w'])!\w"), "'!' negation (use '~')"),
    (re.compile(r"\bunwind_protect\b"), "unwind_protect (use onCleanup)"),
]


def _code_part(line):
    """The line without its trailing comment, with string literals blanked."""
    out, in_str, i = [], False, 0
    while i < len(line):
        c = line[i]
        if in_str:
            if c == "'":
                if i + 1 < len(line) and line[i + 1] == "'":
                    i += 2
                    continue
                in_str = False
            i += 1
            continue
        if c == "%":
            break
        if c == "'":
            prev = line[:i].rstrip()[-1:] if line[:i].strip() else ""
            # A quote after an identifier, ')' ']' '}' or '.' is a transpose.
            if prev and (prev.isalnum() or prev in ")]}_.'"):
                out.append(c)
            else:
                in_str = True
                out.append("''")
            i += 1
            continue
        out.append(c)
        i += 1
    return "".join(out)


@pytest.mark.parametrize("path", _m_files(), ids=lambda p: str(p.relative_to(OCT_PKG)))
def test_m_files_use_shared_matlab_octave_syntax(path):
    problems = []
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        if re.match(r"^\s*#", line):
            problems.append(f"{lineno}: '#' comment")
            continue
        code = _code_part(line)
        if '"' in code:
            problems.append(f"{lineno}: double-quoted string (use single quotes)")
        for pattern, what in _OCTAVE_ONLY[1:]:
            if pattern.search(code):
                problems.append(f"{lineno}: {what}")
    assert not problems, f"{path.name}:\n" + "\n".join(problems)
