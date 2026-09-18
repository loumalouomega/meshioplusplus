"""Guards roadmap §1's locale item from regressing.

``std::strtod``/``std::atof``/``std::stod``/``std::stof`` and a ``%f``/``%e``/
``%g`` ``std::snprintf``/``std::sprintf`` conversion all honour the process's
``LC_NUMERIC`` category, which is exactly the bug ``detail/fast_number.hpp``
(``parse_double``/``snprintf_c``) exists to route around -- see that header's
own doc comment. Every call site under ``src/cpp/`` was migrated to it; this
is a locale-free static guard (no comma-decimal locale needs to be installed
for it to run, unlike the header's own gtest suite) so a new call site added
later cannot silently reintroduce the bug.

Integer parsing (``strtol``/``strtoll``) is deliberately not guarded: the C
locale's ``LC_NUMERIC`` affects only the decimal point, which an integer
literal doesn't have -- see ``fast_number.hpp``'s own doc comment for why
that's out of scope.
"""

import pathlib
import re

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
SRC = REPO / "src" / "cpp"

# The header itself is the one place these names are allowed to appear (its
# own implementation, and doc-comment prose explaining the tiered fallback).
_ALLOWED_FILE = SRC / "include" / "meshioplusplus" / "detail" / "fast_number.hpp"

_FLOAT_FORMAT_SPEC = re.compile(r"%[-0-9.#+ *]*[eEfgG]")
_BANNED_CALL = re.compile(r"\b(strtod|atof|stod|stof)\s*\(")
_BANNED_PRINTF = re.compile(r"\b(snprintf|sprintf)\s*\(")


def _iter_cpp_files():
    if not SRC.is_dir():
        return
    for path in SRC.rglob("*.[hc]pp"):
        if "third_party" in path.parts:
            continue
        if path == _ALLOWED_FILE:
            continue
        yield path


def test_no_locale_sensitive_float_parse_outside_fast_number_hpp():
    violations = []
    for path in _iter_cpp_files():
        for lineno, line in enumerate(path.read_text().splitlines(), start=1):
            stripped = line.strip()
            if stripped.startswith("//") or stripped.startswith("*"):
                continue
            for m in _BANNED_CALL.finditer(line):
                # A qualified detail::parse_double-style call, or a call on
                # some other object (e.g. a variable named "stod"), is fine;
                # only a bare/std::-qualified strtod/atof/stod/stof is banned.
                prefix = line[: m.start()]
                if (
                    prefix.endswith(".")
                    or prefix.endswith("::detail::")
                    or "parse_double" in line
                ):
                    continue
                violations.append(f"{path.relative_to(REPO)}:{lineno}: {stripped}")

    assert not violations, (
        "locale-sensitive float parsing outside detail/fast_number.hpp; use "
        "detail::parse_double instead (see that header for why):\n"
        + "\n".join(violations)
    )


def test_no_locale_sensitive_float_format_outside_fast_number_hpp():
    violations = []
    for path in _iter_cpp_files():
        for lineno, line in enumerate(path.read_text().splitlines(), start=1):
            stripped = line.strip()
            if stripped.startswith("//") or stripped.startswith("*"):
                continue
            for m in _BANNED_PRINTF.finditer(line):
                if not _FLOAT_FORMAT_SPEC.search(line):
                    continue  # a pure-integer/string format has nothing to repair
                prefix = line[: m.start()]
                if prefix.endswith(".") or "snprintf_c" in line:
                    continue
                violations.append(f"{path.relative_to(REPO)}:{lineno}: {stripped}")

    assert not violations, (
        "locale-sensitive float formatting outside detail/fast_number.hpp; "
        "use detail::snprintf_c instead (see that header for why):\n"
        + "\n".join(violations)
    )
