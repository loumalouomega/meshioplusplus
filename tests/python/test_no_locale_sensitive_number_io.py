"""Guards locale-independent number I/O in the C++ core from regressing.

Two independent exposures, two headers:

* **The C locale** (``setlocale`` -- a Qt app's ``setlocale(LC_ALL, "")``, or
  Python's ``locale.setlocale(locale.LC_ALL, "")``). ``std::strtod``/``std::atof``/
  ``std::stod``/``std::stof`` and a ``%f``/``%e``/``%g`` ``std::snprintf``/
  ``std::sprintf`` conversion all honour ``LC_NUMERIC``, which is exactly what
  ``detail/fast_number.hpp`` (``parse_double``/``snprintf_c``) routes around --
  see that header's doc comment.
* **The C++ locale** (``std::locale::global``). A stream takes the global locale
  when it is constructed, so ``iss >> x`` misreads ``1.5`` and, worse, ``os << n``
  writes an *integer* with the locale's digit grouping (``1.234.567``).
  ``detail/classic_stream.hpp`` (``make_classic_*``) pins every stream to the
  classic locale -- see that header's doc comment.

Every call site under ``src/cpp/`` was migrated; these are locale-free static
guards (no comma-decimal locale needs to be installed for them to run, unlike the
headers' own gtest suites) so a new call site added later cannot silently
reintroduce either bug.

Integer *parsing* (``strtol``/``strtoll``) is deliberately not guarded: the C
locale's ``LC_NUMERIC`` affects only the radix character, which an integer
literal doesn't have -- see ``fast_number.hpp``'s own doc comment. Integer
*stream* I/O is guarded, because ``std::num_put`` applies digit grouping, which
``snprintf``'s ``%d`` does not; that difference is why the two rules look
asymmetric.
"""

import pathlib
import re

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
SRC = REPO / "src" / "cpp"

# The header itself is the one place these names are allowed to appear (its
# own implementation, and doc-comment prose explaining the tiered fallback).
_ALLOWED_FILE = SRC / "include" / "meshioplusplus" / "detail" / "fast_number.hpp"

_FLOAT_FORMAT_SPEC = re.compile(r"%[-0-9.#+ *]*[eEfgG]")
_BANNED_CALL = re.compile(r"\b(strtod|strtof|atof|stod|stof)\s*\(")
_BANNED_PRINTF = re.compile(r"\b(snprintf|sprintf)\s*\(")

# A stream *declaration* of a concrete standard type: `std::ifstream in(...)`,
# `std::ostringstream os;`, `std::ofstream f = ...`. A reference/pointer
# parameter (`std::ifstream& in`) is not a construction, hence the lookahead.
# The fix is always `auto x = detail::make_classic_<kind>(...)`, which no longer
# matches -- a concrete stream type never appears as a declaration.
_STREAM_DECL = re.compile(
    r"\bstd::(?:i|o)?(?:string|f)stream\s+(?![&*])\**[A-Za-z_]\w*\s*[({=;]"
)

_CLASSIC_STREAM = SRC / "include" / "meshioplusplus" / "detail" / "classic_stream.hpp"
# `format_compat.hpp` is the one deliberate exception: it is an installed header
# whose function is a template, so imbuing inside it would be an ABI break (see
# doc/abi.md, Tier B), and it only renders log lines and exception text -- never
# a byte of a file, and never anything that is re-parsed.
_FORMAT_COMPAT = SRC / "include" / "meshioplusplus" / "detail" / "format_compat.hpp"
_STREAM_ALLOWED = (_CLASSIC_STREAM, _FORMAT_COMPAT)


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


def test_no_bare_stream_construction_outside_classic_stream_hpp():
    violations = []
    for path in SRC.rglob("*.[hc]pp"):
        if "third_party" in path.parts or path in _STREAM_ALLOWED:
            continue
        for lineno, line in enumerate(path.read_text().splitlines(), start=1):
            stripped = line.strip()
            if stripped.startswith(("//", "*", "/*")):
                continue
            if _STREAM_DECL.search(line):
                violations.append(f"{path.relative_to(REPO)}:{lineno}: {stripped}")

    assert not violations, (
        "a stream constructed without pinning the classic locale; use "
        "`auto x = detail::make_classic_<ifstream|ofstream|istringstream|"
        "ostringstream|stringstream>(...)` from detail/classic_stream.hpp "
        "(see that header for why):\n" + "\n".join(violations)
    )


def test_the_stream_guard_actually_sees_the_factories():
    # A guard that matches nothing passes vacuously; pin that the migration is
    # visible to it, and that the regex would have caught a reintroduction.
    factories = 0
    for path in SRC.rglob("*.[hc]pp"):
        if "third_party" in path.parts or path == _CLASSIC_STREAM:
            continue
        factories += len(re.findall(r"\bmake_classic_\w+stream\(", path.read_text()))
    assert factories >= 150, factories

    for reintroduced in (
        "    std::ifstream in(rPath);",
        "    std::ofstream os(rPath, std::ios::binary);",
        "    std::istringstream iss(line);",
        "    std::ostringstream os;",
        "    std::stringstream ss(line);",
        "    std::ofstream f = foam_open(rPath);",
    ):
        assert _STREAM_DECL.search(reintroduced), reintroduced
    for fine in (
        "    auto in = detail::make_classic_ifstream(rPath);",
        "void read_ascii(std::ifstream& rIn) {",
        "    std::ifstream* p = nullptr;",
    ):
        assert not _STREAM_DECL.search(fine), fine
