"""The ghost/owner/global-id array names are one convention in two languages.

``src/cpp/src/detail/ghost_names.hpp`` is the C++ twin of
``meshioplusplus/_ghosts.py``; see ``doc/ghosts.md``. The C++ header is private
(no installed header may name it), so the names are compared by parsing it.
"""

import pathlib
import re

from meshioplusplus import _ghosts

HEADER = (
    pathlib.Path(__file__).parents[2]
    / "src"
    / "cpp"
    / "src"
    / "detail"
    / "ghost_names.hpp"
)


def _cpp_constants():
    text = HEADER.read_text(encoding="utf-8")
    out = {}
    for name, value in re.findall(
        r"inline constexpr const char\*\s+(\w+)\s*=\s*\"([^\"]*)\"\s*;", text
    ):
        out[name] = value
    for name, value in re.findall(
        r"inline constexpr unsigned char\s+(\w+)\s*=\s*(\d+)\s*;", text
    ):
        out[name] = int(value)
    return out


def _snake(name):
    # kPartitionOriginalPointId -> PARTITION_ORIGINAL_POINT_ID
    return re.sub(r"(?<!^)(?=[A-Z])", "_", name[1:]).upper()


def test_constants_match_the_cpp_header():
    cpp = _cpp_constants()
    assert cpp, "no constants parsed from ghost_names.hpp"
    for name, value in cpp.items():
        py_name = _snake(name)
        assert hasattr(_ghosts, py_name), f"{py_name} missing from _ghosts.py"
        assert getattr(_ghosts, py_name) == value, py_name
    public = {n for n in dir(_ghosts) if n.isupper()}
    assert public == {_snake(n) for n in cpp}
