"""The built-in colormap LUTs and their C++/Python parity.

The tables exist twice -- `src/cpp/include/meshioplusplus/detail/colormap.hpp` and
`src/python/meshioplusplus/_colormap.py`, both emitted by `tools/gen_colormaps.py` --
so a hand-edit to either copy has to fail here rather than silently desync the
SVG/TikZ writers' byte-parity guarantee.
"""

import pytest

import meshioplusplus
from meshioplusplus import _colormap

# Canonical matplotlib endpoints and midpoints, as a sanity anchor independent
# of both copies of the table.
KNOWN = {
    "viridis": ((68, 1, 84), (33, 145, 140), (253, 231, 37)),
    "coolwarm": ((59, 76, 192), (221, 220, 220), (180, 4, 38)),
    "turbo": ((48, 18, 59), (164, 252, 60), (122, 4, 3)),
}


@pytest.mark.parametrize("name", ["viridis", "coolwarm", "turbo"])
def test_table_matches_cpp(name):
    # The whole point of the generator: both copies are the same 768 bytes.
    assert _colormap.colormap_table(name) == meshioplusplus._core.colormap_table(name)


def test_names_match_cpp():
    assert _colormap.NAMES == meshioplusplus._core.colormap_names()


@pytest.mark.parametrize("name", ["viridis", "coolwarm", "turbo"])
def test_known_values(name):
    table = _colormap.colormap_table(name)
    assert len(table) == _colormap.COLORMAP_SIZE * 3
    lo, mid, hi = KNOWN[name]
    assert _colormap.colormap_lookup(table, 0.0) == lo
    assert _colormap.colormap_lookup(table, 0.5) == mid
    assert _colormap.colormap_lookup(table, 1.0) == hi


@pytest.mark.parametrize("name", ["viridis", "coolwarm", "turbo"])
def test_lookup_clamps_and_handles_nan(name):
    table = _colormap.colormap_table(name)
    lo = _colormap.colormap_lookup(table, 0.0)
    hi = _colormap.colormap_lookup(table, 1.0)
    assert _colormap.colormap_lookup(table, -5.0) == lo
    assert _colormap.colormap_lookup(table, 5.0) == hi
    # NaN fails every comparison, so the `not (t > 0.0)` guard catches it. The
    # writers route NaN to nan_color before this point; this pins the fallback.
    assert _colormap.colormap_lookup(table, float("nan")) == lo


@pytest.mark.parametrize("k", [1, 3, 5, 7, 9, 11, 255, 509])
def test_tie_prone_indices_agree_with_cpp(k):
    # t = k/510 for odd k makes t * 255.0 + 0.5 land on an exact .5 tie, which
    # is precisely where std::lround (half away from zero) and Python round()
    # (half to even) would disagree. Both implementations truncate instead, so
    # these must match -- this is the test that would catch a "cleanup" to
    # round() on either side.
    t = k / 510.0
    table = _colormap.colormap_table("viridis")
    expected_index = int(t * 255.0 + 0.5)
    assert _colormap.colormap_lookup(table, t) == (
        table[expected_index * 3],
        table[expected_index * 3 + 1],
        table[expected_index * 3 + 2],
    )


def test_unknown_name_raises():
    with pytest.raises(ValueError, match="unknown colormap"):
        _colormap.colormap_table("nope")
    with pytest.raises(Exception, match="unknown colormap"):
        meshioplusplus._core.colormap_table("nope")
