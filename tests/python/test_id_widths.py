"""Integers above 2**31 either round-trip or fail by name (roadmap 3.19)."""

import importlib.util
import pathlib

import pytest

from . import id_widths_spec as spec

REPO = pathlib.Path(__file__).resolve().parents[2]


def _declared():
    return sorted(spec.EXPECTED)


def test_every_conformance_format_is_declared():
    needed = {f for f, s in spec.cs.SPEC.items() if s and "cells" in s}
    assert (
        needed - set(spec.EXPECTED) == set()
    ), "declare these in id_widths_spec.EXPECTED"


@pytest.mark.parametrize("fmt", _declared())
def test_observed_matches_declared(fmt):
    seen = spec.observe(fmt)
    want = spec.EXPECTED[fmt]
    assert seen is not None
    if seen.get("error") == "not observed":
        pytest.skip("an optional library is not installed")
    assert "WRAPPED" not in repr(seen), f"{fmt} silently narrows an integer: {seen}"
    if "error" in seen:
        assert seen.get("named"), f"{fmt}: the refusal does not name the format: {seen}"
    assert {k: v for k, v in seen.items() if k != "named"} == want


def test_the_page_is_fresh():
    gen_path = REPO / "tools" / "gen_id_widths_table.py"
    s = importlib.util.spec_from_file_location("gen_id_widths_table", gen_path)
    gen = importlib.util.module_from_spec(s)
    s.loader.exec_module(gen)
    assert gen.PAGE.read_text(encoding="utf-8") == gen.render(spec.EXPECTED)


def _triangle(**point_data):
    import numpy as np

    import meshioplusplus as mio

    return mio.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]]),
        [("triangle", np.array([[0, 1, 2]]))],
        point_data=point_data,
    )


def test_gifti_refuses_an_integer_it_would_wrap(tmp_path):
    import numpy as np

    import meshioplusplus as mio

    mesh = _triangle(tag=np.array([1, 2, 2**31 + 5], dtype=np.int64))
    with pytest.raises(mio.WriteError, match="gifti.*32-bit"):
        mio.write(tmp_path / "m.gii", mesh, file_format="gifti")
    ok = _triangle(tag=np.array([1, 2, 2**31 - 1], dtype=np.int64))
    mio.write(tmp_path / "ok.gii", ok, file_format="gifti")
    back = mio.read(tmp_path / "ok.gii", file_format="gifti")
    assert back.point_data["tag"].tolist() == [1, 2, 2**31 - 1]


def test_netgen_refuses_a_name_table_it_cannot_write(tmp_path):
    import numpy as np

    import meshioplusplus as mio

    mesh = mio.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        [("tetra", np.array([[0, 1, 2, 3]]))],
        cell_data={"netgen:index": [np.array([2**31 + 5], dtype=np.int64)]},
    )
    with pytest.raises(mio.WriteError, match="Netgen"):
        mio.write(tmp_path / "m.vol", mesh, file_format="netgen")
