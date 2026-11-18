"""Formats that store integers as floats say so (roadmap 3.19)."""

import pathlib
import re

import numpy as np
import pytest

import meshioplusplus as mio
from meshioplusplus import _helpers

from . import conformance_spec as cs

REPO = pathlib.Path(__file__).resolve().parents[2]
FORMATS = sorted(_helpers._FLOAT_STORED_INTEGER_FORMATS)


def _mesh(fmt, value, dtype=np.int64):
    spec = cs.SPEC[fmt]
    ctypes = [c for c in cs.CELL_TYPES if c in spec["cells"]]
    mesh = cs.canonical_mesh(
        ctypes, data=False, regions=False, planar=spec.get("input") == "2d"
    )
    mesh.point_data["tag"] = np.full(len(mesh.points), value, dtype=dtype)
    mesh.cell_data["tag"] = [
        np.full(len(b.data), value, dtype=dtype) for b in mesh.cells
    ]
    return mesh


def _write(fmt, mesh, tmp_path):
    try:
        cs.round_trip(fmt, mesh, tmp_path)
    except ImportError:
        pytest.skip("an optional library is not installed")


@pytest.mark.parametrize("fmt", FORMATS)
def test_a_large_integer_array_warns(fmt, tmp_path, capsys):
    _write(fmt, _mesh(fmt, 2**24 + 1), tmp_path)
    err = capsys.readouterr().err
    assert "'tag' is stored as floating point" in err, err


@pytest.mark.parametrize("fmt", FORMATS)
def test_small_integers_and_floats_are_silent(fmt, tmp_path, capsys):
    for name, mesh in (("a", _mesh(fmt, 2**24)), ("b", _mesh(fmt, 2**40, np.float64))):
        (tmp_path / name).mkdir()
        _write(fmt, mesh, tmp_path / name)
    assert "stored as floating point" not in capsys.readouterr().err


def test_a_format_with_integer_storage_is_silent(tmp_path, capsys):
    mesh = _mesh("vtu", 2**40)
    mio.write(tmp_path / "m.vtu", mesh)
    assert "stored as floating point" not in capsys.readouterr().err


def test_the_native_table_matches_the_python_one():
    text = (REPO / "src/cpp/src/detail/integer_precision.cpp").read_text(
        encoding="utf-8"
    )
    body = text[text.index("ip_formats()") :]
    body = body[body.index("{") : body.index("};")]
    native = {
        name: (point == "true", cell == "true")
        for name, point, cell in re.findall(
            r'\{"(\w+)", \{(true|false), (true|false)\}\}', body
        )
    }
    python = {
        name: ("point" in kinds, "cell" in kinds)
        for name, kinds in _helpers._FLOAT_STORED_INTEGER_FORMATS.items()
    }
    assert native == python
