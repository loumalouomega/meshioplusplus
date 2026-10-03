"""Each hoisted writer preserves bytes/dtype semantics and caller buffers."""

import numpy as np
import pytest

import meshioplusplus as pp

FORMATS = ("off", "ip", "flux", "permas")
DTYPES = (
    "float32",
    "float64",
    "int8",
    "int16",
    "int32",
    "int64",
    "uint8",
    "uint16",
    "uint32",
    "uint64",
)


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("dtype", DTYPES)
@pytest.mark.parametrize("dimension", [2, 3])
def test_hoisted_writer_all_dtypes_match_default_storage(
    tmp_path, monkeypatch, fmt, dtype, dimension
):
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])[:, :dimension].astype(dtype)
    conn = np.array([[0, 1, 2]], dtype=dtype)
    scalar = np.array([1, 2, 3], dtype=dtype)
    vector = np.array([[1, 2], [3, 4], [5, 6]], dtype=dtype)
    refs = np.array([7], dtype=dtype)
    mesh = pp.Mesh(
        points,
        [("triangle", conn)],
        point_data={"s": scalar, "v": vector},
        cell_data={"pf3:ref": [refs]},
    )
    expected = pp.Mesh(
        points.astype("float64"),
        [("triangle", conn.astype("int64"))],
        point_data={"s": scalar.astype("float64"), "v": vector.astype("float64")},
        cell_data={"pf3:ref": [refs.astype("int64")]},
    )
    inputs = [points, conn, scalar, vector, refs]
    before = [a.tobytes() for a in inputs]
    native, canonical = tmp_path / "native", tmp_path / "canonical"
    pp.write(native, mesh, file_format=fmt)
    pp.write(canonical, expected, file_format=fmt)
    assert native.read_bytes() == canonical.read_bytes()
    assert before == [a.tobytes() for a in inputs]


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("dtype", DTYPES)
def test_hoisted_writer_empty_arrays_keep_writer_contract(
    tmp_path, monkeypatch, fmt, dtype
):
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    points = np.empty((0, 3), dtype=dtype)
    conn = np.empty((0, 3), dtype=dtype)
    mesh = pp.Mesh(points, [("triangle", conn)])
    canonical = pp.Mesh(points.astype("float64"), [("triangle", conn.astype("int64"))])
    first, second = tmp_path / "first", tmp_path / "second"
    pp.write(first, mesh, file_format=fmt)
    pp.write(second, canonical, file_format=fmt)
    assert first.read_bytes() == second.read_bytes()
