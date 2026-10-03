"""Each hoisted writer preserves bytes/dtype semantics and caller buffers."""

import numpy as np
import pytest

import meshioplusplus as pp

FORMATS = (
    "off",
    "ip",
    "flux",
    "permas",
    "gmsh22",
    "gmsh",
    "febio",
    "z88",
    "patran",
    "femap",
    "mdpa",
)
# The integer cell data each format reads per cell, as (name, value).
TAGS = {
    "gmsh22": (("gmsh:physical", 3), ("gmsh:geometrical", 5)),
    "gmsh": (("gmsh:physical", 3), ("gmsh:geometrical", 5)),
    "mdpa": (("gmsh:physical", 3),),
    "patran": (("patran:property", 4),),
    "femap": (("femap:property", 4),),
}
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


def _mesh(fmt, real, index, dimension):
    """The mesh and its arrays: points/point data in @p real, connectivity and
    cell data in @p index. Z88 has no linear triangle, so it gets a tetrahedron.
    """
    if fmt == "z88":
        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=real)
        conn = np.array([[0, 1, 2, 3]], dtype=index)
        types = np.array([17], dtype=index)
        mesh = pp.Mesh(points, [("tetra", conn)], cell_data={"z88:type": [types]})
        return mesh, [points, conn, types]
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])[:, :dimension].astype(real)
    conn = np.array([[0, 1, 2]], dtype=index)
    scalar = np.array([1, 2, 3], dtype=real)
    vector = np.array([[1, 2], [3, 4], [5, 6]], dtype=real)
    cell_data = {"pf3:ref": [np.array([7], dtype=index)]}
    for name, value in TAGS.get(fmt, ()):
        cell_data[name] = [np.array([value], dtype=index)]
    mesh = pp.Mesh(
        points,
        [("triangle", conn)],
        point_data={"s": scalar, "v": vector},
        cell_data=cell_data,
    )
    return mesh, [points, conn, scalar, vector] + [v[0] for v in cell_data.values()]


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("dtype", DTYPES)
@pytest.mark.parametrize("dimension", [2, 3])
def test_hoisted_writer_all_dtypes_match_default_storage(
    tmp_path, monkeypatch, fmt, dtype, dimension
):
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    mesh, inputs = _mesh(fmt, dtype, dtype, dimension)
    expected, _ = _mesh(fmt, "float64", "int64", dimension)
    before = [a.tobytes() for a in inputs]
    native, canonical = tmp_path / "native", tmp_path / "canonical"
    pp.write(native, mesh, file_format=fmt)
    pp.write(canonical, expected, file_format=fmt)
    assert native.read_bytes() == canonical.read_bytes()
    assert before == [a.tobytes() for a in inputs]


@pytest.mark.parametrize("fmt", ("off", "ip", "flux", "permas"))
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
