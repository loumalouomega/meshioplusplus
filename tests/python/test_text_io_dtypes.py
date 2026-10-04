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
    "libmesh",
    "gid",
    "gltf",
    "obj",
    "stl",
    "triangle",
    "tetgen",
    "ugrid",
    "dolfin-xml",
    "freefem",
    "avsucd",
    "wkt",
    "svg",
    "tikz",
    "su2",
    "netgen",
    "elmer",
    "mphtxt",
    "code_aster",
    "nastran",
    "tecplot",
    "ansys",
    "radioss",
    "unv",
)
# Formats that record each array's dtype in the file (PCD SIZE/TYPE, HDF5
# datasets), so their bytes legitimately differ from canonical storage: the
# values must still round-trip.
STORED_DTYPE_FORMATS = ("pcd", "cgns", "med", "ply")
# A writer may embed its own file name (glTF's .bin), pick its extension (GiD)
# or write companion files (Triangle and TetGen's .ele), so every output is
# `<dir>/m<suffix>` in a directory of its own and the whole directory is
# compared.
SUFFIX = {
    "gid": ".post.msh",
    "gltf": ".gltf",
    "obj": ".obj",
    "stl": ".stl",
    "ply": ".ply",
    "triangle": ".node",
    "tetgen": ".node",
    "ugrid": ".lb8.ugrid",
    "dolfin-xml": ".xml",
    "freefem": ".msh",
    "avsucd": ".inp",
    "wkt": ".wkt",
    "svg": ".svg",
    "tikz": ".tex",
    "su2": ".su2",
    "netgen": ".vol",
    "elmer": "",
    "mphtxt": ".mphtxt",
    "code_aster": ".mail",
    "nastran": ".bdf",
    "tecplot": ".dat",
    "ansys": ".msh",
    "radioss": ".rad",
    "unv": ".unv",
    "pcd": ".pcd",
    "cgns": ".cgns",
    "med": ".med",
}
# These writers pick their output by whether an array is float or integer
# (UGRID labels, AVS-UCD materials, DOLFIN's `float`/`int` mesh functions, SU2
# and Netgen region tags), so the canonical form is the same class in float64
# or int64.
CLASS_FORMATS = ("ugrid", "avsucd", "dolfin-xml", "su2", "netgen")
# Writers that take only some meshes: Triangle writes 2-D points only.
ONLY_DIMENSION = {"triangle": 2}
# The integer cell data each format reads per cell, as (name, value).
TAGS = {
    "gmsh22": (("gmsh:physical", 3), ("gmsh:geometrical", 5)),
    "gmsh": (("gmsh:physical", 3), ("gmsh:geometrical", 5)),
    "mdpa": (("gmsh:physical", 3),),
    "patran": (("patran:property", 4),),
    "femap": (("femap:property", 4),),
    "unv": (("unv:pid", 3), ("unv:mid", 5)),
    "nastran": (("nastran:ref", 4),),
    "radioss": (("radioss:part", 2),),
    "mphtxt": (("mphtxt:geom", 2),),
    "ansys": (("ansys:zone", 2),),
    "su2": (("gmsh:physical", 3),),
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
    if fmt == "tetgen":
        points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=real)
        conn = np.array([[0, 1, 2, 3]], dtype=index)
        scalar = np.array([1, 2, 3, 4], dtype=real)
        ref = np.array([7], dtype=index)
        mesh = pp.Mesh(
            points,
            [("tetra", conn)],
            point_data={"s": scalar},
            cell_data={"pf3:ref": [ref]},
        )
        return mesh, [points, conn, scalar, ref]
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
    if ONLY_DIMENSION.get(fmt, dimension) != dimension:
        pytest.skip(f"{fmt} writes {ONLY_DIMENSION[fmt]}-D points only")
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    mesh, inputs = _mesh(fmt, dtype, dtype, dimension)
    canonical_dtype = ("float64", "int64")
    if fmt in CLASS_FORMATS:
        kind = "float64" if np.issubdtype(np.dtype(dtype), np.floating) else "int64"
        canonical_dtype = (kind, kind)
    expected, _ = _mesh(fmt, *canonical_dtype, dimension)
    before = [a.tobytes() for a in inputs]
    name = "m" + SUFFIX.get(fmt, "")
    native, canonical = tmp_path / "native", tmp_path / "canonical"
    native.mkdir()
    canonical.mkdir()
    pp.write(native / name, mesh, file_format=fmt)
    pp.write(canonical / name, expected, file_format=fmt)

    def outputs_of(root):  # Elmer writes a directory, so walk recursively
        return {
            str(f.relative_to(root)): f.read_bytes()
            for f in root.rglob("*")
            if f.is_file()
        }

    outputs = outputs_of(native)
    assert outputs  # something was written
    assert outputs == outputs_of(canonical)
    assert before == [a.tobytes() for a in inputs]


@pytest.mark.parametrize("fmt", STORED_DTYPE_FORMATS)
@pytest.mark.parametrize("dtype", DTYPES)
def test_stored_dtype_writer_round_trips_values_from_all_dtypes(
    tmp_path, monkeypatch, fmt, dtype
):
    monkeypatch.setenv("MESHIOPLUSPLUS_STRICT_CORE", "1")
    mesh, inputs = _mesh(fmt, dtype, dtype, 3)
    expected, _ = _mesh(fmt, "float64", "int64", 3)
    before = [a.tobytes() for a in inputs]
    suffix = SUFFIX[fmt]
    back = {}
    for label, source in (("native", mesh), ("canonical", expected)):
        path = tmp_path / (label + suffix)
        pp.write(path, source, file_format=fmt)
        back[label] = pp.read(path, file_format=fmt)
    assert before == [a.tobytes() for a in inputs]
    assert np.array_equal(back["native"].points, back["canonical"].points)
    assert len(back["native"].cells) == len(back["canonical"].cells)
    for a, b in zip(back["native"].cells, back["canonical"].cells):
        assert a.type == b.type and np.array_equal(a.data, b.data)
    for name, values in back["canonical"].point_data.items():
        assert np.array_equal(back["native"].point_data[name], values)


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
