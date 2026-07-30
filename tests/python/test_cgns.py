import pathlib
import shutil
import subprocess

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._exceptions import ReadError, WriteError
from meshioplusplus.cgns import _cgns

from . import helpers

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)

_ROUNDTRIP_MESHES = [
    helpers.line_mesh,
    helpers.tri_mesh,
    helpers.triangle6_mesh,
    helpers.quad_mesh,
    helpers.quad8_mesh,
    helpers.tet_mesh,
    helpers.tet10_mesh,
    helpers.hex_mesh,
    helpers.hex20_mesh,
    helpers.wedge_mesh,
    helpers.wedge15_mesh,
    helpers.pyramid_mesh,
    helpers.tri_quad_mesh,
]

# FlowSolution_t coverage (v9.9.0). `helpers.write_read` compares point_data and
# cell_data, so these verify the field round-trip for free. Only
# single-dimension meshes: a mixed-dimension mesh's cell_data cannot be written
# (see test_mixed_dimension_cell_data_is_skipped).
_DATA_MESHES = [
    helpers.add_point_data(helpers.tri_mesh, 1),
    helpers.add_point_data(helpers.tri_mesh, 3),  # a genuine vector field
    helpers.add_point_data(helpers.tet_mesh, 3),
    helpers.add_cell_data(helpers.tri_mesh, [("a", (), float)]),
    helpers.add_cell_data(helpers.tri_quad_mesh, [("b", (3,), float)]),
    helpers.add_cell_data(helpers.tet_mesh, [("c", (6,), float)]),
]


@pytest.mark.parametrize("mesh", _ROUNDTRIP_MESHES)
def test(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.cgns.write, meshioplusplus.cgns.read, mesh, 1.0e-15
    )


@pytest.mark.parametrize("mesh", _ROUNDTRIP_MESHES)
def test_python_fallback(mesh, tmp_path):
    """The Python/h5py twin, bypassing the C++-preferring shim entirely --
    exercised even on a machine where the C++ core is present."""
    helpers.write_read(tmp_path, _cgns.write, _cgns.read, mesh, 1.0e-15)


@pytest.mark.parametrize("mesh", _DATA_MESHES)
def test_flow_solution(mesh, tmp_path):
    """point_data/cell_data survive a FlowSolution_t round-trip (v9.9.0);
    before that a CGNS export silently dropped every field."""
    helpers.write_read(
        tmp_path, meshioplusplus.cgns.write, meshioplusplus.cgns.read, mesh, 1.0e-15
    )


@pytest.mark.parametrize("mesh", _DATA_MESHES)
def test_flow_solution_python_fallback(mesh, tmp_path):
    helpers.write_read(tmp_path, _cgns.write, _cgns.read, mesh, 1.0e-15)


def test_vector_field_is_split_into_per_component_arrays(tmp_path):
    """CGNS has no NumberOfComponents, so a k>1 array becomes k sibling
    DataArray_t nodes named `<name>_<i>`, re-joined on read."""
    import h5py

    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2]]))],
    )
    mesh.point_data["velocity"] = np.arange(9, dtype=float).reshape(3, 3)
    mesh.point_data["temperature"] = np.arange(3, dtype=float)

    p = tmp_path / "vec.cgns"
    _cgns.write(p, mesh, compression=None)
    with h5py.File(p, "r") as f:
        sol = f["Base"]["Zone1"]["FlowSolution"]
        assert sorted(k for k in sol if k != "GridLocation") == [
            "temperature",
            "velocity_0",
            "velocity_1",
            "velocity_2",
        ]
        loc = bytes(np.asarray(sol["GridLocation"][" data"][()]).astype("i1"))
        assert loc == b"Vertex"
        # Each component node is a real DataArray_t, one value per vertex.
        assert sol["velocity_1"][" data"].shape == (3,)
        np.testing.assert_allclose(sol["velocity_1"][" data"][()], [1.0, 4.0, 7.0])

    back = _cgns.read(p)
    assert back.point_data["velocity"].shape == (3, 3)
    np.testing.assert_allclose(back.point_data["velocity"], mesh.point_data["velocity"])
    np.testing.assert_allclose(
        back.point_data["temperature"], mesh.point_data["temperature"]
    )


def test_cell_data_uses_cellcenter_grid_location(tmp_path):
    import h5py

    mesh = helpers.add_cell_data(helpers.tri_mesh, [("a", (), float)])
    p = tmp_path / "cells.cgns"
    _cgns.write(p, mesh, compression=None)
    with h5py.File(p, "r") as f:
        sol = f["Base"]["Zone1"]["FlowSolutionCells"]
        loc = bytes(np.asarray(sol["GridLocation"][" data"][()]).astype("i1"))
        assert loc == b"CellCenter"
        assert _read_label(sol) == "FlowSolution_t"


def _read_label(obj):
    v = obj.attrs["label"]
    return bytes(v).split(b"\0", 1)[0].decode()


def test_mixed_dimension_cell_data_is_skipped(tmp_path):
    """A zone-wide CellCenter array cannot be distributed back across blocks of
    different topological dimension, so it is warn-and-skipped rather than
    written wrongly. Geometry is unaffected."""
    import h5py

    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 1.0]]),
        [
            ("triangle", np.array([[0, 1, 2]])),
            ("tetra", np.array([[0, 1, 2, 3]])),
        ],
    )
    mesh.cell_data["a"] = [np.array([1.0]), np.array([2.0])]

    p = tmp_path / "mixed.cgns"
    with pytest.warns(UserWarning, match="topological dimensions"):
        _cgns.write(p, mesh, compression=None)
    with h5py.File(p, "r") as f:
        assert "FlowSolutionCells" not in f["Base"]["Zone1"]
    back = _cgns.read(p)
    assert not back.cell_data
    assert len(back.cells) == 2


def test_two_dimensional_points(tmp_path):
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )
    p = tmp_path / "test.cgns"
    _cgns.write(p, mesh, compression=None)
    out = _cgns.read(p)
    assert out.points.shape == (4, 2)
    assert np.allclose(out.points, mesh.points)

    if _HAS_HDF5:
        p2 = tmp_path / "test_cpp.cgns"
        _core.cgns_write(str(p2), mesh, -1)
        out2 = _core.cgns_read(str(p2))
        assert out2.points.shape == (4, 2)
        assert np.allclose(out2.points, mesh.points)


def _hex27_mesh():
    # Unit cube corners, edge midpoints, face centers, body center -- the
    # real vtkTriQuadraticHexahedron/meshio layout (mirrors test_skin.py's
    # _hex27_mesh, kept independent on purpose -- see that file).
    corners = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 1],
        ],
        dtype=float,
    )
    edges = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ]
    faces = [
        (0, 4, 7, 3),
        (1, 2, 6, 5),
        (0, 1, 5, 4),
        (3, 7, 6, 2),
        (0, 3, 2, 1),
        (4, 5, 6, 7),
    ]
    pts = list(corners)
    for a, b in edges:
        pts.append((corners[a] + corners[b]) / 2)
    for f in faces:
        pts.append(corners[list(f)].mean(axis=0))
    pts.append(corners.mean(axis=0))
    return meshioplusplus.Mesh(
        np.array(pts), [("hexahedron27", np.arange(27).reshape(1, 27))]
    )


def test_hexahedron27_round_trip(tmp_path):
    mesh = _hex27_mesh()
    p = tmp_path / "test.cgns"
    _cgns.write(p, mesh, compression=None)
    out = _cgns.read(p)
    assert np.allclose(out.points, mesh.points)
    assert out.cells[0].type == "hexahedron27"
    assert np.array_equal(out.cells[0].data, mesh.cells[0].data)


@pytest.mark.skipif(not _HAS_HDF5, reason="needs the C++ HDF5-enabled core")
def test_structural_parity_with_cpp(tmp_path):
    """The C++ and Python writers must produce structurally identical files
    (same groups/attributes/dataset shapes+dtypes+bytes), excluding
    ' hdf5version' -- that dataset records the *linked* HDF5 library version,
    which legitimately differs between the C++ core and h5py."""
    import h5py

    mesh = helpers.tri_quad_mesh
    p_cpp = tmp_path / "cpp.cgns"
    p_py = tmp_path / "py.cgns"
    _core.cgns_write(str(p_cpp), mesh, -1)
    _cgns.write(p_py, mesh, compression=None)

    def snapshot(g, path=""):
        out = {}
        for name in g:
            full = f"{path}/{name}"
            obj = g[name]
            if isinstance(obj, h5py.Group):
                attrs = {
                    k: (
                        bytes(v)
                        if isinstance(v, (bytes, np.bytes_))
                        else np.asarray(v).tolist()
                    )
                    for k, v in obj.attrs.items()
                }
                out[full] = ("group", attrs)
                out.update(snapshot(obj, full))
            else:
                if name == " hdf5version":
                    continue
                out[full] = ("dataset", obj.shape, str(obj.dtype), obj[()].tolist())
        return out

    with h5py.File(p_cpp, "r") as f_cpp, h5py.File(p_py, "r") as f_py:
        snap_cpp = snapshot(f_cpp)
        snap_py = snapshot(f_py)
        # Root-level attrs.
        for key in ("name", "label", "type"):
            assert bytes(f_cpp.attrs[key]) == bytes(f_py.attrs[key]), key

    assert snap_cpp.keys() == snap_py.keys()
    for key in snap_cpp:
        assert snap_cpp[key] == snap_py[key], key


def test_legacy_layout_still_reads(tmp_path):
    """The pre-v9.8.0 layout (no node attributes at all) must still read,
    through both the Python and C++ paths."""
    import h5py

    p = tmp_path / "legacy.cgns"
    with h5py.File(p, "w") as f:
        base = f.create_group("Base")
        zone1 = base.create_group("Zone1")
        coords = zone1.create_group("GridCoordinates")
        xyz = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        for i, name in enumerate(["CoordinateX", "CoordinateY", "CoordinateZ"]):
            g = coords.create_group(name)
            g.create_dataset(" data", data=xyz[:, i])
        elems = zone1.create_group("GridElements")
        rng = elems.create_group("ElementRange")
        rng.create_dataset(" data", data=[1, 1])
        conn = elems.create_group("ElementConnectivity")
        conn.create_dataset(" data", data=[1, 2, 3, 4])

    out = _cgns.read(p)
    assert out.points.shape == (4, 3)
    assert out.cells[0].type == "tetra"

    if _HAS_HDF5:
        out2 = _core.cgns_read(str(p))
        assert out2.points.shape == (4, 3)
        assert out2.cells[0].type == "tetra"


@pytest.mark.parametrize("write_fn", [_cgns.write, pytest.param(None, id="cpp")])
def test_ragged_block_raises(tmp_path, write_fn):
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]),
        [("polygon", [np.array([0, 1, 2])])],
    )
    p = tmp_path / "test.cgns"
    if write_fn is None:
        if not _HAS_HDF5:
            pytest.skip("needs the C++ HDF5-enabled core")
        with pytest.raises(WriteError, match="ragged block"):
            _core.cgns_write(str(p), mesh, -1)
    else:
        with pytest.raises(WriteError, match="ragged block"):
            write_fn(p, mesh, compression=None)


@pytest.mark.parametrize("write_fn", [_cgns.write, pytest.param(None, id="cpp")])
def test_unverified_ordering_type_raises(tmp_path, write_fn):
    mesh = meshioplusplus.Mesh(
        np.zeros((20, 3)), [("tetra20", np.arange(20).reshape(1, 20))]
    )
    p = tmp_path / "test.cgns"
    if write_fn is None:
        if not _HAS_HDF5:
            pytest.skip("needs the C++ HDF5-enabled core")
        with pytest.raises(WriteError, match="not yet verified"):
            _core.cgns_write(str(p), mesh, -1)
    else:
        with pytest.raises(WriteError, match="not yet verified"):
            write_fn(p, mesh, compression=None)


def test_structured_zone_raises(tmp_path):
    import h5py

    p = tmp_path / "test.cgns"
    _cgns.write(p, helpers.tet_mesh, compression=None)
    with h5py.File(p, "a") as f:
        del f["Base"]["Zone1"]["ZoneType"][" data"]
        f["Base"]["Zone1"]["ZoneType"].create_dataset(
            " data", data=np.frombuffer(b"Structured", dtype="i1")
        )
    with pytest.raises(ReadError, match="Structured"):
        _cgns.read(p)


def test_mixed_section_raises(tmp_path):
    import h5py

    p = tmp_path / "test.cgns"
    _cgns.write(p, helpers.tet_mesh, compression=None)
    with h5py.File(p, "a") as f:
        section = f["Base"]["Zone1"]["TETRA_4_1"]
        del section[" data"]
        section.create_dataset(" data", data=np.array([20, 0], dtype="<i4"))
    with pytest.raises(ReadError, match="MIXED"):
        _cgns.read(p)


# ---- external validation (see tests/python/meshes/cgns/README.md) ----------
#
# Everything above proves only that meshio++ agrees with itself. These two go
# outside: one reads a file real cgnslib wrote, the other runs cgnslib's own
# conformance checker over our output. The first is unconditional (the fixture
# is committed); the second is opt-in, because cgnscheck is not a pip/apt
# dependency of this repo -- it is skipped with an explicit reason rather than
# silently passing.


def test_read_cgnslib_reference_file():
    """A file whose bytes were written end to end by cgnslib 4.5.2, not by us."""
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "cgns" / "tri_quad_fields_cgnslib.cgns"
    mesh = meshioplusplus.cgns.read(filename)

    assert mesh.points.shape == (6, 3)
    assert [(cb.type, len(cb.data)) for cb in mesh.cells] == [
        ("triangle", 2),
        ("quad", 1),
    ]
    # The whole v9.9.0 FlowSolution surface, read back off cgnslib's own bytes.
    assert mesh.point_data["Temperature"].shape == (6,)
    np.testing.assert_allclose(
        mesh.point_data["Temperature"], np.arange(6, dtype=float) * 1.5
    )
    assert mesh.point_data["Velocity"].shape == (6, 3)
    np.testing.assert_allclose(
        mesh.point_data["Velocity"], np.arange(18, dtype=float).reshape(6, 3)
    )
    np.testing.assert_allclose(mesh.cell_data["Density"][0], [1.0, 2.0])
    np.testing.assert_allclose(mesh.cell_data["Density"][1], [3.0])


def test_python_reader_matches_cpp_on_the_reference_file():
    """Both readers must agree on an externally-written file, not just on ours."""
    if not _HAS_HDF5:
        pytest.skip("needs the C++ HDF5-enabled core")
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "cgns" / "tri_quad_fields_cgnslib.cgns"
    a = _cgns.read(filename)
    b = _core.cgns_read(str(filename))
    np.testing.assert_allclose(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for name in a.point_data:
        np.testing.assert_allclose(a.point_data[name], b.point_data[name])


@pytest.mark.skipif(
    shutil.which("cgnscheck") is None,
    reason="cgnscheck (cgnslib) not on PATH; install e.g. "
    "`micromamba create -r .micromamba -n cgns -c conda-forge cgns` and add its "
    "bin/ to PATH to run this",
)
@pytest.mark.parametrize(
    "mesh",
    _ROUNDTRIP_MESHES + [helpers.add_point_data(helpers.tri_mesh, 3)],
)
def test_cgnscheck_accepts_our_output(mesh, tmp_path):
    """cgnslib's own conformance checker must report no ERROR on what we write.

    Warnings are tolerated deliberately: cgnscheck recommends a `Family_t` on
    every zone and a `DataClass_t` on every array, and flags any field name that
    is not one of SIDS's standard identifiers -- meshio++ preserves the caller's
    own field names instead. None of those is a conformance failure.
    """
    p = tmp_path / "out.cgns"
    meshioplusplus.cgns.write(p, mesh)
    r = subprocess.run(
        ["cgnscheck", "-w3", str(p)], capture_output=True, text=True, check=False
    )
    errors = [ln for ln in (r.stdout + r.stderr).splitlines() if "ERROR" in ln.upper()]
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)
    assert not errors, errors
