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


def _polygon_mesh():
    return meshioplusplus.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]),
        [("polygon", [np.array([0, 1, 2])])],
    )


def test_ragged_block_writes_as_ngon_in_the_cpp_core(tmp_path):
    # Until v9.21.0 both paths refused this outright; the C++ core now emits an
    # NGON_n section.
    if not _HAS_HDF5:
        pytest.skip("needs the C++ HDF5-enabled core")
    p = tmp_path / "test.cgns"
    _core.cgns_write(str(p), _polygon_mesh(), -1)
    out = _core.cgns_read(str(p))
    assert len(out.points) == 3
    assert any(c.type.startswith("polygon") for c in out.cells)


def test_ragged_block_names_the_compiled_core_in_the_python_twin(tmp_path):
    # The h5py twin deliberately does not reimplement NGON_n/NFACE_n (the face
    # dedup + winding repair would be a second implementation of a discrete
    # sign branch). It must say which path is needed, not "unsupported".
    p = tmp_path / "test.cgns"
    with pytest.raises(WriteError, match="compiled core"):
        _cgns.write(p, _polygon_mesh(), compression=None)


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


def test_unknown_zone_type_raises(tmp_path):
    import h5py

    p = tmp_path / "test.cgns"
    _cgns.write(p, helpers.tet_mesh, compression=None)
    with h5py.File(p, "a") as f:
        del f["Base"]["Zone1"]["ZoneType"][" data"]
        f["Base"]["Zone1"]["ZoneType"].create_dataset(
            " data", data=np.frombuffer(b"Weird", dtype="i1")
        )
    with pytest.raises(ReadError, match="Weird"):
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


# ---------------------------------------------------------------------------
# NGON_n / NFACE_n -- polyhedral cells (C++ core only; see _cgns.py)
# ---------------------------------------------------------------------------

_POLY_PTS = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
        [2, 0, 0],
        [2, 1, 0],
        [2, 1, 1],
        [2, 0, 1],
    ],
    dtype=float,
)
_CUBE_A = [
    [0, 3, 2, 1],
    [4, 5, 6, 7],
    [0, 1, 5, 4],
    [2, 3, 7, 6],
    [1, 2, 6, 5],
    [0, 4, 7, 3],
]
_CUBE_B = [
    [1, 2, 9, 8],
    [5, 11, 10, 6],
    [1, 8, 11, 5],
    [2, 6, 10, 9],
    [8, 9, 10, 11],
    [1, 5, 6, 2],
]


def _poly_meshes():
    return {
        "one": meshioplusplus.Mesh(_POLY_PTS[:8], [("polyhedron8", [_CUBE_A])]),
        "two": meshioplusplus.Mesh(_POLY_PTS, [("polyhedron8", [_CUBE_A, _CUBE_B])]),
        "mixed": meshioplusplus.Mesh(
            _POLY_PTS,
            [
                ("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]])),
                ("polyhedron8", [_CUBE_B]),
            ],
        ),
    }


@pytest.mark.skipif(not _HAS_HDF5, reason="needs the C++ HDF5-enabled core")
@pytest.mark.parametrize("key", ["one", "two", "mixed"])
def test_polyhedral_round_trip(key, tmp_path):
    m = _poly_meshes()[key]
    p = tmp_path / "poly.cgns"
    _core.cgns_write(str(p), m, -1)
    back = _core.cgns_read(str(p))

    assert len(back.points) == len(m.points)
    n_poly = sum(len(c.data) for c in back.cells if c.type.startswith("polyhedron"))
    expected = sum(len(c.data) for c in m.cells if c.type.startswith("polyhedron"))
    assert n_poly == expected
    if key == "mixed":
        # The hexahedron keeps its own HEXA_8 section rather than being
        # demoted, which is what the block-filtered face table is for.
        assert any(c.type == "hexahedron" for c in back.cells)


@pytest.mark.skipif(not _HAS_HDF5, reason="needs the C++ HDF5-enabled core")
@pytest.mark.skipif(
    shutil.which("cgnscheck") is None,
    reason="cgnscheck (cgnslib) not on PATH",
)
@pytest.mark.parametrize("key", ["one", "two", "mixed"])
def test_cgnscheck_accepts_our_polyhedral_output(key, tmp_path):
    """The oracle that actually found the bugs in this feature.

    cgnscheck reported no ERROR yet *aborted* on the first version of this
    writer -- twice: once because NCell omitted the polyhedral blocks, and once
    because the file declared CGNSLibraryVersion 3.1, under which cgnslib reads
    NGON_n with the 3.x inline-length layout. Both corrupted cgnslib's heap
    rather than producing a diagnostic, so the return code matters as much as
    the absence of ERROR lines.
    """
    p = tmp_path / "poly.cgns"
    _core.cgns_write(str(p), _poly_meshes()[key], -1)
    r = subprocess.run(
        ["cgnscheck", "-w3", str(p)], capture_output=True, text=True, check=False
    )
    errors = [ln for ln in (r.stdout + r.stderr).splitlines() if "ERROR" in ln.upper()]
    assert r.returncode == 0, (r.returncode, r.stdout, r.stderr)
    assert not errors, errors


def _regions_mesh():
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1.0]]),
        [
            ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
            ("quad", np.array([[0, 1, 2, 3]])),
            ("line", np.array([[0, 1]])),
        ],
    )
    mesh.field_data["answer"] = np.array([42.0])
    mesh.field_data["matrix"] = np.arange(6).reshape(2, 3)
    mesh.regions = [
        meshioplusplus.Region("inlet", "point", np.array([0, 3]), dim=0, tag=5),
        meshioplusplus.Region("skin", "cell", np.array([0, 2]), dim=2, tag=7),
        meshioplusplus.Region("a/b", "cell", np.array([3]), dim=1, tag=1),
        meshioplusplus.Region("x" * 40, "point", np.array([1])),
    ]
    return mesh


def _region_rows(mesh):
    return sorted(
        (r.name, r.kind, r.entries.tolist(), r.dim, r.tag) for r in mesh.regions
    )


@pytest.mark.parametrize("writer", ["core", "python"])
@pytest.mark.parametrize("reader", ["core", "python"])
def test_field_data_and_regions_round_trip_across_engines(tmp_path, writer, reader):
    if "core" in (writer, reader) and not _HAS_HDF5:
        pytest.skip("needs the C++ HDF5-enabled core")
    mesh = _regions_mesh()
    path = str(tmp_path / "r.cgns")
    if writer == "core":
        _core.cgns_write(path, mesh, 4)
    else:
        _cgns.write(path, mesh)
    back = _core.cgns_read(path) if reader == "core" else _cgns.read(path)
    assert _region_rows(back) == _region_rows(mesh)
    assert sorted(back.field_data) == ["answer", "matrix"]
    np.testing.assert_array_equal(back.field_data["matrix"], mesh.field_data["matrix"])
    np.testing.assert_array_equal(back.field_data["answer"], [42.0])


@pytest.mark.skipif(not _HAS_HDF5, reason="needs the C++ HDF5-enabled core")
def test_regions_are_zone_sub_regions_in_the_file(tmp_path):
    h5py = pytest.importorskip("h5py")
    path = str(tmp_path / "r.cgns")
    _core.cgns_write(path, _regions_mesh(), 4)
    with h5py.File(path, "r") as f:
        zone = f["Base/Zone1"]
        assert zone["inlet"].attrs["label"].decode().strip("\0") == "ZoneSubRegion_t"
        assert zone["inlet/PointList/ data"][()].reshape(-1).tolist() == [1, 4]
        # the cell region lists element ids: the triangles are 1 and 2, the quad 3
        assert zone["skin/PointList/ data"][()].reshape(-1).tolist() == [1, 3]
        assert (
            zone["FieldData"].attrs["label"].decode().strip("\0") == "UserDefinedData_t"
        )


@pytest.mark.skipif(not _HAS_HDF5, reason="needs the C++ HDF5-enabled core")
def test_side_regions_are_not_written(tmp_path):
    mesh = _regions_mesh()
    mesh.regions.append(
        meshioplusplus.Region("sides", "side", np.array([[0, 1]]), dim=1, tag=2)
    )
    for write in (
        lambda p: _core.cgns_write(p, mesh, 4),
        lambda p: _cgns.write(p, mesh),
    ):
        path = str(tmp_path / "s.cgns")
        write(path)
        assert "sides" not in {r.name for r in _core.cgns_read(path).regions}


@pytest.mark.parametrize("engine", ["python", "core"])
def test_every_zone_is_read_and_a_narrower_zone_is_padded(tmp_path, engine):
    """Roadmap 2.3.1.1: zones concatenate in name order; a zone with fewer
    coordinates than another is padded with zeros rather than misaligned."""
    import h5py

    pts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    mesh = meshioplusplus.Mesh(pts, [("triangle", np.array([[0, 1, 2]]))])
    p = tmp_path / "two_zone.cgns"
    _cgns.write(p, mesh)
    with h5py.File(p, "r+") as f:
        f.copy("Base/Zone1", "Base/Zone2")
        zone2 = f["Base/Zone2"]
        zone2.attrs["name"] = np.bytes_("Zone2".ljust(33, "\0"))
        del zone2["GridCoordinates/CoordinateZ"]
    read = _cgns.read if engine == "python" else meshioplusplus.read
    if engine == "core" and not hasattr(_core, "cgns_read"):
        pytest.skip("no compiled CGNS reader")
    out = read(p)
    assert out.points.shape == (6, 3)
    np.testing.assert_array_equal(out.points[:3], pts)
    np.testing.assert_array_equal(out.points[3:, :2], pts[:, :2])
    np.testing.assert_array_equal(out.points[3:, 2], 0.0)
    assert [c.type for c in out.cells] == ["triangle", "triangle"]
    np.testing.assert_array_equal(out.cells[1].data, [[3, 4, 5]])


def _add_bc(
    zone, name, bc_type, location, point_list=None, point_range=None, family=None
):
    """A ``BC_t`` under ``ZoneBC``, written with the writer's own node helpers."""
    if "ZoneBC" not in zone:
        bc_root = _cgns._create_group(zone, "ZoneBC")
        _cgns._write_node_attrs(bc_root, "ZoneBC", "ZoneBC_t", "MT")
    node = _cgns._c1_node(zone["ZoneBC"], name, "BC_t", bc_type)
    if point_list is not None:
        pl = _cgns._create_group(node, "PointList")
        _cgns._write_node_attrs(pl, "PointList", "IndexArray_t", "I4")
        pl.create_dataset(" data", data=np.array([point_list], dtype="<i4"))
    if point_range is not None:
        pr = _cgns._create_group(node, "PointRange")
        _cgns._write_node_attrs(pr, "PointRange", "IndexRange_t", "I4")
        pr.create_dataset(" data", data=np.array([point_range], dtype="<i4"))
    _cgns._c1_node(node, "GridLocation", "GridLocation_t", location)
    if family:
        _cgns._c1_node(node, "FamilyName", "FamilyName_t", family)


@pytest.mark.parametrize("engine", ["python", "core"])
def test_boundary_conditions_are_read_as_regions(tmp_path, engine):
    h5py = pytest.importorskip("h5py")
    if engine == "core" and not _HAS_HDF5:
        pytest.skip("needs the C++ HDF5-enabled core")
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0], [1, 0], [1, 1], [0, 1.0]]),
        [
            ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
            ("line", np.array([[0, 1], [1, 2], [2, 3], [3, 0]])),
        ],
    )
    path = str(tmp_path / "bc.cgns")
    _cgns.write(path, mesh)
    with h5py.File(path, "r+") as f:
        zone = f["Base/Zone1"]
        # the line section's elements are ids 3 to 6
        _add_bc(zone, "inlet", "BCInflow", "Vertex", point_list=[1, 4])
        _add_bc(
            zone, "wall_a", "BCWall", "EdgeCenter", point_range=[3, 4], family="walls"
        )
        _add_bc(zone, "wall_b", "BCWall", "EdgeCenter", point_list=[6], family="walls")
    read = _cgns.read if engine == "python" else meshioplusplus.cgns.read
    back = read(path)
    by_name = {(r.name, r.kind): sorted(r.entries.tolist()) for r in back.regions}
    assert by_name[("inlet", "point")] == [0, 3]
    # the two BCs of the family merge; the lines start at cell 2
    assert by_name[("walls", "cell")] == [2, 3, 5]
    assert bytes(back.field_data["cgns:bc_type:inlet"]).decode() == "BCInflow"
    assert bytes(back.field_data["cgns:bc_type:walls"]).decode() == "BCWall"


def _structured_cgns(path, dims=(3, 2, 2), with_solution=True, extra_zone=False):
    """A one-base, one-zone *Structured* file written with the writer's own node
    helpers: coordinates ``(nk, nj, ni)`` (``i`` fastest) and two solutions."""
    h5py = pytest.importorskip("h5py")
    ni, nj, nk = dims
    with h5py.File(path, "w") as f:
        _cgns._write_node_attrs(f, "HDF5 MotherNode", "Root Node of HDF5 File", "MT")
        base = _cgns._create_group(f, "Base")
        _cgns._write_node_attrs(base, "Base", "CGNSBase_t", "I4")
        base.create_dataset(" data", data=np.array([3, 3], dtype="<i4"))
        zone = _cgns._create_group(base, "Block1")
        _cgns._write_node_attrs(zone, "Block1", "Zone_t", "I4")
        cells = (ni - 1) * (nj - 1) * (nk - 1)
        zone.create_dataset(
            " data",
            data=np.array(
                [[ni, nj, nk], [ni - 1, nj - 1, nk - 1], [0, 0, 0]], dtype="<i4"
            ),
        )
        _cgns._c1_node(zone, "ZoneType", "ZoneType_t", "Structured")
        coords = _cgns._create_group(zone, "GridCoordinates")
        _cgns._write_node_attrs(coords, "GridCoordinates", "GridCoordinates_t", "MT")
        k, j, i = np.meshgrid(
            np.arange(nk), np.arange(nj), np.arange(ni), indexing="ij"
        )
        for name, values in (
            ("CoordinateX", i * 1.0),
            ("CoordinateY", j * 2.0),
            ("CoordinateZ", k * 3.0),
        ):
            node = _cgns._create_group(coords, name)
            _cgns._write_node_attrs(node, name, "DataArray_t", "R8")
            node.create_dataset(" data", data=values.astype("<f8"))
        if with_solution:
            sol = _cgns._create_group(zone, "Flow")
            _cgns._write_node_attrs(sol, "Flow", "FlowSolution_t", "MT")
            _cgns._c1_node(sol, "GridLocation", "GridLocation_t", "Vertex")
            node = _cgns._create_group(sol, "Temperature")
            _cgns._write_node_attrs(node, "Temperature", "DataArray_t", "R8")
            node.create_dataset(
                " data", data=np.arange(ni * nj * nk, dtype="<f8").reshape(nk, nj, ni)
            )
            csol = _cgns._create_group(zone, "CellFlow")
            _cgns._write_node_attrs(csol, "CellFlow", "FlowSolution_t", "MT")
            _cgns._c1_node(csol, "GridLocation", "GridLocation_t", "CellCenter")
            node = _cgns._create_group(csol, "Pressure")
            _cgns._write_node_attrs(node, "Pressure", "DataArray_t", "R8")
            node.create_dataset(
                " data",
                data=np.arange(cells, dtype="<f8").reshape(nk - 1, nj - 1, ni - 1),
            )
        if extra_zone:
            # a second, 2 x 2 x 2 zone with no solution
            zone2 = _cgns._create_group(base, "Block2")
            _cgns._write_node_attrs(zone2, "Block2", "Zone_t", "I4")
            zone2.create_dataset(
                " data", data=np.array([[2, 2, 2], [1, 1, 1], [0, 0, 0]], dtype="<i4")
            )
            _cgns._c1_node(zone2, "ZoneType", "ZoneType_t", "Structured")
            coords2 = _cgns._create_group(zone2, "GridCoordinates")
            _cgns._write_node_attrs(
                coords2, "GridCoordinates", "GridCoordinates_t", "MT"
            )
            k2, j2, i2 = np.meshgrid(
                np.arange(2), np.arange(2), np.arange(2), indexing="ij"
            )
            for name, values in (
                ("CoordinateX", i2 + 10.0),
                ("CoordinateY", j2 * 1.0),
                ("CoordinateZ", k2 * 1.0),
            ):
                node = _cgns._create_group(coords2, name)
                _cgns._write_node_attrs(node, name, "DataArray_t", "R8")
                node.create_dataset(" data", data=values.astype("<f8"))
            # an interface of the first zone naming its donor
            conn = _cgns._create_group(zone, "ZoneGridConnectivity")
            _cgns._write_node_attrs(
                conn, "ZoneGridConnectivity", "ZoneGridConnectivity_t", "MT"
            )
            link = _cgns._c1_node(conn, "to_block2", "GridConnectivity1to1_t", "Block2")
            rng = _cgns._create_group(link, "PointRange")
            _cgns._write_node_attrs(rng, "PointRange", "IndexRange_t", "I4")
            # the i = ni face of the first zone: i fixed at 3, all j, all k
            rng.create_dataset(
                " data", data=np.array([[3, 3], [1, 2], [1, 2]], dtype="<i4")
            )


@pytest.mark.parametrize("engine", ["python", "core"])
def test_structured_zone_is_expanded_to_a_lattice(tmp_path, engine):
    if engine == "core" and not _HAS_HDF5:
        pytest.skip("needs the C++ HDF5-enabled core")
    path = str(tmp_path / "s.cgns")
    _structured_cgns(path)
    mesh = (_cgns.read if engine == "python" else meshioplusplus.cgns.read)(path)
    assert len(mesh.points) == 12
    assert mesh.points[1].tolist() == [1.0, 0.0, 0.0]
    assert mesh.points[11].tolist() == [2.0, 2.0, 3.0]
    assert [b.type for b in mesh.cells] == ["hexahedron"]
    assert mesh.cells[0].data.shape == (2, 8)
    assert mesh.point_data["Temperature"].tolist() == list(map(float, range(12)))
    assert mesh.cell_data["Pressure"][0].tolist() == [0.0, 1.0]


@pytest.mark.parametrize("engine", ["python", "core"])
def test_rind_planes_of_a_solution_are_stripped(tmp_path, engine):
    h5py = pytest.importorskip("h5py")
    if engine == "core" and not _HAS_HDF5:
        pytest.skip("needs the C++ HDF5-enabled core")
    path = str(tmp_path / "r.cgns")
    _structured_cgns(path, with_solution=False)
    ni, nj, nk = 3, 2, 2
    padded = np.full((nk + 2, nj + 2, ni + 2), -1.0)
    padded[1:-1, 1:-1, 1:-1] = np.arange(ni * nj * nk, dtype=float).reshape(nk, nj, ni)
    with h5py.File(path, "r+") as f:
        zone = f["Base/Block1"]
        sol = _cgns._create_group(zone, "Flow")
        _cgns._write_node_attrs(sol, "Flow", "FlowSolution_t", "MT")
        _cgns._c1_node(sol, "GridLocation", "GridLocation_t", "Vertex")
        rind = _cgns._create_group(sol, "Rind")
        _cgns._write_node_attrs(rind, "Rind", "Rind_t", "I4")
        rind.create_dataset(" data", data=np.ones(6, dtype="<i4"))
        node = _cgns._create_group(sol, "Temperature")
        _cgns._write_node_attrs(node, "Temperature", "DataArray_t", "R8")
        node.create_dataset(" data", data=padded)
    mesh = (_cgns.read if engine == "python" else meshioplusplus.cgns.read)(path)
    assert mesh.point_data["Temperature"].tolist() == list(map(float, range(12)))


@pytest.mark.parametrize("engine", ["python", "core"])
def test_solutions_of_several_zones_are_concatenated(tmp_path, engine):
    if engine == "core" and not _HAS_HDF5:
        pytest.skip("needs the C++ HDF5-enabled core")
    path = str(tmp_path / "two.cgns")
    _structured_cgns(path, extra_zone=True)
    mesh = (_cgns.read if engine == "python" else meshioplusplus.cgns.read)(path)
    assert len(mesh.points) == 12 + 8
    assert [b.type for b in mesh.cells] == ["hexahedron", "hexahedron"]
    temperature = mesh.point_data["Temperature"]
    assert temperature[:12].tolist() == list(map(float, range(12)))
    assert np.isnan(temperature[12:]).all()
    pressure = mesh.cell_data["Pressure"]
    assert pressure[0].tolist() == [0.0, 1.0]
    assert np.isnan(pressure[1]).all() and len(pressure[1]) == 1
    # the first zone's interface: the points of the i = 3 face, by their lattice number
    by_name = {r.name: sorted(r.entries.tolist()) for r in mesh.regions}
    assert by_name["cgns:interface:to_block2"] == [2, 5, 8, 11]
    assert bytes(mesh.field_data["cgns:interface_donor:to_block2"]).decode() == "Block2"


@pytest.mark.parametrize("engine", ["python", "core"])
@pytest.mark.parametrize("which", ["Base2", "1", "Base", "0"])
def test_a_base_of_a_multi_base_file_can_be_chosen(tmp_path, engine, which):
    h5py = pytest.importorskip("h5py")
    if engine == "core" and not _HAS_HDF5:
        pytest.skip("needs the C++ HDF5-enabled core")
    path = str(tmp_path / "bases.cgns")
    _structured_cgns(path)
    with h5py.File(path, "r+") as f:
        f.copy(f["Base"], "Base2")
        node = f["Base2/Block1/GridCoordinates/CoordinateX/ data"]
        node[...] = node[...] + 100.0
    mesh = (
        (lambda p, b: _cgns.read(p, base=b))
        if engine == "python"
        else (lambda p, b: meshioplusplus.cgns.read(p, base=b))
    )(path, which)
    shifted = which in ("Base2", "1")
    assert (mesh.points[:, 0].min() >= 100.0) == shifted


def test_an_unknown_base_is_refused(tmp_path):
    path = str(tmp_path / "one.cgns")
    _structured_cgns(path)
    with pytest.raises(ReadError, match="base"):
        _cgns.read(path, base="nope")


def _series_mesh(scale, cell=False):
    pts = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], float)
    kwargs = (
        {"cell_data": {"p": [scale * np.array([1.0, 2.0])]}}
        if cell
        else {"point_data": {"T": scale * np.arange(4.0)}}
    )
    return meshioplusplus.Mesh(
        pts, [("triangle", np.array([[0, 1, 2], [1, 3, 2]]))], **kwargs
    )


@pytest.mark.parametrize("cell", [False, True])
@pytest.mark.parametrize("engine", ["core", "python"])
def test_series_writer_writes_iterative_data(tmp_path, engine, cell):
    if not _HAS_HDF5:  # reading a time_step back needs the core, for either writer
        pytest.skip("needs the C++ HDF5-enabled core")
    path = tmp_path / "series.cgns"
    cls = meshioplusplus.cgns.SeriesWriter if engine == "core" else _cgns.SeriesWriter
    with cls(str(path)) as writer:
        for t in (0.0, 0.5, 2.0):
            writer.write(t, _series_mesh(1.0 + t, cell))
    for k, t in enumerate((0.0, 0.5, 2.0)):
        mesh = meshioplusplus.cgns.read(str(path), time_step=k)
        if cell:
            np.testing.assert_allclose(
                mesh.cell_data["p"][0], (1.0 + t) * np.array([1.0, 2.0])
            )
        else:
            np.testing.assert_allclose(mesh.point_data["T"], (1.0 + t) * np.arange(4.0))


def test_series_writer_refuses_both_kinds_of_data(tmp_path):
    mesh = _series_mesh(1.0)
    mesh.cell_data["p"] = [np.array([1.0, 2.0])]
    with meshioplusplus.cgns.SeriesWriter(str(tmp_path / "s.cgns")) as writer:
        with pytest.raises(WriteError):
            writer.write(0.0, mesh)
        writer.write(0.0, _series_mesh(1.0))
