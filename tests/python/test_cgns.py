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
