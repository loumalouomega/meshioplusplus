"""libMesh ``.xda``/``.xdr``: both engines, both encodings, the fixtures written
by ``tools/gen_libmesh_fixtures.py`` in libMesh's own node numbering."""

import pathlib
import shutil
import struct

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus._facets import facet_nodes
from meshioplusplus.libmesh import _libmesh as py_libmesh

from .test_node_order import _is_valid

MESHES = pathlib.Path(__file__).parent / "meshes" / "libmesh"
FIXTURES = sorted(MESHES.glob("*.xd[ar]"))


@pytest.fixture(params=["core", "python"])
def read(request):
    if request.param == "core":
        return meshioplusplus.libmesh.read
    return py_libmesh.read


def _regions(mesh):
    return sorted(
        (r.kind, r.name, r.dim, r.tag, tuple(np.asarray(r.entries).ravel()))
        for r in mesh.regions
    )


def _same(a, b):
    np.testing.assert_array_equal(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_array_equal(x, y)
    assert sorted(a.point_data) == sorted(b.point_data)
    for name in a.point_data:
        np.testing.assert_array_equal(a.point_data[name], b.point_data[name])
    assert _regions(a) == _regions(b)


@pytest.mark.parametrize("path", FIXTURES, ids=[p.name for p in FIXTURES])
def test_engines_agree_on_every_fixture(path):
    _same(meshioplusplus.libmesh.read(path), py_libmesh.read(path))


@pytest.mark.parametrize("stem", sorted({p.stem for p in FIXTURES}))
def test_ascii_and_xdr_read_the_same(read, stem):
    _same(read(MESHES / f"{stem}.xda"), read(MESHES / f"{stem}.xdr"))


def _side_points(mesh, region):
    return [
        np.asarray(mesh.points)[facet_nodes(mesh, int(c), int(f))[1]]
        for c, f in np.asarray(region.entries)
    ]


def test_hex27_in_meshio_order_with_regions(read):
    mesh = read(MESHES / "hex27.xdr")
    assert [c.type for c in mesh.cells] == ["hexahedron27"]
    for row in np.asarray(mesh.cells[0].data):
        assert _is_valid(np.asarray(mesh.points)[row], "hexahedron27")
    np.testing.assert_array_equal(mesh.cell_data["libmesh:subdomain"][0], [1, 2])
    regions = {(r.kind, r.name): r for r in mesh.regions}
    assert regions[("cell", "left")].tag == 1
    assert regions[("cell", "subdomain_2")].tag == 2
    inlet = regions[("side", "inlet")]
    outlet = regions[("side", "boundary_2")]
    assert inlet.tag == 1 and outlet.tag == 2 and inlet.dim == 2
    (face,) = _side_points(mesh, inlet)
    np.testing.assert_array_equal(face[:, 0], 0.0)
    (face,) = _side_points(mesh, outlet)
    np.testing.assert_array_equal(face[:, 0], 2.0)
    wall = regions[("point", "wall")]
    assert wall.tag == 5 and len(wall.entries) == 9
    np.testing.assert_array_equal(np.asarray(mesh.points)[wall.entries][:, 0], 0.0)


def test_legacy_header_and_six_sides(read):
    mesh = read(MESHES / "one_hex.xda")
    assert [c.type for c in mesh.cells] == ["hexahedron"]
    assert _is_valid(np.asarray(mesh.points)[mesh.cells[0].data[0]], "hexahedron")
    sides = {r.tag: r for r in mesh.regions if r.kind == "side"}
    assert sorted(sides) == list(range(6))
    facets = {int(r.entries[0, 1]) for r in sides.values()}
    assert len(facets) == 6
    # libMesh side 0 is the bottom face (z = -1), side 5 the top.
    np.testing.assert_array_equal(_side_points(mesh, sides[0])[0][:, 2], -1.0)
    np.testing.assert_array_equal(_side_points(mesh, sides[5])[0][:, 2], 1.0)


def test_refined_mesh_keeps_the_active_children(read):
    mesh = read(MESHES / "amr_quad.xdr")
    assert [(c.type, len(c.data)) for c in mesh.cells] == [("quad", 4)]
    np.testing.assert_array_equal(mesh.cell_data["libmesh:level"][0], [1, 1, 1, 1])
    np.testing.assert_array_equal(mesh.cell_data["libmesh:p_level"][0], [1, 1, 2, 0])
    # Node 9 is unused (NaN coordinates in the file): dropped, ids kept.
    np.testing.assert_array_equal(mesh.point_data["libmesh:id"], np.arange(9))
    assert np.isfinite(mesh.points).all()
    (bottom,) = [r for r in mesh.regions if r.kind == "side"]
    assert (bottom.name, bottom.tag, len(bottom.entries)) == ("bottom", 3, 2)
    for edge in _side_points(mesh, bottom):
        np.testing.assert_array_equal(edge[:, 1], 0.0)


def test_extra_nodes_are_dropped(read):
    mesh = read(MESHES / "tet14_prism18.xda")
    assert [c.type for c in mesh.cells] == ["tetra10", "wedge18"]
    for block in mesh.cells:
        for row in np.asarray(block.data):
            assert _is_valid(np.asarray(mesh.points)[row], block.type)


def test_legacy_files_are_refused(read, tmp_path):
    path = tmp_path / "old.xda"
    path.write_text("DEAL 003:003\n1\t # Num. Elements\n")
    with pytest.raises(meshioplusplus.ReadError, match="legacy"):
        read(path)


def test_truncated_xdr_is_refused(read, tmp_path):
    data = (MESHES / "hex27.xdr").read_bytes()
    path = tmp_path / "cut.xdr"
    path.write_bytes(data[: len(data) // 2])
    with pytest.raises(meshioplusplus.ReadError):
        read(path)


def test_xdr_starts_with_a_padded_version_string():
    data = (MESHES / "hex27.xdr").read_bytes()
    (length,) = struct.unpack(">I", data[:4])
    assert data[4 : 4 + length] == b"libMesh-1.3.0" and length % 4 != 0


@pytest.mark.parametrize("name", ["hex27.xda", "hex27.xdr"])
def test_sniffed_without_the_extension(tmp_path, name):
    path = tmp_path / "mesh.data"
    shutil.copy(MESHES / name, path)
    assert meshioplusplus.sniff_format(path) == "libmesh"
    assert _core.sniff_format(str(path)) == "libmesh"


def test_extension_dispatch():
    mesh = meshioplusplus.read(MESHES / "hex27.xda")
    assert mesh.cells[0].type == "hexahedron27"


def test_edge_and_shellface_sets(read):
    mesh = read(MESHES / "edges_shell.xda")
    assert [c.type for c in mesh.cells] == ["hexahedron20", "quad8", "line3"]
    # The hex's edge 0 and its top edge (named twice, by the hex and the shell).
    assert len(mesh.cells[2].data) == 2
    np.testing.assert_array_equal(mesh.cell_data["libmesh:subdomain"][2], [-1, -1])
    regions = {(r.kind, r.name): r for r in mesh.regions}
    axis = regions[("cell", "axis:edge")]
    top = regions[("cell", "boundary_12:edge")]
    assert (axis.tag, axis.dim, top.tag) == (11, 1, 12)
    for region, z in ((axis, 0.0), (top, 1.0)):
        (cell,) = np.asarray(region.entries)
        line = mesh.cells[2].data[cell - 2]
        np.testing.assert_array_equal(mesh.points[line][:, 1:], [[0, z]] * 3)
        assert mesh.points[line[2], 0] == 0.5  # the mid-edge node
    assert np.asarray(regions[("cell", "front:shellface0")].entries).tolist() == [1]
    assert regions[("cell", "boundary_21:shellface1")].tag == 21


@pytest.mark.parametrize("name", ["hex27.xda.gz", "hex27.xdr.bz2"])
def test_compressed_files(read, name):
    _same(read(MESHES / name), read(MESHES / "hex27.xda"))


def test_core_leaves_bzip2_to_python():
    with pytest.raises(Exception, match="bzip2"):
        _core.libmesh_read(str(MESHES / "hex27.xdr.bz2"))
    assert meshioplusplus.libmesh.read(MESHES / "hex27.xdr.bz2").cells


@pytest.mark.parametrize("stem", ["hex27", "one_hex", "edges_shell"])
@pytest.mark.parametrize("ext", [".xda", ".xdr"])
def test_writer_round_trips(tmp_path, stem, ext):
    mesh = py_libmesh.read(MESHES / f"{stem}{ext}")
    core, python = tmp_path / f"core{ext}", tmp_path / f"python{ext}"
    _core.libmesh_write(str(core), mesh)
    py_libmesh.write(python, mesh)
    assert core.read_bytes() == python.read_bytes()
    _same(py_libmesh.read(core), mesh)
    _same(meshioplusplus.libmesh.read(core), mesh)


def test_writer_matches_libmesh_layout(tmp_path):
    # libMesh-1.8.0's own header lines and legend (XdrIO::write).
    mesh = py_libmesh.read(MESHES / "hex27.xda")
    path = tmp_path / "mesh.xda"
    meshioplusplus.write(path, mesh)
    lines = path.read_text().splitlines()
    assert lines[:4] == [
        "libMesh-1.8.0",
        "2\t # number of elements",
        "45\t # number of nodes",
        ".\t # boundary condition specification file",
    ]
    assert "2\t # n_elem at level 0, [ type sid (n0 ... nN-1) ]" in lines
    assert lines[-1] == "0\t # number of shellface boundary conditions"


def test_writer_maps_foreign_regions(tmp_path, capfd):
    from meshioplusplus._regions import Region

    points = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0], [2, 1, 0]]
    )
    mesh = meshioplusplus.Mesh(
        points.astype(float),
        [("quad", np.array([[0, 1, 2, 3], [1, 4, 5, 2]])), ("polygon", [[0, 1, 2]])],
    )
    mesh.regions = [
        Region("left", "cell", np.array([0]), 2, -1),
        Region("right", "cell", np.array([1]), 2, 3),
        Region("wall", "side", np.array([[0, 3]]), 1, -1),
        Region("pin", "point", np.array([4]), -1, -1),
    ]
    for ext in (".xda", ".xdr"):
        path = tmp_path / f"mesh{ext}"
        py_libmesh.write(path, mesh)
        assert "polygon" in capfd.readouterr().err
        back = meshioplusplus.libmesh.read(path)
        assert [c.type for c in back.cells] == ["quad"]
        np.testing.assert_array_equal(back.cell_data["libmesh:subdomain"][0], [4, 3])
        got = {(r.kind, r.name, r.tag) for r in back.regions}
        assert got == {
            ("cell", "left", 4),
            ("cell", "right", 3),
            ("side", "wall", 0),
            ("point", "pin", 0),
        }
        wall = next(r for r in back.regions if r.name == "wall")
        assert np.asarray(wall.entries).tolist() == [[0, 3]]


@pytest.mark.parametrize("suffix", [".xda.gz", ".xdr.bz2"])
def test_writer_compresses(tmp_path, suffix):
    mesh = py_libmesh.read(MESHES / "hex27.xda")
    path = tmp_path / f"mesh{suffix}"
    meshioplusplus.write(path, mesh)
    first = path.read_bytes()
    meshioplusplus.write(path, mesh)
    assert path.read_bytes() == first  # no timestamp in the gzip header
    _same(meshioplusplus.read(path), mesh)


def test_writer_sparse_node_ids(tmp_path):
    mesh = py_libmesh.read(MESHES / "hex27.xda")
    mesh.point_data["libmesh:id"] = 2 * np.arange(len(mesh.points), dtype=np.int64)
    for ext in (".xda", ".xdr"):
        core, python = tmp_path / f"core{ext}", tmp_path / f"python{ext}"
        _core.libmesh_write(str(core), mesh)
        py_libmesh.write(python, mesh)
        assert core.read_bytes() == python.read_bytes()
        back = meshioplusplus.libmesh.read(core)
        np.testing.assert_array_equal(
            back.point_data["libmesh:id"], mesh.point_data["libmesh:id"]
        )
    # ASCII writes the unused ids as 0 (libMesh cannot read back `nan`).
    assert "nan" not in (tmp_path / "core.xda").read_text()
