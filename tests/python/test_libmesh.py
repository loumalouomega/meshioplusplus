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
