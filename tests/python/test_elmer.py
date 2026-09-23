"""Elmer mesh directories: both engines, the fixtures written by
``tools/gen_elmer_fixtures.py`` in Elmer's own node numbering, partitions, and
the directory sniff that lets a bare directory be read with no format named."""

import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.elmer import _elmer as py_elmer

from .test_node_order import _is_valid

MESHES = pathlib.Path(__file__).parent / "meshes" / "elmer"
SERIAL = ["tet10_two_bodies", "hex27_block", "quad8_2d"]


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read/write signatures."""
    if request.param == "core":
        return meshioplusplus.elmer
    return py_elmer


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
    assert _regions(a) == _regions(b)
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for name in a.cell_data:
        for x, y in zip(a.cell_data[name], b.cell_data[name]):
            np.testing.assert_array_equal(x, y)


@pytest.mark.parametrize("name", SERIAL + ["partitioned"])
def test_engines_read_the_same_mesh(name):
    _same(meshioplusplus.elmer.read(MESHES / name), py_elmer.read(MESHES / name))


@pytest.mark.parametrize("name", SERIAL)
def test_quadratic_cells_are_valid_in_meshio_order(name):
    """Every mid-edge node on its edge, every centre on its face, every cell
    positively oriented -- so the Elmer node order was undone correctly."""
    mesh = meshioplusplus.read(MESHES / name)
    bulk = mesh.cells[0].type  # bulk blocks come first; boundary faces in 3-D
    for block in mesh.cells:  # have no +z normal for _is_valid to check
        if block.type != bulk:
            continue
        for row in block.data:
            assert _is_valid(mesh.points[row], block.type), (name, block.type)


def test_bodies_and_boundaries_are_regions(engine):
    mesh = engine.read(MESHES / "tet10_two_bodies")
    assert [(c.type, len(c.data)) for c in mesh.cells] == [
        ("tetra10", 2),
        ("triangle6", 2),
        ("vertex", 1),
    ]
    got = {r.name: (r.kind, r.dim, r.tag, r.entries.tolist()) for r in mesh.regions}
    assert got == {
        "left": ("cell", 3, 1, [0]),
        "right": ("cell", 3, 2, [1]),
        "outer": ("cell", 2, 2, [2]),
        "interface": ("cell", 2, 3, [3]),
        "corner": ("cell", 0, 4, [4]),
    }
    # Without mesh.names, ids name the regions.
    names = {r.name: r.tag for r in engine.read(MESHES / "quad8_2d").regions}
    assert names == {"body_1": 1, "boundary_1": 1, "boundary_2": 2}


def test_mesh_header_stands_for_its_directory():
    a = meshioplusplus.read(MESHES / "hex27_block")
    b = meshioplusplus.read(MESHES / "hex27_block" / "mesh.header")
    _same(a, b)


def test_partitions_merge_with_labels(engine):
    mesh = engine.read(MESHES / "partitioned")
    assert len(mesh.points) == 4
    # Element 1's halo copy in part 2 is read once, owned by part 1.
    assert [(c.type, len(c.data)) for c in mesh.cells] == [
        ("triangle", 2),
        ("line", 2),
    ]
    assert [a.tolist() for a in mesh.cell_data["partition:part"]] == [[0, 1], [0, 1]]
    # mesh.names next to the partitioning directory names the regions.
    assert sorted(r.name for r in mesh.regions) == ["bottom", "plate", "right"]
    _same(mesh, engine.read(MESHES / "partitioned" / "partitioning.2"))
    second = engine.read(MESHES / "partitioned", piece=-1)
    assert len(second.cells[0].data) == 2  # its own element and the halo copy
    with pytest.raises(meshioplusplus.ReadError):
        engine.read(MESHES / "partitioned", piece=2)


def test_piece_reaches_the_reader_through_read():
    mesh = meshioplusplus.read(MESHES / "partitioned", piece=0)
    assert [(c.type, len(c.data)) for c in mesh.cells] == [
        ("triangle", 1),
        ("line", 1),
    ]


def test_writers_write_the_same_bytes(tmp_path):
    mesh = meshioplusplus.read(MESHES / "tet10_two_bodies")
    _core.elmer_write(str(tmp_path / "cpp"), mesh)
    py_elmer.write(tmp_path / "py", mesh)
    for name in ("mesh.header", "mesh.nodes", "mesh.elements", "mesh.boundary"):
        assert (tmp_path / "cpp" / name).read_bytes() == (
            tmp_path / "py" / name
        ).read_bytes(), name
    # mesh.names differs only in the provenance block's timestamp, if any.
    strip = [
        line
        for line in (tmp_path / "cpp" / "mesh.names").read_text().splitlines()
        if not line.startswith("! ") or "names for" in line
    ]
    assert strip == [
        line
        for line in (tmp_path / "py" / "mesh.names").read_text().splitlines()
        if not line.startswith("! ") or "names for" in line
    ]


@pytest.mark.parametrize("name", SERIAL)
def test_round_trip(engine, name, tmp_path):
    mesh = meshioplusplus.read(MESHES / name)
    engine.write(tmp_path / "out", mesh)
    back = engine.read(tmp_path / "out")
    np.testing.assert_allclose(back.points, mesh.points)
    _same(back, mesh)


def test_parents_are_regenerated(tmp_path):
    mesh = meshioplusplus.read(MESHES / "tet10_two_bodies")
    meshioplusplus.write(tmp_path / "out", mesh, file_format="elmer")
    rows = [
        line.split()[:5]
        for line in (tmp_path / "out" / "mesh.boundary").read_text().splitlines()
    ]
    # outer: one parent; interface: both tets; the corner point: its first tet.
    assert rows == [
        ["1", "2", "1", "0", "306"],
        ["2", "3", "1", "2", "306"],
        ["3", "4", "1", "0", "101"],
    ]


def test_side_regions_become_boundary_elements(engine, tmp_path):
    mesh = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )
    mesh.regions = [meshioplusplus.Region("wall", "side", [[0, 3]], dim=2, tag=5)]
    engine.write(tmp_path / "out", mesh)
    assert (tmp_path / "out" / "mesh.boundary").read_text() == "1 5 1 0 303 1 3 2\n"
    back = engine.read(tmp_path / "out")
    assert {(r.name, r.kind, r.dim, r.tag) for r in back.regions} == {
        ("body_1", "cell", 3, 1),
        ("wall", "cell", 2, 5),
    }


def test_directory_is_sniffed_but_write_needs_the_format(tmp_path):
    for name in SERIAL + ["partitioned"]:
        assert meshioplusplus.sniff_format(MESHES / name) == "elmer"
    mesh = meshioplusplus.read(MESHES / "quad8_2d")
    with pytest.raises((meshioplusplus.ReadError, meshioplusplus.WriteError)):
        meshioplusplus.write(tmp_path / "out", mesh)
    meshioplusplus.write(tmp_path / "out", mesh, file_format="elmer")
    assert meshioplusplus.read(tmp_path / "out").cells[0].type == "quad8"


def test_unsupported_things(engine, tmp_path):
    bad = tmp_path / "bad"
    bad.mkdir()
    (bad / "mesh.header").write_text("2 0 1\n1\n102 1\n")
    (bad / "mesh.nodes").write_text("1 -1 0 0 0\n2 -1 1 0 0\n")
    (bad / "mesh.elements").write_text("")
    (bad / "mesh.boundary").write_text("1 1 0 0 102 1 2\n")
    with pytest.raises(meshioplusplus.ReadError, match="102"):
        engine.read(bad)
    assert len(engine.read(bad, lenient=True).cells) == 0
    binary = tmp_path / "binary"
    binary.mkdir()
    for name in ("mesh.header", "mesh.nodes.bin", "mesh.elements.bin"):
        (binary / name).write_text("")
    with pytest.raises(meshioplusplus.ReadError, match="binary"):
        engine.read(binary)
    poly = meshioplusplus.Mesh(np.zeros((5, 3)), [("polygon", [np.arange(5)])])
    with pytest.raises(meshioplusplus.WriteError):
        engine.write(tmp_path / "poly", poly)


def test_dropped_things_are_noted(tmp_path, capfd):
    mesh = meshioplusplus.read(MESHES / "quad8_2d")
    mesh.point_data["t"] = np.arange(len(mesh.points), dtype=float)
    mesh.regions.append(meshioplusplus.Region("pin", "point", [0]))
    meshioplusplus.write(tmp_path / "out", mesh, file_format="elmer")
    err = capfd.readouterr().err
    assert "point region" in err and "data" in err
    names = (tmp_path / "out" / "mesh.names").read_text()
    assert "regions-dropped" in names and "data-dropped" in names
    assert "$" not in "".join(
        line for line in names.splitlines() if line.startswith("! ")
    )
