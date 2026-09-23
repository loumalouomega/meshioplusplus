import io
import pathlib

import numpy as np
import pytest

import meshioplusplus

from . import helpers


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.empty_mesh,
        helpers.tri_mesh,
        helpers.tri_mesh_2d,
        helpers.triangle6_mesh,
        helpers.quad_mesh,
        helpers.quad8_mesh,
        helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
    ],
)
def test(mesh, tmp_path):
    helpers.write_read(
        tmp_path,
        meshioplusplus.nastran.write,
        meshioplusplus.nastran.read,
        mesh,
        1.0e-13,
    )


@pytest.mark.parametrize("filename", ["cylinder.fem", "cylinder_cells_first.fem"])
def test_reference_file(filename):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "nastran" / filename

    mesh = meshioplusplus.read(filename)

    # points
    assert np.isclose(mesh.points.sum(), 16.5316866)

    # cells
    ref_num_cells = {
        "line": 241,
        "triangle": 171,
        "quad": 721,
        "pyramid": 1180,
        "tetra": 5309,
    }
    assert {
        cell_block.type: cell_block.data.sum() for cell_block in mesh.cells
    } == ref_num_cells


def test_long_format():
    filename = io.StringIO(
        "BEGIN BULK\n"
        "GRID*    43                             1.50000000000000 0.0\n"
        "*        0.\n"
        "ENDDATA\n"
    )

    mesh = meshioplusplus.read(filename, "nastran")

    # points
    assert len(mesh.points) == 1
    assert np.isclose(mesh.points.sum(), 1.5)


# --- OptiStruct / HyperMesh decks, both engines -------------------------------
#
# ``optistruct_mixed.fem`` is written by ``tools/gen_optistruct_fixture.py`` (not
# by meshio++); the cylinder decks were exported by HyperMesh 2017.3 and
# ``composite_plate_2022.fem`` comes from pyNastran (BSD-3-Clause).

from meshioplusplus import _core  # noqa: E402
from meshioplusplus.nastran import _nastran as py_nastran  # noqa: E402

MESHES = pathlib.Path(__file__).resolve().parent / "meshes" / "nastran"
DECKS = [
    "cylinder.fem",
    "cylinder_cells_first.fem",
    "optistruct_mixed.fem",
    "composite_plate_2022.fem",
]


@pytest.fixture(params=["core", "python"])
def read(request):
    return _core.nastran_read if request.param == "core" else py_nastran.read


def _regions(mesh):
    return {(r.kind, r.name, r.dim, r.tag): r.entries.tolist() for r in mesh.regions}


def test_hypermesh_component_is_a_region(read):
    # The done-when of roadmap §1.1: a HyperMesh-exported deck reads with its
    # component as a named region ($HMNAME COMP comes after the elements).
    mesh = read(str(MESHES / "cylinder.fem"))
    num_cells = sum(len(c) for c in mesh.cells)
    assert _regions(mesh) == {("cell", "misc1", -1, 1): list(range(num_cells))}


def test_optistruct_mixed(read, capfd):
    mesh = read(str(MESHES / "optistruct_mixed.fem"))
    assert [c.type for c in mesh.cells] == [
        "hexahedron20",
        "tetra10",
        "quad",
        "triangle",
        "line",
    ]
    assert _regions(mesh) == {
        ("cell", "solids", 3, 10): [0, 1],  # $HMMOVE with THRU
        ("cell", "shells", 2, 20): [2, 3],  # $HMMOVE with single ids
        ("cell", "beams", 1, 30): [4],  # through its property id only
        ("cell", "empty", -1, 40): [],  # named, no members
        ("point", "base_nodes", -1, 1): [0, 1],  # SET GRID with THRU, $HMSET
        ("cell", "set_2", -1, 2): [0, 2, 4],  # free-field SET ELEM, no name
    }
    assert mesh.cell_data["nastran:ref"][0].tolist() == [1]
    assert mesh.cell_data["nastran:ref"][4].tolist() == [7]
    err = " ".join(capfd.readouterr().err.split())
    assert "CONM2 (1), CONTACT (1), DESVAR (1), DRESP1 (2), TIE (1)" in err
    assert "PSOLID" not in err and "MAT1" not in err
    # quadratic node order: every mid-edge node on its edge's midpoint
    p = mesh.points
    hexa = mesh.cells[0].data[0]
    for k, (a, b) in enumerate(
        [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
        + [(0, 4), (1, 5), (2, 6), (3, 7)]
    ):
        assert np.allclose(p[hexa[8 + k]], (p[hexa[a]] + p[hexa[b]]) / 2)
    tet = mesh.cells[1].data[0]
    for k, (a, b) in enumerate([(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]):
        assert np.allclose(p[tet[4 + k]], (p[tet[a]] + p[tet[b]]) / 2)


def test_optistruct_sets_and_property_components(read):
    mesh = read(str(MESHES / "composite_plate_2022.fem"))
    regions = _regions(mesh)
    assert regions[("cell", "Plate", 2, 1)] == list(range(100))
    assert len(regions[("cell", "ALL", 2, 3)]) == 100
    assert len(regions[("point", "hinge_edges", -1, 2)]) == 40
    assert len(regions[("point", "center_point_load", -1, 1)]) == 1


@pytest.mark.parametrize("deck", DECKS)
def test_engines_agree(deck):
    a = _core.nastran_read(str(MESHES / deck))
    b = py_nastran.read(str(MESHES / deck))
    assert np.array_equal(a.points, b.points)
    assert [(c.type, c.data.tolist()) for c in a.cells] == [
        (c.type, c.data.tolist()) for c in b.cells
    ]
    assert _regions(a) == _regions(b)
    for name in set(a.cell_data) | set(b.cell_data):
        assert [x.tolist() for x in a.cell_data[name]] == [
            x.tolist() for x in b.cell_data[name]
        ]
    assert {k: v.tolist() for k, v in a.point_data.items()} == {
        k: v.tolist() for k, v in b.point_data.items()
    }


def test_errors(read, tmp_path):
    path = tmp_path / "bad.bdf"
    path.write_text("GRID    1               0.0     0.0     0.0\n")
    with pytest.raises(meshioplusplus.ReadError, match="BEGIN BULK"):
        read(str(path))
    path.write_text("BEGIN BULK\nCTETRA  1       1       1       2       3\nENDDATA\n")
    with pytest.raises(meshioplusplus.ReadError, match="3 nodes"):
        read(str(path))
    path.write_text("BEGIN BULK\nCTRIA3  1       1       1       2       3\nENDDATA\n")
    with pytest.raises(meshioplusplus.ReadError, match="grid 1"):
        read(str(path))


@pytest.mark.parametrize(
    "write", [_core.nastran_write, py_nastran.write], ids=["core", "python"]
)
def test_regions_round_trip_as_components(write, tmp_path, capfd):
    tet10 = helpers.tet10_mesh.cells[0].data
    mesh = meshioplusplus.Mesh(
        helpers.tet10_mesh.points, [("tetra10", np.vstack([tet10, tet10]))]
    )
    mesh.cell_data["nastran:ref"] = [np.array([3, 0], dtype=np.int64)]
    Region = meshioplusplus.Region
    mesh.regions = [
        Region("both", "cell", np.array([0, 1]), 3, 6),  # overlaps: dropped
        Region("first", "cell", np.array([0]), 3, 5),
        Region("nodes", "point", np.array([1]), -1, 2),  # points: dropped
    ]
    path = tmp_path / "out.bdf"
    write(str(path), mesh)
    text = path.read_text()
    assert "CTETRA_" not in text and "CTETRA  " in text
    # regions are taken in (kind, name, dim, tag) order: "both" wins
    assert '$HMNAME COMP                   6"both"' in text
    err = " ".join(capfd.readouterr().err.split())
    assert "'first' overlaps" in err and "'nodes' dropped" in err
    back = meshioplusplus.nastran.read(str(path))
    assert _regions(back) == {("cell", "both", 3, 6): [0, 1]}
    assert back.cell_data["nastran:ref"][0].tolist() == [3, 0]
    mesh.regions = [Region("bad", "cell", np.array([2]))]
    with pytest.raises(meshioplusplus.WriteError, match="names cell 2 of 2"):
        write(str(path), mesh)


def test_hypermesh_blocks_are_identical(tmp_path):
    mesh = meshioplusplus.Mesh(helpers.tri_mesh.points, helpers.tri_mesh.cells)
    mesh.regions = [meshioplusplus.Region('a "q"', "cell", np.array([0, 1]))]
    _core.nastran_write(str(tmp_path / "c.bdf"), mesh)
    py_nastran.write(str(tmp_path / "p.bdf"), mesh)

    def block(p):
        lines = p.read_text().splitlines()
        return lines[lines.index("$") :]

    assert block(tmp_path / "c.bdf") == block(tmp_path / "p.bdf")
    assert block(tmp_path / "c.bdf")[-2] == "$HMNAME COMP                   1\"a 'q'\""
