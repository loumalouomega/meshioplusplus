import copy
import pathlib
from functools import partial

import numpy as np
import pytest

import meshioplusplus

from . import helpers


def gmsh_periodic():
    mesh = copy.deepcopy(helpers.quad_mesh)
    trns = [0] * 16  # just for io testing
    mesh.gmsh_periodic = [
        [0, (3, 1), None, [[2, 0]]],
        [0, (4, 6), None, [[3, 5]]],
        [1, (2, 1), trns, [[5, 0], [4, 1], [4, 2]]],
    ]
    return mesh


@pytest.mark.parametrize(
    "mesh",
    [
        # helpers.empty_mesh,
        helpers.line_mesh,
        helpers.tri_mesh,
        helpers.triangle6_mesh,
        helpers.quad_mesh,
        helpers.quad8_mesh,
        # helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
        helpers.add_point_data(helpers.tri_mesh, 1),
        helpers.add_point_data(helpers.tri_mesh, 3),
        helpers.add_point_data(helpers.tri_mesh, 9),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (3,), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (9,), np.float64)]),
        helpers.add_field_data(helpers.tri_mesh, [1, 2], int),
        helpers.add_field_data(helpers.tet_mesh, [1, 3], int),
        helpers.add_field_data(helpers.hex_mesh, [1, 3], int),
        gmsh_periodic(),
    ],
)
@pytest.mark.parametrize("binary", [False, True])
def test_gmsh22(mesh, binary, tmp_path):
    writer = partial(meshioplusplus.gmsh.write, fmt_version="2.2", binary=binary)
    helpers.write_read(tmp_path, writer, meshioplusplus.gmsh.read, mesh, 1.0e-15)


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tri_mesh,
        helpers.triangle6_mesh,
        helpers.quad_mesh,
        helpers.quad8_mesh,
        # helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
        helpers.add_point_data(helpers.tri_mesh, 1),
        helpers.add_point_data(helpers.tri_mesh, 3),
        helpers.add_point_data(helpers.tri_mesh, 9),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (3,), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (9,), np.float64)]),
        helpers.add_field_data(helpers.tri_mesh, [1, 2], int),
        helpers.add_field_data(helpers.tet_mesh, [1, 3], int),
        helpers.add_field_data(helpers.hex_mesh, [1, 3], int),
    ],
)
@pytest.mark.parametrize("binary", [False, True])
def test_gmsh40(mesh, binary, tmp_path):
    writer = partial(meshioplusplus.gmsh.write, fmt_version="4.0", binary=binary)

    helpers.write_read(tmp_path, writer, meshioplusplus.gmsh.read, mesh, 1.0e-15)


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tri_mesh,
        helpers.triangle6_mesh,
        helpers.quad_mesh,
        helpers.quad8_mesh,
        # helpers.tri_quad_mesh,
        helpers.tet_mesh,
        helpers.tet10_mesh,
        helpers.hex_mesh,
        helpers.hex20_mesh,
        helpers.add_point_data(helpers.tri_mesh, 1),
        helpers.add_point_data(helpers.tri_mesh, 3),
        helpers.add_point_data(helpers.tri_mesh, 9),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (3,), np.float64)]),
        helpers.add_cell_data(helpers.tri_mesh, [("a", (9,), np.float64)]),
        helpers.add_field_data(helpers.tri_mesh, [1, 2], int),
        helpers.add_field_data(helpers.tet_mesh, [1, 3], int),
        helpers.add_field_data(helpers.hex_mesh, [1, 3], int),
        gmsh_periodic(),
    ],
)
@pytest.mark.parametrize("binary", [False, True])
def test_gmsh41(mesh, binary, tmp_path):
    writer = partial(meshioplusplus.gmsh.write, fmt_version="4.1", binary=binary)
    helpers.write_read(tmp_path, writer, meshioplusplus.gmsh.read, mesh, 1.0e-15)


def test_generic_io(tmp_path):
    helpers.generic_io(tmp_path / "test.msh")
    # With additional, insignificant suffix:
    helpers.generic_io(tmp_path / "test.0.msh")


@pytest.mark.parametrize(
    "filename, ref_sum, ref_num_cells",
    [("insulated-2.2.msh", 2.001762136876221, [21, 111])],
)
@pytest.mark.parametrize("binary", [False, True])
def test_reference_file(filename, ref_sum, ref_num_cells, binary, tmp_path):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "msh" / filename
    mesh = meshioplusplus.read(filename)
    tol = 1.0e-2
    s = mesh.points.sum()
    assert abs(s - ref_sum) < tol * ref_sum
    assert [c.type for c in mesh.cells] == ["line", "triangle"]
    assert [len(c.data) for c in mesh.cells] == ref_num_cells
    assert list(map(len, mesh.cell_data["gmsh:geometrical"])) == ref_num_cells
    assert list(map(len, mesh.cell_data["gmsh:physical"])) == ref_num_cells

    writer = partial(meshioplusplus.gmsh.write, fmt_version="2.2", binary=binary)
    helpers.write_read(tmp_path, writer, meshioplusplus.gmsh.read, mesh, 1.0e-15)


@pytest.mark.parametrize(
    "filename, ref_sum, ref_num_cells, ref_num_cells_in_cell_sets",
    [
        (
            "insulated-4.1.msh",
            2.001762136876221,
            {"line": 21, "triangle": 111},
            {"line": 27, "triangle": 120},
        )
    ],
    # Note that testing on number of cells in
    # cell_sets_dict will count both cells associated with physical tags, and
    # bounding entities.
)
@pytest.mark.parametrize("binary", [False, True])
def test_reference_file_with_entities(
    filename, ref_sum, ref_num_cells, ref_num_cells_in_cell_sets, binary, tmp_path
):
    this_dir = pathlib.Path(__file__).resolve().parent
    filename = this_dir / "meshes" / "msh" / filename

    mesh = meshioplusplus.read(filename)
    tol = 1.0e-2
    s = mesh.points.sum()
    assert abs(s - ref_sum) < tol * ref_sum
    assert {k: len(v) for k, v in mesh.cells_dict.items()} == ref_num_cells
    assert {
        k: len(v) for k, v in mesh.cell_data_dict["gmsh:physical"].items()
    } == ref_num_cells

    writer = partial(meshioplusplus.gmsh.write, fmt_version="4.1", binary=binary)

    num_cells = {k: 0 for k in ref_num_cells_in_cell_sets}
    for vv in mesh.cell_sets_dict.values():
        for k, v in vv.items():
            num_cells[k] += len(v)
    assert num_cells == ref_num_cells_in_cell_sets

    # $Entities is the only place 4.1 records physical-group membership, so the
    # regions it yields carry the group's real dimension and tag -- not the -1
    # placeholders a set-derived region has.
    assert sorted((r.name, r.dim, r.tag, len(r.entries)) for r in mesh.regions) == [
        ("convection", 1, 3, 21),
        ("insulation", 2, 2, 66),
        ("wire", 2, 1, 45),
    ]

    helpers.write_read(tmp_path, writer, meshioplusplus.gmsh.read, mesh, 1.0e-15)


_ENTITY_FILES = ["tests/python/meshes/msh/insulated-4.1.msh", "example/example.msh"]


def _assert_same_mesh(a, b):
    np.testing.assert_allclose(a.points, b.points)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for ca, cb in zip(a.cells, b.cells):
        np.testing.assert_array_equal(ca.data, cb.data)
    assert sorted(a.point_data) == sorted(b.point_data)
    for k in a.point_data:
        np.testing.assert_array_equal(a.point_data[k], b.point_data[k])
    assert sorted(a.cell_data) == sorted(b.cell_data)
    for k in a.cell_data:
        for x, y in zip(a.cell_data[k], b.cell_data[k]):
            np.testing.assert_array_equal(x, y)
    assert sorted(a.field_data) == sorted(b.field_data)
    assert sorted(a.cell_sets) == sorted(b.cell_sets)
    for k in a.cell_sets:
        for x, y in zip(a.cell_sets[k], b.cell_sets[k]):
            xa = np.asarray([] if x is None else x).ravel()
            ya = np.asarray([] if y is None else y).ravel()
            np.testing.assert_array_equal(xa, ya)


@pytest.mark.parametrize("filename", _ENTITY_FILES)
def test_cpp_matches_python_on_entities(filename):
    # Before $Entities landed in the C++ reader these files reached Python only
    # by way of the shim's blanket except -- so this is the gate that the two
    # readers now genuinely agree rather than one of them being unreachable.
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as py_read

    root = pathlib.Path(__file__).resolve().parents[2]
    path = str(root / filename)
    _assert_same_mesh(_core.gmsh_read(path), py_read(path))


@pytest.mark.parametrize("filename", _ENTITY_FILES)
@pytest.mark.parametrize("binary", [False, True])
def test_entities_survive_a_cpp_round_trip(filename, binary, tmp_path):
    from meshioplusplus import _core

    root = pathlib.Path(__file__).resolve().parents[2]
    mesh = _core.gmsh_read(str(root / filename))
    out = str(tmp_path / "rt.msh")
    _core.gmsh41_write(out, mesh, binary, mesh.cell_sets.get("gmsh:bounding_entities"))
    _assert_same_mesh(mesh, _core.gmsh_read(out))


def test_cpp_ascii_41_output_matches_the_python_writer(tmp_path):
    # Byte parity, not just equivalence: the entity records' spacing and line
    # discipline are easy to get subtly wrong and no read-back would notice.
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as py_read
    from meshioplusplus.gmsh.main import write as py_write

    root = pathlib.Path(__file__).resolve().parents[2]
    src = str(root / "tests/python/meshes/msh/insulated-4.1.msh")
    mesh = py_read(src)

    py_path = tmp_path / "py.msh"
    cpp_path = tmp_path / "cpp.msh"
    py_write(str(py_path), mesh, fmt_version="4.1", binary=False)
    _core.gmsh41_write(
        str(cpp_path),
        _core.gmsh_read(src),
        False,
        mesh.cell_sets.get("gmsh:bounding_entities"),
    )
    assert py_path.read_bytes() == cpp_path.read_bytes()


def test_writing_41_without_dim_tags_emits_no_entities(tmp_path):
    # No entity information to describe -> the legacy single-$Nodes-block
    # output, unchanged.
    from meshioplusplus import _core

    mesh = copy.deepcopy(helpers.tet_mesh)
    path = tmp_path / "plain.msh"
    _core.gmsh41_write(str(path), mesh, False)
    text = path.read_text()
    assert "$Entities" not in text
    # A single $Nodes block covering every point, as before.
    assert text.split("$Nodes\n")[1].splitlines()[0].split()[0] == "1"


def test_a_41_file_with_no_physical_groups_has_no_physical_cell_data():
    # example.msh tags no entity at all. Synthesizing an all-zero
    # gmsh:physical there would invent a group the file does not have.
    from meshioplusplus import _core

    root = pathlib.Path(__file__).resolve().parents[2]
    mesh = _core.gmsh_read(str(root / "example/example.msh"))
    assert "gmsh:physical" not in mesh.cell_data
    assert mesh.regions == []
    assert "gmsh:geometrical" in mesh.cell_data
