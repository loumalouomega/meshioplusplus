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
        helpers.pyramid13_mesh,
        helpers.wedge15_mesh,
        helpers.pyramid14_mesh,
        helpers.wedge18_mesh,
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
    """A `.msh` written without a format and without gmsh tags is Fluent's
    (the first `.msh` candidate): it reads back as the same planar triangles,
    rebuilt from their faces, so compare them as point sets."""
    mesh = helpers.tri_mesh
    # With additional, insignificant suffix too:
    for path in (tmp_path / "test.msh", tmp_path / "test.0.msh"):
        meshioplusplus.write_points_cells(path, mesh.points, mesh.cells)
        out = meshioplusplus.read(path)
        dim = out.points.shape[1]
        assert (mesh.points[:, dim:] == 0).all()

        def triangles(m, pts):
            return sorted(
                tuple(sorted(map(tuple, pts[row].tolist())))
                for block in m.cells
                if block.type == "triangle"
                for row in block.data
            )

        assert triangles(out, out.points) == triangles(mesh, mesh.points[:, :dim])


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


def test_read_metadata_reports_both_steps_of_a_transient_file(tmp_path):
    """roadmap §1 tier B1: read_gmsh_metadata (4.1) reports the sorted union
    of every $NodeData/$ElementData section's time value, and ``time_step``
    on a real read picks exactly the one matching that resolved time --
    replacing the pre-v11.3.0 "first section per name wins" rule. A second
    step is appended by hand (ASCII text, no external tool)."""
    from meshioplusplus import _core

    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["u"] = np.array([10.0, 20.0, 30.0, 40.0])
    path = tmp_path / "transient.msh"
    _core.gmsh41_write(str(path), mesh, False)

    n = len(mesh.points)
    lines = [
        "$NodeData",
        "1",
        '"u"',
        "1",
        "2.5",
        "3",
        "1",
        "1",
        str(n),
    ]
    lines += [f"{i + 1} {11.0 + i * 10.0}" for i in range(n)]
    lines.append("$EndNodeData")
    with open(path, "a") as f:
        f.write("\n".join(lines) + "\n")

    meta = meshioplusplus.read_metadata(path, "gmsh")
    assert meta["fell_back_to_full_read"] is False
    assert meta["time_values"] == [0.0, 2.5]
    assert meta["point_data_names"].count("u") == 1

    mesh0 = meshioplusplus.gmsh.read(path, time_step=0)
    assert mesh0.point_data["u"][0] == 10.0
    mesh1 = meshioplusplus.gmsh.read(path, time_step=1)
    assert mesh1.point_data["u"][0] == 11.0
    mesh_last = meshioplusplus.gmsh.read(path, time_step=-1)
    assert mesh_last.point_data["u"][0] == 11.0

    with pytest.raises(meshioplusplus.ReadError):
        meshioplusplus.gmsh.read(path, time_step=5)


def test_malformed_time_value_raises_read_error_not_value_error(tmp_path):
    """A $NodeData section whose time value isn't a valid number used to
    reach std::stod directly, propagating as a bare ValueError instead of a
    proper ReadError -- part of the roadmap's locale item (std::stod's
    exceptions are the reason A4 gave these two call sites an explicit
    check rather than a blind swap to the non-throwing parse_double)."""
    from meshioplusplus import _core

    mesh = copy.deepcopy(helpers.tri_mesh)
    mesh.point_data["u"] = np.array([10.0, 20.0, 30.0, 40.0])
    path = tmp_path / "bad_time.msh"
    _core.gmsh41_write(str(path), mesh, False)

    n = len(mesh.points)
    lines = [
        "$NodeData",
        "1",
        '"u"',
        "1",
        "not-a-number",
        "3",
        "1",
        "1",
        str(n),
    ]
    lines += [f"{i + 1} {11.0 + i * 10.0}" for i in range(n)]
    lines.append("$EndNodeData")
    with open(path, "a") as f:
        f.write("\n".join(lines) + "\n")

    # Straight through the C++ core, not the meshioplusplus.gmsh.read() shim:
    # the shim's `except Exception: if time_step: raise` only re-raises for a
    # truthy time_step, and 0 (the default/first step) is falsy, so it would
    # otherwise swallow this and silently fall back to the Python reader --
    # a separate, already-tracked gap (roadmap §1's other reader-fallback
    # item), not what this test is about.
    with pytest.raises(meshioplusplus.ReadError, match="time value"):
        _core.gmsh_read(str(path), time_step=0)


def test_untagged_region_gets_an_allocated_tag_on_write(tmp_path):
    """Roadmap §1 tier B3 (v11.5.0): a Cell region with no gmsh tag of its
    own (as Abaqus/MED/MDPA produce) gets a freshly allocated one instead of
    being dropped from $PhysicalNames/gmsh:physical. C++-core only -- the
    pure-Python fallback writer keeps its existing field_data-only behaviour.
    """
    mesh = meshioplusplus.Mesh(
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
        [("hexahedron", [[0, 1, 2, 3, 4, 5, 6, 7]])],
        regions=[meshioplusplus.Region("body", "cell", [0])],
    )
    path = tmp_path / "body.msh"
    meshioplusplus.write(path, mesh, file_format="gmsh22")
    back = meshioplusplus.read(path)
    cell_regions = {r.name: r for r in back.regions if r.kind == "cell"}
    assert "body" in cell_regions
    assert cell_regions["body"].tag >= 0, "an allocated tag must not round-trip as -1"
    assert list(cell_regions["body"].entries) == [0]


# --- roadmap §1 "gmsh pyramid14 (and wedge18) node ordering is not
# permuted" ---
#
# A write->read round trip cannot catch a wrong-but-consistently-inverse
# permutation pair, so this checks against gmsh's OWN edge/face tables
# independently: `_meshio_to_gmsh_order` tells us which meshio corner sits
# at each gmsh slot; gmsh's corners 0..k_corner-1 are unchanged from
# meshio's (confirmed against gmsh's own src/geo/MPyramid.h /
# src/geo/MPrism.h edge/face tables), so re-deriving each mid-edge/
# face-centre point from THOSE corners and comparing to what actually
# landed there is a check independent of meshio's own edge/face
# convention (which is what the fixture meshes in helpers.py were built
# from).


def _assert_gmsh_order_matches_gmsh_geometry(
    cell_type, points, num_corners, edges, faces
):
    from meshioplusplus.gmsh.common import _meshio_to_gmsh_order

    idx = np.arange(len(points)).reshape(1, -1)
    gmsh_conn = _meshio_to_gmsh_order(cell_type, idx)[0]
    gmsh_pts = points[gmsh_conn]

    for slot, (a, b) in enumerate(edges, start=num_corners):
        expected = (gmsh_pts[a] + gmsh_pts[b]) / 2.0
        np.testing.assert_allclose(
            gmsh_pts[slot],
            expected,
            atol=1e-12,
            err_msg=f"{cell_type}: gmsh mid-edge slot {slot} (edge {a},{b})",
        )
    edge_end = num_corners + len(edges)
    for slot, corners in enumerate(faces, start=edge_end):
        expected = gmsh_pts[list(corners)].mean(axis=0)
        np.testing.assert_allclose(
            gmsh_pts[slot],
            expected,
            atol=1e-12,
            err_msg=f"{cell_type}: gmsh face-centre slot {slot} (face {corners})",
        )


def test_pyramid14_gmsh_order_matches_gmsh_own_edge_and_face_tables():
    _assert_gmsh_order_matches_gmsh_geometry(
        "pyramid14",
        helpers.pyramid14_mesh.points,
        num_corners=5,
        edges=[(0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (2, 4), (3, 4)],
        faces=[(0, 3, 2, 1)],
    )


def test_wedge18_gmsh_order_matches_gmsh_own_edge_and_face_tables():
    _assert_gmsh_order_matches_gmsh_geometry(
        "wedge18",
        helpers.wedge18_mesh.points,
        num_corners=6,
        edges=[
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 2),
            (1, 4),
            (2, 5),
            (3, 4),
            (3, 5),
            (4, 5),
        ],
        faces=[(0, 1, 4, 3), (0, 3, 5, 2), (1, 2, 5, 4)],
    )


def test_pyramid13_gmsh_order_matches_gmsh_own_edge_table():
    # The already-shipped pyramid13 table, as an in-repo positive control
    # for the helper above.
    _assert_gmsh_order_matches_gmsh_geometry(
        "pyramid13",
        helpers.pyramid13_mesh.points,
        num_corners=5,
        edges=[(0, 1), (0, 3), (0, 4), (1, 2), (1, 4), (2, 3), (2, 4), (3, 4)],
        faces=[],
    )


def test_wedge15_gmsh_order_matches_gmsh_own_edge_table():
    _assert_gmsh_order_matches_gmsh_geometry(
        "wedge15",
        helpers.wedge15_mesh.points,
        num_corners=6,
        edges=[
            (0, 1),
            (0, 2),
            (0, 3),
            (1, 2),
            (1, 4),
            (2, 5),
            (3, 4),
            (3, 5),
            (4, 5),
        ],
        faces=[],
    )


@pytest.mark.parametrize("engine", ["core", "python"])
def test_indented_file_reads(engine):
    # FEconv's samples indent every line (tools/gen_feconv_quirk_fixtures.py);
    # the third element line is past the declared count and ignored.
    from meshioplusplus import _core
    from meshioplusplus.gmsh.main import read as py_read

    path = pathlib.Path(__file__).resolve().parent / "meshes" / "gmsh" / "indented.msh"
    mesh = _core.gmsh_read(str(path)) if engine == "core" else py_read(path)
    assert len(mesh.points) == 6
    assert [c.type for c in mesh.cells] == ["quad"]
    np.testing.assert_array_equal(mesh.cells[0].data, [[0, 1, 2, 3], [1, 4, 5, 2]])
