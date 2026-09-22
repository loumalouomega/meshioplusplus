import pathlib

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _core
from meshioplusplus.su2 import _su2 as py_su2

from . import helpers

test_set = [
    # helpers.empty_mesh,
    helpers.tri_mesh_2d,
    helpers.tet_mesh,
    helpers.hex_mesh,
]
this_dir = pathlib.Path(__file__).resolve().parent


@pytest.mark.parametrize("mesh", test_set)
def test(mesh, tmp_path):
    helpers.write_read(
        tmp_path, meshioplusplus.su2.write, meshioplusplus.su2.read, mesh, 1.0e-15
    )


@pytest.mark.parametrize(
    "filename, ref_num_cells, ref_num_points,ref_num_unique_tags,sum_tags",
    [("square.su2", 16, 9, 4, 20), ("mixgrid.su2", 30, 16, 6, 62)],
)
def test_structured(
    filename, ref_num_cells, ref_num_points, ref_num_unique_tags, sum_tags
):
    filename = this_dir / "meshes" / "su2" / filename

    mesh = meshioplusplus.read(filename)

    assert sum(len(block.data) for block in mesh.cells) == ref_num_cells
    assert len(mesh.points) == ref_num_points

    all_tags = np.concatenate([tags for tags in mesh.cell_data["su2:tag"]])

    assert sum(all_tags) == sum_tags

    all_unique_tags = np.unique(all_tags)

    assert len(all_unique_tags) == ref_num_unique_tags + 1


# --------------------------------------------------------------------------- #
# NZONE / IZONE: single-file multizone (roadmap 1.1)                          #
# --------------------------------------------------------------------------- #

MULTIZONE_TEXT = """\
NZONE= 2

IZONE= 1
NDIME= 2
NPOIN= 4
0.0 0.0
1.0 0.0
1.0 1.0
0.0 1.0
NELEM= 2
5 0 1 2
5 0 2 3
NMARK= 1
MARKER_TAG= wall
MARKER_ELEMS= 2
3 0 1
3 2 3

IZONE= 2
NDIME= 2
NPOIN= 4
2.0 0.0
3.0 0.0
3.0 1.0
2.0 1.0
NELEM= 2
5 0 1 2
5 0 2 3
NMARK= 2
MARKER_TAG= wall
MARKER_ELEMS= 1
3 0 1
MARKER_TAG= inlet
MARKER_ELEMS= 1
3 1 2
"""


@pytest.fixture(params=["core", "python"])
def engine(request):
    """Both engines behind the same read/write signature."""
    if request.param == "core":
        return _core.su2_read, _core.su2_write
    return py_su2.read, py_su2.write


def test_multizone_read_builds_zones_and_marker_regions(engine, tmp_path):
    read, _ = engine
    path = tmp_path / "multizone.su2"
    path.write_text(MULTIZONE_TEXT)
    mesh = read(str(path))

    assert len(mesh.points) == 8  # two disjoint 4-point zones, never welded
    assert sum(len(c.data) for c in mesh.cells) == 8  # 4 triangles + 4 lines

    assert "su2:zone" in mesh.cell_data
    assert "su2:tag" in mesh.cell_data

    regions = {r.name: r for r in mesh.regions}
    assert {"zone_0", "zone_1", "zone_0/wall", "zone_1/wall", "zone_1/inlet"} <= set(
        regions
    )
    assert len(regions["zone_0"].entries) == 4
    assert len(regions["zone_1/inlet"].entries) == 1


def test_multizone_round_trips(engine, tmp_path):
    read, write = engine
    path = tmp_path / "multizone.su2"
    path.write_text(MULTIZONE_TEXT)
    original = read(str(path))

    out = tmp_path / "out.su2"
    write(str(out), original)
    reread = read(str(out))

    np.testing.assert_allclose(original.points, reread.points)
    assert [c.type for c in original.cells] == [c.type for c in reread.cells]
    for c1, c2 in zip(original.cells, reread.cells):
        np.testing.assert_array_equal(c1.data, c2.data)
    for key in original.cell_data:
        for a, b in zip(original.cell_data[key], reread.cell_data[key]):
            np.testing.assert_array_equal(a, b)

    original_regions = {r.name: r.entries.tolist() for r in original.regions}
    reread_regions = {r.name: r.entries.tolist() for r in reread.regions}
    assert original_regions == reread_regions


def test_multizone_engines_agree(tmp_path):
    path = tmp_path / "multizone.su2"
    path.write_text(MULTIZONE_TEXT)
    core_mesh = _core.su2_read(str(path))
    py_mesh = py_su2.read(str(path))

    np.testing.assert_allclose(core_mesh.points, py_mesh.points)
    for c1, c2 in zip(core_mesh.cells, py_mesh.cells):
        assert c1.type == c2.type
        np.testing.assert_array_equal(c1.data, c2.data)
    for key in core_mesh.cell_data:
        for a, b in zip(core_mesh.cell_data[key], py_mesh.cell_data[key]):
            np.testing.assert_array_equal(a, b)
    assert sorted((r.name, r.entries.tolist()) for r in core_mesh.regions) == sorted(
        (r.name, r.entries.tolist()) for r in py_mesh.regions
    )


def test_single_zone_string_marker_names_round_trip_as_regions(engine, tmp_path):
    read, write = engine
    path = tmp_path / "single.su2"
    path.write_text(
        "NDIME= 2\nNPOIN= 4\n0 0\n1 0\n1 1\n0 1\nNELEM= 1\n5 0 1 2\nNMARK= 1\n"
        "MARKER_TAG= inlet\nMARKER_ELEMS= 1\n3 0 1\n"
    )
    mesh = read(str(path))
    region_names = {r.name for r in mesh.regions}
    assert "inlet" in region_names
    assert "zone_0" not in region_names  # single-zone: no zone region

    out = tmp_path / "out.su2"
    write(str(out), mesh)
    text = out.read_text()
    assert "MARKER_TAG= inlet" in text
    assert "NZONE" not in text  # a single zone never writes NZONE


def test_vtm_to_su2_carries_a_boundary_piece_name_as_a_marker(tmp_path):
    """.vtm names each of its pieces as a Cell region ("block_<i>", since the
    writer used here has no custom names to give); su2's own write side turns
    any such region over boundary-typed cells into a named marker -- this
    only exercises the mechanism (a region survives an unrelated format in
    between), not su2:zone, since .vtm has no zone concept of its own."""
    surface = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [("triangle", [[0, 1, 2]])],
    )
    vtm_path = tmp_path / "in.vtm"
    meshioplusplus.write(str(vtm_path), surface, file_format="vtm")
    back = meshioplusplus.read(str(vtm_path))
    assert any(r.kind == "cell" and r.name == "block_0" for r in back.regions)

    su2_path = tmp_path / "out.su2"
    meshioplusplus.write(str(su2_path), back, file_format="su2")
    text = su2_path.read_text()
    assert "MARKER_TAG= block_0" in text
