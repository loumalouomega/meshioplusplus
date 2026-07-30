import copy

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus._skin import _extract_skin_py

from . import helpers


def _cell_counts(mesh):
    counts = {}
    for block in mesh.cells:
        counts[block.type] = counts.get(block.type, 0) + len(block.data)
    return counts


def _hex27_mesh():
    # Unit cube corners, edge midpoints, face centers, body center — the VTK
    # triquadratic hexahedron layout.
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
    # Nodes 20-25 in vtkTriQuadraticHexahedron::Faces order (x-min, x-max,
    # y-min, y-max, bottom, top); corrected in v9.9.0 with _skin.py's table.
    face_centers = [
        (0, 4, 7, 3),
        (1, 2, 6, 5),
        (0, 1, 5, 4),
        (3, 7, 6, 2),
        (0, 1, 2, 3),
        (4, 5, 6, 7),
    ]
    pts = list(corners)
    for a, b in edges:
        pts.append((corners[a] + corners[b]) / 2)
    for f in face_centers:
        pts.append(corners[list(f)].mean(axis=0))
    pts.append(corners.mean(axis=0))
    return meshioplusplus.Mesh(np.array(pts), [("hexahedron27", [np.arange(27)])])


def _pyramid13_mesh(num_nodes=13):
    corners = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0.5, 0.5, 1]], dtype=float
    )
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4)]
    pts = list(corners)
    for a, b in edges:
        pts.append((corners[a] + corners[b]) / 2)
    if num_nodes == 14:
        pts.append(corners[:4].mean(axis=0))
    cell_type = "pyramid13" if num_nodes == 13 else "pyramid14"
    return meshioplusplus.Mesh(np.array(pts), [(cell_type, [np.arange(num_nodes)])])


def _pyramid_cube_mesh():
    # Unit cube split into 6 pyramids sharing the center apex (node 8): the
    # skin is the 6 cube faces and compaction must drop the center node.
    points = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
        [0.5, 0.5, 0.5],
    ]
    cells = [
        [0, 1, 2, 3, 8],
        [7, 6, 5, 4, 8],
        [0, 4, 5, 1, 8],
        [2, 6, 7, 3, 8],
        [0, 3, 7, 4, 8],
        [1, 5, 6, 2, 8],
    ]
    return meshioplusplus.Mesh(points, [("pyramid", cells)])


@pytest.mark.parametrize(
    "mesh, expected",
    [
        (helpers.tet_mesh, {"triangle": 6}),
        (helpers.hex_mesh, {"quad": 6}),
        (helpers.wedge_mesh, {"triangle": 2, "quad": 3}),
        (helpers.pyramid_mesh, {"triangle": 4, "quad": 1}),
        (helpers.tet10_mesh, {"triangle6": 4}),
        (helpers.hex20_mesh, {"quad8": 6}),
        (helpers.wedge15_mesh, {"triangle6": 2, "quad8": 3}),
        (_hex27_mesh(), {"quad9": 6}),
        (_pyramid13_mesh(13), {"triangle6": 4, "quad8": 1}),
        (_pyramid13_mesh(14), {"triangle6": 4, "quad9": 1}),
    ],
)
def test_counts(mesh, expected):
    skin = meshioplusplus.extract_skin(mesh)
    assert _cell_counts(skin) == expected


def test_two_hex_stack():
    # Two stacked cubes: 12 quad faces, 2 interior -> 10 boundary quads.
    points = []
    for layer in range(3):
        for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]:
            points.append([x, y, layer])
    mesh = meshioplusplus.Mesh(
        points, [("hexahedron", [[0, 1, 2, 3, 4, 5, 6, 7], [4, 5, 6, 7, 8, 9, 10, 11]])]
    )
    skin = meshioplusplus.extract_skin(mesh)
    assert _cell_counts(skin) == {"quad": 10}
    assert len(skin.points) == 12


def test_compaction_drops_interior_points():
    skin = meshioplusplus.extract_skin(_pyramid_cube_mesh())
    assert _cell_counts(skin) == {"quad": 6}
    # The center node is gone; compaction keeps ascending original order, so
    # the remaining points are exactly the 8 cube corners in order.
    assert len(skin.points) == 8
    assert np.array_equal(skin.points, np.asarray(_pyramid_cube_mesh().points)[:8])


def test_point_data_subset():
    mesh = _pyramid_cube_mesh()
    mesh.point_data["tag"] = np.arange(9.0) * 10.0
    skin = meshioplusplus.extract_skin(mesh)
    assert np.array_equal(skin.point_data["tag"], np.arange(8.0) * 10.0)


def test_linearize():
    skin = meshioplusplus.extract_skin(helpers.tet10_mesh, linearize=True)
    assert _cell_counts(skin) == {"triangle": 4}
    # The 6 mid-edge nodes compact away.
    assert len(skin.points) == 4


def test_surface_only_raises():
    with pytest.raises(ValueError):
        meshioplusplus.extract_skin(helpers.tri_mesh)
    with pytest.raises(ValueError):
        _extract_skin_py(helpers.tri_mesh)


def test_surface_blocks_ignored():
    # A mesh carrying both volume cells and a pre-existing surface block:
    # the skin comes from the volume cells only.
    mesh = copy.deepcopy(helpers.tet_mesh)
    mesh.cells.append(meshioplusplus.CellBlock("triangle", np.array([[0, 1, 2]])))
    skin = meshioplusplus.extract_skin(mesh)
    assert _cell_counts(skin) == {"triangle": 6}


def test_polyhedron_block_skipped():
    # A polyhedron block next to a tetra block is skipped with a warning; the
    # tetra skin is still extracted (via the numpy fallback, since ragged
    # blocks cannot cross the pybind boundary).
    mesh = copy.deepcopy(helpers.polyhedron_mesh)
    mesh.cells.append(meshioplusplus.CellBlock("tetra", np.array([[0, 1, 3, 4]])))
    skin = meshioplusplus.extract_skin(mesh)
    assert _cell_counts(skin) == {"triangle": 4}


def test_polyhedron_only_raises():
    # A mesh whose only volume blocks are polyhedra has nothing to extract.
    with pytest.raises(ValueError):
        meshioplusplus.extract_skin(copy.deepcopy(helpers.polyhedron_mesh))


@pytest.mark.parametrize(
    "mesh",
    [
        helpers.tet_mesh,
        helpers.hex_mesh,
        helpers.wedge_mesh,
        helpers.pyramid_mesh,
        helpers.tet10_mesh,
        helpers.hex20_mesh,
        helpers.wedge15_mesh,
        _hex27_mesh(),
        _pyramid13_mesh(13),
        _pyramid13_mesh(14),
        _pyramid_cube_mesh(),
    ],
)
@pytest.mark.parametrize("linearize", [False, True])
def test_cpp_matches_python(mesh, linearize):
    # The C++ core and the numpy reference produce identical results —
    # same block order, same connectivity, same compacted point order.
    cpp = meshioplusplus._core.extract_skin(mesh, linearize)
    py = _extract_skin_py(mesh, linearize=linearize)
    assert np.array_equal(cpp.points, py.points)
    assert len(cpp.cells) == len(py.cells)
    for a, b in zip(cpp.cells, py.cells):
        assert a.type == b.type
        assert np.array_equal(a.data, b.data)
