"""Uniform refinement (``meshioplusplus.refine``)."""

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import refine


def _two_triangles():
    """Two triangles sharing edge (0, 2): 4 points, 5 unique edges."""
    return meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]]))],
    )


def _unit_cube():
    return meshioplusplus.Mesh(
        np.array(
            [
                [0.0, 0, 0],
                [1, 0, 0],
                [1, 1, 0],
                [0, 1, 0],
                [0, 0, 1],
                [1, 0, 1],
                [1, 1, 1],
                [0, 1, 1],
            ]
        ),
        [("hexahedron", np.array([[0, 1, 2, 3, 4, 5, 6, 7]]))],
    )


def _unit_tet():
    return meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        [("tetra", np.array([[0, 1, 2, 3]]))],
    )


def _unit_quad():
    return meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [2, 0, 0], [2, 2, 0], [0, 2, 0]]),
        [("quad", np.array([[0, 1, 2, 3]]))],
    )


def _signed_tet_volume(pts, row):
    a, b, c, d = (pts[i] for i in row)
    return np.dot(np.cross(b - a, c - a), d - a) / 6.0


# --- 2D --------------------------------------------------------------------- #


def test_triangle_splits_into_four_with_shared_mid_edge_nodes():
    out = refine(_two_triangles())
    assert len(out.cells) == 1
    assert out.cells[0].type == "triangle"
    assert len(out.cells[0].data) == 8
    # 4 original + 5 unique edges -- NOT 4 + 2*3, which is what an undeduped
    # per-cell midpoint would give.
    assert len(out.points) == 9


def test_triangle_mid_nodes_sit_at_edge_midpoints():
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [2, 0, 0], [0, 2, 0]]),
        [("triangle", np.array([[0, 1, 2]]))],
    )
    out = refine(mesh)
    assert len(out.points) == 6
    # New nodes are appended in edge order (0,1), (1,2), (2,0).
    assert np.allclose(out.points[3], [1, 0, 0])
    assert np.allclose(out.points[4], [1, 1, 0])
    assert np.allclose(out.points[5], [0, 1, 0])


def test_linear_point_data_is_interpolated_exactly():
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [2, 0, 0], [0, 2, 0]]),
        [("triangle", np.array([[0, 1, 2]]))],
    )
    # f(x, y) = 3x + 5y + 1
    mesh.point_data["f"] = np.array([1.0, 7.0, 11.0])
    out = refine(mesh)
    expected = 3 * out.points[:, 0] + 5 * out.points[:, 1] + 1
    assert np.allclose(out.point_data["f"], expected)


def test_quad_splits_into_four_with_centre_at_corner_mean():
    out = refine(_unit_quad())
    assert out.cells[0].type == "quad"
    assert len(out.cells[0].data) == 4
    # 4 corners + 4 edge mids + 1 face centre.
    assert len(out.points) == 9
    assert np.allclose(out.points[8], [1, 1, 0])


# --- 3D --------------------------------------------------------------------- #


def test_tetra_splits_into_eight_conserving_volume_and_orientation():
    mesh = _unit_tet()
    before = meshioplusplus.compute_stats(mesh)
    out = refine(mesh)
    after = meshioplusplus.compute_stats(out)

    assert out.cells[0].type == "tetra"
    assert len(out.cells[0].data) == 8
    assert after["signed_volume"] == pytest.approx(before["signed_volume"])
    assert after["num_inverted"] == 0

    pts = np.asarray(out.points)
    for row in out.cells[0].data:
        assert _signed_tet_volume(pts, row) > 0


def test_hex_splits_into_eight_with_27_nodes():
    mesh = _unit_cube()
    before = meshioplusplus.compute_stats(mesh)
    out = refine(mesh)
    after = meshioplusplus.compute_stats(out)

    assert out.cells[0].type == "hexahedron"
    assert len(out.cells[0].data) == 8
    # 8 corners + 12 edge mids + 6 face centres + 1 body.
    assert len(out.points) == 27
    assert np.allclose(out.points[26], [0.5, 0.5, 0.5])
    assert after["signed_volume"] == pytest.approx(before["signed_volume"])
    assert after["num_inverted"] == 0


def test_wedge_splits_into_eight_with_18_nodes():
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]]),
        [("wedge", np.array([[0, 1, 2, 3, 4, 5]]))],
    )
    before = meshioplusplus.compute_stats(mesh)
    out = refine(mesh)
    after = meshioplusplus.compute_stats(out)

    assert out.cells[0].type == "wedge"
    assert len(out.cells[0].data) == 8
    # 6 corners + 9 edge mids + 3 quad-face centres, no body node.
    assert len(out.points) == 18
    assert after["signed_volume"] == pytest.approx(before["signed_volume"])
    assert after["num_inverted"] == 0


def test_adjacent_hexes_refine_conformingly():
    """Face-centre nodes must be shared, or every interior face becomes boundary."""
    mesh = meshioplusplus.Mesh(
        np.array(
            [
                [0.0, 0, 0],
                [1, 0, 0],
                [1, 1, 0],
                [0, 1, 0],
                [0, 0, 1],
                [1, 0, 1],
                [1, 1, 1],
                [0, 1, 1],
                [0, 0, 2],
                [1, 0, 2],
                [1, 1, 2],
                [0, 1, 2],
            ]
        ),
        [
            (
                "hexahedron",
                np.array([[0, 1, 2, 3, 4, 5, 6, 7], [4, 5, 6, 7, 8, 9, 10, 11]]),
            )
        ],
    )
    out = refine(mesh)
    # A conforming 2x2x4 lattice of hexahedra has 3x3x5 nodes.
    assert len(out.points) == 45

    before = sum(len(cb.data) for cb in meshioplusplus.extract_surface(mesh).cells)
    after = sum(len(cb.data) for cb in meshioplusplus.extract_surface(out).cells)
    # Every boundary facet splits into 4; no interior facet may show up.
    assert after == before * 4


# --- levels, data, maps ------------------------------------------------------ #


def test_levels_two_equals_refining_twice():
    mesh = _two_triangles()
    direct = refine(mesh, levels=2)
    twice = refine(refine(mesh))
    assert np.array_equal(direct.cells[0].data, twice.cells[0].data)
    assert np.allclose(direct.points, twice.points)


def test_zero_levels_is_unchanged():
    mesh = _two_triangles()
    out = refine(mesh, levels=0)
    assert np.array_equal(out.cells[0].data, mesh.cells[0].data)
    assert np.allclose(out.points, mesh.points)


def test_cell_data_is_replicated_to_children():
    mesh = _two_triangles()
    mesh.cell_data["tag"] = [np.array([10, 20])]
    out = refine(mesh)
    assert list(out.cell_data["tag"][0]) == [10, 10, 10, 10, 20, 20, 20, 20]


def test_record_parent_ids_names_the_original_ancestor():
    out = refine(_two_triangles(), levels=2, record_parent_ids=True)
    parents = np.asarray(out.cell_data["refine:parent_cell"][0]).reshape(-1)
    assert len(parents) == 32  # 2 parents x 16 grandchildren
    assert set(parents[:16]) == {0}
    assert set(parents[16:]) == {1}


def test_cell_sets_expand_to_children():
    mesh = _two_triangles()
    mesh.cell_sets = {"left": [np.array([0])]}
    out = refine(mesh)
    assert list(out.cell_sets["left"][0]) == [0, 1, 2, 3]


def test_mixed_type_mesh_refines_each_block_by_its_own_template():
    mesh = meshioplusplus.Mesh(
        np.array([[0.0, 0, 0], [1, 0, 0], [2, 0, 0], [2, 1, 0], [1, 1, 0], [0, 1, 0]]),
        [
            ("triangle", np.array([[0, 1, 4], [0, 4, 5]])),
            ("quad", np.array([[1, 2, 3, 4]])),
        ],
    )
    out = refine(mesh)
    assert [cb.type for cb in out.cells] == ["triangle", "quad"]
    assert len(out.cells[0].data) == 8
    assert len(out.cells[1].data) == 4


# --- rejections, parity, I/O ------------------------------------------------- #


@pytest.mark.parametrize(
    "cell_type, conn",
    [
        ("tetra10", np.arange(10).reshape(1, 10)),
        ("triangle6", np.arange(6).reshape(1, 6)),
        ("pyramid", np.arange(5).reshape(1, 5)),
    ],
)
def test_unsupported_cell_types_raise(cell_type, conn):
    mesh = meshioplusplus.Mesh(np.zeros((16, 3)), [(cell_type, conn)])
    with pytest.raises(ValueError):
        refine(mesh)


def test_ragged_block_raises():
    mesh = meshioplusplus.Mesh(
        np.zeros((5, 3)), [("polygon", [[0, 1, 2, 3], [3, 2, 4]])]
    )
    with pytest.raises(ValueError):
        refine(mesh)


@pytest.mark.parametrize("factory", [_two_triangles, _unit_quad, _unit_tet, _unit_cube])
def test_cpp_matches_python(factory):
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._refine import _refine_py

    mesh = factory()
    got = core.refine(mesh, 2, False)["mesh"]
    ref, _, _ = _refine_py(mesh, 2, False)

    assert [cb.type for cb in got.cells] == [cb.type for cb in ref.cells]
    for a, b in zip(got.cells, ref.cells):
        assert np.array_equal(a.data, b.data)
    assert np.allclose(got.points, ref.points)


def test_roundtrip_write_read(tmp_path):
    out = refine(_unit_cube())
    path = tmp_path / "refined.vtu"
    meshioplusplus.write(path, out)
    back = meshioplusplus.read(path)
    assert len(back.points) == len(out.points)
    assert np.array_equal(back.cells[0].data, out.cells[0].data)


# --------------------------------------------------------------------------- #
# selective / adaptive refinement                                              #
# --------------------------------------------------------------------------- #
def _tri_grid(n):
    pts = [[float(i), float(j), 0.0] for j in range(n + 1) for i in range(n + 1)]
    cells = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            cells.append([a, a + 1, a + n + 2])
            cells.append([a, a + n + 2, a + n + 1])
    return meshioplusplus.Mesh(np.array(pts), [("triangle", np.array(cells))])


def _quad_grid(n):
    pts = [[float(i), float(j), 0.0] for j in range(n + 1) for i in range(n + 1)]
    cells = [
        [
            j * (n + 1) + i,
            j * (n + 1) + i + 1,
            j * (n + 1) + i + n + 2,
            j * (n + 1) + i + n + 1,
        ]
        for j in range(n)
        for i in range(n)
    ]
    return meshioplusplus.Mesh(np.array(pts), [("quad", np.array(cells))])


def _hex_grid(n):
    s = n + 1
    pts = [
        [float(i), float(j), float(k)]
        for k in range(s)
        for j in range(s)
        for i in range(s)
    ]

    def idx(i, j, k):
        return (k * s + j) * s + i

    cells = [
        [
            idx(i, j, k),
            idx(i + 1, j, k),
            idx(i + 1, j + 1, k),
            idx(i, j + 1, k),
            idx(i, j, k + 1),
            idx(i + 1, j, k + 1),
            idx(i + 1, j + 1, k + 1),
            idx(i, j + 1, k + 1),
        ]
        for k in range(n)
        for j in range(n)
        for i in range(n)
    ]
    return meshioplusplus.Mesh(np.array(pts), [("hexahedron", np.array(cells))])


def _count_hanging_nodes(mesh):
    """Output nodes sitting exactly on another cell's edge midpoint.

    The oracle that actually finds a hanging node: facet counting cannot, since
    a hanging node leaves one big facet on one side and two small ones on the
    other, all with a count of 1. New nodes are order-independent corner means,
    so the comparison is exact rather than approximate.
    """
    from meshioplusplus._refine_templates import EDGES

    points = np.asarray(mesh.points)
    at = {tuple(p): i for i, p in enumerate(points)}
    hanging = 0
    for cb in mesh.cells:
        data = np.asarray(cb.data)
        for a, b in EDGES[cb.type]:
            mids = (points[data[:, a]] + points[data[:, b]]) * 0.5
            hanging += sum(1 for m in mids if tuple(m) in at)
    return hanging


def _assert_conforming(mesh):
    assert _count_hanging_nodes(mesh) == 0, "the refined mesh has hanging nodes"


def test_the_hanging_node_oracle_actually_fires():
    # An oracle that cannot fail proves nothing. The left quadrilateral is split
    # across its mid height and its neighbour is not, so node 7 hangs.
    bad = meshioplusplus.Mesh(
        np.array(
            [
                [0.0, 0, 0],
                [1, 0, 0],
                [2, 0, 0],
                [0, 1, 0],
                [1, 1, 0],
                [2, 1, 0],
                [0, 0.5, 0],
                [1, 0.5, 0],
            ]
        ),
        [("quad", np.array([[0, 1, 7, 6], [6, 7, 4, 3], [1, 2, 5, 4]]))],
    )
    assert _count_hanging_nodes(bad) == 1


def test_no_selector_is_the_uniform_refinement():
    mesh = _tri_grid(2)
    assert np.array_equal(
        np.asarray(refine(mesh).cells[0].data),
        np.asarray(refine(mesh, cells=None).cells[0].data),
    )


def test_one_triangle_greens_its_neighbours_conformingly():
    mesh = _tri_grid(3)
    before = len(mesh.cells[0].data)
    out = refine(mesh, cells=[8])
    _assert_conforming(out)
    assert before < len(out.cells[0].data) < before + 12
    assert out.cells[0].type == "triangle"


def test_one_quad_refines_one_row_not_the_whole_grid():
    # The anisotropic rule: a single split edge promotes to its opposite pair,
    # so the bisection travels along one row instead of the whole grid. Full
    # propagation would give 4 * 64 = 256 cells.
    n = 8
    mesh = _quad_grid(n)
    out = refine(mesh, cells=[27])
    _assert_conforming(out)
    assert out.cells[0].type == "quad", "green quads stay quads"
    assert (
        len(mesh.cells[0].data)
        < len(out.cells[0].data)
        < len(mesh.cells[0].data) + 6 * n
    )


def test_one_hex_refines_sheets_not_the_whole_block():
    n = 4
    mesh = _hex_grid(n)
    out = refine(mesh, cells=[21])
    _assert_conforming(out)
    assert out.cells[0].type == "hexahedron"
    assert len(out.cells[0].data) < len(mesh.cells[0].data) + 4 * n * n


@pytest.mark.parametrize("factory,n", [(_tri_grid, 3), (_quad_grid, 3), (_hex_grid, 2)])
def test_propagate_reproduces_uniform_refinement(factory, n):
    # The oracle, and the honest statement of what that mode costs: on a
    # connected mesh propagation reaches every cell.
    mesh = factory(n)
    uniform = refine(mesh)
    propagated = refine(mesh, cells=[0], closure="propagate")
    assert np.array_equal(
        np.asarray(uniform.cells[0].data), np.asarray(propagated.cells[0].data)
    )
    assert np.array_equal(np.asarray(uniform.points), np.asarray(propagated.points))


def test_two_successive_passes_stay_conforming():
    mesh = _tri_grid(4)
    first = refine(mesh, cells=[10])
    _assert_conforming(first)
    _assert_conforming(refine(first, cells=[0]))


def test_selectors_agree_with_each_other():
    mesh = _tri_grid(3)
    by_index = refine(mesh, cells=[4, 9])

    region_mesh = _tri_grid(3)
    region_mesh.cell_sets["hot"] = [np.array([4, 9])]
    by_region = refine(region_mesh, region="hot")
    assert np.array_equal(
        np.asarray(by_index.cells[0].data), np.asarray(by_region.cells[0].data)
    )

    pred_mesh = _tri_grid(3)
    score = np.full(len(pred_mesh.cells[0].data), 0.9)
    score[[4, 9]] = 0.1
    pred_mesh.cell_data["q"] = [score]
    by_pred = refine(pred_mesh, where="q < 0.3")
    assert np.array_equal(
        np.asarray(by_index.cells[0].data), np.asarray(by_pred.cells[0].data)
    )


def test_a_point_region_selects_every_cell_touching_it():
    mesh = _tri_grid(3)
    mesh.point_sets["corner"] = np.array([0])
    out = refine(mesh, region="corner")
    _assert_conforming(out)
    assert len(out.cells[0].data) > len(mesh.cells[0].data)


def test_non_finite_predicate_values_never_match():
    # compute_quality reports NaN where a metric does not apply, so a predicate
    # over `quality:*` on a mixed mesh is the headline use case; rejecting the
    # array outright would break it.
    mesh = _tri_grid(2)
    mesh.cell_data["q"] = [np.full(len(mesh.cells[0].data), np.nan)]
    out = refine(mesh, where="q < 1.0")
    assert len(out.cells[0].data) == len(mesh.cells[0].data)


def test_levels_are_opt_in_and_green_children_inherit():
    mesh = _tri_grid(3)
    assert "refine:level" not in refine(mesh, cells=[8]).cell_data

    out = refine(mesh, cells=[8], record_levels=True)
    levels = np.asarray(out.cell_data["refine:level"][0]).reshape(-1)
    first = np.asarray(_first_child_map(mesh, out, 8))
    assert set(levels[first[0] : first[1]]) == {1}, "the selected cell's children"
    assert 0 in set(levels), "a transitional split is a closure, not a refinement"


def _first_child_map(mesh, out, parent):
    """`[first, next_first)` for one parent, recomputed from a fresh call."""
    from meshioplusplus import _core

    res = _core.refine(
        mesh,
        1,
        False,
        np.array([parent], dtype=np.int64),
        "",
        "",
        "<",
        0.0,
        "redgreen",
        True,
    )
    first = np.asarray(res["cell_maps"][0])
    total = len(res["mesh"].cells[0].data)
    return (first[parent], first[parent + 1] if parent + 1 < len(first) else total)


def test_an_existing_level_array_accumulates_rather_than_replicating():
    mesh = _tri_grid(3)
    first = refine(mesh, cells=[4], record_levels=True)
    # The flag only controls CREATING the array; once present it is maintained.
    again = refine(first, cells=[0])
    levels = np.asarray(again.cell_data["refine:level"][0]).reshape(-1)
    assert levels.max() >= 1


# --------------------------------------------------------------------------- #
# refine:cell_id / refine:parent_id: the persistent hierarchy                  #
# --------------------------------------------------------------------------- #
def test_hierarchy_is_opt_in():
    mesh = _tri_grid(3)
    out = refine(mesh, cells=[8])
    assert "refine:cell_id" not in out.cell_data
    assert "refine:parent_id" not in out.cell_data


def test_an_unsplit_cell_keeps_its_id_and_is_its_own_parent():
    mesh = _tri_grid(3)
    out = refine(mesh, cells=[8], record_hierarchy=True)
    ids = np.asarray(out.cell_data["refine:cell_id"][0]).reshape(-1)
    parents = np.asarray(out.cell_data["refine:parent_id"][0]).reshape(-1)
    first, last = _first_child_map(mesh, out, 8)
    # A cell far from the selection is untouched; it must keep its own
    # implicit id (its own global block-major index).
    untouched = 0
    for c in range(len(ids)):
        if first <= c < last:
            continue
        if ids[c] == parents[c]:
            untouched += 1
    assert untouched > 0


def test_every_parent_id_resolves_in_the_coarse_mesh():
    mesh = _tri_grid(3)
    out = refine(mesh, cells=[8], record_hierarchy=True)
    coarse_n = len(mesh.cells[0].data)
    parents = np.asarray(out.cell_data["refine:parent_id"][0]).reshape(-1)
    assert np.all(parents >= 0) and np.all(parents < coarse_n)


def test_fresh_ids_exceed_every_input_id():
    mesh = _tri_grid(3)
    first = refine(mesh, cells=[8], record_hierarchy=True)
    max_first_id = np.asarray(first.cell_data["refine:cell_id"][0]).max()
    second = refine(first, cells=[0], record_hierarchy=True)
    ids = np.asarray(second.cell_data["refine:cell_id"][0]).reshape(-1)
    parents = np.asarray(second.cell_data["refine:parent_id"][0]).reshape(-1)
    fresh = ids[ids != parents]
    assert len(fresh) > 0 and np.all(fresh > max_first_id)


def test_an_existing_hierarchy_accumulates_rather_than_replicating():
    mesh = _tri_grid(3)
    first = refine(mesh, cells=[4], record_hierarchy=True)
    # The flag only controls CREATING the arrays; once present they are
    # maintained even without the flag.
    again = refine(first, cells=[0])
    assert "refine:cell_id" in again.cell_data
    ids = np.asarray(again.cell_data["refine:cell_id"][0]).reshape(-1)
    assert len(set(ids.tolist())) == len(ids), "ids stay unique across calls"


def test_two_levels_names_the_original_input_like_parent_cell():
    mesh = _tri_grid(3)
    out = refine(mesh, levels=2, record_parent_ids=True, record_hierarchy=True)
    parent_cell = np.asarray(out.cell_data["refine:parent_cell"][0]).reshape(-1)
    parent_id = np.asarray(out.cell_data["refine:parent_id"][0]).reshape(-1)
    # For a mesh with no PRIOR hierarchy the implicit id of an input cell IS
    # its row index, so the two must induce the same partition even across
    # several internal levels.
    assert np.array_equal(parent_cell, parent_id)


def test_record_hierarchy_also_attaches_the_entity_keys():
    # redgreen leaves no hanging nodes, so refine:entity would normally never
    # be attached -- but record_hierarchy needs it as the multigrid
    # prolongation stencil.
    mesh = _tri_grid(3)
    out = refine(mesh, cells=[8], closure="redgreen", record_hierarchy=True)
    assert "refine:entity" in out.point_data
    keys = np.asarray(out.point_data["refine:entity"])
    pts = np.asarray(out.points)
    new_nodes = 0
    for p in range(len(pts)):
        key = keys[p]
        if key[3] < 0:
            continue
        corners = key[2:] if key[0] < 0 else key
        want = pts[list(corners)].mean(axis=0)
        assert np.allclose(pts[p], want)
        new_nodes += 1
    assert new_nodes > 0


def test_split_by_parent_id_groups_siblings():
    mesh = _tri_grid(3)
    out = refine(mesh, cells=[4, 8], record_hierarchy=True)
    parents = np.asarray(out.cell_data["refine:parent_id"][0]).reshape(-1)
    distinct_parents = set(parents.tolist())
    pieces = meshioplusplus.split(out, by="tag", tag="refine:parent_id")
    assert len(pieces) == len(distinct_parents)


def test_hierarchy_classifies_red_green_untouched_from_the_two_meshes():
    mesh = _tri_grid(3)
    out = refine(mesh, cells=[8], record_hierarchy=True, record_levels=True)
    ids = np.asarray(out.cell_data["refine:cell_id"][0]).reshape(-1)
    parents = np.asarray(out.cell_data["refine:parent_id"][0]).reshape(-1)
    levels = np.asarray(out.cell_data["refine:level"][0]).reshape(-1)
    untouched = green = red = 0
    for c in range(len(ids)):
        if ids[c] == parents[c]:
            assert levels[c] == 0
            untouched += 1
        elif levels[c] == 1:
            red += 1
        else:
            assert levels[c] == 0
            green += 1
    assert untouched > 0 and green > 0 and red > 0


def test_a_duplicated_id_is_dropped_rather_than_trusted():
    """A duplicated id means the mesh was merged or the array was replicated by
    another operation. Mirrors test_a_stale_entity_array_is_ignored_rather_than_trusted:
    the warning is asserted on the pure-numpy reader directly (refine() itself
    may take the C++ path, whose equivalent warning goes to meshio++'s own log,
    not a Python warning), and the drop behaviour is asserted through refine()."""
    mesh = _tri_grid(2)
    n = len(mesh.cells[0].data)
    ids = np.arange(n, dtype=np.int64)
    ids[0] = ids[1]  # duplicate
    mesh.cell_data["refine:cell_id"] = [ids.reshape(-1, 1)]
    mesh.cell_data["refine:parent_id"] = [ids.reshape(-1, 1)]

    with pytest.warns(UserWarning, match="not both unique"):
        from meshioplusplus._refine import _read_hierarchy

        blocks = [(cb.type, np.asarray(cb.data)) for cb in mesh.cells]
        got_ids, _ = _read_hierarchy(mesh, blocks)
        assert got_ids is None

    out = refine(mesh, cells=[0])
    assert "refine:cell_id" not in out.cell_data
    assert "refine:parent_id" not in out.cell_data


def test_cell_data_is_gathered_not_blindly_repeated():
    mesh = _tri_grid(2)
    n = len(mesh.cells[0].data)
    mesh.cell_data["tag"] = [np.arange(n)]
    out = refine(mesh, cells=[3])
    tags = np.asarray(out.cell_data["tag"][0]).reshape(-1)
    # Every value still appears, and only values that were in the input.
    assert set(tags) == set(range(n))


def test_cell_sets_expand_under_a_selective_pass():
    mesh = _tri_grid(2)
    mesh.cell_sets["patch"] = [np.array([3])]
    out = refine(mesh, cells=[3])
    assert len(out.cell_sets["patch"][0]) == 4


def test_regions_survive_a_selective_pass():
    mesh = _tri_grid(2)
    mesh.cell_sets["patch"] = [np.array([3])]
    out = refine(mesh, cells=[3], record_levels=True)
    assert "patch" in out.cell_sets


@pytest.mark.parametrize(
    "kwargs",
    [
        {"cells": [0], "region": "anything"},
        {"cells": [1000]},
        {"cells": [-1]},
        {"region": "nope"},
        {"where": "absent < 1"},
        {"where": "not-an-expression"},
        {"closure": "blue"},
    ],
)
def test_bad_selectors_raise(kwargs):
    with pytest.raises(ValueError):
        refine(_tri_grid(2), **kwargs)


def test_selective_numbering_is_stable_across_repeated_runs():
    mesh = _quad_grid(5)
    first = refine(mesh, cells=[7, 18])
    for _ in range(5):
        again = refine(mesh, cells=[7, 18])
        assert np.array_equal(
            np.asarray(first.cells[0].data), np.asarray(again.cells[0].data)
        )
        assert np.array_equal(np.asarray(first.points), np.asarray(again.points))


# --------------------------------------------------------------------------- #
# the tables, and C++/numpy parity for the selective path                      #
# --------------------------------------------------------------------------- #
def test_mask_tables_match_the_cpp_core():
    """The Python tables ARE the C++ tables, checked rather than transcribed.

    Both sides generate the tetrahedron's 27 admissible masks from six orbit
    representatives by the same 12 even permutations, so this also pins that the
    two generators agree on which representative reaches each mask.
    """
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus import _refine_templates as t

    for cell_type in t.SUPPORTED_TYPES:
        cpp = core.refine_mask_table(cell_type)
        assert set(cpp) == set(t.TABLES[cell_type]), cell_type
        for mask, entry in cpp.items():
            children, children_alt, tie_a, tie_b = t.TABLES[cell_type][mask]
            assert [list(r) for r in children] == entry["children"], (cell_type, mask)
            assert [list(r) for r in children_alt] == entry["children_alt"], (
                cell_type,
                mask,
            )
            if children_alt:
                assert (tie_a, tie_b) == (entry["tie_a"], entry["tie_b"])
        for mask in range(t.FULL_MASK[cell_type] + 1):
            for propagate in (False, True):
                assert t.promote_mask(
                    cell_type, mask, propagate
                ) == core.refine_promote_mask(cell_type, mask, propagate), (
                    cell_type,
                    mask,
                    propagate,
                )


@pytest.mark.parametrize(
    "factory,n,selected", [(_tri_grid, 3, 8), (_quad_grid, 4, 5), (_hex_grid, 3, 13)]
)
@pytest.mark.parametrize("closure", ["redgreen", "propagate", "balanced"])
@pytest.mark.parametrize("levels", [1, 2])
def test_cpp_matches_python_selective(factory, n, selected, closure, levels):
    """Byte-identical across the C++-core/numpy-fallback boundary, selective too."""
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._refine import _refine_py, _resolve_selection
    from meshioplusplus._refine_templates import closure_from_name

    mesh = factory(n)
    ids = np.array([selected], dtype=np.int64)
    got = core.refine(mesh, levels, False, ids, "", "", "<", 0.0, closure, True)["mesh"]
    seed = _resolve_selection(mesh, ids, None, "", "<", 0.0)
    want, _, _ = _refine_py(mesh, levels, False, seed, closure_from_name(closure), True)

    assert [cb.type for cb in got.cells] == [cb.type for cb in want.cells]
    for a, b in zip(got.cells, want.cells):
        assert np.array_equal(np.asarray(a.data), np.asarray(b.data))
    assert np.array_equal(np.asarray(got.points), np.asarray(want.points))
    for a, b in zip(got.cell_data["refine:level"], want.cell_data["refine:level"]):
        assert np.array_equal(np.asarray(a).reshape(-1), np.asarray(b).reshape(-1))
    # The two reserved point arrays too: they decide whether a later pass reuses
    # a node or tears the mesh, so the two engines disagreeing about them would
    # make every *subsequent* pass diverge rather than this one.
    assert set(got.point_data) == set(want.point_data)
    for name in ("refine:hanging", "refine:entity"):
        if name in got.point_data:
            a = np.asarray(got.point_data[name])
            b = np.asarray(want.point_data[name])
            assert a.dtype == b.dtype and a.shape == b.shape, name
            assert a.tobytes() == b.tobytes(), name


@pytest.mark.parametrize(
    "factory,n,selected", [(_tri_grid, 3, 8), (_quad_grid, 4, 5), (_hex_grid, 3, 13)]
)
@pytest.mark.parametrize("levels", [1, 2])
def test_cpp_matches_python_hierarchy(factory, n, selected, levels):
    """Byte-identical refine:cell_id/refine:parent_id, uniform and selective.

    The pre-existing parity test above compares only geometry, level and the
    two pre-existing reserved point arrays -- it would not have caught a
    cell_data divergence in the new hierarchy arrays, so this is a dedicated
    element-wise comparison, not an addition to the existing assertions.
    """
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._refine import _refine_py, _resolve_selection

    mesh = factory(n)
    ids = np.array([selected], dtype=np.int64)
    got = core.refine(
        mesh, levels, False, ids, "", "", "<", 0.0, "redgreen", False, True
    )["mesh"]
    seed = _resolve_selection(mesh, ids, None, "", "<", 0.0)
    want, _, _ = _refine_py(mesh, levels, False, seed, "redgreen", False, True)

    for name in ("refine:cell_id", "refine:parent_id"):
        assert name in got.cell_data and name in want.cell_data
        for a, b in zip(got.cell_data[name], want.cell_data[name]):
            a = np.asarray(a)
            b = np.asarray(b)
            assert a.dtype == b.dtype and a.shape == b.shape, name
            assert a.tobytes() == b.tobytes(), name
    assert np.array_equal(
        got.point_data["refine:entity"], want.point_data["refine:entity"]
    )


@pytest.mark.parametrize("factory", [_two_triangles, _unit_quad, _unit_tet, _unit_cube])
def test_cpp_matches_python_hierarchy_uniform(factory):
    """The uniform path too: every cell splits, so rule 1 (keep the id) never
    fires, which is exactly the case
    test_cpp_matches_python_hierarchy's selected-cell fixtures do not exercise
    on their own (some cells there are always left untouched)."""
    core = pytest.importorskip("meshioplusplus._core")
    from meshioplusplus._refine import _refine_py

    mesh = factory()
    got = core.refine(mesh, 2, False, None, "", "", "<", 0.0, "redgreen", False, True)[
        "mesh"
    ]
    want, _, _ = _refine_py(mesh, 2, False, None, "redgreen", False, True)
    for name in ("refine:cell_id", "refine:parent_id"):
        for a, b in zip(got.cell_data[name], want.cell_data[name]):
            assert np.asarray(a).tobytes() == np.asarray(b).tobytes(), name


# --------------------------------------------------------------------------- #
# the balanced closure: 2:1 balance, hanging nodes kept                        #
# --------------------------------------------------------------------------- #
def _max_level_gap(mesh):
    """Largest level difference between two cells sharing a NODE (the 2:1 rule).

    Node adjacency, not edge adjacency: across a hanging interface the coarse
    cell spans a whole edge while the fine cell has only half of it, so the two
    are different entities and an edge-keyed check is blind to exactly the
    coarse/fine adjacency 2:1 balance exists to police.
    """
    levels = np.concatenate(
        [np.asarray(a).reshape(-1) for a in mesh.cell_data["refine:level"]]
    )
    npts = len(mesh.points)
    lo = np.full(npts, 1 << 30, dtype=np.int64)
    hi = np.full(npts, -(1 << 30), dtype=np.int64)
    base = 0
    for cb in mesh.cells:
        data = np.asarray(cb.data)
        for c in range(len(data)):
            np.minimum.at(lo, data[c], levels[base + c])
            np.maximum.at(hi, data[c], levels[base + c])
        base += len(data)
    seen = hi > -(1 << 30)
    return int((hi[seen] - lo[seen]).max())


def _hanging_node_ids(mesh):
    """The DISTINCT constrained nodes: a node on some cell's edge midpoint OR at
    the centre of one of its quadrilateral faces.

    The face centres matter as much as the edge midpoints: a split face leaves
    one on a neighbour that still spans the face whole, and a solver has to
    constrain it too.

    A node the cell itself references is not constrained *by that cell*, which
    only starts to matter once a pass refines more than one level: a child's
    corner can then coincide with a coarser neighbour's edge midpoint.
    """
    from meshioplusplus._refine_templates import EDGES, QUAD_FACES

    points = np.asarray(mesh.points)
    at = {tuple(p): i for i, p in enumerate(points)}
    out = set()
    for cb in mesh.cells:
        data = np.asarray(cb.data)
        for c, row in enumerate(data):
            own = set(int(x) for x in row)
            spots = [
                tuple((points[row[a]] + points[row[b]]) * 0.5)
                for a, b in EDGES[cb.type]
            ]
            spots += [
                tuple(sum(points[row[i]] for i in face) / 4.0)
                for face in QUAD_FACES[cb.type]
            ]
            for m in spots:
                hit = at.get(m)
                if hit is not None and hit not in own:
                    out.add(hit)
    return out


def _torn_positions(mesh):
    """Positions carrying two distinct node ids that cells both reference.

    Distinct from a hanging node: the mesh is torn, and no `refine:hanging` flag
    can describe it because both nodes are genuinely used.
    """
    used = set()
    for cb in mesh.cells:
        used.update(int(x) for x in np.asarray(cb.data).reshape(-1))
    seen = {}
    for i, p in enumerate(np.asarray(mesh.points)):
        if i in used:
            seen.setdefault(tuple(p), []).append(i)
    return [ids for ids in seen.values() if len(ids) > 1]


def test_balanced_keeps_hanging_nodes_instead_of_closing():
    mesh = _hex_grid(4)
    out = refine(mesh, cells=[21], closure="balanced", record_levels=True)
    # One cell split into 8, every other cell untouched: 64 - 1 + 8.
    assert len(out.cells[0].data) == 71
    assert out.cells[0].type == "hexahedron"
    # It is deliberately NOT conforming, and says so in refine:hanging.
    assert _count_hanging_nodes(out) > 0
    assert "refine:hanging" in out.point_data
    # refine:hanging must mark EXACTLY the constrained nodes -- not a superset
    # (which would over-constrain a solver) and not a subset (which would leave
    # a crack).
    flags = np.asarray(out.point_data["refine:hanging"]).reshape(-1)
    assert set(np.flatnonzero(flags).tolist()) == _hanging_node_ids(out)


def test_the_conforming_closures_leave_no_hanging_array():
    mesh = _hex_grid(3)
    for closure in ("redgreen", "propagate"):
        out = refine(mesh, cells=[13], closure=closure, record_levels=True)
        _assert_conforming(out)
        assert "refine:hanging" not in out.point_data


# --------------------------------------------------------------------------- #
# multi-pass balanced: tearing, and reporting every constrained node           #
# --------------------------------------------------------------------------- #
def test_the_tear_oracle_actually_fires():
    """A mesh with two referenced nodes at one position must be caught."""
    torn = meshioplusplus.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [1, 0, 0]], dtype=float),
        [("triangle", np.array([[0, 1, 2], [0, 2, 3], [4, 2, 3]], dtype=np.int64))],
    )
    assert len(_torn_positions(torn)) == 1
    assert _torn_positions(_hex_grid(2)) == []


@pytest.mark.parametrize("levels", [2, 3])
def test_multi_pass_balanced_does_not_tear_the_mesh(levels):
    """Refining one cell twice draws its coarse neighbours in by the balance rule.

    Those already carry hanging nodes from the first pass, and their own
    refinement must reuse them rather than allocate a coincident second node --
    which would leave the two sides of the interface referencing different nodes
    at the same point. That is what ``refine:entity`` is for.
    """
    out = refine(
        _hex_grid(4), cells=[0], levels=levels, closure="balanced", record_levels=True
    )
    assert _torn_positions(out) == []
    assert _max_level_gap(out) <= 1


@pytest.mark.parametrize("levels", [1, 2, 3])
def test_multi_pass_balanced_reports_every_constrained_node(levels):
    """Neither a superset nor a subset, which is what doc/refine.md promises.

    Stating the rule over the *input* cells' entities under-reports here: a cell
    the balance rule draws in has children whose sub-edges the input cell never
    had, so a node the neighbouring refinement places on one of them is
    constrained without any input entity ever naming it.
    """
    out = refine(
        _hex_grid(4), cells=[0], levels=levels, closure="balanced", record_levels=True
    )
    flags = np.asarray(out.point_data["refine:hanging"]).reshape(-1)
    assert set(np.flatnonzero(flags).tolist()) == _hanging_node_ids(out)


def test_entity_keys_describe_the_nodes_they_name():
    out = refine(_hex_grid(3), cells=[13], closure="balanced")
    keys = np.asarray(out.point_data["refine:entity"])
    points = np.asarray(out.points)
    assert keys.shape == (len(points), 4)
    keyed = 0
    for i, key in enumerate(keys):
        if key[3] < 0:
            continue  # sentinel: an original point or a body centre
        keyed += 1
        corners = key[2:] if key[0] < 0 else key
        assert np.allclose(
            points[i], points[list(corners)].mean(axis=0), rtol=0, atol=0
        )
    assert keyed > 0


def test_the_conforming_closures_attach_no_entity_array():
    """They cannot tear, so they must not pay for the bookkeeping that prevents it."""
    mesh = _hex_grid(3)
    for closure in ("redgreen", "propagate"):
        out = refine(mesh, cells=[13], levels=2, closure=closure)
        assert "refine:entity" not in out.point_data
    assert "refine:entity" not in refine(mesh, levels=2).point_data


def test_a_stale_entity_array_is_ignored_rather_than_trusted():
    """The keys name point indices, which refine never renumbers but other
    operations do. A key that no longer reproduces its own point's coordinates
    invalidates the whole array: refine falls back to allocating fresh nodes,
    which is the old behaviour rather than a wrong answer."""
    stale = _hex_grid(3)
    keys = np.full((len(stale.points), 4), -1, dtype=np.int64)
    keys[0] = (-1, -1, 1, 2)  # point 0 is a corner, not that edge's midpoint
    stale.point_data["refine:entity"] = keys

    with pytest.warns(UserWarning, match="refine:entity"):
        from meshioplusplus._refine import _read_entity_keys

        assert _read_entity_keys(stale) is None

    got = refine(stale, cells=[13], closure="balanced")
    want = refine(_hex_grid(3), cells=[13], closure="balanced")
    assert np.array_equal(np.asarray(got.points), np.asarray(want.points))
    assert np.array_equal(np.asarray(got.cells[0].data), np.asarray(want.cells[0].data))


def test_balanced_does_not_propagate_on_a_mesh_of_uniform_level():
    # The whole point: on a mesh where every cell is at the same level, refining
    # one cell puts nothing else out of balance, so nothing else is touched.
    for factory, n, sel in [(_tri_grid, 4, 8), (_quad_grid, 5, 7), (_hex_grid, 4, 21)]:
        mesh = factory(n)
        before = len(mesh.cells[0].data)
        out = refine(mesh, cells=[sel], closure="balanced", record_levels=True)
        children = {"triangle": 4, "quad": 4, "hexahedron": 8}[mesh.cells[0].type]
        assert len(out.cells[0].data) == before - 1 + children
        assert _max_level_gap(out) <= 1


def test_balanced_draws_in_neighbours_that_would_fall_two_levels_behind():
    mesh = _hex_grid(4)
    first = refine(mesh, cells=[21], closure="balanced", record_levels=True)
    levels = np.asarray(first.cell_data["refine:level"][0]).reshape(-1)

    # Refining a level-1 child would leave its level-0 neighbours two levels
    # coarser, so balancing must promote them -- and only them.
    child = int(np.flatnonzero(levels == 1)[0])
    second = refine(first, cells=[child], closure="balanced", record_levels=True)
    grew = len(second.cells[0].data) - len(first.cells[0].data)
    assert grew > 7, "balancing must have drawn in cells beyond the selected one"
    assert grew < 7 * 30, "and must stay local, not reach the whole mesh"
    assert _max_level_gap(second) <= 1
    assert set(np.asarray(second.cell_data["refine:level"][0]).reshape(-1)) == {0, 1, 2}


def test_balanced_is_far_cheaper_than_propagate():
    mesh = _hex_grid(3)
    balanced = refine(mesh, cells=[13], closure="balanced")
    propagated = refine(mesh, cells=[13], closure="propagate")
    assert len(balanced.cells[0].data) < len(propagated.cells[0].data) / 4


def test_balanced_parses_by_name():
    mesh = _tri_grid(2)
    a = refine(mesh, cells=[0], closure="balanced")
    b = refine(mesh, cells=[0], closure="2:1")
    assert np.array_equal(np.asarray(a.cells[0].data), np.asarray(b.cells[0].data))
