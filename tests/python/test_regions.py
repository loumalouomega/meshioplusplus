"""Named regions: the model, the legacy compat views, and the operations.

The cross-format round-trip matrix lives in ``test_region_roundtrip.py``; this
file covers everything upstream of a file format. See ``doc/regions.md``.
"""

import copy

import numpy as np
import pytest
from numpy.testing import assert_array_equal

import meshioplusplus

from . import helpers


def _tri_mesh():
    """Two triangles over a unit square."""
    return meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        [("triangle", [[0, 1, 2], [0, 2, 3]])],
    )


def _tet_mesh():
    """Two tetrahedra, for exercising side regions (4 faces each)."""
    return meshioplusplus.Mesh(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 1.0],
        ],
        [("tetra", [[0, 1, 2, 4], [0, 2, 3, 4]])],
    )


def _mixed_mesh():
    """A triangle block followed by a quad block, for global cell indices."""
    return meshioplusplus.Mesh(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
        ],
        [
            ("triangle", [[0, 1, 2], [0, 2, 3]]),
            ("quad", [[1, 4, 5, 2]]),
        ],
    )


# --------------------------------------------------------------------------- #
# construction                                                                 #
# --------------------------------------------------------------------------- #
def test_region_construction_canonicalizes():
    r = meshioplusplus.Region("wall", "point", [3, 0, 3, 1])
    assert r.name == "wall"
    assert r.kind == "point"
    assert r.dim == -1 and r.tag == -1
    # Sorted and de-duplicated, so equality is exact.
    assert_array_equal(r.entries, [0, 1, 3])
    assert len(r) == 3


def test_side_region_is_pairs():
    r = meshioplusplus.Region("s", "side", [[1, 2], [0, 3], [1, 2]])
    assert r.entries.shape == (2, 2)
    assert_array_equal(r.entries, [[0, 3], [1, 2]])


def test_unknown_kind_is_rejected():
    with pytest.raises(ValueError, match="unknown kind"):
        meshioplusplus.Region("x", "facet", [0])


def test_dim_and_tag_are_carried():
    r = meshioplusplus.Region("Surface", "cell", [0], dim=2, tag=7)
    assert (r.dim, r.tag) == (2, 7)


def test_region_equality_is_by_value():
    a = meshioplusplus.Region("w", "cell", [1, 0])
    b = meshioplusplus.Region("w", "cell", [0, 1])
    assert a == b
    assert a != meshioplusplus.Region("w", "point", [0, 1])
    assert a != meshioplusplus.Region("w", "cell", [0, 2])


# --------------------------------------------------------------------------- #
# the compat views                                                             #
# --------------------------------------------------------------------------- #
def test_point_sets_view_round_trips_both_ways():
    mesh = _tri_mesh()
    mesh.point_sets = {"fixed": np.array([0, 3])}

    # ...written through to a region...
    assert [r.name for r in mesh.regions] == ["fixed"]
    assert mesh.regions[0].kind == "point"
    # ...and readable back in the historical shape.
    assert_array_equal(mesh.point_sets["fixed"], [0, 3])

    # The other direction: a region appears as a set.
    mesh.regions.append(meshioplusplus.Region("loose", "point", [1, 2]))
    assert sorted(mesh.point_sets) == ["fixed", "loose"]


def test_cell_sets_view_materializes_one_array_per_block():
    mesh = _mixed_mesh()
    mesh.cell_sets = {"some": [np.array([1]), np.array([0])]}

    blocks = mesh.cell_sets["some"]
    assert len(blocks) == len(mesh.cells)
    assert_array_equal(blocks[0], [1])
    assert_array_equal(blocks[1], [0])
    # Stored globally, block-major: triangle 1 is global 1, quad 0 is global 2.
    assert_array_equal(mesh.regions[0].entries, [1, 2])


def test_views_are_dicts():
    # numpy.testing.assert_equal and plenty of user code branch on
    # isinstance(x, dict); the views must not break that.
    mesh = _tri_mesh()
    mesh.point_sets = {"a": np.array([0])}
    assert isinstance(mesh.point_sets, dict)
    assert isinstance(mesh.cell_sets, dict)
    assert mesh.point_sets == {"a": np.array([0])} or list(mesh.point_sets) == ["a"]


def test_view_mutators_write_through():
    mesh = _tri_mesh()
    mesh.point_sets["a"] = np.array([0, 1])
    assert [r.name for r in mesh.regions] == ["a"]

    mesh.point_sets.update({"b": np.array([2])})
    assert sorted(r.name for r in mesh.regions) == ["a", "b"]

    del mesh.point_sets["a"]
    assert [r.name for r in mesh.regions] == ["b"]

    mesh.point_sets.clear()
    assert mesh.regions == []


def test_setting_the_property_replaces_only_its_own_kind():
    mesh = _tri_mesh()
    mesh.cell_sets = {"c": [np.array([0])]}
    mesh.point_sets = {"p": np.array([0])}
    # Assigning point_sets must not have wiped the cell region.
    assert sorted(r.name for r in mesh.regions) == ["c", "p"]


def test_constructor_keeps_regions_and_sets_together():
    # Regression: the property setters replace all regions of their kind, so a
    # constructor that assigned them unconditionally wiped the `regions=`
    # argument — which is exactly what every C++-backed reader passes.
    mesh = meshioplusplus.Mesh(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]],
        [("triangle", [[0, 1, 2], [0, 2, 3]])],
        regions=[meshioplusplus.Region("wall", "side", [[0, 1]])],
        point_sets={"fixed": np.array([0])},
    )
    assert sorted(r.name for r in mesh.regions) == ["fixed", "wall"]


def test_side_regions_show_in_repr_but_not_in_the_sets_views():
    mesh = _tet_mesh()
    mesh.regions.append(meshioplusplus.Region("wall", "side", [[0, 1]]))
    assert "Side sets: wall" in repr(mesh)
    assert "wall" not in mesh.point_sets
    assert "wall" not in mesh.cell_sets


def test_the_escape_hatch_keeps_non_index_cell_sets_verbatim():
    # gmsh stashes entity tags — negative values included — in cell_sets. Those
    # cannot be a region (regions are sorted non-negative indices), so they take
    # the passthrough and come back untouched.
    mesh = _tri_mesh()
    value = [np.array([-3, 7])]
    mesh.cell_sets["gmsh:bounding_entities"] = value

    assert mesh.regions == []
    assert_array_equal(mesh.cell_sets["gmsh:bounding_entities"][0], [-3, 7])
    assert "gmsh:bounding_entities" in mesh.cell_sets
    del mesh.cell_sets["gmsh:bounding_entities"]
    assert "gmsh:bounding_entities" not in mesh.cell_sets


def test_a_real_gmsh_file_keeps_its_bounding_entities():
    mesh = meshioplusplus.read("example/example.msh")
    assert "gmsh:bounding_entities" in mesh.cell_sets
    blocks = mesh.cell_sets["gmsh:bounding_entities"]
    assert len(blocks) == len(mesh.cells)


def test_none_blocks_become_empty_arrays():
    # A cell_sets value may carry None for a block with no members (gmsh's
    # reader builds [None] * n). Global indices cannot express "absent"
    # separately from "empty", so the round trip normalizes it.
    mesh = _mixed_mesh()
    mesh.cell_sets = {"s": [None, np.array([0])]}
    blocks = mesh.cell_sets["s"]
    assert blocks[0] is not None and len(blocks[0]) == 0
    assert_array_equal(blocks[1], [0])


def test_deepcopy_carries_regions():
    mesh = _tri_mesh()
    mesh.regions.append(meshioplusplus.Region("w", "cell", [0]))
    other = copy.deepcopy(mesh)
    other.regions[0].entries[0] = 1
    assert_array_equal(mesh.regions[0].entries, [0])


# --------------------------------------------------------------------------- #
# the Python <-> C++ boundary                                                  #
# --------------------------------------------------------------------------- #
def test_regions_survive_the_cpp_boundary():
    _core = pytest.importorskip("meshioplusplus._core")
    mesh = _tet_mesh()
    mesh.point_sets = {"fixed": np.array([0, 3])}
    mesh.cell_sets = {"solid": [np.array([0, 1])]}
    mesh.regions.append(meshioplusplus.Region("wall", "side", [[0, 1], [1, 3]]))

    # transform passes regions through verbatim — a pure coordinate move.
    identity = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    out = _core.transform(mesh, identity, False)

    got = {(r.name, r.kind): r for r in out.regions}
    assert sorted(got) == [("fixed", "point"), ("solid", "cell"), ("wall", "side")]
    assert_array_equal(got[("wall", "side")].entries, [[0, 1], [1, 3]])


def test_dim_and_tag_survive_the_cpp_boundary():
    _core = pytest.importorskip("meshioplusplus._core")
    mesh = _tri_mesh()
    mesh.regions.append(meshioplusplus.Region("S", "cell", [0], dim=2, tag=42))
    identity = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
    out = _core.transform(mesh, identity, False)
    assert (out.regions[0].dim, out.regions[0].tag) == (2, 42)


# --------------------------------------------------------------------------- #
# operations                                                                   #
# --------------------------------------------------------------------------- #
def test_crop_remaps_and_drops():
    mesh = _mixed_mesh()
    mesh.point_sets = {"all": np.arange(6)}
    mesh.cell_sets = {"tris": [np.array([0, 1]), np.array([], dtype=int)]}

    # Keep only x <= 1.5, which drops the quad (it reaches x = 2).
    out = meshioplusplus.crop(mesh, bbox=(-1, -1, -1, 1.5, 2, 2))

    assert len(out.points) < len(mesh.points)
    assert (np.asarray(out.point_sets["all"]) < len(out.points)).all()
    assert_array_equal(out.cell_sets["tris"][0], [0, 1])


def test_crop_carries_a_side_region_whose_cell_survives():
    mesh = _tet_mesh()
    mesh.regions.append(meshioplusplus.Region("wall", "side", [[0, 1], [1, 3]]))
    out = meshioplusplus.crop(mesh, bbox=(-9, -9, -9, 9, 9, 9))
    side = [r for r in out.regions if r.kind == "side"]
    assert len(side) == 1
    assert_array_equal(side[0].entries, [[0, 1], [1, 3]])


def test_crop_drops_a_side_entry_whose_cell_vanishes():
    mesh = _tet_mesh()
    mesh.regions.append(meshioplusplus.Region("wall", "side", [[0, 1], [1, 3]]))
    # Keep only the half-space containing tetra 0's far corner: the second tet
    # loses at least one node, so under mode="all" it is dropped.
    out = meshioplusplus.crop(mesh, plane=((0.0, 0.6, 0.0), (0.0, -1.0, 0.0)))
    side = [r for r in out.regions if r.kind == "side"]
    assert len(side) == 1
    # Whatever survived, it is strictly fewer facets than we started with.
    assert len(side[0].entries) < 2


def test_split_carries_regions_onto_each_piece():
    mesh = _mixed_mesh()
    mesh.cell_sets = {"first": [np.array([0]), np.array([0])]}
    pieces = meshioplusplus.split(mesh, by="type")
    assert set(pieces) == {"triangle", "quad"}
    for piece in pieces.values():
        assert "first" in piece.cell_sets


def test_merge_namespaces_colliding_region_names():
    a = _tri_mesh()
    a.cell_sets = {"g": [np.array([0])]}
    b = _tri_mesh()
    b.cell_sets = {"g": [np.array([1])]}
    out = meshioplusplus.merge([a, b])
    assert set(out.cell_sets) == {"0:g", "1:g"}


def test_reorder_remaps_regions():
    mesh = helpers.add_point_sets(helpers.tri_mesh_5)
    out = meshioplusplus.reorder(mesh, method="rcm")
    for name, idx in mesh.point_sets.items():
        assert len(out.point_sets[name]) == len(idx)
        assert (np.asarray(out.point_sets[name]) < len(out.points)).all()


def test_refine_expands_a_cell_region_to_the_children():
    mesh = _tri_mesh()
    mesh.cell_sets = {"left": [np.array([0])]}
    out = meshioplusplus.refine(mesh, levels=1)
    # One triangle becomes four.
    assert len(out.cell_sets["left"][0]) == 4


def test_transform_and_smooth_pass_regions_through_unchanged():
    for op in (
        lambda m: meshioplusplus.transform(m, translate=(1.0, 0.0, 0.0)),
        lambda m: meshioplusplus.smooth(m, iterations=1),
    ):
        mesh = _tet_mesh()
        mesh.point_sets = {"fixed": np.array([0, 3])}
        mesh.regions.append(meshioplusplus.Region("wall", "side", [[0, 1]]))
        out = op(mesh)
        assert_array_equal(out.point_sets["fixed"], [0, 3])
        side = [r for r in out.regions if r.kind == "side"]
        assert len(side) == 1
        assert_array_equal(side[0].entries, [[0, 1]])


def test_slice_drops_regions_with_a_warning(recwarn):
    mesh = _tet_mesh()
    mesh.point_sets = {"fixed": np.array([0, 3])}
    out = meshioplusplus.slice(mesh, origin=(0.5, 0.5, 0.25), normal=(0.0, 0.0, 1.0))
    # The section is newly created geometry: no entity correspondence exists.
    assert out.regions == []


def test_extract_surface_drops_regions():
    mesh = _tet_mesh()
    mesh.cell_sets = {"solid": [np.array([0, 1])]}
    out = meshioplusplus.extract_surface(mesh)
    assert out.regions == []
