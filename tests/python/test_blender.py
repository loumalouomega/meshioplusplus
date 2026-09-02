"""Tests for the Blender bridge.

Two halves, deliberately — the ``test_interop.py`` split, for the same reason.

The **pure payload** tests run in the default CI matrix, where Blender is
nowhere in sight. That is what ``_blender``'s ``_interop``-style split buys:
the loop/polygon layout, the block-major indexing, the reduction that keeps
n-gons, and the attribute typing are all exercised with no ``bpy`` at all.

The **gated** tests use ``pytest.importorskip("bpy")`` and assert what only a
real Blender data model can settle — chiefly that ``mesh.validate()`` finds
nothing to repair, which is Blender's own topology checker used as an oracle.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _blender, _interop
from meshioplusplus._blender import (
    FIELD_PREFIX,
    _attribute_arrays,
    _mesh_from_blender_arrays,
    _to_blender_payload,
)
from meshioplusplus._mesh import CellBlock, Mesh
from meshioplusplus._regions import Region


# --------------------------------------------------------------------------- #
# fixtures                                                                    #
# --------------------------------------------------------------------------- #
def _quad_mesh():
    return Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
        [("quad", np.array([[0, 1, 2, 3]], dtype=np.int64))],
    )


def _two_tets():
    """Two tetrahedra sharing the face (1, 2, 3), tagged 10 and 20."""
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ]
    )
    return Mesh(
        points,
        [("tetra", np.array([[0, 1, 2, 3], [1, 2, 3, 4]], dtype=np.int64))],
        cell_data={"mat": [np.array([10, 20])]},
    )


def _mixed_surface():
    """triangle + quad + pentagon, with point and cell data and two regions."""
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [3.0, 0.5, 0.0],
            [2.5, 1.5, 0.0],
            [1.5, 1.5, 0.0],
        ]
    )
    return Mesh(
        points,
        [
            ("triangle", np.array([[0, 1, 5]], dtype=np.int64)),
            ("quad", np.array([[1, 2, 3, 4]], dtype=np.int64)),
            CellBlock("polygon", np.array([[3, 6, 7, 8, 4]], dtype=np.int64)),
        ],
        point_data={"T": np.arange(9, dtype=np.float64)},
        cell_data={"c": [np.array([1.0]), np.array([2.0]), np.array([3.0])]},
        field_data={"step": np.array([7])},
        regions=[
            Region("inlet", "point", [0, 1], dim=0, tag=5),
            Region("wall", "cell", [1, 2], dim=2, tag=9),
        ],
    )


# --------------------------------------------------------------------------- #
# pure payload — the reduction                                                #
# --------------------------------------------------------------------------- #
def test_payload_keeps_a_quad_as_a_quad():
    """The headline oracle: Blender holds n-gons, so nothing is simplexified.

    This is the one assertion that distinguishes the Blender reduction from
    ``to_trimesh``'s, whose ``_triangulate`` ends in ``simplexify`` because a
    ``Trimesh`` holds triangles only.
    """
    payload = _to_blender_payload(_quad_mesh())
    assert payload.num_polygons == 1
    assert payload.loop_totals.tolist() == [4]
    assert payload.loop_vertices.tolist() == [0, 1, 2, 3]
    assert "simplexify" not in payload.ops


def test_payload_routes_a_volume_through_extract_surface():
    payload = _to_blender_payload(_two_tets())
    assert payload.ops == ["extract_surface"]
    # Two tets sharing one face: 8 faces, 2 of them interior.
    assert payload.num_polygons == 6
    assert set(payload.loop_totals.tolist()) == {3}


def test_payload_gathers_cell_data_through_parents():
    """Each boundary face carries its owning tet's value, not a default.

    Pins the ``record_parent_ids`` gather in ``_viewer_browser._renderable_surface``
    — without it a solid would arrive in Blender with no material at all.
    """
    payload = _to_blender_payload(_two_tets())
    attr = payload.face_attributes["mat"]
    assert attr.data_type == "INT"
    assert sorted(attr.values.tolist()) == [10, 10, 10, 20, 20, 20]


def test_payload_linearizes_higher_order():
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    )
    mesh = Mesh(points, [("triangle6", np.array([[0, 1, 2, 3, 4, 5]], dtype=np.int64))])
    payload = _to_blender_payload(mesh)
    assert payload.ops == ["linearize"]
    assert payload.loop_totals.tolist() == [3]


def test_payload_ops_ignore_a_supplied_skin_on_a_volume_mesh():
    """extract_surface drops the supplied 2-D block, so nothing is linearized."""
    mesh = Mesh(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.5, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.0],
            ]
        ),
        [
            ("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64)),
            ("triangle6", np.array([[0, 1, 2, 4, 5, 6]], dtype=np.int64)),
        ],
    )
    payload = _to_blender_payload(mesh)
    assert payload.ops == ["extract_surface"]


def test_payload_ragged_polygon_block():
    """A ragged block's ``data`` is a Python list, not an ndarray."""
    mesh = Mesh(
        np.zeros((7, 3)),
        [CellBlock("polygon", [[0, 1, 2, 3], [0, 1, 2, 3, 4], [0, 1, 2, 3, 4, 5]])],
    )
    payload = _to_blender_payload(mesh)
    assert payload.loop_totals.tolist() == [4, 5, 6]
    assert payload.loop_starts.tolist() == [0, 4, 9]
    assert payload.loop_vertices.tolist() == [0, 1, 2, 3] + [0, 1, 2, 3, 4] + [
        0,
        1,
        2,
        3,
        4,
        5,
    ]


@pytest.mark.parametrize(
    "factory",
    [_quad_mesh, _two_tets, _mixed_surface],
    ids=["quad", "two_tets", "mixed_surface"],
)
def test_payload_loop_polygon_consistency(factory):
    """The structural invariant Blender crashes on if it is broken."""
    payload = _to_blender_payload(factory())
    starts, totals = payload.loop_starts, payload.loop_totals
    assert len(starts) == len(totals) == payload.num_polygons
    assert starts[0] == 0
    assert starts[1:].tolist() == np.cumsum(totals)[:-1].tolist()
    assert int(starts[-1]) + int(totals[-1]) == len(payload.loop_vertices)
    assert payload.loop_vertices.min() >= 0
    assert payload.loop_vertices.max() < len(payload.vertices)
    assert payload.vertices.dtype == np.float32
    assert payload.vertices.shape[1] == 3
    for array in (
        payload.vertices,
        payload.loop_vertices,
        payload.loop_starts,
        payload.loop_totals,
    ):
        assert array.flags["C_CONTIGUOUS"]
    for array in (payload.loop_vertices, payload.loop_starts, payload.loop_totals):
        assert array.dtype == np.int32


def test_payload_pads_2d_points():
    mesh = Mesh(
        np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        [("triangle", np.array([[0, 1, 2]], dtype=np.int64))],
    )
    payload = _to_blender_payload(mesh)
    assert payload.vertices.shape == (3, 3)
    assert np.all(payload.vertices[:, 2] == 0.0)
    assert any("2D" in note for note in payload.notes)


def test_payload_line_blocks_become_edges():
    mesh = Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]),
        [
            ("triangle", np.array([[0, 1, 2]], dtype=np.int64)),
            ("line", np.array([[0, 2]], dtype=np.int64)),
        ],
    )
    payload = _to_blender_payload(mesh)
    assert payload.edge_vertices.tolist() == [[0, 2]]
    assert payload.edge_vertices.dtype == np.int32
    assert payload.num_polygons == 1


def test_payload_vertex_blocks_are_noted_not_dropped_silently():
    mesh = Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]),
        [
            ("triangle", np.array([[0, 1, 2]], dtype=np.int64)),
            ("vertex", np.array([[1]], dtype=np.int64)),
        ],
    )
    payload = _to_blender_payload(mesh)
    assert payload.dropped == []
    assert any("no geometry" in note for note in payload.notes)


def test_payload_drops_unsupported_blocks_with_a_note():
    """A block Blender cannot hold costs a note, never an exception.

    ``VTK_LAGRANGE_TRIANGLE`` is the honest fixture here: it is 2-D, so it
    reaches the classifier untouched, and it is absent from ``_LINEAR_BASE``,
    so ``linearize`` cannot reduce it to something Blender holds.
    """
    mesh = Mesh(
        np.zeros((6, 3)),
        [
            ("triangle", np.array([[0, 1, 2]], dtype=np.int64)),
            CellBlock(
                "VTK_LAGRANGE_TRIANGLE",
                np.array([[0, 1, 2, 3, 4, 5]], dtype=np.int64),
            ),
        ],
    )
    payload = _to_blender_payload(mesh)
    assert payload.num_polygons == 1
    assert [t for _, t in payload.dropped] == ["VTK_LAGRANGE_TRIANGLE"]
    assert any("VTK_LAGRANGE_TRIANGLE" in note for note in payload.notes)


# --------------------------------------------------------------------------- #
# pure payload — attributes                                                   #
# --------------------------------------------------------------------------- #
def _attr_mesh(point_data):
    return Mesh(
        np.zeros((3, 3)),
        [("triangle", np.array([[0, 1, 2]], dtype=np.int64))],
        point_data=point_data,
    )


def test_payload_scalar_attribute_types():
    payload = _to_blender_payload(
        _attr_mesh(
            {
                "f": np.zeros(3, dtype=np.float64),
                "i": np.zeros(3, dtype=np.int64),
                "b": np.zeros(3, dtype=bool),
            }
        )
    )
    got = {k: v.data_type for k, v in payload.point_attributes.items()}
    assert got == {"f": "FLOAT", "i": "INT", "b": "BOOLEAN"}
    assert payload.point_attributes["f"].values.dtype == np.float32
    assert payload.point_attributes["i"].values.dtype == np.int32


def test_payload_vector_attributes_use_blenders_own_types():
    payload = _to_blender_payload(
        _attr_mesh({"v2": np.zeros((3, 2)), "v3": np.zeros((3, 3))})
    )
    assert payload.point_attributes["v2"].data_type == "FLOAT2"
    assert payload.point_attributes["v2"].field == "vector"
    assert payload.point_attributes["v2"].values.size == 6
    assert payload.point_attributes["v3"].data_type == "FLOAT_VECTOR"
    assert payload.point_attributes["v3"].values.size == 9


def test_payload_wide_attributes_expand_with_the_pandas_suffix_rule():
    payload = _to_blender_payload(_attr_mesh({"v": np.zeros((3, 5))}))
    assert sorted(payload.point_attributes) == ["v_0", "v_1", "v_2", "v_3", "v_4"]
    assert all(a.data_type == "FLOAT" for a in payload.point_attributes.values())


def test_payload_integer_vectors_expand_rather_than_widening_to_float():
    payload = _to_blender_payload(_attr_mesh({"ij": np.zeros((3, 3), dtype=np.int64)}))
    assert sorted(payload.point_attributes) == ["ij_0", "ij_1", "ij_2"]
    assert all(a.data_type == "INT" for a in payload.point_attributes.values())


def test_payload_reserved_attribute_names_are_renamed():
    payload = _to_blender_payload(
        _attr_mesh({"position": np.zeros(3), ".hidden": np.zeros(3)})
    )
    assert "attr_position" in payload.point_attributes
    assert "attr.hidden" in payload.point_attributes
    assert "position" not in payload.point_attributes
    assert any("reserved by Blender" in note for note in payload.notes)


def test_payload_out_of_range_integers_fall_back_to_float():
    big = np.array([2**40, 0, 1], dtype=np.int64)
    payload = _to_blender_payload(_attr_mesh({"big": big}))
    assert payload.point_attributes["big"].data_type == "FLOAT"
    assert any("32-bit" in note for note in payload.notes)


def test_payload_string_attributes_are_dropped_with_a_note():
    payload = _to_blender_payload(_attr_mesh({"name": np.array(["a", "b", "c"])}))
    assert payload.point_attributes == {}
    assert any("string array" in note for note in payload.notes)


# --------------------------------------------------------------------------- #
# pure payload — regions and field data                                       #
# --------------------------------------------------------------------------- #
def test_payload_regions_become_boolean_masks_and_a_sidecar():
    payload = _to_blender_payload(_mixed_surface())
    point_key = _interop.REGION_PREFIX + "inlet"
    face_key = _interop.REGION_PREFIX + "wall"
    assert payload.point_attributes[point_key].data_type == "BOOLEAN"
    assert payload.point_attributes[point_key].values.tolist() == [
        True,
        True,
        False,
        False,
        False,
        False,
        False,
        False,
        False,
    ]
    assert payload.face_attributes[face_key].values.tolist() == [False, True, True]
    meta = json.loads(payload.custom_properties[_interop.REGION_META_KEY])
    by_name = {(e["kind"], e["name"]): e for e in meta}
    assert by_name[("point", "inlet")]["tag"] == 5
    assert by_name[("cell", "wall")]["dim"] == 2


def test_payload_volume_input_loses_regions():
    """extract_surface drops regions; the Blender bridge gets no second policy.

    The deliberate mirror of ``test_interop.py::test_trimesh_volume_input_loses_regions``.
    """
    mesh = Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
        regions=[Region("blob", "cell", [0])],
    )
    payload = _to_blender_payload(mesh)
    assert not any(
        k.startswith(_interop.REGION_PREFIX) for k in payload.face_attributes
    )


def test_payload_field_data_becomes_custom_properties():
    payload = _to_blender_payload(_mixed_surface())
    assert payload.custom_properties[FIELD_PREFIX + "step"] == [7]


def test_payload_categories_can_be_switched_off():
    payload = _to_blender_payload(
        _mixed_surface(),
        point_data=False,
        cell_data=False,
        regions=False,
        field_data=False,
    )
    assert payload.point_attributes == {}
    assert payload.face_attributes == {}
    assert payload.custom_properties == {}


# --------------------------------------------------------------------------- #
# pure inverse                                                                #
# --------------------------------------------------------------------------- #
def _round_trip(mesh):
    payload = _to_blender_payload(mesh)
    return payload, _mesh_from_blender_arrays(
        payload.vertices,
        payload.loop_vertices,
        payload.loop_starts,
        payload.loop_totals,
        edge_vertices=payload.edge_vertices,
        point_attributes=_attribute_arrays(
            payload.point_attributes, len(payload.vertices)
        ),
        face_attributes=_attribute_arrays(
            payload.face_attributes, payload.num_polygons
        ),
        custom_properties=payload.custom_properties,
    )


def test_pure_round_trip():
    """The strongest oracle available without Blender."""
    mesh = _mixed_surface()
    payload, back = _round_trip(mesh)
    assert [(b.type, len(b)) for b in back.cells] == [
        ("triangle", 1),
        ("quad", 1),
        ("polygon", 1),
    ]
    assert back.cells[0].data.tolist() == [[0, 1, 5]]
    assert back.cells[1].data.tolist() == [[1, 2, 3, 4]]
    assert back.cells[2].data.tolist() == [[3, 6, 7, 8, 4]]
    assert np.allclose(back.points, mesh.points, atol=1e-6)
    assert np.allclose(back.point_data["T"], mesh.point_data["T"])
    assert [a.tolist() for a in back.cell_data["c"]] == [[1.0], [2.0], [3.0]]
    assert back.field_data["step"].tolist() == [7]
    by_key = {(r.kind, r.name): r for r in back.regions}
    assert by_key[("point", "inlet")].entries.tolist() == [0, 1]
    assert by_key[("cell", "wall")].entries.tolist() == [1, 2]
    assert by_key[("cell", "wall")].tag == 9


def test_inverse_groups_interleaved_polygons_into_one_block_per_size():
    """``vtk_cells_from_data`` splits on runs, so the sort is load-bearing.

    Without the stable sort by side count an alternating triangle/quad mesh
    comes back as one ``CellBlock`` per cell.
    """
    vertices = np.zeros((8, 3), dtype=np.float32)
    loop_vertices = np.array([0, 1, 2, 0, 1, 2, 3, 3, 4, 5, 4, 5, 6, 7], np.int32)
    loop_totals = np.array([3, 4, 3, 4], dtype=np.int32)
    loop_starts = np.array([0, 3, 7, 10], dtype=np.int32)
    back = _mesh_from_blender_arrays(
        vertices,
        loop_vertices,
        loop_starts,
        loop_totals,
        face_attributes={"c": np.array([1.0, 2.0, 3.0, 4.0])},
    )
    assert [(b.type, len(b)) for b in back.cells] == [("triangle", 2), ("quad", 2)]
    assert back.cells[0].data.tolist() == [[0, 1, 2], [3, 4, 5]]
    assert back.cells[1].data.tolist() == [[0, 1, 2, 3], [4, 5, 6, 7]]
    # cell_data follows the same permutation as the cells.
    assert [a.tolist() for a in back.cell_data["c"]] == [[1.0, 3.0], [2.0, 4.0]]


def test_inverse_keeps_only_loose_edges():
    """An edge a polygon already implies must not become a second `line` cell."""
    vertices = np.zeros((4, 3), dtype=np.float32)
    back = _mesh_from_blender_arrays(
        vertices,
        np.array([0, 1, 2], dtype=np.int32),
        np.array([0], dtype=np.int32),
        np.array([3], dtype=np.int32),
        edge_vertices=np.array([[0, 1], [1, 2], [2, 0], [0, 3]], dtype=np.int32),
    )
    assert [(b.type, b.data.tolist()) for b in back.cells] == [
        ("triangle", [[0, 1, 2]]),
        ("line", [[0, 3]]),
    ]


def test_inverse_pads_float_cell_data_over_loose_edges_and_refuses_integers():
    """A Blender edge carries no face value: NaN for floats, refusal for ints."""
    vertices = np.zeros((4, 3), dtype=np.float32)
    notes: list = []
    back = _mesh_from_blender_arrays(
        vertices,
        np.array([0, 1, 2], dtype=np.int32),
        np.array([0], dtype=np.int32),
        np.array([3], dtype=np.int32),
        edge_vertices=np.array([[0, 3]], dtype=np.int32),
        face_attributes={"f": np.array([1.5]), "i": np.array([7], dtype=np.int32)},
        notes=notes,
    )
    assert back.cell_data["f"][0].tolist() == [1.5]
    assert np.isnan(back.cell_data["f"][1]).all()
    assert "i" not in back.cell_data
    assert any("no missing-value representation" in note for note in notes)


def test_inverse_handles_a_mesh_with_no_faces():
    """A point cloud or a curve is an ordinary Blender mesh, not an error.

    ``vtk_cells_from_data`` indexes ``types[0]`` unconditionally, so the empty
    run has to be skipped rather than passed through.
    """
    empty = np.empty(0, dtype=np.int32)
    back = _mesh_from_blender_arrays(np.zeros((0, 3), np.float32), empty, empty, empty)
    assert len(back.points) == 0
    assert back.cells == []

    back = _mesh_from_blender_arrays(
        np.zeros((3, 3), np.float32),
        empty,
        empty,
        empty,
        edge_vertices=np.array([[0, 1], [1, 2]], dtype=np.int32),
    )
    assert [(b.type, b.data.tolist()) for b in back.cells] == [
        ("line", [[0, 1], [1, 2]])
    ]


def test_inverse_rebuilds_regions_without_a_sidecar():
    vertices = np.zeros((3, 3), dtype=np.float32)
    notes: list = []
    back = _mesh_from_blender_arrays(
        vertices,
        np.array([0, 1, 2], dtype=np.int32),
        np.array([0], dtype=np.int32),
        np.array([3], dtype=np.int32),
        point_attributes={_interop.REGION_PREFIX + "tip": np.array([1, 0, 0])},
        notes=notes,
    )
    assert back.regions[0].name == "tip"
    assert back.regions[0].dim == -1
    assert any("dim/tag are unknown" in note for note in notes)


# --------------------------------------------------------------------------- #
# the gate                                                                    #
# --------------------------------------------------------------------------- #
def test_has_blender_is_a_boolean():
    assert isinstance(meshioplusplus.has_blender(), bool)


def test_public_api_is_exported():
    for name in ("to_blender", "from_blender", "has_blender"):
        assert name in meshioplusplus.__all__
        assert hasattr(meshioplusplus, name)


@pytest.mark.skipif(_blender.has_blender(), reason="bpy is installed")
def test_missing_dependency_error_names_blender_and_no_extra():
    """There is no ``meshioplusplus[blender]``, so the message must not claim one."""
    with pytest.raises(ImportError) as excinfo:
        meshioplusplus.to_blender(_quad_mesh())
    message = str(excinfo.value)
    assert "bpy" in message
    assert "doc/blender.md" in message
    assert "meshioplusplus[" not in message
    assert "accelerator" not in message


def test_require_framework_default_message_is_unchanged():
    """Widening ``_require_framework`` must not move the torch/JAX wording."""
    from meshioplusplus._gpu import _require_framework

    with pytest.raises(ImportError) as excinfo:
        _require_framework("to_torch", "definitely_not_a_module", "pip install torch")
    assert "(pick the wheel matching your accelerator). See doc/ml.md." in str(
        excinfo.value
    )


# --------------------------------------------------------------------------- #
# gated: a real Blender data model                                            #
# --------------------------------------------------------------------------- #
def _new_mesh(bpy, mesh, **kwargs):
    return meshioplusplus.to_blender(mesh, **kwargs)


def test_blender_validate_reports_nothing_to_fix():
    """Blender's own topology checker as the oracle.

    ``validate()`` returns True when it *had* to repair something, so False is
    a genuine structural proof that the loops and polygons are well-formed.
    """
    bpy = pytest.importorskip("bpy")
    for factory in (_quad_mesh, _two_tets, _mixed_surface):
        me = _new_mesh(bpy, factory(), validate=False)
        assert me.validate(verbose=False) is False, factory.__name__


def test_blender_quad_is_not_triangulated():
    bpy = pytest.importorskip("bpy")
    me = _new_mesh(bpy, _quad_mesh())
    assert len(me.polygons) == 1
    assert me.polygons[0].loop_total == 4


def test_blender_point_data_reaches_point_attributes():
    bpy = pytest.importorskip("bpy")
    me = _new_mesh(bpy, _mixed_surface())
    layer = me.attributes["T"]
    assert layer.domain == "POINT"
    assert layer.data_type == "FLOAT"
    values = np.empty(len(me.vertices), dtype=np.float32)
    layer.data.foreach_get("value", values)
    assert values.tolist() == list(range(9))


def test_blender_cell_data_reaches_face_attributes_with_the_owning_cell_value():
    bpy = pytest.importorskip("bpy")
    me = _new_mesh(bpy, _two_tets())
    layer = me.attributes["mat"]
    assert layer.domain == "FACE"
    values = np.empty(len(me.polygons), dtype=np.int32)
    layer.data.foreach_get("value", values)
    assert sorted(values.tolist()) == [10, 10, 10, 20, 20, 20]


def test_blender_regions_and_field_data_survive():
    bpy = pytest.importorskip("bpy")
    me = _new_mesh(bpy, _mixed_surface())
    assert me.attributes[_interop.REGION_PREFIX + "wall"].domain == "FACE"
    assert me.attributes[_interop.REGION_PREFIX + "inlet"].domain == "POINT"
    # Blender's ID-property arrays come back as a tuple, not a list.
    assert list(me[FIELD_PREFIX + "step"]) == [7]


def test_blender_round_trip():
    bpy = pytest.importorskip("bpy")
    mesh = _mixed_surface()
    me = _new_mesh(bpy, mesh)
    back = meshioplusplus.from_blender(me)
    assert [(b.type, len(b)) for b in back.cells] == [
        ("triangle", 1),
        ("quad", 1),
        ("polygon", 1),
    ]
    assert np.allclose(back.points, mesh.points, atol=1e-6)
    assert np.allclose(back.point_data["T"], mesh.point_data["T"])
    by_key = {(r.kind, r.name): r for r in back.regions}
    assert by_key[("cell", "wall")].entries.tolist() == [1, 2]
    assert by_key[("cell", "wall")].tag == 9
