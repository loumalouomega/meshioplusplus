"""Tests for the in-memory interoperability layer.

The tests are in two halves, deliberately.

The **pure payload** tests run in the default CI matrix, where none of PyVista,
trimesh or pyarrow is installed. That is the whole point of ``_interop``'s
``_viewer``-style split: the subtle parts — block-major indexing, what is shared
and what is copied, the routing through existing operations — are exercised with
no optional dependency present at all.

The **gated** tests use ``pytest.importorskip`` (the HDF5/netCDF precedent) and
assert the things that can only be checked against the real library: that the
round-trips are lossless, that ``np.shares_memory`` really holds where sharing
is claimed, that ``zero_copy_only=True`` raises rather than copying silently,
and that a returned wrapper outlives the mesh it was built from.
"""

from __future__ import annotations

import gc
import json

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _interop
from meshioplusplus._mesh import Mesh
from meshioplusplus._regions import Region, block_bases


# --------------------------------------------------------------------------- #
# fixtures                                                                    #
# --------------------------------------------------------------------------- #
def _mixed_mesh():
    """triangle + quad + tetra, with point/cell data and all three region kinds."""
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    cells = [
        ("triangle", np.array([[0, 1, 2]], dtype=np.int64)),
        ("quad", np.array([[0, 1, 2, 3]], dtype=np.int64)),
        ("tetra", np.array([[0, 1, 2, 4], [0, 2, 3, 4]], dtype=np.int64)),
    ]
    return Mesh(
        points,
        cells,
        point_data={
            "T": np.arange(5, dtype=np.float64),
            "v": np.arange(15, dtype=np.float64).reshape(5, 3),
        },
        cell_data={
            "mat": [
                np.array([10], dtype=np.int64),
                np.array([20], dtype=np.int64),
                np.array([30, 31], dtype=np.int64),
            ]
        },
        regions=[
            Region("inlet", "point", [0, 1], dim=0, tag=5),
            Region("steel", "cell", [1, 3], dim=3, tag=7),
            Region("wall", "side", [[0, 1], [2, 0]], dim=2, tag=9),
        ],
    )


def _triangle_mesh():
    """A single-block, int64/float64 triangle mesh — the shareable-everything case."""
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    return Mesh(
        points,
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64))],
        point_data={
            "T": np.arange(4, dtype=np.float64),
            "v": np.arange(12, dtype=np.float64).reshape(4, 3),
        },
        cell_data={"mat": [np.array([1, 2], dtype=np.int64)]},
        regions=[
            Region("inlet", "point", [0, 1], dim=0, tag=5),
            Region("skin", "cell", [1], dim=2, tag=9),
        ],
    )


# --------------------------------------------------------------------------- #
# pure payload: VTK                                                           #
# --------------------------------------------------------------------------- #
def test_vtk_payload_mixed_types():
    mesh = _mixed_mesh()
    payload = _interop._to_vtk_payload(mesh)

    # VTK type ids: triangle=5, quad=9, tetra=10.
    assert payload.celltypes.tolist() == [5, 9, 10, 10]
    assert payload.celltypes.dtype == np.uint8
    # offsets are the VTK 9 layout: num_cells + 1 entries, starting at 0.
    assert payload.offsets.tolist() == [0, 3, 7, 11, 15]
    assert payload.connectivity.tolist() == [
        0, 1, 2,
        0, 1, 2, 3,
        0, 1, 2, 4,
        0, 2, 3, 4,
    ]  # fmt: skip
    assert payload.num_cells == 4
    assert payload.dropped == []


def test_vtk_payload_cell_data_is_block_major():
    """cell_data concatenation must follow the global block-major cell index."""
    mesh = _mixed_mesh()
    payload = _interop._to_vtk_payload(mesh)

    assert payload.cell_data["mat"].tolist() == [10, 20, 30, 31]
    # ... and that is exactly the numbering block_bases defines.
    bases = block_bases(mesh.cells)
    assert bases.tolist() == [0, 1, 2, 4]
    assert payload.kept_global.tolist() == [0, 1, 2, 3]


def test_vtk_payload_drops_unsupported_blocks():
    """A polyhedron block is dropped by name, and the global index skips it."""
    mesh = Mesh(
        np.zeros((6, 3)),
        [
            ("triangle", np.array([[0, 1, 2]], dtype=np.int64)),
            ("polyhedron4", [[[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]]),
            ("tetra", np.array([[0, 1, 2, 4], [1, 2, 3, 5]], dtype=np.int64)),
        ],
        cell_data={
            "mat": [
                np.array([1], dtype=np.int64),
                np.array([2], dtype=np.int64),
                np.array([3, 4], dtype=np.int64),
            ]
        },
    )
    payload = _interop._to_vtk_payload(mesh)

    assert [t for _, t in payload.dropped] == ["polyhedron4"]
    assert any("polyhedron4" in note for note in payload.notes)
    assert payload.celltypes.tolist() == [5, 10, 10]
    # The dropped block's cell data must not leak in, and the global index must
    # advance over it -- block 2's cells are global 2 and 3, not 1 and 2.
    assert payload.cell_data["mat"].tolist() == [1, 3, 4]
    assert payload.kept_global.tolist() == [0, 2, 3]


def test_vtk_payload_cell_region_entry_in_dropped_block_is_noted():
    mesh = Mesh(
        np.zeros((4, 3)),
        [
            ("triangle", np.array([[0, 1, 2]], dtype=np.int64)),
            ("polyhedron4", [[[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]]),
        ],
        regions=[Region("both", "cell", [0, 1])],
    )
    payload = _interop._to_vtk_payload(mesh)

    assert payload.region_masks[("cell", "both")].tolist() == [1]
    assert any("lie in dropped cell blocks" in note for note in payload.notes)


def test_vtk_payload_region_masks_and_metadata():
    mesh = _mixed_mesh()
    payload = _interop._to_vtk_payload(mesh)

    assert payload.region_masks[("point", "inlet")].tolist() == [1, 1, 0, 0, 0]
    # global cells 1 and 3 out of four.
    assert payload.region_masks[("cell", "steel")].tolist() == [0, 1, 0, 1]
    assert ("side", "wall") not in payload.region_masks
    assert payload.region_meta == [
        {"name": "inlet", "kind": "point", "dim": 0, "tag": 5},
        {"name": "steel", "kind": "cell", "dim": 3, "tag": 7},
    ]


def test_vtk_payload_side_region_dropped_with_a_note():
    mesh = _mixed_mesh()
    payload = _interop._to_vtk_payload(mesh)
    assert any(
        "side region 'wall' was dropped" in note for note in payload.notes
    ), payload.notes


def test_vtk_payload_shares_single_block_connectivity():
    mesh = _triangle_mesh()
    payload = _interop._to_vtk_payload(mesh)

    assert payload.shared["connectivity"] is True
    assert payload.shared["points"] is True
    assert payload.shared["point_data:T"] is True
    assert np.shares_memory(payload.connectivity, mesh.cells[0].data)
    assert np.shares_memory(payload.points, mesh.points)


def test_vtk_payload_derived_arrays_are_not_notes():
    """Deriving offsets/celltypes is not a loss, so it must not warn."""
    payload = _interop._to_vtk_payload(_triangle_mesh())
    assert payload.notes == []
    assert any("offsets/celltypes" in d for d in payload.derived)


def test_vtk_payload_wedge_is_reordered():
    """meshio's wedge node order differs from VTK's, which forces a copy."""
    mesh = Mesh(
        np.zeros((6, 3)),
        [("wedge", np.array([[0, 1, 2, 3, 4, 5]], dtype=np.int64))],
    )
    payload = _interop._to_vtk_payload(mesh)
    assert payload.connectivity.tolist() == [0, 2, 1, 3, 5, 4]
    assert payload.shared["connectivity"] is False


def test_vtk_payload_pads_2d_points():
    mesh = Mesh(
        np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        [("triangle", np.array([[0, 1, 2]], dtype=np.int64))],
    )
    payload = _interop._to_vtk_payload(mesh)
    assert payload.points.shape == (3, 3)
    assert payload.points[:, 2].tolist() == [0.0, 0.0, 0.0]
    assert payload.shared["points"] is False
    assert any("2D" in note for note in payload.notes)


@pytest.mark.parametrize(
    "mesh_factory, expected",
    [
        (
            lambda: Mesh(
                np.zeros((3, 3)),
                [("triangle", np.array([[0, 1, 2]], dtype=np.int32))],
            ),
            "int32",
        ),
        (
            lambda: Mesh(
                np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
                [("triangle", np.array([[0, 1, 2]], dtype=np.int64))],
            ),
            "2D",
        ),
    ],
)
def test_vtk_payload_zero_copy_only_raises(mesh_factory, expected):
    with pytest.raises(_interop.InteropError) as excinfo:
        _interop._to_vtk_payload(mesh_factory(), zero_copy_only=True)
    assert expected in str(excinfo.value)


def test_vtk_payload_zero_copy_only_allows_derived_arrays():
    """A mixed mesh must still convert: offsets/concatenation are exempt."""
    payload = _interop._to_vtk_payload(_mixed_mesh(), zero_copy_only=True)
    assert payload.num_cells == 4
    assert payload.shared["connectivity"] is False


# --------------------------------------------------------------------------- #
# pure payload: triangles                                                     #
# --------------------------------------------------------------------------- #
def test_triangles_payload_shares_a_pure_triangle_mesh():
    mesh = _triangle_mesh()
    payload = _interop._to_triangles_payload(mesh)

    assert payload.ops == []
    assert payload.shared["vertices"] is True
    assert payload.shared["faces"] is True
    assert np.shares_memory(payload.faces, mesh.cells[0].data)
    assert payload.region_masks[("point", "inlet")].tolist() == [1, 1, 0, 0]
    assert payload.region_masks[("cell", "skin")].tolist() == [0, 1]


def test_triangles_payload_routes_a_volume_through_extract_surface():
    mesh = Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
    )
    payload = _interop._to_triangles_payload(mesh)

    assert payload.ops == ["extract_surface"]
    assert payload.faces.shape == (4, 3)
    assert any("boundary was extracted" in note for note in payload.notes)


def test_triangles_payload_simplexifies_quads():
    mesh = Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]),
        [("quad", np.array([[0, 1, 2, 3]], dtype=np.int64))],
    )
    payload = _interop._to_triangles_payload(mesh)

    assert payload.ops == ["simplexify"]
    assert payload.faces.shape == (2, 3)


def test_triangles_payload_linearizes_higher_order():
    mesh = Mesh(
        np.zeros((6, 3)),
        [("triangle6", np.array([[0, 1, 2, 3, 4, 5]], dtype=np.int64))],
    )
    payload = _interop._to_triangles_payload(mesh)

    assert payload.ops == ["linearize"]
    assert payload.faces.shape == (1, 3)


def test_triangles_payload_raises_without_triangles():
    mesh = Mesh(np.zeros((2, 3)), [("line", np.array([[0, 1]], dtype=np.int64))])
    with pytest.raises(ValueError, match="no triangles"):
        _interop._to_triangles_payload(mesh)


# --------------------------------------------------------------------------- #
# pure payload: table                                                         #
# --------------------------------------------------------------------------- #
def test_table_payload_point_columns():
    payload = _interop._to_table_payload(_mixed_mesh(), "point")

    assert payload.names == ["x", "y", "z", "T", "v"]
    components = {c.name: c.components for c in payload.columns}
    # Multi-component arrays keep their shape (a fixed_size_list), rather than
    # being flattened into T_0/T_1/T_2 columns.
    assert components == {"x": 1, "y": 1, "z": 1, "T": 1, "v": 3}
    assert payload.num_rows == 5
    # x/y/z are structural: they are slices of a row-major points array.
    assert [c.derived for c in payload.columns] == [True, True, True, False, False]


def test_table_payload_cell_columns_use_the_global_index():
    mesh = _mixed_mesh()
    payload = _interop._to_table_payload(mesh, "cell")

    assert payload.names == ["block", "cell_type", "cell", "mat"]
    by_name = {c.name: c.values for c in payload.columns}
    assert by_name["block"].tolist() == [0, 1, 2, 2]
    assert list(by_name["cell_type"]) == ["triangle", "quad", "tetra", "tetra"]
    # The `cell` column is the global block-major index -- the numbering regions
    # and partition_labels use.
    assert by_name["cell"].tolist() == [0, 1, 2, 3]
    assert by_name["mat"].tolist() == [10, 20, 30, 31]


def test_table_payload_metadata_is_self_describing():
    payload = _interop._to_table_payload(_mixed_mesh(), "cell")
    meta = payload.metadata

    assert meta["meshioplusplus:location"] == "cell"
    assert meta["meshioplusplus:num_points"] == "5"
    assert meta["meshioplusplus:num_cells"] == "4"
    assert json.loads(meta["meshioplusplus:cell_types"]) == [
        ["triangle", 1],
        ["quad", 1],
        ["tetra", 2],
    ]
    names = {r["name"] for r in json.loads(meta[_interop.REGION_META_KEY])}
    assert names == {"inlet", "steel", "wall"}


def test_table_payload_unknown_location():
    with pytest.raises(ValueError, match="unknown location"):
        _interop._to_table_payload(_mixed_mesh(), "field")


def test_table_payload_unknown_location_names_op():
    with pytest.raises(ValueError, match="to_pandas"):
        _interop._to_table_payload(_mixed_mesh(), "field", op="to_pandas")


# --------------------------------------------------------------------------- #
# pure payload: frame columns (the pandas suffix expansion)                   #
# --------------------------------------------------------------------------- #
def test_frame_columns_expands_components():
    payload = _interop._to_table_payload(_mixed_mesh(), "point")
    columns, groups = _interop._frame_columns(payload)

    assert [c.name for c in columns] == ["x", "y", "z", "T", "v_0", "v_1", "v_2"]
    assert all(c.components == 1 for c in columns)
    assert groups == {"v": ["v_0", "v_1", "v_2"]}
    by_name = {c.name: c for c in columns}
    v = _mixed_mesh().point_data["v"]
    for i in range(3):
        assert by_name[f"v_{i}"].values.tolist() == v[:, i].tolist()
        # The expansion inherits the source's derived flag: v is real mesh
        # data, so its components are copy-noted, unlike x/y/z.
        assert by_name[f"v_{i}"].derived is False


def test_frame_columns_scalar_passthrough_is_identity():
    payload = _interop._to_table_payload(_mixed_mesh(), "point")
    columns, _ = _interop._frame_columns(payload)

    by_name = {c.name: c for c in columns}
    for original in payload.columns:
        if original.components == 1:
            assert by_name[original.name] is original


def test_frame_columns_collision_drops_with_note():
    mesh = _mixed_mesh()
    # A scalar whose name is exactly what expanding `v` produces. Sorted data
    # names put `v` before `v_1`, so the expansion claims the name first.
    mesh.point_data["v_1"] = np.arange(5, dtype=np.float64)
    payload = _interop._to_table_payload(mesh, "point")
    columns, groups = _interop._frame_columns(payload)

    names = [c.name for c in columns]
    assert names == ["x", "y", "z", "T", "v_0", "v_1", "v_2"]
    assert groups == {"v": ["v_0", "v_1", "v_2"]}
    # The v_1 column that survived is the vector component, not the scalar.
    by_name = {c.name: c for c in columns}
    assert by_name["v_1"].values.tolist() == mesh.point_data["v"][:, 1].tolist()
    assert any(
        "column 'v_1' was dropped" in note for note in payload.notes
    ), payload.notes


# --------------------------------------------------------------------------- #
# predicates, named errors, Phase 2                                           #
# --------------------------------------------------------------------------- #
def test_missing_dependency_error_names_the_extra():
    with pytest.raises(ImportError) as excinfo:
        _interop._require("meshioplusplus_no_such_module", "pyvista", "to_pyvista")
    message = str(excinfo.value)
    assert "to_pyvista" in message
    assert "pip install meshioplusplus[pyvista]" in message


def test_predicates_are_booleans():
    for predicate in (
        _interop.has_pyvista,
        _interop.has_trimesh,
        _interop.has_arrow,
        _interop.has_pandas,
        _interop.has_polars,
        _interop.has_open3d,
        _interop.has_dolfinx,
    ):
        assert isinstance(predicate(), bool)


def test_phase_two_targets_are_not_implemented():
    assert _interop.has_open3d() is False
    assert _interop.has_dolfinx() is False
    for fn in (_interop.to_open3d, _interop.to_dolfinx):
        with pytest.raises(NotImplementedError, match="Phase 2"):
            fn(_triangle_mesh())


def test_public_api_is_exported():
    for name in (
        "to_pyvista",
        "from_pyvista",
        "to_trimesh",
        "from_trimesh",
        "to_arrow",
        "from_arrow",
        "write_parquet",
        "read_parquet",
        "to_pandas",
        "to_polars",
        "has_pyvista",
        "has_trimesh",
        "has_arrow",
        "has_pandas",
        "has_polars",
    ):
        assert name in meshioplusplus.__all__
        assert hasattr(meshioplusplus, name)


# --------------------------------------------------------------------------- #
# gated: PyVista                                                              #
# --------------------------------------------------------------------------- #
def test_pyvista_round_trip():
    pytest.importorskip("pyvista")
    mesh = _mixed_mesh()
    grid = meshioplusplus.to_pyvista(mesh)

    assert grid.n_points == 5
    assert grid.n_cells == 4
    back = meshioplusplus.from_pyvista(grid)

    assert [cb.type for cb in back.cells] == ["triangle", "quad", "tetra"]
    assert np.allclose(back.points, mesh.points)
    assert np.allclose(back.point_data["T"], mesh.point_data["T"])
    assert np.array_equal(
        np.concatenate(back.cell_data["mat"]), np.array([10, 20, 30, 31])
    )
    assert meshioplusplus.meshes_equal(mesh, back)


def test_pyvista_round_trip_preserves_regions_with_dim_and_tag():
    pytest.importorskip("pyvista")
    mesh = _mixed_mesh()
    back = meshioplusplus.from_pyvista(meshioplusplus.to_pyvista(mesh))

    by_key = {(r.kind, r.name): r for r in back.regions}
    assert set(by_key) == {("point", "inlet"), ("cell", "steel")}
    assert by_key[("point", "inlet")] == Region("inlet", "point", [0, 1], dim=0, tag=5)
    # The mask alone could not carry dim/tag; the field_data sidecar does.
    assert by_key[("cell", "steel")] == Region("steel", "cell", [1, 3], dim=3, tag=7)


def test_pyvista_side_region_is_dropped():
    pytest.importorskip("pyvista")
    back = meshioplusplus.from_pyvista(meshioplusplus.to_pyvista(_mixed_mesh()))
    assert all(r.kind != "side" for r in back.regions)


def test_pyvista_regions_without_the_sidecar_lose_dim_and_tag():
    pytest.importorskip("pyvista")
    grid = meshioplusplus.to_pyvista(_mixed_mesh())
    del grid.field_data[_interop.REGION_META_KEY]

    back = meshioplusplus.from_pyvista(grid)
    by_key = {(r.kind, r.name): r for r in back.regions}
    assert by_key[("cell", "steel")].entries.tolist() == [1, 3]
    assert by_key[("cell", "steel")].dim == -1
    assert by_key[("cell", "steel")].tag == -1


def test_pyvista_shares_buffers():
    pytest.importorskip("pyvista")
    mesh = _triangle_mesh()
    grid = meshioplusplus.to_pyvista(mesh)

    assert np.shares_memory(np.asarray(grid.points), mesh.points)
    assert np.shares_memory(np.asarray(grid.point_data["T"]), mesh.point_data["T"])
    assert np.shares_memory(np.asarray(grid.point_data["v"]), mesh.point_data["v"])
    assert np.shares_memory(np.asarray(grid.cell_connectivity), mesh.cells[0].data)


def test_pyvista_does_not_share_when_a_conversion_is_needed():
    pytest.importorskip("pyvista")
    mesh = Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2]], dtype=np.int32))],
    )
    grid = meshioplusplus.to_pyvista(mesh)
    assert not np.shares_memory(np.asarray(grid.cell_connectivity), mesh.cells[0].data)


def test_pyvista_zero_copy_only_raises_on_dtype_mismatch():
    pytest.importorskip("pyvista")
    mesh = Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2]], dtype=np.int32))],
    )
    with pytest.raises(_interop.InteropError, match="int32"):
        meshioplusplus.to_pyvista(mesh, zero_copy_only=True)


def test_pyvista_outlives_the_source_mesh():
    """The readers hand back numpy owned by the C++ core -- a shared buffer must
    be kept alive by the wrapper, or this is a use-after-free rather than a
    wrong answer."""
    pytest.importorskip("pyvista")

    def build():
        mesh = Mesh(
            np.arange(12, dtype=np.float64).reshape(4, 3),
            [("triangle", np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64))],
            point_data={"T": np.arange(4, dtype=np.float64) * 3.0},
        )
        return meshioplusplus.to_pyvista(mesh)

    grid = build()
    gc.collect()
    for _ in range(3):
        churn = np.random.rand(1_000_000)
        del churn
    gc.collect()

    assert np.array_equal(
        np.asarray(grid.points), np.arange(12, dtype=np.float64).reshape(4, 3)
    )
    assert np.array_equal(np.asarray(grid.point_data["T"]), [0.0, 3.0, 6.0, 9.0])
    assert np.asarray(grid.cell_connectivity).tolist() == [0, 1, 2, 0, 2, 3]


# --------------------------------------------------------------------------- #
# gated: trimesh                                                              #
# --------------------------------------------------------------------------- #
def test_trimesh_round_trip():
    pytest.importorskip("trimesh")
    mesh = _triangle_mesh()
    tm = meshioplusplus.to_trimesh(mesh)

    assert tm.vertices.shape == (4, 3)
    assert tm.faces.shape == (2, 3)
    back = meshioplusplus.from_trimesh(tm)
    assert np.allclose(back.points, mesh.points)
    assert np.array_equal(back.cells[0].data, mesh.cells[0].data)
    assert np.allclose(back.point_data["T"], mesh.point_data["T"])


def test_trimesh_shares_buffers():
    pytest.importorskip("trimesh")
    mesh = _triangle_mesh()
    tm = meshioplusplus.to_trimesh(mesh)

    assert np.shares_memory(np.asarray(tm.vertices), mesh.points)
    assert np.shares_memory(np.asarray(tm.faces), mesh.cells[0].data)


def test_trimesh_does_not_merge_vertices():
    """process=True would silently renumber, breaking every index-based array."""
    pytest.importorskip("trimesh")
    mesh = Mesh(
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 0.0],
                # a deliberate duplicate of vertex 0
                [0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ]
        ),
        [("triangle", np.array([[0, 1, 2], [3, 4, 5]], dtype=np.int64))],
    )
    tm = meshioplusplus.to_trimesh(mesh)
    assert len(tm.vertices) == 6


def test_trimesh_carries_regions_for_a_triangle_mesh():
    pytest.importorskip("trimesh")
    back = meshioplusplus.from_trimesh(meshioplusplus.to_trimesh(_triangle_mesh()))
    by_key = {(r.kind, r.name): r for r in back.regions}
    assert by_key[("point", "inlet")] == Region("inlet", "point", [0, 1], dim=0, tag=5)
    assert by_key[("cell", "skin")] == Region("skin", "cell", [1], dim=2, tag=9)


def test_trimesh_volume_input_loses_regions():
    """extract_surface drops regions; to_trimesh gets no second policy."""
    pytest.importorskip("trimesh")
    mesh = Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
        [("tetra", np.array([[0, 1, 2, 3]], dtype=np.int64))],
        regions=[Region("blob", "cell", [0])],
    )
    tm = meshioplusplus.to_trimesh(mesh)
    assert not any(k.startswith(_interop.REGION_PREFIX) for k in tm.face_attributes)


def test_trimesh_zero_copy_only_raises_on_dtype_mismatch():
    pytest.importorskip("trimesh")
    mesh = Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float32),
        [("triangle", np.array([[0, 1, 2]], dtype=np.int64))],
    )
    with pytest.raises(_interop.InteropError, match="float32"):
        meshioplusplus.to_trimesh(mesh, zero_copy_only=True)


def test_trimesh_outlives_the_source_mesh():
    pytest.importorskip("trimesh")

    def build():
        mesh = Mesh(
            np.arange(12, dtype=np.float64).reshape(4, 3),
            [("triangle", np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64))],
        )
        return meshioplusplus.to_trimesh(mesh)

    tm = build()
    gc.collect()
    for _ in range(3):
        churn = np.random.rand(1_000_000)
        del churn
    gc.collect()
    assert np.array_equal(
        np.asarray(tm.vertices), np.arange(12, dtype=np.float64).reshape(4, 3)
    )


# --------------------------------------------------------------------------- #
# gated: Arrow / Parquet                                                      #
# --------------------------------------------------------------------------- #
def test_arrow_point_table():
    pa = pytest.importorskip("pyarrow")
    mesh = _mixed_mesh()
    table = meshioplusplus.to_arrow(mesh, "point")

    assert table.schema.names == ["x", "y", "z", "T", "v"]
    assert table.num_rows == 5
    # A multi-component array keeps its shape as a fixed_size_list.
    assert pa.types.is_fixed_size_list(table.schema.field("v").type)
    assert table.schema.field("v").type.list_size == 3


def test_arrow_cell_table():
    pytest.importorskip("pyarrow")
    table = meshioplusplus.to_arrow(_mixed_mesh(), "cell")

    assert table.schema.names == ["block", "cell_type", "cell", "mat"]
    assert table.column("cell").to_pylist() == [0, 1, 2, 3]
    assert table.column("cell_type").to_pylist() == [
        "triangle",
        "quad",
        "tetra",
        "tetra",
    ]


def test_arrow_schema_metadata():
    pytest.importorskip("pyarrow")
    table = meshioplusplus.to_arrow(_mixed_mesh(), "cell")
    meta = {k.decode(): v.decode() for k, v in table.schema.metadata.items()}

    assert meta["meshioplusplus:location"] == "cell"
    assert meta["meshioplusplus:num_cells"] == "4"
    assert "meshioplusplus:version" in meta
    assert json.loads(meta["meshioplusplus:cell_types"])[0] == ["triangle", 1]


def test_arrow_shares_numeric_buffers():
    pytest.importorskip("pyarrow")
    mesh = _triangle_mesh()
    table = meshioplusplus.to_arrow(mesh, "point")

    scalar = _interop._single_array(table.column("T")).to_numpy(zero_copy_only=True)
    assert np.shares_memory(scalar, mesh.point_data["T"])
    listed = (
        _interop._single_array(table.column("v"))
        .flatten()
        .to_numpy(zero_copy_only=True)
    )
    assert np.shares_memory(listed, mesh.point_data["v"])


def test_arrow_round_trip_restores_component_shape():
    pytest.importorskip("pyarrow")
    mesh = _triangle_mesh()
    arrays = meshioplusplus.from_arrow(meshioplusplus.to_arrow(mesh, "point"))

    assert arrays["v"].shape == (4, 3)
    assert np.allclose(arrays["v"], mesh.point_data["v"])
    assert np.allclose(arrays["T"], mesh.point_data["T"])


def test_parquet_round_trip(tmp_path):
    pytest.importorskip("pyarrow")
    mesh = _triangle_mesh()
    path = tmp_path / "data.parquet"
    meshioplusplus.write_parquet(mesh, path, location="point")

    arrays, metadata = meshioplusplus.read_parquet(path, return_metadata=True)
    assert arrays["v"].shape == (4, 3)
    assert arrays["v"].dtype == np.float64
    assert np.allclose(arrays["v"], mesh.point_data["v"])
    assert np.allclose(arrays["x"], mesh.points[:, 0])
    assert metadata["meshioplusplus:location"] == "point"


def test_parquet_cell_round_trip_preserves_dtypes(tmp_path):
    pytest.importorskip("pyarrow")
    mesh = _mixed_mesh()
    path = tmp_path / "cells.parquet"
    meshioplusplus.write_parquet(mesh, path, location="cell")

    arrays = meshioplusplus.read_parquet(path)
    assert arrays["mat"].dtype == np.int64
    assert arrays["mat"].tolist() == [10, 20, 30, 31]
    assert arrays["cell"].tolist() == [0, 1, 2, 3]


def test_parquet_is_not_a_mesh_format():
    """It exports data for analytics; it must not look like a mesh writer."""
    from meshioplusplus._helpers import _writer_map, extension_to_filetypes

    assert "parquet" not in _writer_map
    assert ".parquet" not in extension_to_filetypes


def test_cli_data_export(tmp_path):
    pytest.importorskip("pyarrow")
    from meshioplusplus._cli import main

    mesh = _triangle_mesh()
    infile = tmp_path / "in.vtu"
    outfile = tmp_path / "out.parquet"
    meshioplusplus.write(infile, mesh)

    assert (
        main(["data", "export", str(infile), str(outfile), "--location", "cell"]) == 0
    )
    arrays = meshioplusplus.read_parquet(outfile)
    assert arrays["mat"].tolist() == [1, 2]


# --------------------------------------------------------------------------- #
# gated: pandas                                                               #
# --------------------------------------------------------------------------- #
def test_pandas_point_frame():
    pytest.importorskip("pandas")
    mesh = _mixed_mesh()
    df = meshioplusplus.to_pandas(mesh, "point")

    assert list(df.columns) == ["x", "y", "z", "T", "v_0", "v_1", "v_2"]
    assert len(df) == 5
    assert df["T"].dtype == np.float64
    assert df["v_1"].tolist() == mesh.point_data["v"][:, 1].tolist()


def test_pandas_cell_frame():
    pytest.importorskip("pandas")
    df = meshioplusplus.to_pandas(_mixed_mesh(), "cell")

    assert list(df.columns) == ["block", "cell_type", "cell", "mat"]
    assert df["cell_type"].tolist() == ["triangle", "quad", "tetra", "tetra"]
    assert df["cell"].tolist() == [0, 1, 2, 3]
    assert df["mat"].dtype == np.int64
    assert df["mat"].tolist() == [10, 20, 30, 31]


def test_pandas_attrs_metadata():
    pytest.importorskip("pandas")
    df = meshioplusplus.to_pandas(_mixed_mesh(), "point")

    assert df.attrs["meshioplusplus:location"] == "point"
    assert df.attrs["meshioplusplus:num_points"] == "5"
    assert json.loads(df.attrs[_interop.COMPONENTS_META_KEY]) == {
        "v": ["v_0", "v_1", "v_2"]
    }
    # No expansion -> the key is still present, recording that fact.
    cells = meshioplusplus.to_pandas(_mixed_mesh(), "cell")
    assert json.loads(cells.attrs[_interop.COMPONENTS_META_KEY]) == {}


def test_pandas_expansion_warns_and_zero_copy_raises(monkeypatch):
    pytest.importorskip("pandas")
    recorded = []
    monkeypatch.setattr(_interop, "warn", lambda s, highlight=True: recorded.append(s))
    mesh = _mixed_mesh()

    meshioplusplus.to_pandas(mesh, "point")
    assert any("column 'v'" in s for s in recorded), recorded

    with pytest.raises(_interop.InteropError, match="column 'v'"):
        meshioplusplus.to_pandas(mesh, "point", zero_copy_only=True)


def test_pandas_scalar_sharing_is_honest(monkeypatch):
    """A scalar column is either genuinely shared or its copy was warned."""
    pytest.importorskip("pandas")
    recorded = []
    monkeypatch.setattr(_interop, "warn", lambda s, highlight=True: recorded.append(s))
    mesh = _triangle_mesh()
    del mesh.point_data["v"]  # keep the frame all-scalar

    df = meshioplusplus.to_pandas(mesh, "point")
    shared = np.shares_memory(df["T"].to_numpy(), mesh.point_data["T"])
    warned = any("column 'T'" in s for s in recorded)
    assert shared != warned, (shared, recorded)


# --------------------------------------------------------------------------- #
# gated: polars                                                               #
# --------------------------------------------------------------------------- #
def test_polars_point_frame():
    pl = pytest.importorskip("polars")
    mesh = _mixed_mesh()
    df = meshioplusplus.to_polars(mesh, "point")

    assert df.columns == ["x", "y", "z", "T", "v"]
    assert df.height == 5
    # The (n, 3) shape survives structurally -- no suffix columns.
    assert df.schema["v"] == pl.Array(pl.Float64, 3)
    assert df["v"].to_numpy().shape == (5, 3)
    assert np.allclose(df["v"].to_numpy(), mesh.point_data["v"])


def test_polars_cell_frame():
    pl = pytest.importorskip("polars")
    df = meshioplusplus.to_polars(_mixed_mesh(), "cell")

    assert df.columns == ["block", "cell_type", "cell", "mat"]
    assert df.schema["cell_type"] == pl.String
    assert df["cell_type"].to_list() == ["triangle", "quad", "tetra", "tetra"]
    assert df.schema["mat"] == pl.Int64
    assert df["mat"].to_list() == [10, 20, 30, 31]
