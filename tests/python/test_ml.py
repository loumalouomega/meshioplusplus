"""Tests for the ML data-handling layer (edge_index / feature_matrix / datasets).

Two halves, the ``test_interop.py`` pattern: the graph and feature-matrix
tests are pure numpy and run in the default CI matrix with no optional library
installed; the dataset-export tests gate per backend with
``pytest.importorskip`` (pyarrow / zarr / h5py).
"""

from __future__ import annotations

import json

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import _ml
from meshioplusplus._mesh import Mesh
from meshioplusplus._regions import Region


# --------------------------------------------------------------------------- #
# fixtures                                                                    #
# --------------------------------------------------------------------------- #
def _mixed_mesh():
    """triangle + quad + tetra with point/cell data and point/cell regions."""
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


def _triangle_mesh(temperature=None):
    points = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    return Mesh(
        points,
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64))],
        point_data={
            "T": (
                np.arange(4, dtype=np.float64)
                if temperature is None
                else np.asarray(temperature, dtype=np.float64)
            )
        },
    )


# --------------------------------------------------------------------------- #
# edge_index: node graph                                                      #
# --------------------------------------------------------------------------- #
def test_edge_index_node_unique_pairs():
    ei = meshioplusplus.edge_index(_mixed_mesh(), undirected=False)

    # triangle 012, quad 0123, tetra 0124 and 0234 — nine unique edges.
    expected = [
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 4),
        (1, 2),
        (1, 4),
        (2, 3),
        (2, 4),
        (3, 4),
    ]
    assert ei.shape == (2, 9)
    assert ei.dtype == np.int64
    assert list(zip(ei[0].tolist(), ei[1].tolist())) == expected
    assert (ei[0] < ei[1]).all()


def test_edge_index_undirected_emits_both_directions():
    mesh = _mixed_mesh()
    ei = meshioplusplus.edge_index(mesh)  # undirected=True is the default

    assert ei.shape == (2, 18)
    pairs = set(zip(ei[0].tolist(), ei[1].tolist()))
    assert len(pairs) == 18
    for a, b in list(pairs):
        assert (b, a) in pairs
    # Sorted lexicographically by (source, target) — the recorded order.
    order = np.lexsort((ei[1], ei[0]))
    assert (order == np.arange(18)).all()
    assert ei.flags["C_CONTIGUOUS"]


def test_edge_index_is_deterministic():
    mesh = _mixed_mesh()
    first = meshioplusplus.edge_index(mesh)
    for _ in range(3):
        assert np.array_equal(meshioplusplus.edge_index(mesh), first)


def test_edge_index_ragged_polygon_ring():
    mesh = Mesh(
        np.zeros((5, 3)),
        [("polygon", [[0, 1, 2, 3, 4]])],
    )
    ei = meshioplusplus.edge_index(mesh, undirected=False)
    assert list(zip(ei[0].tolist(), ei[1].tolist())) == [
        (0, 1),
        (0, 4),
        (1, 2),
        (2, 3),
        (3, 4),
    ]


def test_edge_index_warns_on_unknown_topology(monkeypatch):
    recorded = []
    monkeypatch.setattr(_ml, "warn", lambda s, highlight=True: recorded.append(s))
    mesh = Mesh(
        np.zeros((10, 3)),
        [("tetra10", np.arange(10, dtype=np.int64).reshape(1, 10))],
    )
    ei = meshioplusplus.edge_index(mesh)
    assert ei.shape == (2, 0)
    assert any("tetra10" in s for s in recorded), recorded


def test_edge_index_unknown_kind():
    with pytest.raises(ValueError, match="unknown kind"):
        meshioplusplus.edge_index(_mixed_mesh(), kind="face")


# --------------------------------------------------------------------------- #
# edge_index: cell dual                                                       #
# --------------------------------------------------------------------------- #
def test_edge_index_cell_dual():
    ei = meshioplusplus.edge_index(_mixed_mesh(), kind="cell")

    # Only the two tetras (global cells 2 and 3) share a facet (face 0-2-4);
    # the triangle and quad sit below the dual dimension.
    assert ei.tolist() == [[2, 3], [3, 2]]


def test_edge_index_cell_dual_2d():
    ei = meshioplusplus.edge_index(_triangle_mesh(), kind="cell", undirected=False)
    # The two triangles share edge 0-2.
    assert ei.tolist() == [[0], [1]]


# --------------------------------------------------------------------------- #
# feature_matrix                                                              #
# --------------------------------------------------------------------------- #
def test_feature_matrix_point_contract():
    fm = meshioplusplus.feature_matrix(_mixed_mesh(), "point")

    assert fm.columns == ("x", "y", "z", "T", "v_0", "v_1", "v_2", "region:inlet")
    assert fm.matrix.shape == (5, 8)
    assert fm.matrix.dtype == np.float64
    assert fm.matrix.flags["C_CONTIGUOUS"]
    assert fm.matrix[:, 3].tolist() == [0.0, 1.0, 2.0, 3.0, 4.0]
    assert fm.matrix[:, 7].tolist() == [1.0, 1.0, 0.0, 0.0, 0.0]
    assert fm.schema["feature_schema_version"] == _ml.FEATURE_SCHEMA_VERSION
    assert fm.schema["columns"] == list(fm.columns)


def test_feature_matrix_cell_excludes_structural_columns():
    fm = meshioplusplus.feature_matrix(_mixed_mesh(), "cell")

    # block/cell_type/cell are identifiers, not features.
    assert fm.columns == ("mat", "region:steel")
    assert fm.matrix[:, 0].tolist() == [10.0, 20.0, 30.0, 31.0]
    assert fm.matrix[:, 1].tolist() == [0.0, 1.0, 0.0, 1.0]


def test_feature_matrix_fields_order_is_honoured():
    mesh = _mixed_mesh()
    fm = meshioplusplus.feature_matrix(
        mesh, "point", fields=["v", "T"], coords=False, regions=False
    )
    assert fm.columns == ("v_0", "v_1", "v_2", "T")

    with pytest.raises(ValueError, match="no point data array named"):
        meshioplusplus.feature_matrix(mesh, "point", fields=["nope"])


def test_feature_matrix_region_selection():
    mesh = _mixed_mesh()
    none = meshioplusplus.feature_matrix(mesh, "point", regions=False)
    assert "region:inlet" not in none.columns
    picked = meshioplusplus.feature_matrix(mesh, "point", regions=["inlet"])
    assert picked.columns[-1] == "region:inlet"


def test_feature_matrix_sources_describe_every_column():
    fm = meshioplusplus.feature_matrix(_mixed_mesh(), "point")
    by_name = {s["name"]: s for s in fm.schema["sources"]}
    assert set(by_name) == set(fm.columns)
    assert by_name["v_1"] == {
        "name": "v_1",
        "source": "v",
        "kind": "data",
        "component": 1,
    }
    assert by_name["x"]["kind"] == "coords"
    assert by_name["region:inlet"]["kind"] == "region"


# --------------------------------------------------------------------------- #
# write_dataset: validation (no optional dep needed)                          #
# --------------------------------------------------------------------------- #
def test_write_dataset_rejects_unknown_format(tmp_path):
    with pytest.raises(ValueError, match="unknown format"):
        meshioplusplus.write_dataset([], tmp_path / "d", format="csv")


def test_write_dataset_rejects_unknown_mesh_id_rule(tmp_path):
    with pytest.raises(ValueError, match="unknown mesh_id rule"):
        meshioplusplus.write_dataset([], tmp_path / "d", mesh_id="uuid")


def test_mesh_ids_stem_collision_is_named():
    entries = [
        {"path": "/a/mesh.vtu", "step": 0},
        {"path": "/b/mesh.vtu", "step": 0},
    ]
    with pytest.raises(ValueError, match="not unique"):
        _ml._mesh_ids(entries, "stem")
    # index ids never collide
    assert _ml._mesh_ids(entries, "index") == ["0000", "0001"]


def test_mesh_ids_multi_step_file_appends_the_step():
    entries = [
        {"path": "/a/series.xdmf", "step": 0},
        {"path": "/a/series.xdmf", "step": 1},
    ]
    assert _ml._mesh_ids(entries, "stem") == ["series_0", "series_1"]


# --------------------------------------------------------------------------- #
# gated: parquet dataset                                                      #
# --------------------------------------------------------------------------- #
def _write_inputs(tmp_path, count=2):
    paths = []
    for i in range(count):
        p = tmp_path / f"case_{i}.vtu"
        meshioplusplus.write(p, _triangle_mesh(temperature=[i, i, i, i]))
        paths.append(str(p))
    return paths


def test_parquet_dataset_layout_and_manifest(tmp_path):
    pq = pytest.importorskip("pyarrow.parquet")
    paths = _write_inputs(tmp_path)
    out = tmp_path / "dataset"

    manifest = meshioplusplus.write_dataset(paths, out, location="point")

    assert manifest["layout_version"] == _ml.DATASET_LAYOUT_VERSION
    assert manifest["columns"] == ["x", "y", "z", "T"]
    assert [e["mesh_id"] for e in manifest["entries"]] == ["case_0", "case_1"]
    on_disk = json.loads((out / _ml.DATASET_MANIFEST_NAME).read_text())
    assert on_disk == manifest

    table = pq.read_table(str(out / "mesh_id=case_1" / "data.parquet"))
    assert table.schema.names == ["x", "y", "z", "T", "mesh_id"]
    assert table.column("mesh_id").to_pylist() == ["case_1"] * 4
    assert table.column("T").to_pylist() == [1.0, 1.0, 1.0, 1.0]
    meta = {k.decode(): v.decode() for k, v in table.schema.metadata.items()}
    assert meta["meshioplusplus:mesh_id"] == "case_1"


def test_parquet_dataset_glob_source(tmp_path):
    pytest.importorskip("pyarrow")
    _write_inputs(tmp_path, count=3)
    manifest = meshioplusplus.write_dataset(
        str(tmp_path / "case_*.vtu"), tmp_path / "ds", location="point"
    )
    assert [e["mesh_id"] for e in manifest["entries"]] == [
        "case_0",
        "case_1",
        "case_2",
    ]


def test_dataset_schema_is_strict(tmp_path):
    pytest.importorskip("pyarrow")
    paths = _write_inputs(tmp_path)
    odd = _triangle_mesh()
    odd.point_data["extra"] = np.zeros(4)
    odd_path = tmp_path / "case_odd.vtu"
    meshioplusplus.write(odd_path, odd)

    with pytest.raises(ValueError, match="does not match the dataset schema"):
        meshioplusplus.write_dataset(
            paths + [str(odd_path)], tmp_path / "ds", location="point"
        )
    # A failed run leaves no manifest — a partial dataset must not look done.
    assert not (tmp_path / "ds" / _ml.DATASET_MANIFEST_NAME).exists()


def test_cli_data_export_dataset(tmp_path):
    pytest.importorskip("pyarrow")
    from meshioplusplus._cli import main

    paths = _write_inputs(tmp_path)
    out = tmp_path / "cli_ds"
    assert main(["data", "export-dataset", *paths, str(out)]) == 0
    assert (out / _ml.DATASET_MANIFEST_NAME).exists()


# --------------------------------------------------------------------------- #
# gated: zarr / hdf5 datasets                                                 #
# --------------------------------------------------------------------------- #
def test_zarr_dataset_layout(tmp_path):
    zarr = pytest.importorskip("zarr")
    paths = _write_inputs(tmp_path)
    out = tmp_path / "dataset.zarr"

    manifest = meshioplusplus.write_dataset(paths, out, format="zarr")

    root = zarr.open_group(str(out), mode="r")
    assert json.loads(root.attrs["meshioplusplus_dataset"]) == manifest
    group = root["case_0"]
    assert sorted(group.array_keys()) == ["T", "x", "y", "z"]
    assert np.asarray(group["T"]).tolist() == [0.0, 0.0, 0.0, 0.0]
    assert group.attrs["num_rows"] == 4


def test_hdf5_dataset_layout(tmp_path):
    h5py = pytest.importorskip("h5py")
    paths = _write_inputs(tmp_path)
    out = tmp_path / "dataset.h5"

    manifest = meshioplusplus.write_dataset(paths, out, format="hdf5")

    with h5py.File(out, "r") as handle:
        assert json.loads(handle.attrs["meshioplusplus_dataset"]) == manifest
        assert sorted(handle["case_1"].keys()) == ["T", "x", "y", "z"]
        assert handle["case_1"]["T"][...].tolist() == [1.0, 1.0, 1.0, 1.0]
        assert handle["case_1"].attrs["num_rows"] == 4


def test_group_layout_drops_non_numeric_with_warning(tmp_path, monkeypatch):
    h5py = pytest.importorskip("h5py")
    recorded = []
    monkeypatch.setattr(_ml, "warn", lambda s, highlight=True: recorded.append(s))
    paths = _write_inputs(tmp_path)
    out = tmp_path / "cells.h5"

    meshioplusplus.write_dataset(paths, out, format="hdf5", location="cell")

    # The object-dtype cell_type column cannot ride; block/cell (numeric) do.
    with h5py.File(out, "r") as handle:
        assert "cell_type" not in handle["case_0"]
        assert "block" in handle["case_0"]
        assert "cell_types" in handle["case_0"].attrs
    assert any("cell_type" in s for s in recorded), recorded


# --------------------------------------------------------------------------- #
# exports                                                                     #
# --------------------------------------------------------------------------- #
def test_public_api_is_exported():
    for name in (
        "edge_index",
        "feature_matrix",
        "FeatureMatrix",
        "write_dataset",
        "has_zarr",
    ):
        assert name in meshioplusplus.__all__
        assert hasattr(meshioplusplus, name)
    assert isinstance(meshioplusplus.has_zarr(), bool)
