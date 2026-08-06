"""Tests for the PhysicsNeMo adapter (`meshioplusplus.physicsnemo`).

Two halves, the `test_gpu.py` pattern: `graph_sample`, the stats accumulators
and the manifest iteration are pure numpy and run in the default CI matrix
with nothing optional installed; the framework-facing halves gate with
`pytest.importorskip` (torch / torch_geometric / physicsnemo) and are
exercised on a real GPU machine, not by public CI.
"""

from __future__ import annotations

import numpy as np
import pytest

import meshioplusplus
import meshioplusplus.physicsnemo as mpn
from meshioplusplus import DatasetManifest, _gpu
from meshioplusplus._mesh import Mesh
from meshioplusplus._regions import Region


def _mesh(offset=0.0):
    """Two triangles with scalar + vector point data, cell data, regions."""
    points = (
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
        + offset
    )
    return Mesh(
        points,
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64))],
        point_data={
            "T": np.arange(4, dtype=np.float64) + offset,
            "v": np.arange(8, dtype=np.float64).reshape(4, 2),
        },
        cell_data={"mat": [np.array([3.0, 4.0])]},
        regions=[Region("inlet", "point", [0, 1], dim=0, tag=1)],
    )


def _manifest(tmp_path, n=4):
    cases = tmp_path / "cases"
    cases.mkdir(exist_ok=True)
    manifest = DatasetManifest(base_dir=str(tmp_path))
    for i in range(n):
        meshioplusplus.write(str(cases / f"case_{i}.vtu"), _mesh(float(i)))
        manifest.add(
            f"cases/case_{i}.vtu", id=f"c{i}", split="train" if i < n - 1 else "test"
        )
    return manifest


# --------------------------------------------------------------------------- #
# graph_sample: shapes, dtypes, the contracts                                 #
# --------------------------------------------------------------------------- #
def test_graph_sample_keys_shapes_dtypes():
    s = mpn.graph_sample(_mesh(), fields=["T"], target_fields=["v"])
    assert set(s.arrays) == {"pos", "x", "y", "edge_index", "edge_attr"}
    n = 4
    assert s.arrays["pos"].shape == (n, 3)
    assert s.arrays["x"].shape == (n, 2)  # T + region:inlet one-hot
    assert s.arrays["y"].shape == (n, 2)  # v_0, v_1
    e = s.arrays["edge_index"].shape[1]
    assert s.arrays["edge_index"].shape == (2, e)
    assert s.arrays["edge_attr"].shape == (e, 4)  # dx, dy, dz, norm
    for name in ("pos", "x", "y", "edge_attr"):
        assert s.arrays[name].dtype == np.float32, name
    assert s.arrays["edge_index"].dtype == np.int64


def test_graph_sample_x_columns_are_the_feature_matrix_contract():
    mesh = _mesh()
    s = mpn.graph_sample(mesh, fields=["T"])
    fm = meshioplusplus.feature_matrix(
        mesh, "point", fields=["T"], coords=False, regions=True
    )
    assert s.x_columns == fm.columns == ("T", "region:inlet")
    assert s.schema["x_columns"] == list(fm.columns)
    assert s.schema["x_sources"] == fm.schema["sources"]
    # coordinates never enter x — the MGN convention
    assert "x" not in s.x_columns and "y" not in s.x_columns


def test_graph_sample_edge_attr_matches_the_convention_by_hand():
    s = mpn.graph_sample(_mesh(), float32=False)
    pos, ei, attr = s.arrays["pos"], s.arrays["edge_index"], s.arrays["edge_attr"]
    row, col = ei
    disp = pos[row] - pos[col]  # source minus destination — the sign matters
    expected = np.concatenate(
        (disp, np.linalg.norm(disp, axis=-1, keepdims=True)), axis=-1
    )
    np.testing.assert_array_equal(attr, expected)
    assert attr.dtype == np.float64  # float32=False keeps canonical dtypes


def test_graph_sample_flags():
    mesh = _mesh()
    s = mpn.graph_sample(mesh, edge_features=False, regions=False)
    assert "edge_attr" not in s.arrays and "y" not in s.arrays
    assert s.x_columns == ("T", "v_0", "v_1")  # every point array, no one-hots
    assert s.y_columns == () and s.schema["edge_features"] == []
    directed = mpn.graph_sample(mesh, undirected=False)
    assert (
        directed.arrays["edge_index"].shape[1] == s.arrays["edge_index"].shape[1] // 2
    )
    with pytest.raises(ValueError, match="kind must be"):
        mpn.graph_sample(mesh, kind="bogus")


def test_graph_sample_cell_kind_is_coherent():
    """kind='cell' switches pos/x to centroids and cell-located data, so the
    vertices of edge_index and the rows of pos/x are the same objects."""
    mesh = _mesh()
    s = mpn.graph_sample(mesh, kind="cell", fields=["mat"], float32=False)
    assert s.arrays["pos"].shape == (2, 3)
    np.testing.assert_allclose(s.arrays["pos"][0], mesh.points[[0, 1, 2]].mean(axis=0))
    assert s.x_columns == ("mat",)
    assert s.arrays["x"].shape == (2, 1)
    assert s.arrays["edge_index"].max() < 2


def test_graph_sample_schema_is_json_ready():
    import json

    s = mpn.graph_sample(_mesh(), fields=["T"], target_fields=["v"])
    text = json.dumps(s.schema)
    assert json.loads(text)["y_columns"] == ["v_0", "v_1"]
    assert s.schema["edge_features"] == ["dx", "dy", "dz", "norm"]


# --------------------------------------------------------------------------- #
# manifest walking: iter_samples, field_stats, edge_stats                     #
# --------------------------------------------------------------------------- #
def test_iter_samples_streams_and_filters_by_split(tmp_path):
    manifest = _manifest(tmp_path)
    seen = list(mpn.iter_samples(manifest, split="train", fields=["T"]))
    assert [entry_id for entry_id, _, _ in seen] == ["c0", "c1", "c2"]
    for _, _, sample in seen:
        assert isinstance(sample, mpn.GraphSample)
        # VTU carries no regions, so the re-read mesh has no one-hot column —
        # the format's documented behaviour, not the adapter's.
        assert sample.x_columns == ("T",)
    assert len(list(mpn.iter_samples(manifest))) == 4


def test_field_stats_matches_direct_concatenation(tmp_path):
    manifest = _manifest(tmp_path)
    stats = mpn.field_stats(manifest, split="train", fields=["T", "v"])
    all_t = np.concatenate([np.arange(4, dtype=np.float64) + i for i in range(3)])
    assert stats["T_mean"] == pytest.approx([all_t.mean()])
    assert stats["T_std"] == pytest.approx([all_t.std()])
    assert len(stats["v_mean"]) == len(stats["v_std"]) == 2
    v = np.arange(8, dtype=np.float64).reshape(4, 2)
    assert stats["v_mean"] == pytest.approx(np.tile(v, (3, 1)).mean(axis=0).tolist())
    # the keys are exactly the node_stats.json convention
    assert set(stats) == {"T_mean", "T_std", "v_mean", "v_std"}


def test_edge_stats_matches_direct_computation(tmp_path):
    manifest = _manifest(tmp_path, n=2)
    stats = mpn.edge_stats(manifest)
    attrs = []
    for entry in manifest:
        _, mesh = entry.time_series()[0]
        s = mpn.graph_sample(mesh, float32=False)
        attrs.append(s.arrays["edge_attr"])
    stacked = np.concatenate(attrs)
    assert stats["edge_mean"] == pytest.approx(stacked.mean(axis=0).tolist())
    assert stats["edge_std"] == pytest.approx(stacked.std(axis=0).tolist())


def test_manifest_argument_is_tri_modal(tmp_path):
    manifest = _manifest(tmp_path, n=2)
    path = tmp_path / "m.json"
    manifest.save(path)
    from_path = list(mpn.iter_samples(str(path), fields=["T"]))
    from_obj = list(mpn.iter_samples(manifest, fields=["T"]))
    assert [i for i, _, _ in from_path] == [i for i, _, _ in from_obj]


# --------------------------------------------------------------------------- #
# gating: predicates and the named install errors                             #
# --------------------------------------------------------------------------- #
def test_predicates_return_plain_bools():
    assert mpn.has_physicsnemo() in (True, False)
    assert mpn.has_torch_geometric() in (True, False)


def test_install_errors_name_the_command_and_no_extra(monkeypatch, tmp_path):
    manifest = _manifest(tmp_path, n=1)
    monkeypatch.setattr(_gpu, "_importable", lambda module: False)
    with pytest.raises(ImportError, match="pip install nvidia-physicsnemo") as e:
        mpn.make_reader(manifest)
    # There is deliberately no [physicsnemo] extra — the CuPy/torch precedent.
    assert "meshioplusplus[" not in str(e.value)
    assert "doc/physicsnemo.md" in str(e.value)
    with pytest.raises(ImportError, match="pip install torch_geometric") as e:
        mpn.make_dataset(manifest)
    assert "meshioplusplus[" not in str(e.value)


# --------------------------------------------------------------------------- #
# gated halves — dev GPU box only, never public CI                            #
# --------------------------------------------------------------------------- #
def test_pyg_dataset_yields_data_objects(tmp_path):
    pytest.importorskip("torch_geometric")
    import torch

    manifest = _manifest(tmp_path)
    ds = mpn.make_dataset(manifest, split="train", fields=["T"], target_fields=["v"])
    assert len(ds) == 3
    data = ds[0]
    # VTU carries no regions, so x is the T column alone on the re-read mesh
    assert data.x.shape == (4, 1) and data.x.dtype == torch.float32
    assert data.y.shape == (4, 2)
    assert data.pos.shape == (4, 3)
    assert data.edge_index.dtype == torch.int64
    assert data.edge_attr.shape[1] == 4
    assert ds.schema["x_columns"] == ["T"]
    from torch_geometric.loader import DataLoader

    batch = next(iter(DataLoader(ds, batch_size=3)))
    assert batch.num_graphs == 3 and batch.x.shape == (12, 1)


def test_reader_load_sample_returns_tensor_dict(tmp_path):
    pytest.importorskip("physicsnemo")
    import torch

    manifest = _manifest(tmp_path)
    reader = mpn.make_reader(manifest, split="train", fields=["T"])
    assert len(reader) == 3
    sample = reader._load_sample(0)
    assert isinstance(sample, dict)
    assert set(sample) == {"pos", "x", "edge_index", "edge_attr"}
    assert all(isinstance(t, torch.Tensor) for t in sample.values())
    assert sample["pos"].device.type == "cpu"
