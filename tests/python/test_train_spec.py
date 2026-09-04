"""The training spec, the run-directory files and the model card -- the pure
half of the trainer (``meshioplusplus.physicsnemo._train``). Runs in the
default matrix with no framework installed.
"""

from __future__ import annotations

import json
import os

import numpy as np
import pytest

from meshioplusplus.physicsnemo import _train as t

DOC = {
    "Version": 1,
    "Manifest": "dataset_manifest.json",
    "RunDir": "runs/a",
    "Fields": ["q_scaled"],
    "TargetFields": ["T"],
    "Epochs": 3,
    "Model": {"HiddenDim": 32},
    "Graph": {"Regions": True, "TargetOffset": 1, "TargetDelta": True},
    "Tags": ["smoke"],
}


def test_spec_round_trips_and_fills_defaults():
    spec = t.spec_from_dict(DOC)
    assert spec.fields == ("q_scaled",) and spec.target_fields == ("T",)
    assert spec.epochs == 3 and spec.batch_size == 8 and spec.learning_rate == 1e-3
    assert spec.hidden_dim == 32 and spec.processor_size == 8
    assert spec.regions is True and spec.target_offset == 1 and spec.target_delta
    doc = t.spec_to_dict(spec)
    assert doc["Model"] == {
        "Name": "meshgraphnet",
        "ProcessorSize": 8,
        "HiddenDim": 32,
        "Aggregation": "sum",
    }
    assert "Read" not in doc and "Notes" not in doc and doc["Tags"] == ["smoke"]
    assert t.spec_from_dict(doc) == spec
    kwargs = spec.graph_kwargs()
    assert kwargs["fields"] == ["q_scaled"] and kwargs["target_delta"] is True


@pytest.mark.parametrize(
    "doc, message",
    [
        ({**DOC, "Bogus": 1}, "unknown key 'Bogus'"),
        ({**DOC, "Version": 2}, "unsupported Version"),
        ({k: v for k, v in DOC.items() if k != "Manifest"}, "Manifest is required"),
        ({**DOC, "Fields": []}, "Fields must name"),
        ({**DOC, "TargetFields": "T", "Fields": "q"}, None),  # a string is one name
        ({**DOC, "Model": {"Nope": 1}}, "unknown key 'Nope' in Model"),
        ({**DOC, "Model": {"Name": "gpt"}}, "Model.Name"),
        ({**DOC, "Graph": {"Kind": "edge"}}, "Graph.Kind"),
        ({**DOC, "Graph": {"Regions": "yes"}}, "Graph.Regions must be"),
        ({**DOC, "Epochs": 0}, "Epochs must be"),
        ({**DOC, "LearningRate": -1}, "LearningRate"),
        ({**DOC, "Device": "tpu"}, "Device must be"),
        ({**DOC, "Read": 3}, "Read must be"),
        ("[", "not valid JSON"),
    ],
)
def test_strict_validation_names_the_offender(doc, message):
    if message is None:
        spec = t.load_spec(doc)
        assert spec.fields == ("q",) and spec.target_fields == ("T",)
        return
    with pytest.raises(ValueError, match=message):
        t.load_spec(doc)


def test_load_is_tri_modal_and_a_file_anchors_relative_paths(tmp_path):
    spec_path = tmp_path / "campaign" / "spec.json"
    spec_path.parent.mkdir()
    spec_path.write_text(json.dumps(DOC), encoding="utf-8")
    from_file = t.load_spec(str(spec_path))
    assert from_file.base_dir == str(spec_path.parent)
    assert from_file.resolved_manifest() == str(
        spec_path.parent / "dataset_manifest.json"
    )
    assert from_file.resolved_run_dir() == str(spec_path.parent / "runs" / "a")
    from_text = t.load_spec(json.dumps(DOC))
    assert from_text.base_dir is None
    assert from_text.resolved_manifest() == os.path.join(
        os.getcwd(), "dataset_manifest.json"
    )
    assert from_text == from_file  # base_dir is excluded from equality
    assert t.load_spec(from_file) is from_file
    saved = tmp_path / "out.json"
    t.save_spec(from_file, str(saved))
    assert t.load_spec(str(saved)) == from_file


def test_default_spec_applies_overrides_through_validation():
    spec = t.default_spec(
        "m.json", ["q"], "T", run_dir="runs/x", epochs=5, hidden_dim=16
    )
    assert spec.epochs == 5 and spec.hidden_dim == 16 and spec.run_dir == "runs/x"
    with pytest.raises(ValueError, match="unknown spec field"):
        t.default_spec("m.json", ["q"], ["T"], bogus=1)
    with pytest.raises(ValueError, match="Epochs must be"):
        t.default_spec("m.json", ["q"], ["T"], epochs=0)


def test_metrics_and_progress_files(tmp_path):
    run = str(tmp_path)
    assert t.read_metrics(run) == []
    for epoch in range(3):
        t.append_metrics(
            run,
            t.metrics_row(
                epoch,
                1.0 / (epoch + 1),
                None if epoch == 0 else 0.5,
                1e-3,
                epoch * 2.0,
                2.0,
            ),
        )
    # a torn last line (a trainer mid-write) is skipped, never fatal
    with open(os.path.join(run, t.METRICS_FILE), "a", encoding="utf-8") as fh:
        fh.write('{"epoch": 3, "train_lo')
    rows = t.read_metrics(run)
    assert [r["epoch"] for r in rows] == [0, 1, 2]
    assert rows[0]["valid_loss"] is None and rows[1]["valid_loss"] == 0.5
    assert [r["epoch"] for r in t.read_metrics(run, since_epoch=2)] == [2]
    t.write_json_atomic(os.path.join(run, "p.json"), {"a": 1})
    assert t.read_json(os.path.join(run, "p.json")) == {"a": 1}
    assert t.read_json(os.path.join(run, "missing.json"), "dflt") == "dflt"
    assert not [n for n in os.listdir(run) if n.startswith(".tmp-")]


def test_stat_vectors_concatenate_in_order_and_floor_the_std():
    stats = {"q_mean": [1.0], "q_std": [0.0], "v_mean": [2.0, 3.0], "v_std": [0.5, 4.0]}
    assert t.stat_vectors(stats, ["q", "v"], "mean").tolist() == [1.0, 2.0, 3.0]
    std = t.stat_vectors(stats, ["q", "v"], "std")
    assert std[0] == t.STATS_STD_FLOOR and std[1:].tolist() == [0.5, 4.0]
    with pytest.raises(KeyError, match="no 'w_mean'"):
        t.stat_vectors(stats, ["w"], "mean")


def test_card_round_trip_and_checkpoint_listing(tmp_path):
    spec = t.spec_from_dict({**DOC, "Graph": {}})
    schema = {"x_columns": ["q_scaled"], "y_columns": ["T"], "edge_features": True}
    node_stats = {
        "q_scaled_mean": [0.1],
        "q_scaled_std": [0.0],
        "T_mean": [0.5],
        "T_std": [0.25],
    }
    edge_stats = {"edge_mean": [0.0, 0.0, 0.0, 0.1], "edge_std": [0.1, 0.1, 0.0, 0.05]}
    card = t.card_from_run(
        spec,
        schema,
        node_stats,
        edge_stats,
        epoch=7,
        valid_loss=0.01,
        checkpoint="/x/checkpoints/best.mdlus",
    )
    assert (
        card["model"]["input_dim_nodes"] == 1
        and card["model"]["input_dim_edges"] == 4
        and card["model"]["output_dim"] == 1
    )
    assert card["checkpoint"] == "best.mdlus" and card["epoch"] == 7
    assert card["input_normalization"]["std"] == [t.STATS_STD_FLOOR]
    norms = t.normalizers_from_card(card)
    assert norms["y_mean"].tolist() == [0.5] and norms["e_std"][2] == t.STATS_STD_FLOOR
    assert all(isinstance(v, np.ndarray) for v in norms.values())
    # delta targets normalize with the diff stats
    delta = t.spec_from_dict(DOC)
    with pytest.raises(KeyError, match="T_diff_mean"):
        t.card_from_run(
            delta,
            schema,
            node_stats,
            edge_stats,
            epoch=0,
            valid_loss=None,
            checkpoint="c.mdlus",
        )

    ckpt = tmp_path / t.CHECKPOINT_DIR
    ckpt.mkdir()
    for name, epoch in (
        ("Model.0.9.mdlus", 9),
        ("Model.0.19.mdlus", 19),
        ("best.mdlus", 9),
        ("final.mdlus", 29),
    ):
        (ckpt / name).write_bytes(b"x")
        t.write_json_atomic(
            t.card_path(str(ckpt / name)),
            {"epoch": epoch, "valid_loss": 1.0 / (epoch + 1)},
        )
    (ckpt / "orphan.mdlus").write_bytes(b"x")  # no card: epoch unknown, sorted last
    listed = t.list_checkpoints(str(tmp_path), best=str(ckpt / "best.mdlus"))
    assert [c["name"] for c in listed] == [
        "Model.0.9.mdlus",
        "best.mdlus",
        "Model.0.19.mdlus",
        "final.mdlus",
        "orphan.mdlus",
    ]
    assert [c["kind"] for c in listed] == [
        "periodic",
        "best",
        "periodic",
        "final",
        "periodic",
    ]
    assert [c["is_best"] for c in listed] == [False, True, False, False, False]
    assert listed[-1]["epoch"] is None
