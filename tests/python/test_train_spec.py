"""The training spec, the run-directory files and the model card -- the pure
half of the trainer (``meshioplusplus.physicsnemo._train``). Runs in the
default matrix with no framework installed.
"""

from __future__ import annotations

import json
import os
import re

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


# --------------------------------------------------------------------------- #
# the input normalizer must be as wide as x                                   #
# --------------------------------------------------------------------------- #
def test_the_input_normalizer_is_as_wide_as_the_recorded_x_columns():
    """Region one-hots widen `x`, and the normalizer has to widen with it.

    Before this was fixed the normalizer was built field-by-field, so with
    `Graph.Regions` on it was narrower than the batch. That has two faces and
    the first is the dangerous one:

    * one field, N columns -> numpy **broadcasts** the single value across
      every column, so the region one-hots are silently normalized by the
      field's own mean and standard deviation. No error, wrong model.
    * several fields, a different N -> a shape crash somewhere in the loop.

    A one-hot contributes mean 0 / std 1: it is already 0/1, and shifting it
    would stop it meaning "member".
    """
    spec = t.spec_from_dict(
        {
            "Manifest": "m.json",
            "Fields": ["q"],
            "TargetFields": ["T"],
            "Graph": {"Regions": True},
        }
    )
    schema = {
        "x_columns": ["q", "region:inlet", "region:wall"],
        "y_columns": ["T"],
        "x_sources": [
            {"name": "q", "source": "q", "kind": "data", "component": None},
            {
                "name": "region:inlet",
                "source": "inlet",
                "kind": "region",
                "component": None,
            },
            {
                "name": "region:wall",
                "source": "wall",
                "kind": "region",
                "component": None,
            },
        ],
    }
    stats = {"q_mean": [2.0], "q_std": [0.5], "T_mean": [1.0], "T_std": [1.0]}
    card = t.card_from_run(
        spec,
        schema,
        stats,
        {"edge_mean": [0.0], "edge_std": [1.0]},
        epoch=0,
        valid_loss=None,
        checkpoint="best.mdlus",
    )
    norms = t.normalizers_from_card(card)
    assert len(norms["x_mean"]) == len(schema["x_columns"])
    assert norms["x_mean"].tolist() == [2.0, 0.0, 0.0]
    assert norms["x_std"].tolist() == [0.5, 1.0, 1.0]


def test_a_multi_component_field_takes_its_own_component_statistics():
    spec = t.spec_from_dict(
        {"Manifest": "m.json", "Fields": ["v"], "TargetFields": ["T"]}
    )
    schema = {
        "x_columns": ["v_0", "v_1"],
        "y_columns": ["T"],
        "x_sources": [
            {"name": "v_0", "source": "v", "kind": "data", "component": 0},
            {"name": "v_1", "source": "v", "kind": "data", "component": 1},
        ],
    }
    stats = {"v_mean": [1.0, 9.0], "v_std": [2.0, 8.0], "T_mean": [0.0], "T_std": [1.0]}
    card = t.card_from_run(
        spec, schema, stats, {}, epoch=0, valid_loss=None, checkpoint="c.mdlus"
    )
    assert card["input_normalization"]["mean"] == [1.0, 9.0]
    assert card["input_normalization"]["std"] == [2.0, 8.0]


def test_a_card_without_an_edge_block_loads():
    """A grid model has no edges at all, so its card carries no edge block."""
    card = {
        "input_normalization": {"mean": [0.0], "std": [1.0]},
        "output_normalization": {"mean": [0.0], "std": [1.0]},
    }
    norms = t.normalizers_from_card(card)
    assert norms["e_mean"].tolist() == []
    assert norms["e_std"].tolist() == []


# --------------------------------------------------------------------------- #
# the srresnet family                                                         #
# --------------------------------------------------------------------------- #
SR_DOC = {
    "Manifest": "m.json",
    "Fields": ["T"],
    "TargetFields": ["T"],
    "Model": {"Name": "srresnet", "ScalingFactor": 2},
    "Grid": {"Resolution": [8, 8, 8]},
}


def test_srresnet_round_trips_and_emits_only_its_own_blocks():
    spec = t.spec_from_dict(SR_DOC)
    doc = t.spec_to_dict(spec)
    assert doc["Model"] == {
        "Name": "srresnet",
        "ScalingFactor": 2,
        "ConvLayerSize": 32,
        "ResidBlocks": 8,
        "LargeKernelSize": 7,
        "SmallKernelSize": 3,
        "ActivationFn": "prelu",
    }
    assert doc["Grid"]["Resolution"] == [8, 8, 8]
    # a superresolution document that listed graph aggregation would invite
    # someone to change it and wonder why nothing happened
    assert "Graph" not in doc
    assert t.spec_from_dict(doc) == spec


def test_meshgraphnet_still_emits_exactly_what_it_did():
    """The other family's document must be untouched by the second one."""
    doc = t.spec_to_dict(t.spec_from_dict(DOC))
    assert doc["Model"] == {
        "Name": "meshgraphnet",
        "ProcessorSize": 8,
        "HiddenDim": 32,
        "Aggregation": "sum",
    }
    assert "Grid" not in doc


@pytest.mark.parametrize(
    "doc, needle",
    [
        # a hyperparameter meant for the other family is refused, not ignored
        (
            {**SR_DOC, "Model": {"Name": "srresnet", "HiddenDim": 32}},
            "unknown key 'HiddenDim' in Model (with Name 'srresnet')",
        ),
        (
            {**DOC, "Model": {"Name": "meshgraphnet", "ScalingFactor": 2}},
            "unknown key 'ScalingFactor'",
        ),
        ({**DOC, "Grid": {"Resolution": [8, 8, 8]}}, "reads the Graph block, not Grid"),
        ({**SR_DOC, "Graph": {"Regions": True}}, "reads the Grid block, not Graph"),
        (
            {**SR_DOC, "Model": {"Name": "srresnet", "ScalingFactor": 3}},
            "Model.ScalingFactor must be one of 2, 4, 8",
        ),
        (
            {**SR_DOC, "Grid": {"Resolution": [8, 8, 8], "Squeeze": 2}},
            "Grid.Squeeze does not apply to 'srresnet'",
        ),
        (
            {**SR_DOC, "Grid": {"Resolution": [8, 8, 8], "CellSize": 0.1}},
            "exactly one of Grid.Resolution and Grid.CellSize",
        ),
        ({**SR_DOC, "Grid": {}}, "exactly one of Grid.Resolution and Grid.CellSize"),
        ({**SR_DOC, "Grid": {"Bogus": 1}}, "unknown key 'Bogus' in Grid"),
    ],
)
def test_cross_family_strictness_names_the_offender(doc, needle):
    with pytest.raises(ValueError, match=re.escape(needle)):
        t.spec_from_dict(doc)


def test_grid_kwargs_describe_the_pair():
    spec = t.spec_from_dict(SR_DOC)
    kwargs = spec.grid_kwargs()
    assert kwargs["scaling_factor"] == 2
    assert kwargs["coarse"]["resolution"] == [8, 8, 8]
    assert kwargs["fields"] == ["T"] and kwargs["target_fields"] == ["T"]


# --------------------------------------------------------------------------- #
# Graph.Proximity (v10.31.0)                                                  #
# --------------------------------------------------------------------------- #
def test_graph_proximity_round_trips_and_reaches_the_sample_kwargs():
    doc = dict(DOC)
    doc["Graph"] = {"Proximity": {"Method": "radius", "Radius": 0.05}}
    spec = t.spec_from_dict(doc)
    assert spec.proximity == {"method": "radius", "radius": 0.05}
    assert spec.graph_kwargs()["proximity"] == {"method": "radius", "radius": 0.05}
    assert t.spec_to_dict(spec)["Graph"]["Proximity"] == {
        "Method": "radius",
        "Radius": 0.05,
    }


def test_graph_proximity_carries_a_periodic_box():
    doc = dict(DOC)
    doc["Graph"] = {
        "Proximity": {"Method": "knn", "MaxNeighbors": 8, "BoxSize": [1.0, 2.0, 3.0]}
    }
    spec = t.spec_from_dict(doc)
    assert spec.proximity == {
        "method": "knn",
        "max_neighbors": 8,
        "box_size": [1.0, 2.0, 3.0],
    }
    assert t.spec_to_dict(spec)["Graph"]["Proximity"]["BoxSize"] == [1.0, 2.0, 3.0]


def test_an_unset_proximity_leaves_the_document_exactly_as_it_was():
    # The v10.30.0 shape, byte for byte: a key nobody asked for must not
    # appear in a spec that does not use the feature.
    spec = t.spec_from_dict(dict(DOC))
    assert spec.proximity is None
    assert "Proximity" not in t.spec_to_dict(spec)["Graph"]
    assert spec.graph_kwargs()["proximity"] is None


@pytest.mark.parametrize(
    "block, message",
    [
        ({"Method": "tree", "Radius": 1.0}, "Method must be one of"),
        ({"Method": "radius"}, "Radius must be a positive number"),
        ({"Method": "radius", "Radius": 0.0}, "Radius must be a positive number"),
        (
            {"Method": "radius", "Radius": 1.0, "MaxNeighbors": 4},
            "MaxNeighbors belongs to Method 'knn'",
        ),
        ({"Method": "knn", "Radius": 1.0}, "Radius belongs to Method 'radius'"),
        ({"Method": "knn", "MaxNeighbors": 0}, "MaxNeighbors must be an integer"),
        ({"Method": "radius", "Radius": 1.0, "Cutoff": 2}, "unknown key 'Cutoff'"),
        (
            {"Method": "radius", "Radius": 1.0, "BoxSize": "big"},
            "BoxSize must be a positive number",
        ),
        (
            {"Method": "radius", "Radius": 1.0, "BoxSize": [1.0, -2.0]},
            "BoxSize must be a positive number",
        ),
    ],
)
def test_a_malformed_proximity_block_names_the_offender(block, message):
    doc = dict(DOC)
    doc["Graph"] = {"Proximity": block}
    with pytest.raises(ValueError, match=re.escape(message)):
        t.spec_from_dict(doc)


def test_proximity_is_a_graph_key_and_not_a_grid_one():
    doc = dict(DOC)
    doc["Model"] = {"Name": "srresnet"}
    doc["Grid"] = {"Resolution": [4, 4, 4], "Proximity": {"Method": "radius"}}
    with pytest.raises(ValueError, match="unknown key 'Proximity' in Grid"):
        t.spec_from_dict(doc)


# --------------------------------------------------------------------------- #
# the family table (v10.40.0): every family x every block                     #
# --------------------------------------------------------------------------- #
#: A block-legal probe per block. The per-block `_check_keys` runs BEFORE the
#: family's own refusal, so a probe must pass it -- an unknown key would be
#: refused for the wrong reason and the test would prove nothing.
_BLOCK_PROBES = {
    "Graph": {"Regions": True},
    "Grid": {"Resolution": [8, 8, 8]},
    "Operator": {"Parameters": ["a"]},
}
#: Families whose own block needs more than the probe to be valid at all.
_OWN_BLOCKS = {"afno": {"Resolution": [8, 8, 8], "Squeeze": 2}}


def _minimal_doc(fam):
    doc = {"Manifest": "m.json", "TargetFields": ["T"], "Model": {"Name": fam.name}}
    # an operator family's inputs are its parameters, not a data array
    if fam.block != "Operator":
        doc["Fields"] = ["T"]
    doc[fam.block] = _OWN_BLOCKS.get(fam.name, _BLOCK_PROBES[fam.block])
    return doc


def test_the_block_table_covers_every_family():
    assert {f.block for f in t._FAMILIES} <= set(t._BLOCKS)
    assert set(_BLOCK_PROBES) == set(t._BLOCKS)
    assert t._MODELS == tuple(f.name for f in t._FAMILIES)


@pytest.mark.parametrize("fam", t._FAMILIES, ids=lambda f: f.name)
def test_every_family_parses_with_its_own_block(fam):
    spec = t.spec_from_dict(_minimal_doc(fam))
    assert spec.model_name == fam.name
    assert t.spec_from_dict(t.spec_to_dict(spec)) == spec


@pytest.mark.parametrize(
    "fam, block",
    [(f, b) for f in t._FAMILIES for b in t._BLOCKS if b != f.block],
    ids=lambda v: v if isinstance(v, str) else v.name,
)
def test_every_foreign_block_is_refused_by_name(fam, block):
    """The refusal loops over EVERY block, not "the other one": with a
    third block a binary `unwanted = ...` lets a foreign block pass silently."""
    doc = _minimal_doc(fam)
    doc[block] = _BLOCK_PROBES[block]
    needle = f"a '{fam.name}' model reads the {fam.block} block, not {block}"
    with pytest.raises(ValueError, match=re.escape(needle)):
        t.spec_from_dict(doc)


def test_grid_float32_is_read_from_the_grid_block():
    """`Grid.Float32` was accepted by the key check and emitted by
    `spec_to_dict`, but READ from the Graph block -- so `false` on an srresnet
    was silently ignored and did not round-trip."""
    doc = {**SR_DOC, "Grid": {"Resolution": [8, 8, 8], "Float32": False}}
    spec = t.spec_from_dict(doc)
    assert spec.float32 is False
    assert spec.grid_kwargs()["float32"] is False
    emitted = t.spec_to_dict(spec)
    assert emitted["Grid"]["Float32"] is False
    assert t.spec_from_dict(emitted) == spec
    with pytest.raises(ValueError, match="Grid.Float32 must be"):
        t.spec_from_dict({**SR_DOC, "Grid": {"Resolution": [8, 8, 8], "Float32": 1}})


def test_require_frameworks_asks_only_for_the_family_s_own(monkeypatch):
    from meshioplusplus import _gpu

    asked = []
    monkeypatch.setattr(_gpu, "_importable", lambda module: asked.append(module))
    for fam in t._FAMILIES:
        asked.clear()
        with pytest.raises(ImportError) as excinfo:
            t.require_frameworks("op", fam.name)
        assert asked[0] == fam.frameworks[0]
        assert t._FRAMEWORK_HINTS[fam.frameworks[0]] in str(excinfo.value)
        assert "meshioplusplus[" not in str(excinfo.value)
    with pytest.raises(ValueError, match="Model.Name must be one of"):
        t.require_frameworks("op", "gpt")


# --------------------------------------------------------------------------- #
# the fno / afno / deeponet families (v10.40.0)                               #
# --------------------------------------------------------------------------- #
FNO_DOC = {
    "Manifest": "m.json",
    "Fields": ["K"],
    "TargetFields": ["p"],
    "Model": {"Name": "fno", "NumFnoModes": 12},
    "Grid": {"Resolution": [32, 32, 1], "Squeeze": 2, "SqueezeIndex": 0},
}
AFNO_DOC = {
    "Manifest": "m.json",
    "Fields": ["c0", "u", "v"],
    "TargetFields": ["c"],
    "Model": {"Name": "afno", "PatchSize": [8, 8]},
    "Grid": {"Resolution": [63, 63, 1], "Squeeze": 2, "SqueezeIndex": 0},
}
DEEP_DOC = {
    "Manifest": "m.json",
    "TargetFields": ["w"],
    "Model": {"Name": "deeponet"},
    "Operator": {"Parameters": ["Load", "Modulus", "PoissonRatio"]},
}


def test_fno_round_trips_and_emits_only_its_own_blocks():
    spec = t.spec_from_dict(FNO_DOC)
    doc = t.spec_to_dict(spec)
    assert doc["Model"] == {
        "Name": "fno",
        "LatentChannels": 32,
        "NumFnoLayers": 4,
        "NumFnoModes": 12,
        "SpectralPadding": 8,
        "PaddingType": "constant",
        "ActivationFn": "gelu",
        "DecoderLayers": 1,
        "DecoderLayerSize": 32,
        "DecoderActivationFn": "silu",
        "CoordFeatures": True,
    }
    assert doc["Grid"]["Squeeze"] == 2 and doc["Grid"]["SqueezeIndex"] == 0
    assert "Graph" not in doc and "Operator" not in doc
    assert t.spec_from_dict(doc) == spec
    # a 3-D fno needs no squeeze at all
    three = t.spec_from_dict({**FNO_DOC, "Grid": {"Resolution": [8, 8, 8]}})
    assert three.squeeze is None and three.grid_kwargs()["squeeze"] is None


def test_afno_round_trips_and_emits_only_its_own_blocks():
    spec = t.spec_from_dict(AFNO_DOC)
    doc = t.spec_to_dict(spec)
    assert doc["Model"] == {
        "Name": "afno",
        "PatchSize": [8, 8],
        "EmbedDim": 256,
        "Depth": 4,
        "MlpRatio": 4.0,
        "DropRate": 0.0,
        "NumBlocks": 16,
        "SparsityThreshold": 0.01,
        "HardThresholdingFraction": 1.0,
    }
    assert spec.patch_size == (8, 8)
    assert "Graph" not in doc and "Operator" not in doc
    assert t.spec_from_dict(doc) == spec


def test_deeponet_round_trips_and_emits_only_its_own_blocks():
    spec = t.spec_from_dict(DEEP_DOC)
    doc = t.spec_to_dict(spec)
    assert doc["Model"] == {
        "Name": "deeponet",
        "BranchLayers": 4,
        "BranchLayerSize": 128,
        "TrunkLayers": 4,
        "TrunkLayerSize": 128,
        "Width": 64,
        "DecoderType": "mlp",
        "DecoderWidth": 128,
        "DecoderLayers": 2,
        "DecoderActivationFn": "relu",
        "ActivationFn": "silu",
    }
    assert doc["Operator"] == {
        "Parameters": ["Load", "Modulus", "PoissonRatio"],
        "Trunk": "points",
        "TrunkMethod": "farthest",
        "TrunkSeed": 0,
        "Float32": True,
    }
    assert doc["Fields"] == [] and "Grid" not in doc and "Graph" not in doc
    assert t.spec_from_dict(doc) == spec
    budget = t.spec_from_dict(
        {
            **DEEP_DOC,
            "Operator": {"Parameters": ["Load"], "Trunk": "budget", "TrunkCount": 64},
        }
    )
    assert t.spec_to_dict(budget)["Operator"]["TrunkCount"] == 64
    assert budget.operator_kwargs() == {
        "parameter_names": ["Load"],
        "target_fields": ["w"],
        "trunk": "budget",
        "trunk_count": 64,
        "trunk_method": "farthest",
        "trunk_seed": 0,
        "float32": True,
    }


def test_shared_keys_take_the_family_s_own_default():
    """`ActivationFn` and `DecoderLayers` are shared PascalCase keys whose
    constructor defaults differ per family; `default_spec` re-validates
    through the document, so a shared dataclass field would hand fno
    srresnet's prelu. The family-prefixed fields are what prevent that."""
    assert t.spec_from_dict(FNO_DOC).fno_activation_fn == "gelu"
    assert t.spec_from_dict(FNO_DOC).fno_decoder_layers == 1
    assert t.spec_from_dict(DEEP_DOC).deeponet_activation_fn == "silu"
    assert t.spec_from_dict(DEEP_DOC).decoder_layers == 2
    assert t.spec_from_dict(SR_DOC).activation_fn == "prelu"
    fno = t.default_spec("m.json", ["K"], ["p"], model_name="fno", resolution=(8, 8, 8))
    assert fno.fno_activation_fn == "gelu" and fno.activation_fn == "prelu"
    deep = t.default_spec(
        "m.json", [], ["w"], model_name="deeponet", parameters=("Load",)
    )
    assert deep.deeponet_activation_fn == "silu" and deep.decoder_layers == 2


def test_grid_kwargs_of_a_resolution_preserving_family_pair_the_grid_with_itself():
    """fno/afno have no ScalingFactor: `upscale_samples(1)` is the identity,
    so the coarse/fine pairing is reused unchanged with a factor of one."""
    assert t.spec_from_dict(FNO_DOC).grid_kwargs()["scaling_factor"] == 1
    assert t.spec_from_dict(AFNO_DOC).grid_kwargs()["scaling_factor"] == 1
    assert t.spec_from_dict(SR_DOC).grid_kwargs()["scaling_factor"] == 2


@pytest.mark.parametrize(
    "doc, needle",
    [
        (
            {**FNO_DOC, "Model": {"Name": "fno", "HiddenDim": 32}},
            "unknown key 'HiddenDim' in Model (with Name 'fno')",
        ),
        (
            {**FNO_DOC, "Model": {"Name": "fno", "ScalingFactor": 2}},
            "unknown key 'ScalingFactor' in Model (with Name 'fno')",
        ),
        (
            {**AFNO_DOC, "Grid": {"Resolution": [63, 63, 1]}},
            "a 'afno' model is 2-D only (AFNO patches a fixed (H, W) image); give "
            "Grid.Squeeze the world axis to collapse (0, 1 or 2)",
        ),
        (
            {**SR_DOC, "Grid": {"Resolution": [8, 8, 8], "Squeeze": 2}},
            "Grid.Squeeze does not apply to 'srresnet'",
        ),
        (
            {**FNO_DOC, "Operator": {"Parameters": ["a"]}},
            "a 'fno' model reads the Grid block, not Operator",
        ),
        (
            {**DEEP_DOC, "Grid": {"Resolution": [8, 8, 8]}},
            "a 'deeponet' model reads the Operator block, not Grid",
        ),
        (
            {**DEEP_DOC, "Augmentation": {"Seed": 1, "Rotation": {"Axis": "z"}}},
            "Augmentation applies to the graph families; a 'deeponet' shares one "
            "trunk across the batch",
        ),
        (
            {**FNO_DOC, "Augmentation": {"Seed": 1, "Rotation": {"Axis": "z"}}},
            "a 'fno' samples a fixed lattice",
        ),
        (
            {**DEEP_DOC, "Operator": {"Parameters": ["Load"], "Trunk": "budget"}},
            "Operator.TrunkCount is required with Trunk 'budget'",
        ),
        (
            {**DEEP_DOC, "Operator": {"Parameters": ["Load"], "TrunkCount": 8}},
            "Operator.TrunkCount belongs to Trunk 'budget'",
        ),
        (
            {
                **DEEP_DOC,
                "Model": {"Name": "deeponet", "DecoderType": "temporal_projection"},
            },
            "Model.DecoderType must be one of mlp; 'temporal_projection' needs an "
            "output_window",
        ),
        (
            {**DEEP_DOC, "Model": {"Name": "deeponet", "DecoderType": "conv"}},
            "'conv' needs a spatial branch",
        ),
        (
            {**DEEP_DOC, "Fields": ["T"]},
            "a 'deeponet' model takes its inputs from Operator.Parameters, not Fields",
        ),
        ({**DEEP_DOC, "Operator": {}}, "Operator.Parameters must name at least one"),
        (
            {**DEEP_DOC, "Operator": {"Parameters": ["a", "a"]}},
            "Operator.Parameters must not repeat",
        ),
        (
            {**DEEP_DOC, "Operator": {"Parameters": ["a"], "TrunkMethod": "kd"}},
            "Operator.TrunkMethod must be one of farthest, grid, random",
        ),
        (
            {**AFNO_DOC, "Model": {"Name": "afno", "PatchSize": [8]}},
            "Model.PatchSize must be two positive integers",
        ),
        ({**FNO_DOC, "Grid": {"Squeeze": 2}}, "exactly one of Grid.Resolution"),
    ],
)
def test_new_family_strictness_names_the_offender(doc, needle):
    with pytest.raises(ValueError, match=re.escape(needle)):
        t.spec_from_dict(doc)


def test_grid_card_records_the_squeeze_contract():
    """The card of a 2-D grid family carries the layout of the REMAINING axes,
    the spatial rank, the sample shape and how to expand a plane back --
    without them a 2-D checkpoint could neither size its model nor write its
    answer onto the thin lattice."""
    from meshioplusplus.physicsnemo import train as trainer

    schema = {
        "x_channels": ["K"],
        "y_channels": ["p"],
        "squeeze": 2,
        "squeeze_index": 0,
        "x_shape": [16, 16],
        "y_shape": [16, 16],
        "expand_size": 2,
        "coarse": None,
        "fine": None,
    }
    stats = {"x_mean": [0.0], "x_std": [1.0], "y_mean": [0.0], "y_std": [1.0]}
    card = trainer.grid_card_from_run(
        t.spec_from_dict(FNO_DOC),
        schema,
        stats,
        epoch=0,
        valid_loss=None,
        checkpoint="c",
    )
    assert card["layout"] == "channels_first_yx" and card["spatial_ndim"] == 2
    assert card["squeeze"] == 2 and card["squeeze_index"] == 0
    assert card["x_shape"] == [16, 16]
    assert card["expand_axis"] == 2 and card["expand_size"] == 2
    assert card["model"]["name"] == "fno" and card["model"]["dimension"] == 2
    assert card["model"]["num_fno_modes"] == 12
    afno = trainer.grid_card_from_run(
        t.spec_from_dict(AFNO_DOC),
        schema,
        stats,
        epoch=0,
        valid_loss=None,
        checkpoint="c",
    )
    assert afno["model"]["inp_shape"] == [16, 16] and afno["model"]["patch_size"] == [
        8,
        8,
    ]
    # the srresnet card is untouched in meaning: 3-D, the zyx layout
    three = dict(
        schema, squeeze=None, squeeze_index=None, x_shape=[5, 5, 5], expand_size=None
    )
    sr = trainer.grid_card_from_run(
        t.spec_from_dict(SR_DOC), three, stats, epoch=0, valid_loss=None, checkpoint="c"
    )
    assert sr["layout"] == "channels_first_zyx" and sr["spatial_ndim"] == 3
    assert sr["model"]["scaling_factor"] == 2 and "dimension" not in sr["model"]
    assert t.normalizers_from_card(card)["e_mean"].tolist() == []


def test_afno_divisibility_is_checked_before_torch_is_imported():
    """A patch that does not divide the sample shape is refused naming the
    off-by-one (`Grid.Resolution [n, ...]` gives n + 1 samples) BEFORE torch
    is imported -- which is why this runs in the default matrix: were the
    check after the constructor, this box, having no physicsnemo, would see
    an ImportError instead."""
    from meshioplusplus.physicsnemo import train as trainer

    spec = t.spec_from_dict(AFNO_DOC)
    schema = {"x_channels": ["c0"], "y_channels": ["c"], "x_shape": [65, 65]}
    with pytest.raises(ValueError, match=r"n \+ 1 sample points"):
        trainer._build_afno(spec, schema)
    with pytest.raises(ValueError, match="AFNO is 2-D only"):
        trainer._build_afno(spec, {**schema, "x_shape": [8, 8, 8]})
    odd = t.spec_from_dict(
        {
            **AFNO_DOC,
            "Model": {
                "Name": "afno",
                "PatchSize": [8, 8],
                "EmbedDim": 12,
                "NumBlocks": 8,
            },
        }
    )
    with pytest.raises(
        ValueError, match="EmbedDim 12 must be divisible by Model.NumBlocks 8"
    ):
        trainer._build_afno(odd, {**schema, "x_shape": [64, 64]})


def test_operator_card_records_the_parameter_and_trunk_contract():
    from meshioplusplus.physicsnemo import train as trainer

    spec = t.spec_from_dict(
        {
            **DEEP_DOC,
            "Operator": {
                "Parameters": ["Load", "E"],
                "Trunk": "budget",
                "TrunkCount": 3,
            },
        }
    )
    schema = {
        "parameter_columns": ["Load", "E_0", "E_1"],
        "y_columns": ["w"],
        "num_points": 10,
        "num_trunk": 3,
    }
    stats = {
        "params_mean": [1.0, 2.0, 3.0],
        "params_std": [0.0, 1.0, 1.0],
        "trunk_mean": [0.5, 0.5, 0.0],
        "trunk_std": [0.3, 0.3, 0.0],
        "y_mean": [0.1],
        "y_std": [0.2],
    }
    card = trainer.operator_card_from_run(
        spec,
        schema,
        stats,
        [1, 4, 7],
        epoch=2,
        valid_loss=0.5,
        checkpoint="/r/best.mdlus",
    )
    assert card["model"]["name"] == "deeponet"
    assert card["model"]["in_parameters"] == 3 and card["model"]["out_channels"] == 1
    assert card["parameters"] == ["Load", "E_0", "E_1"]
    assert card["parameter_names"] == ["Load", "E"]
    assert card["trunk_indices"] == [1, 4, 7] and card["num_points"] == 10
    assert card["input_normalization"]["std"][0] == t.STATS_STD_FLOOR
    assert card["trunk_normalization"]["std"][2] == t.STATS_STD_FLOOR
    assert card["fields"] == [] and card["operator"]["trunk_count"] == 3
    norms = t.normalizers_from_card(card)
    assert (
        norms["x_mean"].tolist() == [1.0, 2.0, 3.0] and norms["e_mean"].tolist() == []
    )


def test_the_viewer_family_table_matches_the_python_one():
    """`src/viewer/src/dataset/families.ts` drives the launch form; its
    family names are pinned equal, in order, to `_MODELS` (the `manifest.ts`
    parity precedent), so a family added on one side only is a red build."""
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(
        here, "..", "..", "src", "viewer", "src", "dataset", "families.ts"
    )
    with open(path, encoding="utf-8") as fh:
        names = re.findall(r"name: '([a-z]+)'", fh.read())
    assert names == list(t._MODELS)
