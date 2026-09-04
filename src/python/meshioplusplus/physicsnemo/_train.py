"""The training spec and run-directory layout: the pure half of the trainer.

Everything here is stdlib + numpy and Python 3.8 -- the spec parser, the run
layout, the metrics/progress files, the model card -- so it is tested in the
default CI matrix with no framework installed. ``train.py`` (the module run
as ``python -m meshioplusplus.physicsnemo.train``) is the only place torch,
PhysicsNeMo and PyTorch Geometric are imported.

**The spec is a hand-editable input, so it is a settings-family document**
(PascalCase keys, ``"Version": 1``, strict unknown-key refusal -- the
``DatasetManifest``/pipeline convention); everything a run *writes*
(``job.json``, ``progress.json``, ``metrics.jsonl``, ``*.card.json``) is a
machine artefact and is snake_case (the ``write_dataset`` precedent).

Run directory layout (``RunDir``)::

    spec.json           the spec the run was started from (PascalCase)
    job.json            written by the job manager (mcp/_jobs.py), not here
    progress.json       {epoch, epochs, best_epoch, best_valid_loss, ...}
    metrics.jsonl       one JSON object per epoch
    log.txt             stdout+stderr of the trainer (the manager redirects)
    node_stats.json     the PhysicsNeMo-convention normalization stats
    edge_stats.json
    checkpoints/        <Model>.0.<epoch>.mdlus + checkpoint.0.<epoch>.pt
                        (physicsnemo's save_checkpoint), best.mdlus,
                        final.mdlus -- every .mdlus with a .card.json sidecar
    predictions/        <entry_id>.vtu written by predict()

The **model card** (``<checkpoint>.card.json``, the Kratos
PhysicsNeMoApplication convention) records what a checkpoint's channels
mean -- field names, column order, the normalization the model was trained
under -- so a prediction can never silently write normalized numbers onto a
physical field.
"""

from __future__ import annotations

import glob
import json
import os
import tempfile
import time
from dataclasses import dataclass, field, replace
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..__about__ import __version__

SPEC_VERSION = 1
SPEC_FILE = "spec.json"
JOB_FILE = "job.json"
PROGRESS_FILE = "progress.json"
METRICS_FILE = "metrics.jsonl"
LOG_FILE = "log.txt"
CHECKPOINT_DIR = "checkpoints"
NODE_STATS_FILE = "node_stats.json"
EDGE_STATS_FILE = "edge_stats.json"
PREDICTIONS_DIR = "predictions"
BEST_CHECKPOINT = "best.mdlus"
FINAL_CHECKPOINT = "final.mdlus"
CARD_SUFFIX = ".card.json"
STATS_STD_FLOOR = 1e-8

_ERR = "meshio++: train: "
_TOP_KEYS = (
    "Version",
    "Manifest",
    "RunDir",
    "Fields",
    "TargetFields",
    "TrainSplit",
    "ValidSplit",
    "Epochs",
    "BatchSize",
    "LearningRate",
    "Seed",
    "CheckpointEvery",
    "Device",
    "Model",
    "Graph",
    "Read",
    "Notes",
    "Tags",
)
_MODEL_KEYS = ("Name", "ProcessorSize", "HiddenDim", "Aggregation")
_GRAPH_KEYS = (
    "Regions",
    "Kind",
    "Undirected",
    "EdgeFeatures",
    "Float32",
    "TargetOffset",
    "TargetDelta",
)
_MODELS = ("meshgraphnet",)
_KINDS = ("node", "cell")
_AGGREGATIONS = ("sum", "mean")
_DEVICES_PREFIX = ("auto", "cpu", "cuda")


def _check_keys(obj, where, allowed):
    for key in obj:
        if key not in allowed:
            raise ValueError(
                f"{_ERR}unknown key '{key}' in {where} (known: {', '.join(allowed)})"
            )


def _names(value, where):
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, (list, tuple)) or not all(
        isinstance(v, str) and v for v in value
    ):
        raise ValueError(f"{_ERR}{where} must be a list of non-empty strings")
    return tuple(value)


def _int(value, where, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{_ERR}{where} must be an integer >= {minimum}")
    return int(value)


def _bool(value, where):
    if not isinstance(value, bool):
        raise ValueError(f"{_ERR}{where} must be true or false")
    return value


@dataclass(frozen=True)
class TrainSpec:
    """A training run's inputs -- the PascalCase document, typed.

    ``base_dir`` (never serialized, excluded from equality) is the directory
    a spec loaded from a file resolves its relative ``Manifest``/``RunDir``
    against; a dict/JSON-text load resolves against the CWD.
    """

    manifest: str
    fields: Tuple[str, ...]
    target_fields: Tuple[str, ...]
    run_dir: str = "runs/run"
    train_split: str = "train"
    valid_split: str = "valid"
    epochs: int = 100
    batch_size: int = 8
    learning_rate: float = 1e-3
    seed: int = 0
    checkpoint_every: int = 10
    device: str = "auto"
    model_name: str = "meshgraphnet"
    processor_size: int = 8
    hidden_dim: int = 64
    aggregation: str = "sum"
    regions: bool = False
    kind: str = "node"
    undirected: bool = True
    edge_features: bool = True
    float32: bool = True
    target_offset: int = 0
    target_delta: bool = False
    read: Dict[str, Any] = field(default_factory=dict)
    notes: Optional[str] = None
    tags: Tuple[str, ...] = ()
    base_dir: Optional[str] = field(default=None, compare=False, repr=False)

    def resolved_manifest(self) -> str:
        return _resolve_against(self.manifest, self.base_dir)

    def resolved_run_dir(self) -> str:
        return _resolve_against(self.run_dir, self.base_dir)

    def graph_kwargs(self) -> Dict[str, Any]:
        """The :func:`graph_sample`/``make_dataset`` keyword arguments."""
        return {
            "fields": list(self.fields),
            "target_fields": list(self.target_fields),
            "regions": self.regions,
            "kind": self.kind,
            "undirected": self.undirected,
            "edge_features": self.edge_features,
            "float32": self.float32,
            "target_offset": self.target_offset,
            "target_delta": self.target_delta,
        }


def _resolve_against(path, base_dir):
    if os.path.isabs(path):
        return os.path.normpath(path)
    return os.path.normpath(os.path.join(base_dir or os.getcwd(), path))


def spec_from_dict(doc, *, base_dir=None) -> TrainSpec:
    """Parse a PascalCase spec document, strictly."""
    if not isinstance(doc, dict):
        raise ValueError(f"{_ERR}the spec must be a JSON object")
    _check_keys(doc, "the spec", _TOP_KEYS)
    version = doc.get("Version", SPEC_VERSION)
    if version != SPEC_VERSION:
        raise ValueError(
            f"{_ERR}unsupported Version {version!r} (this build knows {SPEC_VERSION})"
        )
    if not isinstance(doc.get("Manifest"), str) or not doc["Manifest"]:
        raise ValueError(f"{_ERR}Manifest is required (a dataset manifest path)")
    fields = _names(doc.get("Fields", ()), "Fields")
    targets = _names(doc.get("TargetFields", ()), "TargetFields")
    if not fields:
        raise ValueError(f"{_ERR}Fields must name at least one input array")
    if not targets:
        raise ValueError(f"{_ERR}TargetFields must name at least one target array")
    model = doc.get("Model", {})
    if not isinstance(model, dict):
        raise ValueError(f"{_ERR}Model must be an object")
    _check_keys(model, "Model", _MODEL_KEYS)
    graph = doc.get("Graph", {})
    if not isinstance(graph, dict):
        raise ValueError(f"{_ERR}Graph must be an object")
    _check_keys(graph, "Graph", _GRAPH_KEYS)
    read = doc.get("Read", {})
    if not isinstance(read, dict):
        raise ValueError(f"{_ERR}Read must be an object of read() keyword arguments")
    name = str(model.get("Name", "meshgraphnet")).lower()
    if name not in _MODELS:
        raise ValueError(f"{_ERR}Model.Name must be one of {', '.join(_MODELS)}")
    aggregation = str(model.get("Aggregation", "sum"))
    if aggregation not in _AGGREGATIONS:
        raise ValueError(
            f"{_ERR}Model.Aggregation must be one of {', '.join(_AGGREGATIONS)}"
        )
    kind = str(graph.get("Kind", "node"))
    if kind not in _KINDS:
        raise ValueError(f"{_ERR}Graph.Kind must be one of {', '.join(_KINDS)}")
    device = str(doc.get("Device", "auto"))
    if not device.startswith(_DEVICES_PREFIX):
        raise ValueError(f"{_ERR}Device must be auto, cpu, cuda or cuda:N")
    lr = doc.get("LearningRate", 1e-3)
    if isinstance(lr, bool) or not isinstance(lr, (int, float)) or lr <= 0:
        raise ValueError(f"{_ERR}LearningRate must be a positive number")
    notes = doc.get("Notes")
    if notes is not None and not isinstance(notes, str):
        raise ValueError(f"{_ERR}Notes must be a string")
    return TrainSpec(
        manifest=doc["Manifest"],
        fields=fields,
        target_fields=targets,
        run_dir=str(doc.get("RunDir", "runs/run")),
        train_split=str(doc.get("TrainSplit", "train")),
        valid_split=str(doc.get("ValidSplit", "valid")),
        epochs=_int(doc.get("Epochs", 100), "Epochs", 1),
        batch_size=_int(doc.get("BatchSize", 8), "BatchSize", 1),
        learning_rate=float(lr),
        seed=_int(doc.get("Seed", 0), "Seed"),
        checkpoint_every=_int(doc.get("CheckpointEvery", 10), "CheckpointEvery"),
        device=device,
        model_name=name,
        processor_size=_int(model.get("ProcessorSize", 8), "Model.ProcessorSize", 1),
        hidden_dim=_int(model.get("HiddenDim", 64), "Model.HiddenDim", 1),
        aggregation=aggregation,
        regions=_bool(graph.get("Regions", False), "Graph.Regions"),
        kind=kind,
        undirected=_bool(graph.get("Undirected", True), "Graph.Undirected"),
        edge_features=_bool(graph.get("EdgeFeatures", True), "Graph.EdgeFeatures"),
        float32=_bool(graph.get("Float32", True), "Graph.Float32"),
        target_offset=_int(graph.get("TargetOffset", 0), "Graph.TargetOffset"),
        target_delta=_bool(graph.get("TargetDelta", False), "Graph.TargetDelta"),
        read=dict(read),
        notes=notes,
        tags=_names(doc.get("Tags", ()), "Tags"),
        base_dir=base_dir,
    )


def spec_to_dict(spec: TrainSpec) -> dict:
    """The PascalCase document (round-trip exact through :func:`spec_from_dict`)."""
    doc = {
        "Version": SPEC_VERSION,
        "Manifest": spec.manifest,
        "RunDir": spec.run_dir,
        "Fields": list(spec.fields),
        "TargetFields": list(spec.target_fields),
        "TrainSplit": spec.train_split,
        "ValidSplit": spec.valid_split,
        "Epochs": spec.epochs,
        "BatchSize": spec.batch_size,
        "LearningRate": spec.learning_rate,
        "Seed": spec.seed,
        "CheckpointEvery": spec.checkpoint_every,
        "Device": spec.device,
        "Model": {
            "Name": spec.model_name,
            "ProcessorSize": spec.processor_size,
            "HiddenDim": spec.hidden_dim,
            "Aggregation": spec.aggregation,
        },
        "Graph": {
            "Regions": spec.regions,
            "Kind": spec.kind,
            "Undirected": spec.undirected,
            "EdgeFeatures": spec.edge_features,
            "Float32": spec.float32,
            "TargetOffset": spec.target_offset,
            "TargetDelta": spec.target_delta,
        },
    }
    if spec.read:
        doc["Read"] = dict(spec.read)
    if spec.notes is not None:
        doc["Notes"] = spec.notes
    if spec.tags:
        doc["Tags"] = list(spec.tags)
    return doc


def load_spec(settings) -> TrainSpec:
    """Load a spec from a dict, JSON text, or a file path (tri-modal, the
    ``DatasetManifest.load`` rule); a file's directory becomes ``base_dir``."""
    if isinstance(settings, TrainSpec):
        return settings
    if isinstance(settings, dict):
        return spec_from_dict(settings)
    text = str(settings)
    if text.lstrip()[:1] in ("{", "["):
        try:
            return spec_from_dict(json.loads(text))
        except json.JSONDecodeError as e:
            raise ValueError(f"{_ERR}the spec is not valid JSON: {e}") from None
    if not os.path.isfile(text):
        raise ValueError(f"{_ERR}spec file not found: '{text}'")
    with open(text, encoding="utf-8") as fh:
        try:
            doc = json.load(fh)
        except json.JSONDecodeError as e:
            raise ValueError(f"{_ERR}'{text}' is not valid JSON: {e}") from None
    return spec_from_dict(doc, base_dir=os.path.dirname(os.path.abspath(text)))


def save_spec(spec: TrainSpec, path) -> str:
    write_json_atomic(path, spec_to_dict(spec))
    return str(path)


def default_spec(manifest, fields, target_fields, *, run_dir="runs/run", **overrides):
    """A spec with every default filled in; ``overrides`` are field names."""
    spec = TrainSpec(
        manifest=str(manifest),
        fields=_names(fields, "Fields"),
        target_fields=_names(target_fields, "TargetFields"),
        run_dir=str(run_dir),
    )
    if overrides:
        unknown = [k for k in overrides if k not in spec.__dataclass_fields__]
        if unknown:
            raise ValueError(f"{_ERR}unknown spec field(s) {unknown}")
        spec = replace(spec, **overrides)
    # Re-validate through the document, so overrides obey the same rules.
    return spec_from_dict(spec_to_dict(spec), base_dir=spec.base_dir)


# --------------------------------------------------------------------------- #
# Run-directory files                                                         #
# --------------------------------------------------------------------------- #
def write_json_atomic(path, obj) -> None:
    """Write JSON via a sibling temp file + rename, so a reader (the job
    manager polling ``progress.json``) never sees a torn file."""
    path = str(path)
    directory = os.path.dirname(path) or "."
    os.makedirs(directory, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp-", dir=directory)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(obj, fh, indent=2)
            fh.write("\n")
        os.replace(tmp, path)
    except BaseException:
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


def read_json(path, default=None):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return default


def append_metrics(run_dir, row: dict) -> None:
    """Append one epoch's row to ``metrics.jsonl``."""
    with open(os.path.join(run_dir, METRICS_FILE), "a", encoding="utf-8") as fh:
        fh.write(json.dumps(row) + "\n")


def read_metrics(run_dir, since_epoch: int = 0) -> List[dict]:
    """Every metrics row with ``epoch >= since_epoch`` (a torn last line, from
    a trainer mid-write, is skipped rather than fatal)."""
    path = os.path.join(run_dir, METRICS_FILE)
    rows = []
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    row = json.loads(line)
                except ValueError:
                    continue
                if int(row.get("epoch", -1)) >= since_epoch:
                    rows.append(row)
    except OSError:
        return []
    return rows


def metrics_row(epoch, train_loss, valid_loss, lr, elapsed, epoch_seconds) -> dict:
    return {
        "epoch": int(epoch),
        "train_loss": float(train_loss),
        "valid_loss": None if valid_loss is None else float(valid_loss),
        "lr": float(lr),
        "elapsed": float(elapsed),
        "epoch_seconds": float(epoch_seconds),
        "timestamp": time.time(),
    }


# --------------------------------------------------------------------------- #
# Normalization stats and the model card                                      #
# --------------------------------------------------------------------------- #
def stat_vectors(stats: dict, names: Sequence[str], suffix: str) -> np.ndarray:
    """Concatenate ``{name}_{suffix}`` vectors (the ``node_stats.json``
    convention) for ``names``, in order; a ``std`` is floored at
    ``STATS_STD_FLOOR`` so a constant field never divides by zero."""
    values: List[float] = []
    for name in names:
        key = f"{name}_{suffix}"
        if key not in stats:
            raise KeyError(f"{_ERR}no '{key}' in the normalization stats")
        values.extend(float(v) for v in stats[key])
    out = np.asarray(values, dtype=np.float64)
    if suffix.endswith("std"):
        out = np.maximum(out, STATS_STD_FLOOR)
    return out


def card_from_run(
    spec: TrainSpec,
    schema: dict,
    node_stats: dict,
    edge_stats: dict,
    *,
    epoch: int,
    valid_loss: Optional[float],
    checkpoint: str,
    input_dims: Optional[Tuple[int, int, int]] = None,
) -> dict:
    """The model card for a checkpoint: what its channels mean and how they
    were normalized. ``input_dims`` = (nodes, edges, output)."""
    x_cols = list(schema["x_columns"])
    y_cols = list(schema["y_columns"])
    y_suffix = "diff_" if spec.target_delta else ""
    return {
        "version": 1,
        "meshioplusplus_version": __version__,
        "model": {
            "name": spec.model_name,
            "processor_size": spec.processor_size,
            "hidden_dim": spec.hidden_dim,
            "aggregation": spec.aggregation,
            "input_dim_nodes": input_dims[0] if input_dims else len(x_cols),
            "input_dim_edges": (
                input_dims[1] if input_dims else len(edge_stats.get("edge_mean", []))
            ),
            "output_dim": input_dims[2] if input_dims else len(y_cols),
        },
        "schema": schema,
        "x_columns": x_cols,
        "y_columns": y_cols,
        "fields": list(spec.fields),
        "target_fields": list(spec.target_fields),
        "graph": spec.graph_kwargs(),
        "read": dict(spec.read),
        "input_normalization": {
            "mean": stat_vectors(node_stats, spec.fields, "mean").tolist(),
            "std": stat_vectors(node_stats, spec.fields, "std").tolist(),
        },
        "output_normalization": {
            "mean": stat_vectors(
                node_stats, spec.target_fields, y_suffix + "mean"
            ).tolist(),
            "std": stat_vectors(
                node_stats, spec.target_fields, y_suffix + "std"
            ).tolist(),
        },
        "edge_normalization": {
            "mean": [float(v) for v in edge_stats.get("edge_mean", [])],
            "std": [
                max(float(v), STATS_STD_FLOOR) for v in edge_stats.get("edge_std", [])
            ],
        },
        "epoch": int(epoch),
        "valid_loss": None if valid_loss is None else float(valid_loss),
        "checkpoint": os.path.basename(checkpoint),
    }


def normalizers_from_card(card: dict) -> Dict[str, np.ndarray]:
    """``x_mean/x_std/y_mean/y_std/e_mean/e_std`` as float64 arrays."""
    out = {}
    for prefix, key in (
        ("x", "input_normalization"),
        ("y", "output_normalization"),
        ("e", "edge_normalization"),
    ):
        block = card[key]
        out[f"{prefix}_mean"] = np.asarray(block["mean"], dtype=np.float64)
        out[f"{prefix}_std"] = np.maximum(
            np.asarray(block["std"], dtype=np.float64), STATS_STD_FLOOR
        )
    return out


def card_path(checkpoint) -> str:
    return str(checkpoint) + CARD_SUFFIX


def list_checkpoints(run_dir, best: Optional[str] = None) -> List[dict]:
    """Every ``.mdlus`` under ``checkpoints/`` with its card's epoch and
    validation loss (never parsed from the file name), oldest epoch first;
    ``best`` names the one the job marks best."""
    directory = os.path.join(run_dir, CHECKPOINT_DIR)
    out = []
    for path in sorted(glob.glob(os.path.join(directory, "*.mdlus"))):
        name = os.path.basename(path)
        card = read_json(card_path(path), {}) or {}
        kind = (
            "best"
            if name == BEST_CHECKPOINT
            else "final" if name == FINAL_CHECKPOINT else "periodic"
        )
        out.append(
            {
                "path": path,
                "name": name,
                "kind": kind,
                "epoch": card.get("epoch"),
                "valid_loss": card.get("valid_loss"),
                "size": os.path.getsize(path),
                "is_best": best is not None and os.path.basename(best) == name,
            }
        )
    out.sort(key=lambda c: (c["epoch"] is None, c["epoch"] or 0, c["name"]))
    return out
