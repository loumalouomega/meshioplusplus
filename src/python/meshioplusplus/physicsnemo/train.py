"""Train PhysicsNeMo's MeshGraphNet from a dataset manifest, and predict
with the result: ``python -m meshioplusplus.physicsnemo.train --spec spec.json``.

The ONLY module in the package importing torch, PyTorch Geometric and
PhysicsNeMo's model/checkpoint utilities -- and only inside functions, so
``import meshioplusplus.physicsnemo`` stays pure. Reached through the lazy
``meshioplusplus.physicsnemo.run_training``/``predict`` factories, by the
job manager (``mcp/_jobs.py``) as a subprocess, and by the worked example.

What a run writes (``physicsnemo/_train.py`` owns the layout): one
``metrics.jsonl`` row and a ``progress.json`` update per epoch, a resumable
training checkpoint every ``CheckpointEvery`` epochs through
``physicsnemo.utils.checkpoint.save_checkpoint`` (model weights as ``.mdlus``
plus optimizer state), ``best.mdlus`` whenever the validation loss improves
and ``final.mdlus`` at the end -- every ``.mdlus`` with a **model card**
sidecar (what its channels mean and how they were normalized). SIGTERM is
honoured: the running epoch finishes, ``final.*`` is written, and the
process exits 143. The training path is PyG's (``make_dataset`` +
``DataLoader``): MeshGraphNet consumes PyG ``Data`` and PyG batches
variable-size graphs natively. See ``doc/physicsnemo.md``.
"""

from __future__ import annotations

import argparse
import functools
import glob
import math
import os
import signal
import sys
import time
import warnings
from typing import List, Optional

import numpy as np

from .. import DatasetManifest, write
from ..__about__ import __version__
from .._regions import block_bases
from . import (
    _DOC,
    _tessellate_for_sample,
    edge_stats,
    field_stats,
    graph_sample,
    make_dataset,
    operator_sample,
)
from ._train import (
    _FAMILY,
    BEST_CHECKPOINT,
    CHECKPOINT_DIR,
    EDGE_STATS_FILE,
    FINAL_CHECKPOINT,
    NODE_STATS_FILE,
    PROGRESS_FILE,
    SPEC_FILE,
    STATS_STD_FLOOR,
    TrainSpec,
    append_metrics,
    card_from_run,
    card_path,
    load_spec,
    metrics_row,
    normalizers_from_card,
    read_json,
    require_frameworks,
    save_spec,
    write_json_atomic,
)

_ERR = "meshio++: train: "
_STOP = {"requested": False}


def _frameworks(op, model_name="meshgraphnet"):
    """Require only what the chosen family actually imports -- the shared
    :func:`~meshioplusplus.physicsnemo._train.require_frameworks`, which reads
    the family table rather than a ``!= "srresnet"`` comparison."""
    require_frameworks(op, model_name, doc=_DOC)


def _family_block(model_name):
    """The spec block a family reads (``Graph``/``Grid``), from the table."""
    fam = _FAMILY.get(str(model_name).lower())
    if fam is None:
        raise ValueError(
            f"{_ERR}the card names model family {model_name!r}, which this build "
            f"does not know (known: {', '.join(_FAMILY)})"
        )
    return fam.block


def _device(name):
    import torch

    if name == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    return name


def build_model(spec: TrainSpec, schema: dict, *, input_dim_edges: int):
    """The MeshGraphNet a spec describes, sized from the recorded schema --
    never from ``len(Fields)``, which is wrong once region one-hots widen
    ``x``."""
    from physicsnemo.models.meshgraphnet import MeshGraphNet

    return MeshGraphNet(
        input_dim_nodes=len(schema["x_columns"]),
        input_dim_edges=input_dim_edges,
        output_dim=len(schema["y_columns"]),
        processor_size=spec.processor_size,
        hidden_dim_processor=spec.hidden_dim,
        hidden_dim_node_encoder=spec.hidden_dim,
        hidden_dim_edge_encoder=spec.hidden_dim,
        hidden_dim_node_decoder=spec.hidden_dim,
        aggregation=spec.aggregation,
    )


def _norm_tensors(card: dict, device):
    import torch

    norms = normalizers_from_card(card)
    return {
        k: torch.tensor(v, dtype=torch.float32, device=device) for k, v in norms.items()
    }


def _run_epoch(model, loader, norms, device, optimizer=None):
    import torch

    total, count = 0.0, 0
    for batch in loader:
        batch = batch.to(device)
        x = (batch.x - norms["x_mean"]) / norms["x_std"]
        edge_attr = (batch.edge_attr - norms["e_mean"]) / norms["e_std"]
        target = (batch.y - norms["y_mean"]) / norms["y_std"]
        pred = model(x, edge_attr, batch)
        loss = torch.nn.functional.mse_loss(pred, target)
        if optimizer is not None:
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        total += float(loss.detach()) * batch.num_graphs
        count += batch.num_graphs
    return total / max(count, 1) if count else None


def _write_card(
    spec, schema, node_stats, e_stats, path, *, epoch, valid_loss, dims, guard=None
):
    card = card_from_run(
        spec,
        schema,
        node_stats,
        e_stats,
        epoch=epoch,
        valid_loss=valid_loss,
        checkpoint=path,
        input_dims=dims,
    )
    if guard is not None:
        card["guard"] = guard.to_dict()
    write_json_atomic(card_path(path), card)


def _save_periodic(run_dir, model, optimizer, epoch, extra):
    """``physicsnemo.utils.checkpoint.save_checkpoint`` (weights + optimizer
    state, resumable) into ``checkpoints/``; returns the ``.mdlus`` it wrote."""
    from physicsnemo.utils.checkpoint import save_checkpoint

    try:  # save_checkpoint initializes it anyway, with a warning otherwise
        from physicsnemo.distributed import DistributedManager

        if not DistributedManager.is_initialized():
            DistributedManager.initialize()
    except Exception:  # noqa: BLE001 - best-effort; save_checkpoint copes
        pass
    directory = os.path.join(run_dir, CHECKPOINT_DIR)
    before = set(glob.glob(os.path.join(directory, "*.mdlus")))
    save_checkpoint(
        directory, models=model, optimizer=optimizer, epoch=epoch, metadata=extra
    )
    new = sorted(set(glob.glob(os.path.join(directory, "*.mdlus"))) - before)
    return new[-1] if new else None


def run(spec, *, log=print) -> dict:
    """Train per ``spec`` (a :class:`TrainSpec`, dict, JSON text or path);
    returns the final ``progress.json`` record."""
    spec = load_spec(spec)
    _frameworks("run_training", spec.model_name)
    return _RUNNERS[spec.model_name](spec, log=log)


def _run_graph(spec, *, log=print) -> dict:
    """Train a MeshGraphNet."""
    import torch
    from torch_geometric.loader import DataLoader

    if not spec.edge_features:
        raise ValueError(f"{_ERR}MeshGraphNet needs edge features (Graph.EdgeFeatures)")
    run_dir = spec.resolved_run_dir()
    manifest_path = spec.resolved_manifest()
    os.makedirs(os.path.join(run_dir, CHECKPOINT_DIR), exist_ok=True)
    if not os.path.isfile(os.path.join(run_dir, SPEC_FILE)):
        save_spec(spec, os.path.join(run_dir, SPEC_FILE))
    torch.manual_seed(spec.seed)
    np.random.seed(spec.seed)
    device = _device(spec.device)
    log(f"device: {device}")

    manifest = DatasetManifest.load(manifest_path)
    read_kwargs = dict(spec.read)
    fields = list(spec.fields)
    targets = list(spec.target_fields)
    node_stats = field_stats(
        manifest, split=spec.train_split, fields=fields + targets, **read_kwargs
    )
    if spec.target_delta:
        node_stats.update(
            field_stats(
                manifest,
                split=spec.train_split,
                fields=targets,
                delta=spec.target_offset or True,
                **read_kwargs,
            )
        )
    # The SAME edge options the samples will use: statistics gathered over a
    # different graph normalize the wrong thing, silently.
    e_stats = edge_stats(
        manifest,
        split=spec.train_split,
        kind=spec.kind,
        undirected=spec.undirected,
        proximity=spec.proximity,
        **read_kwargs,
    )
    write_json_atomic(os.path.join(run_dir, NODE_STATS_FILE), node_stats)
    write_json_atomic(os.path.join(run_dir, EDGE_STATS_FILE), e_stats)

    graph_kwargs = spec.graph_kwargs()
    guard = _fit_guard(spec, manifest, log)
    # Augmentation applies to the TRAIN split only: a validation loss that
    # moves because the poses moved measures nothing.
    augmentation = None
    if spec.augmentation:
        from ._augment import Augmentation

        augmentation = Augmentation.from_dict(spec.augmentation)
        log(f"augmentation: {augmentation.to_dict()}")
    train_ds = make_dataset(
        manifest,
        split=spec.train_split,
        augmentation=augmentation,
        **graph_kwargs,
        **read_kwargs,
    )
    valid_ds = make_dataset(
        manifest, split=spec.valid_split, **graph_kwargs, **read_kwargs
    )
    if len(train_ds) == 0:
        raise ValueError(
            f"{_ERR}split '{spec.train_split}' of {manifest_path} yields no samples"
        )
    log(f"train: {len(train_ds)} graphs, valid: {len(valid_ds)} graphs")
    schema = train_ds.schema
    input_dim_edges = len(e_stats["edge_mean"])
    dims = (len(schema["x_columns"]), input_dim_edges, len(schema["y_columns"]))
    model = build_model(spec, schema, input_dim_edges=input_dim_edges).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=spec.learning_rate)
    norms = _norm_tensors(
        card_from_run(
            spec,
            schema,
            node_stats,
            e_stats,
            epoch=0,
            valid_loss=None,
            checkpoint="",
            input_dims=dims,
        ),
        device,
    )
    train_loader = DataLoader(train_ds, batch_size=spec.batch_size, shuffle=True)
    valid_loader = (
        DataLoader(valid_ds, batch_size=spec.batch_size) if len(valid_ds) else None
    )

    def step(loader, optimizer_or_none):
        return _run_epoch(model, loader, norms, device, optimizer_or_none)

    def write_card(path, epoch, valid_loss):
        _write_card(
            spec,
            schema,
            node_stats,
            e_stats,
            path,
            epoch=epoch,
            valid_loss=valid_loss,
            dims=dims,
            guard=guard,
        )

    return _epoch_loop(
        spec,
        run_dir,
        device,
        model,
        optimizer,
        train_loader,
        valid_loader,
        step=step,
        write_card=write_card,
        log=log,
    )


def _build_srresnet(spec: TrainSpec, schema: dict):
    from physicsnemo.models.srrn import SRResNet

    return SRResNet(
        in_channels=len(schema["x_channels"]),
        out_channels=len(schema["y_channels"]),
        large_kernel_size=spec.large_kernel_size,
        small_kernel_size=spec.small_kernel_size,
        conv_layer_size=spec.conv_layer_size,
        n_resid_blocks=spec.resid_blocks,
        scaling_factor=spec.scaling_factor,
        activation_fn=spec.activation_fn,
    )


def _build_fno(spec: TrainSpec, schema: dict):
    """The FNO a spec describes; ``dimension`` IS the presence of a squeeze."""
    from physicsnemo.models.fno import FNO

    return FNO(
        in_channels=len(schema["x_channels"]),
        out_channels=len(schema["y_channels"]),
        decoder_layers=spec.fno_decoder_layers,
        decoder_layer_size=spec.fno_decoder_layer_size,
        decoder_activation_fn=spec.fno_decoder_activation_fn,
        dimension=2 if schema.get("squeeze") is not None else 3,
        latent_channels=spec.latent_channels,
        num_fno_layers=spec.num_fno_layers,
        num_fno_modes=spec.num_fno_modes,
        padding=spec.spectral_padding,
        padding_type=spec.padding_type,
        activation_fn=spec.fno_activation_fn,
        coord_features=spec.coord_features,
    )


def _build_afno(spec: TrainSpec, schema: dict):
    """The AFNO a spec describes, at the recorded 2-D sample shape.

    AFNO patches a FIXED ``(H, W)`` image, so its shape is taken from the
    recorded ``x_shape`` and the divisibility it needs is checked **before**
    torch is imported, naming the off-by-one everyone hits: ``Grid.Resolution
    [n, ...]`` gives ``n + 1`` sample points per axis.
    """
    inp_shape = [int(v) for v in schema["x_shape"]]
    patch = [int(v) for v in spec.patch_size]
    if len(inp_shape) != 2:
        raise ValueError(
            f"{_ERR}afno: the grid samples are {len(inp_shape)}-D ({inp_shape}); "
            "AFNO is 2-D only -- give Grid.Squeeze the world axis to collapse"
        )
    for n, p in zip(inp_shape, patch):
        if n % p:
            raise ValueError(
                f"{_ERR}afno: the grid's sample shape {inp_shape} is not divisible "
                f"by Model.PatchSize {patch}; note that Grid.Resolution [n, ...] "
                "gives n + 1 sample points per axis, so pick a resolution whose "
                "n + 1 is a multiple of the patch"
            )
    if spec.embed_dim % spec.num_blocks:
        raise ValueError(
            f"{_ERR}afno: Model.EmbedDim {spec.embed_dim} must be divisible by "
            f"Model.NumBlocks {spec.num_blocks}"
        )
    from physicsnemo.models.afno import AFNO

    return AFNO(
        inp_shape=inp_shape,
        in_channels=len(schema["x_channels"]),
        out_channels=len(schema["y_channels"]),
        patch_size=patch,
        embed_dim=spec.embed_dim,
        depth=spec.depth,
        mlp_ratio=spec.mlp_ratio,
        drop_rate=spec.drop_rate,
        num_blocks=spec.num_blocks,
        sparsity_threshold=spec.sparsity_threshold,
        hard_thresholding_fraction=spec.hard_thresholding_fraction,
    )


#: One constructor per grid family, keyed by name (never an else branch).
_GRID_BUILDERS = {"srresnet": _build_srresnet, "fno": _build_fno, "afno": _build_afno}


def build_grid_model(spec: TrainSpec, schema: dict):
    """The grid model a spec describes, sized from the recorded channel
    contract.

    ``in_channels``/``out_channels`` come from ``schema``, never from
    ``len(Fields)`` -- a multi-component array expands into one channel per
    component, so the two are different numbers the moment a vector field is
    involved. Same rule as :func:`build_model`, for the same reason.
    """
    return _GRID_BUILDERS[spec.model_name](spec, schema)


def _grid_model_block(spec: TrainSpec, schema: dict) -> dict:
    """The per-family ``"model"`` block of a grid card."""
    x_channels = list(schema["x_channels"])
    y_channels = list(schema["y_channels"])
    common = {
        "name": spec.model_name,
        "in_channels": len(x_channels),
        "out_channels": len(y_channels),
    }
    if spec.model_name == "srresnet":
        return {
            **common,
            "scaling_factor": spec.scaling_factor,
            "conv_layer_size": spec.conv_layer_size,
            "resid_blocks": spec.resid_blocks,
            "large_kernel_size": spec.large_kernel_size,
            "small_kernel_size": spec.small_kernel_size,
            "activation_fn": spec.activation_fn,
        }
    if spec.model_name == "fno":
        return {
            **common,
            "dimension": 2 if schema.get("squeeze") is not None else 3,
            "latent_channels": spec.latent_channels,
            "num_fno_layers": spec.num_fno_layers,
            "num_fno_modes": spec.num_fno_modes,
            "spectral_padding": spec.spectral_padding,
            "padding_type": spec.padding_type,
            "activation_fn": spec.fno_activation_fn,
            "decoder_layers": spec.fno_decoder_layers,
            "decoder_layer_size": spec.fno_decoder_layer_size,
            "decoder_activation_fn": spec.fno_decoder_activation_fn,
            "coord_features": spec.coord_features,
        }
    if spec.model_name == "afno":
        return {
            **common,
            "inp_shape": [int(v) for v in schema["x_shape"]],
            "patch_size": list(spec.patch_size),
            "embed_dim": spec.embed_dim,
            "depth": spec.depth,
            "mlp_ratio": spec.mlp_ratio,
            "drop_rate": spec.drop_rate,
            "num_blocks": spec.num_blocks,
            "sparsity_threshold": spec.sparsity_threshold,
            "hard_thresholding_fraction": spec.hard_thresholding_fraction,
        }
    raise ValueError(f"{_ERR}{spec.model_name!r} is not a grid family")


def _grid_norm_tensors(card, device):
    """Per-channel normalizers reshaped to broadcast over ``(B, C, *spatial)``.

    The spatial rank comes from the card's ``spatial_ndim`` (3, or 2 under a
    squeeze). It used to be hard-coded to three: against a ``(B, C, H, W)``
    batch a ``(1, C, 1, 1, 1)`` normalizer broadcasts along the BATCH axis with
    no error whenever ``C == 1`` -- so this is asserted on the shape, never on
    the arithmetic.
    """
    import torch

    norms = normalizers_from_card(card)
    ndim = int(card.get("spatial_ndim", 3))
    out = {}
    for key, values in norms.items():
        tensor = torch.as_tensor(values, dtype=torch.float32, device=device)
        # (C,) -> (1, C, 1, ..., 1): a channel statistic applies across the
        # whole volume, and the explicit reshape is what keeps that from
        # silently broadcasting along the wrong axis.
        out[key] = tensor.reshape((1, -1) + (1,) * ndim)
    return out


def _run_grid_epoch(model, loader, norms, device, optimizer=None):
    import torch

    total, count = 0.0, 0
    for x, y in loader:
        x = x.to(device)
        y = y.to(device)
        x = (x - norms["x_mean"]) / norms["x_std"]
        target = (y - norms["y_mean"]) / norms["y_std"]
        pred = model(x)
        loss = torch.nn.functional.mse_loss(pred, target)
        if optimizer is not None:
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        total += float(loss.detach()) * x.shape[0]
        count += x.shape[0]
    return total / max(count, 1) if count else None


def grid_card_from_run(
    spec: TrainSpec,
    schema: dict,
    stats: dict,
    *,
    epoch: int,
    valid_loss,
    checkpoint: str,
) -> dict:
    """The model card for a grid checkpoint.

    :func:`card_from_run`'s counterpart. It records the **grid construction**
    alongside the channel contract, and the array layout as a literal string:
    a checkpoint that does not say which axis is which produces finite,
    plausible, transposed numbers months later, and nothing downstream can tell.

    There is no edge block -- a convolutional model has no edges -- which
    :func:`normalizers_from_card` treats as empty rather than as an error.
    """
    from .._grid_transfer import grid_layout_after_squeeze

    x_channels = list(schema["x_channels"])
    y_channels = list(schema["y_channels"])
    squeeze = schema.get("squeeze")
    return {
        "version": 1,
        "meshioplusplus_version": __version__,
        "model": _grid_model_block(spec, schema),
        "schema": schema,
        # The layout string is a function of WHICH world axis was squeezed
        # (the remaining axes' order depends on it), so it is derived rather
        # than a constant; the shape contract beside it is what a 2-D family
        # builds its model from and writes its answer back through.
        "layout": grid_layout_after_squeeze(squeeze),
        "spatial_ndim": 3 if squeeze is None else 2,
        "squeeze": squeeze,
        "squeeze_index": schema.get("squeeze_index"),
        "x_shape": list(schema.get("x_shape") or []),
        "expand_axis": squeeze,
        "expand_size": schema.get("expand_size"),
        "x_columns": x_channels,
        "y_columns": y_channels,
        "x_channels": x_channels,
        "y_channels": y_channels,
        "fields": list(spec.fields),
        "target_fields": list(spec.target_fields),
        "grid": spec.grid_kwargs(),
        "coarse": schema.get("coarse"),
        "fine": schema.get("fine"),
        "read": dict(spec.read),
        # Informational: a checkpoint says what poses it was shown. Inference
        # never augments, so nothing reads this back.
        "augmentation": spec.augmentation,
        "input_normalization": {
            "mean": [float(v) for v in stats["x_mean"]],
            "std": [max(float(v), STATS_STD_FLOOR) for v in stats["x_std"]],
        },
        "output_normalization": {
            "mean": [float(v) for v in stats["y_mean"]],
            "std": [max(float(v), STATS_STD_FLOOR) for v in stats["y_std"]],
        },
        "epoch": int(epoch),
        "valid_loss": None if valid_loss is None else float(valid_loss),
        "checkpoint": os.path.basename(checkpoint),
    }


def _run_grid(spec, *, log=print) -> dict:
    """Train a grid family (srresnet/fno/afno) on a manifest's grid pairs --
    shared unchanged by all three; only ``build_grid_model`` differs."""
    import torch
    from torch.utils.data import DataLoader

    from . import grid_stats
    from ._torch import make_grid_dataset

    run_dir = spec.resolved_run_dir()
    manifest_path = spec.resolved_manifest()
    os.makedirs(os.path.join(run_dir, CHECKPOINT_DIR), exist_ok=True)
    if not os.path.isfile(os.path.join(run_dir, SPEC_FILE)):
        save_spec(spec, os.path.join(run_dir, SPEC_FILE))
    torch.manual_seed(spec.seed)
    np.random.seed(spec.seed)
    device = _device(spec.device)
    log(f"device: {device}")

    manifest = DatasetManifest.load(manifest_path)
    read_kwargs = dict(spec.read)
    guard = _fit_guard(spec, manifest, log)
    grid_kwargs = spec.grid_kwargs()

    stats = grid_stats(manifest, split=spec.train_split, **grid_kwargs, **read_kwargs)
    write_json_atomic(os.path.join(run_dir, NODE_STATS_FILE), stats)
    coverage = stats.get("coverage")
    if coverage is not None:
        log(f"grid coverage: {coverage:.4f} of points inside the mesh")
        if coverage < 0.5:
            log(
                "warning: more than half of every grid is fill -- a model trained "
                "on this learns the fill (tighten Grid.Bounds or the resolution)"
            )

    train_ds = make_grid_dataset(
        manifest, split=spec.train_split, read_kwargs=read_kwargs, **grid_kwargs
    )
    valid_ds = make_grid_dataset(
        manifest, split=spec.valid_split, read_kwargs=read_kwargs, **grid_kwargs
    )
    if len(train_ds) == 0:
        raise ValueError(
            f"{_ERR}split '{spec.train_split}' of {manifest_path} yields no samples"
        )
    log(f"train: {len(train_ds)} pairs, valid: {len(valid_ds)} pairs")

    schema = train_ds.schema
    model = build_grid_model(spec, schema).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=spec.learning_rate)
    norms = _grid_norm_tensors(
        grid_card_from_run(
            spec, schema, stats, epoch=0, valid_loss=None, checkpoint=""
        ),
        device,
    )
    train_loader = DataLoader(train_ds, batch_size=spec.batch_size, shuffle=True)
    valid_loader = (
        DataLoader(valid_ds, batch_size=spec.batch_size) if len(valid_ds) else None
    )

    def step(loader, optimizer_or_none):
        return _run_grid_epoch(model, loader, norms, device, optimizer_or_none)

    def write_card(path, epoch, valid_loss):
        card = grid_card_from_run(
            spec,
            schema,
            stats,
            epoch=epoch,
            valid_loss=valid_loss,
            checkpoint=path,
        )
        if guard is not None:
            card["guard"] = guard.to_dict()
        write_json_atomic(card_path(path), card)

    return _epoch_loop(
        spec,
        run_dir,
        device,
        model,
        optimizer,
        train_loader,
        valid_loader,
        step=step,
        write_card=write_card,
        log=log,
    )


def _epoch_loop(
    spec,
    run_dir,
    device,
    model,
    optimizer,
    train_loader,
    valid_loader,
    *,
    step,
    write_card,
    log,
):
    """The training loop's bookkeeping, shared by both model families.

    Everything in here is family-independent: per-epoch metrics, best-score
    tracking, ``progress.json``, periodic checkpoints and the SIGTERM handling
    that makes a stopped run leave a usable checkpoint. The **only** thing the
    families differ in is one epoch's forward pass, which arrives as ``step``.

    ``step(loader, optimizer_or_None)`` returns that epoch's mean loss;
    ``write_card(path, epoch, valid_loss)`` writes the card beside a checkpoint.
    """
    import torch

    _STOP["requested"] = False
    try:
        signal.signal(signal.SIGTERM, lambda *_: _STOP.update(requested=True))
    except ValueError:  # not the main thread (e.g. called from a worker)
        pass

    checkpoints = os.path.join(run_dir, CHECKPOINT_DIR)
    best = None
    best_path = None
    start = time.time()
    progress = {}
    for epoch in range(spec.epochs):
        t0 = time.time()
        # A dataset carrying an augmentation redraws its poses per epoch; the
        # loop is shared with the grid family, which has none, so this is a
        # capability check rather than a family branch.
        dataset = getattr(getattr(train_loader, "dataset", None), "set_epoch", None)
        if dataset is not None:
            dataset(epoch)
        model.train()
        train_loss = step(train_loader, optimizer)
        valid_loss = None
        if valid_loader is not None:
            model.eval()
            with torch.no_grad():
                valid_loss = step(valid_loader, None)
        now = time.time()
        append_metrics(
            run_dir,
            metrics_row(
                epoch, train_loss, valid_loss, spec.learning_rate, now - start, now - t0
            ),
        )
        score = valid_loss if valid_loss is not None else train_loss
        if best is None or score < best:
            best = score
            best_path = os.path.join(checkpoints, BEST_CHECKPOINT)
            model.save(best_path)
            write_card(best_path, epoch, valid_loss)
            best_epoch = epoch
        per_epoch = (now - start) / (epoch + 1)
        progress = {
            "epoch": epoch + 1,
            "epochs": spec.epochs,
            "best_epoch": best_epoch,
            "best_valid_loss": best,
            "best_checkpoint": best_path,
            "eta_seconds": per_epoch * (spec.epochs - epoch - 1),
            "device": device,
            "completed": False,
            "stopped": False,
        }
        write_json_atomic(os.path.join(run_dir, PROGRESS_FILE), progress)
        valid_text = f"{valid_loss:.3e}" if valid_loss is not None else "n/a"
        log(
            f"epoch {epoch:4d}  train {train_loss:.3e}  valid {valid_text}  ({now - t0:.1f} s)"
        )
        if spec.checkpoint_every and (epoch + 1) % spec.checkpoint_every == 0:
            written = _save_periodic(
                run_dir,
                model,
                optimizer,
                epoch,
                {"valid_loss": valid_loss, "train_loss": train_loss},
            )
            if written:
                write_card(written, epoch, valid_loss)
        if _STOP["requested"]:
            log("stop requested: finishing after this epoch")
            break

    final_path = os.path.join(checkpoints, FINAL_CHECKPOINT)
    model.save(final_path)
    write_card(
        final_path, progress.get("epoch", 0) - 1, progress.get("best_valid_loss")
    )
    progress["completed"] = not _STOP["requested"]
    progress["stopped"] = bool(_STOP["requested"])
    progress["elapsed"] = time.time() - start
    write_json_atomic(os.path.join(run_dir, PROGRESS_FILE), progress)
    log(
        f"{'stopped' if progress['stopped'] else 'trained'} after {progress['epoch']} epoch(s) in {progress['elapsed']:.1f} s; wrote {final_path}"
    )
    return progress


class _Loaded:
    """One checkpoint, loaded once: the model, its card, its normalizers.

    ``predict`` over a manifest pays for this once and reuses it per entry;
    ``predict_mesh`` pays for it per call, which is what a one-off inference
    is. The family comes from the card, so no caller has to say which
    ``predict`` it wants.
    """

    __slots__ = ("model", "card", "norms", "device", "family", "block")

    def __init__(self, model, card, norms, device, family, block):
        self.model = model
        self.card = card
        self.norms = norms
        self.device = device
        self.family = family
        #: the spec block the family reads -- what every per-mesh body and
        #: normalizer loader is keyed on, so a new family in a known block
        #: needs no new dispatch here
        self.block = block


def _load_checkpoint(checkpoint, device="auto") -> _Loaded:
    """Load a ``.mdlus`` and its card, and require only what its family needs."""
    checkpoint = str(checkpoint)
    card = read_json(card_path(checkpoint))
    if not isinstance(card, dict):
        raise ValueError(
            f"{_ERR}no model card beside '{checkpoint}' ({card_path(checkpoint)}); "
            "only checkpoints written by this trainer carry one"
        )
    family = card.get("model", {}).get("name", "meshgraphnet")
    block = _family_block(family)
    if block == "Grid":
        from .._grid_transfer import grid_layout_after_squeeze

        known = {grid_layout_after_squeeze(axis) for axis in (None, 0, 1, 2)}
        if card.get("layout") not in known | {None}:
            raise ValueError(
                f"{_ERR}the card records layout {card['layout']!r}, which this "
                "build does not know how to read"
            )
    _frameworks("predict", family)
    from physicsnemo.core.module import Module

    device = _device(device)
    with warnings.catch_warnings():
        # the experimental DeepONet announces itself on import; a checkpoint
        # load is not the place to hear it again
        warnings.simplefilter("ignore")
        model = Module.from_checkpoint(checkpoint).to(device)
    model.eval()
    norms = _NORM_LOADERS[block](card, device)
    return _Loaded(model, card, norms, device, family, block)


def _available_targets(mesh, card, location):
    """Whether this mesh actually carries the card's target fields.

    A mesh that was never in a manifest usually has no truth in it, and
    ``graph_sample`` would then quietly take ``y`` from the input's own step --
    an "error" against itself. Absent truth means no ``y``, no ``_error``
    arrays and a ``None`` rmse, which is the honest answer.
    """
    data = mesh.point_data if location == "point" else mesh.cell_data
    return all(name in data for name in card.get("target_fields", ()))


def _predict_graph_mesh(loaded, mesh, target_mesh=None, label="mesh"):
    """One mesh through a graph checkpoint -> ``(mesh, row)``; no file I/O."""
    import torch
    from torch_geometric.data import Data

    card, norms, device = loaded.card, loaded.norms, loaded.device
    graph_kwargs = dict(card["graph"])
    kind = graph_kwargs.get("kind", "node")
    location = "point" if kind == "node" else "cell"
    y_columns = list(card["y_columns"])
    if target_mesh is None and (
        graph_kwargs.get("target_offset", 0)
        or not _available_targets(mesh, card, location)
    ):
        graph_kwargs["target_fields"] = None
        graph_kwargs["target_delta"] = False

    # When tessellation was used at training time (recorded in the card),
    # rebuild it on THIS mesh before sampling -- the same isoparametric
    # linearization the checkpoint's own x_columns/edge structure expect.
    # The ORIGINAL mesh is what gets returned (and what a prediction is
    # written onto): `sample_mesh` is only the tessellated stand-in
    # `graph_sample` reads, never the caller's own mesh object.
    orig_mesh = mesh
    tess = None
    tess_opt = graph_kwargs.pop("tessellate", None)
    forward_kwargs = graph_kwargs
    sample_mesh = mesh
    sample_target = target_mesh
    if tess_opt:
        sample_mesh, tess = _tessellate_for_sample(
            mesh, tess_opt, forward_kwargs.get("fields")
        )
        if target_mesh is not None:
            sample_target, _ = _tessellate_for_sample(
                target_mesh, tess_opt, forward_kwargs.get("target_fields")
            )

    sample = graph_sample(sample_mesh, target_mesh=sample_target, **forward_kwargs)
    if list(sample.x_columns) != list(card["x_columns"]):
        raise ValueError(
            f"{_ERR}feature drift: the checkpoint was trained on x columns "
            f"{card['x_columns']} but {label} yields {list(sample.x_columns)}"
        )
    arrays = sample.arrays
    with torch.no_grad():
        x = (torch.from_numpy(arrays["x"]).to(device) - norms["x_mean"]) / norms[
            "x_std"
        ]
        edge_attr = (
            torch.from_numpy(arrays["edge_attr"]).to(device) - norms["e_mean"]
        ) / norms["e_std"]
        graph = Data(
            x=x,
            edge_index=torch.from_numpy(arrays["edge_index"]).to(device),
            edge_attr=edge_attr,
        )
        pred_norm = loaded.model(x, edge_attr, graph)
        pred = (
            (pred_norm * norms["y_std"] + norms["y_mean"])
            .cpu()
            .numpy()
            .astype(np.float64)
        )
    truth = arrays.get("y")
    truth = None if truth is None else np.asarray(truth, dtype=np.float64)
    for i, column in enumerate(y_columns):
        _attach(orig_mesh, kind, f"{column}_pred", pred[:, i], tess=tess)
        if truth is not None:
            _attach(
                orig_mesh,
                kind,
                f"{column}_error",
                np.abs(pred[:, i] - truth[:, i]),
                tess=tess,
            )
    mesh = orig_mesh
    error = None if truth is None else np.abs(pred - truth)
    row = {
        "num_rows": int(pred.shape[0]),
        "rmse": (None if error is None else float(math.sqrt(float(np.mean(error**2))))),
        "max_error": None if error is None else float(error.max()),
    }
    _check_guard(card, mesh, label, row)
    return mesh, row


def _predict_grid_mesh(loaded, mesh, target_mesh=None, label="mesh"):
    """One mesh through a grid checkpoint -> ``(mesh, row)``; no file I/O."""
    import torch

    from .._grid_transfer import GridArray, GridSpec, expand_grid, scatter_grid
    from . import grid_sample_pair

    card, norms, device = loaded.card, loaded.norms, loaded.device
    grid_kwargs = dict(card["grid"])
    x_channels = list(card["x_channels"])
    y_channels = list(card["y_channels"])
    if target_mesh is None and not _available_targets(mesh, card, "point"):
        grid_kwargs["target_fields"] = None

    sample = grid_sample_pair(mesh, target_mesh=target_mesh, **grid_kwargs)
    if list(sample.x_channels) != x_channels:
        raise ValueError(
            f"{_ERR}{label} produces channels {list(sample.x_channels)} but the "
            f"checkpoint was trained on {x_channels}"
        )
    family = card.get("model", {}).get("name", "srresnet")
    x_np = sample.arrays["x"]
    inp_shape = card.get("model", {}).get("inp_shape")
    if inp_shape is not None and tuple(x_np.shape[1:]) != tuple(inp_shape):
        # AFNO's resolution is fixed at construction; fail by name here rather
        # than inside its patch embedding.
        raise ValueError(
            f"{_ERR}{label} samples to shape {tuple(x_np.shape[1:])} but this "
            f"'{family}' checkpoint was built for {tuple(inp_shape)} (AFNO patches "
            "a fixed image; use the training grid)"
        )
    x = torch.from_numpy(x_np).float().unsqueeze(0).to(device)
    with torch.no_grad():
        pred = loaded.model((x - norms["x_mean"]) / norms["x_std"])
    pred = (pred * norms["y_std"] + norms["y_mean"])[0].cpu().numpy()

    truth = sample.arrays.get("y")
    if truth is not None and pred.shape != truth.shape:
        reason = (
            "the checkpoint's scaling factor does not match the pair this spec "
            "builds"
            if family == "srresnet"
            else f"a '{family}' model is resolution-preserving, so its output "
            "follows the input grid and the target grid must match it"
        )
        raise ValueError(
            f"{_ERR}{label}: the model emitted {pred.shape} but the target grid "
            f"is {truth.shape}; {reason}"
        )
    squeeze = card.get("squeeze")
    if squeeze is not None:
        # A 2-D prediction is written back over every plane of the thin axis.
        pred = expand_grid(pred, card["expand_axis"], card["expand_size"])
        if truth is not None:
            truth = expand_grid(truth, card["expand_axis"], card["expand_size"])
    fine = GridSpec.from_dict(sample.schema["fine"])
    out_mesh = scatter_grid(
        GridArray(pred, fine, tuple(y_channels)),
        target_mesh if target_mesh is not None else mesh,
        names=[f"{c}_pred" for c in y_channels],
    )
    row = {
        "num_rows": int(pred.reshape(pred.shape[0], -1).shape[1]),
        "coverage": sample.schema.get("y_coverage"),
    }
    if truth is not None:
        error = pred - truth
        row["rmse"] = float(np.sqrt(np.mean(np.square(error))))
        row["max_error"] = float(np.max(np.abs(error)))
        # A 2-D power spectrum is a documented follow-up: reported as None
        # under a squeeze rather than computed over duplicated planes.
        row["spectrum_rel_l2"] = (
            None if squeeze is not None else _spectrum_rel_l2(pred, truth, fine)
        )
        truth_mesh = scatter_grid(
            GridArray(truth, fine, tuple(y_channels)),
            out_mesh,
            names=[f"{c}_true" for c in y_channels],
            on_conflict="overwrite",
        )
        out_mesh = scatter_grid(
            GridArray(error, fine, tuple(y_channels)),
            truth_mesh,
            names=[f"{c}_error" for c in y_channels],
            on_conflict="overwrite",
        )
    else:
        row["rmse"] = None
        row["max_error"] = None
    _check_guard(card, mesh, label, row)
    return out_mesh, row


# --------------------------------------------------------------------------- #
# the deeponet family: parameters in, field out                               #
# --------------------------------------------------------------------------- #
def build_operator_model(spec: TrainSpec, schema: dict):
    """The DeepONet a spec describes: an MLP **branch** over the recorded
    parameter columns, an MLP **trunk** over the 3 coordinates, both to
    ``Width``, combined by the installed experimental ``DeepONet``.

    Sized from the recorded schema, never from ``len(Operator.Parameters)``: a
    length-``k`` parameter expands into ``k`` columns.
    """
    from physicsnemo.models.mlp.fully_connected import FullyConnected

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from physicsnemo.experimental.models.xdeeponet import DeepONet

    branch = FullyConnected(
        in_features=len(schema["parameter_columns"]),
        layer_size=spec.branch_layer_size,
        out_features=spec.width,
        num_layers=spec.branch_layers,
        activation_fn=spec.deeponet_activation_fn,
    )
    trunk = FullyConnected(
        in_features=3,
        layer_size=spec.trunk_layer_size,
        out_features=spec.width,
        num_layers=spec.trunk_layers,
        activation_fn=spec.deeponet_activation_fn,
    )
    return DeepONet(
        branch,
        trunk=trunk,
        width=spec.width,
        out_channels=len(schema["y_columns"]),
        decoder_type=spec.decoder_type,
        decoder_width=spec.decoder_width,
        decoder_layers=spec.decoder_layers,
        decoder_activation_fn=spec.decoder_activation_fn,
    )


def _operator_norm_tensors(card, device):
    """The three normalizers of an operator card, shaped for their tensors:
    parameters ``(1, p)`` against a ``(B, p)`` batch, the trunk ``(1, 3)``
    against ``(T, 3)``, the output ``(1, 1, k)`` against ``(B, T, k)``."""
    import torch

    norms = normalizers_from_card(card)
    trunk = card.get("trunk_normalization") or {"mean": [0.0] * 3, "std": [1.0] * 3}

    def tensor(values, shape):
        return torch.as_tensor(
            np.asarray(values, dtype=np.float64), dtype=torch.float32, device=device
        ).reshape(shape)

    return {
        "x_mean": tensor(norms["x_mean"], (1, -1)),
        "x_std": tensor(norms["x_std"], (1, -1)),
        "t_mean": tensor(trunk["mean"], (1, -1)),
        "t_std": tensor(
            np.maximum(np.asarray(trunk["std"], dtype=np.float64), STATS_STD_FLOOR),
            (1, -1),
        ),
        "y_mean": tensor(norms["y_mean"], (1, 1, -1)),
        "y_std": tensor(norms["y_std"], (1, 1, -1)),
    }


def _run_operator_epoch(model, loader, norms, trunk_norm, device, optimizer=None):
    import torch

    total, count = 0.0, 0
    for params, y in loader:
        params = params.to(device)
        y = y.to(device)
        x = (params - norms["x_mean"]) / norms["x_std"]
        target = (y - norms["y_mean"]) / norms["y_std"]
        pred = model(x, trunk_norm)
        loss = torch.nn.functional.mse_loss(pred, target)
        if optimizer is not None:
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        total += float(loss.detach()) * params.shape[0]
        count += params.shape[0]
    return total / max(count, 1) if count else None


def operator_card_from_run(
    spec: TrainSpec,
    schema: dict,
    stats: dict,
    trunk_indices,
    *,
    epoch: int,
    valid_loss,
    checkpoint: str,
) -> dict:
    """The model card for a deeponet checkpoint.

    Records the parameter contract (``parameters`` -- the expanded columns
    the branch was trained on -- and ``parameter_names``, the ``Metadata``
    keys they came from), the trunk construction, and **``trunk_indices``
    when the trunk is a budget**: the trunk is rebuilt from the mesh at
    predict time, so the indices are what detects a budget that selected
    different points (a JSON card is no place for a ``(T, 3)`` array). Three
    normalizations, one per tensor.
    """
    columns = list(schema["parameter_columns"])
    y_columns = list(schema["y_columns"])
    return {
        "version": 1,
        "meshioplusplus_version": __version__,
        "model": {
            "name": spec.model_name,
            "branch_layers": spec.branch_layers,
            "branch_layer_size": spec.branch_layer_size,
            "trunk_layers": spec.trunk_layers,
            "trunk_layer_size": spec.trunk_layer_size,
            "width": spec.width,
            "decoder_type": spec.decoder_type,
            "decoder_width": spec.decoder_width,
            "decoder_layers": spec.decoder_layers,
            "decoder_activation_fn": spec.decoder_activation_fn,
            "activation_fn": spec.deeponet_activation_fn,
            "in_parameters": len(columns),
            "out_channels": len(y_columns),
        },
        "schema": schema,
        "x_columns": columns,
        "y_columns": y_columns,
        "parameters": columns,
        "parameter_names": list(spec.parameters),
        "fields": [],
        "target_fields": list(spec.target_fields),
        "operator": spec.operator_kwargs(),
        "trunk": spec.trunk,
        "trunk_count": spec.trunk_count,
        "trunk_method": spec.trunk_method,
        "trunk_seed": spec.trunk_seed,
        "trunk_indices": (
            None if trunk_indices is None else [int(i) for i in trunk_indices]
        ),
        "num_points": int(schema["num_points"]),
        "read": dict(spec.read),
        "augmentation": None,
        "input_normalization": {
            "mean": [float(v) for v in stats["params_mean"]],
            "std": [max(float(v), STATS_STD_FLOOR) for v in stats["params_std"]],
        },
        "trunk_normalization": {
            "mean": [float(v) for v in stats["trunk_mean"]],
            "std": [max(float(v), STATS_STD_FLOOR) for v in stats["trunk_std"]],
        },
        "output_normalization": {
            "mean": [float(v) for v in stats["y_mean"]],
            "std": [max(float(v), STATS_STD_FLOOR) for v in stats["y_std"]],
        },
        "epoch": int(epoch),
        "valid_loss": None if valid_loss is None else float(valid_loss),
        "checkpoint": os.path.basename(checkpoint),
    }


def _run_operator(spec, *, log=print) -> dict:
    """Train a DeepONet on a manifest's cases: each entry's ``Metadata``
    parameters in, its target fields at the (shared) trunk points out."""
    import torch
    from torch.utils.data import DataLoader

    from . import operator_stats
    from ._torch import make_operator_dataset

    run_dir = spec.resolved_run_dir()
    manifest_path = spec.resolved_manifest()
    os.makedirs(os.path.join(run_dir, CHECKPOINT_DIR), exist_ok=True)
    if not os.path.isfile(os.path.join(run_dir, SPEC_FILE)):
        save_spec(spec, os.path.join(run_dir, SPEC_FILE))
    torch.manual_seed(spec.seed)
    np.random.seed(spec.seed)
    device = _device(spec.device)
    log(f"device: {device}")

    manifest = DatasetManifest.load(manifest_path)
    read_kwargs = dict(spec.read)
    guard = _fit_guard(spec, manifest, log)
    op_kwargs = spec.operator_kwargs()

    stats = operator_stats(manifest, split=spec.train_split, **op_kwargs, **read_kwargs)
    write_json_atomic(os.path.join(run_dir, NODE_STATS_FILE), stats)

    train_ds = make_operator_dataset(
        manifest, split=spec.train_split, read_kwargs=read_kwargs, **op_kwargs
    )
    valid_ds = make_operator_dataset(
        manifest, split=spec.valid_split, read_kwargs=read_kwargs, **op_kwargs
    )
    if len(train_ds) == 0:
        raise ValueError(
            f"{_ERR}split '{spec.train_split}' of {manifest_path} yields no samples"
        )
    log(f"train: {len(train_ds)} cases, valid: {len(valid_ds)} cases")

    schema = train_ds.schema
    trunk_indices = train_ds.trunk_indices
    log(
        f"trunk: {schema['num_trunk']} of {schema['num_points']} points "
        f"({spec.trunk}), {len(schema['parameter_columns'])} parameter column(s)"
    )
    model = build_operator_model(spec, schema).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=spec.learning_rate)
    norms = _operator_norm_tensors(
        operator_card_from_run(
            spec, schema, stats, trunk_indices, epoch=0, valid_loss=None, checkpoint=""
        ),
        device,
    )
    trunk_norm = (train_ds.trunk.to(device) - norms["t_mean"]) / norms["t_std"]
    train_loader = DataLoader(train_ds, batch_size=spec.batch_size, shuffle=True)
    valid_loader = (
        DataLoader(valid_ds, batch_size=spec.batch_size) if len(valid_ds) else None
    )

    def step(loader, optimizer_or_none):
        return _run_operator_epoch(
            model, loader, norms, trunk_norm, device, optimizer_or_none
        )

    def write_card(path, epoch, valid_loss):
        card = operator_card_from_run(
            spec,
            schema,
            stats,
            trunk_indices,
            epoch=epoch,
            valid_loss=valid_loss,
            checkpoint=path,
        )
        if guard is not None:
            card["guard"] = guard.to_dict()
        write_json_atomic(card_path(path), card)

    return _epoch_loop(
        spec,
        run_dir,
        device,
        model,
        optimizer,
        train_loader,
        valid_loader,
        step=step,
        write_card=write_card,
        log=log,
    )


def _predict_operator_mesh(
    loaded, mesh, target_mesh=None, label="mesh", parameters=None
):
    """One case through a deeponet checkpoint -> ``(mesh, row)``; no file I/O.

    ``parameters`` is the case's ``Metadata``-shaped dict (required -- the
    branch has nothing to read from a mesh). The trunk is rebuilt from the
    mesh and checked against the card: the point count must match (fixed
    geometry, named), and a budget must select the SAME points it did at
    training time. Under a budget the write-back fills **NaN, never 0**, at
    the points the model never saw.
    """
    import torch

    card, norms, device = loaded.card, loaded.norms, loaded.device
    family = card.get("model", {}).get("name", "deeponet")
    names = list(card.get("parameter_names", []))
    if parameters is None:
        raise ValueError(
            f"{_ERR}parameters= is required for a '{family}' checkpoint (this one "
            f"expects {names})"
        )
    op_kwargs = dict(card["operator"])
    y_columns = list(card["y_columns"])
    y_mesh = target_mesh if target_mesh is not None else mesh
    if not _available_targets(y_mesh, card, "point"):
        op_kwargs["target_fields"] = None
    sample = operator_sample(y_mesh, parameters, **op_kwargs)
    if list(sample.parameter_columns) != list(card["parameters"]):
        raise ValueError(
            f"{_ERR}parameter drift: the checkpoint was trained on parameter "
            f"columns {card['parameters']} but {label} yields "
            f"{list(sample.parameter_columns)}"
        )
    num_points = int(card["num_points"])
    if sample.schema["num_points"] != num_points:
        raise ValueError(
            f"{_ERR}{label} has {sample.schema['num_points']} points but the "
            f"checkpoint's trunk was trained on {num_points}; a DeepONet is "
            "fixed-geometry, so predict on the training geometry (or use "
            "meshgraphnet for varying geometry)"
        )
    indices = sample.arrays.get("trunk_indices")
    recorded = card.get("trunk_indices")
    if recorded is not None and (
        indices is None or [int(i) for i in indices] != list(recorded)
    ):
        raise ValueError(
            f"{_ERR}{label}: the token budget selected different points than at "
            "training time, so the trunk no longer matches the checkpoint; the "
            "geometry has changed"
        )

    params = torch.from_numpy(np.asarray(sample.arrays["params"], dtype=np.float32))
    trunk = torch.from_numpy(np.asarray(sample.arrays["trunk"], dtype=np.float32))
    with torch.no_grad():
        x = (params.unsqueeze(0).to(device) - norms["x_mean"]) / norms["x_std"]
        t = (trunk.to(device) - norms["t_mean"]) / norms["t_std"]
        pred_norm = loaded.model(x, t)
        pred = (
            (pred_norm * norms["y_std"] + norms["y_mean"])[0]
            .cpu()
            .numpy()
            .astype(np.float64)
        )
    truth = sample.arrays.get("y")
    truth = None if truth is None else np.asarray(truth, dtype=np.float64)

    def whole(values):
        # every point of the mesh: the trunk's rows in place, NaN elsewhere
        if indices is None:
            return values
        out = np.full(num_points, np.nan, dtype=np.float64)
        out[np.asarray(indices, dtype=np.int64)] = values
        return out

    for i, column in enumerate(y_columns):
        _attach(mesh, "node", f"{column}_pred", whole(pred[:, i]))
        if truth is not None:
            _attach(
                mesh, "node", f"{column}_error", whole(np.abs(pred[:, i] - truth[:, i]))
            )
    error = None if truth is None else np.abs(pred - truth)
    row = {
        "num_rows": int(pred.shape[0]),
        "rmse": None if error is None else float(math.sqrt(float(np.mean(error**2)))),
        "max_error": None if error is None else float(error.max()),
    }
    _check_guard(card, mesh, label, row)
    return mesh, row


def predict_operator(
    checkpoint,
    manifest,
    *,
    entry_ids: Optional[List[str]] = None,
    split: Optional[str] = "test",
    step: int = 0,
    output_dir,
    device: str = "auto",
) -> List[dict]:
    """Predict over a manifest's entries with a deeponet checkpoint, each
    entry's ``Metadata`` supplying its parameters -- :func:`predict`'s
    operator counterpart, reached through it."""
    loaded = _load_checkpoint(checkpoint, device)
    read_kwargs = dict(loaded.card.get("read", {}))
    manifest = DatasetManifest.load(manifest)
    entries = list(manifest.entries(split=split))
    if entry_ids:
        wanted = set(entry_ids)
        entries = [e for e in entries if e.id in wanted]
    os.makedirs(output_dir, exist_ok=True)

    rows = []
    for entry in entries:
        series = entry.time_series(**read_kwargs)
        time_value, mesh = series[step]
        mesh, row = _predict_operator_mesh(
            loaded, mesh, None, f"entry '{entry.id}'", parameters=entry.metadata
        )
        suffix = "" if step == 0 else f"_s{step}"
        out_path = os.path.join(output_dir, f"{entry.id}{suffix}.vtu")
        write(out_path, mesh)
        row["entry_id"] = entry.id
        row["time"] = time_value
        row["output_path"] = out_path
        rows.append(row)
        del mesh
    return rows


def _check_guard(card, mesh, label, row):
    """Score a mesh against the card's guardrail, if it carries one.

    Advisory: the prediction is made and written either way. A model cannot
    refuse to answer, so the honest thing is to answer and say loudly that the
    input is unlike anything it was trained on.
    """
    doc = card.get("guard")
    if not doc:
        return
    from .._guard import GeometryGuard

    try:
        report = GeometryGuard.from_dict(doc).check(mesh)
    except Exception as exc:  # noqa: BLE001 -- advisory, never fatal
        warnings.warn(
            f"{_ERR}guard: could not score {label} ({type(exc).__name__})",
            stacklevel=2,
        )
        return
    row["guard"] = report
    if report["verdict"] == "out":
        worst = report["worst"][0]["name"] if report["worst"] else "?"
        warnings.warn(
            f"{_ERR}{label} is out of the distribution this checkpoint was "
            f"trained on (score {report['score']:.3g} > threshold "
            f"{report['threshold']:.3g}, worst descriptor '{worst}'); the "
            "prediction is still written",
            stacklevel=2,
        )


def _body_for(loaded, parameters):
    """The per-mesh body for a loaded checkpoint, with ``parameters`` bound
    for an operator family and refused by name for the others."""
    body = _PREDICT_BODIES[loaded.block]
    if loaded.block == "Operator":
        return functools.partial(body, parameters=parameters)
    if parameters is not None:
        raise ValueError(
            f"{_ERR}parameters= applies to a 'deeponet' checkpoint; this one is a "
            f"'{loaded.family}'"
        )
    return body


def predict_mesh(
    checkpoint,
    mesh,
    *,
    target_mesh=None,
    device="auto",
    label="mesh",
    parameters=None,
):
    """Predict with a checkpoint on ONE in-memory mesh -> ``(mesh, row)``.

    The card says which family wrote the checkpoint, so a caller does not have
    to. Everything the prediction needs -- the sample options, the read
    options, the column contract, the normalization -- comes from the card;
    nothing here consults a manifest. A mesh carrying no truth predicts anyway,
    with ``rmse``/``max_error`` reported as ``None`` rather than measured
    against itself. A ``deeponet`` checkpoint needs ``parameters`` (the case's
    ``Metadata``-shaped dict), the one thing a mesh cannot supply.
    """
    loaded = _load_checkpoint(checkpoint, device)
    return _body_for(loaded, parameters)(loaded, mesh, target_mesh, label)


def predict_file(
    checkpoint,
    input_path,
    output_path,
    *,
    time_step: Optional[int] = None,
    target_path=None,
    input_format: Optional[str] = None,
    output_format: Optional[str] = None,
    device: str = "auto",
    parameters=None,
) -> dict:
    """Predict with a checkpoint on ONE mesh file, writing the result.

    The single-mesh counterpart of :func:`predict`: no manifest, no split, no
    entry -- point a trained checkpoint at a file that was never catalogued
    anywhere and get the prediction written back as ordinary data arrays.

    ``time_step`` selects a step of a multi-step file. ``target_path`` supplies
    the paired mesh a t->t+n or coarse/fine checkpoint compares against; with
    no target and no truth in the file the prediction is still written, and
    ``rmse``/``max_error`` come back ``None``.
    """
    from .._sequence import TimeSeries

    loaded = _load_checkpoint(checkpoint, device)
    read_kwargs = dict(loaded.card.get("read", {}))
    if input_format is not None:
        read_kwargs["file_format"] = input_format
    series = TimeSeries(input_path, **read_kwargs)
    step = 0 if time_step is None else int(time_step)
    if step < 0:
        step += len(series)
    if not 0 <= step < len(series):
        raise ValueError(
            f"{_ERR}'{input_path}' has {len(series)} step(s); step {time_step} "
            "is out of range"
        )
    time_value, mesh = series[step]

    target_mesh = None
    if target_path is not None:
        target_mesh = TimeSeries(target_path, **read_kwargs)[step][1]
    elif loaded.block == "Graph":
        offset = int(dict(loaded.card["graph"]).get("target_offset", 0))
        if offset and step + offset < len(series):
            target_mesh = series[step + offset][1]

    label = f"'{input_path}'"
    out, row = _body_for(loaded, parameters)(loaded, mesh, target_mesh, label)
    write(output_path, out, file_format=output_format)
    row["input_path"] = str(input_path)
    row["output_path"] = str(output_path)
    row["time"] = None if time_value is None else float(time_value)
    return row


def _fit_guard(spec, manifest, log=print):
    """Fit the geometry guardrail on the training split, or ``None``.

    Advisory by construction: a fit that fails is reported and the run goes
    on, because a guardrail is a warning system and refusing to train for
    want of one would be the wrong trade.
    """
    if not spec.guard and spec.guard != {}:
        return None
    from .._guard import GeometryGuard

    options = dict(spec.guard or {})
    try:
        fitted = GeometryGuard.fit(
            manifest, split=spec.train_split, **options, **dict(spec.read)
        )
    except Exception as exc:  # noqa: BLE001 -- advisory, never fatal
        log(f"guard: could not fit ({type(exc).__name__}: {exc}); continuing")
        return None
    log(
        f"guard: fitted over {fitted.schema.get('num_samples')} mesh(es), "
        f"threshold {fitted.threshold:.4g}"
    )
    return fitted


def predict(
    checkpoint,
    manifest,
    *,
    entry_ids: Optional[List[str]] = None,
    split: Optional[str] = "test",
    step: int = 0,
    output_dir,
    device: str = "auto",
) -> List[dict]:
    """Predict with a ``.mdlus`` checkpoint (plus its card) over a manifest's
    entries and write each prediction back as ordinary data arrays --
    ``<column>_pred`` and, where the truth is present, ``<column>_error`` --
    into ``output_dir/<entry_id>.vtu``. Returns one row per entry.

    A loop over :func:`predict_mesh`'s own per-mesh body: the manifest supplies
    which meshes to read and nothing else, which is why
    :func:`predict_file` can do the same job for a file that was never
    catalogued."""
    loaded = _load_checkpoint(checkpoint, device)
    if loaded.block in _MANIFEST_PREDICTORS:
        return _MANIFEST_PREDICTORS[loaded.block](
            checkpoint,
            manifest,
            entry_ids=entry_ids,
            split=split,
            step=step,
            output_dir=output_dir,
            device=device,
        )
    card = loaded.card
    read_kwargs = dict(card.get("read", {}))
    offset = int(dict(card["graph"]).get("target_offset", 0))

    if not isinstance(manifest, DatasetManifest):
        manifest = DatasetManifest.load(manifest)
    entries = manifest.entries(split=split) if split else list(manifest)
    if entry_ids:
        wanted = list(entry_ids)
        known = {e.id for e in entries}
        unknown = [i for i in wanted if i not in known]
        if unknown:
            raise ValueError(f"{_ERR}unknown entry id(s): {unknown}")
        entries = [e for e in entries if e.id in set(wanted)]
    os.makedirs(output_dir, exist_ok=True)

    rows = []
    for entry in entries:
        series = entry.time_series(**read_kwargs)
        if step + offset >= len(series):
            raise ValueError(
                f"{_ERR}entry '{entry.id}' has {len(series)} step(s); step {step}"
                + (f" + offset {offset}" if offset else "")
                + " is out of range"
            )
        time_value, mesh = series[step]
        target_mesh = series[step + offset][1] if offset else None
        mesh, row = _predict_graph_mesh(
            loaded, mesh, target_mesh, f"entry '{entry.id}'"
        )
        name = f"{entry.id}.vtu" if step == 0 else f"{entry.id}_s{step}.vtu"
        out_path = os.path.join(output_dir, name)
        write(out_path, mesh)
        row["entry_id"] = entry.id
        row["time"] = float(time_value)
        row["output_path"] = out_path
        rows.append(row)
        del mesh, target_mesh
    return rows


def predict_grid(
    checkpoint,
    manifest,
    *,
    entry_ids: Optional[List[str]] = None,
    split: Optional[str] = "test",
    step: int = 0,
    output_dir,
    device: str = "auto",
) -> List[dict]:
    """Superresolve a manifest's entries with a grid checkpoint.

    :func:`predict`'s grid counterpart, and reached through it -- the card says
    which family wrote the checkpoint, so a caller does not have to.

    The prediction is scattered back onto the entry's own mesh as ordinary
    ``<column>_pred`` / ``<column>_error`` point data, so the result is a plain
    ``.vtu`` that every viewer already understands and the dashboard's
    prediction preview needs no change to show.
    """
    loaded = _load_checkpoint(checkpoint, device)
    read_kwargs = dict(loaded.card.get("read", {}))

    manifest = DatasetManifest.load(manifest)
    entries = list(manifest.entries(split=split))
    if entry_ids:
        wanted = set(entry_ids)
        entries = [e for e in entries if e.id in wanted]
    os.makedirs(output_dir, exist_ok=True)

    rows = []
    for entry in entries:
        series = entry.time_series(**read_kwargs)
        time_value, mesh = series[step]
        target_mesh = (
            entry.target_time_series(**read_kwargs)[step][1] if entry.target else None
        )
        out_mesh, row = _predict_grid_mesh(
            loaded, mesh, target_mesh, f"entry '{entry.id}'"
        )
        suffix = "" if step == 0 else f"_s{step}"
        out_path = os.path.join(output_dir, f"{entry.id}{suffix}.vtu")
        write(out_path, out_mesh)
        row["entry_id"] = entry.id
        row["time"] = time_value
        row["output_path"] = out_path
        rows.append(row)
    return rows


#: The family -> trainer and block -> per-mesh body / normalizer tables.
#: Keyed lookups, never an `else` branch: an unknown family fails by name
#: instead of silently running some other family's code.
_RUNNERS = {
    "meshgraphnet": _run_graph,
    "srresnet": _run_grid,
    "fno": _run_grid,
    "afno": _run_grid,
    "deeponet": _run_operator,
}
_PREDICT_BODIES = {
    "Graph": _predict_graph_mesh,
    "Grid": _predict_grid_mesh,
    "Operator": _predict_operator_mesh,
}
_NORM_LOADERS = {
    "Graph": _norm_tensors,
    "Grid": _grid_norm_tensors,
    "Operator": _operator_norm_tensors,
}
#: The manifest-walking predictors `predict` hands off to, by block (the
#: graph family's loop is inline in `predict` itself).
_MANIFEST_PREDICTORS = {"Grid": predict_grid, "Operator": predict_operator}


def _spectrum_rel_l2(pred, truth, spec):
    """Relative L2 between the two azimuthally averaged power spectra.

    The number a superresolution result should be reported against: a pointwise
    error cannot tell a field with the right small-scale content from a
    plausible smoothed one, and this can. ``None`` on an anisotropic lattice,
    where a shell average is not meaningful.
    """
    from .._grid_transfer import power_spectrum

    if not spec.is_isotropic:
        return None
    a = power_spectrum(pred, spec).power
    b = power_spectrum(truth, spec).power
    denominator = float(np.sqrt(np.sum(np.square(b))))
    if denominator == 0.0:
        return None
    return float(np.sqrt(np.sum(np.square(a - b))) / denominator)


def _attach(mesh, kind, name, values, tess=None):
    """Write a per-node or per-cell array onto ``mesh``.

    ``tess`` is the :class:`~meshioplusplus.Tessellation` a prediction was
    actually computed through (or ``None``, the ordinary case): when given,
    ``values`` are on the TESSELLATED sample's own points/cells, so they are
    first written back onto ``mesh``'s own (pre-tessellation) points via
    :meth:`~meshioplusplus.Tessellation.scatter` (node) or
    :meth:`~meshioplusplus.Tessellation.aggregate` (cell) -- the map from
    each simplex to the cell it was carved out of. With ``tess=None`` this
    function's behaviour is byte-identical to before tessellation existed.
    """
    if tess is not None:
        values = tess.scatter(values) if kind == "node" else tess.aggregate(values)[1]
    if kind == "node":
        mesh.point_data[name] = values
        return
    bases = block_bases(mesh.cells)
    mesh.cell_data[name] = [
        values[bases[i] : bases[i + 1]] for i in range(len(mesh.cells))
    ]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="python -m meshioplusplus.physicsnemo.train",
        description=__doc__.split("\n\n")[0],
    )
    ap.add_argument(
        "--spec", required=True, help="the PascalCase training spec (JSON file)"
    )
    args = ap.parse_args(argv)
    try:
        progress = run(args.spec)
    except ImportError as e:
        print(str(e), file=sys.stderr)
        return 1
    return 143 if progress.get("stopped") else 0


if __name__ == "__main__":
    sys.exit(main())
