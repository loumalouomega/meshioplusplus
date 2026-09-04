"""PhysicsNeMo adapter: meshes -> the tensors its datapipes expect.

The ``mcp/`` split, applied again: this ``__init__`` is **pure** -- numpy,
stdlib and meshioplusplus only -- and holds all the behaviour (`graph_sample`,
the stats accumulators, the manifest iteration), so it is testable in the
default CI matrix with nothing optional installed. ``_reader.py`` is the only
module importing ``physicsnemo`` (the Gen-2 ``Reader`` subclass), ``_pyg.py``
the only one importing ``torch_geometric`` (the PyG ``Dataset``), and both
are reached through the lazy ``make_reader``/``make_dataset`` factories.

Deliberately **not** imported by ``import meshioplusplus`` (the ``mcp``
subpackage precedent): the gated halves have their own Python floor
(PhysicsNeMo requires >= 3.11 against the package's 3.8) and heavyweight
dependencies, and a class subclassing an optional-dep base cannot live in an
eagerly imported module. Use ``import meshioplusplus.physicsnemo as mpn``.

There is deliberately no ``meshioplusplus[physicsnemo]`` pip extra:
``nvidia-physicsnemo`` hard-depends on ``torch``, whose default Linux wheel
bundles CUDA at multiple GB -- the exact dependency the repo already refuses
to pin for ``to_torch``/``to_jax`` (the CuPy packaging precedent). See
``doc/physicsnemo.md``.
"""

import warnings
from dataclasses import dataclass

import numpy as np

from ..__about__ import __version__
from .._dataset import DatasetManifest
from .._gpu import _require_framework
from .._interop import _emit, _importable
from .._ml import FEATURE_SCHEMA_VERSION, edge_index, feature_matrix
from .._regions import block_bases
from ._train import TrainSpec, default_spec, load_spec

__all__ = [
    "GraphSample",
    "graph_sample",
    "iter_samples",
    "field_stats",
    "edge_stats",
    "make_reader",
    "make_dataset",
    "to_physicsnemo",
    "from_physicsnemo",
    "has_physicsnemo",
    "has_torch_geometric",
    "TrainSpec",
    "default_spec",
    "load_spec",
    "run_training",
    "predict",
]

# Version 2 (v9.30.0): the schema gained `target_offset`/`target_delta` --
# the documented bump rule (a stored v1 schema now compares unequal, which is
# the drift guard doing its job).
GRAPH_SAMPLE_VERSION = 2

_DOC = "doc/physicsnemo.md"

#: keyword arguments consumed by :func:`graph_sample`; everything else a
#: manifest-walking entry point receives is forwarded to ``read()``.
#: ``target_offset`` is also read (without popping) by the iteration layer,
#: which turns it into the per-sample ``target_mesh`` -- ``target_mesh``
#: itself is never an iteration-level kwarg, being derived per step.
_GRAPH_KWARGS = (
    "fields",
    "target_fields",
    "target_delta",
    "target_offset",
    "regions",
    "kind",
    "undirected",
    "edge_features",
    "float32",
)


def has_physicsnemo():
    """Whether PhysicsNeMo is importable. There is deliberately no pip extra
    (install ``nvidia-physicsnemo`` directly; it imports as ``physicsnemo``)."""
    return _importable("physicsnemo")


def has_torch_geometric():
    """Whether PyTorch Geometric is importable. There is deliberately no pip
    extra (install ``torch_geometric`` directly)."""
    return _importable("torch_geometric")


@dataclass(frozen=True)
class GraphSample:
    """One graph's worth of arrays, plus the recorded feature contract.

    ``arrays`` holds plain numpy: ``pos``, ``x``, ``edge_index``, and
    optionally ``y`` and ``edge_attr``. ``x_columns``/``y_columns`` are the
    :func:`~meshioplusplus.feature_matrix` column contract carried through;
    ``schema`` is the JSON-serializable record to store at training time and
    compare at inference time (the feature-drift guard).
    """

    arrays: dict
    x_columns: tuple
    y_columns: tuple
    schema: dict


def _positions(mesh, kind):
    """Graph-vertex positions: mesh points for the node graph, block-major
    cell centroids (the `partition` numbering) for the cell dual."""
    points = np.ascontiguousarray(np.asarray(mesh.points, dtype=np.float64))
    if kind == "node":
        return points
    from .._partition import _centroids

    total = block_bases(mesh.cells)[-1] if mesh.cells else 0
    dim = min(points.shape[1] if points.ndim == 2 else 0, 3)
    return np.ascontiguousarray(_centroids(mesh, total)[:, :dim])


def _edge_attributes(pos, ei):
    """The PhysicsNeMo MeshGraphNet edge-feature convention, exactly:
    relative displacement (source minus destination) plus its norm."""
    row, col = ei
    disp = pos[row] - pos[col]
    norm = np.linalg.norm(disp, axis=-1, keepdims=True)
    return np.concatenate((disp, norm), axis=-1)


def graph_sample(
    mesh,
    *,
    fields=None,
    target_fields=None,
    target_mesh=None,
    target_delta=False,
    target_offset=0,
    regions=True,
    kind="node",
    undirected=True,
    edge_features=True,
    float32=True,
):
    """One mesh as one GNN training sample.

    - ``pos`` -- vertex positions ``(N, dim)``; never folded into ``x`` (the
      MeshGraphNet convention: coordinates enter through edge features only).
    - ``x`` -- :func:`~meshioplusplus.feature_matrix` at the graph's location
      with ``coords=False``: ``fields`` in stated order, then the
      ``region:<name>`` one-hots (the node-type slot); ``regions`` as there.
    - ``y`` -- a second ``feature_matrix`` over ``target_fields``; absent
      when ``None``. Three target modes: same-step (the default, the
      steady-state surrogate shape), **next-step** (``target_mesh`` -- the
      autoregressive t->t+1 pairing; the target's row count must match the
      input's, so a remeshed series fails by name rather than mis-pairing),
      and **delta** (``target_delta=True`` with a ``target_mesh``: y is the
      target's values minus this mesh's own -- the MeshGraphNet increment
      convention; normalize with ``field_stats(delta=...)``).
    - ``edge_index`` -- :func:`~meshioplusplus.edge_index` verbatim
      (``kind``/``undirected``); ``kind="cell"`` switches ``pos``/``x``/``y``
      to cell centroids and cell-located data, so the sample stays coherent.
    - ``edge_attr`` -- ``cat((pos[row] - pos[col], ||.||))``, the stable
      PhysicsNeMo convention; ``edge_features=False`` skips it.

    ``target_offset`` is pure metadata recorded in the schema -- this
    function sees two meshes, not a series, so the iteration layer
    (:func:`iter_samples`/`make_dataset`/`make_reader`) supplies both the
    ``target_mesh`` and the offset it used. ``float32=True`` (default -- the
    training convention) casts the float arrays; ``edge_index`` is int64
    always. Returns a :class:`GraphSample`.
    """
    if kind not in ("node", "cell"):
        raise ValueError(
            f"meshio++: graph_sample: kind must be 'node' or 'cell', not {kind!r}"
        )
    if target_mesh is not None and target_fields is None:
        raise ValueError(
            "meshio++: graph_sample: target_mesh requires target_fields "
            "(there is nothing to pair without target arrays)"
        )
    if target_delta and target_mesh is None:
        raise ValueError(
            "meshio++: graph_sample: target_delta requires target_mesh -- "
            "a delta against the same step is identically zero"
        )
    if target_offset < 0:
        raise ValueError(
            f"meshio++: graph_sample: target_offset must be >= 0, not {target_offset}"
        )
    location = "point" if kind == "node" else "cell"
    float_dtype = np.float32 if float32 else np.float64

    pos = _positions(mesh, kind)
    fmx = feature_matrix(mesh, location, fields=fields, coords=False, regions=regions)
    ei = edge_index(mesh, kind=kind, undirected=undirected)

    arrays = {
        "pos": pos.astype(float_dtype, copy=False),
        "x": fmx.matrix.astype(float_dtype, copy=False),
        "edge_index": ei,
    }
    y_columns = ()
    y_sources = []
    if target_fields is not None:
        y_mesh = mesh if target_mesh is None else target_mesh
        fmy = feature_matrix(
            y_mesh, location, fields=target_fields, coords=False, regions=False
        )
        if fmy.matrix.shape[0] != pos.shape[0]:
            raise ValueError(
                f"meshio++: graph_sample: target_mesh yields "
                f"{fmy.matrix.shape[0]} {location} rows but mesh yields "
                f"{pos.shape[0]} -- t->t+1 pairing requires matching meshes"
            )
        y = fmy.matrix
        if target_delta:
            fmb = feature_matrix(
                mesh, location, fields=target_fields, coords=False, regions=False
            )
            if fmy.columns != fmb.columns:
                raise ValueError(
                    "meshio++: graph_sample: target_delta needs the same "
                    f"columns at both steps, got {fmy.columns} vs {fmb.columns}"
                )
            y = fmy.matrix - fmb.matrix
        arrays["y"] = y.astype(float_dtype, copy=False)
        y_columns = fmy.columns
        y_sources = fmy.schema["sources"]
    if edge_features:
        arrays["edge_attr"] = _edge_attributes(pos, ei).astype(float_dtype, copy=False)

    schema = {
        "graph_sample_version": GRAPH_SAMPLE_VERSION,
        "feature_schema_version": FEATURE_SCHEMA_VERSION,
        "meshioplusplus_version": str(__version__),
        "kind": kind,
        "undirected": bool(undirected),
        "float32": bool(float32),
        "target_offset": int(target_offset),
        "target_delta": bool(target_delta),
        "x_columns": list(fmx.columns),
        "y_columns": list(y_columns),
        "edge_features": (
            ["d" + "xyz"[i] for i in range(pos.shape[1])] + ["norm"]
            if edge_features
            else []
        ),
        "x_sources": fmx.schema["sources"],
        "y_sources": y_sources,
    }
    return GraphSample(
        arrays=arrays,
        x_columns=fmx.columns,
        y_columns=tuple(y_columns),
        schema=schema,
    )


# --------------------------------------------------------------------------- #
# Walking a dataset manifest                                                  #
# --------------------------------------------------------------------------- #
def _as_manifest(manifest):
    if isinstance(manifest, DatasetManifest):
        return manifest
    return DatasetManifest.load(manifest)


def _split_kwargs(kwargs):
    graph_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k in _GRAPH_KWARGS}
    return graph_kwargs, kwargs


def _flat_items(manifest, split, read_kwargs, offset=0):
    """The flat (entry, step) index: one ``(entry_id, TimeSeries, step)`` per
    sample, resolved once from the manifest's plans -- no mesh read.

    With ``offset >= 1`` (t->t+offset pairing) each entry contributes
    ``len(series) - offset`` samples; an entry too short to pair contributes
    nothing, with one warning naming it -- emitted here at index-build time,
    never once per epoch.
    """
    items = []
    for entry in _as_manifest(manifest).entries(split=split):
        series = entry.time_series(**read_kwargs)
        if offset and len(series) <= offset:
            warnings.warn(
                f"meshio++: physicsnemo: entry '{entry.id}' has "
                f"{len(series)} step(s), fewer than target_offset+1="
                f"{offset + 1} -- it contributes no samples",
                stacklevel=2,
            )
            continue
        for step in range(max(len(series) - offset, 0)):
            items.append((entry.id, series, step))
    return items


def _read_sample(series, step, graph_kwargs):
    """One (or two, when pairing) mesh reads -> ``(time, GraphSample)``.

    The single owner of the read-and-sample step, shared by
    :func:`iter_samples`, the Gen-2 ``Reader`` and the PyG ``Dataset`` so
    the three cannot drift. ``TimeSeries`` caches nothing, so a paired
    sample (``target_offset >= 1``) costs exactly two reads -- **at most two
    meshes alive**, the pairing amendment to the streaming invariant.
    """
    offset = int(graph_kwargs.get("target_offset", 0))
    time, mesh = series[step]
    target = series[step + offset][1] if offset else None
    return time, graph_sample(mesh, target_mesh=target, **graph_kwargs)


def iter_samples(manifest, *, split=None, **kwargs):
    """Yield ``(entry_id, time, GraphSample)`` over a manifest's entries.

    A generator honouring the sequence streaming invariant: one mesh is
    alive per yielded sample -- at most two when ``target_offset >= 1``
    pairs step k with step k+offset (the yielded ``time`` is the *input*
    step's). ``kwargs`` split into :func:`graph_sample` parameters and
    ``read()`` kwargs (``arrays=...`` narrowing, ...); ``manifest`` is a
    :class:`~meshioplusplus.DatasetManifest` or anything its ``load``
    accepts.
    """
    graph_kwargs, read_kwargs = _split_kwargs(dict(kwargs))
    offset = int(graph_kwargs.get("target_offset", 0))
    for entry_id, series, step in _flat_items(manifest, split, read_kwargs, offset):
        time, sample = _read_sample(series, step, graph_kwargs)
        yield entry_id, time, sample


class _Moments:
    """Streaming per-component mean/std, the PhysicsNeMo convention:
    ``std = sqrt(mean(x^2) - mean(x)^2)`` accumulated in float64."""

    def __init__(self, width):
        self.count = 0
        self.total = np.zeros(width, dtype=np.float64)
        self.total_sq = np.zeros(width, dtype=np.float64)

    def add(self, rows):
        rows = np.asarray(rows, dtype=np.float64)
        self.count += rows.shape[0]
        self.total += rows.sum(axis=0)
        self.total_sq += np.square(rows).sum(axis=0)

    def mean(self):
        return self.total / max(self.count, 1)

    def std(self):
        variance = self.total_sq / max(self.count, 1) - np.square(self.mean())
        return np.sqrt(np.maximum(variance, 0.0))


def field_stats(
    manifest, *, split=None, location="point", fields=None, delta=False, **read_kwargs
):
    """Per-field normalization stats over a manifest, streaming.

    Returns ``{"<field>_mean": [...], "<field>_std": [...]}`` with one value
    per component -- the Gen-1 ``node_stats.json`` key convention, directly
    usable by PhysicsNeMo's datapipes and by its Gen-2 ``Normalize``
    transform. Plain JSON away from disk: ``json.dump(stats, fh)``. One mesh
    alive at a time; the multi-component grouping comes from the
    ``feature_matrix`` schema, so it is the same rule as the columns.

    ``delta`` accumulates the statistics of **consecutive-step differences**
    instead -- what normalizing ``graph_sample(target_delta=True)`` targets
    needs -- under ``{field}_diff_mean``/``{field}_diff_std`` keys (still
    plain JSON, no new stats file format). ``True`` means lag 1; an int
    ``n >= 1`` is the lag matching ``target_offset=n``. Only the previous
    feature *matrices* (numpy) are held between steps, so one mesh stays
    alive; entries with too few steps contribute nothing.
    """
    # Identity checks, not `in`: 0 == False in Python, and a `delta=0` almost
    # certainly meant something -- name it rather than silently doing nothing.
    if delta is False or delta is None:
        lag = 0
    else:
        lag = 1 if delta is True else int(delta)
        if lag < 1:
            raise ValueError(
                f"meshio++: field_stats: delta must be True or an int >= 1, "
                f"not {delta!r}"
            )
    suffix = "_diff" if lag else ""
    accumulators = {}
    order = []
    for entry in _as_manifest(manifest).entries(split=split):
        previous = []  # the last `lag` feature matrices (numpy, never meshes)
        for _, mesh in entry.time_series(**read_kwargs):
            fm = feature_matrix(
                mesh, location, fields=fields, coords=False, regions=False
            )
            if lag:
                previous.append(fm)
                if len(previous) <= lag:
                    continue
                oldest = previous.pop(0)
                if oldest.columns != fm.columns:
                    raise ValueError(
                        f"meshio++: field_stats: entry '{entry.id}' changes "
                        f"columns across steps ({oldest.columns} vs "
                        f"{fm.columns}) -- delta stats need a stable field set"
                    )
                if oldest.matrix.shape[0] != fm.matrix.shape[0]:
                    raise ValueError(
                        f"meshio++: field_stats: entry '{entry.id}' changes "
                        f"row count across steps ({oldest.matrix.shape[0]} vs "
                        f"{fm.matrix.shape[0]}) -- delta stats need a stable mesh"
                    )
                matrix = fm.matrix - oldest.matrix
            else:
                matrix = fm.matrix
            groups = {}
            for column, source in zip(fm.columns, fm.schema["sources"]):
                groups.setdefault(source["source"], []).append(column)
            for name, columns in groups.items():
                idx = [fm.columns.index(c) for c in columns]
                acc = accumulators.get(name)
                if acc is None:
                    acc = accumulators[name] = _Moments(len(idx))
                    order.append(name)
                elif acc.total.shape[0] != len(idx):
                    raise ValueError(
                        f"meshio++: field_stats: '{name}' changes component "
                        f"count across the dataset ({acc.total.shape[0]} vs "
                        f"{len(idx)})"
                    )
                acc.add(matrix[:, idx])
    stats = {}
    for name in order:
        acc = accumulators[name]
        stats[f"{name}{suffix}_mean"] = acc.mean().tolist()
        stats[f"{name}{suffix}_std"] = acc.std().tolist()
    return stats


def edge_stats(manifest, *, split=None, kind="node", undirected=True, **read_kwargs):
    """Edge-feature normalization stats over a manifest, streaming.

    Returns ``{"edge_mean": [...], "edge_std": [...]}`` per component of the
    ``(disp, norm)`` edge attribute -- the Gen-1 ``edge_stats.json`` key
    convention.
    """
    acc = None
    for entry in _as_manifest(manifest).entries(split=split):
        for _, mesh in entry.time_series(**read_kwargs):
            pos = _positions(mesh, kind)
            ei = edge_index(mesh, kind=kind, undirected=undirected)
            attr = _edge_attributes(pos, ei)
            if acc is None:
                acc = _Moments(attr.shape[1])
            elif acc.total.shape[0] != attr.shape[1]:
                raise ValueError(
                    "meshio++: edge_stats: the edge-attribute width changes "
                    f"across the dataset ({acc.total.shape[0]} vs {attr.shape[1]})"
                )
            acc.add(attr)
    if acc is None:
        return {"edge_mean": [], "edge_std": []}
    return {"edge_mean": acc.mean().tolist(), "edge_std": acc.std().tolist()}


# --------------------------------------------------------------------------- #
# Framework objects (lazy factories over the gated modules)                   #
# --------------------------------------------------------------------------- #
def make_reader(manifest, *, split=None, **kwargs):
    """A PhysicsNeMo Gen-2 ``Reader`` over a manifest: ``__len__`` +
    ``_load_sample(i) -> dict[str, torch.Tensor]`` (the sanctioned datapipe
    extension point). Needs ``nvidia-physicsnemo``; raises a named install
    error otherwise. ``kwargs`` as in :func:`iter_samples`."""
    _require_framework(
        "make_reader", "physicsnemo", "pip install nvidia-physicsnemo", doc=_DOC
    )
    from ._reader import MeshManifestReader

    return MeshManifestReader(manifest, split=split, **kwargs)


# --------------------------------------------------------------------------- #
# The physicsnemo.mesh.Mesh bridge (v9.30.0)                                  #
# --------------------------------------------------------------------------- #
#: physicsnemo.mesh.Mesh holds exactly ONE simplex kind, named by the cells
#: tensor's trailing dim (n_manifold_dims = k - 1).
_SIMPLEX_BY_DIM = {0: "vertex", 1: "line", 2: "triangle", 3: "tetra"}
_TYPE_BY_WIDTH = {1: "vertex", 2: "line", 3: "triangle", 4: "tetra"}


def _to_physicsnemo_payload(mesh, *, manifold_dim="auto", float32=True):
    """The pure half of :func:`to_physicsnemo`: plain numpy, no torch.

    ``physicsnemo.mesh.Mesh`` is single-manifold-dim and simplicial-only (its
    own docs: subdivide non-simplicial elements before use), so this selects
    ONE topological dimension, tessellates it with the existing
    ``convert_cells`` operations (which replicate ``cell_data`` rows to
    children natively), and drops everything the target cannot hold -- each
    drop recorded in ``notes``, never silent. Every returned buffer is
    freshly owned: upstream in-place edits invalidate that Mesh's caches, so
    sharing memory with the source mesh would be a spooky-action bug.
    """
    from .._convert_cells import convert_cells

    float_dtype = np.float32 if float32 else np.float64
    notes = []

    dims = sorted({blk.dim for blk in mesh.cells})
    if manifold_dim == "auto":
        target = dims[-1] if dims else 0
    else:
        target = int(manifold_dim)
        if target not in _SIMPLEX_BY_DIM:
            raise ValueError(
                f"meshio++: to_physicsnemo: manifold_dim must be 'auto' or "
                f"0..3, not {manifold_dim!r}"
            )
        if dims and target not in dims:
            raise ValueError(
                f"meshio++: to_physicsnemo: the mesh has no dimension-{target} "
                f"cells (dimensions present: {dims or 'none'})"
            )
    spatial = mesh.points.shape[1] if mesh.points.ndim == 2 else 0
    if target > spatial:
        raise ValueError(
            f"meshio++: to_physicsnemo: manifold dimension {target} exceeds "
            f"the spatial dimension {spatial} -- a tetra Mesh needs 3-D "
            f"points (physicsnemo.mesh.Mesh's own constraint)"
        )

    simplex = _SIMPLEX_BY_DIM[target]
    converted = mesh
    if dims and any(
        blk.dim == target and blk.type != simplex for blk in converted.cells
    ):
        # Higher-order first, then decompose -- both idempotent on cells they
        # do not apply to, and both replicate cell_data rows to children.
        converted = convert_cells(converted, mode="linearize")
        converted = convert_cells(converted, mode="simplexify")
        notes.append(
            f"non-{simplex} dimension-{target} cells were tessellated "
            f"(linearize + simplexify); cell_data rows replicated to children"
        )

    kept = [
        i
        for i, blk in enumerate(converted.cells)
        if blk.dim == target and blk.type == simplex
    ]
    dropped = [
        (blk.type, len(blk.data))
        for i, blk in enumerate(converted.cells)
        if i not in kept
    ]
    if dropped:
        listing = ", ".join(f"{t} x {n}" for t, n in dropped)
        notes.append(
            f"physicsnemo.mesh.Mesh is single-manifold-dim; dropped {listing} "
            f"-- export each dimension as its own Mesh"
        )

    if kept:
        cells = np.ascontiguousarray(
            np.concatenate([np.asarray(converted.cells[i].data) for i in kept]),
            dtype=np.int64,
        )
    else:
        cells = None
        if dims:
            notes.append("no cells at the selected dimension survive; points only")

    cell_data = {}
    for name, arrays in converted.cell_data.items():
        rows = [np.asarray(arrays[i]) for i in kept]
        if not rows:
            continue
        if not all(r.dtype.kind in "fiub" for r in rows):
            notes.append(f"cell_data '{name}' is non-numeric and was dropped")
            continue
        cell_data[name] = np.ascontiguousarray(np.concatenate(rows))

    point_data = {}
    for name, array in converted.point_data.items():
        array = np.asarray(array)
        if array.dtype.kind not in "fiub":
            notes.append(f"point_data '{name}' is non-numeric and was dropped")
            continue
        point_data[name] = np.array(array, copy=True)

    global_data = {}
    for name, value in converted.field_data.items():
        try:
            array = np.asarray(value)
        except Exception:
            array = np.asarray(None)
        if array.dtype.kind not in "fiub":
            notes.append(f"field_data '{name}' is non-numeric and was dropped")
            continue
        global_data[name] = np.array(array, copy=True)

    if getattr(converted, "regions", None):
        names = sorted({r.name for r in converted.regions})
        notes.append(
            f"regions {names} have no physicsnemo.mesh slot and were dropped "
            f"-- encode membership as features via graph_sample's one-hots"
        )

    return {
        "points": np.array(converted.points, dtype=float_dtype, copy=True),
        "cells": cells,
        "point_data": point_data,
        "cell_data": cell_data,
        "global_data": global_data,
        "manifold_dim": target,
        "notes": notes,
    }


def _from_physicsnemo_payload(payload):
    """The pure half of :func:`from_physicsnemo`: numpy payload -> Mesh."""
    from .._mesh import Mesh

    cells = payload.get("cells")
    blocks = []
    if cells is not None and len(cells):
        cells = np.ascontiguousarray(np.asarray(cells), dtype=np.int64)
        width = cells.shape[1] if cells.ndim == 2 else 0
        if width not in _TYPE_BY_WIDTH:
            raise ValueError(
                f"meshio++: from_physicsnemo: cells with {width} nodes per "
                f"row map to no simplex (vertex/line/triangle/tetra)"
            )
        blocks = [(_TYPE_BY_WIDTH[width], cells)]
    cell_data = {
        name: [np.asarray(array)]
        for name, array in payload.get("cell_data", {}).items()
    }
    return Mesh(
        np.asarray(payload["points"]),
        blocks,
        point_data={
            name: np.asarray(a) for name, a in payload.get("point_data", {}).items()
        },
        cell_data=cell_data if blocks else {},
        field_data={
            name: np.asarray(a) for name, a in payload.get("global_data", {}).items()
        },
    )


def to_physicsnemo(mesh, *, manifold_dim="auto", float32=True):
    """A meshio++ mesh as a ``physicsnemo.mesh.Mesh`` (CPU tensors).

    The target holds exactly one simplex kind: ``manifold_dim`` selects it
    (``"auto"`` = the highest dimension present), non-simplex cells at that
    dimension are tessellated through the existing ``convert_cells``
    operations, and everything the target cannot hold (other-dim blocks,
    regions, non-numeric data) is dropped **with a warning**, never
    silently. Points follow the upstream float32 convention
    (``float32=False`` keeps float64); cells are int64; move to a device
    with the upstream ``Mesh.to(device)``.

    Needs ``nvidia-physicsnemo`` (no pip extra, deliberately); pin
    ``>=2.1,<2.2`` in training projects -- ``physicsnemo.mesh`` is the
    framework's newest surface. See ``doc/physicsnemo.md``.
    """
    _require_framework(
        "to_physicsnemo", "physicsnemo", "pip install nvidia-physicsnemo", doc=_DOC
    )
    payload = _to_physicsnemo_payload(mesh, manifold_dim=manifold_dim, float32=float32)
    _emit("to_physicsnemo", payload["notes"])
    from ._mesh import to_mesh

    return to_mesh(payload)


def from_physicsnemo(pm):
    """A ``physicsnemo.mesh.Mesh`` back as a meshio++ mesh.

    The cells tensor's trailing dim names the one simplex block
    (vertex/line/triangle/tetra); ``point_data``/``cell_data`` leaves come
    back as numpy arrays and ``global_data`` lands in ``field_data``.
    Needs ``nvidia-physicsnemo``; see :func:`to_physicsnemo`.
    """
    _require_framework(
        "from_physicsnemo", "physicsnemo", "pip install nvidia-physicsnemo", doc=_DOC
    )
    from ._mesh import from_mesh

    return _from_physicsnemo_payload(from_mesh(pm))


def make_dataset(manifest, *, split=None, device=None, **kwargs):
    """A ``torch.utils.data.Dataset`` of ``torch_geometric.data.Data``
    (``pos``/``x``/``y``/``edge_index``/``edge_attr``) over a manifest --
    what MeshGraphNet training consumes, batched natively by PyG's
    ``DataLoader``. Needs ``torch_geometric``; raises a named install error
    otherwise. ``device=`` moves each sample's tensors on access."""
    _require_framework(
        "make_dataset", "torch_geometric", "pip install torch_geometric", doc=_DOC
    )
    from ._pyg import MeshManifestDataset

    return MeshManifestDataset(manifest, split=split, device=device, **kwargs)


# --------------------------------------------------------------------------- #
# Training and prediction (v10.24.0; the gated module is ``train.py``)        #
# --------------------------------------------------------------------------- #
def run_training(spec, *, log=print):
    """Train MeshGraphNet per a :class:`TrainSpec` (or its dict / JSON text /
    file path) into the spec's run directory; returns the final progress
    record. The same thing ``python -m meshioplusplus.physicsnemo.train
    --spec`` and the dashboard's ``train_start`` run. Needs ``torch_geometric``
    and ``nvidia-physicsnemo``; raises a named install error otherwise.
    Named ``run_training`` rather than ``train`` so importing the ``train``
    submodule can never shadow it. See ``doc/physicsnemo.md``."""
    _require_framework(
        "run_training", "torch_geometric", "pip install torch_geometric", doc=_DOC
    )
    from .train import run

    return run(spec, log=log)


def predict(
    checkpoint,
    manifest,
    *,
    entry_ids=None,
    split="test",
    step=0,
    output_dir,
    device="auto",
):
    """Predict with a trained ``.mdlus`` checkpoint (and its model card) over
    a manifest's entries, writing ``<column>_pred``/``<column>_error`` back
    as data arrays into ``output_dir/<entry_id>.vtu``; returns one row per
    entry (``entry_id``, ``output_path``, ``rmse``, ``max_error``). Needs
    ``nvidia-physicsnemo`` and ``torch_geometric``."""
    _require_framework(
        "predict", "physicsnemo", "pip install nvidia-physicsnemo", doc=_DOC
    )
    from .train import predict as _predict

    return _predict(
        checkpoint,
        manifest,
        entry_ids=entry_ids,
        split=split,
        step=step,
        output_dir=output_dir,
        device=device,
    )
