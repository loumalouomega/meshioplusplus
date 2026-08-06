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

from dataclasses import dataclass

import numpy as np

from ..__about__ import __version__
from .._dataset import DatasetManifest
from .._gpu import _require_framework
from .._interop import _importable
from .._ml import FEATURE_SCHEMA_VERSION, edge_index, feature_matrix
from .._regions import block_bases

__all__ = [
    "GraphSample",
    "graph_sample",
    "iter_samples",
    "field_stats",
    "edge_stats",
    "make_reader",
    "make_dataset",
    "has_physicsnemo",
    "has_torch_geometric",
]

GRAPH_SAMPLE_VERSION = 1

_DOC = "doc/physicsnemo.md"

#: keyword arguments consumed by :func:`graph_sample`; everything else a
#: manifest-walking entry point receives is forwarded to ``read()``.
_GRAPH_KWARGS = (
    "fields",
    "target_fields",
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
    - ``y`` -- a second ``feature_matrix`` over ``target_fields`` (same-step
      targets, the steady-state surrogate shape); absent when ``None``.
    - ``edge_index`` -- :func:`~meshioplusplus.edge_index` verbatim
      (``kind``/``undirected``); ``kind="cell"`` switches ``pos``/``x``/``y``
      to cell centroids and cell-located data, so the sample stays coherent.
    - ``edge_attr`` -- ``cat((pos[row] - pos[col], ||.||))``, the stable
      PhysicsNeMo convention; ``edge_features=False`` skips it.

    ``float32=True`` (default -- the training convention) casts the float
    arrays; ``edge_index`` is int64 always. Returns a :class:`GraphSample`.
    """
    if kind not in ("node", "cell"):
        raise ValueError(
            f"meshio++: graph_sample: kind must be 'node' or 'cell', not {kind!r}"
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
        fmy = feature_matrix(
            mesh, location, fields=target_fields, coords=False, regions=False
        )
        arrays["y"] = fmy.matrix.astype(float_dtype, copy=False)
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


def _flat_items(manifest, split, read_kwargs):
    """The flat (entry, step) index: one ``(entry_id, TimeSeries, step)`` per
    sample, resolved once from the manifest's plans -- no mesh read."""
    items = []
    for entry in _as_manifest(manifest).entries(split=split):
        series = entry.time_series(**read_kwargs)
        for step in range(len(series)):
            items.append((entry.id, series, step))
    return items


def iter_samples(manifest, *, split=None, **kwargs):
    """Yield ``(entry_id, time, GraphSample)`` over a manifest's entries.

    A generator honouring the sequence streaming invariant: exactly one mesh
    is alive per yielded sample. ``kwargs`` split into :func:`graph_sample`
    parameters and ``read()`` kwargs (``arrays=...`` narrowing, ...);
    ``manifest`` is a :class:`~meshioplusplus.DatasetManifest` or anything
    its ``load`` accepts.
    """
    graph_kwargs, read_kwargs = _split_kwargs(dict(kwargs))
    for entry_id, series, step in _flat_items(manifest, split, read_kwargs):
        time, mesh = series[step]
        yield entry_id, time, graph_sample(mesh, **graph_kwargs)


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


def field_stats(manifest, *, split=None, location="point", fields=None, **read_kwargs):
    """Per-field normalization stats over a manifest, streaming.

    Returns ``{"<field>_mean": [...], "<field>_std": [...]}`` with one value
    per component -- the Gen-1 ``node_stats.json`` key convention, directly
    usable by PhysicsNeMo's datapipes and by its Gen-2 ``Normalize``
    transform. Plain JSON away from disk: ``json.dump(stats, fh)``. One mesh
    alive at a time; the multi-component grouping comes from the
    ``feature_matrix`` schema, so it is the same rule as the columns.
    """
    accumulators = {}
    order = []
    for entry in _as_manifest(manifest).entries(split=split):
        for _, mesh in entry.time_series(**read_kwargs):
            fm = feature_matrix(
                mesh, location, fields=fields, coords=False, regions=False
            )
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
                acc.add(fm.matrix[:, idx])
    stats = {}
    for name in order:
        acc = accumulators[name]
        stats[f"{name}_mean"] = acc.mean().tolist()
        stats[f"{name}_std"] = acc.std().tolist()
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
