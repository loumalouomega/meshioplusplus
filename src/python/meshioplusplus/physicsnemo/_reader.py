"""The Gen-2 PhysicsNeMo ``Reader`` over a dataset manifest.

The ONLY module in the package importing ``physicsnemo`` (and, through it,
``torch``) -- reached exclusively through the lazy ``make_reader`` factory,
never at package import time. The ``Reader`` ABC is the framework's sanctioned
datapipe extension point: two abstract methods, everything else (TensorDict
wrapping, pin_memory, per-sample RNG) supplied by the base class.
"""

import torch
from physicsnemo.datapipes.readers.base import Reader

from . import _flat_items, _split_kwargs, graph_sample


class MeshManifestReader(Reader):
    """``__len__`` + ``_load_sample(i) -> dict[str, torch.Tensor]`` over the
    flat (entry, step) index of a :class:`~meshioplusplus.DatasetManifest`.

    The index is resolved once from the manifest's plans (no mesh read);
    ``_load_sample`` reads exactly one mesh and returns
    :func:`~meshioplusplus.physicsnemo.graph_sample`'s arrays as CPU tensors.
    Batching variable-size graphs across meshes is the PyG path's job
    (``make_dataset``); TensorDict has no ragged-graph batching convention.
    """

    def __init__(self, manifest, *, split=None, **kwargs):
        super().__init__()
        graph_kwargs, read_kwargs = _split_kwargs(dict(kwargs))
        self._graph_kwargs = graph_kwargs
        self._items = _flat_items(manifest, split, read_kwargs)

    def __len__(self):
        return len(self._items)

    def _load_sample(self, index):
        entry_id, series, step = self._items[index]
        _, mesh = series[step]
        sample = graph_sample(mesh, **self._graph_kwargs)
        return {name: torch.from_numpy(array) for name, array in sample.arrays.items()}


# Hydra addressability (`_target_: ${dp:MeshManifestReader}`) is a bonus, not
# a pinned contract: the registry moved once already in the Gen-1 -> Gen-2
# transition, so a missing/reshaped register() must not break the reader.
try:
    from physicsnemo.datapipes.registry import register

    MeshManifestReader = register()(MeshManifestReader)
except Exception:  # pragma: no cover - registry layout is framework-version-bound
    pass
