"""The PyTorch Geometric ``Dataset`` over a dataset manifest.

The ONLY module in the package importing ``torch_geometric`` (and, through
it, ``torch``) -- reached exclusively through the lazy ``make_dataset``
factory. This is the training path: MeshGraphNet consumes PyG ``Data``
objects, and PyG's ``DataLoader`` batches variable-size graphs natively.
"""

import torch
from torch_geometric.data import Data

from . import _flat_items, _read_sample, _split_kwargs


class MeshManifestDataset(torch.utils.data.Dataset):
    """``torch.utils.data.Dataset`` of ``torch_geometric.data.Data``
    (``pos``/``x``/``y``/``edge_index``/``edge_attr``) over the flat
    (entry, step) index of a :class:`~meshioplusplus.DatasetManifest`.

    One mesh is read per ``__getitem__`` -- two when ``target_offset >= 1``
    pairs step k with step k+offset -- and nothing is cached across accesses
    (the sequence streaming invariant). ``.schema`` reads sample 0 once and
    reports the recorded feature contract for the whole dataset -- it is
    deliberately not attached to each ``Data`` (a non-tensor attribute would
    muddy PyG's collate).
    """

    def __init__(self, manifest, *, split=None, device=None, **kwargs):
        graph_kwargs, read_kwargs = _split_kwargs(dict(kwargs))
        self._graph_kwargs = graph_kwargs
        offset = int(graph_kwargs.get("target_offset", 0))
        self._items = _flat_items(manifest, split, read_kwargs, offset)
        self._device = device
        self._schema = None

    def __len__(self):
        return len(self._items)

    def __getitem__(self, index):
        entry_id, series, step = self._items[index]
        _, sample = _read_sample(series, step, self._graph_kwargs)
        if self._schema is None:
            self._schema = sample.schema
        tensors = {
            name: torch.from_numpy(array) for name, array in sample.arrays.items()
        }
        data = Data(**tensors)
        return data.to(self._device) if self._device is not None else data

    @property
    def schema(self):
        """The recorded feature contract (from the first sample; one read)."""
        if self._schema is None and self._items:
            self[0]
        return self._schema
