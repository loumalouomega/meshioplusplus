"""A plain torch ``Dataset`` over a manifest's grid pairs.

The grid counterpart of :mod:`._pyg`, and deliberately a separate module for one
reason: a convolutional model needs **torch only**. PyTorch Geometric exists to
batch ragged graphs, which a dense ``(C, D, H, W)`` tensor is not, so requiring
it for a superresolution run would be demanding a large dependency the model
never touches -- and, since its prebuilt wheels lag torch releases, one that is
often simply unavailable.

Imports ``torch`` at call time, never at module scope.
"""

from __future__ import annotations

from .._dataset import DatasetManifest


def make_grid_dataset(manifest, *, split=None, **kwargs):
    """A ``torch.utils.data.Dataset`` of ``(x, y)`` grid pairs.

    ``x`` is the coarse grid and ``y`` the fine one, both float tensors shaped
    ``(C, D, H, W)`` -- exactly what ``SRResNet.forward`` consumes and produces,
    so the default collate function batches them with no custom collator.

    The index is built once from the manifest's plans (no mesh read); each
    ``__getitem__`` reads one mesh, or two for a paired entry, so the streaming
    invariant holds through the DataLoader.
    """
    import torch

    from . import _grid_flat_items, _read_grid_sample

    grid_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k != "read_kwargs"}
    read_kwargs = kwargs.pop("read_kwargs", {})
    if isinstance(manifest, str) or not isinstance(manifest, DatasetManifest):
        manifest = DatasetManifest.load(manifest)
    items = _grid_flat_items(manifest, split, read_kwargs)

    class GridPairDataset(torch.utils.data.Dataset):
        #: The recorded contract of sample 0, read once. `run()` sizes the model
        #: from this rather than from `len(Fields)`, which is not the same
        #: number once a multi-component array expands into its columns.
        _schema = None

        def __len__(self):
            return len(items)

        def __getitem__(self, index):
            _, series, target, step = items[index]
            _, sample = _read_grid_sample(series, target, step, grid_kwargs)
            return (
                torch.from_numpy(sample.arrays["x"]).float(),
                torch.from_numpy(sample.arrays["y"]).float(),
            )

        @property
        def schema(self):
            if self._schema is None:
                if not items:
                    raise ValueError(
                        "meshio++: train: the split yields no samples, so there "
                        "is no schema to read"
                    )
                _, series, target, step = items[0]
                _, sample = _read_grid_sample(series, target, step, grid_kwargs)
                type(self)._schema = sample.schema
            return self._schema

    return GridPairDataset()
