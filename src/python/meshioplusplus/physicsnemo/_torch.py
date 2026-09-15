"""Plain torch ``Dataset``s over a manifest: grid pairs and operator samples.

The grid and operator counterparts of :mod:`._pyg`, and deliberately a separate
module for one reason: a convolutional model, a neural operator on a grid and a
DeepONet need **torch only**. PyTorch Geometric exists to batch ragged graphs,
which a dense ``(C, D, H, W)`` tensor or a ``(p,)``/``(T, k)`` pair is not, so
requiring it for these runs would be demanding a large dependency the model
never touches -- and, since its prebuilt wheels lag torch releases, one that is
often simply unavailable.

Imports ``torch`` at call time, never at module scope: the dataset classes are
built inside the factories so this module stays importable with no torch.
"""

from __future__ import annotations

from .._dataset import DatasetManifest


def _as_manifest(manifest):
    if isinstance(manifest, str) or not isinstance(manifest, DatasetManifest):
        return DatasetManifest.load(manifest)
    return manifest


def _no_samples():
    return ValueError(
        "meshio++: train: the split yields no samples, so there is no schema to " "read"
    )


def make_grid_dataset(manifest, *, split=None, **kwargs):
    """A ``torch.utils.data.Dataset`` of ``(x, y)`` grid pairs.

    ``x`` is the coarse grid and ``y`` the fine one, both float tensors shaped
    ``(C, D, H, W)`` -- or ``(C, H, W)`` under a squeeze -- exactly what the
    grid families' ``forward`` consumes and produces, so the default collate
    function batches them with no custom collator.

    The index is built once from the manifest's plans (no mesh read); each
    ``__getitem__`` reads one mesh, or two for a paired entry, so the streaming
    invariant holds through the DataLoader.
    """
    import torch

    from . import _grid_flat_items, _read_grid_sample

    grid_kwargs = {k: kwargs.pop(k) for k in list(kwargs) if k != "read_kwargs"}
    read_kwargs = kwargs.pop("read_kwargs", {})
    items = _grid_flat_items(_as_manifest(manifest), split, read_kwargs)

    class GridPairDataset(torch.utils.data.Dataset):
        def __init__(self):
            # The recorded contract of sample 0, read once PER INSTANCE. It
            # used to be cached on the class, so the second dataset built in a
            # process (the validation split) inherited the first one's schema.
            self._schema = None

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
            """`run()` sizes the model from this rather than from
            `len(Fields)`, which is not the same number once a multi-component
            array expands into its columns."""
            if self._schema is None:
                if not items:
                    raise _no_samples()
                _, series, target, step = items[0]
                _, sample = _read_grid_sample(series, target, step, grid_kwargs)
                self._schema = sample.schema
            return self._schema

    return GridPairDataset()


def make_operator_dataset(manifest, *, split=None, read_kwargs=None, **kwargs):
    """A ``torch.utils.data.Dataset`` of ``(params, y)`` pairs for a DeepONet.

    ``params`` is the case's parameter vector ``(p,)`` and ``y`` the target
    field at the trunk points ``(T, k)``; the **trunk itself is shared** across
    the dataset (fixed geometry, checked by name at index-build time) and is
    exposed once as ``.trunk``, a ``(T, 3)`` float tensor read from sample 0,
    so the model is called as ``model(params_batch, trunk)`` with the trunk
    broadcast over the batch -- the installed ``DeepONet``'s own contract.
    ``kwargs`` are :func:`operator_sample`'s.
    """
    import torch

    from . import _operator_flat_items, operator_sample

    read_kwargs = dict(read_kwargs or {})
    if "parameter_names" not in kwargs:
        raise ValueError("meshio++: train: parameter_names is required")
    items = _operator_flat_items(
        _as_manifest(manifest), split, read_kwargs, kwargs["parameter_names"]
    )

    class OperatorDataset(torch.utils.data.Dataset):
        def __init__(self):
            self._first = None

        def __len__(self):
            return len(items)

        def _sample(self, index):
            _, series, step, metadata = items[index]
            _, mesh = series[step]
            return operator_sample(mesh, metadata, **kwargs)

        def __getitem__(self, index):
            sample = self._sample(index)
            return (
                torch.from_numpy(sample.arrays["params"]).float(),
                torch.from_numpy(sample.arrays["y"]).float(),
            )

        def _first_sample(self):
            if self._first is None:
                if not items:
                    raise _no_samples()
                self._first = self._sample(0)
            return self._first

        @property
        def schema(self):
            return self._first_sample().schema

        @property
        def trunk(self):
            """The shared ``(T, 3)`` trunk, as a float tensor."""
            return torch.from_numpy(self._first_sample().arrays["trunk"]).float()

        @property
        def trunk_indices(self):
            """The budget's selected point indices (``None`` for ``"points"``)."""
            arrays = self._first_sample().arrays
            return None if "trunk_indices" not in arrays else arrays["trunk_indices"]

    return OperatorDataset()
