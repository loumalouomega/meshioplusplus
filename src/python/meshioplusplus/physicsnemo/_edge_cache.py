"""The opt-in edge-index cache of the sample path (roadmap 5.10.4).

A series whose topology does not change between steps -- the common case:
one mesh, many time steps of fields -- rebuilds the same ``edge_index`` at
every draw. ``edge_cache=N`` keeps the last ``N`` distinct edge sets, keyed by
the topology (the point count and every cell block's type, dtype, shape and a
sample of its connectivity, then the whole connectivity compared) and the
graph options that shape them (``kind``, ``undirected``). Only the edge arrays
and a copy of the connectivity are held, never a mesh, so the streaming
invariant -- at most one mesh alive per sample, two when pairing -- is
untouched.

Edges built by ``proximity`` or ``world_edges`` depend on positions and are
never cached. A mesh with a ragged or object-dtype block has no cheap digest
and is never cached either. Each worker process of a ``DataLoader`` holds its
own cache.
"""

import hashlib
from collections import OrderedDict

import numpy as np

#: Connectivity values sampled per block for the lookup key; the key only
#: finds a candidate, which is then compared whole.
_KEY_SAMPLES = 1024


def _topology(mesh):
    """``(key, blocks)`` for ``mesh``'s topology, or ``None`` when it has no
    cheap one: the key holds the point count, every block's type, dtype and
    shape and a strided sample of its connectivity; ``blocks`` are the
    connectivity arrays themselves, which a candidate must equal."""
    parts = [len(mesh.points)]
    blocks = []
    h = hashlib.blake2b(digest_size=16)
    for block in mesh.cells:
        data = block.data
        if not isinstance(data, np.ndarray) or data.dtype == object:
            return None
        parts.append((str(block.type), data.dtype.str, data.shape))
        flat = data.reshape(-1)
        step = max(1, flat.size // _KEY_SAMPLES)
        h.update(np.ascontiguousarray(flat[::step]).view(np.uint8).data)
        blocks.append(data)
    return (tuple(parts), h.digest()), blocks


class EdgeCache:
    """A bounded LRU of edge-index arrays keyed by topology.

    An entry holds the edge array and a copy of the connectivity it was built
    from (no points, no data), so a hit is a whole-array comparison rather
    than a rebuild; a key collision with other connectivity is a miss."""

    def __init__(self, maxsize):
        maxsize = int(maxsize)
        if maxsize < 1:
            raise ValueError(
                f"meshio++: physicsnemo: edge_cache must be a positive "
                f"number of entries, got {maxsize}"
            )
        self.maxsize = maxsize
        self.hits = 0
        self.misses = 0
        self._entries = OrderedDict()

    def __len__(self):
        return len(self._entries)

    def edges(self, mesh, kind, undirected, build):
        """``build()``'s edges for ``mesh``, from the cache when its topology
        was seen; a fresh, writable copy either way (torch tensors share
        memory with the arrays they wrap, so a caller's edit must not reach
        the cached one)."""
        topology = _topology(mesh)
        if topology is None:
            return build()
        digest, blocks = topology
        key = (kind, bool(undirected), digest)
        cached = self._entries.get(key)
        if cached is not None and all(
            np.array_equal(a, b) for a, b in zip(cached[0], blocks)
        ):
            self._entries.move_to_end(key)
            self.hits += 1
            return cached[1].copy()
        self.misses += 1
        ei = build()
        self._entries[key] = ([b.copy() for b in blocks], ei.copy())
        self._entries.move_to_end(key)
        if len(self._entries) > self.maxsize:
            self._entries.popitem(last=False)
        return ei


def pop_edge_cache(kwargs):
    """Remove ``edge_cache`` from ``kwargs`` and return its cache (or None).

    ``edge_cache`` is ``None``/``0`` (off, the default), a positive number of
    entries, or an :class:`EdgeCache` to share between readers.
    """
    value = kwargs.pop("edge_cache", None)
    if value is None or value is False:
        return None
    if isinstance(value, EdgeCache):
        return value
    if value is True:
        return EdgeCache(1)
    if int(value) == 0:
        return None
    return EdgeCache(value)
