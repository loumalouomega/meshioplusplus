"""Partition a mesh into N balanced pieces for domain decomposition.

The count-driven complement to the criterion-driven :func:`meshioplusplus.split`.
Two methods: ``"sfc"`` (Hilbert space-filling-curve cut of the cell centroids —
always available, deterministic, byte-identical between the C++ core and this
numpy fallback) and ``"kahip"`` (KaHIP's serial ``kaffpa()`` on the shared-face
dual graph — the quality path, available when the C++ core was built with
``-DMESHIOPLUSPLUS_WITH_KAHIP=ON`` or when the MIT ``kahip`` PyPI package is
installed). ``"auto"`` picks KaHIP when available, else SFC. Requesting
``"kahip"`` when neither exists raises by name — never a silent downgrade.

Public API:

* :func:`partition` — decompose; returns a list of exactly ``nparts`` meshes.
* :func:`partition_labels` — the per-cell part assignment only (one Int64 array
  per cell block, block-aligned like a ``cell_data`` entry).

Pieces keep the input block structure 1:1 (empty blocks included, unlike
``split``), so concatenating them reproduces the input mesh: every cell lands
in exactly one piece when ``ghost_layers == 0``; a positive ``ghost_layers``
grows each piece by that many shared-node BFS layers of other parts' cells,
tagged ``partition:ghost`` (0 = owned, L = reached at layer L), so the pieces
then overlap.
"""

from __future__ import annotations

import numpy as np

from ._mesh import Mesh, topological_dimension

_SFC_BITS = 21


# --------------------------------------------------------------------------- #
# numpy SFC path (must stay label-identical to src/cpp/src/operations/partition.cpp)
# --------------------------------------------------------------------------- #
def _hilbert_keys(q):
    """Skilling AxesToTranspose Hilbert distance of quantized (n, 3) points.

    The vectorized twin of ``detail::sfc_hilbert_key`` — every bitwise step is
    lane-wise exact, so the keys (and therefore the labels) match the C++ core
    bit for bit.
    """
    X = [q[:, i].astype(np.uint32) for i in range(3)]
    bits = _SFC_BITS
    M = 1 << (bits - 1)
    # Inverse undo excess work.
    Q = M
    while Q > 1:
        P = np.uint32(Q - 1)
        for i in range(3):
            mask = (X[i] & np.uint32(Q)) != 0
            X[0] = np.where(mask, X[0] ^ P, X[0]).astype(np.uint32)
            t = np.where(mask, np.uint32(0), (X[0] ^ X[i]) & P).astype(np.uint32)
            X[0] ^= t
            X[i] ^= t
        Q >>= 1
    # Gray encode.
    for i in range(1, 3):
        X[i] ^= X[i - 1]
    t = np.zeros_like(X[0])
    Q = M
    while Q > 1:
        t = np.where((X[2] & np.uint32(Q)) != 0, t ^ np.uint32(Q - 1), t).astype(
            np.uint32
        )
        Q >>= 1
    for i in range(3):
        X[i] ^= t
    # Interleave transpose columns, most-significant bit first.
    d = np.zeros(len(q), dtype=np.uint64)
    for b in range(bits - 1, -1, -1):
        for i in range(3):
            d = (d << np.uint64(1)) | ((X[i] >> np.uint32(b)) & np.uint32(1)).astype(
                np.uint64
            )
    return d


def _centroids(mesh, total):
    """Per-cell centroids in global cell order, (total, 3), unused axes zero.

    Accumulates node coordinates in node order then divides — the same
    sequential summation the C++ core uses, so the doubles match exactly.
    """
    pts = np.asarray(mesh.points, dtype=np.float64)
    dim = min(pts.shape[1] if pts.ndim == 2 else 0, 3)
    cent = np.zeros((total, 3), dtype=np.float64)
    base = 0
    for cb in mesh.cells:
        data = np.asarray(cb.data)
        if data.ndim != 2:
            raise NotImplementedError(
                "partition numpy fallback: ragged/polyhedron blocks are C++-core only"
            )
        ncells, npc = data.shape
        acc = np.zeros((ncells, dim), dtype=np.float64)
        for k in range(npc):
            acc += pts[data[:, k], :dim]
        if npc > 0:
            acc /= float(npc)
        cent[base : base + ncells, :dim] = acc
        base += ncells
    return cent


def _sfc_labels(mesh, nparts, weights, total):
    labels = np.zeros(total, dtype=np.int64)
    if total == 0 or nparts == 1:
        return labels
    cent = _centroids(mesh, total)
    pts = np.asarray(mesh.points)
    dim = min(pts.shape[1] if pts.ndim == 2 else 0, 3)

    qmax = float((1 << _SFC_BITS) - 1)
    q = np.zeros((total, 3), dtype=np.uint32)
    for d in range(dim):
        lo = float(np.min(cent[:, d]))
        hi = float(np.max(cent[:, d]))
        span = hi - lo
        scale = qmax / span if span > 0.0 else 0.0
        t = (cent[:, d] - lo) * scale
        qi = (t + 0.5).astype(np.int64)
        np.clip(qi, 0, int(qmax), out=qi)
        q[:, d] = qi.astype(np.uint32)
    keys = _hilbert_keys(q)
    order = np.argsort(keys, kind="stable")

    if weights is None:
        ranks = np.empty(total, dtype=np.int64)
        ranks[order] = np.arange(total, dtype=np.int64)
        labels = (ranks * np.int64(nparts)) // np.int64(total)
    else:
        # Sequential prefix accumulation, matching the C++ core exactly.
        wsum = 0.0
        for w in weights:
            wsum += float(w)
        prefix = 0.0
        for r in range(total):
            w = float(weights[order[r]])
            p = int(np.floor((prefix + 0.5 * w) * nparts / wsum))
            labels[order[r]] = min(nparts - 1, max(0, p))
            prefix += w
    return labels


# --------------------------------------------------------------------------- #
# numpy KaHIP path (dual graph + the MIT `kahip` PyPI package)                 #
# --------------------------------------------------------------------------- #
def _dual_graph(mesh, total):
    """CSR dual graph: cells are vertices, edges where two cells share a face
    (3D) / an edge (2D). Cells without a facet table are isolated vertices."""
    from ._skin import _CELL_FACES
    from ._surface import _CELL_EDGES

    dual_dim = 0
    for cb in mesh.cells:
        dual_dim = max(dual_dim, topological_dimension.get(cb.type, 0))

    keys = []
    parents = []
    base = 0
    for cb in mesh.cells:
        data = np.asarray(cb.data)
        if data.ndim != 2:
            raise NotImplementedError(
                "partition numpy fallback: ragged/polyhedron blocks are C++-core only"
            )
        ncells = data.shape[0]
        if topological_dimension.get(cb.type, 0) == dual_dim:
            if dual_dim == 3:
                facets = _CELL_FACES.get(cb.type, [])
            elif dual_dim == 2:
                facets = _CELL_EDGES.get(cb.type, [])
            else:
                facets = []
            for _ftype, ncorners, nodes in facets:
                corners = np.sort(data[:, list(nodes[:ncorners])], axis=1)
                padded = np.full((ncells, 4), -1, dtype=np.int64)
                padded[:, :ncorners] = corners
                keys.append(padded)
                parents.append(np.arange(base, base + ncells, dtype=np.int64))
        base += ncells

    adj = [[] for _ in range(total)]
    if keys:
        all_keys = np.concatenate(keys)
        all_parents = np.concatenate(parents)
        _uniq, inverse = np.unique(all_keys, axis=0, return_inverse=True)
        first_owner = {}
        for rec in range(len(all_parents)):
            k = int(inverse[rec])
            p = int(all_parents[rec])
            if k not in first_owner:
                first_owner[k] = p
            elif first_owner[k] != p:
                adj[first_owner[k]].append(p)
                adj[p].append(first_owner[k])

    xadj = [0]
    adjncy = []
    for i in range(total):
        neigh = sorted(set(adj[i]))
        adjncy.extend(neigh)
        xadj.append(len(adjncy))
    return xadj, adjncy


_KAHIP_MODES = {"fast": 0, "eco": 1, "strong": 2}


def _kahip_labels(mesh, nparts, imbalance, mode, seed, weights, total):
    try:
        import kahip
    except ImportError:
        raise ValueError(
            "partition: method 'kahip' requires a meshioplusplus build with "
            "-DMESHIOPLUSPLUS_WITH_KAHIP=ON or the 'kahip' package "
            "(pip install kahip)"
        ) from None

    if total == 0:
        return np.zeros(0, dtype=np.int64)
    if nparts == 1:
        return np.zeros(total, dtype=np.int64)

    xadj, adjncy = _dual_graph(mesh, total)
    if weights is None:
        vwgt = [1] * total
    else:
        wsum = float(np.sum(weights))
        scale = float(1 << 30) / wsum
        vwgt = [max(1, int(round(float(w) * scale))) for w in weights]
    adjcwgt = [1] * len(adjncy)
    _edgecut, labels = kahip.kaffpa(
        vwgt,
        xadj,
        adjcwgt,
        adjncy,
        int(nparts),
        float(imbalance),
        True,  # suppress_output
        int(seed),
        _KAHIP_MODES[mode],
    )
    return np.asarray(labels, dtype=np.int64)


# --------------------------------------------------------------------------- #
# numpy piece builder (blocks kept 1:1, unlike split's subset builder)         #
# --------------------------------------------------------------------------- #
def _build_piece(mesh, flat, part, record_ids):
    n = len(mesh.points)
    used = np.zeros(n, dtype=bool)
    kept_per_block = []
    base = 0
    for cb in mesh.cells:
        data = np.asarray(cb.data)
        ncells = data.shape[0]
        kept = np.nonzero(flat[base : base + ncells] == part)[0]
        kept_per_block.append(kept)
        if len(kept):
            used[data[kept].reshape(-1)] = True
        base += ncells

    point_map = np.full(n, -1, dtype=np.int64)
    point_src = np.nonzero(used)[0]
    point_map[point_src] = np.arange(len(point_src), dtype=np.int64)

    pts = np.asarray(mesh.points)
    out_points = pts[point_src] if len(point_src) else pts[:0]

    new_cells = []
    cell_maps = []
    for cb, kept in zip(mesh.cells, kept_per_block):
        data = np.asarray(cb.data)
        cm = np.full(data.shape[0], -1, dtype=np.int64)
        cm[kept] = np.arange(len(kept), dtype=np.int64)
        conn = point_map[data[kept]] if len(kept) else data[:0].astype(np.int64)
        new_cells.append((cb.type, conn.astype(np.int64)))
        cell_maps.append(cm)

    point_data = {}
    for key, value in mesh.point_data.items():
        value = np.asarray(value)
        point_data[key] = value[point_src] if value.shape[:1] == (n,) else value.copy()
    if record_ids:
        point_data["partition:original_point_id"] = point_src.astype(np.int64)

    cell_data = {}
    for key, blk in mesh.cell_data.items():
        out_blk = []
        for b, kept in enumerate(kept_per_block):
            value = blk[b] if b < len(blk) else None
            out_blk.append(None if value is None else np.asarray(value)[kept])
        cell_data[key] = out_blk
    if record_ids:
        cell_data["partition:original_cell_id"] = [
            kept.astype(np.int64) for kept in kept_per_block
        ]

    field_data = {
        k: (v.copy() if isinstance(v, np.ndarray) else v)
        for k, v in mesh.field_data.items()
    }
    out = Mesh(
        out_points,
        new_cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
    )
    return out, point_map, cell_maps


def _cell_set_is_remappable(mesh, blocks):
    """Whether a ``cell_sets`` entry is per-block cell-index lists.

    Formats can stash *other* per-block side data in ``cell_sets`` — gmsh's
    ``gmsh:bounding_entities`` holds entity tags (negative values included),
    not cell indices — and remapping those as indices would crash or corrupt
    them.
    """
    for b, cb in enumerate(mesh.cells):
        if b >= len(blocks) or blocks[b] is None:
            continue
        idx = np.asarray(blocks[b]).reshape(-1)
        if len(idx) == 0:
            continue
        if idx.dtype.kind not in "iu":
            return False
        if idx.min() < 0 or idx.max() >= len(cb.data):
            return False
    return True


def _remap_piece_sets(src, out, point_map, cell_maps):
    """Remap the shim-only point_sets/cell_sets into a piece (all blocks kept).

    Entries that are not actually cell/point indices are copied through
    unchanged (the piece keeps the input's block structure 1:1, so per-block
    side data stays aligned).
    """
    point_sets = getattr(src, "point_sets", None)
    if point_sets:
        pm = np.asarray(point_map, dtype=np.int64)
        new_ps = {}
        for name, idx in point_sets.items():
            idx = np.asarray(idx)
            if idx.dtype.kind in "iu" and (
                len(idx) == 0 or (idx.min() >= 0 and idx.max() < len(pm))
            ):
                mapped = pm[idx.astype(np.int64)]
                new_ps[name] = mapped[mapped >= 0]
            else:
                new_ps[name] = idx.copy()
        out.point_sets = new_ps
    cell_sets = getattr(src, "cell_sets", None)
    if cell_sets and cell_maps is not None:
        new_cs = {}
        for name, blocks in cell_sets.items():
            if not _cell_set_is_remappable(src, blocks):
                new_cs[name] = [
                    None if b is None else np.asarray(b).copy() for b in blocks
                ]
                continue
            out_blocks = []
            for b in range(len(cell_maps)):
                if b < len(blocks) and blocks[b] is not None:
                    mapped = np.asarray(cell_maps[b], dtype=np.int64)[
                        np.asarray(blocks[b], dtype=np.int64).reshape(-1)
                    ]
                    out_blocks.append(mapped[mapped >= 0])
                else:
                    out_blocks.append(np.empty(0, dtype=np.int64))
            new_cs[name] = out_blocks
        if new_cs:
            out.cell_sets = new_cs


# --------------------------------------------------------------------------- #
# shared validation + label computation                                        #
# --------------------------------------------------------------------------- #
def _validate(mesh, nparts, method, mode, ghost_layers):
    if int(nparts) < 1:
        raise ValueError(f"partition: nparts must be >= 1, got {nparts}")
    if method not in ("auto", "sfc", "kahip"):
        raise ValueError(
            f"partition: unknown method '{method}' (expected 'sfc', 'kahip', or 'auto')"
        )
    if mode not in _KAHIP_MODES:
        raise ValueError(
            f"partition: unknown mode '{mode}' (expected 'fast', 'eco', or 'strong')"
        )
    if int(ghost_layers) < 0:
        raise ValueError(f"partition: ghost_layers must be >= 0, got {ghost_layers}")


def _weights_vector(mesh, weights, total):
    if weights is None:
        return None
    name = str(weights)
    if name not in mesh.cell_data:
        raise ValueError(f"partition: no cell_data array named '{name}'")
    blk = mesh.cell_data[name]
    if len(blk) != len(mesh.cells):
        raise ValueError(
            f"partition: weights cell_data '{name}' does not cover every cell block"
        )
    values = []
    for b, cb in enumerate(mesh.cells):
        a = np.asarray(blk[b], dtype=np.float64)
        if a.ndim > 1 and np.prod(a.shape[1:]) != 1:
            raise ValueError(
                f"partition: weights cell_data '{name}' must be scalar "
                "(one component per cell)"
            )
        a = a.reshape(-1)
        if len(a) != len(cb.data):
            raise ValueError(
                f"partition: weights cell_data '{name}' has {len(a)} rows on a "
                f"block of {len(cb.data)} cells"
            )
        values.append(a)
    w = np.concatenate(values) if values else np.zeros(0)
    if not np.all(np.isfinite(w)):
        raise ValueError(
            f"partition: weights cell_data '{name}' contains a non-finite value"
        )
    if np.any(w < 0):
        raise ValueError(
            f"partition: weights cell_data '{name}' contains a negative value"
        )
    if not float(np.sum(w)) > 0:
        raise ValueError(f"partition: weights cell_data '{name}' has zero total weight")
    return w


def _pip_kahip_available():
    try:
        import kahip  # noqa: F401

        return True
    except ImportError:
        return False


def _use_python_kahip(method):
    """Whether to bypass ``_core`` and run KaHIP via the ``kahip`` PyPI package.

    The C++ core wins when it has KaHIP compiled in. When it does not, a
    ``method="kahip"`` request (or an ``"auto"`` that could be upgraded) is
    served by the pure-Python dual graph + ``kahip.kaffpa`` — so pip-only users
    still get the quality path. Without the package either, ``"kahip"`` falls
    through to ``_core`` whose error names the build option.
    """
    if method not in ("kahip", "auto"):
        return False
    try:
        from . import _core

        if getattr(_core, "__has_kahip__", False):
            return False
    except Exception:
        pass  # no _core at all: the fallback handles every method anyway
    return _pip_kahip_available()


def _labels_py(mesh, nparts, method, imbalance, mode, seed, weights):
    total = sum(len(cb.data) for cb in mesh.cells)
    w = _weights_vector(mesh, weights, total)
    if method == "auto":
        method = "kahip" if _pip_kahip_available() else "sfc"
    if method == "kahip":
        if not 0.0 < float(imbalance) < 1.0:
            raise ValueError(f"partition: imbalance must be in (0, 1), got {imbalance}")
        return _kahip_labels(mesh, nparts, imbalance, mode, seed, w, total)
    return _sfc_labels(mesh, nparts, w, total)


def _split_flat_labels(mesh, flat):
    out = []
    base = 0
    for cb in mesh.cells:
        ncells = len(cb.data)
        out.append(np.asarray(flat[base : base + ncells], dtype=np.int64))
        base += ncells
    return out


# --------------------------------------------------------------------------- #
# public API                                                                   #
# --------------------------------------------------------------------------- #
def partition_labels(
    mesh,
    nparts: int,
    method: str = "auto",
    imbalance: float = 0.03,
    mode: str = "eco",
    seed: int = 0,
    weights: str | None = None,
) -> list:
    """Compute the per-cell part assignment without building the pieces.

    Parameters
    ----------
    mesh :
        the mesh to partition (unchanged).
    nparts :
        the number of parts (>= 1).
    method :
        ``"sfc"`` (Hilbert curve cut, always available), ``"kahip"`` (KaHIP
        ``kaffpa`` on the dual graph; needs a KaHIP-enabled C++ core or the
        ``kahip`` package), or ``"auto"`` (KaHIP when available, else SFC).
    imbalance :
        KaHIP only: allowed imbalance fraction in (0, 1).
    mode :
        KaHIP only: ``"fast"``/``"eco"``/``"strong"`` (default ``"eco"`` — per
        the Kratos KaHIPApplication benchmarks, eco/strong carry the quality).
    seed :
        KaHIP only: random seed (``kaffpa`` is deterministic per seed).
    weights :
        optional name of a scalar numeric ``cell_data`` array of per-cell
        weights; SFC cuts the weight prefix sum, KaHIP receives vertex weights.

    Returns
    -------
    list
        one Int64 numpy array per cell block (block-aligned, like a
        ``cell_data`` entry), values in ``[0, nparts)``.
    """
    _validate(mesh, nparts, method, mode, 0)
    if not _use_python_kahip(method):
        try:
            from . import _core

            return [
                np.asarray(a)
                for a in _core.partition_labels(
                    mesh,
                    int(nparts),
                    method,
                    float(imbalance),
                    mode,
                    int(seed),
                    weights or "",
                )
            ]
        except ValueError:
            raise
        except Exception:
            pass
    flat = _labels_py(mesh, int(nparts), method, imbalance, mode, seed, weights)
    return _split_flat_labels(mesh, flat)


def partition(
    mesh,
    nparts: int,
    method: str = "auto",
    imbalance: float = 0.03,
    mode: str = "eco",
    seed: int = 0,
    record_ids: bool = False,
    ghost_layers: int = 0,
    weights: str | None = None,
) -> list:
    """Decompose a mesh into exactly ``nparts`` balanced pieces.

    Same parameters as :func:`partition_labels`, plus:

    Parameters
    ----------
    record_ids :
        attach Int64 ``partition:original_point_id`` ``point_data`` and
        ``partition:original_cell_id`` ``cell_data`` (original input indices)
        to every piece.
    ghost_layers :
        grow each piece by this many shared-node BFS layers of other parts'
        cells (a halo), tagged ``partition:ghost``. Needs the C++ core.

    Returns
    -------
    list
        exactly ``nparts`` :class:`Mesh` pieces, part id = list index. Every
        piece keeps the input's cell-block structure 1:1 (empty blocks
        included, unlike :func:`split`), so concatenating the pieces
        reproduces the input. ``point_sets``/``cell_sets`` are remapped into
        each piece.

    Raises
    ------
    ValueError
        on ``nparts < 1``, ``ghost_layers < 0``, a bad weights array, an
        out-of-range imbalance, or ``method="kahip"`` when no KaHIP backend
        exists (the error names the build option).
    """
    _validate(mesh, nparts, method, mode, ghost_layers)
    pieces = None
    use_py_kahip = _use_python_kahip(method)
    if use_py_kahip and int(ghost_layers) != 0:
        # The pure-Python KaHIP route computes labels and builds disjoint
        # pieces; it cannot grow a halo. An EXPLICIT method="kahip" request
        # fails by name (the documented no-silent-SFC-downgrade rule), while
        # "auto" — a preference order, not a method — routes to the C++ core,
        # the one ghost-capable path.
        if method == "kahip":
            raise NotImplementedError(
                "partition: ghost_layers with method='kahip' needs KaHIP "
                "compiled into the C++ core (-DMESHIOPLUSPLUS_WITH_KAHIP=ON); "
                "the kahip PyPI route cannot build ghosted pieces"
            )
        use_py_kahip = False
    if not use_py_kahip:
        try:
            from . import _core

            raw = _core.partition(
                mesh,
                int(nparts),
                method,
                float(imbalance),
                mode,
                int(seed),
                bool(record_ids),
                int(ghost_layers),
                weights or "",
            )
            pieces = [
                (
                    p["mesh"],
                    np.asarray(p["point_map"]),
                    [np.asarray(c) for c in p["cell_maps"]],
                )
                for p in raw
            ]
        except ValueError:
            raise
        except Exception:
            pieces = None
    used_cpp = pieces is not None
    if pieces is None:
        # The numpy fallback builds disjoint pieces only. Ghost growth is a
        # second, different traversal (shared-node BFS over the whole mesh), and
        # silently returning unghosted pieces would be a wrong answer rather
        # than a slower one -- so say so instead, as _smooth.py does for its
        # inversion guard.
        if int(ghost_layers) != 0:
            raise NotImplementedError(
                "partition: ghost_layers needs the C++ core (the numpy fallback "
                "builds disjoint pieces only)"
            )
        flat = _labels_py(mesh, int(nparts), method, imbalance, mode, seed, weights)
        pieces = []
        for part in range(int(nparts)):
            out, pm, cm = _build_piece(mesh, flat, part, bool(record_ids))
            pieces.append((out, pm, cm))

    result = []
    for out, pm, cm in pieces:
        # The C++ core carries named regions onto each piece itself (and
        # therefore `point_sets`/`cell_sets`, which are views over them);
        # remapping again here would apply the maps twice.
        if not used_cpp:
            _remap_piece_sets(mesh, out, pm, cm)
        result.append(out)
    return result
