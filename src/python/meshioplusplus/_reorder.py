"""Mesh renumbering / reordering (RCM + Morton/Hilbert space-filling curves).

A dependency-free mesh *operation* (not a file format): it permutes node and
element indices to reduce sparse-matrix bandwidth and improve cache locality,
as a pure permutation — no nodes or cells added or removed, geometry and all
data preserved. The C++ core (``_core.reorder``) does the work; this module is
the thin shim (try C++, fall back to a pure-numpy reference) plus the
``point_sets``/``cell_sets`` remapping the core cannot see across the Python
boundary.

Public API:

* :func:`reorder` — renumber a mesh, optionally returning the applied
  node/cell permutations (old -> new).
* :func:`compute_bandwidth` — connectivity bandwidth, a before/after measure.
"""

from __future__ import annotations

import numpy as np

from ._mesh import Mesh

_SFC_BITS = 21


# --------------------------------------------------------------------------- #
# cell-node helpers                                                           #
# --------------------------------------------------------------------------- #
def _flatten_cell_nodes(cell):
    """Return the node ids of one (possibly ragged/polyhedron) cell as a list."""
    arr = np.asarray(cell, dtype=object) if isinstance(cell, list) else cell
    if (
        isinstance(cell, list)
        and cell
        and isinstance(cell[0], (list, tuple, np.ndarray))
    ):
        # polyhedron: list of faces (each a list of node ids)
        out = []
        for face in cell:
            out.extend(int(x) for x in np.asarray(face).reshape(-1))
        return out
    return [int(x) for x in np.asarray(arr).reshape(-1)]


def _iter_cell_node_lists(blk, n):
    """Yield the in-range node-id list of every cell of a block."""
    data = blk.data
    if isinstance(data, np.ndarray) and data.ndim == 2:
        for row in data:
            yield [int(x) for x in row if 0 <= int(x) < n]
    else:
        for cell in data:
            yield [x for x in _flatten_cell_nodes(cell) if 0 <= x < n]


# --------------------------------------------------------------------------- #
# RCM                                                                          #
# --------------------------------------------------------------------------- #
def _build_adjacency(mesh, n):
    adj = [set() for _ in range(n)]
    for blk in mesh.cells:
        for nodes in _iter_cell_node_lists(blk, n):
            for a in nodes:
                sa = adj[a]
                for b in nodes:
                    if a != b:
                        sa.add(b)
    return [sorted(s) for s in adj]


def _rls(adj, root, level):
    """Rooted level structure (BFS): return (eccentricity, sorted last level).

    ``level`` is a length-n list of -1; it is reset to all -1 before return.
    """
    from collections import deque

    touched = [root]
    level[root] = 0
    q = deque([root])
    ecc = 0
    while q:
        u = q.popleft()
        lu = level[u]
        if lu > ecc:
            ecc = lu
        for v in adj[u]:
            if level[v] == -1:
                level[v] = lu + 1
                touched.append(v)
                q.append(v)
    last = sorted(t for t in touched if level[t] == ecc)
    for t in touched:
        level[t] = -1
    return ecc, last


def _pseudo_peripheral(adj, start, degree, level):
    root = start
    ecc, last = _rls(adj, root, level)
    for _ in range(50):
        cand = last[0] if last else root
        best = degree[cand]
        for v in last:
            if degree[v] < best:
                best = degree[v]
                cand = v
        ecc2, last2 = _rls(adj, cand, level)
        if ecc2 > ecc:
            root, ecc, last = cand, ecc2, last2
        else:
            break
    return root


def _rcm_perm(mesh, n):
    if n == 0:
        return np.empty(0, dtype=np.int64)
    from collections import deque

    adj = _build_adjacency(mesh, n)
    degree = [len(a) for a in adj]
    visited = bytearray(n)
    order = []
    level = [-1] * n
    for s in range(n):
        if visited[s]:
            continue
        root = _pseudo_peripheral(adj, s, degree, level)
        visited[root] = 1
        q = deque([root])
        while q:
            u = q.popleft()
            order.append(u)
            nbrs = [v for v in adj[u] if not visited[v]]
            nbrs.sort(key=lambda v: (degree[v], v))
            for v in nbrs:
                if not visited[v]:
                    visited[v] = 1
                    q.append(v)
    order.reverse()  # Cuthill-McKee -> reverse
    perm = np.empty(n, dtype=np.int64)
    perm[np.asarray(order, dtype=np.int64)] = np.arange(n, dtype=np.int64)
    return perm


# --------------------------------------------------------------------------- #
# space-filling curves                                                         #
# --------------------------------------------------------------------------- #
def _part1by2(x):
    x &= 0x1FFFFF
    x = (x | (x << 32)) & 0x1F00000000FFFF
    x = (x | (x << 16)) & 0x1F0000FF0000FF
    x = (x | (x << 8)) & 0x100F00F00F00F00F
    x = (x | (x << 4)) & 0x10C30C30C30C30C3
    x = (x | (x << 2)) & 0x1249249249249249
    return x


def _morton_key(q):
    return (
        _part1by2(int(q[0])) | (_part1by2(int(q[1])) << 1) | (_part1by2(int(q[2])) << 2)
    )


def _hilbert_key(q, bits):
    X = [int(q[0]), int(q[1]), int(q[2])]
    n = 3
    M = 1 << (bits - 1)
    Q = M
    while Q > 1:
        P = Q - 1
        for i in range(n):
            if X[i] & Q:
                X[0] ^= P
            else:
                t = (X[0] ^ X[i]) & P
                X[0] ^= t
                X[i] ^= t
        Q >>= 1
    for i in range(1, n):
        X[i] ^= X[i - 1]
    t = 0
    Q = M
    while Q > 1:
        if X[n - 1] & Q:
            t ^= Q - 1
        Q >>= 1
    for i in range(n):
        X[i] ^= t
    d = 0
    for b in range(bits - 1, -1, -1):
        for i in range(n):
            d = (d << 1) | ((X[i] >> b) & 1)
    return d


def _sfc_perm(mesh, n, hilbert):
    if n == 0:
        return np.empty(0, dtype=np.int64)
    pts = np.asarray(mesh.points, dtype=np.float64)
    if pts.ndim == 1:
        pts = pts.reshape(n, -1)
    pdim = pts.shape[1] if pts.ndim == 2 else 1
    dim = min(pdim, 3)
    bits = _SFC_BITS
    q = np.zeros((n, 3), dtype=np.uint64)
    maxq = (1 << bits) - 1
    for d in range(dim):
        col = pts[:, d]
        lo = float(col.min())
        hi = float(col.max())
        span = hi - lo
        scale = maxq / span if span > 0 else 0.0
        qi = np.floor((col - lo) * scale + 0.5).astype(np.int64)
        qi = np.clip(qi, 0, maxq)
        q[:, d] = qi.astype(np.uint64)
    if hilbert:
        keys = np.array([_hilbert_key(q[i], bits) for i in range(n)], dtype=object)
    else:
        keys = np.array([_morton_key(q[i]) for i in range(n)], dtype=object)
    order = np.argsort(keys, kind="stable")
    perm = np.empty(n, dtype=np.int64)
    perm[order] = np.arange(n, dtype=np.int64)
    return perm


# --------------------------------------------------------------------------- #
# applying a permutation                                                       #
# --------------------------------------------------------------------------- #
def _scatter(arr, perm):
    """``out[perm[i]] = arr[i]`` along the first axis, preserving dtype/shape."""
    out = np.empty_like(arr)
    out[perm] = arr
    return out


def _remap_ragged_cell(cell, node_perm):
    if (
        isinstance(cell, list)
        and cell
        and isinstance(cell[0], (list, tuple, np.ndarray))
    ):
        return [
            [int(node_perm[int(x)]) for x in np.asarray(face).reshape(-1)]
            for face in cell
        ]
    return [int(node_perm[int(x)]) for x in np.asarray(cell).reshape(-1)]


def _apply_perm(mesh, node_perm):
    n = len(mesh.points)
    node_perm = np.asarray(node_perm, dtype=np.int64)

    new_points = (
        _scatter(np.asarray(mesh.points), node_perm) if n else np.asarray(mesh.points)
    )

    new_cells = []
    cell_perms = []
    for blk in mesh.cells:
        data = blk.data
        if isinstance(data, np.ndarray) and data.ndim == 2:
            nc = data.shape[0]
            remapped = node_perm[data] if nc else data.astype(np.int64)
            if nc:
                key = remapped.min(axis=1)
                order = np.argsort(key, kind="stable")
            else:
                order = np.empty(0, dtype=np.int64)
            cell_perm = np.empty(nc, dtype=np.int64)
            cell_perm[order] = np.arange(nc, dtype=np.int64)
            new_cells.append((blk.type, remapped[order].astype(np.int64)))
        else:
            cells_list = list(data)
            nc = len(cells_list)
            remapped_cells = [_remap_ragged_cell(c, node_perm) for c in cells_list]
            keys = np.array(
                [
                    min(_flatten_cell_nodes(rc)) if _flatten_cell_nodes(rc) else 0
                    for rc in remapped_cells
                ],
                dtype=np.int64,
            )
            order = (
                np.argsort(keys, kind="stable") if nc else np.empty(0, dtype=np.int64)
            )
            cell_perm = np.empty(nc, dtype=np.int64)
            cell_perm[order] = np.arange(nc, dtype=np.int64)
            new_cells.append((blk.type, [remapped_cells[o] for o in order]))
        cell_perms.append(cell_perm)

    point_data = {}
    for key, value in mesh.point_data.items():
        value = np.asarray(value)
        if value.shape[:1] == (n,):
            point_data[key] = _scatter(value, node_perm)
        else:
            point_data[key] = value.copy()

    cell_data = {}
    for key, blocks in mesh.cell_data.items():
        out_blocks = []
        for b, value in enumerate(blocks):
            if value is None:
                out_blocks.append(None)
                continue
            value = np.asarray(value)
            if b < len(cell_perms) and value.shape[:1] == (len(cell_perms[b]),):
                out_blocks.append(_scatter(value, cell_perms[b]))
            else:
                out_blocks.append(value.copy())
        cell_data[key] = out_blocks

    field_data = {
        k: (v.copy() if isinstance(v, np.ndarray) else v)
        for k, v in mesh.field_data.items()
    }

    out = Mesh(
        new_points,
        new_cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
    )
    return out, node_perm, cell_perms


def _reorder_py(mesh, method):
    m = str(method).lower()
    n = len(mesh.points)
    if m == "rcm":
        perm = _rcm_perm(mesh, n)
    elif m in ("morton", "z", "zorder", "z-order"):
        perm = _sfc_perm(mesh, n, hilbert=False)
    elif m == "hilbert":
        perm = _sfc_perm(mesh, n, hilbert=True)
    else:
        raise ValueError(
            f"reorder: unknown method '{method}' (expected 'rcm', 'morton', or 'hilbert')"
        )
    return _apply_perm(mesh, perm)


# --------------------------------------------------------------------------- #
# set remapping (Python-side side-channels the C++ core cannot see)            #
# --------------------------------------------------------------------------- #
def _remap_sets(src, out, node_perm, cell_perms):
    point_sets = getattr(src, "point_sets", None)
    if point_sets:
        np_perm = np.asarray(node_perm, dtype=np.int64)
        out.point_sets = {
            name: np_perm[np.asarray(idx, dtype=np.int64)]
            for name, idx in point_sets.items()
        }
    cell_sets = getattr(src, "cell_sets", None)
    if cell_sets and cell_perms is not None:
        new_cs = {}
        for name, blocks in cell_sets.items():
            remapped = []
            for b, idx in enumerate(blocks):
                if idx is None:
                    remapped.append(None)
                elif b < len(cell_perms):
                    remapped.append(
                        np.asarray(cell_perms[b])[np.asarray(idx, dtype=np.int64)]
                    )
                else:
                    remapped.append(np.asarray(idx))
            new_cs[name] = remapped
        out.cell_sets = new_cs


# --------------------------------------------------------------------------- #
# public API                                                                   #
# --------------------------------------------------------------------------- #
def reorder(mesh, method: str = "rcm", return_permutation: bool = False):
    """Renumber ``mesh`` to reduce bandwidth / improve locality.

    Parameters
    ----------
    mesh :
        the mesh to renumber.
    method :
        ``"rcm"`` (Reverse Cuthill-McKee, the default), ``"morton"`` (Z-order),
        or ``"hilbert"`` (Hilbert space-filling curve).
    return_permutation :
        when true, also return the applied node/cell permutations (old -> new).

    Returns
    -------
    Mesh, or (Mesh, node_permutation, cell_permutations)
        The renumbered mesh (a pure permutation of the input). ``point_data``
        follows its nodes, ``cell_data`` follows its cells, ``field_data`` is
        preserved, and ``point_sets``/``cell_sets`` indices are remapped.
        ``node_permutation`` is an int64 array (old index -> new index);
        ``cell_permutations`` is one int64 array per cell block.

    Raises
    ------
    ValueError
        if ``method`` is not one of ``"rcm"``, ``"morton"``, ``"hilbert"``.
    """
    out = None
    node_perm = None
    cell_perms = None
    try:
        from . import _core

        res = _core.reorder(mesh, method)
        out = res["mesh"]
        node_perm = np.asarray(res["node_permutation"])
        cell_perms = [np.asarray(c) for c in res["cell_permutations"]]
    except ValueError:
        raise
    except Exception:
        out = None
    used_cpp = out is not None
    if out is None:
        out, node_perm, cell_perms = _reorder_py(mesh, method)

    # The C++ core carries named regions across itself (and therefore
    # `point_sets`/`cell_sets`, which are views over them), so remapping again
    # here would apply the maps twice. Only the numpy fallback needs it.
    if not used_cpp:
        _remap_sets(mesh, out, node_perm, cell_perms)

    if return_permutation:
        return out, node_perm, cell_perms
    return out


def _compute_bandwidth_py(mesh):
    bw = 0
    for blk in mesh.cells:
        data = blk.data
        if isinstance(data, np.ndarray) and data.ndim == 2 and data.size:
            bw = max(bw, int((data.max(axis=1) - data.min(axis=1)).max()))
        else:
            for cell in data:
                nodes = _flatten_cell_nodes(cell)
                if nodes:
                    bw = max(bw, max(nodes) - min(nodes))
    return int(bw)


def compute_bandwidth(mesh) -> int:
    """Connectivity bandwidth: max over all cells of (max - min) node index."""
    try:
        from . import _core

        return int(_core.compute_bandwidth(mesh))
    except Exception:
        pass
    return _compute_bandwidth_py(mesh)
