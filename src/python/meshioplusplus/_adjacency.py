"""A mesh's node or cell graph in compressed-sparse-row form, and the METIS and
Scotch graph files written from it (roadmap 7.2.2).

``adjacency(mesh, kind)`` returns ``(indptr, indices)``: the neighbours of
vertex ``i`` are ``indices[indptr[i]:indptr[i + 1]]``, ascending, with no
duplicates and no self-loops, and every undirected edge stored in both
directions. The kinds are ``operations/adjacency.hpp``'s:

- ``"node"``: nodes joined by a cell edge -- :func:`edge_index`'s graph in CSR
  form;
- ``"clique"``: nodes sharing any cell -- the sparsity pattern of an assembled
  finite-element matrix, without its diagonal;
- ``"cell"``: cells sharing a facet (a face in 3-D, an edge in 2-D), in the
  global block-major numbering -- the dual graph partitioners take.

The C++ core builds the graph; the numpy twins below reproduce it for a
pure-Python install, except the cell graph of a polyhedral mesh, which needs
the core. ``write_graph`` writes a graph (a mesh's, or any ``(indptr,
indices)`` pair) as a METIS ``.graph`` or a Scotch ``.grf`` file.
"""

from __future__ import annotations

import os

import numpy as np

from ._fallback import core_op_declined
from ._mesh import Mesh, topological_dimension

_KINDS = ("node", "clique", "cell")
_FORMATS = ("metis", "scotch")


def _check_kind(kind):
    if kind not in _KINDS:
        raise ValueError(
            f"meshio++: adjacency: unknown kind '{kind}' (expected 'node', 'clique' or 'cell')"
        )


def _is_rectangular(data):
    return isinstance(data, np.ndarray) and data.ndim == 2


def _csr_from_pairs(a, b, n):
    """Directed pairs -> CSR over ``n`` rows: in range, no loops, deduplicated,
    each row ascending."""
    a = np.asarray(a, dtype=np.int64)
    b = np.asarray(b, dtype=np.int64)
    keep = (a >= 0) & (a < n) & (b >= 0) & (b < n) & (a != b)
    a, b = a[keep], b[keep]
    order = np.lexsort((b, a))
    a, b = a[order], b[order]
    if a.size:
        uniq = np.ones(a.size, dtype=bool)
        uniq[1:] = (a[1:] != a[:-1]) | (b[1:] != b[:-1])
        a, b = a[uniq], b[uniq]
    indptr = np.zeros(n + 1, dtype=np.int64)
    indptr[1:] = np.cumsum(np.bincount(a, minlength=n))
    return indptr, np.ascontiguousarray(b, dtype=np.int64)


def _rings(block):
    """The node rings of a polygon block's rows, or of a polyhedron's faces."""
    if block.type.startswith("polyhedron"):
        for cell in block.data:
            for face in cell:
                yield np.asarray(face, dtype=np.int64).ravel()
    else:
        for row in block.data:
            yield np.asarray(row, dtype=np.int64).ravel()


def _node_graph(mesh):
    """The twin of ``detail::build_node_adjacency(..., Edge)`` plus the rings
    of uniform polygon blocks (``operations/adjacency.cpp``)."""
    from ._smooth import _REFINE_EDGES

    n = len(mesh.points)
    src, dst = [], []
    for block in mesh.cells:
        data = block.data
        if block.type in _REFINE_EDGES and _is_rectangular(data):
            data = np.asarray(data, dtype=np.int64)
            for a, b in _REFINE_EDGES[block.type]:
                src += [data[:, a], data[:, b]]
                dst += [data[:, b], data[:, a]]
        elif block.type.startswith(("polygon", "polyhedron")) or not _is_rectangular(
            data
        ):
            for ring in _rings(block):
                if ring.size >= 2:
                    nxt = np.roll(ring, -1)
                    src += [ring, nxt]
                    dst += [nxt, ring]
        # else: no tabulated edge topology (higher-order, Lagrange, custom).
    if not src:
        return np.zeros(n + 1, dtype=np.int64), np.zeros(0, dtype=np.int64)
    return _csr_from_pairs(np.concatenate(src), np.concatenate(dst), n)


def _clique_graph(mesh):
    """The twin of ``detail::build_node_adjacency(..., Clique)``."""
    n = len(mesh.points)
    src, dst = [], []
    for block in mesh.cells:
        data = block.data
        if _is_rectangular(data) and not block.type.startswith("polyhedron"):
            data = np.asarray(data, dtype=np.int64)
            k = data.shape[1]
            for i in range(k):
                for j in range(k):
                    if i != j:
                        src.append(data[:, i])
                        dst.append(data[:, j])
            continue
        for cell in data:
            if block.type.startswith("polyhedron"):
                nodes = np.unique(
                    np.concatenate(
                        [np.asarray(f, dtype=np.int64).ravel() for f in cell]
                    )
                    if len(cell)
                    else np.zeros(0, dtype=np.int64)
                )
            else:
                nodes = np.unique(np.asarray(cell, dtype=np.int64).ravel())
            if nodes.size >= 2:
                src.append(np.repeat(nodes, nodes.size))
                dst.append(np.tile(nodes, nodes.size))
    if not src:
        return np.zeros(n + 1, dtype=np.int64), np.zeros(0, dtype=np.int64)
    return _csr_from_pairs(np.concatenate(src), np.concatenate(dst), n)


def _block_dim(block):
    if block.type.startswith("polyhedron"):
        return 3
    if block.type.startswith("polygon") or not _is_rectangular(block.data):
        return 2
    return topological_dimension.get(block.type, 0)


def _cell_graph(mesh):
    """The twin of ``detail::build_facet_dual_graph(..., PolygonRings=true)``.

    A facet's owner is its first occurrence in block, cell, facet order -- the
    lowest global cell id -- and every other occurrence joins it, which is the
    core's rule for a non-manifold facet too."""
    from ._skin import _CELL_FACES
    from ._surface import _CELL_EDGES

    total = int(sum(len(block.data) for block in mesh.cells))
    dual_dim = max([_block_dim(block) for block in mesh.cells] + [0])
    keys, parents = [], []  # (rows, width) sorted corner ids, -1 padded
    base = 0
    for block in mesh.cells:
        ncells = len(block.data)
        if _block_dim(block) == dual_dim and dual_dim in (2, 3):
            if block.type.startswith("polyhedron"):
                raise NotImplementedError(
                    "meshio++: adjacency: the cell graph of a polyhedral mesh needs the C++ core"
                )
            if block.type.startswith("polygon") or not _is_rectangular(block.data):
                for c, ring in enumerate(_rings(block)):
                    if ring.size >= 1:
                        pairs = np.sort(
                            np.column_stack([ring, np.roll(ring, -1)]), axis=1
                        )
                        keys.append(pairs)
                        parents.append(np.full(len(pairs), base + c, dtype=np.int64))
            else:
                table = (_CELL_FACES if dual_dim == 3 else _CELL_EDGES).get(
                    block.type, []
                )
                if table:
                    data = np.asarray(block.data, dtype=np.int64)
                    # Cell-major, facet-minor: the core's record order.
                    per_facet = [
                        np.sort(data[:, list(nodes[:nc])], axis=1)
                        for _t, nc, nodes in table
                    ]
                    width = max(p.shape[1] for p in per_facet)
                    padded = np.full((ncells, len(table), width), -1, dtype=np.int64)
                    for f, p in enumerate(per_facet):
                        padded[:, f, : p.shape[1]] = p
                    keys.append(padded.reshape(ncells * len(table), width))
                    parents.append(
                        np.repeat(
                            np.arange(base, base + ncells, dtype=np.int64), len(table)
                        )
                    )
        base += ncells
    if not keys:
        return np.zeros(total + 1, dtype=np.int64), np.zeros(0, dtype=np.int64)
    width = max(k.shape[1] for k in keys)
    rows = np.concatenate(
        [np.pad(k, ((0, 0), (0, width - k.shape[1])), constant_values=-1) for k in keys]
    )
    parent = np.concatenate(parents)
    _, inverse = np.unique(rows, axis=0, return_inverse=True)
    inverse = np.asarray(inverse).ravel()
    owner = np.full(int(inverse.max()) + 1, np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(owner, inverse, parent)
    own = owner[inverse]
    joined = own != parent
    a, b = own[joined], parent[joined]
    return _csr_from_pairs(np.concatenate([a, b]), np.concatenate([b, a]), total)


def _adjacency_py(mesh, kind):
    if kind == "node":
        return _node_graph(mesh)
    if kind == "clique":
        return _clique_graph(mesh)
    return _cell_graph(mesh)


def adjacency(mesh, kind: str = "node"):
    """The mesh's graph in compressed-sparse-row form.

    :param mesh: the mesh (never modified).
    :param kind: ``"node"`` (nodes joined by a cell edge), ``"clique"`` (nodes
        sharing a cell: the FEM sparsity pattern) or ``"cell"`` (cells sharing
        a facet, numbered block-major: the dual graph).
    :returns: ``(indptr, indices)``, int64 arrays: the neighbours of vertex
        ``i`` are ``indices[indptr[i]:indptr[i + 1]]``, ascending, with every
        undirected edge stored both ways, so ``indptr[-1]`` is twice the edge
        count. ``scipy.sparse.csr_matrix((np.ones(len(indices)), indices,
        indptr))`` is the adjacency matrix.
    """
    _check_kind(kind)
    try:
        from . import _core

        indptr, indices = _core.adjacency(mesh, kind)
        return np.asarray(indptr, dtype=np.int64), np.asarray(indices, dtype=np.int64)
    except Exception as exc:
        if not core_op_declined(exc, "adjacency"):
            raise
    return _adjacency_py(mesh, kind)


def _format_for(path, file_format):
    if file_format is not None:
        if file_format not in _FORMATS:
            raise ValueError(
                f"meshio++: write_graph: unknown format '{file_format}' "
                "(expected 'metis' or 'scotch')"
            )
        return file_format
    ext = os.path.splitext(os.fspath(path))[1].lower()
    if ext in (".graph", ".metis"):
        return "metis"
    if ext == ".grf":
        return "scotch"
    raise ValueError(
        f"meshio++: write_graph: cannot infer the format of '{os.fspath(path)}' "
        "(use .graph or .metis for METIS, .grf for Scotch, or name the format)"
    )


def _check_graph(indptr, indices):
    """The core's ``adjgraph_check``: an undirected simple graph."""
    if indptr.ndim != 1 or indices.ndim != 1:
        raise ValueError("meshio++: write_graph: indptr and indices must be 1-D arrays")
    if indptr.size == 0:
        raise ValueError(
            "meshio++: write_graph: indptr is empty (it holds n + 1 offsets)"
        )
    n = indptr.size - 1
    if indptr[0] != 0 or indptr[-1] != indices.size:
        raise ValueError(
            "meshio++: write_graph: indptr must start at 0 and end at the number of "
            f"indices ({indices.size})"
        )
    degree = np.diff(indptr)
    if np.any(degree < 0):
        raise ValueError(
            f"meshio++: write_graph: indptr decreases at row {int(np.argmax(degree < 0))}"
        )
    rows = np.repeat(np.arange(n, dtype=np.int64), degree)
    if indices.size and (indices.min() < 0 or indices.max() >= n):
        raise ValueError(
            f"meshio++: write_graph: a row names a vertex outside the {n} vertices"
        )
    if np.any(rows == indices):
        raise ValueError(
            "meshio++: write_graph: a row has a self-loop, which a graph file cannot hold"
        )
    fwd = np.lexsort((indices, rows))
    a, b = rows[fwd], indices[fwd]
    if a.size > 1 and np.any((a[1:] == a[:-1]) & (b[1:] == b[:-1])):
        raise ValueError("meshio++: write_graph: a row lists a neighbour twice")
    rev = np.lexsort((rows, indices))
    if not (np.array_equal(a, indices[rev]) and np.array_equal(b, rows[rev])):
        raise ValueError(
            "meshio++: write_graph: the graph is not symmetric (METIS and Scotch files "
            "are undirected)"
        )
    return n, int(degree.max()) if n else 0


def _write_graph_py(path, indptr, indices, file_format):
    n, max_degree = _check_graph(indptr, indices)
    nnz = int(indices.size)
    shift = 1 if file_format == "metis" else 0
    lines = (
        [f"{n} {nnz // 2}"] if file_format == "metis" else ["0", f"{n} {nnz}", "0 000"]
    )
    for i in range(n):
        row = indices[indptr[i] : indptr[i + 1]] + shift
        words = [str(int(v)) for v in row]
        if file_format == "scotch":
            words.insert(0, str(len(row)))
        lines.append(" ".join(words))
    with open(path, "w", encoding="ascii", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")
    return {
        "format": file_format,
        "num_vertices": n,
        "num_edges": nnz // 2,
        "max_degree": max_degree,
    }


def write_graph(path, graph, kind: str = "node", *, file_format=None):
    """Write a graph as a METIS or Scotch graph file.

    :param path: the file to write.
    :param graph: a :class:`Mesh`, whose ``kind`` graph is written (see
        :func:`adjacency`), or an ``(indptr, indices)`` pair, which must be an
        undirected simple graph: symmetric, without self-loops or repeats.
    :param kind: the graph of a mesh: ``"node"``, ``"clique"`` or ``"cell"``.
    :param file_format: ``"metis"`` (``.graph``: a header ``n m``, then each
        vertex's 1-based neighbours) or ``"scotch"`` (``.grf``: ``0``, ``n
        2m``, ``0 000``, then each vertex's degree and 0-based neighbours);
        inferred from the extension (``.graph``/``.metis``, ``.grf``) when
        omitted.
    :returns: ``{"format", "num_vertices", "num_edges", "max_degree"}``.
    """
    fmt = _format_for(path, file_format)
    if isinstance(graph, Mesh):
        indptr, indices = adjacency(graph, kind)
    else:
        try:
            indptr, indices = graph
        except (TypeError, ValueError):
            raise TypeError(
                "meshio++: write_graph: graph must be a Mesh or an (indptr, indices) pair"
            ) from None
    indptr = np.ascontiguousarray(indptr, dtype=np.int64)
    indices = np.ascontiguousarray(indices, dtype=np.int64)
    path = os.fspath(path)
    try:
        from . import _core

        return dict(_core.write_graph(path, indptr, indices, fmt))
    except Exception as exc:
        if not core_op_declined(exc, "write_graph"):
            raise
    return _write_graph_py(path, indptr, indices, fmt)
