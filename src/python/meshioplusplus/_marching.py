"""The marching-tetrahedra cutter shared by ``_slice`` and ``_isosurface``.

Private module -- the pure-numpy twin of ``cpp/include/meshioplusplus/detail/marching.hpp``
(hoisted out of ``_slice.py`` in v7.15.0 alongside the C++ hoist; slice's output
is unchanged, byte for byte). Nothing about marching tetrahedra is specific to a
plane: the algorithm cuts wherever a **per-node signed scalar** changes sign.
``slice`` supplies ``d_i = dot(x_i - origin, normal)``, ``isosurface`` supplies
``d_i = f_i - isovalue``.

Written to reproduce the C++ result **byte for byte** -- the traversal order, the
``t = d_lo/(d_lo - d_hi)`` crossing formula from the sorted edge endpoints, the
Newell-normal orientation flip and the bbox-relative degeneracy tolerance --
pinned by ``tests/python/test_slice.py::test_cpp_matches_python`` and
``tests/python/test_isosurface.py::test_cpp_matches_python``.

Two phases, mirroring the C++ header:

* :func:`_marching_prepare` -- everything independent of the cut scalar
  (simplexify, per-block input cell bases, the max topological dimension). The
  caller resolves its scalar against the **simplexified** mesh, which is what
  lets ``isosurface`` contour a ``point_data`` array that rode through the
  linearization.
* :func:`_marching_cut` -- one cut, returning the section mesh. ``isosurface``
  calls it once per isovalue against the same prepared input.
"""

from __future__ import annotations

import numpy as np

from ._mesh import Mesh

# Marching-tetrahedra / -triangle tables (twins of the C++ constants in
# detail/marching.cpp). Local tetra edges e0..e5, local triangle edges e0..e2.
_TET_EDGE = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
_TRI_EDGE = ((0, 1), (1, 2), (0, 2))

# Crossing-edge ring per below-zero sign mask (bit i set iff d_i < 0). A single
# isolated vertex crosses 3 edges (triangle); a 2-2 split crosses 4 (quad).
_TET_RING = (
    (),  # 0x0
    (0, 1, 2),  # 0x1
    (0, 3, 4),  # 0x2
    (2, 1, 3, 4),  # 0x3
    (1, 3, 5),  # 0x4
    (0, 2, 5, 3),  # 0x5
    (0, 4, 5, 1),  # 0x6
    (2, 4, 5),  # 0x7
    (2, 4, 5),  # 0x8
    (0, 4, 5, 1),  # 0x9
    (0, 2, 5, 3),  # 0xA
    (1, 3, 5),  # 0xB
    (2, 1, 3, 4),  # 0xC
    (0, 3, 4),  # 0xD
    (0, 1, 2),  # 0xE
    (),  # 0xF
)
_TRI_RING = (
    (),  # 0x0
    (0, 2),  # 0x1
    (0, 1),  # 0x2
    (1, 2),  # 0x3
    (1, 2),  # 0x4
    (0, 1),  # 0x5
    (0, 2),  # 0x6
    (),  # 0x7
)

_DIM_OF = {
    "tetra": 3,
    "hexahedron": 3,
    "wedge": 3,
    "pyramid": 3,
    "triangle": 2,
    "quad": 2,
}

_REL_TOL = 1e-10

#: Wind every section face toward a fixed direction (slice's rule).
FIXED_DIRECTION = "fixed"

#: Wind every section face toward increasing field (isosurface's rule).
FIELD_GRADIENT = "gradient"


def _cell_dim(cell_type):
    return _DIM_OF.get(cell_type, 0)


class _MarchingInput:
    """Everything the cutter needs that does not depend on the cut scalar."""

    __slots__ = ("simp", "orig_block_base", "maxdim", "points", "dim")

    def __init__(self, simp, orig_block_base, maxdim, points, dim):
        self.simp = simp
        self.orig_block_base = orig_block_base
        self.maxdim = maxdim
        self.points = points
        self.dim = dim


def _marching_prepare(mesh):
    """Simplexify ``mesh`` and collect everything independent of the cut."""
    from ._convert_cells import _simplexify_py

    # Simplexify first: every 3D cell -> tetra, every 2D cell -> triangle. This
    # is byte-identical to the C++ convert_cells(simplexify) the core calls.
    simp, _, _ = _simplexify_py(mesh, record_parent_ids=True)
    pts = np.asarray(simp.points, dtype=np.float64)
    dim = pts.shape[1] if pts.ndim == 2 else 0

    # Original (input) global cell base per block (blocks map 1:1).
    orig_block_base = []
    base = 0
    for cb in mesh.cells:
        orig_block_base.append(base)
        base += len(cb.data)

    # Max topological dimension decides the cut: volume -> surface, 2D -> lines.
    maxdim = 0
    for cb in simp.cells:
        maxdim = max(maxdim, _cell_dim(cb.type))

    return _MarchingInput(simp, orig_block_base, maxdim, pts, dim)


def _pad3(pts, node, dim):
    """A node's coordinates padded to three components (z = 0 in 2D)."""
    return np.array(
        (
            [pts[node, 0], pts[node, 1], pts[node, 2]]
            if dim == 3
            else [pts[node, 0], pts[node, 1], 0.0]
        ),
        dtype=np.float64,
    )


def _marching_cut(
    prep,
    dist,
    provenance_name,
    record_parent_ids,
    orientation=FIXED_DIRECTION,
    direction=(0.0, 0.0, 1.0),
):
    """Cut the prepared mesh where ``dist`` changes sign.

    ``dist`` carries one signed value per point of ``prep.simp``; the cut is its
    zero level set, with zero itself counting as **positive** (so the sign mask
    is total and a grazing cut is emitted exactly once).
    """
    simp = prep.simp
    pts = prep.points
    dim = prep.dim
    maxdim = prep.maxdim
    nblocks = len(simp.cells)
    direction = np.asarray(direction, dtype=np.float64).reshape(3)
    field_gradient = orientation == FIELD_GRADIENT

    has_parent = "convert:parent_cell" in simp.cell_data

    edge_tbl = _TET_EDGE if maxdim == 3 else _TRI_EDGE
    ring_tbl = _TET_RING if maxdim == 3 else _TRI_RING
    ncorners = 4 if maxdim == 3 else 3

    # --- serial cut + edge-key dedup (the determinism pin) -------------------
    node_id = {}
    node_edges = []  # section node id -> (lo, hi) simplexified global node ids
    faces = []  # (nverts, [node ids], parent_block, parent_local, parent_global, ref)

    if maxdim in (2, 3):
        want = "tetra" if maxdim == 3 else "triangle"
        for block, cb in enumerate(simp.cells):
            if cb.type != want:
                continue
            data = np.asarray(cb.data)
            if data.ndim != 2:
                continue
            parent_arr = (
                np.asarray(simp.cell_data["convert:parent_cell"][block])
                if has_parent
                else None
            )
            for c in range(len(data)):
                corners = [int(data[c, k]) for k in range(ncorners)]
                mask = 0
                for k in range(ncorners):
                    if dist[corners[k]] < 0.0:
                        mask |= 1 << k
                ring = ring_tbl[mask]
                if not ring:
                    continue
                within = int(parent_arr[c]) if parent_arr is not None else c
                ref = None
                if field_gradient:
                    # Reference direction = centroid(corners with d >= 0) minus
                    # centroid(corners with d < 0), i.e. "toward increasing
                    # field". Summed in ascending corner index and divided once,
                    # so this matches the C++ arithmetic bit for bit.
                    pos = np.zeros(3, dtype=np.float64)
                    neg = np.zeros(3, dtype=np.float64)
                    npos = 0
                    nneg = 0
                    for k in range(ncorners):
                        x = _pad3(pts, corners[k], dim)
                        if mask & (1 << k):
                            neg = neg + x
                            nneg += 1
                        else:
                            pos = pos + x
                            npos += 1
                    # Both are non-empty whenever the ring is non-empty.
                    ref = pos * (1.0 / float(npos)) - neg * (1.0 / float(nneg))
                nodes = []
                for ei in ring:
                    ga = corners[edge_tbl[ei][0]]
                    gb = corners[edge_tbl[ei][1]]
                    key = (ga, gb) if ga < gb else (gb, ga)
                    nid = node_id.get(key)
                    if nid is None:
                        nid = len(node_edges)
                        node_id[key] = nid
                        node_edges.append(key)
                    nodes.append(nid)
                faces.append(
                    (
                        len(ring),
                        nodes,
                        block,
                        c,
                        prep.orig_block_base[block] + within,
                        ref,
                    )
                )

    nnodes = len(node_edges)

    # --- crossing-point coordinates ------------------------------------------
    out_points = np.zeros((nnodes, dim), dtype=np.float64)
    for i in range(nnodes):
        lo, hi = node_edges[i]
        dl = dist[lo]
        dh = dist[hi]
        t = dl / (dl - dh)
        for k in range(dim):
            pl = pts[lo, k]
            ph = pts[hi, k]
            out_points[i, k] = pl + t * (ph - pl)

    # --- point_data: interpolate every source array at the same t (Float64) --
    point_data = {}
    for name in sorted(simp.point_data):
        a = np.asarray(simp.point_data[name])
        src_rows = a.shape[0] if a.ndim >= 1 else 0
        if src_rows == 0:
            continue
        comp_shape = a.shape[1:]
        flat = a.reshape(src_rows, -1).astype(np.float64)
        comp = flat.shape[1]
        out_arr = np.zeros((nnodes, comp), dtype=np.float64)
        for i in range(nnodes):
            lo, hi = node_edges[i]
            dl = dist[lo]
            dh = dist[hi]
            t = dl / (dl - dh)
            for k in range(comp):
                vl = flat[lo, k]
                vh = flat[hi, k]
                out_arr[i, k] = vl + t * (vh - vl)
        point_data[name] = out_arr.reshape((nnodes,) + comp_shape)

    if not faces:
        return Mesh(out_points, [], point_data=point_data)

    # --- bbox-relative degeneracy tolerance ----------------------------------
    lo = out_points.min(axis=0)
    hi = out_points.max(axis=0)
    bbdiag = float(np.sqrt(float(np.sum((hi - lo) * (hi - lo)))))
    area_tol = _REL_TOL * bbdiag * bbdiag
    len_tol = _REL_TOL * bbdiag

    # --- stage output blocks (fixed order: triangle, quad, line) -------------
    blocks = {
        "triangle": {"npc": 3, "conn": [], "pb": [], "pl": [], "pg": []},
        "quad": {"npc": 4, "conn": [], "pb": [], "pl": [], "pg": []},
        "line": {"npc": 2, "conn": [], "pb": [], "pl": [], "pg": []},
    }
    for nverts, nodes, pblk, plocal, pglobal, ref in faces:
        if nverts >= 3:
            ring = [_pad3(out_points, nodes[v], dim) for v in range(nverts)]
            nrm = np.zeros(3, dtype=np.float64)
            for i in range(nverts):
                a = ring[i]
                b = ring[(i + 1) % nverts]
                nrm[0] += (a[1] - b[1]) * (a[2] + b[2])
                nrm[1] += (a[2] - b[2]) * (a[0] + b[0])
                nrm[2] += (a[0] - b[0]) * (a[1] + b[1])
            if (
                float(np.sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]))
                < area_tol
            ):
                continue
            d = ref if field_gradient else direction
            flip = float(nrm[0] * d[0] + nrm[1] * d[1] + nrm[2] * d[2]) < 0.0
            blk = blocks["triangle"] if nverts == 3 else blocks["quad"]
            order = range(nverts - 1, -1, -1) if flip else range(nverts)
            for v in order:
                blk["conn"].append(nodes[v])
            blk["pb"].append(pblk)
            blk["pl"].append(plocal)
            blk["pg"].append(pglobal)
        else:  # line segment
            a = nodes[0]
            b = nodes[1]
            d2 = 0.0
            for k in range(dim):
                dd = out_points[a, k] - out_points[b, k]
                d2 += dd * dd
            if float(np.sqrt(d2)) < len_tol:
                continue
            blk = blocks["line"]
            blk["conn"].append(nodes[0])
            blk["conn"].append(nodes[1])
            blk["pb"].append(pblk)
            blk["pl"].append(plocal)
            blk["pg"].append(pglobal)

    order = [name for name in ("triangle", "quad", "line") if blocks[name]["conn"]]
    cells = []
    for name in order:
        blk = blocks[name]
        conn = np.array(blk["conn"], dtype=np.int64).reshape(-1, blk["npc"])
        cells.append((name, conn))

    out = Mesh(out_points, cells, point_data=point_data)

    # --- cell_data: replicate each parent's row to its section cells ----------
    cell_data = {}
    for name in sorted(simp.cell_data):
        if name == "convert:parent_cell":
            continue
        blk_list = simp.cell_data[name]
        if len(blk_list) != nblocks:
            continue
        blk0 = np.asarray(blk_list[0])
        out_blocks = []
        for name_out in order:
            b = blocks[name_out]
            ncells = len(b["pg"])
            arr = np.empty((ncells,) + blk0.shape[1:], dtype=blk0.dtype)
            for c in range(ncells):
                arr[c] = np.asarray(blk_list[b["pb"][c]])[b["pl"][c]]
            out_blocks.append(arr)
        cell_data[name] = out_blocks
    if cell_data:
        out.cell_data = cell_data

    if record_parent_ids:
        out.cell_data[provenance_name] = [
            np.array(blocks[name_out]["pg"], dtype=np.int64) for name_out in order
        ]

    return out
