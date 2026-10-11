"""Split a surface into smooth patches by region growing.

The numpy twin of ``operations/segment.hpp``. A patch grows from its seed --
the lowest-numbered surface cell not yet in a patch -- breadth first across
every edge shared by exactly two cells that is smooth (dihedral at most
``feature_angle``, after reorienting an inconsistent pair) and, when enabled,
keeps the neighbour's normal within ``normal_tolerance`` of the seed's and its
curvature estimate within ``curvature_tolerance`` of the patch's mean so far.
Patches
smaller than ``min_patch_cells`` are merged, in one ascending pass, into the
neighbour they share the most edges with.

The result is byte-identical to the core's.
"""

from __future__ import annotations

import numpy as np

from ._fallback import core_op_declined
from ._feature_edges import _coords3, ring_unit_normals
from ._regions import block_bases

_PI = 3.14159265358979323846
_PREFIX = "meshio++: segment: "

#: Cell data (Int64): the patch of each cell, or -1.
PATCH_NAME = "segment:patch"
#: Field data (Int64, ``(m, 3)``): ``(a, b, shared edges)`` per adjacent pair.
ADJACENCY_NAME = "segment:adjacency"


def _selection(mesh, region, total):
    if not region:
        return np.ones(total, dtype=bool)
    for r in getattr(mesh, "regions", ()) or ():
        if r.name == region and r.kind == "cell":
            mask = np.zeros(total, dtype=bool)
            e = np.asarray(r.entries, dtype=np.int64).ravel()
            mask[e[(e >= 0) & (e < total)]] = True
            return mask
    names = sorted({r.name for r in getattr(mesh, "regions", ()) or ()})
    raise ValueError(
        f"{_PREFIX}no cell region named '{region}' "
        f"(available: {', '.join(names) if names else 'none'})"
    )


def surface_rings(mesh, selected, prefix=_PREFIX):
    """``(start, nodes, cell)``: the 2-D cells as face rings, with global cell ids.

    Volume blocks raise ``ValueError``; lines, vertices and cells naming a
    missing point are skipped.
    """
    from ._mesh import topological_dimension
    from ._surface import _CELL_EDGES

    n = len(mesh.points)
    bases = block_bases(mesh.cells)
    rings, cells = [], []
    for base, block in zip(bases, mesh.cells):
        ctype = block.type
        if ctype.startswith("polyhedron") or topological_dimension.get(ctype) == 3:
            raise ValueError(
                f"{prefix}block '{ctype}' is a volume block; segment a surface "
                "(take the skin first with extract_skin)"
            )
        if ctype.startswith("polygon"):
            rows = [np.asarray(r, dtype=np.int64) for r in block.data]
        else:
            defs = _CELL_EDGES.get(ctype)
            if topological_dimension.get(ctype) != 2 or not defs or len(defs) < 3:
                continue
            data = np.asarray(block.data, dtype=np.int64)
            if data.ndim != 2:
                continue
            corners = [nodes[0] for _name, _k, nodes in defs]
            rows = list(data[:, corners])
        for c, row in enumerate(rows):
            if not selected[base + c] or len(row) < 3:
                continue
            if row.min() < 0 or row.max() >= n:
                continue
            rings.append(row)
            cells.append(base + c)
    start = np.zeros(len(rings) + 1, dtype=np.int64)
    if rings:
        start[1:] = np.cumsum([len(r) for r in rings])
        nodes = np.concatenate(rings).astype(np.int64)
    else:
        nodes = np.empty(0, dtype=np.int64)
    return start, nodes, np.asarray(cells, dtype=np.int64)


def _edge_table(start, nodes, normals, feature_angle):
    """Sorted edge uses -> boundary/non-manifold counts and the paired edges."""
    nf = len(start) - 1
    sizes = np.diff(start)
    face_of = np.repeat(np.arange(nf, dtype=np.int64), sizes)
    pos = np.arange(len(nodes), dtype=np.int64)
    nxt = pos + 1
    if nf > 0:
        last = start[1:] - 1
        nxt[last[sizes > 0]] = start[:-1][sizes > 0]
    u = nodes
    v = nodes[nxt] if len(nodes) else nodes
    keep = u != v
    u, v, face = u[keep], v[keep], face_of[keep]
    lo, hi, fwd = np.minimum(u, v), np.maximum(u, v), u < v
    order = np.lexsort((fwd, face, hi, lo))
    lo, hi, face, fwd = lo[order], hi[order], face[order], fwd[order]
    pairs = {
        "a": np.empty(0, dtype=np.int64),
        "b": np.empty(0, dtype=np.int64),
        "inconsistent": np.empty(0, dtype=bool),
        "smooth": np.empty(0, dtype=bool),
        "angle": np.empty(0),
    }
    if len(lo) == 0:
        return 0, 0, 0, pairs
    head = np.ones(len(lo), dtype=bool)
    head[1:] = (lo[1:] != lo[:-1]) | (hi[1:] != hi[:-1])
    first = np.nonzero(head)[0]
    uses = np.diff(np.append(first, len(lo)))
    num_boundary = int((uses == 1).sum())
    num_non_manifold = int((uses >= 3).sum())
    two = first[uses == 2]
    two = two[face[two] != face[two + 1]]
    a, b = face[two], face[two + 1]
    inconsistent = fwd[two] == fwd[two + 1]
    nx, ny = normals[a], normals[b]
    zx = (nx[:, 0] == 0.0) & (nx[:, 1] == 0.0) & (nx[:, 2] == 0.0)
    zy = (ny[:, 0] == 0.0) & (ny[:, 1] == 0.0) & (ny[:, 2] == 0.0)
    dot = nx[:, 0] * ny[:, 0] + nx[:, 1] * ny[:, 1] + nx[:, 2] * ny[:, 2]
    dot = np.where(inconsistent, -dot, dot)
    defined = ~zx & ~zy
    cos_feature = np.cos(feature_angle * _PI / 180.0)
    smooth = ~defined | ~(dot < cos_feature)
    angle = np.where(defined, np.arccos(np.clip(dot, -1.0, 1.0)), np.nan)
    pairs = {
        "a": a,
        "b": b,
        "inconsistent": inconsistent,
        "smooth": smooth,
        "angle": angle,
    }
    num_feature = int((~smooth).sum())
    return num_boundary, num_non_manifold, num_feature, pairs


def _segment_py(mesh, opts):
    from ._regions import Region

    bases = block_bases(mesh.cells)
    total = int(bases[-1]) + len(mesh.cells[-1].data) if len(mesh.cells) else 0
    selected = _selection(mesh, opts["region"], total)
    start, nodes, cell_of = surface_rings(mesh, selected)
    nf = len(cell_of)
    xyz = _coords3(mesh)
    normals = ring_unit_normals(xyz, start, nodes)
    centroid = np.zeros((nf, 3))
    for f in range(nf):
        ring = nodes[start[f] : start[f + 1]]
        c = np.zeros(3)
        for k in ring:
            c = c + xyz[k]
        centroid[f] = c / float(len(ring))

    nb, nnm, nfeat, pairs = _edge_table(start, nodes, normals, opts["feature_angle"])
    pa, pb = pairs["a"].tolist(), pairs["b"].tolist()
    inc = pairs["inconsistent"].tolist()
    smooth = pairs["smooth"].tolist()
    angle = pairs["angle"].tolist()

    adj = [[] for _ in range(nf)]
    for e, (a, b) in enumerate(zip(pa, pb)):
        adj[a].append((b, e))
        adj[b].append((a, e))
    for lst in adj:
        lst.sort()

    curv = None
    if opts["curvature_tolerance"] > 0.0:
        curv = [0.0] * nf
        for e, (a, b) in enumerate(zip(pa, pb)):
            if not smooth[e] or angle[e] != angle[e]:
                continue
            d = centroid[a] - centroid[b]
            dist = float(np.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]))
            if not dist > 0.0:
                continue
            k = angle[e] / dist
            curv[a] = k if k > curv[a] else curv[a]
            curv[b] = k if k > curv[b] else curv[b]

    use_normal = opts["normal_tolerance"] > 0.0
    cos_normal = float(np.cos(opts["normal_tolerance"] * _PI / 180.0))
    ctol = opts["curvature_tolerance"]
    nrm = normals.tolist()
    patch = [-1] * nf
    flip = [0] * nf
    npatches = 0
    for seed in range(nf):
        if patch[seed] != -1:
            continue
        p = npatches
        npatches += 1
        patch[seed] = p
        flip[seed] = 0
        sn = nrm[seed]
        curv_sum = curv[seed] if curv is not None else 0.0
        curv_count = 1.0
        queue = [seed]
        q = 0
        while q < len(queue):
            g = queue[q]
            q += 1
            for nbr, e in adj[g]:
                if patch[nbr] != -1 or not smooth[e]:
                    continue
                fl = flip[g] ^ (1 if inc[e] else 0)
                if use_normal:
                    nn = nrm[nbr]
                    dot = sn[0] * nn[0] + sn[1] * nn[1] + sn[2] * nn[2]
                    if fl:
                        dot = -dot
                    if dot < cos_normal:
                        continue
                if curv is not None:
                    if abs(curv[nbr] - curv_sum / curv_count) > ctol:
                        continue
                    curv_sum += curv[nbr]
                    curv_count += 1.0
                patch[nbr] = p
                flip[nbr] = fl
                queue.append(nbr)

    label = list(range(npatches))
    num_merged = 0
    if opts["min_patch_cells"] > 1 and npatches > 1:
        size = [0] * npatches
        for f in range(nf):
            size[patch[f]] += 1
        shared = [{} for _ in range(npatches)]
        for a, b in zip(pa, pb):
            x, y = patch[a], patch[b]
            if x == y:
                continue
            shared[x][y] = shared[x].get(y, 0) + 1
            shared[y][x] = shared[y].get(x, 0) + 1
        for p in range(npatches):
            if label[p] != p or size[p] >= opts["min_patch_cells"] or not shared[p]:
                continue
            best, best_count = -1, -1
            for r in sorted(shared[p]):
                if shared[p][r] > best_count:
                    best, best_count = r, shared[p][r]
            for r, c in shared[p].items():
                shared[r].pop(p, None)
                if r == best:
                    continue
                shared[best][r] = shared[best].get(r, 0) + c
                shared[r][best] = shared[r].get(best, 0) + c
            shared[p] = {}
            size[best] += size[p]
            size[p] = 0
            label[p] = best
            num_merged += 1
        for p in range(npatches):
            r = p
            while label[r] != r:
                r = label[r]
            label[p] = r

    final = [-1] * npatches
    nfinal = 0
    for f in range(nf):
        r = label[patch[f]]
        if final[r] == -1:
            final[r] = nfinal
            nfinal += 1
        patch[f] = final[r]

    cell_patch = np.full(total, -1, dtype=np.int64)
    if nf:
        cell_patch[cell_of] = np.asarray(patch, dtype=np.int64)

    adjacency = {}
    for a, b in zip(pa, pb):
        x, y = patch[a], patch[b]
        if x == y:
            continue
        key = (min(x, y), max(x, y))
        adjacency[key] = adjacency.get(key, 0) + 1
    table = np.array(
        [[k[0], k[1], adjacency[k]] for k in sorted(adjacency)], dtype=np.int64
    ).reshape(-1, 3)

    out = mesh.copy()
    out.cell_data = dict(out.cell_data)
    out.cell_data[PATCH_NAME] = [
        cell_patch[base : base + len(block.data)].copy()
        for base, block in zip(bases, mesh.cells)
    ]
    out.field_data = dict(out.field_data)
    out.field_data[ADJACENCY_NAME] = table

    head = opts["prefix"] + ":"

    def _stale(r):
        tail = r.name[len(head) :]
        return r.kind == "cell" and r.name.startswith(head) and tail and tail.isdigit()

    regions = [r for r in (getattr(out, "regions", None) or []) if not _stale(r)]
    for p in range(nfinal):
        regions.append(
            Region(
                f"{head}{p}",
                "cell",
                np.nonzero(cell_patch == p)[0].astype(np.int64),
                2,
                -1,
            )
        )
    out.regions = regions
    report = {
        "num_patches": nfinal,
        "num_merged": num_merged,
        "num_feature_edges": nfeat,
        "num_boundary_edges": nb,
        "num_non_manifold_edges": nnm,
    }
    return out, report


_REPORT_KEYS = (
    "num_patches",
    "num_merged",
    "num_feature_edges",
    "num_boundary_edges",
    "num_non_manifold_edges",
)


def segment(
    mesh,
    feature_angle: float = 30.0,
    normal_tolerance: float = 0.0,
    curvature_tolerance: float = 0.0,
    min_patch_cells: int = 1,
    region: str = "",
    prefix: str = "patch",
    return_report: bool = False,
):
    """Split a surface into smooth patches by region growing.

    :param mesh: a surface mesh (triangles, quads, polygons and their quadratic
        variants, which contribute their corners); lines and vertices are
        carried and labelled -1. Volume blocks raise ``ValueError``.
    :param feature_angle: the largest dihedral angle, in degrees, still crossed;
        in ``[0, 180]``.
    :param normal_tolerance: the largest angle, in degrees, between a cell's
        normal and its patch seed's; ``0`` disables the test.
    :param curvature_tolerance: the largest difference between a cell's
        curvature estimate and its patch's mean so far, in inverse length units; ``0``
        disables the test.
    :param min_patch_cells: patches with fewer cells are merged into the
        neighbour they share the most edges with.
    :param region: restrict to the cells of this named cell region.
    :param prefix: the patch regions are named ``<prefix>:<k>``.
    :param return_report: also return the counts.
    :returns: the input with cell data ``segment:patch``, field data
        ``segment:adjacency`` (rows ``(a, b, shared edges)``) and one cell
        region per patch. With ``return_report``, ``(mesh, report)``.
    """
    opts = {
        "feature_angle": float(feature_angle),
        "normal_tolerance": float(normal_tolerance),
        "curvature_tolerance": float(curvature_tolerance),
        "min_patch_cells": int(min_patch_cells),
        "region": str(region),
        "prefix": str(prefix),
    }
    if not (0.0 <= opts["feature_angle"] <= 180.0):
        raise ValueError(f"{_PREFIX}the feature angle must lie in [0, 180] degrees")
    if not (0.0 <= opts["normal_tolerance"] <= 180.0):
        raise ValueError(f"{_PREFIX}the normal tolerance must lie in [0, 180] degrees")
    ctol = opts["curvature_tolerance"]
    if not (ctol >= 0.0 and np.isfinite(ctol)):
        raise ValueError(f"{_PREFIX}the curvature tolerance must be finite and >= 0")
    if opts["min_patch_cells"] < 1:
        raise ValueError(f"{_PREFIX}min_patch_cells must be >= 1")
    if not opts["prefix"]:
        raise ValueError(f"{_PREFIX}the region prefix must not be empty")

    out = None
    try:
        from . import _core

        res = _core.segment(
            mesh,
            opts["feature_angle"],
            opts["normal_tolerance"],
            opts["curvature_tolerance"],
            opts["min_patch_cells"],
            opts["region"],
            opts["prefix"],
        )
        out = res["mesh"]
        report = {k: int(res[k]) for k in _REPORT_KEYS}
    except Exception as exc:
        if not core_op_declined(exc, "segment"):
            raise
        out = None
    if out is None:
        out, report = _segment_py(mesh, opts)
    return (out, report) if return_report else out


__all__ = ["segment"]
