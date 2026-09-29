"""Carry ``side`` regions through an operation: the numpy twin of the
"Side facets by containment" half of ``detail/region_remap.cpp``.

A side entry ``(global cell, local facet)`` whose cell keeps its identity -- one
output cell of the same type, reached by no other input cell, whose facet at
that number still has the same nodes -- keeps its number. Any other entry is
found again by what the facet is made of:

- **refined** (several children): every child facet lying within it -- each of
  the child facet's nodes the image of one of the facet's, or a point the
  operation created that lies on it;
- **merged or retyped** (one child that others share, or of another type): the
  output facet containing it -- each surviving node one of that facet's, and a
  node the operation removed lying on it (a second pass also accepts a
  surviving node that only lies on it).

When neither finds it, a facet whose corners all survive is looked up across
the whole output by those corners (a decimation collapse hands a boundary edge
to a neighbouring cell) -- for a ``FirstChild`` map only, never a merge.

Node identity decides wherever it can; geometry (a relative 1e-9 plane or
segment test) only places points an operation created or removed. KEEP IN SYNC
with the C++ core: the fallbacks that use this must agree with it entry for
entry.
"""

from __future__ import annotations

import numpy as np

from ._common import warn
from ._facets import FacetIndex
from ._regions import Region, block_bases
from ._skin import _CELL_FACES
from ._surface import _CELL_EDGES


def _block_row(bases, g):
    b = int(np.searchsorted(bases, g, side="right")) - 1
    if b < 0 or b >= len(bases) - 1 or g >= bases[-1]:
        return None, None
    return b, g - int(bases[b])


def _facet_count(mesh, bases, g):
    b, row = _block_row(bases, g)
    if b is None:
        return 0
    block = mesh.cells[b]
    if isinstance(block.data, list):
        cell = block.data[row]
        return len(cell)  # a polygon's edges, or a polyhedron's faces
    table = _CELL_FACES.get(block.type) or _CELL_EDGES.get(block.type)
    return len(table) if table else 0


def _facet(mesh, bases, g, f):
    """``(nodes, num_corners)`` of facet ``f`` of global cell ``g``, or None."""
    b, row = _block_row(bases, g)
    if b is None or f < 0:
        return None
    block = mesh.cells[b]
    if isinstance(block.data, list):
        cell = block.data[row]
        if f >= len(cell):
            return None
        if block.type.startswith("polyhedron"):
            nodes = [int(v) for v in cell[f]]
        else:
            nodes = [int(cell[f]), int(cell[(f + 1) % len(cell)])]
        return (nodes, len(nodes)) if len(nodes) >= 2 else None
    table = _CELL_FACES.get(block.type) or _CELL_EDGES.get(block.type)
    if not table or f >= len(table):
        return None
    _, ncorner, local = table[f]
    data = block.data[row]
    return [int(data[i]) for i in local], ncorner


def _coords(mesh, p):
    x = np.zeros(3)
    pts = np.asarray(mesh.points)
    d = min(pts.shape[1], 3)
    x[:d] = pts[p, :d]
    return x


def _on_facet(q, corners):
    n = len(corners)
    if n < 2:
        return False
    size = max(float(np.linalg.norm(c - corners[0])) for c in corners[1:])
    if size == 0.0:
        return False
    tol = 1e-9 * size
    if n == 2:
        e = corners[1] - corners[0]
        w = q - corners[0]
        t = float(np.dot(w, e) / np.dot(e, e))
        if t < -1e-9 or t > 1.0 + 1e-9:
            return False
        return float(np.linalg.norm(w - t * e)) <= tol
    normal = np.zeros(3)
    for i in range(n):
        a, b = corners[i], corners[(i + 1) % n]
        normal[0] += (a[1] - b[1]) * (a[2] + b[2])
        normal[1] += (a[2] - b[2]) * (a[0] + b[0])
        normal[2] += (a[0] - b[0]) * (a[1] + b[1])
    length = float(np.linalg.norm(normal))
    if length == 0.0:
        return False
    centre = sum(corners) / n
    return abs(float(np.dot(q - centre, normal)) / length) <= tol


def _within(small_mesh, small, big_mesh, big, small_to_big, geometric):
    big_nodes, big_corners = big
    corner_xyz = None
    for node in small[0]:
        there = small_to_big(node)
        if there >= 0 and there in big_nodes:
            continue
        if there >= 0 and not geometric:
            return False
        if corner_xyz is None:
            corner_xyz = [_coords(big_mesh, v) for v in big_nodes[:big_corners]]
        if not _on_facet(_coords(small_mesh, node), corner_xyz):
            return False
    return True


def first_child_children(src, out, cell_maps):
    """``children(g)`` for a ``FirstChild`` map (a parent's children are a run
    ending at the next non-negative entry, or at the block's end)."""
    in_bases = block_bases(src.cells)
    out_bases = block_bases(out.cells)

    def children(g):
        b, row = _block_row(in_bases, g)
        if b is None or b >= len(cell_maps) or b >= len(out.cells):
            return []
        m = np.asarray(cell_maps[b], dtype=np.int64)
        first = int(m[row]) if row < len(m) else -1
        if first < 0:
            return []
        n_out = int(out_bases[b + 1] - out_bases[b])
        later = m[row + 1 :]
        later = later[later >= 0]
        end = int(later[0]) if len(later) else n_out
        return [int(out_bases[b]) + k for k in range(first, min(end, n_out))]

    return children


def direct_children(src, out, cell_maps):
    """``children(g)`` for a ``Direct`` map (``map[c]`` is the output cell)."""
    in_bases = block_bases(src.cells)
    out_bases = block_bases(out.cells)

    def children(g):
        b, row = _block_row(in_bases, g)
        if b is None or b >= len(cell_maps) or b >= len(out.cells):
            return []
        m = np.asarray(cell_maps[b], dtype=np.int64)
        c = int(m[row]) if row < len(m) else -1
        n_out = int(out_bases[b + 1] - out_bases[b])
        return [int(out_bases[b]) + c] if 0 <= c < n_out else []

    return children


def carry_side_regions(
    src, out, children, point_map=None, op_name="operation", many_to_one=False
):
    """Append ``src``'s side regions, carried, to ``out.regions``.

    ``children(g)`` gives the output global cells of input global cell ``g``;
    ``point_map`` maps input points to output points (-1 dropped; ``None`` the
    identity). ``many_to_one`` says several input cells may share an output
    cell (a ``Direct`` merge such as ``undo_green``).
    """
    sides = [r for r in getattr(src, "regions", []) if r.kind == "side"]
    if not sides:
        return
    in_bases = block_bases(src.cells)
    out_bases = block_bases(out.cells)
    n_out_points = len(out.points)
    pm = None if point_map is None else np.asarray(point_map, dtype=np.int64)

    def in_to_out(p):
        q = p if pm is None else (int(pm[p]) if 0 <= p < len(pm) else -1)
        return q if q < n_out_points else -1

    inv = None

    def out_to_in(q):
        return int(inv[q]) if 0 <= q < len(inv) else -1

    preimages = None
    out_facets = None
    new = []
    for region in sides:
        entries = np.asarray(region.entries, dtype=np.int64).reshape(-1, 2)
        out_entries = []
        lost = 0
        for cell, facet in entries.tolist():
            before = len(out_entries)
            kids = children(cell)
            source = _facet(src, in_bases, cell, facet)
            if source is None:
                lost += 1
                continue
            if not kids:
                pass  # the cell is gone; only the last resort below can place it
            elif len(kids) == 1:
                child = kids[0]
                if preimages is None and many_to_one:
                    preimages = np.zeros(int(out_bases[-1]), dtype=np.int64)
                    for g in range(int(in_bases[-1])):
                        for k in children(g):
                            preimages[k] += 1
                one_to_one = preimages is None or preimages[child] == 1
                in_type = src.cells[_block_row(in_bases, cell)[0]].type
                out_type = out.cells[_block_row(out_bases, child)[0]].type
                target = _facet(out, out_bases, child, facet)
                if (
                    one_to_one
                    and in_type == out_type
                    and target is not None
                    and len(target[0]) == len(source[0])
                    and _within(src, source, out, target, in_to_out, False)
                ):
                    out_entries.append((child, facet))
                    continue
                nf = _facet_count(out, out_bases, child)
                for geometric in (False, True):
                    for j in range(nf):
                        target = _facet(out, out_bases, child, j)
                        if target is not None and _within(
                            src, source, out, target, in_to_out, geometric
                        ):
                            out_entries.append((child, j))
                    if len(out_entries) > before:
                        break
            else:
                if inv is None:
                    inv = np.full(n_out_points, -1, dtype=np.int64)
                    for p in range(len(src.points)):
                        q = in_to_out(p)
                        if q >= 0 and inv[q] < 0:
                            inv[q] = p
                for child in kids:
                    for j in range(_facet_count(out, out_bases, child)):
                        target = _facet(out, out_bases, child, j)
                        if target is not None and _within(
                            out, target, src, source, out_to_in, False
                        ):
                            out_entries.append((child, j))
            if len(out_entries) == before:
                # Last resort: the facet still exists, by its surviving
                # corners, on another output cell.
                corners = [in_to_out(v) for v in source[0][: source[1]]]
                if (
                    not many_to_one
                    and all(q >= 0 for q in corners)
                    and len(corners) <= 4
                ):
                    if out_facets is None:
                        out_facets = FacetIndex(out)
                    hit = out_facets.find(corners)
                    if hit is not None:
                        out_entries.append(hit.first)
            if len(out_entries) == before:
                lost += 1
        if lost:
            warn(
                f"{op_name}: side region '{region.name}' lost {lost} of its "
                f"{len(entries)} facet(s) -- they have no counterpart in the output"
            )
        new.append(
            Region(
                region.name,
                "side",
                np.asarray(out_entries, dtype=np.int64).reshape(-1, 2),
                region.dim,
                region.tag,
            )
        )
    out.regions = list(out.regions) + new


def carry_regions_to_facet_mesh(src, out, out_parent, out_facet, point_map, op_name):
    """The twin of ``detail::carry_regions_to_facet_mesh``: a side region
    becomes a cell region of the facet mesh, a point region follows
    ``point_map``, a cell region is dropped with a warning."""
    regions = list(getattr(src, "regions", []))
    if not regions:
        return
    where = {
        (int(c), int(f)): i
        for i, (c, f) in enumerate(zip(np.asarray(out_parent), np.asarray(out_facet)))
    }
    pm = np.asarray(point_map, dtype=np.int64)
    new = []
    for region in regions:
        if region.kind == "point":
            e = np.asarray(region.entries, dtype=np.int64)
            e = e[(e >= 0) & (e < len(pm))]
            mapped = pm[e]
            new.append(
                Region(
                    region.name, "point", mapped[mapped >= 0], region.dim, region.tag
                )
            )
        elif region.kind == "side":
            pairs = np.asarray(region.entries, dtype=np.int64).reshape(-1, 2).tolist()
            hits = [where.get((c, f), -1) for c, f in pairs]
            lost = sum(1 for h in hits if h < 0)
            if lost:
                warn(
                    f"{op_name}: side region '{region.name}' lost {lost} of its "
                    f"{len(pairs)} facet(s) -- they are not on the extracted boundary"
                )
            new.append(
                Region(
                    region.name,
                    "cell",
                    np.asarray([h for h in hits if h >= 0], dtype=np.int64),
                    region.dim,
                    region.tag,
                )
            )
        else:
            warn(
                f"{op_name}: cell region '{region.name}' dropped -- it names input "
                "cells, and the output holds only their facets"
            )
    out.regions = new
