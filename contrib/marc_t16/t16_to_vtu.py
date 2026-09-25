# MIT License -- part of meshio++ (https://github.com/loumalouomega/meshioplusplus).
"""Export an MSC Marc binary post file (``.t16``) to a ParaView series meshio++
reads.

Run with the Python that ships with Marc/Mentat, whose PyPost module
(``py_post``) reads post files::

    python t16_to_vtu.py job.t16 [job.pvd] [--increments 1,5,10]

``mio_vtu_series.py`` (from ``contrib/common``) must sit next to this script.
The output is ``job.pvd`` plus a ``job/`` directory of ``.vtu`` files, one per
increment, at the increment's time (``p.time``; the increment number where two
increments share a time):

- the undeformed coordinates, with the node vectors (``Displacement``, ...) as
  point data and the node scalars too;
- element scalars and tensors averaged over the element's nodes as cell data (a
  tensor as six components 11, 22, 33, 12, 23, 13);
- node and element ids as ``marc:node_id`` / ``marc:element_id``, and every set
  as a ``set:<name>`` 0/1 mask.

Element types map as the native ``.t19`` reader maps them (Volume B; the node
lists are in meshio++'s order, the extra nodes of Herrmann, bubble and
generalized-plane-strain types dropped); others are skipped with a count. The
formatted ``.t19`` twin of a post file is read natively (``meshioplusplus.read``).
See ``doc/routes/marc_t16.md``.
"""

from __future__ import division, print_function

import argparse
import os
import sys

import numpy as np

# mio_vtu_series.py: copied next to this script, or in the repository's
# contrib/common.
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_HERE, os.path.join(os.path.dirname(_HERE), "common")]

from mio_vtu_series import PvdSeries  # noqa: E402

# Marc element type -> (VTK cell type, the cell's node count).
MARC_TYPES = {}
for _types, _cell in (
    ((3, 10, 11, 18, 75, 139, 140, 143, 144, 145, 147, 151, 152), (9, 4)),
    ((2, 6, 138, 158, 201), (5, 3)),
    (
        (22, 26, 27, 28, 30, 32, 33, 46, 48, 53, 54, 55, 58, 59, 63, 66, 142, 148)
        + (153, 154),
        (23, 8),
    ),
    ((124, 125, 126, 128, 129, 200), (22, 6)),
    ((7, 43, 117, 123, 146, 149), (12, 8)),
    ((21, 23, 35, 44, 57, 61, 150), (25, 20)),
    ((134, 135), (10, 4)),
    ((127, 130, 133), (24, 10)),
    ((136, 137), (13, 6)),
    ((9, 31, 52, 98, 165, 166, 167), (3, 2)),
    ((64, 168, 169, 170), (21, 3)),
    # More nodes than the cell keeps: the leading geometric nodes are the cell.
    ((80, 81, 82, 83, 118, 119), (9, 4)),
    ((34, 47, 60), (23, 8)),
    ((84, 120), (12, 8)),
    ((155, 156), (5, 3)),
    ((157,), (10, 4)),
):
    for _t in _types:
        MARC_TYPES[_t] = _cell


def increment_mesh(p):
    """Points, cell blocks and ids of the current increment."""
    nn = p.nodes()
    node_ids = np.array([p.node_id(i) for i in range(nn)], dtype=np.int64)
    index = dict((int(n), i) for i, n in enumerate(node_ids))
    points = np.zeros((nn, 3))
    for i in range(nn):
        node = p.node(i)
        points[i] = (node.x, node.y, getattr(node, "z", 0.0))
    by_type, order, element_rows, skipped = {}, [], {}, {}
    for e in range(p.elements()):
        element = p.element(e)
        if element.type not in MARC_TYPES:
            skipped[element.type] = skipped.get(element.type, 0) + 1
            continue
        vtk, n = MARC_TYPES[element.type]
        if vtk not in by_type:
            by_type[vtk] = []
            order.append(vtk)
        by_type[vtk].append(e)
        element_rows[e] = [index[int(node)] for node in list(element.items)[:n]]
    for etype, count in sorted(skipped.items()):
        print(
            "t16_to_vtu: %d elements of type %d have no cell, skipped" % (count, etype)
        )
    blocks, kept = [], []
    for vtk in order:
        rows = by_type[vtk]
        blocks.append((vtk, np.array([element_rows[e] for e in rows], dtype=np.int64)))
        kept.extend(rows)
    element_ids = np.array([p.element_id(e) for e in kept], dtype=np.int64)
    return points, blocks, node_ids, kept, element_ids


def _node_average(values):
    return float(np.mean([v.value for v in values])) if values else np.nan


def _tensor_average(tensors):
    if not tensors:
        return [np.nan] * 6
    names = ("t11", "t22", "t33", "t12", "t23", "t13")
    return [float(np.mean([getattr(t, k) for t in tensors])) for k in names]


def increment_data(p, nn, kept):
    point, cell = {}, {}
    for k in range(p.node_scalars()):
        point[p.node_scalar_label(k)] = np.array(
            [p.node_scalar(i, k) for i in range(nn)], dtype=np.float64
        )
    for k in range(p.node_vectors()):
        rows = []
        for i in range(nn):
            v = p.node_vector(i, k)
            rows.append((v.x, v.y, getattr(v, "z", 0.0)))
        point[p.node_vector_label(k)] = np.array(rows, dtype=np.float64)
    for k in range(p.element_scalars()):
        cell[p.element_scalar_label(k)] = np.array(
            [_node_average(p.element_scalar(e, k)) for e in kept], dtype=np.float64
        )
    for k in range(p.element_tensors()):
        cell[p.element_tensor_label(k)] = np.array(
            [_tensor_average(p.element_tensor(e, k)) for e in kept], dtype=np.float64
        )
    return point, cell


def set_masks(p, node_ids, element_ids):
    point, cell = {}, {}
    count = p.sets() if hasattr(p, "sets") else 0
    for k in range(count):
        s = p.set(k)
        kind = str(s.type).lower()
        members = set(int(i) for i in s.items)
        if kind.startswith("node"):
            point["set:" + s.name] = np.array(
                [1 if int(n) in members else 0 for n in node_ids], dtype=np.uint8
            )
        elif kind.startswith("element"):
            cell["set:" + s.name] = np.array(
                [1 if int(e) in members else 0 for e in element_ids], dtype=np.uint8
            )
    return point, cell


def export(p, out, increments=None):
    count = p.increments()
    # Position 0 holds the model; the results start at 1 (a post file with a
    # single position holds only the model, exported as it is).
    positions = list(range(1, count)) if count > 1 else [0]
    if increments:
        positions = [i for i in positions if i in increments]
    last = None
    with PvdSeries(out) as series:
        for position in positions:
            p.moveto(position)
            points, blocks, node_ids, kept, element_ids = increment_mesh(p)
            point, cell = increment_data(p, len(points), kept)
            point["marc:node_id"] = node_ids
            cell["marc:element_id"] = element_ids
            point_sets, cell_sets = set_masks(p, node_ids, element_ids)
            point.update(point_sets)
            cell.update(cell_sets)
            t = float(getattr(p, "time", position))
            if last is not None and t <= last:
                t = float(position) if position > last else last + 1.0
            last = t
            series.add(t, points, blocks, point_data=point, cell_data=cell)
    return len(positions)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("t16")
    parser.add_argument("out", nargs="?", help="the .pvd to write (default: <t16>.pvd)")
    parser.add_argument(
        "--increments", help="comma-separated post file positions (default: all)"
    )
    args = parser.parse_args(argv)
    out = args.out or os.path.splitext(args.t16)[0] + ".pvd"
    wanted = (
        set(int(i) for i in args.increments.split(",")) if args.increments else None
    )

    import py_post

    p = py_post.post_open(args.t16)
    try:
        count = export(p, out, wanted)
    finally:
        p.close()
    print("t16_to_vtu: wrote %d increments to %s" % (count, out))


if __name__ == "__main__":
    main()
