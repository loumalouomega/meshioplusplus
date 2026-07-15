"""
I/O for the I-DEAS Universal file format (``.unv``), following the datasets used
by FEconv <https://github.com/victorsndvg/FEconv>: 2411 (nodes), 2412
(elements, with the FE-descriptor -> element-type map and the Salome/Code-Aster
mid-node "sandwich" ordering for parabolic elements) and 2467 (permanent
groups, mapped to point/cell sets).

A UNV file is a sequence of datasets, each delimited by a line whose content is
``-1``, followed by a dataset-id line and the dataset records.
"""

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._files import open_file
from .._mesh import CellBlock, Mesh

__all__ = ["read", "write"]

# Salome/UNV parabolic node order -> meshio position (0-based).  For UNV node i
# (0-based within the element), meshio_conn[nd[i]] = unv_conn[i].
_ND = {
    "line3": [0, 2, 1],
    "triangle6": [0, 3, 1, 4, 2, 5],
    "quad8": [0, 4, 1, 5, 2, 6, 3, 7],
    "tetra10": [0, 4, 1, 5, 2, 6, 7, 8, 9, 3],
    "hexahedron20": [
        0,
        8,
        1,
        9,
        2,
        10,
        3,
        11,
        12,
        13,
        14,
        15,
        4,
        16,
        5,
        17,
        6,
        18,
        7,
        19,
    ],
}

# UNV FE descriptor id -> meshio cell type
_unv_to_meshio_type = {
    11: "line",
    21: "line",
    22: "line3",
    24: "line3",
    41: "triangle",
    81: "triangle",
    91: "triangle",
    42: "triangle6",
    82: "triangle6",
    92: "triangle6",
    44: "quad",
    84: "quad",
    94: "quad",
    122: "quad",
    45: "quad8",
    85: "quad8",
    95: "quad8",
    111: "tetra",
    118: "tetra10",
    112: "wedge",
    115: "hexahedron",
    116: "hexahedron20",
}

# meshio type -> (descriptor id, is_beam) used on write
_meshio_to_unv = {
    "line": (21, True),
    "line3": (24, True),
    "triangle": (91, False),
    "triangle6": (92, False),
    "quad": (94, False),
    "quad8": (95, False),
    "tetra": (111, False),
    "tetra10": (118, False),
    "wedge": (112, False),
    "hexahedron": (115, False),
    "hexahedron20": (116, False),
}
_beam_descriptors = {11, 21, 22, 24}


def _split_datasets(lines):
    """Yield (dataset_id, body_lines) for each -1 delimited dataset."""
    i, n = 0, len(lines)
    while i < n:
        if lines[i].strip() == "-1":
            i += 1
            if i >= n:
                break
            ds_id = int(lines[i].strip())
            i += 1
            body = []
            while i < n and lines[i].strip() != "-1":
                body.append(lines[i])
                i += 1
            i += 1  # skip closing -1
            yield ds_id, body
        else:
            i += 1


def read(filename):
    with open_file(filename, "r") as f:
        lines = f.read().splitlines()

    points = []
    point_label_to_index = {}
    # cells grouped by type, in first-appearance order
    cell_groups = {}  # meshio_type -> list of connectivity rows
    cell_pid = {}  # meshio_type -> list of property ids
    elem_label_to_ref = {}  # element label -> (meshio_type, index within type)
    node_label_to_index = {}
    groups = []  # (name, entity_type, [tags])

    for ds_id, body in _split_datasets(lines):
        if ds_id in (2411, 781):
            k = 0
            idx = 0
            while k + 1 < len(body):
                rec1 = body[k].split()
                if not rec1:
                    k += 1
                    continue
                label = int(rec1[0])
                coords = body[k + 1].replace("D", "E").replace("d", "e").split()
                points.append([float(c) for c in coords])
                point_label_to_index[label] = idx
                node_label_to_index[label] = idx
                idx += 1
                k += 2
        elif ds_id == 2412:
            _read_2412(
                body, cell_groups, cell_pid, elem_label_to_ref, point_label_to_index
            )
        elif ds_id in (2467, 2477):
            groups.extend(_read_2467(body))
        # other datasets (fields 2414, ...) are ignored on read for now

    points = np.array(points, dtype=float) if points else np.empty((0, 3))

    cells = []
    cell_data = {"unv:pid": []}
    type_order = list(cell_groups.keys())
    for t in type_order:
        cells.append(CellBlock(t, np.array(cell_groups[t], dtype=int)))
        cell_data["unv:pid"].append(np.array(cell_pid[t], dtype=int))

    # groups -> point_sets / cell_sets
    point_sets = {}
    cell_sets = {}
    # map element label -> (type index in cells list, local index)
    for name, etype, tags in groups:
        if etype == 8:  # nodes
            point_sets[name] = np.array(
                [node_label_to_index[t] for t in tags if t in node_label_to_index],
                dtype=int,
            )
        elif etype == 7:  # elements
            blocks = [np.array([], dtype=int) for _ in cells]
            for t in tags:
                if t not in elem_label_to_ref:
                    continue
                mtype, local = elem_label_to_ref[t]
                bi = type_order.index(mtype)
                blocks[bi] = np.append(blocks[bi], local)
            cell_sets[name] = blocks

    mesh = Mesh(points, cells, cell_data=cell_data if cells else {})
    if point_sets:
        mesh.point_sets = point_sets
    if cell_sets:
        mesh.cell_sets = cell_sets
    return mesh


def _read_2412(body, cell_groups, cell_pid, elem_label_to_ref, node_map):
    k = 0
    while k < len(body):
        rec1 = body[k].split()
        if len(rec1) < 6:
            break
        label, fedesc, pid = int(rec1[0]), int(rec1[1]), int(rec1[2])
        num_nodes = int(rec1[5])
        k += 1
        if fedesc in _beam_descriptors:
            k += 1  # skip beam orientation record
        # gather num_nodes node labels
        node_labels = []
        while len(node_labels) < num_nodes and k < len(body):
            node_labels += [int(v) for v in body[k].split()]
            k += 1
        node_labels = node_labels[:num_nodes]

        if fedesc not in _unv_to_meshio_type:
            raise ReadError(f"UNV: FE descriptor {fedesc} not supported")
        mtype = _unv_to_meshio_type[fedesc]
        unv_conn = [node_map[n] for n in node_labels]
        nd = _ND.get(mtype)
        if nd is None:
            conn = unv_conn
        else:
            conn = [0] * len(nd)
            for i, pos in enumerate(nd):
                conn[pos] = unv_conn[i]

        if mtype not in cell_groups:
            cell_groups[mtype] = []
            cell_pid[mtype] = []
        elem_label_to_ref[label] = (mtype, len(cell_groups[mtype]))
        cell_groups[mtype].append(conn)
        cell_pid[mtype].append(pid)


def _read_2467(body):
    out = []
    k = 0
    while k < len(body):
        rec1 = body[k].split()
        if len(rec1) < 8:
            break
        n_entities = int(rec1[7])
        k += 1
        name = body[k].strip() if k < len(body) else ""
        k += 1
        # entities: 4 ints each, 2 per line
        vals = []
        while len(vals) < 4 * n_entities and k < len(body):
            vals += [int(v) for v in body[k].split()]
            k += 1
        # split by entity type (all-node vs all-element groups are the common case)
        by_type = {}
        for e in range(n_entities):
            etype = vals[4 * e]
            tag = vals[4 * e + 1]
            by_type.setdefault(etype, []).append(tag)
        for etype, tags in by_type.items():
            out.append((name, etype, tags))
    return out


def write(filename, mesh):
    if mesh.points.shape[1] == 2:
        points = np.column_stack([mesh.points, np.zeros(len(mesh.points))])
    else:
        points = mesh.points

    with open_file(filename, "w") as f:
        # 2411 nodes
        f.write("    -1\n  2411\n")
        for k, pt in enumerate(points):
            f.write(f"{k + 1:10d}{1:10d}{1:10d}{11:10d}\n")
            f.write("".join(f"{x:25.16E}" for x in pt) + "\n")
        f.write("    -1\n")

        # 2412 elements
        f.write("    -1\n  2412\n")
        label = 0
        # map (block index, local) -> label for group writing
        elem_labels = []
        for cell_block in mesh.cells:
            t = cell_block.type
            if t not in _meshio_to_unv:
                warn(f"UNV does not support '{t}' cells. Skipping.")
                elem_labels.append(None)
                continue
            descriptor, is_beam = _meshio_to_unv[t]
            nd = _ND.get(t)
            nnodes = cell_block.data.shape[1]
            labels_this = []
            for row in cell_block.data:
                label += 1
                labels_this.append(label)
                if nd is None:
                    unv_row = row
                else:
                    unv_row = [0] * len(nd)
                    for i, pos in enumerate(nd):
                        unv_row[i] = row[pos]
                f.write(
                    f"{label:10d}{descriptor:10d}{1:10d}{1:10d}{11:10d}{nnodes:10d}\n"
                )
                if is_beam:
                    f.write(f"{0:10d}{0:10d}{0:10d}\n")
                # node labels, 8 per line, 1-based
                ints = [int(v) + 1 for v in unv_row]
                for i in range(0, len(ints), 8):
                    f.write("".join(f"{v:10d}" for v in ints[i : i + 8]) + "\n")
            elem_labels.append(labels_this)
        f.write("    -1\n")

        # 2467 groups from point_sets / cell_sets
        psets = getattr(mesh, "point_sets", {}) or {}
        csets = getattr(mesh, "cell_sets", {}) or {}
        if psets or csets:
            f.write("    -1\n  2467\n")
            gid = 0
            for name, ids in psets.items():
                gid += 1
                _write_group(f, gid, name, 8, [int(i) + 1 for i in ids])
            for name, blocks in csets.items():
                gid += 1
                tags = []
                for bi, sel in enumerate(blocks):
                    if elem_labels[bi] is None:
                        continue
                    for local in np.asarray(sel, dtype=int):
                        tags.append(elem_labels[bi][int(local)])
                _write_group(f, gid, name, 7, tags)
            f.write("    -1\n")


def _write_group(f, gid, name, entity_type, tags):
    f.write(f"{gid:10d}{0:10d}{0:10d}{0:10d}{0:10d}{0:10d}{0:10d}{len(tags):10d}\n")
    f.write(f"{name}\n")
    row = []
    for t in tags:
        row += [entity_type, t, 0, 0]
        if len(row) == 8:
            f.write("".join(f"{v:10d}" for v in row) + "\n")
            row = []
    if row:
        f.write("".join(f"{v:10d}" for v in row) + "\n")
