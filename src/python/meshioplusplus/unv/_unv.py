"""
I/O for the I-DEAS Universal file format (``.unv``), following the datasets used
by FEconv <https://github.com/victorsndvg/FEconv>: 2411/781 (nodes), 2412
(elements, with the FE-descriptor -> element-type map and the Salome/Code-Aster
mid-node "sandwich" ordering for parabolic elements), 2467/2477/2452/2435/2432/
2430 (permanent groups, mapped to point/cell sets) and the field/results
datasets 2414 (default) and legacy 55/57 (Code-Aster mode), mapped to
``point_data`` (data at nodes) / ``cell_data`` (data on elements).

A UNV file is a sequence of datasets, each delimited by a line whose content is
``-1``, followed by a dataset-id line and the dataset records.
"""

import numpy as np

from .. import _provenance
from .._common import warn
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
    "wedge15": [0, 6, 1, 7, 2, 8, 9, 10, 11, 3, 12, 4, 13, 5, 14],
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
    113: "wedge15",
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
    "wedge15": (113, False),
    "hexahedron": (115, False),
    "hexahedron20": (116, False),
}
_beam_descriptors = {11, 21, 22, 24}

# Group datasets that FEconv reads (all share the 2467 record layout); 2467 is
# used for writing.
_group_datasets = {2467, 2477, 2452, 2435, 2432, 2430}

# Field datasets: 2414 (modern, node/element/point location) + legacy 55 (data
# at nodes) / 56 (data at nodes on elements) / 57 (data at elements).  On write
# the default is 2414; Code-Aster mode emits 55 (node data) / 57 (element data).
_field_datasets = {2414, 55, 56, 57}

# UNV data-characteristic code -> number of components.  NDV (record-9 field 6)
# is authoritative when present; this maps the write direction.
_ncomp_to_char = {1: 1, 3: 2, 6: 4, 9: 5}
# Data-type codes considered real (2 = single, 4 = double precision float).
_real_data_types = {2, 4}


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
    node_fields = []  # (key, ncomp, {node_label: [values]})
    elem_fields = []  # (key, ncomp, {elem_label: [values]})
    used_keys = set()

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
        elif ds_id in _group_datasets:
            groups.extend(_read_2467(body))
        elif ds_id in _field_datasets:
            parsed = _read_field(ds_id, body)
            if parsed is None:
                continue
            location, ncomp, name, values = parsed
            key = _unique_key(name, used_keys)
            if location == 1:  # data at nodes
                node_fields.append((key, ncomp, values))
            elif location == 2:  # data on elements
                elem_fields.append((key, ncomp, values))
            # location 3 (nodes-on-elements) is warned+skipped in _read_field
        # other datasets are ignored on read

    points = np.array(points, dtype=float) if points else np.empty((0, 3))

    cells = []
    cell_data = {"unv:pid": []}
    type_order = list(cell_groups.keys())
    for t in type_order:
        cells.append(CellBlock(t, np.array(cell_groups[t], dtype=int)))
        cell_data["unv:pid"].append(np.array(cell_pid[t], dtype=int))

    # node fields -> point_data (indexed by node label -> point index)
    point_data = {}
    for key, ncomp, values in node_fields:
        arr = np.zeros((len(points), ncomp), dtype=float)
        for label, vals in values.items():
            if label in node_label_to_index:
                arr[node_label_to_index[label]] = vals[:ncomp]
        point_data[key] = arr[:, 0] if ncomp == 1 else arr

    # element fields -> cell_data (one array per block, aligned by element label)
    for key, ncomp, values in elem_fields:
        blocks = [np.zeros((len(cb.data), ncomp), dtype=float) for cb in cells]
        for label, vals in values.items():
            ref = elem_label_to_ref.get(label)
            if ref is None:
                continue
            mtype, local = ref
            bi = type_order.index(mtype)
            blocks[bi][local] = vals[:ncomp]
        cell_data[key] = [b[:, 0] if ncomp == 1 else b for b in blocks]

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

    if not cells:
        cell_data = {}
    mesh = Mesh(
        points,
        cells,
        point_data=point_data if point_data else {},
        cell_data=cell_data,
    )
    if point_sets:
        mesh.point_sets = point_sets
    if cell_sets:
        mesh.cell_sets = cell_sets
    return mesh


def _unique_key(name, used_keys):
    """Return a collision-free data key derived from a UNV field name."""
    base = name.strip() or "unv:field"
    key = base
    n = 1
    while key in used_keys:
        n += 1
        key = f"{base}_{n}"
    used_keys.add(key)
    return key


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
            warn(f"UNV: FE descriptor {fedesc} not supported; skipping element.")
            continue
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


def _read_field(ds_id, body):
    """Parse a field dataset (2414 or legacy 55/56/57).

    Returns ``(location, ncomp, name, {entity_label: [values]})`` or ``None``
    when the dataset is unsupported (complex data, or nodes-on-elements data,
    which is warned about and skipped).  ``location`` is 1 (nodes) or 2
    (elements).
    """
    lines = [ln for ln in body]
    if ds_id == 2414:
        # record1 label, record2 name, record3 location, records 4-8 id lines,
        # record9 six ints (model, analysis, char, spec, data_type, ndv),
        # records 10-13 header, then per-entity (label line + value line(s)).
        if len(lines) < 13:
            return None
        name = lines[1].strip()
        location = int(lines[2].split()[0]) if lines[2].split() else 1
        rec9 = lines[8].split()
        if len(rec9) < 6:
            return None
        data_type = int(rec9[4])
        ndv = int(rec9[5])
        data_start = 13
    else:
        # Legacy 55/56/57: records 1-5 id lines, record6 six ints
        # (model, analysis, char, spec, data_type, ndv), records 7-10 header,
        # then per-entity (label line + value line(s)).  55 = nodes,
        # 56 = nodes-on-elements, 57 = elements.
        if len(lines) < 10:
            return None
        name = lines[0].strip()
        rec6 = lines[5].split()
        if len(rec6) < 6:
            return None
        data_type = int(rec6[4])
        ndv = int(rec6[5])
        location = {55: 1, 56: 3, 57: 2}.get(ds_id, 1)
        data_start = 10

    if data_type not in _real_data_types:
        warn(f"UNV: skipping complex field '{name}' (dataset {ds_id}).")
        return None
    if location == 3:
        warn(
            f"UNV: field '{name}' is at nodes-on-elements (dataset {ds_id}); "
            "averaging to element barycenters is not implemented, skipping."
        )
        return None
    if ndv <= 0:
        return None

    values = {}
    k = data_start
    while k < len(lines):
        rec = lines[k].split()
        if not rec:
            k += 1
            continue
        label = int(rec[0])
        k += 1
        vals = []
        while len(vals) < ndv and k < len(lines):
            for v in lines[k].replace("D", "E").replace("d", "e").split():
                vals.append(float(v))
            k += 1
        values[label] = vals[:ndv]
    return location, ndv, name, values


def write(filename, mesh, code_aster=False, node_dataset=2411):
    if mesh.points.shape[1] == 2:
        points = np.column_stack([mesh.points, np.zeros(len(mesh.points))])
    else:
        points = mesh.points

    if node_dataset not in (2411, 781):
        node_dataset = 2411

    with open_file(filename, "w") as f:
        # 2411 / 781 nodes
        f.write(f"    -1\n{node_dataset:6d}\n")
        for k, pt in enumerate(points):
            f.write(f"{k + 1:10d}{1:10d}{1:10d}{11:10d}\n")
            f.write("".join(f"{x:25.16E}" for x in pt) + "\n")
        f.write("    -1\n")

        # 2412 elements
        f.write("    -1\n  2412\n")
        label = 0
        pid_blocks = (getattr(mesh, "cell_data", None) or {}).get("unv:pid")
        # map (block index, local) -> label for group writing
        elem_labels = []
        for bi, cell_block in enumerate(mesh.cells):
            t = cell_block.type
            if t not in _meshio_to_unv:
                warn(f"UNV does not support '{t}' cells. Skipping.")
                _provenance.note(
                    "cells-dropped", f"cell block(s) of type {t} have no UNV equivalent"
                )
                elem_labels.append(None)
                continue
            descriptor, is_beam = _meshio_to_unv[t]
            nd = _ND.get(t)
            nnodes = cell_block.data.shape[1]
            pids = (
                np.asarray(pid_blocks[bi], dtype=int)
                if pid_blocks is not None and bi < len(pid_blocks)
                else None
            )
            labels_this = []
            for li, row in enumerate(cell_block.data):
                label += 1
                labels_this.append(label)
                pid = int(pids[li]) if pids is not None and li < len(pids) else 1
                if nd is None:
                    unv_row = row
                else:
                    unv_row = [0] * len(nd)
                    for i, pos in enumerate(nd):
                        unv_row[i] = row[pos]
                f.write(
                    f"{label:10d}{descriptor:10d}{pid:10d}{pid:10d}{11:10d}{nnodes:10d}\n"
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

        # field datasets from point_data / cell_data
        _write_fields(f, mesh, elem_labels, code_aster)


def _field_char(ncomp):
    """Map a component count to a UNV (data_characteristic, ndv) pair."""
    return _ncomp_to_char.get(ncomp, 0), ncomp


def _as_2d(arr):
    arr = np.asarray(arr, dtype=float)
    return arr.reshape(len(arr), 1) if arr.ndim == 1 else arr


def _write_fields(f, mesh, elem_labels, code_aster):
    field_id = 0
    # data at nodes -> location 1 (2414) or dataset 55 (Code Aster)
    for name, arr in (getattr(mesh, "point_data", None) or {}).items():
        arr = _as_2d(arr)
        ncomp = arr.shape[1]
        field_id += 1
        node_labels = range(1, len(arr) + 1)
        if code_aster:
            _write_field_55_57(f, 55, name, ncomp, node_labels, arr)
        else:
            _write_field_2414(f, field_id, name, 1, ncomp, node_labels, arr)

    # data on elements -> location 2 (2414) or dataset 57 (Code Aster)
    for name, blocks in (getattr(mesh, "cell_data", None) or {}).items():
        if name == "unv:pid":  # element property id, carried by dataset 2412
            continue
        # gather (element_label, values) across blocks
        labels = []
        rows = []
        ncomp = None
        for bi, blk in enumerate(blocks):
            if elem_labels[bi] is None:
                continue
            blk = _as_2d(blk)
            if ncomp is None:
                ncomp = blk.shape[1]
            for local in range(len(blk)):
                labels.append(elem_labels[bi][local])
                rows.append(blk[local])
        if ncomp is None or not rows:
            continue
        field_id += 1
        data = np.array(rows, dtype=float)
        if code_aster:
            _write_field_55_57(f, 57, name, ncomp, labels, data)
        else:
            _write_field_2414(f, field_id, name, 2, ncomp, labels, data)


def _write_values(f, labels, data):
    for label, vals in zip(labels, data):
        f.write(f"{int(label):10d}\n")
        f.write("".join(f"{v:13.5E}" for v in vals) + "\n")


def _write_field_2414(f, field_id, name, location, ncomp, labels, data):
    char, ndv = _field_char(ncomp)
    f.write("    -1\n  2414\n")
    f.write(f"{field_id:10d}\n")  # record 1: analysis dataset label
    f.write(f"{name}\n")  # record 2: analysis dataset name
    f.write(f"{location:10d}\n")  # record 3: 1=nodes, 2=elements
    for _ in range(5):  # records 4-8: ID lines
        f.write("meshioplusplus\n")
    # record 9: model, analysis, data_char, spec, data_type(4=double), ndv
    f.write(f"{1:10d}{0:10d}{char:10d}{0:10d}{4:10d}{ndv:10d}\n")
    f.write(f"{0:10d}" * 8 + "\n")  # record 10
    f.write(f"{0:10d}" * 2 + "\n")  # record 11
    f.write("".join(f"{0.0:13.5E}" for _ in range(6)) + "\n")  # record 12
    f.write("".join(f"{0.0:13.5E}" for _ in range(6)) + "\n")  # record 13
    _write_values(f, labels, data)
    f.write("    -1\n")


def _write_field_55_57(f, ds_id, name, ncomp, labels, data):
    char, ndv = _field_char(ncomp)
    f.write(f"    -1\n{ds_id:6d}\n")
    for _ in range(5):  # records 1-5: ID lines
        f.write(f"{name}\n")
    # record 6: model, analysis, data_char, spec, data_type(4=double), ndv
    f.write(f"{1:10d}{0:10d}{char:10d}{0:10d}{4:10d}{ndv:10d}\n")
    f.write(f"{0:10d}" * 8 + "\n")  # record 7
    f.write(f"{0:10d}" * 2 + "\n")  # record 8
    f.write("".join(f"{0.0:13.5E}" for _ in range(6)) + "\n")  # record 9
    f.write("".join(f"{0.0:13.5E}" for _ in range(6)) + "\n")  # record 10
    _write_values(f, labels, data)
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
