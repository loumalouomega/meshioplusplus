"""
I/O for Exodus II.

See
<https://src.fedoraproject.org/repo/pkgs/exodusii/922137.pdf/a45d67f4a1a8762bcf66af2ec6eb35f9/922137.pdf>,
in particular Appendix A (page 171, Implementation of EXODUS II with netCDF).
"""

import datetime
import re

import numpy as np

from ..__about__ import __version__
from .._exceptions import ReadError
from .._mesh import Mesh, topological_dimension

exodus_to_meshio_type = {
    "SPHERE": "vertex",
    # curves
    "BEAM": "line",
    "BEAM2": "line",
    "BEAM3": "line3",
    "BAR2": "line",
    # surfaces
    "SHELL": "quad",
    "SHELL4": "quad",
    "SHELL8": "quad8",
    "SHELL9": "quad9",
    "QUAD": "quad",
    "QUAD4": "quad",
    "QUAD5": "quad5",
    "QUAD8": "quad8",
    "QUAD9": "quad9",
    #
    "TRI": "triangle",
    "TRIANGLE": "triangle",
    "TRI3": "triangle",
    "TRI6": "triangle6",
    "TRI7": "triangle7",
    # 'TRISHELL': 'triangle',
    # 'TRISHELL3': 'triangle',
    # 'TRISHELL6': 'triangle6',
    # 'TRISHELL7': 'triangle',
    #
    # volumes
    "HEX": "hexahedron",
    "HEXAHEDRON": "hexahedron",
    "HEX8": "hexahedron",
    "HEX9": "hexahedron9",
    "HEX20": "hexahedron20",
    "HEX27": "hexahedron27",
    #
    "TETRA": "tetra",
    "TETRA4": "tetra4",
    "TET4": "tetra4",
    "TETRA8": "tetra8",
    "TETRA10": "tetra10",
    "TETRA14": "tetra14",
    #
    "PYRAMID": "pyramid",
    "WEDGE": "wedge",
}
meshio_to_exodus_type = {v: k for k, v in exodus_to_meshio_type.items()}


# Exodus numbers an element's sides in its own order, which is NOT
# ``detail/cell_faces.hpp``'s. This is the Python twin of ``exo_face_index`` in
# ``src/cpp/src/formats/exodus.cpp`` — KEEP THE TWO IN SYNC. Each entry is the
# meshio++ facet index whose corner-node set equals that Exodus side's; the
# derivations are spelled out on the C++ side.
_EXODUS_SIDE_TO_FACET = {
    "tetra": (0, 1, 2, 3),
    "hexahedron": (2, 1, 3, 0, 4, 5),
    "wedge": (2, 3, 4, 0, 1),
    "pyramid": (1, 2, 3, 4, 0),
    # 2-D elements: Exodus walks the edges 1-2, 2-3, ... in the same order
    # ``detail/cell_edges.hpp`` does, so the mapping is the identity.
    "triangle": (0, 1, 2),
    "quad": (0, 1, 2, 3),
}

# Higher-order variants share their linear base's facet ordering.
_EXODUS_FACET_BASE = {
    "tetra4": "tetra",
    "tetra10": "tetra",
    "tetra14": "tetra",
    "hexahedron20": "hexahedron",
    "hexahedron27": "hexahedron",
    "wedge15": "wedge",
    "wedge18": "wedge",
    "pyramid13": "pyramid",
    "pyramid14": "pyramid",
    "triangle6": "triangle",
    "triangle7": "triangle",
    "quad8": "quad",
    "quad9": "quad",
}


def exodus_face_index(cell_type, exodus_side):
    """Map a 1-based Exodus side number to a meshio++ local facet index.

    Returns ``-1`` when the pair has no mapping, which the caller skips rather
    than storing a facet pointing at the wrong face.
    """
    base = _EXODUS_FACET_BASE.get(cell_type, cell_type)
    table = _EXODUS_SIDE_TO_FACET.get(base)
    if table is None or exodus_side < 1 or exodus_side > len(table):
        return -1
    return table[exodus_side - 1]


def _group_name(names, ids, index, prefix):
    """The name for set/block ``index``, or a stable synthetic one.

    SEACAS writes blank names for unnamed groups, but a group with no name is
    still a group — dropping it would lose the only handle a consumer has on it.
    """
    if index < len(names) and names[index]:
        return names[index]
    gid = ids[index] if index < len(ids) else index + 1
    return f"{prefix} {int(gid)}"


def _resolve_time_step(time_step, num_steps):
    """Resolve ``time_step`` against an actual step count.

    ``0`` is the first step (the historical behaviour); negative counts from the
    end. Out of range is an error naming the available count, never a silent
    clamp — quietly returning step 0 when step 7 was asked for is exactly the
    failure this option exists to remove. Twin of
    ``ReadOptions::ResolveTimeStep``.
    """
    if num_steps == 0:
        if time_step in (0, -1):
            return 0
        raise ReadError(
            f"meshio++: time step {time_step} requested, "
            "but this file carries no time steps"
        )
    resolved = num_steps + time_step if time_step < 0 else time_step
    if resolved < 0 or resolved >= num_steps:
        raise ReadError(
            f"meshio++: time step {time_step} is out of range: this file has "
            f"{num_steps} step{'' if num_steps == 1 else 's'}"
        )
    return resolved


def read(filename, time_step=0):  # noqa: C901
    import netCDF4

    with netCDF4.Dataset(filename) as nc:
        # assert nc.version == np.float32(5.1)
        # assert nc.api_version == np.float32(5.1)
        # assert nc.floating_point_word_size == 8

        # assert b''.join(nc.variables['coor_names'][0]) == b'X'
        # assert b''.join(nc.variables['coor_names'][1]) == b'Y'
        # assert b''.join(nc.variables['coor_names'][2]) == b'Z'

        points = np.zeros((len(nc.dimensions["num_nodes"]), 3))
        point_data_names = []
        cell_data_names = []
        pd = {}
        cd = {}
        cells = []
        ns_names = []
        eb_names = []
        ss_names = []
        eb_ids = []
        ns_ids = []
        ss_ids = []
        ns = []
        side_elems = {}
        side_sides = {}
        info = []

        # Resolve the requested step up front, so an out-of-range request fails
        # before any heavy array is sliced rather than midway through.
        time_values = (
            np.asarray(nc.variables["time_whole"][:]).ravel()
            if "time_whole" in nc.variables
            else np.zeros(0)
        )
        step = _resolve_time_step(time_step, len(time_values))

        def _slice_step(value, key):
            # `time_whole` and a `vals_*` array can legitimately disagree on
            # length (a writer that died mid-step leaves a short array), so the
            # resolved step is re-checked against the array actually sliced.
            if step >= len(value):
                raise ReadError(
                    f"Exodus: time step {step} is out of range for '{key}', "
                    f"which has {len(value)} step{'' if len(value) == 1 else 's'}"
                )
            return value[step]

        for key, value in nc.variables.items():
            if key == "info_records":
                value.set_auto_mask(False)
                for c in value[:]:
                    try:
                        info += [b"".join(c).decode("UTF-8")]
                    except UnicodeDecodeError:
                        # https://github.com/nschloe/meshio/issues/983
                        pass
            elif key == "qa_records":
                value.set_auto_mask(False)
                for val in value:
                    info += [b"".join(c).decode("UTF-8") for c in val[:]]
            elif key[:7] == "connect":
                meshio_type = exodus_to_meshio_type[value.elem_type.upper()]
                cells.append((meshio_type, value[:] - 1))
            elif key == "coord":
                points = nc.variables["coord"][:].T
            elif key == "coordx":
                points[:, 0] = value[:]
            elif key == "coordy":
                points[:, 1] = value[:]
            elif key == "coordz":
                points[:, 2] = value[:]
            elif key == "name_nod_var":
                value.set_auto_mask(False)
                point_data_names = [b"".join(c).decode("UTF-8") for c in value[:]]
            elif key[:12] == "vals_nod_var":
                idx = 0 if len(key) == 12 else int(key[12:]) - 1
                value.set_auto_mask(False)
                pd[idx] = _slice_step(value, key)
            elif key == "name_elem_var":
                value.set_auto_mask(False)
                cell_data_names = [b"".join(c).decode("UTF-8") for c in value[:]]
            elif key[:13] == "vals_elem_var":
                # eb: element block
                m = re.match("vals_elem_var(\\d+)?(?:eb(\\d+))?", key)
                idx = 0 if m.group(1) is None else int(m.group(1)) - 1
                block = 0 if m.group(2) is None else int(m.group(2)) - 1

                value.set_auto_mask(False)
                if idx not in cd:
                    cd[idx] = {}
                cd[idx][block] = _slice_step(value, key)
            elif key == "ns_names":
                value.set_auto_mask(False)
                ns_names = [b"".join(c).decode("UTF-8") for c in value[:]]
            elif key == "eb_names":
                value.set_auto_mask(False)
                eb_names = [b"".join(c).decode("UTF-8") for c in value[:]]
            elif key == "ss_names":
                value.set_auto_mask(False)
                ss_names = [b"".join(c).decode("UTF-8") for c in value[:]]
            elif key == "eb_prop1":
                eb_ids = np.asarray(value[:]).ravel().tolist()
            elif key == "ns_prop1":
                ns_ids = np.asarray(value[:]).ravel().tolist()
            elif key == "ss_prop1":
                ss_ids = np.asarray(value[:]).ravel().tolist()
            elif key.startswith("node_ns"):  # Expected keys: node_ns1, node_ns2
                ns.append(np.asarray(value[:]) - 1)  # Exodus is 1-based
            elif key.startswith("elem_ss"):  # side set: owning element ids
                side_elems[int(key[7:] or 1)] = np.asarray(value[:])
            elif key.startswith("side_ss"):  # side set: Exodus side numbers
                side_sides[int(key[7:] or 1)] = np.asarray(value[:])

        # merge element block data; can't handle blocks yet
        for k, value in cd.items():
            cd[k] = np.concatenate(list(value.values()))

        # Check if there are any <name>R, <name>Z tuples or <name>X, <name>Y, <name>Z
        # triplets in the point data. If yes, they belong together.
        single, double, triple = categorize(point_data_names)

        point_data = {}
        for name, idx in single:
            point_data[name] = pd[idx]
        for name, idx0, idx1 in double:
            point_data[name] = np.column_stack([pd[idx0], pd[idx1]])
        for name, idx0, idx1, idx2 in triple:
            point_data[name] = np.column_stack([pd[idx0], pd[idx1], pd[idx2]])

        cell_data = {}
        k = 0
        for _, cell in cells:
            n = len(cell)
            for name, data in zip(cell_data_names, cd.values()):
                if name not in cell_data:
                    cell_data[name] = []
                cell_data[name].append(data[k : k + n])
            k += n

    mesh = Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        info=info,
    )
    # Attached after construction rather than passed as `regions=`, so
    # `block_bases` sees real CellBlocks and the global block-major numbering
    # comes from its one owner instead of being re-derived from raw tuples.
    # `point_sets` is a write-through view over these, so it needs no separate
    # argument -- and passing both would have the set setter wipe the regions.
    mesh.regions = _build_regions(
        mesh.cells,
        eb_names,
        eb_ids,
        ns_names,
        ns_ids,
        ns,
        ss_names,
        ss_ids,
        side_elems,
        side_sides,
        len(mesh.points),
    )
    # Side channel like `.regions`/`.info`: `read_metadata`'s generic
    # (no-native-core) fallback builds its summary from an already-read Mesh
    # and has nowhere else to learn the file's full time-series -- the
    # returned mesh itself only ever holds the one requested step.
    mesh.time_values = [float(v) for v in time_values]
    return mesh


def _build_regions(  # noqa: C901
    cells,
    eb_names,
    eb_ids,
    ns_names,
    ns_ids,
    ns,
    ss_names,
    ss_ids,
    side_elems,
    side_sides,
    num_points,
):
    """Turn Exodus element blocks, node sets and side sets into Regions.

    Twin of ``exo_add_regions`` in ``src/cpp/src/formats/exodus.cpp``. Cell/Side
    entries are global block-major cell indices, taken from ``block_bases`` —
    the Python owner of ``detail/cell_index.hpp`` — rather than re-derived.
    """
    from .._regions import Region, block_bases

    bases = block_bases(cells)
    total_cells = int(bases[-1]) if len(bases) else 0
    regions = []

    # Element blocks -> Cell regions. The name and tag come from per-block
    # arrays, which is what keeps two blocks of the SAME element type
    # distinguishable — exactly what a materials assignment depends on.
    for b, cb in enumerate(cells):
        tag = eb_ids[b] if b < len(eb_ids) else -1
        regions.append(
            Region(
                _group_name(eb_names, eb_ids, b, "Block"),
                "cell",
                np.arange(bases[b], bases[b + 1], dtype=np.int64),
                dim=topological_dimension.get(cb.type, -1),
                tag=tag,
            )
        )

    # Node sets -> Point regions.
    for k, nodes in enumerate(ns):
        nodes = np.asarray(nodes, dtype=np.int64).ravel()
        nodes = nodes[(nodes >= 0) & (nodes < num_points)]
        regions.append(
            Region(
                _group_name(ns_names, ns_ids, k, "Nodeset"),
                "point",
                nodes,
                tag=ns_ids[k] if k < len(ns_ids) else -1,
            )
        )

    # Side sets -> Side regions: (global cell, local facet) pairs. The facet goes
    # through exodus_face_index; an unmappable pair is skipped rather than stored
    # pointing at the wrong face.
    for k, set_id in enumerate(sorted(side_elems)):
        elems = np.asarray(side_elems[set_id], dtype=np.int64).ravel() - 1
        sides = np.asarray(side_sides.get(set_id, []), dtype=np.int64).ravel()
        pairs = []
        for g, side in zip(elems, sides):
            if g < 0 or g >= total_cells:
                continue
            b = int(np.searchsorted(bases, g, side="right")) - 1
            if b < 0 or b >= len(cells):
                continue
            facet = exodus_face_index(cells[b].type, int(side))
            if facet < 0:
                continue
            pairs.append((int(g), facet))
        regions.append(
            Region(
                _group_name(ss_names, ss_ids, k, "Sideset"),
                "side",
                np.asarray(pairs, dtype=np.int64).reshape(-1, 2),
                tag=ss_ids[k] if k < len(ss_ids) else -1,
            )
        )

    return regions


def categorize(names):
    # Check if there are any <name>R, <name>Z tuples or <name>X, <name>Y, <name>Z
    # triplets in the point data. If yes, they belong together.
    single = []
    double = []
    triple = []
    is_accounted_for = [False] * len(names)
    k = 0
    while True:
        if k == len(names):
            break
        if is_accounted_for[k]:
            k += 1
            continue
        name = names[k]
        if name[-1] == "X":
            ix = k
            try:
                iy = names.index(name[:-1] + "Y")
            except ValueError:
                iy = None
            try:
                iz = names.index(name[:-1] + "Z")
            except ValueError:
                iz = None
            if iy and iz:
                triple.append((name[:-1], ix, iy, iz))
                is_accounted_for[ix] = True
                is_accounted_for[iy] = True
                is_accounted_for[iz] = True
            else:
                single.append((name, ix))
                is_accounted_for[ix] = True
        elif name[-2:] == "_R":
            ir = k
            try:
                iz = names.index(name[:-2] + "_Z")
            except ValueError:
                iz = None
            if iz:
                double.append((name[:-2], ir, iz))
                is_accounted_for[ir] = True
                is_accounted_for[iz] = True
            else:
                single.append((name, ir))
                is_accounted_for[ir] = True
        else:
            single.append((name, k))
            is_accounted_for[k] = True

        k += 1

    if not all(is_accounted_for):
        raise ReadError()
    return single, double, triple


numpy_to_exodus_dtype = {
    "float32": "f4",
    "float64": "f8",
    "int8": "i1",
    "int16": "i2",
    "int32": "i4",
    "int64": "i8",
    "uint8": "u1",
    "uint16": "u2",
    "uint32": "u4",
    "uint64": "u8",
}


def write(filename, mesh):
    import netCDF4

    with netCDF4.Dataset(filename, "w") as rootgrp:
        # set global data
        now = datetime.datetime.now().isoformat()
        rootgrp.title = f"Created by meshio++ v{__version__}, {now}"
        rootgrp.version = np.float32(5.1)
        rootgrp.api_version = np.float32(5.1)
        rootgrp.floating_point_word_size = 8

        # set dimensions
        total_num_elems = sum(c.data.shape[0] for c in mesh.cells)
        rootgrp.createDimension("num_nodes", len(mesh.points))
        rootgrp.createDimension("num_dim", mesh.points.shape[1])
        rootgrp.createDimension("num_elem", total_num_elems)
        rootgrp.createDimension("num_el_blk", len(mesh.cells))
        rootgrp.createDimension("num_node_sets", len(mesh.point_sets))
        rootgrp.createDimension("len_string", 33)
        rootgrp.createDimension("len_line", 81)
        rootgrp.createDimension("four", 4)
        rootgrp.createDimension("time_step", None)

        # dummy time step
        data = rootgrp.createVariable("time_whole", "f4", ("time_step",))
        data[:] = 0.0

        # points
        coor_names = rootgrp.createVariable(
            "coor_names", "S1", ("num_dim", "len_string")
        )
        coor_names.set_auto_mask(False)
        coor_names[0, 0] = b"X"
        coor_names[1, 0] = b"Y"
        if mesh.points.shape[1] == 3:
            coor_names[2, 0] = b"Z"
        data = rootgrp.createVariable(
            "coord",
            numpy_to_exodus_dtype[mesh.points.dtype.name],
            ("num_dim", "num_nodes"),
        )
        data[:] = mesh.points.T

        # cells
        # ParaView needs eb_prop1 -- some ID. The values don't seem to matter as
        # long as they are different for the for different blocks.
        data = rootgrp.createVariable("eb_prop1", "i4", "num_el_blk")
        for k in range(len(mesh.cells)):
            data[k] = k
        for k, cell_block in enumerate(mesh.cells):
            dim1 = f"num_el_in_blk{k + 1}"
            dim2 = f"num_nod_per_el{k + 1}"
            rootgrp.createDimension(dim1, cell_block.data.shape[0])
            rootgrp.createDimension(dim2, cell_block.data.shape[1])
            dtype = numpy_to_exodus_dtype[cell_block.data.dtype.name]
            data = rootgrp.createVariable(f"connect{k + 1}", dtype, (dim1, dim2))
            data.elem_type = meshio_to_exodus_type[cell_block.type]
            # Exodus is 1-based
            data[:] = cell_block.data + 1

        # point data
        # The variable `name_nod_var` holds the names and indices of the node variables, the
        # variables `vals_nod_var{1,2,...}` hold the actual data.
        num_nod_var = len(mesh.point_data)
        if num_nod_var > 0:
            rootgrp.createDimension("num_nod_var", num_nod_var)
            # set names
            point_data_names = rootgrp.createVariable(
                "name_nod_var", "S1", ("num_nod_var", "len_string")
            )
            point_data_names.set_auto_mask(False)
            for k, name in enumerate(mesh.point_data.keys()):
                for i, letter in enumerate(name):
                    point_data_names[k, i] = letter.encode()

            # Set data. ParaView might have some problems here, see
            # <https://gitlab.kitware.com/paraview/paraview/-/issues/18403>.
            for k, (name, data) in enumerate(mesh.point_data.items()):
                for i, s in enumerate(data.shape):
                    rootgrp.createDimension(f"dim_nod_var{k}{i}", s)
                dims = ["time_step"] + [
                    f"dim_nod_var{k}{i}" for i in range(len(data.shape))
                ]
                node_data = rootgrp.createVariable(
                    f"vals_nod_var{k + 1}",
                    numpy_to_exodus_dtype[data.dtype.name],
                    tuple(dims),
                    fill_value=False,
                )
                node_data[0] = data

        # node sets
        num_point_sets = len(mesh.point_sets)
        if num_point_sets > 0:
            data = rootgrp.createVariable("ns_prop1", "i4", "num_node_sets")
            data_names = rootgrp.createVariable(
                "ns_names", "S1", ("num_node_sets", "len_string")
            )
            for k, name in enumerate(mesh.point_sets.keys()):
                data[k] = k
                for i, letter in enumerate(name):
                    data_names[k, i] = letter.encode()
            for k, (key, values) in enumerate(mesh.point_sets.items()):
                dim1 = f"num_nod_ns{k + 1}"
                rootgrp.createDimension(dim1, values.shape[0])
                dtype = numpy_to_exodus_dtype[values.dtype.name]
                data = rootgrp.createVariable(f"node_ns{k + 1}", dtype, (dim1,))
                # Exodus is 1-based
                data[:] = values + 1
