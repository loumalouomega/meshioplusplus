"""
I/O for Abaqus inp files.
"""

import pathlib
from itertools import count

import numpy as np

from .. import _provenance
from .._common import num_nodes_per_cell
from .._exceptions import ReadError
from .._files import open_file
from .._mesh import CellBlock, Mesh

abaqus_to_meshio_type = {
    # trusses
    "T2D2": "line",
    "T2D2H": "line",
    "T2D3": "line3",
    "T2D3H": "line3",
    "T3D2": "line",
    "T3D2H": "line",
    "T3D3": "line3",
    "T3D3H": "line3",
    # beams
    "B21": "line",
    "B21H": "line",
    "B22": "line3",
    "B22H": "line3",
    "B31": "line",
    "B31H": "line",
    "B32": "line3",
    "B32H": "line3",
    "B33": "line3",
    "B33H": "line3",
    # surfaces
    "CPS4": "quad",
    "CPS4R": "quad",
    "S4": "quad",
    "S4R": "quad",
    "S4RS": "quad",
    "S4RSW": "quad",
    "S4R5": "quad",
    "S8R": "quad8",
    "S8R5": "quad8",
    "S9R5": "quad9",
    # "QUAD": "quad",
    # "QUAD4": "quad",
    # "QUAD5": "quad5",
    # "QUAD8": "quad8",
    # "QUAD9": "quad9",
    #
    "CPS3": "triangle",
    "STRI3": "triangle",
    "S3": "triangle",
    "S3R": "triangle",
    "S3RS": "triangle",
    "R3D3": "triangle",
    # "TRI7": "triangle7",
    # 'TRISHELL': 'triangle',
    # 'TRISHELL3': 'triangle',
    # 'TRISHELL7': 'triangle',
    #
    "STRI65": "triangle6",
    # 'TRISHELL6': 'triangle6',
    # volumes
    "C3D8": "hexahedron",
    "C3D8H": "hexahedron",
    "C3D8I": "hexahedron",
    "C3D8IH": "hexahedron",
    "C3D8R": "hexahedron",
    "C3D8RH": "hexahedron",
    # "HEX9": "hexahedron9",
    "C3D20": "hexahedron20",
    "C3D20H": "hexahedron20",
    "C3D20R": "hexahedron20",
    "C3D20RH": "hexahedron20",
    # "HEX27": "hexahedron27",
    #
    "C3D4": "tetra",
    "C3D4H": "tetra4",
    # "TETRA8": "tetra8",
    "C3D10": "tetra10",
    "C3D10H": "tetra10",
    "C3D10I": "tetra10",
    "C3D10M": "tetra10",
    "C3D10MH": "tetra10",
    # "TETRA14": "tetra14",
    #
    # "PYRAMID": "pyramid",
    "C3D6": "wedge",
    "C3D15": "wedge15",
    #
    # 4-node bilinear displacement and pore pressure
    "CAX4P": "quad",
    # 6-node quadratic
    "CPE6": "triangle6",
}
meshio_to_abaqus_type = {v: k for k, v in abaqus_to_meshio_type.items()}


# Abaqus face identifier (S1..S6) -> the meshio++ local facet index of
# `detail/cell_faces.hpp` (3-D) / `cell_edges.hpp` (2-D). The two numberings
# genuinely differ -- Abaqus C3D8 S1 is the 1-2-3-4 face, i.e. local nodes
# {0,1,2,3}, which is meshio++'s face 4 -- so this is the Python twin of the
# C++ table and the two must stay in step. Shell SPOS/SNEG name a side rather
# than a facet and map to 0/1 (see doc/regions.md).
_ABAQUS_FACE_ORDER = {
    "tetra": [3, 0, 1, 2],
    "tetra10": [3, 0, 1, 2],
    "hexahedron": [4, 5, 2, 1, 3, 0],
    "hexahedron20": [4, 5, 2, 1, 3, 0],
    "wedge": [0, 1, 2, 3, 4],
    "wedge15": [0, 1, 2, 3, 4],
    "triangle": [0, 1, 2],
    "triangle6": [0, 1, 2],
    "quad": [0, 1, 2, 3],
    "quad8": [0, 1, 2, 3],
    "quad9": [0, 1, 2, 3],
}


def _face_index(cell_type, face):
    """Abaqus face identifier -> local facet index, or None."""
    face = face.strip().upper()
    if face == "SPOS":
        return 0
    if face == "SNEG":
        return 1
    if not face.startswith("S") or not face[1:].isdigit():
        return None
    order = _ABAQUS_FACE_ORDER.get(cell_type)
    n = int(face[1:]) - 1
    if order is None or n < 0 or n >= len(order):
        return None
    return order[n]


def _face_name(cell_type, facet):
    """The inverse of :func:`_face_index`, for the writer."""
    order = _ABAQUS_FACE_ORDER.get(cell_type)
    if order is None:
        return None
    for n, f in enumerate(order):
        if f == facet:
            return f"S{n + 1}"
    return None


def read(filename):
    """Reads a Abaqus inp file."""
    with open_file(filename, "r") as f:
        out = read_buffer(f)
    return out


def read_buffer(f):
    # Initialize the optional data fields
    points = []
    cells = []
    cell_ids = []
    point_sets = {}
    cell_sets = {}
    cell_sets_element = {}  # Handle cell sets defined in ELEMENT
    cell_sets_element_order = []  # Order of keys is not preserved in Python 3.5
    surfaces = []  # (name, [(element id | elset name, face identifier)])
    field_data = {}
    cell_data = {}
    point_data = {}
    point_ids = None

    line = f.readline()
    while True:
        if not line:  # EOF
            break

        # Comments
        if line.startswith("**"):
            line = f.readline()
            continue

        keyword = line.partition(",")[0].strip().replace("*", "").upper()
        if keyword == "NODE":
            points, point_ids, line = _read_nodes(f)
        elif keyword == "ELEMENT":
            if point_ids is None:
                raise ReadError("Expected NODE before ELEMENT")
            params_map = get_param_map(line, required_keys=["TYPE"])
            cell_type, cells_data, ids, sets, line = _read_cells(
                f, params_map, point_ids
            )
            cells.append(CellBlock(cell_type, cells_data))
            cell_ids.append(ids)
            if sets:
                cell_sets_element.update(sets)
                cell_sets_element_order += list(sets.keys())
        elif keyword == "NSET":
            params_map = get_param_map(line, required_keys=["NSET"])
            set_ids, _, line = _read_set(f, params_map)
            name = params_map["NSET"]
            point_sets[name] = np.array(
                [point_ids[point_id] for point_id in set_ids], dtype="int32"
            )
        elif keyword == "ELSET":
            params_map = get_param_map(line, required_keys=["ELSET"])
            set_ids, set_names, line = _read_set(f, params_map)
            name = params_map["ELSET"]
            cell_sets[name] = []
            if set_ids.size:
                for cell_ids_ in cell_ids:
                    cell_sets_ = np.array(
                        [
                            cell_ids_[set_id]
                            for set_id in set_ids
                            if set_id in cell_ids_
                        ],
                        dtype="int32",
                    )
                    cell_sets[name].append(cell_sets_)
            elif set_names:
                for set_name in set_names:
                    if set_name in cell_sets.keys():
                        cell_sets[name].append(cell_sets[set_name])
                    elif set_name in cell_sets_element.keys():
                        cell_sets[name].append(cell_sets_element[set_name])
                    else:
                        raise ReadError(f"Unknown cell set '{set_name}'")
        elif keyword == "SURFACE":
            # `*SURFACE, NAME=..., TYPE=ELEMENT`: each data row is
            # `<element id | elset name>, <face identifier>`. Collected raw and
            # resolved once every element and elset is known.
            params_map = get_param_map(line)
            rows, line = _read_surface_rows(f)
            name = params_map.get("NAME")
            stype = (params_map.get("TYPE") or "ELEMENT").upper()
            if name and stype == "ELEMENT":
                surfaces.append((name, rows))
        elif keyword == "INCLUDE":
            # Splitting line to get external input file path (example: *INCLUDE,INPUT=wInclude_bulk.inp)
            ext_input_file = pathlib.Path(line.split("=")[-1].strip())
            if ext_input_file.exists() is False:
                cd = pathlib.Path(f.name).parent
                ext_input_file = cd / ext_input_file

            # Read contents from external input file into mesh object
            out = read(ext_input_file)

            # Merge contents of external file only if it is containing mesh data
            if len(out.points) > 0:
                points, cells = merge(
                    out,
                    points,
                    cells,
                    point_data,
                    cell_data,
                    field_data,
                    point_sets,
                    cell_sets,
                )

            line = f.readline()
        else:
            # There are just too many Abaqus keywords to explicitly skip them.
            line = f.readline()

    # Parse cell sets defined in ELEMENT
    for i, name in enumerate(cell_sets_element_order):
        # Not sure whether this case would ever happen
        if name in cell_sets.keys():
            cell_sets[name][i] = cell_sets_element[name]
        else:
            cell_sets[name] = []
            for ic in range(len(cells)):
                cell_sets[name].append(
                    cell_sets_element[name] if i == ic else np.array([], dtype="int32")
                )

    mesh = Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
        point_sets=point_sets,
        cell_sets=cell_sets,
    )

    if surfaces:
        _attach_surfaces(mesh, surfaces, cells, cell_ids, cell_sets)

    return mesh


def _read_surface_rows(f):
    """The `<who>, <face>` data rows of a `*SURFACE` block."""
    rows = []
    while True:
        line = f.readline()
        if not line or line.startswith("*"):
            break
        stripped = line.strip()
        if not stripped:
            continue
        parts = [p.strip() for p in stripped.split(",")]
        if len(parts) >= 2 and parts[0] and parts[1]:
            rows.append((parts[0], parts[1]))
    return rows, line


def _attach_surfaces(mesh, surfaces, cells, cell_ids, cell_sets):
    """Turn collected `*SURFACE` blocks into `side` regions on ``mesh``.

    Entries are `(global block-major cell index, local facet index)` pairs, the
    same convention the C++ reader produces -- see doc/regions.md.
    """
    from .._regions import Region

    # element id -> global (block-major) cell index
    elem_index = {}
    base = 0
    for block_ids, block in zip(cell_ids, cells):
        for file_id, local in block_ids.items():
            elem_index[file_id] = base + local
        base += len(block.data)

    bases = np.cumsum([0] + [len(c.data) for c in cells])

    def type_of(g):
        b = int(np.searchsorted(bases, g, side="right") - 1)
        return cells[b].type

    for name, rows in surfaces:
        pairs = []
        for who, face in rows:
            if who in cell_sets:
                members = []
                for b, idx in enumerate(cell_sets[who]):
                    if idx is None:
                        continue
                    members += [int(bases[b]) + int(v) for v in np.asarray(idx)]
            elif who.lstrip("+-").isdigit() and int(who) in elem_index:
                members = [elem_index[int(who)]]
            else:
                continue
            for g in members:
                facet = _face_index(type_of(g), face)
                if facet is not None:
                    pairs.append((g, facet))
        entries = (
            np.asarray(pairs, dtype=np.int64) if pairs else np.empty((0, 2), np.int64)
        )
        mesh.regions.append(Region(name, "side", entries))


def _read_nodes(f):
    points = []
    point_ids = {}
    counter = 0
    while True:
        line = f.readline()
        if not line or line.startswith("*"):
            break
        if line.strip() == "":
            continue

        line = line.strip().split(",")
        point_id, coords = line[0], line[1:]
        point_ids[int(point_id)] = counter
        points.append([float(x) for x in coords])
        counter += 1

    return np.array(points, dtype=float), point_ids, line


def _read_cells(f, params_map, point_ids):
    etype = params_map["TYPE"]
    if etype not in abaqus_to_meshio_type.keys():
        raise ReadError(f"Element type not available: {etype}")

    cell_type = abaqus_to_meshio_type[etype]
    # ElementID + NodesIDs
    num_data = num_nodes_per_cell[cell_type] + 1

    idx = []
    while True:
        line = f.readline()
        if not line or line.startswith("*"):
            break
        line = line.strip()
        if line == "":
            continue
        idx += [int(k) for k in filter(None, line.split(","))]

    # Check for expected number of data
    if len(idx) % num_data != 0:
        raise ReadError("Expected number of data items does not match element type")

    idx = np.array(idx).reshape((-1, num_data))
    cell_ids = dict(zip(idx[:, 0], count(0)))
    cells = np.array([[point_ids[node] for node in elem] for elem in idx[:, 1:]])

    cell_sets = (
        {params_map["ELSET"]: np.arange(len(cells), dtype="int32")}
        if "ELSET" in params_map.keys()
        else {}
    )

    return cell_type, cells, cell_ids, cell_sets, line


def merge(
    mesh, points, cells, point_data, cell_data, field_data, point_sets, cell_sets
):
    """
    Merge Mesh object into existing containers for points, cells, sets, etc..

    :param mesh:
    :param points:
    :param cells:
    :param point_data:
    :param cell_data:
    :param field_data:
    :param point_sets:
    :param cell_sets:
    :type mesh: Mesh
    """
    ext_points = np.array([p for p in mesh.points])

    if len(points) > 0:
        new_point_id = points.shape[0]
        # new_cell_id = len(cells) + 1
        points = np.concatenate([points, ext_points])
    else:
        # new_cell_id = 0
        new_point_id = 0
        points = ext_points

    new_block_id = len(cells)
    cnt = 0
    for c in mesh.cells:
        new_data = np.array([d + new_point_id for d in c.data])
        cells.append(CellBlock(c.type, new_data))
        cnt += 1

    # The following aren't currently included in the abaqus parser, and are therefore
    # excluded?
    # point_data.update(mesh.point_data)
    # cell_data.update(mesh.cell_data)
    # field_data.update(mesh.field_data)

    # Update point and cell sets to account for change in cell and point ids
    for key, val in mesh.point_sets.items():
        point_sets[key] = [x + new_point_id for x in val]

    # Cell sets: the included file's blocks were appended after the ones already
    # present, so its per-block lists are padded on the left by that many empty
    # blocks, and every set declared earlier is padded on the right. (This used
    # to be a TODO that silently dropped them; the C++ reader carries them, and
    # the two paths must agree.)
    for key, val in cell_sets.items():
        cell_sets[key] = list(val) + [np.array([], dtype="int32")] * cnt
    pad = [np.array([], dtype="int32")] * new_block_id
    for key, val in mesh.cell_sets.items():
        cell_sets[key] = pad + [np.asarray(v, dtype="int32") for v in val]

    return points, cells


def get_param_map(word, required_keys=None):
    """
    get the optional arguments on a line

    Example
    -------
    >>> word = 'elset,instance=dummy2,generate'
    >>> params = get_param_map(word, required_keys=['instance'])
    params = {
        'elset' : None,
        'instance' : 'dummy2,
        'generate' : None,
    }
    """
    if required_keys is None:
        required_keys = []
    words = word.split(",")
    param_map = {}
    for wordi in words:
        if "=" not in wordi:
            key = wordi.strip().upper()
            value = None
        else:
            sword = wordi.split("=")
            if len(sword) != 2:
                raise ReadError(sword)
            key = sword[0].strip().upper()
            value = sword[1].strip()
        param_map[key] = value

    msg = ""
    for key in required_keys:
        if key not in param_map:
            msg += f"{key} not found in {word}\n"
    if msg:
        raise RuntimeError(msg)
    return param_map


def _read_set(f, params_map):
    set_ids = []
    set_names = []
    while True:
        line = f.readline()
        if not line or line.startswith("*"):
            break
        if line.strip() == "":
            continue

        line = line.strip().strip(",").split(",")
        if line[0].isnumeric():
            set_ids += [int(k) for k in line]
        else:
            set_names.append(line[0])

    set_ids = np.array(set_ids, dtype="int32")
    if "GENERATE" in params_map:
        if len(set_ids) != 3:
            raise ReadError(set_ids)
        set_ids = np.arange(set_ids[0], set_ids[1] + 1, set_ids[2], dtype="int32")
    return set_ids, set_names, line


def write(
    filename, mesh: Mesh, float_fmt: str = ".16e", translate_cell_names: bool = True
) -> None:
    with open_file(filename, "wt") as f:
        f.write("*HEADING\n")
        f.write("Abaqus DataFile Version 6.14\n")
        f.write(_provenance.render_lines(_provenance.SlotTier.BLOCK, ""))
        f.write("*NODE\n")
        fmt = ", ".join(["{}"] + ["{:" + float_fmt + "}"] * mesh.points.shape[1]) + "\n"
        for k, x in enumerate(mesh.points):
            f.write(fmt.format(k + 1, *x))
        eid = 0
        for cell_block in mesh.cells:
            cell_type = cell_block.type
            node_idcs = cell_block.data
            name = (
                meshio_to_abaqus_type[cell_type] if translate_cell_names else cell_type
            )
            f.write(f"*ELEMENT, TYPE={name}\n")
            for row in node_idcs:
                eid += 1
                nids_strs = (str(nid + 1) for nid in row.tolist())
                f.write(str(eid) + "," + ",".join(nids_strs) + "\n")

        nnl = 8
        offset = 0
        for ic in range(len(mesh.cells)):
            for k, v in mesh.cell_sets.items():
                if len(v[ic]) > 0:
                    els = [str(i + 1 + offset) for i in v[ic]]
                    f.write(f"*ELSET, ELSET={k}\n")
                    f.write(
                        ",\n".join(
                            ",".join(els[i : i + nnl]) for i in range(0, len(els), nnl)
                        )
                        + "\n"
                    )
            offset += len(mesh.cells[ic].data)

        for k, v in mesh.point_sets.items():
            nds = [str(i + 1) for i in v]
            f.write(f"*NSET, NSET={k}\n")
            f.write(
                ",\n".join(",".join(nds[i : i + nnl]) for i in range(0, len(nds), nnl))
                + "\n"
            )

        # Side regions -> *SURFACE. They have no point_sets/cell_sets
        # equivalent, so this is the only way they leave the mesh.
        bases = np.cumsum([0] + [len(c.data) for c in mesh.cells])
        for region in getattr(mesh, "regions", []):
            if region.kind != "side":
                continue
            f.write(f"*SURFACE, NAME={region.name}, TYPE=ELEMENT\n")
            for g, facet in np.asarray(region.entries).reshape(-1, 2):
                b = int(np.searchsorted(bases, g, side="right") - 1)
                if b < 0 or b >= len(mesh.cells):
                    continue
                face = _face_name(mesh.cells[b].type, int(facet))
                if face is not None:
                    f.write(f"{int(g) + 1}, {face}\n")

        # https://github.com/nschloe/meshio/issues/747#issuecomment-643479921
        # f.write("*END")


# NOTE: format registration now lives in meshioplusplus/abaqus/__init__.py, which wraps
# the reader/writer above with the C++-backed fast paths.
