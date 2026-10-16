"""
I/O for Tecplot ASCII data format (and reading binary ``.plt``), cf.
<https://github.com/su2code/SU2/raw/master/externals/tecio/360_data_format_guide.pdf>,
<http://paulbourke.net/dataformats/tp/>.
"""

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError, WriteError
from .._files import is_buffer, open_file
from .._regions import block_bases, canonical_regions
from . import _ascii, _plt, _zones
from ._plt import is_plt

meshio_to_tecplot_type = {
    "line": "FELINESEG",
    "triangle": "FETRIANGLE",
    "quad": "FEQUADRILATERAL",
    "tetra": "FETETRAHEDRON",
    "pyramid": "FEBRICK",
    "wedge": "FEBRICK",
    "hexahedron": "FEBRICK",
}


meshio_only = set(meshio_to_tecplot_type.keys())


meshio_to_tecplot_order = {
    "line": [0, 1],
    "triangle": [0, 1, 2],
    "quad": [0, 1, 2, 3],
    "tetra": [0, 1, 2, 3],
    "pyramid": [0, 1, 2, 3, 4, 4, 4, 4],
    "wedge": [0, 1, 4, 3, 2, 2, 5, 5],
    "hexahedron": [0, 1, 2, 3, 4, 5, 6, 7],
}


def _zone_type(cell_type):
    if cell_type in meshio_to_tecplot_type:
        return meshio_to_tecplot_type[cell_type]
    if cell_type.startswith("polygon"):
        return "FEPOLYGON"
    if cell_type.startswith("polyhedron"):
        return "FEPOLYHEDRON"
    return None


def _face_map(cell):
    """``tecplot_write_faces``: a polygon or polyhedron block's faces in
    first-seen order, each shared face once; the first cell to use a face is
    its left element (its faces are outward), the second its right. Elements
    1-based, 0 = none."""
    faces, left, right, seen = [], [], [], {}

    def add(face, c):
        key = tuple(sorted(face))
        k = seen.get(key)
        if k is not None and right[k] == 0:
            right[k] = c
            return
        seen[key] = len(faces)
        faces.append(face)
        left.append(c)
        right.append(0)

    for c, cell_nodes in enumerate(cell.data, start=1):
        if cell.type.startswith("polyhedron"):
            for face in cell_nodes:
                add([int(v) for v in face], c)
        else:
            row = [int(v) for v in cell_nodes]
            for j in range(len(row)):
                add([row[j], row[(j + 1) % len(row)]], c)
    return faces, left, right


def _load(filename):
    if not is_buffer(filename, "r"):
        with open(filename, "rb") as f:
            head = f.read(8)
        if is_plt(head):
            return _plt.load(filename)
    return _ascii.load(filename)


def _resolve_step(time_step, count):
    step = time_step + count if time_step < 0 else time_step
    if not 0 <= step < count:
        raise ReadError(
            f"meshio++: time step {time_step} is out of range: this file has {count} "
            + ("step" if count == 1 else "steps")
        )
    return step


def read(filename, time_step=0):
    """Reads an ASCII (``.dat``/``.tec``) or binary (``.plt``) Tecplot file.

    Every zone of the selected step becomes a cell block and a Cell region
    (``time_step`` picks a SOLUTIONTIME of a transient file).
    """
    variables, zones, source = _load(filename)
    steps = _zones.timeline(zones)
    return _zones.build_step(
        steps[_resolve_step(time_step, len(steps))], zones, variables, source
    )


def time_values(filename):
    """The SOLUTIONTIME of each step (empty for a static file)."""
    variables, zones, source = _load(filename)
    return _zones.metadata(zones, variables, source)[2]


def _zone_title(mesh, bases, block):
    """An exactly-matching Cell region's name, else ``block_<block>``.

    Mirrors the C++ writer: a region counts only when its (sorted, canonical)
    entries are precisely the global cell-index range of this block.
    """
    lo, hi = bases[block], bases[block + 1]
    for region in canonical_regions(mesh):
        if region.kind != "cell" or len(region.entries) == 0:
            continue
        if len(region.entries) != hi - lo:
            continue
        if region.entries[0] == lo and region.entries[-1] == hi - 1:
            return region.name
    return f"block_{block}"


def write_plan(mesh):
    """The zones and variables every Tecplot writer lays out: ``(cell_blocks,
    variables, shared_data, cell_var_blocks, bases)``.

    One Tecplot ZONE per cell block -- no single-type restriction, no 2D/3D
    padding hack. Zone 1 carries the coordinates and nodal fields; later zones
    reuse them through VARSHARELIST rather than duplicating the (unwelded,
    shared) point array. A cell-data array is expected one entry per cell block
    (the meshio convention); a block outside ``cell_blocks`` (an unsupported
    type) simply never claims it.
    """
    cell_blocks = []
    for ic, c in enumerate(mesh.cells):
        if _zone_type(c.type) is not None:
            cell_blocks.append(ic)
        else:
            warn(
                f"Tecplot does not support cell type '{c.type}'. Skipping cell block {ic}."
            )
    if not cell_blocks:
        raise WriteError("No cell type supported by Tecplot in mesh")

    dim = mesh.points.shape[1]

    variables = ["X", "Y"]
    shared_data = [mesh.points[:, 0], mesh.points[:, 1]]
    if dim == 3:
        variables.append("Z")
        shared_data.append(mesh.points[:, 2])

    for k, v in sorted(mesh.point_data.items()):
        if k in {"X", "Y", "Z", "x", "y", "z"}:
            warn(f"Skipping point data '{k}'.")
            continue
        if v.ndim == 1:
            variables.append(k)
            shared_data.append(v)
        elif v.ndim == 2:
            for i, vv in enumerate(v.T):
                variables.append(f"{k}_{i}")
                shared_data.append(vv)

    # Cell-centred variables, component-expanded: `cell_var_blocks[j]` maps
    # block index -> that component's 1-D array, for variable `num_shared+j`.
    cell_var_blocks = []
    for k, v in sorted(mesh.cell_data.items()):
        if k in {"X", "Y", "Z", "x", "y", "z"}:
            warn(f"Skipping cell data '{k}'.")
            continue
        first = next((np.asarray(v[ic]) for ic in cell_blocks if ic < len(v)), None)
        if first is None:
            continue
        ncomp = first.shape[1] if first.ndim == 2 else 1
        for c in range(ncomp):
            variables.append(k if ncomp == 1 else f"{k}_{c}")
            per_block = {}
            for ic in cell_blocks:
                if ic >= len(v):
                    continue
                arr = np.asarray(v[ic])
                per_block[ic] = arr[:, c] if arr.ndim == 2 else arr
            cell_var_blocks.append(per_block)

    bases = block_bases(mesh.cells)
    return cell_blocks, variables, shared_data, cell_var_blocks, bases


def write(filename, mesh):
    cell_blocks, variables, shared_data, cell_var_blocks, bases = write_plan(mesh)
    num_nodes = len(mesh.points)
    num_shared = len(shared_data)

    if str(filename).lower().endswith(".plt"):
        _write_plt(
            filename,
            mesh,
            cell_blocks,
            variables,
            shared_data,
            cell_var_blocks,
            bases,
        )
        return

    with open_file(filename, "w") as f:
        f.write(f'TITLE = "{_provenance.lines(_provenance.SlotTier.SINGLE_LINE)[0]}"\n')
        variables_str = ", ".join(f'"{var}"' for var in variables)
        f.write(f"VARIABLES = {variables_str}\n")

        for bi, ic in enumerate(cell_blocks):
            cell = mesh.cells[ic]
            zone_type = _zone_type(cell.type)
            face_based = zone_type in ("FEPOLYGON", "FEPOLYHEDRON")
            num_cells = len(cell.data)
            title = _zone_title(mesh, bases, ic)

            present = [j for j, pb in enumerate(cell_var_blocks) if ic in pb]
            passive = [j for j, pb in enumerate(cell_var_blocks) if ic not in pb]

            f.write(
                f'ZONE T = "{title}", NODES = {num_nodes}, ELEMENTS = {num_cells},\n'
            )
            if face_based:
                faces, left, right = _face_map(cell)
                f.write(
                    f"FACES = {len(faces)}, TOTALNUMFACENODES = {sum(map(len, faces))},\n"
                    "NUMCONNECTEDBOUNDARYFACES = 0, TOTALNUMBOUNDARYCONNECTIONS = 0,\n"
                )
            f.write(f"DATAPACKING = BLOCK, ZONETYPE = {zone_type}")
            if bi > 0:
                f.write(f",\nVARSHARELIST = ([1-{num_shared}] = 1)")
            if present:
                rng = ",".join(str(num_shared + j + 1) for j in present)
                f.write(f",\nVARLOCATION = ([{rng}] = CELLCENTERED)")
            if passive:
                rng = ",".join(str(num_shared + j + 1) for j in passive)
                f.write(f",\nPASSIVEVARLIST = ([{rng}])")
            f.write("\n")

            if bi == 0:
                for arr in shared_data:
                    _write_table(f, arr)
            for j in present:
                _write_table(f, cell_var_blocks[j][ic])

            if face_based:
                if zone_type == "FEPOLYHEDRON":
                    _write_ints(f, [len(face) for face in faces])
                for face in faces:
                    f.write(" ".join(str(v + 1) for v in face) + "\n")
                _write_ints(f, left)
                _write_ints(f, right)
                continue
            order = meshio_to_tecplot_order[cell.type]
            for row in cell.data[:, order]:
                f.write(" ".join(str(c) for c in row + 1) + "\n")


_PLT_ZONE_TYPES = {
    "FELINESEG": 1,
    "FETRIANGLE": 2,
    "FEQUADRILATERAL": 3,
    "FETETRAHEDRON": 4,
    "FEBRICK": 5,
    "FEPOLYGON": 6,
    "FEPOLYHEDRON": 7,
}


def _write_plt(
    filename, mesh, cell_blocks, variables, shared_data, cell_var_blocks, bases
):
    """The same zones as the ASCII writer, as a binary ``#!TDV112`` file.

    One zone per cell block; zone 1 owns the nodal variables and the later
    zones share them (and carry only their own cell-centred values, the others
    passive). Values are written as doubles, so nothing narrows; connectivity is
    0-based. The layout is the Data Format Guide's appendix A, which
    ``_plt.py`` reads back.
    """
    num_nodes = len(mesh.points)
    nvar = len(variables)
    num_shared = len(shared_data)
    out = bytearray()

    def i32(*values):
        out.extend(np.asarray(values, dtype="<i4").tobytes())

    def f32(value):
        out.extend(np.float32(value).tobytes())

    def f64(*values):
        out.extend(np.asarray(values, dtype="<f8").tobytes())

    def string(text):
        out.extend(np.asarray([ord(c) for c in text] + [0], dtype="<i4").tobytes())

    out.extend(b"#!TDV112")
    i32(1, 0)  # byte order, full file
    string(_provenance.lines(_provenance.SlotTier.SINGLE_LINE)[0])
    i32(nvar)
    for name in variables:
        string(name)

    zones = []
    for bi, ic in enumerate(cell_blocks):
        cell = mesh.cells[ic]
        zone_type = _zone_type(cell.type)
        face_based = zone_type in ("FEPOLYGON", "FEPOLYHEDRON")
        faces = _face_map(cell) if face_based else None
        present = [j for j, pb in enumerate(cell_var_blocks) if ic in pb]
        zones.append((ic, cell, zone_type, faces, present))
        f32(299.0)
        string(_zone_title(mesh, bases, ic))
        i32(-1, -1)  # parent zone, strand
        f64(0.0)  # solution time
        i32(-1, _PLT_ZONE_TYPES[zone_type])  # color, zone type
        if cell_var_blocks:
            located = [0] * num_shared + [1] * len(cell_var_blocks)
            i32(1, *located)
        else:
            i32(0)
        i32(0, 0)  # raw and miscellaneous face neighbours
        i32(num_nodes)
        if face_based:
            nodes = sum(len(face) for face in faces[0])
            i32(len(faces[0]), nodes, 0, 0)
        i32(len(cell.data), 0, 0, 0)  # cells, I/J/K cell dimensions
        i32(0)  # no auxiliary data
    f32(357.0)

    for bi, (ic, cell, zone_type, faces, present) in enumerate(zones):
        f32(299.0)
        i32(*([2] * nvar))  # every variable a double
        passive = [
            num_shared + j for j in range(len(cell_var_blocks)) if j not in present
        ]
        if passive:
            flags = [0] * nvar
            for v in passive:
                flags[v] = 1
            i32(1, *flags)
        else:
            i32(0)
        if bi > 0:
            i32(1, *([0] * num_shared + [-1] * len(cell_var_blocks)))
        else:
            i32(0)
        i32(-1)  # own connectivity
        owned = []
        if bi == 0:
            owned += [np.asarray(col, dtype=np.float64) for col in shared_data]
        owned += [np.asarray(cell_var_blocks[j][ic], dtype=np.float64) for j in present]
        for column in owned:
            finite = column[~np.isnan(column)]
            f64(*((finite.min(), finite.max()) if finite.size else (0.0, 0.0)))
        for column in owned:
            out.extend(np.ascontiguousarray(column, dtype="<f8").tobytes())
        if faces is not None:
            face_list, left, right = faces
            if zone_type == "FEPOLYHEDRON":
                offsets = np.concatenate([[0], np.cumsum([len(f) for f in face_list])])
                i32(*offsets)
            for face in face_list:
                i32(*face)
            i32(*[v - 1 for v in left])
            i32(*[v - 1 for v in right])
            continue
        order = meshio_to_tecplot_order[cell.type]
        out.extend(np.ascontiguousarray(cell.data[:, order], dtype="<i4").tobytes())

    with open_file(filename, "wb") as f:
        f.write(bytes(out))


def _write_ints(f, ints, ncol=20):
    if not len(ints):
        f.write("\n")
    for s in range(0, len(ints), ncol):
        f.write(" ".join(str(v) for v in ints[s : s + ncol]) + "\n")


def _write_table(f, data, ncol=20):
    nrow = len(data) // ncol
    lines = np.split(data, np.full(nrow, ncol).cumsum())
    for line in lines:
        if len(line):
            f.write(" ".join(str(l) for l in line) + "\n")


# NOTE: format registration now lives in meshioplusplus/tecplot/__init__.py, which wraps
# the reader/writer above with the C++-backed fast paths.
