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
from .._regions import block_bases
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
    variables, zones, _ = _load(filename)
    return _zones.metadata(zones, variables)[2]


def _zone_title(mesh, bases, block):
    """An exactly-matching Cell region's name, else ``block_<block>``.

    Mirrors the C++ writer: a region counts only when its (sorted, canonical)
    entries are precisely the global cell-index range of this block.
    """
    lo, hi = bases[block], bases[block + 1]
    for region in mesh.regions:
        if region.kind != "cell" or len(region.entries) == 0:
            continue
        if len(region.entries) != hi - lo:
            continue
        if region.entries[0] == lo and region.entries[-1] == hi - 1:
            return region.name
    return f"block_{block}"


def write(filename, mesh):
    # One Tecplot ZONE per cell block -- no single-type restriction, no
    # 2D/3D padding hack. Zone 1 carries the coordinates and nodal fields;
    # later zones reuse them through VARSHARELIST rather than duplicating
    # the (unwelded, shared) point array. A cell-data array is expected one
    # entry per cell block (the meshio convention); a block outside
    # ``cell_blocks`` (an unsupported type) simply never claims it.
    cell_blocks = []
    for ic, c in enumerate(mesh.cells):
        if c.type in meshio_only:
            cell_blocks.append(ic)
        else:
            warn(
                f"Tecplot does not support cell type '{c.type}'. Skipping cell block {ic}."
            )
    if not cell_blocks:
        raise WriteError("No cell type supported by Tecplot in mesh")

    num_nodes = len(mesh.points)
    dim = mesh.points.shape[1]

    variables = ["X", "Y"]
    shared_data = [mesh.points[:, 0], mesh.points[:, 1]]
    if dim == 3:
        variables.append("Z")
        shared_data.append(mesh.points[:, 2])

    for k, v in mesh.point_data.items():
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
    num_shared = len(variables)

    # Cell-centred variables, component-expanded: `cell_var_blocks[j]` maps
    # block index -> that component's 1-D array, for variable `num_shared+j`.
    cell_var_blocks = []
    for k, v in mesh.cell_data.items():
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

    with open_file(filename, "w") as f:
        f.write(f'TITLE = "{_provenance.lines(_provenance.SlotTier.SINGLE_LINE)[0]}"\n')
        variables_str = ", ".join(f'"{var}"' for var in variables)
        f.write(f"VARIABLES = {variables_str}\n")

        for bi, ic in enumerate(cell_blocks):
            cell = mesh.cells[ic]
            zone_type = meshio_to_tecplot_type[cell.type]
            order = meshio_to_tecplot_order[cell.type]
            num_cells = len(cell.data)
            title = _zone_title(mesh, bases, ic)

            present = [j for j, pb in enumerate(cell_var_blocks) if ic in pb]
            passive = [j for j, pb in enumerate(cell_var_blocks) if ic not in pb]

            f.write(
                f'ZONE T = "{title}", NODES = {num_nodes}, ELEMENTS = {num_cells},\n'
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

            for row in cell.data[:, order]:
                f.write(" ".join(str(c) for c in row + 1) + "\n")


def _write_table(f, data, ncol=20):
    nrow = len(data) // ncol
    lines = np.split(data, np.full(nrow, ncol).cumsum())
    for line in lines:
        if len(line):
            f.write(" ".join(str(l) for l in line) + "\n")


# NOTE: format registration now lives in meshioplusplus/tecplot/__init__.py, which wraps
# the reader/writer above with the C++-backed fast paths.
