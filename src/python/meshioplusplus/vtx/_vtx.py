"""DOLFINx VTX output: an ADIOS2 ``.bp`` directory, read through ``adios2``.

``dolfinx.io.VTXWriter`` writes a ``vtk.xml`` attribute (the schema ParaView's
``ADIOS2VTXReader`` follows) and, per step, the variables it names: ``geometry``
(n x 3), ``connectivity`` (cells x (1 + nodes), each row prefixed with its node
count), ``types``, the point and cell arrays and the ``step`` time. Every MPI
rank writes one block with its own point numbering.

The Python twin of ``src/cpp/src/formats/vtx.cpp``; see ``doc/formats/vtx.md``.
"""

import logging
import xml.etree.ElementTree as ET

import numpy as np

from .._exceptions import ReadError
from .._mesh import Mesh
from .._vtk_common import vtk_cells_from_data

TIME_KEY = "meshio:time"
_GHOST_ARRAYS = ("vtkGhostType", "vtkOriginalPointIds")

# VTK_LAGRANGE_* cells whose node order is a linear or quadratic VTK cell's:
# (Lagrange type, node count) -> VTK type.
_LOWER = {
    (68, 2): 3,
    (68, 3): 21,
    (69, 3): 5,
    (69, 6): 22,
    (70, 4): 9,
    (70, 9): 28,
    (71, 4): 10,
    (71, 10): 24,
    (72, 8): 12,
    (73, 6): 13,
    (74, 5): 14,
}

_log = logging.getLogger("meshioplusplus")


def _adios2():
    from .._interop import _require

    return _require("adios2", "adios2", "reading a DOLFINx VTX (.bp) file")


def _fail(message):
    raise ReadError(f"VTX (.bp): {message}")


def _parse_schema(xml):
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        _fail("the vtk.xml schema is not valid XML")
    if root.tag != "VTKFile" or root.get("type") != "UnstructuredGrid":
        _fail(
            "only an UnstructuredGrid schema is read, this one is "
            f"'{root.get('type', '')}'"
        )
    piece = root.find("UnstructuredGrid/Piece")
    if piece is None:
        _fail("the vtk.xml schema has no Piece")

    def variable(da):
        text = (da.text or "").strip()
        return text or da.get("Name", "")

    schema = {
        "geometry": "geometry",
        "connectivity": "connectivity",
        "types": "types",
        "time": "",
        "arrays": [],
    }
    da = piece.find("Points/DataArray")
    if da is not None:
        schema["geometry"] = variable(da)
    for da in piece.findall("Cells/DataArray"):
        if da.get("Name") in ("connectivity", "types"):
            schema[da.get("Name")] = variable(da)
    for group, is_cell in (("PointData", False), ("CellData", True)):
        for da in piece.findall(f"{group}/DataArray"):
            name = da.get("Name", "")
            if not name:
                continue
            if name == "TIME":
                schema["time"] = variable(da)
                continue
            schema["arrays"].append((name, variable(da), is_cell))
    return schema


def _attribute_string(reader, name):
    value = reader.read_attribute_string(name)
    if isinstance(value, (list, tuple)):
        value = "".join(value)
    return value


def _scan(filename):
    """Schema, the variables each step holds, and the step times."""
    adios2 = _adios2()
    present, times = [], []
    schema = None
    with adios2.Stream(str(filename), "r") as stream:
        for _ in stream.steps():
            if schema is None:
                if "vtk.xml" not in stream.available_attributes():
                    _fail(
                        f"'{filename}' has no vtk.xml attribute: only DOLFINx "
                        "VTXWriter output is read (adios4dolfinx checkpoints and "
                        "Fides output are other layouts)"
                    )
                schema = _parse_schema(_attribute_string(stream, "vtk.xml"))
            names = set(stream.available_variables())
            t = schema["time"]
            times.append(
                float(np.asarray(stream.read(t)).ravel()[0])
                if t and t in names
                else float(len(present))
            )
            present.append(names)
    if schema is None:
        _fail(f"'{filename}' holds no step")
    return schema, present, times


def _relative(present, name, step):
    return sum(1 for s in range(step) if name in present[s])


def _last_with(present, name, step):
    for s in range(step, -1, -1):
        if name in present[s]:
            return s
    return None


def _blocks(reader, name, rel, rows=None):
    """Every block of a variable at one (relative) step, as arrays.

    With ``rows``, only the blocks matched to the rank blocks by length, in
    order, are read: DOLFINx 0.11 writes stray blocks whose count does not
    describe their payload, so an unmatched block is never read."""
    var = reader.inquire_variable(name)
    if var is None:
        return []
    infos = reader.all_blocks_info(name)[rel]
    if var.shape() == [] and all(info.get("IsValue") == "True" for info in infos):
        # A global value (one instance per writer): its block list does not carry
        # the value in every engine (BP4 records none), so it is read; the read
        # gives the first instance.
        var.set_step_selection([rel, 1])
        return [np.atleast_1d(np.asarray(reader.read(var)).ravel()[:1])]
    wanted = [rows is None] * len(infos)
    if rows is not None:
        j = 0
        for index, info in enumerate(infos):
            count = [int(c) for c in info["Count"].split(",") if c.strip()]
            n = 1 if info.get("IsValue") == "True" or not count else count[0]
            if j < len(rows) and n == rows[j]:
                wanted[index] = True
                j += 1
            else:
                _log.debug(
                    "VTX: skipping block %d of '%s': it matches no rank's size",
                    index,
                    name,
                )
    out = []
    # A block is selected by its position in the step's list: BlockID restarts
    # at 0 for every writer (MPI rank).
    for index, info in enumerate(infos):
        if not wanted[index]:
            continue
        if info.get("IsValue") == "True":
            # A local value: the Python bindings report it as Min (= Max).
            value = info.get("Value", info.get("Min"))
            out.append(np.atleast_1d(np.asarray(float(value))))
            continue
        var.set_step_selection([rel, 1])
        var.set_block_selection(index)
        count = [int(c) for c in info["Count"].split(",") if c.strip()]
        data = np.asarray(reader.read(var))
        out.append(data.reshape(count) if count else np.atleast_1d(data))
    return out


def _block_rows(reader, name, rel):
    rows = []
    for info in reader.all_blocks_info(name)[rel]:
        count = [int(c) for c in info["Count"].split(",") if c.strip()]
        rows.append(1 if info.get("IsValue") == "True" or not count else count[0])
    return rows


def _clean(reader, present, name, step, rows):
    return name in present[step] and _block_rows(
        reader, name, _relative(present, name, step)
    ) == list(rows)


def _concat_matching(blocks, counts, name):
    """Blocks matched to the rank blocks by length, in order; None when they do
    not cover every rank."""
    picked = []
    j = 0
    for b in blocks:
        if j < len(counts) and b.shape[0] == counts[j]:
            picked.append(b)
            j += 1
        else:
            _log.debug(
                "VTX: skipping a block of '%s' that matches no rank's size", name
            )
    if j != len(counts) or not picked:
        return None
    width = {b.shape[1:] for b in picked}
    if len(width) != 1:
        return None
    out = np.concatenate(picked, axis=0)
    if out.ndim == 2 and out.shape[1] == 1:
        out = out[:, 0]
    return out


def _weld_ghosts(points, conn, point_data):
    ids = point_data.get("vtkOriginalPointIds")
    if ids is None:
        _log.warning("VTX: no vtkOriginalPointIds, so ghost points cannot be welded")
        return points, conn, point_data
    ghost = point_data.get("vtkGhostType")
    n = len(ids)
    is_ghost = (
        np.zeros(n, dtype=bool) if ghost is None else np.asarray(ghost).ravel() != 0
    )
    owner = {}
    for wanted in (False, True):
        for i in range(n):
            if is_ghost[i] == wanted:
                owner.setdefault(int(ids[i]), i)
    rep = np.array([owner[int(v)] for v in ids], dtype=np.int64)
    keep = np.flatnonzero(rep == np.arange(n))
    new_index = np.full(n, -1, dtype=np.int64)
    new_index[keep] = np.arange(len(keep))
    conn = new_index[rep[conn]]
    point_data = {k: v[keep] for k, v in point_data.items() if k not in _GHOST_ARRAYS}
    return points[keep], conn, point_data


def _read_step(filename, step_request, points_only, arrays, ghosts):
    adios2 = _adios2()
    schema, present, times = _scan(filename)
    nsteps = len(present)
    if not -nsteps <= step_request < nsteps:
        raise ReadError(
            f"meshio++: time step {step_request} requested, but this file has "
            f"{nsteps} time steps"
        )
    step = step_request % nsteps
    geometry, connectivity, types_name = (
        schema["geometry"],
        schema["connectivity"],
        schema["types"],
    )
    mesh_step = _last_with(present, geometry, step)
    if mesh_step is None:
        _fail(f"no step up to {step} holds the mesh ('{geometry}')")
    conn_step = _last_with(present, connectivity, step)
    types_step = _last_with(present, types_name, step)
    if conn_step is None or types_step is None:
        _fail(f"no step up to {step} holds the cells")

    drop = ghosts == "drop"
    with adios2.FileReader(str(filename)) as reader:

        def blocks_at(name, s, rows=None):
            return _blocks(reader, name, _relative(present, name, s), rows)

        geom = blocks_at(geometry, mesh_step)
        point_counts = [g.shape[0] for g in geom]
        points = np.concatenate(geom, axis=0)
        if points.ndim != 2:
            _fail(f"'{geometry}' is not an (n, dim) array")

        conn_blocks = blocks_at(connectivity, conn_step)
        if len(conn_blocks) != len(geom):
            _fail(
                f"'{connectivity}' has {len(conn_blocks)} blocks, "
                f"'{geometry}' {len(geom)}"
            )
        type_blocks = blocks_at(types_name, types_step)
        if not type_blocks:
            _fail(f"'{types_name}' holds no value")
        conn, offsets, types, cell_counts = [], [], [], []
        base, end = 0, 0
        for b, c in enumerate(conn_blocks):
            c = np.asarray(c, dtype=np.int64)
            if c.ndim != 2 or c.shape[0] == 0:
                cell_counts.append(0)
                base += point_counts[b]
                continue
            rows, width = c.shape
            if width < 2:
                _fail(f"'{connectivity}' is not a (cells, 1 + nodes) array")
            if np.any(c[:, 0] != width - 1):
                _fail(f"a row of '{connectivity}' counts other than {width - 1} nodes")
            nodes = c[:, 1:]
            if nodes.size and (nodes.min() < 0 or nodes.max() >= point_counts[b]):
                _fail(
                    f"'{connectivity}' names a point outside a block with "
                    f"{point_counts[b]}"
                )
            cell_counts.append(rows)
            conn.append(nodes.ravel() + base)
            offsets.append(end + (width - 1) * np.arange(1, rows + 1, dtype=np.int64))
            end += (width - 1) * rows
            t = np.asarray(
                type_blocks[b]
                if len(type_blocks) == len(conn_blocks)
                else type_blocks[0]
            ).ravel()
            t = t if t.size == rows else np.full(rows, t[0])
            types.append(
                np.array([_LOWER.get((int(v), width - 1), int(v)) for v in t], np.int64)
            )
            base += point_counts[b]
        conn = np.concatenate(conn) if conn else np.zeros(0, np.int64)
        offsets = np.concatenate(offsets) if offsets else np.zeros(0, np.int64)
        types = np.concatenate(types) if types else np.zeros(0, np.int64)

        point_data, cell_data_raw = {}, {}
        for name, var, is_cell in schema["arrays"]:
            is_ghost_array = name in _GHOST_ARRAYS
            wanted = not points_only and (arrays is None or name in arrays)
            if not wanted and not (drop and is_ghost_array and not is_cell):
                continue
            if var not in present[step]:
                continue
            rows = cell_counts if is_cell else point_counts
            # DOLFINx 0.11 writes the ghost arrays twice in the first step of a
            # `reuse` file, the second time as one-value blocks whose payload is
            # not what their count says (ADIOS2 then misreads the next rank's
            # block too): a ghost array whose blocks do not match the rank blocks
            # one to one is taken from the nearest step where they do; any other
            # array is dropped.
            source = step
            if not _clean(reader, present, var, step, rows):
                source = None
                if is_ghost_array:
                    for d in range(1, nsteps):
                        for s in (step + d, step - d):
                            if 0 <= s < nsteps and _clean(
                                reader, present, var, s, rows
                            ):
                                source = s
                                break
                        if source is not None:
                            break
                if source is None:
                    _log.warning(
                        "VTX: the blocks of '%s' match no rank's %s count and it is "
                        "not read",
                        name,
                        "cell" if is_cell else "point",
                    )
                    continue
            joined = _concat_matching(blocks_at(var, source, rows), rows, name)
            if joined is None:
                _log.warning(
                    "VTX: '%s' matches no rank's %s count and is not read",
                    name,
                    "cell" if is_cell else "point",
                )
                continue
            (cell_data_raw if is_cell else point_data)[name] = joined

    if drop:
        points, conn, point_data = _weld_ghosts(points, conn, point_data)
    cells, cell_data = vtk_cells_from_data(conn, offsets, types, cell_data_raw)
    mesh = Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data={TIME_KEY: np.array([times[step]], dtype=np.float64)},
    )
    mesh.time_values = list(times)  # read_metadata's side channel
    return mesh


def read(filename, points_only=False, arrays=None, time_step=0, ghosts="keep"):
    if ghosts not in ("keep", "drop"):
        raise ValueError(f"meshio++: ghosts must be 'keep' or 'drop', got '{ghosts}'")
    return _read_step(filename, time_step, points_only, arrays, ghosts)


def time_values(filename):
    return _scan(filename)[2]
