"""Write Tecplot SZL (``.szplt``) files through a user-supplied shared TecIO.

The Python twin of ``write_szplt`` in ``src/cpp/src/formats/tecplot.cpp``: both
lay the zones out as the ASCII and ``.plt`` writers do (``tecplot._tecplot``'s
``write_plan``): one zone per cell block, a cell-centred variable passive in a
zone that has no values for it. Wedges and pyramids are degenerate bricks, as
there. Each variable is written in the TecIO type that holds its array's dtype
exactly: single precision stays single, 16-bit and byte integers stay small, and
what no TecIO type holds (64-bit integers, unsigned 32-bit) is a double.

A static file's tetrahedron and brick zones whose cells carry ``partition:part``
values of several parts are written as partitioned zones, one partition per part:
a partition holds the nodes its cells use, the lowest part to use a node owns it
and the others hold it as a ghost node (TecIO takes partitions of 3D FE zones only,
and numbers the zone's nodes by partition, each partition's owned nodes in order,
and its cells by partition too). A zone that shares its connectivity from a
partitioned one would need the partitions again, so a series is written whole.

TecIO shares a variable between zones only when it shares their connectivity
too, so the blocks of a mesh cannot share one point array as the other writers'
zones do: a mesh of one block is one zone with every point, a mesh of several
gives each zone the points its cells use (a point several blocks use is
repeated, one no cell uses is dropped), as TecIO's own files have it. A series
is a strand of zones per step, the later steps sharing the first step's
coordinates and connectivity, which they must therefore have.

SZL is TecIO's to write: ``MESHIOPLUSPLUS_TECIO_LIBRARY`` names the shared
library, as for reading.
"""

from __future__ import annotations

import ctypes

import numpy as np

from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from ..tecplot._tecplot import (
    _zone_title,
    _zone_type,
    meshio_to_tecplot_order,
    write_plan,
)
from ._szplt import library

__all__ = ["write", "SeriesWriter"]

_ZONE_CODES = {
    "FELINESEG": 1,
    "FETRIANGLE": 2,
    "FEQUADRILATERAL": 3,
    "FETETRAHEDRON": 4,
    "FEBRICK": 5,
}
_DOUBLE = 2
# TecIO's variable types: float, double, int32, int16, byte
_FLOAT, _INT32, _INT16, _BYTE = 1, 3, 4, 5
_TYPE_OF_KIND = {
    "float32": _FLOAT,
    "int32": _INT32,
    "uint16": _INT32,
    "int16": _INT16,
    "int8": _INT16,
    "uint8": _BYTE,
    "bool": _BYTE,
}
_WRITE = {
    _FLOAT: ("tecZoneVarWriteFloatValues", np.float32, ctypes.c_float),
    _DOUBLE: ("tecZoneVarWriteDoubleValues", np.float64, ctypes.c_double),
    _INT32: ("tecZoneVarWriteInt32Values", np.int32, ctypes.c_int32),
    _INT16: ("tecZoneVarWriteInt16Values", np.int16, ctypes.c_int16),
    _BYTE: ("tecZoneVarWriteUInt8Values", np.uint8, ctypes.c_uint8),
}


def _check(status, what):
    if status != 0:
        raise WriteError(f"Tecplot .szplt: TecIO failed in {what}")


def _type_of(dtypes):
    """TecIO's type for a variable whose arrays have these dtypes: the one type that holds
    them all exactly, else a double."""
    types = {_TYPE_OF_KIND.get(np.dtype(d).name, _DOUBLE) for d in dtypes}
    return types.pop() if len(types) == 1 else _DOUBLE


def _pointer(array, ctype):
    """A ctypes pointer to an array's data, or None for an empty one."""
    return array.ctypes.data_as(ctypes.POINTER(ctype)) if len(array) else None


def _ints(values):
    return (ctypes.c_int32 * max(len(values), 1))(*[int(v) for v in values])


class _Writer:
    """One SZL file: the zones of every step added, written as they come."""

    def __init__(self, filename):
        self.lib = library()
        self.path = str(filename)
        self.handle = ctypes.c_void_p()
        self.layout = None  # the first step's variables and block types
        self.first = None  # its coordinates and connectivity
        self.first_step_zone = {}  # block position -> zone number (step 0)
        self.steps = 0

    def _open(self, variables):
        names = ",".join(v.replace(",", "_") for v in variables)
        _check(
            self.lib.tecFileWriterOpen(
                self.path.encode(),
                b"meshio++",
                names.encode(),
                ctypes.c_int32(1),
                ctypes.c_int32(0),
                ctypes.c_int32(_DOUBLE),
                None,
                ctypes.byref(self.handle),
            ),
            "tecFileWriterOpen",
        )

    def add_step(self, mesh, time=None):
        cell_blocks, variables, shared_data, cell_var_blocks, bases = write_plan(mesh)
        usable = []
        for ic in cell_blocks:
            if _zone_type(mesh.cells[ic].type) in _ZONE_CODES:
                usable.append(ic)
            else:
                warn(
                    f"Tecplot .szplt: skipping '{mesh.cells[ic].type}' cells (TecIO "
                    "2018.3 refuses polygon and polyhedral zones in SZL files)"
                )
        if not usable:
            raise WriteError("Tecplot .szplt: no cell type SZL can hold in the mesh")
        signature = (
            len(mesh.points),
            variables,
            [(mesh.cells[ic].type, len(mesh.cells[ic].data)) for ic in usable],
        )
        conns = [np.asarray(mesh.cells[ic].data, dtype=np.int64) for ic in usable]
        coords = np.asarray(mesh.points, dtype=np.float64)
        if self.layout is None:
            self.layout = signature
            self.first = (coords, conns)
            self._open(variables)
        elif signature != self.layout:
            raise WriteError(
                f"Tecplot .szplt: series: step {self.steps}'s points, cells or "
                "variables differ from the first step's; write one file per step"
            )
        elif self.steps > 0 and not (
            np.array_equal(coords, self.first[0])
            and all(np.array_equal(a, b) for a, b in zip(conns, self.first[1]))
        ):
            raise WriteError(
                f"Tecplot .szplt: series: step {self.steps}'s coordinates or "
                "connectivity differ from the first step's, which the later steps "
                "share; write one file per step"
            )
        nvar = len(variables)
        num_shared = len(shared_data)
        types = [_type_of([np.asarray(col).dtype]) for col in shared_data]
        for blocks in cell_var_blocks:
            types.append(_type_of([np.asarray(a).dtype for a in blocks.values()]))
        ncoord = 3 if mesh.points.shape[1] == 3 else 2
        compact = len(usable) > 1
        for bi, ic in enumerate(usable):
            cell = mesh.cells[ic]
            ztype = _ZONE_CODES[_zone_type(cell.type)]
            data = conns[bi]
            if compact:
                used, inverse = np.unique(data, return_inverse=True)
                local = inverse.reshape(data.shape)
            else:
                used, local = np.arange(len(mesh.points)), data
            present = [j for j, pb in enumerate(cell_var_blocks) if ic in pb]
            share = [0] * nvar
            location = [1 if v < num_shared else 0 for v in range(nvar)]
            passive = [0] * nvar
            for j in range(len(cell_var_blocks)):
                if j not in present:
                    passive[num_shared + j] = 1
            parts = _parts(mesh, ic) if time is None and ztype in (4, 5) else None
            plan = None
            if parts is not None and len(np.unique(parts)) > 1:
                plan = _partition_plan(parts, len(used), local)
                if plan is None:
                    warn(
                        "Tecplot .szplt: a part owns no node of its own, which TecIO "
                        "cannot partition; the zone is written whole"
                    )
            conn_share = self.first_step_zone.get(bi, 0) if self.steps > 0 else 0
            if conn_share:
                for v in range(ncoord):
                    share[v] = conn_share
            zone = ctypes.c_int32()
            _check(
                self.lib.tecZoneCreateFE(
                    self.handle,
                    _zone_title(mesh, bases, ic).encode(),
                    ctypes.c_int32(ztype),
                    ctypes.c_int64(len(used)),
                    ctypes.c_int64(len(data)),
                    _ints(types),
                    _ints(share),
                    _ints(location),
                    _ints(passive),
                    ctypes.c_int32(conn_share),
                    ctypes.c_int64(0),
                    ctypes.c_int32(0),
                    ctypes.byref(zone),
                ),
                "tecZoneCreateFE",
            )
            number = zone.value
            if self.steps == 0:
                self.first_step_zone[bi] = number
            if time is not None:
                _check(
                    self.lib.tecZoneSetUnsteadyOptions(
                        self.handle,
                        ctypes.c_int32(number),
                        ctypes.c_double(time),
                        ctypes.c_int32(1),
                    ),
                    "tecZoneSetUnsteadyOptions",
                )
            owned = []
            for v in range(num_shared):
                if share[v] == 0:
                    owned.append((v, np.asarray(shared_data[v])[used]))
            for j in present:
                owned.append((num_shared + j, np.asarray(cell_var_blocks[j][ic])))
            order = meshio_to_tecplot_order[cell.type]
            if plan is not None:
                self._write_partitions(
                    number, plan, owned, num_shared, types, local, order
                )
                continue
            for v, values in owned:
                self._write_values(number, v, 0, types[v], values)
            if conn_share == 0:
                nodes = np.ascontiguousarray(local[:, order] + 1, dtype=np.int32)
                _check(
                    self.lib.tecZoneNodeMapWrite32(
                        self.handle,
                        ctypes.c_int32(number),
                        ctypes.c_int32(0),
                        ctypes.c_int32(1),
                        ctypes.c_int64(nodes.size),
                        nodes.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                    ),
                    "tecZoneNodeMapWrite32",
                )
        self.steps += 1

    def _write_values(self, zone, var, partition, kind, values):
        """One variable's values in TecIO's type for it, to a zone's partition (0: the
        whole zone)."""
        name, dtype, ctype = _WRITE[kind]
        array = np.ascontiguousarray(values, dtype=dtype)
        if len(array):
            _check(
                getattr(self.lib, name)(
                    self.handle,
                    ctypes.c_int32(zone),
                    ctypes.c_int32(var + 1),
                    ctypes.c_int32(partition),
                    ctypes.c_int64(len(array)),
                    array.ctypes.data_as(ctypes.POINTER(ctype)),
                ),
                name,
            )

    def _write_partitions(self, zone, plan, owned, num_shared, types, local, order):
        """A 3D FE zone as one partition per part; see the module's notes."""
        cells, nodes, owner = plan
        # a partition lists the nodes it owns, then its ghost nodes, each in node order
        own = [n[owner[n] == p + 1] for p, n in enumerate(nodes)]
        ghost = [n[owner[n] != p + 1] for p, n in enumerate(nodes)]

        def locate(p, found):
            """The 1-based numbers of the nodes ``found`` in partition ``p``."""
            mine = owner[found] == p + 1
            return (
                np.where(
                    mine,
                    np.searchsorted(own[p], found),
                    len(own[p]) + np.searchsorted(ghost[p], found),
                )
                + 1
            )

        for p in range(len(nodes)):
            neighbor = owner[ghost[p]].astype(np.int32)
            neighbor_node = np.array(
                [
                    locate(o - 1, np.array([g]))[0]
                    for o, g in zip(neighbor.tolist(), ghost[p].tolist())
                ],
                dtype=np.int64,
            )
            ghost_local = (len(own[p]) + np.arange(1, len(ghost[p]) + 1)).astype(
                np.int64
            )
            _check(
                self.lib.tecFEPartitionCreate64(
                    self.handle,
                    ctypes.c_int32(zone),
                    ctypes.c_int32(p + 1),
                    ctypes.c_int64(len(nodes[p])),
                    ctypes.c_int64(len(cells[p])),
                    ctypes.c_int64(len(ghost[p])),
                    _pointer(ghost_local, ctypes.c_int64),
                    _pointer(neighbor, ctypes.c_int32),
                    _pointer(neighbor_node, ctypes.c_int64),
                    ctypes.c_int64(0),
                    None,
                ),
                "tecFEPartitionCreate64",
            )
        for p in range(len(nodes)):
            listed = np.concatenate([own[p], ghost[p]])
            for v, values in owned:
                if v < num_shared:
                    self._write_values(zone, v, p + 1, types[v], values[listed])
                else:
                    self._write_values(zone, v, p + 1, types[v], values[cells[p]])
            nodemap = np.ascontiguousarray(
                locate(p, local[cells[p]][:, order]), dtype=np.int32
            )
            _check(
                self.lib.tecZoneNodeMapWrite32(
                    self.handle,
                    ctypes.c_int32(zone),
                    ctypes.c_int32(p + 1),
                    ctypes.c_int32(1),
                    ctypes.c_int64(nodemap.size),
                    nodemap.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
                ),
                "tecZoneNodeMapWrite32",
            )

    def close(self):
        if self.handle:
            _check(
                self.lib.tecFileWriterClose(ctypes.byref(self.handle)),
                "tecFileWriterClose",
            )
            self.handle = ctypes.c_void_p()


def _partition_plan(parts, num_nodes, local):
    """How a 3D FE zone splits into one partition per part: ``(cells, nodes, owner)``, each
    partition's cells and nodes and each node's owning partition (from 1); None when TecIO
    could not take it (a partition owning no node)."""
    labels = np.unique(parts)
    slot = np.searchsorted(labels, parts)  # each cell's partition, from 0
    cells = [np.flatnonzero(slot == p) for p in range(len(labels))]
    nodes = [np.unique(local[c]) for c in cells]
    # a node belongs to the lowest partition that uses it (so a ghost node always names an
    # earlier partition); one no cell uses is the first partition's
    owner = np.zeros(num_nodes, dtype=np.int32)
    for p in range(len(nodes) - 1, -1, -1):
        owner[nodes[p]] = p + 1
    loose = np.flatnonzero(owner == 0)
    owner[loose] = 1
    nodes[0] = np.union1d(nodes[0], loose)
    if not np.all(np.bincount(owner, minlength=len(nodes) + 1)[1:]):
        return None
    return cells, nodes, owner


def _parts(mesh, block):
    """The ``partition:part`` values of a cell block, or None."""
    data = getattr(mesh, "cell_data", {}) or {}
    values = data.get("partition:part")
    if values is None or block >= len(values):
        return None
    a = np.asarray(values[block])
    if (
        a.ndim != 1
        or len(a) != len(mesh.cells[block].data)
        or a.dtype.kind not in "iub"
    ):
        return None
    return a


def _time_of(mesh):
    fd = getattr(mesh, "field_data", {}) or {}
    if "meshio:time" in fd:
        return float(np.asarray(fd["meshio:time"]).reshape(-1)[0])
    return None


def write(filename, mesh):
    """Write one ``.szplt`` file holding the mesh (a static data set)."""
    if is_buffer(filename, "w"):
        raise WriteError("Tecplot .szplt: cannot be written to a stream")
    writer = _Writer(filename)
    try:
        writer.add_step(mesh)
    finally:
        writer.close()


class SeriesWriter:
    """A strand of zones per step; the later steps share the first step's
    coordinates and connectivity, so every step needs the first's points and
    cells."""

    def __init__(self, filename):
        if is_buffer(filename, "w"):
            raise WriteError("Tecplot .szplt: cannot be written to a stream")
        self._writer = _Writer(filename)

    def write(self, time, mesh):
        self._writer.add_step(mesh, float(time))

    def close(self):
        if self._writer is not None:
            self._writer.close()
            self._writer = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
