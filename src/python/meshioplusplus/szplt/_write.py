"""Write Tecplot SZL (``.szplt``) files through a user-supplied shared TecIO.

The Python twin of ``write_szplt`` in ``src/cpp/src/formats/tecplot.cpp``: both
lay the zones out as the ASCII and ``.plt`` writers do (``tecplot._tecplot``'s
``write_plan``): one zone per cell block, a cell-centred variable passive in a
zone that has no values for it. Wedges and pyramids are degenerate bricks, as
there. Values are doubles.

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


def _check(status, what):
    if status != 0:
        raise WriteError(f"Tecplot .szplt: TecIO failed in {what}")


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
                    f"Tecplot .szplt: skipping '{mesh.cells[ic].type}' cells (SZL "
                    "files are written without face-based zones)"
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
                    _ints([_DOUBLE] * nvar),
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
                owned.append((num_shared + j, cell_var_blocks[j][ic]))
            for v, values in owned:
                array = np.ascontiguousarray(values, dtype=np.float64)
                if len(array):
                    _check(
                        self.lib.tecZoneVarWriteDoubleValues(
                            self.handle,
                            ctypes.c_int32(number),
                            ctypes.c_int32(v + 1),
                            ctypes.c_int32(0),
                            ctypes.c_int64(len(array)),
                            array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ),
                        "tecZoneVarWriteDoubleValues",
                    )
            if conn_share == 0:
                order = meshio_to_tecplot_order[cell.type]
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

    def close(self):
        if self.handle:
            _check(
                self.lib.tecFileWriterClose(ctypes.byref(self.handle)),
                "tecFileWriterClose",
            )
            self.handle = ctypes.c_void_p()


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
