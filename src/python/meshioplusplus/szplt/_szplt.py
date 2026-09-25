"""Tecplot SZL (``.szplt``) through a user-supplied shared TecIO, via ctypes.

SZL is undocumented; TecIO, Tecplot's own library, is its only reader, and
meshio++ never redistributes it. This reader loads a shared ``libtecio`` the
user already has -- ``MESHIOPLUSPLUS_TECIO_LIBRARY`` names it, else the system
library path is searched (Tecplot 360 ships ``libtecio.so``/``tecio.dll``) --
and fills the zone model the ASCII and ``.plt`` readers share, so a ``.szplt``
reads exactly as the same data saved as ``.plt``.

The Python twin of the TecIO source in ``src/cpp/src/formats/tecplot.cpp``;
see ``doc/formats/szplt.md``.
"""

import ctypes
import ctypes.util
import os

import numpy as np

from .._exceptions import ReadError
from ..tecplot import _zones
from ..tecplot._zones import Zone

_ZONE_TYPES = (
    "ORDERED",
    "FELINESEG",
    "FETRIANGLE",
    "FEQUADRILATERAL",
    "FETETRAHEDRON",
    "FEBRICK",
    "FEPOLYGON",
    "FEPOLYHEDRON",
)
_NODES_PER_CELL = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8}
# TecIO data types: 1 float, 2 double, 3 int32, 4 int16, 5 byte.
_GETTERS = {
    1: ("tecZoneVarGetFloatValues", ctypes.c_float),
    2: ("tecZoneVarGetDoubleValues", ctypes.c_double),
    3: ("tecZoneVarGetInt32Values", ctypes.c_int32),
    4: ("tecZoneVarGetInt16Values", ctypes.c_int16),
    5: ("tecZoneVarGetUInt8Values", ctypes.c_uint8),
}

_LIB = None


def library():
    """The shared TecIO, loaded once; ``ImportError`` naming the remedies."""
    global _LIB
    if _LIB is not None:
        return _LIB
    path = os.environ.get("MESHIOPLUSPLUS_TECIO_LIBRARY") or ctypes.util.find_library(
        "tecio"
    )
    if not path:
        raise ImportError(
            "meshio++: reading a Tecplot .szplt file needs TecIO, Tecplot's own "
            "library: set MESHIOPLUSPLUS_TECIO_LIBRARY to a shared libtecio (Tecplot "
            "360 ships one), build meshio++ with MESHIOPLUSPLUS_WITH_TECIO=ON, or "
            "save the file as .plt in Tecplot"
        )
    try:
        lib = ctypes.CDLL(path)
    except OSError as e:
        raise ImportError(f"meshio++: cannot load TecIO from '{path}': {e}") from e
    _LIB = lib
    return lib


def _check(status, what):
    if status != 0:
        raise ReadError(f"Tecplot .szplt: TecIO failed in {what}")


class _Handle:
    def __init__(self, filename):
        self.lib = library()
        self.h = ctypes.c_void_p()
        if (
            self.lib.tecFileReaderOpen(str(filename).encode(), ctypes.byref(self.h))
            != 0
            or not self.h
        ):
            raise ReadError(f"Tecplot .szplt: TecIO cannot open '{filename}'")

    def __del__(self):
        if getattr(self, "h", None):
            self.lib.tecFileReaderClose(ctypes.byref(self.h))

    def int32(self, fn, *args):
        out = ctypes.c_int32()
        _check(getattr(self.lib, fn)(self.h, *args, ctypes.byref(out)), fn)
        return out.value

    def int64(self, fn, *args):
        out = ctypes.c_int64()
        _check(getattr(self.lib, fn)(self.h, *args, ctypes.byref(out)), fn)
        return out.value

    def string(self, fn, *args):
        out = ctypes.c_char_p()
        _check(getattr(self.lib, fn)(self.h, *args, ctypes.byref(out)), fn)
        text = out.value.decode("utf-8", "replace") if out.value else ""
        if out:
            self.lib.tecStringFree(ctypes.byref(out))
        return text


class SzlSource:
    def __init__(self, handle, zones, num_variables):
        self.handle = handle
        self.zones = zones
        self.num_variables = num_variables

    def own_data(self, idx):
        z = self.zones[idx]
        h = self.handle
        zone = ctypes.c_int32(idx + 1)
        cols = {}
        for v in range(self.num_variables):
            if not z.owns(v):
                continue
            var = ctypes.c_int32(v + 1)
            n = h.int64("tecZoneVarGetNumValues", zone, var)
            dtype = h.int32("tecZoneVarGetType", zone, var)
            if dtype not in _GETTERS:
                raise ReadError(
                    f"Tecplot .szplt: variable {v + 1} has the unknown data type {dtype}"
                )
            fn, ctype = _GETTERS[dtype]
            buf = (ctype * max(n, 1))()
            if n > 0:
                _check(
                    getattr(h.lib, fn)(
                        h.h, zone, var, ctypes.c_int64(1), ctypes.c_int64(n), buf
                    ),
                    fn,
                )
            values = np.ctypeslib.as_array(buf)[:n].astype(np.float64)
            if len(values) != z.data_length(v):
                # An ordered zone's cell-centred values over the full I x J x K
                # index space: keep i < I-1, j < J-1, k < K-1.
                i, j, k = z.ijk
                if not (z.ordered and z.cell_centered[v] and len(values) == i * j * k):
                    raise ReadError(
                        f"Tecplot .szplt: zone {idx + 1} variable {v + 1} holds "
                        f"{len(values)} values"
                    )
                values = values.reshape(k, j, i)[
                    : max(k - 1, 1), : max(j - 1, 1), : max(i - 1, 1)
                ].ravel()
            cols[v] = values
        conn = None
        if not z.ordered and z.conn_share < 0:
            npc = _NODES_PER_CELL[_ZONE_TYPES.index(z.type_name)]
            count = h.int64(
                "tecZoneNodeMapGetNumValues", zone, ctypes.c_int64(z.num_cells)
            )
            if count != z.num_cells * npc:
                raise ReadError(
                    f"Tecplot .szplt: zone {idx + 1} node map holds {count} values, "
                    f"not {z.num_cells * npc}"
                )
            is64 = h.int32("tecZoneNodeMapIs64Bit", zone)
            fn, ctype = (
                ("tecZoneNodeMapGet64", ctypes.c_int64)
                if is64
                else ("tecZoneNodeMapGet", ctypes.c_int32)
            )
            buf = (ctype * max(count, 1))()
            _check(
                getattr(h.lib, fn)(
                    h.h, zone, ctypes.c_int64(1), ctypes.c_int64(z.num_cells), buf
                ),
                fn,
            )
            # TecIO hands the node map back 1-based.
            conn = np.ctypeslib.as_array(buf)[:count].astype(np.int64).reshape(-1, npc)
            conn -= 1
        return cols, conn


def load(filename):
    """Reads a ``.szplt`` header through TecIO: (variables, zones, source)."""
    h = _Handle(filename)
    nvar = h.int32("tecDataSetGetNumVars")
    nzone = h.int32("tecDataSetGetNumZones")
    variables = [
        h.string("tecVarGetName", ctypes.c_int32(v)) for v in range(1, nvar + 1)
    ]
    if nzone < 1:
        raise ReadError("Tecplot .szplt: no zone")
    zones = []
    for zi in range(1, nzone + 1):
        zone = ctypes.c_int32(zi)
        ztype = h.int32("tecZoneGetType", zone)
        if not 0 <= ztype <= 7:
            raise ReadError(f"Tecplot .szplt: zone {zi} has type {ztype}")
        z = Zone()
        z.type_name = _ZONE_TYPES[ztype]
        if z.is_poly:
            raise ReadError(
                f"Tecplot .szplt: {z.type_name} zones are not read; save the file as "
                ".plt in Tecplot"
            )
        z.ordered = ztype == 0
        z.title = h.string("tecZoneGetTitle", zone)
        ijk = [ctypes.c_int64() for _ in range(3)]
        _check(
            h.lib.tecZoneGetIJK(h.h, zone, *(ctypes.byref(x) for x in ijk)),
            "tecZoneGetIJK",
        )
        i, j, k = (x.value for x in ijk)
        if z.ordered:
            z.ijk = (max(i, 1), max(j, 1), max(k, 1))
            z.finish()
        else:
            # An FE zone's I is its node count, J its cell count.
            z.num_nodes, z.num_cells = i, j
        z.cell_centered = []
        for v in range(nvar):
            var = ctypes.c_int32(v + 1)
            location = h.int32("tecZoneVarGetValueLocation", zone, var)
            passive = h.int32("tecZoneVarIsPassive", zone, var)
            shared = h.int32("tecZoneVarGetSharedZone", zone, var)
            z.cell_centered.append(location == 0)  # 0 cell-centred, 1 nodal
            if passive:
                z.passive.add(v)
            elif shared > 0:
                z.var_share[v] = shared - 1
        conn_share = h.int32("tecZoneConnectivityGetSharedZone", zone)
        z.conn_share = conn_share - 1 if conn_share > 0 else -1
        time = ctypes.c_double()
        _check(
            h.lib.tecZoneGetSolutionTime(h.h, zone, ctypes.byref(time)),
            "tecZoneGetSolutionTime",
        )
        strand = h.int32("tecZoneGetStrandID", zone)
        # TecIO's strands are 1-based (0 static), like the ASCII STRANDID.
        z.has_solution_time = strand != 0 or time.value != 0.0
        z.solution_time = time.value
        z.has_strand = strand > 0
        z.strand = strand
        zones.append(z)
    return variables, zones, SzlSource(h, zones, nvar)


def read(filename, time_step=0):
    variables, zones, source = load(filename)
    steps = _zones.timeline(zones)
    count = len(steps)
    step = time_step + count if time_step < 0 else time_step
    if not 0 <= step < count:
        raise ReadError(
            f"meshio++: time step {time_step} is out of range: this file has {count} "
            + ("step" if count == 1 else "steps")
        )
    mesh = _zones.build_step(steps[step], zones, variables, source)
    # read_metadata's side channel.
    mesh.time_values = _zones.metadata(zones, variables, source)[2]
    return mesh


def time_values(filename):
    variables, zones, source = load(filename)
    return _zones.metadata(zones, variables, source)[2]
