"""
Binary Tecplot ``.plt`` reader: the header and data sections of the Data
Format Guide's appendix A ("Binary Data File Format", version ``#!TDV112``),
decoded into the same zone model as the ASCII reader (``_zones.py``) so both
build one mesh the same way.
"""

import numpy as np

from .._exceptions import ReadError
from .._files import open_file
from ._zones import Zone

ZONE_MARKER = 299.0
GEOMETRY_MARKER = 399.0
TEXT_MARKER = 499.0
LABEL_MARKER = 599.0
USERREC_MARKER = 699.0
DATASETAUX_MARKER = 799.0
VARAUX_MARKER = 899.0
EOH_MARKER = 357.0

ZONE_TYPES = {
    0: "ORDERED",
    1: "FELINESEG",
    2: "FETRIANGLE",
    3: "FEQUADRILATERAL",
    4: "FETETRAHEDRON",
    5: "FEBRICK",
    6: "FEPOLYGON",
    7: "FEPOLYHEDRON",
}
NODES_PER_CELL = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8}
FACES_PER_CELL = {1: 0, 2: 3, 3: 4, 4: 4, 5: 6}
# Variable data formats: 1 float, 2 double, 3 int32, 4 int16, 5 byte, 6 bit.
FORMATS = {1: "f4", 2: "f8", 3: "i4", 4: "i2", 5: "u1"}


def is_plt(head):
    return head[:5] == b"#!TDV"


class _Cursor:
    def __init__(self, data, filename):
        self.data = data
        self.pos = 0
        self.order = "<"
        self.filename = filename

    def need(self, n):
        if self.pos + n > len(self.data):
            raise ReadError(f"Tecplot .plt: {self.filename} is truncated")

    def array(self, fmt, n):
        dt = np.dtype(self.order + fmt)
        self.need(dt.itemsize * n)
        out = np.frombuffer(self.data, dtype=dt, count=n, offset=self.pos)
        self.pos += dt.itemsize * n
        return out

    def i32(self):
        return int(self.array("i4", 1)[0])

    def f32(self):
        return float(self.array("f4", 1)[0])

    def f64(self):
        return float(self.array("f8", 1)[0])

    def string(self):
        chars = []
        while True:
            c = self.i32()
            if c == 0:
                return "".join(chars)
            chars.append(chr(c))


def _skip_geometry(cur):
    coord_sys = cur.i32()  # 4 = Grid3D: polylines carry Z too
    cur.i32()  # scope
    cur.i32()  # draw order
    cur.array("f8", 3)  # anchor (the layout TecIO's own reader expects for V112)
    cur.array("i4", 4)  # zone, color, fill color, is filled
    gtype = cur.i32()
    cur.i32()  # line pattern
    cur.array("f8", 2)  # pattern length, line thickness
    cur.array("i4", 3)  # ellipse points, arrowhead style, attachment
    cur.array("f8", 2)  # arrowhead size, angle
    cur.string()  # macro function command
    field = cur.i32()  # 1 float, 2 double
    cur.i32()  # clipping
    fmt = "f8" if field == 2 else "f4"
    if gtype == 0:  # polylines: per line, its point count then X, Y (, Z) blocks
        for _ in range(cur.i32()):
            n = cur.i32()
            cur.array(fmt, (3 if coord_sys == 4 else 2) * n)
    elif gtype in (1, 4):  # rectangle, ellipse
        cur.array(fmt, 2)
    elif gtype in (2, 3):  # square, circle
        cur.array(fmt, 1)
    else:
        raise ReadError(f"Tecplot .plt: unknown geometry type {gtype}")


def _skip_text(cur):
    cur.i32()  # position coordinate system
    cur.i32()  # scope
    cur.array("f8", 3)
    cur.i32()  # font
    cur.i32()  # height units
    cur.f64()
    cur.i32()  # box type
    cur.array("f8", 2)
    cur.array("i4", 2)
    cur.array("f8", 2)  # angle, line spacing
    cur.i32()  # anchor
    cur.i32()  # zone
    cur.i32()  # color
    cur.string()  # macro function command
    cur.i32()  # clipping
    cur.string()


def _read_header(cur):
    magic = bytes(cur.data[:8])
    if not is_plt(magic):
        raise ReadError("Tecplot .plt: missing #!TDV magic")
    try:
        version = int(magic[5:8].decode("ascii").strip())
    except ValueError:
        version = 0
    if version != 112:
        raise ReadError(
            f"Tecplot .plt: version {magic[5:8].decode('ascii', 'replace')} is not supported "
            "(only #!TDV112, written by Tecplot 360 2009 and later, is read)"
        )
    cur.pos = 8
    if cur.i32() != 1:
        cur.order = ">"
        cur.pos = 8
        if cur.i32() != 1:
            raise ReadError("Tecplot .plt: bad byte-order word")
    file_type = cur.i32()
    title = cur.string()
    nvar = cur.i32()
    variables = [cur.string() for _ in range(nvar)]

    zones = []
    while True:
        marker = cur.f32()
        if marker == ZONE_MARKER:
            z = Zone()
            z.title = cur.string()
            cur.i32()  # parent zone
            strand = cur.i32()
            time = cur.f64()
            cur.i32()  # not used (-1)
            ztype = cur.i32()
            if ztype not in ZONE_TYPES:
                raise ReadError(f"Tecplot .plt: unknown zone type {ztype}")
            z.type_name = ZONE_TYPES[ztype]
            z.type_code = ztype
            z.cell_centered = [0] * nvar
            if cur.i32() == 1:
                z.cell_centered = [int(v) for v in cur.array("i4", nvar)]
            z.raw_face_neighbors = cur.i32()
            z.misc_face_neighbors = cur.i32()
            z.face_neighbor_mode = 0
            if z.misc_face_neighbors:
                z.face_neighbor_mode = cur.i32()
                if ztype != 0:
                    cur.i32()  # FE face neighbours completely specified
            if ztype == 0:
                z.ordered = True
                z.ijk = tuple(int(v) for v in cur.array("i4", 3))
            else:
                z.num_nodes = cur.i32()
                if ztype in (6, 7):
                    cur.array("i4", 4)
                z.num_cells = cur.i32()
                cur.array("i4", 3)  # I/J/K cell dims, unused
            while cur.i32() == 1:  # auxiliary name/value pairs
                cur.string()
                cur.i32()
                cur.string()
            # File strands are 0-based (-1 = static); the ASCII STRANDID is
            # 1-based, but only the grouping by strand matters.
            z.has_solution_time = strand != -1 or time != 0.0
            z.solution_time = time
            z.has_strand = strand >= 0
            z.strand = strand
            zones.append(z)
        elif marker == GEOMETRY_MARKER:
            _skip_geometry(cur)
        elif marker == TEXT_MARKER:
            _skip_text(cur)
        elif marker == LABEL_MARKER:
            for _ in range(cur.i32()):
                cur.string()
        elif marker == USERREC_MARKER:
            cur.string()
        elif marker == DATASETAUX_MARKER:
            cur.string()
            cur.i32()
            cur.string()
        elif marker == VARAUX_MARKER:
            cur.i32()
            cur.string()
            cur.i32()
            cur.string()
        elif marker == EOH_MARKER:
            break
        else:
            raise ReadError(
                f"Tecplot .plt: unexpected header marker {marker} at byte {cur.pos - 4}"
            )
    if not zones:
        raise ReadError("Tecplot: no ZONE")
    return file_type, title, variables, zones


def _ordered_cc_layout(ijk):
    """How a cell-centred variable of an ordered zone is stored: the node
    dimensions with the last one longer than 1 shortened by one (appendix A,
    note 5), so every direction but that last carries a ghost value at its end.
    Returns (stored shape as (K, J, I), the slices that keep the real cells)."""
    dims = list(ijk)
    long = [a for a in range(3) if dims[a] > 1]
    stored = dims[:]
    if long:
        stored[long[-1]] -= 1
    keep = [slice(0, d - 1 if d > 1 else 1) for d in dims]
    return (stored[2], stored[1], stored[0]), (keep[2], keep[1], keep[0])


def _skip_face_neighbors(cur, zone):
    """Skips the user-defined face-neighbour connections (appendix A, note 4)."""
    mode = zone.face_neighbor_mode
    for _ in range(zone.misc_face_neighbors):
        if mode == 0:  # local one-to-one: cz, fz, cz
            cur.array("i4", 3)
        elif mode == 2:  # global one-to-one: cz, fz, ZZ, CZ
            cur.array("i4", 4)
        elif mode in (1, 3):  # one-to-many: cz, fz, oz, nz, then nz (or 2*nz) values
            nz = int(cur.array("i4", 4)[3])
            cur.array("i4", nz if mode == 1 else 2 * nz)
        else:
            raise ReadError(f"Tecplot .plt: unknown face-neighbour mode {mode}")


def _scan_data(cur, variables, zones):
    """Walks every zone's data section, recording where each owned variable
    and the connectivity start, and the sharing that lives here (not in the
    header) in a binary file."""
    nvar = len(variables)
    for index, z in enumerate(zones):
        marker = cur.f32()
        if marker != ZONE_MARKER:
            raise ReadError(
                f"Tecplot .plt: expected the zone {index + 1} data marker, got {marker}"
            )
        formats = [int(f) for f in cur.array("i4", nvar)]
        if cur.i32():
            z.passive = {v for v, p in enumerate(cur.array("i4", nvar)) if p}
        if cur.i32():
            z.var_share = {
                v: int(s) for v, s in enumerate(cur.array("i4", nvar)) if s >= 0
            }
        z.conn_share = cur.i32()
        owned = [v for v in range(nvar) if z.owns(v)]
        cur.array("f8", 2 * len(owned))  # min/max pairs
        z.var_offsets = {}
        z.var_formats = {}
        for v in owned:
            if formats[v] == 6:
                raise ReadError(
                    f"Tecplot .plt: variable '{variables[v]}' is BIT-packed, which is not supported"
                )
            if formats[v] not in FORMATS:
                raise ReadError(f"Tecplot .plt: unknown data format {formats[v]}")
            if z.ordered and z.cell_centered[v]:
                shape, _ = _ordered_cc_layout(z.ijk)
                n = int(np.prod(shape))
            else:
                n = z.data_length(v)
            z.var_offsets[v] = (cur.pos, n)
            z.var_formats[v] = FORMATS[formats[v]]
            cur.array(FORMATS[formats[v]], n)
        z.conn_offset = None
        if z.ordered:
            if z.conn_share < 0 and z.misc_face_neighbors:
                _skip_face_neighbors(cur, z)
            continue
        if z.type_code in (6, 7):
            raise ReadError(
                f"Tecplot: {z.type_name} zones (polygonal/polyhedral) are not supported"
            )
        if z.conn_share < 0:
            z.conn_offset = cur.pos
            cur.array("i4", NODES_PER_CELL[z.type_code] * z.num_cells)
            if z.raw_face_neighbors:
                cur.array("i4", FACES_PER_CELL[z.type_code] * z.num_cells)
            if z.misc_face_neighbors:
                _skip_face_neighbors(cur, z)


class PltSource:
    def __init__(self, cur, zones):
        self.cur = cur
        self.zones = zones

    def own_data(self, idx):
        z = self.zones[idx]
        cur = self.cur
        cols = {}
        for v, (offset, n) in z.var_offsets.items():
            cur.pos = offset
            values = cur.array(z.var_formats[v], n).astype(np.float64)
            if z.ordered and z.cell_centered[v]:
                shape, keep = _ordered_cc_layout(z.ijk)
                values = values.reshape(shape)[keep].ravel()
            cols[v] = values
        conn = None
        if z.conn_offset is not None:
            cur.pos = z.conn_offset
            npc = NODES_PER_CELL[z.type_code]
            conn = cur.array("i4", npc * z.num_cells).astype(np.int64).reshape(-1, npc)
        return cols, conn


def load(filename):
    """Parses a ``.plt`` file: (variables, zones, source)."""
    with open_file(filename, "rb") as f:
        data = f.read()
    cur = _Cursor(data, filename)
    _, _, variables, zones = _read_header(cur)
    for z in zones:
        z.finish()
    _scan_data(cur, variables, zones)
    return variables, zones, PltSource(cur, zones)
