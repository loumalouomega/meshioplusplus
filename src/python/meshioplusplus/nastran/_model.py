"""The Nastran model half shared by the result readers (``nastran_h5``,
``nastran_op2``): element cards to cell blocks, element ids to cells, property
ids to regions. The Python twin of ``src/cpp/src/detail/nastran_model.cpp``.
"""

import math

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._mesh import CellBlock
from .._regions import Region

__all__ = [
    "CARDS",
    "CoordCard",
    "CoordSystems",
    "add_cells",
    "apply_frames",
    "rotate_to_basic",
]

# Nastran numbers the hex20/wedge15 mid-side nodes bottom, vertical, top;
# meshio++ (VTK) numbers them bottom, top, vertical: conn[k] = G[perm[k]].
_HEXA20 = list(range(12)) + [16, 17, 18, 19, 12, 13, 14, 15]
_PENTA15 = list(range(9)) + [12, 13, 14, 9, 10, 11]

# card -> (linear type, nodes, quadratic type or None, nodes, permutation or None)
CARDS = {
    "CBAR": ("line", 2, None, 0, None),
    "CBEAM": ("line", 2, None, 0, None),
    "CBUSH": ("line", 2, None, 0, None),
    "CHEXA": ("hexahedron", 8, "hexahedron20", 20, _HEXA20),
    "CONM2": ("vertex", 1, None, 0, None),
    "CONROD": ("line", 2, None, 0, None),
    "CPENTA": ("wedge", 6, "wedge15", 15, _PENTA15),
    "CPYRAM": ("pyramid", 5, "pyramid13", 13, None),
    "CQUAD": ("quad", 4, "quad9", 9, None),
    "CQUAD4": ("quad", 4, None, 0, None),
    "CQUAD8": ("quad", 4, "quad8", 8, None),
    "CQUADR": ("quad", 4, None, 0, None),
    "CROD": ("line", 2, None, 0, None),
    "CSHEAR": ("quad", 4, None, 0, None),
    "CTETRA": ("tetra", 4, "tetra10", 10, None),
    "CTRIA3": ("triangle", 3, None, 0, None),
    "CTRIA6": ("triangle", 3, "triangle6", 6, None),
    "CTRIAR": ("triangle", 3, None, 0, None),
    "CTUBE": ("line", 2, None, 0, None),
    "CVISC": ("line", 2, None, 0, None),
    "PLOTEL": ("line", 2, None, 0, None),
}

_DIM = {
    "vertex": 0,
    "line": 1,
    "triangle": 2,
    "triangle6": 2,
    "quad": 2,
    "quad8": 2,
    "quad9": 2,
    "tetra": 3,
    "tetra10": 3,
    "pyramid": 3,
    "pyramid13": 3,
    "wedge": 3,
    "wedge15": 3,
    "hexahedron": 3,
    "hexahedron20": 3,
}


def add_cells(cards, grid_index, scalar_points, ptype, who):
    """Cell blocks, cell data and property regions of a Nastran model.

    ``cards`` is a list of ``(card, eid, pid, g)`` with ``g`` an ``(n, width)``
    array of grid ids, in the order the blocks should take (cards without a
    spec in :data:`CARDS` are ignored). Returns ``(cells, cell_data, regions,
    cell_index, offsets, sizes)``: ``cell_index`` maps EID -> global cell
    (without CONM2, which MSC accepts sharing an id with a structural element).
    """
    blocks = []  # [type, conn list, eid list, pid list, card]
    dropped = 0
    for card, eid, pid, g in cards:
        spec = CARDS.get(card)
        if spec is None:
            continue
        lin_type, lin_n, quad_type, quad_n, perm = spec
        g = np.asarray(g, dtype=np.int64)
        width = g.shape[1] if g.ndim == 2 else 0
        linear = [lin_type, [], [], [], card]
        quadratic = [quad_type, [], [], [], card]
        partial = 0
        eid = [int(v) for v in eid]
        pid = [int(v) for v in pid]
        for i, row in enumerate(g.tolist()):
            quad = False
            if quad_type is not None and width >= quad_n:
                given = sum(1 for k in range(lin_n, quad_n) if row[k] != 0)
                quad = given == quad_n - lin_n
                if given and not quad:
                    partial += 1
            nodes = quad_n if quad else lin_n
            src = perm if (quad and perm) else range(nodes)
            conn = []
            for k in src:
                p = grid_index.get(row[k])
                if p is None:
                    if row[k] not in scalar_points:
                        what = "no node" if row[k] == 0 else f"undefined GRID {row[k]}"
                        raise ReadError(
                            f"{who}: {card} {eid[i]} references {what} as its node {k + 1}"
                        )
                    conn = None  # a scalar point, not a GRID
                    break
                conn.append(p)
            if conn is None:
                dropped += 1
                continue
            b = quadratic if quad else linear
            b[1].append(conn)
            b[2].append(eid[i])
            b[3].append(pid[i])
        if partial:
            warn(
                f"{who}: {partial} {card} element(s) have only some mid-side nodes; read as {lin_type}"
            )
        for b in (linear, quadratic):
            if b[2]:
                blocks.append(b)
    if dropped:
        warn(f"{who}: skipped {dropped} element(s) that connect scalar points")

    cells = []
    cell_data = {}
    cell_index = {}
    offsets = []
    sizes = []
    ncells = 0
    shared = 0
    for ctype, conn, eid, pid, card in blocks:
        offsets.append(ncells)
        sizes.append(len(eid))
        if card != "CONM2":
            for i, e in enumerate(eid):
                if e in cell_index:
                    shared += 1
                else:
                    cell_index[e] = ncells + i
        ncells += len(eid)
        cells.append(
            CellBlock(ctype, np.asarray(conn, dtype=np.int64).reshape(len(eid), -1))
        )
    if shared:
        warn(
            f"{who}: {shared} element id(s) are used by more than one card; "
            "their element results go to the first"
        )
    if blocks:
        cell_data["nastran:eid"] = [np.asarray(b[2], dtype=np.int64) for b in blocks]
        cell_data["nastran:pid"] = [np.asarray(b[3], dtype=np.int64) for b in blocks]

    by_pid = {}
    for b, (ctype, _, _, pid, _) in enumerate(blocks):
        dim = _DIM[ctype]
        for i, p in enumerate(pid):
            if p <= 0:
                continue
            entry = by_pid.setdefault(p, [[], dim])
            entry[0].append(offsets[b] + i)
            entry[1] = max(entry[1], dim)
    regions = [
        Region(
            f"{ptype.get(p, 'PID')}_{p}",
            "cell",
            np.asarray(entries, dtype=np.int64),
            dim,
            p,
        )
        for p, (entries, dim) in sorted(by_pid.items())
    ]
    regions.sort(key=lambda r: r.name)
    return cells, cell_data, regions, cell_index, offsets, sizes


# Degrees to radians; one constant so both engines round the same way.
_DEGREE = 3.14159265358979323846 / 180.0


class CoordCard:
    """One CORD1R/C/S (``grids``) or CORD2R/C/S (``rid`` and ``abc``) card."""

    def __init__(self, cid, ctype, rid=0, abc=None, grids=None):
        self.cid = int(cid)
        self.type = int(ctype)  # 1 rectangular, 2 cylindrical, 3 spherical
        self.rid = int(rid)
        self.abc = [float(v) for v in abc] if abc is not None else None
        self.grids = [int(g) for g in grids] if grids is not None else None


def _local_to_cartesian(ctype, p):
    if ctype == 2:
        t = p[1] * _DEGREE
        return [p[0] * math.cos(t), p[0] * math.sin(t), p[2]]
    if ctype == 3:
        t = p[1] * _DEGREE
        f = p[2] * _DEGREE
        return [
            p[0] * math.sin(t) * math.cos(f),
            p[0] * math.sin(t) * math.sin(f),
            p[0] * math.cos(t),
        ]
    return [p[0], p[1], p[2]]


def _system_from_points(ctype, a, b, c):
    z = [b[0] - a[0], b[1] - a[1], b[2] - a[2]]
    nz = math.sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2])
    if not nz > 0.0:
        return None
    ez = [z[0] / nz, z[1] / nz, z[2] / nz]
    v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]]
    y = [
        ez[1] * v[2] - ez[2] * v[1],
        ez[2] * v[0] - ez[0] * v[2],
        ez[0] * v[1] - ez[1] * v[0],
    ]
    ny = math.sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2])
    if not ny > 0.0:
        return None
    ey = [y[0] / ny, y[1] / ny, y[2] / ny]
    ex = [
        ey[1] * ez[2] - ey[2] * ez[1],
        ey[2] * ez[0] - ey[0] * ez[2],
        ey[0] * ez[1] - ey[1] * ez[0],
    ]
    return (ctype, list(a), [ex, ey, ez])


class CoordSystems:
    """Resolved systems: type, origin and axes (rows x, y, z) in basic."""

    def __init__(self):
        self.systems = {
            0: (1, [0.0, 0.0, 0.0], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        }

    def has(self, cid):
        return cid in self.systems

    def to_basic(self, cid, local):
        ctype, o, ax = self.systems[cid]
        c = _local_to_cartesian(ctype, local)
        return [
            o[k] + ax[0][k] * c[0] + ax[1][k] * c[1] + ax[2][k] * c[2] for k in range(3)
        ]

    def vector_to_basic(self, cid, point, v):
        ctype, o, ax = self.systems[cid]
        w = [v[0], v[1], v[2]]
        if ctype in (2, 3):
            d = [point[0] - o[0], point[1] - o[1], point[2] - o[2]]
            L = [ax[j][0] * d[0] + ax[j][1] * d[1] + ax[j][2] * d[2] for j in range(3)]
            ph = math.atan2(L[1], L[0])
            cp = math.cos(ph)
            sp = math.sin(ph)
            if ctype == 2:
                w[0] = v[0] * cp - v[1] * sp
                w[1] = v[0] * sp + v[1] * cp
            else:
                th = math.atan2(math.sqrt(L[0] * L[0] + L[1] * L[1]), L[2])
                ct = math.cos(th)
                st = math.sin(th)
                er = [st * cp, st * sp, ct]
                et = [ct * cp, ct * sp, -st]
                ep = [-sp, cp, 0.0]
                w = [v[0] * er[k] + v[1] * et[k] + v[2] * ep[k] for k in range(3)]
        return [ax[0][k] * w[0] + ax[1][k] * w[1] + ax[2][k] * w[2] for k in range(3)]


def apply_frames(points, point_data, cards, ids, cp, cd, who):
    """``nastran_apply_frames``: move CP != 0 GRIDs to basic (``points`` in
    place), keep ``nastran:cp``/``nastran:cd``, and return the systems."""
    systems = CoordSystems()
    cp = np.asarray(cp, dtype=np.int64)
    cd = np.asarray(cd, dtype=np.int64)
    any_cp = bool(np.count_nonzero(cp))
    any_cd = bool(np.count_nonzero(cd))
    if not any_cp and not any_cd:
        return systems

    pending = {}
    duplicates = 0
    for c in cards:
        if c.cid <= 0:
            continue
        if c.cid in pending:
            duplicates += 1
        else:
            pending[c.cid] = c
    if duplicates:
        warn(
            f"{who}: {duplicates} coordinate system(s) are defined more than once; "
            "the first definition is used"
        )
    index = {int(g): i for i, g in enumerate(ids)}
    resolved = cp == 0
    cpl = cp.tolist()
    degenerate = set()
    progress = True
    while progress:
        progress = False
        for i in range(len(cpl)):
            if resolved[i] or not systems.has(cpl[i]):
                continue
            points[i] = systems.to_basic(cpl[i], points[i].tolist())
            resolved[i] = True
            progress = True
        for cid in sorted(pending):
            c = pending[cid]
            if c.grids is not None:
                rows = [index.get(g) for g in c.grids]
                if any(r is None or not resolved[r] for r in rows):
                    continue
                a, b, p = (points[r].tolist() for r in rows)
            elif systems.has(c.rid):
                a = systems.to_basic(c.rid, c.abc[0:3])
                b = systems.to_basic(c.rid, c.abc[3:6])
                p = systems.to_basic(c.rid, c.abc[6:9])
            else:
                continue
            sys = _system_from_points(c.type, a, b, p)
            if sys is None:
                degenerate.add(cid)
            else:
                systems.systems[cid] = sys
            del pending[cid]
            progress = True
    if degenerate:
        warn(
            f"{who}: {len(degenerate)} coordinate system(s) have coincident or collinear "
            "defining points and are ignored"
        )
    kept = [c for c, r in zip(cpl, resolved) if c != 0 and not r]
    if kept:
        ids_text = ", ".join(str(c) for c in sorted(set(kept)))
        warn(
            f"{who}: {len(kept)} GRID(s) are in coordinate system(s) {ids_text} that cannot "
            "be resolved; their coordinates are kept as written"
        )
    if any_cp:
        point_data["nastran:cp"] = cp
    if any_cd:
        unresolved = sum(1 for c in cd.tolist() if c > 0 and not systems.has(c))
        if unresolved:
            warn(
                f"{who}: {unresolved} GRID(s) have an output system (CD) that cannot be "
                "resolved; their results stay in it"
            )
        point_data["nastran:cd"] = cd
    return systems


def rotate_to_basic(systems, cd, point, values):
    """Rotate consecutive triplets of ``values`` (a list, in place) from system
    ``cd`` at ``point`` to basic; False when there is nothing to do."""
    if cd <= 0 or not systems.has(cd):
        return False
    for k in range(0, len(values), 3):
        values[k : k + 3] = systems.vector_to_basic(cd, point, values[k : k + 3])
    return True
