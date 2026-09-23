"""
I/O for the MSC Nastran HDF5 result database (``.h5``): the Python reference twin
of ``src/cpp/src/formats/nastran_h5.cpp``. Read-only.

The model comes from ``/NASTRAN/INPUT`` (GRID points, element cards as cell blocks,
one cell region per property id), and every result domain an
``/INDEX/NASTRAN/RESULT/...`` table references is one step of a sequence: its
``/NASTRAN/RESULT/NODAL`` tables become point data and its
``/NASTRAN/RESULT/ELEMENTAL`` tables cell data. See ``doc/formats/nastran_h5.md``.
"""

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._mesh import CellBlock, Mesh
from .._regions import Region

__all__ = ["read", "time_values"]

TIME_KEY = "meshio:time"
_NAN = float("nan")

_GRID = "/NASTRAN/INPUT/NODE/GRID"
_ELEMENTS = "/NASTRAN/INPUT/ELEMENT"
_PROPERTIES = "/NASTRAN/INPUT/PROPERTY"
_DOMAINS = "/NASTRAN/RESULT/DOMAINS"
_NODAL = "/NASTRAN/RESULT/NODAL"
_ELEMENTAL = "/NASTRAN/RESULT/ELEMENTAL"

# Nastran numbers the hex20/wedge15 mid-side nodes bottom, vertical, top;
# meshio++ (VTK) numbers them bottom, top, vertical: conn[k] = G[perm[k]].
_HEXA20 = list(range(12)) + [16, 17, 18, 19, 12, 13, 14, 15]
_PENTA15 = list(range(9)) + [12, 13, 14, 9, 10, 11]

# card -> (linear type, nodes, quadratic type or None, nodes, permutation or None)
_CARDS = {
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


def _fail(message):
    raise ReadError(f"MSC Nastran HDF5: {message}")


def _open(filename):
    import h5py

    try:
        return h5py.File(str(filename), "r")
    except OSError:
        _fail(f"'{filename}' is not an HDF5 file that can be opened")


def _children(f, path, groups):
    import h5py

    if path not in f or not isinstance(f[path], h5py.Group):
        return []
    kind = h5py.Group if groups else h5py.Dataset
    return sorted(k for k, v in f[path].items() if isinstance(v, kind))


def _is_dataset(f, path):
    import h5py

    return path in f and isinstance(f[path], h5py.Dataset)


def _names(ds):
    return ds.dtype.names or ()


def _is_float(ds, name):
    t = ds.dtype.fields[name][0]
    return (t.base if t.subdtype else t).kind == "f"


def _first(values):
    """A member's values, reduced to the first entry per row for an array member."""
    a = np.asarray(values, dtype=np.float64)
    return a.reshape(a.shape[0], -1)[:, 0] if a.ndim > 1 else a


class _File:
    """The file's structure: everything but the result payloads."""

    def __init__(self, filename):
        self.f = f = _open(filename)
        if not _is_dataset(f, _GRID):
            _fail(
                f"'{filename}' has no /NASTRAN/INPUT/NODE/GRID table; it is not an MSC "
                "Nastran HDF5 result file"
            )
        version = f["/NASTRAN"].attrs.get("VERSION")
        if version is not None:
            version = np.asarray(version).ravel()
            text = version[0] if version.size else b""
            text = (
                text.decode("ascii", "replace")
                if isinstance(text, bytes)
                else str(text)
            )
            text = text.rstrip(" \x00")
            if not text.lstrip(" \t").lower().startswith("msc"):
                _fail(
                    f"'{filename}' was written by '{text}', not MSC Nastran; other "
                    "vendors' HDF5 schemas are not supported"
                )
        self.grid_ids = np.asarray(f[_GRID]["ID"], dtype=np.int64)
        self.grid_index = {}
        for i, g in enumerate(self.grid_ids.tolist()):
            if g in self.grid_index:
                _fail(f"GRID {g} is defined twice")
            self.grid_index[g] = i
        self.scalar_points = set()  # SPOINT/EPOINT ids
        for table in ("/NASTRAN/INPUT/NODE/SPOINT", "/NASTRAN/INPUT/NODE/EPOINT"):
            if _is_dataset(f, table) and "ID" in _names(f[table]):
                self.scalar_points.update(
                    np.asarray(f[table]["ID"], dtype=np.int64).tolist()
                )
        self.tables = []
        self.steps = []
        self._read_result_tables()

    def _add_table(self, path, group, name, nodal):
        f = self.f
        ds = f[path]
        names = _names(ds)
        if "DOMAIN_ID" not in names:
            return  # a lookup table (ENERGY/IDENT), not a result
        if not nodal and "EID" in names:
            key = "EID"
        elif "ID" in names:
            key = "ID"
        else:
            warn(f"MSC Nastran HDF5: result table {path} has no ID/EID column; skipped")
            return
        index = "/INDEX" + path
        if not _is_dataset(f, index):
            warn(f"MSC Nastran HDF5: result table {path} has no {index} table; skipped")
            return
        idx = f[index][()]
        rows = ds.shape[0]
        table_index = {}
        for d, p, n in zip(
            idx["DOMAIN_ID"].tolist(), idx["POSITION"].tolist(), idx["LENGTH"].tolist()
        ):
            if p < 0 or n < 0 or p + n > rows:
                _fail(f"{index} gives rows outside {path}")
            table_index[d] = (p, n)
        self.tables.append(
            {
                "path": path,
                "group": group,
                "name": name,
                "nodal": nodal,
                "key": key,
                "index": table_index,
            }
        )

    def _read_result_tables(self):
        f = self.f
        for name in _children(f, _NODAL, False):
            self._add_table(f"{_NODAL}/{name}", "", name, True)
        for group in _children(f, _ELEMENTAL, True):
            gpath = f"{_ELEMENTAL}/{group}"
            for name in _children(f, gpath, False):
                self._add_table(f"{gpath}/{name}", group, name, False)
        if not self.tables:
            return
        if not _is_dataset(f, _DOMAINS):
            _fail("the file has result tables but no /NASTRAN/RESULT/DOMAINS table")
        dom = f[_DOMAINS][()]
        names = _names(f[_DOMAINS])
        if "ID" not in names:
            _fail("/NASTRAN/RESULT/DOMAINS has no ID column")
        n = len(dom)

        def col(name, dtype):
            if name in names:
                return np.asarray(dom[name], dtype=dtype).tolist()
            return [dtype(0)] * n

        ids = col("ID", np.int64)
        cols = {
            "subcase": col("SUBCASE", np.int64),
            "step": col("STEP", np.int64),
            "analysis": col("ANALYSIS", np.int64),
            "mode": col("MODE", np.int64),
            "time": col("TIME_FREQ_EIGR", np.float64),
            "eigi": col("EIGI", np.float64),
        }
        used = {d for t in self.tables for d in t["index"]}
        known = set()
        for i in range(n):
            known.add(ids[i])
            if ids[i] not in used:
                continue
            step = {"domain": ids[i]}
            step.update({k: v[i] for k, v in cols.items()})
            self.steps.append(step)
            used.discard(ids[i])  # a repeated DOMAINS row is one step
        for missing in sorted(used):
            if missing not in known:
                _fail(
                    f"a result table names domain {missing}, which "
                    "/NASTRAN/RESULT/DOMAINS does not define"
                )


def _nodal_outputs(table, ds):
    """The point-data arrays a nodal table yields: [(name, [members])]."""
    names = _names(ds)
    floats = [
        m for m in names if m not in (table["key"], "DOMAIN_ID") and _is_float(ds, m)
    ]
    left = set(floats)
    out = []
    name = table["name"]
    cplx = name.endswith("_CPLX") and len(name) > len("_CPLX")
    base = name[: -len("_CPLX")] if cplx else name

    def group(out_name, members):
        if all(m in left for m in members):
            left.difference_update(members)
            out.append((out_name, members))

    if cplx:
        group(base + "_real", ["XR", "YR", "ZR"])
        group(base + "_imag", ["XI", "YI", "ZI"])
        group(base + "_ROT_real", ["RXR", "RYR", "RZR"])
        group(base + "_ROT_imag", ["RXI", "RYI", "RZI"])
    else:
        group(base, ["X", "Y", "Z"])
        group(base + "_ROT", ["RX", "RY", "RZ"])
    if left == {"VALUE"}:
        out.append((name, ["VALUE"]))
        left.clear()
    out.extend((f"{name}:{m}", [m]) for m in floats if m in left)
    return out


def _read_cells(nf):
    f = nf.f
    blocks = []  # [type, conn list, eid list, pid list, card]
    skipped_cards = []
    dropped = 0
    for card in _children(f, _ELEMENTS, False):
        spec = _CARDS.get(card)
        if spec is None:
            skipped_cards.append(card)
            continue
        lin_type, lin_n, quad_type, quad_n, perm = spec
        path = f"{_ELEMENTS}/{card}"
        data = f[path][()]
        names = _names(f[path])
        n = len(data)
        if "EID" not in names:
            _fail(f"{path} has no EID column")
        eid = np.asarray(data["EID"], dtype=np.int64).tolist()
        pid = (
            np.asarray(data["PID"], dtype=np.int64).tolist()
            if "PID" in names
            else [-1] * n
        )
        if "G" in names:
            g = np.asarray(data["G"], dtype=np.int64).reshape(n, -1)
        else:
            a, b = ("GA", "GB") if "GA" in names else ("G1", "G2")
            if a not in names or b not in names:
                _fail(f"{path} has no G, GA/GB or G1/G2 columns")
            g = np.stack(
                [
                    np.asarray(data[a], dtype=np.int64),
                    np.asarray(data[b], dtype=np.int64),
                ],
                axis=1,
            ).reshape(n, 2)
        width = g.shape[1]
        if width < lin_n:
            _fail(f"{path} has {width} node columns, {card} needs {lin_n}")
        linear = [lin_type, [], [], [], card]
        quadratic = [quad_type, [], [], [], card]
        partial = 0
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
                p = nf.grid_index.get(row[k])
                if p is None:
                    if row[k] not in nf.scalar_points:
                        what = "no node" if row[k] == 0 else f"undefined GRID {row[k]}"
                        _fail(f"{card} {eid[i]} references {what} as its node {k + 1}")
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
                f"MSC Nastran HDF5: {partial} {card} element(s) have only some "
                f"mid-side nodes; read as {lin_type}"
            )
        for b in (linear, quadratic):
            if b[2]:
                blocks.append(b)
    if skipped_cards:
        warn(
            "MSC Nastran HDF5: skipped element tables with no cell type: "
            + ", ".join(skipped_cards)
        )
    if dropped:
        warn(
            f"MSC Nastran HDF5: skipped {dropped} element(s) that connect scalar points"
        )
    return blocks


def _property_types(f):
    ptype = {}
    for prop in _children(f, _PROPERTIES, False):
        ds = f[f"{_PROPERTIES}/{prop}"]
        if "PID" in _names(ds):
            for p in np.asarray(ds["PID"], dtype=np.int64).tolist():
                ptype.setdefault(p, prop)
    # Grouped properties (PCOMP/IDENTITY) live one level deeper.
    for prop in _children(f, _PROPERTIES, True):
        path = f"{_PROPERTIES}/{prop}/IDENTITY"
        if _is_dataset(f, path) and "PID" in _names(f[path]):
            for p in np.asarray(f[path]["PID"], dtype=np.int64).tolist():
                ptype.setdefault(p, prop)
    return ptype


def _resolve_step(time_step, count):
    if count == 0:
        if time_step in (0, -1):
            return None
        raise ReadError(
            f"meshio++: time step {time_step} requested, but this file carries no time steps"
        )
    step = time_step + count if time_step < 0 else time_step
    if not 0 <= step < count:
        raise ReadError(
            f"meshio++: time step {time_step} is out of range: this file has {count} "
            + ("step" if count == 1 else "steps")
        )
    return step


def read(filename, points_only=False, arrays=None, time_step=0):
    nf = _File(filename)
    try:
        return _read(nf, points_only, arrays, time_step)
    finally:
        nf.f.close()


def _read(nf, points_only, arrays, time_step):
    f = nf.f
    grid = f[_GRID][()]
    npts = len(grid)
    points = np.asarray(grid["X"], dtype=np.float64).reshape(npts, 3)
    point_data = {}
    for frame in ("CP", "CD"):
        if frame not in grid.dtype.names:
            continue
        v = np.asarray(grid[frame], dtype=np.int64)
        nonzero = int(np.count_nonzero(v))
        if nonzero == 0:
            continue
        if frame == "CP":
            warn(
                f"MSC Nastran HDF5: {nonzero} GRID(s) have CP != 0; their coordinates "
                "are kept in the local system, not transformed"
            )
        else:
            warn(
                f"MSC Nastran HDF5: {nonzero} GRID(s) have CD != 0; their results are "
                "in the local output system"
            )
        point_data["nastran:" + frame.lower()] = v

    blocks = _read_cells(nf)
    cells = []
    cell_data = {}
    # EID -> global cell, for the element results. A CONM2 has none, and MSC
    # accepts one sharing its id with a structural element, so it stays out.
    cell_index = {}
    offsets = []
    ncells = 0
    shared = 0
    for ctype, conn, eid, pid, card in blocks:
        offsets.append(ncells)
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
            f"MSC Nastran HDF5: {shared} element id(s) are used by more than one card; "
            "their element results go to the first"
        )
    if blocks:
        cell_data["nastran:eid"] = [np.asarray(b[2], dtype=np.int64) for b in blocks]
        cell_data["nastran:pid"] = [np.asarray(b[3], dtype=np.int64) for b in blocks]

    ptype = _property_types(f)
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

    mesh = Mesh(points, cells, point_data=point_data, cell_data=cell_data)
    mesh.regions = sorted(regions, key=lambda r: r.name)
    mesh.time_values = [s["time"] for s in nf.steps]

    s = _resolve_step(time_step, len(nf.steps))
    if s is None:
        return mesh
    dom = nf.steps[s]
    mesh.field_data[TIME_KEY] = np.array([dom["time"]], dtype=np.float64)
    for key in ("domain", "subcase", "step", "analysis", "mode"):
        mesh.field_data["nastran:" + key] = np.array([dom[key]], dtype=np.int64)
    mesh.field_data["nastran:eigi"] = np.array([dom["eigi"]], dtype=np.float64)
    if points_only or (arrays is not None and len(arrays) == 0):
        return mesh

    def wants(name):
        return arrays is None or name in arrays

    used_point_names = []
    cell_arrays = {}
    for t in nf.tables:
        rng = t["index"].get(dom["domain"])
        if rng is None:
            continue
        row0, count = rng
        ds = f[t["path"]]
        if t["nodal"]:
            outputs = []
            for name, members in _nodal_outputs(t, ds):
                out_name = name
                k = 2
                while out_name in used_point_names:
                    out_name = f"{name}_{k}"
                    k += 1
                if wants(out_name):
                    outputs.append((out_name, members))
        else:
            outputs = [
                (f"{t['group']}:{m}", [m])
                for m in _names(ds)
                if m not in (t["key"], "DOMAIN_ID")
                and _is_float(ds, m)
                and wants(f"{t['group']}:{m}")
            ]
        if not outputs:
            continue
        rows = ds[row0 : row0 + count]
        keys = np.asarray(rows[t["key"]], dtype=np.int64).tolist()
        if len(set(keys)) != len(keys):
            warn(
                f"MSC Nastran HDF5: {t['path']} has several rows per "
                f"{'node' if t['nodal'] else 'element'}; skipped"
            )
            continue
        lookup = nf.grid_index if t["nodal"] else cell_index
        pairs = [(r, lookup[k]) for r, k in enumerate(keys) if k in lookup]
        if not pairs:
            continue
        src = np.array([p[0] for p in pairs], dtype=np.int64)
        dst = np.array([p[1] for p in pairs], dtype=np.int64)
        for name, members in outputs:
            if t["nodal"]:
                nc = len(members)
                data = np.full((npts, nc), _NAN)
                for c, m in enumerate(members):
                    data[dst, c] = _first(rows[m])[src]
                used_point_names.append(name)
                mesh.point_data[name] = data[:, 0] if nc == 1 else data
            else:
                values = cell_arrays.get(name)
                if values is None:
                    values = cell_arrays[name] = np.full(ncells, _NAN)
                values[dst] = _first(rows[members[0]])[src]
    for name, values in cell_arrays.items():
        mesh.cell_data[name] = [
            values[offsets[b] : offsets[b] + len(blocks[b][2])]
            for b in range(len(blocks))
        ]
    return mesh


def time_values(filename):
    """The ``TIME_FREQ_EIGR`` of every result domain, in step order."""
    nf = _File(filename)
    try:
        return [s["time"] for s in nf.steps]
    finally:
        nf.f.close()
