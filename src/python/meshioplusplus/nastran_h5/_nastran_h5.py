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
from .._mesh import Mesh
from ..nastran._model import CARDS, CoordCard, add_cells, apply_frames, rotate_to_basic

__all__ = ["read", "time_values"]

TIME_KEY = "meshio:time"
_NAN = float("nan")

_GRID = "/NASTRAN/INPUT/NODE/GRID"
_ELEMENTS = "/NASTRAN/INPUT/ELEMENT"
_PROPERTIES = "/NASTRAN/INPUT/PROPERTY"
_DOMAINS = "/NASTRAN/RESULT/DOMAINS"
_NODAL = "/NASTRAN/RESULT/NODAL"
_ELEMENTAL = "/NASTRAN/RESULT/ELEMENTAL"


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


def _read_cards(nf):
    """Every element table with a cell type, as ``(card, eid, pid, g)`` rows."""
    f = nf.f
    cards = []
    skipped_cards = []
    for card in _children(f, _ELEMENTS, False):
        spec = CARDS.get(card)
        if spec is None:
            skipped_cards.append(card)
            continue
        lin_n = spec[1]
        path = f"{_ELEMENTS}/{card}"
        data = f[path][()]
        names = _names(f[path])
        n = len(data)
        if "EID" not in names:
            _fail(f"{path} has no EID column")
        eid = np.asarray(data["EID"], dtype=np.int64)
        pid = (
            np.asarray(data["PID"], dtype=np.int64)
            if "PID" in names
            else np.full(n, -1, dtype=np.int64)
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
        cards.append((card, eid, pid, g))
    if skipped_cards:
        warn(
            "MSC Nastran HDF5: skipped element tables with no cell type: "
            + ", ".join(skipped_cards)
        )
    return cards


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


_COORD_TABLES = (
    ("CORD2R", 1),
    ("CORD2C", 2),
    ("CORD2S", 3),
    ("CORD1R", 1),
    ("CORD1C", 2),
    ("CORD1S", 3),
)


def _coord_cards(f):
    """The CORD1R/C/S and CORD2R/C/S tables of /NASTRAN/INPUT/COORDINATE_SYSTEM."""
    import h5py

    out = []
    for name, ctype in _COORD_TABLES:
        path = "/NASTRAN/INPUT/COORDINATE_SYSTEM/" + name
        if not isinstance(f.get(path), h5py.Dataset):
            continue
        ds = f[path]
        by_grids = name[4] == "1"
        needed = (
            ["CID", "G1", "G2", "G3"]
            if by_grids
            else ["CID", "RID", "A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"]
        )
        if not all(m in ds.dtype.names for m in needed):
            warn(
                f"MSC Nastran HDF5: {path} lacks the expected columns; its systems are not read"
            )
            continue
        rows = ds[()]
        for r in rows:
            if by_grids:
                out.append(
                    CoordCard(r["CID"], ctype, grids=[r["G1"], r["G2"], r["G3"]])
                )
            else:
                out.append(
                    CoordCard(
                        r["CID"], ctype, rid=r["RID"], abc=[r[m] for m in needed[2:]]
                    )
                )
    return out


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
    names = grid.dtype.names
    zeros = np.zeros(npts, dtype=np.int64)
    cp = np.asarray(grid["CP"], dtype=np.int64) if "CP" in names else zeros
    cd = np.asarray(grid["CD"], dtype=np.int64) if "CD" in names else zeros
    ids = np.asarray(grid["ID"], dtype=np.int64)
    systems = apply_frames(
        points, point_data, _coord_cards(f), ids, cp, cd, "MSC Nastran HDF5"
    )
    cdl = cd.tolist()

    cells, cell_data, regions, cell_index, offsets, sizes = add_cells(
        _read_cards(nf),
        nf.grid_index,
        nf.scalar_points,
        _property_types(f),
        "MSC Nastran HDF5",
    )
    ncells = sum(sizes)

    mesh = Mesh(points, cells, point_data=point_data, cell_data=cell_data)
    mesh.regions = regions
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
                if nc == 3:  # vector results are in each GRID's output system (CD)
                    for p in dst.tolist():
                        v = data[p].tolist()
                        if rotate_to_basic(systems, cdl[p], points[p].tolist(), v):
                            data[p] = v
                used_point_names.append(name)
                mesh.point_data[name] = data[:, 0] if nc == 1 else data
            else:
                values = cell_arrays.get(name)
                if values is None:
                    values = cell_arrays[name] = np.full(ncells, _NAN)
                values[dst] = _first(rows[members[0]])[src]
    for name, values in cell_arrays.items():
        mesh.cell_data[name] = [
            values[offsets[b] : offsets[b] + sizes[b]] for b in range(len(sizes))
        ]
    return mesh


def time_values(filename):
    """The ``TIME_FREQ_EIGR`` of every result domain, in step order."""
    nf = _File(filename)
    try:
        return [s["time"] for s in nf.steps]
    finally:
        nf.f.close()
