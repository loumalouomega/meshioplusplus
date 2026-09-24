"""LS-DYNA binary output database (``binout``), read-only.

A binout is an LSDA file (Livermore Software Data Archival): a symbol table of
directories and typed variables whose values sit in DATA records. Each ASCII
database LS-DYNA would otherwise write (``nodout``, ``glstat``, ``matsum``,
``elout``...) is a directory holding ``metadata`` (ids, title...) and one
folder per output time, ``d000001``, ``d000002``..., each with a ``time`` and
the variables of that output.

The steps are ``nodout``'s outputs (else those of the first database with
any): each is a point cloud of the ``nodout`` nodes at their current
coordinates, with their displacements, velocities and accelerations as point
data. Every other database's latest output at or before the step's time is
field data ``binout:<database>:<variable>``, with its own time as
``binout:<database>:time`` (the databases have their own output intervals).

The LSDA layout follows lasso-python's ``lsda_py3`` (BSD-3), itself LSTC's.
"""

import glob
import os
import struct

import numpy as np

from .._exceptions import ReadError
from .._mesh import CellBlock, Mesh

TIME_KEY = "meshio:time"

# LSDA commands and data types
_CD, _DATA, _VARIABLE, _BEGIN, _END, _OFFSET = 2, 3, 4, 5, 6, 7
_LINK = 11
_TYPES = {
    1: "i1",
    2: "i2",
    3: "i4",
    4: "i8",
    5: "u1",
    6: "u2",
    7: "u4",
    8: "u8",
    9: "f4",
    10: "f8",
    11: "u1",
}

# nodout's per-node variables: (point data name, its x, y, z variables)
_NODOUT_VECTORS = (
    ("displacement", ("x_displacement", "y_displacement", "z_displacement")),
    ("rotation", ("rx_displacement", "ry_displacement", "rz_displacement")),
    ("velocity", ("x_velocity", "y_velocity", "z_velocity")),
    ("rotational_velocity", ("rx_velocity", "ry_velocity", "rz_velocity")),
    ("acceleration", ("x_acceleration", "y_acceleration", "z_acceleration")),
    (
        "rotational_acceleration",
        ("rx_acceleration", "ry_acceleration", "rz_acceleration"),
    ),
)
_COORDINATES = ("x_coordinate", "y_coordinate", "z_coordinate")


def _fail(message):
    raise ReadError(f"LS-DYNA binout: {message}")


def is_binout(head):
    """The LSDA header: its size (8 or more), length, offset, command and type
    sizes, the byte order, then the symbol table offset command."""
    if len(head) < 8:
        return False
    size, lsize, osize, csize, tsize, order = head[:6]
    if size < 8 or lsize not in (4, 8) or osize not in (4, 8):
        return False
    if csize != 1 or tsize != 1 or order not in (0, 1):
        return False
    at = size
    if len(head) < at + lsize + csize:
        return True
    return head[at + lsize] == _OFFSET


class _Var:
    __slots__ = ("type", "offset", "length", "file")

    def __init__(self, typ, offset, length, file):
        self.type, self.offset, self.length, self.file = typ, offset, length, file


class _File:
    """One LSDA file: its layout and its bytes."""

    def __init__(self, path):
        with open(path, "rb") as f:
            self.data = f.read()
        d = self.data
        if not is_binout(d[:32]):
            _fail(f"'{path}' is not an LSDA (binout) file")
        size, self.lsize, self.osize, self.csize, self.tsize, order = d[:6]
        self.order = "<" if order == 1 else ">"
        self.start = size
        self.icode = {1: "b", 2: "h", 4: "i", 8: "q"}

    def uint(self, pos, n):
        if pos + n > len(self.data):
            _fail("a record runs past the end of the file")
        return struct.unpack(self.order + self.icode[n], self.data[pos : pos + n])[0]


class _Archive:
    """The symbol table of a binout and its continuation files (``%001``...)."""

    def __init__(self, path):
        paths = [path] + sorted(glob.glob(glob.escape(path) + "%[0-9][0-9]*"))
        self.files = [_File(p) for p in paths]
        self.root = {}
        for f in self.files:
            self._read_symbols(f)

    def _read_symbols(self, f):
        pos = f.start
        lc = f.lsize + f.csize
        _, cmd = f.uint(pos, f.lsize), f.uint(pos + f.lsize, f.csize)
        if cmd != _OFFSET:
            return
        pos += lc
        cwd = [self.root, []]  # the directory dict and its path
        while True:
            offset = f.uint(pos, f.osize)
            if offset == 0:
                return
            pos = offset
            length, cmd = f.uint(pos, f.lsize), f.uint(pos + f.lsize, f.csize)
            pos += lc
            if cmd != _BEGIN:
                return
            while True:
                length, cmd = f.uint(pos, f.lsize), f.uint(pos + f.lsize, f.csize)
                body = length - lc
                pos += lc
                if cmd == _CD:
                    cwd = self._cd(cwd, f.data[pos : pos + body].decode("latin-1"))
                    pos += body
                elif cmd == _VARIABLE:
                    n = body - (f.tsize + f.osize + f.lsize)
                    name = f.data[pos : pos + n].decode("latin-1")
                    p = pos + n
                    typ = f.uint(p, f.tsize)
                    off = f.uint(p + f.tsize, f.osize)
                    count = f.uint(p + f.tsize + f.osize, f.lsize)
                    cwd[0][name] = _Var(typ, off, count, f)
                    pos += body
                else:  # the end of this table: the next table's offset follows
                    break

    def _cd(self, cwd, path):
        node, parts = cwd
        if path.startswith("/"):
            node, parts = self.root, []
            path = path[1:]
        for part in path.rstrip("/").split("/"):
            if not part or part == ".":
                continue
            if part == "..":
                parts = parts[:-1]
                node = self.root
                for p in parts:
                    node = node[p]
                continue
            child = node.get(part)
            if not isinstance(child, dict):
                child = node[part] = {}
            node, parts = child, parts + [part]
        return [node, parts]

    def get(self, path):
        node = self.root
        for part in path:
            if not isinstance(node, dict) or part not in node:
                return None
            node = node[part]
        return node

    def values(self, var):
        """A variable's values as a 1-D array (a link's target's)."""
        seen = 0
        while var.type == _LINK and seen < 16:
            target = bytes(self._raw(var)).decode("latin-1")
            var = self.get([p for p in target.split("/") if p])
            if not isinstance(var, _Var):
                _fail(f"a link points to '{target}', which is not a variable")
            seen += 1
        return self._raw(var)

    def _raw(self, var):
        f = var.file
        name_len = f.data[var.offset + f.lsize + f.csize + f.tsize]
        pos = var.offset + f.lsize + f.csize + f.tsize + 1 + name_len
        dtype = np.dtype(f.order + _TYPES.get(var.type, "u1"))
        end = pos + dtype.itemsize * var.length
        if end > len(f.data):
            _fail("a variable's data runs past the end of the file")
        return np.frombuffer(f.data[pos:end], dtype)


def _is_step(name):
    """Whether a folder name is an output's: ``d`` and digits (``d000001``);
    ``deforc`` is a database."""
    return name[:1] == "d" and name[1:].isdigit()


def _steps(archive, database):
    """The step folders (``d000001``...) of a database, in order."""
    folder = archive.get(database)
    if not isinstance(folder, dict):
        return []
    return sorted(k for k, v in folder.items() if isinstance(v, dict) and _is_step(k))


def _databases(archive):
    """Every database with step folders, as its path: ``nodout``,
    ``elout/beam``..."""
    out = []

    def walk(node, path):
        steps = [k for k, v in node.items() if isinstance(v, dict) and _is_step(k)]
        if steps and path:
            out.append(tuple(path))
        for k, v in sorted(node.items()):
            if isinstance(v, dict) and not _is_step(k) and k != "metadata":
                walk(v, path + [k])

    walk(archive.root, [])
    return out


def _time(archive, database, step):
    var = archive.get(list(database) + [step, "time"])
    if not isinstance(var, _Var):
        return None
    v = archive.values(var)
    return float(v[0]) if len(v) else None


class _Binout:
    def __init__(self, filename):
        path = os.fspath(filename)
        if not os.path.isfile(path):
            _fail(f"'{path}' does not exist")
        self.archive = a = _Archive(path)
        self.databases = _databases(a)
        if not self.databases:
            _fail("the file holds no database with outputs")
        self.main = ("nodout",) if ("nodout",) in self.databases else self.databases[0]
        self.steps = _steps(a, self.main)
        self.times = [_time(a, self.main, s) for s in self.steps]
        self.times = [0.0 if t is None else t for t in self.times]
        # the other databases' outputs: (times, folders) in time order
        self.others = {}
        for db in self.databases:
            if db == self.main:
                continue
            pairs = [(_time(a, db, s), s) for s in _steps(a, db)]
            pairs = sorted((t, s) for t, s in pairs if t is not None)
            self.others[db] = ([t for t, _ in pairs], [s for _, s in pairs])


def time_values(filename):
    return list(_Binout(filename).times)


def _numeric(values):
    """Numbers, not text (titles and dates are 1-byte integers)."""
    return values.dtype.kind in "iuf" and values.dtype.itemsize > 1


def read(filename, points_only=False, arrays=None, time_step=0):
    b = _Binout(filename)
    a = b.archive
    n = len(b.steps)
    index = time_step + n if time_step < 0 else time_step
    if not 0 <= index < n:
        raise ReadError(
            f"time step {time_step} is out of range: the file has {n} step(s)"
        )

    def want(name):
        return arrays is None or name in arrays

    step = b.steps[index]
    time = b.times[index]
    main = list(b.main)
    points = np.zeros((0, 3))
    ids = None
    if b.main == ("nodout",):
        ids_var = a.get(["nodout", "metadata", "ids"])
        ids = a.values(ids_var).astype(np.int64) if ids_var is not None else None
        n_nodes = len(ids) if ids is not None else 0
        folder = a.get(main + [step])
        coords = [folder.get(c) for c in _COORDINATES]
        if n_nodes and all(isinstance(c, _Var) for c in coords):
            points = np.column_stack([a.values(c).astype(np.float64) for c in coords])
        else:
            points = np.full((n_nodes, 3), np.nan)
    mesh = Mesh(points, [CellBlock("vertex", np.arange(len(points)).reshape(-1, 1))])
    mesh.field_data[TIME_KEY] = np.array([time], dtype=np.float64)
    mesh.time_values = list(b.times)
    if ids is not None:
        mesh.point_data["lsdyna:nid"] = ids
    if points_only:
        return mesh

    folder = a.get(main + [step])
    if b.main == ("nodout",):
        used = set(_COORDINATES) | {"time"}
        for name, parts in _NODOUT_VECTORS:
            vs = [folder.get(p) for p in parts]
            if all(isinstance(v, _Var) for v in vs):
                used |= set(parts)
                if want(name):
                    mesh.point_data[name] = np.column_stack(
                        [a.values(v).astype(np.float64) for v in vs]
                    )
        for name, var in sorted(folder.items()):
            if name in used or not isinstance(var, _Var):
                continue
            v = a.values(var)
            if not _numeric(v):
                continue
            key = "binout:nodout:" + name
            if len(v) == len(points) and len(v) > 1:
                if want(name):
                    mesh.point_data[name] = v.astype(np.float64)
            elif want(key):
                mesh.field_data[key] = v.astype(
                    np.float64 if v.dtype.kind == "f" else np.int64
                )
    else:
        _field_data(mesh, a, b.main, step, want)
    # every other database's latest output at or before this time
    for db, (times, folders) in b.others.items():
        k = _latest(times, time)
        if k is not None:
            _field_data(mesh, a, db, folders[k], want, times[k])
    return mesh


def _latest(times, time):
    """The index of the last of ``times`` (sorted) not after ``time`` (their
    float32 roundings compared: a binout writes times in single precision)."""
    t = np.float32(time)
    k = int(np.searchsorted(np.asarray(times, dtype=np.float32), t, side="right")) - 1
    return k if k >= 0 else None


def _field_data(mesh, archive, database, step, want, time=None):
    folder = archive.get(list(database) + [step])
    prefix = "binout:" + "/".join(database) + ":"
    if time is not None and want(prefix + "time"):
        mesh.field_data[prefix + "time"] = np.array([time], dtype=np.float64)
    meta = archive.get(list(database) + ["metadata"])
    if isinstance(meta, dict) and isinstance(meta.get("ids"), _Var):
        ids = archive.values(meta["ids"])
        if _numeric(ids) and want(prefix + "ids"):
            mesh.field_data[prefix + "ids"] = ids.astype(np.int64)
    for name, var in sorted(folder.items()):
        if name == "time" or not isinstance(var, _Var):
            continue
        v = archive.values(var)
        if not _numeric(v) or not want(prefix + name):
            continue
        mesh.field_data[prefix + name] = v.astype(
            np.float64 if v.dtype.kind == "f" else np.int64
        )
