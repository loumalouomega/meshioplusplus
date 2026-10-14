"""LS-DYNA binout files (LSDA), written.

The twin of the writing half of ``src/cpp/src/formats/lsdyna_binout.cpp``: both
engines write the same bytes. The file is the reader's inverse: one LSDA file
(8-byte lengths and offsets, little-endian), the DATA records of every output
and, behind them, the one symbol table of directories and variables.

``nodout`` comes from the mesh: its ids (``lsdyna:nid``), and per step its
coordinates, point data and ``binout:nodout:<variable>`` field data. Every
other database comes from ``binout:<database>:<variable>`` field data, one
output folder whenever the database's own time changes.
"""

import re
import struct

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer

TIME_KEY = "meshio:time"

# LSDA commands and the two data types written.
_CD, _DATA, _VARIABLE, _BEGIN, _END, _OFFSET = 2, 3, 4, 5, 6, 7
_INT32, _FLOAT32 = 3, 9

_HEADER = bytes([8, 8, 8, 1, 1, 1, 0, 0])

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
_VECTOR_NAMES = {name for name, _ in _NODOUT_VECTORS}

_STEP = re.compile(r"^d[0-9]+$")
_INT32_RANGE = (-(2**31), 2**31 - 1)


def _wfail(message):
    raise WriteError(f"LS-DYNA binout: {message}")


def _check_name(name, what):
    if (
        not name
        or len(name.encode("latin-1", "replace")) > 255
        or "/" in name
        or not name.isascii()
    ):
        _wfail(f"{what} '{name}' is not a name a binout can hold")


# -- the container ---------------------------------------------------------------------


class _Lsda:
    """DATA records as they are added, and the entries of the symbol table."""

    def __init__(self):
        self.data = bytearray(_HEADER)
        self.data += struct.pack("<QB", 17, _OFFSET) + struct.pack("<Q", 0)
        self.entries = []  # (directory, name, type, offset, count)

    def var(self, directory, name, values):
        """One variable: float32 for reals, int32 for integers."""
        values = np.asarray(values)
        if values.dtype.kind == "f":
            code, raw = _FLOAT32, _f32(values)
        else:
            code, raw = _INT32, _i32(values, f"'{directory}/{name}'")
        offset = len(self.data)
        encoded = name.encode("latin-1")
        length = 8 + 1 + 1 + 1 + len(encoded) + len(raw)
        self.data += struct.pack("<QBBB", length, _DATA, code, len(encoded))
        self.data += encoded + raw
        self.entries.append((directory, name, code, offset, values.size))

    def finish(self):
        table = len(self.data)
        out = bytearray(self.data)
        out[17:25] = struct.pack("<Q", table)  # after the header, length and command
        out += struct.pack("<QB", 9, _BEGIN)
        cwd = ""
        for directory, name, code, offset, count in self.entries:
            if directory != cwd:
                encoded = directory.encode("latin-1")
                out += struct.pack("<QB", 9 + len(encoded), _CD) + encoded
                cwd = directory
            encoded = name.encode("latin-1")
            out += struct.pack("<QB", 9 + len(encoded) + 1 + 8 + 8, _VARIABLE)
            out += encoded + struct.pack("<BQQ", code, offset, count)
        out += struct.pack("<QB", 17, 6) + struct.pack("<Q", 0)
        return bytes(out)


def _f32(values):
    with np.errstate(over="ignore"):
        return np.ascontiguousarray(values, dtype="<f4").tobytes()


def _i32(values, what):
    wide = np.asarray(values, dtype=np.int64)
    if wide.size and (wide.min() < _INT32_RANGE[0] or wide.max() > _INT32_RANGE[1]):
        _wfail(f"{what} holds an integer outside the 32-bit range")
    return np.ascontiguousarray(wide, dtype="<i4").tobytes()


# -- one step of the mesh --------------------------------------------------------------


def _numeric(array):
    return array.dtype.kind in "biuf"


def _flat(array):
    array = np.asarray(array)
    return array.astype(np.int64) if array.dtype.kind in "bu" else array


class _Output:
    """One database's variables of one step: ``(name, array)`` in file order."""

    def __init__(self, time, variables, ids=None):
        self.time = time
        self.variables = variables
        self.ids = ids

    def layout(self):
        return [(n, "f" if a.dtype.kind == "f" else "i") for n, a in self.variables]


class _Step:
    def __init__(self, mesh, time):
        self.time = float(np.float32(time))
        self.npoints = len(mesh.points)
        self.nodout = None
        self.databases = {}  # name -> _Output
        self.dropped = []
        used = set()
        fields = {}
        for name, array in mesh.field_data.items():
            if name == TIME_KEY:
                continue
            parts = name.split(":")
            if parts[0] != "binout" or len(parts) < 3:
                self.dropped.append(name)
                continue
            database, variable = ":".join(parts[1:-1]), parts[-1]
            fields.setdefault(database, {})[variable] = (name, np.asarray(array))
        if self.npoints:
            self.nodout = self._nodout(mesh, used, fields.get("nodout", {}))
        for database in sorted(fields):
            self._database(database, fields[database])
        for name in mesh.point_data:
            if name not in used:
                self.dropped.append(name)
        self.dropped.extend(mesh.cell_data)
        self.dropped = sorted(set(self.dropped))

    def _nodout(self, mesh, used, items):
        n = self.npoints
        points = np.zeros((n, 3))
        raw = np.asarray(mesh.points, dtype=np.float64).reshape(n, -1)
        points[:, : min(3, raw.shape[1])] = raw[:, :3]
        variables = [(c, points[:, k]) for k, c in enumerate(_COORDINATES)]
        for name, parts in _NODOUT_VECTORS:
            array = mesh.point_data.get(name)
            if array is None:
                continue
            used.add(name)
            array = np.asarray(array)
            if array.size != 3 * n or not _numeric(array):
                self.dropped.append(name)
                continue
            array = array.reshape(n, 3).astype(np.float64)
            variables += [(p, array[:, k]) for k, p in enumerate(parts)]
        ids = mesh.point_data.get("lsdyna:nid")
        used.add("lsdyna:nid")
        if ids is None:
            ids = np.arange(1, n + 1, dtype=np.int64)
        ids = np.asarray(ids)
        if ids.size != n or ids.dtype.kind not in "biu":
            _wfail("'lsdyna:nid' must hold one integer per point")
        for name in sorted(mesh.point_data):
            array = np.asarray(mesh.point_data[name])
            if name in used or name in _VECTOR_NAMES:
                continue
            if array.size != n or array.ndim != 1 or not _numeric(array) or n < 2:
                continue  # dropped below, with the unused
            _check_name(name, "point data")
            if name == "time" or any(name == v for v, _ in variables):
                _wfail(f"point data '{name}' is a variable of nodout of its own")
            used.add(name)
            variables.append((name, _flat(array)))
        for variable in sorted(items):
            name, array = items[variable]
            _check_name(variable, "variable")
            if variable in ("time", "ids") or not _numeric(array):
                self.dropped.append(name)
                continue
            if any(variable == v for v, _ in variables):
                _wfail(f"'{name}' repeats a variable of nodout")
            variables.append((variable, _flat(array).ravel()))
        return _Output(self.time, variables, ids.astype(np.int64))

    def _database(self, database, items):
        if database == "nodout":
            if self.nodout is None:  # no points: nodout is not a database of its own
                self.dropped.extend(name for name, _ in items.values())
            return
        for part in database.split("/"):
            _check_name(part, "database")
            if _STEP.match(part) or part == "metadata":
                _wfail(f"'{database}' is not a name a database can have")
        time = self.time
        ids = None
        variables = []
        for variable in sorted(items):
            name, array = items[variable]
            _check_name(variable, "variable")
            if variable == "time":
                values = array.astype(np.float64).ravel()
                if values.size:
                    time = float(np.float32(values[0]))
                continue
            if variable == "ids":
                if array.dtype.kind in "biu":
                    ids = array.astype(np.int64).ravel()
                else:
                    self.dropped.append(name)
                continue
            if not _numeric(array):
                self.dropped.append(name)
                continue
            variables.append((variable, _flat(array).ravel()))
        self.databases[database] = _Output(time, variables, ids)


# -- the file ----------------------------------------------------------------------------


class _Archive:
    """The whole file: every step added so far."""

    def __init__(self):
        self.lsda = _Lsda()
        self.layout = None  # nodout's variables, and each database's
        self.counts = {}  # database -> outputs written
        self.last = {}  # database -> time of its last output
        self.db_layout = {}
        self.ids = None
        self.npoints = None
        self.steps = 0

    def add(self, step):
        index = self.steps
        if self.steps == 0:
            self.npoints = step.npoints
            self.layout = None if step.nodout is None else step.nodout.layout()
        elif step.npoints != self.npoints:
            _wfail(
                f"series: step {index} has {step.npoints} points, the first "
                f"{self.npoints}; a binout holds one set of nodes -- write one file "
                "per step with '{step}'"
            )
        if step.nodout is not None:
            self._nodout(step, index)
        for database in sorted(step.databases):
            self._output(database, step.databases[database], index)
        self.steps += 1

    def _nodout(self, step, index):
        out = step.nodout
        if self.steps == 0:
            self.ids = out.ids
            self.lsda.var("/nodout/metadata", "ids", out.ids)
        elif out.layout() != self.layout or not np.array_equal(out.ids, self.ids):
            _wfail(
                f"series: step {index}'s nodes or variables differ from the first "
                "step's; a binout holds one set of them"
            )
        self._folder("nodout", out)

    def _output(self, database, out, index):
        layout = out.layout()
        if database not in self.db_layout:
            self.db_layout[database] = (layout, None if out.ids is None else out.ids)
            if out.ids is not None:
                self.lsda.var(f"/{database}/metadata", "ids", out.ids)
        elif layout != self.db_layout[database][0]:
            _wfail(
                f"series: step {index}'s variables of '{database}' differ from "
                "its first output's"
            )
        # An output is new when its time moved: the reader carries the latest.
        if database in self.last and self.last[database] == out.time:
            return
        self._folder(database, out)

    def _folder(self, database, out):
        count = self.counts.get(database, 0) + 1
        self.counts[database] = count
        self.last[database] = out.time
        directory = f"/{database}/d{count:06d}"
        self.lsda.var(directory, "time", np.array([out.time], dtype=np.float32))
        for name, array in out.variables:
            self.lsda.var(directory, name, array)

    def finish(self):
        if not self.counts:
            _wfail("the mesh holds no database: no points and no binout:* field data")
        return self.lsda.finish()


def _warn_dropped(step):
    if step.dropped:
        warn(
            "LS-DYNA binout: arrays without a place in the file are not written: "
            + ", ".join(step.dropped)
        )


def _time_of(mesh):
    array = mesh.field_data.get(TIME_KEY)
    if array is None:
        return 0.0
    values = np.asarray(array, dtype=np.float64).ravel()
    return float(values[0]) if values.size else 0.0


def _save(filename, data):
    if is_buffer(filename, "w"):
        filename.write(data)
    else:
        with open(filename, "wb") as fh:
            fh.write(data)


@_provenance.slotless_writer
def write(filename, mesh):
    step = _Step(mesh, _time_of(mesh))
    _warn_dropped(step)
    archive = _Archive()
    archive.add(step)
    _save(filename, archive.finish())


class SeriesWriter:
    """Write a time series into one binout: ``nodout``'s ids once and an output
    per step, and any other database's output whenever its own
    ``binout:<database>:time`` moves. A step's time is the argument, not its
    ``meshio:time``. Every step must have the first step's points and the
    variables of each database (a database may be absent from a step, and
    appear later). Use as a context manager, or call :meth:`close`."""

    def __init__(self, filename):
        self._filename = filename
        self._archive = _Archive()

    def write(self, time, mesh):
        step = _Step(mesh, float(time))
        if self._archive.steps == 0:
            _warn_dropped(step)
        self._archive.add(step)

    def close(self):
        if self._filename is None:
            return
        filename, self._filename = self._filename, None
        if not self._archive.steps:
            return
        _save(filename, self._archive.finish())

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if exc[0] is None:
            self.close()
        else:
            self._filename = None
        return False
