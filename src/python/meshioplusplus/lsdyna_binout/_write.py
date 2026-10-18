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
_INT8, _INT32, _FLOAT32, _LINK = 1, 3, 9, 11


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

    def __init__(self, word=8, big=False):
        self.word = word
        self.o = ">" if big else "<"
        self.L = "Q" if word == 8 else "I"  # a length or an offset
        self.data = bytearray([8, word, word, 1, 1, 0 if big else 1, 0, 0])
        self.data += struct.pack(self.o + self.L + "B", word + 1 + word, _OFFSET)
        self.data += struct.pack(self.o + self.L, 0)
        self.entries = []  # (directory, name, type, offset, count)

    def var(self, directory, name, values, text=False, link=None):
        """One variable: float32 for reals, int32 for integers, (``text``) the
        1-byte integers LS-DYNA writes titles and dates as, or (``link``) a link
        to the variable at that absolute path."""
        values = np.asarray(values)
        if link is not None:
            code, raw = _LINK, link.encode("latin-1")
            values = np.zeros(len(raw))
        elif text:
            code, raw = _INT8, (values.astype(np.int64) & 0xFF).astype("u1").tobytes()
        elif values.dtype.kind == "f":
            code, raw = _FLOAT32, _f32(values, self.o)
        else:
            code, raw = _INT32, _i32(values, f"'{directory}/{name}'", self.o)
        offset = len(self.data)
        encoded = name.encode("latin-1")
        length = self.word + 1 + 1 + 1 + len(encoded) + len(raw)
        self.data += struct.pack(
            self.o + self.L + "BBB", length, _DATA, code, len(encoded)
        )
        self.data += encoded + raw
        self.entries.append((directory, name, code, offset, values.size))

    def has(self, directory, name):
        return any(d == directory and n == name for d, n, *_ in self.entries)

    def finish(self):
        w, o, L = self.word, self.o, self.L
        table = len(self.data)
        if w == 4 and table > 0xFFFFFFFF:
            _wfail("the file is too large for 4-byte offsets")
        out = bytearray(self.data)
        at = 8 + w + 1  # after the header, the length and the command
        out[at : at + w] = struct.pack(o + L, table)
        out += struct.pack(o + L + "B", w + 1, _BEGIN)
        cwd = ""
        for directory, name, code, offset, count in self.entries:
            if directory != cwd:
                encoded = directory.encode("latin-1")
                out += struct.pack(o + L + "B", w + 1 + len(encoded), _CD) + encoded
                cwd = directory
            encoded = name.encode("latin-1")
            out += struct.pack(o + L + "B", w + 1 + len(encoded) + 1 + w + w, _VARIABLE)
            out += encoded + struct.pack(o + "B" + L + L, code, offset, count)
        out += struct.pack(o + L + "B", w + 1 + w, 6) + struct.pack(o + L, 0)
        return bytes(out)


def _f32(values, order="<"):
    with np.errstate(over="ignore"):
        return np.ascontiguousarray(values, dtype=order + "f4").tobytes()


def _i32(values, what, order="<"):
    wide = np.asarray(values, dtype=np.int64)
    if wide.size and (wide.min() < _INT32_RANGE[0] or wide.max() > _INT32_RANGE[1]):
        _wfail(f"{what} holds an integer outside the 32-bit range")
    return np.ascontiguousarray(wide, dtype=order + "i4").tobytes()


# -- one step of the mesh --------------------------------------------------------------


def _numeric(array):
    return array.dtype.kind in "biuf"


def _flat(array):
    array = np.asarray(array)
    return array.astype(np.int64) if array.dtype.kind in "bu" else array


class _Output:
    """One database's variables of one step: ``(name, array)`` in file order."""

    def __init__(self, time, variables, ids=None, text=(), links=None):
        self.time = time
        self.variables = variables
        self.ids = ids
        self.text = list(text)  # (name, bytes) of the database's metadata
        self.links = links or {}  # variable -> the absolute path it is a link to

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
        text = []
        links = {}
        for variable in sorted(items):
            name, array = items[variable]
            if variable.startswith("@"):
                self._text(text, variable, name, array, links)
                continue
            _check_name(variable, "variable")
            if variable in ("time", "ids") or not _numeric(array):
                self.dropped.append(name)
                continue
            if any(variable == v for v, _ in variables):
                _wfail(f"'{name}' repeats a variable of nodout")
            variables.append((variable, _flat(array).ravel()))
        return _Output(self.time, variables, ids.astype(np.int64), text, links)

    def _text(self, text, variable, name, array, links):
        """``@<name>``: a metadata text variable (bytes), kept in ``text``;
        ``@link=<variable>``: the path that variable links to, kept in ``links``."""
        label = variable[1:]
        if label.startswith("link="):
            _check_name(label[5:], "linked variable")
            if array.dtype.kind not in "iu":
                self.dropped.append(name)
                return
            target = bytes(
                int(v) & 0xFF for v in array.astype(np.int64).ravel()
            ).decode("latin-1")
            if not target.startswith("/"):
                _wfail(f"the link '{name}' must hold an absolute path")
            links[label[5:]] = target
            return
        _check_name(label, "text variable")
        if label == "ids" or array.dtype.kind not in "iu":
            self.dropped.append(name)
            return
        values = array.astype(np.int64).ravel()
        if values.size and (values.min() < -128 or values.max() > 255):
            _wfail(f"the text '{name}' holds a value that is not a byte")
        text.append((label, values))

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
        text = []
        links = {}
        for variable in sorted(items):
            name, array = items[variable]
            if variable.startswith("@"):
                self._text(text, variable, name, array, links)
                continue
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
        self.databases[database] = _Output(time, variables, ids, text, links)


# -- the file ----------------------------------------------------------------------------


class _Archive:
    """The whole file: every step added so far."""

    def __init__(self, split_size=0, word_size=8, big_endian=False):
        if word_size not in (4, 8):
            _wfail("the word size must be 4 or 8")
        self.word, self.big = word_size, big_endian
        self.files = [_Lsda(word_size, big_endian)]
        self.split = int(split_size or 0)
        self.folders = 0  # folders in the current file
        self.layout = None  # nodout's variables, and each database's
        self.counts = {}  # database -> outputs written
        self.last = {}  # database -> time of its last output
        self.db_layout = {}
        self.ids = None
        self.npoints = None
        self.steps = 0
        self.link_targets = []  # every link written, to check it resolves

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
            for label, values in out.text:
                self.lsda.var("/nodout/metadata", label, values, text=True)
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
            for label, values in out.text:
                self.lsda.var(f"/{database}/metadata", label, values, text=True)
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
        time = np.array([out.time], dtype=np.float32)
        # A folder stays whole in one file: a new file starts when it would pass the split.
        if self.split > 0 and self.folders > 0:
            size = _record_size("time", time)
            size += sum(_record_size(n, a) for n, a in out.variables)
            if len(self.lsda.data) + size > self.split:
                self.files.append(_Lsda(self.word, self.big))
                self.folders = 0
        self.folders += 1
        self.lsda.var(directory, "time", time)
        for name, array in out.variables:
            if name in out.links:
                self.lsda.var(directory, name, array, link=out.links[name])
                self.link_targets.append(out.links[name])
            else:
                self.lsda.var(directory, name, array)

    @property
    def lsda(self):
        return self.files[-1]

    def finish(self):
        if not self.counts:
            _wfail("the mesh holds no database: no points and no binout:* field data")
        # A link whose target is not in the file could not be read back.
        for target in self.link_targets:
            directory, _, name = target.rpartition("/")
            if not any(f.has(directory or "/", name) for f in self.files):
                _wfail(
                    f"a link points to '{target}', which is not a variable of the file"
                )
        return [f.finish() for f in self.files]


def _record_size(name, values):
    """The bytes of one DATA record of a variable (4 per value; text is 1)."""
    return 8 + 3 + len(name) + 4 * np.asarray(values).size


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


def _save(filename, files):
    """The binout, then ``<name>%001``, ``%002``... for its continuation files."""
    if is_buffer(filename, "w"):
        if len(files) > 1:
            _wfail("a split binout needs a file path, not a stream")
        filename.write(files[0])
        return
    for k, data in enumerate(files):
        path = str(filename) + (f"%{k:03d}" if k else "")
        with open(path, "wb") as fh:
            fh.write(data)


@_provenance.slotless_writer
def write(filename, mesh, split_size=0, word_size=8, big_endian=False):
    step = _Step(mesh, _time_of(mesh))
    _warn_dropped(step)
    archive = _Archive(split_size, word_size, big_endian)
    archive.add(step)
    _save(filename, archive.finish())


class SeriesWriter:
    """Write a time series into one binout: ``nodout``'s ids once and an output
    per step, and any other database's output whenever its own
    ``binout:<database>:time`` moves. A step's time is the argument, not its
    ``meshio:time``. Every step must have the first step's points and the
    variables of each database (a database may be absent from a step, and
    appear later). Use as a context manager, or call :meth:`close`."""

    def __init__(self, filename, split_size=0, word_size=8, big_endian=False):
        self._filename = filename
        self._archive = _Archive(split_size, word_size, big_endian)

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
