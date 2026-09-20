"""Transient VTKHDF, the pure-Python (h5py) writer and reader.

The Python twin of ``meshioplusplus._core.VtkhdfTimeSeriesWriter`` and the VTKHDF
counterpart of ``meshioplusplus.xdmf.TimeSeriesWriter``: one static grid, then one
step at a time. Geometry x 1, fields x N -- exactly what VTKHDF's ``Steps`` offset
tables are for -- and the two writers produce the same on-disk layout.

Documented API (mirrors the XDMF twin, raw arrays rather than a ``Mesh``)::

    with TimeSeriesWriter("out.vtkhdf") as w:
        w.write_points_cells(points, cells)
        for k, t in enumerate(times):
            w.write_data(t, point_data={"u": u[k]}, cell_data={"p": p[k]})

Nothing is buffered: every step lands in the file as it is written and ``NSteps``
is rewritten after each, so a run that is killed leaves a file ParaView opens. The
set of array names is fixed at the first ``write_data``; a later step that
introduces, drops or reshapes an array raises ``WriteError`` naming it, because VTK
indexes ``PointDataOffsets/<name>`` by step and a short table is a silent truncation.
"""

from __future__ import annotations

import os

import numpy as np

from .. import _provenance
from .._exceptions import ReadError, WriteError
from .._files import is_buffer
from .._mesh import Mesh
from . import _vtkhdf as _v

_BOOKKEEPING_ROWS = 256
_CHUNK_BYTES_CAP = 4 << 20
_ZERO_TABLES = ("PartOffsets", "PointOffsets", "CellOffsets", "ConnectivityIdOffsets")
_POLY_TABLES = (
    "FaceConnectivityOffsets",
    "FaceOffsetsOffsets",
    "PolyhedronToFaceIdOffsets",
)


def _row_shape(arr):
    return (arr.shape[1],) if arr.ndim == 2 else ()


def _chunk_rows(rows, arr):
    """One step's rows per chunk (the unit of random access), capped for huge meshes."""
    row_bytes = max(1, int(np.prod(_row_shape(arr) or (1,))) * arr.dtype.itemsize)
    return max(1, min(int(rows), max(1, _CHUNK_BYTES_CAP // row_bytes)))


class TimeSeriesWriter:
    """Transient VTKHDF ``UnstructuredGrid`` writer (h5py).

    :param compression: ``None`` (default) or ``"gzip"``.
    :param compression_opts: gzip level 0-9.
    :param mode: ``"truncate"`` starts a fresh series; ``"append"`` continues the one
        already at ``filename`` (a path that does not exist yet is just a fresh series).

    Usable only as a context manager. There is also a C++ writer, reachable
    explicitly as ``meshioplusplus._core.VtkhdfTimeSeriesWriter``; it takes whole
    ``Mesh`` objects and writes the same layout.
    """

    def __init__(
        self, filename, compression=None, compression_opts=4, mode="truncate"
    ) -> None:
        if compression not in (None, "gzip"):
            raise WriteError(
                f"meshio++: vtkhdf: compression must be 'gzip' or None, got {compression!r}"
            )
        if mode not in ("truncate", "append"):
            raise WriteError(
                f"meshio++: vtkhdf: mode must be 'truncate' or 'append', got {mode!r}"
            )
        if is_buffer(filename, "w"):
            raise WriteError(
                "meshio++: vtkhdf: a time series needs a path, not a buffer"
            )
        self.filename = str(filename)
        self._gzip = None
        if compression == "gzip":
            self._gzip = 4 if compression_opts is None else int(compression_opts)
            if not 0 <= self._gzip <= 9:
                raise WriteError("meshio++: vtkhdf: gzip level must be 0-9")
        self._mode = mode
        self._file = None
        self._root = None
        self._geometry = False
        self._adopted = False
        self._has_poly = False
        self._num_points = 0
        self._num_cells = 0
        self._num_steps = 0
        self._fixed = False
        self._info = {"PointData": {}, "CellData": {}, "FieldData": {}}
        self._field_tuples = {}
        self.auto_flush = True
        self.finalized = False

    # -- context ------------------------------------------------------------- #
    def __enter__(self):
        import h5py

        if self._mode == "append" and os.path.exists(self.filename):
            self._file = h5py.File(self.filename, "r+")
            self._adopt()
        else:
            self._file = h5py.File(self.filename, "w")
        return self

    def __exit__(self, *_):
        self.finalize()
        return False

    def _require_open(self):
        if self._file is None or self.finalized:
            raise WriteError(
                "meshio++: vtkhdf: the series is not open (use it as a context manager, "
                "and not after finalize())"
            )

    @property
    def num_steps(self):
        return self._num_steps

    # -- appending to an existing series ------------------------------------- #
    def _adopt(self):
        why = (
            f"meshio++: vtkhdf: '{self.filename}' is not a transient UnstructuredGrid "
            "this writer can continue"
        )
        f = self._file
        if _v.ROOT not in f:
            raise WriteError(why)
        root = f[_v.ROOT]
        if _v._text_attr(root, "Type") != _v.UNSTRUCTURED_GRID or "Steps" not in root:
            raise WriteError(why)
        npts, ncells = (
            root["NumberOfPoints"][()].ravel(),
            root["NumberOfCells"][()].ravel(),
        )
        if npts.size != 1 or ncells.size != 1:
            raise WriteError(why + " (it is partitioned)")
        self._num_points, self._num_cells = int(npts[0]), int(ncells[0])
        self._has_poly = "FaceConnectivity" in root
        self._num_steps = int(root["Steps"].attrs["NSteps"])
        for grp in self._info:
            for name in root[grp]:
                ds = root[grp][name]
                self._info[grp][name] = (ds.dtype, tuple(ds.shape[1:]))
        for name, ds in ((n, root["FieldData"][n]) for n in root["FieldData"]):
            self._field_tuples[name] = int(ds.shape[0])
        self._root = root
        self._fixed = self._num_steps > 0
        self._geometry = True
        self._adopted = True

    # -- geometry ------------------------------------------------------------- #
    def write_points_cells(self, points, cells) -> None:
        """Write the static grid. Once, before the first ``write_data``.

        When continuing an existing series this only checks that the counts match.
        """
        self._require_open()
        mesh = Mesh(np.asarray(points), cells)
        n_cells = sum(len(cb) for cb in mesh.cells)
        if self._adopted:
            self._adopted = False
            if len(mesh.points) != self._num_points or n_cells != self._num_cells:
                raise WriteError(
                    f"meshio++: vtkhdf: the appended series has {self._num_points} points "
                    f"and {self._num_cells} cells but the mesh has {len(mesh.points)} "
                    f"and {n_cells}"
                )
            return
        if self._geometry:
            raise WriteError("meshio++: vtkhdf: write_points_cells was already called")
        self._has_poly = _v._has_polyhedra(mesh)
        version = _v._resolve_version(None, _v.UNSTRUCTURED_GRID, self._has_poly)
        root = self._file.create_group(_v.ROOT, track_order=True)
        _v._write_leaf(root, mesh, _v.UNSTRUCTURED_GRID, self._gzip, version)
        provenance = "\n".join(_provenance.lines(_provenance.SlotTier.BLOCK))
        if provenance:
            root.attrs[_v.PROVENANCE_ATTR] = provenance
        self._root = root
        self._num_points, self._num_cells = len(mesh.points), n_cells
        self._geometry = True
        if self.auto_flush:
            self._file.flush()

    # -- steps ---------------------------------------------------------------- #
    def _prep(self, name, arr, what, rows=None):
        a = _v._data_array(name, arr, what)
        if rows is not None and a.shape[0] != rows:
            raise WriteError(
                f"meshio++: vtkhdf: {what} '{name}' has {a.shape[0]} rows; the series "
                f"has {rows} {'points' if what == 'point_data' else 'cells'}"
            )
        return a

    def _check_set(self, group, what, got):
        expected = self._info[group]
        for name, arr in got.items():
            if name not in expected:
                raise WriteError(
                    f"meshio++: vtkhdf: step {self._num_steps} introduces {what} '{name}', "
                    "which the first step did not have; the set of array names is fixed "
                    "at the first write_data"
                )
            if (arr.dtype, _row_shape(arr)) != (
                expected[name][0],
                tuple(expected[name][1]),
            ):
                raise WriteError(
                    f"meshio++: vtkhdf: step {self._num_steps} changes the dtype or "
                    f"component count of {what} '{name}'"
                )
        for name in expected:
            if name not in got:
                raise WriteError(
                    f"meshio++: vtkhdf: step {self._num_steps} is missing {what} '{name}', "
                    "which the first step had"
                )

    def _create_tables(self, point, cell, field):
        root = self._root
        steps = root.create_group("Steps")
        steps.attrs["NSteps"] = np.int64(0)

        def table(grp, name, dtype="i8", row=()):
            grp.create_dataset(
                name,
                shape=(0,) + row,
                maxshape=(None,) + row,
                chunks=(_BOOKKEEPING_ROWS,) + row,
                dtype=dtype,
            )

        table(steps, "Values", "f8")
        for name in _ZERO_TABLES + (("NumberOfParts",)):
            table(steps, name)
        if self._has_poly:
            for name in _POLY_TABLES:
                table(steps, name)
        for grp in (
            "PointDataOffsets",
            "CellDataOffsets",
            "FieldDataOffsets",
            "FieldDataSizes",
        ):
            steps.create_group(grp)
        for group, arrays, rows in (
            ("PointData", point, self._num_points),
            ("CellData", cell, self._num_cells),
            ("FieldData", field, None),
        ):
            g = root[group]
            for name, arr in arrays.items():
                n = rows if rows is not None else arr.shape[0]
                kw = {}
                if self._gzip is not None:
                    kw = {"compression": "gzip", "compression_opts": self._gzip}
                g.create_dataset(
                    name,
                    shape=(0,) + _row_shape(arr),
                    maxshape=(None,) + _row_shape(arr),
                    chunks=(_chunk_rows(n, arr),) + _row_shape(arr),
                    dtype=arr.dtype,
                    **kw,
                )
                table(steps[f"{group}Offsets"], name)
                if group == "FieldData":
                    table(steps["FieldDataSizes"], name, row=(2,))
                    self._field_tuples[name] = 0
                self._info[group][name] = (arr.dtype, _row_shape(arr))
        self._fixed = True

    @staticmethod
    def _append(ds, rows):
        n = ds.shape[0]
        ds.resize(n + rows.shape[0], axis=0)
        ds[n:] = rows

    def write_data(self, t, point_data=None, cell_data=None, field_data=None) -> None:
        """Append one step at simulation time ``t``.

        ``cell_data`` values may be one array or a per-block list (concatenated in
        block order); the reserved ``meshio:time`` field is skipped, since ``t``
        already records it.
        """
        self._require_open()
        if not self._geometry:
            raise WriteError(
                "meshio++: vtkhdf: write_points_cells must be called before write_data"
            )
        point = {
            n: self._prep(n, a, "point_data", self._num_points)
            for n, a in sorted((point_data or {}).items())
        }
        cell = {}
        for n, a in sorted((cell_data or {}).items()):
            if isinstance(a, (list, tuple)):
                a = (
                    np.concatenate([np.asarray(x) for x in a])
                    if len(a)
                    else np.zeros(0)
                )
            cell[n] = self._prep(n, a, "cell_data", self._num_cells)
        field = {
            n: self._prep(n, a, "field_data")
            for n, a in sorted((field_data or {}).items())
            if n != _v.TIME_KEY
        }
        if self._fixed:
            self._check_set("PointData", "point_data", point)
            self._check_set("CellData", "cell_data", cell)
            self._check_set("FieldData", "field_data", field)
        else:
            self._create_tables(point, cell, field)
        root = self._root
        steps = root["Steps"]
        step = self._num_steps
        one = lambda v, dt="i8": np.array([v], dtype=dt)  # noqa: E731
        self._append(steps["Values"], one(t, "f8"))
        for name in _ZERO_TABLES:
            self._append(steps[name], one(0))
        self._append(steps["NumberOfParts"], one(1))
        if self._has_poly:
            for name in _POLY_TABLES:
                self._append(steps[name], one(0))
        for group, arrays, rows in (
            ("PointData", point, self._num_points),
            ("CellData", cell, self._num_cells),
        ):
            for name, arr in arrays.items():
                self._append(root[group][name], arr)
                self._append(steps[f"{group}Offsets"][name], one(step * rows))
        for name, arr in field.items():
            ncomp = arr.shape[1] if arr.ndim == 2 else 1
            self._append(root["FieldData"][name], arr)
            self._append(steps["FieldDataOffsets"][name], one(self._field_tuples[name]))
            self._append(
                steps["FieldDataSizes"][name],
                np.array([[ncomp, arr.shape[0]]], dtype="i8"),
            )
            self._field_tuples[name] += arr.shape[0]
        self._num_steps += 1
        steps.attrs["NSteps"] = np.int64(self._num_steps)
        if self.auto_flush:
            self._file.flush()

    def flush(self) -> None:
        """Make everything written so far durable. Cheap; a no-op once finalized."""
        if self._file is not None and not self.finalized:
            self._file.flush()

    def finalize(self) -> None:
        """Flush and close the file. Idempotent; ``__exit__`` calls it."""
        if self.finalized:
            return
        self.finalized = True
        if self._file is not None:
            self._file.close()
            self._file = None


class TimeSeriesReader:
    """Random access to the steps of a transient VTKHDF file.

    ``num_steps`` and ``times`` come from ``Steps/Values`` without touching the
    geometry; ``read_points_cells()`` and ``read_data(k)`` mirror the XDMF twin.
    """

    def __init__(self, filename):
        import h5py

        self.filename = str(filename)
        try:
            f = h5py.File(self.filename, "r")
        except OSError as exc:
            raise ReadError(
                f"meshio++: vtkhdf: cannot open '{filename}' as an HDF5 file ({exc})"
            ) from exc
        with f:
            if _v.ROOT not in f:
                raise ReadError(
                    f"meshio++: vtkhdf: '{filename}' has no /{_v.ROOT} group"
                )
            steps = _v._Steps(f[_v.ROOT])
        self.num_steps = steps.count
        self.times = [] if steps.values is None else [float(x) for x in steps.values]

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def read_points_cells(self):
        mesh = _v.read(self.filename, points_only=True)
        return mesh.points, mesh.cells

    def read_data(self, k):
        """``(time, point_data, cell_data)`` of step ``k`` (negative counts from the end)."""
        mesh = _v.read(self.filename, time_step=k)
        t = mesh.field_data.get(_v.TIME_KEY)
        return (
            None if t is None else float(np.asarray(t).ravel()[0]),
            mesh.point_data,
            mesh.cell_data,
        )
