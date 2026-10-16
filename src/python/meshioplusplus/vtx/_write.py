"""Write DOLFINx VTX output: an ADIOS2 ``.bp`` directory, through ``adios2``.

The Python twin of ``write_vtx`` in ``src/cpp/src/formats/vtx.cpp``: both write
the layout ``dolfinx.io.VTXWriter`` does and ``_vtx.py`` (and ParaView's
``ADIOS2VTXReader``) read. One rank: a ``vtk.xml`` schema attribute and, per
step, ``step`` (the time), ``geometry`` (points x 3), ``connectivity`` (cells x
(1 + nodes), each row prefixed with its node count), ``types``, the sizes
``NumberOfNodes`` and ``NumberOfEntities`` and one variable per point or cell
array. The mesh is written at every step (``VTXMeshPolicy.update``).

A VTX file holds one cell type, as DOLFINx's does: a mesh with several raises.
"""

from __future__ import annotations

from xml.sax.saxutils import quoteattr

import numpy as np

from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from .._vtk_common import meshio_to_vtk_order, meshio_to_vtk_type

__all__ = ["write", "SeriesWriter"]

TIME_KEY = "meshio:time"
# (VTK type, node count) -> VTK_LAGRANGE_* type, the inverse of the reader's
# lowering: DOLFINx writes Lagrange cells of every degree.
_LAGRANGE = {
    (3, 2): 68,
    (21, 3): 68,
    (5, 3): 69,
    (22, 6): 69,
    (9, 4): 70,
    (28, 9): 70,
    (10, 4): 71,
    (24, 10): 71,
    (12, 8): 72,
    (13, 6): 73,
    (14, 5): 74,
}
_RESERVED = {
    "geometry",
    "connectivity",
    "types",
    "step",
    "TIME",
    "NumberOfNodes",
    "NumberOfEntities",
}


def _adios2():
    from .._interop import _require

    return _require("adios2", "adios2", "writing a DOLFINx VTX (.bp) file")


def _fail(message):
    raise WriteError(f"VTX (.bp): {message}")


class _Model:
    """The mesh, as VTXWriter lays it out, and the layout of its arrays."""

    def __init__(self, mesh):
        points = np.asarray(mesh.points, dtype=np.float64)
        n = len(points)
        points = points.reshape(n, -1)
        self.geometry = np.zeros((n, 3))
        self.geometry[:, : min(3, points.shape[1])] = points[:, :3]
        types = {block.type for block in mesh.cells if len(block.data)}
        if len(types) > 1:
            _fail(
                "a VTX file holds one cell type, as DOLFINx's does; this mesh has "
                + ", ".join(sorted(types))
                + " (write each type to its own file)"
            )
        self.rows = []
        self.vtk_type = 0
        for block in mesh.cells:
            data = np.asarray(block.data, dtype=np.int64)
            if not len(data):
                continue
            if block.type not in meshio_to_vtk_type:
                _fail(f"cell type '{block.type}' has no VTK type")
            if data.min() < 0 or data.max() >= n:
                _fail(f"a '{block.type}' cell refers to a point the mesh does not have")
            order = meshio_to_vtk_order(block.type)
            if order is not None and data.shape[1] == len(order):
                data = data[:, order]
            vtk = meshio_to_vtk_type[block.type]
            self.vtk_type = _LAGRANGE.get((vtk, data.shape[1]), vtk)
            self.rows.append(
                np.hstack(
                    [np.full((len(data), 1), data.shape[1], dtype=np.int64), data]
                )
            )
        self.connectivity = (
            np.concatenate(self.rows) if self.rows else np.zeros((0, 1), np.int64)
        )
        self.ncells = len(self.connectivity)
        if not self.ncells:
            _fail("the mesh has no cells")
        self.point_names, self.cell_names = [], []
        self.signature = (
            n,
            tuple((b.type, np.asarray(b.data).shape) for b in mesh.cells),
        )
        for name in sorted(mesh.point_data):
            if self._usable(name, mesh.point_data[name], n, "point"):
                self.point_names.append(name)
        for name in sorted(mesh.cell_data):
            parts = [np.asarray(a) for a in mesh.cell_data[name]]
            if self._usable(
                name, np.concatenate(parts) if parts else None, self.ncells, "cell"
            ):
                self.cell_names.append(name)
        self.signature += (tuple(self.point_names), tuple(self.cell_names))

    @staticmethod
    def _usable(name, values, rows, what):
        if name in _RESERVED:
            warn(
                f"VTX (.bp): {what} data '{name}' is a reserved name and is not written"
            )
            return False
        a = np.asarray(values)
        if (
            a.dtype.kind not in "fiub"
            or a.ndim == 0
            or a.shape[0] != rows
            or a.ndim > 2
        ):
            warn(f"VTX (.bp): {what} data '{name}' is not written")
            return False
        return True

    def schema(self):
        def arrays(names):
            return "".join(
                f"        <DataArray Name={quoteattr(name)} />\n" for name in names
            )

        point = '        <DataArray Name="TIME">step</DataArray>\n' + arrays(
            self.point_names
        )
        cell = (
            f"      <CellData>\n{arrays(self.cell_names)}      </CellData>\n"
            if self.cell_names
            else ""
        )
        return (
            '<?xml version="1.0"?>\n'
            '<VTKFile type="UnstructuredGrid" version="0.1">\n'
            "  <UnstructuredGrid>\n"
            '    <Piece NumberOfPoints="NumberOfNodes" NumberOfCells="NumberOfEntities">\n'
            "      <Points>\n"
            '        <DataArray Name="geometry" />\n'
            "      </Points>\n"
            "      <Cells>\n"
            '        <DataArray Name="connectivity" />\n'
            '        <DataArray Name="types" />\n'
            "      </Cells>\n"
            f"      <PointData>\n{point}      </PointData>\n"
            f"{cell}"
            "    </Piece>\n"
            "  </UnstructuredGrid>\n"
            "</VTKFile>\n"
        )

    def put_step(self, stream, mesh, time):
        stream.begin_step()
        stream.write("step", float(time))
        stream.write(
            "NumberOfNodes", np.array([len(self.geometry)], np.uint32), [1], [0], [1]
        )
        stream.write(
            "NumberOfEntities", np.array([self.ncells], np.uint32), [1], [0], [1]
        )
        _local(stream, "geometry", self.geometry)
        _local(stream, "connectivity", self.connectivity)
        stream.write("types", np.array(self.vtk_type, dtype=np.uint32))
        for name in self.point_names:
            _local(stream, name, _storable(mesh.point_data[name]))
        for name in self.cell_names:
            values = np.concatenate([np.asarray(a) for a in mesh.cell_data[name]])
            _local(stream, name, _storable(values))
        stream.end_step()


def _local(stream, name, values):
    """Write a rank-local array (no global shape), as VTXWriter does."""
    io = stream._io
    variable = io.inquire_variable(name)
    if not variable:
        variable = io.define_variable(name, values, [], [], list(values.shape))
    stream.write(variable, values)


def _storable(values):
    a = np.ascontiguousarray(values)
    if a.dtype == np.bool_:
        a = a.astype(np.uint8)
    return a


def _time_of(mesh):
    fd = getattr(mesh, "field_data", {}) or {}
    if TIME_KEY in fd:
        return float(np.asarray(fd[TIME_KEY]).reshape(-1)[0])
    return 0.0


def _open(filename):
    adios2 = _adios2()
    return adios2.Stream(str(filename), "w")


def write(filename, mesh):
    """Write one step of a DOLFINx VTX file (an ADIOS2 ``.bp`` directory)."""
    if is_buffer(filename, "w"):
        _fail("cannot be written to a stream")
    model = _Model(mesh)
    with _open(filename) as stream:
        stream.write_attribute("vtk.xml", model.schema())
        model.put_step(stream, mesh, _time_of(mesh))


class SeriesWriter:
    """A VTX file with a step per mesh; every step must have the first step's
    points, cells and array names."""

    def __init__(self, filename):
        if is_buffer(filename, "w"):
            _fail("cannot be written to a stream")
        self._stream = _open(filename)
        self._model = None
        self._steps = 0

    def write(self, time, mesh):
        model = _Model(mesh)
        if self._model is None:
            self._model = model
            self._stream.write_attribute("vtk.xml", model.schema())
        elif model.signature != self._model.signature:
            _fail(
                f"series: step {self._steps}'s points, cells or arrays differ from "
                "the first step's; write one file per step with '{step}'"
            )
        model.put_step(self._stream, mesh, float(time))
        self._steps += 1

    def close(self):
        if self._stream is not None:
            self._stream.close()
            self._stream = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
