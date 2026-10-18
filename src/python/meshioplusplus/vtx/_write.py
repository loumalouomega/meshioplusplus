"""Write DOLFINx VTX output: an ADIOS2 ``.bp`` directory, through ``adios2``.

The Python twin of ``write_vtx`` in ``src/cpp/src/formats/vtx.cpp``: both write
the layout ``dolfinx.io.VTXWriter`` does and ``_vtx.py`` (and ParaView's
``ADIOS2VTXReader``) read. One rank: a ``vtk.xml`` schema attribute and, per
step, ``step`` (the time), ``geometry`` (points x 3), ``connectivity`` (cells x
(1 + nodes), each row prefixed with its node count), ``types``, the sizes
``NumberOfNodes`` and ``NumberOfEntities`` and one variable per point or cell
array. The mesh is written at every step (``VTXMeshPolicy.update``), or at the
first only (``VTXMeshPolicy.reuse``, ``SeriesWriter(mesh_policy="reuse")``).

A VTX file written by DOLFINx holds one cell type on one block per MPI rank. A mesh
with several cell types, or one whose cells carry ``partition:part`` values of
several parts, is written as one block per cell type (and part), every block with
its own points: a point used by several blocks is a ghost copy in all but the first,
marked by ``vtkGhostType`` and ``vtkOriginalPointIds`` and by the ``meshio.shared_points``
attribute, which tells the reader to weld the copies back into the one mesh.
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
PARTITION_KEY = "partition:part"
GHOST_TYPE = "vtkGhostType"
ORIGINAL_IDS = "vtkOriginalPointIds"
SHARED_ATTRIBUTE = "meshio.shared_points"
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


class _Group:
    """One block of the file: a cell type (of one part) and the points it uses."""

    def __init__(self, vtk_type, connectivity, cells):
        self.vtk_type = vtk_type
        self.connectivity = connectivity  # (cells, 1 + nodes), global point numbers
        self.cells = cells  # the cells' numbers in the mesh, None for all in order
        self.points = None  # the mesh points the block holds, None for all
        self.ghost = None  # uint8 per block point: 1 where another block owns it


class _Model:
    """The mesh, as VTXWriter lays it out, and the layout of its arrays."""

    def __init__(self, mesh):
        points = np.asarray(mesh.points, dtype=np.float64)
        n = len(points)
        points = points.reshape(n, -1)
        self.geometry = np.zeros((n, 3))
        self.geometry[:, : min(3, points.shape[1])] = points[:, :3]
        blocks, offset = [], 0
        for k, block in enumerate(mesh.cells):
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
            rows = np.hstack(
                [np.full((len(data), 1), data.shape[1], dtype=np.int64), data]
            )
            blocks.append(
                (k, block.type, _LAGRANGE.get((vtk, data.shape[1]), vtk), rows, offset)
            )
            offset += len(data)
        self.ncells = offset
        if not self.ncells:
            _fail("the mesh has no cells")
        parts = self._parts(mesh, [b[0] for b in blocks])
        labels = sorted({int(v) for p in (parts or []) for v in np.unique(p)})
        if len(labels) > 1:
            # one block per part and cell block: a rank's cells, as DOLFINx writes them
            self.groups = []
            for label in labels:
                for (_k, _t, vtk, rows, first), part in zip(blocks, parts):
                    mask = part == label
                    if mask.any():
                        self.groups.append(
                            _Group(vtk, rows[mask], first + np.flatnonzero(mask))
                        )
        elif len({b[1] for b in blocks}) > 1:
            self.groups = [
                _Group(vtk, rows, first + np.arange(len(rows)))
                for (_k, _t, vtk, rows, first) in blocks
            ]
        else:
            self.groups = [
                _Group(
                    blocks[-1][2],
                    np.concatenate([b[3] for b in blocks]),
                    None,
                )
            ]
        self.shared = len(self.groups) > 1
        if self.shared:
            self._share_points(n)
        reserved = (GHOST_TYPE, ORIGINAL_IDS)
        self.point_names, self.cell_names = [], []
        for name in sorted(mesh.point_data):
            if self.shared and name in reserved:
                warn(
                    f"VTX (.bp): point data '{name}' is replaced by the ghost array of "
                    "the blocks"
                )
            elif self._usable(name, mesh.point_data[name], n, "point"):
                self.point_names.append(name)
        if self.shared:
            self.point_names = sorted(self.point_names + list(reserved))
        for name in sorted(mesh.cell_data):
            arrays = [np.asarray(a) for a in mesh.cell_data[name]]
            if self._usable(
                name,
                np.concatenate(arrays) if arrays else None,
                self.ncells,
                "cell",
            ):
                self.cell_names.append(name)
        self.signature = (
            n,
            tuple(
                (
                    g.vtk_type,
                    g.connectivity.shape,
                    None if g.cells is None else len(g.cells),
                )
                for g in self.groups
            ),
            tuple(self.point_names),
            tuple(self.cell_names),
        )

    @staticmethod
    def _parts(mesh, indices):
        """The ``partition:part`` values of every non-empty cell block, or None."""
        data = getattr(mesh, "cell_data", {}) or {}
        if PARTITION_KEY not in data:
            return None
        arrays = list(data[PARTITION_KEY])
        out = []
        for k in indices:
            rows = len(mesh.cells[k].data)
            if k >= len(arrays):
                return None
            a = np.asarray(arrays[k])
            if a.ndim != 1 or len(a) != rows or a.dtype.kind not in "iub":
                return None
            out.append(a)
        return out

    def _share_points(self, n):
        """Give every block the points it uses: the first block to use a point owns
        it, the others hold a ghost copy; points no cell uses go to the first."""
        owner = np.full(n, -1, dtype=np.int64)
        used = []
        for g, group in enumerate(self.groups):
            nodes = np.unique(group.connectivity[:, 1:])
            used.append(nodes)
            fresh = nodes[owner[nodes] < 0]
            owner[fresh] = g
        owner[owner < 0] = 0
        for g, group in enumerate(self.groups):
            points = np.union1d(used[g], np.flatnonzero(owner == g))
            group.points = points.astype(np.int64)
            group.ghost = (owner[points] != g).astype(np.uint8)
            rows = group.connectivity.copy()
            rows[:, 1:] = np.searchsorted(points, rows[:, 1:])
            group.connectivity = rows

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

    def same_mesh(self, other):
        """True when ``other`` has this model's points and cells exactly."""
        return (
            np.array_equal(self.geometry, other.geometry)
            and len(self.groups) == len(other.groups)
            and all(
                a.vtk_type == b.vtk_type
                and np.array_equal(a.connectivity, b.connectivity)
                and np.array_equal(a.points, b.points)
                for a, b in zip(self.groups, other.groups)
            )
        )

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

    def write_attributes(self, stream):
        stream.write_attribute("vtk.xml", self.schema())
        if self.shared:
            stream.write_attribute(SHARED_ATTRIBUTE, ORIGINAL_IDS)

    def _point_block(self, mesh, name, group):
        if name == GHOST_TYPE and self.shared:
            return group.ghost
        if name == ORIGINAL_IDS and self.shared:
            return group.points
        values = _storable(mesh.point_data[name])
        return values if group.points is None else values[group.points]

    def put_step(self, stream, mesh, time, with_mesh=True):
        """One step. ``with_mesh=False`` is the step of a ``reuse`` file after the
        first: the time and the arrays only."""
        stream.begin_step()
        stream.write("step", float(time))
        if with_mesh:
            for name, count in (
                (
                    "NumberOfNodes",
                    [
                        len(g.points) if g.points is not None else len(self.geometry)
                        for g in self.groups
                    ],
                ),
                ("NumberOfEntities", [len(g.connectivity) for g in self.groups]),
            ):
                for value in count:
                    array = np.array([value], np.uint32)
                    if self.shared:
                        _local(stream, name, array)
                    else:
                        stream.write(name, array, [1], [0], [1])
            for group in self.groups:
                geometry = (
                    self.geometry
                    if group.points is None
                    else self.geometry[group.points]
                )
                _local(stream, "geometry", geometry)
            for group in self.groups:
                _local(stream, "connectivity", group.connectivity)
            if self.shared:
                for group in self.groups:
                    _local(stream, "types", np.array([group.vtk_type], np.uint32))
            else:
                stream.write(
                    "types", np.array(self.groups[0].vtk_type, dtype=np.uint32)
                )
        for name in self.point_names:
            if not with_mesh and name in (GHOST_TYPE, ORIGINAL_IDS) and self.shared:
                continue
            for group in self.groups:
                _local(stream, name, self._point_block(mesh, name, group))
        for name in self.cell_names:
            values = _storable(
                np.concatenate([np.asarray(a) for a in mesh.cell_data[name]])
            )
            for group in self.groups:
                _local(
                    stream,
                    name,
                    values if group.cells is None else values[group.cells],
                )
        stream.end_step()


def _local(stream, name, values):
    """Write a rank-local array (no global shape), as VTXWriter does; a variable
    written again in the same step is a further block."""
    io = stream._io
    variable = io.inquire_variable(name)
    if not variable:
        variable = io.define_variable(name, values, [], [], list(values.shape))
    else:
        variable.set_selection([[], list(values.shape)])
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
        model.write_attributes(stream)
        model.put_step(stream, mesh, _time_of(mesh))


class SeriesWriter:
    """A VTX file with a step per mesh; every step must have the first step's
    points, cells and array names.

    ``mesh_policy="update"`` (the default) writes the mesh at every step;
    ``"reuse"`` writes it at the first step only (DOLFINx's ``VTXMeshPolicy.reuse``),
    so every step must also have the first step's points and connectivity."""

    def __init__(self, filename, mesh_policy="update"):
        if mesh_policy not in ("update", "reuse"):
            raise ValueError(
                f"meshio++: mesh_policy must be 'update' or 'reuse', got '{mesh_policy}'"
            )
        if is_buffer(filename, "w"):
            _fail("cannot be written to a stream")
        self._reuse = mesh_policy == "reuse"
        self._stream = _open(filename)
        self._model = None
        self._steps = 0

    def write(self, time, mesh):
        model = _Model(mesh)
        if self._model is None:
            self._model = model
            model.write_attributes(self._stream)
        elif model.signature != self._model.signature:
            _fail(
                f"series: step {self._steps}'s points, cells or arrays differ from "
                "the first step's; write one file per step with '{step}'"
            )
        elif self._reuse and not self._model.same_mesh(model):
            _fail(
                f"series: step {self._steps}'s points or cells differ from the first "
                "step's, which the reuse policy cannot write; use mesh_policy='update'"
            )
        model.put_step(
            self._stream,
            mesh,
            float(time),
            with_mesh=not (self._reuse and self._steps > 0),
        )
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
