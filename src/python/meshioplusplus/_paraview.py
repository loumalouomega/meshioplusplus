"""ParaView bridge: the conversions behind the ParaView plugin.

The plugin (``tools/paraview-meshioplusplus-plugin.py``, installed with the
wheel under ``share/meshioplusplus/paraview/``) declares ParaView proxies and
nothing else; everything it does to a mesh lives here, so it can be tested
without ParaView. See ``doc/paraview_plugin.md``.

Architecture
------------
Split the way :mod:`meshioplusplus._blender` is, in two layers:

* a **payload** layer that imports neither VTK nor ParaView. The mesh-to-VTK
  direction is :func:`meshioplusplus._interop._to_vtk_payload`, the one
  ``to_pyvista`` uses (vectorised, every ``polygon*`` block mapped to
  ``VTK_POLYGON``, points of any width padded to three columns); the
  VTK-to-mesh direction is
  :func:`meshioplusplus._vtk_common.vtk_cells_from_data`, the VTU reader's
  (cells split into runs in their original order). Node orders need no
  permutation: meshio++'s are VTK's. Neither is copied here.
* a **VTK adapter** layer, :func:`grid_from_mesh` and :func:`mesh_from_grid`,
  which imports ``vtkmodules`` inside each function and moves whole arrays
  through ``numpy_support``: no per-cell Python loop in either direction. It
  needs VTK but not ParaView, so the PyPI ``vtk`` wheel tests it.

Notes about what could not be carried (a cell block VTK has no type for, a
string array) are surfaced as warnings, the interop layer's convention.
"""

from __future__ import annotations

import json

import numpy as np

from ._interop import REGION_META_KEY, REGION_PREFIX

#: The arrays VTK reads a Lagrange or Bézier cell's degrees and a Bézier
#: cell's weights from, once set as the active attribute of that name (the
#: names the VTU readers and writers use).
HIGHER_ORDER_DEGREES = "HigherOrderDegrees"
RATIONAL_WEIGHTS = "RationalWeights"

#: Operation names the notes are reported under.
READER_OP = "ParaView reader"
WRITER_OP = "ParaView writer"


# --------------------------------------------------------------------------- #
# VTK adapter: mesh -> vtkUnstructuredGrid                                    #
# --------------------------------------------------------------------------- #
def _vtk_array(values, name=None):
    """A VTK data array sharing ``values`` (kept alive by the array)."""
    from vtkmodules.util.numpy_support import numpy_to_vtk

    values = np.ascontiguousarray(values)
    array = numpy_to_vtk(values, deep=False)
    # numpy_to_vtk attaches the numpy array to the VTK object, which holds the
    # buffer for as long as VTK holds the array.
    if name is not None:
        array.SetName(name)
    return array


def _id_array(values):
    """A ``vtkIdTypeArray`` over int64 ``values``, shared."""
    from vtkmodules.util.numpy_support import numpy_to_vtkIdTypeArray

    return numpy_to_vtkIdTypeArray(
        np.ascontiguousarray(values, dtype=np.int64), deep=False
    )


def _field_array(name, raw, notes):
    """One field-data array: numbers as a data array, strings as a string array."""
    a = np.asarray(raw)
    if a.dtype.kind in "SU" or a.dtype == object:
        from vtkmodules.vtkCommonCore import vtkStringArray

        strings = vtkStringArray()
        strings.SetName(name)
        for value in a.reshape(-1):
            strings.InsertNextValue(str(value))
        return strings
    if np.issubdtype(a.dtype, np.complexfloating):
        notes.append(f"field_data:{name} has dtype {a.dtype} and was dropped")
        return None
    if a.dtype == bool:
        a = a.astype(np.uint8)
    if a.ndim == 0:
        a = a.reshape(1)
    return _vtk_array(a, name)


def legacy_face_stream(face_connectivity, face_offsets, face_location_offsets):
    """Polyhedron faces in VTK's pre-9.4 layout: ``(faces, face_locations)``.

    ``faces`` holds, per polyhedron, ``nfaces, n0, ids..., n1, ids...``;
    ``face_locations[c]`` is where cell ``c``'s entry starts, ``-1`` for a cell
    that is not a polyhedron. Built from the payload's cell-array layout
    without a loop per cell.
    """
    from ._ragged import ragged_gather

    loc = np.asarray(face_location_offsets, dtype=np.int64)
    per_cell = np.diff(loc)
    poly = per_cell > 0
    sizes = np.diff(np.asarray(face_offsets, dtype=np.int64))
    tokens = np.zeros(len(sizes) + 1, dtype=np.int64)
    np.cumsum(1 + sizes, out=tokens[1:])
    before = np.cumsum(poly) - poly  # polyhedra ahead of each cell
    start = before + tokens[loc[:-1]]
    out = np.empty(int(poly.sum()) + int(tokens[-1]), dtype=np.int64)
    out[start[poly]] = per_cell[poly]
    cell_of_face = np.repeat(np.arange(len(per_cell), dtype=np.int64), per_cell)
    face_at = before[cell_of_face] + 1 + tokens[:-1]
    out[face_at] = sizes
    out[ragged_gather(face_at + 1, sizes)[0]] = face_connectivity
    return out, np.where(poly, start, -1).astype(np.int64)


def faces_from_legacy_stream(faces, face_locations):
    """The inverse of :func:`legacy_face_stream`: ``(face_connectivity,
    face_offsets, face_location_offsets)``. One loop over the polyhedra, which
    only a VTK older than 9.4 needs."""
    faces = np.asarray(faces, dtype=np.int64)
    where = np.asarray(face_locations, dtype=np.int64)
    conn, sizes, per_cell = [], [], np.zeros(len(where), dtype=np.int64)
    for c in np.flatnonzero(where >= 0):
        at = int(where[c])
        per_cell[c] = faces[at]
        at += 1
        for _ in range(int(faces[at - 1])):
            n = int(faces[at])
            conn.append(faces[at + 1 : at + 1 + n])
            sizes.append(n)
            at += 1 + n
    face_offsets = np.zeros(len(sizes) + 1, dtype=np.int64)
    np.cumsum(sizes, out=face_offsets[1:])
    locations = np.zeros(len(where) + 1, dtype=np.int64)
    np.cumsum(per_cell, out=locations[1:])
    flat = np.concatenate(conn) if conn else np.empty(0, dtype=np.int64)
    return flat, face_offsets, locations


def set_cells(grid, payload, legacy=None):
    """Set ``payload``'s cells on the ``vtkUnstructuredGrid`` ``grid``.

    Polyhedra go through ``SetPolyhedralCells`` (VTK 9.4 and later) or, on an
    older VTK or with ``legacy=True``, the face stream ``SetCells`` took
    before. ``to_pyvista`` sets its cells here too.
    """
    from vtkmodules.util.numpy_support import numpy_to_vtk
    from vtkmodules.vtkCommonCore import VTK_UNSIGNED_CHAR
    from vtkmodules.vtkCommonDataModel import vtkCellArray

    cells = vtkCellArray()
    cells.SetData(_id_array(payload.offsets), _id_array(payload.connectivity))
    types = numpy_to_vtk(
        np.ascontiguousarray(payload.celltypes, dtype=np.uint8),
        deep=False,
        array_type=VTK_UNSIGNED_CHAR,
    )
    if not payload.has_polyhedra:
        grid.SetCells(types, cells)
        return
    if legacy is None:
        legacy = not hasattr(grid, "SetPolyhedralCells")
    if legacy:
        stream, where = legacy_face_stream(
            payload.face_connectivity,
            payload.face_offsets,
            payload.face_location_offsets,
        )
        grid.SetCells(types, cells, _id_array(where), _id_array(stream))
        return
    faces = vtkCellArray()
    faces.SetData(_id_array(payload.face_offsets), _id_array(payload.face_connectivity))
    locations = vtkCellArray()
    locations.SetData(
        _id_array(payload.face_location_offsets),
        _id_array(np.arange(len(payload.face_offsets) - 1, dtype=np.int64)),
    )
    grid.SetPolyhedralCells(types, cells, locations, faces)


def fill_grid(grid, mesh, op=READER_OP):
    """Fill the ``vtkUnstructuredGrid`` ``grid`` from ``mesh``, in place.

    The plugin's reader passes the grid ParaView hands it; :func:`grid_from_mesh`
    passes a new one. Returns the payload's notes, already emitted.
    """
    from vtkmodules.vtkCommonCore import vtkPoints

    from ._interop import _emit, _to_vtk_payload

    payload = _to_vtk_payload(mesh, carry_sides=True)
    notes = list(payload.notes)

    points = vtkPoints()
    points.SetData(_vtk_array(payload.points))
    grid.SetPoints(points)

    set_cells(grid, payload)

    for name, values in payload.point_data.items():
        array = _vtk_array(values, name)
        grid.GetPointData().AddArray(array)
        if name == RATIONAL_WEIGHTS:
            grid.GetPointData().SetRationalWeights(array)
    for name, values in payload.cell_data.items():
        array = _vtk_array(values, name)
        grid.GetCellData().AddArray(array)
        if name == HIGHER_ORDER_DEGREES:
            grid.GetCellData().SetHigherOrderDegrees(array)
    for (kind, name), mask in payload.region_masks.items():
        target = grid.GetPointData() if kind == "point" else grid.GetCellData()
        target.AddArray(_vtk_array(mask, REGION_PREFIX + name))
    if payload.region_meta:
        meta = json.dumps(payload.region_meta, separators=(",", ":"))
        grid.GetFieldData().AddArray(_field_array(REGION_META_KEY, [meta], notes))
    for name, entries in side_region_arrays(mesh, payload.kept_global, notes).items():
        grid.GetFieldData().AddArray(_vtk_array(entries, name))
    for name in sorted(mesh.field_data):
        array = _field_array(name, mesh.field_data[name], notes)
        if array is not None:
            grid.GetFieldData().AddArray(array)

    _emit(op, notes)
    return notes


def side_region_arrays(mesh, kept_global, notes):
    """Side regions as the VTU convention's field arrays
    (``region:side:<name>``, ``(cell, facet)`` pairs in the grid's cell
    numbering, and ``region-meta:side:<name>``), which ParaView carries as
    field data and the writer reads back.

    ``kept_global`` maps a grid cell to its global cell; an entry in a cell
    block the grid does not hold is dropped with a note.
    """
    from types import SimpleNamespace

    from ._region_field_data import regions_to_field_arrays
    from ._regions import Region

    sides = [r for r in getattr(mesh, "regions", []) if r.kind == "side"]
    if not sides:
        return {}
    kept_global = np.asarray(kept_global, dtype=np.int64)
    total = int(kept_global.max()) + 1 if len(kept_global) else 0
    total = max(
        [total]
        + [
            int(np.asarray(r.entries).reshape(-1, 2)[:, 0].max()) + 1
            for r in sides
            if np.asarray(r.entries).size
        ]
    )
    to_grid = np.full(total, -1, dtype=np.int64)
    to_grid[kept_global] = np.arange(len(kept_global), dtype=np.int64)

    carried = []
    for r in sides:
        e = np.asarray(r.entries, dtype=np.int64).reshape(-1, 2)
        keep = to_grid[e[:, 0]] >= 0 if len(e) else np.zeros(0, dtype=bool)
        if not bool(keep.all()):
            notes.append(
                f"side region '{r.name}': {int((~keep).sum())} entrie(s) lie in "
                "cell blocks VTK cannot express"
            )
        carried.append(Region(r.name, "side", e[keep], r.dim, r.tag))
    # regions_to_field_arrays reads only a mesh's regions.
    return regions_to_field_arrays(SimpleNamespace(regions=carried), to_grid)


def grid_from_mesh(mesh, op=READER_OP):
    """A new ``vtkUnstructuredGrid`` holding ``mesh``."""
    from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid

    grid = vtkUnstructuredGrid()
    fill_grid(grid, mesh, op)
    return grid


# --------------------------------------------------------------------------- #
# VTK adapter: vtkUnstructuredGrid -> mesh                                    #
# --------------------------------------------------------------------------- #
def _numpy(array):
    """A VTK data array as numpy, through a copy into AOS storage if needed.

    ``vtk_to_numpy`` needs a contiguous array-of-structures buffer; VTK can hand
    back other storage (an implicit or affine array, say) for offsets or
    connectivity, which a deep copy into a plain array of the same type
    converts.
    """
    from vtkmodules.util.numpy_support import vtk_to_numpy
    from vtkmodules.vtkCommonCore import vtkDataArray

    try:
        return vtk_to_numpy(array)
    except Exception:
        plain = vtkDataArray.CreateDataArray(array.GetDataType())
        plain.DeepCopy(array)
        return vtk_to_numpy(plain)


def _arrays(attributes, what, notes):
    """``{name: numpy}`` for the numeric arrays of a ``vtkFieldData``.

    String and other non-numeric arrays (the QA records and information lines
    ParaView's Exodus reader adds as field data) have no meshio++ counterpart
    among point and cell data and are skipped with a note; for field data
    a string array is kept as a numpy string array.
    """
    from vtkmodules.vtkCommonCore import vtkDataArray, vtkStringArray

    out = {}
    for i in range(attributes.GetNumberOfArrays()):
        array = attributes.GetAbstractArray(i)
        if array is None:
            continue
        name = array.GetName() or f"{what}_{i}"
        if isinstance(array, vtkDataArray):
            out[name] = _numpy(array)
        elif what == "field_data" and isinstance(array, vtkStringArray):
            out[name] = np.array(
                [array.GetValue(k) for k in range(array.GetNumberOfValues())]
            )
        else:
            notes.append(
                f"{what}:{name} is a {array.GetClassName()} and was skipped "
                "(meshio++ holds numeric arrays only)"
            )
    return out


def _cell_array(cell_array):
    """``(offsets, connectivity)`` of a ``vtkCellArray``, int64."""
    return (
        _numpy(cell_array.GetOffsetsArray()).astype(np.int64, copy=False),
        _numpy(cell_array.GetConnectivityArray()).astype(np.int64, copy=False),
    )


def _polyhedron_faces(grid):
    """``(face_connectivity, face_offsets, location_offsets, location_ids)`` of
    a grid's polyhedra, or ``None``; either VTK API."""
    if hasattr(grid, "GetPolyhedronFaces"):
        faces = grid.GetPolyhedronFaces()
        locations = grid.GetPolyhedronFaceLocations()
        if faces is None or locations is None:
            return None
        face_offsets, face_conn = _cell_array(faces)
        loc_offsets, loc_ids = _cell_array(locations)
        return face_conn, face_offsets, loc_offsets, loc_ids
    stream, where = grid.GetFaces(), grid.GetFaceLocations()
    if stream is None or where is None:
        return None
    flat, offsets, locations = faces_from_legacy_stream(_numpy(stream), _numpy(where))
    return flat, offsets, locations, np.arange(len(offsets) - 1, dtype=np.int64)


def _polyhedron_blocks(first, last, npts, faces, cell_flat):
    """Cell blocks of the polyhedra ``first:last`` (a run of the grid), one per
    run of equal point count, so the cells keep their order (the VTU reader's
    ``polyhedron<n>`` names)."""
    from ._mesh import CellBlock
    from ._ragged import CsrRows, ragged_gather

    face_conn, face_offsets, loc_offsets, loc_ids = faces
    blocks, data = [], {}
    counts = npts[first:last]
    bounds = np.concatenate(
        [[0], np.flatnonzero(counts[1:] != counts[:-1]) + 1, [len(counts)]]
    )
    for a, b in zip(bounds[:-1], bounds[1:]):
        cells = np.arange(first + a, first + b)
        per_cell = loc_offsets[cells + 1] - loc_offsets[cells]
        face_ids = loc_ids[ragged_gather(loc_offsets[cells], per_cell)[0]]
        sizes = face_offsets[face_ids + 1] - face_offsets[face_ids]
        flat = face_conn[ragged_gather(face_offsets[face_ids], sizes)[0]]
        offsets = np.zeros(len(sizes) + 1, dtype=np.int64)
        np.cumsum(sizes, out=offsets[1:])
        cell_faces = np.zeros(len(cells) + 1, dtype=np.int64)
        np.cumsum(per_cell, out=cell_faces[1:])
        blocks.append(
            CellBlock(
                f"polyhedron{counts[a]}", CsrRows.from_csr(flat, offsets, cell_faces)
            )
        )
        for name, values in cell_flat.items():
            data.setdefault(name, []).append(values[cells])
    return blocks, data


def _cell_types_array(grid):
    """The grid's cell types as a data array: ``GetCellTypes()`` from VTK 9.6
    (where ``GetCellTypesArray`` is deprecated), which before took a
    ``vtkCellTypes`` to fill and returned nothing."""
    try:
        types = grid.GetCellTypes()
    except TypeError:
        types = None
    if types is None or not hasattr(types, "GetNumberOfTuples"):
        types = grid.GetCellTypesArray()
    return types


def cells_from_grid(grid, cell_flat):
    """``(cells, cell_data)`` of a ``vtkUnstructuredGrid``, in its cell order.

    ``cell_flat`` holds the grid's cell arrays, one value per cell. A run of
    consecutive cells of one type (and, for polygons, one size; for polyhedra,
    one point count) becomes a cell block, so the blocks' cells, laid end to
    end, are the grid's in its order. ``from_pyvista`` reads its cells
    here too.
    """
    from ._vtk_common import vtk_cells_from_data

    if not grid.GetNumberOfCells():
        return [], {}
    types = _numpy(_cell_types_array(grid)).astype(np.uint8)
    offsets, connectivity = _cell_array(grid.GetCells())
    poly = types == 42  # VTK_POLYHEDRON
    if not poly.any():
        # vtk_cells_from_data takes END offsets (the VTU convention); a VTK 9
        # cell array's offsets start at 0 and have num_cells + 1 entries.
        return vtk_cells_from_data(connectivity, offsets[1:], types, cell_flat)

    faces = _polyhedron_faces(grid)
    npts = np.diff(offsets)
    bounds = np.concatenate(
        [[0], np.flatnonzero(poly[1:] != poly[:-1]) + 1, [len(types)]]
    )
    cells, cell_data = [], {}
    for first, last in zip(bounds[:-1], bounds[1:]):
        part = {name: v[first:last] for name, v in cell_flat.items()}
        if poly[first]:
            blocks, data = _polyhedron_blocks(first, last, npts, faces, cell_flat)
        else:
            lo, hi = offsets[first], offsets[last]
            blocks, data = vtk_cells_from_data(
                connectivity[lo:hi],
                offsets[first + 1 : last + 1] - lo,
                types[first:last],
                part,
            )
        cells += blocks
        for name in cell_flat:
            cell_data.setdefault(name, []).extend(data.get(name, []))
    return cells, cell_data


def mesh_from_grid(grid, op=WRITER_OP):
    """A mesh from a ``vtkUnstructuredGrid``, cells in the grid's order.

    Consecutive cells of one type (and, for polygons, one size) form a cell
    block, so the grid's cell order is kept: a mixed grid becomes several
    blocks rather than one per type.
    """
    from ._interop import _emit
    from ._mesh import Mesh

    notes: list = []
    if grid.GetPoints() is not None and grid.GetNumberOfPoints():
        points = _numpy(grid.GetPoints().GetData())
    else:
        points = np.zeros((0, 3))

    point_data = _arrays(grid.GetPointData(), "point_data", notes)
    cell_flat = _arrays(grid.GetCellData(), "cell_data", notes)
    field_data = _arrays(grid.GetFieldData(), "field_data", notes)
    regions = _regions_from_grid_arrays(
        point_data, cell_flat, field_data, len(points), notes
    )

    # cells_from_grid keeps the grid's cell order, so a grid cell index is the
    # global (block-major) cell index the regions above are numbered in.
    cells, cell_data = cells_from_grid(grid, cell_flat)

    _emit(op, notes)
    return Mesh(
        points,
        cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
        regions=regions,
    )


def _regions_from_grid_arrays(point_data, cell_flat, field_data, num_points, notes):
    """Regions from a grid's arrays, removing the arrays that carried them:
    ``region:<name>`` masks with the ``meshioplusplus:regions`` sidecar (what
    the reader and ``to_pyvista`` write), and the VTU convention's
    ``region:<kind>:<name>`` field arrays (side regions from the reader, and
    every region of a meshio++ VTU that ParaView's own reader opened)."""
    from ._interop import _regions_from_masks
    from ._region_field_data import regions_from_field_arrays

    meta_by_name = {}
    raw = field_data.pop(REGION_META_KEY, None)
    if raw is not None:
        try:
            entries = json.loads(str(np.asarray(raw).reshape(-1)[0]))
            meta_by_name = {(e["kind"], e["name"]): e for e in entries}
        except (ValueError, KeyError, TypeError, IndexError):
            notes.append(f"the '{REGION_META_KEY}' sidecar is unreadable; ignored")
    regions = _regions_from_masks(point_data, "point", meta_by_name, notes)
    regions += _regions_from_masks(cell_flat, "cell", meta_by_name, notes)
    from_fields, rest = regions_from_field_arrays(
        field_data, num_points, None, WRITER_OP
    )
    field_data.clear()
    field_data.update(rest)
    return regions + from_fields


# --------------------------------------------------------------------------- #
# time (roadmap 6.4)                                                          #
# --------------------------------------------------------------------------- #
#: ``{key: plan}``; a key names each file with its size and modification
#: time, so a file rewritten in place (ParaView's *Reload Files*) is planned
#: again and an unchanged one is not, which matters for the formats whose
#: ``read_metadata`` is a full read.
_PLANS: dict = {}


def _plan_key(paths, file_format):
    import os

    stats = []
    for path in paths:
        try:
            st = os.stat(path)
            stats.append((os.path.abspath(path), st.st_mtime_ns, st.st_size))
        except OSError:
            stats.append((str(path), None, None))
    return (file_format, tuple(stats))


def time_plan(paths, file_format=None):
    """The steps of ``paths``: one entry per step, ``{"path", "step", "time",
    "time_source"}`` (:func:`meshioplusplus._sequence.sequence_entries`).

    One path may hold several steps (the formats ``read_metadata`` reports
    ``time_values`` for); several paths are a file series, each file one step
    or more, in the order given. A step's time is the file's own where it has
    one, else its position. Cached per file identity.
    """
    from ._sequence import _entry_time, sequence_entries

    paths = [str(p) for p in paths]
    key = _plan_key(paths, file_format)
    if key not in _PLANS:
        source = paths[0] if len(paths) == 1 else paths
        plan = sequence_entries(source, file_format=file_format)
        if len(plan) > 1:
            # A step whose only time is inside its file (the meshio:time
            # field data write_sequence stores) is read once to find it, as
            # sequence_times does; the cache keeps it to once per file.
            for entry in plan:
                entry["time"] = _entry_time(entry, file_format, "auto")
        _PLANS[key] = plan
    return _PLANS[key]


def forget_plans():
    """Drop every cached plan (the reader's *Refresh*)."""
    _PLANS.clear()


def time_values(plan):
    """The time steps ParaView should see: ``[]`` for a single step."""
    if len(plan) < 2:
        return []
    return [float(e["time"]) for e in plan]


def step_index(times, t):
    """The step shown at time ``t``: the last one at or before it (the first
    before the first), as VTK's own readers choose."""
    if not times or t is None:
        return 0
    i = int(np.searchsorted(np.asarray(times, dtype=float), t, side="right")) - 1
    return min(max(i, 0), len(times) - 1)


def read_step(
    plan, index, file_format=None, deflection=0.0, angular_deflection=0.5, **read_kwargs
):
    """Read step ``index`` of ``plan`` with ``meshioplusplus.read``'s keywords.

    A STEP or IGES file is triangulated at ``deflection`` and
    ``angular_deflection`` (their readers' defaults are these), then filtered
    as ``read`` filters (``points_only``, ``arrays``).
    """
    from ._helpers import _apply_read_filter, read

    entry = plan[index]
    cad = cad_format(entry["path"], file_format)
    if cad is not None:
        mesh = read_cad(entry["path"], cad, deflection, angular_deflection)
        return _apply_read_filter(
            mesh, read_kwargs.get("points_only", False), read_kwargs.get("arrays")
        )
    return read(
        entry["path"],
        file_format=file_format or None,
        time_step=entry["step"],
        **read_kwargs,
    )


def series_path(path, file_format=None):
    """Where the writer puts several steps: ``path`` itself when its format
    holds a series in one file, else ``path`` with ``_{step}`` before its
    extension, one file per step (:func:`~meshioplusplus._sequence.write_sequence`'s
    two modes)."""
    import pathlib

    from ._exceptions import WriteError
    from ._sequence import _check_series_target

    try:
        _check_series_target(path, file_format)
        return str(path)
    except WriteError:
        p = pathlib.Path(str(path))
        return str(p.with_name(p.stem + "_{step}" + p.suffix))


def write_one_step(template, index, count, time, mesh, file_format=None, **kwargs):
    """Write step ``index`` of ``count`` to its own file of ``template`` (a
    ``{step}`` path), with its time, as ``write_sequence``'s fan-out does; the
    plugin's writer calls it as each step arrives, so no step is held."""
    from ._helpers import write
    from ._sequence import TIME_KEY, _time_array, expand_pattern

    out = expand_pattern(template, index, count)
    mesh.field_data[TIME_KEY] = _time_array(time)
    write(out, mesh, file_format=file_format or None, **kwargs)
    return out


# --------------------------------------------------------------------------- #
# read options and status (roadmap 6.5)                                       #
# --------------------------------------------------------------------------- #
#: The values of the reader's Ghosts property, ``read``'s ``ghosts``.
GHOST_POLICIES = ("keep", "drop")
#: The formats the tessellation options reach (``step.read``, ``iges.read``).
CAD_FORMATS = ("step", "iges")
_METADATA: dict = {}


def metadata(path, file_format=None):
    """``read_metadata(path)``, cached per file identity as the plans are."""
    from ._helpers import read_metadata

    key = _plan_key([str(path)], file_format)
    if key not in _METADATA:
        _METADATA[key] = read_metadata(str(path), file_format)
    return _METADATA[key]


def read_kwargs(
    point_names=(),
    cell_names=(),
    point_enabled=None,
    cell_enabled=None,
    points_only=False,
    piece=-1,
    ghosts="keep",
):
    """``meshioplusplus.read`` keywords for the reader's properties.

    ``point_enabled``/``cell_enabled`` are the ticked names (``None``: all);
    ``arrays=`` is passed only when one is unticked, so a format that cannot
    filter its arrays is not asked to. ``piece`` below 0 reads every piece.
    """
    kwargs = {}
    if points_only:
        kwargs["points_only"] = True
    point_on = list(point_names) if point_enabled is None else list(point_enabled)
    cell_on = list(cell_names) if cell_enabled is None else list(cell_enabled)
    every = set(point_names) | set(cell_names)
    wanted = {n for n in point_on if n in point_names} | {
        n for n in cell_on if n in cell_names
    }
    if wanted != every:
        kwargs["arrays"] = sorted(wanted)
    if piece is not None and piece >= 0:
        kwargs["piece"] = int(piece)
    if ghosts not in GHOST_POLICIES:
        raise ValueError(f"ghosts must be one of {GHOST_POLICIES}, not {ghosts!r}")
    if ghosts != "keep":
        kwargs["ghosts"] = ghosts
    return kwargs


def cad_format(path, file_format=None):
    """``"step"`` or ``"iges"`` when ``path`` is one (by format or extension),
    else ``None``."""
    import pathlib

    from ._helpers import formats

    if file_format:
        return file_format if file_format in CAD_FORMATS else None
    fmts = formats()["extensions"].get(pathlib.Path(str(path)).suffix.lower(), [])
    return next((f for f in fmts if f in CAD_FORMATS), None)


def read_cad(path, fmt, deflection=0.0, angular_deflection=0.5):
    """A STEP or IGES file triangulated at the given deflections."""
    import importlib

    module = importlib.import_module(f"meshioplusplus.{fmt}")
    return module.read(
        str(path), deflection=deflection, angular_deflection=angular_deflection
    )


def addon_status():
    """Which add-ons this meshio++ has: ``{"occt", "netgen", "mmg"}`` -> bool.
    A failed STEP read in ParaView is explained by the first."""
    from ._fit_cad import has_occt
    from ._mmg import has_mmg
    from ._tetrahedralize import has_netgen

    return {"occt": has_occt(), "netgen": has_netgen(), "mmg": has_mmg()}


# --------------------------------------------------------------------------- #
# regions as blocks (roadmap 6.6.2)                                           #
# --------------------------------------------------------------------------- #
#: The assembly's groups, by region kind, in this order.
REGION_GROUPS = {
    "cell": "cell_regions",
    "side": "side_regions",
    "point": "point_regions",
}


def _side_piece(mesh, region):
    """A side region as a mesh of its facets (``facet_nodes``), one block per
    facet type in first-seen order, over the whole mesh's points."""
    from ._facets import facet_nodes
    from ._mesh import Mesh

    rows: dict = {}
    for cell, facet in np.asarray(region.entries, dtype=np.int64).reshape(-1, 2):
        got = facet_nodes(mesh, int(cell), int(facet))
        if got is not None:
            rows.setdefault(got[0], []).append(np.asarray(got[1], dtype=np.int64))
    blocks = [(t, np.vstack(r)) for t, r in rows.items()]
    return Mesh(mesh.points, blocks)


def region_pieces(mesh):
    """``[(group, name, Mesh)]``: one piece per named region, the blocks of the
    reader's *RegionsAsBlocks* output.

    Cell regions are :func:`meshioplusplus.split` ``(by="regions")``'s pieces
    (which may overlap), side regions their facets (``facet_nodes``), point
    regions their points as ``vertex`` cells. A mesh with no cell region is
    also one piece, ``mesh``, so its cells are not lost.
    """
    from ._mesh import Mesh
    from ._split import split

    regions = list(getattr(mesh, "regions", []))
    pieces = []
    cell_regions = [r for r in regions if r.kind == "cell"]
    if cell_regions:
        # Split a view holding the cell regions alone: the others are pieces
        # of their own, and split would warn that each piece lost them.
        view = Mesh(
            mesh.points,
            mesh.cells,
            point_data=mesh.point_data,
            cell_data=mesh.cell_data,
            regions=cell_regions,
        )
        for name, piece in split(view, by="regions").items():
            pieces.append((REGION_GROUPS["cell"], str(name), piece))
    else:
        pieces.append(("mesh", "mesh", mesh))
    for r in regions:
        if r.kind == "side":
            pieces.append((REGION_GROUPS["side"], r.name, _side_piece(mesh, r)))
        elif r.kind == "point":
            ids = np.asarray(r.entries, dtype=np.int64).reshape(-1, 1)
            pieces.append(
                (REGION_GROUPS["point"], r.name, Mesh(mesh.points, [("vertex", ids)]))
            )
    return pieces


def fill_collection(collection, mesh, op=READER_OP):
    """Fill a ``vtkPartitionedDataSetCollection`` with ``mesh``'s regions,
    one partitioned dataset each, named and grouped in a ``vtkDataAssembly``
    (``cell_regions``, ``side_regions``, ``point_regions``), as ParaView's
    own Exodus reader presents its blocks and sets."""
    from vtkmodules.vtkCommonDataModel import (
        vtkCompositeDataSet,
        vtkDataAssembly,
        vtkPartitionedDataSet,
    )

    assembly = vtkDataAssembly()
    assembly.SetRootNodeName("meshioplusplus")
    groups = {}
    collection.SetNumberOfPartitionedDataSets(0)
    for index, (group, name, piece) in enumerate(region_pieces(mesh)):
        partitioned = vtkPartitionedDataSet()
        partitioned.SetPartition(0, grid_from_mesh(piece, op))
        collection.SetPartitionedDataSet(index, partitioned)
        collection.GetMetaData(index).Set(vtkCompositeDataSet.NAME(), name)
        if group not in groups:
            groups[group] = assembly.AddNode(group, 0)
        node = assembly.AddNode(vtkDataAssembly.MakeValidNodeName(name), groups[group])
        assembly.AddDataSetIndex(node, index)
    collection.SetDataAssembly(assembly)
    return collection
