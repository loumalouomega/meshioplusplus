"""Blender bridge: a mesh in and out of Blender's own data model.

Blender ships Python and reads almost no FEA formats — STL, OBJ, PLY and
essentially nothing else. meshio++ reads 43. This module is the bridge, and
``src/blender/`` packages it as a Blender 4.2+ extension whose only job is to
put ``File > Import`` in front of it.

Nothing here is part of the C++ core. It is pure Python over the numpy the
readers already return, and ``bpy`` is imported lazily and named.

Architecture
------------
Split exactly the way :mod:`meshioplusplus._interop` and
:mod:`meshioplusplus.mcp` are, and for the same reason. The whole of the
behaviour lives in a **pure payload layer** — :func:`_to_blender_payload` and
its inverse :func:`_mesh_from_blender_arrays` — which imports no ``bpy`` at
all, does not mutate its input, and speaks in flat numpy arrays. That is what
makes the subtle parts (the loop/polygon layout, the block-major indexing, the
attribute-domain mapping) testable in the default CI matrix with Blender
nowhere in sight. :func:`to_blender` and :func:`from_blender` are thin wrappers
that import ``bpy`` *inside* the function.

n-gons are kept, not triangulated
---------------------------------
This is the one place the design deliberately parts company with
:func:`~meshioplusplus._interop.to_trimesh`. ``_interop._triangulate`` ends in
``convert_cells("simplexify")`` because a ``trimesh.Trimesh`` holds triangles
and nothing else. Blender's mesh is vertices + **loops** + **polygons** — the
same shape as ``vtkPolyData`` — and holds triangles, quads and n-gons natively,
so simplexifying would destroy exactly the quad topology someone importing a
hex-dominant mesh came for. The reduction used here is instead
:func:`~meshioplusplus._viewer_browser._renderable_surface`, which extracts a
volume mesh's boundary, gathers each ``cell_data`` array through the owning
cell, and linearizes — and stops there.

That function is **called, not copied**. Its docstring binds it to the WASM
``convertSurface`` binding ("the two must agree"); a second caller costs
nothing, while a second implementation would put that contract out of reach of
the file documenting it. It also means a mesh looks the same in Blender, in the
browser viewer and in the hosted demo.

Why there is no ``[blender]`` extra
-----------------------------------
``bpy`` on PyPI is a several-hundred-megabyte wheel pinned to one CPython
minor, and inside Blender — the only place this code actually runs — ``bpy`` is
a builtin that pip must never touch; installing it there would install a second
Blender. So the gate is :func:`meshioplusplus._gpu._require_framework` (the
torch / JAX / CuPy shape), not ``_interop._require``, and the error names
Blender rather than a ``pip install meshioplusplus[...]`` that would be wrong
twice over. See ``doc/blender.md``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import NamedTuple

import numpy as np

from ._common import warn
from ._interop import (
    REGION_META_KEY,
    REGION_PREFIX,
    _derived_note,
    _emit,
    _flatten_cell_data,
    _importable,
    _region_payload,
    _regions_from_masks,
)

#: Cell types that become Blender **polygons**. ``polygon`` covers the ragged
#: storage shape and ``polygon<n>`` the uniform one; both are n-gons to Blender.
#: Checked by type *name*, never ``IsRagged()`` — a uniform n-gon block stores
#: rectangularly and is still an n-gon (the ``cgns.cpp:476`` lesson).
_POLYGON_TYPES = ("triangle", "quad", "polygon")

#: Cell types that become Blender **edges**. Higher-order lines are linearized
#: to ``line`` before we ever see them.
_EDGE_TYPES = ("line",)

#: Cell types carrying no Blender geometry of their own. Every point is
#: imported regardless, so a ``vertex`` block is a no-op rather than a loss.
_LOOSE_TYPES = ("vertex",)

#: Attribute names Blender either owns or refuses. A ``.``-prefixed name is
#: Blender's internal namespace; the rest are built-ins whose type and domain
#: are fixed, so writing one would either fail or corrupt the mesh.
#:
#: Deliberately **conservative**: it cannot be exhaustive across Blender
#: versions, so :func:`to_blender` additionally catches a rejected
#: ``attributes.new()`` and turns it into a note. Neither mechanism alone is
#: enough — this one keeps the common cases readable, that one keeps an
#: unforeseen name from aborting an import.
_RESERVED_ATTRIBUTE_NAMES = frozenset(
    (
        "position",
        "material_index",
        "sharp_face",
        "sharp_edge",
        "crease",
        "crease_edge",
        "crease_vert",
        "id",
        "radius",
        "resolution",
        "velocity",
    )
)

#: Blender's own limit on a name, in bytes.
_MAX_ATTRIBUTE_NAME = 64

#: Prefix for a ``field_data`` entry parked in the mesh's custom properties.
FIELD_PREFIX = "meshioplusplus:field:"


class BlenderAttr(NamedTuple):
    """One attribute layer, ready for ``layer.data.foreach_set(field, values)``."""

    #: Blender's ``data_type`` string.
    data_type: str
    #: The property ``foreach_set`` writes: ``"value"`` for scalars,
    #: ``"vector"`` for ``FLOAT_VECTOR``/``FLOAT2``.
    field: str
    #: Flat, already in the dtype Blender stores.
    values: np.ndarray
    #: Components per element (1 for a scalar).
    components: int = 1


@dataclass(frozen=True)
class BlenderPayload:
    """Everything a ``bpy.types.Mesh`` needs, as flat numpy in Blender's dtypes.

    ``loop_starts[i] : loop_starts[i] + loop_totals[i]`` slices polygon ``i``
    out of ``loop_vertices``. That is the same layout as
    :class:`~meshioplusplus._interop.VtkPayload`'s ``offsets``/``connectivity``
    under different names, which is why
    :func:`~meshioplusplus._interop._vtk_block_conn` builds it here.
    """

    #: ``(V, 3)`` float32, C-contiguous — ``MeshVertex.co``.
    vertices: np.ndarray
    #: ``(L,)`` int32 — ``MeshLoop.vertex_index``.
    loop_vertices: np.ndarray
    #: ``(P,)`` int32 — ``MeshPolygon.loop_start``.
    loop_starts: np.ndarray
    #: ``(P,)`` int32 — ``MeshPolygon.loop_total``.
    loop_totals: np.ndarray
    #: ``(E, 2)`` int32 — loose edges, from ``line`` blocks only.
    edge_vertices: np.ndarray
    #: ``name -> BlenderAttr`` on the ``POINT`` domain.
    point_attributes: dict = field(default_factory=dict)
    #: ``name -> BlenderAttr`` on the ``FACE`` domain.
    face_attributes: dict = field(default_factory=dict)
    #: JSON-safe values for the mesh datablock's custom properties.
    custom_properties: dict = field(default_factory=dict)
    #: Ascending global block-major cell index of each emitted polygon. Regions
    #: resolve through this; see ``_interop._region_payload``.
    kept_global: np.ndarray = field(default_factory=lambda: np.empty(0, np.int64))
    #: Operations composed to get here, e.g. ``["extract_surface", "linearize"]``.
    ops: list = field(default_factory=list)
    #: ``(block index, cell type)`` of every block Blender cannot hold.
    dropped: list = field(default_factory=list)
    #: Lossy or surprising steps; surfaced as warnings by the wrappers.
    notes: list = field(default_factory=list)
    #: Arrays with no meshio++ counterpart. See ``_interop._derived_note``.
    derived: list = field(default_factory=list)

    @property
    def num_polygons(self) -> int:
        return len(self.loop_totals)


# --------------------------------------------------------------------------- #
# attribute naming and typing                                                 #
# --------------------------------------------------------------------------- #
def _safe_attribute_name(name, used, notes):
    """A name Blender will accept, unique within ``used``."""
    out = str(name)
    if out.startswith("."):
        # Blender's internal namespace. `.select_vert` and friends live there.
        out = "attr" + out
        notes.append(f"attribute '{name}' was renamed to '{out}' (leading dot)")
    if out in _RESERVED_ATTRIBUTE_NAMES:
        out = f"attr_{out}"
        notes.append(f"attribute '{name}' was renamed to '{out}' (reserved by Blender)")
    encoded = out.encode("utf-8")
    if len(encoded) > _MAX_ATTRIBUTE_NAME:
        out = encoded[:_MAX_ATTRIBUTE_NAME].decode("utf-8", "ignore")
        notes.append(f"attribute '{name}' was truncated to '{out}'")
    if out in used:
        stem = out
        i = 1
        while f"{stem}_{i}" in used:
            i += 1
        out = f"{stem}_{i}"
        notes.append(f"attribute '{name}' was renamed to '{out}' (name already used)")
    used.add(out)
    return out


def _blender_attribute(name, raw, n_expected, used, notes, derived):
    """One data array as one or more ``(name, BlenderAttr)`` pairs.

    Returns ``[]`` — with a note — when the array cannot be represented at all,
    which is a drop rather than an error: the ``_to_vtk_payload`` policy, since
    one unusable array must not cost the user the whole import.

    Multi-component arrays map onto Blender's own vector types where one
    exists, and otherwise expand into ``v_0``..``v_{k-1}`` scalars. That suffix
    spelling is :func:`~meshioplusplus._interop._frame_columns`' rule, which
    ``doc/interop.md`` already calls the one component-flattening convention in
    the repo; only the *name* is shared, since that function operates on a
    ``TablePayload`` whose structural columns have no meaning here.
    """
    a = np.asarray(raw)
    if a.dtype == object or np.issubdtype(a.dtype, np.complexfloating):
        notes.append(f"attribute '{name}' has dtype {a.dtype} and was dropped")
        return []
    if a.dtype.kind in "SU":
        notes.append(f"attribute '{name}' is a string array and was dropped")
        return []
    if a.ndim == 0 or a.shape[0] != n_expected:
        got = "scalar" if a.ndim == 0 else str(a.shape[0])
        notes.append(
            f"attribute '{name}' has {got} row(s) for {n_expected} element(s); dropped"
        )
        return []

    flat = a.reshape(a.shape[0], -1)
    k = int(flat.shape[1])
    is_float = flat.dtype.kind == "f"
    is_bool = flat.dtype == bool

    if k == 1:
        column = flat[:, 0]
        if is_bool:
            return [(_safe_attribute_name(name, used, notes), _attr_bool(column))]
        if is_float:
            return [(_safe_attribute_name(name, used, notes), _attr_float(column))]
        return [
            (_safe_attribute_name(name, used, notes), _attr_int(name, column, notes))
        ]

    # Blender's vector types are float-only, so an integer or boolean vector
    # expands rather than being silently widened to float.
    if is_float and k in (2, 3):
        data_type = "FLOAT2" if k == 2 else "FLOAT_VECTOR"
        values = np.ascontiguousarray(flat, dtype=np.float32).reshape(-1)
        if flat.dtype != np.float32:
            _derived_note(
                f"attribute '{name}'", "Blender stores attributes as float32", derived
            )
        return [
            (
                _safe_attribute_name(name, used, notes),
                BlenderAttr(data_type, "vector", values, k),
            )
        ]

    notes.append(
        f"attribute '{name}' has {k} components and was expanded into "
        f"{k} scalar attributes"
    )
    out = []
    for i in range(k):
        column = np.ascontiguousarray(flat[:, i])
        safe = _safe_attribute_name(f"{name}_{i}", used, notes)
        if is_bool:
            out.append((safe, _attr_bool(column)))
        elif is_float:
            out.append((safe, _attr_float(column)))
        else:
            out.append((safe, _attr_int(name, column, notes)))
    return out


def _attr_float(column):
    return BlenderAttr(
        "FLOAT", "value", np.ascontiguousarray(column, dtype=np.float32), 1
    )


def _attr_bool(column):
    return BlenderAttr("BOOLEAN", "value", np.ascontiguousarray(column, dtype=bool), 1)


def _attr_int(name, column, notes):
    """An integer scalar as Blender's ``INT`` (a C ``int``, i.e. 32-bit)."""
    a = np.ascontiguousarray(column)
    if a.dtype != np.int32:
        info = np.iinfo(np.int32)
        if a.size and (int(a.min()) < info.min or int(a.max()) > info.max):
            notes.append(
                f"attribute '{name}' does not fit Blender's 32-bit integer "
                "attribute and was stored as float"
            )
            return _attr_float(a)
        a = a.astype(np.int32)
    return BlenderAttr("INT", "value", a, 1)


# --------------------------------------------------------------------------- #
# the pure payload                                                            #
# --------------------------------------------------------------------------- #
def _blender_points(mesh, notes, derived):
    """``(V, 3)`` float32 vertices — ``_vtk_points``' shape logic, Blender's dtype."""
    pts = np.asarray(mesh.points)
    if pts.ndim == 1:
        pts = pts.reshape(-1, 1)
    n, dim = pts.shape
    if dim != 3:
        notes.append(f"points are {dim}D and Blender needs 3 columns (Z padded)")
    if dim == 3 and pts.dtype == np.float32 and pts.flags["C_CONTIGUOUS"]:
        return pts
    # Always a copy: Blender's `MeshVertex.co` is float[3], and `foreach_set`
    # only memcpys on an exact dtype match -- otherwise it converts element by
    # element, which is far slower than converting once here.
    _derived_note("vertices", "Blender stores coordinates as float32", derived)
    out = np.zeros((n, 3), dtype=np.float32)
    keep = min(dim, 3)
    out[:, :keep] = pts[:, :keep]
    return out


def _reduce_for_blender(mesh, ops):
    """The surface Blender can draw, and an honest record of how we got it.

    Calls :func:`~meshioplusplus._viewer_browser._renderable_surface`; ``ops``
    is rebuilt from the *same two predicates* that function branches on, so it
    stays accurate without changing a signature the WASM twin depends on.
    """
    from ._convert_cells import _LINEAR_BASE
    from ._viewer_browser import _renderable_surface

    had_volume = any(blk.dim == 3 for blk in mesh.cells)
    if had_volume:
        ops.append("extract_surface")
        # extract_surface preserves element order, so a higher-order volume
        # yields higher-order faces: the linearize predicate is decided by the
        # 3-D blocks alone. A `tetra` + `triangle6` mesh linearizes nothing,
        # because the supplied `triangle6` block never reaches the output.
        linearized = any(blk.type in _LINEAR_BASE for blk in mesh.cells if blk.dim == 3)
    else:
        linearized = any(blk.type in _LINEAR_BASE for blk in mesh.cells)
    if linearized:
        ops.append("linearize")
    # "surface" is a literal: the `kind` parameter exists for the browser
    # viewer's own error message, and Blender never offers the other values, so
    # that branch is unreachable from here.
    return _renderable_surface(mesh, "surface")


def _to_blender_payload(
    mesh,
    *,
    point_data: bool = True,
    cell_data: bool = True,
    regions: bool = True,
    field_data: bool = True,
) -> BlenderPayload:
    """Map a mesh onto the flat arrays a ``bpy.types.Mesh`` is built from.

    Pure: imports no ``bpy``, touches no Blender object, and does not modify
    ``mesh``. It does call meshio++'s own operations, each of which returns a
    new mesh.

    Blocks Blender cannot hold — ``polyhedron``, ``custom``, anything left with
    a dimension of 3 — are **dropped with a note naming them**, never an
    exception. A mixed mesh is the normal case.
    """
    from ._interop import _vtk_block_conn
    from ._regions import block_bases

    notes: list = []
    derived: list = []
    dropped: list = []
    ops: list = []

    surface = _reduce_for_blender(mesh, ops)
    vertices = _blender_points(surface, notes, derived)
    bases = block_bases(surface.cells)

    poly_blocks: list = []
    conn_parts: list = []
    count_parts: list = []
    global_parts: list = []
    edge_parts: list = []
    loose: list = []
    unsupported: list = []

    for b, cb in enumerate(surface.cells):
        base = (
            cb.type.rstrip("0123456789") if cb.type.startswith("polygon") else cb.type
        )
        if base in _POLYGON_TYPES:
            flat, counts, _shareable, _reason = _vtk_block_conn(cb)
            poly_blocks.append(b)
            conn_parts.append(flat)
            count_parts.append(counts)
            global_parts.append(np.arange(bases[b], bases[b] + len(cb), dtype=np.int64))
        elif base in _EDGE_TYPES:
            data = np.asarray(cb.data)
            if data.ndim == 2 and data.shape[1] == 2:
                edge_parts.append(data)
            else:
                unsupported.append(cb.type)
                dropped.append((b, cb.type))
        elif base in _LOOSE_TYPES:
            loose.append(cb.type)
        else:
            unsupported.append(cb.type)
            dropped.append((b, cb.type))

    if unsupported:
        notes.append(
            "cell blocks Blender cannot hold were dropped: "
            + ", ".join(sorted(set(unsupported)))
        )
    if loose:
        notes.append(
            f"'{sorted(set(loose))[0]}' blocks contribute no geometry; their "
            "points are imported like every other point"
        )

    if conn_parts:
        loop_vertices = np.ascontiguousarray(np.concatenate(conn_parts), dtype=np.int32)
        loop_totals = np.ascontiguousarray(np.concatenate(count_parts), dtype=np.int32)
        kept_global = np.concatenate(global_parts)
    else:
        loop_vertices = np.empty(0, dtype=np.int32)
        loop_totals = np.empty(0, dtype=np.int32)
        kept_global = np.empty(0, dtype=np.int64)
    offsets = np.cumsum(loop_totals, dtype=np.int64)
    loop_starts = np.ascontiguousarray(
        np.concatenate([[0], offsets[:-1]]) if len(loop_totals) else offsets,
        dtype=np.int32,
    )
    _derived_note("loop_starts", "Blender indexes polygons by loop offset", derived)

    if edge_parts:
        edge_vertices = np.ascontiguousarray(
            np.concatenate(edge_parts, axis=0), dtype=np.int32
        )
    else:
        edge_vertices = np.empty((0, 2), dtype=np.int32)

    used_names: set = set()
    point_attributes: dict = {}
    face_attributes: dict = {}
    custom_properties: dict = {}

    if point_data:
        for name in sorted(surface.point_data):
            for safe, attr in _blender_attribute(
                name,
                surface.point_data[name],
                len(vertices),
                used_names,
                notes,
                derived,
            ):
                point_attributes[safe] = attr

    if cell_data:
        for name in sorted(surface.cell_data):
            flat = _flatten_cell_data(surface, name, poly_blocks, notes)
            if flat is None:
                continue
            for safe, attr in _blender_attribute(
                name, flat, len(loop_totals), used_names, notes, derived
            ):
                face_attributes[safe] = attr

    if regions:
        masks, meta = _region_payload(surface, len(vertices), kept_global, notes)
        for (kind, name), mask in masks.items():
            target = point_attributes if kind == "point" else face_attributes
            safe = _safe_attribute_name(REGION_PREFIX + name, used_names, notes)
            target[safe] = BlenderAttr(
                "BOOLEAN", "value", np.ascontiguousarray(mask != 0), 1
            )
        if meta:
            custom_properties[REGION_META_KEY] = json.dumps(meta)

    if field_data:
        for name in sorted(surface.field_data):
            value = surface.field_data[name]
            try:
                encoded = np.asarray(value).tolist()
                json.dumps(encoded)
            except (TypeError, ValueError):
                notes.append(f"field_data '{name}' is not JSON-safe and was dropped")
                continue
            custom_properties[FIELD_PREFIX + name] = encoded

    return BlenderPayload(
        vertices=vertices,
        loop_vertices=loop_vertices,
        loop_starts=loop_starts,
        loop_totals=loop_totals,
        edge_vertices=edge_vertices,
        point_attributes=point_attributes,
        face_attributes=face_attributes,
        custom_properties=custom_properties,
        kept_global=kept_global,
        ops=ops,
        dropped=dropped,
        notes=notes,
        derived=derived,
    )


# --------------------------------------------------------------------------- #
# the pure inverse                                                            #
# --------------------------------------------------------------------------- #
def _polygon_edge_keys(loop_vertices, loop_starts, loop_totals):
    """The sorted vertex pairs already implied by the polygons."""
    keys = set()
    for start, total in zip(loop_starts, loop_totals):
        ring = loop_vertices[start : start + total]
        for i in range(total):
            a = int(ring[i])
            b = int(ring[(i + 1) % total])
            keys.add((a, b) if a < b else (b, a))
    return keys


def _pad_cell_data(values, n_extra, name, notes):
    """Extend a per-polygon array to cover a trailing ``line`` block.

    A Blender edge carries no FACE attribute, so those rows genuinely have no
    value. Floats get NaN, which *is* the "no value" a float can hold — the
    ``merge`` operation's ``Fill`` policy. An integer has no such value, so the
    array is dropped rather than given an invented 0: refusing beats guessing,
    the same call ``_interop._flatten_cell_data`` makes.
    """
    if values.dtype.kind != "f":
        notes.append(
            f"cell_data '{name}' was dropped: the mesh has loose edges, which "
            f"carry no value, and dtype {values.dtype} has no missing-value "
            "representation"
        )
        return None
    pad = np.full((n_extra,) + values.shape[1:], np.nan, dtype=values.dtype)
    return np.concatenate([values, pad], axis=0)


def _mesh_from_blender_arrays(
    vertices,
    loop_vertices,
    loop_starts,
    loop_totals,
    edge_vertices=None,
    point_attributes=None,
    face_attributes=None,
    custom_properties=None,
    notes=None,
):
    """Rebuild a :class:`~meshioplusplus._mesh.Mesh` from Blender's own arrays.

    Takes exactly what ``foreach_get`` produces and imports no ``bpy``, which
    is what puts the whole round-trip inside the default CI matrix.
    ``point_attributes``/``face_attributes`` are plain ``{name: (n,) or (n, k)
    array}`` dicts.

    The block reconstruction is
    :func:`~meshioplusplus._vtk_common.vtk_cells_from_data` — the same function
    :func:`~meshioplusplus._interop.from_pyvista` uses — because Blender's
    loop layout and VTK 9's cell array are the same thing. It splits on *runs*
    of equal cell type, so the polygons are first stably sorted by side count:
    without that, an alternating triangle/quad mesh would come back as one
    block per cell.
    """
    from ._mesh import CellBlock, Mesh
    from ._vtk_common import vtk_cells_from_data

    notes = [] if notes is None else notes
    point_attributes = dict(point_attributes or {})
    face_attributes = dict(face_attributes or {})
    custom_properties = dict(custom_properties or {})

    points = np.ascontiguousarray(vertices, dtype=np.float64).reshape(-1, 3)
    loop_vertices = np.asarray(loop_vertices).reshape(-1)
    loop_starts = np.asarray(loop_starts).reshape(-1)
    loop_totals = np.asarray(loop_totals).reshape(-1)

    order = np.argsort(loop_totals, kind="stable")
    totals = loop_totals[order]
    starts = loop_starts[order]
    if len(totals):
        gather = np.concatenate(
            [np.arange(s, s + t, dtype=np.int64) for s, t in zip(starts, totals)]
        )
    else:
        gather = np.empty(0, dtype=np.int64)
    connectivity = np.ascontiguousarray(loop_vertices[gather], dtype=np.int64)
    offsets = np.cumsum(totals, dtype=np.int64)
    types = np.where(totals == 3, 5, np.where(totals == 4, 9, 7)).astype(np.int64)

    face_arrays = {k: np.asarray(v)[order] for k, v in face_attributes.items()}
    region_face = {
        k: face_arrays.pop(k) for k in list(face_arrays) if k.startswith(REGION_PREFIX)
    }
    if len(totals):
        cells, cell_data = vtk_cells_from_data(
            connectivity, offsets, types, face_arrays
        )
    else:
        # A Blender mesh with no faces is ordinary -- a point cloud, or a curve
        # imported as loose edges. `vtk_cells_from_data` indexes `types[0]`
        # unconditionally, so it must not be called with an empty run.
        cells, cell_data = [], {}

    # Loose edges -- those no polygon already implies -- become a `line` block,
    # emitted last so a cell region's global block-major indices, which are
    # positions among the polygons, stay correct.
    if edge_vertices is not None and len(edge_vertices):
        pairs = np.asarray(edge_vertices, dtype=np.int64).reshape(-1, 2)
        implied = _polygon_edge_keys(loop_vertices, loop_starts, loop_totals)
        keep = [
            i
            for i, (a, b) in enumerate(pairs.tolist())
            if ((a, b) if a < b else (b, a)) not in implied
        ]
        if keep:
            n_extra = len(keep)
            cells = list(cells) + [CellBlock("line", pairs[keep])]
            for name in list(cell_data):
                joined = np.concatenate(cell_data[name], axis=0)
                padded = _pad_cell_data(joined, n_extra, name, notes)
                if padded is None:
                    del cell_data[name]
                    continue
                sizes = [len(b) for b in cell_data[name]] + [n_extra]
                bounds = np.cumsum([0] + sizes)
                cell_data[name] = [
                    padded[bounds[i] : bounds[i + 1]] for i in range(len(sizes))
                ]

    point_arrays = dict(point_attributes)
    region_point = {
        k: point_arrays.pop(k)
        for k in list(point_arrays)
        if k.startswith(REGION_PREFIX)
    }

    meta_by_name: dict = {}
    raw_meta = custom_properties.get(REGION_META_KEY)
    if raw_meta:
        try:
            for entry in json.loads(raw_meta):
                meta_by_name[(entry["kind"], entry["name"])] = entry
        except (TypeError, ValueError, KeyError):
            notes.append(f"'{REGION_META_KEY}' is malformed; region dim/tag are lost")

    regions = _regions_from_masks(region_point, "point", meta_by_name, notes)
    regions += _regions_from_masks(region_face, "cell", meta_by_name, notes)

    field_data = {
        key[len(FIELD_PREFIX) :]: np.asarray(value)
        for key, value in custom_properties.items()
        if key.startswith(FIELD_PREFIX)
    }

    return Mesh(
        points,
        cells,
        point_data={k: np.asarray(v) for k, v in point_arrays.items()},
        cell_data=cell_data,
        field_data=field_data,
        regions=regions,
    )


def _attribute_arrays(attributes, count):
    """``{name: BlenderAttr}`` back to the ``{name: (n,) or (n, k)}`` shape.

    The pure counterpart of what :func:`from_blender` gets out of
    ``foreach_get``, so a round-trip test needs no Blender.
    """
    out = {}
    for name, attr in attributes.items():
        values = np.asarray(attr.values)
        out[name] = values if attr.components == 1 else values.reshape(count, -1)
    return out


# --------------------------------------------------------------------------- #
# the gated wrappers                                                          #
# --------------------------------------------------------------------------- #
def has_blender() -> bool:
    """Whether ``bpy`` is importable — Blender's own Python, or the PyPI wheel."""
    return _importable("bpy")


def _require_bpy(op):
    """Import ``bpy`` or raise naming Blender.

    Deliberately not ``_interop._require``: its message is
    ``pip install meshioplusplus[extra]``, which would be wrong twice over
    here. See this module's docstring.
    """
    from ._gpu import _require_framework

    return _require_framework(
        op,
        "bpy",
        "run this inside Blender, or `pip install bpy` for a headless build",
        doc="doc/blender.md",
        qualifier="(the PyPI wheel is pinned to one CPython version)",
    )


def _fill_mesh(me, payload):
    """Populate a ``bpy.types.Mesh`` from ``payload``, fast path first.

    Blender 3.6 moved polygon storage to an offset-indices array, and
    ``MeshPolygon.loop_total`` is derived from it rather than stored — so
    whether it can be written is a **capability**, not a version number. Probe
    it once and fall back to ``from_pydata``, which supports n-gons on every
    supported release at the cost of building Python lists.
    """
    me.vertices.add(len(payload.vertices))
    me.loops.add(len(payload.loop_vertices))
    me.polygons.add(payload.num_polygons)
    me.vertices.foreach_set("co", payload.vertices.reshape(-1))
    me.loops.foreach_set("vertex_index", payload.loop_vertices)
    me.polygons.foreach_set("loop_start", payload.loop_starts)
    try:
        me.polygons.foreach_set("loop_total", payload.loop_totals)
    except (AttributeError, RuntimeError, TypeError):
        return False
    if len(payload.edge_vertices):
        me.edges.add(len(payload.edge_vertices))
        me.edges.foreach_set("vertices", payload.edge_vertices.reshape(-1))
    me.update(calc_edges=True)
    return True


def _fill_mesh_from_pydata(me, payload):
    """The portable fallback: build the polygons as Python lists."""
    starts = payload.loop_starts.tolist()
    totals = payload.loop_totals.tolist()
    loops = payload.loop_vertices.tolist()
    faces = [loops[s : s + t] for s, t in zip(starts, totals)]
    edges = payload.edge_vertices.tolist()
    me.from_pydata(payload.vertices.tolist(), edges, faces)
    me.update(calc_edges=True)


def to_blender(
    mesh,
    name: str = "mesh",
    *,
    point_data: bool = True,
    cell_data: bool = True,
    regions: bool = True,
    field_data: bool = True,
    validate: bool = True,
):
    """Convert a mesh to a ``bpy.types.Mesh`` datablock.

    Returns the **datablock**, not an ``Object``, and links nothing into a
    collection — the same separation as :func:`~meshioplusplus.to_pyvista`
    returning a grid rather than a plotter. The add-on wraps it.

    A mesh with volume cells arrives as its **boundary surface**: Blender has
    no tetrahedron, hexahedron, wedge or pyramid, so no other mapping exists.
    Each ``cell_data`` array is carried through the owning cell, so a solid can
    still be coloured by its per-cell material. Quads and n-gons are **kept**,
    not triangulated. See ``doc/blender.md``.

    There is deliberately no ``zero_copy_only``: Blender stores coordinates as
    float32 and derives loops and polygons from a layout no mesh holds, so
    every array is a copy by construction and the flag could only ever raise
    (the :func:`~meshioplusplus.to_polars` precedent).
    """
    bpy = _require_bpy("to_blender")

    payload = _to_blender_payload(
        mesh,
        point_data=point_data,
        cell_data=cell_data,
        regions=regions,
        field_data=field_data,
    )
    _emit("to_blender", payload.notes)

    me = bpy.data.meshes.new(name)
    if not _fill_mesh(me, payload):
        me.clear_geometry()
        _fill_mesh_from_pydata(me, payload)

    for domain, attributes in (
        ("POINT", payload.point_attributes),
        ("FACE", payload.face_attributes),
    ):
        for attr_name, attr in attributes.items():
            try:
                layer = me.attributes.new(attr_name, attr.data_type, domain)
                layer.data.foreach_set(attr.field, attr.values)
            except (RuntimeError, TypeError, ValueError) as exc:
                # A name or type Blender refuses must cost one array, never the
                # whole import -- the reserved-name list cannot be exhaustive
                # across releases, so this is the second half of that guard.
                warn(f"to_blender: attribute '{attr_name}' was rejected: {exc}")

    for key, value in payload.custom_properties.items():
        me[key] = value

    if validate:
        me.validate(verbose=False)
    return me


def from_blender(source, *, apply_modifiers: bool = False, attributes: bool = True):
    """Convert a ``bpy.types.Mesh`` or ``bpy.types.Object`` back to a mesh.

    Lossy by construction: the result holds ``triangle``/``quad``/``polygon``
    and ``line`` blocks only, because that is everything Blender has. Attributes
    on the ``CORNER`` and ``EDGE`` domains are dropped with a note — meshio++
    has point and cell data locations and no others.
    """
    bpy = _require_bpy("from_blender")

    notes: list = []
    owner = None
    if isinstance(source, bpy.types.Object):
        if apply_modifiers:
            depsgraph = bpy.context.evaluated_depsgraph_get()
            owner = source.evaluated_get(depsgraph)
            me = owner.to_mesh()
        else:
            me = source.data
    else:
        me = source

    try:
        n_v, n_l = len(me.vertices), len(me.loops)
        n_p, n_e = len(me.polygons), len(me.edges)

        vertices = np.empty(n_v * 3, dtype=np.float32)
        me.vertices.foreach_get("co", vertices)
        loop_vertices = np.empty(n_l, dtype=np.int32)
        me.loops.foreach_get("vertex_index", loop_vertices)
        loop_starts = np.empty(n_p, dtype=np.int32)
        me.polygons.foreach_get("loop_start", loop_starts)
        loop_totals = np.empty(n_p, dtype=np.int32)
        me.polygons.foreach_get("loop_total", loop_totals)
        edge_vertices = np.empty(n_e * 2, dtype=np.int32)
        me.edges.foreach_get("vertices", edge_vertices)

        point_attributes: dict = {}
        face_attributes: dict = {}
        if attributes:
            for layer in me.attributes:
                got = _read_attribute(layer, n_v, n_p, notes)
                if got is None:
                    continue
                target, values = got
                target_dict = point_attributes if target == "POINT" else face_attributes
                target_dict[layer.name] = values

        custom_properties = {
            key: me[key]
            for key in me.keys()
            if isinstance(key, str)
            and (key == REGION_META_KEY or key.startswith(FIELD_PREFIX))
        }

        out = _mesh_from_blender_arrays(
            vertices.reshape(-1, 3),
            loop_vertices,
            loop_starts,
            loop_totals,
            edge_vertices=edge_vertices.reshape(-1, 2),
            point_attributes=point_attributes,
            face_attributes=face_attributes,
            custom_properties=custom_properties,
            notes=notes,
        )
    finally:
        if owner is not None:
            owner.to_mesh_clear()

    _emit("from_blender", notes)
    return out


#: ``data_type`` -> ``(foreach_get property, components, numpy dtype)``.
_ATTRIBUTE_READERS = {
    "FLOAT": ("value", 1, np.float32),
    "INT": ("value", 1, np.int32),
    "INT8": ("value", 1, np.int32),
    "BOOLEAN": ("value", 1, bool),
    "FLOAT2": ("vector", 2, np.float32),
    "FLOAT_VECTOR": ("vector", 3, np.float32),
    "FLOAT_COLOR": ("color", 4, np.float32),
}


def _read_attribute(layer, n_points, n_polygons, notes):
    """``(domain, array)`` for one Blender attribute layer, or ``None``."""
    if layer.domain not in ("POINT", "FACE"):
        notes.append(
            f"attribute '{layer.name}' is on the {layer.domain} domain and was "
            "dropped: meshio++ has point and cell data locations only"
        )
        return None
    reader = _ATTRIBUTE_READERS.get(layer.data_type)
    if reader is None:
        notes.append(
            f"attribute '{layer.name}' has type {layer.data_type} and was dropped"
        )
        return None
    prop, components, dtype = reader
    count = n_points if layer.domain == "POINT" else n_polygons
    values = np.empty(count * components, dtype=dtype)
    layer.data.foreach_get(prop, values)
    return layer.domain, (
        values if components == 1 else values.reshape(count, components)
    )
