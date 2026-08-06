"""In-memory interoperability with the wider mesh/analytics ecosystem.

Convert a :class:`~meshioplusplus._mesh.Mesh` to and from the in-memory objects of
the tools people actually use next — **without a file round-trip** — sharing the
underlying numpy buffers wherever the target accepts them as they are.

Nothing here is part of the C++ core. This is pure Python over the numpy the
readers already return, and every third-party import is lazy and named:

* **PyVista** — ``to_pyvista`` / ``from_pyvista``
  (``pip install meshioplusplus[pyvista]``). PyVista ships its own
  ``from_meshio``, but it targets the upstream ``meshio`` package and does not
  recognize ``meshioplusplus``, so it cannot be used here.
* **trimesh** — ``to_trimesh`` / ``from_trimesh``
  (``pip install meshioplusplus[trimesh]``).
* **Apache Arrow / Parquet** — ``to_arrow`` / ``from_arrow`` /
  ``write_parquet`` / ``read_parquet`` (``pip install meshioplusplus[arrow]``).
  This is a **tabular export of data arrays**, not a mesh format: it moves
  ``point_data`` / ``cell_data`` into pandas / polars / DuckDB and deliberately
  does not round-trip geometry. It is not registered in the format registry, so
  ``meshioplusplus convert mesh.vtu out.parquet`` does not and will not work.
* **pandas** — ``to_pandas`` (``pip install meshioplusplus[pandas]``). The same
  tabular export, straight to a ``pandas.DataFrame`` with no pyarrow or Parquet
  detour. pandas columns are one-dimensional, so a multi-component array is
  flattened into suffixed columns (``v`` → ``v_0``/``v_1``/``v_2``) with the
  grouping recorded in ``df.attrs`` — see :func:`to_pandas`.
* **polars** — ``to_polars`` (``pip install meshioplusplus[polars]``). The
  polars counterpart; multi-component arrays keep their true ``(n, k)`` shape
  as ``pl.Array`` columns. There is no ``from_pandas``/``from_polars``:
  :func:`from_arrow` already returns plain arrays, and a frame is one
  ``{c: df[c].to_numpy() for c in df.columns}`` away from that dict.

Architecture
------------
The module is split exactly the way :mod:`meshioplusplus._viewer` is, and for
the same reason. The bulk of it is a **pure payload layer** —
:func:`_to_vtk_payload`, :func:`_to_triangles_payload`,
:func:`_to_table_payload` — which imports no third-party library at all, does
not mutate its input, and returns plain numpy plus dataclasses. That is what
makes the genuinely subtle part (block-major indexing, the routing through
existing operations, what is shared and what is copied) testable in the default
CI matrix with none of the optional libraries installed. The public ``to_*`` /
``from_*`` functions are thin wrappers that import their target *inside* the
function and raise a named install error when it is missing.

The zero-copy contract
----------------------
A buffer is **shared** when the target accepts the array as-is: contiguous,
supported dtype, right shape. Anything else forces a copy, and every ``to_*``
takes ``zero_copy_only`` (the one exception is :func:`to_polars` — polars
always copies into its own Arrow-backed buffers, so a flag that could only
ever raise would be dishonest):

* ``False`` (default) — copies happen and each is recorded in ``notes``, which
  is surfaced as a warning.
* ``True`` — any step that would copy an array **that exists in the Mesh**
  raises :class:`ValueError`, naming the array and the reason (pyarrow's own
  ``zero_copy_only`` is the precedent).

``zero_copy_only`` governs ``points``, each ``point_data`` array, each
``cell_data`` array, and — for a single-block mesh — connectivity. It does
*not* govern **derived** arrays: VTK's ``offsets`` and ``celltypes`` have no
meshio++ counterpart at all, and a multi-block mesh has no single connectivity
array to share, so those are always constructed and merely noted. Making them
fatal would mean the flag rejected every mixed-type mesh, which is the normal
case.

Lifetime is a hard requirement, not a nicety: meshio++'s readers hand back
capsule-backed numpy owned by the C++ core, so a wrapper that shares a buffer
without holding a reference is a use-after-free. VTK's ``numpy_to_vtk(deep=0)``
does attach a ``_numpy_reference``, but this module does not rely on that — every
returned wrapper additionally carries a ``_meshioplusplus_refs`` list holding
each shared array.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np

from ._common import warn

#: dtypes VTK (and Arrow) can hold without a conversion. ``bool`` is absent on
#: purpose: VTK has no boolean array type, so a mask is widened to uint8.
_NATIVE_DTYPES = frozenset(
    (
        np.int8,
        np.uint8,
        np.int16,
        np.uint16,
        np.int32,
        np.uint32,
        np.int64,
        np.uint64,
        np.float32,
        np.float64,
    )
)

#: Prefix of the mask arrays a region is exported as, on both PyVista's
#: point_data/cell_data and trimesh's vertex/face attributes.
REGION_PREFIX = "region:"

#: Key of the JSON sidecar carrying each exported region's name/kind/dim/tag.
#: A mask array alone cannot express ``dim``/``tag``, so without this a gmsh
#: physical group's integer tag would not survive a round-trip.
REGION_META_KEY = "meshioplusplus:regions"


class InteropError(ValueError):
    """Raised when ``zero_copy_only=True`` cannot be honoured."""


# --------------------------------------------------------------------------- #
# shared helpers                                                              #
# --------------------------------------------------------------------------- #
def _copy_note(op, label, reason, notes, shared, zero_copy_only):
    """Record (or refuse) a copy of an array that exists in the mesh."""
    if zero_copy_only:
        raise InteropError(
            f"meshio++: {op}: {label} cannot be shared with the target: {reason}"
        )
    notes.append(f"{label} was copied ({reason})")
    shared[label] = False


def _derived_note(label, reason, derived):
    """Record construction of an array with no meshio++ counterpart.

    Kept out of ``notes`` on purpose: this is not a loss and not a surprise —
    it is how the target's data model differs from meshio++'s — so warning
    about it on every single call would be noise. It stays visible on the
    payload's ``derived`` list for tests and for anyone auditing what is shared.
    Never fatal under ``zero_copy_only`` — see the module docstring.
    """
    derived.append(f"{label} is derived, not shared ({reason})")


def _is_native(a):
    """Whether ``a`` can be handed over without a conversion."""
    return a.dtype.type in _NATIVE_DTYPES and a.flags["C_CONTIGUOUS"]


def _prepare_data_array(op, label, raw, notes, shared, zero_copy_only):
    """One ``point_data``/``cell_data`` array, ready to hand over.

    Returns ``None`` when the array cannot be represented at all (object dtype,
    complex), which is a drop-with-a-note rather than an error — the same policy
    the viewer uses for quantities it cannot render.
    """
    a = np.asarray(raw)
    if a.dtype == object or np.issubdtype(a.dtype, np.complexfloating):
        notes.append(f"{label} has dtype {a.dtype} and was dropped")
        return None
    if a.dtype == bool:
        _copy_note(
            op,
            label,
            "boolean arrays are widened to uint8",
            notes,
            shared,
            zero_copy_only,
        )
        return a.astype(np.uint8)
    if a.dtype.kind in "SU":
        notes.append(f"{label} is a string array and was dropped")
        return None
    if _is_native(a):
        shared[label] = True
        return a
    if a.dtype.type in _NATIVE_DTYPES:
        _copy_note(op, label, "not C-contiguous", notes, shared, zero_copy_only)
        return np.ascontiguousarray(a)
    _copy_note(
        op,
        label,
        f"dtype {a.dtype} has no native equivalent",
        notes,
        shared,
        zero_copy_only,
    )
    return np.ascontiguousarray(
        a, dtype=np.float64 if a.dtype.kind == "f" else np.int64
    )


def _flatten_cell_data(mesh, name, block_indices, notes):
    """``cell_data[name]`` for ``block_indices``, concatenated block-major.

    Returns ``None`` (with a note) when the array cannot be interpreted — one
    sub-array per cell block is the invariant, and a malformed set of blocks has
    no single flat meaning. Refusing beats guessing.
    """
    blocks = mesh.cell_data[name]
    if len(blocks) != len(mesh.cells):
        notes.append(
            f"cell_data '{name}' has {len(blocks)} block(s) for "
            f"{len(mesh.cells)} cell block(s); dropped"
        )
        return None
    parts = []
    for b in block_indices:
        blk = blocks[b]
        if blk is None:
            notes.append(f"cell_data '{name}' has an empty block; dropped")
            return None
        a = np.asarray(blk)
        if a.shape[0] != len(mesh.cells[b]):
            notes.append(
                f"cell_data '{name}' block {b} has {a.shape[0]} row(s) for the "
                f"{len(mesh.cells[b])}-cell '{mesh.cells[b].type}' block; dropped"
            )
            return None
        parts.append(a)
    if not parts:
        return None
    widths = {p.shape[1:] for p in parts}
    if len(widths) != 1:
        notes.append(
            f"cell_data '{name}' blocks disagree on component shape {sorted(map(str, widths))}; dropped"
        )
        return None
    if len(parts) == 1:
        return parts[0]
    return np.concatenate(parts, axis=0)


def _region_payload(mesh, num_points, kept_global, notes):
    """Region masks + the JSON-serializable metadata sidecar.

    ``Point`` and ``Cell`` regions become int8 0/1 masks keyed
    ``(kind, name)``; ``Side`` regions are dropped by name, because neither VTK
    nor trimesh has a ``(cell, local facet)`` concept to map them onto.

    ``kept_global`` is the ascending array of *global block-major* cell indices
    actually emitted (blocks the target cannot express are absent from it), so a
    cell region entry in a dropped block is removed here rather than silently
    landing on the wrong cell.
    """
    masks, meta = {}, []
    n_cells = len(kept_global)
    for region in getattr(mesh, "regions", []):
        if region.kind == "side":
            notes.append(
                f"side region '{region.name}' was dropped: the target has no "
                "(cell, local facet) concept"
            )
            continue
        entries = np.asarray(region.entries, dtype=np.int64).reshape(-1)
        if region.kind == "point":
            mask = np.zeros(num_points, dtype=np.int8)
            valid = entries[(entries >= 0) & (entries < num_points)]
            mask[valid] = 1
        else:
            mask = np.zeros(n_cells, dtype=np.int8)
            if n_cells and entries.size:
                pos = np.searchsorted(kept_global, entries)
                pos = np.clip(pos, 0, n_cells - 1)
                hit = kept_global[pos] == entries
                mask[pos[hit]] = 1
                if not bool(hit.all()):
                    notes.append(
                        f"cell region '{region.name}': "
                        f"{int((~hit).sum())} entrie(s) lie in dropped cell blocks"
                    )
        masks[(region.kind, region.name)] = mask
        meta.append(
            {
                "name": region.name,
                "kind": region.kind,
                "dim": int(region.dim),
                "tag": int(region.tag),
            }
        )
    return masks, meta


def _regions_from_masks(arrays, kind, meta_by_name, notes):
    """Rebuild regions from ``region:<name>`` masks, using the sidecar if present.

    ``arrays`` is consumed: every ``region:`` key is removed from it, so what is
    left is the caller's real data. Without a sidecar entry the region comes back
    with ``dim = tag = -1`` and a note — a grid built by something other than
    :func:`to_pyvista` carries the masks but not the metadata.
    """
    from ._regions import Region

    regions = []
    for key in [k for k in arrays if str(k).startswith(REGION_PREFIX)]:
        mask = np.asarray(arrays.pop(key)).reshape(-1)
        name = str(key)[len(REGION_PREFIX) :]
        entries = np.flatnonzero(mask != 0).astype(np.int64)
        info = meta_by_name.get((kind, name))
        if info is None:
            notes.append(
                f"region '{name}' was rebuilt from its mask alone; dim/tag are unknown"
            )
            info = {"dim": -1, "tag": -1}
        regions.append(
            Region(name, kind, entries, dim=int(info["dim"]), tag=int(info["tag"]))
        )
    return regions


def _emit(op, notes):
    """Surface every note as a warning, the viewer's convention."""
    for note in notes:
        warn(f"{op}: {note}")


# --------------------------------------------------------------------------- #
# payload 1: VTK arrays                                                       #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class VtkPayload:
    """Everything a VTK 9 ``vtkUnstructuredGrid`` needs, as numpy.

    ``connectivity``/``offsets``/``celltypes`` are the VTK 9 cell-array layout:
    ``offsets`` has ``num_cells + 1`` entries and ``offsets[i]:offsets[i+1]``
    slices cell ``i`` out of ``connectivity``.
    """

    #: ``(P, 3)`` float32 or float64.
    points: np.ndarray
    #: Flat int64 node indices.
    connectivity: np.ndarray
    #: int64, ``(num_cells + 1,)``, starting at 0.
    offsets: np.ndarray
    #: uint8 VTK cell type ids, ``(num_cells,)``.
    celltypes: np.ndarray
    point_data: dict = field(default_factory=dict)
    #: Concatenated block-major over the *kept* blocks.
    cell_data: dict = field(default_factory=dict)
    #: ``(kind, name) -> int8 mask``, for kinds VTK can express.
    region_masks: dict = field(default_factory=dict)
    #: ``[{name, kind, dim, tag}]`` — the sidecar.
    region_meta: list = field(default_factory=list)
    #: Ascending global (block-major) cell index of each emitted cell.
    kept_global: np.ndarray = field(default_factory=lambda: np.empty(0, np.int64))
    #: ``[(block_index, cell_type)]`` VTK cannot express.
    dropped: list = field(default_factory=list)
    #: label -> whether the mesh's own buffer is shared.
    shared: dict = field(default_factory=dict)
    #: Lossy or surprising steps, surfaced as warnings.
    notes: list = field(default_factory=list)
    #: Arrays built because the target's data model has no meshio++ counterpart.
    #: Not warnings — see :func:`_derived_note`.
    derived: list = field(default_factory=list)

    @property
    def num_cells(self) -> int:
        return len(self.celltypes)


def _vtk_type_id(cell_type):
    """VTK cell type id for a meshio++ type name, or ``None``."""
    from ._vtk_common import meshio_to_vtk_type

    if cell_type in meshio_to_vtk_type:
        return meshio_to_vtk_type[cell_type]
    # `polygon2`, `polygon7`, ... are all VTK_POLYGON; only the uniform-n-gon
    # storage shape distinguishes them, which VTK does not care about.
    if cell_type.startswith("polygon"):
        return meshio_to_vtk_type["polygon"]
    return None


def _vtk_points(mesh, notes, shared, zero_copy_only):
    pts = np.asarray(mesh.points)
    if pts.ndim == 1:
        pts = pts.reshape(-1, 1)
    n, dim = pts.shape
    if (
        dim == 3
        and pts.dtype.type in (np.float32, np.float64)
        and pts.flags["C_CONTIGUOUS"]
    ):
        shared["points"] = True
        return pts
    if dim != 3:
        reason = f"points are {dim}D and VTK needs 3 columns"
    elif not pts.flags["C_CONTIGUOUS"]:
        reason = "not C-contiguous"
    else:
        reason = f"dtype {pts.dtype} is not float32/float64"
    _copy_note("to_pyvista", "points", reason, notes, shared, zero_copy_only)
    out = np.zeros((n, 3), dtype=np.float64)
    keep = min(dim, 3)
    out[:, :keep] = pts[:, :keep]
    return out


def _vtk_block_conn(cell_block):
    """``(flat int64 connectivity, per-cell counts, shareable, reason)``."""
    from ._vtk_common import meshio_to_vtk_order

    data = cell_block.data
    perm = meshio_to_vtk_order(cell_block.type)
    if isinstance(data, np.ndarray) and data.ndim == 2:
        n, k = data.shape
        counts = np.full(n, k, dtype=np.int64)
        if perm is not None:
            # `wedge` is the one type whose meshio (gmsh Prism) node order
            # differs from VTK's; a permutation is always a copy.
            flat = np.ascontiguousarray(data[:, perm], dtype=np.int64).reshape(-1)
            return (
                flat,
                counts,
                False,
                f"'{cell_block.type}' nodes are reordered for VTK",
            )
        if data.dtype == np.int64 and data.flags["C_CONTIGUOUS"]:
            return data.reshape(-1), counts, True, None
        return (
            np.ascontiguousarray(data, dtype=np.int64).reshape(-1),
            counts,
            False,
            f"connectivity dtype {data.dtype} is not int64 (VTK's vtkIdType)",
        )
    rows = [np.asarray(r, dtype=np.int64).reshape(-1) for r in data]
    counts = np.array([len(r) for r in rows], dtype=np.int64)
    flat = np.concatenate(rows) if rows else np.empty(0, dtype=np.int64)
    return flat, counts, False, "ragged polygon rows are flattened"


def _to_vtk_payload(mesh, zero_copy_only: bool = False) -> VtkPayload:
    """Map a mesh onto the arrays a VTK 9 unstructured grid needs.

    Pure: imports no third-party library, touches no VTK object, and does not
    modify ``mesh``.

    Blocks VTK cannot express — ragged ``polyhedron``, ``custom``, anything
    absent from the meshio↔VTK type map — are **dropped with a note naming
    them**, never an exception. Mixed-type meshes are the normal case.

    Raises
    ------
    InteropError
        under ``zero_copy_only`` when an array of the mesh would be copied.
    """
    from ._regions import block_bases

    notes: list = []
    derived: list = []
    shared: dict = {}
    bases = block_bases(mesh.cells)

    kept, dropped = [], []
    for b, cb in enumerate(mesh.cells):
        type_id = _vtk_type_id(cb.type)
        if type_id is None or cb.type.startswith("polyhedron"):
            dropped.append((b, cb.type))
            continue
        kept.append((b, cb, type_id))
    if dropped:
        names = ", ".join(sorted({t for _, t in dropped}))
        notes.append(f"cell blocks VTK cannot express were dropped: {names}")

    conn_parts, count_parts, type_parts, global_parts = [], [], [], []
    share_conn, conn_reason = True, None
    for b, cb, type_id in kept:
        flat, counts, ok, reason = _vtk_block_conn(cb)
        conn_parts.append(flat)
        count_parts.append(counts)
        type_parts.append(np.full(len(counts), type_id, dtype=np.uint8))
        global_parts.append(np.arange(bases[b], bases[b] + len(cb), dtype=np.int64))
        if not ok:
            share_conn, conn_reason = False, reason

    if len(conn_parts) == 1 and share_conn:
        connectivity = conn_parts[0]
        shared["connectivity"] = True
    elif conn_parts:
        connectivity = np.concatenate(conn_parts)
        if len(conn_parts) > 1:
            _derived_note(
                "connectivity",
                f"{len(conn_parts)} cell blocks concatenated into one VTK cell array",
                derived,
            )
            shared["connectivity"] = False
        else:
            _copy_note(
                "to_pyvista", "connectivity", conn_reason, notes, shared, zero_copy_only
            )
    else:
        connectivity = np.empty(0, dtype=np.int64)
        shared["connectivity"] = False

    counts = np.concatenate(count_parts) if count_parts else np.empty(0, dtype=np.int64)
    offsets = np.zeros(len(counts) + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    celltypes = (
        np.concatenate(type_parts) if type_parts else np.empty(0, dtype=np.uint8)
    )
    kept_global = (
        np.concatenate(global_parts) if global_parts else np.empty(0, dtype=np.int64)
    )
    _derived_note(
        "offsets/celltypes", "VTK's cell array has no meshio++ counterpart", derived
    )

    points = _vtk_points(mesh, notes, shared, zero_copy_only)

    point_data = {}
    for name in sorted(mesh.point_data):
        got = _prepare_data_array(
            "to_pyvista",
            f"point_data:{name}",
            mesh.point_data[name],
            notes,
            shared,
            zero_copy_only,
        )
        if got is not None:
            point_data[name] = got

    block_indices = [b for b, _, _ in kept]
    cell_data = {}
    for name in sorted(mesh.cell_data):
        flat = _flatten_cell_data(mesh, name, block_indices, notes)
        if flat is None:
            continue
        label = f"cell_data:{name}"
        if len(block_indices) > 1:
            _derived_note(label, f"{len(block_indices)} blocks concatenated", derived)
            shared[label] = False
            got = np.ascontiguousarray(flat)
            if got.dtype == bool:
                got = got.astype(np.uint8)
            cell_data[name] = got
            continue
        got = _prepare_data_array(
            "to_pyvista", label, flat, notes, shared, zero_copy_only
        )
        if got is not None:
            cell_data[name] = got

    masks, meta = _region_payload(mesh, len(points), kept_global, notes)

    return VtkPayload(
        points=points,
        connectivity=connectivity,
        offsets=offsets,
        celltypes=celltypes,
        point_data=point_data,
        cell_data=cell_data,
        region_masks=masks,
        region_meta=meta,
        kept_global=kept_global,
        dropped=dropped,
        shared=shared,
        notes=notes,
        derived=derived,
    )


# --------------------------------------------------------------------------- #
# payload 2: triangles                                                        #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class TrianglesPayload:
    """An all-triangle surface, for targets that hold nothing else."""

    #: ``(V, 3)`` float64.
    vertices: np.ndarray
    #: ``(F, 3)`` int64.
    faces: np.ndarray
    vertex_data: dict = field(default_factory=dict)
    face_data: dict = field(default_factory=dict)
    region_masks: dict = field(default_factory=dict)
    region_meta: list = field(default_factory=list)
    #: Operations composed to get here, in order, e.g.
    #: ``["extract_surface", "linearize", "simplexify"]``.
    ops: list = field(default_factory=list)
    shared: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)
    #: See :func:`_derived_note`.
    derived: list = field(default_factory=list)


def _triangulate(mesh, notes, ops):
    """Compose the existing operations until only triangles are left.

    Never a reimplementation: volume meshes go through
    :func:`~meshioplusplus.extract_surface`, higher-order cells through
    ``convert_cells("linearize")`` and quads/polygons through
    ``convert_cells("simplexify")`` — in that order, because each stage's output
    is the next stage's input.
    """
    from ._convert_cells import _LINEAR_BASE, convert_cells
    from ._surface import extract_surface

    if max((blk.dim for blk in mesh.cells), default=0) == 3:
        notes.append(
            "the mesh has volume cells; its boundary was extracted "
            "(extract_surface drops cell_data and regions)"
        )
        mesh = extract_surface(mesh)
        ops.append("extract_surface")

    if any(blk.type in _LINEAR_BASE for blk in mesh.cells if blk.dim == 2):
        notes.append("higher-order surface cells were linearized")
        mesh = convert_cells(mesh, "linearize")
        ops.append("linearize")

    if any(blk.type == "quad" or blk.type.startswith("polygon") for blk in mesh.cells):
        notes.append("quad/polygon cells were split into triangles")
        mesh = convert_cells(mesh, "simplexify")
        ops.append("simplexify")

    return mesh


def _to_triangles_payload(mesh, zero_copy_only: bool = False) -> TrianglesPayload:
    """Map a mesh onto vertices + triangle faces.

    Pure: imports no third-party library and does not modify ``mesh``. It does
    call meshio++'s own operations, each of which returns a new mesh.

    Raises
    ------
    InteropError
        under ``zero_copy_only`` when an array of the (triangulated) mesh would
        be copied.
    ValueError
        when the mesh has no extractable triangles at all.
    """
    from ._regions import block_bases

    notes: list = []
    derived: list = []
    shared: dict = {}
    ops: list = []

    tri = _triangulate(mesh, notes, ops)
    bases = block_bases(tri.cells)

    kept = [(b, cb) for b, cb in enumerate(tri.cells) if cb.type == "triangle"]
    other = sorted({cb.type for b, cb in enumerate(tri.cells) if cb.type != "triangle"})
    if other:
        notes.append(f"non-triangle cell blocks were dropped: {', '.join(other)}")
    if not kept:
        raise ValueError(
            "meshio++: to_trimesh: the mesh has no triangles after "
            f"{' -> '.join(ops) if ops else 'no conversion'}; trimesh holds triangles only"
        )

    face_parts, global_parts = [], []
    share_faces = True
    face_reason = None
    for b, cb in kept:
        data = np.asarray(cb.data)
        if data.dtype != np.int64 or not data.flags["C_CONTIGUOUS"]:
            share_faces = False
            face_reason = f"connectivity dtype {data.dtype} is not int64"
            data = np.ascontiguousarray(data, dtype=np.int64)
        face_parts.append(data)
        global_parts.append(np.arange(bases[b], bases[b] + len(cb), dtype=np.int64))

    if len(face_parts) == 1 and share_faces:
        faces = face_parts[0]
        shared["faces"] = True
    elif len(face_parts) == 1:
        _copy_note("to_trimesh", "faces", face_reason, notes, shared, zero_copy_only)
        faces = face_parts[0]
    else:
        _derived_note(
            "faces", f"{len(face_parts)} triangle blocks concatenated", derived
        )
        shared["faces"] = False
        faces = np.concatenate(face_parts)
    kept_global = np.concatenate(global_parts)

    pts = np.asarray(tri.points)
    if pts.ndim == 1:
        pts = pts.reshape(-1, 1)
    n, dim = pts.shape
    if dim == 3 and pts.dtype == np.float64 and pts.flags["C_CONTIGUOUS"]:
        vertices = pts
        shared["vertices"] = True
    else:
        reason = (
            f"points are {dim}D and trimesh needs 3 columns"
            if dim != 3
            else f"dtype {pts.dtype} is not float64"
        )
        _copy_note("to_trimesh", "vertices", reason, notes, shared, zero_copy_only)
        vertices = np.zeros((n, 3), dtype=np.float64)
        keep = min(dim, 3)
        vertices[:, :keep] = pts[:, :keep]

    vertex_data = {}
    for name in sorted(tri.point_data):
        got = _prepare_data_array(
            "to_trimesh",
            f"point_data:{name}",
            tri.point_data[name],
            notes,
            shared,
            zero_copy_only,
        )
        if got is not None:
            vertex_data[name] = got

    block_indices = [b for b, _ in kept]
    face_data = {}
    for name in sorted(tri.cell_data):
        flat = _flatten_cell_data(tri, name, block_indices, notes)
        if flat is None:
            continue
        label = f"cell_data:{name}"
        if len(block_indices) > 1:
            _derived_note(label, f"{len(block_indices)} blocks concatenated", derived)
            shared[label] = False
            face_data[name] = np.ascontiguousarray(flat)
            continue
        got = _prepare_data_array(
            "to_trimesh", label, flat, notes, shared, zero_copy_only
        )
        if got is not None:
            face_data[name] = got

    masks, meta = _region_payload(tri, len(vertices), kept_global, notes)

    return TrianglesPayload(
        vertices=vertices,
        faces=faces,
        vertex_data=vertex_data,
        face_data=face_data,
        region_masks=masks,
        region_meta=meta,
        ops=ops,
        shared=shared,
        notes=notes,
        derived=derived,
    )


# --------------------------------------------------------------------------- #
# payload 3: a table                                                          #
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Column:
    """One table column.

    ``components > 1`` keeps its shape wherever the target can hold it — an
    Arrow ``fixed_size_list``, a polars ``pl.Array``. Only :func:`to_pandas`
    flattens it into ``_0``/``_1``/``_2`` suffix columns (pandas columns are
    one-dimensional), and there the grouping is recorded in ``df.attrs`` so the
    shape is not lost.
    """

    name: str
    #: ``(N,)`` when ``components == 1``, else ``(N, components)``.
    values: np.ndarray
    components: int = 1
    #: True for the structural columns (``x``/``y``/``z``, ``block``,
    #: ``cell_type``, ``cell``), which are built here rather than taken from the
    #: mesh — a coordinate column is a strided slice of a row-major ``points``
    #: array, so it is a copy by construction. ``zero_copy_only`` exempts them,
    #: exactly as it exempts VTK's ``offsets``/``celltypes``.
    derived: bool = False


@dataclass(frozen=True)
class TablePayload:
    """A tabular view of one data location, plus self-describing metadata."""

    location: str
    columns: list = field(default_factory=list)
    #: str -> str; becomes the Arrow schema metadata.
    metadata: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    @property
    def names(self) -> list:
        return [c.name for c in self.columns]

    @property
    def num_rows(self) -> int:
        return len(self.columns[0].values) if self.columns else 0


def _column(name, raw, notes):
    """One data array as a :class:`Column`, or ``None`` when unusable."""
    a = np.asarray(raw)
    if a.dtype == object or np.issubdtype(a.dtype, np.complexfloating):
        notes.append(f"column '{name}' has dtype {a.dtype} and was dropped")
        return None
    if a.dtype == bool:
        a = a.astype(np.uint8)
    if a.ndim == 1:
        return Column(name, np.ascontiguousarray(a), 1)
    flat = a.reshape(a.shape[0], -1)
    return Column(name, np.ascontiguousarray(flat), int(flat.shape[1]))


def _to_table_payload(
    mesh, location: str = "point", op: str = "to_arrow"
) -> TablePayload:
    """Map one data location of a mesh onto table columns.

    Pure: imports no third-party library and does not modify ``mesh``.

    ``location="point"`` builds ``x``/``y``/``z`` (as many as the mesh has) plus
    one column per ``point_data`` array. ``location="cell"`` builds ``block``,
    ``cell_type`` and ``cell`` — the **global block-major** cell index, the same
    numbering regions use — plus one column per ``cell_data`` array,
    concatenated block-major. ``op`` names the calling entry point in the
    unknown-location error, so a ``to_pandas`` mistake is not reported as a
    ``to_arrow`` one.

    Raises
    ------
    ValueError
        on an unknown ``location``.
    """
    from .__about__ import __version__
    from ._regions import block_bases

    if location not in ("point", "cell"):
        raise ValueError(
            f"meshio++: {op}: unknown location '{location}' (expected point or cell)"
        )

    notes: list = []
    columns: list = []
    bases = block_bases(mesh.cells)
    num_cells = int(bases[-1])

    if location == "point":
        pts = np.asarray(mesh.points)
        if pts.ndim == 1:
            pts = pts.reshape(-1, 1)
        for i, axis in enumerate("xyz"[: min(pts.shape[1], 3)]):
            # A column slice of a row-major array is strided, so this is a copy
            # by construction; the mesh's own layout is what makes it one.
            columns.append(
                Column(axis, np.ascontiguousarray(pts[:, i]), 1, derived=True)
            )
        if pts.shape[1] > 3:
            notes.append(f"points are {pts.shape[1]}D; only x/y/z were exported")
        for name in sorted(mesh.point_data):
            col = _column(name, mesh.point_data[name], notes)
            if col is not None:
                columns.append(col)
    else:
        counts = np.diff(bases)
        columns.append(
            Column(
                "block",
                np.repeat(np.arange(len(mesh.cells), dtype=np.int32), counts),
                1,
                derived=True,
            )
        )
        columns.append(
            Column(
                "cell_type",
                np.repeat(
                    np.array([cb.type for cb in mesh.cells], dtype=object), counts
                ),
                1,
                derived=True,
            )
        )
        columns.append(
            Column("cell", np.arange(num_cells, dtype=np.int64), 1, derived=True)
        )
        block_indices = list(range(len(mesh.cells)))
        for name in sorted(mesh.cell_data):
            flat = _flatten_cell_data(mesh, name, block_indices, notes)
            if flat is None:
                continue
            col = _column(name, flat, notes)
            if col is not None:
                columns.append(col)

    metadata = {
        "meshioplusplus:version": str(__version__),
        "meshioplusplus:location": location,
        "meshioplusplus:num_points": str(len(mesh.points)),
        "meshioplusplus:num_cells": str(num_cells),
        "meshioplusplus:cell_types": json.dumps(
            [[cb.type, len(cb)] for cb in mesh.cells], separators=(",", ":")
        ),
        REGION_META_KEY: json.dumps(
            [
                {
                    "name": r.name,
                    "kind": r.kind,
                    "dim": int(r.dim),
                    "tag": int(r.tag),
                }
                for r in getattr(mesh, "regions", [])
            ],
            separators=(",", ":"),
        ),
    }
    return TablePayload(
        location=location, columns=columns, metadata=metadata, notes=notes
    )


# --------------------------------------------------------------------------- #
# availability predicates                                                     #
# --------------------------------------------------------------------------- #
def _importable(module: str) -> bool:
    try:
        __import__(module)
    except Exception:
        return False
    return True


def has_pyvista() -> bool:
    """Whether PyVista is installed (``pip install meshioplusplus[pyvista]``)."""
    return _importable("pyvista")


def has_trimesh() -> bool:
    """Whether trimesh is installed (``pip install meshioplusplus[trimesh]``)."""
    return _importable("trimesh")


def has_arrow() -> bool:
    """Whether pyarrow is installed (``pip install meshioplusplus[arrow]``)."""
    return _importable("pyarrow")


def has_pandas() -> bool:
    """Whether pandas is installed (``pip install meshioplusplus[pandas]``)."""
    return _importable("pandas")


def has_polars() -> bool:
    """Whether polars is installed (``pip install meshioplusplus[polars]``)."""
    return _importable("polars")


def has_open3d() -> bool:
    """Always ``False``: the Open3D bridge is Phase 2. See ``doc/interop.md``."""
    return False


def has_dolfinx() -> bool:
    """Always ``False``: the DOLFINx bridge is Phase 2. See ``doc/interop.md``."""
    return False


def _require(module: str, extra: str, op: str):
    try:
        return __import__(module)
    except ImportError as e:
        raise ImportError(
            f"meshio++: {op}: {module} is not installed; "
            f"install it with `pip install meshioplusplus[{extra}]`"
        ) from e


# --------------------------------------------------------------------------- #
# PyVista                                                                     #
# --------------------------------------------------------------------------- #
def to_pyvista(mesh, zero_copy_only: bool = False):
    """Convert a mesh to a :class:`pyvista.UnstructuredGrid`, sharing buffers.

    PyVista ships its own ``from_meshio``, but it targets the upstream ``meshio``
    package and does not recognize ``meshioplusplus`` — hence this function.

    ``point_data``/``cell_data`` land on ``grid.point_data``/``grid.cell_data``,
    the latter concatenated block-major (the numbering regions use). Cell blocks
    VTK cannot express (ragged ``polyhedron``, ``custom``, anything absent from
    the meshio↔VTK map) are dropped with a warning naming them.

    ``Point`` and ``Cell`` regions are exported as int8 masks named
    ``region:<name>``, plus a JSON sidecar in
    ``grid.field_data["meshioplusplus:regions"]`` carrying each region's
    ``dim``/``tag`` — a mask alone cannot express those. ``Side`` regions have no
    VTK equivalent (VTK has no ``(cell, local facet)`` concept) and are dropped
    with a warning.

    Parameters
    ----------
    mesh :
        the mesh to convert (unmodified).
    zero_copy_only :
        when true, raise instead of copying any array of the mesh. VTK's
        ``offsets``/``celltypes`` and multi-block concatenation are derived
        arrays and remain exempt — see the module docstring.

    Returns
    -------
    pyvista.UnstructuredGrid
        holding references to every shared numpy array in
        ``grid._meshioplusplus_refs``, so the grid stays valid after the source
        mesh is garbage-collected.

    Raises
    ------
    ImportError
        when PyVista is not installed.
    InteropError
        under ``zero_copy_only`` when an array would be copied.
    """
    pv = _require("pyvista", "pyvista", "to_pyvista")
    import vtk
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

    payload = _to_vtk_payload(mesh, zero_copy_only=zero_copy_only)
    _emit("to_pyvista", payload.notes)

    refs = [payload.points, payload.connectivity, payload.offsets, payload.celltypes]

    # pyvista's UnstructuredGrid(cells, celltypes, points) constructor takes the
    # *legacy* `[n, p0, ..., n, p1, ...]` layout and always copies. Going through
    # vtkCellArray.SetData is what actually shares the VTK 9 offsets/connectivity
    # pair, so this route is mandatory rather than a preference.
    cell_array = vtk.vtkCellArray()
    cell_array.SetData(
        numpy_to_vtkIdTypeArray(payload.offsets, deep=False),
        numpy_to_vtkIdTypeArray(payload.connectivity, deep=False),
    )
    grid = pv.UnstructuredGrid()
    grid.points = payload.points
    grid.SetCells(
        numpy_to_vtk(payload.celltypes, deep=False, array_type=vtk.VTK_UNSIGNED_CHAR),
        cell_array,
    )

    for name, values in payload.point_data.items():
        grid.point_data[name] = values
        refs.append(values)
    for name, values in payload.cell_data.items():
        grid.cell_data[name] = values
        refs.append(values)
    for (kind, name), mask in payload.region_masks.items():
        target = grid.point_data if kind == "point" else grid.cell_data
        target[REGION_PREFIX + name] = mask
        refs.append(mask)
    if payload.region_meta:
        grid.field_data[REGION_META_KEY] = np.array(
            [json.dumps(payload.region_meta, separators=(",", ":"))]
        )

    grid._meshioplusplus_refs = refs
    return grid


def from_pyvista(grid):
    """Convert a :class:`pyvista.UnstructuredGrid` back to a mesh.

    The inverse of :func:`to_pyvista`. ``region:<name>`` mask arrays are turned
    back into regions, using the ``meshioplusplus:regions`` sidecar for each
    region's ``dim``/``tag`` when present; a grid built by something else carries
    the masks but not the sidecar, and its regions come back with
    ``dim = tag = -1`` and a warning.

    Raises
    ------
    ImportError
        when PyVista is not installed.
    """
    _require("pyvista", "pyvista", "from_pyvista")

    from ._mesh import Mesh
    from ._vtk_common import vtk_cells_from_data

    notes: list = []
    connectivity = np.asarray(grid.cell_connectivity)
    # `vtk_cells_from_data` wants END offsets (the VTU convention); a VTK 9 cell
    # array's `offset` starts at 0 and has num_cells + 1 entries.
    offsets = np.asarray(grid.offset)[1:]
    celltypes = np.asarray(grid.celltypes)

    meta_by_name = {}
    try:
        raw_meta = grid.field_data[REGION_META_KEY]
    except KeyError:
        raw_meta = None
    if raw_meta is not None:
        try:
            entries = json.loads(np.asarray(raw_meta).reshape(-1)[0])
            meta_by_name = {(e["kind"], e["name"]): e for e in entries}
        except Exception:
            notes.append(
                f"the '{REGION_META_KEY}' sidecar is unreadable and was ignored"
            )

    # Cell region masks are indexed by VTK cell, which is exactly the global
    # block-major numbering `vtk_cells_from_data` produces (it splits the type
    # array into contiguous runs, in order). So they must be harvested from the
    # FLAT arrays, before that call turns cell data into per-block lists.
    cell_flat = {str(k): np.asarray(v) for k, v in grid.cell_data.items()}
    point_data = {str(k): np.asarray(v) for k, v in grid.point_data.items()}
    regions = _regions_from_masks(point_data, "point", meta_by_name, notes)
    regions += _regions_from_masks(cell_flat, "cell", meta_by_name, notes)

    cells, cell_data = vtk_cells_from_data(connectivity, offsets, celltypes, cell_flat)

    _emit("from_pyvista", notes)
    return Mesh(
        np.asarray(grid.points),
        cells,
        point_data=point_data,
        cell_data=cell_data,
        regions=regions,
    )


# --------------------------------------------------------------------------- #
# trimesh                                                                     #
# --------------------------------------------------------------------------- #
def to_trimesh(mesh, zero_copy_only: bool = False):
    """Convert a mesh to a :class:`trimesh.Trimesh`, sharing buffers.

    trimesh holds **triangles only**, so non-triangle input is routed through
    meshio++'s own operations rather than a reimplementation: volume meshes
    through :func:`~meshioplusplus.extract_surface`, higher-order cells through
    ``convert_cells("linearize")`` and quads/polygons through
    ``convert_cells("simplexify")``. Each step is warned about.

    Regions follow those operations' own documented behaviour and get no second
    policy here: a pure-triangle input keeps its ``Point``/``Cell`` regions
    (as ``region:<name>`` masks on ``vertex_attributes``/``face_attributes``),
    while anything routed through ``extract_surface`` loses them, because
    ``extract_surface`` drops regions. ``Side`` regions are always dropped.

    Parameters
    ----------
    mesh :
        the mesh to convert (unmodified).
    zero_copy_only :
        when true, raise instead of copying any array of the (triangulated)
        mesh.

    Returns
    -------
    trimesh.Trimesh
        holding references to every shared numpy array in
        ``tm._meshioplusplus_refs``.

    Raises
    ------
    ImportError
        when trimesh is not installed.
    InteropError
        under ``zero_copy_only`` when an array would be copied.
    ValueError
        when the mesh yields no triangles at all.
    """
    trimesh = _require("trimesh", "trimesh", "to_trimesh")

    payload = _to_triangles_payload(mesh, zero_copy_only=zero_copy_only)
    _emit("to_trimesh", payload.notes)

    # process=False is mandatory: the default merges vertices, which silently
    # renumbers the mesh and would break every index-based attribute below.
    tm = trimesh.Trimesh(vertices=payload.vertices, faces=payload.faces, process=False)
    refs = [payload.vertices, payload.faces]

    from ._viewer import _COLOR_NAME_RE

    for name, values in payload.vertex_data.items():
        tm.vertex_attributes[name] = values
        refs.append(values)
        # Colours are recognized by NAME, never inferred from a [0, 1] range --
        # the repo-wide rule (see _viewer._looks_like_color).
        if (
            _COLOR_NAME_RE.search(name)
            and values.ndim == 2
            and values.shape[1] in (3, 4)
        ):
            tm.visual = trimesh.visual.ColorVisuals(mesh=tm, vertex_colors=values)
    for name, values in payload.face_data.items():
        tm.face_attributes[name] = values
        refs.append(values)
    for (kind, name), mask in payload.region_masks.items():
        target = tm.vertex_attributes if kind == "point" else tm.face_attributes
        target[REGION_PREFIX + name] = mask
        refs.append(mask)
    if payload.region_meta:
        tm.metadata[REGION_META_KEY] = payload.region_meta

    tm._meshioplusplus_refs = refs
    return tm


def from_trimesh(tm):
    """Convert a :class:`trimesh.Trimesh` back to a mesh.

    Raises
    ------
    ImportError
        when trimesh is not installed.
    """
    _require("trimesh", "trimesh", "from_trimesh")

    from ._mesh import Mesh

    notes: list = []
    point_data = {str(k): np.asarray(v) for k, v in dict(tm.vertex_attributes).items()}
    cell_data_flat = {
        str(k): np.asarray(v) for k, v in dict(tm.face_attributes).items()
    }

    meta_by_name = {}
    raw_meta = (tm.metadata or {}).get(REGION_META_KEY)
    if raw_meta:
        meta_by_name = {(e["kind"], e["name"]): e for e in raw_meta}

    regions = _regions_from_masks(point_data, "point", meta_by_name, notes)
    regions += _regions_from_masks(cell_data_flat, "cell", meta_by_name, notes)

    faces = np.ascontiguousarray(np.asarray(tm.faces), dtype=np.int64)
    cell_data = {name: [values] for name, values in cell_data_flat.items()}

    _emit("from_trimesh", notes)
    return Mesh(
        np.asarray(tm.vertices),
        [("triangle", faces)],
        point_data=point_data,
        cell_data=cell_data,
        regions=regions,
    )


# --------------------------------------------------------------------------- #
# Arrow / Parquet                                                             #
# --------------------------------------------------------------------------- #
def to_arrow(mesh, location: str = "point", zero_copy_only: bool = False):
    """Export one data location of a mesh as a :class:`pyarrow.Table`.

    **This is not a mesh format.** It moves ``point_data``/``cell_data`` into the
    analytics stack (pandas, polars, DuckDB); it does not round-trip geometry,
    and it is deliberately not registered in meshio++'s format registry, so
    ``meshioplusplus convert mesh.vtu out.parquet`` does not work.

    Multi-component arrays become Arrow ``fixed_size_list`` columns rather than
    ``_0``/``_1``/``_2`` suffix columns, which would lose the shape. The mesh's
    counts, cell types, version, location and region names ride along in the
    schema metadata, so a table is self-describing.

    Parameters
    ----------
    mesh :
        the mesh to export (unmodified).
    location :
        ``"point"`` (default) or ``"cell"``.
    zero_copy_only :
        when true, raise instead of copying a numeric buffer. Note that the
        ``x``/``y``/``z`` columns are column slices of a row-major ``points``
        array and are therefore always materialized.

    Raises
    ------
    ImportError
        when pyarrow is not installed.
    """
    pa = _require("pyarrow", "arrow", "to_arrow")

    payload = _to_table_payload(mesh, location, op="to_arrow")

    arrays, names = [], []
    for column in payload.columns:
        values = column.values
        if values.dtype == object:
            array = pa.array(list(values))
        elif column.components > 1:
            array = pa.FixedSizeListArray.from_arrays(
                pa.array(values.reshape(-1)), column.components
            )
        else:
            array = pa.array(values)
        # Measured, not assumed: Arrow's zero-copy for a contiguous numeric
        # buffer is real, but a string column or a masked type quietly is not.
        if not column.derived and not _arrow_shares(pa, array, values):
            _copy_note(
                "to_arrow",
                f"column '{column.name}'",
                f"Arrow materialized its {values.dtype} buffer",
                payload.notes,
                {},
                zero_copy_only,
            )
        arrays.append(array)
        names.append(column.name)
    _emit("to_arrow", payload.notes)

    schema = pa.schema(
        [pa.field(n, a.type) for n, a in zip(names, arrays)], metadata=payload.metadata
    )
    return pa.Table.from_arrays(arrays, schema=schema)


def _arrow_shares(pa, array, values) -> bool:
    """Whether ``array`` actually points at ``values``' buffer."""
    try:
        flat = array.flatten() if pa.types.is_fixed_size_list(array.type) else array
        return bool(np.shares_memory(flat.to_numpy(zero_copy_only=True), values))
    except Exception:
        return False


def from_arrow(table) -> dict:
    """Read a :class:`pyarrow.Table` back into a dict of numpy arrays.

    ``fixed_size_list`` columns come back with their ``(N, components)`` shape.
    This is the inverse of :func:`to_arrow`'s *data*; it does not reconstruct a
    mesh, because a table never carried the geometry.

    Raises
    ------
    ImportError
        when pyarrow is not installed.
    """
    pa = _require("pyarrow", "arrow", "from_arrow")

    out = {}
    for name, column in zip(table.schema.names, table.columns):
        array = _single_array(column)
        if pa.types.is_fixed_size_list(array.type):
            width = array.type.list_size
            out[name] = _to_numpy(array.flatten()).reshape(-1, width)
        else:
            out[name] = _to_numpy(array)
    return out


def _single_array(column):
    """One contiguous Arrow array out of a table column.

    ``ChunkedArray.combine_chunks()`` copies even when there is a single chunk
    (measured against pyarrow 25), which would silently defeat the zero-copy
    claim on the common one-chunk case — so take the chunk directly when we can.
    """
    if not hasattr(column, "num_chunks"):
        return column
    if column.num_chunks == 1:
        return column.chunk(0)
    return column.combine_chunks()


def _to_numpy(array):
    """Numpy view of an Arrow array, sharing its buffer when Arrow allows."""
    try:
        return array.to_numpy(zero_copy_only=True)
    except Exception:
        # Strings, nulls, and anything else with no flat numeric buffer.
        return array.to_numpy(zero_copy_only=False)


def write_parquet(mesh, path, location: str = "point", **kwargs) -> None:
    """Write one data location of a mesh to a Parquet file.

    A thin wrapper over :func:`to_arrow` — same caveat: this exports **data for
    analytics**, not a mesh. Extra keyword arguments go to
    ``pyarrow.parquet.write_table`` (``compression=``, ...).

    Raises
    ------
    ImportError
        when pyarrow is not installed.
    """
    _require("pyarrow", "arrow", "write_parquet")
    import pyarrow.parquet as pq

    pq.write_table(to_arrow(mesh, location=location), str(path), **kwargs)


def read_parquet(path, return_metadata: bool = False):
    """Read a Parquet file written by :func:`write_parquet`.

    Returns a dict of numpy arrays; with ``return_metadata=True`` returns
    ``(arrays, metadata)`` where ``metadata`` is the decoded
    ``meshioplusplus:*`` schema metadata.

    Raises
    ------
    ImportError
        when pyarrow is not installed.
    """
    _require("pyarrow", "arrow", "read_parquet")
    import pyarrow.parquet as pq

    table = pq.read_table(str(path))
    arrays = from_arrow(table)
    if not return_metadata:
        return arrays
    raw = table.schema.metadata or {}
    metadata = {
        k.decode() if isinstance(k, bytes) else k: (
            v.decode() if isinstance(v, bytes) else v
        )
        for k, v in raw.items()
    }
    return arrays, metadata


# --------------------------------------------------------------------------- #
# pandas / polars                                                             #
# --------------------------------------------------------------------------- #
#: ``df.attrs`` key recording the pandas suffix expansion: compact JSON
#: ``{source_name: [suffixed column names, in component order]}``. Always
#: present (``"{}"`` when nothing was expanded), so the grouping — and with it
#: the original ``(n, k)`` shape — is recoverable from the frame alone.
COMPONENTS_META_KEY = "meshioplusplus:components"


def _frame_columns(payload):
    """Expand multi-component columns into suffixed one-component columns.

    Pure: imports no third-party library. This is the part of :func:`to_pandas`
    with actual decisions in it — pandas columns are one-dimensional, so a
    ``(n, k)`` column ``v`` becomes ``v_0`` .. ``v_{k-1}`` — kept separate so it
    is testable in the default CI matrix with pandas absent.

    Returns ``(columns, groups)``: ``columns`` is a list of one-component
    :class:`Column` in payload order (scalar columns pass through by identity,
    no copy at this layer), ``groups`` maps each expanded source name to its
    suffixed column names. Duplicate final names (a scalar ``v_1`` next to a
    vector ``v``, a ``point_data["x"]`` next to the derived ``x``) keep the
    first occurrence and drop the rest with a note appended to
    ``payload.notes`` — pandas cannot hold two columns of one name.
    """
    seen: set = set()
    columns: list = []
    groups: dict = {}
    for column in payload.columns:
        if column.components == 1:
            if column.name in seen:
                payload.notes.append(
                    f"column '{column.name}' was dropped "
                    "(name collides with a column already produced)"
                )
                continue
            seen.add(column.name)
            columns.append(column)
            continue
        emitted = []
        for i in range(column.components):
            name = f"{column.name}_{i}"
            if name in seen:
                payload.notes.append(
                    f"column '{name}' was dropped "
                    "(name collides with a column already produced)"
                )
                continue
            seen.add(name)
            # A component is a strided slice of a row-major array — a copy by
            # construction, the x/y/z argument (but these arrays DO exist in
            # the mesh, so to_pandas records the copy; x/y/z are exempt).
            columns.append(
                Column(
                    name,
                    np.ascontiguousarray(column.values[:, i]),
                    1,
                    derived=column.derived,
                )
            )
            emitted.append(name)
        groups[column.name] = emitted
    return columns, groups


def to_pandas(mesh, location: str = "point", zero_copy_only: bool = False):
    """Export one data location of a mesh as a :class:`pandas.DataFrame`.

    **This is not a mesh format.** Like :func:`to_arrow`, it moves
    ``point_data``/``cell_data`` into the analytics stack without a file
    round-trip; it does not carry geometry, and there is no ``from_pandas`` —
    :func:`from_arrow` already returns plain arrays, and a frame is one
    ``{c: df[c].to_numpy() for c in df.columns}`` away from that dict.

    pandas columns are one-dimensional, so — unlike Arrow's
    ``fixed_size_list`` and polars' ``pl.Array`` — a multi-component array
    ``v`` is flattened into suffixed columns ``v_0``/``v_1``/``v_2``. The
    grouping is recorded in ``df.attrs["meshioplusplus:components"]`` (compact
    JSON, source name → suffixed names) alongside the full ``meshioplusplus:*``
    metadata :func:`to_arrow` puts in the schema, so nothing is lost — but note
    pandas drops ``attrs`` on many operations.

    Parameters
    ----------
    mesh :
        the mesh to export (unmodified).
    location :
        ``"point"`` (default) or ``"cell"``.
    zero_copy_only :
        when true, raise instead of copying a numeric buffer that exists in
        the mesh. A multi-component array **always** copies here (the suffix
        expansion slices a row-major array), so under this flag it raises —
        :func:`to_arrow` is the shared-buffer path for vector data. Scalar
        sharing is measured, never assumed.

    Raises
    ------
    ImportError
        when pandas is not installed.
    InteropError
        under ``zero_copy_only=True`` when a mesh array had to be copied.
    """
    pd = _require("pandas", "pandas", "to_pandas")

    payload = _to_table_payload(mesh, location, op="to_pandas")
    columns, groups = _frame_columns(payload)
    for name, emitted in groups.items():
        _copy_note(
            "to_pandas",
            f"column '{name}'",
            "pandas frames hold 1-D columns only; expanding to "
            f"{', '.join(emitted)} copies (to_arrow keeps the shape)",
            payload.notes,
            {},
            zero_copy_only,
        )
    expanded = {name for emitted in groups.values() for name in emitted}

    # copy=False is load-bearing: pandas >= 2.0 copies a dict of ndarrays by
    # default. No _meshioplusplus_refs list is needed here, unlike PyVista —
    # a shared column is a numpy view whose .base chain keeps the mesh's
    # capsule-backed buffer alive by ordinary refcounting.
    df = pd.DataFrame({c.name: c.values for c in columns}, copy=False)
    for column in columns:
        if column.derived or column.name in expanded:
            continue
        if column.values.dtype == object:
            continue
        # Measured, not assumed — the _arrow_shares precedent.
        if not np.shares_memory(df[column.name].to_numpy(), column.values):
            _copy_note(
                "to_pandas",
                f"column '{column.name}'",
                f"pandas materialized its {column.values.dtype} buffer",
                payload.notes,
                {},
                zero_copy_only,
            )
    _emit("to_pandas", payload.notes)

    df.attrs = dict(payload.metadata)
    df.attrs[COMPONENTS_META_KEY] = json.dumps(groups, separators=(",", ":"))
    return df


def to_polars(mesh, location: str = "point"):
    """Export one data location of a mesh as a :class:`polars.DataFrame`.

    **This is not a mesh format** — :func:`to_arrow`'s caveat verbatim, and
    like pandas there is deliberately no ``from_polars``.

    Multi-component arrays keep their true ``(n, k)`` shape as ``pl.Array``
    columns; the derived ``cell_type`` strings become ``pl.String``. There is
    **no** ``zero_copy_only`` parameter — polars always copies into its own
    Arrow-backed buffers, so the frame is independent of the mesh by
    construction and a flag that could only ever raise would be dishonest
    (use :func:`to_arrow` when you need shared buffers). Polars also has no
    ``attrs``/schema-metadata slot, so the ``meshioplusplus:*`` metadata does
    not ride along; :func:`to_arrow` is the self-describing-table path.

    Raises
    ------
    ImportError
        when polars is not installed.
    """
    pl = _require("polars", "polars", "to_polars")

    payload = _to_table_payload(mesh, location, op="to_polars")
    series = []
    for column in payload.columns:
        if column.values.dtype == object:
            series.append(pl.Series(column.name, list(column.values), dtype=pl.String))
        else:
            # A (n, k) numpy array maps to pl.Array(inner, k) — the shape
            # survives structurally, which is why polars needs no suffix rule.
            series.append(pl.Series(column.name, column.values))
    _emit("to_polars", payload.notes)
    return pl.DataFrame(series)


# --------------------------------------------------------------------------- #
# Phase 2                                                                     #
# --------------------------------------------------------------------------- #
def _phase_two(name: str, reason: str):
    raise NotImplementedError(
        f"meshio++: {name} is Phase 2 and not implemented yet: {reason}. "
        "See doc/interop.md."
    )


def to_open3d(mesh, zero_copy_only: bool = False):
    """Phase 2. See ``doc/interop.md`` for the design and its constraints."""
    _phase_two(
        "to_open3d",
        "Open3D's Vector3dVector typically copies, so the zero-copy contract "
        "differs per structure, and the wheel is ~400 MB",
    )


def to_dolfinx(mesh, **kwargs):
    """Phase 2. See ``doc/interop.md`` for the design and its constraints."""
    _phase_two(
        "to_dolfinx",
        "DOLFINx needs a single-cell-type mesh, a ufl/basix domain, an MPI "
        "communicator and a VTK->basix node-ordering permutation, and is "
        "conda/apt-only so it can never be a pip extra",
    )
