"""Interactive visualization of a mesh, with two interchangeable backends.

A dependency-free *mapping* plus two optional renderers. Nothing here reads or
writes a file format, and nothing here is part of the C++ core — meshio++
already knows how to produce every array a renderer needs, so the job is to
hand those arrays over rather than to draw anything:

* **polyscope** — a native desktop window, and headless screenshots. A
  Python-only optional dependency (``pip install meshioplusplus[viewer]``); it
  is never imported at module scope and never reaches C++/WASM/C/Fortran.
* **browser** — the mesh is written to VTP in memory and rendered by the
  vtk.js app under ``src/viewer/``, either inline in a notebook or in the default
  browser. See :mod:`meshioplusplus._viewer_browser`.

The bulk of this module is :func:`_to_polyscope_payload`, which is *pure*: it
takes a :class:`~meshioplusplus._mesh.Mesh` and returns plain numpy arrays and
dataclasses, importing no renderer at all. That is what makes the routing
logic — which is the only genuinely subtle part — testable without a display.

Two invariants govern that routing, and getting either wrong produces a mesh
that renders but is silently coloured wrong:

**1. polyscope volume meshes hold only tetrahedra and hexahedra.**
``register_volume_mesh`` takes ``tets``/``hexes``/``mixed_cells`` and nothing
else, so wedges, pyramids and the higher-order family must go through
:func:`~meshioplusplus.convert_cells` in ``"simplexify"`` mode first. Hexes
must *not* — decomposing them would be a gratuitous loss of the element
structure the user is usually trying to look at.

**2.** ``convert_cells`` **preserves cell_data block-for-block;**
``extract_surface`` **destroys it.** So after a ``convert_cells`` call, re-enter
the pipeline on the converted mesh and read *its* ``cell_data``; after an
``extract_surface`` call, pass ``record_parent_ids=True`` and gather from the
**original** mesh's flattened ``cell_data`` through ``surface:parent_cell``.

The mechanism tying cell data to rendered primitives is the ``cell_source``
array built by :func:`_dimension_rows`: the index, in the *global* cell
numbering over all blocks, of the input cell behind each emitted primitive.
A gmsh file whose ``triangle`` block is preceded by a boundary-marker ``line``
block is the case that makes this necessary — a naive per-dimension
concatenation offsets every cell quantity by the length of the skipped block.

Public API:

* :func:`view` — show a mesh, dispatching on ``backend``.
* :func:`screenshot` — render a mesh to a PNG without opening a window.
* :func:`has_viewer` — whether the polyscope backend is importable.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ._common import warn

#: Volume cell types polyscope can hold directly, and their corner counts.
#: Everything else 3D is simplexified into tetrahedra first.
_POLYSCOPE_VOLUME_TYPES = {"tetra": 4, "hexahedron": 8}

#: Width of a polyscope ``mixed_cells`` row (a hexahedron); tetrahedra pad the
#: remaining slots with a negative value.
_MIXED_CELL_WIDTH = 8

#: Names that mean "these three components are a colour, not a vector". Matched
#: against the array name, *not* inferred from the value range: a normalized
#: displacement field or a unit normal shifted into [0, 1] would otherwise be
#: silently rendered as RGB with no way for the user to tell.
_COLOR_NAME_RE = re.compile(r"(?:^(?:colou?rs?|rgb)$)|(?:[_:](?:colou?rs?|rgb)$)", re.I)

#: Components above which per-component scalars stop being useful and only a
#: magnitude is emitted.
_MAX_COMPONENT_QUANTITIES = 16


@dataclass(frozen=True)
class Quantity:
    """One renderable data array, already reduced to what polyscope accepts."""

    #: Quantity name, after collision resolution.
    name: str
    #: ``"scalar"``, ``"vector"`` or ``"color"``.
    kind: str
    #: ``(N,)`` for a scalar, ``(N, 3)`` for a vector/colour. Always C-contiguous
    #: float64 — polyscope's binding layer would otherwise copy on every frame.
    values: np.ndarray
    #: ``"vertices"``, ``"faces"``, ``"cells"`` or ``"edges"``.
    defined_on: str
    #: Explicit colormap range. Set only when the array holds non-finite
    #: entries, where polyscope's autoscale would otherwise poison the whole
    #: colormap from a single NaN.
    vminmax: Optional[tuple] = None
    #: ``"point_data"`` or ``"cell_data"`` — provenance, for tests and errors.
    source: str = "point_data"


@dataclass(frozen=True)
class PolyscopePayload:
    """Everything needed to register one polyscope structure.

    Exactly one of ``faces``, ``edges``, ``{tets, hexes}`` and ``mixed_cells``
    is populated, selected by ``kind``.
    """

    #: ``"surface"``, ``"volume"``, ``"curve"`` or ``"points"``.
    kind: str
    #: ``(P, 3)`` float64 — always three columns, whatever the input dimension.
    vertices: np.ndarray
    #: Ragged list of index lists, for ``kind == "surface"``.
    faces: Optional[list] = None
    #: ``(E, 2)`` int32, for ``kind == "curve"``.
    edges: Optional[np.ndarray] = None
    #: ``(T, 4)`` int32.
    tets: Optional[np.ndarray] = None
    #: ``(H, 8)`` int32.
    hexes: Optional[np.ndarray] = None
    #: ``(M, 8)`` int32, negative-padded on tetrahedron rows.
    mixed_cells: Optional[np.ndarray] = None
    quantities: list = field(default_factory=list)
    #: Human-readable log of every lossy or surprising step taken. Empty when
    #: the mesh was rendered exactly as given.
    notes: list = field(default_factory=list)

    @property
    def num_primitives(self) -> int:
        """How many faces/cells/edges this payload draws."""
        for a in (self.faces, self.edges, self.tets, self.hexes, self.mixed_cells):
            if a is not None:
                return len(a)
        return 0


# --- geometry ------------------------------------------------------------- #


def _vertices_of(mesh, notes: list) -> np.ndarray:
    """``mesh.points`` as ``(P, 3)`` float64, padded or truncated as needed."""
    pts = np.asarray(mesh.points, dtype=np.float64)
    if pts.ndim == 1:
        pts = pts.reshape(-1, 1)
    n, dim = pts.shape
    if dim == 3:
        return np.ascontiguousarray(pts)
    if dim < 3:
        out = np.zeros((n, 3), dtype=np.float64)
        out[:, :dim] = pts
        notes.append(f"points are {dim}D; padded to 3D with zeros for display")
        return out
    notes.append(f"points are {dim}D; only the first 3 coordinates are displayed")
    return np.ascontiguousarray(pts[:, :3])


def _dimension_rows(mesh, dim: int):
    """Split ``mesh.cells`` into the blocks of topological dimension ``dim``.

    Returns ``(blocks, cell_source)`` where ``blocks`` are the matching
    :class:`CellBlock`\\ s in mesh order and ``cell_source`` is an int64 array
    giving, per emitted cell, its index in the *global* numbering across all
    blocks. The counter advances over skipped blocks too — that is the whole
    point, and is what keeps cell data aligned when a mesh carries lower-
    dimensional blocks (gmsh boundary markers) alongside the ones being drawn.
    """
    blocks, sources, g = [], [], 0
    for blk in mesh.cells:
        n = len(blk)
        if blk.dim == dim:
            blocks.append(blk)
            sources.append(np.arange(g, g + n, dtype=np.int64))
        g += n
    source = np.concatenate(sources) if sources else np.empty(0, dtype=np.int64)
    return blocks, source


def _faces_list(blocks) -> list:
    """A ragged list-of-lists of face indices, in block order.

    polyscope accepts a ragged list directly, so mixed triangle/quad blocks and
    ragged polygon blocks need no triangulation — which would otherwise add
    edges the user never asked to see.
    """
    faces: list = []
    for blk in blocks:
        data = blk.data
        if isinstance(data, np.ndarray) and data.ndim == 2:
            faces.extend(data.astype(np.int64, copy=False).tolist())
        else:  # ragged polygon block: a list of variable-length rows
            faces.extend([[int(i) for i in row] for row in data])
    return faces


def _mixed_cells(blocks) -> np.ndarray:
    """Pack 3D blocks into a negative-padded ``(M, 8)`` polyscope array.

    Built here rather than by handing polyscope separate ``tets=``/``hexes=``
    arrays, because then the row order — and so the meaning of every per-cell
    quantity — would depend on polyscope's internal merge order. Owning it
    makes the cell-data mapping correct by construction.
    """
    rows = []
    for blk in blocks:
        conn = np.asarray(blk.data, dtype=np.int32)
        pad = np.full((len(conn), _MIXED_CELL_WIDTH), -1, dtype=np.int32)
        pad[:, : conn.shape[1]] = conn
        rows.append(pad)
    if not rows:
        return np.empty((0, _MIXED_CELL_WIDTH), dtype=np.int32)
    return np.ascontiguousarray(np.concatenate(rows))


def _stack_conn(blocks) -> np.ndarray:
    """Concatenate same-type blocks into one int32 connectivity array."""
    if not blocks:
        return np.empty((0, 0), dtype=np.int32)
    return np.ascontiguousarray(
        np.concatenate([np.asarray(b.data, dtype=np.int32) for b in blocks])
    )


# --- data ----------------------------------------------------------------- #


def _flatten_cell_data(mesh, name: str, notes: list):
    """``cell_data[name]`` as one array indexed by *global* cell index.

    Returns ``None`` (with a note) when the array cannot be interpreted: meshio++
    stores one sub-array per cell block, and a malformed or inconsistent set of
    blocks has no single flat meaning. Refusing beats guessing — a wrong answer
    here is a plausible-looking but wrongly-coloured picture.
    """
    blocks = mesh.cell_data[name]
    if len(blocks) != len(mesh.cells):
        notes.append(
            f"cell_data '{name}' has {len(blocks)} block(s) for "
            f"{len(mesh.cells)} cell block(s); dropped"
        )
        return None
    parts = []
    for blk, cb in zip(blocks, mesh.cells):
        if blk is None:
            notes.append(f"cell_data '{name}' has an empty block; dropped")
            return None
        a = np.asarray(blk)
        if a.shape[0] != len(cb):
            notes.append(
                f"cell_data '{name}' block length {a.shape[0]} does not match "
                f"the {len(cb)}-cell '{cb.type}' block; dropped"
            )
            return None
        parts.append(a.reshape(a.shape[0], -1))
    widths = {p.shape[1] for p in parts}
    if len(widths) != 1:
        notes.append(
            f"cell_data '{name}' blocks disagree on component count {sorted(widths)}; dropped"
        )
        return None
    return np.concatenate(parts, axis=0) if parts else np.empty((0, 1))


def _looks_like_color(name: str, values: np.ndarray) -> bool:
    """Whether a 3-wide array should render as RGB rather than as arrows.

    Gated on the *name*; the value range only confirms. See ``_COLOR_NAME_RE``.
    """
    if not _COLOR_NAME_RE.search(name):
        return False
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return False
    lo, hi = float(finite.min()), float(finite.max())
    if values.dtype == np.uint8:
        return 0 <= lo and hi <= 255
    return -1e-9 <= lo and hi <= 1.0 + 1e-9


def _finite_range(values: np.ndarray):
    """``(min, max)`` over finite entries, or ``None`` if all entries are finite.

    A single NaN makes polyscope's autoscale span the whole colormap over
    nothing, so an explicit range is needed exactly when one is present.
    """
    finite = np.isfinite(values)
    if finite.all():
        return None
    if not finite.any():
        return None
    vals = values[finite]
    return (float(vals.min()), float(vals.max()))


def _quantities_from_array(name, raw, defined_on, source, notes) -> list:
    """Map one data array onto zero or more polyscope quantities."""
    arr = np.asarray(raw)
    if arr.dtype == bool:
        arr = arr.astype(np.float64)
    if not np.issubdtype(arr.dtype, np.number) or np.issubdtype(
        arr.dtype, np.complexfloating
    ):
        notes.append(f"{source} '{name}' has non-numeric dtype {arr.dtype}; dropped")
        return []

    was_uint8 = arr.dtype == np.uint8
    flat = arr.reshape(arr.shape[0], -1) if arr.ndim > 1 else arr.reshape(-1, 1)
    width = flat.shape[1]
    vals = np.ascontiguousarray(flat, dtype=np.float64)

    def scalar(qname, column):
        col = np.ascontiguousarray(column, dtype=np.float64)
        rng = _finite_range(col)
        if rng is None and not np.isfinite(col).any() and col.size:
            notes.append(f"{source} '{qname}' has no finite values; dropped")
            return None
        return Quantity(qname, "scalar", col, defined_on, rng, source)

    def vector(qname, cols):
        v = np.ascontiguousarray(cols, dtype=np.float64)
        return Quantity(qname, "vector", v, defined_on, None, source)

    out: list = []
    if width == 1:
        q = scalar(name, vals[:, 0])
        return [q] if q else []
    if width == 2:
        notes.append(f"{source} '{name}' is 2-wide; zero-padded to a 3D vector")
        padded = np.zeros((vals.shape[0], 3), dtype=np.float64)
        padded[:, :2] = vals
        return [vector(name, padded)]
    if width == 3:
        if _looks_like_color(name, arr if was_uint8 else vals):
            rgb = vals / 255.0 if was_uint8 else np.clip(vals, 0.0, 1.0)
            return [
                Quantity(
                    name, "color", np.ascontiguousarray(rgb), defined_on, None, source
                )
            ]
        out.append(vector(name, vals))
        q = scalar(f"{name}:magnitude", np.linalg.norm(vals, axis=1))
        if q:
            out.append(q)
        return out
    if width in (6, 9):
        notes.append(
            f"{source} '{name}' is a {width}-component tensor; polyscope has no "
            "tensor quantity, so components and the Frobenius norm are shown"
        )
    if width > _MAX_COMPONENT_QUANTITIES:
        notes.append(
            f"{source} '{name}' has {width} components; only its magnitude is shown"
        )
        q = scalar(f"{name}:magnitude", np.linalg.norm(vals, axis=1))
        return [q] if q else []
    for c in range(width):
        q = scalar(f"{name}:{c}", vals[:, c])
        if q:
            out.append(q)
    if width in (6, 9):
        q = scalar(f"{name}:norm", np.linalg.norm(vals, axis=1))
        if q:
            out.append(q)
    return out


def _collect_quantities(
    point_data, cell_data_flat, cell_defined_on, notes, point_defined_on="vertices"
) -> list:
    """Build every quantity, resolving point/cell name collisions.

    Both live in one polyscope name space per structure, so a name present on
    both sides would silently shadow. The point one keeps the plain name.

    ``*_defined_on`` are polyscope's own per-structure vocabularies, which are
    not uniform: surface meshes say ``vertices``/``faces``, volume meshes
    ``vertices``/``cells``, curve networks ``nodes``/``edges``, and point
    clouds take no such argument at all (``""`` here — see
    :func:`_add_quantity`).
    """
    quantities: list = []
    for name in sorted(point_data):
        quantities.extend(
            _quantities_from_array(
                name, point_data[name], point_defined_on, "point_data", notes
            )
        )
    point_names = {q.name for q in quantities}
    for name in sorted(cell_data_flat):
        qname = name
        if name in point_names:
            qname = f"{name} (cells)"
            notes.append(
                f"'{name}' exists as both point and cell data; the cell one is "
                f"shown as '{qname}'"
            )
        quantities.extend(
            _quantities_from_array(
                qname, cell_data_flat[name], cell_defined_on, "cell_data", notes
            )
        )
    return quantities


def _gathered_cell_data(mesh, cell_source, notes) -> dict:
    """Every usable ``cell_data`` array, gathered onto the emitted primitives."""
    out = {}
    for name in mesh.cell_data:
        flat = _flatten_cell_data(mesh, name, notes)
        if flat is None:
            continue
        out[name] = flat[cell_source] if cell_source.size else flat[:0]
    return out


# --- routing -------------------------------------------------------------- #


def _max_dim(mesh) -> int:
    return max((blk.dim for blk in mesh.cells), default=0)


def _payload_points(mesh, notes) -> PolyscopePayload:
    return PolyscopePayload(
        kind="points",
        vertices=_vertices_of(mesh, notes),
        # A point cloud has nowhere else data could live, so polyscope's
        # quantity methods take no `defined_on` at all.
        quantities=_collect_quantities(mesh.point_data, {}, "", notes, ""),
        notes=notes,
    )


def _payload_curve(mesh, notes) -> PolyscopePayload:
    blocks, cell_source = _dimension_rows(mesh, 1)
    conn = _stack_conn(blocks)
    if conn.size and conn.shape[1] > 2:
        notes.append("higher-order line cells are drawn as their end segments")
        conn = np.ascontiguousarray(conn[:, :2])
    return PolyscopePayload(
        kind="curve",
        vertices=_vertices_of(mesh, notes),
        edges=conn,
        quantities=_collect_quantities(
            mesh.point_data,
            _gathered_cell_data(mesh, cell_source, notes),
            "edges",
            notes,
            # A curve network calls its points "nodes", not "vertices".
            "nodes",
        ),
        notes=notes,
    )


def _payload_surface(mesh, notes) -> PolyscopePayload:
    from ._convert_cells import _LINEAR_BASE, convert_cells

    higher = sorted(
        {blk.type for blk in mesh.cells if blk.dim == 2 and blk.type in _LINEAR_BASE}
    )
    if higher:
        # Linearize rather than truncating columns by hand: it preserves
        # cell_data 1:1, prunes the now-orphaned mid-side nodes, and remaps
        # point_data — all of which would otherwise have to be redone here.
        notes.append(
            f"higher-order surface cells ({', '.join(higher)}) are drawn as their "
            "linear bases"
        )
        mesh = convert_cells(mesh, "linearize")

    blocks, cell_source = _dimension_rows(mesh, 2)
    return PolyscopePayload(
        kind="surface",
        vertices=_vertices_of(mesh, notes),
        faces=_faces_list(blocks),
        quantities=_collect_quantities(
            mesh.point_data,
            _gathered_cell_data(mesh, cell_source, notes),
            "faces",
            notes,
        ),
        notes=notes,
    )


def _payload_surface_of_volume(mesh, notes) -> PolyscopePayload:
    """The boundary of a volume mesh, with cell data sampled from owning cells."""
    from ._convert_cells import _LINEAR_BASE, convert_cells
    from ._surface import extract_surface

    surf = extract_surface(mesh, record_parent_ids=True)
    if not surf.cells:
        raise ValueError("meshio++: view: the mesh has no extractable boundary")

    # Linearize BEFORE consuming surface:parent_cell — linearize carries it
    # through as ordinary cell data, so the two compose in this order only.
    if any(blk.type in _LINEAR_BASE for blk in surf.cells):
        notes.append("the extracted boundary is drawn with linearized cells")
        surf = convert_cells(surf, "linearize")

    parent = _flatten_cell_data(surf, "surface:parent_cell", notes)
    if parent is None:
        raise ValueError(
            "meshio++: view: could not recover boundary-face provenance; "
            "cannot map cell data onto the surface"
        )
    face_source = parent.reshape(-1).astype(np.int64)

    blocks, _ = _dimension_rows(surf, 2)
    cell_flat = _gathered_cell_data(mesh, face_source, notes)
    if mesh.cell_data:
        notes.append(
            "cell data is sampled from each boundary face's owning volume cell"
        )
    return PolyscopePayload(
        kind="surface",
        vertices=_vertices_of(surf, notes),
        faces=_faces_list(blocks),
        quantities=_collect_quantities(surf.point_data, cell_flat, "faces", notes),
        notes=notes,
    )


def _payload_volume(mesh, notes) -> PolyscopePayload:
    from ._convert_cells import convert_cells

    types = {blk.type for blk in mesh.cells if blk.dim == 3}
    bad = sorted(t for t in types if t.startswith("polyhedron"))
    if bad:
        raise ValueError(
            f"meshio++: view: cannot render polyhedron cells ({', '.join(bad)}) as a "
            "volume mesh; pass kind='surface' to draw the boundary instead"
        )

    if not types <= set(_POLYSCOPE_VOLUME_TYPES):
        # polyscope holds only tets and hexes, so wedges/pyramids/higher-order
        # must be decomposed. Any hexes present are collateral damage.
        other = sorted(types - set(_POLYSCOPE_VOLUME_TYPES))
        note = (
            f"volume cells ({', '.join(other)}) are not renderable directly and were "
            "split into tetrahedra"
        )
        if "hexahedron" in types:
            note += " — hexahedra in the same mesh were split too"
        notes.append(note)
        mesh = convert_cells(mesh, "simplexify")
        types = {blk.type for blk in mesh.cells if blk.dim == 3}

    blocks, cell_source = _dimension_rows(mesh, 3)
    quantities = _collect_quantities(
        mesh.point_data, _gathered_cell_data(mesh, cell_source, notes), "cells", notes
    )
    common = dict(
        kind="volume",
        vertices=_vertices_of(mesh, notes),
        quantities=quantities,
        notes=notes,
    )
    if types == {"tetra"}:
        return PolyscopePayload(tets=_stack_conn(blocks), **common)
    if types == {"hexahedron"}:
        return PolyscopePayload(hexes=_stack_conn(blocks), **common)
    return PolyscopePayload(mixed_cells=_mixed_cells(blocks), **common)


def _to_polyscope_payload(mesh, kind: str = "auto") -> PolyscopePayload:
    """Map a mesh onto the arrays a polyscope structure needs.

    Pure: imports no renderer, touches no display, and does not modify ``mesh``.

    Parameters
    ----------
    mesh :
        the mesh to map.
    kind :
        ``"auto"`` (default) picks by the mesh's highest cell dimension;
        ``"surface"`` draws the boundary of a volume mesh; ``"volume"`` requires
        3D cells; ``"points"`` draws the point cloud alone.

    Returns
    -------
    PolyscopePayload
        vertices, one connectivity array, the mapped quantities, and a
        ``notes`` list recording every lossy step taken.

    Raises
    ------
    ValueError
        on an unknown ``kind``, ``kind="volume"`` without 3D cells, or a
        polyhedron block under ``kind="volume"``.
    """
    if kind not in ("auto", "surface", "volume", "points", "curve"):
        raise ValueError(
            f"meshio++: view: unknown kind '{kind}' "
            "(expected auto, surface, volume, curve or points)"
        )
    notes: list = []
    dim = _max_dim(mesh)

    if kind == "points" or (kind == "auto" and dim == 0):
        return _payload_points(mesh, notes)
    if kind == "curve" or (kind == "auto" and dim == 1):
        return _payload_curve(mesh, notes)
    if kind == "volume":
        if dim != 3:
            raise ValueError(
                "meshio++: view: kind='volume' needs 3D cells, but the mesh has "
                f"none (highest cell dimension is {dim}); use kind='surface'"
            )
        return _payload_volume(mesh, notes)
    if kind == "surface" and dim == 3:
        return _payload_surface_of_volume(mesh, notes)
    if kind == "auto" and dim == 3:
        return _payload_volume(mesh, notes)
    return _payload_surface(mesh, notes)


# --- polyscope backend ---------------------------------------------------- #


def has_viewer() -> bool:
    """Whether the polyscope desktop backend is installed.

    Importing polyscope does not open a window or need a display, so this is
    safe to call anywhere.
    """
    try:
        import polyscope  # noqa: F401
    except Exception:
        return False
    return True


def _import_polyscope():
    try:
        import polyscope as ps
    except ImportError as e:
        raise ImportError(
            "meshio++: view: the polyscope backend is not installed; "
            "install it with `pip install meshioplusplus[viewer]`"
        ) from e
    return ps


def _register(ps, payload: PolyscopePayload, name: str, color_by: Optional[str] = None):
    """Register one payload as a polyscope structure and attach its quantities."""
    if payload.kind == "points":
        struct = ps.register_point_cloud(name, payload.vertices)
    elif payload.kind == "curve":
        struct = ps.register_curve_network(name, payload.vertices, payload.edges)
    elif payload.kind == "surface":
        struct = ps.register_surface_mesh(name, payload.vertices, payload.faces)
    else:
        kwargs: dict = {}
        if payload.tets is not None:
            kwargs["tets"] = payload.tets
        if payload.hexes is not None:
            kwargs["hexes"] = payload.hexes
        if payload.mixed_cells is not None:
            kwargs["mixed_cells"] = payload.mixed_cells
        struct = ps.register_volume_mesh(name, payload.vertices, **kwargs)

    for q in payload.quantities:
        _add_quantity(struct, q, enabled=q.name == color_by)
    return struct


def _add_quantity(struct, q: Quantity, enabled: bool):
    """Attach one quantity, dispatching on its kind.

    An empty ``defined_on`` means the structure has only one place data can
    live (a point cloud), where polyscope's methods take no such argument and
    passing one is a ``TypeError``.
    """
    kwargs = {"enabled": enabled}
    if q.defined_on:
        kwargs["defined_on"] = q.defined_on
    if q.kind == "scalar":
        if q.vminmax is not None:
            kwargs["vminmax"] = q.vminmax
        struct.add_scalar_quantity(q.name, q.values, **kwargs)
    elif q.kind == "vector":
        struct.add_vector_quantity(q.name, q.values, **kwargs)
    else:
        struct.add_color_quantity(q.name, q.values, **kwargs)


def _resolve_color_by(
    payload: PolyscopePayload, color_by: Optional[str]
) -> Optional[str]:
    """Validate ``color_by`` against the mapped quantities, or raise by name."""
    if color_by is None:
        return None
    names = [q.name for q in payload.quantities]
    if color_by not in names:
        available = ", ".join(names) if names else "none"
        raise ValueError(
            f"meshio++: view: no data array named '{color_by}' (available: {available})"
        )
    return color_by


def _emit_notes(payload: PolyscopePayload) -> None:
    for note in payload.notes:
        warn(f"view: {note}")


#: polyscope backends that need no window server, most preferred first. EGL is
#: the real headless path; the mock backend renders nothing but lets the whole
#: registration path run, which is enough to keep a CI job honest.
_HEADLESS_BACKENDS = ("openGL3_egl", "openGL_mock")


def _init_polyscope(ps, headless: bool) -> None:
    """Initialize polyscope once per process, headless if asked.

    polyscope refuses to re-initialize with a *different* backend, so a process
    that took a headless screenshot and then opened a window would otherwise
    raise. The first initialization wins; a later call with a different
    ``headless`` is a no-op rather than an error, because an already-working
    context is a better outcome than a correct-but-fatal complaint.
    """
    if ps.is_initialized():
        return
    if not headless:
        ps.init()
        return
    # Not every polyscope wheel is built with the EGL backend, and when it is
    # absent `set_allow_headless_backends` falls through to GLFW and fails
    # without a display. Ask for each headless backend by name instead, so the
    # failure says which ones were tried.
    errors = []
    for backend in _HEADLESS_BACKENDS:
        try:
            ps.init(backend)
            return
        except Exception as e:
            errors.append(f"{backend}: {e}")
    raise RuntimeError(
        "meshio++: screenshot: no headless rendering backend is available "
        "(tried " + "; ".join(errors) + "). Install EGL or run under xvfb."
    )


def _show_polyscope(
    payload, name, color_by, show, headless=False, window_size=None, **ps_options
):
    """Register ``payload`` and hand control to polyscope."""
    ps = _import_polyscope()
    _emit_notes(payload)
    color_by = _resolve_color_by(payload, color_by)

    _init_polyscope(ps, headless)
    # Not routed through the generic setter loop below: set_window_size takes
    # two positional arguments, not a pair.
    if window_size is not None:
        ps.set_window_size(int(window_size[0]), int(window_size[1]))
    for key, value in ps_options.items():
        setter = getattr(ps, f"set_{key}", None)
        if setter is None:
            raise TypeError(f"meshio++: view: unknown polyscope option '{key}'")
        setter(value)

    struct = _register(ps, payload, name, color_by)
    if show:
        ps.show()
    return struct


def view(
    mesh,
    kind: str = "auto",
    backend: str = "auto",
    color_by: Optional[str] = None,
    name: str = "mesh",
    show: bool = True,
    **backend_options,
):
    """Show a mesh in an interactive viewer.

    Parameters
    ----------
    mesh :
        the mesh to show (unmodified).
    kind :
        ``"auto"`` (default) picks by the mesh's highest cell dimension;
        ``"surface"`` draws the boundary of a volume mesh, which is much
        lighter for a large solid; ``"volume"``, ``"curve"`` and ``"points"``
        force a representation.
    backend :
        ``"auto"`` (default) uses polyscope when it is installed and a display
        is available, else the browser; ``"polyscope"`` opens a native window;
        ``"browser"`` renders with vtk.js, inline in a notebook or in the
        default browser.
    color_by :
        name of a mapped quantity to enable on load. Note that a vector array
        ``v`` also yields ``v:magnitude``, and a tensor its components — the
        error message on a miss lists every available name.
    name :
        structure name shown in the viewer's UI.
    show :
        when false, register the mesh but do not block on the UI. Useful for
        building up several structures, and for tests.
    **backend_options :
        backend-specific extras. For polyscope, forwarded to its
        ``set_<option>`` functions (``ground_plane_mode="none"``,
        ``up_dir="z_up"``). For the browser, ``quality=True`` bakes per-cell
        quality metrics into the page so it can colour by them offline.

    Returns
    -------
    The registered polyscope structure, or for the browser backend the path to
    the generated page (``None`` when displayed inline in a notebook).

    Raises
    ------
    ImportError
        when ``backend="polyscope"`` and polyscope is not installed. Install it
        with ``pip install meshioplusplus[viewer]``.
    ValueError
        on an unknown ``backend``/``kind``, or a ``color_by`` naming no array.
    """
    if backend not in ("auto", "polyscope", "browser"):
        raise ValueError(
            f"meshio++: view: unknown backend '{backend}' "
            "(expected auto, polyscope or browser)"
        )
    if backend == "auto":
        backend = "polyscope" if (has_viewer() and _has_display()) else "browser"

    if backend == "browser":
        from ._viewer_browser import view_browser

        return view_browser(
            mesh, kind=kind, color_by=color_by, title=name, **backend_options
        )

    payload = _to_polyscope_payload(mesh, kind)
    return _show_polyscope(payload, name, color_by, show, **backend_options)


def screenshot(
    mesh,
    path,
    kind: str = "auto",
    color_by: Optional[str] = None,
    name: str = "mesh",
    size=(1280, 960),
    transparent: bool = False,
    camera=None,
    **ps_options,
):
    """Render a mesh to a PNG without opening a window.

    Always uses the polyscope backend, headless — which is what makes it usable
    from CI and from a docs build. A vtk.js/Playwright equivalent would be a
    much heavier way to get the same picture.

    Parameters
    ----------
    mesh :
        the mesh to render (unmodified).
    path :
        output PNG path.
    kind, color_by, name, **ps_options :
        as for :func:`view`.
    size :
        ``(width, height)`` in pixels.
    transparent :
        transparent rather than opaque background.
    camera :
        ``(position, target)``, each a 3-vector, for an explicit viewpoint.
        The default home view looks straight down an axis, which renders a
        box-shaped mesh as a flat rectangle — pass a three-quarter viewpoint
        for anything meant to look like a solid.

    Raises
    ------
    ImportError
        when polyscope is not installed.
    RuntimeError
        when no headless rendering backend is available.
    """
    ps = _import_polyscope()
    payload = _to_polyscope_payload(mesh, kind)
    _show_polyscope(
        payload,
        name,
        color_by,
        show=False,
        headless=True,
        window_size=size,
        **ps_options,
    )
    ps.set_ground_plane_mode("none")
    if camera is None:
        ps.reset_camera_to_home_view()
    else:
        position, target = camera
        ps.look_at(tuple(position), tuple(target))
    ps.screenshot(str(path), transparent_bg=transparent)
    ps.remove_all_structures()
    return str(path)


def _has_display() -> bool:
    """Whether a native window can plausibly be opened.

    Deliberately conservative: on Linux the absence of both ``DISPLAY`` and
    ``WAYLAND_DISPLAY`` is a reliable "no", while macOS and Windows always
    have one. Getting this wrong only changes which backend ``"auto"`` picks.
    """
    import os
    import sys

    if sys.platform.startswith(("win", "darwin")):
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
