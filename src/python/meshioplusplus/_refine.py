"""Mesh refinement: subdivide cells into same-type children.

A dependency-free mesh *operation* (not a file format): it increases a mesh's
resolution while leaving the object it describes intact. The C++ core
(``_core.refine``) does the work; this module is the thin shim (try C++, fall
back to a pure-numpy reference) plus the ``point_sets``/``cell_sets`` remapping
the core cannot see across the Python boundary.

**Uniform** refinement applies a fixed template per cell type: ``line`` -> 2,
``triangle`` -> 4, ``quad`` -> 4, ``tetra`` -> 8, ``wedge`` -> 8, ``hexahedron``
-> 8. New nodes sit at the midpoints of the parent's edges, quad faces and
(hexahedron only) body, and carry the mean of that entity's corner values for
every ``point_data`` array. Mid-edge and quad-face-centre nodes are **shared**
between every cell touching the entity, so the refined mesh has no hanging
nodes.

**Selective** refinement (``cells=`` / ``region=`` / ``where=``) splits only the
chosen cells and resolves the resulting hanging nodes with a conforming closure,
so the output is still a valid mesh. Everything is driven by one set — which
edges carry a new node — and every other decision is derived from it: a quad
face gets a centre iff all four of its edges are split, a hexahedron a body node
iff all twelve are, and a cell's subdivision is the template for its own
split-edge mask, promoted to the smallest admissible one. See
:mod:`._refine_templates` and ``doc/refine.md``.

This reference implementation is a byte-for-byte twin of the C++ core, pinned by
``tests/python/test_refine.py::test_cpp_matches_python`` — selective mode
included. Its tables are the same tables, checked across the ``_core`` boundary
rather than transcribed.

Higher-order cells, pyramids, and ragged polygon/polyhedron blocks have no
same-type subdivision and raise :class:`ValueError`.

Public API:

* :func:`refine` — refine a mesh, uniformly or selectively.
"""

from __future__ import annotations

import warnings

import numpy as np

from ._mesh import Mesh
from ._refine_templates import (
    EDGES,
    FULL_MASK,
    NUM_CORNERS,
    QUAD_FACES,
    SUPPORTED_TYPES,
    closure_from_name,
    compare_from_name,
    face_edge_mask,
    promote_mask,
    template,
)
from ._regions import block_bases

__all__ = ["refine", "SUPPORTED_TYPES"]

#: The ``cell_data`` array :func:`refine` attaches for ``record_parent_ids``.
PARENT_CELL_NAME = "refine:parent_cell"

#: The Int64 ``cell_data`` array recording each cell's refinement depth. The name
#: is **reserved**: an input that already carries it is updated rather than
#: replicated, so successive passes accumulate.
LEVEL_NAME = "refine:level"

#: The Int64 ``point_data`` array the balanced closure attaches: ``1`` for a
#: hanging (constrained) node, ``0`` otherwise.
HANGING_NAME = "refine:hanging"

#: The Int64 ``(num_points, 4)`` ``point_data`` array recording the entity each
#: point was created on -- ``(-1, -1, a, b)`` for the midpoint of edge ``(a, b)``,
#: the sorted ``(p, q, r, s)`` for a quad-face centre, and ``(-1, -1, -1, -1)``
#: for a point ``refine`` did not create. It lets a later pass reuse a node an
#: earlier one already placed on an entity instead of allocating a coincident
#: duplicate. See ``operations/refine.hpp``'s ``kRefineEntityName``.
ENTITY_NAME = "refine:entity"

#: The Int64 ``cell_data`` array ``record_hierarchy`` attaches: a stable,
#: never-reused id per cell -- the persistent parent/child hierarchy
#: ``refine:level`` alone cannot give. See ``kRefineCellIdName``.
CELL_ID_NAME = "refine:cell_id"

#: The Int64 ``cell_data`` array naming, per output cell, the ``CELL_ID_NAME``
#: of the cell in the mesh passed to *this call* that it came from. See
#: ``kRefineParentIdName``.
PARENT_ID_NAME = "refine:parent_id"


# --------------------------------------------------------------------------- #
# validation                                                                   #
# --------------------------------------------------------------------------- #
def _check_type(cell_type, data):
    if data.ndim != 2:
        raise ValueError(
            f"refine: cannot refine ragged cell block '{cell_type}' "
            "(polygon/polyhedron blocks have no same-type subdivision template)"
        )
    if cell_type in SUPPORTED_TYPES:
        return
    if cell_type == "pyramid":
        raise ValueError(
            "refine: cannot refine cell type 'pyramid' into same-type children "
            "(its uniform refinement is 6 pyramids + 4 tetrahedra; simplexify "
            "the mesh first)"
        )
    raise ValueError(
        f"refine: cannot refine cell type '{cell_type}' into same-type children "
        "(higher-order cells have no same-type subdivision template; linearize "
        "the mesh first)"
    )


# --------------------------------------------------------------------------- #
# selection                                                                    #
# --------------------------------------------------------------------------- #
def _resolve_selection(
    mesh, cells, region, predicate_array, predicate_op, predicate_value
):
    """The per-global-cell red flag, mirroring the C++ resolver's checks.

    Only called when a selector was given, so "no selector" and "a selector that
    matched nothing" stay distinguishable — the latter is a legitimate request
    for an unchanged mesh, not a silent fall-through to refining everything.
    """
    num_set = (cells is not None) + bool(region) + bool(predicate_array)
    if num_set > 1:
        raise ValueError(
            "refine: set at most one cell selector -- cells, region or predicate "
            f"array -- but {num_set} were given"
        )

    bases = block_bases(mesh.cells)
    ncells = int(bases[-1])
    red = np.zeros(ncells, dtype=bool)

    if cells is not None:
        ids = np.asarray(cells, dtype=np.int64).reshape(-1)
        for g in ids:
            if g < 0 or g >= ncells:
                raise ValueError(
                    f"refine: cell index {int(g)} is out of range; the mesh has "
                    f"{ncells} cells"
                )
        red[ids] = True  # duplicates collapse here
        return red

    if region:
        found = {r.kind: r for r in getattr(mesh, "regions", []) if r.name == region}
        if "cell" in found:
            entries = np.asarray(found["cell"].entries, dtype=np.int64)
            entries = entries[(entries >= 0) & (entries < ncells)]
            red[entries] = True
            return red
        if "point" in found:
            npoints = len(mesh.points)
            in_region = np.zeros(npoints, dtype=bool)
            entries = np.asarray(found["point"].entries, dtype=np.int64)
            in_region[entries[(entries >= 0) & (entries < npoints)]] = True
            # ANY node in the region, not all: "all" would select nothing at all
            # for a point set describing a surface.
            for b, cb in enumerate(mesh.cells):
                data = np.asarray(cb.data)
                if data.ndim != 2:
                    continue  # rejected by _check_type later, by name
                base = int(bases[b])
                red[base : base + len(data)] = in_region[data].any(axis=1)
            return red
        if "side" in found:
            raise ValueError(
                f"refine: region '{region}' is a side region; a facet is not a "
                "cell, so it cannot select what to refine (use a cell or point "
                "region)"
            )
        names = ", ".join(repr(r.name) for r in getattr(mesh, "regions", []))
        raise ValueError(
            f"refine: no region named '{region}' (available: {names or 'none'})"
        )

    blk = mesh.cell_data.get(predicate_array)
    if blk is None:
        available = ", ".join(sorted(mesh.cell_data)) or "none"
        raise ValueError(
            f"refine: no cell_data array named '{predicate_array}' "
            f"(available: {available})"
        )
    if len(blk) != len(mesh.cells):
        raise ValueError(
            f"refine: cell_data '{predicate_array}' does not cover every cell block"
        )
    op = compare_from_name(predicate_op)
    for b, cb in enumerate(mesh.cells):
        value = np.asarray(blk[b])
        n = len(np.asarray(cb.data))
        if len(value) != n:
            raise ValueError(
                f"refine: cell_data '{predicate_array}' has {len(value)} rows on a "
                f"block of {n} cells"
            )
        if n and value.reshape(n, -1).shape[1] != 1:
            raise ValueError(
                f"refine: cell_data '{predicate_array}' must be scalar (one value "
                "per cell) to be used as a refinement predicate"
            )
        flat = value.reshape(-1).astype(np.float64, copy=False)
        # A non-finite value never matches. compute_quality deliberately reports
        # NaN where a metric does not apply, so a predicate over `quality:*` on a
        # mixed mesh is the headline use case; rejecting the array would break it.
        with np.errstate(invalid="ignore"):
            hit = np.isfinite(flat) & op(flat, predicate_value)
        base = int(bases[b])
        red[base : base + n] = hit
    return red


# --------------------------------------------------------------------------- #
# pure-numpy reference (rectangular blocks only; ragged is C++-core only)      #
# --------------------------------------------------------------------------- #
def _read_levels(mesh, blocks):
    """The per-cell refinement depth already on a mesh, or ``None``.

    A malformed array is ignored rather than rejected: it is this operation's own
    bookkeeping, not user input.
    """
    blk = mesh.cell_data.get(LEVEL_NAME)
    if blk is None or len(blk) != len(blocks):
        return None
    out = []
    for b, (_, data) in enumerate(blocks):
        value = np.asarray(blk[b]) if blk[b] is not None else None
        if value is None or len(value) != len(data):
            return None
        if len(data) and value.reshape(len(data), -1).shape[1] != 1:
            return None
        out.append(value.reshape(-1).astype(np.int64, copy=False))
    return np.concatenate(out) if out else np.empty(0, dtype=np.int64)


def _read_entity_keys(mesh):
    """The per-point entity keys already on a mesh, or ``None``.

    Mirrors ``refine_read_entity_keys`` in the C++ core, guard included: the
    entries name point indices, which ``refine`` never renumbers but
    ``reorder``/``clean``/``crop``/``merge`` do and ``transform``/``smooth``
    move, so every key must still reproduce its own point's coordinates by the
    same arithmetic that wrote them. A stale array is warned about and dropped
    rather than trusted, since trusting one grafts an entity's node onto another.
    """
    value = mesh.point_data.get(ENTITY_NAME)
    if value is None:
        return None
    n = len(mesh.points)
    value = np.asarray(value)
    if value.shape[:1] != (n,) or (n and value.reshape(n, -1).shape[1] != 4):
        warnings.warn(
            f"refine: ignoring {ENTITY_NAME!r}: it is not four values per point.",
            stacklevel=2,
        )
        return None
    keys = value.reshape(n, 4).astype(np.int64, copy=False)

    points = np.asarray(mesh.points)
    for p in range(n):
        key = keys[p]
        if key[3] < 0:
            continue  # the (-1, -1, -1, -1) sentinel
        corners = key[2:] if key[0] < 0 else key
        if np.any(corners < 0) or np.any(corners >= n):
            ok = False
        else:
            # _extend's arithmetic, verbatim: the same mean, cast back through
            # the points' own dtype.
            want = points[list(corners)].mean(axis=0).astype(points.dtype)
            ok = bool(np.array_equal(want, points[p]))
        if not ok:
            warnings.warn(
                f"refine: ignoring {ENTITY_NAME!r}: point {p} no longer sits on the "
                "entity it records, so the mesh was renumbered or moved since it "
                "was written.",
                stacklevel=2,
            )
            return None
    return keys


def _read_hierarchy(mesh, blocks):
    """Read the input's ``refine:cell_id``/``refine:parent_id``, if any.

    The numpy twin of ``refine_read_hierarchy`` in the C++ core, including its
    warn-and-drop policy: this is the operation's own bookkeeping, not user
    input, so a malformed or non-unique pair is dropped with a warning rather
    than rejected outright.

    Returns ``(ids, id_base)`` -- one id per global (block-major) input cell,
    and one past the largest id in use anywhere, taken over BOTH arrays
    defensively (a parent's id can outlive its own row once the cell is
    split) -- or ``(None, 0)`` when the input carries no valid pair.
    """
    ids_blk = mesh.cell_data.get(CELL_ID_NAME)
    if ids_blk is None:
        return None, 0
    parent_blk = mesh.cell_data.get(PARENT_ID_NAME)
    if len(ids_blk) != len(blocks) or (
        parent_blk is not None and len(parent_blk) != len(blocks)
    ):
        warnings.warn(
            f"refine: ignoring {CELL_ID_NAME!r}/{PARENT_ID_NAME!r}: they do not cover every "
            "cell block.",
            stacklevel=2,
        )
        return None, 0

    ids_parts = []
    parent_parts = [] if parent_blk is not None else None
    for b, (_, data) in enumerate(blocks):
        value = np.asarray(ids_blk[b]) if ids_blk[b] is not None else None
        if (
            value is None
            or len(value) != len(data)
            or (len(data) and value.reshape(len(data), -1).shape[1] != 1)
        ):
            warnings.warn(
                f"refine: ignoring {CELL_ID_NAME!r}: block {b} is not one scalar value per "
                "cell.",
                stacklevel=2,
            )
            return None, 0
        ids_parts.append(value.reshape(-1).astype(np.int64, copy=False))
        if parent_parts is not None:
            pvalue = np.asarray(parent_blk[b]) if parent_blk[b] is not None else None
            if (
                pvalue is None
                or len(pvalue) != len(data)
                or (len(data) and pvalue.reshape(len(data), -1).shape[1] != 1)
            ):
                warnings.warn(
                    f"refine: ignoring {PARENT_ID_NAME!r}: block {b} is not one scalar value "
                    "per cell.",
                    stacklevel=2,
                )
                return None, 0
            parent_parts.append(pvalue.reshape(-1).astype(np.int64, copy=False))

    ids = np.concatenate(ids_parts) if ids_parts else np.empty(0, dtype=np.int64)
    if np.any(ids < 0) or len(np.unique(ids)) != len(ids):
        warnings.warn(
            f"refine: ignoring {CELL_ID_NAME!r}/{PARENT_ID_NAME!r}: the ids are not both "
            "unique and non-negative, so the mesh was merged or the array was replicated by "
            "another operation.",
            stacklevel=2,
        )
        return None, 0
    max_id = int(ids.max()) if len(ids) else -1
    if parent_parts is not None:
        parents = (
            np.concatenate(parent_parts)
            if parent_parts
            else np.empty(0, dtype=np.int64)
        )
        if len(parents):
            max_id = max(max_id, int(parents.max()))
    return ids, max_id + 1


def _refine_once_py(
    mesh, red_seed=None, closure="redgreen", record_levels=False, record_hierarchy=False
):
    n = len(mesh.points)
    blocks = [(cb.type, np.asarray(cb.data)) for cb in mesh.cells]
    for cell_type, data in blocks:
        _check_type(cell_type, data)

    # --- phase 0: per-block plumbing
    first_slot, first_cell = [], []
    total_slots = total_cells = 0
    for cell_type, data in blocks:
        first_slot.append(total_slots)
        first_cell.append(total_cells)
        total_slots += len(data) * (len(EDGES[cell_type]) + len(QUAD_FACES[cell_type]))
        total_cells += len(data)

    # --- phases 1-2: entity keys in the same (block, cell, slot) order the C++
    # core uses, then a serial first-seen dedup. The id order is a pure function
    # of the key order, which is where determinism comes from.
    keys = []
    for cell_type, data in blocks:
        edges = EDGES[cell_type]
        faces = QUAD_FACES[cell_type]
        for row in data:
            for a, b in edges:
                x, y = int(row[a]), int(row[b])
                keys.append((-1, -1, x, y) if x < y else (-1, -1, y, x))
            for f in faces:
                keys.append(tuple(sorted(int(row[i]) for i in f)))

    entity_id = {}
    entities = []
    entity_of_slot = np.empty(len(keys), dtype=np.int64)
    for i, key in enumerate(keys):
        found = entity_id.get(key)
        if found is None:
            found = len(entities)
            entity_id[key] = found
            entities.append(key)
        entity_of_slot[i] = found

    # --- phase 3: the split-edge set and its closure
    masks = np.zeros(total_cells, dtype=np.int64)
    split = np.zeros(len(entities), dtype=bool)
    base_levels = _read_levels(mesh, blocks)
    base_entities = _read_entity_keys(mesh)
    if red_seed is not None and closure == "balanced":
        # 2:1 balance. A cell is split fully or not at all -- no transitional
        # templates -- and is drawn in only when a neighbour would otherwise end
        # up more than one level finer. Adjacency is by shared NODE, not by
        # shared edge: across a hanging interface the coarse cell spans a whole
        # edge while the fine cell has only half of it, so the two are different
        # entities and an edge-keyed rule is blind to exactly the coarse/fine
        # adjacency it exists to police.
        red = np.asarray(red_seed, dtype=bool).copy()
        levels = (
            np.zeros(total_cells, dtype=np.int64)
            if base_levels is None
            else base_levels
        )
        while True:
            node_max = np.full(len(mesh.points), np.iinfo(np.int64).min, dtype=np.int64)
            for b, (_, data) in enumerate(blocks):
                lo, hi = first_cell[b], first_cell[b] + len(data)
                post = levels[lo:hi] + red[lo:hi].astype(np.int64)
                np.maximum.at(node_max, data, post[:, None])
            changed = False
            for b, (_, data) in enumerate(blocks):
                lo = first_cell[b]
                for c in range(len(data)):
                    g = lo + c
                    if red[g]:
                        continue
                    if node_max[data[c]].max() > levels[g] + 1:
                        red[g] = True
                        changed = True
            if not changed:
                break
        # Only the fully split cells create nodes; an unrefined neighbour keeps
        # its whole shape and simply acquires hanging nodes on its edges.
        for b, (cell_type, data) in enumerate(blocks):
            nedges = len(EDGES[cell_type])
            nslots = nedges + len(QUAD_FACES[cell_type])
            for c in range(len(data)):
                g = first_cell[b] + c
                if not red[g]:
                    continue
                masks[g] = FULL_MASK[cell_type]
                s0 = first_slot[b] + c * nslots
                split[entity_of_slot[s0 : s0 + nedges]] = True
    elif red_seed is None:
        for b, (cell_type, data) in enumerate(blocks):
            masks[first_cell[b] : first_cell[b] + len(data)] = FULL_MASK[cell_type]
        split[:] = True
    else:
        for b, (cell_type, data) in enumerate(blocks):
            nedges = len(EDGES[cell_type])
            nslots = nedges + len(QUAD_FACES[cell_type])
            for c in range(len(data)):
                if not red_seed[first_cell[b] + c]:
                    continue
                s0 = first_slot[b] + c * nslots
                split[entity_of_slot[s0 : s0 + nedges]] = True
        # `promote_mask` is monotone and idempotent, so this converges to a unique
        # fixed point no matter what order cells are visited in.
        while True:
            promoted = np.zeros(total_cells, dtype=np.int64)
            for b, (cell_type, data) in enumerate(blocks):
                nedges = len(EDGES[cell_type])
                nslots = nedges + len(QUAD_FACES[cell_type])
                for c in range(len(data)):
                    s0 = first_slot[b] + c * nslots
                    mask = 0
                    for k in range(nedges):
                        if split[entity_of_slot[s0 + k]]:
                            mask |= 1 << k
                    promoted[first_cell[b] + c] = promote_mask(
                        cell_type, mask, closure == "propagate"
                    )
            changed = False
            for b, (cell_type, data) in enumerate(blocks):
                nedges = len(EDGES[cell_type])
                nslots = nedges + len(QUAD_FACES[cell_type])
                for c in range(len(data)):
                    mask = int(promoted[first_cell[b] + c])
                    if not mask:
                        continue
                    s0 = first_slot[b] + c * nslots
                    for k in range(nedges):
                        if mask & (1 << k):
                            e = entity_of_slot[s0 + k]
                            if not split[e]:
                                split[e] = True
                                changed = True
            if not changed:
                break
        masks = promoted

    # --- phase 4: which entities earn a node, and their point ids
    needs = np.array(
        [split[e] if key[0] < 0 else False for e, key in enumerate(entities)]
    )
    for b, (cell_type, data) in enumerate(blocks):
        faces = QUAD_FACES[cell_type]
        if not faces or not len(data):
            continue
        nedges = len(EDGES[cell_type])
        nslots = nedges + len(faces)
        face_masks = [face_edge_mask(cell_type, f) for f in faces]
        for c in range(len(data)):
            mask = int(masks[first_cell[b] + c])
            s0 = first_slot[b] + c * nslots
            for f, need in enumerate(face_masks):
                if mask & need == need:
                    needs[entity_of_slot[s0 + nedges + f]] = True

    # An entity whose key the input already recorded keeps the node it already
    # has, whether or not this pass splits it: reusing a split entity's node is
    # what stops the balance rule from tearing the mesh, and resolving an unsplit
    # one keeps the hanging rule a single statement. The (-1, -1, -1, -1)
    # sentinel is skipped -- it is what every original point and body centre
    # carries, so admitting it would collide them all on one key.
    persisted_node = {}
    if base_entities is not None:
        for p in range(n):
            key = tuple(int(x) for x in base_entities[p])
            if key[3] < 0:
                continue
            persisted_node.setdefault(key, p)

    node_of_entity = np.full(len(entities), -1, dtype=np.int64)
    if persisted_node:
        for e, key in enumerate(entities):
            found = persisted_node.get(key)
            if found is not None:
                node_of_entity[e] = found
    kept = [e for e in range(len(entities)) if needs[e] and node_of_entity[e] < 0]
    for i, e in enumerate(kept):
        node_of_entity[e] = n + i

    body_base = n + len(kept)
    body_of_cell = np.full(total_cells, -1, dtype=np.int64)
    n_bodies = 0
    for b, (cell_type, data) in enumerate(blocks):
        if cell_type != "hexahedron":
            continue
        for c in range(len(data)):
            if masks[first_cell[b] + c] != FULL_MASK[cell_type]:
                continue
            body_of_cell[first_cell[b] + c] = body_base + n_bodies
            n_bodies += 1
    n_out = body_base + n_bodies

    # --- phases 5 and 7: points and point_data
    def _extend(values):
        """Append the deduped-entity and body means to a per-point array.

        The input dtype is preserved, matching the C++ core (which writes the
        mean back through the destination array's own dtype).
        """
        out = np.empty((n_out,) + values.shape[1:], dtype=values.dtype)
        out[:n] = values
        for i, e in enumerate(kept):
            key = entities[e]
            corners = key[2:] if key[0] < 0 else key
            out[n + i] = values[list(corners)].mean(axis=0)
        for b, (cell_type, data) in enumerate(blocks):
            if cell_type != "hexahedron":
                continue
            nc = NUM_CORNERS[cell_type]
            for c, row in enumerate(data):
                body = body_of_cell[first_cell[b] + c]
                if body >= 0:
                    out[body] = values[[int(x) for x in row[:nc]]].mean(axis=0)
        return out

    out_points = _extend(np.asarray(mesh.points))

    # --- phase 6: child counts, then connectivity
    new_cells = []
    cell_maps = []
    parent_of_child = []
    red_child = []
    for b, (cell_type, data) in enumerate(blocks):
        edges = EDGES[cell_type]
        faces = QUAD_FACES[cell_type]
        nslots = len(edges) + len(faces)
        ncorners = NUM_CORNERS[cell_type]
        full = FULL_MASK[cell_type]
        ncells = len(data)

        counts = np.empty(ncells, dtype=np.int64)
        for c in range(ncells):
            counts[c] = len(template(cell_type, int(masks[first_cell[b] + c]))[0])
        first_child = np.zeros(ncells, dtype=np.int64)
        if ncells:
            first_child[1:] = np.cumsum(counts)[:-1]
        total_children = int(counts.sum()) if ncells else 0

        conn = np.empty((total_children, ncorners), dtype=np.int64)
        parents = np.empty(total_children, dtype=np.int64)
        reds = np.zeros(total_children, dtype=bool)
        for c in range(ncells):
            # A slot whose entity earned no node stays -1; a template never
            # references one, by construction of the node rules in phase 4.
            local = np.empty(ncorners + nslots + 1, dtype=np.int64)
            local[:ncorners] = data[c, :ncorners]
            s0 = first_slot[b] + c * nslots
            if nslots:
                local[ncorners : ncorners + nslots] = node_of_entity[
                    entity_of_slot[s0 : s0 + nslots]
                ]
            local[ncorners + nslots] = body_of_cell[first_cell[b] + c]

            mask = int(masks[first_cell[b] + c])
            children, children_alt, tie_a, tie_b = template(cell_type, mask)
            # The one choice made from GLOBAL node ids rather than the template's
            # local numbering, so that the cell across the affected face resolves
            # it the same way.
            if children_alt and local[tie_b] < local[tie_a]:
                children = children_alt
            lo = int(first_child[c])
            for k, child in enumerate(children):
                conn[lo + k] = local[list(child)]
                parents[lo + k] = c
                reds[lo + k] = mask == full and mask != 0
        new_cells.append((cell_type, conn))
        cell_maps.append(first_child)
        parent_of_child.append(parents)
        red_child.append(reds)

    # --- phase 8: the output's own entity keys, and the hanging set
    # Body centres and original points keep the sentinel: interior to one cell or
    # created by nothing, so nothing can ever look one up.
    out_entities = np.full((n_out, 4), -1, dtype=np.int64)
    if base_entities is not None:
        out_entities[:n] = base_entities
    for i, e in enumerate(kept):
        out_entities[n + i] = entities[e]

    # The rule is "a node sits on an entity of a cell that does not reference
    # it", evaluated over the EMITTED cells rather than over the input cells'
    # entities. Those differ once a pass splits more than one level: a cell the
    # balance rule draws in has children whose sub-edges the input cell never
    # had, and a node the neighbouring refinement places on one of them is
    # constrained without any input entity naming it.
    hanging = np.zeros(n_out, dtype=bool)
    any_partial = any(
        int(masks[first_cell[b] + c]) != FULL_MASK[cell_type]
        for b, (cell_type, data) in enumerate(blocks)
        for c in range(len(data))
    )
    if any_partial or base_entities is not None:
        node_at_entity = {}
        for p in range(n_out):
            key = tuple(int(x) for x in out_entities[p])
            if key[3] < 0:
                continue
            node_at_entity.setdefault(key, p)
        if node_at_entity:
            for b, (cell_type, conn) in enumerate(new_cells):
                edges = EDGES[cell_type]
                faces = QUAD_FACES[cell_type]
                for row in conn:
                    own = set(int(x) for x in row)
                    keys_here = []
                    for a, bb in edges:
                        x, y = int(row[a]), int(row[bb])
                        keys_here.append((-1, -1, x, y) if x < y else (-1, -1, y, x))
                    for f in faces:
                        keys_here.append(tuple(sorted(int(row[i]) for i in f)))
                    for key in keys_here:
                        node = node_at_entity.get(key)
                        if node is not None and node not in own:
                            hanging[node] = True

    point_data = {}
    for key, value in mesh.point_data.items():
        # Both reserved point arrays are this operation's own bookkeeping and are
        # rebuilt below. Carrying them through here would interpolate them: a
        # hanging flag averaged to 0.5 and truncated, an entity key to the mean
        # of two unrelated node ids.
        if key in (HANGING_NAME, ENTITY_NAME):
            continue
        value = np.asarray(value)
        if value.shape[:1] != (n,):
            point_data[key] = value.copy()
            continue
        point_data[key] = _extend(value)

    # --- phase 8: cell_data gathered parent -> children
    # base_levels was hoisted above: the balanced closure compares levels.
    cell_data = {}
    for key, blk in mesh.cell_data.items():
        if key in (LEVEL_NAME, CELL_ID_NAME, PARENT_ID_NAME):
            continue  # this operation's own bookkeeping, rebuilt separately
        out_blk = []
        for b, value in enumerate(blk):
            if value is None or b >= len(blocks):
                out_blk.append(None if value is None else np.asarray(value).copy())
                continue
            value = np.asarray(value)
            if value.shape[:1] != (len(blocks[b][1]),):
                out_blk.append(value.copy())
                continue
            out_blk.append(value[parent_of_child[b]])
        cell_data[key] = out_blk

    # refine:level. Red children increment; green and untouched cells inherit,
    # because a green split is a closure, not a refinement.
    if record_levels or base_levels is not None:
        levels_out = []
        for b in range(len(blocks)):
            parents = parent_of_child[b]
            base = (
                np.zeros(len(parents), dtype=np.int64)
                if base_levels is None
                else base_levels[first_cell[b] + parents]
            )
            levels_out.append(base + red_child[b].astype(np.int64))
        cell_data[LEVEL_NAME] = levels_out

    if hanging.any():
        point_data[HANGING_NAME] = hanging.astype(np.int64).reshape(-1, 1)

    # Attached whenever this pass leaves hanging nodes, and maintained thereafter
    # once present, exactly as refine:level is. A pass that leaves none and
    # inherits none writes nothing, which keeps every conforming closure's output
    # byte-identical to what it was before the array existed -- UNLESS
    # record_hierarchy, which forces it: it already records the coarse corners
    # each new fine node is the mean of, the multigrid prolongation stencil
    # redgreen/propagate would otherwise never expose.
    if hanging.any() or base_entities is not None or record_hierarchy:
        point_data[ENTITY_NAME] = out_entities

    out = Mesh(
        out_points,
        new_cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data={
            k: (v.copy() if isinstance(v, np.ndarray) else v)
            for k, v in mesh.field_data.items()
        },
    )
    out_red = np.concatenate(red_child) if red_child else np.zeros(0, dtype=bool)
    return out, np.arange(n, dtype=np.int64), cell_maps, out_red


def _refine_py(
    mesh,
    levels,
    record_parent_ids,
    red_seed=None,
    closure="redgreen",
    record_levels=False,
    record_hierarchy=False,
):
    if levels <= 0:
        out = Mesh(
            np.asarray(mesh.points).copy(),
            [(cb.type, np.asarray(cb.data).copy()) for cb in mesh.cells],
            point_data={k: np.asarray(v).copy() for k, v in mesh.point_data.items()},
            cell_data={
                k: [None if v is None else np.asarray(v).copy() for v in blk]
                for k, blk in mesh.cell_data.items()
            },
            field_data=dict(mesh.field_data),
        )
        point_map = np.arange(len(mesh.points), dtype=np.int64)
        cell_maps = [np.arange(len(cb.data), dtype=np.int64) for cb in mesh.cells]
    else:
        out, point_map, cell_maps, next_seed = _refine_once_py(
            mesh, red_seed, closure, record_levels, record_hierarchy
        )
        for _ in range(levels - 1):
            # Level k > 1 refines the children of level k - 1's red cells; green
            # and untouched cells are not re-refined.
            out, next_pm, next_cm, next_seed = _refine_once_py(
                out,
                next_seed if red_seed is not None else None,
                closure,
                record_levels,
                record_hierarchy,
            )
            point_map = next_pm[point_map]
            cell_maps = [next_cm[b][cm] for b, cm in enumerate(cell_maps)]

    if record_parent_ids:
        out.cell_data[PARENT_CELL_NAME] = _parent_ids(out, cell_maps)
    # Reads the ORIGINAL mesh, never an intermediate: cell_maps already
    # composes through every internal level, so "the mesh passed to this
    # call" is unambiguous regardless of levels -- see _attach_hierarchy.
    _attach_hierarchy(mesh, out, cell_maps, record_hierarchy)
    return out, point_map, cell_maps


def _parent_ids(out, cell_maps):
    """Per output cell, the index of its ORIGINAL ancestor within its block."""
    blocks = []
    for b, cb in enumerate(out.cells):
        ncells = len(cb.data)
        ids = np.full(ncells, -1, dtype=np.int64)
        first = np.asarray(cell_maps[b], dtype=np.int64) if b < len(cell_maps) else None
        if first is not None and len(first):
            ends = np.append(first[1:], ncells)
            for p, (lo, hi) in enumerate(zip(first, ends)):
                ids[lo:hi] = p
        blocks.append(ids)
    return blocks


def _attach_hierarchy(input_mesh, out, cell_maps, record_hierarchy):
    """Attach ``refine:cell_id``/``refine:parent_id`` from the FINAL composed
    cell maps -- the numpy twin of ``refine_attach_hierarchy``. One serial pass
    over the whole call's result, exactly parity with :func:`_parent_ids`: a
    composed run of length one is a cell no level ever split, so it keeps its
    own identity and is its own parent; every other cell in a longer run gets
    one fresh id from a single counter walking the output's own block/cell
    order, and carries the run's original (in THIS call's input) ancestor as
    its parent.
    """
    input_ids, id_base = _read_hierarchy(
        input_mesh, [(cb.type, np.asarray(cb.data)) for cb in input_mesh.cells]
    )
    maintained = input_ids is not None
    if not maintained and not record_hierarchy:
        return
    if not maintained:
        # Absent, or dropped as malformed/non-unique (warned already): start a
        # fresh id space from the mesh's implicit global block-major ids.
        total = sum(len(cb.data) for cb in input_mesh.cells)
        input_ids = np.arange(total, dtype=np.int64)
        id_base = total

    bases = block_bases(input_mesh.cells)
    next_id = id_base
    id_blocks = []
    parent_blocks = []
    for b, cb in enumerate(out.cells):
        ncells = len(cb.data)
        ids = np.full(ncells, -1, dtype=np.int64)
        parents = np.full(ncells, -1, dtype=np.int64)
        first = np.asarray(cell_maps[b], dtype=np.int64) if b < len(cell_maps) else None
        if first is not None and len(first):
            ends = np.append(first[1:], ncells)
            base = int(bases[b]) if b < len(bases) else 0
            for p, (lo, hi) in enumerate(zip(first, ends)):
                original_id = int(input_ids[base + p])
                if hi - lo == 1:
                    # A run of length one is a cell no level of this call ever
                    # split (every admissible non-empty mask yields more than
                    # one child), so it keeps its own identity.
                    ids[lo] = original_id
                    parents[lo] = original_id
                else:
                    for c in range(lo, hi):
                        ids[c] = next_id
                        next_id += 1
                        parents[c] = original_id
        id_blocks.append(ids.reshape(-1, 1))
        parent_blocks.append(parents.reshape(-1, 1))
    out.cell_data[CELL_ID_NAME] = id_blocks
    out.cell_data[PARENT_ID_NAME] = parent_blocks


# --------------------------------------------------------------------------- #
# set remapping (never reaches the C++ core)                                   #
# --------------------------------------------------------------------------- #
def _remap_sets(src, out, point_map, cell_maps):
    point_sets = getattr(src, "point_sets", None)
    if point_sets:
        pm = np.asarray(point_map, dtype=np.int64)
        new_ps = {}
        for name, idx in point_sets.items():
            mapped = pm[np.asarray(idx, dtype=np.int64)]
            new_ps[name] = mapped[mapped >= 0]
        out.point_sets = new_ps

    cell_sets = getattr(src, "cell_sets", None)
    if cell_sets and cell_maps is not None:
        new_cs = {}
        for name, blocks in cell_sets.items():
            remapped = []
            for b, idx in enumerate(blocks):
                if idx is None or b >= len(cell_maps):
                    remapped.append(None)
                    continue
                first = np.asarray(cell_maps[b], dtype=np.int64)
                total = len(out.cells[b].data) if b < len(out.cells) else 0
                # The end of parent c's child range is the next parent's start.
                ends = np.append(first[1:], total) if len(first) else first
                expanded = []
                for c in np.asarray(idx, dtype=np.int64):
                    if c < 0 or c >= len(first):
                        continue
                    expanded.append(np.arange(first[c], ends[c], dtype=np.int64))
                remapped.append(
                    np.concatenate(expanded)
                    if expanded
                    else np.empty(0, dtype=np.int64)
                )
            new_cs[name] = remapped
        out.cell_sets = new_cs


# --------------------------------------------------------------------------- #
# public API                                                                   #
# --------------------------------------------------------------------------- #
def refine(
    mesh,
    levels: int = 1,
    record_parent_ids: bool = False,
    *,
    cells=None,
    region: str | None = None,
    where: str | None = None,
    closure: str = "redgreen",
    record_levels: bool = False,
    record_hierarchy: bool = False,
) -> Mesh:
    """Refine a mesh, subdividing cells into same-type children.

    With no selector every cell is refined (the uniform case). With one, only the
    selected cells are, and the hanging nodes that leaves are resolved by a
    conforming closure — at most **one** selector may be given.

    :param mesh: the mesh to refine (unchanged).
    :param levels: how many times to apply the subdivision templates. ``0`` (or
        less) returns an unchanged copy. With a selector, level ``k > 1`` refines
        the children of level ``k - 1``'s fully split cells.
    :param record_parent_ids: attach an int64 ``refine:parent_cell`` ``cell_data``
        array naming, per output cell, the **original** input cell it descends
        from within its own block.
    :param cells: global (block-major) indices of the cells to refine.
    :param region: name of a region to refine. A ``cell`` region selects its own
        cells; a ``point`` region selects every cell with **any** node in it; a
        ``side`` region is an error.
    :param where: a predicate over a scalar ``cell_data`` array, spelled
        ``"NAME OP VALUE"`` (e.g. ``"quality:scaled_jacobian < 0.3"``). ``OP`` is
        one of ``<``, ``<=``, ``>``, ``>=``, ``==``, ``!=``. A non-finite cell
        value never matches.
    :param closure: ``"redgreen"`` (default) promotes an affected cell's
        split-edge mask to the smallest admissible superset, keeping the extra
        refinement local; ``"propagate"`` promotes it straight to a full split,
        which is always conforming but converges to uniform refinement of the
        whole connected component; ``"balanced"`` does not close at all, keeping
        the hanging nodes and only enforcing 2:1 balance, which makes the output
        **non-conforming** but bounds the cost by the selection (the constrained
        nodes are reported in ``refine:hanging``).
    :param record_levels: attach the int64 ``refine:level`` ``cell_data`` array
        (0 for a cell no full split touched, +1 per full split; a transitional
        child inherits its parent's level). An input already carrying it is
        updated whatever this says — the flag only controls *creating* it.
    :param record_hierarchy: attach the int64 ``refine:cell_id``/
        ``refine:parent_id`` ``cell_data`` arrays -- the persistent parent/child
        hierarchy a multigrid caller resolves across the sequence of meshes it
        keeps ("a link between two meshes, not a tree inside one"): an
        unsplit cell keeps its id and is its own parent; a split cell's
        children each get a fresh id and carry the parent's id. An input
        already carrying ``refine:cell_id`` is updated whatever this says.
        Also forces ``refine:entity`` to be attached even when the closure
        leaves no hanging node, since it already records the coarse corners
        each new fine node is the mean of -- the multigrid prolongation
        weights, which ``redgreen``/``propagate`` would otherwise never
        expose. Interacts with :func:`meshioplusplus.split`: an unqualified
        ``split(mesh, by="tag")`` auto-picks the first sorted integer
        ``cell_data`` array, and ``"refine:cell_id"`` sorts ahead of ordinary
        material tags -- name the array explicitly on a mesh carrying it.
    :returns: the refined mesh.
    :raises ValueError: on a higher-order, ``pyramid``, or ragged cell block; on
        more than one selector; on an out-of-range cell index; on an unknown
        region or a ``side`` region used as a selector; and on a predicate array
        that is missing, does not cover every block, or is not scalar.
    """
    levels = int(levels)
    predicate_array, predicate_op, predicate_value = _parse_where(where)
    closure = closure_from_name(closure)
    selective = cells is not None or bool(region) or bool(predicate_array)

    out = None
    point_map = None
    cell_maps = None
    try:
        from . import _core

        res = _core.refine(
            mesh,
            levels,
            bool(record_parent_ids),
            cells,
            region or "",
            predicate_array or "",
            predicate_op,
            predicate_value,
            closure,
            bool(record_levels),
            bool(record_hierarchy),
        )
        out = res["mesh"]
        point_map = np.asarray(res["point_map"])
        cell_maps = [np.asarray(a) for a in res["cell_maps"]]
    except ValueError:
        # A genuine user error (a block with no same-type subdivision, two
        # selectors, an unknown region) must not fall through to the numpy path
        # and silently succeed.
        raise
    except Exception:
        out = None
    used_cpp = out is not None
    if out is None:
        red_seed = (
            _resolve_selection(
                mesh, cells, region, predicate_array, predicate_op, predicate_value
            )
            if selective
            else None
        )
        out, point_map, cell_maps = _refine_py(
            mesh,
            levels,
            record_parent_ids,
            red_seed,
            closure,
            record_levels,
            record_hierarchy,
        )

    # The C++ core carries named regions across itself (and therefore
    # `point_sets`/`cell_sets`, which are views over them), so remapping again
    # here would apply the maps twice. Only the numpy fallback needs it.
    if not used_cpp:
        _remap_sets(mesh, out, point_map, cell_maps)
    return out


def _parse_where(where):
    """Split ``"NAME OP VALUE"`` into its three parts.

    The operator is the maximal run of ``<>=!`` characters, located from the
    first of them — array names routinely contain ``:`` and ``.`` but never a
    comparison character. Deliberately not a second ``data_calc`` grammar.
    """
    if where is None:
        return "", "<", 0.0
    text = str(where)
    start = next((i for i, ch in enumerate(text) if ch in "<>=!"), -1)
    if start < 0:
        raise ValueError(
            f"refine: cannot parse predicate {where!r} "
            "(expected 'NAME OP VALUE', e.g. 'quality:scaled_jacobian < 0.3')"
        )
    end = start
    while end < len(text) and text[end] in "<>=!":
        end += 1
    name = text[:start].strip()
    op = text[start:end]
    rhs = text[end:].strip()
    if not name or not rhs:
        raise ValueError(
            f"refine: cannot parse predicate {where!r} "
            "(expected 'NAME OP VALUE', e.g. 'quality:scaled_jacobian < 0.3')"
        )
    try:
        value = float(rhs)
    except ValueError:
        raise ValueError(
            f"refine: predicate {where!r} has a non-numeric threshold {rhs!r}"
        ) from None
    return name, op, value
