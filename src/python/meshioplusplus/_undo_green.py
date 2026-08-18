"""Green-element undo: restore ``refine``'s transitional (closure-only) cells
back to their original parent.

A dependency-free mesh *operation* (not a file format): the standard rule for
selective refinement is to restore a transitional ("green") cell to its
parent and re-split from scratch before a new refinement pass touches the
same region; :func:`meshioplusplus.refine` refines the transitional children
directly instead, so repeated selective passes over one region degrade
element quality without bound. :func:`undo_green` is the missing half.

**A TWO-MESH operation**, ``undo_green(coarse, fine)`` -- the "link between
two meshes, not a tree inside one" framing ``refine``'s own persistent
hierarchy (``record_hierarchy=True``) already uses, taken literally: ``coarse``
is the mesh a prior ``refine(coarse, ..., record_hierarchy=True,
record_levels=True)`` call was run on, ``fine`` is that call's output.

**Design: lookup and substitution, not reconstruction.** ``refine`` never
renumbers or prunes points, so a green parent's *exact* original connectivity
and cell_data are already sitting, byte-for-byte, in ``coarse`` at the row
``fine``'s ``refine:parent_id`` names (resolved against ``coarse``'s own
``refine:cell_id``, or its implicit global-block-major id when it carries
none). This needs no per-type subdivision-table inversion and no winding
repair or other discrete sign branch -- unlike :func:`meshioplusplus.subdivide`
and :func:`meshioplusplus.agglomerate`, this operation is pure array
bookkeeping and row copies, so **it has a full numpy twin** rather than being
C++-core only. The C++ core (``_core.undo_green``) does the work; this module
is the thin shim (try C++, fall back to the pure-numpy reference below) plus
the reference itself.

A cell's mask -- and hence its red/green status -- is uniform across every one
of its children, so classification is per SIBLING GROUP (cells sharing one
``refine:parent_id``), not per cell: a singleton group
(``refine:cell_id == refine:parent_id``) is **untouched**, kept verbatim; a
group whose ``refine:level`` is one more than its coarse parent's own level is
**red** (a genuine, wanted refinement -- passed through unchanged); a group
whose level equals its coarse parent's own level is **green** (a closure
artefact -- the whole group is replaced by ONE cell copied verbatim from
``coarse``).

Points are never pruned or renumbered -- this is what makes the substitution a
zero-translation row copy. :func:`meshioplusplus.clean` with
``remove_orphans=True`` is the documented follow-up for a caller wanting the
orphaned mid-edge nodes (left behind by a substituted green group) pruned.

The six reserved ``refine:*`` cell_data arrays are unconditionally dropped
from the output -- they describe a hierarchy relationship that is now stale
after the undo. Every other ``cell_data`` array on ``fine`` is carried:
untouched/red rows copied from ``fine`` as usual, a green group's one output
row copied from the same-named array on ``coarse`` -- if ``coarse`` lacks
that array, or its shape/dtype does not match, the WHOLE array is dropped
with a warning rather than guessing a value for the substituted row.

Named **Side** regions do not survive at all (the ``subdivide``/
``agglomerate`` precedent): a removed green child's local facet numbering has
no correspondence to the substituted parent's own facets. Point and Cell
regions do survive -- Cell regions through a genuinely non-injective remap
(several fine cells collapsing onto one output row), deduplicated the same
way :class:`meshioplusplus.Region` always canonicalizes its entries.

Two honest limitations, not gaps: it can only undo the LAST generation
relative to the specific ``coarse`` mesh passed in, and it needs the caller
to hold both meshes -- there is no single-mesh fallback. It also only
supports a single-pass (``levels=1``) hierarchy; a multi-level ``refine()``
call's deeper branches are refused by name.

Public API:

* :func:`undo_green` -- restore transitional cells to their coarse parent.
"""

from __future__ import annotations

import warnings

import numpy as np

from ._mesh import Mesh
from ._refine import (
    CELL_ID_NAME,
    ENTITY_NAME,
    HANGING_NAME,
    LEVEL_NAME,
    PARENT_CELL_NAME,
    PARENT_ID_NAME,
    _read_hierarchy,
)
from ._regions import Region, block_bases

__all__ = ["undo_green"]

_RESERVED_NAMES = frozenset(
    {
        PARENT_CELL_NAME,
        LEVEL_NAME,
        HANGING_NAME,
        ENTITY_NAME,
        CELL_ID_NAME,
        PARENT_ID_NAME,
    }
)

_PREFIX = "meshio++: undo_green: "


def _global_to_block_row(bases, global_idx):
    """The inverse of `block_bases`: (block, row) for a global cell index, or
    `(None, 0)` when out of range. A linear scan -- meshes have a handful of
    blocks -- mirroring `detail::global_to_block_row` in the C++ core.
    """
    if global_idx < 0 or global_idx >= bases[-1]:
        return None, 0
    for b in range(len(bases) - 1):
        if global_idx < bases[b + 1]:
            return b, int(global_idx - bases[b])
    return None, 0


def _read_required(mesh, blocks, name):
    """A required int64 scalar cell_data array (one value per cell, covering
    every block), flattened to global-cell order. Raises `ValueError` naming
    `name` on any mismatch -- undo_green's fine-mesh preconditions are hard
    requirements, unlike `_read_hierarchy`'s warn-and-fall-back contract.
    """
    arr_blocks = mesh.cell_data.get(name)
    if arr_blocks is None or len(arr_blocks) != len(blocks):
        raise ValueError(
            f"{_PREFIX}the fine mesh has no {name!r} cell_data covering every block; run "
            "refine(..., record_hierarchy=True, record_levels=True) first"
        )
    parts = []
    for b, (_, data) in enumerate(blocks):
        value = arr_blocks[b]
        value = None if value is None else np.asarray(value)
        if (
            value is None
            or len(value) != len(data)
            or (len(data) and value.reshape(len(data), -1).shape[1] != 1)
        ):
            raise ValueError(
                f"{_PREFIX}{name!r} block {b} is not one scalar value per cell"
            )
        parts.append(value.reshape(-1).astype(np.int64, copy=False))
    return np.concatenate(parts) if parts else np.empty(0, dtype=np.int64)


def _read_optional(mesh, blocks, name):
    """Same shape check as `_read_required`, but returns `None` on absence or
    mismatch rather than raising -- used for the coarse mesh's OPTIONAL
    refine:level, where absence means "never refined" (implicit level 0).
    """
    arr_blocks = mesh.cell_data.get(name)
    if arr_blocks is None or len(arr_blocks) != len(blocks):
        return None
    parts = []
    for b, (_, data) in enumerate(blocks):
        value = arr_blocks[b]
        value = None if value is None else np.asarray(value)
        if (
            value is None
            or len(value) != len(data)
            or (len(data) and value.reshape(len(data), -1).shape[1] != 1)
        ):
            return None
        parts.append(value.reshape(-1).astype(np.int64, copy=False))
    return np.concatenate(parts) if parts else np.empty(0, dtype=np.int64)


def _undo_green_py(coarse, fine):
    if len(coarse.points) > len(fine.points):
        raise ValueError(
            f"{_PREFIX}the coarse mesh has more points than the fine mesh, so they cannot be "
            "the coarse/fine pair of one refine() call"
        )

    fine_blocks = [(cb.type, cb.data) for cb in fine.cells]
    coarse_blocks = [(cb.type, cb.data) for cb in coarse.cells]

    fine_id = _read_required(fine, fine_blocks, CELL_ID_NAME)
    fine_parent = _read_required(fine, fine_blocks, PARENT_ID_NAME)
    fine_level = _read_required(fine, fine_blocks, LEVEL_NAME)

    coarse_bases = block_bases(coarse.cells)
    total_coarse = int(coarse_bases[-1])
    coarse_ids, _ = _read_hierarchy(coarse, coarse_blocks)
    if coarse_ids is None:
        coarse_ids = np.arange(total_coarse, dtype=np.int64)
    id_to_coarse_row = {int(i): row for row, i in enumerate(coarse_ids)}

    coarse_level = _read_optional(coarse, coarse_blocks, LEVEL_NAME)

    def coarse_level_at(row):
        return 0 if coarse_level is None else int(coarse_level[row])

    total_fine = len(fine_id)
    groups_by_parent = {}
    for g in range(total_fine):
        groups_by_parent.setdefault(int(fine_parent[g]), []).append(g)

    KEEP, ANCHOR, SUPPRESSED = 0, 1, 2
    role = np.zeros(total_fine, dtype=np.int8)
    group_id_of = np.full(total_fine, -1, dtype=np.int64)
    group_coarse_row = []
    group_size = []

    for parent_id, members in groups_by_parent.items():
        members.sort()
        if len(members) == 1:
            g = members[0]
            if int(fine_id[g]) != parent_id:
                raise ValueError(
                    f"{_PREFIX}malformed hierarchy: a singleton sibling group's "
                    f"refine:cell_id does not equal its refine:parent_id (cell {g})"
                )
            continue  # untouched: role stays KEEP

        coarse_row = id_to_coarse_row.get(parent_id)
        if coarse_row is None:
            raise ValueError(
                f"{_PREFIX}refine:parent_id {parent_id} does not resolve in the coarse mesh's "
                "id space -- these two meshes are not the coarse/fine pair of one refine() call"
            )
        clevel = coarse_level_at(coarse_row)
        flevel = int(fine_level[members[0]])
        if any(int(fine_level[g]) != flevel for g in members):
            raise ValueError(
                f"{_PREFIX}malformed hierarchy: the sibling group under refine:parent_id "
                f"{parent_id} does not agree on refine:level"
            )

        if flevel == clevel + 1:
            continue  # red: a genuine refinement, kept unchanged

        if flevel > clevel + 1:
            raise ValueError(
                f"{_PREFIX}the sibling group under refine:parent_id {parent_id} has "
                f"refine:level {flevel}, more than one deeper than its coarse parent's own "
                f"level ({clevel}); undo_green only supports a single-pass (levels=1) "
                "hierarchy"
            )
        if flevel != clevel:
            raise ValueError(
                f"{_PREFIX}the sibling group under refine:parent_id {parent_id} has "
                f"refine:level {flevel}, which is neither its coarse parent's own level "
                f"({clevel}, green) nor one more ({clevel + 1}, red)"
            )

        # green: substitute the whole group with one row from the coarse mesh
        gid = len(group_coarse_row)
        group_coarse_row.append(coarse_row)
        group_size.append(len(members))
        role[members[0]] = ANCHOR
        group_id_of[members[0]] = gid
        for g in members[1:]:
            role[g] = SUPPRESSED
            group_id_of[g] = gid

    # --- build the output mesh: fine's own block structure, unchanged types
    # and order, rows compacted/substituted -----------------------------------
    fine_bases = block_bases(fine.cells)
    nblocks = len(fine.cells)
    out_ncells = [
        int(np.count_nonzero(role[fine_bases[b] : fine_bases[b + 1]] != SUPPRESSED))
        for b in range(nblocks)
    ]

    cell_maps = [np.empty(len(cb.data), dtype=np.int64) for cb in fine.cells]
    out_conn = []
    group_output_global = [-1] * len(group_coarse_row)
    group_fine_block = [0] * len(group_coarse_row)

    out_block_base = 0
    for b, cb in enumerate(fine.cells):
        n = len(cb.data)
        conn = np.empty((out_ncells[b],) + cb.data.shape[1:], dtype=cb.data.dtype)
        cm = cell_maps[b]
        out_row = 0
        for r in range(n):
            g = fine_bases[b] + r
            if role[g] == SUPPRESSED:
                cm[r] = group_output_global[group_id_of[g]]
                continue
            if role[g] == KEEP:
                conn[out_row] = cb.data[r]
            else:  # ANCHOR
                gid = group_id_of[g]
                blk, row = _global_to_block_row(coarse_bases, group_coarse_row[gid])
                conn[out_row] = coarse.cells[blk].data[row]
                group_output_global[gid] = out_block_base + out_row
                group_fine_block[gid] = b
            cm[r] = out_block_base + out_row
            out_row += 1
        out_conn.append(conn)
        out_block_base += out_ncells[b]

    out_cells = [(cb.type, out_conn[b]) for b, cb in enumerate(fine.cells)]

    # --- point_data / field_data: unchanged, no new or pruned points --------
    # (refine:entity / refine:hanging are POINT data, not cell_data, but are
    # just as much stale hierarchy bookkeeping as the four cell_data ones)
    point_data = {
        name: np.array(arr).copy()
        for name, arr in fine.point_data.items()
        if name not in _RESERVED_NAMES
    }
    field_data = {
        name: (arr.copy() if isinstance(arr, np.ndarray) else arr)
        for name, arr in fine.field_data.items()
    }

    # --- cell_data: reserved refine:* arrays dropped; everything else kept
    # (fine's own row for KEEP, coarse's for ANCHOR) --------------------------
    cell_data = {}
    for name, arr_blocks in fine.cell_data.items():
        if name in _RESERVED_NAMES:
            continue
        if len(arr_blocks) != nblocks or any(b is None for b in arr_blocks):
            warnings.warn(
                f"{_PREFIX}cell_data {name!r} does not have one array per fine block; dropped "
                "rather than guessed at",
                stacklevel=2,
            )
            continue

        ok = True
        if group_coarse_row:
            coarse_arr_blocks = coarse.cell_data.get(name)
            if coarse_arr_blocks is None or len(coarse_arr_blocks) != len(coarse.cells):
                ok = False
            else:
                for gid in range(len(group_coarse_row)):
                    blk, _row = _global_to_block_row(
                        coarse_bases, group_coarse_row[gid]
                    )
                    if coarse_arr_blocks[blk] is None:
                        ok = False
                        break
                    fsrc = np.asarray(arr_blocks[group_fine_block[gid]])
                    csrc = np.asarray(coarse_arr_blocks[blk])
                    if fsrc.dtype != csrc.dtype or fsrc.shape[1:] != csrc.shape[1:]:
                        ok = False
                        break
        if not ok:
            warnings.warn(
                f"{_PREFIX}cell_data {name!r} cannot be honestly restored for a substituted "
                "cell (missing, incomplete or a different shape/dtype on the coarse mesh); "
                "dropped",
                stacklevel=2,
            )
            continue

        out_blocks = []
        for b, cb in enumerate(fine.cells):
            fsrc = np.asarray(arr_blocks[b])
            dst = np.empty((out_ncells[b],) + fsrc.shape[1:], dtype=fsrc.dtype)
            out_row = 0
            for r in range(len(cb.data)):
                g = fine_bases[b] + r
                if role[g] == SUPPRESSED:
                    continue
                if role[g] == KEEP:
                    dst[out_row] = fsrc[r]
                else:  # ANCHOR
                    gid = group_id_of[g]
                    blk, row = _global_to_block_row(coarse_bases, group_coarse_row[gid])
                    dst[out_row] = np.asarray(coarse.cell_data[name][blk])[row]
                out_row += 1
            out_blocks.append(dst)
        cell_data[name] = out_blocks

    out = Mesh(
        np.array(fine.points).copy(),
        out_cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data=field_data,
    )

    # --- regions: non-injective remap + dedup (Region auto-canonicalizes);
    # named Side regions do not survive at all --------------------------------
    new_regions = []
    for region in fine.regions:
        if region.kind == "side":
            continue
        if region.kind == "point":
            new_regions.append(region.copy())
            continue
        mapped = []
        for entry in region.entries:
            blk, row = _global_to_block_row(fine_bases, int(entry))
            if blk is None:
                continue
            mapped.append(int(cell_maps[blk][row]))
        new_regions.append(
            Region(
                region.name,
                "cell",
                np.array(mapped, dtype=np.int64),
                region.dim,
                region.tag,
            )
        )
    out.regions = new_regions

    report = {
        "num_groups_undone": len(group_coarse_row),
        "num_cells_removed": int(sum(sz - 1 for sz in group_size)),
    }
    return out, report


def undo_green(coarse, fine, return_report: bool = False):
    """Restore ``fine``'s transitional (green) cells back to their original
    parent, read verbatim from ``coarse``.

    :param coarse: the mesh a prior ``refine(coarse, ...)`` call was run on
        (never modified).
    :param fine: that call's output -- must carry one int64 scalar
        ``cell_data`` array per block for each of ``refine:cell_id``,
        ``refine:parent_id`` and ``refine:level`` (i.e. that call must have
        used ``record_hierarchy=True, record_levels=True``); raises
        otherwise.
    :param return_report: also return ``{num_groups_undone,
        num_cells_removed}``.
    :returns: the undone mesh, or ``(mesh, report)`` when ``return_report``
        is set.
    :raises ValueError: when ``fine`` lacks the required hierarchy arrays,
        when a ``refine:parent_id`` value cannot be resolved against
        ``coarse``'s id space (the two meshes are not the input/output pair
        of one ``refine()`` call), or when a sibling group's
        ``refine:level`` matches neither the red nor the green relationship
        to its coarse parent's own level.
    """
    out = None
    report = None
    try:
        from . import _core

        res = _core.undo_green(coarse, fine)
        out = res["mesh"]
        report = {
            "num_groups_undone": res["num_groups_undone"],
            "num_cells_removed": res["num_cells_removed"],
        }
    except (ValueError, TypeError):
        raise
    except Exception:
        out = None

    if out is None:
        out, report = _undo_green_py(coarse, fine)

    return (out, report) if return_report else out
