"""Carry the shim-only ``point_sets``/``cell_sets`` through an operation.

The C++ core never sees a mesh's compatibility sets, so the Python shims of
``clean``, ``crop``, ``convert_cells`` and ``refine`` remap them onto the
result from the operation's point and cell maps.
"""

import numpy as np


def remap_sets(src, out, point_map, cell_maps, *, expand=False):
    """Remap ``src``'s ``point_sets``/``cell_sets`` onto ``out``.

    ``point_map[p]`` is old point ``p``'s new index (negative: dropped).

    With ``expand=False``, ``cell_maps[b][c]`` is old cell ``c``'s new index in
    block ``b`` (negative: dropped). With ``expand=True`` it is the index of
    parent ``c``'s *first* child, and a parent's children are contiguous, so a
    cell set entry expands to the whole ``[first[c], first[c + 1])`` range.
    """
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
                if expand:
                    remapped.append(_expand_children(out, b, cell_maps[b], idx))
                    continue
                mapped = np.asarray(cell_maps[b], dtype=np.int64)[
                    np.asarray(idx, dtype=np.int64)
                ]
                remapped.append(mapped[mapped >= 0])
            new_cs[name] = remapped
        out.cell_sets = new_cs


def _expand_children(out, b, cell_map, idx):
    first = np.asarray(cell_map, dtype=np.int64)
    total = len(out.cells[b].data) if b < len(out.cells) else 0
    # The end of parent c's child range is the next parent's start.
    ends = np.append(first[1:], total) if len(first) else first
    expanded = []
    for c in np.asarray(idx, dtype=np.int64):
        if c < 0 or c >= len(first):
            continue
        expanded.append(np.arange(first[c], ends[c], dtype=np.int64))
    return np.concatenate(expanded) if expanded else np.empty(0, dtype=np.int64)


# --------------------------------------------------------------------------- #
# pieces: split and partition                                                  #
# --------------------------------------------------------------------------- #
def point_set_is_remappable(idx, num_points):
    """Whether a ``point_sets`` entry is point indices into ``num_points`` points."""
    idx = np.asarray(idx)
    return idx.dtype.kind in "iu" and (
        len(idx) == 0 or (idx.min() >= 0 and idx.max() < num_points)
    )


def cell_set_is_remappable(cells, blocks):
    """Whether a ``cell_sets`` entry is per-block cell-index lists into ``cells``.

    Formats can stash other per-block side data in ``cell_sets`` -- gmsh's
    ``gmsh:bounding_entities`` holds entity tags (negative values included),
    not cell indices -- and remapping those as indices would crash or corrupt
    them. Stricter than ``_regions.cell_set_is_index_data``: an index past its
    block's end also fails.
    """
    for b, cb in enumerate(cells):
        if b >= len(blocks) or blocks[b] is None:
            continue
        idx = np.asarray(blocks[b]).reshape(-1)
        if len(idx) == 0:
            continue
        if idx.dtype.kind not in "iu":
            return False
        if idx.min() < 0 or idx.max() >= len(cb.data):
            return False
    return True


def remap_piece_sets(
    src, out, point_map, cell_maps, *, keep_invalid, drop_empty_blocks
):
    """Remap ``src``'s ``point_sets``/``cell_sets`` onto one piece ``out``.

    ``point_map``/``cell_maps`` map old indices to the piece's (negative:
    not in the piece). An entry that is not indices (see
    ``cell_set_is_remappable``) is copied through unchanged with
    ``keep_invalid`` -- ``partition``, whose pieces keep every block, so
    per-block side data stays aligned -- and skipped otherwise (``split``,
    which drops blocks). ``drop_empty_blocks`` leaves out the blocks no cell
    of the piece comes from, as the piece itself does.
    """
    point_sets = getattr(src, "point_sets", None)
    if point_sets:
        pm = np.asarray(point_map, dtype=np.int64)
        new_ps = {}
        for name, idx in point_sets.items():
            idx = np.asarray(idx)
            if point_set_is_remappable(idx, len(pm)):
                mapped = pm[idx.astype(np.int64)]
                new_ps[name] = mapped[mapped >= 0]
            elif keep_invalid:
                new_ps[name] = idx.copy()
        if new_ps:
            out.point_sets = new_ps
    cell_sets = getattr(src, "cell_sets", None)
    if cell_sets and cell_maps is not None:
        emit_blocks = range(len(cell_maps))
        if drop_empty_blocks:
            emit_blocks = [
                b for b in emit_blocks if np.any(np.asarray(cell_maps[b]) >= 0)
            ]
        new_cs = {}
        for name, blocks in cell_sets.items():
            if not cell_set_is_remappable(src.cells, blocks):
                if keep_invalid:
                    new_cs[name] = [
                        None if b is None else np.asarray(b).copy() for b in blocks
                    ]
                continue
            out_blocks = []
            for b in emit_blocks:
                if b < len(blocks) and blocks[b] is not None:
                    mapped = np.asarray(cell_maps[b], dtype=np.int64)[
                        np.asarray(blocks[b], dtype=np.int64).reshape(-1)
                    ]
                    out_blocks.append(mapped[mapped >= 0])
                else:
                    out_blocks.append(np.empty(0, dtype=np.int64))
            new_cs[name] = out_blocks
        if new_cs:
            out.cell_sets = new_cs
