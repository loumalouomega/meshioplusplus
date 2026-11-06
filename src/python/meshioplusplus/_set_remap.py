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
