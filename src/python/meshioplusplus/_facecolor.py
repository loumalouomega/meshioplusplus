"""Resolve a named data array into one color per drawn face.

The twin of ``src/cpp/include/meshioplusplus/detail/face_color.hpp`` -- see that
header for the resolution rules. The SVG and TikZ writers are pinned
byte-for-byte against the C++ core (``tests/python/test_svg.py``,
``tests/python/test_tikz.py``), so every arithmetic step here mirrors the C++ one
expression for expression: sums run left to right in index order, values are
cast to float before any arithmetic, and no NumPy reduction (``np.mean``,
``np.linalg.norm``, ``np.nanmin``) is used, since their pairwise summation
rounds differently from a sequential loop.
"""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np

from ._colormap import colormap_lookup, colormap_table

# Cell types the flat 2D path draws, in the writers' own enumeration order.
_FLAT_TYPES = ("line", "triangle", "quad")


class ColorSpec(NamedTuple):
    """What the writer was asked to color by."""

    color_by: str | None = None
    component: int | None = None
    cmap: str = "viridis"
    vmin: float | None = None
    vmax: float | None = None
    colorbar: bool = False

    @property
    def active(self) -> bool:
        return bool(self.color_by)


class ColorFace(NamedTuple):
    """One drawn face, as the resolver needs to see it."""

    nodes: object  # corner node ids into the drawn mesh
    source_cell: int  # cell index in the drawn mesh, block-major


class FaceColors(NamedTuple):
    """Per-face scalars and the resolved range, ready for lookup."""

    active: bool
    values: list  # one per face; non-finite => nan_color
    vmin: float
    vmax: float
    table: bytes | None

    def color(self, face: int) -> tuple[int, int, int] | None:
        """The color for one face, or ``None`` to use the writer's nan_color."""
        if not self.active:
            return None
        v = self.values[face]
        if not math.isfinite(v):
            return None
        return colormap_lookup(self.table, color_param(v, self.vmin, self.vmax))


INACTIVE = FaceColors(False, [], 0.0, 1.0, None)


def color_param(v: float, vmin: float, vmax: float) -> float:
    """Normalize a value into ``[0, 1]`` over the range, clamping.

    A degenerate range (``vmin == vmax``, which an all-constant array produces)
    maps everything to the middle of the colormap rather than dividing by zero.
    """
    if not (vmax > vmin):
        return 0.5
    t = (v - vmin) / (vmax - vmin)
    if not (t > 0.0):
        return 0.0
    if t > 1.0:
        return 1.0
    return t


def _num_components(array, num_entries: int) -> int:
    """Components = product of the trailing dimensions."""
    if num_entries == 0:
        return 0
    return int(np.asarray(array).size // num_entries)


def _scalarize(flat, row: int, ncomp: int, component: int | None) -> float:
    """Reduce one row to a scalar: the requested component, or the magnitude.

    The magnitude sums left to right and takes one sqrt at the end -- never
    ``np.linalg.norm`` -- so it rounds exactly as the C++ core does.
    """
    if ncomp == 0:
        return float("nan")
    base = row * ncomp
    if component is not None:
        if component < 0 or component >= ncomp:
            raise ValueError(
                f"meshio++: component {component} is out of range for an array "
                f"with {ncomp} component(s)"
            )
        return float(flat[base + component])
    if ncomp == 1:
        return float(flat[base])
    total = 0.0
    for k in range(ncomp):
        v = float(flat[base + k])
        total += v * v
    return math.sqrt(total)


def _flatten_cell_data(mesh, name: str, component: int | None) -> list:
    """One scalar per global cell, block-major -- the mSourceCell order."""
    blocks = mesh.cell_data[name]
    if len(blocks) != len(mesh.cells):
        raise ValueError(
            f"meshio++: cell_data array {name!r} has {len(blocks)} block(s) but "
            f"the mesh has {len(mesh.cells)}"
        )
    out = []
    for block, values in zip(mesh.cells, blocks):
        n = len(block.data)
        arr = np.asarray(values)
        ncomp = _num_components(arr, n)
        flat = arr.reshape(-1)
        for r in range(n):
            out.append(_scalarize(flat, r, ncomp, component))
    return out


def _flatten_parent_ids(mesh) -> list:
    """Per-cell parent ids of the drawn mesh, or ``[]`` when it has none."""
    key = "surface:parent_cell"
    blocks = mesh.cell_data.get(key)
    if blocks is None or len(blocks) != len(mesh.cells):
        return []
    out = []
    for values in blocks:
        out.extend(int(v) for v in np.asarray(values).reshape(-1))
    return out


def faces_from_projection(faces) -> list:
    """The faces a projected surface draws, in emission order."""
    return [ColorFace(nodes, source) for nodes, _, source in faces]


def faces_flat(mesh) -> list:
    """The faces the flat 2D path draws, in its emission order.

    Enumerates ``line``/``triangle``/``quad`` blocks block-major then
    cell-major -- the same rule, and the same order, both flat writers loop
    over -- with the global cell counter advancing over skipped blocks too.
    """
    out = []
    base = 0
    for block in mesh.cells:
        block_base = base
        base += len(block.data)
        if block.type not in _FLAT_TYPES:
            continue
        conn = np.asarray(block.data, dtype=np.int64)
        for r in range(len(conn)):
            out.append(ColorFace(conn[r][:4], block_base + r))
    return out


def resolve_face_colors(spec: ColorSpec, source, draw, faces) -> FaceColors:
    """Resolve ``spec`` into one scalar per face plus the mapped range.

    ``source`` is the writer's input mesh -- the one the array name refers to;
    ``draw`` is the mesh actually being drawn, equal to ``source`` unless a
    boundary was extracted, in which case its ``"surface:parent_cell"`` carries
    the provenance back to ``source``'s cells.
    """
    if not spec.active:
        return INACTIVE

    if spec.vmin is not None and spec.vmax is not None and spec.vmin > spec.vmax:
        raise ValueError("meshio++: vmin must not exceed vmax")

    # Resolve the colormap first: an unknown name should fail before any work.
    table = colormap_table(spec.cmap)

    name = spec.color_by
    is_point = name in source.point_data
    if not is_point and name not in source.cell_data:
        available = [f"point:{n}" for n in sorted(source.point_data)] + [
            f"cell:{n}" for n in sorted(source.cell_data)
        ]
        detail = (
            f" (available: {', '.join(available)})"
            if available
            else " (the mesh carries no data arrays)"
        )
        raise ValueError(
            f"meshio++: no point_data or cell_data array named {name!r}{detail}"
        )

    values = []
    if is_point:
        # Point data lives on the drawn mesh: a boundary extraction carries
        # point_data through its point compaction, so the drawn mesh's node ids
        # index the already-subset array.
        if name not in draw.point_data:
            raise ValueError(
                f"meshio++: point_data array {name!r} did not survive boundary "
                "extraction"
            )
        arr = np.asarray(draw.point_data[name])
        ncomp = _num_components(arr, len(draw.points))
        flat = arr.reshape(-1)
        for face in faces:
            # Mean of the corner values, summed in node order then divided
            # once -- the same shape as the projection's depth accumulation.
            nodes = face.nodes
            total = 0.0
            count = 0
            for p in nodes:
                total += _scalarize(flat, int(p), ncomp, spec.component)
                count += 1
            values.append(float("nan") if count == 0 else total / float(count))
    else:
        cell_values = _flatten_cell_data(source, name, spec.component)
        parents = _flatten_parent_ids(draw)
        for face in faces:
            cell = face.source_cell
            if parents:
                if cell < 0 or cell >= len(parents):
                    values.append(float("nan"))
                    continue
                cell = parents[cell]
            if cell < 0 or cell >= len(cell_values):
                values.append(float("nan"))
            else:
                values.append(cell_values[cell])

    # Auto range over the DRAWN faces, not the whole array: the visible figure
    # then spans the whole colorbar. (ParaView's default is the whole array;
    # this is documented as differing.) Non-finite values are excluded, with
    # FiniteStats' seed-on-first-finite rule rather than +/-inf sentinels.
    lo = 0.0
    hi = 0.0
    seen = False
    for v in values:
        if not math.isfinite(v):
            continue
        if not seen:
            lo = v
            hi = v
            seen = True
        else:
            if v < lo:
                lo = v
            if v > hi:
                hi = v
    vmin = lo if spec.vmin is None else float(spec.vmin)
    vmax = hi if spec.vmax is None else float(spec.vmax)
    if vmin > vmax:
        raise ValueError("meshio++: vmin must not exceed vmax")
    return FaceColors(True, values, vmin, vmax, table)
