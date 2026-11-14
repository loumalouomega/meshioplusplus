"""Bake a point or cell array into a texture image through ``uv`` (roadmap 2.4.8).

The core's ``bake_texture`` (``operations/bake_texture.hpp``) first; the numpy twin
below draws the same pixels with the same arithmetic when the core declines.
"""

from __future__ import annotations

import copy
import math

import numpy as np

from ._colormap import colormap_lookup, colormap_table
from ._png import encode


def _values(array, rows, component):
    a = np.asarray(array, dtype=np.float64)
    width = a.size // rows if rows else 0
    if width == 0 or a.size != rows * width:
        raise ValueError(
            "meshio++: bake_texture: the array does not hold one row per entry"
        )
    if component >= width:
        raise ValueError(
            f"meshio++: bake_texture: component {component} is out of range"
        )
    a = a.reshape(rows, width)
    if width == 1:
        return a[:, 0].copy()
    if component >= 0:
        return a[:, component].copy()
    out = np.zeros(rows)
    for c in range(width):
        out = out + a[:, c] * a[:, c]
    return np.sqrt(out)


def _twin(mesh, array, width, height, cmap, vmin, vmax, component):
    if not (1 <= width <= 16384 and 1 <= height <= 16384):
        raise ValueError(
            "meshio++: bake_texture: width and height must lie in [1, 16384]"
        )
    np_ = len(mesh.points)
    uv = mesh.point_data.get("uv")
    if uv is None or np.asarray(uv).size != 2 * np_:
        raise ValueError(
            "meshio++: bake_texture: the mesh has no point_data['uv'] (two per point); "
            "run parameterize first"
        )
    uv = np.asarray(uv, dtype=np.float64).reshape(np_, 2)
    table = colormap_table(cmap)
    on_points = array in mesh.point_data
    if on_points:
        values = _values(mesh.point_data[array], np_, component)
    elif array in mesh.cell_data and len(mesh.cell_data[array]) == len(mesh.cells):
        values = np.concatenate(
            [
                _values(a, len(b.data), component)
                for a, b in zip(mesh.cell_data[array], mesh.cells)
            ]
        )
    else:
        raise ValueError(
            f"meshio++: bake_texture: no point or cell data array '{array}'"
        )
    finite = values[np.isfinite(values)]
    lo, hi = (float(finite.min()), float(finite.max())) if finite.size else (0.0, 0.0)
    lo = lo if vmin is None else float(vmin)
    hi = hi if vmax is None else float(vmax)
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    W, H = float(width), float(height)

    def draw(a, b, c, cell):
        xa, ya = uv[a, 0] * W, (1.0 - uv[a, 1]) * H
        xb, yb = uv[b, 0] * W, (1.0 - uv[b, 1]) * H
        xc, yc = uv[c, 0] * W, (1.0 - uv[c, 1]) * H
        area = (xb - xa) * (yc - ya) - (xc - xa) * (yb - ya)
        if not (area != 0.0) or not math.isfinite(area):
            return
        x_lo = max(0.0, math.ceil(min(xa, xb, xc) - 0.5))
        x_hi = min(W - 1.0, math.floor(max(xa, xb, xc) - 0.5))
        y_lo = max(0.0, math.ceil(min(ya, yb, yc) - 0.5))
        y_hi = min(H - 1.0, math.floor(max(ya, yb, yc) - 0.5))
        if not (x_lo <= x_hi) or not (y_lo <= y_hi):
            return
        px = np.arange(int(x_lo), int(x_hi) + 1)
        py = np.arange(int(y_lo), int(y_hi) + 1)
        cx = px.astype(np.float64)[None, :] + 0.5
        cy = py.astype(np.float64)[:, None] + 0.5
        w0 = ((xb - cx) * (yc - cy) - (xc - cx) * (yb - cy)) / area
        w1 = ((xc - cx) * (ya - cy) - (xa - cx) * (yc - cy)) / area
        w2 = ((xa - cx) * (yb - cy) - (xb - cx) * (ya - cy)) / area
        inside = (w0 >= 0.0) & (w1 >= 0.0) & (w2 >= 0.0)
        if on_points:
            value = w0 * values[a] + w1 * values[b] + w2 * values[c]
        else:
            value = np.full(w0.shape, values[cell])
        for j, i in zip(*np.nonzero(inside)):
            v = float(value[j, i])
            pixel = rgba[py[j], px[i]]
            if not math.isfinite(v):
                pixel[:] = 0
                continue
            t = (v - lo) / (hi - lo) if hi > lo else 0.5
            pixel[:3] = colormap_lookup(table, t)
            pixel[3] = 255

    cell = 0
    for block in mesh.cells:
        face = block.type in ("triangle", "quad", "polygon")
        for row in block.data:
            if face:
                row = [int(v) for v in row]
                for k in range(1, len(row) - 1):
                    draw(row[0], row[k], row[k + 1], cell)
            cell += 1
    out = copy.copy(mesh)
    out.field_data = dict(mesh.field_data)
    out.field_data["texture:png"] = np.frombuffer(encode(rgba), np.uint8).copy()
    covered = int(np.count_nonzero(rgba[..., 3]))
    warnings = (
        ["bake_texture: no 2-D cell covers a pixel of the texture"]
        if covered == 0
        else []
    )
    return {
        "mesh": out,
        "image": rgba,
        "vmin": lo,
        "vmax": hi,
        "covered_pixels": covered,
        "warnings": warnings,
    }


def bake_texture(
    mesh,
    array,
    *,
    size=(512, 512),
    cmap="viridis",
    vmin=None,
    vmax=None,
    component=None,
    return_report=False,
):
    """Bake ``array`` (point or cell data) into an RGBA texture through ``uv``.

    Each 2-D cell is drawn in texture space (``x = u * width``, ``y = (1 - v) *
    height``, top row first); a pixel inside a cell takes the array's barycentric
    value (a point array) or the cell's value (a cell array), mapped through
    ``cmap`` over ``vmin``..``vmax`` (the finite range by default). Uncovered and
    non-finite pixels stay transparent. A vector array reduces to ``component``,
    else to its magnitude.

    :param size: ``(width, height)`` in pixels.
    :returns: ``(mesh, image)``: the mesh with ``field_data["texture:png"]`` (what
        the glTF and OBJ writers embed) and the ``(height, width, 4)`` uint8 image;
        with ``return_report`` also ``{"vmin", "vmax", "covered_pixels",
        "warnings"}``.
    """
    width, height = (int(v) for v in size)
    comp = -1 if component is None else int(component)
    try:
        from . import _core
    except ImportError:
        _core = None
    res = None
    if _core is not None and hasattr(_core, "bake_texture"):
        try:
            res = _core.bake_texture(mesh, array, width, height, cmap, vmin, vmax, comp)
        except Exception as exc:
            from ._fallback import core_op_declined

            if not core_op_declined(exc, "bake_texture"):
                raise
    if res is None:
        res = _twin(mesh, array, width, height, cmap, vmin, vmax, comp)
    out, image = res.pop("mesh"), res.pop("image")
    return (out, image, res) if return_report else (out, image)
