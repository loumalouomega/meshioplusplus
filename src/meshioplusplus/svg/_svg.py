from __future__ import annotations

from typing import Union
from xml.etree import ElementTree as ET

import numpy as np

from .._projection import ISO_AZIMUTH, ISO_ELEVATION, project_surface
from .._skin import _extract_skin_py, _has_skinnable_cells


def write(
    filename,
    mesh,
    float_fmt: str = ".3f",
    stroke_width: Union[str, None] = None,
    # Use a default image_width (not None). If set to None, images will come out at the
    # width of the mesh (which is okay). Some viewers (e.g., eog) have problems
    # displaying SVGs of width around 1 since they interpret it as the width in pixels.
    image_width: Union[int, float, None] = 100,
    # ParaView's default colors
    fill: str = "#c8c5bd",
    stroke: str = "#000080",
    # Orthographic camera for genuinely 3D input (ignored for flat meshes);
    # the default is the classic CAD isometric view.
    azimuth: float = ISO_AZIMUTH,
    elevation: float = ISO_ELEVATION,
    roll: float = 0.0,
):
    if mesh.points.shape[1] == 3 and not np.allclose(
        mesh.points[:, 2], 0.0, rtol=0.0, atol=1.0e-14
    ):
        # Genuinely non-flat 3D mesh: render the projected skin instead.
        return _write_projected(
            filename,
            mesh,
            float_fmt,
            stroke_width,
            image_width,
            fill,
            stroke,
            azimuth,
            elevation,
            roll,
        )

    pts = mesh.points[:, :2].copy()
    min_x = np.min(pts[:, 0]) if len(pts) > 0 else 0.0
    max_x = np.max(pts[:, 0]) if len(pts) > 0 else 0.0
    min_y = np.min(pts[:, 1]) if len(pts) > 0 else 0.0
    max_y = np.max(pts[:, 1]) if len(pts) > 0 else 0.0

    pts[:, 1] = max_y + min_y - pts[:, 1]

    width = max_x - min_x
    height = max_y - min_y

    if image_width is not None and width != 0:
        scaling_factor = image_width / width
        min_x *= scaling_factor
        min_y *= scaling_factor
        width *= scaling_factor
        height *= scaling_factor
        pts *= scaling_factor

    if stroke_width is None:
        stroke_width = str(width / 100)

    fmt = " ".join(4 * [f"{{:{float_fmt}}}"])
    svg = ET.Element(
        "svg",
        xmlns="http://www.w3.org/2000/svg",
        version="1.1",
        viewBox=fmt.format(min_x, min_y, width, height),
    )

    style = ET.SubElement(svg, "style")
    opts = [
        f"fill: {fill}",
        f"stroke: {stroke}",
        f"stroke-width: {stroke_width}",
        "stroke-linejoin:bevel",
    ]
    # Use path, not polygon, because svgo converts polygons to paths and doesn't convert
    # the style alongside. No problem if it's paths all the way.
    style.text = "path {" + "; ".join(opts) + "}"

    for cell_block in mesh.cells:
        if cell_block.type not in ["line", "triangle", "quad"]:
            continue

        if cell_block.type == "line":
            fmt = (
                f"M {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
            )
        elif cell_block.type == "triangle":
            fmt = (
                f"M {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
                + "Z"
            )
        elif cell_block.type == "quad":
            fmt = (
                f"M {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
                + f"L {{:{float_fmt}}} {{:{float_fmt}}}"
                + "Z"
            )
        for cell in cell_block.data:
            ET.SubElement(
                svg,
                "path",
                d=fmt.format(*pts[cell].flatten()),
            )

    tree = ET.ElementTree(svg)
    tree.write(filename)


def _write_projected(
    filename,
    mesh,
    float_fmt,
    stroke_width,
    image_width,
    fill,
    stroke,
    azimuth,
    elevation,
    roll,
):
    # 3D rendering path: extract the skin of any volume cells, project with
    # the orthographic camera, and paint the faces back-to-front. Mirrors the
    # C++ core's svg_proj_write (bbox -> y-flip -> scaling -> emission).
    draw_mesh = (
        _extract_skin_py(mesh, linearize=True) if _has_skinnable_cells(mesh) else mesh
    )
    x, y, faces = project_surface(draw_mesh, azimuth, elevation, roll)

    min_x = np.min(x) if len(x) > 0 else 0.0
    max_x = np.max(x) if len(x) > 0 else 0.0
    min_y = np.min(y) if len(y) > 0 else 0.0
    max_y = np.max(y) if len(y) > 0 else 0.0

    # Flip y (projected math convention y-up -> SVG screen convention y-down).
    y = max_y + min_y - y

    width = max_x - min_x
    height = max_y - min_y

    if image_width is not None and width != 0:
        scaling_factor = image_width / width
        min_x *= scaling_factor
        min_y *= scaling_factor
        width *= scaling_factor
        height *= scaling_factor
        x = x * scaling_factor
        y = y * scaling_factor

    if stroke_width is None:
        stroke_width = str(width / 100)

    fmt = " ".join(4 * [f"{{:{float_fmt}}}"])
    svg = ET.Element(
        "svg",
        xmlns="http://www.w3.org/2000/svg",
        version="1.1",
        viewBox=fmt.format(min_x, min_y, width, height),
    )

    style = ET.SubElement(svg, "style")
    opts = [
        f"fill: {fill}",
        f"stroke: {stroke}",
        f"stroke-width: {stroke_width}",
        "stroke-linejoin:bevel",
    ]
    style.text = "path {" + "; ".join(opts) + "}"

    point_fmt = f"{{:{float_fmt}}} {{:{float_fmt}}}"
    for nodes, is_line in faces:
        parts = []
        for k, p in enumerate(nodes):
            parts.append(("M " if k == 0 else "L ") + point_fmt.format(x[p], y[p]))
        d = "".join(parts)
        if not is_line:
            d += "Z"
        ET.SubElement(svg, "path", d=d)

    tree = ET.ElementTree(svg)
    tree.write(filename)
