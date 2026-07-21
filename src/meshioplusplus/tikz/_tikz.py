from __future__ import annotations

from typing import Union

import numpy as np

from .._colormap import colormap_lookup
from .._facecolor import (
    ColorSpec,
    faces_flat,
    faces_from_projection,
    resolve_face_colors,
)
from .._files import open_file
from .._projection import ISO_AZIMUTH, ISO_ELEVATION, project_surface
from .._skin import _extract_skin_py, _has_skinnable_cells

# Number of gradient swatches the colorbar is built from (twin of
# kTikzColorbarSegments).
_COLORBAR_SEGMENTS = 32


def _color_rgb(rgb) -> str:
    """TikZ's inline colour vocabulary.

    The braces are required: without them the commas inside would split the
    surrounding key=value option list.
    """
    r, g, b = rgb
    return f"{{rgb,255:red,{r};green,{g};blue,{b}}}"


def _fill_style(fill_color: str, draw: str, line_width) -> str:
    """The per-path option list for a filled face.

    Factored out so the constant-fill and per-face-colour paths build it
    identically; with colouring off this reproduces the previous string
    exactly (twin of ``tikz_fill_style``).
    """
    opts = [f"fill={fill_color}", f"draw={draw}"]
    if line_width is not None:
        opts.append(f"line width={line_width}")
    return ", ".join(opts)


def _append_colorbar(lines, colors, min_x, max_x, min_y, max_y, float_fmt):
    """A vertical gradient bar to the right of the drawing, plus min/max labels.

    Twin of ``tikz_append_colorbar``; emitted only when asked for, so the
    default output is untouched.
    """
    num = f"{{:{float_fmt}}}"
    w = max_x - min_x
    h = max_y - min_y
    x0 = max_x + w * 0.08
    x1 = x0 + w * 0.06
    for i in range(_COLORBAR_SEGMENTS):
        y0 = min_y + h * (float(i) / float(_COLORBAR_SEGMENTS))
        y1 = min_y + h * (float(i + 1) / float(_COLORBAR_SEGMENTS))
        c = colormap_lookup(colors.table, float(i) / float(_COLORBAR_SEGMENTS - 1))
        lines.append(
            f"  \\fill[{_color_rgb(c)}] ({num.format(x0)},{num.format(y0)}) "
            f"rectangle ({num.format(x1)},{num.format(y1)});"
        )
    lines.append(
        f"  \\node[anchor=west, font=\\tiny] at ({num.format(x1)},"
        f"{num.format(min_y)}) {{{num.format(colors.vmin)}}};"
    )
    lines.append(
        f"  \\node[anchor=west, font=\\tiny] at ({num.format(x1)},"
        f"{num.format(max_y)}) {{{num.format(colors.vmax)}}};"
    )


def write(
    filename,
    mesh,
    float_fmt: str = ".6f",
    # Emit a full, directly-compilable LaTeX document by default. Set to False to
    # emit only the bare `tikzpicture` environment for \input into a larger document.
    standalone: bool = True,
    # TikZ line width for the edges, e.g. "0.4pt". None -> use TikZ's default.
    line_width: Union[str, None] = None,
    # xcolor spec for the filled faces (triangles/quads) and the edge stroke.
    fill: str = "gray!30",
    draw: str = "black",
    # Optional \begin{tikzpicture}[scale=...]. None -> omit (coordinates verbatim).
    scale: Union[int, float, None] = None,
    # Orthographic camera for genuinely 3D input (ignored for flat meshes);
    # the default is the classic CAD isometric view.
    azimuth: float = ISO_AZIMUTH,
    elevation: float = ISO_ELEVATION,
    roll: float = 0.0,
    # Data-driven colouring. `color_by` unset keeps the flat `fill` and leaves
    # the output byte-identical to previous releases.
    color_by: Union[str, None] = None,
    component: Union[int, None] = None,
    cmap: str = "viridis",
    vmin: Union[float, None] = None,
    vmax: Union[float, None] = None,
    nan_color: str = "gray",
    colorbar: bool = False,
):
    spec = ColorSpec(color_by, component, cmap, vmin, vmax, colorbar)

    if mesh.points.shape[1] == 3 and not np.allclose(
        mesh.points[:, 2], 0.0, rtol=0.0, atol=1.0e-14
    ):
        # Genuinely non-flat 3D mesh: render the projected skin instead.
        return _write_projected(
            filename,
            mesh,
            float_fmt,
            standalone,
            line_width,
            fill,
            draw,
            scale,
            azimuth,
            elevation,
            roll,
            spec,
            nan_color,
        )

    # TikZ/PGF uses the math convention (y grows upward), so unlike the SVG
    # writer there is no y-flip: the first two columns map straight to (x, y).
    pts = mesh.points[:, :2]

    coord = f"({{:{float_fmt}}},{{:{float_fmt}}})"

    # Per-path style options.
    fill_style = _fill_style(fill, draw, line_width)
    line_opts = [f"draw={draw}"]
    if line_width is not None:
        line_opts.append(f"line width={line_width}")
    line_style = ", ".join(line_opts)

    # On the flat path the drawn mesh IS the source mesh, so a face's cell
    # index needs no provenance indirection.
    colors = resolve_face_colors(spec, mesh, mesh, faces_flat(mesh))

    lines: list[str] = []
    face_index = 0
    for cell_block in mesh.cells:
        if cell_block.type not in ["line", "triangle", "quad"]:
            continue

        for cell in cell_block.data:
            path = " -- ".join(coord.format(x, y) for x, y in pts[cell])
            f = face_index
            face_index += 1
            if cell_block.type == "line":
                lines.append(f"  \\draw[{line_style}] {path};")
            elif colors.active:
                c = colors.color(f)
                fill_color = _color_rgb(c) if c is not None else nan_color
                style = _fill_style(fill_color, draw, line_width)
                lines.append(f"  \\draw[{style}] {path} -- cycle;")
            else:
                # triangle / quad: closed, filled face
                lines.append(f"  \\draw[{fill_style}] {path} -- cycle;")

    if colors.active and colorbar and len(mesh.points) > 0:
        _append_colorbar(
            lines,
            colors,
            float(np.min(pts[:, 0])),
            float(np.max(pts[:, 0])),
            float(np.min(pts[:, 1])),
            float(np.max(pts[:, 1])),
            float_fmt,
        )

    pic_opts = []
    if scale is not None:
        pic_opts.append(f"scale={scale}")
    if line_width is not None:
        pic_opts.append(f"line width={line_width}")
    pic_opt_str = f"[{', '.join(pic_opts)}]" if pic_opts else ""

    body = [f"\\begin{{tikzpicture}}{pic_opt_str}", *lines, "\\end{tikzpicture}"]

    if standalone:
        out = [
            "\\documentclass{standalone}",
            "\\usepackage{tikz}",
            "\\begin{document}",
            *body,
            "\\end{document}",
        ]
    else:
        out = body

    with open_file(filename, "w", newline="\n") as f:
        f.write("\n".join(out) + "\n")


def _write_projected(
    filename,
    mesh,
    float_fmt,
    standalone,
    line_width,
    fill,
    draw,
    scale,
    azimuth,
    elevation,
    roll,
    spec,
    nan_color,
):
    # 3D rendering path: extract the skin of any volume cells, project with
    # the orthographic camera, and draw the faces back-to-front. Mirrors the
    # C++ core's tikz_proj_write byte for byte.
    if _has_skinnable_cells(mesh):
        # Colouring by cell data needs to know which input cell each skin facet
        # came from; only ask for it when it will actually be read.
        need_parents = spec.active and spec.color_by in mesh.cell_data
        draw_mesh = _extract_skin_py(
            mesh, linearize=True, record_parent_ids=need_parents
        )
    else:
        draw_mesh = mesh
    x, y, faces = project_surface(draw_mesh, azimuth, elevation, roll)

    colors = resolve_face_colors(spec, mesh, draw_mesh, faces_from_projection(faces))

    coord = f"({{:{float_fmt}}},{{:{float_fmt}}})"

    fill_style = _fill_style(fill, draw, line_width)
    line_opts = [f"draw={draw}"]
    if line_width is not None:
        line_opts.append(f"line width={line_width}")
    line_style = ", ".join(line_opts)

    lines: list[str] = []
    for f, (nodes, is_line, _) in enumerate(faces):
        path = " -- ".join(coord.format(x[p], y[p]) for p in nodes)
        if is_line:
            lines.append(f"  \\draw[{line_style}] {path};")
        elif colors.active:
            c = colors.color(f)
            fill_color = _color_rgb(c) if c is not None else nan_color
            lines.append(
                f"  \\draw[{_fill_style(fill_color, draw, line_width)}] {path} -- cycle;"
            )
        else:
            lines.append(f"  \\draw[{fill_style}] {path} -- cycle;")

    if colors.active and spec.colorbar and len(x) > 0:
        _append_colorbar(
            lines,
            colors,
            float(np.min(x)),
            float(np.max(x)),
            float(np.min(y)),
            float(np.max(y)),
            float_fmt,
        )

    pic_opts = []
    if scale is not None:
        pic_opts.append(f"scale={scale}")
    if line_width is not None:
        pic_opts.append(f"line width={line_width}")
    pic_opt_str = f"[{', '.join(pic_opts)}]" if pic_opts else ""

    body = [f"\\begin{{tikzpicture}}{pic_opt_str}", *lines, "\\end{tikzpicture}"]

    if standalone:
        out = [
            "\\documentclass{standalone}",
            "\\usepackage{tikz}",
            "\\begin{document}",
            *body,
            "\\end{document}",
        ]
    else:
        out = body

    with open_file(filename, "w", newline="\n") as f:
        f.write("\n".join(out) + "\n")
