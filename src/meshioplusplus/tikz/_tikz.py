from __future__ import annotations

from typing import Union

import numpy as np

from .._files import open_file
from .._projection import ISO_AZIMUTH, ISO_ELEVATION, project_surface
from .._skin import _extract_skin_py, _has_skinnable_cells


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
):
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
        )

    # TikZ/PGF uses the math convention (y grows upward), so unlike the SVG
    # writer there is no y-flip: the first two columns map straight to (x, y).
    pts = mesh.points[:, :2]

    coord = f"({{:{float_fmt}}},{{:{float_fmt}}})"

    # Per-path style options.
    fill_opts = [f"fill={fill}", f"draw={draw}"]
    line_opts = [f"draw={draw}"]
    if line_width is not None:
        fill_opts.append(f"line width={line_width}")
        line_opts.append(f"line width={line_width}")
    fill_style = ", ".join(fill_opts)
    line_style = ", ".join(line_opts)

    lines: list[str] = []
    for cell_block in mesh.cells:
        if cell_block.type not in ["line", "triangle", "quad"]:
            continue

        for cell in cell_block.data:
            path = " -- ".join(coord.format(x, y) for x, y in pts[cell])
            if cell_block.type == "line":
                lines.append(f"  \\draw[{line_style}] {path};")
            else:
                # triangle / quad: closed, filled face
                lines.append(f"  \\draw[{fill_style}] {path} -- cycle;")

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

    with open_file(filename, "w") as f:
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
):
    # 3D rendering path: extract the skin of any volume cells, project with
    # the orthographic camera, and draw the faces back-to-front. Mirrors the
    # C++ core's tikz_proj_write byte for byte.
    draw_mesh = (
        _extract_skin_py(mesh, linearize=True) if _has_skinnable_cells(mesh) else mesh
    )
    x, y, faces = project_surface(draw_mesh, azimuth, elevation, roll)

    coord = f"({{:{float_fmt}}},{{:{float_fmt}}})"

    fill_opts = [f"fill={fill}", f"draw={draw}"]
    line_opts = [f"draw={draw}"]
    if line_width is not None:
        fill_opts.append(f"line width={line_width}")
        line_opts.append(f"line width={line_width}")
    fill_style = ", ".join(fill_opts)
    line_style = ", ".join(line_opts)

    lines: list[str] = []
    for nodes, is_line in faces:
        path = " -- ".join(coord.format(x[p], y[p]) for p in nodes)
        if is_line:
            lines.append(f"  \\draw[{line_style}] {path};")
        else:
            lines.append(f"  \\draw[{fill_style}] {path} -- cycle;")

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

    with open_file(filename, "w") as f:
        f.write("\n".join(out) + "\n")
