"""Software rendering: a mesh to an RGBA image, terminal text or a file.

This is the Python face of the C++ rasterizer (``operations/render.hpp``):
deterministic, dependency-free, and the same pixels on every platform and
thread count. It needs no display, no GPU and no browser, so it works over SSH,
in a container or in CI, where :func:`view` and :func:`screenshot` cannot. See
``doc/tui.md``.

It is core-only by decision (roadmap 7.1.2): a NumPy twin of an integer
edge-function rasterizer would be slow and would double the surface that has to
stay byte-identical, so without the compiled extension these functions raise
``NotImplementedError`` rather than drawing something else.

Every function takes the render options as keyword arguments:

``width``, ``height``
    image size in pixels (``render_image`` and ``.png`` snapshots only; the
    text forms size the frame from ``cols`` and ``rows``).
``supersample``
    1, 2 or 4 samples per pixel along each axis, resolved by a box filter.
``view``
    ``"iso"``, ``"+x"``, ``"-x"``, ``"+y"``, ``"-y"``, ``"+z"`` or ``"-z"``
    (the camera sits on that side); overrides ``azimuth``/``elevation``.
``azimuth``, ``elevation``, ``roll``
    camera direction in degrees, the SVG writer's convention (the default is
    the isometric view).
``projection``, ``fov``
    ``"orthographic"`` (default) or ``"perspective"`` with a vertical field of
    view in degrees.
``zoom``, ``pan``
    magnification over the fit-to-frame view, and an ``(x, y)`` shift in
    fractions of the frame.
``shading``
    ``"none"``, ``"flat"`` (default) or ``"smooth"``; with ``two_sided``,
    ``ambient``, ``light_dir`` (camera space) and ``split_angle``.
``edges``
    ``"none"``, ``"all"`` or ``"feature"`` (with ``feature_angle``) drawn in
    ``edge_color``.
``fill_color``, ``line_color``, ``background``, ``nan_color``
    ``(r, g, b)`` or ``(r, g, b, a)`` tuples of 0-255; a background alpha of 0
    (the default) is transparent in a PNG and the terminal's own background
    in text.
``point_radius``
    disc radius of drawn points, in pixels.
``color_by``, ``component``, ``cmap``, ``vmin``, ``vmax``, ``colorbar``
    colour by a point or cell data array, as the SVG writer does.
``axes``, ``scale_bar``
    draw the world axes and a scale bar in the corners.
"""

from __future__ import annotations

import shutil
from typing import Optional

__all__ = ["render_image", "render_text", "snapshot"]

# Options that belong to the text encoding, not to the camera or the colours.
_TEXT_KEYS = (
    "encoding",
    "color_depth",
    "format",
    "cols",
    "rows",
    "cell_aspect",
    "cell_pixels",
    "tmux",
    "notes",
)
_SNAPSHOT_KEYS = ("png_compress", "cast_frames", "cast_fps", "cast_degrees")


def parse_color(text):
    """``#rrggbb``, ``#rrggbbaa`` or ``none``/``transparent`` as an RGBA tuple."""
    if text in ("none", "transparent"):
        return (0, 0, 0, 0)
    if not text.startswith("#") or len(text) not in (7, 9):
        raise ValueError(f"expected #rrggbb, #rrggbbaa or none, not '{text}'")
    try:
        values = [int(text[i : i + 2], 16) for i in range(1, len(text), 2)]
    except ValueError:
        raise ValueError(f"expected #rrggbb, #rrggbbaa or none, not '{text}'") from None
    return tuple(values) if len(values) == 4 else tuple(values) + (255,)


def _core_module(name):
    try:
        from . import _core
    except ImportError:
        _core = None
    if _core is None:
        raise NotImplementedError(
            f"meshio++: {name}: the software rasterizer is C++-core only and has "
            "no pure-Python fallback (its pixels are pinned byte for byte across "
            "platforms, which a second implementation could not keep) -- install "
            "a build with the compiled meshioplusplus._core extension"
        )
    return _core


def _render_dict(options):
    out = {}
    for key, value in options.items():
        if value is None:
            continue
        if key in ("edge_color", "fill_color", "line_color", "background", "nan_color"):
            value = [int(c) for c in value]
        elif key in ("pan", "light_dir"):
            value = [float(c) for c in value]
        out[key] = value
    return out


def _split(options):
    render = {}
    text = {}
    snap = {}
    for key, value in options.items():
        if key in _TEXT_KEYS:
            if value is not None:
                text[key] = value
        elif key in _SNAPSHOT_KEYS:
            if value is not None:
                snap[key] = value
        else:
            render[key] = value
    return _render_dict(render), text, snap


def _environment_color_depth(core):
    import os

    depth = core.detect_color_depth(
        os.environ.get("NO_COLOR"), os.environ.get("COLORTERM"), os.environ.get("TERM")
    )
    if depth == "16" and os.environ.get("WT_SESSION"):
        depth = "truecolor"  # Windows Terminal speaks 24-bit colour without saying so
    return depth


def render_image(
    mesh, width: int = 640, height: int = 480, return_info: bool = False, **options
):
    """Render a mesh into an RGBA image.

    :param mesh: the mesh to draw (unmodified). Volume cells are drawn through
        their boundary skin, 2-D cells as they are, lines as lines and vertices
        (or a cell-less point cloud) as discs.
    :param width: image width in pixels.
    :param height: image height in pixels.
    :param return_info: also return a dict with ``cell_ids`` (an ``(H, W)``
        int64 array: the input cell drawn at each pixel, global block-major,
        -1 for the background), ``colored``, ``vmin``, ``vmax`` and ``notes``
        (the text lines a terminal rendering prints under the picture).
    :param options: the render options listed in this module's docstring.
    :returns: an ``(H, W, 4)`` uint8 array with straight alpha, or
        ``(image, info)`` when ``return_info`` is set.
    :raises ValueError: on an invalid option or an unknown ``color_by`` array.
    :raises NotImplementedError: when the compiled core is unavailable.
    """
    core = _core_module("render_image")
    render, text, snap = _split(options)
    if text or snap:
        raise TypeError(
            "meshio++: render_image: "
            + ", ".join(sorted(list(text) + list(snap)))
            + " belong to render_text/snapshot, not to an image"
        )
    render["width"] = int(width)
    render["height"] = int(height)
    result = core.render(mesh, render)
    if not return_info:
        return result["image"]
    info = {k: result[k] for k in ("cell_ids", "colored", "vmin", "vmax", "notes")}
    return result["image"], info


def render_text(
    mesh,
    cols: Optional[int] = None,
    rows: Optional[int] = None,
    encoding: str = "halfblock",
    color_depth: Optional[str] = None,
    format: str = "ansi",
    **options,
):
    """Render a mesh as text for a terminal.

    :param mesh: the mesh to draw (unmodified).
    :param cols: width in terminal cells; the current terminal's by default
        (80 when there is none).
    :param rows: height in terminal cells; the current terminal's less one row
        by default (24 when there is none).
    :param encoding: ``"halfblock"`` (1x2 pixels per cell), ``"quadrant"``
        (2x2), ``"sextant"`` (2x3, Unicode 13), ``"braille"`` (2x4),
        ``"ascii"`` (a luminance ramp), or a graphics protocol: ``"kitty"``,
        ``"iterm2"`` or ``"sixel"``.
    :param color_depth: ``"truecolor"``, ``"256"``, ``"16"`` or ``"mono"``; by
        default what the environment advertises (``NO_COLOR``, ``COLORTERM``,
        ``TERM``).
    :param format: ``"ansi"`` (SGR colour escapes), ``"plain"`` (glyphs only)
        or ``"html"`` (a self-contained page).
    :param options: the render options listed in this module's docstring,
        plus ``cell_aspect`` (cell height over width, default 2),
        ``cell_pixels`` (``(w, h)`` of a cell for the graphics protocols),
        ``tmux`` (wrap graphics protocols for tmux passthrough) and ``notes``
        (print the colour range and scale under the picture, default true).
    :returns: the text, ready to print.
    :raises ValueError: on an invalid option or an unknown ``color_by`` array.
    :raises NotImplementedError: when the compiled core is unavailable.
    """
    core = _core_module("render_text")
    render, text, snap = _split(options)
    if snap:
        raise TypeError(
            "meshio++: render_text: " + ", ".join(sorted(snap)) + " belong to snapshot"
        )
    if cols is None or rows is None:
        size = shutil.get_terminal_size((80, 25))
        cols = size.columns if cols is None else cols
        rows = max(1, size.lines - 1) if rows is None else rows
    text.update(
        {
            "cols": int(cols),
            "rows": int(rows),
            "encoding": encoding,
            "color_depth": color_depth or _environment_color_depth(core),
            "format": format,
        }
    )
    return core.render_text(mesh, render, text)


def snapshot(
    mesh,
    path,
    width: int = 800,
    height: int = 600,
    cols: int = 100,
    rows: int = 40,
    encoding: str = "halfblock",
    color_depth: str = "truecolor",
    **options,
):
    """Render a mesh to a file chosen by its extension.

    ``.png`` writes an RGBA image of ``width`` x ``height`` pixels with no
    dependency (stored deflate blocks, the same bytes everywhere; pass
    ``png_compress=1..9`` to compress through zlib). ``.txt`` (plain cells),
    ``.ansi`` (cells with SGR colour, for ``cat``) and ``.html`` (a
    self-contained page) draw ``cols`` x ``rows`` terminal cells. ``.cast``
    records an asciinema v2 orbit around the mesh (``cast_frames``,
    ``cast_fps``, ``cast_degrees``), timed by frame index, never the clock.

    :param mesh: the mesh to draw (unmodified).
    :param path: the output path.
    :param options: the render options listed in this module's docstring, the
        text options of :func:`render_text`, and the ``png_compress`` /
        ``cast_*`` options above.
    :returns: ``str(path)``.
    :raises ValueError: on an unknown extension, an invalid option or an
        unknown ``color_by`` array.
    :raises NotImplementedError: when the compiled core is unavailable.
    """
    import os

    core = _core_module("snapshot")
    render, text, snap = _split(options)
    render["width"] = int(width)
    render["height"] = int(height)
    text.update(
        {
            "cols": int(cols),
            "rows": int(rows),
            "encoding": encoding,
            "color_depth": color_depth,
        }
    )
    core.write_snapshot(os.fspath(path), mesh, render, text, snap)
    return str(path)
