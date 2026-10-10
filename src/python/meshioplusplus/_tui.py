"""The interactive terminal viewer.

:func:`tui` takes over the terminal and lets you orbit, zoom and pan a mesh
with the mouse and the keyboard: the software rasterizer of
:mod:`meshioplusplus._render`, driven by a loop that redraws only the terminal
cells that changed. It is the viewer for a machine with no display: an SSH
session, a container, a text console. ``view()`` picks it automatically when
there is no window system but a terminal; see :doc:`/tui`.

The loop is C++ (the same code as the native ``meshioplusplus tui``), so it
needs the compiled ``meshioplusplus._core`` extension and has no pure-Python
fallback.
"""

from __future__ import annotations

import os
import sys
from typing import Optional

from ._render import _core_module, _environment_color_depth, _split

__all__ = ["tui"]

_SIGINT_EXIT = 128 + 2


def tui(
    mesh,
    *,
    encoding: str = "halfblock",
    color_depth: Optional[str] = None,
    title: Optional[str] = None,
    replay: Optional[bytes] = None,
    cols: int = 100,
    rows: int = 40,
    **options,
) -> dict:
    """Browse a mesh in the terminal until the user quits.

    Keys: drag with the mouse (or the arrow keys) to orbit and pan, the wheel
    or ``+``/``-`` to zoom, ``1``..``7`` for the views ``iso``, ``+x``, ``-x``,
    ``+y``, ``-y``, ``+z``, ``-z``, ``p`` for perspective, ``e`` for edges,
    ``s`` for shading, ``a``/``b``/``c`` for the axes, scale bar and colour
    bar, ``r`` to reset, ``?`` for help, ``q`` (or Escape, or Ctrl-C) to quit.
    The terminal is restored however the viewer ends.

    :param mesh: a mesh, or the path of a file to read.
    :param encoding: how the picture becomes cells: ``"halfblock"`` (default),
        ``"quadrant"``, ``"sextant"``, ``"braille"``, ``"ascii"``, or a graphics
        protocol (``"kitty"``, ``"iterm2"``, ``"sixel"``), which redraws the whole
        image on every change.
    :param color_depth: ``"truecolor"``, ``"256"``, ``"16"`` or ``"mono"``; by
        default what the environment advertises.
    :param title: shown in the status line (the file name for a path).
    :param replay: a recorded input stream (bytes) to play on a ``cols`` x
        ``rows`` screen instead of the terminal; the result then carries
        ``"output"``, everything the loop wrote. The tests and the
        documentation figures use it, and it needs no terminal.
    :param options: the render options of :func:`~meshioplusplus.render_image`
        (``view``, ``color_by``, ``cmap``, ``edges``, ``shading``, ...); the
        frame size comes from the terminal.
    :returns: a dict: ``exit`` (0 for a normal quit), ``frames`` drawn, the
        final ``azimuth``, ``elevation``, ``zoom``, ``pan_x``, ``pan_y``, the last
        ``status`` line, and ``output`` for a replay.
    :raises RuntimeError: when standard input or output is not a terminal
        (use :func:`~meshioplusplus.snapshot` or
        :func:`~meshioplusplus.render_text` for one frame), or the terminal
        cannot be switched to raw mode.
    :raises KeyboardInterrupt: when the process received SIGINT.
    :raises ValueError: on an invalid option or an unknown ``color_by`` array.
    :raises NotImplementedError: when the compiled core is unavailable.
    """
    core = _core_module("tui")
    render, text, snap = _split(options)
    if snap:
        raise TypeError(
            "meshio++: tui: " + ", ".join(sorted(snap)) + " belong to snapshot"
        )
    for key in ("cols", "rows", "format", "notes"):
        text.pop(key, None)  # the screen sets the size; the loop prints the notes
    if isinstance(mesh, (str, os.PathLike)):
        from ._helpers import read as _read

        title = title or os.path.basename(os.fspath(mesh))
        mesh = _read(os.fspath(mesh))
    if replay is None and not (sys.stdin.isatty() and sys.stdout.isatty()):
        raise RuntimeError(
            "meshio++: tui: needs a terminal on standard input and output; "
            "use snapshot() or render_text() for one frame to a file or a string"
        )
    text.update(
        {
            "encoding": encoding,
            "color_depth": color_depth or _environment_color_depth(core),
        }
    )
    result = core.tui(
        mesh,
        render,
        text,
        title or "",
        None if replay is None else bytes(replay),
        int(cols),
        int(rows),
    )
    if result["exit"] == _SIGINT_EXIT:
        raise KeyboardInterrupt
    return result
