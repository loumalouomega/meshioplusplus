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
    mesh=None,
    *,
    encoding: str = "halfblock",
    color_depth: Optional[str] = None,
    title: Optional[str] = None,
    replay: Optional[bytes] = None,
    cols: int = 100,
    rows: int = 40,
    compare=None,
    diff: bool = False,
    shared_range: bool = True,
    series=None,
    input_format: Optional[str] = None,
    follow=False,
    settle: int = 300,
    fps: float = 4.0,
    session: Optional[str] = None,
    **options,
) -> dict:
    """Browse a mesh in the terminal until the user quits.

    Keys: drag with the mouse (or the arrow keys) to orbit and pan, the wheel
    or ``+``/``-`` to zoom, ``1``..``7`` for the views ``iso``, ``+x``, ``-x``,
    ``+y``, ``-y``, ``+z``, ``-z``, ``p`` for perspective, ``e`` for edges,
    ``s`` for shading, ``a``/``b``/``c`` for the axes, scale bar and colour
    bar, ``r`` to reset, ``?`` for help, ``q`` (or Escape, or Ctrl-C) to quit.
    A click probes the cell under the pointer, ``i`` pins it (two pins show
    their difference) and ``0`` clears; ``x``/``y``/``z`` cut a plane away at
    the middle of the model, ``,`` and ``.`` slide it; ``[`` ``]`` step through
    a series, ``{`` ``}`` by ten, space plays; ``:`` opens a command line that
    takes the flags of ``snapshot`` (``:cmap turbo``, ``:range 0 1``,
    ``:clip +x 0.5``, ``:w shot.png``, ``:session save view.json``, ``:help``).
    The terminal is restored however the viewer ends.

    :param mesh: a mesh, or the path of a file to read; ``None`` with ``series``.
    :param encoding: how the picture becomes cells: ``"halfblock"`` (default),
        ``"quadrant"``, ``"sextant"``, ``"braille"``, ``"ascii"``, or a graphics
        protocol (``"kitty"``, ``"iterm2"``, ``"sixel"``), which redraws the whole
        image on every change and cannot compare two meshes.
    :param color_depth: ``"truecolor"``, ``"256"``, ``"16"`` or ``"mono"``; by
        default what the environment advertises.
    :param title: shown in the status line (the file name for a path).
    :param replay: a recorded input stream (bytes) to play on a ``cols`` x
        ``rows`` screen instead of the terminal; the result then carries
        ``"output"``, everything the loop wrote. The tests and the
        documentation figures use it, and it needs no terminal.
    :param compare: a second mesh (or path) drawn beside the first under one
        camera, with one colour range over both unless ``shared_range`` is
        false; ``diff`` draws ``|compare - mesh|`` of the ``color_by`` point
        array instead of the second mesh (the two share their nodes).
    :param series: a list of paths, or a glob (``"out_*.vtu"``, natural-numeric
        order), shown as a time series one step at a time with a colour range
        fixed from the first step; ``input_format`` forces the format.
    :param follow: poll the series for a new step every this many milliseconds
        (``True`` is 500) and show it once its file has stopped changing for
        ``settle`` ms; a half-written file keeps the last picture.
    :param fps: steps per second when playing.
    :param session: a JSON session file, read at the start if it exists and
        written when the viewer ends.
    :param options: the render options of :func:`~meshioplusplus.render_image`
        (``view``, ``color_by``, ``cmap``, ``edges``, ``shading``, ``cutaway``,
        ...); the frame size comes from the terminal.
    :returns: a dict: ``exit`` (0 for a normal quit), ``frames`` drawn, the
        final ``azimuth``, ``elevation``, ``zoom``, ``pan_x``, ``pan_y``, the
        series ``step``, the last ``status`` line and ``probe`` lines, and
        ``output`` for a replay.
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
    from ._helpers import read as _read

    if isinstance(mesh, (str, os.PathLike)):
        title = title or os.path.basename(os.fspath(mesh))
        mesh = _read(os.fspath(mesh), file_format=input_format)
    paths: list = []
    pattern = ""
    if series is not None:
        if isinstance(series, (str, os.PathLike)):
            pattern = os.fspath(series)
            title = title or pattern
        else:
            paths = [os.fspath(p) for p in series]
            title = title or (os.path.basename(paths[0]) if paths else "")
        mesh = None
    elif mesh is None:
        raise TypeError("meshio++: tui: give a mesh, a path, or a series")
    compare_title = ""
    if isinstance(compare, (str, os.PathLike)):
        compare_title = os.path.basename(os.fspath(compare))
        compare = _read(os.fspath(compare), file_format=input_format)
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
    follow_ms = 500 if follow is True else int(follow or 0)
    result = core.tui(
        mesh,
        render,
        text,
        title or "",
        None if replay is None else bytes(replay),
        int(cols),
        int(rows),
        compare,
        compare_title,
        bool(diff),
        bool(shared_range),
        paths,
        pattern,
        input_format or "",
        follow_ms,
        int(settle),
        float(fps),
        os.fspath(session) if session else "",
    )
    if result["exit"] == _SIGINT_EXIT:
        raise KeyboardInterrupt
    return result
