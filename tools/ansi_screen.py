#!/usr/bin/env python3
"""A tiny terminal emulator for what ``meshioplusplus tui`` writes.

It understands exactly the escape sequences the viewer and ``snapshot`` emit --
cursor addressing (``CSI row;col H``), clear (``CSI 2 J``), the SGR colours
(truecolor, the 256-colour palette, the 16 ANSI colours, reverse video) and
text -- and ignores the rest (mode switches such as the alternate screen,
mouse reporting and bracketed paste, graphics-protocol payloads).

Two users:

* the tests, which replay a recorded input stream through the loop and assert
  on the final screen (``Screen.text()``, ``Screen.cell()``), without a real
  terminal; and
* ``tools/gen_doc_images.py``, which turns that screen into the documentation
  screenshots (``Screen.to_svg()``): block, quadrant, sextant and Braille glyphs
  are drawn as shapes, so the figure looks the same in every viewer and needs no
  font.

No dependency beyond the standard library.
"""

from __future__ import annotations

import html
from dataclasses import dataclass
from typing import List, Optional, Tuple

RGB = Tuple[int, int, int]

# The default xterm 16-colour palette.
_ANSI16: List[RGB] = [
    (0, 0, 0),
    (205, 49, 49),
    (13, 188, 121),
    (229, 229, 16),
    (36, 114, 200),
    (188, 63, 188),
    (17, 168, 205),
    (229, 229, 229),
    (102, 102, 102),
    (241, 76, 76),
    (35, 209, 139),
    (245, 245, 67),
    (59, 142, 234),
    (214, 112, 214),
    (41, 184, 219),
    (255, 255, 255),
]


def _xterm256(n: int) -> RGB:
    if n < 16:
        return _ANSI16[n]
    if n < 232:
        n -= 16
        levels = [0, 95, 135, 175, 215, 255]
        return levels[n // 36], levels[(n // 6) % 6], levels[n % 6]
    v = 8 + 10 * (n - 232)
    return v, v, v


@dataclass
class Cell:
    char: str = " "
    fg: Optional[RGB] = None  # None: the terminal's own foreground
    bg: Optional[RGB] = None
    reverse: bool = False

    def colors(self, default_fg: RGB, default_bg: RGB) -> Tuple[RGB, RGB]:
        fg = self.fg if self.fg is not None else default_fg
        bg = self.bg if self.bg is not None else default_bg
        return (bg, fg) if self.reverse else (fg, bg)


class Screen:
    """A ``cols`` x ``rows`` screen that bytes can be fed to."""

    def __init__(self, cols: int, rows: int):
        self.cols = cols
        self.rows = rows
        self.cells = [[Cell() for _ in range(cols)] for _ in range(rows)]
        self.row = 0
        self.col = 0
        self._fg: Optional[RGB] = None
        self._bg: Optional[RGB] = None
        self._reverse = False
        self._pending = b""

    # ------------------------------------------------------------------ input
    def feed(self, data: bytes) -> "Screen":
        data = self._pending + data
        self._pending = b""
        text = data.decode("utf-8", errors="replace")
        i = 0
        n = len(text)
        while i < n:
            ch = text[i]
            if ch == "\x1b":
                if i + 1 >= n:
                    self._pending = text[i:].encode("utf-8")
                    break
                nxt = text[i + 1]
                if nxt == "[":
                    j = i + 2
                    while j < n and not ("@" <= text[j] <= "~"):
                        j += 1
                    if j >= n:
                        self._pending = text[i:].encode("utf-8")
                        break
                    self._csi(text[i + 2 : j], text[j])
                    i = j + 1
                    continue
                if nxt in "_]P^X":  # APC/OSC/DCS...: skip to the string terminator
                    j = i + 2
                    while j < n and text[j] != "\x07" and text[j : j + 2] != "\x1b\\":
                        j += 1
                    if j >= n:
                        self._pending = text[i:].encode("utf-8")
                        break
                    i = j + (1 if text[j] == "\x07" else 2)
                    continue
                i += 2
                continue
            if ch == "\r":
                self.col = 0
            elif ch == "\n":
                self.row = min(self.row + 1, self.rows - 1)
            elif ch >= " ":
                self._put(ch)
            i += 1
        return self

    def _put(self, ch: str) -> None:
        if 0 <= self.row < self.rows and 0 <= self.col < self.cols:
            self.cells[self.row][self.col] = Cell(ch, self._fg, self._bg, self._reverse)
        self.col += 1

    def _csi(self, body: str, final: str) -> None:
        if body.startswith("?") or body.startswith("<") or body.startswith(">"):
            return  # a private mode switch or a report: nothing to draw
        params = [int(p) if p.isdigit() else 0 for p in body.split(";")] if body else []
        if final in "Hf":
            row = params[0] if len(params) > 0 and params[0] else 1
            col = params[1] if len(params) > 1 and params[1] else 1
            self.row = min(max(row - 1, 0), self.rows - 1)
            self.col = min(max(col - 1, 0), self.cols - 1)
        elif final == "J":
            if (params[0] if params else 0) == 2:
                self.cells = [
                    [Cell() for _ in range(self.cols)] for _ in range(self.rows)
                ]
        elif final == "K":
            mode = params[0] if params else 0
            line = self.cells[self.row]
            if mode == 2:
                lo, hi = 0, self.cols
            elif mode == 1:
                lo, hi = 0, self.col + 1
            else:
                lo, hi = self.col, self.cols
            for c in range(lo, hi):
                line[c] = Cell(" ", None, self._bg, False)
        elif final == "m":
            self._sgr(params or [0])

    def _sgr(self, params: List[int]) -> None:
        i = 0
        while i < len(params):
            p = params[i]
            if p == 0:
                self._fg = self._bg = None
                self._reverse = False
            elif p == 7:
                self._reverse = True
            elif p == 27:
                self._reverse = False
            elif 30 <= p <= 37:
                self._fg = _ANSI16[p - 30]
            elif 90 <= p <= 97:
                self._fg = _ANSI16[p - 90 + 8]
            elif 40 <= p <= 47:
                self._bg = _ANSI16[p - 40]
            elif 100 <= p <= 107:
                self._bg = _ANSI16[p - 100 + 8]
            elif p == 39:
                self._fg = None
            elif p == 49:
                self._bg = None
            elif p in (38, 48):
                color = None
                if i + 1 < len(params) and params[i + 1] == 5 and i + 2 < len(params):
                    color = _xterm256(params[i + 2] & 255)
                    i += 2
                elif i + 1 < len(params) and params[i + 1] == 2 and i + 4 < len(params):
                    color = tuple(c & 255 for c in params[i + 2 : i + 5])
                    i += 4
                if color is not None:
                    if p == 38:
                        self._fg = color
                    else:
                        self._bg = color
            i += 1

    # ----------------------------------------------------------------- output
    def cell(self, row: int, col: int) -> Cell:
        """The cell at 0-based (row, col)."""
        return self.cells[row][col]

    def text(self) -> List[str]:
        """The characters of each row, trailing blanks removed."""
        return ["".join(c.char for c in line).rstrip() for line in self.cells]

    def non_blank_rows(self) -> int:
        return sum(1 for line in self.text() if line)

    def to_svg(
        self,
        cell_w: int = 8,
        cell_h: int = 16,
        default_fg: RGB = (212, 212, 212),
        default_bg: RGB = (24, 24, 28),
        title: str = "",
    ) -> str:
        """The screen as a self-contained SVG: block, quadrant, sextant and
        Braille glyphs as shapes, other text as ``<text>``."""
        width = self.cols * cell_w
        height = self.rows * cell_h
        out = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" role="img"'
            + (f' aria-label="{html.escape(title, quote=True)}"' if title else "")
            + ">",
            f'<rect width="{width}" height="{height}" fill="{_hex(default_bg)}"/>',
        ]
        if title:
            out.insert(1, f"<title>{html.escape(title)}</title>")
        text_runs: List[str] = []
        pending: List[Tuple[int, str, RGB]] = []
        for r, line in enumerate(self.cells):
            # Backgrounds first, merged across neighbours of the same colour.
            c = 0
            while c < len(line):
                bg = line[c].colors(default_fg, default_bg)[1]
                run = 1
                while (
                    c + run < len(line)
                    and line[c + run].colors(default_fg, default_bg)[1] == bg
                ):
                    run += 1
                if bg != default_bg:
                    out.append(_rect(c * cell_w, r * cell_h, cell_w * run, cell_h, bg))
                c += run
            c = 0
            while c < len(line):
                cell = line[c]
                # A run of identical cells draws as one wide shape where the glyph
                # spans the full cell width (blocks, half blocks, blanks).
                run = 1
                while c + run < len(line) and line[c + run] == cell:
                    run += 1
                fg, bg = cell.colors(default_fg, default_bg)
                x = c * cell_w
                y = r * cell_h
                shapes = _glyph_shapes(cell.char, cell_w, cell_h)
                if shapes is not None:
                    for sx, sy, sw, sh in shapes:
                        if sx == 0 and sw == cell_w:
                            out.append(_rect(x, y + sy, cell_w * run, sh, fg))
                        else:
                            for k in range(run):
                                out.append(
                                    _rect(x + k * cell_w + sx, y + sy, sw, sh, fg)
                                )
                elif cell.char.strip() == "" or cell.char == "\u2800":
                    pass
                elif "\u2801" <= cell.char <= "\u28ff":
                    for k in range(run):
                        for cx, cy, rad in _braille_dots(
                            ord(cell.char) - 0x2800, cell_w, cell_h
                        ):
                            out.append(
                                f'<circle cx="{x + k * cell_w + cx:.2f}" cy="{y + cy:.2f}" '
                                f'r="{rad:.2f}" fill="{_hex(fg)}"/>'
                            )
                else:
                    for k in range(run):
                        pending.append((c + k, cell.char, fg))
                c += run
            text_runs.extend(
                _text_runs(
                    self.cells[r], pending, r, cell_w, cell_h, default_fg, default_bg
                )
            )
            pending = []
        if text_runs:
            out.append(
                f'<g font-family="DejaVu Sans Mono, Menlo, Consolas, monospace" '
                f'font-size="{cell_h - 4}" xml:space="preserve">'
            )
            out.extend(text_runs)
            out.append("</g>")
        out.append("</svg>")
        return "\n".join(out) + "\n"


def _text_runs(line, pending, row, cw, ch, default_fg, default_bg):
    """One ``<text>`` per run of neighbouring characters of the same colour
    (spaces between words stay in the run), forced onto the cell grid with
    ``textLength`` so any monospace font lines up."""
    out = []
    i = 0
    while i < len(pending):
        col, _, fg = pending[i]
        j = i
        while (
            j + 1 < len(pending)
            and pending[j + 1][2] == fg
            and pending[j + 1][0] - pending[j][0] <= 2
            and all(
                line[k].char == " " and _glyph_shapes(line[k].char, cw, ch) is None
                for k in range(pending[j][0] + 1, pending[j + 1][0])
            )
        ):
            j += 1
        first, last = pending[i][0], pending[j][0]
        chars = [" "] * (last - first + 1)
        for k in range(i, j + 1):
            chars[pending[k][0] - first] = pending[k][1]
        out.append(
            f'<text x="{first * cw}" y="{row * ch + ch - 4}" fill="{_hex(fg)}" '
            f'textLength="{len(chars) * cw}" lengthAdjust="spacing">'
            f"{html.escape(''.join(chars))}</text>"
        )
        i = j + 1
    return out


def _hex(color: RGB) -> str:
    return "#%02x%02x%02x" % tuple(color)


def _rect(x: float, y: float, w: float, h: float, color: RGB) -> str:
    return (
        f'<rect x="{_num(x)}" y="{_num(y)}" width="{_num(w)}" height="{_num(h)}" '
        f'fill="{_hex(color)}"/>'
    )


def _num(v: float) -> str:
    return f"{v:g}"


def _grid_shapes(
    mask_cells: List[Tuple[int, int]], cols: int, rows: int, cw: int, ch: int
):
    """The sub-rectangles of a ``cols`` x ``rows`` grid of subpixels."""
    w = cw / cols
    h = ch / rows
    return [(c * w, r * h, w, h) for r, c in mask_cells]


def _glyph_shapes(char: str, cw: int, ch: int):
    """The rectangles that draw a block-element, quadrant or sextant glyph, or
    ``None`` for any other character."""
    cp = ord(char)
    if cp == 0x2588:
        return [(0, 0, cw, ch)]
    if cp == 0x2580:
        return [(0, 0, cw, ch / 2)]
    if cp == 0x2584:
        return [(0, ch / 2, cw, ch / 2)]
    if cp == 0x258C:
        return [(0, 0, cw / 2, ch)]
    if cp == 0x2590:
        return [(cw / 2, 0, cw / 2, ch)]
    quadrants = {
        0x2596: [(1, 0)],
        0x2597: [(1, 1)],
        0x2598: [(0, 0)],
        0x2599: [(0, 0), (1, 0), (1, 1)],
        0x259A: [(0, 0), (1, 1)],
        0x259B: [(0, 0), (0, 1), (1, 0)],
        0x259C: [(0, 0), (0, 1), (1, 1)],
        0x259D: [(0, 1)],
        0x259E: [(0, 1), (1, 0)],
        0x259F: [(0, 1), (1, 0), (1, 1)],
    }
    if cp in quadrants:
        return _grid_shapes(quadrants[cp], 2, 2, cw, ch)
    if 0x1FB00 <= cp <= 0x1FB3B:
        mask = cp - 0x1FB00 + 1
        if mask >= 21:
            mask += 1
        if mask >= 42:
            mask += 1
        # bit i: row i // 2, column i % 2 of a 2x3 grid
        cells = [(i // 2, i % 2) for i in range(6) if mask & (1 << i)]
        return _grid_shapes(cells, 2, 3, cw, ch)
    return None


def _braille_dots(bits: int, cw: int, ch: int):
    # dot number -> (column, row) in the 2x4 cell
    position = {
        1: (0, 0),
        2: (0, 1),
        3: (0, 2),
        4: (1, 0),
        5: (1, 1),
        6: (1, 2),
        7: (0, 3),
        8: (1, 3),
    }
    rad = min(cw / 4.5, ch / 9.0)
    dots = []
    for dot in range(1, 9):
        if bits & (1 << (dot - 1)):
            col, row = position[dot]
            dots.append(((col + 0.5) * cw / 2, (row + 0.5) * ch / 4, rad))
    return dots


def replay(output: bytes, cols: int, rows: int) -> Screen:
    """The screen after ``output`` (what the viewer wrote) on a terminal of that size."""
    return Screen(cols, rows).feed(output)


if __name__ == "__main__":  # pragma: no cover - a convenience for ad-hoc figures
    import sys

    if len(sys.argv) != 4:
        raise SystemExit("usage: ansi_screen.py COLS ROWS OUT.svg  (reads stdin)")
    scr = replay(sys.stdin.buffer.read(), int(sys.argv[1]), int(sys.argv[2]))
    with open(sys.argv[3], "w", encoding="utf-8") as handle:
        handle.write(scr.to_svg())
