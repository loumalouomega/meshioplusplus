"""
Fixed-width card tokenizer shared by the LS-DYNA keyword reader and writer.

An LS-DYNA card is one line whose fields sit in fixed columns. Four layouts exist
and are told apart per card:

* **standard**: the widths a card declares (``I8``, ``E16``, ``I10`` ...);
* **long** (``*KEYWORD LONG=Y``, or a ``+`` suffix on one keyword): every field,
  integer and real, is 20 columns wide;
* **i10** (``*KEYWORD I10=Y``): only the 8-column integer fields widen to 10, which
  is a different mechanism from long;
* **free**: a line that contains a comma is split on commas, whatever the mode.

A layout is a sequence of ``(kind, width)`` pairs, ``kind`` being ``"i"`` or ``"r"``.
The C++ twin is ``detail/keyword_card.hpp``; keep the two in step.
"""

import re

from .._exceptions import ReadError

STD = "std"
LONG = "long"
I10 = "i10"

# The layouts of the cards the reader looks at.
NODE = (("i", 8), ("r", 16), ("r", 16), ("r", 16), ("r", 8), ("r", 8))
ELEMENT = tuple(("i", 8) for _ in range(10))
ELEMENT_MASS = (("i", 8), ("i", 8), ("r", 16), ("i", 8))
IDS10 = tuple(("i", 10) for _ in range(8))
SEGMENT = (("i", 10),) * 4 + (("r", 10),) * 4

_EXPONENT_WITHOUT_LETTER = re.compile(r"^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d+)$")


def field_width(kind, width, mode):
    """The column width of one field under ``mode``."""
    if mode == LONG:
        return 20
    if mode == I10 and kind == "i" and width == 8:
        return 10
    return width


def is_free(line):
    """A comma anywhere on a data line switches that card to free format."""
    return "," in line


def split_card(line, layout, mode):
    """The stripped text of each field of ``line``; a missing field is ``""``.

    Free-format lines return every comma-separated entry, so a card with a variable
    number of fields (a node list) is not truncated to the layout.
    """
    if is_free(line):
        parts = [p.strip() for p in line.split(",")]
        if len(parts) < len(layout):
            parts.extend([""] * (len(layout) - len(parts)))
        return parts
    out = []
    pos = 0
    for kind, width in layout:
        w = field_width(kind, width, mode)
        out.append(line[pos : pos + w].strip())
        pos += w
    return out


def to_int(text, where=""):
    """An integer field; blank is 0."""
    if not text:
        return 0
    try:
        return int(text)
    except ValueError:
        raise ReadError(f"LS-DYNA: invalid integer field {text!r}{where}") from None


def to_float(text, where=""):
    """A real field; blank is 0. Accepts ``D`` exponents and ``1.5-3`` (no letter)."""
    if not text:
        return 0.0
    s = text.replace("D", "E").replace("d", "e")
    m = _EXPONENT_WITHOUT_LETTER.match(s)
    if m:
        s = m.group(1) + "e" + m.group(2)
    try:
        return float(s)
    except ValueError:
        raise ReadError(f"LS-DYNA: invalid real field {text!r}{where}") from None


def format_real16(x):
    """``x`` in at most 16 columns: the shortest scientific string that round-trips,
    else as many digits as fit. Twin of ``kwc_format_real16``."""
    x = float(x)
    if x == 0.0:
        return "0.0"
    neg = 1 if x < 0.0 else 0
    e3 = 1 if abs(x) >= 1e100 or abs(x) < 1e-99 else 0
    pmax = 10 - neg - e3
    s = ""
    for p in range(1, pmax + 1):
        s = f"{x:.{p}e}"
        if float(s) == x:
            return s
    return s


def pack_card(values, layout):
    """Standard-format line for ``values`` (ints or preformatted strings)."""
    out = []
    for (_, width), v in zip(layout, values):
        out.append(f"{v:>{width}}")
    return "".join(out)


class _FormatParser:
    """Recursive-descent parse of a Fortran edit-descriptor list (the twin of
    keyword_card.cpp's ``KcFormatParser``)."""

    def __init__(self, text):
        self.text = text
        self.pos = 0

    def fail(self):
        raise ReadError(f"cannot parse Fortran format '{self.text}'")

    def skip(self):
        while self.pos < len(self.text) and self.text[self.pos] in " \t":
            self.pos += 1

    def number(self):
        self.skip()
        start = self.pos
        while self.pos < len(self.text) and self.text[self.pos].isdigit():
            self.pos += 1
        if start == self.pos:
            return -1
        value = int(self.text[start : self.pos])
        if value > 100000:
            self.fail()
        return value

    def items(self, out, nested):
        while True:
            self.skip()
            if self.pos >= len(self.text):
                if nested:
                    self.fail()
                return
            if self.text[self.pos] == ")":
                if not nested:
                    self.fail()
                self.pos += 1
                return
            self.item(out)
            self.skip()
            if self.pos < len(self.text) and self.text[self.pos] == ",":
                self.pos += 1

    def item(self, out):
        repeat = self.number()
        count = 1 if repeat < 0 else repeat
        self.skip()
        if self.pos < len(self.text) and self.text[self.pos] == "(":
            self.pos += 1
            group = []
            self.items(group, True)
            out.extend(group * count)
            return
        self.skip()
        if self.pos >= len(self.text):
            self.fail()
        c = self.text[self.pos].lower()
        self.pos += 1
        if c == "p":
            if repeat < 0:
                self.fail()
            return
        if c == "x":
            out.append(("x", count))
            return
        kind = {"i": "i", "a": "a", "e": "r", "d": "r", "f": "r", "g": "r"}.get(c)
        if kind is None:
            self.fail()
        if (
            c == "e"
            and self.pos < len(self.text)
            and self.text[self.pos].lower() in "sn"
        ):
            self.pos += 1
        width = self.number()
        if width <= 0:
            self.fail()
        self.skip()
        if self.pos < len(self.text) and self.text[self.pos] == ".":
            self.pos += 1
            if self.number() < 0:
                self.fail()
            self.skip()
            if (
                self.pos < len(self.text)
                and self.text[self.pos].lower() == "e"
                and kind == "r"
            ):
                self.pos += 1
                if self.number() < 0:
                    self.fail()
        out.extend([(kind, width)] * count)


def parse_fortran_format(text):
    """The ``(kind, width)`` fields of a Fortran edit-descriptor list such as
    ``"(1i7,2i9,6e21.13e3)"``: ``i`` integer, ``r`` real, ``a`` characters, ``x``
    skipped columns. Repeat counts and groups expand; a scale factor is ignored."""
    parser = _FormatParser(text)
    parser.skip()
    out = []
    if parser.pos < len(text) and text[parser.pos] == "(":
        parser.pos += 1
        parser.items(out, True)
        parser.skip()
        if parser.pos != len(text):
            parser.fail()
    else:
        parser.items(out, False)
    if not out:
        parser.fail()
    return out


def split_fixed(line, fields):
    """The stripped text of each field ``fields`` lays out on ``line``, in fixed
    columns; stops at the end of the line and returns no ``x`` fields."""
    line = line.rstrip("\r\n")
    out = []
    col = 0
    for kind, width in fields:
        if col >= len(line):
            break
        if kind != "x":
            out.append(line[col : col + width].strip(" \t"))
        col += width
    return out
