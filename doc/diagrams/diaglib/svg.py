"""A small deterministic SVG builder.

Every primitive formats its coordinates with two decimals and emits its
attributes in a fixed order, so the same figure code always yields the same
bytes -- that is what lets the committed SVGs be checked for staleness.
"""

from xml.sax.saxutils import escape

from . import palette as P


def fmt(v):
    """Two decimals, trailing zeros trimmed: stable and readable."""
    s = f"{float(v):.2f}"
    s = s.rstrip("0").rstrip(".") if "." in s else s
    return "0" if s in ("-0", "") else s


def _attrs(pairs):
    return "".join(f' {k}="{v}"' for k, v in pairs if v is not None)


class Canvas:
    """Accumulates SVG elements; ``render()`` returns the document."""

    def __init__(self, width, height, title):
        self.w = float(width)
        self.h = float(height)
        self.title = title
        self._parts = []
        self._markers = {}
        self.paper()

    # -- structure -----------------------------------------------------------
    def paper(self):
        """The light frame every figure sits on (the dark-mode rule)."""
        self.rrect(
            0.75,
            0.75,
            self.w - 1.5,
            self.h - 1.5,
            fill=P.PAPER,
            stroke=P.PAPER_BORDER,
            rx=8,
            sw=1.5,
        )

    def _marker_id(self, color):
        key = color.lstrip("#")
        if key not in self._markers:
            self._markers[key] = (
                f'<marker id="arrow-{key}" viewBox="0 0 10 10" refX="9" refY="5"'
                f' markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
                f'<path d="M 0 0 L 10 5 L 0 10 z" fill="{color}"/></marker>'
            )
        return f"arrow-{key}"

    def raw(self, text):
        self._parts.append(text)

    def group(self, opacity=None, transform=None):
        """Context manager wrapping subsequent elements in a ``<g>``."""
        canvas = self

        class _G:
            def __enter__(self_):
                canvas._parts.append(
                    "<g"
                    + _attrs([("opacity", opacity), ("transform", transform)])
                    + ">"
                )
                return canvas

            def __exit__(self_, *exc):
                canvas._parts.append("</g>")
                return False

        return _G()

    # -- primitives ----------------------------------------------------------
    def rrect(
        self,
        x,
        y,
        w,
        h,
        fill=P.WHITE,
        stroke=P.INK_2,
        rx=6,
        sw=P.STROKE_BOX,
        dash=None,
        fill_opacity=None,
        opacity=None,
    ):
        self._parts.append(
            "<rect"
            + _attrs(
                [
                    ("x", fmt(x)),
                    ("y", fmt(y)),
                    ("width", fmt(w)),
                    ("height", fmt(h)),
                    ("rx", fmt(rx) if rx else None),
                    ("fill", fill),
                    ("fill-opacity", fmt(fill_opacity) if fill_opacity else None),
                    ("stroke", stroke),
                    ("stroke-width", fmt(sw) if stroke and stroke != "none" else None),
                    ("stroke-dasharray", dash),
                    ("opacity", fmt(opacity) if opacity else None),
                ]
            )
            + "/>"
        )

    def box(
        self,
        x,
        y,
        w,
        h,
        label,
        color=P.CORE,
        sub=None,
        size=P.SIZE_LABEL,
        sub_size=P.SIZE_SMALL,
        mono=False,
        dash=None,
        fill=None,
        fill_opacity=P.FILL_OPACITY,
        rx=6,
        weight="600",
        text_color=P.INK,
        sw=P.STROKE_BOX,
    ):
        """A tinted, outlined box with a centred label and optional sublabel."""
        self.rrect(
            x,
            y,
            w,
            h,
            fill=fill or color,
            fill_opacity=None if fill else fill_opacity,
            stroke=color,
            rx=rx,
            dash=dash,
            sw=sw,
        )
        lines = label if isinstance(label, (list, tuple)) else [label]
        subs = [] if sub is None else (sub if isinstance(sub, (list, tuple)) else [sub])
        lh = size * 1.25
        sh = sub_size * 1.25
        total = lh * len(lines) + sh * len(subs)
        cy = y + h / 2 - total / 2 + lh * 0.78
        for line in lines:
            self.text(
                x + w / 2,
                cy,
                line,
                size=size,
                weight=weight,
                mono=mono,
                fill=text_color,
            )
            cy += lh
        cy += (sh - lh) * 0.2
        for line in subs:
            self.text(x + w / 2, cy, line, size=sub_size, fill=P.INK_2, mono=mono)
            cy += sh

    def text(
        self,
        x,
        y,
        s,
        size=P.SIZE_LABEL,
        anchor="middle",
        weight=None,
        fill=P.INK,
        mono=False,
        italic=False,
        rotate=None,
        opacity=None,
    ):
        transform = f"rotate({fmt(rotate)} {fmt(x)} {fmt(y)})" if rotate else None
        self._parts.append(
            "<text"
            + _attrs(
                [
                    ("x", fmt(x)),
                    ("y", fmt(y)),
                    ("font-family", P.MONO if mono else P.FONT),
                    ("font-size", fmt(size)),
                    ("font-weight", weight),
                    ("font-style", "italic" if italic else None),
                    ("text-anchor", anchor if anchor != "start" else None),
                    ("fill", fill),
                    ("opacity", fmt(opacity) if opacity else None),
                    ("transform", transform),
                ]
            )
            + f">{escape(str(s))}</text>"
        )

    def label(self, x, y, s, **kw):
        """Small secondary-ink annotation."""
        kw.setdefault("size", P.SIZE_SMALL)
        kw.setdefault("fill", P.INK_2)
        self.text(x, y, s, **kw)

    def line(
        self,
        x1,
        y1,
        x2,
        y2,
        stroke=P.INK_2,
        sw=P.STROKE_HAIR,
        dash=None,
        opacity=None,
        cap=None,
    ):
        self._parts.append(
            "<line"
            + _attrs(
                [
                    ("x1", fmt(x1)),
                    ("y1", fmt(y1)),
                    ("x2", fmt(x2)),
                    ("y2", fmt(y2)),
                    ("stroke", stroke),
                    ("stroke-width", fmt(sw)),
                    ("stroke-dasharray", dash),
                    ("stroke-linecap", cap),
                    ("opacity", fmt(opacity) if opacity else None),
                ]
            )
            + "/>"
        )

    def arrow(
        self,
        x1,
        y1,
        x2,
        y2,
        stroke=P.INK_2,
        sw=P.STROKE_ARROW,
        dash=None,
        both=False,
        opacity=None,
    ):
        mid = self._marker_id(stroke)
        self._parts.append(
            "<line"
            + _attrs(
                [
                    ("x1", fmt(x1)),
                    ("y1", fmt(y1)),
                    ("x2", fmt(x2)),
                    ("y2", fmt(y2)),
                    ("stroke", stroke),
                    ("stroke-width", fmt(sw)),
                    ("stroke-dasharray", dash),
                    ("marker-end", f"url(#{mid})"),
                    ("marker-start", f"url(#{mid})" if both else None),
                    ("opacity", fmt(opacity) if opacity else None),
                ]
            )
            + "/>"
        )

    def path(
        self,
        d,
        stroke=P.INK_2,
        sw=P.STROKE_ARROW,
        fill="none",
        dash=None,
        arrow=False,
        opacity=None,
        fill_opacity=None,
    ):
        mid = self._marker_id(stroke) if arrow else None
        self._parts.append(
            "<path"
            + _attrs(
                [
                    ("d", d),
                    ("fill", fill),
                    ("fill-opacity", fmt(fill_opacity) if fill_opacity else None),
                    ("stroke", stroke),
                    ("stroke-width", fmt(sw) if stroke != "none" else None),
                    ("stroke-dasharray", dash),
                    ("stroke-linejoin", "round" if stroke != "none" else None),
                    ("marker-end", f"url(#{mid})" if mid else None),
                    ("opacity", fmt(opacity) if opacity else None),
                ]
            )
            + "/>"
        )

    def elbow(
        self,
        x1,
        y1,
        x2,
        y2,
        stroke=P.INK_2,
        sw=P.STROKE_ARROW,
        dash=None,
        arrow=True,
        via="v",
    ):
        """An orthogonal connector: vertical-then-horizontal (``via="v"``) or the reverse."""
        if via == "v":
            d = f"M {fmt(x1)} {fmt(y1)} V {fmt(y2)} H {fmt(x2)}"
        elif via == "h":
            d = f"M {fmt(x1)} {fmt(y1)} H {fmt(x2)} V {fmt(y2)}"
        else:
            ym = (y1 + y2) / 2
            d = f"M {fmt(x1)} {fmt(y1)} V {fmt(ym)} H {fmt(x2)} V {fmt(y2)}"
        self.path(d, stroke=stroke, sw=sw, dash=dash, arrow=arrow)

    def polygon(
        self,
        pts,
        fill="none",
        stroke=P.INK,
        sw=P.STROKE_BOX,
        fill_opacity=None,
        dash=None,
        opacity=None,
    ):
        d = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in pts)
        self._parts.append(
            "<polygon"
            + _attrs(
                [
                    ("points", d),
                    ("fill", fill),
                    ("fill-opacity", fmt(fill_opacity) if fill_opacity else None),
                    ("stroke", stroke),
                    ("stroke-width", fmt(sw) if stroke != "none" else None),
                    ("stroke-linejoin", "round"),
                    ("stroke-dasharray", dash),
                    ("opacity", fmt(opacity) if opacity else None),
                ]
            )
            + "/>"
        )

    def polyline(
        self, pts, stroke=P.INK, sw=P.STROKE_BOX, dash=None, arrow=False, opacity=None
    ):
        d = " ".join(f"{fmt(x)},{fmt(y)}" for x, y in pts)
        mid = self._marker_id(stroke) if arrow else None
        self._parts.append(
            "<polyline"
            + _attrs(
                [
                    ("points", d),
                    ("fill", "none"),
                    ("stroke", stroke),
                    ("stroke-width", fmt(sw)),
                    ("stroke-linejoin", "round"),
                    ("stroke-dasharray", dash),
                    ("marker-end", f"url(#{mid})" if mid else None),
                    ("opacity", fmt(opacity) if opacity else None),
                ]
            )
            + "/>"
        )

    def circle(
        self,
        cx,
        cy,
        r,
        fill=P.INK,
        stroke=None,
        sw=P.STROKE_BOX,
        opacity=None,
        fill_opacity=None,
    ):
        self._parts.append(
            "<circle"
            + _attrs(
                [
                    ("cx", fmt(cx)),
                    ("cy", fmt(cy)),
                    ("r", fmt(r)),
                    ("fill", fill),
                    ("fill-opacity", fmt(fill_opacity) if fill_opacity else None),
                    ("stroke", stroke),
                    ("stroke-width", fmt(sw) if stroke else None),
                    ("opacity", fmt(opacity) if opacity else None),
                ]
            )
            + "/>"
        )

    def node(
        self,
        cx,
        cy,
        index,
        kind="corner",
        r=4.5,
        label_dx=7,
        label_dy=-6,
        label=True,
        anchor="start",
    ):
        """A numbered mesh node: corners solid ink, higher-order nodes tinted."""
        color = {"corner": P.INK, "edge": P.CORE, "face": P.PYTHON, "body": P.CABI}[
            kind
        ]
        if kind == "corner":
            self.circle(cx, cy, r, fill=color)
        else:
            self.circle(cx, cy, r, fill=P.WHITE, stroke=color, sw=2)
        if label:
            self.text(
                cx + label_dx,
                cy + label_dy,
                str(index),
                size=P.SIZE_SMALL,
                anchor=anchor,
                fill=color if kind != "corner" else P.INK,
                weight="600",
            )

    def caption(self, x, y, s, anchor="middle", size=P.SIZE_LABEL):
        self.text(x, y, s, size=size, anchor=anchor, weight="600")

    def legend(self, x, y, items, swatch=12, gap=18, size=P.SIZE_SMALL, mono=False):
        """Horizontal legend: ``items`` are ``(colour, label)`` pairs."""
        cx = x
        for color, label in items:
            self.rrect(
                cx,
                y - swatch + 2,
                swatch,
                swatch,
                fill=color,
                fill_opacity=P.FILL_OPACITY * 3,
                stroke=color,
                rx=2,
                sw=1,
            )
            self.text(
                cx + swatch + 5,
                y,
                label,
                size=size,
                anchor="start",
                fill=P.INK_2,
                mono=mono,
            )
            cx += swatch + 5 + len(label) * size * 0.56 + gap

    # -- output --------------------------------------------------------------
    def render(self):
        defs = "".join(self._markers[k] for k in sorted(self._markers))
        head = (
            '<svg xmlns="http://www.w3.org/2000/svg"'
            f' viewBox="0 0 {fmt(self.w)} {fmt(self.h)}"'
            f' width="{fmt(self.w)}" height="{fmt(self.h)}"'
            f' font-family="{P.FONT}" role="img" aria-label="{escape(self.title)}">\n'
            "<!-- generated by doc/diagrams/gen_diagrams.py; do not edit by hand -->\n"
            f"<title>{escape(self.title)}</title>\n"
        )
        body = "\n".join(self._parts)
        return head + (f"<defs>{defs}</defs>\n" if defs else "") + body + "\n</svg>\n"


def text_width(s, size=P.SIZE_LABEL, mono=False):
    """A rough width estimate for laying out labels (no font metrics needed)."""
    return len(str(s)) * size * (0.6 if mono else 0.55)
