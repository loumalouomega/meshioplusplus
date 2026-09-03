#!/usr/bin/env python3
"""Generate the documentation diagrams under ``doc/public/diagrams/``.

Usage::

    python3 doc/diagrams/gen_diagrams.py            # write every SVG + PNG twin
    python3 doc/diagrams/gen_diagrams.py --check    # CI gate: stale/missing/unused?
    python3 doc/diagrams/gen_diagrams.py --list     # names only, touches nothing
    python3 doc/diagrams/gen_diagrams.py --only architecture --no-png

Standard library only. PNG twins are rasterised at 2x with ``rsvg-convert``
(falling back to ``inkscape``, then ImageMagick); their bytes are not pinned
because rasterisers differ per machine, only their existence and size.
"""

import argparse
import pathlib
import re
import shutil
import struct
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent
OUT = REPO / "doc" / "public" / "diagrams"
DOC = REPO / "doc"
README = REPO / "README.md"
SCALE = 2

sys.path.insert(0, str(HERE))
from diaglib import palette  # noqa: E402
from figures import REGISTRY  # noqa: E402


def render_all(names=None):
    """Return ``{name: svg_text}`` for the requested (default: all) figures."""
    names = list(REGISTRY) if not names else list(names)
    unknown = [n for n in names if n not in REGISTRY]
    if unknown:
        raise SystemExit(f"unknown diagram(s): {', '.join(unknown)}")
    return {n: REGISTRY[n]() for n in names}


def viewbox_size(svg_text):
    m = re.search(r'viewBox="0 0 (\d+(?:\.\d+)?) (\d+(?:\.\d+)?)"', svg_text)
    return float(m.group(1)), float(m.group(2))


def png_size(path):
    data = path.read_bytes()[:24]
    if len(data) < 24 or data[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return struct.unpack(">II", data[16:24])


def write_svgs(out, names=None):
    out.mkdir(parents=True, exist_ok=True)
    written = []
    for name, text in render_all(names).items():
        path = out / f"{name}.svg"
        path.write_text(text, encoding="utf-8", newline="\n")
        written.append(path)
    return written


def rasterise(svg_path, tool="auto", scale=SCALE):
    png_path = svg_path.with_suffix(".png")
    w, h = viewbox_size(svg_path.read_text(encoding="utf-8"))
    width = int(round(w * scale))
    have = {
        t: shutil.which(t) for t in ("rsvg-convert", "inkscape", "magick", "convert")
    }
    order = {
        "auto": ["rsvg-convert", "inkscape", "magick", "convert"],
        "rsvg": ["rsvg-convert"],
        "inkscape": ["inkscape"],
        "magick": ["magick", "convert"],
    }[tool]
    for t in order:
        if not have[t]:
            continue
        if t == "rsvg-convert":
            cmd = [have[t], "--zoom", str(scale), "-o", str(png_path), str(svg_path)]
        elif t == "inkscape":
            cmd = [
                have[t],
                "--export-type=png",
                f"--export-width={width}",
                f"--export-filename={png_path}",
                str(svg_path),
            ]
        else:
            cmd = [
                have[t],
                "-density",
                str(96 * scale),
                "-background",
                "none",
                str(svg_path),
                str(png_path),
            ]
        subprocess.run(
            cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        return png_path
    print(
        f"warning: no rasteriser found for {svg_path.name} (install librsvg, inkscape or ImageMagick)",
        file=sys.stderr,
    )
    return None


def referenced_names():
    """Diagram names embedded somewhere under doc/ or in README.md."""
    pattern = re.compile(r"/diagrams/([A-Za-z0-9_]+)\.(?:svg|png)")
    found = set()
    for md in list(DOC.rglob("*.md")) + [README]:
        if "node_modules" in md.parts or ".vitepress" in md.parts:
            continue
        found.update(pattern.findall(md.read_text(encoding="utf-8")))
    return found


def check(out):
    problems = []
    rendered = render_all()
    committed = {p.stem: p for p in out.glob("*.svg")} if out.exists() else {}
    for name, text in rendered.items():
        path = committed.get(name)
        if path is None:
            problems.append(f"missing: {name}.svg (run gen_diagrams.py)")
            continue
        if path.read_text(encoding="utf-8") != text:
            problems.append(
                f"stale: {name}.svg differs from the generator (run gen_diagrams.py)"
            )
        png = path.with_suffix(".png")
        size = png_size(png) if png.exists() else None
        if size is None:
            problems.append(
                f"missing or invalid PNG twin: {png.name} (run gen_diagrams.py)"
            )
        else:
            w, h = viewbox_size(text)
            if size != (int(round(w * SCALE)), int(round(h * SCALE))):
                problems.append(
                    f"PNG twin has the wrong size: {png.name} is {size[0]}x{size[1]}, expected {int(round(w * SCALE))}x{int(round(h * SCALE))}"
                )
        for colour in set(re.findall(r"#[0-9a-fA-F]{6}\b", text)):
            if colour.lower() not in palette.ALL:
                problems.append(
                    f"{name}.svg uses a colour outside diaglib/palette.py: {colour}"
                )
    for name in sorted(committed):
        if name not in rendered:
            problems.append(
                f"stray: {name}.svg is not produced by any registered figure"
            )
    if out.exists():
        for p in sorted(out.iterdir()):
            if p.suffix not in (".svg", ".png") or (
                p.suffix == ".png" and p.stem not in rendered
            ):
                problems.append(f"stray file in {out.relative_to(REPO)}: {p.name}")
    used = referenced_names()
    for name in rendered:
        if name not in used:
            problems.append(
                f"unused: {name} is embedded by no page under doc/ nor by README.md"
            )
    for name in sorted(used - set(rendered)):
        problems.append(
            f"dangling: a page references /diagrams/{name} but no figure produces it"
        )
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--list", action="store_true", help="print the figure names and exit"
    )
    ap.add_argument(
        "--check",
        action="store_true",
        help="verify the committed outputs are current, twinned and used",
    )
    ap.add_argument("--only", nargs="+", metavar="NAME", help="only these figures")
    ap.add_argument(
        "--out",
        type=pathlib.Path,
        default=OUT,
        help=f"output directory (default {OUT.relative_to(REPO)})",
    )
    ap.add_argument("--no-png", action="store_true", help="skip the PNG twins")
    ap.add_argument(
        "--png-tool", choices=["auto", "rsvg", "inkscape", "magick"], default="auto"
    )
    args = ap.parse_args(argv)

    if args.list:
        for name in REGISTRY:
            print(name)
        return 0
    if args.check:
        problems = check(args.out)
        if problems:
            print(f"doc/diagrams: {len(problems)} problem(s)", file=sys.stderr)
            for p in problems:
                print(f"  - {p}", file=sys.stderr)
            return 1
        print(f"doc/diagrams: {len(REGISTRY)} figures current, twinned and referenced")
        return 0
    written = write_svgs(args.out, args.only)
    for path in written:
        print(f"wrote {path.relative_to(REPO) if path.is_relative_to(REPO) else path}")
        if not args.no_png:
            rasterise(path, args.png_tool)
    return 0


if __name__ == "__main__":
    sys.exit(main())
