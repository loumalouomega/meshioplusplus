#!/usr/bin/env python3
"""Derive small square/landscape raster assets from the compiled logo PDFs.

Run after `logo.pdf`/`logo-icon.pdf` exist (build.sh calls this once the two
`compile_one` steps finish). Produces, all committed:

- ``logo-icon-square.png`` -- the bunny icon padded (transparent) into a
  square canvas, for small inline use (README badge, favicon source).
- ``favicon.ico`` -- multi-resolution (16/32/48/64 px) favicon derived from
  the square icon, copied into ``doc/public/`` for the docs site.
- ``social-preview.png`` -- the full banner (icon + wordmark) centered on a
  white 1280x640 canvas, sized for GitHub's repo social-preview image
  (Settings -> General -> Social preview; GitHub has no API/config file for
  this, it must be uploaded by hand).

Requires PyMuPDF (fitz) + Pillow, both already used by build.sh's PNG step.
"""

from __future__ import annotations

import os

import fitz  # PyMuPDF
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
DOC_PUBLIC = os.path.join(HERE, "..", "public")

ICON_MASTER_SIZE = 512  # px, square, transparent
FAVICON_SIZES = [16, 32, 48, 64]
SOCIAL_PREVIEW_SIZE = (1280, 640)  # GitHub's ideal social-preview dimensions
# GitHub's own "repository-open-graph-template.png" (Settings -> General ->
# Social preview -> "learn more") draws a safe-area guide ~80px in from each
# edge of the 1280x640 canvas ("leave a 40pt border ... so nothing gets
# cropped"); measured directly off that template, kept here as the
# authoritative margin so the banner is sized to fill the safe area rather
# than an arbitrary fraction of the canvas.
SOCIAL_PREVIEW_SAFE_MARGIN = 80


def render_pdf_page(path: str, target_long_edge: int) -> Image.Image:
    """Rasterize a single-page PDF's page to an RGBA PIL image, transparent
    background, scaled so its longer edge is ``target_long_edge`` pixels."""
    doc = fitz.open(path)
    page = doc[0]
    scale = target_long_edge / max(page.rect.width, page.rect.height)
    pix = page.get_pixmap(matrix=fitz.Matrix(scale, scale), alpha=True)
    return Image.frombytes("RGBA", (pix.width, pix.height), pix.samples)


def paste_centered(canvas: Image.Image, art: Image.Image) -> None:
    x = (canvas.width - art.width) // 2
    y = (canvas.height - art.height) // 2
    canvas.paste(art, (x, y), art)


def make_square_icon() -> Image.Image:
    art = render_pdf_page(os.path.join(HERE, "logo-icon.pdf"), ICON_MASTER_SIZE)
    side = max(art.width, art.height)
    square = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    paste_centered(square, art)
    return square


def make_social_preview() -> Image.Image:
    canvas = Image.new("RGBA", SOCIAL_PREVIEW_SIZE, (255, 255, 255, 255))
    # Fill the safe area GitHub's template guides against, so the banner
    # doesn't touch the crop edges but still uses the available space.
    safe_width = SOCIAL_PREVIEW_SIZE[0] - 2 * SOCIAL_PREVIEW_SAFE_MARGIN
    safe_height = SOCIAL_PREVIEW_SIZE[1] - 2 * SOCIAL_PREVIEW_SAFE_MARGIN
    art = render_pdf_page(os.path.join(HERE, "logo.pdf"), safe_width)
    if art.height > safe_height:
        scale = safe_height / art.height
        art = art.resize(
            (int(art.width * scale), int(art.height * scale)), Image.LANCZOS
        )
    paste_centered(canvas, art)
    return canvas.convert("RGB")


def main() -> None:
    square = make_square_icon()
    square_path = os.path.join(HERE, "logo-icon-square.png")
    square.save(square_path)
    print(f"wrote {square_path}: {square.size}")

    favicon_frames = [square.resize((s, s), Image.LANCZOS) for s in FAVICON_SIZES]
    favicon_path = os.path.join(DOC_PUBLIC, "favicon.ico")
    favicon_frames[-1].save(
        favicon_path, format="ICO", sizes=[(s, s) for s in FAVICON_SIZES]
    )
    print(f"wrote {favicon_path}: sizes {FAVICON_SIZES}")

    preview = make_social_preview()
    preview_path = os.path.join(HERE, "social-preview.png")
    preview.save(preview_path)
    print(f"wrote {preview_path}: {preview.size}")


if __name__ == "__main__":
    main()
