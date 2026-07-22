#!/usr/bin/env python3
"""Regenerate the committed doc/README figures from the bundled sample meshes.

One command, two kinds of output:

* **Vector** (``doc/public/images/*.svg``) -- the SVG writer's data-driven
  colouring. These are the canonical figures for anything that is really a
  picture *of a mesh coloured by a field*: they are dependency-free, they scale,
  and they are produced by exactly the code path they document.
* **Raster** (``doc/public/viewer/*.png``) -- ``screenshot()``, i.e. the
  Polyscope ``[viewer]`` extra, for the shaded 3D views a flat vector
  projection cannot convey.

This *complements* the PyVista notebook convention rather than replacing it:
``example/*.ipynb`` still owns the executable, re-run-with-outputs
demonstrations. What lives here is the small set of static figures the docs and
README embed, so they can be rebuilt reproducibly instead of by hand.

Usage::

    python tools/gen_doc_images.py            # everything available
    python tools/gen_doc_images.py --list     # just say what would be written

The raster half is skipped with a warning when polyscope is not installed
(``pip install meshioplusplus[viewer]``); the vector half needs nothing beyond
meshio++ itself.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import meshioplusplus  # noqa: E402

SAMPLE = REPO / "example" / "example.msh"
IMAGES = REPO / "doc" / "public" / "images"
VIEWER = REPO / "doc" / "public" / "viewer"

# Rendered small: these are figures, not full-resolution plots, and a smaller
# viewBox keeps the committed SVGs a sane size.
IMAGE_WIDTH = 800.0

# The bundled bracket has ~298k cells, whose skin is ~58k facets -- a 4.5 MB SVG
# that is also illegible at figure size, since every element is sub-pixel. The
# figures therefore use a corner of it: real bracket geometry, ~1.1k facets, and
# individual elements you can actually see. (Decimating the whole skin instead
# was tried when `decimate` landed: feature pinning floors the collapse at
# ~6.6k facets even at feature_angle=60 -- the bracket is mostly creases and
# fillets at this resolution -- which is still a ~700 kB SVG, so the crop
# stays.)
CROP_FRACTION = 0.40


def _display(path: pathlib.Path) -> str:
    """Repo-relative when it can be (the tests redirect output elsewhere)."""
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def _bracket():
    """A legible corner of the bundled Gmsh bracket, the repo's demo mesh."""
    full = meshioplusplus.read(SAMPLE)
    # Drop the sets before cropping: crop's set remapping raises on this file
    # (an unrelated pre-existing bug), and a figure has no use for them anyway.
    bare = meshioplusplus.Mesh(full.points, full.cells)
    lo = bare.points.min(axis=0)
    hi = bare.points.max(axis=0)
    cropped = meshioplusplus.crop(bare, bbox=[*lo, *(lo + (hi - lo) * CROP_FRACTION)])
    # crop keeps the input's block structure 1:1, so blocks fully outside the
    # box survive as empty ones. Drop them: they contribute nothing to a figure,
    # and _viewer.py raises on the zero-length cell_data array that attaching a
    # per-cell quantity to an empty block produces (an unrelated pre-existing
    # bug, hit here only because this is the rare mesh with empty blocks).
    return meshioplusplus.Mesh(
        cropped.points,
        [b for b in cropped.cells if len(b.data) > 0],
    )


def vector_figures(dry_run: bool) -> list[pathlib.Path]:
    """The colour-by-a-field SVGs."""
    written = []

    def emit(name, mesh, **kwargs):
        path = IMAGES / name
        written.append(path)
        if dry_run:
            return
        meshioplusplus.svg.write(
            path, mesh, image_width=IMAGE_WIDTH, colorbar=True, **kwargs
        )
        print(f"  wrote {_display(path)}")

    mesh = _bracket()

    # 1. Element quality. attach_quality writes one cell_data array per metric,
    #    and the projected skin picks up each facet's owning cell's value.
    quality = meshioplusplus.attach_quality(mesh)
    emit(
        "color_by_quality.svg",
        quality,
        color_by="quality:scaled_jacobian",
        cmap="viridis",
    )

    # 2. Domain decomposition. partition_labels attaches the "partition:part"
    #    contract array, which is exactly a per-cell integer to colour by.
    labels = meshioplusplus.partition_labels(mesh, 6)
    tagged = meshioplusplus.Mesh(
        mesh.points,
        mesh.cells,
        cell_data={"partition:part": [np.asarray(b) for b in labels]},
    )
    emit("color_by_partition.svg", tagged, color_by="partition:part", cmap="turbo")

    # 3. Per-point error against a perturbed copy -- the kind of field a
    #    regression diff produces. coolwarm because the quantity is signed.
    rng = np.random.default_rng(0)
    moved = meshioplusplus.Mesh(
        mesh.points + rng.normal(scale=0.002, size=mesh.points.shape),
        mesh.cells,
    )
    error = np.linalg.norm(moved.points - mesh.points, axis=1)
    with_error = meshioplusplus.Mesh(
        mesh.points, mesh.cells, point_data={"error": error}
    )
    emit("color_by_diff_error.svg", with_error, color_by="error", cmap="coolwarm")

    return written


def raster_figures(dry_run: bool) -> list[pathlib.Path]:
    """The shaded Polyscope screenshots."""
    path = VIEWER / "desktop-viewer.png"
    if dry_run:
        return [path]
    try:
        import polyscope  # noqa: F401
    except ImportError:
        print(
            "  skipping the raster figures: polyscope is not installed "
            "(pip install meshioplusplus[viewer])",
            file=sys.stderr,
        )
        return []
    quality = meshioplusplus.attach_quality(_bracket())
    meshioplusplus.screenshot(
        quality, path, color_by="quality:scaled_jacobian", size=(1600, 1200)
    )
    print(f"  wrote {_display(path)}")
    return [path]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--list",
        action="store_true",
        help="print the files that would be written and exit",
    )
    parser.add_argument(
        "--vector-only",
        action="store_true",
        help="skip the polyscope raster figures",
    )
    args = parser.parse_args()

    if not SAMPLE.exists():
        print(f"missing sample mesh {SAMPLE}", file=sys.stderr)
        return 1

    if not args.list:
        IMAGES.mkdir(parents=True, exist_ok=True)
        VIEWER.mkdir(parents=True, exist_ok=True)
        print("regenerating doc figures from", SAMPLE.relative_to(REPO))

    written = vector_figures(args.list)
    if not args.vector_only:
        written += raster_figures(args.list)

    if args.list:
        for path in written:
            print(_display(path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
