"""``grid-sample`` / ``grid-scatter`` / ``grid-resample`` / ``grid-spectrum``.

The mesh-to-grid data path a convolutional or superresolution model needs. Write
grids as ``.vti``: it stores the lattice as origin/spacing/extent, so the spec is
recovered exactly on the way back in.
"""

import numpy as np

from .._grid_transfer import (
    GridArray,
    GridSpec,
    power_spectrum,
    resample_grid,
    sample_grid,
    scatter_grid,
)
from .._helpers import _writer_map, read, reader_map, write


def _fmt_args(parser, output=True):
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    if output:
        parser.add_argument(
            "--output-format",
            "-o",
            type=str,
            choices=sorted(list(_writer_map.keys())),
            help="output file format",
            default=None,
        )


def _triple(text, flag, verb):
    values = [int(x) for x in text.split(",")]
    if len(values) != 3:
        raise ValueError(f"{verb}: {flag} expects 'nx,ny,nz'")
    return values


def _report_spec(spec):
    d, s = spec.dims, spec.spacing
    print(f"  grid:           {int(d[0])} x {int(d[1])} x {int(d[2])} cells")
    print(f"  shape (D,H,W):  {spec.shape}")
    print(f"  cell size:      {s[0]:g}, {s[1]:g}, {s[2]:g}")


# --------------------------------------------------------------------------- #
# grid-sample                                                                 #
# --------------------------------------------------------------------------- #
def add_sample_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument("outfile", type=str, help="grid to be written to (.vti)")
    _fmt_args(parser)
    parser.add_argument(
        "--resolution",
        type=str,
        default=None,
        help="cell counts 'nx,ny,nz' (exactly one of --resolution/--cell-size)",
    )
    parser.add_argument(
        "--cell-size",
        type=float,
        default=None,
        help="cubic cell size (exactly one of --resolution/--cell-size)",
    )
    parser.add_argument(
        "--bounds",
        type=str,
        default=None,
        help="'xlo,ylo,zlo,xhi,yhi,zhi'; the mesh's bounding box by default. "
        "For a coarse/fine pair take the box from the FINE mesh, so every fine "
        "node is inside the coarse grid (negatives need --bounds=)",
    )
    parser.add_argument(
        "--padding", type=float, default=0.0, help="grow the box by this on every side"
    )
    parser.add_argument(
        "--padding-relative",
        type=float,
        default=0.0,
        help="grow the box by this fraction of its diagonal on every side",
    )
    parser.add_argument(
        "--fields",
        type=str,
        default=None,
        help="comma-separated point_data names; every one by default",
    )
    parser.add_argument(
        "--extrapolate",
        action="store_true",
        help="give a grid point outside the mesh its nearest value instead of the fill",
    )
    parser.add_argument(
        "--fill-value",
        type=float,
        default=0.0,
        help="what an outside point gets (negatives need --fill-value=; nan makes "
        "the fill visible downstream)",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=20000000,
        help="refuse above this many cells (default ~256^3)",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def grid_sample_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    resolution = (
        None
        if args.resolution is None
        else _triple(args.resolution, "--resolution", "grid-sample")
    )
    bounds = None
    if args.bounds is not None:
        bounds = [float(x) for x in args.bounds.split(",")]
        if len(bounds) != 6:
            raise ValueError("grid-sample: --bounds expects 'xlo,ylo,zlo,xhi,yhi,zhi'")
    fields = None if args.fields is None else args.fields.split(",")

    spec = GridSpec.from_mesh(
        mesh,
        resolution=resolution,
        cell_size=args.cell_size,
        bounds=bounds,
        padding=args.padding,
        padding_relative=args.padding_relative,
        max_cells=args.max_cells,
    )
    array = sample_grid(
        mesh,
        spec,
        fields=fields,
        extrapolate=args.extrapolate,
        fill_value=args.fill_value,
    )
    if not args.quiet:
        print("sampled onto a grid")
        _report_spec(spec)
        print(f"  channels:       {', '.join(array.channels)}")
        if array.coverage is not None:
            print(f"  coverage:       {array.coverage:.4f} of grid points inside")
    write(args.outfile, array.to_mesh(), file_format=args.output_format)
    return 0


# --------------------------------------------------------------------------- #
# grid-scatter                                                                #
# --------------------------------------------------------------------------- #
def add_scatter_args(parser):
    parser.add_argument("gridfile", type=str, help="grid to read (.vti)")
    parser.add_argument("targetfile", type=str, help="mesh to write the fields onto")
    parser.add_argument("outfile", type=str, help="mesh to be written to")
    parser.add_argument(
        "--grid-format", type=str, default=None, help="grid file format"
    )
    parser.add_argument(
        "--target-format", type=str, default=None, help="target file format"
    )
    parser.add_argument(
        "--output-format",
        "-o",
        type=str,
        choices=sorted(list(_writer_map.keys())),
        help="output file format",
        default=None,
    )
    parser.add_argument(
        "--fields",
        type=str,
        default=None,
        help="comma-separated grid point_data names; every one by default",
    )
    parser.add_argument(
        "--on-conflict",
        type=str,
        default="error",
        choices=["error", "overwrite", "suffix"],
        help="what to do when the target already has an array of that name",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def grid_scatter_cmd(args):
    fields = None if args.fields is None else args.fields.split(",")
    array = GridArray.from_mesh(
        read(args.gridfile, file_format=args.grid_format), fields=fields
    )
    target = read(args.targetfile, file_format=args.target_format)
    out = scatter_grid(array, target, on_conflict=args.on_conflict)
    if not args.quiet:
        print("scattered a grid onto a mesh")
        _report_spec(array.spec)
        print(f"  channels:       {', '.join(array.channels)}")
        print(f"  target points:  {len(target.points)}")
    write(args.outfile, out, file_format=args.output_format)
    return 0


# --------------------------------------------------------------------------- #
# grid-resample                                                               #
# --------------------------------------------------------------------------- #
def add_resample_args(parser):
    parser.add_argument("infile", type=str, help="grid to read (.vti)")
    parser.add_argument("outfile", type=str, help="grid to be written to (.vti)")
    _fmt_args(parser)
    parser.add_argument(
        "--factor",
        type=int,
        default=None,
        help="integer upscale, same box (exactly one of --factor/--resolution)",
    )
    parser.add_argument(
        "--resolution",
        type=str,
        default=None,
        help="target cell counts 'nx,ny,nz' over the same box",
    )
    parser.add_argument(
        "--fields",
        type=str,
        default=None,
        help="comma-separated point_data names; every one by default",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", help="suppress the summary"
    )


def grid_resample_cmd(args):
    fields = None if args.fields is None else args.fields.split(",")
    array = GridArray.from_mesh(
        read(args.infile, file_format=args.input_format), fields=fields
    )
    if (args.factor is None) == (args.resolution is None):
        raise ValueError("grid-resample: give exactly one of --factor and --resolution")
    if args.factor is not None:
        target = array.spec.upscale(args.factor)
    else:
        dims = np.asarray(
            _triple(args.resolution, "--resolution", "grid-resample"), dtype=np.int64
        )
        lo, hi = array.spec.bounds
        target = GridSpec(origin=lo, spacing=(hi - lo) / dims.astype(float), dims=dims)
    values = resample_grid(array.values, array.spec, target)
    out = GridArray(values, target, array.channels, dict(array.schema))
    if not args.quiet:
        print("resampled a grid")
        _report_spec(target)
        factor = array.spec.scaling_factor(target)
        print(f"  factor:         {factor if factor else 'not a whole number'}")
    write(args.outfile, out.to_mesh(), file_format=args.output_format)
    return 0


# --------------------------------------------------------------------------- #
# grid-spectrum                                                               #
# --------------------------------------------------------------------------- #
def add_spectrum_args(parser):
    parser.add_argument("infile", type=str, help="grid to read (.vti)")
    _fmt_args(parser, output=False)
    parser.add_argument("--field", type=str, required=True, help="point_data name")
    parser.add_argument(
        "--max-bins", type=int, default=0, help="print at most this many bins (0 = all)"
    )
    parser.add_argument("--json", action="store_true", help="emit JSON")


def grid_spectrum_cmd(args):
    array = GridArray.from_mesh(read(args.infile, file_format=args.input_format))
    names = [
        c for c in array.channels if c == args.field or c.startswith(args.field + "_")
    ]
    if not names:
        raise ValueError(
            f"grid-spectrum: no channel named {args.field!r} "
            f"(have {list(array.channels)})"
        )
    values = np.stack([array.channel(n) for n in names], axis=0)
    ps = power_spectrum(values, array.spec)
    keep = min(args.max_bins, len(ps.power)) if args.max_bins > 0 else len(ps.power)

    if args.json:
        import json

        print(
            json.dumps(
                {
                    "field": args.field,
                    "channels": names,
                    "units": ps.units,
                    "wavenumber": ps.wavenumber[:keep].tolist(),
                    "power": ps.power[:keep].tolist(),
                    "counts": ps.counts[:keep].tolist(),
                    "num_bins": int(len(ps.power)),
                    "total_power": float(ps.power.sum()),
                },
                indent=2,
            )
        )
        return 0

    print(f"power spectrum of {args.field} ({', '.join(names)})")
    print(f"  bins:           {len(ps.power)}   units: {ps.units}")
    print(f"  total power:    {float(ps.power.sum()):.6g}  (= mean of the square)")
    print(f"  {'wavenumber':>14}  {'power':>14}  {'modes':>7}")
    for k, pw, n in zip(ps.wavenumber[:keep], ps.power[:keep], ps.counts[:keep]):
        print(f"  {k:14.6g}  {pw:14.6g}  {int(n):7d}")
    return 0
