import os
import pathlib

from .._common import error
from .._helpers import _filetypes_from_path, read
from ._common import add_io_args


def add_args(parser):
    add_io_args(parser, in_help="mesh file to decompress")


def decompress(args):
    # Imported when the verb runs, so building the parser loads no format.
    from .. import cgns, h5m, pcd, vtu, xdmf

    if args.input_format:
        fmts = [args.input_format]
    else:
        fmts = _filetypes_from_path(pathlib.Path(args.infile))
    # pick the first
    fmt = fmts[0]

    size = os.stat(args.infile).st_size
    print(f"File size before: {size / 1024 ** 2:.2f} MB")
    mesh = read(args.infile, file_format=args.input_format)

    # # Some converters (like VTK) require `points` to be contiguous.
    # mesh.points = np.ascontiguousarray(mesh.points)

    # write it out
    if fmt == "cgns":
        cgns.write(args.infile, mesh, compression=None)
    elif fmt == "h5m":
        h5m.write(args.infile, mesh, compression=None)
    elif fmt == "pcd":
        pcd.write(args.infile, mesh, data="binary", point_dtype="keep")
    elif fmt == "vtu":
        vtu.write(args.infile, mesh, binary=True, compression=None)
    elif fmt == "xdmf":
        xdmf.write(args.infile, mesh, data_format="HDF", compression=None)
    else:
        error(f"Don't know how to decompress {args.infile}.")
        exit(1)

    size = os.stat(args.infile).st_size
    print(f"File size after:  {size / 1024 ** 2:.2f} MB")
