import os
import pathlib

from .._common import error
from .._helpers import _filetypes_from_path, read
from ._common import add_io_args


def add_args(parser):
    add_io_args(parser, in_help="mesh file to convert")


def ascii(args):
    # Imported when the verb runs, so building the parser loads no format.
    from .. import (
        ansys,
        flac3d,
        frd,
        gmsh,
        mdpa,
        openfoam,
        pcd,
        ply,
        stl,
        vtk,
        vtu,
        xdmf,
    )

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
    if fmt == "ansys":
        ansys.write(args.infile, mesh, binary=False)
    elif fmt == "flac3d":
        flac3d.write(args.infile, mesh, binary=False)
    elif fmt == "frd":
        frd.write(args.infile, mesh, binary=False)
    elif fmt == "gmsh":
        gmsh.write(args.infile, mesh, binary=False)
    elif fmt == "mdpa":
        mdpa.write(args.infile, mesh, binary=False)
    elif fmt == "openfoam":
        openfoam.write(args.infile, mesh, binary=False)
    elif fmt == "pcd":
        pcd.write(args.infile, mesh, binary=False, point_dtype="keep")
    elif fmt == "ply":
        ply.write(args.infile, mesh, binary=False)
    elif fmt == "stl":
        stl.write(args.infile, mesh, binary=False)
    elif fmt == "vtk":
        vtk.write(args.infile, mesh, binary=False)
    elif fmt == "vtu":
        vtu.write(args.infile, mesh, binary=False)
    elif fmt == "xdmf":
        xdmf.write(args.infile, mesh, data_format="XML")
    else:
        error(f"Don't know how to convert {args.infile} to ASCII format.")
        return 1

    size = os.stat(args.infile).st_size
    print(f"File size after: {size / 1024 ** 2:.2f} MB")
    return 0
