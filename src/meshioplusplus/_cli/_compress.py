import os
import pathlib

from .. import ansys, cgns, gmsh, h5m, mdpa, ply, stl, vtk, vtp, vtu, xdmf
from .._common import error
from .._helpers import _filetypes_from_path, read, reader_map


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to compress")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--codec",
        type=str,
        choices=["zlib", "lz4", "zstd"],
        default=None,
        help=(
            "block-compression codec for VTU/VTP (default: zlib). lz4 is a real "
            "VTK compressor; zstd is a meshio++ extension that ParaView cannot "
            "read. Rejected for formats with no such choice rather than being "
            "silently ignored."
        ),
    )
    parser.add_argument(
        "--max",
        "-max",
        action="store_true",
        help="maximum compression",
        default=False,
    )


def compress(args):
    if args.input_format:
        fmts = [args.input_format]
    else:
        fmts = _filetypes_from_path(pathlib.Path(args.infile))
    # pick the first
    fmt = fmts[0]

    # --codec is only meaningful where a block codec is actually chosen.
    # Accepting and ignoring it elsewhere would be the worst outcome: the user
    # would believe they got zstd and silently get gzip (or plain binary).
    if args.codec is not None and fmt not in ("vtu", "vtp"):
        error(
            f"--codec is not applicable to '{fmt}'; it selects the VTK XML "
            "block codec and only vtu/vtp have one."
        )
        exit(1)

    size = os.stat(args.infile).st_size
    print(f"File size before: {size / 1024 ** 2:.2f} MB")
    mesh = read(args.infile, file_format=args.input_format)

    # # Some converters (like VTK) require `points` to be contiguous.
    # mesh.points = np.ascontiguousarray(mesh.points)

    # write it out
    if fmt == "ansys":
        ansys.write(args.infile, mesh, binary=True)
    elif fmt == "cgns":
        cgns.write(
            args.infile, mesh, compression="gzip", compression_opts=9 if args.max else 4
        )
    elif fmt == "gmsh":
        gmsh.write(args.infile, mesh, binary=True)
    elif fmt == "h5m":
        h5m.write(
            args.infile, mesh, compression="gzip", compression_opts=9 if args.max else 4
        )
    elif fmt == "mdpa":
        mdpa.write(args.infile, mesh, binary=True)
    elif fmt == "ply":
        ply.write(args.infile, mesh, binary=True)
    elif fmt == "stl":
        stl.write(args.infile, mesh, binary=True)
    elif fmt == "vtk":
        vtk.write(args.infile, mesh, binary=True)
    elif fmt in ("vtu", "vtp"):
        if args.codec is not None:
            compression = args.codec
        else:
            # Unchanged default: --max still means lzma, so no existing
            # invocation behaves differently.
            compression = "lzma" if args.max else "zlib"
        writer = vtu if fmt == "vtu" else vtp
        writer.write(args.infile, mesh, binary=True, compression=compression)
    elif fmt == "xdmf":
        xdmf.write(
            args.infile,
            mesh,
            data_format="HDF",
            compression="gzip",
            compression_opts=9 if args.max else 4,
        )
    else:
        error(f"Don't know how to compress {args.infile}.")
        exit(1)

    size = os.stat(args.infile).st_size
    print(f"File size after: {size / 1024 ** 2:.2f} MB")
