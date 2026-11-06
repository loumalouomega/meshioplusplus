"""The ``read``/``write``/``write_pieces`` shims of the parallel VTK XML formats.

``pvtu`` and ``pvtp`` differ only in their format token, the dataset their
index names and the extension of their pieces, so both build their public
functions here: the C++ core first, the Python reference when the core
declines (``_fallback.core_declined``).
"""

from . import _core
from ._fallback import core_declined
from ._files import is_buffer

_CPP_CODECS = {
    None: "none",
    "zlib": "zlib",
    "lz4": "lz4",
    "zstd": "zstd",
    "lzma": "lzma",
}


def _cpp_codec_ok(compression, header_type):
    if header_type is not None or compression not in _CPP_CODECS:
        return False
    if compression == "lz4" and not getattr(_core, "__has_lz4__", False):
        return False
    if compression == "zstd" and not getattr(_core, "__has_zstd__", False):
        return False
    if compression == "lzma" and not getattr(_core, "__has_lzma__", False):
        return False
    return True


def make_parallel_format(fmt, dataset, piece_ext, py_read, py_write, py_write_pieces):
    """``(read, write, write_pieces)`` for one parallel format.

    ``fmt`` is the format token (``"pvtu"``), ``dataset`` the index's dataset
    type for the docstring (``"PUnstructuredGrid"``) and ``piece_ext`` the
    pieces' extension (``".vtu"``)."""

    def read(filename, points_only=False, arrays=None, piece=None, ghosts="keep"):
        if not is_buffer(filename, "r"):
            try:
                return getattr(_core, f"{fmt}_read")(
                    str(filename), points_only, arrays, piece, ghosts
                )
            except Exception as exc:
                if not core_declined(exc, fmt, "read", filename):
                    raise
        return py_read(filename, piece=piece, ghosts=ghosts)

    def write(
        filename,
        mesh,
        binary=True,
        compression="zlib",
        header_type=None,
        part_key="partition:part",
    ):
        if _cpp_codec_ok(compression, header_type) and not is_buffer(filename, "w"):
            try:
                getattr(_core, f"{fmt}_write_codec")(
                    str(filename), mesh, binary, _CPP_CODECS[compression], part_key
                )
                return
            except Exception as exc:
                if not core_declined(exc, fmt, "write", filename):
                    raise
        return py_write(
            filename,
            mesh,
            binary=binary,
            compression=compression,
            header_type=header_type,
            part_key=part_key,
        )

    def write_pieces(
        filename, pieces, binary=True, compression="zlib", header_type=None
    ):
        pieces = list(pieces)
        if _cpp_codec_ok(compression, header_type) and not is_buffer(filename, "w"):
            try:
                getattr(_core, f"{fmt}_write_pieces_codec")(
                    str(filename), pieces, binary, _CPP_CODECS[compression]
                )
                return
            except Exception as exc:
                if not core_declined(exc, fmt, "write", filename):
                    raise
        return py_write_pieces(
            filename,
            pieces,
            binary=binary,
            compression=compression,
            header_type=header_type,
        )

    read.__doc__ = f"""Read a VTK XML parallel index ({dataset}): every piece merged into one mesh
    with one ``"cell"`` region per piece (``piece_0``, ``piece_1``, ...), or
    only piece ``piece`` (negative counts from the end; no regions).

    ``ghosts="keep"`` (the default) leaves ghost cells in, ``"drop"`` removes
    every cell with a ``vtkGhostType`` bit set, and the points only they used,
    before merging. Piece paths are resolved against the index's own directory."""
    write.__doc__ = f"""Write a mesh as a VTK XML parallel index plus one ``{piece_ext}`` piece per part.

    A part is one value of the integer ``cell_data[part_key]`` (what
    ``partition_labels`` produces); without that array the whole mesh is one
    piece. To write ``partition`` output with its halo layers use
    ``write_pieces``."""
    write_pieces.__doc__ = """Write already-carved pieces (e.g. the list ``partition`` returns) as an
    index plus one file each. Every piece must declare identical arrays; that is
    checked before anything is written. A ``partition:ghost`` array becomes
    ``vtkGhostType`` on cells and points and sets the index's ``GhostLevel``."""
    module = f"meshioplusplus.{fmt}"
    for fn in (read, write, write_pieces):
        fn.__module__ = module
        fn.__qualname__ = fn.__name__
    return read, write, write_pieces
