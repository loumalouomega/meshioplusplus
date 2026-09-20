from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._pvtu import read as _py_read
from ._pvtu import write as _py_write
from ._pvtu import write_pieces as _py_write_pieces

_CPP_CODECS = {None: "none", "zlib": "zlib", "lz4": "lz4", "zstd": "zstd"}


def _cpp_codec_ok(compression, header_type):
    if header_type is not None or compression not in _CPP_CODECS:
        return False
    if compression == "lz4" and not getattr(_core, "__has_lz4__", False):
        return False
    if compression == "zstd" and not getattr(_core, "__has_zstd__", False):
        return False
    return True


def read(filename, points_only=False, arrays=None, piece=None, ghosts="keep"):
    """Read a VTK XML parallel index (PUnstructuredGrid): every piece merged into one mesh
    with one ``"cell"`` region per piece (``piece_0``, ``piece_1``, ...), or
    only piece ``piece`` (negative counts from the end; no regions).

    ``ghosts="keep"`` (the default) leaves ghost cells in, ``"drop"`` removes
    every cell with a ``vtkGhostType`` bit set, and the points only they used,
    before merging. Piece paths are resolved against the index's own directory."""
    if not is_buffer(filename, "r"):
        try:
            return _core.pvtu_read(str(filename), points_only, arrays, piece, ghosts)
        except Exception as exc:
            if not core_declined(exc, "pvtu", "read", filename):
                raise
    return _py_read(filename, piece=piece, ghosts=ghosts)


def write(
    filename,
    mesh,
    binary=True,
    compression="zlib",
    header_type=None,
    part_key="partition:part",
):
    """Write a mesh as a VTK XML parallel index plus one ``.vtu`` piece per part.

    A part is one value of the integer ``cell_data[part_key]`` (what
    ``partition_labels`` produces); without that array the whole mesh is one
    piece. To write ``partition`` output with its halo layers use
    ``write_pieces``."""
    if _cpp_codec_ok(compression, header_type) and not is_buffer(filename, "w"):
        try:
            _core.pvtu_write_codec(
                str(filename), mesh, binary, _CPP_CODECS[compression], part_key
            )
            return
        except Exception as exc:
            if not core_declined(exc, "pvtu", "write", filename):
                raise
    return _py_write(
        filename,
        mesh,
        binary=binary,
        compression=compression,
        header_type=header_type,
        part_key=part_key,
    )


def write_pieces(filename, pieces, binary=True, compression="zlib", header_type=None):
    """Write already-carved pieces (e.g. the list ``partition`` returns) as an
    index plus one file each. Every piece must declare identical arrays; that is
    checked before anything is written. A ``partition:ghost`` array becomes
    ``vtkGhostType`` on cells and points and sets the index's ``GhostLevel``."""
    pieces = list(pieces)
    if _cpp_codec_ok(compression, header_type) and not is_buffer(filename, "w"):
        try:
            _core.pvtu_write_pieces_codec(
                str(filename), pieces, binary, _CPP_CODECS[compression]
            )
            return
        except Exception as exc:
            if not core_declined(exc, "pvtu", "write", filename):
                raise
    return _py_write_pieces(
        filename,
        pieces,
        binary=binary,
        compression=compression,
        header_type=header_type,
    )


register_format("pvtu", [".pvtu"], read, {"pvtu": write})

__all__ = ["read", "write", "write_pieces"]
