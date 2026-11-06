from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from .._vtk_codecs import CPP_CODECS, cpp_codec_ok
from ._pvd import SeriesWriter
from ._pvd import read as _py_read
from ._pvd import write as _py_write


def read(
    filename,
    points_only=False,
    arrays=None,
    time_step=0,
    piece=None,
    ghosts="keep",
):
    """Read a ParaView ``.pvd`` collection: one step, its parts merged.

    ``time_step`` selects the step (the distinct ``timestep=`` values, ascending;
    negative counts from the end) and ``piece`` selects one ``part`` of that step
    instead of merging them. A step that is a single ``.pvtu`` is returned as
    that file reads, keeping its own ``piece_<i>`` regions. Every step's time is
    in ``read_metadata(path)["time_values"]``, read without opening a piece.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.pvd_read(
                str(filename), points_only, arrays, time_step, piece, ghosts
            )
        except Exception as exc:
            if not core_declined(exc, "pvd", "read", filename):
                raise
    return _py_read(filename, time_step=time_step, piece=piece, ghosts=ghosts)


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    """Write a one-step ``.pvd`` collection plus its ``.vtu`` piece.

    The step's time is ``mesh.field_data["meshio:time"]`` when present, else 0.
    For many steps use ``write_sequence`` with a ``.pvd`` path."""
    if cpp_codec_ok(compression, header_type) and not is_buffer(filename, "w"):
        try:
            _core.pvd_write_codec(str(filename), mesh, binary, CPP_CODECS[compression])
            return
        except Exception as exc:
            if not core_declined(exc, "pvd", "write", filename):
                raise
    return _py_write(
        filename, mesh, binary=binary, compression=compression, header_type=header_type
    )


register_format("pvd", [".pvd"], read, {"pvd": write})

__all__ = ["SeriesWriter", "read", "write"]
