from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._abaqus import read as _py_read
from ._abaqus import write as _py_write


def read(filename):
    """Read an Abaqus .inp file.

    The C++ core handles the whole format the Python reference does — ``*NODE``,
    ``*ELEMENT`` (including a trailing ``ELSET=``), ``*NSET``/``*ELSET`` (with
    ``GENERATE`` and set-of-set references) and ``*INCLUDE`` — plus ``*SURFACE``,
    which becomes a ``side`` region. The Python reader stays as the fallback for
    anything the C++ path raises on.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.abaqus_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, float_fmt=".16e", translate_cell_names=True):
    """Write an Abaqus .inp file.

    Named regions are written as ``*NSET`` / ``*ELSET`` / ``*SURFACE``, so a mesh
    carrying sets no longer needs the Python writer.
    """
    if float_fmt == ".16e" and translate_cell_names and not is_buffer(filename, "w"):
        try:
            _core.abaqus_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(
        filename, mesh, float_fmt=float_fmt, translate_cell_names=translate_cell_names
    )


register_format("abaqus", [".inp"], read, {"abaqus": write})

__all__ = ["read", "write"]
