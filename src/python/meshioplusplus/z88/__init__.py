from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._z88 import read as _py_read
from ._z88 import write as _py_write


def read(filename, results=True):
    """Read a Z88 structure file (``z88i1.txt``, ``z88structure.txt``).

    Elements become cells (the Z88 type is the ``z88:type`` cell data). With
    ``results``, ``z88o2.txt`` next to it becomes ``point_data["U"]`` and the
    solid and plane-stress blocks of ``z88o3.txt`` the element mean stress
    ``cell_data["SIG"]`` (plus ``"SIGV"``). A ``z88o2.txt``/``z88o3.txt`` path
    reads the structure file next to it. The C++ core reads the whole format;
    the Python reader is the fallback for buffers and for anything the C++ path
    raises on.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.z88_read(str(filename), results)
        except Exception as exc:
            if not core_declined(exc, "z88", "read", filename):
                raise
    return _py_read(filename, results=results)


def write(filename, mesh, stubs=False):
    """Write a Z88OS v15 structure file (``z88i1.txt``).

    Element types come from ``z88:type`` (else from the cell type); cells with no
    Z88 type, regions and other data arrays are dropped with a warning. With
    ``stubs``, empty ``z88i2.txt`` and ``z88i5.txt`` are written next to it.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.z88_write(str(filename), mesh, stubs)
            return
        except Exception as exc:
            if not core_declined(exc, "z88", "write", filename):
                raise
    return _py_write(filename, mesh, stubs=stubs)


# Z88's files have fixed names (z88i1.txt ...): _helpers matches the basename;
# ".txt" itself stays the xyz reader's.
register_format("z88", [], read, {"z88": write})

__all__ = ["read", "write"]
