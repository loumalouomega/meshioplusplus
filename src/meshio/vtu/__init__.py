from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._vtu import read
from ._vtu import write as _py_write


def _has_polyhedron(mesh):
    return any(c.type.startswith("polyhedron") for c in mesh.cells)


def write(filename, mesh, binary=True, compression="zlib", header_type=None):
    """Write a VTU file.

    Uses the C++ core for the cases it currently supports (ASCII output of
    non-polyhedron meshes to a real file path) and otherwise falls back to the
    reference Python writer. The fallback also catches any limitation hit by
    the C++ path, so behaviour is identical to the pure-Python implementation.
    """
    if (
        not binary
        and header_type is None
        and not is_buffer(filename, "w")
        and not _has_polyhedron(mesh)
    ):
        try:
            _core.vtu_write_ascii(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(
        filename, mesh, binary=binary, compression=compression, header_type=header_type
    )


register_format("vtu", [".vtu"], read, {"vtu": write})

__all__ = ["read", "write"]
