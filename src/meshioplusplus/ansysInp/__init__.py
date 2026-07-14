from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._ansysInp import read as _py_read
from ._ansysInp import write as _py_write


def read(filename):
    """Read an Ansys MAPDL .cdb/.inp file (C++ core, Python fallback)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.ansysinp_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write an Ansys MAPDL .cdb/.inp file (C++ core, Python fallback)."""
    if not is_buffer(filename, "w"):
        point_sets = dict(mesh.point_sets)
        cell_sets = {k: list(v) for k, v in mesh.cell_sets.items()}
        try:
            _core.ansysinp_write(str(filename), mesh, point_sets, cell_sets)
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("ansysInp", [".cdb", ".inp"], read, {"ansysInp": write})

__all__ = ["read", "write"]
