"""3MF packages (roadmap 2.4.5): C++ core for real file paths, the Python reference
otherwise. The format is named ``3mf``; the package ``threemf``."""

from .. import _core
from .._fallback import core_declined
from .._helpers import register_format
from ._threemf import read as _py_read
from ._threemf import write as _py_write


def read(filename):
    """Read a 3MF package's build: every item's object (meshes and components, with
    their transforms) in one triangle block, a cell region per item, the unit as
    ``meshio:length_unit`` and base materials as the ``3mf:color`` RGBA cell data."""
    try:
        return _core.threemf_read(str(filename))
    except Exception as exc:
        if not core_declined(exc, "3mf", "read", filename):
            raise
    return _py_read(filename)


def write(filename, mesh):
    """Write a triangle mesh as a 3MF package: an object per cell region, the unit
    from ``meshio:length_unit`` and base materials from ``3mf:color``."""
    try:
        _core.threemf_write(str(filename), mesh)
        return
    except Exception as exc:
        if not core_declined(exc, "3mf", "write", filename):
            raise
    _py_write(filename, mesh)


register_format("3mf", [".3mf"], read, {"3mf": write})

__all__ = ["read", "write"]
