from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import vtk_xml_shims
from ._vti import read as _py_read  # noqa: F401
from ._vti import write as _py_write  # noqa: F401

read, write = vtk_xml_shims(
    __name__,
    "vti",
    """Read a VTK XML ImageData file, expanding the extent into explicit cells.""",
    """Write a mesh as VTK XML ImageData. The mesh must be a dense lattice.""",
)


register_format("vti", [".vti"], read, {"vti": write})

__all__ = ["read", "write"]
