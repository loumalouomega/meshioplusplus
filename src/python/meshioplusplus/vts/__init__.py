from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import vtk_xml_shims
from ._vts import read as _py_read  # noqa: F401
from ._vts import write as _py_write  # noqa: F401

read, write = vtk_xml_shims(
    __name__,
    "vts",
    """Read a VTK XML StructuredGrid file.""",
    """Write a mesh as VTK XML StructuredGrid. The mesh must be a dense lattice.""",
)


register_format("vts", [".vts"], read, {"vts": write})

__all__ = ["read", "write"]
