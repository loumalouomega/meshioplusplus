from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import vtk_xml_shims
from ._vtr import read as _py_read  # noqa: F401
from ._vtr import write as _py_write  # noqa: F401

read, write = vtk_xml_shims(
    __name__,
    "vtr",
    """Read a VTK XML RectilinearGrid file. Fully general: coordinates need
    not be uniformly spaced (a graded mesh reads correctly).""",
    """Write a mesh as VTK XML RectilinearGrid.

    Unlike the reader, the mesh must be a UNIFORM dense lattice today: a
    genuinely graded mesh cannot be written as ``.vtr`` yet.
    """,
)


register_format("vtr", [".vtr"], read, {"vtr": write})

__all__ = ["read", "write"]
