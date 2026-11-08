from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import vtk_xml_shims
from ._vtm import read as _py_read  # noqa: F401
from ._vtm import write as _py_write  # noqa: F401

read, write = vtk_xml_shims(
    __name__,
    "vtm",
    """Read a VTK XML MultiBlock file: an index plus one piece per cell
    block, combined with ``merge()`` (no welding) into one mesh with one
    ``"cell"`` region per piece, named from the index's ``name=`` attribute.""",
    """Write a mesh as VTK XML MultiBlock: one ``.vtu`` piece per cell block,
    pruned to only the points it references, plus an index file.""",
)


register_format("vtm", [".vtm"], read, {"vtm": write})

__all__ = ["read", "write"]
