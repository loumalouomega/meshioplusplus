from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from . import _vtk_42
from ._main import read as _py_read
from ._main import write as _main_write


def _cpp_ok(mesh):
    if any(c.type.startswith("polyhedron") for c in mesh.cells):
        return False
    # The Python writer pads 2-component vectors to 3 (mutating the mesh); the
    # C++ path doesn't, so defer those to the Python writer.
    for v in mesh.point_data.values():
        if v.ndim == 2 and v.shape[1] == 2:
            return False
    for blocks in mesh.cell_data.values():
        for v in blocks:
            if getattr(v, "ndim", 1) == 2 and v.shape[1] == 2:
                return False
    return True


def read(filename):
    """Read a VTK legacy file.

    Uses the C++ core for version 5.1 UNSTRUCTURED_GRID files (ascii or
    big-endian binary), falling back to the reference Python reader otherwise
    (version 4.2, structured grids, SCALARS/VECTORS sections, polyhedron).
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.vtk_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, fmt_version="5.1", binary=True, **kwargs):
    """Write a VTK legacy file.

    Uses the C++ core for versions 5.1 (default) and 4.2, ascii or big-endian
    binary, on supported meshes; otherwise falls back to the Python writer.
    """
    if fmt_version in ("5.1", "4.2") and not is_buffer(filename, "w") and _cpp_ok(mesh):
        try:
            _core.vtk_write(str(filename), mesh, binary, fmt_version == "5.1")
            return
        except Exception:
            pass
    return _main_write(filename, mesh, fmt_version=fmt_version, binary=binary, **kwargs)


register_format(
    "vtk",
    [".vtk"],
    read,
    {
        "vtk42": _vtk_42.write,
        "vtk51": _vtk_42.write,
        "vtk": write,
    },
)

__all__ = ["read", "write"]
