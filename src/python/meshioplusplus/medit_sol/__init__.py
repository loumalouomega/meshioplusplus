from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._medit_sol import attach as _py_attach
from ._medit_sol import read as _py_read
from ._medit_sol import read_sections as _py_read_sections
from ._medit_sol import write as _py_write


def read(filename):
    """Read a Medit solution: ``.sol`` ascii, or ``.solb`` binary when the name
    ends in ``b``. The result is a geometry-less mesh carrying the vertex fields
    as ``point_data["medit:sol:<i>"]``; element fields need their mesh, through
    ``medit.read(mesh_path, sol=...)``."""
    if not is_buffer(filename, "r"):
        try:
            return _core.medit_sol_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "medit_sol", "read", filename):
                raise
    return _py_read(filename)


def attach(mesh, filename):
    """Put the fields of the solution ``filename`` on ``mesh``: vertex fields as
    point data, element fields as cell data. Returns the mesh with them."""
    if not is_buffer(filename, "r"):
        try:
            attached = _core.medit_sol_read_into(str(filename), mesh)
        except Exception as exc:
            if not core_declined(exc, "medit_sol", "read", filename):
                raise
        else:
            # Only the solution's arrays cross back: the caller keeps its mesh.
            for name, values in attached.point_data.items():
                if name.startswith("medit:sol:"):
                    mesh.point_data[name] = values
            for name, values in attached.cell_data.items():
                if name.startswith("medit:sol:"):
                    mesh.cell_data[name] = values
            return mesh
    dim, sections = _py_read_sections(filename)
    return _py_attach(mesh, dim, sections)


def write(filename, mesh, float_fmt=".16e"):
    """Write the mesh's floating-point point and cell data as a Medit solution:
    ascii ``.sol``, or little-endian binary ``.solb`` when the name ends in ``b``."""
    if not is_buffer(filename, "w") and (
        str(filename).endswith("b") or float_fmt == ".16e"
    ):
        try:
            _core.medit_sol_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "medit_sol", "write", filename):
                raise
    return _py_write(filename, mesh, float_fmt=float_fmt)


register_format("medit_sol", [".sol", ".solb"], read, {"medit_sol": write})

__all__ = ["read", "write", "attach"]
