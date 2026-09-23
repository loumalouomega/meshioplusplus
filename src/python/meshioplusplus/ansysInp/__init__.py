from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._ansysInp import read as _py_read
from ._ansysInp import write as _py_write


def read(filename, lenient=False):
    """Read an Ansys MAPDL coded database (``.cdb``).

    Blocks are sliced by their own Fortran format lines. Elements become cells
    by their element type's category (degenerate bricks as wedges, pyramids and
    tetrahedra; shells with K == L as triangles), with ``ansys:element``,
    ``ansys:type``, ``ansys:mat``, ``ansys:real`` and ``ansys:secnum`` cell data;
    ``CMBLOCK`` components become point and cell regions (so also
    ``point_sets``/``cell_sets``). ``lenient`` skips elements whose type has no
    meshio++ cell instead of failing.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.ansysinp_read(str(filename), lenient=lenient)
        except Exception as exc:
            if not core_declined(exc, "ansysInp", "read", filename):
                raise
    return _py_read(filename, lenient=lenient)


def write(filename, mesh):
    """Write an Ansys MAPDL coded database (``.cdb``).

    Every cell type in the layout MAPDL expects (wedges, pyramids and tetrahedra
    as degenerate bricks, triangles as degenerate quads), ``ansys:*`` cell data
    kept where it fits, point and cell regions as ``CMBLOCK`` components.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.ansysinp_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "ansysInp", "write", filename):
                raise
    return _py_write(filename, mesh)


register_format("ansysInp", [".cdb", ".inp"], read, {"ansysInp": write})

__all__ = ["read", "write"]
