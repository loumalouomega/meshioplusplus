"""UGRID-conventions NetCDF meshes (roadmap 2.4.3.2): C++ core for real file paths,
the ``netCDF4`` reference otherwise."""

from .. import _core
from .._fallback import core_declined
from .._helpers import register_format
from ._ugrid_nc import read as _py_read
from ._ugrid_nc import time_values as _py_time_values
from ._ugrid_nc import write as _py_write


def read(filename, time_step=0):
    """Read a UGRID NetCDF mesh: the first ``mesh_topology``'s nodes, faces (grouped
    by node count) and edges, its ``node``/``face``/``edge`` data, and one step of a
    time series (``time_step``; 0 is the first, negatives count from the end, its
    time in ``field_data["meshio:time"]``)."""
    try:
        return _core.ugrid_nc_read(str(filename), int(time_step))
    except Exception as exc:
        if not core_declined(exc, "ugrid_nc", "read", filename):
            raise
    return _py_read(filename, time_step=time_step)


def time_values(filename):
    """Every step's time (empty for a static file)."""
    try:
        return list(_core.ugrid_nc_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "ugrid_nc", "read", filename):
            raise
    return _py_time_values(filename)


def write(filename, mesh):
    """Write a 1-D or 2-D mesh (triangles, quads, polygons, lines) as UGRID NetCDF."""
    try:
        _core.ugrid_nc_write(str(filename), mesh)
        return
    except Exception as exc:
        if not core_declined(exc, "ugrid_nc", "write", filename):
            raise
    _py_write(filename, mesh)


register_format("ugrid_nc", [".nc"], read, {"ugrid_nc": write})

__all__ = ["read", "write", "time_values"]
