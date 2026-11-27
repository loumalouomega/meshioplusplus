"""UGRID-conventions NetCDF meshes (roadmap 2.4.3.2): C++ core for real file paths,
the ``netCDF4`` reference otherwise."""

from .. import _core  # noqa: F401
from .._helpers import register_format
from .._shims import series_shims
from ._ugrid_nc import read as _py_read  # noqa: F401
from ._ugrid_nc import time_values as _py_time_values  # noqa: F401
from ._ugrid_nc import write as _py_write  # noqa: F401

read, write, time_values = series_shims(
    __name__,
    "ugrid_nc",
    """Read a UGRID NetCDF mesh: the first ``mesh_topology``'s nodes, faces (grouped
    by node count) and edges, its ``node``/``face``/``edge`` data, and one step of a
    time series (``time_step``; 0 is the first, negatives count from the end, its
    time in ``field_data["meshio:time"]``).""",
    """Write a 1-D or 2-D mesh (triangles, quads, polygons, lines) as UGRID NetCDF.""",
    """Every step's time (empty for a static file).""",
    selective=False,
)

register_format("ugrid_nc", [".nc"], read, {"ugrid_nc": write})

__all__ = ["read", "write", "time_values"]
