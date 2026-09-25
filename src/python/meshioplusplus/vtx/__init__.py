from .. import _core
from .._helpers import register_format
from ._vtx import read as _py_read
from ._vtx import time_values as _py_time_values

_HAS_ADIOS2 = getattr(_core, "__has_adios2__", False)


def read(filename, points_only=False, arrays=None, time_step=0, ghosts="keep"):
    """Read one step of a DOLFINx VTX file (an ADIOS2 ``.bp`` directory).

    The steps are the file's ADIOS2 steps, each with its time in
    ``field_data["meshio:time"]``; a step without mesh variables (a
    ``VTXMeshPolicy.reuse`` file) takes the mesh of the last step that has one.
    The MPI ranks' blocks are concatenated as ParaView merges them;
    ``ghosts="drop"`` instead welds every ghost point onto its owner through
    ``vtkOriginalPointIds``, which gives the serial mesh.

    Needs a core built with ``MESHIOPLUSPLUS_WITH_ADIOS2=ON`` or the ``adios2``
    Python package (``pip install meshioplusplus[adios2]``).
    """
    # No fallback from a core built with ADIOS2: the native reader is the whole
    # of this format, and the Python twin would load the `adios2` package's own
    # ADIOS2 libraries into a process that already holds the core's.
    if _HAS_ADIOS2:
        return _core.vtx_read(
            str(filename),
            points_only=points_only,
            arrays=arrays,
            time_step=time_step,
            ghosts=ghosts,
        )
    return _py_read(filename, points_only, arrays, time_step, ghosts)


def time_values(filename):
    """The ``step`` time of every step (the step index where it is absent)."""
    if _HAS_ADIOS2:
        return list(_core.vtx_time_values(str(filename)))
    return _py_time_values(filename)


register_format("vtx", [".bp"], read, {})

__all__ = ["read", "time_values"]
