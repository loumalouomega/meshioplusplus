from .. import _core
from .._helpers import register_format
from ._vtx import read as _py_read
from ._vtx import time_values as _py_time_values
from ._write import SeriesWriter as _PySeriesWriter
from ._write import write as _py_write

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


def write(filename, mesh):
    """Write a DOLFINx VTX file (an ADIOS2 ``.bp`` directory) holding one step:
    the mesh at ``field_data["meshio:time"]`` in ``VTXWriter``'s layout (a
    ``vtk.xml`` schema, ``geometry``, ``connectivity``, ``types``, one variable
    per point or cell array). The file holds one cell type, as DOLFINx's does;
    cells are written as the ``VTK_LAGRANGE_*`` types DOLFINx uses where the
    node order allows. Needs a core built with ``MESHIOPLUSPLUS_WITH_ADIOS2=ON``
    or the ``adios2`` Python package (``pip install meshioplusplus[adios2]``).
    """
    if _HAS_ADIOS2:
        return _core.vtx_write(str(filename), mesh)
    return _py_write(filename, mesh)


class SeriesWriter:
    """Write a time series into one VTX file, a step per mesh;
    ``write_sequence(path, steps, file_format="vtx")`` uses it.

    >>> with meshioplusplus.vtx.SeriesWriter("u.bp") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)
    """

    def __init__(self, filename):
        if _HAS_ADIOS2:
            self._core, self._py = _core.VtxSeriesWriter(str(filename)), None
        else:
            self._core, self._py = None, _PySeriesWriter(filename)

    def write(self, time, mesh):
        if self._core is not None:
            self._core.write(float(time), mesh)
        else:
            self._py.write(time, mesh)

    def close(self):
        if self._core is not None:
            self._core.finalize()
            self._core = None
        if self._py is not None:
            self._py.close()
            self._py = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


register_format("vtx", [".bp"], read, {"vtx": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
