from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._d3plot import read as _py_read
from ._d3plot import time_values as _py_time_values
from ._write import SeriesWriter as _PySeriesWriter
from ._write import write as _py_write


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read one state of an LS-DYNA d3plot family.

    Open the base file (``d3plot``): its numbered members (``d3plot01``...)
    are found beside it. ``time_step`` picks the state (negative counts from
    the end); its time is ``field_data["meshio:time"]``. Nodal variables are
    point data (``displacement`` is relative to the initial coordinates, which
    stay the points), element variables cell data (NaN on blocks without them;
    per integration point or layer as ``(cells, points * components)`` with
    ``field_data["lsdyna_d3plot:layout:<name>"]``), the deletion flags the int8
    mask ``lsdyna:alive`` and the globals field data. Parts are cell regions.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.lsdyna_d3plot_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
            )
        except Exception as exc:
            if not core_declined(exc, "lsdyna_d3plot", "read", filename):
                raise
    return _py_read(
        filename, points_only=points_only, arrays=arrays, time_step=time_step
    )


def time_values(filename):
    """The times of the family's states."""
    try:
        return list(_core.lsdyna_d3plot_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "lsdyna_d3plot", "read", filename):
            raise
    return _py_time_values(filename)


# Found by its file name (``d3plot``) or its control block, never an extension.
def write(filename, mesh, double=False):
    """Write one state as an LS-DYNA d3plot (``double=True``: 8-byte words).

    Hexahedra, wedges, pyramids, tetrahedra and ``tetra10`` are solids, ``line``
    beams, ``triangle`` and ``quad`` shells; other cells are not written. The
    variables go by the names the reader returns: point data ``displacement``,
    ``velocity``, ``acceleration``, ``temperature``, ``heat_flux``,
    ``mass_scaling``, ``temperature_gradient``, ``residual_forces``,
    ``residual_moments`` and ``lsdyna:alive``; cell data ``stress``,
    ``effective_plastic_strain``, ``strain``, ``plastic_strain_tensor``,
    ``thermal_strain_tensor``, ``history_variables``, ``thermal_variables``, the
    shell resultants and extras, the ``beam_*`` resultants and layers and
    ``lsdyna:alive``; the ``global_*`` and ``part_*`` field data; ``lsdyna:nid``,
    ``lsdyna:eid`` and ``lsdyna:part`` the user ids, cell regions the parts and
    their titles; ``field_data["meshio:time"]`` the state's time.
    """
    if not is_buffer(filename, "w"):
        try:
            return _core.lsdyna_d3plot_write(str(filename), mesh, double=double)
        except Exception as exc:
            if not core_declined(exc, "lsdyna_d3plot", "write", filename):
                raise
    return _py_write(filename, mesh, double=double)


class SeriesWriter:
    """Write a time series as one d3plot family: the geometry once and a state
    per step, split over ``d3plot``, ``d3plot01``... once a file would pass
    ``max_file_bytes`` (default 1 GiB). Files of the same name already beside it
    (``d3plot01``...) are replaced.

    >>> with meshioplusplus.lsdyna_d3plot.SeriesWriter("run/d3plot") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)
    """

    def __init__(self, filename, max_file_bytes=1 << 30, double=False):
        self._core = self._py = None
        if not is_buffer(filename, "w"):
            try:
                self._core = _core.LsdynaD3plotSeriesWriter(
                    str(filename), int(max_file_bytes), bool(double)
                )
            except Exception as exc:
                if not core_declined(exc, "lsdyna_d3plot", "write", filename):
                    raise
        if self._core is None:
            self._py = _PySeriesWriter(filename, max_file_bytes, double)

    def write(self, time, mesh):
        if self._core is not None:
            self._core.write(float(time), mesh)
        else:
            self._py.write(time, mesh)

    def close(self):
        if self._core is not None:
            if self._core.num_steps():
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


register_format("lsdyna_d3plot", [], read, {"lsdyna_d3plot": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
