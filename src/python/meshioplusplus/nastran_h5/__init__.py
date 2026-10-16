from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._nastran_h5 import read as _py_read
from ._nastran_h5 import time_values
from ._write import SeriesWriter as _PySeriesWriter
from ._write import write as _py_write

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read an MSC Nastran HDF5 result database (``.h5``, ``MDLPRM,HDF5``).

    The model comes from ``/NASTRAN/INPUT``: GRID points (coordinates as written; a
    ``CP != 0`` is warned about and kept as ``point_data["nastran:cp"]``), one cell
    block per element card and order, ``cell_data["nastran:eid"]`` and
    ``["nastran:pid"]``, and one cell region per property id named after its
    property card (``PSHELL_4``). Every result domain is a step: ``time_step``
    selects one (0 = first, negative counts from the end), its ``TIME_FREQ_EIGR`` is
    ``field_data["meshio:time"]`` and ``nastran:subcase``/``mode``/``analysis``...
    sit beside it. Nodal tables become ``point_data`` (``DISPLACEMENT``,
    ``DISPLACEMENT_ROT``, ``EIGENVECTOR_real``...), element tables ``cell_data``
    (``STRESS:X``...). Other vendors' HDF5 schemas are refused. The C++ core reads
    the file; the Python reader (h5py) is the fallback.
    """
    if _HAS_HDF5 and not is_buffer(filename, "r"):
        try:
            return _core.nastran_h5_read(str(filename), points_only, arrays, time_step)
        except Exception as exc:
            if not core_declined(exc, "nastran_h5", "read", filename):
                raise
    return _py_read(
        filename, points_only=points_only, arrays=arrays, time_step=time_step
    )


def write(filename, mesh):
    """Write a mesh as an MSC Nastran HDF5 result database (one domain).

    The model is ``/NASTRAN/INPUT`` (GRIDs, one table per element card, the
    property and material tables the elements name); ``nastran:eid`` and
    ``nastran:pid`` cell data give the ids. Point data becomes the nodal tables
    (``DISPLACEMENT`` and ``DISPLACEMENT_ROT`` of 3 components, ``T`` of one,
    ``T:M`` members, ``_real``/``_imag`` pairs as ``_CPLX``) and ``<G>:<M>``
    cell data the element tables of the centre values, each with its
    ``/INDEX`` twin; ``meshio:time`` and ``nastran:*`` field data fill the
    ``DOMAINS`` row. Needs ``h5py`` (the Python twin) or a core built with HDF5.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.nastran_h5_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "nastran_h5", "write", filename):
                raise
    return _py_write(filename, mesh)


class SeriesWriter:
    """Write a time series into one results file: the model once and a domain
    per step; ``write_sequence(path, steps, file_format="nastran_h5")`` uses
    it. The C++ core writes the file; the Python twin, which writes the same
    tables, when the core declines to.

    >>> with meshioplusplus.nastran_h5.SeriesWriter("job.h5") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)
    """

    def __init__(self, filename):
        self._core = self._py = None
        if _HAS_HDF5 and not is_buffer(filename, "w"):
            try:
                self._core = _core.NastranH5SeriesWriter(str(filename))
            except Exception as exc:
                if not core_declined(exc, "nastran_h5", "write", filename):
                    raise
        if self._core is None:
            self._py = _PySeriesWriter(filename)

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


register_format("nastran_h5", [".h5"], read, {"nastran_h5": write})

__all__ = ["read", "write", "time_values", "SeriesWriter"]
