from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._nastran_h5 import read as _py_read
from ._nastran_h5 import time_values

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


register_format("nastran_h5", [".h5"], read, {})

__all__ = ["read", "time_values"]
