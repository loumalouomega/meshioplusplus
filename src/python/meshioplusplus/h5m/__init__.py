from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._h5m import read as _py_read
from ._h5m import write as _py_write

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)


def read(filename):
    """Read a MOAB h5m file (C++ core when built with HDF5, Python/h5py fallback)."""
    if _HAS_HDF5 and not is_buffer(filename, "r"):
        try:
            return _core.h5m_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, add_global_ids=True, compression="gzip", compression_opts=4):
    """Write a MOAB h5m file (C++ core when built with HDF5, Python/h5py fallback)."""
    if (
        _HAS_HDF5
        and compression in (None, "gzip")
        and not mesh.cell_data
        and not is_buffer(filename, "w")
    ):
        gzip_level = -1 if compression is None else int(compression_opts or 4)
        try:
            _core.h5m_write(str(filename), mesh, bool(add_global_ids), gzip_level)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, add_global_ids, compression, compression_opts)


register_format("h5m", [".h5m"], read, {"h5m": write})

__all__ = ["read", "write"]
