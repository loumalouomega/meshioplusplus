from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._cgns import read as _py_read
from ._cgns import write as _py_write

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)


def read(filename):
    """Read a CGNS file (C++ core when built with HDF5, Python/h5py fallback)."""
    if _HAS_HDF5 and not is_buffer(filename, "r"):
        try:
            return _core.cgns_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, compression="gzip", compression_opts=4):
    """Write a CGNS file (C++ core when built with HDF5, Python/h5py fallback)."""
    if _HAS_HDF5 and compression in (None, "gzip") and not is_buffer(filename, "w"):
        gzip_level = -1 if compression is None else int(compression_opts or 4)
        try:
            _core.cgns_write(str(filename), mesh, gzip_level)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, compression, compression_opts)


register_format("cgns", [".cgns"], read, {"cgns": write})

__all__ = ["read", "write"]
