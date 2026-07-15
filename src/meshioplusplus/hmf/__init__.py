from .. import _core
from .._common import warn
from .._files import is_buffer
from .._helpers import register_format
from ._hmf import read as _py_read
from ._hmf import write as _py_write

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)


def read(filename):
    """Read an HMF file (C++ core when built with HDF5, Python/h5py fallback)."""
    if _HAS_HDF5 and not is_buffer(filename, "r"):
        try:
            return _core.hmf_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, compression="gzip", compression_opts=4):
    """Write an HMF file (C++ core when built with HDF5, Python/h5py fallback)."""
    if _HAS_HDF5 and compression in (None, "gzip") and not is_buffer(filename, "w"):
        warn("Experimental file format. Format can change at any time.")
        gzip_level = -1 if compression is None else int(compression_opts or 4)
        try:
            _core.hmf_write(str(filename), mesh, gzip_level)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, compression, compression_opts)


register_format("hmf", [".hmf"], read, {"hmf": write})

__all__ = ["read", "write"]
