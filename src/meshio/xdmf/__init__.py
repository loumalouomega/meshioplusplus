"""
I/O for XDMF.
https://xdmf.org/index.php/XDMF_Model_and_Format
"""

from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from .main import read as _py_read
from .main import write as _py_write
from .time_series import TimeSeriesReader, TimeSeriesWriter


def read(filename):
    """Read an XDMF file (C++ core for XML/Binary data, Python/h5py for HDF)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.xdmf_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, data_format="HDF", **kwargs):
    """Write an XDMF file (C++ core for XML/Binary data, Python/h5py for HDF)."""
    if data_format in ("XML", "Binary") and not is_buffer(filename, "w"):
        try:
            _core.xdmf_write(str(filename), mesh, data_format)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, data_format=data_format, **kwargs)


register_format("xdmf", [".xdmf", ".xmf"], read, {"xdmf": write})

__all__ = ["read", "write", "TimeSeriesWriter", "TimeSeriesReader"]
