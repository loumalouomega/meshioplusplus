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

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read an XDMF file (C++ core; HDF DataItems need an HDF5-enabled build).

    ``time_step`` selects one step of a temporal collection (0 = the first,
    negative counts from the end). Unlike ``points_only``/``arrays`` it is
    *not* a narrowing option that a post-filter could emulate, so a non-default
    step must not silently fall back to the pure-Python reader, which always
    returns step 0 -- see the guard below.
    """
    # points_only/arrays reach the C++ reader, which skips the unwanted
    # <DataArray>/section bodies outright. The Python fallback below has no
    # selective support, so _helpers.read trims its result instead -- same
    # answer, just without the saving.
    if not is_buffer(filename, "r"):
        try:
            return _core.xdmf_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
            )
        except Exception:
            if time_step:
                # Falling back here would quietly hand back step 0 under the
                # name of step N -- a wrong answer, not a slower one.
                raise
    return _py_read(filename)


def write(filename, mesh, data_format="HDF", **kwargs):
    """Write an XDMF file (C++ core; HDF needs an HDF5-enabled build)."""
    compression = kwargs.get("compression", "gzip")
    compression_opts = kwargs.get("compression_opts", 4)
    cpp_ok = data_format in ("XML", "Binary") or (
        data_format == "HDF" and _HAS_HDF5 and compression in (None, "gzip")
    )
    if cpp_ok and not is_buffer(filename, "w"):
        gzip_level = -1 if compression is None else int(compression_opts or 4)
        try:
            _core.xdmf_write(str(filename), mesh, data_format, gzip_level)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, data_format=data_format, **kwargs)


register_format("xdmf", [".xdmf", ".xmf"], read, {"xdmf": write})

__all__ = ["read", "write", "TimeSeriesWriter", "TimeSeriesReader"]
