from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._cgns import read as _py_read
from ._cgns import write as _py_write

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)
# The optional cgnslib backend. When present the C++ read path additionally
# handles ADF-container files and NGON_n/NFACE_n polyhedral sections; the
# Python reference below is h5py-based and can do neither, which is why a
# failure on a non-HDF5 file is re-raised rather than falling back into a
# confusing h5py error (see `read`).
_HAS_CGNSLIB = getattr(_core, "__has_cgnslib__", False)


def _is_hdf5(filename):
    """Whether the file starts with the HDF5 signature.

    `.cgns` has two on-disk containers, HDF5 and ADF. The Python reference
    reader is h5py-based, so an ADF file is unreadable there no matter what --
    worth saying by name instead of surfacing "file signature not found".
    """
    try:
        with open(filename, "rb") as fh:
            return fh.read(8) == b"\x89HDF\r\n\x1a\n"
    except OSError:
        return True  # not our problem to diagnose here; let the reader report


def read(filename):
    """Read a CGNS file (C++ core when built with HDF5, Python/h5py fallback)."""
    if _HAS_HDF5 and not is_buffer(filename, "r"):
        try:
            return _core.cgns_read(str(filename))
        except Exception:
            # Falling back is only sound when the fallback answers the same
            # question. It cannot for an ADF-container file: the reference
            # reader is h5py-based and would report a confusing signature
            # error rather than the real one. (The xdmf `time_step` precedent.)
            if not is_buffer(filename, "r") and not _is_hdf5(filename):
                raise
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
