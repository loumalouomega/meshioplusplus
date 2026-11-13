from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from .._series_writer import CoreSeriesWriter
from ._cgns import SeriesWriter as _PySeriesWriter
from ._cgns import read as _py_read
from ._cgns import write as _py_write

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)
# The optional cgnslib backend. When present the C++ read path additionally
# handles ADF-container files and the CGNS 3.x NGON_n/NFACE_n layout (the
# 4.x layout is the raw HDF5 path's own, in the core without cgnslib); the
# Python reference below is h5py-based and can read neither, which is why a
# failure on a non-HDF5 file is re-raised rather than falling back into a
# confusing h5py error (see `read`). ADF itself stays out of scope: it is a
# bespoke binary container, and cgnslib is the route.
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


def read(filename, time_step: int = 0, base=None):
    """Read a CGNS file (C++ core when built with HDF5, Python/h5py fallback).

    ``time_step`` selects one step of a transient (`BaseIterativeData_t`/
    `ZoneIterativeData_t`) file (0 = first, negative counts from the end),
    resolved the same way the C API/Fortran/Julia/R/WASM surfaces do -- see
    :func:`meshioplusplus.cgns.read`'s C++ counterpart, ``read_cgns``. A
    non-default value forces the C++ path (the Python h5py reference has no
    transient-step support) and re-raises rather than silently falling back.

    ``base`` selects one base of a multi-base file, by name or by position
    (0 is the first); by default the first is read, with a warning when there
    are more. Reading every base into one mesh is not offered: bases may differ
    in dimension.
    """
    if _HAS_HDF5 and not is_buffer(filename, "r"):
        try:
            if base is not None:
                return _core.cgns_read_base(str(filename), str(base), time_step)
            return _core.cgns_read(str(filename), time_step)
        except Exception as exc:
            # Falling back is only sound when the fallback answers the same
            # question. It cannot for an ADF-container file: the reference
            # reader is h5py-based and would report a confusing signature
            # error rather than the real one. (The xdmf `time_step` precedent.)
            if not is_buffer(filename, "r") and not _is_hdf5(filename):
                raise
            if (
                time_step
                or base is not None
                or not core_declined(exc, "cgns", "read", filename)
            ):
                raise
    return _py_read(filename, base=base)


class SeriesWriter(CoreSeriesWriter):
    """Write a time series into one CGNS file: the mesh with the first step, then one
    ``FlowSolution_t`` per later step, with ``BaseIterativeData_t`` (``TimeValues``) and
    ``ZoneIterativeData_t`` (``FlowSolutionPointers``), which ``read(..., time_step=k)``
    resolves a step through. A step carries point data or cell data, not both; every step
    must have the first step's points and cells.

    >>> with meshioplusplus.cgns.SeriesWriter("run.cgns") as w:
    ...     for t, mesh in steps:
    ...         w.write(t, mesh)
    """

    def __init__(self, filename, compression="gzip", compression_opts=4):
        self._core = self._py = None
        if _HAS_HDF5 and not is_buffer(filename, "w") and compression == "gzip":
            try:
                self._core = _core.CgnsSeriesWriter(str(filename), compression_opts)
                return
            except Exception as exc:
                if not core_declined(exc, "cgns", "write", filename):
                    raise
        self._py = _PySeriesWriter(filename, compression, compression_opts)


def write(filename, mesh, compression="gzip", compression_opts=4):
    """Write a CGNS file (C++ core when built with HDF5, Python/h5py fallback)."""
    if _HAS_HDF5 and compression in (None, "gzip") and not is_buffer(filename, "w"):
        gzip_level = -1 if compression is None else int(compression_opts or 4)
        try:
            _core.cgns_write(str(filename), mesh, gzip_level)
            return
        except Exception as exc:
            if not core_declined(exc, "cgns", "write", filename):
                raise
    return _py_write(filename, mesh, compression, compression_opts)


register_format("cgns", [".cgns"], read, {"cgns": write})

__all__ = ["read", "write", "SeriesWriter"]
