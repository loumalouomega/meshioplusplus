"""VTKHDF (``.vtkhdf``) -- Kitware's HDF5-based VTK format.

Registers ``vtkhdf`` for ``.vtkhdf`` and ``.hdf``. ``.hdf`` is a generic HDF5
extension: a file with no ``/VTKHDF`` group is refused with a message telling the
caller to name the format explicitly.

Both directions use the C++ core when it is built with HDF5 and fall back to the
h5py reference (``_vtkhdf.py``) otherwise, or when the core declines a file.
"""

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._vtkhdf import read as _py_read
from ._vtkhdf import write as _py_write
from .time_series import TimeSeriesReader, TimeSeriesWriter

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)


def read(
    filename,
    points_only=False,
    arrays=None,
    time_step=0,
    piece=None,
    lenient=False,
):
    """Read a VTKHDF file (C++ core when built with HDF5, h5py fallback).

    ``piece=None`` merges every piece of a partitioned file or composite into one
    mesh with one cell region per piece; an int keeps only that piece (negative
    counts from the end). ``time_step`` picks a step of a transient file.
    ``lenient`` skips cells with no meshio++ type (poly-vertex, poly-line,
    triangle strips) with a warning instead of raising.
    """
    if _HAS_HDF5 and not is_buffer(filename, "r"):
        try:
            return _core.vtkhdf_read(
                str(filename), points_only, arrays, time_step, piece, lenient
            )
        except Exception as exc:
            if not core_declined(exc, "vtkhdf", "read", filename):
                raise
    return _py_read(filename, points_only, arrays, time_step, piece, lenient)


def write(
    filename,
    mesh,
    compression="gzip",
    compression_opts=4,
    dataset_type="UnstructuredGrid",
    version=None,
):
    """Write a static VTKHDF file (C++ core when built with HDF5, h5py fallback).

    ``dataset_type`` is ``"UnstructuredGrid"`` (default), ``"PolyData"``,
    ``"PartitionedDataSetCollection"`` or ``"MultiBlockDataSet"``; ``version`` is
    ``None`` (the oldest version covering the content) or ``(major, minor)``.
    """
    if _HAS_HDF5 and compression in (None, "gzip") and not is_buffer(filename, "w"):
        if compression is None:
            gzip_level = -1
        else:
            gzip_level = 4 if compression_opts is None else int(compression_opts)
        try:
            _core.vtkhdf_write(
                str(filename),
                mesh,
                gzip_level,
                dataset_type,
                None if version is None else tuple(int(v) for v in version),
            )
            return
        except Exception as exc:
            if not core_declined(exc, "vtkhdf", "write", filename):
                raise
    return _py_write(
        filename, mesh, compression, compression_opts, dataset_type, version
    )


register_format("vtkhdf", [".vtkhdf", ".hdf"], read, {"vtkhdf": write})

__all__ = ["read", "write", "TimeSeriesReader", "TimeSeriesWriter"]
