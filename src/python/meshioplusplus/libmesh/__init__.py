import bz2
import gzip
import os
import tempfile

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._libmesh import read as _py_read
from ._libmesh import write as _py_write

_EXTENSIONS = [".xda", ".xdr", ".xda.gz", ".xdr.gz", ".xda.bz2", ".xdr.bz2"]


def _is_bzip2(filename):
    # A file that cannot be opened is left to the reader, which says why.
    try:
        with open(filename, "rb") as fh:
            return fh.read(3) == b"BZh"
    except OSError:
        return False


def read(filename):
    """Read a libMesh ``.xda`` (ASCII) or ``.xdr`` (XDR binary) mesh.

    Active (leaf) elements become the cells, the subdomain id the
    ``libmesh:subdomain`` cell data and a cell region per subdomain (named by the
    file's subdomain map, else ``subdomain_<id>``). Side sets become side regions
    (``boundary_<id>`` when unnamed), node sets point regions (``nodeset_<id>``),
    edge sets ``line`` cells in a ``<name>:edge`` cell region and shell-face sets
    ``<name>:shellface0``/``1`` cell regions. The encoding is detected from the
    content, and gzip/bzip2 compression is inflated. The C++ core reads the
    whole format except bzip2, which the Python reader inflates; the Python
    reader is also the fallback for buffers and for anything the C++ path
    raises on.
    """
    if not is_buffer(filename, "r") and not _is_bzip2(filename):
        try:
            return _core.libmesh_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "libmesh", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh):
    """Write a libMesh ``.xda`` (ASCII) or ``.xdr`` (XDR binary) mesh.

    The libMesh-1.8.0 layout, encoded by the extension; a trailing ``.gz`` or
    ``.bz2`` compresses it. Every cell with a libMesh type becomes a level-0
    element (polygons, polyhedra and Lagrange cells are dropped with a
    warning). Subdomain ids come from ``libmesh:subdomain``, else the cell
    regions; side regions become side sets, point regions node sets, and the
    ``:edge``/``:shellface<k>`` cell regions edge and shell-face sets.
    """
    if is_buffer(filename, "w"):
        _py_write(filename, mesh)
        return
    path = str(filename)
    lower = path.lower()
    compress = (
        gzip if lower.endswith(".gz") else bz2 if lower.endswith(".bz2") else None
    )
    if compress is None:
        try:
            _core.libmesh_write(path, mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "libmesh", "write", filename):
                raise
        _py_write(filename, mesh)
        return
    # The core writes the plain stream; Python compresses it.
    plain = path[: path.rfind(".")]
    fd, tmp = tempfile.mkstemp(
        suffix=os.path.splitext(plain)[1], dir=os.path.dirname(os.path.abspath(path))
    )
    os.close(fd)
    try:
        try:
            _core.libmesh_write(tmp, mesh)
        except Exception as exc:
            if not core_declined(exc, "libmesh", "write", filename):
                raise
            _py_write(tmp, mesh)
        with open(tmp, "rb") as src:
            data = src.read()
        with open(path, "wb") as dst:
            dst.write(
                gzip.compress(data, mtime=0) if compress is gzip else bz2.compress(data)
            )
    finally:
        os.remove(tmp)


register_format("libmesh", _EXTENSIONS, read, {"libmesh": write})

__all__ = ["read", "write"]
