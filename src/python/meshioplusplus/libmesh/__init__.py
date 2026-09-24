from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._libmesh import read as _py_read


def read(filename):
    """Read a libMesh ``.xda`` (ASCII) or ``.xdr`` (XDR binary) mesh.

    Active (leaf) elements become the cells, the subdomain id the
    ``libmesh:subdomain`` cell data and a cell region per subdomain (named by the
    file's subdomain map, else ``subdomain_<id>``). Side sets become side regions
    (``boundary_<id>`` when unnamed), node sets point regions (``nodeset_<id>``).
    The encoding is detected from the content. The C++ core reads the whole
    format; the Python reader is the fallback for buffers and for anything the
    C++ path raises on.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.libmesh_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "libmesh", "read", filename):
                raise
    return _py_read(filename)


register_format("libmesh", [".xda", ".xdr"], read, {})

__all__ = ["read"]
