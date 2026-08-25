from .. import _core
from .._helpers import register_format


def write(filename, mesh, mode="auto", analysis_name="meshio++", step=1.0):
    """Write a mesh as a GiD postprocess file (write-only).

    ``mode`` is one of ``"auto"`` (infer the flavour from ``filename``'s
    extension: ``.post.bin`` -> binary, ``.post.h5`` -> hdf5, anything else,
    including ``.post.msh``/``.post.res`` -> ascii), ``"ascii"`` (two
    sibling files, ``<stem>.post.msh`` + ``<stem>.post.res``), ``"binary"``
    (one deflated ``<stem>.post.bin``), or ``"hdf5"`` (one
    ``<stem>.post.h5``, needing a build with ``MESHIOPLUSPLUS_WITH_HDF5=ON``
    in addition to gidpost itself).

    ``analysis_name``/``step`` are the GiD "analysis name" and time/load
    step every written result is grouped under.

    There is deliberately **no pure-Python fallback** for this format:
    gidpost (the vendored C library this writer is built on) cannot be
    cheaply reimplemented, and a second implementation of its node-ordering
    permutations risks silently disagreeing with the first. A build without
    gidpost (``-DMESHIOPLUSPLUS_WITH_GIDPOST=OFF`` or a missing zlib) still
    exposes this function -- it raises a :class:`meshioplusplus.WriteError`
    naming the missing build flags rather than falling through to another
    format, which is why this writer stays registered unconditionally
    (deregistering it, the way the write-optional ``openfoam`` shim does,
    would let ``.post.msh`` fall through to a generic "unknown format"
    error instead of this format's specific, actionable one).
    """
    _core.gid_write(str(filename), mesh, mode, analysis_name, step)


# Registered unconditionally (never gated on a "_HAS_WRITE" probe): see
# write()'s docstring for why. Reading is out of scope -- gidpost's own
# public API has zero read functions (doc/roadmap.md section 1's remaining
# item) -- so there is no reader here and no `read` in __all__.
register_format(
    "gid", [".post.msh", ".post.res", ".post.bin", ".post.h5"], None, {"gid": write}
)

__all__ = ["write"]
