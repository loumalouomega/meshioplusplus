from .. import _core
from .._helpers import register_format


def read(filename, time_step=0):
    """Read a GiD postprocess file.

    Accepts any of the four ``.post.*`` spellings. The flavour is resolved from
    the extension and confirmed against the leading bytes, so a gzipped
    ``.post.msh`` (gidpost's ``GiD_PostAsciiZipped``, which is the same ASCII
    text through ``gzprintf``) reads correctly too.

    For the ASCII flavour the geometry file ``<stem>.post.msh`` is mandatory
    and the results file ``<stem>.post.res`` is **optional** -- a mesh with no
    results reads back as geometry only. Passing the ``.post.res`` path
    directly derives and reads the ``.post.msh``; *its* absence is an error,
    since results alone carry no geometry.

    ``time_step`` selects one step of a multi-step results file (0 = the first,
    negative counts from the end); it is honoured natively, because no
    caller-side filter can recover a step that was never decoded.

    Unlike :func:`write`, reading needs **no gidpost and no zlib** for the
    ASCII flavour, so it works in builds that cannot write GiD at all -- see
    the module docstring in ``gid.hpp``.
    """
    return _core.gid_read(str(filename), time_step)


def write(filename, mesh, mode="auto", analysis_name="meshio++", step=1.0):
    """Write a mesh as a GiD postprocess file.

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
# write()'s docstring for why.
#
# `read` is passed directly rather than wrapped in the usual
# try-C++/except-fall-back-to-Python shim, because there is no Python engine to
# fall back to. Swallowing the exception would hide a core gap from the only
# tests most people run -- the hazard CLAUDE.md records for gmsh.
register_format(
    "gid", [".post.msh", ".post.res", ".post.bin", ".post.h5"], read, {"gid": write}
)

__all__ = ["read", "write"]
