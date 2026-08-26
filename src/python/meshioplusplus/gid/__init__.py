import enum

from .. import _core
from .._helpers import register_format

#: ``field_data`` key prefix declaring an array's :class:`ResultType`.
#:
#: The full key is this prefix plus the array's own name, and its value is a
#: single integer :class:`ResultType`. ``field_data`` is global rather than
#: per-location, so one key covers an array of that name wherever it appears:
#: if the same name is in both ``point_data`` and ``cell_data``, the single
#: declaration applies to both and must be legal for both component counts.
RESULT_TYPE_PREFIX = "gid:result_type:"


class ResultType(enum.IntEnum):
    """What kind of quantity a GiD result array holds.

    meshio++'s ``Mesh`` cannot say "this array is a symmetric tensor" or "this
    array is complex" -- numpy's ``complex128`` has no meshio++ dtype -- so a
    caller declares it out of band::

        mesh.field_data[gid.RESULT_TYPE_PREFIX + "stress"] = [gid.ResultType.MATRIX]

    An array with no declaration keeps the historical inference (1 component ->
    ``Scalar``, 2 or 3 -> ``Vector``, anything else split into that many named
    scalars), so output for undeclared arrays is byte-identical to before this
    existed. A 6-component array is deliberately *not* inferred as ``MATRIX``:
    ``(n, 6)`` could equally be ``Matrix:6``, ``ComplexMatrix:3`` or
    ``ComplexVector:6``, so inferring would silently pick one meaning.

    Values are stored **verbatim in GiD's own component order**, which meshio++
    does not reinterpret. Two orders are worth knowing: ``MATRIX`` with 6
    components is ``Sxx Syy Szz Sxy Syz Sxz``, which is already meshio/VTK's
    symmetric-tensor order (so a stress tensor needs no permutation); and
    ``COMPLEX_VECTOR`` *interleaves* real and imaginary parts per component
    (``x_re x_im y_re y_im z_re z_im``) while ``COMPLEX_MATRIX`` *blocks* them
    (every real, then every imaginary) -- the same family, opposite
    conventions.

    Legal component counts, which the writer enforces (an illegal count is a
    :class:`meshioplusplus.WriteError` naming the array, never a silent
    fallback): ``SCALAR`` 1; ``VECTOR`` 2/3/4; ``MATRIX`` 3/6;
    ``PLAIN_DEFORMATION_MATRIX`` 4; ``MAIN_MATRIX`` 12; ``LOCAL_AXES`` 3;
    ``COMPLEX_SCALAR`` 2; ``COMPLEX_VECTOR`` 4/6; ``COMPLEX_MATRIX`` 6/12.
    """

    SCALAR = 0
    VECTOR = 1
    MATRIX = 2
    PLAIN_DEFORMATION_MATRIX = 3
    MAIN_MATRIX = 4
    LOCAL_AXES = 5
    COMPLEX_SCALAR = 6
    COMPLEX_VECTOR = 7
    COMPLEX_MATRIX = 8


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

__all__ = ["read", "write", "ResultType", "RESULT_TYPE_PREFIX"]
