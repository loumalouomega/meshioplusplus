from .. import _core
from .._fallback import core_declined
from .._helpers import register_format
from ._openfoam import read as _py_read

# The writer lives only in the compiled core. That is deliberate rather than an
# omission: a pure-Python twin would have to re-implement the per-cell winding
# repair, which is a discrete branch on the sign of an enclosed volume, and two
# implementations of such a branch can land on opposite sides for a
# near-degenerate cell and then diverge macroscopically -- the same reasoning
# that keeps `smooth`'s inversion guard out of its numpy fallback. It is also
# dead code in practice: `openfoam_write` needs no optional dependency, so it is
# present in every wheel and in every source build that has Python bindings at
# all, including the Windows CI leg that runs the format fallbacks.
_HAS_WRITE = hasattr(_core, "openfoam_write")


def read(filename, region="", points_only=False, arrays=None, time_step=0):
    """Read an OpenFOAM polyMesh case (C++ core, Python fallback).

    ``region`` selects one region of a multi-region case (v11.4.0, roadmap §1
    tier B2): ``<case>/constant/<region>/polyMesh``. ``time_step``/``arrays``
    (same tier) select a time directory (``<case>/<time>/``) and which of its
    fields to attach as ``point_data``/``cell_data`` -- also reachable through
    ``meshioplusplus.read(path, time_step=..., arrays=...)``, since this
    shim's signature is what lets the generic dispatcher forward them. The
    Python fallback reader has none of this -- no multi-region concept, no
    field reading at all -- so a ``region``/non-default ``time_step``
    request, or a case the compiled core recognised as multi-region, re-raises
    rather than silently falling back to a worse or plain-wrong result.
    """
    try:
        return _core.openfoam_read(
            str(filename), region, points_only, arrays, time_step
        )
    except Exception as exc:
        if region or time_step or "multi-region" in str(exc):
            raise
        if not core_declined(exc, "openfoam", "read", filename):
            raise
    return _py_read(filename)


def write(filename, mesh, binary=False, label_bits=32, scalar_bits=64):
    """Write an OpenFOAM polyMesh case.

    ``filename`` may be a ``.foam`` marker file, a directory named ``polyMesh``,
    or any other directory taken as the case root; ``<case>/constant/polyMesh/``
    is created as needed. This is the only meshio++ writer that produces a
    directory rather than a file.

    Patch names and types are taken from ``mesh.cell_tags`` and
    ``mesh.openfoam_patch_types`` when present (both are set by :func:`read`); a
    mesh carrying neither -- anything converted from another format -- gets a
    single ``defaultFaces`` patch, which is what ``blockMesh`` itself produces.

    ``binary`` (v15.5.0, roadmap §1.1) writes ``points``/``owner``/``neighbour``/
    zone label lists as raw little-endian bytes and ``faces`` the same
    length-prefixed-per-face way the reader already understood; ASCII (the
    default) is unchanged. ``label_bits``/``scalar_bits`` (32 or 64) set the
    integer/floating-point width and are only meaningful with ``binary=True`` --
    reachable from Python and C++ only, not from the flat bindings or either
    CLI's ``--binary``/``--ascii`` (which always use the 32/64 default). A
    binary request on a big-endian host raises by name.

    Note there is deliberately no Python fallback: a swallowed write error would
    mean a different implementation silently produced a different case.
    """
    if not _HAS_WRITE:
        raise RuntimeError(
            "meshio++: the OpenFOAM writer requires the compiled core "
            "(meshioplusplus._core.openfoam_write), which this installation "
            "does not provide"
        )
    _core.openfoam_write(
        str(filename),
        mesh,
        getattr(mesh, "cell_tags", None) or {},
        getattr(mesh, "openfoam_patch_types", None) or {},
        binary,
        label_bits,
        scalar_bits,
    )


# Advertise the writer only when it can actually run, so `meshioplusplus.write`
# does not offer a format that always raises.
register_format("openfoam", [".foam"], read, {"openfoam": write} if _HAS_WRITE else {})

__all__ = ["read", "write"]
