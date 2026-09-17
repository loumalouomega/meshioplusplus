from .. import _core
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


def read(filename, region=""):
    """Read an OpenFOAM polyMesh case (C++ core, Python fallback).

    ``region`` selects one region of a multi-region case (v11.4.0, roadmap §1
    tier B2): ``<case>/constant/<region>/polyMesh``. The Python fallback
    reader has no multi-region concept at all, so a ``region`` request, or a
    case the compiled core recognised as multi-region (and asked for one),
    re-raises rather than silently falling back to a worse, generic
    "no polyMesh" error from a reader that never looked for one.
    """
    try:
        return _core.openfoam_read(str(filename), region)
    except Exception as exc:
        if region or "multi-region" in str(exc):
            raise
    return _py_read(filename)


def write(filename, mesh):
    """Write an OpenFOAM polyMesh case.

    ``filename`` may be a ``.foam`` marker file, a directory named ``polyMesh``,
    or any other directory taken as the case root; ``<case>/constant/polyMesh/``
    is created as needed. This is the only meshio++ writer that produces a
    directory rather than a file.

    Patch names and types are taken from ``mesh.cell_tags`` and
    ``mesh.openfoam_patch_types`` when present (both are set by :func:`read`); a
    mesh carrying neither -- anything converted from another format -- gets a
    single ``defaultFaces`` patch, which is what ``blockMesh`` itself produces.

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
    )


# Advertise the writer only when it can actually run, so `meshioplusplus.write`
# does not offer a format that always raises.
register_format("openfoam", [".foam"], read, {"openfoam": write} if _HAS_WRITE else {})

__all__ = ["read", "write"]
