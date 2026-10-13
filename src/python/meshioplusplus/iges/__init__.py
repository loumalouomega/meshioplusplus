"""IGES (``.igs``, ``.iges``): read through the optional OpenCASCADE backend.

Read-only, the twin of :func:`meshioplusplus.step.read`: a core built with
``-DMESHIOPLUSPLUS_WITH_OCCT=ON`` triangulates the file's shape and keeps its
NURBS faces, trimming loops and edges in ``mesh.info["brep"]``; without it,
:func:`read` raises ``NotImplementedError`` naming the flag.
"""

from .._helpers import register_format
from ..step import _read_cad


def read(filename, deflection=0.0, angular_deflection=0.5):
    """Read an IGES file as its triangulation, through OpenCASCADE.

    See :func:`meshioplusplus.step.read` for ``deflection`` and
    ``angular_deflection`` and what the mesh carries.

    :raises NotImplementedError: in a build without OpenCASCADE.
    """
    return _read_cad("iges", filename, deflection, angular_deflection)


register_format("iges", [".igs", ".iges"], read, {})

__all__ = ["read"]
