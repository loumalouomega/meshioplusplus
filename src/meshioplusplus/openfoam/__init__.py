from .. import _core
from .._helpers import register_format
from ._openfoam import read as _py_read


def read(filename):
    """Read an OpenFOAM polyMesh case (C++ core, Python fallback)."""
    try:
        return _core.openfoam_read(str(filename))
    except Exception:
        pass
    return _py_read(filename)


register_format("openfoam", [".foam"], read, {})

__all__ = ["read"]
