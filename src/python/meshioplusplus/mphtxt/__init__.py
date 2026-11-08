from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._mphtxt import read as _py_read  # noqa: F401
from ._mphtxt import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "mphtxt",
    """Read a COMSOL text mesh (.mphtxt) file.

    Uses the C++ core for real file paths, the Python reference otherwise.
    Selections become ``mesh.regions``; entity indices ``cell_data["mphtxt:geom"]``.
    """,
    """Write a COMSOL text mesh (.mphtxt) file (Mesh version 4).

    Uses the C++ core for real file paths, the Python reference otherwise. Cell
    regions that are unions of whole geometric entities become Selections.
    """,
)


register_format("mphtxt", [".mphtxt"], read, {"mphtxt": write})

__all__ = ["read", "write"]
