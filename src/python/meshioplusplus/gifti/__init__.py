from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._gifti import read as _py_read  # noqa: F401
from ._gifti import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "gifti",
    """Read a GIFTI surface: the POINTSET, the TRIANGLE array, every per-vertex array as
point data and the LabelTable as point regions (C++ core for real file paths, Python
fallback).""",
    """Write a triangle surface and its point data as GIFTI, little-endian Base64Binary
(C++ core for real file paths, Python fallback).""",
)


register_format("gifti", [".gii"], read, {"gifti": write})

__all__ = ["read", "write"]
