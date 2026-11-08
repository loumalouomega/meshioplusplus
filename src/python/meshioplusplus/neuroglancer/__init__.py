from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._neuroglancer import read as _py_read  # noqa: F401
from ._neuroglancer import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "neuroglancer",
    """Read a Neuroglancer precomputed mesh (C++ core for file paths, Python fallback).""",
    """Write a Neuroglancer precomputed mesh (C++ core for file paths, Python fallback).""",
)


register_format("neuroglancer", [], read, {"neuroglancer": write})

__all__ = ["read", "write"]
