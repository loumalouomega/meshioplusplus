from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ._triangle import read as _py_read  # noqa: F401
from ._triangle import write as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "triangle",
    """Read Triangle files (.node/.ele pair or .poly; C++ core for real file
    paths, Python fallback).""",
    """Write Triangle files (.node/.ele pair or .poly; C++ core for real file
    paths, Python fallback).""",
)


# .node/.ele are also claimed by tetgen (registered first): the dispatcher
# falls through to this format when tetgen rejects a 2D file.
register_format("triangle", [".node", ".ele", ".poly"], read, {"triangle": write})

__all__ = ["read", "write"]
