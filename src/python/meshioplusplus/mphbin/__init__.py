from .. import _core  # noqa: F401 (the shims look it up at call time)
from .._helpers import register_format
from .._shims import plain_shims
from ..mphtxt._mphtxt import read_binary as _py_read  # noqa: F401
from ..mphtxt._mphtxt import write_binary as _py_write  # noqa: F401

read, write = plain_shims(
    __name__,
    "mphbin",
    """Read a COMSOL binary mesh (.mphbin) file: the content of a .mphtxt file,
    little-endian int32/float64, strings as one int32 per character.

    Uses the C++ core for real file paths, the Python reference otherwise.
    """,
    """Write a COMSOL binary mesh (.mphbin) file; see ``mphtxt.write``.

    The format has no comment slot, so no provenance is written.
    """,
    read_mode="rb",
    write_mode="wb",
)


register_format("mphbin", [".mphbin"], read, {"mphbin": write})

__all__ = ["read", "write"]
