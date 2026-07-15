from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._netgen import read as _py_read
from ._netgen import write as _py_write


def read(filename):
    """Read a Netgen .vol file (C++ core for the common path, Python fallback)."""
    if not is_buffer(filename, "r") and not str(filename).endswith(".vol.gz"):
        try:
            return _core.netgen_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, float_fmt=".16e"):
    """Write a Netgen .vol file (C++ core for the common path, Python fallback)."""
    # The C++ writer covers points + cells + the single integer cell index.
    # Identifications (stored in mesh.info) and field_data (codim material/bc
    # names) and the gzip container are left to the reference Python writer.
    info = getattr(mesh, "info", None)
    has_ident = isinstance(info, dict) and (
        info.get("netgen:identifications") is not None
    )
    if (
        not is_buffer(filename, "w")
        and not str(filename).endswith(".vol.gz")
        and not mesh.field_data
        and not has_ident
    ):
        try:
            _core.netgen_write(str(filename), mesh, float_fmt)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, float_fmt)


register_format("netgen", [".vol", ".vol.gz"], read, {"netgen": write})

__all__ = ["read", "write"]
