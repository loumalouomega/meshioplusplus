from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._mdpa import read as _py_read
from ._mdpa import write as _py_write


def _has_mdpa_extras(mesh):
    """Whether ``mesh`` carries MDPA data only the Python reference can write.

    ``misc_data`` (SubModelPart hierarchy, original Kratos ids, ``Mesh``
    blocks), ``geometries_block`` and the ``field_data`` conventions
    (``properties_<id>``/``table_<id>`` dicts, gmsh physical-group tuples) have
    no counterpart in the C++ ``Mesh``, so a mesh holding any of them keeps the
    pure-Python writer.
    """
    if getattr(mesh, "misc_data", None):
        return True
    if getattr(mesh, "geometries_block", None):
        return True
    return bool(getattr(mesh, "field_data", None))


def read(filename):
    """Read a Kratos MDPA mesh file.

    This deliberately stays on the pure-Python reference reader
    (``meshioplusplus.mdpa._mdpa``) rather than trying ``_core.mdpa_read``
    first, unlike the other format shims. The C++ core *does* implement the
    format (``meshioplusplus::read_mdpa``, reachable from the C API, Fortran,
    Julia, R, WebAssembly and the native CLI through the shared registry), but
    it necessarily returns a repo-standard ``Mesh``: ``cell_data`` keyed by
    array name with one array per cell block, no ``misc_data`` and no
    ``geometries_block``. The reference reader instead keys ``cell_data`` by
    cell type and carries MDPA's ``Tables``/``Properties``/``Geometries``/
    ``Mesh`` blocks and the SubModelPart hierarchy in ``mesh.misc_data`` — the
    documented behaviour of this module. Preferring the C++ reader here would
    silently change all of that, so it is reached explicitly
    (``_core.mdpa_read``) or through the flat bindings instead.
    """
    return _py_read(filename)


def write(filename, mesh, float_fmt=".16e", binary=False):
    """Write a Kratos MDPA mesh file (C++ core when possible, Python fallback).

    The C++ writer is used for real file paths and meshes carrying none of the
    MDPA-specific extras listed in :func:`_has_mdpa_extras`; anything else (and
    any failure in the C++ path) falls back to the pure-Python reference
    writer, which is the only one that can reproduce them.
    """
    if not binary and not is_buffer(filename, "w") and not _has_mdpa_extras(mesh):
        try:
            _core.mdpa_write(str(filename), mesh)
            return
        except Exception:
            pass
    return _py_write(filename, mesh, float_fmt=float_fmt, binary=binary)


register_format("mdpa", [".mdpa"], read, {"mdpa": write})

__all__ = ["read", "write"]
