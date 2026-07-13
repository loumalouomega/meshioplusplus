from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._med import read as _py_read
from ._med import write as _py_write

_HAS_HDF5 = getattr(_core, "__has_hdf5__", False)


def read(filename):
    """Read a MED file (C++ core when built with HDF5, Python/h5py fallback)."""
    if _HAS_HDF5 and not is_buffer(filename, "r"):
        try:
            return _core.med_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh):
    """Write a MED file (C++ core when built with HDF5, Python/h5py fallback)."""
    if _HAS_HDF5 and not is_buffer(filename, "w"):
        point_tags = getattr(mesh, "point_tags", None) or {}
        cell_tags = getattr(mesh, "cell_tags", None) or {}
        med_nom = mesh.field_data.get("med:nom", [])
        try:
            _core.med_write(
                str(filename), mesh, dict(point_tags), dict(cell_tags), list(med_nom)
            )
            return
        except Exception:
            pass
    return _py_write(filename, mesh)


register_format("med", [".med"], read, {"med": write})

__all__ = ["read", "write"]
