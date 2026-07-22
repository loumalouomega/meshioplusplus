from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from ._unv import read as _py_read
from ._unv import write as _py_write


def read(filename):
    """Read an I-DEAS Universal file (C++ core; Python fallback for buffers)."""
    if not is_buffer(filename, "r"):
        try:
            return _core.unv_read(str(filename))
        except Exception:
            pass
    return _py_read(filename)


def write(filename, mesh, code_aster=False, node_dataset=2411):
    """Write an I-DEAS Universal file.

    The C++ core handles nodes, elements, field data (``point_data`` ->
    datasets 2414/55, ``cell_data`` -> 2414/57) and permanent groups
    (``point_sets``/``cell_sets`` -> dataset 2467); it falls back to the
    Python reference for buffer targets.  ``code_aster`` emits legacy field
    datasets 55/57 instead of 2414; ``node_dataset`` selects ``2411``
    (default) or ``781``.
    """
    if not is_buffer(filename, "w"):
        point_sets = dict(getattr(mesh, "point_sets", None) or {})
        cell_sets = {
            k: list(v) for k, v in (getattr(mesh, "cell_sets", None) or {}).items()
        }
        try:
            _core.unv_write(
                str(filename),
                mesh,
                point_sets,
                cell_sets,
                code_aster=code_aster,
                node_dataset=node_dataset,
            )
            return
        except Exception:
            pass
    return _py_write(filename, mesh, code_aster=code_aster, node_dataset=node_dataset)


register_format("unv", [".unv"], read, {"unv": write})

__all__ = ["read", "write"]
