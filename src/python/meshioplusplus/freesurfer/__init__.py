"""FreeSurfer surfaces (roadmap 2.4.4): C++ core for real file paths, the Python
reference otherwise."""

from .. import _core
from .._fallback import core_declined
from .._helpers import register_format
from ._freesurfer import read as _py_read
from ._freesurfer import write as _py_write


def read(filename, curv=None, annot=None):
    """Read a FreeSurfer triangle surface (``lh.pial``, ``rh.white``, ...).

    ``curv`` names curv files (``lh.thickness``...), each read as point data named
    after the file without its ``lh.``/``rh.`` prefix; ``annot`` names an annot file,
    read as the ``freesurfer:annot`` point data (structure index, -1 for none) and
    one point region per colour-table entry.
    """
    if isinstance(curv, (str, bytes)) or hasattr(curv, "__fspath__"):
        curv = [curv]
    try:
        return _core.freesurfer_read(
            str(filename),
            [str(c) for c in curv or []],
            "" if annot is None else str(annot),
        )
    except Exception as exc:
        if not core_declined(exc, "freesurfer", "read", filename):
            raise
    return _py_read(filename, curv=curv, annot=annot)


def write(filename, mesh):
    """Write one triangle block as a FreeSurfer triangle surface."""
    try:
        _core.freesurfer_write(str(filename), mesh)
        return
    except Exception as exc:
        if not core_declined(exc, "freesurfer", "write", filename):
            raise
    _py_write(filename, mesh)


register_format(
    "freesurfer",
    [".pial", ".white", ".inflated", ".orig", ".smoothwm", ".sphere"],
    read,
    {"freesurfer": write},
)

__all__ = ["read", "write"]
