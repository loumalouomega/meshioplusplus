import pathlib

from .. import _core, medit_sol
from .._exceptions import ReadError, WriteError
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._medit import read as _py_read
from ._medit import write as _py_write


def _is_binary_name(filename):
    return str(filename).endswith("b")


def _companion(filename):
    """``<stem>.sol`` beside an ascii mesh, ``<stem>.solb`` beside a binary one."""
    path = pathlib.Path(filename)
    return path.with_suffix(".solb" if _is_binary_name(filename) else ".sol")


def _sol_path(filename, sol):
    if sol is True:
        if is_buffer(filename, "r"):
            raise ReadError("Medit: sol=True needs a mesh path; pass the solution's")
        path = pathlib.Path(filename)
        for candidate in (path.with_suffix(".sol"), path.with_suffix(".solb")):
            if candidate.exists():
                return candidate
        raise ReadError(f"Medit: no {path.stem}.sol or {path.stem}.solb beside {path}")
    return sol


def read(filename, sol=None):
    """Read a Medit file: ``.mesh`` ascii, or ``.meshb`` binary (libMeshb GMF)
    when the name ends in ``b``. The C++ core reads both; the reference Python
    reader is the fallback.

    ``sol`` attaches a Medit solution (``.sol``/``.solb``) to the mesh: a path,
    or ``True`` for ``<stem>.sol`` (else ``<stem>.solb``) beside the mesh. Its
    vertex fields become ``point_data["medit:sol:<i>"]`` and its element fields
    ``cell_data["medit:sol:<i>"]``."""
    mesh = _read_mesh(filename)
    if sol is not None and sol is not False:
        mesh = medit_sol.attach(mesh, _sol_path(filename, sol))
    return mesh


def _read_mesh(filename):
    if not is_buffer(filename, "r"):
        try:
            if _is_binary_name(filename):
                return _core.medit_read_binary(str(filename))
            return _core.medit_read_ascii(str(filename))
        except Exception as exc:
            if not core_declined(exc, "medit", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh, float_fmt=".16e", sol=False):
    """Write a Medit file: ascii ``.mesh``, or binary ``.meshb`` (little-endian
    float64, version 3 or, for 64-bit indices, 4) when the name ends in ``b``.

    ``sol=True`` also writes the floating-point point and cell data as a Medit
    solution beside it: ``<stem>.sol`` for an ascii mesh, ``<stem>.solb`` for a
    binary one."""
    if sol:
        if is_buffer(filename, "w"):
            raise WriteError("Medit: sol=True needs a mesh path")
        medit_sol.write(_companion(filename), mesh, float_fmt=float_fmt)
    if not is_buffer(filename, "w") and (
        _is_binary_name(filename) or float_fmt == ".16e"
    ):
        try:
            if _is_binary_name(filename):
                _core.medit_write_binary(str(filename), mesh)
            else:
                _core.medit_write_ascii(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "medit", "write", filename):
                raise
    return _py_write(filename, mesh, float_fmt=float_fmt)


register_format("medit", [".mesh", ".meshb"], read, {"medit": write})

__all__ = ["read", "write"]
