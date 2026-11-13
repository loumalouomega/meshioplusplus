"""PLOT3D multi-block structured grids (roadmap 2.4.2): C++ core for real file paths,
the Python reference otherwise."""

from .. import _core
from .._fallback import core_declined
from .._helpers import register_format
from ._plot3d import _option
from ._plot3d import read as _py_read
from ._plot3d import write as _py_write


def read(
    filename,
    solution=None,
    function=None,
    precision=None,
    multiblock=None,
    iblank=None,
    dimensions=None,
):
    """Read a PLOT3D grid, its variant inferred from its size.

    ASCII, Fortran unformatted or raw binary; single or multi-block; 2-D or 3-D;
    single or double precision; with or without ``iblank``; either byte order. A
    file that reads several ways names them: pick one with ``precision``
    (``"single"``/``"double"``), ``multiblock``, ``iblank`` or ``dimensions`` (2/3).
    Each block is a cell block and a ``block_<b>`` cell region; the dimensions go to
    ``field_data["plot3d:dims"]`` and ``iblank`` to the ``plot3d:iblank`` point data.
    ``solution=`` reads a ``.q`` file (``density``, ``momentum``, ``energy``, and
    ``field_data["plot3d:q_reference"]``), ``function=`` a ``.f`` file
    (``plot3d:f``).
    """
    try:
        w = _option(precision, "precision", {"single": 4, "double": 8, 4: 4, 8: 8})
        dims = _option(dimensions, "dimensions", {2: 2, 3: 3})
        return _core.plot3d_read(
            str(filename),
            "" if solution is None else str(solution),
            "" if function is None else str(function),
            w or 0,
            -1 if multiblock is None else int(bool(multiblock)),
            -1 if iblank is None else int(bool(iblank)),
            dims or 0,
        )
    except Exception as exc:
        if not core_declined(exc, "plot3d", "read", filename):
            raise
    return _py_read(
        filename,
        solution=solution,
        function=function,
        precision=precision,
        multiblock=multiblock,
        iblank=iblank,
        dimensions=dimensions,
    )


def write(filename, mesh):
    """Write a mesh carrying ``plot3d:dims`` as a Fortran unformatted multi-block
    double-precision grid (with ``iblank`` when the mesh has ``plot3d:iblank``)."""
    try:
        _core.plot3d_write(str(filename), mesh)
        return
    except Exception as exc:
        if not core_declined(exc, "plot3d", "write", filename):
            raise
    _py_write(filename, mesh)


register_format("plot3d", [".p3d", ".p3dfmt", ".x", ".g"], read, {"plot3d": write})

__all__ = ["read", "write"]
