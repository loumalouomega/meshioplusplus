import pathlib

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._mfem import read as _py_read
from ._mfem import write as _py_write


def _grid_function_items(grid_functions):
    if not grid_functions:
        return []
    if isinstance(grid_functions, dict):
        return [(str(k), str(v)) for k, v in grid_functions.items()]
    return [(pathlib.Path(p).stem, str(p)) for p in grid_functions]


def read(filename, grid_functions=None, piece=None):
    """Read an MFEM mesh (``.mesh``), optionally with grid functions (``.gf``).

    Elements and boundary elements become cell blocks with an ``mfem:attribute``
    cell data array and ``attribute_<n>``/``boundary_<n>`` cell regions; v1.3
    attribute sets become named cell regions. Order-2 ``H1`` nodes give quadratic
    cells whose points are MFEM's degrees of freedom; order 3 and up (nodal
    ``H1`` spaces and the legacy ``Cubic``) give VTK Lagrange cells. An ``MFEM
    NC mesh`` is read as its leaf elements.

    ``grid_functions`` is a ``{name: path}`` dict, or a list of paths named by
    their stem. ``H1`` fields become point data (an order-2 field on a linear
    mesh makes the cells quadratic, a higher one VTK Lagrange) and ``L2``
    order-0 fields cell data.

    A rank file of a parallel mesh (``<prefix>.000000``, ...) is read with its
    siblings as one mesh, ``cell_data["partition:part"]`` holding each cell's
    rank; ``piece=k`` reads rank ``k`` alone. A parallel mesh's grid function is
    named by one of its rank files, ``<name>.000000``.
    """
    items = _grid_function_items(grid_functions)
    if not is_buffer(filename, "r"):
        try:
            return _core.mfem_read(str(filename), items, piece)
        except Exception as exc:
            if not core_declined(exc, "mfem", "read", filename):
                raise
    return _py_read(filename, dict(items) if items else None, piece=piece)


def write(filename, mesh, grid_functions=False):
    """Write an MFEM mesh (``.mesh``).

    The highest-dimensional cells are the elements and the cells one dimension
    lower the boundary; side regions add boundary elements. Attributes come from
    ``mfem:attribute``, else cell regions, else 1; named regions become v1.3
    attribute sets. Quadratic cells are written as an ``H1_<d>D_P2`` nodes space.

    With ``grid_functions=True`` every point data array is also written as an
    ``H1`` grid function ``<stem>.<name>.gf`` next to the mesh (and cell data as
    ``L2`` order 0); otherwise data arrays are dropped with a warning.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.mfem_write(str(filename), mesh, bool(grid_functions))
            return
        except Exception as exc:
            if not core_declined(exc, "mfem", "write", filename):
                raise
    return _py_write(filename, mesh, grid_functions=grid_functions)


register_format("mfem", [".mesh"], read, {"mfem": write})

__all__ = ["read", "write"]
