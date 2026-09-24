from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._elmer import read as _py_read
from ._elmer import write as _py_write


def read(filename, points_only=False, arrays=None, piece=None, lenient=False):
    """Read an Elmer mesh directory (ElmerSolver's native mesh).

    ``filename`` is the mesh directory, a ``partitioning.N`` directory, or the
    ``mesh.header`` file inside one. Bulk and boundary elements become cell
    blocks; body and boundary ids become cell regions, named from ``mesh.names``
    (or ``body_<id>``/``boundary_<id>``) with the id as the tag. A partitioned
    mesh is merged, each cell's part in ``cell_data["partition:part"]``;
    ``piece`` reads one part alone. ``lenient`` skips element types with no
    meshio++ cell type instead of failing.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.elmer_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                piece=piece,
                lenient=lenient,
            )
        except Exception as exc:
            if not core_declined(exc, "elmer", "read", filename):
                raise
    return _py_read(
        filename, points_only=points_only, arrays=arrays, piece=piece, lenient=lenient
    )


def write(filename, mesh, halo: bool = False):
    """Write a serial Elmer mesh directory, creating it if needed.

    The cells of the highest dimension are the bulk elements, every
    lower-dimensional cell (and every side region's facets) a boundary element
    with its parents regenerated; regions give the body and boundary ids and
    names. Point regions and data arrays are dropped with a warning.

    A ``partition:part`` cell array also writes ElmerGrid's
    ``partitioning.<N>`` directory; ``halo=True`` adds ElmerGrid's ``-halo``
    layer to it (each bulk element with a whole side in another part copied
    there as ``id/owner``), which discontinuous Galerkin solvers need.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.elmer_write(str(filename), mesh, halo)
            return
        except Exception as exc:
            if not core_declined(exc, "elmer", "write", filename):
                raise
    return _py_write(filename, mesh, halo=halo)


# A directory, not a file: no extension maps to it; sniff_format finds it.
register_format("elmer", [], read, {"elmer": write})

__all__ = ["read", "write"]
