from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._febio import read as _py_read
from ._febio import write as _py_write


def read(filename, points_only=False, arrays=None, lenient=False):
    """Read the mesh of an FEBio input file (``.feb``, febio_spec 2.5, 3.0 or 4.0).

    ``<Nodes>`` give the points and each ``<Elements>`` block a cell block and a
    cell region of its name; node sets, element sets and surfaces become point,
    cell and side regions, ``<Edge>``/``<DiscreteSet>`` line blocks, and
    ``<MeshData>`` node and element data (NaN outside its set). Materials,
    loads and steps are not read. ``lenient`` downgrades ``tet5``/``tet15``
    blocks to ``tetra``/``tetra10`` instead of failing.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.febio_read(
                str(filename), points_only=points_only, arrays=arrays, lenient=lenient
            )
        except Exception as exc:
            if not core_declined(exc, "febio", "read", filename):
                raise
    return _py_read(filename, points_only=points_only, arrays=arrays, lenient=lenient)


def write(filename, mesh):
    """Write an FEBio spec-4.0 ``.feb`` mesh.

    ``<Mesh>`` and ``<MeshDomains>`` plus a placeholder material per domain.
    2-D blocks lying on solid faces become ``<Surface>``s and line blocks
    ``<Edge>``s; regions become node sets, element sets and surfaces. Vertex
    cells and data arrays are dropped with a warning.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.febio_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "febio", "write", filename):
                raise
    return _py_write(filename, mesh)


register_format("febio", [".feb"], read, {"febio": write})

__all__ = ["read", "write"]
