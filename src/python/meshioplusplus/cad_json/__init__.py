"""Kratos CAD JSON (``*.cad.json``): the NURBS B-rep of Kratos Multiphysics' IGA.

Faces, their trimming loops, edges and vertices are read into
``mesh.info["brep"]`` and sampled into Lagrange cells; the writer writes the
B-rep back (or, for a plain mesh, its quads as bilinear faces). C++ core for
real file paths; the Python reference otherwise. See doc/formats/cad_json.md.
"""

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._cad_json import from_mesh, parse_brep
from ._cad_json import read as _py_read
from ._cad_json import read_brep, render
from ._cad_json import sample as _py_sample
from ._cad_json import write as _py_write


def read(filename, subdivisions=0, trim_samples=16, edges=True, trim=True):
    """Read a Kratos CAD JSON file as its sampled mesh.

    Every knot span of every face becomes one Lagrange quadrilateral of the
    face's degree (``quad``, ``quad9``, else ``VTK_LAGRANGE_QUADRILATERAL``),
    ``subdivisions`` times finer per direction (``0``, the default, refines only
    a trimmed face, to about 16 cells per direction); a trimmed face keeps the cells
    whose parametric centre lies inside its boundary loops (sampled at
    ``trim_samples`` points per knot span; ``trim=False`` keeps every span).
    Edges with a 3-D curve become Lagrange line cells (``edges=False`` skips
    them). Cell data ``cad:kind``, ``cad:index``, ``cad:id``, ``cad:span`` and
    ``cad:trimmed`` link each cell back to ``mesh.info["brep"]``, which holds
    the exact faces, loops, edges and vertices. The ``Surfaces``/``Lines``
    dialect of Kratos's ``ImportNurbsSbmModeler`` is read too.
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.cad_json_read(
                str(filename), subdivisions, trim_samples, edges, trim
            )
        except Exception as exc:
            if not core_declined(exc, "cad_json", "read", filename):
                raise
    return _py_read(
        filename,
        subdivisions=subdivisions,
        trim_samples=trim_samples,
        edges=edges,
        trim=trim,
    )


def write(filename, mesh, brep=None):
    """Write Kratos CAD JSON.

    ``brep`` (default ``mesh.info["brep"]``) is written exactly, in the layout
    of Kratos's own ``CadJsonOutput``; without one, every ``quad`` cell becomes
    an untrimmed bilinear face and every ``line`` cell an edge, other cells
    being dropped with a warning.
    """
    if brep is None and isinstance(getattr(mesh, "info", None), dict):
        brep = mesh.info.get("brep")
    if not is_buffer(filename, "w"):
        try:
            _core.cad_json_write(str(filename), mesh, brep)
            return
        except Exception as exc:
            if not core_declined(exc, "cad_json", "write", filename):
                raise
    return _py_write(filename, mesh, brep=brep)


def sample(brep, subdivisions=0, trim_samples=16, edges=True, trim=True):
    """Sample a ``brep`` dict (``mesh.info["brep"]``) into a mesh, as ``read`` does."""
    try:
        return _core.brep_sample(brep, subdivisions, trim_samples, edges, trim)
    except Exception as exc:
        if not core_declined(exc, "cad_json", "read", None):
            raise
    return _py_sample(
        brep,
        subdivisions=subdivisions,
        trim_samples=trim_samples,
        edges=edges,
        trim=trim,
    )


register_format("cad_json", [".cad.json"], read, {"cad_json": write})

__all__ = ["read", "write", "sample", "read_brep", "parse_brep", "render", "from_mesh"]
