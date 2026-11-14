from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from .._series_writer import CoreSeriesWriter
from .common import _gmsh_to_meshio_type as gmsh_to_meshioplusplus_type
from .common import _meshio_to_gmsh_type as meshioplusplus_to_gmsh_type
from .common import _validate_periodic
from .main import SeriesWriter as _PySeriesWriter
from .main import read as _py_read
from .main import write as _py_write


def read(filename, points_only=False, arrays=None, time_step: int = 0):
    """Read a Gmsh .msh file.

    Uses the C++ core for format versions 2.2, 4.0 and 4.1 (ASCII or binary),
    including periodic links, falling back to the reference Python reader for
    constructs the C++ reader doesn't handle.

    ``time_step`` selects one step of a `$NodeData`/`$ElementData` timeline
    (0 = first, negative counts from the end), resolved the same way the C
    API/Fortran/Julia/R/WASM surfaces do -- see
    :func:`meshioplusplus.gmsh.read`'s C++ counterpart, ``read_gmsh``. A
    non-default value forces the C++ path (the Python reference has no
    transient-step support) and re-raises rather than silently falling back.
    """
    # points_only/arrays reach the C++ reader, which skips the unwanted
    # <DataArray>/section bodies outright. The Python fallback below has no
    # selective support, so _helpers.read trims its result instead -- same
    # answer, just without the saving.
    if not is_buffer(filename, "r"):
        try:
            return _core.gmsh_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
            )
        except Exception as exc:
            if time_step:
                raise
            if not core_declined(exc, "gmsh", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh, fmt_version="4.1", binary=True, float_fmt=".16e"):
    """Write a Gmsh .msh file.

    Uses the C++ core for format versions 2.2 and 4.1 (ASCII or binary), including
    periodic links; otherwise falls back to the reference Python writer.
    """
    _validate_periodic(getattr(mesh, "gmsh_periodic", None), len(mesh.points))
    if float_fmt == ".16e" and not is_buffer(filename, "w"):
        if fmt_version == "2.2":
            try:
                _core.gmsh22_write(str(filename), mesh, binary)
                return
            except Exception as exc:
                if not core_declined(exc, "gmsh", "write", filename):
                    raise
        elif fmt_version == "4.1":
            try:
                # Bounding entities are signed entity tags, so they live in
                # cell_sets rather than on the C++ Mesh and are handed over
                # separately (the read path's GmshInfo channel, in reverse).
                _core.gmsh41_write(
                    str(filename),
                    mesh,
                    binary,
                    mesh.cell_sets.get("gmsh:bounding_entities"),
                )
                return
            except Exception as exc:
                if not core_declined(exc, "gmsh", "write", filename):
                    raise
    return _py_write(
        filename, mesh, fmt_version=fmt_version, binary=binary, float_fmt=float_fmt
    )


class SeriesWriter(CoreSeriesWriter):
    """Write a time series into one Gmsh file (2.2 or 4.1): the mesh with the first
    step, then each later step's ``$NodeData`` and ``$ElementData`` sections, with the
    step's time as their real tag and its position as their first integer tag (what
    ``read(..., time_step=k)`` selects by). Every step must have the first step's
    points and cells.

    >>> with meshioplusplus.gmsh.SeriesWriter("run.msh") as w:
    ...     for t, mesh in steps:
    ...         w.write(t, mesh)
    """

    def __init__(self, filename, fmt_version="4.1", binary=False):
        self._core = self._py = None
        if not is_buffer(filename, "w") and fmt_version in ("2.2", "4.1"):
            try:
                self._core = _core.GmshSeriesWriter(
                    str(filename), fmt_version == "4.1", binary
                )
                return
            except Exception as exc:
                if not core_declined(exc, "gmsh", "write", filename):
                    raise
        self._py = _PySeriesWriter(filename, fmt_version, binary)


register_format(
    "gmsh",
    [".msh"],
    read,
    {
        "gmsh22": lambda f, m, **kwargs: write(f, m, "2.2", **kwargs),
        "gmsh": lambda f, m, **kwargs: write(f, m, "4.1", **kwargs),
    },
)

__all__ = [
    "read",
    "write",
    "SeriesWriter",
    "gmsh_to_meshioplusplus_type",
    "meshioplusplus_to_gmsh_type",
]
