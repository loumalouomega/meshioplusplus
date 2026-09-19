from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from .common import _gmsh_to_meshio_type as gmsh_to_meshioplusplus_type
from .common import _meshio_to_gmsh_type as meshioplusplus_to_gmsh_type
from .main import read as _py_read
from .main import write as _py_write


def read(filename, points_only=False, arrays=None, time_step: int = 0):
    """Read a Gmsh .msh file.

    Uses the C++ core for format version 2.2 (ascii or binary), falling back to
    the reference Python reader for versions 4.0/4.1, periodic meshes, and
    anything else the C++ reader doesn't handle.

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

    Uses the C++ core for format versions 2.2 and 4.1 (ascii or binary) on
    non-periodic meshes; otherwise falls back to the reference Python writer.
    """
    if (
        float_fmt == ".16e"
        and getattr(mesh, "gmsh_periodic", None) is None
        and not is_buffer(filename, "w")
    ):
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
    "gmsh_to_meshioplusplus_type",
    "meshioplusplus_to_gmsh_type",
]
