from .. import _core
from .._files import is_buffer
from .._helpers import register_format
from .common import _gmsh_to_meshio_type as gmsh_to_meshioplusplus_type
from .common import _meshio_to_gmsh_type as meshioplusplus_to_gmsh_type
from .main import read as _py_read
from .main import write as _py_write


def read(filename, points_only=False, arrays=None):
    """Read a Gmsh .msh file.

    Uses the C++ core for format version 2.2 (ascii or binary), falling back to
    the reference Python reader for versions 4.0/4.1, periodic meshes, and
    anything else the C++ reader doesn't handle.
    """
    # points_only/arrays reach the C++ reader, which skips the unwanted
    # <DataArray>/section bodies outright. The Python fallback below has no
    # selective support, so _helpers.read trims its result instead -- same
    # answer, just without the saving.
    if not is_buffer(filename, "r"):
        try:
            return _core.gmsh_read(
                str(filename), points_only=points_only, arrays=arrays
            )
        except Exception:
            pass
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
            except Exception:
                pass
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
            except Exception:
                pass
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
