from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._abaqus_fil import read as _py_read
from ._abaqus_fil import time_values as _py_time_values


def read(filename, points_only=False, arrays=None, time_step=0):
    """Read an Abaqus results file (``.fil``), binary or ASCII.

    Nodes (1901) and elements (1900/1990) become the mesh, with their labels as
    ``abaqus:id`` point and cell data; node and element sets (1931-1934) become
    point and cell regions, long names resolved through the 1940 labels. Every
    increment is a step: ``time_step`` selects one (0 = first, negative counts
    from the end), its total time is ``field_data["meshio:time"]``, and its
    nodal and element records become point and cell data named by their Abaqus
    identifier (``U``, ``RF``, ``S``, ``E``, ...).
    """
    if not is_buffer(filename, "r"):
        try:
            return _core.abaqus_fil_read(str(filename), points_only, arrays, time_step)
        except Exception as exc:
            if not core_declined(exc, "abaqus_fil", "read", filename):
                raise
    return _py_read(
        filename, points_only=points_only, arrays=arrays, time_step=time_step
    )


def time_values(filename):
    """The total time of every increment, in file order."""
    if not is_buffer(filename, "r"):
        try:
            return list(_core.abaqus_fil_time_values(str(filename)))
        except Exception as exc:
            if not core_declined(exc, "abaqus_fil", "read", filename):
                raise
    return _py_time_values(filename)


register_format("abaqus_fil", [".fil"], read, {})

__all__ = ["read", "time_values"]
