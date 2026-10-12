"""STEP (ISO 10303-21, AP214) writer: a closed surface as a faceted B-rep solid.

Write-only: reading STEP means a CAD kernel (roadmap §6.2.1).
"""

from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._step import write as _py_write


def write(filename, mesh):
    """Write a closed surface (or a volume mesh's skin) as a STEP faceted solid.

    C++ core for real file paths; the byte-identical Python twin otherwise.
    """
    if not is_buffer(filename, "w"):
        try:
            _core.step_write(str(filename), mesh)
            return
        except Exception as exc:
            if not core_declined(exc, "step", "write", filename):
                raise
    return _py_write(filename, mesh)


register_format("step", [".step", ".stp"], None, {"step": write})

__all__ = ["write"]
