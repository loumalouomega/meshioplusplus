from .. import _core
from .._exceptions import ReadError
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._anim import read as _py_read
from ._write import SeriesWriter as _PySeriesWriter
from ._write import numbered_names
from ._write import write as _py_write


def read(filename, time_step=0):
    """Read one OpenRadioss / Radioss animation file (``<run>A001``...).

    Each file is one state of the run: its time is ``field_data["meshio:time"]``,
    so a glob over the run's files (``runA*``) is a transient sequence. The
    points are the state's coordinates; nodal scalars and vectors (velocity,
    displacement, contact forces...) are point data; element scalars, tensors
    (xx yy zz xy yz zx) and 1-D force/moment sets are cell data, NaN where an
    element family has none. Parts are cell regions; ``radioss:part``,
    ``radioss:alive`` (0 once an element is deleted) and, when the file has them, ids, masses, materials and
    properties are cell data. The file holds one step: ``time_step`` may be 0
    or -1.
    """
    if time_step not in (0, -1):
        raise ReadError(
            f"Radioss animation: time_step {time_step} is out of range; "
            f"'{filename}' holds 1 step"
        )
    if not is_buffer(filename, "r"):
        try:
            return _core.radioss_anim_read(str(filename))
        except Exception as exc:
            if not core_declined(exc, "radioss_anim", "read", filename):
                raise
    return _py_read(filename)


def write(filename, mesh):
    """Write one OpenRadioss animation file: the state of ``mesh`` (its time is
    ``field_data["meshio:time"]``), big-endian, magic ``0x542C``. The name is
    ``<run>A001``-shaped, with no extension. Points, nodal scalars and vectors,
    element scalars, 6-wide tensors (xx yy zz xy yz zx), 9-wide 1-D force sets,
    parts (``radioss:part`` or the cell regions), ids, masses, materials and
    properties are written; vertices become SPH particles, lines 1-D elements,
    triangles and quads facets, and tetrahedra, pyramids, wedges and hexahedra
    bricks."""
    if not is_buffer(filename, "w"):
        try:
            return _core.radioss_anim_write(str(filename), mesh)
        except Exception as exc:
            if not core_declined(exc, "radioss_anim", "write", filename):
                raise
    return _py_write(filename, mesh)


class SeriesWriter:
    """One animation file per step, numbered from the number in ``filename``
    (``runA001``, ``runA002``...); ``write_sequence(path, steps,
    file_format="radioss_anim")`` uses it.

    >>> with meshioplusplus.radioss_anim.SeriesWriter("runA001") as w:
    ...     for time, mesh in steps:
    ...         w.write(time, mesh)
    """

    def __init__(self, filename):
        self._core = self._py = None
        self._names = numbered_names(filename)
        self.paths = []
        if not is_buffer(filename, "w"):
            try:
                self._core = _core.RadiossAnimSeriesWriter(str(filename))
            except Exception as exc:
                if not core_declined(exc, "radioss_anim", "write", filename):
                    raise
        if self._core is None:
            self._py = _PySeriesWriter(filename)

    def write(self, time, mesh):
        self.paths.append(next(self._names))
        if self._core is not None:
            self._core.write(float(time), mesh)
        else:
            self._py.write(time, mesh)

    def close(self):
        self._core = None
        if self._py is not None:
            self._py.close()
            self._py = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


# Found by its file name (a stem, then A and three or more digits) or its magic
# number, never an extension.
register_format("radioss_anim", [], read, {"radioss_anim": write})

__all__ = ["read", "write", "SeriesWriter"]
