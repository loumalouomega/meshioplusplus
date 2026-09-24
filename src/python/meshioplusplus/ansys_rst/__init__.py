from .. import _core
from .._fallback import core_declined
from .._files import is_buffer
from .._helpers import register_format
from ._ansys_rst import read as _py_read
from ._ansys_rst import time_values as _py_time_values


def read(filename, points_only=False, arrays=None, time_step=0, lenient=False):
    """Read one result set of an Ansys MAPDL results file (``.rst``, ``.rth``).

    The geometry records become cells exactly as a ``.cdb`` deck's elements do
    (degenerate shapes resolved, ``ansys:*`` cell data, components as point and
    cell regions). ``time_step`` picks the result set (negative counts from the
    end); its time, or frequency in a modal analysis, is
    ``field_data["meshio:time"]``, with ``ansys:load_step``, ``ansys:substep``
    and ``ansys:cumulative``. The set's nodal DOF solution is point data: ``U``,
    ``ROT``, ``A`` and ``V`` as vectors rotated to the global axes, other DOFs
    (``TEMP``, ``PRES`` ...) as scalars; NaN where a node has no value.
    Reaction forces are ``RF`` and ``RMOM`` (global axes) and ``RF_<DOF>``.
    Element nodal stresses and strains (``S``, ``EPEL``, ``EPPL``, ``EPCR``,
    ``EPTH``: ``xx yy zz xy yz xz``, rotated to the global axes; a layered
    shell's top surface as ``<name>@top``) are cell data per element node,
    ``(cells, nodes * components)`` flattened point-major with the layout in
    ``field_data["ansys:layout:<name>"]``, and point data averaged over the
    elements at each corner node; element nodal forces are cell data ``ENF``. The main file of a distributed solve
    (``file0.rst``) reads its partial files with it. ``lenient`` skips elements
    whose type has no meshio++ cell.
    """
    return _read(filename, points_only, arrays, time_step, lenient, False)


def read_cyclic(filename, points_only=False, arrays=None, time_step=0, lenient=False):
    """Read one result set of a static cyclic-symmetry model as the full rotor.

    The base sector (element numbers up to the model's ``csEls``) is repeated
    round the cyclic axis, its nodal and element results rotated with it;
    ``ansys:sector`` numbers the copies. Coincident nodes on the sector
    boundaries are not merged. Otherwise as :func:`read`.
    """
    return _read(filename, points_only, arrays, time_step, lenient, True)


def _read(filename, points_only, arrays, time_step, lenient, cyclic):
    fmt = "ansys_rst_cyclic" if cyclic else "ansys_rst"
    if not is_buffer(filename, "r"):
        try:
            return _core.ansys_rst_read(
                str(filename),
                points_only=points_only,
                arrays=arrays,
                time_step=time_step,
                lenient=lenient,
                cyclic=cyclic,
            )
        except Exception as exc:
            if not core_declined(exc, fmt, "read", filename):
                raise
    return _py_read(
        filename,
        points_only=points_only,
        arrays=arrays,
        time_step=time_step,
        lenient=lenient,
        cyclic=cyclic,
    )


def time_values(filename):
    """The time (or frequency) of every result set."""
    try:
        return list(_core.ansys_rst_time_values(str(filename)))
    except Exception as exc:
        if not core_declined(exc, "ansys_rst", "read", filename):
            raise
    return _py_time_values(filename)


register_format("ansys_rst", [".rst", ".rth"], read, {})
register_format("ansys_rst_cyclic", [], read_cyclic, {})

__all__ = ["read", "read_cyclic", "time_values"]
