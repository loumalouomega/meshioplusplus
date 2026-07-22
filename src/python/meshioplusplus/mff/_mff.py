"""
I/O for the Modulef Formatted Field (``.mff``) format, the field companion to
the Modulef Formatted Mesh (``.mfm``), following FEconv
<https://github.com/victorsndvg/FEconv>.

An MFF file stores a *single* field over the whole mesh as an integer value
count followed by a flat list of double-precision floats.  It carries no
geometry and no component/location metadata: the value count is a multiple of
the number of nodes (or elements) of the companion mesh, and that ratio is the
number of field components.  Read standalone here, the values are exposed as a
geometry-less :class:`Mesh` (no cells, ``points`` with zero columns) carrying
``point_data["mff:field"]``; scalar values round-trip exactly, but the
component count and node-vs-element location cannot be recovered without the
companion mesh.
"""

import numpy as np

from .._files import open_file
from .._mesh import Mesh

__all__ = ["read", "write"]


def read(filename):
    with open_file(filename, "r") as f:
        tokens = f.read().replace("D", "E").replace("d", "e").split()
    if not tokens:
        return Mesh(np.empty((0, 0)), [])
    count = int(tokens[0])
    values = np.array([float(t) for t in tokens[1 : 1 + count]], dtype=float)
    return Mesh(
        np.empty((len(values), 0)),
        [],
        point_data={"mff:field": values},
    )


def _first_field(mesh):
    """Pick the flat value vector to write: first point_data, else cell_data."""
    for arr in (getattr(mesh, "point_data", None) or {}).values():
        return np.asarray(arr, dtype=float).reshape(-1)
    for name, blks in (getattr(mesh, "cell_data", None) or {}).items():
        if name == "unv:pid":
            continue
        return np.concatenate([np.asarray(b, dtype=float).reshape(-1) for b in blks])
    return np.empty(0, dtype=float)


def write(filename, mesh, float_fmt=".16e"):
    values = _first_field(mesh)
    with open_file(filename, "w") as f:
        f.write(f"{len(values)}\n")
        for v in values:
            f.write(f"{v:{float_fmt}}\n")
