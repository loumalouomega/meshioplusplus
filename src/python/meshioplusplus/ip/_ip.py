"""
I/O for the ANSYS Fluent interpolation file (``.ip``), following FEconv
<https://github.com/victorsndvg/FEconv>.

An IP file stores one or more fields over a set of points.  Layout: version,
spatial dimension, point count, component count, then the component names, then
a section of all values for each coordinate (x, y, z), then a section of all
values for each field component -- in version 3 each section is wrapped in
``(``/``)``.  Only text files (versions 2 and 3) are supported; read here as a
geometry-less :class:`Mesh` (no cells) with ``points`` from the coordinate
sections and one ``point_data`` entry per field component.  Written as a
version-3 file.
"""

import numpy as np

from .._files import open_file
from .._mesh import Mesh

__all__ = ["read", "write"]


def read(filename):
    with open_file(filename, "r") as f:
        raw = f.read()
    lines = raw.splitlines()

    # header: first four non-empty lines are version, dim, npoint, ncomp
    ints = []
    idx = 0
    while len(ints) < 4 and idx < len(lines):
        s = lines[idx].strip()
        idx += 1
        if s:
            ints.append(int(float(s.split()[0])))
    dim, npoint, ncomp = ints[1], ints[2], ints[3]

    # next ncomp non-empty lines are the field names
    names = []
    while len(names) < ncomp and idx < len(lines):
        s = lines[idx].strip()
        idx += 1
        if s:
            names.append(s)

    # the rest is (dim + ncomp) sections of npoint reals each, column-major;
    # '(' and ')' are section delimiters -> treat as whitespace.
    rest = " ".join(lines[idx:]).replace("(", " ").replace(")", " ")
    rest = rest.replace("D", "E").replace("d", "e")
    flat = [float(t) for t in rest.split()]

    nsec = dim + ncomp
    need = nsec * npoint
    flat = flat[:need]
    sections = [flat[s * npoint : (s + 1) * npoint] for s in range(nsec)]

    coords = np.array(sections[:dim], dtype=float).T if dim else np.empty((npoint, 0))
    if coords.shape[0] != npoint:
        coords = np.empty((npoint, dim))
    point_data = {}
    for c in range(ncomp):
        vals = np.array(sections[dim + c], dtype=float)
        point_data[names[c]] = vals
    return Mesh(coords, [], point_data=point_data)


def write(filename, mesh, float_fmt=".16g"):
    points = mesh.points
    npoint = len(points)
    dim = points.shape[1] if points.ndim == 2 else 0
    pd = getattr(mesh, "point_data", None) or {}
    # flatten any multi-component point_data into scalar component columns
    names = []
    columns = []
    for name, arr in pd.items():
        arr = np.asarray(arr, dtype=float)
        if arr.ndim == 1:
            names.append(name)
            columns.append(arr)
        else:
            for c in range(arr.shape[1]):
                names.append(f"{name}_{c}")
                columns.append(arr[:, c])
    ncomp = len(columns)

    with open_file(filename, "w") as f:
        f.write("3\n")
        f.write(f"{dim}\n")
        f.write(f"{npoint}\n")
        f.write(f"{ncomp}\n")
        for name in names:
            f.write(f"{name}\n")
        for d in range(dim):
            f.write("(")
            f.write("\n".join(f"{x:{float_fmt}}" for x in points[:, d]))
            f.write("\n)\n")
        for col in columns:
            f.write("(")
            f.write("\n".join(f"{x:{float_fmt}}" for x in col))
            f.write("\n)\n")
