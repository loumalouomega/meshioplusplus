"""
I/O for the FLUX field file (``.dex``), the field companion to the FLUX mesh
(``.pf3``), used by Altair/CEDRAT FLUX and following FEconv
<https://github.com/victorsndvg/FEconv>.

A DEX file stores a single nodal field: a two-line header delimited by ``#``
giving the piece and field names and the ``NB_REAL``/``NB_COMP``/``NB_POINT``
counts, then one row per point holding the point's coordinates (x y z) followed
by its ``NB_COMP`` field values.  Read here as a geometry-less :class:`Mesh`
(no cells) whose ``points`` come from the coordinates and whose
``point_data[<field>]`` holds the values.
"""

import re

import numpy as np

from .._files import open_file
from .._mesh import Mesh

__all__ = ["read", "write"]

_DIM = 3  # DEX coordinates are always written as x y z


def _header_value(text, key, cast=str):
    m = re.search(rf"{key}\s*=\s*(\S+)", text)
    return cast(m.group(1)) if m else None


def read(filename):
    with open_file(filename, "r") as f:
        lines = f.read().splitlines()

    # header: the first two non-empty lines (the second ends with '#')
    header = []
    body_start = 0
    for i, ln in enumerate(lines):
        if ln.strip():
            header.append(ln)
        if len(header) == 2:
            body_start = i + 1
            break
    head_text = " ".join(header)
    field = _header_value(head_text, "FORMULA") or "dex:field"
    ncomp = _header_value(head_text, "NB_COMP", int) or 1
    npoint = _header_value(head_text, "NB_POINT", int) or 0

    rows = []
    for ln in lines[body_start:]:
        toks = ln.replace("D", "E").replace("d", "e").split()
        if toks:
            rows.append([float(t) for t in toks])
        if npoint and len(rows) >= npoint:
            break
    data = np.array(rows, dtype=float) if rows else np.empty((0, _DIM + ncomp))

    points = data[:, :_DIM] if data.shape[1] >= _DIM else data
    values = data[:, _DIM : _DIM + ncomp]
    point_data = {field: values[:, 0] if ncomp == 1 else values}
    return Mesh(points, [], point_data=point_data)


def write(filename, mesh, piece="PIECE", float_fmt=".16g"):
    points = mesh.points
    if points.shape[1] < _DIM:
        points = np.column_stack(
            [points, np.zeros((len(points), _DIM - points.shape[1]))]
        )
    pd = getattr(mesh, "point_data", None) or {}
    if not pd:
        raise ValueError("DEX write needs a nodal field in point_data")
    field, arr = next(iter(pd.items()))
    arr = np.asarray(arr, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    ncomp = arr.shape[1]
    n = len(points)

    with open_file(filename, "w") as f:
        f.write(f"# NAME = {piece} FORMULA = {field}\n")
        f.write(f"NB_REAL = 1 NB_COMP = {ncomp} NB_POINT = {n} #\n")
        for i in range(n):
            coords = " ".join(f"{x:{float_fmt}}" for x in points[i, :_DIM])
            vals = " ".join(f"{x:{float_fmt}}" for x in arr[i])
            f.write(f"{coords} {vals}\n")
