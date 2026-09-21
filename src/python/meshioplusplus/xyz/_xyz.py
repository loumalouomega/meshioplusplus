"""
I/O for headerless ASCII point clouds (``.xyz``, ``.xyzn``, ``.xyzrgb``, ``.asc``,
``.pts``, ``.txt``): one point per line, columns separated by whitespace, commas or
semicolons.

XYZ is a convention, not a specification, so the column meaning is resolved in
this order: an explicit ``columns=`` list; a header comment naming the columns
(``# x y z nx ny nz``, what the writer emits); the column count and the file
extension; the value ranges (six columns: unit vectors are normals, byte-valued
integers are colours). Anything still ambiguous is an error that asks for
``columns=`` rather than a guess. Chemistry XYZ (atom count, comment, ``element x
y z``) shares the extension and is refused by name.
"""

import re

import numpy as np

from .. import _provenance
from .._common import warn
from .._exceptions import ReadError
from .._files import open_file
from .._mesh import Mesh

_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_ELEMENT = re.compile(r"^[A-Za-z]{1,3}[0-9]*$")
_ROLES = {
    "x": "x", "y": "y", "z": "z",
    "nx": "nx", "ny": "ny", "nz": "nz",
    "normal_x": "nx", "normal_y": "ny", "normal_z": "nz",
    "r": "r", "g": "g", "b": "b", "a": "a",
    "red": "r", "green": "g", "blue": "b", "alpha": "a",
}  # fmt: skip
_SKIP = ("_", "skip", "ignore")
_CHEMISTRY = (
    "XYZ: this looks like a chemistry/molecular XYZ file (atom count, a comment "
    "line, then 'element x y z' rows), not a point cloud; meshio++ has no reader "
    "for it"
)


def _suffix(filename):
    name = filename if isinstance(filename, str) else getattr(filename, "name", "")
    name = str(name).lower()
    return name[name.rfind(".") :] if "." in name else ""


def _is_comment(line):
    return line.startswith("#") or line.startswith("//")


def _split(line, delimiter):
    if delimiter is None:
        return line.split()
    tokens = [t.strip() for t in line.split(delimiter)]
    if tokens and tokens[-1] == "":
        tokens.pop()
    return tokens


def _header_names(comments):
    for line in reversed(comments):
        body = line[2:] if line.startswith("//") else line.lstrip("#")
        tokens = [t for t in re.split(r"[\s,;]+", body.strip()) if t]
        if (
            len(tokens) >= 3
            and all(_NAME.match(t) for t in tokens)
            and [t.lower() for t in tokens[:3]] == ["x", "y", "z"]
        ):
            return tokens
    return None


def _colours(values):
    """(n, k) float columns -> uint8: integral bytes as they are, or unit floats * 255."""
    if np.all(values == np.rint(values)) and values.min() >= 0 and values.max() <= 255:
        return values.astype(np.uint8)
    if values.min() >= 0 and values.max() <= 1:
        return np.rint(values * 255).astype(np.uint8)
    raise ReadError("XYZ: colour columns must be bytes (0..255) or unit floats (0..1)")


def _looks_like_normals(block):
    norms = np.linalg.norm(block, axis=1)
    unit = np.abs(norms - 1.0) < 1e-2
    return bool(np.all(unit | (norms == 0)) and unit.any())


def _looks_like_colours(block):
    integral = (
        np.all(block == np.rint(block)) and block.min() >= 0 and block.max() <= 255
    )
    return bool(integral or (block.min() >= 0 and block.max() <= 1))


def _default_columns(table, suffix):
    ncols = table.shape[1]
    if ncols == 3:
        return ["x", "y", "z"]
    if ncols == 4:
        return ["x", "y", "z", "intensity" if suffix == ".pts" else "scalar"]
    if ncols == 6:
        if suffix == ".xyzn":
            return ["x", "y", "z", "nx", "ny", "nz"]
        if suffix == ".xyzrgb":
            return ["x", "y", "z", "r", "g", "b"]
        if _looks_like_normals(table[:, 3:6]):
            return ["x", "y", "z", "nx", "ny", "nz"]
        if _looks_like_colours(table[:, 3:6]):
            return ["x", "y", "z", "r", "g", "b"]
    if ncols == 7 and suffix == ".pts":
        return ["x", "y", "z", "intensity", "r", "g", "b"]
    raise ReadError(
        f"XYZ: cannot tell what the {ncols} columns are; pass columns=[...] "
        "(names x y z nx ny nz r g b a, any other name is a scalar, '_' skips a column)"
    )


def read(filename, columns=None, delimiter=None):
    suffix = _suffix(filename)
    with open_file(filename, "r") as f:
        text = f.read()
    if isinstance(text, (bytes, bytearray)):
        text = bytes(text).decode("latin-1")
    lines = [line.strip() for line in text.splitlines()]

    if len(lines) >= 3 and re.fullmatch(r"\d+", lines[0]):
        atom = lines[2].split()
        if (
            len(atom) >= 4
            and _ELEMENT.match(atom[0])
            and atom[0].lower() not in ("nan", "inf")
        ):
            raise ReadError(_CHEMISTRY)

    comments, rows, numbers = [], [], []
    for number, line in enumerate(lines, 1):
        if not line:
            continue
        if _is_comment(line):
            if not rows:
                comments.append(line)
            continue
        rows.append(line)
        numbers.append(number)

    declared = None
    if suffix == ".pts" and rows and re.fullmatch(r"\d+", rows[0]):
        declared = int(rows[0])
        rows, numbers = rows[1:], numbers[1:]

    if delimiter is None and rows:
        delimiter = ";" if ";" in rows[0] else "," if "," in rows[0] else None
    if delimiter is not None and delimiter.strip() == "":
        delimiter = None

    split = [_split(r, delimiter) for r in rows]
    if not split:
        empty = np.empty((0, 3))
        return Mesh(empty, [("vertex", np.empty((0, 1), dtype=np.int64))])
    ncols = len(split[0])
    for tokens, number in zip(split, numbers):
        if len(tokens) != ncols:
            raise ReadError(
                f"XYZ: line {number}: expected {ncols} columns, found {len(tokens)}"
            )
    try:
        table = np.array(split, dtype=np.float64)
    except ValueError:
        for tokens, number in zip(split, numbers):
            try:
                [float(t) for t in tokens]
            except ValueError:
                raise ReadError(
                    f"XYZ: line {number}: '{' '.join(tokens)}' is not numeric"
                )
        raise
    if declared is not None and declared != len(table):
        raise ReadError(
            f"XYZ: the .pts header declares {declared} points, found {len(table)}"
        )

    if columns is not None:
        names = list(columns)
        if len(names) != ncols:
            raise ReadError(
                f"XYZ: columns= names {len(names)} columns, the file has {ncols}"
            )
    else:
        names = _header_names(comments)
        if names is not None and len(names) != ncols:
            warn(
                f"XYZ: header names {len(names)} columns, the file has {ncols}; ignoring it"
            )
            names = None
        if names is None:
            names = _default_columns(table, suffix)
    return _to_mesh(table, names)


def _to_mesh(table, names):
    roles = {}
    scalars = []
    for j, name in enumerate(names):
        if name.lower() in _SKIP:
            continue
        role = _ROLES.get(name.lower())
        if role is None:
            scalars.append((name, j))
        elif role in roles:
            raise ReadError(f"XYZ: column '{name}' appears twice")
        else:
            roles[role] = j
    if not all(a in roles for a in "xyz"):
        raise ReadError("XYZ: the columns must include x, y and z")
    point_data = {}
    if any(a in roles for a in ("nx", "ny", "nz")):
        if not all(a in roles for a in ("nx", "ny", "nz")):
            raise ReadError("XYZ: normals need all of nx, ny, nz")
        point_data["normals"] = table[:, [roles["nx"], roles["ny"], roles["nz"]]].copy()
    if any(a in roles for a in ("r", "g", "b")):
        if not all(a in roles for a in "rgb"):
            raise ReadError("XYZ: colours need all of r, g, b")
        idx = [roles["r"], roles["g"], roles["b"]] + (
            [roles["a"]] if "a" in roles else []
        )
        point_data["rgba" if "a" in roles else "rgb"] = _colours(table[:, idx])
    elif "a" in roles:
        raise ReadError("XYZ: an alpha column needs r, g and b")
    for name, j in scalars:
        if name in point_data:
            raise ReadError(f"XYZ: column '{name}' appears twice")
        point_data[name] = table[:, j].copy()
    points = table[:, [roles["x"], roles["y"], roles["z"]]].copy()
    cells = [("vertex", np.arange(len(points), dtype=np.int64).reshape(-1, 1))]
    return Mesh(points, cells, point_data=point_data)


def _clean(name):
    return re.sub(r"[\s,;]+", "_", name) or "field"


def write(filename, mesh, float_fmt=None):
    points = np.asarray(mesh.points)
    if points.shape[1] < 3:
        warn("XYZ requires 3D points; padding with zeros.")
        _provenance.note("point-padding", "points padded with zero coordinates to 3D")
        pad = np.zeros((len(points), 3 - points.shape[1]), dtype=points.dtype)
        points = np.concatenate([points, pad], axis=1)
    skipped = [c.type for c in mesh.cells if c.type != "vertex"]
    if skipped:
        string = ", ".join(skipped)
        warn(f"XYZ holds points only. Skipping {string} cells.")
        _provenance.note(
            "cells-dropped", f"cell block(s) of type {string} have no XYZ equivalent"
        )
    if mesh.cell_data:
        _provenance.note("data-dropped", "cell data has no XYZ equivalent")

    n = len(points)
    names = ["x", "y", "z"]
    columns = [points[:, i] for i in range(3)]
    dropped = []
    for name in sorted(mesh.point_data):
        array = np.asarray(mesh.point_data[name])
        if array.dtype.kind not in "fiub":
            dropped.append(name)
            continue
        if array.dtype.kind == "b":
            array = array.astype(np.uint8)
        flat = array.reshape(n, -1) if n else array.reshape(0, 1)
        width = flat.shape[1]
        if name == "normals" and width == 3:
            labels = ["nx", "ny", "nz"]
        elif name == "rgb" and width == 3:
            labels = ["r", "g", "b"]
        elif name == "rgba" and width == 4:
            labels = ["r", "g", "b", "a"]
        elif width == 1:
            labels = [_clean(name)]
        else:
            labels = [f"{_clean(name)}_{k}" for k in range(width)]
        if any(label in names for label in labels):
            # e.g. rgb and rgba together, or point data called "x": the reader would
            # see a repeated role, so fall back to indexed names.
            labels = [f"{_clean(name)}_{k}" for k in range(width)]
        names += labels
        columns += [flat[:, k] for k in range(width)]
    if dropped:
        string = ", ".join(dropped)
        warn(f"XYZ cannot store point data {string}; skipping.")
        _provenance.note(
            "data-dropped", f"point data not representable in XYZ: {string}"
        )

    def render(value):
        if value.dtype.kind in "iu":
            return str(int(value))
        v = float(value)
        if v != v:
            return "nan"
        # 9 significant digits round-trip a float32; anything wider needs 17.
        return format(v, float_fmt or (".9g" if value.dtype.itemsize <= 4 else ".17g"))

    lines = [" ".join(render(col[i]) for col in columns) for i in range(n)]
    with open_file(filename, "wb") as fh:
        fh.write(_provenance.render_lines(_provenance.SlotTier.BLOCK, "# ").encode())
        fh.write(("# " + " ".join(names) + "\n").encode())
        if lines:
            fh.write(("\n".join(lines) + "\n").encode())
