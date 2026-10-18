"""Write OpenRadioss / Radioss animation files (``<run>A001``...).

The pure-Python twin of ``write_radioss_anim`` in
``src/cpp/src/formats/radioss_anim.cpp``: both engines write the same bytes.
One mesh is one state, big-endian, magic ``0x542C``, the layout the reader
(``_anim.py``) and OpenRadioss's ``anim_to_vtk`` read. A series is one file per
step, numbered from the number in the first file's name.

Element families, in file order: ``vertex`` (SPH particles), ``line`` (1-D),
``triangle`` and ``quad`` (2-D facets, a triangle with its last node repeated),
and ``tetra``, ``pyramid``, ``wedge`` and ``hexahedron`` (3-D bricks, the
degenerate ones with repeated nodes). Other cell types are skipped with a
warning.
"""

from __future__ import annotations

import os
import re

import numpy as np

from .._common import warn
from .._exceptions import WriteError
from .._files import is_buffer
from .._regions import canonical_regions
from ..lsdyna._lsdyna import _expand_solid
from ._anim import MAGIC, TH_PREFIX

__all__ = ["write", "SeriesWriter"]

_DIM_TYPES = {
    "vertex": "sph",
    "line": "1d",
    "triangle": "2d",
    "quad": "2d",
    "tetra": "3d",
    "pyramid": "3d",
    "wedge": "3d",
    "hexahedron": "3d",
}
_FAMILIES = ("2d", "3d", "1d", "sph")  # the order the sections follow the nodes
_SPECIAL = {
    "radioss:part",
    "radioss:alive",
    "radioss:element_id",
    "radioss:mass",
    "radioss:material",
    "radioss:property",
}


def _i4(values):
    return np.asarray(values, dtype=">i4").tobytes()


def _f4(values):
    return np.asarray(values, dtype=">f4").tobytes()


def _text(value, width):
    raw = str(value).encode("latin-1", "replace")[: width - 1]
    return raw + b"\0" * (width - len(raw))


def _alive_bytes(alive, family):
    on = 0xFF if family == "2d" else 1
    return np.where(alive, on, 0).astype(np.uint8).tobytes()


def _concat(mesh, name):
    """The block-major cell array `name`, or None."""
    if name not in mesh.cell_data:
        return None
    parts = [np.asarray(a) for a in mesh.cell_data[name]]
    return np.concatenate(parts) if parts else np.zeros(0)


class _Family:
    def __init__(self, key):
        self.key = key
        self.index = []  # global cell numbers, ordered by part
        self.rows = []  # file node rows
        self.parts = []  # part id per element

    def __len__(self):
        return len(self.index)


def _classify(mesh):
    """The element families (global cell numbers, file rows) and the number of
    cells that have no place in the file."""
    families = {key: _Family(key) for key in _FAMILIES}
    skipped = {}
    base = 0
    for block in mesh.cells:
        data = np.asarray(block.data, dtype=np.int64)
        key = _DIM_TYPES.get(block.type)
        if key is None:
            skipped[block.type] = skipped.get(block.type, 0) + len(data)
        else:
            fam = families[key]
            for r, row in enumerate(data.tolist()):
                fam.index.append(base + r)
                if key == "sph":
                    fam.rows.append([row[0]])
                elif key == "2d":
                    fam.rows.append(row if len(row) == 4 else row + [row[-1]])
                elif key == "3d":
                    fam.rows.append(_expand_solid(block.type, row))
                else:
                    fam.rows.append(row[:2])
        base += len(data)
    return families, skipped, base


def _region_parts(mesh, ncells):
    """(part id, name) per cell from the cell regions, or None."""
    regions = [
        r
        for r in canonical_regions(mesh)
        if r.kind == "cell" and not r.name.startswith(TH_PREFIX)
    ]
    if not regions:
        return None
    ids = np.full(ncells, -1, dtype=np.int64)
    names = {}
    for k, region in enumerate(regions):
        pid = region.tag if region.tag >= 0 else k + 1
        names.setdefault(pid, region.name)
        for c in region.entries.tolist():
            if 0 <= c < ncells and ids[c] < 0:
                ids[c] = pid
    return ids, names


def _bytes_state(mesh, time):
    points = np.asarray(mesh.points, dtype=np.float64)
    nn = len(points)
    if points.ndim != 2 or points.shape[1] not in (2, 3):
        points = points.reshape(nn, -1)[:, :3]
    if points.shape[1] == 2:
        points = np.hstack([points, np.zeros((nn, 1))])
    families, skipped, ncells = _classify(mesh)
    for ctype, count in skipped.items():
        warn(f"Radioss animation: {count} '{ctype}' cell(s) are not written")

    # Per-cell attributes, in block-major order.
    part_arr = _concat(mesh, "radioss:part")
    region_parts = None
    if part_arr is None:
        region_parts = _region_parts(mesh, ncells)
    names_by_pid = region_parts[1] if region_parts else {}
    if part_arr is not None:
        for r in canonical_regions(mesh):
            if r.kind == "cell" and r.tag >= 0 and not r.name.startswith(TH_PREFIX):
                names_by_pid.setdefault(r.tag, r.name)
    alive = _concat(mesh, "radioss:alive")
    eids = _concat(mesh, "radioss:element_id")
    emass = _concat(mesh, "radioss:mass")
    emat = _concat(mesh, "radioss:material")
    eprop = _concat(mesh, "radioss:property")
    node_ids = mesh.point_data.get("radioss:node_id")
    nmass = mesh.point_data.get("radioss:mass")

    def pid_of(c):
        if part_arr is not None:
            return int(part_arr[c])
        if region_parts is not None and region_parts[0][c] >= 0:
            return int(region_parts[0][c])
        return 1

    for fam in families.values():
        order = sorted(range(len(fam.index)), key=lambda k: pid_of(fam.index[k]))
        fam.index = [fam.index[k] for k in order]
        fam.rows = [fam.rows[k] for k in order]
        fam.parts = [pid_of(c) for c in fam.index]

    # Time-history items: each entry of a `radioss_th:<label>` region, as the index
    # of a node or of an element in its family (the order just sorted into).
    position = {}
    for key, fam in families.items():
        for k, c in enumerate(fam.index):
            position[c] = (key, k)
    th_items = [[], [], [], []]
    for region in canonical_regions(mesh):
        if not region.name.startswith(TH_PREFIX):
            continue
        label = region.name[len(TH_PREFIX) :]
        for v in region.entries.tolist():
            if region.kind == "point":
                if 0 <= v < nn:
                    th_items[0].append((v, label))
            elif region.kind == "cell" and v in position:
                key, k = position[v]
                if key != "sph":
                    th_items[{"2d": 1, "3d": 2, "1d": 3}[key]].append((k, label))

    # Nodal and element results.
    scalars, vectors = [], []
    for name, values in sorted(mesh.point_data.items()):
        if name in ("radioss:node_id", "radioss:mass"):
            continue
        a = np.asarray(values)
        if a.ndim == 2 and a.shape[1] == 1:
            a = a[:, 0]
        if a.ndim == 1 and a.dtype.kind in "fiu":
            scalars.append((name, a.astype(np.float64)))
        elif a.ndim == 2 and a.shape[1] == 3 and a.dtype.kind in "fiu":
            vectors.append((name, a.astype(np.float64)))
        else:
            warn(f"Radioss animation: point data '{name}' is not written")
    escalars, tensors, forces = [], [], []
    for name in sorted(mesh.cell_data):
        if name in _SPECIAL:
            continue
        a = _concat(mesh, name)
        if a.ndim == 2 and a.shape[1] == 1:
            a = a[:, 0]
        if a.dtype.kind not in "fiu":
            warn(f"Radioss animation: cell data '{name}' is not written")
        elif a.ndim == 1:
            escalars.append((name, a.astype(np.float64)))
        elif a.ndim == 2 and a.shape[1] == 6:
            tensors.append((name, a.astype(np.float64)))
        elif a.ndim == 2 and a.shape[1] == 9:
            forces.append((name, a.astype(np.float64)))
        else:
            warn(f"Radioss animation: cell data '{name}' is not written")

    have_mass = emass is not None or nmass is not None
    have_ids = eids is not None or node_ids is not None
    have_hier = emat is not None or eprop is not None
    flags = [
        1 if have_mass else 0,
        1 if have_ids else 0,
        1 if len(families["3d"]) else 0,
        1 if len(families["1d"]) else 0,
        1 if have_hier else 0,
        1 if any(th_items) else 0,
        0,
        1 if len(families["sph"]) else 0,
        0,
        0,
    ]

    def part_table(fam):
        """Part ends, titles and (subset, material, property) per part."""
        ends, texts, mats, props = [], [], [], []
        for k, pid in enumerate(fam.parts):
            last = k + 1 == len(fam.parts) or fam.parts[k + 1] != pid
            if not last:
                continue
            ends.append(k + 1)
            name = names_by_pid.get(pid, f"Part {pid}")
            texts.append(f"{pid:9d}:{name}")
            c = fam.index[k]
            mats.append(int(emat[c]) if emat is not None else 0)
            props.append(int(eprop[c]) if eprop is not None else 0)
        return ends, texts, mats, props

    out = bytearray()
    out += _i4([MAGIC]) + _f4([time])
    for t in ("Time=", "ModAnim", "Radioss Run="):
        out += _text(t, 81)
    out += _i4(flags)

    # 2-D: the nodes and the facets.
    f2 = families["2d"]
    ends2, texts2, mats2, props2 = part_table(f2)
    tens2 = tensors
    out += _i4(
        [
            nn,
            len(f2),
            len(ends2),
            len(scalars),
            len(escalars),
            len(vectors),
            len(tens2),
            0,
        ]
    )
    out += _f4(points.reshape(-1))
    out += _i4(np.asarray(f2.rows, dtype=np.int64).reshape(-1))
    out += _alive_bytes(
        np.array(
            [bool(alive[c]) if alive is not None else True for c in f2.index],
            dtype=bool,
        ),
        "2d",
    )
    out += _i4(ends2)
    for t in texts2:
        out += _text(t, 50)
    out += np.zeros((nn, 3), dtype=">u2").tobytes()  # nodal normals
    if scalars or escalars:
        for name, _ in scalars + escalars:
            out += _text(name, 81)
        for _, a in scalars:
            out += _f4(a)
        for _, a in escalars:
            out += _f4(a[f2.index])
    for name, _ in vectors:
        out += _text(name, 81)
    for _, a in vectors:
        out += _f4(a.reshape(-1))
    if tens2:
        for name, _ in tens2:
            out += _text(name, 81)
        for _, a in tens2:
            out += _f4(a[f2.index][:, [0, 1, 3]].reshape(-1))
    if flags[0]:
        out += _f4([emass[c] if emass is not None else 0.0 for c in f2.index])
        out += _f4(nmass if nmass is not None else np.zeros(nn))
    if flags[1]:
        out += _i4(node_ids if node_ids is not None else np.arange(1, nn + 1))
        out += _i4([eids[c] if eids is not None else c + 1 for c in f2.index])
    if flags[4]:
        out += _i4([1] * len(ends2)) + _i4(mats2) + _i4(props2)

    hierarchy = {}

    def section(fam, with_counts, scal, tens):
        ends, texts, mats, props = part_table(fam)
        hierarchy[fam.key] = len(ends)
        sec = bytearray()
        sec += _i4(with_counts(len(ends)))
        sec += _i4(np.asarray(fam.rows, dtype=np.int64).reshape(-1))
        sec += _alive_bytes(
            np.array(
                [bool(alive[c]) if alive is not None else True for c in fam.index],
                dtype=bool,
            ),
            fam.key,
        )
        sec += _i4(ends)
        for t in texts:
            sec += _text(t, 50)
        if scal:
            for name, _ in scal:
                sec += _text(name, 81)
            for _, a in scal:
                sec += _f4(a[fam.index])
        if tens:
            for name, _ in tens:
                sec += _text(name, 81)
            for _, a in tens:
                sec += _f4(a[fam.index].reshape(-1))
        if flags[0]:
            sec += _f4([emass[c] if emass is not None else 0.0 for c in fam.index])
        if flags[1]:
            sec += _i4([eids[c] if eids is not None else c + 1 for c in fam.index])
        if flags[4]:
            sec += _i4([1] * len(ends)) + _i4(mats) + _i4(props)
        return sec

    f3 = families["3d"]
    if flags[2]:
        # The counts precede the connectivity: n, parts, scalars, tensors.
        out += section(
            f3,
            lambda p: [len(f3), p, len(escalars), len(tensors)],
            escalars,
            tensors,
        )
    f1 = families["1d"]
    if flags[3]:
        out += section(
            f1,
            lambda p: [len(f1), p, len(escalars), len(forces), 0],
            escalars,
            forces,
        )
    if flags[4]:
        out += _hierarchy(families, hierarchy, part_table)
    if flags[5]:
        out += _i4([len(items) for items in th_items])
        for items in th_items:
            out += _i4([i for i, _ in items])
            for _, label in items:
                out += _text(label, 50)
    fs = families["sph"]
    if flags[7]:
        out += section(
            fs,
            lambda p: [len(fs), p, len(escalars), len(tensors)],
            escalars,
            tensors,
        )
    return bytes(out)


def _hierarchy(families, counts, part_table):
    """One subset holding every part, and the material and property tables."""
    out = bytearray()
    out += _i4([1]) + _text("model", 50) + _i4([0, 0])
    for key in ("2d", "3d", "1d"):
        n = counts.get(key, 0) if key != "2d" else len(part_table(families["2d"])[0])
        out += _i4([n]) + _i4(range(n))
    mats = sorted({m for k in families for m in part_table(families[k])[2]})
    props = sorted({p for k in families for p in part_table(families[k])[3]})
    out += _i4([len(mats), len(props)])
    for m in mats:
        out += _text(f"Material {m}", 50)
    out += _i4(mats)
    for p in props:
        out += _text(f"Property {p}", 50)
    out += _i4(props)
    return out


def write(filename, mesh):
    """Write one OpenRadioss animation file (the state of ``mesh``): the name
    is ``<run>A001``-shaped, with no extension."""
    if is_buffer(filename, "w"):
        raise WriteError("Radioss animation: cannot be written to a stream")
    time = 0.0
    fd = getattr(mesh, "field_data", {}) or {}
    if "meshio:time" in fd:
        time = float(np.asarray(fd["meshio:time"]).reshape(-1)[0])
    with open(filename, "wb") as fh:
        fh.write(_bytes_state(mesh, time))


def numbered_names(filename):
    """The file names ``<stem>A<n>``, ``<stem>A<n+1>``... from the first one.
    A name that is not animation-shaped is a stem: ``<name>A001``, ..."""
    path = os.fspath(filename)
    m = re.match(r"^(.*A)(\d{3,})$", os.path.basename(path))
    directory = os.path.dirname(path)
    if m is None:
        stem, width, first = os.path.basename(path) + "A", 3, 1
    else:
        stem, width, first = m.group(1), len(m.group(2)), int(m.group(2))

    def names(number):
        while True:
            yield os.path.join(directory, f"{stem}{number:0{width}d}")
            number += 1

    return names(first)


class SeriesWriter:
    """One animation file per step, numbered from the number in ``filename``
    (``runA001``, ``runA002``...)."""

    def __init__(self, filename):
        self._names = numbered_names(filename)
        self.paths = []

    def write(self, time, mesh):
        path = next(self._names)
        with open(path, "wb") as fh:
            fh.write(_bytes_state(mesh, float(time)))
        self.paths.append(path)

    def close(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False
