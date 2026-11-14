"""ADCIRC ``fort.14`` grids with ``fort.63``/``fort.64`` results: the Python reference
of ``formats/adcirc.cpp`` (roadmap 2.4.3.3). See that header for the layout."""

from __future__ import annotations

import os

import numpy as np

from .._exceptions import ReadError, WriteError
from .._mesh import Mesh
from .._regions import Region

BOUNDARY_PREFIX = "adcirc:boundary:"
_NAMES = {
    "fort.63": "zeta",
    "fort.64": "velocity",
    "fort.73": "pressure",
    "fort.74": "wind_velocity",
}


class _Lines:
    def __init__(self, path):
        self.path = str(path)
        try:
            with open(self.path, "rb") as f:
                text = f.read().decode("latin-1")
        except OSError:
            raise ReadError(f"meshio++: adcirc: could not open file: {path}") from None
        self.lines = text.split("\n")
        if self.lines and self.lines[-1] == "":
            self.lines.pop()
        self.lines = [ln[:-1] if ln.endswith("\r") else ln for ln in self.lines]
        self.at = 0

    def eof(self):
        return self.at >= len(self.lines)

    def next(self, what):
        if self.at >= len(self.lines):
            raise ReadError(f"meshio++: adcirc: '{self.path}' ends before its {what}")
        self.at += 1
        return self.lines[self.at - 1]


def _numbers(line):
    out = []
    for tok in line.replace(",", " ").split():
        if not (tok[0].isdigit() or tok[0] in "-+."):
            break
        out.append(tok)
    return out


def _need(lines, what, count):
    t = _numbers(lines.next(what))
    if len(t) < count:
        raise ReadError(f"meshio++: adcirc: {what} needs {count} numbers")
    return t


def _int(tok):
    try:
        return int(tok)
    except ValueError:
        try:
            v = float(tok)
        except ValueError:
            raise ReadError(f"meshio++: adcirc: '{tok}' is not an integer") from None
        return int(v)


def _real(tok):
    try:
        return float(tok.replace("d", "e").replace("D", "e"))
    except ValueError:
        raise ReadError(f"meshio++: adcirc: '{tok}' is not a number") from None


def _boundary_key(region):
    """Boundaries in their number's order (``open_boundary_10`` after ``_2``)."""
    tail = region.name.rsplit("_", 1)[-1]
    return (int(tail) if tail.isdigit() else 1 << 62, region.name)


def _is_barrier(kind):
    return kind % 10 in (3, 4, 5)


def _read_grid(path):
    lines = _Lines(path)
    title = lines.next("title")
    counts = _need(lines, "element and node counts (NE NP)", 2)
    limit = len(lines.lines)
    ne, np_ = _int(counts[0]), _int(counts[1])
    for what, n in (("element", ne), ("node", np_)):
        if n < 0 or n > limit:
            raise ReadError(
                f"meshio++: adcirc: {what} count {n} exceeds what the file holds ({limit})"
            )
    index = {}
    points = np.zeros((np_, 3))
    depth = np.zeros(np_)
    for i in range(np_):
        t = _need(lines, "node", 4)
        key = _int(t[0])
        if key in index:
            raise ReadError(f"meshio++: adcirc: duplicate node id {t[0]}")
        index[key] = i
        points[i, 0], points[i, 1] = _real(t[1]), _real(t[2])
        depth[i] = _real(t[3])

    def node(tok):
        try:
            return index[_int(tok)]
        except KeyError:
            raise ReadError(
                f"meshio++: adcirc: node {tok} is not in the grid"
            ) from None

    tri = np.zeros((ne, 3), dtype=np.int64)
    for e in range(ne):
        t = _need(lines, "element", 5)
        if _int(t[1]) != 3:
            raise ReadError(
                f"meshio++: adcirc: element {t[0]} has {t[1]} nodes; ADCIRC grids are "
                "triangles"
            )
        tri[e] = [node(x) for x in t[2:5]]
    regions, field_data = [], {
        "adcirc:title": np.frombuffer(title.encode("latin-1"), dtype=np.uint8).copy()
    }

    def skip_blank():
        while not lines.eof() and not _numbers(lines.lines[lines.at]):
            lines.at += 1

    def strings(kind, is_open):
        count = _int(_need(lines, "NOPE" if is_open else "NBOU", 1)[0])
        if count < 0 or count > limit:
            raise ReadError(f"meshio++: adcirc: boundary count {count} is out of range")
        _need(lines, "NETA" if is_open else "NVEL", 1)
        for k in range(count):
            h = _need(lines, "boundary header", 1)
            n = _int(h[0])
            if n < 0 or n > limit:
                raise ReadError(f"meshio++: adcirc: node count {n} is out of range")
            kind_tag = _int(h[1]) if len(h) > 1 else 0
            ids = np.asarray(
                [node(_need(lines, "boundary node", 1)[0]) for _ in range(n)],
                dtype=np.int64,
            )
            name = f"{kind}_{k}"
            field_data[BOUNDARY_PREFIX + name] = ids
            regions.append(Region(name, "point", ids, tag=kind_tag))

    skip_blank()
    if not lines.eof():
        strings("open_boundary", True)
        skip_blank()
        if not lines.eof():
            strings("land_boundary", False)
    mesh = Mesh(
        points,
        [("triangle", tri)] if ne else [],
        point_data={"depth": depth},
        field_data=field_data,
        regions=regions,
    )
    return mesh, index


def _scan_result(lines, num_points):
    base = os.path.basename(lines.path)
    name = _NAMES.get(base, f"adcirc:{base}")
    lines.next("title")
    head = _need(lines, "dataset header (NDSETS NP DT NSPOOL IRTYPE)", 5)
    if _int(head[1]) != num_points:
        raise ReadError(
            f"meshio++: adcirc: '{lines.path}' has {head[1]} nodes but the grid has "
            f"{num_points}"
        )
    width = _int(head[4])
    if width not in (1, 2, 3):
        raise ReadError(f"meshio++: adcirc: '{lines.path}' has record type {head[4]}")
    times, starts = [], []
    while not lines.eof():
        t = _numbers(lines.next("dataset"))
        if not t:
            continue
        if len(t) < 2:
            raise ReadError("meshio++: adcirc: a dataset header needs TIME IT")
        starts.append(lines.at - 1)
        times.append(_real(t[0]))
        rows = _int(t[2]) if len(t) >= 4 else num_points
        if rows > num_points:
            raise ReadError(
                "meshio++: adcirc: a sparse dataset lists more nodes than the grid"
            )
        lines.at += rows
    return name, width, times, starts


def _read_dataset(lines, width, start, index, num_points):
    lines.at = start
    t = _numbers(lines.next("dataset header"))
    sparse = len(t) >= 4
    rows = _int(t[2]) if sparse else num_points
    fill = _real(t[3]) if sparse else 0.0
    out = np.full((num_points, width), fill)
    for _ in range(rows):
        v = _need(lines, "dataset value", 1 + width)
        try:
            i = index[_int(v[0])]
        except KeyError:
            raise ReadError(
                f"meshio++: adcirc: a dataset names node {v[0]}, which the grid does "
                "not have"
            ) from None
        out[i] = [_real(x) for x in v[1 : 1 + width]]
    return out[:, 0] if width == 1 else out


def _result_paths(path, results, siblings):
    if results:
        return [str(r) for r in results]
    if not siblings or os.path.basename(str(path)) != "fort.14":
        return []
    folder = os.path.dirname(str(path))
    return [
        os.path.join(folder, n)
        for n in ("fort.63", "fort.64")
        if os.path.isfile(os.path.join(folder, n))
    ]


def _resolve(times, time_step):
    n = len(times)
    if n == 0:
        if time_step in (0, -1):
            return 0
        raise ReadError(
            f"meshio++: time step {time_step} requested, but this file carries no time steps"
        )
    k = time_step + n if time_step < 0 else time_step
    if not 0 <= k < n:
        raise ReadError(
            f"meshio++: time step {time_step} is out of range: this file has {n} "
            f"step{'' if n == 1 else 's'}"
        )
    return k


def read(filename, results=None, siblings=True, time_step=0):
    mesh, index = _read_grid(filename)
    num_points = len(mesh.points)
    paths = _result_paths(filename, results, siblings)
    times, step = None, 0
    for path in paths:
        lines = _Lines(path)
        name, width, ftimes, starts = _scan_result(lines, num_points)
        if times is None:
            times = ftimes
            step = _resolve(times, time_step)
        elif len(ftimes) != len(times):
            raise ReadError(
                f"meshio++: adcirc: '{path}' has {len(ftimes)} datasets but "
                f"'{paths[0]}' has {len(times)}"
            )
        if not ftimes:
            continue
        mesh.point_data[name] = _read_dataset(
            lines, width, starts[step], index, num_points
        )
    if times:
        mesh.field_data["meshio:time"] = np.array([times[step]])
    return mesh


def time_values(filename):
    mesh, _ = _read_grid(filename)
    paths = _result_paths(filename, None, True)
    if not paths:
        return []
    return _scan_result(_Lines(paths[0]), len(mesh.points))[2]


def write(filename, mesh):
    for block in mesh.cells:
        if block.type != "triangle":
            raise WriteError(
                f"meshio++: adcirc: a '{block.type}' block; ADCIRC grids hold triangles "
                "only"
            )
    opened, land = [], []
    for r in sorted(mesh.regions, key=lambda r: (r.kind, r.name, r.dim, r.tag)):
        if r.kind != "point":
            continue
        if r.name.startswith("open_boundary"):
            opened.append(r)
        elif r.name.startswith("land_boundary"):
            if _is_barrier(r.tag):
                raise WriteError(
                    f"meshio++: adcirc: land boundary '{r.name}' has type {r.tag}, whose "
                    "records need barrier columns meshio++ does not keep"
                )
            land.append(r)
    title = "meshio++"
    if "adcirc:title" in mesh.field_data:
        title = np.asarray(mesh.field_data["adcirc:title"], dtype=np.uint8).tobytes()
        title = title.decode("latin-1")
    points = np.asarray(mesh.points, dtype=np.float64)
    depth = mesh.point_data.get("depth")
    depth = np.zeros(len(points)) if depth is None else np.asarray(depth, dtype=float)
    tris = (
        np.concatenate([np.asarray(b.data, dtype=np.int64) for b in mesh.cells])
        if mesh.cells
        else np.zeros((0, 3), dtype=np.int64)
    )
    out = [f"{title}\n{len(tris)} {len(points)}\n"]
    for i, p in enumerate(points):
        y = p[1] if len(p) > 1 else 0.0
        out.append("%d %.17g %.17g %.17g\n" % (i + 1, p[0], y, depth[i]))
    for e, t in enumerate(tris):
        out.append(f"{e + 1} 3 {t[0] + 1} {t[1] + 1} {t[2] + 1}\n")

    def strings(regions, is_open):
        total = sum(len(r.entries) for r in regions)
        kind = "open" if is_open else "land"
        out.append(f"{len(regions)} = Number of {kind} boundaries\n")
        out.append(f"{total} = Total number of {kind} boundary nodes\n")
        for k, r in enumerate(regions):
            head = f"{len(r.entries)}"
            if not is_open or r.tag >= 0:
                head += f" {max(int(r.tag), 0)}"
            out.append(f"{head} = Number of nodes for {kind} boundary {k + 1}\n")
            ids = np.asarray(r.entries, dtype=np.int64)
            order = mesh.field_data.get(BOUNDARY_PREFIX + r.name)
            if order is not None:
                order = np.asarray(order, dtype=np.int64).reshape(-1)
                if len(order) == len(ids) and np.array_equal(np.sort(order), ids):
                    ids = order
            out.extend(f"{int(i) + 1}\n" for i in ids)

    strings(sorted(opened, key=_boundary_key), True)
    strings(sorted(land, key=_boundary_key), False)
    with open(filename, "w", encoding="latin-1", newline="\n") as f:
        f.write("".join(out))
