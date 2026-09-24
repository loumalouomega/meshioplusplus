"""I/O for OpenRadioss / Radioss animation files (``<run>A001``, ``A002``...).

The pure-Python twin of ``src/cpp/src/formats/radioss_anim.cpp``: both engines
read the same meshes.

Each animation file is one state of the run, big-endian, magic ``0x542C``: its
time, three titles and ten flags, then a 2-D section (the nodes and the 4-node
facets of shells, with nodal scalars, vectors and 2-D tensors) and, as the
flags say, 3-D (8-node bricks), 1-D (2-node elements) and SPH sections, masses,
node and element ids, a part/subset/material hierarchy and time-history lists.
The layout follows OpenRadioss's ``anim_to_vtk`` converter (MIT).
"""

import numpy as np

from .._exceptions import ReadError
from .._mesh import Mesh
from .._regions import Region
from ..lsdyna._lsdyna import _collapse_solid

__all__ = ["read", "MAGIC"]

MAGIC = 0x542C


class _Cursor:
    def __init__(self, data, path):
        self.data = data
        self.pos = 0
        self.path = path

    def take(self, n):
        if self.pos + n > len(self.data):
            raise ReadError(
                f"Radioss animation: '{self.path}' is truncated (needs {n} more bytes "
                f"at offset {self.pos} of {len(self.data)})"
            )
        chunk = self.data[self.pos : self.pos + n]
        self.pos += n
        return chunk

    def ints(self, n):
        if n < 0:
            raise ReadError(f"Radioss animation: '{self.path}' has a negative count")
        return np.frombuffer(self.take(4 * n), dtype=">i4").astype(np.int64)

    def int(self):
        return int(self.ints(1)[0])

    def floats(self, n):
        return np.frombuffer(self.take(4 * n), dtype=">f4").astype(np.float64)

    def bytes(self, n):
        return np.frombuffer(self.take(n), dtype=np.uint8)

    def text(self, n):
        raw = self.take(n).split(b"\0", 1)[0]
        return raw.decode("latin-1").strip()

    def texts(self, count, n):
        return [self.text(n) for _ in range(count)]


def _parts(ends, texts, count):
    """Per element, the (id, name) of its part: part k ends before `ends[k]`."""
    ids = np.zeros(count, dtype=np.int64)
    names = []
    start = 0
    for end, text in zip(ends.tolist(), texts):
        head, _, tail = text.partition(":")
        try:
            pid = int(head)
        except ValueError:
            pid = 0
        ids[start : max(start, end)] = pid
        names.append(
            (pid, tail.strip() if ":" in text else text, start, max(start, end))
        )
        start = max(start, end)
    return ids, names


def _section(c, flags, n_elem, n_parts, n_scalars, tensor_width, n_tensors, extra):
    """The shared tail of the 3-D, 1-D and SPH sections, after their counts and
    connectivity: alive flags, parts, scalars, tensors, then `extra` (the 1-D
    skews), masses, ids and the part hierarchy."""
    s = {"alive": c.bytes(n_elem) != 0}
    s["part_ends"] = c.ints(n_parts)
    s["part_texts"] = c.texts(n_parts, 50)
    s["scalars"] = []
    if n_scalars:
        names = c.texts(n_scalars, 81)
        values = c.floats(n_scalars * n_elem).reshape(n_scalars, n_elem)
        s["scalars"] = list(zip(names, values))
    s["tensors"] = []
    if n_tensors:
        names = c.texts(n_tensors, 81)
        values = c.floats(n_elem * tensor_width * n_tensors)
        values = values.reshape(n_tensors, n_elem, tensor_width)
        s["tensors"] = list(zip(names, values))
    if extra:
        c.ints(n_elem)
    s["mass"] = c.floats(n_elem) if flags[0] == 1 else None
    s["ids"] = c.ints(n_elem) if flags[1] == 1 else None
    s["hierarchy"] = (
        [c.ints(n_parts) for _ in range(3)] if flags[4] else None
    )  # subset, material, property of each part
    return s


def _parse(data, path):
    c = _Cursor(data, path)
    magic = c.int()
    if magic != MAGIC:
        raise ReadError(
            f"Radioss animation: '{path}' has magic {magic:#x}; only the current "
            f"layout ({MAGIC:#x}) is read"
        )
    f = {"time": float(c.floats(1)[0])}
    c.texts(3, 81)  # "Time=", "ModAnim", "Radioss Run="
    flags = c.ints(10).tolist()
    nn, nf, np2, nfun, nefun, nvec, nten, nskew = c.ints(8).tolist()
    c.take(2 * 6 * nskew)
    f["points"] = c.floats(3 * nn).reshape(nn, 3)
    f["2d"] = {"conn": c.ints(4 * nf).reshape(nf, 4), "alive": c.bytes(nf) != 0}
    two = f["2d"]
    two["part_ends"] = c.ints(np2) if np2 else np.zeros(0, dtype=np.int64)
    two["part_texts"] = c.texts(np2, 50) if np2 else []
    c.take(2 * 3 * nn)  # nodal normals
    f["point_scalars"] = []
    two["scalars"] = []
    if nfun + nefun:
        names = c.texts(nfun + nefun, 81)
        values = c.floats(nn * nfun).reshape(nfun, nn)
        f["point_scalars"] = list(zip(names[:nfun], values))
        evalues = c.floats(nf * nefun).reshape(nefun, nf)
        two["scalars"] = list(zip(names[nfun:], evalues))
    vnames = c.texts(nvec, 81)
    vvalues = c.floats(3 * nn * nvec).reshape(nvec, nn, 3)
    f["vectors"] = list(zip(vnames, vvalues))
    two["tensors"] = []
    if nten:
        names = c.texts(nten, 81)
        values = c.floats(nf * 3 * nten).reshape(nten, nf, 3)
        two["tensors"] = list(zip(names, values))
    two["mass"] = f["point_mass"] = None
    if flags[0] == 1:
        two["mass"] = c.floats(nf)
        f["point_mass"] = c.floats(nn)
    f["node_ids"] = two["ids"] = None
    if flags[1]:
        f["node_ids"] = c.ints(nn)
        two["ids"] = c.ints(nf)
    two["hierarchy"] = [c.ints(np2) for _ in range(3)] if flags[4] else None

    f["3d"] = f["1d"] = f["sph"] = None
    if flags[2]:
        n3, p3, ef3, t3 = c.ints(4).tolist()
        conn = c.ints(8 * n3).reshape(n3, 8)
        f["3d"] = _section(c, flags, n3, p3, ef3, 6, t3, False)
        f["3d"]["conn"] = conn
    if flags[3]:
        n1, p1, ef1, to1, sk1 = c.ints(5).tolist()
        conn = c.ints(2 * n1).reshape(n1, 2)
        f["1d"] = _section(c, flags, n1, p1, ef1, 9, to1, sk1 != 0)
        f["1d"]["conn"] = conn
    if flags[4]:
        for _ in range(c.int()):  # subsets
            c.text(50)
            c.int()  # parent
            c.ints(c.int())  # child subsets
            for _ in range(3):  # 2-D, 3-D and 1-D parts
                c.ints(c.int())
        nmat, nprop = c.int(), c.int()
        c.texts(nmat, 50)
        c.ints(nmat)
        c.texts(nprop, 50)
        c.ints(nprop)
    if flags[5]:
        counts = c.ints(4).tolist()  # nodes, 2-D, 3-D, 1-D elements
        for count in counts:
            c.ints(count)
            c.texts(count, 50)
    if flags[7]:
        ns, ps, efs, ts = c.ints(4).tolist()
        conn = c.ints(ns) if ns else np.zeros(0, dtype=np.int64)
        alive = c.bytes(ns) != 0 if ns else np.zeros(0, dtype=bool)
        sph = {"conn": conn, "alive": alive}
        sph["part_ends"] = c.ints(ps) if ps else np.zeros(0, dtype=np.int64)
        sph["part_texts"] = c.texts(ps, 50) if ps else []
        sph["scalars"] = []
        if efs:
            names = c.texts(efs, 81)
            sph["scalars"] = list(zip(names, c.floats(efs * ns).reshape(efs, ns)))
        sph["tensors"] = []
        if ts:
            names = c.texts(ts, 81)
            sph["tensors"] = list(zip(names, c.floats(ns * ts * 6).reshape(ts, ns, 6)))
        sph["mass"] = c.floats(ns) if flags[0] == 1 else None
        sph["ids"] = c.ints(ns) if flags[1] == 1 else None
        sph["hierarchy"] = [c.ints(ps) for _ in range(3)] if flags[4] else None
        f["sph"] = sph
    return f


def _tensor6(values):
    """A 2-D tensor (xx, yy, xy) as the 3-D (xx yy zz xy yz zx): the rest NaN."""
    if values.shape[1] == 6:
        return values
    out = np.full((len(values), 6), np.nan)
    out[:, 0], out[:, 1], out[:, 3] = values[:, 0], values[:, 1], values[:, 2]
    return out


def read(filename):
    path = str(filename)
    with open(filename, "rb") as fh:
        data = fh.read()
    f = _parse(data, path)
    points = f["points"]
    nn = len(points)

    # The element families in file order: 1-D, 2-D, 3-D, SPH.
    families = []
    for key in ("1d", "2d", "3d", "sph"):
        s = f[key]
        if s is None or len(s["conn"]) == 0:
            continue
        conn = s["conn"]
        n = len(conn)
        if np.any((conn < 0) | (conn >= nn)):
            raise ReadError(f"Radioss animation: '{path}' names a node out of range")
        types, rows = [], []
        for r in conn.tolist():
            if key == "1d":
                t, row = "line", r
            elif key == "sph":
                t, row = "vertex", [r]
            elif key == "2d":
                row = list(dict.fromkeys(r))
                t = "quad" if len(row) == 4 else "triangle" if len(row) == 3 else None
            else:
                t, row = _collapse_solid(r)
            types.append(t)
            rows.append(row)
        pids, part_names = _parts(s["part_ends"], s["part_texts"], n)
        families.append((key, s, types, rows, pids, part_names))

    block_types = []
    members = {}  # type -> [(family index, element index)]
    for fi, (_, _, types, _, _, _) in enumerate(families):
        for e, t in enumerate(types):
            if t is None:
                continue
            if t not in members:
                members[t] = []
                block_types.append(t)
            members[t].append((fi, e))
    skipped = sum(t is None for fam in families for t in fam[2])
    if skipped:
        from .._common import warn

        warn(f"Radioss animation: {skipped} degenerate facet(s) skipped")

    cells = []
    cell_of = {}
    for t in block_types:
        rows = [families[fi][3][e] for fi, e in members[t]]
        for fi, e in members[t]:
            cell_of[(fi, e)] = len(cell_of)
        cells.append((t, np.array(rows, dtype=np.int64)))
    mesh = Mesh(points, cells)
    mesh.field_data["meshio:time"] = np.array([f["time"]])
    # The one step's time, for `read_metadata` (and so a sequence) to report.
    mesh.time_values = [f["time"]]

    for name, values in f["point_scalars"]:
        mesh.point_data[name] = values
    for name, values in f["vectors"]:
        mesh.point_data[name] = values
    if f["node_ids"] is not None:
        mesh.point_data["radioss:node_id"] = f["node_ids"]
    if f["point_mass"] is not None:
        mesh.point_data["radioss:mass"] = f["point_mass"]
    if not cells:
        return mesh

    def per_block(get, width, fill, dtype):
        """One array per block from `get(family_index, element_index)`."""
        out = []
        for t in block_types:
            shape = (len(members[t]), width) if width else (len(members[t]),)
            a = np.full(shape, fill, dtype=dtype)
            for r, (fi, e) in enumerate(members[t]):
                v = get(fi, e)
                if v is not None:
                    a[r] = v
            out.append(a)
        return out

    mesh.cell_data["radioss:part"] = per_block(
        lambda fi, e: families[fi][4][e], 0, 0, np.int64
    )
    # A deleted (eroded) element's flag is 0; any other byte is alive (the
    # facets write 0xFF, bricks and beams 1).
    mesh.cell_data["radioss:alive"] = per_block(
        lambda fi, e: families[fi][1]["alive"][e], 0, 0, np.int8
    )
    if all(fam[1]["ids"] is not None for fam in families):
        mesh.cell_data["radioss:element_id"] = per_block(
            lambda fi, e: families[fi][1]["ids"][e], 0, 0, np.int64
        )
    if all(fam[1]["mass"] is not None for fam in families):
        mesh.cell_data["radioss:mass"] = per_block(
            lambda fi, e: families[fi][1]["mass"][e], 0, 0.0, np.float64
        )
    if all(fam[1]["hierarchy"] is not None for fam in families):
        for k, name in ((1, "radioss:material"), (2, "radioss:property")):
            lookups = []
            for fam in families:
                ends = fam[1]["part_ends"].tolist()
                values = fam[1]["hierarchy"][k].tolist()
                per = np.zeros(len(fam[2]), dtype=np.int64)
                start = 0
                for end, v in zip(ends, values):
                    per[start : max(start, end)] = v
                    start = max(start, end)
                lookups.append(per)
            mesh.cell_data[name] = per_block(
                lambda fi, e: lookups[fi][e], 0, 0, np.int64
            )

    # Element results: the union of the families' names, NaN where a family has
    # none; tensors as xx yy zz xy yz zx (a 2-D one's out-of-plane parts NaN),
    # 1-D force/moment sets with their nine components.
    scalars, tensors = {}, {}
    for fi, fam in enumerate(families):
        for name, values in fam[1]["scalars"]:
            scalars.setdefault(name, {})[fi] = values
        for name, values in fam[1]["tensors"]:
            tensors.setdefault(name, {})[fi] = values
    for name, by_family in scalars.items():
        mesh.cell_data[name] = per_block(
            lambda fi, e: by_family[fi][e] if fi in by_family else None,
            0,
            np.nan,
            np.float64,
        )
    for name, by_family in tensors.items():
        width = 9 if any(v.shape[1] == 9 for v in by_family.values()) else 6
        converted = {
            fi: (v if width == 9 else _tensor6(v)) for fi, v in by_family.items()
        }
        if width == 9 and any(v.shape[1] != 9 for v in by_family.values()):
            raise ReadError(
                f"Radioss animation: '{path}' names both a 1-D force set and a tensor "
                f"'{name}'"
            )
        mesh.cell_data[name] = per_block(
            lambda fi, e: converted[fi][e] if fi in converted else None,
            width,
            np.nan,
            np.float64,
        )

    # Parts -> cell regions (a part id can name a 1-D and a 2-D part).
    regions = {}
    for fi, fam in enumerate(families):
        for pid, name, start, end in fam[5]:
            key = (name or f"Part {pid}", pid)
            cells_in = [
                cell_of[(fi, e)] for e in range(start, end) if (fi, e) in cell_of
            ]
            regions.setdefault(key, []).extend(cells_in)
    dims = {"vertex": 0, "line": 1, "triangle": 2, "quad": 2}
    cell_dim = []
    for t in block_types:
        cell_dim += [dims.get(t, 3)] * len(members[t])
    mesh.regions = [
        Region(
            name,
            "cell",
            np.array(sorted(ids), dtype=np.int64),
            max((cell_dim[c] for c in ids), default=-1),
            pid,
        )
        for (name, pid), ids in regions.items()
    ]
    return mesh
