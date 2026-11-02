"""
I/O for the Wavefront .obj file format, cf.
<https://en.wikipedia.org/wiki/Wavefront_.obj_file>.
"""

import numpy as np

from .. import _provenance
from .._exceptions import WriteError
from .._files import open_file
from .._mesh import CellBlock, Mesh


def read(filename):
    with open_file(filename, "r") as f:
        mesh = read_buffer(f)
    return mesh


def read_buffer(f):
    points = []
    vertex_normals = []
    texture_coords = []
    face_groups = []
    face_vts = []
    face_group_ids = []
    face_group_id = -1
    while True:
        line = f.readline()

        if not line:
            # EOF
            break

        strip = line.strip()

        if len(strip) == 0 or strip[0] == "#":
            continue

        split = strip.split()

        if split[0] == "v":
            points.append([float(item) for item in split[1:]])
        elif split[0] == "vn":
            vertex_normals.append([float(item) for item in split[1:]])
        elif split[0] == "vt":
            texture_coords.append([float(item) for item in split[1:]])
        elif split[0] == "s":
            # "s 1" or "s off" controls smooth shading
            pass
        elif split[0] == "f":
            dat = [int(item.split("/")[0]) for item in split[1:]]
            vts = [_vt_index(item) for item in split[1:]]
            if len(face_groups) == 0 or (
                len(face_groups[-1]) > 0 and len(face_groups[-1][-1]) != len(dat)
            ):
                face_groups.append([])
                face_vts.append([])
                face_group_ids.append([])

            face_groups[-1].append(dat)
            face_vts[-1].append(vts)
            face_group_ids[-1].append(face_group_id)
        elif split[0] == "g":
            # new group
            face_groups.append([])
            face_vts.append([])
            face_group_ids.append([])
            face_group_id += 1
        else:
            # who knows
            pass

    # There may be empty groups, too. <https://github.com/nschloe/meshio/issues/770>
    # Remove them.
    face_vts = [v for f, v in zip(face_groups, face_vts) if len(f) > 0]
    face_groups = [f for f in face_groups if len(f) > 0]
    face_group_ids = [g for g in face_group_ids if len(g) > 0]

    points = np.array(points)
    texture_coords = np.array(texture_coords)
    vertex_normals = np.array(vertex_normals)
    num_in = len(points)

    # convert to numpy arrays
    face_groups = [np.array(f) for f in face_groups]
    uv, parent = _corner_uv(points, texture_coords, face_groups, face_vts)
    if uv is not None and len(parent) > num_in:
        points = points[parent]

    point_data = {}

    # The raw rows; one per input point is carried to the seam copies. Once
    # ``uv`` holds the texture coordinates, ``vt`` rows that are not one per
    # point (a seamed file's) say nothing more and are left out.
    def raw_rows(rows):
        if uv is None or len(parent) == num_in or len(rows) != num_in:
            return rows
        return rows[parent]

    if len(texture_coords) > 0 and (uv is None or len(texture_coords) == num_in):
        point_data["obj:vt"] = raw_rows(texture_coords)
    if len(vertex_normals) > 0:
        point_data["obj:vn"] = raw_rows(vertex_normals)
    if uv is not None:
        point_data["uv"] = uv
        if len(parent) > num_in:
            point_data["obj:parent_point"] = np.asarray(parent, dtype=np.int64)
    cell_data = {"obj:group_ids": []}
    cells = []
    for f, gid in zip(face_groups, face_group_ids):
        conn = f if uv is not None else f - 1
        if f.shape[1] == 3:
            cells.append(CellBlock("triangle", conn))
        elif f.shape[1] == 4:
            cells.append(CellBlock("quad", conn))
        else:
            cells.append(CellBlock("polygon", conn))
        cell_data["obj:group_ids"].append(gid)

    return Mesh(points, cells, point_data=point_data, cell_data=cell_data)


def _vt_index(token):
    """The 0-based ``vt`` index of a face corner (``v/vt``, ``v/vt/vn``), or -1."""
    parts = token.split("/")
    if len(parts) < 2 or not parts[1]:
        return -1
    try:
        return int(parts[1]) - 1
    except ValueError:
        return -1


def _corner_uv(points, texture_coords, face_groups, face_vts):
    """``uv`` from the corners' ``vt`` indices, the twin of ``read_obj``'s.

    When every corner names a ``vt``, each point gets one ``uv``; a point whose
    corners name different ones (a seam) is split, the first copy keeping its
    index and the others appended in the order met. Rewrites ``face_groups`` to
    0-based ids into the split points and returns ``(uv, parent)``, or
    ``(None, None)`` (the faces left 1-based) without texture coordinates.
    """
    num_in = len(points)
    nvt = len(texture_coords)
    if nvt == 0 or texture_coords.ndim != 2 or texture_coords.shape[1] < 2:
        return None, None
    if not face_groups:
        return None, None
    for f, v in zip(face_groups, face_vts):
        vts = np.asarray(v)
        if np.any(vts < 0) or np.any(vts >= nvt) or np.any(f < 1) or np.any(f > num_in):
            return None, None
    copies = [[] for _ in range(num_in)]
    parent = list(range(num_in))
    uv = [[np.nan, np.nan] for _ in range(num_in)]
    for gi, (f, v) in enumerate(zip(face_groups, face_vts)):
        out = np.empty_like(f)
        for r in range(f.shape[0]):
            for k in range(f.shape[1]):
                p = int(f[r, k]) - 1
                t = int(v[r][k])
                ident = next((o for vt_at, o in copies[p] if vt_at == t), -1)
                if ident < 0:
                    if copies[p]:
                        ident = len(parent)
                        parent.append(p)
                        uv.append([0.0, 0.0])
                    else:
                        ident = p
                    copies[p].append((t, ident))
                    uv[ident] = [
                        float(texture_coords[t, 0]),
                        float(texture_coords[t, 1]),
                    ]
                out[r, k] = ident
        face_groups[gi] = out
    return np.array(uv, dtype=float).reshape(-1, 2), np.array(parent, dtype=np.int64)


def write(filename, mesh):
    for c in mesh.cells:
        if c.type not in ["triangle", "quad", "polygon"]:
            raise WriteError(
                "Wavefront .obj files can only contain triangle or quad cells."
            )

    with open_file(filename, "w") as f:
        f.write(_provenance.render_lines(_provenance.SlotTier.BLOCK, "# "))
        for p in mesh.points:
            f.write(f"v {p[0]} {p[1]} {p[2]}\n")

        if "obj:vn" in mesh.point_data:
            dat = mesh.point_data["obj:vn"]
            fmt = "vn " + " ".join(["{}"] * dat.shape[1]) + "\n"
            for vn in dat:
                f.write(fmt.format(*vn))

        # ``uv`` (one pair per point) becomes one ``vt`` per point, which every
        # face corner names (``f v/vt``); it replaces the raw ``obj:vt`` rows.
        uv = mesh.point_data.get("uv")
        with_uv = (
            uv is not None
            and np.ndim(uv) == 2
            and len(uv) == len(mesh.points)
            and np.shape(uv)[1] >= 2
        )
        if with_uv:
            for u, v in np.asarray(uv, dtype=float)[:, :2]:
                if not (np.isfinite(u) and np.isfinite(v)):
                    u = v = 0.0  # a point no face uses
                f.write(f"vt {u} {v}\n")
        elif "obj:vt" in mesh.point_data:
            dat = mesh.point_data["obj:vt"]
            fmt = "vt " + " ".join(["{}"] * dat.shape[1]) + "\n"
            for vt in dat:
                f.write(fmt.format(*vt))

        for cell_block in mesh.cells:
            if with_uv:
                for c in cell_block.data:
                    f.write("f " + " ".join(f"{i}/{i}" for i in c + 1) + "\n")
                continue
            fmt = "f " + " ".join(["{}"] * cell_block.data.shape[1]) + "\n"
            for c in cell_block.data:
                f.write(fmt.format(*(c + 1)))


# NOTE: format registration now lives in meshioplusplus/obj/__init__.py, which wraps the
# reader/writer below with the C++-backed fast paths.
