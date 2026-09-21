"""glTF 2.0 writer -- the Python reference, the twin of ``formats/gltf.cpp``.

A web-native surface export (``.glb`` binary container, or ``.gltf`` with a
``.bin`` beside it). It exports the *surface* of whatever it is given:

* the skin of volume blocks (linearised), with a 2-D cell that coincides with a
  skin facet winning over it, so a boundary patch such as ``wall`` keeps its own
  cells instead of z-fighting with the skin;
* 2-D cells, fanned like ``convert_cells(simplexify)`` and linearised;
* ``line`` cells as ``LINES`` and ``vertex`` cells (or a cell-less mesh, which is
  how a ``.pcd``/``.xyz`` point cloud arrives) as ``POINTS``.

glTF normals are per vertex, so a crease needs its position twice: triangles are
grouped into smooth fans exactly as :func:`meshioplusplus.compute_normals` does
(``_normals.vertex_normal_groups``) and each fan becomes its own vertex.

Everything follows the C++ writer step for step with the same arithmetic in the
same order, and the JSON is written by a hand-rolled serialiser with a fixed key
order and ``%.17g`` numbers, so the two engines write identical bytes
(``tests/python/test_gltf.py`` pins that).

Conventions: ``float32``, metres, right-handed, Y-up. Source data is usually
Z-up, so the axis change is a rotation on the root node rather than a rewrite of
the coordinates, and so is the recentring offset: the transform is exactly
reversible (a permutation with signs) and the coordinates stay small.

Fields: every one-to-four component ``point_data`` array is exported raw as an
underscore-prefixed custom attribute (``temperature`` -> ``_TEMPERATURE``),
cast to ``float32``. Naming one with ``color_by`` also bakes it through a
colormap into ``COLOR_0`` (linear, so the colormap's sRGB bytes go through the
inverse transfer function) and marks the material ``KHR_materials_unlit``.
"""

from __future__ import annotations

import math
import struct

import numpy as np

from .. import _provenance
from .._colormap import SRGB_TO_LINEAR_BITS, colormap_lookup, colormap_table
from .._common import warn
from .._exceptions import WriteError
from .._facecolor import _scalarize, color_param
from .._mesh import Mesh, topological_dimension
from .._normals import vertex_normal_groups
from .._regions import block_bases
from .._skin import _extract_skin_py, _has_skinnable_cells

_PREFIX = "meshio++: gltf: "
_SOURCE_POINT = "gltf:source_point"
_HALF_SQRT2 = 0.70710678118654757

_UINT32 = 5125
_FLOAT32 = 5126
_TARGET_ARRAY = 34962
_TARGET_ELEMENT_ARRAY = 34963
_MODE_POINTS, _MODE_LINES, _MODE_TRIANGLES = 0, 1, 4

_CONTAINERS = ("auto", "glb", "binary", "gltf", "json")
_UP_AXES = ("auto", "x", "y", "z")
_WEIGHTS = ("angle", "area")

_LINEAR = np.array(SRGB_TO_LINEAR_BITS, dtype=np.uint32).view(np.float32)


# --------------------------------------------------------------------------- #
# the JSON serialiser                                                          #
# --------------------------------------------------------------------------- #
_ESCAPES = {
    '"': '\\"',
    "\\": "\\\\",
    "\b": "\\b",
    "\f": "\\f",
    "\n": "\\n",
    "\r": "\\r",
    "\t": "\\t",
}


def _string(text: str) -> str:
    out = ['"']
    for ch in text:
        if ch in _ESCAPES:
            out.append(_ESCAPES[ch])
        elif ord(ch) < 0x20:
            out.append("\\u%04x" % ord(ch))
        else:
            out.append(ch)
    out.append('"')
    return "".join(out)


def _dumps(value) -> str:
    """Compact JSON in insertion order, floats as ``%.17g`` -- the C++ ``GltfJson``."""
    if value is True:
        return "true"
    if value is False:
        return "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return "%.17g" % value
    if isinstance(value, str):
        return _string(value)
    if isinstance(value, dict):
        return (
            "{" + ",".join(_string(k) + ":" + _dumps(v) for k, v in value.items()) + "}"
        )
    return "[" + ",".join(_dumps(v) for v in value) + "]"


# --------------------------------------------------------------------------- #
# options                                                                      #
# --------------------------------------------------------------------------- #
def _validate(
    container,
    up_axis,
    normal_weight,
    split_angle,
    normals,
    scale,
    color_by,
    cmap,
    vmin,
    vmax,
    nan_color,
):
    if container not in _CONTAINERS:
        raise ValueError(
            f"{_PREFIX}unknown container '{container}' (expected 'auto', 'glb' or 'gltf')"
        )
    if up_axis not in _UP_AXES:
        raise ValueError(
            f"{_PREFIX}unknown up axis '{up_axis}' (expected 'auto', 'x', 'y' or 'z')"
        )
    if normal_weight not in _WEIGHTS:
        raise ValueError(
            f"{_PREFIX}unknown weight '{normal_weight}' (expected 'angle' or 'area')"
        )
    if normals and not (0.0 <= float(split_angle) <= 180.0):
        raise ValueError(f"{_PREFIX}split_angle must lie in [0, 180] degrees")
    if not (float(scale) > 0.0) or not math.isfinite(float(scale)):
        raise ValueError(f"{_PREFIX}scale must be a positive number")
    if color_by:
        colormap_table(cmap)  # raises ValueError for an unknown name
        _parse_color(nan_color)
        if vmin is not None and vmax is not None and vmin > vmax:
            raise ValueError(f"{_PREFIX}vmin must not exceed vmax")


def _parse_color(text: str):
    ok = len(text) == 7 and text[0] == "#"
    if ok:
        try:
            int(text[1:], 16)
        except ValueError:
            ok = False
    if not ok or any(c not in "0123456789abcdefABCDEF" for c in text[1:]):
        raise ValueError(f"{_PREFIX}nan_color '{text}' is not a '#rrggbb' colour")
    return int(text[1:3], 16), int(text[3:5], 16), int(text[5:7], 16)


# --------------------------------------------------------------------------- #
# collecting the surface                                                       #
# --------------------------------------------------------------------------- #
def _cell_ids(data, c):
    return [int(x) for x in data[c]]


def _collect(mesh, bases):
    """The twin of ``gltf_collect``: triangles, lines and vertex cells."""
    tri_v, tri_facet, tri_cell = [], [], []
    line_v, line_cell = [], []
    vert_p, vert_cell = [], []
    flat = []  # (cell, ids)
    keys2d = set()
    skipped = []
    has_volume = False

    for b, cb in enumerate(mesh.cells):
        base = int(bases[b])
        name = cb.type
        if name.startswith("polyhedron") or topological_dimension.get(name) == 3:
            has_volume = True
            continue
        kind, corners = 0, 0
        if name.startswith("triangle"):
            corners = 3
        elif name.startswith("quad"):
            corners = 4
        elif name.startswith("polygon"):
            corners = 0
        elif name.startswith("line"):
            kind = 1
        elif name == "vertex":
            kind = 2
        else:
            skipped.append(name)
            continue
        for c in range(len(cb.data)):
            ids = _cell_ids(cb.data, c)
            cell = base + c
            if kind == 0:
                k = len(ids) if corners == 0 else min(corners, len(ids))
                if k < 3:
                    continue
                flat.append((cell, ids[:k]))
                if k in (3, 4):
                    keys2d.add(tuple(sorted(ids[:k])) + (-1,) * (4 - k))
            elif kind == 1:
                if len(ids) < 2 or ids[0] == ids[1]:
                    continue
                line_v.append((ids[0], ids[1]))
                line_cell.append(cell)
            else:
                if not ids:
                    continue
                vert_p.append(ids[0])
                vert_cell.append(cell)
    if skipped:
        joined = ", ".join(skipped)
        warn(f"gltf: cell block(s) of type {joined} have no glTF equivalent; skipping.")
        _provenance.note(
            "cells-dropped", f"cell block(s) of type {joined} have no glTF equivalent"
        )

    facet = 0
    if has_volume:
        if _has_skinnable_cells(mesh):
            vol = Mesh(
                mesh.points,
                list(mesh.cells),
                point_data={_SOURCE_POINT: np.arange(len(mesh.points), dtype=np.int64)},
            )
            skin = _extract_skin_py(vol, linearize=True, record_parent_ids=True)
            source = skin.point_data[_SOURCE_POINT]
            parents = np.concatenate(
                [
                    np.asarray(p).reshape(-1)
                    for p in skin.cell_data["surface:parent_cell"]
                ]
            )
            cell_index = 0
            for cb in skin.cells:
                data = np.asarray(cb.data)
                if cb.type not in ("triangle", "quad"):
                    cell_index += len(data)
                    facet += len(data)
                    continue
                npc = 4 if cb.type == "quad" else 3
                for row in data:
                    corner = [int(source[int(i)]) for i in row[:npc]]
                    parent = int(parents[cell_index])
                    cell_index += 1
                    this_facet = facet
                    facet += 1
                    key = tuple(sorted(corner)) + (-1,) * (4 - npc)
                    if key in keys2d:
                        continue
                    for k in range(1, npc - 1):
                        tri = (corner[0], corner[k], corner[k + 1])
                        if tri[0] == tri[1] or tri[1] == tri[2] or tri[0] == tri[2]:
                            continue
                        tri_v.append(tri)
                        tri_facet.append(this_facet)
                        tri_cell.append(parent)
        else:
            warn(
                "gltf: the volume cell blocks are not supported by the skin "
                "extractor; skipping."
            )

    for f, (cell, ids) in enumerate(flat):
        this_facet = facet + f
        for j in range(1, len(ids) - 1):
            tri = (ids[0], ids[j], ids[j + 1])
            if tri[0] == tri[1] or tri[1] == tri[2] or tri[0] == tri[2]:
                continue
            tri_v.append(tri)
            tri_facet.append(this_facet)
            tri_cell.append(cell)
    return tri_v, tri_facet, tri_cell, line_v, line_cell, vert_p, vert_cell


def _assign_nodes(mesh, total_cells, by_region):
    """The twin of ``gltf_assign_nodes``: (node names, node index per cell)."""
    kind_order = {"point": 0, "cell": 1, "side": 2}
    regions = sorted(
        (
            r
            for r in (getattr(mesh, "regions", ()) or ())
            if r.kind == "cell" and len(r.entries) > 0
        ),
        key=lambda r: (kind_order[r.kind], r.name, r.dim, r.tag),
    )
    if not by_region or not regions:
        return ["mesh"], np.zeros(total_cells, dtype=np.int64)

    names = []
    region_node = []
    for r in regions:
        if r.name not in names:
            names.append(r.name)
        region_node.append(names.index(r.name))

    best = np.full(total_cells, -1, dtype=np.int64)
    best_size = np.zeros(total_cells, dtype=np.int64)
    members = np.zeros(total_cells, dtype=np.int64)
    for i, r in enumerate(regions):
        entries = np.asarray(r.entries, dtype=np.int64).reshape(-1)
        entries = entries[(entries >= 0) & (entries < total_cells)]
        size = len(r.entries)
        members[entries] += 1
        better = (best[entries] < 0) | (size < best_size[entries])
        best[entries[better]] = i
        best_size[entries[better]] = size
    overlapping = int(np.count_nonzero(members > 1))
    if np.any(best < 0):
        names.append("unassigned")
    node = np.array(
        [len(names) - 1 if b < 0 else region_node[b] for b in best], dtype=np.int64
    )
    if overlapping > 0:
        _provenance.note(
            "regions-overlap",
            f"{overlapping} cell(s) are in more than one region and were exported "
            "under the smallest",
        )
    return names, node


# --------------------------------------------------------------------------- #
# primitives                                                                   #
# --------------------------------------------------------------------------- #
class _Prim:
    """Vertices numbered by first use, the twin of ``GltfPrimBuilder``."""

    def __init__(self, mode):
        self.mode = mode
        self.indices = []
        self.vertices = []
        self._lookup = {}

    def add(self, key):
        idx = self._lookup.get(key)
        if idx is None:
            idx = len(self.vertices)
            if idx >= 0xFFFFFFFF - 1:
                raise WriteError(
                    f"{_PREFIX}a primitive has more vertices than a 32-bit index can address"
                )
            self._lookup[key] = idx
            self.vertices.append(key)
        if self.mode != _MODE_POINTS:
            self.indices.append(idx)

    @property
    def empty(self):
        return not self.vertices


def _attr_name(name: str, taken: list) -> str:
    base = ["_"]
    for c in name.encode("utf-8").upper():
        ok = 65 <= c <= 90 or 48 <= c <= 57 or c == 95
        base.append(chr(c) if ok else "_")
    base = "".join(base)
    attr = base
    suffix = 2
    while attr in taken:
        attr = f"{base}_{suffix}"
        suffix += 1
    taken.append(attr)
    return attr


_TYPE_NAMES = {1: "SCALAR", 2: "VEC2", 3: "VEC3", 4: "VEC4"}


def _linear(byte: int) -> float:
    return float(_LINEAR[byte])


def _render(
    mesh,
    container_binary,
    bin_uri,
    up_axis,
    normal_weight,
    normals,
    fields_on,
    recenter,
    by_region,
    unlit,
    split_angle,
    scale,
    color_by,
    component,
    cmap,
    vmin,
    vmax,
    nan_color,
):
    points = np.asarray(mesh.points, dtype=np.float64)
    if points.ndim == 2 and points.shape[1] == 2:
        points = np.column_stack([points, np.zeros(len(points))])
    if points.ndim != 2 or points.shape[1] != 3:
        points = points.reshape(len(points), -1)
    n = len(points)
    dim = np.asarray(mesh.points).shape[1] if n else 3

    color_active = bool(color_by)
    table = None
    nan_rgb = (128, 128, 128)
    if color_active:
        table = colormap_table(cmap)
        nan_rgb = _parse_color(nan_color)
    color_point = color_active and color_by in mesh.point_data
    color_cell = color_active and not color_point and color_by in mesh.cell_data
    if color_active and not color_point and not color_cell:
        avail = sorted(set(mesh.point_data) | set(mesh.cell_data))
        raise ValueError(
            f"{_PREFIX}color_by array '{color_by}' is in neither point_data nor "
            f"cell_data (available: {', '.join(avail) if avail else 'none'})"
        )

    bases = block_bases(mesh.cells)
    total_cells = int(bases[-1]) + len(mesh.cells[-1].data) if len(mesh.cells) else 0
    tri_v, tri_facet, tri_cell, line_v, line_cell, vert_p, vert_cell = _collect(
        mesh, bases
    )
    node_names, cell_node = _assign_nodes(mesh, total_cells, by_region)

    # A non-finite coordinate on an exported point is an error; look before any
    # arithmetic touches it.
    exported = set(v for tri in tri_v for v in tri)
    exported.update(v for seg in line_v for v in seg)
    exported.update(vert_p)
    if len(mesh.cells) == 0:
        exported.update(range(n))
    if exported:
        idx = np.fromiter(sorted(exported), dtype=np.int64)
        bad = np.flatnonzero(~np.isfinite(points[idx]).all(axis=1))
        if len(bad):
            raise WriteError(
                f"{_PREFIX}point {int(idx[bad[0]])} has a non-finite coordinate"
            )

    ntri = len(tri_v)
    groups = None
    if normals and ntri:
        verts = np.asarray(tri_v, dtype=np.int64).reshape(-1, 3)
        corners = points[verts]
        groups = vertex_normal_groups(
            verts,
            corners,
            np.asarray(tri_facet, dtype=np.int64),
            n,
            normal_weight,
            float(split_angle),
        )
        corner_group = groups["corner_group"]

    builders = [
        {
            "tri": _Prim(_MODE_TRIANGLES),
            "line": _Prim(_MODE_LINES),
            "point": _Prim(_MODE_POINTS),
        }
        for _ in node_names
    ]
    for t in range(ntri):
        b = builders[int(cell_node[tri_cell[t]])]["tri"]
        for i in range(3):
            b.add(
                (
                    tri_v[t][i],
                    int(corner_group[t * 3 + i]) if normals else -1,
                    tri_cell[t] if color_cell else -1,
                )
            )
    for s in range(len(line_v)):
        b = builders[int(cell_node[line_cell[s]])]["line"]
        for v in line_v[s]:
            b.add((v, -1, line_cell[s] if color_cell else -1))
    for s in range(len(vert_p)):
        builders[int(cell_node[vert_cell[s]])]["point"].add(
            (vert_p[s], -1, vert_cell[s] if color_cell else -1)
        )
    if len(mesh.cells) == 0:
        for p in range(n):
            builders[0]["point"].add((p, -1, -1))

    prims_of = lambda nb: (nb["tri"], nb["line"], nb["point"])  # noqa: E731
    used = np.zeros(n, dtype=bool)
    for nb in builders:
        for pb in prims_of(nb):
            for key in pb.vertices:
                used[key[0]] = True

    any_used = bool(used.any())
    flat = dim < 3
    if any_used:
        up_pts = points[used]
        if not np.all(np.isfinite(up_pts)):
            bad = int(np.flatnonzero(used & ~np.isfinite(points).all(axis=1))[0])
            raise WriteError(f"{_PREFIX}point {bad} has a non-finite coordinate")
        lo = up_pts.min(axis=0)
        hi = up_pts.max(axis=0)
        if np.any(np.abs(up_pts[:, 2]) > 1e-14):
            flat = False
    center = [0.0, 0.0, 0.0]
    if recenter and any_used:
        center = [0.5 * (float(lo[k]) + float(hi[k])) for k in range(3)]

    # colour values and their range over the exported vertices
    values = None
    lo_v = hi_v = 0.0
    if color_active:
        if color_point:
            arr = np.asarray(mesh.point_data[color_by])
            ncomp = 0 if n == 0 else arr.size // n
            flatarr = arr.reshape(-1).tolist()
            values = [float("nan")] * n
            for p in range(n):
                if used[p]:
                    values[p] = _scalarize(flatarr, p, ncomp, component)
        else:
            blocks = mesh.cell_data[color_by]
            if len(blocks) != len(mesh.cells):
                raise ValueError(
                    f"{_PREFIX}cell_data array '{color_by}' has {len(blocks)} block(s) "
                    f"but the mesh has {len(mesh.cells)}"
                )
            values = []
            for blk in blocks:
                a = np.asarray(blk)
                nc = len(a)
                ncomp = 0 if nc == 0 else a.size // nc
                fl = a.reshape(-1).tolist()
                for r in range(nc):
                    values.append(_scalarize(fl, r, ncomp, component))
        seen = False
        for nb in builders:
            for pb in prims_of(nb):
                for key in pb.vertices:
                    v = values[key[2] if color_cell else key[0]]
                    if not math.isfinite(v):
                        continue
                    if not seen:
                        lo_v = hi_v = v
                        seen = True
                    else:
                        lo_v = v if v < lo_v else lo_v
                        hi_v = v if v > hi_v else hi_v
        lo_v = float(vmin) if vmin is not None else lo_v
        hi_v = float(vmax) if vmax is not None else hi_v
        if lo_v > hi_v:
            raise ValueError(f"{_PREFIX}vmin must not exceed vmax")

    def vertex_color(key):
        v = values[key[2] if color_cell else key[0]]
        rgb = nan_rgb
        if math.isfinite(v):
            rgb = colormap_lookup(table, color_param(v, lo_v, hi_v))
        return _linear(rgb[0]), _linear(rgb[1]), _linear(rgb[2])

    # raw fields
    fields = []
    if fields_on:
        taken = []
        for name in sorted(mesh.point_data):
            if name == "normals":
                continue
            arr = np.asarray(mesh.point_data[name])
            ncomp = 0
            if n > 0 and arr.ndim in (1, 2) and arr.shape[0] == n:
                ncomp = 1 if arr.ndim == 1 else arr.shape[1]
            if ncomp < 1 or ncomp > 4:
                warn(
                    f"gltf: point_data '{name}' has no glTF attribute type (1 to 4 "
                    "components per point); skipping."
                )
                _provenance.note(
                    "fields-dropped",
                    f"point_data '{name}' is not a 1 to 4 component array",
                )
                continue
            data = (
                np.asarray(arr, dtype=np.float64).reshape(n, ncomp).astype(np.float32)
            )
            if not np.all(np.isfinite(data[used])):
                warn(
                    f"gltf: point_data '{name}' has non-finite values, which glTF "
                    "cannot hold; skipping."
                )
                _provenance.note(
                    "fields-dropped", f"point_data '{name}' has non-finite values"
                )
                continue
            fields.append((name, _attr_name(name, taken), ncomp, data))

    # point normals for a POINTS primitive
    point_normals = None
    if normals and "normals" in mesh.point_data:
        arr = np.asarray(mesh.point_data["normals"])
        if n > 0 and arr.ndim == 2 and arr.shape == (n, 3):
            arr = arr.astype(np.float64)
            length = np.sqrt(
                arr[:, 0] * arr[:, 0] + arr[:, 1] * arr[:, 1] + arr[:, 2] * arr[:, 2]
            )
            ok = np.isfinite(length) & (length > 0.0)
            unit = np.zeros_like(arr)
            unit[ok] = arr[ok] * (1.0 / length[ok])[:, None]
            point_normals = (unit, ok)

    # ---- the binary chunk and the accessors ----
    bin_parts = []
    bin_size = 0
    accessors = []

    def add_floats(data, ncomp, minmax, name):
        nonlocal bin_size
        data = np.ascontiguousarray(data, dtype="<f4").reshape(-1, ncomp)
        acc = {
            "ct": _FLOAT32,
            "count": len(data),
            "ncomp": ncomp,
            "name": name,
            "target": _TARGET_ARRAY,
            "offset": bin_size,
            "length": data.size * 4,
        }
        if minmax and len(data):
            acc["max"] = [float(x) for x in data.max(axis=0)]
            acc["min"] = [float(x) for x in data.min(axis=0)]
        bin_parts.append(data.tobytes())
        bin_size += data.size * 4
        accessors.append(acc)
        return len(accessors) - 1

    def add_indices(idx):
        nonlocal bin_size
        data = np.asarray(idx, dtype="<u4")
        accessors.append(
            {
                "ct": _UINT32,
                "count": len(data),
                "ncomp": 1,
                "name": "",
                "target": _TARGET_ELEMENT_ARRAY,
                "offset": bin_size,
                "length": data.size * 4,
            }
        )
        bin_parts.append(data.tobytes())
        bin_size += data.size * 4
        return len(accessors) - 1

    out_nodes = []
    center_arr = np.array(center)
    for ni, nb in enumerate(builders):
        prim_outs = []
        for pb in prims_of(nb):
            if pb.empty:
                continue
            nv = len(pb.vertices)
            pidx = np.array([k[0] for k in pb.vertices], dtype=np.int64)
            po = {"mode": pb.mode, "indices": -1, "attrs": []}
            if pb.mode != _MODE_POINTS:
                po["indices"] = add_indices(pb.indices)
            position = (points[pidx] - center_arr).astype(np.float32)
            po["attrs"].append(("POSITION", add_floats(position, 3, True, "")))

            if pb.mode == _MODE_TRIANGLES and normals:
                gid = np.array([k[1] for k in pb.vertices], dtype=np.int64)
                gn = groups["group_normal"][gid]
                defined = (gn != 0.0).any(axis=1)
                nrm = np.where(defined[:, None], gn, np.array([0.0, 0.0, 1.0]))
                po["attrs"].append(
                    ("NORMAL", add_floats(nrm.astype(np.float32), 3, False, ""))
                )
            elif pb.mode == _MODE_POINTS and point_normals is not None:
                unit, ok = point_normals
                if bool(np.all(ok[pidx])):
                    po["attrs"].append(
                        (
                            "NORMAL",
                            add_floats(unit[pidx].astype(np.float32), 3, False, ""),
                        )
                    )
                else:
                    warn(
                        "gltf: point_data 'normals' has a zero or non-finite row on the "
                        "exported points; NORMAL omitted."
                    )
            if color_active:
                rgb = np.array([vertex_color(k) for k in pb.vertices], dtype=np.float32)
                po["attrs"].append(("COLOR_0", add_floats(rgb, 3, False, "")))
            for name, attr, ncomp, data in fields:
                po["attrs"].append((attr, add_floats(data[pidx], ncomp, False, name)))
            prim_outs.append(po)
            del nv
        if prim_outs:
            out_nodes.append((node_names[ni], prim_outs))
    if not out_nodes:
        warn(
            "gltf: nothing to export (no surface, line or vertex cells); writing an empty scene."
        )

    # ---- the root transform ----
    up = up_axis
    if up == "auto":
        up = "y" if flat else "z"
    s = float(scale)
    sx, sy, sz = s * center[0], s * center[1], s * center[2]
    if up == "z":
        translation = [sx, sz, -sy]
    elif up == "x":
        translation = [-sy, sx, sz]
    else:
        translation = [sx, sy, sz]
    translation = [t + 0.0 for t in translation]

    # ---- the JSON ----
    provenance = _provenance.lines(_provenance.SlotTier.BLOCK)
    asset = {
        "generator": provenance[0] if provenance else _provenance.TAG,
        "version": "2.0",
    }
    if len(provenance) > 1:
        asset["extras"] = {"meshioplusplus:provenance": list(provenance)}
    doc = {"asset": asset}
    has_unlit = color_active and unlit
    if has_unlit:
        doc["extensionsUsed"] = ["KHR_materials_unlit"]
    doc["scene"] = 0
    doc["scenes"] = [{"nodes": [0]} if out_nodes else {}]
    if out_nodes:
        root = {"name": "meshio++", "children": list(range(1, len(out_nodes) + 1))}
        if up == "z":
            root["rotation"] = [-_HALF_SQRT2, 0.0, 0.0, _HALF_SQRT2]
        elif up == "x":
            root["rotation"] = [0.0, 0.0, _HALF_SQRT2, _HALF_SQRT2]
        if s != 1.0:
            root["scale"] = [s, s, s]
        root["translation"] = translation
        doc["nodes"] = [root] + [
            {"name": nm, "mesh": i} for i, (nm, _p) in enumerate(out_nodes)
        ]
        meshes = []
        for nm, prim_outs in out_nodes:
            prims = []
            for po in prim_outs:
                prim = {"attributes": {sem: idx for sem, idx in po["attrs"]}}
                if po["indices"] >= 0:
                    prim["indices"] = po["indices"]
                prim["material"] = 0
                prim["mode"] = po["mode"]
                prims.append(prim)
            meshes.append({"name": nm, "primitives": prims})
        doc["meshes"] = meshes
        base_c = 1.0 if color_active else 0.8
        material = {
            "name": "meshio++",
            "pbrMetallicRoughness": {
                "baseColorFactor": [base_c, base_c, base_c, 1.0],
                "metallicFactor": 0.0,
                "roughnessFactor": 0.5,
            },
        }
        if has_unlit:
            material["extensions"] = {"KHR_materials_unlit": {}}
        material["doubleSided"] = True
        doc["materials"] = [material]
        acc_json = []
        for i, a in enumerate(accessors):
            aj = {
                "bufferView": i,
                "componentType": a["ct"],
                "count": a["count"],
                "type": _TYPE_NAMES[a["ncomp"]],
            }
            if "max" in a:
                aj["max"] = a["max"]
                aj["min"] = a["min"]
            if a["name"]:
                aj["name"] = a["name"]
            acc_json.append(aj)
        doc["accessors"] = acc_json
        doc["bufferViews"] = [
            {
                "buffer": 0,
                "byteOffset": a["offset"],
                "byteLength": a["length"],
                "target": a["target"],
            }
            for a in accessors
        ]
        buf = {"byteLength": bin_size}
        if not container_binary:
            buf["uri"] = bin_uri
        doc["buffers"] = [buf]
    return _dumps(doc).encode("utf-8"), b"".join(bin_parts)


def _glb(json_bytes: bytes, bin_bytes: bytes) -> bytes:
    """The GLB container: header, a space-padded JSON chunk, a zero-padded BIN chunk."""
    json_bytes += b" " * (-len(json_bytes) % 4)
    bin_bytes += b"\0" * (-len(bin_bytes) % 4)
    total = 12 + 8 + len(json_bytes) + (8 + len(bin_bytes) if bin_bytes else 0)
    if total > 0xFFFFFFFF:
        raise WriteError(f"{_PREFIX}the file is larger than a GLB can hold (4 GiB)")
    out = struct.pack("<III", 0x46546C67, 2, total)
    out += struct.pack("<II", len(json_bytes), 0x4E4F534A) + json_bytes
    if bin_bytes:
        out += struct.pack("<II", len(bin_bytes), 0x004E4942) + bin_bytes
    return out


def _quote_uri(name: str) -> str:
    out = []
    for c in name.encode("utf-8"):
        if 65 <= c <= 90 or 97 <= c <= 122 or 48 <= c <= 57 or c in (45, 46, 95, 126):
            out.append(chr(c))
        else:
            out.append("%%%02X" % c)
    return "".join(out)


def write(
    filename,
    mesh,
    container="auto",
    up_axis="auto",
    normal_weight="angle",
    normals=True,
    fields=True,
    recenter=True,
    by_region=True,
    unlit=True,
    split_angle=30.0,
    scale=1.0,
    color_by=None,
    component=None,
    cmap="viridis",
    vmin=None,
    vmax=None,
    nan_color="#808080",
):
    """Write a glTF 2.0 file; see the module docstring and :func:`gltf.write`."""
    _validate(
        container,
        up_axis,
        normal_weight,
        split_angle,
        normals,
        scale,
        color_by,
        cmap,
        vmin,
        vmax,
        nan_color,
    )
    path = str(filename)
    lower = path.lower()
    binary = container in ("glb", "binary") or (
        container == "auto" and lower.endswith(".glb")
    )
    bin_path = ""
    if not binary:
        bin_path = (path[:-5] if lower.endswith(".gltf") else path) + ".bin"
    bin_name = bin_path.replace("\\", "/").rsplit("/", 1)[-1]
    json_bytes, bin_bytes = _render(
        mesh,
        binary,
        _quote_uri(bin_name),
        up_axis,
        normal_weight,
        normals,
        fields,
        recenter,
        by_region,
        unlit,
        split_angle,
        scale,
        color_by or "",
        component,
        cmap,
        vmin,
        vmax,
        nan_color,
    )
    if binary:
        with open(path, "wb") as f:
            f.write(_glb(json_bytes, bin_bytes))
        return
    with open(path, "wb") as f:
        f.write(json_bytes)
    if bin_bytes:
        with open(bin_path, "wb") as f:
            f.write(bin_bytes)


def write_buffer(
    buffer,
    mesh,
    up_axis="auto",
    normal_weight="angle",
    normals=True,
    fields=True,
    recenter=True,
    by_region=True,
    unlit=True,
    split_angle=30.0,
    scale=1.0,
    color_by=None,
    component=None,
    cmap="viridis",
    vmin=None,
    vmax=None,
    nan_color="#808080",
):
    """Write the GLB container into a binary file-like object."""
    _validate(
        "glb",
        up_axis,
        normal_weight,
        split_angle,
        normals,
        scale,
        color_by,
        cmap,
        vmin,
        vmax,
        nan_color,
    )
    json_bytes, bin_bytes = _render(
        mesh,
        True,
        "",
        up_axis,
        normal_weight,
        normals,
        fields,
        recenter,
        by_region,
        unlit,
        split_angle,
        scale,
        color_by or "",
        component,
        cmap,
        vmin,
        vmax,
        nan_color,
    )
    buffer.write(_glb(json_bytes, bin_bytes))
