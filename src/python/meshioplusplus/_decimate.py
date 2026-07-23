"""Surface decimation (``meshioplusplus.decimate``).

A dependency-free mesh *operation* (not a file format): reduce a surface mesh's
face count by greedy quadric-error-metric (Garland-Heckbert) edge collapse,
preserving shape, boundaries and features. The resolution-*reducing* inverse of
``refine``.

The C++ core (``_core.decimate``) is used when available and this module's
reference (``_decimate_py``) is the fallback. The fallback is a **byte-parity
twin**: the setup is vectorized numpy with the accumulation order pinned
(``np.add.at``, whose unbuffered index-order scatter reproduces the C++
per-vertex ascending-face quadric sums exactly), and the greedy loop is pure
Python driven by ``heapq`` on the same ``(error, lower id, higher id,
versions)`` tuples the C++ priority queue orders by. Every floating-point
expression that feeds a branch (plane, quadric, 3x3 solve, error form, blend
parameter) is transcribed token for token from
``src/cpp/src/operations/decimate.cpp`` -- KEEP THEM IN SYNC -- so both sides take
the same side of every discrete decision and the output is byte-identical
(pinned by ``tests/python/test_decimate.py::test_cpp_matches_python``).

``point_sets`` are remapped through the returned point map (each entry lands on
its survivor; duplicates collapse) and ``cell_sets`` through the per-block cell
maps (removed faces drop out); neither ever reaches the C++ core.
"""

import heapq
import math
from bisect import bisect_left, insort

import numpy as np

from ._convert_cells import _simplexify_py
from ._mesh import Mesh, topological_dimension

_PLACEMENTS = ("optimal", "midpoint", "endpoint")


# --------------------------------------------------------------------------- #
# validation (mirrors decim_resolve_stop / decim_check_blocks)                 #
# --------------------------------------------------------------------------- #
def _check_criteria(target_ratio, target_faces, max_error):
    num_set = (
        (1 if target_ratio >= 0.0 else 0)
        + (1 if target_faces >= 0 else 0)
        + (1 if max_error >= 0.0 else 0)
    )
    if num_set != 1:
        raise ValueError(
            "meshio++: decimate: exactly one of target_ratio, target_faces, "
            f"max_error must be set; got {num_set}"
        )
    if target_ratio >= 0.0 and not 0.0 < target_ratio <= 1.0:
        raise ValueError(
            "meshio++: decimate: target_ratio must lie in (0, 1]; got "
            f"{target_ratio:f}"
        )


def _check_blocks(mesh):
    has_surface = False
    for cb in mesh.cells:
        if cb.type.startswith("polyhedron"):
            raise ValueError(
                "meshio++: decimate: mesh contains a polyhedron cell block; "
                "decimate operates on surface meshes (extract_surface first)"
            )
        dim = topological_dimension.get(cb.type, -1)
        if dim == 3:
            raise ValueError(
                f"meshio++: decimate: mesh contains 3D volume cell block "
                f"'{cb.type}'; decimate operates on surface meshes "
                "(extract_surface first)"
            )
        if dim == 2:
            has_surface = True
    if not has_surface:
        raise ValueError("meshio++: decimate: mesh contains no surface (2D) cell block")
    for cb in mesh.cells:
        dim = topological_dimension.get(cb.type, -1)
        if dim <= 1:
            raise ValueError(
                f"meshio++: decimate: cannot decimate a mesh containing cell "
                f"block '{cb.type}' alongside its surface blocks (its nodes "
                "would dangle after the collapse; drop it first, e.g. via split)"
            )
        data = np.asarray(cb.data)
        if data.dtype == object:
            raise ValueError(
                f"meshio++: decimate: cannot decimate ragged cell block "
                f"'{cb.type}' (triangulate/clean it first)"
            )
        if cb.type not in ("triangle", "quad", "polygon"):
            raise ValueError(
                f"meshio++: decimate: cannot decimate higher-order cell block "
                f"'{cb.type}' (linearize the mesh first)"
            )


# --------------------------------------------------------------------------- #
# pinned floating-point kernels (token-for-token twins of decimate.cpp)        #
# --------------------------------------------------------------------------- #
def _quadric_error(q, x, y, z):
    return (
        q[0] * x * x
        + q[4] * y * y
        + q[7] * z * z
        + 2.0 * (q[1] * x * y + q[2] * x * z + q[5] * y * z)
        + 2.0 * (q[3] * x + q[6] * y + q[8] * z)
        + q[9]
    )


def _place(xyz, quads, pinned, placement, a, b):
    """Twin of decim_place: (x, y, z, err) of collapsing edge (a, b), a < b."""
    ao = a * 3
    bo = b * 3
    qa = quads[a * 10 : a * 10 + 10]
    qb = quads[b * 10 : b * 10 + 10]
    q = [qa[i] + qb[i] for i in range(10)]

    if pinned[a]:
        x, y, z = xyz[ao], xyz[ao + 1], xyz[ao + 2]
        return x, y, z, _quadric_error(q, x, y, z)
    if pinned[b]:
        x, y, z = xyz[bo], xyz[bo + 1], xyz[bo + 2]
        return x, y, z, _quadric_error(q, x, y, z)

    if placement == "endpoint":
        err_a = _quadric_error(q, xyz[ao], xyz[ao + 1], xyz[ao + 2])
        err_b = _quadric_error(q, xyz[bo], xyz[bo + 1], xyz[bo + 2])
        if err_b < err_a:  # tie -> a, the lower id
            return xyz[bo], xyz[bo + 1], xyz[bo + 2], err_b
        return xyz[ao], xyz[ao + 1], xyz[ao + 2], err_a

    if placement == "optimal":
        c00 = q[4] * q[7] - q[5] * q[5]
        c01 = q[2] * q[5] - q[1] * q[7]
        c02 = q[1] * q[5] - q[2] * q[4]
        det = q[0] * c00 + q[1] * c01 + q[2] * c02
        scale = max(abs(q[0]), abs(q[1]), abs(q[2]), abs(q[4]), abs(q[5]), abs(q[7]))
        if abs(det) > 1e-12 * (scale * scale * scale):
            c11 = q[0] * q[7] - q[2] * q[2]
            c12 = q[1] * q[2] - q[0] * q[5]
            c22 = q[0] * q[4] - q[1] * q[1]
            inv = 1.0 / det
            x = -(c00 * q[3] + c01 * q[6] + c02 * q[8]) * inv
            y = -(c01 * q[3] + c11 * q[6] + c12 * q[8]) * inv
            z = -(c02 * q[3] + c12 * q[6] + c22 * q[8]) * inv
            return x, y, z, _quadric_error(q, x, y, z)

    # midpoint, or optimal's ill-conditioned fallback
    x = (xyz[ao] + xyz[bo]) * 0.5
    y = (xyz[ao + 1] + xyz[bo + 1]) * 0.5
    z = (xyz[ao + 2] + xyz[bo + 2]) * 0.5
    return x, y, z, _quadric_error(q, x, y, z)


def _face_normal(xyz, corners, a, b, sub):
    """Twin of decim_face_normal: unnormalized normal, endpoints read as sub."""
    p = []
    for cid in corners:
        if sub is not None and (cid == a or cid == b):
            p.append(sub)
        else:
            o = cid * 3
            p.append((xyz[o], xyz[o + 1], xyz[o + 2]))
    e1x = p[1][0] - p[0][0]
    e1y = p[1][1] - p[0][1]
    e1z = p[1][2] - p[0][2]
    e2x = p[2][0] - p[0][0]
    e2y = p[2][1] - p[0][1]
    e2z = p[2][2] - p[0][2]
    return (
        e1y * e2z - e1z * e2y,
        e1z * e2x - e1x * e2z,
        e1x * e2y - e1y * e2x,
    )


def _vertex_ring(vfaces, corners, v):
    ring = set()
    for f in vfaces:
        o = f * 3
        for cid in (corners[o], corners[o + 1], corners[o + 2]):
            if cid != v:
                ring.add(cid)
    return sorted(ring)


# --------------------------------------------------------------------------- #
# the numpy/pure-Python reference                                              #
# --------------------------------------------------------------------------- #
def _decimate_py(
    mesh,
    target_ratio,
    target_faces,
    max_error,
    placement,
    preserve_boundary,
    preserve_features,
    feature_angle,
    frozen_ids,
):
    _check_criteria(target_ratio, target_faces, max_error)
    if placement not in _PLACEMENTS:
        raise ValueError(
            f"meshio++: decimate: unknown placement '{placement}' "
            "(expected 'optimal', 'midpoint' or 'endpoint')"
        )
    n = len(mesh.points)
    _check_blocks(mesh)

    # Triangulate. Higher-order input was rejected, so points are untouched.
    simp, _, _ = _simplexify_py(mesh, record_parent_ids=True)

    dim = np.asarray(simp.points).shape[1] if n else 3
    pts3 = np.zeros((n, 3), dtype=np.float64)
    if n:
        pts3[:, : min(dim, 3)] = np.asarray(simp.points, dtype=np.float64)[
            :, : min(dim, 3)
        ]

    # Flat block-major face table.
    block_cells = [np.asarray(cb.data, dtype=np.int64) for cb in simp.cells]
    block_base = []
    base = 0
    for data in block_cells:
        block_base.append(base)
        base += len(data)
    corners_np = (
        np.concatenate([d.reshape(-1, 3) for d in block_cells])
        if base
        else np.empty((0, 3), dtype=np.int64)
    )
    if corners_np.size and (corners_np.min() < 0 or corners_np.max() >= n):
        bad = corners_np.min() if corners_np.min() < 0 else corners_np.max()
        raise ValueError(
            f"meshio++: decimate: connectivity references point {bad} "
            f"outside [0, {n})"
        )
    num_faces = len(corners_np)

    # Vertex -> incident faces, ascending face index per vertex (counting sort).
    flat_verts = corners_np.reshape(-1)
    order = np.argsort(flat_verts, kind="stable")
    faces_of = order // 3
    counts = np.bincount(flat_verts, minlength=n) if num_faces else np.zeros(n, int)
    xadj = np.concatenate([[0], np.cumsum(counts)])

    # Edges + boundary marks (once-used-edge test). The edge SET is what
    # matters -- pop order comes from the heap's total order, not push order.
    if num_faces:
        ek = np.stack(
            [
                corners_np[:, [0, 1]],
                corners_np[:, [1, 2]],
                corners_np[:, [2, 0]],
            ],
            axis=1,
        ).reshape(-1, 2)
        ek = np.sort(ek, axis=1)
        ek = ek[ek[:, 0] != ek[:, 1]]
        edges, edge_counts = np.unique(ek, axis=0, return_counts=True)
    else:
        edges = np.empty((0, 2), dtype=np.int64)
        edge_counts = np.empty(0, dtype=np.int64)
    boundary = np.zeros(n, dtype=bool)
    once = edges[edge_counts == 1]
    boundary[once.reshape(-1)] = True

    # Per-face area-weighted plane quadrics + unit normals; expressions are the
    # elementwise twins of decim_face_planes.
    if num_faces:
        p0 = pts3[corners_np[:, 0]]
        p1 = pts3[corners_np[:, 1]]
        p2 = pts3[corners_np[:, 2]]
        e1 = p1 - p0
        e2 = p2 - p0
        nx = e1[:, 1] * e2[:, 2] - e1[:, 2] * e2[:, 1]
        ny = e1[:, 2] * e2[:, 0] - e1[:, 0] * e2[:, 2]
        nz = e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0]
        len2 = nx * nx + ny * ny + nz * nz
        nonzero = len2 != 0.0
        length = np.zeros(num_faces)
        length[nonzero] = np.sqrt(len2[nonzero])
        a_n = np.zeros(num_faces)
        b_n = np.zeros(num_faces)
        c_n = np.zeros(num_faces)
        a_n[nonzero] = nx[nonzero] / length[nonzero]
        b_n[nonzero] = ny[nonzero] / length[nonzero]
        c_n[nonzero] = nz[nonzero] / length[nonzero]
        d_n = -(a_n * p0[:, 0] + b_n * p0[:, 1] + c_n * p0[:, 2])
        d_n[~nonzero] = 0.0
        w = 0.5 * length
        wa = w * a_n
        wb = w * b_n
        wc = w * c_n
        wd = w * d_n
        K = np.stack(
            [
                wa * a_n,
                wa * b_n,
                wa * c_n,
                wa * d_n,
                wb * b_n,
                wb * c_n,
                wb * d_n,
                wc * c_n,
                wc * d_n,
                wd * d_n,
            ],
            axis=1,
        )
        normals = np.stack([a_n, b_n, c_n], axis=1)
    else:
        K = np.empty((0, 10))
        normals = np.empty((0, 3))

    # Pin mask: caller | boundary | features.
    pinned = np.zeros(n, dtype=bool)
    if frozen_ids is not None:
        ids = np.asarray(frozen_ids, dtype=np.int64).reshape(-1)
        if ids.size and (ids.min() < 0 or ids.max() >= n):
            bad = ids.min() if ids.min() < 0 else ids.max()
            raise ValueError(
                f"meshio++: decimate: frozen node id {bad} is out of range"
            )
        pinned[ids] = True
    if preserve_boundary:
        pinned |= boundary
    if preserve_features:
        cos_thr = math.cos(feature_angle * 3.14159265358979323846 / 180.0)
        for v in range(n):
            row = faces_of[xadj[v] : xadj[v + 1]]
            nrm = [normals[f] for f in row if not (normals[f] == 0.0).all()]
            done = False
            for i in range(len(nrm)):
                for j in range(i + 1, len(nrm)):
                    dot = (
                        nrm[i][0] * nrm[j][0]
                        + nrm[i][1] * nrm[j][1]
                        + nrm[i][2] * nrm[j][2]
                    )
                    if dot < cos_thr:
                        pinned[v] = True
                        done = True
                        break
                if done:
                    break

    # Per-vertex quadrics: np.add.at is unbuffered and processes indices in
    # order, so each vertex accumulates its faces ascending -- the same FP
    # order as the C++ per-vertex sums. Never a pairwise reduction here.
    Q = np.zeros((n, 10), dtype=np.float64)
    if num_faces:
        np.add.at(Q, flat_verts, np.repeat(K, 3, axis=0))

    # From here on: pure Python state, scalar IEEE doubles.
    xyz = [float(v) for v in pts3.reshape(-1)]
    quads = [float(v) for v in Q.reshape(-1)]
    pinned_l = [bool(v) for v in pinned]
    boundary_l = [bool(v) for v in boundary]
    corners = [int(v) for v in corners_np.reshape(-1)]
    face_alive = [True] * num_faces
    successor = [-1] * n
    version = [0] * n
    vfaces = []
    for v in range(n):
        row = []
        for f in faces_of[xadj[v] : xadj[v + 1]]:
            f = int(f)
            if row and row[-1] == f:
                continue  # a degenerate face repeating this corner
            row.append(f)
        vfaces.append(row)

    float_data = []  # (name, ncomp, flat float64 list)
    for name in sorted(simp.point_data):
        arr = np.asarray(simp.point_data[name])
        if not np.issubdtype(arr.dtype, np.floating) or len(arr) != n or n == 0:
            continue
        ncomp = arr.size // n
        float_data.append(
            (name, ncomp, [float(v) for v in arr.astype(np.float64).reshape(-1)])
        )

    use_target = target_ratio >= 0.0 or target_faces >= 0
    if target_ratio >= 0.0:
        stop_faces = max(1, int(target_ratio * float(num_faces) + 0.5))
    else:
        stop_faces = int(target_faces)

    heap = []
    for a, b in edges:
        a = int(a)
        b = int(b)
        if pinned_l[a] and pinned_l[b]:
            continue
        err = _place(xyz, quads, pinned_l, placement, a, b)[3]
        if not math.isfinite(err):
            continue
        heap.append((err, a, b, 0, 0))
    heapq.heapify(heap)

    faces_removed = 0
    collapses_rejected = 0
    max_error_applied = 0.0
    alive_faces = num_faces
    while True:
        if use_target and alive_faces <= stop_faces:
            break
        if not heap:
            break
        err, a, b, va, vb = heapq.heappop(heap)
        if max_error >= 0.0 and err > max_error:
            break
        if va != version[a] or vb != version[b]:
            continue  # stale
        if pinned_l[a] and pinned_l[b]:
            continue  # defensive

        # Guard 1: link condition (twin of the C++ guard block).
        fa = vfaces[a]
        fb = vfaces[b]
        shared = []
        i = j = 0
        while i < len(fa) and j < len(fb):
            if fa[i] < fb[j]:
                i += 1
            elif fb[j] < fa[i]:
                j += 1
            else:
                shared.append(fa[i])
                i += 1
                j += 1
        if not shared:
            continue
        if len(shared) > 2:
            collapses_rejected += 1
            continue
        if boundary_l[a] and boundary_l[b] and len(shared) == 2:
            collapses_rejected += 1
            continue
        ring_a = _vertex_ring(fa, corners, a)
        ring_b = _vertex_ring(fb, corners, b)
        common = 0
        i = j = 0
        while i < len(ring_a) and j < len(ring_b):
            if ring_a[i] < ring_b[j]:
                i += 1
            elif ring_b[j] < ring_a[i]:
                j += 1
            else:
                common += 1
                i += 1
                j += 1
        if common != len(shared):
            collapses_rejected += 1
            continue

        # Guard 2: normal flip / degenerate survivor, ascending face index.
        x, y, z, _ = _place(xyz, quads, pinned_l, placement, a, b)
        sub = (x, y, z)
        flips = False
        merged = sorted(fa + fb)
        si = 0
        for f in merged:
            while si < len(shared) and shared[si] < f:
                si += 1
            if si < len(shared) and shared[si] == f:
                continue  # dies in the collapse: no constraint
            o = f * 3
            tri = (corners[o], corners[o + 1], corners[o + 2])
            n0 = _face_normal(xyz, tri, a, b, None)
            if n0[0] == 0.0 and n0[1] == 0.0 and n0[2] == 0.0:
                continue
            n1 = _face_normal(xyz, tri, a, b, sub)
            if n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2] <= 0.0:
                flips = True
                break
        if flips:
            collapses_rejected += 1
            continue

        # Commit: r merges into the survivor s.
        s = b if pinned_l[b] else a
        r = b if s == a else a

        so = s * 3
        ro = r * 3
        dx = xyz[ro] - xyz[so]
        dy = xyz[ro + 1] - xyz[so + 1]
        dz = xyz[ro + 2] - xyz[so + 2]
        den = dx * dx + dy * dy + dz * dz
        if den == 0.0:
            t = 0.5
        else:
            t = (
                (x - xyz[so]) * dx + (y - xyz[so + 1]) * dy + (z - xyz[so + 2]) * dz
            ) / den
            if t < 0.0:
                t = 0.0
            if t > 1.0:
                t = 1.0
        for _, ncomp, vals in float_data:
            sd = s * ncomp
            rd = r * ncomp
            for k in range(ncomp):
                vals[sd + k] = vals[sd + k] + t * (vals[rd + k] - vals[sd + k])

        xyz[so] = x
        xyz[so + 1] = y
        xyz[so + 2] = z
        sq = s * 10
        rq = r * 10
        for i in range(10):
            quads[sq + i] += quads[rq + i]

        for f in shared:
            face_alive[f] = False
            alive_faces -= 1
            faces_removed += 1
            o = f * 3
            for cid in (corners[o], corners[o + 1], corners[o + 2]):
                if cid != s and cid != r:
                    row = vfaces[cid]
                    k = bisect_left(row, f)
                    if k < len(row) and row[k] == f:
                        del row[k]
            row = vfaces[s]
            k = bisect_left(row, f)
            if k < len(row) and row[k] == f:
                del row[k]
        for f in vfaces[r]:
            if not face_alive[f]:
                continue  # was in shared
            o = f * 3
            for k in range(3):
                if corners[o + k] == r:
                    corners[o + k] = s
            insort(vfaces[s], f)
        vfaces[r] = []
        successor[r] = s
        version[s] += 1
        version[r] += 1
        max_error_applied = max(max_error_applied, err)

        for w2 in _vertex_ring(vfaces[s], corners, s):
            if pinned_l[s] and pinned_l[w2]:
                continue
            lo, hi = (s, w2) if s < w2 else (w2, s)
            e2 = _place(xyz, quads, pinned_l, placement, lo, hi)[3]
            if not math.isfinite(e2):
                continue
            heapq.heappush(heap, (e2, lo, hi, version[lo], version[hi]))

    # ---- emission (twin of decimate.cpp's phase 7) --------------------------
    used = np.zeros(n, dtype=bool)
    corners_out = np.asarray(corners, dtype=np.int64).reshape(-1, 3)
    alive_mask = np.asarray(face_alive, dtype=bool)
    if num_faces:
        used[corners_out[alive_mask].reshape(-1)] = True
    remap = np.full(n, -1, dtype=np.int64)
    remap[used] = np.arange(int(used.sum()), dtype=np.int64)
    num_used = int(used.sum())

    src_dtype = np.asarray(simp.points).dtype
    xyz_arr = np.asarray(xyz, dtype=np.float64).reshape(-1, 3)
    out_points = xyz_arr[used][:, : min(dim, 3)].astype(src_dtype)

    new_cells = []
    cell_maps = []
    parent_blocks = simp.cell_data.get("convert:parent_cell")
    in_cells = [len(np.asarray(cb.data)) for cb in mesh.cells]
    for bi, data in enumerate(block_cells):
        nc = len(data)
        lo = block_base[bi]
        blk_alive = alive_mask[lo : lo + nc]
        conn = remap[corners_out[lo : lo + nc][blk_alive]]
        new_cells.append(("triangle", conn.astype(np.int64)))
        nparents = in_cells[bi] if bi < len(in_cells) else nc
        cmap = np.full(nparents, -1, dtype=np.int64)
        parents = (
            np.asarray(parent_blocks[bi], dtype=np.int64)
            if parent_blocks is not None and bi < len(parent_blocks)
            else np.arange(nc, dtype=np.int64)
        )
        w = 0
        for c in range(nc):
            if not blk_alive[c]:
                continue
            p = int(parents[c])
            if 0 <= p < nparents and cmap[p] < 0:
                cmap[p] = w
            w += 1
        cell_maps.append(cmap)

    cell_data = {}
    for name in sorted(simp.cell_data):
        if name == "convert:parent_cell":
            continue
        blocks = []
        for bi, value in enumerate(simp.cell_data[name]):
            if value is None or bi >= len(block_cells):
                blocks.append(None if value is None else np.asarray(value).copy())
                continue
            value = np.asarray(value)
            nc = len(block_cells[bi])
            if len(value) != nc or nc == 0:
                blocks.append(value.copy())
                continue
            lo = block_base[bi]
            blocks.append(value[alive_mask[lo : lo + nc]])
        cell_data[name] = blocks

    point_data = {}
    for name in sorted(simp.point_data):
        arr = np.asarray(simp.point_data[name])
        if len(arr) != n:
            point_data[name] = arr.copy()
            continue
        if np.issubdtype(arr.dtype, np.floating):
            blended = None
            for fname, ncomp, vals in float_data:
                if fname == name:
                    blended = np.asarray(vals, dtype=np.float64).reshape(n, ncomp)
                    break
            sub = blended[used].astype(arr.dtype)
            point_data[name] = sub.reshape((num_used,) + arr.shape[1:])
        else:
            point_data[name] = arr[used]

    out = Mesh(
        out_points,
        new_cells,
        point_data=point_data,
        cell_data=cell_data,
        field_data={k: np.asarray(v).copy() for k, v in simp.field_data.items()},
    )

    point_map = np.empty(n, dtype=np.int64)
    for i in range(n):
        v = i
        while successor[v] >= 0:
            v = successor[v]
        point_map[i] = remap[v]

    report = {
        "faces_removed": int(faces_removed),
        "points_removed": int(n - num_used),
        "collapses_rejected": int(collapses_rejected),
        "max_error_applied": float(max_error_applied),
    }
    return out, point_map, cell_maps, report


# --------------------------------------------------------------------------- #
# set remapping (never reaches the C++ core)                                   #
# --------------------------------------------------------------------------- #
def _remap_sets(src, out, point_map, cell_maps):
    point_sets = getattr(src, "point_sets", None)
    if point_sets:
        pm = np.asarray(point_map, dtype=np.int64)
        new_ps = {}
        for name, idx in point_sets.items():
            mapped = pm[np.asarray(idx, dtype=np.int64)]
            # Several set members may collapse onto the same survivor.
            new_ps[name] = np.unique(mapped[mapped >= 0])
        out.point_sets = new_ps

    cell_sets = getattr(src, "cell_sets", None)
    if cell_sets and cell_maps is not None:
        new_cs = {}
        for name, blocks in cell_sets.items():
            remapped = []
            for b, idx in enumerate(blocks):
                if idx is None or b >= len(cell_maps):
                    remapped.append(None)
                    continue
                cm = np.asarray(cell_maps[b], dtype=np.int64)
                idx = np.asarray(idx, dtype=np.int64)
                idx = idx[(idx >= 0) & (idx < len(cm))]
                mapped = cm[idx]
                remapped.append(np.unique(mapped[mapped >= 0]))
            new_cs[name] = remapped
        out.cell_sets = new_cs


# --------------------------------------------------------------------------- #
# public API                                                                   #
# --------------------------------------------------------------------------- #
def decimate(
    mesh,
    ratio=None,
    target_faces=None,
    max_error=None,
    placement="optimal",
    preserve_boundary=True,
    preserve_features=True,
    feature_angle=30.0,
    frozen=None,
    return_report=False,
):
    """Decimate a surface mesh by quadric-error-metric edge collapse.

    Exactly one of ``ratio``, ``target_faces``, ``max_error`` must be given.
    The output is all-triangle: ``quad``/``polygon`` blocks are triangulated
    first, and the block structure stays 1:1 with the input.

    :param mesh: the surface mesh to decimate (never modified). Volume meshes
        raise -- run :func:`extract_surface` first.
    :param ratio: fraction of the (triangulated) faces to KEEP, in ``(0, 1]``.
    :param target_faces: absolute face count to stop at. A collapse removes one
        or two faces, so the result lands within one collapse of the target.
    :param max_error: collapse only while the cheapest candidate's quadric
        error is at most this (squared mesh units).
    :param placement: where the surviving vertex goes -- ``"optimal"`` (the
        quadric minimizer, midpoint when ill-conditioned), ``"midpoint"``, or
        ``"endpoint"`` (the endpoint with the lower error).
    :param preserve_boundary: pin boundary vertices (once-used-edge test), so
        the outline of an open patch is preserved exactly.
    :param preserve_features: pin vertices whose incident face normals differ
        by more than ``feature_angle``, so corners and creases survive.
    :param feature_angle: that angle, in degrees, between face normals.
    :param frozen: optional extra vertex ids to pin, as an index array or the
        name of one of ``mesh.point_sets``.
    :param return_report: also return the run summary dict
        (``faces_removed``, ``points_removed``, ``collapses_rejected``,
        ``max_error_applied``).
    :returns: the decimated mesh, or ``(mesh, report)`` if ``return_report``.
    :raises ValueError: when not exactly one stopping criterion is given, and
        on any block outside the surface scope (volume, higher-order, ragged,
        ``line``/``vertex``).

    ``point_data`` at a surviving vertex is blended between the collapse's two
    endpoints (clamped edge parameter; exact for midpoint/endpoint placement,
    an approximation for optimal) -- except **integer** arrays, which keep the
    survivor's own value, a blended material id being meaningless. Blended
    arrays keep their float dtype; the run accumulates in float64.
    """
    if isinstance(frozen, str):
        try:
            frozen = np.asarray(mesh.point_sets[frozen], dtype=np.int64)
        except KeyError:
            raise ValueError(
                f"meshio++: decimate: no point_set named {frozen!r} "
                f"(available: {sorted(mesh.point_sets)})"
            ) from None

    target_ratio = -1.0 if ratio is None else float(ratio)
    faces = -1 if target_faces is None else int(target_faces)
    error = -1.0 if max_error is None else float(max_error)

    out = None
    point_map = None
    cell_maps = None
    report = None
    try:
        from . import _core

        res = _core.decimate(
            mesh,
            target_ratio,
            faces,
            error,
            placement,
            bool(preserve_boundary),
            bool(preserve_features),
            float(feature_angle),
            None if frozen is None else np.asarray(frozen, dtype=np.int64),
        )
        out = res["mesh"]
        point_map = np.asarray(res["point_map"])
        cell_maps = [np.asarray(a) for a in res["cell_maps"]]
        report = {
            "faces_removed": res["faces_removed"],
            "points_removed": res["points_removed"],
            "collapses_rejected": res["collapses_rejected"],
            "max_error_applied": res["max_error_applied"],
        }
    except (ValueError, TypeError):
        # A genuine user error (no criterion, a volume mesh, a bad placement, a
        # frozen id out of bounds) must not fall through to the numpy path and
        # there either succeed silently or raise something less precise.
        raise
    except Exception:
        out = None
    used_cpp = out is not None
    if out is None:
        out, point_map, cell_maps, report = _decimate_py(
            mesh,
            target_ratio,
            faces,
            error,
            placement,
            bool(preserve_boundary),
            bool(preserve_features),
            float(feature_angle),
            frozen,
        )

    # The C++ core carries named regions across itself (and therefore
    # `point_sets`/`cell_sets`, which are views over them), so remapping again
    # here would apply the maps twice. Only the numpy fallback needs it.
    if not used_cpp:
        _remap_sets(mesh, out, point_map, cell_maps)
    return (out, report) if return_report else out
