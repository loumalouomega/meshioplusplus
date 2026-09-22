"""Point and cell normals of a surface mesh, optionally splitting vertices at
creases.

A *vertex* normal is a property of a smooth patch, not of a position: at the
edge of a cube the one position has three normals. With ``split_angle=None`` the
operation returns one normal per point, the angle- (or area-) weighted mean of
the incident faces, and smooths the crease over. With a split angle it does what
a renderer needs: the corners around a vertex are grouped into smooth fans and
every fan beyond the first gets its own copy of the point, so each point carries
exactly one normal.

The fans are found with a union-find over triangle corners. Two corners of one
vertex belong to the same fan when the triangles that own them are joined along
an edge through that vertex; an edge joins its two triangles only when they are a
proper manifold pair (exactly two users, walked in opposite directions) whose
dihedral angle is within the split angle. Boundary, non-manifold and
inconsistently wound edges therefore always cut, and the triangles fanned from
one polygon are always joined, so a non-planar polygon never splits along its own
diagonals. The union-find root is always the smallest corner index of its set and
groups are numbered by ascending root, so the partition depends only on the
soup and the angle. This module is the numpy twin of
``detail/surface_normals.cpp`` and ``operations/normals.cpp``; the glTF writer's
Python reference is built on :func:`vertex_normal_groups`.

The split layout is the one :func:`meshioplusplus.repair` uses for its bowtie
copies: originals keep their indices, copies are appended, ``point_data`` gathers
by row, and a copy joins its source's point regions. Cell numbering is untouched.

**It never reorients.** Two triangles that disagree about which side is out
cannot both be right, and averaging them gives a wrong normal that looks
plausible. ``report["quality"]["inconsistent_pairs"]`` reports the count; a split
always cuts at such an edge; the fix is ``repair(mesh, fix_orientation=True)``.

The arrays are called ``normals`` -- the name the PCD and XYZ formats already use
-- so ``compute_normals`` followed by a write to ``.pcd`` or ``.xyz`` emits the
normal columns with no further step. An existing array of that name is replaced.

Public API:

* :func:`compute_normals` -- the operation.
"""

from __future__ import annotations

import math

import numpy as np

from ._regions import Region, block_bases

_PREFIX = "meshio++: normals: "

#: Point / cell data written by :func:`compute_normals`. Twins of
#: ``normals.hpp``'s ``kNormals*Name`` constants -- keep the two in step.
NORMALS_NAME = "normals"
PARENT_POINT_NAME = "normals:parent_point"

_WEIGHTS = ("angle", "area")


# --------------------------------------------------------------------------- #
# the triangle soup                                                            #
# --------------------------------------------------------------------------- #
def _dot(u, v):
    """``vec3_dot``'s own left-to-right order, never ``np.sum``/``einsum``."""
    return u[:, 0] * v[:, 0] + u[:, 1] * v[:, 1] + u[:, 2] * v[:, 2]


def _cross(u, v):
    """``vec3_cross``'s own component order."""
    return np.column_stack(
        [
            u[:, 1] * v[:, 2] - u[:, 2] * v[:, 1],
            u[:, 2] * v[:, 0] - u[:, 0] * v[:, 2],
            u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0],
        ]
    )


def _check_surface(mesh, region):
    """Refuse what a normal is not defined on, naming the fix."""
    from ._mesh import topological_dimension

    for cb in mesh.cells:
        name = cb.type
        polygon = name.startswith("polygon")
        if name.startswith("polyhedron") or (
            not polygon and topological_dimension.get(name) == 3
        ):
            raise ValueError(
                f"{_PREFIX}cell block '{name}' is a volume; normals are defined on "
                "a surface (run extract_surface first)"
            )
        if (
            not polygon
            and topological_dimension.get(name) == 2
            and name not in ("triangle", "quad")
        ):
            raise ValueError(
                f"{_PREFIX}cell block '{name}' is a higher-order surface cell "
                "(run linearize first)"
            )
    if region:
        names = {r.name for r in (getattr(mesh, "regions", ()) or ())}
        if not any(
            r.name == region and r.kind == "cell"
            for r in (getattr(mesh, "regions", ()) or ())
        ):
            raise ValueError(
                f"{_PREFIX}no cell region named '{region}' "
                f"(available: {', '.join(sorted(names)) if names else 'none'})"
            )


def _region_mask(mesh, region):
    if not region:
        return None
    bases = block_bases(mesh.cells)
    total = int(bases[-1]) + len(mesh.cells[-1].data) if len(mesh.cells) else 0
    mask = np.zeros(total, dtype=bool)
    for r in getattr(mesh, "regions", ()) or ():
        if r.name == region and r.kind == "cell":
            e = np.asarray(r.entries, dtype=np.int64).ravel()
            mask[e[(e >= 0) & (e < total)]] = True
    return mask


def soup(mesh, region=""):
    """The triangle soup of a surface mesh, region filter included.

    :returns: ``(points, verts, corners, source)``: the ``(n, 3)`` float64
        points, the ``(t, 3)`` vertex ids and ``(t, 3, 3)`` corner coordinates
        of each triangle, and the global cell index each triangle came from.
    """
    from ._sdf import _soup as sdf_soup

    _check_surface(mesh, region)
    points, verts, corners, source = sdf_soup(mesh)
    mask = _region_mask(mesh, region)
    if mask is not None and len(verts):
        keep = mask[source]
        verts, corners, source = verts[keep], corners[keep], source[keep]
    return points, verts, corners, source


def _quality(verts, corners):
    """``soup_quality``: the four edge defect counts and the verdict."""
    degenerate = 0
    if len(verts):
        n = _cross(corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0])
        degenerate = int(np.count_nonzero(~(_dot(n, n) > 0.0)))
    counts = {}
    for tri in verts:
        for e in range(3):
            u, w = int(tri[e]), int(tri[(e + 1) % 3])
            key = (u, w) if u < w else (w, u)
            used, forward = counts.get(key, (0, 0))
            counts[key] = (used + 1, forward + (1 if u < w else 0))
    boundary = non_manifold = inconsistent = 0
    for used, forward in counts.values():
        if used == 1:
            boundary += 1
        elif used > 2:
            non_manifold += 1
        elif used == 2 and forward != 1:
            inconsistent += 1
    return {
        "boundary_edges": boundary,
        "non_manifold_edges": non_manifold,
        "inconsistent_pairs": inconsistent,
        "degenerate_triangles": degenerate,
        "watertight": boundary == 0
        and non_manifold == 0
        and inconsistent == 0
        and degenerate == 0,
    }


# --------------------------------------------------------------------------- #
# the corner grouping                                                          #
# --------------------------------------------------------------------------- #
def _find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def _unite(parent, a, b):
    """Union with the smaller root winning, as ``SnUnionFind::Unite`` does."""
    ra, rb = _find(parent, a), _find(parent, b)
    if ra == rb:
        return
    if ra < rb:
        parent[rb] = ra
    else:
        parent[ra] = rb


def face_normals(corners):
    """``soup_face_normals``: ``cross(b - a, c - a)`` per triangle."""
    if not len(corners):
        return np.zeros((0, 3))
    return _cross(corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0])


def _corner_angles(corners):
    """``corner_angle`` at each of the three corners, per triangle.

    ``math.acos`` per element rather than ``np.arccos``: numpy's vectorised
    ``arccos`` is free to use a SIMD implementation that differs from libm in
    the last bit, and the C++ core calls ``std::acos``.
    """
    ntri = len(corners)
    out = np.zeros((ntri, 3))
    for i in range(3):
        a = corners[:, i]
        u = corners[:, (i + 1) % 3] - a
        v = corners[:, (i + 2) % 3] - a
        nu = np.sqrt(_dot(u, u))
        nv = np.sqrt(_dot(v, v))
        ok = (nu > 0.0) & (nv > 0.0)
        c = np.zeros(ntri)
        c[ok] = _dot(u, v)[ok] / (nu[ok] * nv[ok])
        c = np.clip(c, -1.0, 1.0)
        ang = np.array([math.acos(x) for x in c.tolist()]) if ntri else np.zeros(0)
        out[:, i] = np.where(ok, ang, 0.0)
    return out


def vertex_normal_groups(
    verts, corners, source, npts, weight="angle", split_angle=-1.0
):
    """Group the corners of a soup into smooth fans and compute a unit normal each.

    The numpy twin of ``detail::vertex_normal_groups``.

    :param split_angle: degrees; an edge joins its two triangles only when the
        angle between their face normals does not exceed it. Negative disables
        splitting (one group per touched point); 180 and above keep every proper
        manifold pair joined.
    :returns: a dict with ``corner_group`` ``(3t,)``, ``group_point``,
        ``group_root``, ``group_normal`` ``(g, 3)`` (zeros when undefined) and
        ``num_degenerate``.
    """
    ntri = len(verts)
    ncorner = ntri * 3
    face = face_normals(corners)
    length = np.sqrt(_dot(face, face)) if ntri else np.zeros(0)
    ok = length > 0.0
    unit = np.zeros((ntri, 3))
    if ntri:
        unit[ok] = face[ok] * (1.0 / length[ok])[:, None]
    num_degenerate = int(np.count_nonzero(~ok))

    flat_pts = verts.reshape(-1)
    parent = list(range(ncorner))
    if split_angle < 0.0:
        if ncorner:
            uniq, first = np.unique(flat_pts, return_index=True)
            pos = np.searchsorted(uniq, flat_pts)
            root_flat = first[pos]
            for c in np.flatnonzero(root_flat != np.arange(ncorner)).tolist():
                parent[c] = int(root_flat[c])
    else:
        always = split_angle >= 180.0
        cos_threshold = math.cos(split_angle * (math.pi / 180.0))
        tri_idx = np.arange(ntri, dtype=np.int64)
        recs = []
        for i in range(3):
            j = (i + 1) % 3
            u, w = verts[:, i], verts[:, j]
            keep = u != w
            fwd = u < w
            lo = np.where(fwd, u, w)[keep]
            hi = np.where(fwd, w, u)[keep]
            ci = (tri_idx * 3 + i)[keep]
            cj = (tri_idx * 3 + j)[keep]
            recs.append(
                (
                    lo,
                    hi,
                    tri_idx[keep],
                    np.where(fwd[keep], ci, cj),
                    np.where(fwd[keep], cj, ci),
                    fwd[keep],
                )
            )
        if recs and sum(len(r[0]) for r in recs):
            lo, hi, tri, clo, chi, fwd = (np.concatenate(x) for x in zip(*recs))
            order = np.lexsort((clo, tri, hi, lo))
            lo, hi, tri, clo, chi, fwd = (
                lo[order],
                hi[order],
                tri[order],
                clo[order],
                chi[order],
                fwd[order],
            )
            new_edge = np.ones(len(lo), dtype=bool)
            new_edge[1:] = (lo[1:] != lo[:-1]) | (hi[1:] != hi[:-1])
            starts = np.flatnonzero(new_edge).tolist() + [len(lo)]
            unit_l = unit.tolist()
            length_l = length.tolist()
            source_l = source.tolist()
            tri_l, clo_l, chi_l, fwd_l = (
                tri.tolist(),
                clo.tolist(),
                chi.tolist(),
                fwd.tolist(),
            )
            for s in range(len(starts) - 1):
                b, e = starts[s], starts[s + 1]
                proper = e - b == 2
                for x in range(b, e):
                    for y in range(x + 1, e):
                        tx, ty = tri_l[x], tri_l[y]
                        if tx == ty:
                            continue
                        join = source_l[tx] == source_l[ty]
                        if (
                            not join
                            and proper
                            and fwd_l[x] != fwd_l[y]
                            and length_l[tx] > 0.0
                            and length_l[ty] > 0.0
                        ):
                            ux, uy = unit_l[tx], unit_l[ty]
                            join = always or (
                                ux[0] * uy[0] + ux[1] * uy[1] + ux[2] * uy[2]
                                >= cos_threshold
                            )
                        if join:
                            _unite(parent, clo_l[x], clo_l[y])
                            _unite(parent, chi_l[x], chi_l[y])

    root = [_find(parent, c) for c in range(ncorner)]
    gid = {}
    group_root = []
    group_point = []
    corner_group = np.zeros(ncorner, dtype=np.int64)
    for c in range(ncorner):
        if root[c] == c:
            gid[c] = len(group_root)
            group_root.append(c)
            group_point.append(int(flat_pts[c]))
        corner_group[c] = gid[root[c]]

    ngroups = len(group_root)
    acc = np.zeros((ngroups, 3))
    if ntri:
        if weight == "angle":
            w = _corner_angles(corners)
        else:
            w = np.repeat(length[:, None], 3, axis=1)
        idx = corner_group.reshape(ntri, 3)[ok].ravel()
        val = (unit[:, None, :] * w[:, :, None])[ok].reshape(-1, 3)
        np.add.at(acc, idx, val)
    norm = np.sqrt(_dot(acc, acc)) if ngroups else np.zeros(0)
    good = norm > 0.0
    normal = np.zeros((ngroups, 3))
    if ngroups:
        normal[good] = acc[good] * (1.0 / norm[good])[:, None]
    return {
        "corner_group": corner_group,
        "group_point": np.asarray(group_point, dtype=np.int64),
        "group_root": np.asarray(group_root, dtype=np.int64),
        "group_normal": normal,
        "num_degenerate": num_degenerate,
    }


# --------------------------------------------------------------------------- #
# the operation                                                                #
# --------------------------------------------------------------------------- #
def _validate(weight, split_angle):
    if weight not in _WEIGHTS:
        raise ValueError(
            f"{_PREFIX}unknown weight '{weight}' (expected 'angle' or 'area')"
        )
    if split_angle is not None and not (0.0 <= float(split_angle) <= 180.0):
        raise ValueError(f"{_PREFIX}the split angle must lie in [0, 180] degrees")


def _compute_normals_py(
    mesh,
    point_normals,
    cell_normals,
    weight,
    split_angle,
    record_parent_ids,
    region,
):
    """The numpy twin of ``compute_normals``."""
    points, verts, corners, source = soup(mesh, region)
    n = len(points)
    ntri = len(verts)
    nan = float("nan")
    quality = _quality(verts, corners)

    g = vertex_normal_groups(
        verts,
        corners,
        source,
        n,
        weight,
        float(split_angle) if split_angle is not None else -1.0,
    )
    ngroups = len(g["group_point"])
    gpoint = g["group_point"].tolist()
    gnormal = g["group_normal"]
    defined = (gnormal != 0.0).any(axis=1).tolist() if ngroups else []

    # A group with no direction folds into its point's primary group.
    primary = [-1] * n
    for k in range(ngroups):
        if primary[gpoint[k]] < 0 and defined[k]:
            primary[gpoint[k]] = k
    for k in range(ngroups):
        if primary[gpoint[k]] < 0:
            primary[gpoint[k]] = k

    group_out = [0] * ngroups
    parent_of_new = []
    has_copy = set()
    for k in range(ngroups):
        p = gpoint[k]
        if k == primary[p] or not defined[k]:
            group_out[k] = p
        else:
            group_out[k] = n + len(parent_of_new)
            parent_of_new.append(p)
            has_copy.add(p)
    n_added = len(parent_of_new)
    n_out = n + n_added

    pn = np.full((n_out, 3), nan)
    for k in range(ngroups):
        p = gpoint[k]
        if k != primary[p] and group_out[k] == p:
            continue
        if not defined[k]:
            continue
        pn[group_out[k]] = gnormal[k]
    num_isolated = sum(1 for p in range(n) if primary[p] < 0)
    num_undefined = sum(
        1 for p in range(n) if primary[p] >= 0 and not defined[primary[p]]
    )

    bases = block_bases(mesh.cells)
    total = int(bases[-1]) + len(mesh.cells[-1].data) if len(mesh.cells) else 0
    first_tri = np.full(total, -1, dtype=np.int64)
    if ntri:
        _uniq, first = np.unique(source, return_index=True)
        first_tri[_uniq] = first
    cell_sum = np.zeros((total, 3))
    if ntri:
        np.add.at(cell_sum, source, face_normals(corners))

    out = mesh.copy()
    if n_added:
        parent_arr = np.asarray(parent_of_new, dtype=np.int64)
        cg = g["corner_group"]
        go = np.asarray(group_out, dtype=np.int64)

        def corner_point(t, i):
            return int(go[cg[t * 3 + i]])

        def rewrite(base, c, ids):
            t0 = int(first_tri[base + c])
            if t0 < 0 or len(ids) < 3:
                return ids
            nv = len(ids)
            ids = list(ids)
            ids[0] = corner_point(t0, 0)
            for k in range(1, nv - 1):
                ids[k] = corner_point(t0 + k - 1, 1)
            ids[nv - 1] = corner_point(t0 + nv - 3, 2)
            return ids

        out.points = np.concatenate(
            [np.asarray(mesh.points), np.asarray(mesh.points)[parent_arr]]
        )
        for b, cb in enumerate(out.cells):
            base = int(bases[b])
            if isinstance(cb.data, list):
                cb.data = [rewrite(base, c, row) for c, row in enumerate(cb.data)]
            else:
                data = np.array(cb.data, dtype=np.int64)
                for c in range(len(data)):
                    data[c] = rewrite(base, c, data[c].tolist())
                cb.data = data
        for name in list(out.point_data):
            a = np.asarray(mesh.point_data[name])
            out.point_data[name] = (
                np.concatenate([a, a[parent_arr]]) if len(a) == n and n else a.copy()
            )
        new_regions = []
        for r in out.regions:
            if r.kind != "point":
                new_regions.append(r)
                continue
            entries = np.asarray(r.entries, dtype=np.int64)
            members = set(entries.tolist())
            extra = [n + k for k, p in enumerate(parent_of_new) if p in members]
            new_regions.append(
                Region(
                    r.name,
                    r.kind,
                    np.concatenate([entries, np.asarray(extra, dtype=np.int64)]),
                    r.dim,
                    r.tag,
                )
            )
        out.regions = new_regions

    if point_normals:
        out.point_data[NORMALS_NAME] = pn
    if cell_normals:
        blocks = []
        for b, cb in enumerate(mesh.cells):
            base = int(bases[b])
            arr = np.full((len(cb.data), 3), nan)
            for c in range(len(cb.data)):
                gc = base + c
                s = cell_sum[gc : gc + 1]
                ln = float(np.sqrt(_dot(s, s))[0])
                if first_tri[gc] >= 0 and ln > 0.0:
                    arr[c] = s[0] * (1.0 / ln)
            blocks.append(arr)
        out.cell_data[NORMALS_NAME] = blocks
    if record_parent_ids:
        out.point_data[PARENT_POINT_NAME] = np.concatenate(
            [np.arange(n, dtype=np.int64), np.asarray(parent_of_new, dtype=np.int64)]
        )

    report = {
        "num_isolated": num_isolated,
        "num_undefined": num_undefined,
        "num_degenerate": g["num_degenerate"],
        "num_split_points": len(has_copy),
        "num_added_points": n_added,
        "quality": quality,
    }
    return out, report


def compute_normals(
    mesh,
    point_normals=True,
    cell_normals=False,
    weight="angle",
    split_angle=None,
    region="",
    record_parent_ids=False,
    return_report=False,
):
    """Point and/or cell normals of a surface mesh.

    Triangles come from the same fan :func:`meshioplusplus.convert_cells`'s
    ``simplexify`` uses, so a quad mesh's normals are those of its canonical
    triangulation; a polygon's cell normal is the sum of its fan (Newell's
    normal), which does not depend on which corner the fan starts from. Lines and
    vertices are ignored and get NaN. A volume or polyhedron block is refused by
    name pointing at ``extract_surface``, a higher-order one pointing at
    ``linearize``.

    :param mesh: a surface mesh (never modified).
    :param point_normals: attach the unit point normals as ``point_data["normals"]``,
        shape ``(n, 3)``.
    :param cell_normals: attach the unit vector area of each cell as
        ``cell_data["normals"]``.
    :param weight: ``"angle"`` (default) weights each incident face by the angle
        it subtends at the vertex; ``"area"`` by its area.
    :param split_angle: ``None`` (default) gives one smooth normal per point. A
        number of degrees in ``[0, 180]`` duplicates points wherever the surface
        creases by more than that, appending the copies after the original
        points.
    :param region: restrict to this named ``Cell`` region; ``""`` takes every
        surface cell.
    :param record_parent_ids: also attach ``point_data["normals:parent_point"]``,
        for each output point the input point it came from.
    :param return_report: also return the counters and the input's surface
        ``quality``.
    :returns: the mesh, or ``(mesh, report)`` when ``return_report`` is set. The
        report has ``num_isolated`` and ``num_undefined`` (points whose normal is
        NaN), ``num_degenerate`` (triangles with no area), ``num_split_points``
        and ``num_added_points``, and the input's ``quality``.
    :raises ValueError: on a non-surface input (naming the fix), an unknown
        region name, or a split angle outside ``[0, 180]``.

    It never reorients: check ``report["quality"]["inconsistent_pairs"]`` and run
    ``repair(mesh, fix_orientation=True)`` first if it is non-zero.
    """
    _validate(weight, split_angle)

    out = None
    report = None
    try:
        from . import _core

        res = _core.compute_normals(
            mesh,
            bool(point_normals),
            bool(cell_normals),
            weight,
            split_angle is not None,
            float(split_angle) if split_angle is not None else 30.0,
            bool(record_parent_ids),
            str(region),
        )
        out = res["mesh"]
        report = {
            k: res[k]
            for k in (
                "num_isolated",
                "num_undefined",
                "num_degenerate",
                "num_split_points",
                "num_added_points",
                "quality",
            )
        }
    except (ValueError, TypeError):
        # A genuine user error must not fall through to the numpy path.
        raise
    except Exception:
        out = None

    if out is None:
        out, report = _compute_normals_py(
            mesh,
            bool(point_normals),
            bool(cell_normals),
            weight,
            split_angle,
            bool(record_parent_ids),
            str(region),
        )

    if report["quality"]["inconsistent_pairs"] and split_angle is None:
        from ._common import warn

        warn(
            f"normals: {report['quality']['inconsistent_pairs']} edge pair(s) wind the "
            "same way, so the normals there average faces that disagree about which "
            "side is out; run repair(mesh, fix_orientation=True) first"
        )

    return (out, report) if return_report else out


__all__ = ["compute_normals"]
