"""Analytic surfaces for the segment / fit_primitives / step tests.

Quads on a unit-cube skin and on surfaces of revolution, built from indices
so every point lies on its surface (up to rounding). The C++ tests build the
same shapes in ``tests/cpp/analytic_surfaces.hpp``.
"""

import math

import numpy as np

import meshioplusplus as mio


def box_skin(n=3):
    """The closed, outward-wound skin of the unit cube, ``n x n`` quads per face."""
    ids, pts, quads = {}, [], []

    def pid(x):
        key = tuple(int(v) for v in x)
        if key not in ids:
            ids[key] = len(pts)
            pts.append([v / n for v in key])
        return ids[key]

    for axis in range(3):
        for side in (0, 1):
            u, v = (axis + 1) % 3, (axis + 2) % 3
            for a in range(n):
                for b in range(n):
                    ring = []
                    for cu, cv in ((a, b), (a + 1, b), (a + 1, b + 1), (a, b + 1)):
                        x = [0, 0, 0]
                        x[axis], x[u], x[v] = side * n, cu, cv
                        ring.append(pid(x))
                    if side == 0:
                        ring = [ring[0], ring[3], ring[2], ring[1]]
                    quads.append(ring)
    return mio.Mesh(np.array(pts, dtype=float), [("quad", np.array(quads))])


def revolve(profile, n=32):
    """Quads on the surface of revolution about z through the ``(r, z)`` profile."""
    pts = []
    for r, z in profile:
        for k in range(n):
            a = 2.0 * np.pi * k / n
            pts.append([r * np.cos(a), r * np.sin(a), z])
    quads = []
    for i in range(len(profile) - 1):
        for k in range(n):
            k1 = (k + 1) % n
            quads.append([i * n + k, i * n + k1, (i + 1) * n + k1, (i + 1) * n + k])
    return mio.Mesh(np.array(pts, dtype=float), [("quad", np.array(quads))])


def cylinder(r=1.0, h=2.0, m=6, n=32):
    return revolve([(r, h * i / m) for i in range(m + 1)], n)


def cone(half_angle=20.0, m=6, n=32):
    t = np.tan(np.radians(half_angle))
    zs = [1.0 + 2.0 * i / m for i in range(m + 1)]
    return revolve([(t * z, z) for z in zs], n)


def sphere_band(r=1.5, m=12, n=32):
    ths = [0.3 + (np.pi - 0.6) * i / m for i in range(m + 1)]
    return revolve([(r * np.sin(t), -r * np.cos(t)) for t in ths], n)


def torus(major=3.0, minor=1.0, n=32, m=24):
    pts, quads = [], []
    for i in range(m):
        b = 2.0 * np.pi * i / m
        for k in range(n):
            a = 2.0 * np.pi * k / n
            rr = major + minor * np.cos(b)
            pts.append([rr * np.cos(a), rr * np.sin(a), minor * np.sin(b)])
    for i in range(m):
        for k in range(n):
            i1, k1 = (i + 1) % m, (k + 1) % n
            quads.append([i * n + k, i * n + k1, i1 * n + k1, i1 * n + k])
    return mio.Mesh(np.array(pts, dtype=float), [("quad", np.array(quads))])


def jitter(mesh, amplitude, seed=7):
    """The mesh with every coordinate moved by up to ``amplitude`` (seeded)."""
    rng = np.random.default_rng(seed)
    pts = mesh.points + amplitude * rng.uniform(-1.0, 1.0, mesh.points.shape)
    return mio.Mesh(pts, [(b.type, b.data) for b in mesh.cells])


def turned_part(n=64, noise=0.0, seed=7):
    """A closed turned part: a flat base, a cylinder wall, a toroidal fillet
    tangent to the wall and to a flat shoulder, a cone and a spherical dome --
    every primitive type, with creases between some faces and only curvature
    between others. Triangle fans close the two poles."""
    prof = [(2.0 * i / 4, 0.0) for i in range(1, 5)]
    prof += [(2.0, 1.5 * i / 8) for i in range(1, 9)]
    prof += [
        (1.5 + 0.5 * np.cos(a), 1.5 + 0.5 * np.sin(a))
        for a in np.linspace(0.0, np.pi / 2, 9)[1:]
    ]
    prof += [(1.5 - 0.5 * i / 3, 2.0) for i in range(1, 4)]
    prof += [(1.0 - 0.4 * i / 6, 2.0 + 1.0 * i / 6) for i in range(1, 7)]
    zc, big_r = 1.625, 1.5
    a0 = np.arcsin(0.6 / big_r)
    prof += [
        (big_r * np.sin(a), zc + big_r * np.cos(a))
        for a in np.linspace(a0, 0.0, 7)[1:-1]
    ]
    pts = [[0.0, 0.0, 0.0]]
    for r, z in prof:
        for k in range(n):
            t = 2.0 * np.pi * k / n
            pts.append([r * np.cos(t), r * np.sin(t), z])
    top = len(pts)
    pts.append([0.0, 0.0, zc + big_r])

    def ring(i, k):
        return 1 + i * n + (k % n)

    last = len(prof) - 1
    tris = [[0, ring(0, k + 1), ring(0, k)] for k in range(n)]
    tris += [[top, ring(last, k), ring(last, k + 1)] for k in range(n)]
    quads = [
        [ring(i, k), ring(i, k + 1), ring(i + 1, k + 1), ring(i + 1, k)]
        for i in range(last)
        for k in range(n)
    ]
    pts = np.array(pts)
    if noise:
        pts = pts + noise * np.random.default_rng(seed).uniform(-1.0, 1.0, pts.shape)
    return mio.Mesh(pts, [("triangle", np.array(tris)), ("quad", np.array(quads))])


def _closed_revolve(prof, z_bottom, z_top, n):
    """A closed surface of revolution about z, triangle fans at both ends."""
    pts = [[0.0, 0.0, z_bottom]]
    for r, z in prof:
        for k in range(n):
            t = 2.0 * np.pi * k / n
            pts.append([r * np.cos(t), r * np.sin(t), z])
    top = len(pts)
    pts.append([0.0, 0.0, z_top])

    def ring(i, k):
        return 1 + i * n + (k % n)

    last = len(prof) - 1
    tris = [[0, ring(0, k + 1), ring(0, k)] for k in range(n)]
    tris += [[top, ring(last, k), ring(last, k + 1)] for k in range(n)]
    for i in range(last):
        for k in range(n):
            tris.append([ring(i, k), ring(i, k + 1), ring(i + 1, k + 1)])
            tris.append([ring(i, k), ring(i + 1, k + 1), ring(i + 1, k)])
    return mio.Mesh(np.array(pts), [("triangle", np.array(tris))])


def wavy_capped(n=96, m=64):
    """A capped band of revolution whose radius waves along z, ``r = 1 + 0.15
    sin(2 pi z)`` for z in [0, 2]: no primitive explains the wall, and, closed
    around the axis, neither its plane nor an unrolled cylinder parametrises
    it (``analytic::wavy_capped``)."""
    prof = [(i / 3.0, 0.0) for i in range(1, 4)]
    prof += [
        (1.0 + 0.15 * np.sin(2.0 * np.pi * 2.0 * i / m), 2.0 * i / m)
        for i in range(1, m + 1)
    ]
    prof += [(i / 3.0, 2.0) for i in (2, 1)]
    return _closed_revolve(prof, 0.0, 2.0, n)


class _Skin:
    """Triangles wound along an expected outward direction."""

    def __init__(self):
        self.pts, self.tris = [], []

    def add(self, p):
        self.pts.append([float(v) for v in p])
        return len(self.pts) - 1

    def tri(self, a, b, c, out):
        p, q, r = (np.asarray(self.pts[i]) for i in (a, b, c))
        n = np.cross(q - p, r - p)
        self.tris.append([a, b, c] if np.dot(n, out) >= 0.0 else [a, c, b])

    def quad(self, a, b, c, d, out):
        self.tri(a, b, c, out)
        self.tri(a, c, d, out)

    def mesh(self):
        return mio.Mesh(np.array(self.pts), [("triangle", np.array(self.tris))])


def inclined_boss(alpha=20.0, r=0.6, rf=0.2, h=1.5, plate=1.6, t=0.4, n=96, k=8):
    """A cylinder boss leaning ``alpha`` degrees on a round plate, its foot
    rounded by a concave rolling-ball fillet of radius ``rf``: a canal surface
    along an ellipse that no primitive explains (``analytic::inclined_boss``)."""
    al = np.radians(alpha)
    ax = np.array([np.sin(al), 0.0, np.cos(al)])
    e1 = np.array([np.cos(al), 0.0, -np.sin(al)])
    b = _Skin()
    fillet = [[] for _ in range(k + 1)]
    centre, unit, foot = [], [], []
    for i in range(n):
        th = 2.0 * np.pi * i / n
        u = np.cos(th) * e1 + np.array([0.0, np.sin(th), 0.0])
        hh = (rf + (r + rf) * np.sin(al) * np.cos(th)) / np.cos(al)
        c = hh * ax + (r + rf) * u
        centre.append(c)
        unit.append(u)
        foot.append(hh)
        w = np.arccos(np.clip(u[2], -1.0, 1.0))
        for j in range(k + 1):
            s = j / k
            f0, f1 = np.sin((1.0 - s) * w) / np.sin(w), np.sin(s * w) / np.sin(w)
            d = -f1 * u - f0 * np.array([0.0, 0.0, 1.0])
            fillet[j].append(b.add(c + rf * d))

    def ring_quads(ra, rb, out):
        for i in range(n):
            j = (i + 1) % n
            b.quad(ra[i], ra[j], rb[j], rb[i], out(i))

    for j in range(k):
        ring_quads(
            fillet[j],
            fillet[j + 1],
            lambda i, j=j: centre[i] - np.asarray(b.pts[fillet[j][i]]),
        )
    prev = fillet[k]
    for j in range(1, 9):
        row = [
            b.add((foot[i] + (h - foot[i]) * j / 8) * ax + r * unit[i])
            for i in range(n)
        ]
        ring_quads(prev, row, lambda i: unit[i])
        prev = row
    top = b.add(h * ax)
    for i in range(n):
        b.tri(top, prev[i], prev[(i + 1) % n], ax)
    inner = fillet[0]
    psi = [np.arctan2(b.pts[v][1], b.pts[v][0]) for v in inner]
    rho0 = [np.hypot(b.pts[v][0], b.pts[v][1]) for v in inner]
    for j in range(1, 7):
        row = []
        for i in range(n):
            rr = rho0[i] + (plate - rho0[i]) * j / 6
            row.append(b.add([rr * np.cos(psi[i]), rr * np.sin(psi[i]), 0.0]))
        ring_quads(inner, row, lambda i: np.array([0.0, 0.0, 1.0]))
        inner = row
    lower = [b.add([plate * np.cos(p), plate * np.sin(p), -t]) for p in psi]
    ring_quads(inner, lower, lambda i: np.array([np.cos(psi[i]), np.sin(psi[i]), 0.0]))
    bottom = b.add([0.0, 0.0, -t])
    for i in range(n):
        b.tri(bottom, lower[i], lower[(i + 1) % n], np.array([0.0, 0.0, -1.0]))
    return b.mesh()


def _helical_rod(w, rad, pitch, turns, n, depth=1.0, taper=0.0, pitch_rate=0.0):
    """A rod about z whose surface is ``r = f(z - K theta)``, ``K = pitch /
    (2 pi)``, one period of ``f`` given as nodes ``(w, rad)``; the grid
    follows the helices through every node, and each end is closed by a fan
    to a point on the axis ``depth`` pitches inside the rod.

    ``taper`` and ``pitch_rate`` generalise the screw: a node at the unwrapped
    angle ``T`` sits at ``z = w + g(T)``, ``g(T) = K T + pitch_rate T^2 / (8
    pi^2)`` (the pitch grows by ``pitch_rate`` per turn), and ``r = f + taper
    (w + g(T))`` (``analytic::helical_rod``)."""
    k = pitch / (2.0 * np.pi)
    per = len(w)
    rows = per * turns
    b = _Skin()
    ids = []
    for i in range(n):
        th = 2.0 * np.pi * i / n
        col = []
        for j in range(rows + 1):
            if taper == 0.0 and pitch_rate == 0.0:
                wj = w[j % per] + pitch * (j // per)
                rj = rad[j % per]
                col.append(b.add([rj * np.cos(th), rj * np.sin(th), wj + k * th]))
                continue
            big = th + 2.0 * np.pi * (j // per)
            g = k * big + pitch_rate * big * big / (8.0 * np.pi * np.pi)
            rj = rad[j % per] + taper * (w[j % per] + g)
            col.append(b.add([rj * np.cos(th), rj * np.sin(th), w[j % per] + g]))
        ids.append(col)

    def node(i, j):
        if i < n:
            return ids[i][j]
        return ids[0][j + per] if j + per <= rows else -1

    for i in range(n):
        for j in range(rows):
            a, bb, c, d = node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)
            if bb < 0 or c < 0:
                continue
            pa, pd = np.asarray(b.pts[a]), np.asarray(b.pts[d])
            ra, rd = np.hypot(pa[0], pa[1]), np.hypot(pd[0], pd[1])
            radial = np.array([pa[0], pa[1], 0.0]) / ra
            radial[2] = -1.0 if rd > ra + 1e-12 else (1.0 if rd < ra - 1e-12 else 0.0)
            b.quad(a, bb, c, d, radial)

    def fan(loop, reverse):
        zs = [b.pts[v][2] for v in loop]
        c = b.add(
            [0.0, 0.0, max(zs) + depth * pitch if reverse else min(zs) - depth * pitch]
        )
        for i in range(len(loop)):
            p, q = loop[i], loop[(i + 1) % len(loop)]
            b.tris.append([c, q, p] if reverse else [c, p, q])

    low = [node(i, 0) for i in range(n)] + [node(0, j) for j in range(per, 0, -1)]
    high = [node(i, rows) for i in range(n)]
    high += [node(n - 1, j) for j in range(rows - 1, rows - per - 1, -1)]
    fan(low, True)
    fan(high, False)
    return b.mesh()


def threaded_rod(pitch=0.5, r0=0.8, r1=1.0, turns=4, n=64):
    """A rod with a trapezoidal thread: ``r = f(z - K theta)``, root and crest
    helical strips of cylinders, both flanks ruled helicoids, closed by two
    shallow fans a pitch inside the rod (``analytic::threaded_rod``)."""
    seg = [0.26 * pitch, 0.24 * pitch, 0.26 * pitch, 0.24 * pitch]
    sub = [2, 3, 2, 3]
    w, rad, at = [], [], 0.0
    for s in range(4):
        for i in range(sub[s]):
            tt = i / sub[s]
            w.append(at + seg[s] * tt)
            rad.append(
                [r0, r0 + (r1 - r0) * tt, r1, r1 - (r1 - r0) * tt][s],
            )
        at += seg[s]
    return _helical_rod(w, rad, pitch, turns, n)


def rounded_thread_rod(pitch=0.5, r0=0.8, r1=1.0, rr=0.2, turns=4, n=64):
    """A rod with a thread like ``threaded_rod``'s but no cylinder anywhere:
    the root and the crest are circular arcs of radius ``rr`` across the
    trapezoid's flats, the crest bulging out and the root in, meeting the
    straight flanks at a crease (``analytic::rounded_thread_rod``)."""
    seg = [0.26 * pitch, 0.24 * pitch, 0.26 * pitch, 0.24 * pitch]
    sub = [4, 3, 4, 3]
    half = 0.13 * pitch
    d = np.sqrt(rr * rr - half * half)
    w, rad, at = [], [], 0.0
    for s in range(4):
        for i in range(sub[s]):
            tt = i / sub[s]
            w.append(at + seg[s] * tt)
            x = seg[s] * tt - half  # from the flat's middle
            if s == 0:
                rad.append(r0 + d - np.sqrt(rr * rr - x * x))
            elif s == 1:
                rad.append(r0 + (r1 - r0) * tt)
            elif s == 2:
                rad.append(r1 - d + np.sqrt(rr * rr - x * x))
            else:
                rad.append(r1 - (r1 - r0) * tt)
        at += seg[s]
    return _helical_rod(w, rad, pitch, turns, n)


def _trapezoid_profile(pitch, r0, r1, rounded=None, split_root=False):
    """One period of ``threaded_rod``'s profile as ``(w, rad)``; ``rounded``
    (the arc radius) rounds the root and crest as ``rounded_thread_rod``
    does, and ``split_root`` runs the period from the middle of one root to
    the middle of the next."""
    if split_root:
        seg = [0.13 * pitch, 0.24 * pitch, 0.26 * pitch, 0.24 * pitch, 0.13 * pitch]
        sub = [1, 3, 2, 3, 1]
        kinds = [0, 1, 2, 3, 0]
    else:
        seg = [0.26 * pitch, 0.24 * pitch, 0.26 * pitch, 0.24 * pitch]
        sub = [4, 3, 4, 3] if rounded else [2, 3, 2, 3]
        kinds = [0, 1, 2, 3]
    half = 0.13 * pitch
    d = np.sqrt(rounded * rounded - half * half) if rounded else 0.0
    w, rad, at = [], [], 0.0
    for s, kind in enumerate(kinds):
        for i in range(sub[s]):
            tt = i / sub[s]
            w.append(at + seg[s] * tt)
            x = seg[s] * tt - half
            if kind == 0:
                rad.append(
                    r0 + d - np.sqrt(rounded * rounded - x * x) if rounded else r0
                )
            elif kind == 1:
                rad.append(r0 + (r1 - r0) * tt)
            elif kind == 2:
                rad.append(
                    r1 - d + np.sqrt(rounded * rounded - x * x) if rounded else r1
                )
            else:
                rad.append(r1 - (r1 - r0) * tt)
        at += seg[s]
    return w, rad


def tapered_thread_rod(pitch=0.5, r0=0.8, r1=1.0, taper=0.0625, turns=4, n=64):
    """``threaded_rod``'s thread cut on a cone whose radius grows by ``taper``
    per unit length: the root and crest are helical strips of cones, the
    flanks ruled helicoids (``analytic::tapered_thread_rod``)."""
    w, rad = _trapezoid_profile(pitch, r0, r1)
    return _helical_rod(w, rad, pitch, turns, n, taper=taper)


def tapered_rounded_thread_rod(
    pitch=0.5, r0=0.8, r1=1.0, rr=0.2, taper=0.0625, turns=4, n=64
):
    """``rounded_thread_rod``'s thread cut on a cone, with no cylinder or cone
    anywhere (``analytic::tapered_rounded_thread_rod``)."""
    w, rad = _trapezoid_profile(pitch, r0, r1, rounded=rr)
    return _helical_rod(w, rad, pitch, turns, n, taper=taper)


def variable_pitch_thread_rod(
    pitch=0.5, r0=0.8, r1=1.0, pitch_rate=0.05, turns=4, n=64
):
    """``threaded_rod``'s thread with a pitch growing by ``pitch_rate`` per
    turn, each period running from the middle of one root to the middle of the
    next so the root takes up the growing gap
    (``analytic::variable_pitch_thread_rod``)."""
    w, rad = _trapezoid_profile(pitch, r0, r1, split_root=True)
    return _helical_rod(w, rad, pitch, turns, n, pitch_rate=pitch_rate)


def section_blend_block(
    weight=0.0,
    length=2.0,
    width=1.0,
    height=1.0,
    d1=0.3,
    d2a=0.2,
    d2b=0.4,
    n=40,
    k=12,
    m=6,
):
    """The block ``[0, length] x [0, width] x [0, height]`` whose edge at
    ``y = width``, ``z = height`` is cut by a blend from the top, ``d1`` in from
    the edge, to the side, ``d2(x)`` down from it, growing linearly from
    ``d2a`` to ``d2b``. Each cross-section is the rational quadratic from the
    side foot through the corner to the top foot with weight ``weight``: 0 a
    chamfer of varying width, more a conic blend tangent to both faces. The
    triangles carry their faces as ``segment:patch`` (0 bottom, 1 side, 2
    top, 3 back, 4 and 5 the ends, 6 the blend)
    (``analytic::section_blend_block``)."""
    W, H, L = width, height, length
    b = _Skin()
    loops = []
    centres = []
    for i in range(n + 1):
        x = L * i / n
        d2 = d2a + (d2b - d2a) * x / L
        top, side, apex = (
            np.array(v) for v in ([x, W - d1, H], [x, W, H - d2], [x, W, H])
        )
        loop = []
        for j in range(m):  # the bottom, y from 0 to W
            loop.append(([x, W * j / m, 0.0], [0.0, 0.0, -1.0]))
        for j in range(m):  # the side, z from 0 to H - d2
            loop.append(([x, W, (H - d2) * j / m], [0.0, 1.0, 0.0]))
        for j in range(k):  # the blend, side foot to top foot
            u = j / k  # evenly along a chamfer's chord, by the conic's parameter else
            if weight > 0.0:
                b0, b1, b2 = (1.0 - u) ** 2, 2.0 * weight * u * (1.0 - u), u * u
            else:
                b0, b1, b2 = 1.0 - u, 0.0, u
            loop.append(((b0 * side + b1 * apex + b2 * top) / (b0 + b1 + b2), None))
        for j in range(m):  # the top, y from W - d1 to 0
            loop.append(([x, (W - d1) * (1.0 - j / m), H], [0.0, 0.0, 1.0]))
        for j in range(m):  # the back, z from H to 0
            loop.append(([x, 0.0, H * (1.0 - j / m)], [0.0, -1.0, 0.0]))
        loops.append([(b.add(p), o) for p, o in loop])
        centres.append(np.array([x, W - d1, H - d2]))
    count = len(loops[0])
    for i in range(n):
        for j in range(count):
            jn = (j + 1) % count
            a0, out = loops[i][j]
            if out is None:
                out = np.asarray(b.pts[a0]) - centres[i]
            b.quad(a0, loops[i + 1][j][0], loops[i + 1][jn][0], loops[i][jn][0], out)
    for i, sign in ((0, -1.0), (n, 1.0)):
        ring = [v for v, _ in loops[i]]
        mid = b.add(np.mean([b.pts[v] for v in ring], axis=0))
        for j in range(count):
            b.tri(mid, ring[j], ring[(j + 1) % count], np.array([sign, 0.0, 0.0]))
    mesh = b.mesh()
    m = mesh.points[mesh.cells[0].data].mean(axis=1)
    e = 1e-9
    label = np.full(len(m), 6, dtype=np.int64)
    for code, hit in reversed(
        [
            (0, np.abs(m[:, 2]) < e),
            (1, np.abs(m[:, 1] - W) < e),
            (2, np.abs(m[:, 2] - H) < e),
            (3, np.abs(m[:, 1]) < e),
            (4, np.abs(m[:, 0]) < e),
            (5, np.abs(m[:, 0] - L) < e),
        ]
    ):
        label[hit] = code
    mesh.cell_data["segment:patch"] = [label]
    return mesh


def rounded_box(r=0.3, sx=1.0, sy=1.3, sz=0.8, nf=4, nb=6):
    """The box ``[-sx, sx] x [-sy, sy] x [-sz, sz]`` with every edge and corner
    rounded: the cube ``[-1, 1]^3`` rounded by a ball of radius ``r``, scaled.
    Scaled unequally, its edge blends are conic blends and its corners
    ellipsoid patches, each a vertex blend tangent to three edge blends; the
    triangles carry the patches as ``segment:patch`` (``analytic::rounded_box``)."""
    inner = 1.0 - r
    g = [-(inner + r * np.tan(0.25 * np.pi * k / nb)) for k in range(nb, 0, -1)]
    g += [-inner + 2.0 * inner * k / nf for k in range(nf + 1)]
    g += [inner + r * np.tan(0.25 * np.pi * k / nb) for k in range(1, nb + 1)]
    scale = np.array([sx, sy, sz])
    b = _Skin()
    ids = {}

    def vertex(x):
        key = tuple(int(round(v * 1e9)) for v in x)
        if key not in ids:
            x = np.asarray(x, dtype=float)
            c = np.clip(x, -inner, inner)
            d = x - c
            ids[key] = b.add(scale * (c + r * d / np.linalg.norm(d)))
        return ids[key]

    def patch_of(c):
        code = 0
        for x in c:
            code = code * 3 + (2 if x > inner else 0 if x < -inner else 1)
        return code

    label = []
    n = len(g)
    for a in range(3):
        for sign in (-1.0, 1.0):
            u, v = (a + 1) % 3, (a + 2) % 3

            def at(i, j):
                x = [0.0, 0.0, 0.0]
                x[a], x[u], x[v] = sign, g[i], g[j]
                return x

            out = np.zeros(3)
            out[a] = sign
            for i in range(n - 1):
                for j in range(n - 1):
                    c = 0.5 * (np.asarray(at(i, j)) + np.asarray(at(i + 1, j + 1)))
                    b.quad(
                        vertex(at(i, j)),
                        vertex(at(i + 1, j)),
                        vertex(at(i + 1, j + 1)),
                        vertex(at(i, j + 1)),
                        out,
                    )
                    label += [patch_of(c)] * 2
    mesh = b.mesh()
    mesh.cell_data["segment:patch"] = [np.array(label, dtype=np.int64)]
    return mesh


def variable_fillet_block(
    length=2.0, width=1.0, height=1.0, rm=0.25, ra=0.05, n=48, k=16, m=6
):
    """The block ``[0, length] x [0, width] x [0, height]`` whose edge at
    ``y = width``, ``z = height`` is rounded by a rolling ball of radius
    ``r(x) = rm + ra cos(pi x / length)``, varying along it. Each row of the
    fillet is the exact characteristic circle of the ball family, from its
    foot on the side to its foot on the top; ``r'`` vanishes at both ends, so
    the end caps are planes (``analytic::variable_fillet_block``)."""
    W, H, L = width, height, length
    b = _Skin()
    loops = []
    centres = []
    for i in range(n + 1):
        x = L * i / n
        r = rm + ra * np.cos(np.pi * x / L)
        dr = -ra * np.pi / L * np.sin(np.pi * x / L)
        c = np.array([x, W - r, H - r])
        cd = np.array([1.0, -dr, -dr])
        t = cd / np.linalg.norm(cd)
        delta = -r * dr / np.linalg.norm(cd)
        cc = c + delta * t
        rho = np.sqrt(r * r - delta * delta)
        f_side, f_top = np.array([x, W, H - r]), np.array([x, W - r, H])
        a = (f_side - cc) / rho
        b2 = f_top - cc
        bv = b2 - np.dot(b2, a) * a
        bv /= np.linalg.norm(bv)
        om = np.arctan2(np.dot(b2, bv), np.dot(b2, a))
        loop = []
        for j in range(m):  # the bottom, y from 0 to W
            loop.append(([x, W * j / m, 0.0], [0.0, 0.0, -1.0]))
        for j in range(m):  # the side, z from 0 to H - r
            loop.append(([x, W, (H - r) * j / m], [0.0, 1.0, 0.0]))
        for j in range(k):  # the fillet, side foot to top foot
            p = cc + rho * (np.cos(om * j / k) * a + np.sin(om * j / k) * bv)
            loop.append((p, None))
        for j in range(m):  # the top, y from W - r to 0
            loop.append(([x, (W - r) * (1.0 - j / m), H], [0.0, 0.0, 1.0]))
        for j in range(m):  # the back, z from H to 0
            loop.append(([x, 0.0, H * (1.0 - j / m)], [0.0, -1.0, 0.0]))
        loops.append([(b.add(p), o) for p, o in loop])
        centres.append(c)
    count = len(loops[0])
    for i in range(n):
        for j in range(count):
            jn = (j + 1) % count
            a0, out = loops[i][j]
            if out is None:
                out = np.asarray(b.pts[a0]) - centres[i]
            b.quad(a0, loops[i + 1][j][0], loops[i + 1][jn][0], loops[i][jn][0], out)
    for i, sign in ((0, -1.0), (n, 1.0)):
        ring = [v for v, _ in loops[i]]
        mid = b.add(np.mean([b.pts[v] for v in ring], axis=0))
        for j in range(count):
            b.tri(mid, ring[j], ring[(j + 1) % count], np.array([sign, 0.0, 0.0]))
    return b.mesh()


def wavy_bar(
    length=3.0, wave=0.06, twist=3.0, r=0.2, back=1.0, depth=1.0, n=96, k=8, m=8
):
    """A bar along x whose top and front faces are both freeform, meeting at an
    edge rounded by a ball of radius ``r`` rolling along the spine ``(x, 0,
    wave cos(4 pi x / L))``: the top waves with the spine, and the front leans
    in and out by up to ``twist`` degrees, so the angle between the two faces
    changes along the edge. Both faces are ruled surfaces through the ball's
    feet, their rulings normal to its radius there, so tangent to every ball
    by construction; the back, the bottom and the ends are planes
    (``analytic::wavy_bar``)."""
    L = length
    b = _Skin()
    loops = []
    centres = []
    for i in range(n + 1):
        x = L * i / n
        c = np.array([x, 0.0, wave * np.cos(4 * np.pi * x / L)])
        cd = np.array([1.0, 0.0, -wave * 4 * np.pi / L * np.sin(4 * np.pi * x / L)])
        t = cd / np.linalg.norm(cd)
        side = np.array([0.0, 1.0, 0.0])  # normal to t, horizontal
        up = np.cross(t, side)  # normal to t, upwards
        psi = np.radians(twist) * np.sin(2 * np.pi * x / L)
        m1 = up
        m2 = np.cos(psi) * side + np.sin(psi) * up
        q1, q2 = c + r * m1, c + r * m2
        # The rulings, away from the edge, normal to the ball's radius at the
        # foot and tilted along the spine so they keep their x.
        w1 = np.cross(m1, t)
        w1 = w1 if w1[1] < 0.0 else -w1
        w1 = w1 - w1[0] / t[0] * t
        w2 = np.cross(m2, t)
        w2 = w2 if w2[2] < 0.0 else -w2
        w2 = w2 - w2[0] / t[0] * t
        top_far = q1 + (q1[1] + back) / -w1[1] * w1
        front_low = q2 + (q2[2] + depth) / -w2[2] * w2
        corner = np.array([top_far[0], -back, -depth])
        om = np.arccos(np.clip(np.dot(m1, m2), -1.0, 1.0))
        loop = []
        for j in range(m):  # the back, up to the top's far edge
            loop.append((corner + (top_far - corner) * j / m, [0.0, -1.0, 0.0]))
        for j in range(m):  # the top, towards the edge
            loop.append((top_far + (q1 - top_far) * j / m, m1))
        for j in range(k):  # the blend, top foot to front foot
            f0 = np.sin((1.0 - j / k) * om) / np.sin(om)
            f1 = np.sin(j / k * om) / np.sin(om)
            loop.append((c + r * (f0 * m1 + f1 * m2), None))
        for j in range(m):  # the front, down to the bottom
            loop.append((q2 + (front_low - q2) * j / m, m2))
        for j in range(m):  # the bottom, back to the corner
            loop.append((front_low + (corner - front_low) * j / m, [0.0, 0.0, -1.0]))
        loops.append([(b.add(p), o) for p, o in loop])
        centres.append(c)
    count = len(loops[0])
    for i in range(n):
        for j in range(count):
            jn = (j + 1) % count
            a0, out = loops[i][j]
            if out is None:
                out = np.asarray(b.pts[a0]) - centres[i]
            b.quad(a0, loops[i + 1][j][0], loops[i + 1][jn][0], loops[i][jn][0], out)
    for i, sign in ((0, -1.0), (n, 1.0)):
        ring = [v for v, _ in loops[i]]
        mid = b.add(np.mean([b.pts[v] for v in ring], axis=0))
        for j in range(count):
            b.tri(mid, ring[j], ring[(j + 1) % count], np.array([sign, 0.0, 0.0]))
    return b.mesh()


def filleted_plate(a=2.0, h=0.5, r=0.6, rc=0.8, n=48, k=8):
    """The plate ``[-a, a]^2 x [0, h]`` with a through-hole of radius ``r`` at
    its centre and its corner ``(a, a)`` rounded by a quarter circle of radius
    ``rc`` tangent to both sides: a CAD part tessellated as CAD exporters do,
    rings of vertices on the hole and on the outline (``k`` segments along the
    fillet) joined by strips, every vertex on its surfaces, the fillet meeting
    its sides at smooth edges. Six planes and two cylinders
    (``analytic::filleted_plate``)."""
    c = a - rc
    angles = {2.0 * np.pi * i / n for i in range(n)}
    angles |= {np.pi * (0.75 + 0.5 * j) for j in range(3)}  # the sharp corners
    for j in range(k + 1):  # along the fillet, from tangent point to tangent point
        phi = 0.5 * np.pi * j / k
        angles.add(np.arctan2(c + rc * np.sin(phi), c + rc * np.cos(phi)))
    angles = sorted(angles)
    lo_t, hi_t = np.arctan2(c, a), np.arctan2(a, c)

    def outline(th):
        d = np.array([np.cos(th), np.sin(th)])
        if lo_t <= th <= hi_t:
            dc = d @ np.array([c, c])
            return (dc + np.sqrt(dc * dc - 2.0 * c * c + rc * rc)) * d
        return a / max(abs(d[0]), abs(d[1])) * d

    b = _Skin()
    rings = []
    for th in angles:
        d = np.array([np.cos(th), np.sin(th)])
        o = outline(th)
        rings.append(
            (
                b.add([r * d[0], r * d[1], 0.0]),
                b.add([r * d[0], r * d[1], h]),
                b.add([o[0], o[1], 0.0]),
                b.add([o[0], o[1], h]),
            )
        )
    m = len(rings)
    for i in range(m):
        hb, ht, ob, ot = rings[i]
        hb2, ht2, ob2, ot2 = rings[(i + 1) % m]
        mid = 0.5 * (
            angles[i] + (angles[(i + 1) % m] + (2.0 * np.pi if i + 1 == m else 0.0))
        )
        d = np.array([np.cos(mid), np.sin(mid), 0.0])
        b.quad(ht, ot, ot2, ht2, [0.0, 0.0, 1.0])
        b.quad(hb, ob, ob2, hb2, [0.0, 0.0, -1.0])
        b.quad(ob, ob2, ot2, ot, d)
        b.quad(hb, hb2, ht2, ht, -d)
    return b.mesh()


# ---------------------------------------------------------------------------
# Freeform fixtures (v17.12.0): the shapes ``analytic_surfaces.hpp`` builds for
# fit_cad's freeform modes, with the same arithmetic in the same order on
# Python floats (``math``, not NumPy, whose SIMD kernels may differ in the last
# bits), so the points match the C++ ones bit for bit on one platform.
# ---------------------------------------------------------------------------


def triangle_box(n, warp=False):
    """The unit cube's skin as triangles (``analytic::triangle_box``)."""
    quads = box_skin(n)
    pts = []
    for x, y, z in quads.points.tolist():
        if warp:
            x += 0.35 * x * (1.0 - x) * (1.0 - 2.0 * x)
            y += 0.35 * y * (1.0 - y) * (1.0 - 2.0 * y)
        pts.append([x, y, z])
    tris = []
    for a, b, c, d in quads.cells[0].data.tolist():
        tris.append([a, b, c])
        tris.append([a, c, d])
    return mio.Mesh(
        np.array(pts, dtype=float), [("triangle", np.array(tris, dtype=np.int64))]
    )


def heightfield_box(n, height, warp=False):
    """A box whose top (z = 1) is lifted to ``1 + height(x, y)``, its sides'
    points following in proportion (``z (1 + height)``) so they stay planar
    and never fold (``analytic::heightfield_box``)."""
    m = triangle_box(n, warp)
    pts = m.points.copy()
    for i in range(len(pts)):
        if pts[i, 2] > 0.0:
            pts[i, 2] = float(pts[i, 2]) * (
                1.0 + height(float(pts[i, 0]), float(pts[i, 1]))
            )
    m.points = pts
    return m


def monkey_saddle_box(n=24, warp=False):
    def f(x, y):
        u, v = 2.0 * x - 1.0, 2.0 * y - 1.0
        return 0.15 * (u * u * u - 3.0 * u * v * v)

    return heightfield_box(n, f, warp)


def gaussian_bump_box(n=24, warp=False):
    def f(x, y):
        r2 = (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5)
        return 0.3 * math.exp(-r2 / 0.03)

    return heightfield_box(n, f, warp)


def s_band_box(n=24, warp=False):
    def f(x, y):
        s = x * x * x * (x * (6.0 * x - 15.0) + 10.0)
        return 0.25 * (s - 0.5) * (1.0 + 0.5 * (y - 0.5))

    return heightfield_box(n, f, warp)


def hood_box(n=24, warp=False):
    def f(x, y):
        return 0.35 * 16.0 * x * (1.0 - x) * y * (1.0 - y) * (1.0 + 0.4 * x)

    return heightfield_box(n, f, warp)


def organic_torus(n=48, m=24):
    """A genus-1 blob: a torus whose tube radius ripples around it."""
    pts, tris = [], []
    for i in range(m):
        t = 2.0 * math.pi * i / m
        for k in range(n):
            a = 2.0 * math.pi * k / n
            r = 0.35 * (1.0 + 0.2 * math.cos(2.0 * t) * math.sin(3.0 * a))
            rr = 1.0 + r * math.cos(t)
            pts.append([rr * math.cos(a), rr * math.sin(a), 0.7 * r * math.sin(t)])
    for i in range(m):
        for k in range(n):
            i1, k1 = (i + 1) % m, (k + 1) % n
            a, b, c, d = i * n + k, i * n + k1, i1 * n + k1, i1 * n + k
            tris.append([a, b, c])
            tris.append([a, c, d])
    return mio.Mesh(np.array(pts), [("triangle", np.array(tris, dtype=np.int64))])


def implicit_skin(field, lo, hi, n):
    """The zero level of ``field`` (negative inside) by marching tetrahedra on
    the six Kuhn tetrahedra of each grid cube (``analytic::implicit_skin``)."""
    m = n + 1

    def node(i, j, k):
        return (k * m + j) * m + i

    xyz = [None] * (m * m * m)
    f = [0.0] * (m * m * m)
    for k in range(m):
        for j in range(m):
            for i in range(m):
                p = [
                    lo[0] + (hi[0] - lo[0]) * float(i) / float(n),
                    lo[1] + (hi[1] - lo[1]) * float(j) / float(n),
                    lo[2] + (hi[2] - lo[2]) * float(k) / float(n),
                ]
                xyz[node(i, j, k)] = p
                f[node(i, j, k)] = field(*p)
    b = _Skin()
    crossing = {}

    def cross(a, c):
        key = (min(a, c), max(a, c))
        if key not in crossing:
            inside = key[0] if f[key[0]] < 0.0 else key[1]
            outside = key[1] if inside == key[0] else key[0]
            fi, fo = f[inside], f[outside]
            t = fi / (fi - fo)
            p, q = xyz[inside], xyz[outside]
            crossing[key] = b.add([p[a2] + t * (q[a2] - p[a2]) for a2 in range(3)])
        return crossing[key]

    orders = ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0))
    for k in range(n):
        for j in range(n):
            for i in range(n):
                for o in orders:
                    c = [i, j, k]
                    v = [node(*c)]
                    for axis in o:
                        c[axis] += 1
                        v.append(node(*c))
                    ins = [x for x in v if f[x] < 0.0]
                    outs = [x for x in v if f[x] >= 0.0]
                    if not ins or not outs:
                        continue
                    out = [0.0, 0.0, 0.0]
                    for a in range(3):
                        out[o[a]] = (f[v[a + 1]] - f[v[a]]) / (hi[o[a]] - lo[o[a]])
                    if len(ins) == 1 or len(outs) == 1:
                        apex, rest = (ins[0], outs) if len(ins) == 1 else (outs[0], ins)
                        b.tri(
                            cross(apex, rest[0]),
                            cross(apex, rest[1]),
                            cross(apex, rest[2]),
                            out,
                        )
                    else:
                        b.quad(
                            cross(ins[0], outs[0]),
                            cross(ins[0], outs[1]),
                            cross(ins[1], outs[1]),
                            cross(ins[1], outs[0]),
                            out,
                        )
    return b.mesh()


def sphere_with_handle(n=32):
    def f(x, y, z):
        sphere = math.sqrt(x * x + y * y + z * z) - 0.6
        q = math.sqrt((x - 0.6) * (x - 0.6) + z * z) - 0.42
        ring = math.sqrt(q * q + y * y) - 0.11
        k = 0.08
        h = max(k - abs(sphere - ring), 0.0) / k
        return min(sphere, ring) - h * h * k * 0.25

    return implicit_skin(f, (-0.8123, -0.8077, -0.7951), (1.2111, 0.8173, 0.8029), n)


def double_torus(n=40):
    def f(x, y, z):
        r2 = x * x + y * y
        lem = r2 * r2 - x * x + y * y
        return lem * lem + 4.0 * z * z - 0.02

    return implicit_skin(f, (-1.3117, -0.6531, -0.4073), (1.3089, 0.6563, 0.4049), n)


def catmull_clark_limit(points, faces, levels):
    """Catmull-Clark steps, then the limit stencil (``analytic::catmull_clark_limit``)."""
    P = [list(map(float, p)) for p in points]
    F = [list(f) for f in faces]
    for _ in range(levels):
        nv = len(P)
        fp = []
        for f in F:
            c = [0.0, 0.0, 0.0]
            for v in f:
                for j in range(3):
                    c[j] += P[v][j]
            fp.append([c[j] / len(f) for j in range(3)])
        edges = {}
        for fi, f in enumerate(F):
            for k in range(len(f)):
                a, b = f[k], f[(k + 1) % len(f)]
                edges.setdefault((min(a, b), max(a, b)), []).append(fi)
        nxt = [None] * (nv + len(F) + len(edges))
        for fi in range(len(F)):
            nxt[nv + fi] = fp[fi]
        edge_id = {}
        e = nv + len(F)
        fsum = [[0.0, 0.0, 0.0] for _ in range(nv)]
        rsum = [[0.0, 0.0, 0.0] for _ in range(nv)]
        valence, nfaces = [0] * nv, [0] * nv
        for key in sorted(edges):
            fa, fb = edges[key]
            a, b = P[key[0]], P[key[1]]
            nxt[e] = [0.25 * (a[j] + b[j] + fp[fa][j] + fp[fb][j]) for j in range(3)]
            edge_id[key] = e
            e += 1
            for v in key:
                valence[v] += 1
                for j in range(3):
                    rsum[v][j] += 0.5 * (a[j] + b[j])
        for fi, f in enumerate(F):
            for v in f:
                nfaces[v] += 1
                for j in range(3):
                    fsum[v][j] += fp[fi][j]
        for v in range(nv):
            n = float(valence[v])
            nxt[v] = [
                (fsum[v][j] / nfaces[v] + 2.0 * rsum[v][j] / n + (n - 3.0) * P[v][j])
                / n
                for j in range(3)
            ]
        quads = []
        for fi, f in enumerate(F):
            k = len(f)
            for i in range(k):
                prev, cur, nx = f[(i + k - 1) % k], f[i], f[(i + 1) % k]
                quads.append(
                    [
                        cur,
                        edge_id[(min(cur, nx), max(cur, nx))],
                        nv + fi,
                        edge_id[(min(prev, cur), max(prev, cur))],
                    ]
                )
        P, F = nxt, quads
    esum = [[0.0, 0.0, 0.0] for _ in P]
    dsum = [[0.0, 0.0, 0.0] for _ in P]
    val = [0] * len(P)
    for q in F:
        for i in range(4):
            v, nb, dg = q[i], P[q[(i + 1) % 4]], P[q[(i + 2) % 4]]
            val[v] += 1
            for j in range(3):
                esum[v][j] += nb[j]
                dsum[v][j] += dg[j]
    pts = []
    for v in range(len(P)):
        n = float(val[v])
        pts.append(
            [
                (n * n * P[v][j] + 4.0 * esum[v][j] + dsum[v][j]) / (n * (n + 5.0))
                for j in range(3)
            ]
        )
    tris = []
    for q in F:
        tris.append([q[0], q[1], q[2]])
        tris.append([q[0], q[2], q[3]])
    return mio.Mesh(np.array(pts), [("triangle", np.array(tris, dtype=np.int64))])


def cc_cage(valence, levels=4):
    """A Catmull-Clark limit surface with extraordinary vertices of ``valence``."""
    if valence == 3:
        P = [
            [-1, -1, -1],
            [1, -1, -1],
            [1, 1, -1],
            [-1, 1, -1],
            [-1, -1, 1],
            [1, -1, 1],
            [1, 1, 1],
            [-1, 1, 1],
        ]
        F = [
            [0, 3, 2, 1],
            [4, 5, 6, 7],
            [0, 1, 5, 4],
            [1, 2, 6, 5],
            [2, 3, 7, 6],
            [3, 0, 4, 7],
        ]
    else:
        n = valence
        P = []
        for z in (0, 1):
            for k in range(n):
                a = 2.0 * math.pi * float(k) / float(n)
                P.append([math.cos(a), math.sin(a), 0.8 if z else -0.8])
        F = [[n - 1 - k for k in range(n)], [n + k for k in range(n)]]
        F += [[k, (k + 1) % n, n + (k + 1) % n, n + k] for k in range(n)]
    return catmull_clark_limit(P, F, levels)


def g1_loaf(n=32, length=24):
    """A superellipse loaf: top and curved sides meeting G1 (``analytic::g1_loaf``)."""
    b = _Skin()
    rings = []
    for j in range(length + 1):
        y = float(j) / length
        swell = 1.0 + 0.25 * 4.0 * y * (1.0 - y)
        ring = []
        for i in range(n + 1):
            t = math.pi * i / n
            c, s = math.cos(t), math.sin(t)
            x = (-1.0 if c < 0.0 else 1.0) * math.sqrt(abs(c))
            z = math.sqrt(abs(s))
            ring.append(b.add([0.5 * x, 2.0 * y, 0.6 * swell * z]))
        rings.append(ring)
    for j in range(length):
        for i in range(n):
            t = math.pi * (i + 0.5) / n
            b.quad(
                rings[j][i],
                rings[j][i + 1],
                rings[j + 1][i + 1],
                rings[j + 1][i],
                [math.cos(t), 0.0, math.sin(t)],
            )
    for j in range(length):
        b.quad(
            rings[j][0], rings[j][n], rings[j + 1][n], rings[j + 1][0], [0.0, 0.0, -1.0]
        )
    for e in (0, 1):
        j = length if e else 0
        y = 2.0 * j / length
        frac = float(j) / length
        swell = 1.0 + 0.25 * 4.0 * frac * (1.0 - frac)
        c = b.add([0.0, y, 0.3 * swell])
        out = [0.0, 1.0 if e else -1.0, 0.0]
        for i in range(n):
            b.tri(c, rings[j][i], rings[j][i + 1], out)
        b.tri(c, rings[j][n], rings[j][0], out)
    return b.mesh()


_MASK64 = (1 << 64) - 1


def _splitmix64(state):
    """One splitmix64 step: ``(new_state, value)``."""
    state = (state + 0x9E3779B97F4A7C15) & _MASK64
    z = state
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _MASK64
    return state, z ^ (z >> 31)


def gaussian_noise(mesh, sigma, seed=7):
    """Irwin-Hall Gaussian noise on every coordinate (``analytic::gaussian_noise``)."""
    pts = np.asarray(mesh.points, dtype=float).reshape(-1).tolist()
    state = seed
    for i in range(len(pts)):
        g = -6.0
        for _ in range(12):
            state, z = _splitmix64(state)
            g += float(z >> 11) * 2.0**-53
        pts[i] += sigma * g
    out = mio.Mesh(np.array(pts).reshape(-1, 3), [(b.type, b.data) for b in mesh.cells])
    return out
