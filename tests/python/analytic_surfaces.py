"""Analytic surfaces for the segment / fit_primitives / step tests.

Quads on a unit-cube skin and on surfaces of revolution, built from indices
so every point lies on its surface (up to rounding). The C++ tests build the
same shapes in ``tests/cpp/analytic_surfaces.hpp``.
"""

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


def threaded_rod(pitch=0.5, r0=0.8, r1=1.0, turns=4, n=64):
    """A rod with a trapezoidal thread: ``r = f(z - K theta)``, root and crest
    helical strips of cylinders, both flanks ruled helicoids, closed by two
    shallow fans a pitch inside the rod (``analytic::threaded_rod``)."""
    k = pitch / (2.0 * np.pi)
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
    per = len(w)
    rows = per * turns
    b = _Skin()
    ids = []
    for i in range(n):
        th = 2.0 * np.pi * i / n
        col = []
        for j in range(rows + 1):
            wj = w[j % per] + pitch * (j // per)
            rj = rad[j % per]
            col.append(b.add([rj * np.cos(th), rj * np.sin(th), wj + k * th]))
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
        c = b.add([0.0, 0.0, max(zs) + pitch if reverse else min(zs) - pitch])
        for i in range(len(loop)):
            p, q = loop[i], loop[(i + 1) % len(loop)]
            b.tris.append([c, q, p] if reverse else [c, p, q])

    low = [node(i, 0) for i in range(n)] + [node(0, j) for j in range(per, 0, -1)]
    high = [node(i, rows) for i in range(n)]
    high += [node(n - 1, j) for j in range(rows - 1, rows - per - 1, -1)]
    fan(low, True)
    fan(high, False)
    return b.mesh()
