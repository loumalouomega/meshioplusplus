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


def _helical_rod(w, rad, pitch, turns, n, depth=1.0):
    """A rod about z whose surface is ``r = f(z - K theta)``, ``K = pitch /
    (2 pi)``, one period of ``f`` given as nodes ``(w, rad)``; the grid
    follows the helices through every node, and each end is closed by a fan
    to a point on the axis ``depth`` pitches inside the rod."""
    k = pitch / (2.0 * np.pi)
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
