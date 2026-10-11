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
