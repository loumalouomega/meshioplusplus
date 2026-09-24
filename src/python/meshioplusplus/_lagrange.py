"""
Arbitrary-order Lagrange cells: the ``VTK_LAGRANGE_*`` node layout and nodal
interpolation between point sets of one polynomial space (the Python twin of
``src/cpp/src/formats/lagrange_common.hpp``).

Shapes are ``"line"``, ``"triangle"``, ``"quad"``, ``"tetra"``, ``"hexahedron"``
and ``"wedge"``, on the reference cells whose corners are meshio++'s (and
MFEM's). ``vtk_lattice`` lists a VTK Lagrange cell's nodes as integer lattice
coordinates in VTK's order (ported from MFEM's ``mesh/vtk.cpp``, BSD-3-Clause,
Lawrence Livermore National Security); ``interpolation_matrix`` maps values at
one unisolvent node set of the order-``p`` space to other points through an
orthogonal basis; ``lattice_weights`` names a node by the corners it combines,
so cells sharing an edge or face share its nodes.
"""

import math

import numpy as np

DIM = {"line": 1, "triangle": 2, "quad": 2, "tetra": 3, "hexahedron": 3, "wedge": 3}


def num_nodes(shape, p):
    return {
        "line": p + 1,
        "triangle": (p + 1) * (p + 2) // 2,
        "quad": (p + 1) ** 2,
        "tetra": (p + 1) * (p + 2) * (p + 3) // 6,
        "hexahedron": (p + 1) ** 3,
        "wedge": (p + 1) * (p + 1) * (p + 2) // 2,
    }[shape]


def gll_points(p):
    """The Gauss-Lobatto points of order ``p`` on [0, 1], ascending."""
    x = [0.0] * (p + 1)
    x[-1] = 1.0
    for k in range(1, p // 2 + 1):
        t = -math.cos(math.pi * k / p)
        for _ in range(100):
            p0, p1 = 1.0, t
            for n in range(1, p):
                p0, p1 = p1, ((2 * n + 1) * t * p1 - n * p0) / (n + 1)
            d1 = p * (t * p1 - p0) / (t * t - 1.0)
            d2 = (2.0 * t * d1 - p * (p + 1) * p1) / (1.0 - t * t)
            step = d1 / d2
            t -= step
            if abs(step) < 1e-19:
                break
        v = (1.0 + t) / 2.0
        x[k] = v
        x[p - k] = 1.0 - v
    if p % 2 == 0 and p > 0:
        x[p // 2] = 0.5
    return x


def uniform_points(p):
    return [i / p for i in range(p + 1)]


def _triangle_index(b, ref):
    mx, mn = ref, 0
    bmin = min(b[0], b[1], b[2])
    idx = 0
    while bmin > mn:
        idx += 3 * ref
        mx -= 2
        mn += 1
        ref -= 3
    for d in range(3):
        if b[(d + 2) % 3] == mx:
            return idx
        idx += 1
    for d in range(3):
        if b[(d + 1) % 3] == mn:
            return idx + b[d] - (mn + 1)
        idx += mx - (mn + 1)
    return idx


def _tetra_index(b, ref):
    idx, mx, mn = 0, ref, 0
    bmin = min(b)
    while bmin > mn:
        idx += 2 * (ref * ref + 1)
        mx -= 3
        mn += 1
        ref -= 4
    vertex_max = (3, 0, 1, 2)
    edge_min = ((1, 2), (2, 3), (0, 2), (0, 1), (1, 3), (0, 3))
    edge_count = (0, 1, 3, 2, 2, 2)
    face_min = (1, 3, 0, 2)
    face_b = ((0, 2, 3), (2, 0, 1), (2, 1, 3), (1, 0, 3))
    for v in range(4):
        if b[vertex_max[v]] == mx:
            return idx
        idx += 1
    for e in range(6):
        if b[edge_min[e][0]] == mn and b[edge_min[e][1]] == mn:
            return idx + b[edge_count[e]] - (mn + 1)
        idx += mx - (mn + 1)
    for f in range(4):
        if b[face_min[f]] == mn:
            projected = [b[face_b[f][i]] - mn for i in range(3)]
            return idx + _triangle_index(projected, ref) - 3 * ref
        idx += (ref + 1) * (ref + 2) // 2 - 3 * ref
    return idx


def _triangle_offset(ref, i, j):
    return i + ref * (j - 1) - (j * (j + 1)) // 2


def _wedge_index(i, j, k, ref):
    om1 = ref - 1
    ibdr, jbdr, ijbdr = int(i == 0), int(j == 0), int(i + j == ref)
    kbdr = int(k == 0 or k == ref)
    nbdr = ibdr + jbdr + ijbdr + kbdr
    if nbdr == 3:
        return (0 if ibdr and jbdr else (1 if jbdr and ijbdr else 2)) + (3 if k else 0)
    offset = 6
    if nbdr == 2:
        if not kbdr:
            offset += om1 * 6
            return (
                offset
                + (k - 1)
                + (0 if ibdr and jbdr else (1 if jbdr and ijbdr else 2)) * om1
            )
        offset += 3 * om1 if k == ref else 0
        if jbdr:
            return offset + i - 1
        offset += om1
        if ijbdr:
            return offset + j - 1
        offset += om1
        return offset + (ref - j - 1)
    offset += 9 * om1
    ntf = (om1 - 1) * om1 // 2
    nqf = om1 * om1
    if nbdr == 1:
        if kbdr:
            if k > 0:
                offset += ntf
            return offset + _triangle_offset(ref, i, j)
        offset += 2 * ntf
        if jbdr:
            return offset + (i - 1) + om1 * (k - 1)
        offset += nqf
        if ijbdr:
            return offset + (j - 1) + om1 * (k - 1)
        offset += nqf
        return offset + (ref - j - 1) + om1 * (k - 1)
    offset += 2 * ntf + 3 * nqf
    return offset + _triangle_offset(ref, i, j) + ntf * (k - 1)


def _tensor_index(idx, ref, dim):
    n = ref + 1
    if dim == 1:
        if idx in (0, ref):
            return 1 if idx else 0
        return idx + 1
    if dim == 2:
        i, j = idx % n, idx // n
        ibdr, jbdr = i in (0, ref), j in (0, ref)
        if ibdr and jbdr:
            return (2 if j else 1) if i else (3 if j else 0)
        offset = 4
        if jbdr:
            return (i - 1) + (2 * (ref - 1) if j else 0) + offset
        if ibdr:
            return (j - 1) + ((ref - 1) if i else 3 * (ref - 1)) + offset
        offset += 4 * (ref - 1)
        return offset + (i - 1) + (ref - 1) * (j - 1)
    i, j, k = idx % n, (idx // n) % n, idx // (n * n)
    ibdr, jbdr, kbdr = i in (0, ref), j in (0, ref), k in (0, ref)
    nbdr = int(ibdr) + int(jbdr) + int(kbdr)
    if nbdr == 3:
        return ((2 if j else 1) if i else (3 if j else 0)) + (4 if k else 0)
    offset = 8
    if nbdr == 2:
        if not ibdr:
            return (
                (i - 1)
                + (2 * (ref - 1) if j else 0)
                + (4 * (ref - 1) if k else 0)
                + offset
            )
        if not jbdr:
            return (
                (j - 1)
                + ((ref - 1) if i else 3 * (ref - 1))
                + (4 * (ref - 1) if k else 0)
                + offset
            )
        offset += 8 * (ref - 1)
        return (
            (k - 1) + (ref - 1) * ((2 if j else 1) if i else (3 if j else 0)) + offset
        )
    offset += 12 * (ref - 1)
    q = (ref - 1) * (ref - 1)
    if nbdr == 1:
        if ibdr:
            return (j - 1) + (ref - 1) * (k - 1) + (q if i else 0) + offset
        offset += 2 * q
        if jbdr:
            return (i - 1) + (ref - 1) * (k - 1) + (q if j else 0) + offset
        offset += 2 * q
        return (i - 1) + (ref - 1) * (j - 1) + (q if k else 0) + offset
    offset += 6 * q
    return offset + (i - 1) + (ref - 1) * ((j - 1) + (ref - 1) * (k - 1))


def vtk_lattice(shape, p):
    """The nodes of an order-``p`` VTK Lagrange cell, in VTK order, as integer
    lattice coordinates ``(i, j, k)`` (``x = i / p``)."""
    out = [None] * num_nodes(shape, p)
    if shape == "triangle":
        for b1 in range(p + 1):
            for b0 in range(p - b1 + 1):
                out[_triangle_index([b0, b1, p - b0 - b1], p)] = (b0, b1, 0)
    elif shape == "tetra":
        for b2 in range(p + 1):
            for b1 in range(p - b2 + 1):
                for b0 in range(p - b1 - b2 + 1):
                    b = [b0, b1, b2, p - b0 - b1 - b2]
                    out[_tetra_index(b, p)] = (b0, b1, b2)
    elif shape == "wedge":
        for k in range(p + 1):
            for j in range(p + 1):
                for i in range(p - j + 1):
                    out[_wedge_index(i, j, k, p)] = (i, j, k)
    else:
        dim = DIM[shape]
        n = p + 1
        for idx in range(n**dim):
            ijk = (
                idx % n,
                (idx // n) % n if dim > 1 else 0,
                idx // (n * n) if dim > 2 else 0,
            )
            out[_tensor_index(idx, p, dim)] = ijk
    return out


def lattice_weights(shape, p, ijk):
    """``(corner, weight)`` of every corner lattice node ``ijk`` combines, over
    the common denominator ``p**3``."""
    i, j, k = ijk
    if shape == "line":
        w = [(p - i) * p * p, i * p * p]
    elif shape == "triangle":
        w = [(p - i - j) * p * p, i * p * p, j * p * p]
    elif shape == "quad":
        w = [(p - i) * (p - j) * p, i * (p - j) * p, i * j * p, (p - i) * j * p]
    elif shape == "tetra":
        w = [(p - i - j - k) * p * p, i * p * p, j * p * p, k * p * p]
    elif shape == "hexahedron":
        w = [
            (p - i) * (p - j) * (p - k),
            i * (p - j) * (p - k),
            i * j * (p - k),
            (p - i) * j * (p - k),
            (p - i) * (p - j) * k,
            i * (p - j) * k,
            i * j * k,
            (p - i) * j * k,
        ]
    else:
        w = [
            (p - i - j) * (p - k) * p,
            i * (p - k) * p,
            j * (p - k) * p,
            (p - i - j) * k * p,
            i * k * p,
            j * k * p,
        ]
    return [(c, v) for c, v in enumerate(w) if v]


def _jacobi(n, a, x):
    v = [1.0]
    if n >= 1:
        v.append(((a + 2.0) * x + a) / 2.0)
    for m in range(1, n):
        a1 = 2.0 * (m + 1) * (m + a + 1) * (2 * m + a)
        a2 = (2 * m + a + 1) * a * a
        a3 = (2 * m + a) * (2 * m + a + 1) * (2 * m + a + 2)
        a4 = 2.0 * (m + a) * m * (2 * m + a + 2)
        v.append(((a2 + a3 * x) * v[m] - a4 * v[m - 1]) / a1)
    return v


def _basis(shape, p, x):
    r, s, t = 2.0 * x[0] - 1.0, 2.0 * x[1] - 1.0, 2.0 * x[2] - 1.0
    if shape == "line":
        return _jacobi(p, 0.0, r)
    if shape == "quad":
        a, b = _jacobi(p, 0.0, r), _jacobi(p, 0.0, s)
        return [a[i] * b[j] for j in range(p + 1) for i in range(p + 1)]
    if shape == "hexahedron":
        a, b, c = _jacobi(p, 0.0, r), _jacobi(p, 0.0, s), _jacobi(p, 0.0, t)
        return [
            a[i] * b[j] * c[k]
            for k in range(p + 1)
            for j in range(p + 1)
            for i in range(p + 1)
        ]
    if shape in ("triangle", "wedge"):
        a = 2.0 * (1.0 + r) / (1.0 - s) - 1.0 if abs(1.0 - s) > 1e-14 else -1.0
        pa = _jacobi(p, 0.0, a)
        tri = []
        for i in range(p + 1):
            pb = _jacobi(p - i, 2.0 * i + 1.0, s)
            f = ((1.0 - s) / 2.0) ** i
            for j in range(p - i + 1):
                tri.append(pa[i] * f * pb[j])
        if shape == "triangle":
            return tri
        c = _jacobi(p, 0.0, t)
        return [v * c[k] for k in range(p + 1) for v in tri]
    a = 2.0 * (1.0 + r) / (-s - t) - 1.0 if abs(s + t) > 1e-14 else -1.0
    b = 2.0 * (1.0 + s) / (1.0 - t) - 1.0 if abs(1.0 - t) > 1e-14 else -1.0
    pa = _jacobi(p, 0.0, a)
    out = []
    for i in range(p + 1):
        pb = _jacobi(p - i, 2.0 * i + 1.0, b)
        fb = ((1.0 - b) / 2.0) ** i
        for j in range(p - i + 1):
            pc = _jacobi(p - i - j, 2.0 * (i + j) + 2.0, t)
            fc = ((1.0 - t) / 2.0) ** (i + j)
            for k in range(p - i - j + 1):
                out.append(pa[i] * fb * pb[j] * fc * pc[k])
    return out


def interpolation_matrix(shape, p, nodes, targets):
    """``(targets, nodes)`` matrix mapping values at ``nodes`` (a unisolvent set
    of the order-``p`` space) to values at ``targets``."""
    if len(nodes) != num_nodes(shape, p):
        raise ValueError("lagrange: node count does not match the space")
    v = np.array([_basis(shape, p, x) for x in nodes])  # nodes x basis
    phi = np.array([_basis(shape, p, x) for x in targets])  # targets x basis
    out = np.linalg.solve(v.T, phi.T).T
    out[np.abs(out) < 1e-14] = 0.0
    out[np.abs(out - 1.0) < 1e-14] = 1.0
    return out
