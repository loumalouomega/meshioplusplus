"""
Point evaluation of B-spline and NURBS curves and tensor-product surfaces at
arbitrary parameters (The NURBS Book, algorithms A2.1, A2.2, A4.1, A4.3): the
Python twin of ``src/cpp/src/detail/nurbs_eval.hpp``.

Knot vectors are in the full form (``n + p + 1`` knots), control nets list
direction 0 fastest, and ``weights=None`` evaluates the polynomial B-spline
with no division. Every operation runs in the C++ kernel's order on Python
floats, so the two give the same bits.
"""


def find_span(p, knots, n, u):
    """The knot span ``i`` (``p <= i <= n - 1``) with ``U[i] <= u < U[i+1]``."""
    last = n - 1
    if u >= knots[last + 1]:
        s = last
        while s > p and knots[s] == knots[s + 1]:
            s -= 1
        return s
    if u <= knots[p]:
        s = p
        while s < last and knots[s] == knots[s + 1]:
            s += 1
        return s
    low, high = p, last + 1
    mid = (low + high) // 2
    while u < knots[mid] or u >= knots[mid + 1]:
        if u < knots[mid]:
            high = mid
        else:
            low = mid
        mid = (low + high) // 2
    return mid


def basis(span, u, p, knots):
    """The ``p + 1`` nonzero basis functions of ``span`` at ``u``."""
    n = [0.0] * (p + 1)
    left = [0.0] * (p + 1)
    right = [0.0] * (p + 1)
    n[0] = 1.0
    for j in range(1, p + 1):
        left[j] = u - knots[span + 1 - j]
        right[j] = knots[span + j] - u
        saved = 0.0
        for r in range(j):
            temp = n[r] / (right[r + 1] + left[j - r])
            n[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        n[j] = saved
    return n


def spans(p, knots, n):
    """Knot indices ``i`` (``p .. n - 1``) of the spans of nonzero length."""
    return [i for i in range(p, n) if knots[i] < knots[i + 1]]


def curve_point_in_span(p, knots, cps, weights, span, u):
    """A curve point at ``u`` in ``span``; ``cps`` rows of floats."""
    nb = basis(span, u, p, knots)
    dim = len(cps[0])
    out = [0.0] * dim
    w = 0.0
    for k in range(p + 1):
        i = span - p + k
        c = nb[k] * weights[i] if weights is not None else nb[k]
        row = cps[i]
        for d in range(dim):
            out[d] += c * row[d]
        w += c
    if weights is not None:
        out = [x / w for x in out]
    return out


def surface_point_in_span(
    p, q, knots_u, knots_v, cps, weights, nu, span_u, span_v, u, v
):
    """A surface point at ``(u, v)`` in spans ``(span_u, span_v)``."""
    bu = basis(span_u, u, p, knots_u)
    bv = basis(span_v, v, q, knots_v)
    dim = len(cps[0])
    out = [0.0] * dim
    w = 0.0
    for ll in range(q + 1):
        temp = [0.0] * (dim + 1)
        j = span_v - q + ll
        for k in range(p + 1):
            i = span_u - p + k
            idx = j * nu + i
            c = bu[k] * weights[idx] if weights is not None else bu[k]
            row = cps[idx]
            for d in range(dim):
                temp[d] += c * row[d]
            temp[dim] += c
        for d in range(dim):
            out[d] += bv[ll] * temp[d]
        w += bv[ll] * temp[dim]
    if weights is not None:
        out = [x / w for x in out]
    return out
