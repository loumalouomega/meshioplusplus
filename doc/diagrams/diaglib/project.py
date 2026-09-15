"""An orthographic projector for the 3-D cell figures.

The view is tikz-3dplot's ``\\tdplotsetmaincoords{70}{110}`` -- the same one
the old ``doc/cell_types.tex`` used -- so z points up, x comes down-left toward
the viewer and y runs right. The matrix is transcribed from tikz-3dplot's own
``\\tdplotsetmaincoords`` definition.
"""

import math

THETA = math.radians(70.0)
PHI = math.radians(110.0)

# Screen-x and screen-up components of the three axes.
_XX, _XY = math.cos(PHI), -math.cos(THETA) * math.sin(PHI)
_YX, _YY = math.sin(PHI), math.cos(THETA) * math.cos(PHI)
_ZX, _ZY = 0.0, math.sin(THETA)

#: Unit vector pointing from the scene toward the viewer (screen-x cross
#: screen-up), so ``dot(p, TOWARD_VIEWER)`` is larger for closer points.
TOWARD_VIEWER = (
    math.sin(PHI) * math.sin(THETA),
    -math.cos(PHI) * math.sin(THETA),
    math.cos(THETA),
)


def project(p, scale=1.0, origin=(0.0, 0.0)):
    """Project a 3-D point to SVG coordinates (y grows downward)."""
    x, y, z = p
    sx = _XX * x + _YX * y + _ZX * z
    sy = _XY * x + _YY * y + _ZY * z
    return origin[0] + scale * sx, origin[1] - scale * sy


def depth(p):
    """Larger is closer to the viewer."""
    return sum(a * b for a, b in zip(p, TOWARD_VIEWER))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def newell_normal(points):
    """Newell's normal of a polygon given as a list of 3-D points."""
    nx = ny = nz = 0.0
    n = len(points)
    for i in range(n):
        x0, y0, z0 = points[i]
        x1, y1, z1 = points[(i + 1) % n]
        nx += (y0 - y1) * (z0 + z1)
        ny += (z0 - z1) * (x0 + x1)
        nz += (x0 - x1) * (y0 + y1)
    return (nx, ny, nz)


def face_is_visible(points):
    """True when an outward-wound face turns toward the viewer."""
    return dot(newell_normal(points), TOWARD_VIEWER) > 1e-12


def centroid(points):
    n = float(len(points))
    return tuple(sum(p[i] for p in points) / n for i in range(len(points[0])))
