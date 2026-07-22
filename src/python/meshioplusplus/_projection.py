"""Orthographic camera projection + painter's-algorithm face ordering.

Shared by the pure-Python SVG/TikZ 3D writer paths (and the logo generator).
This is the Python twin of ``src/cpp/include/meshioplusplus/detail/projection.hpp``
— KEEP THE ARITHMETIC IN SYNC, down to expression order: the TikZ test suite
asserts the C++ core and this fallback produce byte-identical output, which
requires identical floating-point rounding.

The camera is parameterized by azimuth/elevation/roll in degrees: the view
direction (scene toward camera) is
``w = (cos el * cos az, cos el * sin az, sin el)``, screen-right is
``u = normalize(cross((0,0,1), w))`` (y-axis fallback when looking straight
along z), screen-up is ``v = cross(w, u)``, and ``roll`` rotates ``(u, v)``
about ``w``. The defaults (45, atan(1/sqrt(2))) give ``w`` proportional to
``(1,1,1)`` — the classic CAD isometric view. Projection is orthographic;
``depth = p . w`` orders faces back-to-front (larger depth = closer).
"""

from __future__ import annotations

import math

import numpy as np

ISO_AZIMUTH = 45.0
ISO_ELEVATION = 35.264389682754654  # atan(1/sqrt(2))

# Drawable surface types: name -> (corner count, is_line). Higher-order
# types are corner-linearized (their corner nodes come first).
_DRAWABLE = {
    "line": (2, True),
    "triangle": (3, False),
    "triangle6": (3, False),
    "quad": (4, False),
    "quad8": (4, False),
    "quad9": (4, False),
}


def camera_basis(azimuth, elevation, roll):
    """The orthonormal camera basis ``(u, v, w)`` as three float triples."""
    deg = 3.141592653589793 / 180.0
    az = azimuth * deg
    el = elevation * deg
    caz = math.cos(az)
    saz = math.sin(az)
    cel = math.cos(el)
    sel = math.sin(el)
    w0 = cel * caz
    w1 = cel * saz
    w2 = sel

    if abs(w2) > 1.0 - 1.0e-12:
        # Looking straight along z: use y as the up reference.
        n = math.sqrt(w2 * w2 + w0 * w0)
        u0, u1, u2 = w2 / n, 0.0, -w0 / n
    else:
        n = math.sqrt(w1 * w1 + w0 * w0)
        u0, u1, u2 = -w1 / n, w0 / n, 0.0
    v0 = w1 * u2 - w2 * u1
    v1 = w2 * u0 - w0 * u2
    v2 = w0 * u1 - w1 * u0

    if roll != 0.0:
        cr = math.cos(roll * deg)
        sr = math.sin(roll * deg)
        u0, u1, u2, v0, v1, v2 = (
            cr * u0 + sr * v0,
            cr * u1 + sr * v1,
            cr * u2 + sr * v2,
            -sr * u0 + cr * v0,
            -sr * u1 + cr * v1,
            -sr * u2 + cr * v2,
        )
    return (u0, u1, u2), (v0, v1, v2), (w0, w1, w2)


def project_surface(mesh, azimuth, elevation, roll):
    """Project a surface mesh's points and gather its faces back-to-front.

    Returns ``(x, y, faces)`` where ``x``/``y`` are the projected screen
    coordinates per point and ``faces`` is a list of
    ``(node_ids, is_line, source_cell)`` tuples sorted back-to-front by
    face-centroid depth (painter's algorithm; ties keep enumeration order,
    matching the C++ core's stable sort).

    ``source_cell`` is the index of the cell the face came from, counted
    block-major over every cell of every block -- including blocks skipped
    here -- which is the same convention ``"surface:parent_cell"`` uses, so
    the two indices compose. It rides on each face rather than in a parallel
    array so the sort below cannot desync it (twin of
    ``detail::ProjectedFace::mSourceCell``).
    """
    u, v, w = camera_basis(azimuth, elevation, roll)
    points = np.asarray(mesh.points)
    dim = points.shape[1] if points.ndim == 2 else 0
    n = len(points)
    zeros = np.zeros(n, dtype=np.float64)
    # Cast each column to float64 first (the C++ core reads every coordinate
    # as double before the arithmetic).
    px = points[:, 0].astype(np.float64) if dim > 0 else zeros
    py = points[:, 1].astype(np.float64) if dim > 1 else zeros
    pz = points[:, 2].astype(np.float64) if dim > 2 else zeros
    x = px * u[0] + py * u[1] + pz * u[2]
    y = px * v[0] + py * v[1] + pz * v[2]
    depth = px * w[0] + py * w[1] + pz * w[2]

    conns = []
    is_lines = []
    depths = []
    srcs = []
    # Advances over EVERY block, including the skipped ones (mirrors the C++
    # core hoisting global_cell_base above its `continue`).
    base = 0
    for block in mesh.cells:
        block_base = base
        base += len(block.data)
        drawable = _DRAWABLE.get(block.type)
        if drawable is None:
            continue
        num_corners, is_line = drawable
        conn = np.asarray(block.data, dtype=np.int64)[:, :num_corners]
        # Face depth = mean of the corner depths, summed in node order
        # (mirrors the C++ accumulation order).
        d = depth[conn[:, 0]]
        for k in range(1, num_corners):
            d = d + depth[conn[:, k]]
        d = d / float(num_corners)
        conns.append(conn)
        is_lines.append(np.full(len(conn), is_line))
        depths.append(d)
        srcs.append(block_base + np.arange(len(conn), dtype=np.int64))

    faces = []
    if conns:
        all_depth = np.concatenate(depths)
        all_line = np.concatenate(is_lines)
        all_src = np.concatenate(srcs)
        flat = [row for conn in conns for row in conn]
        order = np.argsort(all_depth, kind="stable")
        faces = [(flat[i], bool(all_line[i]), int(all_src[i])) for i in order]
    return x, y, faces
