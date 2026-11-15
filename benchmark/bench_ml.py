#  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
# meshio++ — MIT License (see LICENSE). Main authors: Vicente Mataix Ferrandiz
"""Timings of the physics-ML Python layers and the converters (roadmap 5.1.2).

``bench_ops`` times the C++ operations; the layers a training loop actually
calls sit on top of them in Python and had no number: the lattice sampler
(``sample_grid``), the point-budget selector (``select_points``), the
proximity graph (``proximity_graph``) and a dataset write (a manifest of mesh
files, saved and loaded back). The ParaView and Blender converters had none
either: their payload layers (the arrays a ``vtkUnstructuredGrid`` and a
``bpy.types.Mesh`` are built from) import neither VTK nor Blender, so they are
timed here without the host application. This script builds a tetrahedral cube of
``n`` cells a side with a smooth point field and times each, best of
``--repeats`` after a warm-up, in milliseconds.

Usage:  python bench_ml.py [--n 30] [--repeats 5]
"""

from __future__ import annotations

import argparse
import os
import tempfile
import time

import numpy as np

import meshioplusplus as mio


def _cube(n):
    """The Kuhn tetrahedralisation of the unit cube, ``6 * n**3`` tets."""
    ax = np.arange(n + 1) / n
    k, j, i = np.meshgrid(ax, ax, ax, indexing="ij")
    pts = np.stack([i.ravel(), j.ravel(), k.ravel()], axis=1)
    idx = np.arange((n + 1) ** 3).reshape(n + 1, n + 1, n + 1)  # [k, j, i]
    base = idx[:-1, :-1, :-1].ravel()
    di, dj, dk = 1, n + 1, (n + 1) ** 2
    corner = [
        (0, di, di + dj, di + dj + dk),
        (0, di, di + dk, di + dj + dk),
        (0, dj, di + dj, di + dj + dk),
        (0, dj, dj + dk, di + dj + dk),
        (0, dk, di + dk, di + dj + dk),
        (0, dk, dj + dk, di + dj + dk),
    ]
    tets = np.concatenate(
        [np.stack([base + c for c in cs], axis=1) for cs in corner], axis=0
    )
    p = pts[tets]
    flip = np.linalg.det(p[:, 1:] - p[:, :1]) < 0
    tets[flip, 2], tets[flip, 3] = tets[flip, 3].copy(), tets[flip, 2].copy()
    mesh = mio.Mesh(pts, [("tetra", tets)])
    mesh.point_data["u"] = pts[:, 0] ** 2 + pts[:, 1] * pts[:, 2]
    return mesh


def _best(fn, repeats):
    fn()  # warm-up
    best = float("inf")
    for _ in range(repeats):
        t = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t)
    return best * 1e3


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--n", type=int, default=30, help="cells a side of the cube")
    ap.add_argument("--repeats", type=int, default=5)
    args = ap.parse_args()

    mesh = _cube(args.n)
    npts, ncells = len(mesh.points), len(mesh.cells[0].data)
    print(f"# tetra cube n={args.n}: {npts} points, {ncells} cells")
    print("layer,case,ms")

    def row(layer, case, fn):
        print(f"{layer},{case},{_best(fn, args.repeats):.2f}")

    # The lattice sampler: the coarse and the fine grid of a pair.
    for res in (32, 64):
        spec = mio.GridSpec.from_mesh(mesh, resolution=(res, res, res))
        row("sample_grid", f"{res}^3", lambda spec=spec: mio.sample_grid(mesh, spec))

    # The point-budget selector, all three methods.
    count = min(2048, npts // 4)
    for method in ("farthest", "grid", "random"):
        row(
            "select_points",
            f"{method} {count}",
            lambda method=method: mio.select_points(mesh, count, method=method),
        )

    # The proximity graph, over the mesh's points.
    radius = 2.5 / args.n
    row(
        "proximity_graph",
        f"radius {radius:.3f}",
        lambda: mio.proximity_graph(mesh, method="radius", radius=radius),
    )
    row(
        "proximity_graph",
        "knn 8",
        lambda: mio.proximity_graph(mesh, method="knn", max_neighbors=8),
    )

    # The converters' payload layers, without VTK or Blender: the mesh to the
    # arrays the host object is filled from, and (Blender) back.
    from meshioplusplus._blender import (
        _attribute_arrays,
        _mesh_from_blender_arrays,
        _to_blender_payload,
    )
    from meshioplusplus._interop import _to_vtk_payload

    row("paraview", "mesh -> VTK payload", lambda: _to_vtk_payload(mesh))
    row("blender", "mesh -> payload", lambda: _to_blender_payload(mesh))
    surface = mio.extract_surface(mesh)
    payload = _to_blender_payload(surface)
    row(
        "blender",
        "surface payload -> mesh",
        lambda: _mesh_from_blender_arrays(
            payload.vertices,
            payload.loop_vertices,
            payload.loop_starts,
            payload.loop_totals,
            edge_vertices=payload.edge_vertices,
            point_attributes=_attribute_arrays(
                payload.point_attributes, len(payload.vertices)
            ),
            face_attributes=_attribute_arrays(
                payload.face_attributes, payload.num_polygons
            ),
            custom_properties=payload.custom_properties,
        ),
    )

    # A dataset write: eight mesh files in a manifest, saved and loaded back.
    with tempfile.TemporaryDirectory() as d:
        for k in range(8):
            mio.write(os.path.join(d, f"case_{k}.vtu"), mesh)

        def dataset():
            manifest = mio.DatasetManifest(base_dir=d)
            for k in range(8):
                manifest.add(f"case_{k}.vtu", id=f"case_{k}", validate_source=False)
            path = os.path.join(d, "manifest.json")
            manifest.save(path)
            mio.DatasetManifest.load(path)

        row("dataset", "8 entries, save + load", dataset)


if __name__ == "__main__":
    main()
