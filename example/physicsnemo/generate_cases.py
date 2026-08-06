"""Generate the synthetic steady-heat dataset for the PhysicsNeMo example.

~200 cases of a manufactured Poisson pair on a jittered triangle mesh of the
unit square: the target field ``T`` is a Gaussian bump and the input field
``q = -laplacian(T)`` is its exact analytic source, so ``q -> T`` is a real
(if easy) PDE mapping for MeshGraphNet to learn. Every case varies the bump
centre/width, the node jitter and a small affine transform, all seeded.

Everything geometric goes through meshio++ operations (``convert_cells`` to
simplexify the quad grid, ``transform`` for the per-case affine), and the raw
cases land in ``cases_raw/`` — run the ``preprocess.json`` settings pipeline
next (see README.md).
"""

import argparse
import os

import numpy as np

import meshioplusplus as mio

BASE_RESOLUTION = 16  # cells per side of the unit-square grid


def base_mesh(n):
    """A triangle mesh of the unit square: quad grid -> simplexify."""
    xs = np.linspace(0.0, 1.0, n + 1)
    xg, yg = np.meshgrid(xs, xs, indexing="xy")
    points = np.column_stack([xg.ravel(), yg.ravel(), np.zeros((n + 1) * (n + 1))])
    quads = []
    for j in range(n):
        for i in range(n):
            v = j * (n + 1) + i
            quads.append([v, v + 1, v + n + 2, v + n + 1])
    mesh = mio.Mesh(points, [("quad", np.array(quads, dtype=np.int64))])
    return mio.convert_cells(mesh, mode="simplexify")


def manufactured_pair(points, cx, cy, s):
    """T = exp(-r^2 / (2 s^2)) and its exact 2-D source q = -laplacian(T)."""
    r2 = (points[:, 0] - cx) ** 2 + (points[:, 1] - cy) ** 2
    temperature = np.exp(-r2 / (2.0 * s * s))
    source = temperature * (2.0 / s**2 - r2 / s**4)
    return temperature, source


def make_case(rng, n):
    mesh = base_mesh(n)
    # jitter the interior nodes (boundary stays put so the domain is stable)
    pts = mesh.points.copy()
    interior = (
        (pts[:, 0] > 1e-9)
        & (pts[:, 0] < 1 - 1e-9)
        & (pts[:, 1] > 1e-9)
        & (pts[:, 1] < 1 - 1e-9)
    )
    amplitude = 0.25 / n
    pts[interior, :2] += rng.uniform(-amplitude, amplitude, (interior.sum(), 2))
    mesh = mio.Mesh(pts, mesh.cells)
    # a small per-case affine, through the real transform operation
    angle = float(rng.uniform(-10.0, 10.0))
    scale = float(rng.uniform(0.9, 1.1))
    mesh = mio.transform(mesh, rotate=([0, 0, 1], angle), scale=scale)
    # the manufactured pair, evaluated on the final geometry
    cx, cy = rng.uniform(0.3, 0.7, 2) * scale
    s = float(rng.uniform(0.08, 0.18)) * scale
    temperature, source = manufactured_pair(mesh.points, float(cx), float(cy), s)
    mesh.point_data["T"] = temperature
    mesh.point_data["q"] = source
    return mesh, {"cx": float(cx), "cy": float(cy), "s": s, "angle": angle}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases", type=int, default=200)
    ap.add_argument("--resolution", type=int, default=BASE_RESOLUTION)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=str, default="cases_raw")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    os.makedirs(args.out, exist_ok=True)
    os.makedirs("cases", exist_ok=True)  # the preprocess pipeline's output dir
    width = max(4, len(str(args.cases - 1)))
    params = []
    for i in range(args.cases):
        mesh, meta = make_case(rng, args.resolution)
        mio.write(os.path.join(args.out, f"case_{str(i).zfill(width)}.vtu"), mesh)
        params.append(meta)
    np.save(os.path.join(args.out, "params.npy"), np.array(params))
    print(f"wrote {args.cases} cases to {args.out}/")


if __name__ == "__main__":
    main()
