"""Darcy flow worked example: permeability field -> FNO -> pressure field.

The first of the three neural-operator examples that closed
[roadmap section 1](../../doc/roadmap.md)'s last bullet: a **2-D FNO** through
the thin-axis squeeze idiom. Each case is a random log-normal permeability
`K` on the unit square; the truth is the pressure `p` solving

    -div(K grad p) = 1,   p = 0 on the boundary

by a five-point finite-difference scheme (harmonic-mean interface
conductivity). `Fields: ["K"]`, `TargetFields: ["p"]`.

**The mesh is a thin 3-D lattice, and the squeeze is honest.** `grid((32, 32,
1))` has two z-planes, both carrying the same planar field, so `Grid.Squeeze:
2, SqueezeIndex: 0` keeps plane 0 and loses nothing. `Grid.Resolution [32,
32, 1]` gives 33x33 sample points, which is exactly the lattice's own point
set, so sampling is exact too.

The baseline is the same solve with the **effective medium** `K_eff =
exp(mean(log K))` -- what you get by ignoring the heterogeneity -- so the
number to beat is how much of the pressure field's structure comes from the
permeability's own.

Run:

    python fno_darcy.py --cases 200 --epochs 100

Needs torch + nvidia-physicsnemo (no torch_geometric).
"""

import argparse
import os
import time

import numpy as np
from _common import (
    RENDERS_DIR,
    environment_note,
    loss_curve,
    panels,
    pick_device,
    print_scores,
    relative_l2,
    rmse,
    split_for,
    write_scores,
)

import meshioplusplus as mio
import meshioplusplus.physicsnemo as mpn

N = 32  # cells per side -> N + 1 = 33 sample points per axis
FIELDS = ["K"]
TARGETS = ["p"]


def lattice():
    """The thin lattice every case lives on: two identical z-planes."""
    return mio.convert_cells(
        mio.grid((N, N, 1), spacing=(1.0 / N, 1.0 / N, 0.05)), mode="simplexify"
    )


def plane_to_points(field):
    """A (N+1, N+1) plane, indexed [j, i] (y, x), onto BOTH planes of the
    lattice's x-fastest point numbering."""
    return np.concatenate([field.reshape(-1), field.reshape(-1)])


def points_to_plane(values):
    return np.asarray(values)[: (N + 1) ** 2].reshape(N + 1, N + 1)


def log_normal_field(rng, correlation=0.15, sigma=1.0):
    """A Gaussian random field with a Gaussian covariance, by FFT filtering,
    exponentiated into a permeability."""
    n = N + 1
    white = rng.standard_normal((n, n))
    k = np.fft.fftfreq(n, d=1.0 / n)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    # a Gaussian covariance of length `correlation`: angular wavenumber 2*pi*k
    kernel = np.exp(-0.5 * (correlation * 2.0 * np.pi) ** 2 * (kx**2 + ky**2))
    smooth = np.real(np.fft.ifft2(np.fft.fft2(white) * kernel))
    smooth = (smooth - smooth.mean()) / smooth.std()
    return np.exp(sigma * smooth)


def solve_darcy(K):
    """-div(K grad p) = 1 on the unit square, p = 0 on the boundary, on the
    (N+1)x(N+1) node grid, harmonic-mean conductivities at the interfaces."""
    n = N + 1
    h = 1.0 / N
    interior = np.arange(1, n - 1)
    m = len(interior)
    index = -np.ones((n, n), dtype=np.int64)
    index[np.ix_(interior, interior)] = np.arange(m * m).reshape(m, m)
    A = np.zeros((m * m, m * m))
    b = np.full(m * m, h * h)

    def harmonic(a, c):
        return 2.0 * a * c / (a + c)

    for j in interior:
        for i in interior:
            row = index[j, i]
            kw = harmonic(K[j, i], K[j, i - 1])
            ke = harmonic(K[j, i], K[j, i + 1])
            ks = harmonic(K[j, i], K[j - 1, i])
            kn = harmonic(K[j, i], K[j + 1, i])
            A[row, row] = kw + ke + ks + kn
            for jj, ii, coeff in (
                (j, i - 1, kw),
                (j, i + 1, ke),
                (j - 1, i, ks),
                (j + 1, i, kn),
            ):
                col = index[jj, ii]
                if col >= 0:
                    A[row, col] = -coeff
    p = np.zeros((n, n))
    p[np.ix_(interior, interior)] = np.linalg.solve(A, b).reshape(m, m)
    return p


def build_dataset(directory, cases, seed):
    os.makedirs(directory, exist_ok=True)
    rng = np.random.default_rng(seed)
    manifest = mio.DatasetManifest(base_dir=directory)
    base = lattice()
    for case in range(cases):
        K = log_normal_field(rng)
        p = solve_darcy(K)
        mesh = mio.Mesh(base.points, base.cells)
        mesh.point_data["K"] = plane_to_points(K)
        mesh.point_data["p"] = plane_to_points(p)
        mio.write(os.path.join(directory, f"case_{case:04d}.vtu"), mesh)
        manifest.add(
            f"case_{case:04d}.vtu", id=f"case_{case:04d}", split=split_for(case, cases)
        )
    path = os.path.join(directory, "manifest.json")
    manifest.save(path)
    return path


def score(checkpoint, manifest_path):
    """The trained FNO (through the shipped `predict_mesh`) against the
    effective-medium baseline, on the test split."""
    scores = {
        "effective_medium": {"rmse": [], "rel_l2": []},
        "model": {"rmse": [], "rel_l2": []},
    }
    example = None
    for entry in mio.DatasetManifest.load(manifest_path).entries(split="test"):
        _, mesh = entry.time_series()[0]
        K = points_to_plane(mesh.point_data["K"])
        truth = points_to_plane(mesh.point_data["p"])
        base = solve_darcy(np.full_like(K, np.exp(np.mean(np.log(K)))))
        predicted, row = mpn.predict_mesh(checkpoint, mesh, label=entry.id)
        pred = points_to_plane(predicted.point_data["p_pred"])
        for name, values in (("effective_medium", base), ("model", pred)):
            scores[name]["rmse"].append(rmse(values, truth))
            scores[name]["rel_l2"].append(relative_l2(values, truth))
        if example is None:
            example = (entry.id, K, truth, pred, base)
    averaged = {
        name: {k: float(np.mean(v)) for k, v in metrics.items()}
        for name, metrics in scores.items()
    }
    return averaged, example


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases", type=int, default=200)
    ap.add_argument("--data-dir", type=str, default="fno_cases")
    ap.add_argument("--run-dir", type=str, default="runs/fno_darcy")
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", type=str, default="auto")
    ap.add_argument(
        "--skip-data", action="store_true", help="reuse an existing dataset"
    )
    args = ap.parse_args()

    manifest_path = os.path.join(args.data_dir, "manifest.json")
    if not args.skip_data or not os.path.isfile(manifest_path):
        manifest_path = build_dataset(args.data_dir, args.cases, args.seed)
    print(f"dataset: {manifest_path}")

    spec = mpn.default_spec(
        manifest_path,
        FIELDS,
        TARGETS,
        run_dir=args.run_dir,
        model_name="fno",
        resolution=(N, N, 1),
        squeeze=2,
        squeeze_index=0,
        num_fno_modes=12,
        latent_channels=32,
        num_fno_layers=4,
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.lr,
        seed=args.seed,
        device=args.device,
        checkpoint_every=50,
    )
    started = time.time()
    progress = mpn.run_training(spec)
    elapsed = time.time() - started
    run_dir = spec.resolved_run_dir()
    checkpoint = progress["best_checkpoint"]

    scores, example = score(checkpoint, manifest_path)
    print_scores(scores)
    out = write_scores(
        run_dir,
        "fno_darcy",
        scores,
        extra={
            "cases": args.cases,
            "epochs": progress["epoch"],
            "best_epoch": progress["best_epoch"],
            "best_valid_loss": progress["best_valid_loss"],
            "train_seconds": elapsed,
            "device": pick_device(args.device),
            **environment_note(),
        },
    )
    print(f"scores -> {out}")

    loss_curve(
        run_dir, os.path.join(RENDERS_DIR, "fno_loss_curve.png"), "FNO on Darcy flow"
    )
    entry_id, K, truth, pred, base = example
    lo, hi = float(truth.min()), float(truth.max())
    panels(
        os.path.join(RENDERS_DIR, "fno_panels.png"),
        [
            ("log K (input)", np.log(K), True),
            ("effective medium", base, False),
            ("FNO", pred, False),
            ("truth", truth, False),
            ("|FNO - truth|", np.abs(pred - truth), True),
        ],
        f"{entry_id}: pressure on the unit square",
        vmin=lo,
        vmax=hi,
    )
    print(f"wrote {RENDERS_DIR}/fno_loss_curve.png, fno_panels.png")
    print(
        f"\nbest checkpoint {checkpoint} (epoch {progress['best_epoch']}, "
        f"validation loss {progress['best_valid_loss']:.3e}, {elapsed:.1f} s)"
    )


if __name__ == "__main__":
    main()
