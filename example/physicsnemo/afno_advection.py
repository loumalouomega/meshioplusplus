"""Advection-diffusion worked example: (c0, u, v) -> AFNO -> c(T).

The second neural-operator example: an **AFNO**, which is 2-D only and
patches a *fixed* `(H, W)` image, so the spec's `Grid.Squeeze` is required
and the grid's sample shape must divide the patch. This script prints the
arithmetic that everyone trips on once: `Grid.Resolution [63, 63, 1]` gives
**64** sample points per axis (n cells have n + 1 corners), and `PatchSize
[8, 8]` divides 64.

Each case: a Gaussian blob `c0` at a random centre and width on the periodic
unit square, plus a constant velocity `(u, v)` drawn per case and carried as
two constant fields (`Fields: ["c0", "u", "v"]`). The truth is the exact
spectral solution of

    c_t + u c_x + v c_y = D lap c,   periodic,   at t = T

so there is no discretization error in the target at all. Two baselines: a
first-order **upwind** scheme with the true velocity (it knows the physics
and pays numerical diffusion for it), and **persistence** (`c(T) = c0`, the
floor any model must clear).

Run:

    python afno_advection.py --cases 200 --epochs 100

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

CELLS = 63  # -> 64 sample points per axis, which PatchSize [8, 8] divides
POINTS = CELLS + 1
PATCH = (8, 8)
DIFFUSION = 2e-3
T_FINAL = 0.4
FIELDS = ["c0", "u", "v"]
TARGETS = ["c"]


def lattice():
    """64 points per axis at spacing 1/64: x_i = i/64, so the FFT's period is
    the unit square exactly (the last point sits at 63/64, not at 1)."""
    return mio.convert_cells(
        mio.grid((CELLS, CELLS, 1), spacing=(1.0 / POINTS, 1.0 / POINTS, 0.05)),
        mode="simplexify",
    )


def plane_to_points(field):
    return np.concatenate([field.reshape(-1), field.reshape(-1)])


def points_to_plane(values):
    return np.asarray(values)[: POINTS**2].reshape(POINTS, POINTS)


def blob(rng):
    """A periodic Gaussian blob, indexed [j, i] = (y, x)."""
    x = np.arange(POINTS) / POINTS
    X, Y = np.meshgrid(x, x)  # X varies along axis 1 (i), Y along axis 0 (j)
    cx, cy = rng.random(2)
    width = 0.06 + 0.06 * rng.random()
    dx = (X - cx + 0.5) % 1.0 - 0.5
    dy = (Y - cy + 0.5) % 1.0 - 0.5
    return np.exp(-(dx**2 + dy**2) / (2 * width**2))


def spectral_solution(c0, u, v, t):
    """Exact: each mode advects with phase exp(-i k.u t) and decays with
    exp(-D |k|^2 t)."""
    k = 2.0 * np.pi * np.fft.fftfreq(POINTS, d=1.0 / POINTS)
    kx, ky = np.meshgrid(k, k)  # kx along axis 1 (x), ky along axis 0 (y)
    factor = np.exp(-1j * (kx * u + ky * v) * t - DIFFUSION * (kx**2 + ky**2) * t)
    return np.real(np.fft.ifft2(np.fft.fft2(c0) * factor))


def upwind(c0, u, v, t, cfl=0.4):
    """First-order upwind advection + explicit diffusion with the TRUE
    velocity: the baseline that knows the physics."""
    h = 1.0 / POINTS
    speed = max(abs(u), abs(v), 1e-12)
    dt = min(cfl * h / speed, 0.2 * h * h / DIFFUSION)
    steps = int(np.ceil(t / dt))
    dt = t / steps
    c = c0.copy()
    for _ in range(steps):
        # axis 1 is x, axis 0 is y; roll(+1) brings the upstream neighbour
        dcdx = (
            (c - np.roll(c, 1, axis=1)) / h
            if u >= 0
            else (np.roll(c, -1, axis=1) - c) / h
        )
        dcdy = (
            (c - np.roll(c, 1, axis=0)) / h
            if v >= 0
            else (np.roll(c, -1, axis=0) - c) / h
        )
        lap = (
            np.roll(c, 1, axis=1)
            + np.roll(c, -1, axis=1)
            + np.roll(c, 1, axis=0)
            + np.roll(c, -1, axis=0)
            - 4 * c
        ) / (h * h)
        c = c + dt * (-u * dcdx - v * dcdy + DIFFUSION * lap)
    return c


def build_dataset(directory, cases, seed):
    os.makedirs(directory, exist_ok=True)
    rng = np.random.default_rng(seed)
    manifest = mio.DatasetManifest(base_dir=directory)
    base = lattice()
    for case in range(cases):
        c0 = blob(rng)
        u, v = rng.uniform(-0.6, 0.6, size=2)
        c = spectral_solution(c0, u, v, T_FINAL)
        mesh = mio.Mesh(base.points, base.cells)
        mesh.point_data["c0"] = plane_to_points(c0)
        mesh.point_data["u"] = np.full(len(base.points), u)
        mesh.point_data["v"] = np.full(len(base.points), v)
        mesh.point_data["c"] = plane_to_points(c)
        mio.write(os.path.join(directory, f"case_{case:04d}.vtu"), mesh)
        manifest.add(
            f"case_{case:04d}.vtu",
            id=f"case_{case:04d}",
            split=split_for(case, cases),
            metadata={"u": float(u), "v": float(v)},
        )
    path = os.path.join(directory, "manifest.json")
    manifest.save(path)
    return path


def score(checkpoint, manifest_path):
    scores = {
        name: {"rmse": [], "rel_l2": []} for name in ("persistence", "upwind", "model")
    }
    example = None
    for entry in mio.DatasetManifest.load(manifest_path).entries(split="test"):
        _, mesh = entry.time_series()[0]
        c0 = points_to_plane(mesh.point_data["c0"])
        truth = points_to_plane(mesh.point_data["c"])
        u, v = entry.metadata["u"], entry.metadata["v"]
        base = upwind(c0, u, v, T_FINAL)
        predicted, _ = mpn.predict_mesh(checkpoint, mesh, label=entry.id)
        pred = points_to_plane(predicted.point_data["c_pred"])
        for name, values in (("persistence", c0), ("upwind", base), ("model", pred)):
            scores[name]["rmse"].append(rmse(values, truth))
            scores[name]["rel_l2"].append(relative_l2(values, truth))
        if example is None:
            example = (entry.id, c0, truth, pred, base)
    averaged = {
        name: {k: float(np.mean(v)) for k, v in metrics.items()}
        for name, metrics in scores.items()
    }
    return averaged, example


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases", type=int, default=200)
    ap.add_argument("--data-dir", type=str, default="afno_cases")
    ap.add_argument("--run-dir", type=str, default="runs/afno_advection")
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--lr", type=float, default=5e-4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--device", type=str, default="auto")
    ap.add_argument(
        "--skip-data", action="store_true", help="reuse an existing dataset"
    )
    args = ap.parse_args()

    print(
        f"Grid.Resolution [{CELLS}, {CELLS}, 1] gives {POINTS} sample points per axis "
        f"(n cells have n + 1 corners); PatchSize {list(PATCH)} divides {POINTS}: "
        f"{POINTS % PATCH[0] == 0 and POINTS % PATCH[1] == 0}"
    )
    manifest_path = os.path.join(args.data_dir, "manifest.json")
    if not args.skip_data or not os.path.isfile(manifest_path):
        manifest_path = build_dataset(args.data_dir, args.cases, args.seed)
    print(f"dataset: {manifest_path}")

    spec = mpn.default_spec(
        manifest_path,
        FIELDS,
        TARGETS,
        run_dir=args.run_dir,
        model_name="afno",
        resolution=(CELLS, CELLS, 1),
        squeeze=2,
        squeeze_index=0,
        patch_size=PATCH,
        embed_dim=256,
        depth=4,
        num_blocks=8,
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
        "afno_advection",
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
        run_dir,
        os.path.join(RENDERS_DIR, "afno_loss_curve.png"),
        "AFNO on advection-diffusion",
    )
    entry_id, c0, truth, pred, base = example
    panels(
        os.path.join(RENDERS_DIR, "afno_panels.png"),
        [
            ("c0 (input)", c0, False),
            ("upwind", base, False),
            ("AFNO", pred, False),
            ("truth", truth, False),
            ("|AFNO - truth|", np.abs(pred - truth), True),
        ],
        f"{entry_id}: concentration at t = {T_FINAL}",
        vmin=0.0,
        vmax=float(truth.max()),
    )
    print(f"wrote {RENDERS_DIR}/afno_loss_curve.png, afno_panels.png")
    print(
        f"\nbest checkpoint {checkpoint} (epoch {progress['best_epoch']}, "
        f"validation loss {progress['best_valid_loss']:.3e}, {elapsed:.1f} s)"
    )


if __name__ == "__main__":
    main()
