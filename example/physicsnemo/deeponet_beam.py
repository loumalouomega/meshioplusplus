"""Cantilever worked example: (Load, Modulus, PoissonRatio) -> DeepONet -> w.

The third neural-operator example, and the one with a different *shape* of
problem: **parameters in, field out**. ONE cantilever mesh (`grid((20, 4,
4))`, simplexified; L = 1, h = 0.2) serves every case; the cases differ only
in their `Metadata` -- the tip load `P`, Young's modulus `E` and Poisson's
ratio `nu`, drawn log-uniformly -- and that is what `TrainSpec`'s `Operator`
block reads: `Parameters: ["Load", "Modulus", "PoissonRatio"]`, with the trunk
over every mesh point. No `Fields` at all.

The truth is the Euler-Bernoulli deflection plus a Timoshenko shear
correction,

    w(x) = P x^2 (3L - x) / (6 E I)  +  P x / (k G A),   G = E / (2 (1 + nu))

so `w = (P / E) * [f1(x) + f2(x) (1 + nu)]`: the map from the parameter
vector to the field is **nonlinear** (P/E, and a product with nu). The
baseline is a per-node least-squares fit on `[1, P, E, nu]` over the training
cases -- exact if the map were linear, so the whole of the DeepONet's margin
over it IS the nonlinearity.

Run:

    python deeponet_beam.py --cases 200 --epochs 150

Needs torch + nvidia-physicsnemo with the experimental DeepONet (2.2.0
carries it; `mpn.has_deeponet()` says).
"""

import argparse
import os
import time

import numpy as np
from _common import (
    RENDERS_DIR,
    environment_note,
    loss_curve,
    pick_device,
    print_scores,
    relative_l2,
    rmse,
    split_for,
    write_scores,
)

import meshioplusplus as mio
import meshioplusplus.physicsnemo as mpn

LENGTH = 1.0
HEIGHT = 0.2
CELLS = (20, 4, 4)
SHEAR_FACTOR = 5.0 / 6.0
PARAMETERS = ["Load", "Modulus", "PoissonRatio"]
TARGETS = ["w"]


def beam():
    return mio.convert_cells(
        mio.grid(
            CELLS, spacing=(LENGTH / CELLS[0], HEIGHT / CELLS[1], HEIGHT / CELLS[2])
        ),
        mode="simplexify",
    )


def deflection(x, load, modulus, nu):
    inertia = HEIGHT**4 / 12.0
    area = HEIGHT * HEIGHT
    shear_modulus = modulus / (2.0 * (1.0 + nu))
    bending = load * x**2 * (3.0 * LENGTH - x) / (6.0 * modulus * inertia)
    shear = load * x / (SHEAR_FACTOR * shear_modulus * area)
    return bending + shear


def draw_parameters(rng):
    return {
        "Load": float(np.exp(rng.uniform(np.log(500.0), np.log(5000.0)))),
        "Modulus": float(np.exp(rng.uniform(np.log(5e10), np.log(2.5e11)))),
        "PoissonRatio": float(rng.uniform(0.2, 0.4)),
    }


def build_dataset(directory, cases, seed):
    os.makedirs(directory, exist_ok=True)
    rng = np.random.default_rng(seed)
    manifest = mio.DatasetManifest(base_dir=directory)
    base = beam()
    x = base.points[:, 0]
    for case in range(cases):
        params = draw_parameters(rng)
        mesh = mio.Mesh(base.points, base.cells)
        mesh.point_data["w"] = deflection(
            x, params["Load"], params["Modulus"], params["PoissonRatio"]
        )
        mio.write(os.path.join(directory, f"case_{case:04d}.vtu"), mesh)
        manifest.add(
            f"case_{case:04d}.vtu",
            id=f"case_{case:04d}",
            split=split_for(case, cases),
            metadata=params,
        )
    path = os.path.join(directory, "manifest.json")
    manifest.save(path)
    return path


def design_matrix(entries):
    return np.array(
        [
            [1.0, e.metadata["Load"], e.metadata["Modulus"], e.metadata["PoissonRatio"]]
            for e in entries
        ]
    )


def score(checkpoint, manifest_path, output_dir):
    """The DeepONet (through the shipped manifest `predict`, whose rows take
    each entry's parameters from its Metadata) against the per-node linear
    least-squares baseline fitted on the training cases."""
    manifest = mio.DatasetManifest.load(manifest_path)
    train = list(manifest.entries(split="train"))
    test = list(manifest.entries(split="test"))
    # the baseline: one least-squares fit per node on [1, P, E, nu]
    W = np.stack([e.time_series()[0][1].point_data["w"] for e in train])  # (cases, N)
    coefficients, *_ = np.linalg.lstsq(design_matrix(train), W, rcond=None)
    baseline = design_matrix(test) @ coefficients  # (test cases, N)

    rows = mpn.predict(checkpoint, manifest_path, split="test", output_dir=output_dir)
    by_id = {row["entry_id"]: row for row in rows}
    scores = {
        "least_squares": {"rmse": [], "rel_l2": []},
        "model": {"rmse": [], "rel_l2": []},
    }
    example = None
    for k, entry in enumerate(test):
        truth = entry.time_series()[0][1].point_data["w"]
        predicted = mio.read(by_id[entry.id]["output_path"]).point_data["w_pred"]
        for name, values in (("least_squares", baseline[k]), ("model", predicted)):
            scores[name]["rmse"].append(rmse(values, truth))
            scores[name]["rel_l2"].append(relative_l2(values, truth))
        if example is None:
            example = (entry, truth, predicted, baseline[k])
    averaged = {
        name: {k: float(np.mean(v)) for k, v in metrics.items()}
        for name, metrics in scores.items()
    }
    return averaged, example


def tip_figure(path, example):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    entry, truth, predicted, baseline = example
    mesh = entry.time_series()[0][1]
    p = mesh.points
    top = np.flatnonzero(
        (np.abs(p[:, 1] - HEIGHT) < 1e-9) & (np.abs(p[:, 2] - HEIGHT) < 1e-9)
    )
    order = top[np.argsort(p[top, 0])]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(
        p[order, 0],
        truth[order] * 1e3,
        lw=5,
        alpha=0.35,
        color="#333333",
        label="truth",
    )
    ax.plot(p[order, 0], predicted[order] * 1e3, lw=1.8, label="DeepONet")
    ax.plot(
        p[order, 0],
        baseline[order] * 1e3,
        lw=1.8,
        ls="--",
        label="least squares on [1, P, E, nu]",
    )
    ax.set_xlabel("x along the beam")
    ax.set_ylabel("deflection w (mm)")
    m = entry.metadata
    ax.set_title(
        f"{entry.id}: P = {m['Load']:.0f} N, E = {m['Modulus'] / 1e9:.0f} GPa, nu = {m['PoissonRatio']:.2f}"
    )
    ax.legend()
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases", type=int, default=200)
    ap.add_argument("--data-dir", type=str, default="deeponet_cases")
    ap.add_argument("--run-dir", type=str, default="runs/deeponet_beam")
    ap.add_argument("--epochs", type=int, default=150)
    ap.add_argument("--batch-size", type=int, default=16)
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
        [],  # a DeepONet's inputs are its parameters, not a data array
        TARGETS,
        run_dir=args.run_dir,
        model_name="deeponet",
        parameters=tuple(PARAMETERS),
        trunk="points",
        width=64,
        branch_layers=4,
        branch_layer_size=128,
        trunk_layers=4,
        trunk_layer_size=128,
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

    scores, example = score(
        checkpoint, manifest_path, os.path.join(run_dir, "predictions")
    )
    print_scores(scores)
    out = write_scores(
        run_dir,
        "deeponet_beam",
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
        os.path.join(RENDERS_DIR, "deeponet_loss_curve.png"),
        "DeepONet on the cantilever",
    )
    tip_figure(os.path.join(RENDERS_DIR, "deeponet_deflection.png"), example)
    print(f"wrote {RENDERS_DIR}/deeponet_loss_curve.png, deeponet_deflection.png")
    print(
        f"\nbest checkpoint {checkpoint} (epoch {progress['best_epoch']}, "
        f"validation loss {progress['best_valid_loss']:.3e}, {elapsed:.1f} s)"
    )


if __name__ == "__main__":
    main()
