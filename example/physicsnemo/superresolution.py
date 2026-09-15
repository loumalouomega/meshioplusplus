"""Superresolution worked example: coarse grid -> SRResNet -> fine mesh.

The grid counterpart of this directory's MeshGraphNet example, and the thing
[roadmap section 1](../../doc/roadmap.md)'s grid bullet asked for, executed for
real: synthetic 3-D solution fields -> a coarse/fine grid pair sampled off each
mesh -> `physicsnemo.models.srrn.SRResNet` -> the prediction scattered back onto
the mesh, and scored against the trilinear baseline it has to beat.

**Self-supervised, which is the ordinary shape of a superresolution dataset.**
Each case is one high-resolution solve; the coarse input is made by sampling it
onto a coarser lattice, so no paired second run is needed and the manifest
carries no `Target`. (An entry *can* carry one -- see
[paired cases](../../doc/datasets.md#paired-cases) -- for when the coarse data
is a genuinely separate simulation.)

The field is deliberately multi-scale: a smooth large-scale part plus a
higher-wavenumber part whose amplitude varies per case. A model that only
learns the smooth part scores well on RMSE and badly on the power spectrum,
which is exactly the distinction `spectrum_rel_l2` exists to make.

Run:

    python superresolution.py --cases 60 --epochs 150

Everything a run writes lands in one run directory; this script adds the
figures. Needs torch + nvidia-physicsnemo (no `torch_geometric` -- a
convolutional model never touches it).
"""

import argparse
import json
import os

import numpy as np

import meshioplusplus as mio
import meshioplusplus.physicsnemo as mpn

FIELDS = ["T"]
TARGETS = ["T"]
MESH_RESOLUTION = 12  # cells per side of the unit-cube tet mesh
COARSE = (7, 7, 7)  # coarse grid, in cells -> 8^3 sample points
SCALE = 2  # -> 16^3 fine sample points


def case_mesh(rng):
    """A tetrahedral unit cube carrying a multi-scale scalar field."""
    mesh = mio.convert_cells(
        mio.grid((MESH_RESOLUTION,) * 3, spacing=(1.0 / MESH_RESOLUTION,) * 3),
        mode="simplexify",
    )
    p = mesh.points
    # a smooth part every method gets right, plus a fine part only a trained
    # model recovers -- the amplitude of the second is what varies per case
    low = np.sin(2.0 * np.pi * p[:, 0]) * np.cos(2.0 * np.pi * p[:, 1])
    amplitude = 0.3 + 0.4 * rng.random()
    phase = 2.0 * np.pi * rng.random()
    high = (
        amplitude
        * np.sin(6.0 * np.pi * p[:, 0] + phase)
        * np.sin(5.0 * np.pi * p[:, 2])
    )
    mesh.point_data["T"] = low + high + 0.25 * p[:, 2]
    return mesh


def build_dataset(directory, cases, seed):
    os.makedirs(directory, exist_ok=True)
    rng = np.random.default_rng(seed)
    manifest = mio.DatasetManifest(base_dir=directory)
    for case in range(cases):
        mio.write(os.path.join(directory, f"case_{case:04d}.vtu"), case_mesh(rng))
        split = (
            "train"
            if case < 0.7 * cases
            else ("valid" if case < 0.85 * cases else "test")
        )
        manifest.add(f"case_{case:04d}.vtu", id=f"case_{case:04d}", split=split)
    path = os.path.join(directory, "manifest.json")
    manifest.save(path)
    return path


def baseline_and_model_scores(checkpoint, manifest_path, device):
    """Score the trained model and the trilinear baseline on the test split.

    The baseline is the honest comparison: `resample_grid` is what you get for
    free, so a model that does not beat it has not earned its training time.
    """
    import torch
    from physicsnemo.core.module import Module

    from meshioplusplus._grid_transfer import power_spectrum, resample_grid
    from meshioplusplus.physicsnemo._train import card_path, normalizers_from_card

    card = json.loads(open(card_path(str(checkpoint))).read())
    norms = normalizers_from_card(card)
    model = Module.from_checkpoint(str(checkpoint)).to(device).eval()
    shape = (-1, 1, 1, 1)

    scores = {
        "model": {"rmse": [], "spectrum": []},
        "trilinear": {"rmse": [], "spectrum": []},
    }
    manifest = mio.DatasetManifest.load(manifest_path)
    for entry in manifest.entries(split="test"):
        _, mesh = entry.time_series()[0]
        sample = mpn.grid_sample_pair(mesh, **card["grid"])
        x, y = sample.arrays["x"], sample.arrays["y"]
        coarse = mio.GridSpec.from_dict(sample.schema["coarse"])
        fine = mio.GridSpec.from_dict(sample.schema["fine"])

        xt = (
            torch.from_numpy(
                (x - norms["x_mean"].reshape(shape)) / norms["x_std"].reshape(shape)
            )
            .float()
            .unsqueeze(0)
            .to(device)
        )
        with torch.no_grad():
            pred = model(xt)[0].cpu().numpy()
        pred = pred * norms["y_std"].reshape(shape) + norms["y_mean"].reshape(shape)
        base = resample_grid(x, coarse, fine)

        truth_spectrum = power_spectrum(y, fine).power
        denominator = float(np.sqrt(np.sum(truth_spectrum**2)))
        for name, values in (("model", pred), ("trilinear", base)):
            scores[name]["rmse"].append(float(np.sqrt(np.mean((values - y) ** 2))))
            difference = power_spectrum(values, fine).power - truth_spectrum
            scores[name]["spectrum"].append(
                float(np.sqrt(np.sum(difference**2)) / denominator)
            )
    return {
        name: {k: float(np.mean(v)) for k, v in metrics.items()}
        for name, metrics in scores.items()
    }


def figures(run_dir, checkpoint, manifest_path, device, out_dir="renders"):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import torch
    from physicsnemo.core.module import Module

    from meshioplusplus._grid_transfer import power_spectrum, resample_grid
    from meshioplusplus.physicsnemo._train import (
        card_path,
        normalizers_from_card,
        read_metrics,
    )

    os.makedirs(out_dir, exist_ok=True)
    rows = read_metrics(run_dir)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.semilogy(
        [r["epoch"] for r in rows], [r["train_loss"] for r in rows], label="train"
    )
    valid = [(r["epoch"], r["valid_loss"]) for r in rows if r["valid_loss"] is not None]
    if valid:
        ax.semilogy([e for e, _ in valid], [v for _, v in valid], label="valid")
    ax.set_xlabel("epoch")
    ax.set_ylabel("MSE (normalized)")
    ax.legend()
    ax.set_title("SRResNet on the synthetic multi-scale field")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "sr_loss_curve.png"), dpi=150)
    plt.close(fig)

    card = json.loads(open(card_path(str(checkpoint))).read())
    norms = normalizers_from_card(card)
    model = Module.from_checkpoint(str(checkpoint)).to(device).eval()
    shape = (-1, 1, 1, 1)
    manifest = mio.DatasetManifest.load(manifest_path)
    entry = list(manifest.entries(split="test"))[0]
    _, mesh = entry.time_series()[0]
    sample = mpn.grid_sample_pair(mesh, **card["grid"])
    x, y = sample.arrays["x"], sample.arrays["y"]
    coarse = mio.GridSpec.from_dict(sample.schema["coarse"])
    fine = mio.GridSpec.from_dict(sample.schema["fine"])
    xt = (
        torch.from_numpy(
            (x - norms["x_mean"].reshape(shape)) / norms["x_std"].reshape(shape)
        )
        .float()
        .unsqueeze(0)
        .to(device)
    )
    with torch.no_grad():
        pred = model(xt)[0].cpu().numpy()
    pred = pred * norms["y_std"].reshape(shape) + norms["y_mean"].reshape(shape)
    base = resample_grid(x, coarse, fine)

    mid = y.shape[1] // 2
    panels = [
        ("coarse input", x[0, x.shape[1] // 2]),
        ("trilinear", base[0, mid]),
        ("SRResNet", pred[0, mid]),
        ("truth", y[0, mid]),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.6))
    limits = dict(vmin=float(y.min()), vmax=float(y.max()), cmap="viridis")
    for axis, (title, panel) in zip(axes, panels):
        axis.imshow(panel, origin="lower", **limits)
        axis.set_title(title)
        axis.set_xticks([])
        axis.set_yticks([])
    fig.suptitle(f"{entry.id}: a z-slice through the field")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "sr_slices.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    # Truth first and thick, so the two others reading on top of it show
    # agreement rather than hiding it.
    styles = (
        ("truth", y, {"lw": 5, "alpha": 0.35, "color": "#333333"}),
        ("SRResNet", pred, {"lw": 1.8}),
        ("trilinear", base, {"lw": 1.8}),
    )
    peak = 0.0
    for label, values, style in styles:
        spectrum = power_spectrum(values, fine)
        keep = spectrum.counts > 0
        ax.loglog(
            spectrum.wavenumber[keep][1:],
            spectrum.power[keep][1:],
            label=label,
            **style,
        )
        peak = max(peak, float(spectrum.power[keep][1:].max()))
    # The field is band-limited, so the truth's own tail is numerical noise
    # around 1e-30; plotting it compresses every meaningful decade into a
    # sliver. Eight decades below the peak is well past anything physical.
    ax.set_ylim(bottom=peak * 1e-8, top=peak * 4)
    ax.set_xlabel("wavenumber (cycles per unit length)")
    ax.set_ylabel("power")
    ax.legend()
    ax.set_title("Where the trilinear baseline loses its small scales")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "sr_spectrum.png"), dpi=150)
    plt.close(fig)
    print(f"wrote {out_dir}/sr_loss_curve.png, sr_slices.png, sr_spectrum.png")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cases", type=int, default=60)
    ap.add_argument("--data-dir", type=str, default="sr_cases")
    ap.add_argument("--run-dir", type=str, default="runs/superresolution")
    ap.add_argument("--epochs", type=int, default=150)
    ap.add_argument("--batch-size", type=int, default=4)
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
        model_name="srresnet",
        scaling_factor=SCALE,
        conv_layer_size=32,
        resid_blocks=4,
        resolution=COARSE,
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.lr,
        seed=args.seed,
        device=args.device,
        checkpoint_every=50,
    )
    progress = mpn.run_training(spec)
    run_dir = spec.resolved_run_dir()
    checkpoint = os.path.join(run_dir, "checkpoints", "best.mdlus")

    device = "cuda" if args.device == "auto" else args.device
    try:
        import torch

        if not torch.cuda.is_available():
            device = "cpu"
    except Exception:
        device = "cpu"

    scores = baseline_and_model_scores(checkpoint, manifest_path, device)
    print("\ntest split, averaged:")
    print(f"  {'':10s} {'RMSE':>10} {'spectrum_rel_l2':>18}")
    for name in ("trilinear", "model"):
        print(
            f"  {name:10s} {scores[name]['rmse']:10.5f} "
            f"{scores[name]['spectrum']:18.5f}"
        )
    ratio = scores["trilinear"]["spectrum"] / max(scores["model"]["spectrum"], 1e-12)
    print(f"\n  the model's spectrum is {ratio:.1f}x closer to the truth's")
    with open(os.path.join(run_dir, "scores.json"), "w") as handle:
        json.dump(scores, handle, indent=2)

    figures(run_dir, checkpoint, manifest_path, device)
    print(
        f"\nbest checkpoint {progress['best_checkpoint']} "
        f"(epoch {progress['best_epoch']}, validation loss "
        f"{progress['best_valid_loss']:.3e})"
    )


if __name__ == "__main__":
    main()
