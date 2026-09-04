"""Train PhysicsNeMo's MeshGraphNet on the manifest's train split.

A thin wrapper over the shipped trainer: this script builds a
:class:`~meshioplusplus.physicsnemo.TrainSpec` from its arguments and calls
``meshioplusplus.physicsnemo.run_training``, which owns the training loop,
the normalization stats, the checkpoints and their model cards. The same
entry point backs ``python -m meshioplusplus.physicsnemo.train --spec`` and
the dataset dashboard's *Start training* button, so the example cannot drift
from what the library actually does.

Everything a run writes lands in one run directory (default ``runs/example``);
this script adds only the loss-curve figure, read back from the run's
``metrics.jsonl``. Run ``infer.py`` afterwards.
"""

import argparse
import os

import meshioplusplus.physicsnemo as mpn
from meshioplusplus.physicsnemo._train import read_metrics

FIELDS = ["q_scaled"]
TARGETS = ["T"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=str, default="dataset_manifest.json")
    ap.add_argument("--run-dir", type=str, default="runs/example")
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    spec = mpn.default_spec(
        args.manifest,
        FIELDS,
        TARGETS,
        run_dir=args.run_dir,
        epochs=args.epochs,
        batch_size=args.batch_size,
        learning_rate=args.lr,
        seed=args.seed,
        checkpoint_every=25,
    )
    progress = mpn.run_training(spec)
    run_dir = spec.resolved_run_dir()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = read_metrics(run_dir)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.semilogy(
        [r["epoch"] for r in rows], [r["train_loss"] for r in rows], label="train"
    )
    ax.semilogy(
        [r["epoch"] for r in rows if r["valid_loss"] is not None],
        [r["valid_loss"] for r in rows if r["valid_loss"] is not None],
        label="valid",
    )
    ax.set_xlabel("epoch")
    ax.set_ylabel("MSE (normalized)")
    ax.legend()
    ax.set_title("MeshGraphNet on the synthetic heat dataset")
    fig.tight_layout()
    os.makedirs("renders", exist_ok=True)
    fig.savefig("renders/loss_curve.png", dpi=150)
    print(
        f"best checkpoint {progress['best_checkpoint']} "
        f"(epoch {progress['best_epoch']}, validation loss {progress['best_valid_loss']:.3e}); "
        "wrote renders/loss_curve.png"
    )


if __name__ == "__main__":
    main()
