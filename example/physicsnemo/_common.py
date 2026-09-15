"""What the three neural-operator worked examples share.

`fno_darcy.py`, `afno_advection.py` and `deeponet_beam.py` mirror
`superresolution.py`'s sections -- build a dataset, train through the shipped
`run_training`, score the trained model against a baseline it has to beat,
render -- and this module holds the four things they would otherwise copy:
the split assignment, the device pick, the `scores.json` writer and the
table printer. The scores land in the run directory (gitignored, like every
run's own files) AND in the committed `scores/` directory beside these
scripts, so the numbers quoted in the docs are the numbers a reader can
diff.
"""

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCORES_DIR = os.path.join(HERE, "scores")
RENDERS_DIR = os.path.join(HERE, "renders")


def split_for(case, cases):
    """70 / 15 / 15 train / valid / test, by case index."""
    if case < 0.7 * cases:
        return "train"
    if case < 0.85 * cases:
        return "valid"
    return "test"


def pick_device(arg):
    """`"auto"` -> cuda when torch sees one, else cpu; anything else verbatim."""
    if arg != "auto":
        return arg
    try:
        import torch

        return "cuda" if torch.cuda.is_available() else "cpu"
    except Exception:  # noqa: BLE001 -- no torch at all: the trainer will say so
        return "cpu"


def relative_l2(pred, truth):
    """`||pred - truth|| / ||truth||`, the scale-free companion of the RMSE."""
    pred = np.asarray(pred, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    denominator = float(np.sqrt(np.sum(truth**2)))
    return (
        float(np.sqrt(np.sum((pred - truth) ** 2)) / denominator)
        if denominator
        else 0.0
    )


def rmse(pred, truth):
    return float(np.sqrt(np.mean((np.asarray(pred) - np.asarray(truth)) ** 2)))


def write_scores(run_dir, name, scores, extra=None):
    """`scores.json` into the run directory and `scores/<name>.json` here."""
    document = {"scores": scores}
    if extra:
        document.update(extra)
    os.makedirs(SCORES_DIR, exist_ok=True)
    for path in (
        os.path.join(run_dir, "scores.json"),
        os.path.join(SCORES_DIR, f"{name}.json"),
    ):
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(document, handle, indent=2)
            handle.write("\n")
    return os.path.join(SCORES_DIR, f"{name}.json")


def print_scores(scores, columns=("rmse", "rel_l2"), model="model"):
    """A small table, the model's row last, with the improvement ratio."""
    width = max(len(name) for name in scores) + 2
    print("\ntest split, averaged:")
    print("  " + "".ljust(width) + "".join(f"{c:>14}" for c in columns))
    for name, metrics in scores.items():
        print(
            "  " + name.ljust(width) + "".join(f"{metrics[c]:14.5f}" for c in columns)
        )
    for name, metrics in scores.items():
        if name == model:
            continue
        ratio = metrics[columns[0]] / max(scores[model][columns[0]], 1e-12)
        print(f"\n  {model} is {ratio:.1f}x better than {name} on {columns[0]}")


def loss_curve(run_dir, path, title):
    """The train/valid curve read back from `metrics.jsonl`."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from meshioplusplus.physicsnemo._train import read_metrics

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
    ax.set_title(title)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def panels(path, items, title, *, cmap="viridis", vmin=None, vmax=None):
    """A row of images sharing one colour range, plus an |error| panel."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(items), figsize=(3.4 * len(items), 3.4))
    if len(items) == 1:
        axes = [axes]
    for axis, (label, image, own) in zip(axes, items):
        limits = {} if own else {"vmin": vmin, "vmax": vmax, "cmap": cmap}
        if own:
            limits = {"cmap": "magma"}
        im = axis.imshow(image, origin="lower", **limits)
        axis.set_title(label)
        axis.set_xticks([])
        axis.set_yticks([])
        if own:
            fig.colorbar(im, ax=axis, fraction=0.046)
    fig.suptitle(title)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def environment_note():
    """One line saying what ran this, for the scores file."""
    note = {}
    try:
        import torch

        note["torch"] = torch.__version__
        if torch.cuda.is_available():
            note["device"] = torch.cuda.get_device_name(0)
    except Exception:  # noqa: BLE001
        pass
    try:
        import physicsnemo

        note["physicsnemo"] = physicsnemo.__version__
    except Exception:  # noqa: BLE001
        pass
    try:
        # the compiled core's own version: what actually ran, not whichever
        # dist-info an interpreter happens to resolve `meshioplusplus` to
        from meshioplusplus._core import __version__

        note["meshioplusplus"] = __version__
    except Exception:  # noqa: BLE001
        pass
    return note
