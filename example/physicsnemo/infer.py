"""Inference on the manifest's test split: predictions back onto the mesh.

A thin wrapper over ``meshioplusplus.physicsnemo.predict``, which loads the
checkpoint together with its **model card** (the record of what each channel
means and how it was normalized), guards against feature drift, and writes
``T_pred``/``T_error`` back as ordinary ``point_data`` into
``predictions/<id>.vtu`` — readable by ParaView or any meshio++ consumer.
This script adds only the truth/prediction/error panel for the first case.
"""

import argparse
import os

import numpy as np

import meshioplusplus as mio
import meshioplusplus.physicsnemo as mpn


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=str, default="dataset_manifest.json")
    ap.add_argument("--run-dir", type=str, default="runs/example")
    ap.add_argument("--checkpoint", type=str, default=None)
    args = ap.parse_args()

    checkpoint = args.checkpoint or os.path.join(
        args.run_dir, "checkpoints", "best.mdlus"
    )
    rows = mpn.predict(
        checkpoint, args.manifest, split="test", output_dir="predictions"
    )
    for row in rows:
        print(f"  {row['entry_id']}: RMSE {row['rmse']:.4f}")
    print(
        f"mean RMSE over {len(rows)} test cases: "
        f"{np.mean([r['rmse'] for r in rows]):.4f}"
    )

    # render truth / prediction / error for the first test case
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    first = rows[0]
    mesh = mio.read(first["output_path"])
    truth = mesh.point_data["T"]
    pred = mesh.point_data["T_pred"]
    tri = mtri.Triangulation(mesh.points[:, 0], mesh.points[:, 1], mesh.cells[0].data)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for ax, values, title, cmap in (
        (axes[0], truth, "T (truth)", "viridis"),
        (axes[1], pred, "T (MeshGraphNet)", "viridis"),
        (axes[2], mesh.point_data["T_error"], "|error|", "magma"),
    ):
        m = ax.tripcolor(tri, values, cmap=cmap, shading="gouraud")
        ax.triplot(tri, color="white", linewidth=0.15, alpha=0.5)
        ax.set_aspect("equal")
        ax.set_title(title)
        fig.colorbar(m, ax=ax, shrink=0.85)
    fig.suptitle(f"test case '{first['entry_id']}'")
    os.makedirs("renders", exist_ok=True)
    fig.savefig("renders/prediction.png", dpi=150)
    print("wrote predictions/*.vtu and renders/prediction.png")


if __name__ == "__main__":
    main()
