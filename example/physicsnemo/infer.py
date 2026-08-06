"""Inference on the manifest's test split: predictions back onto the mesh.

For every test case the trained MeshGraphNet predicts ``T`` from ``q``, the
prediction and its error are written back as ordinary ``point_data``
(``T_pred``/``T_error``) into ``predictions/<id>.vtu`` — readable by ParaView
or any meshio++ consumer — and the first case is rendered as a
truth/prediction/error panel into ``renders/prediction.png``.
"""

import argparse
import json
import os

import numpy as np
import torch
from train import FIELDS, TARGETS, build_model, stat_tensor

import meshioplusplus as mio
import meshioplusplus.physicsnemo as mpn


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=str, default="dataset_manifest.json")
    ap.add_argument("--checkpoint", type=str, default="checkpoint.pt")
    args = ap.parse_args()

    device = "cuda" if torch.cuda.is_available() else "cpu"
    checkpoint = torch.load(args.checkpoint, map_location=device, weights_only=False)
    model = build_model(device)
    model.load_state_dict(checkpoint["model"])
    model.eval()

    with open("node_stats.json", encoding="utf-8") as fh:
        node_stats = json.load(fh)
    with open("edge_stats.json", encoding="utf-8") as fh:
        edge_stats = json.load(fh)
    x_mean = stat_tensor(node_stats, FIELDS, "mean", device)
    x_std = stat_tensor(node_stats, FIELDS, "std", device)
    y_mean = stat_tensor(node_stats, TARGETS, "mean", device)
    y_std = stat_tensor(node_stats, TARGETS, "std", device)
    e_mean = torch.tensor(edge_stats["edge_mean"], dtype=torch.float32, device=device)
    e_std = torch.tensor(
        edge_stats["edge_std"], dtype=torch.float32, device=device
    ).clamp_min(1e-8)

    manifest = mio.DatasetManifest.load(args.manifest)
    os.makedirs("predictions", exist_ok=True)
    os.makedirs("renders", exist_ok=True)

    # the feature-drift guard: the schema recorded at training time must match
    trained = checkpoint["schema"]

    first_panel = None
    errors = []
    for entry in manifest.entries(split="test"):
        _, mesh = entry.time_series()[0]
        sample = mpn.graph_sample(
            mesh, fields=FIELDS, target_fields=TARGETS, regions=False
        )
        assert sample.schema["x_columns"] == trained["x_columns"], "feature drift"
        arrays = sample.arrays
        with torch.no_grad():
            x = (torch.from_numpy(arrays["x"]).to(device) - x_mean) / x_std
            edge_attr = (
                torch.from_numpy(arrays["edge_attr"]).to(device) - e_mean
            ) / e_std
            from torch_geometric.data import Data

            graph = Data(
                x=x,
                edge_index=torch.from_numpy(arrays["edge_index"]).to(device),
                edge_attr=edge_attr,
            )
            pred_norm = model(x, edge_attr, graph)
            pred = (pred_norm * y_std + y_mean).cpu().numpy()[:, 0]
        truth = arrays["y"][:, 0]
        mesh.point_data["T_pred"] = pred.astype(np.float64)
        mesh.point_data["T_error"] = np.abs(pred - truth).astype(np.float64)
        mio.write(os.path.join("predictions", f"{entry.id}.vtu"), mesh)
        rmse = float(np.sqrt(np.mean((pred - truth) ** 2)))
        errors.append((entry.id, rmse))
        if first_panel is None:
            first_panel = (entry.id, mesh, truth, pred)

    for entry_id, rmse in errors:
        print(f"  {entry_id}: RMSE {rmse:.4f}")
    print(
        f"mean RMSE over {len(errors)} test cases: "
        f"{np.mean([e[1] for e in errors]):.4f}"
    )

    # render truth / prediction / error for the first test case
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    entry_id, mesh, truth, pred = first_panel
    tri = mtri.Triangulation(mesh.points[:, 0], mesh.points[:, 1], mesh.cells[0].data)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    for ax, values, title, cmap in (
        (axes[0], truth, "T (truth)", "viridis"),
        (axes[1], pred, "T (MeshGraphNet)", "viridis"),
        (axes[2], np.abs(pred - truth), "|error|", "magma"),
    ):
        m = ax.tripcolor(tri, values, cmap=cmap, shading="gouraud")
        ax.triplot(tri, color="white", linewidth=0.15, alpha=0.5)
        ax.set_aspect("equal")
        ax.set_title(title)
        fig.colorbar(m, ax=ax, shrink=0.85)
    fig.suptitle(f"test case '{entry_id}'")
    fig.savefig("renders/prediction.png", dpi=150)
    print("wrote predictions/*.vtu and renders/prediction.png")


if __name__ == "__main__":
    main()
