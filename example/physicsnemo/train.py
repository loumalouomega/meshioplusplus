"""Train PhysicsNeMo's MeshGraphNet on the manifest's train split.

The training path is the PyG one (`make_dataset` -> `torch_geometric`
DataLoader), because MeshGraphNet consumes PyG `Data` and PyG batches
variable-size graphs natively; the Gen-2 `Reader` is smoke-run once so the
Hydra-facing surface is exercised too. Normalization stats stream from the
manifest through `field_stats`/`edge_stats` and are written in the
`node_stats.json`/`edge_stats.json` convention PhysicsNeMo's own datapipes
use. Everything is seeded; run `infer.py` afterwards.
"""

import argparse
import json
import time

import torch

import meshioplusplus.physicsnemo as mpn

FIELDS = ["q_scaled"]
TARGETS = ["T"]


def stat_tensor(stats, names, suffix, device):
    values = []
    for name in names:
        values.extend(stats[f"{name}_{suffix}"])
    t = torch.tensor(values, dtype=torch.float32, device=device)
    return t if suffix == "mean" else t.clamp_min(1e-8)


def build_model(device):
    from physicsnemo.models.meshgraphnet import MeshGraphNet

    return MeshGraphNet(
        input_dim_nodes=len(FIELDS),
        input_dim_edges=4,
        output_dim=len(TARGETS),
        processor_size=8,
        hidden_dim_processor=64,
        hidden_dim_node_encoder=64,
        hidden_dim_edge_encoder=64,
        hidden_dim_node_decoder=64,
    ).to(device)


def run_epoch(model, loader, norms, device, optimizer=None):
    x_mean, x_std, y_mean, y_std, e_mean, e_std = norms
    total, count = 0.0, 0
    for batch in loader:
        batch = batch.to(device)
        x = (batch.x - x_mean) / x_std
        edge_attr = (batch.edge_attr - e_mean) / e_std
        target = (batch.y - y_mean) / y_std
        pred = model(x, edge_attr, batch)
        loss = torch.nn.functional.mse_loss(pred, target)
        if optimizer is not None:
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
        total += float(loss.detach()) * batch.num_graphs
        count += batch.num_graphs
    return total / max(count, 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=str, default="dataset_manifest.json")
    ap.add_argument("--epochs", type=int, default=100)
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(
        f"device: {device} ({torch.cuda.get_device_name(0) if device == 'cuda' else 'cpu'})"
    )

    # normalization stats: streamed from the manifest, stored in the
    # node_stats.json / edge_stats.json convention PhysicsNeMo datapipes use
    node_stats = mpn.field_stats(args.manifest, split="train", fields=FIELDS + TARGETS)
    edge_stats = mpn.edge_stats(args.manifest, split="train")
    with open("node_stats.json", "w", encoding="utf-8") as fh:
        json.dump(node_stats, fh, indent=2)
    with open("edge_stats.json", "w", encoding="utf-8") as fh:
        json.dump(edge_stats, fh, indent=2)

    kwargs = dict(fields=FIELDS, target_fields=TARGETS, regions=False)
    train_ds = mpn.make_dataset(args.manifest, split="train", **kwargs)
    valid_ds = mpn.make_dataset(args.manifest, split="valid", **kwargs)
    print(f"train: {len(train_ds)} graphs, valid: {len(valid_ds)} graphs")

    # smoke the Gen-2 Reader once, so the Hydra-facing surface is exercised
    reader = mpn.make_reader(args.manifest, split="valid", **kwargs)
    sample = reader._load_sample(0)
    assert set(sample) >= {"pos", "x", "y", "edge_index", "edge_attr"}
    print(f"reader: {len(reader)} samples, sample 0 keys: {sorted(sample)}")

    from torch_geometric.loader import DataLoader

    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True)
    valid_loader = DataLoader(valid_ds, batch_size=args.batch_size)

    norms = (
        stat_tensor(node_stats, FIELDS, "mean", device),
        stat_tensor(node_stats, FIELDS, "std", device),
        stat_tensor(node_stats, TARGETS, "mean", device),
        stat_tensor(node_stats, TARGETS, "std", device),
        torch.tensor(edge_stats["edge_mean"], dtype=torch.float32, device=device),
        torch.tensor(
            edge_stats["edge_std"], dtype=torch.float32, device=device
        ).clamp_min(1e-8),
    )

    model = build_model(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    history = []
    start = time.time()
    for epoch in range(args.epochs):
        model.train()
        train_loss = run_epoch(model, train_loader, norms, device, optimizer)
        model.eval()
        with torch.no_grad():
            valid_loss = run_epoch(model, valid_loader, norms, device)
        history.append((train_loss, valid_loss))
        if epoch % 10 == 0 or epoch == args.epochs - 1:
            print(f"epoch {epoch:4d}  train {train_loss:.3e}  valid {valid_loss:.3e}")
    print(f"trained {args.epochs} epochs in {time.time() - start:.1f} s")

    torch.save(
        {"model": model.state_dict(), "schema": train_ds.schema},
        "checkpoint.pt",
    )

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.semilogy([h[0] for h in history], label="train")
    ax.semilogy([h[1] for h in history], label="valid")
    ax.set_xlabel("epoch")
    ax.set_ylabel("MSE (normalized)")
    ax.legend()
    ax.set_title("MeshGraphNet on the synthetic heat dataset")
    fig.tight_layout()
    fig.savefig("renders/loss_curve.png", dpi=150)
    print(
        "wrote checkpoint.pt, node_stats.json, edge_stats.json, "
        "renders/loss_curve.png"
    )


if __name__ == "__main__":
    import os

    os.makedirs("renders", exist_ok=True)
    main()
