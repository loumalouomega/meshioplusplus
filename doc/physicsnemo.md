# PhysicsNeMo integration

::: tip New to PhysicsNeMo?
[**PhysicsNeMo basics**](./physicsnemo/overview.md) is a fourteen-page map of the framework itself — what a `Module` and a `.mdlus` checkpoint are, which of the 25 architecture families fits the shape of your data, how simulation output becomes batched tensors, and where meshio++ ends and the framework begins. This page assumes all of that.
:::

Feed meshes into [NVIDIA PhysicsNeMo](https://github.com/NVIDIA/physicsnemo) training pipelines — MeshGraphNet-style GNNs first of all — without writing the bespoke ingestion glue every project otherwise re-invents:

```python
import meshioplusplus as mio
import meshioplusplus.physicsnemo as mpn

manifest = mio.DatasetManifest.load("dataset_manifest.json")

sample = mpn.graph_sample(mio.read("case_0042.vtu"),
                          fields=["q"], target_fields=["T"])
stats = mpn.field_stats(manifest, split="train", fields=["q", "T"])

ds = mpn.make_dataset(manifest, split="train",
                      fields=["q"], target_fields=["T"])   # PyG Dataset
reader = mpn.make_reader(manifest, split="train")          # Gen-2 Reader
```

Everything here is pure Python over existing machinery — the [`feature_matrix` contract and `edge_index`](./ml), the [sequence machinery](./sequences), the [dataset manifests](./datasets) — and the C++/WASM/C/Fortran core is untouched. The adapter lives in its own subpackage (`meshioplusplus.physicsnemo`), deliberately **not** imported by `import meshioplusplus`: like the [MCP server](./mcp), its gated halves have their own Python floor and heavyweight dependencies, so the top-level package must work without them.

## Installation: deliberately no `[physicsnemo]` extra

```bash
pip install nvidia-physicsnemo     # the framework (imports as `physicsnemo`)
pip install torch_geometric        # for the PyG Dataset path
```

One more, found the hard way while executing the worked example: **MeshGraphNet's GNN layers import `torch_scatter` directly** (it is part of PhysicsNeMo's `gnns` extra, not its core dependencies), and `torch_scatter` ships prebuilt wheels only per torch/CUDA pair, lagging torch releases — at the time of writing there is none for torch 2.13, so the example pins torch 2.12:

```bash
pip install "torch==2.12.0" --index-url https://download.pytorch.org/whl/cu130
pip install torch-scatter -f https://data.pyg.org/whl/torch-2.12.0+cu130.html
```

The adapter itself never imports `torch_scatter` — only training a MeshGraphNet does.

The CuPy/torch packaging precedent, applied a third time: `nvidia-physicsnemo` hard-depends on `torch>=2.10`, whose default Linux wheel bundles CUDA at multiple GB — a `meshioplusplus[physicsnemo]` extra would transitively pin exactly the dependency the repo already [refuses to pin for `to_torch`](./ml#packaging-deliberately-no-torch-jax-extra). Install the framework directly; `mpn.has_physicsnemo()` / `mpn.has_torch_geometric()` answer availability without raising, and a missing install raises a named error quoting the command above.

Two floors to know about: meshio++ supports Python ≥ 3.8, PhysicsNeMo requires **≥ 3.11**. The mcp-SDK precedent applies — everything testable (`graph_sample`, the stats, the manifest iteration) is pure numpy and runs on every supported interpreter; only constructing the framework-facing `Reader`/`Dataset` objects needs the framework.

## Reconnaissance (verified 2026-08-06, physicsnemo 2.1.1)

Written down before the adapter was built, per the roadmap's own rule; these are the facts the design stands on. Package: `nvidia-physicsnemo` on PyPI (there is no `physicsnemo` package), Apache-2.0, pure-Python wheel, quarterly minor releases with real API churn — pin `nvidia-physicsnemo>=2.1,<2.3` in training projects.

- **DGL is gone.** The DGL→PyTorch Geometric migration completed across 25.08–25.11; `main` has zero DGL references. PyG is the only graph backend, and `torch_geometric` is an *extra* of PhysicsNeMo (`gnns`), not a core dependency. Nothing here targets DGL.
- **Two datapipe generations coexist.** Gen-1 (`physicsnemo.datapipes.gnn`) is plain `torch.utils.data.Dataset` subclasses yielding PyG `Data` objects (`.pos`, `.x`, `.y`, `.edge_index`, `.edge_attr`) — this is what the MeshGraphNet examples actually train on today. Gen-2 (`physicsnemo.datapipes`) is a GPU-centric `Reader → Dataset → DataLoader` pipeline over `TensorDict`, whose stated user extension points are **Readers and Transforms, not dataset subclasses**; the `Reader` ABC is two abstract methods (`_load_sample(index) -> dict[str, torch.Tensor]`, `__len__`), and an optional `@register()` decorator makes a reader addressable from Hydra configs. The `Datapipe` base class carries only metadata and imposes no contract — the adapter does not target it.
- **The edge-feature convention is stable across examples**: with `row, col = edge_index`, the edge attribute is `cat((pos[row] - pos[col], ‖pos[row] - pos[col]‖))` — a 4-vector in 3D. Note the sign: source minus destination. `graph_sample` reproduces it exactly, pinned by a numeric test.
- **Normalization stats are plain JSON by convention**: Gen-1 datasets read CWD-relative `node_stats.json` / `edge_stats.json` with `{field}_mean` / `{field}_std` keys; Gen-2 has a `Normalize` transform taking means/stds. `field_stats`/`edge_stats` produce dicts in that key convention; meshio++ defines no new stats file format.
- **`physicsnemo.mesh.Mesh`** (new in 2.1.0) is **simplicial-only** — points, lines, triangles, tetrahedra, exactly one kind per `Mesh` (mixed dimensions are "separate Mesh objects", its own docs) — and its only file on-ramp is a PyVista-backed `from_pyvista`/`to_pyvista` pair; its `.pmsh` save format is self-declared unstable, and its `io` module is a hardcoded PEP-562 set with no plugin registry, so a bridge can only live meshio++-side. (v10.35.0 added [`pmsh`](./formats/pmsh.md) and [`zarr`](./formats/zarr.md) as ordinary meshio++ formats, written in pure numpy against those layouts and pinned against upstream in both directions — so a *file* on-ramp now exists without torch, while the in-memory bridge below stays the way to hand a live `Mesh` over.) v9.28.0 deferred that bridge; **v9.30.0 builds it anyway** ([`to_physicsnemo`/`from_physicsnemo`](#the-physicsnemomesh-bridge)) with the risk stated rather than avoided: `physicsnemo.mesh` is the framework's newest and least stable surface, so the bridge touches only the tensorclass constructor and public attributes, never `.pmsh`, and training projects should pin `nvidia-physicsnemo>=2.1,<2.3`.
- **Every mesh file enters PhysicsNeMo through PyVista/VTK today** (the `VTKReader` hard-requires pyvista and reads `.stl`/`.vtp`/`.vtu`/`.vtk` only). That is the slot this adapter fills: meshio++'s 42 readers, regions, and preprocessing operations, ending in the tensors the datapipes expect.

## The adapter: `meshioplusplus.physicsnemo`

The `mcp/` split, applied again: the subpackage's `__init__` is pure numpy/stdlib/meshioplusplus and holds **all** the behaviour; `_reader.py` is the only module importing `physicsnemo`, `_pyg.py` the only one importing `torch_geometric`, and both are reached through lazy factories.

### One graph's worth of tensors: `graph_sample`

```python
s = mpn.graph_sample(mesh, fields=["q"], target_fields=["T"])
s.arrays        # {"pos", "x", "y", "edge_index", "edge_attr"} — numpy
s.x_columns     # ('q',)            the recorded feature contract
s.y_columns     # ('T',)
s.schema        # JSON-ready: versions, columns, edge-feature rule
```

| key | shape / dtype | rule |
|---|---|---|
| `pos` | `(N, dim)` float32¹ | `mesh.points`, never folded into `x` (the MGN convention — coordinates enter through edge features) |
| `x` | `(N, Fx)` float32¹ | [`feature_matrix(mesh, "point", fields=fields, coords=False, regions=regions)`](./ml#the-feature-matrix-and-its-contract-feature_matrix) — fields in stated order, then `region:<name>` one-hots (the node-type slot) |
| `y` | `(N, Fy)` float32¹ | a second `feature_matrix` over `target_fields`; absent when `target_fields=None`. Same-step by default; [`target_mesh`/`target_delta` pair steps](#autoregressive-t-t1-pairing) |
| `edge_index` | `(2, E)` **int64** | [`edge_index(mesh, kind=, undirected=)`](./ml#graphs-for-gnns-edge_index) verbatim |
| `edge_attr` | `(E, dim+1)` float32¹ | `cat((pos[row] − pos[col], ‖·‖))` — the PhysicsNeMo convention above |
| `world_edge_index` / `world_edge_attr` | `(2, Ew)` int64 / `(Ew, dim+1)` float32¹ | present only with `world_edges`, and kept **separate** rather than concatenated |

¹ `float32=True` is the default (the training convention); `float32=False` keeps meshio++'s canonical float64. `edge_index` is int64 always.

The column contract is `feature_matrix`'s, carried through: store `s.schema` at training time, compare `x_columns` at inference time, and feature drift is a named assertion instead of a silently wrong prediction.

Two keyword arguments change where the edges come from, both dictionaries in [`proximity_graph`](./proximity_graphs)'s vocabulary. **`proximity`** replaces the mesh edges with a radius or k-nearest-neighbour neighbourhood — the particle case, where a state has positions and no connectivity, and where the interaction radius is what makes two points neighbours. **`world_edges`** adds a second edge set beside the mesh edges, for contact between surfaces that are near but not connected.

```python
s = mpn.graph_sample(cloud, fields=["v"], proximity={"method": "radius", "radius": 0.015})
s = mpn.graph_sample(mesh, fields=["u"], world_edges={"method": "knn", "max_neighbors": 8})
```

The world set stays in its own arrays because PyTorch Geometric increments any attribute whose name contains `index` when it batches, so two sets batch correctly for free while a concatenated one could not be split apart afterwards. A model wanting HybridMeshGraphNet's single edge set concatenates them itself, **mesh edges first** — that model splits its edge features positionally.

How the edges were built is recorded in the schema, and `graph_sample_version` moved to **3** to say so: a schema stored under an earlier release now compares unequal, which is the drift guard working. A checkpoint trained on a proximity graph must not be replayed on a mesh one.

### Curved-cell tessellation

A curved-cell dataset (`triangle6`/`quad8`/`quad9`/`tetra10`/`hexahedron27`) feeds a GNN a linear mesh without losing its curvature via [`tessellate`](./tessellation):

```python
time, s = mpn._read_sample(series, step, {"fields": ["q"], "tessellate": True})
```

`Graph.Tessellate` is the spec key (`true` for every default, or `{"Levels": 2, "Curved": true, "Fields": true}` to narrow it), applied inside `_read_sample` before `graph_sample` ever sees the mesh, scoped to `Kind: "node"` graphs with `Graph.Regions: false` (tessellation's synthetic points carry no region membership of their own, so both are refused by name rather than silently producing an incomplete result). The resulting `Tessellation` lives only for that one sample's lifetime, matching the streaming invariant. At prediction time `predict_mesh`/`predict_file`/`_attach` take a `tess=` parameter: passing the sample's own `Tessellation` scatters/aggregates the prediction back onto the *original* mesh's own points/cells, and `tess=None` — the default — is byte-identical to every path that predated tessellation support.

### Streaming over a dataset: `iter_samples`, `field_stats`, `edge_stats`

```python
for entry_id, time, sample in mpn.iter_samples(manifest, split="train",
                                               fields=["q"], target_fields=["T"]):
    ...   # one mesh alive at a time — the sequence streaming invariant

stats = mpn.field_stats(manifest, split="train", fields=["q", "T"])
# {"q_mean": [...], "q_std": [...], "T_mean": [...], "T_std": [...]}
estats = mpn.edge_stats(manifest, split="train")
# {"edge_attr_mean": [...], "edge_attr_std": [...]}

# The SAME edge options training will use, or the statistics normalize the wrong graph
estats = mpn.edge_stats(manifest, split="train", proximity={"method": "radius", "radius": 0.015})
```

`edge_stats` takes `proximity`/`world_edges` too, and **must be given the same ones training will use**: it rebuilds the edges to gather its statistics, and statistics gathered over a different graph normalize the wrong thing with nothing downstream to report the discrepancy. Both go through one shared builder so they cannot drift, and `world_edges` adds `world_edge_mean`/`world_edge_std` to the result.

All three walk a [`DatasetManifest`](./datasets) (or anything `DatasetManifest.load` accepts), resolve each entry through the sequence plan, and read one mesh at a time. The stats are per-component streaming moments in the `{field}_mean`/`{field}_std` key convention — write them with `json.dump` as `node_stats.json`/`edge_stats.json` and Gen-1 datapipes read them as-is; pass them to Gen-2's `Normalize` transform directly.

### Framework objects: `make_dataset`, `make_reader`

```python
ds = mpn.make_dataset(manifest, split="train", fields=["q"], target_fields=["T"])
# torch.utils.data.Dataset of torch_geometric.data.Data(pos=, x=, y=, edge_index=, edge_attr=)

reader = mpn.make_reader(manifest, split="train", fields=["q"])
# physicsnemo Gen-2 Reader: __len__ + _load_sample(i) -> dict[str, torch.Tensor]
```

`make_dataset` is the training path: PyG's `DataLoader` batches variable-size graphs natively, which is what MeshGraphNet consumes. `make_reader` is the Gen-2 path for Hydra-driven pipelines (the reader is `@register()`ed when the registry is importable — addressability is a bonus, not a pinned contract); note TensorDict has no ragged-graph batching convention, so batching across meshes of different sizes remains the PyG path's job. Both index the flat `(entry, step)` sequence resolved once from the manifest's plans, and read one mesh per `__getitem__`/`_load_sample`.

Both of v9.28.0's recorded deferrals shipped in v9.30.0 — the [t→t+1 pairing](#autoregressive-t-t1-pairing) and the [`physicsnemo.mesh` bridge](#the-physicsnemomesh-bridge) below.

### Autoregressive t→t+1 pairing

```python
# steady-state surrogate (v1's shape, still the default): y from the SAME step
ds = mpn.make_dataset(manifest, fields=["q"], target_fields=["T"])

# autoregressive: step k's inputs paired with step k+1's targets
ds = mpn.make_dataset(manifest, target_fields=["v"], target_offset=1)

# the MeshGraphNet increment convention: y = v_{k+1} − v_k
ds = mpn.make_dataset(manifest, target_fields=["v"],
                      target_offset=1, target_delta=True)
stats = mpn.field_stats(manifest, fields=["v"], delta=True)
# {"v_diff_mean": [...], "v_diff_std": [...]} — normalize the delta targets
```

`target_offset=n` pairs step k with step k+n across `iter_samples`, `make_dataset` and `make_reader`: each entry contributes `len(series) − n` samples (an entry too short to pair contributes nothing, with one warning naming it, emitted at index-build time), and the yielded `time` is the *input* step's. `target_delta=True` makes `y` the target's values minus the same step's own — normalize with `field_stats(delta=...)` (`delta=True` is lag 1; an int matches a larger offset), whose `{field}_diff_mean`/`{field}_diff_std` keys stay plain JSON like everything else here.

Two honesty notes. A paired sample costs **two reads** — `TimeSeries` caches nothing, so the streaming invariant becomes "at most two meshes alive"; that is the price of never holding a dataset in memory, stated rather than cached away. And a **remeshed series cannot pair**: the target step's row count must match the input's, and a mismatch is a named error, never a silent mis-pairing. At the `graph_sample` level the same machinery is `target_mesh=` (the iteration layer derives it from the series); the schema records `target_offset`/`target_delta`, which is why `graph_sample_version` bumped to 2 — a stored v1 schema now compares unequal, exactly the drift guard doing its job.

### The `physicsnemo.mesh` bridge

```python
pm = mpn.to_physicsnemo(mesh)              # physicsnemo.mesh.Mesh, CPU
pm = mpn.to_physicsnemo(mesh, manifold_dim=2)   # pick the surface, not the volume
mesh_again = mpn.from_physicsnemo(pm)
pm.to("cuda")                              # device moves are upstream's job
```

`physicsnemo.mesh.Mesh` holds exactly **one simplex kind** — the cells tensor's trailing dim names it — so `to_physicsnemo` selects one topological dimension (`"auto"` = the highest present), tessellates non-simplex cells there through the existing [`convert_cells`](./convert_cells) operations (linearize, then simplexify — `cell_data` rows replicate to the children natively), and **drops with a warning** everything the target cannot hold: blocks of other dimensions, regions (encode membership as `graph_sample` one-hots instead), non-numeric data. Points follow the upstream float32 convention (`float32=False` keeps float64), cells are int64, numeric `field_data` lands in `global_data`, and every buffer is freshly owned — upstream in-place edits invalidate that `Mesh`'s caches, so sharing memory with the source would be a spooky-action bug. `from_physicsnemo` inverts the mapping (vertex/line/triangle/tetra by trailing dim, `global_data` → `field_data`).

**The risk, stated**: this rides the framework's newest surface. The bridge therefore touches only the tensorclass constructor and public attributes; the `.pmsh` and Zarr layouts are handled separately by the [`pmsh`](./formats/pmsh.md) and [`zarr`](./formats/zarr.md) formats, whose parity with upstream is pinned by tests in both directions — its gated module is the only one importing `physicsnemo.mesh` (an import that pulls in NVIDIA Warp, ~1.5 s, which is why it stays lazy), and training projects should pin `nvidia-physicsnemo>=2.1,<2.3`. If upstream's `Mesh` stabilizes, an `io_meshio.py` contributed upstream (mirroring `io_pyvista`'s shape) is still the natural end state; today its `io` module has no plugin registry, so the bridge lives here.

## Grid samples

A convolutional model does not take a graph. `grid_sample_pair` is `graph_sample`'s counterpart for one: it samples a mesh onto a **coarse** lattice and its targets onto a **fine** one, which is the pair a superresolution model trains on.

```python
import meshioplusplus as mio
import meshioplusplus.physicsnemo as mpn

mesh   = mio.read("case_0042.vtu")
coarse = mio.GridSpec.from_mesh(mesh, resolution=(32, 32, 32))

sample = mpn.grid_sample_pair(mesh, coarse, scaling_factor=2,
                              fields=["T"], target_fields=["T"])
sample.arrays["x"]      # (C, D, H, W)      -- the model's input
sample.arrays["y"]      # (C, 2D, 2H, 2W)   -- its target
sample.schema["scaling_factor"]
```

**`scaling_factor` goes through `GridSpec.upscale_samples`, not `upscale`, and the difference is an off-by-one that a model turns into a shape error deep inside a loss.** A convolutional upsampler multiplies *sample* counts: `SRResNet(scaling_factor=2)` maps `(B, C, 5, 5, 5)` to `(B, C, 10, 10, 10)`. `upscale` multiplies *cells*, so on a 4×4×4-cell grid it yields 9×9×9 points — right for resampling, where every coarse point is then also a fine point, and wrong here. Both preserve the box exactly; only the sample counts differ. Give exactly one of `scaling_factor` and an explicit `fine` spec.

`target_mesh=None` means **self-supervised** — the same mesh supplying both sides — and that is the ordinary case, not a fallback. `iter_grid_samples` walks a manifest the way `iter_samples` does, honouring the same streaming invariant (one mesh alive per sample, two for a [paired entry](./datasets#paired-cases)), and checks each entry's pairing once at index-build time rather than at the first epoch that reaches it.

`squeeze=<world axis>` collapses a thin axis for a 2-D operator: an integer `squeeze_index` keeps that plane, and omitting it averages over the axis.

**`grid_stats` is deliberately separate from `field_stats`.** A grid's statistics include the *fill* wherever the lattice reaches outside the mesh, so they are a different number from the nodal ones; normalizing a grid with node stats would be silently wrong. It reports `x_mean`/`x_std`/`y_mean`/`y_std` per channel plus the mean `coverage`, because a dataset whose grids are mostly fill is a dataset a model will learn the fill from.

See [mesh and regular grids](./grids) for the transfer itself.

## Training and prediction

The adapter builds the tensors; `run_training` is the loop over them — MeshGraphNet through the PyG path, driven by a **training spec** and writing everything a run produces into one directory:

```python
spec = mpn.default_spec("dataset_manifest.json", ["q_scaled"], ["T"],
                        run_dir="runs/example", epochs=100)
progress = mpn.run_training(spec)          # or: python -m meshioplusplus.physicsnemo.train --spec spec.json
rows = mpn.predict(progress["best_checkpoint"], "dataset_manifest.json",
                   split="test", output_dir="predictions")
```

The spec is a **hand-editable settings document** — PascalCase keys, `"Version": 1`, strict unknown-key refusal, exactly like a [dataset manifest](./datasets) or a [pipeline](./pipeline) — because a training run is something you write down, review and re-run:

```jsonc
{
  "Version": 1,
  "Manifest": "dataset_manifest.json",   // required; relative to this file
  "RunDir": "runs/example",
  "Fields": ["q_scaled"],                // required: the input arrays
  "TargetFields": ["T"],                 // required: what the model predicts
  "TrainSplit": "train", "ValidSplit": "valid",
  "Epochs": 100, "BatchSize": 8, "LearningRate": 0.001, "Seed": 0,
  "CheckpointEvery": 10, "Device": "auto",
  "Model": { "Name": "meshgraphnet", "ProcessorSize": 8, "HiddenDim": 64, "Aggregation": "sum" },
  "Graph": { "Regions": false, "Kind": "node", "Undirected": true,
             "EdgeFeatures": true, "Float32": true,
             "TargetOffset": 0, "TargetDelta": false,     // graph_sample's own options
             "Proximity": null,                           // or {"Method": "radius", "Radius": 0.015}
             "Tessellate": null },                        // or true, or {"Levels": 2}
  "Read": {}, "Notes": null, "Tags": []
}
```

Everything the run *writes* is a machine artefact and is snake_case (the [`write_dataset`](./ml#dataset-export-write_dataset) convention), in the run directory:

| file | written by | holds |
|---|---|---|
| `spec.json` | the launcher | the spec the run was started from |
| `metrics.jsonl` | the trainer | one row per epoch: `epoch`, `train_loss`, `valid_loss`, `lr`, `elapsed` |
| `progress.json` | the trainer | the live record: epoch, best epoch and loss, ETA, device, `completed` |
| `node_stats.json`, `edge_stats.json` | the trainer | the normalization stats, in PhysicsNeMo's own key convention |
| `checkpoints/` | the trainer | `<Model>.0.<epoch>.mdlus` + `checkpoint.0.<epoch>.pt` (`physicsnemo.utils.checkpoint.save_checkpoint`, resumable), plus `best.mdlus` and `final.mdlus` |
| `predictions/` | `predict` | `<entry_id>.vtu` with the predictions written back |
| `log.txt`, `job.json` | the [job manager](./dashboard#the-companion-process) | only when the run was launched as a job |

**Every `.mdlus` carries a model card** (`<checkpoint>.card.json`): the field names, the `x_columns`/`y_columns` contract, and the input/output/edge normalization the model was trained under. This is the [Kratos PhysicsNeMo application](https://github.com/KratosMultiphysics/Kratos)'s convention, adopted for the same reason it exists there — a checkpoint that does not say what its channels mean is a checkpoint you can silently misuse, and writing a model's normalized output onto a physical field produces finite, plausible, completely wrong numbers. `predict` reads the card rather than being told again, and refuses by name when the columns it recomputes differ from the ones recorded (the feature-drift guard).

`Graph.Proximity` builds the edges from geometry instead of connectivity, in [`proximity_graph`](./proximity_graphs)'s vocabulary (`Method` `"radius"` or `"knn"`, plus `Radius`/`MaxNeighbors` and an optional periodic `BoxSize`). A malformed block, or the other method's key, is refused when the document is read rather than an epoch later; the block is written back only when set, so a spec that does not use it is unchanged. There is deliberately no `Graph.WorldEdges` and no bistride block — no shipped model family reads a second edge set or a hierarchy, and a key a run would silently ignore is exactly what the strict unknown-key refusal exists to prevent.

Two details worth knowing. The model is sized from the **recorded schema**, not from `len(Fields)`: with `Graph.Regions` on, the region one-hots widen `x`, and sizing from the field count alone would build the wrong first layer. And `SIGTERM` is honoured — the running epoch finishes, `final.mdlus` and its card are written, and the process exits 143 — which is what makes the dashboard's *Stop* button leave a usable checkpoint rather than a truncated one.

### Single-mesh inference

`predict` needs a manifest, a split and an entry. Once a model is trained, applying it to *one* mesh that was never catalogued anywhere should not require inventing a manifest for it, and three forms do that:

```python
row = mpn.predict_file(checkpoint, "part.vtu", "part_pred.vtu")   # a file
mesh, row = mpn.predict_mesh(checkpoint, mesh)                    # in memory
```

```bash
meshioplusplus predict runs/example/checkpoints/best.mdlus part.vtu part_pred.vtu
```

Everything the prediction needs comes from the checkpoint's own **card**: which family wrote it (so a grid checkpoint needs no different call), the sample options, the read options, the column contract and the normalization. Nothing consults a manifest — which is exactly why `predict` itself is now a loop over the same per-mesh body, and why the two cannot drift.

`time_step` picks a step of a multi-step input and `target_path` supplies the paired mesh a t→t+n or coarse/fine checkpoint compares against. **A mesh carrying no truth predicts anyway**: the `<column>_pred` arrays are written, no `<column>_error` arrays are, and `rmse`/`max_error` come back `None`. That rule matters more than it looks — with the target fields absent, `graph_sample` would otherwise quietly take `y` from the input's own step and report an "error" of a prediction against itself.

The MCP `predict_file` tool is the same function, and `predict` inherits the rule: an entry whose meshes lack the target field now predicts with a `None` rmse instead of raising.

A checkpoint whose card records `Graph.Tessellate` rebuilds the same `Tessellation` from the input mesh and predicts on it, then `scatter`s/`aggregate`s the result back onto the mesh's own points/cells before writing it out — so a tessellated checkpoint's `predict_mesh`/`predict_file` call looks identical to an ordinary one from the outside.

## Superresolution: the `srresnet` family

`TrainSpec` knows five model families — `meshgraphnet`, `srresnet`, `fno`, `afno` and `deeponet`, one table (`_FAMILIES` in `physicsnemo/_train.py`) that the spec parser, the trainer, prediction, the MCP `train_start` tool and the dashboard's launch form all read — and each reads its own block. `Model.Name: "srresnet"` trains `physicsnemo.models.srrn.SRResNet` on the coarse/fine grid pairs [`grid_sample_pair`](#grid-samples) produces:

```json
{
  "Manifest": "dataset_manifest.json",
  "RunDir": "runs/sr",
  "Fields": ["T"],
  "TargetFields": ["T"],
  "Epochs": 150,
  "Model": { "Name": "srresnet", "ScalingFactor": 2, "ConvLayerSize": 32, "ResidBlocks": 4 },
  "Grid":  { "Resolution": [7, 7, 7] }
}
```

**A hyperparameter meant for another family is refused, not ignored.** `HiddenDim` on an srresnet, `ScalingFactor` on a meshgraphnet or an fno, a `Grid` block on a graph model, a `Graph` or `Operator` block on a CNN, a `ScalingFactor` outside {2, 4, 8}, or `Grid.Squeeze` on an srresnet (which is `Conv3d` throughout) each raise by name — and the refusal loops over *every* block rather than "the other one", so a third block cannot slip past it. A silently dropped key is how a run ends up training a model nobody asked for. `Grid.Float32` is read from the family's own block (it used to be accepted and emitted but read from `Graph`, so `false` on an srresnet was silently ignored), and `run_training` loads the spec *before* asking for frameworks, so a grid run through the public API no longer demands `torch_geometric`.

`Grid` takes the same lattice vocabulary as [`GridSpec.from_mesh`](./grids.md) — exactly one of `Resolution` and `CellSize`, plus `Bounds`, `Padding`, `PaddingRelative`, `Extrapolate`, `FillValue` and `MaxCells`. The **fine** grid is the coarse one through `upscale_samples(ScalingFactor)`, which is what makes the target's shape equal the model's output shape; see [that method's note](./grids.md#pairing-a-coarse-grid-with-a-fine-one) for why `upscale` is the wrong one here.

**A grid run needs `torch` and `nvidia-physicsnemo` but not `torch_geometric`.** PyTorch Geometric exists to batch ragged graphs and a dense `(C, D, H, W)` tensor is not one, so requiring it would refuse a runnable job over a dependency whose prebuilt wheels routinely lag torch releases. The grid dataset is a plain `torch.utils.data.Dataset`.

### The card, and what it must record

Every `.mdlus` gets a `*.card.json` as before, and the grid card records three things the graph card has no equivalent for: the **layout** as a literal string (`"channels_first_zyx"`), and both grid specs (`coarse`, `fine`). The layout matters most — a checkpoint that does not say which axis is which produces finite, plausible, transposed numbers, and nothing downstream can tell.

`predict` reads the card and dispatches on it, so a grid checkpoint needs no special call. It writes `<field>_pred`, `<field>_true` and `<field>_error` into an ordinary `.vtu`, which is why the dashboard's prediction preview works on a superresolution run with no browser change. Each row carries `rmse`, `max_error`, `coverage` and **`spectrum_rel_l2`** — the relative L2 between the predicted and true azimuthally averaged power spectra.

That last number is the one to report. A pointwise error cannot distinguish a field with the right small-scale content from a plausible smoothed one, and on the worked example the gap is stark — over 60 cases of a multi-scale field, 150 epochs in 149 s on one GPU:

| | RMSE | `spectrum_rel_l2` |
|---|---|---|
| trilinear baseline | 0.1556 | 0.2493 |
| SRResNet | **0.0081** | **0.0027** |

Nineteen times better pointwise, and **ninety-three times** closer in the spectrum. The two ratios differ that much because they measure different things: the baseline gets the large scales roughly right and loses the small ones entirely, which is precisely what a pointwise error under-reports. See [`example/physicsnemo/superresolution.py`](https://github.com/loumalouomega/meshioplusplus/blob/main/example/physicsnemo/superresolution.py).

## Neural operators on a grid: the `fno` and `afno` families

`Model.Name: "fno"` and `"afno"` train `physicsnemo.models.fno.FNO` and `physicsnemo.models.afno.AFNO` on the same [grid samples](#grid-samples) the `srresnet` family reads — and they are the families [roadmap section 1](./roadmap)'s last bullet asked for, closing it. Both are **resolution-preserving**: the model answers on the grid it was given, so neither has a `ScalingFactor` and `grid_kwargs()` pairs the coarse grid with *itself* (`upscale_samples(1)` is the identity, which is what lets the coarse/fine machinery serve both shapes with no second code path).

```jsonc
{
  "Manifest": "fno_cases/manifest.json", "RunDir": "runs/fno",
  "Fields": ["K"], "TargetFields": ["p"], "Epochs": 100,
  "Model": { "Name": "fno", "NumFnoModes": 12, "LatentChannels": 32, "NumFnoLayers": 4 },
  "Grid":  { "Resolution": [32, 32, 1], "Squeeze": 2, "SqueezeIndex": 0 }
}
```

**`Grid.Squeeze` is what makes an operator 2-D.** A lattice always has at least two planes on every axis, so a planar problem lives on a *thin* 3-D grid and the squeeze collapses one world axis — `SqueezeIndex` keeps that plane, omitting it averages over the axis (the [thin-axis idiom](./grids#two-dimensional-operators)). For `fno` the squeeze is optional and its presence *is* the FNO's `dimension` (3-D without it); for `afno` it is **required**, since AFNO patches a fixed `(H, W)` image, and the family refuses a spec without one by name. The FNO's `SpectralPadding` maps to the model's own `padding` — the name avoids colliding with `Grid.Padding`, the lattice's bounding-box padding, which is a different thing.

**AFNO's shape is fixed at construction, and the arithmetic everyone trips on once is stated by the error.** `Grid.Resolution [n, ...]` gives `n + 1` sample points per axis (n cells have n + 1 corners), and `PatchSize` must divide that: `[63, 63, 1]` gives 64 points, which `[8, 8]` divides; `[64, 64, 1]` gives 65, which nothing does. The check runs **before torch is imported** — a non-dividing patch, a 3-D grid, or an `EmbedDim` not divisible by `NumBlocks` is refused naming the spec key, never surfacing from inside the model's patch embedding. At inference the card's recorded `inp_shape` is checked against the mesh's sample shape the same way, before the forward.

**The card records the squeeze contract, and every part of it is load-bearing.** A 2-D grid card carries `layout` (the layout of the axes that *remain* — `channels_first_yx` for a squeezed z, `channels_first_zx` for y, `channels_first_zy` for x; `grid_layout_after_squeeze(axis)` in `meshioplusplus._grid_transfer` is the single owner, and a constant would be wrong because which axes remain depends on which one went), `spatial_ndim` (2 or 3), `squeeze`/`squeeze_index`, `x_shape` (the sample shape after the squeeze, what AFNO is built from) and `expand_axis`/`expand_size` (how a plane is duplicated back over the thin axis when the prediction is written onto the mesh). The grid-sample schema gained `x_shape`/`y_shape`/`expand_size` for this, so **`GRID_SAMPLE_VERSION` is now 2** — a stored v1 schema compares unequal, which is the drift guard doing its job. The normalizers are shaped from `spatial_ndim` rather than hard-coded to five dimensions, and that fix is asserted on the *shape*: against a `(B, C, H, W)` batch a `(1, C, 1, 1, 1)` normalizer broadcasts along the batch axis silently whenever `C == 1`, so the arithmetic alone could not have discriminated.

`predict`, `predict_mesh` and `predict_file` need no new call — the card says which family wrote the checkpoint — and write `<field>_pred`/`_true`/`_error` onto the mesh's own points, expanding the 2-D answer over every plane of the thin axis first. `spectrum_rel_l2` is reported as **`None`** for a squeezed card: a 2-D power spectrum is a documented follow-up, and computing the 3-D one over duplicated planes would be a number that means nothing. Neither family needs `torch_geometric`, and neither accepts `Augmentation` (a grid sampled on a fixed lattice has its own coverage changed by rotating the mesh under it — the `srresnet` rule).

Two worked examples, both executed on a GPU, in [`example/physicsnemo/`](https://github.com/loumalouomega/meshioplusplus/tree/main/example/physicsnemo): [`fno_darcy.py`](https://github.com/loumalouomega/meshioplusplus/blob/main/example/physicsnemo/fno_darcy.py) (a log-normal permeability field in, the Darcy pressure out, against the effective-medium solve) and [`afno_advection.py`](https://github.com/loumalouomega/meshioplusplus/blob/main/example/physicsnemo/afno_advection.py) (an initial blob and a velocity in, the exact spectral advection–diffusion solution out, against first-order upwind and persistence). Over 200 Darcy cases, 100 epochs in 103.1 s:

| | RMSE | relative L2 |
|---|---|---|
| effective-medium solve | 1.39e-2 | 0.373 |
| FNO | **2.02e-3** | **0.0535** |

and over 200 advection cases, 100 epochs in 255.7 s:

| | RMSE | relative L2 |
|---|---|---|
| persistence (`c(T) = c0`) | 1.51e-1 | 1.03 |
| first-order upwind, true velocity | **1.11e-2** | **0.0785** |
| AFNO | 1.85e-2 | 0.130 |

The two operators land on opposite sides of their stronger baseline, and both results are reported as run. The FNO beats the effective-medium solve by about 7 times in RMSE and relative L2: a constant permeability cannot place the pressure maximum where the low-permeability pockets trap it, and the operator learns exactly that. The AFNO beats persistence by about 8 times but does **not** beat first-order upwind — its RMSE is about 1.7 times upwind's. That baseline is strong on purpose: it is handed the true velocity and integrates the transport equation, whereas the AFNO has to infer the displacement from the `u`/`v` channels. Two things say the gap is not the family's ceiling, and neither was tuned away here. The run had not converged, its best validation epoch being 98 of 100. And the error map shows the 8x8 patch seams of the patch embedding along the blob's edge, where most of the residual sits — a smaller `PatchSize` or a longer run is the obvious next experiment, left to the reader rather than folded into the committed numbers.

## Parameters in, field out: the `deeponet` family

A great many engineering problems have no *field* as input at all — a load, a modulus, an inlet speed, a handful of numbers per case — and want a field out, and that is the shape people reach for a neural operator for by mistake. `Model.Name: "deeponet"` trains the installed PhysicsNeMo's experimental `DeepONet` (`physicsnemo.experimental.models.xdeeponet`, pinned to 2.2's keyword-only constructor; `mpn.has_deeponet()` says whether the framework has it, quietly): an MLP **branch** over the case's parameters, an MLP **trunk** over the mesh's own points, combined into the field at every trunk point. It reads its own block, `Operator`, and **no `Fields` at all** — its inputs come from each entry's [`Metadata`](./datasets):

```jsonc
{
  "Manifest": "deeponet_cases/manifest.json", "RunDir": "runs/deeponet",
  "TargetFields": ["w"], "Epochs": 150,
  "Model": { "Name": "deeponet", "Width": 64, "BranchLayers": 4, "TrunkLayers": 4 },
  "Operator": { "Parameters": ["Load", "Modulus", "PoissonRatio"], "Trunk": "points" }
}
```

| `Operator` key | meaning |
|---|---|
| `Parameters` | required — the ordered `Metadata` keys that become the branch input, on **every** entry. A scalar is one column under its own name; a list of numbers expands to `name_0`…`name_{k-1}` (the `feature_matrix`/pandas suffix rule, one rule repo-wide); a bool, a string or anything nested is refused by name; a missing key names the entry |
| `Trunk` | `"points"` (every mesh point, the default) or `"budget"` (a token budget through [`select_points`](./point_budgets), in ascending index order) |
| `TrunkCount`, `TrunkMethod`, `TrunkSeed` | the budget's size (required with `"budget"`, refused with `"points"`), `select_points`' method (`farthest`/`grid`/`random`) and seed |
| `Float32` | the sample dtype, as for the other blocks |

`Model` takes `BranchLayers`/`BranchLayerSize`, `TrunkLayers`/`TrunkLayerSize`, `Width` (the branch–trunk product dimension), `DecoderType`/`DecoderWidth`/`DecoderLayers`/`DecoderActivationFn` and `ActivationFn`. **`DecoderType` accepts only `"mlp"`**: the installed `DeepONet` refuses `"conv"` for an MLP branch itself ("pass a SpatialBranch as branch1") and `"temporal_projection"` needs an `output_window`, so both are refused by name here rather than from inside the constructor. Where a PascalCase key is shared with another family but the constructor default differs (`ActivationFn` is prelu/gelu/silu across srresnet/fno/deeponet, `DecoderLayers` is 1 for FNO and 2 for DeepONet), the spec holds a family-prefixed field, because `default_spec` re-validates through the document and a shared field would hand one family the other's default.

**Fixed geometry in v1, checked and named.** A DeepONet's trunk is *shared across the batch* — the model is called as `model(params_batch, trunk)` with one `(T, 3)` trunk broadcast over every case — so every mesh must present the same points. `iter_operator_samples`/`operator_stats`/`make_operator_dataset` read each entry's first step once at index-build time and refuse a point count or a parameter shape that differs from the first entry's, naming both: `entry 'b' yields 1234 point rows but entry 'a' yields 1000 -- a DeepONet trunk is shared across the batch, so every mesh must have the same points; use meshgraphnet for varying geometry`. The streaming invariant holds regardless: one mesh alive per sample.

The pure data path is `operator_sample(mesh, metadata, parameter_names=..., target_fields=..., trunk=...)` → an `OperatorSample` of `params` `(p,)`, `trunk` `(T, 3)` and `y` `(T, k)` (through `feature_matrix`'s column contract), `iter_operator_samples` over a manifest, and `operator_stats` (three normalizers — parameters, trunk coordinates, output — streamed). The **card** records the expanded `parameters` and their source `parameter_names`, the trunk construction and, for a budget, the **`trunk_indices`** it selected: the trunk is rebuilt from the mesh at predict time, and the indices are what detects a budget that selected different points (a JSON card is no place for a `(T, 3)` array). Prediction through the manifest takes each entry's parameters from its `Metadata`; `predict_mesh`/`predict_file` (and the MCP `predict_file` tool) take `parameters=` — a `deeponet` checkpoint without them is refused naming the keys it expects, and any other family given them is refused too. Under a budget the write-back fills **NaN, never 0**, at the points the model never saw. Neither `torch_geometric` nor `Augmentation` applies.

The worked example, [`deeponet_beam.py`](https://github.com/loumalouomega/meshioplusplus/blob/main/example/physicsnemo/deeponet_beam.py): ONE simplexified cantilever (`grid((20, 4, 4))`, L = 1, h = 0.2) and 200 cases differing only in `Metadata: {Load, Modulus, PoissonRatio}`, the truth an Euler–Bernoulli deflection plus a Timoshenko shear term so that `w = (P/E)·[f1(x) + f2(x)(1 + ν)]` is genuinely nonlinear in the parameters. The baseline is a per-node least-squares fit on `[1, P, E, ν]` over the training cases — exact if the map were linear, so the DeepONet's whole margin over it is the nonlinearity. 150 epochs in 21.7 s:

| | RMSE | relative L2 |
|---|---|---|
| least squares on `[1, P, E, ν]` | 5.20e-6 | 0.364 |
| DeepONet | **3.13e-7** | **0.0158** |

The DeepONet's error is about 17 times lower in RMSE and 23 times lower in relative L2 than the least-squares fit, and that margin is the nonlinearity and nothing else: a per-node linear map on `[1, P, E, ν]` cannot represent the `P/E` ratio or its product with `ν`, so the baseline's 0.364 relative error is where any linear surrogate stops on this problem. The run had converged — its best validation epoch was 108 of 150 — so the margin is not an artefact of stopping early.

## Temporal windows and rollout

A steady-state surrogate maps one mesh to one answer. A *transient* one maps a short history to the next state, and the shape it wants is a **window**: the last K states of every node, laid side by side along the channel axis.

```python
for entry_id, time, sample in mpn.iter_windows(manifest, fields=["T"], history_size=3):
    sample.arrays["x"]      # (N, K*W) — the window
    sample.arrays["y"]      # (N, W)   — the state that followed
    sample.x_columns        # ('T@-2', 'T@-1', 'T@-0') — the age is part of the contract
```

**Oldest first**, everywhere: the last `W` columns are always the most recent state. `make_window` is the single owner of that layout, shared by the three schemes and by `rollout`'s own history, because a model trained on one ordering and fed the other produces plausible, wrong numbers with nothing to flag them. The recorded `x_columns` name the age explicitly, so a stored contract catches the disagreement rather than leaving it to a convention nobody re-reads.

| scheme | samples per series | what it learns |
| --- | --- | --- |
| `single_step` | `T − K` | the autoregressive step: each window paired with the state that followed |
| `one_shot` | 1 | the first window paired with the last state — a model that jumps to the end |
| `time_conditional` | `T − K` | the first window plus a normalized time channel, so error does not compound |

The streaming invariant holds: `window_samples` keeps at most K **feature matrices** and one mesh alive, never K meshes, and reads each step exactly once.

**Rollout** is the measurement a one-step validation loss cannot make — feed the model its own prediction and watch the error grow:

```python
result = mpn.rollout(predict_fn, series, history_size=3, fields=["T"])
result.errors        # (T − K,) per-step RMSE against the truth
```

`predict_fn` takes one `(N, K*W)` window and returns the next `(N, W)` state. Injecting it keeps this module torch-free: the caller owns the framework, normalization and device, and meshio++ owns the loop. It is seeded with K *true* states and everything after that is its own output, which is what makes the drift real — the same evaluation under teacher forcing reports a flat error and tells you nothing.

## Dataset augmentation

A surrogate trained on a hundred solves of the same part in the same pose learns the pose. The hard half of fixing that was already here: [`transform(rotate_vector_data=True)`](transform.md) rotates vector and rank-2 tensor fields **coherently with the geometry**, point-located and (since v10.33.0) cell-located alike — a rotated mesh carrying an unrotated velocity field is a physically impossible sample a model will happily learn from.

What was missing is the wrapper:

```python
augmentation = mpn.Augmentation(rotation=True, scale={"low": 0.9, "high": 1.1}, seed=0)

for entry_id, time, sample in mpn.iter_samples(manifest, fields=["T"],
                                               augmentation=augmentation, epoch=3):
    ...
```

Every draw is a pure function of `(seed, epoch, index)` — no global generator, nothing carried between calls — so a run is reproducible, two processes of a distributed job agree, and re-running epoch 3 gives epoch 3's meshes rather than the next ones in a stream. Scaling is isotropic on purpose: an anisotropic scale changes what a vector field *means*.

**A paired sample gets one draw, replayed.** With `target_offset` or a `Target`, a sample is two meshes, and augmenting them independently would teach the model that a part rotates between one step and the next.

In a spec it is a top-level block, and it is refused by name for the `srresnet` family — a grid sampled on a fixed lattice would have its own coverage changed by rotating the mesh under it:

```jsonc
"Augmentation": { "Rotation": {"Axis": "z", "MaxDegrees": 45.0},
                  "Scale": {"Low": 0.9, "High": 1.1}, "Seed": 0 }
```

The trainer applies it to the **train** split only (a validation loss that moves because the poses moved measures nothing), calls the dataset's `set_epoch` between epochs, and records the description in the model card. Statistics are gathered un-augmented; edge norms are rotation-invariant, and a scale range only widens them slightly.

## Dataset manifests

The object the adapter iterates is a [`DatasetManifest`](./datasets) — a hand-editable JSON cataloguing many cases (each possibly a time series) with splits, tags, groups and notes, curated by hand, by the [`dataset` CLI group](./cli#meshioplusplus-dataset), or by the [MCP tools](./mcp), all reading and writing the same file.

## The worked example

[`example/physicsnemo/`](https://github.com/loumalouomega/meshioplusplus/tree/main/example/physicsnemo) is the end-to-end path, executed for real: 200 self-generated steady-heat cases (a manufactured Poisson pair `q = -∆T` on jittered, transformed triangle meshes — `convert_cells(simplexify)` + `transform`), catalogued and split with the `dataset` CLI, preprocessed by a [settings-document pipeline](./pipeline), trained with PhysicsNeMo's MeshGraphNet through `make_dataset`, and the predictions written back as ordinary `point_data` (`T_pred`/`T_error`) into `.vtu` files and rendered. The committed README, stats files and renders are the outputs of a real GPU run: 100 epochs over 160 training graphs in 73.8 s on an RTX 2000 Ada (8 GB, WSL2), validation MSE 7.9×10⁻¹ → 4.2×10⁻⁴, mean test RMSE 0.0040 on a field of amplitude 1 (2026-08-06; physicsnemo 2.1.1, torch 2.12.0+cu130).

The trainer is reachable from every surface the adapter is: `mpn.run_training`/`mpn.predict` in Python, `python -m meshioplusplus.physicsnemo.train --spec` as a process, and the `train_*` [MCP tools](./mcp#training) the [dashboard](./dashboard#launching-and-monitoring-a-run) drives.

**CI note:** public runners install neither torch nor PhysicsNeMo (the [GPU-handoff precedent](./gpu#testing-and-ci)). The pure half — `graph_sample`, the stats accumulators, manifest iteration, the install-error messages — runs in the default CI matrix with nothing optional installed; the gated halves `importorskip` and were exercised on a real GPU machine, which is stated here rather than implied by a green badge.
