# PhysicsNeMo integration

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

Written down before the adapter was built, per the roadmap's own rule; these are the facts the design stands on. Package: `nvidia-physicsnemo` on PyPI (there is no `physicsnemo` package), Apache-2.0, pure-Python wheel, quarterly minor releases with real API churn — pin `nvidia-physicsnemo>=2.1,<2.2` in training projects.

- **DGL is gone.** The DGL→PyTorch Geometric migration completed across 25.08–25.11; `main` has zero DGL references. PyG is the only graph backend, and `torch_geometric` is an *extra* of PhysicsNeMo (`gnns`), not a core dependency. Nothing here targets DGL.
- **Two datapipe generations coexist.** Gen-1 (`physicsnemo.datapipes.gnn`) is plain `torch.utils.data.Dataset` subclasses yielding PyG `Data` objects (`.pos`, `.x`, `.y`, `.edge_index`, `.edge_attr`) — this is what the MeshGraphNet examples actually train on today. Gen-2 (`physicsnemo.datapipes`) is a GPU-centric `Reader → Dataset → DataLoader` pipeline over `TensorDict`, whose stated user extension points are **Readers and Transforms, not dataset subclasses**; the `Reader` ABC is two abstract methods (`_load_sample(index) -> dict[str, torch.Tensor]`, `__len__`), and an optional `@register()` decorator makes a reader addressable from Hydra configs. The `Datapipe` base class carries only metadata and imposes no contract — the adapter does not target it.
- **The edge-feature convention is stable across examples**: with `row, col = edge_index`, the edge attribute is `cat((pos[row] - pos[col], ‖pos[row] - pos[col]‖))` — a 4-vector in 3D. Note the sign: source minus destination. `graph_sample` reproduces it exactly, pinned by a numeric test.
- **Normalization stats are plain JSON by convention**: Gen-1 datasets read CWD-relative `node_stats.json` / `edge_stats.json` with `{field}_mean` / `{field}_std` keys; Gen-2 has a `Normalize` transform taking means/stds. `field_stats`/`edge_stats` produce dicts in that key convention; meshio++ defines no new stats file format.
- **`physicsnemo.mesh.Mesh`** (new in 2.1.0) is **simplicial-only** — points, lines, triangles, tetrahedra, exactly one kind per `Mesh` (mixed dimensions are "separate Mesh objects", its own docs) — and its only file on-ramp is a PyVista-backed `from_pyvista`/`to_pyvista` pair; its `.pmsh` save format is self-declared unstable, and its `io` module is a hardcoded PEP-562 set with no plugin registry, so a bridge can only live meshio++-side. v9.28.0 deferred that bridge; **v9.30.0 builds it anyway** ([`to_physicsnemo`/`from_physicsnemo`](#the-physicsnemomesh-bridge)) with the risk stated rather than avoided: `physicsnemo.mesh` is the framework's newest and least stable surface, so the bridge touches only the tensorclass constructor and public attributes, never `.pmsh`, and training projects should pin `nvidia-physicsnemo>=2.1,<2.2`.
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

¹ `float32=True` is the default (the training convention); `float32=False` keeps meshio++'s canonical float64. `edge_index` is int64 always.

The column contract is `feature_matrix`'s, carried through: store `s.schema` at training time, compare `x_columns` at inference time, and feature drift is a named assertion instead of a silently wrong prediction.

### Streaming over a dataset: `iter_samples`, `field_stats`, `edge_stats`

```python
for entry_id, time, sample in mpn.iter_samples(manifest, split="train",
                                               fields=["q"], target_fields=["T"]):
    ...   # one mesh alive at a time — the sequence streaming invariant

stats = mpn.field_stats(manifest, split="train", fields=["q", "T"])
# {"q_mean": [...], "q_std": [...], "T_mean": [...], "T_std": [...]}
estats = mpn.edge_stats(manifest, split="train")
# {"edge_attr_mean": [...], "edge_attr_std": [...]}
```

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

**The risk, stated**: this rides the framework's newest surface. The bridge therefore touches only the tensorclass constructor and public attributes — never the self-declared-unstable `.pmsh` format — its gated module is the only one importing `physicsnemo.mesh` (an import that pulls in NVIDIA Warp, ~1.5 s, which is why it stays lazy), and training projects should pin `nvidia-physicsnemo>=2.1,<2.2`. If upstream's `Mesh` stabilizes, an `io_meshio.py` contributed upstream (mirroring `io_pyvista`'s shape) is still the natural end state; today its `io` module has no plugin registry, so the bridge lives here.

## Dataset manifests

The object the adapter iterates is a [`DatasetManifest`](./datasets) — a hand-editable JSON cataloguing many cases (each possibly a time series) with splits, tags, groups and notes, curated by hand, by the [`dataset` CLI group](./cli#meshioplusplus-dataset), or by the [MCP tools](./mcp), all reading and writing the same file.

## The worked example

[`example/physicsnemo/`](https://github.com/loumalouomega/meshioplusplus/tree/main/example/physicsnemo) is the end-to-end path, executed for real: 200 self-generated steady-heat cases (a manufactured Poisson pair `q = -∆T` on jittered, transformed triangle meshes — `convert_cells(simplexify)` + `transform`), catalogued and split with the `dataset` CLI, preprocessed by a [settings-document pipeline](./pipeline), trained with PhysicsNeMo's MeshGraphNet through `make_dataset`, and the predictions written back as ordinary `point_data` (`T_pred`/`T_error`) into `.vtu` files and rendered. The committed README, stats files and renders are the outputs of a real GPU run: 100 epochs over 160 training graphs in 73.8 s on an RTX 2000 Ada (8 GB, WSL2), validation MSE 7.9×10⁻¹ → 4.2×10⁻⁴, mean test RMSE 0.0040 on a field of amplitude 1 (2026-08-06; physicsnemo 2.1.1, torch 2.12.0+cu130).

**CI note:** public runners install neither torch nor PhysicsNeMo (the [GPU-handoff precedent](./gpu#testing-and-ci)). The pure half — `graph_sample`, the stats accumulators, manifest iteration, the install-error messages — runs in the default CI matrix with nothing optional installed; the gated halves `importorskip` and were exercised on a real GPU machine, which is stated here rather than implied by a green badge.
