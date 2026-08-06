# PhysicsNeMo integration

Feed meshes into [NVIDIA PhysicsNeMo](https://github.com/NVIDIA/physicsnemo)
training pipelines — MeshGraphNet-style GNNs first of all — without writing
the bespoke ingestion glue every project otherwise re-invents:

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

Everything here is pure Python over existing machinery — the
[`feature_matrix` contract and `edge_index`](./ml), the
[sequence machinery](./sequences), the [dataset manifests](./datasets) — and
the C++/WASM/C/Fortran core is untouched. The adapter lives in its own
subpackage (`meshioplusplus.physicsnemo`), deliberately **not** imported by
`import meshioplusplus`: like the [MCP server](./mcp), its gated halves have
their own Python floor and heavyweight dependencies, so the top-level package
must work without them.

## Installation: deliberately no `[physicsnemo]` extra

```bash
pip install nvidia-physicsnemo     # the framework (imports as `physicsnemo`)
pip install torch_geometric        # for the PyG Dataset path
```

One more, found the hard way while executing the worked example:
**MeshGraphNet's GNN layers import `torch_scatter` directly** (it is part of
PhysicsNeMo's `gnns` extra, not its core dependencies), and `torch_scatter`
ships prebuilt wheels only per torch/CUDA pair, lagging torch releases — at
the time of writing there is none for torch 2.13, so the example pins torch
2.12:

```bash
pip install "torch==2.12.0" --index-url https://download.pytorch.org/whl/cu130
pip install torch-scatter -f https://data.pyg.org/whl/torch-2.12.0+cu130.html
```

The adapter itself never imports `torch_scatter` — only training a
MeshGraphNet does.

The CuPy/torch packaging precedent, applied a third time:
`nvidia-physicsnemo` hard-depends on `torch>=2.10`, whose default Linux wheel
bundles CUDA at multiple GB — a `meshioplusplus[physicsnemo]` extra would
transitively pin exactly the dependency the repo already
[refuses to pin for `to_torch`](./ml#packaging-deliberately-no-torch-jax-extra).
Install the framework directly; `mpn.has_physicsnemo()` /
`mpn.has_torch_geometric()` answer availability without raising, and a
missing install raises a named error quoting the command above.

Two floors to know about: meshio++ supports Python ≥ 3.8, PhysicsNeMo
requires **≥ 3.11**. The mcp-SDK precedent applies — everything testable
(`graph_sample`, the stats, the manifest iteration) is pure numpy and runs on
every supported interpreter; only constructing the framework-facing
`Reader`/`Dataset` objects needs the framework.

## Reconnaissance (verified 2026-08-06, physicsnemo 2.1.1)

Written down before the adapter was built, per the roadmap's own rule; these
are the facts the design stands on. Package: `nvidia-physicsnemo` on PyPI
(there is no `physicsnemo` package), Apache-2.0, pure-Python wheel, quarterly
minor releases with real API churn — pin `nvidia-physicsnemo>=2.1,<2.2` in
training projects.

- **DGL is gone.** The DGL→PyTorch Geometric migration completed across
  25.08–25.11; `main` has zero DGL references. PyG is the only graph backend,
  and `torch_geometric` is an *extra* of PhysicsNeMo (`gnns`), not a core
  dependency. Nothing here targets DGL.
- **Two datapipe generations coexist.** Gen-1 (`physicsnemo.datapipes.gnn`)
  is plain `torch.utils.data.Dataset` subclasses yielding PyG `Data` objects
  (`.pos`, `.x`, `.y`, `.edge_index`, `.edge_attr`) — this is what the
  MeshGraphNet examples actually train on today. Gen-2
  (`physicsnemo.datapipes`) is a GPU-centric `Reader → Dataset → DataLoader`
  pipeline over `TensorDict`, whose stated user extension points are
  **Readers and Transforms, not dataset subclasses**; the `Reader` ABC is two
  abstract methods (`_load_sample(index) -> dict[str, torch.Tensor]`,
  `__len__`), and an optional `@register()` decorator makes a reader
  addressable from Hydra configs. The `Datapipe` base class carries only
  metadata and imposes no contract — the adapter does not target it.
- **The edge-feature convention is stable across examples**: with
  `row, col = edge_index`, the edge attribute is
  `cat((pos[row] - pos[col], ‖pos[row] - pos[col]‖))` — a 4-vector in 3D.
  Note the sign: source minus destination. `graph_sample` reproduces it
  exactly, pinned by a numeric test.
- **Normalization stats are plain JSON by convention**: Gen-1 datasets read
  CWD-relative `node_stats.json` / `edge_stats.json` with `{field}_mean` /
  `{field}_std` keys; Gen-2 has a `Normalize` transform taking means/stds.
  `field_stats`/`edge_stats` produce dicts in that key convention; meshio++
  defines no new stats file format.
- **`physicsnemo.mesh.Mesh`** (new in 2.1.0) is **simplicial-only** —
  points, lines, triangles, tetrahedra; hexes/wedges/pyramids are tessellated
  on ingest with data-loss warnings — and its only file on-ramp is a
  PyVista-backed `from_pyvista`/`to_pyvista` pair; its `.pmsh` save format is
  self-declared unstable. A `from_meshio` bridge is therefore **deferred**:
  it would buy nothing over the Reader/PyG path for training, and would pin
  meshio++ to the least stable surface in the framework. If PhysicsNeMo's
  `Mesh` stabilizes, an upstream `io_meshio.py` mirroring `io_pyvista`'s
  shape is the natural follow-up.
- **Every mesh file enters PhysicsNeMo through PyVista/VTK today** (the
  `VTKReader` hard-requires pyvista and reads `.stl`/`.vtp`/`.vtu`/`.vtk`
  only). That is the slot this adapter fills: meshio++'s 42 readers, regions,
  and preprocessing operations, ending in the tensors the datapipes expect.

## The adapter: `meshioplusplus.physicsnemo`

The `mcp/` split, applied again: the subpackage's `__init__` is pure
numpy/stdlib/meshioplusplus and holds **all** the behaviour; `_reader.py` is
the only module importing `physicsnemo`, `_pyg.py` the only one importing
`torch_geometric`, and both are reached through lazy factories.

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
| `y` | `(N, Fy)` float32¹ | a second `feature_matrix` over `target_fields`; absent when `target_fields=None` |
| `edge_index` | `(2, E)` **int64** | [`edge_index(mesh, kind=, undirected=)`](./ml#graphs-for-gnns-edge_index) verbatim |
| `edge_attr` | `(E, dim+1)` float32¹ | `cat((pos[row] − pos[col], ‖·‖))` — the PhysicsNeMo convention above |

¹ `float32=True` is the default (the training convention); `float32=False`
keeps meshio++'s canonical float64. `edge_index` is int64 always.

The column contract is `feature_matrix`'s, carried through: store
`s.schema` at training time, compare `x_columns` at inference time, and
feature drift is a named assertion instead of a silently wrong prediction.

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

All three walk a [`DatasetManifest`](./datasets) (or anything
`DatasetManifest.load` accepts), resolve each entry through the sequence
plan, and read one mesh at a time. The stats are per-component streaming
moments in the `{field}_mean`/`{field}_std` key convention — write them with
`json.dump` as `node_stats.json`/`edge_stats.json` and Gen-1 datapipes read
them as-is; pass them to Gen-2's `Normalize` transform directly.

### Framework objects: `make_dataset`, `make_reader`

```python
ds = mpn.make_dataset(manifest, split="train", fields=["q"], target_fields=["T"])
# torch.utils.data.Dataset of torch_geometric.data.Data(pos=, x=, y=, edge_index=, edge_attr=)

reader = mpn.make_reader(manifest, split="train", fields=["q"])
# physicsnemo Gen-2 Reader: __len__ + _load_sample(i) -> dict[str, torch.Tensor]
```

`make_dataset` is the training path: PyG's `DataLoader` batches
variable-size graphs natively, which is what MeshGraphNet consumes.
`make_reader` is the Gen-2 path for Hydra-driven pipelines (the reader is
`@register()`ed when the registry is importable — addressability is a bonus,
not a pinned contract); note TensorDict has no ragged-graph batching
convention, so batching across meshes of different sizes remains the PyG
path's job. Both index the flat `(entry, step)` sequence resolved once from
the manifest's plans, and read one mesh per `__getitem__`/`_load_sample`.

Deferred, deliberately, and recorded here rather than half-built:
autoregressive t→t+1 target pairing (`target_fields` is same-step in v1 —
the steady-state surrogate shape), and the `physicsnemo.mesh.Mesh` bridge
(see the reconnaissance section for why).

## Dataset manifests

The object the adapter iterates is a [`DatasetManifest`](./datasets) — a
hand-editable JSON cataloguing many cases (each possibly a time series) with
splits, tags, groups and notes, curated by hand, by the
[`dataset` CLI group](./cli#meshioplusplus-dataset), or by the
[MCP tools](./mcp), all reading and writing the same file.

## The worked example

[`example/physicsnemo/`](https://github.com/loumalouomega/meshioplusplus/tree/main/example/physicsnemo)
is the end-to-end path, executed for real: 200 self-generated steady-heat
cases (a manufactured Poisson pair `q = -∆T` on jittered, transformed
triangle meshes — `convert_cells(simplexify)` + `transform`), catalogued and
split with the `dataset` CLI, preprocessed by a
[settings-document pipeline](./pipeline), trained with PhysicsNeMo's
MeshGraphNet through `make_dataset`, and the predictions written back as
ordinary `point_data` (`T_pred`/`T_error`) into `.vtu` files and rendered.
The committed README, stats files and renders are the outputs of a real GPU
run: 100 epochs over 160 training graphs in 73.8 s on an RTX 2000 Ada (8 GB,
WSL2), validation MSE 7.9×10⁻¹ → 4.2×10⁻⁴, mean test RMSE 0.0040 on a field
of amplitude 1 (2026-08-06; physicsnemo 2.1.1, torch 2.12.0+cu130).

**CI note:** public runners install neither torch nor PhysicsNeMo (the
[GPU-handoff precedent](./gpu#testing-and-ci)). The pure half — `graph_sample`,
the stats accumulators, manifest iteration, the install-error messages — runs
in the default CI matrix with nothing optional installed; the gated halves
`importorskip` and were exercised on a real GPU machine, which is stated here
rather than implied by a green badge.
