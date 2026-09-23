# Partitioning

`meshioplusplus.partition(mesh, nparts)` decomposes a mesh into exactly `nparts` balanced pieces for domain decomposition — the **count-driven** complement to the criterion-driven [`split`](/split). It is a mesh **operation**, not a file format, and runs under every mesh backend.

![A mesh decomposed into 6 balanced parts (SFC / Hilbert curve)](/images/partition_pieces.png)

```python
import meshioplusplus

mesh = meshioplusplus.read("domain.msh")

# four balanced pieces (part id = list index)
pieces = meshioplusplus.partition(mesh, 4)
for part, piece in enumerate(pieces):
    meshioplusplus.write(f"domain_{part}.vtu", piece)

# just the assignment, one Int64 array per cell block
labels = meshioplusplus.partition_labels(mesh, 4)
mesh.cell_data["partition:part"] = labels

# weighted by a per-cell cost, tagged with the original indices
pieces = meshioplusplus.partition(mesh, 4, weights="cost", record_ids=True)
```

Both CLIs expose the same thing:

```bash
meshioplusplus partition domain.msh 'domain_{part}.vtu' --nparts 4
meshioplusplus partition domain.msh labelled.vtu --nparts 4 --labels-only
meshioplusplus partition domain.msh domain.pvtu --nparts 4 --ghost-layers 1
```

## Writing a partition as one dataset

An Elmer mesh partitioned by ElmerGrid (`partitioning.N/part.n.*`) reads back merged with the same `partition:part` labels, one per cell (v16.2.0, see [Elmer](formats/elmer.md#partitioned-meshes)).

A `.pvtu` (or `.pvtp`) is the ParaView index of a partition: one file that declares the arrays and names one piece per part, so the decomposition opens as a single dataset. `pvtu.write_pieces` takes exactly what `partition` returns — halo layers included, which a single mesh with a `partition:part` label cannot express, since a cell can be a ghost of several parts:

```python
pieces = meshioplusplus.partition(mesh, 4, ghost_layers=1)
meshioplusplus.pvtu.write_pieces("domain.pvtu", pieces)   # domain.pvtu + domain/domain_0000.vtu ... domain_0003.vtu

meshioplusplus.pvtu.write("labelled.pvtu", mesh_with_partition_part)   # or carve one mesh by its partition:part labels
```

`partition:ghost` becomes `vtkGhostType` on cells (`DUPLICATECELL` for any layer) and on points (`DUPLICATEPOINT` when no owned cell of the piece uses the point), and sets the index's `GhostLevel` to the deepest layer; `partition:ghost` itself is kept too, since `vtkGhostType` collapses layer 2 onto layer 1. Reading the index merges the pieces (one `piece_<i>` region each); `ghosts="drop"` (on `meshioplusplus.read`, or `--drop-ghosts` on `convert`) removes the halo and returns the partition of unity, so `partition` → `.pvtu` → read is the original mesh up to point ordering (`clean` with `weld=True` fuses the interface points). See [PVTU](/formats/pvtu).

## Two methods

### SFC (always available, the default fallback)

Cells are ordered along a **Hilbert space-filling curve** of their centroids (bounding box quantized to 21 bits per axis — the same key transforms [`reorder`](/reorder) uses, shared via `detail/space_filling.hpp`) and the curve is cut into `nparts` contiguous ranges:

- **Unweighted**: the cell at curve rank `r` of `n` goes to part `(r * nparts) / n`, so part sizes differ by **at most one cell**.
- **Weighted** (`weights=<cell_data name>`, a scalar numeric array): each cell goes to the ideal interval containing its weight midpoint, `min(nparts − 1, ⌊(S_r + w_r/2) · nparts / W⌋)` over the prefix sum `S_r`. Monotone along the curve, hence still contiguous; a part can end up empty when a single weight exceeds `W / nparts`.

The SFC path is **deterministic and byte-identical** across the three mesh backends, any thread count, and the C++-core/numpy-fallback boundary (a test pins the labels as identical). Curve locality keeps parts spatially compact, though it optimizes no explicit edge-cut metric.

### KaHIP (optional, the quality path)

With the optional [KaHIP](https://github.com/KaHIP/KaHIP) backend, the mesh's **dual graph** — cells are vertices, with an edge wherever two cells share a face (3D) or an edge (2D) — is partitioned by KaHIP's serial `kaffpa()` interface, which actively minimizes the edge cut (the halo communication volume of an MPI solver):

```python
pieces = meshioplusplus.partition(mesh, 16, method="kahip",
                                  imbalance=0.03, mode="eco", seed=0)
```

- `imbalance` — allowed part-size imbalance fraction, in (0, 1); default 0.03.
- `mode` — `"fast"`, `"eco"` (default) or `"strong"`. From the benchmarks of the Kratos KaHIPApplication this port derives from ([KratosMultiphysics/Kratos#14453](https://github.com/KratosMultiphysics/Kratos/pull/14453)): `fast` is not dramatically better than Metis — **`eco`/`strong` are where the edge-cut wins are**, at the price of runtime.
- `seed` — `kaffpa` is deterministic per seed, but the assignment may change between KaHIP versions; treat exact labels as unstable, balance/coverage as the contract.

Cells whose block has no facet table (ragged polygons, polyhedra, cells of lower dimension in a mixed mesh) become isolated graph vertices: they are still assigned, just without adjacency preferences.

`method="auto"` (the default) resolves to KaHIP when a backend is available, else SFC. Requesting `method="kahip"` without one **raises an error naming the option** — never a silent downgrade.

## Getting KaHIP

KaHIP is **MIT-licensed, like meshio++ itself**, so unlike GPL-licensed partitioners enabling it has no licensing implication whatsoever. It is a **bring-your-own** dependency: meshio++ never vendors or downloads it (the same policy the Kratos reviewers required of the original integration).

- **C++ core** (the fast path): build/install KaHIP, then configure with the option and prefix:

  ```bash
  git clone --depth 1 --branch v3.25 https://github.com/KaHIP/KaHIP.git
  cmake -S KaHIP -B KaHIP/build -DCMAKE_BUILD_TYPE=Release -DPARHIP=OFF -DNOMPI=ON \
        -DCMAKE_INSTALL_PREFIX="$HOME/kahip-install"
  cmake --build KaHIP/build -j && cmake --install KaHIP/build

  CMAKE_ARGS="-DMESHIOPLUSPLUS_WITH_KAHIP=ON -DKAHIP_ROOT=$HOME/kahip-install" \
      pip install meshioplusplus --no-binary meshioplusplus
  ```

`-DNOMPI=ON` is required, not just `-DPARHIP=OFF`: KaHIP's own CMakeLists.txt does an unconditional `find_package(MPI REQUIRED)` unless `NOMPI` is set, since `kaffpaE` (not only ParHIP) also depends on MPI — meshio++ only ever links the serial `kaffpa` interface, so this build has no MPI dependency at all despite what KaHIP's own default configure implies.

  ::: warning KaHIP must be findable at *runtime* too
KaHIP is a shared library, and the wheel build strips the RPATH from `_core` so the wheel stays relocatable. If the prefix is not on the default loader path, importing meshio++ then fails with

  ```
  ImportError: libkahip.so: cannot open shared object file: No such file or directory
  ```

even though the build succeeded. Put the prefix on the loader path:

  ```bash
  export LD_LIBRARY_PATH="$HOME/kahip-install/lib:$LD_LIBRARY_PATH"   # Linux
  export DYLD_LIBRARY_PATH="$HOME/kahip-install/lib:$DYLD_LIBRARY_PATH"  # macOS
  ```

or install KaHIP to a system prefix. This is the ordinary cost of a bring-your-own shared dependency; if you would rather not manage it, the Python-only route below needs no such setup.
  :::

The in-repo `cmake/FindKaHIP.cmake` honours `KAHIP_ROOT`/`KAHIP_DIR` (environment or cache), standard prefixes, and pkg-config; an in-source KaHIP build (`interface/` + `build/`) works as a prefix too. On macOS, `brew install KaHIP/kahip/kahip` provides an installable prefix. Package managers: Conan option `with_kahip`, vcpkg feature `kahip` (both map the flag only — you still supply the install; KaHIP is on neither registry).

- **Python-only** (no rebuild): `pip install meshioplusplus[kahip]` pulls the MIT [`kahip` wheel](https://pypi.org/project/kahip/); the pure-Python fallback builds the dual graph in numpy and calls `kahip.kaffpa` — same quality, Python-speed graph assembly.

Only the **serial** `kaffpa` interface is ever linked — never ParHIP — so enabling KaHIP adds **no MPI dependency**. `_core.__has_kahip__` reports the C++ backend at runtime.

### 32/64-bit indices

KaHIP built with `-D64BITMODE=On` uses 64-bit graph indices. meshio++ queries the installed library's width at runtime (`kahip_sizeof_idx()`) and builds the CSR arrays to match, so either build works unmodified. Against a 32-bit KaHIP, a dual graph whose directed-edge count overflows `int32` fails with an error telling you to rebuild KaHIP with `-D64BITMODE=On`; the cell count itself is limited to `int32` by the `kaffpa` interface in both builds.

## Output contract

- Exactly `nparts` pieces are returned, part id ascending; pieces may hold empty blocks. Every piece keeps the input's cell-block structure **1:1** (deliberately unlike `split`, which drops empty blocks), so the per-block cell maps index by input block and **concatenating the pieces reproduces the input** — every cell lands in exactly one piece (partition of unity).
- `point_data`/`cell_data` rows are carried into the owning piece; `field_data` is copied to every piece. `point_sets`/`cell_sets` are remapped per piece (Python layer only, as with `split`).
- `record_ids=True` attaches Int64 `partition:original_point_id` `point_data` and `partition:original_cell_id` `cell_data` with the original input indices.
- `partition_labels` returns one Int64 array per cell block (block-aligned, ready to attach as `cell_data`), values in `[0, nparts)`.
- `ghost_layers=N` grows each piece by `N` **shared-node BFS layers** of cells owned by other parts — the halo an MPI-style domain decomposition needs. Each ghosted piece carries an Int64 `partition:ghost` `cell_data`: `0` for a cell this part owns, `L` for one reached at layer `L`. Shared-node (rather than shared-facet) adjacency is the conservative choice: it is what a node-based assembly actually needs.

The pieces then **overlap**, so the partition-of-unity property above holds only for `ghost_layers=0` (the default). `partition_labels` therefore does not accept it at all — a flat per-cell label array is the *ownership* map, and a cell can be a ghost of several parts at once. The numpy fallback builds disjoint pieces only and raises `NotImplementedError` rather than silently returning unghosted ones.

The [pure-Python KaHIP route](#kahip-optional-the-quality-path) (the `kahip` PyPI wheel against a KaHIP-less `_core`) cannot ghost either: with `ghost_layers > 0`, `method="auto"` routes to the C++ core instead (auto is a preference order, and the core is the one ghost-capable path), while an explicit `method="kahip"` raises naming `-DMESHIOPLUSPLUS_WITH_KAHIP=ON` — never a silent SFC downgrade.

  ```python
  pieces = mp.partition(mesh, 4, ghost_layers=1)
  owned = pieces[0].cell_data["partition:ghost"][0] == 0
  ```
- v3.25's `--connected_blocks` (force each part to be a connected subgraph) is not reachable through the `kaffpa` C interface, so it is not exposed — a documented TODO should a later KaHIP add it to the API.

## Other language surfaces

The operation exists on every binding surface: C (`mio_partition`/`mio_partition_labels`, see [the C API](/c_api)), Fortran (`m%partition`/`m%partition_labels`, see [Fortran](/fortran)), WebAssembly (`partition`/`partitionLabels`, SFC only — KaHIP is never part of the WASM build, see [WASM](/wasm); `partition`'s `returnMaps` option gives each piece a `pointMap`/`cellMaps`), and the `partition` verb in both [CLIs](/cli).

## meshio++ is serial

`partition` computes a decomposition; it does not run one. There is no MPI in the library, no distributed reader or writer and no communicator anywhere in the API. The intended workflow is to partition on one rank (or in a pre-step) and have each rank read or receive its own piece; `ghost_layers > 0` produces the shared-node halo an MPI assembly needs. See `doc/cpp_api.md`.
