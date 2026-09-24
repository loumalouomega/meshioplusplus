# Elmer mesh directory

An Elmer mesh is what [ElmerSolver](https://www.elmerfem.org) reads with `Mesh DB`: a **directory**, not a file, holding plain-text `mesh.header`, `mesh.nodes`, `mesh.elements`, `mesh.boundary` and optionally `mesh.names` (or, with `ElmerGrid -bin`, binary `mesh.*.bin` files). ElmerGrid writes it from gmsh, UNV, Abaqus and other meshes, and its `-partition`/`-metis` options add a `partitioning.N` subdirectory for parallel runs.

| | |
|---|---|
| **Format name** | `elmer` |
| **Extensions** | none: a directory is recognised by content (a `mesh.header`, a `partitioning.N` directory of `part.n.*` files, or a directory holding one) |
| **Read / Write** | ✓ / ✓ (text files; binary files are read only) |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("box")                     # the directory; sniffed, no format needed
mesh = meshioplusplus.read("box/mesh.header")         # the same mesh
part = meshioplusplus.read("box", piece=0)            # one part of a partitioned mesh
meshioplusplus.write("out", mesh, file_format="elmer")  # a directory has no extension to infer from

mesh.cell_data["partition:part"] = meshioplusplus.partition_labels(mesh, 4)
meshioplusplus.write("out", mesh, file_format="elmer")  # also writes out/partitioning.4
```

```bash
meshioplusplus convert box box.vtu          # read: the directory is sniffed
meshioplusplus convert box.msh out -o elmer # write: name the format
```

`read` also takes `lenient=True`, which skips element types with no meshio++ cell type (with a warning) instead of failing. Both engines (the C++ core and the pure-Python reference) read the same meshes and write the same bytes.

## Files

The layouts follow ElmerGrid's own writer and ElmerSolver's reader:

| File | Lines |
|---|---|
| `mesh.header` | `nodes elements boundary-elements`, then the number of element types, then one `type count` line per type |
| `mesh.nodes` | `id part x y z` (`part` is `-1` in a serial mesh) |
| `mesh.elements` | `id body type node…` |
| `mesh.boundary` | `id boundary parent1 parent2 type node…` (a missing parent is `0`) |
| `mesh.names` | `$ name = id` lines, after a `names for bodies` or `names for boundaries` heading |

Element types are family × 100 + node count:

| Code | Cell type | Code | Cell type |
|---|---|---|---|
| `101` | `vertex` | `504` / `510` | `tetra` / `tetra10` |
| `202` / `203` / `204` | `line` / `line3` / `line4` | `605` / `613` | `pyramid` / `pyramid13` |
| `303` / `306` / `310` | `triangle` / `triangle6` / `triangle10` | `706` / `715` / `718` | `wedge` / `wedge15` / `wedge18` |
| `404` / `408` / `409` | `quad` / `quad8` / `quad9` | `808` / `820` / `827` | `hexahedron` / `hexahedron20` / `hexahedron27` |

### Binary files

`ElmerGrid -bin` writes `mesh.nodes.bin`, `mesh.elements.bin` and `mesh.boundary.bin` (and `-sbin` a single-precision `mesh.nodes.sbin`) in place of the text files, which ElmerSolver reads as stream files. Since v16.10.0 they are read, in a serial mesh and in the parts of a partitioned one, each file in whichever form the directory holds. Every record is native-endian 32-bit integers followed, for nodes, by 64-bit (`.bin`) or 32-bit (`.sbin`) reals, with no separators:

| File | Record |
|---|---|
| `mesh.nodes.bin` / `.sbin` | `id x y z` |
| `mesh.elements.bin` | `id part body type node…` |
| `mesh.boundary.bin` | `id part boundary parent1 parent2 type node…` |

The part is the owner a halo copy names, else 0. Nothing marks the byte order, so it is taken from the first record's id: of the two readings of it, the file is read in the order that gives the smaller positive id. `mesh.header` and `mesh.names` stay text. The fixtures are ElmerGrid's own `-bin` and `-sbin` output of text meshes in the suite, and they read identically to those meshes. The writer writes text only.

Codes with no meshio++ cell of the same nodes (`102` periodic pairs, `412`/`416` quads, the `614` pyramid and higher orders) are a `ReadError` naming the code, or skipped under `lenient=True`.

## Mapping

- **Bulk and boundary elements are separate cell blocks**, one per type code in file order, bulk first.
- **Bodies and boundaries are `cell` regions.** The region's tag is the body or boundary id, and its name comes from `mesh.names`, or is `body_<id>` / `boundary_<id>` when the file does not name it. A name used for both a body and a boundary becomes `body:<name>` and `boundary:<name>`.
- Boundary elements stay cells, not `side` regions: Elmer allows boundary elements with no parent (point conditions, `101`) and between two bodies, which a side region cannot express. The parent references are read and not kept; the writer regenerates them.

## Node order

Only the `820` and `827` bricks differ from meshio++'s (VTK's) order, and the registry's `"elmer"` tables handle them (see [node ordering](../node_ordering.md)). Elmer lists the four vertical mid-edge nodes before the top ring, and puts the `827` mid-height face centres at `v=−1, u=+1, v=+1, u=−1`. The tables come from the reference coordinates in ElmerSolver's `elements.def`, and they match the permutation ElmerSolver's own VTU writer applies (`Elmer2VtkIndexes`). Every other code, including `715`, `718` and `613`, is already in VTK order.

This was checked outside the repository against ElmerGrid built from source: gmsh `tet10`, `hex20`, `hex27` and `wedge15` meshes converted by `ElmerGrid 14 2` read back with every cell's nodes at the same coordinates, in the same order, as meshio++'s own gmsh reader gives.

## Partitioned meshes

ElmerGrid's `-partition`/`-metis` options write `partitioning.N/part.n.{header,nodes,elements,boundary,shared}`, with shared nodes listed by every part that uses them and, with `-halo`, halo copies of elements written as `id/owner`.

- **Reading merges the parts** into one mesh. Each shared node and each halo element is kept once, and each cell's 0-based part is in `cell_data["partition:part"]`, the key the [partition](../partition.md) operation uses.
- The path may be the `partitioning.N` directory itself, or a mesh directory holding exactly one. A serial mesh next to exactly one `partitioning.N` is read serially and labelled from it. With several `partitioning.*` directories, the serial mesh is read without labels (with a warning) and `piece` is an error until you name one.
- `piece` (C++ `ReadOptions::mPiece`) reads one part alone, halo copies included.
- `part.n.shared` is not read.

Merged reads of ElmerGrid partitions, with and without `-halo`, were checked against the serial mesh: the same cells, and each cell's label matching its part. ElmerGrid's `-bin -partition` output was not usable as a fixture (it leaves out the part headers), so binary parts are only covered by the parser they share with serial binary files.

## Writing

- The writer creates the directory if needed and writes `mesh.header`, `mesh.nodes`, `mesh.elements`, `mesh.boundary` and `mesh.names`. A mesh with `cell_data["partition:part"]` is also written partitioned; see [Partitioned write](#partitioned-write). Without it, an existing `partitioning.*` directory is left in place, with a warning, since it no longer matches.
- **Bulk and boundary.** The cells of the highest dimension are the bulk elements. Every lower-dimensional cell is a boundary element, and so is every facet of a `side` region.
- **Parents are regenerated.** A boundary facet's parents are the bulk cells that share it (two for an interface). An edge of a solid, or a point, gets the first bulk cell that holds all its nodes. A boundary element on no bulk cell is written with parent `0`, with a warning.
- **Ids come from regions**, in name order: a `cell` region over bulk cells is a body and one over lower-dimensional cells a boundary, and a `side` region is a boundary. A region's tag is its id when positive and unused; otherwise it gets the next free id. Cells in no region get one fresh id per block. A cell in two regions goes to the first.
- `mesh.names` names every region written. It is always written, even with no names, and it carries the `!` provenance block, which holds no `$` so ElmerSolver never mistakes it for a name. A boundary sharing a body's name is written as `<name>_boundary`, because ElmerSolver matches `mesh.names` against bodies and boundaries alike.
- Coordinates are written as `%.17g`, which round-trips a double exactly. Node and element ids are renumbered from 1.
- **Dropped, with a warning and a provenance note:** `point` regions (Elmer has no node set) and every data array (an Elmer mesh holds none) other than `partition:part`. A cell type with no Elmer code (polygons, polyhedra, higher-order Lagrange cells) is a `WriteError`, and so is a mesh with no cells of dimension 1 or more.

The writer was checked outside the repository: ElmerGrid re-reads every directory it writes without a warning, and ElmerSolver's heat equation on a written `tet10` mesh gives the same temperatures (to 1e-12) as on ElmerGrid's own conversion of the same gmsh mesh, and on a written `hex27` block the exact solution.

### Partitioned write

Since v16.10.0 a mesh carrying `cell_data["partition:part"]` (0-based parts, as the [partition](../partition.md) operation and the partitioned readers give) is written as the serial mesh plus `partitioning.N/part.k.{header,nodes,elements,boundary,shared}`, `N` being the largest part + 1, in ElmerGrid's layout:

- Each bulk element goes to its part, keeping its serial id, body and nodes. The global node and element ids are the serial mesh's.
- A part's nodes are those its elements use. A node used by several parts is shared: it is listed by each of them, and `part.k.shared` gives `id count owner others…`, the owner being the lowest part (1-based in the file).
- A boundary element goes to every part holding one of its parents, with the parent in another part written as 0.
- `part.k.header` counts nodes, elements, boundary elements and the element types, then the shared nodes.
- Halo elements (`ElmerGrid -halo`) are not written. A negative part is a `WriteError`, and an empty part is a warning, since ElmerSolver needs every part populated.

Reading the result merges back the input mesh with its labels. The layout was checked outside the repository against elmerfem (built from source at `a8a13b5`): a heat equation run with `mpirun -np 2 ElmerSolver_mpi` on a written two-part mesh gives the serial run's temperatures to 2e-14.

## Notes

- **Not read:** `.sif` files and results. ElmerSolver writes results as VTU (`Post File = x.vtu`), but in VTK's *raw* appended encoding, which meshio++'s VTU reader does not yet read.
- An Elmer directory cannot be one step of a [sequence](../sequences.md) found by glob: list its path explicitly.
- A body and a boundary share one id space in `mesh.names` only by name: the ids of bodies and of boundaries are separate, so `body_1` and `boundary_1` are different groups.

The layouts above follow ElmerGrid's `SaveElmerInput`/`LoadElmerInput` (`elmergrid/src/egnative.c`) and ElmerSolver's `ReadTargetNames` (`fem/src/MeshIO.F90`) in the [elmerfem source](https://github.com/ElmerCSC/elmerfem); see also the [ElmerSolver manual](https://www.nic.funet.fi/index/elmer/doc/ElmerSolverManual.pdf) and the [ElmerGrid manual](https://www.nic.funet.fi/index/elmer/doc/ElmerGridManual.pdf).
