# GiD postprocess (`.post.msh` / `.post.res` / `.post.bin` / `.post.h5`)

The [GiD](https://www.gidsimulation.com/) postprocess format. **Writing** goes through a vendored hardcopy of CIMNE's [gidpost 2.14](https://www.gidsimulation.com/downloads/gidpost-2-14-library-to-write-postprocess-results-for-gid-in-ascii-binary-or-hdf5-format/) (`src/cpp/third_party/gidpost/`, BSD-2-Clause-Views). gidpost's public API has **zero read functions**, so **reading** (added in v10.19.0) is meshio++'s own code against the on-disk grammar, independent of the vendored library — which is why the two directions have different build requirements, below.

| | |
|---|---|
| **Format name** | `gid` |
| **Extensions** | `.post.msh` / `.post.res`, `.post.bin`, `.post.h5` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | *Writing*: zlib (hard — gidpost deflates unconditionally), plus HDF5 for the `hdf5` flavour. *Reading*: nothing for ascii, zlib for binary, HDF5 for `hdf5`. |

## Writing

```python
import meshioplusplus

meshioplusplus.write("out.post.msh", mesh)                       # ascii, inferred
meshioplusplus.gid.write("out.post.bin", mesh, mode="binary")
meshioplusplus.gid.write("out.post.h5", mesh, mode="hdf5")
```

- **`mode`** — `"auto"` (default, inferred from the extension: `.post.bin` → binary, `.post.h5` → hdf5, anything else including `.post.msh`/`.post.res` → ascii), `"ascii"`, `"binary"`, `"hdf5"`.
- **`analysis_name`** — the GiD "analysis name" every result is grouped under (default `"meshio++"`).
- **`step`** — the single time/load step every result is written at (default `1.0`).

Three on-disk flavours:

- **ascii** — two sibling files, `<stem>.post.msh` (geometry, human-readable) + `<stem>.post.res` (results).
- **binary** — one deflated file, `<stem>.post.bin`.
- **hdf5** — one HDF5 file, `<stem>.post.h5`; needs a build with `MESHIOPLUSPLUS_WITH_HDF5=ON` in addition to gidpost itself (gidpost's own HDF5 flavour additionally needs the HDF5 *high-level* library, `libhdf5_hl`, alongside the core C API every other HDF5-backed format here already links).

A build without gidpost (`-DMESHIOPLUSPLUS_WITH_GIDPOST=OFF`, or gidpost on but zlib off) still exposes `gid.write` — it raises a `WriteError` naming the missing CMake flags rather than the path silently falling through to another format. **The statically-linked release CLI binaries and the Windows wheels build with zlib off and therefore do not carry `gid`** (documented, not a bug — see `CLAUDE.md`'s "GiD postprocess" note).

## Reading

```python
mesh = meshioplusplus.read("out.post.msh")              # any of the four spellings
mesh = meshioplusplus.gid.read("out.post.msh", time_step=1)
```

- **`time_step`** — selects one step of a multi-step results file (`0` = the first, negative counts from the end). Honoured natively, because no caller-side filter can recover a step that was never decoded; out of range is an error naming the available count.

The flavour is resolved from the extension and then **confirmed against the leading bytes**, so a gzipped `.post.msh` (gidpost's `GiD_PostAsciiZipped`, which is the same ASCII text through `gzprintf`) reads correctly even though its extension cannot say so.

**Sibling policy** (the `tetgen`/`triangle` `.node`/`.ele` precedent): the geometry file `<stem>.post.msh` is **mandatory**, the results file `<stem>.post.res` is **optional** — a mesh with no results reads back as geometry only. Passing the `.post.res` path directly derives and reads the `.post.msh`; *its* absence is an error, since results alone carry no geometry.

**Reading needs no gidpost at all**, and that has a visible consequence: `gid` is readable in **strictly more build configurations than it is writable**. The statically-linked release CLI binaries and the Windows wheels build with zlib off, so they cannot write GiD — but they read the ASCII flavour fine. `gid_available` reports the write side and `gid_readable` the read side; they genuinely differ, which is why there are two.

### Real-world variants the reader handles

Files in the wild differ from meshio++'s own output in two ways, both observed in a genuine Kratos-produced file and neither reproducible through our writer:

- **The full node table may be repeated in every `MESH` block** (rather than written once, with empty `Coordinates` pairs thereafter). Node ids are therefore accumulated into one global table de-duplicated by id, so a repeated-but-identical table is a no-op. Gapped and non-contiguous ids are supported.
- **Element ids may restart at 1 in each block** (rather than being globally unique, as our writer makes them). Element ids are tracked per block, so results resolve through the Gauss-point set's mesh name rather than a global id map.

A trailing material column reads back as `cell_data["gmsh:physical"]` — the exact inverse of the key the writer consumes. `Nnode` from the `MESH` header is the **only** disambiguator for that column, since GiD writes no separator before it. An all-zero column means "no materials" and is deliberately *not* surfaced: the binary and HDF5 writers always emit one, so materializing it would invent data and make the three flavours disagree on a round trip.

### What reading does not recover

- **Gauss-point results with more than one point per element** are dropped with a warning. meshio++'s `cell_data` is `(n,)`/`(n,k)`, never per-node-within-cell — the same structural limit MED's ELNO/ELGA documents. Averaging or taking the first point would invent data.
- **Tensor result types** (`Matrix`, `MainMatrix`, `Complex*`) are read as a single k-component array under the declared name; the declared *type* is dropped rather than reinterpreted, because GiD's symmetric-tensor component order differs from meshio/VTK's. Note the resulting **round-trip asymmetry**: the writer splits a k∉{1,2,3} array into k scalars, and the reader does not re-join them.
- `ResultGroup`, `OnNurbs*` locations, mesh groups (`Group`/`End Group`) and range tables are skipped or refused by name rather than guessed at.

## Cell types

GiD has exactly ten element types; higher-order variants share a type with a larger node count. One meshio++ cell block becomes one named GiD mesh (`"<celltype>_<blockindex>"`, since GiD results reference meshes by name).

| meshio++ type | GiD element type | Nnode |
|---|---|---|
| `vertex` | `Point` | 1 |
| `line` / `line3` | `Linear` | 2 / 3 |
| `triangle` / `triangle6` | `Triangle` | 3 / 6 |
| `quad` / `quad8` / `quad9` | `Quadrilateral` | 4 / 8 / 9 |
| `tetra` / `tetra10` | `Tetrahedra` | 4 / 10 |
| `hexahedron` / `hexahedron20` | `Hexahedra` | 8 / 20 |
| `wedge` | `Prism` | 6 |
| `pyramid` | `Pyramid` | 5 |

Every entry is **identity** (no node permutation) — independently cross-checked against Kratos Multiphysics's own production GiD writer (`kratos/includes/gid_mesh_container.h`), whose only reordering is for `hexahedron20`, and whose own internal hexahedra20 convention turns out to differ from GiD's for a reason unrelated to meshio++: Kratos's own node order (corners, bottom ring, verticals, top ring) differs from GiD's (corners, bottom ring, top ring, verticals) — and GiD's convention is, edge for edge, identical to meshio++'s own `hexahedron20` table. Kratos's own geometry classes for `triangle6`/`tetra10`/`quad8` apply no reorder at all when writing to GiD, and their edge tables match meshio++'s own edge-for-edge too. Pinned by `tests/cpp/test_gid.cpp`'s `GidOrdering` suite, which checks the raw written file against GiD's own geometry — not a round trip through a reader that does not exist.

**Not yet supported** — throws a `WriteError` naming the type, never a guess: `hexahedron27`, `wedge15`, `pyramid13` (orderings not independently verified); `polygon`/`polyhedron` (GiD has no such type); every `VTK_LAGRANGE_*` and higher-degree Lagrange type.

## Geometry and ids

One shared node table is written on the **first** mesh only (`GiD_fWriteCoordinatesBlock`); every subsequent mesh gets an empty `Coordinates`/`End Coordinates` pair (gidpost's own state machine rejects a mesh with neither). Points are always 3-D, z padded with 0 for a 2-D mesh. Element ids are globally unique 1-based integers across every block — gidpost's own bulk writer numbers `1..n` *per mesh*, which would collide the moment a mesh has more than one block, so meshio++ assigns a running id across all blocks instead.

An integral `cell_data` array named `"gmsh:physical"` is written as each element's material id (`GiD_fWriteElementsIdMatBlock`) and excluded from the result output, since it is already in the geometry file. No other key is consulted for material ids.

## Data mapping

`point_data` is written `GiD_OnNodes`. `GiD_ResultLocation` has no "on cells" concept, so `cell_data` is written `GiD_OnGaussPoints` against a synthetic one-point Gauss-point set declared once per block (`"gp_<mesh_name>"`) — the standard GiD idiom for a per-element field, and how Kratos's own GiD writer represents one. An array spanning several blocks becomes several result blocks sharing one result name but different Gauss-point sets.

`GiD_ResultType`'s valid component counts are irregular (Scalar 1; Vector 2/3/4; Matrix 3/6; MainMatrix 12; …), and gidpost does not validate an unsupported count itself — it silently emits a malformed file. meshio++ validates instead: 1 component → `GiD_Scalar`; 2 or 3 → `GiD_Vector`; anything else splits into that many named `GiD_Scalar` results (`"<name>_1"` … `"<name>_k"`, recorded as a provenance note). A 6-component array is **deliberately not** mapped to `GiD_Matrix`, even though stress tensors are GiD's canonical use case: meshio++'s `(n,6)` carries no declaration that it *is* a symmetric tensor, and GiD's own component order (`xx,yy,xy,zz,xz,yz`) differs from meshio/VTK's (`xx,yy,zz,xy,yz,xz`) — mapping on shape alone would silently permute six possibly-unrelated scalars. Splitting is lossless and unambiguous; a tensor-aware side channel is a documented follow-up.

Named regions, `field_data`, and multi-step results are **not carried** (each dropped with a provenance note where applicable); every write emits exactly one step (`step`).

## Provenance

`SlotTier::Block`: the provenance block renders as one GiD "user attribute" per line (`GiD_fWriteMeshUserAttribute` on the geometry file, `GiD_fWriteResultUserAttribute` on the results side — `"meshio++"` for line 0, `"meshio++_<n>"` for any further line under an open scope), which gidpost itself renders as a `# Name: value` comment in ascii/binary and as an HDF5 group attribute otherwise. Confirmed via the ascii flavour and via `h5py` inspection of the hdf5 flavour's `Meshes/1` group attributes; the *result*-file attribute's exact HDF5 placement (which "current result" it associates with, absent an open result block at the time it is written) is not independently verified and may not surface identically across every GiD HDF5 reader — a known limitation of the least-tested flavour, not a claimed guarantee.

There is **no pure-Python reference writer** for this format: gidpost cannot be cheaply reimplemented in pure Python, and a second implementation of its node-ordering permutations risks disagreeing with the first near exactly the cases the bytes tests exist to pin. `meshioplusplus.gid` is therefore a C++-core-only surface, the `openfoam` writer's precedent generalized (here there is no fallback engine at all, not merely one that goes unused).

## Quirks & limitations

- Reading is supported for all three flavours as of v10.19.0. `ResultGroup` blocks (several results packed into one wide `Values` row) are refused by name — meshio++ has never written one, and unpacking needs the per-member `ResultDescription` dimensions.
- The ascii flavour spans two files and **cannot** be read from or written to a buffer; neither can the single-file binary/hdf5 flavours, since the buffer path cannot express a flavour choice.
- Compiled out (no gidpost, or gidpost without zlib), the writer still exists and raises naming both `-DMESHIOPLUSPLUS_WITH_GIDPOST=ON` and `-DMESHIOPLUSPLUS_WITH_ZLIB=ON` — `.post.msh` can never silently resolve to another format.
- Meshes exceeding `INT_MAX` points, cells, or node indices raise a `WriteError` naming the count (gidpost's connectivity API is 32-bit).

## Manual validation checklist (GiD itself)

No headless GiD invocation is known to exist, so this format has no automated `cgnscheck`/`checkMesh`-style external-validation test. Everything above is verified against gidpost's own documented behaviour and by independently re-deriving the node orderings (see "Cell types"), but GiD itself is the only real oracle for whether the *application* actually opens and renders a written file correctly. A maintainer with a GiD installation should periodically check:

1. Write a mesh with mixed cell types, `point_data`, `cell_data`, and a `gmsh:physical`-tagged material id, in each of the three flavours (`mode="ascii"`, `"binary"`, `"hdf5"`).
2. Open each file in GiD's postprocessor; confirm the geometry renders correctly (no inverted/degenerate elements, correct node positions for any quadratic block present).
3. Confirm each written result (`point_data` on nodes, `cell_data` via the Gauss-point-set idiom) displays with the correct values and location.
4. Confirm the material-id column (if present) is readable as GiD material groups.
5. Confirm the provenance line appears as expected (a `# meshio++: ...` comment in ascii/binary; a `meshio++` attribute on the `Meshes/1`/`Results/...` group in hdf5 — see "Provenance" above for the known hdf5 result-side caveat).

## Notes

- Implemented entirely in the C++ core (`src/cpp/src/formats/gid.cpp`); no Python reference exists (see "Provenance" above).
- `tests/cpp/test_gid.cpp`'s `GidOrdering` suite is the written-bytes ordering oracle for `hexahedron20`/`tetra10`/`triangle6`/`quad8`; `GidWrite` covers multi-block element-id uniqueness, the empty-coordinates-block requirement, unsupported-type errors, and Gauss-point-set referencing.
- `src/cpp/third_party/gidpost/README.meshioplusplus.md` records exactly what was vendored from upstream gidpost 2.14 and why (notably, gidpost's own Fortran binding and its `cfortran.h`-dependent file are excluded — a licence carve-out, and meshio++ has its own Fortran module over its own C API regardless).
