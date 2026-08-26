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

- **`mode`** — `"auto"` (default, inferred from the extension: `.post.bin` → binary, `.post.h5` → hdf5, anything else including `.post.msh`/`.post.res` → ascii), `"ascii"`, `"ascii_zipped"`, `"binary"`, `"hdf5"`.
- **`analysis_name`** — the GiD "analysis name" every result is grouped under (default `"meshio++"`).
- **`step`** — the single time/load step every result is written at (default `1.0`).

Three on-disk flavours:

- **ascii** — two sibling files, `<stem>.post.msh` (geometry, human-readable) + `<stem>.post.res` (results).
- **binary** — one deflated file, `<stem>.post.bin`.
- **ascii_zipped** — the same two sibling files as `ascii`, gzipped (gidpost's `GiD_PostAsciiZipped`, which is the identical text through `gzprintf`). The extension is unchanged — a gzipped file still ends `.post.msh` — so **`auto` never resolves to it**: inferring it would change what every existing `.post.msh` write produces. Reading needs no flag at all, since the reader sniffs the gzip magic.
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

### Grammar conformance

Three properties come from CIMNE's published grammar rather than from gidpost's behaviour. None is reachable through meshio++'s own writer — gidpost emits one fixed casing, always writes a mesh name, and always spells the 1-D type `Linear` — so a round trip cannot exercise any of them, and each is pinned by a hand-authored fixture instead.

- **Keywords are case-insensitive.** GiD states this explicitly for `MESH`/`dimension`/`ElemType`/`Nnode` and for the `coordinates`/`end coordinates` and `elements`/`end elements` pairs. It is not a tolerance meshio++ adds: the manual's own worked example opens a block with `Coordinates` and closes it with `end coordinates`, so a case-sensitive reader rejects the specification's own example file.
- **The mesh name is optional.** `MESH dimension 3 ElemType Linear Nnode 2` is legal (and is the manual's own second example). The name is taken only when the token after `MESH` is not the `dimension` keyword — reading it unconditionally names such a mesh `"dimension"`, which then lets a `GaussPoints` set declared `OnMesh "dimension"` bind to it and attach a result that does not belong to it.
- **The 1-D type has two spellings.** gidpost — CIMNE's own writer, vendored here — emits `Linear`, while CIMNE's current published grammar names the type `Line`. Both are read; the writer is unaffected, since it goes through gidpost.

### What reading does not recover

- **Tensor and complex result types** (`Matrix`, `MainMatrix`, `Complex*`) round-trip through the `field_data` declaration below (see [Result types](#result-types)).
- **`ResultGroup`** is read in the **ascii** flavour (see below); the binary flavour refuses it by name. `OnNurbs*` locations, mesh groups (`Group`/`End Group` — see the note under Gauss points) and range tables are skipped or refused by name rather than guessed at.

## Cell types

GiD has exactly ten element types; higher-order variants share a type with a larger node count. One meshio++ cell block becomes one named GiD mesh (`"<celltype>_<blockindex>"`, since GiD results reference meshes by name).

| meshio++ type | GiD element type | Nnode | Node permutation |
|---|---|---|---|
| `vertex` | `Point` | 1 | identity |
| `line` / `line3` | `Linear` | 2 / 3 | identity |
| `triangle` / `triangle6` | `Triangle` | 3 / 6 | identity |
| `quad` / `quad8` / `quad9` | `Quadrilateral` | 4 / 8 / 9 | identity |
| `tetra` / `tetra10` | `Tetrahedra` | 4 / 10 | identity |
| `hexahedron` / `hexahedron20` | `Hexahedra` | 8 / 20 | identity |
| `hexahedron27` | `Hexahedra` | 27 | **top ring ↔ verticals**, 3 face-centre swaps |
| `wedge` | `Prism` | 6 | identity |
| `wedge15` | `Prism` | 15 | **top triangle ↔ verticals** |
| `pyramid` | `Pyramid` | 5 | identity |
| `pyramid13` | `Pyramid` | 13 | identity |

Every type's ordering is cross-checked against Kratos Multiphysics's own geometry classes (`kratos/geometries/`). For most types this is **identity** — no node permutation — because Kratos applies no reorder at all when writing them (`triangle6`/`tetra10`/`quad8`'s edge tables match meshio++'s edge-for-edge too). `hexahedron20`, `hexahedron27` and `wedge15` are the exceptions: Kratos's own *internal* node order splits the mid-edge nodes into "bottom ring, verticals, top ring" (or the triangular-prism equivalent), while meshio++'s own table (mirroring VTK) splits them "bottom ring, top ring, verticals" — the reverse pairing of the last two blocks. `hexahedron27` additionally permutes its six face-centre nodes; its body centre and every corner/bottom-ring node are unaffected. `wedge15`'s permutation is the same two-block swap, one tier smaller. Every permutation here is **self-inverse** (an involution — `dst[c] = src[p[c]]`, the same convention `med.cpp`'s `med_node_perm()` already uses in this repo), so the identical table serves both the writer and the reader; see `gid_common.hpp`'s `gid_cell_perm_table()` for the derivation and the literal arrays. `pyramid13` needs no permutation at all — Kratos's own `Pyramid3D13` order already matches meshio++'s.

::: tip `hexahedron20` — a documentary conflict, resolved
CIMNE's own GiD 6-era figure for the 20-node hexahedron (`hexa20.gif`) numbers the mid-edge nodes **bottom ring, verticals, top ring** — that is, exactly Kratos's *internal* order, the one Kratos's own GiD writer permutes *away* from before emitting a file. Taken at face value, the figure said the identity mapping above was wrong.

Confirmed: GiD's actual expected order is the one Kratos **writes**, i.e. the post-swap order — which is meshio++'s own `hexahedron20` table. Identity is therefore correct, and the figure is outdated; CIMNE's current published grammar dropped the mid-edge figures entirely for exactly this kind of staleness, saying only *"hierarchical order … vertex nodes first, then the middle ones"* with no order given. `hexahedron8` and every lower-order type were never in question — they have no mid-edge nodes — nor were `tetra10`/`triangle6`/`quad8`/`quad9`.

A precision on the evidence, found while deriving `hexahedron27`/`wedge15`'s own orderings: the specific reorder this conclusion cites (`gid_mesh_container.h`) lives **only in that file's *Conditions*-writing path** — the Elements path (the one a volume cell type like a hexahedron actually goes through) writes no reorder at all. The conclusion is unaffected: an Element-agnostic Kratos source (`kratos/input_output/vtk_output.cpp`'s general Kratos-to-VTK conversion, which every `Hexahedra3D20` element *or* condition goes through for VTK/EnSight output, mirrored in `ensight_output.cpp`) independently reproduces the identical swap. That broader, more precise source is what `hexahedron27`/`wedge15`'s own permutations below actually lean on.
:::

The `GidOrdering` suite in `tests/cpp/test_gid.cpp` pins that the writer applies **exactly the permutation documented above** — never more, never less — by writing a cell at geometrically-known corner/edge-midpoint/face-centre/body-centre positions and re-parsing the raw file to check where each slot's value actually landed, with the expected permutation written out literally in the test rather than by calling back into the implementation's own table. For the identity types this pins that no permutation is applied; it deliberately does *not* claim to pin that meshio++'s identity mapping *is* GiD's own convention, since (absent a real GiD-written oracle) there is no way to fully rule that out. For `hexahedron27`/`wedge15`, the permutation *is* the thing being pinned, cross-checked against two independent Kratos sources per the derivation above.

`hexahedron27` and `wedge15` were previously refused in both directions as "not independently verified"; they are now supported, cross-checked against Kratos's own geometry classes as described above. Still **not supported** — throws a `WriteError`/`ReadError` naming the type, never a guess: `polygon`/`polyhedron` (GiD has no such type); every `VTK_LAGRANGE_*` and higher-degree Lagrange type.

## Geometry and ids

One shared node table is written on the **first** mesh only (`GiD_fWriteCoordinatesBlock`); every subsequent mesh gets an empty `Coordinates`/`End Coordinates` pair (gidpost's own state machine rejects a mesh with neither). Points are always 3-D, z padded with 0 for a 2-D mesh. Element ids are globally unique 1-based integers across every block — gidpost's own bulk writer numbers `1..n` *per mesh*, which would collide the moment a mesh has more than one block, so meshio++ assigns a running id across all blocks instead.

An integral `cell_data` array named `"gmsh:physical"` is written as each element's material id (`GiD_fWriteElementsIdMatBlock`) and excluded from the result output, since it is already in the geometry file. No other key is consulted for material ids.

## Data mapping

`point_data` is written `GiD_OnNodes`. `GiD_ResultLocation` has no "on cells" concept, so `cell_data` is written `GiD_OnGaussPoints` against a Gauss-point set declared per block — by default a synthetic **one-point** set named `"gp_<mesh_name>"`, the standard GiD idiom for a per-element field and how Kratos's own GiD writer represents one. An array spanning several blocks becomes several result blocks sharing one result name but different Gauss-point sets. See [Gauss points](#gauss-points) for arrays carrying more than one value per element.

`GiD_ResultType`'s valid component counts are irregular (Scalar 1; Vector 2/3/4; Matrix 3/6; MainMatrix 12; …), and gidpost does not validate an unsupported count itself — it silently emits a malformed file. With no declaration, meshio++ validates and infers: 1 component → `GiD_Scalar`; 2 or 3 → `GiD_Vector`; anything else splits into that many named `GiD_Scalar` results (`"<name>_1"` … `"<name>_k"`, recorded as a provenance note). See [Result types](#result-types) below for declaring `Matrix`/`Complex*` explicitly instead of splitting.

Named regions and multi-step results are **not carried** (each dropped with a provenance note where applicable); every write emits exactly one step (`step`). `field_data` is not *written*, but the `"gid:result_type:*"` keys in it *are* read on write — see below.

## Result types

meshio++'s `Mesh` has no way to say "this array is a symmetric tensor" or "this array is complex" — `NDArray` has neither a complex dtype nor a string dtype for a type name — so a caller **declares** it out of band, through a `field_data` entry:

```python
mesh.field_data["gid:result_type:stress"] = [2]  # GidResultType.MATRIX
```

(`meshioplusplus.gid.ResultType` is the Python `IntEnum`; `meshioplusplus.gid.RESULT_TYPE_PREFIX` is the key prefix.) `field_data` is global rather than per-location, so one key covers an array of that name wherever it appears — a name present in both `point_data` and `cell_data` is declared once, and must be a legal count for both.

All nine `GiD_ResultType`s are supported. Counts are gidpost's own (`_ResultTypeInfo`); orders are quoted from CIMNE's Customization Manual, since meshio++ stores values **verbatim in GiD's order** rather than reinterpreting them:

| Type | Legal counts | Order |
|---|---|---|
| `Scalar` | 1 | the value |
| `Vector` | 2, 3, 4 | X, Y, Z, \|V\| (4th = signed modulus) |
| `Matrix` | 3, 6 | 3: Sxx Syy Sxy — 6: Sxx Syy Szz Sxy Syz Sxz |
| `PlainDeformationMatrix` | 4 | Sxx Syy Sxy Szz |
| `MainMatrix` | 12 | Si Sii Siii, then the three eigenvectors |
| `LocalAxes` | 3 | euler_ang_1..3 |
| `ComplexScalar` | 2 | real, imag |
| `ComplexVector` | 4, 6 | **interleaved**: x_re x_im y_re y_im [z_re z_im] |
| `ComplexMatrix` | 6, 12 | **blocked**: every real, then every imaginary |

Two facts worth not rediscovering, both taken from the manual rather than assumed: **`Matrix:6` is already meshio/VTK's symmetric-tensor order** (`Sxx Syy Szz Sxy Syz Sxz`), so a stress tensor needs no permutation in either direction — a claim earlier drafts of this page stated backwards, as the stated (and false) reason for refusing `Matrix` altogether. And **`ComplexVector` interleaves real/imaginary parts per component while `ComplexMatrix` blocks them** (every real, then every imaginary) — the same family, opposite conventions, so neither may be inferred from the other.

An illegal count for a declared type is a `WriteError` naming the array, its count, and the legal counts — never a silent fallback to splitting. On read, a declaration is recorded **only when it carries information**, i.e. when it differs from what the inference above would have produced for that count — so an ordinary scalar or 2/3-component vector round trip adds no `gid:result_type:*` key, and a mesh with no declarations at all writes byte-identical output to a build without this feature.

## Gauss points

GiD attaches a per-element result to a **Gauss-point set** of G points and emits G rows per element. meshio++'s `cell_data` has no per-point-within-cell axis (the same limit MED's ELNO/ELGA documents), so a G-point, k-component array is stored **flat**:

```python
mesh.cell_data["stress"]                          # (ncells, G*k)
mesh.field_data["gid:gauss_points:stress"] = [G]  # only when G != 1
```

The row is **Gauss-point-major** — `[gp0_c0…gp0_ck-1, gp1_c0…gp1_ck-1, …]` — which is the order GiD's own `Values` rows arrive in, so neither direction re-packs anything.

`G == 1` is the historical layout and **declares nothing**, so an ordinary per-element array is a plain `(ncells, k)` and its bytes are unchanged by this mechanism's existence. The declaration is required rather than inferred because a bare `(ncells, 3)` is genuinely ambiguous: a 3-component vector at one Gauss point, or a scalar at three. It also interacts with [result types](#result-types): a type's legal component counts are checked against **k**, never `G*k`, so a `Matrix` (k=6) at G=3 is an `(n, 18)` array still validated as 6 components.

### Natural coordinates

GiD places the points itself (`Natural Coordinates: Internal`) only for specific counts per element family. Any other count must supply them explicitly, keyed by `(cell type, G)` — that being what a Gauss-point set actually depends on, so two arrays sharing a block and a G share one set rather than repeating identical coordinates:

```python
mesh.field_data["gid:gauss_coords:triangle:5"] = [...]   # G*dim doubles, point-major
```

| Cell family | Internal-valid G | Given range | `dim` |
|---|---|---|---|
| `triangle` | 1, 3, 6 | 0…1 | 2 |
| `quad` | 1, 4, 9 | −1…1 | 2 |
| `tetra` | 1, 4, 10 | 0…1 | 3 |
| `hexahedron` | 1, 8, 27 | −1…1 | 3 |
| `wedge` | 1, 6 | 0…1 | 3 |
| `pyramid` | 1, 5 | −1…1 | 3 |
| `line` | **any** | **forbidden** | — |

Higher-order variants share their base family's rules (`hexahedron27` follows `hexahedron`). Supplying coordinates for a count GiD *could* have placed is honoured — a solver's quadrature need not match GiD's — but omitting them for one it cannot is a `WriteError` naming the type and its legal counts, never a file GiD would silently reject. Line elements are the one family GiD forbids `Given` for outright, which costs nothing since they already accept any count.

Reading captures a `Given` set's coordinates into the same key, so a non-standard count read from a file can be written back; without that the writer would refuse for want of the very coordinates the file supplied. **This capture is ASCII-only**: in the binary flavour the coordinates are raw double records the string-record scanner cannot reach, and the HDF5 attribute layout is unverified. The *count* is read in all three flavours, so G>1 values round-trip everywhere — only re-writing a `Given` set read from a binary or HDF5 file is affected.

::: warning Binary flavour: component counts
Unrelated to Gauss points, but easily mistaken for them: the binary stream carries neither a row width nor a component count, so its reader can only use the declared type's canonical width (Scalar 1, Vector 3, Matrix 6). A 2-component `Vector` is legal GiD and round-trips through ASCII, but is unrecoverable from binary at **any** G, including the default 1. Pinned by `test_binary_cannot_recover_a_non_canonical_component_count`.
:::

## ResultGroup

A `ResultGroup` packs several results sharing one analysis, step and location into a single wide `Values` row:

```
ResultGroup "analysis" step OnNodes|OnGaussPoints ["gauss set"]
ResultDescription "name" Type[:N]        <- one per member
[ResultRangesTable / ComponentNames / Unit]
Values
  <id>  v v v v v v v v v v v            <- every member's components, one line
End Values                                <- there is no "End ResultGroup"
```

Each member is unpacked into an **ordinary result**, so everything else — `point_data`/`cell_data` routing, the flat `(ncells, G*k)` [Gauss-point](#gauss-points) layout, result-type recording, `time_step` selection — applies to it unchanged, with no second code path to drift.

**Member widths** come from `:N` when present (a `ResultDescription` is the one place GiD states a width; a plain `Result` header does not, which is why *that* path infers from the row). Absent it, the manual fixes the widths inside a group: a `Vector` is always 3 components and a `Matrix` always 6. A type admitting several widths with no `:N` — the `Complex*` family — is **refused by name** rather than guessed, since picking one would silently mis-associate every column after it. A row whose width disagrees with the members' total is refused for the same reason.

**Not supported in the binary or HDF5 flavours** — refused by name there. The binary path previously *skipped* the record silently, losing the whole group with no diagnostic; that is fixed. Unpacking it would need the binary record layout of `ResultDescription`/`Values`, which no available file exercises. meshio++ never *writes* a `ResultGroup` in any flavour.

::: tip Malformed results files warn rather than fail
The `.post.res` sibling is optional, so an **absent** one is normal and silent — a mesh with no results reads back as geometry only. A **malformed** one is also not fatal (the geometry is still good, and the alternative is refusing a file GiD itself opens), but it now emits a warning naming the problem. Before this, the reader's own refusals were unreachable from `read`: they were caught by the optional-sibling handler and every result vanished with no diagnostic at all — including results that had already parsed cleanly.
:::

::: warning `Group` / `OnGroup` is not a region concept
Easy to misread from the name: GiD's `Group`/`End Group` wraps the `MESH` blocks belonging to **one time step**, so GiD can swap meshes as the step changes — it exists for re-meshing and adaptive analyses. It is not a named entity set, and the postprocess format has **no node/element set concept at all**; the optional material column (already round-tripping as `gmsh:physical`) is its only grouping mechanism. Its real meshio++ analogue is the transient axis, not `Region`. `Group`/`OnGroup` are skipped on read.
:::

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
