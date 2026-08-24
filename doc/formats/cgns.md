# CGNS (`.cgns`)

The [CGNS](https://cgns.github.io/) (CFD General Notation System) format: since v9.8.0, a genuine **CGNS/SIDS-compliant unstructured-mesh subset** stored in HDF5 (the real ADF-over-HDF5 mapping every CGNS tool uses), readable by cgnslib/ParaView/VTK. Before v9.8.0 this was a private, tetrahedra-only encoding whose writer emitted **only the first `tetra` block it found** — a mesh with any other cell type (including every surface/2-D mesh) wrote a file with empty `ElementRange`/`ElementConnectivity` groups, unreadable even by this library's own reader.

| | |
|---|---|
| **Format name** | `cgns` |
| **Extensions** | `.cgns` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | `h5py` (or a C++ build with HDF5) |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("mesh.cgns")
meshioplusplus.cgns.write("out.cgns", mesh, compression="gzip", compression_opts=4)
```

- **`compression`** / **`compression_opts`** — HDF5 gzip filter and level.

## File structure

**The `" data"` leading-space dataset name is not an ad hoc convention** — it is cgnslib's own ADF-over-HDF5 mapping (`#define D_DATA " data"` in `ADFH.c`); every node with a payload stores it in a child dataset by that exact name. This library's pre-v9.8.0 docs claimed otherwise; that was wrong.

```
/                                    # name="HDF5 MotherNode" label="Root Node of HDF5 File" type=MT
  " format"                          # int8[15] "IEEE_LITTLE_32\0"
  " hdf5version"                     # int8[33] "HDF5 Version x.y.z" zero-padded
  CGNSLibraryVersion/                # CGNSLibraryVersion_t, R4, " data"=float32[1]={3.1}
  Base/                              # CGNSBase_t, I4, " data"=int32[2]={CellDim, PhysDim}
    Zone1/                           # Zone_t, I4/I8, " data" shape (3,1)={NVertex, NCell, 0}
      ZoneType/                      # ZoneType_t, C1, " data"=int8[12] "Unstructured" (no NUL)
      GridCoordinates/               # GridCoordinates_t, MT
        CoordinateX/" data"          # DataArray_t, R4/R8, float[NVertex]
        CoordinateY/" data"
        CoordinateZ/" data"          # only when PhysDim == 3
      <CGNS type>_<n>/               # Elements_t, one section PER MESHIO++ CELL BLOCK
        " data"                      # int32[2] = {ElementType_t code, 0}
        ElementRange/" data"         # IndexRange_t, {first, last} 1-based inclusive
        ElementConnectivity/" data"  # DataArray_t, flat 1-based, ROW-major (element-major)
      FlowSolution/                  # FlowSolution_t, point_data (see "Data mapping")
      FlowSolutionCells/             # FlowSolution_t, cell_data
```

Every node carries CGNS's four standard attributes — `name`/`label`/`type` (fixed-length 33/33/3-byte NULTERM strings) and `flags` (a 1-element int32 array, not a scalar) — and the file and every group are created with HDF5 **link and attribute creation-order tracking**. This is load-bearing, not cosmetic: cgnslib's `has_child`/`has_data` (`ADFH.c`) iterate creation order with **no name-order fallback**, so a file written with the library's un-tracked defaults is structurally invisible to a real CGNS reader even though every node is individually well-formed. The Python writer sets this via h5py's `track_order=True` on `File(...)`/`create_group(...)`; the C++ writer sets `H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED` on the FCPL/GCPL directly.

Unlike MED, **sections of the same type are not consolidated** — every meshio++ cell block becomes its own `Elements_t` section, in mesh order, named `"<CGNS type>_<n>"` with `n` a 1-based per-type counter (so `tri_quad_mesh`'s three blocks become `TRI_3_1`, `QUAD_4_1`, `TRI_3_2`, keeping block count and order exact on read).

`CellDim` is the max topological dimension over the mesh's cell blocks (or `PhysDim` for a point cloud); `PhysDim` is `clamp(max(PointDim(), CellDim), 1, 3)`. `NCell` in the Zone's `" data"` counts only cells at `CellDim` — boundary/lower-dimensional sections are not zone cells, per SIDS. A 2-D-points mesh (`PhysDim == 2`) writes only `CoordinateX`/`CoordinateY`; reading back honours however many of `CoordinateX/Y/Z` a file actually has (minimum 2), so point shape round-trips exactly rather than always coming back `(n, 3)`.

## Cell types

| meshio++ | CGNS | code | permutation |
|---|---|---|---|
| `vertex` | `NODE` | 2 | identity |
| `line` / `line3` | `BAR_2` / `BAR_3` | 3 / 4 | identity |
| `triangle` / `triangle6` | `TRI_3` / `TRI_6` | 5 / 6 | identity |
| `quad` / `quad8` / `quad9` | `QUAD_4` / `QUAD_8` / `QUAD_9` | 7 / 8 / 9 | identity |
| `tetra` / `tetra10` | `TETRA_4` / `TETRA_10` | 10 / 11 | identity |
| `pyramid` / `pyramid13` / `pyramid14` | `PYRA_5` / `PYRA_13` / `PYRA_14` | 12 / **21** / 13 | identity |
| `wedge` | `PENTA_6` | 14 | identity |
| `wedge15` | `PENTA_15` | 15 | `{0..8,12,13,14,9,10,11}` |
| `wedge18` | `PENTA_18` | 16 | `{0..8,12,13,14,9,10,11,15,16,17}` |
| `hexahedron` | `HEXA_8` | 17 | identity |
| `hexahedron20` | `HEXA_20` | 18 | `{0..11,16,17,18,19,12,13,14,15}` |
| `hexahedron27` | `HEXA_27` | 19 | `{0..11,16,17,18,19,12,13,14,15,24,22,21,23,20,25,26}` |

`PYRA_13`'s code (21) is non-monotonic — appended after `MIXED` (20) in the real CGNS `ElementType_t` enum, not a typo here.

Permutations are derived from the SIDS element-numbering-conventions edge/face tables, cross-checked against VTK's own CGNS translator tables (`vtkCGNSReaderInternal.cxx`) and, for the hexahedron/wedge families, against the real `vtkWedge`/`vtkTriQuadraticHexahedron` node-ordering source directly — not against meshio++'s own `detail/cell_faces.cpp`, which at the time carried a `hexahedron27` face-centre defect this work found (its 20/22/23 were a permuted 3-cycle of the real VTK table). That defect was **fixed in v9.9.0**, in lockstep with `cell_refine_quad_faces` and the `refine` templates plus their Python twins, so `cell_faces.cpp` is now on the same convention as this table. Two findings worth recording:

- **`wedge` is identity** — meshio++'s wedge node order (0,1,2 bottom triangle; 3,4,5 top triangle) already matches SIDS's `PENTA_6` face-for-face. Do **not** apply the `{0,2,1,3,5,4}` flip `vtk_common.hpp`'s `meshio_to_vtk_order("wedge")` uses for VTK/VTU interop — that permutation exists for a *different* target convention (`vtkWedge`) and is wrong for CGNS.
- **`hexahedron27`**'s mid-face order follows the real `vtkTriQuadraticHexahedron::Faces` table (`20=(0,4,7,3)`, `21=(1,2,6,5)`, `22=(0,1,5,4)`, `23=(3,7,6,2)`, `24`=bottom, `25`=top, `26`=body centre) mapped onto the SIDS face order `F1..F6 = (0,3,2,1)(0,1,5,4)(1,2,6,5)(2,3,7,6)(0,4,7,3)(4,5,6,7)`.

**Deliberately unsupported, by name, with a distinct error rather than a guessed ordering:** the cubic/quartic Lagrange families — `line4`/`BAR_4`, `triangle10`, `quad16`, `tetra20`, `hexahedron64`, and similar — have a CGNS `ElementType_t` code but their SIDS interior-node order is shown only in the (image-only) SIDS convention figures, with no text source to verify against; VTK's own translator tables for these target `VTK_LAGRANGE_*`, which is not the ordering meshio++'s same-named types use. Writing one raises `WriteError` naming the type ("its CGNS node ordering is not yet verified in meshio++; refusing to write a guessed ordering"); reading a section with one of these `ElementType_t` codes raises the equivalent `ReadError`. Also unsupported: `MIXED` sections, and any type with no CGNS `ElementType_t` at all. (`polygon*`/`polyhedron*` **are** supported since v9.21.0 — see [Polyhedral cells](#polyhedral-cells-ngon_n-nface_n) below.)

## Data mapping

Since v9.9.0, `point_data` and `cell_data` round-trip through `FlowSolution_t` nodes; before that a CGNS export silently dropped every field.

| meshio++ | CGNS |
|---|---|
| `point_data` | `FlowSolution` (`FlowSolution_t`), `GridLocation` = `Vertex` |
| `cell_data` | `FlowSolutionCells` (`FlowSolution_t`), `GridLocation` = `CellCenter` |
| `field_data` | — (neither per-vertex nor per-cell; not written) |

Each array becomes one `DataArray_t` child, `R8`, with one value per vertex or per zone cell. On read, `GridLocation` is honoured (absent ⇒ `Vertex`, the SIDS default) and an array whose length disagrees with `NVertex`/`NCell` is a `ReadError` naming it.

**Multi-component arrays are split, because CGNS has no component concept.** There is no `NumberOfComponents` anywhere in the SIDS — one `DataArray_t` is one scalar. A k-component meshio++ array is therefore written as k sibling nodes named `<name>_0 .. <name>_{k-1}` and re-joined on read from a *contiguous* `0..k-1` run. This is a documented meshio++ convention, not something SIDS specifies (the same class of deliberate extension as `zstd` for VTU). Anything that is not a contiguous run — a lone `foo_7`, a gap — stays a scalar under its own literal name, since guessing would invent components.

```
Zone1/
  FlowSolution/            FlowSolution_t, MT
    GridLocation/" data"   int8[6]  "Vertex"        (no trailing NUL)
    temperature/" data"    R8, float64[NVertex]     (scalar: its own name)
    velocity_0/" data"     R8, float64[NVertex]     (k=3: split into three)
    velocity_1/" data"
    velocity_2/" data"
  FlowSolutionCells/       FlowSolution_t, MT
    GridLocation/" data"   int8[10] "CellCenter"
    material/" data"       R8, float64[NCell]
```

**`cell_data` needs a single-dimension mesh.** A `CellCenter` array is per-*zone*, but meshio++'s `cell_data` is per-*block* — and only blocks at `CellDim` are zone cells. Concatenating block-major is therefore only well defined when every block is at `CellDim`; for a mixed-dimension mesh (tets plus boundary triangles, say) there is no way to distribute the zone-wide array back across the blocks on read without inventing values, so `cell_data` is skipped with a warning rather than written wrongly. Geometry is unaffected. On read the array is split back by each section's cell count in `ElementRange` order.

**`FlowSolution_t` is only read for a single-zone file.** Across several zones the arrays would have to be concatenated in whatever order the zones happen to be listed in, and a solution present on only some zones has no defensible filler; a multi-zone file's solutions are skipped with a warning. meshio++'s own writer always emits one zone.

## Polyhedral cells (`NGON_n` / `NFACE_n`)

Since v9.21.0 meshio++ reads and writes CGNS's face-based sections itself, with no optional dependency — so polyhedral CGNS works in the default build, the PyPI wheels and the WASM artifact. A `polyhedron<N>` block becomes an `NGON_n` face list plus an `NFACE_n` cell list of **signed** face element ids (the sign meaning "traverse this face reversed"); a jagged `polygon<N>` block is itself a face list, so it becomes an `NGON_n` on its own.

Three rules are worth knowing:

- **`NGON_n` is written before the `NFACE_n` sections that reference it**, and its faces are deduplicated across the *polyhedral* blocks only. A mesh mixing hexahedra with polyhedra keeps ordinary `HEXA_8` sections for the former — putting their faces in the pool would leave `NGON_n` elements that no cell references.
- **On read, an `NGON_n` becomes `polygon<N>` cells only when no `NFACE_n` references it.** Otherwise it is the shared face pool, not a set of cells, and emitting it as cells would duplicate every polyhedron's geometry.
- **Only the CGNS ≥ 4.0 layout is read** (`ElementStartOffset` beside `ElementConnectivity`). A CGNS 3.x file prefixes each row with its own length inline instead; normalising the two is exactly what `cg_poly_elements_read` exists for, so such a file is refused **by name**, pointing at the [cgnslib backend](#the-optional-cgnslib-backend), rather than misread. A file meshio++ writes with a face-based section declares `CGNSLibraryVersion` 4.0 for the same reason — below that, cgnslib itself reads `NGON_n` the 3.x way.

::: warning No pure-Python path
The h5py twin deliberately does not implement this; it raises naming the compiled core. The writer deduplicates faces and repairs each cell's winding, and that repair is a discrete branch on the sign of an enclosed volume — two implementations of such a branch can land on opposite sides for a near-degenerate cell and then diverge macroscopically, the same reasoning that keeps `smooth`'s inversion guard out of its numpy fallback. `_core` ships in every wheel, so this is not a practical restriction.
:::

## The optional cgnslib backend

Since v9.18.0 meshio++ can additionally read CGNS through the **official CGNS library** (cgnslib, the Mid-Level Library), behind `-DMESHIOPLUSPLUS_WITH_CGNSLIB=ON`. It is **OFF by default and bring-your-own** (`CGNS_ROOT`) — never vendored, never downloaded, exactly the policy [KaHIP](/partition) follows. cgnslib is Zlib-licensed and ships its own CMake config package, so no Find module is needed and the exported meshio++ package stays relocatable.

The backend is **additive**: everything above still works without it, byte for byte. What it adds is two things the raw-HDF5 reader cannot have at all:

- **ADF-container files.** `.cgns` has two on-disk containers, HDF5 and ADF. The reader above speaks HDF5 directly, so an ADF file is unreachable by construction rather than merely unimplemented — and much of the real-world corpus is ADF.
- **`NGON_n` / `NFACE_n` polyhedral sections.** `NGON_n` lists faces as node lists; `NFACE_n` lists each cell as **signed** face ids, where the sign means "traverse this face reversed" (CGNS's way of orienting a shared face outward from each of the two cells that use it). They map to `polyhedron<N>` blocks, grouped by unique node count like the OpenFOAM and EnSight readers. `cg_poly_elements_read` also absorbs the CGNS 3.x-vs-4.0 `ElementStartOffset` split, which is the most error-prone part of the encoding and the single strongest reason to use the MLL here.

**Routing.** Read goes through the MLL whenever it is built — the input is not ours and the MLL is strictly more capable — with one narrow exception: the pre-v9.8.0 legacy layout, which has no ADF node attributes and which the MLL rejects, falls through to the hand-rolled path. That is a *specific* fallback, not a blanket catch, so a genuine MLL error still surfaces. **Write is untouched** and stays on the hand-rolled path: switching engines for meshes it already handles would churn bytes for no benefit and would cost the C++/Python byte-parity oracle. A useful consequence is that on a cgnslib build the whole CGNS test suite becomes a cross-engine check for free.

**The Python reference reader is h5py-based** and so has neither capability. It is not a fallback for them either: `meshioplusplus.cgns.read` re-raises rather than falling back when the file is not HDF5, because the reference reader would report a confusing signature error instead of the real one.

`_core.__has_cgnslib__` (Python) and `cgns_has_cgnslib()` (C++) report whether the backend is present; `read_cgns_mll` always exists and throws naming the CMake flag when it is not.

### Getting cgnslib

It is plain C with a CMake build and works on every platform meshio++ targets, Windows/MSVC included — there is nothing Unix-specific about it.

| | |
| --- | --- |
| vcpkg | the `cgnslib` feature of this port, which depends on upstream's `cgns` port — the usual route on Windows |
| Conan | the `with_cgnslib` option. cgnslib is **not** on ConanCenter, so this only flips the CMake flag; supply an install and set `CGNS_ROOT` (the same status KaHIP has) |
| Debian / Ubuntu | `libcgns-dev` |
| conda-forge | `cgns` — note the package omits the `libz.so` dev symlink its own exported CMake target names, and its cgnslib is linked against conda's HDF5, so take **both** from the same prefix |
| From source | `cmake -DCGNS_ENABLE_HDF5=ON`; point meshio++ at it with `CGNS_ROOT` or `CGNS_DIR` |

meshio++'s own Windows CI leg builds every native path off and exercises the Python fallbacks, so no optional C++ dependency is tested there — that is a CI policy of this repository, not a limitation of cgnslib.

## Quirks & limitations

- **Backward compatible with the pre-v9.8.0 layout.** The reader's spec/legacy discriminator is *structural*: a root child group whose `label` attribute is `"CGNSBase_t"` selects the new path; its absence falls back to a near-verbatim copy of the old reader (same messages, same `tetra`-only, `Base`/`Zone1`/`GridElements` layout) — the old writer emitted no node attributes at all, so the absence is unambiguous. This also covers files written by upstream `meshio`, whose CGNS writer uses the same historical layout. Nothing writes the legacy layout anymore.
- Indices are 1-based in the file; `-1`/`+1` conversions are applied on read/write. Integer width is `I4` when the source connectivity is 32-bit and the point count fits `int32`, else `I8`.
- Multiple `Unstructured` zones under one base are concatenated (points offset, blocks stay per-zone-per-section); a `Structured` zone raises `ReadError` rather than being silently dropped (the exact failure mode this rewrite otherwise removes); more than one `CGNSBase_t` node reads the first and `log::warn`s about the rest.
- `CGNSLibraryVersion` is written as `3.1` — a judgement call, not a spec requirement: cgnslib rejects a file whose declared version exceeds the reading library's, so writing a low, self-consistent number (3.1 is when HDF5 became CGNS's default storage and `" hdf5version"` replaced the older `" version"` root dataset) is the compatible direction.
- `" hdf5version"` records the **linked** HDF5 library's own version string (`H5get_libversion`) — this is why a byte-for-byte comparison between the C++ and h5py writers' output must exclude that one dataset (they routinely link different HDF5 builds); every other node is asserted byte-identical.
- This remains a subset of full CGNS/SIDS: `FlowSolution_t` covers `Vertex` and `CellCenter` only (no `IFaceCenter`/`EdgeCenter`, no `Rind` planes, no `DataClass_t`/`DimensionalUnits_t` metadata), and there is no `BC_t` (boundary conditions), no structured zones, no `ZoneGridConnectivity_t`, no `Family_t`. What SIDS-compliance buys is that what this format *does* write is genuinely readable by other tools, not that the format is feature-complete.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_HDF5`, otherwise through `h5py` — the Python module is a **structural twin** of the C++ writer (same node tree, same attribute encoding, same type/permutation table), not merely another implementation of the same idea; the two are compared byte-for-byte (except `" hdf5version"`) in `tests/python/test_cgns.py::test_structural_parity_with_cpp`.
- **External validation.** Four layers, added in v9.9.0 (the v9.8.0 rewrite had only the first two, which prove self-consistency and nothing more):
  1. **Internal round-trip**: write→read reproduces the input across every supported type, multi-block, 2-D points and field data.
  2. **Cross-writer parity**: the C++ and Python writers are compared byte-for-byte (excluding `" hdf5version"`, which records the *linked* HDF5 version and legitimately differs).
  3. **Geometric ordering oracle** (`tests/cpp/test_cgns.cpp`'s `CgnsOrdering` suite): builds a reference element with nodes at their true SIDS-defined positions — computed independently of the permutation tables — and asserts the raw written bytes place each node where SIDS says it belongs, plus a pure-computation check that every table entry is a genuine involution. This is what catches a wrong permutation; a round trip through only our own reader and writer cannot, since a self-inverse permutation makes `read(write(m)) == m` even when the table is wrong.
  4. **Real cgnslib**: a reference `.cgns` written end to end by **cgnslib 4.5.2** is committed under [`tests/python/meshes/cgns/`](../../tests/python/meshes/cgns/README.md) (Git LFS) and read unconditionally by both readers; and `test_cgnscheck_accepts_our_output` runs cgnslib's own conformance checker over what we write for every supported type. `cgnscheck` reports **zero errors**. It does report style *warnings* — no `Family_t` on the zone, no `DataClass_t` on the arrays, and "not a CGNS data-name identifier" for any field whose name is not one of SIDS's standard names (`Density`, `VelocityX`, …) — which are recommendations, not conformance failures; meshio++ preserves the caller's own field names rather than renaming them to fit SIDS's vocabulary.

Layer 4 is **opt-in**: `cgnscheck` is not a pip/apt dependency of this repo, so that test skips (with an explicit, actionable reason naming the conda-forge install) where the binary is absent — including in CI. The committed fixture in the same layer needs nothing external and always runs. What still cannot be verified anywhere: that ParaView actually *renders* the output, and a semantically-plausible-but-wrong ordering for a type no geometric test covers.
