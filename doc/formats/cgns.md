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

Permutations are derived from the SIDS element-numbering-conventions edge/face tables, cross-checked against VTK's own CGNS translator tables (`vtkCGNSReaderInternal.cxx`) and, for the hexahedron/wedge families, against the real `vtkWedge`/`vtkTriQuadraticHexahedron` node-ordering source directly — not against meshio++'s own `detail/cell_faces.cpp` (used for `extract_surface`/skin, a different concern with its own pre-existing, separately-tracked `hexahedron27` face-centre defect that this work found but did not fix, to avoid an unrelated risk to the `refine` operation's internal indexing scheme; see that file's comment). Two findings worth recording:

- **`wedge` is identity** — meshio++'s wedge node order (0,1,2 bottom triangle; 3,4,5 top triangle) already matches SIDS's `PENTA_6` face-for-face. Do **not** apply the `{0,2,1,3,5,4}` flip `vtk_common.hpp`'s `meshio_to_vtk_order("wedge")` uses for VTK/VTU interop — that permutation exists for a *different* target convention (`vtkWedge`) and is wrong for CGNS.
- **`hexahedron27`**'s mid-face order follows the real `vtkTriQuadraticHexahedron::Faces` table (`20=(0,4,7,3)`, `21=(1,2,6,5)`, `22=(0,1,5,4)`, `23=(3,7,6,2)`, `24`=bottom, `25`=top, `26`=body centre) mapped onto the SIDS face order `F1..F6 = (0,3,2,1)(0,1,5,4)(1,2,6,5)(2,3,7,6)(0,4,7,3)(4,5,6,7)`.

**Deliberately unsupported, by name, with a distinct error rather than a guessed ordering:** the cubic/quartic Lagrange families — `line4`/`BAR_4`, `triangle10`, `quad16`, `tetra20`, `hexahedron64`, and similar — have a CGNS `ElementType_t` code but their SIDS interior-node order is shown only in the (image-only) SIDS convention figures, with no text source to verify against; VTK's own translator tables for these target `VTK_LAGRANGE_*`, which is not the ordering meshio++'s same-named types use. Writing one raises `WriteError` naming the type ("its CGNS node ordering is not yet verified in meshio++; refusing to write a guessed ordering"); reading a section with one of these `ElementType_t` codes raises the equivalent `ReadError`. Also unsupported: `polygon`/`polygon2`/`polyhedron*` (ragged blocks — CGNS's `MIXED`/`NGON_n`/`NFACE_n` sections are not written by meshio++, and are a `ReadError` naming the code on read) and any type with no CGNS `ElementType_t` at all.

## Data mapping

None — no point_data, cell_data, or field_data is read or written by this format.

## Quirks & limitations

- **Backward compatible with the pre-v9.8.0 layout.** The reader's spec/legacy discriminator is *structural*: a root child group whose `label` attribute is `"CGNSBase_t"` selects the new path; its absence falls back to a near-verbatim copy of the old reader (same messages, same `tetra`-only, `Base`/`Zone1`/`GridElements` layout) — the old writer emitted no node attributes at all, so the absence is unambiguous. This also covers files written by upstream `meshio`, whose CGNS writer uses the same historical layout. Nothing writes the legacy layout anymore.
- Indices are 1-based in the file; `-1`/`+1` conversions are applied on read/write. Integer width is `I4` when the source connectivity is 32-bit and the point count fits `int32`, else `I8`.
- Multiple `Unstructured` zones under one base are concatenated (points offset, blocks stay per-zone-per-section); a `Structured` zone raises `ReadError` rather than being silently dropped (the exact failure mode this rewrite otherwise removes); more than one `CGNSBase_t` node reads the first and `log::warn`s about the rest.
- `CGNSLibraryVersion` is written as `3.1` — a judgement call, not a spec requirement: cgnslib rejects a file whose declared version exceeds the reading library's, so writing a low, self-consistent number (3.1 is when HDF5 became CGNS's default storage and `" hdf5version"` replaced the older `" version"` root dataset) is the compatible direction.
- `" hdf5version"` records the **linked** HDF5 library's own version string (`H5get_libversion`) — this is why a byte-for-byte comparison between the C++ and h5py writers' output must exclude that one dataset (they routinely link different HDF5 builds); every other node is asserted byte-identical.
- This remains a comparatively small subset of full CGNS/SIDS: no `FlowSolution_t` (point/cell data), no `BC_t` (boundary conditions), no structured zones, no `ZoneGridConnectivity_t`. What SIDS-compliance buys is that the geometry+topology this format *does* write is genuinely readable by other tools, not that the format is feature-complete.

## Notes

- Read/written through the C++ core when built with `MESHIO_WITH_HDF5`, otherwise through `h5py` — the Python module is a **structural twin** of the C++ writer (same node tree, same attribute encoding, same type/permutation table), not merely another implementation of the same idea; the two are compared byte-for-byte (except `" hdf5version"`) in `tests/python/test_cgns.py::test_structural_parity_with_cpp`.
- **External validation, and its honest limits.** No `cgnslib`/`cgnscheck` binary was available in the environment this rewrite was developed in, so there is currently no reference `.cgns` file (from real cgnslib) committed under `tests/python/meshes/cgns/`, and no `cgnscheck`-gated test. What CI *can* verify: internal round-trip fidelity (write→read reproduces the input, across every supported type, multi-block, and 2-D points), that the two writers (C++/Python) agree byte-for-byte, and — the test that actually matters for the node-ordering tables — a **geometric** check (`tests/cpp/test_cgns.cpp`'s `CgnsOrdering` suite) that builds a reference element with nodes at their true SIDS-defined positions (corner/edge-midpoint/face-centre coordinates computed independently of the permutation tables) and asserts the raw file bytes place each node where SIDS says it belongs, plus a pure-computation check that every permutation table entry is a genuine mathematical involution. What CI **cannot** verify: that ParaView/cgnslib/VTK actually accept and render the output, or catch a semantically-plausible-but-wrong node ordering for a type no geometric test covers. Adding the missing reference-fixture and `cgnscheck` layers (committing a tiny cgnslib-produced `.cgns` under Git LFS, and an opt-in `cgnscheck`-gated test) is a documented follow-up, not done in this pass.
