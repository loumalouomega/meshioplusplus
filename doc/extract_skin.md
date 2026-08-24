# Skin extraction

`meshioplusplus.extract_skin(mesh, linearize=False)` derives the **boundary surface ("skin") of a 3D volume mesh** as a new surface `Mesh`. It implements the face-hashing algorithm of Kratos Multiphysics' [`SkinDetectionProcess`](https://github.com/KratosMultiphysics/Kratos): every face of every volume cell is keyed by its sorted corner-node ids, and a face whose key occurs **exactly once** is a boundary face — faces shared by two cells are interior and cancel out.

```python
import meshioplusplus

vol = meshioplusplus.read("part.msh")        # tetra/hexa/wedge/pyramid mesh
skin = meshioplusplus.extract_skin(vol)      # triangle/quad/... surface mesh
meshioplusplus.write("skin.vtu", skin)
```

The [STL](./formats/stl.md) and [PLY](./formats/ply.md) writers call this automatically (by default) when handed a volume mesh, and the [SVG](./formats/svg.md)/[TikZ](./formats/tikz.md) writers use it to render 3D meshes.

## Supported cell types

| volume type | boundary faces |
|---|---|
| `tetra` | 4 × `triangle` |
| `tetra10` | 4 × `triangle6` |
| `hexahedron` | 6 × `quad` |
| `hexahedron20` | 6 × `quad8` |
| `hexahedron27` | 6 × `quad9` |
| `wedge` | 2 × `triangle` + 3 × `quad` |
| `wedge15` | 2 × `triangle6` + 3 × `quad8` |
| `pyramid` | 4 × `triangle` + 1 × `quad` |
| `pyramid13` | 4 × `triangle6` + 1 × `quad8` |
| `pyramid14` | 4 × `triangle6` + 1 × `quad9` |

Face matching always uses the **corner nodes only**, so a mixed linear/quadratic mesh still cancels shared faces correctly. Boundary faces are wound **outward** (normals point away from the element interior), which is what makes the STL normals come out right.

## Semantics

- **Output block order** is fixed: `triangle`, `triangle6`, `quad`, `quad8`, `quad9` (only non-empty blocks appear), with faces in volume-cell enumeration order within each block.
- **Points are compacted**: only nodes referenced by a boundary face are kept, in ascending original-index order. `point_data` rows are subset the same way; `cell_data`, `field_data`, sets, and `info` are **dropped** (a boundary face has no canonical 1:1 source cell).
- **`linearize=True`** emits only the corner nodes of each face (`triangle`/`quad` output even for higher-order input); unused mid-edge/face-center nodes then compact away. This is what the STL/PLY/SVG/TikZ writers use.
- **Unsupported volume blocks** (`polyhedron*`, ragged blocks, Lagrange and very-high-order types) are skipped with a warning. If the mesh has **no** supported volume block at all, `ValueError` is raised.
- Existing 0/1/2-D blocks are ignored — the skin is derived from the volume cells only. (This also means the writers' `skin=True` mode drops pre-existing surface blocks, with a warning: a volume mesh's surface blocks are usually its boundary already, and writing both would duplicate every facet.)
- Non-conforming interfaces (a face appearing more than twice) are treated as interior — only count == 1 survives.

For the general form that also handles 2D→edge extraction and records the parent cell id, see [surface extraction](./extract_surface.md).

## Implementation

The C++ core implements the extractor in `src/cpp/src/operations/surface.cpp` — shared with [`extract_surface`](./extract_surface.md), which `extract_skin` calls in volume-only (face) mode — over the uniform mesh API (so it works under all three mesh backends), with the per-cell-type face topology tables in `src/cpp/include/meshioplusplus/detail/cell_faces.hpp`; a pure-numpy fallback (`meshioplusplus._skin._extract_skin_py`) carries a twin of the tables and produces **identical output** (same block order, connectivity, and point order — asserted by `tests/python/test_skin.py`). The public `extract_skin` tries the C++ core first and falls back to numpy, like the format shims. The face tables' outward winding is enforced by a gtest invariant (`tests/cpp/test_skin.cpp`): on each reference element, the Newell normal of every face must point away from the cell centroid.

`extract_skin` is exposed through every binding surface — Python, the C API (`mio_extract_skin`), Fortran (`mesh%extract_skin()`), and WASM (`extractSkin`) — alongside the STL/PLY writers' skin-by-default behavior.
