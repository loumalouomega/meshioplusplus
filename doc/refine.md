# Uniform refinement

`meshioplusplus.refine(mesh, levels=1)` subdivides **every cell into congruent
children of the same cell type**, increasing a mesh's resolution while leaving
the object it describes intact. It is a mesh **operation** (like
[cell conversion](/convert_cells) and [crop](/crop)), not a file format, and uses
only standard C++/numpy, so it runs under every mesh backend.

```python
import meshioplusplus

mesh = meshioplusplus.read("domain.msh")

# one level: every tetra -> 8 tetra, every quad -> 4 quad, ...
fine = meshioplusplus.refine(mesh)

# two levels: 64x the cells for a volume mesh, 16x for a surface mesh
finer = meshioplusplus.refine(mesh, levels=2)

# tag each output cell with the original cell it descends from
tagged = meshioplusplus.refine(mesh, record_parent_ids=True)

meshioplusplus.write("fine.vtu", fine)
```

## Templates

| cell type | children | new nodes per cell |
|---|---|---|
| `line` | 2 × `line` | 1 edge midpoint |
| `triangle` | 4 × `triangle` | 3 edge midpoints (the standard 1-to-4 split) |
| `quad` | 4 × `quad` | 4 edge midpoints + 1 face centre |
| `tetra` | 8 × `tetra` | 6 edge midpoints |
| `wedge` | 8 × `wedge` | 9 edge midpoints + 3 quad-face centres |
| `hexahedron` | 8 × `hexahedron` | 12 edge midpoints + 6 face centres + 1 body centre |

Each block is refined by its own template, so a mixed-type mesh is fine. New
nodes are numbered to coincide with the type's own full-Lagrange layout
(`line3`, `triangle6`, `quad9`, `tetra10`, `wedge18`, `hexahedron27`), which is
what makes the templates readable against the reference elements.

A wedge refines as a *triangle 1-to-4 split × 2 vertical levels*, and the
mid-level triangle's three edge midpoints **are** its three quad-face centres —
which is why a wedge needs no body node while a hexahedron does.

Three constructs raise rather than guess:

- **higher-order** cells (`tetra10`, `quad8`, …) — linearize the mesh first
  (`convert_cells(mesh, mode="linearize")`);
- **`pyramid`**, whose uniform refinement is 6 pyramids + 4 tetrahedra and so
  cannot keep the same-type contract — simplexify the mesh first;
- **ragged** polygon/polyhedron blocks, which have no subdivision template.

Passing an unsupported block through unchanged would leave hanging nodes at
every interface with a refined neighbour, so these fail loudly instead.

## Conformity

**The refined mesh has no hanging nodes.** Mid-edge nodes *and* quad-face-centre
nodes are shared between every cell touching the entity; only the hexahedron body
node is per-cell. Sharing the face centres is not optional — with a per-cell copy,
two hexahedra meeting at a face would reference distinct coincident nodes and the
mesh would be topologically split along every interior face, which
[`extract_surface`](/extract_surface) would then report as boundary.

Two neighbours compute a shared node's coordinate as the *mean* of the same
corner set. The mean is order-independent, so both arrive at bit-identical
coordinates with no tie-breaking rule.

The tetrahedron's interior diagonal is fixed at the opposite-edge pair
`(0,1)`–`(2,3)` for **determinism only**. Being strictly interior — no face of a
tetrahedron contains both opposite edges — it never affects conformity: a face's
subdivision is fixed by that face's own mid-edge nodes whatever the neighbour
does. This is the opposite of `convert_cells`' hex-simplexify diagonal 0–6, whose
endpoints lie on the boundary and which therefore *is* conformity-critical.

## What changes

- **Block structure is preserved 1:1**: the output has exactly as many cell
  blocks as the input, in the same order and of the same types. That is what
  keeps the `cell_data` correspondence (one array per block) trivially correct
  under every backend.
- **Points**: originals keep their indices; new nodes are appended.
- **`point_data`**: a new node gets the **mean of its entity's corner values**
  (2 for an edge, 4 for a quad face, 8 for a hexahedron body), so a linear field
  is interpolated exactly. The input dtype is preserved.
- **`cell_data`**: each parent's row is **replicated** to its children.
  `field_data` passes through.
- **`point_sets` / `cell_sets`** are remapped in the Python layer; a cell-set
  entry expands to the parent's children (which are contiguous).
- `record_parent_ids=True` attaches an Int64 `refine:parent_cell` `cell_data`
  array naming, per output cell, the **original** input cell it descends from
  within its own block — the original ancestor, not the immediate parent, even
  at several levels.
- `mesh.info` and `gmsh_periodic` are not carried through (as for every other
  operation).

Output is **deterministic**: the templates are fixed and the new-node numbering
is assigned by a serial pass over a parallel-filled buffer, never a concurrent
hash insert. Results are byte-identical across the MESHIO/NATIVE/KRATOS backends
and any thread count.

## Volume and orientation

Children inherit the parent's orientation, so a well-oriented input refines to an
output with **zero newly-inverted cells**.

Volume is conserved **exactly** for `line`/`triangle`/`quad`/`tetra` always, and
for `wedge`/`hexahedron` when the parent is affine (a right prism /
parallelepiped). For a general trilinear hexahedron the eight children's volumes
do not sum to the parent's, because the parent's bilinear faces are replaced by
four different bilinear patches — that is a property of the geometry, not of this
implementation.

Refinement is exponential in `levels`: one level multiplies a surface mesh's cell
count by 4 and a volume mesh's by 8. `refine` logs a warning when a level would
produce more than ~20 M cells, since the failure mode at depth is exhausting
memory rather than returning a wrong answer.

The numpy fallback handles rectangular cell blocks, which is every type `refine`
supports.

## CLI

```bash
meshioplusplus refine in.msh out.vtu
meshioplusplus refine in.msh out.vtu --levels 2
meshioplusplus refine in.msh out.vtu --record-parent-ids
```

See the [CLI reference](/cli).

## Other languages

- **C API** — `mio_refine(mesh, levels, record_parent_ids)` returns an opaque
  `mio_refine_result` (`_mesh` borrow, `_take_mesh`, zero-copy `_point_map` /
  `_cell_map`, `_free`). See the [C API reference](/c_api).
- **Fortran** — `mesh%refine(levels, record_parent_ids=..., point_map=...)`. See
  the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `refine(mesh, levels, recordParentIds)`. See the
  [WebAssembly reference](/wasm).
