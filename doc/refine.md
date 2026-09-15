# Refinement

`meshioplusplus.refine(mesh, levels=1)` subdivides cells into congruent **children of the same cell type**, increasing a mesh's resolution while leaving the object it describes intact. It is a mesh **operation** (like [cell conversion](/convert_cells) and [crop](/crop)), not a file format, and uses only standard C++/numpy, so it runs under every mesh backend.

With no selector it refines **every** cell (the uniform case, below). With one — a cell list, a region, or a `cell_data` threshold — it refines only those cells and resolves the hanging nodes that leaves, so the output is still a valid conforming mesh: see [Selective (adaptive) refinement](#selective-adaptive-refinement).

![Each hexahedron refined into 8, volume and orientation preserved](/images/refine_hex.png)

```python
import meshioplusplus

mesh = meshioplusplus.read("domain.msh")

# one level: every tetra -> 8 tetra, every quad -> 4 quad, ...
fine = meshioplusplus.refine(mesh)

# two levels: 64x the cells for a volume mesh, 16x for a surface mesh
finer = meshioplusplus.refine(mesh, levels=2)

# tag each output cell with the original cell it descends from
tagged = meshioplusplus.refine(mesh, record_parent_ids=True)

# adaptive: refine only the worst cells, and close up conformingly around them
graded = meshioplusplus.refine(
    meshioplusplus.attach_quality(mesh),
    where="quality:scaled_jacobian < 0.3",
    record_levels=True,
)

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

Each block is refined by its own template, so a mixed-type mesh is fine. New nodes are numbered to coincide with the type's own full-Lagrange layout (`line3`, `triangle6`, `quad9`, `tetra10`, `wedge18`, `hexahedron27`), which is what makes the templates readable against the reference elements.

A wedge refines as a *triangle 1-to-4 split × 2 vertical levels*, and the mid-level triangle's three edge midpoints **are** its three quad-face centres — which is why a wedge needs no body node while a hexahedron does.

Three constructs raise rather than guess:

- **higher-order** cells (`tetra10`, `quad8`, …) — linearize the mesh first (`convert_cells(mesh, mode="linearize")`);
- **`pyramid`**, whose uniform refinement is 6 pyramids + 4 tetrahedra and so cannot keep the same-type contract — simplexify the mesh first;
- **ragged** polygon/polyhedron blocks, which have no subdivision template.

Passing an unsupported block through unchanged would leave hanging nodes at every interface with a refined neighbour, so these fail loudly instead.

## Conformity

**The refined mesh has no hanging nodes.** Mid-edge nodes *and* quad-face-centre nodes are shared between every cell touching the entity; only the hexahedron body node is per-cell. Sharing the face centres is not optional — with a per-cell copy, two hexahedra meeting at a face would reference distinct coincident nodes and the mesh would be topologically split along every interior face, which [`extract_surface`](/extract_surface) would then report as boundary.

Two neighbours compute a shared node's coordinate as the *mean* of the same corner set. The mean is order-independent, so both arrive at bit-identical coordinates with no tie-breaking rule.

The tetrahedron's interior diagonal is fixed at the opposite-edge pair `(0,1)`–`(2,3)` for **determinism only**. Being strictly interior — no face of a tetrahedron contains both opposite edges — it never affects conformity: a face's subdivision is fixed by that face's own mid-edge nodes whatever the neighbour does. This is the opposite of `convert_cells`' hex-simplexify diagonal 0–6, whose endpoints lie on the boundary and which therefore *is* conformity-critical.

## What changes

- **Block structure is preserved 1:1**: the output has exactly as many cell blocks as the input, in the same order and of the same types. That is what keeps the `cell_data` correspondence (one array per block) trivially correct under every backend.
- **Points**: originals keep their indices; new nodes are appended.
- **`point_data`**: a new node gets the **mean of its entity's corner values** (2 for an edge, 4 for a quad face, 8 for a hexahedron body), so a linear field is interpolated exactly. The input dtype is preserved.
- **`cell_data`**: each parent's row is **replicated** to its children. `field_data` passes through.
- **`point_sets` / `cell_sets`** are remapped in the Python layer; a cell-set entry expands to the parent's children (which are contiguous).
- `record_parent_ids=True` attaches an Int64 `refine:parent_cell` `cell_data` array naming, per output cell, the **original** input cell it descends from within its own block — the original ancestor, not the immediate parent, even at several levels.
- `mesh.info` and `gmsh_periodic` are not carried through (as for every other operation).

Output is **deterministic**: the templates are fixed and the new-node numbering is assigned by a serial pass over a parallel-filled buffer, never a concurrent hash insert. Results are byte-identical across the MESHIO/NATIVE/KRATOS backends and any thread count.

## Selective (adaptive) refinement

Refining a *subset* of the cells leaves hanging nodes on the interface, so the operation also has to close the mesh back up. `refine` does that from one set — which **edges** carry a new node — and derives everything else from it.

![One hexahedron of a 4x4x4 block refined adaptively, coloured by refine:level, beside the same request under propagate](/images/refine_adaptive.png)

Refining one corner cell of a 4×4×4 block takes it from 64 to **125** cells under the default closure, against **512** — the full uniform refinement — under propagation.

### Selectors

At most **one** may be given; two is an error rather than a precedence rule, because silently ignoring a selector the caller typed is the failure mode worth refusing. With none, every cell is refined.

| selector | meaning |
|---|---|
| `cells=[3, 7, 12]` | global **block-major** cell indices — the numbering [named regions](/regions) and `partition_labels` use |
| `region="hot"` | a `Cell` region selects its own cells; a `Point` region selects every cell with **any** node in it; a `Side` region is an error, since a facet is not a cell |
| `where="quality:skewness > 0.7"` | a threshold on a scalar `cell_data` array, spelled `NAME OP VALUE` with `OP` one of `<`, `<=`, `>`, `>=`, `==`, `!=` |

The **closure** then decides what happens to everything the selection touched: `"redgreen"` (default) and `"propagate"` both return a conforming mesh; `"balanced"` keeps the hanging nodes and only enforces 2:1 balance.

A **non-finite** cell value never matches a predicate. That is deliberate rather than incidental: [`compute_quality`](/quality) reports `NaN` where a metric does not apply — `scaled_jacobian` is `NaN` for a quadrilateral in 3-D, for instance — so rejecting such an array outright would break the case the feature exists for. The flip side is that a predicate over a metric that is N/A for a cell type simply selects nothing there.

### How the closure works

![Every triangle mask has a template; a quad admits only none, an opposite pair or all four split edges, and any other mask is promoted](/diagrams/refine_masks_2d.svg)

![Representative admissible tetrahedron masks by number of split edges, and the three-edges-at-one-corner mask that is excluded on purpose](/diagrams/refine_masks_tetra.svg)

A cell's state is a bitmask over its edges. Not every mask has a same-type subdivision, so a mask that has none is **promoted** to the smallest *admissible* superset — which splits more edges, which may promote further, and so on to a fixed point. Because each type's admissible masks are closed under intersection, "the smallest admissible superset" is well defined and the promotion is a monotone idempotent closure operator; the mesh-wide fixed point is therefore unique and **independent of the order cells are visited in**. Determinism here is a property of the formulation, not a convention about traversal.

Whether a new node exists is likewise *derived*, never tabulated:

- an **edge** carries a node iff it is split;
- a **quad face** carries a centre iff **all four** of its edges are split;
- a **hexahedron** carries a body node iff **all twelve** are.

Two cells sharing an entity read the same edges and so reach the same answer, which is what makes conformity structural rather than something the tests merely sample. When every edge is split these rules collapse into the uniform templates above — which is why uniform output is byte-identical to a build without any of this.

### Support matrix

`closure="redgreen"` (the default) promotes to the smallest admissible mask:

| cell type | admissible masks | children | what propagates |
|---|---|---|---|
| `line` | all | 1, 2 | nothing |
| `triangle` | all | 1, 2, 3, 4 | nothing |
| `quad` | none, either **opposite** edge pair, all | 1, 2, 2, 4 | one row of a structured grid |
| `tetra` | none, any 1 or 2 edges, the 4 face-triples, all | 1, 2, 3, 4, 8 | nothing, for those masks |
| `wedge` | none, the 6 triangle edges, the 3 verticals, both | 1, 4, 2, 8 | one layer / one column |
| `hexahedron` | unions of the 3 parallel edge classes | 1, 2, 2, 2, 4, 4, 4, 8 | one dual sheet |

The quadrilateral row is **forced, not chosen**: a quadrangulation of an `n`-gon satisfies `4Q = B + 2I`, so an odd boundary count is impossible and a quad with one or three split edges has *no* all-quad subdivision at any number of interior nodes. Promoting it to the opposite pair is the finest type-preserving answer there is, and the reason a single refined quad costs a row rather than the whole grid. The tetrahedron's three-edges-at-a-common-vertex mask is excluded deliberately: it would put the ambiguous two-edge case on all three incident faces at once.

Every green child has the **parent's own cell type**, which is what keeps output blocks 1:1 with input blocks — and therefore keeps `mCellMaps`, `refine:parent_cell`, the `cell_data` slicing and every binding's map shape exactly as they are in the uniform case.

`closure="propagate"` instead promotes any non-empty mask straight to a full split. It is always conforming and defined for every cell type, but it is **not local**: a fully split cell splits all of its edges, so every edge-neighbour is promoted in turn and the cascade reaches the whole edge-connected component. On a connected mesh it *is* the uniform refinement. It ships as the always-works baseline and as the test oracle — not as an adaptivity mode.

### `closure="balanced"` — keep the hanging nodes

The third option does not close at all. It splits a cell fully or not at all — there are no transitional templates — and only enforces **2:1 balance**, drawing a cell in when a neighbour would otherwise end up more than one level finer:

> refine `C` ⟹ `C`'s level rises by one `D` must refine ⟺ some cell sharing a **node** with `D` would end up more than one level above `D`

Adjacency is by shared *node*, not by shared edge, and that is not a detail: across a hanging interface the coarse cell spans a whole edge while the fine cell has only half of it, so the two are *different* entities and an edge-keyed rule would be blind to exactly the coarse/fine adjacency it exists to police. (It is also the stronger, standard "corner balance", so a diagonal neighbour counts.)

On a mesh of uniform level that condition holds nowhere, so refining one cell propagates to **nothing**:

| 4×4×4 block, one cell selected | cells |
|---|---|
| `redgreen` | 125 |
| `propagate` | 512 (= uniform) |
| **`balanced`** | **71** (= 64 − 1 + 8) |

Balancing only bites once levels differ — from the second adaptive pass onwards — and even then it reaches one ring of neighbours rather than the whole mesh. This is what an adaptive-mesh-refinement code normally means by "propagate", and the only mode whose cost is bounded by the selection rather than by the mesh.

The price is stated rather than hidden: the result is **1-irregular and not conforming**. Every constrained node is reported in the Int64 `refine:hanging` `point_data` array (`1` = hanging), so a solver can eliminate it — the array marks *exactly* the constrained nodes, edge midpoints and quad-face centres alike, neither a superset nor a subset. The conformity guarantees below apply to the other two closures only, and `extract_surface`, `decimate` and anything else assuming a conforming mesh will treat a hanging node as a genuine boundary.

1-irregular is the *only* thing given up. The mesh is never **torn**: no two distinct nodes ever occupy the same position with cells referencing both. That is what [`refine:entity`](#refine-entity) is for.

`refine:level` is what makes the rule well defined across passes, which is another reason that array is *maintained* rather than replicated.

### `refine:entity` {#refine-entity}

The Int64 `(num_points, 4)` `point_data` array recording, per point, the **entity it was created on**:

| key | meaning |
|---|---|
| `(-1, -1, a, b)` | the midpoint of edge `(a, b)`, `a < b` |
| `(p, q, r, s)` sorted | the centre of the quadrilateral face with those corners |
| `(-1, -1, -1, -1)` | a point `refine` did not create — an original point, or a hexahedron body centre, which is interior to one cell and can never be shared |

It exists so that a **later** pass can recognise that an entity it is about to split already carries a node and reuse it, instead of allocating a coincident second one. Without it, a second balanced pass tore the mesh: when the 2:1 balance rule draws in a coarse cell whose edge `(a, b)` already holds a hanging node, that cell's own refinement keys the edge as `(a, b)`, allocates a fresh midpoint, and leaves the two sides of the interface referencing distinct — but exactly coincident — nodes. On a 4×4×4 hex block, `levels=2` produced 12 such positions and `levels=3` produced 96.

Like `refine:level` the name is **reserved** and the array is *maintained* rather than replicated: it is attached whenever a pass leaves hanging nodes and kept thereafter. The conforming closures leave none and attach nothing, so their output is byte-identical to what it was before the array existed.

Entries name **point indices**, which `refine` never renumbers. An operation that does renumber points (`reorder`, `clean`, `crop`, `merge`) or moves them (`transform`, `smooth`) invalidates them. `refine` detects that — every entry must still reproduce its own point's coordinates — and warns and ignores the array rather than trusting it, falling back to allocating fresh nodes.

::: tip Why the hanging rule is evaluated over the output
A node is constrained when it sits on an entity of a cell that does not reference it, and that has to be checked against the **emitted** cells rather than the input cells' entities. The two differ as soon as a pass refines more than one level: a cell the balance rule draws in is split into children whose own sub-edges the input cell never had, so a node the neighbouring refinement places on one of them is constrained without any input entity ever naming it. Checking the input left 42 of 84 constrained nodes unreported on a two-level balanced refinement of a 4×4×4 block.
:::

### One choice made from global node ids

A cell with two *adjacent* split edges on one face leaves a quadrilateral remnant there that needs a diagonal, and the cell on the other side of that face must pick the same one. The choice is therefore made from the **global node ids** — the diagonal starting at whichever of the two surviving corners has the smaller id — never from the template's local numbering, which two differently-oriented neighbours would disagree about.

### `refine:level`

`record_levels=True` attaches an Int64 `refine:level` `cell_data` array: `0` for a cell no full split ever touched, incremented once per full split. A **transitional child inherits its parent's level unchanged**, because a green split is a closure, not a refinement. Colour by it to see the grading.

The name is **reserved**. An input that already carries it is *updated* rather than replicated, whatever the flag says — the flag only controls creating it — so successive passes accumulate rather than reporting a stale depth. A mesh without the array and without the flag is unaffected.

With `levels > 1` and a selector, level *k* refines the children of level *k − 1*'s fully split cells; green and untouched cells are not re-refined.

### `refine:cell_id` and `refine:parent_id` {#refinecell_id-and-refineparent_id}

![The refine:cell_id and refine:parent_id arrays linking a fine mesh back to the coarse mesh it was refined from](/images/refine_hierarchy.png)

`record_hierarchy=True` attaches two Int64 `cell_data` arrays that persist a parent/child hierarchy the caller can resolve *across separate calls*, unlike `mCellMaps` (which dies with the call) or `refine:parent_cell` (a row index into a mesh the caller may since have renumbered).

**The hierarchy is a link between two meshes, not a tree inside one.** A multigrid — or future green-undo — caller keeps every pass's output mesh; the fine mesh's `refine:parent_id` resolves against the coarse mesh's `refine:cell_id`:

- An unsplit cell **keeps its own id** and is its own parent (`parent_id == cell_id`).
- A split cell's children each get a **fresh id** and carry the parent's id.
- When the coarse mesh carries no `refine:cell_id` of its own, a cell's implicit id is its **global block-major index** — the same numbering `detail::block_bases`/`partition_labels` use.
- `parent_id` names a cell of **the mesh passed to this call** — the same rule `refine:parent_cell` follows: two `levels=1` calls give immediate parents, one `levels=n` call gives the original ancestor, because only the mesh actually named in the call is one the caller could have kept.

Given both meshes, red/green/untouched classify with **no further storage**: an untouched cell has `parent_id == cell_id`; a green (transitional) child has `parent_id != cell_id` and the same `refine:level` as its parent; a red child's `refine:level` is one more than its parent's — which is why no third `refine:green` array exists.

`record_hierarchy=True` also **forces `refine:entity` to be attached**, even when the closure leaves no hanging node (the normal `redgreen`/`propagate` case, where it would otherwise never appear at all — see above). That array already records exactly which 2 or 4 coarse nodes each new fine node is the mean of, which is the multigrid prolongation operator's stencil and weights: combined with `mPointMap` being the identity (a coarse node keeps its index in the fine mesh), it completes the coarse→fine transfer operator for free.

The name is **reserved**, on the same policy as `refine:level`: an input that already carries a *valid* `refine:cell_id` is updated whatever the flag says. A **duplicated or malformed id** is dropped with a warning rather than trusted — ids are values, not indices, so `reorder`/`crop`/`clean` carry them correctly with no coordinate check, but a repeat means the mesh was `merge`d with another hierarchy or the array was replicated without being updated.

**Interacts with `split`:** an unqualified `split(mesh, by="tag")` (no explicit array name) auto-picks the first sorted integer `cell_data` array, and `"refine:cell_id"` sorts ahead of ordinary material tags. On a mesh carrying this array, always name the array explicitly — and note that `split(mesh, by="tag", tag="refine:parent_id")` is itself the sibling-grouping recipe: it yields one piece per distinct parent (a singleton for every untouched cell).

### What is not implemented

- **Green elements are not undone** before a later refinement. The standard rule is to restore a transitional cell to its parent and re-split from scratch; this implementation refines the transitional children directly, so *repeated* selective passes over the same region degrade element quality without bound. `refine:cell_id`/`refine:parent_id` plus `refine:level` is exactly the hierarchy a future green-undo would need, which is why they can now be recorded.
- **Mixed-type green** (a quad closing into triangles, a hex into pyramids) is not offered. It would localize the quad's one- and three-edge cases further, at the cost of one input block emitting several output blocks — which the `mCellMaps` "first child within the corresponding output block" contract, and therefore every binding's map shape, cannot express.

### Invariants

- The output is **conforming** under `redgreen` and `propagate`: no hanging nodes, and no facet shared by more than two cells. Under `balanced` it is 1-irregular by design, and `refine:hanging` marks every constrained node.
- The output is **never torn**, under any closure and any number of passes: two distinct nodes never occupy the same position with cells referencing both. See [`refine:entity`](#refine-entity).
- Every selected cell is fully split; every other cell is either untouched or minimally split by the closure.
- New nodes are the **mean of their entity's corners**, which is order-independent, so neighbours agree bit-for-bit with no tie-break.
- `cell_data` is replicated parent → children, transitional children included; `point_data` is interpolated onto the new nodes exactly as in the uniform case.
- Named regions (and so `point_sets`/`cell_sets`) are remapped, with a refined cell's region membership expanding to its children.
- Output is byte-identical across the MESHIO/NATIVE/KRATOS backends, across thread counts, and across the C++-core/numpy-fallback boundary — pinned by `tests/python/test_refine.py::test_cpp_matches_python_selective`.

## Volume and orientation

Children inherit the parent's orientation, so a well-oriented input refines to an output with **zero newly-inverted cells**.

Volume is conserved **exactly** for `line`/`triangle`/`quad`/`tetra` always, and for `wedge`/`hexahedron` when the parent is affine (a right prism / parallelepiped). For a general trilinear hexahedron the eight children's volumes do not sum to the parent's, because the parent's bilinear faces are replaced by four different bilinear patches — that is a property of the geometry, not of this implementation.

Refinement is exponential in `levels`: one level multiplies a surface mesh's cell count by 4 and a volume mesh's by 8. `refine` logs a warning when a level would produce more than ~20 M cells, since the failure mode at depth is exhausting memory rather than returning a wrong answer.

The numpy fallback handles rectangular cell blocks, which is every type `refine` supports.

## CLI

```bash
meshioplusplus refine in.msh out.vtu
meshioplusplus refine in.msh out.vtu --levels 2
meshioplusplus refine in.msh out.vtu --record-parent-ids

# adaptive
meshioplusplus refine in.msh out.vtu --cells 12,13,44 --record-levels
meshioplusplus refine in.msh out.vtu --region hot
meshioplusplus refine in.msh out.vtu --where "quality:scaled_jacobian < 0.3"
meshioplusplus refine in.msh out.vtu --cells 12 --closure propagate
meshioplusplus refine in.msh out.vtu --cells 12 --closure balanced --record-levels
meshioplusplus refine in.msh out.vtu --cells 12,13 --record-hierarchy
```

See the [CLI reference](/cli).

## Other languages

- **C API** — `mio_refine(mesh, levels, record_parent_ids)` returns an opaque `mio_refine_result` (`_mesh` borrow, `_take_mesh`, zero-copy `_point_map` / `_cell_map`, `_free`). `mio_refine_ex(mesh, &opts)` takes a `mio_refine_opts` with the selection, closure and `record_hierarchy`; zero-initialize it through `mio_refine_opts_init` and it reproduces `mio_refine` exactly. See the [C API reference](/c_api).
- **Fortran** — `mesh%refine(levels, record_parent_ids=..., point_map=..., cells=..., region=..., where_array=..., where_op=..., where_value=..., closure=..., record_levels=..., record_hierarchy=...)`; `cells` is 1-based here. See the [Fortran reference](/fortran).
- **Julia** — `refine(mesh; levels, record_parent_ids, cells, region, where_array, where_op, where_value, closure, record_levels, record_hierarchy)`; `cells` is 1-based. `refine:cell_id`/`refine:parent_id` ride the raw C-side numbering unshifted, like `partition_labels`' part ids — they are stable identifiers, not index maps. See the [Julia reference](/julia).
- **R** — `mio_refine(mesh, levels, record_parent_ids, cells, region, where_array, where_op, where_value, closure, record_levels, record_hierarchy)`; `cells` is 1-based. See the [R reference](/r).
- **WebAssembly / JavaScript** — `refine(mesh, levels, recordParentIds, options)` where `options` is `{cells, region, array, compare, value, closure, recordLevels, recordHierarchy}`, and the `convertSurfaceOps` pipeline op `{op: 'refine', ...}` takes the same fields. The comparison is spelled `compare` there because `op` is the pipeline step's own discriminant. See the [WebAssembly reference](/wasm).
