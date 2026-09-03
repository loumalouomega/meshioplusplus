# Green-element undo (undo_green)

`meshioplusplus.undo_green(coarse, fine)` restores [refine](/refine)'s **transitional ("green") cells** back to their original parent cell — the standard rule for selective refinement: before a new refinement pass touches a region that was already closed up, restore the transitional cell to its parent and re-split from scratch. `refine` refines the transitional children directly instead, so repeated selective passes over the same region degrade element quality without bound; `undo_green` is the missing half.

It is a mesh **operation**, not a file format, uses only standard C++/numpy, and runs under every mesh backend — and, like [interpolate](/interpolate), it is the repo's second **two-mesh** operation: `coarse` is the mesh a prior `refine(coarse, ..., record_hierarchy=True, record_levels=True)` call was run on, `fine` is that call's output.

```python
import meshioplusplus as mp

coarse = mp.read("bracket.msh")
fine = mp.refine(
    coarse, cells=[123], record_hierarchy=True, record_levels=True
)
# ... later, decide to refine a DIFFERENT / larger region ...
undone = mp.undo_green(coarse, fine)
redone = mp.refine(undone, cells=[456, 789], record_hierarchy=True, record_levels=True)
```

## Design: lookup and substitution, not reconstruction

The obvious approach — invert `refine`'s per-type subdivision tables against a green group's children to reconstruct the parent's connectivity — is genuinely hard: it needs graph-isomorphism-style matching against a few dozen `(type, mask, variant)` templates. It turns out to be unnecessary. `refine()` **never renumbers or prunes points** — its point map is always the identity — so a green parent's *exact* original connectivity and cell_data are already sitting, byte-for-byte, in the `coarse` mesh at the row `fine`'s `refine:parent_id` names (resolved against `coarse`'s own `refine:cell_id`, or its implicit global-block-major id when it carries none — the same fallback [`refine`'s persistent hierarchy](/refine#refinecell_id-and-refineparent_id) itself uses when starting a fresh id space). So "undo" is a **lookup and substitution**: no template inversion, no winding repair, no discrete sign branch anywhere in it. Unlike [subdivide](/subdivide) and [agglomerate](/agglomerate), this operation has a **full numpy twin**.

## Classification: red, green, untouched

![A coarse grid, the fine mesh a selective refinement produces with its red and green cells, and the result of undo_green collapsing the green groups](/diagrams/undo_green_timeline.svg)

A cell's split mask — and hence its red/green status — is uniform across every one of its children (`refine` sets it once per parent), so classification is per **sibling group** (cells sharing one `refine:parent_id`), not per cell:

- a singleton group (`refine:cell_id == refine:parent_id`) is **untouched** — kept verbatim;
- a group whose `refine:level` is one more than its coarse parent's own level is **red** — a genuine, wanted refinement, passed through unchanged;
- a group whose level equals its coarse parent's own level is **green** — a closure artefact, replaced by **one** cell copied verbatim from `coarse`.

## Preconditions

`fine` must carry one Int64 scalar `cell_data` array per block for each of `refine:cell_id`, `refine:parent_id` and `refine:level` — i.e. it must come from a `refine()` call with `record_hierarchy=True, record_levels=True` (both; `record_hierarchy` alone does not force `record_levels`, they are independent flags). `coarse` must be the exact mesh that call was run on (or an equally-shaped mesh sharing its point/id space); a `refine:parent_id` that does not resolve against `coarse` fails by name rather than guessing, as does a sibling group whose level matches neither the red nor the green relationship to its coarse parent's own level.

## Output structure

`fine`'s own block structure is preserved exactly — same types, same order — unlike `subdivide`/`agglomerate`, which restructure into polyhedron blocks. A green substitution keeps the parent's cell type (green children always share their parent's type), so it only shrinks a block's row count, never its type or position.

Points are **never pruned or renumbered**: this is what makes the substitution a zero-translation row copy in the first place, since a coarse cell's node ids are already valid indices into `fine`'s own point array. [`clean(mesh, remove_orphans=True)`](/clean) is the documented follow-up for a caller wanting the orphaned mid-edge nodes (left behind by a substituted green group) pruned.

The six reserved `refine:*` arrays (`refine:parent_cell`, `refine:level`, `refine:cell_id`, `refine:parent_id` in `cell_data`; `refine:entity`, `refine:hanging` in `point_data`) are **unconditionally dropped** from the output — they describe a hierarchy relationship that is now stale; a subsequent `refine(..., record_hierarchy=True)` call rebuilds them fresh. Every other `cell_data` array on `fine` is carried: untouched/red rows copied from `fine` as usual, a green group's one output row copied from the **same-named array on `coarse`** — if `coarse` lacks that array, or its shape/dtype does not match, the whole array is dropped with a warning rather than guessing a value for the substituted row. `point_data` and `field_data` carry through from `fine` unchanged.

## Regions

Named **Side** regions do not survive at all (the `subdivide`/`agglomerate` precedent): a removed green child's local facet numbering has no correspondence to the substituted parent's own facets, even though the cell type is unchanged. Point regions survive trivially (points are never renumbered). Cell regions survive through the first genuinely **non-injective** `CellMapKind::Direct` remap in the C++ core — several fine cells collapsing onto the same output row — deduplicated the same way `Region`'s entries are always sorted and de-duplicated.

## Limitations

Two honest limitations, not gaps:

- it can only undo the **last generation** relative to the specific `coarse` mesh passed in — the hierarchy sets `parent_id == cell_id` for untouched cells, so an older green closure becomes indistinguishable from an original once a later pass has run over it;
- it needs the caller to hold **both** meshes — there is no single-mesh fallback.

It also supports only a **single-pass** (`levels=1`) hierarchy: a sibling group whose level is more than one deeper than its coarse parent's own (possible after an `mLevels > 1` call, where different branches of one original parent's descendant tree can reach different final depths) is refused by name rather than silently mishandled.

## CLI

```sh
meshioplusplus undo-green COARSE FINE OUT
```

Available in both the Python CLI and the native `meshioplusplus` binary — see the [CLI reference](/cli#meshioplusplus-undo-green).

## Other languages

The operation is exposed on every binding surface — module-level (not type-bound), the same shape [`interpolate`](/interpolate) uses for its own two-mesh input:

```c
/* C API: plain mesh + nullable counters, like mio_smooth */
int64_t num_groups_undone, num_cells_removed;
mio_mesh* undone = mio_undo_green(coarse, fine, &num_groups_undone, &num_cells_removed);
```

```fortran
! Fortran: module-level function (two-mesh input, like mio_interpolate)
undone = mio_undo_green(coarse, fine, num_groups_undone=n, stat=st)
```

```julia
u = undo_green(coarse, fine)
u.mesh, u.num_groups_undone, u.num_cells_removed
```

```r
u <- mio_undo_green(coarse, fine)
u$mesh; u$num_groups_undone; u$num_cells_removed
```

```js
const u = m.undoGreen(coarse, fine);
// u.mesh, u.numGroupsUndone, u.numCellsRemoved
```

Being a two-mesh operation, `undo_green` is deliberately **not** reachable as a settings-pipeline step or a `convertSurfaceOps` chain entry — the same exclusion `Merge`/`Interpolate`/`Split`/`Diff` already have, naming the `undo-green` CLI verb instead.
