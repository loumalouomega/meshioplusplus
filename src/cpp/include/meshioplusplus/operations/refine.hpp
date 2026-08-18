//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file refine.hpp
 * @brief Mesh refinement: subdivide cells into congruent children of the
 * **same** cell type, interpolating `point_data` onto the new nodes. Either
 * every cell (uniform) or a selected subset with a conforming closure
 * (selective / adaptive).
 *
 * This is the resolution-increasing counterpart to the resolution-preserving
 * `convert_cells` and the resolution-reducing `crop`/`clean`. One level applies
 * a fixed template per cell type:
 *
 *  - `line` -> 2 `line` (1 new mid-edge node)
 *  - `triangle` -> 4 `triangle` (3 mid-edge nodes; the standard 1-to-4 split)
 *  - `quad` -> 4 `quad` (4 mid-edge + 1 face-centre node)
 *  - `tetra` -> 8 `tetra` (6 mid-edge nodes)
 *  - `wedge` -> 8 `wedge` (9 mid-edge + 3 quad-face-centre nodes)
 *  - `hexahedron` -> 8 `hexahedron` (12 mid-edge + 6 face-centre + 1 body node)
 *
 * New nodes are placed at the midpoint of the entity that defines them and
 * carry the mean of that entity's corner values for every `point_data` array.
 * Original points keep their indices; new nodes are appended. Each parent's
 * `cell_data` row is replicated to its children.
 *
 * **Conformity.** Mid-edge nodes and quad-face-centre nodes are *shared*
 * between every cell that touches the entity, so the refined mesh has no
 * hanging nodes. Only the hexahedron body node is per-cell. (Sharing the face
 * centres is not optional: with a per-cell copy, two hexahedra meeting at a
 * face would reference distinct coincident nodes and the mesh would be
 * topologically split along every interior face.)
 *
 * Note that the tetrahedron's interior diagonal -- fixed at the opposite-edge
 * pair `(0,1)`-`(2,3)`, i.e. `tetra10` nodes 4-9 -- is chosen for determinism
 * only, **not** for conformity: it is strictly interior, so a face's
 * subdivision is fixed by that face's own mid-edge nodes whatever the
 * neighbour does. This is the opposite of `convert_cells`' hex-simplexify
 * diagonal 0-6, whose endpoints lie on the boundary and which therefore *is*
 * conformity-critical.
 *
 * **Volume and orientation.** Children inherit the parent's orientation, so a
 * well-oriented input refines to a well-oriented output (zero newly-inverted
 * cells). Volume is conserved exactly for `line`/`triangle`/`quad`/`tetra`
 * always, and for `wedge`/`hexahedron` when the parent is affine (a right
 * prism / parallelepiped). For a general trilinear hexahedron the eight
 * children's volumes do not sum to the parent's, because the parent's bilinear
 * faces are replaced by four different bilinear patches -- that is a property
 * of the geometry, not of this implementation.
 *
 * **Block structure is preserved 1:1**: the output has exactly
 * `NumCellBlocks()` blocks, in input order, and each keeps its input cell
 * type, which is what keeps the one-array-per-block `cell_data` invariant
 * trivially correct. That holds in selective mode too -- see below.
 *
 * ### Selective (adaptive) refinement
 *
 * With a selector set (`mCells`, `mRegion` or the `mPredicateArray` predicate)
 * only the selected cells get the full template above -- they are *red* -- and
 * the hanging nodes that leaves on the interface are resolved by a **closure**,
 * so the output is still a valid conforming mesh. Setting no selector is the
 * uniform behaviour, byte-identical to a build without this feature.
 *
 * The whole thing is one table-driven rule rather than three loosely coupled
 * features. Each cell type has a set of **admissible** split-edge masks -- the
 * subsets of its edges for which a same-type subdivision template exists (see
 * `detail/refine_templates.hpp`). A cell's mask is whatever its neighbours have
 * bisected; if that mask is not admissible it is *promoted* to the smallest
 * admissible superset, which bisects more edges, which may promote further. The
 * admissible sets are closed under intersection, so "smallest admissible
 * superset" is well defined and the promotion is a monotone idempotent closure
 * operator -- and the global fixed point is therefore **unique and independent
 * of iteration order**, which is why determinism across backends and thread
 * counts is a property of the formulation rather than a convention.
 *
 * Every new node's existence is likewise *derived* from the split-edge set,
 * never tabulated: an edge carries a node iff it is split, a quad face carries
 * a centre iff **all four** of its edges are split, and a hexahedron carries a
 * body node iff **all twelve** are. Two cells sharing an entity see the same
 * edges and so reach the same answer, which is what makes conformity structural
 * rather than something the tests merely sample. When every edge is split these
 * rules reduce exactly to the uniform templates above.
 *
 * `RefineClosure::RedGreen` (the default) promotes to the smallest admissible
 * superset, which keeps the extra refinement local. Per type:
 *
 *  - `line`, `triangle`: every mask is admissible -- a triangle with 1 bisected
 *    edge splits into 2, with 2 into 3, with 3 into 4. Nothing propagates.
 *  - `quad`: the admissible masks are the two *opposite* edge pairs and the
 *    full split, so a single bisected edge promotes to its opposite pair and
 *    splits into 2 quads. A pentagon or heptagon cannot be partitioned into
 *    quadrilaterals at all (`4Q = B + 2I` forbids an odd boundary count), so
 *    this is the finest possible type-preserving answer; the bisection then
 *    travels along one row of a structured grid and stops at the boundary.
 *  - `hexahedron`: the admissible masks are unions of its three parallel edge
 *    classes -- 1, 2, 2, 2, 4, 4, 4 or 8 children -- so refinement propagates
 *    through one dual sheet rather than the whole block.
 *  - `wedge`: its six triangle edges form one class and its three verticals
 *    another -- 1, 2, 4 or 8 children.
 *  - `tetra`: every mask up to two edges is admissible, as are the four
 *    face-triples; anything else promotes to the full 8-way split. (Three
 *    edges meeting at a common vertex are deliberately *not* admissible: each
 *    of the three incident faces would then hit the ambiguous two-edge case.)
 *
 * `RefineClosure::Propagate` instead promotes any non-empty mask straight to a
 * full split. It is always conforming and defined for every cell type, but it
 * is **not local** -- every edge-neighbour of a red cell becomes red in turn,
 * so it converges to uniform refinement of the whole edge-connected component.
 * It is the always-works baseline and the test oracle, not the adaptivity mode.
 *
 * `RefineClosure::Balanced` does not close at all: it **keeps the hanging
 * nodes** and only enforces 2:1 balance, which is what an adaptive-mesh-
 * refinement code normally means by "propagate". A cell is split fully or not
 * at all -- there are no transitional templates -- and a cell is drawn in only
 * when a neighbour would otherwise end up more than one level finer than it:
 *
 *     refine C  =>  C's level rises by one
 *     D must refine  <=>  some entity D shares has an incident cell whose
 *                         post-refinement level exceeds D's by more than one
 *
 * On a mesh of uniform level that condition is satisfied nowhere, so refining
 * one cell propagates to **nothing** -- 64 hexahedra become 71, against 125
 * under `RedGreen` and 512 under `Propagate`. Balancing only bites once levels
 * differ, i.e. from the second adaptive pass onwards, and even then it reaches
 * one level-ring rather than the whole mesh. The `refine:level` array is what
 * makes that well defined across passes, and reading it back is why the array
 * is *maintained* rather than replicated.
 *
 * The price is stated rather than hidden: the result is **1-irregular and not
 * conforming**. Every constrained node is reported in the `refine:hanging`
 * `point_data` array (see `kRefineHangingName`) so a solver can eliminate it;
 * the conformity guarantees below apply to the other two closures only, and
 * `extract_surface`, `decimate` and anything else assuming a conforming mesh
 * will treat a hanging node as a genuine boundary.
 *
 * A cell with two *adjacent* bisected edges on one face has a remnant
 * quadrilateral there that needs a diagonal, and the neighbour across that face
 * must choose the same one. The choice is therefore made from the **global node
 * ids** -- the diagonal starting at whichever of the two surviving corners has
 * the smaller id -- never from the template's local numbering, which two
 * differently-oriented neighbours would disagree about.
 *
 * Green cells are **not undone** before a later red refinement, so repeated
 * selective passes over the same region degrade element quality without bound.
 * `refine:level` plus `mCellMaps` is the hierarchy a future green-undo needs.
 *
 * Constructs that would break the same-type contract raise rather than guess:
 * higher-order cells (`tetra10`, ...; linearize first), `pyramid` (whose
 * uniform refinement is 6 pyramids + 4 tetrahedra; simplexify first), and
 * ragged polygon/polyhedron blocks.
 *
 * Determinism: the templates are fixed and the new-node numbering comes from a
 * **serial** dedup pass over a `parallel_for`-filled disjoint-slot buffer
 * (`src/cpp/src/operations/surface.cpp`'s phase-split idiom), never from a
 * concurrent hash insert. Output is byte-identical across mesh backends and
 * thread counts.
 *
 * Standard C++ and the uniform mesh API only, so it compiles under every mesh
 * backend. This is an operation, not a file format -- it is not in the format
 * registry.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// The `cell_data` array `RefineOptions::mRecordParentIds` attaches.
inline constexpr const char* kRefineParentCellName = "refine:parent_cell";

/// The Int64 `cell_data` array recording each cell's refinement depth: `0` for
/// a cell no red split ever touched, incremented once per red split. A green
/// (transitional) child inherits its parent's level unchanged, because a green
/// split is a closure, not a refinement. The name is **reserved**: if the input
/// already carries it, `refine` updates it rather than replicating it, so
/// successive passes accumulate.
inline constexpr const char* kRefineLevelName = "refine:level";

/// The Int64 `point_data` array `RefineClosure::Balanced` attaches: `1` for a
/// **hanging** (constrained) node, `0` otherwise. A hanging node is one that
/// exists on an entity of a cell that does not reference it -- the mid-edge node
/// a refined cell created on an edge its unrefined neighbour still spans whole.
/// Only `Balanced` produces any; the other closures leave none by construction
/// and do not attach the array. The name is **reserved**: the array is recomputed
/// from the output's own topology every pass rather than carried through the
/// generic `point_data` path, which would interpolate a flag.
inline constexpr const char* kRefineHangingName = "refine:hanging";

/// The Int64 `(num_points, 4)` `point_data` array recording, per point, the
/// **entity it was created on** -- `{-1, -1, a, b}` for the midpoint of edge
/// `(a, b)`, the sorted `{p, q, r, s}` for a quad-face centre, and
/// `{-1, -1, -1, -1}` for a point `refine` did not create (an original point, or
/// a hexahedron body centre, which is interior to one cell and can never be
/// shared).
///
/// It exists so that a *later* pass can recognise that an entity it is about to
/// split already carries a node, and reuse it instead of allocating a coincident
/// second one. Without it, `RefineClosure::Balanced` tore the mesh: when the 2:1
/// balance rule drew in a coarse cell whose edge `(a, b)` already held a hanging
/// node, that cell's refinement keyed the edge as `(a, b)`, allocated a *new*
/// midpoint and left the two sides of the interface referencing distinct,
/// exactly-coincident nodes.
///
/// The name is **reserved** and the array is *maintained* rather than replicated:
/// it is attached whenever a pass leaves hanging nodes, and kept thereafter. The
/// other closures leave none and do not attach it, so their output is unchanged.
///
/// Entries name **input point indices**, which `refine` never renumbers. An
/// operation that *does* renumber points (`reorder`, `clean`, `crop`, `merge`) or
/// that moves them (`transform`, `smooth`) invalidates them; `refine` detects
/// that -- each entry must still reproduce its own point's coordinates -- and
/// warns and ignores the array rather than trusting it.
inline constexpr const char* kRefineEntityName = "refine:entity";

/// The Int64 `cell_data` array `RefineOptions::mRecordHierarchy` attaches: a
/// stable, never-reused id per cell. Together with `kRefineParentIdName` this is
/// the persistent parent/child hierarchy `refine:level` alone cannot give --
/// **a link between two meshes, not a tree inside one**. A multigrid or
/// green-undo caller keeps every pass's output mesh; the fine mesh's
/// `refine:parent_id` resolves against the coarse mesh's `refine:cell_id`, so no
/// ancestry path or packed location code is needed. An untouched cell KEEPS its
/// id (see `kRefineParentIdName`).
///
/// The name is **reserved**: an input that already carries it is *updated*
/// rather than replicated whatever the flag says -- the `kRefineLevelName` rule
/// -- so ids stay unique across repeated calls. When the input carries none,
/// ids are implicit: a cell's id is its global block-major index (the same
/// numbering `detail::block_bases`/`partition_labels` use). Computed once, in a
/// single serial pass over the FINAL mesh after every internal level has run --
/// like `kRefineParentCellName`, never per intermediate level -- so it costs
/// nothing when `mLevels == 1` and needs no per-pass bookkeeping at all: an
/// untouched original cell (one whose composed cell map run has length one) is
/// its own parent, and every other cell in a split run gets one fresh id from a
/// single counter walking the final mesh's own block/cell order.
///
/// A multigrid interpolation stencil rides alongside for free: setting this
/// flag also forces `kRefineEntityName` to be attached even when no pass left a
/// hanging node (the `redgreen`/`propagate` closures normally never attach it),
/// since it already records exactly which coarse nodes each new fine node is
/// the corner mean of -- the prolongation operator's weights.
///
/// **Interacts with `split`:** an unqualified `split(mesh, by="tag")` (no
/// explicit array name) auto-picks the first sorted integer `cell_data` array,
/// and `"refine:cell_id"` sorts ahead of ordinary material tags. On a mesh
/// carrying this array, always name the array explicitly
/// (`split(mesh, by="tag", tag="my_material")`).
inline constexpr const char* kRefineCellIdName = "refine:cell_id";

/// The Int64 `cell_data` array naming, per output cell, the `kRefineCellIdName`
/// of the cell in the mesh **passed to this call** that it came from -- the
/// same rule `kRefineParentCellName` follows: two `levels=1` calls give
/// immediate parents (the mesh the caller passed to the second call), one
/// `levels=n` call gives the original ancestor (the mesh the caller passed to
/// that one call), because only the input actually named in the call is a mesh
/// the caller could have kept. An untouched cell is its own parent
/// (`parent_id == cell_id`); every cell of a split run carries the run's
/// original-in-this-call ancestor's id.
///
/// Given both meshes, red/green/untouched classify with no further storage:
/// untouched has `parent_id == cell_id`; a green (transitional) child has
/// `parent_id != cell_id` and the same `refine:level` as its parent; a red
/// child's `refine:level` is one more than its parent's. See `doc/refine.md`.
///
/// Attaching a duplicate `kRefineCellIdName` value (e.g. by `merge`-ing two
/// hierarchied meshes) is detected and both reserved arrays are dropped with a
/// warning, never trusted -- the `kRefineEntityName` staleness policy.
inline constexpr const char* kRefineParentIdName = "refine:parent_id";

/// How `refine` resolves the hanging nodes a partial refinement leaves behind.
enum class RefineClosure {
    /// Promote a cell's split-edge mask to the smallest *admissible* superset,
    /// so an affected neighbour is split transitionally rather than fully. Keeps
    /// the extra refinement local, and the output is conforming. The default.
    RedGreen = 0,
    /// Promote any non-empty mask straight to a full split. Conforming and
    /// defined for every cell type, but **not local**: it converges to uniform
    /// refinement of the whole edge-connected component.
    Propagate = 1,
    /// Do not close at all: **keep the hanging nodes** and merely enforce 2:1
    /// balance, refining a cell only when a neighbour would otherwise end up
    /// more than one level finer. The output is 1-irregular and **NOT
    /// conforming** -- the constrained nodes are reported in `refine:hanging`
    /// for a solver to eliminate. This is the classic adaptive-mesh-refinement
    /// meaning of "propagate", and the only mode whose cost is bounded by the
    /// selection rather than by the mesh.
    Balanced = 2,
};

/// The comparison in `RefineOptions`' `cell_data` predicate selector.
enum class RefineCompare {
    Less = 0,
    LessEqual = 1,
    Greater = 2,
    GreaterEqual = 3,
    Equal = 4,
    NotEqual = 5,
};

/**
 * @brief Parse a closure name: `"redgreen"` / `"red-green"` / `"green"`,
 * `"propagate"` / `"red"`, or `"balanced"` / `"2:1"`.
 * @param rName The name; empty means the default (`RedGreen`).
 * @throws std::invalid_argument naming every accepted value.
 */
MESHIOPLUSPLUS_API RefineClosure refine_closure_from_name(const std::string& rName);

/**
 * @brief Parse a comparison operator (`"<"`, `"<="`, `">"`, `">="`, `"=="`,
 * `"!="`; `"="` is accepted as `"=="`).
 * @throws std::invalid_argument naming every accepted value.
 */
MESHIOPLUSPLUS_API RefineCompare refine_compare_from_name(const std::string& rName);

/**
 * @brief Evaluate one `RefineCompare` against a value.
 *
 * Exposed because `crop_predicate` selects cells by exactly this rule, and a
 * second transcription of a discrete branch is precisely the kind of thing that
 * drifts silently -- the two operations would then disagree only on the boundary
 * cases (`==` on a tie, and a NaN), which is where disagreement is hardest to
 * notice.
 *
 * @param Value the cell's value.
 * @param Op the comparison.
 * @param Rhs the threshold.
 * @return whether the comparison holds. **A non-finite @p Value never matches**,
 *   whatever @p Op is: `compute_quality` deliberately reports NaN where a metric
 *   does not apply, and a predicate over `quality:*` on a mixed mesh is the
 *   headline use case, so rejecting such an array outright would break it.
 */
MESHIOPLUSPLUS_API bool refine_compare_value(double Value, RefineCompare Op, double Rhs);

/// Options for `refine`.
struct RefineOptions {
    /// How many times to apply the subdivision templates. `0` (or less) returns
    /// an unchanged clone; `n` multiplies the cell count of a supported block
    /// by `children_per_cell^n` when no selector is set. With a selector, level
    /// `k > 1` refines the children of level `k - 1`'s red cells; green and
    /// untouched cells are not re-refined.
    int mLevels = 1;
    /// Attach an Int64 `refine:parent_cell` `cell_data` array recording, per
    /// output cell, the index of the **original** input cell it descends from
    /// *within its own block* (blocks correspond 1:1). Across several levels
    /// this is the original ancestor, not the immediate parent.
    bool mRecordParentIds = false;

    // --- selective refinement ------------------------------------------------
    // At most ONE of the three selectors below may be set. Two is an error
    // rather than a precedence rule: silently ignoring a selector the caller
    // asked for is the failure mode worth refusing. All empty = uniform.

    /// Explicit **global block-major** cell indices to refine (the numbering
    /// `detail/cell_index.hpp` owns and named `Cell` regions use). Canonicalized
    /// (sorted, de-duplicated) before use; an out-of-range index is an error.
    std::vector<std::int64_t> mCells;
    /// Name of a region to refine. A `Cell` region selects its own cells; a
    /// `Point` region selects every cell with **any** node in it. A `Side`
    /// region is not a selector and is an error -- it names facets, and turning
    /// facets into cells is a policy decision this operation does not make
    /// silently.
    std::string mRegion;
    /// Name of a scalar numeric `cell_data` array to threshold (`""` = unused).
    /// Composes directly with `attach_quality`, e.g.
    /// `quality:scaled_jacobian < 0.3`. Deliberately a single comparison and
    /// not a second `data_calc` expression grammar.
    std::string mPredicateArray;
    /// The predicate's comparison.
    RefineCompare mPredicateOp = RefineCompare::Less;
    /// The predicate's right-hand side. A non-finite cell value never matches.
    double mPredicateValue = 0.0;
    /// How to resolve hanging nodes. Ignored when no selector is set (every
    /// cell is then red and no closure is needed).
    RefineClosure mClosure = RefineClosure::RedGreen;
    /// Attach the Int64 `refine:level` `cell_data` array (see
    /// `kRefineLevelName`). An input that already carries it is updated
    /// whatever this flag says; the flag only controls *creating* it.
    bool mRecordLevels = false;
    /// Attach the `refine:cell_id`/`refine:parent_id` `cell_data` arrays (see
    /// `kRefineCellIdName`/`kRefineParentIdName`) -- the persistent parent/child
    /// hierarchy a multigrid caller (or a future green-undo) resolves across
    /// the sequence of meshes it keeps. An input that already carries
    /// `refine:cell_id` is updated whatever this flag says; the flag only
    /// controls *creating* it. Also forces `refine:entity` (`kRefineEntityName`)
    /// to be attached even when the closure leaves no hanging node, since it
    /// already records the coarse corners each new fine node is the mean of --
    /// the prolongation operator's weights, which `redgreen`/`propagate` would
    /// otherwise never expose.
    bool mRecordHierarchy = false;
};

/// The result of `refine`: the refined mesh plus the index maps.
struct RefineResult {
    /// The refined mesh.
    Mesh mMesh;
    /// Int64 shape `(num_points_in,)`, input point index -> output point index.
    /// Refinement never prunes, so this is the identity; it is returned in full
    /// so callers need not depend on that.
    NDArray mPointMap;
    /// Per input block, Int64 shape `(num_cells_in_block,)`, input cell -> the
    /// index of its **first** child in the corresponding output block. A cell's
    /// children are contiguous, so cell `c` owns
    /// `[map[c], c + 1 < n ? map[c + 1] : num_cells_out)`. In selective mode a
    /// cell may have a single child (itself, unrefined), so the run length
    /// varies -- the map is still monotone and every entry non-negative.
    std::vector<NDArray> mCellMaps;
};

/**
 * @brief Refine a mesh, subdividing cells into same-type children.
 * @param rMesh The mesh to refine (unchanged).
 * @param rOptions Level count, cell selection, closure and bookkeeping. With no
 *   selector set this is the uniform refinement of every cell.
 * @return The refined mesh and its point/cell index maps.
 * @throws std::invalid_argument on a higher-order, `pyramid`, or ragged cell
 *   block, none of which can be subdivided into same-type children; on more
 *   than one selector being set; on an out-of-range cell index; on an unknown
 *   region, or a `Side` region used as a selector; and on a predicate array
 *   that is missing, does not cover every block, is not scalar, or holds a
 *   non-finite value.
 */
MESHIOPLUSPLUS_API RefineResult refine(const Mesh& rMesh, const RefineOptions& rOptions = {});

}  // namespace meshioplusplus
