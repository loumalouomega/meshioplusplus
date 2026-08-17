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
 * @file undo_green.hpp
 * @brief Green-element undo: restore `refine`'s *transitional* (closure-only)
 * cells back to their original parent, so a subsequent selective `refine` call
 * re-splits the parent from scratch rather than compounding an already
 * irregular shape. See `operations/refine.hpp`'s own doc comment for the
 * red/green/untouched vocabulary this builds on.
 *
 * A TWO-MESH operation, `undo_green(coarse, fine)` — the "link between two
 * meshes, not a tree inside one" framing `refine.hpp` already uses for its
 * persistent hierarchy, taken literally: `coarse` is the mesh a prior
 * `refine(coarse, ..., record_hierarchy=True, record_levels=True)` call was
 * run on, `fine` is that call's output.
 *
 * **Design: lookup and substitution, not reconstruction.** `refine()` never
 * renumbers or prunes points -- its point map is always the identity -- so a
 * green parent's *exact* original connectivity and cell_data are already
 * sitting, byte-for-byte, in `coarse` at the row `fine`'s `refine:parent_id`
 * names (resolved against `coarse`'s own `refine:cell_id`, or its implicit
 * global-block-major id when it carries none -- the same fallback
 * `refine_attach_hierarchy` itself uses when starting a fresh id space). This
 * needs no per-type subdivision-table inversion, no graph matching against
 * `detail::refine_templates.hpp`'s admissible masks, and consequently no
 * winding repair or other discrete sign branch -- unlike `subdivide` and
 * `agglomerate`, this operation is pure array bookkeeping and byte copies, and
 * so (unlike those two) it has a full numpy twin rather than being
 * twin-exempt.
 *
 * A cell's mask -- and hence its red/green status -- is uniform across every
 * one of its children (`refine_once` sets it once per parent), so
 * classification is per SIBLING GROUP (cells sharing one `refine:parent_id`),
 * not per cell: a singleton group (`refine:cell_id == refine:parent_id`) is
 * **untouched**, kept verbatim; a group whose `refine:level` is one more than
 * its coarse parent's own level is **red** (a genuine, wanted refinement --
 * passed through unchanged); a group whose level equals its coarse parent's
 * own level is **green** (a closure artefact -- the whole group is replaced
 * by ONE cell copied verbatim from `coarse`).
 *
 * **Points are never pruned or renumbered** -- this is what makes the
 * substitution a zero-translation byte copy: a coarse cell's node ids are
 * valid indices into `fine`'s own point array with no offset, precisely
 * because `refine`'s own point map is always the identity. `clean(mesh,
 * remove_orphans=True)` is the documented follow-up for a caller wanting the
 * orphaned mid-edge nodes (left behind by a substituted green group) pruned.
 *
 * The six reserved `refine:*` cell_data arrays (`parent_cell`, `level`,
 * `hanging`, `entity`, `cell_id`, `parent_id`) are unconditionally dropped
 * from the output -- they describe a hierarchy relationship that is now stale
 * after the undo; a subsequent `refine(..., record_hierarchy=True)` call
 * rebuilds them fresh. Every other `cell_data` array on `fine` is carried:
 * untouched/red rows byte-copied from `fine` as usual, a green group's one
 * output row byte-copied from the SAME-NAMED array on `coarse` -- if `coarse`
 * lacks that array, or its shape doesn't match, the WHOLE array is dropped
 * with a warning rather than guessing a value for the substituted row.
 *
 * **Named Side regions do not survive** (the `subdivide`/`agglomerate`
 * precedent): a removed green child's local facet numbering has no
 * correspondence to the substituted parent's own facets, even though the
 * cell type is unchanged. Point and Cell regions do survive -- Cell regions
 * through the first genuinely non-injective `CellMapKind::Direct` use in the
 * repo (several fine cells collapsing onto one output row), relying on
 * `Region::Canonicalize`'s existing sort+dedup.
 *
 * **Two honest limitations, not gaps**: it can only undo the LAST generation
 * relative to the specific `coarse` mesh passed in (an untouched cell's
 * `parent_id == cell_id`, so an older green closure becomes indistinguishable
 * from an original once a later pass has run over it), and it needs the
 * caller to hold both meshes -- there is no single-mesh fallback.
 */

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/// The result of `undo_green`.
struct UndoGreenResult {
    Mesh mMesh;
    /// Per FINE input block, Int64 shape `(num_cells_in_block,)`,
    /// `CellMapKind::Direct` shape: fine cell -> output cell within the
    /// corresponding output block. Every member of a green sibling group maps
    /// to the SAME output index (the group's one substituted row).
    std::vector<NDArray> mCellMaps;
    /// Number of green sibling groups substituted (i.e. number of NEW cells
    /// created by this call).
    std::int64_t mNumGroupsUndone = 0;
    /// Total number of green child cells removed (mNumGroupsUndone counted
    /// once each; this is the sum of each undone group's original size).
    std::int64_t mNumCellsRemoved = 0;
};

/**
 * @brief Restore `fine`'s transitional (green) cells back to their original
 * parent, read verbatim from `coarse`.
 * @param rCoarse the mesh a prior `refine(rCoarse, ...)` call was run on.
 * @param rFine that call's output -- must carry one Int64 scalar `cell_data`
 *        array per block for each of `refine:cell_id`, `refine:parent_id` and
 *        `refine:level` (i.e. that call must have used
 *        `record_hierarchy=True, record_levels=True`); throws by name
 *        otherwise.
 * @return the undone mesh (fine's own block structure, unchanged types and
 *         order; some blocks' row counts shrink), plus the cell maps and
 *         counters above.
 * @throws std::invalid_argument if `rFine` lacks the required hierarchy
 *         arrays, if a `refine:parent_id` value cannot be resolved against
 *         `rCoarse`'s id space (the two meshes are not the input/output pair
 *         of one `refine()` call), or if a sibling group's `refine:level`
 *         matches neither the red nor the green relationship to its coarse
 *         parent's own level.
 */
MESHIOPLUSPLUS_API UndoGreenResult undo_green(const Mesh& rCoarse, const Mesh& rFine);

}  // namespace meshioplusplus
