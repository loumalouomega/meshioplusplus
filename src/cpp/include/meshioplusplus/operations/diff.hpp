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
 * @file operations/diff.hpp
 * @brief Dependency-free mesh comparison ("diff"): compare two meshes and
 * report whether they are equivalent within a tolerance, plus a structured
 * description of any differences.
 *
 * `diff` compares points, cell blocks (types + counts + connectivity), and
 * `point_data`/`cell_data`/`field_data` (keys, shapes, values), and returns a
 * `DiffReport` with an overall verdict — `Identical`, `EqualWithinTolerance`,
 * or `Different` — plus a per-section breakdown (max abs/rel error, worst
 * offender index, mismatch counts, keys present in only one mesh, first
 * mismatching cells). The elementwise tolerance rule is
 * `abs_err <= atol + rtol*|expected|` with `expected` taken from the first
 * mesh.
 *
 * Two modes:
 *  - **Ordered** (default): points and cells are compared index-by-index — the
 *    fast O(N), parallel path for the common regression-test case where only
 *    numerical values may have drifted.
 *  - **Unordered** (`DiffOptions::unordered`): points are matched between meshes
 *    by spatial proximity within tolerance via a bucket-grid coordinate hash
 *    (no O(N^2) scan), then connectivity is compared up to that node
 *    correspondence. It is conservative: if a robust one-to-one correspondence
 *    cannot be established it says so (`mCorrespondenceFailed`) rather than
 *    guessing.
 *
 * Everything uses only standard C++ and the uniform mesh API, so it runs under
 * every mesh backend. This is an operation, not a file format — it is not in
 * the format registry. Named `point_sets`/`cell_sets` are not visible to the
 * C++ core; the Python shim compares those on top of this report.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Overall outcome of comparing two meshes, in increasing order of severity.
enum class DiffVerdict {
    Identical,             ///< Structurally equal and every value bitwise-equal.
    EqualWithinTolerance,  ///< Structurally equal; values differ but within tolerance.
    Different,             ///< A structural difference or a value beyond tolerance.
};

/// Human-readable name of a verdict (`"identical"` / `"equal within tolerance"`
/// / `"different"`).
MESHIOPLUSPLUS_API const std::string& diff_verdict_name(DiffVerdict verdict);

/**
 * @brief Elementwise comparison summary for one array (the points table or one
 * named data array; cell_data folds all its blocks into one summary).
 */
struct ArrayDiff {
    /// Name of the data key; empty for the points table.
    std::string mName;
    /// The two arrays could not be compared elementwise (differing shape /
    /// column count / row count with no correspondence).
    bool mShapeMismatch = false;
    /// Element counts on each side (product of shape).
    std::int64_t mSizeA = 0;
    std::int64_t mSizeB = 0;
    /// Largest absolute error `|a - b|` over all compared elements.
    double mMaxAbsError = 0.0;
    /// Largest relative error `|a - b| / |a|` (0 where `a == 0`).
    double mMaxRelError = 0.0;
    /// Flat index of the worst-offending element (largest absolute error), or
    /// -1 if nothing was compared.
    std::int64_t mWorstIndex = -1;
    /// Number of elements exceeding `atol + rtol*|a|`.
    std::int64_t mNumExceeding = 0;
    /// Every compared element was bitwise-equal.
    bool mExact = true;
};

/**
 * @brief Comparison of one cell block (matched by position across the two
 * meshes).
 */
struct BlockDiff {
    /// Block index (position in cell-block order).
    std::int64_t mBlock = -1;
    /// meshio cell-type names on each side.
    std::string mTypeA;
    std::string mTypeB;
    /// Cell counts on each side.
    std::int64_t mCountA = 0;
    std::int64_t mCountB = 0;
    /// The two blocks have different cell types.
    bool mTypeMismatch = false;
    /// The two blocks have different cell counts (or shapes not comparable).
    bool mCountMismatch = false;
    /// Cells whose connectivity differs (after any node correspondence).
    std::int64_t mConnMismatchCount = 0;
    /// The first few mismatching cell indices (capped at `max_reported`).
    std::vector<std::int64_t> mFirstMismatches;
};

/**
 * @brief Comparison of one data section (`point_data`, `cell_data`, or
 * `field_data`).
 */
struct DataDiff {
    /// Keys present in mesh A but not mesh B (sorted).
    std::vector<std::string> mOnlyInA;
    /// Keys present in mesh B but not mesh A (sorted).
    std::vector<std::string> mOnlyInB;
    /// Per shared key, its elementwise comparison summary (sorted by name).
    std::vector<ArrayDiff> mShared;
};

/**
 * @brief Comparison of two meshes' named regions (`region.hpp`).
 *
 * A region is identified by its name *and* kind, since the same name may
 * legitimately label a node set and an element set. Membership is compared
 * exactly: entries are canonical (sorted, de-duplicated), so two regions
 * describing the same group are bitwise equal and no tolerance applies.
 */
struct RegionDiff {
    /// `"name (kind)"` for every region in A but not B (sorted).
    std::vector<std::string> mOnlyInA;
    /// `"name (kind)"` for every region in B but not A (sorted).
    std::vector<std::string> mOnlyInB;
    /// `"name (kind)"` for every region present in both whose membership,
    /// dimension or tag differs (sorted).
    std::vector<std::string> mChanged;

    /// Whether any region differs at all.
    bool Differs() const { return !mOnlyInA.empty() || !mOnlyInB.empty() || !mChanged.empty(); }
};

/**
 * @brief The full structured result of `diff`.
 */
struct DiffReport {
    /// Overall verdict (the maximum severity over every section).
    DiffVerdict mVerdict = DiffVerdict::Identical;
    /// The same verdict with `mRegions` left out. `meshes_equal` uses this:
    /// that helper is documented to compare geometry and data, not named
    /// groups, and the Python shim has always recomputed a sets-agnostic
    /// verdict for exactly this reason.
    DiffVerdict mVerdictIgnoringRegions = DiffVerdict::Identical;
    /// The unordered (proximity) node-matching mode was used.
    bool mUnordered = false;
    /// Unordered mode could not build a robust one-to-one node correspondence.
    bool mCorrespondenceFailed = false;

    // --- points ---
    /// Point counts / dimensions differ, so points are not elementwise-comparable.
    bool mPointCountMismatch = false;
    std::int64_t mNumPointsA = 0;
    std::int64_t mNumPointsB = 0;
    /// Coordinate comparison (empty/`mShapeMismatch` when counts differ).
    ArrayDiff mPoints;

    // --- cells ---
    /// The two meshes have different numbers of cell blocks.
    bool mBlockCountMismatch = false;
    std::int64_t mNumBlocksA = 0;
    std::int64_t mNumBlocksB = 0;
    /// Per-block comparison, over `min(numBlocksA, numBlocksB)` blocks.
    std::vector<BlockDiff> mBlocks;

    // --- data ---
    DataDiff mPointData;
    DataDiff mCellData;
    DataDiff mFieldData;

    // --- named regions ---
    /// Named-group comparison. Before meshio++ 8.1 sets never reached the C++
    /// core, so only the Python shim could compare them; this is the core's
    /// own answer, which is what the native CLI and the flat bindings report.
    RegionDiff mRegions;

    /// Human-readable notes (e.g. "connectivity not comparable: ragged").
    std::vector<std::string> mMessages;
};

/// Tuning for `diff`.
struct DiffOptions {
    /// Absolute tolerance in `abs_err <= atol + rtol*|expected|`.
    double atol = 1e-12;
    /// Relative tolerance in `abs_err <= atol + rtol*|expected|`.
    double rtol = 1e-9;
    /// Match points by spatial proximity instead of index (see the file docs).
    bool unordered = false;
    /// Cap on the number of first-mismatching cell indices recorded per block.
    std::int64_t max_reported = 10;
};

/**
 * @brief Compare two meshes and return a structured report.
 * @param rA the first ("expected") mesh
 * @param rB the second mesh
 * @param rOpts tolerances and mode (defaults: `atol=1e-12`, `rtol=1e-9`,
 *              ordered)
 * @return the filled `DiffReport`
 */
MESHIOPLUSPLUS_API DiffReport diff(const Mesh& rA, const Mesh& rB, const DiffOptions& rOpts = {});

/**
 * @brief Convenience boolean: are two meshes equal within tolerance?
 *
 * Equivalent to `diff(rA, rB, {atol, rtol}).mVerdict != DiffVerdict::Different`
 * in the ordered mode — handy in test suites.
 * @param rA the first mesh
 * @param rB the second mesh
 * @param atol absolute tolerance (default `1e-12`)
 * @param rtol relative tolerance (default `1e-9`)
 * @return `true` when the verdict is `Identical` or `EqualWithinTolerance`
 */
MESHIOPLUSPLUS_API bool meshes_equal(const Mesh& rA, const Mesh& rB, double atol = 1e-12, double rtol = 1e-9);

}  // namespace meshioplusplus
