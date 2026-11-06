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
 * @file detail/decimate_shared.hpp
 * @brief The collapse-loop and emission plumbing surface `decimate` and volume
 * `decimate_volume` share, beyond the installed `decimate_common.hpp`.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI. Every body was
 * moved verbatim out of the two operations -- same floating-point expressions,
 * same evaluation order, same iteration order -- so their outputs (pinned by
 * the determinism digests) are unchanged. Where the two copies differed only
 * in corners per cell (3 for a triangle, 4 for a tet), the helper is a
 * template on that count.
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/decimate_common.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace detail {

/// A Float64 working copy of a float-kind point_data array, blended in place
/// at each collapse commit.
struct DecimFloatData {
    std::string mName;
    std::size_t mNumComponents = 0;
    std::vector<double> mValues;  // n * mNumComponents
};

/// A vector of Int64 as a 1-D NDArray.
NDArray decim_int64_vector(const std::vector<std::int64_t>& rValues);

/// Float64 working copies of `rMesh`'s float-kind point_data with `n` rows,
/// in point-data name order.
std::vector<DecimFloatData> decim_float_working_copies(const Mesh& rMesh, std::size_t n);

/// The pin mask: caller's `rFrozen` | `rBoundary` (under `PreserveBoundary`)
/// | crease endpoints of the triangle mesh `rCorners`/`rNormals` (under
/// `PreserveFeatures`).
std::vector<std::uint8_t> decim_pin_mask(const std::vector<std::uint8_t>& rFrozen,
                                         const std::vector<std::uint8_t>& rBoundary,
                                         bool PreserveBoundary, bool PreserveFeatures,
                                         double FeatureAngleDeg,
                                         const std::vector<std::int64_t>& rCorners,
                                         const std::vector<double>& rNormals, std::size_t n);

/**
 * @brief The normal-flip / degenerate-survivor sweep over the faces incident
 * to `a` or `b` (the sorted rows `rFa`, `rFb`), in ascending face index.
 *
 * Faces in the sorted `rShared` die in the collapse and impose no constraint;
 * an already degenerate face (zero normal) imposes none either. `rMerged` is
 * scratch storage. Returns `true` when some surviving face would flip.
 */
bool decim_any_normal_flip(const std::vector<double>& rXyz,
                           const std::vector<std::int64_t>& rCorners,
                           const std::vector<std::int64_t>& rFa,
                           const std::vector<std::int64_t>& rFb,
                           const std::vector<std::int64_t>& rShared, std::int64_t a, std::int64_t b,
                           const double* pPlaced, std::vector<std::int64_t>& rMerged);

/// Commit's FP part: blends the float data of survivor `s` toward removed `r`
/// by the placed point's projection onto the edge, moves `s` to `pPlaced`,
/// and adds `r`'s quadric into `s`'s.
void decim_commit_blend(std::vector<double>& rXyz, std::vector<double>& rQuads,
                        std::vector<DecimFloatData>& rFloatData, const double* pPlaced,
                        std::int64_t s, std::int64_t r);

/// The used input points' working coordinates, cast back once to `rPoints`'
/// dtype.
NDArray decim_emit_points(const NDArray& rPoints, const std::vector<double>& rXyz,
                          const std::vector<char>& rUsed, std::size_t n, std::size_t num_used,
                          std::size_t dim);

/// Copies `rSrc`'s cell_data (alive rows only, `rBlockBase` giving each
/// block's first global cell id; `pSkipName`, when non-null, is not
/// forwarded), point_data (float-kind from `rFloatData`, others as the used
/// rows) and field_data into `rOut`.
void decim_emit_data(const Mesh& rSrc, const std::vector<std::uint8_t>& rAlive,
                     const std::vector<std::int64_t>& rBlockBase, const char* pSkipName,
                     const std::vector<DecimFloatData>& rFloatData, const std::vector<char>& rUsed,
                     std::size_t n, std::size_t num_used, Mesh& rOut);

/// Every input point's survivor's output index (follows `rSuccessor` chains).
NDArray decim_point_map(const std::vector<std::int64_t>& rSuccessor,
                        const std::vector<std::int64_t>& rRemap, std::size_t n);

/// Sorted unique corners of the cells `rRows` (indices into the flat
/// `K`-per-cell `rConn`), minus `v` itself.
template <int K>
std::vector<std::int64_t> decim_corner_ring(const std::vector<std::int64_t>& rRows,
                                            const std::vector<std::int64_t>& rConn,
                                            std::int64_t v) {
    std::vector<std::int64_t> ring;
    ring.reserve(rRows.size() * (K - 1));
    for (std::int64_t f : rRows) {
        const std::int64_t* c = rConn.data() + static_cast<std::size_t>(f) * K;
        for (int k = 0; k < K; ++k)
            if (c[k] != v)
                ring.push_back(c[k]);
    }
    std::sort(ring.begin(), ring.end());
    ring.erase(std::unique(ring.begin(), ring.end()), ring.end());
    return ring;
}

/// Kill the `rShared` cells (removing them from every corner's row), then
/// migrate removed vertex `r`'s remaining alive cells onto survivor `s` and
/// clear `r`'s row. `rConn` is flat `K`-per-cell.
template <int K>
void decim_kill_and_migrate(std::vector<std::int64_t>& rConn,
                            std::vector<std::vector<std::int64_t>>& rRows,
                            std::vector<std::uint8_t>& rAlive,
                            const std::vector<std::int64_t>& rShared, std::int64_t s,
                            std::int64_t r) {
    for (const std::int64_t f : rShared) {
        rAlive[static_cast<std::size_t>(f)] = 0;
        const std::int64_t* c = rConn.data() + static_cast<std::size_t>(f) * K;
        for (int k = 0; k < K; ++k)
            if (c[k] != s && c[k] != r)
                decim_sorted_erase(rRows[static_cast<std::size_t>(c[k])], f);
        decim_sorted_erase(rRows[static_cast<std::size_t>(s)], f);
    }
    for (const std::int64_t f : rRows[static_cast<std::size_t>(r)]) {
        if (!rAlive[static_cast<std::size_t>(f)])
            continue;  // was in shared
        std::int64_t* c = rConn.data() + static_cast<std::size_t>(f) * K;
        for (int k = 0; k < K; ++k)
            if (c[k] == r)
                c[k] = s;
        decim_sorted_insert(rRows[static_cast<std::size_t>(s)], f);
    }
    rRows[static_cast<std::size_t>(r)].clear();
}

/// Marks in `rUsed` every corner of an alive cell of the flat `K`-per-cell
/// `rConn`, numbers the used points ascending into `rRemap` (-1 elsewhere),
/// and returns how many are used.
template <int K>
std::size_t decim_compact_points(const std::vector<std::uint8_t>& rAlive,
                                 const std::vector<std::int64_t>& rConn, std::size_t n,
                                 std::vector<char>& rUsed, std::vector<std::int64_t>& rRemap) {
    rUsed.assign(n, 0);
    for (std::size_t f = 0; f < rAlive.size(); ++f)
        if (rAlive[f])
            for (int k = 0; k < K; ++k)
                rUsed[static_cast<std::size_t>(rConn[f * K + k])] = 1;
    rRemap.assign(n, -1);
    std::size_t num_used = 0;
    for (std::size_t i = 0; i < n; ++i)
        if (rUsed[i])
            rRemap[i] = static_cast<std::int64_t>(num_used++);
    return num_used;
}

}  // namespace detail
}  // namespace meshioplusplus
