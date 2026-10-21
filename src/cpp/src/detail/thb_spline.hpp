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
 * @file detail/thb_spline.hpp
 * @brief Truncated hierarchical B-splines (THB) on the unit square, and the
 * adaptive freeform fit with them behind `local_refine="thb"` (v17.12.0,
 * experimental), free of OpenCASCADE.
 *
 * The hierarchy (Giannelli, Jüttler and Speleers, "THB-splines: the
 * truncated basis for hierarchical splines", CAGD 2012; the construction is
 * the paper's, reimplemented): level `l` is the clamped uniform cubic
 * B-spline space with `s0 2^l` spans per direction, `Omega^0` the square and
 * `Omega^{l+1}` a union of level-`l` cells inside `Omega^l` (the cells marked
 * refined). A level-`l` B-spline is active when its support lies in
 * `Omega^l` but not in `Omega^{l+1}`; it is truncated by dropping, level by
 * level, its components along the finer B-splines whose support lies in the
 * next domain. The truncated functions are a basis, they sum to one, and
 * every one is a combination of the finest level's B-splines -- so a THB
 * surface is **flattened** exactly into one tensor-product B-spline over the
 * finest knots, which is what a face is written as.
 *
 * The fit (`thb_fit_patch`) takes `bspline_fit_patch`'s inputs and options:
 * least squares to the vertices and facet centres, thin-plate-like smoothing
 * of the finest coefficients (strong where no sample is), parameter
 * correction by projection, and, where samples stay beyond the target, the
 * level-`l` cells holding them (and their neighbours) refined -- cells, not
 * whole knot lines, so the degrees of freedom grow where the shape needs
 * them. `mMaxTotalPoles` caps the hierarchical degrees of freedom, and the
 * finest level stays within `mMaxPoles` per direction.
 *
 * A **core-private** header: nothing installed names it. Serial and
 * deterministic.
 */

// System includes
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// Project includes (private, not installed)
#include "bspline_fit.hpp"
#include "bspline_lsq.hpp"

namespace meshioplusplus {
namespace detail {

/// A THB hierarchy on the unit square and its truncated basis.
struct ThbBasis {
    int mDegree = 3;
    std::size_t mSpans0 = 1;  ///< level-0 spans per direction
    /// Per level `l`, its `(s0 2^l)^2` cells, row by row in u: 1 where the
    /// cell is refined into level `l + 1`. The last level refines none.
    std::vector<std::vector<std::uint8_t>> mRefined;
    /// Per truncated function: its level, its B-spline's indices at that
    /// level, and its coefficients on the finest level's B-splines (pole
    /// index `i + n j`), sparse.
    struct Function {
        int mLevel = 0;
        std::size_t mI = 0, mJ = 0;
        std::vector<std::pair<std::size_t, double>> mFine;
    };
    std::vector<Function> mFunctions;

    int NumLevels() const { return static_cast<int>(mRefined.size()); }
    std::size_t Spans(int Level) const { return mSpans0 << Level; }
    /// The finest level's B-splines per direction.
    std::size_t FineSize() const {
        return Spans(NumLevels() - 1) + static_cast<std::size_t>(mDegree);
    }
};

/// A hierarchy of `Levels` levels with no cell refined yet.
ThbBasis thb_hierarchy(int Degree, std::size_t Spans0, int Levels);

/// Fill `mFunctions` from the refined cells (active and truncated).
void thb_build(ThbBasis& rB);

/// The flattened surface: the finest level's tensor-product B-spline on the
/// unit square whose poles are `sum_k rCoef[k] (fine coefficients of k)`.
BsplineSurface thb_flatten(const ThbBasis& rB, const std::vector<PatchVec>& rCoef);

/// What a THB fit adds to `BsplineFitResult`.
struct ThbFitInfo {
    std::size_t mNumDofs = 0;  ///< the hierarchical degrees of freedom
    int mNumLevels = 1;        ///< levels the hierarchy reached
};

/**
 * @brief Fit a patch with a THB surface (see the file comment), flattened:
 * `bspline_fit_patch`'s contract, the surface on the unit square.
 */
BsplineFitResult thb_fit_patch(const std::vector<PatchVec>& rPts,
                               const std::vector<PatchTri>& rTris, const std::vector<double>& rU,
                               const std::vector<double>& rV, const BsplineFitOptions& rOptions,
                               ThbFitInfo* pInfo = nullptr);

}  // namespace detail
}  // namespace meshioplusplus
