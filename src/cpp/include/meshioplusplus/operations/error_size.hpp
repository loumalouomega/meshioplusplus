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
 * @file error_size.hpp
 * @brief A size field from a per-cell error estimate (roadmap 6.10.2).
 *
 * `estimate_error` measures where a field is badly resolved (`error:zz`, one
 * value per cell); `mesh_cad_surface`'s `background` takes a size per point.
 * `error_size` joins them by the adaptive procedure of O. C. Zienkiewicz and
 * J. Z. Zhu, "A simple error estimator and adaptive procedure for practical
 * engineering analysis", Int. J. Numer. Meth. Engng. 24 (1987): with `N`
 * cells evaluated and a target global error `eta`, the error each cell may
 * carry when the error is spread evenly is `e = eta / sqrt(N)`; a cell of
 * size `h_K` (its diameter, the largest distance between two of its nodes)
 * whose indicator is `eta_K` asks for `h_K (e / eta_K)^(1 / p)`, where `p` is
 * the rate at which the indicator falls with `h` (the element order for the
 * recovered-gradient error, 1 for linear elements). The size is clamped to
 * `[mMinSize, mMaxSize]` (0: no bound) and written per point as the least of
 * its cells' (the finer request wins, as every size source of
 * `mesh_cad_surface` refines); a point with no evaluated cell reads NaN,
 * which `background` treats as no size there.
 *
 * Geometry, connectivity, regions, property sets and every existing array
 * pass through. Byte-identical across mesh backends and thread counts.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// Default output array name: what `mesh_cad_surface`'s `background` reads.
inline constexpr const char* kErrorSizeName = "size";

/// Options for `error_size`.
struct ErrorSizeOptions {
    /// The one-component cell_data indicator (`estimate_error` writes
    /// `error:zz`).
    std::string mErrorArray = "error:zz";
    /// The target global error `eta`, in the indicator's norm (positive).
    double mTargetError = 0.0;
    /// The indicator's rate `p` in the cell size (positive; 1 for linear
    /// elements).
    double mOrder = 1.0;
    /// The sizes' bounds; 0 is none.
    double mMinSize = 0.0;
    double mMaxSize = 0.0;
    /// The point_data array written.
    std::string mOutputName = kErrorSizeName;
    /// When false, an output name that already exists is an error.
    bool mOverwrite = false;
};

/// The mesh with the size array, and what the sizes cover.
struct ErrorSizeResult {
    Mesh mMesh;
    std::int64_t mNumCells = 0;    ///< cells evaluated
    std::int64_t mNumSkipped = 0;  ///< cells whose indicator or size is not finite
    std::int64_t mNumClamped = 0;  ///< cell sizes the bounds clamped
    double mMinSize = 0.0;         ///< the least and greatest point size written
    double mMaxSize = 0.0;
};

/**
 * @brief A per-point size field from a per-cell error indicator (see the file
 * comment).
 * @throws std::invalid_argument on a missing or multi-component indicator, a
 *         target error or order that is not finite and positive, a bound that
 *         is negative or not finite, `mMinSize > mMaxSize`, or an output-name
 *         collision when `mOverwrite` is false.
 */
MESHIOPLUSPLUS_API ErrorSizeResult error_size(const Mesh& rMesh, const ErrorSizeOptions& rOptions);

}  // namespace meshioplusplus
