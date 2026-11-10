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
// Predicate-free ODT remeshing of a tetrahedral mesh: alternate ODT vertex
// relocation (reused verbatim from `smooth`'s SmoothMethod::Odt) with
// quality-improving 2-3 / 3-2 topological flips. See
// operations/optimize_volume.hpp for the contract, and doc/optimize_volume.md
// for the algorithm and attribution (Freitag & Ollivier-Gooch's improvement
// rule; no in-sphere/Delaunay predicate).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/optimize_volume.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/fem_kernels.hpp"
#include "../detail/tet_flips.hpp"
#include "../detail/tet_scope.hpp"
#include "../detail/typed_view.hpp"

// Project includes (private, not installed)
#include "smooth_odt.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

constexpr const char* kOvPrefix = "meshio++: optimize_volume: ";

// --- tet-only scope check (shared with adapt_volume) --------------------------

void optvol_check_blocks(const Mesh& rMesh) {
    detail::require_tet_only_blocks(rMesh, kOvPrefix);
}

// --- the shared moves (detail/tet_flips.hpp, moved there verbatim) -----------

using Tet = detail::TetFlipTet;
using detail::tet_flip_boundary_mask;
using detail::tet_flip_guarded_odt;
using detail::tet_flip_min_dihedral_tet;
using detail::tet_flip_orient_positive;
using detail::tet_flip_pass_23;
using detail::tet_flip_pass_32;
using detail::tet_flip_pass_edge_removal;
using detail::tet_flip_pass_multi_face;
using detail::tet_flip_quality_tet;
using detail::tet_flip_smooth_worst;
using detail::tet_flip_worst_q;

// --- ODT relocation half (smooth's SmoothMethod::Odt) -----------------------
// The options of the one-iteration ODT pass each sweep runs: smooth's own
// boundary/feature/frozen pinning and inversion guard.
SmoothOptions optvol_odt_options(const OptimizeVolumeOptions& rOptions) {
    SmoothOptions so;
    so.mMethod = SmoothMethod::Odt;
    so.mIterations = 1;
    so.mFixBoundary = rOptions.mPreserveBoundary;
    so.mGuardInversion = true;
    so.mFrozen = rOptions.mFrozen;
    return so;
}

// The nodes the relocation holds still, computed once per call: the caller's
// frozen mask and (with mPreserveBoundary) the boundary. Neither flip changes
// the boundary -- a 2-3 flip replaces two tets across an interior face and a
// 3-2 flip three tets around a closed edge ring by tets with the same outer
// faces -- so the per-sweep boundary rebuild smooth() used to do is hoisted
// here. The mesh is the single-tetra-block Float64 one smooth() was given.
std::vector<std::uint8_t> optvol_pin_mask(const std::vector<double>& rXyz,
                                          const std::vector<Tet>& rTets,
                                          const SmoothOptions& rOdt) {
    const std::size_t n = rXyz.size() / 3;
    Mesh tmp;
    NDArray pts = NDArray::Uninit(DType::Float64, {n, 3});
    std::memcpy(pts.Data(), rXyz.data(), rXyz.size() * sizeof(double));
    tmp.AssignPoints(std::move(pts));
    NDArray conn = NDArray::Uninit(DType::Int64, {rTets.size(), 4});
    static_assert(sizeof(Tet) == 4 * sizeof(std::int64_t), "Tet rows are copied as a flat block");
    std::memcpy(conn.Data(), rTets.data(), rTets.size() * sizeof(Tet));
    tmp.AddCellBlock(cell_type_name(CellType::Tetra), std::move(conn));
    return detail::smooth_odt_pin_mask(tmp, rOdt);
}

}  // namespace

OptimizeVolumeObjective optimize_volume_objective_from_name(const std::string& rName) {
    if (rName == "worst")
        return OptimizeVolumeObjective::Worst;
    if (rName == "energy")
        return OptimizeVolumeObjective::Energy;
    throw std::invalid_argument(std::string(kOvPrefix) + "unknown objective '" + rName +
                                "' (expected 'worst' or 'energy')");
}

OptimizeVolumeResult optimize_volume(const Mesh& rMesh, const OptimizeVolumeOptions& rOptions) {
    optvol_check_blocks(rMesh);

    OptimizeVolumeResult result;
    const std::size_t n = static_cast<std::size_t>(rMesh.NumPoints());
    const std::size_t dim = rMesh.PointDim();

    // --- working buffers ------------------------------------------------------
    std::vector<double> xyz(n * 3, 0.0);
    {
        const NDArray& points = rMesh.Points();
        const detail::DoubleView points_v(points);
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t d = 0; d < dim && d < 3; ++d)
                xyz[i * 3 + d] = points_v[i * dim + d];
        });
    }
    const std::vector<double> xyz0 = xyz;  // for the moved-count at the end

    std::vector<Tet> tets;
    {
        const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
        tets.resize(static_cast<std::size_t>(detail::total_cells(bases)));
        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            const NDArray& conn = cb.Conn();
            const detail::Int64View conn_v(conn);
            const std::size_t nc = cb.NumCells();
            const std::size_t base = static_cast<std::size_t>(bases[bi]);
            for (std::size_t c = 0; c < nc; ++c)
                for (std::size_t k = 0; k < 4; ++k)
                    tets[base + c][k] = conn_v[c * 4 + k];
            ++bi;
        }
    }

    // Scale-relative epsilons (bbox diagonal), so validity is scale-invariant.
    double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    if (n > 0) {
        for (int d = 0; d < 3; ++d)
            lo[d] = hi[d] = xyz[d];
        for (std::size_t i = 1; i < n; ++i)
            for (int d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], xyz[i * 3 + d]);
                hi[d] = std::max(hi[d], xyz[i * 3 + d]);
            }
    }
    const double L =
        std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                  (hi[2] - lo[2]) * (hi[2] - lo[2]));
    const double eps_len = (L > 0.0 ? L : 1.0) * 1e-14;
    const double eps_vol6 = (L > 0.0 ? L * L * L : 1.0) * 1e-14;

    auto min_quality = [&](const std::vector<Tet>& ts) {
        double m = 1.0;
        for (const Tet& t : ts)
            m = std::min(m, tet_flip_quality_tet(xyz, t, eps_len));
        return ts.empty() ? 0.0 : m;
    };
    result.mMinQualityBefore = min_quality(tets);

    auto min_dihedral = [&](const std::vector<Tet>& ts) {
        double m = 180.0;
        for (const Tet& t : ts)
            m = std::min(m, tet_flip_min_dihedral_tet(xyz, t));
        return ts.empty() ? 0.0 : m;
    };
    result.mMinDihedralBefore = min_dihedral(tets);

    // --- the optimisation loop ------------------------------------------------
    const SmoothOptions odt = optvol_odt_options(rOptions);
    std::vector<std::uint8_t> pinned;
    if (rOptions.mRelocate && rOptions.mMaxIterations > 0)
        pinned = optvol_pin_mask(xyz, tets, odt);
    if (rOptions.mObjective == OptimizeVolumeObjective::Energy) {
        const auto quality = [&](const Tet& t) { return tet_flip_quality_tet(xyz, t, eps_len); };
        for (int sweep = 0; sweep < rOptions.mMaxIterations; ++sweep) {
            std::int64_t moved_before = 0;
            if (rOptions.mRelocate) {
                const std::vector<double> before = xyz;
                detail::smooth_odt_pass(xyz, tets, pinned, odt);
                for (std::size_t i = 0; i < n; ++i) {
                    const double dx = xyz[i * 3] - before[i * 3];
                    const double dy = xyz[i * 3 + 1] - before[i * 3 + 1];
                    const double dz = xyz[i * 3 + 2] - before[i * 3 + 2];
                    if (std::sqrt(dx * dx + dy * dy + dz * dz) > 1e-12 * (L > 0.0 ? L : 1.0))
                        ++moved_before;
                }
            }
            std::int64_t flips_this_sweep = 0;
            if (rOptions.mFlip) {
                const std::int64_t n23 =
                    tet_flip_pass_23(xyz, tets, rOptions.mMinImprovement, eps_vol6, quality);
                const std::int64_t n32 =
                    tet_flip_pass_32(xyz, tets, rOptions.mMinImprovement, eps_len, eps_vol6);
                result.mNum23Flips += n23;
                result.mNum32Flips += n32;
                flips_this_sweep = n23 + n32;
            }
            result.mNumFlips += flips_this_sweep;
            if (moved_before == 0 && flips_this_sweep == 0)
                break;  // fixed point
        }
    } else {
        // The worst-cell schedule: guarded relocation, then the topological
        // moves, then worst-cell smoothing; every step keeps the least angle.
        if (rOptions.mMaxRingSize < 3)
            throw std::invalid_argument(std::string(kOvPrefix) + "max_ring_size must be >= 3");
        for (Tet& t : tets)
            t = tet_flip_orient_positive(xyz, t);
        const auto quality = [&](const Tet& t) { return tet_flip_worst_q(xyz, t, eps_vol6); };
        std::vector<std::uint8_t> fixed;
        if (rOptions.mRelocate && rOptions.mMaxIterations > 0) {
            fixed = tet_flip_boundary_mask(n, tets);
            for (std::size_t i = 0; i < n; ++i)
                fixed[i] = static_cast<std::uint8_t>(fixed[i] || pinned[i]);
        }
        constexpr double kWorstThresholdDeg = 40.0;
        for (int sweep = 0; sweep < rOptions.mMaxIterations; ++sweep) {
            std::int64_t moved = 0;
            if (rOptions.mRelocate)
                moved = tet_flip_guarded_odt(xyz, tets, pinned, odt, eps_vol6);
            std::int64_t flips = 0;
            if (rOptions.mFlip) {
                const std::int64_t n23 = tet_flip_pass_23(xyz, tets, rOptions.mMinImprovement,
                                                          eps_vol6, quality, kWorstThresholdDeg);
                std::int64_t n32 = 0;
                const std::int64_t nedge =
                    tet_flip_pass_edge_removal(xyz, tets, rOptions.mMinImprovement, eps_vol6,
                                               rOptions.mMaxRingSize, kWorstThresholdDeg, n32);
                const std::int64_t nmulti =
                    tet_flip_pass_multi_face(xyz, tets, rOptions.mMinImprovement, eps_vol6,
                                             rOptions.mMaxRingSize, kWorstThresholdDeg);
                result.mNum23Flips += n23;
                result.mNum32Flips += n32;
                result.mNumEdgeRemovals += nedge - n32;
                result.mNumMultiFaceRemovals += nmulti;
                flips = n23 + nedge + nmulti;
            }
            result.mNumFlips += flips;
            std::int64_t smoothed = 0;
            if (rOptions.mRelocate)
                smoothed = tet_flip_smooth_worst(xyz, tets, fixed, eps_vol6, kWorstThresholdDeg);
            result.mNumWorstSmoothed += smoothed;
            if (moved == 0 && flips == 0 && smoothed == 0)
                break;  // fixed point
        }
    }

    result.mMinQualityAfter = min_quality(tets);
    result.mMinDihedralAfter = min_dihedral(tets);
    result.mNumTets = static_cast<std::int64_t>(tets.size());
    for (std::size_t i = 0; i < n; ++i) {
        const double dx = xyz[i * 3] - xyz0[i * 3];
        const double dy = xyz[i * 3 + 1] - xyz0[i * 3 + 1];
        const double dz = xyz[i * 3 + 2] - xyz0[i * 3 + 2];
        if (std::sqrt(dx * dx + dy * dy + dz * dz) > 1e-9 * (L > 0.0 ? L : 1.0))
            ++result.mNumVerticesMoved;
    }

    // --- emit: single tetra block, points invariant --------------------------
    Mesh& out = result.mMesh;
    {
        const NDArray& points = rMesh.Points();
        const detail::DoubleView points_v(points);
        NDArray np = NDArray::Uninit(points.Dtype(), {n, dim});
        detail::dispatch_dtype(points.Dtype(), [&]<class T>() {
            T* dst = np.As<T>();
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t d = 0; d < dim; ++d)
                    dst[i * dim + d] = d < 3 ? static_cast<T>(xyz[i * 3 + d])
                                             : static_cast<T>(points_v[i * dim + d]);
        });
        out.AssignPoints(std::move(np));
    }
    {
        NDArray block = NDArray::Uninit(DType::Int64, {tets.size(), 4});
        std::int64_t* dst = block.As<std::int64_t>();
        for (std::size_t t = 0; t < tets.size(); ++t) {
            const Tet ot = tet_flip_orient_positive(xyz, tets[t]);
            for (int k = 0; k < 4; ++k)
                dst[t * 4 + k] = ot[k];
        }
        out.AddCellBlock(cell_type_name(CellType::Tetra), std::move(block));
    }

    // point_data + field_data + Point regions carry (point set is invariant);
    // cell_data + Cell/Side regions are dropped -- flips have no cell map.
    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::data_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    std::size_t dropped_cell_regions = 0;
    for (std::size_t i = 0; i < static_cast<std::size_t>(rMesh.NumRegions()); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind == RegionKind::Point)
            out.AddRegion(
                Region(r.mName, r.mKind, r.mDim, r.mTag, detail::data_owned_copy(r.mEntries)));
        else
            ++dropped_cell_regions;
    }
    if (!rMesh.CellDataNames().empty())
        log::warn("{}dropped {} cell_data array(s): a flip has no cell correspondence", kOvPrefix,
                  rMesh.CellDataNames().size());
    if (dropped_cell_regions > 0)
        log::warn("{}dropped {} Cell/Side region(s): a flip has no cell correspondence", kOvPrefix,
                  dropped_cell_regions);

    return result;
}

}  // namespace meshioplusplus
