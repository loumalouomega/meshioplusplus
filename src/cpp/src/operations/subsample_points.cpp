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

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/operations/subsample_points.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"

// Project includes (private, not installed)
#include "../detail/farthest_points.hpp"
#include "../detail/subset_gather.hpp"

namespace meshioplusplus {

namespace {

constexpr std::int64_t kSpGridOversample = 2;
constexpr int kSpGridMaxRefinements = 12;

struct SpGridEntry {
    std::int64_t mCell = 0;
    double mDist2 = 0.0;
    std::int64_t mLocal = 0;

    bool operator<(const SpGridEntry& r) const {
        if (mCell != r.mCell)
            return mCell < r.mCell;
        if (mDist2 != r.mDist2)
            return mDist2 < r.mDist2;
        return mLocal < r.mLocal;
    }
};

std::vector<std::int64_t> sp_run_fps(const std::vector<double>& rXyz, std::size_t N,
                                     std::size_t Count, std::size_t First,
                                     const std::string& rPrefix) {
    const detail::FarthestPoints fp = detail::farthest_point_order(rXyz.data(), N, Count, First);
    if (fp.mSelected < Count)
        throw std::invalid_argument(rPrefix + "the cloud has only " + std::to_string(fp.mSelected) +
                                    " distinct positions among its " + std::to_string(N) +
                                    " candidates, fewer than the requested " +
                                    std::to_string(Count) +
                                    "; lower the count or weld the duplicates first");
    return fp.mOrder;
}

std::vector<std::int64_t> sp_grid_select(const std::vector<double>& rLocal, std::size_t M,
                                         std::size_t Count, std::size_t FirstLocal,
                                         std::int64_t MaxCells, const std::string& rPrefix) {
    double lo[3] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[3] = {-std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity()};
    for (std::size_t i = 0; i < M; ++i) {
        for (int d = 0; d < 3; ++d) {
            const double v = rLocal[3 * i + static_cast<std::size_t>(d)];
            if (v < lo[d])
                lo[d] = v;
            if (v > hi[d])
                hi[d] = v;
        }
    }
    double extent[3] = {hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
    bool live[3] = {extent[0] > 0.0, extent[1] > 0.0, extent[2] > 0.0};
    int nlive = (live[0] ? 1 : 0) + (live[1] ? 1 : 0) + (live[2] ? 1 : 0);
    if (nlive == 0)
        throw std::invalid_argument(
            rPrefix +
            "every candidate sits at the same position; there is nothing to select among");

    double measure = 1.0;
    double max_ext = 0.0;
    for (int d = 0; d < 3; ++d) {
        if (live[d]) {
            measure *= extent[d];
            if (extent[d] > max_ext)
                max_ext = extent[d];
        }
    }
    const double target = static_cast<double>(kSpGridOversample) * static_cast<double>(Count);
    double h = std::pow(measure / target, 1.0 / static_cast<double>(nlive));
    h = std::max(h, max_ext * 1e-12);

    std::vector<std::int64_t> reps;
    std::vector<std::int64_t> cell_of(M, 0);
    std::array<std::int64_t, 3> dims{1, 1, 1};
    for (int ref = 0; ref <= kSpGridMaxRefinements; ++ref) {
        std::array<std::int64_t, 3> cur_dims{1, 1, 1};
        double prod = 1.0;
        for (int d = 0; d < 3; ++d) {
            cur_dims[static_cast<std::size_t>(d)] =
                live[d]
                    ? std::max<std::int64_t>(1, static_cast<std::int64_t>(std::ceil(extent[d] / h)))
                    : 1;
            prod *= static_cast<double>(cur_dims[static_cast<std::size_t>(d)]);
        }
        if (prod > static_cast<double>(MaxCells))
            break;
        dims = cur_dims;
        double spacing[3] = {live[0] ? h : 1.0, live[1] ? h : 1.0, live[2] ? h : 1.0};
        std::vector<SpGridEntry> entries(M);
        for (std::size_t i = 0; i < M; ++i) {
            std::int64_t ijk[3];
            double d2 = 0.0;
            for (int d = 0; d < 3; ++d) {
                const double v = rLocal[3 * i + static_cast<std::size_t>(d)];
                std::int64_t k = static_cast<std::int64_t>(std::floor((v - lo[d]) / spacing[d]));
                k = std::min(std::max<std::int64_t>(k, 0), dims[static_cast<std::size_t>(d)] - 1);
                ijk[d] = k;
                const double centre = lo[d] + (static_cast<double>(k) + 0.5) * spacing[d];
                const double diff = v - centre;
                d2 += diff * diff;
            }
            const std::int64_t cid = ijk[0] + dims[0] * (ijk[1] + dims[1] * ijk[2]);
            cell_of[i] = cid;
            entries[i] = SpGridEntry{cid, d2, static_cast<std::int64_t>(i)};
        }
        std::sort(entries.begin(), entries.end());
        reps.clear();
        for (std::size_t i = 0; i < M; ++i) {
            if (i == 0 || entries[i].mCell != entries[i - 1].mCell)
                reps.push_back(entries[i].mLocal);
        }
        if (static_cast<std::int64_t>(reps.size()) >=
                kSpGridOversample * static_cast<std::int64_t>(Count) ||
            reps.size() == M)
            break;
        h *= 0.5;
    }

    if (reps.size() < Count) {
        throw std::invalid_argument(
            rPrefix + "the cloud occupies only " + std::to_string(reps.size()) +
            " distinct lattice cells at the finest lattice tried (" + std::to_string(dims[0]) +
            "x" + std::to_string(dims[1]) + "x" + std::to_string(dims[2]) +
            "), fewer than the requested " + std::to_string(Count) +
            "; lower the count, raise max_cells, or use method='farthest'");
    }

    std::unordered_map<std::int64_t, std::size_t> rep_of_cell;
    rep_of_cell.reserve(reps.size() * 2);
    for (std::size_t i = 0; i < reps.size(); ++i)
        rep_of_cell[cell_of[static_cast<std::size_t>(reps[i])]] = i;
    const std::size_t first_rep = rep_of_cell.at(cell_of[FirstLocal]);

    std::vector<double> rep_xyz(3 * reps.size());
    for (std::size_t i = 0; i < reps.size(); ++i) {
        const auto r = static_cast<std::size_t>(reps[i]);
        rep_xyz[3 * i] = rLocal[3 * r];
        rep_xyz[3 * i + 1] = rLocal[3 * r + 1];
        rep_xyz[3 * i + 2] = rLocal[3 * r + 2];
    }
    const std::vector<std::int64_t> rep_order =
        sp_run_fps(rep_xyz, reps.size(), Count, first_rep, rPrefix);
    std::vector<std::int64_t> order(Count);
    for (std::size_t i = 0; i < Count; ++i)
        order[i] = reps[static_cast<std::size_t>(rep_order[i])];
    return order;
}

}  // namespace

SubsamplePointsResult subsample_points(const Mesh& rMesh, const SubsamplePointsOptions& rOptions) {
    const std::string prefix = "meshio++: subsample_points: ";
    if (rOptions.mCount < 1)
        throw std::invalid_argument(prefix + "count must be at least 1, got " +
                                    std::to_string(rOptions.mCount));
    if (rOptions.mMethod != SubsampleMethod::Farthest &&
        rOptions.mMethod != SubsampleMethod::Grid && rOptions.mMethod != SubsampleMethod::Random)
        throw std::invalid_argument(prefix + "method must be 'farthest', 'grid' or 'random'");
    if (rOptions.mMaxCells < 1)
        throw std::invalid_argument(prefix + "max_cells must be at least 1");

    const std::size_t n = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();
    std::vector<double> xyz(3 * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto pt = detail::read_point(rMesh.Points(), pdim, static_cast<std::int64_t>(i));
        if (!std::isfinite(pt[0] + pt[1] + pt[2]))
            throw std::invalid_argument(prefix + "point " + std::to_string(i) + " is not finite");
        xyz[3 * i] = pt[0];
        xyz[3 * i + 1] = pt[1];
        xyz[3 * i + 2] = pt[2];
    }

    std::vector<std::int64_t> candidates;
    if (!rOptions.mHasBounds) {
        candidates.resize(n);
        std::iota(candidates.begin(), candidates.end(), std::int64_t{0});
    } else {
        for (int d = 0; d < 3; ++d) {
            const double lo = rOptions.mBounds[static_cast<std::size_t>(d)];
            const double hi = rOptions.mBounds[static_cast<std::size_t>(d + 3)];
            if (!std::isfinite(lo) || !std::isfinite(hi) || hi < lo)
                throw std::invalid_argument(prefix + "bounds have hi < lo or non-finite value");
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (xyz[3 * i] >= rOptions.mBounds[0] && xyz[3 * i] <= rOptions.mBounds[3] &&
                xyz[3 * i + 1] >= rOptions.mBounds[1] && xyz[3 * i + 1] <= rOptions.mBounds[4] &&
                xyz[3 * i + 2] >= rOptions.mBounds[2] && xyz[3 * i + 2] <= rOptions.mBounds[5])
                candidates.push_back(static_cast<std::int64_t>(i));
        }
    }

    const std::size_t m = candidates.size();
    const auto count = static_cast<std::size_t>(rOptions.mCount);
    if (count > m)
        throw std::invalid_argument(prefix + "count is " + std::to_string(count) + " but only " +
                                    std::to_string(m) + " points are available");

    std::vector<double> local(3 * m);
    for (std::size_t i = 0; i < m; ++i) {
        const auto c = static_cast<std::size_t>(candidates[i]);
        local[3 * i] = xyz[3 * c];
        local[3 * i + 1] = xyz[3 * c + 1];
        local[3 * i + 2] = xyz[3 * c + 2];
    }

    std::mt19937_64 rng(static_cast<std::uint64_t>(rOptions.mSeed));
    std::size_t first_local = 0;
    if (rOptions.mMethod != SubsampleMethod::Random) {
        if (rOptions.mStart < 0) {
            first_local = static_cast<std::size_t>(rng() % m);
        } else {
            const auto it = std::lower_bound(candidates.begin(), candidates.end(), rOptions.mStart);
            if (it == candidates.end() || *it != rOptions.mStart)
                throw std::invalid_argument(prefix + "start=" + std::to_string(rOptions.mStart) +
                                            " is not a candidate point index");
            first_local = static_cast<std::size_t>(it - candidates.begin());
        }
    }

    std::vector<std::int64_t> local_order;
    if (rOptions.mMethod == SubsampleMethod::Farthest) {
        local_order = sp_run_fps(local, m, count, first_local, prefix);
    } else if (rOptions.mMethod == SubsampleMethod::Grid) {
        local_order = sp_grid_select(local, m, count, first_local, rOptions.mMaxCells, prefix);
    } else {
        std::vector<std::int64_t> pool(m);
        std::iota(pool.begin(), pool.end(), std::int64_t{0});
        local_order.resize(count);
        for (std::size_t i = 0; i < count; ++i) {
            const std::size_t pick = i + static_cast<std::size_t>(rng() % (m - i));
            std::swap(pool[i], pool[pick]);
            local_order[i] = pool[i];
        }
    }

    std::vector<std::int64_t> selected(count);
    for (std::size_t i = 0; i < count; ++i)
        selected[i] = candidates[static_cast<std::size_t>(local_order[i])];

    SubsamplePointsResult result;
    Mesh out;
    detail::subset_gather_points(rMesh, selected, out);
    if (rOptions.mRecordIds) {
        NDArray orig = NDArray::Uninit(DType::Int64, {count});
        std::copy(selected.begin(), selected.end(), orig.As<std::int64_t>());
        out.AddPointData(kBudgetOriginalPointIdName, std::move(orig));
    }
    NDArray conn = NDArray::Uninit(DType::Int64, {count, std::size_t{1}});
    for (std::size_t i = 0; i < count; ++i)
        conn.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    out.AddCellBlock("vertex", std::move(conn));

    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    // Warnings for dropped non-vertex cells, cell_data, and Cell/Side regions.
    std::size_t dropped_cells = 0;
    std::size_t dropped_blocks = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.Type() != "vertex") {
            dropped_cells += cb.NumCells();
            ++dropped_blocks;
        }
    }
    if (dropped_cells > 0)
        result.mWarnings.push_back("dropped " + std::to_string(dropped_cells) + " cells in " +
                                   std::to_string(dropped_blocks) +
                                   " blocks (a point cloud has no cells)");
    if (!rMesh.CellDataNames().empty())
        result.mWarnings.push_back("dropped cell_data with the cells");

    std::vector<std::int64_t> inverse(n, -1);
    for (std::size_t j = 0; j < count; ++j)
        inverse[static_cast<std::size_t>(selected[j])] = static_cast<std::int64_t>(j);
    std::vector<std::string> dropped_regions;
    for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
        const Region& region = rMesh.Region(r);
        if (region.mKind != RegionKind::Point) {
            dropped_regions.push_back(region.mName);
            continue;
        }
        std::vector<std::int64_t> mapped;
        const std::int64_t* entries = region.Entries();
        for (std::size_t e = 0; e < region.NumEntries(); ++e) {
            const std::int64_t p = entries[e];
            if (p >= 0 && static_cast<std::size_t>(p) < n &&
                inverse[static_cast<std::size_t>(p)] >= 0)
                mapped.push_back(inverse[static_cast<std::size_t>(p)]);
        }
        NDArray arr = NDArray::Uninit(DType::Int64, {mapped.size()});
        std::copy(mapped.begin(), mapped.end(), arr.As<std::int64_t>());
        out.AddRegion(
            Region(region.mName, RegionKind::Point, region.mDim, region.mTag, std::move(arr)));
    }
    if (!dropped_regions.empty()) {
        std::string names;
        for (std::size_t i = 0; i < dropped_regions.size(); ++i) {
            if (i > 0)
                names += ", ";
            names += dropped_regions[i];
        }
        result.mWarnings.push_back("dropped " + std::to_string(dropped_regions.size()) +
                                   " cell/side region(s) (" + names + ") with the cells");
    }

    result.mMesh = std::move(out);
    result.mPointMap = std::move(selected);
    return result;
}

}  // namespace meshioplusplus
