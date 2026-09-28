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

// match_periodic_nodes. See operations/periodic.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/periodic.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/node_adjacency.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kPrPrefix = "meshio++: match_periodic_nodes: ";

// The nodes a region names, ascending and unique.
std::vector<std::int64_t> pr_region_points(const Mesh& rMesh, const meshioplusplus::Region& rR) {
    const std::size_t n = rMesh.NumPoints();
    std::vector<std::int64_t> out;
    const std::size_t count = rR.NumEntries();
    if (rR.mKind == RegionKind::Point) {
        for (std::size_t i = 0; i < count; ++i)
            out.push_back(detail::read_int(rR.mEntries, i));
    } else if (rR.mKind == RegionKind::Cell) {
        const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
        std::vector<std::int64_t> ids;
        for (std::size_t i = 0; i < count; ++i) {
            const auto [block, row] =
                detail::global_to_block_row(bases, detail::read_int(rR.mEntries, i));
            if (block == static_cast<std::size_t>(-1))
                continue;
            detail::cell_node_ids(rMesh.Cells(block), static_cast<std::size_t>(row), n, ids);
            out.insert(out.end(), ids.begin(), ids.end());
        }
    } else {
        CellType type;
        std::vector<std::int64_t> ids;
        for (std::size_t i = 0; i < count; ++i)
            if (detail::facet_nodes(rMesh, detail::read_int(rR.mEntries, i * 2),
                                    detail::read_int(rR.mEntries, i * 2 + 1), type, ids))
                out.insert(out.end(), ids.begin(), ids.end());
    }
    out.erase(
        std::remove_if(out.begin(), out.end(),
                       [n](std::int64_t v) { return v < 0 || static_cast<std::size_t>(v) >= n; }),
        out.end());
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

NDArray pr_ids(const std::vector<std::int64_t>& rIds) {
    NDArray a = NDArray::Uninit(DType::Int64, {rIds.size()});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    return a;
}

std::string pr_first_ids(const std::vector<std::int64_t>& rIds) {
    std::string s;
    for (std::size_t i = 0; i < rIds.size() && i < 8; ++i)
        s += (i ? ", " : "") + std::to_string(rIds[i]);
    if (rIds.size() > 8)
        s += ", ...";
    return s;
}

}  // namespace

PeriodicPairs match_periodic_nodes(const Mesh& rMesh, const RegionSelector& rSlave,
                                   const RegionSelector& rMaster, const PeriodicOptions& rOptions) {
    if (!(rOptions.mAtol > 0.0) || !std::isfinite(rOptions.mAtol))
        throw std::invalid_argument(std::string(kPrPrefix) + "the tolerance must be positive");
    const std::vector<std::int64_t> slave =
        pr_region_points(rMesh, rMesh.Region(find_region(rMesh, rSlave)));
    const std::vector<std::int64_t> master =
        pr_region_points(rMesh, rMesh.Region(find_region(rMesh, rMaster)));

    const std::size_t dim = rMesh.PointDim();
    const NDArray& pts = rMesh.Points();
    auto point = [&](std::int64_t Id) {
        detail::Vec3 p{0.0, 0.0, 0.0};
        for (std::size_t d = 0; d < dim && d < 3; ++d)
            p[d] = detail::read_double(pts, static_cast<std::size_t>(Id) * dim + d);
        return p;
    };

    // A grid over the master nodes, inserted in ascending id so every bucket is
    // ascending. The cell is at least the tolerance, so a match within it is
    // always in the 3x3x3 neighbourhood; it is floored relative to the extent
    // so a tiny tolerance cannot overflow the integer keys.
    double extent = 0.0;
    for (std::int64_t m : master) {
        const detail::Vec3 p = point(m);
        for (double c : p)
            extent = std::max(extent, std::fabs(c));
    }
    const double cell = std::max(rOptions.mAtol, extent * 1e-12);
    detail::SpatialGrid grid(cell);
    std::vector<detail::Vec3> master_xyz(master.size());
    for (std::size_t k = 0; k < master.size(); ++k) {
        master_xyz[k] = point(master[k]);
        grid.Insert(grid.KeyOf(master_xyz[k].data()), static_cast<std::int64_t>(k));
    }

    const std::array<double, 16>& m = rOptions.mTransform.mMatrix;
    const double atol2 = rOptions.mAtol * rOptions.mAtol;
    std::vector<std::int64_t> match(slave.size(), -1);  // index into `master`
    std::vector<double> residual(slave.size(), 0.0);
    parallel_for(slave.size(), [&](std::size_t i) {
        const detail::Vec3 p = point(slave[i]);
        detail::Vec3 q;
        for (std::size_t r = 0; r < 3; ++r)
            q[r] = m[r * 4 + 0] * p[0] + m[r * 4 + 1] * p[1] + m[r * 4 + 2] * p[2] + m[r * 4 + 3];
        double best = std::numeric_limits<double>::infinity();
        std::int64_t hit = -1;
        grid.ForEachIn27(grid.KeyOf(q.data()), [&](const std::vector<std::int64_t>& rBucket) {
            for (std::int64_t k : rBucket) {
                const detail::Vec3& c = master_xyz[static_cast<std::size_t>(k)];
                const double d2 = (c[0] - q[0]) * (c[0] - q[0]) + (c[1] - q[1]) * (c[1] - q[1]) +
                                  (c[2] - q[2]) * (c[2] - q[2]);
                if (d2 <= atol2 &&
                    (d2 < best || (d2 == best && master[static_cast<std::size_t>(k)] <
                                                     master[static_cast<std::size_t>(hit)]))) {
                    best = d2;
                    hit = k;
                }
            }
            return true;
        });
        match[i] = hit;
        residual[i] = hit >= 0 ? std::sqrt(best) : 0.0;
    });

    PeriodicPairs out;
    std::vector<std::int64_t> s_ids;
    std::vector<std::int64_t> m_ids;
    std::vector<std::int64_t> unmatched;
    std::vector<std::int64_t> claimed(master.size(), -1);
    for (std::size_t i = 0; i < slave.size(); ++i) {
        if (match[i] < 0) {
            unmatched.push_back(slave[i]);
            continue;
        }
        const std::size_t k = static_cast<std::size_t>(match[i]);
        if (master[k] == slave[i]) {
            ++out.mNumFixed;
            continue;
        }
        if (claimed[k] >= 0)
            throw std::invalid_argument(std::string(kPrPrefix) + "slave nodes " +
                                        std::to_string(claimed[k]) + " and " +
                                        std::to_string(slave[i]) + " both map onto master node " +
                                        std::to_string(master[k]) +
                                        " (the tolerance is too loose, or the regions overlap)");
        claimed[k] = slave[i];
        s_ids.push_back(slave[i]);
        m_ids.push_back(master[k]);
        out.mMaxResidual = std::max(out.mMaxResidual, residual[i]);
    }
    char atol[32];
    detail::snprintf_c(atol, sizeof(atol), "%g", rOptions.mAtol);
    if (!unmatched.empty() && rOptions.mRequireComplete)
        throw std::invalid_argument(std::string(kPrPrefix) + std::to_string(unmatched.size()) +
                                    " slave node(s) have no master node within " +
                                    std::string(atol) + " (" + pr_first_ids(unmatched) +
                                    "); check the transform, or pass require_complete=false");
    out.mSlave = pr_ids(s_ids);
    out.mMaster = pr_ids(m_ids);
    out.mUnmatched = pr_ids(unmatched);
    return out;
}

}  // namespace meshioplusplus
