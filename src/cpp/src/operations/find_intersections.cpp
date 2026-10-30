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
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/find_intersections.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/ndarray.hpp"
// Project includes (private, not installed)
#include "../detail/tri_pairs.hpp"
#include "../detail/tri_tri.hpp"

namespace meshioplusplus {

namespace {

struct FiHit {
    std::int64_t mJ = 0;
    detail::TriTriSegment mSegment;
};

// The input, bit-identical, with the count per cell.
Mesh fi_flagged(const Mesh& rMesh, const std::vector<std::int64_t>& rCount, bool Attach) {
    Mesh out = detail::clone_mesh(
        rMesh, [&](DataLocation Location, const std::string& rName, std::string&) {
            return !(Attach && Location == DataLocation::Cell && rName == kIntersectionsCountName);
        });
    if (Attach) {
        const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
        std::vector<NDArray> blocks;
        for (std::size_t b = 0; b + 1 < bases.size(); ++b) {
            const auto n = static_cast<std::size_t>(bases[b + 1] - bases[b]);
            NDArray a = NDArray::Uninit(DType::Int64, {n});
            std::copy(rCount.begin() + bases[b], rCount.begin() + bases[b + 1],
                      a.As<std::int64_t>());
            blocks.push_back(std::move(a));
        }
        out.AddCellData(kIntersectionsCountName, std::move(blocks));
    }
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        out.AddRegion(rMesh.Region(i));
    return out;
}

FindIntersectionsResult fi_run(const Mesh& rA, const Mesh* pB,
                               const FindIntersectionsOptions& rOptions) {
    if (rOptions.mMaxPairs < 0)
        throw std::invalid_argument("meshio++: find_intersections: max_pairs must be >= 0");
    const bool self = pB == nullptr;
    const detail::TriangleSoup soup_a = detail::build_triangle_soup(rA, rOptions.mRegion);
    const detail::TriangleSoup soup_b =
        self ? detail::TriangleSoup{} : detail::build_triangle_soup(*pB, rOptions.mRegionB);
    const detail::TriangleSoup& r_other = self ? soup_a : soup_b;

    FindIntersectionsResult result;
    if (soup_a.NumTriangles() == 0)
        result.mWarnings.emplace_back("the first mesh has no 2-D cells to test");
    if (!self && soup_b.NumTriangles() == 0)
        result.mWarnings.emplace_back("the second mesh has no 2-D cells to test");

    // Per triangle of A, the triangles it meets, ascending.
    std::vector<std::vector<FiHit>> hits(soup_a.NumTriangles());
    detail::soup_box_pairs(
        soup_a, self ? nullptr : &soup_b, [&](std::size_t I, const std::vector<std::int64_t>& rJs) {
            const detail::Vec3* a = &soup_a.mCorners[3 * I];
            for (std::int64_t j : rJs) {
                const auto uj = static_cast<std::size_t>(j);
                if (self && soup_a.mSourceCell[uj] == soup_a.mSourceCell[I])
                    continue;
                const detail::Vec3* b = &r_other.mCorners[3 * uj];
                const bool meet =
                    self ? detail::tri_tri_intersect_unshared(a[0].data(), a[1].data(), a[2].data(),
                                                              b[0].data(), b[1].data(), b[2].data())
                         : detail::tri_tri_intersect(a[0].data(), a[1].data(), a[2].data(),
                                                     b[0].data(), b[1].data(), b[2].data());
                if (!meet)
                    continue;
                FiHit hit;
                hit.mJ = j;
                hit.mSegment = detail::tri_tri_segment(a[0].data(), a[1].data(), a[2].data(),
                                                       b[0].data(), b[1].data(), b[2].data());
                hits[I].push_back(hit);
            }
        });

    // Reduce in triangle order.
    std::vector<std::array<std::int64_t, 2>> pairs;
    std::vector<double> seg_points;
    std::vector<std::int64_t> seg_a, seg_b;
    for (std::size_t i = 0; i < hits.size(); ++i)
        for (const FiHit& r_hit : hits[i]) {
            const std::int64_t ca = soup_a.mSourceCell[i];
            const std::int64_t cb = r_other.mSourceCell[static_cast<std::size_t>(r_hit.mJ)];
            pairs.push_back(self ? std::array<std::int64_t, 2>{{std::min(ca, cb), std::max(ca, cb)}}
                                 : std::array<std::int64_t, 2>{{ca, cb}});
            ++result.mNumTrianglePairs;
            if (r_hit.mSegment.mKind == detail::TriTriContact::Coplanar)
                ++result.mNumCoplanar;
            if (rOptions.mSegments && r_hit.mSegment.mKind == detail::TriTriContact::Segment) {
                seg_points.insert(seg_points.end(), r_hit.mSegment.mP, r_hit.mSegment.mP + 3);
                seg_points.insert(seg_points.end(), r_hit.mSegment.mQ, r_hit.mSegment.mQ + 3);
                seg_a.push_back(ca);
                seg_b.push_back(cb);
            }
        }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    result.mNumPairs = static_cast<std::int64_t>(pairs.size());
    result.mClean = pairs.empty();

    // Partners per cell.
    const std::vector<std::int64_t> bases_a = detail::block_bases(rA);
    std::vector<std::int64_t> count_a(static_cast<std::size_t>(detail::total_cells(bases_a)), 0);
    std::vector<std::int64_t> count_b;
    if (!self)
        count_b.assign(static_cast<std::size_t>(detail::total_cells(detail::block_bases(*pB))), 0);
    for (const std::array<std::int64_t, 2>& r_pair : pairs) {
        ++count_a[static_cast<std::size_t>(r_pair[0])];
        ++(self ? count_a : count_b)[static_cast<std::size_t>(r_pair[1])];
    }
    for (std::int64_t c : count_a)
        result.mNumCellsA += c > 0;
    for (std::int64_t c : count_b)
        result.mNumCellsB += c > 0;

    if (rOptions.mMaxPairs > 0 && result.mNumPairs > rOptions.mMaxPairs) {
        pairs.resize(static_cast<std::size_t>(rOptions.mMaxPairs));
        result.mTruncated = true;
        result.mWarnings.emplace_back(
            "the pair list is capped at max_pairs = " + std::to_string(rOptions.mMaxPairs) +
            " of " + std::to_string(result.mNumPairs) + " pairs");
    }
    result.mPairs = std::move(pairs);

    result.mMesh = fi_flagged(rA, count_a, rOptions.mAttachFlags);
    if (!self)
        result.mMeshB = fi_flagged(*pB, count_b, rOptions.mAttachFlags);

    if (rOptions.mSegments) {
        const std::size_t ns = seg_a.size();
        result.mNumSegments = static_cast<std::int64_t>(ns);
        NDArray points = NDArray::Uninit(DType::Float64, {2 * ns, std::size_t{3}});
        std::copy(seg_points.begin(), seg_points.end(), points.As<double>());
        result.mSegments.AssignPoints(std::move(points));
        NDArray line = NDArray::Uninit(DType::Int64, {ns, std::size_t{2}});
        std::int64_t* p_line = line.As<std::int64_t>();
        for (std::size_t s = 0; s < 2 * ns; ++s)
            p_line[s] = static_cast<std::int64_t>(s);
        result.mSegments.AddCellBlock("line", std::move(line));
        for (const auto& [name, values] : {std::pair<const char*, const std::vector<std::int64_t>*>{
                                               kIntersectionsCellAName, &seg_a},
                                           std::pair<const char*, const std::vector<std::int64_t>*>{
                                               kIntersectionsCellBName, &seg_b}}) {
            NDArray a = NDArray::Uninit(DType::Int64, {ns});
            std::copy(values->begin(), values->end(), a.As<std::int64_t>());
            std::vector<NDArray> blocks;
            blocks.push_back(std::move(a));
            result.mSegments.AddCellData(name, std::move(blocks));
        }
    }
    return result;
}

}  // namespace

FindIntersectionsResult find_intersections(const Mesh& rMesh,
                                           const FindIntersectionsOptions& rOptions) {
    return fi_run(rMesh, nullptr, rOptions);
}

FindIntersectionsResult find_intersections(const Mesh& rA, const Mesh& rB,
                                           const FindIntersectionsOptions& rOptions) {
    return fi_run(rA, &rB, rOptions);
}

}  // namespace meshioplusplus
