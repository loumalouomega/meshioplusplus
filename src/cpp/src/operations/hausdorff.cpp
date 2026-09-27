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

// hausdorff_distance. See operations/hausdorff.hpp for the contract.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/surface_edge_runs.hpp"

namespace meshioplusplus {
namespace {

using detail::Vec3;

constexpr const char* kHdPrefix = "meshio++: hausdorff_distance: ";

bool hd_has_volume(const Mesh& rMesh) {
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            return true;
        if (!cb.IsRagged() && cell_type_dimension(cell_type_from_name(cb.Type())) == 3)
            return true;
    }
    return false;
}

// The triangle soup of one side: the mesh itself, or its linear skin.
detail::TriangleSoup hd_soup(const Mesh& rMesh, const std::string& rRegion, const char* pSide) {
    detail::TriangleSoup soup;
    if (hd_has_volume(rMesh)) {
        if (!rRegion.empty())
            throw std::invalid_argument(std::string(kHdPrefix) + "mesh " + pSide +
                                        " is a volume mesh, whose skin is compared; a region "
                                        "selects surface cells (run extract_surface first)");
        const Mesh skin = detail::surface_extract(rMesh, /*forceFaceMode=*/true,
                                                  /*linearize=*/true, /*recordParentIds=*/false,
                                                  "hausdorff_distance");
        soup = detail::build_triangle_soup(skin, "");
    } else {
        soup = detail::build_triangle_soup(rMesh, rRegion);
    }
    if (soup.NumTriangles() == 0)
        throw std::invalid_argument(std::string(kHdPrefix) + "mesh " + pSide +
                                    " has no surface triangles to measure");
    return soup;
}

// The sample points of a soup: its vertices (the points its triangles use, in
// ascending id), then per triangle the centroids of its `s * s` sub-triangles.
std::vector<Vec3> hd_samples(const detail::TriangleSoup& rSoup, std::int64_t S) {
    std::vector<std::uint8_t> used(rSoup.mPoints.size(), 0);
    for (const auto& tri : rSoup.mVertices)
        for (std::int64_t v : tri)
            used[static_cast<std::size_t>(v)] = 1;
    std::vector<Vec3> out;
    for (std::size_t p = 0; p < used.size(); ++p)
        if (used[p])
            out.push_back(rSoup.mPoints[p]);
    if (S <= 0)
        return out;
    const double s = static_cast<double>(S);
    const std::size_t per = static_cast<std::size_t>(S * S);
    const std::size_t base = out.size();
    out.resize(base + rSoup.NumTriangles() * per);
    parallel_for(rSoup.NumTriangles(), [&](std::size_t t) {
        const Vec3& a = rSoup.mCorners[t * 3 + 0];
        const Vec3 ab = detail::vec3_sub(rSoup.mCorners[t * 3 + 1], a);
        const Vec3 ac = detail::vec3_sub(rSoup.mCorners[t * 3 + 2], a);
        std::size_t k = base + t * per;
        auto at = [&](double U, double V) {
            return detail::vec3_add(
                a, detail::vec3_add(detail::vec3_scale(ab, U / s), detail::vec3_scale(ac, V / s)));
        };
        for (std::int64_t i = 0; i < S; ++i)
            for (std::int64_t j = 0; i + j < S; ++j) {
                out[k++] =
                    at(static_cast<double>(i) + 1.0 / 3.0, static_cast<double>(j) + 1.0 / 3.0);
                if (i + j + 1 < S)
                    out[k++] =
                        at(static_cast<double>(i) + 2.0 / 3.0, static_cast<double>(j) + 2.0 / 3.0);
            }
    });
    return out;
}

struct HdOneSided {
    double mMax = 0.0;
    double mMean = 0.0;
    double mRms = 0.0;
    Vec3 mWorst{0.0, 0.0, 0.0};
};

// Distances from @p rSamples to the surface @p rQuery was built from, reduced
// serially in sample order (a fixed order, so the sums are reproducible).
HdOneSided hd_reduce(const detail::DistanceQuery& rQuery, const std::vector<Vec3>& rSamples) {
    const std::vector<detail::ClosestPointHit> hits =
        detail::query_closest_points(rQuery, rSamples);
    HdOneSided r;
    double sum = 0.0;
    double sum_sq = 0.0;
    for (std::size_t i = 0; i < hits.size(); ++i) {
        const double d = hits[i].mDistance;
        sum += d;
        sum_sq += d * d;
        if (d > r.mMax || i == 0) {
            r.mMax = d;
            r.mWorst = rSamples[i];
        }
    }
    const double n = static_cast<double>(hits.size());
    r.mMean = hits.empty() ? 0.0 : sum / n;
    r.mRms = hits.empty() ? 0.0 : std::sqrt(sum_sq / n);
    return r;
}

detail::DistanceQuery hd_query(const detail::TriangleSoup& rSoup, double CellSize) {
    SurfaceDistanceOptions o;
    o.mSign = SdfSign::Unsigned;
    o.mWatertightCheck = SdfWatertightCheck::Off;
    o.mGridCellSize = CellSize;
    return detail::build_distance_query_from_runs(rSoup, o, detail::surface_edge_runs(rSoup));
}

}  // namespace

HausdorffResult hausdorff_distance(const Mesh& rA, const Mesh& rB,
                                   const HausdorffOptions& rOptions) {
    if (rOptions.mFaceSamples < 0)
        throw std::invalid_argument(std::string(kHdPrefix) + "face samples must be >= 0");
    if (rOptions.mGridCellSize < 0.0)
        throw std::invalid_argument(std::string(kHdPrefix) + "the grid cell size must be >= 0");
    const detail::TriangleSoup sa = hd_soup(rA, rOptions.mRegionA, "A");
    const detail::TriangleSoup sb = hd_soup(rB, rOptions.mRegionB, "B");
    const std::vector<Vec3> pa = hd_samples(sa, rOptions.mFaceSamples);
    const std::vector<Vec3> pb = hd_samples(sb, rOptions.mFaceSamples);

    const HdOneSided ab = hd_reduce(hd_query(sb, rOptions.mGridCellSize), pa);
    const HdOneSided ba = hd_reduce(hd_query(sa, rOptions.mGridCellSize), pb);

    HausdorffResult r;
    r.mAtoB = ab.mMax;
    r.mBtoA = ba.mMax;
    r.mDistance = ab.mMax >= ba.mMax ? ab.mMax : ba.mMax;
    r.mMeanAtoB = ab.mMean;
    r.mRmsAtoB = ab.mRms;
    r.mMeanBtoA = ba.mMean;
    r.mRmsBtoA = ba.mRms;
    r.mNumSamplesA = static_cast<std::int64_t>(pa.size());
    r.mNumSamplesB = static_cast<std::int64_t>(pb.size());
    r.mWorstPointA = {ab.mWorst[0], ab.mWorst[1], ab.mWorst[2]};
    r.mWorstPointB = {ba.mWorst[0], ba.mWorst[1], ba.mWorst[2]};
    return r;
}

}  // namespace meshioplusplus
