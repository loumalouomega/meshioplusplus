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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

// Project includes
#include "mesh_cad_size.hpp"
#include "tri_bvh.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

// The cell of a coordinate on a grid of edge `Cell` from `Lo`.
std::int64_t mcsz_cell(double X, double Lo, double Cell) {
    return static_cast<std::int64_t>(std::floor((X - Lo) / Cell));
}

}  // namespace

MeshCadSizeField::MeshCadSizeField(double H, double Floor, double Slope, double Reach,
                                   const std::vector<Vec3>& rX, const std::vector<double>& rSize)
    : mH(H), mFloor(Floor), mSlope(Slope), mReach(Reach), mSmallest(H) {
    for (std::size_t i = 0; i < rX.size(); ++i) {
        const double s = std::max(Floor, rSize[i]);
        if (!(s < H))
            continue;
        mX.push_back(rX[i]);
        mSize.push_back(s);
        mSmallest = std::min(mSmallest, s);
    }
    if (mX.empty())
        return;
    // Cells of edge H: a sample farther than r + (H - h_i) / s from a point
    // cannot lower it below H, so only a few rings around it are visited.
    mCell = H;
    Vec3 hi = mX.front();
    mLo = mX.front();
    for (const Vec3& r_x : mX)
        for (int d = 0; d < 3; ++d) {
            mLo[d] = std::min(mLo[d], r_x[d]);
            hi[d] = std::max(hi[d], r_x[d]);
        }
    for (int d = 0; d < 3; ++d)
        mDims[d] = mcsz_cell(hi[d], mLo[d], mCell) + 1;
    for (std::size_t i = 0; i < mX.size(); ++i)
        mCells[Key(mcsz_cell(mX[i][0], mLo[0], mCell), mcsz_cell(mX[i][1], mLo[1], mCell),
                   mcsz_cell(mX[i][2], mLo[2], mCell))]
            .push_back(static_cast<std::int32_t>(i));
}

std::int64_t MeshCadSizeField::Key(std::int64_t I, std::int64_t J, std::int64_t K) const {
    return I + mDims[0] * (J + mDims[1] * K);
}

double MeshCadSizeField::At(const Vec3& rX) const {
    double best = mH;
    if (mpBackground)
        best = std::max(mFloor, std::min(best, mpBackground->At(rX)));
    if (mX.empty())
        return best;
    std::array<std::int64_t, 3> c{};
    for (int d = 0; d < 3; ++d)
        c[d] = mcsz_cell(rX[d], mLo[d], mCell);
    // Ring r holds the cells at Chebyshev distance r, all farther than
    // (r - 1) cells; stop once no sample there could beat the best.
    const std::int64_t reach = mDims[0] + mDims[1] + mDims[2];
    for (std::int64_t r = 0; r <= reach; ++r) {
        if (r > 0 && mSmallest + mSlope * std::max(0.0, double(r - 1) * mCell - mReach) >= best)
            break;
        for (std::int64_t k = c[2] - r; k <= c[2] + r; ++k) {
            if (k < 0 || k >= mDims[2])
                continue;
            for (std::int64_t j = c[1] - r; j <= c[1] + r; ++j) {
                if (j < 0 || j >= mDims[1])
                    continue;
                for (std::int64_t i = c[0] - r; i <= c[0] + r; ++i) {
                    if (i < 0 || i >= mDims[0])
                        continue;
                    if (std::max({std::abs(i - c[0]), std::abs(j - c[1]), std::abs(k - c[2])}) != r)
                        continue;
                    const auto it = mCells.find(Key(i, j, k));
                    if (it == mCells.end())
                        continue;
                    for (const std::int32_t s : it->second) {
                        const double d = std::max(
                            0.0, vec3_norm(vec3_sub(rX, mX[static_cast<std::size_t>(s)])) - mReach);
                        best = std::min(best, mSize[static_cast<std::size_t>(s)] + mSlope * d);
                    }
                }
            }
        }
    }
    return std::max(mFloor, best);
}

std::vector<double> mesh_cad_sample_sizes(const std::vector<MeshCadSizeSample>& rSamples,
                                          const Mesh& rFirst,
                                          const std::vector<std::vector<std::int64_t>>& rAdjacent,
                                          double H, double Floor, double Deviation,
                                          double CellsAcrossGaps) {
    const std::size_t n = rSamples.size();
    std::vector<double> out(n, H);
    // Curvature (roadmap 6.3.2): a chord of length l on a circle of radius
    // 1 / k sits d = R - sqrt(R^2 - l^2 / 4) from it, so l = 2 sqrt(2 d R - d^2).
    if (Deviation > 0.0)
        for (std::size_t i = 0; i < n; ++i) {
            const double k = rSamples[i].mCurvature;
            if (!(k > 0.0) || !std::isfinite(k))
                continue;
            const double q = 2.0 * Deviation / k - Deviation * Deviation;
            if (q > 0.0)
                out[i] = std::min(out[i], 2.0 * std::sqrt(q) / kMeshCadSizeSlack);
        }
    // Local feature size (roadmap 6.3.7): along the normal, inward to the far
    // wall of a thin part, outward to the near wall across a gap, skipping
    // the point's own face and its neighbours.
    if (CellsAcrossGaps > 0.0 && rFirst.NumCellBlocks() > 0) {
        const TriangleSoup soup = build_triangle_soup(rFirst, "");
        const TriBvh tree = build_tri_bvh(soup);
        std::vector<std::int64_t> cell_face;
        for (std::size_t b = 0; b < rFirst.NumCellBlocks(); ++b) {
            const std::size_t nc = rFirst.Cells(b).NumCells();
            const bool tagged =
                rFirst.HasCellData("cad:face") && rFirst.CellData("cad:face", b).Size() == nc;
            for (std::size_t c = 0; c < nc; ++c)
                cell_face.push_back(tagged ? read_int(rFirst.CellData("cad:face", b), c) : -1);
        }
        const auto excluded = [&](std::int64_t Face, std::int64_t Other) {
            if (Other < 0 || Other == Face)
                return true;
            const auto& r_adj = rAdjacent[static_cast<std::size_t>(Face)];
            return std::binary_search(r_adj.begin(), r_adj.end(), Other);
        };
        const double tmin = 1e-9 * H;
        const double inf = std::numeric_limits<double>::infinity();
        std::vector<double> gap(n, inf);
        parallel_for(
            n,
            [&](std::size_t i) {
                const MeshCadSizeSample& r_s = rSamples[i];
                if (r_s.mFace < 0 || vec3_norm_sq(r_s.mNormal) == 0.0)
                    return;
                for (const double side : {-1.0, 1.0}) {
                    const Vec3 dir = vec3_scale(r_s.mNormal, side);
                    double from = tmin;
                    for (int hop = 0; hop < 64; ++hop) {
                        const TriRayHit hit = tri_bvh_ray(tree, r_s.mX, dir, from, inf);
                        if (hit.mTriangle < 0)
                            break;
                        const std::int64_t face = cell_face[static_cast<std::size_t>(
                            soup.mSourceCell[static_cast<std::size_t>(hit.mTriangle)])];
                        if (excluded(r_s.mFace, face)) {
                            from = hit.mT + std::max(tmin, 1e-12 * hit.mT);
                            continue;
                        }
                        // Inward, the far wall faces along the ray; outward,
                        // the wall across the gap faces back.
                        const double facing = vec3_dot(hit.mNormal, dir);
                        if ((side < 0.0 && facing > 0.0) || (side > 0.0 && facing < 0.0))
                            gap[i] = std::min(gap[i], hit.mT);
                        break;
                    }
                }
            },
            256);
        for (std::size_t i = 0; i < n; ++i)
            if (std::isfinite(gap[i]))
                out[i] = std::min(out[i], gap[i] / (CellsAcrossGaps * kMeshCadSizeSlack));
    }
    for (double& r_h : out)
        r_h = std::max(Floor, std::min(H, r_h));
    return out;
}

struct MeshCadBackgroundSize::Impl {
    TriangleSoup mSoup;
    DistanceQuery mQuery;
    std::vector<double> mSize;  ///< per point of the surface
};

MeshCadBackgroundSize::MeshCadBackgroundSize(const Mesh& rBackground, const std::string& rArray)
    : mpImpl(std::make_unique<Impl>()) {
    bool volume = false;
    for (std::size_t b = 0; b < rBackground.NumCellBlocks(); ++b)
        volume =
            volume || cell_type_dimension(cell_type_from_name(rBackground.Cells(b).Type())) == 3;
    // A skin keeps the point_data of the points it keeps.
    const Mesh skin = volume ? extract_surface(rBackground) : Mesh{};
    const Mesh& r_surface = volume ? skin : rBackground;
    const NDArray& r_size = r_surface.PointData(rArray);
    mpImpl->mSize.resize(r_surface.NumPoints());
    for (std::size_t i = 0; i < mpImpl->mSize.size(); ++i)
        mpImpl->mSize[i] = read_double(r_size, i);
    mpImpl->mSoup = build_triangle_soup(r_surface, "");
    if (mpImpl->mSoup.NumTriangles() == 0)
        return;
    SurfaceDistanceOptions opt;
    opt.mSign = SdfSign::Unsigned;
    opt.mWatertightCheck = SdfWatertightCheck::Off;
    mpImpl->mQuery = build_distance_query(mpImpl->mSoup, opt);
}

MeshCadBackgroundSize::~MeshCadBackgroundSize() = default;

double MeshCadBackgroundSize::At(const Vec3& rX) const {
    const double none = std::numeric_limits<double>::infinity();
    const TriangleSoup& r_soup = mpImpl->mSoup;
    if (r_soup.NumTriangles() == 0)
        return none;
    const SurfaceProjection hit = query_surface_projections(mpImpl->mQuery, {rX}).front();
    if (!hit.mFound || hit.mTriangle < 0)
        return none;
    const auto& r_v = r_soup.mVertices[static_cast<std::size_t>(hit.mTriangle)];
    const Vec3& a = r_soup.mPoints[static_cast<std::size_t>(r_v[0])];
    const Vec3& b = r_soup.mPoints[static_cast<std::size_t>(r_v[1])];
    const Vec3& c = r_soup.mPoints[static_cast<std::size_t>(r_v[2])];
    // The hit's barycentric coordinates in its triangle.
    const Vec3 e1 = vec3_sub(b, a), e2 = vec3_sub(c, a), ep = vec3_sub(hit.mPoint, a);
    const double d11 = vec3_dot(e1, e1), d12 = vec3_dot(e1, e2), d22 = vec3_dot(e2, e2);
    const double dp1 = vec3_dot(ep, e1), dp2 = vec3_dot(ep, e2);
    const double den = d11 * d22 - d12 * d12;
    double wb = 1.0 / 3.0, wc = 1.0 / 3.0;
    if (den > 0.0) {
        wb = std::clamp((d22 * dp1 - d12 * dp2) / den, 0.0, 1.0);
        wc = std::clamp((d11 * dp2 - d12 * dp1) / den, 0.0, 1.0 - wb);
    }
    const auto& r_s = mpImpl->mSize;
    const double s = (1.0 - wb - wc) * r_s[static_cast<std::size_t>(r_v[0])] +
                     wb * r_s[static_cast<std::size_t>(r_v[1])] +
                     wc * r_s[static_cast<std::size_t>(r_v[2])];
    return std::isfinite(s) && s > 0.0 ? s : none;
}

}  // namespace detail
}  // namespace meshioplusplus
