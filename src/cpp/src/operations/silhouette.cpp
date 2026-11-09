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
// The outline of a part seen along a direction (roadmap 7.6.3).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/silhouette.hpp"
#include "../detail/lsq_fit.hpp"
#include "../detail/planar_union.hpp"
#include "../detail/predicates.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/skin.hpp"

namespace meshioplusplus {

namespace {

using SilVec2 = std::array<double, 2>;
using SilTri = std::array<std::int64_t, 3>;

constexpr const char* kSilPrefix = "meshio++: silhouette: ";
constexpr double kSilPi = 3.14159265358979323846;
constexpr double kSilMitreLimit = 4.0;  // a mitre longer than this many clearances is bevelled

[[noreturn]] void sil_refuse(const std::string& rWhy) {
    throw std::invalid_argument(std::string(kSilPrefix) + rWhy);
}

SilVec2 sil_rotate(const SilVec2& rV, double Angle) {
    const double c = std::cos(Angle), s = std::sin(Angle);
    return {rV[0] * c - rV[1] * s, rV[0] * s + rV[1] * c};
}

SilVec2 sil_at(const SilVec2& rP, const SilVec2& rDir, double R) {
    return {rP[0] + R * rDir[0], rP[1] + R * rDir[1]};
}

// Appends a triangle over three new points.
void sil_push(detail::PlanarSoup& rSoup, const SilVec2& rA, const SilVec2& rB, const SilVec2& rC) {
    const auto n = static_cast<std::int64_t>(rSoup.mPoints.size());
    rSoup.mPoints.push_back(rA);
    rSoup.mPoints.push_back(rB);
    rSoup.mPoints.push_back(rC);
    rSoup.mTriangles.push_back(SilTri{n, n + 1, n + 2});
}

// The strip a clearance of |R| sweeps along the region's boundary on one
// side (Sigma +1: outside, the right of each loop edge; -1: inside): a
// rectangle per edge and, at each corner the offset opens, an arc fan whose
// sides are tangent to the true arc, or a mitre (bevelled past the limit).
detail::PlanarSoup sil_offset_strip(const detail::PlanarRegion& rRegion, double R, double Sigma,
                                    double RoundRad, std::int64_t ArcSegments) {
    detail::PlanarSoup strip;
    const auto& pts = rRegion.mPoints;
    for (const auto& r_loop : rRegion.mLoops) {
        const std::size_t n = r_loop.size();
        if (n < 3)
            continue;
        for (std::size_t i = 0; i < n; ++i) {
            const SilVec2& a = pts[static_cast<std::size_t>(r_loop[(i + n - 1) % n])];
            const SilVec2& b = pts[static_cast<std::size_t>(r_loop[i])];
            const SilVec2& c = pts[static_cast<std::size_t>(r_loop[(i + 1) % n])];
            const double l1 = std::hypot(b[0] - a[0], b[1] - a[1]);
            const double l2 = std::hypot(c[0] - b[0], c[1] - b[1]);
            if (!(l1 > 0.0) || !(l2 > 0.0))
                continue;
            const SilVec2 e1{(b[0] - a[0]) / l1, (b[1] - a[1]) / l1};
            const SilVec2 e2{(c[0] - b[0]) / l2, (c[1] - b[1]) / l2};
            const SilVec2 s1{Sigma * e1[1], -Sigma * e1[0]};
            const SilVec2 s2{Sigma * e2[1], -Sigma * e2[0]};
            // The edge b -> c's rectangle.
            sil_push(strip, b, c, sil_at(c, s2, R));
            sil_push(strip, b, sil_at(c, s2, R), sil_at(b, s2, R));
            // The corner at b opens on the offset side when the loop turns
            // towards the other side there.
            const double cross = e1[0] * e2[1] - e1[1] * e2[0];
            const double dot = e1[0] * e2[0] + e1[1] * e2[1];
            if (!(Sigma * cross > 0.0))
                continue;
            const double phi = std::atan2(std::fabs(cross), dot);
            const SilVec2 p0 = sil_at(b, s1, R), p1 = sil_at(b, s2, R);
            if (phi < RoundRad) {
                const double turn = cross > 0.0 ? 1.0 : -1.0;
                const auto k = std::max<std::int64_t>(
                    1, static_cast<std::int64_t>(
                           std::ceil(phi / (2.0 * kSilPi) * static_cast<double>(ArcSegments))));
                const double step = phi / static_cast<double>(k);
                const double rr = R / std::cos(step / 2.0);
                SilVec2 prev = p0;
                for (std::int64_t j = 0; j < k; ++j) {
                    const SilVec2 q =
                        sil_at(b, sil_rotate(s1, turn * (static_cast<double>(j) + 0.5) * step), rr);
                    sil_push(strip, b, prev, q);
                    prev = q;
                }
                sil_push(strip, b, prev, p1);
            } else if (1.0 / std::cos(phi / 2.0) <= kSilMitreLimit) {
                const double d = 1.0 + s1[0] * s2[0] + s1[1] * s2[1];
                const SilVec2 m{b[0] + R * (s1[0] + s2[0]) / d, b[1] + R * (s1[1] + s2[1]) / d};
                sil_push(strip, b, p0, m);
                sil_push(strip, b, m, p1);
            } else {
                sil_push(strip, b, p0, p1);
            }
        }
    }
    return strip;
}

detail::PlanarSoup sil_soup_of(const detail::PlanarRegion& rRegion) {
    return detail::PlanarSoup{rRegion.mPoints, rRegion.mTriangles};
}

void sil_append(detail::PlanarSoup& rInto, const detail::PlanarSoup& rFrom) {
    const auto base = static_cast<std::int64_t>(rInto.mPoints.size());
    rInto.mPoints.insert(rInto.mPoints.end(), rFrom.mPoints.begin(), rFrom.mPoints.end());
    for (const SilTri& r_t : rFrom.mTriangles)
        rInto.mTriangles.push_back(SilTri{r_t[0] + base, r_t[1] + base, r_t[2] + base});
}

NDArray sil_int64(const std::vector<std::int64_t>& rValues) {
    NDArray out = NDArray::Uninit(DType::Int64, {rValues.size()});
    std::copy(rValues.begin(), rValues.end(), out.As<std::int64_t>());
    return out;
}

NDArray sil_doubles(const std::vector<double>& rValues) {
    NDArray out = NDArray::Uninit(DType::Float64, {rValues.size()});
    std::copy(rValues.begin(), rValues.end(), out.As<double>());
    return out;
}

}  // namespace

SilhouetteResult silhouette(const Mesh& rMesh, const SilhouetteOptions& rOptions) {
    const double nn = detail::vec3_norm(rOptions.mNormal);
    if (!(nn > 0.0) || !std::isfinite(nn))
        sil_refuse("the normal must be a non-zero finite vector");
    if (!std::isfinite(rOptions.mClearance))
        sil_refuse("the clearance must be finite");
    if (!std::isfinite(rOptions.mRoundAngleDeg) || rOptions.mRoundAngleDeg < 0.0)
        sil_refuse("round_angle must be a finite angle of at least 0 degrees");
    if (rOptions.mArcSegments < 4)
        sil_refuse("arc_segments must be at least 4");
    if (!std::isfinite(rOptions.mMergeTolerance) || rOptions.mMergeTolerance < 0.0)
        sil_refuse("merge_tolerance must be finite and not negative");

    bool volume = false;
    for (const auto cb : rMesh.CellRange())
        volume = volume || cb.IsPolyhedron() ||
                 cell_type_dimension(cell_type_from_name(std::string(cb.Type()))) == 3;
    if (volume && !rOptions.mRegion.empty())
        sil_refuse(
            "a region selects surface cells; a volume mesh is outlined through its skin "
            "(extract_skin first, then the region)");
    const detail::TriangleSoup soup = volume ? detail::build_triangle_soup(extract_skin(rMesh), "")
                                             : detail::build_triangle_soup(rMesh, rOptions.mRegion);
    if (soup.NumTriangles() == 0)
        sil_refuse("the mesh has no surface cells to project (triangles, quads or polygons)");

    const detail::Vec3 n = detail::vec3_scale(rOptions.mNormal, 1.0 / nn);
    detail::Vec3 e1{}, e2{};
    detail::lsq_frame(n.data(), e1.data(), e2.data());

    // Project the points the triangles use.
    detail::PlanarSoup add;
    add.mPoints.resize(soup.mPoints.size());
    std::vector<char> used(soup.mPoints.size(), 0);
    for (const auto& r_t : soup.mVertices)
        for (const std::int64_t v : r_t)
            used[static_cast<std::size_t>(v)] = 1;
    SilhouetteResult result;
    bool first = true;
    double lo_u = 0, lo_v = 0, hi_u = 0, hi_v = 0;
    for (std::size_t i = 0; i < soup.mPoints.size(); ++i) {
        const detail::Vec3& r_p = soup.mPoints[i];
        add.mPoints[i] = {detail::vec3_dot(r_p, e1), detail::vec3_dot(r_p, e2)};
        if (!used[i])
            continue;
        const double h = detail::vec3_dot(r_p, n);
        if (first) {
            result.mMinHeight = result.mMaxHeight = h;
            lo_u = hi_u = add.mPoints[i][0];
            lo_v = hi_v = add.mPoints[i][1];
            first = false;
        }
        result.mMinHeight = std::min(result.mMinHeight, h);
        result.mMaxHeight = std::max(result.mMaxHeight, h);
        lo_u = std::min(lo_u, add.mPoints[i][0]);
        hi_u = std::max(hi_u, add.mPoints[i][0]);
        lo_v = std::min(lo_v, add.mPoints[i][1]);
        hi_v = std::max(hi_v, add.mPoints[i][1]);
    }
    if (!std::isfinite(result.mMinHeight) || !std::isfinite(result.mMaxHeight))
        sil_refuse("the mesh has non-finite coordinates");
    add.mTriangles.assign(soup.mVertices.begin(), soup.mVertices.end());
    for (const SilTri& r_t : add.mTriangles)
        if (detail::orient2d(add.mPoints[static_cast<std::size_t>(r_t[0])].data(),
                             add.mPoints[static_cast<std::size_t>(r_t[1])].data(),
                             add.mPoints[static_cast<std::size_t>(r_t[2])].data()) != 0.0)
            ++result.mNumTriangles;
    const double diagonal = std::hypot(hi_u - lo_u, hi_v - lo_v);
    const double tol = rOptions.mMergeTolerance * diagonal;

    detail::PlanarRegion region = detail::planar_region(add, detail::PlanarSoup{}, tol);
    if (!region.mOk)
        sil_refuse(region.mWhy);
    if (rOptions.mClearance != 0.0 && !region.mLoops.empty()) {
        const double r = std::fabs(rOptions.mClearance);
        const double sigma = rOptions.mClearance > 0.0 ? 1.0 : -1.0;
        const detail::PlanarSoup strip = sil_offset_strip(
            region, r, sigma, rOptions.mRoundAngleDeg * kSilPi / 180.0, rOptions.mArcSegments);
        detail::PlanarSoup base = sil_soup_of(region);
        if (sigma > 0.0) {
            sil_append(base, strip);
            region = detail::planar_region(base, detail::PlanarSoup{}, tol);
        } else {
            region = detail::planar_region(base, strip, tol);
        }
        if (!region.mOk)
            sil_refuse("offsetting by the clearance: " + region.mWhy);
        if (region.mLoops.empty())
            result.mWarnings.push_back("the clearance removes the whole outline");
    }
    if (region.mTriangles.empty() && rOptions.mClearance == 0.0)
        result.mWarnings.push_back("the surface projects to no area along the normal");

    for (const SilTri& r_t : region.mTriangles)
        result.mArea +=
            0.5 * detail::orient2d(region.mPoints[static_cast<std::size_t>(r_t[0])].data(),
                                   region.mPoints[static_cast<std::size_t>(r_t[1])].data(),
                                   region.mPoints[static_cast<std::size_t>(r_t[2])].data());
    for (const auto& r_loop : region.mLoops) {
        ++result.mNumLoops;
        if (detail::planar_loop_twice_area(region.mPoints, r_loop) < 0.0)
            ++result.mNumHoles;
        else
            ++result.mNumOuter;
        for (std::size_t i = 0; i < r_loop.size(); ++i) {
            const SilVec2& a = region.mPoints[static_cast<std::size_t>(r_loop[i])];
            const SilVec2& b =
                region.mPoints[static_cast<std::size_t>(r_loop[(i + 1) % r_loop.size()])];
            result.mPerimeter += std::hypot(b[0] - a[0], b[1] - a[1]);
        }
    }

    // The output: the outline at the lowest height, or the prism.
    auto lift = [&](const SilVec2& rP, double H) {
        detail::Vec3 q{};
        for (std::size_t k = 0; k < 3; ++k)
            q[k] = rP[0] * e1[k] + rP[1] * e2[k] + H * n[k];
        return q;
    };
    Mesh& out = result.mMesh;
    std::vector<double> xyz;
    if (!rOptions.mExtrude) {
        std::unordered_map<std::int64_t, std::int64_t> id_of;
        std::vector<std::int64_t> conn, chain, closed, hole;
        std::int64_t k = 0;
        for (const auto& r_loop : region.mLoops) {
            const std::int64_t is_hole =
                detail::planar_loop_twice_area(region.mPoints, r_loop) < 0.0 ? 1 : 0;
            auto id = [&](std::int64_t V) {
                const auto it = id_of.find(V);
                if (it != id_of.end())
                    return it->second;
                const auto nid = static_cast<std::int64_t>(id_of.size());
                id_of.emplace(V, nid);
                const detail::Vec3 q =
                    lift(region.mPoints[static_cast<std::size_t>(V)], result.mMinHeight);
                xyz.insert(xyz.end(), q.begin(), q.end());
                return nid;
            };
            for (std::size_t i = 0; i < r_loop.size(); ++i) {
                const std::int64_t a = id(r_loop[i]);
                const std::int64_t b = id(r_loop[(i + 1) % r_loop.size()]);
                conn.push_back(a);
                conn.push_back(b);
                chain.push_back(k);
                closed.push_back(1);
                hole.push_back(is_hole);
            }
            ++k;
        }
        NDArray p = NDArray::Uninit(DType::Float64, {xyz.size() / 3, std::size_t{3}});
        std::copy(xyz.begin(), xyz.end(), p.As<double>());
        out.AssignPoints(std::move(p));
        NDArray c = NDArray::Uninit(DType::Int64, {conn.size() / 2, std::size_t{2}});
        std::copy(conn.begin(), conn.end(), c.As<std::int64_t>());
        out.AddCellBlock("line", std::move(c));
        std::vector<NDArray> b_chain, b_closed, b_offset, b_hole;
        b_chain.push_back(sil_int64(chain));
        b_closed.push_back(sil_int64(closed));
        b_offset.push_back(sil_doubles(std::vector<double>(chain.size(), result.mMinHeight)));
        b_hole.push_back(sil_int64(hole));
        out.AddCellData(kSectionChainName, std::move(b_chain));
        out.AddCellData(kSectionClosedName, std::move(b_closed));
        out.AddCellData(kSectionOffsetName, std::move(b_offset));
        out.AddCellData(kSilhouetteHoleName, std::move(b_hole));
    } else {
        if (!(result.mMaxHeight > result.mMinHeight))
            result.mWarnings.push_back(
                "the part is flat along the normal: the prism has no height");
        const auto m = static_cast<std::int64_t>(region.mPoints.size());
        for (const double h : {result.mMinHeight, result.mMaxHeight})
            for (const SilVec2& r_p : region.mPoints) {
                const detail::Vec3 q = lift(r_p, h);
                xyz.insert(xyz.end(), q.begin(), q.end());
            }
        std::vector<std::int64_t> conn;
        for (const SilTri& r_t : region.mTriangles) {
            conn.insert(conn.end(), {r_t[0], r_t[2], r_t[1]});
            conn.insert(conn.end(), {r_t[0] + m, r_t[1] + m, r_t[2] + m});
        }
        for (const auto& r_loop : region.mLoops)
            for (std::size_t i = 0; i < r_loop.size(); ++i) {
                const std::int64_t a = r_loop[i], b = r_loop[(i + 1) % r_loop.size()];
                conn.insert(conn.end(), {a, b, b + m});
                conn.insert(conn.end(), {a, b + m, a + m});
            }
        NDArray p = NDArray::Uninit(DType::Float64, {xyz.size() / 3, std::size_t{3}});
        std::copy(xyz.begin(), xyz.end(), p.As<double>());
        out.AssignPoints(std::move(p));
        NDArray c = NDArray::Uninit(DType::Int64, {conn.size() / 3, std::size_t{3}});
        std::copy(conn.begin(), conn.end(), c.As<std::int64_t>());
        out.AddCellBlock("triangle", std::move(c));
    }
    out.AddFieldData(kSilhouetteNormalName, sil_doubles({n[0], n[1], n[2]}));
    out.AddFieldData(kSilhouetteExtentName, sil_doubles({result.mMinHeight, result.mMaxHeight}));
    if (rMesh.HasFieldData(detail::kLengthUnitName))
        out.AddFieldData(detail::kLengthUnitName, rMesh.FieldData(detail::kLengthUnitName));
    return result;
}

}  // namespace meshioplusplus
