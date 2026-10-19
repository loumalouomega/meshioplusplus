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
// The planar skin shared by voxelize (pixels) and octree (quadtrees). See
// detail/planar_skin.hpp.

// System includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

// Project includes
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/planar_skin.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {
namespace detail {

PlanarPlane planar_plane_from_name(const std::string& rName) {
    if (rName == "auto" || rName.empty())
        return PlanarPlane::Auto;
    if (rName == "xy")
        return PlanarPlane::XY;
    if (rName == "xz")
        return PlanarPlane::XZ;
    if (rName == "yz")
        return PlanarPlane::YZ;
    throw std::invalid_argument("meshio++: unknown plane '" + rName +
                                "' (expected one of: auto, xy, xz, yz)");
}

PlanarFrame planar_frame(const Mesh& rMesh, PlanarPlane Plane) {
    PlanarFrame frame;
    if (Plane == PlanarPlane::Auto) {
        std::array<double, 3> lo{}, hi{};
        if (point_bbox(rMesh, lo, hi)) {
            for (int k : {2, 1, 0}) {
                const std::size_t kk = static_cast<std::size_t>(k);
                if (hi[kk] == lo[kk]) {
                    Plane = k == 2 ? PlanarPlane::XY : (k == 1 ? PlanarPlane::XZ : PlanarPlane::YZ);
                    break;
                }
            }
        }
        if (Plane == PlanarPlane::Auto)
            Plane = PlanarPlane::XY;
    }
    if (Plane == PlanarPlane::XZ)
        frame.mAxis = {{0, 2, 1}};
    else if (Plane == PlanarPlane::YZ)
        frame.mAxis = {{1, 2, 0}};
    return frame;
}

PlanarSkin collect_planar_skin(const Mesh& rMesh, const char* pPrefix, const PlanarFrame& rFrame) {
    PlanarSkin skin;
    const std::size_t dim = rMesh.PointDim();
    const NDArray& points = rMesh.Points();
    const std::size_t npts = rMesh.NumPoints();
    std::map<std::array<std::int64_t, 2>, std::pair<std::int64_t, std::array<std::int64_t, 2>>>
        boundary;
    std::vector<std::array<std::int64_t, 2>> order;
    for (const auto cb : rMesh.CellRange()) {
        const std::string type(cb.Type());
        const bool line = type == "line";
        const bool surf = type == "triangle" || type == "quad";
        if (!line && !surf)
            continue;
        const std::size_t n = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::array<std::int64_t, 4> v{};
            for (std::size_t i = 0; i < n; ++i) {
                v[i] = read_int(cb.Conn(), c * n + i);
                if (v[i] < 0 || static_cast<std::size_t>(v[i]) >= npts)
                    throw std::invalid_argument(std::string(pPrefix) +
                                                "a cell references a point the mesh does not have");
            }
            if (line) {
                skin.mHasLines = true;
                skin.mIds.push_back({v[0], v[1]});
                continue;
            }
            for (std::size_t e = 0; e < n; ++e) {
                const std::int64_t a = v[e], b = v[(e + 1) % n];
                const std::array<std::int64_t, 2> key{std::min(a, b), std::max(a, b)};
                auto it = boundary.find(key);
                if (it == boundary.end()) {
                    boundary.emplace(
                        key, std::make_pair(std::int64_t{1}, std::array<std::int64_t, 2>{a, b}));
                    order.push_back(key);
                } else {
                    ++it->second.first;
                }
            }
        }
    }
    for (const auto& key : order) {
        const auto& rec = boundary.at(key);
        if (rec.first == 1)
            skin.mIds.push_back(rec.second);
    }
    for (const auto& id : skin.mIds)
        skin.mSegments.push_back(
            {rFrame.Permute(read_point(points, dim, id[0])),
             rFrame.Permute(read_point(points, dim, id[1]))});
    return skin;
}

bool mesh_is_polyline_only(const Mesh& rMesh) {
    bool any_line = false;
    for (const auto cb : rMesh.CellRange()) {
        const std::string type(cb.Type());
        if (type == "line")
            any_line = true;
        else if (type != "vertex")
            return false;
    }
    return any_line;
}

void require_closed_polyline(const PlanarSkin& rSkin, const char* pPrefix) {
    if (!rSkin.mHasLines)
        return;
    std::map<std::int64_t, int> degree;
    for (const auto& id : rSkin.mIds) {
        ++degree[id[0]];
        ++degree[id[1]];
    }
    std::string ends;
    int count = 0;
    for (const auto& kv : degree)
        if (kv.second % 2 != 0) {
            if (count < 8)
                ends += (ends.empty() ? "" : ", ") + std::to_string(kv.first);
            ++count;
        }
    if (count > 0)
        throw std::invalid_argument(std::string(pPrefix) +
                                    "the polyline does not close, so it has no inside: open ends "
                                    "at points " +
                                    ends + (count > 8 ? ", ..." : ""));
}

bool segment_box_overlap_2d(const Vec3& rCentre, const Vec3& rHalf, const Vec3& rA,
                            const Vec3& rB) {
    const double ax = rA[0] - rCentre[0], ay = rA[1] - rCentre[1];
    const double bx = rB[0] - rCentre[0], by = rB[1] - rCentre[1];
    if (std::min(ax, bx) > rHalf[0] || std::max(ax, bx) < -rHalf[0])
        return false;
    if (std::min(ay, by) > rHalf[1] || std::max(ay, by) < -rHalf[1])
        return false;
    const double dx = bx - ax, dy = by - ay;
    const double nx = -dy, ny = dx;
    const double proj = nx * ax + ny * ay;
    return std::fabs(proj) <= std::fabs(nx) * rHalf[0] + std::fabs(ny) * rHalf[1];
}

Vec3 polyline_closest_point_2d(const PlanarSkin& rSkin, const Vec3& rP) {
    Vec3 best = rP;
    double best_d2 = std::numeric_limits<double>::infinity();
    for (const auto& s : rSkin.mSegments) {
        const double ax = s[0][0], ay = s[0][1];
        const double dx = s[1][0] - ax, dy = s[1][1] - ay;
        const double len2 = dx * dx + dy * dy;
        double t = len2 > 0.0 ? ((rP[0] - ax) * dx + (rP[1] - ay) * dy) / len2 : 0.0;
        t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
        const double qx = ax + t * dx, qy = ay + t * dy;
        const double d2 = (qx - rP[0]) * (qx - rP[0]) + (qy - rP[1]) * (qy - rP[1]);
        if (d2 < best_d2) {
            best_d2 = d2;
            best = Vec3{{qx, qy, rP[2]}};
        }
    }
    return best;
}

double polyline_signed_distance_2d(const PlanarSkin& rSkin, const Vec3& rP, bool EvenOdd) {
    double best = std::numeric_limits<double>::infinity();
    int winding = 0;
    int crossings = 0;
    for (const auto& s : rSkin.mSegments) {
        const double ax = s[0][0], ay = s[0][1], bx = s[1][0], by = s[1][1];
        const double dx = bx - ax, dy = by - ay;
        const double len2 = dx * dx + dy * dy;
        double t = len2 > 0.0 ? ((rP[0] - ax) * dx + (rP[1] - ay) * dy) / len2 : 0.0;
        t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
        const double qx = ax + t * dx - rP[0], qy = ay + t * dy - rP[1];
        const double d = std::sqrt(qx * qx + qy * qy);
        if (d < best)
            best = d;
        const double left = dx * (rP[1] - ay) - (rP[0] - ax) * dy;
        if (ay <= rP[1]) {
            if (by > rP[1] && left > 0.0) {
                ++winding;
                ++crossings;
            }
        } else if (by <= rP[1] && left < 0.0) {
            --winding;
            ++crossings;
        }
    }
    const bool inside = EvenOdd ? (crossings % 2 != 0) : (winding != 0);
    return inside ? -best : best;
}

}  // namespace detail
}  // namespace meshioplusplus
