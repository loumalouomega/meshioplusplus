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
// Offset, thicken and hollow (roadmap 7.8.4). See operations/offset_surface.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/offset_surface.hpp"
#include "../detail/math_utils.hpp"
#include "../detail/offset_normals.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {

namespace {

using OfsV3 = std::array<double, 3>;
using OfsFace = std::array<std::int64_t, 4>;

[[noreturn]] void ofs_fail(const char* pOp, const std::string& rWhy) {
    throw std::invalid_argument(std::string("meshio++: ") + pOp + ": " + rWhy);
}

// The surface's points (packed xyz) and faces (a triangle's fourth entry -1),
// with the block each face came from.
struct OfsSurface {
    std::vector<double> mX;
    std::vector<OfsFace> mFaces;
};

OfsSurface ofs_read(const Mesh& rMesh, const char* pOp) {
    OfsSurface s;
    const std::size_t np = rMesh.NumPoints();
    const NDArray& r_pts = rMesh.Points();
    if (np > 0 && (r_pts.Ndim() != 2 || r_pts.Shape()[1] != 3))
        ofs_fail(pOp, "the points must be three-dimensional");
    s.mX.resize(3 * np);
    for (std::size_t i = 0; i < 3 * np; ++i)
        s.mX[i] = detail::read_double(r_pts, i);
    for (const auto cb : rMesh.CellRange()) {
        const std::string type(cb.Type());
        if ((type != "triangle" && type != "quad") || cb.IsRagged())
            ofs_fail(pOp, "a '" + type + "' block: only triangles and quads are offset " +
                              "(extract_surface or linearize first)");
        const std::size_t k = type == "triangle" ? 3 : 4;
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            OfsFace f{-1, -1, -1, -1};
            for (std::size_t i = 0; i < k; ++i) {
                f[i] = detail::read_int(cb.Conn(), k * c + i);
                if (f[i] < 0 || static_cast<std::size_t>(f[i]) >= np)
                    ofs_fail(pOp, "a cell names a point out of range");
            }
            s.mFaces.push_back(f);
        }
    }
    if (s.mFaces.empty())
        ofs_fail(pOp, "the mesh has no triangle or quad");
    return s;
}

// A face's Newell normal (not unit) over the points X.
OfsV3 ofs_face_normal(const std::vector<double>& rX, const OfsFace& rF) {
    const int k = rF[3] >= 0 ? 4 : 3;
    OfsV3 n{0.0, 0.0, 0.0};
    for (int i = 0; i < k; ++i) {
        const auto a = static_cast<std::size_t>(rF[static_cast<std::size_t>(i)]);
        const auto b = static_cast<std::size_t>(rF[static_cast<std::size_t>((i + 1) % k)]);
        n[0] += (rX[3 * a + 1] - rX[3 * b + 1]) * (rX[3 * a + 2] + rX[3 * b + 2]);
        n[1] += (rX[3 * a + 2] - rX[3 * b + 2]) * (rX[3 * a] + rX[3 * b]);
        n[2] += (rX[3 * a] - rX[3 * b]) * (rX[3 * a + 1] + rX[3 * b + 1]);
    }
    return n;
}

// The points moved by Distance along the offset directions; the counts go
// into rResult.
std::vector<double> ofs_move(const OfsSurface& rS, const detail::OffsetNormals& rN, double Distance,
                             OffsetSurfaceResult& rResult) {
    std::vector<double> y = rS.mX;
    const std::size_t np = rS.mX.size() / 3;
    for (std::size_t p = 0; p < np; ++p) {
        if (!rN.mUsed[p])
            continue;
        double stretch = rN.mStretch[p];
        if (!(stretch > 0.0)) {
            stretch = 1.0;
            ++rResult.mNumNotVisible;
        } else if (stretch >= detail::kOffsetMaxStretch) {
            ++rResult.mNumClamped;
        }
        const double move = Distance * stretch;
        for (std::size_t c = 0; c < 3; ++c)
            y[3 * p + c] += move * rN.mDirections[p][c];
        rResult.mMaxMove = std::max(rResult.mMaxMove, std::abs(move));
    }
    for (std::size_t f = 0; f < rS.mFaces.size(); ++f)
        if (detail::vec3_dot(ofs_face_normal(rS.mX, rS.mFaces[f]),
                             ofs_face_normal(y, rS.mFaces[f])) <= 0.0)
            ++rResult.mNumFolded;
    return y;
}

void ofs_warn(OffsetSurfaceResult& rResult, const char* pOp) {
    if (rResult.mNumNotVisible > 0)
        rResult.mWarnings.push_back(std::string(pOp) + ": " +
                                    std::to_string(rResult.mNumNotVisible) +
                                    " point(s) have no direction every face around them sees; "
                                    "they moved by the distance alone");
    if (rResult.mNumFolded > 0)
        rResult.mWarnings.push_back(std::string(pOp) + ": " + std::to_string(rResult.mNumFolded) +
                                    " face(s) turned over: the offset is deeper than a concave "
                                    "region's radius");
}

NDArray ofs_points(const std::vector<double>& rX) {
    NDArray p = NDArray::Uninit(DType::Float64, {rX.size() / 3, std::size_t{3}});
    std::copy(rX.begin(), rX.end(), p.As<double>());
    return p;
}

// A shell: the upper layer U (kept), the lower layer L (reversed), with the
// rim quads when Rim; its point data duplicated, thicken:side attached.
void ofs_shell(const Mesh& rMesh, const OfsSurface& rS, const std::vector<double>& rUpper,
               const std::vector<double>& rLower, bool Rim, const char* pOp,
               OffsetSurfaceResult& rResult) {
    const auto np = static_cast<std::int64_t>(rS.mX.size() / 3);
    std::vector<double> xyz = rUpper;
    xyz.insert(xyz.end(), rLower.begin(), rLower.end());
    std::vector<std::int64_t> tris, quads, tri_side, quad_side;
    for (int layer = 0; layer < 2; ++layer)
        for (const OfsFace& r_f : rS.mFaces) {
            const std::int64_t off = layer == 0 ? 0 : np;
            if (r_f[3] < 0) {
                if (layer == 0)
                    tris.insert(tris.end(), {r_f[0], r_f[1], r_f[2]});
                else
                    tris.insert(tris.end(), {r_f[0] + off, r_f[2] + off, r_f[1] + off});
                tri_side.push_back(layer);
            } else {
                if (layer == 0)
                    quads.insert(quads.end(), {r_f[0], r_f[1], r_f[2], r_f[3]});
                else
                    quads.insert(quads.end(),
                                 {r_f[0] + off, r_f[3] + off, r_f[2] + off, r_f[1] + off});
                quad_side.push_back(layer);
            }
        }
    // The boundary edges, in their face's direction, in face order.
    std::map<std::pair<std::int64_t, std::int64_t>, int> uses;
    for (const OfsFace& r_f : rS.mFaces) {
        const int k = r_f[3] >= 0 ? 4 : 3;
        for (int i = 0; i < k; ++i) {
            const std::int64_t a = r_f[static_cast<std::size_t>(i)];
            const std::int64_t b = r_f[static_cast<std::size_t>((i + 1) % k)];
            ++uses[{std::min(a, b), std::max(a, b)}];
        }
    }
    std::int64_t num_boundary = 0;
    for (const OfsFace& r_f : rS.mFaces) {
        const int k = r_f[3] >= 0 ? 4 : 3;
        for (int i = 0; i < k; ++i) {
            const std::int64_t a = r_f[static_cast<std::size_t>(i)];
            const std::int64_t b = r_f[static_cast<std::size_t>((i + 1) % k)];
            if (uses[{std::min(a, b), std::max(a, b)}] != 1)
                continue;
            ++num_boundary;
            if (Rim) {
                quads.insert(quads.end(), {b, a, a + np, b + np});
                quad_side.push_back(2);
            }
        }
    }
    if (!Rim && num_boundary > 0)
        ofs_fail(pOp, "the surface has " + std::to_string(num_boundary) +
                          " boundary edge(s): hollow takes a closed surface (thicken an open one)");
    rResult.mNumRimQuads = Rim ? num_boundary : 0;

    Mesh& out = rResult.mMesh;
    out.AssignPoints(ofs_points(xyz));
    std::vector<NDArray> side;
    const auto add_block = [&](const char* pType, const std::vector<std::int64_t>& rConn,
                               const std::vector<std::int64_t>& rSide, std::size_t K) {
        if (rSide.empty())
            return;
        NDArray c = NDArray::Uninit(DType::Int64, {rSide.size(), K});
        std::copy(rConn.begin(), rConn.end(), c.As<std::int64_t>());
        out.AddCellBlock(pType, std::move(c));
        NDArray s = NDArray::Uninit(DType::Int64, {rSide.size()});
        std::copy(rSide.begin(), rSide.end(), s.As<std::int64_t>());
        side.push_back(std::move(s));
    };
    add_block("triangle", tris, tri_side, 3);
    add_block("quad", quads, quad_side, 4);
    std::vector<std::int64_t> twice(2 * static_cast<std::size_t>(np));
    for (std::size_t i = 0; i < twice.size(); ++i)
        twice[i] = static_cast<std::int64_t>(i % static_cast<std::size_t>(np));
    for (const std::string& r_name : rMesh.PointDataNames())
        out.AddPointData(r_name, detail::subset_gather_rows(rMesh.PointData(r_name), twice));
    out.AddCellData("thicken:side", std::move(side));
    for (const std::string& r_name : rMesh.FieldDataNames())
        out.AddFieldData(r_name, rMesh.FieldData(r_name));
    if (!rMesh.CellDataNames().empty() || rMesh.NumRegions() > 0)
        rResult.mWarnings.push_back(std::string(pOp) +
                                    ": the input's cell data and regions are not carried "
                                    "onto the shell");
}

}  // namespace

ThickenSide thicken_side_from_name(std::string_view Name) {
    if (Name == "inward")
        return ThickenSide::Inward;
    if (Name == "outward")
        return ThickenSide::Outward;
    if (Name == "both")
        return ThickenSide::Both;
    throw std::invalid_argument("meshio++: thicken: side must be inward, outward or both, not '" +
                                std::string(Name) + "'");
}

OffsetSurfaceResult offset_surface(const Mesh& rMesh, const OffsetSurfaceOptions& rOptions) {
    if (!std::isfinite(rOptions.mDistance))
        ofs_fail("offset_surface", "distance must be finite");
    if (rOptions.mSmoothIterations < 0 || rOptions.mSmoothIterations > 1000)
        ofs_fail("offset_surface", "smooth_iterations must be between 0 and 1000");
    const OfsSurface s = ofs_read(rMesh, "offset_surface");
    const detail::OffsetNormals n =
        detail::offset_normals(s.mX, s.mFaces, 1.0, rOptions.mSmoothIterations);
    OffsetSurfaceResult result;
    const std::vector<double> y = ofs_move(s, n, rOptions.mDistance, result);
    result.mMesh = detail::clone_mesh(rMesh);
    result.mMesh.AssignPoints(ofs_points(y));
    ofs_warn(result, "offset_surface");
    return result;
}

OffsetSurfaceResult thicken(const Mesh& rMesh, const ThickenOptions& rOptions) {
    const double t = rOptions.mThickness;
    if (!(std::isfinite(t) && t > 0.0))
        ofs_fail("thicken", "thickness must be positive and finite");
    if (rOptions.mSmoothIterations < 0 || rOptions.mSmoothIterations > 1000)
        ofs_fail("thicken", "smooth_iterations must be between 0 and 1000");
    if (rOptions.mSide != ThickenSide::Inward && rOptions.mSide != ThickenSide::Outward &&
        rOptions.mSide != ThickenSide::Both)
        ofs_fail("thicken", "side must be inward, outward or both");
    const OfsSurface s = ofs_read(rMesh, "thicken");
    const detail::OffsetNormals n =
        detail::offset_normals(s.mX, s.mFaces, 1.0, rOptions.mSmoothIterations);
    OffsetSurfaceResult result;
    const double hi = rOptions.mSide == ThickenSide::Outward  ? t
                      : rOptions.mSide == ThickenSide::Inward ? 0.0
                                                              : 0.5 * t;
    const double lo = hi - t;
    const std::vector<double> upper = hi == 0.0 ? s.mX : ofs_move(s, n, hi, result);
    const std::vector<double> lower = lo == 0.0 ? s.mX : ofs_move(s, n, lo, result);
    ofs_shell(rMesh, s, upper, lower, true, "thicken", result);
    ofs_warn(result, "thicken");
    return result;
}

OffsetSurfaceResult hollow(const Mesh& rMesh, double Wall, int SmoothIterations) {
    if (!(std::isfinite(Wall) && Wall > 0.0))
        ofs_fail("hollow", "wall must be positive and finite");
    if (SmoothIterations < 0 || SmoothIterations > 1000)
        ofs_fail("hollow", "smooth_iterations must be between 0 and 1000");
    const OfsSurface s = ofs_read(rMesh, "hollow");
    const detail::OffsetNormals n = detail::offset_normals(s.mX, s.mFaces, 1.0, SmoothIterations);
    OffsetSurfaceResult result;
    const std::vector<double> inner = ofs_move(s, n, -Wall, result);
    ofs_shell(rMesh, s, s.mX, inner, false, "hollow", result);
    ofs_warn(result, "hollow");
    return result;
}

}  // namespace meshioplusplus
