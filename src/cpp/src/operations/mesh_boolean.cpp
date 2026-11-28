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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/mesh_boolean.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/providers.hpp"

// Core-private includes
#include "../detail/surface_rings.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kBoolPrefix = "meshio++: mesh_boolean: ";

std::string bool_lower(std::string s) {
    for (char& c : s)
        if (c >= 'A' && c <= 'Z')
            c = static_cast<char>(c - 'A' + 'a');
    return s;
}

void bool_extract_triangles(const Mesh& rMesh, std::vector<std::array<double, 3>>& rPts,
                            std::vector<std::array<std::int64_t, 3>>& rTris) {
    const std::vector<double> xyz = detail::surface_xyz(rMesh);
    const std::size_t npts = xyz.size() / 3;
    rPts.resize(npts);
    for (std::size_t i = 0; i < npts; ++i)
        rPts[i] = {xyz[3 * i + 0], xyz[3 * i + 1], xyz[3 * i + 2]};

    std::size_t total = 0;
    for (const auto cb : rMesh.CellRange())
        total += cb.NumCells();
    std::vector<std::uint8_t> sel(total, 1);
    const detail::SurfaceRings rings = detail::surface_rings(rMesh, sel, kBoolPrefix);
    if (rings.NumFaces() == 0)
        throw std::invalid_argument(std::string(kBoolPrefix) + "mesh contains no surface cells");
    for (std::size_t f = 0; f < rings.NumFaces(); ++f) {
        const std::int64_t s0 = rings.mStart[f];
        const std::int64_t s1 = rings.mStart[f + 1];
        const std::size_t m = static_cast<std::size_t>(s1 - s0);
        if (m < 3)
            continue;
        const std::int64_t v0 = rings.mNodes[static_cast<std::size_t>(s0)];
        for (std::size_t k = 1; k + 1 < m; ++k) {
            rTris.push_back({v0, rings.mNodes[static_cast<std::size_t>(s0) + k],
                             rings.mNodes[static_cast<std::size_t>(s0) + k + 1]});
        }
    }
}

}  // namespace

BooleanOp boolean_op_from_name(const std::string& rName) {
    const std::string low = bool_lower(rName);
    if (low == "union" || low == "or" || low == "+")
        return BooleanOp::Union;
    if (low == "difference" || low == "subtract" || low == "minus" || low == "-")
        return BooleanOp::Difference;
    if (low == "intersection" || low == "intersect" || low == "and" || low == "&")
        return BooleanOp::Intersection;
    throw std::invalid_argument(std::string(kBoolPrefix) + "unknown op '" + rName +
                                "' (expected 'union', 'difference', or 'intersection')");
}

bool has_boolean_provider() noexcept {
    return boolean_provider() != nullptr;
}

MeshBooleanResult mesh_boolean(const Mesh& rA, const Mesh& rB, const MeshBooleanOptions& rOptions) {
    auto prov = boolean_provider();
    if (!prov || !prov->mBoolean)
        throw Unsupported(
            std::string(kBoolPrefix) +
            "no BooleanProvider is registered (build with -DMESHIOPLUSPLUS_WITH_MANIFOLD=ON "
            "or register a BooleanProvider plugin)");
    if (!rOptions.mBackend.empty() && rOptions.mBackend != "auto" &&
        rOptions.mBackend != prov->mName)
        throw std::invalid_argument(std::string(kBoolPrefix) + "requested backend '" +
                                    rOptions.mBackend + "' is not registered (active: '" +
                                    prov->mName + "')");

    BooleanRequest req;
    req.mOp = static_cast<std::int32_t>(rOptions.mOp);
    bool_extract_triangles(rA, req.mPointsA, req.mTrianglesA);
    bool_extract_triangles(rB, req.mPointsB, req.mTrianglesB);

    const BooleanReply rep = prov->mBoolean(req);

    const std::size_t npts = rep.mPoints.size();
    const std::size_t ntris = rep.mTriangles.size();
    NDArray pts = NDArray::Uninit(DType::Float64, {npts, std::size_t{3}});
    auto* pp = pts.As<double>();
    for (std::size_t i = 0; i < npts; ++i) {
        pp[3 * i + 0] = rep.mPoints[i][0];
        pp[3 * i + 1] = rep.mPoints[i][1];
        pp[3 * i + 2] = rep.mPoints[i][2];
    }
    NDArray conn = NDArray::Uninit(DType::Int64, {ntris, std::size_t{3}});
    auto* cp = conn.As<std::int64_t>();
    for (std::size_t i = 0; i < ntris; ++i) {
        cp[3 * i + 0] = rep.mTriangles[i][0];
        cp[3 * i + 1] = rep.mTriangles[i][1];
        cp[3 * i + 2] = rep.mTriangles[i][2];
    }

    Mesh out;
    out.AssignPoints(std::move(pts));
    if (ntris > 0)
        out.AddCellBlock("triangle", std::move(conn));

    MeshBooleanResult res;
    res.mMesh = std::move(out);
    res.mBackend = prov->mName;
    res.mNumTriangles = static_cast<std::int64_t>(ntris);
    return res;
}

}  // namespace meshioplusplus
