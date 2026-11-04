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
// Field-aligned quad remeshing through a quad remesher provider. See
// operations/remesh_quads.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/remesh_quads.hpp"
#include "meshioplusplus/providers.hpp"

// Project includes (private, not installed)
#include "../detail/quad_recombine.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kRqPrefix = "meshio++: remesh_quads: ";

[[noreturn]] void rq_fail(const std::string& rWhy) {
    throw std::invalid_argument(std::string(kRqPrefix) + rWhy);
}

std::string rq_join(const std::vector<std::string>& rNames) {
    std::string out;
    for (const std::string& r_name : rNames)
        out += (out.empty() ? "" : ", ") + r_name;
    return out.empty() ? std::string("none") : out;
}

}  // namespace

std::vector<std::string> remesh_quads_backends() {
    return quad_remesher_names();
}

RemeshQuadsResult remesh_quads(const Mesh& rSurface, const RemeshQuadsOptions& rOptions) {
    const RemeshQuadsOptions& o = rOptions;
    if (o.mMaxCells < 1)
        rq_fail("max_cells must be positive");
    if (o.mTargetFaces < 0)
        rq_fail("target_faces must not be negative");
    if (o.mTargetFaces > o.mMaxCells)
        rq_fail("target_faces " + std::to_string(o.mTargetFaces) + " is past max_cells (" +
                std::to_string(o.mMaxCells) + ")");

    // The triangles, block by block.
    QuadRemeshRequest request;
    const std::size_t np = rSurface.NumPoints();
    const NDArray& r_pts = rSurface.Points();
    if (np > 0 && (r_pts.Ndim() != 2 || r_pts.Shape()[1] != 3))
        rq_fail("the surface's points must be three-dimensional");
    request.mPoints.resize(np);
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t k = 0; k < 3; ++k) {
            request.mPoints[i][k] = detail::read_double(r_pts, 3 * i + k);
            if (!std::isfinite(request.mPoints[i][k]))
                rq_fail("the surface has a non-finite point");
        }
    for (std::size_t b = 0; b < rSurface.NumCellBlocks(); ++b) {
        const auto cells = rSurface.Cells(b);
        const std::string type(cells.Type());
        if (type.rfind("line", 0) == 0 || type == "vertex")
            continue;
        if (type != "triangle" || cells.IsRagged())
            rq_fail("a '" + type +
                    "' block: the surface must be linear triangles (convert_cells, linearize or "
                    "extract_surface first)");
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            std::array<std::int64_t, 3> t{};
            for (std::size_t k = 0; k < 3; ++k) {
                t[k] = detail::read_int(cells.Conn(), 3 * c + k);
                if (t[k] < 0 || static_cast<std::size_t>(t[k]) >= np)
                    rq_fail("a triangle names a point out of range");
            }
            request.mTriangles.push_back(t);
        }
    }
    if (request.mTriangles.empty())
        rq_fail("the surface has no triangle");
    request.mTargetFaces =
        o.mTargetFaces > 0
            ? o.mTargetFaces
            : std::max<std::int64_t>(1, static_cast<std::int64_t>(request.mTriangles.size() / 2));
    request.mSharp = o.mSharp;
    request.mBoundary = o.mBoundary;
    request.mAdaptive = o.mAdaptive;
    request.mDominant = o.mDominant;
    request.mSeed = o.mSeed;

    // The remesher.
    std::shared_ptr<const QuadRemesherProvider> p_remesher;
    if (o.mBackend == "auto" || o.mBackend.empty()) {
        for (const char* p_name : {"quadriflow", "instant_meshes"})
            if (!p_remesher)
                p_remesher = quad_remesher(p_name);
        if (!p_remesher) {
            const std::vector<std::string> names = quad_remesher_names();
            if (!names.empty())
                p_remesher = quad_remesher(names.front());
        }
        if (!p_remesher)
            throw Unsupported(std::string(kRqPrefix) +
                              "no quad remesher is installed: install QuadriFlow (`quadriflow`) "
                              "or Instant Meshes (`Instant Meshes`) on PATH, or name them "
                              "through MESHIOPLUSPLUS_QUADRIFLOW and "
                              "MESHIOPLUSPLUS_INSTANT_MESHES (roadmap 7.7.3)");
    } else {
        p_remesher = quad_remesher(o.mBackend);
        if (!p_remesher)
            throw Unsupported(std::string(kRqPrefix) + "no quad remesher named '" + o.mBackend +
                              "' (available: auto, " + rq_join(quad_remesher_names()) +
                              "; quadriflow and instant_meshes are found on PATH, or through "
                              "MESHIOPLUSPLUS_QUADRIFLOW and MESHIOPLUSPLUS_INSTANT_MESHES)");
    }
    QuadRemeshResult out;
    try {
        out = p_remesher->mRemesh(request);
    } catch (const std::runtime_error& rErr) {
        throw std::runtime_error(std::string(kRqPrefix) + p_remesher->mName + ": " + rErr.what());
    }
    const std::size_t nq = out.mQuads.size(), nt = out.mTriangles.size(), nop = out.mPoints.size();
    if (nq + nt == 0)
        throw std::runtime_error(std::string(kRqPrefix) + p_remesher->mName + " returned no face");
    if (static_cast<std::int64_t>(nq + nt) > o.mMaxCells)
        throw std::length_error(std::string(kRqPrefix) + "the result would have more than " +
                                std::to_string(o.mMaxCells) + " faces (max_cells)");
    const auto check = [&](std::int64_t V) {
        if (V < 0 || static_cast<std::size_t>(V) >= nop)
            throw std::runtime_error(std::string(kRqPrefix) + p_remesher->mName +
                                     " returned a face naming no point");
    };
    for (const auto& r_q : out.mQuads)
        for (const std::int64_t v : r_q)
            check(v);
    for (const auto& r_t : out.mTriangles)
        for (const std::int64_t v : r_t)
            check(v);
    for (const auto& r_p : out.mPoints)
        for (const double c : r_p)
            if (!std::isfinite(c))
                throw std::runtime_error(std::string(kRqPrefix) + p_remesher->mName +
                                         " returned a non-finite point");

    // The mesh.
    RemeshQuadsResult result;
    result.mBackend = p_remesher->mName;
    NDArray p = NDArray::Uninit(DType::Float64, {nop, std::size_t{3}});
    double* const p_x = p.As<double>();
    for (std::size_t i = 0; i < nop; ++i)
        for (std::size_t k = 0; k < 3; ++k)
            p_x[3 * i + k] = out.mPoints[i][k];
    result.mMesh.AssignPoints(std::move(p));
    if (nq) {
        NDArray q = NDArray::Uninit(DType::Int64, {nq, std::size_t{4}});
        for (std::size_t c = 0; c < nq; ++c)
            for (std::size_t k = 0; k < 4; ++k)
                q.As<std::int64_t>()[4 * c + k] = out.mQuads[c][k];
        result.mMesh.AddCellBlock("quad", std::move(q));
    }
    if (nt) {
        NDArray t = NDArray::Uninit(DType::Int64, {nt, std::size_t{3}});
        for (std::size_t c = 0; c < nt; ++c)
            for (std::size_t k = 0; k < 3; ++k)
                t.As<std::int64_t>()[3 * c + k] = out.mTriangles[c][k];
        result.mMesh.AddCellBlock("triangle", std::move(t));
        result.mWarnings.push_back(std::to_string(nt) + " triangle(s) left among the quads");
    }
    const double unit = detail::mesh_length_unit(rSurface);
    if (unit > 0.0)
        detail::set_mesh_length_unit(result.mMesh, unit);

    // The measures.
    result.mNumQuads = static_cast<std::int64_t>(nq);
    result.mNumTriangles = static_cast<std::int64_t>(nt);
    result.mNumPoints = static_cast<std::int64_t>(nop);
    if (nq) {
        double lo = 1.0, sum = 0.0;
        for (const auto& r_q : out.mQuads) {
            const double j = detail::quad_scaled_jacobian(p_x, r_q);
            lo = std::min(lo, j);
            sum += j;
        }
        result.mMinScaledJacobian = lo;
        result.mMeanScaledJacobian = sum / static_cast<double>(nq);
    }
    Mesh input;
    {
        NDArray ip = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t k = 0; k < 3; ++k)
                ip.As<double>()[3 * i + k] = request.mPoints[i][k];
        input.AssignPoints(std::move(ip));
        NDArray it = NDArray::Uninit(DType::Int64, {request.mTriangles.size(), std::size_t{3}});
        for (std::size_t c = 0; c < request.mTriangles.size(); ++c)
            for (std::size_t k = 0; k < 3; ++k)
                it.As<std::int64_t>()[3 * c + k] = request.mTriangles[c][k];
        input.AddCellBlock("triangle", std::move(it));
    }
    HausdorffOptions h;
    h.mFaceSamples = 2;
    result.mHausdorffToSurface = hausdorff_distance(result.mMesh, input, h).mDistance;
    for (const std::string& r_w : result.mWarnings)
        log::warn("{}{}", kRqPrefix, r_w);
    return result;
}

}  // namespace meshioplusplus
