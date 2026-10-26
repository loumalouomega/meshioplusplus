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
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "../detail/brep_trim.hpp"
#include "../detail/edge_sampler.hpp"
#include "../detail/surface_project.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/project_to_cad.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kPtcPrefix = "meshio++: project_to_cad: ";

[[noreturn]] void ptc_fail(const std::string& rWhy) {
    throw std::invalid_argument(std::string(kPtcPrefix) + rWhy);
}

bool ptc_is_line(const std::string& rType) {
    return rType.rfind("line", 0) == 0;
}

bool ptc_is_face(const std::string& rType) {
    return rType.rfind("triangle", 0) == 0 || rType.rfind("quad", 0) == 0 || rType == "polygon";
}

// What a point lies on: the lowest face and edge its cells name, and
// whether they name more than one.
struct PtcOn {
    std::int64_t mFace = -1;
    std::int64_t mEdge = -1;
    bool mFaces = false;  ///< more than one face
    bool mEdges = false;  ///< more than one edge
};

void ptc_note(std::int64_t& rFirst, bool& rMore, std::int64_t K) {
    if (K < 0)
        return;
    if (rFirst < 0)
        rFirst = K;
    else if (rFirst != K) {
        rMore = true;
        rFirst = std::min(rFirst, K);
    }
}

// The face a topology reference names, by id or by name.
std::optional<std::size_t> ptc_face_of(const BrepInfo& rBrep, const BrepRef& rRef) {
    for (std::size_t f = 0; f < rBrep.mFaces.size(); ++f) {
        const BrepFace& r_f = rBrep.mFaces[f];
        if ((rRef.mId >= 0 && r_f.mId == rRef.mId) ||
            (rRef.mId < 0 && !rRef.mName.empty() && r_f.mName == rRef.mName))
            return f;
    }
    return std::nullopt;
}

// The trim (on a loop or embedded) of face F with index T.
const BrepTrim* ptc_trim(const BrepFace& rFace, std::int64_t T) {
    for (const BrepLoop& r_l : rFace.mLoops)
        for (const BrepTrim& r_t : r_l.mTrims)
            if (r_t.mTrimIndex == T)
                return &r_t;
    for (const BrepTrim& r_t : rFace.mEmbeddedEdges)
        if (r_t.mTrimIndex == T)
            return &r_t;
    return nullptr;
}

}  // namespace

ProjectToCadResult project_to_cad(const Mesh& rMesh, const BrepInfo& rBrep,
                                  const ProjectToCadOptions& rOptions) {
    if (!(std::isfinite(rOptions.mMaxDistance) && rOptions.mMaxDistance >= 0.0))
        ptc_fail("max_distance must be finite and not negative");
    if (rOptions.mTrimSamples < 2 || rOptions.mTrimSamples > 1024)
        ptc_fail("trim_samples must be between 2 and 1024");
    if (!rMesh.HasCellData(rOptions.mFaceArray))
        ptc_fail("the mesh has no cell data '" + rOptions.mFaceArray +
                 "' naming each cell's CAD face (mesh_cad_surface writes it)");
    const NDArray& r_pts = rMesh.Points();
    const std::size_t np = rMesh.NumPoints();
    if (np > 0 && (r_pts.Ndim() != 2 || r_pts.Shape()[1] != 3))
        ptc_fail("the mesh's points must be three-dimensional");
    const bool has_edges = rMesh.HasCellData(rOptions.mEdgeArray);
    const std::size_t nf = rBrep.mFaces.size();
    const std::size_t ne = rBrep.mEdges.size();

    // What each point lies on, from the cells that use it.
    std::vector<PtcOn> on(np);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cells = rMesh.Cells(b);
        const std::string type(cells.Type());
        const bool line = ptc_is_line(type);
        if (!(line || ptc_is_face(type)) || (line && !has_edges))
            continue;
        const NDArray& r_tag = rMesh.CellData(line ? rOptions.mEdgeArray : rOptions.mFaceArray, b);
        const std::size_t limit = line ? ne : nf;
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            const std::int64_t k = detail::read_int(r_tag, c);
            if (k >= 0 && static_cast<std::size_t>(k) >= limit)
                ptc_fail("cell " + std::to_string(c) + " of block " + std::to_string(b) +
                         " names " + (line ? "edge " : "face ") + std::to_string(k) +
                         ", out of range");
            const auto visit = [&](std::int64_t P) {
                PtcOn& r_o = on[static_cast<std::size_t>(P)];
                if (line)
                    ptc_note(r_o.mEdge, r_o.mEdges, k);
                else
                    ptc_note(r_o.mFace, r_o.mFaces, k);
            };
            if (cells.IsRagged()) {
                if (cells.IsPolyhedron())
                    continue;
                for (std::size_t i = 0; i < cells.RowSize(c); ++i)
                    visit(cells.Row(c)[i]);
            } else {
                const std::size_t n = cells.NodesPerCell();
                for (std::size_t i = 0; i < n; ++i)
                    visit(detail::read_int(cells.Conn(), n * c + i));
            }
        }
    }

    // The surfaces and edge traces the points need, built once.
    std::vector<char> face_used(nf, 0), edge_used(ne, 0);
    for (const PtcOn& r_o : on) {
        if (r_o.mEdge >= 0)
            edge_used[static_cast<std::size_t>(r_o.mEdge)] = 1;
        else if (r_o.mFace >= 0 && !r_o.mFaces)
            face_used[static_cast<std::size_t>(r_o.mFace)] = 1;
    }
    ProjectToCadResult result;
    const auto warn = [&](const std::string& rMsg) {
        log::warn("project_to_cad: {}", rMsg);
        result.mWarnings.push_back(rMsg);
    };
    std::vector<std::unique_ptr<detail::BrepSurfaceEval>> surface(nf);
    const auto surface_of = [&](std::size_t F) -> const detail::BrepSurfaceEval* {
        if (surface[F])
            return surface[F].get();
        const BrepFace& r_f = rBrep.mFaces[F];
        if (r_f.mPatch < 0 || static_cast<std::size_t>(r_f.mPatch) >= rBrep.mPatches.size())
            return nullptr;
        try {
            surface[F] = std::make_unique<detail::BrepSurfaceEval>(
                rBrep.mPatches[static_cast<std::size_t>(r_f.mPatch)], r_f.mRational);
        } catch (const std::invalid_argument&) {
            return nullptr;
        }
        return surface[F].get();
    };
    std::vector<std::optional<detail::SurfaceProjector>> projector(nf);
    std::vector<std::vector<std::vector<std::array<double, 2>>>> loops(nf);
    for (std::size_t f = 0; f < nf; ++f) {
        if (!face_used[f])
            continue;
        const detail::BrepSurfaceEval* p_s = surface_of(f);
        if (!p_s) {
            warn("face " + std::to_string(f) + " has no usable surface");
            continue;
        }
        projector[f].emplace(*p_s,
                             rBrep.mPatches[static_cast<std::size_t>(rBrep.mFaces[f].mPatch)]);
        if (rBrep.mFaces[f].mTrimmed)
            for (const BrepLoop& r_l : rBrep.mFaces[f].mLoops)
                loops[f].push_back(detail::brep_loop_polygon(r_l, rOptions.mTrimSamples));
    }
    std::vector<std::optional<detail::TrimTrace>> trace(ne);
    for (std::size_t e = 0; e < ne; ++e) {
        if (!edge_used[e])
            continue;
        for (const BrepRef& r_r : rBrep.mEdges[e].mTopology) {
            const auto f = ptc_face_of(rBrep, r_r);
            if (!f)
                continue;
            const BrepTrim* p_t = ptc_trim(rBrep.mFaces[*f], r_r.mTrimIndex);
            const detail::BrepSurfaceEval* p_s = surface_of(*f);
            if (!p_t || !p_s)
                continue;
            try {
                trace[e].emplace(p_s, detail::BrepCurveEval(p_t->mCurve), true);
                break;
            } catch (const std::invalid_argument&) {
            }
        }
        if (!trace[e])
            warn("edge " + std::to_string(e) + " has no usable trim");
    }

    // Each point, independently.
    enum class Fate : std::uint8_t { Stay, Face, Edge, Corner, Ambiguous, TooFar, Failed };
    std::vector<Fate> fate(np, Fate::Stay);
    std::vector<char> outside(np, 0);
    std::vector<double> moved(np, 0.0);
    NDArray points = detail::data_owned_copy(r_pts);
    double* px = points.As<double>();
    const double* p0 = r_pts.Size() ? r_pts.As<double>() : nullptr;
    NDArray input;
    if (np > 0 && r_pts.Dtype() != DType::Float64) {
        input = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
        for (std::size_t i = 0; i < 3 * np; ++i)
            input.As<double>()[i] = detail::read_double(r_pts, i);
        p0 = input.As<double>();
        points = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
        px = points.As<double>();
        std::copy(p0, p0 + 3 * np, px);
    }
    parallel_for(
        np,
        [&](std::size_t i) {
            const PtcOn& r_o = on[i];
            const detail::EdgeVec3 p{p0[3 * i], p0[3 * i + 1], p0[3 * i + 2]};
            detail::EdgeVec3 x{};
            bool on_face = false;
            std::size_t f = 0;
            std::array<double, 2> uv{0.0, 0.0};
            if (r_o.mEdge >= 0) {
                const auto& r_t = trace[static_cast<std::size_t>(r_o.mEdge)];
                if (!r_t) {
                    fate[i] = Fate::Failed;
                    return;
                }
                x = r_t->Point(r_t->ParamAtLength(detail::edge_closest_length(*r_t, p)));
                fate[i] = r_o.mEdges ? Fate::Corner : Fate::Edge;
            } else if (r_o.mFace >= 0) {
                if (r_o.mFaces) {
                    fate[i] = Fate::Ambiguous;
                    return;
                }
                f = static_cast<std::size_t>(r_o.mFace);
                if (!projector[f]) {
                    fate[i] = Fate::Failed;
                    return;
                }
                uv = projector[f]->Project(p, x);
                fate[i] = Fate::Face;
                on_face = true;
            } else {
                return;
            }
            const double d =
                std::sqrt((x[0] - p[0]) * (x[0] - p[0]) + (x[1] - p[1]) * (x[1] - p[1]) +
                          (x[2] - p[2]) * (x[2] - p[2]));
            if (rOptions.mMaxDistance > 0.0 && d > rOptions.mMaxDistance) {
                fate[i] = Fate::TooFar;
                return;
            }
            for (std::size_t k = 0; k < 3; ++k)
                px[3 * i + k] = x[k];
            moved[i] = d;
            if (on_face && !loops[f].empty() && !detail::brep_inside(loops[f], uv[0], uv[1]))
                outside[i] = 1;
        },
        64);
    for (std::size_t i = 0; i < np; ++i) {
        switch (fate[i]) {
            case Fate::Face:
            case Fate::Edge:
            case Fate::Corner:
                ++result.mNumProjected;
                result.mNumOnEdges += fate[i] != Fate::Face ? 1 : 0;
                result.mNumCorners += fate[i] == Fate::Corner ? 1 : 0;
                result.mMaxMove = std::max(result.mMaxMove, moved[i]);
                break;
            case Fate::Ambiguous:
                ++result.mNumAmbiguous;
                break;
            case Fate::TooFar:
                ++result.mNumTooFar;
                break;
            case Fate::Failed:
                ++result.mNumFailed;
                break;
            case Fate::Stay:
                break;
        }
        if (outside[i]) {
            ++result.mNumOutsideTrim;
            result.mOutsideTrim.push_back(static_cast<std::int64_t>(i));
        }
    }
    if (result.mNumOutsideTrim > 0)
        warn(std::to_string(result.mNumOutsideTrim) +
             " points projected outside their face's trimmed region");

    // Only the coordinates change.
    Mesh& r_out = result.mMesh;
    r_out = detail::clone_mesh(rMesh);
    r_out.AssignPoints(std::move(points));
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        r_out.AddRegion(rMesh.Region(i));
    for (std::size_t i = 0; i < rMesh.NumPropertySets(); ++i)
        r_out.AddPropertySet(rMesh.GetPropertySet(i));
    return result;
}

}  // namespace meshioplusplus
