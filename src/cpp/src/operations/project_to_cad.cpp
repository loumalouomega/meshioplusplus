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

// A point on a trim lies on it to rounding, which the even-odd test cannot
// decide, so the trim's tolerance is a relative 1e-9 of the loops' extent.
double ptc_trim_tolerance(const std::vector<std::vector<std::array<double, 2>>>& rLoops) {
    constexpr double kInf = std::numeric_limits<double>::infinity();
    double lo[2] = {kInf, kInf};
    double hi[2] = {-kInf, -kInf};
    for (const auto& poly : rLoops)
        for (const auto& q : poly)
            for (std::size_t k = 0; k < 2; ++k) {
                lo[k] = std::min(lo[k], q[k]);
                hi[k] = std::max(hi[k], q[k]);
            }
    if (hi[0] < lo[0])
        return 0.0;
    return 1e-9 * std::hypot(hi[0] - lo[0], hi[1] - lo[1]);
}

// Whether (U, V) lies within `Tol` of an edge of one of the loops.
bool ptc_on_trim(const std::vector<std::vector<std::array<double, 2>>>& rLoops, double U,
                 double V, double Tol) {
    for (const auto& poly : rLoops) {
        const std::size_t n = poly.size();
        for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
            const double ex = poly[i][0] - poly[j][0];
            const double ey = poly[i][1] - poly[j][1];
            const double len2 = ex * ex + ey * ey;
            const double s = (U - poly[j][0]) * ex + (V - poly[j][1]) * ey;
            const double t = len2 > 0.0 ? std::clamp(s / len2, 0.0, 1.0) : 0.0;
            const double dx = poly[j][0] + t * ex - U;
            const double dy = poly[j][1] + t * ey - V;
            if (dx * dx + dy * dy <= Tol * Tol)
                return true;
        }
    }
    return false;
}

// The nearest point of a 3-D curve to `rP`: the best of 64 equal steps of
// its parameter, then Gauss-Newton on `(C - p) . C' = 0` within the steps
// beside it.
detail::EdgeVec3 ptc_curve_closest(const detail::BrepCurveEval& rC, const detail::EdgeVec3& rP) {
    constexpr int kSteps = 64;
    const double t0 = rC.T0(), t1 = rC.T1();
    const auto dist2 = [&](double T, detail::EdgeVec3& rX) {
        double d[3];
        rC.Eval(T, rX.data(), d);
        return (rX[0] - rP[0]) * (rX[0] - rP[0]) + (rX[1] - rP[1]) * (rX[1] - rP[1]) +
               (rX[2] - rP[2]) * (rX[2] - rP[2]);
    };
    double best_t = t0, best = std::numeric_limits<double>::infinity();
    detail::EdgeVec3 x{};
    for (int i = 0; i <= kSteps; ++i) {
        const double t = t0 + (t1 - t0) * double(i) / kSteps;
        const double d = dist2(t, x);
        if (d < best) {
            best = d;
            best_t = t;
        }
    }
    const double step = (t1 - t0) / kSteps;
    const double lo = std::max(t0, best_t - step), hi = std::min(t1, best_t + step);
    double t = best_t;
    for (int it = 0; it < 30; ++it) {
        double d[3];
        rC.Eval(t, x.data(), d);
        const double g = (x[0] - rP[0]) * d[0] + (x[1] - rP[1]) * d[1] + (x[2] - rP[2]) * d[2];
        const double h = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        if (!(h > 0.0))
            break;
        const double next = std::clamp(t - g / h, lo, hi);
        if (std::fabs(next - t) <= 1e-15 * (1.0 + std::fabs(t)))
            break;
        t = next;
    }
    detail::EdgeVec3 y{};
    if (dist2(t, y) <= best)
        return y;
    dist2(best_t, y);
    return y;
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
    std::vector<double> trim_tol(nf, 0.0);
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
        trim_tol[f] = ptc_trim_tolerance(loops[f]);
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
    // An edge's own 3-D curve, where the model gives one that its trace
    // follows (within a millionth of its length at nine points, or a
    // billionth of the model's size): a point on the edge goes onto it, the
    // trace's surface error left out. A curve the trace leaves is warned and
    // the trace kept.
    std::vector<std::optional<detail::BrepCurveEval>> curve(ne);
    for (std::size_t e = 0; e < ne; ++e) {
        const BrepEdge& r_e = rBrep.mEdges[e];
        if (!edge_used[e] || !trace[e] || !r_e.mHasCurve || r_e.mCurve.mControlPoints.Ndim() != 2 ||
            r_e.mCurve.mControlPoints.Shape()[1] != 3)
            continue;
        try {
            curve[e].emplace(r_e.mCurve);
        } catch (const std::invalid_argument&) {
            continue;
        }
        const double length = trace[e]->Length();
        double gap = 0.0, scale = 0.0;
        for (int k = 0; k <= 8; ++k) {
            const detail::EdgeVec3 q = trace[e]->Point(trace[e]->ParamAtLength(length * k / 8.0));
            const detail::EdgeVec3 c = ptc_curve_closest(*curve[e], q);
            gap = std::max(gap,
                           std::sqrt((c[0] - q[0]) * (c[0] - q[0]) + (c[1] - q[1]) * (c[1] - q[1]) +
                                     (c[2] - q[2]) * (c[2] - q[2])));
            for (const double v : q)
                scale = std::max(scale, std::fabs(v));
        }
        if (!(gap <= std::max(1e-6 * length, 1e-9 * scale))) {
            warn("edge " + std::to_string(e) + ": its 3-D curve lies up to " + std::to_string(gap) +
                 " from its trim; its points go onto the trim");
            curve[e].reset();
        }
    }

    // Each point, independently.
    enum class Fate : std::uint8_t { Stay, Face, Edge, Corner, Ambiguous, TooFar, Failed };
    std::vector<Fate> fate(np, Fate::Stay);
    std::vector<char> outside(np, 0);
    std::vector<char> on_curve(np, 0);
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
                const auto& r_c = curve[static_cast<std::size_t>(r_o.mEdge)];
                if (r_c) {
                    x = ptc_curve_closest(*r_c, p);
                    on_curve[i] = 1;
                } else {
                    x = r_t->Point(r_t->ParamAtLength(detail::edge_closest_length(*r_t, p)));
                }
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
            if (on_face && !loops[f].empty() && !detail::brep_inside(loops[f], uv[0], uv[1]) &&
                !ptc_on_trim(loops[f], uv[0], uv[1], trim_tol[f]))
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
                result.mNumOnCurves += on_curve[i] ? 1 : 0;
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
