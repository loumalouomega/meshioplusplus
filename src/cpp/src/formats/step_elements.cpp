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
// Fitted elements to STEP (roadmap 7.6.5). See formats/step.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

// Project includes (private, not installed)
#include "../detail/closed_facets.hpp"
#include "../detail/lsq_fit.hpp"
#include "../detail/predicates.hpp"
#include "../detail/step_brep.hpp"
#include "../detail/step_stream.hpp"

namespace meshioplusplus {

namespace {

using SeV3 = std::array<double, 3>;

constexpr const char* kSePrefix = "meshio++: step: write_step_elements: ";
constexpr double kSePi = 3.14159265358979323846;

SeV3 se_add(const SeV3& rA, const SeV3& rB, double S) {
    return {rA[0] + S * rB[0], rA[1] + S * rB[1], rA[2] + S * rB[2]};
}

double se_dot(const SeV3& rA, const SeV3& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

SeV3 se_sub(const SeV3& rA, const SeV3& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}

// A bounded plane: its corners, counter-clockwise about the normal.
struct SePlane {
    SeV3 mPoint{}, mNormal{}, mE1{};
    std::vector<SeV3> mCorners;
};

// The closed tessellated solids, labelled by the rows of their fit table.
struct SeSolids {
    std::vector<double> mXyz;
    std::vector<std::int64_t> mTri;
    std::vector<std::int64_t> mLabel;
    std::vector<std::int64_t> mType;
    std::vector<double> mParams;

    std::int64_t Row(std::int64_t Type, const std::array<double, 10>& rParams) {
        mType.push_back(Type);
        mParams.insert(mParams.end(), rParams.begin(), rParams.end());
        return static_cast<std::int64_t>(mType.size()) - 1;
    }
    std::int64_t Point(const SeV3& rP) {
        mXyz.insert(mXyz.end(), rP.begin(), rP.end());
        return static_cast<std::int64_t>(mXyz.size() / 3) - 1;
    }
    void Tri(std::int64_t A, std::int64_t B, std::int64_t C, std::int64_t Label) {
        mTri.insert(mTri.end(), {A, B, C});
        mLabel.push_back(Label);
    }
    // A flat lid at `rCentre` with outward normal `rNormal` over a ring.
    void Lid(const SeV3& rCentre, const SeV3& rNormal, const std::vector<std::int64_t>& rRing,
             bool Reverse) {
        const std::int64_t row = Row(1, {rCentre[0], rCentre[1], rCentre[2], rNormal[0], rNormal[1],
                                         rNormal[2], 0, 0, 0, 0});
        const std::int64_t c = Point(rCentre);
        const std::size_t n = rRing.size();
        for (std::size_t i = 0; i < n; ++i) {
            const std::int64_t a = rRing[i], b = rRing[(i + 1) % n];
            if (Reverse)
                Tri(c, b, a, row);
            else
                Tri(c, a, b, row);
        }
    }
};

// A ring of N points at height T along the axis, radius R.
std::vector<std::int64_t> se_ring(SeSolids& rS, const SeV3& rP, const SeV3& rA, const SeV3& rU,
                                  const SeV3& rV, double T, double R, std::int64_t N) {
    std::vector<std::int64_t> ring;
    for (std::int64_t i = 0; i < N; ++i) {
        const double th = 2.0 * kSePi * static_cast<double>(i) / static_cast<double>(N);
        SeV3 q = se_add(rP, rA, T);
        q = se_add(q, rU, R * std::cos(th));
        q = se_add(q, rV, R * std::sin(th));
        ring.push_back(rS.Point(q));
    }
    return ring;
}

// The side between two rings (the second further along the axis), outward.
void se_side(SeSolids& rS, const std::vector<std::int64_t>& rLo,
             const std::vector<std::int64_t>& rHi, std::int64_t Label) {
    const std::size_t n = rHi.size();
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        if (rLo.size() == 1) {  // a cone's apex
            rS.Tri(rLo[0], rHi[j], rHi[i], Label);
            continue;
        }
        rS.Tri(rLo[i], rLo[j], rHi[j], Label);
        rS.Tri(rLo[i], rHi[j], rHi[i], Label);
    }
}

// The convex hull of 2-D points, counter-clockwise (Andrew's monotone chain).
std::vector<std::array<double, 2>> se_hull(std::vector<std::array<double, 2>> Pts) {
    std::sort(Pts.begin(), Pts.end());
    Pts.erase(std::unique(Pts.begin(), Pts.end()), Pts.end());
    if (Pts.size() < 3)
        return {};
    std::vector<std::array<double, 2>> h(2 * Pts.size());
    std::size_t k = 0;
    for (std::size_t i = 0; i < Pts.size(); ++i) {
        while (k >= 2 && detail::orient2d(h[k - 2].data(), h[k - 1].data(), Pts[i].data()) <= 0.0)
            --k;
        h[k++] = Pts[i];
    }
    for (std::size_t i = Pts.size() - 1, t = k + 1; i > 0; --i) {
        while (k >= t &&
               detail::orient2d(h[k - 2].data(), h[k - 1].data(), Pts[i - 1].data()) <= 0.0)
            --k;
        h[k++] = Pts[i - 1];
    }
    h.resize(k - 1);
    return h.size() >= 3 ? h : std::vector<std::array<double, 2>>{};
}

// Each bounded plane as a SHELL_BASED_SURFACE_MODEL of one ADVANCED_FACE.
std::vector<std::int64_t> se_write_planes(detail::StepStream& rSt,
                                          const std::vector<SePlane>& rPlanes) {
    using detail::StepStream;
    std::vector<std::int64_t> models;
    auto point = [&](const SeV3& rP) {
        return rSt.Ent("CARTESIAN_POINT(''," + StepStream::Triple(rP.data()) + ")");
    };
    auto direction = [&](const SeV3& rD) {
        return rSt.Ent("DIRECTION(''," + StepStream::Triple(rD.data()) + ")");
    };
    for (std::size_t k = 0; k < rPlanes.size(); ++k) {
        const SePlane& r_pl = rPlanes[k];
        const std::size_t n = r_pl.mCorners.size();
        std::vector<std::int64_t> verts;
        for (const SeV3& r_c : r_pl.mCorners)
            verts.push_back(rSt.Ent("VERTEX_POINT(''," + StepStream::Ref(point(r_c)) + ")"));
        std::string edges;
        for (std::size_t i = 0; i < n; ++i) {
            const SeV3& a = r_pl.mCorners[i];
            const SeV3& b = r_pl.mCorners[(i + 1) % n];
            SeV3 d = se_sub(b, a);
            const double len = std::sqrt(se_dot(d, d));
            for (double& r_x : d)
                r_x /= len;
            const std::int64_t vec = rSt.Ent("VECTOR(''," + StepStream::Ref(direction(d)) + "," +
                                             detail::step_real(len) + ")");
            const std::int64_t line =
                rSt.Ent("LINE(''," + StepStream::Ref(point(a)) + "," + StepStream::Ref(vec) + ")");
            const std::int64_t ec = rSt.Ent("EDGE_CURVE(''," + StepStream::Ref(verts[i]) + "," +
                                            StepStream::Ref(verts[(i + 1) % n]) + "," +
                                            StepStream::Ref(line) + ",.T.)");
            edges +=
                (edges.empty() ? "" : ",") +
                StepStream::Ref(rSt.Ent("ORIENTED_EDGE('',*,*," + StepStream::Ref(ec) + ",.T.)"));
        }
        const std::int64_t loop = rSt.Ent("EDGE_LOOP('',(" + edges + "))");
        const std::int64_t bound =
            rSt.Ent("FACE_OUTER_BOUND(''," + StepStream::Ref(loop) + ",.T.)");
        const std::int64_t place =
            rSt.Ent("AXIS2_PLACEMENT_3D(''," + StepStream::Ref(point(r_pl.mPoint)) + "," +
                    StepStream::Ref(direction(r_pl.mNormal)) + "," +
                    StepStream::Ref(direction(r_pl.mE1)) + ")");
        const std::int64_t plane = rSt.Ent("PLANE(''," + StepStream::Ref(place) + ")");
        const std::int64_t face =
            rSt.Ent("ADVANCED_FACE('plane" + std::to_string(k + 1) + "',(" +
                    StepStream::Ref(bound) + ")," + StepStream::Ref(plane) + ",.T.)");
        const std::int64_t shell = rSt.Ent("OPEN_SHELL('',(" + StepStream::Ref(face) + "))");
        models.push_back(rSt.Ent("SHELL_BASED_SURFACE_MODEL('',(" + StepStream::Ref(shell) + "))"));
    }
    return models;
}

// Each ball as one SPHERICAL_SURFACE face bounded by a vertex loop at its
// south pole (the form OpenCASCADE writes): exact, with no seam.
std::vector<std::int64_t> se_write_balls(detail::StepStream& rSt,
                                         const std::vector<std::array<double, 4>>& rBalls) {
    using detail::StepStream;
    std::vector<std::int64_t> out;
    for (std::size_t k = 0; k < rBalls.size(); ++k) {
        const auto& r_b = rBalls[k];
        const double pole[3] = {r_b[0], r_b[1], r_b[2] - r_b[3]};
        const double centre[3] = {r_b[0], r_b[1], r_b[2]};
        const double z[3] = {0.0, 0.0, 1.0}, x[3] = {1.0, 0.0, 0.0};
        const std::int64_t vp = rSt.Ent(
            "VERTEX_POINT(''," +
            StepStream::Ref(rSt.Ent("CARTESIAN_POINT(''," + StepStream::Triple(pole) + ")")) + ")");
        const std::int64_t loop = rSt.Ent("VERTEX_LOOP(''," + StepStream::Ref(vp) + ")");
        const std::int64_t bound = rSt.Ent("FACE_BOUND(''," + StepStream::Ref(loop) + ",.T.)");
        const std::int64_t place = rSt.Ent(
            "AXIS2_PLACEMENT_3D(''," +
            StepStream::Ref(rSt.Ent("CARTESIAN_POINT(''," + StepStream::Triple(centre) + ")")) +
            "," + StepStream::Ref(rSt.Ent("DIRECTION(''," + StepStream::Triple(z) + ")")) + "," +
            StepStream::Ref(rSt.Ent("DIRECTION(''," + StepStream::Triple(x) + ")")) + ")");
        const std::int64_t sphere = rSt.Ent("SPHERICAL_SURFACE(''," + StepStream::Ref(place) + "," +
                                            detail::step_real(r_b[3]) + ")");
        const std::int64_t face = rSt.Ent("ADVANCED_FACE('',(" + StepStream::Ref(bound) + ")," +
                                          StepStream::Ref(sphere) + ",.T.)");
        const std::int64_t shell = rSt.Ent("CLOSED_SHELL('',(" + StepStream::Ref(face) + "))");
        out.push_back(rSt.Ent("MANIFOLD_SOLID_BREP('ball" + std::to_string(k + 1) + "'," +
                              StepStream::Ref(shell) + ")"));
    }
    return out;
}

NDArray se_int64(const std::vector<std::int64_t>& rV, std::vector<std::size_t> Shape) {
    NDArray a = NDArray::Uninit(DType::Int64, std::move(Shape));
    std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
    return a;
}

NDArray se_doubles(const std::vector<double>& rV, std::vector<std::size_t> Shape) {
    NDArray a = NDArray::Uninit(DType::Float64, std::move(Shape));
    std::copy(rV.begin(), rV.end(), a.As<double>());
    return a;
}

}  // namespace

StepElementsResult write_step_elements(const std::string& rPath, const Mesh& rFitted,
                                       const StepElementsOptions& rOptions) {
    if (rOptions.mSegments < 8)
        throw WriteError(std::string(kSePrefix) + "segments must be at least 8");
    if (!rFitted.HasCellData("segment:patch") || !rFitted.HasFieldData("fit:patch_id") ||
        !rFitted.HasFieldData("fit:patch_type") || !rFitted.HasFieldData("fit:patch_params"))
        throw WriteError(std::string(kSePrefix) +
                         "the mesh carries no fitted patches (segment, then fit_primitives)");
    const NDArray& ids = rFitted.FieldData("fit:patch_id");
    const NDArray& types = rFitted.FieldData("fit:patch_type");
    const NDArray& params = rFitted.FieldData("fit:patch_params");
    const std::size_t m = ids.Size();
    if (types.Size() != m || params.Size() != 10 * m)
        throw WriteError(std::string(kSePrefix) + "the fit:patch_* arrays disagree");

    // The points of each patch.
    std::map<std::int64_t, std::set<std::int64_t>> points_of;
    for (std::size_t b = 0; b < rFitted.NumCellBlocks(); ++b) {
        const auto cb = rFitted.Cells(b);
        const NDArray& labels = rFitted.CellData("segment:patch", b);
        for (std::size_t c = 0; c < cb.NumCells() && c < labels.Size(); ++c) {
            const std::int64_t l = detail::read_int(labels, c);
            const std::size_t k = cb.IsRagged() ? cb.RowSize(c) : cb.NodesPerCell();
            auto& r_set = points_of[l];
            for (std::size_t i = 0; i < k; ++i)
                r_set.insert(cb.IsRagged()
                                 ? cb.Row(c)[i]
                                 : detail::read_int(cb.Conn(), c * cb.NodesPerCell() + i));
        }
    }
    const NDArray& pts = rFitted.Points();
    const std::size_t dim = rFitted.PointDim();
    auto point = [&](std::int64_t V) {
        SeV3 p{};
        for (std::size_t k = 0; k < dim && k < 3; ++k)
            p[k] = detail::read_double(pts, static_cast<std::size_t>(V) * dim + k);
        return p;
    };

    std::vector<std::size_t> rows;
    if (rOptions.mPatches.empty()) {
        for (std::size_t r = 0; r < m; ++r)
            rows.push_back(r);
    } else {
        for (const std::int64_t want : rOptions.mPatches) {
            std::size_t r = 0;
            while (r < m && detail::read_int(ids, r) != want)
                ++r;
            if (r == m)
                throw WriteError(std::string(kSePrefix) + "no fitted patch " +
                                 std::to_string(want));
            rows.push_back(r);
        }
    }

    StepElementsResult result;
    SeSolids solids;
    std::vector<SePlane> planes;
    std::vector<std::array<double, 4>> balls;  // centre, radius
    const std::int64_t n_seg = rOptions.mSegments;
    double scale = 0.0;
    for (const std::size_t r : rows) {
        const std::int64_t label = detail::read_int(ids, r);
        const std::int64_t type = detail::read_int(types, r);
        std::array<double, 10> p{};
        bool finite = true;
        for (std::size_t q = 0; q < 10; ++q) {
            p[q] = detail::read_double(params, 10 * r + q);
            finite = finite && std::isfinite(p[q]);
        }
        const auto it = points_of.find(label);
        if (type < 1 || type > 5 || !finite || it == points_of.end() || it->second.size() < 3) {
            ++result.mNumSkipped;
            if (!rOptions.mPatches.empty())
                result.mWarnings.push_back("patch " + std::to_string(label) +
                                           " has no fitted primitive to write");
            continue;
        }
        const SeV3 o{p[0], p[1], p[2]};
        SeV3 a{p[3], p[4], p[5]};
        const double an = std::sqrt(se_dot(a, a));
        if (type == 2) {
            a = {0.0, 0.0, 1.0};  // a sphere has no axis
        } else if (!(an > 0.0)) {
            ++result.mNumSkipped;
            continue;
        } else {
            for (double& r_x : a)
                r_x /= an;
        }
        SeV3 u{}, v{};
        detail::lsq_frame(a.data(), u.data(), v.data());
        double t_lo = 0.0, t_hi = 0.0;
        bool first = true;
        for (const std::int64_t q : it->second) {
            const SeV3 x = point(q);
            const double t = se_dot(se_sub(x, o), a);
            t_lo = first ? t : std::min(t_lo, t);
            t_hi = first ? t : std::max(t_hi, t);
            first = false;
            scale = std::max(scale, std::sqrt(se_dot(se_sub(x, o), se_sub(x, o))));
        }
        if (type == 1) {
            SePlane pl;
            pl.mPoint = o;
            pl.mNormal = a;
            pl.mE1 = u;
            std::vector<std::array<double, 2>> uv;
            for (const std::int64_t q : it->second) {
                const SeV3 d = se_sub(point(q), o);
                uv.push_back({se_dot(d, u), se_dot(d, v)});
            }
            for (const auto& r_h : se_hull(uv))
                pl.mCorners.push_back(se_add(se_add(o, u, r_h[0]), v, r_h[1]));
            if (pl.mCorners.size() < 3) {
                ++result.mNumSkipped;
                continue;
            }
            planes.push_back(std::move(pl));
            ++result.mNumPlanes;
        } else if (type == 2) {
            if (!(p[6] > 0.0)) {
                ++result.mNumSkipped;
                continue;
            }
            balls.push_back({p[0], p[1], p[2], p[6]});
            ++result.mNumSpheres;
        } else if (type == 3 || type == 4) {
            const double tan_a = type == 4 ? std::tan(p[8] * kSePi / 180.0) : 0.0;
            if (type == 4)
                t_lo = std::max(t_lo, 0.0);
            if (!(t_hi > t_lo) || (type == 4 && !(tan_a > 0.0))) {
                ++result.mNumSkipped;
                continue;
            }
            auto radius = [&](double T) { return type == 3 ? p[6] : T * tan_a; };
            const std::int64_t row = solids.Row(type, p);
            const bool apex = type == 4 && t_lo <= 1e-9 * t_hi;
            const std::vector<std::int64_t> lo =
                apex ? std::vector<std::int64_t>{solids.Point(o)}
                     : se_ring(solids, o, a, u, v, t_lo, radius(t_lo), n_seg);
            const std::vector<std::int64_t> hi =
                se_ring(solids, o, a, u, v, t_hi, radius(t_hi), n_seg);
            se_side(solids, lo, hi, row);
            if (!apex)
                solids.Lid(se_add(o, a, t_lo), {-a[0], -a[1], -a[2]}, lo, true);
            solids.Lid(se_add(o, a, t_hi), a, hi, false);
            if (type == 3)
                ++result.mNumCylinders;
            else
                ++result.mNumCones;
        } else {
            const double big = p[6], small = p[7];
            if (!(small > 0.0) || !(big > small)) {
                ++result.mNumSkipped;
                result.mWarnings.push_back("patch " + std::to_string(label) +
                                           ": a torus whose tube reaches its axis is not written");
                continue;
            }
            const std::int64_t row = solids.Row(5, p);
            const std::int64_t n_minor = std::max<std::int64_t>(8, n_seg / 2);
            std::vector<std::int64_t> grid;
            for (std::int64_t i = 0; i < n_seg; ++i) {
                const double phi =
                    2.0 * kSePi * static_cast<double>(i) / static_cast<double>(n_seg);
                SeV3 d = se_add(se_add(SeV3{0, 0, 0}, u, std::cos(phi)), v, std::sin(phi));
                for (std::int64_t j = 0; j < n_minor; ++j) {
                    const double psi =
                        2.0 * kSePi * static_cast<double>(j) / static_cast<double>(n_minor);
                    SeV3 q = se_add(o, d, big + small * std::cos(psi));
                    q = se_add(q, a, small * std::sin(psi));
                    grid.push_back(solids.Point(q));
                }
            }
            auto at = [&](std::int64_t I, std::int64_t J) {
                return grid[static_cast<std::size_t>((I % n_seg) * n_minor + (J % n_minor))];
            };
            for (std::int64_t i = 0; i < n_seg; ++i)
                for (std::int64_t j = 0; j < n_minor; ++j) {
                    solids.Tri(at(i, j), at(i + 1, j), at(i + 1, j + 1), row);
                    solids.Tri(at(i, j), at(i + 1, j + 1), at(i, j + 1), row);
                }
            ++result.mNumTori;
        }
    }
    if (solids.mTri.empty() && planes.empty() && balls.empty())
        throw WriteError(std::string(kSePrefix) + "no fitted primitive to write");

    const double unit = detail::mesh_length_unit(rFitted);
    auto extra = [&](detail::StepStream& rSt, std::vector<std::int64_t>& rSolids,
                     std::vector<std::int64_t>& rSurfaces) {
        const std::vector<std::int64_t> b = se_write_balls(rSt, balls);
        rSolids.insert(rSolids.end(), b.begin(), b.end());
        rSurfaces = se_write_planes(rSt, planes);
    };
    if (solids.mTri.empty()) {
        detail::StepStream st(rPath, "meshio++ elements", unit);
        std::vector<std::int64_t> solid_ids, surface_ids;
        extra(st, solid_ids, surface_ids);
        st.FinishMixed(solid_ids, surface_ids, {});
        return result;
    }
    Mesh skin;
    skin.AssignPoints(se_doubles(solids.mXyz, {solids.mXyz.size() / 3, std::size_t{3}}));
    skin.AddCellBlock("triangle", se_int64(solids.mTri, {solids.mTri.size() / 3, std::size_t{3}}));
    std::vector<NDArray> labels;
    labels.push_back(se_int64(solids.mLabel, {solids.mLabel.size()}));
    skin.AddCellData("segment:patch", std::move(labels));
    const std::size_t nrow = solids.mType.size();
    std::vector<std::int64_t> row_ids(nrow);
    for (std::size_t i = 0; i < nrow; ++i)
        row_ids[i] = static_cast<std::int64_t>(i);
    skin.AddFieldData("fit:patch_id", se_int64(row_ids, {nrow}));
    skin.AddFieldData("fit:patch_type", se_int64(solids.mType, {nrow}));
    skin.AddFieldData("fit:patch_params", se_doubles(solids.mParams, {nrow, std::size_t{10}}));
    skin.AddFieldData("fit:patch_tolerance",
                      se_doubles(std::vector<double>(nrow, 1e-9 * std::max(scale, 1.0)), {nrow}));
    const detail::ClosedFacets shell = detail::closed_facets(skin, kSePrefix);
    detail::StepBrepOptions brep;
    brep.mFitted = true;
    brep.mLengthUnit = unit;
    // The tessellation's chords sag off a curved rim by about R (pi/N)^2 / 2:
    // the rims are circles within that.
    double max_r = 0.0;
    for (std::size_t i = 0; i < nrow; ++i)
        if (solids.mType[i] != 1)
            max_r = std::max(max_r, solids.mParams[10 * i + 6]);
    brep.mTolerance =
        std::max(1e-9 * std::max(scale, 1.0), max_r * (kSePi / static_cast<double>(n_seg)) *
                                                  (kSePi / static_cast<double>(n_seg)));
    if (!planes.empty() || !balls.empty())
        brep.mExtraItems = extra;
    const detail::StepBrepResult r = detail::step_brep_write(rPath, skin, shell, brep);
    if (r.mCounts.mNumDemoted > 0)
        result.mWarnings.push_back(std::to_string(r.mCounts.mNumDemoted) +
                                   " element face(s) written as facets");
    return result;
}

}  // namespace meshioplusplus
