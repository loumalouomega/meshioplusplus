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
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// External includes
#ifdef MESHIOPLUSPLUS_HAS_NETGEN
namespace nglib {
#include <nglib.h>
}  // namespace nglib
#endif

// Project includes
#include "../detail/tri_bvh.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/tetrahedralize.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kTetPrefix = "meshio++: tetrahedralize: ";
constexpr double kTetPi = 3.14159265358979323846;

using TetV3 = std::array<double, 3>;

[[noreturn]] void tet_fail(const std::string& rWhy) {
    throw std::invalid_argument(std::string(kTetPrefix) + rWhy);
}

TetV3 tet_sub(const TetV3& rA, const TetV3& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}

TetV3 tet_cross(const TetV3& rA, const TetV3& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}

double tet_dot(const TetV3& rA, const TetV3& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

// Six times the signed volume of (A, B, C, D): positive when D lies on the
// side of (A, B, C) its right-hand normal points away from, meshio's order.
double tet_volume6(const TetV3& rA, const TetV3& rB, const TetV3& rC, const TetV3& rD) {
    return tet_dot(tet_sub(rB, rA), tet_cross(tet_sub(rC, rA), tet_sub(rD, rA)));
}

// The six dihedral angles of a tetrahedron, in degrees: at each edge, pi
// less the angle between the outward normals of the two faces there.
void tet_dihedrals(const std::array<TetV3, 4>& rX, double& rLo, double& rHi) {
    static constexpr int kFaces[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
    std::array<TetV3, 4> n{};
    for (int f = 0; f < 4; ++f) {
        const TetV3& a = rX[static_cast<std::size_t>(kFaces[f][0])];
        const TetV3& b = rX[static_cast<std::size_t>(kFaces[f][1])];
        const TetV3& c = rX[static_cast<std::size_t>(kFaces[f][2])];
        TetV3 v = tet_cross(tet_sub(b, a), tet_sub(c, a));
        const double l = std::sqrt(tet_dot(v, v));
        if (l > 0.0)
            v = {v[0] / l, v[1] / l, v[2] / l};
        n[static_cast<std::size_t>(f)] = v;
    }
    // Faces f and g (the ones opposite corners f and g) meet at the edge of
    // the other two corners.
    for (int f = 0; f < 4; ++f)
        for (int g = f + 1; g < 4; ++g) {
            const double c =
                std::clamp(-tet_dot(n[static_cast<std::size_t>(f)], n[static_cast<std::size_t>(g)]),
                           -1.0, 1.0);
            const double angle = std::acos(c) * (180.0 / kTetPi);
            rLo = std::min(rLo, angle);
            rHi = std::max(rHi, angle);
        }
}

#ifdef MESHIOPLUSPLUS_HAS_NETGEN
// nglib's state is global: one call at a time, and one Ng_Init.
std::mutex& tet_netgen_mutex() {
    static std::mutex m;
    return m;
}
#endif

}  // namespace

bool tetrahedralize_has_netgen() {
#ifdef MESHIOPLUSPLUS_HAS_NETGEN
    return true;
#else
    return false;
#endif
}

TetrahedralizeResult tetrahedralize(const Mesh& rSurface, const TetrahedralizeOptions& rOptions) {
    const TetrahedralizeOptions& o = rOptions;
    if (o.mBackend != TetrahedralizeBackend::Auto && o.mBackend != TetrahedralizeBackend::Netgen &&
        o.mBackend != TetrahedralizeBackend::Octree)
        tet_fail("backend must be auto, netgen or octree");
    if (!(std::isfinite(o.mMaxSize) && o.mMaxSize >= 0.0))
        tet_fail("max_size must be finite and not negative");
    if (!(std::isfinite(o.mMaxVolume) && o.mMaxVolume >= 0.0))
        tet_fail("max_volume must be finite and not negative");
    if (!(std::isfinite(o.mGrading) && o.mGrading >= 0.0 && o.mGrading <= 1.0))
        tet_fail("grading must lie in [0, 1]");
    if (o.mOptimizeSteps < 0 || o.mOptimizeSteps > 100)
        tet_fail("optimize_steps must be between 0 and 100");
    if (o.mMaxCells < 1)
        tet_fail("max_cells must be positive");

    // The triangles, block by block.
    const std::size_t np = rSurface.NumPoints();
    const NDArray& r_pts = rSurface.Points();
    if (np > 0 && (r_pts.Ndim() != 2 || r_pts.Shape()[1] != 3))
        tet_fail("the surface's points must be three-dimensional");
    std::vector<TetV3> x(np);
    for (std::size_t i = 0; i < np; ++i)
        for (std::size_t k = 0; k < 3; ++k)
            x[i][k] = detail::read_double(r_pts, 3 * i + k);
    std::vector<std::array<std::int64_t, 3>> tris;
    std::vector<std::int64_t> tri_face;  // mFaceArray, or -1
    bool has_face = false;
    for (std::size_t b = 0; b < rSurface.NumCellBlocks(); ++b) {
        const auto cells = rSurface.Cells(b);
        const std::string type(cells.Type());
        if (type.rfind("line", 0) == 0 || type == "vertex")
            continue;
        if (type == "quad" || type == "polygon")
            tet_fail("a '" + type +
                     "' block: the surface must be triangles (split them with convert_cells)");
        if (type != "triangle" || cells.IsRagged())
            tet_fail("a '" + type +
                     "' block: the surface must be linear triangles (linearize or extract_surface "
                     "first)");
        const bool face = !o.mFaceArray.empty() && rSurface.HasCellData(o.mFaceArray);
        has_face = has_face || face;
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            std::array<std::int64_t, 3> t{};
            for (std::size_t k = 0; k < 3; ++k) {
                t[k] = detail::read_int(cells.Conn(), 3 * c + k);
                if (t[k] < 0 || static_cast<std::size_t>(t[k]) >= np)
                    tet_fail("a triangle names a point out of range");
            }
            tris.push_back(t);
            tri_face.push_back(face ? detail::read_int(rSurface.CellData(o.mFaceArray, b), c) : -1);
        }
    }
    if (tris.empty())
        tet_fail("the surface has no triangle");
    {
        Mesh check;
        NDArray p = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
        for (std::size_t i = 0; i < np; ++i)
            for (std::size_t k = 0; k < 3; ++k)
                p.As<double>()[3 * i + k] = x[i][k];
        check.AssignPoints(std::move(p));
        NDArray c = NDArray::Uninit(DType::Int64, {tris.size(), std::size_t{3}});
        for (std::size_t t = 0; t < tris.size(); ++t)
            for (std::size_t k = 0; k < 3; ++k)
                c.As<std::int64_t>()[3 * t + k] = tris[t][k];
        check.AddCellBlock("triangle", std::move(c));
        const SurfaceQuality q = surface_watertight_check(check);
        if (!q.mWatertight)
            tet_fail(
                "the surface must be closed, manifold, consistently oriented and without "
                "degenerate triangles: it has " +
                std::to_string(q.mBoundaryEdges) + " open edges, " +
                std::to_string(q.mNonManifoldEdges) + " non-manifold edges, " +
                std::to_string(q.mInconsistentPairs) + " misoriented edge pairs and " +
                std::to_string(q.mDegenerateTriangles) + " degenerate triangles (repair it first)");
    }

    // The backend: Netgen when asked for or, under Auto, when built in; the
    // native octree route otherwise.
    const bool has_netgen = tetrahedralize_has_netgen();
    if (o.mBackend == TetrahedralizeBackend::Netgen && !has_netgen)
        throw Unsupported(std::string(kTetPrefix) +
                          "the netgen backend is not built in: rebuild with "
                          "-DMESHIOPLUSPLUS_WITH_NETGEN=ON and Netgen's nglib installed, or use "
                          "backend=octree (roadmap 6.7)");
    const bool use_netgen = o.mBackend == TetrahedralizeBackend::Netgen ||
                            (o.mBackend == TetrahedralizeBackend::Auto && has_netgen);
    const std::size_t nt = tris.size();
    // Shells: the triangles connected across edges, in the order of their
    // first triangle.
    std::vector<std::size_t> parent(nt);
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    const auto find = [&](std::size_t A) {
        while (parent[A] != A)
            A = parent[A] = parent[parent[A]];
        return A;
    };
    {
        std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> first;
        for (std::size_t t = 0; t < nt; ++t)
            for (std::size_t k = 0; k < 3; ++k) {
                const auto key = std::minmax(tris[t][k], tris[t][(k + 1) % 3]);
                const auto [it, fresh] = first.emplace(key, t);
                if (!fresh) {
                    const std::size_t a = find(it->second), b = find(t);
                    if (a != b)
                        parent[std::max(a, b)] = std::min(a, b);
                }
            }
    }
    std::vector<std::int64_t> shell_of(nt, -1);
    std::vector<std::vector<std::size_t>> shells;
    {
        std::map<std::size_t, std::size_t> index;
        for (std::size_t t = 0; t < nt; ++t) {
            const auto [it, fresh] = index.emplace(find(t), shells.size());
            if (fresh)
                shells.emplace_back();
            shells[it->second].push_back(t);
            shell_of[t] = static_cast<std::int64_t>(it->second);
        }
    }
    const std::size_t ns = shells.size();
    // Each shell's signed volume; one facing in is turned outward.
    std::vector<double> volume(ns, 0.0);
    std::vector<char> turned(ns, 0);
    for (std::size_t s = 0; s < ns; ++s) {
        for (const std::size_t t : shells[s])
            volume[s] += tet_volume6({0.0, 0.0, 0.0}, x[static_cast<std::size_t>(tris[t][0])],
                                     x[static_cast<std::size_t>(tris[t][1])],
                                     x[static_cast<std::size_t>(tris[t][2])]) /
                         6.0;
        if (volume[s] < 0.0) {
            turned[s] = 1;
            volume[s] = -volume[s];
        }
        if (!(volume[s] > 0.0))
            tet_fail("shell " + std::to_string(s) + " encloses no volume");
    }
    const auto oriented = [&](std::size_t T) {
        std::array<std::int64_t, 3> t = tris[T];
        if (turned[static_cast<std::size_t>(shell_of[T])])
            std::swap(t[1], t[2]);
        return t;
    };
    // Nesting: shell s lies in shell r when the winding number of r at one of
    // s's points is 1; s's region is in its innermost container's.
    std::vector<std::int64_t> container(ns, -1);
    for (std::size_t s = 0; s < ns; ++s) {
        const detail::Vec3 q = x[static_cast<std::size_t>(tris[shells[s][0]][0])];
        for (std::size_t r = 0; r < ns; ++r) {
            if (r == s || !(volume[r] > volume[s]))
                continue;
            double w = 0.0;
            for (const std::size_t t : shells[r]) {
                const auto v = oriented(t);
                w += detail::tri_solid_angle(x[static_cast<std::size_t>(v[0])],
                                             x[static_cast<std::size_t>(v[1])],
                                             x[static_cast<std::size_t>(v[2])], q);
            }
            if (std::fabs(w) / (4.0 * kTetPi) > 0.5 &&
                (container[s] < 0 || volume[r] < volume[static_cast<std::size_t>(container[s])]))
                container[s] = static_cast<std::int64_t>(r);
        }
    }

    double maxh = o.mMaxSize > 0.0 ? o.mMaxSize : 0.0;
    if (o.mMaxVolume > 0.0) {
        const double h = std::cbrt(6.0 * std::sqrt(2.0) * o.mMaxVolume);
        maxh = maxh > 0.0 ? std::min(maxh, h) : h;
    }
    // The region of shell s: its triangles outward, its children's inward.
    const auto region_faces = [&](std::size_t S) {
        std::vector<std::array<std::int64_t, 3>> faces;
        for (const std::size_t t : shells[S])
            faces.push_back(oriented(t));
        for (std::size_t c = 0; c < ns; ++c)
            if (container[c] == static_cast<std::int64_t>(S))
                for (const std::size_t t : shells[c]) {
                    auto v = oriented(t);
                    std::swap(v[1], v[2]);
                    faces.push_back(v);
                }
        return faces;
    };
    TetrahedralizeResult result;
    result.mNumRegions = static_cast<std::int64_t>(ns);
    std::vector<TetV3> points;
    std::vector<std::array<std::int64_t, 4>> tets;
    std::vector<std::int64_t> tet_region;
    // The result's triangles and their face_array values: the surface's own
    // (Netgen), or the tetrahedra's boundary facets (octree).
    std::vector<std::array<std::int64_t, 3>> bnd;
    std::vector<std::int64_t> bnd_face;
    if (use_netgen) {
#ifdef MESHIOPLUSPLUS_HAS_NETGEN
        result.mBackend = "netgen";
        result.mBoundaryPreserved = true;
        points.assign(x.begin(), x.end());
        {
            std::lock_guard<std::mutex> lock(tet_netgen_mutex());
            static std::once_flag s_init;
            std::call_once(s_init, [] { nglib::Ng_Init(); });
            for (std::size_t s = 0; s < ns; ++s) {
                // The region's triangles; the points they use, ascending.
                const std::vector<std::array<std::int64_t, 3>> faces = region_faces(s);
                std::vector<std::int64_t> used;
                for (const auto& r_f : faces)
                    used.insert(used.end(), r_f.begin(), r_f.end());
                std::sort(used.begin(), used.end());
                used.erase(std::unique(used.begin(), used.end()), used.end());
                std::map<std::int64_t, int> local;
                for (std::size_t i = 0; i < used.size(); ++i)
                    local[used[i]] = static_cast<int>(i) + 1;
                nglib::Ng_Mesh* p_m = nglib::Ng_NewMesh();
                for (const std::int64_t g : used) {
                    double p[3] = {x[static_cast<std::size_t>(g)][0],
                                   x[static_cast<std::size_t>(g)][1],
                                   x[static_cast<std::size_t>(g)][2]};
                    nglib::Ng_AddPoint(p_m, p);
                }
                for (const auto& r_f : faces) {
                    int pi[3] = {local.at(r_f[0]), local.at(r_f[1]), local.at(r_f[2])};
                    nglib::Ng_AddSurfaceElement(p_m, nglib::NG_TRIG, pi);
                }
                nglib::Ng_Meshing_Parameters mp;
                if (maxh > 0.0)
                    mp.maxh = maxh;
                mp.grading = o.mGrading;
                mp.optsteps_3d = o.mOptimizeSteps;
                mp.optvolmeshenable = o.mOptimizeSteps > 0 ? 1 : 0;
                mp.optsurfmeshenable = 0;
                const nglib::Ng_Result status = nglib::Ng_GenerateVolumeMesh(p_m, &mp);
                const int ne = nglib::Ng_GetNE(p_m);
                const int nlocal = nglib::Ng_GetNP(p_m);
                if (status != nglib::NG_OK || ne <= 0) {
                    nglib::Ng_DeleteMesh(p_m);
                    throw std::runtime_error(std::string(kTetPrefix) + "Netgen failed on region " +
                                             std::to_string(s) + " (status " +
                                             std::to_string(static_cast<int>(status)) +
                                             "); is the surface free of self-intersections?");
                }
                if (static_cast<std::int64_t>(tets.size()) + ne > o.mMaxCells) {
                    nglib::Ng_DeleteMesh(p_m);
                    throw std::length_error(
                        std::string(kTetPrefix) + "the mesh would have more than " +
                        std::to_string(o.mMaxCells) + " tetrahedra (max_cells)");
                }
                // The surface kept: its points where they were, its triangles as
                // given.
                for (std::size_t i = 0; i < used.size(); ++i) {
                    double p[3];
                    nglib::Ng_GetPoint(p_m, static_cast<int>(i) + 1, p);
                    const TetV3& r_g = x[static_cast<std::size_t>(used[i])];
                    if (p[0] != r_g[0] || p[1] != r_g[1] || p[2] != r_g[2])
                        result.mBoundaryPreserved = false;
                }
                if (nglib::Ng_GetNSE(p_m) != static_cast<int>(faces.size())) {
                    result.mBoundaryPreserved = false;
                } else {
                    std::vector<std::array<std::int64_t, 3>> want, got;
                    for (const auto& r_f : faces) {
                        std::array<std::int64_t, 3> k = r_f;
                        std::sort(k.begin(), k.end());
                        want.push_back(k);
                    }
                    for (int e = 1; e <= static_cast<int>(faces.size()); ++e) {
                        int pi[NG_SURFACE_ELEMENT_MAXPOINTS];
                        nglib::Ng_GetSurfaceElement(p_m, e, pi);
                        std::array<std::int64_t, 3> k{};
                        for (std::size_t j = 0; j < 3; ++j)
                            k[j] = pi[j] >= 1 && pi[j] <= static_cast<int>(used.size())
                                       ? used[static_cast<std::size_t>(pi[j] - 1)]
                                       : -1;
                        std::sort(k.begin(), k.end());
                        got.push_back(k);
                    }
                    std::sort(want.begin(), want.end());
                    std::sort(got.begin(), got.end());
                    if (want != got)
                        result.mBoundaryPreserved = false;
                }
                // New points, after every earlier region's.
                std::vector<std::int64_t> global(static_cast<std::size_t>(nlocal) + 1, -1);
                for (std::size_t i = 0; i < used.size(); ++i)
                    global[i + 1] = used[i];
                for (int i = static_cast<int>(used.size()) + 1; i <= nlocal; ++i) {
                    double p[3];
                    nglib::Ng_GetPoint(p_m, i, p);
                    global[static_cast<std::size_t>(i)] = static_cast<std::int64_t>(points.size());
                    points.push_back({p[0], p[1], p[2]});
                }
                for (int e = 1; e <= ne; ++e) {
                    int pi[NG_VOLUME_ELEMENT_MAXPOINTS];
                    if (nglib::Ng_GetVolumeElement(p_m, e, pi) != nglib::NG_TET) {
                        nglib::Ng_DeleteMesh(p_m);
                        throw std::runtime_error(std::string(kTetPrefix) +
                                                 "Netgen returned an element that is not a "
                                                 "tetrahedron");
                    }
                    std::array<std::int64_t, 4> t{};
                    for (std::size_t j = 0; j < 4; ++j)
                        t[j] = global[static_cast<std::size_t>(pi[j])];
                    if (tet_volume6(points[static_cast<std::size_t>(t[0])],
                                    points[static_cast<std::size_t>(t[1])],
                                    points[static_cast<std::size_t>(t[2])],
                                    points[static_cast<std::size_t>(t[3])]) < 0.0)
                        std::swap(t[1], t[2]);
                    tets.push_back(t);
                    tet_region.push_back(static_cast<std::int64_t>(s));
                }
                nglib::Ng_DeleteMesh(p_m);
            }
        }
        if (!result.mBoundaryPreserved)
            throw std::runtime_error(std::string(kTetPrefix) +
                                     "Netgen moved the surface's points or triangles");
        bnd = tris;
        bnd_face = tri_face;
#endif
    } else {
        // The native route (octree_mesh's isosurface stuffing, Labelle and
        // Shewchuk 2007), region by region from the same faces Netgen takes:
        // the winding number of a region's faces is 1 inside it and 0 in its
        // holes, so the stuffing cuts its outer and its inner surfaces. The
        // boundary is fitted to the surface within the distance reported,
        // not kept, and two regions do not share their interface's points.
        result.mBackend = "octree";
        result.mBoundaryPreserved = false;
        double edge_sum = 0.0;
        for (const auto& r_t : tris)
            for (std::size_t k = 0; k < 3; ++k) {
                const TetV3 e = tet_sub(x[static_cast<std::size_t>(r_t[(k + 1) % 3])],
                                        x[static_cast<std::size_t>(r_t[k])]);
                edge_sum += std::sqrt(tet_dot(e, e));
            }
        const double mean_edge = edge_sum / (3.0 * static_cast<double>(nt));
        const double skin_size = maxh > 0.0 ? std::min(maxh, mean_edge) : mean_edge;
        for (std::size_t s = 0; s < ns; ++s) {
            const std::vector<std::array<std::int64_t, 3>> faces = region_faces(s);
            std::map<std::int64_t, std::int64_t> local;
            std::vector<std::int64_t> used;
            for (const auto& r_f : faces)
                for (const std::int64_t v : r_f)
                    if (local.emplace(v, 0).second)
                        used.push_back(v);
            std::sort(used.begin(), used.end());
            for (std::size_t i = 0; i < used.size(); ++i)
                local[used[i]] = static_cast<std::int64_t>(i);
            Mesh skin;
            NDArray sp = NDArray::Uninit(DType::Float64, {used.size(), std::size_t{3}});
            for (std::size_t i = 0; i < used.size(); ++i)
                for (std::size_t k = 0; k < 3; ++k)
                    sp.As<double>()[3 * i + k] = x[static_cast<std::size_t>(used[i])][k];
            skin.AssignPoints(std::move(sp));
            NDArray sc = NDArray::Uninit(DType::Int64, {faces.size(), std::size_t{3}});
            for (std::size_t f = 0; f < faces.size(); ++f)
                for (std::size_t k = 0; k < 3; ++k)
                    sc.As<std::int64_t>()[3 * f + k] = local.at(faces[f][k]);
            skin.AddCellBlock("triangle", std::move(sc));
            OctreeMeshOptions oo;
            oo.mOutput = OctreeMeshOutput::Tetra;
            oo.mSkinGrading = OctreeMeshSkinGrading::Uniform;
            oo.mSurfaceSize = skin_size;
            oo.mMinSize = skin_size;
            // The largest leaf: max_size, or octree_mesh's default (an eighth
            // of the region's diagonal), never under the skin's size.
            double lo3[3] = {std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity()};
            double hi3[3] = {-lo3[0], -lo3[1], -lo3[2]};
            for (const std::int64_t v : used)
                for (std::size_t k = 0; k < 3; ++k) {
                    lo3[k] = std::min(lo3[k], x[static_cast<std::size_t>(v)][k]);
                    hi3[k] = std::max(hi3[k], x[static_cast<std::size_t>(v)][k]);
                }
            const double diagonal = std::sqrt((hi3[0] - lo3[0]) * (hi3[0] - lo3[0]) +
                                              (hi3[1] - lo3[1]) * (hi3[1] - lo3[1]) +
                                              (hi3[2] - lo3[2]) * (hi3[2] - lo3[2]));
            oo.mMaxSize = std::max(maxh > 0.0 ? maxh : diagonal / 8.0, skin_size);
            oo.mFitIterations = o.mOptimizeSteps;
            oo.mMaxTets =
                std::max<std::int64_t>(1, o.mMaxCells - static_cast<std::int64_t>(tets.size()));
            oo.mAttachDistance = false;
            OctreeMeshResult om;
            try {
                om = octree_mesh(skin, oo);
            } catch (const std::invalid_argument& rE) {
                // The lattice's own budget, which max_cells sets.
                if (std::string(rE.what()).find("max_tets") == std::string::npos)
                    throw;
                throw std::length_error(std::string(kTetPrefix) + "the mesh would have more than " +
                                        std::to_string(o.mMaxCells) + " tetrahedra (max_cells)");
            }
            result.mHausdorffToSurface = std::max(result.mHausdorffToSurface, om.mHausdorffToSkin);
            const std::size_t base = points.size();
            const NDArray& r_op = om.mMesh.Points();
            for (std::size_t i = 0; i < om.mMesh.NumPoints(); ++i)
                points.push_back({detail::read_double(r_op, 3 * i),
                                  detail::read_double(r_op, 3 * i + 1),
                                  detail::read_double(r_op, 3 * i + 2)});
            for (std::size_t b = 0; b < om.mMesh.NumCellBlocks(); ++b) {
                const auto cells = om.mMesh.Cells(b);
                if (std::string(cells.Type()) != "tetra")
                    continue;
                for (std::size_t c = 0; c < cells.NumCells(); ++c) {
                    std::array<std::int64_t, 4> t{};
                    for (std::size_t k = 0; k < 4; ++k)
                        t[k] = static_cast<std::int64_t>(base) +
                               detail::read_int(cells.Conn(), 4 * c + k);
                    if (tet_volume6(points[static_cast<std::size_t>(t[0])],
                                    points[static_cast<std::size_t>(t[1])],
                                    points[static_cast<std::size_t>(t[2])],
                                    points[static_cast<std::size_t>(t[3])]) < 0.0)
                        std::swap(t[1], t[2]);
                    tets.push_back(t);
                    tet_region.push_back(static_cast<std::int64_t>(s));
                }
            }
            if (static_cast<std::int64_t>(tets.size()) > o.mMaxCells)
                throw std::length_error(std::string(kTetPrefix) + "the mesh would have more than " +
                                        std::to_string(o.mMaxCells) + " tetrahedra (max_cells)");
        }
        // The boundary: the facets one tetrahedron has, outward, in the order
        // of their tetrahedra.
        static constexpr int kSide[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
        std::map<std::array<std::int64_t, 3>, std::int64_t> uses;
        for (const auto& r_t : tets)
            for (const auto& r_k : kSide) {
                std::array<std::int64_t, 3> k{r_t[static_cast<std::size_t>(r_k[0])],
                                              r_t[static_cast<std::size_t>(r_k[1])],
                                              r_t[static_cast<std::size_t>(r_k[2])]};
                std::sort(k.begin(), k.end());
                ++uses[k];
            }
        for (const auto& r_t : tets)
            for (int f = 0; f < 4; ++f) {
                std::array<std::int64_t, 3> v{r_t[static_cast<std::size_t>(kSide[f][0])],
                                              r_t[static_cast<std::size_t>(kSide[f][1])],
                                              r_t[static_cast<std::size_t>(kSide[f][2])]};
                std::array<std::int64_t, 3> k = v;
                std::sort(k.begin(), k.end());
                if (uses.at(k) != 1)
                    continue;
                // Facing away from the tetrahedron's fourth corner.
                const TetV3& r_a = points[static_cast<std::size_t>(v[0])];
                const TetV3 n = tet_cross(tet_sub(points[static_cast<std::size_t>(v[1])], r_a),
                                          tet_sub(points[static_cast<std::size_t>(v[2])], r_a));
                if (tet_dot(
                        n,
                        tet_sub(points[static_cast<std::size_t>(r_t[static_cast<std::size_t>(f)])],
                                r_a)) > 0.0)
                    std::swap(v[1], v[2]);
                bnd.push_back(v);
            }
        // Each boundary triangle takes the face_array value of the surface
        // triangle nearest its centroid (ties to the lower index).
        if (has_face && !bnd.empty()) {
            Mesh surface;
            NDArray p = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
            for (std::size_t i = 0; i < np; ++i)
                for (std::size_t k = 0; k < 3; ++k)
                    p.As<double>()[3 * i + k] = x[i][k];
            surface.AssignPoints(std::move(p));
            NDArray c = NDArray::Uninit(DType::Int64, {nt, std::size_t{3}});
            for (std::size_t t = 0; t < nt; ++t)
                for (std::size_t k = 0; k < 3; ++k)
                    c.As<std::int64_t>()[3 * t + k] = tris[t][k];
            surface.AddCellBlock("triangle", std::move(c));
            const detail::TriangleSoup soup = detail::build_triangle_soup(surface, "");
            const detail::DistanceQuery query =
                detail::build_distance_query(soup, SurfaceDistanceOptions{});
            std::vector<detail::Vec3> centres(bnd.size());
            for (std::size_t t = 0; t < bnd.size(); ++t)
                for (std::size_t k = 0; k < 3; ++k)
                    centres[t][k] = (points[static_cast<std::size_t>(bnd[t][0])][k] +
                                     points[static_cast<std::size_t>(bnd[t][1])][k] +
                                     points[static_cast<std::size_t>(bnd[t][2])][k]) /
                                    3.0;
            const auto hits = detail::query_closest_points(query, centres);
            bnd_face.resize(bnd.size(), -1);
            for (std::size_t t = 0; t < bnd.size(); ++t)
                if (hits[t].mFound && hits[t].mSourceCell >= 0)
                    bnd_face[t] = tri_face[static_cast<std::size_t>(hits[t].mSourceCell)];
        } else {
            bnd_face.assign(bnd.size(), -1);
        }
        if (o.mBackend == TetrahedralizeBackend::Auto)
            result.mWarnings.push_back(
                "no Netgen backend is built in: the native octree backend fits the boundary "
                "within " +
                std::to_string(result.mHausdorffToSurface) +
                " of the surface instead of keeping its triangles (build with "
                "-DMESHIOPLUSPLUS_WITH_NETGEN=ON to keep them)");
        if (ns > 1)
            result.mWarnings.push_back(
                "the octree backend meshes each of the " + std::to_string(ns) +
                " regions on its own: they do not share points across their interfaces");
        for (const std::string& r_w : result.mWarnings)
            log::warn("{}{}", kTetPrefix, r_w);
    }
    result.mNumBoundaryTriangles = static_cast<std::int64_t>(bnd.size());

    // The measures.
    result.mNumTetrahedra = static_cast<std::int64_t>(tets.size());
    result.mNumPoints = static_cast<std::int64_t>(points.size());
    result.mNumNewPoints =
        static_cast<std::int64_t>(points.size()) - (use_netgen ? static_cast<std::int64_t>(np) : 0);
    double lo = 180.0, hi = 0.0;
    for (const auto& r_t : tets) {
        const std::array<TetV3, 4> c{
            points[static_cast<std::size_t>(r_t[0])], points[static_cast<std::size_t>(r_t[1])],
            points[static_cast<std::size_t>(r_t[2])], points[static_cast<std::size_t>(r_t[3])]};
        result.mVolume += tet_volume6(c[0], c[1], c[2], c[3]) / 6.0;
        tet_dihedrals(c, lo, hi);
    }
    result.mMinDihedralDeg = tets.empty() ? 0.0 : lo;
    result.mMaxDihedralDeg = tets.empty() ? 0.0 : hi;

    // The mesh: the points, the tetrahedra, the surface's triangles.
    Mesh& r_m = result.mMesh;
    NDArray pts = NDArray::Uninit(DType::Float64, {points.size(), std::size_t{3}});
    for (std::size_t i = 0; i < points.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            pts.As<double>()[3 * i + k] = points[i][k];
    r_m.AssignPoints(std::move(pts));
    NDArray tc = NDArray::Uninit(DType::Int64, {tets.size(), std::size_t{4}});
    for (std::size_t t = 0; t < tets.size(); ++t)
        for (std::size_t k = 0; k < 4; ++k)
            tc.As<std::int64_t>()[4 * t + k] = tets[t][k];
    r_m.AddCellBlock("tetra", std::move(tc));
    const std::size_t nb = bnd.size();
    NDArray fc = NDArray::Uninit(DType::Int64, {nb, std::size_t{3}});
    for (std::size_t t = 0; t < nb; ++t)
        for (std::size_t k = 0; k < 3; ++k)
            fc.As<std::int64_t>()[3 * t + k] = bnd[t][k];
    r_m.AddCellBlock("triangle", std::move(fc));
    const auto column = [](std::size_t N, const std::int64_t* pFrom, std::int64_t Fill) {
        NDArray a = NDArray::Uninit(DType::Int64, {N});
        for (std::size_t i = 0; i < N; ++i)
            a.As<std::int64_t>()[i] = pFrom ? pFrom[i] : Fill;
        return a;
    };
    {
        std::vector<NDArray> blocks;
        blocks.push_back(column(tets.size(), tet_region.data(), 0));
        blocks.push_back(column(nb, nullptr, -1));
        r_m.AddCellData(kTetrahedralizeRegion, std::move(blocks));
    }
    if (has_face) {
        std::vector<NDArray> blocks;
        blocks.push_back(column(tets.size(), nullptr, -1));
        blocks.push_back(column(nb, bnd_face.data(), 0));
        r_m.AddCellData(o.mFaceArray, std::move(blocks));
    }
    for (std::size_t s = 0; s < ns; ++s) {
        std::vector<std::int64_t> ids;
        for (std::size_t t = 0; t < tets.size(); ++t)
            if (tet_region[t] == static_cast<std::int64_t>(s))
                ids.push_back(static_cast<std::int64_t>(t));
        r_m.AddRegion(Region("region:" + std::to_string(s), RegionKind::Cell, 3, -1,
                             column(ids.size(), ids.data(), 0)));
    }
    if (has_face) {
        // Each surface triangle's tetrahedron facet, the one outside it
        // first found (a shell between two regions is a facet of both; the
        // lower region's is taken).
        static constexpr int kFacet[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
        std::map<std::array<std::int64_t, 3>, std::pair<std::int64_t, std::int64_t>> facet_of;
        for (std::size_t t = 0; t < tets.size(); ++t)
            for (int f = 0; f < 4; ++f) {
                std::array<std::int64_t, 3> k{tets[t][static_cast<std::size_t>(kFacet[f][0])],
                                              tets[t][static_cast<std::size_t>(kFacet[f][1])],
                                              tets[t][static_cast<std::size_t>(kFacet[f][2])]};
                std::sort(k.begin(), k.end());
                facet_of.emplace(
                    k, std::make_pair(static_cast<std::int64_t>(t), static_cast<std::int64_t>(f)));
            }
        std::map<std::int64_t, std::vector<std::int64_t>> sides;
        for (std::size_t t = 0; t < nb; ++t) {
            if (bnd_face[t] < 0)
                continue;
            std::array<std::int64_t, 3> k = bnd[t];
            std::sort(k.begin(), k.end());
            const auto it = facet_of.find(k);
            if (it == facet_of.end())
                continue;
            auto& r_s = sides[bnd_face[t]];
            r_s.push_back(it->second.first);
            r_s.push_back(it->second.second);
        }
        for (const auto& [face, r_e] : sides) {
            NDArray e = NDArray::Uninit(DType::Int64, {r_e.size() / 2, std::size_t{2}});
            std::copy(r_e.begin(), r_e.end(), e.As<std::int64_t>());
            r_m.AddRegion(
                Region("face:" + std::to_string(face), RegionKind::Side, 2, -1, std::move(e)));
        }
    }
    log::debug("tetrahedralize: {} tetrahedra in {} regions, dihedral angles {:.1f} to {:.1f}",
               result.mNumTetrahedra, result.mNumRegions, result.mMinDihedralDeg,
               result.mMaxDihedralDeg);
    return result;
}

}  // namespace meshioplusplus
