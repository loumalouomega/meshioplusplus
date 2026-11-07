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
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/hull.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "../detail/cdt2d.hpp"
#include "../detail/geometry_utils.hpp"
#include "../detail/predicates.hpp"

namespace meshioplusplus {

namespace {

using HullVec = std::array<double, 3>;

/// Every point as three doubles (a 2-D point gets z = 0); refuses non-finite
/// coordinates and 1-D points by name.
std::vector<HullVec> hull_points(const Mesh& rMesh, const char* pOp) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    if (n > 0 && dim < 2)
        throw std::invalid_argument(std::string("meshio++: ") + pOp +
                                    ": the points must have 2 or 3 coordinates, not " +
                                    std::to_string(dim));
    std::vector<HullVec> out(n, HullVec{0.0, 0.0, 0.0});
    const NDArray& points = rMesh.Points();
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < dim && k < 3; ++k) {
            const double v = detail::read_double(points, i * dim + k);
            if (!std::isfinite(v))
                throw std::invalid_argument(std::string("meshio++: ") + pOp + ": point " +
                                            std::to_string(i) + " has a non-finite coordinate");
            out[i][k] = v;
        }
    return out;
}

/// Andrew's monotone chain over (x, y): the strictly convex hull,
/// counter-clockwise from the lowest-leftmost point, as input indices (the
/// first of any coincident points). Fewer than three entries: collinear.
std::vector<std::int64_t> hull_chain_2d(const std::vector<HullVec>& rP) {
    std::vector<std::int64_t> order(rP.size());
    for (std::size_t i = 0; i < order.size(); ++i)
        order[i] = static_cast<std::int64_t>(i);
    std::sort(order.begin(), order.end(), [&](std::int64_t A, std::int64_t B) {
        const HullVec& a = rP[static_cast<std::size_t>(A)];
        const HullVec& b = rP[static_cast<std::size_t>(B)];
        if (a[0] != b[0])
            return a[0] < b[0];
        if (a[1] != b[1])
            return a[1] < b[1];
        return A < B;
    });
    std::vector<std::int64_t> unique;
    unique.reserve(order.size());
    for (const std::int64_t i : order) {
        if (!unique.empty()) {
            const HullVec& last = rP[static_cast<std::size_t>(unique.back())];
            const HullVec& p = rP[static_cast<std::size_t>(i)];
            if (last[0] == p[0] && last[1] == p[1])
                continue;
        }
        unique.push_back(i);
    }
    if (unique.size() < 3)
        return unique;
    const auto turn = [&](std::int64_t A, std::int64_t B, std::int64_t C) {
        return detail::orient2d(rP[static_cast<std::size_t>(A)].data(),
                                rP[static_cast<std::size_t>(B)].data(),
                                rP[static_cast<std::size_t>(C)].data());
    };
    std::vector<std::int64_t> hull(2 * unique.size());
    std::size_t k = 0;
    for (const std::int64_t i : unique) {  // lower chain
        while (k >= 2 && turn(hull[k - 2], hull[k - 1], i) <= 0.0)
            --k;
        hull[k++] = i;
    }
    for (std::size_t j = unique.size() - 1, lower = k + 1; j-- > 0;) {  // upper chain
        const std::int64_t i = unique[j];
        while (k >= lower && turn(hull[k - 2], hull[k - 1], i) <= 0.0)
            --k;
        hull[k++] = i;
    }
    hull.resize(k - 1);  // the last point repeats the first
    return hull;
}

/// A hull face: counter-clockwise seen from outside; `mN[e]` the face across
/// edge `e` = (`mV[e]`, `mV[(e + 1) % 3]`).
struct HullFace {
    std::array<std::int64_t, 3> mV{};
    std::array<std::int64_t, 3> mN{-1, -1, -1};
    std::vector<std::int64_t> mOutside;
    bool mAlive = true;
    /// The eye that last tested this face, and whether it sees it (a stamp,
    /// so no per-step array over every face ever made).
    std::int64_t mTestedBy = -1;
    bool mVisible = false;
};

class HullBuilder {
public:
    explicit HullBuilder(const std::vector<HullVec>& rP) : mP(rP) {}

    /// Strictly outside the face's plane (exact).
    bool Sees(const HullFace& rF, std::int64_t P) const {
        return detail::orient3d(At(rF.mV[0]), At(rF.mV[1]), At(rF.mV[2]), At(P)) < 0.0;
    }

    void Build() {
        std::array<std::int64_t, 4> s = InitialSimplex();
        // Each face of the tetrahedron, wound so the opposite vertex is
        // inside (orient3d > 0: below the face).
        const std::array<std::array<int, 4>, 4> kFaces{
            {{0, 1, 2, 3}, {0, 3, 1, 2}, {1, 3, 2, 0}, {0, 2, 3, 1}}};
        for (const auto& f : kFaces) {
            HullFace face;
            face.mV = {s[static_cast<std::size_t>(f[0])], s[static_cast<std::size_t>(f[1])],
                       s[static_cast<std::size_t>(f[2])]};
            if (detail::orient3d(At(face.mV[0]), At(face.mV[1]), At(face.mV[2]),
                                 At(s[static_cast<std::size_t>(f[3])])) < 0.0)
                std::swap(face.mV[1], face.mV[2]);
            mFaces.push_back(std::move(face));
        }
        LinkAll();
        std::vector<char> in_simplex(mP.size(), 0);
        for (const std::int64_t v : s)
            in_simplex[static_cast<std::size_t>(v)] = 1;
        for (std::size_t i = 0; i < mP.size(); ++i) {
            if (in_simplex[i])
                continue;
            for (auto& face : mFaces)
                if (Sees(face, static_cast<std::int64_t>(i))) {
                    face.mOutside.push_back(static_cast<std::int64_t>(i));
                    break;
                }
        }
        std::deque<std::int64_t> queue;
        for (std::size_t f = 0; f < mFaces.size(); ++f)
            queue.push_back(static_cast<std::int64_t>(f));
        while (!queue.empty()) {
            const std::int64_t f = queue.front();
            queue.pop_front();
            if (!mFaces[static_cast<std::size_t>(f)].mAlive ||
                mFaces[static_cast<std::size_t>(f)].mOutside.empty())
                continue;
            for (const std::int64_t created : AddPoint(f))
                queue.push_back(created);
        }
    }

    const std::vector<HullFace>& Faces() const { return mFaces; }

private:
    const std::vector<HullVec>& mP;
    std::vector<HullFace> mFaces;

    const double* At(std::int64_t I) const { return mP[static_cast<std::size_t>(I)].data(); }

    bool Collinear(std::int64_t A, std::int64_t B, std::int64_t C) const {
        // Collinear in 3-D exactly when every coordinate projection is.
        for (int drop = 0; drop < 3; ++drop) {
            const int u = (drop + 1) % 3, v = (drop + 2) % 3;
            const double a[2] = {At(A)[u], At(A)[v]};
            const double b[2] = {At(B)[u], At(B)[v]};
            const double c[2] = {At(C)[u], At(C)[v]};
            if (detail::orient2d(a, b, c) != 0.0)
                return false;
        }
        return true;
    }

    std::array<std::int64_t, 4> InitialSimplex() const {
        const std::size_t n = mP.size();
        const auto less = [&](std::size_t A, std::size_t B) { return mP[A] < mP[B]; };
        std::size_t lo = 0, hi = 0;
        for (std::size_t i = 1; i < n; ++i) {
            if (less(i, lo))
                lo = i;
            if (less(hi, i))
                hi = i;
        }
        if (mP[lo] == mP[hi])
            throw std::invalid_argument(
                "meshio++: convex_hull: the points all coincide; they span no hull");
        const auto a = static_cast<std::int64_t>(lo), b = static_cast<std::int64_t>(hi);
        // The point farthest from the line (a, b) that is exactly off it.
        const HullVec d = detail::MathUtils::Sub(mP[hi], mP[lo]);
        std::int64_t c = -1;
        double best = -1.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double dist = detail::MathUtils::NormSq(
                detail::MathUtils::Cross(d, detail::MathUtils::Sub(mP[i], mP[lo])));
            if (dist > best && !Collinear(a, b, static_cast<std::int64_t>(i))) {
                best = dist;
                c = static_cast<std::int64_t>(i);
            }
        }
        if (c < 0)
            throw std::invalid_argument(
                "meshio++: convex_hull: the points are collinear; they span no hull");
        // The point farthest from the plane (a, b, c) that is exactly off it.
        std::int64_t e = -1;
        best = -1.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double o = detail::orient3d(At(a), At(b), At(c), mP[i].data());
            if (o != 0.0 && std::abs(o) > best) {
                best = std::abs(o);
                e = static_cast<std::int64_t>(i);
            }
        }
        if (e < 0)
            throw std::invalid_argument(
                "meshio++: convex_hull: the points are coplanar, so they span no 3-D hull; "
                "pass dimension=2 for points in the x, y plane");
        return {a, b, c, e};
    }

    /// Neighbour links of the initial tetrahedron (four faces, so a plain search).
    void LinkAll() {
        for (std::size_t f = 0; f < mFaces.size(); ++f)
            for (std::size_t e = 0; e < 3; ++e) {
                const std::int64_t from = mFaces[f].mV[e], to = mFaces[f].mV[(e + 1) % 3];
                for (std::size_t g = 0; g < mFaces.size(); ++g)
                    for (std::size_t h = 0; h < 3; ++h)
                        if (mFaces[g].mV[h] == to && mFaces[g].mV[(h + 1) % 3] == from)
                            mFaces[f].mN[e] = static_cast<std::int64_t>(g);
            }
    }

    /// Adds the farthest outside point of face `F`; returns the new faces.
    std::vector<std::int64_t> AddPoint(std::int64_t F) {
        HullFace& start = mFaces[static_cast<std::size_t>(F)];
        // The farthest outside point (floating distance; the visibility
        // decisions below are exact), the lowest index on a tie.
        const double* a = At(start.mV[0]);
        const double* b = At(start.mV[1]);
        const double* c = At(start.mV[2]);
        const HullVec nrm = detail::GeometryUtils::TriangleNormal(a, b, c);
        std::int64_t eye = -1;
        double best = -std::numeric_limits<double>::infinity();
        for (const std::int64_t p : start.mOutside) {
            const double* q = At(p);
            const double dist = detail::MathUtils::Dot(nrm, detail::MathUtils::Sub(q, a));
            if (dist > best || (dist == best && p < eye)) {
                best = dist;
                eye = p;
            }
        }
        // The faces the eye sees: connected, found from F across edges.
        std::vector<std::int64_t> visible{F};
        start.mTestedBy = eye;
        start.mVisible = true;
        for (std::size_t k = 0; k < visible.size(); ++k) {
            const std::array<std::int64_t, 3> around =
                mFaces[static_cast<std::size_t>(visible[k])].mN;
            for (const std::int64_t nb : around) {
                HullFace& other = mFaces[static_cast<std::size_t>(nb)];
                if (other.mTestedBy == eye)
                    continue;
                other.mTestedBy = eye;
                other.mVisible = Sees(other, eye);
                if (other.mVisible)
                    visible.push_back(nb);
            }
        }
        // The horizon: edges of visible faces whose neighbour is not; a new
        // face (a, b, eye) on each keeps the visible face's winding.
        std::vector<std::int64_t> created;
        std::unordered_map<std::int64_t, std::int64_t> by_start, by_end;
        for (const std::int64_t vf : visible) {
            for (int e = 0; e < 3; ++e) {
                const HullFace& face = mFaces[static_cast<std::size_t>(vf)];
                const std::int64_t nb = face.mN[static_cast<std::size_t>(e)];
                if (mFaces[static_cast<std::size_t>(nb)].mVisible &&
                    mFaces[static_cast<std::size_t>(nb)].mTestedBy == eye)
                    continue;
                HullFace nf;
                nf.mV = {face.mV[static_cast<std::size_t>(e)],
                         face.mV[static_cast<std::size_t>((e + 1) % 3)], eye};
                nf.mN[0] = nb;
                nf.mTestedBy = eye;
                const auto id = static_cast<std::int64_t>(mFaces.size());
                HullFace& other = mFaces[static_cast<std::size_t>(nb)];
                for (int oe = 0; oe < 3; ++oe)
                    if (other.mN[static_cast<std::size_t>(oe)] == vf &&
                        other.mV[static_cast<std::size_t>(oe)] == nf.mV[1])
                        other.mN[static_cast<std::size_t>(oe)] = id;
                if (!by_start.emplace(nf.mV[0], id).second || !by_end.emplace(nf.mV[1], id).second)
                    throw std::runtime_error(
                        "meshio++: convex_hull: internal error: the horizon is not a simple "
                        "cycle");
                mFaces.push_back(std::move(nf));
                created.push_back(id);
            }
        }
        for (const std::int64_t id : created) {
            HullFace& nf = mFaces[static_cast<std::size_t>(id)];
            nf.mN[1] = by_start.at(nf.mV[1]);  // across (b, eye): the face starting at b
            nf.mN[2] = by_end.at(nf.mV[0]);    // across (eye, a): the face ending at a
        }
        // The visible faces' outside points go to the first new face that
        // sees them; the rest are inside now.
        for (const std::int64_t vf : visible) {
            HullFace& face = mFaces[static_cast<std::size_t>(vf)];
            face.mAlive = false;
            std::vector<std::int64_t> outside;
            outside.swap(face.mOutside);
            for (const std::int64_t p : outside) {
                if (p == eye)
                    continue;
                for (const std::int64_t id : created) {
                    HullFace& nf = mFaces[static_cast<std::size_t>(id)];
                    if (Sees(nf, p)) {
                        nf.mOutside.push_back(p);
                        break;
                    }
                }
            }
        }
        return created;
    }
};

/// The hull's faces regrouped into its facets: adjacent faces in one plane
/// (exact `orient3d` == 0) form a convex polygon, reduced to its strict
/// corners (a point inside it or on its edges, which quickhull may have made
/// a vertex while the hull was still growing, is not extreme) and fanned
/// into triangles, wound outward.
std::vector<std::array<std::int64_t, 3>> hull_facets(const std::vector<HullVec>& rP,
                                                     const std::vector<HullFace>& rFaces) {
    std::vector<std::int64_t> parent(rFaces.size());
    for (std::size_t f = 0; f < parent.size(); ++f)
        parent[f] = static_cast<std::int64_t>(f);
    const auto find = [&](std::int64_t F) {
        while (parent[static_cast<std::size_t>(F)] != F) {
            parent[static_cast<std::size_t>(F)] =
                parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(F)])];
            F = parent[static_cast<std::size_t>(F)];
        }
        return F;
    };
    for (std::size_t f = 0; f < rFaces.size(); ++f) {
        const HullFace& face = rFaces[f];
        if (!face.mAlive)
            continue;
        for (std::size_t e = 0; e < 3; ++e) {
            const HullFace& other = rFaces[static_cast<std::size_t>(face.mN[e])];
            std::int64_t opposite = -1;
            for (const std::int64_t v : other.mV)
                if (v != face.mV[e] && v != face.mV[(e + 1) % 3])
                    opposite = v;
            if (detail::orient3d(rP[static_cast<std::size_t>(face.mV[0])].data(),
                                 rP[static_cast<std::size_t>(face.mV[1])].data(),
                                 rP[static_cast<std::size_t>(face.mV[2])].data(),
                                 rP[static_cast<std::size_t>(opposite)].data()) == 0.0) {
                const std::int64_t a = find(static_cast<std::int64_t>(f)), b = find(face.mN[e]);
                parent[static_cast<std::size_t>(std::max(a, b))] = std::min(a, b);
            }
        }
    }
    // Facets in the order of their first face; their vertices, ascending.
    std::vector<std::int64_t> roots;
    std::unordered_map<std::int64_t, std::vector<std::int64_t>> members;
    for (std::size_t f = 0; f < rFaces.size(); ++f) {
        if (!rFaces[f].mAlive)
            continue;
        const std::int64_t r = find(static_cast<std::int64_t>(f));
        auto& list = members[r];
        if (list.empty())
            roots.push_back(r);
        list.insert(list.end(), rFaces[f].mV.begin(), rFaces[f].mV.end());
    }
    std::vector<std::array<std::int64_t, 3>> out;
    for (const std::int64_t r : roots) {
        std::vector<std::int64_t>& verts = members[r];
        std::sort(verts.begin(), verts.end());
        verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
        const HullFace& first = rFaces[static_cast<std::size_t>(r)];
        if (verts.size() == 3) {
            out.push_back(first.mV);
            continue;
        }
        // Project along the normal's largest component, keeping the winding
        // (counter-clockwise from outside is counter-clockwise in (u, v)).
        const double* a = rP[static_cast<std::size_t>(first.mV[0])].data();
        const double* b = rP[static_cast<std::size_t>(first.mV[1])].data();
        const double* c = rP[static_cast<std::size_t>(first.mV[2])].data();
        const HullVec n{(b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]),
                        (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]),
                        (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])};
        std::size_t k = 0;
        for (std::size_t j = 1; j < 3; ++j)
            if (std::abs(n[j]) > std::abs(n[k]))
                k = j;
        std::size_t u = (k + 1) % 3, v = (k + 2) % 3;
        if (n[k] < 0.0)
            std::swap(u, v);
        std::vector<HullVec> plane(verts.size());
        for (std::size_t i = 0; i < verts.size(); ++i) {
            const HullVec& q = rP[static_cast<std::size_t>(verts[i])];
            plane[i] = {q[u], q[v], 0.0};
        }
        const std::vector<std::int64_t> chain = hull_chain_2d(plane);
        for (std::size_t i = 1; i + 1 < chain.size(); ++i)
            out.push_back({verts[static_cast<std::size_t>(chain[0])],
                           verts[static_cast<std::size_t>(chain[i])],
                           verts[static_cast<std::size_t>(chain[i + 1])]});
    }
    return out;
}

/// The hull vertices' output mesh: their points, gathered point data,
/// `hull:point_id`, field data.
Mesh hull_vertex_mesh(const Mesh& rMesh, const std::vector<std::int64_t>& rIds) {
    Mesh out;
    out.AssignPoints(detail::subset_gather_rows(rMesh.Points(), rIds));
    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::subset_gather_rows(rMesh.PointData(name), rIds));
    NDArray ids(DType::Int64, {rIds.size()});
    std::copy(rIds.begin(), rIds.end(), ids.As<std::int64_t>());
    out.AddPointData(kHullPointIdName, std::move(ids));
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));
    return out;
}

}  // namespace

ConvexHullResult convex_hull(const Mesh& rMesh, const ConvexHullOptions& rOptions) {
    if (rOptions.mDimension != 0 && rOptions.mDimension != 2 && rOptions.mDimension != 3)
        throw std::invalid_argument(
            "meshio++: convex_hull: dimension must be 0 (automatic), 2 "
            "or 3, not " +
            std::to_string(rOptions.mDimension));
    const std::vector<HullVec> p = hull_points(rMesh, "convex_hull");
    if (p.size() < 3)
        throw std::invalid_argument("meshio++: convex_hull: needs at least three points, not " +
                                    std::to_string(p.size()));
    int dimension = rOptions.mDimension;
    if (dimension == 0) {
        dimension = 2;
        if (rMesh.PointDim() >= 3)
            for (const HullVec& q : p)
                if (q[2] != p[0][2]) {
                    dimension = 3;
                    break;
                }
    }
    if (dimension == 3 && rMesh.PointDim() < 3)
        throw std::invalid_argument(
            "meshio++: convex_hull: dimension 3 needs points with three coordinates");

    ConvexHullResult result;
    result.mDimension = dimension;
    if (dimension == 2) {
        const std::vector<std::int64_t> chain = hull_chain_2d(p);
        if (chain.size() < 3)
            throw std::invalid_argument(
                "meshio++: convex_hull: the points are collinear in x, y; they span no hull");
        // Vertices in ascending input order; the loop follows the chain.
        std::vector<std::int64_t> ids = chain;
        std::sort(ids.begin(), ids.end());
        std::unordered_map<std::int64_t, std::int64_t> local;
        for (std::size_t i = 0; i < ids.size(); ++i)
            local[ids[i]] = static_cast<std::int64_t>(i);
        const std::size_t m = chain.size();
        NDArray conn(DType::Int64, {m, std::size_t{2}});
        std::int64_t* pc = conn.As<std::int64_t>();
        double twice_area = 0.0;
        for (std::size_t i = 0; i < m; ++i) {
            const std::int64_t a = chain[i], b = chain[(i + 1) % m];
            pc[2 * i] = local.at(a);
            pc[2 * i + 1] = local.at(b);
            const HullVec& pa = p[static_cast<std::size_t>(a)];
            const HullVec& pb = p[static_cast<std::size_t>(b)];
            result.mPerimeter += std::hypot(pb[0] - pa[0], pb[1] - pa[1]);
            twice_area += pa[0] * pb[1] - pb[0] * pa[1];
        }
        result.mArea = 0.5 * twice_area;
        result.mNumVertices = static_cast<std::int64_t>(m);
        result.mNumFacets = static_cast<std::int64_t>(m);
        result.mMesh = hull_vertex_mesh(rMesh, ids);
        result.mMesh.AddCellBlock("line", std::move(conn));
        return result;
    }

    HullBuilder builder(p);
    builder.Build();
    const std::vector<std::array<std::int64_t, 3>> tris = hull_facets(p, builder.Faces());
    std::vector<std::int64_t> ids;
    for (const auto& t : tris)
        ids.insert(ids.end(), t.begin(), t.end());
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    std::unordered_map<std::int64_t, std::int64_t> local;
    for (std::size_t i = 0; i < ids.size(); ++i)
        local[ids[i]] = static_cast<std::int64_t>(i);
    NDArray conn(DType::Int64, {tris.size(), std::size_t{3}});
    std::int64_t* pc = conn.As<std::int64_t>();
    const HullVec& o = p[static_cast<std::size_t>(ids.front())];
    double six_volume = 0.0;
    for (std::size_t t = 0; t < tris.size(); ++t) {
        const HullVec& a = p[static_cast<std::size_t>(tris[t][0])];
        const HullVec& b = p[static_cast<std::size_t>(tris[t][1])];
        const HullVec& c = p[static_cast<std::size_t>(tris[t][2])];
        for (int k = 0; k < 3; ++k)
            pc[3 * t + static_cast<std::size_t>(k)] =
                local.at(tris[t][static_cast<std::size_t>(k)]);
        const HullVec x = detail::GeometryUtils::TriangleNormal(a, b, c);
        result.mArea += 0.5 * detail::MathUtils::Norm(x);
        six_volume += detail::MathUtils::Dot(x, detail::MathUtils::Sub(a, o));
    }
    result.mVolume = six_volume / 6.0;
    result.mNumVertices = static_cast<std::int64_t>(ids.size());
    result.mNumFacets = static_cast<std::int64_t>(tris.size());
    result.mMesh = hull_vertex_mesh(rMesh, ids);
    result.mMesh.AddCellBlock("triangle", std::move(conn));
    return result;
}

Delaunay2dResult delaunay_2d(const Mesh& rMesh) {
    const std::vector<HullVec> p = hull_points(rMesh, "delaunay_2d");
    double lo[2] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[2] = {-lo[0], -lo[1]};
    for (std::size_t i = 0; i < p.size(); ++i)
        for (int k = 0; k < 2; ++k) {
            if (!detail::cdt2d_in_range(p[i][static_cast<std::size_t>(k)]))
                throw std::invalid_argument(
                    "meshio++: delaunay_2d: point " + std::to_string(i) +
                    " has a coordinate outside the exact range (zero, or a magnitude in "
                    "[2^-128, 2^128]); rescale the points first");
            lo[k] = std::min(lo[k], p[i][static_cast<std::size_t>(k)]);
            hi[k] = std::max(hi[k], p[i][static_cast<std::size_t>(k)]);
        }
    const std::vector<std::int64_t> chain = hull_chain_2d(p);
    if (chain.size() < 3)
        throw std::invalid_argument(
            "meshio++: delaunay_2d: needs at least three distinct points that are not collinear "
            "in x, y");

    detail::Cdt2d cdt(lo[0], lo[1], hi[0], hi[1]);
    std::vector<detail::Cdt2d::Vec2> xy(p.size());
    for (std::size_t i = 0; i < p.size(); ++i)
        xy[i] = {p[i][0], p[i][1]};
    const std::vector<std::int64_t> vid = cdt.InsertBatch(xy);
    // Vertex id -> the first input point at it.
    std::vector<std::int64_t> input_of(static_cast<std::size_t>(cdt.NumVertices()), -1);
    Delaunay2dResult result;
    for (std::size_t i = 0; i < vid.size(); ++i) {
        std::int64_t& slot = input_of[static_cast<std::size_t>(vid[i])];
        if (slot < 0)
            slot = static_cast<std::int64_t>(i);
        else
            ++result.mNumDuplicates;
    }
    for (std::size_t i = 0; i < chain.size(); ++i) {
        const detail::Cdt2dStatus status =
            cdt.InsertSegment(vid[static_cast<std::size_t>(chain[i])],
                              vid[static_cast<std::size_t>(chain[(i + 1) % chain.size()])]);
        if (status != detail::Cdt2dStatus::Ok)
            throw std::runtime_error(std::string("meshio++: delaunay_2d: internal error: a hull "
                                                 "edge was refused (") +
                                     detail::cdt2d_status_name(status) + ")");
    }
    const detail::Cdt2dStatus status = cdt.Classify();
    if (status != detail::Cdt2dStatus::Ok)
        throw std::runtime_error(std::string("meshio++: delaunay_2d: internal error: ") +
                                 detail::cdt2d_status_name(status));
    const std::vector<detail::Cdt2d::Tri> tris = cdt.InsideTriangles();
    NDArray conn(DType::Int64, {tris.size(), std::size_t{3}});
    std::int64_t* pc = conn.As<std::int64_t>();
    for (std::size_t t = 0; t < tris.size(); ++t)
        for (std::size_t k = 0; k < 3; ++k) {
            const std::int64_t v = input_of[static_cast<std::size_t>(tris[t][k])];
            if (v < 0)
                throw std::runtime_error(
                    "meshio++: delaunay_2d: internal error: a triangle uses a point not in the "
                    "input");
            pc[3 * t + k] = v;
        }
    // Boundary points: the ends of the constrained (hull) edges, which the
    // triangulation split at every point on them.
    std::vector<char> on_boundary(p.size(), 0);
    for (const auto& e : cdt.ConstrainedEdges())
        for (const std::int64_t v : e)
            on_boundary[static_cast<std::size_t>(input_of[static_cast<std::size_t>(v)])] = 1;
    for (const char b : on_boundary)
        result.mNumBoundary += b;
    result.mNumTriangles = static_cast<std::int64_t>(tris.size());

    Mesh out;
    out.AssignPoints(detail::data_owned_copy(rMesh.Points()));
    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::data_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));
    out.AddCellBlock("triangle", std::move(conn));
    result.mMesh = std::move(out);
    return result;
}

}  // namespace meshioplusplus
