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
// Vertex normals of a triangle soup and the corner-group split. See
// detail/surface_normals.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

// Project includes
#include "meshioplusplus/detail/surface_normals.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr double kSnPi = 3.14159265358979323846;

// One edge use of one triangle: the sorted endpoint pair, the triangle, the
// corner index (3 * triangle + corner) of each endpoint, and whether the
// triangle walks it low -> high.
struct SnEdgeUse {
    std::int64_t mLo;
    std::int64_t mHi;
    std::int64_t mTri;
    std::int64_t mCornerLo;
    std::int64_t mCornerHi;
    bool mForward;
};

bool sn_edge_use_less(const SnEdgeUse& rA, const SnEdgeUse& rB) {
    if (rA.mLo != rB.mLo)
        return rA.mLo < rB.mLo;
    if (rA.mHi != rB.mHi)
        return rA.mHi < rB.mHi;
    if (rA.mTri != rB.mTri)
        return rA.mTri < rB.mTri;
    return rA.mCornerLo < rB.mCornerLo;
}

// Union-find whose root is always the smallest index of its set. That is what
// makes the partition, and the group numbering derived from it, independent of
// the order the joins are applied in.
class SnUnionFind {
public:
    explicit SnUnionFind(std::size_t N) : mParent(N) {
        std::iota(mParent.begin(), mParent.end(), std::int64_t{0});
    }
    std::int64_t Find(std::int64_t X) {
        while (mParent[static_cast<std::size_t>(X)] != X) {
            std::int64_t& p = mParent[static_cast<std::size_t>(X)];
            p = mParent[static_cast<std::size_t>(p)];
            X = p;
        }
        return X;
    }
    void Unite(std::int64_t A, std::int64_t B) {
        const std::int64_t ra = Find(A);
        const std::int64_t rb = Find(B);
        if (ra == rb)
            return;
        if (ra < rb)
            mParent[static_cast<std::size_t>(rb)] = ra;
        else
            mParent[static_cast<std::size_t>(ra)] = rb;
    }

private:
    std::vector<std::int64_t> mParent;
};

double sn_corner_weight(const TriangleSoup& rSoup, std::size_t Tri, std::size_t Corner, double Len,
                        SdfPseudonormalWeight Weight) {
    if (Weight != SdfPseudonormalWeight::Angle)
        return Len;  // area weighting: |cross| is twice the area, a positive scale
    const Vec3* pCorner = &rSoup.mCorners[Tri * 3];
    return corner_angle(pCorner[Corner], pCorner[(Corner + 1) % 3], pCorner[(Corner + 2) % 3]);
}

}  // namespace

double corner_angle(const Vec3& rA, const Vec3& rB, const Vec3& rC) {
    const Vec3 u = vec3_sub(rB, rA);
    const Vec3 v = vec3_sub(rC, rA);
    const double nu = vec3_norm(u);
    const double nv = vec3_norm(v);
    if (!(nu > 0.0) || !(nv > 0.0))
        return 0.0;
    double c = vec3_dot(u, v) / (nu * nv);
    c = c < -1.0 ? -1.0 : (c > 1.0 ? 1.0 : c);
    return std::acos(c);
}

std::vector<Vec3> soup_face_normals(const TriangleSoup& rSoup) {
    std::vector<Vec3> normals(rSoup.NumTriangles());
    parallel_for(normals.size(), [&](std::size_t t) {
        const Vec3& a = rSoup.mCorners[t * 3 + 0];
        const Vec3& b = rSoup.mCorners[t * 3 + 1];
        const Vec3& c = rSoup.mCorners[t * 3 + 2];
        normals[t] = vec3_cross(vec3_sub(b, a), vec3_sub(c, a));
    });
    return normals;
}

std::vector<Vec3> accumulate_vertex_normals(const TriangleSoup& rSoup,
                                            const std::vector<Vec3>& rFaceNormal,
                                            SdfPseudonormalWeight Weight) {
    // Gather form: each vertex sums its own corners, in ascending (triangle,
    // corner) order -- the order the serial scatter added them in, so the bits
    // are the same -- found through a counting sort of the corners by vertex.
    const std::size_t npts = rSoup.mPoints.size();
    const std::size_t ncorner = rSoup.NumTriangles() * 3;
    std::vector<std::int64_t> start(npts + 1, 0);
    for (std::size_t c = 0; c < ncorner; ++c)
        ++start[static_cast<std::size_t>(rSoup.mVertices[c / 3][c % 3]) + 1];
    for (std::size_t p = 0; p < npts; ++p)
        start[p + 1] += start[p];
    std::vector<std::int64_t> corners(ncorner);
    {
        std::vector<std::int64_t> cursor(start.begin(), start.end() - 1);
        for (std::size_t c = 0; c < ncorner; ++c)
            corners[static_cast<std::size_t>(
                cursor[static_cast<std::size_t>(rSoup.mVertices[c / 3][c % 3])]++)] =
                static_cast<std::int64_t>(c);
    }
    std::vector<Vec3> sums(npts, Vec3{0.0, 0.0, 0.0});
    parallel_for(npts, [&](std::size_t p) {
        Vec3 acc{0.0, 0.0, 0.0};
        for (std::int64_t k = start[p]; k < start[p + 1]; ++k) {
            const std::size_t c = static_cast<std::size_t>(corners[static_cast<std::size_t>(k)]);
            const std::size_t t = c / 3;
            const Vec3& n = rFaceNormal[t];
            const double len = vec3_norm(n);
            if (!(len > 0.0))
                continue;  // degenerate: no direction to contribute
            const Vec3 unit = vec3_scale(n, 1.0 / len);
            acc = vec3_add(acc, vec3_scale(unit, sn_corner_weight(rSoup, t, c % 3, len, Weight)));
        }
        sums[p] = acc;
    });
    return sums;
}

VertexNormalGroups vertex_normal_groups(const TriangleSoup& rSoup, SdfPseudonormalWeight Weight,
                                        double SplitAngleDeg) {
    const std::size_t ntri = rSoup.NumTriangles();
    const std::size_t ncorner = ntri * 3;
    const std::vector<Vec3> face = soup_face_normals(rSoup);

    std::vector<double> length(ntri, 0.0);
    std::vector<Vec3> unit(ntri, Vec3{0.0, 0.0, 0.0});
    VertexNormalGroups out;
    for (std::size_t t = 0; t < ntri; ++t) {
        length[t] = vec3_norm(face[t]);
        if (length[t] > 0.0)
            unit[t] = vec3_scale(face[t], 1.0 / length[t]);
        else
            ++out.mNumDegenerate;
    }

    SnUnionFind sets(ncorner);
    if (SplitAngleDeg < 0.0) {
        // No split: every corner of a vertex is one group.
        std::vector<std::int64_t> first(rSoup.mPoints.size(), -1);
        for (std::size_t c = 0; c < ncorner; ++c) {
            const std::size_t p = static_cast<std::size_t>(rSoup.mVertices[c / 3][c % 3]);
            if (first[p] < 0)
                first[p] = static_cast<std::int64_t>(c);
            else
                sets.Unite(first[p], static_cast<std::int64_t>(c));
        }
    } else {
        const bool always = SplitAngleDeg >= 180.0;
        const double cos_threshold = std::cos(SplitAngleDeg * (kSnPi / 180.0));

        std::vector<SnEdgeUse> uses;
        uses.reserve(ncorner);
        for (std::size_t t = 0; t < ntri; ++t) {
            const std::array<std::int64_t, 3>& v = rSoup.mVertices[t];
            for (std::size_t i = 0; i < 3; ++i) {
                const std::size_t j = (i + 1) % 3;
                const std::int64_t u = v[i];
                const std::int64_t w = v[j];
                if (u == w)
                    continue;  // a collapsed edge joins nothing
                const std::int64_t ci = static_cast<std::int64_t>(t * 3 + i);
                const std::int64_t cj = static_cast<std::int64_t>(t * 3 + j);
                if (u < w)
                    uses.push_back({u, w, static_cast<std::int64_t>(t), ci, cj, true});
                else
                    uses.push_back({w, u, static_cast<std::int64_t>(t), cj, ci, false});
            }
        }
        std::sort(uses.begin(), uses.end(), sn_edge_use_less);

        for (std::size_t b = 0; b < uses.size();) {
            std::size_t e = b + 1;
            while (e < uses.size() && uses[e].mLo == uses[b].mLo && uses[e].mHi == uses[b].mHi)
                ++e;
            const bool proper_pair = e - b == 2;
            for (std::size_t x = b; x < e; ++x) {
                for (std::size_t y = x + 1; y < e; ++y) {
                    const SnEdgeUse& ux = uses[x];
                    const SnEdgeUse& uy = uses[y];
                    if (ux.mTri == uy.mTri)
                        continue;
                    const std::size_t tx = static_cast<std::size_t>(ux.mTri);
                    const std::size_t ty = static_cast<std::size_t>(uy.mTri);
                    bool join = rSoup.mSourceCell[tx] == rSoup.mSourceCell[ty];
                    if (!join && proper_pair && ux.mForward != uy.mForward && length[tx] > 0.0 &&
                        length[ty] > 0.0)
                        join = always || vec3_dot(unit[tx], unit[ty]) >= cos_threshold;
                    if (join) {
                        sets.Unite(ux.mCornerLo, uy.mCornerLo);
                        sets.Unite(ux.mCornerHi, uy.mCornerHi);
                    }
                }
            }
            b = e;
        }
    }

    // Number the groups by ascending root. A root is the smallest corner of its
    // set, so it has always been numbered by the time a later corner asks.
    std::vector<std::int64_t> root(ncorner);
    std::vector<std::int64_t> id(ncorner, -1);
    out.mCornerGroup.assign(ncorner, 0);
    for (std::size_t c = 0; c < ncorner; ++c) {
        root[c] = sets.Find(static_cast<std::int64_t>(c));
        if (root[c] == static_cast<std::int64_t>(c)) {
            id[c] = static_cast<std::int64_t>(out.mGroupPoint.size());
            out.mGroupRoot.push_back(static_cast<std::int64_t>(c));
            out.mGroupPoint.push_back(rSoup.mVertices[c / 3][c % 3]);
        }
        out.mCornerGroup[c] = id[static_cast<std::size_t>(root[c])];
    }

    // Serial scatter in ascending (triangle, corner) order, then normalise.
    out.mGroupNormal.assign(out.NumGroups(), Vec3{0.0, 0.0, 0.0});
    for (std::size_t t = 0; t < ntri; ++t) {
        if (!(length[t] > 0.0))
            continue;
        for (std::size_t i = 0; i < 3; ++i) {
            Vec3& acc = out.mGroupNormal[static_cast<std::size_t>(out.mCornerGroup[t * 3 + i])];
            acc = vec3_add(acc,
                           vec3_scale(unit[t], sn_corner_weight(rSoup, t, i, length[t], Weight)));
        }
    }
    for (Vec3& n : out.mGroupNormal) {
        const double len = vec3_norm(n);
        n = len > 0.0 ? vec3_scale(n, 1.0 / len) : Vec3{0.0, 0.0, 0.0};
    }
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
