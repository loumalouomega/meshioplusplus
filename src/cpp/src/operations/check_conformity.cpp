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
// check_conformity. See operations/check_conformity.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/check_conformity.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_edges.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/point_triangle.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/surface_rings.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kCcPrefix = "meshio++: check_conformity: ";

// The unmatched facets of the examined cells, as corner rings (CSR), each with
// the global cell that owns it.
struct CcFacets {
    std::vector<std::int64_t> mStart{0};
    std::vector<std::int64_t> mNodes;
    std::vector<std::int64_t> mOwner;
    std::int64_t mNumNonManifold = 0;

    std::size_t Size() const { return mOwner.size(); }
    std::size_t Corners(std::size_t Facet) const {
        return static_cast<std::size_t>(mStart[Facet + 1] - mStart[Facet]);
    }
    const std::int64_t* Ring(std::size_t Facet) const { return mNodes.data() + mStart[Facet]; }
    void Push(const std::int64_t* pIds, std::size_t Size, std::int64_t Owner) {
        mNodes.insert(mNodes.end(), pIds, pIds + Size);
        mStart.push_back(static_cast<std::int64_t>(mNodes.size()));
        mOwner.push_back(Owner);
    }
};

// The topological dimension of a block: 3, 2, or -1 for anything else.
template <class TView>
int cc_block_dimension(const TView& rBlock) {
    if (rBlock.IsPolyhedron())
        return 3;
    const std::string type(rBlock.Type());
    if (rBlock.IsRagged() || type.rfind("polygon", 0) == 0)
        return 2;
    const int dim = cell_type_dimension(cell_type_from_name(type));
    return dim == 2 || dim == 3 ? dim : -1;
}

// Every node of a cell (all faces' nodes of a polyhedron), unsorted.
template <class TView>
void cc_cell_nodes(const TView& rBlock, std::size_t Row, std::vector<std::int64_t>& rOut) {
    rOut.clear();
    if (rBlock.IsPolyhedron()) {
        for (std::size_t f = 0; f < rBlock.NumFaces(Row); ++f) {
            const auto face = rBlock.Face(Row, f);
            rOut.insert(rOut.end(), face.first, face.first + face.second);
        }
    } else if (rBlock.IsRagged()) {
        rOut.assign(rBlock.Row(Row), rBlock.Row(Row) + rBlock.RowSize(Row));
    } else {
        const std::size_t npc = rBlock.NodesPerCell();
        for (std::size_t k = 0; k < npc; ++k)
            rOut.push_back(detail::read_int(rBlock.Conn(), Row * npc + k));
    }
}

// 3-D: the faces of the volume cells that only one cell uses.
CcFacets cc_volume_facets(const Mesh& rMesh, std::vector<std::string>& rWarnings) {
    CcFacets out;
    const detail::GlobalFaces gf = detail::build_global_faces(rMesh);
    out.mNumNonManifold = gf.mNumNonManifold;
    for (std::size_t f = 0; f < gf.NumFaces(); ++f)
        if (gf.mNeighbour[f] < 0)
            out.Push(gf.Face(f), gf.FaceSize(f),
                     gf.mCellToGlobal[static_cast<std::size_t>(gf.mOwner[f])]);
    std::size_t b = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (cc_block_dimension(cb) == 3 &&
            std::find(gf.mNonCellBlocks.begin(), gf.mNonCellBlocks.end(), b) !=
                gf.mNonCellBlocks.end() &&
            cb.NumCells() > 0)
            rWarnings.push_back("block " + std::to_string(b) + " ('" + std::string(cb.Type()) +
                                "') has no known face topology and was not examined");
        ++b;
    }
    if (gf.mNumUnorientable > 0)
        rWarnings.push_back(std::to_string(gf.mNumUnorientable) +
                            " polyhedra are not closed and orientable; their faces were "
                            "examined as given");
    return out;
}

// 2-D: the edges of the 2-D cells that only one cell uses, in ascending
// (low, high) endpoint order.
CcFacets cc_surface_facets(const Mesh& rMesh, std::vector<std::string>& rWarnings) {
    struct Use {
        std::int64_t mLo, mHi, mCell, mA, mB;
    };
    std::vector<Use> uses;
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t npts = rMesh.NumPoints();
    std::vector<std::int64_t> ring;
    std::size_t b = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t block = b++;
        if (cc_block_dimension(cb) != 2)
            continue;
        const std::string type(cb.Type());
        const bool polygon = cb.IsRagged() || type.rfind("polygon", 0) == 0;
        std::vector<std::size_t> corners;
        if (!polygon) {
            for (const detail::CellEdgeDef& ed : detail::cell_edges(cell_type_from_name(type)))
                corners.push_back(ed.mNodes[0]);
            if (corners.size() < 3) {
                if (cb.NumCells() > 0)
                    rWarnings.push_back("block " + std::to_string(block) + " ('" + type +
                                        "') has no known edge topology and was not examined");
                continue;
            }
        }
        const std::size_t npc = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            ring.clear();
            if (cb.IsRagged())
                ring.assign(cb.Row(c), cb.Row(c) + cb.RowSize(c));
            else if (polygon)
                for (std::size_t k = 0; k < npc; ++k)
                    ring.push_back(detail::read_int(cb.Conn(), c * npc + k));
            else
                for (std::size_t k : corners)
                    ring.push_back(detail::read_int(cb.Conn(), c * npc + k));
            const std::int64_t g = bases[block] + static_cast<std::int64_t>(c);
            const std::size_t n = ring.size();
            for (std::size_t k = 0; k < n; ++k) {
                const std::int64_t a = ring[k], bb = ring[(k + 1) % n];
                if (a < 0 || bb < 0 || static_cast<std::size_t>(a) >= npts ||
                    static_cast<std::size_t>(bb) >= npts)
                    throw std::invalid_argument(std::string(kCcPrefix) +
                                                "a cell references a point the mesh does not have");
                if (a == bb)
                    continue;
                uses.push_back(Use{std::min(a, bb), std::max(a, bb), g, a, bb});
            }
        }
    }
    std::sort(uses.begin(), uses.end(), [](const Use& rX, const Use& rY) {
        if (rX.mLo != rY.mLo)
            return rX.mLo < rY.mLo;
        if (rX.mHi != rY.mHi)
            return rX.mHi < rY.mHi;
        return rX.mCell < rY.mCell;
    });
    CcFacets out;
    for (std::size_t i = 0; i < uses.size();) {
        std::size_t j = i + 1;
        while (j < uses.size() && uses[j].mLo == uses[i].mLo && uses[j].mHi == uses[i].mHi)
            ++j;
        if (j - i == 1) {
            const std::int64_t ids[2] = {uses[i].mA, uses[i].mB};
            out.Push(ids, 2, uses[i].mCell);
        } else if (j - i >= 3) {
            ++out.mNumNonManifold;
        }
        i = j;
    }
    return out;
}

detail::Vec3 cc_point(const std::vector<double>& rXyz, std::int64_t Id) {
    const std::size_t i = static_cast<std::size_t>(Id) * 3;
    return {{rXyz[i], rXyz[i + 1], rXyz[i + 2]}};
}

// The squared distance from P to segment AB.
double cc_segment_distance_sq(const detail::Vec3& rP, const detail::Vec3& rA,
                              const detail::Vec3& rB) {
    const detail::Vec3 ab = detail::vec3_sub(rB, rA);
    const double len_sq = detail::vec3_norm_sq(ab);
    double t = len_sq > 0.0 ? detail::vec3_dot(detail::vec3_sub(rP, rA), ab) / len_sq : 0.0;
    t = std::min(1.0, std::max(0.0, t));
    return detail::vec3_norm_sq(
        detail::vec3_sub(rP, detail::vec3_add(rA, detail::vec3_scale(ab, t))));
}

// The squared distance from P to a facet: a segment, or a polygon as the fan
// of triangles from its first corner.
double cc_facet_distance_sq(const detail::Vec3& rP, const std::vector<double>& rXyz,
                            const std::int64_t* pRing, std::size_t Corners) {
    if (Corners == 2)
        return cc_segment_distance_sq(rP, cc_point(rXyz, pRing[0]), cc_point(rXyz, pRing[1]));
    double best = std::numeric_limits<double>::infinity();
    const detail::Vec3 a = cc_point(rXyz, pRing[0]);
    for (std::size_t k = 1; k + 1 < Corners; ++k)
        best = std::min(best, detail::closest_point_on_triangle(rP, a, cc_point(rXyz, pRing[k]),
                                                                cc_point(rXyz, pRing[k + 1]))
                                  .mDistanceSq);
    return best;
}

// A facet's plane: its centroid, unit (Newell) normal and an in-plane frame.
struct CcPlane {
    detail::Vec3 mCentre{{0.0, 0.0, 0.0}};
    detail::Vec3 mNormal{{0.0, 0.0, 0.0}};
    detail::Vec3 mU{{0.0, 0.0, 0.0}};
    detail::Vec3 mV{{0.0, 0.0, 0.0}};
    bool mValid = false;
};

CcPlane cc_plane(const std::vector<double>& rXyz, const std::int64_t* pRing, std::size_t Corners) {
    CcPlane plane;
    for (std::size_t k = 0; k < Corners; ++k)
        plane.mCentre = detail::vec3_add(plane.mCentre, cc_point(rXyz, pRing[k]));
    plane.mCentre = detail::vec3_scale(plane.mCentre, 1.0 / static_cast<double>(Corners));
    detail::Vec3 n{{0.0, 0.0, 0.0}};
    for (std::size_t k = 0; k < Corners; ++k) {
        const detail::Vec3 a = cc_point(rXyz, pRing[k]);
        const detail::Vec3 b = cc_point(rXyz, pRing[(k + 1) % Corners]);
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    const double len = detail::vec3_norm(n);
    if (!(len > 0.0))
        return plane;
    plane.mNormal = detail::vec3_scale(n, 1.0 / len);
    const detail::Vec3 e = detail::vec3_sub(cc_point(rXyz, pRing[1]), cc_point(rXyz, pRing[0]));
    const double elen = detail::vec3_norm(e);
    if (!(elen > 0.0))
        return plane;
    plane.mU = detail::vec3_scale(e, 1.0 / elen);
    plane.mV = detail::vec3_cross(plane.mNormal, plane.mU);
    plane.mValid = true;
    return plane;
}

using CcUv = std::array<double, 2>;

double cc_area(const std::vector<CcUv>& rPoly) {
    double a = 0.0;
    for (std::size_t i = 0, n = rPoly.size(); i < n; ++i) {
        const CcUv& p = rPoly[i];
        const CcUv& q = rPoly[(i + 1) % n];
        a += p[0] * q[1] - q[0] * p[1];
    }
    return 0.5 * a;
}

// The area of the intersection of two triangles, both counter-clockwise, by
// clipping the first against the second's three edges (Sutherland and Hodgman).
double cc_triangle_overlap(const std::array<CcUv, 3>& rA, const std::array<CcUv, 3>& rB) {
    std::vector<CcUv> poly(rA.begin(), rA.end()), next;
    for (std::size_t e = 0; e < 3 && !poly.empty(); ++e) {
        const CcUv& p = rB[e];
        const CcUv& q = rB[(e + 1) % 3];
        auto side = [&](const CcUv& rX) {
            return (q[0] - p[0]) * (rX[1] - p[1]) - (q[1] - p[1]) * (rX[0] - p[0]);
        };
        next.clear();
        for (std::size_t i = 0, n = poly.size(); i < n; ++i) {
            const CcUv& a = poly[i];
            const CcUv& b = poly[(i + 1) % n];
            const double sa = side(a), sb = side(b);
            if (sa >= 0.0)
                next.push_back(a);
            if ((sa >= 0.0) != (sb >= 0.0)) {
                const double t = sa / (sa - sb);
                next.push_back({{a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])}});
            }
        }
        poly.swap(next);
    }
    return poly.size() < 3 ? 0.0 : std::abs(cc_area(poly));
}

// A ring as the counter-clockwise fan of triangles from its first corner.
std::vector<std::array<CcUv, 3>> cc_fan(const std::vector<CcUv>& rRing) {
    std::vector<std::array<CcUv, 3>> fan;
    for (std::size_t k = 1; k + 1 < rRing.size(); ++k) {
        std::array<CcUv, 3> t{{rRing[0], rRing[k], rRing[k + 1]}};
        const double a =
            (t[1][0] - t[0][0]) * (t[2][1] - t[0][1]) - (t[1][1] - t[0][1]) * (t[2][0] - t[0][0]);
        if (a < 0.0)
            std::swap(t[1], t[2]);
        if (a != 0.0)
            fan.push_back(t);
    }
    return fan;
}

// How two unmatched facets relate.
enum class CcPair { Apart, Unwelded, Overlap };

CcPair cc_classify_polygons(const std::vector<double>& rXyz, const std::int64_t* pF, std::size_t Nf,
                            const CcPlane& rPf, const std::int64_t* pG, std::size_t Ng,
                            const CcPlane& rPg, double Tol) {
    if (!rPf.mValid || !rPg.mValid)
        return CcPair::Apart;
    for (std::size_t k = 0; k < Ng; ++k)
        if (std::abs(detail::vec3_dot(rPf.mNormal,
                                      detail::vec3_sub(cc_point(rXyz, pG[k]), rPf.mCentre))) > Tol)
            return CcPair::Apart;
    for (std::size_t k = 0; k < Nf; ++k)
        if (std::abs(detail::vec3_dot(rPg.mNormal,
                                      detail::vec3_sub(cc_point(rXyz, pF[k]), rPg.mCentre))) > Tol)
            return CcPair::Apart;
    if (Nf == Ng) {
        bool same = true;
        for (std::size_t i = 0; i < Nf && same; ++i) {
            bool found = false;
            for (std::size_t j = 0; j < Ng && !found; ++j)
                found = detail::vec3_norm_sq(detail::vec3_sub(cc_point(rXyz, pF[i]),
                                                              cc_point(rXyz, pG[j]))) <= Tol * Tol;
            same = found;
        }
        if (same)
            return CcPair::Unwelded;
    }
    auto uv = [&](std::int64_t Id) {
        const detail::Vec3 d = detail::vec3_sub(cc_point(rXyz, Id), rPf.mCentre);
        return CcUv{{detail::vec3_dot(d, rPf.mU), detail::vec3_dot(d, rPf.mV)}};
    };
    std::vector<CcUv> f(Nf), g(Ng);
    for (std::size_t k = 0; k < Nf; ++k)
        f[k] = uv(pF[k]);
    for (std::size_t k = 0; k < Ng; ++k)
        g[k] = uv(pG[k]);
    // The facets overlap when their intersection is a strip wider than the
    // tolerance: its area beyond the tolerance times the smaller facet's
    // diameter. Facets that only share a boundary intersect in a segment.
    double area = 0.0;
    const std::vector<std::array<CcUv, 3>> ff = cc_fan(f), gf = cc_fan(g);
    for (const auto& a : ff)
        for (const auto& b : gf)
            area += cc_triangle_overlap(a, b);
    auto diameter = [](const std::vector<CcUv>& rRing) {
        double d = 0.0;
        for (const CcUv& p : rRing)
            for (const CcUv& q : rRing)
                d = std::max(d, std::hypot(p[0] - q[0], p[1] - q[1]));
        return d;
    };
    return area > Tol * std::min(diameter(f), diameter(g)) ? CcPair::Overlap : CcPair::Apart;
}

CcPair cc_classify_segments(const std::vector<double>& rXyz, const std::int64_t* pF,
                            const std::int64_t* pG, double Tol) {
    const detail::Vec3 a = cc_point(rXyz, pF[0]), b = cc_point(rXyz, pF[1]);
    const detail::Vec3 c = cc_point(rXyz, pG[0]), d = cc_point(rXyz, pG[1]);
    const detail::Vec3 ab = detail::vec3_sub(b, a), cd = detail::vec3_sub(d, c);
    const double lf = detail::vec3_norm(ab), lg = detail::vec3_norm(cd);
    if (!(lf > 0.0) || !(lg > 0.0))
        return CcPair::Apart;
    // Distance of a point from the infinite line through P0 along unit U.
    auto off_line = [&](const detail::Vec3& rP, const detail::Vec3& rP0, const detail::Vec3& rU) {
        const detail::Vec3 w = detail::vec3_sub(rP, rP0);
        return detail::vec3_norm(detail::vec3_cross(w, rU));
    };
    const detail::Vec3 uf = detail::vec3_scale(ab, 1.0 / lf), ug = detail::vec3_scale(cd, 1.0 / lg);
    if (off_line(c, a, uf) > Tol || off_line(d, a, uf) > Tol || off_line(a, c, ug) > Tol ||
        off_line(b, c, ug) > Tol)
        return CcPair::Apart;
    const double tol_sq = Tol * Tol;
    auto near = [&](const detail::Vec3& rX, const detail::Vec3& rY) {
        return detail::vec3_norm_sq(detail::vec3_sub(rX, rY)) <= tol_sq;
    };
    if ((near(a, c) && near(b, d)) || (near(a, d) && near(b, c)))
        return CcPair::Unwelded;
    double t0 = detail::vec3_dot(detail::vec3_sub(c, a), uf);
    double t1 = detail::vec3_dot(detail::vec3_sub(d, a), uf);
    if (t0 > t1)
        std::swap(t0, t1);
    return std::min(lf, t1) - std::max(0.0, t0) > Tol ? CcPair::Overlap : CcPair::Apart;
}

}  // namespace

ConformityResult check_conformity(const Mesh& rMesh, const ConformityOptions& rOptions) {
    if (!(rOptions.mTolerance >= 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kCcPrefix) +
                                    "the tolerance must be finite and non-negative");
    ConformityResult result;

    int dim = -1;
    for (const auto cb : rMesh.CellRange())
        if (cb.NumCells() > 0)
            dim = std::max(dim, cc_block_dimension(cb));
    CcFacets facets;
    if (dim == 3)
        facets = cc_volume_facets(rMesh, result.mWarnings);
    else if (dim == 2)
        facets = cc_surface_facets(rMesh, result.mWarnings);
    else
        result.mWarnings.push_back("the mesh has no 2-D or 3-D cells; nothing was examined");
    result.mDimension = dim < 2 ? 0 : dim;
    result.mNumNonManifoldFacets = facets.mNumNonManifold;
    const std::size_t nf = facets.Size();
    result.mNumBoundaryFacets = static_cast<std::int64_t>(nf);

    const std::size_t npts = rMesh.NumPoints();
    std::vector<std::int64_t> flags(npts, 0);
    if (nf > 0) {
        const std::vector<double> xyz = detail::surface_xyz(rMesh);

        // The facets' boxes and the tolerance, from the box of their corners.
        std::vector<std::array<double, 6>> box(nf);
        const double inf = std::numeric_limits<double>::infinity();
        std::array<double, 3> lo{{inf, inf, inf}};
        std::array<double, 3> hi{{-inf, -inf, -inf}};
        double extent_sum = 0.0;
        for (std::size_t f = 0; f < nf; ++f) {
            std::array<double, 6>& r_box = box[f];
            r_box = {{inf, inf, inf, -inf, -inf, -inf}};
            const std::int64_t* ring = facets.Ring(f);
            for (std::size_t k = 0; k < facets.Corners(f); ++k) {
                const detail::Vec3 p = cc_point(xyz, ring[k]);
                for (std::size_t d = 0; d < 3; ++d) {
                    r_box[d] = std::min(r_box[d], p[d]);
                    r_box[3 + d] = std::max(r_box[3 + d], p[d]);
                }
            }
            for (std::size_t d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], r_box[d]);
                hi[d] = std::max(hi[d], r_box[3 + d]);
            }
            extent_sum += std::max({r_box[3] - r_box[0], r_box[4] - r_box[1], r_box[5] - r_box[2]});
        }
        const double diag =
            std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                      (hi[2] - lo[2]) * (hi[2] - lo[2]));
        const double tol = rOptions.mTolerance * diag;
        for (std::array<double, 6>& r_box : box)
            for (std::size_t d = 0; d < 3; ++d) {
                r_box[d] -= tol;
                r_box[3 + d] += tol;
            }

        // A bucket grid over the facets' boxes, sized by their mean extent and
        // never finer than 1/256 of the diagonal, so a large facet among small
        // ones spans a bounded number of buckets.
        double cell = std::max(extent_sum / static_cast<double>(nf), diag / 256.0);
        if (!(cell > 0.0))
            cell = 1.0;
        detail::SpatialGrid grid(cell);
        for (std::size_t f = 0; f < nf; ++f)
            grid.InsertBox(grid.KeyOf(box[f].data()), grid.KeyOf(box[f].data() + 3),
                           static_cast<std::int64_t>(f));

        // One view per block, and the nodes of a facet's owner on demand.
        const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
        auto owner_has = [&](std::size_t Facet, std::int64_t Point,
                             std::vector<std::int64_t>& rScratch) {
            const auto br = detail::global_to_block_row(bases, facets.mOwner[Facet]);
            cc_cell_nodes(rMesh.Cells(br.first), static_cast<std::size_t>(br.second), rScratch);
            return std::find(rScratch.begin(), rScratch.end(), Point) != rScratch.end();
        };

        // Candidates: the nodes of the cells that own an unmatched facet.
        std::vector<std::int64_t> owners(facets.mOwner);
        std::sort(owners.begin(), owners.end());
        owners.erase(std::unique(owners.begin(), owners.end()), owners.end());
        std::vector<std::int64_t> candidates;
        {
            std::vector<std::int64_t> nodes;
            for (std::int64_t g : owners) {
                const auto br = detail::global_to_block_row(bases, g);
                cc_cell_nodes(rMesh.Cells(br.first), static_cast<std::size_t>(br.second), nodes);
                for (std::int64_t p : nodes)
                    if (p >= 0 && static_cast<std::size_t>(p) < npts)
                        candidates.push_back(p);
            }
        }
        std::sort(candidates.begin(), candidates.end());
        candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());

        const std::int64_t flag = static_cast<std::int64_t>(dim == 3 ? ConformityFlag::HangingNode
                                                                     : ConformityFlag::TJunction);
        const double tol_sq = tol * tol;
        parallel_for(candidates.size(), [&](std::size_t i) {
            const std::int64_t p = candidates[i];
            const detail::Vec3 x = cc_point(xyz, p);
            const std::vector<std::int64_t>* bucket = grid.Find(grid.KeyOf(x.data()));
            if (bucket == nullptr)
                return;
            std::vector<std::int64_t> scratch;
            for (std::int64_t id : *bucket) {
                const std::size_t f = static_cast<std::size_t>(id);
                const std::array<double, 6>& b = box[f];
                if (x[0] < b[0] || x[1] < b[1] || x[2] < b[2] || x[0] > b[3] || x[1] > b[4] ||
                    x[2] > b[5])
                    continue;
                const std::int64_t* ring = facets.Ring(f);
                const std::size_t nc = facets.Corners(f);
                if (std::find(ring, ring + nc, p) != ring + nc)
                    continue;
                if (cc_facet_distance_sq(x, xyz, ring, nc) > tol_sq)
                    continue;
                // On a corner's location: an unwelded duplicate, not a hanging node.
                bool on_corner = false;
                for (std::size_t k = 0; k < nc && !on_corner; ++k)
                    on_corner =
                        detail::vec3_norm_sq(detail::vec3_sub(x, cc_point(xyz, ring[k]))) <= tol_sq;
                if (on_corner || owner_has(f, p, scratch))
                    continue;
                flags[static_cast<std::size_t>(p)] = flag;
                return;
            }
        });

        // Facet pairs, each counted once at its lower facet.
        std::vector<CcPlane> planes(nf);
        if (dim == 3)
            parallel_for(nf, [&](std::size_t f) {
                planes[f] = cc_plane(xyz, facets.Ring(f), facets.Corners(f));
            });
        std::vector<std::int64_t> partial(nf, 0), unwelded(nf, 0);
        parallel_for(nf, [&](std::size_t f) {
            std::vector<std::int64_t> near;
            grid.ForEachInBox(grid.KeyOf(box[f].data()), grid.KeyOf(box[f].data() + 3),
                              [&](const std::vector<std::int64_t>& rIds) {
                                  for (std::int64_t id : rIds)
                                      if (static_cast<std::size_t>(id) > f)
                                          near.push_back(id);
                              });
            std::sort(near.begin(), near.end());
            near.erase(std::unique(near.begin(), near.end()), near.end());
            const std::array<double, 6>& bf = box[f];
            for (std::int64_t id : near) {
                const std::size_t g = static_cast<std::size_t>(id);
                const std::array<double, 6>& bg = box[g];
                if (bf[0] > bg[3] || bg[0] > bf[3] || bf[1] > bg[4] || bg[1] > bf[4] ||
                    bf[2] > bg[5] || bg[2] > bf[5])
                    continue;
                const CcPair kind =
                    dim == 3
                        ? cc_classify_polygons(xyz, facets.Ring(f), facets.Corners(f), planes[f],
                                               facets.Ring(g), facets.Corners(g), planes[g], tol)
                        : cc_classify_segments(xyz, facets.Ring(f), facets.Ring(g), tol);
                if (kind == CcPair::Overlap)
                    ++partial[f];
                else if (kind == CcPair::Unwelded)
                    ++unwelded[f];
            }
        });
        for (std::size_t f = 0; f < nf; ++f) {
            result.mNumPartialOverlaps += partial[f];
            result.mNumUnweldedFacets += unwelded[f];
        }
    }

    for (std::int64_t v : flags) {
        result.mNumHangingNodes += v == static_cast<std::int64_t>(ConformityFlag::HangingNode);
        result.mNumTJunctions += v == static_cast<std::int64_t>(ConformityFlag::TJunction);
    }
    result.mConforming = result.mNumHangingNodes == 0 && result.mNumTJunctions == 0 &&
                         result.mNumPartialOverlaps == 0 && result.mNumUnweldedFacets == 0;
    if (result.mNumPartialOverlaps > 0 || result.mNumHangingNodes > 0 || result.mNumTJunctions > 0)
        result.mInterfaceMode = "proximity";
    else if (result.mNumUnweldedFacets > 0)
        result.mInterfaceMode = "weld";

    result.mMesh = detail::clone_mesh(
        rMesh, [&](DataLocation Location, const std::string& rName, std::string&) {
            return !(rOptions.mAttachFlags && Location == DataLocation::Point &&
                     rName == kConformityHangingName);
        });
    if (rOptions.mAttachFlags) {
        NDArray out = NDArray::Uninit(DType::Int64, {npts});
        std::copy(flags.begin(), flags.end(), out.As<std::int64_t>());
        result.mMesh.AddPointData(kConformityHangingName, std::move(out));
    }
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        result.mMesh.AddRegion(rMesh.Region(i));
    return result;
}

}  // namespace meshioplusplus
