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
// Isosurface stuffing (BCC lattice, warp, cut). See operations/remesh_volume.hpp
// for the full contract and the algorithm's design reasoning. Anonymous-
// namespace helpers are prefixed `rvol_`.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/remesh_volume.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/slot_runs.hpp"
#include "../detail/surface_edge_runs.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

constexpr const char* kRvolPrefix = "meshio++: remesh_volume: ";

// --- exact-position weld key -------------------------------------------------
//
// Multiple crossing edges incident to the SAME warped lattice vertex each get
// their OWN fresh point (Resolve()'s own design, see its doc comment) that
// nonetheless lands at the IDENTICAL warped-target position -- confirmed by
// RemeshVolume.OutputIsWatertightAndPositivelyOriented, whose non-manifold
// edges traced back to exactly this. Welding those AFTER every tet's winding
// is already fixed (remesh_volume's own weld pass, below) is what makes the
// fix safe: merging two ids changes only which id a face references, never
// any tet's corner order, so it cannot reopen the winding-consistency problem
// a vertex-level SHARED id (tried and reverted -- see git history) did. Keyed
// on the exact bit pattern (the `SurfaceEdgeKeyHash` mixing constant and
// shape), not a tolerance: these positions are not merely close, they are the
// SAME `mFinalPosition[...]` value read twice, so exact equality is both
// sufficient and the more conservative choice (a tolerance could weld two
// genuinely distinct nearby crossings that happen to be close by coincidence).
struct RvolPosKey {
    std::int64_t mBits[3];
    bool operator==(const RvolPosKey& rOther) const {
        return mBits[0] == rOther.mBits[0] && mBits[1] == rOther.mBits[1] &&
               mBits[2] == rOther.mBits[2];
    }
};
struct RvolPosKeyHash {
    std::size_t operator()(const RvolPosKey& rKey) const {
        std::size_t h = 0;
        for (std::int64_t b : rKey.mBits)
            h ^= std::hash<std::int64_t>{}(b) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};
RvolPosKey rvol_pos_key(const Vec3& rP) {
    RvolPosKey key;
    static_assert(sizeof(double) == sizeof(std::int64_t), "double/int64 size mismatch");
    for (int i = 0; i < 3; ++i) {
        double v = rP[static_cast<std::size_t>(i)];
        std::memcpy(&key.mBits[i], &v, sizeof(v));
    }
    return key;
}

// --- the BCC 12-tet local table ---------------------------------------------
//
// Local corner indices 0-7 are the standard hexahedron corner layout
// (0=(0,0,0), 1=(1,0,0), 2=(1,1,0), 3=(0,1,0), 4=(0,0,1), 5=(1,0,1),
// 6=(1,1,1), 7=(0,1,1)); index 8 is the cell's own body center (0.5,0.5,0.5).
// Each of the cell's 6 faces is split into 2 triangles by ONE fixed diagonal
// (the one connecting its two lattice-index-lowest and lattice-index-highest
// corners, in the face's own local (u,v) parametrization -- a rule stated
// purely in terms of the FACE's own global coordinates, so two cells sharing
// a face independently compute the identical diagonal with no adjacency
// bookkeeping), and each resulting triangle is connected to the center.
//
// DERIVED and algebraically pinned (never transcribed from anywhere): the 12
// entries tile the unit cell exactly (sum of |volume| == 1, every point of
// the cell belongs to exactly one tet -- RemeshVolume.BccTableTilesTheUnitCellWithNoGapOrOverlap),
// every entry is positively oriented (RemeshVolume.BccTableIsPositivelyOriented),
// and every tet's 6 dihedral angles are exactly {45, 60, 60, 90, 120}
// degrees, independent of cell size (RemeshVolume.BccTableDihedralAnglesAreFixed)
// -- NOT Labelle & Shewchuk's own reported {60, 90, 120} set, since this is
// an independent derivation of the same general idea rather than a
// transcription of their (unread) construction; see remesh_volume.hpp's file
// doc comment.
constexpr std::array<std::array<int, 3>, 12> kBccLocalTets = {{
    {0, 7, 3},
    {0, 4, 7},
    {1, 2, 6},
    {1, 6, 5},
    {0, 1, 5},
    {0, 5, 4},
    {3, 6, 2},
    {3, 7, 6},
    {0, 2, 1},
    {0, 3, 2},
    {4, 5, 6},
    {4, 6, 7},
}};

// --- lattice generation ------------------------------------------------------

// The root BCC lattice: every point (A-lattice corners, then B-lattice
// centers), and every root tet as 4 point ids, already positively oriented
// (kBccLocalTets's own guarantee, since every cell is a translate of the same
// unit construction).
struct RvolLattice {
    std::vector<Vec3> mPoints;
    std::vector<std::array<std::int64_t, 4>> mTets;
    std::int64_t mNumLatticeIds() const { return static_cast<std::int64_t>(mPoints.size()); }
};

RvolLattice rvol_build_lattice(const detail::LatticeSpec& rSpec) {
    RvolLattice lat;
    const std::int64_t nx = rSpec.mDims[0];
    const std::int64_t ny = rSpec.mDims[1];
    const std::int64_t nz = rSpec.mDims[2];
    if (nx <= 0 || ny <= 0 || nz <= 0)
        return lat;

    const std::array<double, 3>& o = rSpec.mOrigin;
    const std::array<double, 3>& h = rSpec.mSpacing;

    const std::int64_t nax = nx + 1, nay = ny + 1, naz = nz + 1;
    const std::int64_t numA = nax * nay * naz;
    const std::int64_t numB = nx * ny * nz;

    auto a_id = [&](std::int64_t i, std::int64_t j, std::int64_t k) {
        return (i * nay + j) * naz + k;
    };
    auto b_id = [&](std::int64_t i, std::int64_t j, std::int64_t k) {
        return numA + (i * ny + j) * nz + k;
    };

    // Every point and tet has a closed-form index, so both fill in parallel:
    // A point (i, j, k) at a_id, B point at b_id, and cell (i, j, k)'s twelve
    // tets at 12 * ((i * ny + j) * nz + k), in kBccLocalTets order.
    lat.mPoints.resize(static_cast<std::size_t>(numA + numB));
    parallel_for(static_cast<std::size_t>(numA), [&](std::size_t id) {
        const std::int64_t a = static_cast<std::int64_t>(id);
        const std::int64_t k = a % naz;
        const std::int64_t j = (a / naz) % nay;
        const std::int64_t i = a / (naz * nay);
        lat.mPoints[id] =
            Vec3{o[0] + static_cast<double>(i) * h[0], o[1] + static_cast<double>(j) * h[1],
                 o[2] + static_cast<double>(k) * h[2]};
    });
    parallel_for(static_cast<std::size_t>(numB), [&](std::size_t id) {
        const std::int64_t b = static_cast<std::int64_t>(id);
        const std::int64_t k = b % nz;
        const std::int64_t j = (b / nz) % ny;
        const std::int64_t i = b / (nz * ny);
        lat.mPoints[static_cast<std::size_t>(numA) + id] =
            Vec3{o[0] + (static_cast<double>(i) + 0.5) * h[0],
                 o[1] + (static_cast<double>(j) + 0.5) * h[1],
                 o[2] + (static_cast<double>(k) + 0.5) * h[2]};
    });

    lat.mTets.resize(static_cast<std::size_t>(numB) * 12);
    parallel_for(static_cast<std::size_t>(numB), [&](std::size_t id) {
        const std::int64_t b = static_cast<std::int64_t>(id);
        const std::int64_t k = b % nz;
        const std::int64_t j = (b / nz) % ny;
        const std::int64_t i = b / (nz * ny);
        const std::int64_t corner[8] = {a_id(i, j, k),
                                        a_id(i + 1, j, k),
                                        a_id(i + 1, j + 1, k),
                                        a_id(i, j + 1, k),
                                        a_id(i, j, k + 1),
                                        a_id(i + 1, j, k + 1),
                                        a_id(i + 1, j + 1, k + 1),
                                        a_id(i, j + 1, k + 1)};
        const std::int64_t centre = b_id(i, j, k);
        std::size_t out = id * 12;
        for (const auto& t : kBccLocalTets)
            lat.mTets[out++] = {centre, corner[t[0]], corner[t[1]], corner[t[2]]};
    });
    return lat;
}

// --- classify + warp ---------------------------------------------------------

// Every lattice vertex's inside/outside LABEL (from its ORIGINAL signed
// distance, fixed for the rest of the algorithm), its raw signed distance
// (needed by a fresh cut point's interpolation parameter), whether it was
// warped, and its FINAL (possibly-warped) position.
struct RvolClassification {
    std::vector<std::uint8_t> mInside;
    std::vector<std::uint8_t> mWarped;
    std::vector<double> mDistance;
    std::vector<Vec3> mFinalPosition;
};

RvolClassification rvol_classify_and_warp(const RvolLattice& rLattice,
                                          const detail::DistanceQuery& rQuery,
                                          const SurfaceDistanceOptions& rDistOpts,
                                          double WarpThreshold, std::int64_t& rNumWarped) {
    const std::size_t n = rLattice.mPoints.size();
    RvolClassification c;
    c.mInside.resize(n);
    c.mWarped.assign(n, 0);
    c.mDistance.resize(n);
    c.mFinalPosition = rLattice.mPoints;
    rNumWarped = 0;
    if (n == 0)
        return c;

    const std::vector<detail::DistanceHit> hits =
        detail::query_distances(rQuery, rLattice.mPoints, rDistOpts);
    parallel_for_bw(n, [&](std::size_t i) {
        c.mDistance[i] = hits[i].mSignedDistance;
        c.mInside[i] = hits[i].mSignedDistance <= 0.0 ? 1 : 0;
    });

    if (!(WarpThreshold > 0.0))
        return c;

    // The candidates in ascending id: flags, then a scan.
    std::vector<std::uint64_t> near(n);
    parallel_for_bw(
        n, [&](std::size_t i) { near[i] = std::fabs(c.mDistance[i]) <= WarpThreshold ? 1 : 0; });
    std::vector<std::uint64_t> at(n);
    const std::uint64_t ncand =
        parallel_exclusive_scan(near.data(), n, at.data(), std::uint64_t{0});
    if (ncand == 0)
        return c;
    std::vector<std::size_t> candidates(ncand);
    parallel_for_bw(n, [&](std::size_t i) {
        if (near[i])
            candidates[at[i]] = i;
    });

    std::vector<Vec3> qpts(candidates.size());
    parallel_for_bw(candidates.size(),
                    [&](std::size_t ci) { qpts[ci] = rLattice.mPoints[candidates[ci]]; });
    const std::vector<detail::ClosestPointHit> chits = detail::query_closest_points(rQuery, qpts);

    std::int64_t warped = 0;
    for (std::size_t ci = 0; ci < candidates.size(); ++ci) {
        if (!chits[ci].mFound)
            continue;
        const std::size_t i = candidates[ci];
        c.mWarped[i] = 1;
        c.mFinalPosition[i] = chits[ci].mPoint;
        ++warped;
    }
    rNumWarped = warped;
    return c;
}

// --- cut: resolve one crossing edge to a point id ----------------------------

// Every crossing edge (one endpoint inside, one outside -- the caller's
// contract, never checked here) gets exactly one FRESH point, in the combined
// id space: [0, numLatticeIds) are original lattice vertices, ids at or beyond
// numLatticeIds fresh points, one per distinct crossing edge, deduped by the
// sorted endpoint-id pair -- the plain marching-tetrahedra rule, which is what
// makes watertightness a structural consequence of global, edge-keyed dedup
// rather than something this code has to reason about per case. Warping
// affects only WHERE that point sits (reusing an incident endpoint's own warp
// target when available), never WHETHER it is a distinct point --
// deliberately a genuinely FRESH point every time, even when it coincides
// exactly with another edge's own fresh point (both incident to the same
// warped vertex): any scheme that instead SHARES one id across multiple edges
// was tried and rejected here (see git history) -- sharing by the warped
// vertex's own id collapsed a straddling tet's kept corner onto its own cut
// point when that same vertex was both kept and a cut-point source; sharing by
// a SEPARATE id per warped vertex avoided that but broke winding CONSISTENCY
// between neighbouring tets instead
// (RemeshVolume.OutputIsWatertightAndPositivelyOriented caught both, first as
// non-manifold edges, then as inconsistently wound pairs).
// Genuinely-coincident fresh points are instead welded AFTER every tet is
// built (see remesh_volume's own weld pass below), which touches only WHICH id
// a face uses, never any tet's already-decided winding.
//
// Fresh ids are numbered in the order a serial sweep over the root tets, in
// ascending order and each tet's own Resolve order, first meets each edge
// (first-seen numbering over the sort-based table, detail/slot_runs.hpp);
// each root tet is cut independently, in parallel, against a tet-local
// resolver whose placeholders are mapped to those ids afterwards.

// Where a crossing edge's fresh point sits: a function of the edge alone.
Vec3 rvol_fresh_position(const RvolLattice& rLattice, const RvolClassification& rClass,
                         std::int64_t Inside, std::int64_t Outside) {
    if (rClass.mWarped[static_cast<std::size_t>(Inside)])
        return rClass.mFinalPosition[static_cast<std::size_t>(Inside)];
    if (rClass.mWarped[static_cast<std::size_t>(Outside)])
        return rClass.mFinalPosition[static_cast<std::size_t>(Outside)];
    const double du = rClass.mDistance[static_cast<std::size_t>(Inside)];
    const double dv = rClass.mDistance[static_cast<std::size_t>(Outside)];
    double t = (du - dv) != 0.0 ? du / (du - dv) : 0.5;
    if (!(t > 0.0 && t < 1.0))
        t = 0.5;  // a near-tangent crossing: the safe interior fallback
    const Vec3& pu = rLattice.mPoints[static_cast<std::size_t>(Inside)];
    const Vec3& pv = rLattice.mPoints[static_cast<std::size_t>(Outside)];
    return {pu[0] + t * (pv[0] - pu[0]), pu[1] + t * (pv[1] - pu[1]), pu[2] + t * (pv[2] - pu[2])};
}

// One root tet's resolver: its crossing edges in first-Resolve order, as
// placeholders numLatticeIds + k (a tet cuts at most four edges).
class RvolLocalResolver {
public:
    RvolLocalResolver(const RvolLattice& rLattice, const RvolClassification& rClass)
        : mrLattice(rLattice), mrClass(rClass) {}

    void Reset() { mCount = 0; }

    std::int64_t Resolve(std::int64_t Inside, std::int64_t Outside) {
        const std::int64_t lo = Inside < Outside ? Inside : Outside;
        const std::int64_t hi = Inside < Outside ? Outside : Inside;
        const detail::SurfaceEdgeKey key{lo, hi};
        for (std::size_t k = 0; k < mCount; ++k)
            if (mKeys[k] == key)
                return mrLattice.mNumLatticeIds() + static_cast<std::int64_t>(k);
        mKeys[mCount] = key;
        mPositions[mCount] = rvol_fresh_position(mrLattice, mrClass, Inside, Outside);
        return mrLattice.mNumLatticeIds() + static_cast<std::int64_t>(mCount++);
    }

    const Vec3& PositionOf(std::int64_t Id) const {
        return Id < mrLattice.mNumLatticeIds()
                   ? mrClass.mFinalPosition[static_cast<std::size_t>(Id)]
                   : mPositions[static_cast<std::size_t>(Id - mrLattice.mNumLatticeIds())];
    }

    std::size_t NumKeys() const { return mCount; }
    const detail::SurfaceEdgeKey& Key(std::size_t K) const { return mKeys[K]; }

private:
    const RvolLattice& mrLattice;
    const RvolClassification& mrClass;
    std::array<detail::SurfaceEdgeKey, 4> mKeys{};
    std::array<Vec3, 4> mPositions{};
    std::size_t mCount = 0;
};

// --- cut: the sign-mask case table -------------------------------------------
//
// Every sub-case's TOPOLOGY (which corners survive, which edges are cut, how
// the cut region decomposes into tets) is derived and numerically pinned --
// tested by comparing the emitted sub-tets' total signed volume against an
// independently-computed expected value (RemeshVolume.CutCaseTopologyIsExact)
// -- never transcribed. Two facts are load-bearing and both are pinned by a
// dedicated test rather than trusted from the derivation alone:
//
// 1-inside / 3-inside: with `k` the single distinguished corner's index into
// the root tet's own 4-entry array and "others" the remaining 3 in their
// ORIGINAL array order, that natural order gives a POSITIVELY oriented result
// only for EVEN k (0 or 2); ODD k (1 or 3) needs the last two of "others"
// swapped. This is a parity fact about the tet corner numbering (a cyclic
// rotation by one position is an odd permutation of 4 elements), not an
// approximation -- RemeshVolume.CutCaseTableIsCorrectlyOrientedForEveryVertexPosition
// checks all 4 positions explicitly.
//
// 2-inside: with the inside pair's ascending indices (i0, i1) and the outside
// pair's ascending indices (j0, j1), the outside pair's order must be
// SWAPPED exactly when (i0, i1) is (0, 2) or (1, 3) -- the tet's two
// "diagonal" (non-adjacent-in-array) vertex pairs -- and used as-is for the
// other 4 of the 6 possible pairings. Same test covers all 6.
// Rejects a degenerate (near-zero or negative) tet -- possible only from an
// ill-conditioned cut near a tangency, never from a whole (count==0/4) tet,
// whose sign is already guaranteed by kBccLocalTets.
bool rvol_tet_ok(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3, double Eps) {
    const Vec3 d1{p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    const Vec3 d2{p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    const Vec3 d3{p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]};
    const double v = d1[0] * (d2[1] * d3[2] - d2[2] * d3[1]) -
                     d1[1] * (d2[0] * d3[2] - d2[2] * d3[0]) +
                     d1[2] * (d2[0] * d3[1] - d2[1] * d3[0]);
    return v > Eps;
}

template <class TResolver, class TOut>
void rvol_cut_tet(const std::array<std::int64_t, 4>& rV, const RvolClassification& rClass,
                  TResolver& rResolver, TOut& rOut) {
    bool lab[4];
    int count = 0;
    for (int i = 0; i < 4; ++i) {
        lab[i] = rClass.mInside[static_cast<std::size_t>(rV[static_cast<std::size_t>(i)])] != 0;
        if (lab[i])
            ++count;
    }
    if (count == 0)
        return;
    if (count == 4) {
        rOut.push_back(rV);
        return;
    }

    // Checks each candidate tet's sign and swaps its own last two corners if
    // negative -- a purely WITHIN-tet reordering that cannot change which
    // vertices any face uses, so it cannot reopen a cross-tet consistency
    // problem. Declared once here so both the count==3 and count==2 branches
    // below share it.
    auto push_oriented = [&](std::int64_t p0, std::int64_t p1, std::int64_t p2, std::int64_t p3) {
        if (rvol_tet_ok(rResolver.PositionOf(p0), rResolver.PositionOf(p1),
                        rResolver.PositionOf(p2), rResolver.PositionOf(p3), 0.0))
            rOut.push_back({p0, p1, p2, p3});
        else
            rOut.push_back({p0, p1, p3, p2});
    };

    if (count == 1) {
        // The single inside vertex is unambiguous (there is only one
        // candidate), so unlike count==2/3 below, no GLOBAL-vertex-id rule is
        // needed here: which 3 outside neighbours it pairs with -- and hence
        // which 3 faces get emitted -- is fixed by the LABELS alone, and
        // their internal order affects only this one tet's own winding
        // (fixed by push_oriented), never which faces exist.
        int k = -1;
        for (int i = 0; i < 4; ++i)
            if (lab[i]) {
                k = i;
                break;
            }
        const std::int64_t pk = rV[static_cast<std::size_t>(k)];
        const std::int64_t pb = rV[static_cast<std::size_t>((k + 1) % 4)];
        const std::int64_t pc = rV[static_cast<std::size_t>((k + 2) % 4)];
        const std::int64_t pd = rV[static_cast<std::size_t>((k + 3) % 4)];
        const std::int64_t e1 = rResolver.Resolve(pk, pb);
        const std::int64_t e2 = rResolver.Resolve(pk, pc);
        const std::int64_t e3 = rResolver.Resolve(pk, pd);
        push_oriented(pk, e1, e2, e3);
        return;
    }

    if (count == 3) {
        // The single outside vertex pk is unambiguous, but the 3 INSIDE
        // vertices (b, c, d) are not: which one plays the "far corner" that
        // recurs in every sub-tet is exactly the same kind of ambiguity
        // count==2's hub choice has below, and for the identical reason
        // (a neighbouring tet sharing a 2-inside face sees the same two
        // global ids at different local slots) must be resolved by GLOBAL
        // vertex id, not local array position. Sorted ascending by global id
        // (B0, B1, B2), the "inside" region is a triangular prism with ends
        // (B0, B1, B2) and (cut(B0,pk), cut(B1,pk), cut(B2,pk)) -- the exact
        // topological shape count==2's own prism has (2 triangular ends, 3
        // quad sides), split by the identical canonical "staircase" formula,
        // confirmed to agree with count==2's own diagonal choice on every
        // shared 2-inside-type face (RemeshVolume.OutputIsWatertightAndPositivelyOriented
        // is what actually caught this cross-case disagreement -- a
        // same-case self-consistency check alone, which an earlier version
        // of this branch already passed, could not).
        std::int64_t inside_ids[3];
        std::int64_t pk = -1;
        int ii = 0;
        for (int i = 0; i < 4; ++i)
            if (lab[i])
                inside_ids[ii++] = rV[static_cast<std::size_t>(i)];
            else
                pk = rV[static_cast<std::size_t>(i)];
        std::sort(inside_ids, inside_ids + 3);
        const std::int64_t b0 = inside_ids[0], b1 = inside_ids[1], b2 = inside_ids[2];
        const std::int64_t e0 = rResolver.Resolve(b0, pk);
        const std::int64_t e1 = rResolver.Resolve(b1, pk);
        const std::int64_t e2 = rResolver.Resolve(b2, pk);
        push_oriented(b0, b1, b2, e2);
        push_oriented(b0, b1, e2, e1);
        push_oriented(b0, e1, e2, e0);
        return;
    }

    // count == 2. Which of the 2 inside vertices becomes the "hub" (the one
    // appearing in all 3 sub-tets) is genuinely ambiguous from this tet's
    // own local vertex order alone -- and MUST be resolved by a rule that
    // depends only on the two GLOBAL vertex ids, never on which of the four
    // local array slots (0-3) they happen to occupy in THIS root tet.
    // A neighbouring lattice tet sharing this same face sees the identical
    // two global ids but very likely at DIFFERENT local slots (each root
    // tet's own 4-vertex order comes from kBccLocalTets independently), so a
    // local-index-based hub choice picks OPPOSITE diagonals for the two
    // tets' shared quad face -- caught by
    // RemeshVolume.OutputIsWatertightAndPositivelyOriented's non-manifold-edge
    // count staying stubbornly nonzero even after the face-sharing bug above
    // was fixed. Smaller global id = hub ("a"), for both the inside and the
    // outside pair -- an arbitrary but, critically, GLOBALLY CONSISTENT rule.
    std::int64_t inside_ids[2], outside_ids[2];
    int ii = 0, jj = 0;
    for (int i = 0; i < 4; ++i)
        if (lab[i])
            inside_ids[ii++] = rV[static_cast<std::size_t>(i)];
        else
            outside_ids[jj++] = rV[static_cast<std::size_t>(i)];
    if (inside_ids[0] > inside_ids[1])
        std::swap(inside_ids[0], inside_ids[1]);
    if (outside_ids[0] > outside_ids[1])
        std::swap(outside_ids[0], outside_ids[1]);
    const std::int64_t pa = inside_ids[0], pb = inside_ids[1];
    const std::int64_t pc = outside_ids[0], pd = outside_ids[1];
    const std::int64_t eac = rResolver.Resolve(pa, pc);
    const std::int64_t ead = rResolver.Resolve(pa, pd);
    const std::int64_t ebc = rResolver.Resolve(pb, pc);
    const std::int64_t ebd = rResolver.Resolve(pb, pd);
    // The standard "staircase fan" split of a triangular prism whose two
    // triangular ends are (pa, eac, ead) and (pb, ebc, ebd) (corresponding
    // vertex-for-vertex: pa<->pb, eac<->ebc, ead<->ebd) -- the identical
    // canonical formula count==3 above uses, verified against an
    // independently-computed convex-hull volume AND checked for no face
    // shared by more than 2 of the 3 tets (a stronger, and as it turned out
    // necessary, check than volume-sum alone; see the git history for the
    // face-sharing bug this replaced). Global-id ordering fixes WHICH points
    // are used (for cross-tet consistency); push_oriented (declared above)
    // fixes each candidate's winding independently.
    push_oriented(pa, eac, ead, ebd);
    push_oriented(pa, eac, ebd, ebc);
    push_oriented(pa, ebc, ebd, pb);
}

// Does the mesh contain a 3D or polyhedron block? Selects whether the input is
// treated as a volume (its boundary is extracted) or a surface (used as-is).
bool rvol_has_volume_cells(const Mesh& rMesh) {
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            return true;
        if (!cb.IsRagged() && cell_type_dimension(cell_type_from_name(cb.Type())) == 3)
            return true;
    }
    return false;
}

}  // namespace

RemeshVolumeResult remesh_volume(const Mesh& rMesh, const RemeshVolumeOptions& rOptions) {
    if (!(rOptions.mWarpFraction >= 0.0))
        throw std::invalid_argument(std::string(kRvolPrefix) +
                                    "warp_fraction must not be "
                                    "negative, got " +
                                    std::to_string(rOptions.mWarpFraction));

    const bool is_volume = rvol_has_volume_cells(rMesh);
    const Mesh surface_owner = is_volume ? extract_surface(rMesh) : Mesh{};
    const Mesh& surface = is_volume ? surface_owner : rMesh;

    const detail::TriangleSoup soup =
        detail::build_triangle_soup(surface, rOptions.mDistance.mSurfaceRegion);

    detail::LatticeRequest req;
    req.mResolution = rOptions.mResolution;
    req.mCellSize = rOptions.mCellSize;
    req.mBounds = rOptions.mBounds;
    req.mPadding = rOptions.mPadding;
    req.mPaddingRelative = rOptions.mPaddingRelative;
    req.mMaxCells = rOptions.mMaxCells;
    const detail::LatticeSpec spec = detail::lattice_resolve(surface, req, kRvolPrefix);

    // The BCC construction assumes a genuinely cubic cell; a per-axis
    // mResolution over a non-cube bounding box can resolve to a rectangular
    // one, which would silently distort every dihedral-angle guarantee this
    // operation makes. Named rather than allowed through.
    if (spec.mDims[0] > 0 && spec.mDims[1] > 0 && spec.mDims[2] > 0) {
        const double h0 = spec.mSpacing[0];
        const double tol = 1e-9 * (h0 > 0.0 ? h0 : 1.0);
        if (std::fabs(spec.mSpacing[1] - h0) > tol || std::fabs(spec.mSpacing[2] - h0) > tol)
            throw std::invalid_argument(
                std::string(kRvolPrefix) + "the resolved lattice cell is not cubic (spacing " +
                std::to_string(spec.mSpacing[0]) + ", " + std::to_string(spec.mSpacing[1]) + ", " +
                std::to_string(spec.mSpacing[2]) +
                ") -- a per-axis resolution over a non-cube bounding box does not fit this "
                "operation's fixed BCC cell; use cell_size instead");
    }

    detail::warn_regions_dropped(rMesh, "remesh_volume");

    RemeshVolumeResult result;
    const detail::SurfaceEdgeRuns edge_runs = detail::surface_edge_runs(soup);
    result.mQuality =
        detail::soup_quality(soup, detail::build_surface_edges_from_runs(soup, edge_runs));
    if (rOptions.mDistance.mWatertightCheck == SdfWatertightCheck::Error &&
        !result.mQuality.mWatertight)
        throw std::invalid_argument(
            std::string(kRvolPrefix) +
            "the surface is not watertight "
            "(boundary_edges=" +
            std::to_string(result.mQuality.mBoundaryEdges) +
            ", non_manifold_edges=" + std::to_string(result.mQuality.mNonManifoldEdges) +
            ", inconsistent_pairs=" + std::to_string(result.mQuality.mInconsistentPairs) + ")");
    else if (rOptions.mDistance.mWatertightCheck == SdfWatertightCheck::Warn &&
             !result.mQuality.mWatertight)
        log::warn(
            "remesh_volume: the surface is not watertight; signs may be unreliable "
            "near the defect");

    const detail::DistanceQuery query =
        detail::build_distance_query_from_runs(soup, rOptions.mDistance, edge_runs);
    const RvolLattice lattice = rvol_build_lattice(spec);

    const double h = spec.mSpacing[0];
    std::int64_t num_warped = 0;
    const RvolClassification cls = rvol_classify_and_warp(lattice, query, rOptions.mDistance,
                                                          rOptions.mWarpFraction * h, num_warped);
    result.mNumVerticesWarped = num_warped;

    // --- cut: every root tet in parallel against a tet-local resolver (at
    // most three sub-tets and four crossing edges each), then the fresh
    // points numbered first-seen over (root tet, Resolve order) -- the ids the
    // serial sweep handed out.
    const std::int64_t num_lattice = lattice.mNumLatticeIds();
    const std::size_t nroot = lattice.mTets.size();
    struct RvolTetCut {
        std::uint8_t mNumSub = 0;
        std::uint8_t mNumKeys = 0;
    };
    struct RvolSubs {
        std::array<std::array<std::int64_t, 4>, 3> mTet{};
        std::size_t mCount = 0;
        void push_back(const std::array<std::int64_t, 4>& rT) { mTet[mCount++] = rT; }
    };
    std::vector<RvolTetCut> per_root(nroot);
    parallel_for(nroot, [&](std::size_t r) {
        RvolLocalResolver local(lattice, cls);
        RvolSubs subs;
        rvol_cut_tet(lattice.mTets[r], cls, local, subs);
        per_root[r].mNumSub = static_cast<std::uint8_t>(subs.mCount);
        per_root[r].mNumKeys = static_cast<std::uint8_t>(local.NumKeys());
    });
    std::vector<std::uint64_t> nsub(nroot), nkey(nroot);
    parallel_for_bw(nroot, [&](std::size_t r) {
        nsub[r] = per_root[r].mNumSub;
        nkey[r] = per_root[r].mNumKeys;
    });
    std::vector<std::uint64_t> sub_at(nroot), key_at(nroot);
    const std::uint64_t total_sub =
        parallel_exclusive_scan(nsub.data(), nroot, sub_at.data(), std::uint64_t{0});
    const std::uint64_t total_key =
        parallel_exclusive_scan(nkey.data(), nroot, key_at.data(), std::uint64_t{0});
    std::vector<std::array<std::int64_t, 4>> sub_tets(total_sub);
    std::vector<std::uint64_t> sub_root(total_sub);
    std::vector<detail::SurfaceEdgeKey> keys(total_key);
    parallel_for(nroot, [&](std::size_t r) {
        if (nsub[r] == 0)
            return;
        RvolLocalResolver local(lattice, cls);
        RvolSubs subs;
        rvol_cut_tet(lattice.mTets[r], cls, local, subs);
        for (std::size_t k = 0; k < local.NumKeys(); ++k)
            keys[key_at[r] + k] = local.Key(k);
        for (std::size_t q = 0; q < subs.mCount; ++q) {
            sub_tets[sub_at[r] + q] = subs.mTet[q];
            sub_root[sub_at[r] + q] = r;
        }
    });
    const detail::SlotRuns key_runs = detail::group_slots(
        keys, static_cast<std::size_t>(num_lattice) + 1, [&](const detail::SurfaceEdgeKey& rK) {
            return rK[0] >= 0 && rK[0] < num_lattice ? static_cast<std::size_t>(rK[0])
                                                     : static_cast<std::size_t>(num_lattice);
        });
    const detail::FirstSeen fresh_ids = detail::number_first_seen(key_runs, total_key);
    // Each fresh point at its edge's position (the inside endpoint is the one
    // labelled inside: every Resolve call names it first).
    std::vector<Vec3> fresh(fresh_ids.NumIds());
    parallel_for(fresh.size(), [&](std::size_t f) {
        const detail::SurfaceEdgeKey& k = keys[key_runs.Head(fresh_ids.mRunOfId[f])];
        const bool lo_inside = cls.mInside[static_cast<std::size_t>(k[0])] != 0;
        fresh[f] =
            rvol_fresh_position(lattice, cls, lo_inside ? k[0] : k[1], lo_inside ? k[1] : k[0]);
    });
    auto point_at = [&](std::int64_t id) -> const Vec3& {
        return id < num_lattice ? cls.mFinalPosition[static_cast<std::size_t>(id)]
                                : fresh[static_cast<std::size_t>(id - num_lattice)];
    };
    // Placeholders to fresh ids, and the degenerate test, per sub-tet.
    const double eps = 1e-15 * h * h * h;
    std::vector<std::uint64_t> ok(total_sub);
    parallel_for(total_sub, [&](std::size_t q) {
        std::array<std::int64_t, 4>& tt = sub_tets[q];
        for (std::int64_t& id : tt)
            if (id >= num_lattice)
                id = num_lattice + fresh_ids.mIdOfSlot[key_at[sub_root[q]] +
                                                       static_cast<std::size_t>(id - num_lattice)];
        ok[q] = rvol_tet_ok(point_at(tt[0]), point_at(tt[1]), point_at(tt[2]), point_at(tt[3]), eps)
                    ? 1
                    : 0;
    });
    // The kept sub-tets in (root tet, sub-tet) order.
    const auto compact = [](const std::vector<std::array<std::int64_t, 4>>& rAll,
                            const std::vector<std::uint64_t>& rKeep) {
        std::vector<std::uint64_t> at(rAll.size());
        const std::uint64_t nkept =
            parallel_exclusive_scan(rKeep.data(), rAll.size(), at.data(), std::uint64_t{0});
        std::vector<std::array<std::int64_t, 4>> kept(nkept);
        parallel_for_bw(rAll.size(), [&](std::size_t q) {
            if (rKeep[q])
                kept[at[q]] = rAll[q];
        });
        return kept;
    };
    std::vector<std::array<std::int64_t, 4>> raw_tets = compact(sub_tets, ok);
    std::int64_t num_rejected = static_cast<std::int64_t>(total_sub - raw_tets.size());

    // --- weld: exact-position dedup of fresh points, then re-check for any
    // tet that welding itself made degenerate (two of its own corners
    // collapsing onto one welded id) -- see RvolPosKey's own doc comment for
    // why this is the safe place to fix duplicate coincident points, and
    // why it is exact rather than tolerance-based. Each fresh point maps to
    // the first fresh point at its exact position: the head of its run.
    {
        std::vector<RvolPosKey> pos_keys(fresh.size());
        parallel_for_bw(fresh.size(), [&](std::size_t i) { pos_keys[i] = rvol_pos_key(fresh[i]); });
        const std::size_t nb = fresh.size() / 4 + 1;
        const detail::SlotRuns runs = detail::group_slots(
            pos_keys, nb, [nb](const RvolPosKey& rK) { return RvolPosKeyHash{}(rK) % nb; },
            [](const RvolPosKey& rA, const RvolPosKey& rB) {
                return std::lexicographical_compare(rA.mBits, rA.mBits + 3, rB.mBits, rB.mBits + 3);
            });
        std::vector<std::int64_t> fresh_remap(fresh.size());
        parallel_for(runs.NumRuns(), [&](std::size_t r) {
            const std::int64_t head = static_cast<std::int64_t>(runs.Head(r));
            for (const std::uint64_t* q = runs.Begin(r); q != runs.End(r); ++q)
                fresh_remap[*q] = head;
        });
        std::vector<std::uint64_t> keep(raw_tets.size());
        parallel_for(raw_tets.size(), [&](std::size_t q) {
            std::array<std::int64_t, 4>& tt = raw_tets[q];
            for (std::int64_t& id : tt)
                if (id >= num_lattice)
                    id = num_lattice + fresh_remap[static_cast<std::size_t>(id - num_lattice)];
            keep[q] =
                rvol_tet_ok(point_at(tt[0]), point_at(tt[1]), point_at(tt[2]), point_at(tt[3]), eps)
                    ? 1
                    : 0;
        });
        std::vector<std::array<std::int64_t, 4>> welded = compact(raw_tets, keep);
        num_rejected += static_cast<std::int64_t>(raw_tets.size() - welded.size());
        raw_tets = std::move(welded);
    }
    result.mNumTetsRejected = num_rejected;

    if (static_cast<std::int64_t>(raw_tets.size()) > rOptions.mMaxTets)
        throw std::invalid_argument(
            std::string(kRvolPrefix) + "the cut output has " + std::to_string(raw_tets.size()) +
            " tets, exceeding max_tets (" + std::to_string(rOptions.mMaxTets) +
            "); coarsen the lattice or raise max_tets");

    // --- compact: used lattice vertices (in ascending original id) + the
    // fresh points (in their first-seen order) ------------------------------
    std::vector<std::uint64_t> used(static_cast<std::size_t>(num_lattice), 0);
    for (const auto& t : raw_tets)
        for (std::int64_t id : t)
            if (id < num_lattice)
                used[static_cast<std::size_t>(id)] = 1;
    std::vector<std::uint64_t> used_at(static_cast<std::size_t>(num_lattice));
    const std::uint64_t num_used = parallel_exclusive_scan(
        used.data(), static_cast<std::size_t>(num_lattice), used_at.data(), std::uint64_t{0});
    const std::int64_t new_base = static_cast<std::int64_t>(num_used);
    auto final_id = [&](std::int64_t id) -> std::int64_t {
        if (id < num_lattice)
            return used[static_cast<std::size_t>(id)]
                       ? static_cast<std::int64_t>(used_at[static_cast<std::size_t>(id)])
                       : -1;
        return new_base + (id - num_lattice);
    };

    result.mNumTets = static_cast<std::int64_t>(raw_tets.size());

    Mesh& out = result.mMesh;
    {
        const std::size_t npts = static_cast<std::size_t>(num_used) + fresh.size();
        NDArray pts = NDArray::Uninit(DType::Float64, {npts, std::size_t{3}});
        double* dst = pts.As<double>();
        parallel_for_bw(static_cast<std::size_t>(num_lattice), [&](std::size_t i) {
            if (!used[i])
                return;
            for (std::size_t d = 0; d < 3; ++d)
                dst[used_at[i] * 3 + d] = cls.mFinalPosition[i][d];
        });
        parallel_for_bw(fresh.size(), [&](std::size_t f) {
            for (std::size_t d = 0; d < 3; ++d)
                dst[(static_cast<std::size_t>(num_used) + f) * 3 + d] = fresh[f][d];
        });
        out.AssignPoints(std::move(pts));
    }
    {
        NDArray conn = NDArray::Uninit(DType::Int64, {raw_tets.size(), std::size_t{4}});
        std::int64_t* dst = conn.As<std::int64_t>();
        parallel_for_bw(raw_tets.size(), [&](std::size_t i) {
            for (std::size_t c = 0; c < 4; ++c)
                dst[i * 4 + c] = final_id(raw_tets[i][c]);
        });
        if (!raw_tets.empty())
            out.AddCellBlock(cell_type_name(CellType::Tetra), std::move(conn));
    }
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, rMesh.FieldData(name));

    // Measure, don't assume: report the OUTPUT mesh's own boundary
    // manifoldness rather than trust the warp/cut construction blindly (see
    // RemeshVolumeResult::mNumNonManifoldEdges' doc comment). Skipped for an
    // empty result, which extract_surface refuses.
    if (!raw_tets.empty())
        result.mNumNonManifoldEdges =
            surface_watertight_check(extract_surface(out)).mNonManifoldEdges;

    return result;
}

}  // namespace meshioplusplus
