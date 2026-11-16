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
// octree_mesh, tetrahedral output: fitting and optimisation after stuffing and
// feature snapping. Each round relocates the interior nodes by an ODT pass
// (`smooth_odt_pass`), relaxes the nodes on the skin along it (the mean of
// their neighbours on the skin, projected back onto the skin, kept only where
// it raises the worst dihedral angle around the node), puts back any node
// whose worst dihedral fell below the mesh's worst before the round, and
// (one material) tries 2-3 and 3-2 flips (`optimize_volume`), kept only when
// they do not lower that worst angle either. Features, the lattice's box and
// the interfaces between materials stay fixed. Anonymous-namespace helpers
// are prefixed `octmo_` (the amalgamation concatenates every translation
// unit into one).

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

// Project includes
#include "../detail/bcc_octree_tets.hpp"
#include "../detail/mesh_from_soup.hpp"
#include "octree_mesh_internal.hpp"
#include "smooth_odt.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/optimize_volume.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr int kOctmoFace[4][3] = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 2, 1}};

// The smallest and the largest dihedral angle of a tetrahedron, from one
// tet_dihedral_angles call (its angles are ascending); -1 and 181 for one
// that is not positive.
std::pair<double, double> octmo_min_max_dihedral(const std::vector<Vec3>& rP,
                                                 const std::array<std::int64_t, 4>& rT) {
    const Vec3& a = rP[static_cast<std::size_t>(rT[0])];
    const Vec3& b = rP[static_cast<std::size_t>(rT[1])];
    const Vec3& c = rP[static_cast<std::size_t>(rT[2])];
    const Vec3& d = rP[static_cast<std::size_t>(rT[3])];
    if (triple_product(vec3_sub(b, a), vec3_sub(c, a), vec3_sub(d, a)) <= 0.0)
        return {-1.0, 181.0};
    const std::array<double, 6> angles = tet_dihedral_angles(a, b, c, d);
    return {angles[0], angles[5]};
}

// The smallest dihedral angle of every tetrahedron into @p rLo and the largest
// into @p rHi, recomputed only where @p pDirty is set (every tetrahedron when it
// is null, or when the arrays do not have the tetrahedra's size); returns the
// minimum of @p rLo and the maximum of @p rHi.
std::pair<double, double> octmo_angles(const std::vector<Vec3>& rP,
                                       const std::vector<std::array<std::int64_t, 4>>& rT,
                                       std::vector<double>& rLo, std::vector<double>& rHi,
                                       const std::vector<std::uint8_t>* pDirty = nullptr) {
    const bool all = !pDirty || rLo.size() != rT.size() || rHi.size() != rT.size();
    rLo.resize(rT.size());
    rHi.resize(rT.size());
    parallel_for(rT.size(), [&](std::size_t t) {
        if (all || (*pDirty)[t])
            std::tie(rLo[t], rHi[t]) = octmo_min_max_dihedral(rP, rT[t]);
    });
    double lo = 180.0;
    for (double a : rLo)
        lo = std::min(lo, a);
    double hi = 0.0;
    for (double a : rHi)
        hi = std::max(hi, a);
    return {lo, hi};
}

std::vector<std::vector<std::int64_t>> octmo_incidence(
    std::size_t Np, const std::vector<std::array<std::int64_t, 4>>& rT) {
    std::vector<std::vector<std::int64_t>> inc(Np);
    for (std::size_t t = 0; t < rT.size(); ++t)
        for (std::int64_t v : rT[t])
            inc[static_cast<std::size_t>(v)].push_back(static_cast<std::int64_t>(t));
    return inc;
}

}  // namespace

void octree_mesh_optimize(OctreeMeshSetup& rS, const OctreeMeshOptions& rO, OctreeMeshResult& rR,
                          OctreeMeshTets& rTets) {
    if (rO.mFitIterations <= 0)
        return;
    // The tetrahedra that are kept, as their own mesh.
    std::vector<std::size_t> kept;
    std::int32_t lo_mat = 1 << 30, hi_mat = -1;
    for (std::size_t t = 0; t < rTets.mTets.size(); ++t)
        if (rTets.mMaterial[t] > 0 || rO.mKeepOutside) {
            kept.push_back(t);
            lo_mat = std::min(lo_mat, rTets.mMaterial[t]);
            hi_mat = std::max(hi_mat, rTets.mMaterial[t]);
        }
    if (kept.empty())
        return;
    const bool one_material = lo_mat == hi_mat;
    std::vector<std::array<std::int64_t, 4>> tets(kept.size());
    for (std::size_t n = 0; n < kept.size(); ++n)
        tets[n] = rTets.mTets[kept[n]];
    std::vector<Vec3>& pts = rTets.mPoints;
    const std::size_t np = pts.size();
    // The faces used once (the mesh's boundary: the skin, and the box with
    // keep_outside), and the skin's surface graph.
    std::vector<std::pair<std::array<std::int64_t, 3>, std::int64_t>> faces;
    for (std::size_t n = 0; n < tets.size(); ++n)
        for (const auto& f : kOctmoFace) {
            std::array<std::int64_t, 3> k{tets[n][f[0]], tets[n][f[1]], tets[n][f[2]]};
            std::sort(k.begin(), k.end());
            faces.emplace_back(k, rTets.mMaterial[kept[n]]);
        }
    std::sort(faces.begin(), faces.end());
    std::vector<std::uint8_t> boundary(np, 0);
    std::vector<std::vector<std::int64_t>> surf_adj(np);
    auto link = [&](const std::array<std::int64_t, 3>& rF) {
        for (int e = 0; e < 3; ++e) {
            surf_adj[static_cast<std::size_t>(rF[e])].push_back(rF[(e + 1) % 3]);
            surf_adj[static_cast<std::size_t>(rF[(e + 1) % 3])].push_back(rF[e]);
        }
    };
    for (std::size_t i = 0; i < faces.size();) {
        std::size_t j = i;
        while (j < faces.size() && faces[j].first == faces[i].first)
            ++j;
        const bool once = j - i == 1;
        const bool interface = j - i == 2 && faces[i].second != faces[i + 1].second;
        if (once)
            for (std::int64_t v : faces[i].first)
                boundary[static_cast<std::size_t>(v)] = 1;
        // The skin: an inside face used once, or one between two materials.
        if ((once && faces[i].second > 0) || interface)
            link(faces[i].first);
        i = j;
    }
    for (auto& a : surf_adj) {
        std::sort(a.begin(), a.end());
        a.erase(std::unique(a.begin(), a.end()), a.end());
    }
    // What never moves: features, the box, the interfaces of several materials,
    // and (for the interior pass) every node on the skin.
    std::vector<std::uint8_t> frozen(np, 0), slide(np, 0);
    for (std::size_t p = 0; p < np; ++p) {
        const bool on_skin = !surf_adj[p].empty();
        const bool feature = !rTets.mFeature.empty() && rTets.mFeature[p];
        frozen[p] = boundary[p] || on_skin || feature || rTets.mOnSkin[p];
        slide[p] = on_skin && !feature && one_material && rTets.mOnSkin[p];
        if (boundary[p] && !on_skin)
            slide[p] = 0;  // the box
    }
    std::vector<std::uint8_t> moved(np, 0);
    SmoothOptions so;
    so.mMethod = SmoothMethod::Odt;
    // The dihedral floor below already puts back every node that inverts a
    // tetrahedron, so the pass's own guard (and its log line per sweep) is
    // only needed without it.
    so.mGuardInversion = !rO.mDihedralFloorGuard;
    so.mLambda = 0.9;
    std::int64_t flips = 0;
    for (std::int64_t round = 0; round < rO.mFitIterations; ++round) {
        std::vector<double> angle, big;
        const auto [floor, ceiling] = octmo_angles(pts, tets, angle, big);
        const std::vector<Vec3> before = pts;
        // 1. ODT relocation of the interior nodes.
        std::vector<double> xyz(np * 3);
        for (std::size_t p = 0; p < np; ++p)
            for (std::size_t k = 0; k < 3; ++k)
                xyz[p * 3 + k] = pts[p][k];
        smooth_odt_pass(xyz, tets, frozen, so);
        for (std::size_t p = 0; p < np; ++p)
            for (std::size_t k = 0; k < 3; ++k)
                pts[p][k] = xyz[p * 3 + k];
        // 2. The skin's nodes slide along it, where that helps.
        const auto inc = octmo_incidence(np, tets);
        std::vector<std::size_t> sliders;
        for (std::size_t p = 0; p < np; ++p)
            if (slide[p])
                sliders.push_back(p);
        std::vector<Vec3> target(sliders.size());
        for (std::size_t k = 0; k < sliders.size(); ++k) {
            Vec3 c{0.0, 0.0, 0.0};
            for (std::int64_t q : surf_adj[sliders[k]])
                c = vec3_add(c, pts[static_cast<std::size_t>(q)]);
            target[k] = vec3_scale(c, 1.0 / static_cast<double>(surf_adj[sliders[k]].size()));
        }
        const auto hits = query_closest_points(rS.mQuery, target);
        for (std::size_t k = 0; k < sliders.size(); ++k) {
            const std::size_t p = sliders[k];
            // The node's worst angle: the smaller of its smallest dihedral
            // and 180 minus its largest.
            auto worst_at = [&]() {
                double w = 180.0;
                for (std::int64_t t : inc[p]) {
                    const auto [lo, hi] =
                        octmo_min_max_dihedral(pts, tets[static_cast<std::size_t>(t)]);
                    w = std::min({w, lo, 180.0 - hi});
                }
                return w;
            };
            const double worst = worst_at();
            const Vec3 keep = pts[p];
            pts[p] = hits[k].mPoint;
            const double after = worst_at();
            if (after <= worst) {
                pts[p] = keep;
                continue;
            }
        }
        // 3. The dihedral floor: put back every node whose worst angle fell
        // below the mesh's worst before the round, or that inverted or flattened
        // a tetrahedron, until none does.
        if (rO.mDihedralFloorGuard) {
            // Clearly positive before the round, inverted or within the zero
            // tolerance now (the dihedral angles do not see the sign).
            // Per tetrahedron, whether it is fine now (1), zero or inverted (0), or not looked at
            // since a node of it last changed (-1): each of its nodes asks within a sweep.
            std::vector<std::int8_t> fine_now(tets.size(), -1);
            const auto spoiled_by_move = [&](const std::size_t T) {
                const std::array<std::int64_t, 4>& rT = tets[T];
                if (fine_now[T] < 0) {
                    const auto [volume_now, tolerance_now] =
                        octree_mesh_volume_and_tolerance(pts, rT, rO.mZeroVolumeFactor);
                    fine_now[T] = volume_now > tolerance_now ? 1 : 0;
                }
                // Most tetrahedra are fine now: the state before the round is only looked at
                // for the few that are not.
                if (fine_now[T] == 1)
                    return false;
                const auto [volume_before, tolerance_before] =
                    octree_mesh_volume_and_tolerance(before, rT, rO.mZeroVolumeFactor);
                return volume_before > tolerance_before;
            };
            // The angles of a sweep: all of them on the first, then only those of the
            // tetrahedra around a node the previous sweep put back (the others' corners
            // did not move, so their angles are what a full recomputation gives).
            std::vector<std::uint8_t> dirty(tets.size(), 1);
            for (int sweep = 0; sweep < 64; ++sweep) {
                octmo_angles(pts, tets, angle, big, sweep == 0 ? nullptr : &dirty);
                std::fill(dirty.begin(), dirty.end(), std::uint8_t{0});
                std::fill(fine_now.begin(), fine_now.end(), std::int8_t{-1});
                bool any = false;
                for (std::size_t p = 0; p < np; ++p) {
                    if (pts[p] == before[p])
                        continue;
                    for (std::int64_t t : inc[p])
                        if (angle[static_cast<std::size_t>(t)] < floor ||
                            big[static_cast<std::size_t>(t)] > ceiling ||
                            spoiled_by_move(static_cast<std::size_t>(t))) {
                            pts[p] = before[p];
                            for (std::int64_t u : inc[p]) {
                                fine_now[static_cast<std::size_t>(u)] = -1;
                                dirty[static_cast<std::size_t>(u)] = 1;
                            }
                            any = true;
                            break;
                        }
                }
                if (!any)
                    break;
                if (sweep == 63)
                    pts = before;  // never left below the floor
            }
        }
        for (std::size_t p = 0; p < np; ++p)
            if (pts[p] != before[p])
                moved[p] = 1;
        // 4. Flips, for one material, kept when the worst angle holds.
        if (one_material) {
            const Mesh m = mesh_from_soup(pts, tets, "tetra");
            OptimizeVolumeOptions oo;
            oo.mMaxIterations = 1;
            oo.mRelocate = false;
            oo.mFlip = true;
            oo.mPreserveBoundary = true;
            oo.mFrozen = frozen;
            // The flips by scaled Jacobian this mesher was tuned with; it keeps
            // its own dihedral floor below.
            oo.mObjective = OptimizeVolumeObjective::Energy;
            const OptimizeVolumeResult r = optimize_volume(m, oo);
            if (r.mNumFlips > 0 && r.mMesh.NumPoints() == np) {
                const auto block = r.mMesh.Cells(0);
                const std::int64_t* c = block.Conn().As<std::int64_t>();
                std::vector<std::array<std::int64_t, 4>> flipped(block.NumCells());
                for (std::size_t n = 0; n < flipped.size(); ++n)
                    for (std::size_t k = 0; k < 4; ++k)
                        flipped[n][k] = c[n * 4 + k];
                std::vector<double> fa, fb;
                const auto [now, now_hi] = octmo_angles(pts, tets, angle, big);
                const auto [flipped_lo, flipped_hi] = octmo_angles(pts, flipped, fa, fb);
                if (flipped_lo >= std::min(now, floor) && flipped_hi <= std::max(now_hi, ceiling)) {
                    tets.swap(flipped);
                    flips += r.mNumFlips;
                }
            }
        }
    }
    // Write back: tetrahedra that kept their vertices and did not move keep
    // their lattice shape; the rest are optimised shapes.
    std::map<std::array<std::int64_t, 4>, std::size_t> old;
    for (std::size_t n = 0; n < kept.size(); ++n) {
        auto k = rTets.mTets[kept[n]];
        std::sort(k.begin(), k.end());
        old.emplace(k, kept[n]);
    }
    OctreeMeshTets out;
    out.mPoints = rTets.mPoints;
    out.mOnSkin = rTets.mOnSkin;
    out.mFeature = rTets.mFeature;
    out.mChainGaps = rTets.mChainGaps;
    out.mNumMaterials = rTets.mNumMaterials;
    const LinearOctree& tree = rS.mOctree.mTree;
    for (const auto& t : tets) {
        auto k = t;
        std::sort(k.begin(), k.end());
        const auto it = old.find(k);
        bool still = it != old.end();
        for (std::int64_t v : t)
            still = still && !moved[static_cast<std::size_t>(v)];
        out.mTets.push_back(t);
        if (it != old.end()) {
            out.mMaterial.push_back(rTets.mMaterial[it->second]);
            out.mLevel.push_back(rTets.mLevel[it->second]);
        } else {
            out.mMaterial.push_back(lo_mat);
            // The level of the leaf holding the centroid.
            Vec3 c{0.0, 0.0, 0.0};
            for (std::int64_t v : t)
                c = vec3_add(c, vec3_scale(pts[static_cast<std::size_t>(v)], 0.25));
            const double fine = rS.mFrame.Fine();
            std::int64_t x[3];
            for (std::size_t q = 0; q < 3; ++q)
                x[q] = std::clamp<std::int64_t>(
                    static_cast<std::int64_t>(std::floor((c[q] - rS.mFrame.mOrigin[q]) / fine)), 0,
                    tree.Extent(q) - 1);
            const std::int64_t j = tree.Locate(x[0], x[1], x[2]);
            out.mLevel.push_back(j < 0 ? 0 : tree.mLevel[static_cast<std::size_t>(j)]);
        }
        out.mShape.push_back(still ? rTets.mShape[it->second] : 0);
        out.mCut.push_back(still ? rTets.mCut[it->second] : std::uint8_t{1});
    }
    std::int64_t nmoved = 0;
    for (std::uint8_t m : moved)
        nmoved += m;
    rR.mNumOptimizedMoves = nmoved;
    rR.mNumFlips = flips;
    rTets = std::move(out);
}

}  // namespace detail
}  // namespace meshioplusplus
