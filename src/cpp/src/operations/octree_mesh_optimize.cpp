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
#include <utility>
#include <vector>

// Project includes
#include "../detail/bcc_octree_tets.hpp"
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

double octmo_min_dihedral(const std::vector<Vec3>& rP, const std::array<std::int64_t, 4>& rT) {
    const Vec3& a = rP[static_cast<std::size_t>(rT[0])];
    const Vec3& b = rP[static_cast<std::size_t>(rT[1])];
    const Vec3& c = rP[static_cast<std::size_t>(rT[2])];
    const Vec3& d = rP[static_cast<std::size_t>(rT[3])];
    if (triple_product(vec3_sub(b, a), vec3_sub(c, a), vec3_sub(d, a)) <= 0.0)
        return -1.0;
    return tet_dihedral_angles(a, b, c, d)[0];
}

double octmo_max_dihedral(const std::vector<Vec3>& rP, const std::array<std::int64_t, 4>& rT) {
    const Vec3& a = rP[static_cast<std::size_t>(rT[0])];
    const Vec3& b = rP[static_cast<std::size_t>(rT[1])];
    const Vec3& c = rP[static_cast<std::size_t>(rT[2])];
    const Vec3& d = rP[static_cast<std::size_t>(rT[3])];
    if (triple_product(vec3_sub(b, a), vec3_sub(c, a), vec3_sub(d, a)) <= 0.0)
        return 181.0;
    return tet_dihedral_angles(a, b, c, d)[5];
}

// The worst dihedral angle of every tetrahedron, and its minimum.
double octmo_angles(const std::vector<Vec3>& rP, const std::vector<std::array<std::int64_t, 4>>& rT,
                    std::vector<double>& rOut) {
    rOut.resize(rT.size());
    parallel_for(rT.size(), [&](std::size_t t) { rOut[t] = octmo_min_dihedral(rP, rT[t]); });
    double lo = 180.0;
    for (double a : rOut)
        lo = std::min(lo, a);
    return lo;
}

// The largest dihedral angle of every tetrahedron, and its maximum.
double octmo_max_angles(const std::vector<Vec3>& rP,
                        const std::vector<std::array<std::int64_t, 4>>& rT,
                        std::vector<double>& rOut) {
    rOut.resize(rT.size());
    parallel_for(rT.size(), [&](std::size_t t) { rOut[t] = octmo_max_dihedral(rP, rT[t]); });
    double hi = 0.0;
    for (double a : rOut)
        hi = std::max(hi, a);
    return hi;
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
        const double floor = octmo_angles(pts, tets, angle);
        const double ceiling = octmo_max_angles(pts, tets, big);
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
        octmo_angles(pts, tets, angle);
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
                    const auto& tt = tets[static_cast<std::size_t>(t)];
                    w = std::min(
                        {w, octmo_min_dihedral(pts, tt), 180.0 - octmo_max_dihedral(pts, tt)});
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
        // below the mesh's worst before the round, until none does.
        if (rO.mDihedralFloorGuard) {
            for (int sweep = 0; sweep < 64; ++sweep) {
                octmo_angles(pts, tets, angle);
                octmo_max_angles(pts, tets, big);
                bool any = false;
                for (std::size_t p = 0; p < np; ++p) {
                    if (pts[p] == before[p])
                        continue;
                    for (std::int64_t t : inc[p])
                        if (angle[static_cast<std::size_t>(t)] < floor ||
                            big[static_cast<std::size_t>(t)] > ceiling) {
                            pts[p] = before[p];
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
            Mesh m;
            NDArray xyzs = NDArray::Uninit(DType::Float64, {np, std::size_t{3}});
            for (std::size_t p = 0; p < np; ++p)
                for (std::size_t k = 0; k < 3; ++k)
                    xyzs.As<double>()[p * 3 + k] = pts[p][k];
            NDArray conn = NDArray::Uninit(DType::Int64, {tets.size(), std::size_t{4}});
            for (std::size_t n = 0; n < tets.size(); ++n)
                for (std::size_t c = 0; c < 4; ++c)
                    conn.As<std::int64_t>()[n * 4 + c] = tets[n][c];
            m.AssignPoints(std::move(xyzs));
            m.AddCellBlock("tetra", std::move(conn));
            OptimizeVolumeOptions oo;
            oo.mMaxIterations = 1;
            oo.mRelocate = false;
            oo.mFlip = true;
            oo.mPreserveBoundary = true;
            oo.mFrozen = frozen;
            const OptimizeVolumeResult r = optimize_volume(m, oo);
            if (r.mNumFlips > 0 && r.mMesh.NumPoints() == np) {
                const auto block = r.mMesh.Cells(0);
                const std::int64_t* c = block.Conn().As<std::int64_t>();
                std::vector<std::array<std::int64_t, 4>> flipped(block.NumCells());
                for (std::size_t n = 0; n < flipped.size(); ++n)
                    for (std::size_t k = 0; k < 4; ++k)
                        flipped[n][k] = c[n * 4 + k];
                std::vector<double> fa, fb;
                const double now = octmo_angles(pts, tets, angle);
                const double now_hi = octmo_max_angles(pts, tets, big);
                if (octmo_angles(pts, flipped, fa) >= std::min(now, floor) &&
                    octmo_max_angles(pts, flipped, fb) <= std::max(now_hi, ceiling)) {
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
