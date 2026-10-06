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
// Polyhedral coarsening: greedy seed-and-grow over the shared-face dual
// (detail::build_global_faces), one merged polyhedron per group -- see
// operations/agglomerate.hpp for the contract.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kAggPrefix = "meshio++: agglomerate: ";

/// A staged output block: either an unchanged copy of a pass-through input
/// block (any of its three storage shapes), or the one merged polyhedron
/// block. Mirrors subdivide.cpp's own staging shape; not shared with it,
/// since that one is file-private there too.
struct AggOutBlock {
    std::string mType;
    std::size_t mNodesPerCell = 0;  // rectangular pass-through only
    std::vector<std::int64_t> mConn;
    std::vector<std::vector<std::int64_t>> mPolygonRows;
    std::vector<std::vector<std::vector<std::int64_t>>> mPolyhedronCells;
    bool mIsRagged = false;
    bool mIsPolyhedron = false;
};

AggOutBlock agg_stage_passthrough(const Mesh::CellView& rBlock) {
    AggOutBlock out;
    out.mType = std::string(rBlock.Type());
    if (rBlock.IsPolyhedron()) {
        out.mIsRagged = true;
        out.mIsPolyhedron = true;
        out.mPolyhedronCells.resize(rBlock.NumCells());
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c) {
            for (std::size_t f = 0; f < rBlock.NumFaces(c); ++f) {
                const auto face = rBlock.Face(c, f);
                out.mPolyhedronCells[c].emplace_back(face.first, face.first + face.second);
            }
        }
    } else if (rBlock.IsRagged()) {
        out.mIsRagged = true;
        out.mPolygonRows.resize(rBlock.NumCells());
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c) {
            const std::int64_t* row = rBlock.Row(c);
            out.mPolygonRows[c].assign(row, row + rBlock.RowSize(c));
        }
    } else {
        const NDArray& conn = rBlock.Conn();
        const std::size_t npc = rBlock.NodesPerCell();
        out.mNodesPerCell = npc;
        out.mConn.resize(rBlock.NumCells() * npc);
        for (std::size_t c = 0; c < rBlock.NumCells(); ++c)
            for (std::size_t k = 0; k < npc; ++k)
                out.mConn[c * npc + k] = detail::read_int(conn, c * npc + k);
    }
    return out;
}

void agg_emit_block(Mesh& rOut, AggOutBlock& rBlock) {
    if (rBlock.mIsPolyhedron) {
        rOut.AddPolyhedronBlock(rBlock.mType, std::move(rBlock.mPolyhedronCells));
    } else if (rBlock.mIsRagged) {
        rOut.AddPolygonBlock(rBlock.mType, std::move(rBlock.mPolygonRows));
    } else {
        const std::size_t npc = rBlock.mNodesPerCell;
        const std::size_t ncells = npc == 0 ? 0 : rBlock.mConn.size() / npc;
        NDArray conn = NDArray::Uninit(DType::Int64, {ncells, npc});
        std::int64_t* dst = conn.As<std::int64_t>();
        for (std::size_t i = 0; i < rBlock.mConn.size(); ++i)
            dst[i] = rBlock.mConn[i];
        rOut.AddCellBlock(rBlock.mType, std::move(conn));
    }
}

/// Area of global face `f`, gathering its corner coordinates into a local
/// buffer first -- GlobalFaces::Face() returns GLOBAL node ids, while
/// polygon_area's ring-taking overload expects local indices, so the plain
/// (no-ring) overload over a freshly gathered buffer is the fit.
double agg_face_area(const detail::GlobalFaces& rFaces, std::size_t f, const NDArray& rPoints,
                     std::size_t PointDim) {
    const std::size_t n = rFaces.FaceSize(f);
    const std::int64_t* ring = rFaces.Face(f);
    std::vector<detail::Vec3> coords(n);
    for (std::size_t k = 0; k < n; ++k) {
        const auto pid = static_cast<std::size_t>(ring[k]);
        for (std::size_t d = 0; d < 3; ++d)
            coords[k][d] = d < PointDim ? detail::read_double(rPoints, pid * PointDim + d) : 0.0;
    }
    return detail::polygon_area(coords.data(), n);
}

/// Newell area vector and vertex centroid of global face `f`, as stored
/// (wound out of its owner).
void agg_face_geometry(const detail::GlobalFaces& rFaces, std::size_t f, const NDArray& rPoints,
                       std::size_t PointDim, detail::Vec3& rArea, detail::Vec3& rCentroid) {
    const std::size_t n = rFaces.FaceSize(f);
    const std::int64_t* ring = rFaces.Face(f);
    std::vector<detail::Vec3> p(n);
    for (std::size_t k = 0; k < n; ++k) {
        const auto pid = static_cast<std::size_t>(ring[k]);
        for (std::size_t d = 0; d < 3; ++d)
            p[k][d] = d < PointDim ? detail::read_double(rPoints, pid * PointDim + d) : 0.0;
    }
    rArea = {0.0, 0.0, 0.0};
    rCentroid = {0.0, 0.0, 0.0};
    for (std::size_t k = 0; k < n; ++k) {
        rArea = detail::vec3_add(rArea, detail::vec3_cross(p[k], p[(k + 1) % n]));
        rCentroid = detail::vec3_add(rCentroid, p[k]);
    }
    rArea = detail::vec3_scale(rArea, 0.5);
    rCentroid = detail::vec3_scale(rCentroid, n ? 1.0 / static_cast<double>(n) : 0.0);
}

/// Sphericity `pi^(1/3) (6V)^(2/3) / A`: 1 for a ball.
double agg_sphericity(double Volume, double Area) {
    if (!(Area > 0.0) || !(Volume > 0.0))
        return 0.0;
    return std::cbrt(3.14159265358979323846) * std::pow(6.0 * Volume, 2.0 / 3.0) / Area;
}

/// Whether a polyhedron's faces close up: every undirected edge used by
/// exactly two of its faces.
bool agg_closed(const std::vector<std::vector<std::int64_t>>& rFaces) {
    std::vector<std::array<std::int64_t, 2>> edges;
    for (const auto& face : rFaces)
        for (std::size_t k = 0; k < face.size(); ++k) {
            std::int64_t u = face[k];
            std::int64_t v = face[(k + 1) % face.size()];
            if (u == v)
                continue;
            if (u > v)
                std::swap(u, v);
            edges.push_back({u, v});
        }
    std::sort(edges.begin(), edges.end());
    for (std::size_t i = 0; i < edges.size();) {
        std::size_t j = i + 1;
        while (j < edges.size() && edges[j] == edges[i])
            ++j;
        if (j - i != 2)
            return false;
        i = j;
    }
    return true;
}

/// Coplanar patches of external faces. Faces separating the same two sides
/// (group, group-or-boundary) that share an edge and lie on one plane fuse
/// into one ring, wound out of `mSideA` (the lower group; the group itself
/// against the boundary).
struct AggPatches {
    std::vector<std::int64_t> mPatchOfFace;  // per global face, -1 when unfused
    std::vector<std::vector<std::int64_t>> mRing;
    std::vector<std::int64_t> mSideA;
    std::vector<std::int64_t> mSideB;  // -1 for the boundary
    std::vector<std::size_t> mNumFaces;
};

AggPatches agg_coplanar_patches(const detail::GlobalFaces& rFaces,
                                const std::vector<std::int64_t>& rGroupOf,
                                const std::vector<detail::Vec3>& rAreaVec, double CosTol,
                                const std::vector<std::uint8_t>& rRejectedFace) {
    const std::size_t nf = rFaces.NumFaces();
    AggPatches out;
    out.mPatchOfFace.assign(nf, -1);
    // Side keys and orientation of every external face.
    std::vector<std::int64_t> side_a(nf, -2);
    std::vector<std::int64_t> side_b(nf, -2);
    std::vector<std::uint8_t> flip(nf, 0);
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t go = rGroupOf[static_cast<std::size_t>(rFaces.mOwner[f])];
        const std::int64_t nb = rFaces.mNeighbour[f];
        const std::int64_t gn = nb >= 0 ? rGroupOf[static_cast<std::size_t>(nb)] : -1;
        if (go == gn || rRejectedFace[f])
            continue;
        side_a[f] = gn < 0 ? go : std::min(go, gn);
        side_b[f] = gn < 0 ? -1 : std::max(go, gn);
        flip[f] = go == side_a[f] ? 0 : 1;
    }
    auto unit = [&](std::size_t f) {
        detail::Vec3 n = detail::vec3_normalize(rAreaVec[f]);
        return flip[f] ? detail::vec3_scale(n, -1.0) : n;
    };
    // Faces grouped by (side a, side b, undirected edge).
    struct EdgeUse {
        std::int64_t mA, mB, mLo, mHi, mFace;
        bool operator<(const EdgeUse& o) const {
            return std::tie(mA, mB, mLo, mHi, mFace) < std::tie(o.mA, o.mB, o.mLo, o.mHi, o.mFace);
        }
    };
    std::vector<EdgeUse> uses;
    for (std::size_t f = 0; f < nf; ++f) {
        if (side_a[f] == -2)
            continue;
        const std::size_t n = rFaces.FaceSize(f);
        const std::int64_t* ring = rFaces.Face(f);
        for (std::size_t k = 0; k < n; ++k) {
            const std::int64_t u = ring[k];
            const std::int64_t v = ring[(k + 1) % n];
            uses.push_back({side_a[f], side_b[f], std::min(u, v), std::max(u, v),
                            static_cast<std::int64_t>(f)});
        }
    }
    std::sort(uses.begin(), uses.end());
    std::vector<std::int64_t> parent(nf);
    for (std::size_t f = 0; f < nf; ++f)
        parent[f] = static_cast<std::int64_t>(f);
    auto find = [&](std::int64_t x) {
        while (parent[static_cast<std::size_t>(x)] != x) {
            parent[static_cast<std::size_t>(x)] =
                parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
            x = parent[static_cast<std::size_t>(x)];
        }
        return x;
    };
    for (std::size_t i = 0; i < uses.size();) {
        std::size_t j = i + 1;
        while (j < uses.size() && uses[j].mA == uses[i].mA && uses[j].mB == uses[i].mB &&
               uses[j].mLo == uses[i].mLo && uses[j].mHi == uses[i].mHi)
            ++j;
        if (j - i == 2) {
            const auto fa = static_cast<std::size_t>(uses[i].mFace);
            const auto fb = static_cast<std::size_t>(uses[i + 1].mFace);
            if (detail::vec3_dot(unit(fa), unit(fb)) >= CosTol) {
                const std::int64_t ra = find(static_cast<std::int64_t>(fa));
                const std::int64_t rb = find(static_cast<std::int64_t>(fb));
                if (ra != rb)
                    parent[static_cast<std::size_t>(std::max(ra, rb))] = std::min(ra, rb);
            }
        }
        i = j;
    }
    // Components with two or more faces, in ascending root order.
    std::vector<std::vector<std::size_t>> members(nf);
    for (std::size_t f = 0; f < nf; ++f)
        if (side_a[f] != -2)
            members[static_cast<std::size_t>(find(static_cast<std::int64_t>(f)))].push_back(f);
    for (std::size_t root = 0; root < nf; ++root) {
        const std::vector<std::size_t>& comp = members[root];
        if (comp.size() < 2)
            continue;
        const detail::Vec3 n0 = unit(comp.front());
        bool planar = true;
        for (std::size_t f : comp)
            planar = planar && detail::vec3_dot(unit(f), n0) >= CosTol;
        if (!planar)
            continue;
        // The outline: directed edges (wound out of side a) used once.
        std::size_t num_directed = 0;
        for (std::size_t f : comp)
            num_directed += rFaces.FaceSize(f);
        std::vector<std::array<std::int64_t, 2>> dir;
        dir.reserve(num_directed);
        for (std::size_t f : comp) {
            const std::size_t n = rFaces.FaceSize(f);
            const std::int64_t* ring = rFaces.Face(f);
            for (std::size_t k = 0; k < n; ++k) {
                std::int64_t u = ring[k];
                std::int64_t v = ring[(k + 1) % n];
                if (flip[f])
                    std::swap(u, v);
                if (u != v)
                    dir.push_back({u, v});
            }
        }
        std::vector<std::array<std::int64_t, 2>> sorted = dir;
        std::sort(sorted.begin(), sorted.end());
        std::vector<std::array<std::int64_t, 2>> outline;
        outline.reserve(dir.size());
        for (const auto& e : dir)
            if (!std::binary_search(sorted.begin(), sorted.end(),
                                    std::array<std::int64_t, 2>{e[1], e[0]}))
                outline.push_back(e);
        std::sort(outline.begin(), outline.end());
        // A single simple loop: every vertex leaves once and is reached once.
        bool simple = outline.size() >= 3;
        for (std::size_t k = 1; simple && k < outline.size(); ++k)
            simple = outline[k][0] != outline[k - 1][0];
        std::vector<std::int64_t> heads;
        heads.reserve(outline.size());
        for (const auto& e : outline)
            heads.push_back(e[1]);
        std::sort(heads.begin(), heads.end());
        for (std::size_t k = 1; simple && k < heads.size(); ++k)
            simple = heads[k] != heads[k - 1];
        std::vector<std::int64_t> ring;
        if (simple) {
            ring.reserve(outline.size());
            std::int64_t at = outline.front()[0];
            for (std::size_t steps = 0; steps < outline.size(); ++steps) {
                ring.push_back(at);
                const auto it = std::lower_bound(outline.begin(), outline.end(),
                                                 std::array<std::int64_t, 2>{at, INT64_MIN});
                if (it == outline.end() || (*it)[0] != at) {
                    simple = false;
                    break;
                }
                at = (*it)[1];
            }
            simple = simple && at == outline.front()[0];
        }
        if (!simple)
            continue;
        const auto pid = static_cast<std::int64_t>(out.mRing.size());
        for (std::size_t f : comp)
            out.mPatchOfFace[f] = pid;
        out.mRing.push_back(std::move(ring));
        out.mSideA.push_back(side_a[comp.front()]);
        out.mSideB.push_back(side_b[comp.front()]);
        out.mNumFaces.push_back(comp.size());
    }
    return out;
}

/// Frontier entry ordered by DESCENDING accumulated shared-face area, ties
/// broken by ASCENDING compact cell id -- storing the negated area keeps a
/// plain ascending std::set a max-by-area, min-by-id priority structure.
struct FrontierKey {
    double mNegArea;
    std::int64_t mId;
    bool operator<(const FrontierKey& rOther) const {
        if (mNegArea != rOther.mNegArea)
            return mNegArea < rOther.mNegArea;
        return mId < rOther.mId;
    }
};

}  // namespace

AgglomerateResult agglomerate(const Mesh& rMesh, const AgglomerateOptions& rOptions) {
    if (rOptions.mTargetGroupSize == 0)
        throw std::invalid_argument(std::string(kAggPrefix) + "mTargetGroupSize must be >= 1");
    if (!(rOptions.mCoplanarAngleDeg >= 0.0 && rOptions.mCoplanarAngleDeg < 90.0))
        throw std::invalid_argument(std::string(kAggPrefix) +
                                    "mCoplanarAngleDeg must lie in [0, 90)");
    if (!(rOptions.mMinSphericity >= 0.0 && rOptions.mMinSphericity <= 1.0))
        throw std::invalid_argument(std::string(kAggPrefix) + "mMinSphericity must lie in [0, 1]");

    const detail::GlobalFaces gf = detail::build_global_faces(rMesh);
    if (gf.mNumNonManifold > 0)
        throw std::invalid_argument(
            std::string(kAggPrefix) + "mesh contains " + std::to_string(gf.mNumNonManifold) +
            " face(s) shared by three or more cells (non-manifold); refusing rather than "
            "guessing a boundary classification");

    const std::size_t n_compact = gf.NumCells();
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();

    // Every face's area, once and in parallel: the grow loop below used to
    // recompute one each time a face was pushed.
    std::vector<double> face_area(gf.NumFaces());
    parallel_for(gf.NumFaces(),
                 [&](std::size_t f) { face_area[f] = agg_face_area(gf, f, points, pdim); });

    // Face area vectors and centroids, for the cells' volumes (the gate) and
    // the face normals (coplanar merging); only when either is asked for.
    const bool gate = rOptions.mMinSphericity > 0.0;
    std::vector<detail::Vec3> area_vec;
    std::vector<double> cell_volume;
    std::vector<double> cell_area;
    if (gate || rOptions.mMergeCoplanarFaces) {
        area_vec.resize(gf.NumFaces());
        std::vector<detail::Vec3> centroid(gf.NumFaces());
        parallel_for(gf.NumFaces(), [&](std::size_t f) {
            agg_face_geometry(gf, f, points, pdim, area_vec[f], centroid[f]);
        });
        if (gate) {
            // Divergence theorem over each cell's outward faces.
            cell_volume.assign(n_compact, 0.0);
            cell_area.assign(n_compact, 0.0);
            parallel_for(n_compact, [&](std::size_t c) {
                const std::size_t nfc = gf.NumCellFaces(c);
                const std::int64_t* row = gf.CellFaces(c);
                double v = 0.0;
                double a = 0.0;
                for (std::size_t k = 0; k < nfc; ++k) {
                    const std::int64_t sid = row[k];
                    const auto f = static_cast<std::size_t>((sid > 0 ? sid : -sid) - 1);
                    const double d = detail::vec3_dot(centroid[f], area_vec[f]) / 3.0;
                    v += sid > 0 ? d : -d;
                    a += face_area[f];
                }
                cell_volume[c] = v;
                cell_area[c] = a;
            });
        }
    }
    std::int64_t num_rejected = 0;

    // --- greedy seed-and-grow over the face dual --------------------------
    std::vector<std::int64_t> group_of(n_compact, -1);
    std::vector<std::vector<std::int64_t>> groups;

    // A candidate's accumulated shared area, as dense per-cell arrays reset
    // through the ids a seed touched -- no hash map allocated per seed. Same
    // sums, in the same order.
    std::vector<double> pending(n_compact, 0.0);
    // 0 not on the frontier, 1 on it, kRefused refused by the gate for this
    // seed's group (never pushed again: a re-push would restart its shared
    // area from one face and misjudge the union).
    constexpr std::uint8_t kRefused = 2;
    std::vector<std::uint8_t> is_pending(n_compact, 0);
    std::vector<std::int64_t> touched;
    std::set<FrontierKey> frontier;

    for (std::size_t seed = 0; seed < n_compact; ++seed) {
        if (group_of[seed] != -1)
            continue;
        const auto gid = static_cast<std::int64_t>(groups.size());
        std::vector<std::int64_t> members{static_cast<std::int64_t>(seed)};
        group_of[seed] = gid;
        for (std::int64_t t : touched)
            is_pending[static_cast<std::size_t>(t)] = 0;
        touched.clear();
        frontier.clear();

        auto push_neighbours = [&](std::int64_t c) {
            const std::size_t nf = gf.NumCellFaces(static_cast<std::size_t>(c));
            const std::int64_t* row = gf.CellFaces(static_cast<std::size_t>(c));
            for (std::size_t k = 0; k < nf; ++k) {
                const std::int64_t sid = row[k];
                const auto f = static_cast<std::size_t>((sid > 0 ? sid : -sid) - 1);
                const std::int64_t owner = gf.mOwner[f];
                const std::int64_t neigh = gf.mNeighbour[f];
                const std::int64_t other = (owner == c) ? neigh : owner;
                if (other < 0)
                    continue;  // mesh boundary
                if (group_of[static_cast<std::size_t>(other)] != -1)
                    continue;  // already claimed (by this group or would be a bug otherwise)
                if (is_pending[static_cast<std::size_t>(other)] == kRefused)
                    continue;  // refused for this group
                const double a = face_area[f];
                double& acc = pending[static_cast<std::size_t>(other)];
                if (is_pending[static_cast<std::size_t>(other)]) {
                    frontier.erase(FrontierKey{-acc, other});
                    acc += a;
                } else {
                    acc = a;
                    is_pending[static_cast<std::size_t>(other)] = 1;
                    touched.push_back(other);
                }
                frontier.insert(FrontierKey{-acc, other});
            }
        };

        push_neighbours(static_cast<std::int64_t>(seed));
        double group_volume = gate ? cell_volume[seed] : 0.0;
        double group_area = gate ? cell_area[seed] : 0.0;

        while (members.size() < rOptions.mTargetGroupSize && !frontier.empty()) {
            const auto fit = frontier.begin();
            const std::int64_t c = fit->mId;
            frontier.erase(fit);
            is_pending[static_cast<std::size_t>(c)] = 0;
            if (group_of[static_cast<std::size_t>(c)] != -1)
                continue;  // defensive; unreachable given the push-time check above
            if (gate) {
                // The union's volume and external area: the shared faces leave
                // the surface from both sides.
                const auto ci = static_cast<std::size_t>(c);
                const double v = group_volume + cell_volume[ci];
                const double a = group_area + cell_area[ci] - 2.0 * pending[ci];
                if (agg_sphericity(v, a) < rOptions.mMinSphericity) {
                    ++num_rejected;
                    is_pending[ci] = kRefused;  // skipped for this group only
                    continue;
                }
                group_volume = v;
                group_area = a;
            }
            group_of[static_cast<std::size_t>(c)] = gid;
            members.push_back(c);
            push_neighbours(c);
        }

        std::sort(members.begin(), members.end());
        groups.push_back(std::move(members));
    }

    // --- coplanar patches (optional) -----------------------------------------
    AggPatches patches;
    patches.mPatchOfFace.assign(gf.NumFaces(), -1);
    std::vector<std::uint8_t> rejected_face(gf.NumFaces(), 0);
    const double cos_tol = std::cos(rOptions.mCoplanarAngleDeg * 3.14159265358979323846 / 180.0);
    std::int64_t num_faces_merged = 0;

    // --- emit: one polyhedron cell per group, external faces only ---------
    std::vector<std::vector<std::vector<std::int64_t>>> merged_cells(groups.size());
    for (int attempt = 0; attempt < 3; ++attempt) {
        if (rOptions.mMergeCoplanarFaces)
            patches = agg_coplanar_patches(gf, group_of, area_vec, cos_tol, rejected_face);
        num_faces_merged = 0;
        for (std::size_t g = 0; g < groups.size(); ++g) {
            auto& faces_out = merged_cells[g];
            faces_out.clear();
            std::vector<std::uint8_t> patch_done(patches.mRing.size(), 0);
            for (std::int64_t c : groups[g]) {
                const std::size_t nf = gf.NumCellFaces(static_cast<std::size_t>(c));
                const std::int64_t* row = gf.CellFaces(static_cast<std::size_t>(c));
                for (std::size_t k = 0; k < nf; ++k) {
                    const std::int64_t sid = row[k];
                    const auto f = static_cast<std::size_t>((sid > 0 ? sid : -sid) - 1);
                    const std::int64_t owner = gf.mOwner[f];
                    const std::int64_t neigh = gf.mNeighbour[f];
                    const std::int64_t other = (owner == c) ? neigh : owner;
                    if (other >= 0 && group_of[static_cast<std::size_t>(other)] ==
                                          group_of[static_cast<std::size_t>(c)])
                        continue;  // internal to the group: dropped from both sides
                    const std::int64_t pid = patches.mPatchOfFace[f];
                    if (pid >= 0) {
                        // A fused patch: its ring once, wound out of this group.
                        const auto pi = static_cast<std::size_t>(pid);
                        if (!patch_done[pi]) {
                            patch_done[pi] = 1;
                            std::vector<std::int64_t> ring_nodes = patches.mRing[pi];
                            if (patches.mSideA[pi] != static_cast<std::int64_t>(g))
                                std::reverse(ring_nodes.begin(), ring_nodes.end());
                            faces_out.push_back(std::move(ring_nodes));
                            num_faces_merged +=
                                static_cast<std::int64_t>(patches.mNumFaces[pi]) - 1;
                        }
                        continue;
                    }
                    const std::size_t n = gf.FaceSize(f);
                    const std::int64_t* ring = gf.Face(f);
                    std::vector<std::int64_t> face_nodes(ring, ring + n);
                    if (sid < 0)
                        std::reverse(face_nodes.begin(), face_nodes.end());
                    faces_out.push_back(std::move(face_nodes));
                }
            }
        }
        if (patches.mRing.empty())
            break;
        // A fusion must leave every polyhedron it touches closed; a patch of a
        // group that no longer closes is withdrawn (both sides) and the emit rerun.
        bool all_closed = true;
        for (std::size_t g = 0; g < groups.size(); ++g) {
            if (agg_closed(merged_cells[g]))
                continue;
            all_closed = false;
            for (std::size_t f = 0; f < gf.NumFaces(); ++f) {
                const std::int64_t pid = patches.mPatchOfFace[f];
                if (pid < 0)
                    continue;
                const auto pi = static_cast<std::size_t>(pid);
                if (patches.mSideA[pi] == static_cast<std::int64_t>(g) ||
                    patches.mSideB[pi] == static_cast<std::int64_t>(g))
                    rejected_face[f] = 1;
            }
        }
        if (all_closed)
            break;
        if (attempt == 2) {
            // Give up on fusion rather than emit an open polyhedron.
            log::warn("{}coplanar merging left an open polyhedron; faces kept unmerged",
                      kAggPrefix);
            std::fill(rejected_face.begin(), rejected_face.end(), std::uint8_t{1});
            attempt = 1;  // one more pass, with every face rejected
        }
    }

    // --- which original blocks carry volume, and the compact<->global bridge
    std::vector<bool> is_volume(rMesh.NumCellBlocks(), true);
    for (std::size_t b : gf.mNonCellBlocks)
        if (b < is_volume.size())
            is_volume[b] = false;

    const std::vector<std::int64_t> in_bases = detail::block_bases(rMesh);
    std::unordered_map<std::int64_t, std::int64_t> global_to_compact;
    global_to_compact.reserve(n_compact);
    for (std::size_t c = 0; c < n_compact; ++c)
        global_to_compact[gf.mCellToGlobal[c]] = static_cast<std::int64_t>(c);

    // --- build the output mesh's blocks, tracking where each original block
    // landed (pass-through) or whether it fed the merged block -----------
    Mesh out;
    out.AssignPoints(detail::data_owned_copy(points));

    const std::size_t nblocks_in = rMesh.NumCellBlocks();
    std::vector<std::int64_t> pass_out_block(nblocks_in, -1);
    std::int64_t merged_out_block = -1;
    bool merged_emitted = false;
    std::size_t bi = 0;
    std::int64_t out_block_idx = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (is_volume[bi]) {
            if (!merged_emitted) {
                AggOutBlock mb;
                mb.mType = "polyhedron";
                mb.mIsRagged = true;
                mb.mIsPolyhedron = true;
                mb.mPolyhedronCells = std::move(merged_cells);
                if (!mb.mPolyhedronCells.empty()) {
                    agg_emit_block(out, mb);
                    merged_out_block = out_block_idx++;
                }
                merged_emitted = true;
            }
        } else {
            AggOutBlock pb = agg_stage_passthrough(cb);
            agg_emit_block(out, pb);
            pass_out_block[bi] = out_block_idx++;
        }
        ++bi;
    }

    const std::vector<std::int64_t> out_bases = detail::block_bases(out);

    // --- the flat cell map --------------------------------------------------
    NDArray cell_map =
        NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(detail::total_cells(in_bases))});
    std::int64_t* cm = cell_map.As<std::int64_t>();
    bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t base_in = in_bases[bi];
        if (is_volume[bi]) {
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                const std::int64_t global_in = base_in + static_cast<std::int64_t>(r);
                const std::int64_t compact = global_to_compact.at(global_in);
                cm[global_in] = merged_out_block < 0
                                    ? -1
                                    : out_bases[static_cast<std::size_t>(merged_out_block)] +
                                          group_of[static_cast<std::size_t>(compact)];
            }
        } else {
            const std::int64_t ob = pass_out_block[bi];
            for (std::size_t r = 0; r < cb.NumCells(); ++r) {
                const std::int64_t global_in = base_in + static_cast<std::int64_t>(r);
                cm[global_in] =
                    out_bases[static_cast<std::size_t>(ob)] + static_cast<std::int64_t>(r);
            }
        }
        ++bi;
    }

    // --- point_data / field_data: unchanged, no new or pruned points -------
    for (const std::string& name : rMesh.PointDataNames())
        out.AddPointData(name, detail::data_owned_copy(rMesh.PointData(name)));
    for (const std::string& name : rMesh.FieldDataNames())
        out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

    // --- cell_data: pass-through blocks copied verbatim; the merged block's
    // row for group g is its first (ascending compact id) member's row -----
    for (const std::string& name : rMesh.CellDataNames()) {
        const std::size_t ndata = rMesh.CellDataNumBlocks(name);
        if (ndata != nblocks_in) {
            log::warn(
                "{}cell_data '{}' does not have one array per input block; dropped rather "
                "than guessed at",
                kAggPrefix, name);
            continue;
        }
        std::vector<NDArray> out_blocks(static_cast<std::size_t>(out_block_idx));
        bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            (void)cb;
            if (!is_volume[bi] && pass_out_block[bi] >= 0)
                out_blocks[static_cast<std::size_t>(pass_out_block[bi])] =
                    detail::data_owned_copy(rMesh.CellData(name, bi));
            ++bi;
        }
        if (merged_out_block >= 0) {
            // Resolve each group's first member's (block, row) and gather.
            const NDArray& first_src = rMesh.CellData(name, 0);
            std::vector<std::size_t> shape = first_src.Shape();
            shape[0] = groups.size();
            NDArray dst = NDArray::Uninit(first_src.Dtype(), std::move(shape));
            const std::size_t row_bytes =
                first_src.Nbytes() / std::max<std::size_t>(1, detail::rows(first_src));
            std::byte* p_dst = dst.Data();
            for (std::size_t g = 0; g < groups.size(); ++g) {
                const std::int64_t compact_first = groups[g].front();
                const std::int64_t global =
                    gf.mCellToGlobal[static_cast<std::size_t>(compact_first)];
                const auto [blk, row] = detail::global_to_block_row(in_bases, global);
                const NDArray& src = rMesh.CellData(name, blk);
                const std::byte* p_src = src.Data();
                std::memcpy(p_dst + g * row_bytes,
                            p_src + static_cast<std::size_t>(row) * row_bytes, row_bytes);
            }
            out_blocks[static_cast<std::size_t>(merged_out_block)] = std::move(dst);
        }
        out.AddCellData(name, std::move(out_blocks));
    }

    AgglomerateResult res;
    res.mMesh = std::move(out);
    res.mCellMap = std::move(cell_map);
    res.mNumFacesMerged = num_faces_merged;
    res.mNumRejected = num_rejected;

    detail::RegionRemap rmap;
    rmap.mCellMapKind = detail::CellMapKind::Global;
    rmap.pGlobalCellMap = &res.mCellMap;
    rmap.mOpName = "agglomerate";
    detail::remap_regions(rMesh, res.mMesh, rmap);

    return res;
}

}  // namespace meshioplusplus
