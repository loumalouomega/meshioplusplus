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
// The signed-distance kernel. See detail/surface_distance.hpp for the contract,
// in particular why the accelerator's visit order is unobservable and why the
// normal tables are accumulated serially.

// System includes
#include <algorithm>
#ifdef MESHIOPLUSPLUS_SD_COUNT_PROBES
#include <atomic>
#endif
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/point_triangle.hpp"
#include "meshioplusplus/detail/space_filling.hpp"
#include "meshioplusplus/detail/surface_normals.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#ifdef MESHIOPLUSPLUS_SD_COUNT_PROBES
#include "meshioplusplus/log.hpp"
#endif
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "slot_runs.hpp"
#include "surface_edge_runs.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr const char* kSdPrefix = "meshio++: surface distance: ";

// The cells a named Cell region selects, as a per-global-cell flag. An empty
// name means "everything", which is reported as an empty vector rather than an
// all-true one so the caller can skip the test entirely.
std::vector<char> sd_region_mask(const Mesh& rMesh, const std::string& rRegion) {
    if (rRegion.empty())
        return {};
    const std::size_t idx = rMesh.FindRegion(rRegion, RegionKind::Cell);
    if (idx == Mesh::npos) {
        std::string names;
        for (const std::string& n : rMesh.RegionNames())
            names += (names.empty() ? "" : ", ") + n;
        throw std::invalid_argument(std::string(kSdPrefix) + "no cell region named '" + rRegion +
                                    "' (available: " + (names.empty() ? "none" : names) + ")");
    }
    const std::vector<std::int64_t> bases = block_bases(rMesh);
    std::vector<char> mask(static_cast<std::size_t>(total_cells(bases)), 0);
    const NDArray& entries = rMesh.Region(idx).mEntries;
    for (std::size_t e = 0; e < entries.Size(); ++e) {
        const std::int64_t g = read_int(entries, e);
        if (g >= 0 && static_cast<std::size_t>(g) < mask.size())
            mask[static_cast<std::size_t>(g)] = 1;
    }
    return mask;
}

}  // namespace

TriangleSoup build_triangle_soup(const Mesh& rSurface, const std::string& rRegion) {
    TriangleSoup soup;
    const std::size_t dim = rSurface.PointDim();
    const NDArray& points = rSurface.Points();
    const std::size_t npts = rSurface.NumPoints();
    soup.mPoints.resize(npts);
    parallel_for_bw(npts, [&](std::size_t p) {
        soup.mPoints[p] = read_point(points, dim, static_cast<std::int64_t>(p));
    });

    const std::vector<char> mask = sd_region_mask(rSurface, rRegion);
    const std::vector<std::int64_t> bases = block_bases(rSurface);

    // Pass 1 (serial over blocks, parallel over cells): validate each block --
    // in block order, so the first bad block is the one reported -- and count
    // each cell's fan triangles.
    struct Blk {
        Mesh::CellView mCells;
        std::int64_t mBase;
        std::vector<std::int64_t> mFirstTri;  // ncells + 1 prefix
    };
    std::vector<Blk> blocks;
    std::size_t bi = 0;
    for (const auto cb : rSurface.CellRange()) {
        const std::int64_t base = bases[bi++];
        const std::string type(cb.Type());
        const CellType ct = cell_type_from_name(type);

        if (cb.IsPolyhedron() || cell_type_dimension(ct) == 3)
            throw std::invalid_argument(
                std::string(kSdPrefix) + "cell block '" + type +
                "' is a volume; distance is measured to a surface (run extract_surface first)");

        // A block whose cells happen to share a node count stores rectangularly
        // and so is not IsRagged(), which is why the type name is what decides
        // whether it is a polygon -- the trap cgns.cpp records.
        const bool polygon = type.rfind("polygon", 0) == 0;
        if (!polygon && ct != CellType::Triangle && ct != CellType::Quad) {
            if (cell_type_dimension(ct) < 2)
                continue;  // lines and vertices carry no area; silently skipped
            throw std::invalid_argument(std::string(kSdPrefix) + "cell block '" + type +
                                        "' is not a linear surface cell (linearize the mesh "
                                        "first, then run extract_surface if needed)");
        }
        const std::size_t ncells = cb.NumCells();
        std::vector<std::int64_t> ntri(ncells + 1, 0);
        parallel_for(ncells, [&](std::size_t c) {
            const std::int64_t global = base + static_cast<std::int64_t>(c);
            if (!mask.empty() && !mask[static_cast<std::size_t>(global)])
                return;
            const std::size_t n = cb.IsRagged() ? cb.RowSize(c) : cb.NodesPerCell();
            ntri[c + 1] = n < 3 ? 0 : static_cast<std::int64_t>(n - 2);
        });
        for (std::size_t c = 0; c < ncells; ++c)
            ntri[c + 1] += ntri[c];
        blocks.push_back(Blk{cb, base, std::move(ntri)});
    }

    // Pass 2: every cell writes its fan at its own offset -- the triangles
    // come out in (block, cell, fan) order, as the serial append gave them.
    std::size_t total = 0;
    std::vector<std::size_t> block_first(blocks.size());
    for (std::size_t k = 0; k < blocks.size(); ++k) {
        block_first[k] = total;
        total += static_cast<std::size_t>(blocks[k].mFirstTri.back());
    }
    soup.mVertices.resize(total);
    soup.mSourceCell.resize(total);
    soup.mCorners.resize(total * 3);
    for (std::size_t k = 0; k < blocks.size(); ++k) {
        const Blk& b = blocks[k];
        const auto& cb = b.mCells;
        parallel_for(cb.NumCells(), [&](std::size_t c) {
            const std::size_t first = block_first[k] + static_cast<std::size_t>(b.mFirstTri[c]);
            const std::size_t count = static_cast<std::size_t>(b.mFirstTri[c + 1] - b.mFirstTri[c]);
            if (count == 0)
                return;
            // The same fan convert_cells(Simplexify) uses: corner 0 to every
            // non-adjacent edge. Transcribing a different fan here would make
            // the two disagree about which diagonal a quad is split on.
            const auto id = [&](std::size_t i) -> std::int64_t {
                const std::int64_t v =
                    cb.IsRagged() ? cb.Row(c)[i] : read_int(cb.Conn(), c * cb.NodesPerCell() + i);
                if (v < 0 || static_cast<std::size_t>(v) >= npts)
                    throw std::invalid_argument(std::string(kSdPrefix) +
                                                "a cell references a point the mesh does not have");
                return v;
            };
            const std::int64_t a = id(0);
            for (std::size_t j = 0; j < count; ++j) {
                const std::array<std::int64_t, 3> tri{a, id(j + 1), id(j + 2)};
                const std::size_t t = first + j;
                soup.mVertices[t] = tri;
                soup.mSourceCell[t] = b.mBase + static_cast<std::int64_t>(c);
                for (std::size_t i = 0; i < 3; ++i)
                    soup.mCorners[t * 3 + i] = soup.mPoints[static_cast<std::size_t>(tri[i])];
            }
        });
    }
    return soup;
}

SurfaceEdgeRuns surface_edge_runs(const TriangleSoup& rSoup) {
    // One record per (triangle, corner), grouped by a parallel sort
    // (slot_runs.hpp): a run is one edge, its slots in ascending (triangle,
    // corner) order, so its head names the first triangle, as the serial map
    // insert did.
    const std::size_t ntri = rSoup.NumTriangles();
    SurfaceEdgeRuns out;
    out.mKeys.resize(ntri * 3);
    parallel_for(ntri, [&](std::size_t t) {
        const std::array<std::int64_t, 3>& v = rSoup.mVertices[t];
        for (std::size_t e = 0; e < 3; ++e) {
            const std::int64_t u = v[e];
            const std::int64_t w = v[(e + 1) % 3];
            out.mKeys[t * 3 + e] = SurfaceEdgeKey{u < w ? u : w, u < w ? w : u};
        }
    });
    // Bucketed by the lower endpoint, the key's leading component, so the runs
    // come out in ascending key order: the edge map is sorted, as documented.
    const std::size_t npts = rSoup.mPoints.size();
    out.mRuns = group_slots(out.mKeys, npts + 1, [npts](const SurfaceEdgeKey& rK) {
        return rK[0] >= 0 && static_cast<std::size_t>(rK[0]) < npts
                   ? static_cast<std::size_t>(rK[0])
                   : npts;
    });
    return out;
}

SurfaceEdgeMap build_surface_edges(const TriangleSoup& rSoup) {
    return build_surface_edges_from_runs(rSoup, surface_edge_runs(rSoup));
}

SurfaceEdgeMap build_surface_edges_from_runs(const TriangleSoup& rSoup,
                                             const SurfaceEdgeRuns& rRuns) {
    // Per undirected edge: how many triangles use it, and how many use it in the
    // low->high direction. A consistently wound closed surface has every edge
    // used exactly twice, once in each direction.
    const SlotRuns& runs = rRuns.mRuns;
    SurfaceEdgeMap edges(runs.NumRuns());
    parallel_for(runs.NumRuns(), [&](std::size_t r) {
        SurfaceEdgeRecord rec;
        rec.mUsed = static_cast<std::int64_t>(runs.Size(r));
        for (const std::uint64_t* p = runs.Begin(r); p != runs.End(r); ++p) {
            const std::array<std::int64_t, 3>& v = rSoup.mVertices[*p / 3];
            const std::size_t e = static_cast<std::size_t>(*p % 3);
            rec.mForward += v[e] < v[(e + 1) % 3] ? 1 : 0;
        }
        rec.mFirstTriangle = static_cast<std::int64_t>(runs.Head(r) / 3);
        edges[r] = {rRuns.mKeys[runs.Head(r)], rec};
    });
    return edges;
}

SurfaceQuality soup_quality(const TriangleSoup& rSoup) {
    return soup_quality(rSoup, build_surface_edges(rSoup));
}

SurfaceQuality soup_quality(const TriangleSoup& rSoup, const SurfaceEdgeMap& rEdges) {
    SurfaceQuality q;
    const std::size_t ntri = rSoup.NumTriangles();

    for (std::size_t t = 0; t < ntri; ++t) {
        const Vec3& a = rSoup.mCorners[t * 3 + 0];
        const Vec3& b = rSoup.mCorners[t * 3 + 1];
        const Vec3& c = rSoup.mCorners[t * 3 + 2];
        if (!(vec3_norm_sq(vec3_cross(vec3_sub(b, a), vec3_sub(c, a))) > 0.0))
            ++q.mDegenerateTriangles;
    }

    for (const auto& kv : rEdges) {
        const std::int64_t used = kv.second.mUsed;
        const std::int64_t forward = kv.second.mForward;
        if (used == 1)
            ++q.mBoundaryEdges;
        else if (used > 2)
            ++q.mNonManifoldEdges;
        else if (used == 2 && forward != 1)
            ++q.mInconsistentPairs;  // both the same way round: they disagree on "out"
    }
    q.mWatertight = q.mBoundaryEdges == 0 && q.mNonManifoldEdges == 0 &&
                    q.mInconsistentPairs == 0 && q.mDegenerateTriangles == 0;
    return q;
}

namespace {

// The median split of `rBvh.mOrder[Begin, End)` in preorder, leaves of at
// most @p LeafSize triangles; returns the node's index. The depth is
// log2(n / LeafSize), so the recursion is shallow.
std::int64_t sd_bvh_split(SurfaceBvh& rBvh, std::int64_t Begin, std::int64_t End,
                          std::int64_t LeafSize) {
    const auto at = static_cast<std::int64_t>(rBvh.mNodes.size());
    rBvh.mNodes.emplace_back();
    rBvh.mNodes.back().mBegin = Begin;
    rBvh.mNodes.back().mCount = End - Begin;
    if (End - Begin <= LeafSize)
        return at;
    const std::int64_t mid = Begin + (End - Begin) / 2;
    sd_bvh_split(rBvh, Begin, mid, LeafSize);
    const std::int64_t right = sd_bvh_split(rBvh, mid, End, LeafSize);
    rBvh.mNodes[static_cast<std::size_t>(at)].mRight = right;
    return at;
}

// The accelerator: the triangles ordered by the Morton key of their box
// centres (a stable radix argsort, so ties stay in triangle order), split at
// the median count, then every box from its triangles' (leaves) or its
// children's (inner nodes, which come before them in preorder). A function of
// the triangles alone, so the tree is the same on every thread count.
SurfaceBvh sd_build_bvh(const std::vector<Vec3>& rLo, const std::vector<Vec3>& rHi) {
    constexpr std::int64_t kLeafSize = 4;
    SurfaceBvh bvh;
    const std::size_t n = rLo.size();
    if (n == 0)
        return bvh;
    Vec3 lo = rLo[0];
    Vec3 hi = rHi[0];
    for (std::size_t t = 0; t < n; ++t)
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], rLo[t][k]);
            hi[k] = std::max(hi[k], rHi[t][k]);
        }
    const double scale = static_cast<double>((std::uint32_t{1} << sfc_bits) - 1);
    std::vector<std::uint64_t> key(n);
    parallel_for(n, [&](std::size_t t) {
        std::uint32_t q[3];
        for (std::size_t k = 0; k < 3; ++k) {
            const double ext = hi[k] - lo[k];
            const double c = 0.5 * (rLo[t][k] + rHi[t][k]);
            double f = ext > 0.0 ? (c - lo[k]) / ext : 0.0;
            if (!(f >= 0.0))
                f = 0.0;  // also a NaN, which a non-finite corner gives
            q[k] = static_cast<std::uint32_t>(std::min(f, 1.0) * scale);
        }
        key[t] = sfc_morton_key(q);
    });
    bvh.mOrder = sfc_stable_argsort(key);
    sd_bvh_split(bvh, 0, static_cast<std::int64_t>(n), kLeafSize);
    for (std::size_t i = bvh.mNodes.size(); i-- > 0;) {
        SurfaceBvhNode& rNode = bvh.mNodes[i];
        if (rNode.mRight < 0) {
            const auto first =
                static_cast<std::size_t>(bvh.mOrder[static_cast<std::size_t>(rNode.mBegin)]);
            rNode.mLo = rLo[first];
            rNode.mHi = rHi[first];
            for (std::int64_t k = 1; k < rNode.mCount; ++k) {
                const auto t = static_cast<std::size_t>(
                    bvh.mOrder[static_cast<std::size_t>(rNode.mBegin + k)]);
                for (std::size_t d = 0; d < 3; ++d) {
                    rNode.mLo[d] = std::min(rNode.mLo[d], rLo[t][d]);
                    rNode.mHi[d] = std::max(rNode.mHi[d], rHi[t][d]);
                }
            }
            continue;
        }
        const SurfaceBvhNode& rLeft = bvh.mNodes[i + 1];
        const SurfaceBvhNode& rRight = bvh.mNodes[static_cast<std::size_t>(rNode.mRight)];
        for (std::size_t d = 0; d < 3; ++d) {
            rNode.mLo[d] = std::min(rLeft.mLo[d], rRight.mLo[d]);
            rNode.mHi[d] = std::max(rLeft.mHi[d], rRight.mHi[d]);
        }
    }
    return bvh;
}

}  // namespace

DistanceQuery build_distance_query(const TriangleSoup& rSoup,
                                   const SurfaceDistanceOptions& rOptions) {
    if (rSoup.NumTriangles() == 0)
        throw std::invalid_argument(std::string(kSdPrefix) +
                                    "the surface has no triangles to measure against");
    return build_distance_query_from_runs(rSoup, rOptions, surface_edge_runs(rSoup));
}

DistanceQuery build_distance_query_from_runs(const TriangleSoup& rSoup,
                                             const SurfaceDistanceOptions& rOptions,
                                             const SurfaceEdgeRuns& rRuns) {
    const std::size_t ntri = rSoup.NumTriangles();
    if (ntri == 0)
        throw std::invalid_argument(std::string(kSdPrefix) +
                                    "the surface has no triangles to measure against");

    DistanceQuery q;
    q.mpSoup = &rSoup;

    // Every triangle's bounding box, in parallel: the tree is built over them.
    std::vector<Vec3> tri_lo(ntri);
    std::vector<Vec3> tri_hi(ntri);
    parallel_for(ntri, [&](std::size_t t) {
        Vec3 tlo = rSoup.mCorners[t * 3];
        Vec3 thi = tlo;
        for (std::size_t i = 1; i < 3; ++i)
            for (std::size_t k = 0; k < 3; ++k) {
                const double v = rSoup.mCorners[t * 3 + i][k];
                tlo[k] = tlo[k] < v ? tlo[k] : v;
                thi[k] = thi[k] > v ? thi[k] : v;
            }
        tri_lo[t] = tlo;
        tri_hi[t] = thi;
    });
    // `rOptions.mGridCellSize` sized the bucket grid this tree replaced
    // (ABI 40) and is ignored.
    q.mBvh = sd_build_bvh(tri_lo, tri_hi);

    // Face normals, then the vertex and edge tables. Every sum runs in
    // ascending (triangle, corner) order: summing unit normals in a different
    // order changes the last bits, and a last-bit change can flip the sign of a
    // query point sitting almost exactly on the surface.
    q.mFaceNormal = soup_face_normals(rSoup);
    q.mVertexNormal = accumulate_vertex_normals(rSoup, q.mFaceNormal, rOptions.mWeight);
    // Edge normals: the (triangle, corner) records grouped by edge with a
    // parallel sort (slot_runs.hpp), each run's unit normals summed in its
    // slot order -- ascending (triangle, corner), the serial order -- from
    // zero. A degenerate triangle contributes nothing, and an edge only
    // degenerate triangles use has no normal (-1), as it had no map entry.
    std::vector<Vec3> unit(ntri);
    std::vector<std::uint8_t> has_unit(ntri, 0);
    parallel_for(ntri, [&](std::size_t t) {
        const Vec3 n = q.mFaceNormal[t];
        const double len = vec3_norm(n);
        if (!(len > 0.0))
            return;  // degenerate: no direction to contribute
        unit[t] = vec3_scale(n, 1.0 / len);
        has_unit[t] = 1;
    });
    const SlotRuns& runs = rRuns.mRuns;
    std::vector<std::uint8_t> used(runs.NumRuns(), 0);
    parallel_for(runs.NumRuns(), [&](std::size_t r) {
        for (const std::uint64_t* p = runs.Begin(r); p != runs.End(r) && !used[r]; ++p)
            used[r] = has_unit[*p / 3];
    });
    std::vector<std::int64_t> index_of_run(runs.NumRuns());
    const std::int64_t nedges =
        parallel_exclusive_scan(used.data(), runs.NumRuns(), index_of_run.data(), std::int64_t{0});
    q.mEdgeNormals.assign(static_cast<std::size_t>(nedges), Vec3{0.0, 0.0, 0.0});
    q.mEdgeOfCorner.assign(ntri * 3, -1);
    parallel_for(runs.NumRuns(), [&](std::size_t r) {
        if (!used[r])
            return;
        const std::int64_t id = index_of_run[r];
        Vec3 e{0.0, 0.0, 0.0};
        for (const std::uint64_t* p = runs.Begin(r); p != runs.End(r); ++p) {
            if (has_unit[*p / 3])
                e = vec3_add(e, unit[*p / 3]);
            q.mEdgeOfCorner[*p] = id;
        }
        q.mEdgeNormals[static_cast<std::size_t>(id)] = e;
    });
    return q;
}

namespace {

// The nearest-triangle search shared by query_distances,
// query_closest_points and query_surface_projections, so the three cannot
// disagree about which triangle is nearest.
//
// It answers the minimum, under the total order (distance squared, triangle
// id), over every triangle whose distance is at most the best found -- the
// global minimum -- so which candidates it visits, and in what order, is
// unobservable. Depth first through the tree, the nearer child first; a node
// is opened only while its box's lower bound does not exceed the best
// distance. The bound gives up a slack for rounding and the cut is strict, so
// no triangle that could win or tie is ever pruned. Until ABI 40 this walked
// a bucket grid's Chebyshev shells outwards from the query's own bucket,
// O((d / cell)^3) probes for a point at distance d, which for a point deep
// inside a volume is of the order of the triangle count (roadmap 5.3). @p
// Limit2 is a squared distance past which a hit is of no use to the caller (a
// band), +inf for none.
struct SdNearestTriangle {
    std::int64_t mTri = -1;
    PointTriangleHit mHit;
};

#ifdef MESHIOPLUSPLUS_SD_COUNT_PROBES
// Probe counts for the benchmarks page (roadmap 5.3), compiled in only on
// request: queries, nodes opened, leaves opened, triangles measured.
std::atomic<std::uint64_t> g_sd_queries{0}, g_sd_nodes{0}, g_sd_leaves{0}, g_sd_triangles{0};
void sd_report_probes(const char* pWhat) {
    log::info(
        "surface distance: {}: {} queries opened {} nodes and {} leaves and measured {} "
        "triangles",
        pWhat, g_sd_queries.exchange(0), g_sd_nodes.exchange(0), g_sd_leaves.exchange(0),
        g_sd_triangles.exchange(0));
}
#define MIO_SD_COUNT(counter, n) counter.fetch_add(n, std::memory_order_relaxed)
#else
#define MIO_SD_COUNT(counter, n) ((void)0)
#endif

// The squared distance from @p rQ to box [rLo, rHi], less a slack on each
// axis's gap of 2^-40 of the magnitudes involved: far more than the rounding
// of the gap and of the triangle distances it is compared with. A non-finite
// gap bounds nothing.
double sd_box_lb2(const Vec3& rLo, const Vec3& rHi, const Vec3& rQ) {
    double sum = 0.0;
    for (std::size_t a = 0; a < 3; ++a) {
        double gap = rLo[a] - rQ[a] > rQ[a] - rHi[a] ? rLo[a] - rQ[a] : rQ[a] - rHi[a];
        if (!(gap > 0.0) || !std::isfinite(gap))
            continue;
        gap -= (std::fabs(rQ[a]) + std::fabs(rLo[a]) + std::fabs(rHi[a])) * 0x1p-40;
        if (gap > 0.0)
            sum += gap * gap;
    }
    return sum;
}

// The walk, for every triangle @p rAccept accepts (a filter only drops
// candidates, so every pruning bound stays valid).
template <class Accept>
SdNearestTriangle sd_nearest_triangle_if(const DistanceQuery& rQuery, const TriangleSoup& rSoup,
                                         const Vec3& rQueryPoint, double Limit2,
                                         const Accept& rAccept) {
    double best_d2 = std::numeric_limits<double>::infinity();
    std::int64_t best_tri = -1;
    PointTriangleHit best_hit;
    const std::vector<SurfaceBvhNode>& nodes = rQuery.mBvh.mNodes;
    const std::vector<std::int64_t>& order = rQuery.mBvh.mOrder;
    if (nodes.empty())
        return {best_tri, best_hit};
    MIO_SD_COUNT(g_sd_queries, 1);
    // (node, its bound); the tree's depth is log2 of the leaf count, far
    // below the stack's size.
    std::pair<std::size_t, double> stack[128];
    int top = 0;
    const double root = sd_box_lb2(nodes[0].mLo, nodes[0].mHi, rQueryPoint);
    if (!(root > Limit2))
        stack[top++] = {0, root};
    while (top > 0) {
        const auto [index, lb2] = stack[--top];
        if (lb2 > best_d2)
            continue;
        const SurfaceBvhNode& rNode = nodes[index];
        if (rNode.mRight < 0) {
            MIO_SD_COUNT(g_sd_leaves, 1);
            MIO_SD_COUNT(g_sd_triangles, static_cast<std::uint64_t>(rNode.mCount));
            for (std::int64_t k = 0; k < rNode.mCount; ++k) {
                const std::int64_t t = order[static_cast<std::size_t>(rNode.mBegin + k)];
                if (!rAccept(t))
                    continue;
                const std::size_t ti = static_cast<std::size_t>(t);
                const PointTriangleHit hit = closest_point_on_triangle(
                    rQueryPoint, rSoup.mCorners[ti * 3 + 0], rSoup.mCorners[ti * 3 + 1],
                    rSoup.mCorners[ti * 3 + 2]);
                // The total order that makes the accelerator unobservable:
                // distance first, then the triangle id, so two equidistant
                // triangles always resolve the same way whatever node was
                // opened first.
                if (hit.mDistanceSq < best_d2 || (hit.mDistanceSq == best_d2 && t < best_tri)) {
                    best_d2 = hit.mDistanceSq;
                    best_tri = t;
                    best_hit = hit;
                }
            }
            continue;
        }
        MIO_SD_COUNT(g_sd_nodes, 1);
        const std::size_t left = index + 1;
        const auto right = static_cast<std::size_t>(rNode.mRight);
        const double lb_left = sd_box_lb2(nodes[left].mLo, nodes[left].mHi, rQueryPoint);
        const double lb_right = sd_box_lb2(nodes[right].mLo, nodes[right].mHi, rQueryPoint);
        const double cut = best_d2 < Limit2 ? best_d2 : Limit2;
        // The nearer child is pushed last, so it is opened first.
        const bool left_first = lb_left <= lb_right;
        const std::size_t far = left_first ? right : left;
        const std::size_t near = left_first ? left : right;
        const double lb_far = left_first ? lb_right : lb_left;
        const double lb_near = left_first ? lb_left : lb_right;
        if (!(lb_far > cut))
            stack[top++] = {far, lb_far};
        if (!(lb_near > cut))
            stack[top++] = {near, lb_near};
    }
    // A cut at Limit2 still lets a farther leaf's triangle win; the callers
    // that pass a finite limit check the distance themselves.
    return {best_tri, best_hit};
}

struct SdAcceptAll {
    constexpr bool operator()(std::int64_t) const { return true; }
};

SdNearestTriangle sd_nearest_triangle(const DistanceQuery& rQuery, const TriangleSoup& rSoup,
                                      const Vec3& rQueryPoint, double Limit2) {
    return sd_nearest_triangle_if(rQuery, rSoup, rQueryPoint, Limit2, SdAcceptAll{});
}

// The pseudonormal of the FEATURE a hit landed on, not of the nearest
// triangle: using the triangle's own normal is right on convex geometry and
// wrong on the concave side of every crease. Hoisted verbatim out of
// query_distances (a pure refactor -- that function's own suite is the
// regression guard) so query_surface_projections reads the same tables the
// same way. Returned unnormalized.
Vec3 sd_feature_normal(const DistanceQuery& rQuery, const TriangleSoup& rSoup, std::int64_t Tri,
                       TriangleFeature Feature) {
    const std::size_t ti = static_cast<std::size_t>(Tri);
    const std::array<std::int64_t, 3>& v = rSoup.mVertices[ti];
    Vec3 normal = rQuery.mFaceNormal[ti];
    switch (Feature) {
        case TriangleFeature::VertexA:
            normal = rQuery.mVertexNormal[static_cast<std::size_t>(v[0])];
            break;
        case TriangleFeature::VertexB:
            normal = rQuery.mVertexNormal[static_cast<std::size_t>(v[1])];
            break;
        case TriangleFeature::VertexC:
            normal = rQuery.mVertexNormal[static_cast<std::size_t>(v[2])];
            break;
        case TriangleFeature::EdgeAB:
        case TriangleFeature::EdgeBC:
        case TriangleFeature::EdgeCA: {
            const std::size_t e = Feature == TriangleFeature::EdgeAB
                                      ? 0
                                      : (Feature == TriangleFeature::EdgeBC ? 1 : 2);
            const std::int64_t id = rQuery.mEdgeOfCorner[ti * 3 + e];
            if (id >= 0)
                normal = rQuery.mEdgeNormals[static_cast<std::size_t>(id)];
            break;
        }
        case TriangleFeature::Face:
        default:
            break;
    }
    return normal;
}

}  // namespace

std::vector<DistanceHit> query_distances(const DistanceQuery& rQuery,
                                         const std::vector<Vec3>& rPoints,
                                         const SurfaceDistanceOptions& rOptions) {
    const TriangleSoup& soup = *rQuery.mpSoup;
    const std::size_t ntri = soup.NumTriangles();
    const std::size_t n = rPoints.size();
    std::vector<DistanceHit> out(n);

    if (rOptions.mSign == SdfSign::WindingNumber && rOptions.mMaxWindingWork > 0.0) {
        const double work = static_cast<double>(n) * static_cast<double>(ntri);
        if (work > rOptions.mMaxWindingWork)
            throw std::invalid_argument(
                std::string(kSdPrefix) + "sign='winding-number' is O(triangles) per query and " +
                std::to_string(n) + " queries x " + std::to_string(ntri) + " triangles exceeds " +
                "max_winding_work (raise it, use a band, or use sign='pseudonormal')");
    }

    const double band = rOptions.mBand;
    const bool banded = band > 0.0;
    const double band2 = band * band;
    // With a band there is no point opening a block beyond it: a hit found
    // there would be clamped anyway. A hit within the band is the global
    // nearest either way.
    const double limit2 = banded ? band2 : std::numeric_limits<double>::infinity();

    parallel_for(n, [&](std::size_t p) {
        const Vec3& query = rPoints[p];
        const SdNearestTriangle found = sd_nearest_triangle(rQuery, soup, query, limit2);
        const std::int64_t best_tri = found.mTri;
        const PointTriangleHit& best_hit = found.mHit;
        const double best_d2 = best_hit.mDistanceSq;

        DistanceHit& res = out[p];
        if (best_tri < 0) {
            // Nothing within reach: only possible under a band.
            res.mSignedDistance = band;
            res.mSourceCell = -1;
            res.mInBand = false;
            return;
        }
        const double dist = std::sqrt(best_d2);
        if (banded && best_d2 > band2) {
            res.mSignedDistance = band;
            res.mSourceCell = -1;
            res.mInBand = false;
            return;
        }
        res.mSourceCell = soup.mSourceCell[static_cast<std::size_t>(best_tri)];
        res.mInBand = true;

        if (rOptions.mSign == SdfSign::Unsigned) {
            res.mSignedDistance = dist;
            return;
        }
        if (rOptions.mSign == SdfSign::WindingNumber) {
            // Van Oosterom-Strackee solid angle, summed in ascending triangle
            // order. O(triangles) per query, which is why it is guarded above.
            double w = 0.0;
            for (std::size_t t = 0; t < ntri; ++t) {
                const Vec3 a = vec3_sub(soup.mCorners[t * 3 + 0], query);
                const Vec3 b = vec3_sub(soup.mCorners[t * 3 + 1], query);
                const Vec3 c = vec3_sub(soup.mCorners[t * 3 + 2], query);
                const double la = vec3_norm(a);
                const double lb = vec3_norm(b);
                const double lc = vec3_norm(c);
                const double num = triple_product(a, b, c);
                const double den =
                    la * lb * lc + vec3_dot(a, b) * lc + vec3_dot(b, c) * la + vec3_dot(c, a) * lb;
                w += 2.0 * std::atan2(num, den);
            }
            const bool inside = w / (4.0 * 3.14159265358979323846) > 0.5;
            res.mSignedDistance = inside ? -dist : dist;
            return;
        }

        // Pseudonormal: the normal of the nearest FEATURE, not of the nearest
        // triangle (sd_feature_normal says why).
        const Vec3 normal = sd_feature_normal(rQuery, soup, best_tri, best_hit.mFeature);
        const double side = vec3_dot(vec3_sub(query, best_hit.mPoint), normal);
        res.mSignedDistance = side < 0.0 ? -dist : dist;
    });

#ifdef MESHIOPLUSPLUS_SD_COUNT_PROBES
    sd_report_probes("query_distances");
#endif
    return out;
}

std::vector<ClosestPointHit> query_closest_points(const DistanceQuery& rQuery,
                                                  const std::vector<Vec3>& rPoints) {
    const TriangleSoup& soup = *rQuery.mpSoup;
    const std::size_t n = rPoints.size();
    std::vector<ClosestPointHit> out(n);

    parallel_for(n, [&](std::size_t p) {
        const SdNearestTriangle found =
            sd_nearest_triangle(rQuery, soup, rPoints[p], std::numeric_limits<double>::infinity());
        ClosestPointHit& res = out[p];
        if (found.mTri < 0) {
            res.mFound = false;
            return;
        }
        res.mFound = true;
        res.mPoint = found.mHit.mPoint;
        res.mDistance = std::sqrt(found.mHit.mDistanceSq);
        res.mSourceCell = soup.mSourceCell[static_cast<std::size_t>(found.mTri)];
    });

#ifdef MESHIOPLUSPLUS_SD_COUNT_PROBES
    sd_report_probes("query_closest_points");
#endif
    return out;
}

std::vector<SurfaceProjection> query_surface_projections(const DistanceQuery& rQuery,
                                                         const std::vector<Vec3>& rPoints) {
    const TriangleSoup& soup = *rQuery.mpSoup;
    const std::size_t n = rPoints.size();
    std::vector<SurfaceProjection> out(n);

    parallel_for(n, [&](std::size_t p) {
        const SdNearestTriangle found =
            sd_nearest_triangle(rQuery, soup, rPoints[p], std::numeric_limits<double>::infinity());
        SurfaceProjection& res = out[p];
        if (found.mTri < 0) {
            res.mFound = false;
            return;
        }
        res.mFound = true;
        res.mPoint = found.mHit.mPoint;
        res.mDistance = std::sqrt(found.mHit.mDistanceSq);
        res.mTriangle = found.mTri;
        res.mSourceCell = soup.mSourceCell[static_cast<std::size_t>(found.mTri)];
        res.mFeature = found.mHit.mFeature;
        res.mNormal = sd_feature_normal(rQuery, soup, found.mTri, found.mHit.mFeature);
    });

#ifdef MESHIOPLUSPLUS_SD_COUNT_PROBES
    sd_report_probes("query_surface_projections");
#endif
    return out;
}

std::vector<SurfaceProjection> query_surface_projections_filtered(const DistanceQuery& rQuery,
                                                                  const std::vector<Vec3>& rPoints,
                                                                  const TriangleFilter& rAccept,
                                                                  double MaxDistance) {
    const TriangleSoup& soup = *rQuery.mpSoup;
    const std::size_t n = rPoints.size();
    std::vector<SurfaceProjection> out(n);
    const double limit2 = MaxDistance < std::numeric_limits<double>::infinity()
                              ? MaxDistance * MaxDistance
                              : std::numeric_limits<double>::infinity();
    parallel_for(n, [&](std::size_t p) {
        const auto accept = [&](std::int64_t T) { return rAccept(p, T); };
        const SdNearestTriangle found =
            sd_nearest_triangle_if(rQuery, soup, rPoints[p], limit2, accept);
        SurfaceProjection& res = out[p];
        if (found.mTri < 0 || found.mHit.mDistanceSq > limit2) {
            res.mFound = false;
            return;
        }
        res.mFound = true;
        res.mPoint = found.mHit.mPoint;
        res.mDistance = std::sqrt(found.mHit.mDistanceSq);
        res.mTriangle = found.mTri;
        res.mSourceCell = soup.mSourceCell[static_cast<std::size_t>(found.mTri)];
        res.mFeature = found.mHit.mFeature;
        res.mNormal = sd_feature_normal(rQuery, soup, found.mTri, found.mHit.mFeature);
    });
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
