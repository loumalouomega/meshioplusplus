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

// segment. See operations/segment.hpp for the contract.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/segment.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "../detail/crease_edges.hpp"
#include "../detail/surface_rings.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kSegPrefix = "meshio++: segment: ";
constexpr double kSegPi = 3.14159265358979323846;

// Per global cell, whether the region selects it (all cells when unnamed).
std::vector<std::uint8_t> seg_selection(const Mesh& rMesh, const std::string& rRegion,
                                        std::size_t NumCells) {
    if (rRegion.empty())
        return std::vector<std::uint8_t>(NumCells, 1);
    const std::size_t r = rMesh.FindRegion(rRegion, RegionKind::Cell);
    if (r == Mesh::npos) {
        std::string names;
        for (const std::string& n : rMesh.RegionNames())
            names += (names.empty() ? "" : ", ") + n;
        throw std::invalid_argument(std::string(kSegPrefix) + "no cell region named '" + rRegion +
                                    "' (available: " + (names.empty() ? "none" : names) + ")");
    }
    std::vector<std::uint8_t> selected(NumCells, 0);
    const NDArray& entries = rMesh.Region(r).mEntries;
    for (std::size_t e = 0; e < entries.Size(); ++e) {
        const std::int64_t c = detail::read_int(entries, e);
        if (c >= 0 && static_cast<std::size_t>(c) < NumCells)
            selected[static_cast<std::size_t>(c)] = 1;
    }
    return selected;
}

// One use of an undirected edge by one face ring.
struct SegEdgeUse {
    std::int64_t mLo;
    std::int64_t mHi;
    std::int64_t mFace;
    bool mForward;
};

bool seg_use_less(const SegEdgeUse& rA, const SegEdgeUse& rB) {
    if (rA.mLo != rB.mLo)
        return rA.mLo < rB.mLo;
    if (rA.mHi != rB.mHi)
        return rA.mHi < rB.mHi;
    if (rA.mFace != rB.mFace)
        return rA.mFace < rB.mFace;
    return rA.mForward < rB.mForward;
}

// An edge shared by exactly two distinct faces.
struct SegPair {
    std::int64_t mA;  // the lower face
    std::int64_t mB;  // the higher face
    bool mInconsistent;
    bool mSmooth;
    double mAngleRad;  // NaN when either face is degenerate
};

}  // namespace

SegmentResult segment(const Mesh& rMesh, const SegmentOptions& rOptions) {
    if (!(rOptions.mFeatureAngleDeg >= 0.0 && rOptions.mFeatureAngleDeg <= 180.0))
        throw std::invalid_argument(std::string(kSegPrefix) +
                                    "the feature angle must lie in [0, 180] degrees");
    if (!(rOptions.mNormalToleranceDeg >= 0.0 && rOptions.mNormalToleranceDeg <= 180.0))
        throw std::invalid_argument(std::string(kSegPrefix) +
                                    "the normal tolerance must lie in [0, 180] degrees");
    if (!(rOptions.mCurvatureTolerance >= 0.0) || !std::isfinite(rOptions.mCurvatureTolerance))
        throw std::invalid_argument(std::string(kSegPrefix) +
                                    "the curvature tolerance must be finite and >= 0");
    if (rOptions.mMinPatchCells < 1)
        throw std::invalid_argument(std::string(kSegPrefix) + "min_patch_cells must be >= 1");
    if (rOptions.mPrefix.empty())
        throw std::invalid_argument(std::string(kSegPrefix) +
                                    "the region prefix must not be empty");

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t ncells = static_cast<std::size_t>(detail::total_cells(bases));
    const std::vector<std::uint8_t> selected = seg_selection(rMesh, rOptions.mRegion, ncells);

    const detail::SurfaceRings rings = detail::surface_rings(rMesh, selected, kSegPrefix);
    const std::vector<std::int64_t>& start = rings.mStart;
    const std::vector<std::int64_t>& nodes = rings.mNodes;
    const std::vector<std::int64_t>& cell_of = rings.mCell;
    const std::size_t nf = rings.NumFaces();

    // Flat coordinates padded to 3-D, unit normals and centroids per face.
    const std::vector<double> xyz = detail::surface_xyz(rMesh);
    std::vector<double> normals(nf * 3);
    std::vector<double> centroid(nf * 3);
    parallel_for(nf, [&](std::size_t f) {
        const std::int64_t* ring = nodes.data() + start[f];
        const std::size_t size = static_cast<std::size_t>(start[f + 1] - start[f]);
        detail::ring_unit_normal(xyz.data(), ring, size, &normals[f * 3]);
        double c[3] = {0.0, 0.0, 0.0};
        for (std::size_t k = 0; k < size; ++k)
            for (std::size_t d = 0; d < 3; ++d)
                c[d] += xyz[static_cast<std::size_t>(ring[k]) * 3 + d];
        for (std::size_t d = 0; d < 3; ++d)
            centroid[f * 3 + d] = c[d] / static_cast<double>(size);
    });

    // Every non-collapsed ring edge is one use; sorted, each run is one edge.
    std::vector<std::uint64_t> count(nf, 0);
    parallel_for(nf, [&](std::size_t f) {
        const std::int64_t b = start[f];
        const std::int64_t s = start[f + 1] - b;
        std::uint64_t c = 0;
        for (std::int64_t k = 0; k < s; ++k)
            c += nodes[static_cast<std::size_t>(b + k)] !=
                         nodes[static_cast<std::size_t>(b + (k + 1) % s)]
                     ? 1
                     : 0;
        count[f] = c;
    });
    std::vector<std::uint64_t> at(nf);
    const std::uint64_t nuses =
        parallel_exclusive_scan(count.data(), nf, at.data(), std::uint64_t{0});
    std::vector<SegEdgeUse> uses(nuses);
    parallel_for(nf, [&](std::size_t f) {
        const std::int64_t b = start[f];
        const std::int64_t s = start[f + 1] - b;
        std::uint64_t w = at[f];
        for (std::int64_t k = 0; k < s; ++k) {
            const std::int64_t u = nodes[static_cast<std::size_t>(b + k)];
            const std::int64_t v = nodes[static_cast<std::size_t>(b + (k + 1) % s)];
            if (u == v)
                continue;
            const std::int64_t face = static_cast<std::int64_t>(f);
            uses[w++] = u < v ? SegEdgeUse{u, v, face, true} : SegEdgeUse{v, u, face, false};
        }
    });
    parallel_sort(uses.begin(), uses.end(), seg_use_less);

    SegmentResult result;
    const double cos_feature = std::cos(rOptions.mFeatureAngleDeg * kSegPi / 180.0);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<SegPair> pairs;
    for (std::size_t b = 0; b < uses.size();) {
        std::size_t e = b + 1;
        while (e < uses.size() && uses[e].mLo == uses[b].mLo && uses[e].mHi == uses[b].mHi)
            ++e;
        const std::size_t u = e - b;
        if (u == 1) {
            ++result.mNumBoundaryEdges;
        } else if (u >= 3) {
            ++result.mNumNonManifoldEdges;
        } else if (uses[b].mFace != uses[b + 1].mFace) {
            SegPair p;
            p.mA = uses[b].mFace;
            p.mB = uses[b + 1].mFace;
            p.mInconsistent = uses[b].mForward == uses[b + 1].mForward;
            const double* nx = &normals[static_cast<std::size_t>(p.mA) * 3];
            const double* ny = &normals[static_cast<std::size_t>(p.mB) * 3];
            const bool zx = nx[0] == 0.0 && nx[1] == 0.0 && nx[2] == 0.0;
            const bool zy = ny[0] == 0.0 && ny[1] == 0.0 && ny[2] == 0.0;
            p.mSmooth = true;
            p.mAngleRad = nan;
            if (!zx && !zy) {
                double dot = nx[0] * ny[0] + nx[1] * ny[1] + nx[2] * ny[2];
                if (p.mInconsistent)
                    dot = -dot;
                p.mSmooth = !(dot < cos_feature);
                const double c = dot < -1.0 ? -1.0 : (dot > 1.0 ? 1.0 : dot);
                p.mAngleRad = std::acos(c);
            }
            if (!p.mSmooth)
                ++result.mNumFeatureEdges;
            pairs.push_back(p);
        }
        b = e;
    }

    // Per face, its neighbours across paired edges, ascending by face then edge.
    std::vector<std::int64_t> adj_start(nf + 1, 0);
    for (const SegPair& p : pairs) {
        ++adj_start[static_cast<std::size_t>(p.mA) + 1];
        ++adj_start[static_cast<std::size_t>(p.mB) + 1];
    }
    for (std::size_t f = 0; f < nf; ++f)
        adj_start[f + 1] += adj_start[f];
    std::vector<std::pair<std::int64_t, std::int64_t>> adj(static_cast<std::size_t>(adj_start[nf]));
    {
        std::vector<std::int64_t> fill(adj_start.begin(), adj_start.end() - 1);
        for (std::size_t e = 0; e < pairs.size(); ++e) {
            const std::int64_t ei = static_cast<std::int64_t>(e);
            adj[static_cast<std::size_t>(fill[static_cast<std::size_t>(pairs[e].mA)]++)] = {
                pairs[e].mB, ei};
            adj[static_cast<std::size_t>(fill[static_cast<std::size_t>(pairs[e].mB)]++)] = {
                pairs[e].mA, ei};
        }
        parallel_for(nf, [&](std::size_t f) {
            std::sort(adj.begin() + adj_start[f], adj.begin() + adj_start[f + 1]);
        });
    }

    // The curvature estimate, only when it is asked for.
    const bool use_curv = rOptions.mCurvatureTolerance > 0.0;
    std::vector<double> curv;
    if (use_curv) {
        curv.assign(nf, 0.0);
        for (const SegPair& p : pairs) {
            if (!p.mSmooth || std::isnan(p.mAngleRad))
                continue;
            const double* ca = &centroid[static_cast<std::size_t>(p.mA) * 3];
            const double* cb = &centroid[static_cast<std::size_t>(p.mB) * 3];
            const double dx = ca[0] - cb[0], dy = ca[1] - cb[1], dz = ca[2] - cb[2];
            const double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (!(dist > 0.0))
                continue;
            const double k = p.mAngleRad / dist;
            double& ka = curv[static_cast<std::size_t>(p.mA)];
            double& kb = curv[static_cast<std::size_t>(p.mB)];
            ka = k > ka ? k : ka;
            kb = k > kb ? k : kb;
        }
    }

    // Region growing: a serial breadth-first walk from each unclaimed face.
    const bool use_normal = rOptions.mNormalToleranceDeg > 0.0;
    const double cos_normal = std::cos(rOptions.mNormalToleranceDeg * kSegPi / 180.0);
    std::vector<std::int64_t> patch(nf, -1);
    std::vector<std::uint8_t> flip(nf, 0);
    std::vector<std::int64_t> queue;
    queue.reserve(nf);
    std::int64_t npatches = 0;
    for (std::size_t seed = 0; seed < nf; ++seed) {
        if (patch[seed] != -1)
            continue;
        const std::int64_t p = npatches++;
        patch[seed] = p;
        flip[seed] = 0;
        const double* sn = &normals[seed * 3];
        // The patch's mean curvature estimate so far, in walk order.
        double curv_sum = use_curv ? curv[seed] : 0.0;
        double curv_count = 1.0;
        queue.clear();
        queue.push_back(static_cast<std::int64_t>(seed));
        for (std::size_t q = 0; q < queue.size(); ++q) {
            const std::size_t g = static_cast<std::size_t>(queue[q]);
            for (std::int64_t k = adj_start[g]; k < adj_start[g + 1]; ++k) {
                const std::size_t nb =
                    static_cast<std::size_t>(adj[static_cast<std::size_t>(k)].first);
                const SegPair& e =
                    pairs[static_cast<std::size_t>(adj[static_cast<std::size_t>(k)].second)];
                if (patch[nb] != -1 || !e.mSmooth)
                    continue;
                const std::uint8_t fl =
                    static_cast<std::uint8_t>(flip[g] ^ (e.mInconsistent ? 1 : 0));
                if (use_normal) {
                    const double* nn = &normals[nb * 3];
                    double dot = sn[0] * nn[0] + sn[1] * nn[1] + sn[2] * nn[2];
                    if (fl)
                        dot = -dot;
                    if (dot < cos_normal)
                        continue;
                }
                if (use_curv &&
                    std::fabs(curv[nb] - curv_sum / curv_count) > rOptions.mCurvatureTolerance)
                    continue;
                if (use_curv) {
                    curv_sum += curv[nb];
                    curv_count += 1.0;
                }
                patch[nb] = p;
                flip[nb] = fl;
                queue.push_back(static_cast<std::int64_t>(nb));
            }
        }
    }

    // Merge the patches below the minimum size, in one ascending pass.
    std::vector<std::int64_t> label(static_cast<std::size_t>(npatches));
    for (std::int64_t p = 0; p < npatches; ++p)
        label[static_cast<std::size_t>(p)] = p;
    if (rOptions.mMinPatchCells > 1 && npatches > 1) {
        std::vector<std::int64_t> size(static_cast<std::size_t>(npatches), 0);
        for (std::size_t f = 0; f < nf; ++f)
            ++size[static_cast<std::size_t>(patch[f])];
        std::vector<std::map<std::int64_t, std::int64_t>> shared(
            static_cast<std::size_t>(npatches));
        for (const SegPair& e : pairs) {
            const std::int64_t a = patch[static_cast<std::size_t>(e.mA)];
            const std::int64_t b = patch[static_cast<std::size_t>(e.mB)];
            if (a == b)
                continue;
            ++shared[static_cast<std::size_t>(a)][b];
            ++shared[static_cast<std::size_t>(b)][a];
        }
        for (std::int64_t p = 0; p < npatches; ++p) {
            const std::size_t ps = static_cast<std::size_t>(p);
            if (label[ps] != p || size[ps] >= rOptions.mMinPatchCells || shared[ps].empty())
                continue;
            std::int64_t best = -1, best_count = -1;
            for (const auto& kv : shared[ps])
                if (kv.second > best_count) {
                    best = kv.first;
                    best_count = kv.second;
                }
            const std::size_t bs = static_cast<std::size_t>(best);
            for (const auto& kv : shared[ps]) {
                const std::size_t rs = static_cast<std::size_t>(kv.first);
                shared[rs].erase(p);
                if (kv.first == best)
                    continue;
                shared[bs][kv.first] += kv.second;
                shared[rs][best] += kv.second;
            }
            shared[ps].clear();
            size[bs] += size[ps];
            size[ps] = 0;
            label[ps] = best;
            ++result.mNumMerged;
        }
        // Resolve chains of merges to their surviving patch.
        for (std::int64_t p = 0; p < npatches; ++p) {
            std::int64_t r = p;
            while (label[static_cast<std::size_t>(r)] != r)
                r = label[static_cast<std::size_t>(r)];
            label[static_cast<std::size_t>(p)] = r;
        }
    }

    // Final ids in the order of each patch's lowest face.
    std::vector<std::int64_t> final_id(static_cast<std::size_t>(npatches), -1);
    std::int64_t nfinal = 0;
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t r = static_cast<std::size_t>(label[static_cast<std::size_t>(patch[f])]);
        if (final_id[r] == -1)
            final_id[r] = nfinal++;
        patch[f] = final_id[r];
    }
    result.mNumPatches = nfinal;

    // Per global cell, its patch; then the regions and the adjacency table.
    std::vector<std::int64_t> cell_patch(ncells, -1);
    for (std::size_t f = 0; f < nf; ++f)
        cell_patch[static_cast<std::size_t>(cell_of[f])] = patch[f];

    std::map<std::pair<std::int64_t, std::int64_t>, std::int64_t> adjacency;
    for (const SegPair& e : pairs) {
        std::int64_t a = patch[static_cast<std::size_t>(e.mA)];
        std::int64_t b = patch[static_cast<std::size_t>(e.mB)];
        if (a == b)
            continue;
        if (a > b)
            std::swap(a, b);
        ++adjacency[{a, b}];
    }

    Mesh& out = result.mMesh;
    out = detail::clone_mesh(rMesh,
                             [](DataLocation, const std::string&, std::string&) { return true; });

    std::vector<NDArray> blocks;
    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const std::size_t base = static_cast<std::size_t>(bases[bi]);
        const std::size_t nc = rMesh.Cells(bi).NumCells();
        NDArray a = NDArray::Uninit(DType::Int64, {nc});
        std::copy(cell_patch.begin() + static_cast<std::ptrdiff_t>(base),
                  cell_patch.begin() + static_cast<std::ptrdiff_t>(base + nc),
                  a.As<std::int64_t>());
        blocks.push_back(std::move(a));
    }
    out.AddCellData(kSegmentPatchName, std::move(blocks));

    NDArray table = NDArray::Uninit(DType::Int64, {adjacency.size(), std::size_t{3}});
    {
        std::int64_t* t = table.As<std::int64_t>();
        std::size_t r = 0;
        for (const auto& kv : adjacency) {
            t[r * 3 + 0] = kv.first.first;
            t[r * 3 + 1] = kv.first.second;
            t[r * 3 + 2] = kv.second;
            ++r;
        }
    }
    out.AddFieldData(kSegmentAdjacencyName, std::move(table));

    // Replace an earlier run's `<prefix>:<integer>` cell regions.
    const std::string head = rOptions.mPrefix + ":";
    for (std::size_t i = out.NumRegions(); i-- > 0;) {
        const meshioplusplus::Region& r = out.Region(i);
        if (r.mKind != RegionKind::Cell || r.mName.size() <= head.size() ||
            r.mName.compare(0, head.size(), head) != 0)
            continue;
        const bool digits = std::all_of(r.mName.begin() + static_cast<std::ptrdiff_t>(head.size()),
                                        r.mName.end(), [](char c) { return c >= '0' && c <= '9'; });
        if (digits)
            out.RemoveRegion(i);
    }
    std::vector<std::vector<std::int64_t>> members(static_cast<std::size_t>(nfinal));
    for (std::size_t c = 0; c < ncells; ++c)
        if (cell_patch[c] >= 0)
            members[static_cast<std::size_t>(cell_patch[c])].push_back(
                static_cast<std::int64_t>(c));
    for (std::int64_t p = 0; p < nfinal; ++p) {
        const std::vector<std::int64_t>& m = members[static_cast<std::size_t>(p)];
        NDArray entries = NDArray::Uninit(DType::Int64, {m.size()});
        std::copy(m.begin(), m.end(), entries.As<std::int64_t>());
        out.AddRegion(meshioplusplus::Region(head + std::to_string(p), RegionKind::Cell, 2, -1,
                                             std::move(entries)));
    }
    return result;
}

}  // namespace meshioplusplus
