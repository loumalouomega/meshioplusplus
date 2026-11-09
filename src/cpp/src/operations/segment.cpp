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
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/segment.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"

// Project includes (private, not installed)
#include "../detail/cell_region_mask.hpp"
#include "../detail/patch_output.hpp"
#include "../detail/primitive_geom.hpp"
#include "../detail/primitive_grow.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kSegPrefix = "meshio++: segment: ";
constexpr double kSegPi = 3.14159265358979323846;

// The primitive method's patches (steps 1 and 2 of segment.hpp's "Primitive
// method"): sets rPatch for the faces a primitive claims, returns the count.
std::int64_t seg_primitive_patches(detail::PrimitiveGrowInput& rIn, std::size_t NumFaces,
                                   std::vector<std::int64_t>& rPatch) {
    const std::size_t nf = NumFaces;
    const std::vector<std::vector<std::int64_t>> areas = detail::primitive_smooth_areas(rIn, nf);

    std::int64_t npatches = 0;
    std::vector<std::uint8_t> mark(nf, 0), tried(nf, 0);
    constexpr std::uint32_t kAll = (1u << 1) | (1u << 2) | (1u << 3) | (1u << 4) | (1u << 5);
    constexpr std::uint32_t kCurved = (1u << 2) | (1u << 3) | (1u << 4) | (1u << 5);
    for (const std::vector<std::int64_t>& area : areas) {
        // 1. One primitive for the whole smooth area.
        const detail::PrimitiveFit whole = detail::primitive_grow_fit(rIn, area, kAll);
        if (whole.mType > 0) {
            for (std::int64_t f : area)
                rPatch[static_cast<std::size_t>(f)] = npatches;
            ++npatches;
            continue;
        }
        // 2. Grown regions, largest seed first; claimed for the area only once
        // it is clear the primitives explain it (below).
        const double area_total = detail::primitive_grow_area_of(rIn, area);
        const std::int64_t first_patch = npatches;
        std::vector<std::pair<std::vector<std::int64_t>, bool>> grown;  // faces, curved
        std::vector<std::int64_t> seeds = area;
        std::stable_sort(seeds.begin(), seeds.end(), [&](std::int64_t A, std::int64_t B) {
            return rIn.mArea[static_cast<std::size_t>(A)] > rIn.mArea[static_cast<std::size_t>(B)];
        });
        for (std::int64_t s : seeds) {
            const std::size_t ss = static_cast<std::size_t>(s);
            if (rPatch[ss] != -1 || tried[ss])
                continue;
            const detail::PrimitiveGrowRegion region =
                detail::primitive_grow_from_seed(rIn, s, true, kCurved, area_total, 0.02,
                                                 detail::kPrimitiveGrowExactShare, rPatch, mark);
            if (region.mFaces.empty()) {
                // Neither surface holds here: no cell of this neighbourhood
                // seeds again (a freeform area would cost a fit per cell);
                // another seed's region may still grow over them.
                for (std::int64_t f : region.mHood)
                    tried[static_cast<std::size_t>(f)] = 1;
                continue;
            }
            const std::vector<std::int64_t>* chosen = &region.mFaces;
            // Claim the faces now, so later seeds grow around them.
            for (std::int64_t f : *chosen)
                rPatch[static_cast<std::size_t>(f)] = npatches;
            ++npatches;
            grown.emplace_back(*chosen, region.mCurved);
        }
        // A freeform surface fits small spheres and tori anywhere: a curved
        // region under 1% of its smooth area is dropped, and an area the kept
        // regions explain less than half of is left whole to the angle
        // growth, as the angle method segments it.
        double explained = 0.0;
        std::vector<std::uint8_t> keep(grown.size(), 1);
        for (std::size_t i = 0; i < grown.size(); ++i) {
            const double a = detail::primitive_grow_area_of(rIn, grown[i].first);
            if (grown[i].second && a < 0.01 * area_total)
                keep[i] = 0;
            else
                explained += a;
        }
        const bool trusted = explained >= 0.5 * area_total;
        npatches = first_patch;
        for (std::size_t i = 0; i < grown.size(); ++i) {
            const bool kept = trusted && keep[i];
            for (std::int64_t f : grown[i].first)
                rPatch[static_cast<std::size_t>(f)] = kept ? npatches : -1;
            if (kept)
                ++npatches;
        }
    }
    return npatches;
}

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
    if (rOptions.mMethod != SegmentMethod::Angle && rOptions.mMethod != SegmentMethod::Primitive)
        throw std::invalid_argument(std::string(kSegPrefix) +
                                    "the method must be 'angle' or 'primitive'");
    if (!(rOptions.mTolerance > 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kSegPrefix) +
                                    "the tolerance must be finite and > 0");

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t ncells = static_cast<std::size_t>(detail::total_cells(bases));
    const std::vector<std::uint8_t> selected =
        detail::cell_region_mask(rMesh, rOptions.mRegion, ncells, kSegPrefix);

    const detail::PrimitiveGrowSurface surf = detail::build_primitive_grow_surface(
        rMesh, selected, rOptions.mFeatureAngleDeg, kSegPrefix);
    const std::vector<std::int64_t>& start = surf.mRings.mStart;
    const std::vector<std::int64_t>& nodes = surf.mRings.mNodes;
    const std::vector<std::int64_t>& cell_of = surf.mRings.mCell;
    const std::size_t nf = surf.NumFaces();
    const std::vector<double>& xyz = surf.mXyz;
    const std::vector<double>& normals = surf.mNormals;
    const std::vector<double>& centroid = surf.mCentroid;
    const std::vector<detail::PrimitiveGrowPair>& pairs = surf.mPairs;
    const std::vector<std::int64_t>& adj_start = surf.mAdjStart;
    const std::vector<std::pair<std::int64_t, std::int64_t>>& adj = surf.mAdj;

    SegmentResult result;
    result.mNumBoundaryEdges = surf.mNumBoundaryEdges;
    result.mNumNonManifoldEdges = surf.mNumNonManifoldEdges;
    result.mNumFeatureEdges = surf.mNumFeatureEdges;

    // The curvature estimate, only when it is asked for.
    const bool use_curv = rOptions.mCurvatureTolerance > 0.0;
    std::vector<double> curv;
    if (use_curv) {
        curv.assign(nf, 0.0);
        for (const detail::PrimitiveGrowPair& p : pairs) {
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
    if (rOptions.mMethod == SegmentMethod::Primitive && nf > 0) {
        detail::PrimitiveGrowInput in = surf.Input();
        detail::prepare_primitive_grow(in, nf, rOptions.mTolerance);
        if (in.mTol > 0.0)
            npatches = seg_primitive_patches(in, nf, patch);
    }
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
                const detail::PrimitiveGrowPair& e =
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
        for (const detail::PrimitiveGrowPair& e : pairs) {
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
    for (const detail::PrimitiveGrowPair& e : pairs) {
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

    detail::patch_output_write(rMesh, out, bases, cell_patch, static_cast<std::size_t>(nfinal),
                               adjacency, rOptions.mPrefix);
    return result;
}

}  // namespace meshioplusplus
