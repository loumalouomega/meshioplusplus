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
#include "../detail/primitive_geom.hpp"
#include "../detail/surface_rings.hpp"

namespace meshioplusplus {
namespace {

constexpr const char* kSegPrefix = "meshio++: segment: ";
// The primitive method's largest RMS of a grown region, as a share of the
// tolerance.
constexpr double kSegExactShare = 0.02;
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

// The primitive method's view of the surface: per face its ring, area,
// normal and centroid, its smooth neighbours, and the smooth area it lies in.
struct SegPrimitiveInput {
    const std::vector<std::int64_t>* mStart;
    const std::vector<std::int64_t>* mNodes;
    const std::vector<double>* mXyz;
    const std::vector<double>* mNormals;
    const std::vector<double>* mCentroid;
    const std::vector<std::int64_t>* mAdjStart;
    const std::vector<std::pair<std::int64_t, std::int64_t>>* mAdj;
    const std::vector<SegPair>* mPairs;
    std::vector<double> mArea;
    std::vector<std::int64_t> mSmoothArea;
    double mTol = 0.0;
    double mGrowTol = 0.0;  // how far a corner may lie from a growing region's primitive
    double mDiag = 0.0;
};

// The samples fit_primitives fits for a set of faces: their distinct corners,
// normals and centroids.
void seg_samples(const SegPrimitiveInput& rIn, const std::vector<std::int64_t>& rFaces,
                 std::vector<double>& rPts, std::vector<double>& rNormals,
                 std::vector<double>& rCentres) {
    std::vector<std::int64_t> verts;
    for (std::int64_t f : rFaces)
        for (std::int64_t k = (*rIn.mStart)[static_cast<std::size_t>(f)];
             k < (*rIn.mStart)[static_cast<std::size_t>(f) + 1]; ++k)
            verts.push_back((*rIn.mNodes)[static_cast<std::size_t>(k)]);
    std::sort(verts.begin(), verts.end());
    verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
    rPts.resize(verts.size() * 3);
    for (std::size_t i = 0; i < verts.size(); ++i)
        for (std::size_t d = 0; d < 3; ++d)
            rPts[i * 3 + d] = (*rIn.mXyz)[static_cast<std::size_t>(verts[i]) * 3 + d];
    rNormals.clear();
    rCentres.clear();
    for (std::int64_t f : rFaces)
        for (std::size_t d = 0; d < 3; ++d) {
            rNormals.push_back((*rIn.mNormals)[static_cast<std::size_t>(f) * 3 + d]);
            rCentres.push_back((*rIn.mCentroid)[static_cast<std::size_t>(f) * 3 + d]);
        }
}

// A fit accepted with every point within the tolerance, and not a curved
// surface so flat it stands in for a plane.
detail::PrimitiveFit seg_fit(const SegPrimitiveInput& rIn, const std::vector<std::int64_t>& rFaces,
                             std::uint32_t Types) {
    std::vector<double> pts, nrm, ctr;
    seg_samples(rIn, rFaces, pts, nrm, ctr);
    detail::PrimitiveFit fit =
        detail::fit_primitive_samples(pts, nrm, ctr, rIn.mTol, Types, 1.0, 50);
    if (fit.mType >= 2) {
        bool sane = true;
        for (int k = 0; k < detail::kPrimitiveParams; ++k)
            sane = sane && std::isfinite(fit.mParams[k]);
        if (fit.mType == 4)
            sane = sane && fit.mParams[8] < 85.0;
        else
            sane = sane && fit.mParams[6] <= rIn.mDiag && fit.mParams[7] <= rIn.mDiag;
        if (!sane)
            fit.mType = 0;
    }
    return fit;
}

bool seg_face_fits(const SegPrimitiveInput& rIn, std::int64_t F, const detail::PrimitiveFit& rFit,
                   double Within) {
    for (std::int64_t k = (*rIn.mStart)[static_cast<std::size_t>(F)];
         k < (*rIn.mStart)[static_cast<std::size_t>(F) + 1]; ++k) {
        const double* p =
            &(*rIn.mXyz)[static_cast<std::size_t>((*rIn.mNodes)[static_cast<std::size_t>(k)]) * 3];
        if (!(detail::primitive_point_distance(rFit.mType, rFit.mParams, p) <= Within))
            return false;
    }
    return true;
}

// Grow a region from rFaces across smooth edges, within the faces' smooth
// area, over faces not yet claimed, while every corner lies within the
// growth tolerance of the primitive. The primitive was fitted to FittedOn faces; it
// is refitted (same type) each time the region doubles past that, never on
// fewer faces, where a curved fit is ill-conditioned.
std::vector<std::int64_t> seg_grow(const SegPrimitiveInput& rIn, std::vector<std::int64_t> Faces,
                                   detail::PrimitiveFit Fit, std::size_t FittedOn,
                                   const std::vector<std::int64_t>& rPatch,
                                   std::vector<std::uint8_t>& rMark) {
    for (std::int64_t f : Faces)
        rMark[static_cast<std::size_t>(f)] = 1;
    const std::int64_t area = rIn.mSmoothArea[static_cast<std::size_t>(Faces.front())];
    std::size_t fitted_at = std::max(Faces.size(), FittedOn);
    for (std::size_t q = 0; q < Faces.size(); ++q) {
        const std::size_t g = static_cast<std::size_t>(Faces[q]);
        for (std::int64_t k = (*rIn.mAdjStart)[g]; k < (*rIn.mAdjStart)[g + 1]; ++k) {
            const auto& a = (*rIn.mAdj)[static_cast<std::size_t>(k)];
            const std::size_t nb = static_cast<std::size_t>(a.first);
            if (rMark[nb] || rPatch[nb] != -1 ||
                !(*rIn.mPairs)[static_cast<std::size_t>(a.second)].mSmooth ||
                rIn.mSmoothArea[nb] != area || !seg_face_fits(rIn, a.first, Fit, rIn.mGrowTol))
                continue;
            rMark[nb] = 1;
            Faces.push_back(a.first);
        }
        if (Faces.size() >= 2 * fitted_at) {
            fitted_at = Faces.size();
            const detail::PrimitiveFit refit =
                seg_fit(rIn, Faces, 1u << static_cast<unsigned>(Fit.mType));
            if (refit.mType == Fit.mType)
                Fit = refit;
        }
    }
    for (std::int64_t f : Faces)
        rMark[static_cast<std::size_t>(f)] = 0;
    return Faces;
}

// Trim a grown region to the cells its own refit explains within Within,
// keeping the part edge-connected to its first cell, three times over: a
// plane grown to the tolerance takes the first strip of a fillet tangent to
// it, which tilts its fit; the refit without it leaves that strip out.
void seg_trim(const SegPrimitiveInput& rIn, std::vector<std::int64_t>& rRegion, std::uint32_t Types,
              double Within, std::vector<std::uint8_t>& rMark) {
    for (int pass = 0; pass < 3 && !rRegion.empty(); ++pass) {
        const detail::PrimitiveFit fit = seg_fit(rIn, rRegion, Types);
        if (fit.mType == 0) {
            rRegion.clear();
            return;
        }
        std::size_t kept = 0;
        for (std::int64_t f : rRegion) {
            bool on = true;
            for (std::int64_t k = (*rIn.mStart)[static_cast<std::size_t>(f)];
                 k < (*rIn.mStart)[static_cast<std::size_t>(f) + 1] && on; ++k)
                on = detail::primitive_point_distance(
                         fit.mType, fit.mParams,
                         &(*rIn.mXyz)[static_cast<std::size_t>(
                                          (*rIn.mNodes)[static_cast<std::size_t>(k)]) *
                                      3]) <= Within;
            rMark[static_cast<std::size_t>(f)] = on ? 1 : 0;
            kept += on ? 1 : 0;
        }
        if (kept == rRegion.size()) {
            for (std::int64_t f : rRegion)
                rMark[static_cast<std::size_t>(f)] = 0;
            return;
        }
        // The part connected to the seed (the first cell), walk order kept.
        std::vector<std::int64_t> part;
        if (rMark[static_cast<std::size_t>(rRegion.front())]) {
            part.push_back(rRegion.front());
            rMark[static_cast<std::size_t>(rRegion.front())] = 2;
            for (std::size_t q = 0; q < part.size(); ++q) {
                const std::size_t g = static_cast<std::size_t>(part[q]);
                for (std::int64_t k = (*rIn.mAdjStart)[g]; k < (*rIn.mAdjStart)[g + 1]; ++k) {
                    const std::size_t nb =
                        static_cast<std::size_t>((*rIn.mAdj)[static_cast<std::size_t>(k)].first);
                    if (rMark[nb] != 1)
                        continue;
                    rMark[nb] = 2;
                    part.push_back(static_cast<std::int64_t>(nb));
                }
            }
        }
        for (std::int64_t f : rRegion)
            rMark[static_cast<std::size_t>(f)] = 0;
        rRegion.swap(part);
    }
}

double seg_area_of(const SegPrimitiveInput& rIn, const std::vector<std::int64_t>& rFaces) {
    double a = 0.0;
    for (std::int64_t f : rFaces)
        a += rIn.mArea[static_cast<std::size_t>(f)];
    return a;
}

// The primitive method's patches (steps 1 and 2 of segment.hpp's "Primitive
// method"): sets rPatch for the faces a primitive claims, returns the count.
std::int64_t seg_primitive_patches(SegPrimitiveInput& rIn, std::size_t NumFaces,
                                   std::vector<std::int64_t>& rPatch) {
    const std::size_t nf = NumFaces;
    // Smooth areas, numbered in the order of their lowest face.
    rIn.mSmoothArea.assign(nf, -1);
    std::vector<std::vector<std::int64_t>> areas;
    for (std::size_t seed = 0; seed < nf; ++seed) {
        if (rIn.mSmoothArea[seed] != -1)
            continue;
        const std::int64_t id = static_cast<std::int64_t>(areas.size());
        std::vector<std::int64_t> faces{static_cast<std::int64_t>(seed)};
        rIn.mSmoothArea[seed] = id;
        for (std::size_t q = 0; q < faces.size(); ++q) {
            const std::size_t g = static_cast<std::size_t>(faces[q]);
            for (std::int64_t k = (*rIn.mAdjStart)[g]; k < (*rIn.mAdjStart)[g + 1]; ++k) {
                const auto& a = (*rIn.mAdj)[static_cast<std::size_t>(k)];
                const std::size_t nb = static_cast<std::size_t>(a.first);
                if (rIn.mSmoothArea[nb] != -1 ||
                    !(*rIn.mPairs)[static_cast<std::size_t>(a.second)].mSmooth)
                    continue;
                rIn.mSmoothArea[nb] = id;
                faces.push_back(a.first);
            }
        }
        std::sort(faces.begin(), faces.end());
        areas.push_back(std::move(faces));
    }

    std::int64_t npatches = 0;
    std::vector<std::uint8_t> mark(nf, 0), tried(nf, 0);
    constexpr std::uint32_t kAll = (1u << 1) | (1u << 2) | (1u << 3) | (1u << 4) | (1u << 5);
    constexpr std::uint32_t kCurved = (1u << 2) | (1u << 3) | (1u << 4) | (1u << 5);
    for (const std::vector<std::int64_t>& area : areas) {
        // 1. One primitive for the whole smooth area.
        const detail::PrimitiveFit whole = seg_fit(rIn, area, kAll);
        if (whole.mType > 0) {
            for (std::int64_t f : area)
                rPatch[static_cast<std::size_t>(f)] = npatches;
            ++npatches;
            continue;
        }
        // 2. Grown regions, largest seed first; claimed for the area only once
        // it is clear the primitives explain it (below).
        const double area_total = seg_area_of(rIn, area);
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
            // A plane through the seed.
            std::vector<std::int64_t> plane_region;
            {
                detail::PrimitiveFit plane;
                plane.mType = 1;
                for (std::size_t d = 0; d < 3; ++d) {
                    plane.mParams[d] = (*rIn.mCentroid)[ss * 3 + d];
                    plane.mParams[3 + d] = (*rIn.mNormals)[ss * 3 + d];
                }
                if (seg_face_fits(rIn, s, plane, rIn.mGrowTol))
                    plane_region = seg_grow(rIn, {s}, plane, 1, rPatch, mark);
            }
            // A curved surface through the seed's neighbourhood.
            std::vector<std::int64_t> curved_region;
            std::vector<std::int64_t> hood{s};
            {
                mark[ss] = 1;
                std::size_t ring_end = 1;
                for (int ring = 0; ring < 2 && hood.size() < 32; ++ring) {
                    const std::size_t from = ring == 0 ? 0 : ring_end;
                    ring_end = hood.size();
                    for (std::size_t q = from; q < ring_end && hood.size() < 32; ++q) {
                        const std::size_t g = static_cast<std::size_t>(hood[q]);
                        for (std::int64_t k = (*rIn.mAdjStart)[g];
                             k < (*rIn.mAdjStart)[g + 1] && hood.size() < 32; ++k) {
                            const auto& a = (*rIn.mAdj)[static_cast<std::size_t>(k)];
                            const std::size_t nb = static_cast<std::size_t>(a.first);
                            if (mark[nb] || rPatch[nb] != -1 ||
                                !(*rIn.mPairs)[static_cast<std::size_t>(a.second)].mSmooth ||
                                rIn.mSmoothArea[nb] != rIn.mSmoothArea[ss])
                                continue;
                            mark[nb] = 1;
                            hood.push_back(a.first);
                        }
                    }
                }
                for (std::int64_t f : hood)
                    mark[static_cast<std::size_t>(f)] = 0;
                if (hood.size() >= 3) {
                    const detail::PrimitiveFit fit = seg_fit(rIn, hood, kCurved);
                    if (fit.mType >= 2)
                        curved_region = seg_grow(rIn, {s}, fit, hood.size(), rPatch, mark);
                }
            }
            // A region grown to the tolerance is kept only when its refit
            // leaves an RMS under 2% of it: a CAD tessellation lies on its
            // surfaces exactly, while a sphere or cone fitted to a patch of a
            // freeform or helical surface grows until it strays the tolerance.
            auto exact = [&](std::vector<std::int64_t>& rRegion, std::uint32_t Types) {
                if (rRegion.empty())
                    return;
                seg_trim(rIn, rRegion, Types, kSegExactShare * rIn.mTol, mark);
                if (rRegion.empty())
                    return;
                const detail::PrimitiveFit refit = seg_fit(rIn, rRegion, Types);
                if (refit.mType == 0 || !(refit.mRms <= kSegExactShare * rIn.mTol)) {
                    rRegion.clear();
                    return;
                }
                // A curved surface needs three points per parameter (4 for a
                // sphere up to 7 for a torus): a handful of points lies
                // exactly on some sphere or cylinder whatever the surface.
                if (refit.mType >= 2) {
                    std::vector<double> pts, nrm, ctr;
                    seg_samples(rIn, rRegion, pts, nrm, ctr);
                    if (pts.size() / 3 < static_cast<std::size_t>(3 * (refit.mType + 2)))
                        rRegion.clear();
                }
            };
            exact(plane_region, 1u << 1);
            exact(curved_region, kCurved);
            const double plane_area = seg_area_of(rIn, plane_region);
            const double curved_area = seg_area_of(rIn, curved_region);
            const std::vector<std::int64_t>* chosen = nullptr;
            if (curved_region.size() >= 3 && curved_area > 1.5 * plane_area)
                chosen = &curved_region;
            else if (plane_region.size() >= 2 && plane_area >= 0.02 * area_total)
                chosen = &plane_region;
            else if (curved_region.size() >= 3)
                chosen = &curved_region;
            if (chosen == nullptr) {
                // Neither surface holds here: no cell of this neighbourhood
                // seeds again (a freeform area would cost a fit per cell);
                // another seed's region may still grow over them.
                for (std::int64_t f : hood)
                    tried[static_cast<std::size_t>(f)] = 1;
                continue;
            }
            // Claim the faces now, so later seeds grow around them.
            for (std::int64_t f : *chosen)
                rPatch[static_cast<std::size_t>(f)] = npatches;
            ++npatches;
            grown.emplace_back(*chosen, chosen == &curved_region);
        }
        // A freeform surface fits small spheres and tori anywhere: a curved
        // region under 1% of its smooth area is dropped, and an area the kept
        // regions explain less than half of is left whole to the angle
        // growth, as the angle method segments it.
        double explained = 0.0;
        std::vector<std::uint8_t> keep(grown.size(), 1);
        for (std::size_t i = 0; i < grown.size(); ++i) {
            const double a = seg_area_of(rIn, grown[i].first);
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
    if (rOptions.mMethod == SegmentMethod::Primitive && nf > 0) {
        SegPrimitiveInput in;
        in.mStart = &start;
        in.mNodes = &nodes;
        in.mXyz = &xyz;
        in.mNormals = &normals;
        in.mCentroid = &centroid;
        in.mAdjStart = &adj_start;
        in.mAdj = &adj;
        in.mPairs = &pairs;
        in.mArea.assign(nf, 0.0);
        double lo[3], hi[3];
        for (int d = 0; d < 3; ++d) {
            lo[d] = std::numeric_limits<double>::infinity();
            hi[d] = -lo[d];
        }
        for (std::size_t f = 0; f < nf; ++f) {
            const std::int64_t* ring = nodes.data() + start[f];
            const std::size_t size = static_cast<std::size_t>(start[f + 1] - start[f]);
            double twice[3] = {0.0, 0.0, 0.0};
            for (std::size_t k = 0; k < size; ++k) {
                const double* a = &xyz[static_cast<std::size_t>(ring[k]) * 3];
                const double* b = &xyz[static_cast<std::size_t>(ring[(k + 1) % size]) * 3];
                twice[0] += a[1] * b[2] - a[2] * b[1];
                twice[1] += a[2] * b[0] - a[0] * b[2];
                twice[2] += a[0] * b[1] - a[1] * b[0];
                for (int d = 0; d < 3; ++d) {
                    lo[d] = std::min(lo[d], a[d]);
                    hi[d] = std::max(hi[d], a[d]);
                }
            }
            in.mArea[f] =
                0.5 * std::sqrt(twice[0] * twice[0] + twice[1] * twice[1] + twice[2] * twice[2]);
        }
        in.mDiag = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                             (hi[2] - lo[2]) * (hi[2] - lo[2]));
        in.mTol = rOptions.mTolerance * in.mDiag;
        in.mGrowTol = kSegExactShare * in.mTol;
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
