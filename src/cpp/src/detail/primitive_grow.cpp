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
// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/parallel.hpp"
// Project includes (private, not installed)
#include "crease_edges.hpp"
#include "primitive_grow.hpp"

namespace meshioplusplus {
namespace detail {

namespace {
constexpr double kPgrPi = 3.14159265358979323846;
}  // namespace

PrimitiveGrowInput PrimitiveGrowSurface::Input() const {
    PrimitiveGrowInput in;
    in.mStart = &mRings.mStart;
    in.mNodes = &mRings.mNodes;
    in.mXyz = &mXyz;
    in.mNormals = &mNormals;
    in.mCentroid = &mCentroid;
    in.mAdjStart = &mAdjStart;
    in.mAdj = &mAdj;
    in.mPairs = &mPairs;
    return in;
}

PrimitiveGrowSurface build_primitive_grow_surface(const Mesh& rMesh,
                                                  const std::vector<std::uint8_t>& rSelected,
                                                  double FeatureAngleDeg, const char* pPrefix) {
    PrimitiveGrowSurface out;
    out.mRings = surface_rings(rMesh, rSelected, pPrefix);
    const std::vector<std::int64_t>& start = out.mRings.mStart;
    const std::vector<std::int64_t>& nodes = out.mRings.mNodes;
    const std::size_t nf = out.mRings.NumFaces();

    // Flat coordinates padded to 3-D, unit normals and centroids per face.
    out.mXyz = surface_xyz(rMesh);
    const std::vector<double>& xyz = out.mXyz;
    std::vector<double>& normals = out.mNormals;
    std::vector<double>& centroid = out.mCentroid;
    normals.resize(nf * 3);
    centroid.resize(nf * 3);
    parallel_for(nf, [&](std::size_t f) {
        const std::int64_t* ring = nodes.data() + start[f];
        const std::size_t size = static_cast<std::size_t>(start[f + 1] - start[f]);
        ring_unit_normal(xyz.data(), ring, size, &normals[f * 3]);
        double c[3] = {0.0, 0.0, 0.0};
        for (std::size_t k = 0; k < size; ++k)
            for (std::size_t d = 0; d < 3; ++d)
                c[d] += xyz[static_cast<std::size_t>(ring[k]) * 3 + d];
        for (std::size_t d = 0; d < 3; ++d)
            centroid[f * 3 + d] = c[d] / static_cast<double>(size);
    });

    // Every non-collapsed ring edge is one use; sorted, each run is one edge.
    const std::vector<CreaseEdgeUse> uses = crease_edge_uses(start, nodes);

    const double cos_feature = std::cos(FeatureAngleDeg * kPgrPi / 180.0);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<PrimitiveGrowPair>& pairs = out.mPairs;
    for (std::size_t b = 0; b < uses.size();) {
        std::size_t e = b + 1;
        while (e < uses.size() && uses[e].mLo == uses[b].mLo && uses[e].mHi == uses[b].mHi)
            ++e;
        const std::size_t u = e - b;
        if (u == 1) {
            ++out.mNumBoundaryEdges;
        } else if (u >= 3) {
            ++out.mNumNonManifoldEdges;
        } else if (uses[b].mFace != uses[b + 1].mFace) {
            PrimitiveGrowPair p;
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
                ++out.mNumFeatureEdges;
            pairs.push_back(p);
        }
        b = e;
    }

    // Per face, its neighbours across paired edges, ascending by face then edge.
    std::vector<std::int64_t>& adj_start = out.mAdjStart;
    adj_start.assign(nf + 1, 0);
    for (const PrimitiveGrowPair& p : pairs) {
        ++adj_start[static_cast<std::size_t>(p.mA) + 1];
        ++adj_start[static_cast<std::size_t>(p.mB) + 1];
    }
    for (std::size_t f = 0; f < nf; ++f)
        adj_start[f + 1] += adj_start[f];
    std::vector<std::pair<std::int64_t, std::int64_t>>& adj = out.mAdj;
    adj.resize(static_cast<std::size_t>(adj_start[nf]));
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
    return out;
}

void prepare_primitive_grow(PrimitiveGrowInput& rIn, std::size_t NumFaces, double RelTolerance) {
    const std::size_t nf = NumFaces;
    rIn.mArea.assign(nf, 0.0);
    double lo[3], hi[3];
    for (int d = 0; d < 3; ++d) {
        lo[d] = std::numeric_limits<double>::infinity();
        hi[d] = -lo[d];
    }
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t* ring = rIn.mNodes->data() + (*rIn.mStart)[f];
        const std::size_t size = static_cast<std::size_t>((*rIn.mStart)[f + 1] - (*rIn.mStart)[f]);
        double twice[3] = {0.0, 0.0, 0.0};
        for (std::size_t k = 0; k < size; ++k) {
            const double* a = &(*rIn.mXyz)[static_cast<std::size_t>(ring[k]) * 3];
            const double* b = &(*rIn.mXyz)[static_cast<std::size_t>(ring[(k + 1) % size]) * 3];
            twice[0] += a[1] * b[2] - a[2] * b[1];
            twice[1] += a[2] * b[0] - a[0] * b[2];
            twice[2] += a[0] * b[1] - a[1] * b[0];
            for (int d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], a[d]);
                hi[d] = std::max(hi[d], a[d]);
            }
        }
        rIn.mArea[f] =
            0.5 * std::sqrt(twice[0] * twice[0] + twice[1] * twice[1] + twice[2] * twice[2]);
    }
    rIn.mDiag = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                          (hi[2] - lo[2]) * (hi[2] - lo[2]));
    rIn.mTol = RelTolerance * rIn.mDiag;
    rIn.mGrowTol = kPrimitiveGrowExactShare * rIn.mTol;
}

std::vector<std::vector<std::int64_t>> primitive_smooth_areas(PrimitiveGrowInput& rIn,
                                                              std::size_t NumFaces) {
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
    return areas;
}

void primitive_grow_samples(const PrimitiveGrowInput& rIn, const std::vector<std::int64_t>& rFaces,
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
PrimitiveFit primitive_grow_fit(const PrimitiveGrowInput& rIn,
                                const std::vector<std::int64_t>& rFaces, std::uint32_t Types) {
    std::vector<double> pts, nrm, ctr;
    primitive_grow_samples(rIn, rFaces, pts, nrm, ctr);
    PrimitiveFit fit = fit_primitive_samples(pts, nrm, ctr, rIn.mTol, Types, 1.0, 50);
    if (fit.mType >= 2) {
        bool sane = true;
        for (int k = 0; k < kPrimitiveParams; ++k)
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

bool primitive_grow_face_fits(const PrimitiveGrowInput& rIn, std::int64_t F,
                              const PrimitiveFit& rFit, double Within) {
    for (std::int64_t k = (*rIn.mStart)[static_cast<std::size_t>(F)];
         k < (*rIn.mStart)[static_cast<std::size_t>(F) + 1]; ++k) {
        const double* p =
            &(*rIn.mXyz)[static_cast<std::size_t>((*rIn.mNodes)[static_cast<std::size_t>(k)]) * 3];
        if (!(primitive_point_distance(rFit.mType, rFit.mParams, p) <= Within))
            return false;
    }
    return true;
}

// Grow a region from rFaces across smooth edges, within the faces' smooth
// area, over faces not yet claimed, while every corner lies within the
// growth tolerance of the primitive. The primitive was fitted to FittedOn faces; it
// is refitted (same type) each time the region doubles past that, never on
// fewer faces, where a curved fit is ill-conditioned.
std::vector<std::int64_t> primitive_grow(const PrimitiveGrowInput& rIn,
                                         std::vector<std::int64_t> Faces, PrimitiveFit Fit,
                                         std::size_t FittedOn,
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
                rIn.mSmoothArea[nb] != area ||
                !primitive_grow_face_fits(rIn, a.first, Fit, rIn.mGrowTol))
                continue;
            rMark[nb] = 1;
            Faces.push_back(a.first);
        }
        if (Faces.size() >= 2 * fitted_at) {
            fitted_at = Faces.size();
            const PrimitiveFit refit =
                primitive_grow_fit(rIn, Faces, 1u << static_cast<unsigned>(Fit.mType));
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
void primitive_grow_trim(const PrimitiveGrowInput& rIn, std::vector<std::int64_t>& rRegion,
                         std::uint32_t Types, double Within, std::vector<std::uint8_t>& rMark) {
    for (int pass = 0; pass < 3 && !rRegion.empty(); ++pass) {
        const PrimitiveFit fit = primitive_grow_fit(rIn, rRegion, Types);
        if (fit.mType == 0) {
            rRegion.clear();
            return;
        }
        std::size_t kept = 0;
        for (std::int64_t f : rRegion) {
            bool on = true;
            for (std::int64_t k = (*rIn.mStart)[static_cast<std::size_t>(f)];
                 k < (*rIn.mStart)[static_cast<std::size_t>(f) + 1] && on; ++k)
                on = primitive_point_distance(
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

double primitive_grow_area_of(const PrimitiveGrowInput& rIn,
                              const std::vector<std::int64_t>& rFaces) {
    double a = 0.0;
    for (std::int64_t f : rFaces)
        a += rIn.mArea[static_cast<std::size_t>(f)];
    return a;
}

PrimitiveGrowRegion primitive_grow_from_seed(const PrimitiveGrowInput& rIn, std::int64_t Seed,
                                             bool TryPlane, std::uint32_t CurvedTypes,
                                             double AreaTotal, double MinPlaneShare,
                                             double ExactShare,
                                             const std::vector<std::int64_t>& rPatch,
                                             std::vector<std::uint8_t>& rMark) {
    const std::int64_t s = Seed;
    const std::size_t ss = static_cast<std::size_t>(s);
    std::vector<std::uint8_t>& mark = rMark;
    // A plane through the seed.
    std::vector<std::int64_t> plane_region;
    {
        PrimitiveFit plane;
        plane.mType = 1;
        for (std::size_t d = 0; d < 3; ++d) {
            plane.mParams[d] = (*rIn.mCentroid)[ss * 3 + d];
            plane.mParams[3 + d] = (*rIn.mNormals)[ss * 3 + d];
        }
        if (TryPlane && primitive_grow_face_fits(rIn, s, plane, rIn.mGrowTol))
            plane_region = primitive_grow(rIn, {s}, plane, 1, rPatch, mark);
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
            const PrimitiveFit fit = primitive_grow_fit(rIn, hood, CurvedTypes);
            if (fit.mType >= 2)
                curved_region = primitive_grow(rIn, {s}, fit, hood.size(), rPatch, mark);
        }
    }
    // A region grown to the tolerance is kept only when its refit
    // leaves an RMS under 2% of it: a CAD tessellation lies on its
    // surfaces exactly, while a sphere or cone fitted to a patch of a
    // freeform or helical surface grows until it strays the tolerance.
    auto exact = [&](std::vector<std::int64_t>& rRegion, std::uint32_t Types) {
        if (rRegion.empty())
            return;
        primitive_grow_trim(rIn, rRegion, Types, ExactShare * rIn.mTol, mark);
        if (rRegion.empty())
            return;
        const PrimitiveFit refit = primitive_grow_fit(rIn, rRegion, Types);
        if (refit.mType == 0 || !(refit.mRms <= ExactShare * rIn.mTol)) {
            rRegion.clear();
            return;
        }
        // A curved surface needs three points per parameter (4 for a
        // sphere up to 7 for a torus): a handful of points lies
        // exactly on some sphere or cylinder whatever the surface.
        if (refit.mType >= 2) {
            std::vector<double> pts, nrm, ctr;
            primitive_grow_samples(rIn, rRegion, pts, nrm, ctr);
            if (pts.size() / 3 < static_cast<std::size_t>(3 * (refit.mType + 2)))
                rRegion.clear();
        }
    };
    exact(plane_region, 1u << 1);
    exact(curved_region, CurvedTypes);
    const double plane_area = primitive_grow_area_of(rIn, plane_region);
    const double curved_area = primitive_grow_area_of(rIn, curved_region);
    const std::vector<std::int64_t>* chosen = nullptr;
    if (curved_region.size() >= 3 && curved_area > 1.5 * plane_area)
        chosen = &curved_region;
    else if (plane_region.size() >= 2 && plane_area >= MinPlaneShare * AreaTotal)
        chosen = &plane_region;
    else if (curved_region.size() >= 3)
        chosen = &curved_region;
    PrimitiveGrowRegion out;
    out.mHood = std::move(hood);
    if (chosen != nullptr) {
        out.mCurved = chosen == &curved_region;
        out.mFaces = *chosen;
    }
    return out;
}

std::int64_t nearest_centroid_face(const std::vector<double>& rCentroid, const double* pPoint) {
    std::int64_t best = -1;
    double best_d2 = std::numeric_limits<double>::infinity();
    const std::size_t nf = rCentroid.size() / 3;
    for (std::size_t f = 0; f < nf; ++f) {
        const double dx = rCentroid[f * 3] - pPoint[0];
        const double dy = rCentroid[f * 3 + 1] - pPoint[1];
        const double dz = rCentroid[f * 3 + 2] - pPoint[2];
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < best_d2) {
            best_d2 = d2;
            best = static_cast<std::int64_t>(f);
        }
    }
    return best;
}

}  // namespace detail
}  // namespace meshioplusplus
