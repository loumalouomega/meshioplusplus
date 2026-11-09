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
#include "meshioplusplus/operations/deviation.hpp"
#include "meshioplusplus/operations/registration.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/surface_normals.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"
// Project includes (private, not installed)
#include "../detail/dual_graph.hpp"
#include "../detail/geometry_utils.hpp"
#include "../detail/primitive_geom.hpp"
#include "../detail/render_field.hpp"
#include "../detail/surface_edge_runs.hpp"
#include "../detail/surface_rings.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kDevPrefix = "meshio++: deviation: ";
constexpr double kDevPi = 3.14159265358979323846;

void dev_check_options(const DeviationOptions& rOptions) {
    if (!(rOptions.mTolerance >= 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kDevPrefix) + "tolerance must be finite and >= 0");
    if (rOptions.mBins < 1)
        throw std::invalid_argument(std::string(kDevPrefix) + "bins must be >= 1");
    if (!std::isfinite(rOptions.mRangeMin) || !std::isfinite(rOptions.mRangeMax))
        throw std::invalid_argument(std::string(kDevPrefix) + "the range must be finite");
    if (!(rOptions.mMaxNormalAngleDeg >= 0.0 && rOptions.mMaxNormalAngleDeg <= 180.0))
        throw std::invalid_argument(std::string(kDevPrefix) +
                                    "max_normal_angle must lie in [0, 180] degrees");
}

std::vector<detail::Vec3> dev_points(const Mesh& rScan) {
    const NDArray& pts = rScan.Points();
    const std::size_t n = rScan.NumPoints();
    const std::size_t dim = n > 0 ? pts.Size() / n : 3;
    std::vector<detail::Vec3> out(n);
    parallel_for(n, [&](std::size_t p) {
        out[p] = detail::read_point(pts, dim, static_cast<std::int64_t>(p));
    });
    return out;
}

// One unit normal per scan point: its `normals` point data, else the
// angle-weighted normal of the scan's own surface cells; zero where neither.
std::vector<detail::Vec3> dev_point_normals(const Mesh& rScan) {
    const std::size_t n = rScan.NumPoints();
    std::vector<detail::Vec3> out(n, detail::Vec3{0.0, 0.0, 0.0});
    if (rScan.HasPointData("normals")) {
        const NDArray& a = rScan.PointData("normals");
        if (a.Ndim() == 2 && detail::cols(a) == 3 && detail::rows(a) == n) {
            for (std::size_t p = 0; p < n; ++p) {
                const detail::Vec3 v{detail::read_double(a, p * 3),
                                     detail::read_double(a, p * 3 + 1),
                                     detail::read_double(a, p * 3 + 2)};
                const double len = detail::vec3_norm(v);
                if (len > 0.0 && std::isfinite(len))
                    out[p] = detail::vec3_scale(v, 1.0 / len);
            }
            return out;
        }
    }
    bool has_surface = false;
    for (std::size_t b = 0; b < rScan.NumCellBlocks(); ++b) {
        const std::string t = rScan.Cells(b).Type();
        has_surface = has_surface || t == "triangle" || t == "quad" || t.rfind("polygon", 0) == 0;
    }
    if (!has_surface)
        return out;
    detail::TriangleSoup soup;
    try {
        soup = detail::build_triangle_soup(rScan, "");
    } catch (const std::invalid_argument&) {
        return out;  // a volume or higher-order block: no surface normals
    }
    const std::vector<detail::Vec3> fn = detail::soup_face_normals(soup);
    const std::vector<detail::Vec3> vn =
        detail::accumulate_vertex_normals(soup, fn, SdfPseudonormalWeight::Angle);
    for (std::size_t p = 0; p < n && p < vn.size(); ++p) {
        const double len = detail::vec3_norm(vn[p]);
        if (len > 0.0 && std::isfinite(len))
            out[p] = detail::vec3_scale(vn[p], 1.0 / len);
    }
    return out;
}

// The 95th percentile rounded up to two significant figures.
double dev_round_up_2sig(double V) {
    if (!(V > 0.0) || !std::isfinite(V))
        return 0.0;
    const double mag = std::pow(10.0, std::floor(std::log10(V)) - 1.0);
    return std::ceil(V / mag - 1e-9) * mag;
}

// The statistics of the per-point distances (NaN = not measured), in point order.
void dev_summarize(const std::vector<double>& rDist, const DeviationOptions& rOptions,
                   DeviationResult& rOut) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    rOut.mNumPoints = static_cast<std::int64_t>(rDist.size());
    std::vector<double> d;
    d.reserve(rDist.size());
    for (double v : rDist)
        if (std::isfinite(v))
            d.push_back(v);
    const std::size_t m = d.size();
    rOut.mNumMeasured = static_cast<std::int64_t>(m);
    rOut.mHistogramEdges.assign(static_cast<std::size_t>(rOptions.mBins) + 1, 0.0);
    rOut.mHistogramCounts.assign(static_cast<std::size_t>(rOptions.mBins), 0);
    if (m == 0) {
        rOut.mMin = rOut.mMax = rOut.mMean = rOut.mRms = rOut.mSigma = nan;
        rOut.mPercentiles.fill(nan);
        rOut.mShareWithin = nan;
        rOut.mSuggestedRange = 0.0;
        rOut.mWarnings.emplace_back("no scan point was measured");
        return;
    }
    double lo = d[0], hi = d[0], sum = 0.0, sum2 = 0.0, amax = 0.0;
    std::size_t within = 0;
    for (double v : d) {
        lo = v < lo ? v : lo;
        hi = v > hi ? v : hi;
        sum += v;
        sum2 += v * v;
        amax = std::fabs(v) > amax ? std::fabs(v) : amax;
        within += std::fabs(v) <= rOptions.mTolerance ? 1 : 0;
    }
    const double mean = sum / static_cast<double>(m);
    double var = 0.0;
    for (double v : d)
        var += (v - mean) * (v - mean);
    rOut.mMin = lo;
    rOut.mMax = hi;
    rOut.mMean = mean;
    rOut.mRms = std::sqrt(sum2 / static_cast<double>(m));
    rOut.mSigma = std::sqrt(var / static_cast<double>(m));
    rOut.mShareWithin =
        rOptions.mTolerance > 0.0 ? static_cast<double>(within) / static_cast<double>(m) : nan;
    std::vector<double> sorted = d;
    std::sort(sorted.begin(), sorted.end());
    const double ps[5] = {1.0, 5.0, 50.0, 95.0, 99.0};
    for (std::size_t k = 0; k < 5; ++k)
        rOut.mPercentiles[k] = detail::render_percentile(sorted, ps[k]);
    std::vector<double> absd(m);
    for (std::size_t i = 0; i < m; ++i)
        absd[i] = std::fabs(d[i]);
    rOut.mSuggestedRange = dev_round_up_2sig(detail::render_percentile(absd, 95.0));

    double r0 = rOptions.mRangeMin, r1 = rOptions.mRangeMax;
    if (!(r0 < r1)) {
        r0 = -amax;
        r1 = amax;
        if (!(r1 > r0)) {  // every distance is zero
            r0 = -1.0;
            r1 = 1.0;
        }
    }
    const std::size_t bins = static_cast<std::size_t>(rOptions.mBins);
    const double width = (r1 - r0) / static_cast<double>(bins);
    for (std::size_t k = 0; k <= bins; ++k)
        rOut.mHistogramEdges[k] = k == bins ? r1 : r0 + width * static_cast<double>(k);
    for (double v : d) {
        if (v < r0 || v > r1) {
            ++rOut.mNumOutsideRange;
            continue;
        }
        std::size_t k = static_cast<std::size_t>((v - r0) / width);
        k = k >= bins ? bins - 1 : k;  // the last edge belongs to the last bin
        ++rOut.mHistogramCounts[k];
    }
}

void dev_attach(const Mesh& rScan, const std::vector<double>& rDist,
                const DeviationOptions& rOptions, DeviationResult& rOut) {
    rOut.mMesh = detail::clone_mesh(rScan);
    if (!rOptions.mAttach)
        return;
    NDArray a = NDArray::Uninit(DType::Float64, {rDist.size()});
    std::copy(rDist.begin(), rDist.end(), a.As<double>());
    rOut.mMesh.AddPointData(kDeviationDistanceName, std::move(a));
}

}  // namespace

DeviationResult deviation(const Mesh& rScan, const Mesh& rReference,
                          const DeviationOptions& rOptions) {
    dev_check_options(rOptions);
    if (rOptions.mRegister) {
        RegistrationOptions reg;
        reg.mMethod = RegistrationMethod::PointToPlane;
        reg.mInitialMode = RegistrationInitial::Auto;
        reg.mAdaptiveTrim = 3.0;
        reg.mTargetRegion = rOptions.mReferenceRegion;
        const RegistrationResult r = register_rigid(rScan, rReference, reg);
        const Mesh moved = transform(rScan, transform_from_matrix(r.mMatrix.data()), true);
        DeviationOptions plain = rOptions;
        plain.mRegister = false;
        DeviationResult out = deviation(moved, rReference, plain);
        out.mRegistered = true;
        out.mMatrix = r.mMatrix;
        out.mRegistrationRms = r.mRms;
        if (!r.mConverged)
            out.mWarnings.push_back("the registration did not converge in " +
                                    std::to_string(r.mIterations) + " iterations");
        for (const std::string& w : r.mWarnings)
            out.mWarnings.push_back("registration: " + w);
        return out;
    }
    const detail::TriangleSoup soup =
        detail::build_triangle_soup(rReference, rOptions.mReferenceRegion);
    DeviationResult result;
    if (soup.NumTriangles() == 0)
        throw std::invalid_argument(std::string(kDevPrefix) +
                                    "the reference has no surface cells to measure against");
    const detail::SurfaceEdgeRuns runs = detail::surface_edge_runs(soup);
    SurfaceDistanceOptions sd;
    const detail::DistanceQuery query = detail::build_distance_query_from_runs(soup, sd, runs);

    const std::vector<detail::Vec3> pts = dev_points(rScan);
    const std::size_t n = pts.size();
    const bool rule = rOptions.mMaxNormalAngleDeg < 180.0;
    std::vector<detail::Vec3> qn;
    std::vector<detail::Vec3> tn;
    if (rule) {
        qn = dev_point_normals(rScan);
        tn.resize(soup.NumTriangles());
        for (std::size_t t = 0; t < tn.size(); ++t) {
            const double len = detail::vec3_norm(query.mFaceNormal[t]);
            tn[t] = len > 0.0 ? detail::vec3_scale(query.mFaceNormal[t], 1.0 / len)
                              : detail::Vec3{0.0, 0.0, 0.0};
        }
        for (const detail::Vec3& v : qn)
            result.mNumWithoutNormal += (v[0] == 0.0 && v[1] == 0.0 && v[2] == 0.0) ? 1 : 0;
        if (result.mNumWithoutNormal > 0)
            result.mWarnings.push_back(std::to_string(result.mNumWithoutNormal) +
                                       " scan point(s) have no normal and were measured "
                                       "without the far-side rule");
    }
    const double cos_max = std::cos(rOptions.mMaxNormalAngleDeg * kDevPi / 180.0);
    const detail::TriangleFilter accept = [&](std::size_t Q, std::int64_t T) {
        if (!rule)
            return true;
        const detail::Vec3& a = qn[Q];
        if (a[0] == 0.0 && a[1] == 0.0 && a[2] == 0.0)
            return true;
        const detail::Vec3& b = tn[static_cast<std::size_t>(T)];
        if (b[0] == 0.0 && b[1] == 0.0 && b[2] == 0.0)
            return true;
        return detail::vec3_dot(a, b) >= cos_max;
    };
    const std::vector<detail::SurfaceProjection> hits = detail::query_surface_projections_filtered(
        query, pts, accept, std::numeric_limits<double>::infinity());
    std::vector<double> dist(n, std::numeric_limits<double>::quiet_NaN());
    parallel_for(n, [&](std::size_t p) {
        const detail::SurfaceProjection& h = hits[p];
        if (!h.mFound)
            return;
        const double side = detail::vec3_dot(detail::vec3_sub(pts[p], h.mPoint), h.mNormal);
        dist[p] = side < 0.0 ? -h.mDistance : h.mDistance;
    });
    dev_summarize(dist, rOptions, result);
    dev_attach(rScan, dist, rOptions, result);
    return result;
}

DeviationResult deviation_to_primitive(const Mesh& rScan, std::int64_t Type,
                                       const std::array<double, 10>& rParams,
                                       const DeviationOptions& rOptions) {
    dev_check_options(rOptions);
    if (rOptions.mRegister)
        throw std::invalid_argument(std::string(kDevPrefix) +
                                    "register needs a reference surface, not a primitive");
    if (Type < 1 || Type > 5)
        throw std::invalid_argument(std::string(kDevPrefix) + "unknown primitive type " +
                                    std::to_string(Type) +
                                    " (1 plane, 2 sphere, 3 cylinder, 4 cone, 5 torus)");
    for (double v : rParams)
        if (!std::isfinite(v))
            throw std::invalid_argument(std::string(kDevPrefix) +
                                        "the primitive's parameters must be finite");
    const std::vector<detail::Vec3> pts = dev_points(rScan);
    std::vector<double> dist(pts.size());
    parallel_for(pts.size(), [&](std::size_t p) {
        dist[p] = detail::primitive_signed_distance(Type, rParams.data(), pts[p].data());
    });
    DeviationResult result;
    dev_summarize(dist, rOptions, result);
    dev_attach(rScan, dist, rOptions, result);
    return result;
}

DeviationRegionsResult deviation_regions(const Mesh& rMesh,
                                         const DeviationRegionsOptions& rOptions) {
    if (!(rOptions.mTolerance > 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kDevPrefix) +
                                    "deviation_regions: tolerance must be finite and > 0");
    if (!(rOptions.mMinArea >= 0.0) || !std::isfinite(rOptions.mMinArea))
        throw std::invalid_argument(std::string(kDevPrefix) +
                                    "deviation_regions: min_area must be finite and >= 0");
    if (!rMesh.HasPointData(rOptions.mField))
        throw std::invalid_argument(std::string(kDevPrefix) + "deviation_regions: no point data '" +
                                    rOptions.mField + "' (run deviation first)");
    const NDArray& field = rMesh.PointData(rOptions.mField);
    const std::size_t np = rMesh.NumPoints();
    if (field.Size() != np)
        throw std::invalid_argument(std::string(kDevPrefix) + "deviation_regions: '" +
                                    rOptions.mField + "' must hold one value per point");
    std::vector<double> value(np);
    for (std::size_t p = 0; p < np; ++p)
        value[p] = detail::read_double(field, p);

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t ncells = static_cast<std::size_t>(detail::total_cells(bases));
    std::vector<std::uint8_t> all(ncells, 1);
    const detail::SurfaceRings rings =
        detail::surface_rings(rMesh, all, "meshio++: deviation: deviation_regions: ");
    const std::vector<double> xyz = detail::surface_xyz(rMesh);

    // Per cell: its class (+1, -1, 0), value, area and area-weighted centroid.
    std::vector<std::int8_t> cls(ncells, 0);
    std::vector<double> cval(ncells, 0.0), carea(ncells, 0.0);
    std::vector<detail::Vec3> ccen(ncells, detail::Vec3{0.0, 0.0, 0.0});
    std::vector<std::int64_t> ring_of(ncells, -1);
    for (std::size_t f = 0; f < rings.NumFaces(); ++f)
        ring_of[static_cast<std::size_t>(rings.mCell[f])] = static_cast<std::int64_t>(f);
    parallel_for(rings.NumFaces(), [&](std::size_t f) {
        const std::size_t c = static_cast<std::size_t>(rings.mCell[f]);
        const std::int64_t* ring = rings.mNodes.data() + rings.mStart[f];
        const std::size_t size = static_cast<std::size_t>(rings.mStart[f + 1] - rings.mStart[f]);
        double sum = 0.0;
        std::size_t finite = 0;
        for (std::size_t k = 0; k < size; ++k) {
            const double v = value[static_cast<std::size_t>(ring[k])];
            if (std::isfinite(v)) {
                sum += v;
                ++finite;
            }
        }
        const double v = finite ? sum / static_cast<double>(finite) : 0.0;
        cval[c] = v;
        cls[c] = !finite ? 0 : (v > rOptions.mTolerance ? 1 : (v < -rOptions.mTolerance ? -1 : 0));
        double area = 0.0;
        detail::Vec3 cen{0.0, 0.0, 0.0};
        const double* a = &xyz[static_cast<std::size_t>(ring[0]) * 3];
        for (std::size_t k = 1; k + 1 < size; ++k) {
            const double* b = &xyz[static_cast<std::size_t>(ring[k]) * 3];
            const double* d = &xyz[static_cast<std::size_t>(ring[k + 1]) * 3];
            const double t = detail::GeometryUtils::TriangleArea(a, b, d);
            area += t;
            for (std::size_t i = 0; i < 3; ++i)
                cen[i] += t * (a[i] + b[i] + d[i]) / 3.0;
        }
        carea[c] = area;
        ccen[c] = area > 0.0 ? detail::vec3_scale(cen, 1.0 / area) : cen;
    });

    // Components of each class over the edge-sharing dual graph, in ascending
    // cell order.
    const detail::DualGraphCsr graph = detail::build_facet_dual_graph(rMesh, ncells, true);
    std::vector<std::int64_t> comp(ncells, -1);
    std::vector<std::vector<std::int64_t>> members;
    for (std::size_t seed = 0; seed < ncells; ++seed) {
        if (cls[seed] == 0 || comp[seed] >= 0)
            continue;
        const std::int64_t id = static_cast<std::int64_t>(members.size());
        members.emplace_back();
        std::vector<std::int64_t> stack{static_cast<std::int64_t>(seed)};
        comp[seed] = id;
        while (!stack.empty()) {
            const std::int64_t c = stack.back();
            stack.pop_back();
            members.back().push_back(c);
            for (std::int64_t e = graph.mXadj[static_cast<std::size_t>(c)];
                 e < graph.mXadj[static_cast<std::size_t>(c) + 1]; ++e) {
                const std::int64_t n = graph.mAdjncy[static_cast<std::size_t>(e)];
                if (comp[static_cast<std::size_t>(n)] < 0 &&
                    cls[static_cast<std::size_t>(n)] == cls[seed]) {
                    comp[static_cast<std::size_t>(n)] = id;
                    stack.push_back(n);
                }
            }
        }
        std::sort(members.back().begin(), members.back().end());
    }

    DeviationRegionsResult result;
    struct DevComp {
        DeviationRegionRecord mRec;
        std::int64_t mFirst;
        std::size_t mIndex;
    };
    std::vector<DevComp> comps;
    for (std::size_t k = 0; k < members.size(); ++k) {
        DeviationRegionRecord r;
        r.mSign = cls[static_cast<std::size_t>(members[k][0])];
        r.mNumCells = static_cast<std::int64_t>(members[k].size());
        double lo[3] = {std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
        double hi[3] = {-lo[0], -lo[1], -lo[2]};
        r.mExtreme = cval[static_cast<std::size_t>(members[k][0])];
        detail::Vec3 cen{0.0, 0.0, 0.0};
        for (std::int64_t c : members[k]) {
            const std::size_t ci = static_cast<std::size_t>(c);
            r.mArea += carea[ci];
            r.mVolume += cval[ci] * carea[ci];
            cen = detail::vec3_add(cen, detail::vec3_scale(ccen[ci], carea[ci]));
            r.mExtreme =
                r.mSign > 0 ? std::max(r.mExtreme, cval[ci]) : std::min(r.mExtreme, cval[ci]);
            const std::size_t f = static_cast<std::size_t>(ring_of[ci]);
            for (std::int64_t q = rings.mStart[f]; q < rings.mStart[f + 1]; ++q)
                for (std::size_t i = 0; i < 3; ++i) {
                    const double x =
                        xyz[static_cast<std::size_t>(rings.mNodes[static_cast<std::size_t>(q)]) *
                                3 +
                            i];
                    lo[i] = std::min(lo[i], x);
                    hi[i] = std::max(hi[i], x);
                }
        }
        r.mMean = r.mArea > 0.0 ? r.mVolume / r.mArea : 0.0;
        if (r.mArea > 0.0)
            cen = detail::vec3_scale(cen, 1.0 / r.mArea);
        r.mCentroid = cen;
        r.mBounds = {lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]};
        if (r.mArea < rOptions.mMinArea) {
            ++result.mNumDropped;
            continue;
        }
        comps.push_back({std::move(r), members[k][0], k});
    }
    std::sort(comps.begin(), comps.end(), [](const DevComp& rA, const DevComp& rB) {
        const double a = std::fabs(rA.mRec.mVolume), b = std::fabs(rB.mRec.mVolume);
        return a != b ? a > b : rA.mFirst < rB.mFirst;
    });

    Mesh& out = result.mMesh;
    out = detail::clone_mesh(rMesh);
    std::vector<std::int64_t> region_of(ncells, -1);
    for (std::size_t i = 0; i < comps.size(); ++i) {
        DeviationRegionRecord& r = comps[i].mRec;
        const bool excess = r.mSign > 0;
        r.mName = rOptions.mPrefix + (excess ? "excess_" : "missing_") +
                  std::to_string(excess ? result.mNumExcess++ : result.mNumMissing++);
        const std::vector<std::int64_t>& cells = members[comps[i].mIndex];
        NDArray entries = NDArray::Uninit(DType::Int64, {cells.size()});
        std::copy(cells.begin(), cells.end(), entries.As<std::int64_t>());
        for (std::int64_t c : cells)
            region_of[static_cast<std::size_t>(c)] = static_cast<std::int64_t>(i);
        for (std::size_t k = out.NumRegions(); k-- > 0;)
            if (out.Region(k).mKind == RegionKind::Cell && out.Region(k).mName == r.mName)
                out.RemoveRegion(k);
        meshioplusplus::Region region(r.mName, RegionKind::Cell, std::move(entries));
        region.Canonicalize();
        out.AddRegion(std::move(region));
        result.mRegions.push_back(std::move(r));
    }
    std::vector<NDArray> blocks;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const std::size_t base = static_cast<std::size_t>(bases[b]);
        const std::size_t nc = rMesh.Cells(b).NumCells();
        NDArray a = NDArray::Uninit(DType::Int64, {nc});
        std::copy(region_of.begin() + static_cast<std::ptrdiff_t>(base),
                  region_of.begin() + static_cast<std::ptrdiff_t>(base + nc), a.As<std::int64_t>());
        blocks.push_back(std::move(a));
    }
    out.AddCellData(kDeviationRegionName, std::move(blocks));
    if (comps.empty())
        result.mWarnings.emplace_back("no cell lies beyond the tolerance");
    return result;
}

ExaggerateResult exaggerate(const Mesh& rMesh, const ExaggerateOptions& rOptions) {
    if (!std::isfinite(rOptions.mScale))
        throw std::invalid_argument(std::string(kDevPrefix) + "exaggerate: scale must be finite");
    if (!(rOptions.mClamp >= 0.0) || !std::isfinite(rOptions.mClamp))
        throw std::invalid_argument(std::string(kDevPrefix) +
                                    "exaggerate: clamp must be finite and >= 0");
    if (!(rOptions.mSmoothFraction >= 0.0) || !std::isfinite(rOptions.mSmoothFraction))
        throw std::invalid_argument(std::string(kDevPrefix) +
                                    "exaggerate: smooth_fraction must be finite and >= 0");
    if (!rMesh.HasPointData(rOptions.mField))
        throw std::invalid_argument(std::string(kDevPrefix) + "exaggerate: no point data '" +
                                    rOptions.mField + "' (run deviation first)");
    const NDArray& field = rMesh.PointData(rOptions.mField);
    const std::size_t n = rMesh.NumPoints();
    if (field.Size() != n)
        throw std::invalid_argument(std::string(kDevPrefix) + "exaggerate: '" + rOptions.mField +
                                    "' must hold one value per point");
    ExaggerateResult result;
    std::vector<double> d(n);
    std::vector<double> absd;
    for (std::size_t p = 0; p < n; ++p) {
        d[p] = detail::read_double(field, p);
        if (std::isfinite(d[p]))
            absd.push_back(std::fabs(d[p]));
    }
    const double clamp =
        rOptions.mClamp > 0.0
            ? rOptions.mClamp
            : (absd.empty() ? 0.0 : dev_round_up_2sig(detail::render_percentile(absd, 95.0)));
    result.mClamp = clamp;
    for (double& v : d)
        if (std::isfinite(v))
            v = std::max(-clamp, std::min(clamp, v));

    const std::vector<detail::Vec3> pts = dev_points(rMesh);
    // The smoothing radius: a fraction of the bounding-box diagonal.
    detail::Vec3 lo{0.0, 0.0, 0.0}, hi{0.0, 0.0, 0.0};
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = p == 0 ? pts[p][k] : std::min(lo[k], pts[p][k]);
            hi[k] = p == 0 ? pts[p][k] : std::max(hi[k], pts[p][k]);
        }
    const double radius = rOptions.mSmoothFraction * detail::vec3_norm(detail::vec3_sub(hi, lo));
    result.mRadius = radius;
    std::vector<double> smooth = d;
    if (radius > 0.0 && n > 0) {
        detail::SpatialGrid grid(radius);
        for (std::size_t p = 0; p < n; ++p)
            if (std::isfinite(d[p]))
                grid.Insert(grid.KeyOf(pts[p].data()), static_cast<std::int64_t>(p));
        const double r2 = radius * radius;
        parallel_for(n, [&](std::size_t p) {
            if (!std::isfinite(d[p]))
                return;
            double sw = 0.0, sv = 0.0;
            grid.ForEachIn27(grid.KeyOf(pts[p].data()), [&](const std::vector<std::int64_t>& rIds) {
                for (std::int64_t q : rIds) {
                    const double dist2 = detail::vec3_norm_sq(
                        detail::vec3_sub(pts[static_cast<std::size_t>(q)], pts[p]));
                    if (dist2 >= r2)
                        continue;
                    const double w = 1.0 - dist2 / r2;
                    sw += w;
                    sv += w * d[static_cast<std::size_t>(q)];
                }
                return true;
            });
            smooth[p] = sw > 0.0 ? sv / sw : d[p];
        });
    }

    const std::vector<detail::Vec3> nrm = dev_point_normals(rMesh);
    std::vector<detail::Vec3> moved = pts;
    for (std::size_t p = 0; p < n; ++p) {
        const detail::Vec3& u = nrm[p];
        if (!std::isfinite(smooth[p]) || (u[0] == 0.0 && u[1] == 0.0 && u[2] == 0.0))
            continue;
        moved[p] = detail::vec3_add(pts[p], detail::vec3_scale(u, rOptions.mScale * smooth[p]));
        ++result.mNumMoved;
    }
    if (result.mNumMoved < static_cast<std::int64_t>(n))
        result.mWarnings.push_back(std::to_string(static_cast<std::int64_t>(n) - result.mNumMoved) +
                                   " point(s) without a normal or a deviation were not moved");
    result.mMesh = detail::clone_mesh(rMesh);
    const std::size_t dim = n > 0 ? rMesh.Points().Size() / n : 3;
    NDArray out = NDArray::Uninit(DType::Float64, {n, dim});
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t k = 0; k < dim; ++k)
            out.As<double>()[p * dim + k] =
                k < 3 ? moved[p][k] : detail::read_double(rMesh.Points(), p * dim + k);
    result.mMesh.AssignPoints(std::move(out));
    return result;
}

}  // namespace meshioplusplus
