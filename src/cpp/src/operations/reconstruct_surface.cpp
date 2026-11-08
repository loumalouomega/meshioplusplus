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
// A triangle surface from an oriented point cloud: Hoppe et al.'s signed
// distance to the samples' tangent planes, evaluated on a band of lattice
// cells about the samples and contoured by the marching-tetrahedra cutter.
// See operations/reconstruct_surface.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/reconstruct_surface.hpp"
#include "../detail/point_knn.hpp"
#include "../detail/bounding_box_utils.hpp"
#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/detail/marching.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/operations/point_normals.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/detail/geometry.hpp"

namespace meshioplusplus {

namespace {

using RcsVec = std::array<double, 3>;

constexpr const char* kRcsPrefix = "meshio++: reconstruct_surface: ";

[[noreturn]] void rcs_fail(const std::string& rWhat) {
    throw std::invalid_argument(kRcsPrefix + rWhat);
}

/// Sorts and deduplicates @p rIds.
void rcs_unique(std::vector<std::int64_t>& rIds) {
    std::sort(rIds.begin(), rIds.end());
    rIds.erase(std::unique(rIds.begin(), rIds.end()), rIds.end());
}

/// The lattice index set @p rIds (linear, x fastest, dims @p rDims) grown by
/// @p Lo cells down and @p Hi up along @p Axis, clamped to the lattice.
void rcs_dilate(std::vector<std::int64_t>& rIds, const std::array<std::int64_t, 3>& rDims, int Axis,
                std::int64_t Lo, std::int64_t Hi) {
    const std::int64_t stride = Axis == 0 ? 1 : (Axis == 1 ? rDims[0] : rDims[0] * rDims[1]);
    std::vector<std::int64_t> out;
    out.reserve(rIds.size() * static_cast<std::size_t>(Lo + Hi + 1));
    for (const std::int64_t id : rIds) {
        const std::int64_t c = (id / stride) % rDims[static_cast<std::size_t>(Axis)];
        for (std::int64_t d = -Lo; d <= Hi; ++d)
            if (c + d >= 0 && c + d < rDims[static_cast<std::size_t>(Axis)])
                out.push_back(id + d * stride);
    }
    rcs_unique(out);
    rIds = std::move(out);
}

/// Union-find root with path halving.
std::int64_t rcs_find(std::vector<std::int64_t>& rParent, std::int64_t I) {
    while (rParent[static_cast<std::size_t>(I)] != I) {
        rParent[static_cast<std::size_t>(I)] =
            rParent[static_cast<std::size_t>(rParent[static_cast<std::size_t>(I)])];
        I = rParent[static_cast<std::size_t>(I)];
    }
    return I;
}

}  // namespace

ReconstructSurfaceResult reconstruct_surface(const Mesh& rMesh,
                                             const ReconstructSurfaceOptions& rOptions) {
    if (rOptions.mK < 1)
        rcs_fail("k must be at least 1, not " + std::to_string(rOptions.mK));
    if (!std::isfinite(rOptions.mCellSize) || rOptions.mCellSize < 0.0)
        rcs_fail("cell_size must be a non-negative number");
    if (rOptions.mResolution < 0)
        rcs_fail("resolution must be non-negative");
    if (rOptions.mCellSize > 0.0 && rOptions.mResolution > 0)
        rcs_fail("give cell_size or resolution, not both");
    if (!std::isfinite(rOptions.mMaxDistance) || rOptions.mMaxDistance < 0.0)
        rcs_fail("max_distance must be a non-negative number");
    if (rOptions.mMaxCells < 1)
        rcs_fail("max_cells must be positive");

    const std::size_t n = rMesh.NumPoints();
    if (n > 0 && rMesh.PointDim() != 3)
        rcs_fail("the points must have three coordinates");

    ReconstructSurfaceResult result;
    // The normals: the named array, or estimated.
    PointNormalsResult estimated;
    const NDArray* p_normals = nullptr;
    NDArray assembled;
    const bool ply_columns = rOptions.mNormalsName == kNormalsName &&
                             !rMesh.HasPointData(kNormalsName) && rMesh.HasPointData("nx") &&
                             rMesh.HasPointData("ny") && rMesh.HasPointData("nz");
    if (rMesh.HasPointData(rOptions.mNormalsName)) {
        p_normals = &rMesh.PointData(rOptions.mNormalsName);
    } else if (ply_columns) {
        // A PLY scan's normals are three vertex properties, read as scalars.
        assembled = NDArray(DType::Float64, {n, std::size_t{3}});
        double* pa = assembled.As<double>();
        const char* names[3] = {"nx", "ny", "nz"};
        for (std::size_t a = 0; a < 3; ++a) {
            const NDArray& col = rMesh.PointData(names[a]);
            if (!detail::is_float_dtype(col.Dtype()) || col.Size() != n)
                rcs_fail(std::string("'") + names[a] + "' must be a float (n,) array");
            for (std::size_t i = 0; i < n; ++i)
                pa[i * 3 + a] = detail::read_double(col, i);
        }
        p_normals = &assembled;
    } else {
        if (!rOptions.mEstimateNormals)
            rcs_fail("the mesh has no '" + rOptions.mNormalsName +
                     "' point data (run estimate_normals, or allow estimate_normals)");
        if (n < 3)
            rcs_fail("needs at least three samples, not " + std::to_string(n));
        PointNormalsOptions po;
        po.mK = rOptions.mNormalsK;
        estimated = estimate_normals(rMesh, po);
        p_normals = &estimated.mMesh.PointData(kNormalsName);
        result.mEstimatedNormals = true;
    }
    if (!detail::is_float_dtype(p_normals->Dtype()) || p_normals->Size() != n * 3)
        rcs_fail("'" + rOptions.mNormalsName + "' must be a float (n, 3) array");

    // The samples: points with a finite, non-zero normal.
    std::vector<RcsVec> pts;
    std::vector<RcsVec> nrm;
    pts.reserve(n);
    nrm.reserve(n);
    {
        const NDArray& points = rMesh.Points();
        for (std::size_t i = 0; i < n; ++i) {
            RcsVec p{}, v{};
            for (std::size_t a = 0; a < 3; ++a) {
                p[a] = detail::read_double(points, i * 3 + a);
                v[a] = detail::read_double(*p_normals, i * 3 + a);
            }
            if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]))
                rcs_fail("point " + std::to_string(i) + " has a non-finite coordinate");
            const double len = std::sqrt(detail::vec3_dot(v, v));
            if (!std::isfinite(len) || !(len > 0.0))
                continue;
            pts.push_back(p);
            nrm.push_back({v[0] / len, v[1] / len, v[2] / len});
        }
    }
    const std::size_t m = pts.size();
    result.mNumSamples = static_cast<std::int64_t>(m);
    if (m < 3)
        rcs_fail("needs at least three samples with a normal, not " + std::to_string(m));

    const detail::BoundingBox3 box = detail::BoundingBoxUtils::FromPoints(pts);
    const RcsVec& lo = box.mMin;
    const RcsVec ext = detail::BoundingBoxUtils::Extent(box);
    const double extent = detail::BoundingBoxUtils::MaxExtent(box);
    if (!(extent > 0.0))
        rcs_fail("the samples all coincide");

    // The nearest-sample index: buckets of about four samples over the
    // cloud's surface (half its box's area); a speed choice, never the answer.
    double bucket;
    {
        const double area = ext[0] * ext[1] + ext[1] * ext[2] + ext[2] * ext[0];
        bucket = area > 0.0 ? std::sqrt(4.0 * area / static_cast<double>(m)) : extent;
        bucket = std::max(bucket, extent * 1e-9);
    }
    const detail::PointKnn knn(pts, bucket);

    // The spacing: the median distance from a sample to its fourth-nearest
    // other sample (the farthest found, in a smaller cloud). On a regular
    // sampling that is about the nearest-neighbour distance; on a random one,
    // whose nearest neighbours clump, it follows the gaps, which the trim and
    // the band have to bridge.
    {
        constexpr std::size_t kRcsSpacingRank = 4;
        std::vector<double> nn(m, 0.0);
        parallel_for(m, [&](std::size_t i) {
            std::pair<double, std::int64_t> near[kRcsSpacingRank + 1];
            const std::size_t got = knn.KNearest(pts[i], kRcsSpacingRank + 1, near);
            std::size_t others = 0;
            for (std::size_t j = 0; j < got; ++j)
                if (near[j].second != static_cast<std::int64_t>(i)) {
                    nn[i] = std::sqrt(near[j].first);
                    if (++others == kRcsSpacingRank)
                        break;
                }
        });
        nn.erase(std::remove_if(nn.begin(), nn.end(), [](double d) { return !(d > 0.0); }),
                 nn.end());
        if (nn.empty())
            rcs_fail("the samples all coincide");
        const auto mid = nn.begin() + static_cast<std::ptrdiff_t>(nn.size() / 2);
        std::nth_element(nn.begin(), mid, nn.end());
        result.mSpacing = *mid;
    }
    const double h = rOptions.mCellSize > 0.0   ? rOptions.mCellSize
                     : rOptions.mResolution > 0 ? extent / static_cast<double>(rOptions.mResolution)
                                                : result.mSpacing;
    const double sigma = std::max(result.mSpacing, h);
    const double max_dist = rOptions.mMaxDistance > 0.0 ? rOptions.mMaxDistance : 2.0 * sigma;
    result.mCellSize = h;
    result.mMaxDistance = max_dist;

    // The lattice: the samples' box grown by the band's reach and one cell.
    const auto reach = static_cast<std::int64_t>(std::ceil(sigma / h)) + 1;
    RcsVec origin{};
    std::array<std::int64_t, 3> dims{};
    for (std::size_t a = 0; a < 3; ++a) {
        const double cells = std::ceil(ext[a] / h) + 2.0 * static_cast<double>(reach + 1);
        if (!(cells < 2097152.0))  // 2^21 per axis keeps the linear ids in an int64
            rcs_fail("the cell size " + std::to_string(h) +
                     " is too small for the cloud's extent " + std::to_string(extent));
        dims[a] = static_cast<std::int64_t>(cells);
        origin[a] = lo[a] - static_cast<double>(reach + 1) * h;
    }
    const auto cell_of = [&](const RcsVec& rP, std::size_t A) {
        const auto c = static_cast<std::int64_t>(std::floor((rP[A] - origin[A]) / h));
        return std::clamp<std::int64_t>(c, 0, dims[A] - 1);
    };

    // The band: every sample's cell, grown by the reach along each axis.
    std::vector<std::int64_t> band(m);
    for (std::size_t i = 0; i < m; ++i)
        band[i] =
            (cell_of(pts[i], 2) * dims[1] + cell_of(pts[i], 1)) * dims[0] + cell_of(pts[i], 0);
    rcs_unique(band);
    for (int axis = 0; axis < 3; ++axis) {
        rcs_dilate(band, dims, axis, reach, reach);
        if (static_cast<std::int64_t>(band.size()) > rOptions.mMaxCells)
            rcs_fail("the band about the samples holds over " + std::to_string(band.size()) +
                     " cells, more than max_cells (" + std::to_string(rOptions.mMaxCells) +
                     "); raise the cell size (now " + std::to_string(h) + ") or max_cells");
    }

    // Its nodes: each cell's eight corners, by growing the cells' lo corners
    // one node up along each axis on the node lattice.
    const std::array<std::int64_t, 3> ndims{dims[0] + 1, dims[1] + 1, dims[2] + 1};
    const auto node_of_cell = [&](std::int64_t C) {
        const std::int64_t i = C % dims[0], j = (C / dims[0]) % dims[1],
                           k = C / (dims[0] * dims[1]);
        return (k * ndims[1] + j) * ndims[0] + i;
    };
    std::vector<std::int64_t> nodes(band.size());
    for (std::size_t c = 0; c < band.size(); ++c)
        nodes[c] = node_of_cell(band[c]);
    for (int axis = 0; axis < 3; ++axis)
        rcs_dilate(nodes, ndims, axis, 0, 1);
    const std::size_t nn = nodes.size();
    result.mNumLatticePoints = static_cast<std::int64_t>(nn);
    const auto node_pos = [&](std::int64_t V) {
        const std::int64_t i = V % ndims[0], j = (V / ndims[0]) % ndims[1],
                           k = V / (ndims[0] * ndims[1]);
        return RcsVec{origin[0] + static_cast<double>(i) * h,
                      origin[1] + static_cast<double>(j) * h,
                      origin[2] + static_cast<double>(k) * h};
    };

    // The field at each node, or undefined.
    const std::size_t k =
        static_cast<std::size_t>(std::min<std::int64_t>(rOptions.mK, static_cast<std::int64_t>(m)));
    std::vector<double> value(nn, 0.0);
    std::vector<char> defined(nn, 0);
    const double inv_s2 = 1.0 / (sigma * sigma);
    const double max_d2 = max_dist * max_dist;
    // The cutter drops primitives below 1e-10 of its box's diagonal squared;
    // crossings at least this far from every node keep a triangle above it.
    const double nudge = std::max(1e-3 * h, 1e-4 * extent);
    parallel_for(nn, [&](std::size_t v) {
        const RcsVec x = node_pos(nodes[v]);
        std::vector<std::pair<double, std::int64_t>> near(k);
        const std::size_t got = knn.KNearest(x, k, near.data());
        const RcsVec& p0 = pts[static_cast<std::size_t>(near[0].second)];
        const RcsVec& n0 = nrm[static_cast<std::size_t>(near[0].second)];
        const RcsVec d0{x[0] - p0[0], x[1] - p0[1], x[2] - p0[2]};
        const double f0 = detail::vec3_dot(d0, n0);
        // Hoppe's test: the projection onto the nearest plane, x - f0 n0,
        // lies within max_dist of p0 (its offset from p0 is d0 - f0 n0).
        const RcsVec t{d0[0] - f0 * n0[0], d0[1] - f0 * n0[1], d0[2] - f0 * n0[2]};
        if (detail::vec3_dot(t, t) > max_d2)
            return;
        double num = 0.0, den = 0.0;
        for (std::size_t j = 0; j < got; ++j) {
            const auto id = static_cast<std::size_t>(near[j].second);
            const RcsVec d{x[0] - pts[id][0], x[1] - pts[id][1], x[2] - pts[id][2]};
            const double w = std::exp(-near[j].first * inv_s2);
            num += w * detail::vec3_dot(d, nrm[id]);
            den += w;
        }
        // Keep the level set off the node itself: a value near zero would
        // gather every crossing about the node into a speck the cutter drops
        // as collapsed, leaving a hole.
        const double f = den > 0.0 ? num / den : f0;
        value[v] = std::abs(f) >= nudge ? f : (f < 0.0 ? -nudge : nudge);
        defined[v] = 1;
    });
    result.mNumUndefined =
        static_cast<std::int64_t>(std::count(defined.begin(), defined.end(), char{0}));

    // The cells the surface crosses: every corner defined and not all on one
    // side (zero counting as positive, as the cutter does). Each keeps its
    // corners as node indices, in the lattice's hexahedron winding.
    const std::array<std::int64_t, 8> corner{0,
                                             1,
                                             ndims[0] + 1,
                                             ndims[0],
                                             ndims[0] * ndims[1],
                                             ndims[0] * ndims[1] + 1,
                                             ndims[0] * ndims[1] + ndims[0] + 1,
                                             ndims[0] * ndims[1] + ndims[0]};
    std::vector<std::array<std::int64_t, 8>> conn(band.size());
    std::vector<char> active(band.size(), 0);
    parallel_for(band.size(), [&](std::size_t c) {
        const std::int64_t base = node_of_cell(band[c]);
        int pos = 0;
        for (std::size_t q = 0; q < 8; ++q) {
            const auto it = std::lower_bound(nodes.begin(), nodes.end(), base + corner[q]);
            const auto id = static_cast<std::size_t>(it - nodes.begin());
            if (!defined[id])
                return;
            conn[c][q] = static_cast<std::int64_t>(id);
            pos += value[id] >= 0.0 ? 1 : 0;
        }
        active[c] = pos > 0 && pos < 8 ? 1 : 0;
    });

    // The crossed cells as a hexahedron mesh over their own nodes.
    std::vector<std::int64_t> remap(nn, -1);
    std::vector<std::int64_t> used;
    std::size_t num_active = 0;
    for (std::size_t c = 0; c < band.size(); ++c) {
        if (!active[c])
            continue;
        ++num_active;
        for (const std::int64_t id : conn[c])
            if (remap[static_cast<std::size_t>(id)] < 0) {
                remap[static_cast<std::size_t>(id)] = static_cast<std::int64_t>(used.size());
                used.push_back(id);
            }
    }
    Mesh out;
    std::vector<RcsVec> tri_pts;
    std::vector<std::array<std::int64_t, 3>> tris;
    if (num_active > 0) {
        Mesh hexes;
        NDArray hp(DType::Float64, {used.size(), std::size_t{3}});
        double* php = hp.As<double>();
        std::vector<double> hv(used.size());
        for (std::size_t u = 0; u < used.size(); ++u) {
            const RcsVec x = node_pos(nodes[static_cast<std::size_t>(used[u])]);
            std::copy(x.begin(), x.end(), php + u * 3);
            hv[u] = value[static_cast<std::size_t>(used[u])];
        }
        hexes.AssignPoints(std::move(hp));
        NDArray hc(DType::Int64, {num_active, std::size_t{8}});
        std::int64_t* phc = hc.As<std::int64_t>();
        std::size_t row = 0;
        for (std::size_t c = 0; c < band.size(); ++c) {
            if (!active[c])
                continue;
            for (std::size_t q = 0; q < 8; ++q)
                phc[row * 8 + q] = remap[static_cast<std::size_t>(conn[c][q])];
            ++row;
        }
        hexes.AddCellBlock("hexahedron", std::move(hc));

        const detail::MarchingInput input = detail::marching_prepare(hexes);
        if (input.Simp().NumPoints() != used.size())
            throw std::runtime_error(std::string(kRcsPrefix) +
                                     "internal error: the split moved the nodes");
        detail::MarchingOptions mo;
        mo.mOrientation = detail::MarchingOrientation::FieldGradient;
        const Mesh cut = detail::marching_cut(input, hv, mo);

        // Weld coincident crossings (a node exactly on the level set is
        // reached from several edges): equal coordinates take the lowest
        // index, and the survivors keep their order.
        const std::size_t cp = cut.NumPoints();
        std::vector<RcsVec> cpts(cp);
        {
            const NDArray& p = cut.Points();
            for (std::size_t i = 0; i < cp; ++i)
                for (std::size_t a = 0; a < 3; ++a)
                    cpts[i][a] = detail::read_double(p, i * 3 + a);
        }
        std::vector<std::int64_t> order(cp);
        std::iota(order.begin(), order.end(), std::int64_t{0});
        std::sort(order.begin(), order.end(), [&](std::int64_t A, std::int64_t B) {
            const RcsVec& a = cpts[static_cast<std::size_t>(A)];
            const RcsVec& b = cpts[static_cast<std::size_t>(B)];
            return a != b ? a < b : A < B;
        });
        std::vector<std::int64_t> rep(cp);
        for (std::size_t r = 0; r < cp; ++r) {
            const auto i = static_cast<std::size_t>(order[r]);
            rep[i] = r > 0 && cpts[static_cast<std::size_t>(order[r - 1])] == cpts[i]
                         ? rep[static_cast<std::size_t>(order[r - 1])]
                         : static_cast<std::int64_t>(i);
        }
        std::vector<std::int64_t> renum(cp, -1);
        for (std::size_t i = 0; i < cp; ++i)
            if (rep[i] == static_cast<std::int64_t>(i)) {
                renum[i] = static_cast<std::int64_t>(tri_pts.size());
                tri_pts.push_back(cpts[i]);
            }
        const auto id_of = [&](const NDArray& rC, std::size_t I) {
            return renum[static_cast<std::size_t>(
                rep[static_cast<std::size_t>(detail::read_int(rC, I))])];
        };
        const auto emit = [&](std::int64_t A, std::int64_t B, std::int64_t C) {
            if (A != B && B != C && C != A)
                tris.push_back({A, B, C});
        };
        for (std::size_t b = 0; b < cut.NumCellBlocks(); ++b) {
            const auto block = cut.Cells(b);
            const NDArray& c = block.Conn();
            const std::size_t rows = block.NumCells();
            if (block.Type() == "triangle") {
                for (std::size_t r = 0; r < rows; ++r)
                    emit(id_of(c, r * 3), id_of(c, r * 3 + 1), id_of(c, r * 3 + 2));
            } else if (block.Type() == "quad") {
                for (std::size_t r = 0; r < rows; ++r) {
                    const std::int64_t q0 = id_of(c, r * 4), q1 = id_of(c, r * 4 + 1),
                                       q2 = id_of(c, r * 4 + 2), q3 = id_of(c, r * 4 + 3);
                    emit(q0, q1, q2);
                    emit(q0, q2, q3);
                }
            }
        }
    }

    // Points no triangle uses (a welded-away sliver's) are dropped.
    std::vector<std::int64_t> keep(tri_pts.size(), -1);
    std::vector<RcsVec> final_pts;
    for (auto& t : tris)
        for (std::int64_t& id : t) {
            auto& r = keep[static_cast<std::size_t>(id)];
            if (r < 0) {
                r = static_cast<std::int64_t>(final_pts.size());
                final_pts.push_back(tri_pts[static_cast<std::size_t>(id)]);
            }
            id = r;
        }
    NDArray op(DType::Float64, {final_pts.size(), std::size_t{3}});
    double* pop = op.As<double>();
    for (std::size_t i = 0; i < final_pts.size(); ++i)
        std::copy(final_pts[i].begin(), final_pts[i].end(), pop + i * 3);
    out.AssignPoints(std::move(op));
    NDArray oc(DType::Int64, {tris.size(), std::size_t{3}});
    std::int64_t* poc = oc.As<std::int64_t>();
    for (std::size_t t = 0; t < tris.size(); ++t)
        std::copy(tris[t].begin(), tris[t].end(), poc + t * 3);
    out.AddCellBlock("triangle", std::move(oc));
    const double unit = detail::mesh_length_unit(rMesh);
    if (unit > 0.0)
        detail::set_mesh_length_unit(out, unit);

    // Edge-connected pieces.
    std::vector<std::int64_t> parent(final_pts.size());
    std::iota(parent.begin(), parent.end(), std::int64_t{0});
    for (const auto& t : tris)
        for (std::size_t e = 1; e < 3; ++e) {
            const std::int64_t a = rcs_find(parent, t[0]), b = rcs_find(parent, t[e]);
            if (a != b)
                parent[static_cast<std::size_t>(std::max(a, b))] = std::min(a, b);
        }
    for (std::size_t i = 0; i < final_pts.size(); ++i)
        if (rcs_find(parent, static_cast<std::int64_t>(i)) == static_cast<std::int64_t>(i))
            ++result.mNumComponents;
    result.mNumTriangles = static_cast<std::int64_t>(tris.size());
    result.mMesh = std::move(out);
    return result;
}

}  // namespace meshioplusplus
