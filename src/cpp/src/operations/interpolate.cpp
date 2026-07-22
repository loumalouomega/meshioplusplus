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
// Cross-mesh field transfer: sample a source mesh's data arrays onto a target
// mesh (nearest source point, or barycentric in a simplexified source), built
// entirely through the uniform mesh API so it compiles under every mesh
// backend. See operations/interpolate.hpp for the contract, in particular the
// determinism rules the pure-numpy twin (_interpolate.py) replicates
// expression for expression.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// Store `v` at flat element index `flat` of `rOut`, casting to `rOut`'s dtype
// (merge_store_double's twin; anon-namespace helpers cannot be shared).
void interp_store_double(NDArray& rOut, std::size_t flat, double v) {
    switch (rOut.Dtype()) {
        case DType::Float32:
            rOut.As<float>()[flat] = static_cast<float>(v);
            return;
        case DType::Float64:
            rOut.As<double>()[flat] = v;
            return;
        case DType::Int8:
            rOut.As<std::int8_t>()[flat] = static_cast<std::int8_t>(v);
            return;
        case DType::Int16:
            rOut.As<std::int16_t>()[flat] = static_cast<std::int16_t>(v);
            return;
        case DType::Int32:
            rOut.As<std::int32_t>()[flat] = static_cast<std::int32_t>(v);
            return;
        case DType::Int64:
            rOut.As<std::int64_t>()[flat] = static_cast<std::int64_t>(v);
            return;
        case DType::UInt8:
            rOut.As<std::uint8_t>()[flat] = static_cast<std::uint8_t>(v);
            return;
        case DType::UInt16:
            rOut.As<std::uint16_t>()[flat] = static_cast<std::uint16_t>(v);
            return;
        case DType::UInt32:
            rOut.As<std::uint32_t>()[flat] = static_cast<std::uint32_t>(v);
            return;
        case DType::UInt64:
            rOut.As<std::uint64_t>()[flat] = static_cast<std::uint64_t>(v);
            return;
    }
}

// Total cell count across all blocks (the global cell numbering size).
std::size_t interp_total_cells(const Mesh& rMesh) {
    std::size_t total = 0;
    for (const auto cb : rMesh.CellRange())
        total += cb.NumCells();
    return total;
}

// Per-cell centroids in global cell order, always 3 doubles per cell (unused
// axes zero) — the exact accumulation order of partition.cpp's
// partition_centroids (anon-namespace, so duplicated here; keep the two and
// _partition.py's _centroids in step, the numpy parity test pins this).
std::vector<double> interp_centroids(const Mesh& rMesh, std::size_t total) {
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    const std::size_t dim = std::min<std::size_t>(pdim, 3);

    std::vector<double> cent(3 * total, 0.0);
    std::size_t base = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::size_t ncells = cb.NumCells();
        if (cb.IsPolyhedron()) {
            parallel_for(ncells, [&, base](std::size_t c) {
                std::vector<std::int64_t> nodes;
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    const auto face = cb.Face(c, f);
                    nodes.insert(nodes.end(), face.first, face.first + face.second);
                }
                std::sort(nodes.begin(), nodes.end());
                nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
                double* out = cent.data() + 3 * (base + c);
                for (std::int64_t node : nodes)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] +=
                            detail::read_double(points, static_cast<std::size_t>(node) * pdim + d);
                if (!nodes.empty())
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] /= static_cast<double>(nodes.size());
            });
        } else if (cb.IsRagged()) {
            parallel_for(ncells, [&, base](std::size_t c) {
                const std::int64_t* row = cb.Row(c);
                const std::size_t nn = cb.RowSize(c);
                double* out = cent.data() + 3 * (base + c);
                for (std::size_t k = 0; k < nn; ++k)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] += detail::read_double(points,
                                                      static_cast<std::size_t>(row[k]) * pdim + d);
                if (nn > 0)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] /= static_cast<double>(nn);
            });
        } else {
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            parallel_for(ncells, [&, base](std::size_t c) {
                double* out = cent.data() + 3 * (base + c);
                for (std::size_t k = 0; k < npc; ++k) {
                    const std::size_t node =
                        static_cast<std::size_t>(detail::read_int(conn, c * npc + k));
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] += detail::read_double(points, node * pdim + d);
                }
                if (npc > 0)
                    for (std::size_t d = 0; d < dim; ++d)
                        out[d] /= static_cast<double>(npc);
            });
        }
        base += ncells;
    }
    return cent;
}

// The grid cell size for `n` items spread over the bbox `[pLo, pHi]`: the
// largest axis extent divided by the smallest integer R with R*R*R >= n.
// Deliberately no cbrt/pow — the integer loop is exact in both C++ and Python,
// and the single division is IEEE-identical, which is what keeps the numpy
// fallback's grid byte-compatible.
double interp_cell_size(const double* pLo, const double* pHi, std::size_t n) {
    double max_extent = 0.0;
    for (int d = 0; d < 3; ++d) {
        const double e = pHi[d] - pLo[d];
        if (e > max_extent)
            max_extent = e;
    }
    std::int64_t r = 1;
    while (r * r * r < static_cast<std::int64_t>(n))
        ++r;
    return max_extent > 0.0 ? max_extent / static_cast<double>(r) : 1.0;
}

// Serial bbox over `n` coordinates supplied by `coord(i) -> Vec3` (min/max are
// pure comparisons, so the order cannot matter). Callers guarantee n >= 1.
template <class TCoord>
void interp_bbox(std::size_t n, TCoord&& coord, double* pLo, double* pHi) {
    const detail::Vec3 first = coord(0);
    for (int d = 0; d < 3; ++d)
        pLo[d] = pHi[d] = first[d];
    for (std::size_t i = 1; i < n; ++i) {
        const detail::Vec3 c = coord(i);
        for (int d = 0; d < 3; ++d) {
            if (c[d] < pLo[d])
                pLo[d] = c[d];
            if (c[d] > pHi[d])
                pHi[d] = c[d];
        }
    }
}

// Fill `rGrid` with one id per coordinate: parallel key fill into disjoint
// slots, then a serial ascending insert (surface.cpp's phase-split idiom, so
// every bucket vector is ascending and scans over it are deterministic).
template <class TCoord>
void interp_fill_grid(detail::SpatialGrid& rGrid, std::size_t n, TCoord&& coord) {
    std::vector<detail::GridKey> keys(n);
    parallel_for(n, [&](std::size_t i) {
        const detail::Vec3 c = coord(i);
        keys[i] = rGrid.KeyOf(c.data());
    });
    for (std::size_t i = 0; i < n; ++i)
        rGrid.Insert(keys[i], static_cast<std::int64_t>(i));
}

// Nearest id in `rGrid` to `rP` by expanding-shell search: any candidate in a
// cell at Chebyshev key radius r >= 1 lies at least (r-1)*cell away, so once
// ((r-1)*cell)^2 exceeds the best squared distance no farther shell can win;
// the occupied-key bbox bounds the expansion (the grid is never empty here).
// Tie-break: smallest squared distance, then lowest id — enforced by the
// strict comparison plus the ascending id order within every bucket.
template <class TCoord>
std::int64_t interp_nearest(const detail::SpatialGrid& rGrid, TCoord&& coord,
                            const detail::Vec3& rP) {
    const detail::GridKey k = rGrid.KeyOf(rP.data());
    const detail::GridKey& klo = rGrid.OccupiedLo();
    const detail::GridKey& khi = rGrid.OccupiedHi();
    auto axis_cap = [](std::int64_t v, std::int64_t lo, std::int64_t hi) {
        const std::int64_t a = v > lo ? v - lo : lo - v;
        const std::int64_t b = v > hi ? v - hi : hi - v;
        return a > b ? a : b;
    };
    std::int64_t rcap = axis_cap(k.x, klo.x, khi.x);
    rcap = std::max(rcap, axis_cap(k.y, klo.y, khi.y));
    rcap = std::max(rcap, axis_cap(k.z, klo.z, khi.z));
    // Shells below the Chebyshev distance to the occupied bbox are provably
    // empty; starting there skips them without affecting any candidate.
    auto axis_start = [](std::int64_t v, std::int64_t lo, std::int64_t hi) {
        if (v < lo)
            return lo - v;
        if (v > hi)
            return v - hi;
        return std::int64_t{0};
    };
    std::int64_t r0 = axis_start(k.x, klo.x, khi.x);
    r0 = std::max(r0, axis_start(k.y, klo.y, khi.y));
    r0 = std::max(r0, axis_start(k.z, klo.z, khi.z));

    const double cell = rGrid.CellSize();
    double best_d2 = std::numeric_limits<double>::infinity();
    std::int64_t best = -1;
    for (std::int64_t r = r0;; ++r) {
        if (best >= 0 && r >= 1) {
            const double ring = static_cast<double>(r - 1) * cell;
            if (ring * ring > best_d2)
                break;
        }
        if (r > rcap)
            break;
        rGrid.ForEachInShell(k, r, [&](const std::vector<std::int64_t>& rIds) {
            for (std::int64_t id : rIds) {
                const detail::Vec3 q = coord(id);
                double d2 = 0.0;
                for (int d = 0; d < 3; ++d) {
                    const double delta = rP[d] - q[d];
                    d2 += delta * delta;
                }
                if (d2 < best_d2 || (d2 == best_d2 && id < best)) {
                    best_d2 = d2;
                    best = id;
                }
            }
        });
    }
    return best;
}

// Barycentric containment tolerance: dimensionless (the weights are volume
// *ratios*), hence scale-invariant. The same literal appears in the numpy
// twin.
constexpr double INTERP_BARY_EPS = 1e-12;

// Barycentric weights of `rP` in the tetra `pNodes[0..3]`; true iff contained
// (all weights >= -INTERP_BARY_EPS). A degenerate tetra (exactly zero volume)
// contains nothing.
bool interp_tet_contains(const NDArray& rPts, std::size_t pdim, const std::int64_t* pNodes,
                         const detail::Vec3& rP, double* pW) {
    const detail::Vec3 a = detail::read_point(rPts, pdim, pNodes[0]);
    const detail::Vec3 b = detail::read_point(rPts, pdim, pNodes[1]);
    const detail::Vec3 c = detail::read_point(rPts, pdim, pNodes[2]);
    const detail::Vec3 d = detail::read_point(rPts, pdim, pNodes[3]);
    const double v = detail::triple_product(detail::vec3_sub(b, a), detail::vec3_sub(c, a),
                                            detail::vec3_sub(d, a));
    if (v == 0.0)
        return false;
    const double wa = detail::triple_product(detail::vec3_sub(b, rP), detail::vec3_sub(c, rP),
                                             detail::vec3_sub(d, rP)) /
                      v;
    const double wb = detail::triple_product(detail::vec3_sub(rP, a), detail::vec3_sub(c, a),
                                             detail::vec3_sub(d, a)) /
                      v;
    const double wc = detail::triple_product(detail::vec3_sub(b, a), detail::vec3_sub(rP, a),
                                             detail::vec3_sub(d, a)) /
                      v;
    const double wd = detail::triple_product(detail::vec3_sub(b, a), detail::vec3_sub(c, a),
                                             detail::vec3_sub(rP, a)) /
                      v;
    if (wa < -INTERP_BARY_EPS || wb < -INTERP_BARY_EPS || wc < -INTERP_BARY_EPS ||
        wd < -INTERP_BARY_EPS)
        return false;
    pW[0] = wa;
    pW[1] = wb;
    pW[2] = wc;
    pW[3] = wd;
    return true;
}

// Twice the signed area of the triangle p0-p1-p2 **in the xy-plane** (the
// documented planar assumption of the triangle path; z is ignored).
double interp_area2(const detail::Vec3& rP0, const detail::Vec3& rP1, const detail::Vec3& rP2) {
    return (rP1[0] - rP0[0]) * (rP2[1] - rP0[1]) - (rP1[1] - rP0[1]) * (rP2[0] - rP0[0]);
}

// Barycentric weights of `rP` in the (xy-plane) triangle `pNodes[0..2]`; true
// iff contained.
bool interp_tri_contains(const NDArray& rPts, std::size_t pdim, const std::int64_t* pNodes,
                         const detail::Vec3& rP, double* pW) {
    const detail::Vec3 a = detail::read_point(rPts, pdim, pNodes[0]);
    const detail::Vec3 b = detail::read_point(rPts, pdim, pNodes[1]);
    const detail::Vec3 c = detail::read_point(rPts, pdim, pNodes[2]);
    const double v = interp_area2(a, b, c);
    if (v == 0.0)
        return false;
    const double wa = interp_area2(rP, b, c) / v;
    const double wb = interp_area2(a, rP, c) / v;
    const double wc = interp_area2(a, b, rP) / v;
    if (wa < -INTERP_BARY_EPS || wb < -INTERP_BARY_EPS || wc < -INTERP_BARY_EPS)
        return false;
    pW[0] = wa;
    pW[1] = wb;
    pW[2] = wc;
    return true;
}

// The output shape for `nrows` sampled rows mirroring `rSrc`'s trailing
// dimensions (a 1-D scalar array stays 1-D, the merge convention).
std::vector<std::size_t> interp_out_shape(const NDArray& rSrc, std::size_t nrows) {
    std::vector<std::size_t> shape = rSrc.Shape();
    if (shape.empty())
        shape = {nrows};
    else
        shape[0] = nrows;
    return shape;
}

// Gather rows of `rSrc` by index, dtype-preserving (a raw byte gather — never
// through double, so int64 payloads beyond 2^53 survive bit-for-bit).
NDArray interp_gather_rows(const NDArray& rSrc, const std::vector<std::int64_t>& rIdx) {
    NDArray dst = NDArray::Uninit(rSrc.Dtype(), interp_out_shape(rSrc, rIdx.size()));
    const std::size_t src_rows = detail::rows(rSrc);
    if (rIdx.empty() || src_rows == 0)
        return dst;
    const std::size_t rb = rSrc.Nbytes() / src_rows;
    const std::byte* src = rSrc.Data();
    std::byte* out = dst.Data();
    parallel_for_bw(rIdx.size(), [&](std::size_t n) {
        std::memcpy(out + n * rb, src + static_cast<std::size_t>(rIdx[n]) * rb, rb);
    });
    return dst;
}

// The standard unknown-key diagnostic, covering both source locations.
std::string interp_unknown_array_message(const Mesh& rSource, const std::string& rName) {
    auto join = [](const std::vector<std::string>& rNames) -> std::string {
        if (rNames.empty())
            return "none";
        std::string s;
        for (std::size_t i = 0; i < rNames.size(); ++i) {
            if (i > 0)
                s += ", ";
            s += rNames[i];
        }
        return s;
    };
    return "meshio++: interpolate: no source point_data or cell_data array named '" + rName +
           "' (available point_data: " + join(data_names(rSource, DataLocation::Point)) +
           "; cell_data: " + join(data_names(rSource, DataLocation::Cell)) + ")";
}

// The output name for a transferred array under the conflict policy, resolved
// against the *original* target before any sampling happens.
std::string interp_resolve_name(const Mesh& rTarget, bool isPointData, const std::string& rName,
                                InterpolateConflict Conflict) {
    const bool exists = isPointData ? rTarget.HasPointData(rName) : rTarget.HasCellData(rName);
    if (!exists)
        return rName;
    switch (Conflict) {
        case InterpolateConflict::Overwrite:
            return rName;
        case InterpolateConflict::Suffix: {
            const std::string suffixed = rName + "_interp";
            const bool taken =
                isPointData ? rTarget.HasPointData(suffixed) : rTarget.HasCellData(suffixed);
            if (taken)
                throw std::invalid_argument("meshio++: interpolate: array '" + rName +
                                            "' already exists on the target and so does '" +
                                            suffixed + "' (on_conflict='suffix')");
            return suffixed;
        }
        case InterpolateConflict::Error:
        default:
            throw std::invalid_argument("meshio++: interpolate: array '" + rName +
                                        "' already exists on the target (on_conflict='error'; "
                                        "pass 'overwrite' or 'suffix')");
    }
}

}  // namespace

InterpolateMethod interpolate_method_from_name(const std::string& rName) {
    if (rName == "nearest")
        return InterpolateMethod::Nearest;
    if (rName == "barycentric")
        return InterpolateMethod::Barycentric;
    throw std::invalid_argument("meshio++: interpolate: unknown method '" + rName +
                                "' (expected 'nearest' or 'barycentric')");
}

InterpolateConflict interpolate_conflict_from_name(const std::string& rName) {
    if (rName == "error")
        return InterpolateConflict::Error;
    if (rName == "overwrite")
        return InterpolateConflict::Overwrite;
    if (rName == "suffix")
        return InterpolateConflict::Suffix;
    throw std::invalid_argument("meshio++: interpolate: unknown on_conflict '" + rName +
                                "' (expected 'error', 'overwrite' or 'suffix')");
}

Mesh interpolate(const Mesh& rSource, const Mesh& rTarget, const InterpolateOptions& rOptions) {
    if (rSource.NumPoints() == 0)
        throw std::invalid_argument("meshio++: interpolate: the source mesh has no points");

    // --- resolve which arrays transfer (validation before any work) ---------
    std::vector<std::string> point_names;
    std::vector<std::string> cell_names;
    if (rOptions.mArrays.empty()) {
        point_names = rSource.PointDataNames();  // sorted on every backend
    } else {
        for (const std::string& name : rOptions.mArrays) {
            bool any = false;
            if (rSource.HasPointData(name)) {
                point_names.push_back(name);
                any = true;
            }
            if (rSource.HasCellData(name)) {
                cell_names.push_back(name);
                any = true;
            }
            if (!any)
                throw std::invalid_argument(interp_unknown_array_message(rSource, name));
        }
    }

    for (const std::string& name : cell_names) {
        if (rSource.CellDataNumBlocks(name) != rSource.NumCellBlocks())
            throw std::invalid_argument("meshio++: interpolate: source cell_data '" + name +
                                        "' does not have one array per cell block");
        const std::size_t nblocks = rSource.CellDataNumBlocks(name);
        for (std::size_t b = 1; b < nblocks; ++b)
            if (data_num_components(rSource.CellData(name, b)) !=
                data_num_components(rSource.CellData(name, 0)))
                throw std::invalid_argument("meshio++: interpolate: source cell_data '" + name +
                                            "' has inconsistent components across blocks");
    }

    std::vector<std::string> point_out(point_names.size());
    for (std::size_t i = 0; i < point_names.size(); ++i)
        point_out[i] = interp_resolve_name(rTarget, true, point_names[i], rOptions.mOnConflict);
    std::vector<std::string> cell_out(cell_names.size());
    for (std::size_t i = 0; i < cell_names.size(); ++i)
        cell_out[i] = interp_resolve_name(rTarget, false, cell_names[i], rOptions.mOnConflict);

    // KRATOS's Mesh is not copy-constructible — clone through the uniform API.
    Mesh out = detail::clone_mesh(rTarget);

    const NDArray& tpts = rTarget.Points();
    const std::size_t tdim = rTarget.PointDim();
    const std::size_t nt = rTarget.NumPoints();

    // --- point_data: sampled at the target's points -------------------------
    if (!point_names.empty()) {
        if (rOptions.mMethod == InterpolateMethod::Nearest) {
            const NDArray& spts = rSource.Points();
            const std::size_t sdim = rSource.PointDim();
            const std::size_t ns = rSource.NumPoints();
            auto scoord = [&](std::int64_t i) { return detail::read_point(spts, sdim, i); };
            double lo[3], hi[3];
            interp_bbox(
                ns, [&](std::size_t i) { return scoord(static_cast<std::int64_t>(i)); }, lo, hi);
            detail::SpatialGrid grid(interp_cell_size(lo, hi, ns));
            interp_fill_grid(grid, ns,
                             [&](std::size_t i) { return scoord(static_cast<std::int64_t>(i)); });

            std::vector<std::int64_t> nearest(nt);
            parallel_for(nt, [&](std::size_t t) {
                const detail::Vec3 p = detail::read_point(tpts, tdim, static_cast<std::int64_t>(t));
                nearest[t] = interp_nearest(grid, scoord, p);
            });
            for (std::size_t i = 0; i < point_names.size(); ++i)
                out.AddPointData(point_out[i],
                                 interp_gather_rows(rSource.PointData(point_names[i]), nearest));
        } else {
            // Simplexify so containment/weights are exact and uniform; points
            // are untouched for a linear source, and point_data rides through
            // (row-subset under the higher-order linearization), so sampling
            // the simplexified mesh's own arrays is always index-consistent.
            ConvertCellsOptions cco;
            cco.mMode = ConvertCellsMode::Simplexify;
            const ConvertCellsResult simp = convert_cells(rSource, cco);
            const Mesh& smesh = simp.mMesh;

            bool has_tet = false;
            bool has_tri = false;
            for (const auto cb : smesh.CellRange()) {
                if (cb.NumCells() == 0)
                    continue;
                if (cb.Type() == "tetra")
                    has_tet = true;
                else if (cb.Type() == "triangle")
                    has_tri = true;
            }
            if (!has_tet && !has_tri)
                throw std::invalid_argument(
                    "meshio++: interpolate: barycentric needs a source with triangle or "
                    "tetrahedron cells after simplexification");
            const std::size_t nv = has_tet ? 4 : 3;
            const std::string want = has_tet ? "tetra" : "triangle";

            // Flat simplex connectivity, ascending block-then-cell order (the
            // simplex index every tie-break below refers to).
            std::vector<std::int64_t> sconn;
            for (const auto cb : smesh.CellRange()) {
                if (cb.Type() != want || cb.NumCells() == 0)
                    continue;
                const NDArray& conn = cb.Conn();
                const std::size_t n = cb.NumCells() * nv;
                const std::size_t base = sconn.size();
                sconn.resize(base + n);
                for (std::size_t j = 0; j < n; ++j)
                    sconn[base + j] = detail::read_int(conn, j);
            }
            const std::size_t nsimp = sconn.size() / nv;

            const NDArray& spts = smesh.Points();
            const std::size_t sdim = smesh.PointDim();
            const std::size_t ns = smesh.NumPoints();
            auto scoord = [&](std::int64_t i) { return detail::read_point(spts, sdim, i); };

            double lo[3], hi[3];
            interp_bbox(
                ns, [&](std::size_t i) { return scoord(static_cast<std::int64_t>(i)); }, lo, hi);
            detail::SpatialGrid grid(interp_cell_size(lo, hi, nsimp));

            // Per-simplex quantized bbox keys: parallel fill, serial insert.
            std::vector<std::pair<detail::GridKey, detail::GridKey>> boxes(nsimp);
            parallel_for(nsimp, [&](std::size_t s) {
                double blo[3], bhi[3];
                interp_bbox(nv, [&](std::size_t k) { return scoord(sconn[s * nv + k]); }, blo, bhi);
                boxes[s] = {grid.KeyOf(blo), grid.KeyOf(bhi)};
            });
            for (std::size_t s = 0; s < nsimp; ++s)
                grid.InsertBox(boxes[s].first, boxes[s].second, static_cast<std::int64_t>(s));

            // Containment: first containing simplex in ascending index (the
            // bucket vectors are ascending by construction).
            std::vector<std::int64_t> simplex_of(nt, -1);
            std::vector<std::array<double, 4>> wts(nt);
            parallel_for(nt, [&](std::size_t t) {
                const detail::Vec3 p = detail::read_point(tpts, tdim, static_cast<std::int64_t>(t));
                const std::vector<std::int64_t>* p_cand = grid.Find(grid.KeyOf(p.data()));
                if (p_cand == nullptr)
                    return;
                double w[4] = {0.0, 0.0, 0.0, 0.0};
                for (std::int64_t s : *p_cand) {
                    const std::int64_t* nodes = &sconn[static_cast<std::size_t>(s) * nv];
                    const bool inside = nv == 4 ? interp_tet_contains(spts, sdim, nodes, p, w)
                                                : interp_tri_contains(spts, sdim, nodes, p, w);
                    if (inside) {
                        simplex_of[t] = s;
                        for (std::size_t k = 0; k < nv; ++k)
                            wts[t][k] = w[k];
                        break;
                    }
                }
            });

            std::size_t n_outside = 0;
            for (std::size_t t = 0; t < nt; ++t)
                if (simplex_of[t] < 0)
                    ++n_outside;

            std::vector<std::int64_t> near_of;
            if (n_outside > 0 && rOptions.mExtrapolate) {
                detail::SpatialGrid pgrid(interp_cell_size(lo, hi, ns));
                interp_fill_grid(
                    pgrid, ns, [&](std::size_t i) { return scoord(static_cast<std::int64_t>(i)); });
                near_of.assign(nt, -1);
                parallel_for(nt, [&](std::size_t t) {
                    if (simplex_of[t] >= 0)
                        return;
                    const detail::Vec3 p =
                        detail::read_point(tpts, tdim, static_cast<std::int64_t>(t));
                    near_of[t] = interp_nearest(pgrid, scoord, p);
                });
            } else if (n_outside > 0) {
                log::warn(
                    "meshio++: interpolate: {} of {} target points lie outside the source "
                    "domain; filled with default_value = {}",
                    n_outside, nt, rOptions.mDefaultValue);
            }

            for (std::size_t i = 0; i < point_names.size(); ++i) {
                const NDArray& src = smesh.PointData(point_names[i]);
                const std::size_t comps = data_num_components(src);
                NDArray arr = NDArray::Uninit(DType::Float64, interp_out_shape(src, nt));
                double* o = arr.As<double>();
                parallel_for(nt, [&](std::size_t t) {
                    if (simplex_of[t] >= 0) {
                        const std::int64_t* nodes =
                            &sconn[static_cast<std::size_t>(simplex_of[t]) * nv];
                        for (std::size_t c = 0; c < comps; ++c) {
                            // Corner-by-corner in connectivity order (fixed —
                            // the numpy twin sums in the same order).
                            double acc = wts[t][0] *
                                         detail::read_double(
                                             src, static_cast<std::size_t>(nodes[0]) * comps + c);
                            for (std::size_t k = 1; k < nv; ++k)
                                acc += wts[t][k] *
                                       detail::read_double(
                                           src, static_cast<std::size_t>(nodes[k]) * comps + c);
                            o[t * comps + c] = acc;
                        }
                    } else if (rOptions.mExtrapolate) {
                        for (std::size_t c = 0; c < comps; ++c)
                            o[t * comps + c] = detail::read_double(
                                src, static_cast<std::size_t>(near_of[t]) * comps + c);
                    } else {
                        for (std::size_t c = 0; c < comps; ++c)
                            o[t * comps + c] = rOptions.mDefaultValue;
                    }
                });
                out.AddPointData(point_out[i], std::move(arr));
            }
        }
    }

    // --- cell_data: always nearest source-cell centroid ---------------------
    if (!cell_names.empty()) {
        const std::size_t src_total = interp_total_cells(rSource);
        if (src_total == 0)
            throw std::invalid_argument(
                "meshio++: interpolate: cell_data was requested but the source mesh has no cells");
        if (rTarget.NumCellBlocks() == 0) {
            log::warn(
                "meshio++: interpolate: cell_data requested but the target has no cell blocks; "
                "skipping {} array(s)",
                cell_names.size());
        } else {
            const std::size_t tgt_total = interp_total_cells(rTarget);
            const std::vector<double> scent = interp_centroids(rSource, src_total);
            const std::vector<double> tcent = interp_centroids(rTarget, tgt_total);
            auto ccoord = [&](std::int64_t i) {
                const std::size_t j = static_cast<std::size_t>(i);
                return detail::Vec3{scent[3 * j], scent[3 * j + 1], scent[3 * j + 2]};
            };
            double lo[3], hi[3];
            interp_bbox(
                src_total, [&](std::size_t i) { return ccoord(static_cast<std::int64_t>(i)); }, lo,
                hi);
            detail::SpatialGrid grid(interp_cell_size(lo, hi, src_total));
            interp_fill_grid(grid, src_total,
                             [&](std::size_t i) { return ccoord(static_cast<std::int64_t>(i)); });

            std::vector<std::int64_t> nearest_cell(tgt_total);
            parallel_for(tgt_total, [&](std::size_t t) {
                const detail::Vec3 p{tcent[3 * t], tcent[3 * t + 1], tcent[3 * t + 2]};
                nearest_cell[t] = interp_nearest(grid, ccoord, p);
            });

            // Source block offsets in the global cell numbering.
            std::vector<std::size_t> soff;
            soff.reserve(rSource.NumCellBlocks() + 1);
            soff.push_back(0);
            for (const auto cb : rSource.CellRange())
                soff.push_back(soff.back() + cb.NumCells());

            for (std::size_t i = 0; i < cell_names.size(); ++i) {
                const std::string& name = cell_names[i];
                const std::size_t nblocks = rSource.CellDataNumBlocks(name);
                const std::size_t comps = data_num_components(rSource.CellData(name, 0));
                DType dt = rSource.CellData(name, 0).Dtype();
                for (std::size_t b = 1; b < nblocks; ++b)
                    if (rSource.CellData(name, b).Dtype() != dt)
                        dt = DType::Float64;  // mixed dtypes widen

                std::vector<NDArray> blocks_out;
                blocks_out.reserve(rTarget.NumCellBlocks());
                std::size_t tbase = 0;
                for (const auto cb : rTarget.CellRange()) {
                    const std::size_t ncells = cb.NumCells();
                    NDArray arr =
                        NDArray::Uninit(dt, interp_out_shape(rSource.CellData(name, 0), ncells));
                    parallel_for_bw(ncells, [&, tbase](std::size_t c) {
                        const std::size_t g = static_cast<std::size_t>(nearest_cell[tbase + c]);
                        std::size_t sb = 0;
                        while (g >= soff[sb + 1])
                            ++sb;
                        const std::size_t sl = g - soff[sb];
                        const NDArray& src = rSource.CellData(name, sb);
                        if (src.Dtype() == dt) {
                            const std::size_t rb = comps * dtype_size(dt);
                            std::memcpy(arr.Data() + c * rb, src.Data() + sl * rb, rb);
                        } else {
                            for (std::size_t k = 0; k < comps; ++k)
                                interp_store_double(arr, c * comps + k,
                                                    detail::read_double(src, sl * comps + k));
                        }
                    });
                    blocks_out.push_back(std::move(arr));
                    tbase += ncells;
                }
                out.AddCellData(cell_out[i], std::move(blocks_out));
            }
        }
    }

    return out;
}

}  // namespace meshioplusplus
