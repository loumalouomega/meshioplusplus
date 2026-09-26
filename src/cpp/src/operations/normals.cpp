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
// compute_normals. See operations/normals.hpp for the contract and
// detail/surface_normals.hpp for the corner grouping it is built on.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/surface_normals.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "../detail/typed_view.hpp"

namespace meshioplusplus {
namespace {

using detail::Vec3;

constexpr const char* kNrmPrefix = "meshio++: normals: ";

// Refuse what a normal is not defined on, naming the fix. build_triangle_soup
// refuses the same blocks, but with a message about distances.
void nrm_check_surface(const Mesh& rMesh, const std::string& rRegion) {
    for (const auto cb : rMesh.CellRange()) {
        const std::string type(cb.Type());
        const CellType ct = cell_type_from_name(type);
        if (cb.IsPolyhedron() || cell_type_dimension(ct) == 3)
            throw std::invalid_argument(std::string(kNrmPrefix) + "cell block '" + type +
                                        "' is a volume; normals are defined on a surface (run "
                                        "extract_surface first)");
        const bool polygon = type.rfind("polygon", 0) == 0;
        if (!polygon && cell_type_dimension(ct) == 2 && ct != CellType::Triangle &&
            ct != CellType::Quad)
            throw std::invalid_argument(std::string(kNrmPrefix) + "cell block '" + type +
                                        "' is a higher-order surface cell (run linearize first)");
    }
    if (!rRegion.empty() && rMesh.FindRegion(rRegion, RegionKind::Cell) == Mesh::npos) {
        std::string names;
        for (const std::string& n : rMesh.RegionNames())
            names += (names.empty() ? "" : ", ") + n;
        throw std::invalid_argument(std::string(kNrmPrefix) + "no cell region named '" + rRegion +
                                    "' (available: " + (names.empty() ? "none" : names) + ")");
    }
}

// The point data of the input, gathered onto the split layout: originals
// verbatim, each copy from its parent's row.
NDArray nrm_gather_rows(const NDArray& rArray, std::size_t N, std::size_t NOut,
                        const std::vector<std::int64_t>& rParentOfNew) {
    std::vector<std::size_t> shape = rArray.Shape();
    shape[0] = NOut;
    NDArray out = NDArray::Uninit(rArray.Dtype(), std::move(shape));
    std::memcpy(out.Data(), rArray.Data(), rArray.Nbytes());
    const std::size_t row_bytes = rArray.Nbytes() / N;
    for (std::size_t k = 0; k < rParentOfNew.size(); ++k)
        std::memcpy(out.Data() + (N + k) * row_bytes,
                    rArray.Data() + static_cast<std::size_t>(rParentOfNew[k]) * row_bytes,
                    row_bytes);
    return out;
}

}  // namespace

NormalsResult compute_normals(const Mesh& rMesh, const NormalsOptions& rOptions) {
    if (rOptions.mSplit && !(rOptions.mSplitAngle >= 0.0 && rOptions.mSplitAngle <= 180.0))
        throw std::invalid_argument(std::string(kNrmPrefix) +
                                    "the split angle must lie in [0, 180] degrees");
    nrm_check_surface(rMesh, rOptions.mRegion);

    const detail::TriangleSoup soup = detail::build_triangle_soup(rMesh, rOptions.mRegion);
    NormalsResult result;
    result.mQuality = detail::soup_quality(soup);

    const std::size_t n = rMesh.NumPoints();
    const std::size_t ntri = soup.NumTriangles();
    const detail::VertexNormalGroups groups = detail::vertex_normal_groups(
        soup, rOptions.mWeight, rOptions.mSplit ? rOptions.mSplitAngle : -1.0);
    result.mNumDegenerate = groups.mNumDegenerate;

    // A group with no direction (only degenerate triangles) cannot be given a
    // normal, so it is not worth a copy of its point: it folds into the point's
    // primary group, the first defined one.
    const std::size_t ngroups = groups.NumGroups();
    auto defined = [&](std::size_t g) {
        const Vec3& v = groups.mGroupNormal[g];
        return v[0] != 0.0 || v[1] != 0.0 || v[2] != 0.0;
    };
    std::vector<std::int64_t> primary(n, -1);
    for (std::size_t g = 0; g < ngroups; ++g) {
        std::int64_t& slot = primary[static_cast<std::size_t>(groups.mGroupPoint[g])];
        if (slot < 0 && defined(g))
            slot = static_cast<std::int64_t>(g);
    }
    for (std::size_t g = 0; g < ngroups; ++g) {
        std::int64_t& slot = primary[static_cast<std::size_t>(groups.mGroupPoint[g])];
        if (slot < 0)
            slot = static_cast<std::int64_t>(g);
    }

    // Output index of every group: its own point for the primary group (and for
    // an undefined one), a new appended point for any other defined group.
    std::vector<std::int64_t> group_point(ngroups);
    std::vector<std::int64_t> parent_of_new;
    std::vector<char> has_copy(n, 0);
    for (std::size_t g = 0; g < ngroups; ++g) {
        const std::int64_t p = groups.mGroupPoint[g];
        if (static_cast<std::int64_t>(g) == primary[static_cast<std::size_t>(p)] || !defined(g)) {
            group_point[g] = p;
        } else {
            group_point[g] = static_cast<std::int64_t>(n + parent_of_new.size());
            parent_of_new.push_back(p);
            has_copy[static_cast<std::size_t>(p)] = 1;
        }
    }
    const std::size_t n_added = parent_of_new.size();
    const std::size_t n_out = n + n_added;
    result.mNumAddedPoints = static_cast<std::int64_t>(n_added);
    result.mNumSplitPoints =
        static_cast<std::int64_t>(std::count(has_copy.begin(), has_copy.end(), char{1}));

    // The point normals.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    NDArray point_normals(DType::Float64, {n_out, std::size_t{3}});
    double* pn = point_normals.As<double>();
    parallel_for_bw(n_out * 3, [&](std::size_t i) { pn[i] = nan; });
    // Each written row belongs to exactly one group (its primary, or its own
    // appended copy), so the groups fill their rows independently.
    parallel_for(ngroups, [&](std::size_t g) {
        const std::int64_t p = groups.mGroupPoint[g];
        const bool is_primary =
            static_cast<std::int64_t>(g) == primary[static_cast<std::size_t>(p)];
        if (!is_primary && group_point[g] == p)
            return;  // an undefined group folded into the primary one
        if (!defined(g))
            return;  // stays NaN
        const std::size_t row = static_cast<std::size_t>(group_point[g]);
        for (std::size_t k = 0; k < 3; ++k)
            pn[row * 3 + k] = groups.mGroupNormal[g][k];
    });
    for (std::size_t p = 0; p < n; ++p) {
        if (primary[p] < 0)
            ++result.mNumIsolated;
        else if (!defined(static_cast<std::size_t>(primary[p])))
            ++result.mNumUndefined;
    }

    // Cell normals: the sum of each cell's fan crosses (Newell), normalised.
    std::vector<Vec3> cell_sum;
    std::vector<std::int64_t> first_tri;
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t ncells_total = static_cast<std::size_t>(detail::total_cells(bases));
    if (rOptions.mCellNormals || n_added > 0) {
        first_tri.assign(ncells_total, -1);
        cell_sum.assign(ncells_total, Vec3{0.0, 0.0, 0.0});
        const std::vector<Vec3> face = detail::soup_face_normals(soup);
        // A cell's fan triangles are consecutive in the soup, in fan order:
        // each run's first triangle sums the run, in the serial order.
        parallel_for(ntri, [&](std::size_t t) {
            const std::int64_t c = soup.mSourceCell[t];
            if (t > 0 && soup.mSourceCell[t - 1] == c)
                return;
            Vec3 acc = cell_sum[static_cast<std::size_t>(c)];
            for (std::size_t u = t; u < ntri && soup.mSourceCell[u] == c; ++u)
                acc = detail::vec3_add(acc, face[u]);
            first_tri[static_cast<std::size_t>(c)] = static_cast<std::int64_t>(t);
            cell_sum[static_cast<std::size_t>(c)] = acc;
        });
    }

    if (n_added == 0) {
        result.mMesh = detail::clone_mesh(rMesh);
    } else {
        // Rebuild with the corners of every selected cell pointed at their
        // group's point. Cells the soup did not take (lines, vertices, other
        // regions, degenerate rows) keep their original ids.
        Mesh out;
        const NDArray& pts = rMesh.Points();
        const std::size_t dim = rMesh.PointDim();
        NDArray new_pts = NDArray::Uninit(pts.Dtype(), {n_out, dim});
        std::memcpy(new_pts.Data(), pts.Data(), pts.Nbytes());
        const std::size_t row_bytes = pts.Nbytes() / n;
        for (std::size_t k = 0; k < n_added; ++k)
            std::memcpy(new_pts.Data() + (n + k) * row_bytes,
                        pts.Data() + static_cast<std::size_t>(parent_of_new[k]) * row_bytes,
                        row_bytes);
        out.AssignPoints(std::move(new_pts));

        std::size_t bi = 0;
        for (const auto cb : rMesh.CellRange()) {
            const std::int64_t base = bases[bi++];
            const std::size_t ncells = cb.NumCells();
            auto corner_point = [&](std::size_t t, std::size_t i) {
                return group_point[static_cast<std::size_t>(groups.mCornerGroup[t * 3 + i])];
            };
            // The rewritten ids of one cell, given its original ids.
            auto rewrite = [&](std::size_t c, std::vector<std::int64_t>& rIds) {
                const std::int64_t t0 = first_tri[static_cast<std::size_t>(base) + c];
                if (t0 < 0 || rIds.size() < 3)
                    return;
                const std::size_t t = static_cast<std::size_t>(t0);
                const std::size_t nv = rIds.size();
                rIds[0] = corner_point(t, 0);
                for (std::size_t k = 1; k + 1 < nv; ++k)
                    rIds[k] = corner_point(t + k - 1, 1);
                rIds[nv - 1] = corner_point(t + nv - 3, 2);
            };
            if (cb.IsRagged()) {
                // Straight into CSR: the rows keep their sizes, so the offsets
                // are the input's.
                std::vector<std::int64_t> offsets(ncells + 1, 0);
                for (std::size_t c = 0; c < ncells; ++c)
                    offsets[c + 1] = offsets[c] + static_cast<std::int64_t>(cb.RowSize(c));
                std::vector<std::int64_t> flat(static_cast<std::size_t>(offsets[ncells]));
                parallel_for(ncells, [&](std::size_t c) {
                    thread_local std::vector<std::int64_t> ids;
                    ids.assign(cb.Row(c), cb.Row(c) + cb.RowSize(c));
                    rewrite(c, ids);
                    std::copy(ids.begin(), ids.end(),
                              flat.begin() + static_cast<std::ptrdiff_t>(offsets[c]));
                });
                out.AddPolygonBlock(std::string(cb.Type()), std::move(flat), std::move(offsets));
            } else {
                const detail::Int64View conn(cb.Conn());
                const std::size_t npc = cb.NodesPerCell();
                NDArray new_conn = NDArray::Uninit(DType::Int64, {ncells, npc});
                std::int64_t* dst = new_conn.As<std::int64_t>();
                parallel_for(ncells, [&](std::size_t c) {
                    thread_local std::vector<std::int64_t> ids;
                    ids.assign(conn.Data() + c * npc, conn.Data() + (c + 1) * npc);
                    rewrite(c, ids);
                    std::copy(ids.begin(), ids.end(), dst + c * npc);
                });
                out.AddCellBlock(std::string(cb.Type()), std::move(new_conn));
            }
        }

        for (const std::string& name : rMesh.PointDataNames()) {
            const NDArray& a = rMesh.PointData(name);
            out.AddPointData(name, (detail::rows(a) == n && n > 0)
                                       ? nrm_gather_rows(a, n, n_out, parent_of_new)
                                       : detail::data_owned_copy(a));
        }
        for (const std::string& name : rMesh.CellDataNames()) {
            std::vector<NDArray> blocks;
            for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name); ++b)
                blocks.push_back(detail::data_owned_copy(rMesh.CellData(name, b)));
            out.AddCellData(name, std::move(blocks));
        }
        for (const std::string& name : rMesh.FieldDataNames())
            out.AddFieldData(name, detail::data_owned_copy(rMesh.FieldData(name)));

        // Regions: cell numbering is untouched, so everything rides through,
        // and a copy joins its source's Point regions.
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const meshioplusplus::Region& r = rMesh.Region(i);
            if (r.mKind != RegionKind::Point) {
                out.AddRegion(r);
                continue;
            }
            std::vector<std::int64_t> entries(r.mEntries.Size());
            for (std::size_t e = 0; e < entries.size(); ++e)
                entries[e] = detail::read_int(r.mEntries, e);
            const std::size_t original = entries.size();
            for (std::size_t k = 0; k < n_added; ++k)
                if (std::binary_search(entries.begin(), entries.begin() + original,
                                       parent_of_new[k]))
                    entries.push_back(static_cast<std::int64_t>(n + k));
            meshioplusplus::Region nr = r;
            nr.mEntries = NDArray::Uninit(DType::Int64, {entries.size()});
            std::memcpy(nr.mEntries.Data(), entries.data(), entries.size() * sizeof(std::int64_t));
            out.AddRegion(std::move(nr));
        }
        for (std::size_t i = 0; i < rMesh.NumPropertySets(); ++i)
            out.AddPropertySet(rMesh.GetPropertySet(i));
        result.mMesh = std::move(out);
    }

    if (rOptions.mPointNormals)
        result.mMesh.AddPointData(kNormalsName, std::move(point_normals));

    if (rOptions.mCellNormals) {
        std::vector<NDArray> blocks;
        for (const auto cb : rMesh.CellRange()) {
            const std::size_t bi = blocks.size();
            const std::size_t ncells = cb.NumCells();
            NDArray a(DType::Float64, {ncells, std::size_t{3}});
            double* d = a.As<double>();
            parallel_for(ncells, [&](std::size_t c) {
                const std::size_t global = static_cast<std::size_t>(bases[bi]) + c;
                const Vec3& s = cell_sum[global];
                const double len = detail::vec3_norm(s);
                for (std::size_t k = 0; k < 3; ++k)
                    d[c * 3 + k] = (first_tri[global] >= 0 && len > 0.0) ? s[k] * (1.0 / len) : nan;
            });
            blocks.push_back(std::move(a));
        }
        result.mMesh.AddCellData(kNormalsName, std::move(blocks));
    }

    if (rOptions.mRecordParentIds) {
        NDArray parent = NDArray::Uninit(DType::Int64, {n_out});
        std::int64_t* dst = parent.As<std::int64_t>();
        for (std::size_t i = 0; i < n; ++i)
            dst[i] = static_cast<std::int64_t>(i);
        for (std::size_t k = 0; k < n_added; ++k)
            dst[n + k] = parent_of_new[k];
        result.mMesh.AddPointData(kNormalsParentPointName, std::move(parent));
    }

    if (result.mQuality.mInconsistentPairs != 0 && !rOptions.mSplit)
        log::warn(
            "normals: {} edge pair(s) wind the same way, so the normals there average faces "
            "that disagree about which side is out; run repair(mesh, {{.mFixOrientation = true}}) "
            "first",
            result.mQuality.mInconsistentPairs);
    return result;
}

}  // namespace meshioplusplus
