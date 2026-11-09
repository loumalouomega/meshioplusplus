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
// Planar cross-section of a mesh: the signed distance to the plane, handed to
// the shared marching-tetrahedra cutter. Everything below the distance is in
// `detail/marching.hpp`, hoisted out of this file in v7.15.0 (output unchanged,
// byte for byte) so `operations/isosurface.cpp` can cut on a field instead. See
// operations/slice.hpp for the contract, in particular the degeneracy/winding/
// determinism rules the pure-numpy twin (_slice.py) replicates expression for
// expression.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/marching.hpp"
#include "meshioplusplus/detail/subset.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/parallel.hpp"
#include "../detail/polyline_chain.hpp"

namespace meshioplusplus {

namespace {

NDArray slc_int64(const std::vector<std::int64_t>& rValues) {
    NDArray out = NDArray::Uninit(DType::Int64, {rValues.size()});
    std::copy(rValues.begin(), rValues.end(), out.As<std::int64_t>());
    return out;
}

// The line section @p rSection walked into chains (slice.hpp's contract).
Mesh slc_chain(const Mesh& rSection, double MinLength) {
    if (rSection.NumCellBlocks() == 0)
        return detail::clone_mesh(rSection);
    if (rSection.NumCellBlocks() != 1 || rSection.Cells(0).Type() != "line")
        throw std::invalid_argument(
            "meshio++: slice: chain needs a line section; slice a surface mesh (a volume's "
            "section is a surface)");
    const auto cells = rSection.Cells(0);
    const std::size_t nseg = cells.NumCells();
    const NDArray& conn = cells.Conn();
    const std::size_t npts = rSection.NumPoints();
    const std::size_t dim = rSection.PointDim();
    // Points at exactly one position (a plane through a vertex cuts each of
    // its edges there) are one point of the chain: the lowest id stands in.
    std::vector<std::int64_t> same(npts);
    {
        std::map<std::array<double, 3>, std::int64_t> first;
        for (std::size_t i = 0; i < npts; ++i) {
            std::array<double, 3> key{0.0, 0.0, 0.0};
            for (std::size_t a = 0; a < dim && a < 3; ++a)
                key[a] = detail::read_double(rSection.Points(), i * dim + a);
            same[i] = first.emplace(key, static_cast<std::int64_t>(i)).first->second;
        }
    }
    std::vector<std::array<std::int64_t, 2>> segs(nseg);
    for (std::size_t s = 0; s < nseg; ++s)
        segs[s] = {same[static_cast<std::size_t>(detail::read_int(conn, 2 * s))],
                   same[static_cast<std::size_t>(detail::read_int(conn, 2 * s + 1))]};
    const std::vector<detail::PolylineChain> chains = detail::polyline_chains(segs, npts);
    std::vector<std::int64_t> new_of(npts, -1), old_points, seg_order, conn_out, chain_id, closed;
    std::int64_t kept = 0;
    for (const detail::PolylineChain& r_c : chains) {
        double length = 0.0;
        for (std::size_t k = 0; k + 1 < r_c.mPoints.size(); ++k) {
            double d2 = 0.0;
            for (std::size_t a = 0; a < dim; ++a) {
                const double d =
                    detail::read_double(rSection.Points(),
                                        static_cast<std::size_t>(r_c.mPoints[k + 1]) * dim + a) -
                    detail::read_double(rSection.Points(),
                                        static_cast<std::size_t>(r_c.mPoints[k]) * dim + a);
                d2 += d * d;
            }
            length += std::sqrt(d2);
        }
        if (length < MinLength)
            continue;
        for (const std::int64_t p : r_c.mPoints)
            if (new_of[static_cast<std::size_t>(p)] < 0) {
                new_of[static_cast<std::size_t>(p)] = static_cast<std::int64_t>(old_points.size());
                old_points.push_back(p);
            }
        for (std::size_t k = 0; k < r_c.mSegments.size(); ++k) {
            seg_order.push_back(r_c.mSegments[k]);
            conn_out.push_back(new_of[static_cast<std::size_t>(r_c.mPoints[k])]);
            conn_out.push_back(new_of[static_cast<std::size_t>(r_c.mPoints[k + 1])]);
            chain_id.push_back(kept);
            closed.push_back(r_c.mClosed ? 1 : 0);
        }
        ++kept;
    }
    Mesh out;
    out.AssignPoints(detail::subset_gather_rows(rSection.Points(), old_points));
    for (const std::string& r_name : rSection.PointDataNames())
        out.AddPointData(r_name,
                         detail::subset_gather_rows(rSection.PointData(r_name), old_points));
    NDArray c = NDArray::Uninit(DType::Int64, {seg_order.size(), std::size_t{2}});
    std::copy(conn_out.begin(), conn_out.end(), c.As<std::int64_t>());
    out.AddCellBlock("line", std::move(c));
    for (const std::string& r_name : rSection.CellDataNames()) {
        std::vector<NDArray> blocks;
        blocks.push_back(detail::subset_gather_rows(rSection.CellData(r_name, 0), seg_order));
        out.AddCellData(r_name, std::move(blocks));
    }
    std::vector<NDArray> chain_blocks, closed_blocks;
    chain_blocks.push_back(slc_int64(chain_id));
    closed_blocks.push_back(slc_int64(closed));
    out.AddCellData(kSectionChainName, std::move(chain_blocks));
    out.AddCellData(kSectionClosedName, std::move(closed_blocks));
    return out;
}

}  // namespace

Mesh slice(const Mesh& rMesh, const SliceOptions& rOptions) {
    if (!(rOptions.mMinLength >= 0.0) || !std::isfinite(rOptions.mMinLength))
        throw std::invalid_argument("meshio++: slice: min_length must be finite and >= 0");
    using detail::Vec3;

    const Vec3 origin = {rOptions.mOrigin[0], rOptions.mOrigin[1], rOptions.mOrigin[2]};
    const Vec3 normal = {rOptions.mNormal[0], rOptions.mNormal[1], rOptions.mNormal[2]};
    if (detail::vec3_norm_sq(normal) <= 0.0)
        throw std::invalid_argument("meshio++: slice: normal must be non-zero");

    // Simplexify: every 3D cell becomes a tetra, every 2D cell a triangle.
    const detail::MarchingInput input = detail::marching_prepare(rMesh);
    const Mesh& simp = input.Simp();
    const std::size_t dim = simp.PointDim();

    // Signed distance of every simplexified node to the plane (the hot loop).
    const std::size_t nnodes_in = simp.NumPoints();
    std::vector<double> dist(nnodes_in, 0.0);
    {
        const NDArray& points = simp.Points();
        parallel_for(nnodes_in, [&](std::size_t g) {
            dist[g] = detail::vec3_dot(
                detail::vec3_sub(detail::read_point(points, dim, static_cast<std::int64_t>(g)),
                                 origin),
                normal);
        });
    }

    detail::MarchingOptions options;
    options.mProvenanceName = "slice:parent_cell";
    options.mRecordParentIds = rOptions.mRecordParentIds;
    // A plane has a global orientation, so every section face is wound toward
    // the user's +normal side (an isosurface has none and uses FieldGradient).
    options.mOrientation = detail::MarchingOrientation::FixedDirection;
    options.mDirection = rOptions.mNormal;

    detail::warn_regions_dropped(rMesh, "slice");
    if (rOptions.mChain)
        return slc_chain(detail::marching_cut(input, dist, options), rOptions.mMinLength);
    return detail::marching_cut(input, dist, options);
}

Mesh sections(const Mesh& rMesh, const SectionsOptions& rOptions) {
    const double len = std::sqrt(rOptions.mDirection[0] * rOptions.mDirection[0] +
                                 rOptions.mDirection[1] * rOptions.mDirection[1] +
                                 rOptions.mDirection[2] * rOptions.mDirection[2]);
    if (!(len > 0.0) || !std::isfinite(len))
        throw std::invalid_argument("meshio++: sections: direction must be non-zero");
    if (rOptions.mOffsets.empty())
        throw std::invalid_argument("meshio++: sections: give at least one offset");
    std::vector<Mesh> parts;
    std::int64_t chain_base = 0;
    for (const double off : rOptions.mOffsets) {
        if (!std::isfinite(off))
            throw std::invalid_argument("meshio++: sections: offsets must be finite");
        SliceOptions so;
        for (std::size_t a = 0; a < 3; ++a) {
            so.mNormal[a] = rOptions.mDirection[a] / len;
            so.mOrigin[a] = off * so.mNormal[a];
        }
        so.mRecordParentIds = rOptions.mRecordParentIds;
        so.mChain = true;
        so.mMinLength = rOptions.mMinLength;
        Mesh part = slice(rMesh, so);
        if (part.NumCellBlocks() == 0 || part.Cells(0).NumCells() == 0)
            continue;
        const std::size_t n = part.Cells(0).NumCells();
        std::vector<std::int64_t> chains(n);
        std::int64_t top = -1;
        for (std::size_t k = 0; k < n; ++k) {
            const std::int64_t c = detail::read_int(part.CellData(kSectionChainName, 0), k);
            chains[k] = c + chain_base;
            top = std::max(top, c);
        }
        chain_base += top + 1;
        NDArray offs = NDArray::Uninit(DType::Float64, {n});
        std::fill(offs.As<double>(), offs.As<double>() + n, off);
        std::vector<NDArray> cb, ob;
        cb.push_back(slc_int64(chains));
        ob.push_back(std::move(offs));
        part.AddCellData(kSectionChainName, std::move(cb));
        part.AddCellData(kSectionOffsetName, std::move(ob));
        parts.push_back(std::move(part));
    }
    if (parts.empty())
        return Mesh{};
    if (parts.size() == 1)
        return std::move(parts[0]);
    std::vector<const Mesh*> ptrs;
    for (const Mesh& r_p : parts)
        ptrs.push_back(&r_p);
    MergeOptions mo;
    mo.source_tag = false;
    return merge(ptrs, mo).mMesh;
}

}  // namespace meshioplusplus
