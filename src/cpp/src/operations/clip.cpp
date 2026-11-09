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
// Cut a surface by a plane, keep one side and cap it (roadmap 7.8.3). See
// operations/clip.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/clip.hpp"
#include "../detail/cdt2d.hpp"
#include "../detail/lsq_fit.hpp"
#include "../detail/polyline_chain.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/marching.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kClPrefix = "meshio++: clip: ";

std::uint64_t cl_edge_key(std::int64_t A, std::int64_t B) {
    const std::int64_t lo = std::min(A, B), hi = std::max(A, B);
    return (static_cast<std::uint64_t>(lo) << 32) ^ static_cast<std::uint64_t>(hi);
}

}  // namespace

ClipResult clip(const Mesh& rMesh, const ClipOptions& rOptions) {
    const double nn = std::sqrt(rOptions.mNormal[0] * rOptions.mNormal[0] +
                                rOptions.mNormal[1] * rOptions.mNormal[1] +
                                rOptions.mNormal[2] * rOptions.mNormal[2]);
    if (!(nn > 0.0) || !std::isfinite(nn))
        throw std::invalid_argument(std::string(kClPrefix) +
                                    "the normal must be a non-zero finite vector");
    for (const auto cb : rMesh.CellRange())
        if (cb.IsPolyhedron() ||
            cell_type_dimension(cell_type_from_name(std::string(cb.Type()))) == 3)
            throw std::invalid_argument(std::string(kClPrefix) + "cell block '" +
                                        std::string(cb.Type()) +
                                        "' is a volume; clip cuts surfaces (extract_skin first)");
    const detail::MarchingInput inp = detail::marching_prepare(rMesh);
    const Mesh& simp = inp.Simp();
    const std::size_t np = simp.NumPoints();
    const std::size_t dim = simp.PointDim();
    const NDArray& pts = simp.Points();
    auto coord = [&](std::int64_t V, std::size_t K) {
        return K < dim ? detail::read_double(pts, static_cast<std::size_t>(V) * dim + K) : 0.0;
    };
    // Slice's signed distance, d = 0 counting as above.
    std::vector<double> d(np);
    for (std::size_t i = 0; i < np; ++i) {
        double s = 0.0;
        for (std::size_t k = 0; k < 3; ++k)
            s += (coord(static_cast<std::int64_t>(i), k) - rOptions.mOrigin[k]) *
                 rOptions.mNormal[k];
        d[i] = s;
    }
    auto kept_node = [&](std::int64_t V) {
        const bool below = d[static_cast<std::size_t>(V)] < 0.0;
        return rOptions.mKeepBelow ? below : !below;
    };

    ClipResult result;
    // Output points: the kept input points, then the crossings (one per cut edge).
    std::vector<std::int64_t> out_of(np, -1);
    std::vector<std::array<double, 3>> xyz;
    struct ClCross {
        std::int64_t mLo, mHi;
        double mT;
    };
    std::vector<ClCross> crossings;    // per output point beyond the kept ones: its edge
    std::vector<std::int64_t> source;  // per output point: its input point, or -1
    std::unordered_map<std::uint64_t, std::int64_t> cross_of;
    auto keep_point = [&](std::int64_t V) {
        std::int64_t& r_o = out_of[static_cast<std::size_t>(V)];
        if (r_o < 0) {
            r_o = static_cast<std::int64_t>(xyz.size());
            xyz.push_back({coord(V, 0), coord(V, 1), coord(V, 2)});
            source.push_back(V);
            crossings.push_back({-1, -1, 0.0});
        }
        return r_o;
    };
    auto cross_point = [&](std::int64_t A, std::int64_t B) {
        const std::uint64_t key = cl_edge_key(A, B);
        const auto it = cross_of.find(key);
        if (it != cross_of.end())
            return it->second;
        const std::int64_t lo = std::min(A, B), hi = std::max(A, B);
        const double dlo = d[static_cast<std::size_t>(lo)], dhi = d[static_cast<std::size_t>(hi)];
        const double t = dlo / (dlo - dhi);
        std::array<double, 3> x{};
        for (std::size_t k = 0; k < 3; ++k)
            x[k] = coord(lo, k) + t * (coord(hi, k) - coord(lo, k));
        const auto id = static_cast<std::int64_t>(xyz.size());
        xyz.push_back(x);
        source.push_back(-1);
        crossings.push_back({lo, hi, t});
        cross_of.emplace(key, id);
        return id;
    };

    std::vector<std::int64_t> tris;                           // output triangles
    std::vector<std::pair<std::size_t, std::size_t>> parent;  // (simp block, cell) per triangle
    std::vector<std::array<std::int64_t, 2>> cut_segments;    // between crossings
    for (std::size_t b = 0; b < simp.NumCellBlocks(); ++b) {
        const auto cb = simp.Cells(b);
        if (cell_type_from_name(std::string(cb.Type())) != CellType::Triangle)
            continue;  // lines and vertices are not clipped
        const NDArray& conn = cb.Conn();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const std::array<std::int64_t, 3> v{detail::read_int(conn, 3 * c),
                                                detail::read_int(conn, 3 * c + 1),
                                                detail::read_int(conn, 3 * c + 2)};
            const int nk =
                (kept_node(v[0]) ? 1 : 0) + (kept_node(v[1]) ? 1 : 0) + (kept_node(v[2]) ? 1 : 0);
            if (nk == 0)
                continue;
            if (nk == 3) {
                for (const std::int64_t x : v)
                    tris.push_back(keep_point(x));
                parent.emplace_back(b, c);
                ++result.mNumKept;
                continue;
            }
            ++result.mNumSplit;
            // Rotate so the odd corner (the lone kept or the lone dropped one) is first.
            std::size_t r = 0;
            for (std::size_t k = 0; k < 3; ++k)
                if (kept_node(v[k]) == (nk == 1))
                    r = k;
            const std::int64_t a = v[r], bnode = v[(r + 1) % 3], cnode = v[(r + 2) % 3];
            const std::int64_t xab = cross_point(a, bnode), xac = cross_point(a, cnode);
            if (nk == 1) {
                tris.insert(tris.end(), {keep_point(a), xab, xac});
                parent.emplace_back(b, c);
                cut_segments.push_back({xab, xac});
            } else {
                const std::int64_t pb = keep_point(bnode), pc = keep_point(cnode);
                tris.insert(tris.end(), {xab, pb, pc});
                parent.emplace_back(b, c);
                tris.insert(tris.end(), {xab, pc, xac});
                parent.emplace_back(b, c);
                cut_segments.push_back({xab, xac});
            }
        }
    }
    const std::size_t num_body = tris.size() / 3;

    // The cap: the cut's closed loops, triangulated in the plane.
    if (rOptions.mCap && !cut_segments.empty()) {
        const std::vector<detail::PolylineChain> chains =
            detail::polyline_chains(cut_segments, xyz.size());
        std::array<double, 3> n{rOptions.mNormal[0] / nn, rOptions.mNormal[1] / nn,
                                rOptions.mNormal[2] / nn};
        std::array<double, 3> e1{}, e2{};
        detail::lsq_frame(n.data(), e1.data(), e2.data());
        std::vector<std::vector<detail::Cdt2d::Vec2>> loops;
        std::map<std::array<double, 2>, std::int64_t> id_of;
        for (const detail::PolylineChain& r_c : chains) {
            if (!r_c.mClosed) {
                ++result.mNumOpenChains;
                continue;
            }
            std::vector<detail::Cdt2d::Vec2> loop;
            for (std::size_t i = 0; i + 1 < r_c.mPoints.size(); ++i) {
                const auto& x = xyz[static_cast<std::size_t>(r_c.mPoints[i])];
                const detail::Cdt2d::Vec2 uv{x[0] * e1[0] + x[1] * e1[1] + x[2] * e1[2],
                                             x[0] * e2[0] + x[1] * e2[1] + x[2] * e2[2]};
                id_of.emplace(uv, r_c.mPoints[i]);
                loop.push_back(uv);
            }
            if (loop.size() >= 3) {
                loops.push_back(std::move(loop));
                ++result.mNumLoops;
            }
        }
        if (result.mNumOpenChains > 0)
            result.mWarnings.push_back(
                std::to_string(result.mNumOpenChains) +
                " open chain(s) of the cut are not capped (the surface is open there)");
        if (!loops.empty()) {
            const detail::Cdt2dPolygon poly = detail::cdt2d_triangulate_polygon(loops, 0.0, 0.0);
            if (poly.mStatus != detail::Cdt2dStatus::Ok) {
                result.mWarnings.push_back(std::string("the cap could not be triangulated (") +
                                           detail::cdt2d_status_name(poly.mStatus) +
                                           "): " + poly.mWhy);
            } else {
                std::vector<std::int64_t> map(poly.mPoints.size(), -1);
                bool ok = true;
                for (std::size_t i = 0; i < poly.mPoints.size(); ++i) {
                    const auto it = id_of.find(poly.mPoints[i]);
                    if (it == id_of.end())
                        ok = false;
                    else
                        map[i] = it->second;
                }
                if (!ok) {
                    result.mWarnings.push_back("the cap's triangulation added points; not capped");
                } else {
                    for (const auto& r_t : poly.mTriangles) {
                        const std::int64_t a = map[static_cast<std::size_t>(r_t[0])];
                        const std::int64_t b = map[static_cast<std::size_t>(r_t[1])];
                        const std::int64_t c = map[static_cast<std::size_t>(r_t[2])];
                        // Counter-clockwise about the normal faces +n: out of
                        // the kept side when keeping below.
                        if (rOptions.mKeepBelow)
                            tris.insert(tris.end(), {a, b, c});
                        else
                            tris.insert(tris.end(), {a, c, b});
                        ++result.mNumCapTriangles;
                    }
                }
            }
        }
    }

    // Assemble.
    Mesh& out = result.mMesh;
    NDArray p = NDArray::Uninit(DType::Float64, {xyz.size(), std::size_t{3}});
    for (std::size_t i = 0; i < xyz.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            p.As<double>()[3 * i + k] = xyz[i][k];
    out.AssignPoints(std::move(p));
    const std::size_t nt = tris.size() / 3;
    NDArray c = NDArray::Uninit(DType::Int64, {nt, std::size_t{3}});
    std::copy(tris.begin(), tris.end(), c.As<std::int64_t>());
    out.AddCellBlock("triangle", std::move(c));
    for (const std::string& r_name : simp.PointDataNames()) {
        const NDArray& src = simp.PointData(r_name);
        const std::size_t comps = src.Size() / std::max<std::size_t>(np, 1);
        std::vector<std::size_t> shape = src.Shape();
        shape[0] = xyz.size();
        NDArray dst = NDArray::Uninit(DType::Float64, shape);
        for (std::size_t i = 0; i < xyz.size(); ++i)
            for (std::size_t k = 0; k < comps; ++k) {
                double v;
                if (source[i] >= 0) {
                    v = detail::read_double(src, static_cast<std::size_t>(source[i]) * comps + k);
                } else {
                    const ClCross& r_x = crossings[i];
                    const double lo =
                        detail::read_double(src, static_cast<std::size_t>(r_x.mLo) * comps + k);
                    const double hi =
                        detail::read_double(src, static_cast<std::size_t>(r_x.mHi) * comps + k);
                    v = lo + r_x.mT * (hi - lo);
                }
                dst.As<double>()[i * comps + k] = v;
            }
        out.AddPointData(r_name, std::move(dst));
    }
    for (const std::string& r_name : simp.CellDataNames()) {
        if (r_name == "convert:parent_cell")
            continue;
        const NDArray& first = simp.CellData(r_name, parent.empty() ? 0 : parent[0].first);
        const std::size_t comps =
            first.Shape().empty() || first.Shape()[0] == 0 ? 1 : first.Size() / first.Shape()[0];
        const bool floating = first.Dtype() == DType::Float64 || first.Dtype() == DType::Float32;
        std::vector<std::size_t> shape = first.Shape();
        if (shape.empty())
            shape = {0};
        shape[0] = nt;
        NDArray dst = NDArray::Uninit(floating ? DType::Float64 : DType::Int64, shape);
        for (std::size_t t = 0; t < nt; ++t)
            for (std::size_t k = 0; k < comps; ++k) {
                if (t < num_body) {
                    const NDArray& src = simp.CellData(r_name, parent[t].first);
                    if (floating)
                        dst.As<double>()[t * comps + k] =
                            detail::read_double(src, parent[t].second * comps + k);
                    else
                        dst.As<std::int64_t>()[t * comps + k] =
                            detail::read_int(src, parent[t].second * comps + k);
                } else if (floating) {
                    dst.As<double>()[t * comps + k] = std::numeric_limits<double>::quiet_NaN();
                } else {
                    dst.As<std::int64_t>()[t * comps + k] = -1;
                }
            }
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(dst));
        out.AddCellData(r_name, std::move(blocks));
    }
    std::vector<std::int64_t> cap(nt, 0);
    for (std::size_t t = num_body; t < nt; ++t)
        cap[t] = 1;
    NDArray cap_arr = NDArray::Uninit(DType::Int64, {nt});
    std::copy(cap.begin(), cap.end(), cap_arr.As<std::int64_t>());
    std::vector<NDArray> cap_blocks;
    cap_blocks.push_back(std::move(cap_arr));
    out.AddCellData("clip:cap", std::move(cap_blocks));
    for (const std::string& r_name : rMesh.FieldDataNames())
        out.AddFieldData(r_name, rMesh.FieldData(r_name));
    return result;
}

}  // namespace meshioplusplus
