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
// One-pass mesh cleanup: weld coincident points (spatial-hash bucket grid, never
// O(N^2)), drop degenerate/duplicate cells, remove orphaned points, remapping
// connectivity and data. Built entirely through the uniform mesh API so it
// compiles under every mesh backend. See operations/clean.hpp for the contract.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Project includes
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/polyhedron.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

// Row-preserving gather: out[j] = src[idx[j]] (whole trailing dims), dtype kept.
NDArray clean_gather_rows(const NDArray& rSrc, const std::vector<std::int64_t>& rIdx) {
    const std::vector<std::size_t>& shp = rSrc.Shape();
    std::size_t cols = 1;
    for (std::size_t d = 1; d < shp.size(); ++d)
        cols *= shp[d];
    std::vector<std::size_t> out_shape = shp;
    if (out_shape.empty())
        out_shape = {rIdx.size()};
    else
        out_shape[0] = rIdx.size();
    NDArray out = NDArray::Uninit(rSrc.Dtype(), out_shape);
    const std::size_t rb = cols * dtype_size(rSrc.Dtype());
    if (rb == 0)
        return out;
    const std::byte* s = rSrc.Data();
    std::byte* d = out.Data();
    parallel_for_bw(rIdx.size(), [&](std::size_t j) {
        std::memcpy(d + j * rb, s + static_cast<std::size_t>(rIdx[j]) * rb, rb);
    });
    return out;
}

std::int64_t clean_quantize(double v, double atol) {
    return static_cast<std::int64_t>(std::floor(v / atol));
}

struct CleanKey {
    std::int64_t x, y, z;
    bool operator==(const CleanKey& r) const { return x == r.x && y == r.y && z == r.z; }
};

struct CleanKeyHash {
    std::size_t operator()(const CleanKey& k) const {
        std::uint64_t h = 1469598103934665603ULL;
        for (std::int64_t v : {k.x, k.y, k.z}) {
            h ^= static_cast<std::uint64_t>(v);
            h *= 1099511628211ULL;
        }
        return static_cast<std::size_t>(h);
    }
};

// Spatial-hash weld (keep-first, deterministic). Fills rWeldRep[g] in [0, W) and
// rRepSource[r] = the first global point that created representative r; returns W.
std::int64_t clean_build_weld_map(const NDArray& rPts, std::size_t n, std::size_t dim, double atol,
                                  std::vector<std::int64_t>& rWeldRep,
                                  std::vector<std::int64_t>& rRepSource) {
    rWeldRep.assign(n, 0);
    rRepSource.clear();
    const std::size_t ddim = std::min<std::size_t>(dim, 3);
    const double atol2 = atol * atol;
    std::unordered_map<CleanKey, std::vector<std::int64_t>, CleanKeyHash> grid;
    for (std::size_t g = 0; g < n; ++g) {
        double p[3] = {0, 0, 0};
        for (std::size_t d = 0; d < ddim; ++d)
            p[d] = detail::read_double(rPts, g * dim + d);
        const CleanKey k{clean_quantize(p[0], atol), clean_quantize(p[1], atol),
                         clean_quantize(p[2], atol)};
        std::int64_t found = -1;
        for (std::int64_t dz = -1; dz <= 1 && found < 0; ++dz)
            for (std::int64_t dy = -1; dy <= 1 && found < 0; ++dy)
                for (std::int64_t dx = -1; dx <= 1 && found < 0; ++dx) {
                    auto it = grid.find({k.x + dx, k.y + dy, k.z + dz});
                    if (it == grid.end())
                        continue;
                    for (std::int64_t rep : it->second) {
                        double q[3] = {0, 0, 0};
                        const std::size_t rg = static_cast<std::size_t>(rRepSource[rep]);
                        for (std::size_t d = 0; d < ddim; ++d)
                            q[d] = detail::read_double(rPts, rg * dim + d);
                        double d2 = 0;
                        for (std::size_t d = 0; d < ddim; ++d)
                            d2 += (p[d] - q[d]) * (p[d] - q[d]);
                        if (d2 <= atol2) {
                            found = rep;
                            break;
                        }
                    }
                }
        if (found >= 0) {
            rWeldRep[g] = found;
        } else {
            const std::int64_t r = static_cast<std::int64_t>(rRepSource.size());
            rWeldRep[g] = r;
            rRepSource.push_back(static_cast<std::int64_t>(g));
            grid[k].push_back(r);
        }
    }
    return static_cast<std::int64_t>(rRepSource.size());
}

}  // namespace

CleanResult clean(const Mesh& rMesh, const CleanOptions& rOpts) {
    const double eps = 1e-14;
    CleanResult res;

    const NDArray& points = rMesh.Points();
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    // --- weld (or identity) -> welded representatives -----------------------
    std::vector<std::int64_t> weld_rep;    // global point -> rep in [0, W)
    std::vector<std::int64_t> rep_source;  // rep -> first contributing global point
    std::int64_t W = 0;
    if (rOpts.weld && n > 0 && dim > 0 && rOpts.atol > 0) {
        W = clean_build_weld_map(points, n, dim, rOpts.atol, weld_rep, rep_source);
    } else {
        weld_rep.resize(n);
        rep_source.resize(n);
        for (std::size_t g = 0; g < n; ++g) {
            weld_rep[g] = static_cast<std::int64_t>(g);
            rep_source[g] = static_cast<std::int64_t>(g);
        }
        W = static_cast<std::int64_t>(n);
    }
    res.mPointsWelded = static_cast<std::int64_t>(n) - W;

    // --- per-block cell filtering (rep-space connectivity) ------------------
    struct BlockOut {
        std::string type;
        int kind;  // 0 rect, 1 polygon, 2 polyhedron
        std::size_t npc = 0;
        std::vector<std::int64_t> rect_conn;                        // kind 0, flat (kept*npc)
        std::vector<std::vector<std::int64_t>> poly_rows;           // kind 1
        std::vector<std::vector<std::vector<std::int64_t>>> polyh;  // kind 2
        std::vector<std::int64_t> kept_cells;                       // old local indices kept
    };
    std::vector<BlockOut> blocks;
    std::vector<char> rep_used(static_cast<std::size_t>(W), 0);

    for (const auto cb : rMesh.CellRange()) {
        BlockOut bo;
        bo.type = std::string(cb.Type());
        const CellType ct = cell_type_from_name(bo.type);
        const std::size_t nc = cb.NumCells();

        if (cb.IsPolyhedron()) {
            bo.kind = 2;
            // Degenerate/duplicate detection for polyhedra (v9.16.1). Before
            // that these branches kept every cell unconditionally, so a welded
            // polyhedron that had collapsed to nothing survived `clean` while
            // the equivalent hexahedron did not.
            //
            // "Degenerate" here is: a face that lost corners to the weld and is
            // no longer a polygon, a face set that is no longer a closed
            // orientable surface, or a volume that is negligible next to the
            // cell's own size. "Duplicate" keys on the SET OF FACES (each face
            // a sorted node set, the whole cell then sorted), so two cells that
            // list the same faces in a different order or with different
            // per-face rotations still collide -- a plain sorted node list
            // could not tell a cube from a differently-connected solid on the
            // same eight nodes.
            std::unordered_set<std::string> seen_poly;
            for (std::size_t c = 0; c < nc; ++c) {
                std::vector<std::vector<std::int64_t>> cell(cb.NumFaces(c));
                bool degenerate = false;
                for (std::size_t f = 0; f < cb.NumFaces(c); ++f) {
                    auto face = cb.Face(c, f);
                    cell[f].reserve(face.second);
                    for (std::size_t k = 0; k < face.second; ++k)
                        cell[f].push_back(weld_rep[static_cast<std::size_t>(face.first[k])]);
                    if (rOpts.drop_degenerate) {
                        std::vector<std::int64_t> u(cell[f]);
                        std::sort(u.begin(), u.end());
                        if (static_cast<std::size_t>(std::unique(u.begin(), u.end()) - u.begin()) <
                            3)
                            degenerate = true;  // the face collapsed below a triangle
                    }
                }
                if (rOpts.drop_degenerate && !degenerate) {
                    // Measure it through the shared kernel, on the WELDED nodes.
                    detail::CellRings rings;
                    std::vector<Vec3> coords;
                    rings.mFaceStart.push_back(0);
                    for (const std::vector<std::int64_t>& face : cell) {
                        for (std::int64_t id : face) {
                            std::uint32_t local = 0;
                            while (local < rings.mNodes.size() && rings.mNodes[local] != id)
                                ++local;
                            if (local == rings.mNodes.size()) {
                                rings.mNodes.push_back(id);
                                coords.push_back(detail::read_point(
                                    points, dim, rep_source[static_cast<std::size_t>(id)]));
                            }
                            rings.mFaceNodes.push_back(local);
                        }
                        rings.mFaceStart.push_back(
                            static_cast<std::uint32_t>(rings.mFaceNodes.size()));
                    }
                    if (detail::orient_rings(rings, coords.data()) ==
                        detail::RingOrientation::Unorientable) {
                        degenerate = true;
                    } else {
                        const detail::PolyMeasure pm = detail::poly_measure(rings, coords.data());
                        const double scale = pm.mSurfaceArea * std::sqrt(pm.mSurfaceArea);
                        if (!(std::abs(pm.mVolume) > eps * scale))
                            degenerate = true;
                    }
                }
                if (degenerate) {
                    ++res.mCellsDroppedDegenerate;
                    continue;
                }
                if (rOpts.drop_duplicate_cells) {
                    std::vector<std::string> face_keys;
                    face_keys.reserve(cell.size());
                    for (const std::vector<std::int64_t>& face : cell) {
                        std::vector<std::int64_t> sorted(face);
                        std::sort(sorted.begin(), sorted.end());
                        sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
                        std::string k;
                        for (std::int64_t v : sorted) {
                            k.append(reinterpret_cast<const char*>(&v), sizeof(v));
                            k.push_back(',');
                        }
                        face_keys.push_back(std::move(k));
                    }
                    std::sort(face_keys.begin(), face_keys.end());
                    std::string key;
                    for (const std::string& k : face_keys) {
                        key += k;
                        key.push_back(';');
                    }
                    if (!seen_poly.insert(std::move(key)).second) {
                        ++res.mCellsDroppedDuplicate;
                        continue;
                    }
                }
                for (const std::vector<std::int64_t>& face : cell)
                    for (std::int64_t r : face)
                        rep_used[static_cast<std::size_t>(r)] = 1;
                bo.polyh.push_back(std::move(cell));
                bo.kept_cells.push_back(static_cast<std::int64_t>(c));
            }
        } else if (cb.IsRagged()) {
            bo.kind = 1;
            std::unordered_set<std::string> seen_rows;
            for (std::size_t c = 0; c < nc; ++c) {
                std::vector<std::int64_t> row(cb.RowSize(c));
                for (std::size_t k = 0; k < cb.RowSize(c); ++k)
                    row[k] = weld_rep[static_cast<std::size_t>(cb.Row(c)[k])];
                if (rOpts.drop_degenerate) {
                    std::vector<std::int64_t> u(row);
                    std::sort(u.begin(), u.end());
                    if (static_cast<std::size_t>(std::unique(u.begin(), u.end()) - u.begin()) < 3) {
                        ++res.mCellsDroppedDegenerate;
                        continue;  // fewer than three distinct nodes is not a polygon
                    }
                    std::vector<Vec3> coords(row.size());
                    for (std::size_t k = 0; k < row.size(); ++k)
                        coords[k] = detail::read_point(
                            points, dim, rep_source[static_cast<std::size_t>(row[k])]);
                    if (!(detail::polygon_area(coords.data(), coords.size()) > eps)) {
                        ++res.mCellsDroppedDegenerate;
                        continue;
                    }
                }
                if (rOpts.drop_duplicate_cells) {
                    std::vector<std::int64_t> sorted(row);
                    std::sort(sorted.begin(), sorted.end());
                    std::string key(reinterpret_cast<const char*>(sorted.data()),
                                    sorted.size() * sizeof(std::int64_t));
                    if (!seen_rows.insert(std::move(key)).second) {
                        ++res.mCellsDroppedDuplicate;
                        continue;
                    }
                }
                for (std::int64_t r : row)
                    rep_used[static_cast<std::size_t>(r)] = 1;
                bo.poly_rows.push_back(std::move(row));
                bo.kept_cells.push_back(static_cast<std::int64_t>(c));
            }
        } else {
            bo.kind = 0;
            const NDArray& conn = cb.Conn();
            const std::size_t npc = cb.NodesPerCell();
            bo.npc = npc;
            const int corner_count = detail::cell_corner_count(ct);
            const int cdim = cell_type_dimension(ct);
            std::unordered_set<std::string> seen;
            std::vector<std::int64_t> row(npc);
            for (std::size_t c = 0; c < nc; ++c) {
                for (std::size_t k = 0; k < npc; ++k)
                    row[k] =
                        weld_rep[static_cast<std::size_t>(detail::read_int(conn, c * npc + k))];

                // degenerate: repeated corner node, or near-zero measure.
                bool degenerate = false;
                if (rOpts.drop_degenerate) {
                    const int cc = corner_count > 0 ? corner_count : static_cast<int>(npc);
                    for (int i = 0; i < cc && !degenerate; ++i)
                        for (int j = i + 1; j < cc; ++j)
                            if (row[i] == row[j]) {
                                degenerate = true;
                                break;
                            }
                    if (!degenerate && corner_count > 0) {
                        std::vector<Vec3> coords(corner_count);
                        for (int i = 0; i < corner_count; ++i)
                            coords[i] = detail::read_point(
                                points, dim, rep_source[static_cast<std::size_t>(row[i])]);
                        double measure = std::nan("");
                        if (cdim == 2)
                            measure = detail::polygon_area(coords.data(),
                                                           static_cast<std::size_t>(corner_count));
                        else if (cdim == 3)
                            measure = std::abs(detail::cell_volume_from_corners(coords.data(), ct));
                        if (!std::isnan(measure) && measure < eps)
                            degenerate = true;
                    }
                }
                if (degenerate) {
                    ++res.mCellsDroppedDegenerate;
                    continue;
                }

                // exact duplicate: identical sorted connectivity, keep-first.
                if (rOpts.drop_duplicate_cells) {
                    std::vector<std::int64_t> sorted(row);
                    std::sort(sorted.begin(), sorted.end());
                    std::string key(reinterpret_cast<const char*>(sorted.data()),
                                    sorted.size() * sizeof(std::int64_t));
                    if (!seen.insert(std::move(key)).second) {
                        ++res.mCellsDroppedDuplicate;
                        continue;
                    }
                }

                for (std::size_t k = 0; k < npc; ++k) {
                    bo.rect_conn.push_back(row[k]);
                    rep_used[static_cast<std::size_t>(row[k])] = 1;
                }
                bo.kept_cells.push_back(static_cast<std::int64_t>(c));
            }
        }
        blocks.push_back(std::move(bo));
    }

    // --- point compaction (orphan removal) ----------------------------------
    std::vector<std::int64_t> new_point(static_cast<std::size_t>(W), -1);
    std::vector<std::int64_t> kept_reps;  // reps surviving, ascending
    for (std::int64_t r = 0; r < W; ++r) {
        if (!rOpts.remove_orphans || rep_used[static_cast<std::size_t>(r)]) {
            new_point[static_cast<std::size_t>(r)] = static_cast<std::int64_t>(kept_reps.size());
            kept_reps.push_back(r);
        }
    }
    const std::int64_t P = static_cast<std::int64_t>(kept_reps.size());
    res.mPointsRemovedOrphan = W - P;

    // point gather indices (output row -> source global point, keep-first).
    std::vector<std::int64_t> point_src(kept_reps.size());
    for (std::size_t j = 0; j < kept_reps.size(); ++j)
        point_src[j] = rep_source[static_cast<std::size_t>(kept_reps[j])];

    // --- emit points + point_data ------------------------------------------
    Mesh& out = res.mMesh;
    {
        NDArray newpts = NDArray::Uninit(points.Dtype(), {static_cast<std::size_t>(P), dim});
        if (P > 0 && dim > 0) {
            const std::size_t rb = dim * dtype_size(points.Dtype());
            const std::byte* s = points.Data();
            std::byte* d = newpts.Data();
            parallel_for_bw(point_src.size(), [&](std::size_t j) {
                std::memcpy(d + j * rb, s + static_cast<std::size_t>(point_src[j]) * rb, rb);
            });
        }
        out.AssignPoints(std::move(newpts));
    }
    for (const std::string& name : rMesh.PointDataNames()) {
        const NDArray& a = rMesh.PointData(name);
        const std::size_t rows = a.Shape().empty() ? 0 : a.Shape()[0];
        if (rows == n)
            out.AddPointData(name, clean_gather_rows(a, point_src));
        else {
            NDArray c = a;
            c.MakeOwned();
            out.AddPointData(name, std::move(c));
        }
    }

    // --- emit cells (remap rep-space -> final point index) ------------------
    res.mCellMaps.clear();
    std::size_t bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        BlockOut& bo = blocks[bi];
        // cell map: old local -> new local (or -1).
        NDArray cmap = NDArray::Uninit(DType::Int64, {cb.NumCells()});
        std::int64_t* cm = cmap.As<std::int64_t>();
        parallel_for_bw(cb.NumCells(), [&](std::size_t c) { cm[c] = -1; });
        parallel_for_bw(bo.kept_cells.size(), [&](std::size_t j) {
            cm[static_cast<std::size_t>(bo.kept_cells[j])] = static_cast<std::int64_t>(j);
        });
        res.mCellMaps.push_back(std::move(cmap));

        if (bo.kind == 0) {
            const std::size_t kept = bo.kept_cells.size();
            NDArray conn = NDArray::Uninit(DType::Int64, {kept, bo.npc});
            std::int64_t* cd = conn.As<std::int64_t>();
            parallel_for_bw(bo.rect_conn.size(), [&](std::size_t i) {
                cd[i] = new_point[static_cast<std::size_t>(bo.rect_conn[i])];
            });
            out.AddCellBlock(bo.type, std::move(conn));
        } else if (bo.kind == 1) {
            parallel_for_bw(bo.poly_rows.size(), [&](std::size_t r) {
                for (std::int64_t& v : bo.poly_rows[r])
                    v = new_point[static_cast<std::size_t>(v)];
            });
            out.AddPolygonBlock(bo.type, std::move(bo.poly_rows));
        } else {
            parallel_for_bw(bo.polyh.size(), [&](std::size_t ci) {
                for (auto& face : bo.polyh[ci])
                    for (std::int64_t& v : face)
                        v = new_point[static_cast<std::size_t>(v)];
            });
            out.AddPolyhedronBlock(bo.type, std::move(bo.polyh));
        }
        ++bi;
    }

    // --- emit cell_data (gather surviving cells per block) ------------------
    for (const std::string& name : rMesh.CellDataNames()) {
        std::vector<NDArray> outblocks;
        std::size_t b = 0;
        for (const auto cb : rMesh.CellRange()) {
            const NDArray& src = rMesh.CellData(name, b);
            const std::size_t rows = src.Shape().empty() ? 0 : src.Shape()[0];
            if (rows == cb.NumCells())
                outblocks.push_back(clean_gather_rows(src, blocks[b].kept_cells));
            else {
                NDArray c = src;
                c.MakeOwned();
                outblocks.push_back(std::move(c));
            }
            ++b;
        }
        out.AddCellData(name, std::move(outblocks));
    }
    for (const std::string& name : rMesh.FieldDataNames()) {
        NDArray c = rMesh.FieldData(name);
        c.MakeOwned();
        out.AddFieldData(name, std::move(c));
    }

    // --- point map (input global -> output point, or -1) --------------------
    res.mPointMap = NDArray::Uninit(DType::Int64, {n});
    std::int64_t* pm = res.mPointMap.As<std::int64_t>();
    parallel_for_bw(
        n, [&](std::size_t g) { pm[g] = new_point[static_cast<std::size_t>(weld_rep[g])]; });

    // --- named regions ------------------------------------------------------
    // Block structure is 1:1 and cell types are unchanged, so this is a Direct
    // map and side-set facets survive along with their cells.
    {
        detail::RegionRemap rmap;
        rmap.pPointMap = &res.mPointMap;
        rmap.mCellMapKind = detail::CellMapKind::Direct;
        rmap.pCellMaps = &res.mCellMaps;
        rmap.mOpName = "clean";
        detail::remap_regions(rMesh, out, rmap);
    }

    return res;
}

}  // namespace meshioplusplus
