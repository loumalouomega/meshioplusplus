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
// Mesh merge/weld: concatenate (append points, offset connectivity, merge cell
// blocks by type, concatenate data, tag provenance) with optional welding of
// coincident nodes via a spatial-hash bucket grid (never O(N^2)). Built entirely
// through the uniform mesh API so it compiles under every mesh backend. See
// operations/merge.hpp for the contract.

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/detail/spatial_hash.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// A deep copy of an NDArray that always owns its buffer (the source may be a
// view over foreign memory, e.g. numpy memory on the write path).
NDArray merge_owned_copy(const NDArray& rArr) {
    NDArray c = rArr;
    c.MakeOwned();
    return c;
}

// Store `v` at flat element index `flat` of `rOut`, casting to `rOut`'s dtype.
void merge_store_double(NDArray& rOut, std::size_t flat, double v) {
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

// Copy `nrows` rows (`ncols` wide) from `rSrc` (starting at `srcRow0`) into
// `rOut` (starting at `outRow0`), memcpy when the dtypes match, else casting
// element-wise. `rOut`/`rSrc` are treated as row-major (ncols columns).
void merge_copy_rows(NDArray& rOut, std::size_t outRow0, const NDArray& rSrc, std::size_t srcRow0,
                     std::size_t nrows, std::size_t ncols) {
    if (nrows == 0 || ncols == 0)
        return;
    if (rOut.Dtype() == rSrc.Dtype()) {
        const std::size_t rb = ncols * dtype_size(rOut.Dtype());
        std::memcpy(rOut.Data() + outRow0 * rb, rSrc.Data() + srcRow0 * rb, nrows * rb);
    } else {
        parallel_for(nrows, [&](std::size_t r) {
            for (std::size_t k = 0; k < ncols; ++k)
                merge_store_double(rOut, (outRow0 + r) * ncols + k,
                                   detail::read_double(rSrc, (srcRow0 + r) * ncols + k));
        });
    }
}

// Fill `nrows` rows (`ncols` wide) of `rOut` starting at `outRow0` with NaN.
void merge_fill_nan(NDArray& rOut, std::size_t outRow0, std::size_t nrows, std::size_t ncols) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    parallel_for(nrows, [&](std::size_t r) {
        for (std::size_t k = 0; k < ncols; ++k)
            merge_store_double(rOut, (outRow0 + r) * ncols + k, nan);
    });
}

// --- spatial hash (weld) ----------------------------------------------------
// The bucket grid itself lives in detail/spatial_hash.hpp (hoisted from here
// in v7.13.0, shared with operations/interpolate.cpp); cell size = atol, so
// points within atol of each other fall in the same or an adjacent cell (the
// 3x3x3 neighbourhood searched during dedup).

// Build the weld remap over the (Float64, global-order) `rPoints`: fills
// `rRemap` (global point index -> output point index) and returns, per output
// point, the global index of its representative (first contributing point,
// keep-first). `dim` is the coordinate dimension (<= 3 used for distance).
std::vector<std::int64_t> merge_build_weld_map(const NDArray& rPoints, std::size_t total,
                                               std::size_t dim, double atol,
                                               std::vector<std::int64_t>& rRemap) {
    const double* pts = rPoints.As<double>();
    const std::size_t ddim = std::min<std::size_t>(dim, 3);

    // Phase 1 (parallel): quantized bucket key per point.
    std::vector<detail::GridKey> keys(total);
    parallel_for(total, [&](std::size_t g) {
        double c[3] = {0.0, 0.0, 0.0};
        for (std::size_t d = 0; d < ddim; ++d)
            c[d] = pts[g * dim + d];
        keys[g] =
            detail::GridKey{detail::grid_quantize(c[0], atol), detail::grid_quantize(c[1], atol),
                            detail::grid_quantize(c[2], atol)};
    });

    // Phase 2 (serial, deterministic): first occurrence in a bucket wins; a
    // later point within atol of an existing representative welds onto it.
    rRemap.assign(total, -1);
    std::vector<std::int64_t> rep_global;
    rep_global.reserve(total);
    detail::SpatialGrid grid(atol);
    const double atol2 = atol * atol;

    for (std::size_t g = 0; g < total; ++g) {
        const detail::GridKey k = keys[g];
        std::int64_t found = -1;
        grid.ForEachIn27(k, [&](const std::vector<std::int64_t>& rReps) {
            for (std::int64_t rep : rReps) {
                const std::int64_t rg = rep_global[static_cast<std::size_t>(rep)];
                double d2 = 0.0;
                for (std::size_t d = 0; d < ddim; ++d) {
                    const double delta =
                        pts[g * dim + d] - pts[static_cast<std::size_t>(rg) * dim + d];
                    d2 += delta * delta;
                }
                if (d2 <= atol2) {
                    found = rep;
                    return false;
                }
            }
            return true;
        });
        if (found >= 0) {
            rRemap[g] = found;
        } else {
            const std::int64_t nidx = static_cast<std::int64_t>(rep_global.size());
            rRemap[g] = nidx;
            rep_global.push_back(static_cast<std::int64_t>(g));
            grid.Insert(k, nidx);
        }
    }
    return rep_global;
}

// Gather `new_count` rows from `rSrc` by `rRepGlobal[n]` -> output row n (a
// bandwidth-bound byte gather). The output mirrors `rSrc`'s shape with the
// first dimension replaced by `new_count` (so 1-D scalar data stays 1-D).
NDArray merge_gather_rows(const NDArray& rSrc, const std::vector<std::int64_t>& rRepGlobal,
                          std::size_t new_count) {
    std::vector<std::size_t> shape = rSrc.Shape();
    if (shape.empty())
        shape = {new_count};
    else
        shape[0] = new_count;
    NDArray dst = NDArray::Uninit(rSrc.Dtype(), shape);
    const std::size_t src_rows = detail::rows(rSrc);
    if (new_count == 0 || src_rows == 0)
        return dst;
    const std::size_t rb = rSrc.Nbytes() / src_rows;
    const std::byte* src = rSrc.Data();
    std::byte* out = dst.Data();
    parallel_for_bw(new_count, [&](std::size_t n) {
        std::memcpy(out + n * rb, src + static_cast<std::size_t>(rRepGlobal[n]) * rb, rb);
    });
    return dst;
}

// --- data-key planning ------------------------------------------------------

// Resolve the output dtype for a key: the common dtype if every contributing
// array shares it (and no NaN fill is needed), else Float64.
DType merge_common_dtype(const std::vector<const NDArray*>& rArrays, bool any_missing) {
    if (any_missing)
        return DType::Float64;
    DType dt = DType::Float64;
    bool first = true;
    for (const NDArray* a : rArrays) {
        if (!a)
            continue;
        if (first) {
            dt = a->Dtype();
            first = false;
        } else if (a->Dtype() != dt) {
            return DType::Float64;
        }
    }
    return dt;
}

}  // namespace

// --- public API -------------------------------------------------------------

MergeResult merge(const std::vector<const Mesh*>& rMeshes, const MergeOptions& rOpts) {
    if (rMeshes.empty())
        throw std::invalid_argument("merge: at least one input mesh is required");
    for (const Mesh* m : rMeshes)
        if (m == nullptr)
            throw std::invalid_argument("merge: input mesh pointer is null");

    const std::size_t nm = rMeshes.size();
    MergeResult res;
    Mesh& out = res.mMesh;

    // --- point offsets + output dimension -----------------------------------
    std::vector<std::size_t> np(nm), poff(nm + 1, 0);
    std::size_t dim = 0;
    for (std::size_t m = 0; m < nm; ++m) {
        np[m] = rMeshes[m]->NumPoints();
        poff[m + 1] = poff[m] + np[m];
        dim = std::max(dim, rMeshes[m]->PointDim());
    }
    const std::size_t total_pts = poff[nm];

    // --- concatenated points (global order, Float64) ------------------------
    NDArray concat_pts;
    if (total_pts > 0 && dim > 0) {
        concat_pts = NDArray::Uninit(DType::Float64, {total_pts, dim});
        double* dst = concat_pts.As<double>();
        for (std::size_t m = 0; m < nm; ++m) {
            const NDArray& src = rMeshes[m]->Points();
            const std::size_t pdim = rMeshes[m]->PointDim();
            const std::size_t base = poff[m];
            parallel_for(np[m], [&](std::size_t i) {
                for (std::size_t d = 0; d < dim; ++d)
                    dst[(base + i) * dim + d] =
                        (d < pdim) ? detail::read_double(src, i * pdim + d) : 0.0;
            });
        }
    }

    // --- weld remap (or identity) -------------------------------------------
    const bool do_weld = rOpts.weld && total_pts > 0 && dim > 0 && rOpts.atol > 0.0;
    std::vector<std::int64_t> remap;       // global point index -> output point index
    std::vector<std::int64_t> rep_global;  // output point index -> representative global index
    std::size_t out_pts = total_pts;
    if (do_weld) {
        rep_global = merge_build_weld_map(concat_pts, total_pts, dim, rOpts.atol, remap);
        out_pts = rep_global.size();
    }
    auto new_index = [&](std::size_t global) -> std::int64_t {
        return do_weld ? remap[global] : static_cast<std::int64_t>(global);
    };

    // Output points: gather representatives when welding, else the concatenation.
    if (total_pts > 0 && dim > 0) {
        if (do_weld)
            out.AssignPoints(merge_gather_rows(concat_pts, rep_global, out_pts));
        else
            out.AssignPoints(std::move(concat_pts));
    }

    // Per-input point map (local point index -> output point index).
    res.mPointMaps.reserve(nm);
    for (std::size_t m = 0; m < nm; ++m) {
        NDArray pm = NDArray::Uninit(DType::Int64, {np[m]});
        std::int64_t* d = pm.As<std::int64_t>();
        parallel_for_bw(np[m], [&](std::size_t i) { d[i] = new_index(poff[m] + i); });
        res.mPointMaps.push_back(std::move(pm));
    }

    // --- point_data ---------------------------------------------------------
    {
        // Union of names, sorted (each mesh already returns sorted names).
        std::map<std::string, std::vector<const NDArray*>>
            by_name;  // name -> per-mesh array/nullptr
        for (std::size_t m = 0; m < nm; ++m) {
            for (const std::string& name : rMeshes[m]->PointDataNames()) {
                auto& v = by_name[name];
                if (v.empty())
                    v.assign(nm, nullptr);
                v[m] = &rMeshes[m]->PointData(name);
            }
        }
        for (auto& [name, arrays] : by_name) {
            // Decide cols from the first present, valid (rows == np[m]) mesh.
            std::size_t cols = 0;
            bool cols_set = false, cols_ok = true, any_missing = false;
            for (std::size_t m = 0; m < nm; ++m) {
                const NDArray* a = arrays[m];
                if (!a || detail::rows(*a) != np[m]) {
                    if (a && detail::rows(*a) != np[m])
                        cols_ok = false;  // present but misaligned -> cannot merge
                    arrays[m] = nullptr;
                    any_missing = true;
                    continue;
                }
                const std::size_t c = detail::cols(*a);
                if (!cols_set) {
                    cols = c;
                    cols_set = true;
                } else if (c != cols) {
                    cols_ok = false;
                }
            }
            if (!cols_set || !cols_ok) {
                log::warn("merge: dropping point_data '{}' (missing/misaligned/shape mismatch)",
                          name);
                continue;
            }
            if (any_missing && rOpts.data_policy == MergeDataPolicy::Intersection) {
                log::warn("merge: dropping partial point_data '{}' (intersection policy)", name);
                continue;
            }
            // Emit 1-D scalar data as 1-D (mirror the inputs), 2-D otherwise.
            bool one_d = true;
            for (std::size_t m = 0; m < nm; ++m)
                if (arrays[m] && arrays[m]->Ndim() > 1)
                    one_d = false;
            const DType dt = merge_common_dtype(arrays, any_missing);
            std::vector<std::size_t> shape = (one_d && cols == 1)
                                                 ? std::vector<std::size_t>{total_pts}
                                                 : std::vector<std::size_t>{total_pts, cols};
            NDArray concat = NDArray::Uninit(dt, shape);
            for (std::size_t m = 0; m < nm; ++m) {
                if (arrays[m])
                    merge_copy_rows(concat, poff[m], *arrays[m], 0, np[m], cols);
                else
                    merge_fill_nan(concat, poff[m], np[m], cols);
            }
            if (do_weld)
                out.AddPointData(name, merge_gather_rows(concat, rep_global, out_pts));
            else
                out.AddPointData(name, std::move(concat));
        }
    }

    // --- cell blocks: group by type, record contributions -------------------
    struct Contribution {
        std::size_t mesh;
        std::size_t in_block;
        std::size_t ncells;
        std::size_t start;  // pre-dedup start row within the output block
    };
    struct OutBlock {
        std::string type;
        bool poly = false;
        bool ragged = false;
        std::size_t npc = 0;
        std::vector<Contribution> contribs;
        std::size_t pre_count = 0;
    };
    std::vector<OutBlock> blocks;
    std::unordered_map<std::string, std::size_t> type_index;
    // Per-mesh input-global-cell base per block (for the cell maps).
    std::vector<std::vector<std::size_t>> in_block_base(nm);
    std::vector<std::size_t> in_ncells(nm, 0);
    for (std::size_t m = 0; m < nm; ++m) {
        const std::size_t nb = rMeshes[m]->NumCellBlocks();
        in_block_base[m].resize(nb + 1, 0);
        for (std::size_t b = 0; b < nb; ++b) {
            Mesh::CellView cv = rMeshes[m]->Cells(b);
            in_block_base[m][b + 1] = in_block_base[m][b] + cv.NumCells();
            const std::string type(cv.Type());
            std::size_t oi;
            auto it = type_index.find(type);
            if (it == type_index.end()) {
                oi = blocks.size();
                type_index.emplace(type, oi);
                OutBlock ob;
                ob.type = type;
                ob.poly = cv.IsPolyhedron();
                ob.ragged = !ob.poly && cv.IsRagged();
                ob.npc = cv.NodesPerCell();
                blocks.push_back(std::move(ob));
            } else {
                oi = it->second;
            }
            OutBlock& ob = blocks[oi];
            Contribution con{m, b, cv.NumCells(), ob.pre_count};
            ob.pre_count += cv.NumCells();
            ob.contribs.push_back(con);
        }
        in_ncells[m] = in_block_base[m][nb];
    }

    // --- build each output block's remapped connectivity + source ids -------
    struct BuiltBlock {
        NDArray conn;                                                   // rectangular
        std::vector<std::vector<std::int64_t>> polyrows;                // 1-level ragged
        std::vector<std::vector<std::vector<std::int64_t>>> polycells;  // 2-level ragged
        std::vector<std::int64_t> source_ids;                           // pre-dedup, size pre_count
        std::vector<std::int64_t> finalpos;      // pre-local -> final-local or -1
        std::vector<std::int64_t> final_to_pre;  // final-local -> kept pre-local
        std::size_t final_count = 0;
    };
    std::vector<BuiltBlock> built(blocks.size());

    for (std::size_t oi = 0; oi < blocks.size(); ++oi) {
        const OutBlock& ob = blocks[oi];
        BuiltBlock& bb = built[oi];
        bb.source_ids.resize(ob.pre_count);

        auto remap_node = [&](std::size_t mesh, std::int64_t local) -> std::int64_t {
            if (local < 0 || static_cast<std::size_t>(local) >= np[mesh])
                return local;
            return new_index(poff[mesh] + static_cast<std::size_t>(local));
        };

        if (ob.poly) {
            bb.polycells.resize(ob.pre_count);
            for (const Contribution& con : ob.contribs) {
                Mesh::CellView cv = rMeshes[con.mesh]->Cells(con.in_block);
                parallel_for_bw(con.ncells, [&](std::size_t c) {
                    const std::size_t p = con.start + c;
                    bb.source_ids[p] = static_cast<std::int64_t>(con.mesh);
                    bb.polycells[p].resize(cv.NumFaces(c));
                    for (std::size_t f = 0; f < cv.NumFaces(c); ++f) {
                        auto face = cv.Face(c, f);
                        bb.polycells[p][f].reserve(face.second);
                        for (std::size_t k = 0; k < face.second; ++k)
                            bb.polycells[p][f].push_back(remap_node(con.mesh, face.first[k]));
                    }
                });
            }
        } else if (ob.ragged) {
            bb.polyrows.resize(ob.pre_count);
            for (const Contribution& con : ob.contribs) {
                Mesh::CellView cv = rMeshes[con.mesh]->Cells(con.in_block);
                parallel_for_bw(con.ncells, [&](std::size_t c) {
                    const std::size_t p = con.start + c;
                    bb.source_ids[p] = static_cast<std::int64_t>(con.mesh);
                    const std::int64_t* row = cv.Row(c);
                    const std::size_t sz = cv.RowSize(c);
                    bb.polyrows[p].reserve(sz);
                    for (std::size_t k = 0; k < sz; ++k)
                        bb.polyrows[p].push_back(remap_node(con.mesh, row[k]));
                });
            }
        } else {
            const std::size_t npc = ob.npc;
            bb.conn = NDArray::Uninit(DType::Int64, {ob.pre_count, npc});
            std::int64_t* dst = bb.conn.As<std::int64_t>();
            for (const Contribution& con : ob.contribs) {
                Mesh::CellView cv = rMeshes[con.mesh]->Cells(con.in_block);
                const NDArray& conn = cv.Conn();
                const std::size_t start = con.start;
                const std::size_t mesh = con.mesh;
                parallel_for_bw(con.ncells, [&](std::size_t c) {
                    const std::size_t p = start + c;
                    bb.source_ids[p] = static_cast<std::int64_t>(mesh);
                    for (std::size_t k = 0; k < npc; ++k)
                        dst[p * npc + k] = remap_node(mesh, detail::read_int(conn, c * npc + k));
                });
            }
        }

        // Dedup (keep-first) if requested; else identity.
        bb.finalpos.assign(ob.pre_count, 0);
        if (rOpts.weld && rOpts.drop_duplicate_cells) {
            std::map<std::vector<std::int64_t>, std::int64_t> seen;
            for (std::size_t p = 0; p < ob.pre_count; ++p) {
                std::vector<std::int64_t> key;
                if (ob.poly) {
                    for (const auto& face : bb.polycells[p])
                        key.insert(key.end(), face.begin(), face.end());
                } else if (ob.ragged) {
                    key = bb.polyrows[p];
                } else {
                    const std::int64_t* dst = bb.conn.As<std::int64_t>();
                    key.assign(dst + p * ob.npc, dst + (p + 1) * ob.npc);
                }
                std::sort(key.begin(), key.end());
                auto it = seen.find(key);
                if (it != seen.end()) {
                    bb.finalpos[p] = -1;  // dropped duplicate
                } else {
                    const std::int64_t fidx = static_cast<std::int64_t>(bb.final_to_pre.size());
                    seen.emplace(std::move(key), fidx);
                    bb.finalpos[p] = fidx;
                    bb.final_to_pre.push_back(static_cast<std::int64_t>(p));
                }
            }
        } else {
            bb.final_to_pre.resize(ob.pre_count);
            for (std::size_t p = 0; p < ob.pre_count; ++p) {
                bb.finalpos[p] = static_cast<std::int64_t>(p);
                bb.final_to_pre[p] = static_cast<std::int64_t>(p);
            }
        }
        bb.final_count = bb.final_to_pre.size();
    }

    // Output block base offsets (global output cell index).
    std::vector<std::size_t> out_block_base(blocks.size() + 1, 0);
    for (std::size_t oi = 0; oi < blocks.size(); ++oi)
        out_block_base[oi + 1] = out_block_base[oi] + built[oi].final_count;

    // --- emit cell blocks (filtered by dedup) -------------------------------
    for (std::size_t oi = 0; oi < blocks.size(); ++oi) {
        const OutBlock& ob = blocks[oi];
        BuiltBlock& bb = built[oi];
        const std::size_t fc = bb.final_count;
        if (ob.poly) {
            std::vector<std::vector<std::vector<std::int64_t>>> cells(fc);
            parallel_for_bw(fc, [&](std::size_t f) {
                cells[f] = std::move(bb.polycells[static_cast<std::size_t>(bb.final_to_pre[f])]);
            });
            out.AddPolyhedronBlock(ob.type, std::move(cells));
        } else if (ob.ragged) {
            std::vector<std::vector<std::int64_t>> rows(fc);
            parallel_for_bw(fc, [&](std::size_t f) {
                rows[f] = std::move(bb.polyrows[static_cast<std::size_t>(bb.final_to_pre[f])]);
            });
            out.AddPolygonBlock(ob.type, std::move(rows));
        } else {
            NDArray conn = NDArray::Uninit(DType::Int64, {fc, ob.npc});
            std::int64_t* dst = conn.As<std::int64_t>();
            const std::int64_t* src = bb.conn.As<std::int64_t>();
            parallel_for_bw(fc, [&](std::size_t f) {
                std::memcpy(dst + f * ob.npc,
                            src + static_cast<std::size_t>(bb.final_to_pre[f]) * ob.npc,
                            ob.npc * sizeof(std::int64_t));
            });
            out.AddCellBlock(ob.type, std::move(conn));
        }
    }

    // --- source_mesh_id cell_data -------------------------------------------
    if (rOpts.source_tag && !blocks.empty()) {
        std::vector<NDArray> sblocks;
        sblocks.reserve(blocks.size());
        for (std::size_t oi = 0; oi < blocks.size(); ++oi) {
            BuiltBlock& bb = built[oi];
            const std::size_t fc = bb.final_count;
            NDArray a = NDArray::Uninit(DType::Int64, {fc});
            std::int64_t* d = a.As<std::int64_t>();
            parallel_for_bw(fc, [&](std::size_t f) {
                d[f] = bb.source_ids[static_cast<std::size_t>(bb.final_to_pre[f])];
            });
            sblocks.push_back(std::move(a));
        }
        out.AddCellData("source_mesh_id", std::move(sblocks));
    }

    // --- cell_data ----------------------------------------------------------
    {
        std::map<std::string, char> names;  // union of cell_data names
        for (std::size_t m = 0; m < nm; ++m)
            for (const std::string& name : rMeshes[m]->CellDataNames())
                names[name] = 1;
        for (const auto& [name, _] : names) {
            // Per contribution (in emit order = block-major, contrib order), the
            // source array (or nullptr) for this key, plus its cols check.
            bool any_missing = false, cols_ok = true, cols_set = false;
            std::size_t cols = 0;
            std::vector<const NDArray*> present;  // for common-dtype resolution
            for (std::size_t oi = 0; oi < blocks.size(); ++oi) {
                for (const Contribution& con : blocks[oi].contribs) {
                    const Mesh* mesh = rMeshes[con.mesh];
                    const NDArray* a = nullptr;
                    if (mesh->HasCellData(name) && con.in_block < mesh->CellDataNumBlocks(name)) {
                        const NDArray& cand = mesh->CellData(name, con.in_block);
                        if (detail::rows(cand) == con.ncells)
                            a = &cand;
                    }
                    if (!a) {
                        any_missing = true;
                        continue;
                    }
                    present.push_back(a);
                    const std::size_t c = detail::cols(*a);
                    if (!cols_set) {
                        cols = c;
                        cols_set = true;
                    } else if (c != cols) {
                        cols_ok = false;
                    }
                }
            }
            if (!cols_set || !cols_ok) {
                log::warn("merge: dropping cell_data '{}' (missing/misaligned/shape mismatch)",
                          name);
                continue;
            }
            if (any_missing && rOpts.data_policy == MergeDataPolicy::Intersection) {
                log::warn("merge: dropping partial cell_data '{}' (intersection policy)", name);
                continue;
            }
            bool one_d = true;
            for (const NDArray* a : present)
                if (a->Ndim() > 1)
                    one_d = false;
            const DType dt = merge_common_dtype(present, any_missing);
            auto make_shape = [&](std::size_t n) {
                return (one_d && cols == 1) ? std::vector<std::size_t>{n}
                                            : std::vector<std::size_t>{n, cols};
            };

            std::vector<NDArray> outblocks;
            outblocks.reserve(blocks.size());
            for (std::size_t oi = 0; oi < blocks.size(); ++oi) {
                const OutBlock& ob = blocks[oi];
                BuiltBlock& bb = built[oi];
                NDArray pre = NDArray::Uninit(dt, make_shape(ob.pre_count));
                for (const Contribution& con : ob.contribs) {
                    const Mesh* mesh = rMeshes[con.mesh];
                    const NDArray* a = nullptr;
                    if (mesh->HasCellData(name) && con.in_block < mesh->CellDataNumBlocks(name)) {
                        const NDArray& cand = mesh->CellData(name, con.in_block);
                        if (detail::rows(cand) == con.ncells)
                            a = &cand;
                    }
                    if (a)
                        merge_copy_rows(pre, con.start, *a, 0, con.ncells, cols);
                    else
                        merge_fill_nan(pre, con.start, con.ncells, cols);
                }
                // Filter by dedup keep-first.
                NDArray fin = NDArray::Uninit(dt, make_shape(bb.final_count));
                parallel_for_bw(bb.final_count, [&](std::size_t f) {
                    merge_copy_rows(fin, f, pre, static_cast<std::size_t>(bb.final_to_pre[f]), 1,
                                    cols);
                });
                outblocks.push_back(std::move(fin));
            }
            out.AddCellData(name, std::move(outblocks));
        }
    }

    // --- field_data: namespace on collision ---------------------------------
    {
        std::map<std::string, int> counts;
        for (std::size_t m = 0; m < nm; ++m)
            for (const std::string& name : rMeshes[m]->FieldDataNames())
                counts[name] += 1;
        for (std::size_t m = 0; m < nm; ++m) {
            for (const std::string& name : rMeshes[m]->FieldDataNames()) {
                const std::string key = counts[name] > 1 ? (std::to_string(m) + ":" + name) : name;
                out.AddFieldData(key, merge_owned_copy(rMeshes[m]->FieldData(name)));
            }
        }
    }

    // --- per-input cell maps (input global cell -> output global cell) ------
    res.mCellMaps.reserve(nm);
    for (std::size_t m = 0; m < nm; ++m)
        res.mCellMaps.emplace_back(NDArray::Uninit(DType::Int64, {in_ncells[m]}));
    for (std::size_t m = 0; m < nm; ++m) {
        std::int64_t* d = res.mCellMaps[m].As<std::int64_t>();
        parallel_for_bw(in_ncells[m], [&](std::size_t i) { d[i] = -1; });
    }
    for (std::size_t oi = 0; oi < blocks.size(); ++oi) {
        const OutBlock& ob = blocks[oi];
        const BuiltBlock& bb = built[oi];
        for (const Contribution& con : ob.contribs) {
            std::int64_t* d = res.mCellMaps[con.mesh].As<std::int64_t>();
            const std::size_t gbase = in_block_base[con.mesh][con.in_block];
            parallel_for_bw(con.ncells, [&](std::size_t c) {
                const std::int64_t fl = bb.finalpos[con.start + c];
                d[gbase + c] = (fl < 0) ? -1 : static_cast<std::int64_t>(out_block_base[oi]) + fl;
            });
        }
    }

    // --- named regions: remap per input, namespacing on collision -----------
    // Same rule as field_data above: a name carried by more than one input is
    // prefixed with its source index, so "0:wall" and "1:wall" stay distinct
    // groups instead of the second silently replacing the first (they would
    // share a `(kind, name, dim, tag)` region key). Renaming has to happen
    // between remapping and storing, which is why this uses the per-region
    // `remap_region` rather than the whole-mesh `remap_regions`. The map is
    // Global: merge numbers cells across each whole input mesh, not per block.
    {
        std::map<std::string, int> region_counts;
        for (std::size_t m = 0; m < nm; ++m)
            for (const std::string& name : rMeshes[m]->RegionNames())
                region_counts[name] += 1;
        for (std::size_t m = 0; m < nm; ++m) {
            detail::RegionRemap rmap;
            rmap.pPointMap = &res.mPointMaps[m];
            rmap.mCellMapKind = detail::CellMapKind::Global;
            rmap.pGlobalCellMap = &res.mCellMaps[m];
            rmap.mOpName = "merge";
            for (std::size_t i = 0; i < rMeshes[m]->NumRegions(); ++i) {
                const Region& r_src = rMeshes[m]->Region(i);
                Region carried;
                if (!detail::remap_region(*rMeshes[m], out, r_src, rmap, carried))
                    continue;
                if (region_counts[r_src.mName] > 1)
                    carried.mName = std::to_string(m) + ":" + r_src.mName;
                out.AddRegion(std::move(carried));
            }
        }
    }

    return res;
}

}  // namespace meshioplusplus
