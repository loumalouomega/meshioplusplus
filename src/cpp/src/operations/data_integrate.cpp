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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// Cell-measure-weighted field integration: see operations/data_integrate.hpp
// for the contract. Built entirely through the uniform mesh API, so it
// compiles under every mesh backend.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/data_integrate.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

// Project includes (private, not installed)
#include "../detail/quadrature.hpp"
#include "../detail/reference_element.hpp"

namespace meshioplusplus {

namespace {

// The cells belonging to one named Cell region, resolved to (block, row)
// pairs once and reused across every requested array.
struct DintRegionCells {
    std::string mName;
    std::vector<std::pair<std::size_t, std::int64_t>> mCells;
};

// Turns per-component accumulators into the public report shape.
FieldIntegralRegion dint_build_region(std::string name,
                                      const std::vector<detail::WeightedSum>& rStats,
                                      std::int64_t NumCells, std::int64_t NumSkipped,
                                      std::size_t Comps) {
    FieldIntegralRegion out;
    out.mName = std::move(name);
    out.mNumCells = NumCells;
    out.mNumSkipped = NumSkipped;
    out.mDomainMeasurePerComponent.resize(Comps);
    out.mTotalPerComponent.resize(Comps);
    out.mMeanPerComponent.resize(Comps);
    out.mNumNanPerComponent.resize(Comps);
    for (std::size_t k = 0; k < Comps; ++k) {
        out.mDomainMeasurePerComponent[k] = rStats[k].mSumW;
        out.mTotalPerComponent[k] = rStats[k].mSumVW;
        out.mMeanPerComponent[k] = rStats[k].Mean();
        out.mNumNanPerComponent[k] = rStats[k].mNumNan;
    }
    return out;
}

std::string dint_unknown_or_point_message(const Mesh& rMesh, const std::string& rName) {
    if (rMesh.HasPointData(rName))
        return "meshio++: data_integrate: '" + rName +
               "' is a point_data array; convert it first with point_data_to_cell_data (CLI: "
               "`data to-cell`)";
    return data_unknown_key_message(rMesh, DataLocation::Cell, rName);
}

}  // namespace

static DataIntegrateReport dint_core(const Mesh& rMesh, const DataIntegrateOptions& rOptions) {
    // --- resolve which arrays to integrate (validation before any work) -----
    std::vector<std::string> names = rOptions.mArrayNames;
    if (names.empty()) {
        names = rMesh.CellDataNames();  // sorted on every backend
    } else {
        for (const std::string& name : names)
            if (!rMesh.HasCellData(name))
                throw std::invalid_argument(dint_unknown_or_point_message(rMesh, name));
    }
    for (const std::string& name : names)
        if (rMesh.CellDataNumBlocks(name) != rMesh.NumCellBlocks())
            throw std::invalid_argument("meshio++: data_integrate: cell_data '" + name +
                                        "' does not have one array per cell block");

    DataIntegrateReport report;
    if (names.empty())
        return report;

    const std::size_t nblocks = rMesh.NumCellBlocks();

    // --- per-cell measures, computed once and shared by every array/region --
    // A negative sentinel marks a cell whose measure could not be computed
    // (ragged/unsupported/degenerate), so every array's accumulation excludes
    // it identically without recomputing anything geometric.
    std::vector<std::vector<double>> measures(nblocks);
    std::int64_t domain_num_cells = 0;
    std::int64_t domain_num_skipped = 0;
    if (nblocks > 0) {
        const NDArray& points = rMesh.Points();
        const std::size_t pdim = rMesh.PointDim();
        for (std::size_t b = 0; b < nblocks; ++b) {
            const auto cb = rMesh.Cells(b);
            std::vector<double>& w = measures[b];
            w.assign(cb.NumCells(), 0.0);
            parallel_for(cb.NumCells(), [&](std::size_t c) {
                w[c] = std::abs(detail::cell_measure(points, pdim, cb, c));
            });
            for (double& v : w) {
                if (std::isfinite(v) && v > 0.0)
                    ++domain_num_cells;
                else {
                    v = -1.0;  // sentinel: unmeasurable
                    ++domain_num_skipped;
                }
            }
        }
    }

    // --- named Cell regions, resolved once -----------------------------------
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    std::vector<DintRegionCells> region_defs;
    for (const std::string& rname : rMesh.RegionNames()) {
        const std::size_t idx = rMesh.FindRegion(rname, RegionKind::Cell);
        if (idx == Mesh::npos)
            continue;
        const Region& region = rMesh.Region(idx);
        DintRegionCells rc;
        rc.mName = rname;
        const std::int64_t* entries = region.Entries();
        const std::size_t n = region.NumEntries();
        rc.mCells.reserve(n);
        for (std::size_t e = 0; e < n; ++e) {
            const auto [b, row] = detail::global_to_block_row(bases, entries[e]);
            if (b != static_cast<std::size_t>(-1))
                rc.mCells.push_back({b, row});
        }
        region_defs.push_back(std::move(rc));
    }
    std::vector<std::int64_t> region_num_cells(region_defs.size(), 0);
    std::vector<std::int64_t> region_num_skipped(region_defs.size(), 0);
    for (std::size_t r = 0; r < region_defs.size(); ++r)
        for (const auto& [b, row] : region_defs[r].mCells) {
            if (measures[b][static_cast<std::size_t>(row)] > 0.0)
                ++region_num_cells[r];
            else
                ++region_num_skipped[r];
        }

    // --- per array: whole-mesh, then every region ----------------------------
    report.mArrays.reserve(names.size());
    for (const std::string& name : names) {
        FieldIntegralArray arr;
        arr.mName = name;
        const std::size_t comps = data_num_components(rMesh.CellData(name, 0));
        arr.mNumComponents = static_cast<std::int64_t>(comps);

        std::vector<detail::WeightedSum> domain_stats(comps);
        for (std::size_t b = 0; b < nblocks; ++b)
            detail::accumulate_weighted(rMesh.CellData(name, b), comps, measures[b], domain_stats);
        arr.mDomain =
            dint_build_region("", domain_stats, domain_num_cells, domain_num_skipped, comps);

        arr.mRegions.reserve(region_defs.size());
        for (std::size_t r = 0; r < region_defs.size(); ++r) {
            // Regions are typically few and small relative to the whole mesh,
            // so a plain serial fold (in the region's own canonical, sorted
            // entry order) is simpler than chunked parallelism here and is
            // trivially deterministic.
            std::vector<detail::WeightedSum> region_stats(comps);
            for (const auto& [b, row] : region_defs[r].mCells) {
                const double w = measures[b][static_cast<std::size_t>(row)];
                if (!(w > 0.0))
                    continue;
                const NDArray& block_arr = rMesh.CellData(name, b);
                for (std::size_t k = 0; k < comps; ++k)
                    region_stats[k].Add(
                        detail::read_double(block_arr, static_cast<std::size_t>(row) * comps + k),
                        w);
            }
            arr.mRegions.push_back(dint_build_region(region_defs[r].mName, region_stats,
                                                     region_num_cells[r], region_num_skipped[r],
                                                     comps));
        }

        report.mArrays.push_back(std::move(arr));
    }

    return report;
}

// ---------------------------------------------------------------------------
// Integration modes carried by the array name
// ---------------------------------------------------------------------------
//
// An array name may carry a mode prefix, so no signature (and no binding)
// changes: `point:<a>` integrates a point_data array (its cell means, exact
// for a linear field on a simplex); `flux:<a>` integrates the flux `u . n` of a
// 3-component cell_data or point_data array through the surface cells (the
// unit cell normals of `compute_normals`); `boundary:` before either takes the
// extracted boundary of a volume mesh instead of the mesh. A cell_data array
// whose name is the whole string wins, so no existing name changes meaning.

namespace {

struct DintMode {
    bool mBoundary = false;
    bool mPoint = false;
    bool mFlux = false;
    bool mQuad = false;
    std::string mBase;
    bool Plain() const { return !mBoundary && !mPoint && !mFlux && !mQuad; }
};

bool dint_starts_with(const std::string& s, const char* p) {
    return s.rfind(p, 0) == 0;
}

DintMode dint_parse(const Mesh& rMesh, const std::string& rName) {
    DintMode m;
    m.mBase = rName;
    if (rMesh.HasCellData(rName))
        return m;  // an array of exactly this name wins
    std::string s = rName;
    if (dint_starts_with(s, "boundary:")) {
        m.mBoundary = true;
        s = s.substr(9);
    }
    if (dint_starts_with(s, "point:")) {
        m.mPoint = true;
        s = s.substr(6);
    } else if (dint_starts_with(s, "flux:")) {
        m.mFlux = true;
        s = s.substr(5);
    } else if (dint_starts_with(s, "quad:")) {
        m.mQuad = true;
        s = s.substr(5);
    }
    if (m.Plain())
        return m;  // no prefix: the unknown-name error follows
    if (m.mBoundary && !m.mPoint && !m.mFlux && !m.mQuad)
        throw std::invalid_argument(
            "meshio++: data_integrate: 'boundary:' needs 'point:', 'quad:' or 'flux:' after it "
            "(the boundary surface carries no cell_data)");
    m.mBase = s;
    return m;
}

[[noreturn]] void dint_fail(const std::string& rMessage) {
    throw std::invalid_argument("meshio++: data_integrate: " + rMessage);
}

// Integrate a point_data array by Gauss quadrature on each cell's shape functions.
FieldIntegralArray dint_quad_array(const Mesh& rMesh, const std::string& rBaseName,
                                   const std::string& rReportName) {
    if (!rMesh.HasPointData(rBaseName))
        dint_fail("'quad:" + rBaseName + "' names no point_data array");
    const NDArray& pdata = rMesh.PointData(rBaseName);
    const std::size_t comps = data_num_components(pdata);
    const std::size_t nblocks = rMesh.NumCellBlocks();
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();

    // Per-cell quadrature measure and integral (block-major).
    std::vector<std::vector<double>> cell_vol(nblocks);
    std::vector<std::vector<double>> cell_int(nblocks);
    std::vector<std::vector<std::uint8_t>> cell_nan(nblocks);
    std::int64_t domain_num_cells = 0;
    std::int64_t domain_num_skipped = 0;

    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto cb = rMesh.Cells(b);
        const std::size_t nc = cb.NumCells();
        cell_vol[b].assign(nc, -1.0);
        cell_int[b].assign(nc * comps, 0.0);
        cell_nan[b].assign(nc * comps, 0);

        detail::ReferenceElementInfo info;
        const CellType ct = cell_type_from_name(cb.Type());
        if (cb.IsRagged() || !detail::reference_element_info(ct, cb.NodesPerCell(), info)) {
            domain_num_skipped += static_cast<std::int64_t>(nc);
            continue;
        }
        const detail::QuadratureRule rule =
            detail::quadrature_rule(ct, std::max(2 * info.mOrder + 1, 2));
        if (rule.mPoints.empty()) {
            domain_num_skipped += static_cast<std::int64_t>(nc);
            continue;
        }
        const std::size_t n_nodes = cb.NodesPerCell();
        const NDArray& conn = cb.Conn();

        parallel_for(nc, [&](std::size_t c) {
            std::vector<detail::Vec3> cpts(n_nodes);
            std::vector<std::int64_t> cids(n_nodes);
            for (std::size_t i = 0; i < n_nodes; ++i) {
                const std::int64_t nid = detail::read_int(conn, c * n_nodes + i);
                cids[i] = nid;
                const std::size_t u = static_cast<std::size_t>(nid);
                cpts[i] = {detail::read_double(points, u * pdim + 0),
                           pdim > 1 ? detail::read_double(points, u * pdim + 1) : 0.0,
                           pdim > 2 ? detail::read_double(points, u * pdim + 2) : 0.0};
            }
            std::vector<double> vals;
            std::vector<detail::Vec3> grads;
            double vol = 0.0;
            std::vector<double> acc(comps, 0.0);
            std::vector<std::uint8_t> has_nan(comps, 0);
            for (std::size_t k = 0; k < comps; ++k) {
                for (std::size_t i = 0; i < n_nodes; ++i) {
                    const double fv =
                        detail::read_double(pdata, static_cast<std::size_t>(cids[i]) * comps + k);
                    if (!std::isfinite(fv)) {
                        has_nan[k] = 1;
                        break;
                    }
                }
            }
            for (std::size_t q = 0; q < rule.mPoints.size(); ++q) {
                double meas = 0.0;
                if (!detail::shape_eval_and_jacobian(ct, cpts, rule.mPoints[q], vals, grads,
                                                     meas)) {
                    vol = -1.0;
                    break;
                }
                const double w = rule.mWeights[q] * meas;
                vol += w;
                for (std::size_t k = 0; k < comps; ++k) {
                    if (has_nan[k])
                        continue;
                    double fq = 0.0;
                    for (std::size_t i = 0; i < n_nodes; ++i)
                        fq += vals[i] * detail::read_double(
                                            pdata, static_cast<std::size_t>(cids[i]) * comps + k);
                    acc[k] += w * fq;
                }
            }
            if (std::isfinite(vol) && vol > 0.0) {
                cell_vol[b][c] = vol;
                for (std::size_t k = 0; k < comps; ++k) {
                    cell_int[b][c * comps + k] = acc[k];
                    cell_nan[b][c * comps + k] = has_nan[k];
                }
            }
        });

        for (std::size_t c = 0; c < nc; ++c) {
            if (cell_vol[b][c] > 0.0)
                ++domain_num_cells;
            else
                ++domain_num_skipped;
        }
    }

    FieldIntegralArray arr;
    arr.mName = rReportName;
    arr.mNumComponents = static_cast<std::int64_t>(comps);

    std::vector<detail::WeightedSum> domain_stats(comps);
    for (std::size_t b = 0; b < nblocks; ++b) {
        const std::size_t nc = cell_vol[b].size();
        for (std::size_t c = 0; c < nc; ++c) {
            const double w = cell_vol[b][c];
            if (!(w > 0.0))
                continue;
            for (std::size_t k = 0; k < comps; ++k) {
                if (cell_nan[b][c * comps + k])
                    domain_stats[k].Add(std::numeric_limits<double>::quiet_NaN(), w);
                else
                    domain_stats[k].Add(cell_int[b][c * comps + k] / w, w);
            }
        }
    }
    arr.mDomain = dint_build_region("", domain_stats, domain_num_cells, domain_num_skipped, comps);

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    for (const std::string& rname : rMesh.RegionNames()) {
        const std::size_t idx = rMesh.FindRegion(rname, RegionKind::Cell);
        if (idx == Mesh::npos)
            continue;
        const Region& region = rMesh.Region(idx);
        std::vector<detail::WeightedSum> rstats(comps);
        std::int64_t rcells = 0;
        std::int64_t rskip = 0;
        const std::int64_t* entries = region.Entries();
        for (std::size_t e = 0; e < region.NumEntries(); ++e) {
            const auto [b, row] = detail::global_to_block_row(bases, entries[e]);
            if (b == static_cast<std::size_t>(-1))
                continue;
            const std::size_t c = static_cast<std::size_t>(row);
            const double w = cell_vol[b][c];
            if (!(w > 0.0)) {
                ++rskip;
                continue;
            }
            ++rcells;
            for (std::size_t k = 0; k < comps; ++k) {
                if (cell_nan[b][c * comps + k])
                    rstats[k].Add(std::numeric_limits<double>::quiet_NaN(), w);
                else
                    rstats[k].Add(cell_int[b][c * comps + k] / w, w);
            }
        }
        arr.mRegions.push_back(dint_build_region(rname, rstats, rcells, rskip, comps));
    }
    return arr;
}

}  // namespace

DataIntegrateReport data_integrate(const Mesh& rMesh, const DataIntegrateOptions& rOptions) {
    std::vector<DintMode> modes;
    bool any_mode = false;
    for (const std::string& name : rOptions.mArrayNames) {
        modes.push_back(dint_parse(rMesh, name));
        any_mode = any_mode || !modes.back().Plain();
    }
    if (!any_mode)
        return dint_core(rMesh, rOptions);

    // Working meshes, built only when a mode needs them: the mesh itself and its
    // boundary, each gaining one temporary cell_data array per mode entry.
    Mesh plain = detail::clone_mesh(rMesh);
    Mesh bound;
    bool have_bound = false;
    std::vector<std::string> temp_names(modes.size());
    std::vector<FieldIntegralArray> quad_results(modes.size());
    for (std::size_t i = 0; i < modes.size(); ++i) {
        const DintMode& m = modes[i];
        if (m.Plain())
            continue;
        Mesh* wk = &plain;
        if (m.mBoundary) {
            if (!have_bound) {
                if (!has_surface_extractable_cells(rMesh))
                    dint_fail("'boundary:' needs a mesh with cells extract_surface can take");
                bound = extract_surface(rMesh);
                have_bound = true;
            }
            wk = &bound;
        }
        if (m.mQuad) {
            quad_results[i] = dint_quad_array(*wk, m.mBase, rOptions.mArrayNames[i]);
            continue;
        }
        const std::string temp = "\x01integrate" + std::to_string(i);
        temp_names[i] = temp;
        const bool as_point = wk->HasPointData(m.mBase);
        const bool as_cell = !m.mPoint && wk->HasCellData(m.mBase);
        if (m.mPoint && !as_point)
            dint_fail("'point:" + m.mBase + "' names no point_data array");
        if (!as_point && !as_cell)
            dint_fail(data_unknown_key_message(*wk, DataLocation::Point, m.mBase));
        // The cell values: the cell_data array itself, or the cell means of the point_data one.
        std::string src = m.mBase;
        if (!(as_cell && !m.mPoint)) {
            DataAverageOptions avg;
            avg.names = {m.mBase};
            avg.prefix = "\x01u" + std::to_string(i) + ":";
            *wk = point_data_to_cell_data(*wk, avg);
            src = avg.prefix + m.mBase;
        }
        if (!m.mFlux) {
            // point:<a> -- the cell means, under the temporary name.
            std::vector<NDArray> blocks;
            for (std::size_t b = 0; b < wk->NumCellBlocks(); ++b)
                blocks.push_back(wk->CellData(src, b));
            wk->AddCellData(temp, std::move(blocks));
            continue;
        }
        // flux:<a> -- u . n per cell with the unit cell normals.
        NormalsOptions no;
        no.mPointNormals = false;
        no.mCellNormals = true;
        NormalsResult nr;
        try {
            nr = compute_normals(*wk, no);
        } catch (const std::exception& e) {
            dint_fail(std::string("flux needs a surface mesh: ") + e.what());
        }
        std::vector<NDArray> blocks;
        for (std::size_t b = 0; b < wk->NumCellBlocks(); ++b) {
            const NDArray& u = wk->CellData(src, b);
            if (data_num_components(u) != 3)
                dint_fail("flux:" + m.mBase + " needs a 3-component array (got " +
                          std::to_string(data_num_components(u)) + ")");
            const NDArray& nrm = nr.mMesh.CellData(kNormalsName, b);
            const std::size_t rows = detail::rows(u);
            NDArray dens = NDArray::Uninit(DType::Float64, {rows});
            double* d = dens.As<double>();
            for (std::size_t c = 0; c < rows; ++c) {
                double acc = 0.0;
                for (std::size_t k = 0; k < 3; ++k)
                    acc += detail::read_double(u, c * 3 + k) * detail::read_double(nrm, c * 3 + k);
                d[c] = acc;
            }
            blocks.push_back(std::move(dens));
        }
        wk->AddCellData(temp, std::move(blocks));
    }

    // Integrate each working mesh once over its temporary arrays and the plain
    // names, then restore the requested names in the requested order.
    DataIntegrateOptions plain_opts;
    DataIntegrateOptions bound_opts;
    for (std::size_t i = 0; i < modes.size(); ++i) {
        if (modes[i].mQuad)
            continue;
        if (modes[i].Plain())
            plain_opts.mArrayNames.push_back(rOptions.mArrayNames[i]);
        else if (modes[i].mBoundary)
            bound_opts.mArrayNames.push_back(temp_names[i]);
        else
            plain_opts.mArrayNames.push_back(temp_names[i]);
    }
    const DataIntegrateReport pr =
        plain_opts.mArrayNames.empty() ? DataIntegrateReport{} : dint_core(plain, plain_opts);
    const DataIntegrateReport br =
        bound_opts.mArrayNames.empty() ? DataIntegrateReport{} : dint_core(bound, bound_opts);
    DataIntegrateReport out;
    std::size_t pi = 0;
    std::size_t bi = 0;
    for (std::size_t i = 0; i < modes.size(); ++i) {
        if (modes[i].mQuad) {
            out.mArrays.push_back(std::move(quad_results[i]));
            continue;
        }
        FieldIntegralArray arr =
            (!modes[i].Plain() && modes[i].mBoundary) ? br.mArrays[bi++] : pr.mArrays[pi++];
        arr.mName = rOptions.mArrayNames[i];
        out.mArrays.push_back(std::move(arr));
    }
    return out;
}

}  // namespace meshioplusplus
