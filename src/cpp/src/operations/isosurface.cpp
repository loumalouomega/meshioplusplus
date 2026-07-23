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
// Isosurfaces / contours: the level set of a scalar point_data field, cut with
// the same marching-tetrahedra machinery as the planar cross-section
// (detail/marching.hpp). This file owns only what is specific to a *field* cut:
// resolving and reducing the array, looping the isovalues, forcing the contoured
// value exact, and concatenating the contours into one tagged mesh. See
// operations/isosurface.hpp for the contract the pure-numpy twin
// (_isosurface.py) replicates expression for expression.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/marching.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// Reduce one row to a scalar: the requested component, or the magnitude.
//
// Deliberately the same arithmetic as face_color.cpp's fc_scalarize (which is
// file-private there) -- the magnitude sums left to right and takes one sqrt at
// the end, so the numpy twin rounds alike rather than using np.linalg.norm.
double iso_scalarize(const NDArray& rArray, std::size_t Row, std::size_t NumComponents,
                     const std::optional<int>& rComponent) {
    const std::size_t base = Row * NumComponents;
    if (rComponent.has_value()) {
        const int c = *rComponent;
        if (c < 0 || static_cast<std::size_t>(c) >= NumComponents)
            throw std::invalid_argument("meshio++: isosurface: component " + std::to_string(c) +
                                        " is out of range for an array with " +
                                        std::to_string(NumComponents) + " component(s)");
        return detail::read_double(rArray, base + static_cast<std::size_t>(c));
    }
    if (NumComponents == 1)
        return detail::read_double(rArray, base);
    double sum = 0.0;
    for (std::size_t k = 0; k < NumComponents; ++k) {
        const double v = detail::read_double(rArray, base + k);
        sum += v * v;
    }
    return std::sqrt(sum);
}

// The index of the block of type `rType` in `rMesh`, if it has one. The cutter
// emits at most one block per type, in the fixed triangle/quad/line order.
std::optional<std::size_t> iso_type_block(const Mesh& rMesh, const std::string& rType) {
    std::size_t bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        if (std::string(cb.Type()) == rType)
            return bi;
        ++bi;
    }
    return std::nullopt;
}

// One kept contour: the section mesh plus the isovalue that produced it.
struct IsoContour {
    Mesh mMesh;
    double mValue;
    std::int64_t mIndex;
};

// Concatenate the contours into one mesh: points and point_data stacked in
// ascending-isovalue order, cell blocks merged by type in the canonical order,
// per-cell cell_data carried alongside, and the two contour tags attached.
//
// `rExactName`/`rExactComponent` name the contoured array whose value must read
// back as exactly the isovalue (empty = leave everything interpolated, which is
// the magnitude case).
Mesh iso_concat_contours(std::vector<IsoContour>& rContours, std::size_t Dim,
                         const std::string& rExactName, const std::optional<int>& rExactComponent) {
    Mesh out;

    // --- points --------------------------------------------------------------
    std::vector<std::size_t> point_base(rContours.size(), 0);
    std::size_t total_points = 0;
    for (std::size_t k = 0; k < rContours.size(); ++k) {
        point_base[k] = total_points;
        total_points += rContours[k].mMesh.NumPoints();
    }
    {
        NDArray points = NDArray::Uninit(DType::Float64, {total_points, Dim});
        double* dst = points.As<double>();
        for (std::size_t k = 0; k < rContours.size(); ++k) {
            const NDArray& src = rContours[k].mMesh.Points();
            const std::size_t n = rContours[k].mMesh.NumPoints();
            if (n != 0)
                std::memcpy(dst + point_base[k] * Dim, src.Data(), n * Dim * sizeof(double));
        }
        out.AssignPoints(std::move(points));
    }

    // --- point_data (every contour carries the same names, from one simp) ----
    for (const std::string& name : rContours[0].mMesh.PointDataNames()) {
        const NDArray& a0 = rContours[0].mMesh.PointData(name);
        const std::size_t rows0 = detail::rows(a0);
        const std::size_t comp = rows0 == 0 ? a0.Size() : a0.Size() / rows0;
        std::vector<std::size_t> shape = a0.Shape();
        if (shape.empty())
            shape = {total_points};
        else
            shape[0] = total_points;
        NDArray dst_arr = NDArray::Uninit(DType::Float64, std::move(shape));
        double* dst = dst_arr.As<double>();
        for (std::size_t k = 0; k < rContours.size(); ++k) {
            const NDArray& src = rContours[k].mMesh.PointData(name);
            const std::size_t n = rContours[k].mMesh.NumPoints();
            if (n != 0)
                std::memcpy(dst + point_base[k] * comp, src.Data(), n * comp * sizeof(double));
            // The contoured field reads back as *exactly* the isovalue: linear
            // interpolation only reaches it to within round-off, and the value
            // is known analytically. Not applied to a magnitude reduction --
            // |lerp(v)| != lerp(|v|) mathematically, so there is nothing exact
            // to write (see isosurface.hpp).
            if (name == rExactName && comp != 0) {
                const double v = rContours[k].mValue;
                if (rExactComponent.has_value()) {
                    const std::size_t c = static_cast<std::size_t>(*rExactComponent);
                    for (std::size_t i = 0; i < n; ++i)
                        dst[(point_base[k] + i) * comp + c] = v;
                } else if (comp == 1) {
                    for (std::size_t i = 0; i < n; ++i)
                        dst[point_base[k] + i] = v;
                }
            }
        }
        out.AddPointData(name, std::move(dst_arr));
    }

    // --- cells: merge by type, in the cutter's canonical order ---------------
    const std::string types[3] = {cell_type_name(CellType::Triangle),
                                  cell_type_name(CellType::Quad), cell_type_name(CellType::Line)};
    // Per output block: the (contour, source block) pairs feeding it.
    std::vector<std::string> out_types;
    std::vector<std::vector<std::pair<std::size_t, std::size_t>>> out_sources;
    for (const std::string& type : types) {
        std::vector<std::pair<std::size_t, std::size_t>> src;
        for (std::size_t k = 0; k < rContours.size(); ++k) {
            const std::optional<std::size_t> bi = iso_type_block(rContours[k].mMesh, type);
            if (bi.has_value())
                src.emplace_back(k, *bi);
        }
        if (!src.empty()) {
            out_types.push_back(type);
            out_sources.push_back(std::move(src));
        }
    }

    for (std::size_t ob = 0; ob < out_types.size(); ++ob) {
        std::size_t npc = 0;
        std::size_t ncells = 0;
        for (const auto& s : out_sources[ob]) {
            const auto cb = rContours[s.first].mMesh.Cells(s.second);
            npc = cb.NodesPerCell();
            ncells += cb.NumCells();
        }
        NDArray conn = NDArray::Uninit(DType::Int64, {ncells, npc});
        std::int64_t* dst = conn.As<std::int64_t>();
        std::size_t at = 0;
        for (const auto& s : out_sources[ob]) {
            const auto cb = rContours[s.first].mMesh.Cells(s.second);
            const NDArray& src = cb.Conn();
            const std::int64_t off = static_cast<std::int64_t>(point_base[s.first]);
            const std::size_t n = cb.NumCells() * npc;
            for (std::size_t i = 0; i < n; ++i)
                dst[at + i] = detail::read_int(src, i) + off;
            at += n;
        }
        out.AddCellBlock(out_types[ob], std::move(conn));
    }

    // --- cell_data: concatenate the parents' replicated rows -----------------
    for (const std::string& name : rContours[0].mMesh.CellDataNames()) {
        const NDArray& a0 = rContours[0].mMesh.CellData(name, 0);
        const DType dt = a0.Dtype();
        std::size_t comp = 1;
        for (std::size_t d = 1; d < a0.Shape().size(); ++d)
            comp *= a0.Shape()[d];
        const std::size_t row_bytes = comp * dtype_size(dt);
        std::vector<NDArray> blocks;
        blocks.reserve(out_types.size());
        for (std::size_t ob = 0; ob < out_types.size(); ++ob) {
            std::size_t ncells = 0;
            for (const auto& s : out_sources[ob])
                ncells += rContours[s.first].mMesh.Cells(s.second).NumCells();
            std::vector<std::size_t> shape = a0.Shape();
            if (shape.empty())
                shape = {ncells};
            else
                shape[0] = ncells;
            NDArray dst_arr = NDArray::Uninit(dt, std::move(shape));
            std::byte* dst = dst_arr.Data();
            std::size_t at = 0;
            for (const auto& s : out_sources[ob]) {
                const std::size_t n = rContours[s.first].mMesh.Cells(s.second).NumCells();
                const NDArray& src = rContours[s.first].mMesh.CellData(name, s.second);
                if (n != 0)
                    std::memcpy(dst + at * row_bytes, src.Data(), n * row_bytes);
                at += n;
            }
            blocks.push_back(std::move(dst_arr));
        }
        out.AddCellData(name, std::move(blocks));
    }

    // --- the contour tags ----------------------------------------------------
    // iso:value is the real level; iso:index is its ordinal, and is what
    // split(SplitBy::Tag) needs -- that criterion only accepts an integer array.
    {
        std::vector<NDArray> values;
        std::vector<NDArray> indices;
        values.reserve(out_types.size());
        indices.reserve(out_types.size());
        for (std::size_t ob = 0; ob < out_types.size(); ++ob) {
            std::size_t ncells = 0;
            for (const auto& s : out_sources[ob])
                ncells += rContours[s.first].mMesh.Cells(s.second).NumCells();
            NDArray v = NDArray::Uninit(DType::Float64, {ncells});
            NDArray i = NDArray::Uninit(DType::Int64, {ncells});
            double* vp = v.As<double>();
            std::int64_t* ip = i.As<std::int64_t>();
            std::size_t at = 0;
            for (const auto& s : out_sources[ob]) {
                const std::size_t n = rContours[s.first].mMesh.Cells(s.second).NumCells();
                for (std::size_t c = 0; c < n; ++c) {
                    vp[at + c] = rContours[s.first].mValue;
                    ip[at + c] = rContours[s.first].mIndex;
                }
                at += n;
            }
            values.push_back(std::move(v));
            indices.push_back(std::move(i));
        }
        out.AddCellData("iso:value", std::move(values));
        out.AddCellData("iso:index", std::move(indices));
    }

    return out;
}

}  // namespace

Mesh isosurface(const Mesh& rMesh, const IsosurfaceOptions& rOptions) {
    if (rOptions.mArrayName.empty())
        throw std::invalid_argument("meshio++: isosurface: an array name is required");
    if (rOptions.mIsovalues.empty())
        throw std::invalid_argument("meshio++: isosurface: at least one isovalue is required");
    for (double v : rOptions.mIsovalues)
        if (!std::isfinite(v))
            throw std::invalid_argument("meshio++: isosurface: isovalues must be finite");

    // A cell field is piecewise constant: there is no crossing to locate inside
    // a cell, so there is no level set to draw. Say so by name rather than
    // reporting "unknown array", and name the conversion that fixes it.
    if (!rMesh.HasPointData(rOptions.mArrayName)) {
        if (rMesh.HasCellData(rOptions.mArrayName))
            throw std::invalid_argument(
                "meshio++: isosurface: '" + rOptions.mArrayName +
                "' is cell_data, which is piecewise constant and has no level set; convert it "
                "first with cell_data_to_point_data (CLI: meshioplusplus data to-point)");
        throw std::invalid_argument(
            data_unknown_key_message(rMesh, DataLocation::Point, rOptions.mArrayName));
    }

    // Ascending, de-duplicated: contours are emitted in ascending order, and a
    // repeated isovalue would otherwise emit the same contour twice.
    std::vector<double> isovalues = rOptions.mIsovalues;
    std::sort(isovalues.begin(), isovalues.end());
    isovalues.erase(std::unique(isovalues.begin(), isovalues.end()), isovalues.end());

    // Simplexify once; every isovalue is cut against the same prepared input.
    const detail::MarchingInput input = detail::marching_prepare(rMesh);
    const Mesh& simp = input.Simp();
    const std::size_t nnodes = simp.NumPoints();
    const std::size_t dim = simp.PointDim();

    // The field rides through simplexify (points are untouched for a linear
    // input; a higher-order one takes the linearization's row subset), so it is
    // resolved against the *simplexified* mesh -- interpolate's rule.
    std::vector<double> field(nnodes, 0.0);
    {
        const NDArray& a = simp.PointData(rOptions.mArrayName);
        const std::size_t src_rows = detail::rows(a);
        const std::size_t comp = src_rows == 0 ? 0 : a.Size() / src_rows;
        if (comp == 0)
            throw std::invalid_argument("meshio++: isosurface: '" + rOptions.mArrayName +
                                        "' carries no values");
        // Validate the component once, outside the loop, so a bad request fails
        // by name rather than n times from inside a parallel region.
        iso_scalarize(a, 0, comp, rOptions.mComponent);
        parallel_for(nnodes, [&](std::size_t g) {
            field[g] = iso_scalarize(a, g, comp, rOptions.mComponent);
        });
    }

    detail::MarchingOptions cut_options;
    cut_options.mProvenanceName = "iso:parent_cell";
    cut_options.mRecordParentIds = rOptions.mRecordParentIds;
    // A level set has no global direction, so faces are wound toward increasing
    // field (slice, which does have one, uses FixedDirection).
    cut_options.mOrientation = detail::MarchingOrientation::FieldGradient;

    std::vector<IsoContour> contours;
    contours.reserve(isovalues.size());
    Mesh empty_like;
    bool have_empty_like = false;
    std::vector<double> dist(nnodes, 0.0);
    for (std::size_t k = 0; k < isovalues.size(); ++k) {
        const double iso = isovalues[k];
        parallel_for(nnodes, [&](std::size_t g) { dist[g] = field[g] - iso; });
        Mesh section = detail::marching_cut(input, dist, cut_options);
        if (section.NumCellBlocks() == 0) {
            // Out of range (or fully collapsed): no contour at this level. Keep
            // the first such result so an all-empty request still returns a mesh
            // carrying the (zero-row) point_data names, exactly like slice.
            if (!have_empty_like) {
                empty_like = std::move(section);
                have_empty_like = true;
            }
            continue;
        }
        contours.push_back(IsoContour{std::move(section), iso, static_cast<std::int64_t>(k)});
    }

    // A level set is newly created geometry: no input point, cell or facet has
    // a counterpart in it, so no named region can be carried across.
    detail::warn_regions_dropped(rMesh, "isosurface");

    if (contours.empty())
        return have_empty_like ? std::move(empty_like) : Mesh{};

    // The contoured array is written exactly, except under a magnitude
    // reduction, where no exact value exists (see isosurface.hpp).
    const std::size_t comp = [&]() {
        const NDArray& a = simp.PointData(rOptions.mArrayName);
        const std::size_t r = detail::rows(a);
        return r == 0 ? std::size_t{0} : a.Size() / r;
    }();
    const bool exact = rOptions.mComponent.has_value() || comp == 1;
    return iso_concat_contours(contours, dim, exact ? rOptions.mArrayName : std::string(),
                               rOptions.mComponent);
}

}  // namespace meshioplusplus
