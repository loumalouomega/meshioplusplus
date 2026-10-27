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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/error_size.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kEsPrefix = "meshio++: error_size: ";

[[noreturn]] void es_fail(const std::string& rWhy) {
    throw std::invalid_argument(std::string(kEsPrefix) + rWhy);
}

}  // namespace

ErrorSizeResult error_size(const Mesh& rMesh, const ErrorSizeOptions& rOptions) {
    if (!(std::isfinite(rOptions.mTargetError) && rOptions.mTargetError > 0.0))
        es_fail("target_error must be finite and positive");
    if (!(std::isfinite(rOptions.mOrder) && rOptions.mOrder > 0.0))
        es_fail("order must be finite and positive");
    if (!(std::isfinite(rOptions.mMinSize) && rOptions.mMinSize >= 0.0) ||
        !(std::isfinite(rOptions.mMaxSize) && rOptions.mMaxSize >= 0.0))
        es_fail("min_size and max_size must be finite and not negative");
    if (rOptions.mMinSize > 0.0 && rOptions.mMaxSize > 0.0 && rOptions.mMinSize > rOptions.mMaxSize)
        es_fail("min_size must not exceed max_size");
    if (rOptions.mOutputName.empty())
        es_fail("an output name is required");
    if (!rMesh.HasCellData(rOptions.mErrorArray)) {
        if (rMesh.HasPointData(rOptions.mErrorArray))
            es_fail("'" + rOptions.mErrorArray +
                    "' is point_data; the indicator is per cell (estimate_error writes error:zz)");
        throw std::invalid_argument(
            data_unknown_key_message(rMesh, DataLocation::Cell, rOptions.mErrorArray));
    }
    const std::size_t nb = rMesh.NumCellBlocks();
    for (std::size_t b = 0; b < nb; ++b)
        if (data_num_components(rMesh.CellData(rOptions.mErrorArray, b)) != 1)
            es_fail("'" + rOptions.mErrorArray + "' must have one component");
    if (!rOptions.mOverwrite && rMesh.HasPointData(rOptions.mOutputName))
        es_fail("'" + rOptions.mOutputName +
                "' already exists in point_data (pass overwrite=true to replace it)");

    const std::size_t np = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();
    const NDArray& r_pts = rMesh.Points();
    // Each cell's diameter and indicator, block by block, in block order.
    std::vector<std::size_t> first(nb + 1, 0);
    for (std::size_t b = 0; b < nb; ++b)
        first[b + 1] = first[b] + rMesh.Cells(b).NumCells();
    const std::size_t nc = first[nb];
    std::vector<double> diam(nc, 0.0), eta(nc, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t b = 0; b < nb; ++b) {
        const auto cells = rMesh.Cells(b);
        const NDArray& r_e = rMesh.CellData(rOptions.mErrorArray, b);
        parallel_for(cells.NumCells(), [&](std::size_t c) {
            const std::size_t k = first[b] + c;
            eta[k] = detail::read_double(r_e, c);
            std::vector<std::int64_t> ids;
            if (cells.IsPolyhedron()) {
                for (std::size_t f = 0; f < cells.NumFaces(c); ++f) {
                    const auto face = cells.Face(c, f);
                    ids.insert(ids.end(), face.first, face.first + face.second);
                }
            } else if (cells.IsRagged()) {
                const std::int64_t* p_row = cells.Row(c);
                ids.assign(p_row, p_row + cells.RowSize(c));
            } else {
                const std::size_t n = cells.NodesPerCell();
                for (std::size_t i = 0; i < n; ++i)
                    ids.push_back(detail::read_int(cells.Conn(), n * c + i));
            }
            double d2 = 0.0;
            for (std::size_t i = 0; i < ids.size(); ++i)
                for (std::size_t j = i + 1; j < ids.size(); ++j) {
                    double s = 0.0;
                    for (std::size_t a = 0; a < dim; ++a) {
                        const double t =
                            detail::read_double(r_pts, dim * static_cast<std::size_t>(ids[i]) + a) -
                            detail::read_double(r_pts, dim * static_cast<std::size_t>(ids[j]) + a);
                        s += t * t;
                    }
                    d2 = std::max(d2, s);
                }
            diam[k] = std::sqrt(d2);
        });
    }
    ErrorSizeResult result;
    std::vector<char> usable(nc, 0);
    for (std::size_t k = 0; k < nc; ++k)
        if (std::isfinite(eta[k]) && eta[k] >= 0.0 && diam[k] > 0.0) {
            usable[k] = 1;
            ++result.mNumCells;
        } else {
            ++result.mNumSkipped;
        }
    // The evenly spread share (Zienkiewicz and Zhu 1987), and each cell's
    // size: a cell with no error asks for the greatest size allowed.
    const double share =
        result.mNumCells > 0 ? rOptions.mTargetError / std::sqrt(double(result.mNumCells)) : 0.0;
    std::vector<double> size(nc, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t k = 0; k < nc; ++k) {
        if (!usable[k])
            continue;
        double h = eta[k] > 0.0 ? diam[k] * std::pow(share / eta[k], 1.0 / rOptions.mOrder)
                                : std::numeric_limits<double>::infinity();
        const double lo = rOptions.mMinSize, hi = rOptions.mMaxSize;
        if (lo > 0.0 && h < lo) {
            h = lo;
            ++result.mNumClamped;
        } else if (hi > 0.0 && h > hi) {
            h = hi;
            ++result.mNumClamped;
        }
        if (!std::isfinite(h)) {
            --result.mNumCells;
            ++result.mNumSkipped;
            continue;
        }
        size[k] = h;
    }
    // Per point, the least of its cells' (block order, cell order: the
    // minimum does not depend on it).
    std::vector<double> at(np, std::numeric_limits<double>::infinity());
    for (std::size_t b = 0; b < nb; ++b) {
        const auto cells = rMesh.Cells(b);
        for (std::size_t c = 0; c < cells.NumCells(); ++c) {
            const double h = size[first[b] + c];
            if (!std::isfinite(h))
                continue;
            const auto note = [&](std::int64_t P) {
                double& r_a = at[static_cast<std::size_t>(P)];
                r_a = std::min(r_a, h);
            };
            if (cells.IsPolyhedron()) {
                for (std::size_t f = 0; f < cells.NumFaces(c); ++f) {
                    const auto face = cells.Face(c, f);
                    for (std::size_t i = 0; i < face.second; ++i)
                        note(face.first[i]);
                }
            } else if (cells.IsRagged()) {
                const std::int64_t* p_row = cells.Row(c);
                for (std::size_t i = 0; i < cells.RowSize(c); ++i)
                    note(p_row[i]);
            } else {
                const std::size_t n = cells.NodesPerCell();
                for (std::size_t i = 0; i < n; ++i)
                    note(detail::read_int(cells.Conn(), n * c + i));
            }
        }
    }
    NDArray out = NDArray::Uninit(DType::Float64, {np});
    double* po = out.As<double>();
    result.mMinSize = std::numeric_limits<double>::infinity();
    result.mMaxSize = 0.0;
    for (std::size_t i = 0; i < np; ++i) {
        po[i] = std::isfinite(at[i]) ? at[i] : std::numeric_limits<double>::quiet_NaN();
        if (std::isfinite(at[i])) {
            result.mMinSize = std::min(result.mMinSize, at[i]);
            result.mMaxSize = std::max(result.mMaxSize, at[i]);
        }
    }
    if (!std::isfinite(result.mMinSize))
        result.mMinSize = 0.0;
    result.mMesh = detail::clone_mesh(
        rMesh, [&](DataLocation Location, const std::string& rName, std::string&) {
            return !(Location == DataLocation::Point && rName == rOptions.mOutputName);
        });
    result.mMesh.AddPointData(rOptions.mOutputName, std::move(out));
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        result.mMesh.AddRegion(rMesh.Region(i));
    for (std::size_t i = 0; i < rMesh.NumPropertySets(); ++i)
        result.mMesh.AddPropertySet(rMesh.GetPropertySet(i));
    return result;
}

}  // namespace meshioplusplus
