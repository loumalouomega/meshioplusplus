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
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/colorize.hpp"
#include "meshioplusplus/detail/colormap.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/face_color.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/ndarray.hpp"
// Project includes (private, not installed)
#include "../detail/hex_color.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kColPrefix = "meshio++: colorize: ";

// One scalar per row: a component, or the magnitude of the row.
double col_scalar(const NDArray& rA, std::size_t Row, std::size_t NComp, std::int64_t Component) {
    if (NComp <= 1)
        return detail::read_double(rA, Row);
    if (Component >= 0)
        return detail::read_double(rA, Row * NComp + static_cast<std::size_t>(Component));
    double sum = 0.0;
    for (std::size_t k = 0; k < NComp; ++k) {
        const double v = detail::read_double(rA, Row * NComp + k);
        sum += v * v;
    }
    return std::sqrt(sum);
}

// The values of one array (a point array, or one cell block), scalarized.
void col_values(const NDArray& rA, std::size_t Rows, std::int64_t Component,
                std::vector<double>& rOut) {
    const std::size_t ncomp = Rows == 0 ? 1 : rA.Size() / Rows;
    if (Component >= 0 && static_cast<std::size_t>(Component) >= ncomp && ncomp > 1)
        throw std::invalid_argument(std::string(kColPrefix) + "component " +
                                    std::to_string(Component) + " is out of range (the field has " +
                                    std::to_string(ncomp) + ")");
    for (std::size_t r = 0; r < Rows; ++r)
        rOut.push_back(col_scalar(rA, r, ncomp, Component));
}

}  // namespace

ColorizeResult colorize(const Mesh& rMesh, const ColorizeOptions& rOptions) {
    if (rOptions.mField.empty())
        throw std::invalid_argument(std::string(kColPrefix) + "name the field to colour by");
    if (rOptions.mSteps < 0)
        throw std::invalid_argument(std::string(kColPrefix) + "steps must be >= 0");
    if (!std::isfinite(rOptions.mVMin) || !std::isfinite(rOptions.mVMax))
        throw std::invalid_argument(std::string(kColPrefix) + "vmin and vmax must be finite");
    const std::uint8_t* table = detail::colormap_table(rOptions.mCmap);
    const detail::Rgb nan_color = detail::parse_hex_color(rOptions.mNanColor, kColPrefix);
    const bool point = rMesh.HasPointData(rOptions.mField);
    if (!point && !rMesh.HasCellData(rOptions.mField))
        throw std::invalid_argument(std::string(kColPrefix) + "no point or cell data named '" +
                                    rOptions.mField + "'");

    // Every value, point order or block-major cell order.
    std::vector<double> values;
    std::vector<std::size_t> block_rows;
    if (point) {
        col_values(rMesh.PointData(rOptions.mField), rMesh.NumPoints(), rOptions.mComponent,
                   values);
    } else {
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const std::size_t nc = rMesh.Cells(b).NumCells();
            col_values(rMesh.CellData(rOptions.mField, b), nc, rOptions.mComponent, values);
            block_rows.push_back(nc);
        }
    }

    ColorizeResult result;
    double lo = rOptions.mVMin, hi = rOptions.mVMax;
    if (!(lo < hi)) {
        bool seen = false;
        double amax = 0.0;
        for (double v : values) {
            if (!std::isfinite(v))
                continue;
            lo = seen ? std::min(lo, v) : v;
            hi = seen ? std::max(hi, v) : v;
            amax = std::max(amax, std::fabs(v));
            seen = true;
        }
        if (!seen)
            lo = hi = 0.0;
        if (rOptions.mSymmetric) {
            lo = -amax;
            hi = amax;
        }
    }
    result.mVMin = lo;
    result.mVMax = hi;

    const std::size_t n = values.size();
    NDArray red = NDArray::Uninit(DType::UInt8, {n});
    NDArray green = NDArray::Uninit(DType::UInt8, {n});
    NDArray blue = NDArray::Uninit(DType::UInt8, {n});
    const std::int64_t steps = rOptions.mSteps;
    for (std::size_t i = 0; i < n; ++i) {
        detail::Rgb c = nan_color;
        if (std::isfinite(values[i])) {
            double t = detail::color_param(values[i], lo, hi);
            if (steps >= 2) {
                const double s = static_cast<double>(steps);
                double band = std::floor(t * s);
                band = band >= s ? s - 1.0 : band;
                t = (band + 0.5) / s;
            }
            c = detail::colormap_lookup(table, t);
        }
        red.As<std::uint8_t>()[i] = c.mR;
        green.As<std::uint8_t>()[i] = c.mG;
        blue.As<std::uint8_t>()[i] = c.mB;
    }

    result.mMesh = detail::clone_mesh(rMesh);
    if (point) {
        result.mLocation = "point";
        result.mMesh.AddPointData(kColorRedName, std::move(red));
        result.mMesh.AddPointData(kColorGreenName, std::move(green));
        result.mMesh.AddPointData(kColorBlueName, std::move(blue));
        return result;
    }
    result.mLocation = "cell";
    const NDArray* channels[3] = {&red, &green, &blue};
    const char* names[3] = {kColorRedName, kColorGreenName, kColorBlueName};
    for (int ch = 0; ch < 3; ++ch) {
        std::vector<NDArray> blocks;
        std::size_t offset = 0;
        for (std::size_t nc : block_rows) {
            NDArray a = NDArray::Uninit(DType::UInt8, {nc});
            for (std::size_t r = 0; r < nc; ++r)
                a.As<std::uint8_t>()[r] = channels[ch]->As<std::uint8_t>()[offset + r];
            offset += nc;
            blocks.push_back(std::move(a));
        }
        result.mMesh.AddCellData(names[ch], std::move(blocks));
    }
    return result;
}

}  // namespace meshioplusplus
