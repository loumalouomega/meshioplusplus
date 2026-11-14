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
// Bake a point or cell array into a texture through `uv`. See
// operations/bake_texture.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/operations/bake_texture.hpp"
#include "meshioplusplus/detail/colormap.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/zlib_inflate.hpp"
#include "../detail/png_write.hpp"

namespace meshioplusplus {

namespace {

/// One value per row of @p rArray: the component, or the magnitude.
std::vector<double> bake_values(const NDArray& rArray, std::size_t Rows, std::int64_t Component) {
    const std::size_t width = Rows ? rArray.Size() / Rows : 0;
    if (width == 0 || rArray.Size() != Rows * width)
        throw std::invalid_argument(
            "meshio++: bake_texture: the array does not hold one row per "
            "entry");
    if (Component >= static_cast<std::int64_t>(width))
        throw std::invalid_argument("meshio++: bake_texture: component " +
                                    std::to_string(Component) + " is out of range");
    std::vector<double> out(Rows);
    for (std::size_t r = 0; r < Rows; ++r) {
        if (width == 1) {
            out[r] = detail::read_double(rArray, r);
        } else if (Component >= 0) {
            out[r] = detail::read_double(rArray, r * width + static_cast<std::size_t>(Component));
        } else {
            double s = 0.0;
            for (std::size_t c = 0; c < width; ++c) {
                const double v = detail::read_double(rArray, r * width + c);
                s += v * v;
            }
            out[r] = std::sqrt(s);
        }
    }
    return out;
}

}  // namespace

BakeTextureResult bake_texture(const Mesh& rMesh, const BakeTextureOptions& rOptions) {
    if (rOptions.mWidth < 1 || rOptions.mHeight < 1 || rOptions.mWidth > 16384 ||
        rOptions.mHeight > 16384)
        throw std::invalid_argument(
            "meshio++: bake_texture: width and height must lie in [1, 16384]");
    const std::size_t np = rMesh.NumPoints();
    if (!rMesh.HasPointData("uv") || rMesh.PointData("uv").Size() != 2 * np)
        throw std::invalid_argument(
            "meshio++: bake_texture: the mesh has no point_data['uv'] "
            "(two per point); run parameterize first");
    const NDArray& uv = rMesh.PointData("uv");
    const std::uint8_t* table = detail::colormap_table(rOptions.mCmap);
    const bool on_points = rMesh.HasPointData(rOptions.mArray);
    std::size_t total_cells = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
        total_cells += rMesh.Cells(b).NumCells();
    std::vector<double> values;
    if (on_points) {
        values = bake_values(rMesh.PointData(rOptions.mArray), np, rOptions.mComponent);
    } else if (rMesh.HasCellData(rOptions.mArray) &&
               rMesh.CellDataNumBlocks(rOptions.mArray) == rMesh.NumCellBlocks()) {
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const auto v = bake_values(rMesh.CellData(rOptions.mArray, b),
                                       rMesh.Cells(b).NumCells(), rOptions.mComponent);
            values.insert(values.end(), v.begin(), v.end());
        }
    } else {
        throw std::invalid_argument("meshio++: bake_texture: no point or cell data array '" +
                                    rOptions.mArray + "'");
    }
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (const double v : values)
        if (std::isfinite(v)) {
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
    if (!std::isfinite(lo))
        lo = hi = 0.0;
    BakeTextureResult out;
    out.mVMin = rOptions.mVMin.value_or(lo);
    out.mVMax = rOptions.mVMax.value_or(hi);
    const std::size_t w = static_cast<std::size_t>(rOptions.mWidth);
    const std::size_t h = static_cast<std::size_t>(rOptions.mHeight);
    std::vector<std::uint8_t> rgba(w * h * 4, 0);
    const double W = static_cast<double>(w), H = static_cast<double>(h);
    const auto draw = [&](std::int64_t A, std::int64_t B, std::int64_t C, std::size_t Cell) {
        const double xa = detail::read_double(uv, 2 * static_cast<std::size_t>(A)) * W;
        const double ya = (1.0 - detail::read_double(uv, 2 * static_cast<std::size_t>(A) + 1)) * H;
        const double xb = detail::read_double(uv, 2 * static_cast<std::size_t>(B)) * W;
        const double yb = (1.0 - detail::read_double(uv, 2 * static_cast<std::size_t>(B) + 1)) * H;
        const double xc = detail::read_double(uv, 2 * static_cast<std::size_t>(C)) * W;
        const double yc = (1.0 - detail::read_double(uv, 2 * static_cast<std::size_t>(C) + 1)) * H;
        const double area = (xb - xa) * (yc - ya) - (xc - xa) * (yb - ya);
        if (!(area != 0.0) || !std::isfinite(area))
            return;
        const double x_lo = std::max(0.0, std::ceil(std::min({xa, xb, xc}) - 0.5));
        const double x_hi = std::min(W - 1.0, std::floor(std::max({xa, xb, xc}) - 0.5));
        const double y_lo = std::max(0.0, std::ceil(std::min({ya, yb, yc}) - 0.5));
        const double y_hi = std::min(H - 1.0, std::floor(std::max({ya, yb, yc}) - 0.5));
        if (!(x_lo <= x_hi) || !(y_lo <= y_hi))
            return;
        for (std::size_t py = static_cast<std::size_t>(y_lo); py <= static_cast<std::size_t>(y_hi);
             ++py) {
            const double cy = static_cast<double>(py) + 0.5;
            for (std::size_t px = static_cast<std::size_t>(x_lo);
                 px <= static_cast<std::size_t>(x_hi); ++px) {
                const double cx = static_cast<double>(px) + 0.5;
                const double w0 = ((xb - cx) * (yc - cy) - (xc - cx) * (yb - cy)) / area;
                const double w1 = ((xc - cx) * (ya - cy) - (xa - cx) * (yc - cy)) / area;
                const double w2 = ((xa - cx) * (yb - cy) - (xb - cx) * (ya - cy)) / area;
                if (!(w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0))
                    continue;
                const double value = on_points ? w0 * values[static_cast<std::size_t>(A)] +
                                                     w1 * values[static_cast<std::size_t>(B)] +
                                                     w2 * values[static_cast<std::size_t>(C)]
                                               : values[Cell];
                std::uint8_t* p = rgba.data() + (py * w + px) * 4;
                if (!std::isfinite(value)) {
                    p[0] = p[1] = p[2] = p[3] = 0;
                    continue;
                }
                const double t =
                    out.mVMax > out.mVMin ? (value - out.mVMin) / (out.mVMax - out.mVMin) : 0.5;
                const detail::Rgb c = detail::colormap_lookup(table, t);
                p[0] = c.mR;
                p[1] = c.mG;
                p[2] = c.mB;
                p[3] = 255;
            }
        }
    };
    std::size_t cell = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type = cb.Type();
        const bool face = type == "triangle" || type == "quad" || type == "polygon";
        for (std::size_t r = 0; r < cb.NumCells(); ++r, ++cell) {
            if (!face)
                continue;
            const std::size_t n = cb.IsRagged() ? cb.RowSize(r) : cb.NodesPerCell();
            const auto node = [&](std::size_t k) -> std::int64_t {
                return cb.IsRagged() ? cb.Row(r)[k] : detail::read_int(cb.Conn(), r * n + k);
            };
            for (std::size_t k = 1; k + 1 < n; ++k)
                draw(node(0), node(k), node(k + 1), cell);
        }
    }
    for (std::size_t k = 3; k < rgba.size(); k += 4)
        out.mCoveredPixels += rgba[k] != 0;
    out.mImage = NDArray(DType::UInt8, {h, w, 4});
    std::memcpy(out.mImage.Data(), rgba.data(), rgba.size());
    const std::string png =
        detail::png_encode_rgba(rgba.data(), w, h, detail::zlib_available() ? 9 : 0);
    out.mMesh = detail::clone_mesh(rMesh);
    NDArray bytes(DType::UInt8, {png.size()});
    std::memcpy(bytes.Data(), png.data(), png.size());
    out.mMesh.AddFieldData("texture:png", std::move(bytes));
    if (out.mCoveredPixels == 0)
        out.mWarnings.push_back("bake_texture: no 2-D cell covers a pixel of the texture");
    return out;
}

}  // namespace meshioplusplus
