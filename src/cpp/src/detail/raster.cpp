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

// The tile rasterizer. See detail/raster.hpp for the determinism contract.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

// Project includes
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "raster.hpp"

namespace meshioplusplus {
namespace detail {
namespace {

using RstColor = std::array<std::uint8_t, 4>;

constexpr std::int64_t kRstOne = std::int64_t{1} << kRasterSubBits;
constexpr std::int64_t kRstHalf = kRstOne / 2;

// A screen coordinate in fixed point, clamped so every edge-function product
// fits in 64 bits (|coordinate| <= 2^26, differences <= 2^27, products <= 2^54).
std::int64_t rst_fixed(double v) {
    if (!(v == v))
        v = 0.0;
    if (v > kRasterLimitPx)
        v = kRasterLimitPx;
    if (v < -kRasterLimitPx)
        v = -kRasterLimitPx;
    return std::llround(v * static_cast<double>(kRstOne));
}

std::int64_t rst_floor_div(std::int64_t a, std::int64_t b) {
    std::int64_t q = a / b;
    if ((a % b) != 0 && a < 0)
        --q;
    return q;
}

// The first and last pixel whose centre (p * one + half) lies in [lo, hi].
std::int64_t rst_first_px(std::int64_t Lo) {
    return -rst_floor_div(-(Lo - kRstHalf), kRstOne);
}
std::int64_t rst_last_px(std::int64_t Hi) {
    return rst_floor_div(Hi - kRstHalf, kRstOne);
}

std::uint8_t rst_scale(std::uint8_t c, double Intensity) {
    const double v = static_cast<double>(c) * Intensity + 0.5;
    if (!(v > 0.0))
        return 0;
    if (v >= 255.0)
        return 255;
    return static_cast<std::uint8_t>(v);
}

// The top-left rule: an edge a->b of a positively oriented triangle (interior
// where the edge function grows; y points down) is a left edge when it runs
// upward, and a top edge when it is horizontal and runs right.
bool rst_top_left(std::int64_t Dx, std::int64_t Dy) {
    return Dy < 0 || (Dy == 0 && Dx > 0);
}

struct RstTri {
    std::int64_t mX[3];
    std::int64_t mY[3];
    double mZ[3];
    double mI[3];
    std::int64_t mArea = 0;
    std::int64_t mBias[3];  // 0 on a top-left edge, -1 otherwise
    std::int64_t mPx0 = 0, mPx1 = -1, mPy0 = 0, mPy1 = -1;
    RstColor mColor;
    std::int64_t mId = -1;
};

struct RstSeg {
    double mAx, mAy, mAz, mBx, mBy, mBz;
    std::int64_t mSteps = 1;
    std::int64_t mPx0 = 0, mPx1 = -1, mPy0 = 0, mPy1 = -1;
};

// Edge k of a triangle runs from corner k+1 to corner k+2; its function is
// the barycentric weight of corner k scaled by the doubled area.
void rst_prepare_tri(const RasterScene& rScene, const RasterTri& rTri, int Width, int Height,
                     RstTri& rOut) {
    std::int64_t v[3] = {rTri.mV[0], rTri.mV[1], rTri.mV[2]};
    double intensity[3] = {rTri.mIntensity[0], rTri.mIntensity[1], rTri.mIntensity[2]};
    for (int k = 0; k < 3; ++k) {
        const RasterVertex& p = rScene.mVertices[static_cast<std::size_t>(v[k])];
        rOut.mX[k] = rst_fixed(p.mX);
        rOut.mY[k] = rst_fixed(p.mY);
        rOut.mZ[k] = p.mDepth;
        rOut.mI[k] = intensity[k];
    }
    std::int64_t area = (rOut.mX[1] - rOut.mX[0]) * (rOut.mY[2] - rOut.mY[0]) -
                        (rOut.mY[1] - rOut.mY[0]) * (rOut.mX[2] - rOut.mX[0]);
    if (area == 0)
        return;  // degenerate: covers nothing
    if (area < 0) {
        std::swap(rOut.mX[1], rOut.mX[2]);
        std::swap(rOut.mY[1], rOut.mY[2]);
        std::swap(rOut.mZ[1], rOut.mZ[2]);
        std::swap(rOut.mI[1], rOut.mI[2]);
        area = -area;
    }
    rOut.mArea = area;
    for (int k = 0; k < 3; ++k) {
        const int a = (k + 1) % 3;
        const int b = (k + 2) % 3;
        rOut.mBias[k] = rst_top_left(rOut.mX[b] - rOut.mX[a], rOut.mY[b] - rOut.mY[a]) ? 0 : -1;
    }
    const std::int64_t min_x = std::min({rOut.mX[0], rOut.mX[1], rOut.mX[2]});
    const std::int64_t max_x = std::max({rOut.mX[0], rOut.mX[1], rOut.mX[2]});
    const std::int64_t min_y = std::min({rOut.mY[0], rOut.mY[1], rOut.mY[2]});
    const std::int64_t max_y = std::max({rOut.mY[0], rOut.mY[1], rOut.mY[2]});
    rOut.mPx0 = std::max<std::int64_t>(rst_first_px(min_x), 0);
    rOut.mPx1 = std::min<std::int64_t>(rst_last_px(max_x), Width - 1);
    rOut.mPy0 = std::max<std::int64_t>(rst_first_px(min_y), 0);
    rOut.mPy1 = std::min<std::int64_t>(rst_last_px(max_y), Height - 1);
    rOut.mColor = rTri.mColor;
    rOut.mId = rTri.mId;
}

void rst_prepare_seg(const RasterScene& rScene, const RasterLine& rLine, int Width, int Height,
                     RstSeg& rOut) {
    const RasterVertex& a = rScene.mVertices[static_cast<std::size_t>(rLine.mA)];
    const RasterVertex& b = rScene.mVertices[static_cast<std::size_t>(rLine.mB)];
    // The same clamp as a triangle's corners, in pixels rather than fixed point.
    auto clamp = [](double v) {
        if (!(v == v))
            return 0.0;
        return std::min(std::max(v, -kRasterLimitPx), kRasterLimitPx);
    };
    rOut.mAx = clamp(a.mX);
    rOut.mAy = clamp(a.mY);
    rOut.mAz = a.mDepth;
    rOut.mBx = clamp(b.mX);
    rOut.mBy = clamp(b.mY);
    rOut.mBz = b.mDepth;
    const double span = std::max(std::fabs(rOut.mBx - rOut.mAx), std::fabs(rOut.mBy - rOut.mAy));
    rOut.mSteps = std::max<std::int64_t>(1, static_cast<std::int64_t>(std::ceil(span)));
    const std::int64_t w = rScene.mLineWidth;
    rOut.mPx0 = std::max<std::int64_t>(
        static_cast<std::int64_t>(std::floor(std::min(rOut.mAx, rOut.mBx))) - w, 0);
    rOut.mPx1 = std::min<std::int64_t>(
        static_cast<std::int64_t>(std::floor(std::max(rOut.mAx, rOut.mBx))) + w, Width - 1);
    rOut.mPy0 = std::max<std::int64_t>(
        static_cast<std::int64_t>(std::floor(std::min(rOut.mAy, rOut.mBy))) - w, 0);
    rOut.mPy1 = std::min<std::int64_t>(
        static_cast<std::int64_t>(std::floor(std::max(rOut.mAy, rOut.mBy))) + w, Height - 1);
}

struct RstTile {
    std::int64_t mX0, mX1, mY0, mY1;  // inclusive pixel range
};

void rst_draw_tri(const RstTri& rT, const RstTile& rTile, RasterTarget& rOut) {
    const std::int64_t px0 = std::max(rT.mPx0, rTile.mX0);
    const std::int64_t px1 = std::min(rT.mPx1, rTile.mX1);
    const std::int64_t py0 = std::max(rT.mPy0, rTile.mY0);
    const std::int64_t py1 = std::min(rT.mPy1, rTile.mY1);
    if (px0 > px1 || py0 > py1)
        return;
    std::int64_t dx[3], dy[3];
    for (int k = 0; k < 3; ++k) {
        const int a = (k + 1) % 3;
        const int b = (k + 2) % 3;
        dx[k] = rT.mX[b] - rT.mX[a];
        dy[k] = rT.mY[b] - rT.mY[a];
    }
    const double area = static_cast<double>(rT.mArea);
    const std::size_t width = static_cast<std::size_t>(rOut.mWidth);
    const std::int64_t start_x = px0 * kRstOne + kRstHalf;
    for (std::int64_t py = py0; py <= py1; ++py) {
        const std::int64_t sy = py * kRstOne + kRstHalf;
        std::int64_t e[3];
        for (int k = 0; k < 3; ++k) {
            const int a = (k + 1) % 3;
            e[k] = dx[k] * (sy - rT.mY[a]) - dy[k] * (start_x - rT.mX[a]);
        }
        for (std::int64_t px = px0; px <= px1; ++px) {
            if (e[0] + rT.mBias[0] >= 0 && e[1] + rT.mBias[1] >= 0 && e[2] + rT.mBias[2] >= 0) {
                const double w0 = static_cast<double>(e[0]);
                const double w1 = static_cast<double>(e[1]);
                const double w2 = static_cast<double>(e[2]);
                const double z = (w0 * rT.mZ[0] + w1 * rT.mZ[1] + w2 * rT.mZ[2]) / area;
                const float zf = static_cast<float>(z);
                const std::size_t idx =
                    static_cast<std::size_t>(py) * width + static_cast<std::size_t>(px);
                if (zf > rOut.mDepth[idx]) {
                    const double i = (w0 * rT.mI[0] + w1 * rT.mI[1] + w2 * rT.mI[2]) / area;
                    rOut.mDepth[idx] = zf;
                    std::uint8_t* pc = rOut.mRgba.data() + idx * 4;
                    pc[0] = rst_scale(rT.mColor[0], i);
                    pc[1] = rst_scale(rT.mColor[1], i);
                    pc[2] = rst_scale(rT.mColor[2], i);
                    pc[3] = rT.mColor[3];
                    rOut.mIds[idx] = rT.mId;
                }
            }
            for (int k = 0; k < 3; ++k)
                e[k] -= dy[k] * kRstOne;
        }
    }
}

void rst_plot(std::int64_t Px, std::int64_t Py, float Z, const RstColor& rColor, std::int64_t Id,
              const RstTile& rTile, RasterTarget& rOut) {
    if (Px < rTile.mX0 || Px > rTile.mX1 || Py < rTile.mY0 || Py > rTile.mY1)
        return;
    const std::size_t idx = static_cast<std::size_t>(Py) * static_cast<std::size_t>(rOut.mWidth) +
                            static_cast<std::size_t>(Px);
    if (Z > rOut.mDepth[idx]) {
        rOut.mDepth[idx] = Z;
        std::copy(rColor.begin(), rColor.end(), rOut.mRgba.data() + idx * 4);
        rOut.mIds[idx] = Id;
    }
}

void rst_draw_seg(const RstSeg& rS, const RasterLine& rLine, const RasterScene& rScene,
                  const RstTile& rTile, RasterTarget& rOut) {
    if (std::max(rS.mPx0, rTile.mX0) > std::min(rS.mPx1, rTile.mX1) ||
        std::max(rS.mPy0, rTile.mY0) > std::min(rS.mPy1, rTile.mY1))
        return;
    const std::int64_t w = rScene.mLineWidth;
    // Restrict the steps to the part of the segment near this tile (Liang-Barsky
    // against the tile grown by the width). The sample positions themselves
    // depend only on the step index, so every tile agrees on them.
    double t0 = 0.0;
    double t1 = 1.0;
    const double ddx = rS.mBx - rS.mAx;
    const double ddy = rS.mBy - rS.mAy;
    const double box[4] = {
        static_cast<double>(rTile.mX0 - w - 1), static_cast<double>(rTile.mX1 + w + 2),
        static_cast<double>(rTile.mY0 - w - 1), static_cast<double>(rTile.mY1 + w + 2)};
    const double p[4] = {-ddx, ddx, -ddy, ddy};
    const double q[4] = {rS.mAx - box[0], box[1] - rS.mAx, rS.mAy - box[2], box[3] - rS.mAy};
    for (int k = 0; k < 4; ++k) {
        if (p[k] == 0.0) {
            if (q[k] < 0.0)
                return;
            continue;
        }
        const double r = q[k] / p[k];
        if (p[k] < 0.0)
            t0 = std::max(t0, r);
        else
            t1 = std::min(t1, r);
    }
    if (t0 > t1)
        return;
    const double n = static_cast<double>(rS.mSteps);
    const std::int64_t i0 =
        std::max<std::int64_t>(0, static_cast<std::int64_t>(std::floor(t0 * n)) - 1);
    const std::int64_t i1 =
        std::min<std::int64_t>(rS.mSteps, static_cast<std::int64_t>(std::ceil(t1 * n)) + 1);
    const std::int64_t lo = (w - 1) / 2;
    const std::int64_t hi = w / 2;
    for (std::int64_t i = i0; i <= i1; ++i) {
        const double t = static_cast<double>(i) / n;
        const double x = rS.mAx + t * ddx;
        const double y = rS.mAy + t * ddy;
        const double z = rS.mAz + t * (rS.mBz - rS.mAz) + rScene.mDepthBias;
        const std::int64_t cx = static_cast<std::int64_t>(std::floor(x));
        const std::int64_t cy = static_cast<std::int64_t>(std::floor(y));
        for (std::int64_t py = cy - lo; py <= cy + hi; ++py)
            for (std::int64_t px = cx - lo; px <= cx + hi; ++px)
                rst_plot(px, py, static_cast<float>(z), rLine.mColor, rLine.mId, rTile, rOut);
    }
}

void rst_draw_point(const RasterPoint& rPoint, const RasterScene& rScene, const RstTile& rTile,
                    RasterTarget& rOut) {
    const RasterVertex& v = rScene.mVertices[static_cast<std::size_t>(rPoint.mV)];
    if (!(v.mX == v.mX) || !(v.mY == v.mY))
        return;
    const double r = rScene.mPointRadius;
    const double x = std::min(std::max(v.mX, -kRasterLimitPx), kRasterLimitPx);
    const double y = std::min(std::max(v.mY, -kRasterLimitPx), kRasterLimitPx);
    const std::int64_t px0 = std::max(static_cast<std::int64_t>(std::floor(x - r)), rTile.mX0);
    const std::int64_t px1 = std::min(static_cast<std::int64_t>(std::floor(x + r)), rTile.mX1);
    const std::int64_t py0 = std::max(static_cast<std::int64_t>(std::floor(y - r)), rTile.mY0);
    const std::int64_t py1 = std::min(static_cast<std::int64_t>(std::floor(y + r)), rTile.mY1);
    const float z = static_cast<float>(v.mDepth + rScene.mDepthBias);
    for (std::int64_t py = py0; py <= py1; ++py) {
        for (std::int64_t px = px0; px <= px1; ++px) {
            const double ex = static_cast<double>(px) + 0.5 - x;
            const double ey = static_cast<double>(py) + 0.5 - y;
            if (ex * ex + ey * ey <= r * r)
                rst_plot(px, py, z, rPoint.mColor, rPoint.mId, rTile, rOut);
        }
    }
}

// Append every primitive whose pixel range overlaps a tile to that tile's list,
// in primitive order.
void rst_bin(std::int64_t Px0, std::int64_t Px1, std::int64_t Py0, std::int64_t Py1,
             std::size_t Index, std::int64_t TilesX, std::vector<std::vector<std::size_t>>& rBins) {
    if (Px0 > Px1 || Py0 > Py1)
        return;
    for (std::int64_t ty = Py0 / kRasterTile; ty <= Py1 / kRasterTile; ++ty)
        for (std::int64_t tx = Px0 / kRasterTile; tx <= Px1 / kRasterTile; ++tx)
            rBins[static_cast<std::size_t>(ty * TilesX + tx)].push_back(Index);
}

}  // namespace

RasterTarget rasterize(const RasterScene& rScene, int Width, int Height,
                       const std::array<std::uint8_t, 4>& rBackground) {
    RasterTarget out;
    out.mWidth = Width;
    out.mHeight = Height;
    const std::size_t num_px = static_cast<std::size_t>(Width) * static_cast<std::size_t>(Height);
    out.mRgba.resize(num_px * 4);
    for (std::size_t i = 0; i < num_px; ++i)
        std::copy(rBackground.begin(), rBackground.end(), out.mRgba.data() + i * 4);
    out.mDepth.assign(num_px, -std::numeric_limits<float>::infinity());
    out.mIds.assign(num_px, -1);
    if (num_px == 0)
        return out;

    std::vector<RstTri> tris(rScene.mTris.size());
    parallel_for(rScene.mTris.size(), [&](std::size_t i) {
        rst_prepare_tri(rScene, rScene.mTris[i], Width, Height, tris[i]);
    });
    std::vector<RstSeg> segs(rScene.mLines.size());
    parallel_for(rScene.mLines.size(), [&](std::size_t i) {
        rst_prepare_seg(rScene, rScene.mLines[i], Width, Height, segs[i]);
    });

    const std::int64_t tiles_x = (Width + kRasterTile - 1) / kRasterTile;
    const std::int64_t tiles_y = (Height + kRasterTile - 1) / kRasterTile;
    const std::size_t num_tiles = static_cast<std::size_t>(tiles_x * tiles_y);
    std::vector<std::vector<std::size_t>> tri_bins(num_tiles), seg_bins(num_tiles),
        pt_bins(num_tiles);
    for (std::size_t i = 0; i < tris.size(); ++i)
        if (tris[i].mArea > 0)
            rst_bin(tris[i].mPx0, tris[i].mPx1, tris[i].mPy0, tris[i].mPy1, i, tiles_x, tri_bins);
    for (std::size_t i = 0; i < segs.size(); ++i)
        rst_bin(segs[i].mPx0, segs[i].mPx1, segs[i].mPy0, segs[i].mPy1, i, tiles_x, seg_bins);
    const double r = rScene.mPointRadius;
    for (std::size_t i = 0; i < rScene.mPoints.size(); ++i) {
        const RasterVertex& v = rScene.mVertices[static_cast<std::size_t>(rScene.mPoints[i].mV)];
        if (!(v.mX == v.mX) || !(v.mY == v.mY))
            continue;
        const double x = std::min(std::max(v.mX, -kRasterLimitPx), kRasterLimitPx);
        const double y = std::min(std::max(v.mY, -kRasterLimitPx), kRasterLimitPx);
        rst_bin(std::max<std::int64_t>(static_cast<std::int64_t>(std::floor(x - r)), 0),
                std::min<std::int64_t>(static_cast<std::int64_t>(std::floor(x + r)), Width - 1),
                std::max<std::int64_t>(static_cast<std::int64_t>(std::floor(y - r)), 0),
                std::min<std::int64_t>(static_cast<std::int64_t>(std::floor(y + r)), Height - 1), i,
                tiles_x, pt_bins);
    }

    parallel_for(
        num_tiles,
        [&](std::size_t t) {
            const std::int64_t tx = static_cast<std::int64_t>(t) % tiles_x;
            const std::int64_t ty = static_cast<std::int64_t>(t) / tiles_x;
            RstTile tile;
            tile.mX0 = tx * kRasterTile;
            tile.mX1 = std::min<std::int64_t>(tile.mX0 + kRasterTile - 1, Width - 1);
            tile.mY0 = ty * kRasterTile;
            tile.mY1 = std::min<std::int64_t>(tile.mY0 + kRasterTile - 1, Height - 1);
            for (std::size_t i : tri_bins[t])
                rst_draw_tri(tris[i], tile, out);
            for (std::size_t i : seg_bins[t])
                rst_draw_seg(segs[i], rScene.mLines[i], rScene, tile, out);
            for (std::size_t i : pt_bins[t])
                rst_draw_point(rScene.mPoints[i], rScene, tile, out);
        },
        1);
    return out;
}

void raster_downsample(const RasterTarget& rTarget, int Factor, std::vector<std::uint8_t>& rRgba,
                       std::vector<std::int64_t>& rIds) {
    if (Factor <= 1) {
        rRgba = rTarget.mRgba;
        rIds = rTarget.mIds;
        return;
    }
    const std::size_t f = static_cast<std::size_t>(Factor);
    const std::size_t out_w = static_cast<std::size_t>(rTarget.mWidth) / f;
    const std::size_t out_h = static_cast<std::size_t>(rTarget.mHeight) / f;
    const std::size_t in_w = static_cast<std::size_t>(rTarget.mWidth);
    rRgba.assign(out_w * out_h * 4, 0);
    rIds.assign(out_w * out_h, -1);
    const std::uint64_t n = static_cast<std::uint64_t>(f * f);
    parallel_for(
        out_h,
        [&](std::size_t oy) {
            for (std::size_t ox = 0; ox < out_w; ++ox) {
                std::uint64_t sum_a = 0;
                std::uint64_t sum_c[3] = {0, 0, 0};
                for (std::size_t sy = 0; sy < f; ++sy) {
                    const std::uint8_t* row =
                        rTarget.mRgba.data() + ((oy * f + sy) * in_w + ox * f) * 4;
                    for (std::size_t sx = 0; sx < f; ++sx) {
                        const std::uint64_t a = row[sx * 4 + 3];
                        sum_a += a;
                        for (int c = 0; c < 3; ++c)
                            sum_c[c] += static_cast<std::uint64_t>(row[sx * 4 + c]) * a;
                    }
                }
                std::uint8_t* out = rRgba.data() + (oy * out_w + ox) * 4;
                if (sum_a > 0)
                    for (int c = 0; c < 3; ++c)
                        out[c] = static_cast<std::uint8_t>((sum_c[c] + sum_a / 2) / sum_a);
                out[3] = static_cast<std::uint8_t>((sum_a + n / 2) / n);
                rIds[oy * out_w + ox] = rTarget.mIds[(oy * f + f / 2) * in_w + ox * f + f / 2];
            }
        },
        1);
}

}  // namespace detail
}  // namespace meshioplusplus
