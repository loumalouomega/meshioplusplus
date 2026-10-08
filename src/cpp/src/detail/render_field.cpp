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

// The pure helpers of field rendering. See detail/render_field.hpp.

// System includes
#include <algorithm>
#include <cmath>
#include <limits>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"

// Project includes (private, not installed)
#include "render_field.hpp"

namespace meshioplusplus {
namespace detail {

double render_scale_forward(RenderScale Scale, double Threshold, double v) {
    if (!std::isfinite(v))
        return std::numeric_limits<double>::quiet_NaN();
    switch (Scale) {
        case RenderScale::Log:
            return v > 0.0 ? std::log10(v) : std::numeric_limits<double>::quiet_NaN();
        case RenderScale::Symlog: {
            const double t = Threshold > 0.0 ? Threshold : 1.0;
            const double a = std::fabs(v);
            // log10(1 + |v|/t) is odd, continuous and linear (slope 1/(t ln 10)) near 0.
            const double m = std::log10(1.0 + a / t);
            return v < 0.0 ? -m : m;
        }
        default:
            return v;
    }
}

double render_percentile(std::vector<double>& rValues, double P) {
    if (rValues.empty())
        return std::numeric_limits<double>::quiet_NaN();
    std::sort(rValues.begin(), rValues.end());
    const double pos =
        std::min(std::max(P, 0.0), 100.0) / 100.0 * static_cast<double>(rValues.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(lo + 1, rValues.size() - 1);
    const double frac = pos - static_cast<double>(lo);
    return rValues[lo] + frac * (rValues[hi] - rValues[lo]);
}

namespace {

// The 1, 2 or 5 times a power of ten step that gives about Target ticks.
double rf_nice_step(double Range, int Target) {
    const double raw = Range / static_cast<double>(std::max(Target, 1));
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double f = raw / mag;
    const double nice = f <= 1.0 ? 1.0 : (f <= 2.0 ? 2.0 : (f <= 5.0 ? 5.0 : 10.0));
    return nice * mag;
}

}  // namespace

std::vector<double> render_legend_ticks(RenderScale Scale, double Threshold, double VMin,
                                        double VMax, int Target) {
    std::vector<double> ticks;
    if (!(VMax > VMin) || !std::isfinite(VMin) || !std::isfinite(VMax))
        return ticks;
    if (Scale == RenderScale::Log) {
        if (!(VMin > 0.0))
            return ticks;
        const int lo = static_cast<int>(std::ceil(std::log10(VMin) - 1e-12));
        const int hi = static_cast<int>(std::floor(std::log10(VMax) + 1e-12));
        for (int e = lo; e <= hi; ++e)
            ticks.push_back(std::pow(10.0, e));
        return ticks;
    }
    if (Scale == RenderScale::Symlog) {
        const double t = Threshold > 0.0 ? Threshold : 1.0;
        if (VMin < 0.0 && VMax > 0.0)
            ticks.push_back(0.0);
        const double top = std::max(std::fabs(VMin), std::fabs(VMax));
        for (double e = t; e <= top * (1.0 + 1e-12); e *= 10.0) {
            if (e <= VMax)
                ticks.push_back(e);
            if (-e >= VMin)
                ticks.push_back(-e);
        }
        std::sort(ticks.begin(), ticks.end());
        return ticks;
    }
    const double step = rf_nice_step(VMax - VMin, Target);
    const double first = std::ceil(VMin / step - 1e-9) * step;
    for (double v = first; v <= VMax + step * 1e-9; v += step) {
        // Snap to the grid so 0.30000000000000004 is 0.3, and -0 is 0.
        char buf[48];
        snprintf_c(buf, sizeof(buf), "%.12g", std::round(v / step) * step);
        const char* end = nullptr;
        const double snapped = parse_double(buf, end);
        ticks.push_back(snapped == 0.0 ? 0.0 : snapped);
        if (ticks.size() > 64)
            break;
    }
    return ticks;
}

const std::array<std::array<std::uint8_t, 4>, 10>& render_category_palette() {
    static const std::array<std::array<std::uint8_t, 4>, 10> palette = {{
        {0, 114, 178, 255},    // blue
        {230, 159, 0, 255},    // orange
        {0, 158, 115, 255},    // bluish green
        {204, 121, 167, 255},  // reddish purple
        {86, 180, 233, 255},   // sky blue
        {213, 94, 0, 255},     // vermillion
        {240, 228, 66, 255},   // yellow
        {0, 0, 0, 255},        // black
        {153, 153, 153, 255},  // grey
        {204, 204, 204, 255},  // light grey
    }};
    return palette;
}

bool render_triangle_contour(const double* pXyz, const double* pValue, double Level, double* pOut) {
    // Corners at or above the level are "high"; an edge with one of each is cut.
    bool high[3];
    for (int k = 0; k < 3; ++k) {
        if (!std::isfinite(pValue[k]))
            return false;
        high[k] = pValue[k] >= Level;
    }
    if (high[0] == high[1] && high[1] == high[2])
        return false;
    int n = 0;
    for (int k = 0; k < 3; ++k) {
        const int a = k;
        const int b = (k + 1) % 3;
        if (high[a] == high[b])
            continue;
        const double t = (Level - pValue[a]) / (pValue[b] - pValue[a]);
        for (int c = 0; c < 3; ++c)
            pOut[3 * n + c] = pXyz[3 * a + c] + t * (pXyz[3 * b + c] - pXyz[3 * a + c]);
        ++n;
    }
    return n == 2;
}

}  // namespace detail
}  // namespace meshioplusplus
