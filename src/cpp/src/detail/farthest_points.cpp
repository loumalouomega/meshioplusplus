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
#include <vector>

// Project includes
#include "meshioplusplus/parallel.hpp"

#include "farthest_points.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

/// `np.argmax`'s order: a NaN beats any number, and an equal value loses to
/// the lower index.
struct FpsBest {
    double mValue = 0.0;
    std::size_t mIndex = static_cast<std::size_t>(-1);

    bool Beats(const FpsBest& rOther) const {
        if (rOther.mIndex == static_cast<std::size_t>(-1))
            return true;
        if (mIndex == static_cast<std::size_t>(-1))
            return false;
        const bool nan = std::isnan(mValue);
        const bool other_nan = std::isnan(rOther.mValue);
        if (nan != other_nan)
            return nan;
        if (!nan && mValue != rOther.mValue)
            return mValue > rOther.mValue;
        return mIndex < rOther.mIndex;
    }
};

}  // namespace

FarthestPoints farthest_point_order(const double* pXyz, std::size_t N, std::size_t Count,
                                    std::size_t First) {
    FarthestPoints out;
    out.mOrder.reserve(Count);
    out.mOrder.push_back(static_cast<std::int64_t>(First));
    std::vector<double> d2(N);
    // Updates the distances to the newest selection and returns the farthest,
    // in one sweep; the reduction's combine is a total order, so every
    // backend and chunking returns the same index.
    const auto sweep = [&](std::size_t From, bool Init) {
        const double x0 = pXyz[3 * From];
        const double y0 = pXyz[3 * From + 1];
        const double z0 = pXyz[3 * From + 2];
        return parallel_reduce(
            N, 4096, FpsBest{},
            [&](std::size_t Begin, std::size_t End) {
                // A strict `>` keeps the first maximum; a NaN never passes it
                // and is tracked apart, the first one winning.
                double best = -1.0;
                std::size_t at = static_cast<std::size_t>(-1);
                std::size_t first_nan = static_cast<std::size_t>(-1);
                for (std::size_t k = Begin; k < End; ++k) {
                    const double dx = pXyz[3 * k] - x0;
                    const double dy = pXyz[3 * k + 1] - y0;
                    const double dz = pXyz[3 * k + 2] - z0;
                    double d = dx * dx;
                    d = d + dy * dy;
                    d = d + dz * dz;
                    double v = d;
                    if (!Init) {
                        // np.minimum: a NaN on either side is the result.
                        const double cur = d2[k];
                        v = (d < cur || d != d) ? d : cur;
                    }
                    d2[k] = v;
                    if (v > best) {
                        best = v;
                        at = k;
                    } else if (v != v && first_nan == static_cast<std::size_t>(-1)) {
                        first_nan = k;
                    }
                }
                if (first_nan != static_cast<std::size_t>(-1))
                    return FpsBest{std::nan(""), first_nan};
                return FpsBest{best, at};
            },
            [](const FpsBest& rA, const FpsBest& rB) { return rB.Beats(rA) ? rB : rA; });
    };
    FpsBest next = sweep(First, true);
    for (std::size_t i = 1; i < Count; ++i) {
        // `d2[nxt] <= 0.0` in the twin: false for a NaN, which is selected.
        if (next.mValue <= 0.0) {
            out.mSelected = i;
            return out;
        }
        out.mOrder.push_back(static_cast<std::int64_t>(next.mIndex));
        if (i + 1 < Count)
            next = sweep(next.mIndex, false);
    }
    out.mSelected = Count;
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
