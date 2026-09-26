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
// The cell-lattice neighbour search behind proximity_graph; see
// operations/neighbors.hpp for the contract and the parity rules.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/neighbors.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kNbPrefix = "meshio++: neighbor_pairs: ";

/// The points bucketed onto a lattice, sorted by cell.
struct NbLattice {
    std::size_t mDim = 0;
    bool mPeriodic = false;
    std::array<std::int64_t, 3> mCells{1, 1, 1};  // per-axis cell count
    double mMinSide = 0.0;
    std::vector<std::array<std::int64_t, 3>> mCellOf;  // per point
    std::vector<std::array<std::int64_t, 3>> mKeys;    // occupied cells, ascending
    std::vector<std::int64_t> mStart;                  // mKeys.size() + 1 slots into mOrder
    std::vector<std::int64_t> mOrder;                  // point ids grouped by cell

    /// The slice of mOrder holding cell @p rKey, or an empty one.
    std::pair<const std::int64_t*, const std::int64_t*> Find(
        const std::array<std::int64_t, 3>& rKey) const {
        const auto it = std::lower_bound(mKeys.begin(), mKeys.end(), rKey);
        if (it == mKeys.end() || *it != rKey)
            return {nullptr, nullptr};
        const std::size_t c = static_cast<std::size_t>(it - mKeys.begin());
        return {mOrder.data() + mStart[c], mOrder.data() + mStart[c + 1]};
    }
};

NbLattice nb_build_lattice(const std::vector<double>& rXyz, std::size_t n, std::size_t dim,
                           double cell, const std::vector<double>& rBox) {
    NbLattice lat;
    lat.mDim = dim;
    lat.mPeriodic = !rBox.empty();
    std::array<double, 3> lo{0.0, 0.0, 0.0}, side{1.0, 1.0, 1.0};
    if (lat.mPeriodic) {
        // A whole number of cells per periodic axis, each at least `cell`
        // wide, so no short cell sits at the seam.
        for (std::size_t a = 0; a < dim; ++a) {
            lat.mCells[a] =
                std::max<std::int64_t>(1, static_cast<std::int64_t>(std::floor(rBox[a] / cell)));
            side[a] = rBox[a] / static_cast<double>(lat.mCells[a]);
        }
    } else {
        for (std::size_t a = 0; a < dim; ++a) {
            double m = rXyz[a];
            for (std::size_t i = 1; i < n; ++i)
                m = std::min(m, rXyz[i * 3 + a]);
            lo[a] = m;
            side[a] = cell;
        }
    }
    lat.mMinSide = side[0];
    for (std::size_t a = 1; a < dim; ++a)
        lat.mMinSide = std::min(lat.mMinSide, side[a]);
    lat.mCellOf.resize(n);
    parallel_for(n, [&](std::size_t i) {
        std::array<std::int64_t, 3> k{0, 0, 0};
        for (std::size_t a = 0; a < dim; ++a) {
            std::int64_t c =
                static_cast<std::int64_t>(std::floor((rXyz[i * 3 + a] - lo[a]) / side[a]));
            if (lat.mPeriodic)
                c = ((c % lat.mCells[a]) + lat.mCells[a]) % lat.mCells[a];
            k[a] = c;
        }
        lat.mCellOf[i] = k;
    });
    if (!lat.mPeriodic)
        for (std::size_t a = 0; a < dim; ++a) {
            std::int64_t m = 0;
            for (std::size_t i = 0; i < n; ++i)
                m = std::max(m, lat.mCellOf[i][a]);
            lat.mCells[a] = m + 1;
        }
    lat.mOrder.resize(n);
    for (std::size_t i = 0; i < n; ++i)
        lat.mOrder[i] = static_cast<std::int64_t>(i);
    parallel_sort(lat.mOrder.begin(), lat.mOrder.end(), [&](std::int64_t x, std::int64_t y) {
        const auto& kx = lat.mCellOf[static_cast<std::size_t>(x)];
        const auto& ky = lat.mCellOf[static_cast<std::size_t>(y)];
        return kx != ky ? kx < ky : x < y;
    });
    lat.mStart.push_back(0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& k = lat.mCellOf[static_cast<std::size_t>(lat.mOrder[i])];
        if (lat.mKeys.empty() || lat.mKeys.back() != k) {
            if (!lat.mKeys.empty())
                lat.mStart.push_back(static_cast<std::int64_t>(i));
            lat.mKeys.push_back(k);
        }
    }
    lat.mStart.push_back(static_cast<std::int64_t>(n));
    return lat;
}

/// numpy's squared distance: dx*dx + dy*dy + dz*dz in axis order, each
/// displacement reduced to its minimum image (np.round: half to even).
inline double nb_d2(const std::vector<double>& rXyz, std::size_t i, std::size_t j, std::size_t dim,
                    const std::vector<double>& rBox) {
    double s = 0.0;
    for (std::size_t a = 0; a < dim; ++a) {
        double disp = rXyz[i * 3 + a] - rXyz[j * 3 + a];
        if (!rBox.empty())
            disp = disp - rBox[a] * std::nearbyint(disp / rBox[a]);
        const double sq = disp * disp;
        s = a == 0 ? sq : s + sq;
    }
    return s;
}

/// Visits the lattice cells at Chebyshev distance exactly @p R from @p rCenter
/// (wrapped on a periodic lattice), each at most once per call sequence the
/// caller tracks through @p rSeen.
template <class F>
void nb_ring(const NbLattice& rLat, const std::array<std::int64_t, 3>& rCenter, std::int64_t R,
             std::vector<std::array<std::int64_t, 3>>& rSeen, F&& fn) {
    const std::size_t dim = rLat.mDim;
    const std::int64_t rz = dim > 2 ? R : 0, ry = dim > 1 ? R : 0;
    for (std::int64_t dz = -rz; dz <= rz; ++dz)
        for (std::int64_t dy = -ry; dy <= ry; ++dy)
            for (std::int64_t dx = -R; dx <= R; ++dx) {
                if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) != R)
                    continue;
                std::array<std::int64_t, 3> k{rCenter[0] + dx, rCenter[1] + dy, rCenter[2] + dz};
                bool inside = true;
                for (std::size_t a = 0; a < dim; ++a) {
                    if (rLat.mPeriodic)
                        k[a] = ((k[a] % rLat.mCells[a]) + rLat.mCells[a]) % rLat.mCells[a];
                    else if (k[a] < 0 || k[a] >= rLat.mCells[a])
                        inside = false;
                }
                if (!inside)
                    continue;
                if (rLat.mPeriodic) {
                    // A small lattice wraps a ring onto cells already visited.
                    if (std::find(rSeen.begin(), rSeen.end(), k) != rSeen.end())
                        continue;
                    rSeen.push_back(k);
                }
                const auto cell = rLat.Find(k);
                if (cell.first)
                    fn(cell.first, cell.second);
            }
}

}  // namespace

NeighborPairs neighbor_pairs(const NDArray& rPoints, const NeighborOptions& rOptions) {
    if (rPoints.Ndim() != 2 || rPoints.Shape()[1] < 1 || rPoints.Shape()[1] > 3)
        throw std::invalid_argument(std::string(kNbPrefix) +
                                    "points must be an (N, d) array with d in 1..3");
    const std::size_t n = rPoints.Shape()[0];
    const std::size_t dim = rPoints.Shape()[1];
    std::vector<double> xyz(n * 3, 0.0);
    detail::dispatch_dtype(rPoints.Dtype(), [&]<class T>() {
        const T* src = rPoints.As<T>();
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t a = 0; a < dim; ++a)
                xyz[i * 3 + a] = static_cast<double>(src[i * dim + a]);
        });
    });
    for (double v : xyz)
        if (!std::isfinite(v))
            throw std::invalid_argument(std::string(kNbPrefix) +
                                        "the positions contain a non-finite coordinate");
    const std::vector<double>& box = rOptions.mBox;
    if (!box.empty()) {
        if (box.size() != dim)
            throw std::invalid_argument(std::string(kNbPrefix) +
                                        "the box needs one side per coordinate");
        for (double b : box)
            if (!(b > 0.0) || !std::isfinite(b))
                throw std::invalid_argument(std::string(kNbPrefix) +
                                            "box sides must be positive and finite");
    }

    NeighborPairs out;
    const auto finish = [&](std::vector<std::int64_t>& rA, std::vector<std::int64_t>& rB) {
        out.mSource = NDArray::Uninit(DType::Int64, {rA.size()});
        out.mTarget = NDArray::Uninit(DType::Int64, {rB.size()});
        std::copy(rA.begin(), rA.end(), out.mSource.As<std::int64_t>());
        std::copy(rB.begin(), rB.end(), out.mTarget.As<std::int64_t>());
    };
    std::vector<std::int64_t> src, dst;

    if (rOptions.mMethod == NeighborMethod::Radius) {
        const double r = rOptions.mRadius;
        if (!(r > 0.0) || !std::isfinite(r))
            throw std::invalid_argument(std::string(kNbPrefix) + "the radius must be positive");
        if (!box.empty() && r > 0.5 * *std::min_element(box.begin(), box.end()))
            throw std::invalid_argument(std::string(kNbPrefix) +
                                        "the radius exceeds half the smallest box side");
        if (n == 0) {
            finish(src, dst);
            return out;
        }
        // Cells a hair wider than the radius, so rounding in the cell
        // assignment cannot put a pair within the radius two cells apart.
        const double cell = std::max(rOptions.mCellSize, r) * (1.0 + 1e-6);
        const NbLattice lat = nb_build_lattice(xyz, n, dim, cell, box);
        const double r2 = r * r;
        std::vector<std::int64_t> count(n + 1, 0);
        const auto visit = [&](std::size_t i, auto&& emit) {
            std::vector<std::array<std::int64_t, 3>> seen;
            for (std::int64_t R = 0; R <= 1; ++R)
                nb_ring(lat, lat.mCellOf[i], R, seen,
                        [&](const std::int64_t* first, const std::int64_t* last) {
                            for (const std::int64_t* p = first; p != last; ++p) {
                                const std::size_t j = static_cast<std::size_t>(*p);
                                if (j > i && nb_d2(xyz, i, j, dim, box) <= r2)
                                    emit(j);
                            }
                        });
        };
        parallel_for(n, [&](std::size_t i) {
            std::int64_t c = 0;
            visit(i, [&](std::size_t) { ++c; });
            count[i + 1] = c;
        });
        for (std::size_t i = 0; i < n; ++i)
            count[i + 1] += count[i];
        src.resize(static_cast<std::size_t>(count[n]));
        dst.resize(src.size());
        parallel_for(n, [&](std::size_t i) {
            std::size_t at = static_cast<std::size_t>(count[i]);
            visit(i, [&](std::size_t j) {
                src[at] = static_cast<std::int64_t>(i);
                dst[at++] = static_cast<std::int64_t>(j);
            });
        });
        finish(src, dst);
        return out;
    }

    // k nearest: grow a Chebyshev ring of cells around each point until its
    // k-th candidate is strictly closer than anything outside the rings can be.
    const std::int64_t k =
        std::min<std::int64_t>(rOptions.mK, n > 0 ? static_cast<std::int64_t>(n) - 1 : 0);
    if (k <= 0) {
        finish(src, dst);
        return out;
    }
    double cell = rOptions.mCellSize;
    if (!(cell > 0.0)) {
        // About k/2 points per cell (proximity_graph's own sizing rule).
        std::array<double, 3> lo{0, 0, 0}, hi{0, 0, 0};
        for (std::size_t a = 0; a < dim; ++a) {
            lo[a] = hi[a] = xyz[a];
            for (std::size_t i = 1; i < n; ++i) {
                lo[a] = std::min(lo[a], xyz[i * 3 + a]);
                hi[a] = std::max(hi[a], xyz[i * 3 + a]);
            }
        }
        double measure = 1.0, extent = 0.0;
        int live = 0;
        for (std::size_t a = 0; a < dim; ++a)
            if (hi[a] > lo[a]) {
                measure *= hi[a] - lo[a];
                extent = std::max(extent, hi[a] - lo[a]);
                ++live;
            }
        cell =
            live == 0
                ? 1.0
                : std::max(std::pow(measure / static_cast<double>(n) * static_cast<double>(k) * 0.5,
                                    1.0 / live),
                           extent * 1e-9);
    }
    const NbLattice lat = nb_build_lattice(xyz, n, dim, cell, box);
    std::int64_t max_ring = 0;
    for (std::size_t a = 0; a < dim; ++a)
        max_ring = std::max(max_ring, lat.mCells[a]);
    const std::size_t kk = static_cast<std::size_t>(k);
    src.resize(n * kk);
    dst.resize(n * kk);
    parallel_for(n, [&](std::size_t q) {
        static thread_local std::vector<std::pair<double, std::int64_t>> cand;
        cand.clear();
        std::vector<std::array<std::int64_t, 3>> seen;
        const auto less = [](const std::pair<double, std::int64_t>& x,
                             const std::pair<double, std::int64_t>& y) {
            return x.first != y.first ? x.first < y.first : x.second < y.second;
        };
        for (std::int64_t R = 0;; ++R) {
            nb_ring(lat, lat.mCellOf[q], R, seen,
                    [&](const std::int64_t* first, const std::int64_t* last) {
                        for (const std::int64_t* p = first; p != last; ++p)
                            if (static_cast<std::size_t>(*p) != q)
                                cand.emplace_back(
                                    nb_d2(xyz, q, static_cast<std::size_t>(*p), dim, box), *p);
                    });
            if (R >= max_ring)
                break;  // every cell has been visited: exact by exhaustion
            if (cand.size() >= kk) {
                std::nth_element(cand.begin(), cand.begin() + (k - 1), cand.end(), less);
                // Anything outside rings 0..R is at least R cells' width away;
                // the margin covers rounding in the cell assignment.
                const double covered = static_cast<double>(R) * lat.mMinSide * (1.0 - 1e-6);
                if (cand[kk - 1].first < covered * covered)
                    break;
            }
        }
        std::partial_sort(cand.begin(), cand.begin() + k, cand.end(), less);
        for (std::size_t t = 0; t < kk; ++t) {
            src[q * kk + t] = static_cast<std::int64_t>(q);
            dst[q * kk + t] = cand[t].second;
        }
    });
    finish(src, dst);
    return out;
}

}  // namespace meshioplusplus
