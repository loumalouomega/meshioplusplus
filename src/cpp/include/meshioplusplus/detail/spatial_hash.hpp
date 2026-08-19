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
#pragma once

/**
 * @file spatial_hash.hpp
 * @brief The integer bucket-grid spatial hash shared by `operations/merge.cpp`
 * (point welding) and `operations/interpolate.cpp` (nearest-point and
 * simplex-candidate lookup) — hoisted verbatim from merge's weld in v7.13.0,
 * the way `space_filling.hpp` was hoisted from reorder in v7.6.0.
 *
 * A `SpatialGrid` maps a quantized cell key (`floor(coord / cell)` per axis)
 * to the ids inserted into that cell, in insertion order. Both consumers
 * insert ids in ascending order from a **serial** pass, which is what makes
 * every downstream scan deterministic; the hash function only buckets — the
 * map is key-looked-up, never iterated — so the hash cannot affect any
 * consumer's output.
 *
 * Merge's contract (unchanged by the hoist, pinned by its test suites): cell
 * size is `atol`, so two points within `atol` of each other land in the same
 * or an adjacent cell and the fixed 3x3x3 neighbourhood (`ForEachIn27`, in
 * ascending dz -> dy -> dx order with early exit) finds every weld candidate.
 *
 * Interpolate additionally needs an *expanding* search: `ForEachInShell`
 * visits the cells at Chebyshev radius exactly `r` (r = 0 is the centre
 * cell), and the occupied-key bounding box (`OccupiedLo`/`OccupiedHi`,
 * maintained on insert) bounds how far a nearest-point search can ever need
 * to expand. Any point in a cell at Chebyshev key distance `r >= 1` from a
 * query's cell lies at least `(r - 1) * cell` away from the query, which is
 * the search's stopping rule.
 *
 * `operations/conservative_interpolate.cpp` is a third consumer, and adds
 * `ForEachInBox`: unlike a point query (one cell) or a shell (a fixed
 * Chebyshev ring), it needs every id inserted anywhere inside a *query box*
 * (a target simplex's own quantized bbox) — the read-side twin of
 * `InsertBox`, whose nested loop it mirrors exactly. Because `InsertBox`
 * inserts a source simplex's id into *every* cell its bbox spans, and the
 * query box can itself span several cells, `ForEachInBox` will yield the same
 * id more than once whenever the two multi-cell boxes share more than one
 * grid cell — callers must sort-and-deduplicate before treating the result as
 * a candidate set, exactly as `conservative_interpolate.cpp` does.
 *
 * `detail/` header: exempt from the operations layer's anon-namespace prefix
 * rule. `clean.cpp` keeps its own (serial, any-dtype, FNV-hashed) welder —
 * structurally different enough that sharing would only relocate code.
 */

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// Integer bucket-grid cell key: the quantized coordinate of a point.
struct GridKey {
    std::int64_t x, y, z;
    bool operator==(const GridKey& rOther) const {
        return x == rOther.x && y == rOther.y && z == rOther.z;
    }
};

/// A simple, well-mixed hash combine over the three signed axes (moved
/// verbatim from merge's weld; buckets only, never observable in output).
struct GridKeyHash {
    std::size_t operator()(const GridKey& rKey) const {
        std::uint64_t h = static_cast<std::uint64_t>(rKey.x) * 0x9e3779b97f4a7c15ULL;
        h ^= static_cast<std::uint64_t>(rKey.y) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= static_cast<std::uint64_t>(rKey.z) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return static_cast<std::size_t>(h);
    }
};

/// Quantize one coordinate to its bucket index (`floor(coord / cell)`).
inline std::int64_t grid_quantize(double coord, double cell) {
    return static_cast<std::int64_t>(std::floor(coord / cell));
}

/// The bucket grid: cell key -> ids in insertion order. Callers insert ids in
/// ascending order from a serial pass, so every bucket vector is ascending and
/// every scan over it is deterministic.
class SpatialGrid {
public:
    explicit SpatialGrid(double cellSize) : mCell(cellSize) {}

    double CellSize() const { return mCell; }

    /// The cell key of a (z-padded, 3-component) coordinate.
    GridKey KeyOf(const double* pCoords) const {
        return GridKey{grid_quantize(pCoords[0], mCell), grid_quantize(pCoords[1], mCell),
                       grid_quantize(pCoords[2], mCell)};
    }

    /// Append `id` to the cell at `rKey`.
    void Insert(const GridKey& rKey, std::int64_t id) {
        mCells[rKey].push_back(id);
        Cover(rKey);
    }

    /// Append `id` to every cell in the inclusive key box `[rLo, rHi]`
    /// (a simplex's quantized bounding box), in ascending z -> y -> x order.
    void InsertBox(const GridKey& rLo, const GridKey& rHi, std::int64_t id) {
        for (std::int64_t z = rLo.z; z <= rHi.z; ++z)
            for (std::int64_t y = rLo.y; y <= rHi.y; ++y)
                for (std::int64_t x = rLo.x; x <= rHi.x; ++x)
                    mCells[GridKey{x, y, z}].push_back(id);
        Cover(rLo);
        Cover(rHi);
    }

    /// The ids in the cell at `rKey`, or nullptr if the cell is empty.
    const std::vector<std::int64_t>* Find(const GridKey& rKey) const {
        auto it = mCells.find(rKey);
        return it == mCells.end() ? nullptr : &it->second;
    }

    /// Visits the 3x3x3 neighbourhood of `rCenter` in ascending dz -> dy -> dx
    /// order (merge's exact weld-scan order), calling `fn(bucket)` for each
    /// non-empty cell. `fn` returns false to stop early (candidate found).
    template <class F>
    void ForEachIn27(const GridKey& rCenter, F&& fn) const {
        for (std::int64_t dz = -1; dz <= 1; ++dz)
            for (std::int64_t dy = -1; dy <= 1; ++dy)
                for (std::int64_t dx = -1; dx <= 1; ++dx) {
                    const std::vector<std::int64_t>* p_ids =
                        Find(GridKey{rCenter.x + dx, rCenter.y + dy, rCenter.z + dz});
                    if (p_ids != nullptr && !fn(*p_ids))
                        return;
                }
    }

    /// Visits the cells at Chebyshev radius exactly `r` around `rCenter`
    /// (r = 0 is the centre cell alone), calling `fn(bucket)` for each
    /// non-empty cell, in ascending dz -> dy -> dx order. Enumerates only the
    /// shell (never the full cube), clamped to the occupied-key bounding box —
    /// cells outside it are empty by construction, so skipping them cannot
    /// change any consumer's result, and a far-away query stays O(bbox
    /// surface) instead of O(r^2) per shell.
    template <class F>
    void ForEachInShell(const GridKey& rCenter, std::int64_t r, F&& fn) const {
        if (!mCovered)
            return;
        const std::int64_t zl = std::max(-r, mLo.z - rCenter.z);
        const std::int64_t zh = std::min(r, mHi.z - rCenter.z);
        const std::int64_t yl = std::max(-r, mLo.y - rCenter.y);
        const std::int64_t yh = std::min(r, mHi.y - rCenter.y);
        const std::int64_t xl = std::max(-r, mLo.x - rCenter.x);
        const std::int64_t xh = std::min(r, mHi.x - rCenter.x);
        if (zl > zh || yl > yh || xl > xh)
            return;
        auto visit = [&](std::int64_t dx, std::int64_t dy, std::int64_t dz) {
            const std::vector<std::int64_t>* p_ids =
                Find(GridKey{rCenter.x + dx, rCenter.y + dy, rCenter.z + dz});
            if (p_ids != nullptr)
                fn(*p_ids);
        };
        if (r == 0) {
            visit(0, 0, 0);
            return;
        }
        for (std::int64_t dz = zl; dz <= zh; ++dz) {
            if (dz == -r || dz == r) {
                for (std::int64_t dy = yl; dy <= yh; ++dy)
                    for (std::int64_t dx = xl; dx <= xh; ++dx)
                        visit(dx, dy, dz);
            } else {
                for (std::int64_t dy = yl; dy <= yh; ++dy) {
                    if (dy == -r || dy == r) {
                        for (std::int64_t dx = xl; dx <= xh; ++dx)
                            visit(dx, dy, dz);
                    } else {
                        if (xl == -r)
                            visit(-r, dy, dz);
                        if (xh == r)
                            visit(r, dy, dz);
                    }
                }
            }
        }
    }

    /// Visits every non-empty cell in the inclusive key box `[rLo, rHi]`,
    /// calling `fn(bucket)` for each — the read-side twin of `InsertBox`, same
    /// ascending z -> y -> x traversal. A bucket may be visited once per
    /// id-in-that-bucket occurrence from the *inserting* side, so an id
    /// inserted via `InsertBox` over a multi-cell box can surface more than
    /// once here when the query box shares several cells with it; callers must
    /// deduplicate (see the file doc comment).
    template <class F>
    void ForEachInBox(const GridKey& rLo, const GridKey& rHi, F&& fn) const {
        for (std::int64_t z = rLo.z; z <= rHi.z; ++z)
            for (std::int64_t y = rLo.y; y <= rHi.y; ++y)
                for (std::int64_t x = rLo.x; x <= rHi.x; ++x) {
                    const std::vector<std::int64_t>* p_ids = Find(GridKey{x, y, z});
                    if (p_ids != nullptr)
                        fn(*p_ids);
                }
    }

    bool Empty() const { return mCells.empty(); }

    /// The occupied-key bounding box (valid only when not `Empty()`); bounds
    /// the maximum radius an expanding shell search can ever need.
    const GridKey& OccupiedLo() const { return mLo; }
    const GridKey& OccupiedHi() const { return mHi; }

private:
    void Cover(const GridKey& rKey) {
        if (!mCovered) {
            mLo = mHi = rKey;
            mCovered = true;
            return;
        }
        if (rKey.x < mLo.x)
            mLo.x = rKey.x;
        if (rKey.y < mLo.y)
            mLo.y = rKey.y;
        if (rKey.z < mLo.z)
            mLo.z = rKey.z;
        if (rKey.x > mHi.x)
            mHi.x = rKey.x;
        if (rKey.y > mHi.y)
            mHi.y = rKey.y;
        if (rKey.z > mHi.z)
            mHi.z = rKey.z;
    }

    double mCell;
    std::unordered_map<GridKey, std::vector<std::int64_t>, GridKeyHash> mCells;
    bool mCovered = false;
    GridKey mLo{0, 0, 0};
    GridKey mHi{0, 0, 0};
};

}  // namespace detail
}  // namespace meshioplusplus
