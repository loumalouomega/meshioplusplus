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
 * @file region.hpp
 * @brief `Region`: the unified, first-class model of a *named group of mesh
 * entities* — the thing gmsh calls a physical group, Exodus a block / node set
 * / side set, Abaqus an `*NSET` / `*ELSET` / `*SURFACE`, MED a family or group,
 * UNV a group, Ansys a component, OpenFOAM a boundary patch and Kratos a
 * SubModelPart.
 *
 * Before this type existed, named groups lived only on the Python `Mesh` as
 * `point_sets` / `cell_sets`, never crossed the pybind11 boundary, and were
 * therefore invisible to the C API, Fortran, WebAssembly and the native CLI.
 * Each set-capable format smuggled them through its own side-channel struct.
 * `Region` replaces all of that with one representation every mesh backend
 * carries and every binding surface can see.
 *
 * ## The three kinds
 *
 *  - `RegionKind::Point` — entries are **point indices**.
 *  - `RegionKind::Cell`  — entries are **global cell indices**, numbered
 *    *block-major*: block 0's cells first, then block 1's, and so on. This is
 *    the same numbering `partition_labels` uses for its flat buffer;
 *    `detail/cell_index.hpp` is the single owner of the conversion to and from
 *    `(block, row)` so nothing re-derives it.
 *  - `RegionKind::Side`  — entries are `(global cell index, local facet index)`
 *    **pairs**, stored as an `(n, 2)` array. Facets are numbered exactly as
 *    `detail/cell_faces.hpp` numbers the faces of a 3-D cell and
 *    `detail/cell_edges.hpp` the edges of a 2-D one. This kind has no
 *    equivalent in the legacy `point_sets`/`cell_sets` model at all.
 *
 * ## `mDim` and `mTag`
 *
 * `mDim` is the topological dimension the group was declared for, or `-1` when
 * the format does not say. Gmsh needs it: a physical group is per-dimension,
 * and two groups of different dimensions may share a name. `mTag` is the
 * format-native integer id (gmsh physical tag, MED family id, Exodus set id),
 * or `-1` when there is none. Both are carried through conversion so a
 * round-trip back to the originating format can restore them.
 *
 * ## Canonical entry ordering
 *
 * `Region::Canonicalize()` sorts entries ascending and removes duplicates
 * (lexicographically on the pair, for `Side`). Every mesh backend calls it from
 * `AddRegion`, so two regions built by different code paths — or by different
 * backends, or by different thread counts — compare exactly equal whenever they
 * describe the same group. That is what makes the cross-format round-trip
 * matrix an equality assertion rather than a set-comparison heuristic.
 *
 * Entries are held in an `NDArray` rather than a `std::vector<std::int64_t>` so
 * that the read path can hand the buffer to numpy through the existing capsule
 * adoption in `numpy_from_ndarray` without copying. The write path copies once,
 * because canonicalization must sort in place and the incoming buffer may be a
 * view over Python-owned memory.
 */

// System includes
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/**
 * @brief What kind of entity a `Region` groups.
 *
 * The enumerator values are part of the C ABI (`mio_region_kind`) and of the
 * WebAssembly / Python string vocabulary; `c_api.cpp` static_asserts them.
 */
enum class RegionKind {
    Point = 0,  ///< entries are point indices
    Cell = 1,   ///< entries are global (block-major) cell indices
    Side = 2,   ///< entries are (global cell index, local facet index) pairs
};

/**
 * @brief The lower-case spelling of a `RegionKind`, as used by the CLIs, the
 * WebAssembly bindings and the Python API.
 * @param Kind The kind to name.
 * @return `"point"`, `"cell"` or `"side"`.
 */
MESHIOPLUSPLUS_API const char* region_kind_name(RegionKind Kind);

/**
 * @brief Parse a `RegionKind` from its lower-case spelling.
 * @param rName One of `"point"`, `"cell"`, `"side"`.
 * @return The matching enumerator.
 * @throws std::invalid_argument naming the accepted values.
 */
MESHIOPLUSPLUS_API RegionKind region_kind_from_name(const std::string& rName);

/**
 * @brief A named group of points, cells or cell facets.
 *
 * See the file-level documentation for the entry conventions of each kind and
 * for why entries are canonicalized. Construct one directly and hand it to
 * `Mesh::AddRegion`, which canonicalizes and stores it.
 */
struct Region {
    /// The group's name, as the originating format spelled it.
    std::string mName;
    /// Which kind of entity the entries index.
    RegionKind mKind = RegionKind::Point;
    /// Topological dimension the group was declared for, or -1 if unspecified.
    int mDim = -1;
    /// Format-native integer id (gmsh physical tag, MED family id, ...), or -1.
    std::int64_t mTag = -1;
    /// Int64 entries: shape `(n,)` for Point/Cell, `(n, 2)` for Side.
    NDArray mEntries;

    Region() = default;
    Region(std::string name, RegionKind kind, NDArray entries)
        : mName(std::move(name)), mKind(kind), mEntries(std::move(entries)) {}
    Region(std::string name, RegionKind kind, int dim, std::int64_t tag, NDArray entries)
        : mName(std::move(name)), mKind(kind), mDim(dim), mTag(tag), mEntries(std::move(entries)) {}

    /// Number of grouped entities (rows of `mEntries`).
    std::size_t NumEntries() const {
        const std::size_t stride = Stride();
        return stride == 0 ? 0 : mEntries.Size() / stride;
    }

    /// Values per entry: 2 for `Side` (cell, facet), 1 otherwise.
    std::size_t Stride() const { return mKind == RegionKind::Side ? 2u : 1u; }

    /// Read-only view of the canonical Int64 entry buffer (may be null if empty).
    const std::int64_t* Entries() const {
        return mEntries.Size() ? mEntries.As<std::int64_t>() : nullptr;
    }

    /**
     * @brief Sort ascending, remove duplicates and normalize dtype/shape.
     *
     * Converts `mEntries` to owned Int64 storage of shape `(n,)` (Point/Cell)
     * or `(n, 2)` (Side), sorts it (lexicographically on the pair for `Side`)
     * and drops duplicate entries. Idempotent. Called by every backend's
     * `AddRegion`, so stored regions are always canonical.
     */
    MESHIOPLUSPLUS_API void Canonicalize();

    /**
     * @brief The identity of a region: two regions with the same key describe
     * the same group and replace one another in `AddRegion`.
     * @return `(kind, name, dim, tag)`, which is also the storage sort order.
     */
    std::tuple<int, const std::string&, int, std::int64_t> Key() const {
        return {static_cast<int>(mKind), mName, mDim, mTag};
    }
};

/**
 * @brief Whether two regions describe exactly the same group.
 *
 * Compares the key *and* the entries elementwise. Because stored regions are
 * canonical, this is a true set equality and not an ordering-sensitive one.
 * @param rA First region.
 * @param rB Second region.
 * @return `true` when name, kind, dim, tag and every entry match.
 */
MESHIOPLUSPLUS_API bool regions_equal(const Region& rA, const Region& rB);

/**
 * @brief Strict-weak ordering on `Region::Key()`, used to keep a mesh's region
 * list deterministically sorted across backends.
 * @param rA First region.
 * @param rB Second region.
 * @return `true` when `rA` sorts before `rB`.
 */
inline bool region_key_less(const Region& rA, const Region& rB) {
    return rA.Key() < rB.Key();
}

namespace detail {

/**
 * @brief The shared region-list storage every mesh backend delegates to.
 *
 * Kept here rather than duplicated in `meshio_mesh.hpp` / `native_mesh.hpp` /
 * `kratos_mesh.hpp` so the three backends cannot drift on the canonicalization
 * and ordering rules that the cross-backend determinism guarantee rests on.
 */
class RegionList {
public:
    /**
     * @brief Insert a region, or replace the one with the same
     * `(kind, name, dim, tag)` key.
     * @param region The region to store; canonicalized on the way in.
     */
    void Add(Region region) {
        region.Canonicalize();
        auto it = std::lower_bound(mRegions.begin(), mRegions.end(), region, region_key_less);
        if (it != mRegions.end() && it->Key() == region.Key())
            *it = std::move(region);
        else
            mRegions.insert(it, std::move(region));
    }

    /** @brief Number of stored regions. */
    std::size_t Size() const { return mRegions.size(); }
    /** @brief Region @p i in `(kind, name, dim, tag)` order. */
    const Region& At(std::size_t i) const { return mRegions[i]; }
    /** @brief The whole list, for backends that need to forward it wholesale. */
    const std::vector<Region>& All() const { return mRegions; }
    /** @brief Drop every stored region. */
    void Clear() { mRegions.clear(); }

    /**
     * @brief Sorted, de-duplicated region names.
     *
     * Matches the sorted-names guarantee `PointDataNames()` and friends make:
     * a name shared by several kinds or dimensions is reported once.
     */
    std::vector<std::string> Names() const {
        std::vector<std::string> out;
        out.reserve(mRegions.size());
        for (const auto& r_region : mRegions)
            out.push_back(r_region.mName);
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

    /** @brief Whether any region carries this name, regardless of kind. */
    bool Has(const std::string& rName) const {
        for (const auto& r_region : mRegions)
            if (r_region.mName == rName)
                return true;
        return false;
    }

    /** @brief Whether a region of this name *and* kind exists. */
    bool Has(const std::string& rName, RegionKind kind) const { return Find(rName, kind) != npos; }

    /**
     * @brief Index of the first region with this name and kind.
     * @param rName The region name to look for.
     * @param kind The kind to look for.
     * @return The index, or `npos` when absent.
     */
    std::size_t Find(const std::string& rName, RegionKind kind) const {
        for (std::size_t i = 0; i < mRegions.size(); ++i)
            if (mRegions[i].mKind == kind && mRegions[i].mName == rName)
                return i;
        return npos;
    }

    /// Sentinel returned by `Find` when no region matches.
    static constexpr std::size_t npos = static_cast<std::size_t>(-1);

private:
    std::vector<Region> mRegions;  // kept sorted by Region::Key()
};

}  // namespace detail
}  // namespace meshioplusplus
