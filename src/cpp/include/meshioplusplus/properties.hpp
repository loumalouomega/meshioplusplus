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
 * @file properties.hpp
 * @brief `PropertyValue` / `PropertySet`: the key/value content of a Kratos
 * `Properties` block (material data), independent of both the mesh backend and
 * the MDPA format.
 *
 * This exists because two unrelated places need the *same* representation and
 * would otherwise each invent one: `formats/mdpa.hpp` (which reads and writes
 * `Begin Properties <id>` bodies through `MdpaInfo`) and
 * `backends/model_part.hpp` (whose `ModelPart` stores them so material data
 * survives the Kratos bridge instead of being reduced to a bare id). Keeping
 * one owner is what lets a deck round-trip
 * file -> `MdpaInfo` -> `ModelPart` -> `MdpaInfo` -> file with nothing lost.
 *
 * It deliberately does **not** model Kratos's typed `Variable<T>` system. A
 * `Variable<double>` cannot be looked up without Kratos's own component
 * registry (`KratosComponents<Variable<double>>::Get`), which is not linked
 * here and never will be -- see the `to_model_part` overload taking an
 * "apply property" callback in `kratos_bridge.hpp`, which is how a real Kratos
 * consumer turns these key/value pairs into typed variables.
 *
 * @note These do not go through `field_data`: `NDArray` has ten numeric dtypes
 *       and no string or bytes dtype, so a `CONSTITUTIVE_LAW LinearElastic3DLaw`
 *       line has no representation there at all. Since v9.2.0 the `Mesh` carries
 *       them **as themselves** instead, through the uniform API's
 *       `AddPropertySet`/`GetPropertySet` (see `mesh_api.hpp`) -- a
 *       `PropertySet` is an ordinary struct with `std::string` members, so
 *       storing it needs no dtype. Before that they could only ride the
 *       `MdpaInfo` side channel, which no registry-based consumer could reach.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {

/**
 * @brief One `KEY value` entry of a properties block.
 *
 * Exactly one of `mValues` and `mText` carries the value:
 *
 *  - a plain number becomes a one-element Float64 `mValues`;
 *  - an inline `Begin Table` becomes an `(n, k)` Float64 `mValues` with
 *    `mIsTable` set, `mKey` holding the table header's arguments verbatim
 *    (e.g. `"1 TEMPERATURE YOUNG_MODULUS"`) so it re-emits unchanged;
 *  - anything else -- a constitutive-law name, a bracketed vector or matrix
 *    such as `[3] (1.0, 0.0, 0.0)` -- is kept verbatim in `mText`.
 *
 * The text fallback is not a gap to close later: it is what makes an
 * unrecognized value **lossless**, since it is re-emitted byte for byte. It is
 * also what the pure-Python reference does (`_mdpa.py` falls back to the raw
 * string when `float()` fails), so the two readers agree.
 */
struct PropertyValue {
    /** @brief The variable name, or a table's header arguments when `mIsTable`. */
    std::string mKey;
    /** @brief Numeric value: `{1}` for a scalar, `(n, k)` for a table. */
    NDArray mValues;
    /** @brief The value verbatim, when it is not numeric. Empty otherwise. */
    std::string mText;
    /** @brief Whether this entry is an inline `Begin Table` rather than a value. */
    bool mIsTable = false;

    /** @brief Whether the value lives in `mText` rather than `mValues`. */
    bool IsText() const { return mValues.Size() == 0; }
};

/** @brief One `Begin Properties <id>` block: an id plus its entries, in file order. */
struct PropertySet {
    std::int64_t mId = 0;
    std::vector<PropertyValue> mValues;
};

namespace detail {

/**
 * @brief Every properties block of a mesh, kept ascending by `mId`.
 *
 * The structural twin of `detail::RegionList`, and for the same reason: the
 * three mesh backends each hold one of these, so a shared container is what
 * stops them drifting on ordering or on the replace-by-key rule.
 *
 * Unlike regions, property sets are keyed by **id, not by entity index**, so
 * nothing here ever needs remapping when an operation renumbers cells or
 * points -- there is no `detail/region_remap.hpp` counterpart to look for.
 */
class PropertySetList {
public:
    /**
     * @brief Insert a set, or replace the one with the same `mId`.
     * @param propertySet The set to store.
     */
    void Add(PropertySet propertySet) {
        auto it = std::lower_bound(
            mSets.begin(), mSets.end(), propertySet.mId,
            [](const PropertySet& rLhs, std::int64_t id) { return rLhs.mId < id; });
        if (it != mSets.end() && it->mId == propertySet.mId)
            *it = std::move(propertySet);
        else
            mSets.insert(it, std::move(propertySet));
    }

    /** @brief Number of stored property sets. */
    std::size_t Size() const { return mSets.size(); }
    /** @brief Set @p Index, in ascending-`mId` order. */
    const PropertySet& At(std::size_t Index) const { return mSets[Index]; }
    /** @brief The whole list, for backends that forward it wholesale. */
    const std::vector<PropertySet>& All() const { return mSets; }
    /** @brief Drop every stored set. */
    void Clear() { mSets.clear(); }

    /** @brief Whether a set with this id exists. */
    bool Has(std::int64_t Id) const { return Find(Id) != npos; }

    /**
     * @brief Index of the set with this id.
     * @param Id The properties id to look for.
     * @return The index, or `npos` when absent.
     */
    std::size_t Find(std::int64_t Id) const {
        auto it = std::lower_bound(
            mSets.begin(), mSets.end(), Id,
            [](const PropertySet& rLhs, std::int64_t id) { return rLhs.mId < id; });
        if (it == mSets.end() || it->mId != Id)
            return npos;
        return static_cast<std::size_t>(it - mSets.begin());
    }

    /// Sentinel returned by `Find` when no set matches.
    static constexpr std::size_t npos = static_cast<std::size_t>(-1);

private:
    std::vector<PropertySet> mSets;  // kept sorted by mId
};

}  // namespace detail
}  // namespace meshioplusplus
