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
 * @file formats/time_series_common.hpp
 * @brief What the transient writers of the format headers share: the series
 * mode and the raw solver array of `XdmfTimeSeriesWriter` and
 * `VtkhdfTimeSeriesWriter`, and the special members every move-only series
 * writer declares.
 *
 * Every series writer (`PvdSeriesWriter`, `FrdSeriesWriter`,
 * `XdmfTimeSeriesWriter` and the rest) is a pimpl class: one
 * `std::unique_ptr<Impl>`, an out-of-line destructor that finalizes, no copy,
 * and a `noexcept` move that leaves the source empty.
 * `MESHIOPLUSPLUS_SERIES_WRITER_SPECIAL_MEMBERS(Class)` declares those five
 * members once; each class still declares its own constructors, its `Write`
 * (or `WritePointsCells` / `WriteData`), its observers and its `Impl`, with
 * their own documentation.
 */

// System includes
#include <cstddef>
#include <string>
#include <vector>

namespace meshioplusplus {

/**
 * @brief Whether a series starts fresh or continues an existing one.
 *
 * `Append` is what a restartable analysis needs: a run that was killed, hit a
 * node failure, or simply finished a stage leaves its output, and the next run
 * must continue that series rather than overwrite it.
 *
 * Appending to a path that does **not** exist yet is not an error -- it is
 * exactly `Truncate`. That is deliberate: a solver can then pass `Append`
 * unconditionally instead of having to probe the filesystem and branch.
 */
enum class SeriesMode { Truncate, Append };

/**
 * @brief One solver array for a series writer's raw `WriteData` overload.
 *
 * Values are row-major, `mNumComponents` per entity: `{n}` for a scalar
 * field, `{n, 3}` for a vector, and so on.
 */
struct SeriesNamedArray {
    std::string mName;
    std::size_t mNumComponents = 1;
    std::vector<double> mValues;
};

}  // namespace meshioplusplus

/**
 * @brief Declares the special members of a move-only pimpl series writer:
 * the out-of-line destructor (which finalizes the series if it has not been
 * finalized already), the deleted copy constructor and copy assignment, and
 * the `noexcept` move constructor and move assignment, which leave the source
 * empty. Use it in the class's `public:` section; the `.cpp` defines them.
 */
#define MESHIOPLUSPLUS_SERIES_WRITER_SPECIAL_MEMBERS(Class) \
    ~Class();                                               \
    Class(const Class&) = delete;                           \
    Class& operator=(const Class&) = delete;                \
    Class(Class&&) noexcept;                                \
    Class& operator=(Class&&) noexcept
