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
 * @file detail/library_lock.hpp
 * @brief One process-wide lock serialising every call into a third-party
 * library that is not thread-safe (roadmap 3.4.1.1).
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * HDF5, netCDF-4, CGNS, MED and ADIOS2 share one HDF5 (and TecIO and gidpost
 * keep global state), so a lock per library would still race; there is one.
 * It is recursive because entry points nest (a series writer calls the HDF5
 * helpers, a sequence reads through a format reader).
 *
 * Take a `LibraryLock` as the first statement of every public function that
 * can reach such a library, before any handle whose destructor calls into it,
 * so the handles close under the lock.
 *
 * Lock order: a binding releases the GIL *before* it takes this lock; nothing
 * waits for the lock while holding the GIL (the only path back into Python
 * from under it, the gid series callback, re-acquires the GIL itself).
 */

#include <mutex>

namespace meshioplusplus {
namespace detail {

/// The process-wide mutex; prefer `LibraryLock`.
std::recursive_mutex& library_mutex();

/// RAII guard over `library_mutex()`.
class LibraryLock {
public:
    LibraryLock() : mLock(library_mutex()) {}
    LibraryLock(const LibraryLock&) = delete;
    LibraryLock& operator=(const LibraryLock&) = delete;

private:
    std::lock_guard<std::recursive_mutex> mLock;
};

}  // namespace detail
}  // namespace meshioplusplus
