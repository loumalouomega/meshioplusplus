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
 * helpers, a sequence reads through a format reader, `read_exodus(path, opts)`
 * locks and calls the overload that locks again) and because the gid series
 * callback runs caller code, which may read an HDF5 file, under the lock.
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

/// The process-wide mutex; prefer `LibraryLock`. Recursive: see above.
std::recursive_mutex& library_mutex();  // NOSONAR(cpp:S8462): entry points nest

/// How many `LibraryLock`s the calling thread holds (0 when none).
int& library_lock_depth();

/// RAII guard over `library_mutex()`.
class LibraryLock {
public:
    LibraryLock() {
        library_mutex().lock();
        ++library_lock_depth();
    }
    ~LibraryLock() {
        --library_lock_depth();
        library_mutex().unlock();
    }
    LibraryLock(const LibraryLock&) = delete;
    LibraryLock& operator=(const LibraryLock&) = delete;
};

/**
 * @brief Gives the lock up for a stretch of pure computation (roadmap 5.8.2).
 *
 * Releases every level the calling thread holds and takes them all back when
 * it goes out of scope, so another thread's HDF5 read or write runs while
 * this one compresses or inflates. A no-op on a thread that holds no lock.
 * Inside it, make **no** library call and touch no library handle: the
 * invariant -- no two threads inside the libraries -- rests on that. Handles
 * opened before it stay open and valid; only their use waits.
 */
class LibraryUnlock {
public:
    LibraryUnlock() : mDepth(library_lock_depth()) {
        for (int i = 0; i < mDepth; ++i)
            library_mutex().unlock();
        library_lock_depth() = 0;
    }
    ~LibraryUnlock() {
        for (int i = 0; i < mDepth; ++i)
            library_mutex().lock();
        library_lock_depth() = mDepth;
    }
    LibraryUnlock(const LibraryUnlock&) = delete;
    LibraryUnlock& operator=(const LibraryUnlock&) = delete;

private:
    int mDepth;
};

}  // namespace detail
}  // namespace meshioplusplus
