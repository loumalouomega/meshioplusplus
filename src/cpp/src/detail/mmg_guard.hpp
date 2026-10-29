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
 * @file detail/mmg_guard.hpp
 * @brief The embedding hazards of MMG (roadmap 3.2.5), handled once for every
 * MMG operation.
 *
 * A **core-private** header: no installed header names it. It needs no MMG
 * header, so it is built and tested in every configuration.
 *
 * - **Signal handlers (3.2.5.1).** Every MMG 5.8 library entry point --
 *   `MMG3D_mmg3dlib`, `MMG3D_mmg3dls`, `MMG3D_mmg3dmov`, `MMGS_mmgslib`,
 *   `MMGS_mmgsls`, `MMG2D_mmg2dlib`, `MMG2D_mmg2dmesh`, `MMG2D_mmg2dls`,
 *   `MMG2D_mmg2dmov`, and `MMG3D_mmg3dcheck` -- calls `signal()` for SIGABRT,
 *   SIGFPE, SIGILL, SIGSEGV, SIGTERM and SIGINT, installing a handler that
 *   prints a message and calls `exit(EXIT_FAILURE)`. It never puts the old
 *   handlers back, so after one call Python's Ctrl-C handler (and any crash
 *   reporter of the host) is gone for good. `MmgCallGuard` saves the six
 *   handlers before the call and restores them when it leaves scope, on every
 *   path, an exception included. During the call MMG's handlers stand: a
 *   Ctrl-C or a crash inside MMG ends the process, which `doc/mmg.md` states.
 * - **Threads (3.2.5.3).** Whether two meshes can be adapted at once in one
 *   process is not established, so every call is serialised behind one mutex,
 *   taken before the handlers are saved (a second thread would otherwise save
 *   MMG's handlers as the host's).
 * - **Standard output.** MMG prints nothing at `verbose = -1`, but ELAS,
 *   which `move` calls, compares the verbosity MMG hands it as a number
 *   against the character `'0'` and so always prints. `MmgStdoutSilencer`
 *   points file descriptor 1 at the null device for the call (POSIX; a no-op
 *   elsewhere), process-wide: anything another thread prints meanwhile is
 *   lost too.
 * - **Index width (3.2.5.2).** `MMG5_int` is `int32_t` unless MMG was built
 *   with `-DMMG5_INT=int64_t`; `mmg_check_index_range` refuses a mesh beyond
 *   that range by name rather than let it truncate.
 */

#include <csignal>
#include <cstddef>
#include <mutex>

namespace meshioplusplus {
namespace detail {

/// The signals MMG's entry points take over.
inline constexpr int kMmgSignals[] = {SIGABRT, SIGFPE, SIGILL, SIGSEGV, SIGTERM, SIGINT};
inline constexpr std::size_t kMmgNumSignals = sizeof(kMmgSignals) / sizeof(kMmgSignals[0]);

/// The process-wide mutex every MMG call holds; prefer `MmgCallGuard`.
std::mutex& mmg_mutex();

/// Holds `mmg_mutex()` and the host's six signal handlers for one MMG call.
class MmgCallGuard {
public:
    MmgCallGuard();
    ~MmgCallGuard();
    MmgCallGuard(const MmgCallGuard&) = delete;
    MmgCallGuard& operator=(const MmgCallGuard&) = delete;

private:
    std::lock_guard<std::mutex> mLock;
    // Opaque storage for the saved handlers (a `struct sigaction` each on
    // POSIX, a handler pointer on Windows), sized for the larger.
    alignas(std::max_align_t) unsigned char mSaved[kMmgNumSignals * 256];
};

/// Points file descriptor 1 at the null device while it lives, when `Active`.
class MmgStdoutSilencer {
public:
    explicit MmgStdoutSilencer(bool Active);
    ~MmgStdoutSilencer();
    MmgStdoutSilencer(const MmgStdoutSilencer&) = delete;
    MmgStdoutSilencer& operator=(const MmgStdoutSilencer&) = delete;

private:
    int mSaved = -1;
};

/**
 * @brief Refuse `Count` entities of `pWhat` when they overflow `MMG5_int`.
 *
 * @param Count The number of points, cells or the largest index to pass.
 * @param Int64 Whether MMG was built with 64-bit `MMG5_int`.
 * @param pWhat What is counted ("points", "tetrahedra", ...).
 * @param pOp The operation, for the message.
 * @throws Unsupported naming `-DMMG5_INT=int64_t` when `Count` is out of range.
 */
void mmg_check_index_range(std::size_t Count, bool Int64, const char* pWhat, const char* pOp);

}  // namespace detail
}  // namespace meshioplusplus
