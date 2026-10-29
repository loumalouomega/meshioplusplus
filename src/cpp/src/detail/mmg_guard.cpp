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
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "mmg_guard.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

#if defined(_WIN32)
using MmgSavedHandler = void (*)(int);
#else
using MmgSavedHandler = struct sigaction;
#endif
static_assert(sizeof(MmgSavedHandler) <= 256, "mmg_guard: grow MmgCallGuard::mSaved");

MmgSavedHandler* mmg_saved(unsigned char* pStorage) {
    return reinterpret_cast<MmgSavedHandler*>(pStorage);
}

}  // namespace

std::mutex& mmg_mutex() {
    static std::mutex sMutex;
    return sMutex;
}

MmgCallGuard::MmgCallGuard() : mLock(mmg_mutex()) {
    std::memset(mSaved, 0, sizeof(mSaved));
    MmgSavedHandler* saved = mmg_saved(mSaved);
    for (std::size_t i = 0; i < kMmgNumSignals; ++i) {
#if defined(_WIN32)
        // signal() is the only query Windows has: set, read the old, put it back.
        saved[i] = std::signal(kMmgSignals[i], SIG_DFL);
        std::signal(kMmgSignals[i], saved[i]);
#else
        sigaction(kMmgSignals[i], nullptr, &saved[i]);
#endif
    }
}

MmgCallGuard::~MmgCallGuard() {
    MmgSavedHandler* saved = mmg_saved(mSaved);
    for (std::size_t i = 0; i < kMmgNumSignals; ++i) {
#if defined(_WIN32)
        if (saved[i] != SIG_ERR)
            std::signal(kMmgSignals[i], saved[i]);
#else
        sigaction(kMmgSignals[i], &saved[i], nullptr);
#endif
    }
}

void mmg_check_index_range(std::size_t Count, bool Int64, const char* pWhat, const char* pOp) {
    // MMG numbers from 1 and keeps room for the entities it inserts: refuse
    // well before the type's end rather than at it.
    const std::size_t limit =
        Int64 ? static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max() / 2)
              : static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max() / 2);
    if (Count <= limit)
        return;
    if (Int64)
        throw Unsupported(std::string(pOp) + ": " + std::to_string(Count) + " " + pWhat +
                          " are more than MMG can index");
    throw Unsupported(std::string(pOp) + ": " + std::to_string(Count) + " " + pWhat +
                      " overflow MMG's 32-bit MMG5_int; rebuild MMG with -DMMG5_INT=int64_t");
}

}  // namespace detail
}  // namespace meshioplusplus
