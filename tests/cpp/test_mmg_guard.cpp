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
// The MMG backend's embedding hazards (roadmap 3.2.5): the call guard that puts
// the host's signal handlers back, the index-width check and the capability
// queries. None needs MMG, so they run in every build.

// System includes
#include <csignal>
#include <cstddef>
#include <stdexcept>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/mmg_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/mmg.hpp"

namespace {

using HandlerFn = void (*)(int);

void mmg_test_host_handler(int) {}
void mmg_test_intruder_handler(int) {}

HandlerFn mmg_test_current(int Signal) {
    // signal() is the portable query: set, read the old, put it back.
    const HandlerFn previous = std::signal(Signal, SIG_DFL);
    std::signal(Signal, previous);
    return previous;
}

/// What an MMG entry point does on every call: take over the six signals.
void mmg_test_like_mmg() {
    for (const int s : meshioplusplus::detail::kMmgSignals)
        std::signal(s, mmg_test_intruder_handler);
}

}  // namespace

TEST(MmgGuard, PutsTheHostsSignalHandlersBack) {
    HandlerFn before[meshioplusplus::detail::kMmgNumSignals];
    for (std::size_t i = 0; i < meshioplusplus::detail::kMmgNumSignals; ++i)
        before[i] = std::signal(meshioplusplus::detail::kMmgSignals[i], mmg_test_host_handler);
    {
        meshioplusplus::detail::MmgCallGuard guard;
        mmg_test_like_mmg();
        EXPECT_EQ(mmg_test_current(SIGINT), mmg_test_intruder_handler);
    }
    for (const int s : meshioplusplus::detail::kMmgSignals)
        EXPECT_EQ(mmg_test_current(s), mmg_test_host_handler) << "signal " << s;

    // ... and on the way out of an exception, a strong failure's path.
    try {
        meshioplusplus::detail::MmgCallGuard guard;
        mmg_test_like_mmg();
        throw std::runtime_error("MMG5_STRONGFAILURE");
    } catch (const std::runtime_error&) {
    }
    for (const int s : meshioplusplus::detail::kMmgSignals)
        EXPECT_EQ(mmg_test_current(s), mmg_test_host_handler) << "signal " << s;

    for (std::size_t i = 0; i < meshioplusplus::detail::kMmgNumSignals; ++i)
        std::signal(meshioplusplus::detail::kMmgSignals[i], before[i]);
}

TEST(MmgGuard, RefusesAMeshBeyondTheIndexWidth) {
    using meshioplusplus::detail::mmg_check_index_range;
    EXPECT_NO_THROW(mmg_check_index_range(1000, false, "points", "adapt"));
    const std::size_t big = std::size_t{1} << 31;
    EXPECT_THROW(mmg_check_index_range(big, false, "points", "adapt"), meshioplusplus::Unsupported);
    EXPECT_NO_THROW(mmg_check_index_range(big, true, "points", "adapt"));
    try {
        mmg_check_index_range(big, false, "tetrahedra", "adapt");
    } catch (const meshioplusplus::Unsupported& e) {
        EXPECT_NE(std::string(e.what()).find("MMG5_INT=int64_t"), std::string::npos);
    }
}

TEST(MmgGuard, CapabilitiesAgree) {
    const bool has = meshioplusplus::mmg_available();
    EXPECT_EQ(meshioplusplus::mmg_index_bits() != 0, has);
    if (has)
        EXPECT_TRUE(meshioplusplus::mmg_index_bits() == 32 ||
                    meshioplusplus::mmg_index_bits() == 64);
    else
        EXPECT_FALSE(meshioplusplus::mmg_has_elas());
}
