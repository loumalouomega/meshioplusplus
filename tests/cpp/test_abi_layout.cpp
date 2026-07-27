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
/**
 * @file test_abi_layout.cpp
 * @brief The layout half of the ABI gate: a committed snapshot of the size and
 * alignment of every type that crosses the installed-library boundary.
 *
 * `doc/abi.md` splits ABI-relevant changes into three tiers. This file makes
 * **Tier A** (layout) fail mechanically and objectively: add, remove, reorder or
 * retype a data member of any type below and this test stops compiling, naming
 * the type. The reviewer then either reverts it or bumps
 * `MESHIOPLUSPLUS_ABI_VERSION` in `abi_version.hpp` and updates the number here.
 *
 * It is deliberately only half the gate. **Tier B** -- editing the body of an
 * existing inline function or template -- changes no size at all and is
 * invisible here; `tools/check-abi-version.sh` is what catches that, by diffing
 * the installed headers against the previous release tag. Neither tool alone is
 * sufficient, which is why both exist.
 *
 * ### Why the numbers are pinned on one configuration only
 *
 * `sizeof(std::string)` is 32 with libstdc++ and 24 with libc++; `std::vector`
 * and `std::unordered_map` differ likewise, and MSVC differs from both. Hard
 * numbers are therefore only meaningful against a named reference
 * configuration -- Linux x86_64 + libstdc++ + LP64, which is what CI's
 * `cpp-tests` matrix runs on every backend. Everywhere else the test reports the
 * layout and skips rather than failing for a reason that has nothing to do with
 * meshio++. That is not a loophole: the gate runs on every PR in CI, and a
 * developer on macOS seeing a skip has lost nothing, whereas a hard failure
 * there would train people to ignore this file.
 *
 * A snapshot mismatch is NOT automatically a bug -- it is a question. "Did you
 * mean to change the ABI?" The answer is often yes.
 */

// System includes
#include <cstddef>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/abi_version.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/properties.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/write_options.hpp"

namespace {

// The reference configuration the pinned numbers below describe. Anything else
// (libc++, MSVC, 32-bit, non-x86_64) gets the reporting path instead.
#if defined(__linux__) && defined(__x86_64__) && defined(__GLIBCXX__) && !defined(_MSC_VER)
#define MIO_ABI_LAYOUT_PINNED 1
#else
#define MIO_ABI_LAYOUT_PINNED 0
#endif

#if MIO_ABI_LAYOUT_PINNED

// One macro so a failure names the type and both numbers rather than pointing at
// an anonymous static_assert. sizeof AND alignof: a member reordering that
// happens to preserve size can still change alignment, and a change to an
// enum's underlying type changes alignment while size may survive.
#define MIO_ABI_LAYOUT(Type, Size, Align)                                             \
    static_assert(sizeof(Type) == (Size),                                             \
                  "meshio++ ABI: sizeof(" #Type                                       \
                  ") changed. This is a Tier A layout "                               \
                  "break (doc/abi.md). Bump MESHIOPLUSPLUS_ABI_VERSION in "           \
                  "abi_version.hpp and update this snapshot, or revert the member."); \
    static_assert(alignof(Type) == (Align),                                           \
                  "meshio++ ABI: alignof(" #Type                                      \
                  ") changed. This is a Tier A layout "                               \
                  "break (doc/abi.md). Bump MESHIOPLUSPLUS_ABI_VERSION in "           \
                  "abi_version.hpp and update this snapshot, or revert the member.")

// --- Backend-independent boundary types -------------------------------------
// Every one of these is passed by value or reference through an exported
// function, so a consumer and the library must agree on all of them regardless
// of which backend either was built with.

MIO_ABI_LAYOUT(meshioplusplus::NDArray, 72, 8);
MIO_ABI_LAYOUT(meshioplusplus::Region, 120, 8);
MIO_ABI_LAYOUT(meshioplusplus::ReadOptions, 56, 8);
MIO_ABI_LAYOUT(meshioplusplus::WriteOptions, 48, 8);
MIO_ABI_LAYOUT(meshioplusplus::PropertyValue, 144, 8);
MIO_ABI_LAYOUT(meshioplusplus::PropertySet, 32, 8);
MIO_ABI_LAYOUT(meshioplusplus::MeshMetadata, 256, 8);

// CellType is stored inside cell blocks on the NATIVE and KRATOS backends, so
// its width is structural, not cosmetic. Appending an enumerator is fine (and
// deliberately not caught here); widening the underlying type is not.
static_assert(sizeof(meshioplusplus::CellType) == 2,
              "meshio++ ABI: CellType's underlying type changed width. Appending "
              "enumerators is safe and expected; changing `: std::uint16_t` is a "
              "Tier A break (doc/abi.md).");

// --- The mesh itself, which IS the backend ----------------------------------
// Each backend gets its own line because Mesh is a different type per backend;
// only the one this TU was compiled for is checked, and the cpp-tests matrix
// runs this file once per backend so all three are covered across CI.
#if defined(MESHIOPLUSPLUS_MESH_BACKEND_NATIVE)
MIO_ABI_LAYOUT(meshioplusplus::Mesh, 560, 8);
#elif defined(MESHIOPLUSPLUS_MESH_BACKEND_KRATOS)
MIO_ABI_LAYOUT(meshioplusplus::Mesh, 680, 8);
#else
MIO_ABI_LAYOUT(meshioplusplus::Mesh, 392, 8);
#endif

#endif  // MIO_ABI_LAYOUT_PINNED

/// Reports one type's layout, so a mismatch on an unpinned platform is
/// actionable rather than merely skipped.
template <class T>
void report(const char* name) {
    std::printf("  %-42s sizeof=%3zu alignof=%2zu\n", name, sizeof(T), alignof(T));
}

}  // namespace

TEST(AbiLayout, SnapshotIsPinnedOnTheReferenceConfiguration) {
    std::printf("meshio++ ABI version %d, layout on this platform:\n", MESHIOPLUSPLUS_ABI_VERSION);
    report<meshioplusplus::NDArray>("NDArray");
    report<meshioplusplus::Region>("Region");
    report<meshioplusplus::ReadOptions>("ReadOptions");
    report<meshioplusplus::WriteOptions>("WriteOptions");
    report<meshioplusplus::PropertyValue>("PropertyValue");
    report<meshioplusplus::PropertySet>("PropertySet");
    report<meshioplusplus::MeshMetadata>("MeshMetadata");
    report<meshioplusplus::Mesh>("Mesh (this backend)");

#if MIO_ABI_LAYOUT_PINNED
    // The static_asserts above already ran at compile time; reaching here means
    // the snapshot held. Nothing further to assert -- this exists so the test
    // reports rather than silently passing.
    SUCCEED() << "layout snapshot pinned for ABI version " << MESHIOPLUSPLUS_ABI_VERSION;
#else
    GTEST_SKIP() << "layout numbers are pinned only on Linux/x86_64/libstdc++ "
                    "(sizeof(std::string) alone differs across standard libraries); "
                    "CI's cpp-tests matrix is that configuration.";
#endif
}

TEST(AbiLayout, GuardSymbolIsReferencedByThisTranslationUnit) {
    // The v9.1.0 lesson in executable form: the backend guard shipped inert for a
    // full release because an unused `inline` variable is emitted lazily, so
    // nothing ever forced the symbol to resolve. This TU links, and it includes
    // mesh.hpp, so both sentinels must have resolved for this binary to exist at
    // all -- taking their addresses here makes that dependency explicit rather
    // than incidental.
    EXPECT_NE(&meshioplusplus::detail::MESHIOPLUSPLUS_ABI_SYM(MESHIOPLUSPLUS_ABI_VERSION), nullptr);
    EXPECT_NE(&meshioplusplus::detail::MESHIOPLUSPLUS_BACKEND_SYM(MESHIOPLUSPLUS_ACTIVE_BACKEND),
              nullptr);
}
