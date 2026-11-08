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
// Unit tests for detail/facet_slots.hpp -- FacetSlotUtils, the slot numbering
// extract_surface and smooth's boundary pass share. Their outputs are pinned
// by test_facet_tables.cpp; these tests pin the helpers themselves.

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/facet_slots.hpp"

using meshioplusplus::detail::FacetSlotUtils;

namespace {

// A stand-in polyhedron block: only NumFaces is read.
struct FakeBlock {
    std::vector<std::size_t> mFaces;
    std::size_t NumFaces(std::size_t C) const { return mFaces[C]; }
};

}  // namespace

TEST(FacetSlotUtils, FaceStartsArePrefixSumsAndInvert) {
    const FakeBlock block{{4, 6, 5, 0, 7}};
    std::vector<std::size_t> start;
    const std::size_t total = FacetSlotUtils::BuildFaceStarts(block, 5, start);
    EXPECT_EQ(total, 22u);
    EXPECT_EQ(start, (std::vector<std::size_t>{0, 4, 10, 15, 15, 22}));

    std::vector<std::uint32_t> owner;
    FacetSlotUtils::FillCellOfFace(start, 5, owner);
    ASSERT_EQ(owner.size(), total);
    for (std::size_t c = 0; c < 5; ++c)
        for (std::size_t f = start[c]; f < start[c + 1]; ++f)
            EXPECT_EQ(owner[f], c) << "slot " << f;
}

TEST(FacetSlotUtils, AnEmptyBlockHasOneStartAndNoSlots) {
    const FakeBlock block{{}};
    std::vector<std::size_t> start;
    EXPECT_EQ(FacetSlotUtils::BuildFaceStarts(block, 0, start), 0u);
    EXPECT_EQ(start, (std::vector<std::size_t>{0}));
    std::vector<std::uint32_t> owner;
    FacetSlotUtils::FillCellOfFace(start, 0, owner);
    EXPECT_TRUE(owner.empty());
}

TEST(FacetSlotUtils, SlotOwnerFindsTheBlockOfAGlobalSlot) {
    // Blocks of 10, 0 (empty), 5 and 3 slots.
    const std::vector<std::size_t> first{0, 10, 10, 15};
    EXPECT_EQ(FacetSlotUtils::SlotOwner(first, 0), 0u);
    EXPECT_EQ(FacetSlotUtils::SlotOwner(first, 9), 0u);
    // An empty block owns nothing: slot 10 belongs to the next non-empty one.
    EXPECT_EQ(FacetSlotUtils::SlotOwner(first, 10), 2u);
    EXPECT_EQ(FacetSlotUtils::SlotOwner(first, 14), 2u);
    EXPECT_EQ(FacetSlotUtils::SlotOwner(first, 15), 3u);
    EXPECT_EQ(FacetSlotUtils::SlotOwner(first, 17), 3u);
    // One block: no search at all.
    EXPECT_EQ(FacetSlotUtils::SlotOwner({0}, 123), 0u);
}
