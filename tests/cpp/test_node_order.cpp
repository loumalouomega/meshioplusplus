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
 * @file test_node_order.cpp
 * @brief The node-ordering registry: every table is a permutation of its cell
 *        type's width and the two directions are inverses, and the tables the
 *        formats used to keep privately are unchanged.
 *
 * The geometric checks (mid-edge nodes land on midpoints, positive volume, and
 * agreement with Code_Aster's own MED and gmsh readers) live in
 * tests/python/test_node_order.py, which also compares these tables with the
 * Python twin.
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <algorithm>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/node_order.hpp"

using namespace meshioplusplus;

TEST(NodeOrder, EveryTableIsAPermutationOfItsCellTypeWidth) {
    const auto keys = detail::node_order_keys();
    ASSERT_FALSE(keys.empty());
    for (const auto& [format, type] : keys) {
        SCOPED_TRACE(std::string(format) + "/" + std::string(type));
        const detail::NodeOrder* order = detail::node_order(format, type);
        ASSERT_NE(order, nullptr);
        const int n = cell_type_num_nodes(cell_type_from_name(std::string(type)));
        ASSERT_GT(n, 0);
        ASSERT_EQ(order->mToMeshio.size(), static_cast<std::size_t>(n));
        ASSERT_EQ(order->mFromMeshio.size(), static_cast<std::size_t>(n));
        std::vector<int> sorted = order->mToMeshio;
        std::sort(sorted.begin(), sorted.end());
        for (int i = 0; i < n; ++i)
            EXPECT_EQ(sorted[static_cast<std::size_t>(i)], i);
        for (int k = 0; k < n; ++k)
            EXPECT_EQ(order->mFromMeshio[static_cast<std::size_t>(
                          order->mToMeshio[static_cast<std::size_t>(k)])],
                      k);
    }
}

TEST(NodeOrder, UnknownKeysAreTheIdentity) {
    EXPECT_EQ(detail::node_order("med", "triangle"), nullptr);
    EXPECT_EQ(detail::node_order("code_aster", "tetra10"), nullptr);
    EXPECT_EQ(detail::node_order("no_such_format", "hexahedron20"), nullptr);
}

TEST(NodeOrder, KeepsTheTablesTheFormatsUsedToHold) {
    // med.cpp's former med_node_perm() (read direction).
    EXPECT_EQ(
        detail::node_order("med", "hexahedron20")->mToMeshio,
        (std::vector<int>{4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11, 16, 17, 18, 19}));
    // frd.cpp's former frd_perm_pe15.
    EXPECT_EQ(detail::node_order("frd", "wedge15")->mToMeshio,
              (std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11}));
    // unv.cpp's former unv_perm(): meshio_conn[perm[i]] = unv_conn[i].
    EXPECT_EQ(detail::node_order("unv", "tetra10")->mFromMeshio,
              (std::vector<int>{0, 4, 1, 5, 2, 6, 7, 8, 9, 3}));
    // mphtxt.cpp's former perm_of(): the linear tensor-order swaps.
    EXPECT_EQ(detail::node_order("mphtxt", "hexahedron")->mToMeshio,
              (std::vector<int>{0, 1, 3, 2, 4, 5, 7, 6}));
}

TEST(NodeOrder, ComsolMatchesPalaceComposedWithGmsh) {
    // AWS Palace (palace/utils/meshio.cpp, Apache-2.0) scatters COMSOL nodes
    // into gmsh order, gmsh[P[j]] = comsol[j]; meshio++ reads gmsh as
    // meshio[k] = gmsh[G[k]] (gmsh/common.py). So meshio[k] = comsol[P^-1[G[k]]].
    struct Case {
        const char* mType;
        std::vector<int> mPalace;
        std::vector<int> mGmsh;
    };
    const Case cases[] = {
        {"tetra10", {0, 1, 2, 3, 4, 6, 5, 7, 9, 8}, {0, 1, 2, 3, 4, 5, 6, 7, 9, 8}},
        {"quad9", {0, 1, 3, 2, 4, 7, 8, 5, 6}, {0, 1, 2, 3, 4, 5, 6, 7, 8}},
        {"wedge18",
         {0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 15, 10, 16, 17, 11, 12, 13, 14},
         {0, 1, 2, 3, 4, 5, 6, 9, 7, 12, 14, 13, 8, 10, 11, 15, 17, 16}},
        {"pyramid14",
         {0, 1, 3, 2, 4, 5, 6, 13, 8, 10, 7, 9, 12, 11},
         {0, 1, 2, 3, 4, 5, 8, 10, 6, 7, 9, 11, 12, 13}},
        {"hexahedron27",
         {0, 1, 3, 2, 4, 5, 7, 6, 8, 9, 20, 11, 13, 10, 21, 12, 22, 26, 23, 15, 24, 14, 16, 17, 25,
          18, 19},
         {0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15, 22, 23, 21, 24, 20,
          25, 26}},
    };
    for (const Case& c : cases) {
        SCOPED_TRACE(c.mType);
        std::vector<int> inverse(c.mPalace.size());
        for (std::size_t j = 0; j < c.mPalace.size(); ++j)
            inverse[static_cast<std::size_t>(c.mPalace[j])] = static_cast<int>(j);
        std::vector<int> expected;
        for (int g : c.mGmsh)
            expected.push_back(inverse[static_cast<std::size_t>(g)]);
        EXPECT_EQ(detail::node_order("mphtxt", c.mType)->mToMeshio, expected);
    }
}

TEST(NodeOrder, MedHexahedron27IsNotItsOwnInverse) {
    // Why the registry stores both directions: a writer that reused the read
    // table would scramble the face centres.
    const detail::NodeOrder* order = detail::node_order("med", "hexahedron27");
    ASSERT_NE(order, nullptr);
    EXPECT_NE(order->mToMeshio, order->mFromMeshio);
}
