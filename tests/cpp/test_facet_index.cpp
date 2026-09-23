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
// Tests for detail/facet_index: the sorted-corner facet lookup the formats that
// name a facet by its nodes share (LS-DYNA, FEBio, Elmer).

// External includes
#include <gtest/gtest.h>

// System includes
#include <cstdint>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/facet_index.hpp"

using namespace meshioplusplus;
using namespace mt;

TEST(FacetIndex, TwoTetsShareOneFace) {
    // Two tets glued on the face (1, 2, 3).
    Mesh mesh = make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}, "tetra",
                          {{0, 1, 2, 3}, {1, 2, 3, 4}});
    const detail::FacetIndex index(mesh);
    EXPECT_EQ(index.Size(), 7u);  // 8 faces, one shared
    const std::int64_t shared[] = {3, 1, 2};
    const detail::FacetHit* hit = index.Find(shared, 3);
    ASSERT_NE(hit, nullptr);
    EXPECT_EQ(hit->mCount, 2u);
    EXPECT_EQ(hit->mFirst.mCell, 0);
    EXPECT_EQ(hit->mFirst.mFacet, 1);  // tetra face 1 = (1, 2, 3)
    EXPECT_EQ(hit->mSecond.mCell, 1);
    const std::int64_t outer[] = {0, 1, 3};
    hit = index.Find(outer, 3);
    ASSERT_NE(hit, nullptr);
    EXPECT_EQ(hit->mCount, 1u);
    EXPECT_EQ(hit->mFirst.mFacet, 0);
    const std::int64_t missing[] = {0, 1, 4};
    EXPECT_EQ(index.Find(missing, 3), nullptr);
}

TEST(FacetIndex, SurfaceCellsIndexTheirEdgesAndOptionallyThemselves) {
    Mesh mesh = make_mesh({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, "quad", {{0, 1, 2, 3}});
    const detail::FacetIndex edges(mesh);
    const std::int64_t edge[] = {2, 1};
    const detail::FacetHit* hit = edges.Find(edge, 2);
    ASSERT_NE(hit, nullptr);
    EXPECT_EQ(hit->mFirst.mFacet, 1);  // quad edge 1 = (1, 2)
    const std::int64_t face[] = {3, 2, 1, 0};
    EXPECT_EQ(edges.Find(face, 4), nullptr);

    detail::FacetIndexOptions options;
    options.mSurfaceEdges = false;
    options.mSurfaceSelf = true;
    const detail::FacetIndex self(mesh, options);
    hit = self.Find(face, 4);
    ASSERT_NE(hit, nullptr);
    EXPECT_EQ(hit->mFirst.mFacet, 0);
    EXPECT_EQ(self.Find(edge, 2), nullptr);
}

TEST(FacetIndex, FacetNodesIncludesMidSideNodes) {
    const Mesh mesh = tet10_mesh();
    CellType type = CellType::Custom;
    std::vector<std::int64_t> nodes;
    ASSERT_TRUE(detail::facet_nodes(mesh, 0, 0, type, nodes));
    EXPECT_EQ(type, CellType::Triangle6);
    EXPECT_EQ(nodes, (std::vector<std::int64_t>{0, 1, 3, 4, 8, 7}));
    EXPECT_FALSE(detail::facet_nodes(mesh, 0, 4, type, nodes));
    EXPECT_FALSE(detail::facet_nodes(mesh, 1, 0, type, nodes));

    const Mesh quads = quad8_mesh();
    ASSERT_TRUE(detail::facet_nodes(quads, 0, 2, type, nodes));
    EXPECT_EQ(type, CellType::Line3);
    EXPECT_EQ(nodes.size(), 3u);
}
