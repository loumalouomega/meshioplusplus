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
// Characterisation tests for param_boundary_loops (detail/param_disc.hpp): the
// loops, their order and the `ok` flag, on inputs that are a disc, a closed
// surface and three ways to be neither. They pin what the directed-edge count
// did before it moved to EdgeUseTable (roadmap 6.2).

// System includes
#include <cstddef>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/param_disc.hpp"

using meshioplusplus::detail::PatchTri;
using meshioplusplus::detail::param_boundary_loops;

namespace {

using Loops = std::vector<std::vector<std::size_t>>;

Loops loops_of(const std::vector<PatchTri>& rTris, bool& rOk) {
    return param_boundary_loops(rTris, rOk);
}

}  // namespace

TEST(ParamBoundaryLoops, ASingleTriangleHasOneLoop) {
    bool ok = false;
    EXPECT_EQ(loops_of({{0, 1, 2}}, ok), (Loops{{0, 1, 2}}));
    EXPECT_TRUE(ok);
}

TEST(ParamBoundaryLoops, AFanAroundAnInteriorVertexIsADisc) {
    bool ok = false;
    EXPECT_EQ(loops_of({{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 1}}, ok), (Loops{{1, 2, 3, 4}}));
    EXPECT_TRUE(ok);
}

TEST(ParamBoundaryLoops, AClosedSurfaceHasNoLoop) {
    bool ok = false;
    EXPECT_TRUE(loops_of({{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}}, ok).empty());
    EXPECT_TRUE(ok);
}

TEST(ParamBoundaryLoops, AnEdgeUsedTwiceTheSameWayIsNotOk) {
    bool ok = true;
    // Two triangles both walking 0 -> 1: two boundary edges leave vertex 1 as
    // well, and the later one in directed-key order wins.
    EXPECT_EQ(loops_of({{0, 1, 2}, {0, 1, 3}}, ok), (Loops{{0, 1, 3}, {2}}));
    EXPECT_FALSE(ok);
}

TEST(ParamBoundaryLoops, APinchedBoundaryVertexIsNotOkAndKeepsTheLastEdge) {
    bool ok = true;
    // Two triangles sharing only vertex 0 (a bowtie).
    EXPECT_EQ(loops_of({{0, 1, 2}, {0, 3, 4}}, ok), (Loops{{0, 3, 4}, {1, 2}}));
    EXPECT_FALSE(ok);
}
