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
// DXF R2000 drawings (roadmap 7.6.4).
// remove_islands (roadmap 7.8.1).

#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "../../src/cpp/src/detail/cell_components.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/remove_islands.hpp"

using namespace meshioplusplus;

namespace {

// A sphere of radius 1 plus ten crumbs of radius 0.004 scattered around it.
Mesh ri_scan() {
    std::vector<Mesh> parts;
    parts.push_back(icosa_sphere(1.0, 3));
    for (int k = 0; k < 10; ++k)
        parts.push_back(icosa_sphere(0.004, 0, {1.5 + 0.1 * k, 0.2 * (k % 3), -0.3}));
    std::vector<const Mesh*> ptrs;
    for (const Mesh& r_m : parts)
        ptrs.push_back(&r_m);
    return merge(ptrs).mMesh;
}

}  // namespace

TEST(RemoveIslands, ASphereKeepsItselfAndLosesItsCrumbs) {
    const Mesh scan = ri_scan();
    const Mesh sphere = icosa_sphere(1.0, 3);
    RemoveIslandsOptions o;
    const RemoveIslandsResult r = remove_islands(scan, o);
    EXPECT_EQ(r.mNumComponents, 11);
    EXPECT_EQ(r.mNumRemoved, 10);
    EXPECT_EQ(r.mMesh.NumPoints(), sphere.NumPoints());
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), sphere.Cells(0).NumCells());
    // A dry run removes nothing and labels the cells.
    o.mDryRun = true;
    const RemoveIslandsResult d = remove_islands(scan, o);
    EXPECT_EQ(d.mMesh.NumPoints(), scan.NumPoints());
    EXPECT_TRUE(d.mMesh.HasCellData("islands:removed"));
    EXPECT_EQ(d.mNumRemoved, 10);
    // By cell count, absolute.
    RemoveIslandsOptions c;
    c.mMinDiagonal = 0.0;
    c.mMinCells = 100;
    EXPECT_EQ(remove_islands(scan, c).mNumRemoved, 10);
    c.mMinCells = 0;
    c.mMinArea = 0.5;  // half the whole area: only the crumbs fall short
    EXPECT_EQ(remove_islands(scan, c).mNumRemoved, 10);
}

TEST(RemoveIslands, EdgeAndPointConnectivityDiffer) {
    // Two triangles touching at one point: one point component, two edge ones.
    const Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}, {0, -1, 0}},
                                 "triangle", {{0, 1, 2}, {0, 3, 4}});
    EXPECT_EQ(detail::cell_components_by_edge(m).mCount, 2);
    EXPECT_EQ(detail::cell_components(m).mCount, 1);
    RemoveIslandsOptions o;
    o.mConnectivity = IslandConnectivity::Point;
    EXPECT_EQ(remove_islands(m, o).mNumComponents, 1);
    o.mMinDiagonal = -1.0;
    EXPECT_THROW(remove_islands(m, o), std::invalid_argument);
    EXPECT_THROW(island_connectivity_from_name("face"), std::invalid_argument);
}
