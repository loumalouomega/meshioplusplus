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
// flip_normals (roadmap 7.8.2).

#include <cstdint>
#include <stdexcept>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/flip_normals.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/region.hpp"

using namespace meshioplusplus;

namespace {

bool fn_same(const NDArray& rA, const NDArray& rB) {
    if (rA.Size() != rB.Size() || rA.Dtype() != rB.Dtype())
        return false;
    for (std::size_t i = 0; i < rA.Size(); ++i)
        if (rA.As<std::int64_t>()[i] != rB.As<std::int64_t>()[i])
            return false;
    return true;
}

}  // namespace

TEST(FlipNormals, TwiceIsTheInputAndOnceReversesEveryTriangle) {
    const Mesh s = icosa_sphere(1.0, 1);
    const FlipNormalsResult once = flip_normals(s);
    EXPECT_EQ(once.mNumFlipped, static_cast<std::int64_t>(s.Cells(0).NumCells()));
    const std::int64_t* a = s.Cells(0).Conn().As<std::int64_t>();
    const std::int64_t* b = once.mMesh.Cells(0).Conn().As<std::int64_t>();
    EXPECT_EQ(b[0], a[0]);
    EXPECT_EQ(b[1], a[2]);
    EXPECT_EQ(b[2], a[1]);
    const FlipNormalsResult twice = flip_normals(once.mMesh);
    EXPECT_TRUE(fn_same(twice.mMesh.Cells(0).Conn(), s.Cells(0).Conn()));
}

TEST(FlipNormals, ARegionAndItsSides) {
    Mesh m = mt::make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}, "tetra",
                           {{0, 1, 2, 3}, {1, 2, 3, 4}});
    NDArray cells = NDArray::Uninit(DType::Int64, {std::size_t{1}});
    cells.As<std::int64_t>()[0] = 1;
    m.AddRegion(Region("second", RegionKind::Cell, std::move(cells)));
    NDArray sides = NDArray::Uninit(DType::Int64, {std::size_t{2}, std::size_t{2}});
    std::int64_t* p = sides.As<std::int64_t>();
    p[0] = 0;  // cell 0, facet 1 (not flipped)
    p[1] = 1;
    p[2] = 1;  // cell 1, facet 1 (flipped)
    p[3] = 1;
    m.AddRegion(Region("wall", RegionKind::Side, std::move(sides)));
    const FlipNormalsResult r = flip_normals(m, "second");
    EXPECT_EQ(r.mNumFlipped, 1);
    const std::int64_t* c = r.mMesh.Cells(0).Conn().As<std::int64_t>();
    EXPECT_EQ(c[0], 0);  // cell 0 untouched
    EXPECT_EQ(c[1], 1);
    const std::size_t wall = r.mMesh.FindRegion("wall", RegionKind::Side);
    ASSERT_NE(wall, Mesh::npos);
    // Flipping again restores every entry.
    const FlipNormalsResult back = flip_normals(r.mMesh, "second");
    EXPECT_TRUE(fn_same(back.mMesh.Cells(0).Conn(), m.Cells(0).Conn()));
    EXPECT_TRUE(fn_same(back.mMesh.Region(back.mMesh.FindRegion("wall", RegionKind::Side)).mEntries,
                        m.Region(m.FindRegion("wall", RegionKind::Side)).mEntries));
    EXPECT_THROW(flip_normals(m, "nope"), std::invalid_argument);
}
