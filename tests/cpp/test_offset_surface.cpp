// offset_surface, thicken and hollow (roadmap 7.8.4).

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/boundary_layer.hpp"
#include "meshioplusplus/operations/clip.hpp"
#include "meshioplusplus/operations/offset_surface.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/thickness.hpp"

using namespace meshioplusplus;

namespace {

// The unit cube's surface, six outward quads.
Mesh ofs_cube() {
    Mesh m;
    NDArray p = NDArray::Uninit(DType::Float64, {8, 3});
    const double xyz[] = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    std::copy(xyz, xyz + 24, p.As<double>());
    m.AssignPoints(std::move(p));
    NDArray c = NDArray::Uninit(DType::Int64, {6, 4});
    const std::int64_t q[] = {0, 3, 2, 1, 4, 5, 6, 7, 0, 4, 7, 3,
                              1, 2, 6, 5, 0, 1, 5, 4, 3, 7, 6, 2};
    std::copy(q, q + 24, c.As<std::int64_t>());
    m.AddCellBlock("quad", std::move(c));
    return m;
}

// The signed volume of a surface of triangles and quads, and whether every
// directed edge has its reverse once.
double ofs_volume(const Mesh& rMesh, bool& rClosed) {
    const double* p = rMesh.Points().As<double>();
    std::map<std::array<std::int64_t, 2>, int> uses;
    double v = 0.0;
    const auto tet = [&](std::int64_t A, std::int64_t B, std::int64_t C) {
        const double* a = p + 3 * A;
        const double* b = p + 3 * B;
        const double* d = p + 3 * C;
        v += (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) +
              a[2] * (b[0] * d[1] - b[1] * d[0])) /
             6.0;
    };
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::int64_t* t = cb.Conn().As<std::int64_t>();
        const std::size_t k = std::string(cb.Type()) == "triangle" ? 3 : 4;
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const std::int64_t* f = t + k * c;
            tet(f[0], f[1], f[2]);
            if (k == 4)
                tet(f[0], f[2], f[3]);
            for (std::size_t i = 0; i < k; ++i)
                ++uses[{f[i], f[(i + 1) % k]}];
        }
    }
    rClosed = true;
    for (const auto& [e, n] : uses)
        rClosed = rClosed && n == 1 && uses.count({e[1], e[0]}) == 1;
    return v;
}

}  // namespace

TEST(OffsetSurface, ACubeGrowsByTheDistanceOnEveryFace) {
    const OffsetSurfaceResult r = offset_surface(ofs_cube(), {0.1, 0});
    const double* p = r.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < 24; ++i)
        EXPECT_NEAR(std::abs(p[i] - 0.5), 0.6, 1e-12) << i;
    EXPECT_EQ(r.mNumFolded, 0);
    EXPECT_NEAR(r.mMaxMove, 0.1 * std::sqrt(3.0), 1e-12);
    // Inward by 0.1: the cube shrinks to [0.1, 0.9].
    const OffsetSurfaceResult in = offset_surface(ofs_cube(), {-0.1, 0});
    const double* q = in.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < 24; ++i)
        EXPECT_NEAR(std::abs(q[i] - 0.5), 0.4, 1e-12) << i;
}

TEST(OffsetSurface, ASphereKeepsTheDistanceToItsFaces) {
    const Mesh s = icosa_sphere(1.0, 3);
    const OffsetSurfaceResult r = offset_surface(s, {0.05, 0});
    const double* p = r.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
        const double rad = std::hypot(p[3 * i], p[3 * i + 1], p[3 * i + 2]);
        EXPECT_GE(rad, 1.05 - 1e-12);
        EXPECT_LE(rad, 1.05 * 1.01);
    }
    EXPECT_EQ(r.mNumNotVisible, 0);
    EXPECT_TRUE(r.mWarnings.empty());
    EXPECT_THROW(offset_surface(s, {std::nan(""), 0}), std::invalid_argument);
}

TEST(OffsetSurface, AThickenedHemisphereIsClosedAndReadsItsWall) {
    ClipOptions co;
    co.mKeepBelow = false;
    co.mCap = false;
    const Mesh hemi = clip(icosa_sphere(1.0, 4), co).mMesh;
    for (const ThickenSide side : {ThickenSide::Inward, ThickenSide::Outward, ThickenSide::Both}) {
        ThickenOptions o;
        o.mThickness = 0.05;
        o.mSide = side;
        const OffsetSurfaceResult r = thicken(hemi, o);
        bool closed = false;
        const double v = ofs_volume(r.mMesh, closed);
        EXPECT_TRUE(closed);
        EXPECT_GT(v, 0.0);
        EXPECT_GT(r.mNumRimQuads, 0);
        // The two half shells' volume, (2 pi / 3) (R^3 - r^3).
        const double hi = side == ThickenSide::Outward  ? 1.05
                          : side == ThickenSide::Inward ? 1.0
                                                        : 1.025;
        const double want = 2.0 * M_PI / 3.0 * (hi * hi * hi - std::pow(hi - 0.05, 3));
        EXPECT_NEAR(v / want, 1.0, 0.02);
        ThicknessOptions to;
        const ThicknessResult t = thickness(r.mMesh, to);
        EXPECT_NEAR(t.mP2, 0.05, 0.0005) << int(side);
    }
    ThickenOptions bad;
    EXPECT_THROW(thicken(hemi, bad), std::invalid_argument);
    EXPECT_THROW(thicken_side_from_name("up"), std::invalid_argument);
    EXPECT_EQ(thicken_side_from_name("both"), ThickenSide::Both);
}

TEST(OffsetSurface, HollowGivesAClosedWallAndRefusesAnOpenSurface) {
    const OffsetSurfaceResult r = hollow(ofs_cube(), 0.1);
    bool closed = false;
    const double v = ofs_volume(r.mMesh, closed);
    EXPECT_TRUE(closed);
    EXPECT_NEAR(v, 1.0 - 0.8 * 0.8 * 0.8, 1e-12);
    EXPECT_EQ(r.mNumRimQuads, 0);
    ClipOptions co;
    co.mCap = false;
    EXPECT_THROW(hollow(clip(icosa_sphere(1.0, 2), co).mMesh, 0.1), std::invalid_argument);
}

TEST(OffsetSurface, BoundaryLayerSharesTheDirections) {
    // The first layer of a boundary layer and the offset by its height agree.
    const Mesh s = icosa_sphere(1.0, 2);
    BoundaryLayerOptions bo;
    bo.mFirstHeight = 0.02;
    bo.mNumLayers = 1;
    const BoundaryLayerResult bl = boundary_layer(s, bo);
    const OffsetSurfaceResult r = offset_surface(s, {0.02, 0});
    const double* a = r.mMesh.Points().As<double>();
    const double* b = bl.mMesh.Points().As<double>();
    const std::size_t np = s.NumPoints();
    ASSERT_GE(bl.mMesh.NumPoints(), 2 * np);
    double worst = 0.0;
    for (std::size_t i = 0; i < 3 * np; ++i)
        worst = std::max(worst, std::abs(a[i] - b[3 * np + i]));
    EXPECT_LT(worst, 1e-12);
}
