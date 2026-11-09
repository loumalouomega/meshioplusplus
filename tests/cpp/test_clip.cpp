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
// clip: cut by a plane, keep a side, cap it (roadmap 7.8.3).

#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/clip.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/slice.hpp"

using namespace meshioplusplus;

namespace {

// The signed volume of a closed triangle surface, and whether it is closed
// and consistently wound (each directed edge once, its reverse once).
double cl_volume(const Mesh& rMesh, bool& rClosed) {
    const double* p = rMesh.Points().As<double>();
    const std::int64_t* t = rMesh.Cells(0).Conn().As<std::int64_t>();
    std::map<std::array<std::int64_t, 2>, int> uses;
    double v = 0.0;
    for (std::size_t c = 0; c < rMesh.Cells(0).NumCells(); ++c) {
        const double* a = p + 3 * t[3 * c];
        const double* b = p + 3 * t[3 * c + 1];
        const double* d = p + 3 * t[3 * c + 2];
        v += (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0]) +
              a[2] * (b[0] * d[1] - b[1] * d[0])) /
             6.0;
        for (int k = 0; k < 3; ++k)
            ++uses[{t[3 * c + k], t[3 * c + (k + 1) % 3]}];
    }
    rClosed = true;
    for (const auto& [e, n] : uses)
        rClosed = rClosed && n == 1 && uses.count({e[1], e[0]}) == 1;
    return v;
}

}  // namespace

TEST(Clip, AClosedSphereCutAtItsEquatorIsAClosedHemisphere) {
    const Mesh s = icosa_sphere(1.0, 4);
    bool closed = false;
    const double whole = cl_volume(s, closed);
    ASSERT_TRUE(closed);
    for (const bool below : {true, false}) {
        ClipOptions o;
        o.mOrigin = {0.0, 0.0, 0.013};  // off the vertices
        o.mKeepBelow = below;
        const ClipResult r = clip(s, o);
        EXPECT_EQ(r.mNumLoops, 1);
        EXPECT_GT(r.mNumCapTriangles, 0);
        const double half = cl_volume(r.mMesh, closed);
        EXPECT_TRUE(closed);
        EXPECT_GT(half, 0.0);  // outward: positive volume
        // A sphere cap of height 1 -+ 0.013, scaled by the tessellation.
        const double h = below ? 1.013 : 0.987;
        const double cap = 3.14159265358979323846 * h * h * (3.0 - h) / 3.0;
        EXPECT_NEAR(half / whole, cap / (4.0 * 3.14159265358979323846 / 3.0), 2e-3);
    }
}

TEST(Clip, TheBoundaryIsTheSliceAndOpenCutsAreNotCapped) {
    const Mesh s = icosa_sphere(1.0, 3);
    ClipOptions o;
    o.mOrigin = {0.1, 0.0, 0.0};
    o.mNormal = {1.0, 0.0, 0.0};
    o.mCap = false;
    const ClipResult r = clip(s, o);
    SliceOptions so;
    so.mOrigin = o.mOrigin;
    so.mNormal = o.mNormal;
    const Mesh sec = slice(s, so);
    // Every section point is a point of the clipped surface, bit for bit.
    const double* ps = sec.Points().As<double>();
    const double* pc = r.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < sec.NumPoints(); ++i) {
        bool found = false;
        for (std::size_t j = 0; j < r.mMesh.NumPoints() && !found; ++j)
            found = ps[3 * i] == pc[3 * j] && ps[3 * i + 1] == pc[3 * j + 1] &&
                    ps[3 * i + 2] == pc[3 * j + 2];
        EXPECT_TRUE(found);
    }
    // An open disk cut through: one open chain, nothing capped.
    const Mesh d = disk(1.0, 6, {0.0, 0.0, 0.0}, "triangle");
    ClipOptions od;
    od.mOrigin = {0.05, 0.0, 0.0};
    od.mNormal = {1.0, 0.0, 0.0};
    const ClipResult rd = clip(d, od);
    EXPECT_EQ(rd.mNumOpenChains, 1);
    EXPECT_EQ(rd.mNumCapTriangles, 0);
    o.mNormal = {0.0, 0.0, 0.0};
    EXPECT_THROW(clip(s, o), std::invalid_argument);
    EXPECT_THROW(clip(box({{1, 1, 1}}), ClipOptions{}), std::invalid_argument);
}
