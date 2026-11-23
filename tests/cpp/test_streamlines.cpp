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
//
//
// Streamlines of a point vector field as a mesh of line segments (roadmap 15.4).

#include <cmath>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/streamlines.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;
using meshioplusplus::detail::read_double;

namespace {

constexpr double kPi = 3.14159265358979323846;

// An n x n grid of unit-spaced cells over [-1, 1]^2, each split into two
// triangles, carrying `field(x, y)` as a 3-component point array.
template <class TField>
Mesh grid(int N, TField Field) {
    std::vector<std::vector<double>> pts;
    std::vector<double> data;
    for (int j = 0; j <= N; ++j)
        for (int i = 0; i <= N; ++i) {
            const double x = -1.0 + 2.0 * i / N;
            const double y = -1.0 + 2.0 * j / N;
            pts.push_back({x, y, 0.0});
            const std::array<double, 3> v = Field(x, y);
            data.insert(data.end(), v.begin(), v.end());
        }
    std::vector<std::vector<std::int64_t>> tris;
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
            const std::int64_t a = j * (N + 1) + i;
            tris.push_back({a, a + 1, a + N + 2});
            tris.push_back({a, a + N + 2, a + N + 1});
        }
    Mesh m = mt::make_mesh(pts, "triangle", tris);
    m.AddPointData("u", mt::data_array(data, 3));
    return m;
}

Mesh rotation() {
    return grid(20, [](double x, double y) { return std::array<double, 3>{-y, x, 0.0}; });
}

StreamlinesOptions seeded(const std::vector<double>& Seed) {
    StreamlinesOptions o;
    o.mArray = "u";
    o.mSeedPoints = Seed;
    return o;
}

double radius(const Mesh& rMesh, std::size_t I) {
    const double x = read_double(rMesh.Points(), I * 3);
    const double y = read_double(rMesh.Points(), I * 3 + 1);
    return std::sqrt(x * x + y * y);
}

}  // namespace

TEST(Streamlines, ARotationTracesCirclesOfConstantRadius) {
    StreamlinesOptions o = seeded({0.5, 0.0, 0.0});
    o.mDirection = StreamlineDirection::Forward;
    o.mMaxLength = 1.0;  // a fraction of the 2 sqrt 2 diagonal: about 2.83 of the 3.14 circle
    const StreamlinesResult r = streamlines(rotation(), o);
    ASSERT_EQ(r.mNumLines, 1);
    ASSERT_GT(r.mMesh.NumPoints(), 50u);
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i)
        EXPECT_NEAR(radius(r.mMesh, i), 0.5, 2e-3) << i;
    // The speed is |v| = r and the time is the length over it: t = angle.
    const std::size_t last = r.mMesh.NumPoints() - 1;
    EXPECT_NEAR(read_double(r.mMesh.PointData(kStreamlineSpeedName), last), 0.5, 2e-3);
    const double angle = std::atan2(read_double(r.mMesh.Points(), last * 3 + 1),
                                    read_double(r.mMesh.Points(), last * 3));
    // (the angle wraps: compare modulo a full turn)
    const double turned = read_double(r.mMesh.PointData(kStreamlineTimeName), last) - angle;
    EXPECT_NEAR(std::remainder(turned, 2.0 * kPi), 0.0, 5e-3);
    EXPECT_GT(read_double(r.mMesh.PointData(kStreamlineTimeName), last), 5.0);
    EXPECT_EQ(r.mEndForward[0], static_cast<std::int32_t>(StreamlineStop::MaxLength));
    EXPECT_EQ(r.mDim, 2);
}

TEST(Streamlines, BothDirectionsGiveNegativeAndPositiveTime) {
    StreamlinesOptions o = seeded({0.5, 0.0, 0.0});
    o.mMaxLength = 0.3;
    const StreamlinesResult r = streamlines(rotation(), o);
    ASSERT_EQ(r.mNumLines, 1);
    const NDArray& t = r.mMesh.PointData(kStreamlineTimeName);
    const std::size_t n = r.mMesh.NumPoints();
    EXPECT_LT(read_double(t, 0), 0.0);
    EXPECT_GT(read_double(t, n - 1), 0.0);
    for (std::size_t i = 1; i < n; ++i)
        EXPECT_GT(read_double(t, i), read_double(t, i - 1));
    // One segment less than the points, each naming its seed.
    ASSERT_EQ(r.mMesh.NumCellBlocks(), 1u);
    EXPECT_EQ(r.mMesh.Cells(0).NumCells(), n - 1);
}

TEST(Streamlines, TheAdaptiveIntegratorAgreesWithRK4) {
    StreamlinesOptions a = seeded({0.4, 0.1, 0.0});
    a.mDirection = StreamlineDirection::Forward;
    a.mMaxLength = 0.4;
    StreamlinesOptions b = a;
    b.mIntegrator = StreamlineIntegrator::RK45;
    const StreamlinesResult ra = streamlines(rotation(), a);
    const StreamlinesResult rb = streamlines(rotation(), b);
    ASSERT_EQ(ra.mNumLines, 1);
    ASSERT_EQ(rb.mNumLines, 1);
    // The end points are at the same arc length, so they sit on the same circle
    // to within the first-order difference in where the length limit cuts.
    const std::size_t la = ra.mMesh.NumPoints() - 1;
    const std::size_t lb = rb.mMesh.NumPoints() - 1;
    EXPECT_NEAR(radius(ra.mMesh, la), radius(rb.mMesh, lb), 2e-3);
    EXPECT_NEAR(radius(rb.mMesh, lb), std::sqrt(0.4 * 0.4 + 0.1 * 0.1), 2e-3);
}

TEST(Streamlines, TheReasonsALineStops) {
    // A uniform field: the line leaves the mesh, in both directions.
    const Mesh uniform =
        grid(8, [](double, double) { return std::array<double, 3>{1.0, 0.0, 0.0}; });
    StreamlinesOptions o = seeded({0.0, 0.0, 0.0});
    o.mMaxLength = 5.0;
    StreamlinesResult r = streamlines(uniform, o);
    ASSERT_EQ(r.mNumLines, 1);
    EXPECT_EQ(r.mEndBackward[0], static_cast<std::int32_t>(StreamlineStop::LeftDomain));
    EXPECT_EQ(r.mEndForward[0], static_cast<std::int32_t>(StreamlineStop::LeftDomain));
    // The step limit and the length limit.
    o.mMaxSteps = 3;
    r = streamlines(uniform, o);
    EXPECT_EQ(r.mEndForward[0], static_cast<std::int32_t>(StreamlineStop::MaxSteps));
    EXPECT_EQ(r.mMesh.NumPoints(), 7u);  // the seed and 3 steps each way
    o.mMaxSteps = 0;
    o.mMaxLength = 0.1;
    r = streamlines(uniform, o);
    EXPECT_EQ(r.mEndForward[0], static_cast<std::int32_t>(StreamlineStop::MaxLength));
    // A field that dies away: stagnation (the speed falls below the terminal speed).
    const Mesh slow =
        grid(8, [](double x, double) { return std::array<double, 3>{0.5 * (1.0 - x), 0.0, 0.0}; });
    StreamlinesOptions s = seeded({-0.5, 0.0, 0.0});
    s.mDirection = StreamlineDirection::Forward;
    s.mMaxLength = 5.0;
    s.mTerminalSpeed = 0.05;
    r = streamlines(slow, s);
    ASSERT_EQ(r.mNumLines, 1);
    EXPECT_EQ(r.mEndForward[0], static_cast<std::int32_t>(StreamlineStop::Stagnation));
}

TEST(Streamlines, ASeedOutsideTheMeshMakesNoLine) {
    StreamlinesOptions o = seeded({0.5, 0.0, 0.0, 5.0, 5.0, 0.0});
    o.mMaxLength = 0.1;
    const StreamlinesResult r = streamlines(rotation(), o);
    EXPECT_EQ(r.mNumLines, 1);
    EXPECT_EQ(r.mNumSeedsOutside, 1);
    EXPECT_EQ(r.mSeedIds[0], 0);
}

TEST(Streamlines, ACountOfSeedsAndTheFieldRidesAlongUnderItsOwnName) {
    StreamlinesOptions o;
    o.mArray = "u";
    o.mSeeds = 12;
    o.mMaxLength = 0.2;
    const StreamlinesResult r = streamlines(rotation(), o);
    EXPECT_GT(r.mNumLines, 6);
    EXPECT_LE(r.mNumLines, 12);
    EXPECT_TRUE(r.mMesh.HasPointData("u"));
    EXPECT_EQ(data_num_components(r.mMesh.PointData("u")), 3u);
    EXPECT_TRUE(r.mMesh.HasPointData(kStreamlineSeedIdName));
    EXPECT_TRUE(r.mMesh.HasCellData(kStreamlineSeedIdName));
}

TEST(Streamlines, AVolumeIsTracedInItsVolume) {
    // The unit cube as six tetrahedra around a diagonal, with a uniform +z field.
    Mesh m = mt::make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "tetra",
        {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}});
    std::vector<double> v;
    for (int i = 0; i < 8; ++i)
        v.insert(v.end(), {0.0, 0.0, 2.0});
    m.AddPointData("u", mt::data_array(v, 3));
    StreamlinesOptions o = seeded({0.3, 0.4, 0.5});
    o.mMaxLength = 2.0;
    const StreamlinesResult r = streamlines(m, o);
    ASSERT_EQ(r.mNumLines, 1);
    EXPECT_EQ(r.mDim, 3);
    for (std::size_t i = 0; i < r.mMesh.NumPoints(); ++i) {
        EXPECT_NEAR(read_double(r.mMesh.Points(), i * 3), 0.3, 1e-9);
        EXPECT_NEAR(read_double(r.mMesh.Points(), i * 3 + 1), 0.4, 1e-9);
    }
    // It ends at the faces z = 0 and z = 1, with time = distance / 2.
    const std::size_t n = r.mMesh.NumPoints();
    EXPECT_NEAR(
        read_double(r.mMesh.PointData(kStreamlineTimeName), n - 1) -
            read_double(r.mMesh.PointData(kStreamlineTimeName), 0),
        (read_double(r.mMesh.Points(), (n - 1) * 3 + 2) - read_double(r.mMesh.Points(), 2)) / 2.0,
        1e-9);
}

TEST(Streamlines, TheResultDoesNotDependOnTheRun) {
    StreamlinesOptions o;
    o.mArray = "u";
    o.mSeeds = 30;
    o.mIntegrator = StreamlineIntegrator::RK45;
    const Mesh m = rotation();
    const StreamlinesResult a = streamlines(m, o);
    const StreamlinesResult b = streamlines(m, o);
    ASSERT_EQ(a.mMesh.NumPoints(), b.mMesh.NumPoints());
    for (std::size_t i = 0; i < a.mMesh.NumPoints() * 3; ++i)
        EXPECT_EQ(read_double(a.mMesh.Points(), i), read_double(b.mMesh.Points(), i));
}

TEST(Streamlines, Refusals) {
    const Mesh m = rotation();
    StreamlinesOptions o;
    EXPECT_THROW(streamlines(m, o), std::invalid_argument);  // no array
    o.mArray = "nothing";
    EXPECT_THROW(streamlines(m, o), std::invalid_argument);
    o.mArray = "u";
    o.mSeedPoints = {1.0, 2.0};
    EXPECT_THROW(streamlines(m, o), std::invalid_argument);  // not 3 per seed
    o.mSeedPoints.clear();
    o.mMaxLength = 0.0;
    EXPECT_THROW(streamlines(m, o), std::invalid_argument);
    // A scalar field.
    Mesh s = rotation();
    s.AddPointData("t", mt::data_array(std::vector<double>(s.NumPoints(), 1.0)));
    StreamlinesOptions t;
    t.mArray = "t";
    EXPECT_THROW(streamlines(s, t), std::invalid_argument);
    // A mesh with no cell the tracer takes.
    Mesh lines = mt::make_mesh({{0, 0, 0}, {1, 0, 0}}, "line", {{0, 1}});
    lines.AddPointData("u", mt::data_array({1, 0, 0, 1, 0, 0}, 3));
    StreamlinesOptions l;
    l.mArray = "u";
    EXPECT_THROW(streamlines(lines, l), std::invalid_argument);
}
