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
// Elements, dimensions and GD&T (roadmap 7.2); datum alignment (7.3.4).

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/measure.hpp"
#include "meshioplusplus/region.hpp"

using namespace meshioplusplus;

namespace {

constexpr double kMsTestPi = 3.14159265358979323846;

// A block of quads, one cell region per face: bottom (z = 0), top, xmin,
// xmax, ymin, ymax. The xmax face leans by `Tilt` radians about its bottom
// edge; the whole block is then rotated about z by `Turn` and moved by `Shift`.
Mesh ms_block(double H, double Tilt, double Turn = 0.0,
              const std::array<double, 3>& rShift = {0.0, 0.0, 0.0}) {
    const double d = H * std::tan(Tilt);
    std::vector<std::vector<double>> pts = {{0, 0, 0}, {2, 0, 0},     {2, 3, 0},     {0, 3, 0},
                                            {0, 0, H}, {2 + d, 0, H}, {2 + d, 3, H}, {0, 3, H}};
    for (auto& r_p : pts) {
        const double x = r_p[0], y = r_p[1];
        r_p[0] = std::cos(Turn) * x - std::sin(Turn) * y + rShift[0];
        r_p[1] = std::sin(Turn) * x + std::cos(Turn) * y + rShift[1];
        r_p[2] += rShift[2];
    }
    Mesh m = mt::make_mesh(
        pts, "quad",
        {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4}, {3, 7, 6, 2}});
    const char* names[] = {"bottom", "top", "xmin", "xmax", "ymin", "ymax"};
    for (std::int64_t c = 0; c < 6; ++c) {
        NDArray e = NDArray::Uninit(DType::Int64, {std::size_t{1}});
        e.As<std::int64_t>()[0] = c;
        m.AddRegion(Region(names[c], RegionKind::Cell, std::move(e)));
    }
    return m;
}

std::string ms_planes_spec(const std::string& rRest) {
    return R"({"elements": [
        {"name": "A", "fit": "plane", "region": "bottom"},
        {"name": "B", "fit": "plane", "region": "xmax"},
        {"name": "C", "fit": "plane", "region": "ymin"},
        {"name": "T", "fit": "plane", "region": "top"}],)" +
           rRest + "}";
}

}  // namespace

TEST(Measure, ATiltedSideIsOutOfPerpendicularity) {
    const double tilt = 0.05 * kMsTestPi / 180.0;
    const Mesh block = ms_block(20.0, tilt);
    const MeasureResult r = measure(block, ms_planes_spec(R"(
        "tolerances": [{"name": "perp", "type": "perpendicularity", "feature": "B",
                        "datum": "A", "limit": 0.01},
                       {"name": "par", "type": "parallelism", "feature": "T", "datum": "A",
                        "limit": 0.01},
                       {"name": "flat", "type": "flatness", "feature": "A"}],
        "dimensions": [{"name": "h", "type": "distance", "between": ["A", "T"],
                        "nominal": 20, "tolerance": 0.001},
                       {"name": "ang", "type": "angle", "between": ["A", "B"]}])"));
    ASSERT_EQ(r.mTolerances.size(), 3u);
    EXPECT_NEAR(r.mTolerances[0].mValue, 20.0 * std::tan(tilt), 1e-9);
    EXPECT_FALSE(r.mTolerances[0].mPass);
    EXPECT_NEAR(r.mTolerances[1].mValue, 0.0, 1e-12);
    EXPECT_TRUE(r.mTolerances[1].mPass);
    EXPECT_NEAR(r.mTolerances[2].mValue, 0.0, 1e-12);
    EXPECT_NEAR(r.mDimensions[0].mValue, 20.0, 1e-12);
    EXPECT_TRUE(r.mDimensions[0].mPass);
    EXPECT_NEAR(r.mDimensions[1].mValue, 90.0 - 0.05, 1e-9);
    EXPECT_EQ(r.mNumFail, 1);
    EXPECT_TRUE(r.mMesh.HasFieldData(kMeasureTableName));
    EXPECT_NE(r.mMesh.FindRegion("measure:B", RegionKind::Cell), Mesh::npos);
    EXPECT_EQ(measure_table(r.mMesh).size(), 4u);
}

TEST(Measure, ConstructionsAndAngularity) {
    const Mesh block = ms_block(4.0, 0.0);
    const MeasureResult r = measure(block, ms_planes_spec(R"(
        "tolerances": [{"type": "angularity", "feature": "B", "datum": "A", "basic_angle": 90},
                       {"type": "angularity", "feature": "T", "datum": "A", "basic_angle": 0}])"));
    EXPECT_NEAR(r.mTolerances[0].mValue, 0.0, 1e-12);
    EXPECT_NEAR(r.mTolerances[1].mValue, 0.0, 1e-12);
    // The constructions, on their own spec.
    const MeasureResult c = measure(block, R"({"elements": [
        {"name": "A", "fit": "plane", "region": "bottom"},
        {"name": "B", "fit": "plane", "region": "xmax"},
        {"name": "C", "fit": "plane", "region": "ymin"},
        {"name": "T", "fit": "plane", "region": "top"},
        {"name": "edge", "construct": "intersection", "of": ["A", "B"]},
        {"name": "corner", "construct": "intersection", "of": ["edge", "C"]},
        {"name": "mid", "construct": "mid_plane", "of": ["A", "T"]},
        {"name": "origin", "point": [0, 0, 0]},
        {"name": "p", "construct": "projection", "of": ["origin", "T"]}],
      "dimensions": [{"type": "distance", "between": ["corner", "origin"]},
                     {"type": "distance", "between": ["mid", "A"]},
                     {"type": "distance", "between": ["p", "origin"]},
                     {"type": "distance", "between": ["edge", "origin"]}]})");
    EXPECT_NEAR(c.mDimensions[0].mValue, 2.0, 1e-12);  // the corner (2, 0, 0)
    EXPECT_NEAR(c.mDimensions[1].mValue, 2.0, 1e-12);
    EXPECT_NEAR(c.mDimensions[2].mValue, 4.0, 1e-12);
    EXPECT_NEAR(c.mDimensions[3].mValue, 2.0, 1e-12);
}

TEST(Measure, Refusals) {
    const Mesh block = ms_block(4.0, 0.0);
    EXPECT_THROW(measure(block, "not json"), std::invalid_argument);
    EXPECT_THROW(measure(block, R"({"elements": [{"name": "A", "fit": "plane"}]})"),
                 std::invalid_argument);  // no region
    EXPECT_THROW(measure(block, R"({"elements": [{"name": "A", "fit": "plane", "region": "bottom"},
                                                  {"name": "A", "point": [0, 0, 0]}]})"),
                 std::invalid_argument);  // duplicate
    EXPECT_THROW(measure(block, R"({"elements": [{"name": "A", "fit": "plane", "region": "bottom"},
                                                  {"name": "B", "construct": "axis", "of": "A"}]})"),
                 std::invalid_argument);  // a plane has no axis
    EXPECT_THROW(measure(block, R"({"dimensions": [{"type": "distance", "between": ["X", "Y"]}]})"),
                 std::invalid_argument);
}

TEST(AlignDatums, ABlockLiesOnTheCoordinatePlanes) {
    const Mesh block = ms_block(4.0, 0.0, 0.7, {5.0, -3.0, 2.0});
    const MeasureResult r = measure(block, ms_planes_spec(R"("dimensions": [])"));
    const AlignDatumsResult a = align_datums(r.mMesh, "A", "B", "C");
    const NDArray& p = a.mMesh.Points();
    // bottom (z = 0): points 0-3; xmax: points 1, 2, 5, 6 at x = 0 after the
    // move (B sets X); ymin: points 0, 1, 4, 5 at y = 0.
    for (const std::size_t i : {0u, 1u, 2u, 3u})
        EXPECT_NEAR(p.As<double>()[3 * i + 2], 0.0, 1e-12);
    for (const std::size_t i : {1u, 2u, 5u, 6u})
        EXPECT_NEAR(p.As<double>()[3 * i], 0.0, 1e-12);
    for (const std::size_t i : {0u, 1u, 4u, 5u})
        EXPECT_NEAR(p.As<double>()[3 * i + 1], 0.0, 1e-12);
    const std::vector<MeasureElement> t = measure_table(a.mMesh);
    EXPECT_NEAR(std::fabs(t[0].mDir[2]), 1.0, 1e-12);                      // A's normal is now z
    EXPECT_THROW(align_datums(r.mMesh, "A", "T"), std::invalid_argument);  // parallel
    EXPECT_THROW(align_datums(block, "A", "B"), std::invalid_argument);    // no table
}

TEST(Measure, AnElementIsADeviationReference) {
    const Mesh block = ms_block(4.0, 0.0);
    const MeasureResult r = measure(block, ms_planes_spec(R"("dimensions": [])"));
    const auto [type, params] = measure_element_primitive(measure_element(r.mMesh, "T"));
    EXPECT_EQ(type, 1);
    EXPECT_NEAR(std::fabs(params[5]), 1.0, 1e-12);
    EXPECT_THROW(measure_element(r.mMesh, "nope"), std::invalid_argument);
}
