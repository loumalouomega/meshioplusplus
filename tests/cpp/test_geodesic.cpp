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
// Geodesic distance by the heat method (roadmap 7.3.1).

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/geodesic.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

// An icosahedral sphere projected onto the unit sphere.
Mesh geo_sphere(int Subdivisions) {
    Mesh m = icosa_sphere(1.0, Subdivisions);
    NDArray pts = NDArray::Uninit(DType::Float64, {m.NumPoints(), 3});
    for (std::size_t i = 0; i < m.NumPoints(); ++i) {
        double p[3], r = 0.0;
        for (std::size_t k = 0; k < 3; ++k) {
            p[k] = m.Points().As<double>()[i * 3 + k];
            r += p[k] * p[k];
        }
        r = std::sqrt(r);
        for (std::size_t k = 0; k < 3; ++k)
            pts.As<double>()[i * 3 + k] = p[k] / r;
    }
    m.AssignPoints(std::move(pts));
    return m;
}

}  // namespace

TEST(Geodesic, PoleToEquatorOnASphere) {
    double previous = 1.0;
    for (int sub : {3, 4, 5}) {
        const Mesh m = geo_sphere(sub);
        std::size_t pole = 0;
        for (std::size_t i = 0; i < m.NumPoints(); ++i)
            if (m.Points().As<double>()[i * 3 + 2] > m.Points().As<double>()[pole * 3 + 2])
                pole = i;
        GeodesicOptions o;
        o.mSources = {static_cast<std::int64_t>(pole)};
        const GeodesicResult r = geodesic_distance(m, o);
        EXPECT_TRUE(r.mConverged);
        const double* d = r.mMesh.PointData(kGeodesicDistanceName).As<double>();
        double worst = 0.0;
        for (std::size_t i = 0; i < m.NumPoints(); ++i)
            if (std::fabs(m.Points().As<double>()[i * 3 + 2]) < 1e-9)
                worst = std::max(worst, std::fabs(d[i] - M_PI / 2));
        EXPECT_LT(worst, previous);
        previous = worst;
        EXPECT_EQ(d[pole], 0.0);
    }
    EXPECT_LT(previous, 1.3e-2);
}

TEST(Geodesic, ScopeAndArguments) {
    GeodesicOptions o;
    EXPECT_THROW(geodesic_distance(mt::tri_mesh(), o), std::invalid_argument);  // no source
    o.mSources = {0};
    EXPECT_THROW(geodesic_distance(mt::tet_mesh(), o), std::invalid_argument);
    EXPECT_THROW(geodesic_distance(mt::quad_mesh(), o), std::invalid_argument);
    o.mSources = {99};
    EXPECT_THROW(geodesic_distance(mt::tri_mesh(), o), std::invalid_argument);
    o.mSources = {0};
    o.mTimeFactor = -1.0;
    EXPECT_THROW(geodesic_distance(mt::tri_mesh(), o), std::invalid_argument);
    o.mTimeFactor = 1.0;
    const GeodesicResult r = geodesic_distance(mt::tri_mesh(), o);
    EXPECT_EQ(r.mNumSources, 1);
    EXPECT_EQ(r.mNumUnreached, 0);
    // The opposite corner of the unit square is about sqrt(2) away.
    EXPECT_NEAR(r.mMesh.PointData(kGeodesicDistanceName).As<double>()[2], std::sqrt(2.0), 0.3);
}
