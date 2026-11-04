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
// The surface-boolean spike (roadmap 7.8, doc/spikes.md): the expansion
// arithmetic it builds on, exact volumes on boxes with integer-halved
// coordinates, the inclusion-exclusion identities on spheres, the trivial
// nested and disjoint cases, the named coplanar failure and determinism.

// System includes
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/primitives.hpp"

// Project includes (private)
#include "../../src/cpp/src/detail/expansion.hpp"
#include "../../src/cpp/src/detail/mesh_boolean.hpp"

namespace mio = meshioplusplus;
using mio::detail::BooleanOp;

namespace {

/// The box [Lo, Hi] as 12 outward triangles.
mio::Mesh tmb_box(const std::array<double, 3>& rLo, const std::array<double, 3>& rHi) {
    mio::Mesh m;
    mio::NDArray p = mio::NDArray::Uninit(mio::DType::Float64, {std::size_t{8}, std::size_t{3}});
    for (int i = 0; i < 8; ++i)
        for (int k = 0; k < 3; ++k)
            p.As<double>()[3 * i + k] = (i >> k) & 1 ? rHi[k] : rLo[k];
    m.AssignPoints(std::move(p));
    // Vertex i has bit k set when its coordinate k is high.
    const std::int64_t t[12][3] = {{0, 2, 3}, {0, 3, 1}, {4, 5, 7}, {4, 7, 6},
                                   {0, 1, 5}, {0, 5, 4}, {2, 6, 7}, {2, 7, 3},
                                   {0, 4, 6}, {0, 6, 2}, {1, 3, 7}, {1, 7, 5}};
    mio::NDArray c = mio::NDArray::Uninit(mio::DType::Int64, {std::size_t{12}, std::size_t{3}});
    for (int i = 0; i < 12; ++i)
        for (int k = 0; k < 3; ++k)
            c.As<std::int64_t>()[3 * i + k] = t[i][k];
    m.AddCellBlock("triangle", std::move(c));
    return m;
}

double tmb_volume(const mio::Mesh& rA, const mio::Mesh& rB, BooleanOp Op) {
    const auto r = mio::detail::mesh_boolean_spike(rA, rB, Op);
    EXPECT_TRUE(r.mOk) << r.mFailure;
    return r.mOk ? mio::detail::closed_surface_volume(r.mMesh) : -1.0;
}

}  // namespace

TEST(MeshBooleanExpansion, SumsAndProductsAreExact) {
    // 1 + 2^-60 is not a double; the expansion keeps it.
    const auto e =
        mio::detail::exp_sum(mio::detail::exp_of(1.0), mio::detail::exp_of(std::ldexp(1.0, -60)));
    EXPECT_EQ(mio::detail::exp_sign(mio::detail::exp_sub(e, mio::detail::exp_of(1.0))), 1);
    // (1 + 2^-30)^2 - 1 - 2^-29 = 2^-60 exactly.
    const auto a =
        mio::detail::exp_sum(mio::detail::exp_of(1.0), mio::detail::exp_of(std::ldexp(1.0, -30)));
    auto sq = mio::detail::exp_sub(mio::detail::exp_mul(a, a), mio::detail::exp_of(1.0));
    sq = mio::detail::exp_sub(sq, mio::detail::exp_of(std::ldexp(1.0, -29)));
    EXPECT_EQ(mio::detail::exp_estimate(sq), std::ldexp(1.0, -60));
    EXPECT_EQ(mio::detail::exp_sign(mio::detail::exp_sub(sq, sq)), 0);
    EXPECT_EQ(mio::detail::exp_sign(mio::detail::exp_scale(mio::detail::exp_diff(0.1, 0.3), 3.0)),
              -1);
}

TEST(MeshBoolean, OffsetBoxesHaveTheExactVolumes) {
    // The unit box and the box moved by a half along each axis share an
    // eighth of a unit; their faces cross but never lie in one plane.
    const mio::Mesh a = tmb_box({0, 0, 0}, {1, 1, 1});
    const mio::Mesh b = tmb_box({0.5, 0.5, 0.5}, {1.5, 1.5, 1.5});
    EXPECT_NEAR(tmb_volume(a, b, BooleanOp::Union), 1.875, 1e-12);
    EXPECT_NEAR(tmb_volume(a, b, BooleanOp::Intersection), 0.125, 1e-12);
    EXPECT_NEAR(tmb_volume(a, b, BooleanOp::Difference), 0.875, 1e-12);
    const auto r = mio::detail::mesh_boolean_spike(a, b, BooleanOp::Union);
    ASSERT_TRUE(r.mOk) << r.mFailure;
    EXPECT_GT(r.mReport.mNumIntersectingPairs, 0);
    EXPECT_EQ(r.mReport.mNumLabelConflicts, 0);
    EXPECT_EQ(r.mReport.mNumCoplanarPairs, 0);
}

TEST(MeshBoolean, OffsetSpheresMeetTheInclusionExclusionIdentities) {
    const mio::Mesh a = mio::icosa_sphere(1.0, 3);
    const mio::Mesh b = mio::icosa_sphere(0.8, 2, {0.7, 0.31, 0.17});
    const double va = mio::detail::closed_surface_volume(a);
    const double vb = mio::detail::closed_surface_volume(b);
    const double vu = tmb_volume(a, b, BooleanOp::Union);
    const double vi = tmb_volume(a, b, BooleanOp::Intersection);
    const double vd = tmb_volume(a, b, BooleanOp::Difference);
    EXPECT_GT(vi, 0.0);
    EXPECT_LT(vi, vb);
    EXPECT_NEAR(vu + vi, va + vb, 1e-9 * (va + vb));
    EXPECT_NEAR(vd, va - vi, 1e-9 * va);
}

TEST(MeshBoolean, NestedAndDisjointSurfacesSelectWholeSurfaces) {
    const mio::Mesh outer = mio::icosa_sphere(1.0, 2);
    const mio::Mesh inner = mio::icosa_sphere(0.4, 2);
    const mio::Mesh far = mio::icosa_sphere(1.0, 2, {3.0, 0.0, 0.0});
    const double vo = mio::detail::closed_surface_volume(outer);
    const double vn = mio::detail::closed_surface_volume(inner);
    EXPECT_NEAR(tmb_volume(outer, inner, BooleanOp::Union), vo, 1e-12);
    EXPECT_NEAR(tmb_volume(outer, inner, BooleanOp::Intersection), vn, 1e-12);
    EXPECT_NEAR(tmb_volume(outer, inner, BooleanOp::Difference), vo - vn, 1e-12);
    EXPECT_NEAR(tmb_volume(outer, far, BooleanOp::Union), 2.0 * vo, 1e-12);
    const auto r = mio::detail::mesh_boolean_spike(outer, far, BooleanOp::Intersection);
    ASSERT_TRUE(r.mOk) << r.mFailure;
    EXPECT_EQ(r.mReport.mNumOutputTriangles, 0);
    EXPECT_NEAR(tmb_volume(outer, far, BooleanOp::Difference), vo, 1e-12);
}

TEST(MeshBoolean, CoplanarContactFailsByName) {
    const auto r = mio::detail::mesh_boolean_spike(tmb_box({0, 0, 0}, {1, 1, 1}),
                                                   tmb_box({1, 0, 0}, {2, 1, 1}), BooleanOp::Union);
    EXPECT_FALSE(r.mOk);
    EXPECT_EQ(r.mFailure, "coplanar contact");
}

TEST(MeshBoolean, OpenInputFailsByName) {
    mio::Mesh open;
    mio::NDArray p = mio::NDArray::Uninit(mio::DType::Float64, {std::size_t{3}, std::size_t{3}});
    const double xyz[9] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    for (int i = 0; i < 9; ++i)
        p.As<double>()[i] = xyz[i];
    open.AssignPoints(std::move(p));
    mio::NDArray c = mio::NDArray::Uninit(mio::DType::Int64, {std::size_t{1}, std::size_t{3}});
    for (int i = 0; i < 3; ++i)
        c.As<std::int64_t>()[i] = i;
    open.AddCellBlock("triangle", std::move(c));
    const auto r =
        mio::detail::mesh_boolean_spike(open, tmb_box({0, 0, 0}, {1, 1, 1}), BooleanOp::Union);
    EXPECT_FALSE(r.mOk);
    EXPECT_NE(r.mFailure.find("not closed"), std::string::npos) << r.mFailure;
}

TEST(MeshBoolean, IsDeterministic) {
    const mio::Mesh a = mio::icosa_sphere(1.0, 2);
    const mio::Mesh b = mio::icosa_sphere(0.9, 2, {0.6, 0.2, 0.1});
    const auto r1 = mio::detail::mesh_boolean_spike(a, b, BooleanOp::Difference);
    const auto r2 = mio::detail::mesh_boolean_spike(a, b, BooleanOp::Difference);
    ASSERT_TRUE(r1.mOk && r2.mOk);
    const auto& p1 = r1.mMesh.Points();
    const auto& p2 = r2.mMesh.Points();
    ASSERT_EQ(p1.Size(), p2.Size());
    for (std::size_t i = 0; i < p1.Size(); ++i)
        ASSERT_EQ(p1.As<double>()[i], p2.As<double>()[i]);
    EXPECT_EQ(r1.mReport.mNumOutputTriangles, r2.mReport.mNumOutputTriangles);
}
