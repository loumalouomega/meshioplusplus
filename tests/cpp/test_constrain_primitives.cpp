// constrain_primitives under control (roadmap 7.7.3).

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/operations/constrain_primitives.hpp"
#include "analytic_surfaces.hpp"
#include "fit_cad_fixtures.hpp"

#include "../../src/cpp/src/detail/primitive_constraints.hpp"

using namespace meshioplusplus;

namespace {

std::vector<double> cp_row(const Mesh& rMesh, std::int64_t Patch) {
    const NDArray& ids = rMesh.FieldData("fit:patch_id");
    const NDArray& a = rMesh.FieldData("fit:patch_params");
    for (std::size_t r = 0; r < ids.Size(); ++r)
        if (ids.As<std::int64_t>()[r] == Patch)
            return std::vector<double>(a.As<double>() + r * 10, a.As<double>() + r * 10 + 10);
    return {};
}

Mesh cp_box() {
    return fcfix::fitted(fcfix::jittered(analytic::box_skin(4), 2e-5), 0.0, 1e-3);
}

}  // namespace

TEST(ConstrainPrimitives, TheDefaultRunIsFitCadsPass) {
    const Mesh box = cp_box();
    const ConstrainPrimitivesResult r = constrain_primitives(box);
    Mesh direct = detail::clone_mesh(box);
    const detail::PrimitiveConstraintReport rep =
        detail::constrain_primitives(direct, r.mFaceTolerance);
    EXPECT_EQ(r.mNumConstrained, rep.mNumConstrained);
    EXPECT_GT(r.mNumConstrained, 0);
    const NDArray& a = r.mMesh.FieldData("fit:patch_params");
    const NDArray& b = direct.FieldData("fit:patch_params");
    ASSERT_EQ(a.Size(), b.Size());
    for (std::size_t i = 0; i < a.Size(); ++i)
        EXPECT_EQ(a.As<double>()[i], b.As<double>()[i]);
    // Every change is reported, each with a residual.
    EXPECT_FALSE(r.mRelations.empty());
    for (const PrimitiveRelationRecord& rec : r.mRelations) {
        EXPECT_GE(rec.mAngleDeg, 0.0);
        EXPECT_LT(rec.mAngleDeg, 1.0);
        EXPECT_FALSE(rec.mForced);
    }
    EXPECT_TRUE(r.mUnmet.empty());
}

TEST(ConstrainPrimitives, AForbiddenRelationLeavesTheNormalAlone) {
    const Mesh box = cp_box();
    const ConstrainPrimitivesResult first = constrain_primitives(box);
    // A plane the default run snapped, and every relation it took.
    std::int64_t target = -1;
    for (const PrimitiveRelationRecord& rec : first.mRelations)
        if (rec.mKind == PrimitiveRelationKind::Perpendicular ||
            rec.mKind == PrimitiveRelationKind::Parallel)
            target = rec.mPatchB;
    ASSERT_GE(target, 0);
    ConstrainPrimitivesOptions o;
    const NDArray& ids = box.FieldData("fit:patch_id");
    for (std::size_t r = 0; r < ids.Size(); ++r) {
        const std::int64_t other = ids.As<std::int64_t>()[r];
        if (other == target)
            continue;
        o.mForbid.push_back({other, target, PrimitiveRelationKind::Perpendicular});
        o.mForbid.push_back({other, target, PrimitiveRelationKind::Parallel});
    }
    const ConstrainPrimitivesResult r = constrain_primitives(box, o);
    EXPECT_EQ(cp_row(r.mMesh, target), cp_row(box, target));
    for (const PrimitiveRelationRecord& rec : r.mRelations)
        EXPECT_NE(rec.mPatchB, target);
}

TEST(ConstrainPrimitives, ToleranceForceAndRefusals) {
    const Mesh box = cp_box();
    // A zero angle tolerance snaps nothing that is not already exact.
    ConstrainPrimitivesOptions tight;
    tight.mAngleToleranceDeg = 0.0;
    EXPECT_LT(constrain_primitives(box, tight).mNumConstrained,
              constrain_primitives(box).mNumConstrained);
    // A forced relation between two parallel planes that cannot hold (two
    // opposite faces made coaxial) is reported unmet.
    const NDArray& ids = box.FieldData("fit:patch_id");
    ConstrainPrimitivesOptions force;
    force.mForce.push_back(
        {ids.As<std::int64_t>()[0], ids.As<std::int64_t>()[1], PrimitiveRelationKind::Coaxial});
    EXPECT_EQ(constrain_primitives(box, force).mUnmet.size(), 1u);
    ConstrainPrimitivesOptions bad;
    bad.mForbid.push_back({9999, ids.As<std::int64_t>()[0], PrimitiveRelationKind::Parallel});
    EXPECT_THROW(constrain_primitives(box, bad), std::invalid_argument);
    EXPECT_THROW(constrain_primitives(analytic::box_skin(2)), std::invalid_argument);
    EXPECT_EQ(primitive_relation_from_name("equal_radius"), PrimitiveRelationKind::EqualRadius);
    EXPECT_STREQ(primitive_relation_name(PrimitiveRelationKind::Perpendicular), "perpendicular");
    EXPECT_THROW(primitive_relation_from_name("tangent"), std::invalid_argument);
}
