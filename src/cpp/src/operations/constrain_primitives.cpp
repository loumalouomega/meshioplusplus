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
// Constrain fitted primitives under control (roadmap 7.7.3). See
// operations/constrain_primitives.hpp.

// System includes
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>

// Project includes
#include "meshioplusplus/operations/constrain_primitives.hpp"
#include "../detail/primitive_constraints.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/segment.hpp"

namespace meshioplusplus {

namespace {

[[noreturn]] void cpr_fail(const std::string& rWhy) {
    throw std::invalid_argument("meshio++: constrain_primitives: " + rWhy);
}

}  // namespace

PrimitiveRelationKind primitive_relation_from_name(const std::string& rName) {
    if (rName == "coaxial")
        return PrimitiveRelationKind::Coaxial;
    if (rName == "parallel")
        return PrimitiveRelationKind::Parallel;
    if (rName == "perpendicular")
        return PrimitiveRelationKind::Perpendicular;
    if (rName == "equal_radius")
        return PrimitiveRelationKind::EqualRadius;
    cpr_fail("a relation is coaxial, parallel, perpendicular or equal_radius, not '" + rName + "'");
}

const char* primitive_relation_name(PrimitiveRelationKind Kind) noexcept {
    switch (Kind) {
        case PrimitiveRelationKind::Coaxial:
            return "coaxial";
        case PrimitiveRelationKind::Parallel:
            return "parallel";
        case PrimitiveRelationKind::Perpendicular:
            return "perpendicular";
        case PrimitiveRelationKind::EqualRadius:
            return "equal_radius";
    }
    return "";
}

ConstrainPrimitivesResult constrain_primitives(const Mesh& rMesh,
                                               const ConstrainPrimitivesOptions& rOptions) {
    for (const char* name : {"fit:patch_id", "fit:patch_type", "fit:patch_params", "fit:patch_rms",
                             "fit:patch_max", "fit:patch_fraction", "fit:patch_tolerance"})
        if (!rMesh.HasFieldData(name))
            cpr_fail(std::string("the mesh has no '") + name +
                     "' field data (run fit_primitives first)");
    if (!rMesh.HasCellData(kSegmentPatchName))
        cpr_fail("the mesh has no segment:patch cell data (run segment first)");
    const ConstrainPrimitivesOptions& o = rOptions;
    for (const double v :
         {o.mFaceTolerance, o.mAngleToleranceDeg, o.mDistanceTolerance, o.mMinFraction})
        if (!std::isfinite(v))
            cpr_fail("a tolerance is not finite");
    if (o.mMinFraction > 1.0)
        cpr_fail("min_fraction is at most 1");
    if (o.mAngleToleranceDeg > 90.0)
        cpr_fail("angle_tolerance is at most 90 degrees");
    const NDArray& ids = rMesh.FieldData("fit:patch_id");
    std::set<std::int64_t> known;
    for (std::size_t r = 0; r < ids.Size(); ++r)
        known.insert(detail::read_int(ids, r));
    for (const auto* p_list : {&o.mForce, &o.mForbid})
        for (const PrimitiveRelation& r_rel : *p_list)
            if (!known.count(r_rel.mPatchA) || !known.count(r_rel.mPatchB) ||
                r_rel.mPatchA == r_rel.mPatchB)
                cpr_fail("relation (" + std::to_string(r_rel.mPatchA) + ", " +
                         std::to_string(r_rel.mPatchB) + ") names no two fitted patches");
    double face_tol = o.mFaceTolerance;
    if (!(face_tol > 0.0)) {
        const NDArray& tols = rMesh.FieldData("fit:patch_tolerance");
        double most = 0.0;
        // fit_cad's own default (fcn_face_tolerance): three times the largest.
        for (std::size_t r = 0; r < tols.Size(); ++r) {
            const double t = detail::read_double(tols, r);
            if (std::isfinite(t))
                most = std::max(most, t);
        }
        face_tol = 3.0 * most;
    }
    if (!(face_tol > 0.0))
        cpr_fail("no fitted patch has a tolerance; give face_tolerance");

    ConstrainPrimitivesResult result;
    result.mMesh = detail::clone_mesh(rMesh);
    detail::PrimitiveConstraintControl control;
    control.mAngleToleranceDeg = o.mAngleToleranceDeg;
    control.mDistanceTolerance = o.mDistanceTolerance;
    control.mMinFraction = o.mMinFraction;
    control.mForce = o.mForce;
    control.mForbid = o.mForbid;
    detail::PrimitiveConstraintReport rep =
        detail::constrain_primitives(result.mMesh, face_tol, control);
    result.mNumConstrained = rep.mNumConstrained;
    result.mNumGroups = rep.mNumGroups;
    result.mRelations = std::move(rep.mRelations);
    result.mUnmet = std::move(rep.mUnmet);
    result.mFaceTolerance = face_tol;
    return result;
}

}  // namespace meshioplusplus
