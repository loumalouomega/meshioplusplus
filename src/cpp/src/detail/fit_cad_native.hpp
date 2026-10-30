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
/**
 * @file detail/fit_cad_native.hpp
 * @brief The native `fit_cad` backend's helpers that the OpenCASCADE backend
 *        (providers/occt/fit_cad_occt.cpp) builds on.
 *
 * Core-internal and not installed. Both backends share the freeform networks,
 * the refinement and the report; the OpenCASCADE backend lives in the
 * `meshioplusplus-occt` add-on in a split build, so these are exported from
 * the core and valid only within one release (the add-on pins it exactly).
 */
#pragma once

// System includes
#include <cstddef>
#include <cstdint>
#include <vector>

// Project includes
#include "closed_facets.hpp"
#include "patch_flatten.hpp"
#include "patch_network.hpp"
#include "primitive_constraints.hpp"
#include "quad_layout.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/operations/refine_patches.hpp"

namespace meshioplusplus::detail {

/// What the refinement and the primitive constraints did, for the report.
struct FcnRefined {
    bool mRan = false;
    RefinePatchesResult mResult;
    PrimitiveConstraintReport mConstraints;
};

/// A freeform region's patch network (network and auto modes).
struct FcnRegionNetwork {
    std::vector<std::size_t> mItems;
    std::vector<PatchVec> mPts;
    std::vector<PatchTri> mTris;
    std::vector<std::size_t> mTriFacet;
    std::vector<std::int64_t> mGlobal;
    QuadLayout mLayout;
    PatchNetwork mNet;
    std::vector<std::vector<std::size_t>> mPatchFacets;
};

/// The share of the skin's area still facets, from the triangulation.
MESHIOPLUSPLUS_API double native_facet_area_fraction(const std::vector<double>& rXyz,
                                                     const std::vector<std::int64_t>& rTri,
                                                     const std::vector<std::int64_t>& rType);

/// The per-face table, as the field data `fit_cad` writes.
MESHIOPLUSPLUS_API NDArray native_face_table(const std::vector<FitCadFaceRecord>& rFaces);

/// The mesh the fit runs on: segmented and refined as the options ask. Returns
/// `&rMesh` or `&rOwned`.
MESHIOPLUSPLUS_API const Mesh* native_prepare(const Mesh& rMesh, Mesh& rOwned,
                                              const FitCadOptions& rOptions,
                                              FcnRefined& rRefined);

/// What the refinement did, into the report.
MESHIOPLUSPLUS_API void native_report_refined(const FcnRefined& rRefined, const Mesh& rWork,
                                              FitCadResult& rResult);

/// The items grouped into regions, each laid out and fitted as one network.
MESHIOPLUSPLUS_API std::vector<FcnRegionNetwork> native_networks(
    const ClosedFacets& rFc, const std::vector<std::vector<std::size_t>>& rItems,
    const FitCadOptions& rOptions, double Tol);

/// A patch's facets cut into K pieces by their normals.
MESHIOPLUSPLUS_API std::vector<std::vector<std::size_t>> native_resegment(
    const ClosedFacets& rFc, const std::vector<std::size_t>& rFacets, std::size_t K,
    bool OnlyAcrossCreases = false);

/// The OpenCASCADE backend (providers/occt/fit_cad_occt.cpp).
FitCadResult occt_fit_cad(const Mesh& rMesh, const std::string& rStepPath,
                          const FitCadOptions& rOptions);

}  // namespace meshioplusplus::detail
