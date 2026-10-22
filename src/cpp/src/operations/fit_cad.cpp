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
// fit_cad. See operations/fit_cad.hpp for the contract.

// System includes
#include <string>

// System includes
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

// Project includes
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/refine_patches.hpp"
#include "meshioplusplus/operations/segment.hpp"

// Project includes (private, not installed)
#include "../detail/closed_facets.hpp"
#include "../detail/cut_graph.hpp"
#include "../detail/freeform_adaptive.hpp"
#include "../detail/patch_flatten.hpp"
#include "../detail/patch_network.hpp"
#include "../detail/quad_layout.hpp"
#include "../detail/step_brep.hpp"
#include "../detail/surface_rings.hpp"

// System includes
#include <array>
#include <functional>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#ifdef MESHIOPLUSPLUS_HAS_OCCT
// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/segment.hpp"

// Project includes (private, not installed)
#include "../detail/closed_facets.hpp"
#include "../detail/library_lock.hpp"
#include "../detail/lsq_fit.hpp"
#include "../detail/feature_fit.hpp"
#include "../detail/occt_brep.hpp"
#include "../detail/bspline_fit.hpp"
#include "../detail/cut_graph.hpp"
#include "../detail/freeform_adaptive.hpp"
#include "../detail/param_disc.hpp"
#include "../detail/patch_flatten.hpp"
#include "../detail/primitive_geom.hpp"

// External includes (OpenCASCADE)
#include <APIHeaderSection_MakeHeader.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepGProp.hxx>
#include <BRepTools.hxx>
#include <BRepTools_History.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GProp_GProps.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <TColStd_HArray1OfReal.hxx>
#include <TColgp_HArray1OfPnt2d.hxx>
#include <TopoDS_Iterator.hxx>
#include <gp_Pnt2d.hxx>
#include <GeomAPI_IntSS.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <Approx_ParametrizationType.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_BoundedSurface.hxx>
#include <GeomLib.hxx>
#include <GeomFill_Pipe.hxx>
#include <Geom_OffsetSurface.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomLProp_SLProps.hxx>
#include <GeomPlate_BuildPlateSurface.hxx>
#include <GeomPlate_MakeApprox.hxx>
#include <GeomPlate_PointConstraint.hxx>
#include <GeomPlate_Surface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <Poly_Triangulation.hxx>
#include <Precision.hxx>
#include <STEPControl_Writer.hxx>
#include <StepBasic_Product.hxx>
#include <StepData_StepModel.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <ShapeFix_Solid.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TCollection_HAsciiString.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Ax3.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Sphere.hxx>
#include <gp_Torus.hxx>
#endif

namespace meshioplusplus {

bool fit_cad_has_occt() noexcept {
#ifdef MESHIOPLUSPLUS_HAS_OCCT
    return true;
#else
    return false;
#endif
}

const char* fit_cad_occt_version() noexcept {
#ifdef MESHIOPLUSPLUS_HAS_OCCT
    return OCC_VERSION_COMPLETE;
#else
    return "";
#endif
}

namespace {

constexpr const char* kFcnPrefix = "fit_cad: ";

// The share of a triangulated surface's area on facet faces (`fit_cad:type`
// 0): `rTri` three point ids per triangle, `rType` per triangle. Summed in
// triangle order, so the value does not depend on the thread count.
double fcn_facet_area_fraction(const std::vector<double>& rXyz,
                               const std::vector<std::int64_t>& rTri,
                               const std::vector<std::int64_t>& rType) {
    double total = 0.0, facet = 0.0;
    for (std::size_t t = 0; t < rType.size(); ++t) {
        const double* a = &rXyz[static_cast<std::size_t>(rTri[3 * t]) * 3];
        const double* b = &rXyz[static_cast<std::size_t>(rTri[3 * t + 1]) * 3];
        const double* c = &rXyz[static_cast<std::size_t>(rTri[3 * t + 2]) * 3];
        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                             u[0] * v[1] - u[1] * v[0]};
        const double area = 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        total += area;
        if (rType[t] == 0)
            facet += area;
    }
    return total > 0.0 ? facet / total : 0.0;
}

// `fit_cad:face_table`: one row per record, the `FitCadFaceRecord` columns in
// declaration order as float64.
NDArray fcn_face_table(const std::vector<FitCadFaceRecord>& rFaces) {
    constexpr std::size_t kCols = 13;
    NDArray a = NDArray::Uninit(DType::Float64, {rFaces.size(), kCols});
    double* p = a.As<double>();
    for (std::size_t i = 0; i < rFaces.size(); ++i) {
        const FitCadFaceRecord& f = rFaces[i];
        const double row[kCols] = {static_cast<double>(f.mPatch),
                                   static_cast<double>(f.mFace),
                                   static_cast<double>(f.mKind),
                                   static_cast<double>(f.mParamMethod),
                                   static_cast<double>(f.mDemotionReason),
                                   static_cast<double>(f.mNumControl),
                                   static_cast<double>(f.mIterations),
                                   static_cast<double>(f.mRefinements),
                                   static_cast<double>(f.mEscalationSteps),
                                   f.mMaxDist,
                                   f.mRmsDist,
                                   f.mReverseMaxDist,
                                   f.mMaxG1AngleDeg};
        std::copy(row, row + kCols, p + i * kCols);
    }
    return a;
}

// The annotated skin: segmented by the surfaces its facets lie on and fitted
// when `segment:patch` is missing, fitted when the `fit:patch_*` table is.
// Returns a pointer to the mesh to work on, owned by rOwned or the input.
const Mesh* fcn_annotate(const Mesh& rMesh, Mesh& rOwned) {
    if (rMesh.HasCellData(kSegmentPatchName)) {
        if (rMesh.HasFieldData("fit:patch_id"))
            return &rMesh;
        log::info("fit_cad: no fit:patch_* table; running fit_primitives");
        rOwned = fit_primitives(rMesh).mMesh;
        return &rOwned;
    }
    log::info(
        "fit_cad: no segment:patch labels; running segment (method 'primitive') and "
        "fit_primitives");
    SegmentOptions seg;
    seg.mMethod = SegmentMethod::Primitive;
    Mesh segmented = segment(rMesh, seg).mMesh;
    // Fit at the tolerance the patches were grown at: an absolute one, the
    // segmentation's share of the diagonal.
    double lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
        lo[k] = std::numeric_limits<double>::infinity();
        hi[k] = -lo[k];
    }
    const std::vector<double> xyz = detail::surface_xyz(rMesh);
    for (std::size_t i = 0; i + 2 < xyz.size(); i += 3)
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], xyz[i + k]);
            hi[k] = std::max(hi[k], xyz[i + k]);
        }
    double diag = 0.0;
    for (int k = 0; k < 3; ++k)
        diag += hi[k] > lo[k] ? (hi[k] - lo[k]) * (hi[k] - lo[k]) : 0.0;
    diag = std::sqrt(diag);
    FitPrimitivesOptions fit;
    if (diag > 0.0) {
        fit.mTolerance = seg.mTolerance * diag;
        fit.mRelative = false;
    }
    rOwned = fit_primitives(segmented, fit).mMesh;
    return &rOwned;
}

// What the refinement did, for the report.
struct FcnRefined {
    bool mRan = false;
    RefinePatchesResult mResult;
};

// The annotated skin, refined by `refine_patches` when asked.
const Mesh* fcn_prepare(const Mesh& rMesh, Mesh& rOwned, const FitCadOptions& rOptions,
                        FcnRefined& rRefined) {
    const Mesh* work = fcn_annotate(rMesh, rOwned);
    if (!rOptions.mRefineSegmentation)
        return work;
    rRefined.mResult = refine_patches(*work, rOptions.mRefine);
    rRefined.mRan = true;
    log::info("fit_cad: refined the patches: {} -> {} of the area explained by primitives",
              rRefined.mResult.mInitialExplainedAreaFraction,
              rRefined.mResult.mExplainedAreaFraction);
    rOwned = std::move(rRefined.mResult.mMesh);
    return &rOwned;
}

// The refinement's report and its `refine:*` field data on the result.
void fcn_report_refined(const FcnRefined& rRefined, const Mesh& rWork, FitCadResult& rResult) {
    if (!rRefined.mRan)
        return;
    const RefinePatchesResult& r = rRefined.mResult;
    rResult.mRefined = true;
    rResult.mRefinementKept = true;
    rResult.mRefineIterations = r.mNumIterations;
    rResult.mNumRefineSplits = r.mNumSplits;
    rResult.mNumRefineMerges = r.mNumMerges;
    rResult.mInitialExplainedAreaFraction = r.mInitialExplainedAreaFraction;
    rResult.mExplainedAreaFraction = r.mExplainedAreaFraction;
    for (const char* name :
         {kRefinePatchIdName, kRefinePatchProxyName, kRefinePatchIterationsName,
          kRefinePatchRmsName, kRefinePatchMaxName, kRefinePatchCellsName,
          kRefinePatchDemotionName, kRefineExplainedAreaName, kRefineHausdorffName})
        if (rWork.HasFieldData(name))
            rResult.mMesh.AddFieldData(name, rWork.FieldData(name));
}

// The kernel-free backend: the fitted patches as analytic faces and the rest
// as coplanar polygons, written by the B-rep `step` writer.
// What the native backend's adaptive freeform pass (v17.12.0) built: the
// faces for the writer, the mesh to write (its unexplained patches split into
// discs where that was needed), and a record per patch it tried.
struct FcnFreeform {
    std::vector<detail::StepFreeformFace> mFaces;
    std::map<std::int64_t, FitCadFaceRecord> mRecords;  // by label
    Mesh mMesh;  // with the pieces' labels and fit-table rows, when split
    bool mSplit = false;
    std::int64_t mNumSplit = 0;
    std::int64_t mNumResegmented = 0;  // regions the re-segmentation rung cut
    // The patch network's (network and auto modes): its patches' labels, its
    // shared edges and extraordinary vertices, and its largest G1 angle.
    std::set<std::int64_t> mNetworkLabels;
    std::vector<detail::StepEdgeCurve> mEdgeCurves;
    std::int64_t mNumNetworkEdges = 0;
    std::int64_t mNumExtraordinary = 0;
    double mMaxG1AngleDeg = 0.0;
};

// The tolerance the native writer will apply: as `step_brep_write` takes it.
double fcn_tolerance(const Mesh& rMesh, const detail::ClosedFacets& rFc, double Requested) {
    if (Requested > 0.0)
        return Requested;
    double largest = 0.0;
    if (rMesh.HasFieldData("fit:patch_type") && rMesh.HasFieldData("fit:patch_tolerance") &&
        rMesh.HasFieldData("fit:patch_params")) {
        const NDArray& types = rMesh.FieldData("fit:patch_type");
        const NDArray& tols = rMesh.FieldData("fit:patch_tolerance");
        const NDArray& params = rMesh.FieldData("fit:patch_params");
        const std::size_t np = types.Size() > 0 ? params.Size() / types.Size() : 0;
        for (std::size_t r = 0; r < types.Size() && r < tols.Size(); ++r) {
            const std::int64_t t = detail::read_int(types, r);
            bool finite = true;
            for (std::size_t k = 0; k < np; ++k)
                finite = finite && std::isfinite(detail::read_double(params, r * np + k));
            const double pt = detail::read_double(tols, r);
            if (t >= 1 && t <= 5 && finite && std::isfinite(pt))
                largest = std::max(largest, pt);
        }
    }
    if (largest > 0.0)
        return 3.0 * largest;
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (std::int64_t v : rFc.mNodes)
        for (int j = 0; j < 3; ++j) {
            lo[j] = std::min(lo[j], rFc.mXyz[static_cast<std::size_t>(v) * 3 + j]);
            hi[j] = std::max(hi[j], rFc.mXyz[static_cast<std::size_t>(v) * 3 + j]);
        }
    return 1e-3 * std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                            (hi[2] - lo[2]) * (hi[2] - lo[2]));
}

// Facets as local triangles (polygons fanned), points renumbered in the order
// they are first used; with `pGlobal`, each local point's skin point.
void fcn_local_triangles(const detail::ClosedFacets& rFc, const std::vector<std::size_t>& rFacets,
                         std::vector<detail::PatchVec>& rPts, std::vector<detail::PatchTri>& rTris,
                         std::vector<std::size_t>& rTriFacet, std::vector<std::int64_t>* pGlobal) {
    std::unordered_map<std::int64_t, std::size_t> local;
    const auto vid = [&](std::int64_t V) {
        auto ins = local.emplace(V, rPts.size());
        if (ins.second) {
            const std::size_t v = static_cast<std::size_t>(V);
            rPts.push_back({rFc.mXyz[v * 3], rFc.mXyz[v * 3 + 1], rFc.mXyz[v * 3 + 2]});
            if (pGlobal)
                pGlobal->push_back(V);
        }
        return ins.first->second;
    };
    for (std::size_t f : rFacets) {
        const std::int64_t b = rFc.mStart[f], n = rFc.mStart[f + 1] - b;
        const std::size_t a0 = vid(rFc.mNodes[static_cast<std::size_t>(b)]);
        for (std::int64_t k = 1; k + 1 < n; ++k) {
            rTris.push_back({a0, vid(rFc.mNodes[static_cast<std::size_t>(b + k)]),
                             vid(rFc.mNodes[static_cast<std::size_t>(b + k + 1)])});
            rTriFacet.push_back(f);
        }
    }
}

// A face's sense: its surface normal against the facets', area-weighted;
// `rU`, `rV` the local points' parameters on it.
bool fcn_sense(const detail::ClosedFacets& rFc, const detail::BsplineSurface& rS,
               const std::vector<detail::PatchTri>& rTris,
               const std::vector<std::size_t>& rTriFacet, const std::vector<double>& rU,
               const std::vector<double>& rV) {
    double agree = 0.0;
    for (std::size_t t = 0; t < rTris.size(); ++t) {
        const std::size_t f = rTriFacet[t];
        double u = 0.0, v = 0.0;
        for (std::size_t c : rTris[t]) {
            u += rU[c] / 3.0;
            v += rV[c] / 3.0;
        }
        detail::BsplineEval e;
        detail::bspline_eval(rS, u, v, 1, e);
        const double n[3] = {e.mSu[1] * e.mSv[2] - e.mSu[2] * e.mSv[1],
                             e.mSu[2] * e.mSv[0] - e.mSu[0] * e.mSv[2],
                             e.mSu[0] * e.mSv[1] - e.mSu[1] * e.mSv[0]};
        agree += n[0] * rFc.mNormals[f * 3] + n[1] * rFc.mNormals[f * 3 + 1] +
                 n[2] * rFc.mNormals[f * 3 + 2];
    }
    return agree >= 0.0;
}

// A freeform region's patch network (network and auto modes): the items it
// joins, its facets as local triangles, its quad layout (`mLayout.mOk` false:
// none, and no network), its fitted network, and per layout patch its facets
// (each facet with its first triangle's patch).
struct FcnRegionNetwork {
    std::vector<std::size_t> mItems;
    std::vector<detail::PatchVec> mPts;
    std::vector<detail::PatchTri> mTris;
    std::vector<std::size_t> mTriFacet;
    std::vector<std::int64_t> mGlobal;
    detail::QuadLayout mLayout;
    detail::PatchNetwork mNet;
    std::vector<std::vector<std::size_t>> mPatchFacets;
};

// The items (patches, as facet lists) grouped into regions -- items sharing a
// facet edge -- and each laid out and fitted as one network, in the order of
// each region's first item. Free of OpenCASCADE: both backends build on it.
std::vector<FcnRegionNetwork> fcn_networks(const detail::ClosedFacets& rFc,
                                           const std::vector<std::vector<std::size_t>>& rItems,
                                           const FitCadOptions& rOptions, double Tol) {
    double diag = 0.0;
    {
        double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
        for (std::int64_t v : rFc.mNodes)
            for (int j = 0; j < 3; ++j) {
                lo[j] = std::min(lo[j], rFc.mXyz[static_cast<std::size_t>(v) * 3 + j]);
                hi[j] = std::max(hi[j], rFc.mXyz[static_cast<std::size_t>(v) * 3 + j]);
            }
        diag = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                         (hi[2] - lo[2]) * (hi[2] - lo[2]));
    }
    std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> edge_item;
    std::vector<std::size_t> parent(rItems.size());
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    const auto find = [&](std::size_t A) {
        while (parent[A] != A)
            A = parent[A] = parent[parent[A]];
        return A;
    };
    for (std::size_t i = 0; i < rItems.size(); ++i)
        for (std::size_t f : rItems[i]) {
            const std::int64_t b = rFc.mStart[f], n = rFc.mStart[f + 1] - b;
            for (std::int64_t k = 0; k < n; ++k) {
                const std::int64_t a = rFc.mNodes[static_cast<std::size_t>(b + k)];
                const std::int64_t c = rFc.mNodes[static_cast<std::size_t>(b + (k + 1) % n)];
                const auto key = std::make_pair(std::min(a, c), std::max(a, c));
                const auto ins = edge_item.emplace(key, i);
                if (!ins.second) {
                    const std::size_t x = find(ins.first->second), y = find(i);
                    if (x != y)
                        parent[std::max(x, y)] = std::min(x, y);
                }
            }
        }
    std::map<std::size_t, std::vector<std::size_t>> regions;
    for (std::size_t i = 0; i < rItems.size(); ++i)
        regions[find(i)].push_back(i);
    std::vector<FcnRegionNetwork> out;
    for (const auto& kv : regions) {
        FcnRegionNetwork reg;
        reg.mItems = kv.second;
        std::vector<std::size_t> facets;
        for (std::size_t i : kv.second)
            facets.insert(facets.end(), rItems[i].begin(), rItems[i].end());
        std::sort(facets.begin(), facets.end());
        fcn_local_triangles(rFc, facets, reg.mPts, reg.mTris, reg.mTriFacet, &reg.mGlobal);
        detail::QuadLayoutOptions lo;
        lo.mTargetSize = rOptions.mLayoutTargetPatchSize * diag;
        lo.mMaxValenceDeviation = static_cast<int>(rOptions.mLayoutMaxValenceDeviation);
        lo.mSnapToFeatures = rOptions.mLayoutSnapToFeatures;
        reg.mLayout = detail::quad_layout(reg.mPts, reg.mTris, lo);
        if (!reg.mLayout.mOk) {
            log::debug("fit_cad: no quad layout of a {}-triangle freeform region (failure {})",
                       reg.mTris.size(), static_cast<int>(reg.mLayout.mFailure));
            out.push_back(std::move(reg));
            continue;
        }
        detail::PatchNetworkOptions no;
        no.mTolerance = Tol;
        no.mFreeformTolerance = rOptions.mFreeformTolerance;
        no.mMaxPolesPerDirection = static_cast<std::size_t>(rOptions.mMaxFreeformPoles);
        no.mMaxTotalPoles = static_cast<std::size_t>(rOptions.mFreeformMaxControlPoints);
        no.mG1ToleranceDeg = rOptions.mFreeformG1ToleranceDeg;
        no.mExtraordinary = rOptions.mExtraordinaryVertex;
        no.mValidate = rOptions.mValidateBySurfaceSampling;
        reg.mNet = detail::patch_network_fit(reg.mPts, reg.mTris, reg.mLayout, no);
        // The tangent-plane treatment of the extraordinary vertices kept only
        // where it does better there than leaving them alone (`g0`) and fits
        // as much: next to a vertex where the patches' corners meet at very
        // unequal angles, its constraint can bend the edges between them.
        if (no.mExtraordinary != FitCadExtraordinaryVertex::G0 && reg.mNet.mNumExtraordinary > 0) {
            detail::PatchNetworkOptions g0 = no;
            g0.mExtraordinary = FitCadExtraordinaryVertex::G0;
            detail::PatchNetwork plain =
                detail::patch_network_fit(reg.mPts, reg.mTris, reg.mLayout, g0);
            const auto fitted = [&](const detail::PatchNetwork& rN) {
                std::size_t n = 0;
                for (std::size_t i = 0; i < reg.mLayout.mPatches.size(); ++i)
                    n += rN.mPatches[i].mOk ? reg.mLayout.mPatches[i].mTris.size() : 0;
                return n;
            };
            if (plain.mMaxG1AngleEvDeg < reg.mNet.mMaxG1AngleEvDeg &&
                fitted(plain) >= fitted(reg.mNet)) {
                log::debug(
                    "fit_cad: the extraordinary vertices of a {}-triangle region left G0: "
                    "{} degrees next to them instead of {}",
                    reg.mTris.size(), plain.mMaxG1AngleEvDeg, reg.mNet.mMaxG1AngleEvDeg);
                reg.mNet = std::move(plain);
            }
        }
        log::debug(
            "fit_cad: patch network of a {}-triangle freeform region: {} patch(es), {} "
            "edge(s), {} extraordinary vertex(es), G1 {} degrees",
            reg.mTris.size(), reg.mLayout.mPatches.size(), reg.mNet.mNumEdges,
            reg.mNet.mNumExtraordinary, reg.mNet.mMaxG1AngleDeg);
        reg.mPatchFacets.assign(reg.mLayout.mPatches.size(), {});
        for (std::size_t t = 0; t < reg.mTris.size(); ++t)
            if (t == 0 || reg.mTriFacet[t] != reg.mTriFacet[t - 1])
                reg.mPatchFacets[static_cast<std::size_t>(reg.mLayout.mPatchOfTri[t])].push_back(
                    reg.mTriFacet[t]);
        out.push_back(std::move(reg));
    }
    return out;
}

// The record of an accepted network patch.
FitCadFaceRecord fcn_network_record(std::int64_t Label, const detail::NetworkPatch& rP) {
    FitCadFaceRecord rec;
    rec.mPatch = Label;
    rec.mParamMethod = FitCadParamMethod::Square;
    rec.mNumControl = static_cast<std::int64_t>(rP.mSurface.NumPoles());
    rec.mIterations = rP.mIterations;
    rec.mRefinements = rP.mRefinements;
    rec.mMaxDist = rP.mMaxDist;
    rec.mRmsDist = rP.mRmsDist;
    rec.mReverseMaxDist = rP.mReverseMaxDist;
    rec.mMaxG1AngleDeg = rP.mMaxG1AngleDeg;
    return rec;
}

// The re-segmentation rung: a region the adaptive faces and the network left
// as facets cut into pieces. Where edges sharper than 30 degrees already cut
// it, the pieces between them; otherwise `K` pieces by variational shape
// approximation on the L2,1 normal metric (Cohen-Steiner, Alliez and Desbrun
// 2004, the energy refine_patches' Lloyd step uses): proxies seeded farthest
// apart in normal, grown together through one priority queue that crosses an
// edge (in steps of ten degrees) only once every gentler one is taken, so the
// pieces meet at the sharpest ridge between their seeds, re-seeded at the
// facet nearest their proxy, for at most ten Lloyd steps, then the cuts
// smoothed. A cell fanned into several facets stays whole. Empty when no cut
// leaves every piece at least four facets.
// `OnlyAcrossCreases`: none for a region with fewer than three sharp edges
// on curves that close or end on its boundary.
std::vector<std::vector<std::size_t>> fcn_resegment(const detail::ClosedFacets& rFc,
                                                    const std::vector<std::size_t>& rFacets,
                                                    std::size_t K, bool OnlyAcrossCreases = false) {
    const std::size_t n = rFacets.size();
    if (K < 2 || n < 4 * K)
        return {};
    const auto normal = [&](std::size_t I) { return &rFc.mNormals[rFacets[I] * 3]; };
    std::vector<double> area(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t f = rFacets[i];
        const std::int64_t b = rFc.mStart[f], m = rFc.mStart[f + 1] - b;
        const double* p0 =
            &rFc.mXyz[static_cast<std::size_t>(rFc.mNodes[static_cast<std::size_t>(b)]) * 3];
        for (std::int64_t k = 1; k + 1 < m; ++k) {
            const double* p1 =
                &rFc.mXyz[static_cast<std::size_t>(rFc.mNodes[static_cast<std::size_t>(b + k)]) *
                          3];
            const double* p2 = &rFc.mXyz[static_cast<std::size_t>(
                                             rFc.mNodes[static_cast<std::size_t>(b + k + 1)]) *
                                         3];
            const double u[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
            const double v[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
            const double c[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                                 u[0] * v[1] - u[1] * v[0]};
            area[i] += 0.5 * std::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
        }
    }
    // Neighbours across shared edges, with the edge's dihedral angle in
    // tens of degrees (from 3, 30 degrees, it is a crease).
    static constexpr double kCosTens[9] = {
        0.98480775301220806, 0.93969262078590838, 0.86602540378443865,
        0.76604444311897804, 0.64278760968653933, 0.5,
        0.34202014332566873, 0.17364817766693035, 0.0};
    std::vector<std::vector<std::pair<std::size_t, std::uint8_t>>> adj(n);
    std::size_t num_creases = 0;
    {
        std::unordered_map<std::uint64_t, std::size_t> first;
        std::vector<std::pair<std::uint64_t, std::uint64_t>> crease_edges;
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t f = rFacets[i];
            const std::int64_t b = rFc.mStart[f], m = rFc.mStart[f + 1] - b;
            for (std::int64_t k = 0; k < m; ++k) {
                const std::uint64_t a =
                    static_cast<std::uint64_t>(rFc.mNodes[static_cast<std::size_t>(b + k)]);
                const std::uint64_t c = static_cast<std::uint64_t>(
                    rFc.mNodes[static_cast<std::size_t>(b + (k + 1) % m)]);
                const std::uint64_t key = (std::min(a, c) << 32) | std::max(a, c);
                const auto ins = first.emplace(key, i);
                if (ins.second || ins.first->second == i)
                    continue;
                const std::size_t j = ins.first->second;
                const double* ni = normal(i);
                const double* nj = normal(j);
                const double dot = ni[0] * nj[0] + ni[1] * nj[1] + ni[2] * nj[2];
                std::uint8_t tens = 0;
                while (tens < 9 && dot < kCosTens[tens])
                    ++tens;
                adj[i].emplace_back(j, tens);
                adj[j].emplace_back(i, tens);
                ins.first->second = n;  // matched: an inner edge
                if (tens >= 3)
                    crease_edges.emplace_back(a, c);
            }
        }
        // A crease must run somewhere: a group of connected sharp edges that
        // neither reaches the region's boundary nor closes (a coarse or noisy
        // skin's scattered ones) is not counted; one that starts at the
        // boundary and fades out inside (a crease partly under the angle) is.
        std::unordered_set<std::uint64_t> on_boundary;
        for (const auto& kv : first)
            if (kv.second != n) {
                on_boundary.insert(kv.first >> 32);
                on_boundary.insert(kv.first & 0xffffffffu);
            }
        std::unordered_map<std::uint64_t, std::uint64_t> parent;
        const auto find = [&](std::uint64_t V) {
            auto it = parent.emplace(V, V).first;
            while (it->second != it->first) {
                auto up = parent.find(it->second);
                it->second = parent.find(up->second)->first;
                it = up;
            }
            return it->first;
        };
        for (const auto& e : crease_edges) {
            const std::uint64_t ra = find(e.first), rb = find(e.second);
            if (ra != rb)
                parent[std::max(ra, rb)] = std::min(ra, rb);
        }
        std::unordered_map<std::uint64_t, std::size_t> edges_of, verts_of;
        std::unordered_set<std::uint64_t> touches;
        for (const auto& e : crease_edges)
            ++edges_of[find(e.first)];
        for (const auto& kv : parent) {
            const std::uint64_t r = find(kv.first);
            ++verts_of[r];
            if (on_boundary.count(kv.first))
                touches.insert(r);
        }
        for (const auto& kv : edges_of)
            if (touches.count(kv.first) || kv.second >= verts_of[kv.first])  // or a cycle
                num_creases += kv.second;
    }
    // A fanned cell takes its first facet's piece.
    const auto pieces = [&](std::vector<std::size_t> Label, std::size_t NumPieces) {
        std::unordered_map<std::int64_t, std::size_t> of_cell;
        for (std::size_t i = 0; i < n; ++i)
            Label[i] = of_cell.emplace(rFc.mCell[rFacets[i]], Label[i]).first->second;
        std::vector<std::vector<std::size_t>> out(NumPieces);
        for (std::size_t i = 0; i < n; ++i)
            out[Label[i]].push_back(rFacets[i]);
        for (const auto& piece : out)
            if (piece.size() < 4)
                return std::vector<std::vector<std::size_t>>{};
        return out;
    };
    // A region whose creases already cut it: its pieces between them.
    {
        std::vector<std::size_t> comp(n, n);
        std::size_t num = 0;
        for (std::size_t i = 0; i < n; ++i) {
            if (comp[i] != n)
                continue;
            std::vector<std::size_t> stack{i};
            comp[i] = num;
            while (!stack.empty()) {
                const std::size_t a = stack.back();
                stack.pop_back();
                for (const auto& [j, tens] : adj[a])
                    if (tens < 3 && comp[j] == n) {
                        comp[j] = num;
                        stack.push_back(j);
                    }
            }
            ++num;
        }
        // A piece under four facets (a corner between two creases) joins the
        // piece it touches first, then the pieces are numbered again.
        std::vector<std::size_t> size(num, 0);
        for (std::size_t c : comp)
            ++size[c];
        for (bool changed = true; changed;) {
            changed = false;
            for (std::size_t i = 0; i < n; ++i)
                if (size[comp[i]] < 4)
                    for (const auto& nb : adj[i])
                        if (comp[nb.first] != comp[i] && size[comp[nb.first]] >= 4) {
                            const std::size_t from = comp[i], to = comp[nb.first];
                            for (std::size_t& c : comp)
                                if (c == from)
                                    c = to;
                            size[to] += size[from];
                            size[from] = 0;
                            changed = true;
                            break;
                        }
        }
        std::vector<std::size_t> renum(num, num);
        std::size_t kept = 0;
        for (std::size_t& c : comp) {
            if (renum[c] == num)
                renum[c] = kept++;
            c = renum[c];
        }
        num = kept;
        if (num >= 2) {
            auto cut = pieces(std::move(comp), num);
            if (!cut.empty())
                return cut;
        }
    }
    if (OnlyAcrossCreases && num_creases < 3)
        return {};
    const auto error = [&](std::size_t I, const std::array<double, 3>& rN) {
        const double* q = normal(I);
        const double d[3] = {q[0] - rN[0], q[1] - rN[1], q[2] - rN[2]};
        return area[I] * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    };
    // Seeds: farthest from the mean normal, then farthest from those taken.
    std::array<double, 3> mean = {0.0, 0.0, 0.0};
    for (std::size_t i = 0; i < n; ++i)
        for (int j = 0; j < 3; ++j)
            mean[j] += area[i] * normal(i)[j];
    const double len = std::sqrt(mean[0] * mean[0] + mean[1] * mean[1] + mean[2] * mean[2]);
    if (!(len > 0.0))
        return {};
    for (double& c : mean)
        c /= len;
    std::vector<std::array<double, 3>> proxy;
    std::vector<std::size_t> seed;
    std::vector<double> nearest(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        nearest[i] = error(i, mean) / std::max(area[i], 1e-300);
    while (seed.size() < K) {
        const std::size_t s = static_cast<std::size_t>(
            std::max_element(nearest.begin(), nearest.end()) - nearest.begin());
        if (!(nearest[s] > 1e-12))
            return {};
        seed.push_back(s);
        const double* q = normal(s);
        proxy.push_back({q[0], q[1], q[2]});
        for (std::size_t i = 0; i < n; ++i)
            nearest[i] = std::min(nearest[i], error(i, proxy.back()) / std::max(area[i], 1e-300));
    }
    std::vector<std::size_t> label(n);
    for (int step = 0; step < 10; ++step) {
        // Grow: (sharpest edge crossed, error, facet, proxy), smallest first,
        // so the pieces meet at the sharpest ridge between their seeds.
        using Entry = std::tuple<std::uint8_t, double, std::size_t, std::size_t>;
        std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> queue;
        std::vector<std::size_t> grown(n, K);
        for (std::size_t p = 0; p < K; ++p)
            queue.emplace(std::uint8_t{0}, 0.0, seed[p], p);
        while (!queue.empty()) {
            const auto [sharpest, e, i, p] = queue.top();
            queue.pop();
            if (grown[i] != K)
                continue;
            grown[i] = p;
            for (const auto& [j, tens] : adj[i])
                if (grown[j] == K)
                    queue.emplace(std::max(sharpest, tens), error(j, proxy[p]), j, p);
        }
        // A facet no seed reaches (a region in pieces) joins the first proxy.
        for (std::size_t& g : grown)
            if (g == K)
                g = 0;
        const bool same = step > 0 && grown == label;
        label = std::move(grown);
        if (same)
            break;
        // Refit the proxies, and re-seed each at its facet nearest to it.
        for (std::size_t p = 0; p < K; ++p) {
            std::array<double, 3> s = {0.0, 0.0, 0.0};
            for (std::size_t i = 0; i < n; ++i)
                if (label[i] == p)
                    for (int j = 0; j < 3; ++j)
                        s[j] += area[i] * normal(i)[j];
            const double l = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
            if (!(l > 0.0))
                return {};
            for (double& c : s)
                c /= l;
            proxy[p] = s;
            double best = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < n; ++i)
                if (label[i] == p && error(i, s) / std::max(area[i], 1e-300) < best) {
                    best = error(i, s) / std::max(area[i], 1e-300);
                    seed[p] = i;
                }
        }
    }
    // Smooth the cuts: a facet with more sides on one other piece than on its
    // own moves there (a few passes), so no piece keeps a jagged spur.
    for (int pass = 0; pass < 4; ++pass) {
        bool changed = false;
        for (std::size_t i = 0; i < n; ++i) {
            std::vector<std::size_t> votes(K, 0);
            for (const auto& nb : adj[i])
                ++votes[label[nb.first]];
            for (std::size_t p = 0; p < K; ++p)
                if (p != label[i] && votes[p] > votes[label[i]] && votes[p] >= 2) {
                    label[i] = p;
                    changed = true;
                    break;
                }
        }
        if (!changed)
            break;
    }
    return pieces(std::move(label), K);
}

// The adaptive freeform pass of the native backend: every patch the fit table
// leaves unexplained (type 0, or a refinement's freeform or facets verdict)
// fitted by detail/freeform_adaptive; a patch that is not a disc cut into
// discs first (detail/cut_graph), one that is but misses halved, up to
// `mMaxFreeformSplits` times. The pieces get new labels after the largest.
FcnFreeform fcn_freeform(const Mesh& rSource, const detail::ClosedFacets& rFc,
                         const std::vector<std::int64_t>& rCandidates,
                         const FitCadOptions& rOptions, double Tol) {
    FcnFreeform out;
    if (!rSource.HasCellData("segment:patch") || !rSource.HasFieldData("fit:patch_id") ||
        !rSource.HasFieldData("fit:patch_type"))
        return out;
    const NDArray& ids = rSource.FieldData("fit:patch_id");
    const NDArray& types = rSource.FieldData("fit:patch_type");
    const std::vector<std::int64_t> bases = detail::block_bases(rSource);
    std::vector<std::int64_t> cell_label(static_cast<std::size_t>(detail::total_cells(bases)), -1);
    for (std::size_t b = 0; b < rSource.NumCellBlocks(); ++b) {
        const NDArray& a = rSource.CellData("segment:patch", b);
        for (std::size_t c = 0; c < rSource.Cells(b).NumCells() && c < a.Size(); ++c)
            cell_label[static_cast<std::size_t>(bases[b]) + c] = detail::read_int(a, c);
    }
    std::map<std::int64_t, std::vector<std::size_t>> facets_of;  // by label
    for (std::size_t f = 0; f < rFc.NumFacets(); ++f)
        facets_of[cell_label[static_cast<std::size_t>(rFc.mCell[f])]].push_back(f);
    std::int64_t next_label = 0;
    for (std::size_t r = 0; r < ids.Size(); ++r)
        next_label = std::max(next_label, detail::read_int(ids, r) + 1);

    // A work list of (label, facets, halvings so far); pieces join its end.
    struct Item {
        std::int64_t mLabel;
        std::vector<std::size_t> mFacets;
        std::int64_t mDepth;
        std::int64_t mEscalations;
        std::int64_t mCuts = 0;  // early re-segmentation cuts in its line
    };
    std::vector<std::pair<std::int64_t, std::vector<std::size_t>>> pieces;  // new labels
    std::map<std::int64_t, double> piece_tol;

    // The patch network (network and auto modes): the given patches grouped
    // into regions, each laid out and fitted as one network; a patch of it
    // within the tolerance becomes a face under a new label, and the rest
    // come back (as new labels and their facets).
    const auto network =
        [&](const std::vector<std::pair<std::int64_t, std::vector<std::size_t>>>& rItems) {
            std::vector<std::pair<std::int64_t, std::vector<std::size_t>>> left;
            std::vector<std::vector<std::size_t>> item_facets;
            for (const auto& it : rItems)
                item_facets.push_back(it.second);
            for (const FcnRegionNetwork& reg : fcn_networks(rFc, item_facets, rOptions, Tol)) {
                if (!reg.mLayout.mOk) {
                    for (std::size_t i : reg.mItems)
                        left.push_back(rItems[i]);
                    continue;
                }
                const detail::QuadLayout& layout = reg.mLayout;
                const detail::PatchNetwork& net = reg.mNet;
                bool any = false;
                for (std::size_t i = 0; i < layout.mPatches.size(); ++i) {
                    const std::vector<std::size_t>& pf = reg.mPatchFacets[i];
                    if (pf.empty())
                        continue;
                    const detail::NetworkPatch& np = net.mPatches[i];
                    const std::int64_t label = next_label++;
                    pieces.emplace_back(label, pf);
                    if (!np.mOk) {
                        left.emplace_back(label, pf);
                        continue;
                    }
                    any = true;
                    std::vector<detail::PatchTri> ptris;
                    std::vector<std::size_t> pfac;
                    for (std::size_t t : layout.mPatches[i].mTris) {
                        ptris.push_back(reg.mTris[t]);
                        pfac.push_back(reg.mTriFacet[t]);
                    }
                    detail::StepFreeformFace face;
                    face.mLabel = label;
                    face.mSurface = np.mSurface;
                    face.mSense = fcn_sense(rFc, np.mSurface, ptris, pfac, np.mU, np.mV);
                    face.mNetwork = true;
                    out.mFaces.push_back(std::move(face));
                    out.mRecords[label] = fcn_network_record(label, np);
                    out.mNetworkLabels.insert(label);
                }
                if (any) {
                    // Its shared edge curves, for the writer.
                    for (std::size_t c = 0; c < layout.mChains.size(); ++c) {
                        if (net.mChains[c].mPoles.size() < 2)
                            continue;
                        detail::StepEdgeCurve ec;
                        for (std::size_t v : layout.mChains[c].mVerts)
                            ec.mVerts.push_back(reg.mGlobal[v]);
                        ec.mParams = net.mChainParams[c];
                        ec.mCurve = net.mChains[c];
                        out.mEdgeCurves.push_back(std::move(ec));
                    }
                    out.mNumNetworkEdges += static_cast<std::int64_t>(net.mNumEdges);
                    out.mNumExtraordinary += static_cast<std::int64_t>(net.mNumExtraordinary);
                    out.mMaxG1AngleDeg = std::max(out.mMaxG1AngleDeg, net.mMaxG1AngleDeg);
                }
            }
            return left;
        };

    std::vector<Item> work;
    if (rOptions.mFreeformMode == FitCadFreeformMode::Network) {
        std::vector<std::pair<std::int64_t, std::vector<std::size_t>>> items;
        for (std::int64_t label : rCandidates) {
            const auto it = facets_of.find(label);
            if (it != facets_of.end())
                items.emplace_back(label, it->second);
        }
        for (auto& it : network(items))
            work.push_back({it.first, std::move(it.second), 0, 1, rOptions.mResegmentRounds});
    } else {
        for (std::int64_t label : rCandidates) {
            const auto it = facets_of.find(label);
            if (it != facets_of.end())
                work.push_back({label, it->second, 0, 0});
        }
    }
    detail::FreeformAdaptiveOptions o;
    o.mTolerance = Tol;
    o.mFreeformTolerance = rOptions.mFreeformTolerance;
    o.mSmoothing = rOptions.mFreeformSmoothing;
    o.mMaxPoles = static_cast<std::size_t>(rOptions.mMaxFreeformPoles);
    o.mMaxIterations = static_cast<int>(rOptions.mMaxFreeformIterations);
    o.mMaxTotalPoles = static_cast<std::size_t>(rOptions.mFreeformMaxControlPoints);
    o.mParameterization = static_cast<detail::ParamLadder>(rOptions.mFreeformParameterization);
    o.mValidate = rOptions.mValidateBySurfaceSampling;
    o.mThb = rOptions.mLocalRefine == FitCadLocalRefine::Thb;
    // The work list from `From` on: each item a face, or cut into pieces that
    // join the list's end, or facets with its reason.
    const auto run = [&](std::size_t From) {
        for (std::size_t w = From; w < work.size(); ++w) {
            const Item item = work[w];
            // The patch as local triangles (polygons fanned), points renumbered.
            std::vector<detail::PatchVec> pts;
            std::vector<detail::PatchTri> tris;
            std::vector<std::size_t> tri_facet;
            fcn_local_triangles(rFc, item.mFacets, pts, tris, tri_facet, nullptr);
            const detail::FreeformAdaptiveResult fit = detail::freeform_adaptive_fit(pts, tris, o);
            log::debug(
                "fit_cad: adaptive freeform fit of patch {} ({} triangles): {}, {} pole(s), vertex "
                "{} / centroid {} / samples {}, allowance {} (noise {}), {} escalation(s)",
                item.mLabel, tris.size(), fit.mOk ? "accepted" : "rejected",
                fit.mSurface.NumPoles(), fit.mMaxDist, fit.mMaxCentroidDist, fit.mReverseMaxDist,
                fit.mAllowance, fit.mSigma, fit.mEscalations);
            FitCadFaceRecord rec;
            rec.mPatch = item.mLabel;
            rec.mParamMethod = fit.mParam;
            rec.mNumControl = static_cast<std::int64_t>(fit.mNumDofs > 0 ? fit.mNumDofs
                                                                         : fit.mSurface.NumPoles());
            rec.mIterations = fit.mIterations;
            rec.mRefinements = fit.mRefinements;
            rec.mEscalationSteps = item.mEscalations + fit.mEscalations;
            rec.mMaxDist = fit.mMaxDist;
            rec.mRmsDist = fit.mRmsDist;
            rec.mReverseMaxDist = fit.mReverseMaxDist;
            if (fit.mOk) {
                detail::StepFreeformFace face;
                face.mLabel = item.mLabel;
                face.mSurface = fit.mSurface;
                face.mSense = fcn_sense(rFc, fit.mSurface, tris, tri_facet, fit.mU, fit.mV);
                out.mFaces.push_back(std::move(face));
                out.mRecords[item.mLabel] = rec;
                continue;
            }
            // Not a face: cut into discs, else halve, else facets.
            std::size_t num = 0;
            std::vector<std::int64_t> piece;
            bool cut_here = false;
            if (item.mCuts < rOptions.mResegmentRounds) {
                // A patch across sharp edges is cut along them first (at most
                // `mResegmentRounds` times in its line): the re-segmentation
                // rung, where the segmentation left a crease. Its pieces keep
                // the halvings left.
                const auto cut = fcn_resegment(rFc, item.mFacets, 2, true);
                if (!cut.empty()) {
                    std::unordered_map<std::size_t, std::int64_t> of_facet;
                    for (std::size_t k = 0; k < cut.size(); ++k)
                        for (std::size_t f : cut[k])
                            of_facet[f] = static_cast<std::int64_t>(k);
                    piece.resize(tris.size());
                    for (std::size_t t = 0; t < tris.size(); ++t)
                        piece[t] = of_facet[tri_facet[t]];
                    num = cut.size();
                    cut_here = true;
                    ++out.mNumResegmented;
                }
            }
            if (num < 2 && item.mDepth < rOptions.mMaxFreeformSplits && tris.size() >= 32) {
                piece = detail::patch_split_to_discs(pts, tris, 24, num);
                if (num < 2)
                    piece = detail::patch_bisect(pts, tris, num);
            }
            if (num < 2) {
                rec.mDemotionReason = fit.mReason == FitCadDemotionReason::None
                                          ? FitCadDemotionReason::Other
                                          : fit.mReason;
                out.mRecords[item.mLabel] = rec;
                continue;
            }
            std::vector<std::vector<std::size_t>> pf(num);
            for (std::size_t t = 0; t < tris.size(); ++t) {
                if (t > 0 && tri_facet[t] == tri_facet[t - 1])
                    continue;
                pf[static_cast<std::size_t>(piece[t])].push_back(tri_facet[t]);
            }
            for (std::vector<std::size_t>& facets : pf) {
                const std::int64_t label = next_label++;
                work.push_back({label, facets, cut_here ? item.mDepth : item.mDepth + 1,
                                rec.mEscalationSteps + 1, item.mCuts + (cut_here ? 1 : 0)});
                pieces.emplace_back(label, std::move(facets));
            }
            out.mNumSplit += static_cast<std::int64_t>(num);
        }
    };
    run(0);
    // Auto: what the adaptive faces leave goes through the network.
    if (rOptions.mFreeformMode == FitCadFreeformMode::Auto) {
        std::vector<std::pair<std::int64_t, std::vector<std::size_t>>> items;
        for (const Item& item : work) {
            const auto r = out.mRecords.find(item.mLabel);
            if (r != out.mRecords.end() && r->second.mDemotionReason != FitCadDemotionReason::None)
                items.emplace_back(item.mLabel, item.mFacets);
        }
        if (!items.empty()) {
            const std::size_t before = out.mNetworkLabels.size();
            const auto left = network(items);
            if (out.mNetworkLabels.size() > before) {
                // The regions the network took are its pieces now; what it
                // could not fit stays facets, with the adaptive face's reason.
                std::set<std::int64_t> kept;
                for (const auto& l : left)
                    kept.insert(l.first);
                for (const auto& it : items)
                    if (!kept.count(it.first))
                        out.mRecords.erase(it.first);
                for (const auto& l : left)
                    if (!out.mRecords.count(l.first)) {
                        FitCadFaceRecord rec;
                        rec.mPatch = l.first;
                        rec.mEscalationSteps = 1;
                        rec.mDemotionReason = FitCadDemotionReason::HausdorffExceeded;
                        out.mRecords[l.first] = rec;
                    }
            }
        }
    }
    // The re-segmentation rung: up to `mResegmentRounds` times, every region
    // still left as facets is cut by its normals (fcn_resegment, two pieces)
    // and the pieces climb the adaptive rungs again; one that fails is cut
    // again in the next round, never halved.
    for (std::int64_t round = 0; round < rOptions.mResegmentRounds; ++round) {
        std::map<std::int64_t, const std::vector<std::size_t>*> facets_by_label;
        for (const auto& it : facets_of)
            facets_by_label[it.first] = &it.second;
        for (const auto& p : pieces)
            facets_by_label[p.first] = &p.second;
        std::vector<std::tuple<std::int64_t, std::int64_t, std::vector<std::vector<std::size_t>>>>
            cuts;  // (label, escalation steps, pieces)
        for (const auto& [label, rec] : out.mRecords) {
            if (rec.mDemotionReason == FitCadDemotionReason::None)
                continue;
            const auto f = facets_by_label.find(label);
            if (f == facets_by_label.end())
                continue;
            auto cut = fcn_resegment(rFc, *f->second, 2);
            if (!cut.empty())
                cuts.emplace_back(label, rec.mEscalationSteps, std::move(cut));
        }
        if (cuts.empty())
            break;
        const std::size_t from = work.size();
        for (auto& [label, steps, cut] : cuts) {
            log::debug("fit_cad: re-segmentation round {} cuts patch {} into {} piece(s)",
                       round + 1, label, cut.size());
            out.mRecords.erase(label);
            ++out.mNumResegmented;
            for (std::vector<std::size_t>& facets : cut) {
                const std::int64_t piece = next_label++;
                work.push_back({piece, facets, rOptions.mMaxFreeformSplits, steps + 1,
                                rOptions.mResegmentRounds});
                pieces.emplace_back(piece, std::move(facets));
            }
        }
        run(from);
    }
    if (pieces.empty())
        return out;

    // The mesh to write: the pieces under their labels, each a fit-table row.
    out.mSplit = true;
    out.mMesh = detail::clone_mesh(
        rSource, [](DataLocation, const std::string&, std::string&) { return true; });
    std::vector<std::int64_t> new_label = cell_label;
    for (const auto& p : pieces)
        for (std::size_t f : p.second)
            new_label[static_cast<std::size_t>(rFc.mCell[f])] = p.first;
    std::vector<NDArray> blocks;
    for (std::size_t b = 0; b < rSource.NumCellBlocks(); ++b) {
        const std::size_t nc = rSource.Cells(b).NumCells();
        NDArray a = NDArray::Uninit(DType::Int64, {nc});
        for (std::size_t c = 0; c < nc; ++c)
            a.As<std::int64_t>()[c] = new_label[static_cast<std::size_t>(bases[b]) + c];
        blocks.push_back(std::move(a));
    }
    out.mMesh.AddCellData("segment:patch", std::move(blocks));
    const NDArray& params = rSource.FieldData("fit:patch_params");
    const NDArray& tols = rSource.FieldData("fit:patch_tolerance");
    const std::size_t m = ids.Size(), mp = m + pieces.size();
    const std::size_t np = m > 0 ? params.Size() / m : 10;
    NDArray nid = NDArray::Uninit(DType::Int64, {mp});
    NDArray ntype = NDArray::Uninit(DType::Int64, {mp});
    NDArray nparam = NDArray::Uninit(DType::Float64, {mp, np});
    NDArray ntol = NDArray::Uninit(DType::Float64, {mp});
    for (std::size_t r = 0; r < m; ++r) {
        nid.As<std::int64_t>()[r] = detail::read_int(ids, r);
        ntype.As<std::int64_t>()[r] = detail::read_int(types, r);
        for (std::size_t k = 0; k < np; ++k)
            nparam.As<double>()[r * np + k] = detail::read_double(params, r * np + k);
        ntol.As<double>()[r] = detail::read_double(tols, r);
    }
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        nid.As<std::int64_t>()[m + i] = pieces[i].first;
        ntype.As<std::int64_t>()[m + i] = 0;
        for (std::size_t k = 0; k < np; ++k)
            nparam.As<double>()[(m + i) * np + k] = std::numeric_limits<double>::quiet_NaN();
        ntol.As<double>()[m + i] = Tol;
    }
    out.mMesh.AddFieldData("fit:patch_id", std::move(nid));
    out.mMesh.AddFieldData("fit:patch_type", std::move(ntype));
    out.mMesh.AddFieldData("fit:patch_params", std::move(nparam));
    out.mMesh.AddFieldData("fit:patch_tolerance", std::move(ntol));
    return out;
}

FitCadResult fcn_native_fit_cad(const Mesh& rMesh, const std::string& rStepPath,
                                const FitCadOptions& rOptions) {
    // A write entry point: the STEP header's provenance block holds the notes
    // this call raises (its own segmentation and fit included), never ones an
    // earlier, unrelated operation on this thread left behind.
    detail::provenance_begin_write();
    Mesh owned;
    FcnRefined refined;
    const Mesh& work = *fcn_prepare(rMesh, owned, rOptions, refined);
    const detail::ClosedFacets fc = detail::closed_facets(work, kFcnPrefix);
    if (fc.NumFacets() == 0)
        throw std::invalid_argument(std::string(kFcnPrefix) + "no surface facets");
    if (fc.mNumSkipped > 0)
        log::warn("fit_cad: {} cell(s) that are not surface facets (lines, vertices) are ignored",
                  fc.mNumSkipped);

    // The patches whose type is left out, or that a refinement found no
    // primitive for, are written as coplanar groups.
    Mesh typed;
    const Mesh* pSource = &work;
    const bool verdicts = refined.mRan && work.HasFieldData(kRefinePatchProxyName) &&
                          work.HasFieldData(kRefinePatchIdName) &&
                          work.HasFieldData("fit:patch_id");
    if (((rOptions.mTypes & kFitAllTypes) != kFitAllTypes || verdicts) &&
        work.HasFieldData("fit:patch_type")) {
        typed = detail::clone_mesh(
            work, [](DataLocation, const std::string&, std::string&) { return true; });
        const NDArray& types = work.FieldData("fit:patch_type");
        std::unordered_map<std::int64_t, std::int64_t> proxy_of;
        if (verdicts) {
            const NDArray& kinds = work.FieldData(kRefinePatchProxyName);
            const NDArray& labels = work.FieldData(kRefinePatchIdName);
            for (std::size_t i = 0; i < std::min(kinds.Size(), labels.Size()); ++i)
                proxy_of.emplace(detail::read_int(labels, i), detail::read_int(kinds, i));
        }
        NDArray kept = NDArray::Uninit(DType::Int64, {types.Size()});
        for (std::size_t r = 0; r < types.Size(); ++r) {
            const std::int64_t t = detail::read_int(types, r);
            bool drop = t >= 1 && t <= 5 && !(rOptions.mTypes & (1u << static_cast<unsigned>(t)));
            if (verdicts && r < work.FieldData("fit:patch_id").Size()) {
                const auto k = proxy_of.find(detail::read_int(work.FieldData("fit:patch_id"), r));
                drop = drop || (k != proxy_of.end() &&
                                (k->second == static_cast<std::int64_t>(RefineProxy::Facets) ||
                                 k->second == static_cast<std::int64_t>(RefineProxy::Freeform)));
            }
            kept.As<std::int64_t>()[r] = drop ? 0 : t;
        }
        typed.AddFieldData("fit:patch_type", std::move(kept));
        pSource = &typed;
    }
    // Outside the legacy mode, the patches the fit table leaves unexplained (its
    // type 0, or a refinement's freeform or facets verdict; not a type left
    // out of `mTypes`) are tried as adaptive freeform faces (v17.12.0).
    FcnFreeform freeform;
    if (rOptions.mFreeformMode != FitCadFreeformMode::Legacy && work.HasFieldData("fit:patch_id") &&
        work.HasFieldData("fit:patch_type")) {
        std::unordered_map<std::int64_t, std::int64_t> proxy_of;
        if (work.HasFieldData(kRefinePatchProxyName) && work.HasFieldData(kRefinePatchIdName)) {
            const NDArray& kinds = work.FieldData(kRefinePatchProxyName);
            const NDArray& labels = work.FieldData(kRefinePatchIdName);
            for (std::size_t i = 0; i < std::min(kinds.Size(), labels.Size()); ++i)
                proxy_of.emplace(detail::read_int(labels, i), detail::read_int(kinds, i));
        }
        const NDArray& ids = work.FieldData("fit:patch_id");
        const NDArray& types = work.FieldData("fit:patch_type");
        std::vector<std::int64_t> candidates;
        for (std::size_t r = 0; r < std::min(ids.Size(), types.Size()); ++r) {
            const std::int64_t label = detail::read_int(ids, r);
            const auto k = proxy_of.find(label);
            if (detail::read_int(types, r) == 0 ||
                (k != proxy_of.end() &&
                 (k->second == static_cast<std::int64_t>(RefineProxy::Facets) ||
                  k->second == static_cast<std::int64_t>(RefineProxy::Freeform))))
                candidates.push_back(label);
        }
        if (rOptions.mFreeform && !candidates.empty()) {
            const double wtol = fcn_tolerance(*pSource, fc, rOptions.mTolerance);
            freeform = fcn_freeform(*pSource, fc, candidates, rOptions, wtol);
            if (freeform.mSplit)
                pSource = &freeform.mMesh;
        }
    }
    StepWriteOptions opts;
    opts.mFaces = StepFaces::Fitted;
    opts.mTolerance = rOptions.mTolerance;
    detail::StepBrepResult brep;
    detail::write_step_report(rStepPath, *pSource, opts, &brep,
                              freeform.mFaces.empty() ? nullptr : &freeform.mFaces,
                              freeform.mEdgeCurves.empty() ? nullptr : &freeform.mEdgeCurves,
                              rOptions.mFreeformMode != FitCadFreeformMode::Legacy);

    FitCadResult result;
    result.mBackend = FitCadBackend::Native;
    // Outside the legacy mode, the solid as NURBS too (v17.12.0), for
    // `write_cad_json`; empty, with a warning, where a face cannot be one.
    if (rOptions.mFreeformMode != FitCadFreeformMode::Legacy) {
        result.mBrep = std::move(brep.mBrep);
        if (!brep.mBrepWhy.empty())
            log::warn("fit_cad: the native solid has no NURBS form ({}); the B-rep is empty",
                      brep.mBrepWhy);
    }
    result.mValid = brep.mValid;
    result.mTolerance = brep.mTolerance;
    result.mSkinVolume = brep.mVolume;
    result.mVolume = std::numeric_limits<double>::quiet_NaN();
    result.mVolumeRatio = std::numeric_limits<double>::quiet_NaN();
    result.mMaxDeviation = brep.mMaxDeviation;
    result.mNumFaces = brep.mNumFaces;
    result.mNumPlane = brep.mNumPlane;
    result.mNumSphere = brep.mNumSphere;
    result.mNumCylinder = brep.mNumCylinder;
    result.mNumCone = brep.mNumCone;
    result.mNumTorus = brep.mNumTorus;
    result.mNumBspline = brep.mNumBspline;
    result.mNumFreeformAdaptive = brep.mNumBspline;
    result.mNumNetworkEdges = freeform.mNumNetworkEdges;
    result.mNumExtraordinaryVertices = freeform.mNumExtraordinary;
    result.mMaxG1AngleDeg = freeform.mMaxG1AngleDeg;
    result.mNumSplitPatches = freeform.mNumSplit;
    result.mNumResegmented = freeform.mNumResegmented;
    result.mNumFacetedFaces = brep.mNumFacetedFaces;
    result.mNumPatches = brep.mNumPatches;
    result.mNumDemoted = brep.mNumDemoted;
    result.mNumExactEdges = brep.mNumLineEdges + brep.mNumCircleEdges;
    result.mNumApproxEdges = brep.mNumSplineEdges;
    result.mNumSolids = brep.mNumSolids;
    result.mClosed = brep.mValid;

    // The skin itself, each cell with the face it was written into.
    result.mMesh = detail::clone_mesh(
        work, [](DataLocation, const std::string&, std::string&) { return true; });
    const std::vector<std::int64_t> bases = detail::block_bases(work);
    const std::size_t ncells = static_cast<std::size_t>(detail::total_cells(bases));
    std::vector<std::int64_t> face(ncells, -1), type(ncells, 0), patch(ncells, -1);
    for (std::size_t f = fc.NumFacets(); f-- > 0;) {
        const std::size_t c = static_cast<std::size_t>(fc.mCell[f]);
        face[c] = brep.mFaceOfFacet[f];
        type[c] = brep.mTypeOfFacet[f];
        patch[c] = brep.mPatchOfFacet[f];
    }
    std::vector<NDArray> fb, tb, pb;
    for (std::size_t b = 0; b < work.NumCellBlocks(); ++b) {
        const std::size_t base = static_cast<std::size_t>(bases[b]);
        const std::size_t nc = work.Cells(b).NumCells();
        NDArray a = NDArray::Uninit(DType::Int64, {nc});
        NDArray t = NDArray::Uninit(DType::Int64, {nc});
        NDArray q = NDArray::Uninit(DType::Int64, {nc});
        for (std::size_t c = 0; c < nc; ++c) {
            a.As<std::int64_t>()[c] = face[base + c];
            t.As<std::int64_t>()[c] = type[base + c];
            q.As<std::int64_t>()[c] = patch[base + c];
        }
        fb.push_back(std::move(a));
        tb.push_back(std::move(t));
        pb.push_back(std::move(q));
    }
    result.mMesh.AddCellData(kFitCadFaceName, std::move(fb));
    result.mMesh.AddCellData(kFitCadTypeName, std::move(tb));
    result.mMesh.AddCellData(kFitCadPatchName, std::move(pb));

    // The freeform-mode report (v17.12.0): the facet area share over the skin
    // fanned into triangles, and a record per patch tried as a surface.
    result.mFreeformMode = rOptions.mFreeformMode;
    {
        std::vector<std::int64_t> tri, tri_type;
        for (std::size_t f = 0; f < fc.NumFacets(); ++f)
            for (std::int64_t k = fc.mStart[f] + 1; k + 1 < fc.mStart[f + 1]; ++k) {
                tri.insert(tri.end(), {fc.mNodes[static_cast<std::size_t>(fc.mStart[f])],
                                       fc.mNodes[static_cast<std::size_t>(k)],
                                       fc.mNodes[static_cast<std::size_t>(k + 1)]});
                tri_type.push_back(brep.mTypeOfFacet[f]);
            }
        result.mFacetAreaFraction = fcn_facet_area_fraction(fc.mXyz, tri, tri_type);
    }
    {
        std::unordered_map<std::int64_t, std::pair<std::int64_t, std::int64_t>> first_face;
        for (std::size_t f = 0; f < fc.NumFacets(); ++f)
            if (brep.mTypeOfFacet[f] != 0)
                first_face.emplace(brep.mPatchOfFacet[f],
                                   std::make_pair(brep.mFaceOfFacet[f], brep.mTypeOfFacet[f]));
        std::unordered_map<std::int64_t, std::size_t> num_facets;
        for (std::size_t f = 0; f < fc.NumFacets(); ++f)
            ++num_facets[brep.mPatchOfFacet[f]];
        const NDArray* pTypes = pSource->HasFieldData("fit:patch_type")
                                    ? &pSource->FieldData("fit:patch_type")
                                    : nullptr;
        double freeform_sq = 0.0;
        std::size_t freeform_n = 0;
        for (std::size_t r = 0; r < brep.mRowLabel.size(); ++r) {
            const std::int64_t lbl = brep.mRowLabel[r];
            const std::int64_t t = pTypes && r < pTypes->Size() ? detail::read_int(*pTypes, r) : 0;
            if (num_facets.find(lbl) == num_facets.end())
                continue;
            const auto ff = freeform.mRecords.find(lbl);
            if (ff != freeform.mRecords.end()) {
                // A patch the adaptive pass tried: its fit's record, and the
                // writer's verdict on it.
                FitCadFaceRecord rec = ff->second;
                const auto it = first_face.find(lbl);
                if (rec.mDemotionReason == FitCadDemotionReason::None && it != first_face.end() &&
                    brep.mRowReason[r] == 0) {
                    rec.mFace = it->second.first;
                    rec.mKind = it->second.second;
                } else if (rec.mDemotionReason == FitCadDemotionReason::None) {
                    rec.mDemotionReason =
                        brep.mRowReason[r] != 0
                            ? static_cast<FitCadDemotionReason>(brep.mRowReason[r])
                            : FitCadDemotionReason::Other;
                }
                result.mNumEscalations += rec.mEscalationSteps;
                if (rec.mKind == kFitCadBsplineType && freeform.mNetworkLabels.count(lbl)) {
                    ++result.mNumNetworkPatches;
                    --result.mNumFreeformAdaptive;
                }
                if (rec.mKind == kFitCadBsplineType) {
                    result.mMaxFreeformDistance =
                        std::max(result.mMaxFreeformDistance, rec.mMaxDist);
                    freeform_sq += rec.mRmsDist * rec.mRmsDist;
                    ++freeform_n;
                }
                result.mFaces.push_back(rec);
                continue;
            }
            if (t < 1 || t > 5)
                continue;
            FitCadFaceRecord rec;
            rec.mPatch = lbl;
            const auto it = first_face.find(lbl);
            if (it != first_face.end() && brep.mRowReason[r] == 0) {
                rec.mFace = it->second.first;
                rec.mKind = it->second.second;
            } else {
                rec.mDemotionReason = static_cast<FitCadDemotionReason>(
                    brep.mRowReason[r] != 0 ? brep.mRowReason[r] : 10);
            }
            result.mFaces.push_back(rec);
        }
        if (freeform_n > 0)
            result.mRmsFreeformDistance = std::sqrt(freeform_sq / static_cast<double>(freeform_n));
    }
    if (rOptions.mFreeformMode != FitCadFreeformMode::Legacy)
        result.mMesh.AddFieldData(kFitCadFaceTableName, fcn_face_table(result.mFaces));
    if (result.mNumDemoted > 0)
        log::info("fit_cad: {} fitted patch(es) written as coplanar facets", result.mNumDemoted);
    fcn_report_refined(refined, work, result);
    return result;
}

}  // namespace

#ifdef MESHIOPLUSPLUS_HAS_OCCT

namespace {

constexpr const char* kFcPrefix = "fit_cad: ";
constexpr std::size_t kFcParams = 10;
constexpr double kFcPi = 3.14159265358979323846;
// Freeform tiles: a patch no single B-spline follows is halved up to this many
// times (at most 2^depth tiles), and a tile with fewer triangles stays facets.
constexpr int kFcTileMaxDepth = 9;
constexpr std::size_t kFcTileMinTriangles = 24;

template <class T>
using FcHandle = opencascade::handle<T>;

// The OpenCASCADE messenger prints transfer statistics to stdout; a library
// must not, so its printers are detached for the duration of the call.
class FcQuietMessenger {
public:
    FcQuietMessenger() : mMessenger(Message::DefaultMessenger()) {
        if (!mMessenger.IsNull()) {
            mSaved = mMessenger->Printers();
            mMessenger->ChangePrinters().Clear();
        }
    }
    FcQuietMessenger(const FcQuietMessenger&) = delete;
    FcQuietMessenger& operator=(const FcQuietMessenger&) = delete;
    ~FcQuietMessenger() {
        if (!mMessenger.IsNull())
            mMessenger->ChangePrinters() = mSaved;
    }

private:
    FcHandle<Message_Messenger> mMessenger;
    Message_SequenceOfPrinters mSaved;
};

gp_Pnt fc_point(const std::vector<double>& rXyz, std::int64_t V) {
    const std::size_t i = static_cast<std::size_t>(V) * 3;
    return gp_Pnt(rXyz[i], rXyz[i + 1], rXyz[i + 2]);
}

std::uint64_t fc_edge_key(std::int64_t U, std::int64_t V) {
    const std::uint64_t lo = static_cast<std::uint64_t>(std::min(U, V));
    const std::uint64_t hi = static_cast<std::uint64_t>(std::max(U, V));
    return (hi << 32) ^ lo ^ (hi >> 32);
}

struct FcEdgeKeyHash {
    std::size_t operator()(const std::pair<std::int64_t, std::int64_t>& rKey) const {
        return static_cast<std::size_t>(fc_edge_key(rKey.first, rKey.second) *
                                        0x9E3779B97F4A7C15ull);
    }
};

// One fitted patch row of the `fit:patch_*` table.
struct FcPatch {
    std::int64_t mLabel = -1;
    std::int32_t mType = 0;
    std::array<double, kFcParams> mParams{};
    double mTolerance = 0.0;
    bool mAccepted = false;  // a type in the mask with finite parameters
};

// The nearest point of @p rS to @p rP; false when the projection fails.
// One projector per surface for the length of a fit_cad call: building one
// samples the surface, which on a B-spline costs far more than a projection,
// and the edge builder projects the same few surfaces tens of thousands of
// times. Keyed by handle, so a surface lives as long as its entry; cleared by
// `FcProjectorScope` when fit_cad returns.
struct FcProjector {
    FcHandle<Geom_Surface> mSurface;
    std::unique_ptr<GeomAPI_ProjectPointOnSurf> mProj;
};
thread_local std::vector<FcProjector> tFcProjectors;

struct FcProjectorScope {
    FcProjectorScope() { tFcProjectors.clear(); }
    ~FcProjectorScope() { tFcProjectors.clear(); }
};

GeomAPI_ProjectPointOnSurf& fc_projector(const FcHandle<Geom_Surface>& rS) {
    for (FcProjector& p : tFcProjectors)
        if (p.mSurface == rS)
            return *p.mProj;
    double u0 = 0.0, u1 = 0.0, v0 = 0.0, v1 = 0.0;
    rS->Bounds(u0, u1, v0, v1);
    auto proj = std::make_unique<GeomAPI_ProjectPointOnSurf>();
    proj->Init(rS, u0, u1, v0, v1);
    tFcProjectors.push_back({rS, std::move(proj)});
    return *tFcProjectors.back().mProj;
}

bool fc_project(const FcHandle<Geom_Surface>& rS, const gp_Pnt& rP, gp_Pnt& rOut,
                double* pU = nullptr, double* pV = nullptr) {
    GeomAPI_ProjectPointOnSurf& proj = fc_projector(rS);
    proj.Perform(rP);
    if (!proj.IsDone() || proj.NbPoints() == 0)
        return false;
    rOut = proj.NearestPoint();
    if (pU && pV)
        proj.LowerDistanceParameters(*pU, *pV);
    return true;
}

// Move @p rP onto every surface of @p rSurfaces at once, by averaged
// alternating projection; returns the largest remaining distance to any of
// them, or +inf when a projection fails. With @p Close, it also stops once
// the point lies within that distance of every surface: between tangent
// surfaces the averaging creeps, and a fillet's or a thread's edge needs no
// more than a small fraction of the tolerance.
double fc_onto_all(const std::vector<FcHandle<Geom_Surface> >& rSurfaces, gp_Pnt& rP,
                   double Scale, double Close = 0.0) {
    const double inf = std::numeric_limits<double>::infinity();
    if (rSurfaces.empty())
        return 0.0;
    for (int it = 0; it < 200; ++it) {
        gp_XYZ sum(0.0, 0.0, 0.0);
        double max_dist = 0.0;
        for (const auto& s : rSurfaces) {
            gp_Pnt q;
            if (!fc_project(s, rP, q))
                return inf;
            sum += q.XYZ();
            max_dist = std::max(max_dist, q.Distance(rP));
        }
        if (Close > 0.0 && max_dist <= Close)
            break;
        const gp_Pnt next(sum / static_cast<double>(rSurfaces.size()));
        const double moved = next.Distance(rP);
        rP = next;
        if (moved <= 1e-13 * Scale)
            break;
    }
    double worst = 0.0;
    for (const auto& s : rSurfaces) {
        gp_Pnt q;
        if (!fc_project(s, rP, q))
            return inf;
        worst = std::max(worst, q.Distance(rP));
    }
    return worst;
}

// The fitted surface of a patch, or a null handle when the parameters do not
// make one (a cone whose patch lies behind its apex, a degenerate torus).
FcHandle<Geom_Surface> fc_surface(const FcPatch& rPatch, const std::vector<gp_Pnt>& rPts) {
    const double* p = rPatch.mParams.data();
    const gp_Pnt o(p[0], p[1], p[2]);
    const gp_Vec a(p[3], p[4], p[5]);
    if (rPatch.mType != 2 && a.Magnitude() < 1e-12)
        return {};
    switch (rPatch.mType) {
        case 1:
            return new Geom_Plane(gp_Pln(o, gp_Dir(a)));
        case 2: {
            // The pole axis is the one, of 64 spread over the sphere, whose
            // poles stay furthest from the patch: a patch smaller than a
            // hemisphere then has loops that close in parameter space and need
            // no seam, and no pole sits in a hole the face does not cover.
            if (!(p[6] > 0.0))
                return {};
            gp_Dir best(0.0, 0.0, 1.0);
            double best_score = -1.0;
            constexpr int kAxes = 64;
            for (int k = 0; k < kAxes; ++k) {
                const double z = 1.0 - (k + 0.5) / kAxes;  // a hemisphere suffices: +-axis
                const double rr = std::sqrt(std::max(0.0, 1.0 - z * z));
                const double phi = k * 2.399963229728653;  // the golden angle
                const gp_Vec ax(rr * std::cos(phi), rr * std::sin(phi), z);
                double worst = 0.0;
                for (const gp_Pnt& q : rPts) {
                    const gp_Vec d(o, q);
                    const double len = d.Magnitude();
                    if (len > 0.0)
                        worst = std::max(worst, std::fabs(d.Dot(ax)) / len);
                }
                if (1.0 - worst > best_score) {
                    best_score = 1.0 - worst;
                    best = gp_Dir(ax);
                }
            }
            return new Geom_SphericalSurface(gp_Sphere(gp_Ax3(o, best), p[6]));
        }
        case 3:
            if (!(p[6] > 0.0))
                return {};
            return new Geom_CylindricalSurface(gp_Cylinder(gp_Ax3(o, gp_Dir(a)), p[6]));
        case 4: {
            // OpenCASCADE's cone is a reference circle and a half-angle: take the
            // circle at the patch's mean height above the apex.
            const double ang = p[8] * kFcPi / 180.0;
            if (!(ang > 0.0 && ang < kFcPi / 2) || rPts.empty())
                return {};
            const gp_Dir ax(a);
            double h = 0.0;
            for (const gp_Pnt& q : rPts)
                h += gp_Vec(o, q).Dot(gp_Vec(ax));
            h /= static_cast<double>(rPts.size());
            if (!(h > 0.0))
                return {};
            const gp_Pnt at = o.Translated(gp_Vec(ax) * h);
            return new Geom_ConicalSurface(gp_Cone(gp_Ax3(at, ax), ang, h * std::tan(ang)));
        }
        case 5:
            if (!(p[6] > 0.0 && p[7] > 0.0 && p[7] < p[6]))
                return {};
            return new Geom_ToroidalSurface(gp_Torus(gp_Ax3(o, gp_Dir(a)), p[6], p[7]));
        default:
            return {};
    }
}

// The surface's own unit normal at the projection of @p rP, or a zero vector.
gp_Vec fc_normal_at(const FcHandle<Geom_Surface>& rS, const gp_Pnt& rP) {
    gp_Pnt q;
    double u = 0.0, v = 0.0;
    if (!fc_project(rS, rP, q, &u, &v))
        return gp_Vec(0.0, 0.0, 0.0);
    GeomLProp_SLProps props(rS, u, v, 1, Precision::Confusion());
    if (!props.IsNormalDefined())
        return gp_Vec(0.0, 0.0, 0.0);
    return gp_Vec(props.Normal());
}

// Give every edge of @p rWires a parameter-space curve on @p rFace that is
// continuous around each loop and across loops, in one period window. Left to
// itself, ShapeFix projects each edge into the surface's own [0, 2 pi) window,
// so a loop running along a seam (a cylinder's boundary through its u = 0
// line, a fillet's tangent circle on a torus's v = 0 equator) comes out in
// pieces a period apart, and the missing-seam repair then drops a whole wire.
void fc_set_pcurves(const TopoDS_Face& rFace, const std::vector<TopoDS_Wire>& rWires,
                    const FcHandle<Geom_Surface>& rS, double Tol) {
    const double uper = rS->IsUPeriodic() ? rS->UPeriod() : 0.0;
    const double vper = rS->IsVPeriodic() ? rS->VPeriod() : 0.0;
    auto nearest_copy = [](double X, double Ref, double Per) {
        return Per > 0.0 ? X + Per * std::round((Ref - X) / Per) : X;
    };
    constexpr int kSamples = 33;
    BRep_Builder builder;
    bool have_prev = false;
    double pu = 0.0, pv = 0.0;
    for (const TopoDS_Wire& w : rWires)
        for (TopoDS_Iterator it(w); it.More(); it.Next()) {
            const TopoDS_Edge e = TopoDS::Edge(it.Value());
            double a = 0.0, b = 0.0;
            const FcHandle<Geom_Curve> c = BRep_Tool::Curve(e, a, b);
            if (c.IsNull())
                throw Standard_Failure("an edge has no 3-D curve");
            const bool rev = e.Orientation() == TopAbs_REVERSED;
            std::vector<gp_Pnt2d> uv(kSamples);
            for (int k = 0; k < kSamples; ++k) {
                const int j = rev ? kSamples - 1 - k : k;  // walk in loop order
                const double t = a + (b - a) * j / (kSamples - 1);
                GeomAPI_ProjectPointOnSurf proj(c->Value(t), rS);
                if (!proj.IsDone() || proj.NbPoints() == 0)
                    throw Standard_Failure("an edge does not project onto its face");
                double u = 0.0, v = 0.0;
                proj.LowerDistanceParameters(u, v);
                if (have_prev) {
                    u = nearest_copy(u, pu, uper);
                    v = nearest_copy(v, pv, vper);
                }
                pu = u;
                pv = v;
                have_prev = true;
                uv[static_cast<std::size_t>(j)] = gp_Pnt2d(u, v);
            }
            FcHandle<TColgp_HArray1OfPnt2d> pts = new TColgp_HArray1OfPnt2d(1, kSamples);
            FcHandle<TColStd_HArray1OfReal> par = new TColStd_HArray1OfReal(1, kSamples);
            for (int k = 0; k < kSamples; ++k) {
                pts->SetValue(k + 1, uv[static_cast<std::size_t>(k)]);
                par->SetValue(k + 1, a + (b - a) * k / (kSamples - 1));
            }
            Geom2dAPI_Interpolate interp(pts, par, Standard_False, 1e-12);
            interp.Perform();
            if (!interp.IsDone())
                throw Standard_Failure("no parameter-space curve through an edge");
            builder.UpdateEdge(e, interp.Curve(), rFace,
                               std::max(BRep_Tool::Tolerance(e), 1e-3 * Tol));
        }
}

// The B-rep check's complaints about @p rShape and its sub-shapes, for the
// debug log; empty when it is valid.
std::string fc_check_report(const TopoDS_Shape& rShape) {
    BRepCheck_Analyzer check(rShape);
    if (check.IsValid())
        return {};
    std::string out;
    for (TopAbs_ShapeEnum kind : {TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_VERTEX})
        for (TopExp_Explorer ex(rShape, kind); ex.More(); ex.Next()) {
            const FcHandle<BRepCheck_Result> res = check.Result(ex.Current());
            if (res.IsNull())
                continue;
            for (BRepCheck_ListIteratorOfListOfStatus it(res->Status()); it.More(); it.Next())
                if (it.Value() != BRepCheck_NoError)
                    out += " " + std::to_string(static_cast<int>(kind)) + ":" +
                           std::to_string(static_cast<int>(it.Value()));
        }
    return out;
}

// The skin, its facets and patch labels: everything a build round reads.
struct FcSkin {
    detail::ClosedFacets mFacets;
    std::vector<std::int64_t> mRow;  // per facet: its patch row, or -1
    std::vector<FcPatch> mPatches;
    std::vector<std::size_t> mPiece;  // per facet: its edge-connected piece
    std::size_t mNumPieces = 0;
    double mScale = 1.0;  // bounding-box diagonal
    // Per undirected edge, the facets on either side: the one walking it lo->hi first.
    std::unordered_map<std::pair<std::int64_t, std::int64_t>, std::array<std::int64_t, 2>,
                       FcEdgeKeyHash>
        mEdgeFacets;
};

// What `fit_cad` learned about one patch row on its way, for its
// `FitCadFaceRecord` (v17.12.0).
struct FcRowInfo {
    std::uint8_t mReason = 0;  // a FitCadDemotionReason, 0 while not demoted
    bool mTried = false;       // tried as a freeform surface
    FitCadParamMethod mParam = FitCadParamMethod::None;
    std::int64_t mNumControl = 0;
    std::int64_t mIterations = 0;
    std::int64_t mRefinements = 0;
    std::int64_t mEscalations = 0;
    double mMaxDist = std::numeric_limits<double>::quiet_NaN();
    double mRmsDist = std::numeric_limits<double>::quiet_NaN();
    double mReverseMaxDist = std::numeric_limits<double>::quiet_NaN();
    double mMaxG1AngleDeg = std::numeric_limits<double>::quiet_NaN();
    // Outside the legacy mode: the distance the face may lie from its patch
    // (the tolerance, or the patch's noise floor when larger); 0 otherwise.
    double mAllowance = 0.0;
};

// `rInfo[R]`, grown to hold it (split pieces append rows).
FcRowInfo& fc_row(std::vector<FcRowInfo>& rInfo, std::size_t R) {
    if (R >= rInfo.size())
        rInfo.resize(R + 1);
    return rInfo[R];
}

// What one build round produced.
struct FcBuild {
    TopoDS_Shape mShape;
    std::vector<std::int64_t> mFaceGroup;  // per final face (explorer order): its group
    std::vector<TopoDS_Face> mFaces;       // final faces, explorer order
    std::vector<std::uint8_t> mBadRow;     // per patch row: demote it
    // Facets of coplanar groups whose faces failed: one face each next time.
    std::vector<std::size_t> mUnmerge;
    std::int64_t mNumExact = 0;
    std::int64_t mNumApprox = 0;
    std::int64_t mNumSnapped = 0;   // vertices the least-squares junction solve placed
    std::int64_t mNumExtended = 0;  // exact edges found on extended surfaces
    bool mValid = false;
};

// A chain of patch-boundary mesh edges between two vertices that become B-rep
// vertices, with the two groups on either side.
struct FcChain {
    std::vector<std::int64_t> mVerts;
    std::int64_t mGroupFwd = -1;  // the group whose facets walk mVerts forwards
    std::int64_t mGroupBwd = -1;
    TopoDS_Edge mEdge;
    bool mFromStart = true;  // the edge's first vertex is mVerts.front()
    bool mOk = false;
};

// Build the solid for one choice of analytic patches. Every failure is
// recorded against the patch rows it touches in `mBadRow`, never thrown.
FcBuild fc_build(const FcSkin& rSkin, const std::vector<FcHandle<Geom_Surface>>& rSurface,
                 const std::vector<std::uint8_t>& rAnalytic, double Tol, bool ExactEdges,
                 bool MergeFaces, const std::vector<std::uint8_t>& rNoMerge, double SewTol = 0.0,
                 bool SnapJunctions = false,
                 const std::vector<FcHandle<Geom_Surface>>* pExtended = nullptr,
                 const std::vector<detail::StepEdgeCurve>* pNetworkCurves = nullptr,
                 const std::vector<std::uint8_t>* pNetworkRow = nullptr,
                 bool FreeformApprox = false) {
    const detail::ClosedFacets& fc = rSkin.mFacets;
    const std::vector<double>& xyz = fc.mXyz;
    const std::size_t nf = fc.NumFacets();
    const std::int64_t np = static_cast<std::int64_t>(rSkin.mPatches.size());
    const std::size_t npts = xyz.size() / 3;
    FcBuild out;
    out.mBadRow.assign(rSkin.mPatches.size(), 0);

    // The facets outside every analytic face, gathered into edge-connected
    // coplanar groups -- largest facet first, every corner within 1e-6 of the
    // skin's diagonal of that facet's plane -- each one planar face named
    // `np + its seed facet`, so a demoted or unexplained flat patch is one
    // polygon (holes and all), not a face per triangle.
    std::vector<std::int64_t> flat_group(nf, -1);
    {
        auto analytic_facet = [&](std::size_t F) {
            const std::int64_t r = rSkin.mRow[F];
            return r >= 0 && rAnalytic[static_cast<std::size_t>(r)];
        };
        std::vector<double> area(nf, 0.0);
        for (std::size_t f = 0; f < nf; ++f) {
            const std::int64_t b = fc.mStart[f], n = fc.mStart[f + 1] - b;
            gp_XYZ twice(0.0, 0.0, 0.0);
            for (std::int64_t k = 0; k < n; ++k)
                twice +=
                    fc_point(xyz, fc.mNodes[static_cast<std::size_t>(b + k)])
                        .XYZ()
                        .Crossed(fc_point(xyz, fc.mNodes[static_cast<std::size_t>(b + (k + 1) % n)])
                                     .XYZ());
            area[f] = 0.5 * twice.Modulus();
        }
        std::vector<std::size_t> order(nf);
        std::iota(order.begin(), order.end(), std::size_t{0});
        std::stable_sort(order.begin(), order.end(),
                         [&](std::size_t A, std::size_t B) { return area[A] > area[B]; });
        const double eps = 1e-6 * rSkin.mScale;
        std::vector<std::size_t> queue;
        for (std::size_t seed : order) {
            if (flat_group[seed] != -1 || analytic_facet(seed))
                continue;
            const std::int64_t id = np + static_cast<std::int64_t>(seed);
            if (rNoMerge[seed]) {
                flat_group[seed] = id;
                continue;
            }
            const gp_XYZ p0 =
                fc_point(xyz, fc.mNodes[static_cast<std::size_t>(fc.mStart[seed])]).XYZ();
            const gp_XYZ n0(fc.mNormals[seed * 3], fc.mNormals[seed * 3 + 1],
                            fc.mNormals[seed * 3 + 2]);
            flat_group[seed] = id;
            queue.assign(1, seed);
            for (std::size_t q = 0; q < queue.size(); ++q) {
                const std::size_t f = queue[q];
                const std::int64_t b = fc.mStart[f], n = fc.mStart[f + 1] - b;
                for (std::int64_t k = 0; k < n; ++k) {
                    const std::int64_t u = fc.mNodes[static_cast<std::size_t>(b + k)];
                    const std::int64_t v = fc.mNodes[static_cast<std::size_t>(b + (k + 1) % n)];
                    auto it = rSkin.mEdgeFacets.find({std::min(u, v), std::max(u, v)});
                    if (it == rSkin.mEdgeFacets.end())
                        continue;
                    const std::int64_t other = it->second[0] == static_cast<std::int64_t>(f)
                                                   ? it->second[1]
                                                   : it->second[0];
                    if (other < 0)
                        continue;
                    const std::size_t o = static_cast<std::size_t>(other);
                    if (flat_group[o] != -1 || analytic_facet(o) || rNoMerge[o])
                        continue;
                    if (gp_XYZ(fc.mNormals[o * 3], fc.mNormals[o * 3 + 1], fc.mNormals[o * 3 + 2])
                            .Dot(n0) <= 0.0)
                        continue;
                    bool flat = true;
                    for (std::int64_t j = fc.mStart[o]; j < fc.mStart[o + 1] && flat; ++j)
                        flat =
                            std::fabs(
                                (fc_point(xyz, fc.mNodes[static_cast<std::size_t>(j)]).XYZ() - p0)
                                    .Dot(n0)) <= eps;
                    if (!flat)
                        continue;
                    flat_group[o] = id;
                    queue.push_back(o);
                }
            }
        }
    }

    auto group_of = [&](std::size_t F) -> std::int64_t {
        const std::int64_t r = rSkin.mRow[F];
        return r >= 0 && rAnalytic[static_cast<std::size_t>(r)] ? r : flat_group[F];
    };
    auto is_analytic = [&](std::int64_t G) { return G >= 0 && G < np; };
    // `mBadRow` holds the first `FitCadDemotionReason` a row is blamed for.
    auto mark_bad = [&](std::int64_t G, FitCadDemotionReason Why, const char* pWhy) {
        if (!is_analytic(G))
            return;
        if (!out.mBadRow[static_cast<std::size_t>(G)]) {
            log::debug("fit_cad: demoting patch {}: {}",
                       rSkin.mPatches[static_cast<std::size_t>(G)].mLabel, pWhy);
            out.mBadRow[static_cast<std::size_t>(G)] = static_cast<std::uint8_t>(Why);
        }
    };

    // The groups around every point, and whether a facet group is among them.
    std::vector<std::pair<std::int64_t, std::int64_t>> vg;
    vg.reserve(fc.mNodes.size());
    for (std::size_t f = 0; f < nf; ++f)
        for (std::int64_t k = fc.mStart[f]; k < fc.mStart[f + 1]; ++k)
            vg.emplace_back(fc.mNodes[static_cast<std::size_t>(k)], group_of(f));
    std::sort(vg.begin(), vg.end());
    vg.erase(std::unique(vg.begin(), vg.end()), vg.end());
    std::vector<std::int64_t> vg_start(npts + 1, 0);
    for (const auto& p : vg)
        ++vg_start[static_cast<std::size_t>(p.first) + 1];
    for (std::size_t v = 0; v < npts; ++v)
        vg_start[v + 1] += vg_start[v];
    std::vector<std::uint8_t> touches_facet(npts, 0);
    for (const auto& p : vg)
        if (!is_analytic(p.second))
            touches_facet[static_cast<std::size_t>(p.first)] = 1;

    // Boundary edges between groups, deterministic: ascending (lo, hi).
    std::vector<std::pair<std::int64_t, std::int64_t>> bedges;
    for (const auto& kv : rSkin.mEdgeFacets) {
        const std::int64_t ga = group_of(static_cast<std::size_t>(kv.second[0]));
        const std::int64_t gb = group_of(static_cast<std::size_t>(kv.second[1]));
        if (ga != gb)
            bedges.push_back(kv.first);
    }
    std::sort(bedges.begin(), bedges.end());
    std::vector<std::vector<std::int64_t>> nbr(npts);
    for (const auto& e : bedges) {
        nbr[static_cast<std::size_t>(e.first)].push_back(e.second);
        nbr[static_cast<std::size_t>(e.second)].push_back(e.first);
    }
    for (auto& n : nbr)
        std::sort(n.begin(), n.end());
    auto junction = [&](std::int64_t V) {
        const std::size_t v = static_cast<std::size_t>(V);
        return touches_facet[v] || vg_start[v + 1] - vg_start[v] >= 3 || nbr[v].size() != 2;
    };
    std::unordered_map<std::pair<std::int64_t, std::int64_t>, std::uint8_t, FcEdgeKeyHash> seen;
    seen.reserve(bedges.size());
    for (const auto& e : bedges)
        seen.emplace(e, 0);
    auto take = [&](std::int64_t U, std::int64_t V) -> bool {
        auto it = seen.find({std::min(U, V), std::max(U, V)});
        if (it == seen.end() || it->second)
            return false;
        it->second = 1;
        return true;
    };

    // Chains: from every junction along unvisited boundary edges, then the
    // closed loops no junction lies on. A chain that closes is split in two.
    std::vector<FcChain> chains;
    auto add_chain = [&](std::vector<std::int64_t> Verts) {
        if (Verts.front() == Verts.back()) {
            const std::size_t mid = Verts.size() / 2;
            FcChain a, b;
            a.mVerts.assign(Verts.begin(), Verts.begin() + static_cast<std::ptrdiff_t>(mid) + 1);
            b.mVerts.assign(Verts.begin() + static_cast<std::ptrdiff_t>(mid), Verts.end());
            chains.push_back(std::move(a));
            chains.push_back(std::move(b));
        } else {
            FcChain c;
            c.mVerts = std::move(Verts);
            chains.push_back(std::move(c));
        }
    };
    auto walk = [&](std::int64_t Start, std::int64_t First, bool StopAtStart) {
        std::vector<std::int64_t> verts{Start, First};
        std::int64_t prev = Start, cur = First;
        while (!(StopAtStart ? cur == Start : junction(cur))) {
            std::int64_t next = -1;
            for (std::int64_t n : nbr[static_cast<std::size_t>(cur)])
                if (n != prev && take(cur, n)) {
                    next = n;
                    break;
                }
            if (next < 0) {
                // Back to the start of a loop through its last edge.
                if (StopAtStart && take(cur, Start)) {
                    verts.push_back(Start);
                }
                break;
            }
            verts.push_back(next);
            prev = cur;
            cur = next;
        }
        return verts;
    };
    for (std::int64_t v = 0; v < static_cast<std::int64_t>(npts); ++v) {
        if (nbr[static_cast<std::size_t>(v)].empty() || !junction(v))
            continue;
        for (std::int64_t n : nbr[static_cast<std::size_t>(v)])
            if (take(v, n))
                add_chain(walk(v, n, false));
    }
    for (const auto& e : bedges) {
        if (!take(e.first, e.second))
            continue;
        add_chain(walk(e.first, e.second, true));
    }

    // The groups either side of each chain, from the facets of its first edge.
    for (FcChain& c : chains) {
        const std::int64_t u = c.mVerts[0], w = c.mVerts[1];
        const auto& two = rSkin.mEdgeFacets.at({std::min(u, w), std::max(u, w)});
        const std::size_t fwd = static_cast<std::size_t>(u < w ? two[0] : two[1]);
        const std::size_t bwd = static_cast<std::size_t>(u < w ? two[1] : two[0]);
        c.mGroupFwd = group_of(fwd);
        c.mGroupBwd = group_of(bwd);
    }

    // B-rep vertices at the chain ends. A point among analytic faces only is
    // moved onto their surfaces; a point beside a facet keeps the mesh position.
    BRep_Builder builder;
    std::unordered_map<std::int64_t, TopoDS_Vertex> vertex;
    auto groups_at = [&](std::int64_t V) {
        std::vector<std::int64_t> g;
        const std::size_t v = static_cast<std::size_t>(V);
        for (std::int64_t k = vg_start[v]; k < vg_start[v + 1]; ++k)
            g.push_back(vg[static_cast<std::size_t>(k)].second);
        return g;
    };
    for (const FcChain& c : chains)
        for (std::int64_t v : {c.mVerts.front(), c.mVerts.back()}) {
            if (vertex.count(v))
                continue;
            gp_Pnt p = fc_point(xyz, v);
            double vtol = Precision::Confusion();
            if (!touches_facet[static_cast<std::size_t>(v)]) {
                std::vector<FcHandle<Geom_Surface>> surfs;
                for (std::int64_t g : groups_at(v))
                    surfs.push_back(rSurface[static_cast<std::size_t>(g)]);
                if (SnapJunctions && surfs.size() >= 2) {
                    // The least-squares junction point first; the averaged
                    // projection below when it misses the same acceptance.
                    std::vector<detail::ClosestPoint> closest;
                    std::vector<detail::NormalAt> normal;
                    for (const auto& srf : surfs) {
                        closest.emplace_back(
                            [srf](const detail::PatchVec& rX, detail::PatchVec& rFoot) {
                                gp_Pnt f;
                                if (!fc_project(srf, gp_Pnt(rX[0], rX[1], rX[2]), f))
                                    return false;
                                rFoot = {f.X(), f.Y(), f.Z()};
                                return true;
                            });
                        normal.emplace_back(
                            [srf](const detail::PatchVec& rX, detail::PatchVec& rN) {
                                const gp_Vec n = fc_normal_at(srf, gp_Pnt(rX[0], rX[1], rX[2]));
                                const double len = n.Magnitude();
                                if (!(len > 0.0))
                                    return false;
                                rN = {n.X() / len, n.Y() / len, n.Z() / len};
                                return true;
                            });
                    }
                    const detail::JunctionPoint jp = detail::junction_closest(
                        closest, normal, {p.X(), p.Y(), p.Z()}, rSkin.mScale);
                    const gp_Pnt q(jp.mPoint[0], jp.mPoint[1], jp.mPoint[2]);
                    if (jp.mOk && jp.mMaxDistance <= Tol && q.Distance(p) <= Tol) {
                        builder.MakeVertex(vertex[v], q, std::max(vtol, 1.5 * jp.mMaxDistance));
                        ++out.mNumSnapped;
                        continue;
                    }
                }
                gp_Pnt q = p;
                const double d = fc_onto_all(surfs, q, rSkin.mScale);
                if (d <= Tol && q.Distance(p) <= Tol) {
                    p = q;
                    vtol = std::max(vtol, 1.5 * d);
                } else {
                    double worst = 0.0;
                    for (const auto& s : surfs) {
                        gp_Pnt r;
                        if (fc_project(s, p, r))
                            worst = std::max(worst, r.Distance(p));
                    }
                    vtol = std::max(vtol, 1.5 * worst);
                }
            }
            TopoDS_Vertex tv;
            builder.MakeVertex(tv, p, vtol);
            vertex.emplace(v, tv);
        }

    // A vertex sits on a curve when its tolerance covers the gap.
    auto fit_vertex = [&](const TopoDS_Vertex& rV, const gp_Pnt& rOn) {
        const double d = BRep_Tool::Pnt(rV).Distance(rOn);
        if (d > BRep_Tool::Tolerance(rV))
            builder.UpdateVertex(rV, 1.01 * d);
    };

    // The patch networks' shared curves, by the mesh edges they run along
    // (v17.12.0): the edge between two network patches is a piece of one.
    std::unordered_map<std::uint64_t, std::pair<std::size_t, std::size_t>> curve_at;
    if (pNetworkCurves != nullptr)
        for (std::size_t k = 0; k < pNetworkCurves->size(); ++k) {
            const detail::StepEdgeCurve& ec = (*pNetworkCurves)[k];
            for (std::size_t i = 0; i + 1 < ec.mVerts.size(); ++i)
                curve_at.emplace(fc_edge_key(ec.mVerts[i], ec.mVerts[i + 1]), std::make_pair(k, i));
        }
    const auto network_group = [&](std::int64_t G) {
        return pNetworkRow != nullptr && G >= 0 && G < np &&
               static_cast<std::size_t>(G) < pNetworkRow->size() &&
               (*pNetworkRow)[static_cast<std::size_t>(G)] != 0;
    };

    // Edges.
    for (FcChain& c : chains) {
        const TopoDS_Vertex& va = vertex.at(c.mVerts.front());
        const TopoDS_Vertex& vb = vertex.at(c.mVerts.back());
        try {
            // Between two network patches along one shared curve: its exact
            // piece, which both patches' boundaries are.
            if (!curve_at.empty() && c.mVerts.front() != c.mVerts.back() &&
                network_group(c.mGroupFwd) && network_group(c.mGroupBwd)) {
                bool one = true;
                std::size_t curve = 0;
                std::vector<double> ts;
                for (std::size_t i = 0; i + 1 < c.mVerts.size() && one; ++i) {
                    const auto it = curve_at.find(fc_edge_key(c.mVerts[i], c.mVerts[i + 1]));
                    if (it == curve_at.end() || (i > 0 && it->second.first != curve)) {
                        one = false;
                        break;
                    }
                    curve = it->second.first;
                    const detail::StepEdgeCurve& ec = (*pNetworkCurves)[curve];
                    const std::size_t at = it->second.second;
                    const bool fwd = ec.mVerts[at] == c.mVerts[i];
                    if (i == 0)
                        ts.push_back(ec.mParams[fwd ? at : at + 1]);
                    ts.push_back(ec.mParams[fwd ? at + 1 : at]);
                }
                for (std::size_t i = 2; one && i < ts.size(); ++i)
                    one = (ts[i] - ts[i - 1]) * (ts[1] - ts[0]) > 0.0;
                if (one && ts.size() >= 2 && ts.front() != ts.back()) {
                    const detail::BsplineCurve piece = detail::bspline_curve_segment(
                        (*pNetworkCurves)[curve].mCurve, ts.front(), ts.back());
                    const int n = static_cast<int>(piece.mPoles.size());
                    TColgp_Array1OfPnt poles(1, n);
                    for (int i = 0; i < n; ++i) {
                        const detail::PatchVec& q = piece.mPoles[static_cast<std::size_t>(i)];
                        poles.SetValue(i + 1, gp_Pnt(q[0], q[1], q[2]));
                    }
                    std::vector<double> kv;
                    std::vector<int> mv;
                    for (double k : piece.mKnots) {
                        if (!kv.empty() && k == kv.back())
                            ++mv.back();
                        else {
                            kv.push_back(k);
                            mv.push_back(1);
                        }
                    }
                    TColStd_Array1OfReal knots(1, static_cast<int>(kv.size()));
                    TColStd_Array1OfInteger mults(1, static_cast<int>(kv.size()));
                    for (std::size_t i = 0; i < kv.size(); ++i) {
                        knots.SetValue(static_cast<int>(i) + 1, kv[i]);
                        mults.SetValue(static_cast<int>(i) + 1, mv[i]);
                    }
                    const FcHandle<Geom_Curve> cv =
                        new Geom_BSplineCurve(poles, knots, mults, piece.mDegree);
                    fit_vertex(va, cv->Value(piece.T0()));
                    fit_vertex(vb, cv->Value(piece.T1()));
                    BRepBuilderAPI_MakeEdge me(cv, va, vb, piece.T0(), piece.T1());
                    if (me.IsDone()) {
                        c.mEdge = me.Edge();
                        c.mFromStart = true;
                        c.mOk = true;
                        ++out.mNumExact;
                        continue;
                    }
                }
            }
            if (!is_analytic(c.mGroupFwd) || !is_analytic(c.mGroupBwd)) {
                if (c.mVerts.size() != 2)
                    throw Standard_Failure("a facet chain spans several mesh edges");
                BRepBuilderAPI_MakeEdge me(va, vb);
                if (!me.IsDone())
                    throw Standard_Failure("degenerate facet edge");
                c.mEdge = me.Edge();
                c.mOk = true;
                continue;
            }
            const auto& s1 = rSurface[static_cast<std::size_t>(c.mGroupFwd)];
            const auto& s2 = rSurface[static_cast<std::size_t>(c.mGroupBwd)];
            std::vector<gp_Pnt> pts;
            for (std::int64_t v : c.mVerts)
                pts.push_back(fc_point(xyz, v));
            const gp_Pnt pa = BRep_Tool::Pnt(va), pb = BRep_Tool::Pnt(vb);

            // A fillet meets its walls tangentially and a thread flank is a
            // long B-spline: intersecting either is slow and unstable, so their
            // edges go straight to the approximation.
            auto swept = [&](std::int64_t G) {
                const std::int32_t t = rSkin.mPatches[static_cast<std::size_t>(G)].mType;
                return t == kFitCadFilletType || t == kFitCadHelicalType ||
                       t == kFitCadVertexBlendType || t == kFitCadChamferType;
            };
            // Two B-splines (a network patch and its freeform neighbours, which
            // meet it tangentially or nearly so, and in the freeform modes any
            // two) intersect slowly and unstably, and their intersection is
            // no truer to the skin than the chain moved onto both: those
            // edges go straight to the approximation.
            const auto bspline = [&](std::int64_t G) {
                return rSkin.mPatches[static_cast<std::size_t>(G)].mType == kFitCadBsplineType;
            };
            const bool network_pair =
                (FreeformApprox || network_group(c.mGroupFwd) || network_group(c.mGroupBwd)) &&
                bspline(c.mGroupFwd) && bspline(c.mGroupBwd);
            if (ExactEdges && !swept(c.mGroupFwd) && !swept(c.mGroupBwd) && !network_pair) {
                FcHandle<Geom_Curve> best;
                double best_dev = std::numeric_limits<double>::infinity();
                // The intersection branch nearest the chain, kept when nearer
                // than the best so far.
                auto intersect = [&](const FcHandle<Geom_Surface>& rA,
                                     const FcHandle<Geom_Surface>& rB) {
                    GeomAPI_IntSS inter(rA, rB, Precision::Confusion());
                    if (!inter.IsDone())
                        return;
                    for (int i = 1; i <= inter.NbLines(); ++i) {
                        FcHandle<Geom_Curve> cv = inter.Line(i);
                        // A trimmed circle or ellipse: its periodic basis, so
                        // either arc of a closed boundary can be cut from it.
                        while (cv->IsKind(STANDARD_TYPE(Geom_TrimmedCurve)))
                            cv = FcHandle<Geom_TrimmedCurve>::DownCast(cv)->BasisCurve();
                        double dev = 0.0;
                        for (const gp_Pnt& q : pts) {
                            GeomAPI_ProjectPointOnCurve pc(q, cv);
                            dev = std::max(dev, pc.NbPoints() > 0
                                                    ? pc.LowerDistance()
                                                    : std::numeric_limits<double>::infinity());
                            if (dev > best_dev)
                                break;
                        }
                        if (dev < best_dev) {
                            best_dev = dev;
                            best = cv;
                        }
                    }
                };
                intersect(s1, s2);
                bool extended_branch = false;
                if ((best.IsNull() || best_dev > Tol) && pExtended) {
                    // No branch along the chain: a freeform surface may stop
                    // short of it. Its extension's branch, when one fits.
                    const auto ext = [&](std::int64_t G) {
                        const std::size_t g = static_cast<std::size_t>(G);
                        return g < pExtended->size() ? (*pExtended)[g] : FcHandle<Geom_Surface>();
                    };
                    const FcHandle<Geom_Surface> e1 = ext(c.mGroupFwd), e2 = ext(c.mGroupBwd);
                    if (!e1.IsNull() || !e2.IsNull()) {
                        intersect(e1.IsNull() ? s1 : e1, e2.IsNull() ? s2 : e2);
                        extended_branch = !best.IsNull() && best_dev <= Tol;
                    }
                }
                if (!best.IsNull() && best_dev <= Tol) {
                    // A point inside the chain picks which arc of a closed
                    // curve it is; a single mesh edge has only its midpoint.
                    const gp_Pnt mid =
                        pts.size() > 2 ? pts[pts.size() / 2] : gp_Pnt((pa.XYZ() + pb.XYZ()) / 2.0);
                    GeomAPI_ProjectPointOnCurve ea(pa, best), eb(pb, best), em(mid, best);
                    if (ea.NbPoints() > 0 && eb.NbPoints() > 0 && em.NbPoints() > 0) {
                        double ua = ea.LowerDistanceParameter(), ub = eb.LowerDistanceParameter();
                        double um = em.LowerDistanceParameter();
                        bool ok = false, from_start = true;
                        double lo = 0.0, hi = 0.0;
                        if (best->IsPeriodic()) {
                            const double t = best->Period();
                            auto wrap = [t](double X) { return X - t * std::floor(X / t); };
                            const double span = wrap(ub - ua), mid = wrap(um - ua);
                            if (span > Precision::PConfusion()) {
                                if (mid < span) {
                                    lo = ua;
                                    hi = ua + span;
                                } else {
                                    lo = ub;
                                    hi = ub + wrap(ua - ub);
                                    from_start = false;
                                }
                                ok = hi - lo > Precision::PConfusion();
                            }
                        } else if (ua < um && um < ub) {
                            lo = ua;
                            hi = ub;
                            ok = true;
                        } else if (ub < um && um < ua) {
                            lo = ub;
                            hi = ua;
                            from_start = false;
                            ok = true;
                        }
                        // A branch of the extensions is kept only where it
                        // stays on both faces' own surfaces: the faces are
                        // trimmed on those, not on the extensions.
                        for (int k = 0; ok && extended_branch && k <= 16; ++k) {
                            const gp_Pnt at = best->Value(lo + (hi - lo) * k / 16.0);
                            for (const auto& srf : {s1, s2}) {
                                gp_Pnt on;
                                if (!fc_project(srf, at, on) || on.Distance(at) > Tol) {
                                    ok = false;
                                    break;
                                }
                            }
                        }
                        if (ok) {
                            const TopoDS_Vertex& v0 = from_start ? va : vb;
                            const TopoDS_Vertex& v1 = from_start ? vb : va;
                            fit_vertex(v0, best->Value(lo));
                            fit_vertex(v1, best->Value(hi));
                            BRepBuilderAPI_MakeEdge me(best, v0, v1, lo, hi);
                            if (me.IsDone()) {
                                c.mEdge = me.Edge();
                                c.mFromStart = from_start;
                                c.mOk = true;
                                ++out.mNumExact;
                                out.mNumExtended += extended_branch ? 1 : 0;
                                continue;
                            }
                        }
                    }
                }
            }

            // Approximate. A chain along a tangent junction (a fillet meeting a
            // wall) zigzags between two rows of mesh points, and the two
            // surfaces meet too shallowly for projection to settle it: smooth
            // the chain first (the binomial filter cancels an alternating
            // zigzag), move it onto both surfaces, fit a low-degree B-spline
            // and keep it only if it stays within the tolerance of both.
            const std::vector<FcHandle<Geom_Surface>> both{s1, s2};
            std::vector<gp_Pnt> q = pts;
            q.front() = pa;
            q.back() = pb;
            for (int pass = 0; pass < 2 && q.size() > 2; ++pass) {
                std::vector<gp_Pnt> r = q;
                for (std::size_t i = 1; i + 1 < q.size(); ++i)
                    r[i] = gp_Pnt((q[i - 1].XYZ() + 2.0 * q[i].XYZ() + q[i + 1].XYZ()) / 4.0);
                q.swap(r);
            }
            for (std::size_t i = 1; i + 1 < q.size(); ++i) {
                gp_Pnt on = q[i];
                const double d =
                    fc_onto_all(both, on, rSkin.mScale,
                                swept(c.mGroupFwd) || swept(c.mGroupBwd) ? 0.1 * Tol : 0.0);
                if (d <= Tol && on.Distance(pts[i]) <= Tol)
                    q[i] = on;
            }
            if (q.size() == 2) {
                BRepBuilderAPI_MakeEdge me(va, vb);
                if (!me.IsDone())
                    throw Standard_Failure("degenerate edge");
                c.mEdge = me.Edge();
            } else {
                TColgp_Array1OfPnt arr(1, static_cast<int>(q.size()));
                for (std::size_t i = 0; i < q.size(); ++i)
                    arr.SetValue(static_cast<int>(i) + 1, q[i]);
                const int n = static_cast<int>(q.size());
                const int deg_max = std::min(5, n - 1);
                const int deg_min = std::min(3, deg_max);
                const GeomAbs_Shape cont = deg_min >= 3   ? GeomAbs_C2
                                           : deg_min == 2 ? GeomAbs_C1
                                                          : GeomAbs_C0;
                GeomAPI_PointsToBSpline fit(arr, deg_min, deg_max, cont, 0.25 * Tol);
                if (!fit.IsDone())
                    throw Standard_Failure("no B-spline through the chain");
                const FcHandle<Geom_Curve> cv(fit.Curve());
                const double f0 = cv->FirstParameter(), f1 = cv->LastParameter();
                for (int k = 0; k <= 64; ++k) {
                    const gp_Pnt at = cv->Value(f0 + (f1 - f0) * k / 64.0);
                    for (const auto& srf : both) {
                        gp_Pnt on;
                        if (!fc_project(srf, at, on) || on.Distance(at) > Tol)
                            throw Standard_Failure("the chain's B-spline leaves a surface");
                    }
                }
                fit_vertex(va, cv->Value(f0));
                fit_vertex(vb, cv->Value(f1));
                BRepBuilderAPI_MakeEdge me(cv, va, vb, f0, f1);
                if (!me.IsDone())
                    throw Standard_Failure("no edge on the chain's B-spline");
                c.mEdge = me.Edge();
            }
            c.mFromStart = true;
            c.mOk = true;
            ++out.mNumApprox;
        } catch (const Standard_Failure& rErr) {
            log::debug("fit_cad: no edge on a chain of {} point(s): {}", c.mVerts.size(),
                       rErr.GetMessageString());
            mark_bad(c.mGroupFwd, FitCadDemotionReason::TrimFailed,
                     "an edge on its boundary failed");
            mark_bad(c.mGroupBwd, FitCadDemotionReason::TrimFailed,
                     "an edge on its boundary failed");
        }
    }

    // Around vertex V, after the boundary edge X -> V of group G (walked that
    // way by a facet of G): the vertex W of G's next boundary edge V -> W,
    // found by turning through G's facets around V. -1 when the fan breaks.
    auto ring_next = [&](std::size_t F, std::int64_t V) -> std::int64_t {
        const std::int64_t b = fc.mStart[F], n = fc.mStart[F + 1] - b;
        for (std::int64_t k = 0; k < n; ++k)
            if (fc.mNodes[static_cast<std::size_t>(b + k)] == V)
                return fc.mNodes[static_cast<std::size_t>(b + (k + 1) % n)];
        return -1;
    };
    auto walker = [&](std::int64_t U, std::int64_t V) -> std::int64_t {  // the facet walking U -> V
        auto it = rSkin.mEdgeFacets.find({std::min(U, V), std::max(U, V)});
        if (it == rSkin.mEdgeFacets.end())
            return -1;
        return U < V ? it->second[0] : it->second[1];
    };
    auto fan_next = [&](std::int64_t G, std::int64_t X, std::int64_t V) -> std::int64_t {
        std::int64_t f = walker(X, V);
        for (std::size_t guard = 0; f >= 0 && guard <= nf; ++guard) {
            const std::int64_t w = ring_next(static_cast<std::size_t>(f), V);
            if (w < 0)
                return -1;
            const std::int64_t across = walker(w, V);
            if (across < 0 || group_of(static_cast<std::size_t>(across)) != G)
                return w;
            f = across;
        }
        return -1;
    };

    // Wires and faces: each group's oriented chains, linked into loops.
    std::map<std::int64_t, std::vector<std::pair<std::size_t, bool>>>
        uses;  // group -> (chain, fwd)
    for (std::size_t i = 0; i < chains.size(); ++i) {
        uses[chains[i].mGroupFwd].emplace_back(i, true);
        uses[chains[i].mGroupBwd].emplace_back(i, false);
    }
    std::vector<std::pair<TopoDS_Face, std::int64_t>> faces;  // with their group
    for (const auto& kv : uses) {
        const std::int64_t g = kv.first;
        bool broken = false;
        for (const auto& u : kv.second)
            broken = broken || !chains[u.first].mOk;
        if (broken) {
            mark_bad(g, FitCadDemotionReason::TrimFailed, "an edge on its boundary failed");
            continue;
        }
        // start vertex -> unused oriented chains
        std::map<std::int64_t, std::vector<std::size_t>> from;
        auto head = [&](std::size_t K) {
            const auto& u = kv.second[K];
            const FcChain& c = chains[u.first];
            return u.second ? c.mVerts.front() : c.mVerts.back();
        };
        auto tail = [&](std::size_t K) {
            const auto& u = kv.second[K];
            const FcChain& c = chains[u.first];
            return u.second ? c.mVerts.back() : c.mVerts.front();
        };
        auto second = [&](std::size_t K) {
            const auto& u = kv.second[K];
            const FcChain& c = chains[u.first];
            return u.second ? c.mVerts[1] : c.mVerts[c.mVerts.size() - 2];
        };
        auto before_tail = [&](std::size_t K) {
            const auto& u = kv.second[K];
            const FcChain& c = chains[u.first];
            return u.second ? c.mVerts[c.mVerts.size() - 2] : c.mVerts[1];
        };
        for (std::size_t k = 0; k < kv.second.size(); ++k)
            from[head(k)].push_back(k);
        std::vector<std::uint8_t> used(kv.second.size(), 0);
        std::vector<TopoDS_Wire> wires;
        for (std::size_t k0 = 0; k0 < kv.second.size() && !broken; ++k0) {
            if (used[k0])
                continue;
            TopoDS_Wire w;
            builder.MakeWire(w);
            std::size_t k = k0;
            const std::int64_t first = head(k0);
            for (;;) {
                used[k] = 1;
                const auto& u = kv.second[k];
                const FcChain& c = chains[u.first];
                // Walked forwards along the mesh chain exactly when the edge's
                // orientation agrees with the walk.
                const bool along = u.second == c.mFromStart;
                builder.Add(w, along ? c.mEdge : TopoDS::Edge(c.mEdge.Reversed()));
                const std::int64_t end = tail(k);
                // The next chain leaves `end` along the next boundary edge of
                // this group around the facet fan, so a loop that touches
                // itself at a vertex (a pinch) stays two loops, never a figure
                // of eight.
                const std::int64_t w = fan_next(g, before_tail(k), end);
                std::size_t next = kv.second.size();
                for (std::size_t cand : from[end])
                    if (second(cand) == w) {
                        next = cand;
                        break;
                    }
                if (next == k0)
                    break;
                if (next == kv.second.size() || used[next] || end == first) {
                    broken = true;
                    break;
                }
                k = next;
            }
            w.Closed(Standard_True);
            wires.push_back(w);
        }
        if (broken) {
            mark_bad(g, FitCadDemotionReason::TrimFailed, "its boundary does not close into loops");
            continue;
        }
        try {
            if (is_analytic(g)) {
                TopoDS_Face f;
                const FcHandle<Geom_Surface>& surf = rSurface[static_cast<std::size_t>(g)];
                builder.MakeFace(f, surf, Precision::Confusion());
                for (const TopoDS_Wire& w : wires)
                    builder.Add(f, w);
                // Planes need no parameter-space curves; a sphere's loop around
                // its pole is left to the missing-seam repair, which closes it
                // with a degenerated edge.
                const std::int32_t type = rSkin.mPatches[static_cast<std::size_t>(g)].mType;
                if (type >= 2)
                    fc_set_pcurves(f, wires, surf, Tol);
                ShapeFix_Face fix(f);
                fix.SetPrecision(Precision::Confusion());
                fix.SetMaxTolerance(Tol);
                fix.FixAddNaturalBoundMode() = 0;
                fix.FixMissingSeamMode() = 1;
                fix.FixOrientationMode() = 1;
                // OpenCASCADE 7.6 can crash fixing the degenerated edge at a
                // cone's apex; a cone face here is trimmed away from it.
                if (type == 4)
                    fix.FixPeriodicDegeneratedMode() = 0;
                fix.Perform();
                faces.emplace_back(fix.Face(), g);
            } else {
                // A coplanar facet group: a plane through its seed facet,
                // bounded by every loop of the group (holes included).
                const std::size_t facet = static_cast<std::size_t>(g - np);
                const gp_Pnt p0 =
                    fc_point(xyz, fc.mNodes[static_cast<std::size_t>(fc.mStart[facet])]);
                const gp_Dir n(fc.mNormals[facet * 3], fc.mNormals[facet * 3 + 1],
                               fc.mNormals[facet * 3 + 2]);
                if (wires.size() == 1) {
                    BRepBuilderAPI_MakeFace mf(gp_Pln(p0, n), wires.front(), Standard_True);
                    if (!mf.IsDone())
                        throw Standard_Failure("no planar facet face");
                    faces.emplace_back(mf.Face(), g);
                } else {
                    TopoDS_Face f;
                    builder.MakeFace(f, new Geom_Plane(p0, n), Precision::Confusion());
                    for (const TopoDS_Wire& w : wires)
                        builder.Add(f, w);
                    ShapeFix_Face fix(f);
                    fix.SetPrecision(Precision::Confusion());
                    fix.SetMaxTolerance(Tol);
                    fix.FixAddNaturalBoundMode() = 0;
                    fix.FixOrientationMode() = 1;
                    fix.Perform();
                    faces.emplace_back(fix.Face(), g);
                }
            }
        } catch (const Standard_Failure&) {
            if (is_analytic(g))
                mark_bad(g, FitCadDemotionReason::TrimFailed,
                         "OpenCASCADE could not make its face");
            else
                for (std::size_t f = 0; f < nf; ++f)
                    if (flat_group[f] == g)
                        out.mUnmerge.push_back(f);
        }
    }
    // A patch with no boundary at all (a whole sphere or torus) is its
    // surface's natural face.
    {
        std::vector<std::uint8_t> bounded(rSkin.mPatches.size(), 0), present(bounded.size(), 0);
        for (const FcChain& c : chains)
            for (std::int64_t g : {c.mGroupFwd, c.mGroupBwd})
                if (is_analytic(g))
                    bounded[static_cast<std::size_t>(g)] = 1;
        for (std::size_t f = 0; f < nf; ++f)
            if (is_analytic(group_of(f)))
                present[static_cast<std::size_t>(group_of(f))] = 1;
        for (std::int64_t g = 0; g < np; ++g) {
            const std::size_t r = static_cast<std::size_t>(g);
            if (!present[r] || bounded[r])
                continue;
            const std::int32_t type = rSkin.mPatches[r].mType;
            if (type != 2 && type != 5) {
                mark_bad(g, FitCadDemotionReason::TrimFailed,
                         "it has no boundary, and its surface is not closed");
                continue;
            }
            BRepBuilderAPI_MakeFace mf(rSurface[r], Precision::Confusion());
            if (!mf.IsDone()) {
                mark_bad(g, FitCadDemotionReason::TrimFailed,
                         "OpenCASCADE could not make its face");
                continue;
            }
            faces.emplace_back(mf.Face(), g);
        }
    }
    for (std::uint8_t b : out.mBadRow)
        if (b)
            return out;
    if (!out.mUnmerge.empty())
        return out;

    // Sew the faces (their edges are shared already; sewing closes what the
    // face repair re-made) and close each shell into a solid.
    BRepBuilderAPI_Sewing sew(SewTol > 0.0 ? SewTol : Tol);
    for (const auto& f : faces)
        sew.Add(f.first);
    sew.Perform();
    const TopoDS_Shape sewed = sew.SewedShape();
    std::unordered_map<const void*, std::int64_t> group_of_face;
    for (const auto& f : faces) {
        const TopoDS_Shape m = sew.IsModified(f.first) ? sew.Modified(f.first) : f.first;
        for (TopExp_Explorer ex(m, TopAbs_FACE); ex.More(); ex.Next())
            group_of_face[ex.Current().TShape().get()] = f.second;
    }
    auto face_group = [&](const TopoDS_Shape& rFace) -> std::int64_t {
        auto it = group_of_face.find(rFace.TShape().get());
        return it == group_of_face.end() ? -1 : it->second;
    };
    log::debug("fit_cad: sewing {} face(s): {} free, {} multiple, {} degenerated edge(s)",
               faces.size(), sew.NbFreeEdges(), sew.NbMultipleEdges(), sew.NbDegeneratedShapes());
    if (sew.NbFreeEdges() > 0 || sew.NbMultipleEdges() > 0) {
        TopTools_IndexedDataMapOfShapeListOfShape anc;
        TopExp::MapShapesAndAncestors(sewed, TopAbs_EDGE, TopAbs_FACE, anc);
        bool any = false;
        auto blame = [&](const TopoDS_Shape& rEdge) {
            const int i = anc.FindIndex(rEdge);
            if (i == 0)
                return;
            for (TopTools_ListOfShape::Iterator it(anc.FindFromIndex(i)); it.More(); it.Next()) {
                const std::int64_t g = face_group(it.Value());
                if (is_analytic(g)) {
                    mark_bad(g, FitCadDemotionReason::SewFailed,
                             "sewing left a free or multiple edge on it");
                    any = true;
                }
            }
        };
        for (int i = 1; i <= sew.NbFreeEdges(); ++i)
            blame(sew.FreeEdge(i));
        for (int i = 1; i <= sew.NbMultipleEdges(); ++i)
            blame(sew.MultipleEdge(i));
        if (!any) {
            // Free edges among facets only: blame the analytic faces that
            // touch their ends, whose trimming moved those vertices.
            TopTools_IndexedDataMapOfShapeListOfShape vanc;
            TopExp::MapShapesAndAncestors(sewed, TopAbs_VERTEX, TopAbs_FACE, vanc);
            auto blame_ends = [&](const TopoDS_Shape& rEdge) {
                for (TopExp_Explorer vx(rEdge, TopAbs_VERTEX); vx.More(); vx.Next()) {
                    const int i = vanc.FindIndex(vx.Current());
                    if (i == 0)
                        continue;
                    for (TopTools_ListOfShape::Iterator it(vanc.FindFromIndex(i)); it.More();
                         it.Next()) {
                        const std::int64_t g = face_group(it.Value());
                        if (is_analytic(g)) {
                            mark_bad(g, FitCadDemotionReason::SewFailed,
                                     "sewing left a free edge beside it");
                            any = true;
                        }
                    }
                }
            };
            for (int i = 1; i <= sew.NbFreeEdges(); ++i)
                blame_ends(sew.FreeEdge(i));
            for (int i = 1; i <= sew.NbMultipleEdges(); ++i)
                blame_ends(sew.MultipleEdge(i));
        }
        if (!any)
            for (std::int64_t r = 0; r < np; ++r)
                if (rAnalytic[static_cast<std::size_t>(r)])
                    mark_bad(r, FitCadDemotionReason::SewFailed,
                             "sewing left free edges among facets");
        for (std::uint8_t b : out.mBadRow)
            if (b)
                return out;
    }

    std::vector<TopoDS_Solid> solids;
    // Sewing hands back a lone face (a whole sphere or torus) as it is: give
    // it a shell of its own.
    std::vector<TopoDS_Shell> loose;
    for (TopExp_Explorer ex(sewed, TopAbs_FACE, TopAbs_SHELL); ex.More(); ex.Next()) {
        TopoDS_Shell sh;
        builder.MakeShell(sh);
        builder.Add(sh, ex.Current());
        sh.Closed(BRep_Tool::IsClosed(sh));
        loose.push_back(sh);
    }
    for (const TopoDS_Shell& sh : loose) {
        ShapeFix_Solid sfs;
        const TopoDS_Solid s = sfs.SolidFromShell(sh);
        if (!s.IsNull())
            solids.push_back(s);
    }
    for (TopExp_Explorer ex(sewed, TopAbs_SHELL); ex.More(); ex.Next()) {
        ShapeFix_Solid sfs;
        const TopoDS_Solid s = sfs.SolidFromShell(TopoDS::Shell(ex.Current()));
        if (!s.IsNull())
            solids.push_back(s);
    }
    if (solids.size() == 1) {
        out.mShape = solids.front();
    } else {
        TopoDS_Compound cmp;
        builder.MakeCompound(cmp);
        for (const TopoDS_Solid& s : solids)
            builder.Add(cmp, s);
        out.mShape = cmp;
    }

    // Faces on one surface merge (a coplanar facet group into the analytic
    // plane beside it, two patches of one cylinder) and collinear edges join:
    // kept only when the solid stays valid with the same volume. A merged
    // face takes an analytic group over a facet group.
    std::unordered_map<const void*, std::int64_t> merged_group;
    bool merged = false;
    if (MergeFaces && !out.mShape.IsNull()) {
        try {
            ShapeUpgrade_UnifySameDomain usd(out.mShape, Standard_True, Standard_True,
                                             Standard_False);
            usd.Build();
            const TopoDS_Shape shape = usd.Shape();
            GProp_GProps before, after;
            BRepGProp::VolumeProperties(out.mShape, before);
            BRepGProp::VolumeProperties(shape, after);
            const double v0 = std::fabs(before.Mass());
            if (!shape.IsNull() && BRepCheck_Analyzer(shape).IsValid() &&
                std::fabs(std::fabs(after.Mass()) - v0) <= 1e-6 * std::max(v0, 1e-300)) {
                const FcHandle<BRepTools_History> history = usd.History();
                auto take = [&](const TopoDS_Shape& rFace, std::int64_t G) {
                    auto ins = merged_group.emplace(rFace.TShape().get(), G);
                    if (!ins.second && !is_analytic(ins.first->second) && is_analytic(G))
                        ins.first->second = G;
                };
                for (TopExp_Explorer ex(out.mShape, TopAbs_FACE); ex.More(); ex.Next()) {
                    const std::int64_t g = face_group(ex.Current());
                    const TopTools_ListOfShape& mods = history->Modified(ex.Current());
                    if (mods.IsEmpty() && !history->IsRemoved(ex.Current()))
                        take(ex.Current(), g);
                    for (TopTools_ListOfShape::Iterator it(mods); it.More(); it.Next())
                        take(it.Value(), g);
                }
                out.mShape = shape;
                merged = true;
            }
        } catch (const Standard_Failure& rErr) {
            log::debug("fit_cad: merging faces on one surface failed: {}", rErr.GetMessageString());
        }
    }

    // Per face, the group: by the face's own shape, else (a face a repair
    // rebuilt) by the facet nearest a point inside it.
    auto nearest_group = [&](const TopoDS_Face& rFace) -> std::int64_t {
        double u0 = 0.0, u1 = 0.0, v0 = 0.0, v1 = 0.0;
        BRepTools::UVBounds(rFace, u0, u1, v0, v1);
        const FcHandle<Geom_Surface> srf = BRep_Tool::Surface(rFace);
        if (srf.IsNull())
            return -1;
        const gp_Pnt at = srf->Value(0.5 * (u0 + u1), 0.5 * (v0 + v1));
        double best = std::numeric_limits<double>::infinity();
        std::int64_t g = -1;
        for (std::size_t f = 0; f < nf; ++f) {
            gp_XYZ c(0.0, 0.0, 0.0);
            for (std::int64_t k = fc.mStart[f]; k < fc.mStart[f + 1]; ++k)
                c += fc_point(xyz, fc.mNodes[static_cast<std::size_t>(k)]).XYZ();
            c /= static_cast<double>(fc.mStart[f + 1] - fc.mStart[f]);
            const double d = (c - at.XYZ()).SquareModulus();
            if (d < best) {
                best = d;
                g = group_of(f);
            }
        }
        return g;
    };
    for (TopExp_Explorer ex(out.mShape, TopAbs_FACE); ex.More(); ex.Next()) {
        const TopoDS_Face face = TopoDS::Face(ex.Current());
        out.mFaces.push_back(face);
        std::int64_t g = -1;
        if (merged) {
            auto it = merged_group.find(face.TShape().get());
            g = it == merged_group.end() ? -1 : it->second;
        } else {
            g = face_group(face);
        }
        out.mFaceGroup.push_back(g >= 0 ? g : nearest_group(face));
    }

    BRepCheck_Analyzer check(out.mShape);
    out.mValid = check.IsValid() && !solids.empty();
    if (!out.mValid) {
        bool any = false;
        // A face can pass the check alone and fail it in the solid (its wires
        // unorientable there): ask the solid's own analysis first.
        auto fails = [&](const TopoDS_Face& rFace) {
            const FcHandle<BRepCheck_Result> res = check.Result(rFace);
            if (!res.IsNull())
                for (BRepCheck_ListIteratorOfListOfStatus it(res->Status()); it.More(); it.Next())
                    if (it.Value() != BRepCheck_NoError)
                        return true;
            return !BRepCheck_Analyzer(rFace).IsValid();
        };
        for (std::size_t i = 0; i < out.mFaces.size(); ++i) {
            const std::int64_t g = out.mFaceGroup[i];
            if (g < 0 || !fails(out.mFaces[i]))
                continue;
            log::debug("fit_cad: face {} fails the B-rep check (shape kind:status){}", i,
                       fc_check_report(out.mFaces[i]));
            if (is_analytic(g)) {
                mark_bad(g, FitCadDemotionReason::BRepCheckFailed,
                         "its face fails the B-rep check");
            } else {
                // A coplanar group's polygon: its facets become faces of
                // their own, as the faceted writer writes them.
                for (std::size_t f = 0; f < nf; ++f)
                    if (flat_group[f] == g)
                        out.mUnmerge.push_back(f);
            }
            any = true;
        }
        if (!out.mUnmerge.empty())
            return out;
        if (!any)
            for (std::int64_t r = 0; r < np; ++r)
                if (rAnalytic[static_cast<std::size_t>(r)])
                    mark_bad(r, FitCadDemotionReason::BRepCheckFailed,
                             "the solid fails the B-rep check");
    }
    return out;
}

// The triangles of every face of a meshed shape, welded where faces share
// points exactly.
struct FcTriangulation {
    std::vector<double> mXyz;
    std::vector<std::int64_t> mTri;   // three per triangle
    std::vector<std::int64_t> mFace;  // per triangle: the face index
    bool mComplete = true;            // every face carried a triangulation
};

FcTriangulation fc_triangulate(const std::vector<TopoDS_Face>& rFaces) {
    FcTriangulation out;
    struct Key {
        std::array<std::uint64_t, 3> mBits;
        bool operator==(const Key& rO) const { return mBits == rO.mBits; }
    };
    struct KeyHash {
        std::size_t operator()(const Key& rK) const {
            std::uint64_t h = 1469598103934665603ull;
            for (std::uint64_t b : rK.mBits)
                h = (h ^ b) * 1099511628211ull;
            return static_cast<std::size_t>(h);
        }
    };
    std::unordered_map<Key, std::int64_t, KeyHash> id;
    for (std::size_t fi = 0; fi < rFaces.size(); ++fi) {
        const TopoDS_Face& face = rFaces[fi];
        TopLoc_Location loc;
        const FcHandle<Poly_Triangulation> tri = BRep_Tool::Triangulation(face, loc);
        if (tri.IsNull()) {
            out.mComplete = false;
            continue;
        }
        const gp_Trsf trsf = loc.Transformation();
        std::vector<std::int64_t> local(static_cast<std::size_t>(tri->NbNodes()) + 1, -1);
        for (int i = 1; i <= tri->NbNodes(); ++i) {
            const gp_Pnt p = tri->Node(i).Transformed(trsf);
            Key k;
            const double c[3] = {p.X() + 0.0, p.Y() + 0.0, p.Z() + 0.0};
            for (int j = 0; j < 3; ++j)
                std::memcpy(&k.mBits[static_cast<std::size_t>(j)], &c[j], sizeof(double));
            auto ins = id.emplace(k, static_cast<std::int64_t>(out.mXyz.size() / 3));
            if (ins.second)
                out.mXyz.insert(out.mXyz.end(), c, c + 3);
            local[static_cast<std::size_t>(i)] = ins.first->second;
        }
        const bool reversed = face.Orientation() == TopAbs_REVERSED;
        for (int t = 1; t <= tri->NbTriangles(); ++t) {
            int a = 0, b = 0, c = 0;
            tri->Triangle(t).Get(a, b, c);
            if (reversed)
                std::swap(b, c);
            const std::int64_t ia = local[static_cast<std::size_t>(a)];
            const std::int64_t ib = local[static_cast<std::size_t>(b)];
            const std::int64_t ic = local[static_cast<std::size_t>(c)];
            if (ia == ib || ib == ic || ia == ic)
                continue;
            out.mTri.insert(out.mTri.end(), {ia, ib, ic});
            out.mFace.push_back(static_cast<std::int64_t>(fi));
        }
    }
    return out;
}

Mesh fc_triangle_mesh(const std::vector<double>& rXyz, const std::vector<std::int64_t>& rTri) {
    Mesh m;
    const std::size_t np = rXyz.size() / 3, nt = rTri.size() / 3;
    NDArray pts = NDArray::Uninit(DType::Float64, {np, 3});
    std::copy(rXyz.begin(), rXyz.end(), pts.As<double>());
    NDArray conn = NDArray::Uninit(DType::Int64, {nt, 3});
    std::copy(rTri.begin(), rTri.end(), conn.As<std::int64_t>());
    m.AssignPoints(std::move(pts));
    m.AddCellBlock("triangle", std::move(conn));
    return m;
}

std::vector<double> fc_unsigned_distances(const Mesh& rTo, const std::vector<double>& rXyz) {
    const detail::TriangleSoup soup = detail::build_triangle_soup(rTo, "");
    SurfaceDistanceOptions opt;
    opt.mSign = SdfSign::Unsigned;
    opt.mWatertightCheck = SdfWatertightCheck::Off;
    const detail::DistanceQuery q = detail::build_distance_query(soup, opt);
    std::vector<detail::Vec3> pts(rXyz.size() / 3);
    for (std::size_t i = 0; i < pts.size(); ++i)
        pts[i] = {rXyz[i * 3], rXyz[i * 3 + 1], rXyz[i * 3 + 2]};
    const std::vector<detail::DistanceHit> hits = detail::query_distances(q, pts, opt);
    std::vector<double> d(hits.size());
    for (std::size_t i = 0; i < hits.size(); ++i)
        d[i] = std::fabs(hits[i].mSignedDistance);
    return d;
}

// Whether every shell of @p rShape is closed: no edge borders one face only
// (a seam, used twice by its face, counts as closed).
bool fc_closed(const TopoDS_Shape& rShape) {
    if (rShape.IsNull())
        return false;
    bool any = false;
    for (TopExp_Explorer sh(rShape, TopAbs_SHELL); sh.More(); sh.Next()) {
        any = true;
        TopTools_IndexedDataMapOfShapeListOfShape faces_of;
        TopExp::MapShapesAndAncestors(sh.Current(), TopAbs_EDGE, TopAbs_FACE, faces_of);
        for (int i = 1; i <= faces_of.Extent(); ++i) {
            const TopoDS_Edge& e = TopoDS::Edge(faces_of.FindKey(i));
            if (BRep_Tool::Degenerated(e))
                continue;
            const TopTools_ListOfShape& fl = faces_of.FindFromIndex(i);
            if (fl.Extent() >= 2)
                continue;
            if (fl.Extent() == 1 && BRep_Tool::IsClosed(e, TopoDS::Face(fl.First())))
                continue;
            return false;
        }
    }
    return any;
}

// The two-sided distance between each face's re-mesh and the skin: from the
// re-mesh's vertices and triangle centres to the skin, and from the skin's
// vertices and triangle centres to the solid, each skin sample charged to the
// face of the triangle nearest it. Field data `fit_cad:face_hausdorff` (per
// face: solid to skin, skin to solid) and `fit_cad:face_rms`, and the global
// maximum and root mean square on @p rResult.
void fc_face_deviation(const std::vector<double>& rSkinXyz,
                       const std::vector<std::int64_t>& rSkinTri, const FcTriangulation& rTri,
                       std::size_t NumFaces, FitCadResult& rResult) {
    SurfaceDistanceOptions opt;
    opt.mSign = SdfSign::Unsigned;
    opt.mWatertightCheck = SdfWatertightCheck::Off;
    // Every vertex, then every triangle's centre.
    auto samples = [](const std::vector<double>& rXyz, const std::vector<std::int64_t>& rTris) {
        std::vector<detail::Vec3> pts(rXyz.size() / 3);
        for (std::size_t i = 0; i < pts.size(); ++i)
            pts[i] = {rXyz[i * 3], rXyz[i * 3 + 1], rXyz[i * 3 + 2]};
        for (std::size_t t = 0; t + 2 < rTris.size(); t += 3) {
            detail::Vec3 c{0.0, 0.0, 0.0};
            for (std::size_t k = 0; k < 3; ++k)
                for (std::size_t j = 0; j < 3; ++j)
                    c[j] += rXyz[static_cast<std::size_t>(rTris[t + k]) * 3 + j] / 3.0;
            pts.push_back(c);
        }
        return pts;
    };
    auto hits_on = [&](const std::vector<double>& rToXyz, const std::vector<std::int64_t>& rToTri,
                       const std::vector<detail::Vec3>& rPts) {
        const Mesh to = fc_triangle_mesh(rToXyz, rToTri);
        const detail::TriangleSoup soup = detail::build_triangle_soup(to, "");
        const detail::DistanceQuery q = detail::build_distance_query(soup, opt);
        return detail::query_distances(q, rPts, opt);
    };
    std::vector<double> hmax(NumFaces * 2, 0.0), sq(NumFaces, 0.0);
    std::vector<std::size_t> cnt(NumFaces, 0);
    double gmax = 0.0, gsq = 0.0;
    std::size_t gcnt = 0;
    auto add = [&](std::int64_t Face, std::size_t Side, double D) {
        if (Face < 0 || static_cast<std::size_t>(Face) >= NumFaces || !std::isfinite(D))
            return;
        const std::size_t f = static_cast<std::size_t>(Face);
        hmax[f * 2 + Side] = std::max(hmax[f * 2 + Side], D);
        sq[f] += D * D;
        ++cnt[f];
        gmax = std::max(gmax, D);
        gsq += D * D;
        ++gcnt;
    };
    const std::size_t nv = rTri.mXyz.size() / 3, nt = rTri.mFace.size();
    // The solid against the skin; a vertex is charged to its first triangle's face.
    {
        const std::vector<detail::DistanceHit> hits =
            hits_on(rSkinXyz, rSkinTri, samples(rTri.mXyz, rTri.mTri));
        std::vector<std::int64_t> face_of(nv, -1);
        for (std::size_t t = 0; t < nt; ++t)
            for (std::size_t k = 0; k < 3; ++k) {
                const std::size_t v = static_cast<std::size_t>(rTri.mTri[t * 3 + k]);
                if (face_of[v] < 0)
                    face_of[v] = rTri.mFace[t];
            }
        for (std::size_t i = 0; i < nv; ++i)
            add(face_of[i], 0, std::fabs(hits[i].mSignedDistance));
        for (std::size_t t = 0; t < nt; ++t)
            add(rTri.mFace[t], 0, std::fabs(hits[nv + t].mSignedDistance));
    }
    // The skin against the solid, each sample charged to its nearest triangle's face.
    {
        const std::vector<detail::DistanceHit> hits =
            hits_on(rTri.mXyz, rTri.mTri, samples(rSkinXyz, rSkinTri));
        for (const detail::DistanceHit& h : hits)
            if (h.mSourceCell >= 0 && static_cast<std::size_t>(h.mSourceCell) < nt)
                add(rTri.mFace[static_cast<std::size_t>(h.mSourceCell)], 1,
                    std::fabs(h.mSignedDistance));
    }
    NDArray fh = NDArray::Uninit(DType::Float64, {NumFaces, 2});
    NDArray fr = NDArray::Uninit(DType::Float64, {NumFaces});
    for (std::size_t f = 0; f < NumFaces; ++f) {
        fh.As<double>()[f * 2] = hmax[f * 2];
        fh.As<double>()[f * 2 + 1] = hmax[f * 2 + 1];
        fr.As<double>()[f] = cnt[f] > 0 ? std::sqrt(sq[f] / static_cast<double>(cnt[f]))
                                        : std::numeric_limits<double>::quiet_NaN();
    }
    rResult.mMesh.AddFieldData(kFitCadFaceHausdorffName, std::move(fh));
    rResult.mMesh.AddFieldData(kFitCadFaceRmsName, std::move(fr));
    rResult.mHausdorffMax = gmax;
    rResult.mHausdorffRms = gcnt > 0 ? std::sqrt(gsq / static_cast<double>(gcnt))
                                     : std::numeric_limits<double>::quiet_NaN();
}

void fc_write_step(const TopoDS_Shape& rShape, const std::string& rPath) {
    const std::string base = rPath.substr(rPath.find_last_of("/\\") + 1);
    STEPControl_Writer writer;
    Interface_Static::SetCVal("write.step.schema", "AP214IS");
    Interface_Static::SetCVal("write.step.unit", "MM");
    if (writer.Transfer(rShape, STEPControl_AsIs) != IFSelect_RetDone)
        throw WriteError(std::string(kFcPrefix) +
                         "OpenCASCADE could not transfer the solid to STEP");
    // The product is otherwise named "Open CASCADE STEP translator <version>
    // <counter>", which makes the bytes depend on how many files this process
    // wrote before.
    const FcHandle<StepData_StepModel> model = writer.Model();
    for (int i = 1; i <= model->NbEntities(); ++i) {
        const FcHandle<StepBasic_Product> product =
            FcHandle<StepBasic_Product>::DownCast(model->Value(i));
        if (product.IsNull())
            continue;
        product->SetId(new TCollection_HAsciiString("meshio++ fit_cad"));
        product->SetName(new TCollection_HAsciiString("meshio++ fit_cad"));
    }
    APIHeaderSection_MakeHeader header(model);
    header.SetName(new TCollection_HAsciiString(base.c_str()));
    header.SetTimeStamp(new TCollection_HAsciiString("1970-01-01T00:00:00"));
    header.SetOriginatingSystem(new TCollection_HAsciiString("meshio++"));
    header.SetAuthorisation(new TCollection_HAsciiString(""));
    header.SetAuthorValue(1, new TCollection_HAsciiString(""));
    header.SetOrganizationValue(1, new TCollection_HAsciiString(""));
    if (writer.Write(rPath.c_str()) != IFSelect_RetDone)
        throw WriteError(std::string(kFcPrefix) + "cannot write '" + rPath + "'");
}

gp_Pnt fc_pnt(const detail::PatchVec& rV) { return gp_Pnt(rV[0], rV[1], rV[2]); }

// A B-spline face surface for a patch no primitive explains, or a null handle.
// @p rTris: the patch's triangles as local vertex indices; @p rPts: its points.
FcHandle<Geom_Surface> fc_freeform(const std::vector<gp_Pnt>& rPts,
                                   const std::vector<std::array<std::size_t, 3>>& rTris, double Tol,
                                   int Degree, std::int64_t MaxControl,
                                   FcRowInfo* pInfo = nullptr) {
    const std::size_t nt = rTris.size();
    if (rPts.size() < 4 || nt < 2) {
        if (pInfo)
            pInfo->mReason = static_cast<std::uint8_t>(FitCadDemotionReason::TooSmall);
        return {};
    }
    std::vector<detail::PatchVec> pts(rPts.size());
    for (std::size_t i = 0; i < rPts.size(); ++i)
        pts[i] = {rPts[i].X(), rPts[i].Y(), rPts[i].Z()};
    std::vector<detail::PatchVec> normals;
    std::vector<double> areas;
    detail::patch_normals(pts, rTris, normals, areas);
    // Boundary edges: those one triangle has.
    std::map<std::pair<std::size_t, std::size_t>, int> edge_count;
    for (const auto& t : rTris)
        for (int k = 0; k < 3; ++k) {
            const std::size_t a = t[static_cast<std::size_t>(k)];
            const std::size_t b = t[static_cast<std::size_t>((k + 1) % 3)];
            ++edge_count[{std::min(a, b), std::max(a, b)}];
        }
    std::vector<std::pair<std::size_t, std::size_t>> boundary;
    for (const auto& kv : edge_count)
        if (kv.second == 1)
            boundary.push_back(kv.first);

    // Why no surface came out, for the patch's record: every flattening
    // folds, until one does not.
    if (pInfo)
        pInfo->mReason = static_cast<std::uint8_t>(FitCadDemotionReason::ParamNotInjective);
    for (int kind = 0; kind < 2; ++kind) {
        detail::PatchFlat flat;
        if (!detail::patch_flatten(kind, pts, normals, areas, flat))
            continue;
        // A fold (triangles of both orientations) means this flattening is not
        // a parametrisation of the patch.
        if (detail::patch_flat_folds(flat, rTris)) {
            log::debug("fit_cad: a freeform patch folds over its {}", kind ? "cylinder" : "plane");
            continue;
        }
        if (pInfo)
            pInfo->mReason = static_cast<std::uint8_t>(FitCadDemotionReason::BudgetExceeded);
        double s0 = flat.mS[0], s1 = s0, t0 = flat.mT[0], t1 = t0;
        for (std::size_t i = 0; i < rPts.size(); ++i) {
            s0 = std::min(s0, flat.mS[i]);
            s1 = std::max(s1, flat.mS[i]);
            t0 = std::min(t0, flat.mT[i]);
            t1 = std::max(t1, flat.mT[i]);
        }
        if (!(s1 > s0) || !(t1 > t0))
            continue;
        // Bins over the triangles in (s, t), for locating grid points.
        const std::size_t nb = std::max<std::size_t>(
            1, static_cast<std::size_t>(std::sqrt(static_cast<double>(nt))));
        std::vector<std::vector<std::size_t>> bins(nb * nb);
        auto bin = [&](double S, double T, std::size_t& rI, std::size_t& rJ) {
            rI = std::min(nb - 1, static_cast<std::size_t>(std::max(0.0, (S - s0) / (s1 - s0) * nb)));
            rJ = std::min(nb - 1, static_cast<std::size_t>(std::max(0.0, (T - t0) / (t1 - t0) * nb)));
        };
        for (std::size_t t = 0; t < nt; ++t) {
            double ls = flat.mS[rTris[t][0]], hs = ls, lt = flat.mT[rTris[t][0]], ht = lt;
            for (std::size_t k = 1; k < 3; ++k) {
                ls = std::min(ls, flat.mS[rTris[t][k]]);
                hs = std::max(hs, flat.mS[rTris[t][k]]);
                lt = std::min(lt, flat.mT[rTris[t][k]]);
                ht = std::max(ht, flat.mT[rTris[t][k]]);
            }
            std::size_t i0, j0, i1, j1;
            bin(ls, lt, i0, j0);
            bin(hs, ht, i1, j1);
            for (std::size_t j = j0; j <= j1; ++j)
                for (std::size_t i = i0; i <= i1; ++i)
                    bins[j * nb + i].push_back(t);
        }
        // The offset at (S, T): interpolated inside the patch, else that of the
        // nearest boundary point.
        auto height = [&](double S, double T) {
            std::size_t bi, bj;
            bin(S, T, bi, bj);
            for (std::size_t t : bins[bj * nb + bi]) {
                const auto& tr = rTris[t];
                const double x0 = flat.mS[tr[0]], y0 = flat.mT[tr[0]];
                const double ax = flat.mS[tr[1]] - x0, ay = flat.mT[tr[1]] - y0;
                const double bx = flat.mS[tr[2]] - x0, by = flat.mT[tr[2]] - y0;
                const double det = ax * by - ay * bx;
                if (det == 0.0)
                    continue;
                const double l1 = ((S - x0) * by - (T - y0) * bx) / det;
                const double l2 = (ax * (T - y0) - ay * (S - x0)) / det;
                const double l0 = 1.0 - l1 - l2;
                const double eps = -1e-12;
                if (l0 >= eps && l1 >= eps && l2 >= eps)
                    return l0 * flat.mH[tr[0]] + l1 * flat.mH[tr[1]] + l2 * flat.mH[tr[2]];
            }
            double best = std::numeric_limits<double>::infinity(), h = 0.0;
            for (const auto& e : boundary) {
                const double x0 = flat.mS[e.first], y0 = flat.mT[e.first];
                const double dx = flat.mS[e.second] - x0, dy = flat.mT[e.second] - y0;
                const double len2 = dx * dx + dy * dy;
                double u = len2 > 0.0 ? ((S - x0) * dx + (T - y0) * dy) / len2 : 0.0;
                u = std::clamp(u, 0.0, 1.0);
                const double qx = x0 + u * dx - S, qy = y0 + u * dy - T;
                const double d2 = qx * qx + qy * qy;
                if (d2 < best) {
                    best = d2;
                    h = (1.0 - u) * flat.mH[e.first] + u * flat.mH[e.second];
                }
            }
            return h;
        };
        const double aspect = (t1 - t0) / (s1 - s0);
        const std::int64_t g_min = std::max<std::int64_t>(Degree + 1, 4);
        std::int64_t g = std::clamp<std::int64_t>(
            static_cast<std::int64_t>(std::sqrt(static_cast<double>(rPts.size()))), g_min,
            MaxControl);
        for (;; g = std::min(MaxControl, 2 * g)) {
            const std::int64_t nu =
                std::clamp<std::int64_t>(static_cast<std::int64_t>(std::lround(
                                             static_cast<double>(g) / std::sqrt(aspect))),
                                         g_min, MaxControl);
            const std::int64_t nv =
                std::clamp<std::int64_t>(static_cast<std::int64_t>(std::lround(
                                             static_cast<double>(g) * std::sqrt(aspect))),
                                         g_min, MaxControl);
            TColgp_Array2OfPnt grid(1, static_cast<int>(nu), 1, static_cast<int>(nv));
            for (std::int64_t j = 0; j < nv; ++j)
                for (std::int64_t i = 0; i < nu; ++i) {
                    const double S = s0 + (s1 - s0) * static_cast<double>(i) / static_cast<double>(nu - 1);
                    const double T = t0 + (t1 - t0) * static_cast<double>(j) / static_cast<double>(nv - 1);
                    grid.SetValue(static_cast<int>(i) + 1, static_cast<int>(j) + 1,
                                  fc_pnt(flat.At(S, T, height(S, T))));
                }
            try {
                // Iso-parametric: the grid's (i, j) are the surface's
                // parameters, so (s, t) maps to (u, v) affinely and a point is
                // checked by one evaluation instead of a projection.
                GeomAPI_PointsToBSplineSurface fit(grid, Approx_IsoParametric, std::min(3, Degree),
                                                   Degree, GeomAbs_C2, 0.1 * Tol);
                if (fit.IsDone()) {
                    FcHandle<Geom_Surface> surf(fit.Surface());
                    double u0 = 0.0, u1 = 0.0, v0 = 0.0, v1 = 0.0;
                    surf->Bounds(u0, u1, v0, v1);
                    // Within half the tolerance at every vertex, and within it
                    // at every facet centre (the surface point at the same
                    // (s, t), no nearer than the closest one): the face is held
                    // to the tolerance between vertices later, where a smooth
                    // surface is a facet's chordal error away from it, so a fit
                    // loose there would only be demoted. `worst` is the
                    // excess over each point's allowance, in tolerances.
                    double worst = 0.0;
                    for (std::size_t k = 0; k < rPts.size() + nt && worst <= 1.0; ++k) {
                        gp_Pnt p;
                        double S = 0.0, T = 0.0;
                        if (k < rPts.size()) {
                            p = rPts[k];
                            S = flat.mS[k];
                            T = flat.mT[k];
                        } else {
                            const auto& tr = rTris[k - rPts.size()];
                            S = (flat.mS[tr[0]] + flat.mS[tr[1]] + flat.mS[tr[2]]) / 3.0;
                            T = (flat.mT[tr[0]] + flat.mT[tr[1]] + flat.mT[tr[2]]) / 3.0;
                            p = gp_Pnt((rPts[tr[0]].XYZ() + rPts[tr[1]].XYZ() + rPts[tr[2]].XYZ()) /
                                       3.0);
                        }
                        const double u = u0 + (S - s0) / (s1 - s0) * (u1 - u0);
                        const double v = v0 + (T - t0) / (t1 - t0) * (v1 - v0);
                        const double allowance = k < rPts.size() ? 0.5 * Tol : Tol;
                        worst = std::max(worst, surf->Value(u, v).Distance(p) / allowance);
                    }
                    if (worst <= 1.0) {
                        if (pInfo) {
                            pInfo->mReason = 0;
                            pInfo->mParam =
                                kind ? FitCadParamMethod::Cylinder : FitCadParamMethod::Plane;
                            pInfo->mNumControl = static_cast<std::int64_t>(
                                fit.Surface()->NbUPoles() * fit.Surface()->NbVPoles());
                        }
                        return surf;
                    }
                    log::debug(
                        "fit_cad: a {}x{} freeform grid leaves the patch {} of its "
                        "allowance away",
                        nu, nv, worst);
                }
            } catch (const Standard_Failure&) {
            }
            if (g >= MaxControl)
                break;
        }
    }
    return {};
}

// Whether @p rS lies within half the tolerance of every vertex of a patch and
// within the tolerance of every facet centre, by projection: the test a
// fillet or thread face passes before it replaces a freeform fit. With
// @p LooseBoundary, a vertex on the patch's boundary may lie the whole
// tolerance away: a blend stays tangent to fitted neighbours, which are
// themselves held only to half of it, and the edge it shares with them is
// moved onto both and checked there. @p Cached projects through the
// surface's cached projector, which a long helical face needs: building a
// projector samples the whole surface.
bool fc_surface_fits(const FcHandle<Geom_Surface>& rS, const std::vector<gp_Pnt>& rPts,
                     const std::vector<std::array<std::size_t, 3>>& rTris, double Tol,
                     bool LooseBoundary = false, bool Cached = false) {
    auto within = [&](const gp_Pnt& rP, double Allowance) {
        if (Cached) {
            gp_Pnt foot;
            return fc_project(rS, rP, foot) && foot.Distance(rP) <= Allowance;
        }
        GeomAPI_ProjectPointOnSurf proj(rP, rS);
        return proj.IsDone() && proj.NbPoints() > 0 && proj.LowerDistance() <= Allowance;
    };
    std::vector<std::uint8_t> boundary(rPts.size(), 0);
    if (LooseBoundary) {
        std::map<std::pair<std::size_t, std::size_t>, int> uses;
        for (const auto& t : rTris)
            for (int k = 0; k < 3; ++k) {
                const std::size_t a = t[static_cast<std::size_t>(k)];
                const std::size_t b = t[static_cast<std::size_t>((k + 1) % 3)];
                ++uses[{std::min(a, b), std::max(a, b)}];
            }
        for (const auto& kv : uses)
            if (kv.second == 1)
                boundary[kv.first.first] = boundary[kv.first.second] = 1;
    }
    for (std::size_t i = 0; i < rPts.size(); ++i)
        if (!within(rPts[i], boundary[i] ? Tol : 0.5 * Tol))
            return false;
    for (const auto& t : rTris)
        if (!within(gp_Pnt((rPts[t[0]].XYZ() + rPts[t[1]].XYZ() + rPts[t[2]].XYZ()) / 3.0), Tol))
            return false;
    return true;
}

// A detail::BsplineSurface as OpenCASCADE's, its full clamped knot vectors as
// distinct knots and multiplicities.
FcHandle<Geom_Surface> fc_occt_bspline(const detail::BsplineSurface& rS) {
    TColgp_Array2OfPnt poles(1, static_cast<int>(rS.mNu), 1, static_cast<int>(rS.mNv));
    for (std::size_t j = 0; j < rS.mNv; ++j)
        for (std::size_t i = 0; i < rS.mNu; ++i) {
            const double* p = &rS.mPoles[3 * (i + rS.mNu * j)];
            poles.SetValue(static_cast<int>(i) + 1, static_cast<int>(j) + 1,
                           gp_Pnt(p[0], p[1], p[2]));
        }
    auto distinct = [](const std::vector<double>& rK, std::vector<double>& rKnots,
                       std::vector<int>& rMults) {
        for (double k : rK) {
            if (!rKnots.empty() && k == rKnots.back()) {
                ++rMults.back();
            } else {
                rKnots.push_back(k);
                rMults.push_back(1);
            }
        }
    };
    std::vector<double> ku, kv;
    std::vector<int> mu, mv;
    distinct(rS.mKnotsU, ku, mu);
    distinct(rS.mKnotsV, kv, mv);
    TColStd_Array1OfReal uk(1, static_cast<int>(ku.size())), vk(1, static_cast<int>(kv.size()));
    TColStd_Array1OfInteger um(1, static_cast<int>(mu.size())), vm(1, static_cast<int>(mv.size()));
    for (std::size_t i = 0; i < ku.size(); ++i) {
        uk.SetValue(static_cast<int>(i) + 1, ku[i]);
        um.SetValue(static_cast<int>(i) + 1, mu[i]);
    }
    for (std::size_t i = 0; i < kv.size(); ++i) {
        vk.SetValue(static_cast<int>(i) + 1, kv[i]);
        vm.SetValue(static_cast<int>(i) + 1, mv[i]);
    }
    return new Geom_BSplineSurface(poles, uk, vk, um, vm, rS.mDegreeU, rS.mDegreeV);
}

// The adaptive freeform face (v17.12.0, outside the legacy mode): the
// OpenCASCADE-free detail/freeform_adaptive fit (the parametrisation ladder,
// the noise allowance, budgeted refinement and re-parametrisation rungs, and
// surface-sampling validation), OpenCASCADE given only its poles and knots.
FcHandle<Geom_Surface> fc_freeform_adaptive(const std::vector<gp_Pnt>& rPts,
                                            const std::vector<std::array<std::size_t, 3>>& rTris,
                                            double Tol, const FitCadOptions& rOptions,
                                            FcRowInfo& rInfo) {
    std::vector<detail::PatchVec> pts(rPts.size());
    for (std::size_t i = 0; i < rPts.size(); ++i)
        pts[i] = {rPts[i].X(), rPts[i].Y(), rPts[i].Z()};
    const std::vector<detail::PatchTri> tris(rTris.begin(), rTris.end());
    detail::FreeformAdaptiveOptions o;
    o.mTolerance = Tol;
    o.mFreeformTolerance = rOptions.mFreeformTolerance;
    o.mSmoothing = rOptions.mFreeformSmoothing;
    o.mMaxPoles = static_cast<std::size_t>(rOptions.mMaxFreeformPoles);
    o.mMaxIterations = static_cast<int>(rOptions.mMaxFreeformIterations);
    o.mMaxTotalPoles = static_cast<std::size_t>(rOptions.mFreeformMaxControlPoints);
    o.mParameterization = static_cast<detail::ParamLadder>(rOptions.mFreeformParameterization);
    o.mValidate = rOptions.mValidateBySurfaceSampling;
    o.mThb = rOptions.mLocalRefine == FitCadLocalRefine::Thb;
    const detail::FreeformAdaptiveResult fit = detail::freeform_adaptive_fit(pts, tris, o);
    rInfo.mParam = fit.mParam;
    rInfo.mIterations = fit.mIterations;
    rInfo.mRefinements = fit.mRefinements;
    rInfo.mEscalations += fit.mEscalations;
    rInfo.mNumControl =
        static_cast<std::int64_t>(fit.mNumDofs > 0 ? fit.mNumDofs : fit.mSurface.NumPoles());
    rInfo.mMaxDist = fit.mMaxDist;
    rInfo.mRmsDist = fit.mRmsDist;
    rInfo.mReverseMaxDist = fit.mReverseMaxDist;
    rInfo.mAllowance = fit.mAllowance;
    rInfo.mReason = static_cast<std::uint8_t>(fit.mReason);
    log::debug(
        "fit_cad: adaptive freeform fit of a {}-triangle patch: {}, {} pole(s), vertex {} / "
        "centroid {} / samples {}, allowance {} (noise {}), {} escalation(s)",
        tris.size(), fit.mOk ? "accepted" : "rejected", fit.mSurface.NumPoles(), fit.mMaxDist,
        fit.mMaxCentroidDist, fit.mReverseMaxDist, fit.mAllowance, fit.mSigma, fit.mEscalations);
    if (!fit.mOk)
        return {};
    return fc_occt_bspline(fit.mSurface);
}

// The general freeform face (v17.11.0): the patch parametrised as a disc
// (detail/param_disc) and fitted by the iterative B-spline fit
// (detail/bspline_fit), OpenCASCADE given only the final poles and knots.
// Null when the patch is no disc, the fit misses its target, or the surface
// fails the test every freeform face passes (half the tolerance at a vertex,
// the tolerance at a facet centre, by projection).
FcHandle<Geom_Surface> fc_freeform_general(const std::vector<gp_Pnt>& rPts,
                                           const std::vector<std::array<std::size_t, 3>>& rTris,
                                           double Tol, const FitCadOptions& rOptions,
                                           std::int64_t& rRefinements, std::size_t& rPoles,
                                           FcRowInfo* pInfo = nullptr) {
    std::vector<detail::PatchVec> pts(rPts.size());
    for (std::size_t i = 0; i < rPts.size(); ++i)
        pts[i] = {rPts[i].X(), rPts[i].Y(), rPts[i].Z()};
    const std::vector<detail::PatchTri> tris(rTris.begin(), rTris.end());
    const detail::DiscParam par = detail::param_disc(pts, tris);
    if (!par.mOk) {
        log::debug("fit_cad: a {}-triangle patch is no disc (or folds in every parametrisation)",
                   tris.size());
        if (pInfo)
            pInfo->mReason =
                static_cast<std::uint8_t>(par.mDisc ? FitCadDemotionReason::ParamNotInjective
                                                    : FitCadDemotionReason::NonDisc);
        return {};
    }
    detail::BsplineFitOptions o;
    // A quarter of the tolerance at the vertices by default, so the face meets
    // the half-tolerance acceptance with room for its edges and the re-mesh;
    // a facet centre is held to the tolerance itself, its chordal sag allowed.
    o.mTolerance = rOptions.mFreeformTolerance > 0.0 ? rOptions.mFreeformTolerance : 0.25 * Tol;
    o.mCentroidTolerance = std::max(Tol, 2.0 * o.mTolerance);
    o.mSmoothing = rOptions.mFreeformSmoothing;
    o.mMaxPoles = static_cast<std::size_t>(rOptions.mMaxFreeformPoles);
    o.mMaxIterations = static_cast<int>(rOptions.mMaxFreeformIterations);
    const detail::BsplineFitResult fit = detail::bspline_fit_patch(pts, tris, par.mU, par.mV, o);
    log::debug(
        "fit_cad: general freeform fit of a {}-triangle patch: {} after {} iteration(s), "
        "{}x{} poles, vertex {} / centroid {} for a target of {}",
        tris.size(), fit.mOk ? "met" : "missed", fit.mIterations, fit.mSurface.mNu,
        fit.mSurface.mNv, fit.mMaxVertexError, fit.mMaxCentroidError, o.mTolerance);
    if (pInfo) {
        pInfo->mParam = par.mMethod == 1 ? FitCadParamMethod::Lscm : FitCadParamMethod::MeanValue;
        pInfo->mIterations = fit.mIterations;
        pInfo->mRefinements = fit.mRefinements;
        pInfo->mNumControl = static_cast<std::int64_t>(fit.mSurface.NumPoles());
        pInfo->mMaxDist = fit.mMaxVertexError;
        pInfo->mRmsDist = fit.mRmsError;
    }
    if (!fit.mOk) {
        if (pInfo)
            pInfo->mReason = static_cast<std::uint8_t>(FitCadDemotionReason::BudgetExceeded);
        return {};
    }
    const detail::BsplineSurface& s = fit.mSurface;
    TColgp_Array2OfPnt poles(1, static_cast<int>(s.mNu), 1, static_cast<int>(s.mNv));
    for (std::size_t j = 0; j < s.mNv; ++j)
        for (std::size_t i = 0; i < s.mNu; ++i) {
            const double* p = &s.mPoles[3 * (i + s.mNu * j)];
            poles.SetValue(static_cast<int>(i) + 1, static_cast<int>(j) + 1,
                           gp_Pnt(p[0], p[1], p[2]));
        }
    // The full clamped knot vectors as distinct knots and multiplicities.
    auto distinct = [](const std::vector<double>& rK, std::vector<double>& rKnots,
                       std::vector<int>& rMults) {
        for (double k : rK) {
            if (!rKnots.empty() && k == rKnots.back()) {
                ++rMults.back();
            } else {
                rKnots.push_back(k);
                rMults.push_back(1);
            }
        }
    };
    std::vector<double> ku, kv;
    std::vector<int> mu, mv;
    distinct(s.mKnotsU, ku, mu);
    distinct(s.mKnotsV, kv, mv);
    TColStd_Array1OfReal uk(1, static_cast<int>(ku.size())), vk(1, static_cast<int>(kv.size()));
    TColStd_Array1OfInteger um(1, static_cast<int>(mu.size())), vm(1, static_cast<int>(mv.size()));
    for (std::size_t i = 0; i < ku.size(); ++i) {
        uk.SetValue(static_cast<int>(i) + 1, ku[i]);
        um.SetValue(static_cast<int>(i) + 1, mu[i]);
    }
    for (std::size_t i = 0; i < kv.size(); ++i) {
        vk.SetValue(static_cast<int>(i) + 1, kv[i]);
        vm.SetValue(static_cast<int>(i) + 1, mv[i]);
    }
    FcHandle<Geom_Surface> surf =
        new Geom_BSplineSurface(poles, uk, vk, um, vm, s.mDegreeU, s.mDegreeV);
    // Through the cached projector: a projector per point samples the whole
    // surface each time.
    if (!fc_surface_fits(surf, rPts, rTris, Tol, false, true)) {
        log::debug("fit_cad: the general freeform surface strays from its patch on projection");
        if (pInfo)
            pInfo->mReason = static_cast<std::uint8_t>(FitCadDemotionReason::HausdorffExceeded);
        return {};
    }
    rRefinements = fit.mRefinements;
    rPoles = s.NumPoles();
    return surf;
}

// @p rS moved by @p Offset along its normal, as the equivalent elementary
// surface when there is one (a plane, or a cylinder of another radius), so
// two offsets of primitives intersect exactly.
FcHandle<Geom_Surface> fc_offset(const FcHandle<Geom_Surface>& rS, double Offset) {
    FcHandle<Geom_OffsetSurface> off = new Geom_OffsetSurface(rS, Offset);
    FcHandle<Geom_Surface> plain = off->Surface();
    return plain.IsNull() ? FcHandle<Geom_Surface>(off) : plain;
}

// The rolling-ball fillet of radius `rFit.mRadius` between @p rS1 and @p rS2,
// whose ball centres `rFit.mCentres` the patch says it rolls through: the pipe
// of that radius around the spine where the two surfaces, each offset by the
// radius towards the centres, meet. Null when no spine branch runs through
// the centres.
FcHandle<Geom_Surface> fc_fillet_surface(const FcHandle<Geom_Surface>& rS1,
                                         const FcHandle<Geom_Surface>& rS2,
                                         const detail::FilletFit& rFit, double Tol) {
    const double r = rFit.mRadius;
    if (rFit.mCentres.empty())
        return {};
    // The side of each surface the centres lie on, by majority.
    auto side = [&](const FcHandle<Geom_Surface>& rS) {
        int votes = 0;
        const std::size_t step = std::max<std::size_t>(1, rFit.mCentres.size() / 16);
        for (std::size_t i = 0; i < rFit.mCentres.size(); i += step) {
            const auto& c = rFit.mCentres[i];
            const gp_Pnt pc(c[0], c[1], c[2]);
            gp_Pnt foot;
            if (!fc_project(rS, pc, foot))
                continue;
            const gp_Vec n = fc_normal_at(rS, foot);
            votes += gp_Vec(foot, pc).Dot(n) >= 0.0 ? 1 : -1;
        }
        return votes >= 0 ? 1.0 : -1.0;
    };
    const FcHandle<Geom_Surface> o1 = fc_offset(rS1, side(rS1) * r);
    const FcHandle<Geom_Surface> o2 = fc_offset(rS2, side(rS2) * r);
    GeomAPI_IntSS inter(o1, o2, 1e-6 * Tol);
    if (!inter.IsDone() || inter.NbLines() == 0)
        return {};
    // The branch the centres lie closest to.
    FcHandle<Geom_Curve> spine;
    double best = std::numeric_limits<double>::infinity();
    std::vector<double> params;
    for (int i = 1; i <= inter.NbLines(); ++i) {
        const FcHandle<Geom_Curve> line = inter.Line(i);
        double sum = 0.0;
        std::vector<double> at;
        bool ok = true;
        for (const auto& c : rFit.mCentres) {
            GeomAPI_ProjectPointOnCurve proj(gp_Pnt(c[0], c[1], c[2]), line);
            if (proj.NbPoints() == 0) {
                ok = false;
                break;
            }
            sum += proj.LowerDistance();
            at.push_back(proj.LowerDistanceParameter());
        }
        const double mean = ok ? sum / static_cast<double>(rFit.mCentres.size()) : best;
        if (ok && mean < best) {
            best = mean;
            spine = line;
            params = std::move(at);
        }
    }
    if (spine.IsNull() || !(best <= Tol))
        return {};
    // The intersection comes back trimmed (an arc of an ellipse that only
    // nearly closes, between fitted surfaces) or as a closed B-spline: take the
    // periodic curve underneath, so a ring of a fillet is seen as one.
    FcHandle<Geom_Curve> periodic;
    if (!spine->IsPeriodic() || !spine->IsClosed()) {
        FcHandle<Geom_TrimmedCurve> tc = FcHandle<Geom_TrimmedCurve>::DownCast(spine);
        if (!tc.IsNull() && tc->BasisCurve()->IsPeriodic()) {
            periodic = tc->BasisCurve();
        } else if (spine->IsClosed()) {
            FcHandle<Geom_BSplineCurve> bc = FcHandle<Geom_BSplineCurve>::DownCast(spine->Copy());
            if (!bc.IsNull()) {
                bc->SetPeriodic();
                periodic = bc;
            }
        }
    }
    if (!periodic.IsNull()) {
        spine = periodic;
        params.clear();
        for (const auto& c : rFit.mCentres) {
            GeomAPI_ProjectPointOnCurve proj(gp_Pnt(c[0], c[1], c[2]), spine);
            if (proj.NbPoints() > 0)
                params.push_back(proj.LowerDistanceParameter());
        }
        if (params.empty())
            return {};
    }
    // Trim the spine to the centres, with a margin; a closed spine they go
    // most of the way round is kept whole.
    std::sort(params.begin(), params.end());
    const double f = spine->FirstParameter(), l = spine->LastParameter();
    bool whole = false;
    double u0 = params.front(), u1 = params.back();
    if (spine->IsPeriodic()) {
        const double per = spine->Period();
        double gap = params.front() + per - params.back(), gap_end = params.front();
        for (std::size_t i = 1; i < params.size(); ++i)
            if (params[i] - params[i - 1] > gap) {
                gap = params[i] - params[i - 1];
                gap_end = params[i];
            }
        if (gap < 0.25 * per) {
            whole = true;
        } else {
            u0 = gap_end;
            u1 = gap_end + per - gap;
        }
    }
    if (!whole) {
        const double margin = 0.1 * (u1 - u0) + 1e-9;
        u0 -= margin;
        u1 += margin;
        if (!spine->IsPeriodic()) {
            u0 = std::max(u0, f);
            u1 = std::min(u1, l);
        }
        if (!(u1 > u0))
            return {};
        spine = new Geom_TrimmedCurve(spine, u0, u1);
    }
    GeomFill_Pipe pipe(spine, r);
    pipe.Perform(0.1 * Tol, Standard_False);
    if (!pipe.IsDone())
        return {};
    FcHandle<Geom_Surface> out = pipe.Surface();
    // Around a whole closed spine the pipe comes back closed along it but not
    // periodic, which would cut a ring-shaped fillet's loops at the seam.
    FcHandle<Geom_BSplineSurface> bs = FcHandle<Geom_BSplineSurface>::DownCast(out);
    if (whole && !bs.IsNull() && bs->IsVClosed() && !bs->IsVPeriodic())
        bs->SetVPeriodic();
    if (whole && !bs.IsNull() && bs->IsUClosed() && !bs->IsUPeriodic())
        bs->SetUPeriodic();
    return out;
}

// The ruled helicoid of a thread flank, as a B-spline through a grid over the
// patch's radial and angular extent, slightly widened so the face trims inside it.
FcHandle<Geom_Surface> fc_helicoid_surface(const detail::HelicoidFit& rH, double Tol) {
    const double span = rH.mThetaMax - rH.mThetaMin;
    const double drho = 0.05 * (rH.mRhoMax - rH.mRhoMin), dth = 0.02 * span;
    const double rho0 = rH.mRhoMin - drho, rho1 = rH.mRhoMax + drho;
    const double th0 = rH.mThetaMin - dth, th1 = rH.mThetaMax + dth;
    const int nrho = 4;
    const int nth = std::clamp(static_cast<int>(std::ceil((th1 - th0) / (2.0 * kFcPi) * 48.0)) + 1,
                               8, 2048);
    TColgp_Array2OfPnt grid(1, nrho, 1, nth);
    for (int j = 0; j < nth; ++j)
        for (int i = 0; i < nrho; ++i) {
            const double rho = rho0 + (rho1 - rho0) * i / (nrho - 1);
            const double th = th0 + (th1 - th0) * j / (nth - 1);
            const detail::PatchVec p = rH.At(rho, th);
            grid.SetValue(i + 1, j + 1, gp_Pnt(p[0], p[1], p[2]));
        }
    GeomAPI_PointsToBSplineSurface fit(grid, Approx_IsoParametric, 1, 3, GeomAbs_C2, 0.1 * Tol);
    if (!fit.IsDone())
        return {};
    return FcHandle<Geom_Surface>(fit.Surface());
}

// The closest point of @p rS, through its cached projector, for the
// OpenCASCADE-free fits of `detail/feature_fit`.
detail::ClosestPoint fc_closest(const FcHandle<Geom_Surface>& rS) {
    return [rS](const detail::PatchVec& rX, detail::PatchVec& rFoot) {
        gp_Pnt f;
        if (!fc_project(rS, gp_Pnt(rX[0], rX[1], rX[2]), f))
            return false;
        rFoot = {f.X(), f.Y(), f.Z()};
        return true;
    };
}

// The unit normal of @p rS at the point of it nearest a point, for the
// OpenCASCADE-free fits of `detail/feature_fit`.
detail::NormalAt fc_normal(const FcHandle<Geom_Surface>& rS) {
    return [rS](const detail::PatchVec& rX, detail::PatchVec& rN) {
        const gp_Vec n = fc_normal_at(rS, gp_Pnt(rX[0], rX[1], rX[2]));
        const double len = n.Magnitude();
        if (!(len > 0.0))
            return false;
        rN = {n.X() / len, n.Y() / len, n.Z() / len};
        return true;
    };
}

// A chamfer or conic blend (roadmap 6.3.2), as a B-spline interpolating its
// cross-sections along the first rail -- straight chords, or the conic arcs
// of its weight law -- each widened by 15% at both ends and the rail carried
// a tenth of its length past the ends of an open one, so the face trims
// inside its surface. A closed rail gives a surface periodic along it. Null
// when a section misses.
FcHandle<Geom_Surface> fc_section_blend_surface(const detail::SectionBlendFit& rFit) {
    const double len = rFit.mLength;
    const bool closed = rFit.mClosed;
    const double t0 = closed ? 0.0 : -0.1 * len, t1 = closed ? len : 1.1 * len;
    const int n = std::clamp(static_cast<int>(2 * rFit.mRail1.size()) + 1, 8, 256);
    const int na = rFit.mChamfer ? 4 : 13;
    TColgp_Array2OfPnt grid(1, n, 1, na);
    for (int i = 0; i < n; ++i) {
        const double t = closed ? len * i / n : t0 + (t1 - t0) * i / (n - 1);
        for (int j = 0; j < na; ++j) {
            detail::PatchVec p{};
            if (!rFit.At(t, -0.15 + 1.3 * j / (na - 1), p))
                return {};
            grid.SetValue(i + 1, j + 1, gp_Pnt(p[0], p[1], p[2]));
        }
    }
    GeomAPI_PointsToBSplineSurface fit;
    fit.Interpolate(grid, Approx_IsoParametric, closed);
    if (!fit.IsDone())
        return {};
    return FcHandle<Geom_Surface>(fit.Surface());
}

// A vertex blend (roadmap 6.3.1): a plate surface through the patch's
// vertices that meets each neighbour tangentially along the boundary loop,
// approximated by a B-spline (widened by a tenth, so the face trims inside
// it). @p rLoop and @p rAcross are the loop's vertices and the neighbour
// across the edge leaving each; @p rSurface the neighbours' surfaces. Null
// when the plate cannot be built or a boundary point does not project.
FcHandle<Geom_Surface> fc_vertex_blend_surface(const std::vector<gp_Pnt>& rPts,
                                               const std::vector<std::size_t>& rLoop,
                                               const std::vector<std::int64_t>& rAcross,
                                               const std::vector<std::uint8_t>& rOnLoop,
                                               const std::vector<FcHandle<Geom_Surface>>& rSurface,
                                               double Tol) {
    GeomPlate_BuildPlateSurface plate(3, 15, 2, 1e-5, 0.1 * Tol, 0.01, 0.1);
    const std::size_t loop_step = std::max<std::size_t>(1, rLoop.size() / 96);
    for (std::size_t i = 0; i < rLoop.size(); i += loop_step) {
        const FcHandle<Geom_Surface>& s = rSurface[static_cast<std::size_t>(rAcross[i])];
        GeomAPI_ProjectPointOnSurf proj(rPts[rLoop[i]], s);
        if (!proj.IsDone() || proj.NbPoints() == 0)
            return {};
        double u = 0.0, v = 0.0;
        proj.LowerDistanceParameters(u, v);
        plate.Add(new GeomPlate_PointConstraint(u, v, s, 1, 0.1 * Tol, 0.01, 0.1));
    }
    std::vector<std::size_t> inner;
    for (std::size_t i = 0; i < rPts.size(); ++i)
        if (!rOnLoop[i])
            inner.push_back(i);
    const std::size_t inner_step = std::max<std::size_t>(1, inner.size() / 256);
    for (std::size_t i = 0; i < inner.size(); i += inner_step)
        plate.Add(new GeomPlate_PointConstraint(rPts[inner[i]], 0, 0.1 * Tol));
    plate.Perform();
    if (!plate.IsDone())
        return {};
    const FcHandle<GeomPlate_Surface> ps = plate.Surface();
    GeomPlate_MakeApprox approx(ps, 0.1 * Tol, 16, 8, 10.0 * Tol, 0, GeomAbs_C1, 1.1);
    const FcHandle<Geom_BSplineSurface> bs = approx.Surface();
    if (bs.IsNull())
        return {};
    return FcHandle<Geom_Surface>(bs);
}

// A rolling-ball blend of the radius `rFit.mRadius` gives along its rail, as
// a B-spline through its cross-sections. The balls touching both surfaces are
// solved along the rail (one that does not solve, where a neighbour ends, is
// left out), and smooth laws of the rail parameter are fitted through their
// centres and the directions to their feet, which irons out where a fitted
// neighbour wavers and carries the sections past the ends of an open rail (by
// a tenth of it, continuing the laws linearly) so the face trims inside its
// surface. Each grid row is the
// arc of the ball's characteristic circle from one foot to the other, widened
// by 15% at each end. A closed rail gives a surface periodic along it. Null
// when too few balls solve or no law follows them.
FcHandle<Geom_Surface> fc_rolling_ball_surface(const detail::RollingBallFit& rFit,
                                               const detail::ClosestPoint& rS1,
                                               const detail::ClosestPoint& rS2, double Tol,
                                               double Scale) {
    using detail::PatchVec;
    const double len = rFit.mLength;
    const bool closed = rFit.mClosed;
    double rmin = std::numeric_limits<double>::infinity();
    for (int i = 0; i <= 16; ++i)
        rmin = std::min(rmin, rFit.mRadius.At(len * i / 16.0));
    if (!(rmin > 0.0) || rFit.mT.empty())
        return {};
    // The balls.
    const int ns = std::clamp(static_cast<int>(std::ceil(len / (0.25 * rmin))) + 1, 16, 256);
    std::vector<double> at;
    std::array<std::vector<double>, 10> val;  // centre, foot directions 1 and 2, radius
    for (int i = 0; i < ns; ++i) {
        const double t = closed ? len * i / ns : len * i / (ns - 1);
        const detail::RollingBallSection s =
            detail::rolling_ball_section(rFit, rS1, rS2, t, rFit.RailTangent(t), Scale);
        if (!s.mOk)
            continue;
        at.push_back(t);
        for (int k = 0; k < 3; ++k) {
            val[static_cast<std::size_t>(k)].push_back(s.mCentre[static_cast<std::size_t>(k)]);
            val[static_cast<std::size_t>(3 + k)].push_back(
                (s.mFoot1[static_cast<std::size_t>(k)] - s.mCentre[static_cast<std::size_t>(k)]) /
                s.mRadius);
            val[static_cast<std::size_t>(6 + k)].push_back(
                (s.mFoot2[static_cast<std::size_t>(k)] - s.mCentre[static_cast<std::size_t>(k)]) /
                s.mRadius);
        }
        val[9].push_back(s.mRadius);
    }
    if (4 * at.size() < 3 * static_cast<std::size_t>(ns) || at.size() < 8)
        return {};
    // The laws: the lowest degree that follows every ball's centre and radius
    // within a quarter of the tolerance (a fitted neighbour, within half of
    // it, moves the balls that much). The directions to the feet only ride
    // along, smoothed: a ball resting on a nearly flat stretch of a fitted
    // face touches it where its normal happens to point, which wanders far
    // more than the face does, and the arcs are widened past their feet.
    const double l0 = closed ? 0.0 : at.front(), l1 = closed ? len : at.back();
    std::array<detail::FeatureLaw, 10> law;
    bool found = false;
    const int max_degree = std::min(closed ? 12 : 16, static_cast<int>(at.size()) / 3);
    for (int d = 1; d <= max_degree && !found; ++d) {
        bool ok = true;
        for (std::size_t q = 0; q < 10 && ok; ++q)
            ok = detail::feature_law_fit(at, val[q], closed, l0, l1, d, law[q]);
        if (!ok)
            break;
        // Across the spine only: along it, a ball's place is the rail's
        // arc length, which a wavering rail stretches, not geometry.
        double worst = 0.0;
        for (std::size_t i = 0; i < at.size(); ++i) {
            gp_Vec e(val[0][i] - law[0].At(at[i]), val[1][i] - law[1].At(at[i]),
                     val[2][i] - law[2].At(at[i]));
            gp_Vec tan(law[0].Slope(at[i]), law[1].Slope(at[i]), law[2].Slope(at[i]));
            if (tan.Magnitude() > 0.0) {
                tan.Normalize();
                e -= tan * e.Dot(tan);
            }
            worst = std::max({worst, e.Magnitude(), std::fabs(law[9].At(at[i]) - val[9][i])});
        }
        found = worst <= 0.25 * Tol;
    }
    if (!found)
        return {};
    // A law at T, continued linearly past the ends of an open rail.
    auto eval = [&](std::size_t Q, double T, double* pSlope) {
        const double c = closed ? T : std::clamp(T, l0, l1);
        const double slope = law[Q].Slope(c);
        if (pSlope)
            *pSlope = slope;
        return law[Q].At(c) + (T - c) * slope;
    };
    const double t0 = closed ? 0.0 : -0.1 * len, t1 = closed ? len : 1.1 * len;
    const int n = std::clamp(static_cast<int>(std::ceil((t1 - t0) / (0.25 * rmin))) + 1, 8, 256);
    const int na = 13;
    TColgp_Array2OfPnt grid(1, n, 1, na);
    for (int i = 0; i < n; ++i) {
        const double t = closed ? len * i / n : t0 + (t1 - t0) * i / (n - 1);
        double dc[3];
        const gp_Pnt c(eval(0, t, &dc[0]), eval(1, t, &dc[1]), eval(2, t, &dc[2]));
        const double r = eval(9, t, nullptr);
        gp_Vec u1(eval(3, t, nullptr), eval(4, t, nullptr), eval(5, t, nullptr));
        gp_Vec u2(eval(6, t, nullptr), eval(7, t, nullptr), eval(8, t, nullptr));
        const gp_Vec tan(dc[0], dc[1], dc[2]);
        if (!(r > 0.0) || u1.Magnitude() <= 0.0 || u2.Magnitude() <= 0.0 || tan.Magnitude() <= 0.0)
            return {};
        u1.Normalize();
        u2.Normalize();
        const gp_Pnt f1 = c.Translated(u1 * r), f2 = c.Translated(u2 * r);
        gp_Vec d(f1, f2);
        if (d.Magnitude() <= 0.0)
            return {};
        d.Normalize();
        // The characteristic circle: the ball cut by the plane through both
        // feet most nearly normal to the spine.
        gp_Vec nrm = tan - d * tan.Dot(d);
        if (nrm.Magnitude() <= 1e-12 * tan.Magnitude())
            return {};
        nrm.Normalize();
        const double delta = gp_Vec(c, f1).Dot(nrm);
        const double rho2 = r * r - delta * delta;
        if (!(rho2 > 0.0))
            return {};
        const double rho = std::sqrt(rho2);
        const gp_Pnt cc = c.Translated(nrm * delta);
        gp_Vec a(cc, f1);
        if (a.Magnitude() <= 0.0)
            return {};
        a.Normalize();
        const gp_Vec b2(cc, f2);
        gp_Vec b = b2 - a * b2.Dot(a);
        if (b.Magnitude() <= 1e-12 * rho)
            return {};
        b.Normalize();
        const double w = std::atan2(b2.Dot(b), b2.Dot(a));
        for (int j = 0; j < na; ++j) {
            const double phi = -0.15 * w + 1.3 * w * j / (na - 1);
            grid.SetValue(i + 1, j + 1,
                          cc.Translated(a * (rho * std::cos(phi)) + b * (rho * std::sin(phi))));
        }
    }
    // Interpolated rather than approximated: the grid points lie on the ball
    // family exactly, and a knot at every row and column keeps OpenCASCADE's
    // mesher (7.6 meshes a long, thin approximated surface with few spans far
    // coarser than the deflection asked of it) dense enough across the arcs.
    GeomAPI_PointsToBSplineSurface fit;
    fit.Interpolate(grid, Approx_IsoParametric, closed);
    if (!fit.IsDone())
        return {};
    return FcHandle<Geom_Surface>(fit.Surface());
}

// A helical surface of any profile, as a B-spline through a grid over the
// patch's profile and angular extent, slightly widened so the face trims
// inside it.
FcHandle<Geom_Surface> fc_helical_surface(const detail::HelicalFit& rH, double Tol) {
    const double span = rH.mThetaMax - rH.mThetaMin;
    const double ds = 0.05 * (rH.mS1 - rH.mS0), dth = 0.02 * span;
    const double s0 = rH.mS0 - ds, s1 = rH.mS1 + ds;
    const double th0 = rH.mThetaMin - dth, th1 = rH.mThetaMax + dth;
    const int ns = rH.mStraight ? 4 : 13;
    const int nth =
        std::clamp(static_cast<int>(std::ceil((th1 - th0) / (2.0 * kFcPi) * 48.0)) + 1, 8, 2048);
    TColgp_Array2OfPnt grid(1, ns, 1, nth);
    for (int j = 0; j < nth; ++j)
        for (int i = 0; i < ns; ++i) {
            const detail::PatchVec p =
                rH.At(s0 + (s1 - s0) * i / (ns - 1), th0 + (th1 - th0) * j / (nth - 1));
            grid.SetValue(i + 1, j + 1, gp_Pnt(p[0], p[1], p[2]));
        }
    GeomAPI_PointsToBSplineSurface fit(grid, Approx_IsoParametric, rH.mStraight ? 1 : 3,
                                       rH.mStraight ? 3 : 5, GeomAbs_C2, 0.1 * Tol);
    if (!fit.IsDone())
        return {};
    return FcHandle<Geom_Surface>(fit.Surface());
}

// One build; the solid also goes to @p pShape when given.
FitCadResult fc_occt_fit_cad_once(const Mesh& rMesh, const std::string& rStepPath,
                                  const FitCadOptions& rOptions, TopoDS_Shape* pShape) {
    const FcProjectorScope projectors;
    if (std::isnan(rOptions.mTolerance) || std::isinf(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kFcPrefix) + "tolerance must be finite");
    if (std::isnan(rOptions.mMeshDeflection) || std::isinf(rOptions.mMeshDeflection))
        throw std::invalid_argument(std::string(kFcPrefix) + "mesh_deflection must be finite");
    if ((rOptions.mTypes & ~kFitAllTypes) != 0)
        throw std::invalid_argument(std::string(kFcPrefix) +
                                    "types may name only plane, sphere, cylinder, cone and torus");
    if (rOptions.mMaxDemotionRounds < 0)
        throw std::invalid_argument(std::string(kFcPrefix) + "max_demotion_rounds must be >= 0");
    if (rOptions.mFreeformDegree < 3 || rOptions.mFreeformDegree > 8)
        throw std::invalid_argument(std::string(kFcPrefix) + "freeform_degree must be 3 to 8");
    if (rOptions.mFreeformMaxControl < 4 || rOptions.mFreeformMaxControl > 1024)
        throw std::invalid_argument(std::string(kFcPrefix) +
                                    "freeform_max_control must be 4 to 1024");

    detail::LibraryLock lock;
    FcQuietMessenger quiet;

    // The annotated skin: segment and fit first when the labels are missing.
    Mesh owned;
    FcnRefined refined;
    const Mesh& work = *fcn_prepare(rMesh, owned, rOptions, refined);

    FcSkin skin;
    skin.mFacets = detail::closed_facets(work, kFcPrefix);
    detail::ClosedFacets& fc = skin.mFacets;
    const std::size_t nf = fc.NumFacets();
    if (nf == 0)
        throw std::invalid_argument(std::string(kFcPrefix) + "no surface facets");
    if (fc.mNumSkipped > 0)
        log::warn("fit_cad: {} cell(s) that are not surface facets (lines, vertices) are ignored",
                  fc.mNumSkipped);

    // The fit table.
    for (const char* name :
         {"fit:patch_id", "fit:patch_type", "fit:patch_params", "fit:patch_tolerance"})
        if (!work.HasFieldData(name))
            throw std::invalid_argument(std::string(kFcPrefix) + "no field data named '" + name +
                                        "' (run fit_primitives first)");
    const NDArray& ids = work.FieldData("fit:patch_id");
    const NDArray& types = work.FieldData("fit:patch_type");
    const NDArray& params = work.FieldData("fit:patch_params");
    const NDArray& tols = work.FieldData("fit:patch_tolerance");
    // `m` grows when freeform patches are split; `num_input_patches` stays.
    std::size_t m = ids.Size();
    const std::size_t num_input_patches = m;
    if (types.Size() != m || tols.Size() != m || params.Size() != m * kFcParams)
        throw std::invalid_argument(std::string(kFcPrefix) +
                                    "the fit:patch_* arrays disagree on the number of patches");
    // After `refine_patches`, its verdict stands: a patch it found no
    // primitive for within the tolerance (a freeform or facets proxy) is
    // unexplained here too, even when `fit_primitives`, by its own measure,
    // accepts one (a torus's inner side for a saddle, say).
    std::unordered_map<std::int64_t, std::int64_t> refine_proxy;
    if (work.HasFieldData(kRefinePatchProxyName) && work.HasFieldData(kRefinePatchIdName)) {
        const NDArray& kinds = work.FieldData(kRefinePatchProxyName);
        const NDArray& labels = work.FieldData(kRefinePatchIdName);
        for (std::size_t i = 0; i < std::min(kinds.Size(), labels.Size()); ++i)
            refine_proxy.emplace(detail::read_int(labels, i), detail::read_int(kinds, i));
    }
    std::unordered_map<std::int64_t, std::int64_t> row_of;
    skin.mPatches.resize(m);
    for (std::size_t r = 0; r < m; ++r) {
        FcPatch& p = skin.mPatches[r];
        p.mLabel = detail::read_int(ids, r);
        p.mType = static_cast<std::int32_t>(detail::read_int(types, r));
        const auto verdict = refine_proxy.find(p.mLabel);
        if (verdict != refine_proxy.end() &&
            (verdict->second == static_cast<std::int64_t>(RefineProxy::Facets) ||
             verdict->second == static_cast<std::int64_t>(RefineProxy::Freeform)))
            p.mType = 0;
        bool finite = true;
        for (std::size_t k = 0; k < kFcParams; ++k) {
            p.mParams[k] = detail::read_double(params, r * kFcParams + k);
            finite = finite && std::isfinite(p.mParams[k]);
        }
        p.mTolerance = detail::read_double(tols, r);
        p.mAccepted = p.mType >= 1 && p.mType <= 5 && finite &&
                      (rOptions.mTypes & (1u << static_cast<unsigned>(p.mType))) != 0;
        row_of.emplace(p.mLabel, static_cast<std::int64_t>(r));
    }

    // Each facet's patch row.
    const std::vector<std::int64_t> bases = detail::block_bases(work);
    std::vector<std::int64_t> label(static_cast<std::size_t>(detail::total_cells(bases)), -1);
    for (std::size_t b = 0; b < work.NumCellBlocks(); ++b) {
        const NDArray& a = work.CellData("segment:patch", b);
        const std::size_t nc = work.Cells(b).NumCells();
        if (a.Size() != nc)
            throw std::invalid_argument(std::string(kFcPrefix) +
                                        "cell data 'segment:patch' must hold one label per cell");
        for (std::size_t c = 0; c < nc; ++c)
            label[static_cast<std::size_t>(bases[b]) + c] = detail::read_int(a, c);
    }
    skin.mRow.assign(nf, -1);
    for (std::size_t f = 0; f < nf; ++f) {
        auto it = row_of.find(label[static_cast<std::size_t>(fc.mCell[f])]);
        if (it != row_of.end())
            skin.mRow[f] = it->second;
    }

    // Edges and pieces.
    std::vector<std::size_t> parent(nf);
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    auto find = [&](std::size_t X) {
        while (parent[X] != X) {
            parent[X] = parent[parent[X]];
            X = parent[X];
        }
        return X;
    };
    skin.mEdgeFacets.reserve(fc.mNodes.size());
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t b = fc.mStart[f], s = fc.mStart[f + 1] - b;
        for (std::int64_t k = 0; k < s; ++k) {
            const std::int64_t u = fc.mNodes[static_cast<std::size_t>(b + k)];
            const std::int64_t v = fc.mNodes[static_cast<std::size_t>(b + (k + 1) % s)];
            auto ins = skin.mEdgeFacets.emplace(std::make_pair(std::min(u, v), std::max(u, v)),
                                                std::array<std::int64_t, 2>{-1, -1});
            ins.first->second[u < v ? 0 : 1] = static_cast<std::int64_t>(f);
            const std::int64_t other = ins.first->second[u < v ? 1 : 0];
            if (other >= 0) {
                const std::size_t ra = find(f), rb = find(static_cast<std::size_t>(other));
                if (ra != rb)
                    parent[std::max(ra, rb)] = std::min(ra, rb);
            }
        }
    }
    skin.mPiece.assign(nf, 0);
    std::vector<std::size_t> roots;
    for (std::size_t f = 0; f < nf; ++f) {
        const std::size_t r = find(f);
        if (r == f)
            roots.push_back(f);
        skin.mPiece[f] = static_cast<std::size_t>(std::lower_bound(roots.begin(), roots.end(), r) -
                                                  roots.begin());
    }
    skin.mNumPieces = roots.size();

    // Scale, skin volume and the skin as triangles.
    double lo[3] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[3] = {-lo[0], -lo[1], -lo[2]};
    std::vector<std::int64_t> skin_tri;
    std::vector<double> piece_volume(skin.mNumPieces, 0.0);
    for (std::size_t f = 0; f < nf; ++f) {
        const std::int64_t b = fc.mStart[f], s = fc.mStart[f + 1] - b;
        for (std::int64_t k = 0; k < s; ++k) {
            const std::size_t v =
                static_cast<std::size_t>(fc.mNodes[static_cast<std::size_t>(b + k)]);
            for (int j = 0; j < 3; ++j) {
                lo[j] = std::min(lo[j], fc.mXyz[v * 3 + static_cast<std::size_t>(j)]);
                hi[j] = std::max(hi[j], fc.mXyz[v * 3 + static_cast<std::size_t>(j)]);
            }
        }
        const std::int64_t a0 = fc.mNodes[static_cast<std::size_t>(b)];
        const gp_Pnt p0 = fc_point(fc.mXyz, a0);
        for (std::int64_t k = 1; k + 1 < s; ++k) {
            const std::int64_t a1 = fc.mNodes[static_cast<std::size_t>(b + k)];
            const std::int64_t a2 = fc.mNodes[static_cast<std::size_t>(b + k + 1)];
            skin_tri.insert(skin_tri.end(), {a0, a1, a2});
            const gp_Pnt p1 = fc_point(fc.mXyz, a1), p2 = fc_point(fc.mXyz, a2);
            piece_volume[skin.mPiece[f]] +=
                gp_Vec(p0.XYZ()).Dot(gp_Vec(p1.XYZ()).Crossed(gp_Vec(p2.XYZ()))) / 6.0;
        }
    }
    skin.mScale = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                            (hi[2] - lo[2]) * (hi[2] - lo[2]));
    const Mesh skin_mesh = fc_triangle_mesh(fc.mXyz, skin_tri);

    FitCadResult result;
    double tol = rOptions.mTolerance;
    if (!(tol > 0.0)) {
        double largest = 0.0;
        for (const FcPatch& p : skin.mPatches)
            if (p.mAccepted && std::isfinite(p.mTolerance))
                largest = std::max(largest, p.mTolerance);
        tol = largest > 0.0 ? 3.0 * largest : 1e-3 * skin.mScale;
    }
    if (!(tol > 0.0))
        throw std::invalid_argument(std::string(kFcPrefix) +
                                    "the skin has no extent to take a tolerance from");
    const double deflection =
        rOptions.mMeshDeflection > 0.0 ? rOptions.mMeshDeflection : 0.25 * tol;
    result.mTolerance = tol;

    // The fitted surfaces, each oriented so its normal points out of the skin.
    std::vector<std::vector<std::size_t>> facets_of(m);
    for (std::size_t f = 0; f < nf; ++f)
        if (skin.mRow[f] >= 0)
            facets_of[static_cast<std::size_t>(skin.mRow[f])].push_back(f);
    for (std::size_t r = 0; r < m; ++r)
        log::debug("fit_cad: patch {}: fitted type {}, {} facet(s)", skin.mPatches[r].mLabel,
                   skin.mPatches[r].mType, facets_of[r].size());
    std::vector<FcHandle<Geom_Surface>> surface(m);
    std::vector<std::uint8_t> analytic(m, 0);
    // What each row went through, for `FitCadResult::mFaces` (v17.12.0).
    std::vector<FcRowInfo> row_info(m);
    for (std::size_t r = 0; r < m; ++r) {
        if (!skin.mPatches[r].mAccepted || facets_of[r].empty())
            continue;
        std::vector<std::int64_t> verts;
        for (std::size_t f : facets_of[r])
            for (std::int64_t k = fc.mStart[f]; k < fc.mStart[f + 1]; ++k)
                verts.push_back(fc.mNodes[static_cast<std::size_t>(k)]);
        std::sort(verts.begin(), verts.end());
        verts.erase(std::unique(verts.begin(), verts.end()), verts.end());
        std::vector<gp_Pnt> pts;
        for (std::int64_t v : verts)
            pts.push_back(fc_point(fc.mXyz, v));
        try {
            FcHandle<Geom_Surface> s = fc_surface(skin.mPatches[r], pts);
            if (s.IsNull())
                continue;
            const std::size_t step = std::max<std::size_t>(1, facets_of[r].size() / 64);
            double agree = 0.0;
            for (std::size_t i = 0; i < facets_of[r].size(); i += step) {
                const std::size_t f = facets_of[r][i];
                gp_XYZ c(0.0, 0.0, 0.0);
                for (std::int64_t k = fc.mStart[f]; k < fc.mStart[f + 1]; ++k)
                    c += fc_point(fc.mXyz, fc.mNodes[static_cast<std::size_t>(k)]).XYZ();
                c /= static_cast<double>(fc.mStart[f + 1] - fc.mStart[f]);
                const gp_Vec n = fc_normal_at(s, gp_Pnt(c));
                agree += n.X() * fc.mNormals[f * 3] + n.Y() * fc.mNormals[f * 3 + 1] +
                         n.Z() * fc.mNormals[f * 3 + 2];
            }
            if (agree < 0.0)
                s = s->UReversed();
            surface[r] = s;
            analytic[r] = 1;
        } catch (const Standard_Failure&) {
        }
    }

    // Fillets and threads (roadmap 6.3.1): a patch no primitive explains that
    // runs between two analytic faces is tried as a constant-radius
    // rolling-ball fillet, then as one whose radius varies along it, and a
    // patch as a straight-flank thread about the axis of a cylinder or cone
    // beside it, then as a helical face of any profile about that axis or the
    // screw it fits itself, before the freeform fallback.
    auto patch_mesh = [&](std::size_t R, std::vector<gp_Pnt>& rPts,
                          std::vector<std::array<std::size_t, 3>>& rTris,
                          std::vector<detail::PatchVec>& rFlat) {
        std::unordered_map<std::int64_t, std::size_t> local;
        rPts.clear();
        rTris.clear();
        auto vid = [&](std::int64_t V) {
            auto ins = local.emplace(V, rPts.size());
            if (ins.second)
                rPts.push_back(fc_point(fc.mXyz, V));
            return ins.first->second;
        };
        for (std::size_t f : facets_of[R]) {
            const std::int64_t b = fc.mStart[f], n = fc.mStart[f + 1] - b;
            const std::size_t a0 = vid(fc.mNodes[static_cast<std::size_t>(b)]);
            for (std::int64_t k = 1; k + 1 < n; ++k)
                rTris.push_back({a0, vid(fc.mNodes[static_cast<std::size_t>(b + k)]),
                                 vid(fc.mNodes[static_cast<std::size_t>(b + k + 1)])});
        }
        rFlat.resize(rPts.size());
        for (std::size_t i = 0; i < rPts.size(); ++i)
            rFlat[i] = {rPts[i].X(), rPts[i].Y(), rPts[i].Z()};
    };
    // The global vertex of each local one of patch_mesh, in its first-seen order.
    auto patch_globals = [&](std::size_t R) {
        std::vector<std::int64_t> global;
        std::unordered_map<std::int64_t, std::size_t> seen;
        for (std::size_t f : facets_of[R])
            for (std::int64_t k = fc.mStart[f]; k < fc.mStart[f + 1]; ++k) {
                const std::int64_t g = fc.mNodes[static_cast<std::size_t>(k)];
                if (seen.emplace(g, global.size()).second)
                    global.push_back(g);
            }
        return global;
    };
    // Across the mesh edge (Ga, Gb) from patch R: the other row, or -1, and
    // that facet and R's own (-1 where there is none).
    auto row_across = [&](std::size_t R, std::int64_t Ga, std::int64_t Gb, std::int64_t* pMine,
                          std::int64_t* pTheirs) -> std::int64_t {
        auto it = skin.mEdgeFacets.find({std::min(Ga, Gb), std::max(Ga, Gb)});
        std::int64_t mine = -1, theirs = -1, row = -1;
        if (it != skin.mEdgeFacets.end())
            for (std::int64_t g : it->second) {
                if (g < 0)
                    continue;
                const std::int64_t q = skin.mRow[static_cast<std::size_t>(g)];
                if (q == static_cast<std::int64_t>(R)) {
                    mine = g;
                } else if (q >= 0) {
                    theirs = g;
                    row = q;
                }
            }
        if (pMine)
            *pMine = mine;
        if (pTheirs)
            *pTheirs = theirs;
        return row;
    };
    auto orient = [&](std::size_t R, FcHandle<Geom_Surface> S) {
        const std::size_t step = std::max<std::size_t>(1, facets_of[R].size() / 64);
        double agree = 0.0;
        for (std::size_t i = 0; i < facets_of[R].size(); i += step) {
            const std::size_t f = facets_of[R][i];
            gp_XYZ c(0.0, 0.0, 0.0);
            for (std::int64_t k = fc.mStart[f]; k < fc.mStart[f + 1]; ++k)
                c += fc_point(fc.mXyz, fc.mNodes[static_cast<std::size_t>(k)]).XYZ();
            c /= static_cast<double>(fc.mStart[f + 1] - fc.mStart[f]);
            const gp_Vec nrm = fc_normal_at(S, gp_Pnt(c));
            agree += nrm.X() * fc.mNormals[f * 3] + nrm.Y() * fc.mNormals[f * 3 + 1] +
                     nrm.Z() * fc.mNormals[f * 3 + 2];
        }
        return agree < 0.0 ? FcHandle<Geom_Surface>(S->UReversed()) : S;
    };
    auto is_primitive = [&](std::size_t R) {
        const std::int32_t t = skin.mPatches[R].mType;
        return analytic[R] && t >= 1 && t <= 5;
    };
    // The rows across a patch's boundary, by the length they share (longest
    // first), and the patch's perimeter.
    auto neighbours = [&](std::size_t R, double& rPerimeter) {
        std::map<std::int64_t, double> shared;
        rPerimeter = 0.0;
        for (std::size_t f : facets_of[R])
            for (std::int64_t k = fc.mStart[f]; k < fc.mStart[f + 1]; ++k) {
                const std::int64_t a = fc.mNodes[static_cast<std::size_t>(k)];
                const std::int64_t b = fc.mNodes[static_cast<std::size_t>(
                    k + 1 < fc.mStart[f + 1] ? k + 1 : fc.mStart[f])];
                auto it = skin.mEdgeFacets.find({std::min(a, b), std::max(a, b)});
                if (it == skin.mEdgeFacets.end())
                    continue;
                for (std::int64_t g : it->second) {
                    if (g < 0 || static_cast<std::size_t>(g) == f)
                        continue;
                    const std::int64_t q = skin.mRow[static_cast<std::size_t>(g)];
                    if (q == static_cast<std::int64_t>(R))
                        continue;
                    const double len = fc_point(fc.mXyz, a).Distance(fc_point(fc.mXyz, b));
                    rPerimeter += len;
                    if (q >= 0)
                        shared[q] += len;
                }
            }
        std::vector<std::pair<double, std::int64_t>> by_length;
        for (const auto& kv : shared)
            by_length.emplace_back(-kv.second, kv.first);
        std::sort(by_length.begin(), by_length.end());
        return by_length;
    };
    // What a fillet or thread face is, beyond its type, for the report.
    constexpr std::uint8_t kFcVariable = 1, kFcBlend = 2, kFcRounded = 4, kFcConic = 8,
                           kFcTapered = 16, kFcVariablePitch = 32, kFcGeneral = 64,
                           kFcAdaptive = 128;
    std::vector<std::uint8_t> feature(m, 0);
    auto accept = [&](std::size_t R, const FcHandle<Geom_Surface>& rS, std::int64_t Type,
                      std::uint8_t Feature) {
        surface[R] = orient(R, rS);
        analytic[R] = 1;
        skin.mPatches[R].mType = static_cast<std::int32_t>(Type);
        skin.mPatches[R].mAccepted = true;
        feature[R] = Feature;
    };
    for (std::size_t r = 0; (rOptions.mFillets || rOptions.mThreads) && r < m; ++r) {
        {
            bool finite = true;
            for (double x : skin.mPatches[r].mParams)
                finite = finite && std::isfinite(x);
            const std::int32_t type = skin.mPatches[r].mType;
            if (facets_of[r].empty() || analytic[r] || (type >= 1 && type <= 5 && finite))
                continue;
        }
        double perimeter = 0.0;
        const std::vector<std::pair<double, std::int64_t>> by_length = neighbours(r, perimeter);
        std::vector<gp_Pnt> pts;
        std::vector<std::array<std::size_t, 3>> tris;
        std::vector<detail::PatchVec> flat;
        patch_mesh(r, pts, tris, flat);
        try {
            if (rOptions.mFillets && by_length.size() >= 2) {
                const std::size_t q1 = static_cast<std::size_t>(by_length[0].second);
                const std::size_t q2 = static_cast<std::size_t>(by_length[1].second);
                const double along = -(by_length[0].first + by_length[1].first);
                if (is_primitive(q1) && is_primitive(q2) && along >= 0.6 * perimeter) {
                    const detail::FilletFit ff = detail::fillet_fit(
                        flat, tris, skin.mPatches[q1].mType, skin.mPatches[q1].mParams.data(),
                        skin.mPatches[q2].mType, skin.mPatches[q2].mParams.data(), 0.5 * tol,
                        skin.mScale);
                    if (ff.mOk && ff.mRms <= tol) {
                        FcHandle<Geom_Surface> fs =
                            fc_fillet_surface(surface[q1], surface[q2], ff, tol);
                        if (!fs.IsNull() && fc_surface_fits(fs, pts, tris, tol)) {
                            accept(r, fs, kFitCadFilletType, 0);
                            log::debug("fit_cad: patch {} is a fillet of radius {}",
                                       skin.mPatches[r].mLabel, ff.mRadius);
                            continue;
                        }
                    }
                    log::debug("fit_cad: patch {} is not a rolling-ball fillet (rms {})",
                               skin.mPatches[r].mLabel, ff.mRms);
                    // No one radius explains it: a radius varying along the
                    // rail, the face built through its cross-sections.
                    const detail::ClosestPoint s1 = fc_closest(surface[q1]);
                    const detail::ClosestPoint s2 = fc_closest(surface[q2]);
                    const detail::RollingBallFit rb =
                        detail::rolling_ball_fit(flat, tris, s1, s2, 0.5 * tol, skin.mScale, tol);
                    if (rb.mOk) {
                        FcHandle<Geom_Surface> fs =
                            fc_rolling_ball_surface(rb, s1, s2, tol, skin.mScale);
                        if (!fs.IsNull() && fc_surface_fits(fs, pts, tris, tol, false, true)) {
                            accept(r, fs, kFitCadFilletType, rb.mVariable ? kFcVariable : 0);
                            log::debug("fit_cad: patch {} is a fillet of radius {} to {}",
                                       skin.mPatches[r].mLabel, rb.mRadius.At(0.0),
                                       rb.mRadius.At(rb.mLength));
                            continue;
                        }
                        log::debug("fit_cad: patch {}: the fillet of varying radius {}",
                                   skin.mPatches[r].mLabel,
                                   fs.IsNull() ? "has no cross-section" : "strays from it");
                    } else {
                        log::debug("fit_cad: patch {} is no fillet of varying radius",
                                   skin.mPatches[r].mLabel);
                    }
                }
            }
            if (rOptions.mThreads) {
                // The axis of a neighbouring cylinder or cone, the one sharing
                // the most boundary first: the first whose axis the flank
                // winds about wins (a neighbour can be a cone of the end cap).
                const double* seed = nullptr;
                std::int32_t seed_type = 0;
                bool threaded = false;
                std::vector<std::pair<const double*, std::int32_t>> axes;
                for (const auto& nb : by_length) {
                    const std::size_t q = static_cast<std::size_t>(nb.second);
                    const std::int32_t t = skin.mPatches[q].mType;
                    if (!is_primitive(q) || (t != 3 && t != 4))
                        continue;
                    const double* pp = skin.mPatches[q].mParams.data();
                    axes.emplace_back(pp, t);
                    if (threaded)
                        continue;
                    if (seed == nullptr) {
                        seed = pp;
                        seed_type = t;
                    }
                    const detail::HelicoidFit hf =
                        detail::helicoid_fit(flat, tris, {pp[0], pp[1], pp[2]}, {pp[3], pp[4], pp[5]});
                    if (hf.mOk && hf.mMax <= 0.5 * tol) {
                        FcHandle<Geom_Surface> hs = fc_helicoid_surface(hf, tol);
                        if (!hs.IsNull() && fc_surface_fits(hs, pts, tris, tol)) {
                            // A straight flank beside a cone is cut on its taper:
                            // the same ruled helicoid, of an apparent pitch.
                            accept(r, hs, kFitCadHelicalType, t == 4 ? kFcTapered : 0);
                            threaded = true;
                            log::debug("fit_cad: patch {} is a thread flank of pitch {}",
                                       skin.mPatches[r].mLabel, hf.Pitch());
                        }
                    } else {
                        log::debug("fit_cad: patch {} is not a thread flank (max {})",
                                   skin.mPatches[r].mLabel, hf.mMax);
                    }
                    if (threaded) {
                        seed = pp;
                        seed_type = t;
                    }
                }
                // Any profile, about one of those axes or the screw the patch fits.
                if (!threaded) {
                    detail::HelicalFit hf;
                    for (const auto& ax : axes) {
                        const double* pp = ax.first;
                        const detail::PatchVec o{pp[0], pp[1], pp[2]};
                        const detail::PatchVec a{pp[3], pp[4], pp[5]};
                        hf = detail::helical_fit(flat, tris, &o, &a, tol);
                        if (hf.mOk) {
                            seed = pp;
                            seed_type = ax.second;
                            break;
                        }
                    }
                    if (!hf.mOk)
                        hf = detail::helical_fit(flat, tris, nullptr, nullptr, tol);
                    if (hf.mOk) {
                        FcHandle<Geom_Surface> hs = fc_helical_surface(hf, tol);
                        if (!hs.IsNull() && fc_surface_fits(hs, pts, tris, tol, false, true)) {
                            accept(r, hs, kFitCadHelicalType, hf.mStraight ? 0 : kFcRounded);
                            log::debug("fit_cad: patch {} is a helical face of pitch {}{}",
                                       skin.mPatches[r].mLabel, hf.Pitch(),
                                       hf.mStraight ? "" : ", rounded");
                        } else {
                            log::debug("fit_cad: patch {}: the helical face strays from it",
                                       skin.mPatches[r].mLabel);
                        }
                    }
                }
                // Tapered and variable-pitch threads (roadmap 6.3.3), tried
                // only where every constant screw failed: a straight flank
                // whose pitch changes along it, then any profile on a screw
                // whose pitch changes and that moves out along a taper.
                if (!analytic[r] && seed != nullptr) {
                    const detail::HelicoidFit hv = detail::helicoid_fit(
                        flat, tris, {seed[0], seed[1], seed[2]}, {seed[3], seed[4], seed[5]}, true);
                    if (hv.mOk && hv.mMax <= 0.5 * tol &&
                        std::fabs(hv.PitchAt(hv.mThetaMax) - hv.PitchAt(hv.mThetaMin)) > tol) {
                        FcHandle<Geom_Surface> hs = fc_helicoid_surface(hv, tol);
                        if (!hs.IsNull() && fc_surface_fits(hs, pts, tris, tol)) {
                            accept(r, hs, kFitCadHelicalType,
                                   static_cast<std::uint8_t>(kFcVariablePitch |
                                                             (seed_type == 4 ? kFcTapered : 0)));
                            log::debug("fit_cad: patch {} is a thread flank of pitch {} to {}",
                                       skin.mPatches[r].mLabel, hv.PitchAt(hv.mThetaMin),
                                       hv.PitchAt(hv.mThetaMax));
                        }
                    }
                }
                if (!analytic[r]) {
                    // The taper starts at the neighbouring cone's, if any: its
                    // axis points into the cone, where the radius grows.
                    const double taper0 = seed_type == 4 ? std::tan(seed[8] * kFcPi / 180.0) : 0.0;
                    detail::HelicalFit hg;
                    if (seed != nullptr) {
                        const detail::PatchVec o{seed[0], seed[1], seed[2]};
                        const detail::PatchVec a{seed[3], seed[4], seed[5]};
                        hg = detail::helical_fit(flat, tris, &o, &a, tol, true, taper0);
                    } else {
                        hg = detail::helical_fit(flat, tris, nullptr, nullptr, tol, true, 0.0);
                    }
                    if (hg.mOk) {
                        const double g0 = hg.Advance(hg.mThetaMin), g1 = hg.Advance(hg.mThetaMax);
                        const bool tapered = std::fabs(hg.mTaper * (g1 - g0)) > tol;
                        const bool varying =
                            std::fabs(hg.PitchAt(hg.mThetaMax) - hg.PitchAt(hg.mThetaMin)) > tol;
                        FcHandle<Geom_Surface> hs;
                        if (tapered || varying)
                            hs = fc_helical_surface(hg, tol);
                        if (!hs.IsNull() && fc_surface_fits(hs, pts, tris, tol, false, true)) {
                            accept(r, hs, kFitCadHelicalType,
                                   static_cast<std::uint8_t>((hg.mStraight ? 0 : kFcRounded) |
                                                             (tapered ? kFcTapered : 0) |
                                                             (varying ? kFcVariablePitch : 0)));
                            log::debug(
                                "fit_cad: patch {} is a helical face of pitch {} to {}, "
                                "taper {}{}",
                                skin.mPatches[r].mLabel, hg.PitchAt(hg.mThetaMin),
                                hg.PitchAt(hg.mThetaMax), hg.mTaper,
                                hg.mStraight ? "" : ", rounded");
                        }
                    }
                }
            }
            // Chamfers and conic blends (roadmap 6.3.2): a patch between two
            // analytic faces that no rolling ball and no thread explains,
            // its cross-sections straight or conic arcs tangent to both.
            if (!analytic[r] && rOptions.mFillets && by_length.size() >= 2) {
                const std::size_t q1 = static_cast<std::size_t>(by_length[0].second);
                const std::size_t q2 = static_cast<std::size_t>(by_length[1].second);
                const double along = -(by_length[0].first + by_length[1].first);
                const std::vector<std::int64_t> global = patch_globals(r);
                if (analytic[q1] && analytic[q2] && along >= 0.6 * perimeter &&
                    global.size() == flat.size()) {
                    // The rails are the boundary edges across from each neighbour.
                    const auto side = [&](std::size_t A, std::size_t B) {
                        const std::int64_t q =
                            row_across(r, global[A], global[B], nullptr, nullptr);
                        return q == static_cast<std::int64_t>(q1)   ? 1
                               : q == static_cast<std::int64_t>(q2) ? 2
                                                                    : 0;
                    };
                    const detail::SectionBlendFit sb = detail::section_blend_fit(
                        flat, tris, fc_closest(surface[q1]), fc_normal(surface[q1]),
                        fc_closest(surface[q2]), fc_normal(surface[q2]), tol, side);
                    if (!sb.mOk)
                        log::debug(
                            "fit_cad: patch {} is no chamfer or conic blend (rails of {} "
                            "and {} points)",
                            skin.mPatches[r].mLabel, sb.mRail1.size(), sb.mRail2.size());
                    if (sb.mOk) {
                        FcHandle<Geom_Surface> fs = fc_section_blend_surface(sb);
                        if (!fs.IsNull() &&
                            fc_surface_fits(fs, pts, tris, tol, !sb.mChamfer, true)) {
                            if (sb.mChamfer)
                                accept(r, fs, kFitCadChamferType, 0);
                            else
                                accept(r, fs, kFitCadFilletType, kFcConic);
                            log::debug("fit_cad: patch {} is a {}", skin.mPatches[r].mLabel,
                                       sb.mChamfer ? "chamfer" : "conic blend");
                        } else {
                            log::debug("fit_cad: patch {}: the {} {}", skin.mPatches[r].mLabel,
                                       sb.mChamfer ? "chamfer" : "conic blend",
                                       fs.IsNull() ? "has no surface" : "strays from it");
                        }
                    }
                }
            }
        } catch (const Standard_Failure& rErr) {
            log::debug("fit_cad: no fillet or thread face for patch {}: {}",
                       skin.mPatches[r].mLabel, rErr.GetMessageString());
        }
    }

    // Freeform faces (roadmap 6.3.2): a B-spline for each patch no primitive
    // explains, oriented like the analytic ones. A patch that folds over both
    // its flattenings is split, and the pieces, appended as new rows, are
    // fitted in turn by this same loop.
    std::int64_t num_freeform_tried = 0;
    std::int64_t next_label = 0;
    for (const FcPatch& p : skin.mPatches)
        next_label = std::max(next_label, p.mLabel + 1);
    std::vector<std::uint8_t> is_piece(m, 0), tile_depth(m, 0);
    // The general freeform stage's halvings so far, and what it built per row.
    std::vector<std::int64_t> general_depth(m, 0);
    std::vector<std::int64_t> general_refinements(m, 0), general_poles(m, 0);
    // The patch network (v17.12.0; network mode before the freeform faces,
    // auto mode for the patches they leave): the rows grouped into regions,
    // each laid out and fitted as one network, as the native backend does
    // (`fcn_networks`). A region with a layout gives every layout patch a row
    // of its own: one within the tolerance with its surface, the rest for the
    // freeform faces (`network`) or as facets with their reason (`auto`).
    // Neighbouring patches share their boundary curves exactly, so the edges
    // between them are built from surfaces that meet.
    std::vector<std::uint8_t> network_row(m, 0), network_tried(m, 0);
    std::vector<detail::StepEdgeCurve> network_curves;
    std::int64_t num_network_edges = 0, num_network_ev = 0;
    double network_g1 = 0.0;
    const auto occt_network = [&](const std::vector<std::size_t>& rRows, bool AsFacets) {
        std::vector<std::size_t> failed;
        std::vector<std::vector<std::size_t>> items;
        for (std::size_t r : rRows)
            items.push_back(facets_of[r]);
        for (const FcnRegionNetwork& reg : fcn_networks(fc, items, rOptions, tol)) {
            if (!reg.mLayout.mOk)
                continue;
            bool any = false;
            for (const detail::NetworkPatch& np : reg.mNet.mPatches)
                any = any || np.mOk;
            // Auto: a region the network fits nowhere keeps the adaptive
            // faces' rows and reasons.
            if (!any && AsFacets)
                continue;
            const std::size_t first = rRows[reg.mItems.front()];
            const double ptol = skin.mPatches[first].mTolerance;
            const std::int64_t escalations = fc_row(row_info, first).mEscalations + 1;
            for (std::size_t i : reg.mItems)
                facets_of[rRows[i]].clear();
            for (std::size_t i = 0; i < reg.mLayout.mPatches.size(); ++i) {
                const std::vector<std::size_t>& pf = reg.mPatchFacets[i];
                if (pf.empty())
                    continue;
                const detail::NetworkPatch& np = reg.mNet.mPatches[i];
                FcPatch piece;
                piece.mLabel = next_label++;
                piece.mType = 0;
                piece.mParams.fill(std::numeric_limits<double>::quiet_NaN());
                piece.mTolerance = ptol;
                const std::size_t row = skin.mPatches.size();
                for (std::size_t f : pf)
                    skin.mRow[f] = static_cast<std::int64_t>(row);
                skin.mPatches.push_back(piece);
                facets_of.push_back(pf);
                surface.emplace_back();
                analytic.push_back(0);
                is_piece.push_back(1);
                tile_depth.push_back(0);
                general_depth.push_back(0);
                network_row.resize(row + 1, 0);
                network_tried.resize(row + 1, 1);
                feature.resize(row + 1, 0);
                FcRowInfo& info = fc_row(row_info, row);
                info.mEscalations = escalations;
                FcHandle<Geom_Surface> s;
                if (np.mOk) {
                    try {
                        s = fc_occt_bspline(np.mSurface);
                    } catch (const Standard_Failure& rErr) {
                        log::debug("fit_cad: network patch {} has no OpenCASCADE surface: {}",
                                   piece.mLabel, rErr.GetMessageString());
                    }
                }
                if (s.IsNull()) {
                    if (AsFacets) {
                        info.mTried = true;
                        info.mReason =
                            static_cast<std::uint8_t>(np.mReason == FitCadDemotionReason::None
                                                          ? FitCadDemotionReason::HausdorffExceeded
                                                          : np.mReason);
                        failed.push_back(row);
                    }
                    continue;
                }
                surface[row] = orient(row, s);
                analytic[row] = 1;
                skin.mPatches[row].mType = static_cast<std::int32_t>(kFitCadBsplineType);
                skin.mPatches[row].mAccepted = true;
                feature[row] = kFcAdaptive;
                network_row[row] = 1;
                info.mTried = true;
                info.mParam = FitCadParamMethod::Square;
                info.mNumControl = static_cast<std::int64_t>(np.mSurface.NumPoles());
                info.mIterations = np.mIterations;
                info.mRefinements = np.mRefinements;
                info.mMaxDist = np.mMaxDist;
                info.mRmsDist = np.mRmsDist;
                info.mReverseMaxDist = np.mReverseMaxDist;
                info.mMaxG1AngleDeg = np.mMaxG1AngleDeg;
                info.mAllowance = np.mAllowance > tol ? np.mAllowance : 0.0;
            }
            if (any) {
                for (std::size_t c = 0; c < reg.mLayout.mChains.size(); ++c) {
                    if (reg.mNet.mChains[c].mPoles.size() < 2)
                        continue;
                    detail::StepEdgeCurve ec;
                    for (std::size_t v : reg.mLayout.mChains[c].mVerts)
                        ec.mVerts.push_back(reg.mGlobal[v]);
                    ec.mParams = reg.mNet.mChainParams[c];
                    ec.mCurve = reg.mNet.mChains[c];
                    network_curves.push_back(std::move(ec));
                }
                num_network_edges += static_cast<std::int64_t>(reg.mNet.mNumEdges);
                num_network_ev += static_cast<std::int64_t>(reg.mNet.mNumExtraordinary);
                network_g1 = std::max(network_g1, reg.mNet.mMaxG1AngleDeg);
            }
        }
        return failed;
    };
    const auto freeform_rows = [&](bool Tried) {
        std::vector<std::size_t> rows;
        for (std::size_t r = 0; r < skin.mPatches.size(); ++r) {
            bool finite = true;
            for (double x : skin.mPatches[r].mParams)
                finite = finite && std::isfinite(x);
            const std::int32_t type = skin.mPatches[r].mType;
            if (facets_of[r].empty() || analytic[r] || (type >= 1 && type <= 5 && finite))
                continue;
            if (!Tried || (r < row_info.size() && row_info[r].mTried))
                rows.push_back(r);
        }
        return rows;
    };
    // A row (a piece of a demoted face, a patch of a network) fitted as an
    // adaptive face; true when it is one.
    const auto refit_row = [&](std::size_t R, std::int64_t Escalations) {
        std::vector<gp_Pnt> ppts;
        std::vector<std::array<std::size_t, 3>> ptris;
        std::vector<detail::PatchVec> pflat;
        patch_mesh(R, ppts, ptris, pflat);
        FcRowInfo info;
        info.mEscalations = Escalations;
        FcHandle<Geom_Surface> srf;
        try {
            srf = fc_freeform_adaptive(ppts, ptris, tol, rOptions, info);
        } catch (const Standard_Failure& rErr) {
            log::debug("fit_cad: no adaptive face for patch {}: {}", skin.mPatches[R].mLabel,
                       rErr.GetMessageString());
        }
        info.mTried = true;
        if (srf.IsNull()) {
            if (!info.mReason)
                info.mReason = static_cast<std::uint8_t>(FitCadDemotionReason::HausdorffExceeded);
            fc_row(row_info, R) = info;
            return false;
        }
        fc_row(row_info, R) = info;
        surface[R] = orient(R, srf);
        analytic[R] = 1;
        skin.mPatches[R].mType = static_cast<std::int32_t>(kFitCadBsplineType);
        skin.mPatches[R].mAccepted = true;
        feature.resize(std::max(feature.size(), R + 1), 0);
        feature[R] = kFcAdaptive;
        return true;
    };
    // Facets halved across their longest extent: the pieces, and the facets
    // of the pieces under the size floor left out (false: no halving).
    const auto bisect_facets = [&](const std::vector<std::size_t>& rFacets,
                                   std::vector<std::vector<std::size_t>>& rPieces,
                                   std::vector<std::size_t>& rKept) {
        std::unordered_map<std::int64_t, std::size_t> local;
        std::vector<detail::PatchVec> flat;
        std::vector<std::array<std::size_t, 3>> tris;
        const auto vid = [&](std::int64_t V) {
            auto ins = local.emplace(V, flat.size());
            if (ins.second) {
                const std::size_t v = static_cast<std::size_t>(V);
                flat.push_back({fc.mXyz[v * 3], fc.mXyz[v * 3 + 1], fc.mXyz[v * 3 + 2]});
            }
            return ins.first->second;
        };
        for (std::size_t f : rFacets) {
            const std::int64_t b = fc.mStart[f], n = fc.mStart[f + 1] - b;
            const std::size_t a0 = vid(fc.mNodes[static_cast<std::size_t>(b)]);
            for (std::int64_t k = 1; k + 1 < n; ++k)
                tris.push_back({a0, vid(fc.mNodes[static_cast<std::size_t>(b + k)]),
                                vid(fc.mNodes[static_cast<std::size_t>(b + k + 1)])});
        }
        std::size_t num = 0;
        const std::vector<std::int64_t> piece =
            detail::patch_bisect(flat, tris, kFcTileMinTriangles, num);
        if (num < 2)
            return false;
        // A facet goes with its first triangle's piece.
        rPieces.assign(num, {});
        std::size_t t = 0;
        for (std::size_t f : rFacets) {
            const std::int64_t n = fc.mStart[f + 1] - fc.mStart[f];
            const std::int64_t at = piece[t];
            if (at >= 0)
                rPieces[static_cast<std::size_t>(at)].push_back(f);
            else
                rKept.push_back(f);
            t += static_cast<std::size_t>(std::max<std::int64_t>(n - 2, 0));
        }
        return true;
    };
    // An adaptive face the checks demote, cut across its longest extent into
    // halves and those again into quarters where they are large enough, each
    // piece fitted in turn (v17.12.0): the tiles of the general faces,
    // reached from the checks rather than the fit, two tile levels a round.
    const auto halve_and_refit = [&](std::size_t R) {
        if (R >= tile_depth.size() || tile_depth[R] >= kFcTileMaxDepth)
            return;
        std::vector<std::vector<std::size_t>> halves, pf;
        std::vector<std::size_t> kept;
        if (!bisect_facets(facets_of[R], halves, kept))
            return;
        std::uint8_t depth = static_cast<std::uint8_t>(tile_depth[R] + 1);
        for (std::vector<std::size_t>& h : halves) {
            std::vector<std::vector<std::size_t>> quarters;
            if (depth < kFcTileMaxDepth && bisect_facets(h, quarters, kept)) {
                for (std::vector<std::size_t>& q : quarters)
                    pf.push_back(std::move(q));
            } else {
                pf.push_back(std::move(h));
            }
        }
        if (pf.size() > halves.size())
            depth = static_cast<std::uint8_t>(depth + 1);
        facets_of[R] = std::move(kept);
        fc_row(row_info, R).mReason = static_cast<std::uint8_t>(FitCadDemotionReason::TooSmall);
        const std::int64_t escalations = fc_row(row_info, R).mEscalations + 1;
        for (std::vector<std::size_t>& facets : pf) {
            if (facets.empty())
                continue;
            FcPatch np;
            np.mLabel = next_label++;
            np.mType = 0;
            np.mParams.fill(std::numeric_limits<double>::quiet_NaN());
            np.mTolerance = skin.mPatches[R].mTolerance;
            const std::size_t row = skin.mPatches.size();
            for (std::size_t f : facets)
                skin.mRow[f] = static_cast<std::int64_t>(row);
            skin.mPatches.push_back(np);
            facets_of.push_back(std::move(facets));
            surface.emplace_back();
            analytic.push_back(0);
            is_piece.push_back(1);
            tile_depth.push_back(depth);
            general_depth.push_back(general_depth[R]);
            feature.resize(row + 1, 0);
            network_row.resize(row + 1, 0);
            network_tried.resize(row + 1, 1);
            ++result.mNumSplitPatches;
            refit_row(row, escalations);
        }
    };
    // A row with no face fitted again as it is; failing that, halved.
    const auto refit_or_halve = [&](std::size_t R) {
        if (!refit_row(R, fc_row(row_info, R).mEscalations + 1))
            halve_and_refit(R);
    };
    // The re-segmentation rung (v17.12.0): a row left as facets cut by its
    // normals (fcn_resegment) and each piece fitted as an adaptive face, at
    // most `mResegmentRounds` cuts deep; true when it was cut.
    std::vector<std::uint8_t> reseg_depth;
    // Per row, the early re-segmentation cuts in its line.
    std::vector<std::uint8_t> early_cut;
    const auto resegment_and_refit = [&](std::size_t R) {
        const std::uint8_t depth = R < reseg_depth.size() ? reseg_depth[R] : 0;
        if (depth >= rOptions.mResegmentRounds || facets_of[R].empty())
            return false;
        std::vector<std::vector<std::size_t>> cut = fcn_resegment(fc, facets_of[R], 2);
        if (cut.empty())
            return false;
        log::debug("fit_cad: patch {} re-segmented into {} piece(s)", skin.mPatches[R].mLabel,
                   cut.size());
        ++result.mNumResegmented;
        facets_of[R].clear();
        const std::int64_t escalations = fc_row(row_info, R).mEscalations + 1;
        for (std::vector<std::size_t>& facets : cut) {
            FcPatch np;
            np.mLabel = next_label++;
            np.mType = 0;
            np.mParams.fill(std::numeric_limits<double>::quiet_NaN());
            np.mTolerance = skin.mPatches[R].mTolerance;
            const std::size_t row = skin.mPatches.size();
            for (std::size_t f : facets)
                skin.mRow[f] = static_cast<std::int64_t>(row);
            skin.mPatches.push_back(np);
            facets_of.push_back(std::move(facets));
            surface.emplace_back();
            analytic.push_back(0);
            is_piece.push_back(1);
            tile_depth.push_back(tile_depth[R]);
            general_depth.push_back(general_depth[R]);
            feature.resize(row + 1, 0);
            network_row.resize(row + 1, 0);
            network_tried.resize(row + 1, 1);
            reseg_depth.resize(row + 1, 0);
            reseg_depth[row] = static_cast<std::uint8_t>(depth + 1);
            refit_row(row, escalations);
        }
        return true;
    };

    if (rOptions.mFreeform && rOptions.mFreeformMode == FitCadFreeformMode::Network)
        occt_network(freeform_rows(false), false);
    for (std::size_t r = 0; rOptions.mFreeform && r < skin.mPatches.size(); ++r) {
        bool finite = true;
        for (double x : skin.mPatches[r].mParams)
            finite = finite && std::isfinite(x);
        const std::int32_t type = skin.mPatches[r].mType;
        if (facets_of[r].empty() || analytic[r] || (type >= 1 && type <= 5 && finite))
            continue;
        ++num_freeform_tried;
        std::unordered_map<std::int64_t, std::size_t> local;
        std::vector<gp_Pnt> pts;
        std::vector<std::array<std::size_t, 3>> tris;
        std::vector<std::size_t> tri_facet;
        auto vid = [&](std::int64_t V) {
            auto ins = local.emplace(V, pts.size());
            if (ins.second)
                pts.push_back(fc_point(fc.mXyz, V));
            return ins.first->second;
        };
        for (std::size_t f : facets_of[r]) {
            const std::int64_t b = fc.mStart[f], n = fc.mStart[f + 1] - b;
            const std::size_t a0 = vid(fc.mNodes[static_cast<std::size_t>(b)]);
            for (std::int64_t k = 1; k + 1 < n; ++k) {
                tris.push_back({a0, vid(fc.mNodes[static_cast<std::size_t>(b + k)]),
                                vid(fc.mNodes[static_cast<std::size_t>(b + k + 1)])});
                tri_facet.push_back(f);
            }
        }
        fc_row(row_info, r).mTried = true;
        try {
            FcRowInfo info;
            info.mEscalations = fc_row(row_info, r).mEscalations;
            FcHandle<Geom_Surface> s = fc_freeform(pts, tris, tol, rOptions.mFreeformDegree,
                                                   rOptions.mFreeformMaxControl, &info);
            if (s.IsNull() && rOptions.mFreeformGeneral) {
                // The general freeform face, only where the plane and the
                // unrolled cylinder fall short, so every earlier face stays;
                // outside the legacy mode, the adaptive one (v17.12.0).
                const bool adaptive = rOptions.mFreeformMode != FitCadFreeformMode::Legacy;
                std::int64_t refinements = 0;
                std::size_t poles = 0;
                if (adaptive) {
                    const std::int64_t escalated = info.mEscalations;
                    info = FcRowInfo{};
                    info.mEscalations = escalated;
                    s = fc_freeform_adaptive(pts, tris, tol, rOptions, info);
                } else {
                    s = fc_freeform_general(pts, tris, tol, rOptions, refinements, poles, &info);
                }
                if (!s.IsNull()) {
                    feature.resize(skin.mPatches.size(), 0);
                    feature[r] = static_cast<std::uint8_t>(feature[r] |
                                                           (adaptive ? kFcAdaptive : kFcGeneral));
                    general_refinements.resize(skin.mPatches.size(), 0);
                    general_poles.resize(skin.mPatches.size(), 0);
                    general_refinements[r] = refinements;
                    general_poles[r] = static_cast<std::int64_t>(poles);
                } else if (general_depth[r] < rOptions.mMaxFreeformSplits && tris.size() >= 32) {
                    // Halve it along its normals' spread and fit the halves;
                    // outside the legacy mode, a patch that is not a disc is
                    // first cut along its topology into discs.
                    std::vector<detail::PatchVec> flat_pts(pts.size());
                    for (std::size_t i = 0; i < pts.size(); ++i)
                        flat_pts[i] = {pts[i].X(), pts[i].Y(), pts[i].Z()};
                    std::size_t num_pieces = 0;
                    std::vector<std::int64_t> piece;
                    const std::uint8_t cuts = r < early_cut.size() ? early_cut[r] : 0;
                    bool cut_here = false;
                    if (adaptive && cuts < rOptions.mResegmentRounds) {
                        // Across sharp edges, cut along them first (at most
                        // `mResegmentRounds` times in its line: the
                        // re-segmentation rung); the pieces keep the halvings
                        // left.
                        const auto cut = fcn_resegment(fc, facets_of[r], 2, true);
                        if (!cut.empty()) {
                            std::unordered_map<std::size_t, std::int64_t> of_facet;
                            for (std::size_t k = 0; k < cut.size(); ++k)
                                for (std::size_t f : cut[k])
                                    of_facet[f] = static_cast<std::int64_t>(k);
                            piece.resize(tris.size());
                            for (std::size_t tt = 0; tt < tris.size(); ++tt)
                                piece[tt] = of_facet[tri_facet[tt]];
                            num_pieces = cut.size();
                            cut_here = true;
                            ++result.mNumResegmented;
                        }
                    }
                    if (adaptive && num_pieces < 2)
                        piece = detail::patch_split_to_discs(flat_pts, tris, kFcTileMinTriangles,
                                                             num_pieces);
                    if (num_pieces < 2)
                        piece = detail::patch_bisect(flat_pts, tris, num_pieces);
                    if (num_pieces >= 2) {
                        // A facet goes with its first triangle's piece.
                        std::vector<std::vector<std::size_t>> piece_facets(num_pieces);
                        for (std::size_t tt = 0; tt < tris.size(); ++tt) {
                            if (tt > 0 && tri_facet[tt] == tri_facet[tt - 1])
                                continue;
                            piece_facets[static_cast<std::size_t>(piece[tt])].push_back(
                                tri_facet[tt]);
                        }
                        log::debug(
                            "fit_cad: patch {} has no general freeform face; halved "
                            "into {} piece(s)",
                            skin.mPatches[r].mLabel, num_pieces);
                        facets_of[r].clear();
                        // The split is a rung of the ladder its pieces inherit.
                        const std::int64_t escalations = info.mEscalations + (adaptive ? 1 : 0);
                        for (std::vector<std::size_t>& pf : piece_facets) {
                            FcPatch np;
                            np.mLabel = next_label++;
                            np.mType = 0;
                            np.mParams.fill(std::numeric_limits<double>::quiet_NaN());
                            np.mTolerance = skin.mPatches[r].mTolerance;
                            const std::int64_t row =
                                static_cast<std::int64_t>(skin.mPatches.size());
                            for (std::size_t f : pf)
                                skin.mRow[f] = row;
                            skin.mPatches.push_back(np);
                            facets_of.push_back(std::move(pf));
                            surface.emplace_back();
                            analytic.push_back(0);
                            is_piece.push_back(1);
                            // A piece is as deep in tiles as its parent (the
                            // depth was missing here before v17.12.0).
                            tile_depth.push_back(tile_depth[r]);
                            general_depth.push_back(general_depth[r] + (cut_here ? 0 : 1));
                            early_cut.resize(static_cast<std::size_t>(row) + 1, 0);
                            early_cut[static_cast<std::size_t>(row)] =
                                static_cast<std::uint8_t>(cuts + (cut_here ? 1 : 0));
                            fc_row(row_info, static_cast<std::size_t>(row)).mEscalations =
                                escalations;
                        }
                        result.mNumSplitPatches += static_cast<std::int64_t>(num_pieces);
                        continue;
                    }
                }
            }
            info.mTried = true;
            fc_row(row_info, r) = info;
            if (s.IsNull()) {
                if (!rOptions.mFreeformSplit)
                    continue;

                std::vector<detail::PatchVec> flat_pts(pts.size());
                for (std::size_t i = 0; i < pts.size(); ++i)
                    flat_pts[i] = {pts[i].X(), pts[i].Y(), pts[i].Z()};
                // A patch that folds over both flattenings is split into
                // pieces one of them parametrises; one that flattens, or a
                // piece, but that no single B-spline follows within the
                // tolerance is halved across its longest extent, and each half
                // is tried in turn, to a depth: a tile.
                std::size_t num_pieces = 0;
                std::vector<std::int64_t> piece;
                std::uint8_t piece_depth = tile_depth[r];
                const bool folds = !is_piece[r] && !detail::patch_flattenable(flat_pts, tris);
                if (folds) {
                    piece = detail::patch_split_until_flat(flat_pts, tris, {}, num_pieces);
                } else if (rOptions.mFreeformGeneral && tile_depth[r] < kFcTileMaxDepth) {
                    piece = detail::patch_bisect(flat_pts, tris, kFcTileMinTriangles, num_pieces);
                    piece_depth = static_cast<std::uint8_t>(tile_depth[r] + 1);
                }
                if (num_pieces == 0)
                    continue;
                // A facet goes with its first triangle's piece.
                std::vector<std::vector<std::size_t>> piece_facets(num_pieces);
                std::vector<std::size_t> kept;
                for (std::size_t t = 0; t < tris.size(); ++t) {
                    if (t > 0 && tri_facet[t] == tri_facet[t - 1])
                        continue;
                    if (piece[t] >= 0)
                        piece_facets[static_cast<std::size_t>(piece[t])].push_back(tri_facet[t]);
                    else
                        kept.push_back(tri_facet[t]);
                }
                log::debug("fit_cad: patch {} {}; split into {} piece(s), {} facet(s) left",
                           skin.mPatches[r].mLabel,
                           folds ? "folds over every flattening"
                                 : "has no single B-spline within the tolerance",
                           num_pieces, kept.size());
                facets_of[r] = std::move(kept);
                // What the pieces leave behind is under their size floor.
                fc_row(row_info, r).mReason =
                    static_cast<std::uint8_t>(FitCadDemotionReason::TooSmall);
                for (std::vector<std::size_t>& pf : piece_facets) {
                    FcPatch np;
                    np.mLabel = next_label++;
                    np.mType = 0;
                    np.mParams.fill(std::numeric_limits<double>::quiet_NaN());
                    np.mTolerance = skin.mPatches[r].mTolerance;
                    const std::int64_t row = static_cast<std::int64_t>(skin.mPatches.size());
                    for (std::size_t f : pf)
                        skin.mRow[f] = row;
                    skin.mPatches.push_back(np);
                    facets_of.push_back(std::move(pf));
                    surface.emplace_back();
                    analytic.push_back(0);
                    is_piece.push_back(1);
                    tile_depth.push_back(piece_depth);
                    general_depth.push_back(general_depth[r]);
                }
                result.mNumSplitPatches += static_cast<std::int64_t>(num_pieces);
                continue;
            }
            const std::size_t step = std::max<std::size_t>(1, facets_of[r].size() / 64);
            double agree = 0.0;
            for (std::size_t i = 0; i < facets_of[r].size(); i += step) {
                const std::size_t f = facets_of[r][i];
                gp_XYZ c(0.0, 0.0, 0.0);
                for (std::int64_t k = fc.mStart[f]; k < fc.mStart[f + 1]; ++k)
                    c += fc_point(fc.mXyz, fc.mNodes[static_cast<std::size_t>(k)]).XYZ();
                c /= static_cast<double>(fc.mStart[f + 1] - fc.mStart[f]);
                const gp_Vec nrm = fc_normal_at(s, gp_Pnt(c));
                agree += nrm.X() * fc.mNormals[f * 3] + nrm.Y() * fc.mNormals[f * 3 + 1] +
                         nrm.Z() * fc.mNormals[f * 3 + 2];
            }
            if (agree < 0.0)
                s = s->UReversed();
            surface[r] = s;
            analytic[r] = 1;
            skin.mPatches[r].mType = static_cast<std::int32_t>(kFitCadBsplineType);
            skin.mPatches[r].mAccepted = true;
        } catch (const Standard_Failure& rErr) {
            log::debug("fit_cad: no freeform face for patch {}: {}", skin.mPatches[r].mLabel,
                       rErr.GetMessageString());
        }
    }
    if (rOptions.mFreeform && rOptions.mFreeformMode == FitCadFreeformMode::Auto) {
        const std::vector<std::size_t> rows = freeform_rows(true);
        network_tried.resize(skin.mPatches.size(), 0);
        for (std::size_t r : rows)
            network_tried[r] = 1;
        // A patch of the network that misses is a cleaner piece than the
        // tile it came from: fitted as an adaptive face, or halved, first.
        if (!rows.empty())
            for (std::size_t r : occt_network(rows, true))
                refit_or_halve(r);
    }
    // What the adaptive faces and the network leave as facets is re-segmented
    // and fitted again; a piece that fails is cut again, up to the rounds.
    if (rOptions.mFreeform && rOptions.mFreeformMode != FitCadFreeformMode::Legacy)
        for (const std::size_t r0 : freeform_rows(true))
            for (std::vector<std::size_t> todo{r0}; !todo.empty();) {
                const std::size_t r = todo.back();
                todo.pop_back();
                if (analytic[r] || (r < network_row.size() && network_row[r]))
                    continue;
                const std::size_t from = skin.mPatches.size();
                if (resegment_and_refit(r))
                    for (std::size_t q = skin.mPatches.size(); q-- > from;)
                        todo.push_back(q);
            }
    m = skin.mPatches.size();
    if (num_freeform_tried > 0)
        log::debug("fit_cad: {} patch(es) tried as freeform B-spline faces", num_freeform_tried);
    feature.resize(m, 0);
    network_row.resize(m, 0);
    network_tried.resize(m, 0);

    // Blends (roadmap 6.3.1): a patch between two fitted faces, one of them a
    // freeform B-spline, tried as a rolling-ball blend of constant or varying
    // radius. The neighbours are taken as the freeform stage left them, so
    // the result does not depend on the order the patches are tried in.
    if (rOptions.mFillets) {
        const std::vector<FcHandle<Geom_Surface>> before = surface;
        const std::vector<std::uint8_t> analytic_before = analytic;
        std::vector<std::int32_t> type_before(m);
        for (std::size_t r = 0; r < m; ++r)
            type_before[r] = skin.mPatches[r].mType;
        auto fitted = [&](std::size_t Q) {
            return analytic_before[Q] && type_before[Q] >= 1 &&
                   type_before[Q] <= static_cast<std::int32_t>(kFitCadBsplineType);
        };
        for (std::size_t r = 0; r < m; ++r) {
            bool finite = true;
            for (double x : skin.mPatches[r].mParams)
                finite = finite && std::isfinite(x);
            const std::int32_t type = type_before[r];
            const bool freeform =
                analytic_before[r] && type == static_cast<std::int32_t>(kFitCadBsplineType);
            const bool unexplained = !analytic_before[r] && !(type >= 1 && type <= 5 && finite);
            // A network patch is a face of its own network, not a blend.
            if (facets_of[r].empty() || (!freeform && !unexplained) || network_row[r])
                continue;
            double perimeter = 0.0;
            const std::vector<std::pair<double, std::int64_t>> by_length = neighbours(r, perimeter);
            if (by_length.size() < 2)
                continue;
            const std::size_t q1 = static_cast<std::size_t>(by_length[0].second);
            const std::size_t q2 = static_cast<std::size_t>(by_length[1].second);
            const double along = -(by_length[0].first + by_length[1].first);
            const std::int32_t bs = static_cast<std::int32_t>(kFitCadBsplineType);
            if (!fitted(q1) || !fitted(q2) || along < 0.6 * perimeter ||
                (type_before[q1] != bs && type_before[q2] != bs))
                continue;
            // Nor is a gap left inside a network, between two of its patches.
            if (network_row[q1] && network_row[q2])
                continue;
            std::vector<gp_Pnt> pts;
            std::vector<std::array<std::size_t, 3>> tris;
            std::vector<detail::PatchVec> flat;
            patch_mesh(r, pts, tris, flat);
            try {
                const detail::ClosestPoint s1 = fc_closest(before[q1]);
                const detail::ClosestPoint s2 = fc_closest(before[q2]);
                const detail::RollingBallFit rb =
                    detail::rolling_ball_fit(flat, tris, s1, s2, 0.5 * tol, skin.mScale, tol);
                if (!rb.mOk) {
                    log::debug("fit_cad: patch {} is no rolling-ball blend",
                               skin.mPatches[r].mLabel);
                    continue;
                }
                FcHandle<Geom_Surface> fs = fc_rolling_ball_surface(rb, s1, s2, tol, skin.mScale);
                if (fs.IsNull() || !fc_surface_fits(fs, pts, tris, tol, true, true)) {
                    log::debug("fit_cad: patch {}: the blend {}", skin.mPatches[r].mLabel,
                               fs.IsNull() ? "has no cross-section" : "strays from it");
                    continue;
                }
                accept(r, fs, kFitCadFilletType,
                       static_cast<std::uint8_t>(kFcBlend | (rb.mVariable ? kFcVariable : 0)));
                log::debug("fit_cad: patch {} is a blend of radius {} to {}",
                           skin.mPatches[r].mLabel, rb.mRadius.At(0.0), rb.mRadius.At(rb.mLength));
            } catch (const Standard_Failure& rErr) {
                log::debug("fit_cad: no blend face for patch {}: {}", skin.mPatches[r].mLabel,
                           rErr.GetMessageString());
            }
        }
    }

    // Vertex blends (roadmap 6.3.1): a patch no primitive explains (or only
    // a freeform face does) whose boundary is one loop around three or more
    // fitted faces, meeting each along one run without a crease -- the
    // corner where blends meet, set back from each other -- tried as a plate
    // surface tangent to every one of them. A neighbour at a crease (a cap,
    // a wall the patch merely ends at) rules it out, so a freeform face
    // between such faces is left as it was.
    if (rOptions.mFillets) {
        const double crease = std::cos(20.0 * kFcPi / 180.0);
        for (std::size_t r = 0; r < m; ++r) {
            bool finite = true;
            for (double x : skin.mPatches[r].mParams)
                finite = finite && std::isfinite(x);
            const std::int32_t type = skin.mPatches[r].mType;
            const bool freeform =
                analytic[r] && type == static_cast<std::int32_t>(kFitCadBsplineType);
            const bool unexplained = !analytic[r] && !(type >= 1 && type <= 5 && finite);
            // A network patch is a face of its own network, not a blend.
            if (facets_of[r].empty() || (!freeform && !unexplained) ||
                (r < network_row.size() && network_row[r]))
                continue;
            std::vector<gp_Pnt> pts;
            std::vector<std::array<std::size_t, 3>> tris;
            std::vector<detail::PatchVec> flat;
            patch_mesh(r, pts, tris, flat);
            const std::vector<std::int64_t> global = patch_globals(r);
            if (global.size() != pts.size())
                continue;
            bool creased = false;
            auto across = [&](std::size_t A, std::size_t B) -> std::int64_t {
                std::int64_t mine = -1, theirs = -1;
                const std::int64_t row = row_across(r, global[A], global[B], &mine, &theirs);
                if (mine >= 0 && theirs >= 0) {
                    const std::size_t a = static_cast<std::size_t>(mine) * 3;
                    const std::size_t b = static_cast<std::size_t>(theirs) * 3;
                    const double c = fc.mNormals[a] * fc.mNormals[b] +
                                     fc.mNormals[a + 1] * fc.mNormals[b + 1] +
                                     fc.mNormals[a + 2] * fc.mNormals[b + 2];
                    creased = creased || c < crease;
                }
                return row;
            };
            const detail::VertexBlendBoundary vb = detail::vertex_blend_boundary(tris, across);
            if (!vb.mOk || creased)
                continue;
            bool fitted_all = true;
            for (std::int64_t q : vb.mNeighbours)
                fitted_all = fitted_all && analytic[static_cast<std::size_t>(q)] &&
                             !surface[static_cast<std::size_t>(q)].IsNull();
            if (!fitted_all)
                continue;
            try {
                std::vector<std::uint8_t> on_loop(pts.size(), 0);
                for (std::size_t v : vb.mLoop)
                    on_loop[v] = 1;
                FcHandle<Geom_Surface> vs =
                    fc_vertex_blend_surface(pts, vb.mLoop, vb.mAcross, on_loop, surface, tol);
                if (vs.IsNull() || !fc_surface_fits(vs, pts, tris, tol, true, true)) {
                    log::debug("fit_cad: patch {}: the vertex blend {}", skin.mPatches[r].mLabel,
                               vs.IsNull() ? "has no plate" : "strays from it");
                    continue;
                }
                // Tangent to every neighbour along the loop, within 3 degrees.
                bool tangent = true;
                const double tol_angle = std::cos(3.0 * kFcPi / 180.0);
                for (std::size_t i = 0; i < vb.mLoop.size() && tangent; ++i) {
                    const gp_Vec a = fc_normal_at(vs, pts[vb.mLoop[i]]);
                    const gp_Vec b = fc_normal_at(surface[static_cast<std::size_t>(vb.mAcross[i])],
                                                  pts[vb.mLoop[i]]);
                    if (a.Magnitude() <= 0.0 || b.Magnitude() <= 0.0)
                        continue;
                    tangent = std::fabs(a.Dot(b)) >= tol_angle * a.Magnitude() * b.Magnitude();
                }
                if (!tangent) {
                    log::debug(
                        "fit_cad: patch {}: the vertex blend is not tangent to its "
                        "neighbours",
                        skin.mPatches[r].mLabel);
                    continue;
                }
                accept(r, vs, kFitCadVertexBlendType, 0);
                log::debug("fit_cad: patch {} is a vertex blend among {} faces",
                           skin.mPatches[r].mLabel, vb.mNeighbours.size());
            } catch (const Standard_Failure& rErr) {
                log::debug("fit_cad: no vertex blend for patch {}: {}", skin.mPatches[r].mLabel,
                           rErr.GetMessageString());
            }
        }
    }

    // How far each face may lie from the skin: the tolerance, or further
    // where the skin's own facets lie further from the fitted surface (a
    // coarse tessellation's chords sag away from its cylinder by their
    // sagitta, measured at each facet's centre and edge midpoints).
    std::vector<double> allowance(m, tol);
    for (std::size_t r = 0; r < m; ++r) {
        const FcPatch& p = skin.mPatches[r];
        if (!p.mAccepted || p.mType < 1 || p.mType > 5 ||
            static_cast<std::size_t>(r) >= facets_of.size())
            continue;
        double sag = 0.0;
        for (std::size_t f : facets_of[r]) {
            const std::int64_t b = fc.mStart[f], n = fc.mStart[f + 1] - b;
            double c[3] = {0.0, 0.0, 0.0};
            for (std::int64_t k = 0; k < n; ++k) {
                const double* a =
                    &fc.mXyz[static_cast<std::size_t>(fc.mNodes[static_cast<std::size_t>(b + k)]) *
                             3];
                const double* e =
                    &fc.mXyz[static_cast<std::size_t>(
                                 fc.mNodes[static_cast<std::size_t>(b + (k + 1) % n)]) *
                             3];
                const double mid[3] = {0.5 * (a[0] + e[0]), 0.5 * (a[1] + e[1]),
                                       0.5 * (a[2] + e[2])};
                sag =
                    std::max(sag, detail::primitive_point_distance(p.mType, p.mParams.data(), mid));
                for (int j = 0; j < 3; ++j)
                    c[j] += a[j] / static_cast<double>(n);
            }
            sag = std::max(sag, detail::primitive_point_distance(p.mType, p.mParams.data(), c));
        }
        if (std::isfinite(sag))
            allowance[r] = std::max(tol, 1.05 * sag);
    }
    // An adaptive freeform face (outside the legacy mode) may lie as far from
    // its patch as the patch's own noise floor allows.
    for (std::size_t r = 0; r < m && r < row_info.size(); ++r)
        if (row_info[r].mAllowance > tol)
            allowance[r] = row_info[r].mAllowance;

    // Freeform surfaces run on past their patches, so an edge with a
    // neighbour just outside the fitted area is still an intersection; the
    // extensions are intersected only where the surfaces themselves find no
    // branch along a chain, and each face keeps its own surface.
    std::vector<FcHandle<Geom_Surface>> extended;
    if (rOptions.mExtendSurfaces)
        for (std::size_t r = 0; r < m; ++r) {
            if (skin.mPatches[r].mType != kFitCadBsplineType || surface[r].IsNull() ||
                facets_of[r].empty())
                continue;
            const FcHandle<Geom_BSplineSurface> bs =
                FcHandle<Geom_BSplineSurface>::DownCast(surface[r]);
            if (bs.IsNull())
                continue;
            double blo[3] = {std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity()};
            double bhi[3] = {-blo[0], -blo[1], -blo[2]};
            for (std::size_t f : facets_of[r])
                for (std::int64_t k = fc.mStart[f]; k < fc.mStart[f + 1]; ++k) {
                    const std::size_t a =
                        static_cast<std::size_t>(fc.mNodes[static_cast<std::size_t>(k)]);
                    for (int j = 0; j < 3; ++j) {
                        blo[j] = std::min(blo[j], fc.mXyz[a * 3 + static_cast<std::size_t>(j)]);
                        bhi[j] = std::max(bhi[j], fc.mXyz[a * 3 + static_cast<std::size_t>(j)]);
                    }
                }
            const double len = 0.1 * std::sqrt((bhi[0] - blo[0]) * (bhi[0] - blo[0]) +
                                               (bhi[1] - blo[1]) * (bhi[1] - blo[1]) +
                                               (bhi[2] - blo[2]) * (bhi[2] - blo[2]));
            if (!(len > 0.0))
                continue;
            try {
                FcHandle<Geom_BoundedSurface> ext =
                    FcHandle<Geom_BoundedSurface>::DownCast(bs->Copy());
                for (const bool in_u : {true, false}) {
                    if (in_u ? (bs->IsUPeriodic() || bs->IsUClosed())
                             : (bs->IsVPeriodic() || bs->IsVClosed()))
                        continue;
                    for (const bool after : {true, false})
                        GeomLib::ExtendSurfByLength(ext, len, 1, in_u, after);
                }
                if (!ext.IsNull()) {
                    extended.resize(m);
                    extended[r] = ext;
                }
            } catch (const Standard_Failure& rErr) {
                log::debug("fit_cad: patch {}'s B-spline surface could not be extended: {}",
                           skin.mPatches[r].mLabel, rErr.GetMessageString());
            }
        }

    // An adaptive face the checks demote, halved across its longest extent and
    // each half fitted in turn (v17.12.0): the tiles of the general faces,
    // reached from the checks rather than the fit.
    // Build, check, demote and rebuild.
    int escalation_rounds = 0;
    FcBuild build;
    FcTriangulation tri;
    std::vector<std::uint8_t> no_merge(nf, 0);
    for (std::int64_t round = 0;; ++round) {
        bool any = false;
        for (std::uint8_t a : analytic)
            any = any || a;
        try {
            build = fc_build(skin, surface, analytic, tol, rOptions.mExactEdges,
                             rOptions.mMergeFaces, no_merge, rOptions.mSewTolerance,
                             rOptions.mSnapJunctions, extended.empty() ? nullptr : &extended,
                             network_curves.empty() ? nullptr : &network_curves, &network_row,
                             rOptions.mFreeformMode != FitCadFreeformMode::Legacy);
        } catch (const Standard_Failure& rErr) {
            build = FcBuild{};
            build.mBadRow.assign(m, 0);
            for (std::size_t r = 0; r < m; ++r)
                build.mBadRow[r] =
                    analytic[r] ? static_cast<std::uint8_t>(FitCadDemotionReason::Other) : 0;
            log::debug("fit_cad: OpenCASCADE failed building the solid: {}",
                       rErr.GetMessageString());
        }
        bool bad = false;
        for (std::uint8_t b : build.mBadRow)
            bad = bad || b;
        if (!build.mUnmerge.empty()) {
            // Rebuild with those coplanar groups split into their facets;
            // not a demotion round, and it ends as each pass splits a group.
            bool split = false;
            for (std::size_t f : build.mUnmerge) {
                split = split || !no_merge[f];
                no_merge[f] = 1;
            }
            if (split) {
                --round;
                continue;
            }
        }
        if (!bad && !build.mShape.IsNull()) {
            BRepMesh_IncrementalMesh mesher(build.mShape, deflection, Standard_False, 0.5,
                                            Standard_False);
            tri = fc_triangulate(build.mFaces);
            if (any) {
                // Faces against the skin, and the skin against the solid.
                if (!tri.mComplete) {
                    for (std::size_t i = 0; i < build.mFaces.size(); ++i) {
                        TopLoc_Location loc;
                        const std::int64_t g = build.mFaceGroup[i];
                        if (BRep_Tool::Triangulation(build.mFaces[i], loc).IsNull() && g >= 0 &&
                            g < static_cast<std::int64_t>(m)) {
                            log::debug("fit_cad: demoting patch {}: its face does not mesh",
                                       skin.mPatches[static_cast<std::size_t>(g)].mLabel);
                            build.mBadRow[static_cast<std::size_t>(g)] =
                                static_cast<std::uint8_t>(FitCadDemotionReason::KernelMeshCoarse);
                        }
                    }
                }
                const std::vector<double> to_skin = fc_unsigned_distances(skin_mesh, tri.mXyz);
                for (std::size_t t = 0; t < tri.mFace.size(); ++t) {
                    const std::int64_t g = build.mFaceGroup[static_cast<std::size_t>(tri.mFace[t])];
                    if (g < 0 || g >= static_cast<std::int64_t>(m))
                        continue;
                    for (int k = 0; k < 3; ++k) {
                        const double d = to_skin[static_cast<std::size_t>(
                            tri.mTri[t * 3 + static_cast<std::size_t>(k)])];
                        if (d > allowance[static_cast<std::size_t>(g)] &&
                            !build.mBadRow[static_cast<std::size_t>(g)]) {
                            const std::size_t v = tri.mTri[t * 3 + static_cast<std::size_t>(k)];
                            log::debug(
                                "fit_cad: demoting patch {}: its face lies {} from the "
                                "skin, at ({}, {}, {})",
                                skin.mPatches[static_cast<std::size_t>(g)].mLabel, d,
                                tri.mXyz[v * 3], tri.mXyz[v * 3 + 1], tri.mXyz[v * 3 + 2]);
                            build.mBadRow[static_cast<std::size_t>(g)] =
                                static_cast<std::uint8_t>(FitCadDemotionReason::HausdorffExceeded);
                        }
                    }
                }
                if (!tri.mTri.empty()) {
                    const Mesh solid = fc_triangle_mesh(tri.mXyz, tri.mTri);
                    const std::vector<double> to_solid = fc_unsigned_distances(solid, fc.mXyz);
                    for (std::size_t f = 0; f < nf; ++f) {
                        const std::int64_t r = skin.mRow[f];
                        if (r < 0 || !analytic[static_cast<std::size_t>(r)])
                            continue;
                        for (std::int64_t k = fc.mStart[f]; k < fc.mStart[f + 1]; ++k) {
                            const double d = to_solid[static_cast<std::size_t>(
                                fc.mNodes[static_cast<std::size_t>(k)])];
                            const double limit =
                                static_cast<std::size_t>(r) < row_info.size() &&
                                        row_info[static_cast<std::size_t>(r)].mAllowance > tol
                                    ? allowance[static_cast<std::size_t>(r)]
                                    : tol;
                            if (d > limit && !build.mBadRow[static_cast<std::size_t>(r)]) {
                                const gp_Pnt at =
                                    fc_point(fc.mXyz, fc.mNodes[static_cast<std::size_t>(k)]);
                                log::debug(
                                    "fit_cad: demoting patch {}: the skin lies {} from it, "
                                    "at ({}, {}, {})",
                                    skin.mPatches[static_cast<std::size_t>(r)].mLabel, d, at.X(),
                                    at.Y(), at.Z());
                                build.mBadRow[static_cast<std::size_t>(r)] =
                                    static_cast<std::uint8_t>(
                                        FitCadDemotionReason::HausdorffExceeded);
                            }
                        }
                    }
                }
                for (std::uint8_t b : build.mBadRow)
                    bad = bad || b;
            }
        }
        if (!bad || !any)
            break;
        // Each round demotes at least one patch, so the loop ends. Past
        // `max_demotion_rounds` a failing patch takes its analytic neighbours
        // with it, which converges faster than one at a time and still keeps
        // the faces far from any failure.
        const bool widen = round >= rOptions.mMaxDemotionRounds;
        std::vector<std::uint8_t> drop(build.mBadRow.begin(), build.mBadRow.end());
        if (widen)
            for (const auto& kv : skin.mEdgeFacets) {
                const std::int64_t ra = skin.mRow[static_cast<std::size_t>(kv.second[0])];
                const std::int64_t rb = skin.mRow[static_cast<std::size_t>(kv.second[1])];
                if (ra < 0 || rb < 0 || ra == rb)
                    continue;
                // A neighbour taken along keeps its own reason, if it has one.
                constexpr auto kOther = static_cast<std::uint8_t>(FitCadDemotionReason::Other);
                if (build.mBadRow[static_cast<std::size_t>(ra)] &&
                    !drop[static_cast<std::size_t>(rb)])
                    drop[static_cast<std::size_t>(rb)] = kOther;
                if (build.mBadRow[static_cast<std::size_t>(rb)] &&
                    !drop[static_cast<std::size_t>(ra)])
                    drop[static_cast<std::size_t>(ra)] = kOther;
            }
        // In the freeform modes, a failing round that involves one of their
        // faces demotes only those first: an analytic neighbour failing with
        // it (on an edge or a sliver they share) usually stands once the
        // freeform face is facets, as it would have in the legacy mode.
        if (rOptions.mFreeformMode != FitCadFreeformMode::Legacy && !widen) {
            const auto freeform = [&](std::size_t R) {
                return R < feature.size() && (feature[R] & (kFcAdaptive | kFcGeneral)) != 0;
            };
            bool any_freeform = false;
            for (std::size_t r = 0; r < m; ++r)
                any_freeform = any_freeform || (drop[r] && analytic[r] && freeform(r));
            if (any_freeform)
                for (std::size_t r = 0; r < m; ++r)
                    if (!freeform(r))
                        drop[r] = 0;
        }
        std::vector<std::size_t> escalate;
        for (std::size_t r = 0; r < m; ++r)
            if (drop[r] && analytic[r]) {
                analytic[r] = 0;
                fc_row(row_info, r).mReason = drop[r];
                if ((feature[r] & kFcAdaptive) != 0 && !network_row[r])
                    escalate.push_back(r);
            }
        // Outside the legacy mode (v17.12.0), an adaptive face the checks
        // demote climbs one more rung before it is facets: under `auto` the
        // patch network of its region, then (and in the other modes) its
        // halves, each fitted in turn, down to the tile depth. Such a round
        // is not a demotion round: each rung is climbed once per row.
        if (rOptions.mFreeform && rOptions.mFreeformMode != FitCadFreeformMode::Legacy &&
            !escalate.empty() && escalation_rounds < rOptions.mMaxDemotionRounds) {
            ++escalation_rounds;
            const std::size_t before = skin.mPatches.size();
            std::vector<std::size_t> halve;
            std::vector<std::size_t> fresh;
            if (rOptions.mFreeformMode == FitCadFreeformMode::Auto)
                for (std::size_t r : escalate)
                    if (!network_tried[r]) {
                        network_tried[r] = 1;
                        fresh.push_back(r);
                    }
            std::vector<std::size_t> missed;
            if (!fresh.empty())
                missed = occt_network(fresh, true);
            for (std::size_t r : escalate)
                if (!facets_of[r].empty())
                    halve.push_back(r);
            const std::size_t networked = skin.mPatches.size();
            for (std::size_t r : halve)
                halve_and_refit(r);
            for (std::size_t r : missed)
                refit_or_halve(r);
            m = skin.mPatches.size();
            allowance.resize(m, tol);
            for (std::size_t r = before; r < m; ++r)
                if (r < row_info.size() && row_info[r].mAllowance > tol)
                    allowance[r] = row_info[r].mAllowance;
            feature.resize(m, 0);
            network_row.resize(m, 0);
            network_tried.resize(m, 1);
            if (m > before) {
                log::debug(
                    "fit_cad: {} demoted adaptive face(s) rebuilt as {} network row(s) and {} "
                    "half(s)",
                    escalate.size(), networked - before, m - networked);
                --round;
            }
        }
    }
    if (build.mShape.IsNull())
        throw std::runtime_error(std::string(kFcPrefix) +
                                 "OpenCASCADE could not build even the faceted solid");

    // The report.
    result.mValid = build.mValid;
    result.mNumPatches = static_cast<std::int64_t>(num_input_patches);
    result.mNumExactEdges = build.mNumExact;
    result.mNumApproxEdges = build.mNumApprox;
    result.mNumSnappedJunctions = build.mNumSnapped;
    result.mNumExtendedEdges = build.mNumExtended;
    result.mClosed = fc_closed(build.mShape);
    for (std::size_t r = 0; r < m; ++r)
        if (skin.mPatches[r].mAccepted && !facets_of[r].empty() && !analytic[r])
            ++result.mNumDemoted;
    for (TopExp_Explorer ex(build.mShape, TopAbs_SOLID); ex.More(); ex.Next())
        ++result.mNumSolids;
    GProp_GProps props;
    BRepGProp::VolumeProperties(build.mShape, props);
    result.mVolume = props.Mass();
    for (double v : piece_volume)
        result.mSkinVolume += std::fabs(v);
    result.mVolumeRatio = result.mSkinVolume > 0.0 ? result.mVolume / result.mSkinVolume
                                                   : std::numeric_limits<double>::quiet_NaN();
    result.mNumFaces = static_cast<std::int64_t>(build.mFaces.size());
    std::vector<std::int64_t> face_type(build.mFaces.size(), 0),
        face_patch(build.mFaces.size(), -1);
    const std::int64_t mi = static_cast<std::int64_t>(m);
    for (std::size_t i = 0; i < build.mFaces.size(); ++i) {
        const std::int64_t g = build.mFaceGroup[i];
        if (g >= 0 && g < mi) {
            const FcPatch& p = skin.mPatches[static_cast<std::size_t>(g)];
            face_type[i] = p.mType;
            face_patch[i] = p.mLabel;
        } else if (g >= mi) {
            const std::int64_t r = skin.mRow[static_cast<std::size_t>(g - mi)];
            face_patch[i] = r >= 0 ? skin.mPatches[static_cast<std::size_t>(r)].mLabel : -1;
        }
        const std::uint8_t feat =
            g >= 0 && g < mi ? feature[static_cast<std::size_t>(g)] : std::uint8_t(0);
        switch (face_type[i]) {
            case 1:
                ++result.mNumPlane;
                break;
            case 2:
                ++result.mNumSphere;
                break;
            case 3:
                ++result.mNumCylinder;
                break;
            case 4:
                ++result.mNumCone;
                break;
            case 5:
                ++result.mNumTorus;
                break;
            case kFitCadBsplineType:
                ++result.mNumBspline;
                if (g >= 0 && g < mi && network_row[static_cast<std::size_t>(g)])
                    ++result.mNumNetworkPatches;
                else if (feat & kFcAdaptive)
                    ++result.mNumFreeformAdaptive;
                if (feat & kFcGeneral) {
                    const std::size_t gr = static_cast<std::size_t>(g);
                    ++result.mNumFreeformGeneral;
                    result.mNumFreeformRefinements +=
                        gr < general_refinements.size() ? general_refinements[gr] : 0;
                    result.mAvgFreeformPoles +=
                        static_cast<double>(gr < general_poles.size() ? general_poles[gr] : 0);
                }
                break;
            case kFitCadFilletType:
                ++result.mNumFillet;
                result.mNumVariableFillet += (feat & kFcVariable) ? 1 : 0;
                result.mNumBlend += (feat & kFcBlend) ? 1 : 0;
                result.mNumConicBlend += (feat & kFcConic) ? 1 : 0;
                break;
            case kFitCadHelicalType:
                ++result.mNumHelical;
                result.mNumRoundedThread += (feat & kFcRounded) ? 1 : 0;
                result.mNumTaperedThread += (feat & kFcTapered) ? 1 : 0;
                result.mNumVariablePitchThread += (feat & kFcVariablePitch) ? 1 : 0;
                break;
            case kFitCadVertexBlendType:
                ++result.mNumVertexBlend;
                break;
            case kFitCadChamferType:
                ++result.mNumChamfer;
                break;
            default:
                ++result.mNumFacetedFaces;
                break;
        }
    }
    if (result.mNumFreeformGeneral > 0)
        result.mAvgFreeformPoles /= static_cast<double>(result.mNumFreeformGeneral);

    // The re-meshed solid.
    result.mMesh = fc_triangle_mesh(tri.mXyz, tri.mTri);
    const std::size_t nt = tri.mFace.size();
    NDArray cface = NDArray::Uninit(DType::Int64, {nt});
    NDArray ctype = NDArray::Uninit(DType::Int64, {nt});
    NDArray cpatch = NDArray::Uninit(DType::Int64, {nt});
    for (std::size_t t = 0; t < nt; ++t) {
        const std::size_t i = static_cast<std::size_t>(tri.mFace[t]);
        cface.As<std::int64_t>()[t] = static_cast<std::int64_t>(i);
        ctype.As<std::int64_t>()[t] = face_type[i];
        cpatch.As<std::int64_t>()[t] = face_patch[i];
    }
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(cface));
    result.mMesh.AddCellData(kFitCadFaceName, std::move(blocks));
    blocks.clear();
    blocks.push_back(std::move(ctype));
    result.mMesh.AddCellData(kFitCadTypeName, std::move(blocks));
    blocks.clear();
    blocks.push_back(std::move(cpatch));
    result.mMesh.AddCellData(kFitCadPatchName, std::move(blocks));
    if (nt > 0)
        result.mMaxDeviation = hausdorff_distance(skin_mesh, result.mMesh).mDistance;
    else
        result.mMaxDeviation = std::numeric_limits<double>::infinity();
    if (rOptions.mFaceDeviation && nt > 0)
        fc_face_deviation(fc.mXyz, skin_tri, tri, build.mFaces.size(), result);

    // The freeform-mode report (v17.12.0): the facet area share of the
    // re-mesh, and a record per patch tried as a surface.
    result.mBrepChecked = true;
    if (result.mNumNetworkPatches > 0) {
        result.mNumNetworkEdges = num_network_edges;
        result.mNumExtraordinaryVertices = num_network_ev;
        result.mMaxG1AngleDeg = network_g1;
    }
    result.mFreeformMode = rOptions.mFreeformMode;
    {
        std::vector<std::int64_t> tri_type(nt);
        for (std::size_t t = 0; t < nt; ++t)
            tri_type[t] = face_type[static_cast<std::size_t>(tri.mFace[t])];
        result.mFacetAreaFraction = fcn_facet_area_fraction(tri.mXyz, tri.mTri, tri_type);
        std::vector<std::int64_t> first_face(m, -1);
        for (std::size_t i = build.mFaces.size(); i-- > 0;) {
            const std::int64_t g = build.mFaceGroup[i];
            if (g >= 0 && g < mi)
                first_face[static_cast<std::size_t>(g)] = static_cast<std::int64_t>(i);
        }
        double sum_sq = 0.0;
        std::size_t num_freeform = 0;
        for (std::size_t r = 0; r < m; ++r) {
            const FcRowInfo& info = fc_row(row_info, r);
            if (facets_of[r].empty() || (surface[r].IsNull() && !info.mTried))
                continue;
            FitCadFaceRecord rec;
            rec.mPatch = skin.mPatches[r].mLabel;
            result.mNumEscalations += info.mEscalations;
            rec.mParamMethod = info.mParam;
            rec.mNumControl = info.mNumControl;
            rec.mIterations = info.mIterations;
            rec.mRefinements = info.mRefinements;
            rec.mEscalationSteps = info.mEscalations;
            rec.mMaxDist = info.mMaxDist;
            rec.mRmsDist = info.mRmsDist;
            rec.mReverseMaxDist = info.mReverseMaxDist;
            rec.mMaxG1AngleDeg = info.mMaxG1AngleDeg;
            if (analytic[r] && first_face[r] >= 0) {
                rec.mFace = first_face[r];
                rec.mKind = skin.mPatches[r].mType;
                if (rec.mKind == kFitCadBsplineType && std::isfinite(info.mMaxDist)) {
                    result.mMaxFreeformDistance =
                        std::max(result.mMaxFreeformDistance, info.mMaxDist);
                    sum_sq += std::isfinite(info.mRmsDist) ? info.mRmsDist * info.mRmsDist : 0.0;
                    ++num_freeform;
                }
            } else {
                rec.mDemotionReason = static_cast<FitCadDemotionReason>(
                    info.mReason ? info.mReason
                                 : static_cast<std::uint8_t>(FitCadDemotionReason::Other));
            }
            result.mFaces.push_back(rec);
        }
        if (num_freeform > 0)
            result.mRmsFreeformDistance = std::sqrt(sum_sq / static_cast<double>(num_freeform));
    }
    if (rOptions.mFreeformMode != FitCadFreeformMode::Legacy)
        result.mMesh.AddFieldData(kFitCadFaceTableName, fcn_face_table(result.mFaces));

    if (!result.mValid)
        log::warn("fit_cad: the solid does not pass OpenCASCADE's B-rep check");
    if (result.mNumDemoted > 0)
        log::info("fit_cad: {} fitted patch(es) demoted to facets", result.mNumDemoted);
    // The solid as NURBS, its faces in the order fit_cad:face numbers them.
    try {
        result.mBrep = detail::occt_nurbs(build.mShape).mInfo;
    } catch (const Standard_Failure& rErr) {
        log::warn("fit_cad: the solid could not be converted to NURBS: {}",
                  rErr.GetMessageString());
    }
    if (!rStepPath.empty())
        fc_write_step(build.mShape, rStepPath);
    if (pShape)
        *pShape = build.mShape;
    fcn_report_refined(refined, work, result);
    return result;
}

// The OpenCASCADE backend. With `mRefineSegmentation`, the solid is built from
// the refined patches and from the patches as given, and the refined one is
// kept unless the other is better -- valid where it is not, then fewer facet
// faces, fewer demotions, fewer faces: the refinement knows primitives only,
// and may cut into strips a blend or a corner the feature fits explain whole.
FitCadResult fc_occt_fit_cad(const Mesh& rMesh, const std::string& rStepPath,
                             const FitCadOptions& rOptions) {
    if (!rOptions.mRefineSegmentation)
        return fc_occt_fit_cad_once(rMesh, rStepPath, rOptions, nullptr);
    TopoDS_Shape refined_shape, plain_shape;
    FitCadResult refined = fc_occt_fit_cad_once(rMesh, "", rOptions, &refined_shape);
    FitCadOptions plain_options = rOptions;
    plain_options.mRefineSegmentation = false;
    FitCadResult plain = fc_occt_fit_cad_once(rMesh, "", plain_options, &plain_shape);
    const auto rank = [](const FitCadResult& rR) {
        return std::make_tuple(rR.mValid ? 0 : 1, rR.mNumFacetedFaces, rR.mNumDemoted,
                               rR.mNumFaces);
    };
    const bool keep = rank(refined) <= rank(plain);
    log::info(
        "fit_cad: the refined patches build {} face(s), {} faceted; the patches as given "
        "{}, {} faceted: {} kept",
        refined.mNumFaces, refined.mNumFacetedFaces, plain.mNumFaces, plain.mNumFacetedFaces,
        keep ? "the refinement" : "the patches as given");
    if (!keep) {
        // The refinement's report stays, its per-patch data does not: it
        // describes patches the solid was not built from.
        plain.mRefined = refined.mRefined;
        plain.mRefineIterations = refined.mRefineIterations;
        plain.mNumRefineSplits = refined.mNumRefineSplits;
        plain.mNumRefineMerges = refined.mNumRefineMerges;
        plain.mInitialExplainedAreaFraction = refined.mInitialExplainedAreaFraction;
        plain.mExplainedAreaFraction = refined.mExplainedAreaFraction;
    }
    FitCadResult& out = keep ? refined : plain;
    out.mRefinementKept = keep;
    if (!rStepPath.empty()) {
        detail::LibraryLock lock;
        FcQuietMessenger quiet;
        fc_write_step(keep ? refined_shape : plain_shape, rStepPath);
    }
    return std::move(out);
}

}  // namespace

#endif  // MESHIOPLUSPLUS_HAS_OCCT

FitCadResult fit_cad(const Mesh& rMesh, const std::string& rStepPath,
                     const FitCadOptions& rOptions) {
    if (std::isnan(rOptions.mTolerance) || std::isinf(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kFcnPrefix) + "tolerance must be finite");
    if ((rOptions.mTypes & ~kFitAllTypes) != 0)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "types may name only plane, sphere, cylinder, cone and torus");
    if (rOptions.mBackend != FitCadBackend::Auto && rOptions.mBackend != FitCadBackend::Occt &&
        rOptions.mBackend != FitCadBackend::Native)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "backend must be 'auto', 'occt' or 'native'");
    // The general freeform options (v17.11.0), checked whether or not they are used.
    if (!(rOptions.mFreeformTolerance >= 0.0) || std::isinf(rOptions.mFreeformTolerance))
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "freeform_tolerance must be finite and >= 0");
    if (!(rOptions.mFreeformSmoothing >= 0.0) || std::isinf(rOptions.mFreeformSmoothing))
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "freeform_smoothing must be finite and >= 0");
    if (rOptions.mMaxFreeformPoles < 4 || rOptions.mMaxFreeformPoles > 1024)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "max_freeform_poles must be 4 to 1024");
    if (rOptions.mMaxFreeformIterations < 1)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "max_freeform_iterations must be >= 1");
    if (rOptions.mMaxFreeformSplits < 0)
        throw std::invalid_argument(std::string(kFcnPrefix) + "max_freeform_splits must be >= 0");
    if (!(rOptions.mSewTolerance >= 0.0) || std::isinf(rOptions.mSewTolerance))
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "sew_tolerance must be finite and >= 0");
    // The freeform modes (v17.12.0), checked whether or not they are used.
    if (static_cast<std::uint32_t>(rOptions.mFreeformMode) > 3u)
        throw std::invalid_argument(
            std::string(kFcnPrefix) +
            "freeform_mode must be 'legacy', 'adaptive', 'network' or 'auto'");
    if (static_cast<std::uint32_t>(rOptions.mFreeformParameterization) > 3u)
        throw std::invalid_argument(
            std::string(kFcnPrefix) +
            "freeform_parameterization must be 'auto', 'projection', 'lscm' or 'harmonic'");
    if (static_cast<std::uint32_t>(rOptions.mExtraordinaryVertex) > 3u)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "extraordinary_vertex must be 'g0', 'g1', 'g1x' or 'g1xx'");
    if (static_cast<std::uint32_t>(rOptions.mLocalRefine) > 1u)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "local_refine must be 'knots' or 'thb'");
    if (!(rOptions.mFreeformG1ToleranceDeg > 0.0) || std::isinf(rOptions.mFreeformG1ToleranceDeg))
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "freeform_g1_tolerance_deg must be finite and > 0");
    if (rOptions.mFreeformMaxControlPoints < 16 || rOptions.mFreeformMaxControlPoints > 1048576)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "freeform_max_control_points must be 16 to 1048576");
    if (!(rOptions.mLayoutTargetPatchSize > 0.0) || std::isinf(rOptions.mLayoutTargetPatchSize))
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "layout_target_patch_size must be finite and > 0");
    if (rOptions.mLayoutMaxValenceDeviation < 1 || rOptions.mLayoutMaxValenceDeviation > 3)
        throw std::invalid_argument(std::string(kFcnPrefix) +
                                    "layout_max_valence_deviation must be 1 to 3");
    if (rOptions.mResegmentRounds < 0)
        throw std::invalid_argument(std::string(kFcnPrefix) + "resegment_rounds must be >= 0");
    // The freeform modes build on the general fit, and their faces run on
    // past their patches so an edge with an analytic neighbour is still an
    // intersection.
    if (rOptions.mFreeformMode != FitCadFreeformMode::Legacy &&
        (!rOptions.mFreeformGeneral || !rOptions.mExtendSurfaces)) {
        FitCadOptions general = rOptions;
        general.mFreeformGeneral = true;
        general.mExtendSurfaces = true;
        return fit_cad(rMesh, rStepPath, general);
    }
#ifdef MESHIOPLUSPLUS_HAS_OCCT
    if (rOptions.mBackend != FitCadBackend::Native)
        return fc_occt_fit_cad(rMesh, rStepPath, rOptions);
#else
    if (rOptions.mBackend == FitCadBackend::Occt)
        throw Unsupported(
            "fit_cad: this build has no OpenCASCADE backend (the PyPI wheels and the release "
            "binaries carry it); rebuild with -DMESHIOPLUSPLUS_WITH_OCCT=ON (OpenCASCADE 7.6 or "
            "later, build/build-occt.sh builds one), or use the native backend");
    // OpenCASCADE is what the result is meant to be built with: say once that
    // this build falls back, and what that costs. An atomic flag rather than
    // std::call_once, whose libstdc++ thread-locals the xeus-cpp JIT running
    // the C++ notebooks on the single header cannot resolve.
    if (rOptions.mBackend == FitCadBackend::Auto) {
        static std::atomic<bool> s_fcn_fallback_warned{false};
        if (!s_fcn_fallback_warned.exchange(true))
            log::warn(
                "fit_cad: this build has no OpenCASCADE, so the kernel-free native backend "
                "is used: no freeform, fillet or thread faces, no B-rep check and no solid "
                "volume. The PyPI wheels and the release binaries carry OpenCASCADE; a "
                "source build finds one with -DMESHIOPLUSPLUS_WITH_OCCT=ON (the default)");
    }
#endif
    return fcn_native_fit_cad(rMesh, rStepPath, rOptions);
}

}  // namespace meshioplusplus
