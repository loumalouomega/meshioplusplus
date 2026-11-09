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
// face_support and fit_cad_iterate (roadmap 7.7.1, 7.7.2). See
// operations/fit_cad_loop.hpp.

// System includes
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

// Project includes
#include "meshioplusplus/operations/fit_cad_loop.hpp"
#include "../detail/bounding_box_utils.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

[[noreturn]] void fcl_fail(const char* pOp, const std::string& rWhy) {
    throw std::invalid_argument(std::string("meshio++: ") + pOp + ": " + rWhy);
}

double fcl_diagonal(const Mesh& rMesh) {
    const NDArray& p = rMesh.Points();
    std::vector<Vec3> pts(rMesh.NumPoints());
    for (std::size_t i = 0; i < pts.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            pts[i][k] = detail::read_double(p, 3 * i + k);
    return detail::BoundingBoxUtils::Diagonal(detail::BoundingBoxUtils::FromPoints(pts));
}

// A per-cell Int64 array, block-major (absent: all -1).
std::vector<std::int64_t> fcl_cell_ints(const Mesh& rMesh, const std::string& rName) {
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    std::vector<std::int64_t> out(static_cast<std::size_t>(detail::total_cells(bases)), -1);
    if (!rMesh.HasCellData(rName))
        return out;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const NDArray& a = rMesh.CellData(rName, b);
        for (std::size_t c = 0; c < a.Size() && c < rMesh.Cells(b).NumCells(); ++c)
            out[static_cast<std::size_t>(bases[b]) + c] = detail::read_int(a, c);
    }
    return out;
}

// A per-cell Int64 array as cell data, block by block.
std::vector<NDArray> fcl_blocks(const Mesh& rMesh, const std::vector<std::int64_t>& rValues) {
    std::vector<NDArray> out;
    std::size_t at = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const std::size_t nc = rMesh.Cells(b).NumCells();
        NDArray a = NDArray::Uninit(DType::Int64, {nc});
        std::copy(rValues.begin() + static_cast<std::ptrdiff_t>(at),
                  rValues.begin() + static_cast<std::ptrdiff_t>(at + nc), a.As<std::int64_t>());
        out.push_back(std::move(a));
        at += nc;
    }
    return out;
}

detail::DistanceQuery fcl_query(const detail::TriangleSoup& rSoup) {
    SurfaceDistanceOptions o;
    o.mSign = SdfSign::Unsigned;
    o.mWatertightCheck = SdfWatertightCheck::Off;
    return detail::build_distance_query(rSoup, o);
}

// The share of a fit_cad re-mesh's area on analytic or feature faces
// (fit_cad:type other than 0, a facet).
double fcl_explained(const Mesh& rSolid) {
    const std::vector<std::int64_t> type = fcl_cell_ints(rSolid, kFitCadTypeName);
    const detail::TriangleSoup soup = detail::build_triangle_soup(rSolid, "");
    double all = 0.0, explained = 0.0;
    for (std::size_t t = 0; t < soup.NumTriangles(); ++t) {
        const Vec3& a = soup.mCorners[3 * t];
        const double area = 0.5 * detail::vec3_norm(detail::vec3_cross(
                                      detail::vec3_sub(soup.mCorners[3 * t + 1], a),
                                      detail::vec3_sub(soup.mCorners[3 * t + 2], a)));
        all += area;
        if (type[static_cast<std::size_t>(soup.mSourceCell[t])] != 0)
            explained += area;
    }
    return all > 0.0 ? explained / all : 0.0;
}

}  // namespace

FaceSupportResult face_support(const Mesh& rSolid, const Mesh& rScan,
                               const FaceSupportOptions& rOptions) {
    const FaceSupportOptions& o = rOptions;
    if (!rSolid.HasCellData(kFitCadFaceName))
        fcl_fail("face_support", "the solid has no fit_cad:face cell data (a fit_cad re-mesh)");
    if (!(std::isfinite(o.mTolerance) && std::isfinite(o.mMeasuredShare) &&
          std::isfinite(o.mInferredShare) && o.mInferredShare >= 0.0 &&
          o.mInferredShare <= o.mMeasuredShare && o.mMeasuredShare <= 1.0))
        fcl_fail("face_support", "the shares must satisfy 0 <= inferred <= measured <= 1");
    FaceSupportResult result;
    result.mTolerance = o.mTolerance > 0.0 ? o.mTolerance : 1e-3 * fcl_diagonal(rScan);
    const double tol = result.mTolerance;

    // The scan without its fills, as a distance query.
    Mesh scan = detail::clone_mesh(rScan);
    const std::vector<std::int64_t> hole = fcl_cell_ints(rScan, o.mHoleData);
    std::vector<std::int64_t> kept;
    for (std::size_t c = 0; c < hole.size(); ++c)
        if (hole[c] < 0)
            kept.push_back(static_cast<std::int64_t>(c));
    NDArray ids = NDArray::Uninit(DType::Int64, {kept.size()});
    std::copy(kept.begin(), kept.end(), ids.As<std::int64_t>());
    scan.AddRegion(Region("fit_cad:measured_cells", RegionKind::Cell, std::move(ids)));
    const detail::TriangleSoup scan_soup =
        detail::build_triangle_soup(scan, "fit_cad:measured_cells");
    if (scan_soup.NumTriangles() == 0)
        fcl_fail("face_support", "the scan has no surface cell outside its fills");
    const detail::DistanceQuery scan_query = fcl_query(scan_soup);

    // Each solid triangle's four sub-triangle centroids, a quarter of its area each.
    const detail::TriangleSoup solid_soup = detail::build_triangle_soup(rSolid, "");
    const std::vector<std::int64_t> face_of_cell = fcl_cell_ints(rSolid, kFitCadFaceName);
    std::int64_t nfaces = 0;
    for (const std::int64_t f : face_of_cell)
        nfaces = std::max(nfaces, f + 1);
    std::vector<Vec3> samples;
    std::vector<double> weight;
    std::vector<std::int64_t> sample_face;
    for (std::size_t t = 0; t < solid_soup.NumTriangles(); ++t) {
        const Vec3& a = solid_soup.mCorners[3 * t];
        const Vec3& b = solid_soup.mCorners[3 * t + 1];
        const Vec3& c = solid_soup.mCorners[3 * t + 2];
        const double area = 0.5 * detail::vec3_norm(detail::vec3_cross(detail::vec3_sub(b, a),
                                                                       detail::vec3_sub(c, a)));
        const std::int64_t f = face_of_cell[static_cast<std::size_t>(solid_soup.mSourceCell[t])];
        const Vec3 ab = detail::vec3_scale(detail::vec3_add(a, b), 0.5);
        const Vec3 bc = detail::vec3_scale(detail::vec3_add(b, c), 0.5);
        const Vec3 ca = detail::vec3_scale(detail::vec3_add(c, a), 0.5);
        const Vec3 corners[4][3] = {{a, ab, ca}, {ab, b, bc}, {ca, bc, c}, {ab, bc, ca}};
        for (const auto& r_q : corners) {
            samples.push_back(detail::vec3_scale(
                detail::vec3_add(detail::vec3_add(r_q[0], r_q[1]), r_q[2]), 1.0 / 3.0));
            weight.push_back(0.25 * area);
            sample_face.push_back(f);
        }
    }
    const std::vector<detail::ClosestPointHit> hits =
        detail::query_closest_points(scan_query, samples);
    std::vector<double> within(static_cast<std::size_t>(nfaces), 0.0);
    std::vector<double> total(static_cast<std::size_t>(nfaces), 0.0);
    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (sample_face[i] < 0)
            continue;
        const auto f = static_cast<std::size_t>(sample_face[i]);
        total[f] += weight[i];
        if (hits[i].mDistance <= tol)
            within[f] += weight[i];
    }
    result.mFaceShare.resize(static_cast<std::size_t>(nfaces), 0.0);
    result.mFaceClass.resize(static_cast<std::size_t>(nfaces), 2);
    for (std::size_t f = 0; f < within.size(); ++f) {
        const double share = total[f] > 0.0 ? within[f] / total[f] : 0.0;
        result.mFaceShare[f] = share;
        const std::int64_t cls = share >= o.mMeasuredShare ? 0 : share >= o.mInferredShare ? 1 : 2;
        result.mFaceClass[f] = cls;
        (cls == 0 ? result.mNumMeasured : cls == 1 ? result.mNumPartial : result.mNumInferred) += 1;
    }

    // The scan's points (on measured cells) against the solid.
    std::set<std::int64_t> scan_pts;
    for (const auto& r_v : scan_soup.mVertices)
        scan_pts.insert(r_v.begin(), r_v.end());
    std::vector<Vec3> qs;
    for (const std::int64_t v : scan_pts)
        qs.push_back(scan_soup.mPoints[static_cast<std::size_t>(v)]);
    const std::vector<detail::ClosestPointHit> back =
        detail::query_closest_points(fcl_query(solid_soup), qs);
    std::size_t ok = 0;
    for (const detail::ClosestPointHit& r_h : back)
        ok += r_h.mDistance <= tol ? 1 : 0;
    result.mShareWithin =
        qs.empty() ? 0.0 : static_cast<double>(ok) / static_cast<double>(qs.size());

    result.mMesh = detail::clone_mesh(rSolid);
    NDArray cls = NDArray::Uninit(DType::Int64, {result.mFaceClass.size()});
    std::copy(result.mFaceClass.begin(), result.mFaceClass.end(), cls.As<std::int64_t>());
    NDArray share = NDArray::Uninit(DType::Float64, {result.mFaceShare.size()});
    std::copy(result.mFaceShare.begin(), result.mFaceShare.end(), share.As<double>());
    result.mMesh.AddFieldData(kFitCadFaceSupportName, std::move(cls));
    result.mMesh.AddFieldData(kFitCadFaceMeasuredName, std::move(share));
    std::vector<std::int64_t> per_cell(face_of_cell.size(), 2);
    for (std::size_t c = 0; c < face_of_cell.size(); ++c)
        if (face_of_cell[c] >= 0)
            per_cell[c] = result.mFaceClass[static_cast<std::size_t>(face_of_cell[c])];
    result.mMesh.AddCellData(kFitCadSupportName, fcl_blocks(result.mMesh, per_cell));
    return result;
}

FitCadPartClass fit_cad_part_class_from_name(const std::string& rName) {
    if (rName == "any")
        return FitCadPartClass::Any;
    if (rName == "analytic")
        return FitCadPartClass::Analytic;
    if (rName == "cast")
        return FitCadPartClass::Cast;
    if (rName == "freeform")
        return FitCadPartClass::Freeform;
    fcl_fail("fit_cad_iterate",
             "part_class must be any, analytic, cast or freeform, not '" + rName + "'");
}

FitCadIterateResult fit_cad_iterate(const Mesh& rScan, const std::string& rStepPath,
                                    const FitCadIterateOptions& rOptions) {
    const FitCadIterateOptions& o = rOptions;
    const char* op = "fit_cad_iterate";
    if (o.mMaxRounds < 1 || o.mWorstFaces < 1 || o.mMinPatchCells < 1)
        fcl_fail(op, "max_rounds, worst_faces and min_patch_cells must be at least 1");
    if (!(o.mRefineFactor > 0.0 && o.mRefineFactor < 1.0))
        fcl_fail(op, "refine_factor must be in (0, 1)");
    for (const double v :
         {o.mSegmentTolerance, o.mFitTolerance, o.mTargetWithin, o.mTargetExplained})
        if (!std::isfinite(v) || v > 1.0)
            fcl_fail(op, "the tolerances and targets are finite and at most 1");
    if (o.mPartClass != FitCadPartClass::Any && o.mPartClass != FitCadPartClass::Analytic &&
        o.mPartClass != FitCadPartClass::Cast && o.mPartClass != FitCadPartClass::Freeform)
        fcl_fail(op, "part_class must be any, analytic, cast or freeform");
    const int cls = static_cast<int>(o.mPartClass);
    const double class_target[4] = {0.9, 0.9, 0.85, 0.75};
    const double class_tol[4] = {1e-3, 1e-3, 2e-3, 4e-3};
    FitCadIterateResult result;
    result.mTargetWithin = o.mTargetWithin > 0.0 ? o.mTargetWithin : class_target[cls];
    const double seg_tol0 = o.mSegmentTolerance > 0.0 ? o.mSegmentTolerance : class_tol[cls];
    const double fit_tol = o.mFitTolerance > 0.0 ? o.mFitTolerance : seg_tol0;

    FitCadOptions fc = o.mFitCad;
    fc.mFaceDeviation = true;
    FitPrimitivesOptions fp;
    fp.mTolerance = fit_tol;
    fp.mExtendedReport = false;

    // Round 0's segmentation.
    SegmentOptions so;
    so.mMethod = SegmentMethod::Primitive;
    so.mTolerance = seg_tol0;
    so.mMinPatchCells = o.mMinPatchCells;
    Mesh segmented = segment(rScan, so).mMesh;
    double seg_tol = seg_tol0;
    double best_share = -1.0, best_explained = -1.0;
    for (std::int64_t round = 0; round < o.mMaxRounds; ++round) {
        Mesh skin = fit_primitives(segmented, fp).mMesh;
        FitCadResult fit = fit_cad(skin, "", fc);
        FaceSupportOptions sup_o = o.mSupport;
        if (!(sup_o.mTolerance > 0.0))
            sup_o.mTolerance = fit.mTolerance;
        FaceSupportResult sup = face_support(fit.mMesh, rScan, sup_o);
        FitCadIterateRound rec;
        rec.mShareWithin = sup.mShareWithin;
        rec.mExplained = fcl_explained(fit.mMesh);
        rec.mHausdorffMax = std::isfinite(fit.mHausdorffMax) ? fit.mHausdorffMax : 0.0;
        rec.mSegmentTolerance = seg_tol;
        {
            std::set<std::int64_t> labels;
            for (const std::int64_t l : fcl_cell_ints(segmented, kSegmentPatchName))
                if (l >= 0)
                    labels.insert(l);
            rec.mNumPatches = static_cast<std::int64_t>(labels.size());
        }
        rec.mNumFaces = fit.mCounts.mNumFaces;
        const bool better = rec.mShareWithin > best_share ||
                            (rec.mShareWithin == best_share && rec.mExplained > best_explained);
        if (better) {
            best_share = rec.mShareWithin;
            best_explained = rec.mExplained;
            result.mBestRound = round;
            result.mFit = std::move(fit);
            result.mSupport = std::move(sup);
            result.mSkin = detail::clone_mesh(skin);
        }
        result.mRounds.push_back(rec);
        const bool met =
            rec.mShareWithin >= result.mTargetWithin && rec.mExplained >= o.mTargetExplained;
        if (met) {
            result.mStopReason = "target";
            break;
        }
        if (round > 0 && !better) {
            result.mStopReason = "no improvement";
            break;
        }
        if (round + 1 == o.mMaxRounds) {
            result.mStopReason = "rounds";
            break;
        }
        // The worst faces: their deviation beyond the tolerance first, then
        // the largest faceted faces.
        const FitCadResult& cur = better ? result.mFit : fit;
        const std::vector<std::int64_t> face = fcl_cell_ints(cur.mMesh, kFitCadFaceName);
        const std::vector<std::int64_t> type = fcl_cell_ints(cur.mMesh, kFitCadTypeName);
        const std::vector<std::int64_t> patch = fcl_cell_ints(cur.mMesh, kFitCadPatchName);
        std::map<std::int64_t, std::set<std::int64_t>> face_patches;
        std::map<std::int64_t, double> face_area;
        std::map<std::int64_t, std::int64_t> face_type;
        {
            const detail::TriangleSoup soup = detail::build_triangle_soup(cur.mMesh, "");
            for (std::size_t t = 0; t < soup.NumTriangles(); ++t) {
                const auto c = static_cast<std::size_t>(soup.mSourceCell[t]);
                if (face[c] < 0)
                    continue;
                const Vec3& a = soup.mCorners[3 * t];
                face_area[face[c]] += 0.5 * detail::vec3_norm(detail::vec3_cross(
                                                detail::vec3_sub(soup.mCorners[3 * t + 1], a),
                                                detail::vec3_sub(soup.mCorners[3 * t + 2], a)));
                face_type[face[c]] = type[c];
                if (patch[c] >= 0)
                    face_patches[face[c]].insert(patch[c]);
            }
        }
        std::vector<std::pair<double, std::int64_t>> ranked;  // (-key, face)
        if (cur.mMesh.HasFieldData(kFitCadFaceHausdorffName)) {
            const NDArray& h = cur.mMesh.FieldData(kFitCadFaceHausdorffName);
            const std::size_t nf = h.Ndim() == 2 ? h.Shape()[0] : 0;
            for (std::size_t f = 0; f < nf; ++f) {
                const double d =
                    std::max(detail::read_double(h, 2 * f), detail::read_double(h, 2 * f + 1));
                if (std::isfinite(d) && d > cur.mTolerance)
                    ranked.push_back({-(1e300 + d), static_cast<std::int64_t>(f)});
            }
        }
        for (const auto& [f, area] : face_area)
            if (face_type[f] == 0)
                ranked.push_back({-area, f});
        std::sort(ranked.begin(), ranked.end());
        std::set<std::int64_t> worst;
        std::set<std::int64_t> seen_faces;
        for (const auto& [key, f] : ranked) {
            if (static_cast<std::int64_t>(seen_faces.size()) >= o.mWorstFaces)
                break;
            if (!seen_faces.insert(f).second)
                continue;
            for (const std::int64_t p : face_patches[f])
                worst.insert(p);
        }
        if (worst.empty()) {
            result.mStopReason = "nothing to refine";
            break;
        }
        result.mRounds.back().mResegmented.assign(worst.begin(), worst.end());
        // Re-segment those patches' cells finer and splice the labels in.
        std::vector<std::int64_t> labels = fcl_cell_ints(segmented, kSegmentPatchName);
        std::int64_t next = 0;
        std::vector<std::int64_t> region;
        for (std::size_t c = 0; c < labels.size(); ++c) {
            next = std::max(next, labels[c] + 1);
            if (labels[c] >= 0 && worst.count(labels[c]))
                region.push_back(static_cast<std::int64_t>(c));
        }
        Mesh marked = detail::clone_mesh(segmented);
        NDArray rid = NDArray::Uninit(DType::Int64, {region.size()});
        std::copy(region.begin(), region.end(), rid.As<std::int64_t>());
        marked.AddRegion(Region("fit_cad:refine", RegionKind::Cell, std::move(rid)));
        seg_tol *= o.mRefineFactor;
        SegmentOptions fine = so;
        fine.mTolerance = seg_tol;
        fine.mRegion = "fit_cad:refine";
        fine.mPrefix = "fit_cad_refine";
        const std::vector<std::int64_t> sub =
            fcl_cell_ints(segment(marked, fine).mMesh, kSegmentPatchName);
        for (const std::int64_t c : region)
            if (sub[static_cast<std::size_t>(c)] >= 0)
                labels[static_cast<std::size_t>(c)] = next + sub[static_cast<std::size_t>(c)];
        segmented = detail::clone_mesh(rScan);
        segmented.AddCellData(kSegmentPatchName, fcl_blocks(segmented, labels));
    }
    result.mTargetMet = result.mRounds[static_cast<std::size_t>(result.mBestRound)].mShareWithin >=
                            result.mTargetWithin &&
                        result.mRounds[static_cast<std::size_t>(result.mBestRound)].mExplained >=
                            o.mTargetExplained;
    if (!rStepPath.empty())
        fit_cad(result.mSkin, rStepPath, fc);
    return result;
}

}  // namespace meshioplusplus
