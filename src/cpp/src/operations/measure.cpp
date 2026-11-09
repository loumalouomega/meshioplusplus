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
// Elements, dimensions and GD&T tolerances (roadmap 7.2); datum alignment
// (roadmap 7.3.4). See operations/measure.hpp and doc/measure.md.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/measure.hpp"
#include "../detail/cell_region_mask.hpp"
#include "../detail/json_value.hpp"
#include "../detail/lsq_fit.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

using MsV3 = std::array<double, 3>;
using detail::JsonValue;

constexpr const char* kMsPrefix = "meshio++: measure: ";
constexpr double kMsPi = 3.14159265358979323846;
constexpr double kMsParallel = 1e-9;  // |sin| below which directions are parallel

[[noreturn]] void ms_refuse(const std::string& rWhy) {
    throw std::invalid_argument(std::string(kMsPrefix) + rWhy);
}

double ms_dot(const MsV3& rA, const MsV3& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}
MsV3 ms_add(const MsV3& rA, const MsV3& rB, double S = 1.0) {
    return {rA[0] + S * rB[0], rA[1] + S * rB[1], rA[2] + S * rB[2]};
}
MsV3 ms_sub(const MsV3& rA, const MsV3& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}
MsV3 ms_scale(const MsV3& rA, double S) {
    return {S * rA[0], S * rA[1], S * rA[2]};
}
MsV3 ms_cross(const MsV3& rA, const MsV3& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}
double ms_norm(const MsV3& rA) {
    return std::sqrt(ms_dot(rA, rA));
}
MsV3 ms_unit(const MsV3& rA) {
    const double n = ms_norm(rA);
    return n > 0.0 ? ms_scale(rA, 1.0 / n) : MsV3{0.0, 0.0, 0.0};
}

const char* ms_kind_name(ElementKind Kind) {
    switch (Kind) {
        case ElementKind::Point:
            return "point";
        case ElementKind::Line:
            return "line";
        case ElementKind::Plane:
            return "plane";
        case ElementKind::Sphere:
            return "sphere";
        case ElementKind::Cylinder:
            return "cylinder";
        default:
            return "cone";
    }
}

ElementKind ms_kind_from(const std::string& rName) {
    for (int k = 1; k <= 6; ++k)
        if (rName == ms_kind_name(static_cast<ElementKind>(k)))
            return static_cast<ElementKind>(k);
    ms_refuse("unknown element kind '" + rName + "'");
}

bool ms_line_like(const MeasureElement& rE) {
    return rE.mKind == ElementKind::Line || rE.mKind == ElementKind::Cylinder ||
           rE.mKind == ElementKind::Cone;
}
bool ms_point_like(const MeasureElement& rE) {
    return rE.mKind == ElementKind::Point || rE.mKind == ElementKind::Sphere;
}

// --- spec reading -----------------------------------------------------------

const JsonValue* ms_find(const JsonValue& rObj, const char* pKey) {
    return rObj.IsObject() ? rObj.Find(pKey) : nullptr;
}

std::string ms_text(const JsonValue& rObj, const char* pKey, const std::string& rWhere,
                    const std::string& rFallback = "") {
    const JsonValue* v = ms_find(rObj, pKey);
    if (!v)
        return rFallback;
    if (!v->IsString())
        ms_refuse(rWhere + ": '" + pKey + "' must be a string");
    return v->mString;
}

double ms_number(const JsonValue& rObj, const char* pKey, const std::string& rWhere,
                 double Fallback) {
    const JsonValue* v = ms_find(rObj, pKey);
    if (!v)
        return Fallback;
    if (!v->IsNumber() || !std::isfinite(v->mNumber))
        ms_refuse(rWhere + ": '" + pKey + "' must be a finite number");
    return v->mNumber;
}

MsV3 ms_vec3(const JsonValue& rV, const std::string& rWhere) {
    if (!rV.IsArray() || rV.mItems.size() != 3)
        ms_refuse(rWhere + ": a point is [x, y, z]");
    MsV3 out{};
    for (std::size_t k = 0; k < 3; ++k) {
        if (!rV.mItems[k].IsNumber())
            ms_refuse(rWhere + ": a point is [x, y, z]");
        out[k] = rV.mItems[k].mNumber;
    }
    return out;
}

std::vector<std::string> ms_names(const JsonValue& rObj, const char* pKey,
                                  const std::string& rWhere) {
    const JsonValue* v = ms_find(rObj, pKey);
    std::vector<std::string> out;
    if (!v)
        return out;
    if (v->IsString()) {
        out.push_back(v->mString);
        return out;
    }
    if (!v->IsArray())
        ms_refuse(rWhere + ": '" + pKey + "' must name elements");
    for (const JsonValue& r_i : v->mItems) {
        if (!r_i.IsString())
            ms_refuse(rWhere + ": '" + pKey + "' must name elements");
        out.push_back(r_i.mString);
    }
    return out;
}

// --- geometry -----------------------------------------------------------------

double ms_point_line(const MsV3& rP, const MsV3& rA, const MsV3& rD) {
    return ms_norm(ms_cross(ms_sub(rP, rA), rD));
}

// The axis' ends: its points' extent along the axis (the point alone without).
std::pair<MsV3, MsV3> ms_axis_ends(const MeasureElement& rE) {
    if (rE.mPoints.empty())
        return {rE.mPoint, rE.mPoint};
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (std::size_t i = 0; i + 2 < rE.mPoints.size(); i += 3) {
        const double t = ms_dot(
            ms_sub({rE.mPoints[i], rE.mPoints[i + 1], rE.mPoints[i + 2]}, rE.mPoint), rE.mDir);
        lo = std::min(lo, t);
        hi = std::max(hi, t);
    }
    return {ms_add(rE.mPoint, rE.mDir, lo), ms_add(rE.mPoint, rE.mDir, hi)};
}

// The spread of an element's points along a unit direction.
double ms_spread(const MeasureElement& rE, const MsV3& rM) {
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (std::size_t i = 0; i + 2 < rE.mPoints.size(); i += 3) {
        const double t = ms_dot({rE.mPoints[i], rE.mPoints[i + 1], rE.mPoints[i + 2]}, rM);
        lo = std::min(lo, t);
        hi = std::max(hi, t);
    }
    return rE.mPoints.empty() ? 0.0 : hi - lo;
}

// A unit vector perpendicular to rD, from rV when it is not parallel.
MsV3 ms_perp(const MsV3& rD, const MsV3& rV) {
    const MsV3 w = ms_sub(rV, ms_scale(rD, ms_dot(rV, rD)));
    if (ms_norm(w) > kMsParallel)
        return ms_unit(w);
    MsV3 e1{}, e2{};
    detail::lsq_frame(rD.data(), e1.data(), e2.data());
    return e1;
}

// --- the run ------------------------------------------------------------------

struct MsRun {
    const Mesh& mMesh;
    FitCriterion mCriterion = FitCriterion::LeastSquares;
    std::vector<MeasureElement> mElements;
    std::map<std::string, std::size_t> mIndex;
    std::vector<std::pair<std::string, std::vector<std::int64_t>>> mRegions;
    std::vector<std::string> mWarnings;

    const MeasureElement& Get(const std::string& rName, const std::string& rWhere) const {
        const auto it = mIndex.find(rName);
        if (it == mIndex.end())
            ms_refuse(rWhere + ": no element named '" + rName + "' before it");
        return mElements[it->second];
    }

    MeasureElement Fit(const JsonValue& rSpec, const std::string& rName);
    MeasureElement Construct(const JsonValue& rSpec, const std::string& rName);
    MeasureCheck Dimension(const JsonValue& rSpec, std::size_t K);
    MeasureCheck Tolerance(const JsonValue& rSpec, std::size_t K);
};

MeasureElement MsRun::Fit(const JsonValue& rSpec, const std::string& rName) {
    const std::string where = "element '" + rName + "'";
    const std::string type = ms_text(rSpec, "fit", where);
    FitPrimitivesOptions o;
    o.mCriterion = mCriterion;
    if (type == "plane")
        o.mTypes = 1u << 1;
    else if (type == "sphere")
        o.mTypes = 1u << 2;
    else if (type == "cylinder")
        o.mTypes = 1u << 3;
    else if (type == "cone")
        o.mTypes = 1u << 4;
    else if (type == "auto")
        o.mTypes = (1u << 1) | (1u << 2) | (1u << 3) | (1u << 4);
    else
        ms_refuse(where + ": 'fit' must be plane, sphere, cylinder, cone or auto");
    o.mTolerance = ms_number(rSpec, "tolerance", where, 1e-3);
    const JsonValue* absolute = ms_find(rSpec, "absolute");
    o.mRelative = !(absolute && absolute->IsBool() && absolute->mBool);
    o.mMinInlierFraction = ms_number(rSpec, "min_inlier_fraction", where, 0.9);
    const std::string region = ms_text(rSpec, "region", where);
    const JsonValue* seeds = ms_find(rSpec, "seeds");
    const JsonValue* patch = ms_find(rSpec, "patch");
    const std::string out_region = std::string("measure:") + rName;
    if (!region.empty()) {
        o.mRegion = region;
    } else if (seeds) {
        if (!seeds->IsArray() || seeds->mItems.empty())
            ms_refuse(where + ": 'seeds' is a list of [x, y, z]");
        for (const JsonValue& r_s : seeds->mItems)
            o.mSeeds.push_back(ms_vec3(r_s, where));
        o.mOutputRegion = out_region;
        o.mFeatureAngleDeg = ms_number(rSpec, "feature_angle", where, 30.0);
        o.mRimAngleDeg = ms_number(rSpec, "rim_angle", where, 15.0);
    } else if (patch) {
        if (!patch->IsNumber())
            ms_refuse(where + ": 'patch' is a segment:patch label");
    } else {
        ms_refuse(where + ": a fitted element needs a 'region', 'seeds' or a 'patch'");
    }
    // The cells of the element.
    std::size_t ncells = 0;
    for (std::size_t b = 0; b < mMesh.NumCellBlocks(); ++b)
        ncells += mMesh.Cells(b).NumCells();
    std::vector<std::uint8_t> mask;
    Mesh fitted;
    if (patch) {
        const std::int64_t label = static_cast<std::int64_t>(patch->mNumber);
        if (!mMesh.HasCellData("segment:patch"))
            ms_refuse(where + ": 'patch' needs segment:patch cell data");
        mask.assign(ncells, 0);
        std::size_t base = 0;
        for (std::size_t b = 0; b < mMesh.NumCellBlocks(); ++b) {
            const NDArray& l = mMesh.CellData("segment:patch", b);
            for (std::size_t c = 0; c < mMesh.Cells(b).NumCells() && c < l.Size(); ++c)
                mask[base + c] = detail::read_int(l, c) == label ? 1 : 0;
            base += mMesh.Cells(b).NumCells();
        }
        // Fit the patch's cells as a region of their own.
        fitted = detail::clone_mesh(
            mMesh, [](DataLocation, const std::string&, std::string&) { return true; });
        std::vector<std::int64_t> ids;
        for (std::size_t c = 0; c < ncells; ++c)
            if (mask[c])
                ids.push_back(static_cast<std::int64_t>(c));
        if (ids.empty())
            ms_refuse(where + ": no cell carries patch " + std::to_string(label));
        NDArray e = NDArray::Uninit(DType::Int64, {ids.size()});
        std::copy(ids.begin(), ids.end(), e.As<std::int64_t>());
        fitted.AddRegion(Region(out_region, RegionKind::Cell, std::move(e)));
        o.mRegion = out_region;
    }
    FitPrimitivesResult r = fit_primitives(patch ? fitted : mMesh, o);
    for (const std::string& w : r.mWarnings)
        mWarnings.push_back(where + ": " + w);
    const Mesh& m = r.mMesh;
    if (!patch)
        mask = detail::cell_region_mask(m, region.empty() ? out_region : region, ncells, kMsPrefix);
    const std::int64_t t = detail::read_int(m.FieldData("fit:patch_type"), 0);
    if (t < 1 || t > 4)
        ms_refuse(where + ": no " + (type == "auto" ? std::string("primitive") : type) +
                  " explains its cells within the tolerance");
    MeasureElement e;
    e.mName = rName;
    e.mMeasured = true;
    e.mSource = "fit:" + type;
    double p[10];
    for (std::size_t q = 0; q < 10; ++q)
        p[q] = detail::read_double(m.FieldData("fit:patch_params"), q);
    e.mPoint = {p[0], p[1], p[2]};
    e.mDir = ms_unit({p[3], p[4], p[5]});
    e.mKind = t == 1   ? ElementKind::Plane
              : t == 2 ? ElementKind::Sphere
              : t == 3 ? ElementKind::Cylinder
                       : ElementKind::Cone;
    if (t == 2)
        e.mDir = {0.0, 0.0, 0.0};
    e.mRadius = (t == 2 || t == 3) ? p[6] : 0.0;
    e.mAngleDeg = t == 4 ? p[8] : 0.0;
    e.mFormError = m.HasFieldData(kFitPatchFormErrorName)
                       ? detail::read_double(m.FieldData(kFitPatchFormErrorName), 0)
                       : std::numeric_limits<double>::quiet_NaN();
    // The points of its cells.
    std::set<std::int64_t> verts;
    std::vector<std::int64_t> cells;
    std::size_t base = 0;
    for (std::size_t b = 0; b < mMesh.NumCellBlocks(); ++b) {
        const auto cb = mMesh.Cells(b);
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            if (!mask[base + c])
                continue;
            cells.push_back(static_cast<std::int64_t>(base + c));
            const std::size_t k = cb.IsRagged() ? cb.RowSize(c) : cb.NodesPerCell();
            for (std::size_t i = 0; i < k; ++i)
                verts.insert(cb.IsRagged()
                                 ? cb.Row(c)[i]
                                 : detail::read_int(cb.Conn(), c * cb.NodesPerCell() + i));
        }
        base += cb.NumCells();
    }
    const NDArray& pts = mMesh.Points();
    const std::size_t dim = mMesh.PointDim();
    for (const std::int64_t v : verts)
        for (std::size_t k = 0; k < 3; ++k)
            e.mPoints.push_back(
                k < dim ? detail::read_double(pts, static_cast<std::size_t>(v) * dim + k) : 0.0);
    mRegions.emplace_back(out_region, std::move(cells));
    return e;
}

MeasureElement MsRun::Construct(const JsonValue& rSpec, const std::string& rName) {
    const std::string where = "element '" + rName + "'";
    const std::string how = ms_text(rSpec, "construct", where);
    const std::vector<std::string> of = ms_names(rSpec, "of", where);
    MeasureElement e;
    e.mName = rName;
    e.mSource = how;
    e.mInputs = of;
    e.mFormError = std::numeric_limits<double>::quiet_NaN();
    auto need = [&](std::size_t N) {
        if (of.size() != N)
            ms_refuse(where + ": '" + how + "' takes " + std::to_string(N) + " element(s)");
    };
    if (how == "intersection") {
        need(2);
        const MeasureElement& a = Get(of[0], where);
        const MeasureElement& b = Get(of[1], where);
        if (a.mKind == ElementKind::Plane && b.mKind == ElementKind::Plane) {
            const MsV3 u = ms_cross(a.mDir, b.mDir);
            const double s = ms_norm(u);
            if (s < kMsParallel)
                ms_refuse(where + ": the planes are parallel");
            const double d1 = ms_dot(a.mDir, a.mPoint), d2 = ms_dot(b.mDir, b.mPoint);
            const MsV3 p = ms_scale(
                ms_add(ms_scale(ms_cross(b.mDir, u), d1), ms_scale(ms_cross(u, a.mDir), d2)),
                1.0 / (s * s));
            e.mKind = ElementKind::Line;
            e.mDir = ms_unit(u);
            // The point of the line nearest the planes' points' midpoint.
            const MsV3 mid = ms_scale(ms_add(a.mPoint, b.mPoint), 0.5);
            e.mPoint = ms_add(p, e.mDir, ms_dot(ms_sub(mid, p), e.mDir));
            return e;
        }
        const MeasureElement* line = ms_line_like(a) ? &a : (ms_line_like(b) ? &b : nullptr);
        const MeasureElement* plane =
            a.mKind == ElementKind::Plane ? &a : (b.mKind == ElementKind::Plane ? &b : nullptr);
        if (line && plane) {
            const double dn = ms_dot(line->mDir, plane->mDir);
            if (std::fabs(dn) < kMsParallel)
                ms_refuse(where + ": the line is parallel to the plane");
            const double t = ms_dot(ms_sub(plane->mPoint, line->mPoint), plane->mDir) / dn;
            e.mKind = ElementKind::Point;
            e.mPoint = ms_add(line->mPoint, line->mDir, t);
            return e;
        }
        if (ms_line_like(a) && ms_line_like(b)) {
            const MsV3 w0 = ms_sub(a.mPoint, b.mPoint);
            const double aa = 1.0, bb = ms_dot(a.mDir, b.mDir), cc = 1.0;
            const double d = ms_dot(a.mDir, w0), f = ms_dot(b.mDir, w0);
            const double den = aa * cc - bb * bb;
            if (den < kMsParallel)
                ms_refuse(where + ": the lines are parallel");
            const double s = (bb * f - cc * d) / den, t = (aa * f - bb * d) / den;
            e.mKind = ElementKind::Point;
            e.mPoint =
                ms_scale(ms_add(ms_add(a.mPoint, a.mDir, s), ms_add(b.mPoint, b.mDir, t)), 0.5);
            return e;
        }
        ms_refuse(where + ": 'intersection' takes two planes, a line and a plane, or two lines");
    }
    if (how == "axis") {
        need(1);
        const MeasureElement& a = Get(of[0], where);
        if (!ms_line_like(a))
            ms_refuse(where + ": 'axis' takes a cylinder, a cone or a line");
        e.mKind = ElementKind::Line;
        e.mPoint = a.mPoint;
        e.mDir = a.mDir;
        e.mPoints = a.mPoints;  // its extent, for the axis' ends
        return e;
    }
    if (how == "centre" || how == "apex") {
        need(1);
        const MeasureElement& a = Get(of[0], where);
        if (how == "apex" ? a.mKind != ElementKind::Cone
                          : !(ms_point_like(a) || a.mKind == ElementKind::Cylinder))
            ms_refuse(where + ": '" + how + "' takes " +
                      (how == "apex" ? "a cone" : "a sphere, a cylinder or a point"));
        e.mKind = ElementKind::Point;
        e.mPoint = a.mPoint;
        return e;
    }
    if (how == "mid_plane") {
        need(2);
        const MeasureElement& a = Get(of[0], where);
        const MeasureElement& b = Get(of[1], where);
        if (a.mKind != ElementKind::Plane || b.mKind != ElementKind::Plane)
            ms_refuse(where + ": 'mid_plane' takes two planes");
        const MsV3 nb = ms_dot(a.mDir, b.mDir) < 0.0 ? ms_scale(b.mDir, -1.0) : b.mDir;
        e.mKind = ElementKind::Plane;
        e.mDir = ms_unit(ms_add(a.mDir, nb));
        e.mPoint = ms_scale(ms_add(a.mPoint, b.mPoint), 0.5);
        return e;
    }
    if (how == "projection") {
        need(2);
        const MeasureElement& a = Get(of[0], where);
        const MeasureElement& pl = Get(of[1], where);
        if (pl.mKind != ElementKind::Plane)
            ms_refuse(where + ": 'projection' projects onto a plane (the second element)");
        e.mPoint =
            ms_sub(a.mPoint, ms_scale(pl.mDir, ms_dot(ms_sub(a.mPoint, pl.mPoint), pl.mDir)));
        if (ms_line_like(a)) {
            const MsV3 d = ms_sub(a.mDir, ms_scale(pl.mDir, ms_dot(a.mDir, pl.mDir)));
            if (ms_norm(d) < kMsParallel)
                ms_refuse(where + ": the line is perpendicular to the plane");
            e.mKind = ElementKind::Line;
            e.mDir = ms_unit(d);
        } else if (ms_point_like(a)) {
            e.mKind = ElementKind::Point;
        } else {
            ms_refuse(where + ": 'projection' takes a point or a line");
        }
        return e;
    }
    if (how == "plane_through") {
        need(3);
        const MsV3 p0 = Get(of[0], where).mPoint, p1 = Get(of[1], where).mPoint,
                   p2 = Get(of[2], where).mPoint;
        const MsV3 n = ms_cross(ms_sub(p1, p0), ms_sub(p2, p0));
        if (ms_norm(n) <
            kMsParallel * std::max(1.0, ms_norm(ms_sub(p1, p0)) * ms_norm(ms_sub(p2, p0))))
            ms_refuse(where + ": the three points are collinear");
        e.mKind = ElementKind::Plane;
        e.mDir = ms_unit(n);
        e.mPoint = p0;
        return e;
    }
    ms_refuse(where + ": unknown construction '" + how +
              "' (intersection, axis, centre, apex, mid_plane, projection, plane_through)");
}

MeasureCheck MsRun::Dimension(const JsonValue& rSpec, std::size_t K) {
    MeasureCheck c;
    c.mName = ms_text(rSpec, "name", "dimension", "dimension" + std::to_string(K + 1));
    const std::string where = "dimension '" + c.mName + "'";
    c.mType = ms_text(rSpec, "type", where);
    c.mOf = ms_names(rSpec, "between", where);
    if (c.mOf.size() != 2)
        ms_refuse(where + ": 'between' names two elements");
    const MeasureElement& a = Get(c.mOf[0], where);
    const MeasureElement& b = Get(c.mOf[1], where);
    if (c.mType == "distance") {
        const bool pa = ms_point_like(a), pb = ms_point_like(b);
        const bool la = ms_line_like(a), lb = ms_line_like(b);
        const bool qa = a.mKind == ElementKind::Plane, qb = b.mKind == ElementKind::Plane;
        if (pa && pb) {
            c.mValue = ms_norm(ms_sub(a.mPoint, b.mPoint));
        } else if ((pa && qb) || (qa && pb)) {
            const MeasureElement& p = pa ? a : b;
            const MeasureElement& q = pa ? b : a;
            c.mValue = std::fabs(ms_dot(ms_sub(p.mPoint, q.mPoint), q.mDir));
        } else if ((pa && lb) || (la && pb)) {
            const MeasureElement& p = pa ? a : b;
            const MeasureElement& l = pa ? b : a;
            c.mValue = ms_point_line(p.mPoint, l.mPoint, l.mDir);
        } else if (la && lb) {
            const MsV3 u = ms_cross(a.mDir, b.mDir);
            c.mValue = ms_norm(u) < kMsParallel
                           ? ms_point_line(b.mPoint, a.mPoint, a.mDir)
                           : std::fabs(ms_dot(ms_sub(b.mPoint, a.mPoint), u)) / ms_norm(u);
        } else if ((la && qb) || (qa && lb)) {
            const MeasureElement& l = la ? a : b;
            const MeasureElement& q = la ? b : a;
            if (std::fabs(ms_dot(l.mDir, q.mDir)) > kMsParallel)
                mWarnings.push_back(where + ": the line meets the plane; the distance is 0");
            c.mValue = std::fabs(ms_dot(l.mDir, q.mDir)) > kMsParallel
                           ? 0.0
                           : std::fabs(ms_dot(ms_sub(l.mPoint, q.mPoint), q.mDir));
        } else {
            if (ms_norm(ms_cross(a.mDir, b.mDir)) > 1e-6)
                mWarnings.push_back(
                    where +
                    ": the planes are not parallel; the distance is from the second's point");
            c.mValue = std::fabs(ms_dot(ms_sub(b.mPoint, a.mPoint), a.mDir));
        }
    } else if (c.mType == "angle") {
        if (ms_point_like(a) || ms_point_like(b))
            ms_refuse(where + ": an angle is between lines and planes");
        const double cosv = std::min(1.0, std::fabs(ms_dot(a.mDir, b.mDir)));
        double deg = std::acos(cosv) * 180.0 / kMsPi;
        const bool mixed = (a.mKind == ElementKind::Plane) != (b.mKind == ElementKind::Plane);
        if (mixed)
            deg = 90.0 - deg;
        c.mValue = deg;
    } else {
        ms_refuse(where + ": 'type' must be distance or angle");
    }
    if (const JsonValue* n = ms_find(rSpec, "nominal")) {
        if (!n->IsNumber())
            ms_refuse(where + ": 'nominal' must be a number");
        c.mHasNominal = true;
        c.mNominal = n->mNumber;
        c.mDeviation = c.mValue - c.mNominal;
    }
    if (const JsonValue* t = ms_find(rSpec, "tolerance")) {
        if (!t->IsNumber() || !(t->mNumber >= 0.0))
            ms_refuse(where + ": 'tolerance' must be a number >= 0");
        if (!c.mHasNominal)
            ms_refuse(where + ": a 'tolerance' needs a 'nominal'");
        c.mHasTolerance = true;
        c.mTolerance = t->mNumber;
        c.mPass = std::fabs(c.mDeviation) <= c.mTolerance;
    }
    return c;
}

MeasureCheck MsRun::Tolerance(const JsonValue& rSpec, std::size_t K) {
    MeasureCheck c;
    c.mName = ms_text(rSpec, "name", "tolerance", "tolerance" + std::to_string(K + 1));
    const std::string where = "tolerance '" + c.mName + "'";
    c.mType = ms_text(rSpec, "type", where);
    const std::string feature = ms_text(rSpec, "feature", where);
    const std::string datum = ms_text(rSpec, "datum", where);
    if (feature.empty())
        ms_refuse(where + ": it needs a 'feature'");
    const MeasureElement& f = Get(feature, where);
    c.mOf.push_back(feature);
    const std::string& t = c.mType;
    auto need_form = [&](ElementKind Kind) {
        if (f.mKind != Kind || !f.mMeasured)
            ms_refuse(where + ": '" + t + "' takes a measured " + ms_kind_name(Kind));
    };
    if (t == "flatness" || t == "cylindricity" || t == "sphericity") {
        need_form(t == "flatness"       ? ElementKind::Plane
                  : t == "cylindricity" ? ElementKind::Cylinder
                                        : ElementKind::Sphere);
        c.mValue = f.mFormError;
    } else if (t == "circularity") {
        // The largest radial spread about the axis over sections along it.
        need_form(ElementKind::Cylinder);
        const auto sections = static_cast<std::int64_t>(ms_number(rSpec, "sections", where, 10.0));
        if (sections < 1)
            ms_refuse(where + ": 'sections' must be at least 1");
        const auto [e0, e1] = ms_axis_ends(f);
        const double len = ms_norm(ms_sub(e1, e0));
        std::vector<double> lo(static_cast<std::size_t>(sections),
                               std::numeric_limits<double>::infinity());
        std::vector<double> hi(static_cast<std::size_t>(sections),
                               -std::numeric_limits<double>::infinity());
        for (std::size_t i = 0; i + 2 < f.mPoints.size(); i += 3) {
            const MsV3 x{f.mPoints[i], f.mPoints[i + 1], f.mPoints[i + 2]};
            const double s = len > 0.0 ? ms_dot(ms_sub(x, e0), f.mDir) / len : 0.0;
            const auto k = static_cast<std::size_t>(std::clamp<std::int64_t>(
                static_cast<std::int64_t>(s * static_cast<double>(sections)), 0, sections - 1));
            const double r = ms_point_line(x, f.mPoint, f.mDir);
            lo[k] = std::min(lo[k], r);
            hi[k] = std::max(hi[k], r);
        }
        c.mValue = 0.0;
        for (std::size_t k = 0; k < lo.size(); ++k)
            if (hi[k] >= lo[k])
                c.mValue = std::max(c.mValue, hi[k] - lo[k]);
    } else if (t == "parallelism" || t == "perpendicularity" || t == "angularity") {
        if (datum.empty())
            ms_refuse(where + ": '" + t + "' needs a 'datum'");
        const MeasureElement& d = Get(datum, where);
        c.mOf.push_back(datum);
        if (!(d.mKind == ElementKind::Plane || ms_line_like(d)))
            ms_refuse(where + ": the datum must be a plane or an axis");
        const double basic =
            t == "parallelism" ? 0.0
            : t == "perpendicularity"
                ? 90.0
                : ms_number(rSpec, "basic_angle", where, std::numeric_limits<double>::quiet_NaN());
        if (!std::isfinite(basic))
            ms_refuse(where + ": 'angularity' needs a 'basic_angle' in degrees");
        const double th = basic * kMsPi / 180.0;
        const MsV3 D = d.mDir;
        if (f.mKind == ElementKind::Plane) {
            if (f.mPoints.empty())
                ms_refuse(where + ": the feature must be a measured plane");
            // The zone's normal at the basic angle to the datum's normal (a
            // plane datum) or to the plane the datum axis is normal to.
            const double beta = d.mKind == ElementKind::Plane ? th : kMsPi / 2.0 - th;
            const MsV3 w = ms_perp(D, f.mDir);
            const MsV3 m = ms_unit(ms_add(ms_scale(D, std::cos(beta)), w, std::sin(beta)));
            c.mValue = ms_spread(f, m);
        } else if (ms_line_like(f)) {
            const auto [e0, e1] = ms_axis_ends(f);
            const MsV3 delta = ms_sub(e1, e0);
            const bool plane_datum = d.mKind == ElementKind::Plane;
            if ((t == "perpendicularity" && plane_datum) || (t == "parallelism" && !plane_datum)) {
                // A cylindrical zone along the datum direction.
                c.mValue = ms_norm(ms_sub(delta, ms_scale(D, ms_dot(delta, D))));
            } else {
                // Two planes: inclined at the basic angle to the datum plane,
                // or perpendicular to... the datum axis' at the basic angle.
                const double beta = plane_datum ? th : kMsPi / 2.0 - th;
                const MsV3 w = ms_perp(D, f.mDir);
                const MsV3 m = ms_unit(ms_add(ms_scale(D, std::cos(beta)), w, -std::sin(beta)));
                c.mValue = std::fabs(ms_dot(m, delta));
            }
        } else {
            ms_refuse(where + ": the feature must be a plane or an axis");
        }
    } else if (t == "coaxiality" || t == "concentricity") {
        if (datum.empty())
            ms_refuse(where + ": '" + t + "' needs a 'datum'");
        const MeasureElement& d = Get(datum, where);
        c.mOf.push_back(datum);
        auto dist = [&](const MsV3& rP) {
            if (ms_line_like(d))
                return ms_point_line(rP, d.mPoint, d.mDir);
            if (ms_point_like(d))
                return ms_norm(ms_sub(rP, d.mPoint));
            ms_refuse(where + ": the datum must be an axis or a point");
        };
        if (t == "coaxiality") {
            if (!ms_line_like(f))
                ms_refuse(where + ": 'coaxiality' takes an axis (a cylinder, a cone or a line)");
            const auto [e0, e1] = ms_axis_ends(f);
            c.mValue = 2.0 * std::max(dist(e0), dist(e1));
        } else {
            if (!(ms_point_like(f) || ms_line_like(f)))
                ms_refuse(where +
                          ": 'concentricity' takes a centre (a sphere, a point or an axis)");
            c.mValue = 2.0 * dist(f.mPoint);
        }
    } else {
        ms_refuse(where + ": unknown tolerance '" + t + "'");
    }
    if (const JsonValue* l = ms_find(rSpec, "limit")) {
        if (!l->IsNumber() || !(l->mNumber >= 0.0))
            ms_refuse(where + ": 'limit' must be a number >= 0");
        c.mHasTolerance = true;
        c.mTolerance = l->mNumber;
        c.mPass = c.mValue <= c.mTolerance;
    }
    return c;
}

JsonValue ms_vec_json(const MsV3& rV) {
    JsonValue a = JsonValue::Array();
    for (const double x : rV)
        a.Push(JsonValue::Number(x));
    return a;
}

JsonValue ms_number_or_null(double V) {
    return std::isfinite(V) ? JsonValue::Number(V) : JsonValue::Null();
}

JsonValue ms_element_json(const MeasureElement& rE) {
    JsonValue o = JsonValue::Object();
    o.Add("name", JsonValue::String(rE.mName));
    o.Add("kind", JsonValue::String(ms_kind_name(rE.mKind)));
    o.Add("point", ms_vec_json(rE.mPoint));
    o.Add("dir", ms_vec_json(rE.mDir));
    o.Add("radius", JsonValue::Number(rE.mRadius));
    o.Add("angle", JsonValue::Number(rE.mAngleDeg));
    o.Add("form_error", ms_number_or_null(rE.mFormError));
    o.Add("measured", JsonValue::Boolean(rE.mMeasured));
    o.Add("source", JsonValue::String(rE.mSource));
    JsonValue in = JsonValue::Array();
    for (const std::string& s : rE.mInputs)
        in.Push(JsonValue::String(s));
    o.Add("inputs", std::move(in));
    o.Add("num_points", JsonValue::Integer(static_cast<std::int64_t>(rE.mPoints.size() / 3)));
    return o;
}

JsonValue ms_check_json(const MeasureCheck& rC) {
    JsonValue o = JsonValue::Object();
    o.Add("name", JsonValue::String(rC.mName));
    o.Add("type", JsonValue::String(rC.mType));
    JsonValue of = JsonValue::Array();
    for (const std::string& s : rC.mOf)
        of.Push(JsonValue::String(s));
    o.Add("of", std::move(of));
    o.Add("value", JsonValue::Number(rC.mValue));
    o.Add("nominal", rC.mHasNominal ? JsonValue::Number(rC.mNominal) : JsonValue::Null());
    o.Add("tolerance", rC.mHasTolerance ? JsonValue::Number(rC.mTolerance) : JsonValue::Null());
    o.Add("deviation", rC.mHasNominal ? JsonValue::Number(rC.mDeviation) : JsonValue::Null());
    o.Add("pass", JsonValue::Boolean(rC.mPass));
    return o;
}

std::string ms_table_json(const std::vector<MeasureElement>& rElements) {
    JsonValue t = JsonValue::Object();
    JsonValue arr = JsonValue::Array();
    for (const MeasureElement& r_e : rElements)
        arr.Push(ms_element_json(r_e));
    t.Add("elements", std::move(arr));
    return detail::write_json(t, 0);
}

void ms_store_table(Mesh& rMesh, const std::vector<MeasureElement>& rElements) {
    const std::string text = ms_table_json(rElements);
    NDArray bytes = NDArray::Uninit(DType::UInt8, {text.size()});
    std::copy(text.begin(), text.end(), bytes.As<std::uint8_t>());
    rMesh.AddFieldData(kMeasureTableName, std::move(bytes));
}

}  // namespace

MeasureResult measure(const Mesh& rMesh, const std::string& rSpecJson) {
    JsonValue spec;
    try {
        spec = detail::parse_json(rSpecJson, "measure spec");
    } catch (const std::exception& rErr) {
        ms_refuse(std::string("the spec is not JSON: ") + rErr.what());
    }
    if (!spec.IsObject())
        ms_refuse("the spec is a JSON object (elements, dimensions, tolerances)");
    MsRun run{rMesh, FitCriterion::LeastSquares, {}, {}, {}, {}};
    const std::string crit = ms_text(spec, "criterion", "the spec", "least_squares");
    run.mCriterion = fit_criterion_from_name(crit);
    MeasureResult result;
    if (const JsonValue* els = ms_find(spec, "elements")) {
        if (!els->IsArray())
            ms_refuse("'elements' is a list");
        for (const JsonValue& r_e : els->mItems) {
            const std::string name = ms_text(r_e, "name", "an element");
            if (name.empty())
                ms_refuse("every element needs a 'name'");
            if (run.mIndex.count(name))
                ms_refuse("two elements are named '" + name + "'");
            MeasureElement e;
            if (ms_find(r_e, "fit")) {
                e = run.Fit(r_e, name);
            } else if (ms_find(r_e, "construct")) {
                e = run.Construct(r_e, name);
            } else if (const JsonValue* p = ms_find(r_e, "point")) {
                e.mName = name;
                e.mKind = ElementKind::Point;
                e.mPoint = ms_vec3(*p, "element '" + name + "'");
                e.mSource = "point";
                e.mFormError = std::numeric_limits<double>::quiet_NaN();
            } else {
                ms_refuse("element '" + name + "': it needs 'fit', 'construct' or 'point'");
            }
            run.mIndex.emplace(name, run.mElements.size());
            run.mElements.push_back(std::move(e));
        }
    }
    if (const JsonValue* dims = ms_find(spec, "dimensions")) {
        if (!dims->IsArray())
            ms_refuse("'dimensions' is a list");
        for (std::size_t k = 0; k < dims->mItems.size(); ++k)
            result.mDimensions.push_back(run.Dimension(dims->mItems[k], k));
    }
    if (const JsonValue* tols = ms_find(spec, "tolerances")) {
        if (!tols->IsArray())
            ms_refuse("'tolerances' is a list");
        for (std::size_t k = 0; k < tols->mItems.size(); ++k)
            result.mTolerances.push_back(run.Tolerance(tols->mItems[k], k));
    }
    for (const MeasureCheck& r_c : result.mDimensions)
        result.mNumFail += r_c.mPass ? 0 : 1;
    for (const MeasureCheck& r_c : result.mTolerances)
        result.mNumFail += r_c.mPass ? 0 : 1;

    result.mMesh = detail::clone_mesh(
        rMesh, [](DataLocation, const std::string&, std::string&) { return true; });
    for (auto& [name, cells] : run.mRegions) {
        NDArray e = NDArray::Uninit(DType::Int64, {cells.size()});
        std::copy(cells.begin(), cells.end(), e.As<std::int64_t>());
        result.mMesh.AddRegion(Region(name, RegionKind::Cell, std::move(e)));
    }
    ms_store_table(result.mMesh, run.mElements);
    result.mTableJson = ms_table_json(run.mElements);
    JsonValue rep = JsonValue::Object();
    JsonValue els = JsonValue::Array(), dims = JsonValue::Array(), tols = JsonValue::Array();
    for (const MeasureElement& r_e : run.mElements)
        els.Push(ms_element_json(r_e));
    for (const MeasureCheck& r_c : result.mDimensions)
        dims.Push(ms_check_json(r_c));
    for (const MeasureCheck& r_c : result.mTolerances)
        tols.Push(ms_check_json(r_c));
    rep.Add("elements", std::move(els));
    rep.Add("dimensions", std::move(dims));
    rep.Add("tolerances", std::move(tols));
    rep.Add("num_fail", JsonValue::Integer(result.mNumFail));
    JsonValue warns = JsonValue::Array();
    for (const std::string& w : run.mWarnings)
        warns.Push(JsonValue::String(w));
    rep.Add("warnings", std::move(warns));
    result.mReportJson = detail::write_json(rep, 0);
    result.mElements = std::move(run.mElements);
    result.mWarnings = std::move(run.mWarnings);
    return result;
}

std::vector<MeasureElement> measure_table(const Mesh& rMesh) {
    std::vector<MeasureElement> out;
    if (!rMesh.HasFieldData(kMeasureTableName))
        return out;
    const NDArray& bytes = rMesh.FieldData(kMeasureTableName);
    std::string text(bytes.Size(), '\0');
    for (std::size_t i = 0; i < bytes.Size(); ++i)
        text[i] = static_cast<char>(detail::read_int(bytes, i));
    const JsonValue t = detail::parse_json(text, "measure:table");
    const JsonValue* els = ms_find(t, "elements");
    if (!els || !els->IsArray())
        ms_refuse("the measure:table field data is not an element table");
    for (const JsonValue& r_e : els->mItems) {
        MeasureElement e;
        e.mName = ms_text(r_e, "name", "measure:table");
        e.mKind = ms_kind_from(ms_text(r_e, "kind", "measure:table"));
        if (const JsonValue* p = ms_find(r_e, "point"))
            e.mPoint = ms_vec3(*p, "measure:table");
        if (const JsonValue* d = ms_find(r_e, "dir"))
            e.mDir = ms_vec3(*d, "measure:table");
        e.mRadius = ms_number(r_e, "radius", "measure:table", 0.0);
        e.mAngleDeg = ms_number(r_e, "angle", "measure:table", 0.0);
        const JsonValue* fe = ms_find(r_e, "form_error");
        e.mFormError =
            fe && fe->IsNumber() ? fe->mNumber : std::numeric_limits<double>::quiet_NaN();
        const JsonValue* m = ms_find(r_e, "measured");
        e.mMeasured = m && m->IsBool() && m->mBool;
        e.mSource = ms_text(r_e, "source", "measure:table");
        e.mInputs = ms_names(r_e, "inputs", "measure:table");
        out.push_back(std::move(e));
    }
    return out;
}

MeasureElement measure_element(const Mesh& rMesh, const std::string& rName) {
    const std::vector<MeasureElement> table = measure_table(rMesh);
    if (table.empty())
        ms_refuse("the mesh carries no element table (run measure first)");
    for (const MeasureElement& r_e : table)
        if (r_e.mName == rName)
            return r_e;
    ms_refuse("no element named '" + rName + "'");
}

std::pair<std::int64_t, std::array<double, 10>> measure_element_primitive(
    const MeasureElement& rElement) {
    std::array<double, 10> p{};
    for (std::size_t k = 0; k < 3; ++k) {
        p[k] = rElement.mPoint[k];
        p[3 + k] = rElement.mDir[k];
    }
    switch (rElement.mKind) {
        case ElementKind::Plane:
            return {1, p};
        case ElementKind::Sphere:
            p[3] = p[4] = p[5] = 0.0;
            p[6] = rElement.mRadius;
            return {2, p};
        case ElementKind::Cylinder:
            p[6] = rElement.mRadius;
            return {3, p};
        case ElementKind::Cone:
            p[8] = rElement.mAngleDeg;
            return {4, p};
        default:
            ms_refuse("element '" + rElement.mName + "' is a " + ms_kind_name(rElement.mKind) +
                      ", which has no surface to measure against");
    }
}

AlignDatumsResult align_datums(const Mesh& rMesh, const std::string& rPrimary,
                               const std::string& rSecondary, const std::string& rOrigin) {
    std::vector<MeasureElement> table = measure_table(rMesh);
    if (table.empty())
        ms_refuse("align_datums: the mesh carries no element table (run measure first)");
    auto find = [&](const std::string& rName) -> const MeasureElement& {
        for (const MeasureElement& r_e : table)
            if (r_e.mName == rName)
                return r_e;
        ms_refuse("align_datums: no element named '" + rName + "'");
    };
    const MeasureElement& p = find(rPrimary);
    const MeasureElement& s = find(rSecondary);
    if (!(p.mKind == ElementKind::Plane || ms_line_like(p)))
        ms_refuse("align_datums: the primary must be a plane or an axis");
    if (!(s.mKind == ElementKind::Plane || ms_line_like(s)))
        ms_refuse("align_datums: the secondary must be a plane or an axis");
    const MsV3 z = ms_unit(p.mDir);
    const MsV3 xs = ms_sub(s.mDir, ms_scale(z, ms_dot(s.mDir, z)));
    if (ms_norm(xs) < 1e-9)
        ms_refuse("align_datums: the secondary is parallel to the primary");
    const MsV3 x = ms_unit(xs);
    const MsV3 y = ms_cross(z, x);
    MsV3 o = p.mPoint;
    if (!rOrigin.empty()) {
        const MeasureElement& q = find(rOrigin);
        if (q.mKind == ElementKind::Plane) {
            if (p.mKind == ElementKind::Plane && s.mKind == ElementKind::Plane) {
                // The meeting point of three planes (Cramer).
                const MsV3 n1 = p.mDir, n2 = s.mDir, n3 = q.mDir;
                const double det = ms_dot(n1, ms_cross(n2, n3));
                if (std::fabs(det) < 1e-12)
                    ms_refuse("align_datums: the three planes do not meet in a point");
                const double d1 = ms_dot(n1, p.mPoint), d2 = ms_dot(n2, s.mPoint),
                             d3 = ms_dot(n3, q.mPoint);
                o = ms_scale(
                    ms_add(ms_add(ms_scale(ms_cross(n2, n3), d1), ms_scale(ms_cross(n3, n1), d2)),
                           ms_scale(ms_cross(n1, n2), d3)),
                    1.0 / det);
            } else {
                const MeasureElement& l = ms_line_like(p) ? p : s;
                const double dn = ms_dot(l.mDir, q.mDir);
                if (std::fabs(dn) < 1e-12)
                    ms_refuse("align_datums: the axis is parallel to the origin plane");
                o = ms_add(l.mPoint, l.mDir, ms_dot(ms_sub(q.mPoint, l.mPoint), q.mDir) / dn);
            }
        } else {
            o = q.mPoint;
        }
    }
    AlignDatumsResult result;
    AffineTransform xf;
    const MsV3 rows[3] = {x, y, z};
    for (std::size_t r = 0; r < 3; ++r) {
        for (std::size_t c = 0; c < 3; ++c)
            xf.mMatrix[4 * r + c] = rows[r][c];
        xf.mMatrix[4 * r + 3] = -ms_dot(rows[r], o);
    }
    xf.mMatrix[12] = xf.mMatrix[13] = xf.mMatrix[14] = 0.0;
    xf.mMatrix[15] = 1.0;
    result.mMatrix = xf.mMatrix;
    result.mMesh = transform(rMesh, xf);
    for (MeasureElement& r_e : table) {
        const MsV3 d = ms_sub(r_e.mPoint, o);
        r_e.mPoint = {ms_dot(x, d), ms_dot(y, d), ms_dot(z, d)};
        r_e.mDir = {ms_dot(x, r_e.mDir), ms_dot(y, r_e.mDir), ms_dot(z, r_e.mDir)};
    }
    ms_store_table(result.mMesh, table);
    return result;
}

}  // namespace meshioplusplus
