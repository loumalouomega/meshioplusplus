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

// Kratos CAD JSON reader, sampler and writer. See formats/cad_json.hpp.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

// Project includes (private, not installed)
#include "../detail/brep_trim.hpp"
#include "../detail/json_value.hpp"
#include "../detail/nurbs_eval.hpp"
#include "lagrange_common.hpp"

namespace meshioplusplus {
namespace {

using detail::JsonValue;

constexpr const char* kCjFormat = "cad_json";
constexpr int kCjMaxDegree = 32;
constexpr double kCjCorner = 0.999999;
// Automatic subdivision: a clipped face gets about this many cells per direction.
constexpr std::size_t kCjAutoCells = 16;

// ---------------------------------------------------------------------------
// Reading: JSON -> BrepInfo
// ---------------------------------------------------------------------------

[[noreturn]] void cj_fail(const std::string& rWhat) {
    throw ReadError(std::string("meshio++: cad_json: ") + rWhat);
}

double cj_number(const JsonValue& rV, const std::string& rWhat) {
    if (!rV.IsNumber() || !std::isfinite(rV.mNumber))
        cj_fail(rWhat + " is not a finite number");
    return rV.mNumber;
}

std::int64_t cj_integer(const JsonValue& rV, const std::string& rWhat) {
    if (!rV.IsNumber())
        cj_fail(rWhat + " is not a number");
    if (rV.mIsInteger)
        return rV.mInteger;
    const double x = rV.mNumber;
    if (!std::isfinite(x) || x != std::floor(x) || std::fabs(x) > 9007199254740992.0)
        cj_fail(rWhat + " is not an integer");
    return static_cast<std::int64_t>(x);
}

bool cj_bool(const JsonValue& rV, const std::string& rWhat) {
    if (!rV.IsBool())
        cj_fail(rWhat + " is not a boolean");
    return rV.mBool;
}

const JsonValue& cj_array(const JsonValue& rV, const std::string& rWhat) {
    if (!rV.IsArray())
        cj_fail(rWhat + " is not an array");
    return rV;
}

const JsonValue& cj_object(const JsonValue& rV, const std::string& rWhat) {
    if (!rV.IsObject())
        cj_fail(rWhat + " is not an object");
    return rV;
}

// The members of @p rObj outside @p rKnown, as compact object text ("" for none).
std::string cj_extra(const JsonValue& rObj, std::initializer_list<std::string_view> Known) {
    JsonValue extra = JsonValue::Object();
    for (const auto& m : rObj.mMembers) {
        bool known = false;
        for (std::string_view k : Known)
            if (m.mKey == k) {
                known = true;
                break;
            }
        if (!known)
            extra.mMembers.push_back(m);
    }
    return extra.mMembers.empty() ? std::string() : detail::write_json(extra, 0);
}

// `brep_id` / `brep_name`.
void cj_identity(const JsonValue& rObj, const std::string& rWhat, std::int64_t& rId,
                 std::string& rName) {
    rId = -1;
    rName.clear();
    if (const JsonValue* id = rObj.Find("brep_id")) {
        rId = cj_integer(*id, rWhat + " brep_id");
        if (rId < 0)
            cj_fail(rWhat + " brep_id is negative");
    }
    if (const JsonValue* name = rObj.Find("brep_name")) {
        if (!name->IsString())
            cj_fail(rWhat + " brep_name is not a string");
        rName = name->mString;
    }
}

std::vector<double> cj_vector(const JsonValue& rV, const std::string& rWhat) {
    cj_array(rV, rWhat);
    std::vector<double> out;
    out.reserve(rV.mItems.size());
    for (std::size_t i = 0; i < rV.mItems.size(); ++i)
        out.push_back(cj_number(rV.mItems[i], rWhat + "[" + std::to_string(i) + "]"));
    return out;
}

void cj_check_knots(const std::vector<double>& rKnots, const std::string& rWhat) {
    for (std::size_t i = 1; i < rKnots.size(); ++i)
        if (rKnots[i] < rKnots[i - 1])
            cj_fail(rWhat + " decreases at entry " + std::to_string(i));
}

// `[[id, [x, y, z, w]], ...]` or `[[x, y, z, w], ...]`.
void cj_control_points(const JsonValue& rV, const std::string& rWhat, NDArray& rXyz,
                       NDArray& rWeights, NDArray& rIds) {
    cj_array(rV, rWhat);
    const std::size_t n = rV.mItems.size();
    rXyz = NDArray(DType::Float64, {n, 3});
    rWeights = NDArray(DType::Float64, {n});
    int form = -1;  // 1 with ids, 0 without
    std::vector<std::int64_t> ids;
    for (std::size_t i = 0; i < n; ++i) {
        const std::string what = rWhat + "[" + std::to_string(i) + "]";
        const JsonValue& cp = cj_array(rV.mItems[i], what);
        const JsonValue* coords = nullptr;
        int this_form;
        if (cp.mItems.size() == 2 && cp.mItems[1].IsArray()) {
            this_form = 1;
            ids.push_back(cj_integer(cp.mItems[0], what + " id"));
            coords = &cp.mItems[1];
        } else {
            this_form = 0;
            coords = &cp;
            if (cp.mItems.size() == 1 && cp.mItems[0].IsArray())
                coords = &cp.mItems[0];
        }
        if (form >= 0 && this_form != form)
            cj_fail(rWhat + " mixes control points with and without ids");
        form = this_form;
        if (coords->mItems.size() != 4)
            cj_fail(what + " needs [x, y, z, weight]");
        for (std::size_t d = 0; d < 3; ++d)
            rXyz.As<double>()[i * 3 + d] = cj_number(coords->mItems[d], what);
        rWeights.As<double>()[i] = cj_number(coords->mItems[3], what + " weight");
    }
    rIds = NDArray(DType::Int64, {ids.size()});
    std::copy(ids.begin(), ids.end(), rIds.As<std::int64_t>());
}

// A full knot vector from one in Kratos's reduced (n + p - 1) or the full form.
void cj_fit_knots(std::vector<double>& rKnots, bool Reduced) {
    if (!Reduced)
        return;
    std::vector<double> full;
    full.reserve(rKnots.size() + 2);
    full.push_back(rKnots.front());
    full.insert(full.end(), rKnots.begin(), rKnots.end());
    full.push_back(rKnots.back());
    rKnots = std::move(full);
}

std::int32_t cj_degree(const JsonValue& rV, const std::string& rWhat) {
    const std::int64_t p = cj_integer(rV, rWhat);
    if (p < 1 || p > kCjMaxDegree)
        cj_fail(rWhat + " " + std::to_string(p) + " is outside 1.." + std::to_string(kCjMaxDegree));
    return static_cast<std::int32_t>(p);
}

void cj_check_weights(const NDArray& rW, bool Rational, const std::string& rWhat) {
    if (!Rational)
        return;
    for (std::size_t i = 0; i < rW.Size(); ++i)
        if (!(rW.As<double>()[i] > 0.0))
            cj_fail(rWhat + " has a weight that is not positive");
}

BrepCurve cj_curve(const JsonValue& rObj, const std::string& rWhat) {
    cj_object(rObj, rWhat);
    BrepCurve c;
    if (const JsonValue* r = rObj.Find("is_rational"))
        c.mRational = cj_bool(*r, rWhat + " is_rational");
    const JsonValue* deg = rObj.Find("degree");
    const JsonValue* knots = rObj.Find("knot_vector");
    const JsonValue* cps = rObj.Find("control_points");
    if (!deg || !knots || !cps)
        cj_fail(rWhat + " needs degree, knot_vector and control_points");
    c.mDegree = cj_degree(*deg, rWhat + " degree");
    c.mKnots = cj_vector(*knots, rWhat + " knot_vector");
    cj_control_points(*cps, rWhat + " control_points", c.mControlPoints, c.mWeights, c.mIds);
    const std::size_t n = c.mWeights.Size();
    const std::size_t p = static_cast<std::size_t>(c.mDegree);
    if (n < p + 1)
        cj_fail(rWhat + " has fewer than degree + 1 control points");
    const std::size_t k = c.mKnots.size();
    if (k + 1 == n + p)
        c.mReducedKnots = true;
    else if (k != n + p + 1)
        cj_fail(rWhat + " has " + std::to_string(k) + " knots for " + std::to_string(n) +
                " control points of degree " + std::to_string(p));
    cj_check_knots(c.mKnots, rWhat + " knot_vector");
    cj_fit_knots(c.mKnots, c.mReducedKnots);
    cj_check_weights(c.mWeights, c.mRational, rWhat);
    if (const JsonValue* range = rObj.Find("active_range")) {
        c.mActiveRange = cj_vector(*range, rWhat + " active_range");
        if (c.mActiveRange.size() != 2)
            cj_fail(rWhat + " active_range needs two values");
    }
    c.mExtra = cj_extra(rObj, {"is_rational", "degree", "knot_vector", "active_range",
                               "control_points"});
    return c;
}

BrepTrim cj_trim(const JsonValue& rObj, const std::string& rWhat) {
    cj_object(rObj, rWhat);
    BrepTrim t;
    if (const JsonValue* idx = rObj.Find("trim_index")) {
        t.mTrimIndex = cj_integer(*idx, rWhat + " trim_index");
        if (t.mTrimIndex < 0)
            cj_fail(rWhat + " trim_index is negative");
    }
    if (const JsonValue* dir = rObj.Find("curve_direction"))
        t.mDirection = cj_bool(*dir, rWhat + " curve_direction");
    const JsonValue* curve = rObj.Find("parameter_curve");
    if (!curve)
        cj_fail(rWhat + " has no parameter_curve");
    t.mCurve = cj_curve(*curve, rWhat + " parameter_curve");
    t.mExtra = cj_extra(rObj, {"trim_index", "curve_direction", "parameter_curve"});
    return t;
}

BrepRef cj_ref(const JsonValue& rObj, const std::string& rWhat) {
    cj_object(rObj, rWhat);
    BrepRef r;
    cj_identity(rObj, rWhat, r.mId, r.mName);
    if (const JsonValue* idx = rObj.Find("trim_index")) {
        r.mTrimIndex = cj_integer(*idx, rWhat + " trim_index");
        if (r.mTrimIndex < 0)
            cj_fail(rWhat + " trim_index is negative");
    }
    if (const JsonValue* dir = rObj.Find("relative_direction"))
        r.mRelativeDirection = cj_bool(*dir, rWhat + " relative_direction") ? 1 : 0;
    if (const JsonValue* lc = rObj.Find("local_coordinates"))
        r.mLocalCoordinates = cj_vector(*lc, rWhat + " local_coordinates");
    r.mExtra = cj_extra(rObj, {"brep_id", "brep_name", "trim_index", "relative_direction",
                               "local_coordinates"});
    return r;
}

std::vector<BrepRef> cj_topology(const JsonValue& rObj, const std::string& rWhat) {
    std::vector<BrepRef> out;
    if (const JsonValue* topo = rObj.Find("topology")) {
        cj_array(*topo, rWhat + " topology");
        for (std::size_t i = 0; i < topo->mItems.size(); ++i)
            out.push_back(cj_ref(topo->mItems[i], rWhat + " topology[" + std::to_string(i) + "]"));
    }
    return out;
}

void cj_face(const JsonValue& rObj, const std::string& rWhat, BrepInfo& rInfo) {
    cj_object(rObj, rWhat);
    BrepFace f;
    cj_identity(rObj, rWhat, f.mId, f.mName);
    if (const JsonValue* sw = rObj.Find("swapped_surface_normal"))
        f.mSwappedNormal = cj_bool(*sw, rWhat + " swapped_surface_normal");
    const JsonValue* surface = rObj.Find("surface");
    if (!surface)
        cj_fail(rWhat + " has no surface");
    const std::string sw = rWhat + " surface";
    cj_object(*surface, sw);
    if (const JsonValue* t = surface->Find("is_trimmed"))
        f.mTrimmed = cj_bool(*t, sw + " is_trimmed");
    if (const JsonValue* r = surface->Find("is_rational"))
        f.mRational = cj_bool(*r, sw + " is_rational");
    const JsonValue* degs = surface->Find("degrees");
    const JsonValue* knots = surface->Find("knot_vectors");
    const JsonValue* cps = surface->Find("control_points");
    if (!degs || !knots || !cps)
        cj_fail(sw + " needs degrees, knot_vectors and control_points");
    cj_array(*degs, sw + " degrees");
    cj_array(*knots, sw + " knot_vectors");
    if (degs->mItems.size() != 2 || knots->mItems.size() != 2)
        cj_fail(sw + " needs two degrees and two knot vectors");
    NurbsPatch patch;
    for (std::size_t d = 0; d < 2; ++d) {
        patch.mDegrees.push_back(cj_degree(degs->mItems[d], sw + " degree"));
        patch.mKnots.push_back(cj_vector(knots->mItems[d], sw + " knot vector"));
        cj_check_knots(patch.mKnots.back(), sw + " knot vector");
    }
    cj_control_points(*cps, sw + " control_points", patch.mControlPoints, patch.mWeights,
                      f.mControlPointIds);
    const std::size_t n = patch.mWeights.Size();
    const std::size_t p = static_cast<std::size_t>(patch.mDegrees[0]);
    const std::size_t q = static_cast<std::size_t>(patch.mDegrees[1]);
    const std::size_t ku = patch.mKnots[0].size(), kv = patch.mKnots[1].size();
    // Kratos's own test (NurbsSurfaceGeometry::CheckAndFitKnotVectors): the
    // reduced form first, then the full one.
    if (ku + 1 > p && kv + 1 > q && (ku + 1 - p) * (kv + 1 - q) == n)
        f.mReducedKnots = true;
    else if (!(ku > p + 1 && kv > q + 1 && (ku - p - 1) * (kv - q - 1) == n))
        cj_fail(sw + ": " + std::to_string(n) + " control points do not match knot vectors of " +
                std::to_string(ku) + " and " + std::to_string(kv) + " knots at degrees " +
                std::to_string(p) + ", " + std::to_string(q));
    for (auto& k : patch.mKnots)
        cj_fit_knots(k, f.mReducedKnots);
    if (patch.mKnots[0].size() - p - 1 < p + 1 || patch.mKnots[1].size() - q - 1 < q + 1)
        cj_fail(sw + " has fewer than degree + 1 control points in a direction");
    cj_check_weights(patch.mWeights, f.mRational, sw);
    f.mSurfaceExtra = cj_extra(*surface, {"is_trimmed", "is_rational", "degrees",
                                          "knot_vectors", "control_points"});
    if (const JsonValue* loops = rObj.Find("boundary_loops")) {
        cj_array(*loops, rWhat + " boundary_loops");
        for (std::size_t l = 0; l < loops->mItems.size(); ++l) {
            const std::string lw = rWhat + " boundary_loops[" + std::to_string(l) + "]";
            const JsonValue& lo = cj_object(loops->mItems[l], lw);
            BrepLoop loop;
            if (const JsonValue* t = lo.Find("loop_type")) {
                if (!t->IsString())
                    cj_fail(lw + " loop_type is not a string");
                loop.mType = t->mString;
            } else {
                cj_fail(lw + " has no loop_type");
            }
            const JsonValue* trims = lo.Find("trimming_curves");
            if (!trims)
                cj_fail(lw + " has no trimming_curves");
            cj_array(*trims, lw + " trimming_curves");
            for (std::size_t t = 0; t < trims->mItems.size(); ++t)
                loop.mTrims.push_back(
                    cj_trim(trims->mItems[t], lw + " trimming_curves[" + std::to_string(t) + "]"));
            loop.mExtra = cj_extra(lo, {"loop_type", "trimming_curves"});
            f.mLoops.push_back(std::move(loop));
        }
    }
    if (const JsonValue* emb = rObj.Find("embedded_edges")) {
        cj_array(*emb, rWhat + " embedded_edges");
        for (std::size_t t = 0; t < emb->mItems.size(); ++t)
            f.mEmbeddedEdges.push_back(
                cj_trim(emb->mItems[t], rWhat + " embedded_edges[" + std::to_string(t) + "]"));
    }
    f.mExtra = cj_extra(rObj, {"brep_id", "brep_name", "swapped_surface_normal", "surface",
                               "boundary_loops", "embedded_edges"});
    f.mPatch = static_cast<std::int64_t>(rInfo.mPatches.size());
    rInfo.mPatches.push_back(std::move(patch));
    rInfo.mFaces.push_back(std::move(f));
}

void cj_read_breps(const JsonValue& rRoot, BrepInfo& rInfo) {
    const JsonValue& breps = cj_array(*rRoot.Find("breps"), "breps");
    for (std::size_t b = 0; b < breps.mItems.size(); ++b) {
        const std::string bw = "breps[" + std::to_string(b) + "]";
        const JsonValue& bo = cj_object(breps.mItems[b], bw);
        Brep brep;
        cj_identity(bo, bw, brep.mId, brep.mName);
        if (const JsonValue* faces = bo.Find("faces")) {
            cj_array(*faces, bw + " faces");
            for (std::size_t i = 0; i < faces->mItems.size(); ++i) {
                brep.mFaces.push_back(static_cast<std::int64_t>(rInfo.mFaces.size()));
                cj_face(faces->mItems[i], bw + " faces[" + std::to_string(i) + "]", rInfo);
            }
        }
        if (const JsonValue* edges = bo.Find("edges")) {
            cj_array(*edges, bw + " edges");
            for (std::size_t i = 0; i < edges->mItems.size(); ++i) {
                const std::string ew = bw + " edges[" + std::to_string(i) + "]";
                const JsonValue& eo = cj_object(edges->mItems[i], ew);
                BrepEdge e;
                cj_identity(eo, ew, e.mId, e.mName);
                if (const JsonValue* c = eo.Find("3d_curve")) {
                    e.mHasCurve = true;
                    e.mCurve = cj_curve(*c, ew + " 3d_curve");
                }
                e.mTopology = cj_topology(eo, ew);
                e.mExtra = cj_extra(eo, {"brep_id", "brep_name", "3d_curve", "topology"});
                brep.mEdges.push_back(static_cast<std::int64_t>(rInfo.mEdges.size()));
                rInfo.mEdges.push_back(std::move(e));
            }
        }
        if (const JsonValue* verts = bo.Find("vertices")) {
            cj_array(*verts, bw + " vertices");
            for (std::size_t i = 0; i < verts->mItems.size(); ++i) {
                const std::string vw = bw + " vertices[" + std::to_string(i) + "]";
                const JsonValue& vo = cj_object(verts->mItems[i], vw);
                BrepVertex v;
                cj_identity(vo, vw, v.mId, v.mName);
                v.mTopology = cj_topology(vo, vw);
                v.mExtra = cj_extra(vo, {"brep_id", "brep_name", "topology"});
                brep.mVertices.push_back(static_cast<std::int64_t>(rInfo.mVertices.size()));
                rInfo.mVertices.push_back(std::move(v));
            }
        }
        brep.mExtra = cj_extra(bo, {"brep_id", "brep_name", "faces", "edges", "vertices"});
        rInfo.mBreps.push_back(std::move(brep));
    }
    rInfo.mExtra = cj_extra(rRoot, {"breps"});
}

// `CPCoordinates` [[x, y, z]...] and `Weights` of the import dialect.
void cj_import_net(const JsonValue& rObj, const std::string& rWhat, NDArray& rXyz,
                   NDArray& rWeights) {
    const JsonValue* coords = rObj.Find("CPCoordinates");
    if (!coords)
        cj_fail(rWhat + " has no CPCoordinates");
    cj_array(*coords, rWhat + " CPCoordinates");
    const std::size_t n = coords->mItems.size();
    rXyz = NDArray(DType::Float64, {n, 3});
    for (std::size_t i = 0; i < n; ++i) {
        const std::vector<double> xyz = cj_vector(coords->mItems[i], rWhat + " CPCoordinates");
        if (xyz.size() < 2 || xyz.size() > 3)
            cj_fail(rWhat + " CPCoordinates need two or three values");
        for (std::size_t d = 0; d < 3; ++d)
            rXyz.As<double>()[i * 3 + d] = d < xyz.size() ? xyz[d] : 0.0;
    }
    rWeights = NDArray(DType::Float64, {n});
    const JsonValue* w = rObj.Find("Weights");
    if (!w) {
        for (std::size_t i = 0; i < n; ++i)
            rWeights.As<double>()[i] = 1.0;
        return;
    }
    const std::vector<double> ws = cj_vector(*w, rWhat + " Weights");
    if (ws.size() != n)
        cj_fail(rWhat + " has " + std::to_string(ws.size()) + " weights for " + std::to_string(n) +
                " control points");
    std::copy(ws.begin(), ws.end(), rWeights.As<double>());
}

bool cj_import_rational(const NDArray& rW) {
    bool all_one = true, all_zero = true;
    for (std::size_t i = 0; i < rW.Size(); ++i) {
        all_one = all_one && rW.As<double>()[i] == 1.0;
        all_zero = all_zero && rW.As<double>()[i] == 0.0;
    }
    return !all_one && !all_zero;
}

void cj_read_import(const JsonValue& rRoot, BrepInfo& rInfo) {
    Brep brep;
    if (const JsonValue* surfaces = rRoot.Find("Surfaces")) {
        cj_array(*surfaces, "Surfaces");
        for (std::size_t i = 0; i < surfaces->mItems.size(); ++i) {
            const std::string sw = "Surfaces[" + std::to_string(i) + "]";
            const JsonValue& so = cj_object(surfaces->mItems[i], sw);
            BrepFace f;
            NurbsPatch patch;
            if (const JsonValue* id = so.Find("id"))
                f.mId = cj_integer(*id, sw + " id");
            const JsonValue* pu = so.Find("pDegreeU");
            const JsonValue* pv = so.Find("pDegreeV");
            const JsonValue* kvs = so.Find("knotVectors");
            if (!pu || !pv || !kvs)
                cj_fail(sw + " needs pDegreeU, pDegreeV and knotVectors");
            patch.mDegrees = {cj_degree(*pu, sw + " pDegreeU"), cj_degree(*pv, sw + " pDegreeV")};
            cj_object(*kvs, sw + " knotVectors");
            const JsonValue* xi = kvs->Find("Xi");
            const JsonValue* eta = kvs->Find("Eta");
            if (!xi || !eta)
                cj_fail(sw + " knotVectors needs Xi and Eta");
            patch.mKnots = {cj_vector(*xi, sw + " Xi"), cj_vector(*eta, sw + " Eta")};
            for (const auto& k : patch.mKnots)
                cj_check_knots(k, sw + " knot vector");
            cj_import_net(so, sw, patch.mControlPoints, patch.mWeights);
            const std::size_t n = patch.mWeights.Size();
            const std::size_t p = static_cast<std::size_t>(patch.mDegrees[0]);
            const std::size_t q = static_cast<std::size_t>(patch.mDegrees[1]);
            const std::size_t ku = patch.mKnots[0].size(), kv = patch.mKnots[1].size();
            if (ku + 1 > p && kv + 1 > q && (ku + 1 - p) * (kv + 1 - q) == n)
                f.mReducedKnots = true;
            else if (!(ku > p + 1 && kv > q + 1 && (ku - p - 1) * (kv - q - 1) == n))
                cj_fail(sw + ": the control points do not match the knot vectors");
            for (auto& k : patch.mKnots)
                cj_fit_knots(k, f.mReducedKnots);
            if (patch.mKnots[0].size() - p - 1 < p + 1 || patch.mKnots[1].size() - q - 1 < q + 1)
                cj_fail(sw + " has fewer than degree + 1 control points in a direction");
            // ImportNurbsSbmModeler reads no surface weights; all zero or all
            // one is a polynomial surface.
            f.mRational = cj_import_rational(patch.mWeights);
            cj_check_weights(patch.mWeights, f.mRational, sw);
            f.mExtra = cj_extra(so, {"id", "pDegreeU", "pDegreeV", "knotVectors",
                                     "CPCoordinates", "Weights"});
            f.mPatch = static_cast<std::int64_t>(rInfo.mPatches.size());
            rInfo.mPatches.push_back(std::move(patch));
            brep.mFaces.push_back(static_cast<std::int64_t>(rInfo.mFaces.size()));
            rInfo.mFaces.push_back(std::move(f));
        }
    }
    if (const JsonValue* lines = rRoot.Find("Lines")) {
        cj_array(*lines, "Lines");
        for (std::size_t i = 0; i < lines->mItems.size(); ++i) {
            const std::string lw = "Lines[" + std::to_string(i) + "]";
            const JsonValue& lo = cj_object(lines->mItems[i], lw);
            BrepEdge e;
            e.mHasCurve = true;
            if (const JsonValue* id = lo.Find("id"))
                e.mId = cj_integer(*id, lw + " id");
            const JsonValue* deg = lo.Find("pDegree");
            const JsonValue* kv = lo.Find("knotVector");
            if (!deg || !kv)
                cj_fail(lw + " needs pDegree and knotVector");
            BrepCurve& c = e.mCurve;
            c.mDegree = cj_degree(*deg, lw + " pDegree");
            c.mKnots = cj_vector(*kv, lw + " knotVector");
            cj_check_knots(c.mKnots, lw + " knotVector");
            cj_import_net(lo, lw, c.mControlPoints, c.mWeights);
            const std::size_t n = c.mWeights.Size();
            const std::size_t p = static_cast<std::size_t>(c.mDegree);
            if (n < p + 1)
                cj_fail(lw + " has fewer than degree + 1 control points");
            if (c.mKnots.size() + 1 == n + p)
                c.mReducedKnots = true;
            else if (c.mKnots.size() != n + p + 1)
                cj_fail(lw + ": the control points do not match the knot vector");
            cj_fit_knots(c.mKnots, c.mReducedKnots);
            // ImportNurbsSbmModeler builds every line as a rational curve.
            c.mRational = true;
            cj_check_weights(c.mWeights, true, lw);
            e.mExtra = cj_extra(lo, {"id", "pDegree", "knotVector", "CPCoordinates", "Weights"});
            brep.mEdges.push_back(static_cast<std::int64_t>(rInfo.mEdges.size()));
            rInfo.mEdges.push_back(std::move(e));
        }
    }
    rInfo.mBreps.push_back(std::move(brep));
    rInfo.mExtra = cj_extra(rRoot, {"Surfaces", "Lines"});
}

BrepInfo cj_parse(const std::string& rPath) {
    const std::size_t bytes = detail::file_bytes(rPath);
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("meshio++: cad_json: cannot open " + rPath);
    std::string text(bytes, '\0');
    in.read(text.data(), static_cast<std::streamsize>(bytes));
    text.resize(static_cast<std::size_t>(in.gcount()));
    const JsonValue root = detail::parse_json(text, kCjFormat);
    if (!root.IsObject())
        cj_fail("the document is not an object");
    BrepInfo info;
    if (root.Find("breps"))
        cj_read_breps(root, info);
    else if (root.Find("Surfaces") || root.Find("Lines"))
        cj_read_import(root, info);
    else
        cj_fail("neither a \"breps\" list (CadJsonInput) nor \"Surfaces\"/\"Lines\" "
                "(ImportNurbsSbmModeler)");
    return info;
}

// ---------------------------------------------------------------------------
// Sampling: BrepInfo -> Mesh
// ---------------------------------------------------------------------------

[[noreturn]] void cj_bad(const std::string& rWhat) {
    throw std::invalid_argument("meshio++: sample_brep: " + rWhat);
}

std::string cj_quad_type(int Order) {
    return Order <= 1 ? "quad" : Order == 2 ? "quad9" : "VTK_LAGRANGE_QUADRILATERAL";
}

std::string cj_line_type(int Order) {
    return Order <= 1 ? "line" : Order == 2 ? "line3" : "VTK_LAGRANGE_CURVE";
}

void cj_check_curve(const BrepCurve& rC, const std::string& rWhat) {
    detail::brep_check_curve(rC, rWhat, "meshio++: sample_brep: ");
}

void cj_check_patch(const NurbsPatch& rP, const std::string& rWhat) {
    if (rP.mDegrees.size() != 2 || rP.mKnots.size() != 2)
        cj_bad(rWhat + " is not a surface patch (two degrees, two knot vectors)");
    std::size_t n = 1;
    for (std::size_t d = 0; d < 2; ++d) {
        if (rP.mDegrees[d] < 1 || rP.mDegrees[d] > kCjMaxDegree)
            cj_bad(rWhat + " degree is out of range");
        const auto p = static_cast<std::size_t>(rP.mDegrees[d]);
        if (rP.mKnots[d].size() < 2 * p + 2)
            cj_bad(rWhat + " knot vector is too short");
        for (std::size_t i = 1; i < rP.mKnots[d].size(); ++i)
            if (!(rP.mKnots[d][i] >= rP.mKnots[d][i - 1]))
                cj_bad(rWhat + " knot vector decreases");
        n *= rP.mKnots[d].size() - p - 1;
    }
    if (rP.mControlPoints.Dtype() != DType::Float64 || rP.mControlPoints.Shape().size() != 2 ||
        rP.mControlPoints.Shape()[0] != n || rP.mControlPoints.Shape()[1] != 3 ||
        rP.mWeights.Dtype() != DType::Float64 || rP.mWeights.Size() != n)
        cj_bad(rWhat + " control net does not match its knot vectors");
}

struct CjBlock {
    std::string mType;
    std::vector<std::int64_t> mConn;
    std::size_t mNodes = 0;
    std::vector<std::int64_t> mKind, mIndex, mId, mSpan, mTrimmed;
};

CjBlock& cj_block(std::vector<CjBlock>& rBlocks, const std::string& rType, std::size_t Nodes) {
    // Keyed by node count too: faces of degrees 3 and 4 are both
    // VTK_LAGRANGE_QUADRILATERAL, with 16 and 25 nodes.
    for (auto& b : rBlocks)
        if (b.mType == rType && b.mNodes == Nodes)
            return b;
    rBlocks.push_back(CjBlock{});
    rBlocks.back().mType = rType;
    rBlocks.back().mNodes = Nodes;
    return rBlocks.back();
}

// What one face or edge contributes before its points are evaluated: the
// lattice points it uses (as lattice indices) and its cells.
struct CjPending {
    bool mFace = true;
    std::size_t mEntity = 0;
    std::vector<double> mU, mV;              // lattice coordinates
    std::vector<std::size_t> mSpanU, mSpanV;  // span per lattice line
    std::size_t mNu = 0;                      // lattice width
    std::vector<std::size_t> mPoints;         // lattice index of each new point
};

void cj_check_options(const CadJsonReadOptions& rOpt) {
    if (rOpt.mSubdivisions < 0 || rOpt.mSubdivisions > 64)
        cj_bad("subdivisions must be 0 (automatic) to 64");
    if (rOpt.mTrimSamples < 2 || rOpt.mTrimSamples > 4096)
        cj_bad("trim_samples must be 2 to 4096");
}

Mesh cj_sample(const BrepInfo& rInfo, const CadJsonReadOptions& rOpt) {
    cj_check_options(rOpt);
    // Automatic (0): knot-span cells, except that a clipped face is refined
    // to about kCjAutoCells cells per direction so its trim shows.
    const int S = std::max(rOpt.mSubdivisions, 1);
    auto auto_split = [&](int P, const std::vector<double>& rKnots, std::size_t N, bool Clip) {
        if (rOpt.mSubdivisions > 0)
            return rOpt.mSubdivisions;
        if (!Clip)
            return 1;
        const std::size_t spans = std::max<std::size_t>(detail::nurbs_spans(P, rKnots, N).size(), 1);
        return static_cast<int>(std::min<std::size_t>((kCjAutoCells + spans - 1) / spans, 64));
    };
    for (std::size_t i = 0; i < rInfo.mPatches.size(); ++i)
        cj_check_patch(rInfo.mPatches[i], "patch " + std::to_string(i));
    std::vector<CjBlock> blocks;
    std::vector<CjPending> pending;
    std::size_t num_points = 0;

    for (std::size_t fi = 0; fi < rInfo.mFaces.size(); ++fi) {
        const BrepFace& face = rInfo.mFaces[fi];
        const std::string fw = "face " + std::to_string(fi);
        if (face.mPatch < 0 || static_cast<std::size_t>(face.mPatch) >= rInfo.mPatches.size())
            cj_bad(fw + " names no patch");
        const NurbsPatch& patch = rInfo.mPatches[static_cast<std::size_t>(face.mPatch)];
        const int p = patch.mDegrees[0], q = patch.mDegrees[1];
        const std::size_t nu = patch.mKnots[0].size() - static_cast<std::size_t>(p) - 1;
        const std::size_t nv = patch.mKnots[1].size() - static_cast<std::size_t>(q) - 1;
        for (const BrepLoop& loop : face.mLoops)
            for (const BrepTrim& t : loop.mTrims)
                cj_check_curve(t.mCurve, fw + " trimming curve");
        const int k = std::max(p, q);
        const bool clip = rOpt.mTrim && face.mTrimmed && !face.mLoops.empty();
        const detail::BrepAxis au = detail::brep_axis(p, patch.mKnots[0], nu, auto_split(p, patch.mKnots[0], nu, clip),
                                  nullptr);
        const detail::BrepAxis av = detail::brep_axis(q, patch.mKnots[1], nv, auto_split(q, patch.mKnots[1], nv, clip),
                                  nullptr);
        CjPending pd;
        pd.mFace = true;
        pd.mEntity = fi;
        pd.mU = detail::brep_lattice(au, k);
        pd.mV = detail::brep_lattice(av, k);
        pd.mNu = pd.mU.size();
        for (std::size_t i = 0; i < pd.mU.size(); ++i)
            pd.mSpanU.push_back(au.mSpan.empty() ? 0 : au.mSpan[std::min(i / k, au.mSpan.size() - 1)]);
        for (std::size_t i = 0; i < pd.mV.size(); ++i)
            pd.mSpanV.push_back(av.mSpan.empty() ? 0 : av.mSpan[std::min(i / k, av.mSpan.size() - 1)]);
        std::vector<std::vector<std::array<double, 2>>> loops;
        if (clip)
            for (const BrepLoop& loop : face.mLoops) {
                auto poly = detail::brep_loop_polygon(loop, rOpt.mTrimSamples);
                if (poly.size() >= 3)
                    loops.push_back(std::move(poly));
            }
        const auto lattice = lagrange::vtk_lattice(lagrange::Shape::Quad, k);
        std::vector<std::int64_t> local(au.mLo.empty() || av.mLo.empty() ? 0 : pd.mU.size() * pd.mV.size(), -1);
        CjBlock& blk = cj_block(blocks, cj_quad_type(k), lattice.size());
        for (std::size_t jv = 0; jv < av.mLo.size(); ++jv)
            for (std::size_t iu = 0; iu < au.mLo.size(); ++iu) {
                std::int64_t trimmed = 0;
                if (clip) {
                    const double u0 = au.mLo[iu], u1 = au.mHi[iu];
                    const double v0 = av.mLo[jv], v1 = av.mHi[jv];
                    const double uc = (u0 + u1) * 0.5, vc = (v0 + v1) * 0.5;
                    if (!detail::brep_inside(loops, uc, vc))
                        continue;
                    // The corners pulled a hair towards the centre, so a cell
                    // whose side lies on a loop is not counted as crossed.
                    const double a0 = uc + (u0 - uc) * kCjCorner, a1 = uc + (u1 - uc) * kCjCorner;
                    const double b0 = vc + (v0 - vc) * kCjCorner, b1 = vc + (v1 - vc) * kCjCorner;
                    if (!detail::brep_inside(loops, a0, b0) || !detail::brep_inside(loops, a1, b0) ||
                        !detail::brep_inside(loops, a1, b1) || !detail::brep_inside(loops, a0, b1))
                        trimmed = 1;
                }
                for (const auto& ijk : lattice) {
                    const std::size_t gi = iu * static_cast<std::size_t>(k) + static_cast<std::size_t>(ijk[0]);
                    const std::size_t gj = jv * static_cast<std::size_t>(k) + static_cast<std::size_t>(ijk[1]);
                    const std::size_t g = gj * pd.mNu + gi;
                    if (local[g] < 0) {
                        local[g] = static_cast<std::int64_t>(num_points + pd.mPoints.size());
                        pd.mPoints.push_back(g);
                    }
                    blk.mConn.push_back(local[g]);
                }
                blk.mKind.push_back(0);
                blk.mIndex.push_back(static_cast<std::int64_t>(fi));
                blk.mId.push_back(face.mId);
                blk.mSpan.push_back(au.mSpanNo[iu]);
                blk.mSpan.push_back(av.mSpanNo[jv]);
                blk.mTrimmed.push_back(trimmed);
            }
        num_points += pd.mPoints.size();
        pending.push_back(std::move(pd));
    }

    if (rOpt.mEdges)
        for (std::size_t ei = 0; ei < rInfo.mEdges.size(); ++ei) {
            const BrepEdge& edge = rInfo.mEdges[ei];
            if (!edge.mHasCurve)
                continue;
            const BrepCurve& c = edge.mCurve;
            cj_check_curve(c, "edge " + std::to_string(ei) + " 3d_curve");
            const int k = c.mDegree;
            const detail::BrepAxis axis = detail::brep_axis(k, c.mKnots, c.mWeights.Size(), S, &c.mActiveRange);
            CjPending pd;
            pd.mFace = false;
            pd.mEntity = ei;
            pd.mU = detail::brep_lattice(axis, k);
            pd.mNu = pd.mU.size();
            for (std::size_t i = 0; i < pd.mU.size(); ++i)
                pd.mSpanU.push_back(axis.mSpan.empty() ? 0
                                                       : axis.mSpan[std::min(i / k, axis.mSpan.size() - 1)]);
            const auto lattice = lagrange::vtk_lattice(lagrange::Shape::Line, k);
            std::vector<std::int64_t> local(pd.mU.size(), -1);
            CjBlock& blk = cj_block(blocks, cj_line_type(k), lattice.size());
            for (std::size_t iu = 0; iu < axis.mLo.size(); ++iu) {
                for (const auto& ijk : lattice) {
                    const std::size_t g = iu * static_cast<std::size_t>(k) + static_cast<std::size_t>(ijk[0]);
                    if (local[g] < 0) {
                        local[g] = static_cast<std::int64_t>(num_points + pd.mPoints.size());
                        pd.mPoints.push_back(g);
                    }
                    blk.mConn.push_back(local[g]);
                }
                blk.mKind.push_back(1);
                blk.mIndex.push_back(static_cast<std::int64_t>(ei));
                blk.mId.push_back(edge.mId);
                blk.mSpan.push_back(axis.mSpanNo[iu]);
                blk.mSpan.push_back(-1);
                blk.mTrimmed.push_back(0);
            }
            num_points += pd.mPoints.size();
            pending.push_back(std::move(pd));
        }

    NDArray points(DType::Float64, {num_points, 3});
    double* xyz = points.As<double>();
    std::size_t offset = 0;
    for (const CjPending& pd : pending) {
        double* out = xyz + offset * 3;
        if (pd.mFace) {
            const BrepFace& face = rInfo.mFaces[pd.mEntity];
            const NurbsPatch& patch = rInfo.mPatches[static_cast<std::size_t>(face.mPatch)];
            const std::size_t nu = patch.mKnots[0].size() - static_cast<std::size_t>(patch.mDegrees[0]) - 1;
            const double* cps = patch.mControlPoints.As<double>();
            const double* w = face.mRational ? patch.mWeights.As<double>() : nullptr;
            parallel_for(pd.mPoints.size(), [&](std::size_t i) {
                const std::size_t g = pd.mPoints[i];
                const std::size_t gi = g % pd.mNu, gj = g / pd.mNu;
                detail::nurbs_surface_point_in_span(patch.mDegrees[0], patch.mDegrees[1],
                                                    patch.mKnots[0], patch.mKnots[1], cps, w, nu, 3,
                                                    pd.mSpanU[gi], pd.mSpanV[gj], pd.mU[gi],
                                                    pd.mV[gj], out + i * 3);
            });
        } else {
            const BrepCurve& c = rInfo.mEdges[pd.mEntity].mCurve;
            const double* cps = c.mControlPoints.As<double>();
            const double* w = c.mRational ? c.mWeights.As<double>() : nullptr;
            parallel_for(pd.mPoints.size(), [&](std::size_t i) {
                const std::size_t g = pd.mPoints[i];
                detail::nurbs_curve_point_in_span(c.mDegree, c.mKnots, cps, w, 3, pd.mSpanU[g],
                                                  pd.mU[g], out + i * 3);
            });
        }
        offset += pd.mPoints.size();
    }

    Mesh mesh;
    mesh.AssignPoints(std::move(points));
    std::vector<NDArray> kind, index, id, span, trimmed;
    // A face trimmed away or an edge outside its range leaves no cells: keep
    // only the blocks that have some.
    blocks.erase(std::remove_if(blocks.begin(), blocks.end(),
                                [](const CjBlock& rB) { return rB.mKind.empty(); }),
                 blocks.end());
    for (CjBlock& b : blocks) {
        const std::size_t n = b.mKind.size();
        NDArray conn(DType::Int64, {n, b.mNodes});
        std::copy(b.mConn.begin(), b.mConn.end(), conn.As<std::int64_t>());
        mesh.AddCellBlock(b.mType, std::move(conn));
        auto column = [n](const std::vector<std::int64_t>& rV, std::size_t Width) {
            NDArray a(DType::Int64, Width == 1 ? std::vector<std::size_t>{n}
                                               : std::vector<std::size_t>{n, Width});
            std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
            return a;
        };
        kind.push_back(column(b.mKind, 1));
        index.push_back(column(b.mIndex, 1));
        id.push_back(column(b.mId, 1));
        span.push_back(column(b.mSpan, 2));
        trimmed.push_back(column(b.mTrimmed, 1));
    }
    if (!blocks.empty()) {
        mesh.AddCellData("cad:id", std::move(id));
        mesh.AddCellData("cad:index", std::move(index));
        mesh.AddCellData("cad:kind", std::move(kind));
        mesh.AddCellData("cad:span", std::move(span));
        mesh.AddCellData("cad:trimmed", std::move(trimmed));
    }
    return mesh;
}

// ---------------------------------------------------------------------------
// Writing: BrepInfo -> JSON
// ---------------------------------------------------------------------------

[[noreturn]] void cj_wfail(const std::string& rWhat) {
    throw WriteError("meshio++: cad_json: " + rWhat);
}

// Append the members of compact object text @p rExtra to @p rObj.
void cj_add_extra(JsonValue& rObj, const std::string& rExtra, const std::string& rWhat) {
    if (rExtra.empty())
        return;
    JsonValue extra;
    try {
        extra = detail::parse_json(rExtra, kCjFormat);
    } catch (const ReadError&) {
        cj_wfail(rWhat + " extra is not JSON object text");
    }
    if (!extra.IsObject())
        cj_wfail(rWhat + " extra is not JSON object text");
    for (auto& m : extra.mMembers)
        rObj.mMembers.push_back(std::move(m));
}

void cj_add_identity(JsonValue& rObj, std::int64_t Id, const std::string& rName) {
    if (Id >= 0)
        rObj.Add("brep_id", JsonValue::Integer(Id));
    if (!rName.empty())
        rObj.Add("brep_name", JsonValue::String(rName));
}

JsonValue cj_numbers(const std::vector<double>& rV) {
    JsonValue a = JsonValue::Array();
    for (double x : rV)
        a.Push(JsonValue::Number(x));
    return a;
}

std::vector<double> cj_written_knots(const std::vector<double>& rKnots, bool Reduced) {
    if (!Reduced || rKnots.size() < 2)
        return rKnots;
    return std::vector<double>(rKnots.begin() + 1, rKnots.end() - 1);
}

JsonValue cj_net(const NDArray& rXyz, const NDArray& rW, const NDArray& rIds,
                 const std::string& rWhat) {
    const std::size_t n = rW.Size();
    if (rXyz.Dtype() != DType::Float64 || rXyz.Shape().size() != 2 || rXyz.Shape()[0] != n ||
        rXyz.Shape()[1] != 3 || rW.Dtype() != DType::Float64)
        cj_wfail(rWhat + " control points are not Float64 (n, 3) beside n weights");
    if (rIds.Size() != 0 && (rIds.Size() != n || rIds.Dtype() != DType::Int64))
        cj_wfail(rWhat + " ids do not match its control points");
    JsonValue a = JsonValue::Array();
    for (std::size_t i = 0; i < n; ++i) {
        JsonValue c = JsonValue::Array();
        for (std::size_t d = 0; d < 3; ++d)
            c.Push(JsonValue::Number(rXyz.As<double>()[i * 3 + d]));
        c.Push(JsonValue::Number(rW.As<double>()[i]));
        if (rIds.Size()) {
            JsonValue e = JsonValue::Array();
            e.Push(JsonValue::Integer(rIds.As<std::int64_t>()[i]));
            e.Push(std::move(c));
            a.Push(std::move(e));
        } else {
            JsonValue e = JsonValue::Array();
            e.Push(std::move(c));
            a.Push(std::move(e));
        }
    }
    return a;
}

JsonValue cj_curve_json(const BrepCurve& rC, const std::string& rWhat) {
    try {
        cj_check_curve(rC, rWhat);
    } catch (const std::invalid_argument& e) {
        cj_wfail(e.what());
    }
    JsonValue o = JsonValue::Object();
    o.Add("is_rational", JsonValue::Boolean(rC.mRational));
    o.Add("degree", JsonValue::Integer(rC.mDegree));
    o.Add("knot_vector", cj_numbers(cj_written_knots(rC.mKnots, rC.mReducedKnots)));
    if (!rC.mActiveRange.empty()) {
        if (rC.mActiveRange.size() != 2)
            cj_wfail(rWhat + " active range needs two values");
        o.Add("active_range", cj_numbers(rC.mActiveRange));
    }
    o.Add("control_points", cj_net(rC.mControlPoints, rC.mWeights, rC.mIds, rWhat));
    cj_add_extra(o, rC.mExtra, rWhat);
    return o;
}

JsonValue cj_trim_json(const BrepTrim& rT, const std::string& rWhat) {
    JsonValue o = JsonValue::Object();
    if (rT.mTrimIndex >= 0)
        o.Add("trim_index", JsonValue::Integer(rT.mTrimIndex));
    o.Add("curve_direction", JsonValue::Boolean(rT.mDirection));
    o.Add("parameter_curve", cj_curve_json(rT.mCurve, rWhat + " parameter curve"));
    cj_add_extra(o, rT.mExtra, rWhat);
    return o;
}

JsonValue cj_topology_json(const std::vector<BrepRef>& rRefs, const std::string& rWhat) {
    JsonValue a = JsonValue::Array();
    for (const BrepRef& r : rRefs) {
        JsonValue o = JsonValue::Object();
        cj_add_identity(o, r.mId, r.mName);
        if (r.mTrimIndex >= 0)
            o.Add("trim_index", JsonValue::Integer(r.mTrimIndex));
        if (r.mRelativeDirection >= 0)
            o.Add("relative_direction", JsonValue::Boolean(r.mRelativeDirection != 0));
        if (!r.mLocalCoordinates.empty())
            o.Add("local_coordinates", cj_numbers(r.mLocalCoordinates));
        cj_add_extra(o, r.mExtra, rWhat);
        a.Push(std::move(o));
    }
    return a;
}

JsonValue cj_face_json(const BrepInfo& rInfo, const BrepFace& rF, const std::string& rWhat) {
    if (rF.mPatch < 0 || static_cast<std::size_t>(rF.mPatch) >= rInfo.mPatches.size())
        cj_wfail(rWhat + " names no patch");
    const NurbsPatch& patch = rInfo.mPatches[static_cast<std::size_t>(rF.mPatch)];
    try {
        cj_check_patch(patch, rWhat + " patch");
    } catch (const std::invalid_argument& e) {
        cj_wfail(e.what());
    }
    JsonValue o = JsonValue::Object();
    cj_add_identity(o, rF.mId, rF.mName);
    o.Add("swapped_surface_normal", JsonValue::Boolean(rF.mSwappedNormal));
    JsonValue s = JsonValue::Object();
    s.Add("is_trimmed", JsonValue::Boolean(rF.mTrimmed));
    s.Add("is_rational", JsonValue::Boolean(rF.mRational));
    JsonValue degs = JsonValue::Array();
    degs.Push(JsonValue::Integer(patch.mDegrees[0]));
    degs.Push(JsonValue::Integer(patch.mDegrees[1]));
    s.Add("degrees", std::move(degs));
    JsonValue kvs = JsonValue::Array();
    kvs.Push(cj_numbers(cj_written_knots(patch.mKnots[0], rF.mReducedKnots)));
    kvs.Push(cj_numbers(cj_written_knots(patch.mKnots[1], rF.mReducedKnots)));
    s.Add("knot_vectors", std::move(kvs));
    s.Add("control_points",
          cj_net(patch.mControlPoints, patch.mWeights, rF.mControlPointIds, rWhat + " surface"));
    cj_add_extra(s, rF.mSurfaceExtra, rWhat + " surface");
    o.Add("surface", std::move(s));
    JsonValue loops = JsonValue::Array();
    for (std::size_t l = 0; l < rF.mLoops.size(); ++l) {
        const BrepLoop& loop = rF.mLoops[l];
        JsonValue lo = JsonValue::Object();
        lo.Add("loop_type", JsonValue::String(loop.mType));
        JsonValue trims = JsonValue::Array();
        for (const BrepTrim& t : loop.mTrims)
            trims.Push(cj_trim_json(t, rWhat + " loop " + std::to_string(l)));
        lo.Add("trimming_curves", std::move(trims));
        cj_add_extra(lo, loop.mExtra, rWhat + " loop");
        loops.Push(std::move(lo));
    }
    o.Add("boundary_loops", std::move(loops));
    JsonValue emb = JsonValue::Array();
    for (const BrepTrim& t : rF.mEmbeddedEdges)
        emb.Push(cj_trim_json(t, rWhat + " embedded edge"));
    o.Add("embedded_edges", std::move(emb));
    cj_add_extra(o, rF.mExtra, rWhat);
    return o;
}

std::string cj_render(const BrepInfo& rInfo) {
    // Kratos makes model-part nodes of surface and 3-D curve control points,
    // which need ids; trimming curves' are points and may go without.
    std::size_t bare = 0;
    for (const auto& f : rInfo.mFaces)
        bare += f.mControlPointIds.Size() == 0 ? 1 : 0;
    for (const auto& e : rInfo.mEdges)
        bare += e.mHasCurve && e.mCurve.mIds.Size() == 0 ? 1 : 0;
    if (bare)
        log::warn("cad_json: {} surface or 3-D curve control nets have no ids; Kratos's "
                  "CadJsonInput needs [id, [x, y, z, w]] for them",
                  bare);
    JsonValue root = JsonValue::Object();
    cj_add_extra(root, rInfo.mExtra, "the file");
    std::vector<Brep> breps = rInfo.mBreps;
    if (breps.empty()) {
        // A B-rep built in memory with no grouping: one entry holding all.
        Brep all;
        std::int64_t next = 0;
        for (const auto& f : rInfo.mFaces)
            next = std::max(next, f.mId);
        for (const auto& e : rInfo.mEdges)
            next = std::max(next, e.mId);
        for (const auto& v : rInfo.mVertices)
            next = std::max(next, v.mId);
        all.mId = next + 1;
        for (std::size_t i = 0; i < rInfo.mFaces.size(); ++i)
            all.mFaces.push_back(static_cast<std::int64_t>(i));
        for (std::size_t i = 0; i < rInfo.mEdges.size(); ++i)
            all.mEdges.push_back(static_cast<std::int64_t>(i));
        for (std::size_t i = 0; i < rInfo.mVertices.size(); ++i)
            all.mVertices.push_back(static_cast<std::int64_t>(i));
        breps.push_back(std::move(all));
    }
    JsonValue list = JsonValue::Array();
    for (std::size_t b = 0; b < breps.size(); ++b) {
        const Brep& brep = breps[b];
        const std::string bw = "brep " + std::to_string(b);
        JsonValue bo = JsonValue::Object();
        cj_add_identity(bo, brep.mId, brep.mName);
        JsonValue faces = JsonValue::Array();
        for (std::int64_t fi : brep.mFaces) {
            if (fi < 0 || static_cast<std::size_t>(fi) >= rInfo.mFaces.size())
                cj_wfail(bw + " names face " + std::to_string(fi) + ", which does not exist");
            faces.Push(cj_face_json(rInfo, rInfo.mFaces[static_cast<std::size_t>(fi)],
                                    "face " + std::to_string(fi)));
        }
        bo.Add("faces", std::move(faces));
        JsonValue edges = JsonValue::Array();
        for (std::int64_t ei : brep.mEdges) {
            if (ei < 0 || static_cast<std::size_t>(ei) >= rInfo.mEdges.size())
                cj_wfail(bw + " names edge " + std::to_string(ei) + ", which does not exist");
            const BrepEdge& e = rInfo.mEdges[static_cast<std::size_t>(ei)];
            const std::string ew = "edge " + std::to_string(ei);
            JsonValue eo = JsonValue::Object();
            cj_add_identity(eo, e.mId, e.mName);
            if (e.mHasCurve)
                eo.Add("3d_curve", cj_curve_json(e.mCurve, ew + " 3d_curve"));
            eo.Add("topology", cj_topology_json(e.mTopology, ew));
            cj_add_extra(eo, e.mExtra, ew);
            edges.Push(std::move(eo));
        }
        bo.Add("edges", std::move(edges));
        JsonValue verts = JsonValue::Array();
        for (std::int64_t vi : brep.mVertices) {
            if (vi < 0 || static_cast<std::size_t>(vi) >= rInfo.mVertices.size())
                cj_wfail(bw + " names vertex " + std::to_string(vi) + ", which does not exist");
            const BrepVertex& v = rInfo.mVertices[static_cast<std::size_t>(vi)];
            const std::string vw = "vertex " + std::to_string(vi);
            JsonValue vo = JsonValue::Object();
            cj_add_identity(vo, v.mId, v.mName);
            vo.Add("topology", cj_topology_json(v.mTopology, vw));
            cj_add_extra(vo, v.mExtra, vw);
            verts.Push(std::move(vo));
        }
        bo.Add("vertices", std::move(verts));
        cj_add_extra(bo, brep.mExtra, bw);
        list.Push(std::move(bo));
    }
    root.Add("breps", std::move(list));
    return detail::write_json(root, 4) + "\n";
}

// A plain mesh as a B-rep: quads as bilinear faces, lines as edges.
BrepInfo cj_from_mesh(const Mesh& rMesh) {
    BrepInfo info;
    info.mExtra = "{\"version_number\":1}";
    const NDArray& pts = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t npts = rMesh.NumPoints();
    auto coord = [&](std::int64_t P, std::size_t D) -> double {
        if (D >= dim)
            return 0.0;
        const std::size_t i = static_cast<std::size_t>(P) * dim + D;
        switch (pts.Dtype()) {
            case DType::Float32:
                return static_cast<double>(pts.As<float>()[i]);
            case DType::Float64:
                return pts.As<double>()[i];
            default:
                cj_wfail("points are not floating point");
        }
    };
    std::int64_t next_id = 1;
    std::size_t dropped = 0;
    Brep brep;
    for (const auto cb : rMesh.CellRange()) {
        const std::string type = cb.Type();
        const bool quad = type == "quad", line = type == "line";
        if (!(quad || line) || cb.IsRagged()) {
            dropped += cb.NumCells();
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::size_t k = cb.NodesPerCell();
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            std::vector<std::int64_t> nodes(k);
            for (std::size_t j = 0; j < k; ++j) {
                const std::int64_t v = conn.Dtype() == DType::Int32
                                           ? conn.As<std::int32_t>()[c * k + j]
                                           : conn.As<std::int64_t>()[c * k + j];
                if (v < 0 || static_cast<std::size_t>(v) >= npts)
                    cj_wfail("a cell names point " + std::to_string(v) + ", which does not exist");
                nodes[j] = v;
            }
            // The net runs u fastest: corners 0, 1, 3, 2 of the quad.
            const std::vector<std::int64_t> order =
                quad ? std::vector<std::int64_t>{nodes[0], nodes[1], nodes[3], nodes[2]} : nodes;
            NDArray xyz(DType::Float64, {order.size(), 3});
            NDArray w(DType::Float64, {order.size()});
            NDArray ids(DType::Int64, {order.size()});
            for (std::size_t j = 0; j < order.size(); ++j) {
                for (std::size_t d = 0; d < 3; ++d)
                    xyz.As<double>()[j * 3 + d] = coord(order[j], d);
                w.As<double>()[j] = 1.0;
                ids.As<std::int64_t>()[j] = order[j] + 1;
            }
            if (quad) {
                NurbsPatch patch;
                patch.mDegrees = {1, 1};
                patch.mKnots = {{0.0, 0.0, 1.0, 1.0}, {0.0, 0.0, 1.0, 1.0}};
                patch.mControlPoints = std::move(xyz);
                patch.mWeights = std::move(w);
                BrepFace f;
                f.mId = next_id++;
                f.mPatch = static_cast<std::int64_t>(info.mPatches.size());
                f.mControlPointIds = std::move(ids);
                info.mPatches.push_back(std::move(patch));
                brep.mFaces.push_back(static_cast<std::int64_t>(info.mFaces.size()));
                info.mFaces.push_back(std::move(f));
            } else {
                BrepEdge e;
                e.mId = next_id++;
                e.mHasCurve = true;
                e.mCurve.mDegree = 1;
                e.mCurve.mKnots = {0.0, 0.0, 1.0, 1.0};
                e.mCurve.mControlPoints = std::move(xyz);
                e.mCurve.mWeights = std::move(w);
                e.mCurve.mIds = std::move(ids);
                brep.mEdges.push_back(static_cast<std::int64_t>(info.mEdges.size()));
                info.mEdges.push_back(std::move(e));
            }
        }
    }
    if (dropped)
        log::warn("cad_json: {} cells that are neither quad nor line were not written", dropped);
    if (info.mFaces.empty() && info.mEdges.empty())
        cj_wfail("the mesh has no quad or line cells to write as faces or edges");
    brep.mId = next_id;
    info.mBreps.push_back(std::move(brep));
    return info;
}

void cj_write_text(const std::string& rPath, const std::string& rText) {
    auto out = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!out)
        throw WriteError("meshio++: cad_json: cannot open " + rPath + " for writing");
    out.write(rText.data(), static_cast<std::streamsize>(rText.size()));
    if (!out)
        throw WriteError("meshio++: cad_json: failed writing " + rPath);
}

}  // namespace

Mesh read_cad_json(const std::string& rPath) {
    BrepInfo info;
    return read_cad_json(rPath, info);
}

Mesh read_cad_json(const std::string& rPath, BrepInfo& rInfo, const CadJsonReadOptions& rOptions) {
    cj_check_options(rOptions);
    BrepInfo info = cj_parse(rPath);
    Mesh mesh;
    try {
        mesh = cj_sample(info, rOptions);
    } catch (const std::invalid_argument& e) {
        // The parse checks what sampling needs; this is a file it let through.
        throw ReadError(std::string("meshio++: cad_json: ") + e.what());
    }
    rInfo = std::move(info);
    return mesh;
}

BrepInfo read_cad_json_info(const std::string& rPath) {
    return cj_parse(rPath);
}

Mesh sample_brep(const BrepInfo& rInfo, const CadJsonReadOptions& rOptions) {
    return cj_sample(rInfo, rOptions);
}

void write_cad_json(const std::string& rPath, const Mesh& rMesh) {
    cj_write_text(rPath, cj_render(cj_from_mesh(rMesh)));
}

void write_cad_json(const std::string& rPath, const Mesh& rMesh, const BrepInfo& rInfo) {
    if (rInfo.mFaces.empty() && rInfo.mEdges.empty() && rInfo.mVertices.empty()) {
        write_cad_json(rPath, rMesh);
        return;
    }
    cj_write_text(rPath, cj_render(rInfo));
}

}  // namespace meshioplusplus
