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

// Project includes
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/exceptions.hpp"

#ifdef MESHIOPLUSPLUS_HAS_OCCT
// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
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

// External includes (OpenCASCADE)
#include <APIHeaderSection_MakeHeader.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepGProp.hxx>
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
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomLProp_SLProps.hxx>
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
#include <ShapeFix_Solid.hxx>
#include <Standard_Failure.hxx>
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

#ifndef MESHIOPLUSPLUS_HAS_OCCT

FitCadResult fit_cad(const Mesh&, const std::string&, const FitCadOptions&) {
    throw Unsupported(
        "fit_cad: this build has no OpenCASCADE backend; rebuild with "
        "-DMESHIOPLUSPLUS_WITH_OCCT=ON (OpenCASCADE 7.6 or later)");
}

#else

namespace {

constexpr const char* kFcPrefix = "fit_cad: ";
constexpr std::size_t kFcParams = 10;
constexpr double kFcPi = 3.14159265358979323846;

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
bool fc_project(const FcHandle<Geom_Surface>& rS, const gp_Pnt& rP, gp_Pnt& rOut,
                double* pU = nullptr, double* pV = nullptr) {
    GeomAPI_ProjectPointOnSurf proj(rP, rS);
    if (!proj.IsDone() || proj.NbPoints() == 0)
        return false;
    rOut = proj.NearestPoint();
    if (pU && pV)
        proj.LowerDistanceParameters(*pU, *pV);
    return true;
}

// Move @p rP onto every surface of @p rSurfaces at once, by averaged
// alternating projection; returns the largest remaining distance to any of
// them, or +inf when a projection fails.
double fc_onto_all(const std::vector<FcHandle<Geom_Surface> >& rSurfaces, gp_Pnt& rP,
                   double Scale) {
    const double inf = std::numeric_limits<double>::infinity();
    if (rSurfaces.empty())
        return 0.0;
    for (int it = 0; it < 200; ++it) {
        gp_XYZ sum(0.0, 0.0, 0.0);
        for (const auto& s : rSurfaces) {
            gp_Pnt q;
            if (!fc_project(s, rP, q))
                return inf;
            sum += q.XYZ();
        }
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
    auto near = [](double X, double Ref, double Per) {
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
                    u = near(u, pu, uper);
                    v = near(v, pv, vper);
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

// What one build round produced.
struct FcBuild {
    TopoDS_Shape mShape;
    std::vector<std::int64_t> mFaceGroup;  // per final face (explorer order): its group
    std::vector<TopoDS_Face> mFaces;       // final faces, explorer order
    std::vector<std::uint8_t> mBadRow;     // per patch row: demote it
    std::int64_t mNumExact = 0;
    std::int64_t mNumApprox = 0;
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
FcBuild fc_build(const FcSkin& rSkin, const std::vector<FcHandle<Geom_Surface> >& rSurface,
                 const std::vector<std::uint8_t>& rAnalytic, double Tol, bool ExactEdges) {
    const detail::ClosedFacets& fc = rSkin.mFacets;
    const std::vector<double>& xyz = fc.mXyz;
    const std::size_t nf = fc.NumFacets();
    const std::int64_t np = static_cast<std::int64_t>(rSkin.mPatches.size());
    const std::size_t npts = xyz.size() / 3;
    FcBuild out;
    out.mBadRow.assign(rSkin.mPatches.size(), 0);

    auto group_of = [&](std::size_t F) -> std::int64_t {
        const std::int64_t r = rSkin.mRow[F];
        return r >= 0 && rAnalytic[static_cast<std::size_t>(r)] ? r
                                                                : np + static_cast<std::int64_t>(F);
    };
    auto is_analytic = [&](std::int64_t G) { return G >= 0 && G < np; };
    auto mark_bad = [&](std::int64_t G, const char* pWhy) {
        if (!is_analytic(G))
            return;
        if (!out.mBadRow[static_cast<std::size_t>(G)])
            log::debug("fit_cad: demoting patch {}: {}",
                       rSkin.mPatches[static_cast<std::size_t>(G)].mLabel, pWhy);
        out.mBadRow[static_cast<std::size_t>(G)] = 1;
    };

    // The groups around every point, and whether a facet group is among them.
    std::vector<std::pair<std::int64_t, std::int64_t> > vg;
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
    std::vector<std::pair<std::int64_t, std::int64_t> > bedges;
    for (const auto& kv : rSkin.mEdgeFacets) {
        const std::int64_t ga = group_of(static_cast<std::size_t>(kv.second[0]));
        const std::int64_t gb = group_of(static_cast<std::size_t>(kv.second[1]));
        if (ga != gb)
            bedges.push_back(kv.first);
    }
    std::sort(bedges.begin(), bedges.end());
    std::vector<std::vector<std::int64_t> > nbr(npts);
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
                std::vector<FcHandle<Geom_Surface> > surfs;
                for (std::int64_t g : groups_at(v))
                    surfs.push_back(rSurface[static_cast<std::size_t>(g)]);
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

    // Edges.
    for (FcChain& c : chains) {
        const TopoDS_Vertex& va = vertex.at(c.mVerts.front());
        const TopoDS_Vertex& vb = vertex.at(c.mVerts.back());
        try {
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

            if (ExactEdges) {
                GeomAPI_IntSS inter(s1, s2, Precision::Confusion());
                FcHandle<Geom_Curve> best;
                double best_dev = std::numeric_limits<double>::infinity();
                if (inter.IsDone())
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
            const std::vector<FcHandle<Geom_Surface> > both{s1, s2};
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
                const double d = fc_onto_all(both, on, rSkin.mScale);
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
            mark_bad(c.mGroupFwd, "an edge on its boundary failed");
            mark_bad(c.mGroupBwd, "an edge on its boundary failed");
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
    std::map<std::int64_t, std::vector<std::pair<std::size_t, bool> > >
        uses;  // group -> (chain, fwd)
    for (std::size_t i = 0; i < chains.size(); ++i) {
        uses[chains[i].mGroupFwd].emplace_back(i, true);
        uses[chains[i].mGroupBwd].emplace_back(i, false);
    }
    std::vector<std::pair<TopoDS_Face, std::int64_t> > faces;  // with their group
    for (const auto& kv : uses) {
        const std::int64_t g = kv.first;
        bool broken = false;
        for (const auto& u : kv.second)
            broken = broken || !chains[u.first].mOk;
        if (broken) {
            mark_bad(g, "an edge on its boundary failed");
            continue;
        }
        // start vertex -> unused oriented chains
        std::map<std::int64_t, std::vector<std::size_t> > from;
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
            mark_bad(g, "its boundary does not close into loops");
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
                fix.Perform();
                faces.emplace_back(fix.Face(), g);
            } else {
                const std::size_t facet = static_cast<std::size_t>(g - np);
                const gp_Pnt p0 =
                    fc_point(xyz, fc.mNodes[static_cast<std::size_t>(fc.mStart[facet])]);
                const gp_Dir n(fc.mNormals[facet * 3], fc.mNormals[facet * 3 + 1],
                               fc.mNormals[facet * 3 + 2]);
                BRepBuilderAPI_MakeFace mf(gp_Pln(p0, n), wires.front(), Standard_True);
                if (!mf.IsDone())
                    throw Standard_Failure("no planar facet face");
                faces.emplace_back(mf.Face(), g);
            }
        } catch (const Standard_Failure&) {
            mark_bad(g, "OpenCASCADE could not make its face");
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
                mark_bad(g, "it has no boundary, and its surface is not closed");
                continue;
            }
            BRepBuilderAPI_MakeFace mf(rSurface[r], Precision::Confusion());
            if (!mf.IsDone()) {
                mark_bad(g, "OpenCASCADE could not make its face");
                continue;
            }
            faces.emplace_back(mf.Face(), g);
        }
    }
    for (std::uint8_t b : out.mBadRow)
        if (b)
            return out;

    // Sew the faces (their edges are shared already; sewing closes what the
    // face repair re-made) and close each shell into a solid.
    BRepBuilderAPI_Sewing sew(Tol);
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
                    mark_bad(g, "sewing left a free or multiple edge on it");
                    any = true;
                }
            }
        };
        for (int i = 1; i <= sew.NbFreeEdges(); ++i)
            blame(sew.FreeEdge(i));
        for (int i = 1; i <= sew.NbMultipleEdges(); ++i)
            blame(sew.MultipleEdge(i));
        if (!any)
            for (std::int64_t r = 0; r < np; ++r)
                if (rAnalytic[static_cast<std::size_t>(r)])
                    mark_bad(r, "sewing left free edges among facets");
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

    // Per face, the group: by the face's own shape, else by the nearest facet.
    for (TopExp_Explorer ex(out.mShape, TopAbs_FACE); ex.More(); ex.Next()) {
        out.mFaces.push_back(TopoDS::Face(ex.Current()));
        out.mFaceGroup.push_back(face_group(ex.Current()));
    }

    BRepCheck_Analyzer check(out.mShape);
    out.mValid = check.IsValid() && !solids.empty();
    if (!out.mValid) {
        bool any = false;
        for (std::size_t i = 0; i < out.mFaces.size(); ++i)
            if (is_analytic(out.mFaceGroup[i]) && !BRepCheck_Analyzer(out.mFaces[i]).IsValid()) {
                log::debug("fit_cad: face {} fails the B-rep check (shape kind:status){}", i,
                           fc_check_report(out.mFaces[i]));
                mark_bad(out.mFaceGroup[i], "its face fails the B-rep check");
                any = true;
            }
        if (!any)
            for (std::int64_t r = 0; r < np; ++r)
                if (rAnalytic[static_cast<std::size_t>(r)])
                    mark_bad(r, "the solid fails the B-rep check");
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

}  // namespace

FitCadResult fit_cad(const Mesh& rMesh, const std::string& rStepPath,
                     const FitCadOptions& rOptions) {
    if (std::isnan(rOptions.mTolerance) || std::isinf(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kFcPrefix) + "tolerance must be finite");
    if (std::isnan(rOptions.mMeshDeflection) || std::isinf(rOptions.mMeshDeflection))
        throw std::invalid_argument(std::string(kFcPrefix) + "mesh_deflection must be finite");
    if ((rOptions.mTypes & ~kFitAllTypes) != 0)
        throw std::invalid_argument(std::string(kFcPrefix) +
                                    "types may name only plane, sphere, cylinder, cone and torus");
    if (rOptions.mMaxDemotionRounds < 0)
        throw std::invalid_argument(std::string(kFcPrefix) + "max_demotion_rounds must be >= 0");

    detail::LibraryLock lock;
    FcQuietMessenger quiet;

    // The annotated skin: segment and fit first when the labels are missing.
    const Mesh* pWork = &rMesh;
    Mesh owned;
    if (!rMesh.HasCellData("segment:patch")) {
        log::info("fit_cad: no segment:patch labels; running segment and fit_primitives");
        owned = fit_primitives(segment(rMesh).mMesh).mMesh;
        pWork = &owned;
    } else if (!rMesh.HasFieldData("fit:patch_id")) {
        log::info("fit_cad: no fit:patch_* table; running fit_primitives");
        owned = fit_primitives(rMesh).mMesh;
        pWork = &owned;
    }
    const Mesh& work = *pWork;

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
    const std::size_t m = ids.Size();
    if (types.Size() != m || tols.Size() != m || params.Size() != m * kFcParams)
        throw std::invalid_argument(std::string(kFcPrefix) +
                                    "the fit:patch_* arrays disagree on the number of patches");
    std::unordered_map<std::int64_t, std::int64_t> row_of;
    skin.mPatches.resize(m);
    for (std::size_t r = 0; r < m; ++r) {
        FcPatch& p = skin.mPatches[r];
        p.mLabel = detail::read_int(ids, r);
        p.mType = static_cast<std::int32_t>(detail::read_int(types, r));
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
    std::vector<std::vector<std::size_t> > facets_of(m);
    for (std::size_t f = 0; f < nf; ++f)
        if (skin.mRow[f] >= 0)
            facets_of[static_cast<std::size_t>(skin.mRow[f])].push_back(f);
    std::vector<FcHandle<Geom_Surface> > surface(m);
    std::vector<std::uint8_t> analytic(m, 0);
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

    // Build, check, demote and rebuild.
    FcBuild build;
    FcTriangulation tri;
    for (std::int64_t round = 0;; ++round) {
        bool any = false;
        for (std::uint8_t a : analytic)
            any = any || a;
        try {
            build = fc_build(skin, surface, analytic, tol, rOptions.mExactEdges);
        } catch (const Standard_Failure& rErr) {
            build = FcBuild{};
            build.mBadRow.assign(m, 0);
            for (std::size_t r = 0; r < m; ++r)
                build.mBadRow[r] = analytic[r];
            log::debug("fit_cad: OpenCASCADE failed building the solid: {}",
                       rErr.GetMessageString());
        }
        bool bad = false;
        for (std::uint8_t b : build.mBadRow)
            bad = bad || b;
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
                            g < static_cast<std::int64_t>(m))
                            build.mBadRow[static_cast<std::size_t>(g)] = 1;
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
                        if (d > tol && !build.mBadRow[static_cast<std::size_t>(g)]) {
                            log::debug("fit_cad: demoting patch {}: its face lies {} from the skin",
                                       skin.mPatches[static_cast<std::size_t>(g)].mLabel, d);
                            build.mBadRow[static_cast<std::size_t>(g)] = 1;
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
                            if (d > tol && !build.mBadRow[static_cast<std::size_t>(r)]) {
                                log::debug("fit_cad: demoting patch {}: the skin lies {} from it",
                                           skin.mPatches[static_cast<std::size_t>(r)].mLabel, d);
                                build.mBadRow[static_cast<std::size_t>(r)] = 1;
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
        if (round >= rOptions.mMaxDemotionRounds) {
            std::fill(analytic.begin(), analytic.end(), std::uint8_t{0});
        } else {
            for (std::size_t r = 0; r < m; ++r)
                if (build.mBadRow[r])
                    analytic[r] = 0;
        }
    }
    if (build.mShape.IsNull())
        throw std::runtime_error(std::string(kFcPrefix) +
                                 "OpenCASCADE could not build even the faceted solid");

    // The report.
    result.mValid = build.mValid;
    result.mNumPatches = static_cast<std::int64_t>(m);
    result.mNumExactEdges = build.mNumExact;
    result.mNumApproxEdges = build.mNumApprox;
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
            default:
                ++result.mNumFacetedFaces;
                break;
        }
    }

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

    if (!result.mValid)
        log::warn("fit_cad: the solid does not pass OpenCASCADE's B-rep check");
    if (result.mNumDemoted > 0)
        log::info("fit_cad: {} fitted patch(es) demoted to facets", result.mNumDemoted);
    if (!rStepPath.empty())
        fc_write_step(build.mShape, rStepPath);
    return result;
}

#endif  // MESHIOPLUSPLUS_HAS_OCCT

}  // namespace meshioplusplus
