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

// OpenCASCADE shapes as BrepInfo and Mesh. See detail/occt_brep.hpp.

#ifdef MESHIOPLUSPLUS_HAS_OCCT

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// External includes (OpenCASCADE)
#include <BRepBndLib.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <Geom2dAPI_PointsToBSpline.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <GeomConvert.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom2d_Line.hxx>
#include <Geom_Line.hxx>
#include <Geom_Plane.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IGESControl_Reader.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_Failure.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"

// Project includes (private, not installed)
#include "json_value.hpp"
#include "nurbs_eval.hpp"
#include "occt_brep.hpp"

namespace meshioplusplus::detail {

namespace {

template <class T>
using ObHandle = opencascade::handle<T>;

std::vector<double> ob_knots(const TColStd_Array1OfReal& rSeq) {
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(rSeq.Length()));
    for (int i = rSeq.Lower(); i <= rSeq.Upper(); ++i)
        out.push_back(rSeq.Value(i));
    return out;
}

NDArray ob_ids(std::int64_t& rNext, std::size_t N) {
    NDArray ids(DType::Int64, {N});
    for (std::size_t i = 0; i < N; ++i)
        ids.As<std::int64_t>()[i] = rNext++;
    return ids;
}

// A 3-D or 2-D B-spline curve as a BrepCurve (2-D points as (u, v, 0)).
template <class TCurve>
BrepCurve ob_curve(const ObHandle<TCurve>& rC, double First, double Last) {
    BrepCurve out;
    out.mDegree = rC->Degree();
    out.mRational = rC->IsRational() == Standard_True;
    out.mKnots = ob_knots(rC->KnotSequence());
    const auto n = static_cast<std::size_t>(rC->NbPoles());
    out.mControlPoints = NDArray(DType::Float64, {n, 3});
    out.mWeights = NDArray(DType::Float64, {n});
    double* xyz = out.mControlPoints.As<double>();
    for (std::size_t i = 0; i < n; ++i) {
        const auto p = rC->Pole(static_cast<int>(i) + 1);
        xyz[i * 3] = p.X();
        xyz[i * 3 + 1] = p.Y();
        if constexpr (std::is_same_v<TCurve, Geom_BSplineCurve>)
            xyz[i * 3 + 2] = p.Z();
        else
            xyz[i * 3 + 2] = 0.0;
        out.mWeights.As<double>()[i] = rC->Weight(static_cast<int>(i) + 1);
    }
    out.mActiveRange = {First, Last};
    return out;
}

// A B-spline through samples of At over [A, B] (3-D or 2-D), its degree at
// most 8, in At's own parameter: the fit comes back over [0, 1], and the
// samples being uniform, its knots move onto [A, B]. Sampled more densely
// until it stays within Tol of At halfway between them.
template <class TPoint, class TCurve, class TApprox, class TPoints, class F>
ObHandle<TCurve> ob_approx_curve(F&& rAt, double A, double B, double Tol) {
    for (int n = 9; n <= 1025; n = 2 * n - 1) {
        TPoints pts(1, n);
        TColStd_Array1OfReal params(1, n);
        for (int i = 0; i < n; ++i) {
            const double t = A + (B - A) * i / (n - 1);
            params.SetValue(i + 1, t);
            pts.SetValue(i + 1, rAt(t));
        }
        TApprox fit(pts, params, 3, 8, GeomAbs_C2, 0.1 * Tol);
        if (!fit.IsDone())
            continue;
        ObHandle<TCurve> c = fit.Curve();
        TColStd_Array1OfReal knots(1, c->NbKnots());
        c->Knots(knots);
        const double k0 = knots(1), k1 = knots(knots.Upper());
        for (int k = 1; k <= knots.Upper(); ++k)
            knots(k) = A + (B - A) * (knots(k) - k0) / (k1 - k0);
        c->SetKnots(knots);
        double worst = 0.0;
        for (int i = 0; i + 1 < n; ++i) {
            const double t = A + (B - A) * (i + 0.5) / (n - 1);
            worst = std::max(worst, c->Value(t).Distance(rAt(t)));
        }
        if (worst <= Tol)
            return c;
    }
    return {};
}

// A non-rational B-spline over [U0, U1] x [V0, V1] in the surface's own
// parameters, within Tol of it: a grid fitted iso-parametrically, its knots
// then moved onto those ranges (the grid is uniform in both), sampled more
// densely until the cell centres stay within Tol.
ObHandle<Geom_BSplineSurface> ob_approx_surface(const ObHandle<Geom_Surface>& rS, double U0,
                                                double U1, double V0, double V1, double Tol) {
    for (int n = 9; n <= 129; n = 2 * n - 1) {
        TColgp_Array2OfPnt grid(1, n, 1, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                grid.SetValue(
                    i + 1, j + 1,
                    rS->Value(U0 + (U1 - U0) * i / (n - 1), V0 + (V1 - V0) * j / (n - 1)));
        GeomAPI_PointsToBSplineSurface fit;
        fit.Init(grid, Approx_IsoParametric, 3, 8, GeomAbs_C2, 0.1 * Tol);
        if (!fit.IsDone())
            continue;
        ObHandle<Geom_BSplineSurface> bs = fit.Surface();
        auto remap = [](TColStd_Array1OfReal& rKnots, double Lo, double Hi) {
            const double a = rKnots(rKnots.Lower()), b = rKnots(rKnots.Upper());
            for (int k = rKnots.Lower(); k <= rKnots.Upper(); ++k)
                rKnots(k) = Lo + (Hi - Lo) * (rKnots(k) - a) / (b - a);
        };
        TColStd_Array1OfReal uk(1, bs->NbUKnots()), vk(1, bs->NbVKnots());
        bs->UKnots(uk);
        bs->VKnots(vk);
        remap(uk, U0, U1);
        remap(vk, V0, V1);
        bs->SetUKnots(uk);
        bs->SetVKnots(vk);
        double worst = 0.0;
        for (int i = 0; i + 1 < n; ++i)
            for (int j = 0; j + 1 < n; ++j) {
                const double u = U0 + (U1 - U0) * (i + 0.5) / (n - 1);
                const double v = V0 + (V1 - V0) * (j + 0.5) / (n - 1);
                worst = std::max(worst, bs->Value(u, v).Distance(rS->Value(u, v)));
            }
        if (worst <= Tol)
            return bs;
    }
    return {};
}

// Making a periodic B-spline non-periodic keeps the period its knots start
// at. When the face's trims (or the edge's range) lie outside it, a whole
// period away, the period is first moved to start where they do: a knot
// inserted there, the origin moved onto it and the knots shifted by whole
// periods, which leaves the geometry as it is.
int ob_knot_at(const TColStd_Array1OfReal& rKnots, double U, double Tol) {
    for (int k = rKnots.Lower(); k <= rKnots.Upper(); ++k)
        if (std::fabs(rKnots(k) - U) <= Tol)
            return k;
    return 0;
}

void ob_surface_origin(Geom_BSplineSurface& rS, bool U, double Lo, double Hi) {
    const double period = U ? rS.UPeriod() : rS.VPeriod();
    const double tol = 1e-9 * period;
    TColStd_Array1OfReal k(1, U ? rS.NbUKnots() : rS.NbVKnots());
    U ? rS.UKnots(k) : rS.VKnots(k);
    if (Lo >= k(1) - tol && Hi <= k(1) + period + tol)
        return;  // already inside
    const double at = Lo - period * std::floor((Lo - k(1)) / period);
    if (ob_knot_at(k, at, tol) == 0) {
        if (U)
            rS.InsertUKnot(at, 1, tol, Standard_False);
        else
            rS.InsertVKnot(at, 1, tol, Standard_False);
    }
    TColStd_Array1OfReal grown(1, U ? rS.NbUKnots() : rS.NbVKnots());
    U ? rS.UKnots(grown) : rS.VKnots(grown);
    const int index = ob_knot_at(grown, at, tol);
    if (index > 1) {
        if (U)
            rS.SetUOrigin(index);
        else
            rS.SetVOrigin(index);
    }
    TColStd_Array1OfReal moved(1, U ? rS.NbUKnots() : rS.NbVKnots());
    U ? rS.UKnots(moved) : rS.VKnots(moved);
    const double shift = Lo - moved(1);
    for (int i = moved.Lower(); i <= moved.Upper(); ++i)
        moved(i) += shift;
    if (U)
        rS.SetUKnots(moved);
    else
        rS.SetVKnots(moved);
}

template <class TCurve>
void ob_curve_origin(TCurve& rC, double Lo, double Hi) {
    const double period = rC.Period();
    const double tol = 1e-9 * period;
    TColStd_Array1OfReal k(1, rC.NbKnots());
    rC.Knots(k);
    if (Lo >= k(1) - tol && Hi <= k(1) + period + tol)
        return;
    const double at = Lo - period * std::floor((Lo - k(1)) / period);
    if (ob_knot_at(k, at, tol) == 0)
        rC.InsertKnot(at, 1, tol);
    TColStd_Array1OfReal grown(1, rC.NbKnots());
    rC.Knots(grown);
    const int index = ob_knot_at(grown, at, tol);
    if (index > 1)
        rC.SetOrigin(index);
    TColStd_Array1OfReal moved(1, rC.NbKnots());
    rC.Knots(moved);
    const double shift = Lo - moved(1);
    for (int i = moved.Lower(); i <= moved.Upper(); ++i)
        moved(i) += shift;
    rC.SetKnots(moved);
}

// The face's surface as a non-periodic B-spline in the world frame. A surface
// that is not one (NurbsConvert leaves nothing else) is, with ApproxTol > 0,
// approximated in its own parameters, so the face's trimming curves stay
// valid; without, converted over the face's bounds.
ObHandle<Geom_BSplineSurface> ob_surface(const TopoDS_Face& rFace, double ApproxTol = 0.0) {
    TopLoc_Location loc;
    ObHandle<Geom_Surface> s = BRep_Tool::Surface(rFace, loc);
    if (s.IsNull())
        return {};
    if (ApproxTol > 0.0)
        for (ObHandle<Geom_RectangularTrimmedSurface> t =
                 ObHandle<Geom_RectangularTrimmedSurface>::DownCast(s);
             !t.IsNull(); t = ObHandle<Geom_RectangularTrimmedSurface>::DownCast(s))
            s = t->BasisSurface();
    ObHandle<Geom_BSplineSurface> bs = ObHandle<Geom_BSplineSurface>::DownCast(s);
    if (bs.IsNull()) {
        double u0 = 0.0, u1 = 0.0, v0 = 0.0, v1 = 0.0;
        BRepTools::UVBounds(rFace, u0, u1, v0, v1);
        if (ApproxTol > 0.0 && ObHandle<Geom_Plane>::DownCast(s).IsNull()) {
            // A margin, so the trimming curves lie inside the parameter range.
            const double du = 0.01 * (u1 - u0), dv = 0.01 * (v1 - v0);
            bs = ob_approx_surface(s, u0 - du, u1 + du, v0 - dv, v1 + dv, ApproxTol);
        } else {
            // A plane converts in its own parameters.
            bs = GeomConvert::SurfaceToBSplineSurface(
                new Geom_RectangularTrimmedSurface(s, u0, u1, v0, v1));
        }
        if (bs.IsNull())
            return {};
    } else {
        bs = ObHandle<Geom_BSplineSurface>::DownCast(bs->Copy());
    }
    if (bs->IsUPeriodic() || bs->IsVPeriodic()) {
        double u0 = 0.0, u1 = 0.0, v0 = 0.0, v1 = 0.0;
        BRepTools::UVBounds(rFace, u0, u1, v0, v1);
        if (bs->IsUPeriodic()) {
            ob_surface_origin(*bs, true, u0, u1);
            bs->SetUNotPeriodic();
        }
        if (bs->IsVPeriodic()) {
            ob_surface_origin(*bs, false, v0, v1);
            bs->SetVNotPeriodic();
        }
    }
    if (!loc.IsIdentity())
        bs->Transform(loc.Transformation());
    return bs;
}

ObHandle<Geom2d_BSplineCurve> ob_pcurve(const TopoDS_Edge& rEdge, const TopoDS_Face& rFace,
                                        double& rFirst, double& rLast, double ApproxTol = 0.0) {
    ObHandle<Geom2d_Curve> c = BRep_Tool::CurveOnSurface(rEdge, rFace, rFirst, rLast);
    if (c.IsNull())
        return {};
    ObHandle<Geom2d_BSplineCurve> bc = ObHandle<Geom2d_BSplineCurve>::DownCast(c);
    if (bc.IsNull() && ApproxTol > 0.0) {
        // In the curve's own parameter: a line exactly, anything else within
        // a tiny fraction of its extent in the face's parameters.
        const double a = rFirst, b = rLast;
        if (!ObHandle<Geom2d_Line>::DownCast(c).IsNull()) {
            TColgp_Array1OfPnt2d poles(1, 2);
            poles.SetValue(1, c->Value(a));
            poles.SetValue(2, c->Value(b));
            TColStd_Array1OfReal knots(1, 2);
            knots.SetValue(1, a);
            knots.SetValue(2, b);
            TColStd_Array1OfInteger mults(1, 2);
            mults.Init(2);
            bc = new Geom2d_BSplineCurve(poles, knots, mults, 1);
        } else {
            const double extent = c->Value(a).Distance(c->Value(0.5 * (a + b))) +
                                  c->Value(0.5 * (a + b)).Distance(c->Value(b));
            bc = ob_approx_curve<gp_Pnt2d, Geom2d_BSplineCurve, Geom2dAPI_PointsToBSpline,
                                 TColgp_Array1OfPnt2d>([&](double T) { return c->Value(T); }, a, b,
                                                       std::max(1e-9, 1e-7 * extent));
        }
        if (bc.IsNull())
            return {};
    } else if (bc.IsNull()) {
        bc = Geom2dConvert::CurveToBSplineCurve(new Geom2d_TrimmedCurve(c, rFirst, rLast));
        if (bc.IsNull())
            return {};
        rFirst = bc->FirstParameter();
        rLast = bc->LastParameter();
    } else {
        bc = ObHandle<Geom2d_BSplineCurve>::DownCast(bc->Copy());
    }
    if (bc->IsPeriodic()) {
        ob_curve_origin(*bc, rFirst, rLast);
        bc->SetNotPeriodic();
    }
    return bc;
}

ObHandle<Geom_BSplineCurve> ob_edge_curve(const TopoDS_Edge& rEdge, double& rFirst, double& rLast,
                                          double ApproxTol = 0.0) {
    if (BRep_Tool::Degenerated(rEdge))
        return {};
    TopLoc_Location loc;
    ObHandle<Geom_Curve> c = BRep_Tool::Curve(rEdge, loc, rFirst, rLast);
    if (c.IsNull())
        return {};
    ObHandle<Geom_BSplineCurve> bc = ObHandle<Geom_BSplineCurve>::DownCast(c);
    if (bc.IsNull() && ApproxTol > 0.0) {
        // In the edge's own parameter, as its trims are: a line exactly.
        const double a = rFirst, b = rLast;
        if (!ObHandle<Geom_Line>::DownCast(c).IsNull()) {
            TColgp_Array1OfPnt poles(1, 2);
            poles.SetValue(1, c->Value(a));
            poles.SetValue(2, c->Value(b));
            TColStd_Array1OfReal knots(1, 2);
            knots.SetValue(1, a);
            knots.SetValue(2, b);
            TColStd_Array1OfInteger mults(1, 2);
            mults.Init(2);
            bc = new Geom_BSplineCurve(poles, knots, mults, 1);
        } else {
            bc = ob_approx_curve<gp_Pnt, Geom_BSplineCurve, GeomAPI_PointsToBSpline,
                                 TColgp_Array1OfPnt>([&](double T) { return c->Value(T); }, a, b,
                                                     ApproxTol);
        }
        if (bc.IsNull())
            return {};
    } else if (bc.IsNull()) {
        bc = GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(c, rFirst, rLast));
        if (bc.IsNull())
            return {};
        rFirst = bc->FirstParameter();
        rLast = bc->LastParameter();
    } else {
        bc = ObHandle<Geom_BSplineCurve>::DownCast(bc->Copy());
    }
    if (bc->IsPeriodic()) {
        ob_curve_origin(*bc, rFirst, rLast);
        bc->SetNotPeriodic();
    }
    if (!loc.IsIdentity())
        bc->Transform(loc.Transformation());
    return bc;
}

// The ordinal, among the spans of nonzero length, of the span holding U.
std::int64_t ob_span_ordinal(int P, const std::vector<double>& rKnots, double U) {
    const std::size_t n = rKnots.size() - static_cast<std::size_t>(P) - 1;
    const std::size_t span = nurbs_find_span(P, rKnots, n, U);
    const std::vector<std::size_t> spans = nurbs_spans(P, rKnots, n);
    return static_cast<std::int64_t>(std::lower_bound(spans.begin(), spans.end(), span) -
                                     spans.begin());
}

}  // namespace

struct OcctQuiet::Impl {
    ObHandle<Message_Messenger> mMessenger;
    Message_SequenceOfPrinters mSaved;
};

OcctQuiet::OcctQuiet() : mpImpl(new Impl) {
    mpImpl->mMessenger = Message::DefaultMessenger();
    if (!mpImpl->mMessenger.IsNull()) {
        mpImpl->mSaved = mpImpl->mMessenger->Printers();
        mpImpl->mMessenger->ChangePrinters().Clear();
    }
}

OcctQuiet::~OcctQuiet() {
    if (!mpImpl->mMessenger.IsNull())
        mpImpl->mMessenger->ChangePrinters() = mpImpl->mSaved;
    delete mpImpl;
}

TopoDS_Shape occt_read_cad_file(const std::string& rPath, bool Iges, const char* pFormat) {
    const std::string prefix = std::string("meshio++: ") + pFormat + ": ";
    try {
        TopoDS_Shape shape;
        if (Iges) {
            IGESControl_Reader reader;
            if (reader.ReadFile(rPath.c_str()) != IFSelect_RetDone)
                throw ReadError(prefix + "OpenCASCADE cannot read '" + rPath + "'");
            reader.TransferRoots();
            shape = reader.OneShape();
        } else {
            STEPControl_Reader reader;
            if (reader.ReadFile(rPath.c_str()) != IFSelect_RetDone)
                throw ReadError(prefix + "OpenCASCADE cannot read '" + rPath + "'");
            reader.TransferRoots();
            shape = reader.OneShape();
        }
        if (shape.IsNull())
            throw ReadError(prefix + "'" + rPath + "' holds no shape OpenCASCADE can transfer");
        return shape;
    } catch (const Standard_Failure& rErr) {
        throw ReadError(prefix + "OpenCASCADE failed reading '" + rPath +
                        "': " + rErr.GetMessageString());
    }
}

namespace {

BrepInfo ob_brep_info(const TopoDS_Shape& rNurbs, double ApproxTol);

// Whether the parameter curves of the face's edges stay in its surface's
// parameter range, along each direction that is not periodic.
bool ob_trims_inside(const TopoDS_Face& rFace, const Geom_BSplineSurface& rS) {
    double u0 = 0.0, u1 = 0.0, v0 = 0.0, v1 = 0.0;
    rS.Bounds(u0, u1, v0, v1);
    const double tu = 1e-7 * (u1 - u0), tv = 1e-7 * (v1 - v0);
    for (TopExp_Explorer ex(rFace, TopAbs_EDGE); ex.More(); ex.Next()) {
        double a = 0.0, b = 0.0;
        const ObHandle<Geom2d_Curve> c =
            BRep_Tool::CurveOnSurface(TopoDS::Edge(ex.Current()), rFace, a, b);
        if (c.IsNull())
            return false;
        for (int i = 0; i <= 8; ++i) {
            const gp_Pnt2d p = c->Value(a + (b - a) * i / 8.0);
            if ((!rS.IsUPeriodic() && (p.X() < u0 - tu || p.X() > u1 + tu)) ||
                (!rS.IsVPeriodic() && (p.Y() < v0 - tv || p.Y() > v1 + tv)))
                return false;
        }
    }
    return true;
}

}  // namespace

OcctNurbs occt_nurbs(const TopoDS_Shape& rShape) {
    OcctNurbs out;
    try {
        // Convert a copy with its own geometry: OpenCASCADE 7.9 crashes
        // converting the patch network's solid in place (in
        // BRepTools_NurbsConvertModification), not a copy of it.
        BRepBuilderAPI_Copy copy(rShape, Standard_True);
        BRepBuilderAPI_NurbsConvert convert(copy.Shape(), Standard_True);
        if (convert.IsDone()) {
            // Complete only when every face keeps a surface, a B-spline, whose
            // parameter range holds the face's trims (along a direction that
            // is not periodic: a periodic one is moved onto them).
            int before = 0, kept = 0;
            for (TopExp_Explorer ex(rShape, TopAbs_FACE); ex.More(); ex.Next())
                ++before;
            for (TopExp_Explorer ex(convert.Shape(), TopAbs_FACE); ex.More(); ex.Next()) {
                const TopoDS_Face face = TopoDS::Face(ex.Current());
                TopLoc_Location loc;
                const ObHandle<Geom_BSplineSurface> bs =
                    ObHandle<Geom_BSplineSurface>::DownCast(BRep_Tool::Surface(face, loc));
                if (!bs.IsNull() && ob_trims_inside(face, *bs))
                    ++kept;
            }
            if (kept == before) {
                out.mShape = convert.Shape();
                out.mInfo = ob_brep_info(out.mShape, 0.0);
                out.mExact = true;
                return out;
            }
            log::debug(
                "occt: NURBS conversion left {} of {} face(s) without a B-spline; "
                "approximating instead",
                before - kept, before);
        }
    } catch (const Standard_Failure& rErr) {
        log::debug("occt: NURBS conversion failed ({}); approximating instead",
                   rErr.GetMessageString());
    }
    Bnd_Box box;
    BRepBndLib::Add(rShape, box);
    const double diagonal = box.IsVoid() ? 1.0 : std::sqrt(box.SquareExtent());
    out.mShape = rShape;
    out.mInfo = ob_brep_info(rShape, 1e-6 * diagonal);
    return out;
}

BrepInfo occt_brep_info(const TopoDS_Shape& rNurbs) {
    return ob_brep_info(rNurbs, 0.0);
}

namespace {

BrepInfo ob_brep_info(const TopoDS_Shape& rNurbs, double ApproxTol) {
    BrepInfo info;
    info.mExtra = "{\"version_number\":1}";
    TopTools_IndexedMapOfShape faces, edges;
    TopExp::MapShapes(rNurbs, TopAbs_FACE, faces);
    TopExp::MapShapes(rNurbs, TopAbs_EDGE, edges);
    const auto nf = static_cast<std::int64_t>(faces.Extent());
    const auto ne = static_cast<std::int64_t>(edges.Extent());
    std::vector<std::vector<BrepRef> > edge_refs(static_cast<std::size_t>(ne));
    std::int64_t next_id = 1;
    Brep brep;
    brep.mId = nf + ne + 1;
    for (int fi = 1; fi <= faces.Extent(); ++fi) {
        // The face's references, kept only once all of it converts: a face
        // OpenCASCADE throws on (a surface trimmed to an empty range, an
        // approximation that fails) is left out, as one with no surface is.
        std::vector<std::pair<int, BrepRef> > refs;
        NurbsPatch patch;
        BrepFace f;
        try {
            const TopoDS_Face face = TopoDS::Face(faces(fi));
            const ObHandle<Geom_BSplineSurface> bs = ob_surface(face, ApproxTol);
            if (bs.IsNull())
                continue;
            patch.mDegrees = {bs->UDegree(), bs->VDegree()};
            patch.mKnots = {ob_knots(bs->UKnotSequence()), ob_knots(bs->VKnotSequence())};
            const auto nu = static_cast<std::size_t>(bs->NbUPoles());
            const auto nv = static_cast<std::size_t>(bs->NbVPoles());
            patch.mControlPoints = NDArray(DType::Float64, {nu * nv, 3});
            patch.mWeights = NDArray(DType::Float64, {nu * nv});
            for (std::size_t j = 0; j < nv; ++j)
                for (std::size_t i = 0; i < nu; ++i) {
                    const gp_Pnt p = bs->Pole(static_cast<int>(i) + 1, static_cast<int>(j) + 1);
                    const std::size_t k = j * nu + i;
                    patch.mControlPoints.As<double>()[k * 3] = p.X();
                    patch.mControlPoints.As<double>()[k * 3 + 1] = p.Y();
                    patch.mControlPoints.As<double>()[k * 3 + 2] = p.Z();
                    patch.mWeights.As<double>()[k] =
                        bs->Weight(static_cast<int>(i) + 1, static_cast<int>(j) + 1);
                }
            f.mId = fi;
            f.mPatch = static_cast<std::int64_t>(info.mPatches.size());
            f.mSwappedNormal = face.Orientation() == TopAbs_REVERSED;
            f.mRational = bs->IsURational() || bs->IsVRational();
            f.mControlPointIds = ob_ids(next_id, nu * nv);
            // Loops in the surface's own orientation: the outer one first.
            const TopoDS_Face fwd = TopoDS::Face(face.Oriented(TopAbs_FORWARD));
            const TopoDS_Wire outer = BRepTools::OuterWire(fwd);
            std::vector<TopoDS_Wire> wires;
            if (!outer.IsNull())
                wires.push_back(outer);
            for (TopExp_Explorer ex(fwd, TopAbs_WIRE); ex.More(); ex.Next())
                if (!ex.Current().IsSame(outer))
                    wires.push_back(TopoDS::Wire(ex.Current()));
            std::int64_t trim_index = 0;
            for (const TopoDS_Wire& w : wires) {
                BrepLoop loop;
                loop.mType = w.IsSame(outer) ? "outer" : "inner";
                for (BRepTools_WireExplorer we(w, fwd); we.More(); we.Next()) {
                    const TopoDS_Edge& e = we.Current();
                    double a = 0.0, b = 0.0;
                    const ObHandle<Geom2d_BSplineCurve> pc = ob_pcurve(e, fwd, a, b, ApproxTol);
                    if (pc.IsNull())
                        continue;
                    BrepTrim t;
                    t.mTrimIndex = trim_index++;
                    t.mDirection = e.Orientation() != TopAbs_REVERSED;
                    t.mCurve = ob_curve(pc, a, b);
                    const int ei = edges.FindIndex(e);
                    if (ei > 0) {
                        BrepRef r;
                        r.mId = f.mId;
                        r.mTrimIndex = t.mTrimIndex;
                        // Edge and trim share their parameter in OpenCASCADE.
                        r.mRelativeDirection = 1;
                        refs.emplace_back(ei, std::move(r));
                    }
                    loop.mTrims.push_back(std::move(t));
                }
                if (!loop.mTrims.empty())
                    f.mLoops.push_back(std::move(loop));
            }
            f.mTrimmed = !f.mLoops.empty();
        } catch (const Standard_Failure& rErr) {
            log::warn("occt: face {} could not be converted ({}); left out", fi,
                      rErr.GetMessageString());
            continue;
        }
        for (auto& [ei, r] : refs)
            edge_refs[static_cast<std::size_t>(ei - 1)].push_back(std::move(r));
        info.mPatches.push_back(std::move(patch));
        brep.mFaces.push_back(static_cast<std::int64_t>(info.mFaces.size()));
        info.mFaces.push_back(std::move(f));
    }
    for (int ei = 1; ei <= edges.Extent(); ++ei) {
        const TopoDS_Edge edge = TopoDS::Edge(edges(ei));
        BrepEdge e;
        e.mId = nf + ei;
        double a = 0.0, b = 0.0;
        ObHandle<Geom_BSplineCurve> c;
        try {
            c = ob_edge_curve(edge, a, b, ApproxTol);
        } catch (const Standard_Failure& rErr) {
            // Its uses still pair through their trims.
            log::warn("occt: edge {} could not be converted ({}); its curve is left out", ei,
                      rErr.GetMessageString());
        }
        if (!c.IsNull()) {
            e.mHasCurve = true;
            e.mCurve = ob_curve(c, a, b);
            e.mCurve.mIds = ob_ids(next_id, e.mCurve.mWeights.Size());
        }
        e.mTopology = std::move(edge_refs[static_cast<std::size_t>(ei - 1)]);
        // The model's own tolerances (roadmap 6.5.1): the edge's and its two
        // vertices', each with the approximation's error on that path, which
        // mesh_cad_surface's model_tolerances accepts as the gaps it declares.
        {
            JsonValue extra = JsonValue::Object();
            extra.Add("tolerance", JsonValue::Number(BRep_Tool::Tolerance(edge) + ApproxTol));
            TopoDS_Vertex v1, v2;
            TopExp::Vertices(edge, v1, v2);
            JsonValue ends = JsonValue::Array();
            for (const TopoDS_Vertex& r_v : {v1, v2})
                if (!r_v.IsNull())
                    ends.Push(JsonValue::Number(BRep_Tool::Tolerance(r_v) + ApproxTol));
            if (!ends.mItems.empty())
                extra.Add("end_tolerances", std::move(ends));
            e.mExtra = write_json(extra, 0);
        }
        brep.mEdges.push_back(static_cast<std::int64_t>(info.mEdges.size()));
        info.mEdges.push_back(std::move(e));
    }
    info.mBreps.push_back(std::move(brep));
    return info;
}

}  // namespace

Mesh occt_tessellate(const TopoDS_Shape& rNurbs, const BrepInfo& rInfo, double Deflection,
                     double Angular) {
    double deflection = Deflection;
    if (!(deflection > 0.0)) {
        Bnd_Box box;
        BRepBndLib::Add(rNurbs, box);
        double x0 = 0.0, y0 = 0.0, z0 = 0.0, x1 = 0.0, y1 = 0.0, z1 = 0.0;
        if (!box.IsVoid())
            box.Get(x0, y0, z0, x1, y1, z1);
        const double diag = std::sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) +
                                      (z1 - z0) * (z1 - z0));
        deflection = diag > 0.0 ? 1e-3 * diag : 1e-3;
    }
    BRepMesh_IncrementalMesh mesher(rNurbs, deflection, Standard_False, Angular, Standard_False);
    TopTools_IndexedMapOfShape faces;
    TopExp::MapShapes(rNurbs, TopAbs_FACE, faces);
    // Faces by id: occt_brep_info numbers them by the same map.
    std::unordered_map<std::int64_t, std::size_t> face_of_id;
    for (std::size_t i = 0; i < rInfo.mFaces.size(); ++i)
        face_of_id.emplace(rInfo.mFaces[i].mId, i);

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
    std::vector<double> xyz;
    std::vector<std::int64_t> tris, index, ids, span;
    for (int fi = 1; fi <= faces.Extent(); ++fi) {
        const TopoDS_Face face = TopoDS::Face(faces(fi));
        TopLoc_Location loc;
        const ObHandle<Poly_Triangulation> tri = BRep_Tool::Triangulation(face, loc);
        if (tri.IsNull())
            continue;
        auto fit = face_of_id.find(fi);
        const BrepFace* pFace = fit == face_of_id.end() ? nullptr : &rInfo.mFaces[fit->second];
        const NurbsPatch* pPatch =
            pFace ? &rInfo.mPatches[static_cast<std::size_t>(pFace->mPatch)] : nullptr;
        const gp_Trsf trsf = loc.Transformation();
        std::vector<std::int64_t> local(static_cast<std::size_t>(tri->NbNodes()) + 1, -1);
        for (int i = 1; i <= tri->NbNodes(); ++i) {
            const gp_Pnt p = tri->Node(i).Transformed(trsf);
            const double c[3] = {p.X() + 0.0, p.Y() + 0.0, p.Z() + 0.0};
            Key k;
            for (int j = 0; j < 3; ++j)
                std::memcpy(&k.mBits[static_cast<std::size_t>(j)], &c[j], sizeof(double));
            auto ins = id.emplace(k, static_cast<std::int64_t>(xyz.size() / 3));
            if (ins.second)
                xyz.insert(xyz.end(), c, c + 3);
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
            tris.insert(tris.end(), {ia, ib, ic});
            index.push_back(pFace ? static_cast<std::int64_t>(fit->second) : -1);
            ids.push_back(pFace ? pFace->mId : -1);
            std::int64_t su = -1, sv = -1;
            if (pPatch && tri->HasUVNodes()) {
                const gp_Pnt2d uv0 = tri->UVNode(a), uv1 = tri->UVNode(b), uv2 = tri->UVNode(c);
                const double u = (uv0.X() + uv1.X() + uv2.X()) / 3.0;
                const double v = (uv0.Y() + uv1.Y() + uv2.Y()) / 3.0;
                su = ob_span_ordinal(pPatch->mDegrees[0], pPatch->mKnots[0], u);
                sv = ob_span_ordinal(pPatch->mDegrees[1], pPatch->mKnots[1], v);
            }
            span.insert(span.end(), {su, sv});
        }
    }
    Mesh mesh;
    const std::size_t np = xyz.size() / 3, nt = tris.size() / 3;
    NDArray pts(DType::Float64, {np, 3});
    std::copy(xyz.begin(), xyz.end(), pts.As<double>());
    mesh.AssignPoints(std::move(pts));
    if (nt == 0)
        return mesh;
    NDArray conn(DType::Int64, {nt, 3});
    std::copy(tris.begin(), tris.end(), conn.As<std::int64_t>());
    mesh.AddCellBlock("triangle", std::move(conn));
    auto column = [nt](const std::vector<std::int64_t>& rV, std::size_t Width) {
        NDArray a(DType::Int64, Width == 1 ? std::vector<std::size_t>{nt}
                                           : std::vector<std::size_t>{nt, Width});
        std::copy(rV.begin(), rV.end(), a.As<std::int64_t>());
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(a));
        return blocks;
    };
    mesh.AddCellData("cad:id", column(ids, 1));
    mesh.AddCellData("cad:index", column(index, 1));
    mesh.AddCellData("cad:kind", column(std::vector<std::int64_t>(nt, 0), 1));
    mesh.AddCellData("cad:span", column(span, 2));
    mesh.AddCellData("cad:trimmed", column(std::vector<std::int64_t>(nt, 0), 1));
    return mesh;
}

}  // namespace meshioplusplus::detail

#endif  // MESHIOPLUSPLUS_HAS_OCCT
