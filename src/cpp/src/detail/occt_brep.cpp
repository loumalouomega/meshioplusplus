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
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRepTools.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomConvert.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IGESControl_Reader.hxx>
#include <Poly_Triangulation.hxx>
#include <STEPControl_Reader.hxx>
#include <Standard_Failure.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/ndarray.hpp"

// Project includes (private, not installed)
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

// The face's surface as a non-periodic B-spline in the world frame.
ObHandle<Geom_BSplineSurface> ob_surface(const TopoDS_Face& rFace) {
    TopLoc_Location loc;
    ObHandle<Geom_Surface> s = BRep_Tool::Surface(rFace, loc);
    if (s.IsNull())
        return {};
    ObHandle<Geom_BSplineSurface> bs = ObHandle<Geom_BSplineSurface>::DownCast(s);
    if (bs.IsNull()) {
        // NurbsConvert leaves nothing else; a surface it could not convert is
        // approximated over the face's bounds.
        double u0 = 0.0, u1 = 0.0, v0 = 0.0, v1 = 0.0;
        BRepTools::UVBounds(rFace, u0, u1, v0, v1);
        bs = GeomConvert::SurfaceToBSplineSurface(
            new Geom_RectangularTrimmedSurface(s, u0, u1, v0, v1));
        if (bs.IsNull())
            return {};
    } else {
        bs = ObHandle<Geom_BSplineSurface>::DownCast(bs->Copy());
    }
    if (bs->IsUPeriodic())
        bs->SetUNotPeriodic();
    if (bs->IsVPeriodic())
        bs->SetVNotPeriodic();
    if (!loc.IsIdentity())
        bs->Transform(loc.Transformation());
    return bs;
}

ObHandle<Geom2d_BSplineCurve> ob_pcurve(const TopoDS_Edge& rEdge, const TopoDS_Face& rFace,
                                        double& rFirst, double& rLast) {
    ObHandle<Geom2d_Curve> c = BRep_Tool::CurveOnSurface(rEdge, rFace, rFirst, rLast);
    if (c.IsNull())
        return {};
    ObHandle<Geom2d_BSplineCurve> bc = ObHandle<Geom2d_BSplineCurve>::DownCast(c);
    if (bc.IsNull()) {
        bc = Geom2dConvert::CurveToBSplineCurve(new Geom2d_TrimmedCurve(c, rFirst, rLast));
        if (bc.IsNull())
            return {};
        rFirst = bc->FirstParameter();
        rLast = bc->LastParameter();
    } else {
        bc = ObHandle<Geom2d_BSplineCurve>::DownCast(bc->Copy());
    }
    if (bc->IsPeriodic())
        bc->SetNotPeriodic();
    return bc;
}

ObHandle<Geom_BSplineCurve> ob_edge_curve(const TopoDS_Edge& rEdge, double& rFirst,
                                          double& rLast) {
    if (BRep_Tool::Degenerated(rEdge))
        return {};
    TopLoc_Location loc;
    ObHandle<Geom_Curve> c = BRep_Tool::Curve(rEdge, loc, rFirst, rLast);
    if (c.IsNull())
        return {};
    ObHandle<Geom_BSplineCurve> bc = ObHandle<Geom_BSplineCurve>::DownCast(c);
    if (bc.IsNull()) {
        bc = GeomConvert::CurveToBSplineCurve(new Geom_TrimmedCurve(c, rFirst, rLast));
        if (bc.IsNull())
            return {};
        rFirst = bc->FirstParameter();
        rLast = bc->LastParameter();
    } else {
        bc = ObHandle<Geom_BSplineCurve>::DownCast(bc->Copy());
    }
    if (bc->IsPeriodic())
        bc->SetNotPeriodic();
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

TopoDS_Shape occt_nurbs_shape(const TopoDS_Shape& rShape) {
    BRepBuilderAPI_NurbsConvert convert(rShape, Standard_True);
    return convert.IsDone() ? convert.Shape() : rShape;
}

BrepInfo occt_brep_info(const TopoDS_Shape& rNurbs) {
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
        const TopoDS_Face face = TopoDS::Face(faces(fi));
        const ObHandle<Geom_BSplineSurface> bs = ob_surface(face);
        if (bs.IsNull())
            continue;
        NurbsPatch patch;
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
        BrepFace f;
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
                const ObHandle<Geom2d_BSplineCurve> pc = ob_pcurve(e, fwd, a, b);
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
                    edge_refs[static_cast<std::size_t>(ei - 1)].push_back(std::move(r));
                }
                loop.mTrims.push_back(std::move(t));
            }
            if (!loop.mTrims.empty())
                f.mLoops.push_back(std::move(loop));
        }
        f.mTrimmed = !f.mLoops.empty();
        info.mPatches.push_back(std::move(patch));
        brep.mFaces.push_back(static_cast<std::int64_t>(info.mFaces.size()));
        info.mFaces.push_back(std::move(f));
    }
    for (int ei = 1; ei <= edges.Extent(); ++ei) {
        const TopoDS_Edge edge = TopoDS::Edge(edges(ei));
        BrepEdge e;
        e.mId = nf + ei;
        double a = 0.0, b = 0.0;
        const ObHandle<Geom_BSplineCurve> c = ob_edge_curve(edge, a, b);
        if (!c.IsNull()) {
            e.mHasCurve = true;
            e.mCurve = ob_curve(c, a, b);
            e.mCurve.mIds = ob_ids(next_id, e.mCurve.mWeights.Size());
        }
        e.mTopology = std::move(edge_refs[static_cast<std::size_t>(ei - 1)]);
        brep.mEdges.push_back(static_cast<std::int64_t>(info.mEdges.size()));
        info.mEdges.push_back(std::move(e));
    }
    info.mBreps.push_back(std::move(brep));
    return info;
}

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
