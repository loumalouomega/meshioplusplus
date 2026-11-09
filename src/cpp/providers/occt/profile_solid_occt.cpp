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
// Profile solids with OpenCASCADE (roadmap 7.7.4): the provider behind
// operations/profile_solid.hpp. Compiled into the core in a monolithic build,
// and into the meshioplusplus-occt add-on in a split one.

#ifdef MESHIOPLUSPLUS_HAS_OCCT

// System includes
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

// External includes (OpenCASCADE)
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepFilletAPI_MakeChamfer.hxx>
#include <BRepFilletAPI_MakeFillet.hxx>
#include <BRepGProp.hxx>
#include <BRepOffsetAPI_MakeThickSolid.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRep_Tool.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <GC_MakeSegment.hxx>
#include <GProp_GProps.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <Standard_Failure.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopLoc_Location.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/profile_solid.hpp"

// Project includes (private, not installed)
#include "detail/library_lock.hpp"
#include "occt_brep.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kPsoPrefix = "meshio++: profile_solid: ";

gp_Pnt pso_pnt(const std::array<double, 3>& rP) {
    return gp_Pnt(rP[0], rP[1], rP[2]);
}

// A loop as an OpenCASCADE wire.
TopoDS_Wire pso_wire(const std::vector<ProfilePiece>& rLoop, const gp_Dir& rNormal) {
    BRepBuilderAPI_MakeWire wire;
    for (const ProfilePiece& r_p : rLoop) {
        TopoDS_Edge edge;
        if (r_p.mType == 1) {
            edge = BRepBuilderAPI_MakeEdge(pso_pnt(r_p.mStart), pso_pnt(r_p.mEnd)).Edge();
        } else if (r_p.mType == 3) {
            const gp_Circ circ(gp_Ax2(pso_pnt(r_p.mCentre), rNormal), r_p.mRadius);
            edge = BRepBuilderAPI_MakeEdge(circ).Edge();
        } else {
            // An arc from start to end, counter-clockwise about the normal for
            // a positive sweep.
            const gp_Dir axis = r_p.mSweepDeg >= 0.0 ? rNormal : rNormal.Reversed();
            const gp_Circ circ(gp_Ax2(pso_pnt(r_p.mCentre), axis), r_p.mRadius);
            edge =
                BRepBuilderAPI_MakeEdge(
                    GC_MakeArcOfCircle(circ, pso_pnt(r_p.mStart), pso_pnt(r_p.mEnd), true).Value())
                    .Edge();
        }
        wire.Add(edge);
        if (!wire.IsDone())
            throw std::invalid_argument(std::string(kPsoPrefix) +
                                        "OpenCASCADE could not chain a loop's edges");
    }
    return wire.Wire();
}

ProfileSolidResult pso_build(const ProfileSolidRequest& rReq) {
    detail::LibraryLock lock;
    detail::OcctQuiet quiet;
    try {
        const gp_Dir normal(rReq.mNormal[0], rReq.mNormal[1], rReq.mNormal[2]);
        // The face: the outer loop with the holes.
        const TopoDS_Wire outer = pso_wire(rReq.mLoops[0], normal);
        BRepBuilderAPI_MakeFace make_face(outer, true);
        for (std::size_t l = 1; l < rReq.mLoops.size(); ++l) {
            TopoDS_Wire hole = pso_wire(rReq.mLoops[l], normal);
            make_face.Add(TopoDS::Wire(hole.Reversed()));
        }
        if (!make_face.IsDone())
            throw std::invalid_argument(std::string(kPsoPrefix) +
                                        "OpenCASCADE could not make the profile's face");
        const TopoDS_Face face = make_face.Face();

        TopoDS_Shape shape;
        TopoDS_Shape end_face;  // the face an extrusion ends at, for a shell
        ProfileSolidResult result;
        for (std::size_t k = 0; k < rReq.mFeatures.size(); ++k) {
            const ProfileFeature& f = rReq.mFeatures[k];
            switch (f.mKind) {
                case ProfileFeatureKind::Extrude: {
                    const double d = f.mDistance;
                    if (f.mSymmetric) {
                        gp_Trsf shift;
                        shift.SetTranslation(gp_Vec(normal) * (-0.5 * d));
                        const TopoDS_Shape moved = face.Moved(TopLoc_Location(shift));
                        BRepPrimAPI_MakePrism prism(moved, gp_Vec(normal) * d);
                        shape = prism.Shape();
                        end_face = prism.LastShape();
                    } else {
                        BRepPrimAPI_MakePrism prism(face, gp_Vec(normal) * d);
                        shape = prism.Shape();
                        end_face = prism.LastShape();
                    }
                    break;
                }
                case ProfileFeatureKind::Revolve: {
                    const gp_Ax1 axis(
                        pso_pnt(f.mAxisOrigin),
                        gp_Dir(f.mAxisDirection[0], f.mAxisDirection[1], f.mAxisDirection[2]));
                    const double angle = f.mAngleDeg * std::acos(-1.0) / 180.0;
                    BRepPrimAPI_MakeRevol revol(face, axis, angle);
                    shape = revol.Shape();
                    end_face = f.mAngleDeg < 360.0 ? revol.LastShape() : TopoDS_Shape();
                    break;
                }
                case ProfileFeatureKind::Fillet: {
                    BRepFilletAPI_MakeFillet fillet(shape);
                    TopTools_IndexedMapOfShape edges;
                    TopExp::MapShapes(shape, TopAbs_EDGE, edges);
                    for (int e = 1; e <= edges.Extent(); ++e)
                        fillet.Add(f.mRadius, TopoDS::Edge(edges(e)));
                    fillet.Build();
                    if (!fillet.IsDone())
                        throw std::invalid_argument(std::string(kPsoPrefix) +
                                                    "OpenCASCADE could not fillet every edge "
                                                    "by that radius");
                    shape = fillet.Shape();
                    end_face = TopoDS_Shape();
                    break;
                }
                case ProfileFeatureKind::Chamfer: {
                    BRepFilletAPI_MakeChamfer chamfer(shape);
                    TopTools_IndexedMapOfShape edges;
                    TopExp::MapShapes(shape, TopAbs_EDGE, edges);
                    for (int e = 1; e <= edges.Extent(); ++e)
                        chamfer.Add(f.mDistance, TopoDS::Edge(edges(e)));
                    chamfer.Build();
                    if (!chamfer.IsDone())
                        throw std::invalid_argument(std::string(kPsoPrefix) +
                                                    "OpenCASCADE could not chamfer every edge "
                                                    "by that distance");
                    shape = chamfer.Shape();
                    end_face = TopoDS_Shape();
                    break;
                }
                case ProfileFeatureKind::Shell: {
                    TopTools_ListOfShape removed;
                    if (!end_face.IsNull())
                        removed.Append(end_face);
                    BRepOffsetAPI_MakeThickSolid thick;
                    thick.MakeThickSolidByJoin(shape, removed, -f.mDistance, 1e-6);
                    if (!thick.IsDone())
                        throw std::invalid_argument(std::string(kPsoPrefix) +
                                                    "OpenCASCADE could not shell the solid to "
                                                    "that wall");
                    shape = thick.Shape();
                    end_face = TopoDS_Shape();
                    break;
                }
            }
        }
        GProp_GProps volume, area;
        BRepGProp::VolumeProperties(shape, volume);
        BRepGProp::SurfaceProperties(shape, area);
        result.mVolume = std::fabs(volume.Mass());
        result.mArea = area.Mass();
        TopTools_IndexedMapOfShape faces;
        TopExp::MapShapes(shape, TopAbs_FACE, faces);
        result.mNumFaces = faces.Extent();
        result.mValid = BRepCheck_Analyzer(shape).IsValid();
        if (!result.mValid)
            result.mWarnings.push_back("OpenCASCADE's check rejects the solid");
        if (!rReq.mStepPath.empty())
            detail::occt_write_step(shape, rReq.mStepPath, "meshio++ profile_solid", kPsoPrefix);
        const detail::OcctNurbs nurbs = detail::occt_nurbs(shape);
        result.mMesh = detail::occt_tessellate(nurbs.mShape, nurbs.mInfo, rReq.mDeflection, 0.3);
        return result;
    } catch (const Standard_Failure& rErr) {
        throw std::invalid_argument(std::string(kPsoPrefix) +
                                    "OpenCASCADE failed: " + rErr.GetMessageString());
    }
}

}  // namespace

namespace providers {

ProfileSolidProvider make_occt_profile_provider() {
    ProfileSolidProvider provider;
    provider.mName = "occt";
    provider.mBuild = pso_build;
    return provider;
}

}  // namespace providers

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_OCCT
