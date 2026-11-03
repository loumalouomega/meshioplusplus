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
#pragma once

/**
 * @file detail/occt_brep.hpp
 * @brief An OpenCASCADE shape as meshio++ sees it: its faces, loops and edges
 * as a `BrepInfo` (NURBS throughout) and its triangulation as a `Mesh`. The
 * STEP and IGES readers and `fit_cad` share it.
 *
 * The shape is first converted with `BRepBuilderAPI_NurbsConvert`, which turns
 * every surface into a `Geom_BSplineSurface` and every parameter curve into a
 * `Geom2d_BSplineCurve` *consistently*: converting a face's surface alone
 * reparametrises a cylinder or a sphere and leaves its trimming curves in the
 * old parameters. OpenCASCADE 7.9 throws converting some trimmed and multi-turn
 * elementary faces, and leaves a face whose surface already is a B-spline with
 * none; then the shape is kept as it is and every surface and curve that is
 * not a B-spline is approximated by one *in its own parameters*, so the
 * trimming curves stay valid. Faces get ids `1 .. F` and edges `F + 1 .. F + E` (Kratos
 * gives every geometry one id space), control points ids from 1, so the
 * `BrepInfo` writes as CAD JSON Kratos loads.
 *
 * Core-private (not installed); compiled only with `MESHIOPLUSPLUS_HAS_OCCT`.
 */

#ifdef MESHIOPLUSPLUS_HAS_OCCT

// System includes
#include <string>
#include <type_traits>
#include <utility>

// External includes (OpenCASCADE)
#include <Message_Messenger.hxx>
#include <TopoDS_Shape.hxx>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus::detail {

/// The shape a STEP (or, with @p Iges, an IGES) file transfers to.
/// @throws ReadError when OpenCASCADE cannot read or transfer it.
TopoDS_Shape occt_read_cad_file(const std::string& rPath, bool Iges, const char* pFormat);

/// What `occt_heal_shape` did to a shape read from a STEP or IGES file.
struct OcctHealReport {
    /// `ShapeFix_Shape` ran (`CadReadOptions::mHeal`).
    bool mHealed = false;
    /// `ShapeFix_Shape` changed the shape.
    bool mFixed = false;
    /// The sewing tolerance used, as a length (0 when no sewing ran).
    double mSewTolerance = 0.0;
    /// Edges with exactly one face, before and after (degenerate edges excluded).
    long long mFreeEdgesBefore = 0;
    long long mFreeEdgesAfter = 0;
    /// The edges sewing reports as shared by more than two faces.
    long long mMultipleEdges = 0;
    /// The degenerate shapes sewing removed.
    long long mDegenerated = 0;
    /// The report as a compact JSON object.
    std::string ToJson() const;
};

/// Edges of the shape that bound exactly one face (degenerate edges excluded).
long long occt_count_free_edges(const TopoDS_Shape& rShape);

/**
 * @brief Heal a shape at read time (`CadReadOptions::mHeal`, `mSewTolerance`):
 * `ShapeFix_Shape` when healing, then `BRepBuilderAPI_Sewing` when healing or
 * when a sewing tolerance is given. A tolerance at most 0 takes 1e-5 of the
 * bounding-box diagonal, never below the shape's own largest tolerance.
 */
TopoDS_Shape occt_heal_shape(const TopoDS_Shape& rShape, const CadReadOptions& rOptions,
                             OcctHealReport& rReport);

/// A shape and its faces, loops and edges as NURBS.
struct OcctNurbs {
    /// The converted shape, or the shape itself when the conversion failed.
    TopoDS_Shape mShape;
    BrepInfo mInfo;
    /// The conversion succeeded on every face: the NURBS are exact.
    bool mExact = false;
};

/**
 * @brief The shape as NURBS: converted exactly by `BRepBuilderAPI_NurbsConvert`
 * when that leaves every face a B-spline surface, otherwise approximated face by
 * face and edge by edge in the shape's own parameters, within 1e-6 of its
 * bounding-box diagonal.
 */
OcctNurbs occt_nurbs(const TopoDS_Shape& rShape);

/// The faces, loops and edges of a NURBS-converted shape.
BrepInfo occt_brep_info(const TopoDS_Shape& rNurbs);

/**
 * @brief The shape triangulated by `BRepMesh_IncrementalMesh`: one `triangle`
 * block, points shared between faces welded where they coincide exactly, and
 * the `cad:*` cell data of `BrepInfo` (`cad:span` from each triangle's centre
 * in the face's parameters, -1 where the triangulation has none).
 * @param Deflection the chordal deflection (a length); at most 0 takes 1e-3 of
 *        the bounding-box diagonal.
 * @param Angular the angular deflection in radians.
 */
Mesh occt_tessellate(const TopoDS_Shape& rNurbs, const BrepInfo& rInfo, double Deflection,
                     double Angular);

/// The printers a messenger holds, named through its accessor: OpenCASCADE 8.0
/// deprecates the typedef `Message_SequenceOfPrinters`, and no header includes it.
using OcctPrinters = std::remove_cvref_t<decltype(std::declval<Message_Messenger&>().Printers())>;

/// Detach OpenCASCADE's message printers (it prints to stdout) while alive.
class OcctQuiet {
public:
    OcctQuiet();
    OcctQuiet(const OcctQuiet&) = delete;
    OcctQuiet& operator=(const OcctQuiet&) = delete;
    ~OcctQuiet();

private:
    struct Impl;
    Impl* mpImpl;
};

}  // namespace meshioplusplus::detail

#endif  // MESHIOPLUSPLUS_HAS_OCCT
