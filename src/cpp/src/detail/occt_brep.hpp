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

// External includes (OpenCASCADE)
#include <TopoDS_Shape.hxx>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus::detail {

/// The shape a STEP (or, with @p Iges, an IGES) file transfers to.
/// @throws ReadError when OpenCASCADE cannot read or transfer it.
TopoDS_Shape occt_read_cad_file(const std::string& rPath, bool Iges, const char* pFormat);

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
