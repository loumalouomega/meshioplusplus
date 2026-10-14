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
 * old parameters. Faces get ids `1 .. F` and edges `F + 1 .. F + E` (Kratos
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

/// The shape with every surface and curve as a NURBS.
TopoDS_Shape occt_nurbs_shape(const TopoDS_Shape& rShape);

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
