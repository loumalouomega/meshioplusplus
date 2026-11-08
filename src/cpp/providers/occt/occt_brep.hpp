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
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// External includes (OpenCASCADE)
#include <Message_Messenger.hxx>
#include <Poly_Triangulation.hxx>
#include <TopoDS_Shape.hxx>
#include <gp_Trsf.hxx>

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
 * bounding-box diagonal, raised to the shape's own largest tolerance but never
 * past 1e-3 of the diagonal.
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

/// A tessellation node's welding key: the bits of its coordinates, `-0` read
/// as `+0`, so nodes weld only where they coincide exactly.
struct OcctNodeKey {
    std::array<std::uint64_t, 3> mBits;
    bool operator==(const OcctNodeKey& rO) const { return mBits == rO.mBits; }
};

/// FNV-1a over the key's three words. Buckets only: the map is looked up,
/// never iterated, so the hash cannot reach the output.
struct OcctNodeKeyHash {
    std::size_t operator()(const OcctNodeKey& rK) const {
        std::uint64_t h = 1469598103934665603ull;
        for (std::uint64_t b : rK.mBits)
            h = (h ^ b) * 1099511628211ull;
        return static_cast<std::size_t>(h);
    }
};

/// The welded id of every node met so far, across the faces of one tessellation.
using OcctNodeIds = std::unordered_map<OcctNodeKey, std::int64_t, OcctNodeKeyHash>;

/**
 * @brief Welds the nodes of one face's triangulation into a shared point list:
 * each node, moved by @p rTrsf, takes the id of an earlier node with exactly
 * its coordinates, or the next new id, appending its `x y z` to @p rXyz (ids in
 * first-seen order).
 * @return The welded id of each node, indexed 1-based as the triangulation's
 *         own (entry 0 unused).
 */
std::vector<std::int64_t> occt_weld_nodes(const Poly_Triangulation& rTri, const gp_Trsf& rTrsf,
                                          OcctNodeIds& rIds, std::vector<double>& rXyz);

/**
 * @brief Visits the triangles of a welded triangulation, in its order, turned
 * over when @p Reversed (the face is `TopAbs_REVERSED`), skipping those the
 * weld collapsed: `rVisit(ia, ib, ic, a, b, c)`, the welded ids and the
 * triangulation's own 1-based node numbers.
 */
template <class TVisit>
void occt_welded_triangles(const Poly_Triangulation& rTri, bool Reversed,
                           const std::vector<std::int64_t>& rLocal, TVisit&& rVisit) {
    for (int t = 1; t <= rTri.NbTriangles(); ++t) {
        int a = 0, b = 0, c = 0;
        rTri.Triangle(t).Get(a, b, c);
        if (Reversed)
            std::swap(b, c);
        const std::int64_t ia = rLocal[static_cast<std::size_t>(a)];
        const std::int64_t ib = rLocal[static_cast<std::size_t>(b)];
        const std::int64_t ic = rLocal[static_cast<std::size_t>(c)];
        if (ia == ib || ib == ic || ia == ic)
            continue;
        rVisit(ia, ib, ic, a, b, c);
    }
}

}  // namespace meshioplusplus::detail

#endif  // MESHIOPLUSPLUS_HAS_OCCT
