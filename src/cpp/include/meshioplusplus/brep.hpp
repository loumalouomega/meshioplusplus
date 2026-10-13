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
 * @file brep.hpp
 * @brief `BrepInfo`: a NURBS boundary representation -- trimmed NURBS faces,
 * their trimming loops, edges and vertices -- kept beside the `Mesh` a reader
 * sampled it into (the `NurbsInfo`/`MdpaInfo` side-channel pattern).
 *
 * A `Mesh` holds linear and Lagrange polytopes only, so a reader of a B-rep
 * file samples it: Kratos CAD JSON (`formats/cad_json.hpp`, the IGA geometry
 * of Kratos Multiphysics) evaluates every knot span of every face, and STEP
 * and IGES (`formats/step.hpp`, `formats/iges.hpp`, through the optional
 * OpenCASCADE backend) are triangulated at a tolerance. The overloads that
 * take a `BrepInfo` also return the exact geometry the sampling drops; the
 * CAD JSON writer writes it back, so a CAD JSON file read and written keeps
 * every face, loop, edge and vertex. The shared registry passes none.
 *
 * Every patch is a `NurbsPatch` of parametric dimension 2 with its knot
 * vectors in the full form (`n + p + 1` knots); a file that writes them in
 * Kratos's reduced form (`n + p - 1`, without the outer knots) sets
 * `mReducedKnots` so the writer can restore it. Weights are carried as the
 * file writes them; a face or curve whose `mRational` is false is evaluated
 * as a polynomial B-spline whatever they are (as Kratos does).
 *
 * The members a format has that meshio++ does not model are kept, per object,
 * in `mExtra` as compact JSON object text (`{}` or empty for none), and
 * written back. Every index is 0-based. See doc/formats/cad_json.md.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/nurbs.hpp"

namespace meshioplusplus {

/// A NURBS curve: a 3-D edge curve or a trimming curve in a face's (u, v).
struct MESHIOPLUSPLUS_API BrepCurve {
    /// Polynomial degree.
    std::int32_t mDegree = 1;
    /// Full knot vector (`n + degree + 1` knots).
    std::vector<double> mKnots;
    /// Whether the file wrote the reduced form (without the two outer knots).
    bool mReducedKnots = false;
    /// Whether the curve is rational (weights used).
    bool mRational = false;
    /// Float64 `(n, 3)`: control points; a trimming curve's are `(u, v, 0)`.
    NDArray mControlPoints{DType::Float64, {0, 3}};
    /// Float64 `(n)`: the weights as written.
    NDArray mWeights{DType::Float64, {0}};
    /// Int64 `(n)`: the control points' ids (Kratos node ids); empty when the
    /// file gives none.
    NDArray mIds{DType::Int64, {0}};
    /// The active parameter range `[t0, t1]`, or empty for none.
    std::vector<double> mActiveRange;
    /// Members not modelled, as compact JSON object text.
    std::string mExtra;
};

/// One trimming curve of a loop (or an embedded edge of a face).
struct MESHIOPLUSPLUS_API BrepTrim {
    /// The trim index the edges refer to; -1 when the file gives none.
    std::int64_t mTrimIndex = -1;
    /// Whether the loop runs along the curve's direction (`curve_direction`).
    bool mDirection = true;
    /// The curve in the face's parameter space.
    BrepCurve mCurve;
    std::string mExtra;
};

/// A boundary loop of a face.
struct MESHIOPLUSPLUS_API BrepLoop {
    /// `"outer"` or `"inner"`.
    std::string mType = "outer";
    std::vector<BrepTrim> mTrims;
    std::string mExtra;
};

/// A reference from an edge or vertex to the face (and trim) it lies on.
struct MESHIOPLUSPLUS_API BrepRef {
    /// The referenced entity's id; -1 when it is named instead.
    std::int64_t mId = -1;
    /// The referenced entity's name (`brep_name`), or empty.
    std::string mName;
    /// The trim index on that face; -1 when absent.
    std::int64_t mTrimIndex = -1;
    /// `relative_direction`: 1 true, 0 false, -1 absent.
    std::int32_t mRelativeDirection = -1;
    /// A vertex's parametric position on the face, or empty.
    std::vector<double> mLocalCoordinates;
    std::string mExtra;
};

/// A trimmed NURBS face.
struct MESHIOPLUSPLUS_API BrepFace {
    /// `brep_id`; -1 when the face is named instead.
    std::int64_t mId = -1;
    /// `brep_name`, or empty.
    std::string mName;
    /// Its surface: an index into `BrepInfo::mPatches`.
    std::int64_t mPatch = -1;
    bool mSwappedNormal = false;
    /// Whether the surface is trimmed by `mLoops` (`is_trimmed`).
    bool mTrimmed = false;
    bool mRational = false;
    bool mReducedKnots = false;
    /// Int64 `(n)`: the surface's control-point ids, or empty.
    NDArray mControlPointIds{DType::Int64, {0}};
    std::vector<BrepLoop> mLoops;
    /// Curves embedded in the face (`embedded_edges`).
    std::vector<BrepTrim> mEmbeddedEdges;
    /// The surface object's members not modelled.
    std::string mSurfaceExtra;
    std::string mExtra;
};

/// An edge: an optional 3-D curve and the trims it is (its topology).
struct MESHIOPLUSPLUS_API BrepEdge {
    std::int64_t mId = -1;
    std::string mName;
    bool mHasCurve = false;
    /// The `3d_curve`, when `mHasCurve`.
    BrepCurve mCurve;
    std::vector<BrepRef> mTopology;
    std::string mExtra;
};

/// A vertex: the points on faces or edges it is.
struct MESHIOPLUSPLUS_API BrepVertex {
    std::int64_t mId = -1;
    std::string mName;
    std::vector<BrepRef> mTopology;
    std::string mExtra;
};

/// One `breps` entry: indices into `BrepInfo`'s face, edge and vertex lists.
struct MESHIOPLUSPLUS_API Brep {
    std::int64_t mId = -1;
    std::string mName;
    std::vector<std::int64_t> mFaces;
    std::vector<std::int64_t> mEdges;
    std::vector<std::int64_t> mVertices;
    std::string mExtra;
};

/**
 * @brief The B-rep a sampled mesh came from. Empty (`mFaces`, `mEdges` and
 * `mVertices` empty) for a mesh that has none.
 *
 * The sampled mesh links its cells back through cell data: `cad:kind` (0 a
 * face cell, 1 an edge cell), `cad:index` (into `mFaces` or `mEdges`),
 * `cad:id` (that entity's id), `cad:span` (the knot span per direction,
 * counted over the spans of nonzero length, -1 for none) and `cad:trimmed`
 * (1 for a cell the trimming crossed).
 */
struct MESHIOPLUSPLUS_API BrepInfo {
    /// The face surfaces.
    std::vector<NurbsPatch> mPatches;
    std::vector<Brep> mBreps;
    std::vector<BrepFace> mFaces;
    std::vector<BrepEdge> mEdges;
    std::vector<BrepVertex> mVertices;
    /// The file's top-level members besides `breps` (Kratos: `tolerances`,
    /// `version_number`), as compact JSON object text.
    std::string mExtra;
};

/// How a STEP or IGES B-rep is triangulated (`read_step`, `read_iges`).
struct CadReadOptions {
    /// The chordal deflection, as a length; at most 0 takes 1e-3 of the
    /// bounding-box diagonal.
    double mDeflection = 0.0;
    /// The angular deflection, in radians.
    double mAngularDeflection = 0.5;
};

/// Whether this build reads STEP and IGES (it carries the OpenCASCADE backend,
/// `-DMESHIOPLUSPLUS_WITH_OCCT=ON`).
MESHIOPLUSPLUS_API bool cad_read_has_occt() noexcept;

}  // namespace meshioplusplus
