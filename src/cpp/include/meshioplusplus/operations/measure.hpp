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
 * @file operations/measure.hpp
 * @brief Elements, dimensions and GD&T tolerances measured from a scan
 * (roadmap 7.2), and the 3-2-1 datum alignment on them (roadmap 7.3.4).
 *
 * `measure(mesh, spec)` reads a JSON spec (doc/measure.md):
 *
 *  - **elements**, each with a unique name: measured from the mesh by
 *    `fit_primitives` (`"fit": "plane" | "sphere" | "cylinder" | "cone" |
 *    "auto"` on a named cell `"region"`, from `"seeds"`, or the row
 *    `"patch"` of a fitted mesh's table), constructed from earlier elements
 *    (`"construct": "intersection" | "axis" | "centre" | "apex" |
 *    "mid_plane" | "projection" | "plane_through"`, `"of": [names]`), or a
 *    literal `"point"`. A constructed element records its construction, so a
 *    run re-evaluates the whole spec in order and a refit propagates;
 *  - **dimensions**: `"distance"` and `"angle"` between two elements, with an
 *    optional `"nominal"` and `"tolerance"` (PASS when within);
 *  - **tolerances**: form (`"flatness"`, `"cylindricity"`, `"sphericity"`,
 *    `"circularity"`), orientation (`"parallelism"`, `"perpendicularity"`,
 *    `"angularity"` at a `"basic_angle"`) and location (`"coaxiality"`,
 *    `"concentricity"`) of a `"feature"` against a `"datum"`, with an optional
 *    `"limit"` (PASS when within);
 *  - `"criterion"`: the reference the form tolerances are measured against
 *    (`fit_primitives`' criterion; least squares by default).
 *
 * The result is the element table, the report (both also as JSON), and the
 * input mesh with the element table as the field data `measure:table` (the
 * JSON text's bytes, UInt8) and each fitted element's cells as the cell
 * region `measure:<name>`.
 */

// System includes
#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// The field data holding a mesh's element table (JSON text, UInt8).
inline constexpr const char* kMeasureTableName = "measure:table";

/// What an element is.
enum class ElementKind : std::int32_t {
    Point = 1,
    Line = 2,
    Plane = 3,
    Sphere = 4,
    Cylinder = 5,
    Cone = 6,
};

/// One element of the table.
struct MeasureElement {
    std::string mName;
    ElementKind mKind = ElementKind::Point;
    /// A point: the point; a line: a point on it; a plane: a point on it; a
    /// sphere: its centre; a cylinder or cone: a point on its axis (a cone's
    /// apex).
    std::array<double, 3> mPoint{};
    /// A line's direction, a plane's normal, an axis' direction (unit).
    std::array<double, 3> mDir{};
    double mRadius = 0.0;    ///< a sphere's or cylinder's
    double mAngleDeg = 0.0;  ///< a cone's half-angle
    /// Measured: the form error under the spec's criterion, else NaN.
    double mFormError = 0.0;
    bool mMeasured = false;
    /// The construction (`"intersection"`, ...) or the fit (`"fit:plane"`, ...).
    std::string mSource;
    std::vector<std::string> mInputs;  ///< a construction's elements
    /// Measured: the points the fit used (x, y, z each).
    std::vector<double> mPoints;
};

/// One dimension or tolerance, as reported.
struct MeasureCheck {
    std::string mName;
    std::string mType;             ///< "distance", "angle", "flatness", ...
    std::vector<std::string> mOf;  ///< the elements: feature then datum
    double mValue = 0.0;
    bool mHasNominal = false;
    double mNominal = 0.0;
    bool mHasTolerance = false;  ///< dimensions: +- tolerance; tolerances: the limit
    double mTolerance = 0.0;
    double mDeviation = 0.0;  ///< dimensions: value - nominal
    bool mPass = true;        ///< true without a nominal or limit
};

struct MeasureResult {
    Mesh mMesh;
    std::vector<MeasureElement> mElements;
    std::vector<MeasureCheck> mDimensions;
    std::vector<MeasureCheck> mTolerances;
    std::int64_t mNumFail = 0;
    std::string mTableJson;   ///< the element table
    std::string mReportJson;  ///< elements, dimensions, tolerances, num_fail
    std::vector<std::string> mWarnings;
};

/**
 * @brief Measure @p rMesh by the JSON spec @p rSpecJson (see the file comment).
 * @throws std::invalid_argument on a malformed spec, an unknown or duplicate
 *         name, an element a construction cannot take, or a fit that fails.
 */
MESHIOPLUSPLUS_API MeasureResult measure(const Mesh& rMesh, const std::string& rSpecJson);

/// The element table a mesh carries (`measure:table`); empty when none.
MESHIOPLUSPLUS_API std::vector<MeasureElement> measure_table(const Mesh& rMesh);

/// The element named @p rName of @p rMesh's table.
/// @throws std::invalid_argument when the mesh has no table or no such element.
MESHIOPLUSPLUS_API MeasureElement measure_element(const Mesh& rMesh, const std::string& rName);

/**
 * @brief A plane, sphere, cylinder or cone element as the fit table's type and
 * ten parameters (`fit:patch_params`), the reference `deviation_to_primitive`
 * takes (roadmap 7.4.2).
 * @throws std::invalid_argument for a point or a line, which have no surface.
 */
MESHIOPLUSPLUS_API std::pair<std::int64_t, std::array<double, 10>> measure_element_primitive(
    const MeasureElement& rElement);

/// What `align_datums` did.
struct AlignDatumsResult {
    Mesh mMesh;                        ///< moved, its element table with it
    std::array<double, 16> mMatrix{};  ///< row-major, world to the datum frame
    std::vector<std::string> mWarnings;
};

/**
 * @brief The 3-2-1 datum alignment (roadmap 7.3.4): the primary element (a
 * plane's normal, or a cylinder's or line's axis) sets Z, the secondary
 * (projected) sets X, and the origin is a point-like element, or the meeting
 * point of the primary, the secondary and a third plane named by @p rOrigin;
 * with no origin, the primary's point projected along the frame. The mesh is
 * moved by the matrix, its element table with it.
 * @param rMesh a mesh carrying an element table (`measure`).
 * @throws std::invalid_argument on an unknown element, parallel primary and
 *         secondary, or planes that do not meet in a point.
 */
MESHIOPLUSPLUS_API AlignDatumsResult align_datums(const Mesh& rMesh, const std::string& rPrimary,
                                                  const std::string& rSecondary,
                                                  const std::string& rOrigin = "");

}  // namespace meshioplusplus
