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
// DXF R2000 drawings (roadmap 7.6.4).
// STEP curves and element solids (roadmap 7.6.5).

#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/operations/fit_curves.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "meshioplusplus/operations/silhouette.hpp"
#include "meshioplusplus/skin.hpp"

using namespace meshioplusplus;

namespace {

std::string se_tmp(const std::string& rName) {
    return (std::filesystem::temp_directory_path() / ("mio_step_el_" + rName)).string();
}

std::string se_slurp(const std::string& rPath) {
    std::ifstream f(rPath, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::size_t se_count(const std::string& rText, const std::string& rWhat) {
    std::size_t n = 0;
    for (std::size_t p = rText.find(rWhat); p != std::string::npos; p = rText.find(rWhat, p + 1))
        ++n;
    return n;
}

Mesh se_fitted(const Mesh& rSurface) {
    SegmentOptions so;
    so.mMethod = SegmentMethod::Primitive;
    return fit_primitives(segment(rSurface, so).mMesh).mMesh;
}

}  // namespace

TEST(StepElements, AFittedCylinderIsALiddedSolid) {
    const Mesh fitted = se_fitted(extract_skin(cylinder(2.0, 5.0, 12, 10)));
    const std::string path = se_tmp("cyl.step");
    const StepElementsResult r = write_step_elements(path, fitted);
    EXPECT_EQ(r.mNumCylinders, 1);
    EXPECT_EQ(r.mNumPlanes, 2);  // the two end faces, bounded
    const std::string text = se_slurp(path);
    EXPECT_EQ(se_count(text, "CYLINDRICAL_SURFACE("), 1u);
    EXPECT_NE(text.find("CYLINDRICAL_SURFACE('',#"), std::string::npos);
    EXPECT_GE(se_count(text, "CIRCLE("), 2u);  // the rims
    EXPECT_EQ(se_count(text, "MANIFOLD_SOLID_BREP("), 1u);
    EXPECT_EQ(se_count(text, "SHELL_BASED_SURFACE_MODEL("), 2u);
    EXPECT_EQ(se_count(text, "SHAPE_REPRESENTATION_RELATIONSHIP("), 2u);
    std::filesystem::remove(path);
}

TEST(StepElements, ASphereIsAnExactBall) {
    const Mesh fitted = se_fitted(icosa_sphere(1.5, 3));
    const std::string path = se_tmp("ball.step");
    const StepElementsResult r = write_step_elements(path, fitted);
    EXPECT_EQ(r.mNumSpheres, 1);
    const std::string text = se_slurp(path);
    EXPECT_NE(text.find("SPHERICAL_SURFACE('',#"), std::string::npos);
    EXPECT_NE(text.find(",1.5)"), std::string::npos);
    EXPECT_EQ(se_count(text, "VERTEX_LOOP("), 1u);
    EXPECT_NE(text.find("ADVANCED_BREP_SHAPE_REPRESENTATION("), std::string::npos);
    std::filesystem::remove(path);
}

TEST(StepElements, SelectionAndRefusals) {
    const Mesh fitted = se_fitted(extract_skin(box({{2, 2, 2}})));
    const std::string path = se_tmp("sel.step");
    StepElementsOptions o;
    o.mPatches = {0};
    EXPECT_EQ(write_step_elements(path, fitted, o).mNumPlanes, 1);
    o.mPatches = {99};
    EXPECT_THROW(write_step_elements(path, fitted, o), WriteError);
    EXPECT_THROW(write_step_elements(path, icosa_sphere()), WriteError);  // no fit table
    StepElementsOptions few;
    few.mSegments = 3;
    EXPECT_THROW(write_step_elements(path, fitted, few), WriteError);
    std::filesystem::remove(path);
}

TEST(StepCurves, ASketchIsTrimmedLinesAndArcs) {
    // The silhouette of a cylinder along its axis is one circle.
    SilhouetteOptions so;
    const Mesh outline = silhouette(extract_skin(cylinder(1.0, 2.0, 12, 2)), so).mMesh;
    FitCurvesOptions fo;
    fo.mTolerance = 1e-3;
    const Mesh sketch = fit_curves(outline, fo).mMesh;
    const std::string path = se_tmp("sketch.step");
    write_step(path, sketch);
    const std::string text = se_slurp(path);
    EXPECT_NE(text.find("GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION("),
              std::string::npos);
    EXPECT_NE(text.find("GEOMETRIC_CURVE_SET('sketch'"), std::string::npos);
    EXPECT_GE(se_count(text, "CIRCLE("), 1u);
    // A plain line mesh is polylines.
    write_step(path, outline);
    const std::string poly = se_slurp(path);
    EXPECT_NE(poly.find("GEOMETRIC_CURVE_SET('sections'"), std::string::npos);
    EXPECT_EQ(se_count(poly, "POLYLINE("), 1u);
    std::filesystem::remove(path);
}
