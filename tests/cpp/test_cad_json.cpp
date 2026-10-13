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

/**
 * @file test_cad_json.cpp
 * @brief Kratos CAD JSON: the sampled mesh of untrimmed, trimmed and rational
 *        faces, the BrepInfo side channel, a byte-stable rewrite, the plain
 *        mesh writer, registry dispatch and the error paths, on the native
 *        path directly (the Python tests cover parity with the reference).
 */

// External includes
#include <gtest/gtest.h>

// System includes
#include <cmath>
#include <filesystem>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/formats/mdpa.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::BrepInfo;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

std::string write_file(const std::string& rBody, const std::string& rSuffix = ".cad.json") {
    const std::string path = mt::temp_path(rSuffix);
    std::ofstream(path, std::ios::binary) << rBody;
    return path;
}

std::string slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// A rational quadratic quarter circle of radius 1 as a 3-D edge, a degree-1
// unit square face trimmed by a square outer loop and a diamond hole, and
// Kratos's reduced knots on the face.
const char* kModel = R"({
    // comments are accepted, as Kratos's Parameters accept them
    "version_number": 1,
    "breps": [{
        "brep_id": 1,
        "faces": [{
            "brep_id": 2,
            "swapped_surface_normal": false,
            "surface": {
                "is_trimmed": true,
                "is_rational": false,
                "degrees": [1, 1],
                "knot_vectors": [[0, 4], [0, 4]],
                "control_points": [[1, [0, 0, 0, 1]], [2, [4, 0, 0, 1]],
                                   [3, [0, 4, 0, 1]], [4, [4, 4, 0, 1]]]
            },
            "boundary_loops": [
                {"loop_type": "outer", "trimming_curves": [
                    {"trim_index": 0, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[0, 0, 0, 1]], [[4, 0, 0, 1]]]}},
                    {"trim_index": 1, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[4, 0, 0, 1]], [[4, 4, 0, 1]]]}},
                    {"trim_index": 2, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[4, 4, 0, 1]], [[0, 4, 0, 1]]]}},
                    {"trim_index": 3, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 1, 1], "active_range": [0, 1],
                        "control_points": [[[0, 4, 0, 1]], [[0, 0, 0, 1]]]}}]},
                {"loop_type": "inner", "trimming_curves": [
                    {"trim_index": 4, "curve_direction": true, "parameter_curve": {
                        "degree": 1, "knot_vector": [0, 0, 0.25, 0.5, 0.75, 1, 1], "active_range": [0, 1],
                        "control_points": [[[2, 0.8, 0, 1]], [[3.2, 2, 0, 1]], [[2, 3.2, 0, 1]],
                                           [[0.8, 2, 0, 1]], [[2, 0.8, 0, 1]]]}}]}],
            "embedded_edges": [],
            "custom_member": {"kept": [1, 2.5]}
        }],
        "edges": [{
            "brep_id": 3,
            "3d_curve": {
                "is_rational": true,
                "degree": 2,
                "knot_vector": [0, 0, 0, 1, 1, 1],
                "active_range": [0, 1],
                "control_points": [[10, [1, 0, 0, 1]], [11, [1, 1, 0, 0.7071067811865476]],
                                   [12, [0, 1, 0, 1]]]
            },
            "topology": [{"brep_id": 2, "trim_index": 0, "relative_direction": true}]
        }],
        "vertices": []
    }]
})";

std::size_t cells_of(const Mesh& rMesh, const std::string& rType) {
    std::size_t n = 0;
    for (const auto cb : rMesh.CellRange())
        if (cb.Type() == rType)
            n += cb.NumCells();
    return n;
}

TEST(CadJson, ReadsTheSideChannel) {
    BrepInfo info;
    meshioplusplus::read_cad_json(write_file(kModel), info);
    ASSERT_EQ(info.mFaces.size(), 1u);
    ASSERT_EQ(info.mEdges.size(), 1u);
    ASSERT_EQ(info.mBreps.size(), 1u);
    const auto& f = info.mFaces[0];
    EXPECT_EQ(f.mId, 2);
    EXPECT_TRUE(f.mTrimmed);
    EXPECT_TRUE(f.mReducedKnots);
    EXPECT_EQ(info.mPatches[0].mKnots[0], (std::vector<double>{0, 0, 4, 4}));
    EXPECT_EQ(f.mControlPointIds.Size(), 4u);
    ASSERT_EQ(f.mLoops.size(), 2u);
    EXPECT_EQ(f.mLoops[1].mType, "inner");
    EXPECT_EQ(f.mLoops[0].mTrims[3].mTrimIndex, 3);
    EXPECT_EQ(f.mLoops[0].mTrims[0].mCurve.mIds.Size(), 0u);
    EXPECT_EQ(f.mExtra, R"({"custom_member":{"kept":[1,2.5]}})");
    EXPECT_EQ(info.mExtra, R"({"version_number":1})");
    const auto& e = info.mEdges[0];
    EXPECT_TRUE(e.mHasCurve);
    EXPECT_TRUE(e.mCurve.mRational);
    ASSERT_EQ(e.mTopology.size(), 1u);
    EXPECT_EQ(e.mTopology[0].mRelativeDirection, 1);
}

TEST(CadJson, SamplesExactGeometryAndClipsTheHole) {
    BrepInfo info;
    meshioplusplus::CadJsonReadOptions opts;
    opts.mSubdivisions = 8;
    const Mesh mesh = meshioplusplus::read_cad_json(write_file(kModel), info, opts);
    const std::size_t kept = cells_of(mesh, "quad");
    // 64 cells, the diamond |x-2| + |y-2| < 1.2 holding the centres of 12.
    EXPECT_EQ(kept, 52u);
    opts.mTrim = false;
    const Mesh full = meshioplusplus::read_cad_json(write_file(kModel), info, opts);
    EXPECT_EQ(cells_of(full, "quad"), 64u);
    // The rational edge is an exact quarter circle.
    const NDArray& pts = mesh.Points();
    bool any = false;
    for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b) {
        const auto cb = mesh.Cells(b);
        if (cb.Type() != "line3")
            continue;
        const auto* c = cb.Conn().As<std::int64_t>();
        for (std::size_t i = 0; i < cb.NumCells() * 3; ++i) {
            const double* p = pts.As<double>() + c[i] * 3;
            EXPECT_NEAR(std::hypot(p[0], p[1]), 1.0, 1e-15);
            any = true;
        }
    }
    EXPECT_TRUE(any);
    ASSERT_TRUE(mesh.HasCellData("cad:trimmed"));
    ASSERT_TRUE(mesh.HasCellData("cad:span"));
}

TEST(CadJson, RewriteIsByteStable) {
    BrepInfo info;
    const Mesh mesh = meshioplusplus::read_cad_json(write_file(kModel), info);
    const std::string a = mt::temp_path(".cad.json");
    meshioplusplus::write_cad_json(a, mesh, info);
    BrepInfo back;
    const Mesh again = meshioplusplus::read_cad_json(a, back);
    const std::string b = mt::temp_path(".cad.json");
    meshioplusplus::write_cad_json(b, again, back);
    EXPECT_EQ(slurp(a), slurp(b));
    // Reduced knots are written reduced; the extra member is kept.
    EXPECT_NE(slurp(a).find("\"custom_member\": {"), std::string::npos);
    EXPECT_EQ(back.mFaces[0].mExtra, info.mFaces[0].mExtra);
    EXPECT_TRUE(back.mFaces[0].mReducedKnots);
}

TEST(CadJson, PlainMeshWritesBilinearFaces) {
    Mesh mesh = mt::quad_mesh();
    const std::string path = mt::temp_path(".cad.json");
    meshioplusplus::write_cad_json(path, mesh);
    const Mesh back = meshioplusplus::read_cad_json(path);
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).Type(), "quad");
    EXPECT_EQ(back.Cells(0).NumCells(), mesh.Cells(0).NumCells());
    EXPECT_THROW(meshioplusplus::write_cad_json(mt::temp_path(".cad.json"), mt::tet_mesh()),
                 meshioplusplus::WriteError);
}

TEST(CadJson, RegistryDispatch) {
    const std::string path = write_file(kModel);
    EXPECT_EQ(meshioplusplus::resolve_format(path, ""), "cad_json");
    const Mesh mesh = meshioplusplus::registry_readers().at("cad_json")(path);
    EXPECT_GT(mesh.NumCellBlocks(), 0u);
    const std::string out = mt::temp_path(".cad.json");
    meshioplusplus::registry_writers().at("cad_json")(out, mt::quad_mesh());
    EXPECT_EQ(meshioplusplus::resolve_format(out, ""), "cad_json");
}

TEST(CadJson, MalformedInputIsReadError) {
    for (const char* body :
         {"", "[1]", "{\"breps\": [", "{\"x\": 1}",
          R"({"breps":[{"faces":[{"surface":{"degrees":[1,1],"knot_vectors":[[0,1],[0,1]],"control_points":[[[0,0,0,1]]]}}]}]})",
          R"({"breps":[{"edges":[{"3d_curve":{"degree":1,"knot_vector":[1,0],"control_points":[[1,[0,0,0,1]],[2,[1,0,0,1]]]}}]}]})",
          R"({"breps":[{"edges":[{"3d_curve":{"is_rational":true,"degree":1,"knot_vector":[0,1],"control_points":[[1,[0,0,0,-1]],[2,[1,0,0,1]]]}}]}]})",
          // Not UTF-8: a stray byte, an overlong form, a surrogate.
          "{\"breps\":[{\"brep_name\":\"\xfd\"}]}", "{\"breps\":[{\"brep_name\":\"\xc0\xaf\"}]}",
          "{\"breps\":[{\"brep_name\":\"\xed\xa0\x80\"}]}"}) {
        BrepInfo info;
        EXPECT_THROW(meshioplusplus::read_cad_json(write_file(body), info),
                     meshioplusplus::ReadError)
            << body;
    }
}

TEST(CadJson, NonAsciiNamesAreKept) {
    BrepInfo info;
    (void)meshioplusplus::read_cad_json(
        write_file(R"({"breps":[{"brep_name":"Öffnung Ω","faces":[]}]})"
                   "\n// \xc3\xa9 in a comment\n"),
        info);
    ASSERT_EQ(info.mBreps.size(), 1u);
    EXPECT_EQ(info.mBreps[0].mName, "\xc3\x96" "ffnung \xe2\x84\xa6");
}

TEST(CadJson, BadOptionsAndInconsistentInfo) {
    BrepInfo info;
    meshioplusplus::CadJsonReadOptions opts;
    opts.mSubdivisions = 65;
    EXPECT_THROW(meshioplusplus::read_cad_json(write_file(kModel), info, opts),
                 std::invalid_argument);
    const Mesh mesh = meshioplusplus::read_cad_json(write_file(kModel), info);
    info.mFaces[0].mPatch = 7;
    EXPECT_THROW(meshioplusplus::write_cad_json(mt::temp_path(".cad.json"), mesh, info),
                 meshioplusplus::WriteError);
    EXPECT_THROW(meshioplusplus::sample_brep(info), std::invalid_argument);
}

TEST(CadJson, MdpaSidecar) {
    namespace fs = std::filesystem;
    const fs::path dir = fs::path(mt::temp_path("_cad_project"));
    fs::remove_all(dir);  // a previous run's sidecar would be found first
    fs::create_directories(dir);
    const std::string deck = (dir / "model.mdpa").string();
    std::ofstream(deck, std::ios::binary)
        << "Begin Properties 0\nEnd Properties\n\nBegin Nodes\n 1 0 0 0\n 2 1 0 0\n 3 1 1 0\n"
           "End Nodes\n\nBegin Elements Element2D3N\n 1 0 1 2 3\nEnd Elements\n";
    EXPECT_EQ(meshioplusplus::mdpa_cad_json_path(deck), "");
    std::ofstream(dir / "geometry.cad.json", std::ios::binary) << kModel;
    std::ofstream(dir / "ProjectParameters.json", std::ios::binary)
        << R"({"modelers": [{"modeler_name": "CadIoModeler", "Parameters": {"geometry_file_name": "geometry"}}]})";
    EXPECT_EQ(meshioplusplus::mdpa_cad_json_path(deck), (dir / "geometry.cad.json").string());
    meshioplusplus::MdpaInfo info;
    BrepInfo brep;
    const Mesh mesh = meshioplusplus::read_mdpa(deck, info, brep);
    EXPECT_EQ(mesh.NumCellBlocks(), 1u);
    ASSERT_EQ(brep.mFaces.size(), 1u);
    const std::string out = (dir / "copy.mdpa").string();
    meshioplusplus::write_mdpa(out, mesh, info, brep);
    BrepInfo back;
    meshioplusplus::read_mdpa(out, info, back);
    EXPECT_EQ(back.mFaces.size(), 1u);
    EXPECT_TRUE(fs::exists(dir / "copy.cad.json"));
    EXPECT_THROW(meshioplusplus::read_mdpa(deck, info, brep, {}, (dir / "missing.cad.json").string()),
                 meshioplusplus::ReadError);
}

}  // namespace
