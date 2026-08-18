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

// System includes
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/pipeline.hpp"
#include "meshioplusplus/registry.hpp"
#include "mesh_fixtures.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::Pipeline;
using meshioplusplus::PipelineReport;
using meshioplusplus::PipelineStep;

PipelineStep step(std::string op,
                  std::map<std::string, meshioplusplus::PipelineValue> params = {}) {
    PipelineStep s;
    s.mOp = std::move(op);
    s.mParams = std::move(params);
    return s;
}

std::size_t total_cells(const Mesh& rMesh) {
    std::size_t n = 0;
    for (const auto& block : rMesh.CellRange())
        n += block.NumCells();
    return n;
}

// ---------------------------------------------------------------------------
// The typed layer (always compiled, JSON-independent).
// ---------------------------------------------------------------------------

TEST(Pipeline, QualityAttachesMetrics) {
    PipelineReport report;
    Mesh out = meshioplusplus::apply_pipeline_step(mt::tet_mesh(), step("Quality"), report);
    EXPECT_TRUE(out.HasCellData("quality:scaled_jacobian"));
    ASSERT_EQ(report.mSteps.size(), 1u);
    EXPECT_EQ(report.mSteps[0].mOp, "Quality");
    EXPECT_TRUE(report.mSteps[0].mCounters.empty());
}

TEST(Pipeline, TransformRotatesAboutZ) {
    PipelineReport report;
    Mesh out = meshioplusplus::apply_pipeline_step(
        mt::tri_mesh(),
        step("Transform", {{"RotateAxis", std::vector<double>{0.0, 0.0, 1.0}},
                           {"RotateDegrees", std::int64_t(90)}}),
        report);
    // (1, 0, 0) -> (0, 1, 0) under a +90 deg rotation about z.
    const double* p = out.Points().As<double>();
    const std::size_t dim = out.PointDim();
    EXPECT_NEAR(p[1 * dim + 0], 0.0, 1e-12);
    EXPECT_NEAR(p[1 * dim + 1], 1.0, 1e-12);
}

TEST(Pipeline, TransformRequiresExactlyOneSource) {
    PipelineReport report;
    EXPECT_THROW(
        meshioplusplus::apply_pipeline_step(
            mt::tri_mesh(),
            step("Transform", {{"Translate", std::vector<double>{1.0, 0.0, 0.0}}, {"Scale", 2.0}}),
            report),
        std::invalid_argument);
    EXPECT_THROW(meshioplusplus::apply_pipeline_step(mt::tri_mesh(), step("Transform"), report),
                 std::invalid_argument);
}

TEST(Pipeline, ConvertCellsSimplexifiesAHex) {
    PipelineReport report;
    Mesh out = meshioplusplus::apply_pipeline_step(
        mt::hex_mesh(), step("ConvertCells", {{"Mode", std::string("simplexify")}}), report);
    ASSERT_EQ(out.NumCellBlocks(), 1u);
    EXPECT_EQ(out.Cells(0).Type(), "tetra");
    EXPECT_EQ(out.Cells(0).NumCells(), 6u);
}

TEST(Pipeline, CleanReportsItsCounters) {
    PipelineReport report;
    Mesh out = meshioplusplus::apply_pipeline_step(mt::tet_mesh(), step("Clean"), report);
    ASSERT_EQ(report.mSteps.size(), 1u);
    const auto& counters = report.mSteps[0].mCounters;
    ASSERT_EQ(counters.size(), 4u);
    EXPECT_EQ(counters[0].first, "PointsWelded");
    EXPECT_EQ(counters[1].first, "PointsRemovedOrphan");
    EXPECT_EQ(counters[2].first, "CellsDroppedDegenerate");
    EXPECT_EQ(counters[3].first, "CellsDroppedDuplicate");
    EXPECT_EQ(total_cells(out), 2u);
}

TEST(Pipeline, PartitionAttachesLabels) {
    PipelineReport report;
    Mesh out = meshioplusplus::apply_pipeline_step(
        mt::tet_mesh(), step("Partition", {{"Nparts", std::int64_t(2)}}), report);
    EXPECT_TRUE(out.HasCellData("partition:part"));
    ASSERT_EQ(report.mSteps.size(), 1u);
    ASSERT_EQ(report.mSteps[0].mCounters.size(), 1u);
    EXPECT_EQ(report.mSteps[0].mCounters[0].first, "Nparts");
    EXPECT_EQ(report.mSteps[0].mCounters[0].second, 2.0);
}

TEST(Pipeline, SliceWarnsWhenThePlaneMissesTheMesh) {
    PipelineReport report;
    Mesh out = meshioplusplus::apply_pipeline_step(
        mt::tet_mesh(),
        step("Slice", {{"Point", std::vector<double>{0.0, 0.0, 100.0}},
                       {"Normal", std::vector<double>{0.0, 0.0, 1.0}}}),
        report);
    EXPECT_EQ(total_cells(out), 0u);
    ASSERT_EQ(report.mWarnings.size(), 1u);
    EXPECT_NE(report.mWarnings[0].find("the section is empty"), std::string::npos);
}

TEST(Pipeline, SectionIsAnAliasOfSlice) {
    PipelineReport ra, rb;
    Mesh a = meshioplusplus::apply_pipeline_step(
        mt::tet_mesh(),
        step("Slice", {{"Point", std::vector<double>{0.0, 0.0, 0.25}},
                       {"Normal", std::vector<double>{0.0, 0.0, 1.0}}}),
        ra);
    Mesh b = meshioplusplus::apply_pipeline_step(
        mt::tet_mesh(),
        step("Section", {{"Point", std::vector<double>{0.0, 0.0, 0.25}},
                         {"Normal", std::vector<double>{0.0, 0.0, 1.0}}}),
        rb);
    EXPECT_EQ(total_cells(a), total_cells(b));
    EXPECT_GT(total_cells(a), 0u);
    EXPECT_EQ(rb.mSteps.at(0).mOp, "Section");
}

TEST(Pipeline, CropByBbox) {
    PipelineReport report;
    Mesh out = meshioplusplus::apply_pipeline_step(
        mt::quad_mesh(),
        step("Crop", {{"Bbox", std::vector<double>{-0.5, -0.5, -0.5, 1.5, 1.5, 0.5}}}), report);
    EXPECT_EQ(total_cells(out), 1u);
    ASSERT_EQ(report.mSteps.at(0).mCounters.size(), 1u);
    EXPECT_EQ(report.mSteps.at(0).mCounters[0].first, "CellsKept");
}

TEST(Pipeline, CropRejectsBothOrNeitherRegion) {
    PipelineReport report;
    EXPECT_THROW(meshioplusplus::apply_pipeline_step(mt::quad_mesh(), step("Crop"), report),
                 std::invalid_argument);
    EXPECT_THROW(meshioplusplus::apply_pipeline_step(
                     mt::quad_mesh(),
                     step("Crop", {{"Bbox", std::vector<double>{0, 0, 0, 1, 1, 1}},
                                   {"Point", std::vector<double>{0, 0, 0}},
                                   {"Normal", std::vector<double>{0, 0, 1}}}),
                     report),
                 std::invalid_argument);
}

TEST(Pipeline, DataCalcSplitsOnTheFirstEquals) {
    Mesh in = mt::tri_mesh();
    {
        meshioplusplus::NDArray a(meshioplusplus::DType::Float64, {4});
        for (int i = 0; i < 4; ++i)
            a.As<double>()[i] = i;
        in.AddPointData("a", std::move(a));
    }
    PipelineReport report;
    // move `in` in: apply_pipeline_step takes Mesh by value, and KRATOS's
    // Mesh is not copy-constructible.
    Mesh out = meshioplusplus::apply_pipeline_step(
        std::move(in), step("DataCalc", {{"Expr", std::string("b = a + 1")}}), report);
    EXPECT_TRUE(out.HasPointData("b"));
}

TEST(Pipeline, DataRenameSplitsOnTheLastColon) {
    Mesh in = mt::tri_mesh();
    {
        meshioplusplus::NDArray a(meshioplusplus::DType::Float64, {4});
        for (int i = 0; i < 4; ++i)
            a.As<double>()[i] = i;
        in.AddPointData("ns:old", std::move(a));
    }
    PipelineReport report;
    // The LAST colon separates OLD from NEW (the CLI's rule: OLD may carry
    // colons -- `gmsh:physical` -- NEW may not), so "ns:old:new" renames the
    // array "ns:old" to "new".
    // move `in` in: apply_pipeline_step takes Mesh by value, and KRATOS's
    // Mesh is not copy-constructible.
    Mesh out = meshioplusplus::apply_pipeline_step(
        std::move(in), step("DataRename", {{"Point", std::vector<std::string>{"ns:old:new"}}}),
        report);
    EXPECT_FALSE(out.HasPointData("ns:old"));
    EXPECT_TRUE(out.HasPointData("new"));
}

TEST(Pipeline, RunPipelineStepsChains) {
    PipelineReport report;
    Mesh out = meshioplusplus::run_pipeline_steps(
        mt::hex_mesh(),
        {step("ConvertCells", {{"Mode", std::string("simplexify")}}), step("Quality"),
         step("Clean")},
        report);
    EXPECT_EQ(report.mSteps.size(), 3u);
    EXPECT_TRUE(out.HasCellData("quality:scaled_jacobian"));
    EXPECT_EQ(out.Cells(0).Type(), "tetra");
}

TEST(Pipeline, ValidateRejectsUnknownOpAndKey) {
    try {
        meshioplusplus::validate_pipeline_step(step("Nope"));
        FAIL() << "expected invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("unknown operation 'Nope'"), std::string::npos);
        EXPECT_NE(std::string(e.what()).find("Quality"), std::string::npos);  // lists known ops
    }
    try {
        meshioplusplus::validate_pipeline_step(step("Clean", {{"Foo", 1.0}}));
        FAIL() << "expected invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("unknown parameter 'Foo'"), std::string::npos);
        EXPECT_NE(std::string(e.what()).find("RemoveOrphans"), std::string::npos);
    }
}

TEST(Pipeline, ExcludedOpsNameTheCliVerb) {
    for (const char* op : {"Merge", "Interpolate", "Split", "Diff", "UndoGreen"}) {
        try {
            meshioplusplus::validate_pipeline_step(step(op));
            FAIL() << "expected invalid_argument for " << op;
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find("CLI verb"), std::string::npos) << op;
        }
    }
}

TEST(Pipeline, MistypedParametersFailByName) {
    PipelineReport report;
    try {
        meshioplusplus::apply_pipeline_step(
            mt::tet_mesh(), step("Clean", {{"Atol", std::string("tight")}}), report);
        FAIL() << "expected invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("'Atol' must be a number"), std::string::npos);
    }
    try {
        meshioplusplus::apply_pipeline_step(mt::tet_mesh(), step("Clean", {{"Weld", 1.0}}), report);
        FAIL() << "expected invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("'Weld' must be a boolean"), std::string::npos);
    }
}

TEST(Pipeline, OpTableExposesTheVocabulary) {
    const auto table = meshioplusplus::pipeline_op_table();
    bool has_slice = false, has_section = false, has_transform = false;
    for (const auto& entry : table) {
        if (entry.first == "Slice")
            has_slice = true;
        if (entry.first == "Section")
            has_section = true;
        if (entry.first == "Transform") {
            has_transform = true;
            EXPECT_NE(std::find(entry.second.begin(), entry.second.end(), "RotateAxis"),
                      entry.second.end());
        }
    }
    EXPECT_TRUE(has_slice);
    EXPECT_TRUE(has_section);
    EXPECT_TRUE(has_transform);
}

TEST(Pipeline, RunPipelineEndToEndTyped) {
    const std::string in_path = mt::temp_path("_pipe_in.vtk");
    const std::string out_path = mt::temp_path("_pipe_out.vtk");
    meshioplusplus::registry_writers().at("vtk")(in_path, mt::tet_mesh());

    Pipeline pipeline;
    pipeline.mInput.mPath = in_path;
    pipeline.mSteps = {step("Quality"), step("Clean")};
    pipeline.mOutput.mPath = out_path;
    PipelineReport report = meshioplusplus::run_pipeline(pipeline);
    EXPECT_EQ(report.mSteps.size(), 2u);
    ASSERT_TRUE(std::filesystem::exists(out_path));

    Mesh back = meshioplusplus::registry_read(out_path, "vtk", {});
    EXPECT_TRUE(back.HasCellData("quality:scaled_jacobian"));
    std::filesystem::remove(in_path);
    std::filesystem::remove(out_path);
}

TEST(Pipeline, RunPipelineValidatesBeforeReading) {
    Pipeline pipeline;
    pipeline.mInput.mPath = "/definitely/not/a/real/file.vtk";
    pipeline.mSteps = {step("Nope")};
    pipeline.mOutput.mPath = "/tmp/never_written.vtk";
    // The bad step must fail first -- the input is never opened, so the error
    // is invalid_argument (schema), not ReadError (missing file).
    EXPECT_THROW(meshioplusplus::run_pipeline(pipeline), std::invalid_argument);
}

TEST(Pipeline, OptionSpellingsAreOwnedHere) {
    EXPECT_EQ(meshioplusplus::pipeline_encoding_from_name(""),
              meshioplusplus::WriteEncoding::Default);
    EXPECT_EQ(meshioplusplus::pipeline_encoding_from_name("ascii"),
              meshioplusplus::WriteEncoding::Ascii);
    EXPECT_EQ(meshioplusplus::pipeline_encoding_from_name("binary"),
              meshioplusplus::WriteEncoding::Binary);
    EXPECT_THROW(meshioplusplus::pipeline_encoding_from_name("base64"), std::invalid_argument);
    EXPECT_EQ(meshioplusplus::pipeline_codec_from_name("zstd"),
              meshioplusplus::detail::VtkCodec::ZSTD);
    EXPECT_THROW(meshioplusplus::pipeline_codec_from_name(""), std::invalid_argument);
    EXPECT_EQ(meshioplusplus::pipeline_mmap_from_name("off"), meshioplusplus::MmapMode::Off);
    EXPECT_THROW(meshioplusplus::pipeline_mmap_from_name("maybe"), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// The JSON front-end. Split per build: with the parser, parse/run round
// trips; without it, every entry point throws NAMING THE FLAG (and this leg
// runs on the -DMESHIOPLUSPLUS_WITH_JSON=OFF CI matrix entry, so the gate is
// verified to fire, not just written).
// ---------------------------------------------------------------------------

#ifdef MESHIOPLUSPLUS_HAS_JSON

TEST(PipelineJson, ReportsTheParser) {
    EXPECT_TRUE(meshioplusplus::pipeline_has_json());
}

TEST(PipelineJson, ParsesAFullDocument) {
    const std::string text = R"({
        "Version": 1,
        "Input": {"Path": "in.msh", "Format": "gmsh",
                  "Options": {"PointsOnly": true, "TimeStep": 2,
                              "DataArrays": ["u", "v"], "Mmap": "off",
                              "Lenient": true}},
        "Operations": [
            {"Op": "Transform", "RotateAxis": [0, 0, 1], "RotateDegrees": 45},
            {"Op": "Clean", "RemoveOrphans": false}
        ],
        "Output": {"Path": "out.vtu", "Encoding": "binary", "Codec": "zstd",
                   "FloatFormat": ".9e"}
    })";
    Pipeline p = meshioplusplus::parse_pipeline_json(text);
    EXPECT_EQ(p.mVersion, 1);
    EXPECT_EQ(p.mInput.mPath, "in.msh");
    EXPECT_EQ(p.mInput.mFormat, "gmsh");
    EXPECT_TRUE(p.mInput.mOptions.mPointsOnly);
    EXPECT_TRUE(p.mInput.mOptions.mLenient);
    EXPECT_EQ(p.mInput.mOptions.mTimeStep, 2);
    ASSERT_TRUE(p.mInput.mOptions.mDataArrays.has_value());
    EXPECT_EQ(p.mInput.mOptions.mDataArrays->size(), 2u);
    EXPECT_EQ(p.mInput.mOptions.mMmap, meshioplusplus::MmapMode::Off);
    ASSERT_EQ(p.mSteps.size(), 2u);
    EXPECT_EQ(p.mSteps[0].mOp, "Transform");
    EXPECT_EQ(p.mSteps[1].mOp, "Clean");
    EXPECT_EQ(p.mOutput.mPath, "out.vtu");
    EXPECT_EQ(p.mOutput.mOptions.mEncoding, meshioplusplus::WriteEncoding::Binary);
    EXPECT_TRUE(p.mOutput.mOptions.mCodecSet);
    EXPECT_EQ(p.mOutput.mOptions.mCodec, meshioplusplus::detail::VtkCodec::ZSTD);
    EXPECT_EQ(p.mOutput.mOptions.mFloatFormat, ".9e");
}

TEST(PipelineJson, StrictSchemaErrorsNameTheOffender) {
    const char* cases[][2] = {
        {R"({"Input": {"Path": "a"}, "Output": {"Path": "b"}, "Bogus": 1})", "Bogus"},
        {R"({"Version": 2, "Input": {"Path": "a"}, "Output": {"Path": "b"}})", "Version"},
        {R"({"Output": {"Path": "b"}})", "Input is required"},
        {R"({"Input": {"Path": "a"}, "Output": {}})", "Output.Path"},
        {R"({"Input": {"Path": "a"}, "Output": {"Path": "b"},
             "Operations": [{"Op": "Nope"}]})",
         "unknown operation"},
        {R"({"Input": {"Path": "a"}, "Output": {"Path": "b"},
             "Operations": [{"Op": "Clean", "Foo": 1}]})",
         "unknown parameter 'Foo'"},
        {R"({"Input": {"Path": "a"}, "Output": {"Path": "b", "Codec": "brotli"}})", "Codec"},
        {R"({"Input": {"Path": "a"}, "Output": {"Path": "b"},
             "Operations": [{"Op": "Refine", "Cells": [1, "x"]}]})",
         "homogeneous"},
    };
    for (const auto& c : cases) {
        try {
            meshioplusplus::parse_pipeline_json(c[0]);
            FAIL() << "expected invalid_argument for: " << c[0];
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find(c[1]), std::string::npos)
                << e.what() << " should mention " << c[1];
        }
    }
    EXPECT_THROW(meshioplusplus::parse_pipeline_json("not json at all"), std::invalid_argument);
}

TEST(PipelineJson, RunPipelineFileEndToEnd) {
    const std::string in_path = mt::temp_path("_pipe_json_in.vtk");
    const std::string out_path = mt::temp_path("_pipe_json_out.vtk");
    const std::string settings_path = mt::temp_path("_pipe_settings.json");
    meshioplusplus::registry_writers().at("vtk")(in_path, mt::hex_mesh());

    {
        std::ofstream settings(settings_path);
        settings << R"({
            "Input": {"Path": ")"
                 << in_path << R"("},
            "Operations": [
                {"Op": "ConvertCells", "Mode": "simplexify"},
                {"Op": "Quality"}
            ],
            "Output": {"Path": ")"
                 << out_path << R"("}
        })";
    }
    PipelineReport report = meshioplusplus::run_pipeline_file(settings_path);
    EXPECT_EQ(report.mSteps.size(), 2u);
    Mesh back = meshioplusplus::registry_read(out_path, "vtk", {});
    EXPECT_EQ(back.Cells(0).Type(), "tetra");
    EXPECT_TRUE(back.HasCellData("quality:scaled_jacobian"));
    std::filesystem::remove(in_path);
    std::filesystem::remove(out_path);
    std::filesystem::remove(settings_path);
}

TEST(PipelineJson, MissingSettingsFileFailsByPath) {
    try {
        meshioplusplus::parse_pipeline_file("/no/such/settings.json");
        FAIL() << "expected ReadError";
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("/no/such/settings.json"), std::string::npos);
    }
}

#else  // !MESHIOPLUSPLUS_HAS_JSON

TEST(PipelineJson, ReportsNoParser) {
    EXPECT_FALSE(meshioplusplus::pipeline_has_json());
}

TEST(PipelineJson, EntryPointsThrowNamingTheFlag) {
    // The typed layer stays fully functional; only the JSON front-end is
    // compiled out, and it must fail naming the option -- never a link error,
    // never a silent no-op.
    for (auto fn : {+[] { meshioplusplus::parse_pipeline_json("{}"); },
                    +[] { meshioplusplus::parse_pipeline_file("x.json"); },
                    +[] { meshioplusplus::run_pipeline_json("{}"); },
                    +[] { meshioplusplus::run_pipeline_file("x.json"); }}) {
        try {
            fn();
            FAIL() << "expected runtime_error";
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("MESHIOPLUSPLUS_WITH_JSON"), std::string::npos);
        }
    }
}

#endif  // MESHIOPLUSPLUS_HAS_JSON

}  // namespace
