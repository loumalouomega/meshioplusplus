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
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/pipeline.hpp"
#include "meshioplusplus/registry.hpp"
#include "analytic_surfaces.hpp"
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

TEST(Pipeline, VersionTwoTypedMergeAndFanoutWorkWithoutJson) {
    using namespace meshioplusplus;
    const auto input = mt::temp_path("_v2_in.vtu");
    const auto extra = mt::temp_path("_v2_extra.vtu");
    const auto pattern = mt::temp_path("_piece_{part}.vtu");
    WriteOptions options;
    options.mCodecSet = true;
    options.mCodec = detail::VtkCodec::None;
    registry_write_ex(input, mt::tet_mesh(), "vtu", options);
    registry_write_ex(extra, mt::tet_mesh(), "vtu", options);
    Pipeline pipeline;
    pipeline.mVersion = 2;
    pipeline.mInput.mPath = input;
    pipeline.mOutput.mPath = pattern;
    pipeline.mOutput.mOptions = options;
    pipeline.mSteps = {
        step("Merge", {{"Inputs", std::vector<std::string>{extra}}}),
        step("Partition", {{"Nparts", std::int64_t(2)}, {"Method", std::string("sfc")}})};
    const auto report = run_pipeline(pipeline);
    ASSERT_EQ(report.mSteps.size(), 2u);
    EXPECT_EQ(report.mSteps[0].mCounters[0].second, 2.0);
    std::size_t cells = 0;
    for (const auto part : {"0", "1"}) {
        auto path = pattern;
        path.replace(path.find("{part}"), 6, part);
        auto piece = registry_read(path, "vtu", ReadOptions{});
        cells += total_cells(piece);
        std::filesystem::remove(path);
    }
    EXPECT_EQ(cells, 2 * total_cells(mt::tet_mesh()));
    std::filesystem::remove(input);
    std::filesystem::remove(extra);
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

TEST(Pipeline, SetsDataStepsAndReports) {
    PipelineReport report;
    Mesh in = mt::data_mesh();
    in.AddRegion(
        meshioplusplus::Region("a", meshioplusplus::RegionKind::Point, mt::int_data_array({0, 1})));
    in.AddRegion(
        meshioplusplus::Region("b", meshioplusplus::RegionKind::Point, mt::int_data_array({1, 2})));
    auto labels = meshioplusplus::apply_pipeline_step(
        std::move(in),
        step("SetsToData",
             {{"Location", std::string("point")}, {"Order", std::vector<std::string>{"b", "a"}}}),
        report);
    EXPECT_EQ(labels.PointData("b-a").As<std::int64_t>()[1], 1);
    EXPECT_EQ(labels.NumRegions(), 0u);
    auto restored = meshioplusplus::apply_pipeline_step(
        std::move(labels),
        step("DataToSets", {{"Location", std::string("point")}, {"Key", std::string("b-a")}}),
        report);
    EXPECT_EQ(restored.NumRegions(), 3u);
    EXPECT_FALSE(restored.HasPointData("b-a"));
    ASSERT_EQ(report.mSteps.size(), 2u);
    EXPECT_EQ(report.mSteps[0].mOp, "SetsToData");
    EXPECT_EQ(report.mSteps[1].mOp, "DataToSets");
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

// The counter `rName` of a step report, or nullptr when the step has none.
const double* pipeline_counter(const meshioplusplus::PipelineStepReport& rStep,
                               const std::string& rName) {
    for (const auto& [name, value] : rStep.mCounters)
        if (name == rName)
            return &value;
    return nullptr;
}

TEST(Pipeline, RefinePatchesStep) {
    // A feature angle above the box's 90 degrees leaves one patch no primitive
    // explains; RefinePatches splits the six sides back out of it.
    PipelineReport report;
    const Mesh out = meshioplusplus::run_pipeline_steps(
        analytic::box_skin(3),
        {step("Segment", {{"FeatureAngle", 100.0}}),
         step("RefinePatches", {{"FeatureAngle", 100.0}, {"MaxIterations", std::int64_t{20}}})},
        report);
    ASSERT_EQ(report.mSteps.size(), 2u);
    const auto& refined = report.mSteps[1];
    EXPECT_EQ(refined.mOp, "RefinePatches");
    for (const char* key : {"NumPatches", "NumIterations", "NumSplits", "NumMerges", "NumPrimitive",
                            "NumFreeform", "ExplainedAreaFraction", "HausdorffMax", "Stop"})
        EXPECT_NE(pipeline_counter(refined, key), nullptr) << key;
    EXPECT_EQ(*pipeline_counter(refined, "NumPatches"), 6.0);
    EXPECT_EQ(*pipeline_counter(refined, "NumPrimitive"), 6.0);
    EXPECT_GT(*pipeline_counter(refined, "NumSplits"), 0.0);
    EXPECT_NEAR(*pipeline_counter(refined, "ExplainedAreaFraction"), 1.0, 1e-12);
    const double stop = *pipeline_counter(refined, "Stop");
    EXPECT_GE(stop, 0.0);
    EXPECT_LE(stop, 3.0);
    EXPECT_TRUE(out.HasCellData("segment:patch"));
    EXPECT_TRUE(out.HasFieldData("refine:explained_area_fraction"));

    // The seed must be a JSON-exact integer: in [0, 2^53].
    for (const meshioplusplus::PipelineValue& bad :
         {meshioplusplus::PipelineValue{std::int64_t{-1}},
          meshioplusplus::PipelineValue{std::int64_t{9007199254740993}},
          meshioplusplus::PipelineValue{1.5}, meshioplusplus::PipelineValue{std::string("7")}}) {
        PipelineReport r;
        try {
            meshioplusplus::run_pipeline_steps(
                analytic::box_skin(2),
                {step("Segment"), step("RefinePatches", {{"RansacSeed", bad}})}, r);
            ADD_FAILURE() << "a bad RansacSeed was accepted";
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find("2^53"), std::string::npos) << e.what();
        }
    }
    PipelineReport ok;
    EXPECT_NO_THROW(meshioplusplus::run_pipeline_steps(
        analytic::box_skin(2),
        {step("Segment"), step("RefinePatches", {{"RansacSeed", std::int64_t{9007199254740992}}})},
        ok));
    // Without `segment:patch` there is nothing to refine.
    PipelineReport none;
    EXPECT_THROW(
        meshioplusplus::run_pipeline_steps(analytic::box_skin(2), {step("RefinePatches")}, none),
        std::invalid_argument);
}

TEST(Pipeline, FitCadRefineKeys) {
    // Only a run that refined reports it, so an existing report stays as it was.
    PipelineReport plain;
    meshioplusplus::run_pipeline_steps(
        analytic::box_skin(2), {step("FitCad", {{"Backend", std::string("native")}})}, plain);
    ASSERT_EQ(plain.mSteps.size(), 1u);
    EXPECT_EQ(pipeline_counter(plain.mSteps[0], "Refined"), nullptr);
    EXPECT_EQ(pipeline_counter(plain.mSteps[0], "ExplainedAreaFraction"), nullptr);

    PipelineReport refined;
    meshioplusplus::run_pipeline_steps(analytic::box_skin(2),
                                       {step("FitCad", {{"Backend", std::string("native")},
                                                        {"RefineSegmentation", true},
                                                        {"MaxRefineIterations", std::int64_t{5}},
                                                        {"RefineEpsilon", 1e-3},
                                                        {"NormalWeight", 1.0},
                                                        {"AllowMerge", true},
                                                        {"AllowSplit", true},
                                                        {"RansacSeed", std::int64_t{42}},
                                                        {"RefineTolerance", 1e-3},
                                                        {"RelativeTolerance", true}})},
                                       refined);
    ASSERT_EQ(refined.mSteps.size(), 1u);
    const auto& entry = refined.mSteps[0];
    ASSERT_NE(pipeline_counter(entry, "Refined"), nullptr);
    EXPECT_EQ(*pipeline_counter(entry, "Refined"), 1.0);
    EXPECT_NE(pipeline_counter(entry, "NumRefineSplits"), nullptr);
    EXPECT_NE(pipeline_counter(entry, "NumRefineMerges"), nullptr);
    ASSERT_NE(pipeline_counter(entry, "ExplainedAreaFraction"), nullptr);
    EXPECT_NEAR(*pipeline_counter(entry, "ExplainedAreaFraction"), 1.0, 1e-12);
    // The existing counters come first, in their old order.
    EXPECT_EQ(entry.mCounters.front().first, "Valid");
    ASSERT_NE(pipeline_counter(entry, "RefinementKept"), nullptr);
    EXPECT_EQ(entry.mCounters.size(), plain.mSteps[0].mCounters.size() + 5);

    PipelineReport bad;
    EXPECT_THROW(
        meshioplusplus::run_pipeline_steps(analytic::box_skin(2),
                                           {step("FitCad", {{"Backend", std::string("native")},
                                                            {"RansacSeed", std::int64_t{-3}}})},
                                           bad),
        std::invalid_argument);
}

TEST(Pipeline, FitCadPhaseThreeKeys) {
    // `Closed` is always reported; each Phase 3 option adds its counters only
    // when it is set, so an existing report gains just the one key.
    PipelineReport plain;
    meshioplusplus::run_pipeline_steps(
        analytic::box_skin(2), {step("FitCad", {{"Backend", std::string("native")}})}, plain);
    ASSERT_EQ(plain.mSteps.size(), 1u);
    ASSERT_NE(pipeline_counter(plain.mSteps[0], "Closed"), nullptr);
    EXPECT_EQ(*pipeline_counter(plain.mSteps[0], "Closed"), 1.0);
    EXPECT_EQ(plain.mSteps[0].mCounters.back().first, "Closed");
    for (const char* key :
         {"NumFreeformGeneral", "NumFreeformRefinements", "AvgFreeformPoles", "NumExtendedEdges",
          "NumSnappedJunctions", "HausdorffMax", "HausdorffRms"})
        EXPECT_EQ(pipeline_counter(plain.mSteps[0], key), nullptr) << key;

    PipelineReport all;
    meshioplusplus::run_pipeline_steps(analytic::box_skin(2),
                                       {step("FitCad", {{"Backend", std::string("native")},
                                                        {"FreeformGeneral", true},
                                                        {"FreeformTolerance", 0.0},
                                                        {"FreeformSmoothing", 1e-10},
                                                        {"MaxFreeformPoles", std::int64_t{20}},
                                                        {"MaxFreeformIterations", std::int64_t{5}},
                                                        {"MaxFreeformSplits", std::int64_t{1}},
                                                        {"SewTolerance", 0.0},
                                                        {"ExtendSurfaces", true},
                                                        {"SnapJunctions", true},
                                                        {"FaceDeviation", true}})},
                                       all);
    ASSERT_EQ(all.mSteps.size(), 1u);
    const auto& entry = all.mSteps[0];
    for (const char* key :
         {"NumFreeformGeneral", "NumFreeformRefinements", "AvgFreeformPoles", "NumExtendedEdges",
          "NumSnappedJunctions", "HausdorffMax", "HausdorffRms"})
        EXPECT_NE(pipeline_counter(entry, key), nullptr) << key;
    EXPECT_EQ(entry.mCounters.size(), plain.mSteps[0].mCounters.size() + 7);
    // The native backend builds no freeform faces and measures nothing.
    EXPECT_EQ(*pipeline_counter(entry, "NumFreeformGeneral"), 0.0);
    EXPECT_TRUE(std::isnan(*pipeline_counter(entry, "HausdorffMax")));

    PipelineReport bad;
    EXPECT_THROW(meshioplusplus::run_pipeline_steps(
                     analytic::box_skin(2),
                     {step("FitCad", {{"Backend", std::string("native")},
                                      {"MaxFreeformPoles", std::int64_t{3}}})},
                     bad),
                 std::invalid_argument);
    EXPECT_THROW(
        meshioplusplus::run_pipeline_steps(
            analytic::box_skin(2),
            {step("FitCad", {{"Backend", std::string("native")}, {"SewTolerance", -1.0}})}, bad),
        std::invalid_argument);
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

TEST(Pipeline, RemeshStepProducesANewMeshWithTheRequestedClusterCount) {
    // A regular octahedron -- the smallest closed 2-manifold triangle mesh.
    std::vector<std::vector<double>> pts = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                                            {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
    std::vector<std::vector<std::int64_t>> conn = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
                                                   {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};
    const Mesh octa = mt::make_mesh(pts, "triangle", conn);

    const std::string in_path = mt::temp_path("_pipe_remesh_in.vtk");
    const std::string out_path = mt::temp_path("_pipe_remesh_out.vtk");
    meshioplusplus::registry_writers().at("vtk")(in_path, octa);

    Pipeline pipeline;
    pipeline.mInput.mPath = in_path;
    pipeline.mSteps = {
        step("Remesh", {{"NumClusters", std::int64_t(30)}, {"Metric", std::string("quadric")}})};
    pipeline.mOutput.mPath = out_path;
    PipelineReport report = meshioplusplus::run_pipeline(pipeline);
    ASSERT_EQ(report.mSteps.size(), 1u);
    ASSERT_FALSE(report.mSteps[0].mCounters.empty());
    EXPECT_EQ(report.mSteps[0].mCounters[0].first, "NumClusters");
    EXPECT_EQ(report.mSteps[0].mCounters[0].second, 30.0);

    Mesh back = meshioplusplus::registry_read(out_path, "vtk", {});
    EXPECT_EQ(back.NumPoints(), 30u);
    std::filesystem::remove(in_path);
    std::filesystem::remove(out_path);
}

TEST(Pipeline, RemeshVolumeStepProducesATetraMesh) {
    std::vector<std::vector<double>> pts = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                                            {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
    std::vector<std::vector<std::int64_t>> conn = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
                                                   {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};
    const Mesh octa = mt::make_mesh(pts, "triangle", conn);

    const std::string in_path = mt::temp_path("_pipe_remesh_volume_in.vtk");
    const std::string out_path = mt::temp_path("_pipe_remesh_volume_out.vtk");
    meshioplusplus::registry_writers().at("vtk")(in_path, octa);

    Pipeline pipeline;
    pipeline.mInput.mPath = in_path;
    pipeline.mSteps = {
        step("RemeshVolume", {{"CellSize", 0.4}, {"WatertightCheck", std::string("off")}})};
    pipeline.mOutput.mPath = out_path;
    PipelineReport report = meshioplusplus::run_pipeline(pipeline);
    ASSERT_EQ(report.mSteps.size(), 1u);
    bool found_num_tets = false;
    for (const auto& counter : report.mSteps[0].mCounters)
        if (counter.first == "NumTets") {
            found_num_tets = true;
            EXPECT_GT(counter.second, 0.0);
        }
    EXPECT_TRUE(found_num_tets);

    Mesh back = meshioplusplus::registry_read(out_path, "vtk", {});
    EXPECT_GT(back.NumPoints(), 0u);
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).Type(), "tetra");
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
    EXPECT_EQ(meshioplusplus::pipeline_encoding_from_name("raw_appended"),
              meshioplusplus::WriteEncoding::RawAppended);
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
        {R"({"Version": 99, "Input": {"Path": "a"}, "Output": {"Path": "b"}})", "Version"},
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

TEST(PipelineJson, ReportSerializationEscapesStringsAndNullsNonFiniteCounters) {
    PipelineReport report;
    report.mSteps.push_back(
        {"quoted\"op", {{"count", 3.0}, {"invalid", std::numeric_limits<double>::quiet_NaN()}}});
    report.mWarnings = {"line\nwith\tcontrols\\and quotes\""};
    const auto text = meshioplusplus::pipeline_report_json(report);
    EXPECT_NE(text.find("quoted\\\"op"), std::string::npos);
    EXPECT_NE(text.find("\"invalid\":null"), std::string::npos);
    EXPECT_NE(text.find("line\\nwith\\tcontrols\\\\and quotes\\\""), std::string::npos);
    EXPECT_EQ(meshioplusplus::pipeline_report_json({}), "{\"steps\":[],\"warnings\":[]}");
}

TEST(PipelineJson, MissingSettingsFileFailsByPath) {
    try {
        meshioplusplus::parse_pipeline_file("/no/such/settings.json");
        FAIL() << "expected ReadError";
    } catch (const meshioplusplus::ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("/no/such/settings.json"), std::string::npos);
    }
}

TEST(PipelineJson, GenerateNeedsNoInput) {
    const std::string out_path = mt::temp_path("_pipe_generate.vtk");
    const std::string text = R"({
        "Version": 1,
        "Operations": [{"Op": "Generate", "Shape": "box", "Dims": [2, 2, 2],
                        "CellType": "tetra", "Order": 2}],
        "Output": {"Path": ")" + out_path + R"("}
    })";
    const PipelineReport report = meshioplusplus::run_pipeline_json(text);
    ASSERT_EQ(report.mSteps.size(), 1u);
    EXPECT_EQ(report.mSteps[0].mOp, "Generate");
    ASSERT_FALSE(report.mSteps[0].mCounters.empty());
    EXPECT_EQ(report.mSteps[0].mCounters[0].first, "NumPoints");
    EXPECT_TRUE(std::filesystem::exists(out_path));
    std::filesystem::remove(out_path);

    // Input stays required unless the first step generates the mesh.
    EXPECT_THROW(meshioplusplus::run_pipeline_json(R"({"Version": 1,
        "Operations": [{"Op": "Quality"}], "Output": {"Path": "x.vtk"}})"),
                 std::exception);
    EXPECT_THROW(meshioplusplus::run_pipeline_json(R"({"Version": 1,
        "Operations": [{"Op": "Generate", "Shape": "blob"}], "Output": {"Path": "x.vtk"}})"),
                 std::invalid_argument);
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
                    +[] { meshioplusplus::run_pipeline_file("x.json"); },
                    +[] { meshioplusplus::pipeline_report_json({}); }}) {
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
