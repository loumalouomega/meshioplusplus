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
#include <fstream>
#include <sstream>
#include <string>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/region_remap.hpp"
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/operations/pipeline.hpp"
#include "meshioplusplus/formats/openfoam.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/vtu.hpp"

// Asserts the one-line provenance credit every writer emits (see
// detail/provenance.hpp and doc/roadmap.md's "audit and normalize" bullet)
// appears in the *written bytes* -- never inferred by reading the file back
// through our own reader, which cannot tell a right tag from a wrong one
// (the CGNS-ordering-suite rationale: a self-consistent round trip is not an
// oracle for the writer alone).

namespace {

std::string slurp(const std::string& rPath) {
    std::ifstream f(rPath, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

}  // namespace

TEST(Provenance, ObjCarriesTheTag) {
    std::string path = mt::temp_path(".obj");
    meshioplusplus::write_obj(path, mt::tri_mesh());
    std::string bytes = slurp(path);
    EXPECT_NE(bytes.find(meshioplusplus::detail::kProvenanceTag), std::string::npos);
    // The tag is the very first content line.
    EXPECT_EQ(bytes.rfind(std::string("# ") + meshioplusplus::detail::kProvenanceTag, 0), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Provenance, OffCarriesTheTag) {
    std::string path = mt::temp_path(".off");
    meshioplusplus::write_off(path, mt::tri_mesh());
    std::string bytes = slurp(path);
    std::string expected = std::string("OFF\n# ") + meshioplusplus::detail::kProvenanceTag;
    EXPECT_EQ(bytes.rfind(expected, 0), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Provenance, VtuCarriesTheTagAsAnXmlComment) {
    std::string path = mt::temp_path(".vtu");
    meshioplusplus::write_vtu(path, mt::tri_mesh(), /*binary=*/false, /*zlib=*/false);
    std::string bytes = slurp(path);
    std::string expected = std::string("<!--") + meshioplusplus::detail::kProvenanceTag + "-->";
    EXPECT_NE(bytes.find(expected), std::string::npos);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// The binary STL header is a fixed 80-byte slot: the tag, then NUL padding,
// with no trailing newline -- unlike every other writer here.
TEST(Provenance, StlBinaryHeaderIsTagThenNulPadded) {
    std::string path = mt::temp_path(".stl");
    meshioplusplus::write_stl(path, mt::tri_mesh(), /*binary=*/true, /*skin=*/true);
    std::ifstream f(path, std::ios::binary);
    char header[80];
    f.read(header, 80);
    ASSERT_EQ(f.gcount(), 80);
    std::string tag = meshioplusplus::detail::kProvenanceTag;
    EXPECT_EQ(std::string(header, tag.size()), tag);
    for (std::size_t i = tag.size(); i < 80; ++i)
        EXPECT_EQ(header[i], '\0') << "byte " << i << " of the STL header is not NUL-padded";
    // Must not start with "solid" -- that is the ASCII/binary sniff heuristic
    // both the C++ and Python STL readers use.
    EXPECT_NE(tag.rfind("solid", 0), 0u);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// nastran is the one documented exception: the C++ reader gates on a
// sentinel comment line the Python writer never emits (doc/formats/nastran.md),
// so the sentinel stays first and the tag follows it as a second `$` line --
// both must be present, in that order, for the reader's own gate to keep
// working exactly as before this change.
TEST(Provenance, NastranSentinelPrecedesTheTag) {
    std::string path = mt::temp_path(".bdf");
    meshioplusplus::write_nastran(path, mt::tri_mesh());
    std::string bytes = slurp(path);
    auto sentinel_pos = bytes.find("meshioplusplus-cpp-nastran");
    auto tag_pos = bytes.find(meshioplusplus::detail::kProvenanceTag);
    ASSERT_NE(sentinel_pos, std::string::npos);
    ASSERT_NE(tag_pos, std::string::npos);
    EXPECT_LT(sentinel_pos, tag_pos);
    // The reader must still accept its own output.
    EXPECT_NO_THROW(meshioplusplus::read_nastran(path));
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// The OpenFOAM banner box is fixed-width; the tag must not push the
// box's closing "|" out of column.
TEST(Provenance, OpenfoamBannerStaysColumnAligned) {
    std::string dir = mt::temp_path("_of");
    meshioplusplus::OpenFoamInfo info;
    meshioplusplus::write_openfoam(dir, mt::hex_mesh(), info);
    std::string points_path = dir + "/constant/polyMesh/points";
    std::string bytes = slurp(points_path);
    auto tag_pos = bytes.find(meshioplusplus::detail::kProvenanceTag);
    ASSERT_NE(tag_pos, std::string::npos);
    auto line_start = bytes.rfind('\n', tag_pos) + 1;
    auto line_end = bytes.find('\n', tag_pos);
    std::string line = bytes.substr(line_start, line_end - line_start);
    // Every banner line (the box border, blank cells, and this one) is the
    // same fixed width, box borders included.
    EXPECT_EQ(line.size(), 79u) << "banner line: " << line;
    EXPECT_EQ(line.front(), '|');
    EXPECT_EQ(line.back(), '|');
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

// ---------------------------------------------------------------------------
// The opt-in record (v10.16.0): Mode/SlotTier semantics, the operation-chain
// bridge, and conversion-assumption capture -- roadmap #1's bullets 2-7.
// ---------------------------------------------------------------------------

namespace {
using meshioplusplus::detail::ProvenanceMode;
using meshioplusplus::detail::ProvenanceRecord;
using meshioplusplus::detail::ProvenanceScope;
using meshioplusplus::detail::SlotTier;

/// An Int64 `(n,)` entry array -- the `test_region_api.cpp` helper,
/// duplicated locally rather than shared, since it is three lines.
meshioplusplus::NDArray region_entries(const std::vector<std::int64_t>& rVals) {
    meshioplusplus::NDArray a =
        meshioplusplus::NDArray::Uninit(meshioplusplus::DType::Int64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<std::int64_t>()[i] = rVals[i];
    return a;
}

}  // namespace

TEST(Provenance, OffModeIgnoresEverySlotTier) {
    // No scope open -- exactly v10.15.0's contract, for every tier.
    for (SlotTier tier :
         {SlotTier::None, SlotTier::Bounded, SlotTier::SingleLine, SlotTier::Block}) {
        auto lines = meshioplusplus::detail::provenance_lines(tier);
        ASSERT_EQ(lines.size(), 1u);
        EXPECT_EQ(lines[0], meshioplusplus::detail::kProvenanceTag);
    }
}

TEST(Provenance, BestEffortRendersTheFullBlockOnlyForBlockTier) {
    ProvenanceScope scope(ProvenanceMode::BestEffort);
    meshioplusplus::detail::provenance_set_source("in.vtu", "vtu");
    meshioplusplus::detail::provenance_add_operation("clean(weld=true)");
    meshioplusplus::detail::provenance_note("regions-dropped", "3 regions dropped");

    auto block = meshioplusplus::detail::provenance_lines(SlotTier::Block);
    ASSERT_EQ(block.size(), 5u);  // tag, source, operation, note, timestamp
    EXPECT_EQ(block[0], meshioplusplus::detail::kProvenanceTag);
    EXPECT_NE(block[1].find("in.vtu"), std::string::npos);
    EXPECT_NE(block[2].find("clean(weld=true)"), std::string::npos);
    EXPECT_NE(block[3].find("3 regions dropped"), std::string::npos);

    // SingleLine/Bounded degrade to the tag alone, silently, even though a
    // record is active -- the slot is the honest maximum, not a preference.
    EXPECT_EQ(meshioplusplus::detail::provenance_lines(SlotTier::SingleLine).size(), 1u);
    EXPECT_EQ(meshioplusplus::detail::provenance_lines(SlotTier::Bounded).size(), 1u);
}

TEST(Provenance, RequiredThrowsOnlyForNoneTier) {
    ProvenanceScope scope(ProvenanceMode::Required);
    EXPECT_THROW(meshioplusplus::detail::provenance_lines(SlotTier::None), std::exception);
    // Degrading to the tag on a structurally smaller slot is not a failure.
    EXPECT_NO_THROW(meshioplusplus::detail::provenance_lines(SlotTier::SingleLine));
    EXPECT_NO_THROW(meshioplusplus::detail::provenance_lines(SlotTier::Bounded));
    EXPECT_NO_THROW(meshioplusplus::detail::provenance_lines(SlotTier::Block));
}

TEST(Provenance, DuplicateNotesAreCollapsed) {
    ProvenanceScope scope(ProvenanceMode::BestEffort);
    meshioplusplus::detail::provenance_note("dtype", "cast to int32");
    meshioplusplus::detail::provenance_note("dtype", "cast to int32");
    meshioplusplus::detail::provenance_note("dtype", "cast to int32");
    EXPECT_EQ(scope.Get().mNotes.size(), 1u);
}

TEST(Provenance, ScopesNestAndRestore) {
    ProvenanceScope outer(ProvenanceMode::BestEffort);
    meshioplusplus::detail::provenance_set_source("outer.vtu", "vtu");
    {
        ProvenanceScope inner(ProvenanceMode::BestEffort);
        meshioplusplus::detail::provenance_set_source("inner.vtu", "vtu");
        EXPECT_EQ(meshioplusplus::detail::current_provenance().mSourcePath, "inner.vtu");
    }
    EXPECT_EQ(meshioplusplus::detail::current_provenance().mSourcePath, "outer.vtu");
}

TEST(Provenance, NoScopeMeansOffAndNoteIsANoOp) {
    EXPECT_EQ(meshioplusplus::detail::current_provenance_mode(), ProvenanceMode::Off);
    meshioplusplus::detail::provenance_note("x", "y");  // must not throw or crash
    EXPECT_TRUE(meshioplusplus::detail::current_provenance().mNotes.empty());
}

TEST(Provenance, TimestampHonoursSourceDateEpochAndTheOffSwitch) {
    setenv("SOURCE_DATE_EPOCH", "1000000000", 1);
    EXPECT_EQ(meshioplusplus::detail::provenance_timestamp(), "2001-09-09T01:46:40Z");
    unsetenv("SOURCE_DATE_EPOCH");

    setenv("MESHIOPLUSPLUS_PROVENANCE_TIMESTAMP", "off", 1);
    EXPECT_EQ(meshioplusplus::detail::provenance_timestamp(), "");
    unsetenv("MESHIOPLUSPLUS_PROVENANCE_TIMESTAMP");
}

TEST(Provenance, PipelineRecordsSourceTargetAndOperations) {
    std::string in_path = mt::temp_path(".obj");
    meshioplusplus::write_obj(in_path, mt::tri_mesh());
    std::string out_path = mt::temp_path("_out.obj");

    meshioplusplus::Pipeline p;
    p.mInput.mPath = in_path;
    meshioplusplus::PipelineStep step;
    step.mOp = "Clean";
    step.mParams["Weld"] = true;
    p.mSteps.push_back(step);
    p.mOutput.mPath = out_path;

    ProvenanceScope scope(ProvenanceMode::BestEffort);
    meshioplusplus::run_pipeline(p);
    const auto& rec = scope.Get();
    EXPECT_EQ(rec.mSourcePath, in_path);
    EXPECT_EQ(rec.mSourceFormat, "obj");
    EXPECT_EQ(rec.mTargetFormat, "obj");
    ASSERT_EQ(rec.mOperations.size(), 1u);
    EXPECT_NE(rec.mOperations[0].find("Clean("), std::string::npos);
    EXPECT_NE(rec.mOperations[0].find("Weld=true"), std::string::npos);

    std::error_code ec;
    std::filesystem::remove(in_path, ec);
    std::filesystem::remove(out_path, ec);
}

TEST(Provenance, WarnRegionsDroppedRecordsANote) {
    meshioplusplus::Mesh mesh = mt::tri_mesh();
    mesh.AddRegion(meshioplusplus::Region("boundary", meshioplusplus::RegionKind::Point,
                                          region_entries({0, 1})));

    ProvenanceScope scope(ProvenanceMode::BestEffort);
    meshioplusplus::detail::warn_regions_dropped(mesh, "isosurface");
    ASSERT_EQ(scope.Get().mNotes.size(), 1u);
    EXPECT_EQ(scope.Get().mNotes[0].mCategory, "regions-dropped");
    EXPECT_NE(scope.Get().mNotes[0].mDetail.find("isosurface"), std::string::npos);
}

TEST(Provenance, OffWriterRecordsDroppedCellTypes) {
    meshioplusplus::Mesh mesh = mt::tri_quad_mesh();
    meshioplusplus::NDArray tet_conn =
        meshioplusplus::NDArray::Uninit(meshioplusplus::DType::Int64, {1, 4});
    for (int i = 0; i < 4; ++i)
        tet_conn.As<std::int64_t>()[i] = i % static_cast<int>(mesh.NumPoints());
    mesh.AddCellBlock("tetra", tet_conn);

    ProvenanceScope scope(ProvenanceMode::BestEffort);
    std::string path = mt::temp_path(".off");
    meshioplusplus::write_off(path, mesh);
    bool found = false;
    for (const auto& n : scope.Get().mNotes)
        if (n.mCategory == "cells-dropped" && n.mDetail.find("tetra") != std::string::npos)
            found = true;
    EXPECT_TRUE(found);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
