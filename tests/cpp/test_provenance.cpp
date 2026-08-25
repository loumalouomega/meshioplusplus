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
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/write_options.hpp"
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

// Provenance is on by default: with nothing scoped a writer still renders the
// assumptions raised while it ran. `current_provenance_mode()` still reports
// Off -- it reports the *scope's* mode, and there is no scope -- but
// `default_provenance_mode()` is what a writer actually consults.
TEST(Provenance, NoScopeStillRecordsBecauseProvenanceIsOnByDefault) {
    EXPECT_EQ(meshioplusplus::detail::current_provenance_mode(), ProvenanceMode::Off);
    EXPECT_EQ(meshioplusplus::detail::default_provenance_mode(), ProvenanceMode::BestEffort);

    meshioplusplus::detail::provenance_begin_write();  // what a real write does
    meshioplusplus::detail::provenance_note("x", "y");
    ASSERT_EQ(meshioplusplus::detail::current_provenance().mNotes.size(), 1u);
    EXPECT_EQ(meshioplusplus::detail::current_provenance().mNotes[0].mCategory, "x");

    // The next write starts clean, so that note cannot attach itself to an
    // unrelated file -- the misattribution this bounding exists to prevent.
    meshioplusplus::detail::provenance_begin_write();
    EXPECT_TRUE(meshioplusplus::detail::current_provenance().mNotes.empty());
}

TEST(Provenance, DefaultModeCanBeTurnedOff) {
    const ProvenanceMode previous = meshioplusplus::detail::default_provenance_mode();
    meshioplusplus::detail::set_default_provenance_mode(ProvenanceMode::Off);
    meshioplusplus::detail::provenance_begin_write();
    meshioplusplus::detail::provenance_note("x", "y");
    const auto lines = meshioplusplus::detail::provenance_lines(SlotTier::Block);
    meshioplusplus::detail::set_default_provenance_mode(previous);

    ASSERT_EQ(lines.size(), 1u);
    EXPECT_EQ(lines[0], meshioplusplus::detail::kProvenanceTag);
}

// A note raised by an operation, with no write of its own, must NOT attach
// itself to the next unrelated file. Reproduced before the bounding existed:
// extract_surface(A) dropping a region put its note into mesh B's header.
TEST(Provenance, AnOperationsNoteDoesNotLeakIntoAnUnrelatedWrite) {
    meshioplusplus::detail::provenance_begin_write();
    meshioplusplus::detail::provenance_note("regions-dropped", "from an earlier operation");

    std::string path = mt::temp_path("_leak.obj");
    // Through the public write path, which is where the bounding lives.
    meshioplusplus::registry_write_ex(path, mt::tri_mesh(), "obj", meshioplusplus::WriteOptions{});
    auto meta = meshioplusplus::registry_read_metadata(path, "obj", meshioplusplus::ReadOptions{});

    for (const auto& l : meta.mProvenance)
        EXPECT_EQ(l.find("from an earlier operation"), std::string::npos)
            << "an earlier operation's note leaked into this file";
    std::error_code ec;
    std::filesystem::remove(path, ec);
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

// ---------------------------------------------------------------------------
// Read-back (v10.17.0): recovering a file's own block. Roadmap #1's last
// bullet -- see doc/provenance.md#reading-a-block-back.
// ---------------------------------------------------------------------------

TEST(Provenance, ScannerHandlesEverySlotShape) {
    // One scanner serves every format because render_block's CONTENT lines are
    // format-independent; only the wrapping differs. Each case below is a real
    // writer's wrapping.
    struct Case {
        const char* mLabel;
        const char* mText;
    };
    const Case cases[] = {
        {"hash", "# Written by meshio++ v1.2.3\n# Converted from in.vtu (vtu)\nv 0 0 0\n"},
        {"bang", "!PERMAS DataFile\n! Written by meshio++ v1.2.3\n$STRUCTURE\n"},
        {"xml-inline", "<?xml version=\"1.0\"?>\n<!--Written by meshio++ v1.2.3-->\n<VTKFile>\n"},
        {"xml-block", "<!--\nWritten by meshio++ v1.2.3\nTimestamp: X\n-->\n"},
        {"tecplot", "TITLE = \"Written by meshio++ v1.2.3\"\nVARIABLES = \"X\"\n"},
        {"ansys", "(1 \"Written by meshio++ v1.2.3\")\n(2 3)\n"},
        {"openfoam", "|   \\\\  /    A nd           | Written by meshio++ v1.2.3            |\n"},
    };
    for (const Case& c : cases) {
        auto r = meshioplusplus::detail::scan_provenance_text(c.mText);
        EXPECT_TRUE(r.mRecognised) << c.mLabel;
        ASSERT_FALSE(r.mLines.empty()) << c.mLabel;
        EXPECT_EQ(r.mLines[0], "Written by meshio++ v1.2.3") << c.mLabel;
    }
}

// The bug this pins: an unconditional trailing-')' strip (added for Ansys's
// `(1 "...")`) silently truncated every line that legitimately ends in one.
TEST(Provenance, ScannerKeepsParenthesesThatBelongToTheContent) {
    auto r = meshioplusplus::detail::scan_provenance_text(
        "# Written by meshio++ v1.2.3\n"
        "# Converted from in.vtu (vtu)\n"
        "# Operation: Clean(Weld=true)\n");
    ASSERT_EQ(r.mLines.size(), 3u);
    EXPECT_EQ(r.mLines[1], "Converted from in.vtu (vtu)");
    EXPECT_EQ(r.mLines[2], "Operation: Clean(Weld=true)");
}

TEST(Provenance, ScannerIsHonestAboutForeignAndAbsentBlocks) {
    EXPECT_FALSE(meshioplusplus::detail::scan_provenance_text("v 0 0 0\nf 1 2 3\n").mRecognised);
    auto foreign = meshioplusplus::detail::scan_provenance_text("# Created by SomeTool 3.2\n");
    EXPECT_FALSE(foreign.mRecognised);
    EXPECT_TRUE(foreign.mLines.empty());
    // A missing file is "nothing found", never a throw -- this enriches a
    // summary and must not be able to fail one.
    EXPECT_NO_THROW(meshioplusplus::detail::read_provenance_lines("/nonexistent/nope.obj"));
}

TEST(Provenance, DefaultWriteRecoversExactlyTheTag) {
    std::string path = mt::temp_path("_rb_off.obj");
    // A direct low-level writer call does not bound notes for itself (that is
    // the public write path's job), so do what a public write would.
    meshioplusplus::detail::provenance_begin_write();
    meshioplusplus::write_obj(path, mt::tri_mesh());
    auto meta = meshioplusplus::registry_read_metadata(path, "obj", meshioplusplus::ReadOptions{});
    EXPECT_TRUE(meta.mProvenanceRecognised);
    ASSERT_EQ(meta.mProvenance.size(), 1u);
    EXPECT_EQ(meta.mProvenance[0], meshioplusplus::detail::kProvenanceTag);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

TEST(Provenance, ScopedWriteRoundTripsTheWholeBlock) {
    std::string path = mt::temp_path("_rb_on.obj");
    {
        ProvenanceScope scope(ProvenanceMode::BestEffort);
        meshioplusplus::detail::provenance_set_source("in.msh", "gmsh");
        meshioplusplus::detail::provenance_note("regions-dropped", "Side regions dropped");
        meshioplusplus::write_obj(path, mt::tri_mesh());
    }
    auto meta = meshioplusplus::registry_read_metadata(path, "obj", meshioplusplus::ReadOptions{});
    EXPECT_TRUE(meta.mProvenanceRecognised);
    ASSERT_GE(meta.mProvenance.size(), 3u);
    EXPECT_EQ(meta.mProvenance[0], meshioplusplus::detail::kProvenanceTag);
    bool saw_source = false, saw_note = false;
    for (const auto& l : meta.mProvenance) {
        if (l.find("in.msh") != std::string::npos)
            saw_source = true;
        if (l.find("regions-dropped") != std::string::npos)
            saw_note = true;
    }
    EXPECT_TRUE(saw_source);
    EXPECT_TRUE(saw_note);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}

// The property that makes "replace, never append" structural rather than a
// rule a writer has to remember: nothing carries a read block into a write.
TEST(Provenance, AReadBlockIsNeverReEmitted) {
    std::string first = mt::temp_path("_rb_1.obj");
    {
        ProvenanceScope scope(ProvenanceMode::BestEffort);
        meshioplusplus::detail::provenance_set_source("original.msh", "gmsh");
        meshioplusplus::write_obj(first, mt::tri_mesh());
    }
    meshioplusplus::Mesh m = meshioplusplus::read_obj(first);

    std::string second = mt::temp_path("_rb_2.obj");
    {
        ProvenanceScope scope(ProvenanceMode::BestEffort);
        meshioplusplus::detail::provenance_set_source("second.msh", "gmsh");
        meshioplusplus::write_obj(second, m);
    }
    auto meta =
        meshioplusplus::registry_read_metadata(second, "obj", meshioplusplus::ReadOptions{});
    std::size_t tags = 0;
    bool saw_original = false;
    for (const auto& l : meta.mProvenance) {
        if (l.rfind("Written by meshio++ v", 0) == 0)
            ++tags;
        if (l.find("original.msh") != std::string::npos)
            saw_original = true;
    }
    EXPECT_EQ(tags, 1u) << "the block accumulated across a convert";
    EXPECT_FALSE(saw_original) << "the first file's source leaked into the second";
    std::error_code ec;
    std::filesystem::remove(first, ec);
    std::filesystem::remove(second, ec);
}

// OpenFOAM is the one slot with trailing structure (a fixed-width cell closed
// by '|'), so it gets its own end-to-end case rather than only the synthetic
// one above.
TEST(Provenance, OpenfoamBannerRoundTrips) {
    std::string dir = mt::temp_path("_rb_of");
    {
        ProvenanceScope scope(ProvenanceMode::BestEffort);
        meshioplusplus::OpenFoamInfo info;
        meshioplusplus::write_openfoam(dir, mt::hex_mesh(), info);
    }
    auto r = meshioplusplus::detail::read_provenance_lines(dir + "/constant/polyMesh/points");
    EXPECT_TRUE(r.mRecognised);
    ASSERT_FALSE(r.mLines.empty());
    EXPECT_EQ(r.mLines[0], meshioplusplus::detail::kProvenanceTag)
        << "the banner's padding or closing '|' leaked into the recovered line";
    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}
