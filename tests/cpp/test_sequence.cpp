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
#include <atomic>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <fstream>
#include <system_error>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/operations/sequence.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::ReadOptions;
using meshioplusplus::SequenceEntry;
using meshioplusplus::SequenceInput;
using meshioplusplus::SequenceMode;
using meshioplusplus::SequenceOutput;
using meshioplusplus::SequencePipeline;
using meshioplusplus::SequenceTimeFrom;
using meshioplusplus::SequenceTimeSource;

namespace {

/// A private directory per test, removed on scope exit.
class SeqTempDir {
public:
    SeqTempDir() {
        static std::atomic<unsigned> counter{0};
        mPath = std::filesystem::temp_directory_path() /
                ("meshio_seq_" + std::to_string(counter++) + "_" +
                 std::to_string(::testing::UnitTest::GetInstance()->random_seed()));
        std::filesystem::create_directories(mPath);
    }
    ~SeqTempDir() {
        std::error_code ec;
        std::filesystem::remove_all(mPath, ec);
    }
    SeqTempDir(const SeqTempDir&) = delete;
    SeqTempDir& operator=(const SeqTempDir&) = delete;

    std::string operator/(const std::string& rName) const { return (mPath / rName).string(); }
    const std::filesystem::path& Path() const { return mPath; }

private:
    std::filesystem::path mPath;
};

/// Write `Count` single-step `.vtu` files named `<stem><i>.vtu` into `rDir`.
std::vector<std::string> seq_write_files(const SeqTempDir& rDir, const std::string& rStem,
                                         std::size_t Count) {
    std::vector<std::string> out;
    for (std::size_t i = 0; i < Count; ++i) {
        const std::string path = rDir / (rStem + std::to_string(i) + ".vtu");
        Mesh m = mt::tri_mesh();
        meshioplusplus::registry_writers().at("vtu")(path, m);
        out.push_back(path);
    }
    return out;
}

}  // namespace

// ---------------------------------------------------------------------------
// The pure units. No I/O -- these are what the whole feature's correctness
// rests on, and they are the first thing to get green.
// ---------------------------------------------------------------------------

TEST(SequenceOrdering, NaturalNumericBeatsLexicographic) {
    using meshioplusplus::sequence_natural_less;

    // The headline case: lexicographically "out_10" < "out_9".
    EXPECT_TRUE(sequence_natural_less("out_9.vtu", "out_10.vtu"));
    EXPECT_FALSE(sequence_natural_less("out_10.vtu", "out_9.vtu"));

    // Zero padding does not change the numeric value.
    EXPECT_TRUE(sequence_natural_less("out_0009.vtu", "out_0010.vtu"));
    EXPECT_TRUE(sequence_natural_less("out_0009.vtu", "out_10.vtu"));
    EXPECT_TRUE(sequence_natural_less("out_9.vtu", "out_0010.vtu"));

    // Two numeric runs: the first run decides.
    EXPECT_TRUE(sequence_natural_less("run2_step10.vtu", "run10_step2.vtu"));
    EXPECT_TRUE(sequence_natural_less("run2_step2.vtu", "run2_step10.vtu"));

    // The whole path is compared, not the basename.
    EXPECT_TRUE(sequence_natural_less("a/out_9.vtu", "a/out_10.vtu"));
    EXPECT_TRUE(sequence_natural_less("a/out_10.vtu", "b/out_1.vtu"));

    // A digit run sorts before a non-digit run at the same position.
    EXPECT_TRUE(sequence_natural_less("1abc", "abc1"));

    // Equal numeric value, different padding: rule 5's tie-break, deterministic
    // in both directions.
    EXPECT_TRUE(sequence_natural_less("out_01", "out_1"));
    EXPECT_FALSE(sequence_natural_less("out_1", "out_01"));

    // Equal strings are never less than each other (irreflexivity).
    EXPECT_FALSE(sequence_natural_less("out_1.vtu", "out_1.vtu"));

    // A 40-digit run must not overflow -- this is why the comparison is done on
    // the digits and never through stoull.
    const std::string big_a = "x" + std::string(39, '9') + "8";
    const std::string big_b = "x" + std::string(40, '9');
    EXPECT_TRUE(sequence_natural_less(big_a, big_b));

    // Empty and digits-only inputs.
    EXPECT_TRUE(sequence_natural_less("", "a"));
    EXPECT_FALSE(sequence_natural_less("a", ""));
    EXPECT_TRUE(sequence_natural_less("2", "10"));

    // High bytes compare as unsigned char, so a UTF-8 name sorts after ASCII
    // rather than before it (which is what a signed char would give).
    EXPECT_TRUE(sequence_natural_less("a", "\xC3\xA9"));
}

TEST(SequenceOrdering, IsAStrictWeakOrdering) {
    using meshioplusplus::sequence_natural_less;
    // Without the rule-5 tie-break, "out_1" and "out_01" would be mutually
    // "not less" yet not equivalent, and std::sort would be undefined
    // behaviour on any directory mixing padded and unpadded names. Brute-force
    // all four axioms rather than trusting the argument.
    const std::vector<std::string> t = {
        "",         "0",        "00",       "1",      "01",       "001",      "2",
        "9",        "10",       "010",      "a",      "a0",       "a1",       "a01",
        "a10",      "a9",       "b",        "0a",     "1a",       "a_1_b",    "a_01_b",
        "a_2_b",    "a_10_b",   "out_9",    "out_10", "out_0009", "out_0010", "x",
        "\xC3\xA9", "run2_s10", "run10_s2", "z999",   "z1000",
    };
    for (const std::string& a : t) {
        EXPECT_FALSE(sequence_natural_less(a, a)) << "irreflexive: " << a;
        for (const std::string& b : t) {
            if (sequence_natural_less(a, b))
                EXPECT_FALSE(sequence_natural_less(b, a)) << "asymmetric: " << a << " / " << b;
            for (const std::string& c : t) {
                if (sequence_natural_less(a, b) && sequence_natural_less(b, c))
                    EXPECT_TRUE(sequence_natural_less(a, c))
                        << "transitive: " << a << " / " << b << " / " << c;
                // Transitivity of incomparability -- the axiom the tie-break
                // is actually there to satisfy.
                const bool ab = !sequence_natural_less(a, b) && !sequence_natural_less(b, a);
                const bool bc = !sequence_natural_less(b, c) && !sequence_natural_less(c, b);
                if (ab && bc) {
                    const bool ac = !sequence_natural_less(a, c) && !sequence_natural_less(c, a);
                    EXPECT_TRUE(ac) << "incomparability: " << a << " / " << b << " / " << c;
                }
            }
        }
    }
}

TEST(SequenceGlob, MatchesStarAndQuestionOnly) {
    using meshioplusplus::sequence_glob_match;

    EXPECT_TRUE(sequence_glob_match("out_*.vtu", "out_0001.vtu"));
    EXPECT_TRUE(sequence_glob_match("*", "anything"));
    EXPECT_TRUE(sequence_glob_match("*", ""));
    EXPECT_TRUE(sequence_glob_match("*.vtu", ".vtu"));
    EXPECT_FALSE(sequence_glob_match("out_*.vtu", "out_0001.vtk"));
    EXPECT_FALSE(sequence_glob_match("out_*.vtu", "in_0001.vtu"));

    // `?` is exactly one character, never zero.
    EXPECT_TRUE(sequence_glob_match("out_?.vtu", "out_1.vtu"));
    EXPECT_FALSE(sequence_glob_match("out_?.vtu", "out_.vtu"));
    EXPECT_FALSE(sequence_glob_match("out_?.vtu", "out_12.vtu"));

    // Multiple and adjacent stars.
    EXPECT_TRUE(sequence_glob_match("a*b*c", "abc"));
    EXPECT_TRUE(sequence_glob_match("a*b*c", "axxbyyc"));
    EXPECT_TRUE(sequence_glob_match("**", "ab"));
    EXPECT_FALSE(sequence_glob_match("a*b*c", "abd"));

    // The pattern language is deliberately narrower than glob(3): `[abc]` is
    // three literal characters, not a set, so that the C++ matcher and its
    // Python twin cannot accept different things.
    EXPECT_FALSE(sequence_glob_match("f[ab].vtu", "fa.vtu"));
    EXPECT_TRUE(sequence_glob_match("f[ab].vtu", "f[ab].vtu"));

    // `**` is not a recursive wildcard -- it is just two stars, and neither
    // crosses a separator differently from one.
    EXPECT_TRUE(sequence_glob_match("**.vtu", "a/b.vtu"));

    // The classic backtracking blow-up: must terminate, and must be false.
    EXPECT_FALSE(sequence_glob_match("a*a*a*a*a*a*b", std::string(40, 'a')));
}

TEST(SequencePattern, ExpandsStepAndIndexWithPadding) {
    using meshioplusplus::sequence_expand_pattern;
    using meshioplusplus::sequence_pattern_has_token;

    EXPECT_TRUE(sequence_pattern_has_token("out_{step}.vtu"));
    EXPECT_TRUE(sequence_pattern_has_token("out_{index}.vtu"));
    EXPECT_FALSE(sequence_pattern_has_token("out.vtu"));

    // `{step}` pads to at least 4; `{index}` never pads.
    EXPECT_EQ(sequence_expand_pattern("out_{step}.vtu", 0, 12), "out_0000.vtu");
    EXPECT_EQ(sequence_expand_pattern("out_{step}.vtu", 11, 12), "out_0011.vtu");
    EXPECT_EQ(sequence_expand_pattern("out_{index}.vtu", 11, 12), "out_11.vtu");

    // The width grows once the count needs it, so a naive `ls` still sorts.
    EXPECT_EQ(sequence_expand_pattern("out_{step}.vtu", 7, 20000), "out_00007.vtu");
    EXPECT_EQ(sequence_expand_pattern("out_{step}.vtu", 0, 1), "out_0000.vtu");

    // Both tokens, and repeats of one.
    EXPECT_EQ(sequence_expand_pattern("{index}/out_{step}.vtu", 3, 5), "3/out_0003.vtu");
    EXPECT_EQ(sequence_expand_pattern("{step}_{step}.vtu", 3, 5), "0003_0003.vtu");

    // Substring replace-all semantics: an unrelated brace is a literal here,
    // deliberately matching the native CLI's {key}/{part} helpers rather than
    // str.format (whose Python twin raises on this). The asymmetry is
    // documented, not accidental.
    EXPECT_EQ(sequence_expand_pattern("o{ther}_{step}.vtu", 1, 5), "o{ther}_0001.vtu");
}

// ---------------------------------------------------------------------------
// Expansion, ordering in practice, and the time-value precedence.
// ---------------------------------------------------------------------------

TEST(SequenceExpand, OrdersNaturallyNotLexicographically) {
    SeqTempDir dir;
    seq_write_files(dir, "out_", 12);

    SequenceInput in;
    in.mPattern = dir / "out_*.vtu";
    const std::vector<SequenceEntry> entries = meshioplusplus::sequence_expand(in);

    ASSERT_EQ(entries.size(), 12u);
    for (std::size_t i = 0; i < entries.size(); ++i) {
        const std::string want = "out_" + std::to_string(i) + ".vtu";
        EXPECT_EQ(std::filesystem::path(entries[i].mPath).filename().string(), want);
    }
}

TEST(SequenceExpand, EmptyPatternIsAnErrorNotAnEmptySequence) {
    SeqTempDir dir;
    SequenceInput in;
    in.mPattern = dir / "nothing_*.vtu";
    EXPECT_THROW(meshioplusplus::sequence_expand(in), meshioplusplus::ReadError);
}

TEST(SequenceExpand, RejectsBothOrNeitherInput) {
    SequenceInput neither;
    EXPECT_THROW(meshioplusplus::sequence_expand(neither), std::invalid_argument);

    SequenceInput both;
    both.mPaths = {"a.vtu"};
    both.mPattern = "*.vtu";
    EXPECT_THROW(meshioplusplus::sequence_expand(both), std::invalid_argument);
}

TEST(SequenceExpand, ExplicitListKeepsItsOrderUnlessAsked) {
    SeqTempDir dir;
    const std::vector<std::string> files = seq_write_files(dir, "out_", 3);

    SequenceInput in;
    in.mPaths = {files[2], files[0], files[1]};  // a caller's stated order
    std::vector<SequenceEntry> entries = meshioplusplus::sequence_expand(in);
    EXPECT_EQ(entries[0].mPath, files[2]);

    in.mSortExplicit = true;
    entries = meshioplusplus::sequence_expand(in);
    EXPECT_EQ(entries[0].mPath, files[0]);
    EXPECT_EQ(entries[2].mPath, files[2]);
}

TEST(SequenceExpand, PatternDirectoryComponentIsLiteral) {
    SequenceInput in;
    in.mPattern = "some*dir/out_*.vtu";
    EXPECT_THROW(meshioplusplus::sequence_expand(in), std::invalid_argument);
}

TEST(SequenceTime, ExplicitBeatsEverythingElse) {
    SeqTempDir dir;
    seq_write_files(dir, "out_", 3);

    SequenceInput in;
    in.mPattern = dir / "out_*.vtu";
    in.mTimes = {10.0, 20.0, 30.0};
    const std::vector<SequenceEntry> entries = meshioplusplus::sequence_expand(in);

    ASSERT_EQ(entries.size(), 3u);
    EXPECT_DOUBLE_EQ(entries[1].mTime, 20.0);
    EXPECT_EQ(entries[1].mTimeSource, SequenceTimeSource::Explicit);
}

TEST(SequenceTime, ExplicitLengthMismatchNamesBothCounts) {
    SeqTempDir dir;
    seq_write_files(dir, "out_", 3);

    SequenceInput in;
    in.mPattern = dir / "out_*.vtu";
    in.mTimes = {10.0, 20.0};
    try {
        meshioplusplus::sequence_expand(in);
        FAIL() << "expected a length mismatch";
    } catch (const std::invalid_argument& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find('2'), std::string::npos);
        EXPECT_NE(msg.find('3'), std::string::npos);
    }
}

TEST(SequenceTime, FilenameIsTheLastDigitRunOfTheStem) {
    SeqTempDir dir;
    // A directory whose own name carries digits, to pin "last run of the STEM".
    std::filesystem::create_directories(dir.Path() / "run17");
    for (int i : {1, 42}) {
        Mesh m = mt::tri_mesh();
        meshioplusplus::registry_writers().at("vtu")(
            (dir.Path() / "run17" / ("out_" + std::string(i == 1 ? "0001" : "0042") + ".vtu"))
                .string(),
            m);
    }

    SequenceInput in;
    in.mPattern = (dir.Path() / "run17" / "out_*.vtu").string();
    const std::vector<SequenceEntry> entries = meshioplusplus::sequence_expand(in);

    ASSERT_EQ(entries.size(), 2u);
    EXPECT_DOUBLE_EQ(entries[0].mTime, 1.0);
    EXPECT_DOUBLE_EQ(entries[1].mTime, 42.0);
    EXPECT_EQ(entries[0].mTimeSource, SequenceTimeSource::Filename);
}

TEST(SequenceTime, IndexIsTheFallbackAndIsRecordedAsSuch) {
    SeqTempDir dir;
    // Names with no digits at all: nothing to parse, so the index it is.
    for (const char* stem : {"alpha", "beta", "gamma"}) {
        Mesh m = mt::tri_mesh();
        meshioplusplus::registry_writers().at("vtu")(dir / (std::string(stem) + ".vtu"), m);
    }

    SequenceInput in;
    in.mPattern = dir / "*.vtu";
    const std::vector<SequenceEntry> entries = meshioplusplus::sequence_expand(in);

    ASSERT_EQ(entries.size(), 3u);
    for (std::size_t i = 0; i < entries.size(); ++i) {
        EXPECT_DOUBLE_EQ(entries[i].mTime, static_cast<double>(i));
        EXPECT_EQ(entries[i].mTimeSource, SequenceTimeSource::Index);
    }
}

TEST(SequenceTime, TimeFromIndexOverridesAParseableFilename) {
    SeqTempDir dir;
    seq_write_files(dir, "out_", 3);

    SequenceInput in;
    in.mPattern = dir / "out_*.vtu";
    in.mTimeFrom = SequenceTimeFrom::Index;
    const std::vector<SequenceEntry> entries = meshioplusplus::sequence_expand(in);

    ASSERT_EQ(entries.size(), 3u);
    EXPECT_EQ(entries[2].mTimeSource, SequenceTimeSource::Index);
    EXPECT_DOUBLE_EQ(entries[2].mTime, 2.0);
}

// ---------------------------------------------------------------------------
// Registry-derived capability queries.
// ---------------------------------------------------------------------------

TEST(SequenceCapability, NumStepsIsAtLeastOneAndComesFromTheRegistry) {
    SeqTempDir dir;
    const std::vector<std::string> files = seq_write_files(dir, "out_", 1);

    // A .vtu carries no time concept, so one step -- never zero.
    EXPECT_EQ(meshioplusplus::sequence_num_steps(files[0], ""), 1u);
    // An unreadable path is a capability question, not a read: it answers 1 and
    // lets the real read produce the diagnostics.
    EXPECT_EQ(meshioplusplus::sequence_num_steps(dir / "missing.vtu", ""), 1u);
}

TEST(SequenceCapability, WriteSupportsTimeAgreesWithReality) {
    // The anti-drift gate: a format that grows a multi-step writer without
    // updating the predicate turns this red naming itself.
    SeqTempDir dir;
    const std::vector<std::string> files = seq_write_files(dir, "in_", 2);

    for (const auto& [fmt, writer] : meshioplusplus::registry_writers()) {
        (void)writer;
        std::string why;
        const bool claimed = meshioplusplus::sequence_write_supports_time(fmt, why);
        EXPECT_EQ(claimed, why.empty()) << fmt << ": rWhy must be set iff unsupported";

        SequenceInput in;
        in.mPaths = files;
        SequenceOutput out;
        out.mPath = dir / ("series_" + fmt + ".out");
        out.mFormat = fmt;
        bool actually_works = true;
        try {
            meshioplusplus::sequence_to_timeseries(in, out);
        } catch (const std::exception&) {
            actually_works = false;
        }
        EXPECT_EQ(claimed, actually_works)
            << "format '" << fmt << "': sequence_write_supports_time says " << claimed
            << " but a real two-step fan-in " << (actually_works ? "succeeded" : "failed");
    }
}

TEST(SequenceCapability, NonTimeCarryingTargetFailsByName) {
    SeqTempDir dir;
    const std::vector<std::string> files = seq_write_files(dir, "in_", 2);

    SequenceInput in;
    in.mPaths = files;
    SequenceOutput out;
    out.mPath = dir / "series.vtu";
    try {
        meshioplusplus::sequence_to_timeseries(in, out);
        FAIL() << "expected a fan-in to .vtu to fail";
    } catch (const meshioplusplus::WriteError& e) {
        const std::string msg = e.what();
        // Names the format, and names the remedy -- never a silent truncation
        // to step 0.
        EXPECT_NE(msg.find("'vtu'"), std::string::npos) << msg;
        EXPECT_NE(msg.find("{step}"), std::string::npos) << msg;
    }
}

TEST(SequenceCapability, TheTransientWriterRejectsOptionsItCannotHonour) {
    // sequence_to_timeseries drives XdmfTimeSeriesWriter directly rather than
    // going through registry_write_ex, so it must enforce write_options.hpp's
    // rule itself: an option the writer cannot honour is an error, never
    // silently ignored. Only Encoding has anywhere to go.
    SeqTempDir dir;
    const std::vector<std::string> files = seq_write_files(dir, "in_", 2);
    SequenceInput in;
    in.mPaths = files;

    SequenceOutput bad_codec;
    bad_codec.mPath = dir / "s1.xdmf";
    bad_codec.mOptions.mCodec = meshioplusplus::detail::VtkCodec::Zlib;
    bad_codec.mOptions.mCodecSet = true;
    EXPECT_THROW(meshioplusplus::sequence_to_timeseries(in, bad_codec), meshioplusplus::WriteError);

    SequenceOutput bad_float;
    bad_float.mPath = dir / "s2.xdmf";
    bad_float.mOptions.mFloatFormat = ".9e";
    EXPECT_THROW(meshioplusplus::sequence_to_timeseries(in, bad_float), meshioplusplus::WriteError);

    // Encoding alone is fine -- it is what selects XML over HDF.
    SequenceOutput ok;
    ok.mPath = dir / "s3.xdmf";
    ok.mOptions.mEncoding = meshioplusplus::WriteEncoding::Ascii;
    EXPECT_NO_THROW(meshioplusplus::sequence_to_timeseries(in, ok));
}

// ---------------------------------------------------------------------------
// Mode inference.
// ---------------------------------------------------------------------------

TEST(SequenceMode, InferenceFollowsTheDocumentedTable) {
    using meshioplusplus::sequence_resolve_mode;

    auto entry = [](const char* path, std::size_t step) {
        SequenceEntry e;
        e.mPath = path;
        e.mStep = step;
        return e;
    };
    SequenceOutput token;
    token.mPath = "out_{step}.vtu";
    SequenceOutput plain;
    plain.mPath = "out.xdmf";

    const std::vector<SequenceEntry> one = {entry("a.vtu", 0)};
    const std::vector<SequenceEntry> many_files = {entry("a.vtu", 0), entry("b.vtu", 0)};
    const std::vector<SequenceEntry> many_steps = {entry("a.xdmf", 0), entry("a.xdmf", 1)};

    EXPECT_EQ(sequence_resolve_mode(one, plain, SequenceMode::Auto), SequenceMode::Sequence);
    EXPECT_EQ(sequence_resolve_mode(one, token, SequenceMode::Auto), SequenceMode::Sequence);
    EXPECT_EQ(sequence_resolve_mode(many_files, plain, SequenceMode::Auto), SequenceMode::FanIn);
    EXPECT_EQ(sequence_resolve_mode(many_files, token, SequenceMode::Auto), SequenceMode::Sequence);
    EXPECT_EQ(sequence_resolve_mode(many_steps, token, SequenceMode::Auto), SequenceMode::FanOut);

    // The refusal-to-truncate case: several steps aimed at one single-step file.
    try {
        sequence_resolve_mode(many_steps, plain, SequenceMode::Auto);
        FAIL() << "expected a refusal to truncate";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("{step}"), std::string::npos);
    }
}

TEST(SequenceMode, AStatedModeAssertsRatherThanSelects) {
    using meshioplusplus::sequence_resolve_mode;
    SequenceEntry a;
    a.mPath = "a.vtu";
    SequenceOutput plain;
    plain.mPath = "out.xdmf";

    // One file inferred as Sequence; asserting FanIn must fail rather than
    // quietly do something else -- that is the whole point of allowing Mode.
    try {
        sequence_resolve_mode({a}, plain, SequenceMode::FanIn);
        FAIL() << "expected a mode mismatch";
    } catch (const std::invalid_argument& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("fan-in"), std::string::npos) << msg;
        EXPECT_NE(msg.find("sequence"), std::string::npos) << msg;
    }
}

TEST(SequenceMode, NamesParseAndRejectByName) {
    EXPECT_EQ(meshioplusplus::sequence_mode_from_name(""), SequenceMode::Auto);
    EXPECT_EQ(meshioplusplus::sequence_mode_from_name("fan-in"), SequenceMode::FanIn);
    EXPECT_STREQ(meshioplusplus::sequence_mode_name(SequenceMode::FanOut), "fan-out");
    EXPECT_THROW(meshioplusplus::sequence_mode_from_name("nope"), std::invalid_argument);

    EXPECT_EQ(meshioplusplus::sequence_time_from_name("filename"), SequenceTimeFrom::Filename);
    EXPECT_THROW(meshioplusplus::sequence_time_from_name("nope"), std::invalid_argument);

    EXPECT_STREQ(meshioplusplus::sequence_time_source_name(SequenceTimeSource::File), "file");
}

// ---------------------------------------------------------------------------
// The drivers, end to end.
// ---------------------------------------------------------------------------

TEST(SequenceDriver, FanOutRequiresAPattern) {
    SeqTempDir dir;
    const std::vector<std::string> files = seq_write_files(dir, "in_", 1);
    SequenceOutput out;
    out.mPath = dir / "plain.vtu";
    EXPECT_THROW(meshioplusplus::timeseries_to_sequence(files[0], "", ReadOptions{}, out),
                 std::invalid_argument);
}

TEST(SequenceDriver, FanOutWritesOneFilePerStep) {
    SeqTempDir dir;
    const std::vector<std::string> files = seq_write_files(dir, "in_", 1);

    SequenceOutput out;
    out.mPath = dir / "step_{step}.vtu";
    meshioplusplus::timeseries_to_sequence(files[0], "", ReadOptions{}, out);

    const std::string written = dir / "step_0000.vtu";
    ASSERT_TRUE(std::filesystem::exists(written));
    mt::expect_same_geometry(mt::tri_mesh(), meshioplusplus::registry_readers().at("vtu")(written));
}

TEST(SequenceDriver, FanOutTimesSurviveOnlyWhereFieldDataDoes) {
    // A documented limitation, pinned so it cannot regress silently in either
    // direction: `timeseries_to_sequence` always attaches
    // `field_data["meshio:time"]`, but only a handful of formats carry
    // field_data at all -- VTU carries none, in either direction. So a fan-out
    // to .vtu followed by a fan-in recovers the step INDEX (from the `{step}`
    // filename), not the original time value, and `mTimeSource` says so rather
    // than letting a caller believe otherwise.
    SeqTempDir dir;
    const std::vector<std::string> files = seq_write_files(dir, "in_", 1);

    SequenceOutput out;
    out.mPath = dir / "step_{step}.vtu";
    meshioplusplus::timeseries_to_sequence(files[0], "", ReadOptions{}, out);

    const Mesh back = meshioplusplus::registry_readers().at("vtu")(dir / "step_0000.vtu");
    EXPECT_FALSE(back.HasFieldData(meshioplusplus::kSequenceTimeKey))
        << "VTU grew field_data support; doc/sequences.md's round-trip note needs revisiting";

    SequenceInput in;
    in.mPattern = dir / "step_*.vtu";
    const std::vector<SequenceEntry> entries = meshioplusplus::sequence_expand(in);
    ASSERT_EQ(entries.size(), 1u);
    EXPECT_EQ(entries[0].mTimeSource, SequenceTimeSource::Filename);
}

TEST(SequenceDriver, ExplicitTimesAreTheRoundTripRecipe) {
    // The documented way to close a fan-out -> fan-in round trip on time for a
    // format that cannot carry it: hand the times back in.
    SeqTempDir dir;
    seq_write_files(dir, "in_", 3);

    SequenceInput in;
    in.mPattern = dir / "in_*.vtu";
    in.mTimes = {0.0, 0.125, 0.25};
    SequenceOutput out;
    out.mPath = dir / "series.xdmf";
    meshioplusplus::sequence_to_timeseries(in, out);

    ASSERT_TRUE(std::filesystem::exists(out.mPath));
    // The series knows its own times, so reading it back needs no help.
    SequenceInput back;
    back.mPaths = {out.mPath};
    const std::vector<SequenceEntry> entries = meshioplusplus::sequence_expand(back);
    ASSERT_EQ(entries.size(), 3u);
    EXPECT_EQ(entries[1].mTimeSource, SequenceTimeSource::File);
    EXPECT_DOUBLE_EQ(entries[1].mTime, 0.125);
    EXPECT_EQ(entries[2].mStep, 2u);
}

TEST(SequenceDriver, FanInThenFanOutReproducesEveryStep) {
    SeqTempDir dir;
    seq_write_files(dir, "in_", 4);

    SequenceInput in;
    in.mPattern = dir / "in_*.vtu";
    SequenceOutput series;
    series.mPath = dir / "series.xdmf";
    meshioplusplus::sequence_to_timeseries(in, series);

    SequenceOutput fanned;
    fanned.mPath = dir / "back_{step}.vtu";
    meshioplusplus::timeseries_to_sequence(series.mPath, "", ReadOptions{}, fanned);

    for (std::size_t i = 0; i < 4; ++i) {
        const std::string path =
            meshioplusplus::sequence_expand_pattern(dir / "back_{step}.vtu", i, 4);
        ASSERT_TRUE(std::filesystem::exists(path)) << path;
        mt::expect_same_geometry(mt::tri_mesh(),
                                 meshioplusplus::registry_readers().at("vtu")(path));
    }
}

TEST(SequenceDriver, SequenceModeAppliesTheChainPerStep) {
    SeqTempDir dir;
    seq_write_files(dir, "in_", 3);

    SequencePipeline p;
    p.mInput.mPattern = dir / "in_*.vtu";
    p.mSteps.push_back({"Quality", {}});
    p.mOutput.mPath = dir / "out_{step}.vtu";

    const meshioplusplus::PipelineReport report = meshioplusplus::run_sequence_pipeline(p);

    // One report entry per (step, op) -- three steps, one op each.
    EXPECT_EQ(report.mSteps.size(), 3u);
    for (std::size_t i = 0; i < 3; ++i) {
        const std::string path =
            meshioplusplus::sequence_expand_pattern(dir / "out_{step}.vtu", i, 3);
        ASSERT_TRUE(std::filesystem::exists(path)) << path;
        const Mesh m = meshioplusplus::registry_readers().at("vtu")(path);
        EXPECT_TRUE(m.HasCellData("quality:scaled_jacobian"));
    }
}

TEST(SequenceDriver, ValidatesTheChainBeforeReadingAnything) {
    SeqTempDir dir;
    seq_write_files(dir, "in_", 2);

    SequencePipeline p;
    p.mInput.mPattern = dir / "in_*.vtu";
    p.mSteps.push_back({"NoSuchOp", {}});
    p.mOutput.mPath = dir / "out_{step}.vtu";

    EXPECT_THROW(meshioplusplus::run_sequence_pipeline(p), std::invalid_argument);
    EXPECT_FALSE(std::filesystem::exists(dir / "out_0000.vtu"));
}

TEST(SequenceDriver, ReadStepRejectsAnIndexPastTheEnd) {
    SeqTempDir dir;
    const std::vector<std::string> files = seq_write_files(dir, "in_", 2);
    SequenceInput in;
    in.mPaths = files;
    const std::vector<SequenceEntry> entries = meshioplusplus::sequence_expand(in);
    EXPECT_THROW(meshioplusplus::sequence_read_step(entries, 5, "", ReadOptions{}),
                 std::out_of_range);
}

// ---------------------------------------------------------------------------
// The streaming invariant (see @ref sequence_streaming). This is a contract,
// not an optimization: the whole feature exists so a 500-step dataset is
// traversable on a laptop, so it is tested rather than asserted in prose.
// ---------------------------------------------------------------------------

namespace {

/// Live-bytes high-water mark over the `BufferAllocator` hook, which every
/// owning `NDArray` allocation routes through. The `test_ndarray_allocator.cpp`
/// machinery plus a peak.
struct SeqCounters {
    std::size_t live = 0;
    std::size_t peak = 0;
};

void* seq_counting_alloc(std::size_t bytes, void* pUser) {
    SeqCounters* c = static_cast<SeqCounters*>(pUser);
    c->live += bytes;
    c->peak = std::max(c->peak, c->live);
    return ::operator new(bytes);
}

void seq_counting_free(void* pPtr, std::size_t bytes, void* pUser) {
    SeqCounters* c = static_cast<SeqCounters*>(pUser);
    c->live -= std::min(c->live, bytes);
    ::operator delete(pPtr);
}

class SeqScopedAllocator {
public:
    explicit SeqScopedAllocator(SeqCounters& rCounters) {
        auto a = std::make_shared<meshioplusplus::BufferAllocator>();
        a->alloc = &seq_counting_alloc;
        a->free = &seq_counting_free;
        a->pUser = &rCounters;
        meshioplusplus::set_buffer_allocator(std::move(a));
    }
    ~SeqScopedAllocator() { meshioplusplus::set_buffer_allocator(nullptr); }
};

std::size_t seq_fan_in_peak(std::size_t Count) {
    SeqTempDir dir;
    seq_write_files(dir, "in_", Count);
    SequenceInput in;
    in.mPattern = dir / "in_*.vtu";
    SequenceOutput out;
    out.mPath = dir / "series.xdmf";

    SeqCounters c;
    {
        SeqScopedAllocator guard(c);
        meshioplusplus::sequence_to_timeseries(in, out);
    }
    return c.peak;
}

}  // namespace

TEST(SequenceStreaming, FanInPeakIsConstantInTheStepCount) {
    // The real assertion is not "the peak is small" but "the peak does not grow
    // with the sequence length" -- nothing that accumulates meshes can pass it.
    const std::size_t peak20 = seq_fan_in_peak(20);
    const std::size_t peak40 = seq_fan_in_peak(40);
    ASSERT_GT(peak20, 0u) << "the allocator hook saw nothing; the test is inert";
    EXPECT_LT(peak40, peak20 * 3 / 2)
        << "fan-in peak grew from " << peak20 << " to " << peak40
        << " bytes when the step count doubled: something is buffering the sequence";
}

TEST(SequenceStreaming, ReadStepHoldsExactlyOneMesh) {
    // The lazy read surface: reading step i must not retain step i-1. Measured
    // through the same hook, over a loop that drops each mesh immediately.
    SeqTempDir dir;
    seq_write_files(dir, "in_", 12);
    SequenceInput in;
    in.mPattern = dir / "in_*.vtu";
    const std::vector<SequenceEntry> entries = meshioplusplus::sequence_expand(in);

    SeqCounters c;
    {
        SeqScopedAllocator guard(c);
        for (std::size_t i = 0; i < entries.size(); ++i) {
            Mesh m = meshioplusplus::sequence_read_step(entries, i, "", ReadOptions{});
            EXPECT_GT(m.NumPoints(), 0u);
        }
    }
    EXPECT_EQ(c.live, 0u) << "every step's buffers must be released before the next";
}

// ---------------------------------------------------------------------------
// The JSON front-end. Guarded exactly like the pipeline's own parser: the
// typed driver above compiles either way, but a settings document cannot be
// read without a parser, and the compiled-out entry points must fail naming the
// option rather than link-erroring or silently doing nothing.
// ---------------------------------------------------------------------------

#ifdef MESHIOPLUSPLUS_HAS_JSON

TEST(SequenceJson, ParsesAFullSequenceDocument) {
    const std::string text = R"({
        "Version": 1,
        "Mode": "sequence",
        "Parallel": true,
        "Workers": 4,
        "Input": {"Pattern": "in/out_*.vtu", "Times": [0.0, 0.5],
                  "TimeFrom": "filename", "Format": "vtu",
                  "Options": {"PointsOnly": true}},
        "Operations": [{"Op": "Clean", "RemoveOrphans": false}],
        "Output": {"Path": "out_{step}.vtu", "Encoding": "binary"}
    })";
    const SequencePipeline p = meshioplusplus::parse_sequence_json(text);
    EXPECT_EQ(p.mMode, SequenceMode::Sequence);
    EXPECT_TRUE(p.mParallel);
    EXPECT_EQ(p.mWorkers, 4);
    EXPECT_EQ(p.mInput.mPattern, "in/out_*.vtu");
    EXPECT_TRUE(p.mInput.mPaths.empty());
    ASSERT_EQ(p.mInput.mTimes.size(), 2u);
    EXPECT_DOUBLE_EQ(p.mInput.mTimes[1], 0.5);
    EXPECT_EQ(p.mInput.mTimeFrom, SequenceTimeFrom::Filename);
    EXPECT_TRUE(p.mInput.mOptions.mPointsOnly);
    ASSERT_EQ(p.mSteps.size(), 1u);
    EXPECT_EQ(p.mSteps[0].mOp, "Clean");
    EXPECT_EQ(p.mOutput.mPath, "out_{step}.vtu");
}

TEST(SequenceJson, ExplicitPathListParses) {
    const SequencePipeline p = meshioplusplus::parse_sequence_json(
        R"({"Input": {"Paths": ["a.vtu", "b.vtu"]}, "Output": {"Path": "o.xdmf"}})");
    ASSERT_EQ(p.mInput.mPaths.size(), 2u);
    EXPECT_EQ(p.mInput.mPaths[1], "b.vtu");
    // An explicit list is a stated order and is not re-sorted.
    EXPECT_FALSE(p.mInput.mSortExplicit);
}

TEST(SequenceJson, APlainDocumentStillParsesAsBoth) {
    // The shared parser must leave the v9.11.0 shape exactly where it was.
    const std::string text =
        R"({"Version": 1, "Input": {"Path": "a.vtu"}, "Output": {"Path": "b.vtu"}})";
    const SequencePipeline seq = meshioplusplus::parse_sequence_json(text);
    ASSERT_EQ(seq.mInput.mPaths.size(), 1u);
    EXPECT_EQ(seq.mInput.mPaths[0], "a.vtu");

    const meshioplusplus::Pipeline single = meshioplusplus::parse_pipeline_json(text);
    EXPECT_EQ(single.mInput.mPath, "a.vtu");
    EXPECT_EQ(single.mOutput.mPath, "b.vtu");
}

TEST(SequenceJson, SequenceKeysAreRefusedByTheSingleFileParser) {
    // Ignoring them instead would quietly run a transient document as its first
    // step, which is the failure this whole feature exists to prevent.
    const char* cases[] = {
        R"({"Mode": "fan-in", "Input": {"Path": "a"}, "Output": {"Path": "b"}})",
        R"({"Input": {"Pattern": "*.vtu"}, "Output": {"Path": "b"}})",
        R"({"Input": {"Paths": ["a"]}, "Output": {"Path": "b"}})",
        R"({"Input": {"Path": "a", "Times": [0]}, "Output": {"Path": "b"}})",
        R"({"Input": {"Path": "a", "TimeFrom": "index"}, "Output": {"Path": "b"}})",
        R"({"Parallel": true, "Input": {"Path": "a"}, "Output": {"Path": "b"}})",
        R"({"Workers": 2, "Input": {"Path": "a"}, "Output": {"Path": "b"}})",
    };
    for (const char* text : cases) {
        try {
            meshioplusplus::parse_pipeline_json(text);
            FAIL() << "expected a refusal for: " << text;
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find("sequence keys"), std::string::npos) << e.what();
        }
    }
}

TEST(SequenceJson, StrictSchemaErrorsNameTheOffender) {
    const char* cases[][2] = {
        {R"({"Input": {"Path": "a", "Pattern": "*"}, "Output": {"Path": "b"}})",
         "more than one source"},
        {R"({"Input": {"Bogus": 1}, "Output": {"Path": "b"}})", "Bogus"},
        {R"({"Input": {"Paths": []}, "Output": {"Path": "b"}})", "non-empty"},
        {R"({"Input": {"Path": "a", "Times": ["x"]}, "Output": {"Path": "b"}})", "Times"},
        {R"({"Mode": "sideways", "Input": {"Path": "a"}, "Output": {"Path": "b"}})", "Mode"},
        {R"({"Input": {"Path": "a", "TimeFrom": "vibes"}, "Output": {"Path": "b"}})", "TimeFrom"},
        {R"({"Workers": -1, "Input": {"Path": "a"}, "Output": {"Path": "b"}})", "Workers"},
        {R"({"Input": {}, "Output": {"Path": "b"}})", "Input.Path is required"},
    };
    for (const auto& c : cases) {
        try {
            meshioplusplus::parse_sequence_json(c[0]);
            FAIL() << "expected a schema error for: " << c[0];
        } catch (const std::invalid_argument& e) {
            EXPECT_NE(std::string(e.what()).find(c[1]), std::string::npos)
                << "message was: " << e.what();
        }
    }
}

TEST(SequenceJson, RunSequenceFileEndToEnd) {
    SeqTempDir dir;
    seq_write_files(dir, "in_", 3);
    const std::string settings = dir / "settings.json";
    {
        std::ofstream f(settings);
        f << R"({"Version": 1, "Input": {"Pattern": ")" << (dir / "in_*.vtu")
          << R"("}, "Operations": [{"Op": "Quality"}], "Output": {"Path": ")"
          << (dir / "out_{step}.vtu") << R"("}})";
    }
    const meshioplusplus::PipelineReport report = meshioplusplus::run_sequence_file(settings);
    EXPECT_EQ(report.mSteps.size(), 3u);
    EXPECT_TRUE(std::filesystem::exists(dir / "out_0002.vtu"));
}

TEST(SequenceJson, APlainDocumentNamingAMultiStepInputStillRefusesToTruncate) {
    // The routing gap this pins: a document that uses NO sequence key and names
    // a plain output would take the single-file path, which reads step 0 --
    // exactly the silent truncation this layer exists to prevent. Found by the
    // WASM smoke test, whose settings surface shares this decision.
    SeqTempDir dir;
    seq_write_files(dir, "in_", 3);
    SequenceInput in;
    in.mPattern = dir / "in_*.vtu";
    SequenceOutput series;
    series.mPath = dir / "series.xdmf";
    meshioplusplus::sequence_to_timeseries(in, series);

    const std::string text = std::string(R"({"Version": 1, "Input": {"Path": ")") + series.mPath +
                             R"("}, "Output": {"Path": ")" + (dir / "one.vtu") + R"("}})";
    try {
        meshioplusplus::run_sequence_json(text);
        FAIL() << "expected a refusal to truncate a 3-step input";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("{step}"), std::string::npos) << e.what();
    }

    // An explicit TimeStep IS a deliberate single-step selection and opts out.
    const std::string picked = std::string(R"({"Version": 1, "Input": {"Path": ")") + series.mPath +
                               R"(", "Options": {"TimeStep": 1}}, "Output": {"Path": ")" +
                               (dir / "one.vtu") + R"("}})";
    EXPECT_NO_THROW(meshioplusplus::run_sequence_json(picked));
    EXPECT_TRUE(std::filesystem::exists(dir / "one.vtu"));
}

TEST(SequenceJson, APlainDocumentTakesTheUnchangedSingleFilePath) {
    // run_sequence_json must delegate to run_pipeline for a document that uses
    // no sequence key and names a plain output, so an existing settings.json
    // produces byte-identical output right down to which function wrote it.
    SeqTempDir dir;
    const std::vector<std::string> files = seq_write_files(dir, "in_", 1);
    const std::string text = std::string(R"({"Version": 1, "Input": {"Path": ")") + files[0] +
                             R"("}, "Operations": [], "Output": {"Path": ")" +
                             (dir / "single.vtu") + R"("}})";

    meshioplusplus::run_sequence_json(text);
    const std::string via_sequence = dir / "single.vtu";
    ASSERT_TRUE(std::filesystem::exists(via_sequence));
    // No `meshio:time` may be attached: that key labels one step OF A SERIES,
    // and a single named output is not one.
    const Mesh back = meshioplusplus::registry_readers().at("vtu")(via_sequence);
    EXPECT_FALSE(back.HasFieldData(meshioplusplus::kSequenceTimeKey));
    mt::expect_same_geometry(mt::tri_mesh(), back);
}

#else  // MESHIOPLUSPLUS_HAS_JSON

TEST(SequenceJson, EntryPointsThrowNamingTheFlag) {
    // The typed sequence driver stays fully functional in this build; only the
    // settings-document front-end is compiled out.
    for (auto fn : {+[] { meshioplusplus::parse_sequence_json("{}"); },
                    +[] { meshioplusplus::parse_sequence_file("x.json"); },
                    +[] { meshioplusplus::run_sequence_json("{}"); },
                    +[] { meshioplusplus::run_sequence_file("x.json"); }}) {
        try {
            fn();
            FAIL() << "expected runtime_error";
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("MESHIOPLUSPLUS_WITH_JSON"), std::string::npos);
        }
    }
}

TEST(SequenceJson, TheTypedDriverStillWorksWithoutAParser) {
    SeqTempDir dir;
    seq_write_files(dir, "in_", 2);
    SequencePipeline p;
    p.mInput.mPattern = dir / "in_*.vtu";
    p.mOutput.mPath = dir / "out_{step}.vtu";
    EXPECT_NO_THROW(meshioplusplus::run_sequence_pipeline(p));
    EXPECT_TRUE(std::filesystem::exists(dir / "out_0001.vtu"));
}

#endif  // MESHIOPLUSPLUS_HAS_JSON
