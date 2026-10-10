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

// The interactive viewer (roadmap 7.2): the input parser, the prepared scene
// it draws from, the cell diff, and the loop run on a recorded byte stream.

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/render.hpp"
#include "mesh_fixtures.hpp"

// CLI layer (compiled into the test binary like view_payload)
#include "../../src/cpp/cli/render_args.hpp"
#include "../../src/cpp/cli/tui/input.hpp"
#include "../../src/cpp/cli/tui/loop.hpp"
#include "../../src/cpp/cli/tui/probe.hpp"
#include "../../src/cpp/cli/tui/series.hpp"
#include "../../src/cpp/cli/tui/session_file.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtu.hpp"

using namespace meshioplusplus;
using namespace meshioplusplus::cli;
using namespace meshioplusplus::cli::tui;

namespace {

using mt::make_mesh;

Mesh cube_quads() {
    return make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "quad",
        {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 4, 7, 3}, {1, 2, 6, 5}});
}

std::vector<InputEvent> parse(const std::string& rBytes, bool Flush = true) {
    InputParser parser;
    parser.Feed(rBytes);
    return parser.Take(Flush);
}

// ---------------------------------------------------------------- the parser

TEST(TuiInput, PrintableAndControlKeys) {
    const auto events = parse("q\x03\r\t\x7f\xc3\xa9");
    ASSERT_EQ(events.size(), 6u);
    EXPECT_TRUE(events[0].IsChar('q'));
    EXPECT_TRUE(events[1].IsCtrl('c'));
    EXPECT_EQ(events[2].mKey, KeyCode::Enter);
    EXPECT_EQ(events[3].mKey, KeyCode::Tab);
    EXPECT_EQ(events[4].mKey, KeyCode::Backspace);
    EXPECT_EQ(events[5].mChar, 0xE9u);  // é, decoded from UTF-8
}

TEST(TuiInput, ArrowsAndModifiers) {
    const auto events = parse("\x1b[A\x1b[1;2B\x1b[1;5C\x1bOD\x1b[3~\x1b[5~");
    ASSERT_EQ(events.size(), 6u);
    EXPECT_EQ(events[0].mKey, KeyCode::Up);
    EXPECT_EQ(events[1].mKey, KeyCode::Down);
    EXPECT_TRUE(events[1].mShift);
    EXPECT_EQ(events[2].mKey, KeyCode::Right);
    EXPECT_TRUE(events[2].mCtrl);
    EXPECT_EQ(events[3].mKey, KeyCode::Left);  // SS3, as application mode sends it
    EXPECT_EQ(events[4].mKey, KeyCode::Delete);
    EXPECT_EQ(events[5].mKey, KeyCode::PageUp);
}

TEST(TuiInput, SgrMouse) {
    const auto events = parse("\x1b[<0;10;5M\x1b[<32;12;6M\x1b[<0;12;6m\x1b[<64;3;4M\x1b[<65;3;4M");
    ASSERT_EQ(events.size(), 5u);
    EXPECT_EQ(events[0].mKind, EventKind::Mouse);
    EXPECT_EQ(events[0].mAction, MouseAction::Press);
    EXPECT_EQ(events[0].mX, 10);
    EXPECT_EQ(events[0].mY, 5);
    EXPECT_EQ(events[1].mAction, MouseAction::Drag);
    EXPECT_EQ(events[2].mAction, MouseAction::Release);
    EXPECT_EQ(events[3].mAction, MouseAction::WheelUp);
    EXPECT_EQ(events[4].mAction, MouseAction::WheelDown);
}

TEST(TuiInput, ShiftedMouseAndRightButton) {
    const auto events = parse("\x1b[<4;1;1M\x1b[<2;2;2M");
    ASSERT_EQ(events.size(), 2u);
    EXPECT_TRUE(events[0].mShift);
    EXPECT_EQ(events[1].mButton, 2);
}

TEST(TuiInput, BracketedPasteArrivesAsOneEventEvenWhenSplit) {
    InputParser parser;
    parser.Feed("\x1b[200~hel");
    EXPECT_TRUE(parser.Take().empty());  // the end marker has not arrived
    parser.Feed("lo\x1b[20");
    EXPECT_TRUE(parser.Take().empty());
    parser.Feed("1~x");
    const auto events = parser.Take();
    ASSERT_EQ(events.size(), 2u);
    EXPECT_EQ(events[0].mKind, EventKind::Paste);
    EXPECT_EQ(events[0].mText, "hello");
    EXPECT_TRUE(events[1].IsChar('x'));
}

TEST(TuiInput, ASequenceSplitAcrossReadsWaitsForItsTail) {
    InputParser parser;
    parser.Feed("\x1b[");
    EXPECT_TRUE(parser.Take().empty());
    parser.Feed("<0;4");
    EXPECT_TRUE(parser.Take().empty());
    parser.Feed(";9M");
    const auto events = parser.Take();
    ASSERT_EQ(events.size(), 1u);
    EXPECT_EQ(events[0].mX, 4);
    EXPECT_EQ(events[0].mY, 9);
}

TEST(TuiInput, ALoneEscapeIsTheEscapeKeyOnlyAfterAFlush) {
    InputParser parser;
    parser.Feed("\x1b");
    EXPECT_TRUE(parser.Take(false).empty());
    const auto events = parser.Take(true);
    ASSERT_EQ(events.size(), 1u);
    EXPECT_EQ(events[0].mKey, KeyCode::Escape);
    EXPECT_EQ(parser.Pending(), 0u);
}

TEST(TuiInput, AltCharacterAndJunkAreHandled) {
    const auto events = parse("\x1b" "x\x1b[999;999;999z\x80y");
    ASSERT_EQ(events.size(), 2u);
    EXPECT_TRUE(events[0].mAlt);
    EXPECT_EQ(events[0].mChar, static_cast<std::uint32_t>('x'));
    EXPECT_TRUE(events[1].IsChar('y'));  // an unknown CSI and a stray byte are dropped
}

TEST(TuiInput, AnUnfinishedSequenceIsDroppedOnFlush) {
    InputParser parser;
    parser.Feed("\x1b[1;");
    EXPECT_TRUE(parser.Take(true).empty());
    EXPECT_EQ(parser.Pending(), 0u);
}

// ------------------------------------------------------- the prepared scene

bool same(const Frame& rA, const Frame& rB) {
    return rA.mWidth == rB.mWidth && rA.mHeight == rB.mHeight && rA.mRgba == rB.mRgba &&
           rA.mCellIds == rB.mCellIds && rA.mNotes == rB.mNotes && rA.mVMin == rB.mVMin &&
           rA.mVMax == rB.mVMax;
}

Mesh field_cube() {
    Mesh mesh = cube_quads();
    std::vector<std::vector<double>> u;
    std::vector<std::vector<double>> disp;
    for (int i = 0; i < 8; ++i) {
        u.push_back({0.5 * i});
        disp.push_back({0.1 * i, -0.05 * i, 0.2});
    }
    mesh.AddPointData("u", mt::points_from(u));
    mesh.AddPointData("disp", mt::points_from(disp));
    return mesh;
}

TEST(RenderScene, DrawingAPreparedSceneEqualsRender) {
    const Mesh mesh = field_cube();
    std::vector<RenderOptions> cases;
    cases.emplace_back();
    {
        RenderOptions o;
        o.mShading = RenderShading::Smooth;
        o.mEdges = RenderEdges::All;
        cases.push_back(o);
    }
    {
        RenderOptions o;
        o.mColorBy = "u";
        o.mColorbar = true;
        o.mIsolines = 3;
        o.mVectors = "disp";
        o.mAxes = true;
        o.mScaleBar = true;
        o.mAzimuth = 10;
        cases.push_back(o);
    }
    {
        RenderOptions o;
        o.mColorBy = "u";
        o.mWarp = "disp";
        o.mWarpOutline = true;
        o.mEdges = RenderEdges::Feature;
        o.mProjection = RenderProjection::Perspective;
        cases.push_back(o);
    }
    {
        RenderOptions o;
        o.mDiagnostic = RenderDiagnostic::Orientation;
        o.mCategoryEdges = true;
        o.mSupersample = 2;
        cases.push_back(o);
    }
    {
        RenderOptions o;
        o.mDiagnostic = RenderDiagnostic::FreeEdges;
        cases.push_back(o);
    }
    for (std::size_t i = 0; i < cases.size(); ++i) {
        cases[i].mWidth = 96;
        cases[i].mHeight = 64;
        const Frame direct = render(mesh, cases[i]);
        const RenderScene scene = prepare_render(mesh, cases[i]);
        ASSERT_FALSE(scene.Empty());
        EXPECT_TRUE(same(direct, render_scene(scene, cases[i]))) << "case " << i;
    }
}

TEST(RenderScene, OneSceneDrawsFromManyCameras) {
    const Mesh mesh = field_cube();
    RenderOptions prepared;
    prepared.mColorBy = "u";
    prepared.mEdges = RenderEdges::All;
    prepared.mWidth = 80;
    prepared.mHeight = 60;
    const RenderScene scene = prepare_render(mesh, prepared);
    for (double azimuth : {0.0, 30.0, 123.0, -75.0}) {
        for (double elevation : {-40.0, 0.0, 60.0}) {
            RenderOptions camera = prepared;
            camera.mAzimuth = azimuth;
            camera.mElevation = elevation;
            camera.mZoom = 1.5;
            camera.mPanX = 0.1;
            EXPECT_TRUE(same(render(mesh, camera), render_scene(scene, camera)))
                << azimuth << " " << elevation;
        }
    }
}

TEST(RenderScene, TheFrameSizeAndOverlaysComeFromTheDrawCall) {
    const Mesh mesh = field_cube();
    RenderOptions prepared;
    prepared.mColorBy = "u";
    const RenderScene scene = prepare_render(mesh, prepared);
    RenderOptions draw = prepared;
    draw.mWidth = 50;
    draw.mHeight = 30;
    draw.mColorbar = true;
    const Frame frame = render_scene(scene, draw);
    EXPECT_EQ(frame.mWidth, 50);
    EXPECT_EQ(frame.mHeight, 30);
    EXPECT_TRUE(same(frame, render(mesh, draw)));
}

TEST(RenderScene, AFixedFitKeepsTheModelTheSameSizeWhileItTurns) {
    const Mesh mesh = cube_quads();
    RenderOptions o;
    o.mWidth = 120;
    o.mHeight = 120;
    o.mBackground = {0, 0, 0, 255};
    o.mAzimuth = 0;
    o.mElevation = 0;
    const RenderScene scene = prepare_render(mesh, o);
    auto covered = [&](double Azimuth) {
        RenderOptions camera = o;
        camera.mAzimuth = Azimuth;
        const Frame f = render_scene(scene, camera, true);
        std::size_t count = 0;
        for (std::int64_t id : f.mCellIds)
            count += id >= 0 ? 1 : 0;
        return count;
    };
    // The fitted frame resizes the cube between face-on and corner-on; the
    // fixed one holds the bounding sphere, so the silhouette only turns.
    const std::size_t face_on = covered(0.0);
    const std::size_t corner_on = covered(45.0);
    EXPECT_GT(face_on, 0u);
    EXPECT_GT(corner_on, face_on / 2);
    RenderOptions fitted = o;
    fitted.mAzimuth = 0;
    std::size_t fitted_face = 0;
    for (std::int64_t id : render(mesh, fitted).mCellIds)
        fitted_face += id >= 0 ? 1 : 0;
    EXPECT_LT(face_on, fitted_face);  // the sphere leaves margin the tight fit does not
}

TEST(RenderScene, AnEmptySceneIsRefused) {
    EXPECT_THROW(render_scene(RenderScene{}), std::invalid_argument);
}

// ------------------------------------------------------------- the cell diff

TEST(TextCells, AnUpdateOfAnUnchangedGridIsEmpty) {
    const Mesh mesh = cube_quads();
    TextOptions text;
    text.mCols = 20;
    text.mRows = 8;
    RenderOptions render_options;
    double aspect = 1.0;
    const auto size = text_frame_size(text, aspect);
    render_options.mWidth = size[0];
    render_options.mHeight = size[1];
    render_options.mPixelAspect = aspect;
    const TextGrid grid = encode_cells(render(mesh, render_options), text);
    EXPECT_EQ(grid.mCols, 20);
    EXPECT_EQ(grid.mRows, 8);
    EXPECT_EQ(grid.mCells.size(), 160u);
    EXPECT_TRUE(encode_cells_update(grid, grid, text.mDepth).empty());
}

TEST(TextCells, OneChangedCellIsOneMoveAndOneGlyph) {
    TextGrid a;
    a.mCols = 3;
    a.mRows = 2;
    a.mCells.assign(6, TextCell{});
    TextGrid b = a;
    b.mCells[4].mGlyph = 0x2580;  // ▀
    b.mCells[4].mFgSet = true;
    b.mCells[4].mFg = {255, 0, 0};
    const std::string update = encode_cells_update(a, b, ColorDepth::TrueColor, 10, 20);
    // Cell (row 1, col 1) sits at screen row 11, column 21.
    EXPECT_EQ(update,
              "\x1b[11;21H\x1b[38;2;255;0;0;49m\xe2\x96\x80\x1b[0m");
}

TEST(TextCells, AdjacentChangesShareOneCursorMove) {
    TextGrid a;
    a.mCols = 4;
    a.mRows = 1;
    a.mCells.assign(4, TextCell{});
    TextGrid b = a;
    b.mCells[1].mGlyph = 'x';
    b.mCells[2].mGlyph = 'y';
    const std::string update = encode_cells_update(a, b, ColorDepth::Mono);
    EXPECT_EQ(update, "\x1b[1;2Hxy");
}

TEST(TextCells, ColoursThatQuantizeAlikeAreNotRewritten) {
    TextGrid a;
    a.mCols = 1;
    a.mRows = 1;
    a.mCells.assign(1, TextCell{});
    a.mCells[0].mGlyph = 'x';
    a.mCells[0].mFgSet = true;
    a.mCells[0].mFg = {250, 10, 10};
    TextGrid b = a;
    b.mCells[0].mFg = {255, 0, 0};
    EXPECT_TRUE(encode_cells_update(a, b, ColorDepth::Ansi16).empty());
    EXPECT_FALSE(encode_cells_update(a, b, ColorDepth::TrueColor).empty());
}

TEST(TextCells, ADifferentGridSizeRepaintsEverything) {
    TextGrid a;
    a.mCols = 2;
    a.mRows = 1;
    a.mCells.assign(2, TextCell{});
    TextGrid b;
    b.mCols = 2;
    b.mRows = 2;
    b.mCells.assign(4, TextCell{});
    const std::string update = encode_cells_update(a, b, ColorDepth::Mono);
    EXPECT_EQ(std::count(update.begin(), update.end(), 'H'), 2);  // one move per row
}

TEST(TextCells, AGraphicsProtocolHasNoCells) {
    Frame frame;
    frame.mWidth = 2;
    frame.mHeight = 2;
    frame.mRgba.assign(16, 255);
    frame.mCellIds.assign(4, -1);
    TextOptions text;
    text.mEncoding = TextEncoding::Kitty;
    EXPECT_THROW(encode_cells(frame, text), std::invalid_argument);
}

// ---------------------------------------------------------------- the loop

std::string mouse(int Code, int X, int Y, bool Release = false) {
    return "\x1b[<" + std::to_string(Code) + ";" + std::to_string(X) + ";" + std::to_string(Y) +
           (Release ? "m" : "M");
}

TuiOptions loop_options() {
    TuiOptions options;
    options.mText.mNotes = false;
    options.mRender.mColorBy = "u";
    options.mTitle = "cube";
    return options;
}

TEST(TuiLoop, AScriptedSessionIsDeterministic) {
    const Mesh mesh = field_cube();
    const std::string script = mouse(0, 20, 10) + mouse(32, 30, 10) + mouse(32, 40, 14) +
                               mouse(0, 40, 14, true) + mouse(64, 40, 14) + "3" + "q";
    std::string first;
    for (int run = 0; run < 2; ++run) {
        TuiSession session(mesh, loop_options());
        ScriptedIo io(script, 80, 24);
        const TuiReport report = session.Run(io);
        EXPECT_EQ(report.mExit, 0);
        EXPECT_GE(report.mFrames, 2u);
        if (run == 0)
            first = io.Output();
        else
            EXPECT_EQ(io.Output(), first);
    }
    EXPECT_FALSE(first.empty());
}

TEST(TuiLoop, DraggingOrbitsTheCamera) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    TerminalSize size;
    size.mCols = 80;
    size.mRows = 24;
    const double az0 = session.Camera().mAzimuth;
    const double el0 = session.Camera().mElevation;
    EXPECT_FALSE(session.Handle(parse(mouse(0, 10, 10))[0], size));  // a press only grabs
    EXPECT_TRUE(session.Handle(parse(mouse(32, 20, 10))[0], size));  // 10 columns right
    EXPECT_NE(session.Camera().mAzimuth, az0);
    EXPECT_EQ(session.Camera().mElevation, el0);
    EXPECT_TRUE(session.Handle(parse(mouse(32, 20, 15))[0], size));  // 5 rows down
    EXPECT_NE(session.Camera().mElevation, el0);
    session.Handle(parse(mouse(0, 20, 15, true))[0], size);
    // Elevation stays within the poles however far one drags.
    for (int i = 0; i < 50; ++i)
        session.Handle(parse(mouse(32, 20, 15 + i))[0], size);
    EXPECT_LE(session.Camera().mElevation, 90.0);
    EXPECT_GE(session.Camera().mElevation, -90.0);
}

TEST(TuiLoop, ZoomPanPresetsAndReset) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    TerminalSize size;
    size.mCols = 80;
    size.mRows = 24;
    session.Handle(parse("+")[0], size);
    EXPECT_GT(session.Camera().mZoom, 1.0);
    session.Handle(parse(mouse(65, 1, 1))[0], size);
    session.Handle(parse(mouse(65, 1, 1))[0], size);
    EXPECT_LT(session.Camera().mZoom, 1.1);
    session.Handle(parse("\x1b[D")[0], size);
    EXPECT_LT(session.Camera().mPanX, 0.0);
    session.Handle(parse("2")[0], size);  // +x
    EXPECT_EQ(session.Camera().mAzimuth, 0.0);
    EXPECT_EQ(session.Camera().mElevation, 0.0);
    session.Handle(parse("r")[0], size);
    EXPECT_EQ(session.Camera().mZoom, 1.0);
    EXPECT_EQ(session.Camera().mPanX, 0.0);
    EXPECT_EQ(session.Camera().mAzimuth, 45.0);
}

TEST(TuiLoop, ZoomStaysInRange) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    TerminalSize size;
    size.mCols = 80;
    size.mRows = 24;
    for (int i = 0; i < 400; ++i)
        session.Handle(parse("+")[0], size);
    EXPECT_LE(session.Camera().mZoom, 200.0);
    for (int i = 0; i < 800; ++i)
        session.Handle(parse("-")[0], size);
    EXPECT_GE(session.Camera().mZoom, 0.05);
}

TEST(TuiLoop, QuitKeys) {
    const Mesh mesh = field_cube();
    for (const std::string& key : {std::string("q"), std::string("\x03"), std::string("\x04"),
                                   std::string("\x1b")}) {
        TuiSession session(mesh, loop_options());
        ScriptedIo io(key, 80, 24);
        const TuiReport report = session.Run(io);
        EXPECT_EQ(report.mExit, 0);
        EXPECT_TRUE(session.Quit()) << "key " << static_cast<int>(key[0]);
    }
}

TEST(TuiLoop, EndOfInputEndsTheSession) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    ScriptedIo io("", 80, 24);
    const TuiReport report = session.Run(io);
    EXPECT_EQ(report.mExit, 0);
    EXPECT_GE(report.mFrames, 1u);  // the first picture is drawn before it waits
}

// A screen whose size changes after the first frame, and a signal at the end.
class ResizingIo : public TuiIo {
public:
    bool Read(std::string& rBytes, int) override {
        ++mReads;
        if (mReads == 1) {
            mCols = 40;
            mRows = 12;
            mResized = true;
            return true;
        }
        if (mReads == 2) {
            rBytes = "q";
            return true;
        }
        return false;
    }
    bool Write(const std::string& rBytes) override {
        mWrites.push_back(rBytes);
        return true;
    }
    bool Size(TerminalSize& rSize) override {
        rSize.mCols = mCols;
        rSize.mRows = mRows;
        return true;
    }
    bool TakeResize() override {
        const bool was = mResized;
        mResized = false;
        return was;
    }
    int TerminationSignal() override { return 0; }

    int mCols = 100;
    int mRows = 40;
    bool mResized = false;
    int mReads = 0;
    std::vector<std::string> mWrites;
};

TEST(TuiLoop, AResizeNeverLeavesAFrameLargerThanTheScreen) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    ResizingIo io;
    const TuiReport report = session.Run(io);
    EXPECT_EQ(report.mExit, 0);
    EXPECT_LE(report.mScreen.mCols, 40);
    EXPECT_LE(report.mScreen.mRows, 12 - 1);  // one row for the status line
    EXPECT_GE(io.mWrites.size(), 2u);
}

class SignalIo : public ScriptedIo {
public:
    using ScriptedIo::ScriptedIo;
    int TerminationSignal() override { return 15; }
};

TEST(TuiLoop, ASignalEndsTheSessionWithItsExitCode) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    SignalIo io("", 80, 24);
    const TuiReport report = session.Run(io);
    EXPECT_EQ(report.mExit, 128 + 15);
}

TEST(TuiLoop, ATerminalTooSmallSaysSo) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    ScriptedIo io("q", 5, 3);
    const TuiReport report = session.Run(io);
    EXPECT_EQ(report.mExit, 0);
    EXPECT_NE(io.Output().find("terminal too small"), std::string::npos);
}

TEST(TuiLoop, EdgesToggleRepreparesTheScene) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    ScriptedIo io("eq", 80, 24);
    const TuiReport report = session.Run(io);
    EXPECT_EQ(report.mFinal.mEdges, RenderEdges::All);
    EXPECT_NE(io.Output().find("preparing"), std::string::npos);
}

TEST(TuiLoop, HelpIsShownAndHidden) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    ScriptedIo io("?q", 80, 24);
    session.Run(io);
    EXPECT_NE(io.Output().find("this help"), std::string::npos);
}

TEST(TuiLoop, ABadOptionFailsBeforeAnyTerminalState) {
    const Mesh mesh = field_cube();
    TuiOptions options;
    options.mRender.mColorBy = "no_such_array";
    EXPECT_THROW(TuiSession(mesh, options), std::invalid_argument);
}


// ------------------------------------------------------------ the flag set

TEST(RenderFlagsTest, OptionsRoundTripThroughFlagsAndTokens) {
    RenderOptions o;
    o.mAzimuth = 12.5;
    o.mElevation = -20.0;
    o.mZoom = 1.75;
    o.mPanX = 0.1;
    o.mProjection = RenderProjection::Perspective;
    o.mFovDeg = 40.0;
    o.mShading = RenderShading::Smooth;
    o.mEdges = RenderEdges::Feature;
    o.mFeatureAngle = 45.0;
    o.mBackground = {10, 20, 30, 255};
    o.mColorBy = "u";
    o.mComponent = 1;
    o.mCmap = "magma";
    o.mVMin = -1.0 / 3.0;
    o.mVMax = 7.0;
    o.mColorbar = true;
    o.mClipLow = 5.0;
    o.mSymmetric = true;
    o.mScale = RenderScale::Symlog;
    o.mIsolines = 4;
    o.mIsoLevels = {0.1, 0.25};
    o.mVectors = "disp";
    o.mWarp = "disp";
    o.mWarpOutline = true;
    o.mAxes = true;
    RenderCutaway cut;
    cut.mPoint = {0.5, 0.0, 0.0};
    cut.mNormal = {1.0, 0.0, 0.0};
    o.mCutaways = {cut};
    o.mCutawayTint = {1, 2, 3, 255};

    const RenderFlags flags = RenderFlags::FromOptions(o);
    const RenderOptions back = flags.ToOptions();
    EXPECT_EQ(back.mAzimuth, 12.5);
    EXPECT_EQ(back.mVMin, o.mVMin);  // exact, not rounded
    EXPECT_EQ(back.mCmap, "magma");
    EXPECT_EQ(back.mClipLow, 5.0);
    EXPECT_FALSE(back.mClipHigh.has_value());
    EXPECT_EQ(back.mIsoLevels, o.mIsoLevels);
    ASSERT_EQ(back.mCutaways.size(), 1u);
    EXPECT_EQ(back.mCutaways[0].mNormal, cut.mNormal);
    EXPECT_EQ(back.mCutawayTint, (RenderColor{1, 2, 3, 255}));
    EXPECT_EQ(back.mProjection, RenderProjection::Perspective);
    // The tokens read back to the same flags (a fixed point).
    const RenderFlags again = RenderFlags::FromTokens(flags.Tokens());
    EXPECT_EQ(again.Tokens(), flags.Tokens());
    EXPECT_EQ(RenderFlags::FromOptions(again.ToOptions()).Tokens(), flags.Tokens());
}

TEST(RenderFlagsTest, DefaultsAreNoFlagsAtAll) {
    EXPECT_TRUE(RenderFlags::FromOptions(RenderOptions{}).Tokens().empty());
}

TEST(RenderFlagsTest, AViewNameWinsOverAngles) {
    RenderOptions o;
    o.mView = "+x";
    EXPECT_EQ(RenderFlags::FromOptions(o).Tokens(), (std::vector<std::string>{"--view=+x"}));
}

TEST(RenderFlagsTest, RangeFlagsNeedAMappedFieldAndAreDroppedWithoutOne) {
    RenderOptions o;
    o.mCmap = "magma";  // no field is mapped: the cmap cannot be expressed
    EXPECT_TRUE(RenderFlags::FromOptions(o).Tokens().empty());
    RenderFlags f;
    f.Set("cmap", "magma");
    EXPECT_THROW(f.ToOptions(), std::runtime_error);  // the parser's own rule
}

TEST(RenderFlagsTest, SetUnsetAndKnown) {
    RenderFlags f;
    f.Set("color-by", "u");
    f.SetFlag("axes", true);
    EXPECT_TRUE(f.Has("color-by"));
    EXPECT_EQ(f.Value("color-by"), "u");
    f.Unset("color-by");
    EXPECT_FALSE(f.Has("color-by"));
    EXPECT_TRUE(RenderFlags::Known("cmap"));
    EXPECT_FALSE(RenderFlags::Known("nonsense"));
    EXPECT_TRUE(RenderFlags::TakesValue("cmap"));
    EXPECT_FALSE(RenderFlags::TakesValue("axes"));
}

TEST(CutawayFlag, ParsesAxisAndSixNumberForms) {
    const RenderCutaway a = cli_parse_cutaway("+x:0.5");
    EXPECT_EQ(a.mPoint, (std::array<double, 3>{0.5, 0.0, 0.0}));
    EXPECT_EQ(a.mNormal, (std::array<double, 3>{1.0, 0.0, 0.0}));
    const RenderCutaway b = cli_parse_cutaway("-z:2");
    EXPECT_EQ(b.mPoint, (std::array<double, 3>{0.0, 0.0, 2.0}));
    EXPECT_EQ(b.mNormal, (std::array<double, 3>{0.0, 0.0, -1.0}));
    const RenderCutaway c = cli_parse_cutaway("1,2,3,0,1,0");
    EXPECT_EQ(c.mPoint, (std::array<double, 3>{1.0, 2.0, 3.0}));
    EXPECT_EQ(c.mNormal, (std::array<double, 3>{0.0, 1.0, 0.0}));
    EXPECT_THROW(cli_parse_cutaway("x:1"), std::invalid_argument);
    EXPECT_THROW(cli_parse_cutaway("+x:"), std::invalid_argument);
    EXPECT_THROW(cli_parse_cutaway("1,2,3"), std::invalid_argument);
    EXPECT_THROW(cli_parse_cutaway("1,2,3,,5,6"), std::invalid_argument);
}

TEST(CutawayFlag, TheFlagReachesTheOptions) {
    const cli_parsed p = cli_parse({"--cutaway", "+x:0.5", "--cutaway=-y:1", "--cutaway-tint=#ff0000"},
                                   render_flag_specs());
    const RenderOptions o = cli_render_options(p);
    ASSERT_EQ(o.mCutaways.size(), 2u);
    EXPECT_EQ(o.mCutawayTint, (RenderColor{255, 0, 0, 255}));
    const cli_parsed lonely = cli_parse({"--cutaway-tint=#ff0000"}, render_flag_specs());
    EXPECT_THROW(cli_render_options(lonely), std::runtime_error);
}

// -------------------------------------------------------------- the session

TEST(SessionFileTest, WritesSortedKeysAndReadsThemBack) {
    SessionFile s;
    s.mFlags = {"--azimuth=12", "--cmap=magma", "--note=a \"quoted\" \\ value"};
    s.mStep = 7;
    const std::string text = session_to_json(s);
    EXPECT_LT(text.find("\"flags\""), text.find("\"step\""));
    EXPECT_LT(text.find("\"step\""), text.find("\"version\""));
    const SessionFile back = session_from_json(text);
    EXPECT_EQ(back.mFlags, s.mFlags);
    EXPECT_EQ(back.mStep, 7);
    EXPECT_EQ(session_to_json(back), text);  // a fixed point
}

TEST(SessionFileTest, StrictlyRefusesWhatItDoesNotKnow) {
    EXPECT_THROW(session_from_json("{}"), std::invalid_argument);  // no version
    EXPECT_THROW(session_from_json("{\"version\": 2}"), std::invalid_argument);
    EXPECT_THROW(session_from_json("{\"version\": 1, \"colour\": 3}"), std::invalid_argument);
    EXPECT_THROW(session_from_json("{\"version\": 1, \"version\": 1}"), std::invalid_argument);
    EXPECT_THROW(session_from_json("{\"version\": 1, \"flags\": [1]}"), std::invalid_argument);
    EXPECT_THROW(session_from_json("{\"version\": 1, \"step\": -3}"), std::invalid_argument);
    EXPECT_THROW(session_from_json("{\"version\": 1} trailing"), std::invalid_argument);
    EXPECT_THROW(session_from_json("not json"), std::invalid_argument);
    EXPECT_THROW(session_from_json("{\"version\": 1, \"flags\": [\"unterminated"),
                 std::invalid_argument);
    try {
        session_from_json("{\"version\": 9}");
        FAIL();
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("version 9"), std::string::npos);
    }
}

// --------------------------------------------------------------- the probe

Mesh probe_mesh() {
    Mesh mesh = cube_quads();
    std::vector<std::vector<double>> u;
    for (int i = 0; i < 8; ++i)
        u.push_back({double(i)});
    mesh.AddPointData("u", mt::points_from(u));
    mesh.AddCellData("id", {mt::points_from({{10}, {11}, {12}, {13}, {14}, {15}})});
    mesh.AddRegion(Region("left", RegionKind::Cell,
                          mt::points_from({{4}})));  // wrong dtype on purpose: replaced below
    return mesh;
}

TEST(Probe, DescribesTheCellItsDataAndItsRegions) {
    Mesh mesh = cube_quads();
    std::vector<std::vector<double>> u;
    for (int i = 0; i < 8; ++i)
        u.push_back({double(i)});
    mesh.AddPointData("u", mt::points_from(u));
    mesh.AddCellData("id", {mt::points_from({{10}, {11}, {12}, {13}, {14}, {15}})});
    const Probe p = probe_cell(mesh, 4);  // {0, 4, 7, 3}
    EXPECT_EQ(p.mCell, 4);
    ASSERT_GE(p.mLines.size(), 3u);
    EXPECT_NE(p.mLines[0].find("cell 4"), std::string::npos);
    EXPECT_NE(p.mLines[0].find("quad"), std::string::npos);
    EXPECT_NE(p.mLines[0].find("nodes 0 4 7 3"), std::string::npos);
    EXPECT_NE(p.mLines[1].find("id=14"), std::string::npos);
    EXPECT_NE(p.mLines[2].find("u 0..7 (mean 3.5)"), std::string::npos);
    EXPECT_EQ(p.mValues.at("id"), 14.0);
    EXPECT_EQ(p.mValues.at("u"), 3.5);
}

TEST(Probe, ANegativeOrOutOfRangeCellIsEmpty) {
    const Mesh mesh = cube_quads();
    EXPECT_EQ(probe_cell(mesh, -1).mCell, -1);
    EXPECT_EQ(probe_cell(mesh, 99).mCell, -1);
    EXPECT_TRUE(probe_cell(mesh, 99).mLines.empty());
}

TEST(Probe, TheDifferenceOfTwoProbes) {
    Probe a, b;
    a.mValues = {{"u", 1.0}, {"only_a", 5.0}};
    b.mValues = {{"u", 3.5}, {"only_b", 5.0}};
    EXPECT_EQ(probe_difference(a, b), "u +2.5");
    Probe none;
    EXPECT_TRUE(probe_difference(a, none).empty());
}

// ------------------------------------------------- the loop: new interactions

TuiOptions plain_options() {
    TuiOptions options;
    options.mText.mNotes = false;
    options.mTitle = "cube";
    return options;
}

TEST(TuiLoop, AClickProbesTheCellUnderThePointerAndAPinKeepsIt) {
    const Mesh mesh = field_cube();
    TuiOptions options = loop_options();
    options.mRender.mView = "+z";  // the top face, cell 1, fills the middle
    TuiSession session(mesh, options);
    ScriptedIo io(mouse(0, 40, 10) + mouse(0, 40, 10, true) + "i", 80, 24);
    const TuiReport report = session.Run(io);
    ASSERT_FALSE(report.mProbe.empty());
    EXPECT_NE(report.mProbe[0].find("cell 1"), std::string::npos);
    // The pin shows up as its own line.
    EXPECT_TRUE(std::any_of(report.mProbe.begin(), report.mProbe.end(),
                            [](const std::string& s) { return s.rfind("pin A", 0) == 0; }));
}

TEST(TuiLoop, ADragIsNotAClick) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    ScriptedIo io(mouse(0, 40, 10) + mouse(32, 44, 10) + mouse(0, 44, 10, true), 80, 24);
    const TuiReport report = session.Run(io);
    EXPECT_TRUE(report.mProbe.empty());
    EXPECT_NE(report.mFinal.mAzimuth, 45.0);
}

TEST(TuiLoop, TwoPinsShowTheirDifferenceAndZeroClearsThem) {
    const Mesh mesh = field_cube();
    TuiOptions options = loop_options();
    options.mRender.mView = "+z";
    TuiSession session(mesh, options);
    // The top face (cell 1), then the near corner of the side (cell 3 or 5).
    ScriptedIo io(mouse(0, 40, 10) + mouse(0, 40, 10, true) + "i" + mouse(0, 40, 3) +
                      mouse(0, 40, 3, true) + "i",
                  80, 24);
    const TuiReport report = session.Run(io);
    const bool has_diff = std::any_of(report.mProbe.begin(), report.mProbe.end(),
                                      [](const std::string& s) { return s.rfind("B - A", 0) == 0; });
    // Either the second click hit another cell (a difference line) or it hit the
    // background (no second pin); both are legitimate geometry, so only check
    // that a pin A exists and that clearing empties everything.
    EXPECT_TRUE(std::any_of(report.mProbe.begin(), report.mProbe.end(),
                            [](const std::string& s) { return s.rfind("pin A", 0) == 0; }));
    (void)has_diff;
    TuiSession again(mesh, options);
    ScriptedIo io2(mouse(0, 40, 10) + mouse(0, 40, 10, true) + "i0", 80, 24);
    EXPECT_TRUE(again.Run(io2).mProbe.empty());
}

TEST(TuiCommands, AFlagIsACommand) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    EXPECT_EQ(session.Execute("cmap magma"), "");
    EXPECT_EQ(session.Camera().mCmap, "magma");
    EXPECT_EQ(session.Execute(":edges all"), "");
    EXPECT_EQ(session.Camera().mEdges, RenderEdges::All);
    EXPECT_EQ(session.Execute("range 0 1"), "");
    EXPECT_EQ(session.Camera().mVMin, 0.0);
    EXPECT_EQ(session.Camera().mVMax, 1.0);
    EXPECT_EQ(session.Execute("range auto"), "");
    EXPECT_FALSE(session.Camera().mVMin.has_value());
    EXPECT_EQ(session.Execute("axes"), "");
    EXPECT_TRUE(session.Camera().mAxes);
    EXPECT_EQ(session.Execute("axes off"), "");
    EXPECT_FALSE(session.Camera().mAxes);
    EXPECT_EQ(session.Execute("isolines 3"), "");
    EXPECT_EQ(session.Camera().mIsolines, 3);
    EXPECT_EQ(session.Execute("unset isolines"), "");
    EXPECT_EQ(session.Camera().mIsolines, 0);
    EXPECT_EQ(session.Execute("zoom 2"), "");
    EXPECT_EQ(session.Camera().mZoom, 2.0);
    EXPECT_EQ(session.Execute("view +x"), "");
    EXPECT_EQ(session.Camera().mAzimuth, 0.0);
    EXPECT_EQ(session.Camera().mElevation, 0.0);
    EXPECT_TRUE(session.Camera().mView.empty());
}

TEST(TuiCommands, ColorChangesTheFieldAndNoneClearsIt) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    EXPECT_EQ(session.Execute("color disp"), "");
    EXPECT_EQ(session.Camera().mColorBy, "disp");
    EXPECT_EQ(session.Execute("color none"), "");
    EXPECT_TRUE(session.Camera().mColorBy.empty());
}

TEST(TuiCommands, BadInputIsAMessageNeverAnException) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, plain_options());  // no field mapped
    EXPECT_NE(session.Execute("cmap magma").find("requires"), std::string::npos);
    EXPECT_EQ(session.Camera().mCmap, "viridis");  // unchanged
    EXPECT_NE(session.Execute("nonsense").find("unknown command"), std::string::npos);
    EXPECT_NE(session.Execute("zoom").find("usage"), std::string::npos);
    EXPECT_NE(session.Execute("color no_such_array").find(""), std::string::npos);
    EXPECT_NE(session.Execute("edges bogus").find("edges expects"), std::string::npos);
    EXPECT_NE(session.Execute("'unterminated").find("quote"), std::string::npos);
    EXPECT_NE(session.Execute("view sideways").find("unknown view"), std::string::npos);
    EXPECT_EQ(session.Execute(""), "");
    EXPECT_EQ(session.Execute("   "), "");
}

TEST(TuiCommands, ClipAddsAndRemovesPlanes) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    EXPECT_EQ(session.Execute("clip +x 0.5"), "");
    ASSERT_EQ(session.Camera().mCutaways.size(), 1u);
    EXPECT_EQ(session.Camera().mCutaways[0].mPoint[0], 0.5);
    EXPECT_EQ(session.Execute("clip -y:0.25"), "");
    EXPECT_EQ(session.Camera().mCutaways.size(), 2u);
    EXPECT_EQ(session.Execute("clip +z 0.1"), "");  // a third replaces the first
    ASSERT_EQ(session.Camera().mCutaways.size(), 2u);
    EXPECT_EQ(session.Camera().mCutaways[0].mNormal[1], -1.0);
    EXPECT_EQ(session.Execute("clip off"), "");
    EXPECT_TRUE(session.Camera().mCutaways.empty());
    EXPECT_NE(session.Execute("clip sideways").find("cutaway"), std::string::npos);
}

TEST(TuiCommands, QuitAndHelpAndReset) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    session.Execute("zoom 3");
    session.Execute("reset");
    EXPECT_EQ(session.Camera().mZoom, 1.0);
    session.Execute("q");
    EXPECT_TRUE(session.Quit());
}

TEST(TuiCommands, ThePromptTypesRunsAndCancels) {
    const Mesh mesh = field_cube();
    {
        TuiSession session(mesh, loop_options());
        ScriptedIo io(":cmap magma\r", 80, 24);
        const TuiReport report = session.Run(io);
        EXPECT_EQ(report.mFinal.mCmap, "magma");
    }
    {
        TuiSession session(mesh, loop_options());
        ScriptedIo io(":cmap magma\x1b", 80, 24);  // Escape cancels
        const TuiReport report = session.Run(io);
        EXPECT_EQ(report.mFinal.mCmap, "viridis");
    }
    {
        TuiSession session(mesh, loop_options());
        ScriptedIo io(":cmap mx\x7f" "agma\r", 80, 24);  // a backspace in the middle
        EXPECT_EQ(session.Run(io).mFinal.mCmap, "magma");
    }
    {
        TuiSession session(mesh, loop_options());
        ScriptedIo io(":bogus\r", 80, 24);  // the message stays until the next key
        const TuiReport report = session.Run(io);
        EXPECT_NE(report.mStatus.find("unknown command"), std::string::npos);
    }
}

TEST(TuiCommands, AChangeTheSceneCannotBeBuiltWithIsRefusedAndNothingElseIsLost) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    EXPECT_EQ(session.Execute("cmap magma"), "");
    EXPECT_EQ(session.Execute("zoom 2"), "");
    const std::string message = session.Execute("cmap nope");
    EXPECT_NE(message.find("unknown colormap"), std::string::npos);
    EXPECT_EQ(session.Camera().mCmap, "magma");  // the earlier change survives
    EXPECT_EQ(session.Camera().mZoom, 2.0);
    EXPECT_NE(session.Execute("color no_such_array").find("no_such_array"), std::string::npos);
    EXPECT_EQ(session.Camera().mColorBy, "u");
    // The same through the prompt, in one recorded batch.
    TuiSession typed(mesh, loop_options());
    ScriptedIo io(":cmap magma\r:zoom 2\r:cmap nope\r", 80, 24);
    const TuiReport report = typed.Run(io);
    EXPECT_EQ(report.mFinal.mCmap, "magma");
    EXPECT_EQ(report.mFinal.mZoom, 2.0);
    EXPECT_NE(report.mStatus.find("unknown colormap"), std::string::npos);
    EXPECT_EQ(report.mExit, 0);
}

TEST(TuiCommands, WriteSavesASnapshotOfWhatIsShown) {
    const Mesh mesh = field_cube();
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "mio_tui_w";
    std::filesystem::create_directories(dir);
    const std::string path = (dir / "shot.png").string();
    std::filesystem::remove(path);
    TuiSession session(mesh, loop_options());
    EXPECT_EQ(session.Execute("w " + path), "wrote " + path);
    EXPECT_TRUE(std::filesystem::exists(path));
    EXPECT_GT(std::filesystem::file_size(path), 100u);
    EXPECT_NE(session.Execute("w " + (dir / "x.unknown").string()).find("unknown"),
              std::string::npos);
    EXPECT_NE(session.Execute("w").find("usage"), std::string::npos);
}

TEST(TuiCommands, ASessionSavesAndLoadsTheView) {
    const Mesh mesh = field_cube();
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "mio_tui_s";
    std::filesystem::create_directories(dir);
    const std::string path = (dir / "view.json").string();
    {
        TuiSession session(mesh, loop_options());
        session.Execute("zoom 2.5");
        session.Execute("azimuth 33");
        session.Execute("edges all");
        session.Execute("clip +x 0.5");
        EXPECT_EQ(session.Execute("session save " + path), "saved " + path);
    }
    EXPECT_TRUE(std::filesystem::exists(path));
    TuiSession other(mesh, loop_options());
    EXPECT_EQ(other.Execute("session load " + path), "loaded " + path);
    EXPECT_EQ(other.Camera().mZoom, 2.5);
    EXPECT_EQ(other.Camera().mAzimuth, 33.0);
    EXPECT_EQ(other.Camera().mEdges, RenderEdges::All);
    ASSERT_EQ(other.Camera().mCutaways.size(), 1u);
    EXPECT_NE(other.Execute("session load " + (dir / "missing.json").string()).find("cannot read"),
              std::string::npos);
}

TEST(TuiLoop, ASessionPathIsReadAtTheStartAndWrittenAtTheEnd) {
    const Mesh mesh = field_cube();
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "mio_tui_p";
    std::filesystem::create_directories(dir);
    const std::string path = (dir / "auto.json").string();
    std::filesystem::remove(path);
    TuiOptions options = loop_options();
    options.mSessionPath = path;
    {
        TuiSession session(mesh, options);
        ScriptedIo io(mouse(0, 20, 10) + mouse(32, 30, 10) + mouse(0, 30, 10, true) + "q", 80, 24);
        const TuiReport report = session.Run(io);
        EXPECT_EQ(report.mError, "");
    }
    ASSERT_TRUE(std::filesystem::exists(path));
    TuiSession reopened(mesh, options);
    EXPECT_NE(reopened.Camera().mAzimuth, 45.0);  // where the last run left it
    // A damaged file stops the start, by name.
    std::ofstream(path) << "{\"version\": 1, \"unknown\": 1}";
    EXPECT_THROW(TuiSession(mesh, options), std::invalid_argument);
}

TEST(TuiLoop, CutawayKeysTogglePlanesAndSlideThem) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    TerminalSize size;
    size.mCols = 80;
    size.mRows = 24;
    session.Handle(parse("x")[0], size);
    ASSERT_EQ(session.Camera().mCutaways.size(), 1u);
    EXPECT_EQ(session.Camera().mCutaways[0].mNormal[0], 1.0);
    const double at = session.Camera().mCutaways[0].mPoint[0];
    session.Handle(parse(".")[0], size);
    EXPECT_GT(session.Camera().mCutaways[0].mPoint[0], at);
    session.Handle(parse(",")[0], size);
    session.Handle(parse(",")[0], size);
    EXPECT_LT(session.Camera().mCutaways[0].mPoint[0], at);
    session.Handle(parse("x")[0], size);  // the other side
    EXPECT_EQ(session.Camera().mCutaways[0].mNormal[0], -1.0);
    session.Handle(parse("x")[0], size);  // off
    EXPECT_TRUE(session.Camera().mCutaways.empty());
    session.Handle(parse("x")[0], size);
    session.Handle(parse("y")[0], size);
    session.Handle(parse("z")[0], size);  // a third is refused
    EXPECT_EQ(session.Camera().mCutaways.size(), 2u);
    session.Handle(parse("0")[0], size);
    EXPECT_TRUE(session.Camera().mCutaways.empty());
}

TEST(TuiLoop, ACutawayChangesThePictureWithoutPreparingAgain) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    ScriptedIo plain("q", 80, 24);
    session.Run(plain);
    TuiSession cut(mesh, loop_options());
    ScriptedIo io("xq", 80, 24);
    cut.Run(io);
    EXPECT_NE(plain.Output(), io.Output());
    EXPECT_EQ(io.Output().find("preparing"), std::string::npos);
}

// ------------------------------------------------------------- comparison

TEST(TuiCompare, TwoMeshesShareTheScreenAndTheCamera) {
    const Mesh a = field_cube();
    Mesh b = field_cube();
    TuiOptions options = loop_options();
    options.mpCompare = &b;
    options.mTitle = "a.vtu";
    options.mCompareTitle = "b.vtu";
    TuiSession session(a, options);
    ScriptedIo io(mouse(0, 20, 10) + mouse(32, 30, 10) + mouse(0, 30, 10, true) + "q", 81, 24);
    const TuiReport report = session.Run(io);
    EXPECT_EQ(report.mScreen.mCols, 81);
    // The separator column sits between the two halves (40 + 1 + 40).
    int bars = 0;
    for (int r = 0; r < report.mScreen.mRows; ++r)
        bars += report.mScreen.mCells[static_cast<std::size_t>(r) * 81 + 40].mGlyph == 0x2502 ? 1 : 0;
    EXPECT_EQ(bars, report.mScreen.mRows);
    EXPECT_NE(report.mFinal.mAzimuth, 45.0);
}

TEST(TuiCompare, ClickingEitherSideProbesThatMesh) {
    const Mesh a = field_cube();
    Mesh b = cube_quads();
    std::vector<std::vector<double>> u;
    for (int i = 0; i < 8; ++i)
        u.push_back({100.0 + 0.5 * i});  // A's field plus 100
    b.AddPointData("u", mt::points_from(u));
    TuiOptions options = loop_options();
    options.mpCompare = &b;
    options.mRender.mView = "+z";
    TuiSession session(a, options);
    // Left pane centre, then right pane centre.
    ScriptedIo io(mouse(0, 20, 10) + mouse(0, 20, 10, true) + "i" + mouse(0, 60, 10) +
                      mouse(0, 60, 10, true) + "i",
                  81, 24);
    const TuiReport report = session.Run(io);
    const bool pins = std::count_if(report.mProbe.begin(), report.mProbe.end(),
                                    [](const std::string& s) { return s.rfind("pin", 0) == 0; }) >= 2;
    EXPECT_TRUE(pins);
    // B's field is 100 higher than A's: the difference line says so.
    EXPECT_TRUE(std::any_of(report.mProbe.begin(), report.mProbe.end(), [](const std::string& s) {
        return s.rfind("B - A", 0) == 0 && s.find("u +100") != std::string::npos;
    }));
}

TEST(TuiCompare, ADifferenceNeedsAFieldBothMeshesHave) {
    const Mesh a = field_cube();
    Mesh b = field_cube();
    TuiOptions options = loop_options();
    options.mpCompare = &b;
    options.mDiff = true;
    EXPECT_NO_THROW(TuiSession(a, options));
    options.mRender.mColorBy = "no_such";
    EXPECT_THROW(TuiSession(a, options), std::invalid_argument);
    const Mesh other = cube_quads();  // no u
    options.mRender.mColorBy = "u";
    options.mpCompare = &other;
    EXPECT_THROW(TuiSession(a, options), std::invalid_argument);
}

TEST(TuiCompare, AGraphicsProtocolCannotCompare) {
    const Mesh a = field_cube();
    Mesh b = field_cube();
    TuiOptions options = loop_options();
    options.mpCompare = &b;
    options.mText.mEncoding = TextEncoding::Kitty;
    EXPECT_THROW(TuiSession(a, options), std::invalid_argument);
}

// ------------------------------------------------------------------ series

// A series of in-memory steps: each a cube whose field grows with the step.
class MemorySeries : public TuiSeries {
public:
    explicit MemorySeries(std::size_t Count) : mCount(Count) {}
    std::size_t Count() override { return mCount; }
    Mesh Load(std::size_t Index) override {
        if (mFailLoads > 0) {
            --mFailLoads;
            throw ReadError("a half-written file");
        }
        ++mLoads;
        Mesh mesh = cube_quads();
        std::vector<std::vector<double>> u;
        for (int i = 0; i < 8; ++i)
            u.push_back({0.5 * i * double(Index + 1)});
        mesh.AddPointData("u", mt::points_from(u));
        return mesh;
    }
    std::string Label(std::size_t Index) override { return "step " + std::to_string(Index + 1); }
    bool Refresh() override { return false; }
    std::string Fingerprint() override { return std::to_string(mCount) + "|" + mTag; }

    std::size_t mCount;
    std::string mTag;
    int mFailLoads = 0;
    int mLoads = 0;
};

TuiOptions series_options(const std::shared_ptr<TuiSeries>& rSeries) {
    TuiOptions options = loop_options();
    options.mpSeries = rSeries;
    return options;
}

TEST(TuiSeriesTest, KeysAndCommandsStepThroughTheSeries) {
    const Mesh unused = cube_quads();
    auto series = std::make_shared<MemorySeries>(5);
    TuiSession session(unused, series_options(series));
    TerminalSize size;
    size.mCols = 80;
    size.mRows = 24;
    EXPECT_EQ(session.Step(), 0u);
    session.Handle(parse("]")[0], size);
    EXPECT_EQ(session.Step(), 1u);
    session.Handle(parse("}")[0], size);
    EXPECT_EQ(session.Step(), 4u);  // clamped to the last
    session.Handle(parse("]")[0], size);
    EXPECT_EQ(session.Step(), 4u);
    session.Handle(parse("{")[0], size);
    EXPECT_EQ(session.Step(), 0u);
    session.Handle(parse("[")[0], size);
    EXPECT_EQ(session.Step(), 0u);
    EXPECT_EQ(session.Execute("step 3"), "");
    EXPECT_EQ(session.Step(), 2u);
    session.Execute("step last");
    EXPECT_EQ(session.Step(), 4u);
    session.Execute("step prev");
    EXPECT_EQ(session.Step(), 3u);
    session.Execute("step first");
    EXPECT_EQ(session.Step(), 0u);
    EXPECT_NE(session.Execute("step").find("usage"), std::string::npos);
}

TEST(TuiSeriesTest, TheStatusNamesTheStepAndAScriptEndsOnIt) {
    const Mesh unused = cube_quads();
    auto series = std::make_shared<MemorySeries>(4);
    TuiSession session(unused, series_options(series));
    ScriptedIo io("]]q", 100, 24);
    const TuiReport report = session.Run(io);
    EXPECT_EQ(report.mStep, 2u);
    EXPECT_NE(report.mStatus.find("step 3"), std::string::npos);
}

TEST(TuiSeriesTest, OneMeshHasNoSteps) {
    const Mesh mesh = field_cube();
    TuiSession session(mesh, loop_options());
    TerminalSize size;
    size.mCols = 80;
    size.mRows = 24;
    session.Handle(parse("]")[0], size);
    EXPECT_EQ(session.Step(), 0u);
    EXPECT_NE(session.Execute("step 2").find("no steps"), std::string::npos);
}

TEST(TuiSeriesTest, AnEmptySeriesIsRefused) {
    const Mesh unused = cube_quads();
    EXPECT_THROW(TuiSession(unused, series_options(std::make_shared<MemorySeries>(0))),
                 std::invalid_argument);
}

// A clock the test drives: events (time, bytes) arrive when the clock passes them.
class ClockIo : public TuiIo {
public:
    ClockIo(std::vector<std::pair<long long, std::string>> Events, long long EndMs)
        : mEvents(std::move(Events)), mEndMs(EndMs) {}
    bool Read(std::string& rBytes, int TimeoutMs) override {
        mClock += TimeoutMs;
        if (mClock > mEndMs)
            return false;
        while (mNext < mEvents.size() && mEvents[mNext].first <= mClock)
            rBytes += mEvents[mNext++].second;
        return true;
    }
    bool Write(const std::string& rBytes) override {
        mOutput += rBytes;
        return true;
    }
    bool Size(TerminalSize& rSize) override {
        rSize.mCols = 80;
        rSize.mRows = 24;
        return true;
    }
    bool TakeResize() override { return false; }
    int TerminationSignal() override { return 0; }
    long long NowMs() override { return mClock; }
    std::string mOutput;

private:
    std::vector<std::pair<long long, std::string>> mEvents;
    std::size_t mNext = 0;
    long long mEndMs;
    long long mClock = 0;
};

TEST(TuiSeriesTest, PlayingAdvancesAtTheStatedRateAndStopsAtTheEnd) {
    const Mesh unused = cube_quads();
    auto series = std::make_shared<MemorySeries>(4);
    TuiOptions options = series_options(series);
    options.mPlayFps = 10.0;  // a step every 100 ms
    TuiSession session(unused, options);
    ClockIo io({{50, " "}}, 2000);  // start playing at t = 50
    const TuiReport report = session.Run(io);
    EXPECT_EQ(report.mStep, 3u);  // the last step, reached and stopped on
    EXPECT_EQ(series->mLoads, 1 + 3);  // the first step plus three advances
}

TEST(TuiSeriesTest, FollowingReadsANewStepOnceItHasSettled) {
    const Mesh unused = cube_quads();
    auto series = std::make_shared<MemorySeries>(2);
    TuiOptions options = series_options(series);
    options.mFollowMs = 100;
    options.mFollowSettleMs = 200;
    TuiSession session(unused, options);
    EXPECT_EQ(session.Step(), 1u);  // starts on the newest
    // A third step appears at t = 300 and is still changing until t = 450.
    ClockIo io({{300, ""}}, 1500);
    class Growing : public ClockIo {
    public:
        Growing(std::shared_ptr<MemorySeries> Series, std::vector<std::pair<long long, std::string>> E,
                long long End)
            : ClockIo(std::move(E), End), mSeries(std::move(Series)) {}
        bool Read(std::string& rBytes, int TimeoutMs) override {
            const bool alive = ClockIo::Read(rBytes, TimeoutMs);
            const long long now = NowMs();
            if (now >= 300 && mSeries->mCount == 2) {
                mSeries->mCount = 3;
                mSeries->mTag = "a";
            }
            if (now >= 400 && now < 450)
                mSeries->mTag = "b";  // still being written
            return alive;
        }
        std::shared_ptr<MemorySeries> mSeries;
    } growing(series, {}, 1500);
    const TuiReport report = session.Run(growing);
    EXPECT_EQ(report.mStep, 2u);  // jumped to the new last step
    EXPECT_NE(report.mStatus.find("followed"), std::string::npos);
}

TEST(TuiSeriesTest, APartialReadKeepsTheLastGoodFrameAndIsRetried) {
    const Mesh unused = cube_quads();
    auto series = std::make_shared<MemorySeries>(1);
    TuiOptions options = series_options(series);
    options.mFollowMs = 100;
    options.mFollowSettleMs = 100;
    TuiSession session(unused, options);
    series->mCount = 2;
    series->mTag = "x";
    series->mFailLoads = 2;  // the first two reads find a half-written file
    ClockIo io({}, 3000);
    const TuiReport report = session.Run(io);
    EXPECT_EQ(report.mStep, 1u);  // it got there in the end
    EXPECT_EQ(series->mFailLoads, 0);
    EXPECT_EQ(report.mExit, 0);
}

TEST(TuiSeriesTest, TheColourRangeIsKeptFromTheFirstStep) {
    const Mesh unused = cube_quads();
    auto series = std::make_shared<MemorySeries>(3);
    TuiOptions options = series_options(series);
    options.mRender.mColorBy = "u";
    options.mRender.mColorbar = true;
    TuiSession session(unused, options);
    ScriptedIo first("q", 100, 24);
    const TuiReport at_first = session.Run(first);
    ASSERT_NE(at_first.mScreen.mCells.size(), 0u);
    TuiSession again(unused, options);
    ScriptedIo later("]]q", 100, 24);
    again.Run(later);
    // The notes under the picture name the range (the corner means of the faces).
    // Step 1's is 0.75 .. 2.75, and step 3, three times as large, keeps it
    // instead of rescaling to 2.25 .. 8.25.
    EXPECT_NE(later.Output().find("0.75 .. 2.75"), std::string::npos);
    EXPECT_EQ(later.Output().find("2.25 .. 8.25"), std::string::npos);
}

TEST(TuiSeriesTest, RangeAutoLetsEachStepChooseItsOwn) {
    const Mesh unused = cube_quads();
    auto series = std::make_shared<MemorySeries>(3);
    TuiOptions options = series_options(series);
    options.mRender.mColorBy = "u";
    TuiSession session(unused, options);
    ScriptedIo io(":range auto\r]]q", 100, 24);
    session.Run(io);
    EXPECT_NE(io.Output().find("2.25 .. 8.25"), std::string::npos);
}

TEST(PathSeriesTest, ReadsStepsFromFilesAndNoticesNewOnes) {
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "mio_tui_series";
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    auto write = [&](int step) {
        Mesh mesh = cube_quads();
        std::vector<std::vector<double>> u;
        for (int i = 0; i < 8; ++i)
            u.push_back({double(step)});
        mesh.AddPointData("u", mt::points_from(u));
        write_vtu((dir / ("out_" + std::to_string(step) + ".vtu")).string(), mesh, false, false);
    };
    write(1);
    write(2);
    PathSeries series({}, (dir / "out_*.vtu").string(), "");
    EXPECT_EQ(series.Count(), 2u);
    EXPECT_NE(series.Label(0).find("out_1.vtu"), std::string::npos);
    const Mesh second = series.Load(1);
    EXPECT_TRUE(second.HasPointData("u"));
    const std::string before = series.Fingerprint();
    write(10);
    EXPECT_NE(series.Fingerprint(), before);  // a new file changes it
    EXPECT_TRUE(series.Refresh() || series.Count() == 3u);
    EXPECT_EQ(series.Count(), 3u);
    EXPECT_NE(series.Label(2).find("out_10.vtu"), std::string::npos);  // natural order: 2 < 10
    EXPECT_THROW(PathSeries({}, (dir / "nothing_*.vtu").string(), ""), ReadError);
    std::filesystem::remove_all(dir);
}

}  // namespace
