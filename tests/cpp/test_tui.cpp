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
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "meshioplusplus/operations/render.hpp"
#include "mesh_fixtures.hpp"

// CLI layer (compiled into the test binary like view_payload)
#include "../../src/cpp/cli/tui/input.hpp"
#include "../../src/cpp/cli/tui/loop.hpp"

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

}  // namespace
