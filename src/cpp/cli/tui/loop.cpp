// SPDX-License-Identifier: MIT
/// @file loop.cpp
/// @brief Implementation of loop.hpp.

#include "loop.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <utility>

#include "meshioplusplus/detail/fast_number.hpp"

namespace meshioplusplus::cli::tui {

namespace {

constexpr double kZoomStep = 1.1;
constexpr double kPanStep = 0.05;
constexpr double kMinZoom = 0.05;
constexpr double kMaxZoom = 200.0;

// The named views, the same camera directions `RenderOptions::mView` selects
// (the library keeps its own copy private); the loop turns a drag into an
// angle change from wherever a preset left the camera.
struct Preset {
    const char* mName;
    double mAzimuth;
    double mElevation;
};
const Preset kPresets[] = {{"iso", 45.0, 35.264389682754654}, {"+x", 0.0, 0.0},
                           {"-x", 180.0, 0.0},                {"+y", 90.0, 0.0},
                           {"-y", -90.0, 0.0},                {"+z", -90.0, 90.0},
                           {"-z", -90.0, -90.0}};

bool is_cell_encoding(TextEncoding Encoding) {
    return Encoding != TextEncoding::Kitty && Encoding != TextEncoding::ITerm2 &&
           Encoding != TextEncoding::Sixel;
}

// `Text` cut or padded to exactly `Cols` columns (one column per code point:
// the status and notes are ASCII, and a wide glyph in a title only shifts the
// padding).
std::string fit_columns(const std::string& rText, int Cols) {
    std::string out;
    int used = 0;
    for (std::size_t i = 0; i < rText.size() && used < Cols;) {
        const unsigned char c = static_cast<unsigned char>(rText[i]);
        const std::size_t len = c < 0x80 ? 1 : c >= 0xF0 ? 4 : c >= 0xE0 ? 3 : 2;
        out.append(rText, i, std::min(len, rText.size() - i));
        i += len;
        ++used;
    }
    if (used < Cols)
        out.append(static_cast<std::size_t>(Cols - used), ' ');
    return out;
}

std::string at_row(int Row) {
    return "\x1b[" + std::to_string(Row) + ";1H";
}

std::string number(double Value, int Decimals) {
    char buffer[32];
    meshioplusplus::detail::snprintf_c(buffer, sizeof(buffer), "%.*f", Decimals, Value);
    return buffer;
}

const char* edges_name(RenderEdges Edges) {
    return Edges == RenderEdges::All ? "all" : Edges == RenderEdges::Feature ? "feature" : "off";
}

const char* const kHelp[] = {
    " meshio++ tui                                          ",
    " drag            orbit            wheel / + -   zoom    ",
    " arrows          pan              r / Home      reset   ",
    " 1 .. 7          iso +x -x +y -y +z -z                  ",
    " p  perspective  a  axes          c  colour bar         ",
    " e  edges        s  shading       b  scale bar          ",
    " ?  this help    q / Esc / Ctrl-C   quit                ",
};

}  // namespace

TuiSession::TuiSession(const Mesh& rMesh, TuiOptions Options)
    : mrMesh(rMesh), mOptions(std::move(Options)) {
    mInitial = mOptions.mRender;
    if (!mInitial.mView.empty()) {
        // Start the camera at the named view, then drag from there.
        for (const Preset& preset : kPresets)
            if (mInitial.mView == preset.mName) {
                mInitial.mAzimuth = preset.mAzimuth;
                mInitial.mElevation = preset.mElevation;
            }
        mInitial.mView.clear();
    }
    mCamera = mInitial;
    mPreparedSmooth = mCamera.mShading == RenderShading::Smooth;
    mScene = prepare_render(mrMesh, mCamera);
}

bool TuiSession::ResetCamera() {
    // The camera returns to where it started; edges and shading keep what the
    // user set, since the prepared scene was built with them.
    const RenderEdges edges = mCamera.mEdges;
    const RenderShading shading = mCamera.mShading;
    mCamera = mInitial;
    mCamera.mEdges = edges;
    mCamera.mShading = shading;
    return true;
}

void TuiSession::Reprepare() {
    mScene = prepare_render(mrMesh, mCamera);
    mPreparedSmooth = mCamera.mShading == RenderShading::Smooth;
    mReprepare = false;
}

bool TuiSession::Handle(const InputEvent& rEvent, const TerminalSize& rSize) {
    const int cols = std::max(1, rSize.mCols);
    const int pic_rows = std::max(1, rSize.mRows - 1 - mNoteRows);
    auto zoom_by = [&](double Factor) {
        mCamera.mZoom = std::clamp(mCamera.mZoom * Factor, kMinZoom, kMaxZoom);
        return true;
    };
    auto orbit = [&](int Dx, int Dy) {
        mCamera.mAzimuth -= Dx * 360.0 / cols;
        mCamera.mAzimuth = std::remainder(mCamera.mAzimuth, 360.0);
        mCamera.mElevation = std::clamp(mCamera.mElevation + Dy * 180.0 / pic_rows, -90.0, 90.0);
        return Dx != 0 || Dy != 0;
    };
    auto pan = [&](double Dx, double Dy) {
        mCamera.mPanX += Dx;
        mCamera.mPanY += Dy;
        return Dx != 0.0 || Dy != 0.0;
    };

    if (rEvent.mKind == EventKind::Mouse) {
        switch (rEvent.mAction) {
            case MouseAction::Press:
                mDragging = true;
                mDragButton = rEvent.mButton;
                mLastX = rEvent.mX;
                mLastY = rEvent.mY;
                return false;
            case MouseAction::Release:
                mDragging = false;
                return false;
            case MouseAction::Drag: {
                if (!mDragging) {
                    mDragging = true;
                    mDragButton = rEvent.mButton;
                    mLastX = rEvent.mX;
                    mLastY = rEvent.mY;
                    return false;
                }
                const int dx = rEvent.mX - mLastX;
                const int dy = rEvent.mY - mLastY;
                mLastX = rEvent.mX;
                mLastY = rEvent.mY;
                if (mDragButton == 0 && !rEvent.mShift)
                    return orbit(dx, dy);
                return pan(static_cast<double>(dx) / cols, -static_cast<double>(dy) / pic_rows);
            }
            case MouseAction::WheelUp:
                return zoom_by(kZoomStep);
            case MouseAction::WheelDown:
                return zoom_by(1.0 / kZoomStep);
            case MouseAction::Move:
                return false;
        }
        return false;
    }
    if (rEvent.mKind == EventKind::Paste)
        return false;

    mMessage.clear();
    if (rEvent.mKey == KeyCode::Escape || rEvent.IsChar('q') || rEvent.IsCtrl('c') ||
        rEvent.IsCtrl('d')) {
        mQuit = true;
        return false;
    }
    switch (rEvent.mKey) {
        case KeyCode::Up:
            return pan(0.0, kPanStep);
        case KeyCode::Down:
            return pan(0.0, -kPanStep);
        case KeyCode::Left:
            return pan(-kPanStep, 0.0);
        case KeyCode::Right:
            return pan(kPanStep, 0.0);
        case KeyCode::PageUp:
            return zoom_by(kZoomStep);
        case KeyCode::PageDown:
            return zoom_by(1.0 / kZoomStep);
        case KeyCode::Home:
            return ResetCamera();
        default:
            break;
    }
    if (rEvent.mKey != KeyCode::Char || rEvent.mCtrl || rEvent.mAlt)
        return false;
    const std::uint32_t ch = rEvent.mChar;
    if (ch >= '1' && ch <= '7') {
        const Preset& preset = kPresets[ch - '1'];
        mCamera.mAzimuth = preset.mAzimuth;
        mCamera.mElevation = preset.mElevation;
        mCamera.mRoll = 0.0;
        mMessage = std::string("view ") + preset.mName;
        return true;
    }
    switch (ch) {
        case '+':
        case '=':
            return zoom_by(kZoomStep);
        case '-':
        case '_':
            return zoom_by(1.0 / kZoomStep);
        case 'r':
            return ResetCamera();
        case 'p':
            mCamera.mProjection = mCamera.mProjection == RenderProjection::Perspective
                                      ? RenderProjection::Orthographic
                                      : RenderProjection::Perspective;
            return true;
        case 'a':
            mCamera.mAxes = !mCamera.mAxes;
            return true;
        case 'b':
            mCamera.mScaleBar = !mCamera.mScaleBar;
            return true;
        case 'c':
            mCamera.mColorbar = !mCamera.mColorbar;
            return true;
        case 'e':
            mCamera.mEdges = mCamera.mEdges == RenderEdges::None  ? RenderEdges::All
                             : mCamera.mEdges == RenderEdges::All ? RenderEdges::Feature
                                                                  : RenderEdges::None;
            mReprepare = true;  // edge lines are part of the prepared scene
            mMessage = std::string("edges ") + edges_name(mCamera.mEdges);
            return true;
        case 's':
            mCamera.mShading = mCamera.mShading == RenderShading::Flat     ? RenderShading::Smooth
                               : mCamera.mShading == RenderShading::Smooth ? RenderShading::None
                                                                           : RenderShading::Flat;
            if (mCamera.mShading == RenderShading::Smooth && !mPreparedSmooth)
                mReprepare = true;  // smooth normals are computed when preparing
            mMessage = "shading " + std::string(mCamera.mShading == RenderShading::Smooth ? "smooth"
                                                : mCamera.mShading == RenderShading::Flat ? "flat"
                                                                                          : "none");
            return true;
        case '?':
        case 'h':
            mHelp = !mHelp;
            return true;
        default:
            return false;
    }
}

void TuiSession::Draw(TuiIo& rIo, const TerminalSize& rSize, bool Full) {
    const int cols = rSize.mCols;
    const int rows = rSize.mRows;
    std::string out;
    if (Full)
        out += "\x1b[0m\x1b[2J";
    if (mReprepare) {
        // The preparation can take seconds on a large volume: say so first.
        rIo.Write(at_row(rows) + "\x1b[7m" + fit_columns(" preparing...", cols) + "\x1b[0m");
        try {
            Reprepare();
        } catch (const std::exception& rErr) {
            mReprepare = false;
            mMessage = rErr.what();
        }
    }
    if (cols < 8 || rows < 4) {
        out += "\x1b[2J" + at_row(1) + "terminal too small";
        rIo.Write(out);
        mPrevious = TextGrid{};
        return;
    }

    const bool cells = is_cell_encoding(mOptions.mText.mEncoding);
    Frame frame;
    TextGrid grid;
    std::string graphics;
    int pic_rows = std::max(1, rows - 1 - mNoteRows);
    for (int pass = 0; pass < 2; ++pass) {
        TextOptions text = mOptions.mText;
        text.mCols = cols;
        text.mRows = pic_rows;
        if (mOptions.mCellAspectFromTerminal && rSize.mPixelWidth > 0 && rSize.mPixelHeight > 0) {
            text.mCellPixelWidth = std::max(1, rSize.mPixelWidth / cols);
            text.mCellPixelHeight = std::max(1, rSize.mPixelHeight / rows);
            text.mCellAspect = static_cast<double>(rSize.mPixelHeight * cols) /
                               static_cast<double>(rSize.mPixelWidth * rows);
        }
        RenderOptions r = mCamera;
        double aspect = 1.0;
        const std::array<int, 2> size = text_frame_size(text, aspect);
        r.mWidth = size[0];
        r.mHeight = size[1];
        r.mPixelAspect = aspect;
        try {
            frame = render_scene(mScene, r, mOptions.mFixedFit);
            if (cells)
                grid = encode_cells(frame, text);
            else {
                TextOptions plain = text;
                plain.mNotes = false;
                graphics = encode_text(frame, plain);
            }
        } catch (const std::exception& rErr) {
            mMessage = rErr.what();
            break;
        }
        // The frame's notes go under the picture; if their number changed, fit
        // the picture to what is left and draw once more.
        const int wanted = static_cast<int>(frame.mNotes.size());
        if (wanted == mNoteRows || pass == 1)
            break;
        mNoteRows = std::min(wanted, std::max(0, rows - 4));
        pic_rows = std::max(1, rows - 1 - mNoteRows);
        Full = true;
        out = "\x1b[0m\x1b[2J";
    }

    if (cells) {
        // After a clear the screen is blank, so a full repaint only writes what is
        // not.
        TextGrid blank;
        if (Full) {
            blank.mCols = grid.mCols;
            blank.mRows = grid.mRows;
            blank.mCells.assign(grid.mCells.size(), TextCell{});
        }
        out += encode_cells_update(Full ? blank : mPrevious, grid, mOptions.mText.mDepth, 1, 1);
        mPrevious = grid;
        mReport.mScreen = grid;
    } else {
        out += "\x1b[H";
        if (mOptions.mText.mEncoding == TextEncoding::Kitty)
            out += "\x1b_Ga=d\x1b\\";  // the previous image goes before the new one
        out += graphics;
    }
    if (mHelp) {
        const int lines = static_cast<int>(sizeof(kHelp) / sizeof(kHelp[0]));
        for (int i = 0; i < lines && i < pic_rows; ++i)
            out += at_row(1 + i) + "\x1b[7m" + fit_columns(kHelp[i], std::min(cols, 58)) + "\x1b[0m";
        mHelpWasShown = true;
    } else if (mHelpWasShown) {
        mHelpWasShown = false;
        mPrevious = TextGrid{};  // the help covered cells the diff thinks are intact
        if (cells) {
            TextGrid blank;
            blank.mCols = grid.mCols;
            blank.mRows = grid.mRows;
            blank.mCells.assign(grid.mCells.size(), TextCell{});
            out = "\x1b[0m\x1b[2J" + encode_cells_update(blank, grid, mOptions.mText.mDepth);
            mPrevious = grid;
        }
    }

    // Notes, then the status line, each only when it changed.
    std::vector<std::string> lines;
    for (int i = 0; i < mNoteRows && i < static_cast<int>(frame.mNotes.size()); ++i)
        lines.push_back(fit_columns(frame.mNotes[static_cast<std::size_t>(i)], cols));
    std::string status = " " + (mOptions.mTitle.empty() ? std::string("meshio++") : mOptions.mTitle);
    status += "  az " + number(mCamera.mAzimuth, 0) + "  el " + number(mCamera.mElevation, 0) +
              "  zoom " + number(mCamera.mZoom, 2) + "  " +
              (mCamera.mProjection == RenderProjection::Perspective ? "persp" : "ortho");
    if (!mMessage.empty())
        status += "  | " + mMessage;
    status += "  | ? help  q quit";
    mReport.mStatus = status;
    lines.push_back("\x1b[7m" + fit_columns(status, cols) + "\x1b[0m");
    const bool relayout = Full;
    for (std::size_t i = 0; i < lines.size(); ++i) {
        const bool is_status = i + 1 == lines.size();
        if (!relayout && i < mPreviousLines.size() && mPreviousLines[i] == lines[i] && !is_status)
            continue;
        const int row = is_status ? rows : pic_rows + 1 + static_cast<int>(i);
        out += at_row(row) + lines[i] + "\x1b[0m";
    }
    mPreviousLines = lines;
    rIo.Write(out);
    ++mReport.mFrames;
}

TuiReport TuiSession::Run(TuiIo& rIo) {
    TerminalSize size;
    if (!rIo.Size(size)) {
        size.mCols = 80;
        size.mRows = 24;
    }
    InputParser parser;
    bool dirty = true;
    bool full = true;
    while (true) {
        if (const int sig = rIo.TerminationSignal(); sig != 0) {
            mReport.mExit = 128 + sig;
            break;
        }
        if (rIo.TakeResize()) {
            rIo.Size(size);
            dirty = full = true;
        }
        if (dirty) {
            Draw(rIo, size, full);
            dirty = full = false;
        }
        std::string bytes;
        const bool alive = rIo.Read(bytes, mOptions.mPollMs);
        parser.Feed(bytes);
        for (const InputEvent& event : parser.Take(bytes.empty())) {
            if (Handle(event, size))
                dirty = true;
            if (mQuit)
                break;
        }
        if (mQuit || !alive)
            break;
    }
    // Keys that arrived together with the one that ended the session (a recorded
    // script, a pasted burst) still count: show where they left the camera.
    if (dirty && mReport.mExit == 0)
        Draw(rIo, size, full);
    mReport.mFinal = mCamera;
    return mReport;
}

TuiReport run_on_terminal(const Mesh& rMesh, TuiOptions Options) {
    TuiReport report;
    try {
        // Prepare first: a bad option or a mesh that cannot be drawn is reported
        // on a normal screen, before anything about the terminal changes.
        TuiSession session(rMesh, std::move(Options));
        RawTerminal terminal;
        std::string error;
        if (!terminal.Enter(error)) {
            report.mExit = 1;
            report.mError = error;
            return report;
        }
        TerminalIo io(terminal);
        report = session.Run(io);
        terminal.Leave();
    } catch (const std::exception& rErr) {
        report.mExit = 1;
        report.mError = rErr.what();
    }
    return report;
}

}  // namespace meshioplusplus::cli::tui
