// SPDX-License-Identifier: MIT
/// @file loop.hpp
/// @brief The interactive viewer's event loop around the software rasterizer.
#ifndef MESHIOPLUSPLUS_CLI_TUI_LOOP_HPP
#define MESHIOPLUSPLUS_CLI_TUI_LOOP_HPP

// `tui`: orbit, zoom and pan a mesh in a terminal. The loop owns the camera
// state and nothing else; the mesh is prepared once (`prepare_render`: the
// skin, the field, the normals), each input batch redraws from the new camera
// (`render_scene`), and only the cells that changed are written back. Input
// and output go through `TuiIo`, so the same loop runs on a real terminal
// (`RawTerminal`), on a recorded byte stream (`ScriptedIo`, which the tests and
// the documentation screenshots use) and inside the Python binding.

#include <cstddef>
#include <string>
#include <vector>

#include "../terminal.hpp"
#include "input.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/render.hpp"

namespace meshioplusplus::cli::tui {

/// Where the loop reads and writes.
class TuiIo {
public:
    virtual ~TuiIo() = default;
    /// Wait up to `TimeoutMs` for input and append it; false at end of input.
    virtual bool Read(std::string& rBytes, int TimeoutMs) = 0;
    /// Write bytes to the screen; false once it is gone.
    virtual bool Write(const std::string& rBytes) = 0;
    /// The screen size in cells (and pixels when known).
    virtual bool Size(TerminalSize& rSize) = 0;
    /// Whether the screen was resized since the last call.
    virtual bool TakeResize() = 0;
    /// The signal that asked the process to stop, or 0.
    virtual int TerminationSignal() = 0;
};

/// A `TuiIo` over a `RawTerminal`.
class TerminalIo : public TuiIo {
public:
    explicit TerminalIo(RawTerminal& rTerminal) : mrTerminal(rTerminal) {}
    bool Read(std::string& rBytes, int TimeoutMs) override {
        return mrTerminal.Read(rBytes, TimeoutMs);
    }
    bool Write(const std::string& rBytes) override { return mrTerminal.Write(rBytes); }
    bool Size(TerminalSize& rSize) override { return terminal_size(rSize); }
    bool TakeResize() override { return mrTerminal.TakeResize(); }
    int TerminationSignal() override { return mrTerminal.TerminationSignal(); }

private:
    RawTerminal& mrTerminal;
};

/// A `TuiIo` that plays a recorded byte stream on a screen of a fixed size and
/// keeps what was written. After the script the input ends, which quits.
class ScriptedIo : public TuiIo {
public:
    ScriptedIo(std::string Script, int Cols, int Rows, int PixelWidth = 0, int PixelHeight = 0)
        : mScript(std::move(Script)), mCols(Cols), mRows(Rows), mPixelWidth(PixelWidth),
          mPixelHeight(PixelHeight) {}
    bool Read(std::string& rBytes, int) override {
        if (mAt >= mScript.size())
            return false;
        rBytes.append(mScript, mAt, mScript.size() - mAt);
        mAt = mScript.size();
        return true;
    }
    bool Write(const std::string& rBytes) override {
        mOutput += rBytes;
        return true;
    }
    bool Size(TerminalSize& rSize) override {
        rSize.mCols = mCols;
        rSize.mRows = mRows;
        rSize.mPixelWidth = mPixelWidth;
        rSize.mPixelHeight = mPixelHeight;
        return true;
    }
    bool TakeResize() override { return false; }
    int TerminationSignal() override { return 0; }
    const std::string& Output() const { return mOutput; }

private:
    std::string mScript;
    std::size_t mAt = 0;
    int mCols;
    int mRows;
    int mPixelWidth;
    int mPixelHeight;
    std::string mOutput;
};

struct TuiOptions {
    /// The starting camera and everything about the picture that is not the
    /// frame size: field, edges, shading, colours. `mWidth`, `mHeight` and
    /// `mPixelAspect` are set from the screen on every frame.
    RenderOptions mRender;
    /// The encoding, colour depth and cell geometry; `mCols`/`mRows` are set
    /// from the screen on every frame.
    TextOptions mText;
    /// Shown at the start of the status line (a file name).
    std::string mTitle;
    /// Keep the model one size and place while the camera moves.
    bool mFixedFit = true;
    /// How often the loop looks at the signal and resize flags, in ms.
    int mPollMs = 50;
    /// Take the cell aspect and the cell's pixel size from the terminal when it
    /// reports its pixel size (`mText` keeps them otherwise).
    bool mCellAspectFromTerminal = true;
};

/// What a run leaves behind.
struct TuiReport {
    /// 0 for a normal quit or end of input, 128 + signal for a signal, 1 for a
    /// failure (see `mError`).
    int mExit = 0;
    std::string mError;
    std::size_t mFrames = 0;
    /// The last picture, as cells, and the status line under it.
    TextGrid mScreen;
    std::string mStatus;
    /// The camera and toggles when the loop ended.
    RenderOptions mFinal;
};

/// One viewing session over a mesh.
class TuiSession {
public:
    /// Prepares the mesh; throws what `prepare_render` throws, before any
    /// terminal state is touched.
    TuiSession(const Mesh& rMesh, TuiOptions Options);

    /// Run until quit, end of input or a signal.
    TuiReport Run(TuiIo& rIo);

    /// Apply one input event to the camera (exposed for tests); returns true
    /// when the screen needs redrawing.
    bool Handle(const InputEvent& rEvent, const TerminalSize& rSize);
    bool Quit() const { return mQuit; }
    const RenderOptions& Camera() const { return mCamera; }

private:
    void Draw(TuiIo& rIo, const TerminalSize& rSize, bool Full);
    void Reprepare();
    bool ResetCamera();

    const Mesh& mrMesh;
    TuiOptions mOptions;
    RenderScene mScene;
    RenderOptions mInitial;
    RenderOptions mCamera;  ///< the live camera and draw-time toggles
    bool mQuit = false;
    bool mHelp = false;
    bool mHelpWasShown = false;
    bool mReprepare = false;
    bool mPreparedSmooth = false;
    int mNoteRows = 0;  ///< rows reserved under the picture for the frame's notes
    std::string mMessage;
    // Mouse drag in progress.
    bool mDragging = false;
    int mDragButton = 0;
    int mLastX = 0;
    int mLastY = 0;
    TuiReport mReport;
    TextGrid mPrevious;
    std::vector<std::string> mPreviousLines;
};

/// Run on the real terminal: prepare the mesh, enter raw mode, run, and leave
/// it again on every path. A failure (not a terminal, a bad option) comes back
/// as `mExit` 1 with `mError`; nothing is thrown and nothing is left changed.
TuiReport run_on_terminal(const Mesh& rMesh, TuiOptions Options);

}  // namespace meshioplusplus::cli::tui

#endif  // MESHIOPLUSPLUS_CLI_TUI_LOOP_HPP
