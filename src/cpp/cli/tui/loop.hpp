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
//
// On top of the camera: probing the cell under the pointer, a `:` command line
// that speaks the render flags of `snapshot`, sessions, two meshes side by side,
// cut-away planes, and a time series that can be stepped, played and followed.

#include <array>
#include <chrono>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "../render_args.hpp"
#include "../terminal.hpp"
#include "input.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/render.hpp"
#include "music.hpp"
#include "probe.hpp"
#include "series.hpp"

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
    /// A monotonic clock in milliseconds, for playing and following. A scripted
    /// input supplies its own so a replay does not depend on the wall clock.
    virtual long long NowMs() {
        using namespace std::chrono;
        return duration_cast<milliseconds>(steady_clock::now().time_since_epoch()).count();
    }
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
/// keeps what was written. After the script the input ends, which quits. Its
/// clock advances by the timeout of each read, so playing and following are
/// deterministic.
class ScriptedIo : public TuiIo {
public:
    ScriptedIo(std::string Script, int Cols, int Rows, int PixelWidth = 0, int PixelHeight = 0)
        : mScript(std::move(Script)), mCols(Cols), mRows(Rows), mPixelWidth(PixelWidth),
          mPixelHeight(PixelHeight) {}
    bool Read(std::string& rBytes, int TimeoutMs) override {
        mClock += TimeoutMs;
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
    long long NowMs() override { return mClock; }
    const std::string& Output() const { return mOutput; }

private:
    std::string mScript;
    std::size_t mAt = 0;
    int mCols;
    int mRows;
    int mPixelWidth;
    int mPixelHeight;
    long long mClock = 0;
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

    /// A second mesh drawn beside the first, under one camera (cell encodings
    /// only). It must outlive the session.
    const Mesh* mpCompare = nullptr;
    std::string mCompareTitle;
    /// With a second mesh: one colour range over both (the default), or each its
    /// own.
    bool mSharedRange = true;
    /// With a second mesh: draw |B - A| of the `color_by` point array, node by
    /// node, in place of B (the meshes must share their nodes).
    bool mDiff = false;

    /// A time series to step through, in place of the mesh the session is built
    /// with (that mesh is then ignored). Start on the newest step with `mFollowMs`.
    std::shared_ptr<TuiSeries> mpSeries;
    /// Poll the series for a new step every this many ms (0: no).
    int mFollowMs = 0;
    /// A new step is read only once its file has stopped changing for this long.
    int mFollowSettleMs = 300;
    /// Steps per second when playing.
    double mPlayFps = 4.0;

    /// A session file: read at the start when it exists, written when the viewer
    /// ends (and by `:session save`).
    std::string mSessionPath;

    /// The optional soundtrack (not owned): polled every turn so the loop
    /// restarts, `m` mutes it and `<` `>` change its volume. Null for none.
    MusicControl* mpMusic = nullptr;
    /// With a theme: the grid's colour and scroll step on every beat of this
    /// tempo, counted from the clock the loop is given (never the wall clock of
    /// a recording). 0 for a still theme. Beats stay under three a second.
    double mPulseTempo = 0.0;
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
    /// The step the series ended on (0 for a single mesh).
    std::size_t mStep = 0;
    /// The text of the last probe and pins, line by line.
    std::vector<std::string> mProbe;
};

/// One viewing session over a mesh (or a series, or two meshes).
class TuiSession {
public:
    /// Prepares the mesh; throws what `prepare_render` throws, before any
    /// terminal state is touched. With `Options.mpSeries` the series is read
    /// instead.
    TuiSession(const Mesh& rMesh, TuiOptions Options);

    /// Run until quit, end of input or a signal.
    TuiReport Run(TuiIo& rIo);

    /// Apply one input event (exposed for tests); returns true when the screen
    /// needs redrawing.
    bool Handle(const InputEvent& rEvent, const TerminalSize& rSize);
    bool Quit() const { return mQuit; }
    const RenderOptions& Camera() const { return mCamera; }

    /// Run one `:` command line (exposed for tests); returns the message to
    /// show, empty for none. Never throws: an error is the message.
    std::string Execute(const std::string& rLine);
    /// The series step shown now.
    std::size_t Step() const { return mStep; }
    /// The status line as of the last frame.
    const std::string& Status() const { return mReport.mStatus; }
    const std::vector<std::string>& ProbeLines() const { return mProbeLines; }

private:
    // One drawn mesh: its prepared scene and what it last showed.
    struct Pane {
        const Mesh* mpMesh = nullptr;
        std::unique_ptr<Mesh> mOwned;  ///< the diff clone, when there is one
        RenderScene mScene;
        std::string mTitle;
        Frame mFrame;
        TextGrid mGrid;
        int mCols = 0;       ///< cells wide, as drawn
        int mColOffset = 0;  ///< columns left of it on the screen
        int mGx = 1;         ///< pixels per cell, horizontally and vertically
        int mGy = 2;
    };

    void Draw(TuiIo& rIo, const TerminalSize& rSize, bool Full);
    void PrepareScenes();
    void RollBackPreparation();
    bool LoadStep(std::size_t Index);
    bool StepBy(long long Delta);
    void Tick();
    bool HandlePromptKey(const InputEvent& rEvent);
    void ApplyFlags(RenderFlags& rFlags, bool Prepare = true);
    bool ResetCamera();
    bool ToggleCutaway(int Axis);
    bool MoveCutaway(double Direction);
    void SyncCutaways();
    void ProbeAt(int X, int Y);
    void RebuildProbeLines();
    void SaveSession(const std::string& rPath);
    void LoadSession(const std::string& rPath);
    std::string WriteSnapshot(const std::string& rPath);

    const Mesh* mpMeshA = nullptr;  ///< the mesh shown (the series' current step)
    Mesh mOwnedA;                   ///< where a series step lives
    TuiOptions mOptions;
    TuiIo* mpIo = nullptr;          ///< where a command may say it is preparing
    std::vector<Pane> mPanes;       ///< one, or two when comparing
    RenderOptions mInitial;
    RenderOptions mCamera;          ///< the live camera and every option (the flags' state)
    RenderOptions mPrepared;        ///< the options the scenes were last prepared with
    bool mQuit = false;
    bool mHelp = false;
    bool mHelpWasShown = false;
    bool mReprepare = false;
    bool mPreparedSmooth = false;
    bool mFull = false;
    int mNoteRows = 0;   ///< rows reserved under the picture for the frame's notes
    int mProbeRows = 0;  ///< rows reserved for the probe and the pins
    std::string mMessage;
    // Mouse drag in progress.
    bool mDragging = false;
    bool mDragMoved = false;
    int mDragButton = 0;
    int mLastX = 0;
    int mLastY = 0;
    // The `:` prompt.
    bool mPrompt = false;
    std::string mPromptText;
    // Probing and pins.
    std::vector<std::string> mProbeLines;
    Probe mProbe;
    std::vector<Probe> mPins;
    // Cut-away planes the keys manage: sign (0 none, +1, -1) and offset per axis.
    std::array<int, 3> mCutSign = {0, 0, 0};
    std::array<double, 3> mCutAt = {0.0, 0.0, 0.0};
    int mCutActive = -1;
    // Series.
    std::shared_ptr<TuiSeries> mpSeries;
    std::size_t mStep = 0;
    bool mPlaying = false;
    long long mNextPlayMs = 0;
    bool mRangeLocked = true;
    bool mHaveLockedRange = false;
    double mLockedMin = 0.0;
    double mLockedMax = 0.0;
    // Following a live run.
    std::string mKnownFingerprint;
    std::string mSettlingFingerprint;
    bool mSettling = false;
    long long mNextFollowMs = 0;
    long long mSettleUntilMs = 0;
    long long mNowMs = 0;
    long long mPulseStartMs = -1;
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
