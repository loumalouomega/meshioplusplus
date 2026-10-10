// SPDX-License-Identifier: MIT
/// @file terminal.hpp
/// @brief Terminal queries for the CLI: tty tests, size, Windows VT output.
#ifndef MESHIOPLUSPLUS_CLI_TERMINAL_HPP
#define MESHIOPLUSPLUS_CLI_TERMINAL_HPP

// The terminal facts the CLI needs, behind one small interface: whether a
// stream is a terminal, its size in cells (and pixels, when it says), and, on
// Windows, switching the console to virtual-terminal output. The library never
// touches a tty; this lives in the CLI layer (roadmap 7.1.4).

#include <string>

#include "meshioplusplus/operations/render.hpp"

namespace meshioplusplus::cli {

/// A terminal's size: columns and rows, and its pixel size when it reports
/// one (0 otherwise).
struct TerminalSize {
    int mCols = 0;
    int mRows = 0;
    int mPixelWidth = 0;
    int mPixelHeight = 0;
};

/// Whether standard output (or input) is a terminal.
bool stdout_is_terminal();
bool stdin_is_terminal();

/// The size of the terminal on standard output; false when it is not one or
/// the size is unknown.
bool terminal_size(TerminalSize& rSize);

/// Prepare standard output for escape sequences and UTF-8 glyphs. A no-op on
/// POSIX; on Windows it sets the console's output code page to UTF-8 and turns
/// on ENABLE_VIRTUAL_TERMINAL_PROCESSING (Windows 10 1511 and later). Returns
/// false when the console refused (an older Windows), so the caller can fall
/// back to plain text.
bool enable_terminal_output();

/// The colour depth the environment advertises (`detect_color_depth` over
/// NO_COLOR, COLORTERM and TERM), with Windows Terminal (WT_SESSION) taken as
/// 24-bit.
ColorDepth environment_color_depth();

/// Whether the process runs inside tmux (TMUX is set).
bool inside_tmux();

/// The interactive session: raw input, the alternate screen, a hidden cursor,
/// SGR mouse reporting and bracketed paste, all undone by `Leave`.
///
/// The terminal is restored on every way out the process controls: `Leave`,
/// the destructor, an exception unwinding through it, `SIGINT`, `SIGTERM` and
/// `SIGHUP` (a flag the loop sees within one poll interval, so it unwinds
/// normally) and `std::exit` (an `atexit` hook). `SIGKILL` cannot be caught.
/// Only one session can be active at a time.
class RawTerminal {
public:
    RawTerminal() = default;
    ~RawTerminal();
    RawTerminal(const RawTerminal&) = delete;
    RawTerminal& operator=(const RawTerminal&) = delete;

    /// Start the session. Fails, with a message naming `snapshot` as the
    /// alternative, when standard input or output is not a terminal.
    bool Enter(std::string& rError);
    /// End the session; safe to call more than once.
    void Leave();
    bool Active() const { return mActive; }

    /// Wait up to `TimeoutMs` for input and append what arrived to `rBytes`.
    /// Returns false at end of input, true otherwise (also after a timeout or a
    /// signal, with nothing appended).
    bool Read(std::string& rBytes, int TimeoutMs);
    /// Write every byte to standard output; false once the output is gone.
    bool Write(const std::string& rBytes);
    /// Whether the terminal was resized since the last call.
    bool TakeResize();
    /// The signal that asked the process to stop (SIGINT, SIGTERM, SIGHUP),
    /// or 0.
    int TerminationSignal() const;

private:
    bool mActive = false;
    struct State;
    State* mpState = nullptr;
};

}  // namespace meshioplusplus::cli

#endif  // MESHIOPLUSPLUS_CLI_TERMINAL_HPP
