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

}  // namespace meshioplusplus::cli

#endif  // MESHIOPLUSPLUS_CLI_TERMINAL_HPP
