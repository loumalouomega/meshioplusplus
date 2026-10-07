// SPDX-License-Identifier: MIT
/// @file terminal.cpp
/// @brief Implementation of terminal.hpp.

#include "terminal.hpp"

#include <cstdio>
#include <cstdlib>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <io.h>
#include <windows.h>
#ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
#endif
#else
#include <sys/ioctl.h>
#include <unistd.h>
#endif

namespace meshioplusplus::cli {

bool stdout_is_terminal() {
#ifdef _WIN32
    return _isatty(_fileno(stdout)) != 0;
#else
    return isatty(STDOUT_FILENO) != 0;
#endif
}

bool stdin_is_terminal() {
#ifdef _WIN32
    return _isatty(_fileno(stdin)) != 0;
#else
    return isatty(STDIN_FILENO) != 0;
#endif
}

bool terminal_size(TerminalSize& rSize) {
#ifdef _WIN32
    CONSOLE_SCREEN_BUFFER_INFO info;
    if (!GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info))
        return false;
    rSize.mCols = info.srWindow.Right - info.srWindow.Left + 1;
    rSize.mRows = info.srWindow.Bottom - info.srWindow.Top + 1;
    rSize.mPixelWidth = 0;
    rSize.mPixelHeight = 0;
    return rSize.mCols > 0 && rSize.mRows > 0;
#else
    struct winsize ws{};
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) != 0)
        return false;
    rSize.mCols = ws.ws_col;
    rSize.mRows = ws.ws_row;
    rSize.mPixelWidth = ws.ws_xpixel;
    rSize.mPixelHeight = ws.ws_ypixel;
    return rSize.mCols > 0 && rSize.mRows > 0;
#endif
}

bool enable_terminal_output() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    if (out == INVALID_HANDLE_VALUE || !GetConsoleMode(out, &mode))
        return true;  // not a console (a pipe or a file): nothing to switch
    if (mode & ENABLE_VIRTUAL_TERMINAL_PROCESSING)
        return true;
    return SetConsoleMode(out, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING) != 0;
#else
    return true;
#endif
}

ColorDepth environment_color_depth() {
    const char* no_color = std::getenv("NO_COLOR");
    const char* color_term = std::getenv("COLORTERM");
    const char* term = std::getenv("TERM");
    ColorDepth depth = detect_color_depth(no_color, color_term, term);
#ifdef _WIN32
    // Windows Terminal speaks 24-bit colour but sets no COLORTERM.
    if (depth == ColorDepth::Ansi16 && std::getenv("WT_SESSION") != nullptr)
        depth = ColorDepth::TrueColor;
#endif
    return depth;
}

bool inside_tmux() {
    const char* tmux = std::getenv("TMUX");
    return tmux != nullptr && tmux[0] != '\0';
}

}  // namespace meshioplusplus::cli
