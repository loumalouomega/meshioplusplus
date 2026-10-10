// SPDX-License-Identifier: MIT
/// @file terminal.cpp
/// @brief Implementation of terminal.hpp.

#include "terminal.hpp"

#include <cerrno>
#include <atomic>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <io.h>
#include <windows.h>
#ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
#endif
#ifndef ENABLE_VIRTUAL_TERMINAL_INPUT
#define ENABLE_VIRTUAL_TERMINAL_INPUT 0x0200
#endif
#else
#include <poll.h>
#include <sys/ioctl.h>
#include <termios.h>
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

// ---------------------------------------------------------------------------
// The interactive session
// ---------------------------------------------------------------------------

namespace {

// Alternate screen, hidden cursor, button-and-drag mouse reporting in SGR
// encoding, bracketed paste -- and their exact inverses, which also reset the
// colours so a killed viewer cannot leave the shell tinted.
const char kEnterSequence[] = "\x1b[?1049h\x1b[?25l\x1b[?1000h\x1b[?1002h\x1b[?1006h\x1b[?2004h";
const char kLeaveSequence[] = "\x1b[?2004l\x1b[?1006l\x1b[?1002l\x1b[?1000l\x1b[0m\x1b[?25h\x1b[?1049l";

volatile std::sig_atomic_t gResize = 0;
volatile std::sig_atomic_t gSignal = 0;
RawTerminal* gActive = nullptr;
bool gAtexitRegistered = false;
std::atomic<long> gChildPid{0};

// Async-signal-safe: kill(2) only (a lock-free atomic load is too).
void kill_child() {
#ifndef _WIN32
    const long pid = gChildPid.load();
    if (pid > 0)
        kill(static_cast<pid_t>(pid), SIGKILL);
#endif
}

void leave_at_exit() {
    kill_child();
    if (gActive != nullptr)
        gActive->Leave();
}

}  // namespace

#ifdef _WIN32

struct RawTerminal::State {
    HANDLE mIn = INVALID_HANDLE_VALUE;
    HANDLE mOut = INVALID_HANDLE_VALUE;
    DWORD mInMode = 0;
    DWORD mOutMode = 0;
    UINT mInCp = 0;
    UINT mOutCp = 0;
    TerminalSize mSize;
};

namespace {

BOOL WINAPI on_console_event(DWORD Type) {
    if (Type == CTRL_C_EVENT || Type == CTRL_BREAK_EVENT || Type == CTRL_CLOSE_EVENT ||
        Type == CTRL_LOGOFF_EVENT || Type == CTRL_SHUTDOWN_EVENT) {
        gSignal = SIGTERM;
        return TRUE;
    }
    return FALSE;
}

}  // namespace

bool RawTerminal::Enter(std::string& rError) {
    if (mActive) {
        rError = "the terminal session is already active";
        return false;
    }
    if (gActive != nullptr) {
        rError = "another terminal session is active";
        return false;
    }
    mpState = new State;
    State& st = *mpState;
    st.mIn = GetStdHandle(STD_INPUT_HANDLE);
    st.mOut = GetStdHandle(STD_OUTPUT_HANDLE);
    if (st.mIn == INVALID_HANDLE_VALUE || st.mOut == INVALID_HANDLE_VALUE ||
        !GetConsoleMode(st.mIn, &st.mInMode) || !GetConsoleMode(st.mOut, &st.mOutMode)) {
        rError =
            "the interactive viewer needs a console on standard input and output; use "
            "`snapshot` for one frame to a file or a pipe";
        delete mpState;
        mpState = nullptr;
        return false;
    }
    st.mInCp = GetConsoleCP();
    st.mOutCp = GetConsoleOutputCP();
    DWORD in_mode = st.mInMode;
    in_mode &= ~static_cast<DWORD>(ENABLE_ECHO_INPUT | ENABLE_LINE_INPUT | ENABLE_PROCESSED_INPUT |
                                   ENABLE_QUICK_EDIT_MODE | ENABLE_MOUSE_INPUT);
    in_mode |= ENABLE_EXTENDED_FLAGS | ENABLE_VIRTUAL_TERMINAL_INPUT;
    if (!SetConsoleMode(st.mIn, in_mode)) {
        rError =
            "this console does not take virtual-terminal input (SetConsoleMode refused "
            "ENABLE_VIRTUAL_TERMINAL_INPUT, as an older Windows does); use `snapshot` for one "
            "frame";
        delete mpState;
        mpState = nullptr;
        return false;
    }
    SetConsoleCP(CP_UTF8);
    SetConsoleOutputCP(CP_UTF8);
    if (!(st.mOutMode & ENABLE_VIRTUAL_TERMINAL_PROCESSING))
        SetConsoleMode(st.mOut, st.mOutMode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
    terminal_size(st.mSize);
    SetConsoleCtrlHandler(on_console_event, TRUE);
    gSignal = 0;
    gResize = 0;
    gActive = this;
    if (!gAtexitRegistered) {
        std::atexit(leave_at_exit);
        gAtexitRegistered = true;
    }
    mActive = true;
    Write(kEnterSequence);
    return true;
}

void RawTerminal::Leave() {
    if (!mActive)
        return;
    mActive = false;
    State& st = *mpState;
    Write(kLeaveSequence);
    SetConsoleCtrlHandler(on_console_event, FALSE);
    SetConsoleMode(st.mIn, st.mInMode);
    SetConsoleMode(st.mOut, st.mOutMode);
    SetConsoleCP(st.mInCp);
    SetConsoleOutputCP(st.mOutCp);
    gActive = nullptr;
    delete mpState;
    mpState = nullptr;
}

bool RawTerminal::Read(std::string& rBytes, int TimeoutMs) {
    if (!mActive)
        return false;
    State& st = *mpState;
    const DWORD wait = WaitForSingleObject(st.mIn, static_cast<DWORD>(TimeoutMs < 0 ? 0 : TimeoutMs));
    if (wait != WAIT_OBJECT_0)
        return true;
    // Anything that is not a key (a focus change, a buffer resize) wakes the wait
    // without bytes to read: drop it rather than block in ReadFile.
    INPUT_RECORD records[64];
    DWORD count = 0;
    if (!PeekConsoleInputW(st.mIn, records, 64, &count) || count == 0)
        return true;
    bool has_key = false;
    for (DWORD i = 0; i < count; ++i)
        if (records[i].EventType == KEY_EVENT && records[i].Event.KeyEvent.bKeyDown)
            has_key = true;
    if (!has_key) {
        ReadConsoleInputW(st.mIn, records, count, &count);
        return true;
    }
    char buffer[4096];
    DWORD got = 0;
    if (!ReadFile(st.mIn, buffer, sizeof(buffer), &got, nullptr))
        return false;
    rBytes.append(buffer, got);
    return true;
}

bool RawTerminal::Write(const std::string& rBytes) {
    HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
    std::size_t done = 0;
    while (done < rBytes.size()) {
        DWORD wrote = 0;
        if (!WriteFile(out, rBytes.data() + done, static_cast<DWORD>(rBytes.size() - done), &wrote,
                       nullptr))
            return false;
        done += wrote;
    }
    return true;
}

bool RawTerminal::TakeResize() {
    if (!mActive)
        return false;
    TerminalSize now;
    if (!terminal_size(now))
        return false;
    TerminalSize& was = mpState->mSize;
    const bool changed = now.mCols != was.mCols || now.mRows != was.mRows;
    was = now;
    return changed;
}

int RawTerminal::TerminationSignal() const {
    return static_cast<int>(gSignal);
}

#else  // POSIX

struct RawTerminal::State {
    struct termios mSaved;
    struct sigaction mOldInt;
    struct sigaction mOldTerm;
    struct sigaction mOldHup;
    struct sigaction mOldWinch;
    struct sigaction mOldPipe;
};

namespace {

struct termios gSavedTermios;
bool gHaveSaved = false;

// Async-signal-safe: write(2) and tcsetattr(3) only.
void restore_from_signal() {
    kill_child();
    if (!gHaveSaved)
        return;
    const ssize_t ignored = write(STDOUT_FILENO, kLeaveSequence, sizeof(kLeaveSequence) - 1);
    (void)ignored;
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &gSavedTermios);
}

void on_signal(int Sig) {
    if (Sig == SIGWINCH) {
        gResize = 1;
        return;
    }
    // The loop polls this flag every few tens of milliseconds and unwinds
    // normally. A second signal means it is not listening (a long frame, a
    // hang): restore the terminal here and go.
    if (gSignal != 0) {
        restore_from_signal();
        _exit(128 + Sig);
    }
    gSignal = Sig;
}

}  // namespace

bool RawTerminal::Enter(std::string& rError) {
    if (mActive) {
        rError = "the terminal session is already active";
        return false;
    }
    if (gActive != nullptr) {
        rError = "another terminal session is active";
        return false;
    }
    if (!stdin_is_terminal() || !stdout_is_terminal()) {
        rError =
            "the interactive viewer needs a terminal on standard input and output; use "
            "`snapshot` for one frame to a file or a pipe";
        return false;
    }
    mpState = new State;
    State& st = *mpState;
    if (tcgetattr(STDIN_FILENO, &st.mSaved) != 0) {
        rError = std::string("cannot read the terminal settings: ") + std::strerror(errno);
        delete mpState;
        mpState = nullptr;
        return false;
    }
    struct termios raw = st.mSaved;
    raw.c_iflag &= ~static_cast<tcflag_t>(BRKINT | ICRNL | INPCK | ISTRIP | IXON);
    raw.c_lflag &= ~static_cast<tcflag_t>(ECHO | ICANON | IEXTEN | ISIG);
    raw.c_cc[VMIN] = 0;
    raw.c_cc[VTIME] = 0;
    if (tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw) != 0) {
        rError = std::string("cannot set raw mode: ") + std::strerror(errno);
        delete mpState;
        mpState = nullptr;
        return false;
    }
    gSavedTermios = st.mSaved;
    gHaveSaved = true;
    gSignal = 0;
    gResize = 0;
    struct sigaction action;
    std::memset(&action, 0, sizeof(action));
    action.sa_handler = on_signal;
    sigemptyset(&action.sa_mask);
    sigaction(SIGINT, &action, &st.mOldInt);
    sigaction(SIGTERM, &action, &st.mOldTerm);
    sigaction(SIGHUP, &action, &st.mOldHup);
    sigaction(SIGWINCH, &action, &st.mOldWinch);
    struct sigaction ignore;
    std::memset(&ignore, 0, sizeof(ignore));
    ignore.sa_handler = SIG_IGN;
    sigemptyset(&ignore.sa_mask);
    sigaction(SIGPIPE, &ignore, &st.mOldPipe);
    gActive = this;
    if (!gAtexitRegistered) {
        std::atexit(leave_at_exit);
        gAtexitRegistered = true;
    }
    mActive = true;
    Write(kEnterSequence);
    return true;
}

void RawTerminal::Leave() {
    if (!mActive)
        return;
    mActive = false;
    State& st = *mpState;
    Write(kLeaveSequence);
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &st.mSaved);
    gHaveSaved = false;
    sigaction(SIGINT, &st.mOldInt, nullptr);
    sigaction(SIGTERM, &st.mOldTerm, nullptr);
    sigaction(SIGHUP, &st.mOldHup, nullptr);
    sigaction(SIGWINCH, &st.mOldWinch, nullptr);
    sigaction(SIGPIPE, &st.mOldPipe, nullptr);
    gActive = nullptr;
    delete mpState;
    mpState = nullptr;
}

bool RawTerminal::Read(std::string& rBytes, int TimeoutMs) {
    if (!mActive)
        return false;
    struct pollfd fd;
    fd.fd = STDIN_FILENO;
    fd.events = POLLIN;
    fd.revents = 0;
    const int ready = poll(&fd, 1, TimeoutMs < 0 ? 0 : TimeoutMs);
    if (ready < 0)
        return errno == EINTR;
    if (ready == 0)
        return true;
    if ((fd.revents & (POLLERR | POLLNVAL)) != 0)
        return false;
    char buffer[4096];
    const ssize_t got = read(STDIN_FILENO, buffer, sizeof(buffer));
    if (got < 0)
        return errno == EINTR || errno == EAGAIN;
    if (got == 0)
        return (fd.revents & POLLHUP) == 0;  // a hung-up terminal ends the session
    rBytes.append(buffer, static_cast<std::size_t>(got));
    return true;
}

bool RawTerminal::Write(const std::string& rBytes) {
    std::size_t done = 0;
    while (done < rBytes.size()) {
        const ssize_t wrote = write(STDOUT_FILENO, rBytes.data() + done, rBytes.size() - done);
        if (wrote < 0) {
            if (errno == EINTR)
                continue;
            if (errno == EAGAIN) {
                struct pollfd fd;
                fd.fd = STDOUT_FILENO;
                fd.events = POLLOUT;
                fd.revents = 0;
                poll(&fd, 1, 100);
                continue;
            }
            return false;
        }
        done += static_cast<std::size_t>(wrote);
    }
    return true;
}

bool RawTerminal::TakeResize() {
    const bool changed = gResize != 0;
    gResize = 0;
    return changed;
}

int RawTerminal::TerminationSignal() const {
    return static_cast<int>(gSignal);
}

#endif

RawTerminal::~RawTerminal() {
    Leave();
}

void terminal_track_child(long Pid) { gChildPid.store(Pid); }

}  // namespace meshioplusplus::cli
