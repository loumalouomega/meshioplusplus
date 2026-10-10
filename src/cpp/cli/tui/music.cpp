// SPDX-License-Identifier: MIT
// The viewer's optional soundtrack. See music.hpp.

#include "music.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "../terminal.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <csignal>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#ifdef __linux__
#include <sys/prctl.h>
#endif
#endif

namespace meshioplusplus::cli::tui {

namespace {

bool env_set(const char* pName) {
    const char* p = std::getenv(pName);
    return p != nullptr && *p != '\0';
}

std::vector<std::string> split_path(const std::string& rPath) {
#ifdef _WIN32
    const char sep = ';';
#else
    const char sep = ':';
#endif
    std::vector<std::string> out;
    std::size_t start = 0;
    while (start <= rPath.size()) {
        std::size_t end = rPath.find(sep, start);
        if (end == std::string::npos)
            end = rPath.size();
        if (end > start)
            out.push_back(rPath.substr(start, end - start));
        start = end + 1;
    }
    return out;
}

bool is_executable(const std::filesystem::path& rPath) {
    std::error_code ec;
    if (!std::filesystem::is_regular_file(rPath, ec))
        return false;
#ifdef _WIN32
    return true;
#else
    return access(rPath.c_str(), X_OK) == 0;
#endif
}

}  // namespace

std::vector<MusicPlayerSpec> music_player_candidates() {
#ifdef _WIN32
    return {{"powershell",
             {"-NoProfile", "-NonInteractive", "-Command",
              "(New-Object System.Media.SoundPlayer '{file}').PlaySync()"}}};
#else
    return {{"afplay", {}},
            {"paplay", {}},
            {"pw-play", {}},
            {"aplay", {"-q"}},
            {"ffplay", {"-nodisp", "-autoexit", "-loglevel", "quiet"}}};
#endif
}

bool music_find_player(MusicPlayerSpec& rSpec, std::string& rPath) {
    const char* hook = std::getenv("MESHIOPLUSPLUS_MUSIC_PLAYERS");
    const char* path = hook != nullptr ? hook : std::getenv("PATH");
    if (path == nullptr)
        return false;
    const std::vector<std::string> dirs = split_path(path);
    for (const MusicPlayerSpec& spec : music_player_candidates())
        for (const std::string& dir : dirs) {
            std::filesystem::path candidate = std::filesystem::path(dir) / spec.mName;
#ifdef _WIN32
            candidate += ".exe";
#endif
            if (is_executable(candidate)) {
                rSpec = spec;
                rPath = candidate.string();
                return true;
            }
        }
    return false;
}

bool music_allowed(bool OverSsh, std::string& rReason) {
    if (env_set("CI")) {
        rReason = "CI is set, so nothing plays (there is nobody to hear it)";
        return false;
    }
    if (!OverSsh && (env_set("SSH_CONNECTION") || env_set("SSH_CLIENT") || env_set("SSH_TTY"))) {
        rReason =
            "this is an SSH session: the music would play on the machine the process runs on, "
            "not where the terminal is (--music-over-ssh plays it there anyway; --music-out "
            "FILE.wav writes it for your own player)";
        return false;
    }
    return true;
}

bool reduced_motion_requested() {
    const char* p = std::getenv("REDUCED_MOTION");
    if (p == nullptr)
        return false;
    const std::string v = p;
    return !(v.empty() || v == "0" || v == "false" || v == "FALSE" || v == "False");
}

void music_write_wav(const std::string& rPath, const meshioplusplus::detail::SynthOptions& rOptions) {
    const std::string wav = meshioplusplus::detail::synth_wav(
        meshioplusplus::detail::synth_synthwave(rOptions), rOptions.mSampleRate);
    auto out = meshioplusplus::detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!out)
        throw std::runtime_error("music: cannot write '" + rPath + "'");
    out.write(wav.data(), static_cast<std::streamsize>(wav.size()));
}

MusicPlayer::MusicPlayer(meshioplusplus::detail::SynthOptions Options, double Volume)
    : mOptions(Options), mVolume(std::clamp(Volume, 0.05, 0.9)) {
    mOptions.mGain = mVolume;
}

MusicPlayer::~MusicPlayer() { Stop(); }

void MusicPlayer::Regenerate() {
    mOptions.mGain = mVolume;
    music_write_wav(mWavPath, mOptions);
}

bool MusicPlayer::Start(std::string& rMessage) {
    std::string tried;
    for (const MusicPlayerSpec& spec : music_player_candidates())
        tried += (tried.empty() ? "" : ", ") + spec.mName;
    if (!music_find_player(mSpec, mPlayerPath)) {
        rMessage = "music: no player found (looked for " + tried +
                   "); the viewer goes on silently. --music-out FILE.wav writes the loop for "
                   "your own player";
        return false;
    }
    std::error_code ec;
    mWavPath = (std::filesystem::temp_directory_path(ec) /
                ("meshioplusplus_music_" + std::to_string(static_cast<long>(
#ifdef _WIN32
                                                             GetCurrentProcessId()
#else
                                                             getpid()
#endif
                                                                 )) +
                 ".wav"))
                   .string();
    Regenerate();
    mStarted = true;
    if (!Launch()) {
        rMessage = "music: could not start " + mSpec.mName;
        return false;
    }
    rMessage = "music: " + mSpec.mName;
    return true;
}

#ifdef _WIN32

bool MusicPlayer::Launch() {
    std::string cmd = "\"" + mPlayerPath + "\"";
    for (std::string arg : mSpec.mArgs) {
        const auto at = arg.find("{file}");
        if (at != std::string::npos)
            arg.replace(at, 6, mWavPath);
        cmd += " \"" + arg + "\"";
    }
    STARTUPINFOA si{};
    si.cb = sizeof si;
    PROCESS_INFORMATION pi{};
    if (!CreateProcessA(nullptr, cmd.data(), nullptr, nullptr, FALSE, CREATE_NO_WINDOW, nullptr,
                        nullptr, &si, &pi))
        return false;
    CloseHandle(pi.hThread);
    mPid = static_cast<long>(pi.dwProcessId);
    mProcess = pi.hProcess;
    return true;
}

bool MusicPlayer::Reap() {
    if (mPid <= 0)
        return true;
    if (WaitForSingleObject(static_cast<HANDLE>(mProcess), 0) == WAIT_OBJECT_0) {
        CloseHandle(static_cast<HANDLE>(mProcess));
        mProcess = nullptr;
        mPid = 0;
        return true;
    }
    return false;
}

void MusicPlayer::Stop() {
    if (mPid > 0 && mProcess != nullptr) {
        TerminateProcess(static_cast<HANDLE>(mProcess), 0);
        WaitForSingleObject(static_cast<HANDLE>(mProcess), 1000);
        CloseHandle(static_cast<HANDLE>(mProcess));
    }
    mProcess = nullptr;
    mPid = 0;
    mStarted = false;
    if (!mWavPath.empty()) {
        std::error_code ec;
        std::filesystem::remove(mWavPath, ec);
    }
}

#else

bool MusicPlayer::Launch() {
    std::vector<std::string> argv_store = {mSpec.mName};
    bool placed = false;
    for (std::string arg : mSpec.mArgs) {
        const auto at = arg.find("{file}");
        if (at != std::string::npos) {
            arg.replace(at, 6, mWavPath);
            placed = true;
        }
        argv_store.push_back(arg);
    }
    if (!placed)
        argv_store.push_back(mWavPath);
    std::vector<char*> argv;
    for (std::string& s : argv_store)
        argv.push_back(s.data());
    argv.push_back(nullptr);
    const pid_t pid = fork();
    if (pid < 0)
        return false;
    if (pid == 0) {
#ifdef __linux__
        // Die with the viewer even if it is killed outright.
        prctl(PR_SET_PDEATHSIG, SIGKILL);
#endif
        const int null = open("/dev/null", O_RDWR);
        if (null >= 0) {
            dup2(null, STDIN_FILENO);
            dup2(null, STDOUT_FILENO);
            dup2(null, STDERR_FILENO);
        }
        execv(mPlayerPath.c_str(), argv.data());
        _exit(127);
    }
    mPid = static_cast<long>(pid);
    meshioplusplus::cli::terminal_track_child(mPid);
    return true;
}

bool MusicPlayer::Reap() {
    if (mPid <= 0)
        return true;
    int status = 0;
    const pid_t done = waitpid(static_cast<pid_t>(mPid), &status, WNOHANG);
    if (done == static_cast<pid_t>(mPid) || done < 0) {
        mPid = 0;
        meshioplusplus::cli::terminal_track_child(0);
        return true;
    }
    return false;
}

void MusicPlayer::Stop() {
    if (mPid > 0) {
        kill(static_cast<pid_t>(mPid), SIGKILL);
        int status = 0;
        waitpid(static_cast<pid_t>(mPid), &status, 0);
        mPid = 0;
        meshioplusplus::cli::terminal_track_child(0);
    }
    mStarted = false;
    if (!mWavPath.empty()) {
        std::error_code ec;
        std::filesystem::remove(mWavPath, ec);
    }
}

#endif

void MusicPlayer::Poll() {
    if (!mStarted || mMuted)
        return;
    if (Reap())
        Launch();  // the loop ended: play it again
}

void MusicPlayer::SetMuted(bool Muted) {
    if (Muted == mMuted || !mStarted)
        return;
    mMuted = Muted;
    if (mMuted) {
#ifndef _WIN32
        if (mPid > 0) {
            kill(static_cast<pid_t>(mPid), SIGKILL);
            int status = 0;
            waitpid(static_cast<pid_t>(mPid), &status, 0);
            mPid = 0;
            meshioplusplus::cli::terminal_track_child(0);
        }
#else
        if (mPid > 0 && mProcess != nullptr) {
            TerminateProcess(static_cast<HANDLE>(mProcess), 0);
            CloseHandle(static_cast<HANDLE>(mProcess));
            mProcess = nullptr;
            mPid = 0;
        }
#endif
    } else {
        Launch();
    }
}

void MusicPlayer::AdjustVolume(int Steps) {
    mVolume = std::clamp(std::round((mVolume + 0.05 * Steps) * 100.0) / 100.0, 0.05, 0.9);
    if (!mStarted)
        return;
    const bool was_muted = mMuted;
#ifndef _WIN32
    if (mPid > 0) {
        kill(static_cast<pid_t>(mPid), SIGKILL);
        int status = 0;
        waitpid(static_cast<pid_t>(mPid), &status, 0);
        mPid = 0;
        meshioplusplus::cli::terminal_track_child(0);
    }
#else
    if (mPid > 0 && mProcess != nullptr) {
        TerminateProcess(static_cast<HANDLE>(mProcess), 0);
        CloseHandle(static_cast<HANDLE>(mProcess));
        mProcess = nullptr;
        mPid = 0;
    }
#endif
    Regenerate();
    if (!was_muted)
        Launch();
}

std::string MusicPlayer::Status() const {
    if (!mStarted)
        return "music: off";
    if (mMuted)
        return "music: muted";
    return "music: " + std::to_string(static_cast<int>(std::lround(mVolume * 100.0))) + "%";
}

}  // namespace meshioplusplus::cli::tui
