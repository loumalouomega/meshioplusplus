// SPDX-License-Identifier: MIT
/// @file music.hpp
/// @brief The viewer's optional soundtrack: a generated loop played by an external player.
#ifndef MESHIOPLUSPLUS_CLI_TUI_MUSIC_HPP
#define MESHIOPLUSPLUS_CLI_TUI_MUSIC_HPP

// The loop is synthesized by `detail/synth.*` (a pure function of its options),
// written as a WAV file and played by whichever player the system has, started as
// a child process, restarted when it ends and killed on every way out the process
// controls (`terminal_track_child` covers the signal and `atexit` paths). No audio
// library is linked. Nothing here is on by default, and `music_allowed` refuses where
// the sound would come out of the wrong machine or where nobody is listening.

#include <string>
#include <vector>

#include "../../src/detail/synth.hpp"

namespace meshioplusplus::cli::tui {

/// What the loop asks of the soundtrack. A test passes a recording stand-in.
class MusicControl {
public:
    virtual ~MusicControl() = default;
    /// Called once per loop turn: restart the loop when the player has ended.
    virtual void Poll() = 0;
    virtual void SetMuted(bool Muted) = 0;
    virtual bool Muted() const = 0;
    /// Change the volume by `Steps` tenths of the range and restart the loop.
    virtual void AdjustVolume(int Steps) = 0;
    virtual double Volume() const = 0;
    /// A few words for the status line ("music 30%" / "music muted").
    virtual std::string Status() const = 0;
};

/// The players tried, in order, and what each needs on its command line.
struct MusicPlayerSpec {
    std::string mName;
    std::vector<std::string> mArgs;  ///< before the file
};
std::vector<MusicPlayerSpec> music_player_candidates();

/// The first candidate found on `PATH` (or in the `MESHIOPLUSPLUS_MUSIC_PLAYERS`
/// search directories, a test hook), with its full path; false when none.
bool music_find_player(MusicPlayerSpec& rSpec, std::string& rPath);

/// Whether sound may play here. False with the reason when `CI` is set or the
/// session is an SSH one (unless `OverSsh`): the music plays where the process
/// runs, not where the terminal is.
bool music_allowed(bool OverSsh, std::string& rReason);

/// Whether the environment asks for reduced motion (`REDUCED_MOTION` set to
/// anything but empty, `0` or `false`).
bool reduced_motion_requested();

/// A looping soundtrack on an external player.
class MusicPlayer : public MusicControl {
public:
    MusicPlayer(meshioplusplus::detail::SynthOptions Options, double Volume);
    ~MusicPlayer() override;
    MusicPlayer(const MusicPlayer&) = delete;
    MusicPlayer& operator=(const MusicPlayer&) = delete;

    /// Find a player, write the loop and start it. False, with a message naming
    /// what was tried, when there is nothing to play it with (the loop goes on
    /// silently).
    bool Start(std::string& rMessage);
    void Stop();

    void Poll() override;
    void SetMuted(bool Muted) override;
    bool Muted() const override { return mMuted; }
    void AdjustVolume(int Steps) override;
    double Volume() const override { return mVolume; }
    std::string Status() const override;
    bool Playing() const { return mPid > 0; }

private:
    bool Launch();
    bool Reap();
    void Regenerate();

    meshioplusplus::detail::SynthOptions mOptions;
    double mVolume;
    bool mMuted = false;
    bool mStarted = false;
    MusicPlayerSpec mSpec;
    std::string mPlayerPath;
    std::string mWavPath;
    long mPid = 0;
    void* mProcess = nullptr;  ///< the Windows process handle
};

/// Write the loop as a WAV file (`--music-out`): no player is needed.
void music_write_wav(const std::string& rPath, const meshioplusplus::detail::SynthOptions& rOptions);

}  // namespace meshioplusplus::cli::tui

#endif  // MESHIOPLUSPLUS_CLI_TUI_MUSIC_HPP
