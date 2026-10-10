// SPDX-License-Identifier: MIT
/// @file session_file.hpp
/// @brief The viewer's session file: a view someone reached, as strict JSON.
#ifndef MESHIOPLUSPLUS_CLI_TUI_SESSION_FILE_HPP
#define MESHIOPLUSPLUS_CLI_TUI_SESSION_FILE_HPP

// A session is the render flags the viewer is showing (the camera, the field and
// its range and scale, the cut-away planes, ...; the same tokens `snapshot` takes)
// and the series step. It is written with sorted keys and a version, and read back
// strictly: an unknown key, a wrong type or another version is an error naming it,
// so a file that someone attaches to an issue means one thing. The reader is a
// small hand-written one for exactly this shape, so the viewer does not depend on
// the optional JSON library.

#include <string>
#include <vector>

namespace meshioplusplus::cli::tui {

inline constexpr int kSessionVersion = 1;

struct SessionFile {
    std::vector<std::string> mFlags;  ///< argv-style tokens, `--name=value`
    long long mStep = 0;
};

/// The JSON text (a trailing newline included).
std::string session_to_json(const SessionFile& rSession);

/// Parse a session. Throws `std::invalid_argument` naming what is wrong.
SessionFile session_from_json(const std::string& rText);

}  // namespace meshioplusplus::cli::tui

#endif  // MESHIOPLUSPLUS_CLI_TUI_SESSION_FILE_HPP
