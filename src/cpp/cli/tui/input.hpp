// SPDX-License-Identifier: MIT
/// @file input.hpp
/// @brief A pure byte-to-event parser for the interactive viewer's input.
#ifndef MESHIOPLUSPLUS_CLI_TUI_INPUT_HPP
#define MESHIOPLUSPLUS_CLI_TUI_INPUT_HPP

// What a terminal sends in raw mode, as events: printable characters (UTF-8),
// control keys, arrows and the other CSI/SS3 keys with their modifiers, SGR
// mouse reports (`CSI ? 1006`) and bracketed paste. No tty, no clock, no
// state beyond the unparsed tail, so it is unit-tested with byte strings.

#include <cstdint>
#include <string>
#include <vector>

namespace meshioplusplus::cli::tui {

enum class EventKind : std::uint8_t { Key, Mouse, Paste };

enum class KeyCode : std::uint8_t {
    Char,  ///< a printable character, in `mChar`
    Up,
    Down,
    Left,
    Right,
    Home,
    End,
    PageUp,
    PageDown,
    Insert,
    Delete,
    Enter,
    Tab,
    Backspace,
    Escape,
    F1,
    F2,
    F3,
    F4,
    F5,
    F6,
    F7,
    F8,
    F9,
    F10,
    F11,
    F12,
};

enum class MouseAction : std::uint8_t { Press, Release, Drag, Move, WheelUp, WheelDown };

struct InputEvent {
    EventKind mKind = EventKind::Key;
    // Key
    KeyCode mKey = KeyCode::Char;
    std::uint32_t mChar = 0;  ///< the code point of a `Char` (the letter for Ctrl-letter)
    bool mCtrl = false;
    bool mAlt = false;
    bool mShift = false;
    // Mouse (cells, 1-based, as the terminal reports them)
    MouseAction mAction = MouseAction::Press;
    int mButton = 0;  ///< 0 left, 1 middle, 2 right
    int mX = 0;
    int mY = 0;
    // Paste
    std::string mText;

    bool IsChar(std::uint32_t Ch) const {
        return mKind == EventKind::Key && mKey == KeyCode::Char && mChar == Ch && !mCtrl &&
               !mAlt;
    }
    bool IsCtrl(std::uint32_t Ch) const {
        return mKind == EventKind::Key && mKey == KeyCode::Char && mChar == Ch && mCtrl;
    }
};

/// Accumulates bytes and yields the events they complete. An escape sequence
/// split across reads waits for its tail; `Flush` (call it when a read timed
/// out with nothing new) turns a lone `ESC` into the Escape key and drops an
/// unfinished sequence.
class InputParser {
public:
    void Feed(const char* pBytes, std::size_t Count);
    void Feed(const std::string& rBytes) { Feed(rBytes.data(), rBytes.size()); }
    /// The events parsed so far (and removed from the parser).
    std::vector<InputEvent> Take(bool Flush = false);
    /// Bytes still waiting for the rest of a sequence.
    std::size_t Pending() const { return mBuffer.size(); }

private:
    std::string mBuffer;
    bool mInPaste = false;
    std::string mPaste;
};

}  // namespace meshioplusplus::cli::tui

#endif  // MESHIOPLUSPLUS_CLI_TUI_INPUT_HPP
