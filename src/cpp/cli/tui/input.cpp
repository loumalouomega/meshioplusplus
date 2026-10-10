// SPDX-License-Identifier: MIT
/// @file input.cpp
/// @brief Implementation of input.hpp.

#include "input.hpp"

#include <algorithm>
#include <cstddef>

namespace meshioplusplus::cli::tui {

namespace {

constexpr char kEsc = '\x1b';

// The number of bytes of the UTF-8 sequence that starts with `Lead`, or 0 for
// a byte that cannot start one.
std::size_t utf8_length(unsigned char Lead) {
    if (Lead < 0x80)
        return 1;
    if (Lead >= 0xC2 && Lead <= 0xDF)
        return 2;
    if (Lead >= 0xE0 && Lead <= 0xEF)
        return 3;
    if (Lead >= 0xF0 && Lead <= 0xF4)
        return 4;
    return 0;
}

std::uint32_t utf8_decode(const std::string& rBytes, std::size_t At, std::size_t Length) {
    const unsigned char lead = static_cast<unsigned char>(rBytes[At]);
    if (Length == 1)
        return lead;
    std::uint32_t cp = lead & (0xFFu >> (Length + 1));
    for (std::size_t k = 1; k < Length; ++k)
        cp = (cp << 6) | (static_cast<unsigned char>(rBytes[At + k]) & 0x3Fu);
    return cp;
}

InputEvent make_key(KeyCode Key, std::uint32_t Ch = 0) {
    InputEvent e;
    e.mKind = EventKind::Key;
    e.mKey = Key;
    e.mChar = Ch;
    return e;
}

// xterm's modifier parameter: 1 + shift(1) + alt(2) + ctrl(4).
void apply_modifier(InputEvent& rEvent, int Param) {
    if (Param < 2)
        return;
    const int bits = Param - 1;
    rEvent.mShift = (bits & 1) != 0;
    rEvent.mAlt = (bits & 2) != 0;
    rEvent.mCtrl = (bits & 4) != 0;
}

bool is_final(char C) {
    return C >= 0x40 && C <= 0x7E;
}

// `~`-terminated keys: CSI <n> ~.
bool tilde_key(int N, KeyCode& rKey) {
    switch (N) {
        case 1:
        case 7:
            rKey = KeyCode::Home;
            return true;
        case 2:
            rKey = KeyCode::Insert;
            return true;
        case 3:
            rKey = KeyCode::Delete;
            return true;
        case 4:
        case 8:
            rKey = KeyCode::End;
            return true;
        case 5:
            rKey = KeyCode::PageUp;
            return true;
        case 6:
            rKey = KeyCode::PageDown;
            return true;
        case 11: rKey = KeyCode::F1; return true;
        case 12: rKey = KeyCode::F2; return true;
        case 13: rKey = KeyCode::F3; return true;
        case 14: rKey = KeyCode::F4; return true;
        case 15: rKey = KeyCode::F5; return true;
        case 17: rKey = KeyCode::F6; return true;
        case 18: rKey = KeyCode::F7; return true;
        case 19: rKey = KeyCode::F8; return true;
        case 20: rKey = KeyCode::F9; return true;
        case 21: rKey = KeyCode::F10; return true;
        case 23: rKey = KeyCode::F11; return true;
        case 24: rKey = KeyCode::F12; return true;
        default:
            return false;
    }
}

bool letter_key(char Final, KeyCode& rKey) {
    switch (Final) {
        case 'A': rKey = KeyCode::Up; return true;
        case 'B': rKey = KeyCode::Down; return true;
        case 'C': rKey = KeyCode::Right; return true;
        case 'D': rKey = KeyCode::Left; return true;
        case 'H': rKey = KeyCode::Home; return true;
        case 'F': rKey = KeyCode::End; return true;
        case 'P': rKey = KeyCode::F1; return true;
        case 'Q': rKey = KeyCode::F2; return true;
        case 'R': rKey = KeyCode::F3; return true;
        case 'S': rKey = KeyCode::F4; return true;
        default:
            return false;
    }
}

// Semicolon-separated numbers; a missing one is -1.
std::vector<int> parse_params(const std::string& rText) {
    std::vector<int> out;
    int cur = -1;
    for (char c : rText) {
        if (c >= '0' && c <= '9') {
            cur = (cur < 0 ? 0 : cur) * 10 + (c - '0');
            if (cur > 100000)
                cur = 100000;
        } else {
            out.push_back(cur);
            cur = -1;
        }
    }
    out.push_back(cur);
    return out;
}

}  // namespace

void InputParser::Feed(const char* pBytes, std::size_t Count) {
    mBuffer.append(pBytes, Count);
}

std::vector<InputEvent> InputParser::Take(bool Flush) {
    std::vector<InputEvent> out;
    std::size_t at = 0;
    const std::string& b = mBuffer;
    while (at < b.size()) {
        if (mInPaste) {
            // Everything up to the end marker is text.
            static const std::string kEnd = "\x1b[201~";
            const std::size_t end = b.find(kEnd, at);
            if (end == std::string::npos) {
                // Hold back a possible partial marker at the tail.
                std::size_t keep = 0;
                for (std::size_t k = std::min(kEnd.size() - 1, b.size() - at); k > 0; --k)
                    if (b.compare(b.size() - k, k, kEnd, 0, k) == 0) {
                        keep = k;
                        break;
                    }
                mPaste.append(b, at, b.size() - at - keep);
                at = b.size() - keep;
                break;
            }
            mPaste.append(b, at, end - at);
            InputEvent e;
            e.mKind = EventKind::Paste;
            e.mText = mPaste;
            out.push_back(e);
            mPaste.clear();
            mInPaste = false;
            at = end + kEnd.size();
            continue;
        }
        const unsigned char c = static_cast<unsigned char>(b[at]);
        if (c == static_cast<unsigned char>(kEsc)) {
            if (at + 1 >= b.size()) {
                if (Flush) {
                    out.push_back(make_key(KeyCode::Escape));
                    ++at;
                }
                break;  // a lone ESC may start a sequence
            }
            const char next = b[at + 1];
            if (next == '[') {
                // CSI: parameters/intermediates, then a final byte.
                std::size_t k = at + 2;
                while (k < b.size() && !is_final(b[k]))
                    ++k;
                if (k >= b.size()) {
                    if (Flush)
                        at = b.size();  // an unfinished sequence is dropped
                    break;
                }
                const char final_byte = b[k];
                const std::string body = b.substr(at + 2, k - (at + 2));
                at = k + 1;
                if (!body.empty() && body[0] == '<' && (final_byte == 'M' || final_byte == 'm')) {
                    const std::vector<int> p = parse_params(body.substr(1));
                    if (p.size() >= 3 && p[0] >= 0 && p[1] > 0 && p[2] > 0) {
                        InputEvent e;
                        e.mKind = EventKind::Mouse;
                        const int code = p[0];
                        e.mShift = (code & 4) != 0;
                        e.mAlt = (code & 8) != 0;
                        e.mCtrl = (code & 16) != 0;
                        e.mX = p[1];
                        e.mY = p[2];
                        e.mButton = code & 3;
                        if ((code & 64) != 0) {
                            e.mAction = (code & 1) != 0 ? MouseAction::WheelDown
                                                        : MouseAction::WheelUp;
                        } else if (final_byte == 'm') {
                            e.mAction = MouseAction::Release;
                        } else if ((code & 32) != 0) {
                            e.mAction = (code & 3) == 3 ? MouseAction::Move : MouseAction::Drag;
                        } else {
                            e.mAction = MouseAction::Press;
                        }
                        out.push_back(e);
                    }
                    continue;
                }
                const std::vector<int> p = parse_params(body);
                if (final_byte == '~' && !p.empty()) {
                    if (p[0] == 200) {
                        mInPaste = true;
                        mPaste.clear();
                        continue;
                    }
                    KeyCode key;
                    if (tilde_key(p[0], key)) {
                        InputEvent e = make_key(key);
                        if (p.size() > 1)
                            apply_modifier(e, p[1]);
                        out.push_back(e);
                    }
                    continue;
                }
                KeyCode key;
                if (letter_key(final_byte, key)) {
                    InputEvent e = make_key(key);
                    if (p.size() > 1)
                        apply_modifier(e, p[1]);
                    out.push_back(e);
                } else if (final_byte == 'Z') {  // back-tab
                    InputEvent e = make_key(KeyCode::Tab);
                    e.mShift = true;
                    out.push_back(e);
                }
                continue;  // an unknown CSI is consumed and ignored
            }
            if (next == 'O') {
                if (at + 2 >= b.size()) {
                    if (Flush)
                        at = b.size();
                    break;
                }
                KeyCode key;
                if (letter_key(b[at + 2], key))
                    out.push_back(make_key(key));
                at += 3;
                continue;
            }
            // ESC followed by a character: Alt+character.
            const std::size_t length = utf8_length(static_cast<unsigned char>(next));
            if (length == 0) {
                ++at;  // not a sequence: drop the ESC
                continue;
            }
            if (at + 1 + length > b.size()) {
                if (Flush)
                    at = b.size();
                break;
            }
            InputEvent e = make_key(KeyCode::Char, utf8_decode(b, at + 1, length));
            e.mAlt = true;
            out.push_back(e);
            at += 1 + length;
            continue;
        }
        // Control keys.
        if (c == '\r' || c == '\n') {
            out.push_back(make_key(KeyCode::Enter));
            ++at;
            continue;
        }
        if (c == '\t') {
            out.push_back(make_key(KeyCode::Tab));
            ++at;
            continue;
        }
        if (c == 0x7F || c == 0x08) {
            out.push_back(make_key(KeyCode::Backspace));
            ++at;
            continue;
        }
        if (c < 0x20) {  // Ctrl+letter (0x01 is Ctrl-A)
            InputEvent e = make_key(KeyCode::Char, static_cast<std::uint32_t>('a' + c - 1));
            if (c == 0x00)
                e.mChar = ' ';
            e.mCtrl = true;
            out.push_back(e);
            ++at;
            continue;
        }
        const std::size_t length = utf8_length(c);
        if (length == 0) {
            ++at;  // a stray continuation byte
            continue;
        }
        if (at + length > b.size()) {
            if (Flush)
                at = b.size();
            break;
        }
        out.push_back(make_key(KeyCode::Char, utf8_decode(b, at, length)));
        at += length;
    }
    mBuffer.erase(0, at);
    return out;
}

}  // namespace meshioplusplus::cli::tui
