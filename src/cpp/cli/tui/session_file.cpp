// SPDX-License-Identifier: MIT
/// @file session_file.cpp
/// @brief Implementation of session_file.hpp.

#include "session_file.hpp"

#include <cctype>
#include <cstddef>
#include <stdexcept>

namespace meshioplusplus::cli::tui {

namespace {

constexpr const char* kPrefix = "meshio++: session: ";

void write_string(std::string& rOut, const std::string& rText) {
    rOut.push_back('"');
    for (unsigned char c : rText) {
        switch (c) {
            case '"':
                rOut += "\\\"";
                break;
            case '\\':
                rOut += "\\\\";
                break;
            case '\n':
                rOut += "\\n";
                break;
            case '\t':
                rOut += "\\t";
                break;
            case '\r':
                rOut += "\\r";
                break;
            default:
                if (c < 0x20) {
                    static const char digits[] = "0123456789abcdef";
                    rOut += "\\u00";
                    rOut.push_back(digits[c >> 4]);
                    rOut.push_back(digits[c & 15]);
                } else {
                    rOut.push_back(static_cast<char>(c));
                }
        }
    }
    rOut.push_back('"');
}

class Reader {
public:
    explicit Reader(const std::string& rText) : mrText(rText) {}

    void SkipSpace() {
        while (mAt < mrText.size() && std::isspace(static_cast<unsigned char>(mrText[mAt])))
            ++mAt;
    }
    bool Peek(char C) {
        SkipSpace();
        return mAt < mrText.size() && mrText[mAt] == C;
    }
    void Expect(char C) {
        SkipSpace();
        if (mAt >= mrText.size() || mrText[mAt] != C)
            Fail(std::string("expected '") + C + "'");
        ++mAt;
    }
    bool AtEnd() {
        SkipSpace();
        return mAt >= mrText.size();
    }

    std::string String() {
        Expect('"');
        std::string out;
        while (true) {
            if (mAt >= mrText.size())
                Fail("unterminated string");
            const char c = mrText[mAt++];
            if (c == '"')
                break;
            if (c != '\\') {
                out.push_back(c);
                continue;
            }
            if (mAt >= mrText.size())
                Fail("unterminated escape");
            const char e = mrText[mAt++];
            switch (e) {
                case '"':
                case '\\':
                case '/':
                    out.push_back(e);
                    break;
                case 'n':
                    out.push_back('\n');
                    break;
                case 't':
                    out.push_back('\t');
                    break;
                case 'r':
                    out.push_back('\r');
                    break;
                case 'u': {
                    if (mAt + 4 > mrText.size())
                        Fail("bad \\u escape");
                    unsigned v = 0;
                    for (int k = 0; k < 4; ++k) {
                        const char h = mrText[mAt++];
                        v = v * 16 +
                            static_cast<unsigned>(h >= '0' && h <= '9'   ? h - '0'
                                                  : h >= 'a' && h <= 'f' ? h - 'a' + 10
                                                  : h >= 'A' && h <= 'F' ? h - 'A' + 10
                                                                         : (Fail("bad \\u escape"), 0));
                    }
                    if (v > 0x7F)
                        Fail("only ASCII \\u escapes are supported");
                    out.push_back(static_cast<char>(v));
                    break;
                }
                default:
                    Fail(std::string("unknown escape \\") + e);
            }
        }
        return out;
    }

    long long Integer() {
        SkipSpace();
        const std::size_t start = mAt;
        if (mAt < mrText.size() && (mrText[mAt] == '-' || mrText[mAt] == '+'))
            ++mAt;
        while (mAt < mrText.size() && std::isdigit(static_cast<unsigned char>(mrText[mAt])))
            ++mAt;
        if (mAt == start || (mAt == start + 1 && !std::isdigit(static_cast<unsigned char>(
                                                      mrText[start]))))
            Fail("expected an integer");
        if (mAt - start > 18)
            Fail("integer out of range");
        return std::stoll(mrText.substr(start, mAt - start));
    }

    [[noreturn]] void Fail(const std::string& rWhat) const {
        throw std::invalid_argument(std::string(kPrefix) + rWhat + " (at byte " +
                                    std::to_string(mAt) + ")");
    }

private:
    const std::string& mrText;
    std::size_t mAt = 0;
};

}  // namespace

std::string session_to_json(const SessionFile& rSession) {
    // Keys in sorted order: flags, step, version.
    std::string out = "{\n  \"flags\": [";
    for (std::size_t i = 0; i < rSession.mFlags.size(); ++i) {
        out += i == 0 ? "\n    " : ",\n    ";
        write_string(out, rSession.mFlags[i]);
    }
    out += rSession.mFlags.empty() ? "]" : "\n  ]";
    out += ",\n  \"step\": " + std::to_string(rSession.mStep);
    out += ",\n  \"version\": " + std::to_string(kSessionVersion) + "\n}\n";
    return out;
}

SessionFile session_from_json(const std::string& rText) {
    Reader in(rText);
    SessionFile session;
    bool have_version = false;
    bool have_flags = false;
    bool have_step = false;
    in.Expect('{');
    if (!in.Peek('}')) {
        while (true) {
            const std::string key = in.String();
            in.Expect(':');
            if (key == "version") {
                if (have_version)
                    in.Fail("duplicate key 'version'");
                have_version = true;
                const long long v = in.Integer();
                if (v != kSessionVersion)
                    throw std::invalid_argument(std::string(kPrefix) + "version " +
                                                std::to_string(v) + " is not supported (this build "
                                                "reads version " +
                                                std::to_string(kSessionVersion) + ")");
            } else if (key == "step") {
                if (have_step)
                    in.Fail("duplicate key 'step'");
                have_step = true;
                session.mStep = in.Integer();
                if (session.mStep < 0)
                    in.Fail("'step' must not be negative");
            } else if (key == "flags") {
                if (have_flags)
                    in.Fail("duplicate key 'flags'");
                have_flags = true;
                in.Expect('[');
                if (!in.Peek(']')) {
                    while (true) {
                        session.mFlags.push_back(in.String());
                        if (in.Peek(',')) {
                            in.Expect(',');
                            continue;
                        }
                        break;
                    }
                }
                in.Expect(']');
            } else {
                in.Fail("unknown key '" + key + "' (expected flags, step and version)");
            }
            if (in.Peek(',')) {
                in.Expect(',');
                continue;
            }
            break;
        }
    }
    in.Expect('}');
    if (!in.AtEnd())
        in.Fail("text after the closing brace");
    if (!have_version)
        in.Fail("missing 'version'");
    return session;
}

}  // namespace meshioplusplus::cli::tui
