//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// The JSON document model of the JSON format readers. See json_value.hpp.

// System includes
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/exceptions.hpp"

// Project includes (private, not installed)
#include "json_value.hpp"

namespace meshioplusplus::detail {

namespace {

constexpr int kJsonMaxDepth = 256;

// The offset of the first byte of @p rText that is not well-formed UTF-8
// (overlong forms, surrogates and code points past U+10FFFF included, as
// Python's strict decoder rejects them), or npos.
std::size_t json_invalid_utf8_at(std::string_view Text) noexcept {
    const auto* p = reinterpret_cast<const unsigned char*>(Text.data());
    const std::size_t n = Text.size();
    for (std::size_t i = 0; i < n;) {
        const unsigned char c = p[i];
        if (c < 0x80) {
            ++i;
            continue;
        }
        std::size_t len = 0;
        unsigned char lo = 0x80, hi = 0xBF;  // the second byte's range
        if (c >= 0xC2 && c <= 0xDF) {
            len = 2;
        } else if (c >= 0xE0 && c <= 0xEF) {
            len = 3;
            if (c == 0xE0)
                lo = 0xA0;  // no overlong forms
            else if (c == 0xED)
                hi = 0x9F;  // no surrogates
        } else if (c >= 0xF0 && c <= 0xF4) {
            len = 4;
            if (c == 0xF0)
                lo = 0x90;
            else if (c == 0xF4)
                hi = 0x8F;  // nothing past U+10FFFF
        } else {
            return i;
        }
        if (i + len > n || p[i + 1] < lo || p[i + 1] > hi)
            return i;
        for (std::size_t k = 2; k < len; ++k)
            if (p[i + k] < 0x80 || p[i + k] > 0xBF)
                return i;
        i += len;
    }
    return std::string_view::npos;
}

class JsonParser {
public:
    JsonParser(std::string_view Text, const char* pFormat) : mText(Text), mpFormat(pFormat) {
        if (mText.size() >= 3 && static_cast<unsigned char>(mText[0]) == 0xEF &&
            static_cast<unsigned char>(mText[1]) == 0xBB &&
            static_cast<unsigned char>(mText[2]) == 0xBF)
            mPos = 3;
        // JSON text is UTF-8 (RFC 8259), and a string read from it reaches
        // callers (Python among them) that must be able to decode it.
        const std::size_t bad = json_invalid_utf8_at(mText);
        if (bad != std::string_view::npos) {
            mPos = bad;
            Fail("the text is not valid UTF-8");
        }
    }

    JsonValue ParseDocument() {
        SkipSpace();
        JsonValue v = ParseValue(0);
        SkipSpace();
        if (mPos != mText.size())
            Fail("unexpected content after the document");
        return v;
    }

private:
    [[noreturn]] void Fail(const std::string& rWhat) const {
        throw ReadError(std::string("meshio++: ") + mpFormat + ": malformed JSON at byte " +
                        std::to_string(mPos) + ": " + rWhat);
    }

    bool AtEnd() const noexcept { return mPos >= mText.size(); }
    char Peek() const noexcept { return AtEnd() ? '\0' : mText[mPos]; }

    void SkipSpace() {
        while (!AtEnd()) {
            const char c = mText[mPos];
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
                ++mPos;
            } else if (c == '/' && mPos + 1 < mText.size() && mText[mPos + 1] == '/') {
                while (!AtEnd() && mText[mPos] != '\n')
                    ++mPos;
            } else if (c == '/' && mPos + 1 < mText.size() && mText[mPos + 1] == '*') {
                const std::size_t close = mText.find("*/", mPos + 2);
                if (close == std::string_view::npos)
                    Fail("unterminated comment");
                mPos = close + 2;
            } else {
                break;
            }
        }
    }

    void Expect(std::string_view Word) {
        if (mText.substr(mPos, Word.size()) != Word)
            Fail("expected '" + std::string(Word) + "'");
        mPos += Word.size();
    }

    JsonValue ParseValue(int Depth) {
        if (Depth > kJsonMaxDepth)
            Fail("nesting deeper than " + std::to_string(kJsonMaxDepth));
        switch (Peek()) {
            case '{':
                return ParseObject(Depth);
            case '[':
                return ParseArray(Depth);
            case '"':
                return JsonValue::String(ParseString());
            case 't':
                Expect("true");
                return JsonValue::Boolean(true);
            case 'f':
                Expect("false");
                return JsonValue::Boolean(false);
            case 'n':
                Expect("null");
                return JsonValue::Null();
            default:
                if (Peek() == '-' || (Peek() >= '0' && Peek() <= '9'))
                    return ParseNumber();
                Fail(AtEnd() ? "unexpected end of input" : "unexpected character");
        }
    }

    JsonValue ParseObject(int Depth) {
        ++mPos;
        JsonValue v = JsonValue::Object();
        SkipSpace();
        if (Peek() == '}') {
            ++mPos;
            return v;
        }
        for (;;) {
            SkipSpace();
            if (Peek() != '"')
                Fail("expected a member name");
            std::string key = ParseString();
            SkipSpace();
            if (Peek() != ':')
                Fail("expected ':'");
            ++mPos;
            SkipSpace();
            JsonValue item = ParseValue(Depth + 1);
            v.mMembers.emplace_back(std::move(key), std::move(item));
            SkipSpace();
            if (Peek() == ',') {
                ++mPos;
                continue;
            }
            if (Peek() == '}') {
                ++mPos;
                return v;
            }
            Fail("expected ',' or '}'");
        }
    }

    JsonValue ParseArray(int Depth) {
        ++mPos;
        JsonValue v = JsonValue::Array();
        SkipSpace();
        if (Peek() == ']') {
            ++mPos;
            return v;
        }
        for (;;) {
            SkipSpace();
            v.mItems.push_back(ParseValue(Depth + 1));
            SkipSpace();
            if (Peek() == ',') {
                ++mPos;
                continue;
            }
            if (Peek() == ']') {
                ++mPos;
                return v;
            }
            Fail("expected ',' or ']'");
        }
    }

    unsigned ParseHex4() {
        if (mPos + 4 > mText.size())
            Fail("truncated \\u escape");
        unsigned code = 0;
        for (int i = 0; i < 4; ++i) {
            const char c = mText[mPos++];
            code <<= 4;
            if (c >= '0' && c <= '9')
                code |= static_cast<unsigned>(c - '0');
            else if (c >= 'a' && c <= 'f')
                code |= static_cast<unsigned>(c - 'a' + 10);
            else if (c >= 'A' && c <= 'F')
                code |= static_cast<unsigned>(c - 'A' + 10);
            else
                Fail("bad \\u escape");
        }
        return code;
    }

    static void AppendUtf8(std::string& rOut, unsigned Code) {
        if (Code < 0x80) {
            rOut += static_cast<char>(Code);
        } else if (Code < 0x800) {
            rOut += static_cast<char>(0xC0 | (Code >> 6));
            rOut += static_cast<char>(0x80 | (Code & 0x3F));
        } else if (Code < 0x10000) {
            rOut += static_cast<char>(0xE0 | (Code >> 12));
            rOut += static_cast<char>(0x80 | ((Code >> 6) & 0x3F));
            rOut += static_cast<char>(0x80 | (Code & 0x3F));
        } else {
            rOut += static_cast<char>(0xF0 | (Code >> 18));
            rOut += static_cast<char>(0x80 | ((Code >> 12) & 0x3F));
            rOut += static_cast<char>(0x80 | ((Code >> 6) & 0x3F));
            rOut += static_cast<char>(0x80 | (Code & 0x3F));
        }
    }

    std::string ParseString() {
        ++mPos;  // the opening quote
        std::string out;
        for (;;) {
            if (AtEnd())
                Fail("unterminated string");
            const char c = mText[mPos++];
            if (c == '"')
                return out;
            if (static_cast<unsigned char>(c) < 0x20)
                Fail("control character in a string");
            if (c != '\\') {
                out += c;
                continue;
            }
            if (AtEnd())
                Fail("unterminated escape");
            const char e = mText[mPos++];
            switch (e) {
                case '"':
                case '\\':
                case '/':
                    out += e;
                    break;
                case 'b':
                    out += '\b';
                    break;
                case 'f':
                    out += '\f';
                    break;
                case 'n':
                    out += '\n';
                    break;
                case 'r':
                    out += '\r';
                    break;
                case 't':
                    out += '\t';
                    break;
                case 'u': {
                    unsigned code = ParseHex4();
                    if (code >= 0xD800 && code <= 0xDBFF) {
                        if (mPos + 2 > mText.size() || mText[mPos] != '\\' || mText[mPos + 1] != 'u')
                            Fail("unpaired surrogate");
                        mPos += 2;
                        const unsigned low = ParseHex4();
                        if (low < 0xDC00 || low > 0xDFFF)
                            Fail("unpaired surrogate");
                        code = 0x10000 + ((code - 0xD800) << 10) + (low - 0xDC00);
                    } else if (code >= 0xDC00 && code <= 0xDFFF) {
                        Fail("unpaired surrogate");
                    }
                    AppendUtf8(out, code);
                    break;
                }
                default:
                    Fail("bad escape");
            }
        }
    }

    JsonValue ParseNumber() {
        const std::size_t start = mPos;
        auto digit = [this]() { return !AtEnd() && mText[mPos] >= '0' && mText[mPos] <= '9'; };
        bool integral = true;
        if (Peek() == '-')
            ++mPos;
        if (Peek() == '0') {
            ++mPos;
        } else if (digit()) {
            while (digit())
                ++mPos;
        } else {
            Fail("bad number");
        }
        if (Peek() == '.') {
            integral = false;
            ++mPos;
            if (!digit())
                Fail("bad number");
            while (digit())
                ++mPos;
        }
        if (Peek() == 'e' || Peek() == 'E') {
            integral = false;
            ++mPos;
            if (Peek() == '+' || Peek() == '-')
                ++mPos;
            if (!digit())
                Fail("bad number");
            while (digit())
                ++mPos;
        }
        const std::string lexeme(mText.substr(start, mPos - start));
        JsonValue v;
        v.mKind = JsonValue::Kind::Number;
        v.mNumber = parse_double(lexeme);
        if (integral) {
            // Exact int64 when it fits; otherwise the double alone.
            const bool negative = lexeme[0] == '-';
            std::uint64_t magnitude = 0;
            bool fits = true;
            for (std::size_t i = negative ? 1 : 0; i < lexeme.size(); ++i) {
                const auto d = static_cast<std::uint64_t>(lexeme[i] - '0');
                if (magnitude > (std::numeric_limits<std::uint64_t>::max() - d) / 10) {
                    fits = false;
                    break;
                }
                magnitude = magnitude * 10 + d;
            }
            const std::uint64_t limit =
                negative ? static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max()) + 1
                         : static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
            if (fits && magnitude <= limit && !(negative && magnitude == 0)) {
                v.mIsInteger = true;
                v.mInteger = negative ? static_cast<std::int64_t>(0 - magnitude)
                                      : static_cast<std::int64_t>(magnitude);
            }
        }
        return v;
    }

    std::string_view mText;
    const char* mpFormat;
    std::size_t mPos = 0;
};

void json_append_string(std::string& rOut, const std::string& rText) {
    static const char kHex[] = "0123456789abcdef";
    rOut += '"';
    for (const char ch : rText) {
        const auto c = static_cast<unsigned char>(ch);
        switch (c) {
            case '"':
                rOut += "\\\"";
                break;
            case '\\':
                rOut += "\\\\";
                break;
            case '\b':
                rOut += "\\b";
                break;
            case '\f':
                rOut += "\\f";
                break;
            case '\n':
                rOut += "\\n";
                break;
            case '\r':
                rOut += "\\r";
                break;
            case '\t':
                rOut += "\\t";
                break;
            default:
                if (c < 0x20) {
                    rOut += "\\u00";
                    rOut += kHex[c >> 4];
                    rOut += kHex[c & 0xF];
                } else {
                    rOut += ch;
                }
        }
    }
    rOut += '"';
}

void json_write(std::string& rOut, const JsonValue& rV, int Indent, int Level) {
    auto newline = [&](int L) {
        if (Indent <= 0)
            return;
        rOut += '\n';
        rOut.append(static_cast<std::size_t>(Indent) * static_cast<std::size_t>(L), ' ');
    };
    switch (rV.mKind) {
        case JsonValue::Kind::Null:
            rOut += "null";
            return;
        case JsonValue::Kind::Bool:
            rOut += rV.mBool ? "true" : "false";
            return;
        case JsonValue::Kind::Number:
            if (rV.mIsInteger)
                rOut += std::to_string(rV.mInteger);
            else
                append_json_number(rOut, rV.mNumber);
            return;
        case JsonValue::Kind::String:
            json_append_string(rOut, rV.mString);
            return;
        case JsonValue::Kind::Array:
            if (rV.mItems.empty()) {
                rOut += "[]";
                return;
            }
            rOut += '[';
            for (std::size_t i = 0; i < rV.mItems.size(); ++i) {
                if (i)
                    rOut += ',';
                newline(Level + 1);
                json_write(rOut, rV.mItems[i], Indent, Level + 1);
            }
            newline(Level);
            rOut += ']';
            return;
        case JsonValue::Kind::Object:
            if (rV.mMembers.empty()) {
                rOut += "{}";
                return;
            }
            rOut += '{';
            for (std::size_t i = 0; i < rV.mMembers.size(); ++i) {
                if (i)
                    rOut += ',';
                newline(Level + 1);
                json_append_string(rOut, rV.mMembers[i].mKey);
                rOut += Indent > 0 ? ": " : ":";
                json_write(rOut, rV.mMembers[i].mValue, Indent, Level + 1);
            }
            newline(Level);
            rOut += '}';
            return;
    }
}

}  // namespace

const JsonValue* JsonValue::Find(std::string_view Key) const noexcept {
    for (const auto& m : mMembers)
        if (m.mKey == Key)
            return &m.mValue;
    return nullptr;
}

JsonValue JsonValue::Null() {
    return {};
}

JsonValue JsonValue::Boolean(bool Value) {
    JsonValue v;
    v.mKind = Kind::Bool;
    v.mBool = Value;
    return v;
}

JsonValue JsonValue::Number(double Value) {
    JsonValue v;
    v.mKind = Kind::Number;
    v.mNumber = Value;
    return v;
}

JsonValue JsonValue::Integer(std::int64_t Value) {
    JsonValue v;
    v.mKind = Kind::Number;
    v.mNumber = static_cast<double>(Value);
    v.mIsInteger = true;
    v.mInteger = Value;
    return v;
}

JsonValue JsonValue::String(std::string Value) {
    JsonValue v;
    v.mKind = Kind::String;
    v.mString = std::move(Value);
    return v;
}

JsonValue JsonValue::Array() {
    JsonValue v;
    v.mKind = Kind::Array;
    return v;
}

JsonValue JsonValue::Object() {
    JsonValue v;
    v.mKind = Kind::Object;
    return v;
}

JsonValue& JsonValue::Push(JsonValue Value) {
    mItems.push_back(std::move(Value));
    return mItems.back();
}

JsonValue& JsonValue::Add(std::string Key, JsonValue Value) {
    mMembers.emplace_back(std::move(Key), std::move(Value));
    return mMembers.back().mValue;
}

bool json_equal(const JsonValue& rA, const JsonValue& rB) noexcept {
    if (rA.mKind != rB.mKind)
        return false;
    switch (rA.mKind) {
        case JsonValue::Kind::Null:
            return true;
        case JsonValue::Kind::Bool:
            return rA.mBool == rB.mBool;
        case JsonValue::Kind::Number:
            if (rA.mIsInteger && rB.mIsInteger)
                return rA.mInteger == rB.mInteger;
            return rA.mNumber == rB.mNumber;
        case JsonValue::Kind::String:
            return rA.mString == rB.mString;
        case JsonValue::Kind::Array:
            if (rA.mItems.size() != rB.mItems.size())
                return false;
            for (std::size_t i = 0; i < rA.mItems.size(); ++i)
                if (!json_equal(rA.mItems[i], rB.mItems[i]))
                    return false;
            return true;
        case JsonValue::Kind::Object:
            if (rA.mMembers.size() != rB.mMembers.size())
                return false;
            for (std::size_t i = 0; i < rA.mMembers.size(); ++i)
                if (rA.mMembers[i].mKey != rB.mMembers[i].mKey ||
                    !json_equal(rA.mMembers[i].mValue, rB.mMembers[i].mValue))
                    return false;
            return true;
    }
    return false;
}

JsonValue parse_json(std::string_view Text, const char* pFormat) {
    return JsonParser(Text, pFormat).ParseDocument();
}

std::string write_json(const JsonValue& rValue, int Indent) {
    std::string out;
    json_write(out, rValue, Indent, 0);
    return out;
}

void append_json_number(std::string& rOut, double Value) {
    if (!std::isfinite(Value)) {
        rOut += "null";  // JSON has no spelling for them (nlohmann's choice too)
        return;
    }
    if (Value == std::floor(Value) && std::fabs(Value) < 9007199254740992.0 &&
        !(Value == 0.0 && std::signbit(Value))) {
        rOut += std::to_string(static_cast<std::int64_t>(Value));
        return;
    }
    char buf[40];
    for (const char* fmt : {"%.15g", "%.16g", "%.17g"}) {
        snprintf_c(buf, sizeof(buf), fmt, Value);
        if (parse_double(std::string(buf)) == Value)
            break;
    }
    rOut += buf;
}

}  // namespace meshioplusplus::detail
