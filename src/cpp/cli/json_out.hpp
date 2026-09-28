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
#pragma once

/**
 * @file cli/json_out.hpp
 * @brief The native CLI's one JSON writer.
 *
 * Every `--json` goes through `JsonOut`, so the output is always valid JSON,
 * laid out the way the Python CLI's `emit_json` (`json.dumps(indent=2)`) lays
 * it out: strings escaped (quotes, backslashes, control characters; UTF-8
 * passes through), non-finite numbers as `null`, floats in the shortest form
 * that reads back to the same double (with a `.0` on an integral value, as
 * Python prints it), and empty containers as `[]` / `{}`.
 * `tests/python/test_cli_json.py` holds the two CLIs to one shape.
 *
 * Numbers are formatted through `detail::snprintf_c`, so a process locale
 * cannot turn `1.5` into `1,5`.
 */

// System includes
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"

namespace meshioplusplus::cli {

/// A streaming, pretty-printing JSON writer (indent 2).
class JsonOut {
public:
    explicit JsonOut(std::ostream& rOs) : mOs(rOs) {}

    ~JsonOut() {
        if (mStack.empty() && mWrote)
            mOs << '\n';
    }

    JsonOut(const JsonOut&) = delete;
    JsonOut& operator=(const JsonOut&) = delete;

    void BeginObject() { Open('{'); }
    void EndObject() { Close('}'); }
    void BeginArray() { Open('['); }
    void EndArray() { Close(']'); }

    /// The key of the next value in the enclosing object.
    void Key(std::string_view Name) {
        Separate();
        WriteString(Name);
        mOs << ": ";
        mAfterKey = true;
    }

    void String(std::string_view Value) {
        Separate();
        WriteString(Value);
    }
    void Int(std::int64_t Value) {
        Separate();
        mOs << Value;
    }
    void Bool(bool Value) {
        Separate();
        mOs << (Value ? "true" : "false");
    }
    void Null() {
        Separate();
        mOs << "null";
    }
    /// A double: `null` when not finite, else the shortest exact form.
    void Number(double Value) {
        Separate();
        if (!std::isfinite(Value)) {
            mOs << "null";
            return;
        }
        mOs << Shortest(Value);
    }

    /// `"name": value` shorthands.
    void Field(std::string_view Name, std::string_view Value) {
        Key(Name);
        String(Value);
    }
    void Field(std::string_view Name, const char* pValue) { Field(Name, std::string_view(pValue)); }
    void Field(std::string_view Name, const std::string& rValue) {
        Field(Name, std::string_view(rValue));
    }
    void Field(std::string_view Name, double Value) {
        Key(Name);
        Number(Value);
    }
    void Field(std::string_view Name, std::int64_t Value) {
        Key(Name);
        Int(Value);
    }
    void Field(std::string_view Name, int Value) { Field(Name, static_cast<std::int64_t>(Value)); }
    void Field(std::string_view Name, std::size_t Value) {
        Field(Name, static_cast<std::int64_t>(Value));
    }
    void Field(std::string_view Name, bool Value) {
        Key(Name);
        Bool(Value);
    }

    /// The shortest `%.Ng` (N = 15..17) that reads back to @p Value exactly.
    static std::string Shortest(double Value) {
        char buf[40];
        for (int digits = 15; digits <= 17; ++digits) {
            detail::snprintf_c(buf, sizeof(buf), "%.*g", digits, Value);
            if (detail::parse_double(std::string(buf)) == Value)
                break;
        }
        std::string s(buf);
        if (s.find_first_of(".eEn") == std::string::npos)
            s += ".0";
        return s;
    }

private:
    struct Level {
        char mClose;
        bool mEmpty;
    };

    void Open(char c) {
        Separate();
        mOs << c;
        mStack.push_back({c == '{' ? '}' : ']', true});
    }

    void Close(char c) {
        const bool empty = mStack.back().mEmpty;
        mStack.pop_back();
        if (!empty)
            Newline();
        mOs << c;
    }

    // Before a value or key: the comma and newline its position needs.
    void Separate() {
        mWrote = true;
        if (mAfterKey) {
            mAfterKey = false;
            return;
        }
        if (mStack.empty())
            return;
        if (!mStack.back().mEmpty)
            mOs << ',';
        mStack.back().mEmpty = false;
        Newline();
    }

    void Newline() {
        mOs << '\n';
        for (std::size_t i = 0; i < mStack.size(); ++i)
            mOs << "  ";
    }

    void WriteString(std::string_view S) {
        mOs << '"';
        for (const char ch : S) {
            const auto u = static_cast<unsigned char>(ch);
            switch (ch) {
                case '"':
                    mOs << "\\\"";
                    break;
                case '\\':
                    mOs << "\\\\";
                    break;
                case '\n':
                    mOs << "\\n";
                    break;
                case '\r':
                    mOs << "\\r";
                    break;
                case '\t':
                    mOs << "\\t";
                    break;
                case '\b':
                    mOs << "\\b";
                    break;
                case '\f':
                    mOs << "\\f";
                    break;
                default:
                    if (u < 0x20) {
                        char esc[8];
                        detail::snprintf_c(esc, sizeof(esc), "\\u%04x", u);
                        mOs << esc;
                    } else {
                        mOs << ch;  // UTF-8 passes through
                    }
            }
        }
        mOs << '"';
    }

    std::ostream& mOs;
    std::vector<Level> mStack;
    bool mAfterKey = false;
    bool mWrote = false;
};

}  // namespace meshioplusplus::cli
