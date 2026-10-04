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
 * @file detail/text_cursor.hpp
 * @brief Text readers' lines, tokens and numbers as views over one buffer.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * Twenty-odd text readers held the whole file as a `std::vector<std::string>`
 * of lines and split each line into `std::string` tokens through a
 * `std::istringstream`: a heap string per line and per token. Here a file is
 * read once (`FileSource`), its lines are `string_view`s into it
 * (`split_lines`, the lines `std::getline` gives), a line's tokens are views
 * into it (`split_blanks`, the tokens `istringstream >> std::string` gives),
 * and a token parses in place with the semantics the readers used before:
 * `parse_double_token` is `parse_double` over the whole token, and
 * `parse_int_token` is `strtoll(…, 10)` over the whole token -- a leading `+`
 * accepted, an out-of-range value saturated. `TextCursor` shares bounded line
 * and prefix-number positioning; `RecordCursor` traverses split records while
 * format adapters retain their comments, quoting and diagnostics. Roadmap §3,
 * "A shared tokenizer and number path".
 */

// System includes
#include <charconv>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cerrno>
#include <limits>
#include <string>
#include <string_view>
#include <span>
#include <system_error>
#include <type_traits>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"

namespace meshioplusplus {
namespace detail {

/// Whether `c` separates tokens: `std::isspace` in the classic locale.
inline bool text_is_blank(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
}

/**
 * @brief Core-private byte cursor shared by format-specific grammar adapters.
 * Lines and prefix numbers stay bounded even when the source has no terminator.
 * A prefix parse advances only through the accepted prefix, not the whole token.
 */
class TextCursor {
public:
    explicit TextCursor(std::string_view Text) : mBuf(Text) {}
    template <class T>
        requires(std::is_same_v<std::remove_cvref_t<T>, std::string> &&
                 !std::is_lvalue_reference_v<T>)
    explicit TextCursor(T&&) = delete;
    bool AtEnd() const { return mPos >= mBuf.size(); }
    std::size_t Pos() const { return mPos; }
    void Seek(std::size_t Pos) { mPos = Pos; }
    std::string_view Data() const { return mBuf; }

    std::string_view Line(bool StripCr = false) {
        if (AtEnd())
            return {};
        const auto eol = mBuf.find('\n', mPos);
        const auto end = eol == std::string_view::npos ? mBuf.size() : eol;
        auto line = mBuf.substr(mPos, end - mPos);
        mPos = eol == std::string_view::npos ? end : end + 1;
        if (StripCr && !line.empty() && line.back() == '\r')
            line.remove_suffix(1);
        return line;
    }

    bool DoublePrefix(double& rValue) {
#ifdef MESHIOPLUSPLUS_HAS_FAST_FROM_CHARS
        auto pos = NumberStart();
        if (pos < mBuf.size() && mBuf[pos] == '+')
            ++pos;
        const auto hex = pos < mBuf.size() && mBuf[pos] == '-' ? pos + 1 : pos;
        const bool hexadecimal = hex + 1 < mBuf.size() && mBuf[hex] == '0' &&
                                 (mBuf[hex + 1] == 'x' || mBuf[hex + 1] == 'X');
        if (pos < mBuf.size() && !hexadecimal) {
            double value = 0.0;
            const auto parsed =
                std::from_chars(mBuf.data() + pos, mBuf.data() + mBuf.size(), value);
            if (parsed.ec == std::errc{}) {
                rValue = value;
                mPos = static_cast<std::size_t>(parsed.ptr - mBuf.data());
                return true;
            }
        }
#endif
        // Hexadecimal, out-of-range and platform fallback paths use the exact
        // locale-independent C parser over a bounded, terminated copy.
        return NumberPrefix(rValue, [](const char* pFirst, const char*& rpEnd) {
            return parse_double(pFirst, rpEnd);
        });
    }
    bool IntPrefix(std::int64_t& rValue) {
        auto pos = NumberStart();
        const bool negative = pos < mBuf.size() && mBuf[pos] == '-';
        if (pos < mBuf.size() && (mBuf[pos] == '+' || negative))
            ++pos;
        if (pos < mBuf.size() && mBuf[pos] >= '0' && mBuf[pos] <= '9') {
            std::uint64_t magnitude = 0;
            const auto parsed =
                std::from_chars(mBuf.data() + pos, mBuf.data() + mBuf.size(), magnitude);
            constexpr auto maximum =
                static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
            if (parsed.ec == std::errc::result_out_of_range ||
                magnitude > maximum + (negative ? 1u : 0u)) {
                rValue = negative ? std::numeric_limits<std::int64_t>::min()
                                  : std::numeric_limits<std::int64_t>::max();
                errno = ERANGE;  // the saturating strtoll contract
            } else {
                rValue = negative ? static_cast<std::int64_t>(0 - magnitude)
                                  : static_cast<std::int64_t>(magnitude);
            }
            mPos = static_cast<std::size_t>(parsed.ptr - mBuf.data());
            return true;
        }
        return NumberPrefix(rValue, [](const char* pFirst, const char*& rpEnd) {
            char* end = nullptr;
            const auto value = std::strtoll(pFirst, &end, 10);
            rpEnd = end;
            return static_cast<std::int64_t>(value);
        });
    }

    // Binary-aware adapters retain direct checked byte positioning.
    std::string_view mBuf;
    std::size_t mPos = 0;

private:
    std::size_t NumberStart() const {
        auto first = mPos;
        while (first < mBuf.size() && text_is_blank(mBuf[first]))
            ++first;
        return first;
    }

    template <class T, class Parser>
    bool NumberPrefix(T& rValue, Parser Parse) {
        const auto first = NumberStart();
        if (first >= mBuf.size())
            return false;
        std::size_t last = first;
        while (last < mBuf.size() && !text_is_blank(mBuf[last]))
            ++last;
        const auto token = mBuf.substr(first, last - first);
        char small[64];
        std::string large;
        const char* start;
        if (token.size() < sizeof small) {
            token.copy(small, token.size());
            small[token.size()] = '\0';
            start = small;
        } else {
            large.assign(token);
            start = large.c_str();
        }
        const char* end = nullptr;
        const auto value = Parse(start, end);
        if (end == start)
            return false;
        rValue = value;
        mPos = first + static_cast<std::size_t>(end - start);
        return true;
    }
};

/// A non-owning cursor over already split records; grammar adapters check EOF.
template <class T>
class RecordCursor {
public:
    RecordCursor() = default;
    explicit RecordCursor(std::span<const T> Records) : mRecords(Records) {}
    bool Done() const { return mPos >= mRecords.size(); }
    std::size_t Remaining() const { return Done() ? 0 : mRecords.size() - mPos; }
    std::size_t Pos() const { return mPos; }
    const T& Peek(std::size_t Ahead = 0) const { return mRecords[mPos + Ahead]; }
    const T& Next() { return mRecords[mPos++]; }

private:
    std::span<const T> mRecords;
    std::size_t mPos = 0;
};

/**
 * @brief The lines `std::getline` reads from @p Text, as views into it: split
 * on `'\n'` (a `'\r'` before it is kept, as `getline` keeps it), with no empty
 * last line after a final newline.
 */
inline std::vector<std::string_view> split_lines(std::string_view Text) {
    std::vector<std::string_view> lines;
    std::size_t pos = 0;
    while (pos < Text.size()) {
        std::size_t eol = Text.find('\n', pos);
        if (eol == std::string_view::npos)
            eol = Text.size();
        lines.push_back(Text.substr(pos, eol - pos));
        pos = eol + 1;
    }
    return lines;
}

/// The tokens `istringstream >> std::string` reads from @p Line, as views,
/// appended to @p rOut (cleared first).
inline void split_blanks(std::string_view Line, std::vector<std::string_view>& rOut) {
    rOut.clear();
    std::size_t i = 0;
    const std::size_t n = Line.size();
    while (true) {
        while (i < n && text_is_blank(Line[i]))
            ++i;
        if (i >= n)
            return;
        const std::size_t b = i;
        while (i < n && !text_is_blank(Line[i]))
            ++i;
        rOut.push_back(Line.substr(b, i - b));
    }
}

/// `split_blanks` returning its tokens.
inline std::vector<std::string_view> split_blanks(std::string_view Line) {
    std::vector<std::string_view> out;
    split_blanks(Line, out);
    return out;
}

/**
 * @brief `parse_double` over the whole of @p Token: true, with the value, when
 * the entire token is one number (what `parse_double(s.c_str(), end)` with
 * `end == s.c_str() + s.size()` accepted), false otherwise.
 */
inline bool parse_double_token(std::string_view Token, double& rOut) {
    if (Token.empty())
        return false;
    // parse_double reads a NUL-terminated string; a token is copied into one.
    // Tokens are short, so this costs little next to the parse itself.
    char small[64];
    std::string large;
    const char* first;
    if (Token.size() < sizeof small) {
        Token.copy(small, Token.size());
        small[Token.size()] = '\0';
        first = small;
    } else {
        large.assign(Token);
        first = large.c_str();
    }
    const char* end = nullptr;
    const double v = parse_double(first, end);
    if (end != first + Token.size())
        return false;
    rOut = v;
    return true;
}

/**
 * @brief `strtoll(…, 10)` over the whole of @p Token: true, with the value,
 * when the entire token is one base-10 integer. Like `strtoll`, a leading `+`
 * is accepted and a value outside `int64` saturates.
 */
inline bool parse_int_token(std::string_view Token, std::int64_t& rOut) {
    std::size_t i = 0;
    while (i < Token.size() && text_is_blank(Token[i]))
        ++i;  // strtoll skips leading blanks
    bool neg = false;
    if (i < Token.size() && (Token[i] == '+' || Token[i] == '-')) {
        neg = Token[i] == '-';
        ++i;
    }
    if (i >= Token.size() || Token[i] < '0' || Token[i] > '9')
        return false;
    std::uint64_t mag = 0;
    const char* p = Token.data() + i;
    const char* last = Token.data() + Token.size();
    const auto r = std::from_chars(p, last, mag);
    if (r.ptr != last)
        return false;
    constexpr std::uint64_t kMax =
        static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
    if (r.ec == std::errc::result_out_of_range || mag > kMax + (neg ? 1 : 0)) {
        rOut = neg ? std::numeric_limits<std::int64_t>::min()
                   : std::numeric_limits<std::int64_t>::max();
        return true;
    }
    rOut = neg ? static_cast<std::int64_t>(0 - mag) : static_cast<std::int64_t>(mag);
    return true;
}

/**
 * @brief `strtoll(token, nullptr, 10)` on @p Token: the longest base-10 prefix
 * (after blanks and a sign), 0 when there is none, saturated when out of
 * range -- the lenient parse, which ignores what follows the number.
 */
inline std::int64_t strtoll_token(std::string_view Token) {
    std::size_t i = 0;
    while (i < Token.size() && text_is_blank(Token[i]))
        ++i;
    bool neg = false;
    if (i < Token.size() && (Token[i] == '+' || Token[i] == '-')) {
        neg = Token[i] == '-';
        ++i;
    }
    std::size_t j = i;
    while (j < Token.size() && Token[j] >= '0' && Token[j] <= '9')
        ++j;
    if (j == i)
        return 0;
    std::uint64_t mag = 0;
    const auto r = std::from_chars(Token.data() + i, Token.data() + j, mag);
    constexpr std::uint64_t kMax =
        static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max());
    if (r.ec == std::errc::result_out_of_range || mag > kMax + (neg ? 1 : 0))
        return neg ? std::numeric_limits<std::int64_t>::min()
                   : std::numeric_limits<std::int64_t>::max();
    return neg ? static_cast<std::int64_t>(0 - mag) : static_cast<std::int64_t>(mag);
}

/**
 * @brief `strtoull(token, nullptr, 10)` on @p Token: the longest base-10
 * prefix after blanks and a sign, 0 when there is none, `UINT64_MAX` when out
 * of range, and -- as `strtoull` does -- a negative value wrapped modulo 2^64.
 */
inline std::uint64_t strtoull_token(std::string_view Token) {
    std::size_t i = 0;
    while (i < Token.size() && text_is_blank(Token[i]))
        ++i;
    bool neg = false;
    if (i < Token.size() && (Token[i] == '+' || Token[i] == '-')) {
        neg = Token[i] == '-';
        ++i;
    }
    std::size_t j = i;
    while (j < Token.size() && Token[j] >= '0' && Token[j] <= '9')
        ++j;
    if (j == i)
        return 0;
    std::uint64_t mag = 0;
    const auto r = std::from_chars(Token.data() + i, Token.data() + j, mag);
    if (r.ec == std::errc::result_out_of_range)
        return std::numeric_limits<std::uint64_t>::max();
    return neg ? 0 - mag : mag;
}

/**
 * @brief `parse_double(token)` on @p Token: the longest numeric prefix, 0.0
 * when there is none -- the lenient parse, which ignores what follows.
 */
inline double parse_double_prefix(std::string_view Token) {
    char small[64];
    std::string large;
    const char* first;
    if (Token.size() < sizeof small) {
        Token.copy(small, Token.size());
        small[Token.size()] = '\0';
        first = small;
    } else {
        large.assign(Token);
        first = large.c_str();
    }
    const char* end = nullptr;
    return parse_double(first, end);
}

/**
 * @brief `std::istringstream` extraction over a view, with no stream: the
 * same values and the same fail state as `operator>>` in the classic locale.
 *
 * Readers that split a line into numbers through a stream built one per line
 * (a copy of the line, a locale-aware parse per value) take this instead,
 * unchanged. A double is read as libstdc++'s `num_get` does: blanks skipped,
 * then a sign, digits with at most one decimal point, and -- after a digit --
 * an exponent (`e`/`E`, a sign, digits) are accumulated and converted as
 * `strtod` would; a conversion that fails stores 0, one that overflows stores
 * the largest finite value of that sign, and either sets the fail state. An
 * integer accumulates a sign and digits, and out of range stores the bound and
 * fails. At the end of the text an extraction fails and leaves the value
 * alone, and once failed every further extraction does nothing, as with a
 * stream. `test_text_cursor.cpp` checks it against `std::istringstream`.
 */
class TextStream {
public:
    explicit TextStream(std::string_view Text) : mText(Text) {}
    // A const temporary cannot be moved into mOwned; do not silently borrow it.
    explicit TextStream(const std::string&&) = delete;
    explicit TextStream(const char* pText) : mText(pText ? pText : "") {}
    explicit TextStream(const std::string& rText) : mText(rText) {}
    /// A temporary is kept, so the view never outlives it.
    explicit TextStream(std::string&& rText) : mOwned(std::move(rText)), mText(mOwned) {}
    TextStream(const TextStream&) = delete;
    TextStream& operator=(const TextStream&) = delete;

    explicit operator bool() const { return !mFail; }
    bool operator!() const { return mFail; }
    bool fail() const { return mFail; }

    TextStream& operator>>(double& rValue) {
        if (!Sentry())
            return *this;
        const std::size_t b = mPos;
        if (mPos < mText.size() && (mText[mPos] == '+' || mText[mPos] == '-'))
            ++mPos;
        bool digits = false, point = false;
        while (mPos < mText.size()) {
            const char c = mText[mPos];
            if (c >= '0' && c <= '9') {
                digits = true;
                ++mPos;
            } else if (c == '.' && !point) {
                point = true;
                ++mPos;
            } else {
                break;
            }
        }
        if (digits && mPos < mText.size() && (mText[mPos] == 'e' || mText[mPos] == 'E')) {
            ++mPos;
            if (mPos < mText.size() && (mText[mPos] == '+' || mText[mPos] == '-'))
                ++mPos;
            while (mPos < mText.size() && mText[mPos] >= '0' && mText[mPos] <= '9')
                ++mPos;
        }
        double v = 0.0;
        if (!parse_double_token(mText.substr(b, mPos - b), v)) {
            rValue = 0.0;
            mFail = true;
        } else if (v == std::numeric_limits<double>::infinity()) {
            rValue = std::numeric_limits<double>::max();
            mFail = true;
        } else if (v == -std::numeric_limits<double>::infinity()) {
            rValue = -std::numeric_limits<double>::max();
            mFail = true;
        } else {
            rValue = v;
        }
        return *this;
    }

    TextStream& operator>>(long long& rValue) { return ExtractInt(rValue); }
    TextStream& operator>>(long& rValue) { return ExtractInt(rValue); }
    TextStream& operator>>(int& rValue) { return ExtractInt(rValue); }
    TextStream& operator>>(unsigned long long& rValue) { return ExtractUnsigned(rValue); }
    TextStream& operator>>(unsigned long& rValue) { return ExtractUnsigned(rValue); }
    TextStream& operator>>(unsigned& rValue) { return ExtractUnsigned(rValue); }

    /// `std::getline(stream, rLine, Delim)`: the characters up to @p Delim
    /// (consumed, not stored) or the end; fails only when none are left.
    friend TextStream& getline(TextStream& rIn, std::string& rLine, char Delim = '\n') {
        if (rIn.mFail)
            return rIn;
        if (rIn.mPos >= rIn.mText.size()) {
            rIn.mFail = true;
            return rIn;
        }
        std::size_t e = rIn.mText.find(Delim, rIn.mPos);
        if (e == std::string_view::npos)
            e = rIn.mText.size();
        rLine.assign(rIn.mText.substr(rIn.mPos, e - rIn.mPos));
        rIn.mPos = e < rIn.mText.size() ? e + 1 : e;
        return rIn;
    }

    TextStream& operator>>(std::string& rValue) {
        if (!Sentry())
            return *this;
        const std::size_t b = mPos;
        while (mPos < mText.size() && !text_is_blank(mText[mPos]))
            ++mPos;
        rValue.assign(mText.substr(b, mPos - b));
        return *this;
    }

private:
    // The stream's sentry: skip blanks; failed already, or at the end, fails
    // and leaves the value alone.
    bool Sentry() {
        if (mFail)
            return false;
        while (mPos < mText.size() && text_is_blank(mText[mPos]))
            ++mPos;
        if (mPos >= mText.size()) {
            mFail = true;
            return false;
        }
        return true;
    }

    template <class T>
    TextStream& ExtractInt(T& rValue) {
        if (!Sentry())
            return *this;
        bool neg = false;
        if (mText[mPos] == '+' || mText[mPos] == '-') {
            neg = mText[mPos] == '-';
            ++mPos;
        }
        const std::size_t d = mPos;
        while (mPos < mText.size() && mText[mPos] >= '0' && mText[mPos] <= '9')
            ++mPos;
        if (mPos == d) {
            rValue = 0;
            mFail = true;
            return *this;
        }
        std::uint64_t mag = 0;
        const auto r = std::from_chars(mText.data() + d, mText.data() + mPos, mag);
        const std::uint64_t bound =
            neg ? static_cast<std::uint64_t>(-(std::numeric_limits<T>::min() + 1)) + 1
                : static_cast<std::uint64_t>(std::numeric_limits<T>::max());
        if (r.ec == std::errc::result_out_of_range || mag > bound) {
            rValue = neg ? std::numeric_limits<T>::min() : std::numeric_limits<T>::max();
            mFail = true;
            return *this;
        }
        rValue = neg ? static_cast<T>(0 - mag) : static_cast<T>(mag);
        return *this;
    }

    // num_get for an unsigned type: a sign is allowed and a negative value
    // wraps, as strtoull does; out of range stores the maximum and fails.
    template <class T>
    TextStream& ExtractUnsigned(T& rValue) {
        if (!Sentry())
            return *this;
        bool neg = false;
        if (mText[mPos] == '+' || mText[mPos] == '-') {
            neg = mText[mPos] == '-';
            ++mPos;
        }
        const std::size_t d = mPos;
        while (mPos < mText.size() && mText[mPos] >= '0' && mText[mPos] <= '9')
            ++mPos;
        if (mPos == d) {
            rValue = 0;
            mFail = true;
            return *this;
        }
        std::uint64_t mag = 0;
        const auto r = std::from_chars(mText.data() + d, mText.data() + mPos, mag);
        if (r.ec == std::errc::result_out_of_range ||
            mag > static_cast<std::uint64_t>(std::numeric_limits<T>::max())) {
            rValue = std::numeric_limits<T>::max();
            mFail = true;
            return *this;
        }
        rValue = neg ? static_cast<T>(0 - static_cast<T>(mag)) : static_cast<T>(mag);
        return *this;
    }

    std::string mOwned;
    std::string_view mText;
    std::size_t mPos = 0;
    bool mFail = false;
};

}  // namespace detail
}  // namespace meshioplusplus
