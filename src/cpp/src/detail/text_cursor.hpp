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
 * accepted, an out-of-range value saturated. Roadmap §4, "A shared tokenizer
 * and number path".
 */

// System includes
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <system_error>
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

}  // namespace detail
}  // namespace meshioplusplus
