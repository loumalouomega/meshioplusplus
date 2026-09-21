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

// System includes
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/detail/keyword_card.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

std::string kwc_strip(std::string_view Text) {
    std::size_t b = 0;
    std::size_t e = Text.size();
    while (b < e && std::isspace(static_cast<unsigned char>(Text[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(Text[e - 1])))
        --e;
    return std::string(Text.substr(b, e - b));
}

bool kwc_is_digit(char c) {
    return c >= '0' && c <= '9';
}

// "1.5-3" -> "1.5e-3": a Fortran exponent written without its letter.
std::string kwc_add_exponent_letter(const std::string& rText) {
    for (char c : rText)
        if (c == 'e' || c == 'E')
            return rText;
    std::size_t k = std::string::npos;
    for (std::size_t i = rText.size(); i-- > 1;) {
        if (rText[i] == '+' || rText[i] == '-') {
            k = i;
            break;
        }
    }
    if (k == std::string::npos || k + 1 >= rText.size())
        return rText;
    for (std::size_t i = k + 1; i < rText.size(); ++i)
        if (!kwc_is_digit(rText[i]))
            return rText;
    std::size_t b = (rText[0] == '+' || rText[0] == '-') ? 1 : 0;
    int digits = 0;
    int dots = 0;
    for (std::size_t i = b; i < k; ++i) {
        if (kwc_is_digit(rText[i]))
            ++digits;
        else if (rText[i] == '.')
            ++dots;
        else
            return rText;
    }
    if (digits == 0 || dots > 1)
        return rText;
    return rText.substr(0, k) + "e" + rText.substr(k);
}

}  // namespace

int card_field_width(const CardField& rField, CardMode Mode) {
    if (Mode == CardMode::Long)
        return 20;
    if (Mode == CardMode::I10 && rField.mKind == 'i' && rField.mWidth == 8)
        return 10;
    return rField.mWidth;
}

std::vector<std::string> split_card(std::string_view Line, const std::vector<CardField>& rLayout,
                                    CardMode Mode) {
    std::vector<std::string> out;
    if (Line.find(',') != std::string_view::npos) {
        std::size_t start = 0;
        while (true) {
            const std::size_t comma = Line.find(',', start);
            if (comma == std::string_view::npos) {
                out.push_back(kwc_strip(Line.substr(start)));
                break;
            }
            out.push_back(kwc_strip(Line.substr(start, comma - start)));
            start = comma + 1;
        }
        while (out.size() < rLayout.size())
            out.emplace_back();
        return out;
    }
    std::size_t pos = 0;
    out.reserve(rLayout.size());
    for (const CardField& field : rLayout) {
        const std::size_t w = static_cast<std::size_t>(card_field_width(field, Mode));
        if (pos >= Line.size())
            out.emplace_back();
        else
            out.push_back(kwc_strip(Line.substr(pos, w)));
        pos += w;
    }
    return out;
}

std::int64_t card_to_int(const std::string& rText, const std::string& rWhere) {
    if (rText.empty())
        return 0;
    errno = 0;
    char* end = nullptr;
    const long long v = std::strtoll(rText.c_str(), &end, 10);
    if (end == rText.c_str() || *end != '\0' || errno == ERANGE)
        throw ReadError("LS-DYNA: invalid integer field '" + rText + "'" + rWhere);
    return static_cast<std::int64_t>(v);
}

double card_to_real(const std::string& rText, const std::string& rWhere) {
    if (rText.empty())
        return 0.0;
    std::string s = rText;
    for (char& c : s) {
        if (c == 'D')
            c = 'E';
        else if (c == 'd')
            c = 'e';
    }
    s = kwc_add_exponent_letter(s);
    const char* end = nullptr;
    const double v = parse_double(s.c_str(), end);
    if (end == s.c_str() || *end != '\0')
        throw ReadError("LS-DYNA: invalid real field '" + rText + "'" + rWhere);
    return v;
}

std::string format_real16(double Value) {
    if (Value == 0.0)
        return "0.0";
    const int neg = Value < 0.0 ? 1 : 0;
    const int e3 = (std::fabs(Value) >= 1e100 || std::fabs(Value) < 1e-99) ? 1 : 0;
    const int pmax = 10 - neg - e3;
    char buf[64];
    for (int p = 1; p <= pmax; ++p) {
        snprintf_c(buf, sizeof(buf), "%.*e", p, Value);
        const char* end = nullptr;
        if (parse_double(buf, end) == Value)
            return std::string(buf);
    }
    return std::string(buf);
}

}  // namespace detail
}  // namespace meshioplusplus
