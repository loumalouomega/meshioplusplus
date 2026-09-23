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
    return card_to_int(rText, rWhere, "LS-DYNA");
}

std::int64_t card_to_int(const std::string& rText, const std::string& rWhere,
                         const std::string& rFormat) {
    if (rText.empty())
        return 0;
    errno = 0;
    char* end = nullptr;
    const long long v = std::strtoll(rText.c_str(), &end, 10);
    if (end == rText.c_str() || *end != '\0' || errno == ERANGE)
        throw ReadError(rFormat + ": invalid integer field '" + rText + "'" + rWhere);
    return static_cast<std::int64_t>(v);
}

double card_to_real(const std::string& rText, const std::string& rWhere) {
    return card_to_real(rText, rWhere, "LS-DYNA");
}

double card_to_real(const std::string& rText, const std::string& rWhere,
                    const std::string& rFormat) {
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
        throw ReadError(rFormat + ": invalid real field '" + rText + "'" + rWhere);
    return v;
}

namespace {

// Recursive-descent parse of a Fortran edit-descriptor list, appending fields.
struct KcFormatParser {
    std::string_view mText;
    std::size_t mPos = 0;

    [[noreturn]] void Fail() const {
        throw ReadError("cannot parse Fortran format '" + std::string(mText) + "'");
    }
    void Skip() {
        while (mPos < mText.size() && (mText[mPos] == ' ' || mText[mPos] == '\t'))
            ++mPos;
    }
    bool Digit() const { return mPos < mText.size() && mText[mPos] >= '0' && mText[mPos] <= '9'; }
    int Number() {
        Skip();
        if (!Digit())
            return -1;
        int v = 0;
        while (Digit()) {
            v = v * 10 + (mText[mPos++] - '0');
            if (v > 100000)
                Fail();
        }
        return v;
    }
    char Letter() {
        Skip();
        if (mPos >= mText.size())
            Fail();
        const char c = mText[mPos];
        return static_cast<char>(c >= 'A' && c <= 'Z' ? c - 'A' + 'a' : c);
    }

    // A comma-separated list up to `Close` (')' or end of text).
    void List(std::vector<CardField>& rOut, bool Nested) {
        while (true) {
            Skip();
            if (mPos >= mText.size()) {
                if (Nested)
                    Fail();
                return;
            }
            if (mText[mPos] == ')') {
                if (!Nested)
                    Fail();
                ++mPos;
                return;
            }
            Item(rOut);
            Skip();
            if (mPos < mText.size() && mText[mPos] == ',')
                ++mPos;
        }
    }

    void Item(std::vector<CardField>& rOut) {
        const int repeat = Number();
        const int count = repeat < 0 ? 1 : repeat;
        Skip();
        if (mPos < mText.size() && mText[mPos] == '(') {
            ++mPos;
            std::vector<CardField> group;
            List(group, true);
            for (int k = 0; k < count; ++k)
                rOut.insert(rOut.end(), group.begin(), group.end());
            return;
        }
        const char c = Letter();
        ++mPos;
        if (c == 'p') {  // scale factor: `1P`
            if (repeat < 0)
                Fail();
            return;
        }
        if (c == 'x') {  // `nX`: skipped columns
            rOut.push_back(CardField{'x', count});
            return;
        }
        char kind = 0;
        if (c == 'i')
            kind = 'i';
        else if (c == 'a')
            kind = 'a';
        else if (c == 'e' || c == 'd' || c == 'f' || c == 'g')
            kind = 'r';
        else
            Fail();
        // ES / EN variants of the E descriptor.
        if (c == 'e' && mPos < mText.size()) {
            const char n = static_cast<char>(mText[mPos] | 0x20);
            if (n == 's' || n == 'n')
                ++mPos;
        }
        const int width = Number();
        if (width <= 0)
            Fail();
        Skip();
        if (mPos < mText.size() && mText[mPos] == '.') {  // .d
            ++mPos;
            if (Number() < 0)
                Fail();
            Skip();
            if (mPos < mText.size() && (mText[mPos] | 0x20) == 'e' && kind == 'r') {  // Ee
                ++mPos;
                if (Number() < 0)
                    Fail();
            }
        }
        for (int k = 0; k < count; ++k)
            rOut.push_back(CardField{kind, width});
    }
};

}  // namespace

std::vector<CardField> parse_fortran_format(std::string_view Format) {
    KcFormatParser parser{Format};
    parser.Skip();
    std::vector<CardField> out;
    if (parser.mPos < Format.size() && Format[parser.mPos] == '(') {
        ++parser.mPos;
        parser.List(out, true);
        parser.Skip();
        if (parser.mPos != Format.size())
            parser.Fail();
    } else {
        parser.List(out, false);
    }
    if (out.empty())
        parser.Fail();
    return out;
}

std::vector<std::string> split_fixed(std::string_view Line, const std::vector<CardField>& rFields) {
    while (!Line.empty() && (Line.back() == '\r' || Line.back() == '\n'))
        Line.remove_suffix(1);
    std::vector<std::string> out;
    std::size_t col = 0;
    for (const CardField& f : rFields) {
        if (col >= Line.size())
            break;
        const std::size_t width = static_cast<std::size_t>(f.mWidth);
        if (f.mKind != 'x') {
            std::string_view text = Line.substr(col, width);
            const std::size_t a = text.find_first_not_of(" \t");
            const std::size_t b = text.find_last_not_of(" \t");
            out.emplace_back(a == std::string_view::npos ? std::string_view()
                                                         : text.substr(a, b - a + 1));
        }
        col += width;
    }
    return out;
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
