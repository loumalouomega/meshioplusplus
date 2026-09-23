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
 * @file keyword_card.hpp
 * @brief Fixed-width card tokenizer for solver keyword decks (LS-DYNA first).
 *
 * A card is one line whose fields sit in fixed columns. Four layouts exist and are
 * told apart per card:
 *
 *  - **Standard**: the widths a card declares (`I8`, `E16`, `I10` ...);
 *  - **Long** (`*KEYWORD LONG=Y`, or a `+` keyword suffix): every field, integer and
 *    real, is 20 columns wide;
 *  - **I10** (`*KEYWORD I10=Y`): only the 8-column integer fields widen to 10, a
 *    different mechanism from long;
 *  - **Free**: a line that contains a comma is split on commas, whatever the mode.
 *
 * The Python twin is `lsdyna/_cards.py`; keep the two in step.
 */

// System includes
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {
namespace detail {

/// How a card's fixed columns are laid out.
enum class CardMode { Standard, Long, I10 };

/// One field of a card layout: `mKind` is `'i'` (integer), `'r'` (real), `'a'`
/// (characters) or `'x'` (skipped columns, from a Fortran format's `nX`).
struct CardField {
    char mKind;
    int mWidth;
};

/// The column width of one field under `Mode`.
MESHIOPLUSPLUS_API int card_field_width(const CardField& rField, CardMode Mode);

/**
 * @brief The stripped text of each field of `Line`; a missing field is `""`.
 *
 * A line containing a comma is split on commas and returns every entry, so a card
 * with a variable number of fields is not truncated to the layout.
 */
MESHIOPLUSPLUS_API std::vector<std::string> split_card(std::string_view Line,
                                                       const std::vector<CardField>& rLayout,
                                                       CardMode Mode);

/// An integer field; blank is 0. @throws ReadError on garbage (`rWhere` is appended).
MESHIOPLUSPLUS_API std::int64_t card_to_int(const std::string& rText, const std::string& rWhere);

/**
 * @brief A real field; blank is 0. Accepts `D` exponents and `1.5-3` (no letter).
 * @throws ReadError on garbage (`rWhere` is appended).
 */
MESHIOPLUSPLUS_API double card_to_real(const std::string& rText, const std::string& rWhere);

/// `card_to_int` whose error names `rFormat` (e.g. `"Nastran"`) instead of LS-DYNA.
MESHIOPLUSPLUS_API std::int64_t card_to_int(const std::string& rText, const std::string& rWhere,
                                            const std::string& rFormat);

/// `card_to_real` whose error names `rFormat` (e.g. `"Nastran"`) instead of LS-DYNA.
MESHIOPLUSPLUS_API double card_to_real(const std::string& rText, const std::string& rWhere,
                                       const std::string& rFormat);

/**
 * @brief The fields of a Fortran edit-descriptor list, such as the format line of
 * an ANSYS `.cdb` block: `"(1i7,2i9,6e21.13e3)"` gives one `i` field of width 7,
 * two of width 9 and six `r` fields of width 21.
 *
 * `Iw` is an integer field; `Ew.d[Ee]`, `ESw.d`, `ENw.d`, `Dw.d`, `Fw.d` and `Gw.d`
 * real fields; `Aw` a character field; `nX` skipped columns. Repeat counts and
 * parenthesised groups (`2(i8,e16.9)`) expand; a scale factor (`1P`) is ignored.
 * Case does not matter, blanks are ignored.
 * @throws ReadError naming the format for anything else.
 */
MESHIOPLUSPLUS_API std::vector<CardField> parse_fortran_format(std::string_view Format);

/**
 * @brief The stripped text of each field `rFields` lays out on `Line`, in fixed
 * columns (values may touch, as in `1.0E+00-2.0E+00`). Slicing stops at the end
 * of the line, so a short line gives fewer fields; `x` fields are not returned.
 */
MESHIOPLUSPLUS_API std::vector<std::string> split_fixed(std::string_view Line,
                                                        const std::vector<CardField>& rFields);

/**
 * @brief `Value` in at most 16 columns: the shortest scientific string that
 * round-trips, else as many digits as fit.
 */
MESHIOPLUSPLUS_API std::string format_real16(double Value);

}  // namespace detail
}  // namespace meshioplusplus
