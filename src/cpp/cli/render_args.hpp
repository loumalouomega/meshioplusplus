// SPDX-License-Identifier: MIT
/// @file render_args.hpp
/// @brief The CLI's argument parser and the software-rendering flag vocabulary.
#ifndef MESHIOPLUSPLUS_CLI_RENDER_ARGS_HPP
#define MESHIOPLUSPLUS_CLI_RENDER_ARGS_HPP

// One parser and one flag vocabulary for every verb that draws: `snapshot`
// and the interactive `tui` (whose `:` prompt feeds its line to the same
// `cli_parse`), in the native CLI and in the Python binding alike. Extracted
// from main.cpp so the loop and the binding share it (roadmap 7.2.5, 7.2.11).

#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "meshioplusplus/operations/render.hpp"

namespace meshioplusplus::cli {

/// One option/flag spec: a canonical long name, its aliases, and whether it
/// consumes a following value.
struct cli_opt_spec {
    std::string canonical;
    std::vector<std::string> aliases;
    bool takes_value;
};

/// Parsed result of a verb's argument list.
struct cli_parsed {
    std::vector<std::string> positionals;
    std::unordered_map<std::string, std::string> values;  ///< canonical -> value (last wins)
    std::unordered_set<std::string> flags;                ///< canonical present
    /// Every occurrence of each value-option, in order. Always populated
    /// alongside `values`, which keeps its last-wins semantics so no existing
    /// verb changes behaviour; the `data` verbs read this instead because they
    /// accept repeated `--point`/`--cell`/`--field`.
    std::unordered_map<std::string, std::vector<std::string>> multi;
    /// Every value-option occurrence as (canonical name, value), in
    /// command-line order across options -- what `regions`' edit flags need.
    std::vector<std::pair<std::string, std::string>> ordered;
};

/// Parse `args` against `specs`. Throws std::runtime_error on an unknown option
/// or a value-option missing its argument. `--` stops option parsing.
cli_parsed cli_parse(const std::vector<std::string>& rArgs,
                     const std::vector<cli_opt_spec>& rSpecs);

std::string opt_value(const cli_parsed& rP, const std::string& rName,
                      const std::string& rDefault = "");
bool has_flag(const cli_parsed& rP, const std::string& rName);
/// Every occurrence of a repeatable value-option, in order (empty if absent).
const std::vector<std::string>& opt_values(const cli_parsed& rP, const std::string& rName);
/// Whether a value-option was supplied at all (distinct from "supplied empty").
bool has_opt(const cli_parsed& rP, const std::string& rName);

/// The numeric half of the `--color-by` family (`--component`, `--vmin`,
/// `--vmax`), shared by `convert` and `snapshot` so both read them alike.
void cli_color_values(const cli_parsed& rP, std::optional<int>& rComponent,
                      std::optional<double>& rVMin, std::optional<double>& rVMax);

/// The flags every software-rendering verb shares: one vocabulary.
std::vector<cli_opt_spec> render_flag_specs();

/// A comma-separated list of numbers; an empty entry is `nullopt` (so
/// `--clip 2,` leaves the high end alone).
std::vector<std::optional<double>> cli_parse_numbers(const std::string& rText, const char* pFlag);

/// `#rrggbb`, `#rrggbbaa`, or `none`/`transparent` (alpha 0).
RenderColor cli_parse_rgba(const std::string& rText, const char* pFlag);

/// A cut-away plane: six numbers `PX,PY,PZ,NX,NY,NZ` (a point and the normal
/// of the side kept), or `AXIS:OFFSET` with AXIS one of `+x -x +y -y +z -z`,
/// which keeps the side that axis points to: `+x:0.5` keeps x >= 0.5 and
/// `-x:0.5` keeps x <= 0.5.
RenderCutaway cli_parse_cutaway(const std::string& rText);

/// Map the shared render flags onto a `RenderOptions`. Throws
/// `std::invalid_argument` / `std::runtime_error` naming a bad or conflicting flag.
RenderOptions cli_render_options(const cli_parsed& rP);

/**
 * The render flags as a mutable set: what `snapshot` takes, held by name so a
 * viewer can change one and re-derive the options. `FromOptions` is the inverse
 * of `cli_render_options` for everything a flag can say (only what differs from
 * the defaults is kept, in the flags' own spellings), so a `:cmap turbo`
 * typed in the viewer, a saved session and a command line are one vocabulary.
 * The frame size is not a render flag and is not kept.
 */
class RenderFlags {
public:
    RenderFlags() = default;
    static RenderFlags FromOptions(const RenderOptions& rOptions);
    /// From argv-style tokens (`--cmap turbo --axes`); throws like `cli_parse`.
    static RenderFlags FromTokens(const std::vector<std::string>& rTokens);

    /// Whether a flag (by canonical name, no dashes) is present.
    bool Has(const std::string& rName) const;
    /// A value flag's last value, or `rDefault`.
    std::string Value(const std::string& rName, const std::string& rDefault = "") const;
    /// Every value of a value flag, in order (a repeatable one has several).
    const std::vector<std::string>& Values(const std::string& rName) const;
    /// Replace a value flag (all earlier values go).
    void Set(const std::string& rName, const std::string& rValue);
    /// Append to a repeatable value flag (`cutaway`).
    void Add(const std::string& rName, const std::string& rValue);
    /// Turn a boolean flag on or off.
    void SetFlag(const std::string& rName, bool On);
    void Unset(const std::string& rName);
    /// Whether the spec table says the flag takes a value.
    static bool TakesValue(const std::string& rName);
    /// Whether `rName` is a render flag at all.
    static bool Known(const std::string& rName);

    /// The options these flags describe; throws like `cli_render_options` (a
    /// range flag without a mapped field, a bad value, ...).
    RenderOptions ToOptions() const;
    /// argv-style tokens, sorted by name, that `FromTokens` reads back.
    std::vector<std::string> Tokens() const;

private:
    std::map<std::string, std::vector<std::string>> mValues;
    std::set<std::string> mFlags;
};

TextEncoding cli_text_encoding(const std::string& rName);
ColorDepth cli_color_depth(const std::string& rName);

}  // namespace meshioplusplus::cli

#endif  // MESHIOPLUSPLUS_CLI_RENDER_ARGS_HPP
