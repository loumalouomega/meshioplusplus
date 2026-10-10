// SPDX-License-Identifier: MIT
/// @file render_args.hpp
/// @brief The CLI's argument parser and the software-rendering flag vocabulary.
#ifndef MESHIOPLUSPLUS_CLI_RENDER_ARGS_HPP
#define MESHIOPLUSPLUS_CLI_RENDER_ARGS_HPP

// One parser and one flag vocabulary for every verb that draws: `snapshot`
// and the interactive `tui` (whose `:` prompt feeds its line to the same
// `cli_parse`), in the native CLI and in the Python binding alike. Extracted
// from main.cpp so the loop and the binding share it (roadmap 7.2.5, 7.2.11).

#include <optional>
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

/// Map the shared render flags onto a `RenderOptions`. Throws
/// `std::invalid_argument` / `std::runtime_error` naming a bad or conflicting flag.
RenderOptions cli_render_options(const cli_parsed& rP);

TextEncoding cli_text_encoding(const std::string& rName);
ColorDepth cli_color_depth(const std::string& rName);

}  // namespace meshioplusplus::cli

#endif  // MESHIOPLUSPLUS_CLI_RENDER_ARGS_HPP
