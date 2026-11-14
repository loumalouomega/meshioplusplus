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
 * @file operations/detect.hpp
 * @brief A ranked format detection report (roadmap 2.4.9).
 *
 * `sniff_format` answers with one name or `""`, so a caller cannot tell a sure
 * answer from a guess, nor that `.msh`, `.mesh`, `.inp`, `.dat`, `.ele` and `.node`
 * name several formats. `detect_format` returns every candidate with how sure it is:
 *
 * - `certain`: a content signature matched (`sniff_format`);
 * - `likely`: the extension names the format and the file's metadata reads as it;
 * - `guess`: the extension alone names it.
 *
 * Candidates are ranked certain, likely, guess, each in the extension's own order.
 * Two `likely` entries and no `certain` one is an ambiguous file. The twin is the
 * Python `meshioplusplus.detect`.
 */

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"

namespace meshioplusplus {

/** @brief One candidate of a detection report. */
struct DetectCandidate {
    std::string mFormat;      ///< a meshio++ format name
    std::string mConfidence;  ///< "certain", "likely" or "guess"
    std::string mReason;      ///< why, in words
};

/**
 * @brief Every format @p rPath could be, ranked (empty when nothing matches).
 * @param rPath an existing file or directory
 */
MESHIOPLUSPLUS_API std::vector<DetectCandidate> detect_format(const std::string& rPath);

/** @brief The report as text: one `format<TAB>confidence<TAB>reason` line per candidate. */
MESHIOPLUSPLUS_API std::string detect_format_text(const std::string& rPath);

}  // namespace meshioplusplus
