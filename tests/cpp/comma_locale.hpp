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
 * @file comma_locale.hpp
 * @brief A comma-decimal `LC_NUMERIC` for the tests that check number I/O
 * does not follow the process locale.
 */

// System includes
#include <clocale>
#include <string>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"

namespace mio_test {

// RAII guard: saves/restores LC_NUMERIC around a test, mandatory since every
// gtest in this binary shares one process -- a leaked comma locale would
// corrupt every later test's floating-point parsing/formatting.
class LcNumericGuard {
public:
    LcNumericGuard() {
        const char* cur = std::setlocale(LC_NUMERIC, nullptr);
        mSaved = cur ? cur : "C";
    }
    ~LcNumericGuard() { std::setlocale(LC_NUMERIC, mSaved.c_str()); }
    LcNumericGuard(const LcNumericGuard&) = delete;
    LcNumericGuard& operator=(const LcNumericGuard&) = delete;

private:
    std::string mSaved;
};

// Sets LC_NUMERIC to the first candidate whose decimal_point is a comma and
// returns its name, or returns empty if none of these locales are installed
// on this machine (`locale-gen de_DE.UTF-8` installs one on Debian/Ubuntu).
inline std::string find_comma_locale() {
    static const char* kCandidates[] = {
        "de_DE.UTF-8", "de_DE.utf8",  "fr_FR.UTF-8", "fr_FR.utf8",          "es_ES.UTF-8",
        "es_ES.utf8",  "nl_NL.UTF-8", "nl_NL.utf8",  "German_Germany.1252", "French_France.1252",
    };
    for (const char* cand : kCandidates) {
        if (std::setlocale(LC_NUMERIC, cand) != nullptr) {
            if (meshioplusplus::detail::is_c_decimal_point(std::localeconv()->decimal_point))
                continue;  // some libc builds accept the name but keep "C"
            return cand;
        }
    }
    return {};
}

}  // namespace mio_test
