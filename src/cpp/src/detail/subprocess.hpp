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
 * @file detail/subprocess.hpp
 * @brief Run an installed program as a separate process (roadmap 7.7.2, the
 * external tetrahedral meshers): find it on `PATH`, run it with an argument
 * vector (no shell), its output sent to a log file, and wait for it.
 *
 * Core-internal. Compiled out under Emscripten, where `subprocess_available()`
 * is false and running throws `Unsupported`.
 */

// System includes
#include <string>
#include <vector>

namespace meshioplusplus::detail {

/// Whether this build can start processes.
bool subprocess_available() noexcept;

/// The full path of the program @p rName: the environment variable
/// @p pOverride when it is set (and names a file), else the first match on
/// `PATH` (with `.exe` on Windows); empty when there is none.
std::string find_executable(const std::string& rName, const char* pOverride = nullptr);

/// Runs @p rArgv (`rArgv[0]` the program's full path) in @p rWorkDir, standard
/// output and error appended to @p rLogPath, and waits. Returns the exit code
/// (-1 when the process ended by a signal). Throws `std::runtime_error` when it
/// cannot start.
int run_process(const std::vector<std::string>& rArgv, const std::string& rWorkDir,
                const std::string& rLogPath);

/// A fresh, empty directory under the system temporary directory, named after
/// @p rTag; the caller removes it.
std::string make_temp_dir(const std::string& rTag);

}  // namespace meshioplusplus::detail
