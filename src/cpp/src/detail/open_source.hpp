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
 * @file detail/open_source.hpp
 * @brief A reader's whole file through `FileSource`, failing with the
 * reader's own message.
 *
 * A **core-private** header (the `slot_runs.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 *
 * Readers that slurped their file through `std::istreambuf_iterator` -- one
 * character at a time into a heap string -- read it through `FileSource`
 * instead: one bulk read, or a mapping above `MESHIOPLUSPLUS_MMAP_THRESHOLD`.
 * `FileSource` reports an unreadable file as "Could not open file"; this keeps
 * the message each reader gave before. Roadmap §4, "Memory and allocation".
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

/// The whole file at @p rPath, or `ReadError(rMessage)` when it cannot be read.
inline FileSource open_source(const std::string& rPath, const std::string& rMessage) {
    try {
        return FileSource(rPath);
    } catch (const ReadError&) {
        throw ReadError(rMessage);
    }
}

}  // namespace detail
}  // namespace meshioplusplus
