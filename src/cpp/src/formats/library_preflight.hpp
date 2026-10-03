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
 * @file library_preflight.hpp
 * @brief Refuse obviously-corrupt containers before the library sees them.
 *
 * A **core-private** header beside the formats (the `vtk_preflight.hpp`
 * precedent): it adds nothing to the installed headers or the ABI.
 *
 * HDF5's own container parser aborts on some corrupted inputs (dense link
 * tables, huge object headers) instead of returning an error, and netCDF-4
 * inherits the dense-link abort through `nc_open`. The fuzzer's fork isolation
 * (tests/fuzz/fuzz_read.cpp) keeps such an abort from corrupting the campaign,
 * but production readers should still refuse what they can recognize without
 * the library: a missing file, an empty one, a wrong magic, or an HDF5
 * superblock version the library never wrote. The checks stay loose on
 * purpose: anything they cannot decide is left for the library and the
 * reader's own validation, whose `ReadError` remains the verdict.
 */

#include <cstdint>
#include <cstdio>
#include <string>

// Project includes
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

/// HDF5's eight-byte magic starting every HDF5 container (including netCDF-4).
inline constexpr char kLibraryPreflightHdf5Magic[] = "\x89HDF\r\n\x1a\n";

/// netCDF classic magic (`CDF\x01`, `CDF\x02`, `CDF\x05` for 64-bit variants).
inline constexpr char kLibraryPreflightCdfMagic[] = "CDF";

/**
 * @brief Throws `ReadError` when `rPath` cannot be an HDF5 container.
 * @param rPath File to inspect.
 * @param pFormat Format name for the message.
 *
 * Checks the file opens, holds at least the eight magic bytes, starts with
 * the HDF5 magic, and names a superblock version 0, 2 or 3 (the versions
 * HDF5 1.8 through 1.14 write). Anything else is left for the library.
 */
inline void library_preflight_hdf5(const std::string& rPath, const char* pFormat) {
    auto in = make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError(std::string("meshio++: ") + pFormat + ": cannot open '" + rPath + "'");
    char head[16] = {0};
    in.read(head, sizeof head);
    const std::streamsize got = in.gcount();
    if (got < 8)
        throw ReadError(std::string("meshio++: ") + pFormat + ": '" + rPath +
                        "' is too short for an HDF5 container");
    for (int i = 0; i < 8; ++i)
        if (head[i] != kLibraryPreflightHdf5Magic[i])
            throw ReadError(std::string("meshio++: ") + pFormat + ": '" + rPath +
                            "' is not an HDF5 container");
    if (got >= 9) {
        const unsigned version = static_cast<unsigned char>(head[8]);
        if (version != 0 && version != 2 && version != 3)
            throw ReadError(std::string("meshio++: ") + pFormat + ": '" + rPath +
                            "' names HDF5 superblock version " + std::to_string(version));
    }
}

/**
 * @brief Throws `ReadError` when `rPath` cannot be a netCDF/Exodus container.
 * @param rPath File to inspect.
 *
 * Accepts either the netCDF classic magic (`CDF\x01`/`\x02`/`\x05`) or the
 * HDF5 magic (netCDF-4, which is an HDF5 container and gets the HDF5
 * superblock check too). Anything else is refused before `nc_open` runs.
 */
inline void library_preflight_netcdf(const std::string& rPath) {
    auto in = make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError(std::string("meshio++: exodus: cannot open '") + rPath + "'");
    char head[16] = {0};
    in.read(head, sizeof head);
    const std::streamsize got = in.gcount();
    if (got < 4)
        throw ReadError(std::string("meshio++: exodus: '") + rPath +
                        "' is too short for a netCDF container");
    const bool is_cdf = head[0] == 'C' && head[1] == 'D' && head[2] == 'F';
    bool is_hdf5 = got >= 8;
    for (int i = 0; is_hdf5 && i < 8; ++i)
        is_hdf5 = head[i] == kLibraryPreflightHdf5Magic[i];
    if (!is_cdf && !is_hdf5)
        throw ReadError(std::string("meshio++: exodus: '") + rPath +
                        "' is neither netCDF classic nor HDF5 (netCDF-4)");
    if (is_hdf5)
        library_preflight_hdf5(rPath, "exodus");
}

}  // namespace detail
}  // namespace meshioplusplus
