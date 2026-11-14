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
 * @file detail/netcdf_util.hpp
 * @brief Core-private netCDF helpers shared by the netCDF-backed formats (Exodus,
 * UGRID NetCDF): status checks, variable reads into an `NDArray`, dimensions, name
 * arrays, text attributes and integer lists (roadmap 2.4.3.1).
 *
 * Every helper takes the format's own error prefix (`pWho`), so a format's messages
 * read as they always did ("Exodus/netCDF: get_vara: ..."). Only compiled with
 * `MESHIOPLUSPLUS_HAS_NETCDF`.
 */

#ifdef MESHIOPLUSPLUS_HAS_NETCDF

// External includes
#include <netcdf.h>

// System includes
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus::detail::nc {

/// Throw `ReadError` (`WriteError` when @p Writing) for a failed netCDF call.
void check(int Status, const char* pWhat, bool Writing = false, const char* pWho = "netCDF");

/// The netCDF type of a dtype, and the dtype of a netCDF type (`ReadError`
/// "<pWhoType>: unsupported netCDF variable type" otherwise).
nc_type nc_type_of(DType Dtype);
DType dtype_of(nc_type Type, const char* pWhoType = "netCDF");

/// A whole variable, or a start/count hyperslab, as an `NDArray` of its own dtype.
NDArray read_var(int Ncid, int Varid, const std::vector<std::size_t>& rStart,
                 const std::vector<std::size_t>& rCount, const char* pWho = "netCDF",
                 const char* pWhoType = "netCDF");

/// The lengths of a variable's dimensions.
std::vector<std::size_t> var_dims(int Ncid, int Varid, const char* pWho = "netCDF");

/// A char variable whose last dimension is the string width, as strings.
std::vector<std::string> read_names(int Ncid, int Varid, const char* pWho = "netCDF");

/// A text attribute with trailing NULs and blanks trimmed (NetCDF.jl stores the NUL).
std::string read_att_text(int Ncid, int Varid, const char* pName, const char* pWho = "netCDF");

/// A 1-D integer variable as int64.
std::vector<std::int64_t> read_ints(int Ncid, int Varid, const char* pWho = "netCDF");

}  // namespace meshioplusplus::detail::nc

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
