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
// Core-private netCDF helpers shared by the netCDF-backed formats. See
// detail/netcdf_util.hpp.

#ifdef MESHIOPLUSPLUS_HAS_NETCDF

// System includes
#include <cstring>
#include <string>
#include <vector>

// Project includes
#include "netcdf_util.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "dtype_table.hpp"

namespace meshioplusplus::detail::nc {

void check(int Status, const char* pWhat, bool Writing, const char* pWho) {
    if (Status != NC_NOERR) {
        std::string msg = std::string(pWho) + ": " + pWhat + ": " + nc_strerror(Status);
        if (Writing)
            throw WriteError(msg);
        throw ReadError(msg);
    }
}

nc_type nc_type_of(DType Dtype) {
    return detail::DTypeUtils::Pick<nc_type>(
        Dtype, {NC_FLOAT, NC_DOUBLE, NC_BYTE, NC_SHORT, NC_INT, NC_INT64, NC_UBYTE, NC_USHORT,
                NC_UINT, NC_UINT64});
}

DType dtype_of(nc_type Type, const char* pWhoType) {
    switch (Type) {
        case NC_FLOAT:
            return DType::Float32;
        case NC_DOUBLE:
            return DType::Float64;
        case NC_BYTE:
            return DType::Int8;
        case NC_SHORT:
            return DType::Int16;
        case NC_INT:
            return DType::Int32;
        case NC_INT64:
            return DType::Int64;
        case NC_UBYTE:
            return DType::UInt8;
        case NC_USHORT:
            return DType::UInt16;
        case NC_UINT:
            return DType::UInt32;
        case NC_UINT64:
            return DType::UInt64;
        default:
            throw ReadError(std::string(pWhoType) + ": unsupported netCDF variable type");
    }
}

NDArray read_var(int Ncid, int Varid, const std::vector<std::size_t>& rStart,
                 const std::vector<std::size_t>& rCount, const char* pWho, const char* pWhoType) {
    nc_type t;
    check(nc_inq_vartype(Ncid, Varid, &t), "inq_vartype", false, pWho);
    const DType dt = dtype_of(t, pWhoType);
    std::vector<std::size_t> shape(rCount.begin(), rCount.end());
    NDArray out(dt, shape);
    if (out.Size() > 0)
        check(nc_get_vara(Ncid, Varid, rStart.data(), rCount.data(), out.Data()), "get_vara", false,
              pWho);
    return out;
}

std::vector<std::size_t> var_dims(int Ncid, int Varid, const char* pWho) {
    int ndims;
    check(nc_inq_varndims(Ncid, Varid, &ndims), "inq_varndims", false, pWho);
    std::vector<int> dimids(static_cast<std::size_t>(ndims));
    check(nc_inq_vardimid(Ncid, Varid, dimids.data()), "inq_vardimid", false, pWho);
    std::vector<std::size_t> out;
    for (int d : dimids) {
        std::size_t len;
        check(nc_inq_dimlen(Ncid, d, &len), "inq_dimlen", false, pWho);
        out.push_back(len);
    }
    return out;
}

// Covers the 2-D `(n, len_string)` name arrays and the 3-D `(num_qa_rec, four,
// len_string)` QA records, which flatten to 4 strings per record.
std::vector<std::string> read_names(int Ncid, int Varid, const char* pWho) {
    const std::vector<std::size_t> dims = var_dims(Ncid, Varid, pWho);
    if (dims.empty())
        return {};
    const std::size_t w = dims.back();
    std::size_t n = 1;
    for (std::size_t k = 0; k + 1 < dims.size(); ++k)
        n *= dims[k];
    std::vector<char> buf(n * w, '\0');
    if (n * w > 0)
        check(nc_get_var_text(Ncid, Varid, buf.data()), "get names", false, pWho);
    std::vector<std::string> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        out.emplace_back(buf.data() + i * w, strnlen(buf.data() + i * w, w));
    return out;
}

// NetCDF.jl stores the C string's terminating NUL inside the attribute (an Exodus
// SPHERE block's `elem_type` arrives as "SPHERE\0"); fixed-width writers pad with
// spaces. Both are trimmed.
std::string read_att_text(int Ncid, int Varid, const char* pName, const char* pWho) {
    std::size_t attlen = 0;
    check(nc_inq_attlen(Ncid, Varid, pName, &attlen), pName, false, pWho);
    std::string out(attlen, '\0');
    if (attlen > 0)
        check(nc_get_att_text(Ncid, Varid, pName, out.data()), pName, false, pWho);
    const std::size_t last = out.find_last_not_of(std::string("\0 \t\r\n", 5));
    out.erase(last == std::string::npos ? 0 : last + 1);
    return out;
}

std::vector<std::int64_t> read_ints(int Ncid, int Varid, const char* pWho) {
    const std::vector<std::size_t> dims = var_dims(Ncid, Varid, pWho);
    const std::size_t n = dims.empty() ? 0 : dims[0];
    std::vector<std::int64_t> out(n, 0);
    if (n > 0)
        check(nc_get_var_longlong(Ncid, Varid, reinterpret_cast<long long*>(out.data())),
              "get ints", false, pWho);
    return out;
}

}  // namespace meshioplusplus::detail::nc

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
