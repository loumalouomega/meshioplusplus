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

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <cstring>

// Project includes
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace h5 {

Hid open_file_read(const std::string& rPath) {
    Hid f(H5Fopen(rPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    if (!f.Valid())
        throw ReadError("HDF5: could not open file " + rPath);
    return f;
}

Hid create_file(const std::string& rPath) {
    Hid f(H5Fcreate(rPath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
    if (!f.Valid())
        throw WriteError("HDF5: could not create file " + rPath);
    return f;
}

Hid open_file_rw(const std::string& rPath) {
    Hid f(H5Fopen(rPath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
    if (!f.Valid())
        throw WriteError("HDF5: could not open file for writing " + rPath);
    return f;
}

void flush_file(Hid& rFile) {
    if (!rFile.Valid())
        return;
    if (H5Fflush(rFile, H5F_SCOPE_GLOBAL) < 0)
        throw WriteError("HDF5: could not flush file");
}

std::vector<std::string> link_names(hid_t loc) {
    H5G_info_t info{};
    std::vector<std::string> out;
    if (H5Gget_info(loc, &info) < 0)
        return out;
    for (hsize_t i = 0; i < info.nlinks; ++i) {
        const ssize_t len =
            H5Lget_name_by_idx(loc, ".", H5_INDEX_NAME, H5_ITER_INC, i, nullptr, 0, H5P_DEFAULT);
        if (len <= 0)
            continue;
        std::string name(static_cast<std::size_t>(len), '\0');
        H5Lget_name_by_idx(loc, ".", H5_INDEX_NAME, H5_ITER_INC, i, name.data(),
                           static_cast<std::size_t>(len) + 1, H5P_DEFAULT);
        out.push_back(std::move(name));
    }
    return out;
}

bool exists(hid_t loc, const std::string& rName) {
    return H5Lexists(loc, rName.c_str(), H5P_DEFAULT) > 0;
}

Hid open_group(hid_t loc, const std::string& rName) {
    Hid g(H5Gopen2(loc, rName.c_str(), H5P_DEFAULT), H5Gclose);
    if (!g.Valid())
        throw ReadError("HDF5: missing group '" + rName + "'");
    return g;
}

Hid create_group(hid_t loc, const std::string& rName) {
    Hid g(H5Gcreate2(loc, rName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
    if (!g.Valid())
        throw WriteError("HDF5: could not create group '" + rName + "'");
    return g;
}

hid_t native_type(DType dt) {
    switch (dt) {
        case DType::Float32:
            return H5T_NATIVE_FLOAT;
        case DType::Float64:
            return H5T_NATIVE_DOUBLE;
        case DType::Int8:
            return H5T_NATIVE_INT8;
        case DType::Int16:
            return H5T_NATIVE_INT16;
        case DType::Int32:
            return H5T_NATIVE_INT32;
        case DType::Int64:
            return H5T_NATIVE_INT64;
        case DType::UInt8:
            return H5T_NATIVE_UINT8;
        case DType::UInt16:
            return H5T_NATIVE_UINT16;
        case DType::UInt32:
            return H5T_NATIVE_UINT32;
        case DType::UInt64:
            return H5T_NATIVE_UINT64;
    }
    return H5T_NATIVE_DOUBLE;
}

hid_t file_type(DType dt) {
    switch (dt) {
        case DType::Float32:
            return H5T_IEEE_F32LE;
        case DType::Float64:
            return H5T_IEEE_F64LE;
        case DType::Int8:
            return H5T_STD_I8LE;
        case DType::Int16:
            return H5T_STD_I16LE;
        case DType::Int32:
            return H5T_STD_I32LE;
        case DType::Int64:
            return H5T_STD_I64LE;
        case DType::UInt8:
            return H5T_STD_U8LE;
        case DType::UInt16:
            return H5T_STD_U16LE;
        case DType::UInt32:
            return H5T_STD_U32LE;
        case DType::UInt64:
            return H5T_STD_U64LE;
    }
    return H5T_IEEE_F64LE;
}

DType dtype_from_h5(hid_t type_id) {
    H5T_class_t cls = H5Tget_class(type_id);
    std::size_t sz = H5Tget_size(type_id);
    if (cls == H5T_FLOAT)
        return sz == 4 ? DType::Float32 : DType::Float64;
    if (cls == H5T_INTEGER) {
        bool is_signed = H5Tget_sign(type_id) != H5T_SGN_NONE;
        switch (sz) {
            case 1:
                return is_signed ? DType::Int8 : DType::UInt8;
            case 2:
                return is_signed ? DType::Int16 : DType::UInt16;
            case 4:
                return is_signed ? DType::Int32 : DType::UInt32;
            default:
                return is_signed ? DType::Int64 : DType::UInt64;
        }
    }
    throw ReadError("HDF5: unsupported datatype class");
}

NDArray read_dataset(hid_t loc, const std::string& rName) {
    Hid d(H5Dopen2(loc, rName.c_str(), H5P_DEFAULT), H5Dclose);
    if (!d.Valid())
        throw ReadError("HDF5: missing dataset '" + rName + "'");
    Hid space(H5Dget_space(d), H5Sclose);
    int ndim = H5Sget_simple_extent_ndims(space);
    std::vector<hsize_t> hdims(ndim > 0 ? ndim : 0);
    if (ndim > 0)
        H5Sget_simple_extent_dims(space, hdims.data(), nullptr);
    Hid dt(H5Dget_type(d), H5Tclose);

    std::vector<std::size_t> shape(hdims.begin(), hdims.end());
    if (shape.empty())
        shape.push_back(1);  // scalar -> length-1

    DType mdt;
    if (H5Tget_class(dt) == H5T_ARRAY) {
        Hid base(H5Tget_super(dt), H5Tclose);
        mdt = dtype_from_h5(base);
        int arank = H5Tget_array_ndims(dt);
        std::vector<hsize_t> adims(arank > 0 ? arank : 0);
        if (arank > 0)
            H5Tget_array_dims2(dt, adims.data());
        for (hsize_t ad : adims)
            shape.push_back(static_cast<std::size_t>(ad));
    } else {
        mdt = dtype_from_h5(dt);
    }

    NDArray out(mdt, shape);
    if (out.Size() > 0) {
        // For ARRAY-typed datasets the memory type must be the matching array
        // type; for scalar types the plain native type suffices.
        if (H5Tget_class(dt) == H5T_ARRAY) {
            int arank = H5Tget_array_ndims(dt);
            std::vector<hsize_t> adims(arank > 0 ? arank : 0);
            if (arank > 0)
                H5Tget_array_dims2(dt, adims.data());
            Hid mem(H5Tarray_create2(native_type(mdt), arank, adims.data()), H5Tclose);
            if (H5Dread(d, mem, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.Data()) < 0)
                throw ReadError("HDF5: failed reading dataset '" + rName + "'");
        } else if (H5Dread(d, native_type(mdt), H5S_ALL, H5S_ALL, H5P_DEFAULT, out.Data()) < 0) {
            throw ReadError("HDF5: failed reading dataset '" + rName + "'");
        }
    }
    return out;
}

void write_dataset(hid_t loc, const std::string& rName, const NDArray& rArr, int gzip_level) {
    std::vector<hsize_t> hdims(rArr.Shape().begin(), rArr.Shape().end());
    if (hdims.empty())
        hdims.push_back(0);
    Hid space(H5Screate_simple(static_cast<int>(hdims.size()), hdims.data(), nullptr), H5Sclose);

    Hid dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose);
    if (gzip_level >= 0 && rArr.Size() > 0) {
        H5Pset_chunk(dcpl, static_cast<int>(hdims.size()), hdims.data());
        H5Pset_deflate(dcpl, static_cast<unsigned>(gzip_level));
    }

    Hid d(H5Dcreate2(loc, rName.c_str(), file_type(rArr.Dtype()), space, H5P_DEFAULT, dcpl,
                     H5P_DEFAULT),
          H5Dclose);
    if (!d.Valid())
        throw WriteError("HDF5: could not create dataset '" + rName + "'");
    if (rArr.Size() > 0) {
        if (H5Dwrite(d, native_type(rArr.Dtype()), H5S_ALL, H5S_ALL, H5P_DEFAULT, rArr.Data()) < 0)
            throw WriteError("HDF5: failed writing dataset '" + rName + "'");
    }
}

bool has_attr(hid_t loc, const std::string& rName) {
    return H5Aexists(loc, rName.c_str()) > 0;
}

std::int64_t read_attr_int(hid_t loc, const std::string& rName) {
    Hid a(H5Aopen(loc, rName.c_str(), H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw ReadError("HDF5: missing attribute '" + rName + "'");
    std::int64_t v = 0;
    if (H5Aread(a, H5T_NATIVE_INT64, &v) < 0)
        throw ReadError("HDF5: failed reading attribute '" + rName + "'");
    return v;
}

void write_attr_int(hid_t loc, const std::string& rName, std::int64_t v, hid_t ftype) {
    Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    Hid a(H5Acreate2(loc, rName.c_str(), ftype, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw WriteError("HDF5: could not create attribute '" + rName + "'");
    H5Awrite(a, H5T_NATIVE_INT64, &v);
}

std::string read_attr_string(hid_t loc, const std::string& rName) {
    Hid a(H5Aopen(loc, rName.c_str(), H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw ReadError("HDF5: missing attribute '" + rName + "'");
    Hid t(H5Aget_type(a), H5Tclose);
    if (H5Tis_variable_str(t) > 0) {
        char* p = nullptr;
        Hid mt(H5Tcopy(H5T_C_S1), H5Tclose);
        H5Tset_size(mt, H5T_VARIABLE);
        H5Tset_cset(mt, H5Tget_cset(t));
        if (H5Aread(a, mt, &p) < 0 || p == nullptr)
            throw ReadError("HDF5: failed reading attribute '" + rName + "'");
        std::string out(p);
        H5free_memory(p);
        return out;
    }
    std::size_t sz = H5Tget_size(t);
    std::vector<char> buf(sz + 1, '\0');
    Hid mt(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(mt, sz);
    H5Tset_cset(mt, H5Tget_cset(t));
    // NULLPAD memory type: converting a NULLPAD file string into a NULLTERM
    // memory string of the same size would truncate the last character to
    // make room for the terminator.
    H5Tset_strpad(mt, H5T_STR_NULLPAD);
    if (H5Aread(a, mt, buf.data()) < 0)
        throw ReadError("HDF5: failed reading attribute '" + rName + "'");
    // trim trailing NULs/spaces
    std::string out(buf.data(), strnlen(buf.data(), sz));
    while (!out.empty() && out.back() == ' ')
        out.pop_back();
    return out;
}

void write_attr_string(hid_t loc, const std::string& rName, const std::string& rValue) {
    Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    Hid t(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(t, H5T_VARIABLE);
    H5Tset_cset(t, H5T_CSET_UTF8);
    Hid a(H5Acreate2(loc, rName.c_str(), t, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw WriteError("HDF5: could not create attribute '" + rName + "'");
    const char* p = rValue.c_str();
    H5Awrite(a, t, &p);
}

std::vector<std::string> group_links(hid_t loc) {
    H5G_info_t info;
    H5Gget_info(loc, &info);
    std::vector<std::string> names;
    names.reserve(info.nlinks);
    for (hsize_t i = 0; i < info.nlinks; ++i) {
        ssize_t len =
            H5Lget_name_by_idx(loc, ".", H5_INDEX_NAME, H5_ITER_INC, i, nullptr, 0, H5P_DEFAULT);
        std::string name(static_cast<std::size_t>(len), '\0');
        H5Lget_name_by_idx(loc, ".", H5_INDEX_NAME, H5_ITER_INC, i, name.data(),
                           static_cast<std::size_t>(len) + 1, H5P_DEFAULT);
        names.push_back(std::move(name));
    }
    return names;
}

std::vector<std::string> group_links_crt(hid_t loc) {
    H5G_info_t info;
    H5Gget_info(loc, &info);
    std::vector<std::string> names;
    names.reserve(info.nlinks);
    for (hsize_t i = 0; i < info.nlinks; ++i) {
        ssize_t len = H5Lget_name_by_idx(loc, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, i, nullptr, 0,
                                         H5P_DEFAULT);
        if (len < 0)
            return group_links(loc);  // creation order not indexed
        std::string name(static_cast<std::size_t>(len), '\0');
        H5Lget_name_by_idx(loc, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, i, name.data(),
                           static_cast<std::size_t>(len) + 1, H5P_DEFAULT);
        names.push_back(std::move(name));
    }
    return names;
}

}  // namespace h5
}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
