#pragma once
//
// Shared HDF5 helpers for the HDF5-backed formats (CGNS, HMF, H5M, MED,
// XDMF-HDF). Only available when the extension is built with MESHIO_HAS_HDF5;
// the format sources are entirely #ifdef-guarded on that macro.
//
// Uses the version-stable classic C API with explicit-version names
// (H5Gcreate2/H5Dcreate2/H5Acreate2, ...) so both HDF5 1.10 and 2.x compile.

#ifdef MESHIO_HAS_HDF5

#include <hdf5.h>

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "meshio/exceptions.hpp"
#include "meshio/ndarray.hpp"

namespace meshio {
namespace h5 {

// RAII wrapper for an hid_t with its closer function.
class Hid {
public:
    using Closer = herr_t (*)(hid_t);
    Hid() = default;
    Hid(hid_t id, Closer closer) : id_(id), closer_(closer) {}
    Hid(Hid&& o) noexcept : id_(o.id_), closer_(o.closer_) { o.id_ = -1; }
    Hid& operator=(Hid&& o) noexcept {
        reset();
        id_ = o.id_;
        closer_ = o.closer_;
        o.id_ = -1;
        return *this;
    }
    Hid(const Hid&) = delete;
    Hid& operator=(const Hid&) = delete;
    ~Hid() { reset(); }

    void reset() {
        if (id_ >= 0 && closer_) closer_(id_);
        id_ = -1;
    }
    bool valid() const { return id_ >= 0; }
    hid_t get() const { return id_; }
    operator hid_t() const { return id_; }

private:
    hid_t id_ = -1;
    Closer closer_ = nullptr;
};

inline Hid open_file_read(const std::string& path) {
    Hid f(H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    if (!f.valid()) throw ReadError("HDF5: could not open file " + path);
    return f;
}

inline Hid create_file(const std::string& path) {
    Hid f(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose);
    if (!f.valid()) throw WriteError("HDF5: could not create file " + path);
    return f;
}

inline bool exists(hid_t loc, const std::string& name) {
    return H5Lexists(loc, name.c_str(), H5P_DEFAULT) > 0;
}

inline Hid open_group(hid_t loc, const std::string& name) {
    Hid g(H5Gopen2(loc, name.c_str(), H5P_DEFAULT), H5Gclose);
    if (!g.valid()) throw ReadError("HDF5: missing group '" + name + "'");
    return g;
}

inline Hid create_group(hid_t loc, const std::string& name) {
    Hid g(H5Gcreate2(loc, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
          H5Gclose);
    if (!g.valid()) throw WriteError("HDF5: could not create group '" + name + "'");
    return g;
}

// meshio DType -> native in-memory HDF5 type (for H5Dread/H5Dwrite).
inline hid_t native_type(DType dt) {
    switch (dt) {
        case DType::Float32: return H5T_NATIVE_FLOAT;
        case DType::Float64: return H5T_NATIVE_DOUBLE;
        case DType::Int8: return H5T_NATIVE_INT8;
        case DType::Int16: return H5T_NATIVE_INT16;
        case DType::Int32: return H5T_NATIVE_INT32;
        case DType::Int64: return H5T_NATIVE_INT64;
        case DType::UInt8: return H5T_NATIVE_UINT8;
        case DType::UInt16: return H5T_NATIVE_UINT16;
        case DType::UInt32: return H5T_NATIVE_UINT32;
        case DType::UInt64: return H5T_NATIVE_UINT64;
    }
    return H5T_NATIVE_DOUBLE;
}

// meshio DType -> little-endian file type (matches what h5py writes on x86).
inline hid_t file_type(DType dt) {
    switch (dt) {
        case DType::Float32: return H5T_IEEE_F32LE;
        case DType::Float64: return H5T_IEEE_F64LE;
        case DType::Int8: return H5T_STD_I8LE;
        case DType::Int16: return H5T_STD_I16LE;
        case DType::Int32: return H5T_STD_I32LE;
        case DType::Int64: return H5T_STD_I64LE;
        case DType::UInt8: return H5T_STD_U8LE;
        case DType::UInt16: return H5T_STD_U16LE;
        case DType::UInt32: return H5T_STD_U32LE;
        case DType::UInt64: return H5T_STD_U64LE;
    }
    return H5T_IEEE_F64LE;
}

// Stored datatype of a dataset/attribute -> meshio DType.
inline DType dtype_from_h5(hid_t type_id) {
    H5T_class_t cls = H5Tget_class(type_id);
    std::size_t sz = H5Tget_size(type_id);
    if (cls == H5T_FLOAT) return sz == 4 ? DType::Float32 : DType::Float64;
    if (cls == H5T_INTEGER) {
        bool is_signed = H5Tget_sign(type_id) != H5T_SGN_NONE;
        switch (sz) {
            case 1: return is_signed ? DType::Int8 : DType::UInt8;
            case 2: return is_signed ? DType::Int16 : DType::UInt16;
            case 4: return is_signed ? DType::Int32 : DType::UInt32;
            default: return is_signed ? DType::Int64 : DType::UInt64;
        }
    }
    throw ReadError("HDF5: unsupported datatype class");
}

// Read a full dataset into an owning NDArray (shape + dtype from the file).
// A dataset whose datatype is an ARRAY of a scalar type (h5py's "(n,) of
// k-tuples" trick, used by H5M) is returned with the array dims appended to
// the shape, i.e. as a plain (n, k) array.
inline NDArray read_dataset(hid_t loc, const std::string& name) {
    Hid d(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), H5Dclose);
    if (!d.valid()) throw ReadError("HDF5: missing dataset '" + name + "'");
    Hid space(H5Dget_space(d), H5Sclose);
    int ndim = H5Sget_simple_extent_ndims(space);
    std::vector<hsize_t> hdims(ndim > 0 ? ndim : 0);
    if (ndim > 0) H5Sget_simple_extent_dims(space, hdims.data(), nullptr);
    Hid dt(H5Dget_type(d), H5Tclose);

    std::vector<std::size_t> shape(hdims.begin(), hdims.end());
    if (shape.empty()) shape.push_back(1);  // scalar -> length-1

    DType mdt;
    if (H5Tget_class(dt) == H5T_ARRAY) {
        Hid base(H5Tget_super(dt), H5Tclose);
        mdt = dtype_from_h5(base);
        int arank = H5Tget_array_ndims(dt);
        std::vector<hsize_t> adims(arank > 0 ? arank : 0);
        if (arank > 0) H5Tget_array_dims2(dt, adims.data());
        for (hsize_t ad : adims) shape.push_back(static_cast<std::size_t>(ad));
    } else {
        mdt = dtype_from_h5(dt);
    }

    NDArray out(mdt, shape);
    if (out.size() > 0) {
        // For ARRAY-typed datasets the memory type must be the matching array
        // type; for scalar types the plain native type suffices.
        if (H5Tget_class(dt) == H5T_ARRAY) {
            int arank = H5Tget_array_ndims(dt);
            std::vector<hsize_t> adims(arank > 0 ? arank : 0);
            if (arank > 0) H5Tget_array_dims2(dt, adims.data());
            Hid mem(H5Tarray_create2(native_type(mdt), arank, adims.data()), H5Tclose);
            if (H5Dread(d, mem, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()) < 0)
                throw ReadError("HDF5: failed reading dataset '" + name + "'");
        } else if (H5Dread(d, native_type(mdt), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           out.data()) < 0) {
            throw ReadError("HDF5: failed reading dataset '" + name + "'");
        }
    }
    return out;
}

// Write a full dataset; gzip-compressed (chunked, chunk = full shape) when
// `gzip_level >= 0` and the data is non-empty.
inline void write_dataset(hid_t loc, const std::string& name, const NDArray& arr,
                          int gzip_level = -1) {
    std::vector<hsize_t> hdims(arr.shape().begin(), arr.shape().end());
    if (hdims.empty()) hdims.push_back(0);
    Hid space(H5Screate_simple(static_cast<int>(hdims.size()), hdims.data(), nullptr),
              H5Sclose);

    Hid dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose);
    if (gzip_level >= 0 && arr.size() > 0) {
        H5Pset_chunk(dcpl, static_cast<int>(hdims.size()), hdims.data());
        H5Pset_deflate(dcpl, static_cast<unsigned>(gzip_level));
    }

    Hid d(H5Dcreate2(loc, name.c_str(), file_type(arr.dtype()), space, H5P_DEFAULT,
                     dcpl, H5P_DEFAULT),
          H5Dclose);
    if (!d.valid()) throw WriteError("HDF5: could not create dataset '" + name + "'");
    if (arr.size() > 0) {
        if (H5Dwrite(d, native_type(arr.dtype()), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     arr.data()) < 0)
            throw WriteError("HDF5: failed writing dataset '" + name + "'");
    }
}

// ---- attribute helpers ----

inline bool has_attr(hid_t loc, const std::string& name) {
    return H5Aexists(loc, name.c_str()) > 0;
}

inline std::int64_t read_attr_int(hid_t loc, const std::string& name) {
    Hid a(H5Aopen(loc, name.c_str(), H5P_DEFAULT), H5Aclose);
    if (!a.valid()) throw ReadError("HDF5: missing attribute '" + name + "'");
    std::int64_t v = 0;
    if (H5Aread(a, H5T_NATIVE_INT64, &v) < 0)
        throw ReadError("HDF5: failed reading attribute '" + name + "'");
    return v;
}

inline void write_attr_int(hid_t loc, const std::string& name, std::int64_t v,
                           hid_t ftype = H5T_STD_I64LE) {
    Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    Hid a(H5Acreate2(loc, name.c_str(), ftype, space, H5P_DEFAULT, H5P_DEFAULT),
          H5Aclose);
    if (!a.valid()) throw WriteError("HDF5: could not create attribute '" + name + "'");
    H5Awrite(a, H5T_NATIVE_INT64, &v);
}

// Read a string attribute (fixed or variable length).
inline std::string read_attr_string(hid_t loc, const std::string& name) {
    Hid a(H5Aopen(loc, name.c_str(), H5P_DEFAULT), H5Aclose);
    if (!a.valid()) throw ReadError("HDF5: missing attribute '" + name + "'");
    Hid t(H5Aget_type(a), H5Tclose);
    if (H5Tis_variable_str(t) > 0) {
        char* p = nullptr;
        Hid mt(H5Tcopy(H5T_C_S1), H5Tclose);
        H5Tset_size(mt, H5T_VARIABLE);
        H5Tset_cset(mt, H5Tget_cset(t));
        if (H5Aread(a, mt, &p) < 0 || p == nullptr)
            throw ReadError("HDF5: failed reading attribute '" + name + "'");
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
        throw ReadError("HDF5: failed reading attribute '" + name + "'");
    // trim trailing NULs/spaces
    std::string out(buf.data(), strnlen(buf.data(), sz));
    while (!out.empty() && out.back() == ' ') out.pop_back();
    return out;
}

// Write a string attribute the way h5py does by default: variable-length UTF-8.
inline void write_attr_string(hid_t loc, const std::string& name,
                              const std::string& value) {
    Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    Hid t(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(t, H5T_VARIABLE);
    H5Tset_cset(t, H5T_CSET_UTF8);
    Hid a(H5Acreate2(loc, name.c_str(), t, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.valid()) throw WriteError("HDF5: could not create attribute '" + name + "'");
    const char* p = value.c_str();
    H5Awrite(a, t, &p);
}

// List the link names of a group, in creation/alphabetical (native) order.
inline std::vector<std::string> group_links(hid_t loc) {
    H5G_info_t info;
    H5Gget_info(loc, &info);
    std::vector<std::string> names;
    names.reserve(info.nlinks);
    for (hsize_t i = 0; i < info.nlinks; ++i) {
        ssize_t len = H5Lget_name_by_idx(loc, ".", H5_INDEX_NAME, H5_ITER_INC, i,
                                         nullptr, 0, H5P_DEFAULT);
        std::string name(static_cast<std::size_t>(len), '\0');
        H5Lget_name_by_idx(loc, ".", H5_INDEX_NAME, H5_ITER_INC, i, name.data(),
                           static_cast<std::size_t>(len) + 1, H5P_DEFAULT);
        names.push_back(std::move(name));
    }
    return names;
}

// Silence HDF5's default stderr error stack (we convert to exceptions).
struct SilenceErrors {
    SilenceErrors() {
        H5Eget_auto2(H5E_DEFAULT, &old_func_, &old_data_);
        H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    }
    ~SilenceErrors() { H5Eset_auto2(H5E_DEFAULT, old_func_, old_data_); }
    H5E_auto2_t old_func_ = nullptr;
    void* old_data_ = nullptr;
};

}  // namespace h5
}  // namespace meshio

#endif  // MESHIO_HAS_HDF5
