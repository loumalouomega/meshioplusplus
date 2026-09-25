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
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

// External includes
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
#include <zlib.h>
#endif

// Project includes
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"

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

namespace {

// Direct chunk I/O (H5Dwrite_chunk / H5Dread_chunk) arrived in HDF5 1.10.2,
// and gzip in parallel needs zlib itself.
#if defined(MESHIOPLUSPLUS_HAS_ZLIB) && H5_VERSION_GE(1, 10, 2)
#define MESHIOPLUSPLUS_H5_PARALLEL_DEFLATE 1
#endif

/// Target bytes per chunk of a compressed dataset (the appendable datasets'
/// policy below): large enough that the chunk index stays small, small
/// enough that the chunks spread over the threads.
constexpr std::size_t kH5uChunkBytes = std::size_t{1} << 20;

/**
 * @brief Rows per chunk for a gzip dataset of `rDims` holding `ElemBytes`-byte
 * elements: every row whose `kH5uChunkBytes` allow, or all of them when the
 * whole dataset fits (one chunk, as every dataset had before chunking).
 */
hsize_t h5u_chunk_rows(const std::vector<hsize_t>& rDims, std::size_t ElemBytes) {
    std::size_t row_bytes = ElemBytes;
    for (std::size_t i = 1; i < rDims.size(); ++i)
        row_bytes *= static_cast<std::size_t>(rDims[i]);
    const hsize_t rows = std::max<hsize_t>(1, kH5uChunkBytes / std::max<std::size_t>(1, row_bytes));
    return std::min(rows, rDims[0]);
}

#ifdef MESHIOPLUSPLUS_H5_PARALLEL_DEFLATE
/**
 * @brief Write `rArr` into the chunked, deflate-filtered dataset `d` chunk by
 * chunk: compress every chunk in parallel exactly as HDF5's deflate filter
 * would (`compress2` at `Level`, the last chunk zero-padded to full size, as
 * the library pads it with the default fill value), then write them in
 * ascending order. Returns false, having written nothing, when the file type
 * is not the memory layout byte for byte.
 */
bool h5u_write_deflate_chunks(hid_t d, const NDArray& rArr, const std::vector<hsize_t>& rDims,
                              hsize_t ChunkRows, int Level) {
    if (H5Tequal(file_type(rArr.Dtype()), native_type(rArr.Dtype())) <= 0)
        return false;
    const std::size_t row_bytes = rArr.Nbytes() / static_cast<std::size_t>(rDims[0]);
    const std::size_t chunk_bytes = static_cast<std::size_t>(ChunkRows) * row_bytes;
    const std::size_t nchunks = static_cast<std::size_t>((rDims[0] + ChunkRows - 1) / ChunkRows);
    const auto* src = reinterpret_cast<const Bytef*>(rArr.Data());
    std::vector<std::vector<Bytef>> packed(nchunks);
    std::vector<int> status(nchunks, Z_OK);
    parallel_for(
        nchunks,
        [&](std::size_t c) {
            const std::size_t first = c * chunk_bytes;
            const std::size_t have = std::min(chunk_bytes, rArr.Nbytes() - first);
            std::vector<Bytef> padded;
            const Bytef* in = src + first;
            if (have < chunk_bytes) {
                padded.assign(chunk_bytes, 0);
                std::memcpy(padded.data(), in, have);
                in = padded.data();
            }
            uLongf n = compressBound(static_cast<uLong>(chunk_bytes));
            packed[c].resize(n);
            status[c] = compress2(packed[c].data(), &n, in, static_cast<uLong>(chunk_bytes), Level);
            packed[c].resize(n);
        },
        1);
    std::vector<hsize_t> offset(rDims.size(), 0);
    for (std::size_t c = 0; c < nchunks; ++c) {
        if (status[c] != Z_OK)
            throw WriteError("HDF5: zlib could not compress a chunk");
        offset[0] = static_cast<hsize_t>(c) * ChunkRows;
        if (H5Dwrite_chunk(d, H5P_DEFAULT, 0, offset.data(), packed[c].size(), packed[c].data()) <
            0)
            throw WriteError("HDF5: failed writing a chunk");
        std::vector<Bytef>().swap(packed[c]);
    }
    return true;
}

/**
 * @brief Read the gzip-compressed, row-chunked dataset `d` straight into
 * `rOut` with the chunks inflated in parallel (the I/O stays serial: HDF5 is
 * not thread-safe). Returns false, having read nothing, for any layout it does
 * not cover -- not chunked, a filter pipeline other than deflate alone, chunks
 * that split a row, one chunk only, a chunk never written, or a file type
 * that is not the memory layout -- so the caller falls back to `H5Dread`.
 */
bool h5u_read_deflate_chunks(hid_t d, hid_t FileType, NDArray& rOut,
                             const std::vector<hsize_t>& rDims) {
    if (rDims.empty() || rOut.Size() == 0 || H5Tequal(FileType, native_type(rOut.Dtype())) <= 0)
        return false;
    Hid dcpl(H5Dget_create_plist(d), H5Pclose);
    if (!dcpl.Valid() || H5Pget_layout(dcpl) != H5D_CHUNKED || H5Pget_nfilters(dcpl) != 1)
        return false;
    unsigned flags = 0;
    std::size_t nelmts = 0;
    unsigned filter_config = 0;
    if (H5Pget_filter2(dcpl, 0, &flags, &nelmts, nullptr, 0, nullptr, &filter_config) !=
        H5Z_FILTER_DEFLATE)
        return false;
    std::vector<hsize_t> chunk(rDims.size(), 0);
    if (H5Pget_chunk(dcpl, static_cast<int>(chunk.size()), chunk.data()) !=
        static_cast<int>(chunk.size()))
        return false;
    for (std::size_t i = 1; i < rDims.size(); ++i)
        if (chunk[i] != rDims[i])
            return false;  // chunks that split a row: not this path's shape
    if (chunk[0] == 0 || chunk[0] >= rDims[0])
        return false;  // one chunk: nothing to spread over threads

    const std::size_t row_bytes = rOut.Nbytes() / static_cast<std::size_t>(rDims[0]);
    const std::size_t chunk_bytes = static_cast<std::size_t>(chunk[0]) * row_bytes;
    const std::size_t nchunks = static_cast<std::size_t>((rDims[0] + chunk[0] - 1) / chunk[0]);
    std::vector<std::vector<Bytef>> packed(nchunks);
    std::vector<std::uint32_t> masks(nchunks, 0);
    std::vector<hsize_t> offset(rDims.size(), 0);
    for (std::size_t c = 0; c < nchunks; ++c) {
        offset[0] = static_cast<hsize_t>(c) * chunk[0];
        hsize_t stored = 0;
        if (H5Dget_chunk_storage_size(d, offset.data(), &stored) < 0 || stored == 0)
            return false;  // unallocated: H5Dread supplies the fill value
        packed[c].resize(static_cast<std::size_t>(stored));
#if H5_VERSION_GE(2, 0, 0)
        std::size_t got = packed[c].size();
        const herr_t rc =
            H5Dread_chunk2(d, H5P_DEFAULT, offset.data(), &masks[c], packed[c].data(), &got);
        if (rc >= 0 && got != packed[c].size())
            return false;
#else
        const herr_t rc = H5Dread_chunk(d, H5P_DEFAULT, offset.data(), &masks[c], packed[c].data());
#endif
        if (rc < 0)
            return false;
    }

    auto* dst = reinterpret_cast<Bytef*>(rOut.Data());
    std::vector<std::uint8_t> ok(nchunks, 0);
    parallel_for(
        nchunks,
        [&](std::size_t c) {
            const std::size_t first = c * chunk_bytes;
            const std::size_t want = std::min(chunk_bytes, rOut.Nbytes() - first);
            if (masks[c] & 1u) {  // the filter was skipped: the chunk is stored raw
                if (packed[c].size() >= want) {
                    std::memcpy(dst + first, packed[c].data(), want);
                    ok[c] = 1;
                }
                return;
            }
            // The last chunk inflates to a whole chunk (HDF5 pads it); only
            // its leading `want` bytes belong to the dataset.
            std::vector<Bytef> tail;
            Bytef* out = dst + first;
            if (want < chunk_bytes) {
                tail.resize(chunk_bytes);
                out = tail.data();
            }
            uLongf n = static_cast<uLongf>(chunk_bytes);
            const int rc =
                uncompress(out, &n, packed[c].data(), static_cast<uLong>(packed[c].size()));
            if (rc != Z_OK || n != chunk_bytes)
                return;
            if (!tail.empty())
                std::memcpy(dst + first, tail.data(), want);
            ok[c] = 1;
        },
        1);
    for (std::uint8_t k : ok)
        if (!k)
            throw ReadError("HDF5: a compressed chunk does not inflate to its declared size");
    return true;
}
#endif

}  // namespace

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

    // Uninitialized: H5Dread overwrites every element (or throws), so zero-
    // filling first would be a wasted pass over the whole dataset.
    NDArray out = NDArray::Uninit(mdt, shape);
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
        } else {
#ifdef MESHIOPLUSPLUS_H5_PARALLEL_DEFLATE
            if (h5u_read_deflate_chunks(d, dt, out, hdims))
                return out;
#endif
            if (H5Dread(d, native_type(mdt), H5S_ALL, H5S_ALL, H5P_DEFAULT, out.Data()) < 0)
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
    const bool compress = gzip_level >= 0 && rArr.Size() > 0;
    hsize_t chunk_rows = 0;
    if (compress) {
        std::vector<hsize_t> chunk = hdims;
        chunk_rows = chunk[0] = h5u_chunk_rows(hdims, dtype_size(rArr.Dtype()));
        H5Pset_chunk(dcpl, static_cast<int>(chunk.size()), chunk.data());
        H5Pset_deflate(dcpl, static_cast<unsigned>(gzip_level));
    }

    Hid d(H5Dcreate2(loc, rName.c_str(), file_type(rArr.Dtype()), space, H5P_DEFAULT, dcpl,
                     H5P_DEFAULT),
          H5Dclose);
    if (!d.Valid())
        throw WriteError("HDF5: could not create dataset '" + rName + "'");
#ifdef MESHIOPLUSPLUS_H5_PARALLEL_DEFLATE
    if (compress && chunk_rows < hdims[0] &&
        h5u_write_deflate_chunks(d, rArr, hdims, chunk_rows, gzip_level))
        return;
#endif
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

void write_attr_string_fixed(hid_t loc, const std::string& rName, const std::string& rValue) {
    if (rValue.empty())
        throw WriteError("HDF5: fixed-length string attribute '" + rName +
                         "' cannot be empty (HDF5 has no zero-size string type)");
    for (unsigned char c : rValue) {
        if (c > 0x7F)
            throw WriteError("HDF5: attribute '" + rName + "' must be ASCII");
    }
    Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    Hid t(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(t, rValue.size());
    H5Tset_cset(t, H5T_CSET_ASCII);
    H5Tset_strpad(t, H5T_STR_NULLPAD);
    Hid a(H5Acreate2(loc, rName.c_str(), t, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw WriteError("HDF5: could not create attribute '" + rName + "'");
    if (H5Awrite(a, t, rValue.data()) < 0)
        throw WriteError("HDF5: failed writing attribute '" + rName + "'");
}

void write_attr_int_array(hid_t loc, const std::string& rName,
                          const std::vector<std::int64_t>& rValues, hid_t ftype) {
    if (rValues.empty())
        throw WriteError("HDF5: integer array attribute '" + rName + "' cannot be empty");
    const hsize_t dim = rValues.size();
    Hid space(H5Screate_simple(1, &dim, nullptr), H5Sclose);
    Hid a(H5Acreate2(loc, rName.c_str(), ftype, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw WriteError("HDF5: could not create attribute '" + rName + "'");
    if (H5Awrite(a, H5T_NATIVE_INT64, rValues.data()) < 0)
        throw WriteError("HDF5: failed writing attribute '" + rName + "'");
}

std::vector<std::int64_t> read_attr_int_array(hid_t loc, const std::string& rName) {
    Hid a(H5Aopen(loc, rName.c_str(), H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw ReadError("HDF5: missing attribute '" + rName + "'");
    Hid t(H5Aget_type(a), H5Tclose);
    if (H5Tget_class(t) != H5T_INTEGER)
        throw ReadError("HDF5: attribute '" + rName + "' is not an integer");
    Hid space(H5Aget_space(a), H5Sclose);
    const hssize_t n = H5Sget_simple_extent_npoints(space);
    if (n < 0)
        throw ReadError("HDF5: could not size attribute '" + rName + "'");
    std::vector<std::int64_t> out(static_cast<std::size_t>(n));
    if (n > 0 && H5Aread(a, H5T_NATIVE_INT64, out.data()) < 0)
        throw ReadError("HDF5: failed reading attribute '" + rName + "'");
    return out;
}

std::vector<std::size_t> dataset_shape(hid_t loc, const std::string& rName) {
    Hid d(H5Dopen2(loc, rName.c_str(), H5P_DEFAULT), H5Dclose);
    if (!d.Valid())
        throw ReadError("HDF5: missing dataset '" + rName + "'");
    Hid space(H5Dget_space(d), H5Sclose);
    const int ndim = H5Sget_simple_extent_ndims(space);
    std::vector<hsize_t> dims(ndim > 0 ? ndim : 0);
    if (ndim > 0)
        H5Sget_simple_extent_dims(space, dims.data(), nullptr);
    std::vector<std::size_t> shape(dims.begin(), dims.end());
    if (shape.empty())
        shape.push_back(1);  // scalar -> length-1, as read_dataset reports it
    return shape;
}

NDArray read_dataset_rows(hid_t loc, const std::string& rName, std::size_t Row0,
                          std::size_t Count) {
    Hid d(H5Dopen2(loc, rName.c_str(), H5P_DEFAULT), H5Dclose);
    if (!d.Valid())
        throw ReadError("HDF5: missing dataset '" + rName + "'");
    Hid space(H5Dget_space(d), H5Sclose);
    const int ndim = H5Sget_simple_extent_ndims(space);
    if (ndim < 1)
        throw ReadError("HDF5: dataset '" + rName + "' is scalar; it has no rows to slice");
    std::vector<hsize_t> dims(static_cast<std::size_t>(ndim));
    H5Sget_simple_extent_dims(space, dims.data(), nullptr);
    if (Row0 > dims[0] || Count > dims[0] - Row0)
        throw ReadError("HDF5: rows [" + std::to_string(Row0) + ", " +
                        std::to_string(Row0 + Count) + ") are outside dataset '" + rName + "' of " +
                        std::to_string(dims[0]) + " rows");

    Hid dt(H5Dget_type(d), H5Tclose);
    const DType mdt = dtype_from_h5(dt);
    std::vector<std::size_t> shape(dims.begin(), dims.end());
    shape[0] = Count;
    NDArray out = NDArray::Uninit(mdt, shape);  // H5Dread fills it, or throws
    if (Count == 0 || out.Size() == 0)
        return out;

    std::vector<hsize_t> start(static_cast<std::size_t>(ndim), 0);
    std::vector<hsize_t> count(dims);
    start[0] = Row0;
    count[0] = Count;
    if (H5Sselect_hyperslab(space, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr) <
        0)
        throw ReadError("HDF5: could not select rows of dataset '" + rName + "'");
    Hid mem(H5Screate_simple(ndim, count.data(), nullptr), H5Sclose);
    if (H5Dread(d, native_type(mdt), mem, space, H5P_DEFAULT, out.Data()) < 0)
        throw ReadError("HDF5: failed reading rows of dataset '" + rName + "'");
    return out;
}

void create_appendable_dataset(hid_t loc, const std::string& rName, DType dt,
                               const std::vector<std::size_t>& rRowShape, std::size_t ChunkRows,
                               int gzip_level) {
    std::size_t row_elems = 1;
    for (std::size_t e : rRowShape) {
        if (e == 0)
            throw WriteError("HDF5: appendable dataset '" + rName +
                             "' cannot have a zero-length row dimension");
        row_elems *= e;
    }
    if (ChunkRows == 0) {
        // Aim for ~1 MiB chunks: small enough that a short series does not pad a
        // whole megabyte per dataset, large enough that a long one is not all metadata.
        constexpr std::size_t kTargetBytes = std::size_t{1} << 20;
        const std::size_t row_bytes = std::max<std::size_t>(1, row_elems * dtype_size(dt));
        ChunkRows = std::max<std::size_t>(1, kTargetBytes / row_bytes);
    }

    const int rank = static_cast<int>(rRowShape.size()) + 1;
    std::vector<hsize_t> dims(static_cast<std::size_t>(rank), 0);
    std::vector<hsize_t> maxdims(static_cast<std::size_t>(rank), 0);
    std::vector<hsize_t> chunk(static_cast<std::size_t>(rank), 0);
    maxdims[0] = H5S_UNLIMITED;
    chunk[0] = ChunkRows;
    for (std::size_t i = 0; i < rRowShape.size(); ++i) {
        dims[i + 1] = maxdims[i + 1] = chunk[i + 1] = rRowShape[i];
    }
    Hid space(H5Screate_simple(rank, dims.data(), maxdims.data()), H5Sclose);
    Hid dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose);
    if (H5Pset_chunk(dcpl, rank, chunk.data()) < 0)
        throw WriteError("HDF5: could not chunk dataset '" + rName + "'");
    if (gzip_level >= 0)
        H5Pset_deflate(dcpl, static_cast<unsigned>(gzip_level));
    Hid d(H5Dcreate2(loc, rName.c_str(), file_type(dt), space, H5P_DEFAULT, dcpl, H5P_DEFAULT),
          H5Dclose);
    if (!d.Valid())
        throw WriteError("HDF5: could not create dataset '" + rName + "'");
}

void append_rows(hid_t loc, const std::string& rName, const NDArray& rRows) {
    Hid d(H5Dopen2(loc, rName.c_str(), H5P_DEFAULT), H5Dclose);
    if (!d.Valid())
        throw WriteError("HDF5: missing dataset '" + rName + "'");
    Hid space(H5Dget_space(d), H5Sclose);
    const int ndim = H5Sget_simple_extent_ndims(space);
    if (ndim < 1 || static_cast<std::size_t>(ndim) != rRows.Shape().size())
        throw WriteError("HDF5: appended rows for '" + rName + "' have rank " +
                         std::to_string(rRows.Shape().size()) + ", the dataset has rank " +
                         std::to_string(ndim));
    std::vector<hsize_t> dims(static_cast<std::size_t>(ndim));
    std::vector<hsize_t> maxdims(static_cast<std::size_t>(ndim));
    H5Sget_simple_extent_dims(space, dims.data(), maxdims.data());
    if (maxdims[0] != H5S_UNLIMITED)
        throw WriteError("HDF5: dataset '" + rName + "' is not extendable");
    for (std::size_t i = 1; i < dims.size(); ++i) {
        if (dims[i] != rRows.Shape()[i])
            throw WriteError("HDF5: appended rows for '" + rName +
                             "' have row shape differing "
                             "from the dataset's in dimension " +
                             std::to_string(i));
    }
    const std::size_t n = rRows.Shape()[0];
    if (n == 0)
        return;

    std::vector<hsize_t> newdims(dims);
    newdims[0] = dims[0] + n;
    if (H5Dset_extent(d, newdims.data()) < 0)
        throw WriteError("HDF5: could not extend dataset '" + rName + "'");
    Hid fspace(H5Dget_space(d), H5Sclose);  // the extent changed: re-fetch
    std::vector<hsize_t> start(static_cast<std::size_t>(ndim), 0);
    std::vector<hsize_t> count(dims);
    start[0] = dims[0];
    count[0] = n;
    if (H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr) <
        0)
        throw WriteError("HDF5: could not select rows of dataset '" + rName + "'");
    Hid mem(H5Screate_simple(ndim, count.data(), nullptr), H5Sclose);
    if (H5Dwrite(d, native_type(rRows.Dtype()), mem, fspace, H5P_DEFAULT, rRows.Data()) < 0)
        throw WriteError("HDF5: failed appending to dataset '" + rName + "'");
}

std::size_t dataset_num_rows(hid_t loc, const std::string& rName) {
    Hid d(H5Dopen2(loc, rName.c_str(), H5P_DEFAULT), H5Dclose);
    if (!d.Valid())
        throw ReadError("HDF5: missing dataset '" + rName + "'");
    Hid space(H5Dget_space(d), H5Sclose);
    const int ndim = H5Sget_simple_extent_ndims(space);
    if (ndim < 1)
        throw ReadError("HDF5: dataset '" + rName + "' is scalar; it has no rows");
    std::vector<hsize_t> dims(static_cast<std::size_t>(ndim));
    H5Sget_simple_extent_dims(space, dims.data(), nullptr);
    return static_cast<std::size_t>(dims[0]);
}

void rewrite_dataset(hid_t loc, const std::string& rName, const NDArray& rArr) {
    Hid d(H5Dopen2(loc, rName.c_str(), H5P_DEFAULT), H5Dclose);
    if (!d.Valid())
        throw WriteError("HDF5: missing dataset '" + rName + "'");
    Hid space(H5Dget_space(d), H5Sclose);
    const int ndim = H5Sget_simple_extent_ndims(space);
    std::vector<hsize_t> dims(static_cast<std::size_t>(ndim > 0 ? ndim : 0));
    if (ndim > 0)
        H5Sget_simple_extent_dims(space, dims.data(), nullptr);
    if (dims.empty())
        dims.push_back(1);
    const std::vector<std::size_t> want(dims.begin(), dims.end());
    if (rArr.Shape() != want)
        throw WriteError("HDF5: cannot rewrite dataset '" + rName + "' with a different shape");
    if (rArr.Size() == 0)
        return;
    if (H5Dwrite(d, native_type(rArr.Dtype()), H5S_ALL, H5S_ALL, H5P_DEFAULT, rArr.Data()) < 0)
        throw WriteError("HDF5: failed rewriting dataset '" + rName + "'");
}

Hid create_group_crt(hid_t loc, const std::string& rName) {
    Hid gcpl(H5Pcreate(H5P_GROUP_CREATE), H5Pclose);
    if (H5Pset_link_creation_order(gcpl, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED) < 0)
        throw WriteError("HDF5: could not request creation-order tracking for group '" + rName +
                         "'");
    Hid g(H5Gcreate2(loc, rName.c_str(), H5P_DEFAULT, gcpl, H5P_DEFAULT), H5Gclose);
    if (!g.Valid())
        throw WriteError("HDF5: could not create group '" + rName + "'");
    return g;
}

void create_soft_link(hid_t loc, const std::string& rName, const std::string& rTargetPath) {
    if (H5Lcreate_soft(rTargetPath.c_str(), loc, rName.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0)
        throw WriteError("HDF5: could not create soft link '" + rName + "' -> '" + rTargetPath +
                         "'");
}

bool is_soft_link(hid_t loc, const std::string& rName) {
    if (H5Lexists(loc, rName.c_str(), H5P_DEFAULT) <= 0)
        return false;
    H5L_info_t info;
    if (H5Lget_info(loc, rName.c_str(), &info, H5P_DEFAULT) < 0)
        return false;
    return info.type == H5L_TYPE_SOFT;
}

std::string soft_link_target(hid_t loc, const std::string& rName) {
    if (!is_soft_link(loc, rName))
        throw ReadError("HDF5: '" + rName + "' is not a soft link");
    H5L_info_t info;
    if (H5Lget_info(loc, rName.c_str(), &info, H5P_DEFAULT) < 0)
        throw ReadError("HDF5: could not inspect link '" + rName + "'");
    std::vector<char> buf(info.u.val_size + 1, '\0');
    if (H5Lget_val(loc, rName.c_str(), buf.data(), buf.size(), H5P_DEFAULT) < 0)
        throw ReadError("HDF5: could not read soft link '" + rName + "'");
    return std::string(buf.data());
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

std::vector<CompoundMember> compound_members(hid_t loc, const std::string& rName) {
    Hid d(H5Dopen2(loc, rName.c_str(), H5P_DEFAULT), H5Dclose);
    if (!d.Valid())
        throw ReadError("HDF5: missing dataset '" + rName + "'");
    Hid dt(H5Dget_type(d), H5Tclose);
    if (H5Tget_class(dt) != H5T_COMPOUND)
        throw ReadError("HDF5: dataset '" + rName + "' is not a compound table");
    const int n = H5Tget_nmembers(dt);
    std::vector<CompoundMember> out;
    out.reserve(n > 0 ? static_cast<std::size_t>(n) : 0u);
    for (int i = 0; i < n; ++i) {
        CompoundMember m;
        char* name = H5Tget_member_name(dt, static_cast<unsigned>(i));
        if (name != nullptr) {
            m.mName = name;
            H5free_memory(name);
        }
        Hid mt(H5Tget_member_type(dt, static_cast<unsigned>(i)), H5Tclose);
        Hid base;
        hid_t elem = mt;
        if (H5Tget_class(mt) == H5T_ARRAY) {
            const int arank = H5Tget_array_ndims(mt);
            std::vector<hsize_t> adims(arank > 0 ? static_cast<std::size_t>(arank) : 0u);
            if (arank > 0)
                H5Tget_array_dims2(mt, adims.data());
            m.mDims.assign(adims.begin(), adims.end());
            base = Hid(H5Tget_super(mt), H5Tclose);
            elem = base;
        }
        const H5T_class_t cls = H5Tget_class(elem);
        if (cls == H5T_INTEGER || cls == H5T_FLOAT) {
            m.mNumeric = true;
            m.mDtype = dtype_from_h5(elem);
        }
        out.push_back(std::move(m));
    }
    return out;
}

NDArray read_compound_member(hid_t loc, const std::string& rName, const std::string& rMember,
                             std::size_t Row0, std::size_t Count, const DType* pAs) {
    const std::vector<CompoundMember> members = compound_members(loc, rName);
    const CompoundMember* found = nullptr;
    for (const CompoundMember& m : members)
        if (m.mName == rMember)
            found = &m;
    if (found == nullptr || !found->mNumeric)
        throw ReadError("HDF5: compound table '" + rName + "' has no numeric member '" + rMember +
                        "'");
    Hid d(H5Dopen2(loc, rName.c_str(), H5P_DEFAULT), H5Dclose);
    Hid space(H5Dget_space(d), H5Sclose);
    if (H5Sget_simple_extent_ndims(space) != 1)
        throw ReadError("HDF5: compound table '" + rName + "' is not one-dimensional");
    hsize_t rows = 0;
    H5Sget_simple_extent_dims(space, &rows, nullptr);
    if (Row0 > rows || Count > rows - Row0)
        throw ReadError("HDF5: rows [" + std::to_string(Row0) + ", " +
                        std::to_string(Row0 + Count) + ") are outside compound table '" + rName +
                        "' of " + std::to_string(rows) + " rows");

    std::vector<std::size_t> shape{Count};
    shape.insert(shape.end(), found->mDims.begin(), found->mDims.end());
    const DType dt = pAs != nullptr ? *pAs : found->mDtype;
    NDArray out(dt, shape);
    if (Count == 0 || out.Size() == 0)
        return out;

    Hid elem;
    if (found->mDims.empty()) {
        elem = Hid(H5Tcopy(native_type(dt)), H5Tclose);
    } else {
        std::vector<hsize_t> adims(found->mDims.begin(), found->mDims.end());
        elem = Hid(
            H5Tarray_create2(native_type(dt), static_cast<unsigned>(adims.size()), adims.data()),
            H5Tclose);
    }
    Hid mem_type(H5Tcreate(H5T_COMPOUND, H5Tget_size(elem)), H5Tclose);
    if (!mem_type.Valid() || H5Tinsert(mem_type, rMember.c_str(), 0, elem) < 0)
        throw ReadError("HDF5: could not build a memory type for '" + rName + "/" + rMember + "'");
    const hsize_t start = Row0;
    const hsize_t count = Count;
    if (H5Sselect_hyperslab(space, H5S_SELECT_SET, &start, nullptr, &count, nullptr) < 0)
        throw ReadError("HDF5: could not select rows of compound table '" + rName + "'");
    Hid mem(H5Screate_simple(1, &count, nullptr), H5Sclose);
    if (H5Dread(d, mem_type, mem, space, H5P_DEFAULT, out.Data()) < 0)
        throw ReadError("HDF5: failed reading member '" + rMember + "' of '" + rName + "'");
    return out;
}

}  // namespace h5
}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
