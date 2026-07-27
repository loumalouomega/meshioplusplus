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
 * @file hdf5_util.hpp
 * @brief Shared low-level HDF5 helpers used by every C++ format that stores
 * data in an HDF5 container (MED, XDMF's `Format="HDF"` DataItems).
 *
 * Provides: an RAII handle wrapper (`Hid`) so every `H5*` resource is closed
 * exactly once even under exceptions; dataset read/write helpers that
 * translate between `meshioplusplus::DType` and HDF5's native/file types
 * (matching what h5py writes on x86: little-endian file types via
 * `file_type()`); scalar/string attribute helpers matching h5py's own
 * variable-length UTF-8 convention; group link listing in both name order
 * and (where the file tracks it) creation order, the latter needed where
 * block order carries meaning (e.g. MED's `MAI` cell blocks, whose order
 * must align with `cell_data`/`cell_sets`); and `SilenceErrors`, which
 * suppresses HDF5's default stderr error-stack printing so failures surface
 * only as the C++ exceptions this codebase converts them to (`ReadError`/
 * `WriteError`). This entire header compiles to nothing when
 * `MESHIOPLUSPLUS_HAS_HDF5` is not defined, i.e. when the build has no HDF5
 * library — the HDF-dependent C++ code paths are then simply absent and
 * callers fall back to the pure-Python (h5py-based) implementation.
 *
 * `Hid`'s methods and `SilenceErrors` are small RAII glue (a handful of lines
 * each), so they stay inline; every free function below is called once per
 * dataset/attribute/group (never once per element), so bodies live in
 * `src/cpp/src/detail/hdf5_util.cpp` (itself `#ifdef`-guarded to compile to an
 * empty TU when HDF5 is absent, matching every HDF5-backed format's own
 * `.cpp`) rather than inline here.
 */

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// External includes
#include <hdf5.h>

// System includes
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace h5 {

/**
 * @brief RAII wrapper for an HDF5 `hid_t` handle, paired with the `H5*Close`
 * function that must release it.
 *
 * Move-only (copying an `hid_t` would double-close it): moving transfers
 * ownership and leaves the source handle invalid (`mId = -1`). Implicitly
 * convertible to `hid_t` so it can be passed straight into `H5*` C API
 * calls. `Valid()` reports whether the handle is currently open (id `>= 0`).
 */
class Hid {
public:
    using Closer = herr_t (*)(hid_t);
    Hid() = default;
    Hid(hid_t id, Closer closer) : mId(id), mCloser(closer) {}
    Hid(Hid&& o) noexcept : mId(o.mId), mCloser(o.mCloser) { o.mId = -1; }
    Hid& operator=(Hid&& o) noexcept {
        Reset();
        mId = o.mId;
        mCloser = o.mCloser;
        o.mId = -1;
        return *this;
    }
    Hid(const Hid&) = delete;
    Hid& operator=(const Hid&) = delete;
    ~Hid() { Reset(); }

    void Reset() {
        if (mId >= 0 && mCloser)
            mCloser(mId);
        mId = -1;
    }
    bool Valid() const { return mId >= 0; }
    hid_t Get() const { return mId; }
    operator hid_t() const { return mId; }

private:
    hid_t mId = -1;
    Closer mCloser = nullptr;
};

/**
 * @brief Opens an existing HDF5 file read-only.
 * @param rPath Filesystem path of the file to open.
 * @return Owning `Hid` for the open file.
 * @throws ReadError if the file cannot be opened.
 */
MESHIOPLUSPLUS_API Hid open_file_read(const std::string& rPath);

/**
 * @brief Creates a new HDF5 file, truncating any existing file at `path`.
 * @param rPath Filesystem path of the file to create.
 * @return Owning `Hid` for the new file.
 * @throws WriteError if the file cannot be created.
 */
MESHIOPLUSPLUS_API Hid create_file(const std::string& rPath);

/**
 * @brief Opens an existing HDF5 file for reading *and* writing.
 *
 * Distinct from both `open_file_read` (read-only) and `create_file` (which
 * truncates): appending to a file that already carries datasets needs a third
 * mode, and reaching for `create_file` would silently destroy the very data an
 * append is meant to continue.
 *
 * @param rPath Filesystem path of the file to open.
 * @return Owning `Hid` for the open file.
 * @throws WriteError if the file cannot be opened read-write.
 */
MESHIOPLUSPLUS_API Hid open_file_rw(const std::string& rPath);

/**
 * @brief Flush a file's buffers to disk without closing it.
 *
 * What makes a partially-written series readable by another process: HDF5
 * buffers aggressively, so an un-flushed `.h5` can be unopenable even though
 * every dataset was "written".
 *
 * @param rFile An open file handle; an invalid handle is a no-op.
 * @throws WriteError if the flush fails.
 */
MESHIOPLUSPLUS_API void flush_file(Hid& rFile);

/**
 * @brief The names of every link directly under a group or file.
 * @param loc Group or file handle to list.
 * @return Link names in the container's own iteration order.
 */
MESHIOPLUSPLUS_API std::vector<std::string> link_names(hid_t loc);

/**
 * @brief Whether a link named `name` exists directly under group/file `loc`.
 * @param loc Group or file handle to look under.
 * @param rName Link name to test.
 * @return `true` if the link exists.
 */
MESHIOPLUSPLUS_API bool exists(hid_t loc, const std::string& rName);

/**
 * @brief Opens an existing HDF5 group.
 * @param loc Parent group or file handle.
 * @param rName Name of the group to open.
 * @return Owning `Hid` for the opened group.
 * @throws ReadError if the group does not exist.
 */
MESHIOPLUSPLUS_API Hid open_group(hid_t loc, const std::string& rName);

/**
 * @brief Creates a new HDF5 group.
 * @param loc Parent group or file handle.
 * @param rName Name of the group to create.
 * @return Owning `Hid` for the new group.
 * @throws WriteError if the group cannot be created.
 */
MESHIOPLUSPLUS_API Hid create_group(hid_t loc, const std::string& rName);

/**
 * @brief Maps a `meshioplusplus::DType` to the native in-memory HDF5 type used
 * for `H5Dread`/`H5Dwrite` (host byte order/representation, not the
 * on-disk file type — see `file_type()` for that).
 * @param dt The dtype to convert.
 * @return The matching `H5T_NATIVE_*` constant (defaults to
 *         `H5T_NATIVE_DOUBLE` for an unrecognized/invalid `dt`).
 */
MESHIOPLUSPLUS_API hid_t native_type(DType dt);

/**
 * @brief Maps a `meshioplusplus::DType` to the on-disk (file) HDF5 type to use
 * when creating a dataset/attribute.
 *
 * Always little-endian (`H5T_*LE`), matching what h5py writes on x86, so
 * files produced by the C++ writer are byte-for-byte compatible with the
 * pure-Python/h5py writer's output.
 * @param dt The dtype to convert.
 * @return The matching `H5T_*LE` constant (defaults to `H5T_IEEE_F64LE`).
 */
MESHIOPLUSPLUS_API hid_t file_type(DType dt);

/**
 * @brief Converts a stored HDF5 datatype (of a dataset or attribute) to the
 * corresponding `meshioplusplus::DType`.
 * @param type_id HDF5 type id, as returned e.g. by `H5Dget_type`.
 * @return The matching `DType`.
 * @throws ReadError if `type_id`'s class is neither float nor integer.
 */
MESHIOPLUSPLUS_API DType dtype_from_h5(hid_t type_id);

/**
 * @brief Reads a full HDF5 dataset into a freshly-allocated, owning `NDArray`.
 *
 * The output shape and dtype are taken from the file. A dataset whose
 * datatype is an `ARRAY` of a scalar type — h5py's "(n,) of k-tuples" trick,
 * used e.g. by MED's `H5M` node coordinates — is unpacked into a plain
 * `(n, k)` `NDArray` by appending the array dimensions to the dataset's
 * shape, rather than exposed as a compound/array-typed element.
 * A scalar (0-dimensional) dataset comes back with shape `{1}`.
 *
 * @param loc Group or file handle the dataset lives under.
 * @param rName Name of the dataset to read.
 * @return A new owning `NDArray` holding the dataset's contents.
 * @throws ReadError if the dataset is missing or the read fails.
 */
MESHIOPLUSPLUS_API NDArray read_dataset(hid_t loc, const std::string& rName);

/**
 * @brief Writes a full dataset in one call, optionally gzip-compressed.
 *
 * When `gzip_level >= 0` and `arr` is non-empty, the dataset is created
 * chunked with a single chunk spanning the whole shape and gzip deflate
 * filtering enabled at that level; otherwise it is a plain contiguous
 * dataset. Uses `file_type(arr.Dtype())` for the on-disk type and
 * `native_type(arr.Dtype())` for the in-memory transfer type.
 *
 * @param loc Group or file handle to create the dataset under.
 * @param rName Name for the new dataset.
 * @param rArr Data to write; its shape and dtype determine the dataset's.
 * @param gzip_level gzip compression level (0-9), or negative to disable
 *                   compression (the default).
 * @throws WriteError if the dataset cannot be created or the write fails.
 */
MESHIOPLUSPLUS_API void write_dataset(hid_t loc, const std::string& rName, const NDArray& rArr,
                                      int gzip_level = -1);

// ---- attribute helpers ----

/**
 * @brief Whether an attribute named `name` exists on `loc`.
 * @param loc Object (group/dataset/file) to check.
 * @param rName Attribute name to test.
 * @return `true` if the attribute exists.
 */
MESHIOPLUSPLUS_API bool has_attr(hid_t loc, const std::string& rName);

/**
 * @brief Reads a scalar integer attribute.
 * @param loc Object the attribute is attached to.
 * @param rName Attribute name.
 * @return The attribute's value as `int64_t`.
 * @throws ReadError if the attribute is missing or unreadable.
 */
MESHIOPLUSPLUS_API std::int64_t read_attr_int(hid_t loc, const std::string& rName);

/**
 * @brief Writes a scalar integer attribute.
 * @param loc Object to attach the attribute to.
 * @param rName Attribute name.
 * @param v Value to write.
 * @param ftype On-disk integer type to store as (default `H5T_STD_I64LE`).
 * @throws WriteError if the attribute cannot be created.
 */
MESHIOPLUSPLUS_API void write_attr_int(hid_t loc, const std::string& rName, std::int64_t v,
                                       hid_t ftype = H5T_STD_I64LE);

/**
 * @brief Reads a string attribute, handling both variable- and fixed-length
 * HDF5 string encodings.
 *
 * For a fixed-length (`NULLPAD`) string, reads into a same-sized buffer
 * (converting to `NULLTERM` would truncate the last character to make room
 * for a terminator) and then trims trailing NUL bytes and spaces.
 * @param loc Object the attribute is attached to.
 * @param rName Attribute name.
 * @return The attribute's value as a `std::string`.
 * @throws ReadError if the attribute is missing or unreadable.
 */
MESHIOPLUSPLUS_API std::string read_attr_string(hid_t loc, const std::string& rName);

/**
 * @brief Writes a string attribute the way h5py does by default:
 * variable-length, UTF-8-tagged.
 *
 * Matching h5py's convention keeps files produced by the C++ writer
 * byte-for-byte compatible with the Python/h5py writer's output.
 * @param loc Object to attach the attribute to.
 * @param rName Attribute name.
 * @param rValue String value to write.
 * @throws WriteError if the attribute cannot be created.
 */
MESHIOPLUSPLUS_API void write_attr_string(hid_t loc, const std::string& rName,
                                          const std::string& rValue);

/**
 * @brief Lists the link (child) names directly under a group, in HDF5's
 * default name-index iteration order.
 * @param loc Group handle to list.
 * @return Child link names, in name order.
 */
MESHIOPLUSPLUS_API std::vector<std::string> group_links(hid_t loc);

/**
 * @brief Like `group_links`, but iterates in HDF5 link *creation* order when
 * the group tracks it (matching h5py's iteration order on
 * `track_order=True` files); silently falls back to name order otherwise.
 *
 * Needed wherever the order children were created in is semantically
 * significant rather than incidental — e.g. MED's `MAI` cell-block groups,
 * whose order must line up with the corresponding entries in `cell_data`/
 * `cell_sets`, which are positional (not keyed by group name).
 * @param loc Group handle to list.
 * @return Child link names, in creation order if indexed, else name order.
 */
MESHIOPLUSPLUS_API std::vector<std::string> group_links_crt(hid_t loc);

/**
 * @brief RAII guard that silences HDF5's default stderr error-stack printing
 * for its lifetime, restoring the previous handler on destruction.
 *
 * The library's own error reporting is redundant here since every failure
 * this codebase cares about is converted to a `ReadError`/`WriteError`
 * exception; without this guard, HDF5 would additionally dump a raw error
 * stack to stderr on every recoverable failure (e.g. a probing "does this
 * attribute exist" call that's expected to fail sometimes).
 */
struct SilenceErrors {
    SilenceErrors() {
        H5Eget_auto2(H5E_DEFAULT, &mOldFunc, &mOldData);
        H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
    }
    ~SilenceErrors() { H5Eset_auto2(H5E_DEFAULT, mOldFunc, mOldData); }
    H5E_auto2_t mOldFunc = nullptr;
    void* mOldData = nullptr;
};

}  // namespace h5
}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
