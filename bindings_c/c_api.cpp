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

/**
 * @file c_api.cpp
 * @brief Implementation of the meshio++ C API (bindings_c/include/
 *        meshioplusplus/meshioplusplus.h), compiled into the installable
 *        `libmeshioplusplus` shared library.
 *
 * Like the WASM binding (and unlike `bindings/_core.cpp`), this is a flat,
 * whole-mesh interface built exclusively on the uniform mesh API
 * (mesh_api.hpp) -- it therefore compiles unchanged under every mesh backend
 * (MESHIO/NATIVE/KRATOS) -- and dispatches formats through the shared
 * registry (registry.hpp). The three ABI rules the implementation enforces:
 *
 *  - No C++ exception ever crosses `extern "C"`: every public function body
 *    runs inside guarded()/guarded_ptr(), which map ReadError/WriteError/
 *    anything else to a `mio_status` plus a thread-local message retrievable
 *    via mio_last_error().
 *  - Setters copy caller memory into owning NDArrays (validated first);
 *    getters hand out pointers into mesh-owned storage, which every backend
 *    keeps stable until the next mutating call (KRATOS serves accessors from
 *    its persistent staging mesh).
 *  - Strings cross via the caller-buffer/required-length protocol -- never a
 *    `c_str()` of a possibly-temporary (CellView::Type()'s return category
 *    differs per backend).
 */

// System includes
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <exception>
#include <new>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/meshioplusplus.h"

#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/skin.hpp"

struct mio_mesh {
    meshioplusplus::Mesh mMesh;
};

struct mio_reorder_result {
    mio_mesh mMesh;  // owns the renumbered mesh; borrowed via mio_reorder_result_mesh
    meshioplusplus::NDArray mNodePerm;
    std::vector<meshioplusplus::NDArray> mCellPerms;
};

namespace {

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

// The C header's MIO_CELL_TYPES list must mirror the C++ enum exactly; any
// added/removed/reordered entry on either side fails to compile right here.
#define MIO_CELL_TYPE_CHECK(Name)                                                              \
    static_assert(                                                                             \
        static_cast<int>(meshioplusplus::CellType::Name) == static_cast<int>(MIO_CELL_##Name), \
        "meshioplusplus.h cell-type list drifted from cell_type.hpp");
MIO_CELL_TYPES(MIO_CELL_TYPE_CHECK)
#undef MIO_CELL_TYPE_CHECK
static_assert(static_cast<int>(meshioplusplus::CellType::Custom) ==
                  static_cast<int>(MIO_CELL_Custom),
              "meshioplusplus.h cell-type list drifted from cell_type.hpp");

#ifndef MIO_VERSION_STRING
#define MIO_VERSION_STRING "unknown"
#endif

thread_local std::string g_last_error;

mio_status fail(mio_status code, std::string msg) {
    g_last_error = std::move(msg);
    return code;
}

// Every fallible extern "C" body runs inside one of these two: no exception
// crosses the ABI, and the thread-local message is always set on failure.
template <class F>
mio_status guarded(F&& f) {
    try {
        return f();
    } catch (const meshioplusplus::ReadError& e) {
        return fail(MIO_ERR_READ, e.what());
    } catch (const meshioplusplus::WriteError& e) {
        return fail(MIO_ERR_WRITE, e.what());
    } catch (const std::bad_alloc&) {
        return fail(MIO_ERR_INTERNAL, "meshio++: out of memory");
    } catch (const std::exception& e) {
        return fail(MIO_ERR_INTERNAL, e.what());
    } catch (...) {
        return fail(MIO_ERR_INTERNAL, "meshio++: unknown error");
    }
}

// Variant for functions returning a pointer/count instead of a status.
template <class T, class F>
T guarded_ptr(T on_error, F&& f) {
    T result = on_error;
    guarded([&]() -> mio_status {
        result = f();
        return MIO_OK;
    });
    return result;
}

bool to_dtype(mio_dtype in, DType& rOut) {
    switch (in) {
        case MIO_FLOAT32:
            rOut = DType::Float32;
            return true;
        case MIO_FLOAT64:
            rOut = DType::Float64;
            return true;
        case MIO_INT8:
            rOut = DType::Int8;
            return true;
        case MIO_INT16:
            rOut = DType::Int16;
            return true;
        case MIO_INT32:
            rOut = DType::Int32;
            return true;
        case MIO_INT64:
            rOut = DType::Int64;
            return true;
        case MIO_UINT8:
            rOut = DType::UInt8;
            return true;
        case MIO_UINT16:
            rOut = DType::UInt16;
            return true;
        case MIO_UINT32:
            rOut = DType::UInt32;
            return true;
        case MIO_UINT64:
            rOut = DType::UInt64;
            return true;
    }
    return false;
}

mio_dtype from_dtype(DType in) {
    switch (in) {
        case DType::Float32:
            return MIO_FLOAT32;
        case DType::Float64:
            return MIO_FLOAT64;
        case DType::Int8:
            return MIO_INT8;
        case DType::Int16:
            return MIO_INT16;
        case DType::Int32:
            return MIO_INT32;
        case DType::Int64:
            return MIO_INT64;
        case DType::UInt8:
            return MIO_UINT8;
        case DType::UInt16:
            return MIO_UINT16;
        case DType::UInt32:
            return MIO_UINT32;
        case DType::UInt64:
            return MIO_UINT64;
    }
    return MIO_FLOAT64;  // unreachable
}

// Validate an (ndim, shape) pair from the caller and return it as the vector
// NDArray wants; MIO_MAX_NDIM bounds the rank in both directions of the API.
bool to_shape(int32_t ndim, const int64_t* pShape, std::vector<std::size_t>& rOut) {
    if (ndim < 1 || ndim > MIO_MAX_NDIM || !pShape)
        return false;
    rOut.assign(static_cast<std::size_t>(ndim), 0);
    for (int32_t i = 0; i < ndim; ++i) {
        if (pShape[i] < 0)
            return false;
        rOut[static_cast<std::size_t>(i)] = static_cast<std::size_t>(pShape[i]);
    }
    return true;
}

// Owning NDArray copied from a caller buffer (setter rule: setters copy).
NDArray copy_in(DType dt, std::vector<std::size_t> shape, const void* pData) {
    NDArray out = NDArray::Uninit(dt, std::move(shape));
    if (out.Nbytes() > 0)
        std::memcpy(out.Data(), pData, out.Nbytes());
    return out;
}

// String-getter protocol (header rule 5): copy what fits, NUL-terminate,
// return the untruncated length.
int64_t copy_string(const std::string& rStr, char* pBuf, int64_t buflen) {
    if (pBuf && buflen > 0) {
        std::size_t n = std::min(rStr.size(), static_cast<std::size_t>(buflen - 1));
        std::memcpy(pBuf, rStr.data(), n);
        pBuf[n] = '\0';
    }
    return static_cast<int64_t>(rStr.size());
}

std::string format_or_empty(const char* pFormat) {
    return pFormat ? std::string(pFormat) : std::string();
}

std::string unknown_format_message(const std::string& rFormat, bool for_write) {
    std::string msg = for_write
                          ? "meshio++: unknown, read-only, or unsupported format '" + rFormat + "'"
                          : "meshio++: unknown or unsupported format '" + rFormat + "'";
    if (const char* dep = meshioplusplus::registry_compiled_out(rFormat))
        msg += " (this build has no " + std::string(dep) + " support)";
    return msg;
}

// Fill the (data, dtype, ndim, shape) out-params from a mesh-owned NDArray;
// any out-param may be NULL.
mio_status array_out(const NDArray& rArr, const void** ppData, mio_dtype* pDtype, int32_t* pNdim,
                     int64_t* pShape) {
    const auto& shape = rArr.Shape();
    if (shape.size() > MIO_MAX_NDIM)
        return fail(MIO_ERR_INTERNAL, "meshio++: array rank exceeds MIO_MAX_NDIM");
    if (ppData)
        *ppData = rArr.Data();
    if (pDtype)
        *pDtype = from_dtype(rArr.Dtype());
    if (pNdim)
        *pNdim = static_cast<int32_t>(shape.size());
    if (pShape) {
        for (std::size_t i = 0; i < shape.size(); ++i)
            pShape[i] = static_cast<int64_t>(shape[i]);
    }
    return MIO_OK;
}

// Shared body of the three named-data setters (they differ only in the
// leading-dimension check and the uniform-API call).
template <class AddFn>
mio_status add_named_array(mio_mesh* pMesh, const char* pName, mio_dtype dtype, int32_t ndim,
                           const int64_t* pShape, const void* pData, const char* pWhat,
                           std::int64_t required_dim0, AddFn&& add) {
    DType dt;
    std::vector<std::size_t> shape;
    if (!pMesh || !pName)
        return fail(MIO_ERR_INVALID_ARG,
                    std::string("meshio++: bad ") + pWhat + " argument (NULL mesh/name)");
    if (!to_dtype(dtype, dt))
        return fail(MIO_ERR_INVALID_ARG, std::string("meshio++: bad ") + pWhat + " dtype");
    if (!to_shape(ndim, pShape, shape))
        return fail(MIO_ERR_INVALID_ARG, std::string("meshio++: bad ") + pWhat +
                                             " shape (rank 1.." + std::to_string(MIO_MAX_NDIM) +
                                             ", non-negative extents)");
    if (required_dim0 >= 0 && shape[0] != static_cast<std::size_t>(required_dim0))
        return fail(MIO_ERR_INVALID_ARG, std::string("meshio++: ") + pWhat + " '" + pName +
                                             "' shape[0] is " + std::to_string(shape[0]) +
                                             ", expected " + std::to_string(required_dim0));
    NDArray arr = NDArray::Uninit(dt, std::move(shape));
    if (arr.Nbytes() > 0) {
        if (!pData)
            return fail(MIO_ERR_INVALID_ARG,
                        std::string("meshio++: bad ") + pWhat + " argument (NULL data)");
        std::memcpy(arr.Data(), pData, arr.Nbytes());
    }
    add(std::move(arr));
    return MIO_OK;
}

bool block_in_range(const mio_mesh* pMesh, int64_t block) {
    return pMesh && block >= 0 && static_cast<std::size_t>(block) < pMesh->mMesh.NumCellBlocks();
}

}  // namespace

extern "C" {

/* ------------------------------------------------------------------ */
/* Version / build introspection                                       */
/* ------------------------------------------------------------------ */

const char* mio_version(void) {
    return MIO_VERSION_STRING;
}

const char* mio_mesh_backend(void) {
    return meshioplusplus::mesh_backend_name();
}

int mio_format_readable(const char* format) {
    return guarded_ptr(0, [&]() -> int {
        return format && meshioplusplus::registry_readers().count(format) ? 1 : 0;
    });
}

int mio_format_writable(const char* format) {
    return guarded_ptr(0, [&]() -> int {
        return format && meshioplusplus::registry_writers().count(format) ? 1 : 0;
    });
}

const char* mio_last_error(void) {
    return g_last_error.c_str();
}

/* ------------------------------------------------------------------ */
/* Cell-type metadata                                                  */
/* ------------------------------------------------------------------ */

const char* mio_cell_type_name(mio_cell_type t) {
    if (t < 0 || t > MIO_CELL_Custom)
        return "";
    return meshioplusplus::cell_type_name(static_cast<meshioplusplus::CellType>(t)).c_str();
}

mio_cell_type mio_cell_type_from_name(const char* name) {
    if (!name)
        return MIO_CELL_Custom;
    return static_cast<mio_cell_type>(meshioplusplus::cell_type_from_name(name));
}

int mio_cell_type_num_nodes(mio_cell_type t) {
    if (t < 0 || t > MIO_CELL_Custom)
        return -1;
    return meshioplusplus::cell_type_num_nodes(static_cast<meshioplusplus::CellType>(t));
}

int mio_cell_type_dimension(mio_cell_type t) {
    if (t < 0 || t > MIO_CELL_Custom)
        return -1;
    return meshioplusplus::cell_type_dimension(static_cast<meshioplusplus::CellType>(t));
}

/* ------------------------------------------------------------------ */
/* Lifecycle & file I/O                                                */
/* ------------------------------------------------------------------ */

mio_mesh* mio_mesh_create(void) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), []() { return new mio_mesh(); });
}

void mio_mesh_free(mio_mesh* mesh) {
    delete mesh;
}

mio_mesh* mio_read(const char* path, const char* format) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!path)
            throw meshioplusplus::ReadError("meshio++: path is NULL");
        std::string fmt;
        try {
            fmt = meshioplusplus::resolve_format(path, format_or_empty(format));
        } catch (const meshioplusplus::ReadError&) {
            // Extension gave nothing: fall back to a conservative content sniff.
            fmt = meshioplusplus::sniff_format(path);
            if (fmt.empty())
                throw;
        }
        auto it = meshioplusplus::registry_readers().find(fmt);
        if (it == meshioplusplus::registry_readers().end())
            throw meshioplusplus::ReadError(unknown_format_message(fmt, /*for_write=*/false));
        return new mio_mesh{it->second(path)};
    });
}

mio_status mio_write(const char* path, const mio_mesh* mesh, const char* format) {
    return guarded([&]() -> mio_status {
        if (!path || !mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: path/mesh is NULL");
        std::string fmt = meshioplusplus::resolve_format(path, format_or_empty(format));
        auto it = meshioplusplus::registry_writers().find(fmt);
        if (it == meshioplusplus::registry_writers().end())
            return fail(MIO_ERR_NOT_FOUND, unknown_format_message(fmt, /*for_write=*/true));
        it->second(path, mesh->mMesh);
        return MIO_OK;
    });
}

mio_status mio_convert(const char* in_path, const char* in_format, const char* out_path,
                       const char* out_format) {
    return guarded([&]() -> mio_status {
        if (!in_path || !out_path)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: path is NULL");
        std::string rfmt = meshioplusplus::resolve_format(in_path, format_or_empty(in_format));
        std::string wfmt = meshioplusplus::resolve_format(out_path, format_or_empty(out_format));
        auto rit = meshioplusplus::registry_readers().find(rfmt);
        if (rit == meshioplusplus::registry_readers().end())
            return fail(MIO_ERR_NOT_FOUND, unknown_format_message(rfmt, /*for_write=*/false));
        auto wit = meshioplusplus::registry_writers().find(wfmt);
        if (wit == meshioplusplus::registry_writers().end())
            return fail(MIO_ERR_NOT_FOUND, unknown_format_message(wfmt, /*for_write=*/true));
        wit->second(out_path, rit->second(in_path));
        return MIO_OK;
    });
}

/* ------------------------------------------------------------------ */
/* Mesh operations                                                     */
/* ------------------------------------------------------------------ */

mio_mesh* mio_extract_surface(const mio_mesh* mesh, int record_parent_ids) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return new mio_mesh{meshioplusplus::extract_surface(mesh->mMesh, record_parent_ids != 0)};
    });
}

mio_mesh* mio_extract_skin(const mio_mesh* mesh, int linearize) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return new mio_mesh{meshioplusplus::extract_skin(mesh->mMesh, linearize != 0)};
    });
}

mio_mesh* mio_attach_quality(const mio_mesh* mesh) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return new mio_mesh{meshioplusplus::attach_quality(mesh->mMesh)};
    });
}

mio_status mio_quality_counts(const mio_mesh* mesh, int64_t* num_cells, int64_t* num_inverted,
                              int64_t* num_degenerate) {
    return guarded([&]() -> mio_status {
        if (!mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh is NULL");
        meshioplusplus::QualityReport rep = meshioplusplus::compute_quality(mesh->mMesh);
        if (num_cells)
            *num_cells = rep.mNumCells;
        if (num_inverted)
            *num_inverted = rep.mNumInverted;
        if (num_degenerate)
            *num_degenerate = rep.mNumDegenerate;
        return MIO_OK;
    });
}

int64_t mio_sniff_format(const char* path, char* buf, int64_t buflen) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!path)
            throw meshioplusplus::ReadError("meshio++: path is NULL");
        return copy_string(meshioplusplus::sniff_format(path), buf, buflen);
    });
}

mio_reorder_result* mio_reorder(const mio_mesh* mesh, const char* method) {
    return guarded_ptr(static_cast<mio_reorder_result*>(nullptr), [&]() -> mio_reorder_result* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        if (!method)
            throw meshioplusplus::ReadError("meshio++: method is NULL");
        meshioplusplus::ReorderResult res =
            meshioplusplus::reorder(mesh->mMesh, meshioplusplus::reorder_method_from_name(method));
        auto* out =
            new mio_reorder_result{mio_mesh{std::move(res.mMesh)}, std::move(res.mNodePermutation),
                                   std::move(res.mCellPermutations)};
        return out;
    });
}

const mio_mesh* mio_reorder_result_mesh(const mio_reorder_result* result) {
    return result ? &result->mMesh : nullptr;
}

mio_status mio_reorder_result_node_perm(const mio_reorder_result* result, const void** data,
                                        mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        const NDArray& a = result->mNodePerm;
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

int64_t mio_reorder_result_num_cell_perms(const mio_reorder_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return static_cast<int64_t>(result->mCellPerms.size());
    });
}

mio_status mio_reorder_result_cell_perm(const mio_reorder_result* result, int64_t block,
                                        const void** data, mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (block < 0 || static_cast<std::size_t>(block) >= result->mCellPerms.size())
            return fail(MIO_ERR_NOT_FOUND, "meshio++: cell-permutation block index out of range");
        const NDArray& a = result->mCellPerms[static_cast<std::size_t>(block)];
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

mio_mesh* mio_reorder_result_take_mesh(mio_reorder_result* result) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return new mio_mesh{std::move(result->mMesh.mMesh)};
    });
}

void mio_reorder_result_free(mio_reorder_result* result) {
    delete result;
}

int64_t mio_compute_bandwidth(const mio_mesh* mesh) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return meshioplusplus::compute_bandwidth(mesh->mMesh);
    });
}

/* ------------------------------------------------------------------ */
/* Building a mesh (setters copy)                                      */
/* ------------------------------------------------------------------ */

mio_status mio_mesh_set_points(mio_mesh* mesh, mio_dtype dtype, int64_t num_points, int64_t dim,
                               const void* xyz) {
    return guarded([&]() -> mio_status {
        DType dt;
        if (!mesh || num_points < 0 || dim <= 0 || (!xyz && num_points > 0))
            return fail(MIO_ERR_INVALID_ARG, "meshio++: bad mio_mesh_set_points argument");
        if (!to_dtype(dtype, dt) || (dt != DType::Float32 && dt != DType::Float64))
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: points dtype must be MIO_FLOAT32 or MIO_FLOAT64");
        mesh->mMesh.AssignPoints(copy_in(
            dt, {static_cast<std::size_t>(num_points), static_cast<std::size_t>(dim)}, xyz));
        return MIO_OK;
    });
}

mio_status mio_mesh_add_cell_block(mio_mesh* mesh, const char* cell_type, int64_t num_cells,
                                   int64_t nodes_per_cell, mio_dtype dtype,
                                   const void* connectivity) {
    return guarded([&]() -> mio_status {
        if (!mesh || !cell_type || num_cells < 0 || nodes_per_cell <= 0 ||
            (!connectivity && num_cells > 0))
            return fail(MIO_ERR_INVALID_ARG, "meshio++: bad mio_mesh_add_cell_block argument");
        if (dtype != MIO_INT32 && dtype != MIO_INT64)
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: connectivity dtype must be MIO_INT32 or MIO_INT64");
        const auto ct = meshioplusplus::cell_type_from_name(cell_type);
        const int fixed = meshioplusplus::cell_type_num_nodes(ct);
        if (ct != meshioplusplus::CellType::Custom && fixed > 0 && fixed != nodes_per_cell)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: cell type '" + std::string(cell_type) +
                                                 "' has " + std::to_string(fixed) +
                                                 " nodes per cell, got " +
                                                 std::to_string(nodes_per_cell));
        // Widen to Int64 during the copy -- the core's connectivity type.
        NDArray conn = NDArray::Uninit(DType::Int64, {static_cast<std::size_t>(num_cells),
                                                      static_cast<std::size_t>(nodes_per_cell)});
        const std::size_t n = conn.Size();
        if (dtype == MIO_INT64) {
            if (n > 0)
                std::memcpy(conn.Data(), connectivity, conn.Nbytes());
        } else {
            const auto* src = static_cast<const std::int32_t*>(connectivity);
            std::int64_t* dst = conn.As<std::int64_t>();
            for (std::size_t i = 0; i < n; ++i)
                dst[i] = src[i];
        }
        mesh->mMesh.AddCellBlock(cell_type, std::move(conn));
        return MIO_OK;
    });
}

mio_status mio_mesh_add_point_data(mio_mesh* mesh, const char* name, mio_dtype dtype, int32_t ndim,
                                   const int64_t* shape, const void* data) {
    return guarded([&]() -> mio_status {
        const std::int64_t npoints = mesh ? static_cast<std::int64_t>(mesh->mMesh.NumPoints()) : -1;
        return add_named_array(mesh, name, dtype, ndim, shape, data, "mio_mesh_add_point_data",
                               npoints,
                               [&](NDArray a) { mesh->mMesh.AddPointData(name, std::move(a)); });
    });
}

mio_status mio_mesh_append_cell_data(mio_mesh* mesh, const char* name, mio_dtype dtype,
                                     int32_t ndim, const int64_t* shape, const void* data) {
    return guarded([&]() -> mio_status {
        // The array being appended belongs to the next block in order: block
        // index = how many arrays this field already has.
        std::int64_t required = -1;
        if (mesh && name) {
            const std::size_t next =
                mesh->mMesh.HasCellData(name) ? mesh->mMesh.CellDataNumBlocks(name) : 0;
            if (next >= mesh->mMesh.NumCellBlocks())
                return fail(MIO_ERR_INVALID_ARG,
                            std::string("meshio++: cell_data '") + name + "' already has one " +
                                "array per cell block (add cell blocks before their data)");
            required = static_cast<std::int64_t>(mesh->mMesh.Cells(next).NumCells());
        }
        return add_named_array(mesh, name, dtype, ndim, shape, data, "mio_mesh_append_cell_data",
                               required,
                               [&](NDArray a) { mesh->mMesh.AppendCellData(name, std::move(a)); });
    });
}

mio_status mio_mesh_add_field_data(mio_mesh* mesh, const char* name, mio_dtype dtype, int32_t ndim,
                                   const int64_t* shape, const void* data) {
    return guarded([&]() -> mio_status {
        return add_named_array(mesh, name, dtype, ndim, shape, data, "mio_mesh_add_field_data",
                               /*required_dim0=*/-1,
                               [&](NDArray a) { mesh->mMesh.AddFieldData(name, std::move(a)); });
    });
}

/* ------------------------------------------------------------------ */
/* Reading a mesh back (getters are zero-copy)                         */
/* ------------------------------------------------------------------ */

int64_t mio_mesh_num_points(const mio_mesh* mesh) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return static_cast<int64_t>(mesh->mMesh.NumPoints());
    });
}

int64_t mio_mesh_point_dim(const mio_mesh* mesh) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return static_cast<int64_t>(mesh->mMesh.PointDim());
    });
}

mio_status mio_mesh_get_points(const mio_mesh* mesh, const void** data, mio_dtype* dtype) {
    return guarded([&]() -> mio_status {
        if (!mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh is NULL");
        const NDArray& points = mesh->mMesh.Points();
        return array_out(points, data, dtype, nullptr, nullptr);
    });
}

int64_t mio_mesh_num_cell_blocks(const mio_mesh* mesh) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return static_cast<int64_t>(mesh->mMesh.NumCellBlocks());
    });
}

mio_status mio_mesh_cell_block_info(const mio_mesh* mesh, int64_t block, int64_t* num_cells,
                                    int64_t* nodes_per_cell, int32_t* is_ragged) {
    return guarded([&]() -> mio_status {
        if (!mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh is NULL");
        if (!block_in_range(mesh, block))
            return fail(MIO_ERR_NOT_FOUND,
                        "meshio++: cell block " + std::to_string(block) + " out of range");
        const auto view = mesh->mMesh.Cells(static_cast<std::size_t>(block));
        if (num_cells)
            *num_cells = static_cast<int64_t>(view.NumCells());
        if (nodes_per_cell)
            *nodes_per_cell = static_cast<int64_t>(view.NodesPerCell());
        if (is_ragged)
            *is_ragged = view.IsRagged() ? 1 : 0;
        return MIO_OK;
    });
}

int64_t mio_mesh_cell_block_type(const mio_mesh* mesh, int64_t block, char* buf, int64_t buflen) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!block_in_range(mesh, block))
            throw meshioplusplus::ReadError("meshio++: cell block " + std::to_string(block) +
                                            " out of range");
        // Bind through a std::string so per-backend Type() return categories
        // (reference vs temporary) both stay valid for the copy.
        const std::string type = mesh->mMesh.Cells(static_cast<std::size_t>(block)).Type();
        return copy_string(type, buf, buflen);
    });
}

mio_status mio_mesh_cell_block_conn(const mio_mesh* mesh, int64_t block, const void** conn,
                                    mio_dtype* dtype) {
    return guarded([&]() -> mio_status {
        if (!mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh is NULL");
        if (!block_in_range(mesh, block))
            return fail(MIO_ERR_NOT_FOUND,
                        "meshio++: cell block " + std::to_string(block) + " out of range");
        const auto view = mesh->mMesh.Cells(static_cast<std::size_t>(block));
        if (view.IsRagged())
            return fail(MIO_ERR_UNSUPPORTED,
                        "meshio++: ragged cell blocks are not accessible through the C API yet");
        const NDArray& c = view.Conn();
        return array_out(c, conn, dtype, nullptr, nullptr);
    });
}

/* Named-data accessors: the three families (point/cell/field) share the same
 * shape; small macros would obscure more than they save, so they are spelled
 * out. Names are returned in ascending lexicographic order -- the uniform
 * API's *DataNames() guarantee, identical on every backend. */

int64_t mio_mesh_num_point_data(const mio_mesh* mesh) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return static_cast<int64_t>(mesh->mMesh.NumPointData());
    });
}

int64_t mio_mesh_point_data_name(const mio_mesh* mesh, int64_t index, char* buf, int64_t buflen) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        const auto names = mesh->mMesh.PointDataNames();
        if (index < 0 || static_cast<std::size_t>(index) >= names.size())
            throw meshioplusplus::ReadError("meshio++: point_data index " + std::to_string(index) +
                                            " out of range");
        return copy_string(names[static_cast<std::size_t>(index)], buf, buflen);
    });
}

mio_status mio_mesh_get_point_data(const mio_mesh* mesh, const char* name, const void** data,
                                   mio_dtype* dtype, int32_t* ndim, int64_t* shape) {
    return guarded([&]() -> mio_status {
        if (!mesh || !name)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh/name is NULL");
        if (!mesh->mMesh.HasPointData(name))
            return fail(MIO_ERR_NOT_FOUND,
                        "meshio++: no point_data named '" + std::string(name) + "'");
        return array_out(mesh->mMesh.PointData(name), data, dtype, ndim, shape);
    });
}

int64_t mio_mesh_num_cell_data(const mio_mesh* mesh) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return static_cast<int64_t>(mesh->mMesh.NumCellData());
    });
}

int64_t mio_mesh_cell_data_name(const mio_mesh* mesh, int64_t index, char* buf, int64_t buflen) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        const auto names = mesh->mMesh.CellDataNames();
        if (index < 0 || static_cast<std::size_t>(index) >= names.size())
            throw meshioplusplus::ReadError("meshio++: cell_data index " + std::to_string(index) +
                                            " out of range");
        return copy_string(names[static_cast<std::size_t>(index)], buf, buflen);
    });
}

int64_t mio_mesh_cell_data_num_blocks(const mio_mesh* mesh, const char* name) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!mesh || !name)
            throw meshioplusplus::ReadError("meshio++: mesh/name is NULL");
        if (!mesh->mMesh.HasCellData(name))
            throw meshioplusplus::ReadError("meshio++: no cell_data named '" + std::string(name) +
                                            "'");
        return static_cast<int64_t>(mesh->mMesh.CellDataNumBlocks(name));
    });
}

mio_status mio_mesh_get_cell_data(const mio_mesh* mesh, const char* name, int64_t block,
                                  const void** data, mio_dtype* dtype, int32_t* ndim,
                                  int64_t* shape) {
    return guarded([&]() -> mio_status {
        if (!mesh || !name)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh/name is NULL");
        if (!mesh->mMesh.HasCellData(name))
            return fail(MIO_ERR_NOT_FOUND,
                        "meshio++: no cell_data named '" + std::string(name) + "'");
        if (block < 0 || static_cast<std::size_t>(block) >= mesh->mMesh.CellDataNumBlocks(name))
            return fail(MIO_ERR_NOT_FOUND, "meshio++: cell_data '" + std::string(name) +
                                               "' block " + std::to_string(block) +
                                               " out of range");
        return array_out(mesh->mMesh.CellData(name, static_cast<std::size_t>(block)), data, dtype,
                         ndim, shape);
    });
}

int64_t mio_mesh_num_field_data(const mio_mesh* mesh) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return static_cast<int64_t>(mesh->mMesh.NumFieldData());
    });
}

int64_t mio_mesh_field_data_name(const mio_mesh* mesh, int64_t index, char* buf, int64_t buflen) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        const auto names = mesh->mMesh.FieldDataNames();
        if (index < 0 || static_cast<std::size_t>(index) >= names.size())
            throw meshioplusplus::ReadError("meshio++: field_data index " + std::to_string(index) +
                                            " out of range");
        return copy_string(names[static_cast<std::size_t>(index)], buf, buflen);
    });
}

mio_status mio_mesh_get_field_data(const mio_mesh* mesh, const char* name, const void** data,
                                   mio_dtype* dtype, int32_t* ndim, int64_t* shape) {
    return guarded([&]() -> mio_status {
        if (!mesh || !name)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh/name is NULL");
        if (!mesh->mMesh.HasFieldData(name))
            return fail(MIO_ERR_NOT_FOUND,
                        "meshio++: no field_data named '" + std::string(name) + "'");
        return array_out(mesh->mMesh.FieldData(name), data, dtype, ndim, shape);
    });
}

}  // extern "C"
