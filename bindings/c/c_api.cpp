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
 * @brief Implementation of the meshio++ C API (bindings/c/include/
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
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <array>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

// Project includes
#include "meshioplusplus/meshioplusplus.h"

#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/mdpa.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/mfem.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/formats/iges.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/formats/xdmf_time_series.hpp"
#include "meshioplusplus/formats/exodus.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/curvature.hpp"
#include "meshioplusplus/operations/normals.hpp"
#include "meshioplusplus/operations/repair.hpp"
#include "meshioplusplus/operations/shrinkwrap.hpp"
#include "meshioplusplus/operations/sobolev_deform.hpp"
#include "meshioplusplus/operations/crop.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/data_calc.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/data_condition.hpp"
#include "meshioplusplus/operations/tensor_invariants.hpp"
#include "meshioplusplus/operations/data_info.hpp"
#include "meshioplusplus/operations/data_integrate.hpp"
#include "meshioplusplus/operations/data_manage.hpp"
#include "meshioplusplus/formats/med.hpp"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/decimate_volume.hpp"
#include "meshioplusplus/operations/conservative_interpolate.hpp"
#include "meshioplusplus/operations/diff.hpp"
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/pipeline.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/tessellate.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/operations/error.hpp"
#include "meshioplusplus/operations/error_size.hpp"
#include "meshioplusplus/operations/check_conformity.hpp"
#include "meshioplusplus/operations/gradient.hpp"
#include "meshioplusplus/operations/hessian.hpp"
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/operations/agglomerate.hpp"
#include "meshioplusplus/operations/subdivide.hpp"
#include "meshioplusplus/operations/voxelize.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/octree.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"
#include "meshioplusplus/operations/octree_mesh.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/operations/blend.hpp"
#include "meshioplusplus/operations/feature_edges.hpp"
#include "meshioplusplus/operations/fit_cad.hpp"
#include "meshioplusplus/operations/refine_patches.hpp"
#include "meshioplusplus/operations/fit_primitives.hpp"
#include "meshioplusplus/operations/iga.hpp"
#include "meshioplusplus/operations/segment.hpp"
#include "meshioplusplus/operations/render.hpp"
#include "meshioplusplus/operations/interfaces.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/periodic.hpp"
#include "meshioplusplus/operations/quality_gate.hpp"
#include "meshioplusplus/operations/region_ops.hpp"
#include "meshioplusplus/operations/boundary_layer.hpp"
#include "meshioplusplus/operations/extrude.hpp"
#include "meshioplusplus/operations/mirror.hpp"
#include "meshioplusplus/operations/tetrahedralize.hpp"
#include "meshioplusplus/operations/project_to_cad.hpp"
#include "meshioplusplus/operations/undo_green.hpp"
#include "meshioplusplus/operations/mmg.hpp"
#include "meshioplusplus/operations/remesh.hpp"
#include "meshioplusplus/operations/remesh_volume.hpp"
#include "meshioplusplus/operations/optimize_volume.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/version.hpp"
#include "meshioplusplus/write_options.hpp"
#include "meshioplusplus/formats/gltf.hpp"
#include "meshioplusplus/skin.hpp"

struct mio_mesh {
    meshioplusplus::Mesh mMesh;
};

struct mio_regions {
    // A snapshot, not a borrow: the KRATOS backend serves Region(i) from lazily
    // rebuilt staging, so a pointer into the mesh would not stay valid.
    std::vector<meshioplusplus::Region> mRegions;
};

struct mio_periodic_pairs {
    meshioplusplus::PeriodicPairs mPairs;
};

struct mio_find_interface_result {
    mio_mesh mMesh;
    meshioplusplus::Region mSideA;
    meshioplusplus::Region mSideB;
    meshioplusplus::InterfaceReport mReport;
};

struct mio_contact_pairs_result {
    meshioplusplus::ContactPairsResult mPairs;
};

struct mio_poly_conn {
    // Likewise a snapshot rather than a rule-3 borrow, and for a stronger
    // reason: the MESHIO backend stores ragged blocks as nested vectors, so no
    // offsets array exists inside the mesh to point at, and only NativeMesh
    // keeps real CSR -- which c_api.cpp cannot reach, since it compiles under
    // all three backends. Built from the uniform API's RowSize/Row/NumFaces/
    // Face exactly as bindings/wasm/js_bindings.cpp does.
    bool mIsPolyhedron = false;
    std::int64_t mNumCells = 0;
    std::vector<std::int64_t> mNodes;
    std::vector<std::int64_t> mFaceOffsets;
    std::vector<std::int64_t> mCellOffsets;  // empty for a 1-level block
};

struct mio_reorder_result {
    mio_mesh mMesh;  // owns the renumbered mesh; borrowed via mio_reorder_result_mesh
    meshioplusplus::NDArray mNodePerm;
    std::vector<meshioplusplus::NDArray> mCellPerms;
};

struct mio_diff_result {
    meshioplusplus::DiffReport mReport;
};

struct mio_pipeline_report {
    std::string mJson;
};

struct mio_frame {
    meshioplusplus::Frame mFrame;
};

struct mio_xdmf_series {
    // Held by value: the writer is already move-only and owns its open
    // heavy-data container, so the handle is just the C-side name for it.
    meshioplusplus::XdmfTimeSeriesWriter mWriter;
};

struct mio_exodus_series {
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
    meshioplusplus::ExodusTimeSeriesWriter mWriter;
    explicit mio_exodus_series(const std::string& rPath) : mWriter(rPath) {}
#endif
};

/// The PLAN for a sequence -- paths, per-file step indices and times. It owns
/// no mesh, deliberately: caching what it handed out is exactly the
/// accumulation the streaming guarantee forbids, which is why
/// mio_sequence_read returns an OWNED mesh rather than a borrow.
struct mio_sequence {
    std::vector<meshioplusplus::SequenceEntry> mEntries;
    std::string mFormat;
    meshioplusplus::ReadOptions mOptions;
};

struct mio_split_result {
    std::vector<mio_mesh> mMeshes;   // owns each piece; borrowed via mio_split_result_mesh
    std::vector<std::string> mKeys;  // per-piece key (type / component / tag value)
};

struct mio_convert_cells_result {
    mio_mesh mMesh;  // owns the converted mesh; borrowed via _result_mesh
    meshioplusplus::NDArray mPointMap;
    std::vector<meshioplusplus::NDArray> mCellMaps;
};

struct mio_subdivide_result {
    mio_mesh mMesh;  // owns the subdivided mesh; borrowed via _result_mesh
    std::vector<meshioplusplus::NDArray> mCellMaps;
};

struct mio_agglomerate_result {
    mio_mesh mMesh;                    // owns the coarsened mesh; borrowed via _result_mesh
    meshioplusplus::NDArray mCellMap;  // flat, unlike mio_subdivide_result's per-block vector
};

struct mio_refine_result {
    mio_mesh mMesh;  // owns the refined mesh; borrowed via _result_mesh
    meshioplusplus::NDArray mPointMap;
    std::vector<meshioplusplus::NDArray> mCellMaps;
};

struct mio_decimate_result {
    mio_mesh mMesh;  // owns the decimated mesh; borrowed via _result_mesh
    meshioplusplus::NDArray mPointMap;
    std::vector<meshioplusplus::NDArray> mCellMaps;
    int64_t mFacesRemoved = 0;
    int64_t mPointsRemoved = 0;
    int64_t mCollapsesRejected = 0;
    double mMaxErrorApplied = 0.0;
};

struct mio_decimate_volume_result {
    mio_mesh mMesh;  // owns the decimated mesh; borrowed via _result_mesh
    meshioplusplus::NDArray mPointMap;
    std::vector<meshioplusplus::NDArray> mCellMaps;
    int64_t mTetsRemoved = 0;
    int64_t mPointsRemoved = 0;
    int64_t mCollapsesRejected = 0;
    double mMaxErrorApplied = 0.0;
};

struct mio_partition_result {
    struct Piece {
        mio_mesh mMesh;  // owns the piece; borrowed via _result_mesh
        meshioplusplus::NDArray mPointMap;
        std::vector<meshioplusplus::NDArray> mCellMaps;
    };
    std::vector<Piece> mPieces;  // exactly nparts entries, part id == index
};

struct mio_data_info {
    meshioplusplus::DataInfoReport mReport;
};

struct mio_data_integrate {
    meshioplusplus::DataIntegrateReport mReport;
};

/// Opaque handle behind mio_read_metadata_*; modelled on mio_data_info above,
/// since the summary is variable-length (cell blocks, name lists).
struct mio_read_metadata {
    meshioplusplus::MeshMetadata mMeta;
};

// The flattened lists the B-rep accessors index (loops, trims and topology
// references are nested per face, edge and vertex in BrepInfo), built once
// when a handle is made; the pointers are into the handle's own BrepInfo.
struct CapiBrepIndex {
    struct Loop {
        std::int64_t mFace;
        std::int64_t mFirstTrim;
        const meshioplusplus::BrepLoop* mpLoop;
    };
    struct Trim {
        std::int64_t mFace;
        std::int64_t mLoop;
        const meshioplusplus::BrepTrim* mpTrim;
    };
    std::vector<Loop> mLoops;
    std::vector<Trim> mTrims;
    std::vector<std::array<std::int64_t, 4>> mFaces;  // first loop, loops, first trim, trims
    std::vector<const meshioplusplus::BrepRef*> mRefs;
    std::vector<std::array<std::int64_t, 2>> mEdgeRefs, mVertexRefs;  // first, count
};

// One side channel, tagged by format. Later formats add alternatives; the
// accessors check the alternative they need, so a handle is never misread.
struct mio_format_info {
    std::string mFormat;
    std::variant<std::monostate, meshioplusplus::MdpaInfo, meshioplusplus::GmshInfo,
                 meshioplusplus::NurbsInfo, meshioplusplus::BrepInfo>
        mInfo;
    std::unique_ptr<CapiBrepIndex> mBrepIndex;
};

// What mio_iga_evaluate returns: the evaluation and whether it came from a B-rep face.
struct mio_iga_evaluation {
    meshioplusplus::IgaEvaluation mEval;
    bool mFace = false;
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
    } catch (const meshioplusplus::Unsupported& e) {
        return fail(MIO_ERR_UNSUPPORTED, e.what());
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

// Helpers for the analysis-and-editing entry points (v16.23.0).
namespace {

meshioplusplus::RegionSelector capi_selector(const mio_region_selector* pSel) {
    if (!pSel || !pSel->name)
        throw std::invalid_argument("meshio++: a region selector or its name is NULL");
    if (pSel->kind < -1 || pSel->kind > 2)
        throw std::invalid_argument(
            "meshio++: a region selector's kind must be -1 or a "
            "mio_region_kind");
    meshioplusplus::RegionSelector s;
    s.mName = pSel->name;
    s.mKind = pSel->kind;
    s.mDim = pSel->dim;
    s.mTag = pSel->tag;
    return s;
}

const int64_t* capi_periodic_ids(const meshioplusplus::NDArray& rIds, int64_t* pCount) {
    if (pCount)
        *pCount = static_cast<int64_t>(rIds.Size());
    return rIds.Size() == 0 ? nullptr : rIds.As<std::int64_t>();
}

}  // namespace

extern "C" {

/* ------------------------------------------------------------------ */
/* Version / build introspection                                       */
/* ------------------------------------------------------------------ */

// The C header repeats the release version as preprocessor macros so consumers
// can feature-detect; these pin them to the C++ ones, so the two copies cannot
// drift. (CMake separately hard-fails if either disagrees with the project
// version.) Same technique as the mio_cell_type / CellType static_asserts.
static_assert(MIO_VERSION_MAJOR == MESHIOPLUSPLUS_VERSION_MAJOR &&
                  MIO_VERSION_MINOR == MESHIOPLUSPLUS_VERSION_MINOR &&
                  MIO_VERSION_PATCH == MESHIOPLUSPLUS_VERSION_PATCH,
              "MIO_VERSION_* drifted from MESHIOPLUSPLUS_VERSION_*");

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

namespace {

/** @brief Resolve a read format, falling back to a content sniff. */
std::string capi_resolve_read_format(const char* pPath, const char* pFormat) {
    try {
        return meshioplusplus::resolve_format(pPath, format_or_empty(pFormat));
    } catch (const meshioplusplus::ReadError&) {
        // Extension gave nothing: fall back to a conservative content sniff.
        std::string sniffed = meshioplusplus::sniff_format(pPath);
        if (sniffed.empty())
            throw;
        return sniffed;
    }
}

/** @brief `mio_read_opts` -> `ReadOptions`; NULL means the defaults. */
meshioplusplus::ReadOptions capi_read_options(const mio_read_opts* pOpts) {
    meshioplusplus::ReadOptions out;
    if (!pOpts)
        return out;
    out.mPointsOnly = pOpts->points_only != 0;
    out.mMetadataOnly = pOpts->metadata_only != 0;
    // A NULL array pointer means "every array"; a non-NULL pointer with count 0
    // means "no arrays". Collapsing the two here would silently turn an
    // explicit "none" into "everything".
    if (pOpts->arrays) {
        std::vector<std::string> names;
        names.reserve(static_cast<std::size_t>(pOpts->num_arrays));
        for (std::int64_t i = 0; i < pOpts->num_arrays; ++i)
            names.emplace_back(pOpts->arrays[i] ? pOpts->arrays[i] : "");
        out.mDataArrays = std::move(names);  // setters copy, per the ABI contract
    }
    switch (pOpts->mmap_mode) {
        case 1:
            out.mMmap = meshioplusplus::MmapMode::On;
            break;
        case 2:
            out.mMmap = meshioplusplus::MmapMode::Off;
            break;
        default:
            out.mMmap = meshioplusplus::MmapMode::Auto;
            break;
    }
    // int64 on the ABI (the reserved slots are int64) but int in the core: a
    // step index beyond int range cannot name a real step, so reject it here
    // rather than let the narrowing wrap into a valid-looking one.
    if (pOpts->time_step > INT_MAX || pOpts->time_step < INT_MIN)
        throw meshioplusplus::ReadError("meshio++: time_step is out of range");
    out.mTimeStep = static_cast<int>(pOpts->time_step);
    out.mLenient = pOpts->lenient != 0;
    out.mPiece = pOpts->piece;
    out.mPieceSet = pOpts->piece_set != 0;
    out.mGhosts = pOpts->drop_ghosts != 0 ? meshioplusplus::GhostPolicy::Drop
                                          : meshioplusplus::GhostPolicy::Keep;
    return out;
}

/** @brief The name list for a `mio_data_location`. */
const std::vector<std::string>& capi_metadata_names(const meshioplusplus::MeshMetadata& rMeta,
                                                    int location) {
    switch (location) {
        case MIO_DATA_POINT:
            return rMeta.mPointDataNames;
        case MIO_DATA_CELL:
            return rMeta.mCellDataNames;
        case MIO_DATA_FIELD:
            return rMeta.mFieldDataNames;
        default:
            throw meshioplusplus::ReadError("meshio++: unknown data location");
    }
}

}  // namespace

// The layout the Fortran (`mio_read_opts_t`) and Julia (`_CReadOpts`, which also
// checks it at load time) mirrors are written against: every field added since
// the struct shipped took a former `reserved` slot, so it has never moved.
static_assert(sizeof(mio_read_opts) == 80, "mio_read_opts grew outside its reserved tail");

void mio_read_opts_init(mio_read_opts* opts) {
    if (!opts)
        return;
    *opts = mio_read_opts{};  // value-initialized: all zero == read everything
}

mio_mesh* mio_read_ex(const char* path, const char* format, const mio_read_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!path)
            throw meshioplusplus::ReadError("meshio++: path is NULL");
        const std::string fmt = capi_resolve_read_format(path, format);
        if (!meshioplusplus::registry_readers().count(fmt) &&
            !meshioplusplus::registry_reader_supports_options(fmt))
            throw meshioplusplus::ReadError(unknown_format_message(fmt, /*for_write=*/false));
        return new mio_mesh{meshioplusplus::registry_read(path, fmt, capi_read_options(opts))};
    });
}

mio_read_metadata* mio_read_metadata_create(const char* path, const char* format) {
    return guarded_ptr(static_cast<mio_read_metadata*>(nullptr), [&]() -> mio_read_metadata* {
        if (!path)
            throw meshioplusplus::ReadError("meshio++: path is NULL");
        const std::string fmt = capi_resolve_read_format(path, format);
        return new mio_read_metadata{
            meshioplusplus::registry_read_metadata(path, fmt, meshioplusplus::ReadOptions{})};
    });
}

int64_t mio_read_metadata_num_points(const mio_read_metadata* meta) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        return static_cast<std::int64_t>(meta->mMeta.mNumPoints);
    });
}

int64_t mio_read_metadata_point_dim(const mio_read_metadata* meta) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        return static_cast<std::int64_t>(meta->mMeta.mPointDim);
    });
}

int64_t mio_read_metadata_num_cells(const mio_read_metadata* meta) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        return static_cast<std::int64_t>(meta->mMeta.NumCells());
    });
}

int64_t mio_read_metadata_num_cell_blocks(const mio_read_metadata* meta) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        return static_cast<std::int64_t>(meta->mMeta.mCellBlocks.size());
    });
}

mio_status mio_read_metadata_cell_block(const mio_read_metadata* meta, int64_t index,
                                        int64_t* num_cells, int64_t* num_nodes_per_cell,
                                        int* is_ragged) {
    return guarded([&]() -> mio_status {
        if (!meta)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: metadata handle is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= meta->mMeta.mCellBlocks.size())
            return fail(MIO_ERR_INVALID_ARG, "meshio++: cell block index out of range");
        const meshioplusplus::CellBlockInfo& block =
            meta->mMeta.mCellBlocks[static_cast<std::size_t>(index)];
        if (num_cells)
            *num_cells = static_cast<std::int64_t>(block.mNumCells);
        if (num_nodes_per_cell)
            *num_nodes_per_cell = static_cast<std::int64_t>(block.mNodesPerCell);
        if (is_ragged)
            *is_ragged = block.mRagged ? 1 : 0;
        return MIO_OK;
    });
}

int64_t mio_read_metadata_cell_block_type(const mio_read_metadata* meta, int64_t index, char* buf,
                                          int64_t buflen) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= meta->mMeta.mCellBlocks.size())
            throw meshioplusplus::ReadError("meshio++: cell block index out of range");
        return copy_string(meta->mMeta.mCellBlocks[static_cast<std::size_t>(index)].mType, buf,
                           buflen);
    });
}

int64_t mio_read_metadata_num_time_values(const mio_read_metadata* meta) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        return static_cast<std::int64_t>(meta->mMeta.mTimeValues.size());
    });
}

int64_t mio_read_metadata_time_values(const mio_read_metadata* meta, double* out, int64_t count) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        const std::vector<double>& values = meta->mMeta.mTimeValues;
        const std::int64_t n =
            std::min<std::int64_t>(count < 0 ? 0 : count, static_cast<std::int64_t>(values.size()));
        if (out)
            for (std::int64_t i = 0; i < n; ++i)
                out[i] = values[static_cast<std::size_t>(i)];
        return n;
    });
}

int64_t mio_read_metadata_num_regions(const mio_read_metadata* meta) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        return static_cast<std::int64_t>(meta->mMeta.mRegions.size());
    });
}

int64_t mio_read_metadata_region_name(const mio_read_metadata* meta, int64_t index, char* buf,
                                      int64_t buflen) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        const std::vector<meshioplusplus::RegionSummary>& regions = meta->mMeta.mRegions;
        if (index < 0 || static_cast<std::size_t>(index) >= regions.size())
            throw meshioplusplus::ReadError("meshio++: region index " + std::to_string(index) +
                                            " out of range");
        return copy_string(regions[static_cast<std::size_t>(index)].mName, buf, buflen);
    });
}

mio_status mio_read_metadata_region_info(const mio_read_metadata* meta, int64_t index,
                                         mio_region_info* out) {
    return guarded([&]() -> mio_status {
        if (!meta || !out)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: metadata/out is NULL");
        const std::vector<meshioplusplus::RegionSummary>& regions = meta->mMeta.mRegions;
        if (index < 0 || static_cast<std::size_t>(index) >= regions.size())
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: region index " + std::to_string(index) + " out of range");
        const meshioplusplus::RegionSummary& r = regions[static_cast<std::size_t>(index)];
        out->kind = static_cast<int32_t>(r.mKind);
        out->dim = static_cast<int32_t>(r.mDim);
        out->tag = r.mTag;
        out->num_entries = static_cast<int64_t>(r.mNumEntries);
        out->stride = 0;  // no entries are carried by a summary
        return MIO_OK;
    });
}

int64_t mio_read_metadata_num_names(const mio_read_metadata* meta, int location) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        return static_cast<std::int64_t>(capi_metadata_names(meta->mMeta, location).size());
    });
}

int64_t mio_read_metadata_name(const mio_read_metadata* meta, int location, int64_t index,
                               char* buf, int64_t buflen) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        const std::vector<std::string>& names = capi_metadata_names(meta->mMeta, location);
        if (index < 0 || static_cast<std::size_t>(index) >= names.size())
            throw meshioplusplus::ReadError("meshio++: name index out of range");
        return copy_string(names[static_cast<std::size_t>(index)], buf, buflen);
    });
}

mio_status mio_read_metadata_bbox(const mio_read_metadata* meta, double* min, double* max) {
    return guarded([&]() -> mio_status {
        if (!meta)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: metadata handle is NULL");
        if (!meta->mMeta.mHasBBox)
            return fail(MIO_ERR_NOT_FOUND,
                        "meshio++: no bounding box in this summary (a native metadata path "
                        "does not decode the point coordinates)");
        for (int d = 0; d < 3; ++d) {
            if (min)
                min[d] = meta->mMeta.mBBoxMin[d];
            if (max)
                max[d] = meta->mMeta.mBBoxMax[d];
        }
        return MIO_OK;
    });
}

int mio_read_metadata_fell_back(const mio_read_metadata* meta) {
    return guarded_ptr(-1, [&]() -> int {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        return meta->mMeta.mFellBackToFullRead ? 1 : 0;
    });
}

int64_t mio_read_metadata_num_provenance_lines(const mio_read_metadata* meta) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        return static_cast<std::int64_t>(meta->mMeta.mProvenance.size());
    });
}

int64_t mio_read_metadata_provenance_line(const mio_read_metadata* meta, int64_t index, char* out,
                                          int64_t cap) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        const std::vector<std::string>& lines = meta->mMeta.mProvenance;
        if (index < 0 || static_cast<std::size_t>(index) >= lines.size())
            throw meshioplusplus::ReadError("meshio++: provenance line index " +
                                            std::to_string(index) + " out of range");
        return copy_string(lines[static_cast<std::size_t>(index)], out, cap);
    });
}

int mio_read_metadata_provenance_recognised(const mio_read_metadata* meta) {
    return guarded_ptr(-1, [&]() -> int {
        if (!meta)
            throw meshioplusplus::ReadError("meshio++: metadata handle is NULL");
        return meta->mMeta.mProvenanceRecognised ? 1 : 0;
    });
}

void mio_read_metadata_free(mio_read_metadata* meta) {
    delete meta;
}

int mio_reader_supports_options(const char* format) {
    return guarded_ptr(-1, [&]() -> int {
        if (!format)
            throw meshioplusplus::ReadError("meshio++: format is NULL");
        return meshioplusplus::registry_reader_supports_options(format) ? 1 : 0;
    });
}

mio_status mio_write(const char* path, const mio_mesh* mesh, const char* format) {
    return guarded([&]() -> mio_status {
        if (!path || !mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: path/mesh is NULL");
        std::string fmt = meshioplusplus::resolve_write_format(path, format_or_empty(format));
        meshioplusplus::detail::provenance_begin_write();
        auto it = meshioplusplus::registry_writers().find(fmt);
        if (it == meshioplusplus::registry_writers().end())
            return fail(MIO_ERR_NOT_FOUND, unknown_format_message(fmt, /*for_write=*/true));
        it->second(path, mesh->mMesh);
        return MIO_OK;
    });
}

void mio_write_opts_init(mio_write_opts* opts) {
    if (!opts)
        return;
    *opts = mio_write_opts{};  // value-initialized: all zero == mio_write()
}

void mio_gltf_opts_init(mio_gltf_opts* opts) {
    if (!opts)
        return;
    *opts = mio_gltf_opts{};
    opts->normals = opts->fields = opts->recenter = opts->by_region = opts->unlit = 1;
    opts->split_angle = 30.0;
    opts->scale = 1.0;
}

mio_status mio_write_gltf(const char* path, const mio_mesh* mesh, const mio_gltf_opts* opts) {
    return guarded([&]() -> mio_status {
        if (!path || !mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: path/mesh is NULL");
        meshioplusplus::GltfWriteOptions out;
        if (opts) {
            if (opts->container < 0 || opts->container > 2 || opts->up_axis < 0 ||
                opts->up_axis > 3 || opts->normal_weight < 0 || opts->normal_weight > 1)
                return fail(MIO_ERR_INVALID_ARG,
                            "meshio++: invalid glTF container, axis or weight");
            out.mContainer = static_cast<meshioplusplus::GltfContainer>(opts->container);
            out.mUpAxis = static_cast<meshioplusplus::GltfUpAxis>(opts->up_axis);
            out.mNormalWeight =
                static_cast<meshioplusplus::SdfPseudonormalWeight>(opts->normal_weight);
            out.mNormals = opts->normals != 0;
            out.mFields = opts->fields != 0;
            out.mRecenter = opts->recenter != 0;
            out.mByRegion = opts->by_region != 0;
            out.mUnlit = opts->unlit != 0;
            out.mSplitAngle = opts->split_angle;
            out.mScale = opts->scale;
            if (opts->component_set)
                out.mComponent = opts->component;
            if (opts->vmin_set)
                out.mVMin = opts->vmin;
            if (opts->vmax_set)
                out.mVMax = opts->vmax;
            if (opts->color_by)
                out.mColorBy = opts->color_by;
            if (opts->cmap)
                out.mCmap = opts->cmap;
            if (opts->nan_color)
                out.mNanColor = opts->nan_color;
        }
        meshioplusplus::detail::provenance_begin_write();
        meshioplusplus::write_gltf(path, mesh->mMesh, out);
        return MIO_OK;
    });
}

static_assert(sizeof(mio_gltf_opts) == 168, "mio_gltf_opts layout changed");

namespace {

/// mio_write_opts -> WriteOptions, shared by mio_write_ex and
/// mio_sequence_to_timeseries_ex. Returns MIO_ERR_INVALID_ARG (via fail(),
/// which sets the thread-local message) on a bad enum value; MIO_OK otherwise.
mio_status write_opts_to_cxx(const mio_write_opts& rOpts, meshioplusplus::WriteOptions& rOut) {
    switch (rOpts.encoding) {
        case MIO_ENCODING_DEFAULT:
            break;
        case MIO_ENCODING_ASCII:
            rOut.mEncoding = meshioplusplus::WriteEncoding::Ascii;
            break;
        case MIO_ENCODING_BINARY:
            rOut.mEncoding = meshioplusplus::WriteEncoding::Binary;
            break;
        case MIO_ENCODING_RAW_APPENDED:
            rOut.mEncoding = meshioplusplus::WriteEncoding::RawAppended;
            break;
        default:
            return fail(MIO_ERR_INVALID_ARG, "meshio++: bad mio_write_opts.encoding");
    }
    switch (rOpts.codec) {
        case MIO_CODEC_DEFAULT:
            break;
        case MIO_CODEC_NONE:
            rOut.mCodec = meshioplusplus::detail::VtkCodec::None;
            rOut.mCodecSet = true;
            break;
        case MIO_CODEC_ZLIB:
            rOut.mCodec = meshioplusplus::detail::VtkCodec::Zlib;
            rOut.mCodecSet = true;
            break;
        case MIO_CODEC_LZ4:
            rOut.mCodec = meshioplusplus::detail::VtkCodec::LZ4;
            rOut.mCodecSet = true;
            break;
        case MIO_CODEC_ZSTD:
            rOut.mCodec = meshioplusplus::detail::VtkCodec::ZSTD;
            rOut.mCodecSet = true;
            break;
        case MIO_CODEC_LZF:
            rOut.mCodec = meshioplusplus::detail::VtkCodec::LZF;
            rOut.mCodecSet = true;
            break;
        default:
            return fail(MIO_ERR_INVALID_ARG, "meshio++: bad mio_write_opts.codec");
    }
    if (rOpts.float_format)
        rOut.mFloatFormat = rOpts.float_format;
    return MIO_OK;
}

}  // namespace

mio_status mio_write_ex(const char* path, const mio_mesh* mesh, const char* format,
                        const mio_write_opts* opts) {
    return guarded([&]() -> mio_status {
        if (!path || !mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: path/mesh is NULL");
        if (!opts)
            return mio_write(path, mesh, format);

        meshioplusplus::WriteOptions w;
        const mio_status opt_status = write_opts_to_cxx(*opts, w);
        if (opt_status != MIO_OK)
            return opt_status;

        meshioplusplus::registry_write_ex(path, mesh->mMesh, format_or_empty(format), w);
        return MIO_OK;
    });
}

namespace {

// Flatten a B-rep handle's nested lists for the accessors (see CapiBrepIndex).
void capi_index_brep(mio_format_info& rHandle) {
    const auto& b = std::get<meshioplusplus::BrepInfo>(rHandle.mInfo);
    auto idx = std::make_unique<CapiBrepIndex>();
    for (std::size_t f = 0; f < b.mFaces.size(); ++f) {
        const auto& face = b.mFaces[f];
        const auto fi = static_cast<std::int64_t>(f);
        const auto first_loop = static_cast<std::int64_t>(idx->mLoops.size());
        const auto first_trim = static_cast<std::int64_t>(idx->mTrims.size());
        for (std::size_t l = 0; l < face.mLoops.size(); ++l) {
            idx->mLoops.push_back(
                {fi, static_cast<std::int64_t>(idx->mTrims.size()), &face.mLoops[l]});
            for (const auto& t : face.mLoops[l].mTrims)
                idx->mTrims.push_back({fi, static_cast<std::int64_t>(idx->mLoops.size() - 1), &t});
        }
        for (const auto& t : face.mEmbeddedEdges)
            idx->mTrims.push_back({fi, -1, &t});
        idx->mFaces.push_back({first_loop, static_cast<std::int64_t>(face.mLoops.size()),
                               first_trim,
                               static_cast<std::int64_t>(idx->mTrims.size()) - first_trim});
    }
    for (const auto& e : b.mEdges) {
        idx->mEdgeRefs.push_back({static_cast<std::int64_t>(idx->mRefs.size()),
                                  static_cast<std::int64_t>(e.mTopology.size())});
        for (const auto& r : e.mTopology)
            idx->mRefs.push_back(&r);
    }
    for (const auto& v : b.mVertices) {
        idx->mVertexRefs.push_back({static_cast<std::int64_t>(idx->mRefs.size()),
                                    static_cast<std::int64_t>(v.mTopology.size())});
        for (const auto& r : v.mTopology)
            idx->mRefs.push_back(&r);
    }
    rHandle.mBrepIndex = std::move(idx);
}

}  // namespace

mio_mesh* mio_read_with_info(const char* path, const char* format, const mio_read_opts* opts,
                             mio_format_info** info) {
    if (info)
        *info = nullptr;
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!path || !info)
            throw std::invalid_argument("meshio++: path/info is NULL");
        const std::string fmt = capi_resolve_read_format(path, format);
        if (fmt == "gmsh") {
            auto handle = std::make_unique<mio_format_info>();
            handle->mFormat = fmt;
            meshioplusplus::GmshInfo gmsh;
            auto mesh = std::make_unique<mio_mesh>(
                mio_mesh{meshioplusplus::read_gmsh(path, gmsh, capi_read_options(opts))});
            handle->mInfo = std::move(gmsh);
            *info = handle.release();
            return mesh.release();
        }
        if (fmt == "mfem") {
            auto handle = std::make_unique<mio_format_info>();
            handle->mFormat = fmt;
            meshioplusplus::NurbsInfo nurbs;
            auto mesh = std::make_unique<mio_mesh>(
                mio_mesh{meshioplusplus::read_mfem(path, {}, capi_read_options(opts), nurbs)});
            handle->mInfo = std::move(nurbs);
            *info = handle.release();
            return mesh.release();
        }
        if (fmt == "cad_json" || fmt == "step" || fmt == "iges") {
            auto handle = std::make_unique<mio_format_info>();
            handle->mFormat = fmt;
            meshioplusplus::BrepInfo brep;
            auto mesh = std::make_unique<mio_mesh>(
                mio_mesh{fmt == "cad_json" ? meshioplusplus::read_cad_json(path, brep)
                         : fmt == "step"   ? meshioplusplus::read_step(path, brep)
                                           : meshioplusplus::read_iges(path, brep)});
            handle->mInfo = std::move(brep);
            capi_index_brep(*handle);
            *info = handle.release();
            return mesh.release();
        }
        if (fmt == "mdpa") {
            auto handle = std::make_unique<mio_format_info>();
            handle->mFormat = fmt;
            meshioplusplus::MdpaInfo mdpa;
            auto mesh = std::make_unique<mio_mesh>(
                mio_mesh{meshioplusplus::read_mdpa(path, mdpa, capi_read_options(opts))});
            handle->mInfo = std::move(mdpa);
            *info = handle.release();
            return mesh.release();
        }
        if (!meshioplusplus::registry_readers().count(fmt) &&
            !meshioplusplus::registry_reader_supports_options(fmt))
            throw meshioplusplus::ReadError(unknown_format_message(fmt, /*for_write=*/false));
        return new mio_mesh{meshioplusplus::registry_read(path, fmt, capi_read_options(opts))};
    });
}

mio_status mio_write_with_info(const char* path, const mio_mesh* mesh, const char* format,
                               const mio_format_info* info) {
    return guarded([&]() -> mio_status {
        if (!path || !mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: path/mesh is NULL");
        if (!info)
            return mio_write(path, mesh, format);
        const std::string fmt = meshioplusplus::resolve_write_format(path, format_or_empty(format));
        const bool brep_as_cad_json =
            fmt == "cad_json" && std::holds_alternative<meshioplusplus::BrepInfo>(info->mInfo);
        if (fmt != info->mFormat && !(fmt == "gmsh22" && info->mFormat == "gmsh") &&
            !brep_as_cad_json)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: a '" + info->mFormat +
                                                 "' side channel cannot be written as '" + fmt +
                                                 "'");
        meshioplusplus::detail::provenance_begin_write();
        if (const auto* p_gmsh = std::get_if<meshioplusplus::GmshInfo>(&info->mInfo)) {
            if (fmt == "gmsh22")
                meshioplusplus::write_gmsh22(path, mesh->mMesh, true, *p_gmsh);
            else
                meshioplusplus::write_gmsh41(path, mesh->mMesh, true, *p_gmsh);
            return MIO_OK;
        }
        if (const auto* p_mdpa = std::get_if<meshioplusplus::MdpaInfo>(&info->mInfo)) {
            meshioplusplus::write_mdpa(path, mesh->mMesh, *p_mdpa);
            return MIO_OK;
        }
        if (const auto* p_nurbs = std::get_if<meshioplusplus::NurbsInfo>(&info->mInfo)) {
            meshioplusplus::write_mfem(path, mesh->mMesh, false, *p_nurbs);
            return MIO_OK;
        }
        if (const auto* p_brep = std::get_if<meshioplusplus::BrepInfo>(&info->mInfo)) {
            if (fmt == "cad_json") {
                meshioplusplus::write_cad_json(path, mesh->mMesh, *p_brep);
                return MIO_OK;
            }
        }
        return mio_write(path, mesh, format);
    });
}

int64_t mio_format_info_format(const mio_format_info* info, char* buf, int64_t buflen) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!info)
            throw std::invalid_argument("meshio++: format info handle is NULL");
        return copy_string(info->mFormat, buf, buflen);
    });
}

void mio_format_info_free(mio_format_info* info) {
    delete info;
}

namespace {
const meshioplusplus::GmshInfo& capi_gmsh(const mio_format_info* pInfo) {
    if (!pInfo)
        throw std::invalid_argument("meshio++: format info handle is NULL");
    const auto* info = std::get_if<meshioplusplus::GmshInfo>(&pInfo->mInfo);
    if (!info)
        throw std::invalid_argument("meshio++: not a gmsh side channel");
    return *info;
}
std::size_t capi_gmsh_count(const meshioplusplus::GmshInfo& rInfo, int section) {
    if (section == MIO_GMSH_BOUNDING_ENTITIES)
        return rInfo.mBoundingEntities.size();
    if (section == MIO_GMSH_PERIODIC)
        return rInfo.mPeriodic.size();
    throw std::invalid_argument("meshio++: unknown gmsh side-channel section");
}
}  // namespace

int64_t mio_gmsh_info_count(const mio_format_info* info, int32_t section) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        return static_cast<std::int64_t>(capi_gmsh_count(capi_gmsh(info), section));
    });
}

mio_status mio_gmsh_info_array(const mio_format_info* info, int32_t section, int64_t index,
                               int32_t field, const void** data, mio_dtype* dtype, int32_t* ndim,
                               int64_t* shape) {
    return guarded([&]() -> mio_status {
        const auto& gmsh = capi_gmsh(info);
        if (index < 0 || static_cast<std::size_t>(index) >= capi_gmsh_count(gmsh, section))
            return fail(MIO_ERR_INVALID_ARG, "meshio++: gmsh side-channel index out of range");
        if (section == MIO_GMSH_BOUNDING_ENTITIES && field == 0) {
            const auto& tags = gmsh.mBoundingEntities[static_cast<std::size_t>(index)];
            return array_out(
                meshioplusplus::NDArray::MakeView(
                    meshioplusplus::DType::Int32, {tags.size()},
                    reinterpret_cast<std::byte*>(const_cast<std::int32_t*>(tags.data()))),
                data, dtype, ndim, shape);
        }
        if (section == MIO_GMSH_PERIODIC) {
            const auto& link = gmsh.mPeriodic[static_cast<std::size_t>(index)];
            if (field == 0)
                return array_out(meshioplusplus::NDArray::MakeView(
                                     meshioplusplus::DType::Int32, {3},
                                     reinterpret_cast<std::byte*>(
                                         const_cast<std::int32_t*>(link.mEntityTags.data()))),
                                 data, dtype, ndim, shape);
            if (field == 1)
                return array_out(
                    meshioplusplus::NDArray::MakeView(
                        meshioplusplus::DType::Float64, {link.mAffine.size()},
                        reinterpret_cast<std::byte*>(const_cast<double*>(link.mAffine.data()))),
                    data, dtype, ndim, shape);
            if (field == 2)
                return array_out(link.mNodePairs, data, dtype, ndim, shape);
        }
        return fail(MIO_ERR_INVALID_ARG, "meshio++: unknown gmsh side-channel array field");
    });
}

namespace {

/// The NurbsInfo behind a handle, or a throw naming what it is instead.
const meshioplusplus::NurbsInfo& capi_nurbs(const mio_format_info* pInfo) {
    if (!pInfo)
        throw std::invalid_argument("meshio++: format info handle is NULL");
    const auto* info = std::get_if<meshioplusplus::NurbsInfo>(&pInfo->mInfo);
    if (!info)
        throw std::invalid_argument("meshio++: not an mfem (NURBS) side channel");
    return *info;
}

std::size_t capi_nurbs_count(const meshioplusplus::NurbsInfo& rInfo, int Section) {
    switch (Section) {
        case MIO_NURBS_PATCHES:
            return rInfo.mPatches.size();
        case MIO_NURBS_FIELDS:
            return rInfo.mFieldNames.size();
        case MIO_NURBS_TOPOLOGY:
            return rInfo.mTopology.mDimension ? 1 : 0;
        case MIO_NURBS_CELL_MAP:
            return rInfo.mPatches.empty() ? 0 : 1;
        default:
            throw std::invalid_argument("meshio++: unknown NURBS side-channel section");
    }
}

bool capi_nurbs_in_range(const meshioplusplus::NurbsInfo& rInfo, int Section, std::int64_t Index) {
    return Section >= MIO_NURBS_PATCHES && Section <= MIO_NURBS_CELL_MAP && Index >= 0 &&
           static_cast<std::size_t>(Index) < capi_nurbs_count(rInfo, Section);
}

std::string capi_nurbs_bad_field(int Section, int Field) {
    return "meshio++: NURBS side-channel section " + std::to_string(Section) +
           " has no such field " + std::to_string(Field);
}

}  // namespace

int64_t mio_nurbs_info_count(const mio_format_info* info, int32_t section) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        return static_cast<std::int64_t>(capi_nurbs_count(capi_nurbs(info), section));
    });
}

mio_status mio_nurbs_info_int(const mio_format_info* info, int32_t section, int64_t index,
                              int32_t field, int64_t* value) {
    return guarded([&]() -> mio_status {
        const meshioplusplus::NurbsInfo& r_info = capi_nurbs(info);
        if (!capi_nurbs_in_range(r_info, section, index))
            return fail(MIO_ERR_INVALID_ARG, "meshio++: NURBS side-channel index out of range");
        const auto i = static_cast<std::size_t>(index);
        std::int64_t v = 0;
        if (section == MIO_NURBS_PATCHES && field == 0)
            v = static_cast<std::int64_t>(r_info.mPatches[i].mDegrees.size());
        else if (section == MIO_NURBS_PATCHES && field == 1)
            v = r_info.mPatches[i].mControlPoints.Shape().size() == 2
                    ? static_cast<std::int64_t>(r_info.mPatches[i].mControlPoints.Shape()[1])
                    : 0;
        else if (section == MIO_NURBS_TOPOLOGY && field == 0)
            v = r_info.mTopology.mDimension;
        else if (section == MIO_NURBS_TOPOLOGY && field == 1)
            v = r_info.mTopology.mNumVertices;
        else
            return fail(MIO_ERR_INVALID_ARG, capi_nurbs_bad_field(section, field));
        if (value)
            *value = v;
        return MIO_OK;
    });
}

mio_status mio_nurbs_info_array(const mio_format_info* info, int32_t section, int64_t index,
                                int32_t field, const void** data, mio_dtype* dtype, int32_t* ndim,
                                int64_t* shape) {
    return guarded([&]() -> mio_status {
        const meshioplusplus::NurbsInfo& r_info = capi_nurbs(info);
        if (!capi_nurbs_in_range(r_info, section, index))
            return fail(MIO_ERR_INVALID_ARG, "meshio++: NURBS side-channel index out of range");
        const auto i = static_cast<std::size_t>(index);
        if (section == MIO_NURBS_PATCHES) {
            const meshioplusplus::NurbsPatch& p = r_info.mPatches[i];
            if (field == 0)
                return array_out(
                    meshioplusplus::NDArray::MakeView(
                        meshioplusplus::DType::Int32, {p.mDegrees.size()},
                        reinterpret_cast<std::byte*>(const_cast<std::int32_t*>(p.mDegrees.data()))),
                    data, dtype, ndim, shape);
            if (field == 1)
                return array_out(p.mControlPoints, data, dtype, ndim, shape);
            if (field == 2)
                return array_out(p.mWeights, data, dtype, ndim, shape);
            if (field >= 3 && field < 6 && static_cast<std::size_t>(field - 3) < p.mKnots.size()) {
                const std::vector<double>& k = p.mKnots[static_cast<std::size_t>(field - 3)];
                return array_out(meshioplusplus::NDArray::MakeView(
                                     meshioplusplus::DType::Float64, {k.size()},
                                     reinterpret_cast<std::byte*>(const_cast<double*>(k.data()))),
                                 data, dtype, ndim, shape);
            }
            if (field >= 6 && static_cast<std::size_t>(field - 6) < p.mFields.size())
                return array_out(p.mFields[static_cast<std::size_t>(field - 6)], data, dtype, ndim,
                                 shape);
        } else if (section == MIO_NURBS_TOPOLOGY) {
            const meshioplusplus::NurbsTopology& t = r_info.mTopology;
            const meshioplusplus::NDArray* arrays[] = {&t.mPatchVertices, &t.mPatchAttributes,
                                                       &t.mBoundaryVertices, &t.mBoundaryAttributes,
                                                       &t.mEdges};
            if (field >= 0 && field < 5)
                return array_out(*arrays[field], data, dtype, ndim, shape);
        } else if (section == MIO_NURBS_CELL_MAP) {
            const meshioplusplus::NDArray* arrays[] = {&r_info.mCellPatch, &r_info.mCellSpan,
                                                       &r_info.mBoundaryCellPatch};
            if (field >= 0 && field < 3)
                return array_out(*arrays[field], data, dtype, ndim, shape);
        }
        return fail(MIO_ERR_INVALID_ARG, capi_nurbs_bad_field(section, field));
    });
}

int64_t mio_nurbs_info_string(const mio_format_info* info, int32_t section, int64_t index,
                              int32_t field, char* buf, int64_t buflen) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        const meshioplusplus::NurbsInfo& r_info = capi_nurbs(info);
        if (!capi_nurbs_in_range(r_info, section, index))
            throw std::invalid_argument("meshio++: NURBS side-channel index out of range");
        if (section == MIO_NURBS_FIELDS && field == 0)
            return copy_string(r_info.mFieldNames[static_cast<std::size_t>(index)], buf, buflen);
        throw std::invalid_argument(capi_nurbs_bad_field(section, field));
    });
}

namespace {

/// The BrepInfo (and its index) behind a handle, or a throw naming what it is.
const meshioplusplus::BrepInfo& capi_brep(const mio_format_info* pInfo) {
    if (!pInfo)
        throw std::invalid_argument("meshio++: format info handle is NULL");
    const auto* info = std::get_if<meshioplusplus::BrepInfo>(&pInfo->mInfo);
    if (!info || !pInfo->mBrepIndex)
        throw std::invalid_argument("meshio++: not a B-rep side channel ('" + pInfo->mFormat +
                                    "')");
    return *info;
}

std::size_t capi_brep_count(const mio_format_info* pInfo, int Section) {
    const meshioplusplus::BrepInfo& b = capi_brep(pInfo);
    const CapiBrepIndex& idx = *pInfo->mBrepIndex;
    switch (Section) {
        case MIO_BREP_FILE:
            return 1;
        case MIO_BREP_BREPS:
            return b.mBreps.size();
        case MIO_BREP_PATCHES:
            return b.mPatches.size();
        case MIO_BREP_FACES:
            return b.mFaces.size();
        case MIO_BREP_LOOPS:
            return idx.mLoops.size();
        case MIO_BREP_TRIMS:
            return idx.mTrims.size();
        case MIO_BREP_EDGES:
            return b.mEdges.size();
        case MIO_BREP_VERTICES:
            return b.mVertices.size();
        case MIO_BREP_REFS:
            return idx.mRefs.size();
        default:
            throw std::invalid_argument("meshio++: unknown B-rep side-channel section");
    }
}

void capi_brep_check(const mio_format_info* pInfo, int Section, std::int64_t Index) {
    if (Index < 0 || static_cast<std::size_t>(Index) >= capi_brep_count(pInfo, Section))
        throw std::invalid_argument("meshio++: B-rep side-channel index out of range");
}

[[noreturn]] void capi_brep_bad_field(int Section, int Field) {
    throw std::invalid_argument("meshio++: B-rep side-channel section " + std::to_string(Section) +
                                " has no such field " + std::to_string(Field));
}

/// The curve an item's curve fields describe (a trim's parameter curve, an
/// edge's 3-D curve), or null.
const meshioplusplus::BrepCurve* capi_brep_curve(const mio_format_info* pInfo, int Section,
                                                 std::size_t Index) {
    const meshioplusplus::BrepInfo& b = capi_brep(pInfo);
    if (Section == MIO_BREP_TRIMS)
        return &pInfo->mBrepIndex->mTrims[Index].mpTrim->mCurve;
    if (Section == MIO_BREP_EDGES && b.mEdges[Index].mHasCurve)
        return &b.mEdges[Index].mCurve;
    return nullptr;
}

meshioplusplus::NDArray capi_doubles_view(const std::vector<double>& rV) {
    return meshioplusplus::NDArray::MakeView(
        meshioplusplus::DType::Float64, {rV.size()},
        reinterpret_cast<std::byte*>(const_cast<double*>(rV.data())));
}

meshioplusplus::NDArray capi_int64s_view(const std::vector<std::int64_t>& rV) {
    return meshioplusplus::NDArray::MakeView(
        meshioplusplus::DType::Int64, {rV.size()},
        reinterpret_cast<std::byte*>(const_cast<std::int64_t*>(rV.data())));
}

// guarded(), with a bad handle, section, index or field reported as
// MIO_ERR_INVALID_ARG rather than an internal error.
extern "C++" {
template <class F>
mio_status capi_brep_guarded(F&& f) {
    return guarded([&]() -> mio_status {
        try {
            return f();
        } catch (const std::invalid_argument& e) {
            return fail(MIO_ERR_INVALID_ARG, e.what());
        }
    });
}
}  // extern "C++"

}  // namespace

int64_t mio_brep_info_count(const mio_format_info* info, int32_t section) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        return static_cast<std::int64_t>(capi_brep_count(info, section));
    });
}

mio_status mio_brep_info_int(const mio_format_info* info, int32_t section, int64_t index,
                             int32_t field, int64_t* value) {
    return capi_brep_guarded([&]() -> mio_status {
        capi_brep_check(info, section, index);
        const meshioplusplus::BrepInfo& b = capi_brep(info);
        const CapiBrepIndex& idx = *info->mBrepIndex;
        const auto i = static_cast<std::size_t>(index);
        std::int64_t v = 0;
        bool found = true;
        if (field >= 10) {
            const meshioplusplus::BrepCurve* c = capi_brep_curve(info, section, i);
            if (!c)
                capi_brep_bad_field(section, field);
            if (field == 10)
                v = c->mDegree;
            else if (field == 11)
                v = c->mRational;
            else if (field == 12)
                v = c->mReducedKnots;
            else
                found = false;
        } else if (section == MIO_BREP_BREPS && field == 0) {
            v = b.mBreps[i].mId;
        } else if (section == MIO_BREP_PATCHES && (field == 0 || field == 1)) {
            const auto& d = b.mPatches[i].mDegrees;
            v = static_cast<std::size_t>(field) < d.size() ? d[static_cast<std::size_t>(field)] : 0;
        } else if (section == MIO_BREP_FACES) {
            const auto& f = b.mFaces[i];
            const std::int64_t fields[] = {f.mId,
                                           f.mPatch,
                                           f.mSwappedNormal,
                                           f.mTrimmed,
                                           f.mRational,
                                           f.mReducedKnots,
                                           idx.mFaces[i][0],
                                           idx.mFaces[i][1],
                                           idx.mFaces[i][2],
                                           idx.mFaces[i][3]};
            if (field >= 0 && field < 10)
                v = fields[field];
            else
                found = false;
        } else if (section == MIO_BREP_LOOPS && field >= 0 && field <= 2) {
            const auto& l = idx.mLoops[i];
            v = field == 0   ? l.mFace
                : field == 1 ? l.mFirstTrim
                             : static_cast<std::int64_t>(l.mpLoop->mTrims.size());
        } else if (section == MIO_BREP_TRIMS && field >= 0 && field <= 3) {
            const auto& t = idx.mTrims[i];
            v = field == 0   ? t.mFace
                : field == 1 ? t.mLoop
                : field == 2 ? t.mpTrim->mTrimIndex
                             : static_cast<std::int64_t>(t.mpTrim->mDirection);
        } else if (section == MIO_BREP_EDGES && field >= 0 && field <= 3) {
            const auto& e = b.mEdges[i];
            v = field == 0   ? e.mId
                : field == 1 ? static_cast<std::int64_t>(e.mHasCurve)
                             : idx.mEdgeRefs[i][static_cast<std::size_t>(field - 2)];
        } else if (section == MIO_BREP_VERTICES && field >= 0 && field <= 2) {
            v = field == 0 ? b.mVertices[i].mId
                           : idx.mVertexRefs[i][static_cast<std::size_t>(field - 1)];
        } else if (section == MIO_BREP_REFS && field >= 0 && field <= 2) {
            const auto* r = idx.mRefs[i];
            v = field == 0 ? r->mId : field == 1 ? r->mTrimIndex : r->mRelativeDirection;
        } else {
            found = false;
        }
        if (!found)
            capi_brep_bad_field(section, field);
        if (value)
            *value = v;
        return MIO_OK;
    });
}

mio_status mio_brep_info_array(const mio_format_info* info, int32_t section, int64_t index,
                               int32_t field, const void** data, mio_dtype* dtype, int32_t* ndim,
                               int64_t* shape) {
    return capi_brep_guarded([&]() -> mio_status {
        capi_brep_check(info, section, index);
        const meshioplusplus::BrepInfo& b = capi_brep(info);
        const auto i = static_cast<std::size_t>(index);
        if (field >= 10) {
            const meshioplusplus::BrepCurve* c = capi_brep_curve(info, section, i);
            if (c) {
                if (field == 10)
                    return array_out(c->mControlPoints, data, dtype, ndim, shape);
                if (field == 11)
                    return array_out(c->mWeights, data, dtype, ndim, shape);
                if (field == 12)
                    return array_out(capi_doubles_view(c->mKnots), data, dtype, ndim, shape);
                if (field == 13)
                    return array_out(c->mIds, data, dtype, ndim, shape);
                if (field == 14)
                    return array_out(capi_doubles_view(c->mActiveRange), data, dtype, ndim, shape);
            }
        } else if (section == MIO_BREP_BREPS && field >= 0 && field <= 2) {
            const auto& br = b.mBreps[i];
            const std::vector<std::int64_t>* lists[] = {&br.mFaces, &br.mEdges, &br.mVertices};
            return array_out(capi_int64s_view(*lists[field]), data, dtype, ndim, shape);
        } else if (section == MIO_BREP_PATCHES) {
            const auto& p = b.mPatches[i];
            if (field == 0)
                return array_out(p.mControlPoints, data, dtype, ndim, shape);
            if (field == 1)
                return array_out(p.mWeights, data, dtype, ndim, shape);
            if ((field == 2 || field == 3) && static_cast<std::size_t>(field - 2) < p.mKnots.size())
                return array_out(capi_doubles_view(p.mKnots[static_cast<std::size_t>(field - 2)]),
                                 data, dtype, ndim, shape);
        } else if (section == MIO_BREP_FACES && field == 0) {
            return array_out(b.mFaces[i].mControlPointIds, data, dtype, ndim, shape);
        } else if (section == MIO_BREP_REFS && field == 0) {
            return array_out(capi_doubles_view(info->mBrepIndex->mRefs[i]->mLocalCoordinates), data,
                             dtype, ndim, shape);
        }
        capi_brep_bad_field(section, field);
    });
}

int64_t mio_brep_info_string(const mio_format_info* info, int32_t section, int64_t index,
                             int32_t field, char* buf, int64_t buflen) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        capi_brep_check(info, section, index);
        const meshioplusplus::BrepInfo& b = capi_brep(info);
        const CapiBrepIndex& idx = *info->mBrepIndex;
        const auto i = static_cast<std::size_t>(index);
        const std::string* text = nullptr;
        if (field == 10) {
            if (const meshioplusplus::BrepCurve* c = capi_brep_curve(info, section, i))
                text = &c->mExtra;
        } else if (section == MIO_BREP_FILE && field == 0) {
            text = &b.mExtra;
        } else if (section == MIO_BREP_BREPS && (field == 0 || field == 1)) {
            text = field == 0 ? &b.mBreps[i].mName : &b.mBreps[i].mExtra;
        } else if (section == MIO_BREP_FACES && field >= 0 && field <= 2) {
            const auto& f = b.mFaces[i];
            text = field == 0 ? &f.mName : field == 1 ? &f.mExtra : &f.mSurfaceExtra;
        } else if (section == MIO_BREP_LOOPS && (field == 0 || field == 1)) {
            text = field == 0 ? &idx.mLoops[i].mpLoop->mType : &idx.mLoops[i].mpLoop->mExtra;
        } else if (section == MIO_BREP_TRIMS && field == 0) {
            text = &idx.mTrims[i].mpTrim->mExtra;
        } else if (section == MIO_BREP_EDGES && (field == 0 || field == 1)) {
            text = field == 0 ? &b.mEdges[i].mName : &b.mEdges[i].mExtra;
        } else if (section == MIO_BREP_VERTICES && (field == 0 || field == 1)) {
            text = field == 0 ? &b.mVertices[i].mName : &b.mVertices[i].mExtra;
        } else if (section == MIO_BREP_REFS && (field == 0 || field == 1)) {
            text = field == 0 ? &idx.mRefs[i]->mName : &idx.mRefs[i]->mExtra;
        }
        if (!text)
            capi_brep_bad_field(section, field);
        return copy_string(*text, buf, buflen);
    });
}

namespace {

/// The MdpaInfo behind a handle, or a throw naming what it is instead.
const meshioplusplus::MdpaInfo& capi_mdpa(const mio_format_info* pInfo) {
    if (!pInfo)
        throw std::invalid_argument("meshio++: format info handle is NULL");
    const auto* p_mdpa = std::get_if<meshioplusplus::MdpaInfo>(&pInfo->mInfo);
    if (!p_mdpa)
        throw std::invalid_argument("meshio++: not an mdpa side channel ('" + pInfo->mFormat +
                                    "')");
    return *p_mdpa;
}

std::size_t capi_mdpa_count(const meshioplusplus::MdpaInfo& rInfo, int32_t Section) {
    switch (Section) {
        case MIO_MDPA_PROPERTIES:
            return rInfo.mProperties.size();
        case MIO_MDPA_ENTITY_NAMES:
            return rInfo.mEntityNames.size();
        case MIO_MDPA_SKIPPED:
            return rInfo.mSkippedConstructs.size();
        case MIO_MDPA_MODEL_PART_DATA:
            return 1;
        case MIO_MDPA_TABLES:
            return rInfo.mTables.size();
        case MIO_MDPA_GEOMETRIES:
            return rInfo.mGeometries.size();
        case MIO_MDPA_MESH_BLOCKS:
            return rInfo.mMeshBlocks.size();
        case MIO_MDPA_SUBMODELPARTS:
            return rInfo.mSubModelParts.size();
        case MIO_MDPA_RAW_BLOCKS:
            return rInfo.mRawBlocks.size();
        default:
            throw std::invalid_argument("meshio++: unknown mdpa section " +
                                        std::to_string(Section));
    }
}

std::size_t capi_mdpa_index(const meshioplusplus::MdpaInfo& rInfo, int32_t Section, int64_t Index) {
    const std::size_t n = capi_mdpa_count(rInfo, Section);
    if (Index < 0 || static_cast<std::size_t>(Index) >= n)
        throw std::out_of_range("meshio++: mdpa section " + std::to_string(Section) + " index " +
                                std::to_string(Index) + " out of range (" + std::to_string(n) +
                                " items)");
    return static_cast<std::size_t>(Index);
}

[[noreturn]] void capi_mdpa_bad_field(int32_t Section, int32_t Field) {
    throw std::invalid_argument("meshio++: mdpa section " + std::to_string(Section) +
                                " has no such field " + std::to_string(Field));
}

/// The key/value entries of item @p Index of a section that has them.
const std::vector<meshioplusplus::PropertyValue>& capi_mdpa_data(
    const meshioplusplus::MdpaInfo& rInfo, int32_t Section, int64_t Index) {
    const std::size_t i = capi_mdpa_index(rInfo, Section, Index);
    switch (Section) {
        case MIO_MDPA_PROPERTIES:
            return rInfo.mProperties[i].mValues;
        case MIO_MDPA_MODEL_PART_DATA:
            return rInfo.mModelPartData;
        case MIO_MDPA_MESH_BLOCKS:
            return rInfo.mMeshBlocks[i].mData;
        case MIO_MDPA_SUBMODELPARTS:
            return rInfo.mSubModelParts[i].mData;
        default:
            throw std::invalid_argument("meshio++: mdpa section " + std::to_string(Section) +
                                        " has no data entries");
    }
}

const meshioplusplus::PropertyValue& capi_mdpa_entry(const meshioplusplus::MdpaInfo& rInfo,
                                                     int32_t Section, int64_t Index,
                                                     int64_t Entry) {
    const auto& data = capi_mdpa_data(rInfo, Section, Index);
    if (Entry < 0 || static_cast<std::size_t>(Entry) >= data.size())
        throw std::out_of_range("meshio++: mdpa data entry " + std::to_string(Entry) +
                                " out of range (" + std::to_string(data.size()) + " entries)");
    return data[static_cast<std::size_t>(Entry)];
}

/// Borrow a std::vector<int64_t> as a 1-D Int64 array.
mio_status capi_ids_out(const std::vector<std::int64_t>& rIds, const void** ppData,
                        mio_dtype* pDtype, int32_t* pNdim, int64_t* pShape) {
    if (ppData)
        *ppData = rIds.data();
    if (pDtype)
        *pDtype = MIO_INT64;
    if (pNdim)
        *pNdim = 1;
    if (pShape)
        pShape[0] = static_cast<int64_t>(rIds.size());
    return MIO_OK;
}

}  // namespace

int64_t mio_mdpa_info_count(const mio_format_info* info, int32_t section) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        return static_cast<std::int64_t>(capi_mdpa_count(capi_mdpa(info), section));
    });
}

int64_t mio_mdpa_info_string(const mio_format_info* info, int32_t section, int64_t index,
                             int32_t field, char* buf, int64_t buflen) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        const meshioplusplus::MdpaInfo& r_info = capi_mdpa(info);
        const std::size_t i = capi_mdpa_index(r_info, section, index);
        switch (section) {
            case MIO_MDPA_ENTITY_NAMES:
                if (field == 0)
                    return copy_string(r_info.mEntityNames[i].mName, buf, buflen);
                break;
            case MIO_MDPA_SKIPPED:
                if (field == 0)
                    return copy_string(r_info.mSkippedConstructs[i], buf, buflen);
                break;
            case MIO_MDPA_TABLES:
                if (field == 0)
                    return copy_string(r_info.mTables[i].mKey, buf, buflen);
                break;
            case MIO_MDPA_GEOMETRIES:
                if (field == 0)
                    return copy_string(r_info.mGeometries[i].mName, buf, buflen);
                if (field == 1)
                    return copy_string(r_info.mGeometries[i].mType, buf, buflen);
                break;
            case MIO_MDPA_SUBMODELPARTS:
                if (field == 0)
                    return copy_string(r_info.mSubModelParts[i].mName, buf, buflen);
                break;
            case MIO_MDPA_RAW_BLOCKS:
                if (field == 0)
                    return copy_string(r_info.mRawBlocks[i].mHeader, buf, buflen);
                if (field == 1)
                    return copy_string(r_info.mRawBlocks[i].mBody, buf, buflen);
                if (field == 2)
                    return copy_string(r_info.mRawBlocks[i].mEnd, buf, buflen);
                break;
            default:
                break;
        }
        capi_mdpa_bad_field(section, field);
    });
}

mio_status mio_mdpa_info_int(const mio_format_info* info, int32_t section, int64_t index,
                             int32_t field, int64_t* value) {
    return guarded([&]() -> mio_status {
        if (!value)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: value is NULL");
        const meshioplusplus::MdpaInfo& r_info = capi_mdpa(info);
        const std::size_t i = capi_mdpa_index(r_info, section, index);
        if (field == 0) {
            switch (section) {
                case MIO_MDPA_PROPERTIES:
                    *value = r_info.mProperties[i].mId;
                    return MIO_OK;
                case MIO_MDPA_ENTITY_NAMES:
                    *value = r_info.mEntityNames[i].mIsCondition ? 1 : 0;
                    return MIO_OK;
                case MIO_MDPA_MESH_BLOCKS:
                    *value = r_info.mMeshBlocks[i].mId;
                    return MIO_OK;
                default:
                    break;
            }
        }
        capi_mdpa_bad_field(section, field);
    });
}

mio_status mio_mdpa_info_array(const mio_format_info* info, int32_t section, int64_t index,
                               int32_t field, const void** data, mio_dtype* dtype, int32_t* ndim,
                               int64_t* shape) {
    return guarded([&]() -> mio_status {
        const meshioplusplus::MdpaInfo& r_info = capi_mdpa(info);
        const std::size_t i = capi_mdpa_index(r_info, section, index);
        switch (section) {
            case MIO_MDPA_TABLES:
                if (field == 0)
                    return array_out(r_info.mTables[i].mValues, data, dtype, ndim, shape);
                break;
            case MIO_MDPA_GEOMETRIES:
                if (field == 0)
                    return array_out(r_info.mGeometries[i].mConn, data, dtype, ndim, shape);
                if (field == 1)
                    return capi_ids_out(r_info.mGeometries[i].mIds, data, dtype, ndim, shape);
                break;
            case MIO_MDPA_MESH_BLOCKS:
                if (field == 0)
                    return capi_ids_out(r_info.mMeshBlocks[i].mNodes, data, dtype, ndim, shape);
                if (field == 1)
                    return capi_ids_out(r_info.mMeshBlocks[i].mElementIds, data, dtype, ndim,
                                        shape);
                if (field == 2)
                    return capi_ids_out(r_info.mMeshBlocks[i].mConditionIds, data, dtype, ndim,
                                        shape);
                break;
            case MIO_MDPA_SUBMODELPARTS:
                if (field == 0)
                    return capi_ids_out(r_info.mSubModelParts[i].mTables, data, dtype, ndim, shape);
                if (field == 1)
                    return capi_ids_out(r_info.mSubModelParts[i].mGeometryIds, data, dtype, ndim,
                                        shape);
                if (field == 2)
                    return capi_ids_out(r_info.mSubModelParts[i].mConstraintIds, data, dtype, ndim,
                                        shape);
                break;
            default:
                break;
        }
        capi_mdpa_bad_field(section, field);
    });
}

int64_t mio_mdpa_info_data_count(const mio_format_info* info, int32_t section, int64_t index) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        return static_cast<std::int64_t>(capi_mdpa_data(capi_mdpa(info), section, index).size());
    });
}

int32_t mio_mdpa_info_data_kind(const mio_format_info* info, int32_t section, int64_t index,
                                int64_t entry) {
    return guarded_ptr(std::int32_t(-1), [&]() -> std::int32_t {
        const auto& r_v = capi_mdpa_entry(capi_mdpa(info), section, index, entry);
        if (r_v.mIsTable)
            return MIO_MDPA_VALUE_TABLE;
        return r_v.IsText() ? MIO_MDPA_VALUE_TEXT : MIO_MDPA_VALUE_NUMBER;
    });
}

int64_t mio_mdpa_info_data_string(const mio_format_info* info, int32_t section, int64_t index,
                                  int64_t entry, int32_t field, char* buf, int64_t buflen) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        const auto& r_v = capi_mdpa_entry(capi_mdpa(info), section, index, entry);
        if (field == 0)
            return copy_string(r_v.mKey, buf, buflen);
        if (field == 1)
            return copy_string(r_v.mText, buf, buflen);
        capi_mdpa_bad_field(section, field);
    });
}

mio_status mio_mdpa_info_data_array(const mio_format_info* info, int32_t section, int64_t index,
                                    int64_t entry, const void** data, mio_dtype* dtype,
                                    int32_t* ndim, int64_t* shape) {
    return guarded([&]() -> mio_status {
        const auto& r_v = capi_mdpa_entry(capi_mdpa(info), section, index, entry);
        if (r_v.mValues.Size() == 0) {
            // Text: an empty Float64 array, so a caller can size-check first.
            static const std::vector<std::int64_t> none;
            capi_ids_out(none, data, dtype, ndim, shape);
            if (dtype)
                *dtype = MIO_FLOAT64;
            return MIO_OK;
        }
        return array_out(r_v.mValues, data, dtype, ndim, shape);
    });
}

namespace {

/// The thread-local stack of open provenance scopes -- the exact twin of the
/// pybind bridge's `provenance_scope_stack()` (bindings/python/_core.cpp),
/// needed for the identical reason: a C caller's begin/end calls are two
/// separate ABI crossings, so something has to keep the scope alive between
/// them rather than relying on C++ RAII across the boundary.
std::vector<std::unique_ptr<meshioplusplus::detail::ProvenanceScope>>& provenance_scope_stack() {
    thread_local std::vector<std::unique_ptr<meshioplusplus::detail::ProvenanceScope>> stack;
    return stack;
}

}  // namespace

mio_status mio_provenance_scope_begin(int mode) {
    return guarded([&]() -> mio_status {
        if (mode < 0 || mode > 2)
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: unknown mio_provenance_mode " + std::to_string(mode));
        provenance_scope_stack().push_back(
            std::make_unique<meshioplusplus::detail::ProvenanceScope>(
                static_cast<meshioplusplus::detail::ProvenanceMode>(mode)));
        return MIO_OK;
    });
}

void mio_provenance_scope_end(void) {
    auto& stack = provenance_scope_stack();
    if (!stack.empty())
        stack.pop_back();
}

void mio_provenance_note(const char* category, const char* detail) {
    if (!category || !detail)
        return;
    meshioplusplus::detail::provenance_note(category, detail);
}

void mio_provenance_set_source(const char* path, const char* format) {
    meshioplusplus::detail::provenance_set_source(path ? path : "", format ? format : "");
}

void mio_provenance_set_target(const char* format, const char* encoding, const char* codec,
                               const char* float_format) {
    meshioplusplus::detail::provenance_set_target(format ? format : "", encoding ? encoding : "",
                                                  codec ? codec : "",
                                                  float_format ? float_format : "");
}

mio_status mio_convert(const char* in_path, const char* in_format, const char* out_path,
                       const char* out_format) {
    return guarded([&]() -> mio_status {
        if (!in_path || !out_path)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: path is NULL");
        std::string rfmt = meshioplusplus::resolve_format(in_path, format_or_empty(in_format));
        std::string wfmt =
            meshioplusplus::resolve_write_format(out_path, format_or_empty(out_format));
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

mio_mesh* mio_merge(const mio_mesh* const* meshes, int64_t count, int weld, double atol,
                    int source_tag, int data_policy, int drop_duplicate_cells) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!meshes || count < 1)
            throw meshioplusplus::ReadError("meshio++: need at least one input mesh");
        std::vector<const meshioplusplus::Mesh*> ptrs;
        ptrs.reserve(static_cast<std::size_t>(count));
        for (int64_t i = 0; i < count; ++i) {
            if (!meshes[i])
                throw meshioplusplus::ReadError("meshio++: input mesh is NULL");
            ptrs.push_back(&meshes[i]->mMesh);
        }
        meshioplusplus::MergeOptions opts;
        opts.weld = weld != 0;
        opts.atol = atol;
        opts.source_tag = source_tag != 0;
        opts.drop_duplicate_cells = drop_duplicate_cells != 0;
        opts.data_policy = (data_policy == 1) ? meshioplusplus::MergeDataPolicy::Fill
                                              : meshioplusplus::MergeDataPolicy::Intersection;
        return new mio_mesh{meshioplusplus::merge(ptrs, opts).mMesh};
    });
}

mio_mesh* mio_transform(const mio_mesh* mesh, const double* matrix, int rotate_vector_data) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        if (!matrix)
            throw meshioplusplus::ReadError("meshio++: matrix is NULL");
        meshioplusplus::AffineTransform xf = meshioplusplus::transform_from_matrix(matrix);
        return new mio_mesh{meshioplusplus::transform(mesh->mMesh, xf, rotate_vector_data != 0)};
    });
}

mio_mesh* mio_clean(const mio_mesh* mesh, int weld, double atol, int remove_orphans,
                    int drop_degenerate, int drop_duplicate_cells, int64_t* points_welded,
                    int64_t* points_removed_orphan, int64_t* cells_dropped_degenerate,
                    int64_t* cells_dropped_duplicate) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::CleanOptions opts;
        opts.weld = weld != 0;
        opts.atol = atol;
        opts.remove_orphans = remove_orphans != 0;
        opts.drop_degenerate = drop_degenerate != 0;
        opts.drop_duplicate_cells = drop_duplicate_cells != 0;
        meshioplusplus::CleanResult r = meshioplusplus::clean(mesh->mMesh, opts);
        if (points_welded)
            *points_welded = r.mPointsWelded;
        if (points_removed_orphan)
            *points_removed_orphan = r.mPointsRemovedOrphan;
        if (cells_dropped_degenerate)
            *cells_dropped_degenerate = r.mCellsDroppedDegenerate;
        if (cells_dropped_duplicate)
            *cells_dropped_duplicate = r.mCellsDroppedDuplicate;
        return new mio_mesh{std::move(r.mMesh)};
    });
}

namespace {

/// 0-based id list -> a per-point pin mask, range-checked against `pMesh`'s
/// own point count. `pFrozen == nullptr` (or `numFrozen <= 0`) means "no
/// extra pins", the empty-mask case every `mFrozen`-accepting options struct
/// treats as unset. Mirrors bindings/python/_core.cpp's identical id-list ->
/// mask conversion.
std::vector<std::uint8_t> capi_frozen_mask(const mio_mesh* pMesh, const int64_t* pFrozen,
                                           int64_t numFrozen, const char* pOp) {
    std::vector<std::uint8_t> mask;
    if (pFrozen == nullptr || numFrozen <= 0)
        return mask;
    const std::size_t num_points = pMesh->mMesh.NumPoints();
    mask.assign(num_points, 0);
    for (int64_t i = 0; i < numFrozen; ++i) {
        const int64_t id = pFrozen[i];
        if (id < 0 || static_cast<std::size_t>(id) >= num_points)
            throw meshioplusplus::ReadError(std::string("meshio++: ") + pOp + ": frozen node id " +
                                            std::to_string(id) + " is out of range");
        mask[static_cast<std::size_t>(id)] = 1;
    }
    return mask;
}

/// The shared MMG parameters of the C API (roadmap 3.2).
meshioplusplus::MmgParameters capi_mmg_params(const mio_mesh* pMesh, const mio_mmg_params& rIn,
                                              const char* pOp) {
    meshioplusplus::MmgParameters p;
    p.mHmin = rIn.hmin;
    p.mHmax = rIn.hmax;
    p.mHausd = rIn.hausd;
    p.mHgrad = rIn.hgrad;
    p.mHgradreq = rIn.hgradreq;
    p.mAngleDetection = rIn.angle_detection;
    p.mOptim = rIn.optim != 0;
    p.mNoInsert = rIn.noinsert != 0;
    p.mNoSwap = rIn.noswap != 0;
    p.mNoMove = rIn.nomove != 0;
    p.mNoSurf = rIn.nosurf != 0;
    p.mNreg = rIn.nreg != 0;
    p.mOpnbdy = rIn.opnbdy != 0;
    p.mVerbose = rIn.verbose;
    p.mMem = rIn.mem;
    p.mFrozen = capi_frozen_mask(pMesh, rIn.frozen, rIn.num_frozen, pOp);
    for (int64_t i = 0; rIn.required && i < rIn.num_required; ++i)
        if (rIn.required[i])
            p.mRequired.emplace_back(rIn.required[i]);
    for (int64_t i = 0; rIn.local && i < rIn.num_local; ++i) {
        if (!rIn.local[i].region)
            throw std::invalid_argument(std::string("meshio++: ") + pOp +
                                        ": a local parameter names no region");
        p.mLocal.push_back(
            {rIn.local[i].region, rIn.local[i].hmin, rIn.local[i].hmax, rIn.local[i].hausd});
    }
    return p;
}

void capi_mmg_report(const meshioplusplus::MmgReport& rR, mio_mmg_report* pOut) {
    if (!pOut)
        return;
    *pOut = mio_mmg_report{};
    pOut->status = static_cast<int32_t>(rR.mStatus);
    std::strncpy(pOut->library, rR.mLibrary.c_str(), sizeof(pOut->library) - 1);
    pOut->num_points_in = rR.mNumPointsIn;
    pOut->num_points_out = rR.mNumPointsOut;
    pOut->num_elements_in = rR.mNumElementsIn;
    pOut->num_elements_out = rR.mNumElementsOut;
    pOut->quality_min = rR.mQuality.mMin;
    pOut->quality_mean = rR.mQuality.mMean;
    pOut->quality_max = rR.mQuality.mMax;
    for (std::size_t k = 0; k < 10; ++k)
        pOut->quality_histogram[k] = rR.mQuality.mHistogram[k];
    pOut->num_warnings = static_cast<int64_t>(rR.mWarnings.size());
}

static_assert(sizeof(mio_mmg_local) == 32, "C ABI: mio_mmg_local layout changed");
static_assert(sizeof(mio_mmg_params) == 168, "C ABI: mio_mmg_params layout changed");
static_assert(sizeof(mio_mmg_report) == 192, "C ABI: mio_mmg_report layout changed");
static_assert(sizeof(mio_adapt_opts) == 224, "C ABI: mio_adapt_opts layout changed");
static_assert(sizeof(mio_mmg_material) == 32, "C ABI: mio_mmg_material layout changed");
static_assert(sizeof(mio_discretize_levelset_opts) == 240,
              "C ABI: mio_discretize_levelset_opts layout changed");
static_assert(sizeof(mio_mesh_domain_opts) == 232, "C ABI: mio_mesh_domain_opts layout changed");
static_assert(sizeof(mio_move_opts) == 232, "C ABI: mio_move_opts layout changed");

mio_mesh* capi_smooth(const mio_mesh* pMesh, const meshioplusplus::SmoothOptions& rOptions,
                      int64_t* pNodesMoved, double* pMaxDisplacement, int64_t* pSkippedInversion) {
    if (!pMesh)
        throw meshioplusplus::ReadError("meshio++: mesh is NULL");
    meshioplusplus::SmoothResult r = meshioplusplus::smooth(pMesh->mMesh, rOptions);
    if (pNodesMoved)
        *pNodesMoved = r.mNumNodesMoved;
    if (pMaxDisplacement)
        *pMaxDisplacement = r.mMaxDisplacement;
    if (pSkippedInversion)
        *pSkippedInversion = r.mNumSkippedInversion;
    return new mio_mesh{std::move(r.mMesh)};
}

}  // namespace

mio_mesh* mio_smooth(const mio_mesh* mesh, const char* method, int iterations, double lambda,
                     double mu, int fix_boundary, int preserve_features, double feature_angle,
                     int guard_inversion, int64_t* nodes_moved, double* max_displacement,
                     int64_t* skipped_inversion) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::SmoothOptions opts;
        // A negative lambda is the "this method's own default" sentinel and is
        // passed through unchanged.
        opts.mMethod = meshioplusplus::smooth_method_from_name(method ? method : "taubin");
        opts.mIterations = iterations;
        opts.mLambda = lambda;
        opts.mMu = mu;
        opts.mFixBoundary = fix_boundary != 0;
        opts.mPreserveFeatures = preserve_features != 0;
        opts.mFeatureAngleDeg = feature_angle;
        opts.mGuardInversion = guard_inversion != 0;
        return capi_smooth(mesh, opts, nodes_moved, max_displacement, skipped_inversion);
    });
}

void mio_smooth_opts_init(mio_smooth_opts* opts) {
    if (!opts)
        return;
    *opts = mio_smooth_opts{};
    opts->method = nullptr;  // "taubin"
    opts->iterations = 10;
    opts->lambda = -1.0;
    opts->mu = -0.34;
    opts->fix_boundary = 1;
    opts->preserve_features = 1;
    opts->feature_angle = 30.0;
    opts->guard_inversion = 1;
}

mio_mesh* mio_smooth_ex(const mio_mesh* mesh, const mio_smooth_opts* opts, int64_t* nodes_moved,
                        double* max_displacement, int64_t* skipped_inversion) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        mio_smooth_opts defaults;
        mio_smooth_opts_init(&defaults);
        const mio_smooth_opts& o = opts ? *opts : defaults;
        meshioplusplus::SmoothOptions options;
        options.mMethod = meshioplusplus::smooth_method_from_name(o.method ? o.method : "taubin");
        options.mIterations = o.iterations;
        options.mLambda = o.lambda;
        options.mMu = o.mu;
        options.mFixBoundary = o.fix_boundary != 0;
        options.mPreserveFeatures = o.preserve_features != 0;
        options.mFeatureAngleDeg = o.feature_angle;
        options.mGuardInversion = o.guard_inversion != 0;
        options.mFrozen = capi_frozen_mask(mesh, o.frozen, o.num_frozen, "smooth");
        return capi_smooth(mesh, options, nodes_moved, max_displacement, skipped_inversion);
    });
}

static_assert(sizeof(mio_smooth_opts) == 104, "mio_smooth_opts grew outside its reserved tail");

mio_mesh* mio_optimize_volume(const mio_mesh* mesh, int max_iterations, int relocate, int flip,
                              int preserve_boundary, double min_improvement, int64_t* num_flips,
                              int64_t* num_23_flips, int64_t* num_32_flips,
                              int64_t* num_vertices_moved, int64_t* num_tets,
                              double* min_quality_before, double* min_quality_after) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::OptimizeVolumeOptions opts;
        opts.mMaxIterations = max_iterations;
        opts.mRelocate = relocate != 0;
        opts.mFlip = flip != 0;
        opts.mPreserveBoundary = preserve_boundary != 0;
        opts.mMinImprovement = min_improvement;
        meshioplusplus::OptimizeVolumeResult r = meshioplusplus::optimize_volume(mesh->mMesh, opts);
        if (num_flips)
            *num_flips = r.mNumFlips;
        if (num_23_flips)
            *num_23_flips = r.mNum23Flips;
        if (num_32_flips)
            *num_32_flips = r.mNum32Flips;
        if (num_vertices_moved)
            *num_vertices_moved = r.mNumVerticesMoved;
        if (num_tets)
            *num_tets = r.mNumTets;
        if (min_quality_before)
            *min_quality_before = r.mMinQualityBefore;
        if (min_quality_after)
            *min_quality_after = r.mMinQualityAfter;
        return new mio_mesh{std::move(r.mMesh)};
    });
}

mio_mesh* mio_undo_green(const mio_mesh* coarse, const mio_mesh* fine, int64_t* num_groups_undone,
                         int64_t* num_cells_removed) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!coarse)
            throw meshioplusplus::ReadError("meshio++: coarse mesh is NULL");
        if (!fine)
            throw meshioplusplus::ReadError("meshio++: fine mesh is NULL");
        meshioplusplus::UndoGreenResult r = meshioplusplus::undo_green(coarse->mMesh, fine->mMesh);
        if (num_groups_undone)
            *num_groups_undone = r.mNumGroupsUndone;
        if (num_cells_removed)
            *num_cells_removed = r.mNumCellsRemoved;
        return new mio_mesh{std::move(r.mMesh)};
    });
}

mio_mesh* mio_interpolate(const mio_mesh* source, const mio_mesh* target, const char* method,
                          const char* const* arrays, int64_t arrays_count, int extrapolate,
                          double default_value, const char* on_conflict) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!source)
            throw meshioplusplus::ReadError("meshio++: source mesh is NULL");
        if (!target)
            throw meshioplusplus::ReadError("meshio++: target mesh is NULL");
        meshioplusplus::InterpolateOptions opts;
        opts.mMethod = meshioplusplus::interpolate_method_from_name(method ? method : "nearest");
        // NULL or a non-positive count means "every source point_data array"
        // (deliberately laxer than data_name_list — empty-means-all is the
        // core's own convention for this option).
        if (arrays && arrays_count > 0) {
            opts.mArrays.reserve(static_cast<std::size_t>(arrays_count));
            for (int64_t i = 0; i < arrays_count; ++i) {
                if (!arrays[i])
                    throw meshioplusplus::ReadError("meshio++: arrays[" + std::to_string(i) +
                                                    "] is NULL");
                opts.mArrays.emplace_back(arrays[i]);
            }
        }
        opts.mExtrapolate = extrapolate != 0;
        opts.mDefaultValue = default_value;
        opts.mOnConflict =
            meshioplusplus::interpolate_conflict_from_name(on_conflict ? on_conflict : "error");
        return new mio_mesh{meshioplusplus::interpolate(source->mMesh, target->mMesh, opts)};
    });
}

mio_mesh* mio_conservative_interpolate(const mio_mesh* source, const mio_mesh* target,
                                       const char* const* arrays, int64_t arrays_count,
                                       double default_value, const char* on_conflict) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!source)
            throw meshioplusplus::ReadError("meshio++: source mesh is NULL");
        if (!target)
            throw meshioplusplus::ReadError("meshio++: target mesh is NULL");
        meshioplusplus::ConservativeInterpolateOptions opts;
        // NULL or a non-positive count means "every source point_data and
        // cell_data array" (the interpolate convention, extended: there is
        // one algorithm regardless of location here).
        if (arrays && arrays_count > 0) {
            opts.mArrays.reserve(static_cast<std::size_t>(arrays_count));
            for (int64_t i = 0; i < arrays_count; ++i) {
                if (!arrays[i])
                    throw meshioplusplus::ReadError("meshio++: arrays[" + std::to_string(i) +
                                                    "] is NULL");
                opts.mArrays.emplace_back(arrays[i]);
            }
        }
        opts.mDefaultValue = default_value;
        opts.mOnConflict = meshioplusplus::conservative_interpolate_conflict_from_name(
            on_conflict ? on_conflict : "error");
        return new mio_mesh{
            meshioplusplus::conservative_interpolate(source->mMesh, target->mMesh, opts)};
    });
}

mio_mesh* mio_slice(const mio_mesh* mesh, const double* origin, const double* normal,
                    int record_parent_ids) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !origin || !normal)
            throw meshioplusplus::ReadError("meshio++: mesh/origin/normal is NULL");
        meshioplusplus::SliceOptions opts;
        opts.mOrigin = {origin[0], origin[1], origin[2]};
        opts.mNormal = {normal[0], normal[1], normal[2]};
        opts.mRecordParentIds = record_parent_ids != 0;
        return new mio_mesh{meshioplusplus::slice(mesh->mMesh, opts)};
    });
}

mio_mesh* mio_isosurface(const mio_mesh* mesh, const char* array_name, const double* isovalues,
                         int n_isovalues, int component, int record_parent_ids) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !array_name || !isovalues)
            throw meshioplusplus::ReadError("meshio++: mesh/array_name/isovalues is NULL");
        if (n_isovalues <= 0)
            throw meshioplusplus::ReadError("meshio++: at least one isovalue is required");
        meshioplusplus::IsosurfaceOptions opts;
        opts.mArrayName = array_name;
        opts.mIsovalues.assign(isovalues, isovalues + n_isovalues);
        if (component >= 0)
            opts.mComponent = component;
        opts.mRecordParentIds = record_parent_ids != 0;
        return new mio_mesh{meshioplusplus::isosurface(mesh->mMesh, opts)};
    });
}

mio_mesh* mio_gradient(const mio_mesh* mesh, const char* array_name, const char* op,
                       const char* method, const char* location, const char* output_name,
                       int component, int overwrite, int64_t* num_skipped, int64_t* num_fallback) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !array_name)
            throw meshioplusplus::ReadError("meshio++: mesh/array_name is NULL");
        meshioplusplus::GradientOptions opts;
        opts.mArrayName = array_name;
        opts.mOperator = meshioplusplus::gradient_operator_from_name(op ? op : "");
        opts.mMethod = meshioplusplus::gradient_method_from_name(method ? method : "");
        opts.mLocation =
            meshioplusplus::data_location_from_name((location && *location) ? location : "cell");
        if (output_name)
            opts.mOutputName = output_name;
        // A negative component means EVERY component here -- deliberately the
        // opposite of mio_isosurface, where it means the row magnitude.
        if (component >= 0)
            opts.mComponent = component;
        opts.mOverwrite = overwrite != 0;
        meshioplusplus::GradientResult r = meshioplusplus::gradient(mesh->mMesh, opts);
        if (num_skipped)
            *num_skipped = r.mNumSkipped;
        if (num_fallback)
            *num_fallback = r.mNumFallback;
        return new mio_mesh{std::move(r.mMesh)};
    });
}

mio_mesh* mio_hessian(const mio_mesh* mesh, const char* array_name, const char* method,
                      const char* location, const char* output_name, int overwrite,
                      int64_t* num_skipped, int64_t* num_fallback) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !array_name)
            throw meshioplusplus::ReadError("meshio++: mesh/array_name is NULL");
        meshioplusplus::HessianOptions opts;
        opts.mArrayName = array_name;
        opts.mMethod = meshioplusplus::gradient_method_from_name(method ? method : "");
        opts.mLocation =
            meshioplusplus::data_location_from_name((location && *location) ? location : "cell");
        if (output_name)
            opts.mOutputName = output_name;
        opts.mOverwrite = overwrite != 0;
        meshioplusplus::HessianResult r = meshioplusplus::hessian(mesh->mMesh, opts);
        if (num_skipped)
            *num_skipped = r.mNumSkipped;
        if (num_fallback)
            *num_fallback = r.mNumFallback;
        return new mio_mesh{std::move(r.mMesh)};
    });
}

mio_mesh* mio_estimate_error(const mio_mesh* mesh, const char* array_name, const char* method,
                             const char* marking, double marking_value, const char* output_name,
                             const char* marked_name, int overwrite, double* global_error,
                             int64_t* num_skipped, int64_t* num_marked) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !array_name)
            throw meshioplusplus::ReadError("meshio++: mesh/array_name is NULL");
        meshioplusplus::ErrorOptions opts;
        opts.mArrayName = array_name;
        opts.mMethod = meshioplusplus::error_method_from_name(method ? method : "");
        opts.mMarking = meshioplusplus::error_marking_from_name(marking ? marking : "");
        opts.mMarkingValue = marking_value;
        if (output_name)
            opts.mOutputName = output_name;
        if (marked_name)
            opts.mMarkedName = marked_name;
        opts.mOverwrite = overwrite != 0;
        meshioplusplus::ErrorResult r = meshioplusplus::estimate_error(mesh->mMesh, opts);
        if (global_error)
            *global_error = r.mGlobalError;
        if (num_skipped)
            *num_skipped = r.mNumSkipped;
        if (num_marked)
            *num_marked = r.mNumMarked;
        return new mio_mesh{std::move(r.mMesh)};
    });
}

namespace {

/// Translate the flat option struct into the core's RemeshOptions -- the
/// mio_refine_opts precedent (capi_refine_options).
meshioplusplus::RemeshOptions capi_remesh_options(const mio_remesh_opts& rOpts) {
    meshioplusplus::RemeshOptions options;
    options.mNumClusters = rOpts.num_clusters;
    options.mSubdivide = rOpts.subdivide;
    options.mSubsampleRatio = rOpts.subsample_ratio;
    options.mMaxSubdivide = rOpts.max_subdivide;
    options.mMaxIterations = rOpts.max_iterations;
    options.mMaxRepairPasses = rOpts.max_repair_passes;
    options.mMetric =
        meshioplusplus::remesh_metric_from_name(rOpts.metric ? rOpts.metric : "isotropic");
    options.mGradation = rOpts.gradation;
    options.mPreserveBoundary = rOpts.preserve_boundary != 0;
    options.mMaxAnisotropy = rOpts.max_anisotropy;
    return options;
}

mio_mesh* capi_remesh(const mio_mesh* pMesh, const meshioplusplus::RemeshOptions& rOptions,
                      mio_remesh_report* pReport) {
    if (!pMesh)
        throw meshioplusplus::ReadError("meshio++: mesh is NULL");
    meshioplusplus::RemeshResult r = meshioplusplus::remesh(pMesh->mMesh, rOptions);
    if (pReport) {
        *pReport = mio_remesh_report{};
        pReport->num_clusters = r.mNumClusters;
        pReport->num_iterations = r.mNumIterations;
        pReport->subdivide_applied = r.mSubdivideApplied;
        pReport->num_isolated_clusters = r.mNumIsolatedClusters;
        pReport->num_non_manifold_vertices = r.mNumNonManifoldVertices;
    }
    return new mio_mesh{std::move(r.mMesh)};
}

}  // namespace

static_assert(sizeof(mio_remesh_opts) == 120, "mio_remesh_opts grew outside its reserved tail");
static_assert(sizeof(mio_remesh_report) == 72, "mio_remesh_report grew outside its reserved tail");

void mio_remesh_opts_init(mio_remesh_opts* opts) {
    if (!opts)
        return;
    *opts = mio_remesh_opts{};  // value-initialized: zero every field first
    opts->subdivide = -1;
    opts->max_subdivide = 4;
    opts->subsample_ratio = 10.0;
    opts->max_iterations = 100;
    opts->max_repair_passes = 10;
    opts->preserve_boundary = 1;
    opts->max_anisotropy = meshioplusplus::kRemeshDefaultMaxAnisotropy;
}

mio_mesh* mio_remesh(const mio_mesh* mesh, int64_t num_clusters, int subdivide,
                     double subsample_ratio, int max_subdivide, int max_iterations,
                     int max_repair_passes, const char* metric, double gradation,
                     int preserve_boundary, int64_t* num_clusters_out, int64_t* num_iterations,
                     int* subdivide_applied, int64_t* num_isolated_clusters,
                     int64_t* num_non_manifold_vertices) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        meshioplusplus::RemeshOptions opts;
        opts.mNumClusters = num_clusters;
        opts.mSubdivide = subdivide;
        opts.mSubsampleRatio = subsample_ratio;
        opts.mMaxSubdivide = max_subdivide;
        opts.mMaxIterations = max_iterations;
        opts.mMaxRepairPasses = max_repair_passes;
        opts.mMetric = meshioplusplus::remesh_metric_from_name(metric ? metric : "isotropic");
        opts.mGradation = gradation;
        opts.mPreserveBoundary = preserve_boundary != 0;
        // max_anisotropy has no flat parameter here -- mio_remesh_ex is
        // required to move it away from the default (RemeshOptions's own
        // member-initializer default, inherited unchanged above).
        mio_remesh_report report{};
        mio_mesh* out = capi_remesh(mesh, opts, &report);
        if (num_clusters_out)
            *num_clusters_out = report.num_clusters;
        if (num_iterations)
            *num_iterations = report.num_iterations;
        if (subdivide_applied)
            *subdivide_applied = report.subdivide_applied;
        if (num_isolated_clusters)
            *num_isolated_clusters = report.num_isolated_clusters;
        if (num_non_manifold_vertices)
            *num_non_manifold_vertices = report.num_non_manifold_vertices;
        return out;
    });
}

mio_mesh* mio_remesh_ex(const mio_mesh* mesh, const mio_remesh_opts* opts,
                        mio_remesh_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (opts == nullptr) {
            mio_remesh_opts defaults;
            mio_remesh_opts_init(&defaults);
            return capi_remesh(mesh, capi_remesh_options(defaults), report);
        }
        return capi_remesh(mesh, capi_remesh_options(*opts), report);
    });
}

mio_mesh* mio_crop_bbox(const mio_mesh* mesh, const double* lo, const double* hi, int mode,
                        int record_ids) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !lo || !hi)
            throw meshioplusplus::ReadError("meshio++: mesh/lo/hi is NULL");
        meshioplusplus::CropMode m =
            (mode == 1) ? meshioplusplus::CropMode::Any : meshioplusplus::CropMode::All;
        return new mio_mesh{
            meshioplusplus::crop_bbox(mesh->mMesh, lo, hi, m, record_ids != 0).mMesh};
    });
}

mio_mesh* mio_crop_plane(const mio_mesh* mesh, const double* point, const double* normal, int mode,
                         int record_ids) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !point || !normal)
            throw meshioplusplus::ReadError("meshio++: mesh/point/normal is NULL");
        meshioplusplus::CropMode m =
            (mode == 1) ? meshioplusplus::CropMode::Any : meshioplusplus::CropMode::All;
        return new mio_mesh{
            meshioplusplus::crop_halfspace(mesh->mMesh, point, normal, m, record_ids != 0).mMesh};
    });
}

mio_mesh* mio_crop_predicate(const mio_mesh* mesh, const char* array, int compare, double value,
                             int record_ids) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !array)
            throw meshioplusplus::ReadError("meshio++: mesh/array is NULL");
        if (compare < 0 || compare > MIO_REFINE_NE)
            throw meshioplusplus::ReadError("meshio++: crop: unknown comparison " +
                                            std::to_string(compare));
        return new mio_mesh{
            meshioplusplus::crop_predicate(mesh->mMesh, array,
                                           static_cast<meshioplusplus::RefineCompare>(compare),
                                           value, record_ids != 0)
                .mMesh};
    });
}

mio_split_result* mio_split(const mio_mesh* mesh, const char* by, const char* tag_name) {
    return guarded_ptr(static_cast<mio_split_result*>(nullptr), [&]() -> mio_split_result* {
        if (!mesh || !by)
            throw meshioplusplus::ReadError("meshio++: mesh/by is NULL");
        meshioplusplus::SplitResult r = meshioplusplus::split(
            mesh->mMesh, meshioplusplus::split_by_from_name(by), tag_name ? tag_name : "");
        auto* out = new mio_split_result{};
        out->mMeshes.reserve(r.mPieces.size());
        out->mKeys.reserve(r.mPieces.size());
        for (meshioplusplus::SplitPiece& p : r.mPieces) {
            out->mMeshes.push_back(mio_mesh{std::move(p.mMesh)});
            out->mKeys.push_back(std::move(p.mKey));
        }
        return out;
    });
}

int64_t mio_split_result_count(const mio_split_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return static_cast<int64_t>(result->mMeshes.size());
    });
}

int64_t mio_split_result_key(const mio_split_result* result, int64_t index, char* buf,
                             int64_t buflen) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mKeys.size())
            throw meshioplusplus::ReadError("meshio++: split piece index out of range");
        return copy_string(result->mKeys[static_cast<std::size_t>(index)], buf, buflen);
    });
}

const mio_mesh* mio_split_result_mesh(const mio_split_result* result, int64_t index) {
    return guarded_ptr(static_cast<const mio_mesh*>(nullptr), [&]() -> const mio_mesh* {
        if (!result || index < 0 || static_cast<std::size_t>(index) >= result->mMeshes.size())
            return nullptr;
        return &result->mMeshes[static_cast<std::size_t>(index)];
    });
}

mio_mesh* mio_split_result_take_mesh(mio_split_result* result, int64_t index) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!result || index < 0 || static_cast<std::size_t>(index) >= result->mMeshes.size())
            throw meshioplusplus::ReadError("meshio++: split piece index out of range");
        return new mio_mesh{std::move(result->mMeshes[static_cast<std::size_t>(index)].mMesh)};
    });
}

void mio_split_result_free(mio_split_result* result) {
    delete result;
}

mio_convert_cells_result* mio_convert_cells(const mio_mesh* mesh, const char* mode,
                                            int record_parent_ids) {
    return guarded_ptr(static_cast<mio_convert_cells_result*>(nullptr),
                       [&]() -> mio_convert_cells_result* {
                           if (!mesh || !mode)
                               throw meshioplusplus::ReadError("meshio++: mesh/mode is NULL");
                           meshioplusplus::ConvertCellsOptions options;
                           options.mMode = meshioplusplus::convert_cells_mode_from_name(mode);
                           options.mRecordParentIds = record_parent_ids != 0;
                           meshioplusplus::ConvertCellsResult r =
                               meshioplusplus::convert_cells(mesh->mMesh, options);
                           auto* out = new mio_convert_cells_result{};
                           out->mMesh = mio_mesh{std::move(r.mMesh)};
                           out->mPointMap = std::move(r.mPointMap);
                           out->mCellMaps = std::move(r.mCellMaps);
                           return out;
                       });
}

const mio_mesh* mio_convert_cells_result_mesh(const mio_convert_cells_result* result) {
    return guarded_ptr(static_cast<const mio_mesh*>(nullptr), [&]() -> const mio_mesh* {
        if (!result)
            return nullptr;
        return &result->mMesh;
    });
}

mio_mesh* mio_convert_cells_result_take_mesh(mio_convert_cells_result* result) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return new mio_mesh{std::move(result->mMesh.mMesh)};
    });
}

mio_status mio_convert_cells_result_point_map(const mio_convert_cells_result* result,
                                              const void** data, mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        const NDArray& a = result->mPointMap;
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

int64_t mio_convert_cells_result_num_cell_maps(const mio_convert_cells_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return static_cast<int64_t>(result->mCellMaps.size());
    });
}

mio_status mio_convert_cells_result_cell_map(const mio_convert_cells_result* result, int64_t block,
                                             const void** data, mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (block < 0 || static_cast<std::size_t>(block) >= result->mCellMaps.size())
            return fail(MIO_ERR_NOT_FOUND, "meshio++: cell-map block index out of range");
        const NDArray& a = result->mCellMaps[static_cast<std::size_t>(block)];
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

void mio_convert_cells_result_free(mio_convert_cells_result* result) {
    delete result;
}

mio_subdivide_result* mio_subdivide(const mio_mesh* mesh, int record_parent_ids) {
    return guarded_ptr(static_cast<mio_subdivide_result*>(nullptr), [&]() -> mio_subdivide_result* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::SubdivideOptions options;
        options.mRecordParentIds = record_parent_ids != 0;
        meshioplusplus::SubdivideResult r = meshioplusplus::subdivide(mesh->mMesh, options);
        auto* out = new mio_subdivide_result{};
        out->mMesh = mio_mesh{std::move(r.mMesh)};
        out->mCellMaps = std::move(r.mCellMaps);
        return out;
    });
}

const mio_mesh* mio_subdivide_result_mesh(const mio_subdivide_result* result) {
    return guarded_ptr(static_cast<const mio_mesh*>(nullptr), [&]() -> const mio_mesh* {
        if (!result)
            return nullptr;
        return &result->mMesh;
    });
}

mio_mesh* mio_subdivide_result_take_mesh(mio_subdivide_result* result) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return new mio_mesh{std::move(result->mMesh.mMesh)};
    });
}

int64_t mio_subdivide_result_num_cell_maps(const mio_subdivide_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return static_cast<int64_t>(result->mCellMaps.size());
    });
}

mio_status mio_subdivide_result_cell_map(const mio_subdivide_result* result, int64_t block,
                                         const void** data, mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (block < 0 || static_cast<std::size_t>(block) >= result->mCellMaps.size())
            return fail(MIO_ERR_NOT_FOUND, "meshio++: cell-map block index out of range");
        const NDArray& a = result->mCellMaps[static_cast<std::size_t>(block)];
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

void mio_subdivide_result_free(mio_subdivide_result* result) {
    delete result;
}

mio_agglomerate_result* mio_agglomerate(const mio_mesh* mesh, int64_t target_group_size) {
    return guarded_ptr(
        static_cast<mio_agglomerate_result*>(nullptr), [&]() -> mio_agglomerate_result* {
            if (!mesh)
                throw meshioplusplus::ReadError("meshio++: mesh is NULL");
            if (target_group_size < 0)
                throw std::invalid_argument(
                    "meshio++: agglomerate: target_group_size must be >= 1");
            meshioplusplus::AgglomerateOptions options;
            options.mTargetGroupSize = static_cast<std::size_t>(target_group_size);
            meshioplusplus::AgglomerateResult r = meshioplusplus::agglomerate(mesh->mMesh, options);
            auto* out = new mio_agglomerate_result{};
            out->mMesh = mio_mesh{std::move(r.mMesh)};
            out->mCellMap = std::move(r.mCellMap);
            return out;
        });
}

void mio_agglomerate_opts_init(mio_agglomerate_opts* opts) {
    if (!opts)
        return;
    *opts = mio_agglomerate_opts{};
    opts->target_group_size = 8;
    opts->coplanar_angle = 1.0;
}

mio_agglomerate_result* mio_agglomerate_ex(const mio_mesh* mesh, const mio_agglomerate_opts* opts,
                                           int64_t* num_faces_merged, int64_t* num_rejected) {
    return guarded_ptr(
        static_cast<mio_agglomerate_result*>(nullptr), [&]() -> mio_agglomerate_result* {
            if (!mesh)
                throw meshioplusplus::ReadError("meshio++: mesh is NULL");
            meshioplusplus::AgglomerateOptions options;
            if (opts) {
                if (opts->target_group_size < 1)
                    throw std::invalid_argument(
                        "meshio++: agglomerate: target_group_size must be >= 1");
                options.mTargetGroupSize = static_cast<std::size_t>(opts->target_group_size);
                options.mMergeCoplanarFaces = opts->merge_coplanar_faces != 0;
                options.mCoplanarAngleDeg = opts->coplanar_angle;
                options.mMinSphericity = opts->min_sphericity;
            }
            meshioplusplus::AgglomerateResult r = meshioplusplus::agglomerate(mesh->mMesh, options);
            if (num_faces_merged)
                *num_faces_merged = r.mNumFacesMerged;
            if (num_rejected)
                *num_rejected = r.mNumRejected;
            auto* out = new mio_agglomerate_result{};
            out->mMesh = mio_mesh{std::move(r.mMesh)};
            out->mCellMap = std::move(r.mCellMap);
            return out;
        });
}

const mio_mesh* mio_agglomerate_result_mesh(const mio_agglomerate_result* result) {
    return guarded_ptr(static_cast<const mio_mesh*>(nullptr), [&]() -> const mio_mesh* {
        if (!result)
            return nullptr;
        return &result->mMesh;
    });
}

mio_mesh* mio_agglomerate_result_take_mesh(mio_agglomerate_result* result) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return new mio_mesh{std::move(result->mMesh.mMesh)};
    });
}

mio_status mio_agglomerate_result_cell_map(const mio_agglomerate_result* result, const void** data,
                                           mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        const NDArray& a = result->mCellMap;
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

void mio_agglomerate_result_free(mio_agglomerate_result* result) {
    delete result;
}

// The two flat enums duplicate the C++ ones across the ABI, so drift is a
// compile error rather than a wrong answer -- the mio_region_kind pattern.
static_assert(static_cast<int>(meshioplusplus::RefineClosure::RedGreen) ==
                  MIO_REFINE_CLOSURE_REDGREEN,
              "mio_refine_closure drifted from RefineClosure");
static_assert(static_cast<int>(meshioplusplus::RefineClosure::Propagate) ==
                  MIO_REFINE_CLOSURE_PROPAGATE,
              "mio_refine_closure drifted from RefineClosure");
static_assert(static_cast<int>(meshioplusplus::RefineClosure::Balanced) ==
                  MIO_REFINE_CLOSURE_BALANCED,
              "mio_refine_closure drifted from RefineClosure");
static_assert(static_cast<int>(meshioplusplus::RefineCompare::Less) == MIO_REFINE_LT &&
                  static_cast<int>(meshioplusplus::RefineCompare::LessEqual) == MIO_REFINE_LE &&
                  static_cast<int>(meshioplusplus::RefineCompare::Greater) == MIO_REFINE_GT &&
                  static_cast<int>(meshioplusplus::RefineCompare::GreaterEqual) == MIO_REFINE_GE &&
                  static_cast<int>(meshioplusplus::RefineCompare::Equal) == MIO_REFINE_EQ &&
                  static_cast<int>(meshioplusplus::RefineCompare::NotEqual) == MIO_REFINE_NE,
              "mio_refine_compare drifted from RefineCompare");

// The size is ABI: the Fortran and Julia mirrors of this struct hard-code it,
// and the Julia binding checks it at load. Pinning it here makes a field added
// outside the `reserved` tail a compile error rather than silent corruption.
static_assert(sizeof(mio_refine_opts) == 112, "mio_refine_opts grew outside its reserved tail");

void mio_refine_opts_init(mio_refine_opts* opts) {
    if (!opts)
        return;
    *opts = mio_refine_opts{};  // value-initialized: all zero == refine everything once
    opts->levels = 1;
}

namespace {

/// Translate the flat option struct into the core's RefineOptions.
meshioplusplus::RefineOptions capi_refine_options(const mio_refine_opts& rOpts) {
    meshioplusplus::RefineOptions options;
    options.mLevels = rOpts.levels;
    options.mRecordParentIds = rOpts.record_parent_ids != 0;
    options.mRecordLevels = rOpts.record_levels != 0;
    options.mRecordHierarchy = rOpts.record_hierarchy != 0;
    if (rOpts.cells != nullptr && rOpts.num_cells > 0)
        options.mCells.assign(rOpts.cells, rOpts.cells + rOpts.num_cells);
    if (rOpts.region != nullptr)
        options.mRegion = rOpts.region;
    if (rOpts.predicate_array != nullptr)
        options.mPredicateArray = rOpts.predicate_array;
    options.mPredicateValue = rOpts.predicate_value;
    if (rOpts.closure < 0 || rOpts.closure > MIO_REFINE_CLOSURE_BALANCED)
        throw meshioplusplus::ReadError("meshio++: refine: unknown closure " +
                                        std::to_string(rOpts.closure));
    options.mClosure = static_cast<meshioplusplus::RefineClosure>(rOpts.closure);
    if (rOpts.predicate_op < 0 || rOpts.predicate_op > MIO_REFINE_NE)
        throw meshioplusplus::ReadError("meshio++: refine: unknown comparison " +
                                        std::to_string(rOpts.predicate_op));
    options.mPredicateOp = static_cast<meshioplusplus::RefineCompare>(rOpts.predicate_op);
    return options;
}

mio_refine_result* capi_refine(const mio_mesh* pMesh,
                               const meshioplusplus::RefineOptions& rOptions) {
    if (!pMesh)
        throw meshioplusplus::ReadError("meshio++: mesh is NULL");
    meshioplusplus::RefineResult r = meshioplusplus::refine(pMesh->mMesh, rOptions);
    auto* out = new mio_refine_result{};
    out->mMesh = mio_mesh{std::move(r.mMesh)};
    out->mPointMap = std::move(r.mPointMap);
    out->mCellMaps = std::move(r.mCellMaps);
    return out;
}

}  // namespace

mio_refine_result* mio_refine(const mio_mesh* mesh, int levels, int record_parent_ids) {
    return guarded_ptr(static_cast<mio_refine_result*>(nullptr), [&]() -> mio_refine_result* {
        meshioplusplus::RefineOptions options;
        options.mLevels = levels;
        options.mRecordParentIds = record_parent_ids != 0;
        return capi_refine(mesh, options);
    });
}

mio_refine_result* mio_refine_ex(const mio_mesh* mesh, const mio_refine_opts* opts) {
    return guarded_ptr(static_cast<mio_refine_result*>(nullptr), [&]() -> mio_refine_result* {
        if (opts == nullptr)
            return capi_refine(mesh, meshioplusplus::RefineOptions{});
        return capi_refine(mesh, capi_refine_options(*opts));
    });
}

const mio_mesh* mio_refine_result_mesh(const mio_refine_result* result) {
    return guarded_ptr(static_cast<const mio_mesh*>(nullptr), [&]() -> const mio_mesh* {
        if (!result)
            return nullptr;
        return &result->mMesh;
    });
}

mio_mesh* mio_refine_result_take_mesh(mio_refine_result* result) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return new mio_mesh{std::move(result->mMesh.mMesh)};
    });
}

mio_status mio_refine_result_point_map(const mio_refine_result* result, const void** data,
                                       mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        const NDArray& a = result->mPointMap;
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

int64_t mio_refine_result_num_cell_maps(const mio_refine_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return static_cast<int64_t>(result->mCellMaps.size());
    });
}

mio_status mio_refine_result_cell_map(const mio_refine_result* result, int64_t block,
                                      const void** data, mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (block < 0 || static_cast<std::size_t>(block) >= result->mCellMaps.size())
            return fail(MIO_ERR_NOT_FOUND, "meshio++: cell-map block index out of range");
        const NDArray& a = result->mCellMaps[static_cast<std::size_t>(block)];
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

void mio_refine_result_free(mio_refine_result* result) {
    delete result;
}

static_assert(sizeof(mio_tessellate_opts) == 72,
              "mio_tessellate_opts grew outside its reserved tail");
static_assert(sizeof(mio_tessellate_report) == 112,
              "mio_tessellate_report grew outside its reserved tail");

void mio_tessellate_opts_init(mio_tessellate_opts* opts) {
    if (!opts)
        return;
    *opts = mio_tessellate_opts{};
    opts->levels = 2;
    opts->curved = 1;
    opts->fields = 1;
}

mio_mesh* mio_tessellate(const mio_mesh* mesh, const mio_tessellate_opts* opts,
                         mio_tessellate_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::TessellateOptions options;
        if (opts) {
            options.mLevels = opts->levels;
            options.mCurved = opts->curved != 0;
            options.mFields = opts->fields != 0;
            options.mRecordStencil = opts->record_stencil != 0;
        }
        meshioplusplus::TessellateResult r = meshioplusplus::tessellate(mesh->mMesh, options);
        if (report) {
            *report = mio_tessellate_report{};
            report->num_source_points = r.mNumSourcePoints;
            report->num_source_cells = r.mNumSourceCells;
            report->num_curved_source_cells = r.mNumCurvedSourceCells;
            report->num_pass_through_source_cells = r.mNumPassThroughSourceCells;
            report->num_points = r.mNumPoints;
            report->num_cells = r.mNumCells;
            report->num_boundary_facets = r.mNumBoundaryFacets;
            report->num_interior_facets = r.mNumInteriorFacets;
            report->num_facets_with_bad_count = r.mNumFacetsWithBadCount;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

namespace {

mio_decimate_result* capi_decimate(const mio_mesh* pMesh,
                                   const meshioplusplus::DecimateOptions& rOptions) {
    if (!pMesh)
        throw meshioplusplus::ReadError("meshio++: mesh is NULL");
    meshioplusplus::DecimateResult r = meshioplusplus::decimate(pMesh->mMesh, rOptions);
    auto* out = new mio_decimate_result{};
    out->mMesh = mio_mesh{std::move(r.mMesh)};
    out->mPointMap = std::move(r.mPointMap);
    out->mCellMaps = std::move(r.mCellMaps);
    out->mFacesRemoved = r.mFacesRemoved;
    out->mPointsRemoved = r.mPointsRemoved;
    out->mCollapsesRejected = r.mCollapsesRejected;
    out->mMaxErrorApplied = r.mMaxErrorApplied;
    return out;
}

}  // namespace

mio_decimate_result* mio_decimate(const mio_mesh* mesh, double target_ratio, int64_t target_faces,
                                  double max_error, const char* placement, int preserve_boundary,
                                  int preserve_features, double feature_angle) {
    return guarded_ptr(static_cast<mio_decimate_result*>(nullptr), [&]() -> mio_decimate_result* {
        meshioplusplus::DecimateOptions options;
        options.mTargetRatio = target_ratio;
        options.mTargetFaces = target_faces;
        options.mMaxError = max_error;
        options.mPlacement =
            meshioplusplus::decimate_placement_from_name(placement ? placement : "optimal");
        options.mPreserveBoundary = preserve_boundary != 0;
        options.mPreserveFeatures = preserve_features != 0;
        options.mFeatureAngleDeg = feature_angle;
        return capi_decimate(mesh, options);
    });
}

void mio_decimate_opts_init(mio_decimate_opts* opts) {
    if (!opts)
        return;
    *opts = mio_decimate_opts{};
    opts->target_ratio = -1.0;
    opts->target_faces = -1;
    opts->max_error = -1.0;
    opts->placement = nullptr;  // "optimal"
    opts->preserve_boundary = 1;
    opts->preserve_features = 1;
    opts->feature_angle = 30.0;
}

mio_decimate_result* mio_decimate_ex(const mio_mesh* mesh, const mio_decimate_opts* opts) {
    return guarded_ptr(static_cast<mio_decimate_result*>(nullptr), [&]() -> mio_decimate_result* {
        if (!opts)
            throw meshioplusplus::ReadError("meshio++: opts is NULL");
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::DecimateOptions options;
        options.mTargetRatio = opts->target_ratio;
        options.mTargetFaces = opts->target_faces;
        options.mMaxError = opts->max_error;
        options.mPlacement = meshioplusplus::decimate_placement_from_name(
            opts->placement ? opts->placement : "optimal");
        options.mPreserveBoundary = opts->preserve_boundary != 0;
        options.mPreserveFeatures = opts->preserve_features != 0;
        options.mFeatureAngleDeg = opts->feature_angle;
        options.mFrozen = capi_frozen_mask(mesh, opts->frozen, opts->num_frozen, "decimate");
        return capi_decimate(mesh, options);
    });
}

static_assert(sizeof(mio_decimate_opts) == 96, "mio_decimate_opts grew outside its reserved tail");

const mio_mesh* mio_decimate_result_mesh(const mio_decimate_result* result) {
    return guarded_ptr(static_cast<const mio_mesh*>(nullptr), [&]() -> const mio_mesh* {
        if (!result)
            return nullptr;
        return &result->mMesh;
    });
}

mio_mesh* mio_decimate_result_take_mesh(mio_decimate_result* result) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return new mio_mesh{std::move(result->mMesh.mMesh)};
    });
}

mio_status mio_decimate_result_point_map(const mio_decimate_result* result, const void** data,
                                         mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        const NDArray& a = result->mPointMap;
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

int64_t mio_decimate_result_num_cell_maps(const mio_decimate_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return static_cast<int64_t>(result->mCellMaps.size());
    });
}

mio_status mio_decimate_result_cell_map(const mio_decimate_result* result, int64_t block,
                                        const void** data, mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (block < 0 || static_cast<std::size_t>(block) >= result->mCellMaps.size())
            return fail(MIO_ERR_NOT_FOUND, "meshio++: cell-map block index out of range");
        const NDArray& a = result->mCellMaps[static_cast<std::size_t>(block)];
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

int64_t mio_decimate_result_faces_removed(const mio_decimate_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return result->mFacesRemoved;
    });
}

int64_t mio_decimate_result_points_removed(const mio_decimate_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return result->mPointsRemoved;
    });
}

int64_t mio_decimate_result_collapses_rejected(const mio_decimate_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return result->mCollapsesRejected;
    });
}

double mio_decimate_result_max_error_applied(const mio_decimate_result* result) {
    return guarded_ptr(-1.0, [&]() -> double {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return result->mMaxErrorApplied;
    });
}

void mio_decimate_result_free(mio_decimate_result* result) {
    delete result;
}

namespace {

mio_decimate_volume_result* capi_decimate_volume(
    const mio_mesh* pMesh, const meshioplusplus::DecimateVolumeOptions& rOptions) {
    if (!pMesh)
        throw meshioplusplus::ReadError("meshio++: mesh is NULL");
    meshioplusplus::DecimateVolumeResult r =
        meshioplusplus::decimate_volume(pMesh->mMesh, rOptions);
    auto* out = new mio_decimate_volume_result{};
    out->mMesh = mio_mesh{std::move(r.mMesh)};
    out->mPointMap = std::move(r.mPointMap);
    out->mCellMaps = std::move(r.mCellMaps);
    out->mTetsRemoved = r.mTetsRemoved;
    out->mPointsRemoved = r.mPointsRemoved;
    out->mCollapsesRejected = r.mCollapsesRejected;
    out->mMaxErrorApplied = r.mMaxErrorApplied;
    return out;
}

}  // namespace

mio_decimate_volume_result* mio_decimate_volume(const mio_mesh* mesh, double target_ratio,
                                                int64_t target_cells, double max_error,
                                                const char* placement, int preserve_boundary,
                                                int preserve_features, double feature_angle) {
    return guarded_ptr(
        static_cast<mio_decimate_volume_result*>(nullptr), [&]() -> mio_decimate_volume_result* {
            meshioplusplus::DecimateVolumeOptions options;
            options.mTargetRatio = target_ratio;
            options.mTargetCells = target_cells;
            options.mMaxError = max_error;
            options.mPlacement =
                meshioplusplus::decimate_placement_from_name(placement ? placement : "optimal");
            options.mPreserveBoundary = preserve_boundary != 0;
            options.mPreserveFeatures = preserve_features != 0;
            options.mFeatureAngleDeg = feature_angle;
            return capi_decimate_volume(mesh, options);
        });
}

void mio_decimate_volume_opts_init(mio_decimate_volume_opts* opts) {
    if (!opts)
        return;
    *opts = mio_decimate_volume_opts{};
    opts->target_ratio = -1.0;
    opts->target_cells = -1;
    opts->max_error = -1.0;
    opts->placement = nullptr;  // "optimal"
    opts->preserve_boundary = 0;
    opts->preserve_features = 1;
    opts->feature_angle = 30.0;
}

mio_decimate_volume_result* mio_decimate_volume_ex(const mio_mesh* mesh,
                                                   const mio_decimate_volume_opts* opts) {
    return guarded_ptr(
        static_cast<mio_decimate_volume_result*>(nullptr), [&]() -> mio_decimate_volume_result* {
            if (!opts)
                throw meshioplusplus::ReadError("meshio++: opts is NULL");
            if (!mesh)
                throw meshioplusplus::ReadError("meshio++: mesh is NULL");
            meshioplusplus::DecimateVolumeOptions options;
            options.mTargetRatio = opts->target_ratio;
            options.mTargetCells = opts->target_cells;
            options.mMaxError = opts->max_error;
            options.mPlacement = meshioplusplus::decimate_placement_from_name(
                opts->placement ? opts->placement : "optimal");
            options.mPreserveBoundary = opts->preserve_boundary != 0;
            options.mPreserveFeatures = opts->preserve_features != 0;
            options.mFeatureAngleDeg = opts->feature_angle;
            options.mFrozen =
                capi_frozen_mask(mesh, opts->frozen, opts->num_frozen, "decimate_volume");
            return capi_decimate_volume(mesh, options);
        });
}

static_assert(sizeof(mio_decimate_volume_opts) == 96,
              "mio_decimate_volume_opts grew outside its reserved tail");

const mio_mesh* mio_decimate_volume_result_mesh(const mio_decimate_volume_result* result) {
    return guarded_ptr(static_cast<const mio_mesh*>(nullptr), [&]() -> const mio_mesh* {
        if (!result)
            return nullptr;
        return &result->mMesh;
    });
}

mio_mesh* mio_decimate_volume_result_take_mesh(mio_decimate_volume_result* result) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return new mio_mesh{std::move(result->mMesh.mMesh)};
    });
}

mio_status mio_decimate_volume_result_point_map(const mio_decimate_volume_result* result,
                                                const void** data, mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        const NDArray& a = result->mPointMap;
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

int64_t mio_decimate_volume_result_num_cell_maps(const mio_decimate_volume_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return static_cast<int64_t>(result->mCellMaps.size());
    });
}

mio_status mio_decimate_volume_result_cell_map(const mio_decimate_volume_result* result,
                                               int64_t block, const void** data, mio_dtype* dtype,
                                               int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (block < 0 || static_cast<std::size_t>(block) >= result->mCellMaps.size())
            return fail(MIO_ERR_NOT_FOUND, "meshio++: cell-map block index out of range");
        const NDArray& a = result->mCellMaps[static_cast<std::size_t>(block)];
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

int64_t mio_decimate_volume_result_tets_removed(const mio_decimate_volume_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return result->mTetsRemoved;
    });
}

int64_t mio_decimate_volume_result_points_removed(const mio_decimate_volume_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return result->mPointsRemoved;
    });
}

int64_t mio_decimate_volume_result_collapses_rejected(const mio_decimate_volume_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return result->mCollapsesRejected;
    });
}

double mio_decimate_volume_result_max_error_applied(const mio_decimate_volume_result* result) {
    return guarded_ptr(-1.0, [&]() -> double {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return result->mMaxErrorApplied;
    });
}

void mio_decimate_volume_result_free(mio_decimate_volume_result* result) {
    delete result;
}

namespace {

// Shared option assembly for mio_partition / mio_partition_labels. NULL strings
// mean "auto" / "eco" / unweighted (the *_from_name parsers accept "").
meshioplusplus::PartitionOptions partition_options_from_c(int nparts, const char* method,
                                                          double imbalance, const char* mode,
                                                          int seed, int record_ids,
                                                          int ghost_layers,
                                                          const char* weights_key) {
    meshioplusplus::PartitionOptions options;
    options.mNParts = nparts;
    options.mMethod = meshioplusplus::partition_method_from_name(method ? method : "");
    options.mImbalance = imbalance;
    options.mMode = meshioplusplus::partition_mode_from_name(mode ? mode : "");
    options.mSeed = seed;
    options.mRecordIds = record_ids != 0;
    options.mGhostLayers = ghost_layers;
    options.mWeightsKey = weights_key ? weights_key : "";
    return options;
}

}  // namespace

mio_partition_result* mio_partition(const mio_mesh* mesh, int nparts, const char* method,
                                    double imbalance, const char* mode, int seed, int record_ids,
                                    int ghost_layers, const char* weights_key) {
    return guarded_ptr(static_cast<mio_partition_result*>(nullptr), [&]() -> mio_partition_result* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::PartitionOptions options = partition_options_from_c(
            nparts, method, imbalance, mode, seed, record_ids, ghost_layers, weights_key);
        meshioplusplus::PartitionResult r = meshioplusplus::partition(mesh->mMesh, options);
        auto* out = new mio_partition_result{};
        out->mPieces.reserve(r.mPieces.size());
        for (meshioplusplus::PartitionPiece& p : r.mPieces) {
            mio_partition_result::Piece piece;
            piece.mMesh = mio_mesh{std::move(p.mMesh)};
            piece.mPointMap = std::move(p.mPointMap);
            piece.mCellMaps = std::move(p.mCellMaps);
            out->mPieces.push_back(std::move(piece));
        }
        return out;
    });
}

int64_t mio_partition_result_num_pieces(const mio_partition_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return static_cast<int64_t>(result->mPieces.size());
    });
}

int mio_partition_result_part_id(const mio_partition_result* result, int64_t index) {
    return guarded_ptr(-1, [&]() -> int {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mPieces.size())
            throw meshioplusplus::ReadError("meshio++: piece index out of range");
        return static_cast<int>(index);
    });
}

const mio_mesh* mio_partition_result_mesh(const mio_partition_result* result, int64_t index) {
    return guarded_ptr(static_cast<const mio_mesh*>(nullptr), [&]() -> const mio_mesh* {
        if (!result)
            return nullptr;
        if (index < 0 || static_cast<std::size_t>(index) >= result->mPieces.size())
            return nullptr;
        return &result->mPieces[static_cast<std::size_t>(index)].mMesh;
    });
}

mio_mesh* mio_partition_result_take_mesh(mio_partition_result* result, int64_t index) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mPieces.size())
            throw meshioplusplus::ReadError("meshio++: piece index out of range");
        return new mio_mesh{
            std::move(result->mPieces[static_cast<std::size_t>(index)].mMesh.mMesh)};
    });
}

mio_status mio_partition_result_point_map(const mio_partition_result* result, int64_t index,
                                          const void** data, mio_dtype* dtype, int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mPieces.size())
            return fail(MIO_ERR_NOT_FOUND, "meshio++: piece index out of range");
        const NDArray& a = result->mPieces[static_cast<std::size_t>(index)].mPointMap;
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

int64_t mio_partition_result_num_cell_maps(const mio_partition_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        if (result->mPieces.empty())
            return 0;
        return static_cast<int64_t>(result->mPieces.front().mCellMaps.size());
    });
}

mio_status mio_partition_result_cell_map(const mio_partition_result* result, int64_t index,
                                         int64_t block, const void** data, mio_dtype* dtype,
                                         int64_t* n) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mPieces.size())
            return fail(MIO_ERR_NOT_FOUND, "meshio++: piece index out of range");
        const mio_partition_result::Piece& piece = result->mPieces[static_cast<std::size_t>(index)];
        if (block < 0 || static_cast<std::size_t>(block) >= piece.mCellMaps.size())
            return fail(MIO_ERR_NOT_FOUND, "meshio++: cell-map block index out of range");
        const NDArray& a = piece.mCellMaps[static_cast<std::size_t>(block)];
        if (data)
            *data = a.Data();
        if (dtype)
            *dtype = from_dtype(a.Dtype());
        if (n)
            *n = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        return MIO_OK;
    });
}

void mio_partition_result_free(mio_partition_result* result) {
    delete result;
}

mio_status mio_partition_labels(const mio_mesh* mesh, int nparts, const char* method,
                                double imbalance, const char* mode, int seed,
                                const char* weights_key, int64_t* labels, int64_t labels_size) {
    return guarded([&]() -> mio_status {
        if (!mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh is NULL");
        if (!labels)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: labels is NULL");
        meshioplusplus::PartitionOptions options =
            partition_options_from_c(nparts, method, imbalance, mode, seed, /*record_ids=*/0,
                                     /*ghost_layers=*/0, weights_key);
        std::vector<NDArray> per_block = meshioplusplus::partition_labels(mesh->mMesh, options);
        int64_t total = 0;
        for (const NDArray& a : per_block)
            total += a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
        if (labels_size != total)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: labels_size (" +
                                                 std::to_string(labels_size) +
                                                 ") does not equal the mesh's total cell count (" +
                                                 std::to_string(total) + ")");
        int64_t offset = 0;
        for (const NDArray& a : per_block) {
            const int64_t rows = a.Shape().empty() ? 0 : static_cast<int64_t>(a.Shape()[0]);
            std::memcpy(labels + offset, a.Data(),
                        static_cast<std::size_t>(rows) * sizeof(int64_t));
            offset += rows;
        }
        return MIO_OK;
    });
}

mio_status mio_stats(const mio_mesh* mesh, mio_stats_report* out) {
    return guarded([&]() -> mio_status {
        if (!mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh is NULL");
        if (!out)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: out is NULL");
        meshioplusplus::StatsReport r = meshioplusplus::compute_stats(mesh->mMesh);
        out->num_points = r.mNumPoints;
        out->num_cells = r.mNumCells;
        for (int d = 0; d < 3; ++d) {
            out->bbox_min[d] = r.mBBoxMin[d];
            out->bbox_max[d] = r.mBBoxMax[d];
            out->extent[d] = r.mExtent[d];
            out->centroid[d] = r.mCentroid[d];
        }
        out->total_area = r.mTotalArea;
        out->signed_volume = r.mSignedVolume;
        out->unsigned_volume = r.mUnsignedVolume;
        out->num_inverted = r.mNumInverted;
        return MIO_OK;
    });
}

/* --- data operations ---------------------------------------------------- */

/* The C enums duplicate the C++ ones, so pin every value: drift becomes a
 * compile error rather than a silently wrong argument (the same discipline the
 * MIO_CELL_TYPES X-macro applies to the cell-type list). */
static_assert(static_cast<int>(meshioplusplus::DataLocation::Point) == MIO_DATA_POINT, "");
static_assert(static_cast<int>(meshioplusplus::DataLocation::Cell) == MIO_DATA_CELL, "");
static_assert(static_cast<int>(meshioplusplus::DataLocation::Field) == MIO_DATA_FIELD, "");
static_assert(static_cast<int>(meshioplusplus::CellPointWeight::Uniform) == MIO_WEIGHT_UNIFORM, "");
static_assert(static_cast<int>(meshioplusplus::CellPointWeight::Measure) == MIO_WEIGHT_MEASURE, "");
static_assert(static_cast<int>(meshioplusplus::ConditionMode::Clamp) == MIO_COND_CLAMP, "");
static_assert(static_cast<int>(meshioplusplus::ConditionMode::Normalize) == MIO_COND_NORMALIZE, "");
static_assert(static_cast<int>(meshioplusplus::ConditionMode::Standardize) == MIO_COND_STANDARDIZE,
              "");
static_assert(static_cast<int>(meshioplusplus::ConditionScope::Component) == MIO_SCOPE_COMPONENT,
              "");
static_assert(static_cast<int>(meshioplusplus::ConditionScope::Magnitude) == MIO_SCOPE_MAGNITUDE,
              "");
static_assert(static_cast<int>(meshioplusplus::NanPolicy::Ignore) == MIO_NAN_IGNORE, "");
static_assert(static_cast<int>(meshioplusplus::NanPolicy::Replace) == MIO_NAN_REPLACE, "");
static_assert(static_cast<int>(meshioplusplus::NanPolicy::Fail) == MIO_NAN_FAIL, "");

namespace {

/* Turn the (names, count) pair the ABI uses into a vector, rejecting a NULL
 * entry rather than dereferencing it. */
std::vector<std::string> data_name_list(const char* const* pNames, int64_t count) {
    if (count < 0)
        throw meshioplusplus::ReadError("meshio++: name count is negative");
    if (count > 0 && !pNames)
        throw meshioplusplus::ReadError("meshio++: names is NULL but count > 0");
    std::vector<std::string> out;
    out.reserve(static_cast<std::size_t>(count));
    for (int64_t i = 0; i < count; ++i) {
        if (!pNames[i])
            throw meshioplusplus::ReadError("meshio++: names[" + std::to_string(i) + "] is NULL");
        out.emplace_back(pNames[i]);
    }
    return out;
}

meshioplusplus::DataLocation data_location_of(mio_data_location loc) {
    switch (loc) {
        case MIO_DATA_POINT:
            return meshioplusplus::DataLocation::Point;
        case MIO_DATA_CELL:
            return meshioplusplus::DataLocation::Cell;
        case MIO_DATA_FIELD:
            return meshioplusplus::DataLocation::Field;
    }
    throw meshioplusplus::ReadError("meshio++: invalid data location");
}

}  // namespace

mio_mesh* mio_data_drop(const mio_mesh* mesh, mio_data_location location, const char* const* names,
                        int64_t count, int ignore_missing) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return new mio_mesh{meshioplusplus::data_drop(mesh->mMesh, data_location_of(location),
                                                      data_name_list(names, count),
                                                      ignore_missing != 0)};
    });
}

mio_mesh* mio_data_keep(const mio_mesh* mesh, mio_data_location location, const char* const* names,
                        int64_t count, int ignore_missing) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return new mio_mesh{meshioplusplus::data_keep(mesh->mMesh, data_location_of(location),
                                                      data_name_list(names, count),
                                                      ignore_missing != 0)};
    });
}

mio_mesh* mio_data_rename(const mio_mesh* mesh, mio_data_location location, const char* from_name,
                          const char* to_name) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !from_name || !to_name)
            throw meshioplusplus::ReadError("meshio++: mesh/from_name/to_name is NULL");
        return new mio_mesh{meshioplusplus::data_rename(mesh->mMesh, data_location_of(location),
                                                        from_name, to_name)};
    });
}

mio_mesh* mio_sets_to_data(const mio_mesh* mesh, mio_data_location location, const char* data_name,
                           const char* join_char, const char* const* order, int64_t count) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw std::invalid_argument("meshio++: sets_to_data: mesh is NULL");
        return new mio_mesh{meshioplusplus::sets_to_data(
            mesh->mMesh, data_location_of(location),
            data_name ? std::optional<std::string>(data_name) : std::nullopt,
            join_char ? join_char : "-", data_name_list(order, count))};
    });
}

mio_mesh* mio_data_to_sets(const mio_mesh* mesh, mio_data_location location, const char* key) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !key)
            throw std::invalid_argument("meshio++: data_to_sets: mesh/key is NULL");
        return new mio_mesh{
            meshioplusplus::data_to_sets(mesh->mMesh, data_location_of(location), key)};
    });
}

mio_mesh* mio_data_point_to_cell(const mio_mesh* mesh, const char* const* names, int64_t count,
                                 const char* suffix) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::DataAverageOptions opts;
        opts.names = data_name_list(names, count);
        opts.suffix = suffix ? suffix : "";
        return new mio_mesh{meshioplusplus::point_data_to_cell_data(mesh->mMesh, opts)};
    });
}

mio_mesh* mio_data_cell_to_point(const mio_mesh* mesh, const char* const* names, int64_t count,
                                 mio_cell_point_weight weight, const char* suffix) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::DataAverageOptions opts;
        opts.names = data_name_list(names, count);
        opts.weight = weight == MIO_WEIGHT_MEASURE ? meshioplusplus::CellPointWeight::Measure
                                                   : meshioplusplus::CellPointWeight::Uniform;
        opts.suffix = suffix ? suffix : "";
        return new mio_mesh{meshioplusplus::cell_data_to_point_data(mesh->mMesh, opts)};
    });
}

mio_mesh* mio_data_calc(const mio_mesh* mesh, const char* expression, mio_data_location location,
                        const char* output_name, int overwrite) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !expression || !output_name)
            throw meshioplusplus::ReadError("meshio++: mesh/expression/output_name is NULL");
        meshioplusplus::DataCalcOptions opts;
        opts.location = data_location_of(location);
        opts.output = output_name;
        opts.overwrite = overwrite != 0;
        return new mio_mesh{meshioplusplus::data_calc(mesh->mMesh, expression, opts)};
    });
}

mio_mesh* mio_data_condition(const mio_mesh* mesh, mio_data_location location,
                             const char* const* names, int64_t count, mio_condition_mode mode,
                             double lo, double hi, mio_condition_scope scope,
                             mio_nan_policy nan_policy, double nan_replacement,
                             const char* suffix) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::DataConditionOptions opts;
        opts.location = data_location_of(location);
        opts.names = data_name_list(names, count);
        opts.mode = static_cast<meshioplusplus::ConditionMode>(mode);
        opts.scope = static_cast<meshioplusplus::ConditionScope>(scope);
        opts.lo = lo;
        opts.hi = hi;
        opts.nan_policy = static_cast<meshioplusplus::NanPolicy>(nan_policy);
        opts.nan_replacement = nan_replacement;
        opts.suffix = suffix ? suffix : "";
        return new mio_mesh{meshioplusplus::data_condition(mesh->mMesh, opts)};
    });
}

mio_mesh* mio_tensor_invariants(const mio_mesh* mesh, mio_data_location location,
                                const char* const* names, int64_t count, int outputs,
                                const char* prefix, const char* suffix, int overwrite) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::TensorInvariantsOptions opts;
        opts.location = data_location_of(location);
        opts.names = data_name_list(names, count);
        opts.outputs =
            static_cast<meshioplusplus::TensorInvariant>(outputs == 0 ? MIO_TINV_ALL : outputs);
        opts.prefix = prefix ? prefix : "";
        opts.suffix = suffix ? suffix : "";
        opts.overwrite = overwrite != 0;
        return new mio_mesh{meshioplusplus::tensor_invariants(mesh->mMesh, opts)};
    });
}

mio_data_info* mio_data_info_create(const mio_mesh* mesh) {
    return guarded_ptr(static_cast<mio_data_info*>(nullptr), [&]() -> mio_data_info* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        return new mio_data_info{meshioplusplus::data_info(mesh->mMesh)};
    });
}

int64_t mio_data_info_count(const mio_data_info* info) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!info)
            throw meshioplusplus::ReadError("meshio++: info is NULL");
        return static_cast<int64_t>(info->mReport.mArrays.size());
    });
}

int64_t mio_data_info_name(const mio_data_info* info, int64_t index, char* buf, int64_t buflen) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!info)
            throw meshioplusplus::ReadError("meshio++: info is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= info->mReport.mArrays.size())
            throw meshioplusplus::ReadError("meshio++: data array index out of range");
        return copy_string(info->mReport.mArrays[static_cast<std::size_t>(index)].mName, buf,
                           buflen);
    });
}

mio_status mio_data_info_entry(const mio_data_info* info, int64_t index, mio_data_array_info* out) {
    return guarded([&]() -> mio_status {
        if (!info)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: info is NULL");
        if (!out)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: out is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= info->mReport.mArrays.size())
            return fail(MIO_ERR_INVALID_ARG, "meshio++: data array index out of range");
        const meshioplusplus::DataArrayInfo& a =
            info->mReport.mArrays[static_cast<std::size_t>(index)];
        out->location = static_cast<int>(a.mLocation);
        out->dtype = static_cast<int>(a.mDtype);
        out->num_blocks = a.mNumBlocks;
        out->num_entries = a.mNumEntries;
        out->num_components = a.mNumComponents;
        out->num_values = a.mNumValues;
        out->min = a.mMin;
        out->max = a.mMax;
        out->mean = a.mMean;
        out->num_nan = a.mNumNan;
        out->num_inf = a.mNumInf;
        out->num_finite = a.mNumFinite;
        out->inconsistent_blocks = a.mInconsistentBlocks ? 1 : 0;
        return MIO_OK;
    });
}

mio_status mio_data_info_component(const mio_data_info* info, int64_t index, int64_t comp,
                                   double* min, double* max, double* mean) {
    return guarded([&]() -> mio_status {
        if (!info)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: info is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= info->mReport.mArrays.size())
            return fail(MIO_ERR_INVALID_ARG, "meshio++: data array index out of range");
        const meshioplusplus::DataArrayInfo& a =
            info->mReport.mArrays[static_cast<std::size_t>(index)];
        if (comp < 0 || static_cast<std::size_t>(comp) >= a.mMinPerComponent.size())
            return fail(MIO_ERR_INVALID_ARG, "meshio++: component index out of range");
        const std::size_t k = static_cast<std::size_t>(comp);
        if (min)
            *min = a.mMinPerComponent[k];
        if (max)
            *max = a.mMaxPerComponent[k];
        if (mean)
            *mean = a.mMeanPerComponent[k];
        return MIO_OK;
    });
}

void mio_data_info_free(mio_data_info* info) {
    delete info;
}

mio_data_integrate* mio_data_integrate_create(const mio_mesh* mesh, const char* const* names,
                                              int64_t count) {
    return guarded_ptr(static_cast<mio_data_integrate*>(nullptr), [&]() -> mio_data_integrate* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::DataIntegrateOptions opts;
        if (names) {
            for (int64_t i = 0; i < count; ++i)
                opts.mArrayNames.emplace_back(names[i] ? names[i] : "");
        }
        return new mio_data_integrate{meshioplusplus::data_integrate(mesh->mMesh, opts)};
    });
}

int64_t mio_data_integrate_count(const mio_data_integrate* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return static_cast<int64_t>(result->mReport.mArrays.size());
    });
}

int64_t mio_data_integrate_name(const mio_data_integrate* result, int64_t index, char* buf,
                                int64_t buflen) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mReport.mArrays.size())
            throw meshioplusplus::ReadError("meshio++: array index out of range");
        return copy_string(result->mReport.mArrays[static_cast<std::size_t>(index)].mName, buf,
                           buflen);
    });
}

namespace {
mio_status fill_field_integral_info(const meshioplusplus::FieldIntegralRegion& rRegion,
                                    mio_field_integral_info* out) {
    if (!out)
        return fail(MIO_ERR_INVALID_ARG, "meshio++: out is NULL");
    out->num_components = static_cast<int64_t>(rRegion.mTotalPerComponent.size());
    out->num_cells = rRegion.mNumCells;
    out->num_skipped = rRegion.mNumSkipped;
    return MIO_OK;
}

mio_status fill_field_integral_component(const meshioplusplus::FieldIntegralRegion& rRegion,
                                         int64_t comp, double* total, double* mean,
                                         double* domain_measure, int64_t* num_nan) {
    if (comp < 0 || static_cast<std::size_t>(comp) >= rRegion.mTotalPerComponent.size())
        return fail(MIO_ERR_INVALID_ARG, "meshio++: component index out of range");
    const std::size_t k = static_cast<std::size_t>(comp);
    if (total)
        *total = rRegion.mTotalPerComponent[k];
    if (mean)
        *mean = rRegion.mMeanPerComponent[k];
    if (domain_measure)
        *domain_measure = rRegion.mDomainMeasurePerComponent[k];
    if (num_nan)
        *num_nan = rRegion.mNumNanPerComponent[k];
    return MIO_OK;
}
}  // namespace

mio_status mio_data_integrate_entry(const mio_data_integrate* result, int64_t index,
                                    mio_field_integral_info* out) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mReport.mArrays.size())
            return fail(MIO_ERR_INVALID_ARG, "meshio++: array index out of range");
        return fill_field_integral_info(
            result->mReport.mArrays[static_cast<std::size_t>(index)].mDomain, out);
    });
}

mio_status mio_data_integrate_component(const mio_data_integrate* result, int64_t index,
                                        int64_t comp, double* total, double* mean,
                                        double* domain_measure, int64_t* num_nan) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mReport.mArrays.size())
            return fail(MIO_ERR_INVALID_ARG, "meshio++: array index out of range");
        return fill_field_integral_component(
            result->mReport.mArrays[static_cast<std::size_t>(index)].mDomain, comp, total, mean,
            domain_measure, num_nan);
    });
}

int64_t mio_data_integrate_region_count(const mio_data_integrate* result, int64_t index) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mReport.mArrays.size())
            throw meshioplusplus::ReadError("meshio++: array index out of range");
        return static_cast<int64_t>(
            result->mReport.mArrays[static_cast<std::size_t>(index)].mRegions.size());
    });
}

int64_t mio_data_integrate_region_name(const mio_data_integrate* result, int64_t index,
                                       int64_t region, char* buf, int64_t buflen) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mReport.mArrays.size())
            throw meshioplusplus::ReadError("meshio++: array index out of range");
        const auto& regions = result->mReport.mArrays[static_cast<std::size_t>(index)].mRegions;
        if (region < 0 || static_cast<std::size_t>(region) >= regions.size())
            throw meshioplusplus::ReadError("meshio++: region index out of range");
        return copy_string(regions[static_cast<std::size_t>(region)].mName, buf, buflen);
    });
}

mio_status mio_data_integrate_region_entry(const mio_data_integrate* result, int64_t index,
                                           int64_t region, mio_field_integral_info* out) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mReport.mArrays.size())
            return fail(MIO_ERR_INVALID_ARG, "meshio++: array index out of range");
        const auto& regions = result->mReport.mArrays[static_cast<std::size_t>(index)].mRegions;
        if (region < 0 || static_cast<std::size_t>(region) >= regions.size())
            return fail(MIO_ERR_INVALID_ARG, "meshio++: region index out of range");
        return fill_field_integral_info(regions[static_cast<std::size_t>(region)], out);
    });
}

mio_status mio_data_integrate_region_component(const mio_data_integrate* result, int64_t index,
                                               int64_t region, int64_t comp, double* total,
                                               double* mean, double* domain_measure,
                                               int64_t* num_nan) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= result->mReport.mArrays.size())
            return fail(MIO_ERR_INVALID_ARG, "meshio++: array index out of range");
        const auto& regions = result->mReport.mArrays[static_cast<std::size_t>(index)].mRegions;
        if (region < 0 || static_cast<std::size_t>(region) >= regions.size())
            return fail(MIO_ERR_INVALID_ARG, "meshio++: region index out of range");
        return fill_field_integral_component(regions[static_cast<std::size_t>(region)], comp, total,
                                             mean, domain_measure, num_nan);
    });
}

void mio_data_integrate_free(mio_data_integrate* result) {
    delete result;
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

mio_status mio_diff(const mio_mesh* a, const mio_mesh* b, double atol, double rtol, int unordered,
                    mio_diff_result** out) {
    return guarded([&]() -> mio_status {
        if (!a || !b)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh is NULL");
        if (!out)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: out is NULL");
        meshioplusplus::DiffOptions opts;
        opts.atol = atol;
        opts.rtol = rtol;
        opts.unordered = unordered != 0;
        *out = new mio_diff_result{meshioplusplus::diff(a->mMesh, b->mMesh, opts)};
        return MIO_OK;
    });
}

mio_status mio_meshes_equal(const mio_mesh* a, const mio_mesh* b, double atol, double rtol,
                            int unordered, int* out_equal) {
    return guarded([&]() -> mio_status {
        if (!a || !b)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh is NULL");
        meshioplusplus::DiffOptions opts;
        opts.atol = atol;
        opts.rtol = rtol;
        opts.unordered = unordered != 0;
        const bool eq = meshioplusplus::diff(a->mMesh, b->mMesh, opts).mVerdict !=
                        meshioplusplus::DiffVerdict::Different;
        if (out_equal)
            *out_equal = eq ? 1 : 0;
        return MIO_OK;
    });
}

mio_diff_verdict mio_diff_result_verdict(const mio_diff_result* result) {
    if (!result)
        return MIO_DIFF_DIFFERENT;
    switch (result->mReport.mVerdict) {
        case meshioplusplus::DiffVerdict::Identical:
            return MIO_DIFF_IDENTICAL;
        case meshioplusplus::DiffVerdict::EqualWithinTolerance:
            return MIO_DIFF_EQUAL_WITHIN_TOLERANCE;
        case meshioplusplus::DiffVerdict::Different:
            return MIO_DIFF_DIFFERENT;
    }
    return MIO_DIFF_DIFFERENT;
}

mio_status mio_diff_result_point_summary(const mio_diff_result* result, double* max_abs_error,
                                         double* max_rel_error, int64_t* worst_index,
                                         int64_t* num_exceeding) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        const meshioplusplus::ArrayDiff& p = result->mReport.mPoints;
        if (max_abs_error)
            *max_abs_error = p.mMaxAbsError;
        if (max_rel_error)
            *max_rel_error = p.mMaxRelError;
        if (worst_index)
            *worst_index = p.mWorstIndex;
        if (num_exceeding)
            *num_exceeding = p.mNumExceeding;
        return MIO_OK;
    });
}

int64_t mio_diff_result_num_block_diffs(const mio_diff_result* result) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: result is NULL");
        return static_cast<int64_t>(result->mReport.mBlocks.size());
    });
}

mio_status mio_diff_result_block(const mio_diff_result* result, int64_t block, int* type_mismatch,
                                 int* count_mismatch, int64_t* conn_mismatch_count) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: result is NULL");
        if (block < 0 || static_cast<std::size_t>(block) >= result->mReport.mBlocks.size())
            return fail(MIO_ERR_NOT_FOUND, "meshio++: block index out of range");
        const meshioplusplus::BlockDiff& bd =
            result->mReport.mBlocks[static_cast<std::size_t>(block)];
        if (type_mismatch)
            *type_mismatch = bd.mTypeMismatch ? 1 : 0;
        if (count_mismatch)
            *count_mismatch = bd.mCountMismatch ? 1 : 0;
        if (conn_mismatch_count)
            *conn_mismatch_count = bd.mConnMismatchCount;
        return MIO_OK;
    });
}

void mio_diff_result_free(mio_diff_result* result) {
    delete result;
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

namespace {

// Validate one CSR offsets array of `count + 1` entries spanning [0, total].
// A malformed offsets array is the one way a caller can make the ragged setters
// read out of bounds, so every property is checked rather than assumed.
mio_status poly_check_offsets(const std::int64_t* pOffsets, std::int64_t Count, std::int64_t Total,
                              const char* pWhat) {
    if (!pOffsets)
        return fail(MIO_ERR_INVALID_ARG, "meshio++: " + std::string(pWhat) + " is NULL (expected " +
                                             std::to_string(Count + 1) + " entries)");
    if (pOffsets[0] != 0)
        return fail(MIO_ERR_INVALID_ARG, "meshio++: " + std::string(pWhat) + "[0] must be 0, got " +
                                             std::to_string(pOffsets[0]));
    for (std::int64_t i = 0; i < Count; ++i)
        if (pOffsets[i + 1] < pOffsets[i])
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: " + std::string(pWhat) + " must be non-decreasing (entry " +
                            std::to_string(i + 1) + " is " + std::to_string(pOffsets[i + 1]) +
                            " after " + std::to_string(pOffsets[i]) + ")");
    if (pOffsets[Count] != Total)
        return fail(MIO_ERR_INVALID_ARG,
                    "meshio++: " + std::string(pWhat) + "[" + std::to_string(Count) + "] must be " +
                        std::to_string(Total) + ", got " + std::to_string(pOffsets[Count]));
    return MIO_OK;
}

// Node ids are checked non-negative but deliberately NOT range-checked against
// NumPoints(): points may legitimately be assigned after cells, and no other
// setter on this ABI range-checks either.
mio_status poly_check_nodes(const std::int64_t* pNodes, std::int64_t Count) {
    if (Count > 0 && !pNodes)
        return fail(MIO_ERR_INVALID_ARG, "meshio++: nodes is NULL");
    for (std::int64_t i = 0; i < Count; ++i)
        if (pNodes[i] < 0)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: node id at index " + std::to_string(i) +
                                                 " is negative (" + std::to_string(pNodes[i]) +
                                                 "); connectivity is 0-based");
    return MIO_OK;
}

}  // namespace

mio_status mio_mesh_add_polygon_block(mio_mesh* mesh, const char* cell_type, int64_t num_cells,
                                      const int64_t* row_offsets, const int64_t* nodes,
                                      int64_t num_nodes) {
    return guarded([&]() -> mio_status {
        if (!mesh || !cell_type || num_cells < 0 || num_nodes < 0)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: bad mio_mesh_add_polygon_block argument");
        if (mio_status s = poly_check_offsets(row_offsets, num_cells, num_nodes, "row_offsets");
            s != MIO_OK)
            return s;
        if (mio_status s = poly_check_nodes(nodes, num_nodes); s != MIO_OK)
            return s;
        // The C arrays are already the CSR the mesh stores: copy them in whole.
        mesh->mMesh.AddPolygonBlock(
            cell_type, std::vector<std::int64_t>(nodes, nodes + num_nodes),
            std::vector<std::int64_t>(row_offsets, row_offsets + num_cells + 1));
        return MIO_OK;
    });
}

mio_status mio_mesh_add_polyhedron_block(mio_mesh* mesh, const char* cell_type, int64_t num_cells,
                                         const int64_t* cell_offsets, int64_t num_faces,
                                         const int64_t* face_offsets, const int64_t* nodes,
                                         int64_t num_nodes) {
    return guarded([&]() -> mio_status {
        if (!mesh || !cell_type || num_cells < 0 || num_faces < 0 || num_nodes < 0)
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: bad mio_mesh_add_polyhedron_block argument");
        if (mio_status s = poly_check_offsets(cell_offsets, num_cells, num_faces, "cell_offsets");
            s != MIO_OK)
            return s;
        if (mio_status s = poly_check_offsets(face_offsets, num_faces, num_nodes, "face_offsets");
            s != MIO_OK)
            return s;
        if (mio_status s = poly_check_nodes(nodes, num_nodes); s != MIO_OK)
            return s;
        // The C arrays are already the CSR the mesh stores (face_offsets over
        // the nodes, cell_offsets over the faces): copy them in whole.
        mesh->mMesh.AddPolyhedronBlock(
            cell_type, std::vector<std::int64_t>(nodes, nodes + num_nodes),
            std::vector<std::int64_t>(face_offsets, face_offsets + num_faces + 1),
            std::vector<std::int64_t>(cell_offsets, cell_offsets + num_cells + 1));
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

namespace {

// The one place a CellView is unrolled into flat CSR. Shared by
// mio_mesh_cell_block_info_ex (which needs the totals) and mio_poly_conn_create
// (which needs the arrays), so the two can never report different shapes.
// Goes through the uniform API only, hence it compiles under all three
// backends -- see the mio_poly_conn comment above.
void poly_flatten(const meshioplusplus::Mesh::CellView& rView, mio_poly_conn& rOut) {
    const std::size_t ncells = rView.NumCells();
    rOut.mIsPolyhedron = rView.IsPolyhedron();
    rOut.mNumCells = static_cast<std::int64_t>(ncells);
    rOut.mFaceOffsets.push_back(0);
    if (rOut.mIsPolyhedron) {
        rOut.mCellOffsets.push_back(0);
        for (std::size_t c = 0; c < ncells; ++c) {
            for (std::size_t f = 0; f < rView.NumFaces(c); ++f) {
                const auto face = rView.Face(c, f);
                rOut.mNodes.insert(rOut.mNodes.end(), face.first, face.first + face.second);
                rOut.mFaceOffsets.push_back(static_cast<std::int64_t>(rOut.mNodes.size()));
            }
            rOut.mCellOffsets.push_back(static_cast<std::int64_t>(rOut.mFaceOffsets.size() - 1));
        }
    } else {
        // 1-level: one face per cell, so face_offsets IS the row-offsets array
        // and cell_offsets stays empty (reported as NULL).
        for (std::size_t c = 0; c < ncells; ++c) {
            const std::int64_t* p_row = rView.Row(c);
            rOut.mNodes.insert(rOut.mNodes.end(), p_row, p_row + rView.RowSize(c));
            rOut.mFaceOffsets.push_back(static_cast<std::int64_t>(rOut.mNodes.size()));
        }
    }
}

}  // namespace

mio_status mio_mesh_cell_block_info_ex(const mio_mesh* mesh, int64_t block,
                                       mio_cell_block_info* out) {
    return guarded([&]() -> mio_status {
        if (!mesh || !out)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh or out is NULL");
        if (!block_in_range(mesh, block))
            return fail(MIO_ERR_NOT_FOUND,
                        "meshio++: cell block " + std::to_string(block) + " out of range");
        const auto view = mesh->mMesh.Cells(static_cast<std::size_t>(block));
        *out = mio_cell_block_info{};
        out->num_cells = static_cast<int64_t>(view.NumCells());
        out->nodes_per_cell = static_cast<int64_t>(view.NodesPerCell());
        out->is_ragged = view.IsRagged() ? 1 : 0;
        out->is_polyhedron = view.IsPolyhedron() ? 1 : 0;
        if (!view.IsRagged()) {
            out->num_faces = out->num_cells;
            out->num_nodes = out->num_cells * out->nodes_per_cell;
        } else if (view.IsPolyhedron()) {
            for (std::size_t c = 0; c < view.NumCells(); ++c) {
                const std::size_t nf = view.NumFaces(c);
                out->num_faces += static_cast<int64_t>(nf);
                for (std::size_t f = 0; f < nf; ++f)
                    out->num_nodes += static_cast<int64_t>(view.Face(c, f).second);
            }
        } else {
            out->num_faces = out->num_cells;
            for (std::size_t c = 0; c < view.NumCells(); ++c)
                out->num_nodes += static_cast<int64_t>(view.RowSize(c));
        }
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
                        "meshio++: cell block " + std::to_string(block) +
                            " is ragged and has no rectangular connectivity to borrow; "
                            "read it with mio_poly_conn_create()");
        const NDArray& c = view.Conn();
        return array_out(c, conn, dtype, nullptr, nullptr);
    });
}

/* Ragged connectivity: an owning snapshot rather than a rule-3 borrow. See the
 * mio_poly_conn struct comment above and doc/polyhedra.md. */

mio_poly_conn* mio_poly_conn_create(const mio_mesh* mesh, int64_t block) {
    return guarded_ptr(static_cast<mio_poly_conn*>(nullptr), [&]() -> mio_poly_conn* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        if (!block_in_range(mesh, block))
            throw meshioplusplus::ReadError("meshio++: cell block " + std::to_string(block) +
                                            " out of range");
        const auto view = mesh->mMesh.Cells(static_cast<std::size_t>(block));
        if (!view.IsRagged())
            throw meshioplusplus::ReadError(
                "meshio++: cell block " + std::to_string(block) +
                " is rectangular; borrow it with mio_mesh_cell_block_conn()");
        // Filled by value first, so a throw mid-flatten cannot leak the handle.
        mio_poly_conn poly;
        poly_flatten(view, poly);
        return new mio_poly_conn(std::move(poly));
    });
}

mio_status mio_poly_conn_get_shape(const mio_poly_conn* poly, mio_poly_conn_shape* out) {
    return guarded([&]() -> mio_status {
        if (!poly || !out)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: poly or out is NULL");
        *out = mio_poly_conn_shape{};
        out->is_polyhedron = poly->mIsPolyhedron ? 1 : 0;
        out->num_cells = poly->mNumCells;
        out->num_faces = static_cast<int64_t>(poly->mFaceOffsets.size()) - 1;
        out->num_nodes = static_cast<int64_t>(poly->mNodes.size());
        return MIO_OK;
    });
}

namespace {

const int64_t* poly_array_out(const mio_poly_conn* pPoly, const std::vector<std::int64_t>& rVec,
                              int64_t* pCount) {
    if (!pPoly)
        throw meshioplusplus::ReadError("meshio++: poly is NULL");
    if (pCount)
        *pCount = static_cast<int64_t>(rVec.size());
    return rVec.empty() ? nullptr : rVec.data();
}

}  // namespace

const int64_t* mio_poly_conn_nodes(const mio_poly_conn* poly, int64_t* count) {
    return guarded_ptr(static_cast<const int64_t*>(nullptr), [&]() -> const int64_t* {
        return poly_array_out(poly, poly->mNodes, count);
    });
}

const int64_t* mio_poly_conn_face_offsets(const mio_poly_conn* poly, int64_t* count) {
    return guarded_ptr(static_cast<const int64_t*>(nullptr), [&]() -> const int64_t* {
        return poly_array_out(poly, poly->mFaceOffsets, count);
    });
}

const int64_t* mio_poly_conn_cell_offsets(const mio_poly_conn* poly, int64_t* count) {
    // Empty (so NULL, count 0) for a 1-level block: a NULL nobody can mistake
    // for data, rather than a synthesized identity that looks like information.
    return guarded_ptr(static_cast<const int64_t*>(nullptr), [&]() -> const int64_t* {
        return poly_array_out(poly, poly->mCellOffsets, count);
    });
}

void mio_poly_conn_free(mio_poly_conn* poly) {
    delete poly;
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

/* ---------------------------------------------------------------------
 * Named regions (see the header, and doc/regions.md)
 * --------------------------------------------------------------------- */

static_assert(static_cast<int>(meshioplusplus::RegionKind::Point) == MIO_REGION_POINT, "");
static_assert(static_cast<int>(meshioplusplus::RegionKind::Cell) == MIO_REGION_CELL, "");
static_assert(static_cast<int>(meshioplusplus::RegionKind::Side) == MIO_REGION_SIDE, "");

mio_regions* mio_regions_create(const mio_mesh* mesh) {
    return guarded_ptr(static_cast<mio_regions*>(nullptr), [&]() -> mio_regions* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        auto* out = new mio_regions{};
        out->mRegions.reserve(mesh->mMesh.NumRegions());
        for (std::size_t i = 0; i < mesh->mMesh.NumRegions(); ++i)
            out->mRegions.push_back(mesh->mMesh.Region(i));
        return out;
    });
}

int64_t mio_regions_count(const mio_regions* regions) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!regions)
            throw meshioplusplus::ReadError("meshio++: regions is NULL");
        return static_cast<int64_t>(regions->mRegions.size());
    });
}

int64_t mio_regions_name(const mio_regions* regions, int64_t index, char* buf, int64_t buflen) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!regions)
            throw meshioplusplus::ReadError("meshio++: regions is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= regions->mRegions.size())
            throw meshioplusplus::ReadError("meshio++: region index " + std::to_string(index) +
                                            " out of range");
        return copy_string(regions->mRegions[static_cast<std::size_t>(index)].mName, buf, buflen);
    });
}

mio_status mio_regions_info(const mio_regions* regions, int64_t index, mio_region_info* out) {
    return guarded([&]() -> mio_status {
        if (!regions || !out)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: regions/out is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= regions->mRegions.size())
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: region index " + std::to_string(index) + " out of range");
        const meshioplusplus::Region& r = regions->mRegions[static_cast<std::size_t>(index)];
        out->kind = static_cast<int32_t>(r.mKind);
        out->dim = static_cast<int32_t>(r.mDim);
        out->tag = r.mTag;
        out->num_entries = static_cast<int64_t>(r.NumEntries());
        out->stride = static_cast<int64_t>(r.Stride());
        return MIO_OK;
    });
}

const int64_t* mio_regions_entries(const mio_regions* regions, int64_t index, int64_t* count) {
    return guarded_ptr(static_cast<const int64_t*>(nullptr), [&]() -> const int64_t* {
        if (!regions)
            throw meshioplusplus::ReadError("meshio++: regions is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= regions->mRegions.size())
            throw meshioplusplus::ReadError("meshio++: region index " + std::to_string(index) +
                                            " out of range");
        const meshioplusplus::Region& r = regions->mRegions[static_cast<std::size_t>(index)];
        if (count)
            *count = static_cast<int64_t>(r.mEntries.Size());
        return r.Entries();
    });
}

void mio_regions_free(mio_regions* regions) {
    delete regions;
}

mio_status mio_mesh_add_region(mio_mesh* mesh, const char* name, mio_region_kind kind, int32_t dim,
                               int64_t tag, const int64_t* entries, int64_t count) {
    return guarded([&]() -> mio_status {
        if (!mesh || !name)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh/name is NULL");
        if (count < 0 || (count > 0 && !entries))
            return fail(MIO_ERR_INVALID_ARG, "meshio++: region entries is NULL");
        const auto region_kind = static_cast<meshioplusplus::RegionKind>(kind);
        if (region_kind == meshioplusplus::RegionKind::Side && count % 2 != 0)
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: a side region needs an even entry count "
                        "((cell, facet) pairs)");
        const auto n = static_cast<std::size_t>(count);
        std::vector<std::size_t> shape;
        if (region_kind == meshioplusplus::RegionKind::Side)
            shape = {n / 2, 2};
        else
            shape = {n};
        meshioplusplus::NDArray arr =
            meshioplusplus::NDArray::Uninit(meshioplusplus::DType::Int64, std::move(shape));
        for (std::size_t i = 0; i < n; ++i)
            arr.As<std::int64_t>()[i] = entries[i];
        mesh->mMesh.AddRegion(meshioplusplus::Region(name, region_kind, dim, tag, std::move(arr)));
        return MIO_OK;
    });
}

// ---- transient (time-series) XDMF ----------------------------------------

mio_exodus_series* mio_exodus_series_create(const char* path) {
    return guarded_ptr(static_cast<mio_exodus_series*>(nullptr), [&]() -> mio_exodus_series* {
        if (!path)
            throw std::invalid_argument("meshio++: Exodus series path is NULL");
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        return new mio_exodus_series(path);
#else
        throw meshioplusplus::WriteError("Exodus series requires -DMESHIOPLUSPLUS_WITH_NETCDF=ON");
#endif
    });
}

mio_status mio_exodus_series_write_points_cells(mio_exodus_series* series, const mio_mesh* mesh) {
    return guarded([&]() -> mio_status {
        if (!series || !mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: series/mesh is NULL");
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        series->mWriter.WritePointsCells(mesh->mMesh);
        return MIO_OK;
#else
        throw meshioplusplus::WriteError("Exodus series requires -DMESHIOPLUSPLUS_WITH_NETCDF=ON");
#endif
    });
}

mio_status mio_exodus_series_write_data(mio_exodus_series* series, double time,
                                        const mio_mesh* mesh) {
    return guarded([&]() -> mio_status {
        if (!series || !mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: series/mesh is NULL");
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        series->mWriter.WriteData(time, mesh->mMesh);
        return MIO_OK;
#else
        (void)time;
        throw meshioplusplus::WriteError("Exodus series requires -DMESHIOPLUSPLUS_WITH_NETCDF=ON");
#endif
    });
}

mio_status mio_exodus_series_flush(mio_exodus_series* series) {
    return guarded([&]() -> mio_status {
        if (!series)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: series is NULL");
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        series->mWriter.Flush();
        return MIO_OK;
#else
        throw meshioplusplus::WriteError("Exodus series requires -DMESHIOPLUSPLUS_WITH_NETCDF=ON");
#endif
    });
}

mio_status mio_exodus_series_finalize(mio_exodus_series* series) {
    return guarded([&]() -> mio_status {
        if (!series)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: series is NULL");
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        series->mWriter.Finalize();
        return MIO_OK;
#else
        throw meshioplusplus::WriteError("Exodus series requires -DMESHIOPLUSPLUS_WITH_NETCDF=ON");
#endif
    });
}

int64_t mio_exodus_series_num_steps(const mio_exodus_series* series) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!series)
            throw std::invalid_argument("meshio++: series is NULL");
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        return static_cast<std::int64_t>(series->mWriter.NumSteps());
#else
        throw meshioplusplus::WriteError("Exodus series requires -DMESHIOPLUSPLUS_WITH_NETCDF=ON");
#endif
    });
}

int32_t mio_exodus_series_finalized(const mio_exodus_series* series) {
    return guarded_ptr(std::int32_t(-1), [&]() -> std::int32_t {
        if (!series)
            throw std::invalid_argument("meshio++: series is NULL");
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
        return series->mWriter.Finalized();
#else
        throw meshioplusplus::WriteError("Exodus series requires -DMESHIOPLUSPLUS_WITH_NETCDF=ON");
#endif
    });
}

void mio_exodus_series_free(mio_exodus_series* series) {
    delete series;
}
//
// The one writer that is a handle rather than a (path, mesh) call: the mesh is
// written once and each step appended, so there is no single call for mio_write
// to be. The C++ object is held by value; every entry point stays inside
// guarded()/guarded_ptr() like the rest of this file.

mio_xdmf_series* mio_xdmf_series_create(const char* path, const char* data_format,
                                        int32_t gzip_level) {
    return guarded_ptr(static_cast<mio_xdmf_series*>(nullptr), [&]() -> mio_xdmf_series* {
        if (!path)
            throw meshioplusplus::WriteError("meshio++: path is NULL");
        const std::string fmt = data_format ? data_format : "HDF";
        return new mio_xdmf_series{
            meshioplusplus::XdmfTimeSeriesWriter(path, fmt, static_cast<int>(gzip_level))};
    });
}

void mio_xdmf_series_opts_init(mio_xdmf_series_opts* opts) {
    if (!opts)
        return;
    *opts = mio_xdmf_series_opts{};  // value-initialized: HDF, truncate, no flush
}

mio_xdmf_series* mio_xdmf_series_create_ex(const char* path, const mio_xdmf_series_opts* opts) {
    return guarded_ptr(static_cast<mio_xdmf_series*>(nullptr), [&]() -> mio_xdmf_series* {
        if (!path)
            throw meshioplusplus::WriteError("meshio++: path is NULL");
        mio_xdmf_series_opts defaults{};
        const mio_xdmf_series_opts& r_o = opts ? *opts : defaults;
        const std::string fmt = r_o.data_format ? r_o.data_format : "HDF";
        const auto mode = r_o.mode == MIO_XDMF_SERIES_APPEND
                              ? meshioplusplus::XdmfSeriesMode::Append
                              : meshioplusplus::XdmfSeriesMode::Truncate;
        auto* p_out = new mio_xdmf_series{meshioplusplus::XdmfTimeSeriesWriter(
            path, fmt, static_cast<int>(r_o.gzip_level), mode)};
        p_out->mWriter.SetAutoFlush(r_o.auto_flush != 0);
        return p_out;
    });
}

mio_status mio_xdmf_series_write_points_cells(mio_xdmf_series* series, const mio_mesh* mesh) {
    return guarded([&]() -> mio_status {
        if (!series || !mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: series/mesh is NULL");
        series->mWriter.WritePointsCells(mesh->mMesh);
        return MIO_OK;
    });
}

mio_status mio_xdmf_series_write_data(mio_xdmf_series* series, double time, const mio_mesh* mesh) {
    return guarded([&]() -> mio_status {
        if (!series || !mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: series/mesh is NULL");
        series->mWriter.WriteData(time, mesh->mMesh);
        return MIO_OK;
    });
}

mio_status mio_xdmf_series_finalize(mio_xdmf_series* series) {
    return guarded([&]() -> mio_status {
        if (!series)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: series is NULL");
        series->mWriter.Finalize();
        return MIO_OK;
    });
}

mio_status mio_xdmf_series_write_data_arrays(mio_xdmf_series* series, double time,
                                             const mio_named_array* point_arrays,
                                             int64_t num_point_arrays,
                                             const mio_named_array* cell_arrays,
                                             int64_t num_cell_arrays) {
    return guarded([&]() -> mio_status {
        if (!series)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: series is NULL");
        const auto convert = [](const mio_named_array* p_src, std::int64_t n) {
            std::vector<meshioplusplus::XdmfTimeSeriesWriter::NamedArray> out;
            for (std::int64_t i = 0; i < n; ++i) {
                if (!p_src)
                    throw meshioplusplus::WriteError("meshio++: array list is NULL");
                meshioplusplus::XdmfTimeSeriesWriter::NamedArray a;
                a.mName = p_src[i].name ? p_src[i].name : "";
                a.mNumComponents = p_src[i].num_components > 0
                                       ? static_cast<std::size_t>(p_src[i].num_components)
                                       : 1u;
                if (p_src[i].num_values < 0)
                    throw meshioplusplus::WriteError("meshio++: num_values is negative");
                // Copied during the call, per the ABI's "setters copy" rule.
                if (p_src[i].values)
                    a.mValues.assign(p_src[i].values, p_src[i].values + p_src[i].num_values);
                out.push_back(std::move(a));
            }
            return out;
        };
        series->mWriter.WriteData(time, convert(point_arrays, num_point_arrays),
                                  convert(cell_arrays, num_cell_arrays));
        return MIO_OK;
    });
}

mio_status mio_xdmf_series_flush(mio_xdmf_series* series) {
    return guarded([&]() -> mio_status {
        if (!series)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: series is NULL");
        series->mWriter.Flush();
        return MIO_OK;
    });
}

int32_t mio_xdmf_series_finalized(const mio_xdmf_series* series) {
    return guarded_ptr(static_cast<std::int32_t>(-1), [&]() -> std::int32_t {
        if (!series)
            throw meshioplusplus::ReadError("meshio++: series is NULL");
        return series->mWriter.Finalized() ? 1 : 0;
    });
}

int64_t mio_xdmf_series_num_steps(const mio_xdmf_series* series) {
    return guarded_ptr(std::int64_t(-1), [&]() -> std::int64_t {
        if (!series)
            throw meshioplusplus::ReadError("meshio++: series is NULL");
        return static_cast<std::int64_t>(series->mWriter.NumSteps());
    });
}

void mio_xdmf_series_free(mio_xdmf_series* series) {
    // ~XdmfTimeSeriesWriter finalizes and swallows any failure (it must not
    // throw); a caller wanting to see one calls mio_xdmf_series_finalize first.
    delete series;
}

// --------------------------------------------------------------------------
// Settings pipeline. JSON text only across the ABI -- the typed step model
// never crosses it. When the build has no JSON parser the core throws naming
// -DMESHIOPLUSPLUS_WITH_JSON=ON, which guarded() surfaces through
// mio_last_error() as usual.
// --------------------------------------------------------------------------

mio_status mio_pipeline_run_file(const char* settings_path) {
    return guarded([&]() -> mio_status {
        if (!settings_path)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: settings_path is NULL");
        meshioplusplus::run_pipeline_file(settings_path);
        return MIO_OK;
    });
}

int64_t mio_med_mesh_count(const char* path) {
    return guarded_ptr(int64_t{-1}, [&]() -> int64_t {
        if (!path)
            throw std::invalid_argument("MED: path is NULL");
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        return static_cast<int64_t>(meshioplusplus::med_mesh_names(path).size());
#else
        throw std::runtime_error("MED requires HDF5");
#endif
    });
}

int64_t mio_med_mesh_name(const char* path, int64_t index, char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() -> int64_t {
        if (!path)
            throw std::invalid_argument("MED: path is NULL");
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        auto names = meshioplusplus::med_mesh_names(path);
        if (index < 0 || static_cast<std::size_t>(index) >= names.size())
            throw std::invalid_argument("MED: mesh index out of range");
        return copy_string(names[static_cast<std::size_t>(index)], buf, buflen);
#else
        throw std::runtime_error("MED requires HDF5");
#endif
    });
}

mio_mesh* mio_med_read_named(const char* path, const char* name, const mio_read_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!path || !name)
            throw std::invalid_argument("MED: path/name is NULL");
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        meshioplusplus::MedInfo info;
        return new mio_mesh{
            meshioplusplus::read_med_named(path, name, info, capi_read_options(opts))};
#else
        throw std::runtime_error("MED requires HDF5");
#endif
    });
}

mio_status mio_med_write_multi(const char* path, const mio_mesh* const* meshes,
                               const char* const* names, int64_t count, const char* version) {
    return guarded([&]() -> mio_status {
        if (!path || !meshes || !names || count <= 0)
            return fail(MIO_ERR_INVALID_ARG, "MED: provide path, meshes, names and positive count");
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        std::vector<const meshioplusplus::Mesh*> inputs;
        std::vector<meshioplusplus::MedInfo> infos;
        for (int64_t i = 0; i < count; ++i) {
            if (!meshes[i] || !names[i])
                throw std::invalid_argument("MED: NULL mesh/name");
            inputs.push_back(&meshes[i]->mMesh);
            infos.emplace_back();
            infos.back().mMeshName = names[i];
        }
        meshioplusplus::write_med_multi(path, inputs, infos, version ? version : "4.1.0");
        return MIO_OK;
#else
        throw std::runtime_error("MED requires HDF5");
#endif
    });
}

mio_status mio_pipeline_run_json(const char* json_text) {
    return guarded([&]() -> mio_status {
        if (!json_text)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: json_text is NULL");
        meshioplusplus::run_pipeline_json(json_text);
        return MIO_OK;
    });
}

mio_pipeline_report* mio_pipeline_run_file_report(const char* settings_path) {
    return guarded_ptr(static_cast<mio_pipeline_report*>(nullptr), [&]() {
        if (!settings_path)
            throw std::invalid_argument("meshio++: settings_path is NULL");
        return new mio_pipeline_report{
            meshioplusplus::pipeline_report_json(meshioplusplus::run_pipeline_file(settings_path))};
    });
}

mio_pipeline_report* mio_pipeline_run_json_report(const char* json_text) {
    return guarded_ptr(static_cast<mio_pipeline_report*>(nullptr), [&]() {
        if (!json_text)
            throw std::invalid_argument("meshio++: json_text is NULL");
        return new mio_pipeline_report{
            meshioplusplus::pipeline_report_json(meshioplusplus::run_pipeline_json(json_text))};
    });
}

mio_pipeline_report* mio_sequence_pipeline_run_file_report(const char* settings_path) {
    return guarded_ptr(static_cast<mio_pipeline_report*>(nullptr), [&]() {
        if (!settings_path)
            throw std::invalid_argument("meshio++: settings_path is NULL");
        return new mio_pipeline_report{
            meshioplusplus::pipeline_report_json(meshioplusplus::run_sequence_file(settings_path))};
    });
}

mio_pipeline_report* mio_sequence_pipeline_run_json_report(const char* json_text) {
    return guarded_ptr(static_cast<mio_pipeline_report*>(nullptr), [&]() {
        if (!json_text)
            throw std::invalid_argument("meshio++: json_text is NULL");
        return new mio_pipeline_report{
            meshioplusplus::pipeline_report_json(meshioplusplus::run_sequence_json(json_text))};
    });
}

int64_t mio_pipeline_report_json(const mio_pipeline_report* report, char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() {
        if (!report)
            throw std::invalid_argument("meshio++: pipeline report is NULL");
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative report buffer length");
        return copy_string(report->mJson, buf, buflen);
    });
}

void mio_pipeline_report_free(mio_pipeline_report* report) {
    delete report;
}

namespace {

/// mio_sequence_opts -> SequenceInput, shared by the three open* entry points.
meshioplusplus::SequenceInput seq_input_from_opts(const mio_sequence_opts* pOpts) {
    meshioplusplus::SequenceInput in;
    if (!pOpts)
        return in;
    if (pOpts->format)
        in.mFormat = pOpts->format;
    if (pOpts->times && pOpts->num_times > 0)
        in.mTimes.assign(pOpts->times, pOpts->times + pOpts->num_times);
    switch (pOpts->time_from) {
        case 1:
            in.mTimeFrom = meshioplusplus::SequenceTimeFrom::File;
            break;
        case 2:
            in.mTimeFrom = meshioplusplus::SequenceTimeFrom::Filename;
            break;
        case 3:
            in.mTimeFrom = meshioplusplus::SequenceTimeFrom::Index;
            break;
        default:
            in.mTimeFrom = meshioplusplus::SequenceTimeFrom::Auto;
            break;
    }
    in.mSortExplicit = pOpts->sort != 0;
    return in;
}

mio_sequence* seq_open(meshioplusplus::SequenceInput in) {
    // Owned until returned: sequence_expand throws on a pattern that matches
    // nothing, and the guard would otherwise leak the half-built handle.
    auto out = std::make_unique<mio_sequence>();
    out->mEntries = meshioplusplus::sequence_expand(in);
    out->mFormat = in.mFormat;
    out->mOptions = in.mOptions;
    return out.release();
}

const meshioplusplus::SequenceEntry& seq_entry(const mio_sequence* pSeq, int64_t Index) {
    if (!pSeq)
        throw meshioplusplus::ReadError("meshio++: sequence is NULL");
    if (Index < 0 || static_cast<std::size_t>(Index) >= pSeq->mEntries.size())
        throw meshioplusplus::ReadError("meshio++: sequence entry index out of range");
    return pSeq->mEntries[static_cast<std::size_t>(Index)];
}

}  // namespace

void mio_sequence_opts_init(mio_sequence_opts* opts) {
    if (opts)
        *opts = mio_sequence_opts{};
}

mio_sequence* mio_sequence_open(const char* pattern) {
    return mio_sequence_open_ex(pattern, nullptr);
}

mio_sequence* mio_sequence_open_ex(const char* pattern, const mio_sequence_opts* opts) {
    return guarded_ptr(static_cast<mio_sequence*>(nullptr), [&]() -> mio_sequence* {
        if (!pattern)
            throw meshioplusplus::ReadError("meshio++: pattern is NULL");
        meshioplusplus::SequenceInput in = seq_input_from_opts(opts);
        in.mPattern = pattern;
        return seq_open(std::move(in));
    });
}

mio_sequence* mio_sequence_open_list(const char* const* paths, int64_t num_paths,
                                     const mio_sequence_opts* opts) {
    return guarded_ptr(static_cast<mio_sequence*>(nullptr), [&]() -> mio_sequence* {
        if (!paths || num_paths <= 0)
            throw meshioplusplus::ReadError("meshio++: paths is NULL or empty");
        meshioplusplus::SequenceInput in = seq_input_from_opts(opts);
        for (int64_t i = 0; i < num_paths; ++i) {
            if (!paths[i])
                throw meshioplusplus::ReadError("meshio++: a path in the list is NULL");
            in.mPaths.emplace_back(paths[i]);
        }
        return seq_open(std::move(in));
    });
}

int64_t mio_sequence_count(const mio_sequence* seq) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        if (!seq)
            throw meshioplusplus::ReadError("meshio++: sequence is NULL");
        return static_cast<int64_t>(seq->mEntries.size());
    });
}

int64_t mio_sequence_path(const mio_sequence* seq, int64_t index, char* buf, int64_t buflen) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        return copy_string(seq_entry(seq, index).mPath, buf, buflen);
    });
}

int64_t mio_sequence_step(const mio_sequence* seq, int64_t index) {
    return guarded_ptr(static_cast<int64_t>(-1), [&]() -> int64_t {
        return static_cast<int64_t>(seq_entry(seq, index).mStep);
    });
}

mio_status mio_sequence_time(const mio_sequence* seq, int64_t index, double* out_time) {
    return guarded([&]() -> mio_status {
        if (!out_time)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: out_time is NULL");
        *out_time = seq_entry(seq, index).mTime;
        return MIO_OK;
    });
}

int32_t mio_sequence_time_source(const mio_sequence* seq, int64_t index) {
    return guarded_ptr(static_cast<int32_t>(-1), [&]() -> int32_t {
        return static_cast<int32_t>(seq_entry(seq, index).mTimeSource);
    });
}

mio_mesh* mio_sequence_read(mio_sequence* seq, int64_t index) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        seq_entry(seq, index);  // bounds/NULL check with the shared message
        return new mio_mesh{meshioplusplus::sequence_read_step(
            seq->mEntries, static_cast<std::size_t>(index), seq->mFormat, seq->mOptions)};
    });
}

void mio_sequence_free(mio_sequence* seq) {
    delete seq;
}

static_assert(sizeof(mio_agglomerate_opts) == 80,
              "mio_agglomerate_opts grew outside its reserved tail");
static_assert(sizeof(mio_resample_opts) == 80, "mio_resample_opts grew outside its reserved tail");

mio_mesh* mio_blend_steps(const mio_mesh* a, const mio_mesh* b, double w, int32_t blend_points) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!a || !b)
            throw meshioplusplus::ReadError("meshio++: blend_steps: a/b is NULL");
        meshioplusplus::BlendOptions options;
        options.mBlendPoints = blend_points != 0;
        return new mio_mesh{meshioplusplus::blend_steps(a->mMesh, b->mMesh, w, options)};
    });
}

void mio_resample_opts_init(mio_resample_opts* opts) {
    if (opts)
        *opts = mio_resample_opts{};
}

mio_status mio_sequence_resample(const mio_sequence* seq, const char* out_path,
                                 const char* out_format, const mio_resample_opts* opts) {
    return guarded([&]() -> mio_status {
        if (!seq || !out_path || !opts)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: resample: sequence/out_path/opts is NULL");
        if (opts->num_times < 1 || !opts->times)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: resample: no target times");
        if (opts->method < 0 || opts->method > 2 || opts->extrapolate < 0 || opts->extrapolate > 1)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: resample: method/extrapolate out of range");
        meshioplusplus::SequencePipeline p;
        for (const meshioplusplus::SequenceEntry& e : seq->mEntries) {
            p.mInput.mPaths.push_back(e.mPath);
            p.mInput.mTimes.push_back(e.mTime);
        }
        p.mInput.mFormat = seq->mFormat;
        p.mInput.mOptions = seq->mOptions;
        p.mOutput.mPath = out_path;
        if (out_format)
            p.mOutput.mFormat = out_format;
        meshioplusplus::SequenceResample rs;
        rs.mTimes.assign(opts->times, opts->times + opts->num_times);
        rs.mMethod = static_cast<meshioplusplus::ResampleMethod>(opts->method);
        rs.mExtrapolate = static_cast<meshioplusplus::ResampleExtrapolate>(opts->extrapolate);
        rs.mBlendPoints = opts->blend_points != 0;
        p.mResample = rs;
        meshioplusplus::run_sequence_pipeline(p);
        return MIO_OK;
    });
}

mio_status mio_sequence_to_timeseries(const mio_sequence* seq, const char* out_path,
                                      const char* out_format) {
    return mio_sequence_to_timeseries_ex(seq, out_path, out_format, nullptr);
}

mio_status mio_sequence_to_timeseries_ex(const mio_sequence* seq, const char* out_path,
                                         const char* out_format, const mio_write_opts* opts) {
    return guarded([&]() -> mio_status {
        if (!seq || !out_path)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: sequence/out_path is NULL");
        // Rebuild the input from the plan rather than re-globbing: the handle's
        // entries (and their resolved times) are what the caller inspected.
        meshioplusplus::SequenceInput in;
        for (const meshioplusplus::SequenceEntry& e : seq->mEntries)
            in.mPaths.push_back(e.mPath);
        in.mFormat = seq->mFormat;
        in.mOptions = seq->mOptions;
        for (const meshioplusplus::SequenceEntry& e : seq->mEntries)
            in.mTimes.push_back(e.mTime);
        meshioplusplus::SequenceOutput out;
        out.mPath = out_path;
        if (out_format)
            out.mFormat = out_format;
        if (opts) {
            const mio_status opt_status = write_opts_to_cxx(*opts, out.mOptions);
            if (opt_status != MIO_OK)
                return opt_status;
        }
        meshioplusplus::sequence_to_timeseries(in, out);
        return MIO_OK;
    });
}

mio_status mio_timeseries_to_sequence(const char* in_path, const char* in_format,
                                      const char* out_pattern, const char* out_format) {
    return guarded([&]() -> mio_status {
        if (!in_path || !out_pattern)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: in_path/out_pattern is NULL");
        meshioplusplus::SequenceOutput out;
        out.mPath = out_pattern;
        if (out_format)
            out.mFormat = out_format;
        meshioplusplus::timeseries_to_sequence(in_path, in_format ? in_format : "",
                                               meshioplusplus::ReadOptions{}, out);
        return MIO_OK;
    });
}

mio_status mio_sequence_pipeline_run_file(const char* settings_path) {
    return guarded([&]() -> mio_status {
        if (!settings_path)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: settings_path is NULL");
        meshioplusplus::run_sequence_file(settings_path);
        return MIO_OK;
    });
}

mio_status mio_sequence_pipeline_run_json(const char* json_text) {
    return guarded([&]() -> mio_status {
        if (!json_text)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: json_text is NULL");
        meshioplusplus::run_sequence_json(json_text);
        return MIO_OK;
    });
}

/* ------------------------------------------------------------------------- */
/* Regular grids and signed distance                                         */
/* ------------------------------------------------------------------------- */

namespace {

// The C enums duplicate the C++ ones across the ABI, so pin every value: drift
// becomes a compile error rather than a silently wrong argument (the
// mio_region_kind / mio_refine_compare pattern).
static_assert(static_cast<int>(meshioplusplus::VoxelFill::All) == MIO_VOXEL_ALL &&
                  static_cast<int>(meshioplusplus::VoxelFill::Surface) == MIO_VOXEL_SURFACE &&
                  static_cast<int>(meshioplusplus::VoxelFill::Inside) == MIO_VOXEL_INSIDE,
              "mio_voxel_fill drifted from VoxelFill");
static_assert(static_cast<int>(meshioplusplus::SdfSign::Unsigned) == MIO_SDF_UNSIGNED &&
                  static_cast<int>(meshioplusplus::SdfSign::Pseudonormal) == MIO_SDF_PSEUDONORMAL &&
                  static_cast<int>(meshioplusplus::SdfSign::WindingNumber) ==
                      MIO_SDF_WINDING_NUMBER,
              "mio_sdf_sign drifted from SdfSign");
static_assert(static_cast<int>(meshioplusplus::SdfPseudonormalWeight::Angle) ==
                      MIO_SDF_WEIGHT_ANGLE &&
                  static_cast<int>(meshioplusplus::SdfPseudonormalWeight::Area) ==
                      MIO_SDF_WEIGHT_AREA,
              "mio_sdf_weight drifted from SdfPseudonormalWeight");
static_assert(static_cast<int>(meshioplusplus::SdfLocation::Corner) == MIO_SDF_CORNER &&
                  static_cast<int>(meshioplusplus::SdfLocation::Center) == MIO_SDF_CENTER,
              "mio_sdf_location drifted from SdfLocation");
static_assert(static_cast<int>(meshioplusplus::SdfWatertightCheck::Off) == MIO_SDF_WATERTIGHT_OFF &&
                  static_cast<int>(meshioplusplus::SdfWatertightCheck::Warn) ==
                      MIO_SDF_WATERTIGHT_WARN &&
                  static_cast<int>(meshioplusplus::SdfWatertightCheck::Error) ==
                      MIO_SDF_WATERTIGHT_ERROR,
              "mio_sdf_watertight_check drifted from SdfWatertightCheck");

// The sizes are ABI: the Fortran and Julia mirrors of these structs hard-code
// them, and the Julia binding checks them at load. Pinning them here makes a
// field added outside the `reserved` tail a compile error rather than silent
// corruption.
//
// Note mio_voxel_opts carries `sign`/`watertight_check` as its OWN int32 fields
// rather than embedding a mio_sdf_opts. That is deliberate: the C++ VoxelOptions
// does embed SurfaceDistanceOptions by value, and doc/abi.md records that
// growing the inner struct is a Tier A break of the outer one. Not repeating
// that nesting across the ABI keeps the two independently extensible.
static_assert(sizeof(mio_sdf_opts) == 104, "mio_sdf_opts grew outside its reserved tail");
static_assert(sizeof(mio_voxel_opts) == 112, "mio_voxel_opts grew outside its reserved tail");
static_assert(sizeof(mio_gen_opts) == 160, "mio_gen_opts grew outside its reserved tail");
static_assert(offsetof(mio_gen_opts, num_twists) == 128 &&
                  offsetof(mio_gen_opts, complete) == 132 && offsetof(mio_gen_opts, ratio) == 136 &&
                  offsetof(mio_gen_opts, seed) == 144,
              "mio_gen_opts layout changed");
static_assert(sizeof(mio_octree_opts) == 176, "mio_octree_opts grew outside its reserved tail");
static_assert(sizeof(mio_octree_report) == 96, "mio_octree_report layout changed");
static_assert(sizeof(mio_octree_mesh_opts) == 264,
              "mio_octree_mesh_opts grew outside its reserved tail");
static_assert(offsetof(mio_octree_mesh_opts, hex_grading) == 228,
              "mio_octree_mesh_opts layout changed");
static_assert(static_cast<int>(meshioplusplus::OctreeMeshHexGrading::Graded) ==
                  MIO_OCTREE_MESH_HEX_GRADED,
              "mio_octree_mesh_hex_grading drifted from OctreeMeshHexGrading");
static_assert(sizeof(mio_octree_mesh_report) == 392, "mio_octree_mesh_report layout changed");
static_assert(offsetof(mio_octree_mesh_report, num_winding_skipped) == 328,
              "mio_octree_mesh_report layout changed");
static_assert(
    static_cast<int>(meshioplusplus::OctreeMeshOutput::Castellated) ==
            MIO_OCTREE_MESH_CASTELLATED &&
        static_cast<int>(meshioplusplus::OctreeMeshOutput::Tetra) == MIO_OCTREE_MESH_TET &&
        static_cast<int>(meshioplusplus::OctreeMeshOutput::Hexa) == MIO_OCTREE_MESH_HEX &&
        static_cast<int>(meshioplusplus::OctreeMeshConformity::Stuffing) ==
            MIO_OCTREE_MESH_STUFFING &&
        static_cast<int>(meshioplusplus::OctreeMeshConformity::Cleaving) ==
            MIO_OCTREE_MESH_CLEAVING &&
        static_cast<int>(meshioplusplus::OctreeMeshSkinGrading::Graded) == MIO_OCTREE_MESH_GRADED &&
        static_cast<int>(meshioplusplus::OctreeMeshSkinGrading::Uniform) ==
            MIO_OCTREE_MESH_UNIFORM &&
        static_cast<int>(meshioplusplus::OctreeMeshFill::All) == MIO_OCTREE_MESH_FILL_ALL &&
        static_cast<int>(meshioplusplus::OctreeMeshFill::Inside) == MIO_OCTREE_MESH_FILL_INSIDE &&
        static_cast<int>(meshioplusplus::OctreeMeshFill::Cover) == MIO_OCTREE_MESH_FILL_COVER,
    "mio_octree_mesh enums drifted from the core's");
static_assert(static_cast<int>(meshioplusplus::OctreeTransition::Hanging) == MIO_OCTREE_HANGING &&
                  static_cast<int>(meshioplusplus::OctreeTransition::Polyhedra) ==
                      MIO_OCTREE_POLYHEDRA &&
                  static_cast<int>(meshioplusplus::OctreeTransition::Templates) ==
                      MIO_OCTREE_TEMPLATES,
              "mio_octree_transition drifted from OctreeTransition");
static_assert(static_cast<int>(meshioplusplus::detail::PlanarPlane::Auto) == MIO_PLANE_AUTO &&
                  static_cast<int>(meshioplusplus::detail::PlanarPlane::XY) == MIO_PLANE_XY &&
                  static_cast<int>(meshioplusplus::detail::PlanarPlane::XZ) == MIO_PLANE_XZ &&
                  static_cast<int>(meshioplusplus::detail::PlanarPlane::YZ) == MIO_PLANE_YZ,
              "mio_planar_plane drifted from PlanarPlane");
static_assert(static_cast<int>(meshioplusplus::OctreeBalance::None) == MIO_OCTREE_BALANCE_NONE &&
                  static_cast<int>(meshioplusplus::OctreeBalance::Face) ==
                      MIO_OCTREE_BALANCE_FACE &&
                  static_cast<int>(meshioplusplus::OctreeBalance::Full) == MIO_OCTREE_BALANCE_FULL,
              "mio_octree_balance drifted from OctreeBalance");
static_assert(sizeof(mio_compute_sdf_opts) == 232,
              "mio_compute_sdf_opts grew outside its reserved tail");
static_assert(sizeof(mio_surface_quality) == 72,
              "mio_surface_quality grew outside its reserved tail");

/// Translate the flat option struct into the core's SurfaceDistanceOptions.
meshioplusplus::SurfaceDistanceOptions capi_sdf_options(const mio_sdf_opts& rOpts) {
    meshioplusplus::SurfaceDistanceOptions options;
    if (rOpts.sign < 0 || rOpts.sign > MIO_SDF_WINDING_NUMBER)
        throw meshioplusplus::ReadError("meshio++: sdf: unknown sign " +
                                        std::to_string(rOpts.sign));
    options.mSign = static_cast<meshioplusplus::SdfSign>(rOpts.sign);
    if (rOpts.weight < 0 || rOpts.weight > MIO_SDF_WEIGHT_AREA)
        throw meshioplusplus::ReadError("meshio++: sdf: unknown weight " +
                                        std::to_string(rOpts.weight));
    options.mWeight = static_cast<meshioplusplus::SdfPseudonormalWeight>(rOpts.weight);
    if (rOpts.location < 0 || rOpts.location > MIO_SDF_CENTER)
        throw meshioplusplus::ReadError("meshio++: sdf: unknown location " +
                                        std::to_string(rOpts.location));
    options.mLocation = static_cast<meshioplusplus::SdfLocation>(rOpts.location);
    if (rOpts.watertight_check < 0 || rOpts.watertight_check > MIO_SDF_WATERTIGHT_ERROR)
        throw meshioplusplus::ReadError("meshio++: sdf: unknown watertight check " +
                                        std::to_string(rOpts.watertight_check));
    options.mWatertightCheck =
        static_cast<meshioplusplus::SdfWatertightCheck>(rOpts.watertight_check);
    options.mBand = rOpts.band;
    options.mRecordClosestCell = rOpts.record_closest_cell != 0;
    options.mRecordInside = rOpts.record_inside != 0;
    if (rOpts.surface_region != nullptr)
        options.mSurfaceRegion = rOpts.surface_region;
    options.mGridCellSize = rOpts.grid_cell_size;
    options.mMaxWindingWork = rOpts.max_winding_work;
    return options;
}

/// Translate the flat option struct into the core's VoxelOptions.
meshioplusplus::VoxelOptions capi_voxel_options(const mio_voxel_opts& rOpts) {
    meshioplusplus::VoxelOptions options;
    if (rOpts.resolution != nullptr)
        options.mResolution = std::array<std::int64_t, 3>{
            {rOpts.resolution[0], rOpts.resolution[1], rOpts.resolution[2]}};
    if (rOpts.cell_size > 0.0)
        options.mCellSize = rOpts.cell_size;
    if (rOpts.bounds != nullptr)
        options.mBounds =
            std::array<double, 6>{{rOpts.bounds[0], rOpts.bounds[1], rOpts.bounds[2],
                                   rOpts.bounds[3], rOpts.bounds[4], rOpts.bounds[5]}};
    options.mPadding = rOpts.padding;
    options.mPaddingRelative = rOpts.padding_relative;
    if (rOpts.fill < 0 || rOpts.fill > MIO_VOXEL_INSIDE)
        throw meshioplusplus::ReadError("meshio++: voxelize: unknown fill " +
                                        std::to_string(rOpts.fill));
    options.mFill = static_cast<meshioplusplus::VoxelFill>(rOpts.fill);
    options.mAttachOccupancy = rOpts.attach_occupancy != 0;
    options.mMaxCells = rOpts.max_cells;
    if (rOpts.sign < 0 || rOpts.sign > MIO_SDF_WINDING_NUMBER)
        throw meshioplusplus::ReadError("meshio++: voxelize: unknown sign " +
                                        std::to_string(rOpts.sign));
    options.mDistance.mSign = static_cast<meshioplusplus::SdfSign>(rOpts.sign);
    if (rOpts.watertight_check < 0 || rOpts.watertight_check > MIO_SDF_WATERTIGHT_ERROR)
        throw meshioplusplus::ReadError("meshio++: voxelize: unknown watertight check " +
                                        std::to_string(rOpts.watertight_check));
    options.mDistance.mWatertightCheck =
        static_cast<meshioplusplus::SdfWatertightCheck>(rOpts.watertight_check);
    if (rOpts.plane < 0 || rOpts.plane > MIO_PLANE_YZ)
        throw meshioplusplus::ReadError("meshio++: voxelize: unknown plane " +
                                        std::to_string(rOpts.plane));
    options.mPlane = static_cast<meshioplusplus::detail::PlanarPlane>(rOpts.plane);
    return options;
}

/// Translate the flat option struct into the core's OctreeOptions.
meshioplusplus::OctreeOptions capi_octree_options(const mio_octree_opts& rOpts) {
    meshioplusplus::OctreeOptions options;
    if (rOpts.root_cells != nullptr)
        options.mRootCells = std::array<std::int64_t, 3>{
            {rOpts.root_cells[0], rOpts.root_cells[1], rOpts.root_cells[2]}};
    if (rOpts.bounds != nullptr)
        options.mBounds =
            std::array<double, 6>{{rOpts.bounds[0], rOpts.bounds[1], rOpts.bounds[2],
                                   rOpts.bounds[3], rOpts.bounds[4], rOpts.bounds[5]}};
    options.mPadding = rOpts.padding;
    options.mPaddingRelative = rOpts.padding_relative;
    options.mSampleDistance = rOpts.sample_distance;
    options.mMinLevel = rOpts.min_level;
    options.mMaxLevel = rOpts.max_level;
    options.mDimension = rOpts.dimension;
    options.mMaxCells = rOpts.max_cells;
    if (rOpts.balance < 0 || rOpts.balance > MIO_OCTREE_BALANCE_FULL)
        throw meshioplusplus::ReadError("meshio++: octree: unknown balance " +
                                        std::to_string(rOpts.balance));
    options.mBalance = static_cast<meshioplusplus::OctreeBalance>(rOpts.balance);
    if (rOpts.fill < 0 || rOpts.fill > MIO_VOXEL_INSIDE)
        throw meshioplusplus::ReadError("meshio++: octree: unknown fill " +
                                        std::to_string(rOpts.fill));
    options.mFill = static_cast<meshioplusplus::VoxelFill>(rOpts.fill);
    options.mAttachDistance = rOpts.attach_distance != 0;
    options.mEvenOdd = rOpts.even_odd != 0;
    if (rOpts.sign < 0 || rOpts.sign > MIO_SDF_WINDING_NUMBER)
        throw meshioplusplus::ReadError("meshio++: octree: unknown sign " +
                                        std::to_string(rOpts.sign));
    options.mDistance.mSign = static_cast<meshioplusplus::SdfSign>(rOpts.sign);
    if (rOpts.watertight_check < 0 || rOpts.watertight_check > MIO_SDF_WATERTIGHT_ERROR)
        throw meshioplusplus::ReadError("meshio++: octree: unknown watertight check " +
                                        std::to_string(rOpts.watertight_check));
    options.mDistance.mWatertightCheck =
        static_cast<meshioplusplus::SdfWatertightCheck>(rOpts.watertight_check);
    if (rOpts.transition < 0 || rOpts.transition > MIO_OCTREE_TEMPLATES)
        throw meshioplusplus::ReadError("meshio++: octree: unknown transition " +
                                        std::to_string(rOpts.transition));
    options.mTransition = static_cast<meshioplusplus::OctreeTransition>(rOpts.transition);
    options.mSnap = rOpts.snap != 0;
    options.mSurfaceLevel = rOpts.surface_level;
    options.mCurvatureThreshold = rOpts.curvature_threshold;
    options.mFeatureAngle = rOpts.feature_angle;
    if (rOpts.size_array != nullptr)
        options.mSizeArray = rOpts.size_array;
    if (rOpts.plane < 0 || rOpts.plane > MIO_PLANE_YZ)
        throw meshioplusplus::ReadError("meshio++: octree: unknown plane " +
                                        std::to_string(rOpts.plane));
    options.mPlane = static_cast<meshioplusplus::detail::PlanarPlane>(rOpts.plane);
    return options;
}

/// Translate the flat option struct into the core's OctreeMeshOptions.
meshioplusplus::OctreeMeshOptions capi_octree_mesh_options(const mio_octree_mesh_opts& rOpts) {
    meshioplusplus::OctreeMeshOptions options;
    auto check = [](std::int32_t Value, std::int32_t Max, const char* pWhat) {
        if (Value < 0 || Value > Max)
            throw meshioplusplus::ReadError(std::string("meshio++: octree_mesh: unknown ") + pWhat +
                                            " " + std::to_string(Value));
        return Value;
    };
    options.mOutput = static_cast<meshioplusplus::OctreeMeshOutput>(
        check(rOpts.output, MIO_OCTREE_MESH_HEX, "output"));
    options.mConformity = static_cast<meshioplusplus::OctreeMeshConformity>(
        check(rOpts.conformity, MIO_OCTREE_MESH_CLEAVING, "conformity"));
    options.mSkinGrading = static_cast<meshioplusplus::OctreeMeshSkinGrading>(
        check(rOpts.skin_grading, MIO_OCTREE_MESH_UNIFORM, "skin grading"));
    options.mFill = static_cast<meshioplusplus::OctreeMeshFill>(
        check(rOpts.fill, MIO_OCTREE_MESH_FILL_COVER, "fill"));
    options.mHexGrading = static_cast<meshioplusplus::OctreeMeshHexGrading>(
        check(rOpts.hex_grading, MIO_OCTREE_MESH_HEX_GRADED, "hex grading"));
    if (rOpts.bounds != nullptr)
        options.mBounds =
            std::array<double, 6>{{rOpts.bounds[0], rOpts.bounds[1], rOpts.bounds[2],
                                   rOpts.bounds[3], rOpts.bounds[4], rOpts.bounds[5]}};
    if (rOpts.size_array != nullptr)
        options.mSizeArray = rOpts.size_array;
    if (rOpts.surface_region != nullptr)
        options.mSurfaceRegion = rOpts.surface_region;
    options.mMinSize = rOpts.min_size;
    options.mMaxSize = rOpts.max_size;
    options.mPadding = rOpts.padding;
    options.mPaddingRelative = rOpts.padding_relative;
    options.mCurvatureAngle = rOpts.curvature_angle;
    options.mFeatureAngle = rOpts.feature_angle;
    options.mCornerAngle = rOpts.corner_angle;
    options.mFeatureSize = rOpts.feature_size;
    options.mGapCells = rOpts.gap_cells;
    options.mSurfaceSize = rOpts.surface_size;
    options.mWindingThreshold = rOpts.winding_threshold;
    options.mWindingBeta = rOpts.winding_beta;
    options.mAlphaLong = rOpts.alpha_long;
    options.mAlphaShort = rOpts.alpha_short;
    options.mMaxLevel = rOpts.max_level;
    options.mGradingLayers = rOpts.grading_layers;
    options.mFitIterations = rOpts.fit_iterations;
    options.mUntangleIterations = rOpts.untangle_iterations;
    options.mMaxCells = rOpts.max_cells;
    options.mMaxTets = rOpts.max_tets;
    options.mGapOutside = rOpts.gap_outside != 0;
    options.mPreserveFeatures = rOpts.preserve_features != 0;
    options.mDihedralFloorGuard = rOpts.dihedral_floor_guard != 0;
    options.mKeepOutside = rOpts.keep_outside != 0;
    options.mAttachDistance = rOpts.attach_distance != 0;
    options.mAttachQuality = rOpts.attach_quality != 0;
    options.mPillow = rOpts.pillow != 0;
    return options;
}

/// The warnings of this thread's last mio_octree_mesh call.
std::string& capi_octree_mesh_warnings() {
    thread_local std::string warnings;
    return warnings;
}

/// Translate the flat option struct into the core's SdfOptions.
meshioplusplus::SdfOptions capi_compute_sdf_options(const mio_compute_sdf_opts& rOpts) {
    meshioplusplus::SdfOptions options;
    if (rOpts.structure < 0 || rOpts.structure > MIO_SDF_OCTREE)
        throw meshioplusplus::ReadError("meshio++: sdf: unknown structure " +
                                        std::to_string(rOpts.structure));
    options.mStructure = static_cast<meshioplusplus::SdfStructure>(rOpts.structure);
    if (rOpts.resolution != nullptr)
        options.mResolution = std::array<std::int64_t, 3>{
            {rOpts.resolution[0], rOpts.resolution[1], rOpts.resolution[2]}};
    if (rOpts.cell_size > 0.0)
        options.mCellSize = rOpts.cell_size;
    if (rOpts.bounds != nullptr)
        options.mBounds =
            std::array<double, 6>{{rOpts.bounds[0], rOpts.bounds[1], rOpts.bounds[2],
                                   rOpts.bounds[3], rOpts.bounds[4], rOpts.bounds[5]}};
    options.mPadding = rOpts.padding;
    options.mPaddingRelative = rOpts.padding_relative;
    options.mRootResolution = rOpts.root_resolution;
    options.mMaxDepth = rOpts.max_depth;
    options.mBandCells = rOpts.band_cells;
    options.mRecordLevels = rOpts.record_levels != 0;
    options.mMaxCells = rOpts.max_cells;
    options.mDistance = capi_sdf_options(rOpts.distance);
    return options;
}

/// Translate the flat option struct into the core's RemeshVolumeOptions.
meshioplusplus::RemeshVolumeOptions capi_remesh_volume_options(
    const mio_remesh_volume_opts& rOpts) {
    meshioplusplus::RemeshVolumeOptions options;
    if (rOpts.resolution != nullptr)
        options.mResolution = std::array<std::int64_t, 3>{
            {rOpts.resolution[0], rOpts.resolution[1], rOpts.resolution[2]}};
    if (rOpts.cell_size > 0.0)
        options.mCellSize = rOpts.cell_size;
    if (rOpts.bounds != nullptr)
        options.mBounds =
            std::array<double, 6>{{rOpts.bounds[0], rOpts.bounds[1], rOpts.bounds[2],
                                   rOpts.bounds[3], rOpts.bounds[4], rOpts.bounds[5]}};
    options.mPadding = rOpts.padding;
    options.mPaddingRelative = rOpts.padding_relative;
    // Verbatim, not gated on > 0: mMaxCells/mMaxTets <= 0 means "unlimited"
    // to the C++ checks themselves (sdf.cpp/remesh_volume.cpp's own
    // `> 0 && ...` guards), so passing it through unconditionally is what
    // makes explicitly requesting "unlimited" (0) actually work, rather than
    // silently falling back to the struct's C++-side default of 20000000.
    options.mMaxCells = rOpts.max_cells;
    options.mMaxTets = rOpts.max_tets;
    options.mWarpFraction = rOpts.warp_fraction;
    options.mDistance = capi_sdf_options(rOpts.distance);
    return options;
}

void capi_fill_quality(const meshioplusplus::SurfaceQuality& rQuality, mio_surface_quality* pOut) {
    if (!pOut)
        return;
    *pOut = mio_surface_quality{};
    pOut->boundary_edges = rQuality.mBoundaryEdges;
    pOut->non_manifold_edges = rQuality.mNonManifoldEdges;
    pOut->inconsistent_pairs = rQuality.mInconsistentPairs;
    pOut->degenerate_triangles = rQuality.mDegenerateTriangles;
    pOut->watertight = rQuality.mWatertight ? 1 : 0;
}

}  // namespace

void mio_sdf_opts_init(mio_sdf_opts* opts) {
    if (!opts)
        return;
    *opts = mio_sdf_opts{};  // value-initialized: all zero == unsigned, full field
    opts->sign = MIO_SDF_PSEUDONORMAL;
    opts->watertight_check = MIO_SDF_WATERTIGHT_WARN;
    opts->max_winding_work = 2.0e9;
}

void mio_voxel_opts_init(mio_voxel_opts* opts) {
    if (!opts)
        return;
    *opts = mio_voxel_opts{};
    opts->max_cells = 20000000;
    opts->sign = MIO_SDF_PSEUDONORMAL;
    opts->watertight_check = MIO_SDF_WATERTIGHT_WARN;
}

void mio_gen_opts_init(mio_gen_opts* opts) {
    if (!opts)
        return;
    *opts = mio_gen_opts{};
    opts->counts[0] = opts->counts[1] = opts->counts[2] = 8;
    opts->size[0] = opts->size[1] = opts->size[2] = 1.0;
    opts->radius = 1.0;
    opts->height = 1.0;
    opts->radius2 = 0.5;
    opts->angle = 6.283185307179586;
    opts->max_cells = 20000000;
    opts->order = 1;
    opts->num_twists = 1;
    opts->ratio = 1.0;
}

void mio_octree_opts_init(mio_octree_opts* opts) {
    if (!opts)
        return;
    *opts = mio_octree_opts{};
    opts->max_level = 4;
    opts->max_cells = 20000000;
    opts->balance = MIO_OCTREE_BALANCE_FACE;
    opts->surface_level = -1;
    opts->sign = MIO_SDF_PSEUDONORMAL;
    opts->watertight_check = MIO_SDF_WATERTIGHT_WARN;
}

void mio_octree_mesh_opts_init(mio_octree_mesh_opts* opts) {
    if (!opts)
        return;
    *opts = mio_octree_mesh_opts{};
    const meshioplusplus::OctreeMeshOptions d;
    opts->padding_relative = d.mPaddingRelative;
    opts->curvature_angle = d.mCurvatureAngle;
    opts->feature_angle = d.mFeatureAngle;
    opts->corner_angle = d.mCornerAngle;
    opts->gap_cells = d.mGapCells;
    opts->winding_threshold = d.mWindingThreshold;
    opts->winding_beta = d.mWindingBeta;
    opts->alpha_long = d.mAlphaLong;
    opts->alpha_short = d.mAlphaShort;
    opts->max_level = d.mMaxLevel;
    opts->grading_layers = d.mGradingLayers;
    opts->fit_iterations = d.mFitIterations;
    opts->untangle_iterations = d.mUntangleIterations;
    opts->max_cells = d.mMaxCells;
    opts->max_tets = d.mMaxTets;
    opts->output = static_cast<std::int32_t>(d.mOutput);
    opts->conformity = static_cast<std::int32_t>(d.mConformity);
    opts->skin_grading = static_cast<std::int32_t>(d.mSkinGrading);
    opts->fill = static_cast<std::int32_t>(d.mFill);
    opts->gap_outside = d.mGapOutside ? 1 : 0;
    opts->preserve_features = d.mPreserveFeatures ? 1 : 0;
    opts->dihedral_floor_guard = d.mDihedralFloorGuard ? 1 : 0;
    opts->keep_outside = d.mKeepOutside ? 1 : 0;
    opts->attach_distance = d.mAttachDistance ? 1 : 0;
    opts->attach_quality = d.mAttachQuality ? 1 : 0;
    opts->pillow = d.mPillow ? 1 : 0;
    opts->hex_grading = static_cast<std::int32_t>(d.mHexGrading);
}

void mio_compute_sdf_opts_init(mio_compute_sdf_opts* opts) {
    if (!opts)
        return;
    *opts = mio_compute_sdf_opts{};
    opts->padding_relative = 0.1;
    opts->band_cells = 1.0;
    opts->max_cells = 20000000;
    opts->root_resolution = 8;
    opts->max_depth = 4;
    opts->structure = MIO_SDF_VOXEL;
    opts->record_levels = 1;
    mio_sdf_opts_init(&opts->distance);
}

mio_mesh* mio_compute_sdf(const mio_mesh* surface, const mio_compute_sdf_opts* opts,
                          int64_t dims_out[3], double origin_out[3], double spacing_out[3],
                          int64_t* max_depth_out, int64_t* num_banded,
                          mio_surface_quality* quality) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!surface)
            throw meshioplusplus::ReadError("meshio++: surface is NULL");
        if (!opts)
            throw meshioplusplus::ReadError(
                "meshio++: sdf: options are NULL, but a grid sizing must be given");
        meshioplusplus::SdfResult r =
            meshioplusplus::compute_sdf(surface->mMesh, capi_compute_sdf_options(*opts));
        for (int k = 0; k < 3; ++k) {
            if (dims_out)
                dims_out[k] = r.mDims[static_cast<std::size_t>(k)];
            if (origin_out)
                origin_out[k] = r.mOrigin[static_cast<std::size_t>(k)];
            if (spacing_out)
                spacing_out[k] = r.mSpacing[static_cast<std::size_t>(k)];
        }
        if (max_depth_out)
            *max_depth_out = r.mMaxDepth;
        if (num_banded)
            *num_banded = r.mNumBanded;
        capi_fill_quality(r.mQuality, quality);
        return new mio_mesh{std::move(r.mMesh)};
    });
}

static_assert(sizeof(mio_remesh_volume_opts) == 216,
              "mio_remesh_volume_opts grew outside its reserved tail");
static_assert(sizeof(mio_remesh_volume_report) == 136,
              "mio_remesh_volume_report grew outside its reserved tail");

void mio_remesh_volume_opts_init(mio_remesh_volume_opts* opts) {
    if (!opts)
        return;
    *opts = mio_remesh_volume_opts{};
    opts->padding_relative = 0.1;
    opts->max_cells = 20000000;
    opts->max_tets = 20000000;
    opts->warp_fraction = 0.35;
    mio_sdf_opts_init(&opts->distance);
}

mio_mesh* mio_remesh_volume_ex(const mio_mesh* mesh, const mio_remesh_volume_opts* opts,
                               mio_remesh_volume_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        if (!opts)
            throw meshioplusplus::ReadError(
                "meshio++: remesh_volume: options are NULL, but exactly one of resolution and "
                "cell_size must be given");
        meshioplusplus::RemeshVolumeResult r =
            meshioplusplus::remesh_volume(mesh->mMesh, capi_remesh_volume_options(*opts));
        if (report) {
            *report = mio_remesh_volume_report{};
            capi_fill_quality(r.mQuality, &report->input_quality);
            report->num_tets = r.mNumTets;
            report->num_vertices_warped = r.mNumVerticesWarped;
            report->num_tets_rejected = r.mNumTetsRejected;
            report->num_non_manifold_edges = r.mNumNonManifoldEdges;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

mio_mesh* mio_grid(const int64_t dims[3], const double origin[3], const double spacing[3],
                   int64_t max_cells) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!dims)
            throw meshioplusplus::ReadError("meshio++: grid: dims is NULL");
        const std::array<std::int64_t, 3> d{{dims[0], dims[1], dims[2]}};
        const std::array<double, 3> o =
            origin ? std::array<double, 3>{{origin[0], origin[1], origin[2]}}
                   : std::array<double, 3>{{0.0, 0.0, 0.0}};
        const std::array<double, 3> s =
            spacing ? std::array<double, 3>{{spacing[0], spacing[1], spacing[2]}}
                    : std::array<double, 3>{{1.0, 1.0, 1.0}};
        return new mio_mesh{meshioplusplus::grid(d, o, s, max_cells)};
    });
}

namespace {

// The spelling of a generator's cell type, or its default when the options leave
// it unset. A value the shape does not build is an error by name, not a silent
// default.
std::string capi_gen_cell(const mio_gen_opts& rOpts, const char* pShape, const char* pDefault,
                          std::initializer_list<std::pair<int, const char*>> Allowed) {
    if (rOpts.cell_type == MIO_GEN_DEFAULT)
        return pDefault;
    for (const auto& [value, name] : Allowed)
        if (rOpts.cell_type == value)
            return name;
    throw meshioplusplus::ReadError(std::string("meshio++: ") + pShape + ": cell_type " +
                                    std::to_string(rOpts.cell_type) +
                                    " is not built by this shape");
}

mio_gen_opts capi_gen_defaults(const mio_gen_opts* pOpts) {
    mio_gen_opts out;
    if (pOpts) {
        out = *pOpts;
        if (out.pad != 0)
            throw meshioplusplus::ReadError("meshio++: mio_gen_opts: pad must be zero");
        for (int64_t r : out.reserved)
            if (r != 0)
                throw meshioplusplus::ReadError(
                    "meshio++: mio_gen_opts: reserved fields must be zero");
    } else {
        mio_gen_opts_init(&out);
    }
    if (out.order == 0)
        out.order = 1;
    return out;
}

}  // namespace

mio_mesh* mio_rectangle(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        const std::string type = capi_gen_cell(
            o, "rectangle", "quad", {{MIO_GEN_QUAD, "quad"}, {MIO_GEN_TRIANGLE, "triangle"}});
        static const char* const splits[] = {"zigzag", "up", "down", "center"};
        if (o.split < 0 || o.split > 3)
            throw meshioplusplus::ReadError("meshio++: rectangle: split is not a mio_gen_split");
        return new mio_mesh{
            meshioplusplus::rectangle({{o.counts[0], o.counts[1]}}, {{o.size[0], o.size[1]}},
                                      {{o.origin[0], o.origin[1], o.origin[2]}}, type,
                                      splits[o.split], o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_box(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        const std::string type = capi_gen_cell(o, "box", "hexahedron",
                                               {{MIO_GEN_HEXAHEDRON, "hexahedron"},
                                                {MIO_GEN_TETRA, "tetra"},
                                                {MIO_GEN_WEDGE, "wedge"}});
        return new mio_mesh{meshioplusplus::box({{o.counts[0], o.counts[1], o.counts[2]}},
                                                {{o.size[0], o.size[1], o.size[2]}},
                                                {{o.origin[0], o.origin[1], o.origin[2]}}, type,
                                                o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_disk(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        const std::string type = capi_gen_cell(
            o, "disk", "quad", {{MIO_GEN_QUAD, "quad"}, {MIO_GEN_TRIANGLE, "triangle"}});
        return new mio_mesh{meshioplusplus::disk(o.radius, o.counts[0],
                                                 {{o.origin[0], o.origin[1], o.origin[2]}}, type,
                                                 o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_uv_sphere(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        const std::string type = capi_gen_cell(
            o, "uv_sphere", "triangle", {{MIO_GEN_TRIANGLE, "triangle"}, {MIO_GEN_QUAD, "quad"}});
        return new mio_mesh{meshioplusplus::uv_sphere(o.radius, o.counts[0], o.counts[1],
                                                      {{o.origin[0], o.origin[1], o.origin[2]}},
                                                      type, o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_icosa_sphere(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        if (o.cell_type != MIO_GEN_DEFAULT && o.cell_type != MIO_GEN_TRIANGLE)
            throw meshioplusplus::ReadError(
                "meshio++: icosa_sphere: cell_type is not built by this shape");
        return new mio_mesh{meshioplusplus::icosa_sphere(o.radius, o.counts[0],
                                                         {{o.origin[0], o.origin[1], o.origin[2]}},
                                                         o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_cylinder(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        const std::string type =
            capi_gen_cell(o, "cylinder", "hexahedron",
                          {{MIO_GEN_HEXAHEDRON, "hexahedron"}, {MIO_GEN_WEDGE, "wedge"}});
        return new mio_mesh{meshioplusplus::cylinder(o.radius, o.height, o.counts[0], o.counts[1],
                                                     {{o.origin[0], o.origin[1], o.origin[2]}},
                                                     type, o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_reference_triangle(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        if (o.cell_type != MIO_GEN_DEFAULT && o.cell_type != MIO_GEN_TRIANGLE)
            throw meshioplusplus::ReadError(
                "meshio++: reference_triangle: cell_type is not built by this shape");
        return new mio_mesh{meshioplusplus::reference_triangle(
            o.counts[0], {{o.size[0], o.size[1]}}, {{o.origin[0], o.origin[1], o.origin[2]}},
            o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_ngon(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        if (o.cell_type != MIO_GEN_DEFAULT && o.cell_type != MIO_GEN_TRIANGLE)
            throw meshioplusplus::ReadError("meshio++: ngon: cell_type is not built by this shape");
        return new mio_mesh{meshioplusplus::ngon(o.counts[0], o.radius, o.counts[1],
                                                 {{o.origin[0], o.origin[1], o.origin[2]}}, o.order,
                                                 o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_annulus(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        const std::string type = capi_gen_cell(
            o, "annulus", "quad", {{MIO_GEN_QUAD, "quad"}, {MIO_GEN_TRIANGLE, "triangle"}});
        return new mio_mesh{meshioplusplus::annulus(
            o.radius2 == 0.0 ? 0.5 : o.radius2, o.radius, o.counts[0], o.counts[1],
            o.angle == 0.0 ? 6.283185307179586 : o.angle, {{o.origin[0], o.origin[1], o.origin[2]}},
            type, o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_lshape(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        const std::string type = capi_gen_cell(
            o, "lshape", "quad", {{MIO_GEN_QUAD, "quad"}, {MIO_GEN_TRIANGLE, "triangle"}});
        return new mio_mesh{meshioplusplus::lshape(o.counts[0], o.size[0],
                                                   {{o.origin[0], o.origin[1], o.origin[2]}}, type,
                                                   o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_plate_with_hole(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        const std::string type = capi_gen_cell(
            o, "plate_with_hole", "quad", {{MIO_GEN_QUAD, "quad"}, {MIO_GEN_TRIANGLE, "triangle"}});
        return new mio_mesh{
            meshioplusplus::plate_with_hole(o.radius, o.size[0], o.counts[0], o.counts[1],
                                            {{o.origin[0], o.origin[1], o.origin[2]}}, type,
                                            o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_honeycomb(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        if (o.cell_type != MIO_GEN_DEFAULT)
            throw meshioplusplus::ReadError(
                "meshio++: honeycomb: cell_type is not built by this shape");
        if (o.order != 1)
            throw meshioplusplus::ReadError(
                "meshio++: honeycomb: order must be 1, a hexagon has no quadratic form");
        return new mio_mesh{meshioplusplus::honeycomb({{o.counts[0], o.counts[1]}}, o.size[0],
                                                      {{o.origin[0], o.origin[1], o.origin[2]}},
                                                      o.max_cells)};
    });
}

namespace {

void capi_gen_linear_only(const mio_gen_opts& rOpts, const char* pShape) {
    if (rOpts.cell_type != MIO_GEN_DEFAULT)
        throw meshioplusplus::ReadError(std::string("meshio++: ") + pShape +
                                        ": cell_type is not built by this shape");
    if (rOpts.order != 1)
        throw meshioplusplus::ReadError(std::string("meshio++: ") + pShape +
                                        ": order must be 1, this shape has no higher order");
}

void capi_gen_no_cell(const mio_gen_opts& rOpts, const char* pShape) {
    if (rOpts.cell_type != MIO_GEN_DEFAULT)
        throw meshioplusplus::ReadError(std::string("meshio++: ") + pShape +
                                        ": cell_type is not built by this shape");
}

std::array<double, 3> capi_gen_origin(const mio_gen_opts& rOpts) {
    return {{rOpts.origin[0], rOpts.origin[1], rOpts.origin[2]}};
}

double capi_gen_or(double Value, double Default) {
    return Value == 0.0 ? Default : Value;
}

std::string capi_gen_surface_cell(const mio_gen_opts& rOpts, const char* pShape) {
    return capi_gen_cell(rOpts, pShape, "quad",
                         {{MIO_GEN_QUAD, "quad"}, {MIO_GEN_TRIANGLE, "triangle"}});
}

std::string capi_gen_solid_cell(const mio_gen_opts& rOpts, const char* pShape) {
    return capi_gen_cell(rOpts, pShape, "hexahedron",
                         {{MIO_GEN_HEXAHEDRON, "hexahedron"}, {MIO_GEN_TETRA, "tetra"}});
}

}  // namespace

mio_mesh* mio_line(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        capi_gen_no_cell(o, "line");
        return new mio_mesh{meshioplusplus::line(o.counts[0], capi_gen_origin(o),
                                                 {{o.size[0], o.size[1], o.size[2]}}, o.order,
                                                 o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_circle(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        capi_gen_no_cell(o, "circle");
        return new mio_mesh{meshioplusplus::circle(o.radius, o.counts[0], capi_gen_origin(o),
                                                   o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_helix(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        capi_gen_no_cell(o, "helix");
        const double turns =
            capi_gen_or(o.angle, 4.0 * 3.141592653589793) / (2.0 * 3.141592653589793);
        return new mio_mesh{meshioplusplus::helix(o.radius, o.height, turns, o.counts[0],
                                                  capi_gen_origin(o), o.order, o.max_cells,
                                                  o.complete != 0)};
    });
}

mio_mesh* mio_sphere_points(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        capi_gen_linear_only(o, "sphere_points");
        return new mio_mesh{
            meshioplusplus::sphere_points(o.counts[0], o.radius, capi_gen_origin(o), o.max_cells)};
    });
}

mio_mesh* mio_box_points(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        capi_gen_linear_only(o, "box_points");
        return new mio_mesh{meshioplusplus::box_points(o.counts[0],
                                                       {{o.size[0], o.size[1], o.size[2]}},
                                                       capi_gen_origin(o), o.seed, o.max_cells)};
    });
}

mio_mesh* mio_tetra_sphere(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        if (o.cell_type != MIO_GEN_DEFAULT && o.cell_type != MIO_GEN_TRIANGLE)
            capi_gen_no_cell(o, "tetra_sphere");
        return new mio_mesh{meshioplusplus::tetra_sphere(o.radius, o.counts[0], capi_gen_origin(o),
                                                         o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_octa_sphere(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        if (o.cell_type != MIO_GEN_DEFAULT && o.cell_type != MIO_GEN_TRIANGLE)
            capi_gen_no_cell(o, "octa_sphere");
        return new mio_mesh{meshioplusplus::octa_sphere(o.radius, o.counts[0], capi_gen_origin(o),
                                                        o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_cube_sphere(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        return new mio_mesh{meshioplusplus::cube_sphere(o.radius, o.counts[0], capi_gen_origin(o),
                                                        capi_gen_surface_cell(o, "cube_sphere"),
                                                        o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_geo_sphere(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        if (o.cell_type != MIO_GEN_DEFAULT && o.cell_type != MIO_GEN_TRIANGLE)
            capi_gen_no_cell(o, "geo_sphere");
        return new mio_mesh{meshioplusplus::geo_sphere(o.radius, o.counts[0], capi_gen_origin(o),
                                                       o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_torus(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        return new mio_mesh{meshioplusplus::torus(
            o.radius, capi_gen_or(o.radius2, 0.25), o.counts[0], o.counts[1], capi_gen_origin(o),
            capi_gen_surface_cell(o, "torus"), o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_moebius(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        return new mio_mesh{meshioplusplus::moebius(
            o.radius, capi_gen_or(o.radius2, 0.5), o.counts[0], o.counts[1], o.num_twists,
            capi_gen_origin(o), capi_gen_surface_cell(o, "moebius"), o.order, o.max_cells,
            o.complete != 0)};
    });
}

mio_mesh* mio_tube(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        return new mio_mesh{meshioplusplus::tube(
            o.radius, o.height, o.counts[0], o.counts[1], capi_gen_origin(o),
            capi_gen_surface_cell(o, "tube"), o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_cone(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        return new mio_mesh{meshioplusplus::cone(
            o.radius, o.height, o.counts[0], o.counts[1], capi_gen_origin(o),
            capi_gen_surface_cell(o, "cone"), o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_capsule(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        return new mio_mesh{meshioplusplus::capsule(
            o.radius, o.height, o.counts[0], o.counts[1], o.counts[2], capi_gen_origin(o),
            capi_gen_surface_cell(o, "capsule"), o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_box_surface(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        return new mio_mesh{meshioplusplus::box_surface(
            {{o.counts[0], o.counts[1], o.counts[2]}}, {{o.size[0], o.size[1], o.size[2]}},
            capi_gen_origin(o), capi_gen_surface_cell(o, "box_surface"), o.order, o.max_cells,
            o.complete != 0)};
    });
}

mio_mesh* mio_platonic(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        capi_gen_linear_only(o, "platonic");
        const char* name = nullptr;
        switch (o.counts[0]) {
            case 4:
                name = "tetrahedron";
                break;
            case 6:
                name = "cube";
                break;
            case 8:
                name = "octahedron";
                break;
            case 12:
                name = "dodecahedron";
                break;
            case 20:
                name = "icosahedron";
                break;
            default:
                throw meshioplusplus::ReadError(
                    "meshio++: platonic: counts[0] must be a face count: 4, 6, 8, 12 or 20, got " +
                    std::to_string(o.counts[0]));
        }
        return new mio_mesh{
            meshioplusplus::platonic(name, o.radius, capi_gen_origin(o), o.max_cells)};
    });
}

mio_mesh* mio_ball(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        return new mio_mesh{meshioplusplus::ball(o.radius, o.counts[0], capi_gen_origin(o),
                                                 capi_gen_solid_cell(o, "ball"), o.order,
                                                 o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_spherical_shell(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        if (o.cell_type != MIO_GEN_DEFAULT && o.cell_type != MIO_GEN_HEXAHEDRON)
            capi_gen_no_cell(o, "spherical_shell");
        return new mio_mesh{meshioplusplus::spherical_shell(
            capi_gen_or(o.radius2, 0.5), o.radius, o.counts[0], o.counts[1],
            capi_gen_or(o.ratio, 1.0), capi_gen_origin(o), o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_pipe(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        if (o.cell_type != MIO_GEN_DEFAULT && o.cell_type != MIO_GEN_HEXAHEDRON)
            capi_gen_no_cell(o, "pipe");
        return new mio_mesh{meshioplusplus::pipe(
            capi_gen_or(o.radius2, 0.5), o.radius, o.height, o.counts[0], o.counts[1], o.counts[2],
            capi_gen_or(o.ratio, 1.0), capi_gen_origin(o), o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_fichera(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        return new mio_mesh{meshioplusplus::fichera(o.counts[0], o.size[0], capi_gen_origin(o),
                                                    capi_gen_solid_cell(o, "fichera"), o.order,
                                                    o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_solid_torus(const mio_gen_opts* opts) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        const mio_gen_opts o = capi_gen_defaults(opts);
        if (o.cell_type != MIO_GEN_DEFAULT && o.cell_type != MIO_GEN_HEXAHEDRON)
            capi_gen_no_cell(o, "solid_torus");
        return new mio_mesh{meshioplusplus::solid_torus(
            o.radius, capi_gen_or(o.radius2, 0.25), o.counts[0], o.counts[1], capi_gen_origin(o),
            o.order, o.max_cells, o.complete != 0)};
    });
}

mio_mesh* mio_reference_cells(const char* types, int32_t pairs) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        std::vector<std::string> names;
        if (types) {
            std::string current;
            for (const char* p = types;; ++p) {
                if (*p == ',' || *p == '\0') {
                    if (!current.empty())
                        names.push_back(current);
                    current.clear();
                    if (*p == '\0')
                        break;
                } else {
                    current.push_back(*p);
                }
            }
        }
        return new mio_mesh{meshioplusplus::reference_cells(names, pairs != 0)};
    });
}

mio_mesh* mio_voxelize(const mio_mesh* mesh, const mio_voxel_opts* opts, int64_t dims_out[3],
                       double origin_out[3], double spacing_out[3], int64_t* num_occupied) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        if (!opts)
            throw meshioplusplus::ReadError(
                "meshio++: voxelize: options are NULL, but exactly one of resolution and "
                "cell_size must be given");
        meshioplusplus::VoxelResult r =
            meshioplusplus::voxelize(mesh->mMesh, capi_voxel_options(*opts));
        for (int k = 0; k < 3; ++k) {
            if (dims_out)
                dims_out[k] = r.mDims[static_cast<std::size_t>(k)];
            if (origin_out)
                origin_out[k] = r.mOrigin[static_cast<std::size_t>(k)];
            if (spacing_out)
                spacing_out[k] = r.mSpacing[static_cast<std::size_t>(k)];
        }
        if (num_occupied)
            *num_occupied = r.mNumOccupied;
        return new mio_mesh{std::move(r.mMesh)};
    });
}

mio_mesh* mio_octree(const mio_mesh* mesh, const mio_octree_opts* opts, mio_octree_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        mio_octree_opts defaults;
        if (!opts) {
            mio_octree_opts_init(&defaults);
            opts = &defaults;
        }
        meshioplusplus::OctreeResult r =
            meshioplusplus::octree(mesh->mMesh, capi_octree_options(*opts));
        if (report) {
            *report = mio_octree_report{};
            report->dimension = r.mDimension;
            report->num_leaves = r.mNumLeaves;
            report->num_hanging = r.mNumHanging;
            report->finest_level = r.mFinestLevel;
            for (int k = 0; k < 3; ++k)
                report->origin[k] = r.mOrigin[static_cast<std::size_t>(k)];
            report->root_size = r.mRootSize;
            report->num_snapped = r.mNumSnapped;
            report->num_refused = r.mNumRefused;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

mio_mesh* mio_octree_mesh(const mio_mesh* mesh, const mio_octree_mesh_opts* opts,
                          mio_octree_mesh_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        capi_octree_mesh_warnings().clear();
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        mio_octree_mesh_opts defaults;
        if (!opts) {
            mio_octree_mesh_opts_init(&defaults);
            opts = &defaults;
        }
        meshioplusplus::OctreeMeshResult r =
            meshioplusplus::octree_mesh(mesh->mMesh, capi_octree_mesh_options(*opts));
        std::string& warnings = capi_octree_mesh_warnings();
        for (const std::string& w : r.mWarnings) {
            if (!warnings.empty())
                warnings += '\n';
            warnings += w;
        }
        if (report) {
            *report = mio_octree_mesh_report{};
            report->output = static_cast<std::int32_t>(r.mOutput);
            report->guarantee_applies = r.mGuaranteeApplies ? 1 : 0;
            for (std::size_t k = 0; k < 3; ++k) {
                report->origin[k] = r.mOrigin[k];
                report->root_cells[k] = r.mRootCells[k];
            }
            report->root_size = r.mRootSize;
            report->max_level = r.mMaxLevel;
            report->finest_level = r.mFinestLevel;
            report->surface_level = r.mSurfaceLevel;
            report->min_size = r.mMinSize;
            report->max_size = r.mMaxSize;
            report->num_leaves = r.mNumLeaves;
            report->num_points = r.mNumPoints;
            report->num_cells = r.mNumCells;
            report->num_hexahedra = r.mNumHexahedra;
            report->num_tetrahedra = r.mNumTetrahedra;
            report->num_inside = r.mNumInside;
            report->num_boundary_cells = r.mNumBoundaryCells;
            report->num_surrogate_faces = r.mNumSurrogateFaces;
            report->num_hanging = r.mNumHanging;
            report->num_cut_points = r.mNumCutPoints;
            report->num_warped = r.mNumWarped;
            report->num_warp_refused = r.mNumWarpRefused;
            report->min_dihedral = r.mMinDihedral;
            report->max_dihedral = r.mMaxDihedral;
            report->min_scaled_jacobian = r.mMinScaledJacobian;
            report->hausdorff_to_skin = r.mHausdorffToSkin;
            report->num_inverted = r.mNumInverted;
            report->num_boundary_edges_open = r.mNumBoundaryEdgesOpen;
            report->num_non_manifold_boundary_edges = r.mNumNonManifoldBoundaryEdges;
            report->num_corners = r.mNumCorners;
            report->num_corners_snapped = r.mNumCornersSnapped;
            report->num_feature_chains = r.mNumFeatureChains;
            report->num_feature_chains_missed = r.mNumFeatureChainsMissed;
            report->num_feature_moves = r.mNumFeatureMoves;
            report->num_optimized_moves = r.mNumOptimizedMoves;
            report->num_flips = r.mNumFlips;
            report->num_winding_rechecks = r.mNumWindingRechecks;
            report->num_winding_skipped = r.mNumWindingSkipped;
            report->num_warnings = static_cast<std::int64_t>(r.mWarnings.size());
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

int64_t mio_octree_mesh_warnings(char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() {
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative warnings buffer length");
        return copy_string(capi_octree_mesh_warnings(), buf, buflen);
    });
}

mio_status mio_surface_watertight_check(const mio_mesh* surface, mio_surface_quality* out) {
    return guarded([&]() -> mio_status {
        if (!surface)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: surface is NULL");
        if (!out)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: out is NULL");
        capi_fill_quality(meshioplusplus::surface_watertight_check(surface->mMesh), out);
        return MIO_OK;
    });
}

mio_status mio_sample_distance(const mio_mesh* surface, const double* points, int64_t n_points,
                               const mio_sdf_opts* opts, double* out) {
    return guarded([&]() -> mio_status {
        if (!surface)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: surface is NULL");
        if (!points || !out)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: points/out is NULL");
        if (n_points < 0)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: n_points is negative");
        if (n_points == 0)
            return MIO_OK;
        // A view, not a copy: sample_distance only reads the query coordinates,
        // and the buffer is the caller's for the duration of the call.
        const meshioplusplus::NDArray query = meshioplusplus::NDArray::MakeView(
            meshioplusplus::DType::Float64, {static_cast<std::size_t>(n_points), std::size_t{3}},
            const_cast<std::byte*>(reinterpret_cast<const std::byte*>(points)));
        const meshioplusplus::SurfaceDistanceOptions options =
            opts ? capi_sdf_options(*opts) : meshioplusplus::SurfaceDistanceOptions{};
        const meshioplusplus::NDArray values =
            meshioplusplus::sample_distance(surface->mMesh, query, options);
        std::memcpy(out, values.Data(), values.Nbytes());
        return MIO_OK;
    });
}

mio_mesh* mio_distance_to_surface(const mio_mesh* query, const mio_mesh* surface,
                                  const mio_sdf_opts* opts, int64_t* num_banded,
                                  mio_surface_quality* quality) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!query || !surface)
            throw meshioplusplus::ReadError("meshio++: query/surface mesh is NULL");
        const meshioplusplus::SurfaceDistanceOptions options =
            opts ? capi_sdf_options(*opts) : meshioplusplus::SurfaceDistanceOptions{};
        meshioplusplus::SurfaceDistanceResult r =
            meshioplusplus::distance_to_surface(query->mMesh, surface->mMesh, options);
        if (num_banded)
            *num_banded = r.mNumBanded;
        capi_fill_quality(r.mQuality, quality);
        return new mio_mesh{std::move(r.mMesh)};
    });
}

static_assert(sizeof(mio_curvature_opts) == 80,
              "mio_curvature_opts grew outside its reserved tail");
static_assert(sizeof(mio_curvature_report) == 136,
              "mio_curvature_report grew outside its reserved tail");

void mio_curvature_opts_init(mio_curvature_opts* opts) {
    if (!opts)
        return;
    *opts = mio_curvature_opts{};  // value-initialized: every opt-in off
    opts->mean = 1;
    opts->gaussian = 1;
    opts->dual_area = MIO_CURVATURE_MIXED_VORONOI;
}

mio_mesh* mio_compute_curvature(const mio_mesh* mesh, const mio_curvature_opts* opts,
                                mio_curvature_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::CurvatureOptions options;
        if (opts) {
            options.mMean = opts->mean != 0;
            options.mGaussian = opts->gaussian != 0;
            if (opts->dual_area != MIO_CURVATURE_MIXED_VORONOI &&
                opts->dual_area != MIO_CURVATURE_BARYCENTRIC)
                throw std::invalid_argument(
                    "meshio++: curvature: dual_area must be MIO_CURVATURE_MIXED_VORONOI or "
                    "MIO_CURVATURE_BARYCENTRIC");
            options.mDualArea = opts->dual_area == MIO_CURVATURE_BARYCENTRIC
                                    ? meshioplusplus::CurvatureDualArea::Barycentric
                                    : meshioplusplus::CurvatureDualArea::MixedVoronoi;
            options.mIncludeBoundary = opts->include_boundary != 0;
            options.mRecordArea = opts->record_area != 0;
            options.mRecordPrincipal = opts->record_principal != 0;
            if (opts->region)
                options.mRegion = opts->region;
        }
        meshioplusplus::CurvatureResult r = meshioplusplus::compute_curvature(mesh->mMesh, options);
        if (report) {
            *report = mio_curvature_report{};
            capi_fill_quality(r.mQuality, &report->quality);
            report->num_boundary = r.mNumBoundary;
            report->num_isolated = r.mNumIsolated;
            report->num_degenerate = r.mNumDegenerate;
            report->total_angle_defect = r.mTotalAngleDefect;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

static_assert(sizeof(mio_normals_opts) == 88, "mio_normals_opts grew outside its reserved tail");
static_assert(sizeof(mio_normals_report) == 144,
              "mio_normals_report grew outside its reserved tail");

void mio_normals_opts_init(mio_normals_opts* opts) {
    if (!opts)
        return;
    *opts = mio_normals_opts{};
    opts->point_normals = 1;
    opts->weight = MIO_SDF_WEIGHT_ANGLE;
    opts->split_angle = 30.0;
}

mio_mesh* mio_compute_normals(const mio_mesh* mesh, const mio_normals_opts* opts,
                              mio_normals_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::NormalsOptions options;
        if (opts) {
            options.mPointNormals = opts->point_normals != 0;
            options.mCellNormals = opts->cell_normals != 0;
            if (opts->weight != MIO_SDF_WEIGHT_ANGLE && opts->weight != MIO_SDF_WEIGHT_AREA)
                throw std::invalid_argument(
                    "meshio++: normals: weight must be MIO_SDF_WEIGHT_ANGLE or "
                    "MIO_SDF_WEIGHT_AREA");
            options.mWeight = opts->weight == MIO_SDF_WEIGHT_AREA
                                  ? meshioplusplus::SdfPseudonormalWeight::Area
                                  : meshioplusplus::SdfPseudonormalWeight::Angle;
            options.mSplit = opts->split != 0;
            options.mSplitAngle = opts->split_angle;
            options.mRecordParentIds = opts->record_parent_ids != 0;
            if (opts->region)
                options.mRegion = opts->region;
        }
        meshioplusplus::NormalsResult r = meshioplusplus::compute_normals(mesh->mMesh, options);
        if (report) {
            *report = mio_normals_report{};
            capi_fill_quality(r.mQuality, &report->quality);
            report->num_isolated = r.mNumIsolated;
            report->num_undefined = r.mNumUndefined;
            report->num_degenerate = r.mNumDegenerate;
            report->num_split_points = r.mNumSplitPoints;
            report->num_added_points = r.mNumAddedPoints;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

static_assert(sizeof(mio_repair_opts) == 80, "mio_repair_opts grew outside its reserved tail");
static_assert(sizeof(mio_repair_report) == 272, "mio_repair_report grew outside its reserved tail");
static_assert(sizeof(mio_shrinkwrap_opts) == 96,
              "mio_shrinkwrap_opts grew outside its reserved tail");
static_assert(sizeof(mio_shrinkwrap_report) == 136,
              "mio_shrinkwrap_report grew outside its reserved tail");
static_assert(sizeof(mio_sobolev_opts) == 88, "mio_sobolev_opts grew outside its reserved tail");
static_assert(sizeof(mio_sobolev_report) == 80,
              "mio_sobolev_report grew outside its reserved tail");

void mio_repair_opts_init(mio_repair_opts* opts) {
    if (!opts)
        return;
    *opts = mio_repair_opts{};
    opts->fix_orientation = 1;
    opts->orient_outward = 1;
    opts->fill_holes = 1;
    opts->split_non_manifold = 1;
    opts->max_hole_edges = 10;
}

mio_mesh* mio_repair(const mio_mesh* mesh, const mio_repair_opts* opts, mio_repair_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::RepairOptions options;
        if (opts) {
            options.mFixOrientation = opts->fix_orientation != 0;
            options.mOrientOutward = opts->orient_outward != 0;
            options.mFillHoles = opts->fill_holes != 0;
            options.mSplitNonManifold = opts->split_non_manifold != 0;
            options.mRecordProvenance = opts->record_provenance != 0;
            options.mMaxHoleEdges = opts->max_hole_edges;
            options.mWeldTolerance = opts->weld_tolerance;
        }
        meshioplusplus::RepairResult r = meshioplusplus::repair(mesh->mMesh, options);
        if (report) {
            *report = mio_repair_report{};
            capi_fill_quality(r.mQualityBefore, &report->quality_before);
            capi_fill_quality(r.mQualityAfter, &report->quality_after);
            report->num_flipped = r.mNumFlipped;
            report->num_components = r.mNumComponents;
            report->largest_component = r.mLargestComponent;
            report->num_oriented_outward = r.mNumOrientedOutward;
            report->num_unorientable = r.mNumUnorientable;
            report->num_vertices_split = r.mNumVerticesSplit;
            report->num_holes_detected = r.mNumHolesDetected;
            report->num_holes_filled = r.mNumHolesFilled;
            report->num_holes_skipped = r.mNumHolesSkipped;
            report->num_faces_added = r.mNumFacesAdded;
            report->num_points_added = r.mNumPointsAdded;
            report->points_welded = r.mPointsWelded;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

void mio_shrinkwrap_opts_init(mio_shrinkwrap_opts* opts) {
    if (!opts)
        return;
    *opts = mio_shrinkwrap_opts{};
    opts->normal_weight = MIO_SDF_WEIGHT_ANGLE;
}

mio_mesh* mio_shrinkwrap(const mio_mesh* mesh, const mio_mesh* target,
                         const mio_shrinkwrap_opts* opts, mio_shrinkwrap_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !target)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::ShrinkwrapOptions options;
        if (opts) {
            if (opts->weights)
                options.mWeights = opts->weights;
            if (opts->target_region)
                options.mTargetRegion = opts->target_region;
            options.mOffset = opts->offset;
            options.mMaxDistance = opts->max_distance;
            options.mGridCellSize = opts->grid_cell_size;
            if (opts->normal_weight < 0 || opts->normal_weight > MIO_SDF_WEIGHT_AREA)
                throw std::invalid_argument(
                    "meshio++: shrinkwrap: normal_weight must be MIO_SDF_WEIGHT_ANGLE or "
                    "MIO_SDF_WEIGHT_AREA");
            options.mNormalWeight = opts->normal_weight == MIO_SDF_WEIGHT_AREA
                                        ? meshioplusplus::SdfPseudonormalWeight::Area
                                        : meshioplusplus::SdfPseudonormalWeight::Angle;
            options.mRecordDistance = opts->record_distance != 0;
            options.mRecordClosestCell = opts->record_closest_cell != 0;
        }
        meshioplusplus::ShrinkwrapResult r =
            meshioplusplus::shrinkwrap(mesh->mMesh, target->mMesh, options);
        if (report) {
            *report = mio_shrinkwrap_report{};
            capi_fill_quality(r.mQuality, &report->quality);
            report->num_projected = r.mNumProjected;
            report->num_missed = r.mNumMissed;
            report->num_skipped = r.mNumSkipped;
            report->max_displacement = r.mMaxDisplacement;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

void mio_sobolev_opts_init(mio_sobolev_opts* opts) {
    if (!opts)
        return;
    *opts = mio_sobolev_opts{};
    opts->max_iterations = 128;
    opts->tolerance = 1e-10;
}

mio_mesh* mio_sobolev_deform(const mio_mesh* mesh, const mio_sobolev_opts* opts,
                             mio_sobolev_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        if (!opts || !opts->array)
            throw std::invalid_argument(
                "meshio++: sobolev_deform: opts->array (the displacement array name) is required");
        meshioplusplus::SobolevOptions options;
        options.mArrayName = opts->array;
        if (opts->fixed_points_array)
            options.mFixedPointsArray = opts->fixed_points_array;
        options.mLengthScale = opts->length_scale;
        options.mTolerance = opts->tolerance;
        options.mMaxIterations = opts->max_iterations;
        options.mFixBoundary = opts->fix_boundary != 0;
        options.mRecordFiltered = opts->record_filtered != 0;
        meshioplusplus::SobolevResult r = meshioplusplus::sobolev_deform(mesh->mMesh, options);
        if (report) {
            *report = mio_sobolev_report{};
            report->num_iterations = r.mNumIterations;
            report->num_fixed = r.mNumFixed;
            report->num_isolated = r.mNumIsolated;
            report->residual = r.mResidual;
            report->max_displacement = r.mMaxDisplacement;
            report->converged = r.mConverged ? 1 : 0;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

int32_t mio_pipeline_has_json(void) {
    return meshioplusplus::pipeline_has_json() ? 1 : 0;
}

// --- analysis and editing (v16.23.0) ----------------------------------------

static_assert(sizeof(mio_feature_edges_opts) == 80,
              "mio_feature_edges_opts grew outside its reserved tail");
static_assert(sizeof(mio_feature_edges_report) == 64,
              "mio_feature_edges_report grew outside its reserved tail");
static_assert(sizeof(mio_segment_opts) == 96, "mio_segment_opts grew outside its reserved tail");
static_assert(sizeof(mio_segment_report) == 64,
              "mio_segment_report grew outside its reserved tail");
static_assert(sizeof(mio_fit_primitives_opts) == 88,
              "mio_fit_primitives_opts grew outside its reserved tail");
static_assert(sizeof(mio_fit_primitives_report) == 96,
              "mio_fit_primitives_report grew outside its reserved tail");
static_assert(sizeof(mio_fit_cad_opts) == 80, "mio_fit_cad_opts grew outside its reserved tail");
static_assert(sizeof(mio_iga_opts) == 96, "mio_iga_opts grew outside its reserved tail");
static_assert(sizeof(mio_mesh_cad_surface_opts) == 104,
              "mio_mesh_cad_surface_opts grew outside its reserved tail");
static_assert(sizeof(mio_mesh_cad_surface_report) == 200,
              "mio_mesh_cad_surface_report grew outside its reserved tail");
static_assert(static_cast<int>(MIO_MESH_CAD_OTHER) ==
                  static_cast<int>(meshioplusplus::MeshCadFaceStatus::Other),
              "mio_mesh_cad_face_status drifted from MeshCadFaceStatus");
static_assert(sizeof(mio_mesh_cad_corner) == 120,
              "mio_mesh_cad_corner grew outside its reserved tail");
static_assert(offsetof(mio_mesh_cad_surface_opts, min_angle_deg) == 56,
              "mio_mesh_cad_surface_opts::min_angle_deg moved");
static_assert(offsetof(mio_mesh_cad_surface_report, num_refinement_points) == 136,
              "mio_mesh_cad_surface_report::num_refinement_points moved");
static_assert(offsetof(mio_mesh_cad_surface_opts, order) == 72,
              "mio_mesh_cad_surface_opts::order moved");
static_assert(offsetof(mio_mesh_cad_surface_opts, pole_rings) == 76 &&
                  offsetof(mio_mesh_cad_surface_opts, grading) == 80 &&
                  offsetof(mio_mesh_cad_surface_opts, curvature_deviation) == 88 &&
                  offsetof(mio_mesh_cad_surface_opts, cells_across_gaps) == 96,
              "mio_mesh_cad_surface_opts: the v17.18.0 fields left their reserved slots");
static_assert(offsetof(mio_mesh_cad_surface_report, num_curved_straightened) == 168,
              "mio_mesh_cad_surface_report::num_curved_straightened moved");
static_assert(sizeof(void*) != 8 ||
                  (sizeof(mio_mesh_cad_surface_opts_ex) == 136 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, background) == 56 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, reparam_condition_number) == 72 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, model_tolerances) == 80 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, virtual_faces) == 84 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, lop_corners) == 88 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, metric_background) == 96 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, metric_array) == 104 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, metric_error) == 112 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, untangle) == 92 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, quads) == 120 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, quad_pairs) == 124 &&
                   offsetof(mio_mesh_cad_surface_opts_ex, reserved) == 128),
              "mio_mesh_cad_surface_opts_ex's layout moved (v17.20.0 took 24 bytes of its "
              "tail, the metric source 24 more)");
static_assert(static_cast<int>(MIO_MESH_CAD_SMOOTH_ODT) ==
                  static_cast<int>(meshioplusplus::MeshCadSmoothing::Odt),
              "mio_mesh_cad_smoothing drifted from MeshCadSmoothing");
static_assert(MIO_MESH_CAD_ANGLE_BINS == meshioplusplus::kMeshCadAngleBins,
              "MIO_MESH_CAD_ANGLE_BINS drifted from kMeshCadAngleBins");
static_assert(sizeof(mio_fit_cad_report) == 192,
              "mio_fit_cad_report grew outside its reserved tail");
static_assert(sizeof(mio_fit_cad_report_ex) == 128,
              "mio_fit_cad_report_ex grew outside its reserved tail");
static_assert(sizeof(mio_refine_patches_opts) == 144,
              "mio_refine_patches_opts grew outside its reserved tail");
static_assert(sizeof(mio_refine_patches_report) == 192,
              "mio_refine_patches_report grew outside its reserved tail");
static_assert(sizeof(mio_fit_cad_opts_ex) == 280,
              "mio_fit_cad_opts_ex grew outside its reserved tail");
static_assert(sizeof(mio_fit_cad_report_ex2) == 176,
              "mio_fit_cad_report_ex2 grew outside its reserved tail");
// The v17.11.0 Phase 3 fields were carved from the front of the reserved tails.
static_assert(offsetof(mio_fit_cad_opts_ex, freeform_general) == 4,
              "mio_fit_cad_opts_ex: freeform_general replaced reserved32");
static_assert(offsetof(mio_fit_cad_opts_ex, freeform_tolerance) == 152,
              "mio_fit_cad_opts_ex: the Phase 3 options start where reserved[16] did");
static_assert(offsetof(mio_fit_cad_opts_ex, sew_tolerance) == 192,
              "mio_fit_cad_opts_ex: sew_tolerance moved");
static_assert(offsetof(mio_fit_cad_opts_ex, reserved) == 216,
              "mio_fit_cad_opts_ex: the reserved tail moved");
static_assert(offsetof(mio_fit_cad_report_ex2, num_freeform_general) == 48,
              "mio_fit_cad_report_ex2: the Phase 3 fields start where reserved[16] did");
static_assert(offsetof(mio_fit_cad_report_ex2, closed) == 104,
              "mio_fit_cad_report_ex2: closed moved");
static_assert(offsetof(mio_fit_cad_report_ex2, refinement_kept) == 112,
              "mio_fit_cad_report_ex2: refinement_kept moved");
static_assert(offsetof(mio_fit_cad_report_ex2, reserved) == 120,
              "mio_fit_cad_report_ex2: the reserved tail moved");
static_assert(sizeof(mio_fit_cad_opts_ex2) == 128,
              "mio_fit_cad_opts_ex2 grew outside its reserved tail");
static_assert(sizeof(mio_fit_cad_report_ex3) == 184,
              "mio_fit_cad_report_ex3 grew outside its reserved tail");
// The constrained primitives (unreleased) were carved from the front of the
// reserved tails.
static_assert(offsetof(mio_fit_cad_opts_ex2, constrain_primitives) == 64,
              "mio_fit_cad_opts_ex2: constrain_primitives left the reserved tail's start");
static_assert(offsetof(mio_fit_cad_opts_ex2, reserved) == 72,
              "mio_fit_cad_opts_ex2: the reserved tail moved");
static_assert(offsetof(mio_fit_cad_report_ex3, num_constrained_primitives) == 104,
              "mio_fit_cad_report_ex3: the counts left the reserved tail's start");
static_assert(offsetof(mio_fit_cad_report_ex3, reserved) == 120,
              "mio_fit_cad_report_ex3: the reserved tail moved");
static_assert(MIO_FREEFORM_AUTO == static_cast<int>(meshioplusplus::FitCadFreeformMode::Auto),
              "mio_fit_cad_freeform_mode must match FitCadFreeformMode");
static_assert(MIO_PARAM_HARMONIC ==
                  static_cast<int>(meshioplusplus::FitCadParameterization::Harmonic),
              "mio_fit_cad_parameterization must match FitCadParameterization");
static_assert(MIO_EV_G1XX == static_cast<int>(meshioplusplus::FitCadExtraordinaryVertex::G1xx),
              "mio_fit_cad_extraordinary_vertex must match FitCadExtraordinaryVertex");
static_assert(MIO_LOCAL_REFINE_THB == static_cast<int>(meshioplusplus::FitCadLocalRefine::Thb),
              "mio_fit_cad_local_refine must match FitCadLocalRefine");
static_assert(MIO_DEMOTION_OTHER == static_cast<int>(meshioplusplus::FitCadDemotionReason::Other),
              "mio_fit_cad_demotion_reason must match FitCadDemotionReason");
static_assert(MIO_REFINE_DEFAULT_SEED == meshioplusplus::kRefineDefaultSeed,
              "MIO_REFINE_DEFAULT_SEED must match kRefineDefaultSeed");
static_assert(MIO_REFINE_ENERGY_ROSE == static_cast<int>(meshioplusplus::RefineStop::EnergyRose),
              "mio_refine_stop must match RefineStop");
static_assert(MIO_FIT_ALL_TYPES == meshioplusplus::kFitAllTypes,
              "MIO_FIT_ALL_TYPES must match kFitAllTypes");
static_assert(MIO_FIT_TORUS == static_cast<int>(meshioplusplus::FitPrimitiveType::Torus),
              "mio_fit_type must match FitPrimitiveType");
static_assert(sizeof(mio_hausdorff_opts) == 80,
              "mio_hausdorff_opts grew outside its reserved tail");
static_assert(sizeof(mio_hausdorff_report) == 152,
              "mio_hausdorff_report grew outside its reserved tail");
static_assert(sizeof(mio_region_selector) == 48,
              "mio_region_selector grew outside its reserved tail");
static_assert(sizeof(mio_periodic_opts) == 192, "mio_periodic_opts grew outside its reserved tail");
static_assert(sizeof(mio_find_interface_opts) == 64,
              "mio_find_interface_opts grew outside its reserved tail");
static_assert(sizeof(mio_contact_pairs_opts) == 48,
              "mio_contact_pairs_opts grew outside its reserved tail");
static_assert(sizeof(mio_split_interface_opts) == 40,
              "mio_split_interface_opts grew outside its reserved tail");

void mio_feature_edges_opts_init(mio_feature_edges_opts* opts) {
    if (!opts)
        return;
    *opts = mio_feature_edges_opts{};
    opts->feature_angle = 30.0;
    opts->feature = 1;
    opts->boundary = 1;
    opts->non_manifold = 1;
    opts->inconsistent = 1;
}

mio_mesh* mio_feature_edges(const mio_mesh* mesh, const mio_feature_edges_opts* opts,
                            mio_feature_edges_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::FeatureEdgeOptions options;
        if (opts) {
            options.mFeatureAngleDeg = opts->feature_angle;
            options.mFeature = opts->feature != 0;
            options.mBoundary = opts->boundary != 0;
            options.mNonManifold = opts->non_manifold != 0;
            options.mInconsistent = opts->inconsistent != 0;
            if (opts->region)
                options.mRegion = opts->region;
        }
        meshioplusplus::FeatureEdgeResult r = meshioplusplus::feature_edges(mesh->mMesh, options);
        if (report) {
            *report = mio_feature_edges_report{};
            report->num_feature = r.mNumFeature;
            report->num_boundary = r.mNumBoundary;
            report->num_non_manifold = r.mNumNonManifold;
            report->num_inconsistent = r.mNumInconsistent;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

void mio_segment_opts_init(mio_segment_opts* opts) {
    if (!opts)
        return;
    *opts = mio_segment_opts{};
    opts->feature_angle = 30.0;
    opts->min_patch_cells = 1;
    opts->method = MIO_SEGMENT_ANGLE;
    opts->tolerance = 1e-3;
}

mio_mesh* mio_segment(const mio_mesh* mesh, const mio_segment_opts* opts,
                      mio_segment_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::SegmentOptions options;
        if (opts) {
            options.mFeatureAngleDeg = opts->feature_angle;
            options.mNormalToleranceDeg = opts->normal_tolerance;
            options.mCurvatureTolerance = opts->curvature_tolerance;
            options.mMinPatchCells = opts->min_patch_cells;
            if (opts->region)
                options.mRegion = opts->region;
            if (opts->prefix && opts->prefix[0] != '\0')
                options.mPrefix = opts->prefix;
            if (opts->method != MIO_SEGMENT_ANGLE && opts->method != MIO_SEGMENT_PRIMITIVE) {
                fail(MIO_ERR_INVALID_ARG,
                     "meshio++: segment: method must be MIO_SEGMENT_ANGLE or "
                     "MIO_SEGMENT_PRIMITIVE");
                return nullptr;
            }
            options.mMethod = static_cast<meshioplusplus::SegmentMethod>(opts->method);
            // A struct zero-filled before the field existed takes the default.
            if (opts->tolerance > 0.0)
                options.mTolerance = opts->tolerance;
        }
        meshioplusplus::SegmentResult r = meshioplusplus::segment(mesh->mMesh, options);
        if (report) {
            *report = mio_segment_report{};
            report->num_patches = r.mNumPatches;
            report->num_merged = r.mNumMerged;
            report->num_feature_edges = r.mNumFeatureEdges;
            report->num_boundary_edges = r.mNumBoundaryEdges;
            report->num_non_manifold_edges = r.mNumNonManifoldEdges;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

void mio_fit_primitives_opts_init(mio_fit_primitives_opts* opts) {
    if (!opts)
        return;
    *opts = mio_fit_primitives_opts{};
    opts->tolerance = 1e-3;
    opts->min_inlier_fraction = 0.9;
    opts->relative = 1;
    opts->types = MIO_FIT_ALL_TYPES;
    opts->max_iterations = 50;
}

mio_mesh* mio_fit_primitives(const mio_mesh* mesh, const mio_fit_primitives_opts* opts,
                             mio_fit_primitives_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::FitPrimitivesOptions options;
        if (opts) {
            if (opts->patch_data)
                options.mPatchData = opts->patch_data;
            options.mTolerance = opts->tolerance;
            options.mMinInlierFraction = opts->min_inlier_fraction;
            options.mRelative = opts->relative != 0;
            options.mTypes = opts->types;
            options.mMaxIterations = opts->max_iterations;
        }
        meshioplusplus::FitPrimitivesResult r =
            meshioplusplus::fit_primitives(mesh->mMesh, options);
        if (report) {
            *report = mio_fit_primitives_report{};
            report->num_patches = r.mNumPatches;
            report->num_plane = r.mNumPlane;
            report->num_sphere = r.mNumSphere;
            report->num_cylinder = r.mNumCylinder;
            report->num_cone = r.mNumCone;
            report->num_torus = r.mNumTorus;
            report->num_unfitted = r.mNumUnfitted;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

void mio_refine_patches_opts_init(mio_refine_patches_opts* opts) {
    if (!opts)
        return;
    *opts = mio_refine_patches_opts{};
    const meshioplusplus::RefinePatchesOptions d;
    opts->tolerance = d.mTolerance;
    opts->normal_weight = d.mNormalWeight;
    opts->epsilon = d.mEpsilon;
    opts->feature_angle = d.mFeatureAngleDeg;
    opts->max_iterations = d.mMaxIterations;
    opts->min_patch_cells = d.mMinPatchCells;
    opts->ransac_seed = d.mRansacSeed;
    opts->types = d.mTypes;
    opts->relative = d.mRelative ? 1 : 0;
    opts->allow_split = d.mAllowSplit ? 1 : 0;
    opts->allow_merge = d.mAllowMerge ? 1 : 0;
    opts->prefix = nullptr;  // "patch"
}

namespace {

// The C refine options as the C++ ones; a NULL or "" prefix keeps "patch".
meshioplusplus::RefinePatchesOptions capi_refine_patches_options(
    const mio_refine_patches_opts& rOpts) {
    meshioplusplus::RefinePatchesOptions options;
    options.mTolerance = rOpts.tolerance;
    options.mNormalWeight = rOpts.normal_weight;
    options.mEpsilon = rOpts.epsilon;
    options.mFeatureAngleDeg = rOpts.feature_angle;
    options.mMaxIterations = rOpts.max_iterations;
    options.mMinPatchCells = rOpts.min_patch_cells;
    options.mRansacSeed = rOpts.ransac_seed;
    options.mTypes = rOpts.types;
    options.mRelative = rOpts.relative != 0;
    options.mAllowSplit = rOpts.allow_split != 0;
    options.mAllowMerge = rOpts.allow_merge != 0;
    if (rOpts.prefix && rOpts.prefix[0] != '\0')
        options.mPrefix = rOpts.prefix;
    return options;
}

}  // namespace

mio_mesh* mio_refine_patches(const mio_mesh* mesh, const mio_refine_patches_opts* opts,
                             mio_refine_patches_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        const meshioplusplus::RefinePatchesOptions options =
            opts ? capi_refine_patches_options(*opts) : meshioplusplus::RefinePatchesOptions{};
        meshioplusplus::RefinePatchesResult r =
            meshioplusplus::refine_patches(mesh->mMesh, options);
        if (report) {
            *report = mio_refine_patches_report{};
            report->num_patches = r.mNumPatches;
            report->num_iterations = r.mNumIterations;
            report->num_splits = r.mNumSplits;
            report->num_merges = r.mNumMerges;
            report->num_absorbed = r.mNumAbsorbed;
            report->num_primitive = r.mNumPrimitive;
            report->num_freeform = r.mNumFreeform;
            report->num_facets = r.mNumFacets;
            report->tolerance = r.mTolerance;
            report->explained_area_fraction = r.mExplainedAreaFraction;
            report->initial_explained_area_fraction = r.mInitialExplainedAreaFraction;
            report->hausdorff_max = r.mHausdorffMax;
            report->hausdorff_rms = r.mHausdorffRms;
            report->initial_energy = r.mInitialEnergy;
            report->energy = r.mEnergy;
            report->stop = static_cast<std::int32_t>(r.mStop);
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

void mio_fit_cad_opts_init(mio_fit_cad_opts* opts) {
    if (!opts)
        return;
    *opts = mio_fit_cad_opts{};
    opts->max_demotion_rounds = 8;
    opts->exact_edges = 1;
    opts->types = MIO_FIT_ALL_TYPES;
    opts->freeform = 1;
    opts->freeform_degree = 3;
    opts->freeform_max_control = 64;
    opts->freeform_split = 1;
    opts->fillets = 1;
    opts->threads = 1;
    opts->merge_faces = 1;
    opts->backend = MIO_FIT_CAD_AUTO;
}

void mio_fit_cad_opts_ex_init(mio_fit_cad_opts_ex* opts) {
    if (!opts)
        return;
    *opts = mio_fit_cad_opts_ex{};
    opts->refine_segmentation = 0;
    mio_refine_patches_opts_init(&opts->refine);
    const meshioplusplus::FitCadOptions defaults;
    opts->freeform_general = defaults.mFreeformGeneral ? 1 : 0;
    opts->freeform_tolerance = defaults.mFreeformTolerance;
    opts->freeform_smoothing = defaults.mFreeformSmoothing;
    opts->max_freeform_poles = defaults.mMaxFreeformPoles;
    opts->max_freeform_iterations = defaults.mMaxFreeformIterations;
    opts->max_freeform_splits = defaults.mMaxFreeformSplits;
    opts->sew_tolerance = defaults.mSewTolerance;
    opts->extend_surfaces = defaults.mExtendSurfaces ? 1 : 0;
    opts->snap_junctions = defaults.mSnapJunctions ? 1 : 0;
    opts->face_deviation = defaults.mFaceDeviation ? 1 : 0;
}

void mio_fit_cad_opts_ex2_init(mio_fit_cad_opts_ex2* opts) {
    if (!opts)
        return;
    *opts = mio_fit_cad_opts_ex2{};
    const meshioplusplus::FitCadOptions defaults;
    opts->freeform_mode = static_cast<std::int32_t>(defaults.mFreeformMode);
    opts->freeform_parameterization = static_cast<std::int32_t>(defaults.mFreeformParameterization);
    opts->extraordinary_vertex = static_cast<std::int32_t>(defaults.mExtraordinaryVertex);
    opts->local_refine = static_cast<std::int32_t>(defaults.mLocalRefine);
    opts->freeform_g1_tolerance_deg = defaults.mFreeformG1ToleranceDeg;
    opts->freeform_max_control_points = defaults.mFreeformMaxControlPoints;
    opts->layout_target_patch_size = defaults.mLayoutTargetPatchSize;
    opts->layout_max_valence_deviation = defaults.mLayoutMaxValenceDeviation;
    opts->resegment_rounds = defaults.mResegmentRounds;
    opts->layout_snap_to_features = defaults.mLayoutSnapToFeatures ? 1 : 0;
    opts->validate_by_surface_sampling = defaults.mValidateBySurfaceSampling ? 1 : 0;
    opts->constrain_primitives = defaults.mConstrainPrimitives ? 1 : 0;
}

void mio_iga_opts_init(mio_iga_opts* opts) {
    if (!opts)
        return;
    *opts = mio_iga_opts{};
    for (int d = 0; d < 3; ++d)
        opts->subdivisions[d] = 1;
    opts->trim_samples = 16;
    opts->edges = 1;
    opts->trim = 1;
}

}  // extern "C"

namespace {

meshioplusplus::IgaRefineOptions capi_iga_refine_options(const mio_iga_opts* pOpts) {
    mio_iga_opts o;
    mio_iga_opts_init(&o);
    if (pOpts)
        o = *pOpts;
    meshioplusplus::IgaRefineOptions r;
    for (std::size_t d = 0; d < 3; ++d) {
        r.mSubdivisions[d] = o.subdivisions[d];
        r.mElevate[d] = o.elevate[d];
    }
    return r;
}

meshioplusplus::IgaSampleOptions capi_iga_sample_options(const mio_iga_opts* pOpts) {
    mio_iga_opts o;
    mio_iga_opts_init(&o);
    if (pOpts)
        o = *pOpts;
    meshioplusplus::IgaSampleOptions s;
    s.mSubdivisions = o.sample_subdivisions;
    s.mTrimSamples = o.trim_samples;
    s.mEdges = o.edges != 0;
    s.mTrim = o.trim != 0;
    return s;
}

// A new handle of `rFormat` holding `rInfo`.
template <class Info>
mio_format_info* capi_iga_handle(const std::string& rFormat, Info&& rInfo) {
    auto handle = std::make_unique<mio_format_info>();
    handle->mFormat = rFormat;
    handle->mInfo = std::forward<Info>(rInfo);
    if (std::holds_alternative<meshioplusplus::BrepInfo>(handle->mInfo))
        capi_index_brep(*handle);
    return handle.release();
}

const char* const kCapiIgaNotNurbs = "meshio++: iga: the side channel is not a NURBS model";

// This thread's last mesh_cad_surface call: its per-face records and warnings.
struct CapiMeshCadState {
    std::vector<meshioplusplus::MeshCadFaceRecord> mFaces;
    std::string mWarnings;
    std::array<std::int64_t, meshioplusplus::kMeshCadAngleBins> mHistogram{};
    std::vector<meshioplusplus::MeshCadCorner> mCorners;
    meshioplusplus::MeshCadSizeReport mSizes;
    meshioplusplus::MeshCadSurfaceExtrasReport mExtras;
};
CapiMeshCadState& capi_mesh_cad_state() {
    thread_local CapiMeshCadState state;
    return state;
}

}  // namespace

extern "C" {

mio_mesh* mio_iga_refine(const mio_format_info* info, const mio_iga_opts* opts,
                         mio_format_info** refined) {
    if (refined)
        *refined = nullptr;
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!info)
            throw std::invalid_argument("meshio++: format info handle is NULL");
        const meshioplusplus::IgaRefineOptions ro = capi_iga_refine_options(opts);
        const meshioplusplus::IgaSampleOptions so = capi_iga_sample_options(opts);
        if (const auto* p_nurbs = std::get_if<meshioplusplus::NurbsInfo>(&info->mInfo)) {
            const meshioplusplus::NurbsInfo r = meshioplusplus::iga_refine(*p_nurbs, ro);
            meshioplusplus::NurbsInfo mapped;
            auto mesh =
                std::make_unique<mio_mesh>(mio_mesh{meshioplusplus::iga_sample(r, so, &mapped)});
            if (refined)
                *refined = capi_iga_handle(info->mFormat, std::move(mapped));
            return mesh.release();
        }
        if (const auto* p_brep = std::get_if<meshioplusplus::BrepInfo>(&info->mInfo)) {
            meshioplusplus::BrepInfo r = meshioplusplus::iga_refine(*p_brep, ro);
            auto mesh = std::make_unique<mio_mesh>(mio_mesh{meshioplusplus::iga_sample(r, so)});
            if (refined)
                *refined = capi_iga_handle(info->mFormat, std::move(r));
            return mesh.release();
        }
        throw std::invalid_argument(kCapiIgaNotNurbs);
    });
}

mio_mesh* mio_iga_sample(const mio_format_info* info, const mio_iga_opts* opts,
                         mio_format_info** mapped) {
    if (mapped)
        *mapped = nullptr;
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!info)
            throw std::invalid_argument("meshio++: format info handle is NULL");
        const meshioplusplus::IgaSampleOptions so = capi_iga_sample_options(opts);
        if (const auto* p_nurbs = std::get_if<meshioplusplus::NurbsInfo>(&info->mInfo)) {
            meshioplusplus::NurbsInfo with_map;
            auto mesh = std::make_unique<mio_mesh>(
                mio_mesh{meshioplusplus::iga_sample(*p_nurbs, so, &with_map)});
            if (mapped)
                *mapped = capi_iga_handle(info->mFormat, std::move(with_map));
            return mesh.release();
        }
        if (const auto* p_brep = std::get_if<meshioplusplus::BrepInfo>(&info->mInfo)) {
            auto mesh =
                std::make_unique<mio_mesh>(mio_mesh{meshioplusplus::iga_sample(*p_brep, so)});
            if (mapped)
                *mapped = capi_iga_handle(info->mFormat, meshioplusplus::BrepInfo(*p_brep));
            return mesh.release();
        }
        throw std::invalid_argument(kCapiIgaNotNurbs);
    });
}

mio_iga_evaluation* mio_iga_evaluate(const mio_format_info* info, int64_t patch,
                                     const double* params, int64_t num_points, int32_t derivatives,
                                     int32_t trim_samples) {
    return guarded_ptr(static_cast<mio_iga_evaluation*>(nullptr), [&]() -> mio_iga_evaluation* {
        if (!info)
            throw std::invalid_argument("meshio++: format info handle is NULL");
        if (!params && num_points > 0)
            throw std::invalid_argument("meshio++: params is NULL");
        if (patch < 0 || num_points < 0)
            throw std::invalid_argument("meshio++: negative patch or point count");
        const auto index = static_cast<std::size_t>(patch);
        const auto n = static_cast<std::size_t>(num_points);
        const auto make_params = [&](std::size_t Dim) {
            meshioplusplus::NDArray a(meshioplusplus::DType::Float64, {n, Dim});
            if (n)
                std::copy(params, params + n * Dim, a.As<double>());
            return a;
        };
        auto out = std::make_unique<mio_iga_evaluation>();
        if (const auto* p_nurbs = std::get_if<meshioplusplus::NurbsInfo>(&info->mInfo)) {
            if (index >= p_nurbs->mPatches.size())
                throw std::invalid_argument("meshio++: iga: patch out of range");
            out->mEval = meshioplusplus::iga_evaluate(
                *p_nurbs, index, make_params(p_nurbs->mPatches[index].mDegrees.size()),
                derivatives);
            return out.release();
        }
        if (const auto* p_brep = std::get_if<meshioplusplus::BrepInfo>(&info->mInfo)) {
            out->mFace = true;
            out->mEval = meshioplusplus::iga_evaluate(*p_brep, index, make_params(2), derivatives,
                                                      trim_samples);
            return out.release();
        }
        throw std::invalid_argument(kCapiIgaNotNurbs);
    });
}

int64_t mio_iga_evaluation_num_fields(const mio_iga_evaluation* evaluation) {
    if (!evaluation)
        return -1;
    return static_cast<int64_t>(evaluation->mEval.mFields.size());
}

mio_status mio_iga_evaluation_array(const mio_iga_evaluation* evaluation, int32_t which,
                                    int64_t index, const void** data, mio_dtype* dtype,
                                    int32_t* ndim, int64_t* shape) {
    return guarded([&]() -> mio_status {
        if (!evaluation)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: evaluation handle is NULL");
        const meshioplusplus::IgaEvaluation& e = evaluation->mEval;
        const meshioplusplus::NDArray* a = nullptr;
        switch (which) {
            case MIO_IGA_POINTS:
                a = &e.mPoints;
                break;
            case MIO_IGA_DERIVATIVES:
                a = e.mDerivatives.Size() ? &e.mDerivatives : nullptr;
                break;
            case MIO_IGA_NORMALS:
                a = e.mNormals.Size() ? &e.mNormals : nullptr;
                break;
            case MIO_IGA_INSIDE:
                a = evaluation->mFace ? &e.mInside : nullptr;
                break;
            case MIO_IGA_FIELD:
                if (index >= 0 && static_cast<std::size_t>(index) < e.mFields.size())
                    a = &e.mFields[static_cast<std::size_t>(index)];
                break;
            default:
                break;
        }
        if (!a)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: iga: the evaluation holds no such array");
        return array_out(*a, data, dtype, ndim, shape);
    });
}

void mio_iga_evaluation_free(mio_iga_evaluation* evaluation) {
    delete evaluation;
}

void mio_mesh_cad_surface_opts_init(mio_mesh_cad_surface_opts* opts) {
    if (!opts)
        return;
    const meshioplusplus::MeshCadSurfaceOptions d;
    *opts = mio_mesh_cad_surface_opts{};
    opts->target_size = d.mTargetSize;
    opts->max_chordal_deviation = d.mMaxChordalDeviation;
    opts->tolerance = d.mGeometricTolerance;
    opts->min_points_per_edge = d.mMinPointsPerEdge;
    opts->max_points = d.mMaxPoints;
    opts->record_face_ids = d.mRecordFaceIds ? 1 : 0;
    opts->edges = d.mEdges ? 1 : 0;
    opts->match_unreferenced_trims = d.mMatchUnreferencedTrims ? 1 : 0;
    opts->strict = d.mStrict ? 1 : 0;
    opts->min_angle_deg = d.mMinAngleDeg;
    opts->smooth_iterations = static_cast<int32_t>(d.mSmoothIterations);
    opts->smooth_method = static_cast<int32_t>(d.mSmoothMethod);
    opts->order = d.mOrder;
    opts->pole_rings = d.mPoleRings;
    opts->grading = d.mGrading;
    opts->curvature_deviation = d.mCurvatureDeviation;
    opts->cells_across_gaps = d.mCellsAcrossGaps;
}

void mio_mesh_cad_surface_opts_ex_init(mio_mesh_cad_surface_opts_ex* opts) {
    if (!opts)
        return;
    *opts = mio_mesh_cad_surface_opts_ex{};
    opts->metric = MIO_MESH_CAD_METRIC_MEAN;
}

mio_mesh* mio_mesh_cad_surface(const mio_format_info* info, const mio_mesh_cad_surface_opts* opts,
                               mio_mesh_cad_surface_report* report) {
    return mio_mesh_cad_surface2(info, opts, nullptr, report);
}

mio_mesh* mio_mesh_cad_surface2(const mio_format_info* info, const mio_mesh_cad_surface_opts* opts,
                                const mio_mesh_cad_surface_opts_ex* opts_ex,
                                mio_mesh_cad_surface_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        CapiMeshCadState& state = capi_mesh_cad_state();
        state.mFaces.clear();
        state.mWarnings.clear();
        state.mHistogram.fill(0);
        state.mCorners.clear();
        state.mSizes = meshioplusplus::MeshCadSizeReport{};
        state.mExtras = meshioplusplus::MeshCadSurfaceExtrasReport{};
        if (!info)
            throw std::invalid_argument("meshio++: format info handle is NULL");
        if (!opts)
            throw std::invalid_argument(
                "meshio++: mesh_cad_surface: options are required (target_size)");
        const auto* p_brep = std::get_if<meshioplusplus::BrepInfo>(&info->mInfo);
        if (!p_brep)
            throw std::invalid_argument(
                "meshio++: mesh_cad_surface: the side channel is not a B-rep (read a "
                "\"cad_json\", \"step\" or \"iges\" file)");
        meshioplusplus::MeshCadSurfaceOptions o;
        o.mTargetSize = opts->target_size;
        o.mMaxChordalDeviation = opts->max_chordal_deviation;
        o.mGeometricTolerance = opts->tolerance;
        o.mMinPointsPerEdge = opts->min_points_per_edge;
        o.mMaxPoints = opts->max_points;
        o.mRecordFaceIds = opts->record_face_ids != 0;
        o.mEdges = opts->edges != 0;
        o.mMatchUnreferencedTrims = opts->match_unreferenced_trims != 0;
        o.mStrict = opts->strict != 0;
        if (opts->smooth_method != MIO_MESH_CAD_SMOOTH_LAPLACIAN &&
            opts->smooth_method != MIO_MESH_CAD_SMOOTH_ODT)
            throw std::invalid_argument(
                "meshio++: mesh_cad_surface: smooth_method must be MIO_MESH_CAD_SMOOTH_LAPLACIAN "
                "or MIO_MESH_CAD_SMOOTH_ODT");
        o.mMinAngleDeg = opts->min_angle_deg;
        o.mSmoothIterations = opts->smooth_iterations;
        o.mSmoothMethod = static_cast<meshioplusplus::MeshCadSmoothing>(opts->smooth_method);
        o.mOrder = opts->order == 0 ? 1 : opts->order;
        o.mPoleRings = opts->pole_rings;
        o.mGrading = opts->grading;
        o.mCurvatureDeviation = opts->curvature_deviation;
        o.mCellsAcrossGaps = opts->cells_across_gaps;
        meshioplusplus::MeshCadSizeSources sources;
        meshioplusplus::MeshCadSurfaceExtras extras;
        if (opts_ex) {
            o.mReparamConditionNumber = opts_ex->reparam_condition_number;
            extras.mModelTolerances = opts_ex->model_tolerances != 0;
            extras.mVirtualFaces = opts_ex->virtual_faces != 0;
            extras.mLopCorners = opts_ex->lop_corners != 0;
            extras.mUntangle = opts_ex->untangle != 0;
            // The core refuses a value out of range by name.
            extras.mQuads = static_cast<meshioplusplus::MeshCadQuads>(opts_ex->quads);
            extras.mQuadPairs = static_cast<meshioplusplus::MeshCadQuadPairs>(opts_ex->quad_pairs);
            if (opts_ex->metric != MIO_MESH_CAD_METRIC_MEAN &&
                opts_ex->metric != MIO_MESH_CAD_METRIC_VARYING)
                throw std::invalid_argument(
                    "meshio++: mesh_cad_surface: metric must be MIO_MESH_CAD_METRIC_MEAN or "
                    "MIO_MESH_CAD_METRIC_VARYING");
            o.mMetric = static_cast<meshioplusplus::MeshCadMetric>(opts_ex->metric);
            sources.mVirtualEdges = opts_ex->virtual_edges != 0;
            const auto pairs = [](const int64_t* pIndex, const double* pSize, int64_t N,
                                  const char* pName) {
                std::vector<std::pair<std::int64_t, double>> v;
                if (N < 0 || (N > 0 && (!pIndex || !pSize)))
                    throw std::invalid_argument(std::string("meshio++: mesh_cad_surface: ") +
                                                pName + ": bad arrays");
                for (int64_t i = 0; i < N; ++i)
                    v.emplace_back(pIndex[i], pSize[i]);
                return v;
            };
            sources.mFaceSizes = pairs(opts_ex->face_size_faces, opts_ex->face_sizes,
                                       opts_ex->num_face_sizes, "face_sizes");
            sources.mEdgeSizes = pairs(opts_ex->edge_size_edges, opts_ex->edge_sizes,
                                       opts_ex->num_edge_sizes, "edge_sizes");
            if (opts_ex->background) {
                sources.pBackground = &opts_ex->background->mMesh;
                sources.mBackgroundArray = opts_ex->background_array && *opts_ex->background_array
                                               ? opts_ex->background_array
                                               : "size";
            }
            if (opts_ex->metric_background)
                sources.pMetricBackground = &opts_ex->metric_background->mMesh;
            if (opts_ex->metric_array)
                sources.mMetricArray = opts_ex->metric_array;
            sources.mMetricError = opts_ex->metric_error;
        }
        meshioplusplus::MeshCadSurfaceResult r = meshioplusplus::mesh_cad_surface(
            *p_brep, o, sources, extras, &state.mSizes, &state.mExtras);
        for (const std::string& w : r.mWarnings) {
            if (!state.mWarnings.empty())
                state.mWarnings += '\n';
            state.mWarnings += w;
        }
        state.mFaces = r.mFaces;
        state.mHistogram = r.mAngleHistogram;
        double min_input = 180.0;
        for (const auto& r_f : r.mFaces) {
            min_input = std::min(min_input, r_f.mMinInputAngleDeg);
            state.mCorners.insert(state.mCorners.end(), r_f.mSmallAngleCorners.begin(),
                                  r_f.mSmallAngleCorners.end());
        }
        if (report) {
            *report = mio_mesh_cad_surface_report{};
            report->num_points = r.mNumPoints;
            report->num_triangles = r.mNumTriangles;
            report->num_edge_segments = r.mNumEdgeSegments;
            report->num_faces = r.mNumFaces;
            report->num_failed = r.mNumFailed;
            report->num_edges = r.mNumEdges;
            report->num_open_edges = r.mNumOpenEdges;
            report->num_non_manifold_edges = r.mNumNonManifoldEdges;
            report->num_misoriented_edges = r.mNumMisorientedEdges;
            report->num_degenerate_dropped = r.mNumDegenerateDropped;
            report->tolerance = r.mTolerance;
            report->min_angle_deg = r.mMinAngleDeg;
            report->max_angle_deg = r.mMaxAngleDeg;
            report->min_edge = r.mMinEdge;
            report->max_edge = r.mMaxEdge;
            report->max_edge_gap = r.mMaxEdgeGap;
            report->num_warnings = static_cast<std::int64_t>(r.mWarnings.size());
            report->num_refinement_points = r.mNumRefinementPoints;
            report->num_faces_below_bound = r.mNumFacesBelowBound;
            report->num_smoothed = r.mNumSmoothed;
            report->min_input_angle_deg = min_input;
            report->num_curved_straightened = r.mNumCurvedStraightened;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

int32_t mio_mesh_cad_surface_face_status(int64_t face) {
    const auto& faces = capi_mesh_cad_state().mFaces;
    if (face < 0 || static_cast<std::size_t>(face) >= faces.size())
        return -1;
    return static_cast<int32_t>(faces[static_cast<std::size_t>(face)].mStatus);
}

int64_t mio_mesh_cad_surface_face_reason(int64_t face, char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() {
        const auto& faces = capi_mesh_cad_state().mFaces;
        if (face < 0 || static_cast<std::size_t>(face) >= faces.size())
            throw std::invalid_argument("meshio++: mesh_cad_surface: face index out of range");
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative reason buffer length");
        return copy_string(faces[static_cast<std::size_t>(face)].mReason, buf, buflen);
    });
}

int64_t mio_mesh_cad_surface_warnings(char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() {
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative warnings buffer length");
        return copy_string(capi_mesh_cad_state().mWarnings, buf, buflen);
    });
}

int64_t mio_mesh_cad_surface_angle_histogram(int64_t* bins, int64_t n) {
    return guarded_ptr(int64_t{-1}, [&]() {
        if (n < 0 || (n > 0 && !bins))
            throw std::invalid_argument("meshio++: mesh_cad_surface: bad histogram buffer");
        const auto& r_h = capi_mesh_cad_state().mHistogram;
        const auto count = std::min<std::int64_t>(n, static_cast<std::int64_t>(r_h.size()));
        for (std::int64_t b = 0; b < count; ++b)
            bins[b] = r_h[static_cast<std::size_t>(b)];
        return static_cast<int64_t>(r_h.size());
    });
}

int64_t mio_mesh_cad_surface_face_size_limits(int64_t face, int64_t* edges, int64_t n) {
    return guarded_ptr(int64_t{-1}, [&]() -> int64_t {
        const auto& r_l = capi_mesh_cad_state().mSizes.mLimitedByEdges;
        if (face < 0 || static_cast<std::size_t>(face) >= r_l.size())
            return -1;
        if (n < 0 || (n > 0 && !edges))
            throw std::invalid_argument("meshio++: mesh_cad_surface: bad edges buffer");
        const auto& r_e = r_l[static_cast<std::size_t>(face)];
        const auto count = std::min<std::int64_t>(n, static_cast<std::int64_t>(r_e.size()));
        std::copy(r_e.begin(), r_e.begin() + count, edges);
        return static_cast<int64_t>(r_e.size());
    });
}

int64_t mio_mesh_cad_surface_num_virtual_edges(void) {
    return static_cast<int64_t>(capi_mesh_cad_state().mSizes.mVirtualEdges.size());
}

int64_t mio_mesh_cad_surface_virtual_edge(int64_t index, int64_t* edges, int64_t n) {
    return guarded_ptr(int64_t{-1}, [&]() -> int64_t {
        const auto& r_v = capi_mesh_cad_state().mSizes.mVirtualEdges;
        if (index < 0 || static_cast<std::size_t>(index) >= r_v.size())
            return -1;
        if (n < 0 || (n > 0 && !edges))
            throw std::invalid_argument("meshio++: mesh_cad_surface: bad edges buffer");
        const auto& r_c = r_v[static_cast<std::size_t>(index)];
        const auto count = std::min<std::int64_t>(n, static_cast<std::int64_t>(r_c.size()));
        std::copy(r_c.begin(), r_c.begin() + count, edges);
        return static_cast<int64_t>(r_c.size());
    });
}

int64_t mio_mesh_cad_surface_num_reparametrised(void) {
    return static_cast<int64_t>(capi_mesh_cad_state().mExtras.mReparametrised.size());
}

mio_status mio_mesh_cad_surface_reparametrised(int64_t index, int64_t* face,
                                               double* condition_number,
                                               double* condition_number_disc, int32_t* method) {
    return guarded([&]() -> mio_status {
        const auto& r_v = capi_mesh_cad_state().mExtras.mReparametrised;
        if (index < 0 || static_cast<std::size_t>(index) >= r_v.size())
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: mesh_cad_surface: reparametrised index out of range");
        const auto& r_f = r_v[static_cast<std::size_t>(index)];
        if (face)
            *face = r_f.mFace;
        if (condition_number)
            *condition_number = r_f.mConditionNumber;
        if (condition_number_disc)
            *condition_number_disc = r_f.mConditionNumberDisc;
        if (method)
            *method = r_f.mTaken ? r_f.mMethod : 0;
        return MIO_OK;
    });
}

int64_t mio_mesh_cad_surface_reparam_reason(int64_t index, char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() {
        const auto& r_v = capi_mesh_cad_state().mExtras.mReparametrised;
        if (index < 0 || static_cast<std::size_t>(index) >= r_v.size())
            throw std::invalid_argument(
                "meshio++: mesh_cad_surface: reparametrised index out of range");
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative reason buffer length");
        return copy_string(r_v[static_cast<std::size_t>(index)].mReason, buf, buflen);
    });
}

int64_t mio_mesh_cad_surface_num_merged_faces(void) {
    return static_cast<int64_t>(capi_mesh_cad_state().mExtras.mMergedFaces.size());
}

mio_status mio_mesh_cad_surface_merged_face(int64_t index, int64_t* sliver, int64_t* host,
                                            double* width, int32_t* merged) {
    return guarded([&]() -> mio_status {
        const auto& r_v = capi_mesh_cad_state().mExtras.mMergedFaces;
        if (index < 0 || static_cast<std::size_t>(index) >= r_v.size())
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: mesh_cad_surface: merged face index out of range");
        const auto& r_m = r_v[static_cast<std::size_t>(index)];
        if (sliver)
            *sliver = r_m.mSliver;
        if (host)
            *host = r_m.mHost;
        if (width)
            *width = r_m.mWidth;
        if (merged)
            *merged = r_m.mMerged ? 1 : 0;
        return MIO_OK;
    });
}

int64_t mio_mesh_cad_surface_merged_face_edges(int64_t index, int64_t* edges, int64_t n) {
    return guarded_ptr(int64_t{-1}, [&]() -> int64_t {
        const auto& r_v = capi_mesh_cad_state().mExtras.mMergedFaces;
        if (index < 0 || static_cast<std::size_t>(index) >= r_v.size())
            return -1;
        if (n < 0 || (n > 0 && !edges))
            throw std::invalid_argument("meshio++: mesh_cad_surface: bad edges buffer");
        const auto& r_e = r_v[static_cast<std::size_t>(index)].mEdges;
        const auto count = std::min<std::int64_t>(n, static_cast<std::int64_t>(r_e.size()));
        std::copy(r_e.begin(), r_e.begin() + count, edges);
        return static_cast<int64_t>(r_e.size());
    });
}

int64_t mio_mesh_cad_surface_merged_face_reason(int64_t index, char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() {
        const auto& r_v = capi_mesh_cad_state().mExtras.mMergedFaces;
        if (index < 0 || static_cast<std::size_t>(index) >= r_v.size())
            throw std::invalid_argument(
                "meshio++: mesh_cad_surface: merged face index out of range");
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative reason buffer length");
        return copy_string(r_v[static_cast<std::size_t>(index)].mReason, buf, buflen);
    });
}

int64_t mio_mesh_cad_surface_num_tolerance_uses(void) {
    return static_cast<int64_t>(capi_mesh_cad_state().mExtras.mToleranceUses.size());
}

mio_status mio_mesh_cad_surface_tolerance_use(int64_t index, int64_t* edge, double* declared,
                                              double* gap, int32_t* where) {
    return guarded([&]() -> mio_status {
        const auto& r_v = capi_mesh_cad_state().mExtras.mToleranceUses;
        if (index < 0 || static_cast<std::size_t>(index) >= r_v.size())
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: mesh_cad_surface: tolerance use index out of range");
        const auto& r_t = r_v[static_cast<std::size_t>(index)];
        if (edge)
            *edge = r_t.mEdge;
        if (declared)
            *declared = r_t.mDeclared;
        if (gap)
            *gap = r_t.mGap;
        if (where)
            *where = r_t.mWhere == "closed_start" ? MIO_MESH_CAD_TOLERANCE_CLOSED_START
                     : r_t.mWhere == "corner"     ? MIO_MESH_CAD_TOLERANCE_CORNER
                     : r_t.mWhere == "pairing"    ? MIO_MESH_CAD_TOLERANCE_PAIRING
                                                  : MIO_MESH_CAD_TOLERANCE_PROJECTION;
        return MIO_OK;
    });
}

static_assert(sizeof(mio_project_to_cad_opts) == 64, "mio_project_to_cad_opts layout");
static_assert(offsetof(mio_project_to_cad_opts, trim_samples) == 24,
              "mio_project_to_cad_opts layout");
static_assert(offsetof(mio_project_to_cad_opts, reserved) == 32, "mio_project_to_cad_opts layout");
static_assert(sizeof(mio_project_to_cad_report) == 96, "mio_project_to_cad_report layout");
static_assert(offsetof(mio_project_to_cad_report, max_move) == 56,
              "mio_project_to_cad_report layout");
static_assert(offsetof(mio_project_to_cad_report, num_on_curves) == 64,
              "mio_project_to_cad_report layout");

namespace {

std::vector<std::int64_t>& capi_project_to_cad_outside() {
    thread_local std::vector<std::int64_t> outside;
    return outside;
}

}  // namespace

void mio_project_to_cad_opts_init(mio_project_to_cad_opts* opts) {
    if (!opts)
        return;
    *opts = mio_project_to_cad_opts{};
    opts->trim_samples = 16;
}

mio_mesh* mio_project_to_cad(const mio_mesh* mesh, const mio_format_info* info,
                             const mio_project_to_cad_opts* opts,
                             mio_project_to_cad_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        capi_project_to_cad_outside().clear();
        if (!mesh)
            throw std::invalid_argument("meshio++: project_to_cad: mesh is NULL");
        if (!info)
            throw std::invalid_argument("meshio++: format info handle is NULL");
        const auto* p_brep = std::get_if<meshioplusplus::BrepInfo>(&info->mInfo);
        if (!p_brep)
            throw std::invalid_argument(
                "meshio++: project_to_cad: the info is not a B-rep (read a CAD JSON, STEP or "
                "IGES file)");
        meshioplusplus::ProjectToCadOptions o;
        if (opts) {
            if (opts->face_array)
                o.mFaceArray = opts->face_array;
            if (opts->edge_array)
                o.mEdgeArray = opts->edge_array;
            o.mMaxDistance = opts->max_distance;
            o.mTrimSamples = opts->trim_samples == 0 ? 16 : opts->trim_samples;
        }
        meshioplusplus::ProjectToCadResult r =
            meshioplusplus::project_to_cad(mesh->mMesh, *p_brep, o);
        if (report) {
            *report = mio_project_to_cad_report{};
            report->num_projected = r.mNumProjected;
            report->num_on_edges = r.mNumOnEdges;
            report->num_corners = r.mNumCorners;
            report->num_ambiguous = r.mNumAmbiguous;
            report->num_too_far = r.mNumTooFar;
            report->num_failed = r.mNumFailed;
            report->num_outside_trim = r.mNumOutsideTrim;
            report->max_move = r.mMaxMove;
            report->num_on_curves = r.mNumOnCurves;
        }
        capi_project_to_cad_outside() = std::move(r.mOutsideTrim);
        return new mio_mesh{std::move(r.mMesh)};
    });
}

static_assert(sizeof(mio_error_size_opts) == 80, "mio_error_size_opts layout");
static_assert(offsetof(mio_error_size_opts, output) == 40, "mio_error_size_opts layout");
static_assert(offsetof(mio_error_size_opts, reserved) == 56, "mio_error_size_opts layout");
static_assert(sizeof(mio_error_size_report) == 64, "mio_error_size_report layout");
static_assert(offsetof(mio_error_size_report, max_size) == 32, "mio_error_size_report layout");

void mio_error_size_opts_init(mio_error_size_opts* opts) {
    if (!opts)
        return;
    *opts = mio_error_size_opts{};
    opts->order = 1.0;
}

mio_mesh* mio_error_size(const mio_mesh* mesh, const mio_error_size_opts* opts,
                         mio_error_size_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw std::invalid_argument("meshio++: error_size: mesh is NULL");
        if (!opts)
            throw std::invalid_argument(
                "meshio++: error_size: opts is NULL (a target_error is "
                "required)");
        meshioplusplus::ErrorSizeOptions o;
        if (opts->error_array)
            o.mErrorArray = opts->error_array;
        o.mTargetError = opts->target_error;
        o.mOrder = opts->order;
        o.mMinSize = opts->min_size;
        o.mMaxSize = opts->max_size;
        if (opts->output)
            o.mOutputName = opts->output;
        o.mOverwrite = opts->overwrite != 0;
        meshioplusplus::ErrorSizeResult r = meshioplusplus::error_size(mesh->mMesh, o);
        if (report) {
            *report = mio_error_size_report{};
            report->num_cells = r.mNumCells;
            report->num_skipped = r.mNumSkipped;
            report->num_clamped = r.mNumClamped;
            report->min_size = r.mMinSize;
            report->max_size = r.mMaxSize;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

static_assert(sizeof(mio_conformity_opts) == 48, "mio_conformity_opts layout");
static_assert(offsetof(mio_conformity_opts, reserved) == 16, "mio_conformity_opts layout");
static_assert(sizeof(mio_conformity_report) == 104, "mio_conformity_report layout");
static_assert(offsetof(mio_conformity_report, conforming) == 56, "mio_conformity_report layout");
static_assert(offsetof(mio_conformity_report, num_warnings) == 64, "mio_conformity_report layout");

}  // extern "C"

namespace {

/// The warnings of this thread's last mio_check_conformity call.
std::string& capi_conformity_warnings() {
    thread_local std::string warnings;
    return warnings;
}

}  // namespace

extern "C" {

void mio_conformity_opts_init(mio_conformity_opts* opts) {
    if (!opts)
        return;
    *opts = mio_conformity_opts{};
    opts->tolerance = 1e-9;
    opts->attach_flags = 1;
}

mio_mesh* mio_check_conformity(const mio_mesh* mesh, const mio_conformity_opts* opts,
                               mio_conformity_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        capi_conformity_warnings().clear();
        if (!mesh)
            throw std::invalid_argument("meshio++: check_conformity: mesh is NULL");
        mio_conformity_opts defaults;
        mio_conformity_opts_init(&defaults);
        const mio_conformity_opts& rOpts = opts ? *opts : defaults;
        meshioplusplus::ConformityOptions o;
        o.mTolerance = rOpts.tolerance;
        o.mAttachFlags = rOpts.attach_flags != 0;
        meshioplusplus::ConformityResult r = meshioplusplus::check_conformity(mesh->mMesh, o);
        std::string& warnings = capi_conformity_warnings();
        for (const std::string& w : r.mWarnings) {
            if (!warnings.empty())
                warnings += '\n';
            warnings += w;
        }
        if (report) {
            *report = mio_conformity_report{};
            report->dimension = r.mDimension;
            report->num_hanging_nodes = r.mNumHangingNodes;
            report->num_t_junctions = r.mNumTJunctions;
            report->num_partial_overlaps = r.mNumPartialOverlaps;
            report->num_unwelded_facets = r.mNumUnweldedFacets;
            report->num_boundary_facets = r.mNumBoundaryFacets;
            report->num_non_manifold_facets = r.mNumNonManifoldFacets;
            report->conforming = r.mConforming ? 1 : 0;
            report->interface_mode = r.mInterfaceMode == "weld" ? MIO_CONFORMITY_WELD
                                                                : (r.mInterfaceMode == "proximity"
                                                                       ? MIO_CONFORMITY_PROXIMITY
                                                                       : MIO_CONFORMITY_CONFORMING);
            report->num_warnings = static_cast<int64_t>(r.mWarnings.size());
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

int64_t mio_check_conformity_warnings(char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() {
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative warnings buffer length");
        return copy_string(capi_conformity_warnings(), buf, buflen);
    });
}

static_assert(sizeof(mio_boundary_layer_opts) == 64, "mio_boundary_layer_opts layout");
static_assert(offsetof(mio_boundary_layer_opts, num_layers) == 24,
              "mio_boundary_layer_opts layout");
static_assert(sizeof(mio_boundary_layer_report) == 128, "mio_boundary_layer_report layout");
static_assert(offsetof(mio_boundary_layer_report, cap_closed) == 88,
              "mio_boundary_layer_report layout");
static_assert(offsetof(mio_boundary_layer_opts, collapse) == 40, "mio_boundary_layer_opts layout");
static_assert(offsetof(mio_boundary_layer_opts, reserved) == 56, "mio_boundary_layer_opts layout");
static_assert(offsetof(mio_boundary_layer_report, num_fewer_layers) == 96,
              "mio_boundary_layer_report layout");
static_assert(offsetof(mio_boundary_layer_report, reserved) == 112,
              "mio_boundary_layer_report layout");

void mio_boundary_layer_opts_init(mio_boundary_layer_opts* opts) {
    if (!opts)
        return;
    *opts = mio_boundary_layer_opts{};
    const meshioplusplus::BoundaryLayerOptions d;
    opts->growth = d.mGrowth;
    opts->collision_fraction = d.mCollisionFraction;
    opts->num_layers = d.mNumLayers;
    opts->smooth_iterations = d.mSmoothIterations;
    opts->max_shrink = d.mMaxShrink;
    opts->collapse = static_cast<int32_t>(d.mCollapse);
    opts->normal_smooth_iterations = d.mNormalSmoothIterations;
    opts->front_check = d.mFrontCheck ? 1 : 0;
}

mio_mesh* mio_boundary_layer(const mio_mesh* mesh, const mio_boundary_layer_opts* opts,
                             mio_boundary_layer_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw std::invalid_argument("meshio++: boundary_layer: mesh is NULL");
        if (!opts)
            throw std::invalid_argument("meshio++: boundary_layer: opts is NULL (first_height)");
        meshioplusplus::BoundaryLayerOptions o;
        o.mFirstHeight = opts->first_height;
        o.mGrowth = opts->growth;
        o.mCollisionFraction = opts->collision_fraction;
        o.mNumLayers = opts->num_layers;
        o.mReverse = opts->reverse != 0;
        o.mSmoothIterations = opts->smooth_iterations;
        o.mMaxShrink = opts->max_shrink;
        // Out-of-range values reach the core, which refuses them by name.
        o.mCollapse = static_cast<meshioplusplus::BoundaryLayerCollapse>(opts->collapse);
        o.mNormalSmoothIterations = opts->normal_smooth_iterations;
        o.mFrontCheck = opts->front_check != 0;
        meshioplusplus::BoundaryLayerResult r = meshioplusplus::boundary_layer(mesh->mMesh, o);
        if (report) {
            *report = mio_boundary_layer_report{};
            report->num_cells = r.mNumCells;
            report->num_wedges = r.mNumWedges;
            report->num_hexahedra = r.mNumHexahedra;
            report->num_pyramids = r.mNumPyramids;
            report->num_tetrahedra = r.mNumTetrahedra;
            report->num_shrunk = r.mNumShrunk;
            report->num_collapsed = r.mNumCollapsed;
            report->num_not_visible = r.mNumNotVisible;
            report->min_scaled_jacobian = r.mMinScaledJacobian;
            report->min_thickness = r.mMinThickness;
            report->max_thickness = r.mMaxThickness;
            report->cap_closed = r.mCapClosed ? 1 : 0;
            report->num_fewer_layers = r.mNumFewerLayers;
            report->num_front_crossings = r.mNumFrontCrossings;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

static_assert(sizeof(mio_extrude_opts) == 96, "mio_extrude_opts layout");
static_assert(offsetof(mio_extrude_opts, offsets) == 32, "mio_extrude_opts layout");
static_assert(offsetof(mio_extrude_opts, num_layers) == 56, "mio_extrude_opts layout");
static_assert(offsetof(mio_extrude_opts, reserved) == 64, "mio_extrude_opts layout");
static_assert(sizeof(mio_revolve_opts) == 112, "mio_revolve_opts layout");
static_assert(offsetof(mio_revolve_opts, angle) == 48, "mio_revolve_opts layout");
static_assert(offsetof(mio_revolve_opts, num_layers) == 72, "mio_revolve_opts layout");
static_assert(offsetof(mio_revolve_opts, reserved) == 80, "mio_revolve_opts layout");
static_assert(sizeof(mio_sweep_report) == 128, "mio_sweep_report layout");
static_assert(offsetof(mio_sweep_report, min_scaled_jacobian) == 80, "mio_sweep_report layout");
static_assert(offsetof(mio_sweep_report, full) == 88, "mio_sweep_report layout");
static_assert(offsetof(mio_sweep_report, reserved) == 96, "mio_sweep_report layout");
static_assert(sizeof(mio_mirror_opts) == 104, "mio_mirror_opts layout");
static_assert(offsetof(mio_mirror_opts, keep_original) == 56, "mio_mirror_opts layout");
static_assert(offsetof(mio_mirror_opts, reserved) == 72, "mio_mirror_opts layout");
static_assert(sizeof(mio_mirror_report) == 72, "mio_mirror_report layout");
static_assert(offsetof(mio_mirror_report, num_warnings) == 32, "mio_mirror_report layout");

}  // extern "C"

namespace {

/// The warnings of this thread's last mio_extrude, mio_revolve and mio_mirror calls.
std::string& capi_extrude_warnings() {
    thread_local std::string warnings;
    return warnings;
}
std::string& capi_revolve_warnings() {
    thread_local std::string warnings;
    return warnings;
}
std::string& capi_mirror_warnings() {
    thread_local std::string warnings;
    return warnings;
}

void capi_join_warnings(const std::vector<std::string>& rFrom, std::string& rInto) {
    rInto.clear();
    for (const std::string& w : rFrom) {
        if (!rInto.empty())
            rInto += '\n';
        rInto += w;
    }
}

void capi_sweep_report(const meshioplusplus::SweepResult& rR, mio_sweep_report* pReport) {
    if (!pReport)
        return;
    *pReport = mio_sweep_report{};
    pReport->num_cells = rR.mNumCells;
    pReport->num_hexahedra = rR.mNumHexahedra;
    pReport->num_wedges = rR.mNumWedges;
    pReport->num_pyramids = rR.mNumPyramids;
    pReport->num_tetrahedra = rR.mNumTetrahedra;
    pReport->num_quads = rR.mNumQuads;
    pReport->num_triangles = rR.mNumTriangles;
    pReport->num_flipped = rR.mNumFlipped;
    pReport->num_on_axis = rR.mNumOnAxis;
    pReport->num_warnings = static_cast<int64_t>(rR.mWarnings.size());
    pReport->min_scaled_jacobian = rR.mMinScaledJacobian;
    pReport->full = rR.mFull ? 1 : 0;
}

}  // namespace

extern "C" {

void mio_extrude_opts_init(mio_extrude_opts* opts) {
    if (!opts)
        return;
    *opts = mio_extrude_opts{};
    const meshioplusplus::ExtrudeOptions d;
    for (int i = 0; i < 3; ++i)
        opts->direction[i] = d.mDirection[static_cast<std::size_t>(i)];
    opts->ratio = d.mRatio;
    opts->max_cells = d.mMaxCells;
    opts->num_layers = d.mNumLayers;
}

void mio_revolve_opts_init(mio_revolve_opts* opts) {
    if (!opts)
        return;
    *opts = mio_revolve_opts{};
    const meshioplusplus::RevolveOptions d;
    for (int i = 0; i < 3; ++i) {
        opts->axis_point[i] = d.mAxisPoint[static_cast<std::size_t>(i)];
        opts->axis[i] = d.mAxis[static_cast<std::size_t>(i)];
    }
    opts->angle = d.mAngle;
    opts->tolerance = d.mTolerance;
    opts->max_cells = d.mMaxCells;
    opts->num_layers = d.mNumLayers;
}

mio_mesh* mio_extrude(const mio_mesh* mesh, const mio_extrude_opts* opts,
                      mio_sweep_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        capi_extrude_warnings().clear();
        if (!mesh)
            throw std::invalid_argument("meshio++: extrude: mesh is NULL");
        mio_extrude_opts defaults;
        mio_extrude_opts_init(&defaults);
        const mio_extrude_opts& rOpts = opts ? *opts : defaults;
        if (rOpts.num_offsets < 0 || (rOpts.num_offsets > 0 && !rOpts.offsets))
            throw std::invalid_argument("meshio++: extrude: offsets is NULL or its count negative");
        meshioplusplus::ExtrudeOptions o;
        for (int i = 0; i < 3; ++i)
            o.mDirection[static_cast<std::size_t>(i)] = rOpts.direction[i];
        o.mRatio = rOpts.ratio;
        o.mOffsets.assign(rOpts.offsets, rOpts.offsets + rOpts.num_offsets);
        o.mMaxCells = rOpts.max_cells;
        o.mNumLayers = rOpts.num_layers;
        o.mRecordLayers = rOpts.record_layers != 0;
        meshioplusplus::SweepResult r = meshioplusplus::extrude(mesh->mMesh, o);
        capi_join_warnings(r.mWarnings, capi_extrude_warnings());
        capi_sweep_report(r, report);
        return new mio_mesh{std::move(r.mMesh)};
    });
}

mio_mesh* mio_revolve(const mio_mesh* mesh, const mio_revolve_opts* opts,
                      mio_sweep_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        capi_revolve_warnings().clear();
        if (!mesh)
            throw std::invalid_argument("meshio++: revolve: mesh is NULL");
        mio_revolve_opts defaults;
        mio_revolve_opts_init(&defaults);
        const mio_revolve_opts& rOpts = opts ? *opts : defaults;
        meshioplusplus::RevolveOptions o;
        for (int i = 0; i < 3; ++i) {
            o.mAxisPoint[static_cast<std::size_t>(i)] = rOpts.axis_point[i];
            o.mAxis[static_cast<std::size_t>(i)] = rOpts.axis[i];
        }
        o.mAngle = rOpts.angle;
        o.mTolerance = rOpts.tolerance;
        o.mMaxCells = rOpts.max_cells;
        o.mNumLayers = rOpts.num_layers;
        o.mRecordLayers = rOpts.record_layers != 0;
        meshioplusplus::SweepResult r = meshioplusplus::revolve(mesh->mMesh, o);
        capi_join_warnings(r.mWarnings, capi_revolve_warnings());
        capi_sweep_report(r, report);
        return new mio_mesh{std::move(r.mMesh)};
    });
}

int64_t mio_extrude_warnings(char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() {
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative warnings buffer length");
        return copy_string(capi_extrude_warnings(), buf, buflen);
    });
}

int64_t mio_revolve_warnings(char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() {
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative warnings buffer length");
        return copy_string(capi_revolve_warnings(), buf, buflen);
    });
}

void mio_mirror_opts_init(mio_mirror_opts* opts) {
    if (!opts)
        return;
    *opts = mio_mirror_opts{};
    const meshioplusplus::MirrorOptions d;
    for (int i = 0; i < 3; ++i) {
        opts->origin[i] = d.mOrigin[static_cast<std::size_t>(i)];
        opts->normal[i] = d.mNormal[static_cast<std::size_t>(i)];
    }
    opts->tolerance = d.mTolerance;
    opts->weld = d.mWeld ? 1 : 0;
}

mio_mesh* mio_mirror(const mio_mesh* mesh, const mio_mirror_opts* opts, mio_mirror_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        capi_mirror_warnings().clear();
        if (!mesh)
            throw std::invalid_argument("meshio++: mirror: mesh is NULL");
        mio_mirror_opts defaults;
        mio_mirror_opts_init(&defaults);
        const mio_mirror_opts& rOpts = opts ? *opts : defaults;
        meshioplusplus::MirrorOptions o;
        for (int i = 0; i < 3; ++i) {
            o.mOrigin[static_cast<std::size_t>(i)] = rOpts.origin[i];
            o.mNormal[static_cast<std::size_t>(i)] = rOpts.normal[i];
        }
        o.mTolerance = rOpts.tolerance;
        o.mKeepOriginal = rOpts.keep_original != 0;
        o.mWeld = rOpts.weld != 0;
        o.mReflectVectorData = rOpts.reflect_vector_data != 0;
        meshioplusplus::MirrorResult r = meshioplusplus::mirror(mesh->mMesh, o);
        capi_join_warnings(r.mWarnings, capi_mirror_warnings());
        if (report) {
            *report = mio_mirror_report{};
            report->num_flipped_cells = r.mNumFlippedCells;
            report->num_points_on_plane = r.mNumPointsOnPlane;
            report->num_dropped_cells = r.mNumDroppedCells;
            report->num_dropped_sides = r.mNumDroppedSides;
            report->num_warnings = static_cast<int64_t>(r.mWarnings.size());
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

int64_t mio_mirror_warnings(char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() {
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative warnings buffer length");
        return copy_string(capi_mirror_warnings(), buf, buflen);
    });
}

static_assert(sizeof(mio_tetrahedralize_opts) == 80, "mio_tetrahedralize_opts layout");
static_assert(offsetof(mio_tetrahedralize_opts, face_array) == 32,
              "mio_tetrahedralize_opts layout");
static_assert(offsetof(mio_tetrahedralize_opts, backend) == 40, "mio_tetrahedralize_opts layout");
static_assert(offsetof(mio_tetrahedralize_opts, reserved) == 48, "mio_tetrahedralize_opts layout");
static_assert(sizeof(mio_tetrahedralize_report) == 128, "mio_tetrahedralize_report layout");
static_assert(offsetof(mio_tetrahedralize_report, boundary_preserved) == 72,
              "mio_tetrahedralize_report layout");
static_assert(offsetof(mio_tetrahedralize_report, reserved) == 80,
              "mio_tetrahedralize_report layout");
static_assert(static_cast<int>(meshioplusplus::TetrahedralizeBackend::Netgen) ==
                      MIO_TETRAHEDRALIZE_BACKEND_NETGEN &&
                  static_cast<int>(meshioplusplus::TetrahedralizeBackend::Octree) ==
                      MIO_TETRAHEDRALIZE_BACKEND_OCTREE,
              "mio_tetrahedralize_backend drifted from TetrahedralizeBackend");

void mio_tetrahedralize_opts_init(mio_tetrahedralize_opts* opts) {
    if (!opts)
        return;
    *opts = mio_tetrahedralize_opts{};
    const meshioplusplus::TetrahedralizeOptions d;
    opts->max_size = d.mMaxSize;
    opts->max_volume = d.mMaxVolume;
    opts->grading = d.mGrading;
    opts->max_cells = d.mMaxCells;
    opts->backend = static_cast<int32_t>(d.mBackend);
    opts->optimize_steps = d.mOptimizeSteps;
}

int mio_has_netgen(void) {
    return meshioplusplus::tetrahedralize_has_netgen() ? 1 : 0;
}

mio_mesh* mio_tetrahedralize(const mio_mesh* mesh, const mio_tetrahedralize_opts* opts,
                             mio_tetrahedralize_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw std::invalid_argument("meshio++: tetrahedralize: mesh is NULL");
        mio_tetrahedralize_opts defaults;
        mio_tetrahedralize_opts_init(&defaults);
        const mio_tetrahedralize_opts& r_in = opts ? *opts : defaults;
        meshioplusplus::TetrahedralizeOptions o;
        // Out-of-range values reach the core, which refuses them by name.
        o.mBackend = static_cast<meshioplusplus::TetrahedralizeBackend>(r_in.backend);
        o.mMaxSize = r_in.max_size;
        o.mMaxVolume = r_in.max_volume;
        o.mGrading = r_in.grading;
        o.mOptimizeSteps = r_in.optimize_steps;
        o.mMaxCells = r_in.max_cells;
        if (r_in.face_array)
            o.mFaceArray = r_in.face_array;
        meshioplusplus::TetrahedralizeResult r = meshioplusplus::tetrahedralize(mesh->mMesh, o);
        if (report) {
            *report = mio_tetrahedralize_report{};
            report->num_tetrahedra = r.mNumTetrahedra;
            report->num_points = r.mNumPoints;
            report->num_new_points = r.mNumNewPoints;
            report->num_regions = r.mNumRegions;
            report->num_boundary_triangles = r.mNumBoundaryTriangles;
            report->hausdorff_to_surface = r.mHausdorffToSurface;
            report->min_dihedral_deg = r.mMinDihedralDeg;
            report->max_dihedral_deg = r.mMaxDihedralDeg;
            report->volume = r.mVolume;
            report->boundary_preserved = r.mBoundaryPreserved ? 1 : 0;
            report->backend = r.mBackend == "netgen" ? MIO_TETRAHEDRALIZE_BACKEND_NETGEN
                                                     : MIO_TETRAHEDRALIZE_BACKEND_OCTREE;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

int64_t mio_project_to_cad_outside_trim(int64_t* points, int64_t n) {
    return guarded_ptr(int64_t{-1}, [&]() -> int64_t {
        if (n < 0 || (n > 0 && !points))
            throw std::invalid_argument("meshio++: project_to_cad: bad points buffer");
        const auto& r_o = capi_project_to_cad_outside();
        const auto count = std::min<std::int64_t>(n, static_cast<std::int64_t>(r_o.size()));
        std::copy(r_o.begin(), r_o.begin() + count, points);
        return static_cast<int64_t>(r_o.size());
    });
}

int64_t mio_mesh_cad_surface_num_untangled(void) {
    return capi_mesh_cad_state().mExtras.mNumUntangled;
}

mio_status mio_mesh_cad_surface_quad_stats(int64_t* num_quads, int64_t* num_leftover_triangles,
                                           double* min_scaled_jacobian) {
    return guarded([&]() -> mio_status {
        const auto& r_x = capi_mesh_cad_state().mExtras;
        if (num_quads)
            *num_quads = r_x.mNumQuads;
        if (num_leftover_triangles)
            *num_leftover_triangles = r_x.mNumLeftoverTriangles;
        if (min_scaled_jacobian)
            *min_scaled_jacobian = r_x.mMinQuadScaledJacobian;
        return MIO_OK;
    });
}

double mio_mesh_cad_surface_min_metric_angle(void) {
    return capi_mesh_cad_state().mExtras.mMinMetricAngleDeg;
}

int64_t mio_mesh_cad_surface_num_lopped_corners(void) {
    return static_cast<int64_t>(capi_mesh_cad_state().mExtras.mLoppedCorners.size());
}

mio_status mio_mesh_cad_surface_lopped_corner(int64_t index, int64_t* face, int64_t* point,
                                              double* angle_deg, double* lfs0, double* radius,
                                              double* turning_deg, double* stated_bound_deg,
                                              int64_t* num_forced, double* min_forced_deg,
                                              int32_t* within_two_radii, int32_t* guaranteed) {
    return guarded([&]() -> mio_status {
        const auto& r_v = capi_mesh_cad_state().mExtras.mLoppedCorners;
        if (index < 0 || static_cast<std::size_t>(index) >= r_v.size())
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: mesh_cad_surface: lopped corner index out of range");
        const auto& r_c = r_v[static_cast<std::size_t>(index)];
        if (face)
            *face = r_c.mFace;
        if (point)
            *point = r_c.mPoint;
        if (angle_deg)
            *angle_deg = r_c.mAngleDeg;
        if (lfs0)
            *lfs0 = r_c.mLfs0;
        if (radius)
            *radius = r_c.mLopRadius;
        if (turning_deg)
            *turning_deg = r_c.mTurningDeg;
        if (stated_bound_deg)
            *stated_bound_deg = r_c.mStatedBoundDeg;
        if (num_forced)
            *num_forced = r_c.mNumForced;
        if (min_forced_deg)
            *min_forced_deg = r_c.mMinForcedDeg;
        if (within_two_radii)
            *within_two_radii = r_c.mWithinTwoRadii ? 1 : 0;
        if (guaranteed)
            *guaranteed = r_c.mGuaranteed ? 1 : 0;
        return MIO_OK;
    });
}

int64_t mio_mesh_cad_surface_lopped_corner_reason(int64_t index, char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() {
        const auto& r_v = capi_mesh_cad_state().mExtras.mLoppedCorners;
        if (index < 0 || static_cast<std::size_t>(index) >= r_v.size())
            throw std::invalid_argument(
                "meshio++: mesh_cad_surface: lopped corner index out of range");
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative reason buffer length");
        return copy_string(r_v[static_cast<std::size_t>(index)].mReason, buf, buflen);
    });
}

int64_t mio_mesh_cad_surface_num_corners(void) {
    return static_cast<int64_t>(capi_mesh_cad_state().mCorners.size());
}

mio_status mio_mesh_cad_surface_corner(int64_t index, mio_mesh_cad_corner* corner) {
    return guarded([&]() -> mio_status {
        const auto& corners = capi_mesh_cad_state().mCorners;
        if (!corner)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh_cad_surface: corner is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= corners.size())
            return fail(MIO_ERR_INVALID_ARG,
                        "meshio++: mesh_cad_surface: corner index out of range");
        const meshioplusplus::MeshCadCorner& c = corners[static_cast<std::size_t>(index)];
        *corner = mio_mesh_cad_corner{};
        corner->face = c.mFace;
        corner->point = c.mPoint;
        corner->u = c.mU;
        corner->v = c.mV;
        for (int k = 0; k < 3; ++k)
            corner->x[k] = c.mX[static_cast<std::size_t>(k)];
        corner->angle_deg = c.mAngleDeg;
        corner->min_forced_deg = c.mMinForcedDeg;
        corner->reach = c.mReach;
        corner->pole = c.mPole ? 1 : 0;
        return MIO_OK;
    });
}

int mio_has_occt(void) {
    return meshioplusplus::fit_cad_has_occt() ? 1 : 0;
}

mio_mesh* mio_fit_cad(const mio_mesh* mesh, const char* step_path, const mio_fit_cad_opts* opts,
                      mio_fit_cad_report* report) {
    return mio_fit_cad2(mesh, step_path, opts, nullptr, report, nullptr, nullptr);
}

mio_mesh* mio_fit_cad_ex(const mio_mesh* mesh, const char* step_path, const mio_fit_cad_opts* opts,
                         mio_fit_cad_report* report, mio_fit_cad_report_ex* report_ex) {
    return mio_fit_cad2(mesh, step_path, opts, nullptr, report, report_ex, nullptr);
}

mio_mesh* mio_fit_cad2(const mio_mesh* mesh, const char* step_path, const mio_fit_cad_opts* opts,
                       const mio_fit_cad_opts_ex* opts_ex, mio_fit_cad_report* report,
                       mio_fit_cad_report_ex* report_ex, mio_fit_cad_report_ex2* report_ex2) {
    return mio_fit_cad3(mesh, step_path, opts, opts_ex, nullptr, report, report_ex, report_ex2,
                        nullptr);
}

mio_mesh* mio_fit_cad3(const mio_mesh* mesh, const char* step_path, const mio_fit_cad_opts* opts,
                       const mio_fit_cad_opts_ex* opts_ex, const mio_fit_cad_opts_ex2* opts_ex2,
                       mio_fit_cad_report* report, mio_fit_cad_report_ex* report_ex,
                       mio_fit_cad_report_ex2* report_ex2, mio_fit_cad_report_ex3* report_ex3) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::FitCadOptions options;
        if (opts) {
            options.mTolerance = opts->tolerance;
            options.mMeshDeflection = opts->mesh_deflection;
            options.mMaxDemotionRounds = opts->max_demotion_rounds;
            options.mExactEdges = opts->exact_edges != 0;
            options.mTypes = opts->types;
            // A struct zero-filled before these fields existed keeps freeform,
            // splitting, fillets and threads off and the defaults for the rest.
            options.mFreeform = opts->freeform != 0;
            options.mFreeformDegree = opts->freeform_degree ? opts->freeform_degree : 3;
            options.mFreeformMaxControl =
                opts->freeform_max_control ? opts->freeform_max_control : 64;
            options.mFreeformSplit = opts->freeform_split != 0;
            options.mFillets = opts->fillets != 0;
            options.mThreads = opts->threads != 0;
            options.mMergeFaces = opts->merge_faces != 0;
            if (opts->backend != MIO_FIT_CAD_AUTO && opts->backend != MIO_FIT_CAD_OCCT &&
                opts->backend != MIO_FIT_CAD_NATIVE) {
                fail(MIO_ERR_INVALID_ARG,
                     "meshio++: fit_cad: backend must be MIO_FIT_CAD_AUTO, "
                     "MIO_FIT_CAD_OCCT or MIO_FIT_CAD_NATIVE");
                return nullptr;
            }
            options.mBackend = static_cast<meshioplusplus::FitCadBackend>(opts->backend);
        }
        if (opts_ex) {
            options.mRefineSegmentation = opts_ex->refine_segmentation != 0;
            options.mRefine = capi_refine_patches_options(opts_ex->refine);
            options.mFreeformGeneral = opts_ex->freeform_general != 0;
            options.mFreeformTolerance = opts_ex->freeform_tolerance;
            options.mFreeformSmoothing = opts_ex->freeform_smoothing;
            // Zero, from a struct filled before these fields existed, reads as
            // the default; any other value is taken as given.
            if (opts_ex->max_freeform_poles != 0)
                options.mMaxFreeformPoles = opts_ex->max_freeform_poles;
            if (opts_ex->max_freeform_iterations != 0)
                options.mMaxFreeformIterations = opts_ex->max_freeform_iterations;
            options.mMaxFreeformSplits = opts_ex->max_freeform_splits;
            options.mSewTolerance = opts_ex->sew_tolerance;
            options.mExtendSurfaces = opts_ex->extend_surfaces != 0;
            options.mSnapJunctions = opts_ex->snap_junctions != 0;
            options.mFaceDeviation = opts_ex->face_deviation != 0;
        }
        if (opts_ex2) {
            // Enumerations out of range are refused by fit_cad, naming them.
            options.mFreeformMode =
                static_cast<meshioplusplus::FitCadFreeformMode>(opts_ex2->freeform_mode);
            options.mFreeformParameterization = static_cast<meshioplusplus::FitCadParameterization>(
                opts_ex2->freeform_parameterization);
            options.mExtraordinaryVertex = static_cast<meshioplusplus::FitCadExtraordinaryVertex>(
                opts_ex2->extraordinary_vertex);
            options.mLocalRefine =
                static_cast<meshioplusplus::FitCadLocalRefine>(opts_ex2->local_refine);
            // Zero reads as the default for the tolerances, budgets and sizes.
            if (opts_ex2->freeform_g1_tolerance_deg != 0.0)
                options.mFreeformG1ToleranceDeg = opts_ex2->freeform_g1_tolerance_deg;
            if (opts_ex2->freeform_max_control_points != 0)
                options.mFreeformMaxControlPoints = opts_ex2->freeform_max_control_points;
            if (opts_ex2->layout_target_patch_size != 0.0)
                options.mLayoutTargetPatchSize = opts_ex2->layout_target_patch_size;
            if (opts_ex2->layout_max_valence_deviation != 0)
                options.mLayoutMaxValenceDeviation = opts_ex2->layout_max_valence_deviation;
            options.mResegmentRounds = opts_ex2->resegment_rounds;
            options.mLayoutSnapToFeatures = opts_ex2->layout_snap_to_features != 0;
            options.mValidateBySurfaceSampling = opts_ex2->validate_by_surface_sampling != 0;
            options.mConstrainPrimitives = opts_ex2->constrain_primitives != 0;
        }
        meshioplusplus::FitCadResult r =
            meshioplusplus::fit_cad(mesh->mMesh, step_path ? step_path : "", options);
        if (report) {
            *report = mio_fit_cad_report{};
            report->tolerance = r.mTolerance;
            report->volume = r.mVolume;
            report->skin_volume = r.mSkinVolume;
            report->volume_ratio = r.mVolumeRatio;
            report->max_deviation = r.mMaxDeviation;
            report->num_faces = r.mNumFaces;
            report->num_plane = r.mNumPlane;
            report->num_sphere = r.mNumSphere;
            report->num_cylinder = r.mNumCylinder;
            report->num_cone = r.mNumCone;
            report->num_torus = r.mNumTorus;
            report->num_faceted_faces = r.mNumFacetedFaces;
            report->num_patches = r.mNumPatches;
            report->num_demoted = r.mNumDemoted;
            report->num_exact_edges = r.mNumExactEdges;
            report->num_approx_edges = r.mNumApproxEdges;
            report->num_solids = r.mNumSolids;
            report->valid = r.mValid ? 1 : 0;
            report->num_bspline = r.mNumBspline;
            report->num_split_patches = r.mNumSplitPatches;
            report->num_fillet = r.mNumFillet;
            report->num_helical = r.mNumHelical;
            report->num_variable_fillet = r.mNumVariableFillet;
            report->num_blend = r.mNumBlend;
            report->num_rounded_thread =
                static_cast<std::int32_t>(std::min<std::int64_t>(r.mNumRoundedThread, INT32_MAX));
        }
        if (report_ex) {
            *report_ex = mio_fit_cad_report_ex{};
            report_ex->num_vertex_blend = r.mNumVertexBlend;
            report_ex->num_chamfer = r.mNumChamfer;
            report_ex->num_conic_blend = r.mNumConicBlend;
            report_ex->num_tapered_thread = r.mNumTaperedThread;
            report_ex->num_variable_pitch_thread = r.mNumVariablePitchThread;
            report_ex->backend = r.mBackend == meshioplusplus::FitCadBackend::Native
                                     ? MIO_FIT_CAD_NATIVE
                                     : MIO_FIT_CAD_OCCT;
        }
        if (report_ex2) {
            *report_ex2 = mio_fit_cad_report_ex2{};
            report_ex2->refined = r.mRefined ? 1 : 0;
            report_ex2->refine_iterations = r.mRefineIterations;
            report_ex2->num_refine_splits = r.mNumRefineSplits;
            report_ex2->num_refine_merges = r.mNumRefineMerges;
            report_ex2->initial_explained_area_fraction = r.mInitialExplainedAreaFraction;
            report_ex2->explained_area_fraction = r.mExplainedAreaFraction;
            report_ex2->num_freeform_general = r.mNumFreeformGeneral;
            report_ex2->num_freeform_refinements = r.mNumFreeformRefinements;
            report_ex2->avg_freeform_poles = r.mAvgFreeformPoles;
            report_ex2->num_extended_edges = r.mNumExtendedEdges;
            report_ex2->num_snapped_junctions = r.mNumSnappedJunctions;
            report_ex2->hausdorff_max = r.mHausdorffMax;
            report_ex2->hausdorff_rms = r.mHausdorffRms;
            report_ex2->closed = r.mClosed ? 1 : 0;
            report_ex2->refinement_kept = r.mRefinementKept ? 1 : 0;
        }
        if (report_ex3) {
            *report_ex3 = mio_fit_cad_report_ex3{};
            report_ex3->brep_checked = r.mBrepChecked ? 1 : 0;
            report_ex3->freeform_mode = static_cast<std::int32_t>(r.mFreeformMode);
            report_ex3->num_freeform_adaptive = r.mNumFreeformAdaptive;
            report_ex3->num_network_patches = r.mNumNetworkPatches;
            report_ex3->num_network_edges = r.mNumNetworkEdges;
            report_ex3->num_extraordinary_vertices = r.mNumExtraordinaryVertices;
            report_ex3->num_resegmented = r.mNumResegmented;
            report_ex3->num_escalations = r.mNumEscalations;
            report_ex3->facet_area_fraction = r.mFacetAreaFraction;
            report_ex3->max_g1_angle_deg = r.mMaxG1AngleDeg;
            report_ex3->max_freeform_distance = r.mMaxFreeformDistance;
            report_ex3->rms_freeform_distance = r.mRmsFreeformDistance;
            report_ex3->num_face_records = static_cast<std::int64_t>(r.mFaces.size());
            report_ex3->num_constrained_primitives = r.mNumConstrainedPrimitives;
            report_ex3->num_constraint_groups = r.mNumConstraintGroups;
        }
        return new mio_mesh{std::move(r.mMesh)};
    });
}

void mio_hausdorff_opts_init(mio_hausdorff_opts* opts) {
    if (opts)
        *opts = mio_hausdorff_opts{};
}

mio_status mio_hausdorff_distance(const mio_mesh* a, const mio_mesh* b,
                                  const mio_hausdorff_opts* opts, mio_hausdorff_report* report) {
    return guarded([&]() -> mio_status {
        if (!a || !b || !report)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: hausdorff_distance: a/b/report is NULL");
        meshioplusplus::HausdorffOptions options;
        if (opts) {
            options.mFaceSamples = opts->face_samples;
            options.mGridCellSize = opts->grid_cell_size;
            if (opts->region_a)
                options.mRegionA = opts->region_a;
            if (opts->region_b)
                options.mRegionB = opts->region_b;
        }
        const meshioplusplus::HausdorffResult r =
            meshioplusplus::hausdorff_distance(a->mMesh, b->mMesh, options);
        *report = mio_hausdorff_report{};
        report->distance = r.mDistance;
        report->a_to_b = r.mAtoB;
        report->b_to_a = r.mBtoA;
        report->mean_a_to_b = r.mMeanAtoB;
        report->rms_a_to_b = r.mRmsAtoB;
        report->mean_b_to_a = r.mMeanBtoA;
        report->rms_b_to_a = r.mRmsBtoA;
        report->num_samples_a = r.mNumSamplesA;
        report->num_samples_b = r.mNumSamplesB;
        for (int k = 0; k < 3; ++k) {
            report->worst_point_a[k] = r.mWorstPointA[static_cast<std::size_t>(k)];
            report->worst_point_b[k] = r.mWorstPointB[static_cast<std::size_t>(k)];
        }
        return MIO_OK;
    });
}

void mio_region_selector_init(mio_region_selector* sel, const char* name) {
    if (!sel)
        return;
    *sel = mio_region_selector{};
    sel->name = name;
    sel->kind = -1;
    sel->dim = MIO_REGION_ANY;
    sel->tag = MIO_REGION_ANY;
}

mio_mesh* mio_edit_regions(const mio_mesh* mesh, int32_t op, const mio_region_selector* inputs,
                           int64_t num_inputs, const char* output, int64_t dim, int64_t tag,
                           int32_t keep_inputs) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        if (op < MIO_REGION_UNION || op > MIO_REGION_DELETE)
            throw std::invalid_argument("meshio++: edit_regions: op is not a mio_region_op");
        if (num_inputs < 0 || (num_inputs > 0 && !inputs))
            throw std::invalid_argument("meshio++: edit_regions: inputs is NULL");
        meshioplusplus::RegionEdit e;
        e.mOp = static_cast<meshioplusplus::RegionOp>(op);
        for (int64_t i = 0; i < num_inputs; ++i)
            e.mInputs.push_back(capi_selector(inputs + i));
        if (output)
            e.mOutputName = output;
        e.mOutputDim = dim;
        e.mOutputTag = tag;
        e.mKeepInputs = keep_inputs != 0;
        return new mio_mesh{meshioplusplus::edit_regions(mesh->mMesh, {e})};
    });
}

mio_mesh* mio_region_adjacency(const mio_mesh* mesh, const mio_region_selector* regions,
                               int64_t num_regions) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        if (num_regions < 0 || (num_regions > 0 && !regions))
            throw std::invalid_argument(
                "meshio++: region_adjacency: invalid selectors pointer/count");
        std::vector<meshioplusplus::RegionSelector> selectors;
        selectors.reserve(static_cast<std::size_t>(num_regions));
        for (int64_t i = 0; i < num_regions; ++i)
            selectors.push_back(capi_selector(regions + i));
        return new mio_mesh{meshioplusplus::region_adjacency(mesh->mMesh, selectors)};
    });
}

void mio_find_interface_opts_init(mio_find_interface_opts* opts) {
    if (opts) {
        *opts = mio_find_interface_opts{};
        opts->angle_tolerance = 30.0;
    }
}

void mio_contact_pairs_opts_init(mio_contact_pairs_opts* opts) {
    if (opts)
        *opts = mio_contact_pairs_opts{};
}

void mio_split_interface_opts_init(mio_split_interface_opts* opts) {
    if (opts)
        *opts = mio_split_interface_opts{};
}

mio_find_interface_result* mio_find_interface(const mio_mesh* mesh_a,
                                              const mio_region_selector* region_a,
                                              const mio_mesh* mesh_b,
                                              const mio_region_selector* region_b,
                                              const mio_find_interface_opts* opts) {
    return guarded_ptr(static_cast<mio_find_interface_result*>(nullptr), [&]() {
        if (!mesh_a)
            throw meshioplusplus::ReadError("meshio++: find_interface: mesh_a is NULL");
        meshioplusplus::FindInterfaceOptions options;
        if (opts) {
            if (opts->mode != MIO_INTERFACE_CONFORMING && opts->mode != MIO_INTERFACE_PROXIMITY)
                throw std::invalid_argument("meshio++: find_interface: invalid mode");
            if (opts->master != MIO_INTERFACE_MASTER_A && opts->master != MIO_INTERFACE_MASTER_B)
                throw std::invalid_argument("meshio++: find_interface: invalid master side");
            options.mMode = static_cast<meshioplusplus::InterfaceMode>(opts->mode);
            options.mMaster = static_cast<meshioplusplus::InterfaceMaster>(opts->master);
            options.mGapTolerance = opts->gap_tolerance;
            options.mAngleTolerance = opts->angle_tolerance;
            options.mOverlapTolerance = opts->overlap_tolerance;
        }
        const meshioplusplus::RegionSelector selector_a = capi_selector(region_a);
        const meshioplusplus::RegionSelector selector_b = capi_selector(region_b);
        meshioplusplus::FindInterfaceResult result =
            mesh_b ? meshioplusplus::find_interface(mesh_a->mMesh, selector_a, mesh_b->mMesh,
                                                    selector_b, options)
                   : meshioplusplus::find_interface(mesh_a->mMesh, selector_a, selector_b, options);
        return new mio_find_interface_result{mio_mesh{std::move(result.mMesh)},
                                             std::move(result.mSideA), std::move(result.mSideB),
                                             result.mReport};
    });
}

const mio_mesh* mio_find_interface_result_mesh(const mio_find_interface_result* result) {
    return result ? &result->mMesh : nullptr;
}

mio_mesh* mio_find_interface_result_take_mesh(mio_find_interface_result* result) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!result)
            throw meshioplusplus::ReadError("meshio++: find_interface result is NULL");
        return new mio_mesh{std::move(result->mMesh.mMesh)};
    });
}

mio_status mio_find_interface_result_report(const mio_find_interface_result* result,
                                            int64_t* num_pairs, double* area, double* max_gap,
                                            int64_t* unmatched_a, int64_t* unmatched_b) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: find_interface result is NULL");
        if (num_pairs)
            *num_pairs = result->mReport.mNumPairs;
        if (area)
            *area = result->mReport.mArea;
        if (max_gap)
            *max_gap = result->mReport.mMaxGap;
        if (unmatched_a)
            *unmatched_a = result->mReport.mUnmatchedA;
        if (unmatched_b)
            *unmatched_b = result->mReport.mUnmatchedB;
        return MIO_OK;
    });
}

const int64_t* mio_find_interface_result_side_a(const mio_find_interface_result* result,
                                                int64_t* count) {
    if (count)
        *count = result ? static_cast<int64_t>(result->mSideA.mEntries.Shape()[0]) : 0;
    return result && result->mSideA.mEntries.Size() > 0 ? result->mSideA.mEntries.As<std::int64_t>()
                                                        : nullptr;
}

const int64_t* mio_find_interface_result_side_b(const mio_find_interface_result* result,
                                                int64_t* count) {
    if (count)
        *count = result ? static_cast<int64_t>(result->mSideB.mEntries.Shape()[0]) : 0;
    return result && result->mSideB.mEntries.Size() > 0 ? result->mSideB.mEntries.As<std::int64_t>()
                                                        : nullptr;
}

void mio_find_interface_result_free(mio_find_interface_result* result) {
    delete result;
}

mio_contact_pairs_result* mio_contact_pairs(const mio_mesh* slave_mesh,
                                            const mio_region_selector* slave_points,
                                            const mio_mesh* master_mesh,
                                            const mio_region_selector* master_cells,
                                            const mio_contact_pairs_opts* opts) {
    return guarded_ptr(static_cast<mio_contact_pairs_result*>(nullptr), [&]() {
        if (!slave_mesh)
            throw meshioplusplus::ReadError("meshio++: contact_pairs: slave mesh is NULL");
        meshioplusplus::ContactPairsOptions options;
        if (opts) {
            options.mTolerance = opts->tolerance;
            options.mRequireComplete = opts->require_complete != 0;
        }
        const meshioplusplus::RegionSelector point_selector = capi_selector(slave_points);
        const meshioplusplus::RegionSelector cell_selector = capi_selector(master_cells);
        const meshioplusplus::Mesh& master = master_mesh ? master_mesh->mMesh : slave_mesh->mMesh;
        return new mio_contact_pairs_result{meshioplusplus::contact_pairs(
            slave_mesh->mMesh, point_selector, master, cell_selector, options)};
    });
}

mio_status mio_contact_pairs_result_info(const mio_contact_pairs_result* result, int64_t* count,
                                         int64_t* unmatched_count) {
    return guarded([&]() -> mio_status {
        if (!result)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: contact_pairs result is NULL");
        if (count)
            *count = static_cast<int64_t>(result->mPairs.mSlavePoint.Size());
        if (unmatched_count)
            *unmatched_count = static_cast<int64_t>(result->mPairs.mUnmatched.Size());
        return MIO_OK;
    });
}

const int64_t* mio_contact_pairs_slave_point(const mio_contact_pairs_result* result) {
    return result && result->mPairs.mSlavePoint.Size() > 0
               ? result->mPairs.mSlavePoint.As<std::int64_t>()
               : nullptr;
}
const int64_t* mio_contact_pairs_master_cell(const mio_contact_pairs_result* result) {
    return result && result->mPairs.mMasterCell.Size() > 0
               ? result->mPairs.mMasterCell.As<std::int64_t>()
               : nullptr;
}
const int64_t* mio_contact_pairs_master_facet(const mio_contact_pairs_result* result) {
    return result && result->mPairs.mMasterFacet.Size() > 0
               ? result->mPairs.mMasterFacet.As<std::int64_t>()
               : nullptr;
}
const int64_t* mio_contact_pairs_master_subfacet(const mio_contact_pairs_result* result) {
    return result && result->mPairs.mMasterSubfacet.Size() > 0
               ? result->mPairs.mMasterSubfacet.As<std::int64_t>()
               : nullptr;
}
const double* mio_contact_pairs_local_coordinates(const mio_contact_pairs_result* result) {
    return result && result->mPairs.mLocalCoordinates.Size() > 0
               ? result->mPairs.mLocalCoordinates.As<double>()
               : nullptr;
}
const double* mio_contact_pairs_closest_point(const mio_contact_pairs_result* result) {
    return result && result->mPairs.mClosestPoint.Size() > 0
               ? result->mPairs.mClosestPoint.As<double>()
               : nullptr;
}
const double* mio_contact_pairs_gap(const mio_contact_pairs_result* result) {
    return result && result->mPairs.mGap.Size() > 0 ? result->mPairs.mGap.As<double>() : nullptr;
}
const double* mio_contact_pairs_normal(const mio_contact_pairs_result* result) {
    return result && result->mPairs.mNormal.Size() > 0 ? result->mPairs.mNormal.As<double>()
                                                       : nullptr;
}
const int64_t* mio_contact_pairs_unmatched(const mio_contact_pairs_result* result) {
    return result && result->mPairs.mUnmatched.Size() > 0
               ? result->mPairs.mUnmatched.As<std::int64_t>()
               : nullptr;
}

void mio_contact_pairs_result_free(mio_contact_pairs_result* result) {
    delete result;
}

mio_mesh* mio_split_interface(const mio_mesh* mesh, const mio_region_selector* side,
                              const mio_split_interface_opts* opts, int64_t* num_duplicated_points,
                              int64_t* num_cohesive_cells) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: split_interface: mesh is NULL");
        meshioplusplus::SplitInterfaceOptions options;
        if (opts)
            options.mAddCohesive = opts->add_cohesive != 0;
        meshioplusplus::SplitInterfaceResult result =
            meshioplusplus::split_interface(mesh->mMesh, capi_selector(side), options);
        if (num_duplicated_points)
            *num_duplicated_points = result.mNumDuplicatedPoints;
        if (num_cohesive_cells)
            *num_cohesive_cells = result.mNumCohesiveCells;
        return new mio_mesh{std::move(result.mMesh)};
    });
}

mio_status mio_mesh_remove_region(mio_mesh* mesh, int64_t index) {
    return guarded([&]() -> mio_status {
        if (!mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: mesh is NULL");
        if (index < 0 || static_cast<std::size_t>(index) >= mesh->mMesh.NumRegions())
            return fail(MIO_ERR_INVALID_ARG, "meshio++: region index out of range");
        mesh->mMesh.RemoveRegion(static_cast<std::size_t>(index));
        return MIO_OK;
    });
}

void mio_periodic_opts_init(mio_periodic_opts* opts) {
    if (!opts)
        return;
    *opts = mio_periodic_opts{};
    for (int k = 0; k < 4; ++k)
        opts->matrix[k * 5] = 1.0;
    opts->atol = 1e-8;
    opts->require_complete = 1;
}

mio_periodic_pairs* mio_match_periodic_nodes(const mio_mesh* mesh, const mio_region_selector* slave,
                                             const mio_region_selector* master,
                                             const mio_periodic_opts* opts) {
    return guarded_ptr(static_cast<mio_periodic_pairs*>(nullptr), [&]() -> mio_periodic_pairs* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        meshioplusplus::PeriodicOptions options;
        if (opts) {
            options.mTransform = meshioplusplus::transform_from_matrix(opts->matrix);
            options.mAtol = opts->atol;
            options.mRequireComplete = opts->require_complete != 0;
        }
        return new mio_periodic_pairs{meshioplusplus::match_periodic_nodes(
            mesh->mMesh, capi_selector(slave), capi_selector(master), options)};
    });
}

mio_status mio_periodic_pairs_info(const mio_periodic_pairs* pairs, int64_t* num_pairs,
                                   int64_t* num_unmatched, int64_t* num_fixed,
                                   double* max_residual) {
    return guarded([&]() -> mio_status {
        if (!pairs)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: pairs is NULL");
        if (num_pairs)
            *num_pairs = static_cast<int64_t>(pairs->mPairs.mSlave.Size());
        if (num_unmatched)
            *num_unmatched = static_cast<int64_t>(pairs->mPairs.mUnmatched.Size());
        if (num_fixed)
            *num_fixed = pairs->mPairs.mNumFixed;
        if (max_residual)
            *max_residual = pairs->mPairs.mMaxResidual;
        return MIO_OK;
    });
}

const int64_t* mio_periodic_pairs_slave(const mio_periodic_pairs* pairs, int64_t* count) {
    if (count)
        *count = 0;
    return pairs ? capi_periodic_ids(pairs->mPairs.mSlave, count) : nullptr;
}

const int64_t* mio_periodic_pairs_master(const mio_periodic_pairs* pairs, int64_t* count) {
    if (count)
        *count = 0;
    return pairs ? capi_periodic_ids(pairs->mPairs.mMaster, count) : nullptr;
}

const int64_t* mio_periodic_pairs_unmatched(const mio_periodic_pairs* pairs, int64_t* count) {
    if (count)
        *count = 0;
    return pairs ? capi_periodic_ids(pairs->mPairs.mUnmatched, count) : nullptr;
}

void mio_periodic_pairs_free(mio_periodic_pairs* pairs) {
    delete pairs;
}

// --- quality gate (v16.24.0) -------------------------------------------------

static_assert(sizeof(mio_quality_gate_report) == 88,
              "mio_quality_gate_report grew outside its reserved tail");

mio_status mio_check_quality(const mio_mesh* mesh, const char* spec, int64_t max_inverted,
                             int64_t max_degenerate, mio_quality_gate_report* report, char* summary,
                             int64_t summary_len) {
    return guarded([&]() -> mio_status {
        if (!mesh || !report)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: check_quality: mesh/report is NULL");
        meshioplusplus::QualityGateOptions options;
        if (spec)
            options.mThresholds = meshioplusplus::parse_quality_thresholds(spec);
        options.mMaxInverted = max_inverted;
        options.mMaxDegenerate = max_degenerate;
        const meshioplusplus::QualityGateResult r =
            meshioplusplus::check_quality(mesh->mMesh, options);
        const std::string text = meshioplusplus::quality_gate_summary(r);
        *report = mio_quality_gate_report{};
        report->passed = r.mPassed ? 1 : 0;
        report->num_checks = static_cast<int64_t>(r.mChecks.size());
        for (const meshioplusplus::QualityCheck& c : r.mChecks)
            report->num_failed += c.mPassed ? 0 : 1;
        report->num_cells = r.mReport.mNumCells;
        report->num_inverted = r.mReport.mNumInverted;
        report->num_degenerate = r.mReport.mNumDegenerate;
        report->summary_length = static_cast<int64_t>(text.size());
        if (summary && summary_len > 0) {
            const std::size_t n = std::min(text.size(), static_cast<std::size_t>(summary_len - 1));
            std::memcpy(summary, text.data(), n);
            summary[n] = '\0';
        }
        return MIO_OK;
    });
}

// --- software rendering (v16.33.0, field rendering v16.34.0) ------------------

static_assert(sizeof(mio_render_opts) == 432, "mio_render_opts grew outside its reserved tail");
static_assert(sizeof(mio_text_opts) == 96, "mio_text_opts grew outside its reserved tail");
static_assert(sizeof(mio_snapshot_opts) == 72, "mio_snapshot_opts grew outside its reserved tail");

namespace {

meshioplusplus::RenderColor render_color(std::uint32_t Rgba) {
    return {static_cast<std::uint8_t>(Rgba >> 24), static_cast<std::uint8_t>(Rgba >> 16),
            static_cast<std::uint8_t>(Rgba >> 8), static_cast<std::uint8_t>(Rgba)};
}

std::string render_string(const char* pText) {
    return pText ? std::string(pText) : std::string();
}

meshioplusplus::RenderOptions render_options_from_c(const mio_render_opts* pOpts) {
    mio_render_opts defaults;
    if (!pOpts) {
        mio_render_opts_init(&defaults);
        pOpts = &defaults;
    }
    for (std::int64_t r : pOpts->reserved)
        if (r != 0)
            throw std::invalid_argument("meshio++: render: reserved option fields must be zero");
    meshioplusplus::RenderOptions o;
    o.mWidth = pOpts->width;
    o.mHeight = pOpts->height;
    o.mPixelAspect = pOpts->pixel_aspect;
    o.mSupersample = pOpts->supersample;
    o.mAzimuth = pOpts->azimuth;
    o.mElevation = pOpts->elevation;
    o.mRoll = pOpts->roll;
    o.mView = render_string(pOpts->view);
    o.mProjection = pOpts->perspective != 0 ? meshioplusplus::RenderProjection::Perspective
                                            : meshioplusplus::RenderProjection::Orthographic;
    o.mFovDeg = pOpts->fov_deg;
    o.mZoom = pOpts->zoom;
    o.mPanX = pOpts->pan_x;
    o.mPanY = pOpts->pan_y;
    switch (pOpts->shading) {
        case MIO_SHADING_NONE:
            o.mShading = meshioplusplus::RenderShading::None;
            break;
        case MIO_SHADING_FLAT:
            o.mShading = meshioplusplus::RenderShading::Flat;
            break;
        case MIO_SHADING_SMOOTH:
            o.mShading = meshioplusplus::RenderShading::Smooth;
            break;
        default:
            throw std::invalid_argument("meshio++: render: unknown shading value");
    }
    o.mTwoSided = pOpts->two_sided != 0;
    o.mAmbient = pOpts->ambient;
    o.mLightDir = {pOpts->light_x, pOpts->light_y, pOpts->light_z};
    o.mSplitAngle = pOpts->split_angle;
    switch (pOpts->edges) {
        case MIO_EDGES_NONE:
            o.mEdges = meshioplusplus::RenderEdges::None;
            break;
        case MIO_EDGES_ALL:
            o.mEdges = meshioplusplus::RenderEdges::All;
            break;
        case MIO_EDGES_FEATURE:
            o.mEdges = meshioplusplus::RenderEdges::Feature;
            break;
        default:
            throw std::invalid_argument("meshio++: render: unknown edges value");
    }
    o.mFeatureAngle = pOpts->feature_angle;
    o.mEdgeColor = render_color(pOpts->edge_color);
    o.mFillColor = render_color(pOpts->fill_color);
    o.mLineColor = render_color(pOpts->line_color);
    o.mBackground = render_color(pOpts->background);
    o.mNanColor = render_color(pOpts->nan_color);
    o.mPointRadius = pOpts->point_radius;
    o.mColorBy = render_string(pOpts->color_by);
    if (pOpts->component >= 0)
        o.mComponent = pOpts->component;
    if (pOpts->cmap && pOpts->cmap[0] != '\0')
        o.mCmap = pOpts->cmap;
    if (pOpts->has_vmin)
        o.mVMin = pOpts->vmin;
    if (pOpts->has_vmax)
        o.mVMax = pOpts->vmax;
    o.mColorbar = pOpts->colorbar != 0;
    o.mAxes = pOpts->axes != 0;
    o.mScaleBar = pOpts->scale_bar != 0;
    o.mReduce = render_string(pOpts->reduce);
    o.mExpr = render_string(pOpts->expr);
    if (pOpts->has_clip_low)
        o.mClipLow = pOpts->clip_low;
    if (pOpts->has_clip_high)
        o.mClipHigh = pOpts->clip_high;
    o.mSymmetric = pOpts->symmetric != 0;
    switch (pOpts->scale) {
        case MIO_SCALE_LINEAR:
            o.mScale = meshioplusplus::RenderScale::Linear;
            break;
        case MIO_SCALE_LOG:
            o.mScale = meshioplusplus::RenderScale::Log;
            break;
        case MIO_SCALE_SYMLOG:
            o.mScale = meshioplusplus::RenderScale::Symlog;
            break;
        default:
            throw std::invalid_argument("meshio++: render: unknown scale value");
    }
    o.mScaleThreshold = pOpts->scale_threshold;
    o.mCategorical = pOpts->categorical != 0;
    o.mColorRegions = pOpts->color_regions != 0;
    o.mCategoryEdges = pOpts->category_edges != 0;
    o.mIsolines = pOpts->isolines;
    if (pOpts->num_iso_levels < 0 || (pOpts->num_iso_levels > 0 && !pOpts->iso_levels))
        throw std::invalid_argument("meshio++: render: iso_levels is NULL or its count negative");
    if (pOpts->num_iso_levels > 0)
        o.mIsoLevels.assign(pOpts->iso_levels, pOpts->iso_levels + pOpts->num_iso_levels);
    o.mIsoColor = render_color(pOpts->iso_color);
    o.mVectors = render_string(pOpts->vectors);
    o.mVectorCount = pOpts->vector_count;
    o.mVectorLength = pOpts->vector_length;
    o.mVectorColor = render_color(pOpts->vector_color);
    o.mWarp = render_string(pOpts->warp);
    o.mWarpScale = pOpts->warp_scale;
    o.mWarpOutline = pOpts->warp_outline != 0;
    o.mOutlineColor = render_color(pOpts->outline_color);
    if (pOpts->diagnostic < MIO_DIAGNOSTIC_NONE || pOpts->diagnostic > MIO_DIAGNOSTIC_EDGE_LENGTH)
        throw std::invalid_argument("meshio++: render: unknown diagnostic value");
    o.mDiagnostic = static_cast<meshioplusplus::RenderDiagnostic>(pOpts->diagnostic);
    o.mQualityMetric = render_string(pOpts->quality_metric);
    if (pOpts->num_cutaways < 0 || pOpts->num_cutaways > 2 ||
        (pOpts->num_cutaways > 0 && !pOpts->cutaways))
        throw std::invalid_argument(
            "meshio++: render: cutaways is NULL, or its count is not 0, 1 or 2");
    for (std::int32_t i = 0; i < pOpts->num_cutaways; ++i) {
        const double* six = pOpts->cutaways + 6 * i;
        meshioplusplus::RenderCutaway plane;
        plane.mPoint = {six[0], six[1], six[2]};
        plane.mNormal = {six[3], six[4], six[5]};
        o.mCutaways.push_back(plane);
    }
    o.mCutawayTint = render_color(pOpts->cutaway_tint);
    o.mStreamlines = render_string(pOpts->streamlines);
    o.mStreamSeeds = pOpts->stream_seeds;
    o.mStreamLength = pOpts->stream_length;
    o.mStreamColor = render_color(pOpts->stream_color);
    return o;
}

meshioplusplus::TextOptions text_options_from_c(const mio_text_opts* pOpts) {
    mio_text_opts defaults;
    if (!pOpts) {
        mio_text_opts_init(&defaults);
        pOpts = &defaults;
    }
    for (std::int64_t r : pOpts->reserved)
        if (r != 0)
            throw std::invalid_argument("meshio++: render: reserved option fields must be zero");
    if (pOpts->flags != 0)
        throw std::invalid_argument("meshio++: render: text option flags must be zero");
    if (pOpts->encoding < MIO_ENCODING_HALFBLOCK || pOpts->encoding > MIO_ENCODING_SIXEL)
        throw std::invalid_argument("meshio++: render: unknown encoding value");
    if (pOpts->color_depth < MIO_COLOR_TRUECOLOR || pOpts->color_depth > MIO_COLOR_MONO)
        throw std::invalid_argument("meshio++: render: unknown color_depth value");
    if (pOpts->format < MIO_TEXT_ANSI || pOpts->format > MIO_TEXT_HTML)
        throw std::invalid_argument("meshio++: render: unknown format value");
    meshioplusplus::TextOptions t;
    t.mEncoding = static_cast<meshioplusplus::TextEncoding>(pOpts->encoding);
    t.mDepth = static_cast<meshioplusplus::ColorDepth>(pOpts->color_depth);
    t.mFormat = static_cast<meshioplusplus::TextFormat>(pOpts->format);
    t.mCols = pOpts->cols;
    t.mRows = pOpts->rows;
    t.mCellAspect = pOpts->cell_aspect;
    t.mCellPixelWidth = pOpts->cell_pixel_width;
    t.mCellPixelHeight = pOpts->cell_pixel_height;
    t.mTmuxPassthrough = pOpts->tmux != 0;
    t.mNotes = pOpts->notes != 0;
    return t;
}

const meshioplusplus::Frame& frame_of(const mio_frame* pFrame) {
    if (!pFrame)
        throw std::invalid_argument("meshio++: frame is NULL");
    return pFrame->mFrame;
}

}  // namespace

void mio_render_opts_init(mio_render_opts* opts) {
    if (!opts)
        return;
    *opts = mio_render_opts{};
    opts->pixel_aspect = 1.0;
    opts->azimuth = 45.0;
    opts->elevation = 35.264389682754654;
    opts->fov_deg = 30.0;
    opts->zoom = 1.0;
    opts->ambient = 0.25;
    opts->light_z = 1.0;
    opts->split_angle = 30.0;
    opts->feature_angle = 30.0;
    opts->point_radius = 1.5;
    opts->scale_threshold = 1.0;
    opts->warp_scale = 1.0;
    opts->width = 320;
    opts->height = 240;
    opts->supersample = 1;
    opts->shading = MIO_SHADING_FLAT;
    opts->two_sided = 1;
    opts->component = -1;
    opts->vector_count = 200;
    opts->edge_color = 0x000000FFu;
    opts->fill_color = 0xC8C5BDFFu;
    opts->line_color = 0x000080FFu;
    opts->background = 0x00000000u;
    opts->nan_color = 0x808080FFu;
    opts->iso_color = 0x1E1E1EFFu;
    opts->vector_color = 0xDC3232FFu;
    opts->outline_color = 0x969696FFu;
    opts->cutaway_tint = 0xE88034FFu;
    opts->stream_seeds = 40;
    opts->stream_length = 0.5;
    opts->stream_color = 0xF050A0FFu;
}

void mio_text_opts_init(mio_text_opts* opts) {
    if (!opts)
        return;
    *opts = mio_text_opts{};
    opts->cell_aspect = 2.0;
    opts->cols = 80;
    opts->rows = 24;
    opts->cell_pixel_width = 8;
    opts->cell_pixel_height = 16;
    opts->notes = 1;
}

void mio_snapshot_opts_init(mio_snapshot_opts* opts) {
    if (!opts)
        return;
    *opts = mio_snapshot_opts{};
    opts->cast_fps = 12.0;
    opts->cast_degrees = 360.0;
    opts->cast_frames = 36;
}

mio_frame* mio_render(const mio_mesh* mesh, const mio_render_opts* opts) {
    return guarded_ptr(static_cast<mio_frame*>(nullptr), [&]() -> mio_frame* {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        const meshioplusplus::RenderOptions options = render_options_from_c(opts);
        return new mio_frame{meshioplusplus::render(mesh->mMesh, options)};
    });
}

void mio_frame_free(mio_frame* frame) {
    delete frame;
}

int32_t mio_frame_width(const mio_frame* frame) {
    return frame ? frame->mFrame.mWidth : 0;
}

int32_t mio_frame_height(const mio_frame* frame) {
    return frame ? frame->mFrame.mHeight : 0;
}

const uint8_t* mio_frame_rgba(const mio_frame* frame) {
    return frame ? frame->mFrame.mRgba.data() : nullptr;
}

const int64_t* mio_frame_cell_ids(const mio_frame* frame) {
    return frame ? frame->mFrame.mCellIds.data() : nullptr;
}

int32_t mio_frame_range(const mio_frame* frame, double* vmin, double* vmax) {
    return guarded_ptr(int32_t{-1}, [&]() -> int32_t {
        const meshioplusplus::Frame& f = frame_of(frame);
        if (vmin)
            *vmin = f.mVMin;
        if (vmax)
            *vmax = f.mVMax;
        return f.mColored ? 1 : 0;
    });
}

int64_t mio_frame_num_notes(const mio_frame* frame) {
    return guarded_ptr(int64_t{-1}, [&]() -> int64_t {
        return static_cast<int64_t>(frame_of(frame).mNotes.size());
    });
}

int64_t mio_frame_note(const mio_frame* frame, int64_t index, char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() -> int64_t {
        const meshioplusplus::Frame& f = frame_of(frame);
        if (index < 0 || static_cast<std::size_t>(index) >= f.mNotes.size())
            throw std::out_of_range("meshio++: frame note index out of range");
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative note buffer length");
        return copy_string(f.mNotes[static_cast<std::size_t>(index)], buf, buflen);
    });
}

int64_t mio_frame_text(const mio_frame* frame, const mio_text_opts* opts, char* buf,
                       int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() -> int64_t {
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative text buffer length");
        const std::string text =
            meshioplusplus::encode_text(frame_of(frame), text_options_from_c(opts));
        return copy_string(text, buf, buflen);
    });
}

int64_t mio_frame_png(const mio_frame* frame, int32_t compress, uint8_t* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() -> int64_t {
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative PNG buffer length");
        const std::string png = meshioplusplus::encode_png(frame_of(frame), compress);
        if (buf && buflen > 0)
            std::memcpy(buf, png.data(), std::min(png.size(), static_cast<std::size_t>(buflen)));
        return static_cast<int64_t>(png.size());
    });
}

int64_t mio_render_text(const mio_mesh* mesh, const mio_render_opts* render,
                        const mio_text_opts* text, char* buf, int64_t buflen) {
    return guarded_ptr(int64_t{-1}, [&]() -> int64_t {
        if (!mesh)
            throw meshioplusplus::ReadError("meshio++: mesh is NULL");
        if (buflen < 0)
            throw std::invalid_argument("meshio++: negative text buffer length");
        const std::string out = meshioplusplus::render_text(
            mesh->mMesh, render_options_from_c(render), text_options_from_c(text));
        return copy_string(out, buf, buflen);
    });
}

mio_status mio_write_snapshot(const char* path, const mio_mesh* mesh, const mio_render_opts* render,
                              const mio_text_opts* text, const mio_snapshot_opts* snapshot) {
    return guarded([&]() -> mio_status {
        if (!path || !mesh)
            return fail(MIO_ERR_INVALID_ARG, "meshio++: write_snapshot: path/mesh is NULL");
        meshioplusplus::SnapshotOptions options;
        if (snapshot) {
            for (std::int64_t r : snapshot->reserved)
                if (r != 0)
                    return fail(MIO_ERR_INVALID_ARG,
                                "meshio++: write_snapshot: reserved fields must be zero");
            options.mPngCompress = snapshot->png_compress;
            options.mCastFrames = snapshot->cast_frames;
            options.mCastFps = snapshot->cast_fps;
            options.mCastDegrees = snapshot->cast_degrees;
        }
        meshioplusplus::write_snapshot(path, mesh->mMesh, render_options_from_c(render),
                                       text_options_from_c(text), options);
        return MIO_OK;
    });
}

int32_t mio_detect_color_depth(const char* no_color, const char* color_term, const char* term) {
    switch (meshioplusplus::detect_color_depth(no_color, color_term, term)) {
        case meshioplusplus::ColorDepth::TrueColor:
            return MIO_COLOR_TRUECOLOR;
        case meshioplusplus::ColorDepth::Palette256:
            return MIO_COLOR_256;
        case meshioplusplus::ColorDepth::Ansi16:
            return MIO_COLOR_16;
        default:
            return MIO_COLOR_MONO;
    }
}

int mio_has_mmg(void) {
    return meshioplusplus::mmg_available() ? 1 : 0;
}

int mio_has_mmg_elas(void) {
    return meshioplusplus::mmg_has_elas() ? 1 : 0;
}

int mio_mmg_index_bits(void) {
    return meshioplusplus::mmg_index_bits();
}

void mio_mmg_params_init(mio_mmg_params* params) {
    if (!params)
        return;
    *params = mio_mmg_params{};
    params->hmin = params->hmax = params->hausd = params->hgrad = params->hgradreq = -1.0;
    params->angle_detection = 45.0;
    params->verbose = -1;
    params->mem = -1;
}

void mio_adapt_opts_init(mio_adapt_opts* opts) {
    if (!opts)
        return;
    *opts = mio_adapt_opts{};
    mio_mmg_params_init(&opts->params);
    opts->hsiz = -1.0;
}

mio_mesh* mio_adapt(const mio_mesh* mesh, const mio_adapt_opts* opts, mio_mmg_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh)
            throw std::invalid_argument("meshio++: adapt: mesh is NULL");
        mio_adapt_opts defaults;
        mio_adapt_opts_init(&defaults);
        const mio_adapt_opts& r_in = opts ? *opts : defaults;
        meshioplusplus::AdaptOptions o;
        o.mParameters = capi_mmg_params(mesh, r_in.params, "adapt");
        o.mHsiz = r_in.hsiz;
        if (r_in.size)
            o.mSize = r_in.size;
        if (r_in.metric)
            o.mMetric = r_in.metric;
        meshioplusplus::AdaptResult r = meshioplusplus::adapt(mesh->mMesh, o);
        capi_mmg_report(r.mReport, report);
        return new mio_mesh{std::move(r.mMesh)};
    });
}

void mio_discretize_levelset_opts_init(mio_discretize_levelset_opts* opts) {
    if (!opts)
        return;
    *opts = mio_discretize_levelset_opts{};
    mio_mmg_params_init(&opts->params);
    opts->remove_small_components = -1.0;
}

mio_mesh* mio_discretize_levelset(const mio_mesh* mesh, const mio_discretize_levelset_opts* opts,
                                  mio_mmg_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !opts || !opts->level_set)
            throw std::invalid_argument(
                "meshio++: discretize_levelset: mesh, opts or level_set is NULL");
        meshioplusplus::DiscretizeLevelsetOptions o;
        o.mParameters = capi_mmg_params(mesh, opts->params, "discretize_levelset");
        o.mLevelSet = opts->level_set;
        o.mIsovalue = opts->isovalue;
        o.mRemoveSmallComponents = opts->remove_small_components;
        for (int64_t i = 0; opts->materials && i < opts->num_materials; ++i)
            o.mMaterials.push_back({opts->materials[i].ref, opts->materials[i].split != 0,
                                    opts->materials[i].ref_inside, opts->materials[i].ref_outside});
        meshioplusplus::DiscretizeLevelsetResult r =
            meshioplusplus::discretize_levelset(mesh->mMesh, o);
        capi_mmg_report(r.mReport, report);
        return new mio_mesh{std::move(r.mMesh)};
    });
}

void mio_mesh_domain_opts_init(mio_mesh_domain_opts* opts) {
    if (!opts)
        return;
    *opts = mio_mesh_domain_opts{};
    mio_mmg_params_init(&opts->params);
    opts->hsiz = -1.0;
}

mio_mesh* mio_mesh_domain(const mio_mesh* boundary, const mio_mesh_domain_opts* opts,
                          mio_mmg_report* report, int64_t* num_subdomains, int64_t* num_holes) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!boundary)
            throw std::invalid_argument("meshio++: mesh_domain: boundary is NULL");
        mio_mesh_domain_opts defaults;
        mio_mesh_domain_opts_init(&defaults);
        const mio_mesh_domain_opts& r_in = opts ? *opts : defaults;
        meshioplusplus::MeshDomainOptions o;
        o.mParameters = capi_mmg_params(boundary, r_in.params, "mesh_domain");
        o.mHsiz = r_in.hsiz;
        if (r_in.size)
            o.mSize = r_in.size;
        if (r_in.holes && r_in.num_holes > 0)
            o.mHoles.assign(r_in.holes, r_in.holes + 2 * r_in.num_holes);
        meshioplusplus::MeshDomainResult r = meshioplusplus::mesh_domain(boundary->mMesh, o);
        capi_mmg_report(r.mReport, report);
        if (num_subdomains)
            *num_subdomains = r.mNumSubdomains;
        if (num_holes)
            *num_holes = r.mNumHoles;
        return new mio_mesh{std::move(r.mMesh)};
    });
}

void mio_move_opts_init(mio_move_opts* opts) {
    if (!opts)
        return;
    *opts = mio_move_opts{};
    mio_mmg_params_init(&opts->params);
    opts->lagrangian = 1;
}

mio_mesh* mio_move(const mio_mesh* mesh, const mio_move_opts* opts, mio_mmg_report* report) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !opts || !opts->displacement)
            throw std::invalid_argument("meshio++: move: mesh, opts or displacement is NULL");
        meshioplusplus::MoveOptions o;
        o.mParameters = capi_mmg_params(mesh, opts->params, "move");
        o.mDisplacement = opts->displacement;
        for (int64_t i = 0; opts->moving && i < opts->num_moving; ++i)
            if (opts->moving[i])
                o.mMoving.emplace_back(opts->moving[i]);
        o.mLagrangian = opts->lagrangian;
        meshioplusplus::AdaptResult r = meshioplusplus::move(mesh->mMesh, o);
        capi_mmg_report(r.mReport, report);
        return new mio_mesh{std::move(r.mMesh)};
    });
}

mio_mesh* mio_metric_from_hessian(const mio_mesh* mesh, const char* hessian, double error,
                                  double hmin, double hmax, const char* output) {
    return guarded_ptr(static_cast<mio_mesh*>(nullptr), [&]() -> mio_mesh* {
        if (!mesh || !hessian)
            throw std::invalid_argument("meshio++: metric_from_hessian: mesh or hessian is NULL");
        meshioplusplus::MetricFromHessianOptions o;
        o.mHessian = hessian;
        o.mError = error;
        o.mHmin = hmin;
        o.mHmax = hmax;
        if (output)
            o.mOutput = output;
        return new mio_mesh{meshioplusplus::metric_from_hessian(mesh->mMesh, o)};
    });
}

}  // extern "C"
