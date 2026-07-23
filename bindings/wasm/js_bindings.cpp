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
 * @file js_bindings.cpp
 * @brief Emscripten/embind entry point for the `@meshioplusplus/wasm` npm
 *        package. Compiled only under `MESHIOPLUSPLUS_BUILD_WASM` (see the
 *        `if(EMSCRIPTEN)` block in the top-level `CMakeLists.txt`).
 *
 * Unlike `bindings/_core.cpp` (which exposes one `<fmt>_read`/`<fmt>_write`
 * pair per format to Python, letting `_helpers.py` own extension-dispatch and
 * `np_conversions.hpp` own zero-copy numpy<->NDArray conversion), this file
 * exposes a small, flat, copy-based JS API: `readMesh`/`writeMesh`/`convert`
 * plus the two pure cell-type metadata tables. There is no zero-copy path
 * here by design -- WASM linear memory and the JS heap are different address
 * spaces, so every value crossing the boundary is copied once (mirroring how
 * the existing pybind11 layer already treats ragged cell blocks: "just copy,
 * it's fine"). `NDArray`/`CellBlock`/`Mesh` are therefore kept entirely
 * internal; JS only ever sees plain objects of typed arrays (see
 * `mesh_to_val`/`val_to_mesh` below for the exact shape). Note: the
 * JS-facing names bound below (`"readMesh"`, `"writeMesh"`, ...) are string
 * literals independent of the C++ function names/symbols on the other side
 * of each `emscripten::function(...)` call, so this file's internal C++
 * identifiers can follow the project's snake_case free-function convention
 * without changing the JS API surface.
 *
 * Format scope: whatever `registry.cpp` was compiled with -- this file has no
 * format table of its own, so it needs no edit when one is added or when an
 * optional dependency comes or goes. Since v8.0.0 that is every format the
 * core has, CGNS/H5M/HMF/MED/Exodus and XDMF's HDF data path included:
 * build/build-wasm-deps.sh produces the wasm32 HDF5/netCDF the build links
 * (see doc/wasm.md). Ambiguous extensions (`.msh` shared
 * by ansys/freefem/gmsh, `.inp` shared by abaqus/ansysinp) require an
 * explicit `format` argument, mirroring Python's `file_format=` kwarg;
 * `.msh` defaults to gmsh and `.inp` to abaqus when `format` is omitted.
 *
 * File I/O goes through Emscripten's virtual filesystem (`Module.FS`,
 * exposed via `-sEXPORTED_RUNTIME_METHODS=['FS']` in CMakeLists.txt) -- JS
 * callers write bytes into the virtual FS themselves before calling
 * `readMesh`, and read them back out after `writeMesh`/`convert`.
 *
 * Exceptions: `meshioplusplus::ReadError`/`WriteError` derive from
 * `std::exception`, so Emscripten's default (exceptions enabled) embind
 * configuration automatically surfaces them as catchable JS errors -- no
 * explicit translator is needed here (unlike `bindings/_core.cpp`'s
 * `py::register_exception_translator`, which exists because CPython has no
 * such automatic C++-exception-to-host-exception bridge).
 */

// System includes
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <type_traits>
#include <vector>

// External includes
#include <emscripten/bind.h>
#include <emscripten/emscripten.h>
#include <emscripten/val.h>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/convert_cells.hpp"
#include "meshioplusplus/operations/crop.hpp"
#include "meshioplusplus/operations/data_average.hpp"
#include "meshioplusplus/operations/data_calc.hpp"
#include "meshioplusplus/operations/data_common.hpp"
#include "meshioplusplus/operations/data_condition.hpp"
#include "meshioplusplus/operations/data_info.hpp"
#include "meshioplusplus/operations/data_manage.hpp"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/diff.hpp"
#include "meshioplusplus/operations/interpolate.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/operations/isosurface.hpp"
#include "meshioplusplus/operations/slice.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/operations/split.hpp"
#include "meshioplusplus/operations/stats.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/operations/transform.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/skin.hpp"
#include "meshioplusplus/types.hpp"

using emscripten::val;

namespace {

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

// ---------------------------------------------------------------------
// Typed-array helpers: copy a C++ buffer into a genuinely JS-owned typed
// array (a transient emscripten::typed_memory_view over Module memory is
// only valid for the duration of this call; `.set()` copies it into `arr`'s
// own backing store immediately, so the returned val outlives the view).
// ---------------------------------------------------------------------

val float64_array_from(const double* pData, std::size_t n) {
    val arr = val::global("Float64Array").new_(n);
    arr.call<void>("set", val(emscripten::typed_memory_view(n, pData)));
    return arr;
}

val int32_array_from(const std::int32_t* pData, std::size_t n) {
    val arr = val::global("Int32Array").new_(n);
    arr.call<void>("set", val(emscripten::typed_memory_view(n, pData)));
    return arr;
}

std::size_t cols_of(const NDArray& rA) {
    return rA.Shape().size() >= 2 ? rA.Shape()[1] : 1;
}

// Any-dtype NDArray -> Float64Array (upcasts ints / lower-precision floats;
// every point/data array in the JS API is double-precision for simplicity).
val ndarray_to_float64_array(const NDArray& rA) {
    return meshioplusplus::detail::dispatch_dtype(rA.Dtype(), [&]<class T>() -> val {
        if constexpr (std::is_same_v<T, double>) {
            return float64_array_from(rA.As<double>(), rA.Size());
        } else {
            std::vector<double> tmp(rA.Size());
            const T* src = rA.As<T>();
            for (std::size_t i = 0; i < rA.Size(); ++i)
                tmp[i] = static_cast<double>(src[i]);
            return float64_array_from(tmp.data(), tmp.size());
        }
    });
}

// Integer-dtype NDArray (mesh connectivity, always Int64 in the C++ core) ->
// Int32Array. Node/point counts for any mesh a browser can reasonably handle
// fit comfortably in 32 bits; Int32Array is far more JS-ergonomic than
// BigInt64Array for typical mesh-processing consumer code.
val ndarray_to_int32_array(const NDArray& rA) {
    std::vector<std::int32_t> tmp(rA.Size());
    meshioplusplus::detail::dispatch_dtype(rA.Dtype(), [&]<class T>() {
        const T* src = rA.As<T>();
        for (std::size_t i = 0; i < rA.Size(); ++i)
            tmp[i] = static_cast<std::int32_t>(src[i]);
    });
    return int32_array_from(tmp.data(), tmp.size());
}

/**
 * @brief Convert a C++ `Mesh` into a plain JS object of typed arrays.
 *
 * Shape: `{ points: Float64Array, dim: number, cells: [{type, data:
 * Int32Array, nodesPerCell}], point_data: {name: Float64Array}, cell_data:
 * {name: Float64Array[]} (one array per cell block, same order as `cells`),
 * field_data: {name: Float64Array} }` -- deliberately mirrors the Python
 * `Mesh`'s structure (points, a list of cell blocks, cell_data as one array
 * per block) for consistency with the rest of meshio++.
 *
 * @throws meshioplusplus::ReadError if any cell block is ragged (polygon/
 *   polyhedron with varying node counts) -- not supported by the v1 JS API.
 */
val mesh_to_val(const Mesh& rMesh) {
    val out = val::object();
    const NDArray& points = rMesh.Points();
    out.set("points", ndarray_to_float64_array(points));
    out.set("dim", static_cast<int>(cols_of(points)));

    val cells = val::array();
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsRagged())
            throw meshioplusplus::ReadError("meshio++ (wasm): ragged cell blocks ('" + cb.Type() +
                                            "') are not supported by the JS API yet");
        const NDArray& conn = cb.Conn();
        val block = val::object();
        block.set("type", cb.Type());
        block.set("data", ndarray_to_int32_array(conn));
        block.set("nodesPerCell", static_cast<int>(cols_of(conn)));
        cells.call<void>("push", block);
    }
    out.set("cells", cells);

    val point_data = val::object();
    for (const auto& name : rMesh.PointDataNames())
        point_data.set(name, ndarray_to_float64_array(rMesh.PointData(name)));
    out.set("point_data", point_data);

    val cell_data = val::object();
    for (const auto& name : rMesh.CellDataNames()) {
        val blocks = val::array();
        for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name); ++b)
            blocks.call<void>("push", ndarray_to_float64_array(rMesh.CellData(name, b)));
        cell_data.set(name, blocks);
    }
    out.set("cell_data", cell_data);

    val field_data = val::object();
    for (const auto& name : rMesh.FieldDataNames())
        field_data.set(name, ndarray_to_float64_array(rMesh.FieldData(name)));
    out.set("field_data", field_data);

    // Named regions (doc/regions.md) ride on the mesh object itself rather than
    // through a function of their own, so `readMesh` / `writeMesh` / `convert`
    // carry them with no extra call. `entries` is flat: `n` indices for
    // point/cell, `n` (cell, facet) pairs for side.
    val regions = val::array();
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        val jr = val::object();
        jr.set("name", r.mName);
        jr.set("kind", std::string(meshioplusplus::region_kind_name(r.mKind)));
        jr.set("dim", r.mDim);
        jr.set("tag", static_cast<double>(r.mTag));
        jr.set("entries", ndarray_to_int32_array(r.mEntries));
        regions.call<void>("push", jr);
    }
    out.set("regions", regions);

    return out;
}

// A JS array-like of numbers -> an owning Float64/Int64 NDArray of `shape`.
// `emscripten::vecFromJSArray` copies once into a std::vector; the second
// copy into the NDArray's own buffer is unavoidable without exposing
// NDArray's internals to JS, which the file-level design deliberately avoids.
NDArray float64_ndarray_from_val(const val& rJsArr, std::vector<std::size_t> shape) {
    std::vector<double> tmp = emscripten::vecFromJSArray<double>(rJsArr);
    NDArray out = NDArray::Uninit(DType::Float64, std::move(shape));
    std::copy(tmp.begin(), tmp.end(), out.As<double>());
    return out;
}

NDArray int64_ndarray_from_val(const val& rJsArr, std::vector<std::size_t> shape) {
    std::vector<std::int64_t> tmp = emscripten::vecFromJSArray<std::int64_t>(rJsArr);
    NDArray out = NDArray::Uninit(DType::Int64, std::move(shape));
    std::copy(tmp.begin(), tmp.end(), out.As<std::int64_t>());
    return out;
}

std::vector<std::string> js_object_keys(const val& rObj) {
    val keys = val::global("Object").call<val>("keys", rObj);
    return emscripten::vecFromJSArray<std::string>(keys);
}

/**
 * @brief Convert a plain JS mesh object (see `mesh_to_val`'s shape) into a C++
 * `Mesh`, for `writeMesh`/`convert`.
 *
 * @throws meshioplusplus::WriteError on malformed input (points/cell-block
 *   lengths not divisible by their declared dim/nodesPerCell). Ragged cell
 *   blocks cannot be constructed through this API at all (there is no way to
 *   express them in the flat typed-array shape) -- mirrors
 *   `py_to_mesh`'s `allow_ragged=false` default.
 */
Mesh val_to_mesh(const val& rObj) {
    Mesh mesh;
    val points_val = rObj["points"];
    auto dim = rObj["dim"].as<std::size_t>();
    auto npts_len = points_val["length"].as<std::size_t>();
    if (dim == 0 || npts_len % dim != 0)
        throw meshioplusplus::WriteError("meshio++ (wasm): points length is not a multiple of dim");
    mesh.AssignPoints(float64_ndarray_from_val(points_val, {npts_len / dim, dim}));

    val cells = rObj["cells"];
    auto ncells = cells["length"].as<unsigned>();
    for (unsigned i = 0; i < ncells; ++i) {
        val block = cells[i];
        std::string type = block["type"].as<std::string>();
        auto nodes_per_cell = block["nodesPerCell"].as<std::size_t>();
        val data_val = block["data"];
        auto data_len = data_val["length"].as<std::size_t>();
        if (nodes_per_cell == 0 || data_len % nodes_per_cell != 0)
            throw meshioplusplus::WriteError("meshio++ (wasm): cell block '" + type +
                                             "' data length is not a multiple of nodesPerCell");
        mesh.AddCellBlock(
            type, int64_ndarray_from_val(data_val, {data_len / nodes_per_cell, nodes_per_cell}));
    }

    if (rObj.hasOwnProperty("point_data")) {
        val pd = rObj["point_data"];
        for (const std::string& name : js_object_keys(pd)) {
            val arr = pd[name];
            mesh.AddPointData(name,
                              float64_ndarray_from_val(arr, {arr["length"].as<std::size_t>()}));
        }
    }
    if (rObj.hasOwnProperty("cell_data")) {
        val cd = rObj["cell_data"];
        for (const std::string& name : js_object_keys(cd)) {
            val blocks_val = cd[name];
            auto nb = blocks_val["length"].as<unsigned>();
            std::vector<NDArray> blocks;
            blocks.reserve(nb);
            for (unsigned b = 0; b < nb; ++b) {
                val arr = blocks_val[b];
                blocks.push_back(float64_ndarray_from_val(arr, {arr["length"].as<std::size_t>()}));
            }
            mesh.AddCellData(name, std::move(blocks));
        }
    }
    if (rObj.hasOwnProperty("field_data")) {
        val fd = rObj["field_data"];
        for (const std::string& name : js_object_keys(fd)) {
            val arr = fd[name];
            mesh.AddFieldData(name,
                              float64_ndarray_from_val(arr, {arr["length"].as<std::size_t>()}));
        }
    }
    // Named regions (see `mesh_to_val`). `dim`/`tag` default to -1, meaning
    // "unspecified", so a caller building a plain group need not supply them.
    if (rObj.hasOwnProperty("regions")) {
        val regions = rObj["regions"];
        const auto nregions = regions["length"].as<unsigned>();
        for (unsigned i = 0; i < nregions; ++i) {
            val jr = regions[i];
            meshioplusplus::Region region;
            region.mName = jr["name"].as<std::string>();
            region.mKind = meshioplusplus::region_kind_from_name(jr["kind"].as<std::string>());
            region.mDim = jr.hasOwnProperty("dim") ? jr["dim"].as<int>() : -1;
            region.mTag =
                jr.hasOwnProperty("tag") ? static_cast<std::int64_t>(jr["tag"].as<double>()) : -1;
            val ent = jr["entries"];
            const auto n = ent["length"].as<std::size_t>();
            std::vector<std::size_t> shape;
            if (region.mKind == meshioplusplus::RegionKind::Side)
                shape = {n / 2, 2};
            else
                shape = {n};
            region.mEntries = int64_ndarray_from_val(ent, std::move(shape));
            mesh.AddRegion(std::move(region));
        }
    }
    return mesh;
}

// ---------------------------------------------------------------------
// Format dispatch goes through the shared registry (registry.hpp), the
// JS-side analogue of Python's extension_to_filetypes in _helpers.py --
// bindings/_core.cpp exposes one function per format and leaves dispatch
// entirely to Python, while the flat bindings (this file and bindings/c/)
// share the C++-level tables in src/cpp/src/registry.cpp. Parameterized writers
// get a fixed default there (documented per-entry); per-call overrides are a
// possible future API addition, deliberately out of scope for v1. Under
// Emscripten the HDF5/netCDF-conditional registry entries are compiled out,
// so the WASM format set is exactly the non-HDF5/netCDF one described in the
// file-level docs above.
// ---------------------------------------------------------------------

using meshioplusplus::registry_readers;
using meshioplusplus::registry_writers;
using meshioplusplus::resolve_format;

// " (this build has no HDF5 support)" for extensions like `.med` whose format
// the registry knows but this build compiled out; "" otherwise.
std::string compiled_out_hint(const std::string& rFormat) {
    const char* dep = meshioplusplus::registry_compiled_out(rFormat);
    return dep ? " (this build has no " + std::string(dep) + " support)" : "";
}

// Throw a genuine, message-carrying JS `Error` from C++. Emscripten's
// -fwasm-exceptions support lets a C++ exception unwind out of an exported
// function without aborting the module, but the resulting JS-side value is a
// bare, message-less `WebAssembly.Exception` -- not useful to a JS caller.
// EM_ASM executes an inline JS snippet synchronously; throwing inside it is a
// plain JS throw at that call site, which propagates normally to the caller
// of the exported function (the standard documented technique for surfacing
// a readable message across the boundary).
[[noreturn]] void throw_js_error(const std::string& rMsg) {
    EM_ASM({ throw new Error(UTF8ToString($0)); }, rMsg.c_str());
    __builtin_unreachable();  // EM_ASM's throw never returns; satisfies -Wreturn-type
}

// Wraps the body of an exported function so any meshioplusplus::ReadError/
// WriteError (or other std::exception) becomes a proper JS Error instead of
// a message-less WebAssembly.Exception at the JS boundary.
template <class F>
auto with_js_errors(F&& f) -> decltype(f()) {
    try {
        return f();
    } catch (const std::exception& e) {
        throw_js_error(e.what());
    } catch (...) {
        throw_js_error("meshio++ (wasm): unknown error");
    }
}

// ---------------------------------------------------------------------
// The JS-facing API (see wasm/src/index.mjs for the ergonomic wrapper that
// defaults `format` to "" and awaits module instantiation).
// ---------------------------------------------------------------------

/**
 * @brief Read a mesh file from the Emscripten virtual filesystem.
 * @param rPath virtual FS path (write the bytes there first via `Module.FS`).
 * @param rFormat explicit format key (see `registry_extension_defaults()`), or "" to
 *   infer from `rPath`'s extension.
 * @return a plain JS mesh object (see `mesh_to_val`).
 * @throws meshioplusplus::ReadError on an unknown/unsupported format or a
 *   malformed file.
 */
/// @return a JS array of `rItems`. Built with the no-arg `val::array()` plus
/// per-index `set`, matching every other array construction in this file --
/// embind's container-taking overloads vary across Emscripten versions and this
/// form is unambiguously available.
template <class T>
val js_array_of(const std::vector<T>& rItems) {
    val out = val::array();
    for (std::size_t i = 0; i < rItems.size(); ++i)
        out.set(i, rItems[i]);
    return out;
}

/// Resolve a read format, falling back to a content sniff (read paths only).
std::string js_resolve_read_format(const std::string& rPath, const std::string& rFormat) {
    try {
        return resolve_format(rPath, rFormat);
    } catch (const meshioplusplus::ReadError&) {
        std::string sniffed = meshioplusplus::sniff_format(rPath);
        if (sniffed.empty())
            throw;
        return sniffed;
    }
}

val read_mesh(const std::string& rPath, const std::string& rFormat) {
    return with_js_errors([&]() -> val {
        const std::string fmt = js_resolve_read_format(rPath, rFormat);
        auto it = registry_readers().find(fmt);
        if (it == registry_readers().end())
            throw meshioplusplus::ReadError("meshio++ (wasm): unknown or unsupported format '" +
                                            fmt + "'" + compiled_out_hint(fmt));
        return mesh_to_val(it->second(rPath));
    });
}

/**
 * @brief Selective read: geometry only, or only the named data arrays.
 * @param rPath virtual FS path to read.
 * @param rFormat explicit format key, or "" to infer from the extension.
 * @param points_only skip every data array.
 * @param rArrays JS array of names to keep; `null`/`undefined` keeps every
 *   array, an empty array keeps none. That distinction is deliberate.
 *
 * Formats without a native selective path are read whole and filtered, so the
 * result is the same either way -- only the cost differs.
 */
val read_mesh_selective(const std::string& rPath, const std::string& rFormat, bool points_only,
                        const val& rArrays) {
    return with_js_errors([&]() -> val {
        const std::string fmt = js_resolve_read_format(rPath, rFormat);
        meshioplusplus::ReadOptions opts;
        opts.mPointsOnly = points_only;
        if (!rArrays.isNull() && !rArrays.isUndefined())
            opts.mDataArrays = emscripten::vecFromJSArray<std::string>(rArrays);
        return mesh_to_val(meshioplusplus::registry_read(rPath, fmt, opts));
    });
}

/**
 * @brief Summarize a file without loading its heavy arrays.
 * @return `{numPoints, pointDim, numCells, cellBlocks: [{type, numCells,
 *   nodesPerCell, ragged}], pointDataNames, cellDataNames, fieldDataNames,
 *   format, fellBackToFullRead, bboxMin?, bboxMax?}`.
 *
 * `bboxMin`/`bboxMax` are present only when a box was computed -- omitted
 * rather than null, so "not computed" cannot read as a box at the origin.
 * `fellBackToFullRead` says whether the summary was actually cheap.
 */
val read_metadata_js(const std::string& rPath, const std::string& rFormat) {
    return with_js_errors([&]() -> val {
        const std::string fmt = js_resolve_read_format(rPath, rFormat);
        const meshioplusplus::MeshMetadata meta =
            meshioplusplus::registry_read_metadata(rPath, fmt, meshioplusplus::ReadOptions{});

        val blocks = val::array();
        for (std::size_t i = 0; i < meta.mCellBlocks.size(); ++i) {
            val entry = val::object();
            entry.set("type", meta.mCellBlocks[i].mType);
            entry.set("numCells", static_cast<double>(meta.mCellBlocks[i].mNumCells));
            entry.set("nodesPerCell", static_cast<double>(meta.mCellBlocks[i].mNodesPerCell));
            entry.set("ragged", meta.mCellBlocks[i].mRagged);
            blocks.set(i, entry);
        }

        val out = val::object();
        out.set("numPoints", static_cast<double>(meta.mNumPoints));
        out.set("pointDim", static_cast<double>(meta.mPointDim));
        out.set("numCells", static_cast<double>(meta.NumCells()));
        out.set("cellBlocks", blocks);
        out.set("pointDataNames", js_array_of(meta.mPointDataNames));
        out.set("cellDataNames", js_array_of(meta.mCellDataNames));
        out.set("fieldDataNames", js_array_of(meta.mFieldDataNames));
        out.set("format", meta.mFormat);
        out.set("fellBackToFullRead", meta.mFellBackToFullRead);
        if (meta.mHasBBox) {
            out.set("bboxMin", js_array_of(std::vector<double>{meta.mBBoxMin[0], meta.mBBoxMin[1],
                                                               meta.mBBoxMin[2]}));
            out.set("bboxMax", js_array_of(std::vector<double>{meta.mBBoxMax[0], meta.mBBoxMax[1],
                                                               meta.mBBoxMax[2]}));
        }
        return out;
    });
}

/// @return whether `rFormat` has a native selective-read path.
bool reader_supports_options_js(const std::string& rFormat) {
    return meshioplusplus::registry_reader_supports_options(rFormat);
}

/**
 * @brief Write a mesh object to the Emscripten virtual filesystem.
 * @param rPath virtual FS path to write (read the bytes back out via
 *   `Module.FS` afterward).
 * @param rMeshObj a plain JS mesh object (see `mesh_to_val`'s shape).
 * @param rFormat explicit format key, or "" to infer from `rPath`'s extension.
 * @throws meshioplusplus::WriteError on an unknown/write-unsupported format
 *   or malformed input.
 */
void write_mesh(const std::string& rPath, const val& rMeshObj, const std::string& rFormat) {
    with_js_errors([&]() {
        std::string fmt = resolve_format(rPath, rFormat);
        auto it = registry_writers().find(fmt);
        if (it == registry_writers().end())
            throw meshioplusplus::WriteError(
                "meshio++ (wasm): unknown, read-only, or unsupported format '" + fmt + "'" +
                compiled_out_hint(fmt));
        it->second(rPath, val_to_mesh(rMeshObj));
    });
}

/**
 * @brief Read `inPath` and immediately write it to `outPath` (both on the
 * virtual FS), without round-tripping through a JS object. Mirrors the CLI's
 * `convert` subcommand.
 */
void convert(const std::string& rInPath, const std::string& rInFormat, const std::string& rOutPath,
             const std::string& rOutFormat) {
    with_js_errors([&]() {
        std::string rfmt = resolve_format(rInPath, rInFormat);
        std::string wfmt = resolve_format(rOutPath, rOutFormat);
        auto rit = registry_readers().find(rfmt);
        auto wit = registry_writers().find(wfmt);
        if (rit == registry_readers().end())
            throw meshioplusplus::ReadError(
                "meshio++ (wasm): unknown or unsupported input format '" + rfmt + "'" +
                compiled_out_hint(rfmt));
        if (wit == registry_writers().end())
            throw meshioplusplus::WriteError(
                "meshio++ (wasm): unknown, read-only, or unsupported output format '" + wfmt + "'" +
                compiled_out_hint(wfmt));
        wit->second(rOutPath, rit->second(rInPath));
    });
}

/**
 * @brief Read `rInPath` and write a *renderable surface* of it to `rOutPath`,
 * without round-tripping through a JS object.
 *
 * The one operation a browser viewer cannot express with the existing
 * bindings. `readMesh` -> `extractSkin` -> `writeMesh` would do the same
 * thing, but the JS mesh representation is flat and drops multi-component
 * (vector/tensor) data, so every such array would be lost on the way to the
 * renderer. Staying inside C++ keeps them.
 *
 * A mesh with skinnable 3D cells becomes its boundary; anything else (a
 * surface, a curve, a point cloud) passes through. Either way the result is
 * linearized, because the target is a triangle/polygon renderer with no
 * concept of a mid-side node -- drawing `triangle6` connectivity verbatim
 * produces visible garbage rather than a curved triangle.
 */
void convert_surface(const std::string& rInPath, const std::string& rInFormat,
                     const std::string& rOutPath, const std::string& rOutFormat) {
    with_js_errors([&]() {
        std::string rfmt = resolve_format(rInPath, rInFormat);
        std::string wfmt = resolve_format(rOutPath, rOutFormat);
        auto rit = registry_readers().find(rfmt);
        auto wit = registry_writers().find(wfmt);
        if (rit == registry_readers().end())
            throw meshioplusplus::ReadError(
                "meshio++ (wasm): unknown or unsupported input format '" + rfmt + "'" +
                compiled_out_hint(rfmt));
        if (wit == registry_writers().end())
            throw meshioplusplus::WriteError(
                "meshio++ (wasm): unknown, read-only, or unsupported output format '" + wfmt + "'" +
                compiled_out_hint(wfmt));

        Mesh mesh = rit->second(rInPath);
        meshioplusplus::ConvertCellsOptions linearize;
        linearize.mMode = meshioplusplus::ConvertCellsMode::Linearize;

        if (meshioplusplus::has_skinnable_cells(mesh)) {
            // record_parent_ids, then gather: extract_surface drops cell data,
            // and colouring a solid by its per-cell tag is the common case.
            Mesh surface = meshioplusplus::extract_surface(mesh, /*recordParentIds=*/true);
            gather_cell_data_onto_surface(mesh, surface);
            // Linearize after gathering; convert_cells carries cell data 1:1.
            surface = meshioplusplus::convert_cells(surface, linearize).mMesh;
            // The parent ids are plumbing, not something to colour by.
            surface = meshioplusplus::data_drop(surface, meshioplusplus::DataLocation::Cell,
                                                {"surface:parent_cell"},
                                                /*ignore_missing=*/true);
            wit->second(rOutPath, surface);
            return;
        }
        wit->second(rOutPath, meshioplusplus::convert_cells(mesh, linearize).mMesh);
    });
}

namespace {

/// Read a 3-vector out of a JS array.
void read_vec3(const val& rArray, double* pOut) {
    for (unsigned i = 0; i < 3; ++i)
        pOut[i] = rArray[i].as<double>();
}

/// Cells across every block. The uniform mesh API counts per block only.
std::size_t total_cells(const Mesh& rMesh) {
    std::size_t n = 0;
    for (const auto& block : rMesh.CellRange())
        n += block.NumCells();
    return n;
}

/**
 * @brief Apply one operation, appending its counters to `rSteps`.
 *
 * Every operation runs on the **full-dimensional** mesh, before any boundary
 * extraction -- smoothing a solid's skin and smoothing a solid are different
 * things, and the second is what the user asked for.
 */
Mesh apply_one_op(Mesh mesh, const val& rSpec, val& rSteps, val& rWarnings) {
    const std::string op = rSpec["op"].as<std::string>();
    val step = val::object();
    step.set("op", op);

    auto number = [&rSpec](const char* key, double fallback) {
        val v = rSpec[key];
        return v.isUndefined() || v.isNull() ? fallback : v.as<double>();
    };
    auto flag = [&rSpec](const char* key, bool fallback) {
        val v = rSpec[key];
        return v.isUndefined() || v.isNull() ? fallback : v.as<bool>();
    };
    auto text = [&rSpec](const char* key, const char* fallback) {
        val v = rSpec[key];
        return v.isUndefined() || v.isNull() ? std::string(fallback) : v.as<std::string>();
    };

    if (op == "quality") {
        Mesh out = meshioplusplus::attach_quality(mesh);
        rSteps.call<void>("push", step);
        return out;
    }
    if (op == "clean") {
        meshioplusplus::CleanOptions opts;
        opts.weld = flag("weld", false);
        opts.atol = number("atol", 1e-8);
        opts.remove_orphans = flag("removeOrphans", true);
        opts.drop_degenerate = flag("dropDegenerate", true);
        opts.drop_duplicate_cells = flag("dropDuplicateCells", true);
        auto result = meshioplusplus::clean(mesh, opts);
        step.set("pointsWelded", static_cast<double>(result.mPointsWelded));
        step.set("pointsRemovedOrphan", static_cast<double>(result.mPointsRemovedOrphan));
        step.set("cellsDroppedDegenerate", static_cast<double>(result.mCellsDroppedDegenerate));
        step.set("cellsDroppedDuplicate", static_cast<double>(result.mCellsDroppedDuplicate));
        rSteps.call<void>("push", step);
        return std::move(result.mMesh);
    }
    if (op == "smooth") {
        meshioplusplus::SmoothOptions opts;
        opts.mMethod = meshioplusplus::smooth_method_from_name(text("method", "taubin"));
        opts.mIterations = static_cast<int>(number("iterations", 10));
        opts.mLambda = number("lambda", -1.0);
        opts.mMu = number("mu", -0.34);
        opts.mFixBoundary = flag("fixBoundary", true);
        auto result = meshioplusplus::smooth(mesh, opts);
        step.set("numNodesMoved", static_cast<double>(result.mNumNodesMoved));
        step.set("maxDisplacement", result.mMaxDisplacement);
        step.set("numSkippedInversion", static_cast<double>(result.mNumSkippedInversion));
        rSteps.call<void>("push", step);
        return std::move(result.mMesh);
    }
    if (op == "refine") {
        meshioplusplus::RefineOptions opts;
        opts.mLevels = static_cast<int>(number("levels", 1));
        auto result = meshioplusplus::refine(mesh, opts);
        rSteps.call<void>("push", step);
        return std::move(result.mMesh);
    }
    if (op == "decimate") {
        meshioplusplus::DecimateOptions opts;
        opts.mTargetRatio = number("ratio", -1.0);
        const double faces = number("targetFaces", -1.0);
        opts.mTargetFaces =
            faces < 0.0 ? static_cast<std::int64_t>(-1) : static_cast<std::int64_t>(faces);
        opts.mMaxError = number("maxError", -1.0);
        if (opts.mTargetRatio < 0.0 && opts.mTargetFaces < 0 && opts.mMaxError < 0.0)
            opts.mTargetRatio = 0.5;  // a usable pipeline-chip default
        opts.mPlacement =
            meshioplusplus::decimate_placement_from_name(text("placement", "optimal"));
        opts.mPreserveBoundary = flag("preserveBoundary", true);
        opts.mPreserveFeatures = flag("preserveFeatures", true);
        opts.mFeatureAngleDeg = number("featureAngle", 30.0);
        auto result = meshioplusplus::decimate(mesh, opts);
        step.set("facesRemoved", static_cast<double>(result.mFacesRemoved));
        step.set("collapsesRejected", static_cast<double>(result.mCollapsesRejected));
        rSteps.call<void>("push", step);
        return std::move(result.mMesh);
    }
    if (op == "partition") {
        meshioplusplus::PartitionOptions opts;
        opts.mNParts = static_cast<int>(number("nparts", 2));
        opts.mMethod = meshioplusplus::partition_method_from_name(text("method", "auto"));
        // Attach the assignment as cell data rather than splitting into
        // pieces: the viewer wants one mesh it can colour by part.
        auto labels = meshioplusplus::partition_labels(mesh, opts);
        mesh.AddCellData("partition:part", std::move(labels));
        step.set("nparts", number("nparts", 2));
        rSteps.call<void>("push", step);
        return mesh;
    }
    if (op == "section") {
        // The planar cross-section: slice() intersects the mesh with the plane
        // and returns a surface (triangle/quad) one dimension lower. Because a
        // surface is not skinnable, convert_surface_ops' shared re-skin tail
        // linearizes and writes it directly -- exactly the section is rendered,
        // not the boundary of a retained half (the old crop_halfspace path).
        meshioplusplus::SliceOptions opts;
        read_vec3(rSpec["point"], opts.mOrigin.data());
        read_vec3(rSpec["normal"], opts.mNormal.data());
        Mesh section = meshioplusplus::slice(mesh, opts);
        step.set("sectionFaces", static_cast<double>(total_cells(section)));
        rSteps.call<void>("push", step);
        if (total_cells(section) == 0)
            rWarnings.call<void>("push", std::string("the section is empty; the plane "
                                                     "misses the mesh -- move or flip it"));
        return section;
    }
    if (op == "isosurface") {
        // The level set of a scalar point_data field -- slice's data-driven
        // sibling, and like it a surface one dimension lower, so the shared
        // re-skin tail linearizes and writes it directly.
        meshioplusplus::IsosurfaceOptions opts;
        opts.mArrayName = text("array", "");
        opts.mIsovalues.push_back(number("isovalue", 0.0));
        const int comp = static_cast<int>(number("component", -1.0));
        if (comp >= 0)
            opts.mComponent = comp;
        Mesh contour = meshioplusplus::isosurface(mesh, opts);
        step.set("contourCells", static_cast<double>(total_cells(contour)));
        rSteps.call<void>("push", step);
        if (total_cells(contour) == 0)
            rWarnings.call<void>("push", std::string("the contour is empty; the isovalue "
                                                     "lies outside the field's range"));
        return contour;
    }
    throw meshioplusplus::ReadError("meshio++ (wasm): unknown operation '" + op + "'");
}

}  // namespace

/**
 * @brief Read `rInPath`, apply an operation pipeline, and write a renderable
 * surface to `rOutPath` -- all inside C++.
 *
 * This exists because every mesh operation in the JS API takes and returns a
 * JS `Mesh`, whose flat representation cannot carry multi-component
 * (vector/tensor) arrays. Chaining operations through that API would silently
 * destroy exactly the data `convertSurface` goes out of its way to preserve.
 * Here no mesh ever crosses the boundary, so nothing is lost.
 *
 * An **empty** pipeline is exactly `convertSurface`, which is deliberate: a
 * viewer can call this for both the plain display and the post-operation
 * display, so the two cannot drift apart. It is also what makes undo exact --
 * replaying a shortened pipeline from the original file needs no inverse
 * operations and no snapshots.
 *
 * @param rOps JS array of `{op, ...params}`; see the `OpSpec` union in
 *   `wasm/index.d.ts`.
 * @param keepProvenance keep `surface:parent_cell` in the output (the picker
 *   needs it; the colour-by menu must filter it out).
 * @return `{steps: [{op, ...counters}], warnings: [string]}`.
 */
val convert_surface_ops(const std::string& rInPath, const std::string& rInFormat,
                        const std::string& rOutPath, const std::string& rOutFormat, const val& rOps,
                        bool keepProvenance) {
    return with_js_errors([&]() -> val {
        std::string rfmt = resolve_format(rInPath, rInFormat);
        std::string wfmt = resolve_format(rOutPath, rOutFormat);
        auto rit = registry_readers().find(rfmt);
        auto wit = registry_writers().find(wfmt);
        if (rit == registry_readers().end())
            throw meshioplusplus::ReadError(
                "meshio++ (wasm): unknown or unsupported input format '" + rfmt + "'" +
                compiled_out_hint(rfmt));
        if (wit == registry_writers().end())
            throw meshioplusplus::WriteError(
                "meshio++ (wasm): unknown, read-only, or unsupported output format '" + wfmt + "'" +
                compiled_out_hint(wfmt));

        val steps = val::array();
        val warnings = val::array();

        Mesh mesh = rit->second(rInPath);
        const unsigned n = rOps["length"].as<unsigned>();
        for (unsigned i = 0; i < n; ++i)
            mesh = apply_one_op(std::move(mesh), rOps[i], steps, warnings);

        // From here on this is `convert_surface`'s tail, unchanged.
        meshioplusplus::ConvertCellsOptions linearize;
        linearize.mMode = meshioplusplus::ConvertCellsMode::Linearize;

        if (meshioplusplus::has_skinnable_cells(mesh)) {
            Mesh surface = meshioplusplus::extract_surface(mesh, /*recordParentIds=*/true);
            gather_cell_data_onto_surface(mesh, surface);
            surface = meshioplusplus::convert_cells(surface, linearize).mMesh;
            if (!keepProvenance)
                surface = meshioplusplus::data_drop(surface, meshioplusplus::DataLocation::Cell,
                                                    {"surface:parent_cell"},
                                                    /*ignore_missing=*/true);
            wit->second(rOutPath, surface);
        } else {
            wit->second(rOutPath, meshioplusplus::convert_cells(mesh, linearize).mMesh);
        }

        val out = val::object();
        out.set("steps", steps);
        out.set("warnings", warnings);
        return out;
    });
}

/** @brief The shared `num_nodes_per_cell` metadata table, as a plain JS object. */
val num_nodes_per_cell_js() {
    val out = val::object();
    for (const auto& kv : meshioplusplus::num_nodes_per_cell())
        out.set(kv.first, kv.second);
    return out;
}

/** @brief The shared `topological_dimension` metadata table, as a plain JS object. */
val topological_dimension_js() {
    val out = val::object();
    for (const auto& kv : meshioplusplus::topological_dimension())
        out.set(kv.first, kv.second);
    return out;
}

/** @brief The compile-time mesh backend ("native" for the shipped wasm build). */
std::string mesh_backend_js() {
    return meshioplusplus::mesh_backend_name();
}

/**
 * @brief The format names this build can actually read and write.
 *
 * `{readers: [...], writers: [...]}`, both sorted (the registry tables are
 * `std::map`). This is what a caller needs to build a file-picker filter or a
 * "convert to" menu without hardcoding a table that silently drifts from the
 * build -- the HDF5/netCDF-backed formats are absent from the wasm build, and
 * a few formats are read-only or write-only.
 */
val available_formats_js() {
    val readers = val::array();
    val writers = val::array();
    unsigned i = 0;
    for (const auto& kv : registry_readers())
        readers.set(i++, kv.first);
    i = 0;
    for (const auto& kv : registry_writers())
        writers.set(i++, kv.first);
    val out = val::object();
    out.set("readers", readers);
    out.set("writers", writers);
    return out;
}

// ---------------------------------------------------------------------
// Mesh operations (computations on a mesh, not file formats).
// ---------------------------------------------------------------------

/** @brief Extract the boundary (volume -> faces, surface -> edges). */
val extract_surface_js(const val& rMeshObj, bool recordParentIds) {
    return with_js_errors([&]() -> val {
        return mesh_to_val(meshioplusplus::extract_surface(val_to_mesh(rMeshObj), recordParentIds));
    });
}

/** @brief Extract the boundary skin of a volume mesh. */
val extract_skin_js(const val& rMeshObj, bool linearize) {
    return with_js_errors([&]() -> val {
        return mesh_to_val(meshioplusplus::extract_skin(val_to_mesh(rMeshObj), linearize));
    });
}

/** @brief Attach per-cell quality metrics as cell_data. */
val attach_quality_js(const val& rMeshObj) {
    return with_js_errors([&]() -> val {
        return mesh_to_val(meshioplusplus::attach_quality(val_to_mesh(rMeshObj)));
    });
}

/** @brief Guess a mesh file's format from its contents ("" if undetermined). */
std::string sniff_format_js(const std::string& rPath) {
    return with_js_errors([&]() -> std::string { return meshioplusplus::sniff_format(rPath); });
}

/**
 * @brief Renumber a mesh (method: "rcm", "morton", or "hilbert"). Returns an
 * object `{mesh, nodePermutation, cellPermutations}` (permutations old->new).
 */
val reorder_js(const val& rMeshObj, const std::string& rMethod) {
    return with_js_errors([&]() -> val {
        meshioplusplus::ReorderResult res = meshioplusplus::reorder(
            val_to_mesh(rMeshObj), meshioplusplus::reorder_method_from_name(rMethod));
        val out = val::object();
        out.set("mesh", mesh_to_val(res.mMesh));
        out.set("nodePermutation", ndarray_to_int32_array(res.mNodePermutation));
        val cell_perms = val::array();
        for (const NDArray& a : res.mCellPermutations)
            cell_perms.call<void>("push", ndarray_to_int32_array(a));
        out.set("cellPermutations", cell_perms);
        return out;
    });
}

/**
 * @brief Merge a JS array of mesh objects into one (concatenate, optional
 * welding). `dataPolicy` is "intersection" (default) or "fill". Returns a plain
 * JS mesh object (point_sets/cell_sets are not carried, as elsewhere in JS).
 */
val merge_js(const val& rMeshes, bool weld, double atol, bool sourceTag,
             const std::string& rDataPolicy, bool dropDuplicateCells) {
    return with_js_errors([&]() -> val {
        const unsigned n = rMeshes["length"].as<unsigned>();
        std::vector<Mesh> owned;
        owned.reserve(n);
        for (unsigned i = 0; i < n; ++i)
            owned.push_back(val_to_mesh(rMeshes[i]));
        std::vector<const Mesh*> ptrs;
        ptrs.reserve(owned.size());
        for (const Mesh& mm : owned)
            ptrs.push_back(&mm);
        meshioplusplus::MergeOptions opts;
        opts.weld = weld;
        opts.atol = atol;
        opts.source_tag = sourceTag;
        opts.drop_duplicate_cells = dropDuplicateCells;
        opts.data_policy = (rDataPolicy == "fill") ? meshioplusplus::MergeDataPolicy::Fill
                                                   : meshioplusplus::MergeDataPolicy::Intersection;
        return mesh_to_val(meshioplusplus::merge(ptrs, opts).mMesh);
    });
}

/**
 * @brief Apply an affine transform to a mesh's point coordinates. `matrix` is a
 * JS array of 16 numbers (row-major 4x4). Returns the transformed mesh object.
 */
val transform_js(const val& rMeshObj, const val& rMatrix, bool rotateVectorData) {
    return with_js_errors([&]() -> val {
        const unsigned n = rMatrix["length"].as<unsigned>();
        if (n != 16)
            throw meshioplusplus::WriteError("transform: matrix must have 16 elements (4x4)");
        double m[16];
        for (unsigned i = 0; i < 16; ++i)
            m[i] = rMatrix[i].as<double>();
        meshioplusplus::AffineTransform xf = meshioplusplus::transform_from_matrix(m);
        return mesh_to_val(meshioplusplus::transform(val_to_mesh(rMeshObj), xf, rotateVectorData));
    });
}

/**
 * @brief Clean a mesh (weld / prune / de-dup). Returns an object
 * `{mesh, pointsWelded, pointsRemovedOrphan, cellsDroppedDegenerate,
 * cellsDroppedDuplicate}` (sets/maps are not carried, as elsewhere in JS).
 */
val clean_js(const val& rMeshObj, bool weld, double atol, bool removeOrphans, bool dropDegenerate,
             bool dropDuplicateCells) {
    return with_js_errors([&]() -> val {
        meshioplusplus::CleanOptions opts;
        opts.weld = weld;
        opts.atol = atol;
        opts.remove_orphans = removeOrphans;
        opts.drop_degenerate = dropDegenerate;
        opts.drop_duplicate_cells = dropDuplicateCells;
        meshioplusplus::CleanResult r = meshioplusplus::clean(val_to_mesh(rMeshObj), opts);
        val out = val::object();
        out.set("mesh", mesh_to_val(r.mMesh));
        out.set("pointsWelded", static_cast<double>(r.mPointsWelded));
        out.set("pointsRemovedOrphan", static_cast<double>(r.mPointsRemovedOrphan));
        out.set("cellsDroppedDegenerate", static_cast<double>(r.mCellsDroppedDegenerate));
        out.set("cellsDroppedDuplicate", static_cast<double>(r.mCellsDroppedDuplicate));
        return out;
    });
}

/**
 * @brief Smooth a mesh's point coordinates (`"laplacian"` / `"taubin"`),
 * leaving connectivity and every data array untouched. A negative `lambda`
 * means "this method's own default" (0.5 Laplacian, 0.33 Taubin) and is passed
 * through unchanged. Returns an object `{mesh, numNodesMoved, maxDisplacement,
 * numSkippedInversion}`. The caller-supplied frozen-node mask is not exposed
 * here, as on the other flat bindings.
 */
val smooth_js(const val& rMeshObj, const std::string& rMethod, int iterations, double lambda,
              double mu, bool fixBoundary, bool preserveFeatures, double featureAngle,
              bool guardInversion) {
    return with_js_errors([&]() -> val {
        meshioplusplus::SmoothOptions options;
        options.mMethod = meshioplusplus::smooth_method_from_name(rMethod);
        options.mIterations = iterations;
        options.mLambda = lambda;
        options.mMu = mu;
        options.mFixBoundary = fixBoundary;
        options.mPreserveFeatures = preserveFeatures;
        options.mFeatureAngleDeg = featureAngle;
        options.mGuardInversion = guardInversion;
        meshioplusplus::SmoothResult r = meshioplusplus::smooth(val_to_mesh(rMeshObj), options);
        val out = val::object();
        out.set("mesh", mesh_to_val(r.mMesh));
        out.set("numNodesMoved", static_cast<double>(r.mNumNodesMoved));
        out.set("maxDisplacement", r.mMaxDisplacement);
        out.set("numSkippedInversion", static_cast<double>(r.mNumSkippedInversion));
        return out;
    });
}

/**
 * @brief Crop a mesh to a bounding box (`lo`/`hi` are 3-element JS arrays).
 * `mode` is "all" (default) or "any". Returns the pruned mesh object.
 */
val crop_bbox_js(const val& rMeshObj, const val& rLo, const val& rHi, const std::string& rMode,
                 bool recordIds) {
    return with_js_errors([&]() -> val {
        double lo[3], hi[3];
        for (unsigned i = 0; i < 3; ++i) {
            lo[i] = rLo[i].as<double>();
            hi[i] = rHi[i].as<double>();
        }
        meshioplusplus::CropMode m =
            (rMode == "any") ? meshioplusplus::CropMode::Any : meshioplusplus::CropMode::All;
        return mesh_to_val(
            meshioplusplus::crop_bbox(val_to_mesh(rMeshObj), lo, hi, m, recordIds).mMesh);
    });
}

/**
 * @brief Crop a mesh to the half-space (p - point) . normal >= 0 (`point`/
 * `normal` are 3-element JS arrays). Returns the pruned mesh object.
 */
val crop_plane_js(const val& rMeshObj, const val& rPoint, const val& rNormal,
                  const std::string& rMode, bool recordIds) {
    return with_js_errors([&]() -> val {
        double point[3], normal[3];
        for (unsigned i = 0; i < 3; ++i) {
            point[i] = rPoint[i].as<double>();
            normal[i] = rNormal[i].as<double>();
        }
        meshioplusplus::CropMode m =
            (rMode == "any") ? meshioplusplus::CropMode::Any : meshioplusplus::CropMode::All;
        return mesh_to_val(
            meshioplusplus::crop_halfspace(val_to_mesh(rMeshObj), point, normal, m, recordIds)
                .mMesh);
    });
}

/**
 * @brief Planar cross-section of a mesh (marching tetrahedra on a simplexified
 * input). Returns a new mesh one dimension below the cut cells: a triangle/quad
 * surface for a volume mesh, a line mesh for a 2D surface mesh.
 */
val slice_js(const val& rMeshObj, const val& rOrigin, const val& rNormal, bool recordParentIds) {
    return with_js_errors([&]() -> val {
        meshioplusplus::SliceOptions options;
        for (unsigned i = 0; i < 3; ++i) {
            options.mOrigin[i] = rOrigin[i].as<double>();
            options.mNormal[i] = rNormal[i].as<double>();
        }
        options.mRecordParentIds = recordParentIds;
        return mesh_to_val(meshioplusplus::slice(val_to_mesh(rMeshObj), options));
    });
}

/**
 * @brief Isosurfaces / contours: the level set of a scalar point_data field,
 * cut with the same marching tetrahedra as slice. `component` is negative for
 * the row magnitude. Returns one mesh holding every contour, tagged per cell
 * with `iso:value` (Float64) and `iso:index` (Int64).
 */
val isosurface_js(const val& rMeshObj, const std::string& rArray, const val& rIsovalues,
                  int component, bool recordParentIds) {
    return with_js_errors([&]() -> val {
        meshioplusplus::IsosurfaceOptions options;
        options.mArrayName = rArray;
        const unsigned n = rIsovalues["length"].as<unsigned>();
        options.mIsovalues.reserve(n);
        for (unsigned i = 0; i < n; ++i)
            options.mIsovalues.push_back(rIsovalues[i].as<double>());
        if (component >= 0)
            options.mComponent = component;
        options.mRecordParentIds = recordParentIds;
        return mesh_to_val(meshioplusplus::isosurface(val_to_mesh(rMeshObj), options));
    });
}

/**
 * @brief Split a mesh into pieces (by "type" / "component" / "region"|"tag").
 * Returns a JS array of `{key, mesh}` objects. `tagName` selects the integer
 * cell_data for the tag criterion (empty = auto-detect).
 */
val split_js(const val& rMeshObj, const std::string& rBy, const std::string& rTagName) {
    return with_js_errors([&]() -> val {
        meshioplusplus::SplitResult r = meshioplusplus::split(
            val_to_mesh(rMeshObj), meshioplusplus::split_by_from_name(rBy), rTagName);
        val out = val::array();
        for (meshioplusplus::SplitPiece& p : r.mPieces) {
            val piece = val::object();
            piece.set("key", p.mKey);
            piece.set("mesh", mesh_to_val(p.mMesh));
            out.call<void>("push", piece);
        }
        return out;
    });
}

/**
 * @brief Convert the element representation of a mesh: `"linearize"`,
 * `"simplexify"`, or `"elevate"`. Returns the converted mesh; the index maps
 * are not carried across the JS boundary (use `recordParentIds` for the
 * `convert:parent_cell` cell_data instead).
 */
val convert_cells_js(const val& rMeshObj, const std::string& rMode, bool recordParentIds) {
    return with_js_errors([&]() -> val {
        meshioplusplus::ConvertCellsOptions options;
        options.mMode = meshioplusplus::convert_cells_mode_from_name(rMode);
        options.mRecordParentIds = recordParentIds;
        return mesh_to_val(meshioplusplus::convert_cells(val_to_mesh(rMeshObj), options).mMesh);
    });
}

/**
 * @brief Uniformly refine a mesh, subdividing every cell into same-type
 * children (line -> 2, triangle -> 4, quad -> 4, tetra -> 8, wedge -> 8,
 * hexahedron -> 8). Returns the refined mesh; the index maps are not carried
 * across the JS boundary (use `recordParentIds` for the `refine:parent_cell`
 * cell_data instead).
 */
val refine_js(const val& rMeshObj, int levels, bool recordParentIds) {
    return with_js_errors([&]() -> val {
        meshioplusplus::RefineOptions options;
        options.mLevels = levels;
        options.mRecordParentIds = recordParentIds;
        return mesh_to_val(meshioplusplus::refine(val_to_mesh(rMeshObj), options).mMesh);
    });
}

/**
 * @brief Decimate a SURFACE mesh by quadric-error-metric edge collapse — the
 * resolution-reducing inverse of `refine`. Exactly one of `ratio` (fraction of
 * faces to keep, in (0, 1]), `targetFaces` and `maxError` must be non-negative.
 * Returns an object `{mesh, facesRemoved, pointsRemoved, collapsesRejected,
 * maxErrorApplied}`; the index maps are not carried across the JS boundary and
 * the frozen mask is not exposed here, as on the other flat bindings.
 */
val decimate_js(const val& rMeshObj, double ratio, double targetFaces, double maxError,
                const std::string& rPlacement, bool preserveBoundary, bool preserveFeatures,
                double featureAngle) {
    return with_js_errors([&]() -> val {
        meshioplusplus::DecimateOptions options;
        options.mTargetRatio = ratio;
        options.mTargetFaces = targetFaces < 0.0 ? static_cast<std::int64_t>(-1)
                                                 : static_cast<std::int64_t>(targetFaces);
        options.mMaxError = maxError;
        options.mPlacement = meshioplusplus::decimate_placement_from_name(rPlacement);
        options.mPreserveBoundary = preserveBoundary;
        options.mPreserveFeatures = preserveFeatures;
        options.mFeatureAngleDeg = featureAngle;
        meshioplusplus::DecimateResult r = meshioplusplus::decimate(val_to_mesh(rMeshObj), options);
        val out = val::object();
        out.set("mesh", mesh_to_val(r.mMesh));
        out.set("facesRemoved", static_cast<double>(r.mFacesRemoved));
        out.set("pointsRemoved", static_cast<double>(r.mPointsRemoved));
        out.set("collapsesRejected", static_cast<double>(r.mCollapsesRejected));
        out.set("maxErrorApplied", r.mMaxErrorApplied);
        return out;
    });
}

namespace {

meshioplusplus::PartitionOptions partition_options_js(int nparts, const std::string& rMethod,
                                                      double imbalance, const std::string& rMode,
                                                      int seed, bool recordIds, int ghostLayers,
                                                      const std::string& rWeightsKey) {
    meshioplusplus::PartitionOptions options;
    options.mNParts = nparts;
    options.mMethod = meshioplusplus::partition_method_from_name(rMethod);
    options.mImbalance = imbalance;
    options.mMode = meshioplusplus::partition_mode_from_name(rMode);
    options.mSeed = seed;
    options.mRecordIds = recordIds;
    options.mGhostLayers = ghostLayers;
    options.mWeightsKey = rWeightsKey;
    return options;
}

}  // namespace

/**
 * @brief Decompose a mesh into `nparts` balanced pieces (`"sfc"` / `"kahip"` /
 * `"auto"`). Returns a JS array of `{partId, mesh}` objects, exactly `nparts`
 * entries (pieces may be empty; blocks kept 1:1 with the input, unlike
 * `split`). The index maps are not carried across the JS boundary (use
 * `recordIds` for the `partition:original_*_id` arrays, or `partitionLabels`
 * for the assignment). KaHIP is never compiled into the WASM build, so
 * `method: "kahip"` throws the error naming `MESHIOPLUSPLUS_WITH_KAHIP`.
 */
val partition_js(const val& rMeshObj, int nparts, const std::string& rMethod, double imbalance,
                 const std::string& rMode, int seed, bool recordIds, int ghostLayers,
                 const std::string& rWeightsKey) {
    return with_js_errors([&]() -> val {
        meshioplusplus::PartitionResult r = meshioplusplus::partition(
            val_to_mesh(rMeshObj), partition_options_js(nparts, rMethod, imbalance, rMode, seed,
                                                        recordIds, ghostLayers, rWeightsKey));
        val out = val::array();
        for (meshioplusplus::PartitionPiece& p : r.mPieces) {
            val piece = val::object();
            piece.set("partId", p.mPartId);
            piece.set("mesh", mesh_to_val(p.mMesh));
            out.call<void>("push", piece);
        }
        return out;
    });
}

/**
 * @brief The per-cell part assignment only: a JS array with one Int32Array-like
 * array per cell block (block-aligned, like a cell_data entry), values in
 * `[0, nparts)`.
 */
val partition_labels_js(const val& rMeshObj, int nparts, const std::string& rMethod,
                        double imbalance, const std::string& rMode, int seed,
                        const std::string& rWeightsKey) {
    return with_js_errors([&]() -> val {
        std::vector<NDArray> labels = meshioplusplus::partition_labels(
            val_to_mesh(rMeshObj),
            partition_options_js(nparts, rMethod, imbalance, rMode, seed,
                                 /*recordIds=*/false, /*ghostLayers=*/0, rWeightsKey));
        val out = val::array();
        for (const NDArray& a : labels)
            out.call<void>("push", ndarray_to_int32_array(a));
        return out;
    });
}

/**
 * @brief Geometric statistics of a mesh (read-only). Returns an object with
 * `numPoints`, `numCells`, `bboxMin`/`bboxMax`/`extent`/`centroid` (3-arrays),
 * `cellTypeCounts` (object), `totalArea`, `signedVolume`, `unsignedVolume`,
 * `numInverted`.
 */
val stats_js(const val& rMeshObj) {
    return with_js_errors([&]() -> val {
        meshioplusplus::StatsReport r = meshioplusplus::compute_stats(val_to_mesh(rMeshObj));
        auto vec3 = [](const double* p) {
            val a = val::array();
            for (int i = 0; i < 3; ++i)
                a.call<void>("push", p[i]);
            return a;
        };
        val out = val::object();
        out.set("numPoints", static_cast<double>(r.mNumPoints));
        out.set("numCells", static_cast<double>(r.mNumCells));
        out.set("bboxMin", vec3(r.mBBoxMin));
        out.set("bboxMax", vec3(r.mBBoxMax));
        out.set("extent", vec3(r.mExtent));
        out.set("centroid", vec3(r.mCentroid));
        val counts = val::object();
        for (const auto& kv : r.mCellTypeCounts)
            counts.set(kv.first, static_cast<double>(kv.second));
        out.set("cellTypeCounts", counts);
        out.set("totalArea", r.mTotalArea);
        out.set("signedVolume", r.mSignedVolume);
        out.set("unsignedVolume", r.mUnsignedVolume);
        out.set("numInverted", static_cast<double>(r.mNumInverted));
        return out;
    });
}

/** @brief Connectivity bandwidth (max |i - j| over node pairs sharing a cell). */
int compute_bandwidth_js(const val& rMeshObj) {
    return with_js_errors([&]() -> int {
        return static_cast<int>(meshioplusplus::compute_bandwidth(val_to_mesh(rMeshObj)));
    });
}

// Convert an ArrayDiff to a JS object.
val array_diff_to_val(const meshioplusplus::ArrayDiff& rAd) {
    val d = val::object();
    d.set("name", rAd.mName);
    d.set("shapeMismatch", rAd.mShapeMismatch);
    d.set("sizeA", static_cast<double>(rAd.mSizeA));
    d.set("sizeB", static_cast<double>(rAd.mSizeB));
    d.set("maxAbsError", rAd.mMaxAbsError);
    d.set("maxRelError", rAd.mMaxRelError);
    d.set("worstIndex", static_cast<double>(rAd.mWorstIndex));
    d.set("numExceeding", static_cast<double>(rAd.mNumExceeding));
    d.set("exact", rAd.mExact);
    return d;
}

val data_diff_to_val(const meshioplusplus::DataDiff& rDd) {
    val d = val::object();
    val only_a = val::array(), only_b = val::array(), shared = val::array();
    for (const std::string& s : rDd.mOnlyInA)
        only_a.call<void>("push", s);
    for (const std::string& s : rDd.mOnlyInB)
        only_b.call<void>("push", s);
    for (const meshioplusplus::ArrayDiff& ad : rDd.mShared)
        shared.call<void>("push", array_diff_to_val(ad));
    d.set("onlyInA", only_a);
    d.set("onlyInB", only_b);
    d.set("shared", shared);
    return d;
}

/**
 * @brief Compare two meshes. Returns a report object with `verdict`
 * ("identical" / "equal within tolerance" / "different"), `points`, `blocks`,
 * `pointData`/`cellData`/`fieldData`, and `messages`. Named point/cell sets are
 * not compared (not visible to the core).
 */
val diff_js(const val& rMeshA, const val& rMeshB, double atol, double rtol, bool unordered) {
    return with_js_errors([&]() -> val {
        meshioplusplus::DiffOptions opts;
        opts.atol = atol;
        opts.rtol = rtol;
        opts.unordered = unordered;
        meshioplusplus::DiffReport rep =
            meshioplusplus::diff(val_to_mesh(rMeshA), val_to_mesh(rMeshB), opts);
        val out = val::object();
        out.set("verdict", meshioplusplus::diff_verdict_name(rep.mVerdict));
        out.set("unordered", rep.mUnordered);
        out.set("correspondenceFailed", rep.mCorrespondenceFailed);
        out.set("pointCountMismatch", rep.mPointCountMismatch);
        out.set("points", array_diff_to_val(rep.mPoints));
        out.set("blockCountMismatch", rep.mBlockCountMismatch);
        val blocks = val::array();
        for (const meshioplusplus::BlockDiff& bd : rep.mBlocks) {
            val d = val::object();
            d.set("block", static_cast<double>(bd.mBlock));
            d.set("typeA", bd.mTypeA);
            d.set("typeB", bd.mTypeB);
            d.set("countA", static_cast<double>(bd.mCountA));
            d.set("countB", static_cast<double>(bd.mCountB));
            d.set("typeMismatch", bd.mTypeMismatch);
            d.set("countMismatch", bd.mCountMismatch);
            d.set("connMismatchCount", static_cast<double>(bd.mConnMismatchCount));
            blocks.call<void>("push", d);
        }
        out.set("blocks", blocks);
        out.set("pointData", data_diff_to_val(rep.mPointData));
        out.set("cellData", data_diff_to_val(rep.mCellData));
        out.set("fieldData", data_diff_to_val(rep.mFieldData));
        val messages = val::array();
        for (const std::string& s : rep.mMessages)
            messages.call<void>("push", s);
        out.set("messages", messages);
        return out;
    });
}

/** @brief Are two meshes equal within tolerance? (verdict != "different") */
bool meshes_equal_js(const val& rMeshA, const val& rMeshB, double atol, double rtol,
                     bool unordered) {
    return with_js_errors([&]() -> bool {
        meshioplusplus::DiffOptions opts;
        opts.atol = atol;
        opts.rtol = rtol;
        opts.unordered = unordered;
        return meshioplusplus::diff(val_to_mesh(rMeshA), val_to_mesh(rMeshB), opts).mVerdict !=
               meshioplusplus::DiffVerdict::Different;
    });
}

// --- data operations -------------------------------------------------------
// These act on a mesh's point/cell/field data arrays; the geometry is never
// modified. Enumerations cross as strings (as `split`'s `by` already does),
// parsed by the same helpers the CLIs use.

/// A JS array of strings -> std::vector. This file deliberately avoids
/// `register_vector`, so the decode is manual.
std::vector<std::string> val_to_string_vector(const val& rArray) {
    std::vector<std::string> out;
    if (rArray.isUndefined() || rArray.isNull())
        return out;
    const unsigned n = rArray["length"].as<unsigned>();
    out.reserve(n);
    for (unsigned i = 0; i < n; ++i)
        out.push_back(rArray[i].as<std::string>());
    return out;
}

/** @brief Drop the named data arrays at `location` ("point"/"cell"/"field"). */
val data_drop_js(const val& rMeshObj, const std::string& rLocation, const val& rNames,
                 bool ignoreMissing) {
    return with_js_errors([&]() -> val {
        return mesh_to_val(meshioplusplus::data_drop(
            val_to_mesh(rMeshObj), meshioplusplus::data_location_from_name(rLocation),
            val_to_string_vector(rNames), ignoreMissing));
    });
}

/** @brief Keep only the named data arrays at `location`, dropping the rest
 *  there. The other locations are untouched. */
val data_keep_js(const val& rMeshObj, const std::string& rLocation, const val& rNames,
                 bool ignoreMissing) {
    return with_js_errors([&]() -> val {
        return mesh_to_val(meshioplusplus::data_keep(
            val_to_mesh(rMeshObj), meshioplusplus::data_location_from_name(rLocation),
            val_to_string_vector(rNames), ignoreMissing));
    });
}

/** @brief Rename one data array, preserving values, dtype and shape. */
val data_rename_js(const val& rMeshObj, const std::string& rLocation, const std::string& rFrom,
                   const std::string& rTo) {
    return with_js_errors([&]() -> val {
        return mesh_to_val(meshioplusplus::data_rename(
            val_to_mesh(rMeshObj), meshioplusplus::data_location_from_name(rLocation), rFrom, rTo));
    });
}

/** @brief Average point_data onto the cells (mean over each cell's nodes).
 *  An empty `names` converts every point_data array. Output is Float64. */
val data_point_to_cell_js(const val& rMeshObj, const val& rNames, const std::string& rSuffix) {
    return with_js_errors([&]() -> val {
        meshioplusplus::DataAverageOptions opts;
        opts.names = val_to_string_vector(rNames);
        opts.suffix = rSuffix;
        return mesh_to_val(meshioplusplus::point_data_to_cell_data(val_to_mesh(rMeshObj), opts));
    });
}

/** @brief Average cell_data onto the points. `weight` is "uniform" or
 *  "measure". Points touched by no cell get NaN. Output is Float64. */
val data_cell_to_point_js(const val& rMeshObj, const val& rNames, const std::string& rWeight,
                          const std::string& rSuffix) {
    return with_js_errors([&]() -> val {
        meshioplusplus::DataAverageOptions opts;
        opts.names = val_to_string_vector(rNames);
        opts.weight = meshioplusplus::cell_point_weight_from_name(rWeight);
        opts.suffix = rSuffix;
        return mesh_to_val(meshioplusplus::cell_data_to_point_data(val_to_mesh(rMeshObj), opts));
    });
}

/** @brief Evaluate an elementwise expression into a new array. The grammar
 *  accepts + - * /, unary minus, parentheses, numbers, array names and the
 *  functions abs/sqrt/min/max/norm -- nothing else is evaluated. */
val data_calc_js(const val& rMeshObj, const std::string& rExpression, const std::string& rLocation,
                 const std::string& rOutputName, bool overwrite) {
    return with_js_errors([&]() -> val {
        meshioplusplus::DataCalcOptions opts;
        opts.location = meshioplusplus::data_location_from_name(rLocation);
        opts.output = rOutputName;
        opts.overwrite = overwrite;
        return mesh_to_val(meshioplusplus::data_calc(val_to_mesh(rMeshObj), rExpression, opts));
    });
}

/** @brief Condition values: `mode` is "clamp"/"normalize"/"standardize",
 *  `scope` is "component"/"magnitude", `nanPolicy` is
 *  "ignore"/"replace"/"fail". */
val data_condition_js(const val& rMeshObj, const std::string& rLocation, const val& rNames,
                      const std::string& rMode, double lo, double hi, const std::string& rScope,
                      const std::string& rNanPolicy, double nanReplacement,
                      const std::string& rSuffix) {
    return with_js_errors([&]() -> val {
        meshioplusplus::DataConditionOptions opts;
        opts.location = meshioplusplus::data_location_from_name(rLocation);
        opts.names = val_to_string_vector(rNames);
        opts.mode = meshioplusplus::condition_mode_from_name(rMode);
        opts.scope = meshioplusplus::condition_scope_from_name(rScope);
        opts.lo = lo;
        opts.hi = hi;
        opts.nan_policy = meshioplusplus::nan_policy_from_name(rNanPolicy);
        opts.nan_replacement = nanReplacement;
        opts.suffix = rSuffix;
        return mesh_to_val(meshioplusplus::data_condition(val_to_mesh(rMeshObj), opts));
    });
}

/** @brief Read-only per-array data summary. Returns a JS array of objects with
 *  `location`, `name`, `dtype`, `shape`, `numBlocks`, `numEntries`,
 *  `numComponents`, `numValues`, `min`, `max`, `mean`, the three
 *  `*PerComponent` arrays, `numNan`, `numInf`, `numFinite` and
 *  `inconsistentBlocks`. */
val data_info_js(const val& rMeshObj) {
    return with_js_errors([&]() -> val {
        meshioplusplus::DataInfoReport r = meshioplusplus::data_info(val_to_mesh(rMeshObj));
        auto to_array = [](const std::vector<double>& rValues) {
            val a = val::array();
            for (double v : rValues)
                a.call<void>("push", v);
            return a;
        };
        val out = val::array();
        for (const meshioplusplus::DataArrayInfo& a : r.mArrays) {
            val o = val::object();
            o.set("location", std::string(meshioplusplus::data_location_name(a.mLocation)));
            o.set("name", a.mName);
            o.set("dtype", std::string(meshioplusplus::dtype_numpy_str(a.mDtype)));
            val shape = val::array();
            for (std::size_t d : a.mShape)
                shape.call<void>("push", static_cast<double>(d));
            o.set("shape", shape);
            o.set("numBlocks", static_cast<double>(a.mNumBlocks));
            o.set("numEntries", static_cast<double>(a.mNumEntries));
            o.set("numComponents", static_cast<double>(a.mNumComponents));
            o.set("numValues", static_cast<double>(a.mNumValues));
            o.set("min", a.mMin);
            o.set("max", a.mMax);
            o.set("mean", a.mMean);
            o.set("minPerComponent", to_array(a.mMinPerComponent));
            o.set("maxPerComponent", to_array(a.mMaxPerComponent));
            o.set("meanPerComponent", to_array(a.mMeanPerComponent));
            o.set("numNan", static_cast<double>(a.mNumNan));
            o.set("numInf", static_cast<double>(a.mNumInf));
            o.set("numFinite", static_cast<double>(a.mNumFinite));
            o.set("inconsistentBlocks", a.mInconsistentBlocks);
            out.call<void>("push", o);
        }
        return out;
    });
}

/**
 * @brief Cross-mesh field transfer: sample the source's data arrays onto the
 * target (source point_data at the target's points, source cell_data by
 * nearest source-cell centroid regardless of the method). `method` is
 * "nearest" or "barycentric"; `arrays` is a JS array of source names
 * (undefined/null/empty = every source point_data array; cell_data only when
 * named); `onConflict` is "error", "overwrite" or "suffix" (name + "_interp").
 * Returns the target copy as a plain mesh object.
 */
val interpolate_js(const val& rSourceObj, const val& rTargetObj, const std::string& rMethod,
                   const val& rArrays, bool extrapolate, double defaultValue,
                   const std::string& rOnConflict) {
    return with_js_errors([&]() -> val {
        meshioplusplus::InterpolateOptions options;
        options.mMethod = meshioplusplus::interpolate_method_from_name(rMethod);
        options.mArrays = val_to_string_vector(rArrays);
        options.mExtrapolate = extrapolate;
        options.mDefaultValue = defaultValue;
        options.mOnConflict = meshioplusplus::interpolate_conflict_from_name(rOnConflict);
        return mesh_to_val(
            meshioplusplus::interpolate(val_to_mesh(rSourceObj), val_to_mesh(rTargetObj), options));
    });
}

}  // namespace

EMSCRIPTEN_BINDINGS(meshioplusplus_wasm) {
    emscripten::function("readMesh", &read_mesh);
    emscripten::function("readMeshSelective", &read_mesh_selective);
    emscripten::function("readMetadata", &read_metadata_js);
    emscripten::function("readerSupportsOptions", &reader_supports_options_js);
    emscripten::function("writeMesh", &write_mesh);
    emscripten::function("convert", &convert);
    emscripten::function("convertSurface", &convert_surface);
    emscripten::function("convertSurfaceOps", &convert_surface_ops);
    emscripten::function("numNodesPerCell", &num_nodes_per_cell_js);
    emscripten::function("topologicalDimension", &topological_dimension_js);
    emscripten::function("meshBackend", &mesh_backend_js);
    emscripten::function("availableFormats", &available_formats_js);
    emscripten::function("extractSurface", &extract_surface_js);
    emscripten::function("extractSkin", &extract_skin_js);
    emscripten::function("attachQuality", &attach_quality_js);
    emscripten::function("sniffFormat", &sniff_format_js);
    emscripten::function("reorder", &reorder_js);
    emscripten::function("computeBandwidth", &compute_bandwidth_js);
    emscripten::function("diff", &diff_js);
    emscripten::function("meshesEqual", &meshes_equal_js);
    emscripten::function("merge", &merge_js);
    emscripten::function("transform", &transform_js);
    emscripten::function("clean", &clean_js);
    emscripten::function("smooth", &smooth_js);
    emscripten::function("interpolate", &interpolate_js);
    emscripten::function("slice", &slice_js);
    emscripten::function("isosurface", &isosurface_js);
    emscripten::function("cropBbox", &crop_bbox_js);
    emscripten::function("cropPlane", &crop_plane_js);
    emscripten::function("split", &split_js);
    emscripten::function("convertCells", &convert_cells_js);
    emscripten::function("refine", &refine_js);
    emscripten::function("decimate", &decimate_js);
    emscripten::function("partition", &partition_js);
    emscripten::function("partitionLabels", &partition_labels_js);
    emscripten::function("stats", &stats_js);
    emscripten::function("dataDrop", &data_drop_js);
    emscripten::function("dataKeep", &data_keep_js);
    emscripten::function("dataRename", &data_rename_js);
    emscripten::function("dataPointToCell", &data_point_to_cell_js);
    emscripten::function("dataCellToPoint", &data_cell_to_point_js);
    emscripten::function("dataCalc", &data_calc_js);
    emscripten::function("dataCondition", &data_condition_js);
    emscripten::function("dataInfo", &data_info_js);
}
