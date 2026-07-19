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
 * Format scope (v1): the 28 formats with no HDF5/netCDF dependency, plus
 * XDMF's XML/Binary data path (not its HDF variant) -- 29 readable, 28
 * writable (`openfoam` is read-only). CGNS/H5M/HMF/MED/Exodus
 * are not registered here -- porting HDF5/netCDF to WASM is a separate,
 * larger undertaking (see doc/wasm.md). Ambiguous extensions (`.msh` shared
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
#include "meshioplusplus/operations/diff.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/quality.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/operations/surface.hpp"
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
    return mesh;
}

// ---------------------------------------------------------------------
// Format dispatch goes through the shared registry (registry.hpp), the
// JS-side analogue of Python's extension_to_filetypes in _helpers.py --
// bindings/_core.cpp exposes one function per format and leaves dispatch
// entirely to Python, while the flat bindings (this file and bindings_c/)
// share the C++-level tables in cpp/src/registry.cpp. Parameterized writers
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
val read_mesh(const std::string& rPath, const std::string& rFormat) {
    return with_js_errors([&]() -> val {
        std::string fmt;
        try {
            fmt = resolve_format(rPath, rFormat);
        } catch (const meshioplusplus::ReadError&) {
            fmt = meshioplusplus::sniff_format(rPath);
            if (fmt.empty())
                throw;
        }
        auto it = registry_readers().find(fmt);
        if (it == registry_readers().end())
            throw meshioplusplus::ReadError("meshio++ (wasm): unknown or unsupported format '" +
                                            fmt + "'" + compiled_out_hint(fmt));
        return mesh_to_val(it->second(rPath));
    });
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

}  // namespace

EMSCRIPTEN_BINDINGS(meshioplusplus_wasm) {
    emscripten::function("readMesh", &read_mesh);
    emscripten::function("writeMesh", &write_mesh);
    emscripten::function("convert", &convert);
    emscripten::function("numNodesPerCell", &num_nodes_per_cell_js);
    emscripten::function("topologicalDimension", &topological_dimension_js);
    emscripten::function("meshBackend", &mesh_backend_js);
    emscripten::function("extractSurface", &extract_surface_js);
    emscripten::function("extractSkin", &extract_skin_js);
    emscripten::function("attachQuality", &attach_quality_js);
    emscripten::function("sniffFormat", &sniff_format_js);
    emscripten::function("reorder", &reorder_js);
    emscripten::function("computeBandwidth", &compute_bandwidth_js);
    emscripten::function("diff", &diff_js);
    emscripten::function("meshesEqual", &meshes_equal_js);
    emscripten::function("merge", &merge_js);
}
