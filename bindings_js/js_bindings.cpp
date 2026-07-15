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
 * Format scope (v1): the 26 formats with no HDF5/netCDF dependency, plus
 * XDMF's XML/Binary data path (not its HDF variant). CGNS/H5M/HMF/MED/Exodus
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
#include <functional>
#include <map>
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
#include "meshioplusplus/formats/abaqus.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/ansysinp.hpp"
#include "meshioplusplus/formats/avsucd.hpp"
#include "meshioplusplus/formats/dolfin.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/flux.hpp"
#include "meshioplusplus/formats/freefem.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/medit.hpp"
#include "meshioplusplus/formats/mfm.hpp"
#include "meshioplusplus/formats/mphtxt.hpp"
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/formats/netgen.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/formats/openfoam.hpp"
#include "meshioplusplus/formats/permas.hpp"
#include "meshioplusplus/formats/ply.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/su2.hpp"
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/formats/tetgen.hpp"
#include "meshioplusplus/formats/ugrid.hpp"
#include "meshioplusplus/formats/unv.hpp"
#include "meshioplusplus/formats/vtk.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/formats/wkt.hpp"
#include "meshioplusplus/formats/xdmf.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/types.hpp"

using emscripten::val;

namespace {

using meshioplusplus::CellBlock;
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
    out.set("points", ndarray_to_float64_array(rMesh.mPoints));
    out.set("dim", static_cast<int>(cols_of(rMesh.mPoints)));

    val cells = val::array();
    for (const auto& cb : rMesh.mCells) {
        if (cb.IsRagged())
            throw meshioplusplus::ReadError("meshio++ (wasm): ragged cell blocks ('" + cb.mType +
                                            "') are not supported by the JS API yet");
        val block = val::object();
        block.set("type", cb.mType);
        block.set("data", ndarray_to_int32_array(cb.mData));
        block.set("nodesPerCell", static_cast<int>(cols_of(cb.mData)));
        cells.call<void>("push", block);
    }
    out.set("cells", cells);

    val point_data = val::object();
    for (const auto& kv : rMesh.mPointData)
        point_data.set(kv.first, ndarray_to_float64_array(kv.second));
    out.set("point_data", point_data);

    val cell_data = val::object();
    for (const auto& kv : rMesh.mCellData) {
        val blocks = val::array();
        for (const auto& arr : kv.second)
            blocks.call<void>("push", ndarray_to_float64_array(arr));
        cell_data.set(kv.first, blocks);
    }
    out.set("cell_data", cell_data);

    val field_data = val::object();
    for (const auto& kv : rMesh.mFieldData)
        field_data.set(kv.first, ndarray_to_float64_array(kv.second));
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
    mesh.mPoints = float64_ndarray_from_val(points_val, {npts_len / dim, dim});

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
        mesh.mCells.emplace_back(
            type, int64_ndarray_from_val(data_val, {data_len / nodes_per_cell, nodes_per_cell}));
    }

    if (rObj.hasOwnProperty("point_data")) {
        val pd = rObj["point_data"];
        for (const std::string& name : js_object_keys(pd)) {
            val arr = pd[name];
            mesh.mPointData.emplace(
                name, float64_ndarray_from_val(arr, {arr["length"].as<std::size_t>()}));
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
            mesh.mCellData.emplace(name, std::move(blocks));
        }
    }
    if (rObj.hasOwnProperty("field_data")) {
        val fd = rObj["field_data"];
        for (const std::string& name : js_object_keys(fd)) {
            val arr = fd[name];
            mesh.mFieldData.emplace(
                name, float64_ndarray_from_val(arr, {arr["length"].as<std::size_t>()}));
        }
    }
    return mesh;
}

// ---------------------------------------------------------------------
// Format dispatch tables (the JS-side analogue of Python's
// extension_to_filetypes in _helpers.py, which has no C++-level
// equivalent -- bindings/_core.cpp exposes one function per format and
// leaves dispatch entirely to Python). Parameterized writers get a fixed
// default here (documented per-entry), matching each format's own Python
// reference default; per-call overrides are a possible future API addition,
// deliberately out of scope for v1.
// ---------------------------------------------------------------------

using ReadFn = std::function<Mesh(const std::string&)>;
using WriteFn = std::function<void(const std::string&, const Mesh&)>;

const std::map<std::string, ReadFn>& readers() {
    static const std::map<std::string, ReadFn> m = {
        {"abaqus", meshioplusplus::read_abaqus},
        {"ansys", meshioplusplus::read_ansys},
        {"avsucd", meshioplusplus::read_avsucd},
        {"dolfin", meshioplusplus::read_dolfin},
        {"flac3d", meshioplusplus::read_flac3d},
        {"flux", meshioplusplus::read_flux},
        {"freefem", meshioplusplus::read_freefem},
        {"gmsh", meshioplusplus::read_gmsh},
        {"medit", meshioplusplus::read_medit_ascii},
        {"mfm", meshioplusplus::read_mfm},
        {"mphtxt", meshioplusplus::read_mphtxt},
        {"nastran", meshioplusplus::read_nastran},
        {"netgen", meshioplusplus::read_netgen},
        {"obj", meshioplusplus::read_obj},
        {"off", meshioplusplus::read_off},
        {"permas", meshioplusplus::read_permas},
        {"ply", meshioplusplus::read_ply},
        {"stl", meshioplusplus::read_stl},
        {"su2", meshioplusplus::read_su2},
        {"tecplot", meshioplusplus::read_tecplot},
        {"tetgen", meshioplusplus::read_tetgen},
        {"ugrid", meshioplusplus::read_ugrid},
        {"unv", meshioplusplus::read_unv},
        {"vtk", meshioplusplus::read_vtk},
        {"vtu", meshioplusplus::read_vtu},
        {"wkt", meshioplusplus::read_wkt},
        {"xdmf", meshioplusplus::read_xdmf},
        // Side-channel info (point_sets/cell_sets, cell-tag family names) is
        // not yet exposed to JS -- v1 limitation, see doc/wasm.md.
        {"ansysinp",
         [](const std::string& path) {
             meshioplusplus::AnsysInfo info;
             return meshioplusplus::read_ansysinp(path, info);
         }},
        {"openfoam",
         [](const std::string& path) {
             meshioplusplus::OpenFoamInfo info;
             return meshioplusplus::read_openfoam(path, info);
         }},
    };
    return m;
}

const std::map<std::string, WriteFn>& writers() {
    static const std::map<std::string, WriteFn> m = {
        {"abaqus", meshioplusplus::write_abaqus},
        {"ansys", [](const std::string& p,
                     const Mesh& mm) { meshioplusplus::write_ansys(p, mm, /*binary=*/true); }},
        {"avsucd", meshioplusplus::write_avsucd},
        {"dolfin", meshioplusplus::write_dolfin},
        {"flac3d",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_flac3d(p, mm, ".16e", /*binary=*/false);
         }},
        {"flux", meshioplusplus::write_flux},
        {"freefem", meshioplusplus::write_freefem},
        {"gmsh", [](const std::string& p,
                    const Mesh& mm) { meshioplusplus::write_gmsh41(p, mm, /*binary=*/true); }},
        {"medit", meshioplusplus::write_medit_ascii},
        {"mfm",
         [](const std::string& p, const Mesh& mm) { meshioplusplus::write_mfm(p, mm, ".16e"); }},
        {"mphtxt", meshioplusplus::write_mphtxt},
        {"nastran", meshioplusplus::write_nastran},
        {"netgen",
         [](const std::string& p, const Mesh& mm) { meshioplusplus::write_netgen(p, mm, ".16e"); }},
        {"obj", meshioplusplus::write_obj},
        {"off", meshioplusplus::write_off},
        {"permas", meshioplusplus::write_permas},
        {"ply", [](const std::string& p,
                   const Mesh& mm) { meshioplusplus::write_ply(p, mm, /*binary=*/true); }},
        {"stl", [](const std::string& p,
                   const Mesh& mm) { meshioplusplus::write_stl(p, mm, /*binary=*/false); }},
        {"su2", meshioplusplus::write_su2},
        {"tecplot", meshioplusplus::write_tecplot},
        {"tetgen", meshioplusplus::write_tetgen},
        {"ugrid", meshioplusplus::write_ugrid},
        {"unv", meshioplusplus::write_unv},
        {"vtk",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_vtk(p, mm, /*binary=*/true, /*v51=*/true);
         }},
        {"vtu",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::write_vtu(p, mm, /*binary=*/true, /*zlib=*/true);
         }},
        {"wkt", meshioplusplus::write_wkt},
        {"xdmf",
         [](const std::string& p, const Mesh& mm) { meshioplusplus::write_xdmf(p, mm, "XML"); }},
        {"ansysinp",
         [](const std::string& p, const Mesh& mm) {
             meshioplusplus::AnsysInfo info;  // no point_sets/cell_sets from JS in v1
             meshioplusplus::write_ansysinp(p, mm, info);
         }},
        // openfoam is read-only in the C++ core (see openfoam.hpp) -> no writer entry.
    };
    return m;
}

// Extension -> canonical format key for the non-ambiguous cases; `.msh`
// defaults to gmsh and `.inp` to abaqus (matching this repo's own import
// order in src/meshioplusplus/__init__.py). Pass an explicit `format` to
// select ansys/freefem (.msh) or ansysinp (.inp) instead.
const std::map<std::string, std::string>& extension_defaults() {
    static const std::map<std::string, std::string> m = {
        {".inp", "abaqus"},  {".avs", "avsucd"},  {".xml", "dolfin"},  {".f3grid", "flac3d"},
        {".pf3", "flux"},    {".mesh", "medit"},  {".mfm", "mfm"},     {".mphtxt", "mphtxt"},
        {".bdf", "nastran"}, {".nas", "nastran"}, {".fem", "nastran"}, {".vol", "netgen"},
        {".obj", "obj"},     {".off", "off"},     {".post", "permas"}, {".dato", "permas"},
        {".ply", "ply"},     {".stl", "stl"},     {".su2", "su2"},     {".dat", "tecplot"},
        {".tec", "tecplot"}, {".ele", "tetgen"},  {".node", "tetgen"}, {".ugrid", "ugrid"},
        {".unv", "unv"},     {".vtk", "vtk"},     {".vtu", "vtu"},     {".wkt", "wkt"},
        {".xdmf", "xdmf"},   {".xmf", "xdmf"},    {".msh", "gmsh"},
    };
    return m;
}

std::string extension_of(const std::string& rPath) {
    auto pos = rPath.find_last_of('.');
    return pos == std::string::npos ? "" : rPath.substr(pos);
}

std::string resolve_format(const std::string& rPath, const std::string& rFormat) {
    if (!rFormat.empty())
        return rFormat;
    auto it = extension_defaults().find(extension_of(rPath));
    if (it == extension_defaults().end())
        throw meshioplusplus::ReadError("meshio++ (wasm): cannot infer format from '" + rPath +
                                        "' -- pass an explicit format argument");
    return it->second;
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
 * @param rFormat explicit format key (see `extension_defaults()`), or "" to
 *   infer from `rPath`'s extension.
 * @return a plain JS mesh object (see `mesh_to_val`).
 * @throws meshioplusplus::ReadError on an unknown/unsupported format or a
 *   malformed file.
 */
val read_mesh(const std::string& rPath, const std::string& rFormat) {
    return with_js_errors([&]() -> val {
        std::string fmt = resolve_format(rPath, rFormat);
        auto it = readers().find(fmt);
        if (it == readers().end())
            throw meshioplusplus::ReadError("meshio++ (wasm): unknown or unsupported format '" +
                                            fmt + "'");
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
        auto it = writers().find(fmt);
        if (it == writers().end())
            throw meshioplusplus::WriteError(
                "meshio++ (wasm): unknown, read-only, or unsupported format '" + fmt + "'");
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
        auto rit = readers().find(rfmt);
        auto wit = writers().find(wfmt);
        if (rit == readers().end())
            throw meshioplusplus::ReadError(
                "meshio++ (wasm): unknown or unsupported input format '" + rfmt + "'");
        if (wit == writers().end())
            throw meshioplusplus::WriteError(
                "meshio++ (wasm): unknown, read-only, or unsupported output format '" + wfmt + "'");
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

}  // namespace

EMSCRIPTEN_BINDINGS(meshioplusplus_wasm) {
    emscripten::function("readMesh", &read_mesh);
    emscripten::function("writeMesh", &write_mesh);
    emscripten::function("convert", &convert);
    emscripten::function("numNodesPerCell", &num_nodes_per_cell_js);
    emscripten::function("topologicalDimension", &topological_dimension_js);
}
