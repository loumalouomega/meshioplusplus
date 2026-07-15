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
 * @file np_conversions.hpp
 * @brief The pybind11 <-> C++ `meshioplusplus::Mesh` conversion boundary.
 *
 * This header is the single choke point every C++-backed format binding goes
 * through to turn a Python `meshioplusplus.Mesh` into a C++
 * `meshioplusplus::Mesh` (for writing) and back (for reading). The design
 * goal is to make the crossing as cheap as possible:
 *
 * - **Python -> C++ (`py_to_mesh`)** is view-based: for every rectangular
 *   numpy array reachable from the mesh (points, per-block cell
 *   connectivity, point_data/cell_data/field_data arrays) it builds a
 *   `meshioplusplus::NDArray` that *points into* the existing numpy buffer
 *   rather than copying it. The numpy `py::array` objects themselves are
 *   kept alive in a `PyMeshRefs` "keepalive" vector for exactly as long as
 *   the C++ views need to remain valid (see `PyMeshRefs` below).
 * - **C++ -> Python (`mesh_to_py`)** is also zero-copy, but in the opposite
 *   direction: each `NDArray` produced by a C++ reader is "adopted" by a
 *   numpy array via `numpy_from_ndarray`, which transfers ownership of the
 *   heap buffer to a `py::capsule` so numpy (not the C++ side) is the last
 *   one to free it.
 * - **Ragged cell blocks** (jagged polygon / polyhedron connectivity, which
 *   cannot be represented as a rectangular `NDArray`) are the one case that
 *   is *not* zero-copy in either direction: they are always materialized as
 *   Python lists of numpy arrays (`ragged_data_to_py`) or parsed from them
 *   (`ragged_cellblock_from_py`), because there is no rectangular buffer to
 *   view into.
 *
 * See the "C++ core" section of the repository's top-level `CLAUDE.md` for
 * the broader architectural picture (side-channel structs such as `MedInfo`
 * for data this layer intentionally does not carry, the `allow_ragged`
 * opt-in policy, etc.).
 */

// System includes
#include <string>
#include <vector>

// External includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"

namespace py = pybind11;

namespace meshioplusplus_py {

/**
 * @brief Map a numpy dtype to the corresponding `meshioplusplus::DType`.
 *
 * Supports the floating-point (`f4`/`f8`), signed integer (`i1`/`i2`/`i4`/
 * `i8`) and unsigned integer (`u1`/`u2`/`u4`/`u8`) kinds/itemsizes that
 * `meshioplusplus::DType` models; anything else (e.g. complex, object,
 * bool, string dtypes) is unsupported by the C++ core.
 *
 * @param dt A numpy dtype, as obtained from `py::array::dtype()`.
 * @return The matching `meshioplusplus::DType` enumerator.
 * @throws meshioplusplus::WriteError if the dtype's (kind, itemsize) pair
 *         has no C++-side representation.
 */
inline meshioplusplus::DType dtype_from_numpy(const py::dtype& dt) {
    const char kind = dt.kind();
    const py::ssize_t isz = dt.itemsize();
    using meshioplusplus::DType;
    if (kind == 'f') {
        if (isz == 4) return DType::Float32;
        if (isz == 8) return DType::Float64;
    } else if (kind == 'i') {
        if (isz == 1) return DType::Int8;
        if (isz == 2) return DType::Int16;
        if (isz == 4) return DType::Int32;
        if (isz == 8) return DType::Int64;
    } else if (kind == 'u') {
        if (isz == 1) return DType::UInt8;
        if (isz == 2) return DType::UInt16;
        if (isz == 4) return DType::UInt32;
        if (isz == 8) return DType::UInt64;
    }
    throw meshioplusplus::WriteError(std::string("Unsupported numpy dtype '") + dt.kind() +
                             std::to_string(isz) + "' for the meshio++ C++ core");
}

/**
 * @brief Keepalive list for the numpy arrays a `meshioplusplus::Mesh` views into.
 *
 * `py_to_mesh` builds `meshioplusplus::NDArray` values that are *views*
 * (non-owning pointers) into existing numpy buffers rather than copies. A
 * view is only valid for as long as the underlying `py::array` it points
 * into is alive; if that `py::array` were a temporary (e.g. the result of
 * `ensure_contiguous`'s byte-order/contiguity coercion, or of pybind11
 * borrowing/casting), it could be destroyed - and its buffer freed - the
 * moment the local C++ variable holding it goes out of scope, silently
 * dangling the view.
 *
 * `PyMeshRefs::keep` exists to prevent exactly that: every `py::array` that
 * ends up backing a view is additionally pushed onto `keep`, which the
 * caller (typically a format-binding function in `bindings/_core.cpp`)
 * keeps alive on its stack for the full duration it uses the resulting
 * `meshioplusplus::Mesh` (i.e. for the whole write call). Once that scope
 * ends, the mesh's views must not be dereferenced anymore.
 *
 * @note This is a value type with a single public member (`keep`) rather
 *       than an opaque handle so callers can simply declare one on the
 *       stack alongside the `Mesh` they build with `py_to_mesh`.
 */
struct PyMeshRefs {
    std::vector<py::array> keep;
};

/**
 * @brief Build a non-owning `meshioplusplus::NDArray` view over a numpy array's buffer.
 *
 * Reads the numpy array's dtype and shape and wraps its raw data pointer in
 * `meshioplusplus::NDArray::make_view`. This performs **no copy**: the
 * returned `NDArray` aliases `a`'s memory directly, so `a` (or whichever
 * `py::array` owns the same buffer) must outlive every use of the view -
 * see `PyMeshRefs`, which is how callers of this function keep that
 * guarantee.
 *
 * @param a A C-contiguous numpy array (callers are expected to have already
 *          run it through `ensure_contiguous`, which also normalizes byte
 *          order, so the raw bytes can be reinterpreted directly by dtype).
 * @return An `NDArray` view (`is_view() == true`) sharing `a`'s buffer.
 */
inline meshioplusplus::NDArray view_from_numpy(const py::array& a) {
    meshioplusplus::DType dt = dtype_from_numpy(a.dtype());
    std::vector<std::size_t> shape(static_cast<std::size_t>(a.ndim()));
    for (py::ssize_t i = 0; i < a.ndim(); ++i)
        shape[static_cast<std::size_t>(i)] = static_cast<std::size_t>(a.shape(i));
    auto* ptr = reinterpret_cast<std::byte*>(const_cast<void*>(a.data()));
    return meshioplusplus::NDArray::make_view(dt, std::move(shape), ptr);
}

/**
 * @brief Coerce a Python object to a C-contiguous, native-byte-order numpy array
 *        and register it in the keepalive list.
 *
 * Two normalizations happen here so that later code (`view_from_numpy`,
 * typed C++ reads) can treat the buffer as a plain native-endian C array:
 *  1. `py::array::ensure(..., c_style | forcecast)` forces C-contiguity,
 *     converting/copying if `obj` is not already an array or is
 *     Fortran-ordered/non-contiguous.
 *  2. **Byte-order normalization**: numpy's `dtype.byteorder` is one of
 *     `'='` (native), `'|'` (not applicable, e.g. single-byte types), `'<'`
 *     (little-endian) or `'>'` (big-endian). The C++ core assumes a
 *     little-endian host (x86/ARM64), so `'='`, `'|'` and `'<'` are all
 *     accepted as-is, but a big-endian (`'>'`) array - which can arise from
 *     e.g. reading a big-endian binary file format into numpy - is byte-
 *     swapped via `dtype.newbyteorder("=")` + `.astype(...)` before use, so
 *     the typed views the C++ side takes read correctly.
 *
 * Either normalization step may allocate a new array; in that case (and in
 * the common case where none was needed) the resulting `py::array` is
 * pushed onto `refs.keep` so it outlives any `NDArray` view built from it.
 *
 * @param obj An array-like Python object (numpy array or anything
 *            `py::array::ensure` can convert).
 * @param refs Keepalive list the resulting array is appended to.
 * @return A C-contiguous, native-byte-order `py::array`.
 * @throws meshioplusplus::WriteError if `obj` cannot be interpreted as an array.
 */
inline py::array ensure_contiguous(py::handle obj, PyMeshRefs& refs) {
    py::array a = py::array::ensure(obj, py::array::c_style | py::array::forcecast);
    if (!a) throw meshioplusplus::WriteError("Expected an array-like object");
    // Normalize to native byte order so the typed views read correctly. numpy
    // dtype.byteorder is '=' native, '|' n/a, '<' little, '>' big. Host is
    // assumed little-endian (x86/ARM64).
    std::string bo = py::cast<std::string>(a.dtype().attr("byteorder"));
    const bool native = (bo == "=" || bo == "|" || bo == "<");
    if (!native) {
        py::object newdt = a.dtype().attr("newbyteorder")("=");
        a = py::array::ensure(a.attr("astype")(newdt), py::array::c_style);
    }
    refs.keep.push_back(a);
    return a;
}

/**
 * @brief Parse a ragged Python cell-block `data` (a list, not an ndarray)
 *        into a C++ `CellBlock`'s ragged members.
 *
 * meshio++ represents cell blocks whose cells don't share a fixed node
 * count as plain Python lists rather than rectangular numpy arrays, so
 * there is no buffer to view into - this function always **copies** the
 * node ids into `std::vector`s owned by the returned `CellBlock`.
 *
 * The nesting depth depends on the cell type name:
 *  - A `"polyhedron"`-prefixed block is 2-level: a list of cells, each a
 *    list of faces, each a sequence of node ids -> populates
 *    `CellBlock::polyhedron_rows`.
 *  - Any other ragged block (a `"polygon"` block with varying node counts
 *    per cell) is 1-level: a list of cells, each a sequence of node ids ->
 *    populates `CellBlock::polygon_rows`.
 *
 * @param type The meshio++ cell type name (e.g. `"polygon"`,
 *             `"polyhedron4"`); consumed by move into the returned block.
 * @param data_obj The Python `cells[i].data` list for this block.
 * @return A `CellBlock` with `type` set and exactly one of
 *         `polygon_rows`/`polyhedron_rows` populated (owned copies).
 */
inline meshioplusplus::CellBlock ragged_cellblock_from_py(std::string type,
                                                  py::handle data_obj) {
    meshioplusplus::CellBlock cb;
    cb.type = std::move(type);
    auto to_ids = [](py::handle seq) {
        std::vector<std::int64_t> ids;
        for (py::handle v : seq) ids.push_back(py::cast<std::int64_t>(v));
        return ids;
    };
    if (cb.type.rfind("polyhedron", 0) == 0) {
        for (py::handle cell : data_obj) {
            std::vector<std::vector<std::int64_t>> faces;
            for (py::handle face : cell) faces.push_back(to_ids(face));
            cb.polyhedron_rows.push_back(std::move(faces));
        }
    } else {
        for (py::handle row : data_obj) cb.polygon_rows.push_back(to_ids(row));
    }
    return cb;
}

/**
 * @brief Convert a Python `meshioplusplus.Mesh` into a C++ `meshioplusplus::Mesh`,
 *        for use by a format writer's C++ binding.
 *
 * This is the write-side half of the conversion boundary: every rectangular
 * numpy array reachable from `pymesh` (points, each cell block's
 * connectivity, and every point_data/cell_data/field_data array) is turned
 * into a **non-owning view** (`view_from_numpy`) over the same memory numpy
 * already holds - no bytes are copied for the rectangular case. Each source
 * `py::array` is first passed through `ensure_contiguous` (which may itself
 * allocate a fresh, C-contiguous/native-byte-order array when the input
 * isn't already one) and is kept alive via `refs` for as long as the
 * returned `Mesh`'s views are used; see `PyMeshRefs`.
 *
 * @param pymesh A Python `meshioplusplus.Mesh` instance (as a `py::handle`;
 *               its `points`, `cells`, `point_data`, `cell_data` and
 *               `field_data` attributes are read).
 * @param refs Keepalive list; every numpy array backing a view built here
 *             is appended to it. Must outlive the returned `Mesh`.
 * @param lenient_field_data When `false` (default), every `field_data`
 *             entry is coerced to a numeric array via `ensure_contiguous`,
 *             and a non-numeric entry throws. When `true`, a `field_data`
 *             entry that fails that coercion (e.g. MED's `"med:nom"`,
 *             which stores a list of strings rather than numbers) is
 *             silently skipped instead of raising - the caller is
 *             responsible for carrying that entry through its own
 *             format-specific side-channel (e.g. `MedInfo`) instead of
 *             through this generic conversion path.
 * @param allow_ragged When `false` (the default), encountering a ragged
 *             cell block - one whose Python `data` is a `list` rather than
 *             an ndarray, i.e. jagged polygon or polyhedron connectivity -
 *             throws `meshioplusplus::WriteError`. This is deliberate
 *             regression-safety: most C++ format writers only handle
 *             rectangular `NDArray` blocks, and by default rejecting ragged
 *             input here means such a writer safely raises and the
 *             caller's Python fallback takes over, rather than silently
 *             mishandling or truncating the mesh. Only the ragged-aware
 *             bindings (currently MED write) pass `allow_ragged=true` to
 *             opt into parsing ragged blocks via `ragged_cellblock_from_py`
 *             (which always copies, since ragged data has no rectangular
 *             buffer to view into).
 * @return A `meshioplusplus::Mesh` whose rectangular array members alias
 *         Python-owned memory (valid only while `refs` is alive) and whose
 *         ragged cell blocks (if any, and if `allow_ragged`) own independent
 *         copies.
 * @throws meshioplusplus::WriteError if a ragged block is encountered while
 *         `allow_ragged` is `false`, or if any array cannot be coerced to a
 *         supported dtype/shape.
 */
inline meshioplusplus::Mesh py_to_mesh(py::handle pymesh, PyMeshRefs& refs,
                               bool lenient_field_data = false,
                               bool allow_ragged = false) {
    meshioplusplus::Mesh m;

    m.points = view_from_numpy(ensure_contiguous(pymesh.attr("points"), refs));

    for (py::handle cb : pymesh.attr("cells")) {
        std::string type = py::cast<std::string>(cb.attr("type"));
        py::object data_obj = py::reinterpret_borrow<py::object>(cb.attr("data"));
        // Ragged polyhedron / jagged polygon data is a Python list, not an
        // ndarray.
        if (py::isinstance<py::list>(data_obj)) {
            if (!allow_ragged) {
                throw meshioplusplus::WriteError(
                    "ragged (polyhedron / jagged polygon) cell blocks are not "
                    "handled by this C++ format");
            }
            m.cells.push_back(ragged_cellblock_from_py(std::move(type), data_obj));
            continue;
        }
        py::array d = ensure_contiguous(data_obj, refs);
        m.cells.emplace_back(std::move(type), view_from_numpy(d));
    }

    for (auto item : pymesh.attr("point_data").cast<py::dict>()) {
        std::string name = py::cast<std::string>(item.first);
        py::array d = ensure_contiguous(item.second, refs);
        m.point_data.emplace(std::move(name), view_from_numpy(d));
    }

    for (auto item : pymesh.attr("cell_data").cast<py::dict>()) {
        std::string name = py::cast<std::string>(item.first);
        std::vector<meshioplusplus::NDArray> blocks;
        for (py::handle a : py::reinterpret_borrow<py::object>(item.second)) {
            py::array d = ensure_contiguous(a, refs);
            blocks.push_back(view_from_numpy(d));
        }
        m.cell_data.emplace(std::move(name), std::move(blocks));
    }

    py::object fd = pymesh.attr("field_data");
    if (!fd.is_none()) {
        for (auto item : fd.cast<py::dict>()) {
            std::string name = py::cast<std::string>(item.first);
            if (lenient_field_data) {
                try {
                    py::array d = ensure_contiguous(item.second, refs);
                    m.field_data.emplace(std::move(name), view_from_numpy(d));
                } catch (...) {
                    // non-numeric entry: handled by the caller's side-channel
                }
            } else {
                py::array d = ensure_contiguous(item.second, refs);
                m.field_data.emplace(std::move(name), view_from_numpy(d));
            }
        }
    }

    return m;
}

/**
 * @brief Adopt an `NDArray`'s buffer into a capsule-backed, writeable numpy array.
 *
 * This is the read-side counterpart of `view_from_numpy`/`ensure_contiguous`:
 * instead of copying `arr`'s data into a fresh numpy array, ownership of the
 * buffer is transferred to Python.
 *
 * Mechanics: `arr` is moved onto the heap, `make_owned()` is called on it so
 * it holds (or already holds, if it wasn't a view) an independently
 * allocated buffer it is responsible for freeing, and that heap-allocated
 * `NDArray*` is wrapped in a `py::capsule` whose destructor `delete`s it.
 * The returned `py::array` is constructed directly over the `NDArray`'s
 * data pointer with the capsule as its owner, so:
 *  - No element is copied by this function itself (`make_owned()` only
 *    copies if `arr` was still a view over someone else's memory - see
 *    `NDArray::make_owned`).
 *  - The numpy array is writeable and remains valid for exactly as long as
 *    Python holds a reference to it; once the last reference drops, the
 *    capsule destructor deletes the heap `NDArray`, freeing the buffer.
 *  - C-contiguous strides are computed from `arr`'s shape (row-major,
 *    itemsize-scaled), treating a zero-length dimension as stride-compatible
 *    with size 1 to avoid a zero multiplier.
 *
 * @param arr An rvalue `NDArray`, consumed by move; the caller must not use
 *            it afterward.
 * @return A new numpy array whose memory is owned (via capsule) by a
 *         heap-allocated copy of `arr`.
 */
inline py::array numpy_from_ndarray(meshioplusplus::NDArray&& arr) {
    auto* heap = new meshioplusplus::NDArray(std::move(arr));
    heap->make_owned();
    py::capsule owner(heap, [](void* p) {
        delete reinterpret_cast<meshioplusplus::NDArray*>(p);
    });

    std::vector<py::ssize_t> shape(heap->shape().begin(), heap->shape().end());
    std::vector<py::ssize_t> strides(shape.size());
    const py::ssize_t itemsize = static_cast<py::ssize_t>(meshioplusplus::dtype_size(heap->dtype()));
    py::ssize_t s = itemsize;
    for (int i = static_cast<int>(shape.size()) - 1; i >= 0; --i) {
        strides[static_cast<std::size_t>(i)] = s;
        s *= (shape[static_cast<std::size_t>(i)] == 0 ? 1 : shape[static_cast<std::size_t>(i)]);
    }
    return py::array(py::dtype(meshioplusplus::dtype_numpy_str(heap->dtype())), shape, strides,
                     heap->data(), owner);
}

/**
 * @brief Build the Python object a ragged `CellBlock` maps to.
 *
 * This is the read-side counterpart of `ragged_cellblock_from_py`. Since a
 * ragged block has no rectangular buffer, its rows/faces are **copied**
 * (via `std::memcpy` into freshly allocated `py::array_t<std::int64_t>`
 * objects) rather than adopted zero-copy the way `numpy_from_ndarray` does
 * for rectangular blocks:
 *  - For a jagged polygon block (`cb.polygon_rows` non-empty): a Python
 *    list of 1-D int64 numpy arrays, one per cell.
 *  - For a polyhedron block (`cb.polyhedron_rows` non-empty): a Python list
 *    of cells, each itself a list of 1-D int64 numpy arrays, one per face.
 *
 * This matches exactly what `meshioplusplus.Mesh`/`CellBlock` expect to
 * store for these cell types on the Python side (kept as a Python list,
 * never coerced into a rectangular ndarray).
 *
 * @param cb A `CellBlock` for which `cb.is_ragged()` is `true`.
 * @return A `py::object` (a `py::list`) as described above.
 */
inline py::object ragged_data_to_py(const meshioplusplus::CellBlock& cb) {
    auto ids_to_arr = [](const std::vector<std::int64_t>& ids) {
        py::array_t<std::int64_t> a(static_cast<py::ssize_t>(ids.size()));
        if (!ids.empty())
            std::memcpy(a.mutable_data(), ids.data(), ids.size() * sizeof(std::int64_t));
        return a;
    };
    py::list out;
    if (!cb.polyhedron_rows.empty()) {
        for (const auto& cell : cb.polyhedron_rows) {
            py::list faces;
            for (const auto& face : cell) faces.append(ids_to_arr(face));
            out.append(faces);
        }
    } else {
        for (const auto& row : cb.polygon_rows) out.append(ids_to_arr(row));
    }
    return out;
}

/**
 * @brief Convert a C++ `meshioplusplus::Mesh` into a Python `meshioplusplus.Mesh`,
 *        for use by a format reader's C++ binding.
 *
 * This is the read-side counterpart of `py_to_mesh`. `m` is consumed by
 * move (`meshioplusplus::Mesh&&`) and every rectangular array member
 * (`points`, each rectangular cell block's `data`, and every
 * point_data/cell_data/field_data array) is handed to Python via
 * `numpy_from_ndarray`, which transfers buffer ownership to numpy through a
 * capsule rather than copying. Ragged cell blocks (`cb.is_ragged()`) are
 * the one exception: they have no rectangular buffer to adopt, so they are
 * copied into Python lists via `ragged_data_to_py`.
 *
 * The resulting `meshioplusplus.Mesh` is constructed by importing the
 * `meshioplusplus` Python module and calling its `Mesh` class directly, so
 * this function has no compile-time dependency on the Python-side `Mesh`
 * definition.
 *
 * @param m The C++ mesh to convert, consumed by move; the caller must not
 *          use it afterward (its `NDArray` members are moved out one by
 *          one as they're adopted).
 * @return A new `py::object` wrapping a `meshioplusplus.Mesh` instance whose
 *         rectangular arrays are zero-copy views owned via capsules, and
 *         whose ragged cell block data (if any) are freshly copied Python
 *         lists.
 * @note This function does not carry over `mesh.info`, `cell_sets` or
 *       `point_sets` - per the conversion-boundary design, formats that
 *       need those attach them out-of-band via a side-channel struct
 *       (e.g. `MedInfo`, `AnsysInfo`, `OpenFoamInfo`) that the binding
 *       `setattr`s onto the returned Python object separately.
 */
inline py::object mesh_to_py(meshioplusplus::Mesh&& m) {
    py::object MeshCls = py::module_::import("meshioplusplus").attr("Mesh");

    py::array points = numpy_from_ndarray(std::move(m.points));

    py::list cells;
    for (auto& cb : m.cells) {
        if (cb.is_ragged()) {
            cells.append(py::make_tuple(py::str(cb.type), ragged_data_to_py(cb)));
        } else {
            cells.append(
                py::make_tuple(py::str(cb.type), numpy_from_ndarray(std::move(cb.data))));
        }
    }

    py::dict point_data;
    for (auto& kv : m.point_data)
        point_data[py::str(kv.first)] = numpy_from_ndarray(std::move(const_cast<meshioplusplus::NDArray&>(kv.second)));

    py::dict cell_data;
    for (auto& kv : m.cell_data) {
        py::list lst;
        for (auto& a : kv.second) lst.append(numpy_from_ndarray(std::move(a)));
        cell_data[py::str(kv.first)] = lst;
    }

    py::dict field_data;
    for (auto& kv : m.field_data)
        field_data[py::str(kv.first)] =
            numpy_from_ndarray(std::move(const_cast<meshioplusplus::NDArray&>(kv.second)));

    return MeshCls(points, cells, py::arg("point_data") = point_data,
                   py::arg("cell_data") = cell_data,
                   py::arg("field_data") = field_data);
}

}  // namespace meshioplusplus_py
