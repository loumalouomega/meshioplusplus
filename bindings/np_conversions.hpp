#pragma once
//
// Conversions between meshioplusplus::Mesh (C++) and the pure-Python
// meshioplusplus.Mesh,
// implementing the "zero-copy at the I/O boundary" strategy:
//
//   * py_to_mesh  : Python mesh -> C++ mesh whose NDArrays are non-owning
//                   *views* over the numpy buffers (kept alive via PyMeshRefs).
//                   Used on the write path -> no input copy.
//   * mesh_to_py  : C++ mesh -> Python mesh; each NDArray's buffer is moved
//                   into a capsule that backs a writeable numpy array.
//                   Used on the read path -> no output copy.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"

namespace py = pybind11;

namespace meshioplusplus_py {

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

// Holds C-contiguous numpy arrays alive for as long as the C++ mesh views
// into them are in use.
struct PyMeshRefs {
    std::vector<py::array> keep;
};

inline meshioplusplus::NDArray view_from_numpy(const py::array& a) {
    meshioplusplus::DType dt = dtype_from_numpy(a.dtype());
    std::vector<std::size_t> shape(static_cast<std::size_t>(a.ndim()));
    for (py::ssize_t i = 0; i < a.ndim(); ++i)
        shape[static_cast<std::size_t>(i)] = static_cast<std::size_t>(a.shape(i));
    auto* ptr = reinterpret_cast<std::byte*>(const_cast<void*>(a.data()));
    return meshioplusplus::NDArray::make_view(dt, std::move(shape), ptr);
}

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

// Parse a ragged Python cell-block `data` (a list, not an ndarray) into a C++
// CellBlock's ragged members. A "polyhedron" block is 2-level (list of cells,
// each a list of faces, each a sequence of node ids); any other ragged block
// (a "polygon" block with varying node counts) is 1-level (list of cells, each
// a sequence of node ids). This copies (ragged data cannot be a zero-copy
// view).
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

// Python meshioplusplus.Mesh -> C++ meshioplusplus::Mesh (views; zero-copy for rectangular
// blocks). By default ragged list data (polyhedron / jagged polygon blocks)
// raises, so every rectangular-only format writer safely defers such meshes to
// its Python fallback. Pass `allow_ragged=true` (only the ragged-aware format
// bindings do) to instead parse ragged blocks into the CellBlock's ragged
// members (a copy). With `lenient_field_data`, non-numeric field_data entries
// (e.g. MED's "med:nom" list of strings) are silently skipped instead of
// raising; the caller carries them through its own side-channel.
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

// Move an NDArray's buffer into a capsule-backed, writeable numpy array.
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

// Build the Python object a ragged CellBlock maps to: for a jagged polygon
// block, a list of 1-D int64 numpy arrays (one per cell); for a polyhedron
// block, a list of cells, each a list of 1-D int64 numpy arrays (one per
// face). This matches exactly what meshioplusplus.Mesh/CellBlock store for these types
// (kept as a Python list, not converted to a rectangular ndarray).
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

// C++ meshioplusplus::Mesh -> Python meshio.Mesh (zero-copy output arrays for
// rectangular blocks; ragged blocks are copied into Python lists).
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
