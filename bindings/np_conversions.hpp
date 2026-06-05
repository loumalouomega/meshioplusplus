#pragma once
//
// Conversions between meshio::Mesh (C++) and the pure-Python meshio.Mesh,
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

#include "meshio/exceptions.hpp"
#include "meshio/mesh.hpp"

namespace py = pybind11;

namespace meshio_py {

inline meshio::DType dtype_from_numpy(const py::dtype& dt) {
    const char kind = dt.kind();
    const py::ssize_t isz = dt.itemsize();
    using meshio::DType;
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
    throw meshio::WriteError(std::string("Unsupported numpy dtype '") + dt.kind() +
                             std::to_string(isz) + "' for the meshio C++ core");
}

// Holds C-contiguous numpy arrays alive for as long as the C++ mesh views
// into them are in use.
struct PyMeshRefs {
    std::vector<py::array> keep;
};

inline meshio::NDArray view_from_numpy(const py::array& a) {
    meshio::DType dt = dtype_from_numpy(a.dtype());
    std::vector<std::size_t> shape(static_cast<std::size_t>(a.ndim()));
    for (py::ssize_t i = 0; i < a.ndim(); ++i)
        shape[static_cast<std::size_t>(i)] = static_cast<std::size_t>(a.shape(i));
    auto* ptr = reinterpret_cast<std::byte*>(const_cast<void*>(a.data()));
    return meshio::NDArray::make_view(dt, std::move(shape), ptr);
}

inline py::array ensure_contiguous(py::handle obj, PyMeshRefs& refs) {
    py::array a = py::array::ensure(obj, py::array::c_style | py::array::forcecast);
    if (!a) throw meshio::WriteError("Expected an array-like object");
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

// Python meshio.Mesh -> C++ meshio::Mesh (views; zero-copy). Polyhedron cell
// blocks (ragged list data) are skipped here and must be handled on the Python
// side by the caller.
inline meshio::Mesh py_to_mesh(py::handle pymesh, PyMeshRefs& refs) {
    meshio::Mesh m;

    m.points = view_from_numpy(ensure_contiguous(pymesh.attr("points"), refs));

    for (py::handle cb : pymesh.attr("cells")) {
        std::string type = py::cast<std::string>(cb.attr("type"));
        py::object data_obj = py::reinterpret_borrow<py::object>(cb.attr("data"));
        // Ragged polyhedron data is a Python list, not an ndarray.
        if (py::isinstance<py::list>(data_obj)) {
            throw meshio::WriteError(
                "polyhedron cell blocks are not yet handled by the C++ core");
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
        std::vector<meshio::NDArray> blocks;
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
            py::array d = ensure_contiguous(item.second, refs);
            m.field_data.emplace(std::move(name), view_from_numpy(d));
        }
    }

    return m;
}

// Move an NDArray's buffer into a capsule-backed, writeable numpy array.
inline py::array numpy_from_ndarray(meshio::NDArray&& arr) {
    auto* heap = new meshio::NDArray(std::move(arr));
    heap->make_owned();
    py::capsule owner(heap, [](void* p) {
        delete reinterpret_cast<meshio::NDArray*>(p);
    });

    std::vector<py::ssize_t> shape(heap->shape().begin(), heap->shape().end());
    std::vector<py::ssize_t> strides(shape.size());
    const py::ssize_t itemsize = static_cast<py::ssize_t>(meshio::dtype_size(heap->dtype()));
    py::ssize_t s = itemsize;
    for (int i = static_cast<int>(shape.size()) - 1; i >= 0; --i) {
        strides[static_cast<std::size_t>(i)] = s;
        s *= (shape[static_cast<std::size_t>(i)] == 0 ? 1 : shape[static_cast<std::size_t>(i)]);
    }
    return py::array(py::dtype(meshio::dtype_numpy_str(heap->dtype())), shape, strides,
                     heap->data(), owner);
}

// C++ meshio::Mesh -> Python meshio.Mesh (zero-copy output arrays).
inline py::object mesh_to_py(meshio::Mesh&& m) {
    py::object MeshCls = py::module_::import("meshio").attr("Mesh");

    py::array points = numpy_from_ndarray(std::move(m.points));

    py::list cells;
    for (auto& cb : m.cells) {
        cells.append(py::make_tuple(py::str(cb.type), numpy_from_ndarray(std::move(cb.data))));
    }

    py::dict point_data;
    for (auto& kv : m.point_data)
        point_data[py::str(kv.first)] = numpy_from_ndarray(std::move(const_cast<meshio::NDArray&>(kv.second)));

    py::dict cell_data;
    for (auto& kv : m.cell_data) {
        py::list lst;
        for (auto& a : kv.second) lst.append(numpy_from_ndarray(std::move(a)));
        cell_data[py::str(kv.first)] = lst;
    }

    py::dict field_data;
    for (auto& kv : m.field_data)
        field_data[py::str(kv.first)] =
            numpy_from_ndarray(std::move(const_cast<meshio::NDArray&>(kv.second)));

    return MeshCls(points, cells, py::arg("point_data") = point_data,
                   py::arg("cell_data") = cell_data,
                   py::arg("field_data") = field_data);
}

}  // namespace meshio_py
