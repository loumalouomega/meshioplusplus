// pybind11 entry point for the meshio C++ core.
//
// Imported as ``meshio._core``. Exposes the cell-type tables, exception
// translation to meshio.ReadError/WriteError, and (for now) a roundtrip
// helper used to validate the zero-copy Mesh conversions. Format readers and
// writers are registered here as they are ported.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "meshio/exceptions.hpp"
#include "meshio/formats/vtk.hpp"
#include "meshio/formats/vtu.hpp"
#include "meshio/types.hpp"
#include "np_conversions.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = "meshio C++ core (pybind11)";
    m.attr("__cpp_core__") = true;

    // Translate C++ I/O errors to the existing Python exception classes.
    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const meshio::ReadError& e) {
            py::object exc = py::module_::import("meshio").attr("ReadError");
            PyErr_SetString(exc.ptr(), e.what());
        } catch (const meshio::WriteError& e) {
            py::object exc = py::module_::import("meshio").attr("WriteError");
            PyErr_SetString(exc.ptr(), e.what());
        }
    });

    // Shared cell-type metadata (single source of truth with Python).
    m.def("num_nodes_per_cell", []() { return meshio::num_nodes_per_cell(); });
    m.def("topological_dimension", []() { return meshio::topological_dimension(); });

    // Debug helper: Python mesh -> C++ mesh (zero-copy views) -> Python mesh
    // (capsule-backed arrays). Exercises both conversion directions.
    m.def("_roundtrip", [](py::object pymesh) {
        meshio_py::PyMeshRefs refs;
        meshio::Mesh cpp = meshio_py::py_to_mesh(pymesh, refs);
        return meshio_py::mesh_to_py(std::move(cpp));
    });

    // VTU writer (ascii / binary / zlib), zero-copy input from the Python mesh.
    m.def("vtu_write",
          [](const std::string& path, py::object pymesh, bool binary, bool zlib) {
              meshio_py::PyMeshRefs refs;
              meshio::Mesh cpp = meshio_py::py_to_mesh(pymesh, refs);
              meshio::write_vtu(path, cpp, binary, zlib);
          });

    // VTU reader -> Python mesh (zero-copy capsule-backed arrays).
    m.def("vtu_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_vtu(path));
    });

    // VTK 5.1 writer (ascii or big-endian binary).
    m.def("vtk_write_51", [](const std::string& path, py::object pymesh, bool binary) {
        meshio_py::PyMeshRefs refs;
        meshio::Mesh cpp = meshio_py::py_to_mesh(pymesh, refs);
        meshio::write_vtk_51(path, cpp, binary);
    });

    // VTK 5.1 reader -> Python mesh (zero-copy capsule-backed arrays).
    m.def("vtk_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_vtk(path));
    });
}
