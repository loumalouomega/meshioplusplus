// pybind11 entry point for the meshio C++ core.
//
// Imported as ``meshio._core``. Exposes the cell-type tables, exception
// translation to meshio.ReadError/WriteError, and (for now) a roundtrip
// helper used to validate the zero-copy Mesh conversions. Format readers and
// writers are registered here as they are ported.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "meshio/exceptions.hpp"
#include "meshio/formats/abaqus.hpp"
#include "meshio/formats/avsucd.hpp"
#include "meshio/formats/gmsh.hpp"
#include "meshio/formats/medit.hpp"
#include "meshio/formats/nastran.hpp"
#include "meshio/formats/obj_off.hpp"
#include "meshio/formats/ply.hpp"
#include "meshio/formats/stl.hpp"
#include "meshio/formats/su2.hpp"
#include "meshio/formats/tecplot.hpp"
#include "meshio/formats/tetgen.hpp"
#include "meshio/formats/ugrid.hpp"
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

    // VTK writer (version 5.1 or 4.2; ascii or big-endian binary).
    m.def("vtk_write",
          [](const std::string& path, py::object pymesh, bool binary, bool v51) {
              meshio_py::PyMeshRefs refs;
              meshio::Mesh cpp = meshio_py::py_to_mesh(pymesh, refs);
              meshio::write_vtk(path, cpp, binary, v51);
          });

    // VTK 5.1 reader -> Python mesh (zero-copy capsule-backed arrays).
    m.def("vtk_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_vtk(path));
    });

    // STL writer / reader (ascii or binary).
    m.def("stl_write", [](const std::string& path, py::object pymesh, bool binary) {
        meshio_py::PyMeshRefs refs;
        meshio::Mesh cpp = meshio_py::py_to_mesh(pymesh, refs);
        meshio::write_stl(path, cpp, binary);
    });
    m.def("stl_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_stl(path));
    });

    // OFF writer / reader.
    m.def("off_write", [](const std::string& path, py::object pymesh) {
        meshio_py::PyMeshRefs refs;
        meshio::write_off(path, meshio_py::py_to_mesh(pymesh, refs));
    });
    m.def("off_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_off(path));
    });

    // OBJ writer / reader.
    m.def("obj_write", [](const std::string& path, py::object pymesh) {
        meshio_py::PyMeshRefs refs;
        meshio::write_obj(path, meshio_py::py_to_mesh(pymesh, refs));
    });
    m.def("obj_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_obj(path));
    });

    // Gmsh 2.2 writer / reader.
    m.def("gmsh22_write", [](const std::string& path, py::object pymesh, bool binary) {
        meshio_py::PyMeshRefs refs;
        meshio::Mesh cpp = meshio_py::py_to_mesh(pymesh, refs);
        meshio::write_gmsh22(path, cpp, binary);
    });
    m.def("gmsh41_write", [](const std::string& path, py::object pymesh, bool binary) {
        meshio_py::PyMeshRefs refs;
        meshio::Mesh cpp = meshio_py::py_to_mesh(pymesh, refs);
        meshio::write_gmsh41(path, cpp, binary);
    });
    m.def("gmsh_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_gmsh(path));
    });

    // PLY writer / reader (ascii or binary).
    m.def("ply_write", [](const std::string& path, py::object pymesh, bool binary) {
        meshio_py::PyMeshRefs refs;
        meshio::Mesh cpp = meshio_py::py_to_mesh(pymesh, refs);
        meshio::write_ply(path, cpp, binary);
    });
    m.def("ply_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_ply(path));
    });

    // Medit ascii writer / reader (.mesh).
    m.def("medit_write_ascii", [](const std::string& path, py::object pymesh) {
        meshio_py::PyMeshRefs refs;
        meshio::write_medit_ascii(path, meshio_py::py_to_mesh(pymesh, refs));
    });
    m.def("medit_read_ascii", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_medit_ascii(path));
    });

    // Abaqus writer / reader (.inp).
    m.def("abaqus_write", [](const std::string& path, py::object pymesh) {
        meshio_py::PyMeshRefs refs;
        meshio::write_abaqus(path, meshio_py::py_to_mesh(pymesh, refs));
    });
    m.def("abaqus_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_abaqus(path));
    });

    // AVS-UCD writer / reader (.avs).
    m.def("avsucd_write", [](const std::string& path, py::object pymesh) {
        meshio_py::PyMeshRefs refs;
        meshio::write_avsucd(path, meshio_py::py_to_mesh(pymesh, refs));
    });
    m.def("avsucd_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_avsucd(path));
    });

    // Nastran writer / reader (.bdf/.fem/.nas) — meshio-C++ files only.
    m.def("nastran_write", [](const std::string& path, py::object pymesh) {
        meshio_py::PyMeshRefs refs;
        meshio::write_nastran(path, meshio_py::py_to_mesh(pymesh, refs));
    });
    m.def("nastran_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_nastran(path));
    });

    // SU2 writer / reader (.su2).
    m.def("su2_write", [](const std::string& path, py::object pymesh) {
        meshio_py::PyMeshRefs refs;
        meshio::write_su2(path, meshio_py::py_to_mesh(pymesh, refs));
    });
    m.def("su2_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_su2(path));
    });

    // Tecplot writer / reader (.dat/.tec).
    m.def("tecplot_write", [](const std::string& path, py::object pymesh) {
        meshio_py::PyMeshRefs refs;
        meshio::write_tecplot(path, meshio_py::py_to_mesh(pymesh, refs));
    });
    m.def("tecplot_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_tecplot(path));
    });

    // UGRID writer / reader (.ugrid, ascii + binary variants).
    m.def("ugrid_write", [](const std::string& path, py::object pymesh) {
        meshio_py::PyMeshRefs refs;
        meshio::write_ugrid(path, meshio_py::py_to_mesh(pymesh, refs));
    });
    m.def("ugrid_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_ugrid(path));
    });

    // TetGen writer / reader (.node/.ele pair).
    m.def("tetgen_write", [](const std::string& path, py::object pymesh) {
        meshio_py::PyMeshRefs refs;
        meshio::write_tetgen(path, meshio_py::py_to_mesh(pymesh, refs));
    });
    m.def("tetgen_read", [](const std::string& path) {
        return meshio_py::mesh_to_py(meshio::read_tetgen(path));
    });
}
