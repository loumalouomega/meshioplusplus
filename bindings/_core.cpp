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
// External includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/abaqus.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/ansysinp.hpp"
#include "meshioplusplus/formats/avsucd.hpp"
#ifdef MESHIOPLUSPLUS_HAS_HDF5
#include "meshioplusplus/formats/cgns.hpp"
#include "meshioplusplus/formats/h5m.hpp"
#include "meshioplusplus/formats/hmf.hpp"
#include "meshioplusplus/formats/med.hpp"
#endif
#include "meshioplusplus/formats/dolfin.hpp"
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
#include "meshioplusplus/formats/exodus.hpp"
#endif
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
#include "meshioplusplus/formats/wkt.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/formats/xdmf.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"
#include "np_conversions.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = "meshio++ C++ core (pybind11)";
    m.attr("__cpp_core__") = true;

    // Build-time capability flags: the HDF5/netCDF formats are optional and
    // compile out when the libraries are absent (Python is the fallback).
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    m.attr("__has_hdf5__") = true;
#else
    m.attr("__has_hdf5__") = false;
#endif
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
    m.attr("__has_netcdf__") = true;
#else
    m.attr("__has_netcdf__") = false;
#endif
    // Active compile-time parallel backend ("seq"/"stl"/"openmp"/"tbb"): lets
    // callers verify that parallel_for actually threads (STL without TBB is
    // effectively sequential).
    m.attr("__parallel_backend__") = meshioplusplus::parallel_backend_name();

    // Translate C++ I/O errors to the existing Python exception classes.
    py::register_exception_translator([](std::exception_ptr p) {
        try {
            if (p) std::rethrow_exception(p);
        } catch (const meshioplusplus::ReadError& e) {
            py::object exc = py::module_::import("meshioplusplus").attr("ReadError");
            PyErr_SetString(exc.ptr(), e.what());
        } catch (const meshioplusplus::WriteError& e) {
            py::object exc = py::module_::import("meshioplusplus").attr("WriteError");
            PyErr_SetString(exc.ptr(), e.what());
        }
    });

    // Shared cell-type metadata (single source of truth with Python).
    m.def("num_nodes_per_cell", []() { return meshioplusplus::num_nodes_per_cell(); });
    m.def("topological_dimension", []() { return meshioplusplus::topological_dimension(); });

    // Debug helper: Python mesh -> C++ mesh (zero-copy views) -> Python mesh
    // (capsule-backed arrays). Exercises both conversion directions. With
    // allow_ragged=True it also round-trips polyhedron / jagged-polygon blocks
    // through the C++ ragged CellBlock representation (a copy).
    m.def(
        "_roundtrip",
        [](py::object pymesh, bool allow_ragged) {
            meshioplusplus_py::PyMeshRefs refs;
            meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs,
                                                     /*lenient_field_data=*/false,
                                                     allow_ragged);
            return meshioplusplus_py::mesh_to_py(std::move(cpp));
        },
        py::arg("pymesh"), py::arg("allow_ragged") = false);

    // VTU writer (ascii / binary / zlib), zero-copy input from the Python mesh.
    m.def("vtu_write",
          [](const std::string& path, py::object pymesh, bool binary, bool zlib) {
              meshioplusplus_py::PyMeshRefs refs;
              meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
              meshioplusplus::write_vtu(path, cpp, binary, zlib);
          });

    // VTU reader -> Python mesh (zero-copy capsule-backed arrays).
    m.def("vtu_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_vtu(path));
    });

    // VTK writer (version 5.1 or 4.2; ascii or big-endian binary).
    m.def("vtk_write",
          [](const std::string& path, py::object pymesh, bool binary, bool v51) {
              meshioplusplus_py::PyMeshRefs refs;
              meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
              meshioplusplus::write_vtk(path, cpp, binary, v51);
          });

    // VTK 5.1 reader -> Python mesh (zero-copy capsule-backed arrays).
    m.def("vtk_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_vtk(path));
    });

    // STL writer / reader (ascii or binary).
    m.def("stl_write", [](const std::string& path, py::object pymesh, bool binary) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
        meshioplusplus::write_stl(path, cpp, binary);
    });
    m.def("stl_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_stl(path));
    });

    // OFF writer / reader.
    m.def("off_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_off(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("off_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_off(path));
    });

    // OBJ writer / reader.
    m.def("obj_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_obj(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("obj_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_obj(path));
    });

    // Gmsh 2.2 writer / reader.
    m.def("gmsh22_write", [](const std::string& path, py::object pymesh, bool binary) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
        meshioplusplus::write_gmsh22(path, cpp, binary);
    });
    m.def("gmsh41_write", [](const std::string& path, py::object pymesh, bool binary) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
        meshioplusplus::write_gmsh41(path, cpp, binary);
    });
    m.def("gmsh_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_gmsh(path));
    });

    // PLY writer / reader (ascii or binary).
    m.def("ply_write", [](const std::string& path, py::object pymesh, bool binary) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
        meshioplusplus::write_ply(path, cpp, binary);
    });
    m.def("ply_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_ply(path));
    });

    // Medit ascii writer / reader (.mesh).
    m.def("medit_write_ascii", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_medit_ascii(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("medit_read_ascii", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_medit_ascii(path));
    });

    // Abaqus writer / reader (.inp).
    m.def("abaqus_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_abaqus(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("abaqus_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_abaqus(path));
    });

    // AVS-UCD writer / reader (.avs).
    m.def("avsucd_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_avsucd(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("avsucd_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_avsucd(path));
    });

    // Nastran writer / reader (.bdf/.fem/.nas) — meshio++-C++ files only.
    m.def("nastran_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_nastran(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("nastran_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_nastran(path));
    });

    // SU2 writer / reader (.su2).
    m.def("su2_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_su2(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("su2_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_su2(path));
    });

    // Tecplot writer / reader (.dat/.tec).
    m.def("tecplot_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_tecplot(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("tecplot_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_tecplot(path));
    });

    // UGRID writer / reader (.ugrid, ascii + binary variants).
    m.def("ugrid_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_ugrid(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("ugrid_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_ugrid(path));
    });

    // UNV (I-DEAS Universal) writer / reader (.unv).
    m.def("unv_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_unv(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("unv_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_unv(path));
    });

    // TetGen writer / reader (.node/.ele pair).
    m.def("tetgen_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_tetgen(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("tetgen_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_tetgen(path));
    });

    // XDMF writer / reader (.xdmf/.xmf) — XML/Binary always; HDF when built
    // with HDF5.
    m.def("xdmf_write",
          [](const std::string& path, py::object pymesh, const std::string& data_format,
             int gzip_level) {
              meshioplusplus_py::PyMeshRefs refs;
              meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
              meshioplusplus::write_xdmf(path, cpp, data_format, gzip_level);
          },
          py::arg("path"), py::arg("mesh"), py::arg("data_format"),
          py::arg("gzip_level") = -1);
    m.def("xdmf_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_xdmf(path));
    });

#ifdef MESHIOPLUSPLUS_HAS_HDF5
    // CGNS writer / reader (.cgns).
    m.def("cgns_write",
          [](const std::string& path, py::object pymesh, int gzip_level) {
              meshioplusplus_py::PyMeshRefs refs;
              meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
              meshioplusplus::write_cgns(path, cpp, gzip_level);
          });
    m.def("cgns_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_cgns(path));
    });

    // HMF writer / reader (.hmf).
    m.def("hmf_write",
          [](const std::string& path, py::object pymesh, int gzip_level) {
              meshioplusplus_py::PyMeshRefs refs;
              meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
              meshioplusplus::write_hmf(path, cpp, gzip_level);
          });
    m.def("hmf_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_hmf(path));
    });

    // MOAB h5m writer / reader (.h5m).
    m.def("h5m_write",
          [](const std::string& path, py::object pymesh, bool add_global_ids,
             int gzip_level) {
              meshioplusplus_py::PyMeshRefs refs;
              meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
              meshioplusplus::write_h5m(path, cpp, add_global_ids, gzip_level);
          });
    m.def("h5m_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_h5m(path));
    });

    // MED/Salome writer / reader (.med). point_tags/cell_tags are custom Mesh
    // attributes and med:nom is a list of string-lists, so they travel outside
    // the Mesh conversion layer.
    m.def(
        "med_write",
        [](const std::string& path, py::object pymesh,
           std::map<std::int64_t, std::vector<std::string>> point_tags,
           std::map<std::int64_t, std::vector<std::string>> cell_tags,
           std::vector<std::vector<std::string>> med_nom, std::string mesh_name,
           std::string description, std::string unit_time, std::string unit_coords,
           std::map<std::int64_t, std::string> point_tag_groups,
           std::map<std::int64_t, std::string> cell_tag_groups,
           std::string med_version) {
            meshioplusplus_py::PyMeshRefs refs;
            meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs,
                                                     /*lenient_field_data=*/true,
                                                     /*allow_ragged=*/true);
            meshioplusplus::MedInfo info;
            info.point_tags = std::move(point_tags);
            info.cell_tags = std::move(cell_tags);
            info.med_nom = std::move(med_nom);
            info.mesh_name = mesh_name.empty() ? "mesh" : std::move(mesh_name);
            info.description = std::move(description);
            info.unit_time = std::move(unit_time);
            info.unit_coords = std::move(unit_coords);
            info.point_tag_groups = std::move(point_tag_groups);
            info.cell_tag_groups = std::move(cell_tag_groups);
            meshioplusplus::write_med(path, cpp, info, med_version);
        });
    m.def("med_read", [](const std::string& path) {
        meshioplusplus::MedInfo info;
        py::object pymesh = meshioplusplus_py::mesh_to_py(meshioplusplus::read_med(path, info));
        py::dict ptags, ctags, pgroups, cgroups;
        for (const auto& kv : info.point_tags) ptags[py::int_(kv.first)] = kv.second;
        for (const auto& kv : info.cell_tags) ctags[py::int_(kv.first)] = kv.second;
        for (const auto& kv : info.point_tag_groups)
            pgroups[py::int_(kv.first)] = kv.second;
        for (const auto& kv : info.cell_tag_groups)
            cgroups[py::int_(kv.first)] = kv.second;
        pymesh.attr("point_tags") = ptags;
        pymesh.attr("cell_tags") = ctags;
        pymesh.attr("point_tag_groups") = pgroups;
        pymesh.attr("cell_tag_groups") = cgroups;
        pymesh.attr("mesh_name") = info.mesh_name;
        pymesh.attr("description") = info.description;
        pymesh.attr("unit_time") = info.unit_time;
        pymesh.attr("unit_coords") = info.unit_coords;
        if (!info.med_nom.empty())
            pymesh.attr("field_data")[py::str("med:nom")] = py::cast(info.med_nom);
        return pymesh;
    });
#endif

#ifdef MESHIOPLUSPLUS_HAS_NETCDF
    // Exodus II writer / reader (.e/.exo/.ex2).
    m.def("exodus_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_exodus(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("exodus_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_exodus(path));
    });
#endif

    // DOLFIN XML writer / reader (.xml).
    m.def("dolfin_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_dolfin(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("dolfin_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_dolfin(path));
    });

    // Ansys/Fluent writer / reader (.msh, ascii + binary).
    m.def("ansys_write", [](const std::string& path, py::object pymesh, bool binary) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
        meshioplusplus::write_ansys(path, cpp, binary);
    });
    m.def("ansys_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_ansys(path));
    });

    // Ansys MAPDL coded database (.cdb/.inp). CMBLOCK components are point_sets
    // / cell_sets, custom Mesh attributes carried through the AnsysInfo
    // side-channel.
    m.def("ansysinp_write",
          [](const std::string& path, py::object pymesh,
             std::map<std::string, std::vector<std::int64_t>> point_sets,
             std::map<std::string, std::vector<std::vector<std::int64_t>>> cell_sets) {
              meshioplusplus_py::PyMeshRefs refs;
              meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
              meshioplusplus::AnsysInfo info;
              info.point_sets = std::move(point_sets);
              info.cell_sets = std::move(cell_sets);
              meshioplusplus::write_ansysinp(path, cpp, info);
          });
    m.def("ansysinp_read", [](const std::string& path) {
        meshioplusplus::AnsysInfo info;
        py::object pymesh = meshioplusplus_py::mesh_to_py(meshioplusplus::read_ansysinp(path, info));
        py::dict psets, csets;
        for (const auto& kv : info.point_sets)
            psets[py::str(kv.first)] = py::array_t<std::int64_t>(
                static_cast<py::ssize_t>(kv.second.size()), kv.second.data());
        for (const auto& kv : info.cell_sets) {
            py::list blocks;
            for (const auto& blk : kv.second)
                blocks.append(py::array_t<std::int64_t>(
                    static_cast<py::ssize_t>(blk.size()), blk.data()));
            csets[py::str(kv.first)] = blocks;
        }
        pymesh.attr("point_sets") = psets;
        pymesh.attr("cell_sets") = csets;
        return pymesh;
    });

    // OpenFOAM polyMesh reader (read-only). Boundary patch names are
    // mesh.cell_tags, carried through the OpenFoamInfo side-channel.
    m.def("openfoam_read", [](const std::string& path) {
        meshioplusplus::OpenFoamInfo info;
        py::object pymesh = meshioplusplus_py::mesh_to_py(meshioplusplus::read_openfoam(path, info));
        py::dict ctags;
        for (const auto& kv : info.cell_tags) ctags[py::int_(kv.first)] = kv.second;
        pymesh.attr("cell_tags") = ctags;
        pymesh.attr("point_tags") = py::dict();
        return pymesh;
    });

    // WKT (TIN) writer / reader (.wkt).
    m.def("wkt_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_wkt(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("wkt_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_wkt(path));
    });

    // PERMAS writer / reader (.post/.dato).
    m.def("permas_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_permas(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("permas_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_permas(path));
    });

    // FLAC3D writer / reader (.f3grid, ascii + binary, common path).
    m.def("flac3d_write",
          [](const std::string& path, py::object pymesh, const std::string& float_fmt,
             bool binary) {
              meshioplusplus_py::PyMeshRefs refs;
              meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
              meshioplusplus::write_flac3d(path, cpp, float_fmt, binary);
          });
    m.def("flac3d_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_flac3d(path));
    });

    // FLUX .pf3 writer / reader.
    m.def("flux_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_flux(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("flux_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_flux(path));
    });

    // COMSOL .mphtxt writer / reader.
    m.def("mphtxt_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_mphtxt(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("mphtxt_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_mphtxt(path));
    });

    // FreeFem++ writer / reader (.msh).
    m.def("freefem_write", [](const std::string& path, py::object pymesh) {
        meshioplusplus_py::PyMeshRefs refs;
        meshioplusplus::write_freefem(path, meshioplusplus_py::py_to_mesh(pymesh, refs));
    });
    m.def("freefem_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_freefem(path));
    });

    // MFM (Modulef Formatted Mesh) writer / reader (.mfm).
    m.def("mfm_write",
          [](const std::string& path, py::object pymesh, const std::string& float_fmt) {
              meshioplusplus_py::PyMeshRefs refs;
              meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
              meshioplusplus::write_mfm(path, cpp, float_fmt);
          });
    m.def("mfm_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_mfm(path));
    });

    // Netgen writer / reader (.vol, common path).
    m.def("netgen_write",
          [](const std::string& path, py::object pymesh, const std::string& float_fmt) {
              meshioplusplus_py::PyMeshRefs refs;
              meshioplusplus::Mesh cpp = meshioplusplus_py::py_to_mesh(pymesh, refs);
              meshioplusplus::write_netgen(path, cpp, float_fmt);
          });
    m.def("netgen_read", [](const std::string& path) {
        return meshioplusplus_py::mesh_to_py(meshioplusplus::read_netgen(path));
    });
}
