// License: MIT License, meshio++ default license: LICENSE
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/formats/avsucd.hpp"
#include "meshioplusplus/formats/cgns.hpp"
#include "meshioplusplus/formats/code_aster.hpp"
#include "meshioplusplus/formats/dex.hpp"
#include "meshioplusplus/formats/dolfin.hpp"
#include "meshioplusplus/formats/elmer.hpp"
#include "meshioplusplus/formats/ensight.hpp"
#include "meshioplusplus/formats/febio.hpp"
#include "meshioplusplus/formats/flac3d.hpp"
#include "meshioplusplus/formats/freefem.hpp"
#include "meshioplusplus/formats/femap.hpp"
#include "meshioplusplus/formats/flux.hpp"
#include "meshioplusplus/formats/gid.hpp"
#include "meshioplusplus/formats/gltf.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/ip.hpp"
#include "meshioplusplus/formats/lsdyna.hpp"
#include "meshioplusplus/formats/libmesh.hpp"
#include "meshioplusplus/formats/med.hpp"
#include "meshioplusplus/formats/marc.hpp"
#include "meshioplusplus/formats/mdpa.hpp"
#include "meshioplusplus/formats/mfem.hpp"
#include "meshioplusplus/formats/mff.hpp"
#include "meshioplusplus/formats/mfm.hpp"
#include "meshioplusplus/formats/mphtxt.hpp"
#include "meshioplusplus/formats/obj_off.hpp"
#include "meshioplusplus/formats/nastran.hpp"
#include "meshioplusplus/formats/netgen.hpp"
#include "meshioplusplus/formats/patran.hpp"
#include "meshioplusplus/formats/ply.hpp"
#include "meshioplusplus/formats/pcd.hpp"
#include "meshioplusplus/formats/permas.hpp"
#include "meshioplusplus/formats/radioss.hpp"
#include "meshioplusplus/formats/stl.hpp"
#include "meshioplusplus/formats/su2.hpp"
#include "meshioplusplus/formats/svg.hpp"
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/formats/tetgen.hpp"
#include "meshioplusplus/formats/tikz.hpp"
#include "meshioplusplus/formats/triangle.hpp"
#include "meshioplusplus/formats/unv.hpp"
#include "meshioplusplus/formats/ugrid.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/wkt.hpp"
#include "meshioplusplus/formats/xyz.hpp"
#include "meshioplusplus/formats/z88.hpp"

namespace {
namespace tid = meshioplusplus;
constexpr tid::DType tid_dtypes[] = {tid::DType::Float32, tid::DType::Float64, tid::DType::Int8,
                                     tid::DType::Int16,   tid::DType::Int32,   tid::DType::Int64,
                                     tid::DType::UInt8,   tid::DType::UInt16,  tid::DType::UInt32,
                                     tid::DType::UInt64};

tid::NDArray tid_array(tid::DType Dtype, const std::vector<std::size_t>& rShape,
                       const std::vector<double>& rValues) {
    tid::NDArray out(Dtype, rShape);
    tid::detail::dispatch_dtype(Dtype, [&]<class T>() {
        for (std::size_t i = 0; i < rValues.size(); ++i)
            out.As<T>()[i] = static_cast<T>(rValues[i]);
    });
    return out;
}

/// The mesh each writer is given: a triangle with point and cell data, plus
/// the integer cell data the format reads per cell (tags, types, properties).
/// Z88 has no linear-triangle element, so it gets a tetrahedron.
tid::Mesh tid_mesh(const std::string& rFormat, tid::DType RealType, tid::DType IndexType,
                   std::size_t Dim) {
    tid::Mesh mesh;
    if (rFormat == "tetgen") {
        mesh.AssignPoints(tid_array(RealType, {4, 3}, {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1}));
        mesh.AddCellBlock("tetra", tid_array(IndexType, {1, 4}, {0, 1, 2, 3}));
        mesh.AddPointData("s", tid_array(RealType, {4}, {1, 2, 3, 4}));
        mesh.AddCellData("pf3:ref", {tid_array(IndexType, {1}, {7})});
        return mesh;
    }
    if (rFormat == "z88") {
        mesh.AssignPoints(tid_array(RealType, {4, 3}, {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1}));
        mesh.AddCellBlock("tetra", tid_array(IndexType, {1, 4}, {0, 1, 2, 3}));
        mesh.AddCellData("z88:type", {tid_array(IndexType, {1}, {17})});
        return mesh;
    }
    mesh.AssignPoints(tid_array(RealType, {3, Dim},
                                Dim == 2 ? std::vector<double>{0, 0, 1, 0, 0, 1}
                                         : std::vector<double>{0, 0, 0, 1, 0, 0, 0, 1, 0}));
    mesh.AddCellBlock("triangle", tid_array(IndexType, {1, 3}, {0, 1, 2}));
    mesh.AddPointData("s", tid_array(RealType, {3}, {1, 2, 3}));
    mesh.AddPointData("v", tid_array(RealType, {3, 2}, {1, 2, 3, 4, 5, 6}));
    mesh.AddCellData("pf3:ref", {tid_array(IndexType, {1}, {7})});
    if (rFormat.rfind("gmsh", 0) == 0) {
        mesh.AddCellData("gmsh:physical", {tid_array(IndexType, {1}, {3})});
        mesh.AddCellData("gmsh:geometrical", {tid_array(IndexType, {1}, {5})});
    } else if (rFormat == "mdpa") {
        mesh.AddCellData("gmsh:physical", {tid_array(IndexType, {1}, {3})});
    } else if (rFormat == "patran") {
        mesh.AddCellData("patran:property", {tid_array(IndexType, {1}, {4})});
    } else if (rFormat == "femap") {
        mesh.AddCellData("femap:property", {tid_array(IndexType, {1}, {4})});
    } else if (rFormat == "unv") {
        mesh.AddCellData("unv:pid", {tid_array(IndexType, {1}, {3})});
        mesh.AddCellData("unv:mid", {tid_array(IndexType, {1}, {5})});
    } else if (rFormat == "nastran") {
        mesh.AddCellData("nastran:ref", {tid_array(IndexType, {1}, {4})});
    } else if (rFormat == "radioss") {
        mesh.AddCellData("radioss:part", {tid_array(IndexType, {1}, {2})});
    } else if (rFormat == "mphtxt") {
        mesh.AddCellData("mphtxt:geom", {tid_array(IndexType, {1}, {2})});
    } else if (rFormat == "ansys") {
        mesh.AddCellData("ansys:zone", {tid_array(IndexType, {1}, {2})});
    } else if (rFormat == "su2") {
        mesh.AddCellData("gmsh:physical", {tid_array(IndexType, {1}, {3})});
    }
    return mesh;
}

std::string tid_file_bytes(const std::filesystem::path& rPath) {
    auto in = tid::detail::make_classic_ifstream(rPath.string(), std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

/// A new, empty directory for one writer run, so the files found in it are
/// exactly what that run wrote (a shared temp directory also holds files left
/// by earlier runs).
std::string tid_fresh_dir() {
    const std::filesystem::path dir = mt::temp_path("_dtype");
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    return dir.string();
}

/// Everything under @p rDir (the file, its companions such as `.post.res` or
/// `.ele`, a DOLFIN mesh-function file or, for a directory writer like Elmer,
/// every file inside it), keyed by path relative to @p rDir.
std::map<std::string, std::string> tid_outputs(const std::string& rDir) {
    std::map<std::string, std::string> out;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(rDir))
        if (entry.is_regular_file())
            out[std::filesystem::relative(entry.path(), rDir).string()] =
                tid_file_bytes(entry.path());
    return out;
}

using TidWriter = std::function<void(const std::string&, const tid::Mesh&)>;

/// A writer under test. Formats that record each array's dtype in the file
/// (PCD, the HDF5 formats) or a file name inside it (glTF's .bin) cannot be
/// compared byte for byte with canonical storage; their values are compared
/// through a read-back in the Python tests, and here they must only leave the
/// caller's arrays untouched.
struct TidFormat {
    std::string mName;
    TidWriter mWrite;
    bool mExact = true;
    std::string mSuffix = ".data";
    std::size_t mOnlyDim = 0;  // 0: any; Triangle writes 2-D points only
    /// The output depends on whether arrays are float or integer (UGRID labels,
    /// AVS-UCD materials, DOLFIN mesh functions): compare with the same class.
    bool mClassCanonical = false;
};

const std::vector<TidFormat>& tid_writers() {
    static const std::vector<TidFormat> writers = {
        {"off", tid::write_off},
        {"ip", tid::write_ip},
        {"flux", tid::write_flux},
        {"permas", tid::write_permas},
        {"gmsh22-ascii", [](const std::string& rPath,
                            const tid::Mesh& rMesh) { tid::write_gmsh22(rPath, rMesh, false); }},
        {"gmsh22-binary", [](const std::string& rPath,
                             const tid::Mesh& rMesh) { tid::write_gmsh22(rPath, rMesh, true); }},
        {"gmsh41-ascii", [](const std::string& rPath,
                            const tid::Mesh& rMesh) { tid::write_gmsh41(rPath, rMesh, false); }},
        {"gmsh41-binary", [](const std::string& rPath,
                             const tid::Mesh& rMesh) { tid::write_gmsh41(rPath, rMesh, true); }},
        {"febio", tid::write_febio},
        {"z88",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_z88(rPath, rMesh); }},
        {"patran", tid::write_patran},
        {"femap", tid::write_femap},
        {"mdpa",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_mdpa(rPath, rMesh); }},
        {"libmesh", tid::write_libmesh},
        {"obj", tid::write_obj, true, ".obj"},
        {"stl-binary",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_stl(rPath, rMesh, true);
         },
         true, ".stl"},
        {"stl-ascii",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_stl(rPath, rMesh, false);
         },
         true, ".stl"},
        {"ply-ascii",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_ply(rPath, rMesh, false);
         },
         false, ".ply"},
        {"ply-binary",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_ply(rPath, rMesh, true);
         },
         false, ".ply"},
        {"triangle", tid::write_triangle, true, ".node", 2},
        {"tetgen", tid::write_tetgen, true, ".node"},
        {"ugrid", tid::write_ugrid, true, ".lb8.ugrid", 0, true},
        {"dolfin", tid::write_dolfin, true, ".xml", 0, true},
        {"freefem", tid::write_freefem, true, ".msh"},
        {"avsucd", tid::write_avsucd, true, ".inp", 0, true},
        {"wkt", tid::write_wkt, true, ".wkt"},
        {"svg",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_svg(rPath, rMesh); },
         true, ".svg"},
        {"su2", tid::write_su2, true, ".su2", 0, true},
        {"netgen",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_netgen(rPath, rMesh, ".16e");
         },
         true, ".vol", 0, true},
        {"elmer",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_elmer(rPath, rMesh); },
         true, ""},
        {"mphtxt", tid::write_mphtxt, true, ".mphtxt"},
        {"code_aster", tid::write_code_aster, true, ".mail"},
        {"nastran", tid::write_nastran, true, ".bdf"},
        {"tecplot", tid::write_tecplot, true, ".dat"},
        {"ansys-ascii",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_ansys(rPath, rMesh, false);
         },
         true, ".msh"},
        {"ansys-binary",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_ansys(rPath, rMesh, true);
         },
         true, ".msh"},
        {"radioss",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_radioss(rPath, rMesh); },
         true, ".rad"},
        {"unv",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_unv(rPath, rMesh); },
         true, ".unv"},
        {"dex", tid::write_dex, true, ".dex"},
        {"mfm",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_mfm(rPath, rMesh, ".16e");
         },
         true, ".mfm"},
        {"flac3d-ascii",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_flac3d(rPath, rMesh, ".16e", false);
         },
         true, ".f3grid"},
        {"flac3d-binary",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_flac3d(rPath, rMesh, ".16e", true);
         },
         true, ".f3grid"},
        {"ensight-ascii",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_ensight(rPath, rMesh, false);
         },
         true, ".case"},
        {"ensight-binary",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_ensight(rPath, rMesh, true);
         },
         true, ".case"},
        {"marc", tid::write_marc, true, ".dat"},
        {"lsdyna", tid::write_lsdyna, true, ".k"},
        {"mfem",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_mfem(rPath, rMesh); },
         true, ".mesh"},
        {"mff", tid::write_mff, true, ".mff"},
        {"xyz",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_xyz(rPath, rMesh); },
         true, ".xyz"},
        {"vtp",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_vtp(rPath, rMesh, true, true);
         },
         false, ".vtp"},
        {"tikz",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_tikz(rPath, rMesh); },
         true, ".tex"},
        {"gid",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_gid(rPath, rMesh, tid::GidMode::Ascii);
         },
         true, ".post.msh"},
        {"pcd-ascii",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_pcd(rPath, rMesh, tid::PcdData::Ascii);
         },
         false},
        {"pcd-binary",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_pcd(rPath, rMesh, tid::PcdData::Binary);
         },
         false},
        {"gltf",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_gltf(rPath, rMesh); },
         false, ".gltf"},
#ifdef MESHIOPLUSPLUS_HAS_HDF5
        {"cgns",
         [](const std::string& rPath, const tid::Mesh& rMesh) { tid::write_cgns(rPath, rMesh, 0); },
         false, ".cgns"},
        {"med",
         [](const std::string& rPath, const tid::Mesh& rMesh) {
             tid::write_med(rPath, rMesh, tid::MedInfo{});
         },
         false, ".med"},
#endif
    };
    return writers;
}
}  // namespace

TEST(TextIoDtypes, HoistedWritersMatchCanonicalStorageOnAllDtypes) {
    for (const auto& fmt : tid_writers()) {
        if (!fmt.mExact)
            continue;
        const std::string& format = fmt.mName;
        const auto& writer = fmt.mWrite;
        SCOPED_TRACE(format);
        for (const auto dtype : tid_dtypes) {
            SCOPED_TRACE(static_cast<int>(dtype));
            for (const std::size_t dim : {2u, 3u}) {
                if (fmt.mOnlyDim != 0 && fmt.mOnlyDim != dim)
                    continue;
                SCOPED_TRACE(dim);
                const auto input = tid_mesh(format, dtype, dtype, dim);
                const bool is_float = dtype == tid::DType::Float32 || dtype == tid::DType::Float64;
                const tid::DType real =
                    fmt.mClassCanonical && !is_float ? tid::DType::Int64 : tid::DType::Float64;
                const tid::DType index =
                    fmt.mClassCanonical && is_float ? tid::DType::Float64 : tid::DType::Int64;
                const auto expected = tid_mesh(format, real, index, dim);
                const std::string first_dir = tid_fresh_dir();
                const std::string second_dir = tid_fresh_dir();
                writer(first_dir + "/m" + fmt.mSuffix, input);
                writer(second_dir + "/m" + fmt.mSuffix, expected);
                const auto first_files = tid_outputs(first_dir);
                EXPECT_FALSE(first_files.empty());
                EXPECT_EQ(first_files, tid_outputs(second_dir));
                std::filesystem::remove_all(first_dir);
                std::filesystem::remove_all(second_dir);
            }
        }
    }
}

TEST(TextIoDtypes, HoistedWritersLeaveCallerArraysUntouched) {
    for (const auto& fmt : tid_writers()) {
        const std::string& format = fmt.mName;
        const auto& writer = fmt.mWrite;
        SCOPED_TRACE(format);
        for (const auto dtype : {tid::DType::Float32, tid::DType::Int32, tid::DType::UInt16}) {
            SCOPED_TRACE(static_cast<int>(dtype));
            const std::size_t dim = fmt.mOnlyDim != 0 ? fmt.mOnlyDim : 3;
            const auto input = tid_mesh(format, dtype, dtype, dim);
            const tid::NDArray& points = input.Points();
            const tid::NDArray& conn = input.Cells(0).Conn();
            const std::vector<std::byte> points_before(points.Data(),
                                                       points.Data() + points.Nbytes());
            const std::vector<std::byte> conn_before(conn.Data(), conn.Data() + conn.Nbytes());
            const std::string dir = tid_fresh_dir();
            writer(dir + "/m" + fmt.mSuffix, input);
            std::filesystem::remove_all(dir);
            EXPECT_EQ(std::vector<std::byte>(points.Data(), points.Data() + points.Nbytes()),
                      points_before);
            EXPECT_EQ(std::vector<std::byte>(conn.Data(), conn.Data() + conn.Nbytes()),
                      conn_before);
        }
    }
}
