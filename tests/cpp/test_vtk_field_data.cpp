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
// `<FieldData>` on VTU / VTP (v15.0.0): mesh-level `field_data` round-trips, and the
// layouts VTK itself writes (grid-level, piece-level, `TimeValue`) are read.

// System includes
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtp.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/read_options.hpp"

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::ReadOptions;

namespace {

std::string fd_slurp(const std::string& rPath) {
    std::ifstream in(rPath);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

NDArray fd_scalar(double V) {
    NDArray a(DType::Float64, {1});
    *a.As<double>() = V;
    return a;
}

NDArray fd_matrix(std::size_t Rows, std::size_t Cols, DType Dt = DType::Float64) {
    NDArray a(Dt, {Rows, Cols});
    for (std::size_t i = 0; i < Rows * Cols; ++i) {
        if (Dt == DType::Float64)
            a.As<double>()[i] = 0.5 * static_cast<double>(i);
        else
            a.As<std::int32_t>()[i] = static_cast<std::int32_t>(i);
    }
    return a;
}

Mesh fd_mesh() {
    Mesh m = mt::tri_mesh();
    m.AddFieldData("TimeValue", fd_scalar(0.25));
    m.AddFieldData("params", fd_matrix(2, 3));
    m.AddFieldData("ids", fd_matrix(4, 1, DType::Int32));
    return m;
}

void fd_expect_same(const Mesh& rBack, const Mesh& rMesh) {
    EXPECT_EQ(rBack.FieldDataNames(), rMesh.FieldDataNames());
    for (const std::string& name : rMesh.FieldDataNames()) {
        const NDArray& a = rMesh.FieldData(name);
        const NDArray& b = rBack.FieldData(name);
        EXPECT_EQ(b.Dtype(), a.Dtype()) << name;
        EXPECT_EQ(b.Size(), a.Size()) << name;
        ASSERT_EQ(b.Nbytes(), a.Nbytes()) << name;
        EXPECT_EQ(std::memcmp(a.Data(), b.Data(), a.Nbytes()), 0) << name;
    }
}

std::string fd_temp(const std::string& rSuffix) {
    return mt::temp_path(rSuffix);
}

}  // namespace

TEST(VtkFieldData, VtuRoundTripsNumericFieldDataInEveryEncoding) {
    const Mesh m = fd_mesh();
    for (const bool binary : {false, true})
        for (const bool zlib : {false, true}) {
            if (!binary && zlib)
                continue;
            const std::string path = fd_temp(".vtu");
            meshioplusplus::write_vtu(path, m, binary, zlib);
            const Mesh back = meshioplusplus::read_vtu(path);
            fd_expect_same(back, m);
            std::filesystem::remove(path);
        }
}

TEST(VtkFieldData, VtpRoundTripsNumericFieldDataInEveryEncoding) {
    const Mesh m = [] {
        Mesh q = mt::quad_mesh();
        q.AddFieldData("TimeValue", fd_scalar(1.5));
        q.AddFieldData("params", fd_matrix(2, 3));
        return q;
    }();
    for (const bool binary : {false, true}) {
        const std::string path = fd_temp(".vtp");
        meshioplusplus::write_vtp(path, m, binary, false);
        const Mesh back = meshioplusplus::read_vtp(path);
        fd_expect_same(back, m);
        std::filesystem::remove(path);
    }
}

TEST(VtkFieldData, TheWrittenLayoutIsWhatVtkReads) {
    const Mesh m = fd_mesh();
    const std::string path = fd_temp(".vtu");
    meshioplusplus::write_vtu(path, m, false, false);
    const std::string text = fd_slurp(path);
    // On the grid, before the piece, one array per name, with the tuple count VTK requires.
    const std::size_t fd = text.find("<FieldData>");
    const std::size_t piece = text.find("<Piece ");
    ASSERT_NE(fd, std::string::npos);
    EXPECT_LT(fd, piece);
    EXPECT_NE(text.find("Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\""),
              std::string::npos)
        << text;
    EXPECT_NE(text.find("Name=\"params\" NumberOfTuples=\"2\" NumberOfComponents=\"3\""),
              std::string::npos)
        << text;
    std::filesystem::remove(path);
}

TEST(VtkFieldData, AMeshWithoutFieldDataWritesNoFieldDataElement) {
    // The guard that keeps every existing file byte-identical (and the pinned
    // baseline hashes valid).
    const Mesh m = mt::tri_mesh();
    const std::string vtu = fd_temp(".vtu"), vtp = fd_temp(".vtp");
    meshioplusplus::write_vtu(vtu, m, false, false);
    meshioplusplus::write_vtp(vtp, mt::quad_mesh(), false, false);
    EXPECT_EQ(fd_slurp(vtu).find("FieldData"), std::string::npos);
    EXPECT_EQ(fd_slurp(vtp).find("FieldData"), std::string::npos);
    std::filesystem::remove(vtu);
    std::filesystem::remove(vtp);
}

namespace {

// A minimal ASCII .vtu with `pGridFd` on the grid and `pPieceFd` inside the piece.
std::string fd_vtu(const char* pGridFd, const char* pPieceFd) {
    return std::string(
               "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" "
               "version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n") +
           pGridFd + "<Piece NumberOfPoints=\"3\" NumberOfCells=\"1\">\n" + pPieceFd +
           "<Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
           "0 0 0 1 0 0 0 1 0</DataArray></Points>\n<Cells>"
           "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">0 1 2</DataArray>"
           "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">3</DataArray>"
           "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">5</DataArray>"
           "</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

std::string fd_time_array(const char* pName, const char* pValue) {
    return std::string("<FieldData><DataArray type=\"Float64\" Name=\"") + pName +
           "\" NumberOfTuples=\"1\" format=\"ascii\">" + pValue + "</DataArray></FieldData>\n";
}

}  // namespace

TEST(VtkFieldData, ReadsVtksTimeValueOnTheGridAndInsideAPieceThePieceWinning) {
    const std::string path = fd_temp(".vtu");
    {
        std::ofstream out(path);
        out << fd_vtu(fd_time_array("TimeValue", "0.25").c_str(),
                      fd_time_array("Other", "9").c_str());
    }
    const Mesh grid_only = [&] {
        std::ofstream(path) << fd_vtu(fd_time_array("TimeValue", "0.25").c_str(), "");
        return meshioplusplus::read_vtu(path);
    }();
    ASSERT_TRUE(grid_only.HasFieldData("TimeValue"));
    EXPECT_EQ(*grid_only.FieldData("TimeValue").As<double>(), 0.25);

    std::ofstream(path) << fd_vtu(fd_time_array("TimeValue", "0.25").c_str(),
                                  fd_time_array("TimeValue", "0.75").c_str());
    const Mesh both = meshioplusplus::read_vtu(path);
    EXPECT_EQ(*both.FieldData("TimeValue").As<double>(), 0.75);  // the piece's overrides
    std::filesystem::remove(path);
}

TEST(VtkFieldData, ANonNumericArrayIsSkippedNotFatal) {
    // `type="String"` (a vtkStringArray) used to be ignored along with the whole
    // section; reading it must not start failing.
    const std::string path = fd_temp(".vtu");
    const std::string fd =
        "<FieldData><DataArray type=\"String\" Name=\"label\" NumberOfTuples=\"1\" "
        "format=\"ascii\">bracket</DataArray>"
        "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" "
        "format=\"ascii\">2</DataArray></FieldData>\n";
    std::ofstream(path) << fd_vtu(fd.c_str(), "");
    const Mesh m = meshioplusplus::read_vtu(path);
    EXPECT_EQ(m.FieldDataNames(), std::vector<std::string>({"TimeValue"}));
    EXPECT_EQ(meshioplusplus::read_vtu_metadata(path).mFieldDataNames,
              std::vector<std::string>({"TimeValue"}));
    std::filesystem::remove(path);
}

TEST(VtkFieldData, MetadataAgreesWithARealRead) {
    const Mesh m = fd_mesh();
    const std::string vtu = fd_temp(".vtu"), vtp = fd_temp(".vtp");
    meshioplusplus::write_vtu(vtu, m, true, false);
    EXPECT_EQ(meshioplusplus::read_vtu_metadata(vtu).mFieldDataNames,
              meshioplusplus::read_vtu(vtu).FieldDataNames());
    Mesh q = mt::quad_mesh();
    q.AddFieldData("TimeValue", fd_scalar(1.0));
    meshioplusplus::write_vtp(vtp, q, true, false);
    EXPECT_EQ(meshioplusplus::read_vtp_metadata(vtp).mFieldDataNames,
              meshioplusplus::read_vtp(vtp).FieldDataNames());
    std::filesystem::remove(vtu);
    std::filesystem::remove(vtp);
}

TEST(VtkFieldData, SelectiveReadsNarrowFieldDataLikeEveryOtherData) {
    const Mesh m = fd_mesh();
    const std::string path = fd_temp(".vtu");
    meshioplusplus::write_vtu(path, m, false, false);

    ReadOptions only;
    only.mDataArrays = std::vector<std::string>{"TimeValue"};
    EXPECT_EQ(meshioplusplus::read_vtu(path, only).FieldDataNames(),
              std::vector<std::string>({"TimeValue"}));

    ReadOptions none;
    none.mPointsOnly = true;
    EXPECT_EQ(meshioplusplus::read_vtu(path, none).NumFieldData(), 0u);
    std::filesystem::remove(path);
}
