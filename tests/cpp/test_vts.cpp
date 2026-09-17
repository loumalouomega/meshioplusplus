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
// VTK XML StructuredGrid (roadmap §1 tier B4, v11.6.0).

// System includes
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vts.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::read_vts;
using meshioplusplus::read_vts_metadata;
using meshioplusplus::ReadError;
using meshioplusplus::write_vts;
using meshioplusplus::WriteError;

namespace {

Mesh vts_lattice(std::int64_t N, double Origin, double Spacing) {
    return meshioplusplus::grid({{N, N, N}}, {{Origin, Origin, Origin}},
                                {{Spacing, Spacing, Spacing}});
}

}  // namespace

TEST(Vts, RoundTripsALatticeExactly) {
    for (bool binary : {false, true}) {
        Mesh m = vts_lattice(3, -0.5, 0.25);
        NDArray pd = NDArray::Uninit(meshioplusplus::DType::Float64, {m.NumPoints()});
        for (std::size_t p = 0; p < m.NumPoints(); ++p)
            pd.As<double>()[p] = static_cast<double>(p) * 0.5;
        m.AddPointData("f", std::move(pd));
        NDArray cd = NDArray::Uninit(meshioplusplus::DType::Int64, {27u});
        for (std::size_t c = 0; c < 27; ++c)
            cd.As<std::int64_t>()[c] = static_cast<std::int64_t>(c);
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(cd));
        m.AddCellData("g", std::move(blocks));

        const std::string path = mt::temp_path(".vts");
        write_vts(path, m, binary, /*zlib=*/false);
        const Mesh back = read_vts(path);

        ASSERT_EQ(back.NumPoints(), m.NumPoints());
        for (std::size_t p = 0; p < m.NumPoints(); ++p)
            for (std::size_t c = 0; c < 3; ++c)
                EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.Points(), p * 3 + c),
                                 meshioplusplus::detail::read_double(m.Points(), p * 3 + c));
        ASSERT_EQ(back.NumCellBlocks(), 1u);
        EXPECT_EQ(back.Cells(0).Type(), "hexahedron");
        EXPECT_EQ(back.Cells(0).NumCells(), 27u);  // 3x3x3 cells
        ASSERT_TRUE(back.HasPointData("f"));
        for (std::size_t p = 0; p < m.NumPoints(); ++p)
            EXPECT_DOUBLE_EQ(back.PointData("f").As<double>()[p], static_cast<double>(p) * 0.5);
        ASSERT_TRUE(back.HasCellData("g"));
        for (std::size_t c = 0; c < 27; ++c)
            EXPECT_EQ(back.CellData("g", 0).As<std::int64_t>()[c], static_cast<std::int64_t>(c));
        std::remove(path.c_str());
    }
}

TEST(Vts, RefusesAMeshThatIsNotALattice) {
    Mesh m = mt::tet_mesh();
    const std::string path = mt::temp_path(".vts");
    EXPECT_THROW(write_vts(path, m, false, false), WriteError);
}

TEST(Vts, MetadataAgreesWithARealRead) {
    Mesh m = vts_lattice(2, 0.0, 1.0);
    const std::string path = mt::temp_path(".vts");
    write_vts(path, m, true, false);

    const meshioplusplus::MeshMetadata meta = read_vts_metadata(path);
    const Mesh back = read_vts(path);
    EXPECT_EQ(meta.mNumPoints, back.NumPoints());
    ASSERT_EQ(meta.mCellBlocks.size(), 1u);
    EXPECT_EQ(meta.mCellBlocks[0].mType, "hexahedron");
    EXPECT_EQ(meta.mCellBlocks[0].mNumCells, back.Cells(0).NumCells());
    std::remove(path.c_str());
}

// A hand-written file with explicit, non-lattice-checked points: the reader
// never verifies uniform spacing (only the writer does, via
// detail::lattice_from_mesh), so a genuinely curved StructuredGrid -- the
// point of the format over .vti -- must still read correctly. One cell
// (WholeExtent 0-1 on every axis, 8 points), the last point moved off the
// regular grid -- something .vti's implicit Origin/Spacing could not state
// at all, and something this scope deliberately does not verify on read
// (only .vts's WRITER requires a lattice; degenerate 2-D/1-D extents are a
// documented remainder, matching .vti's own scope, so this fixture stays
// genuinely 3-D).
TEST(Vts, ReadsAHandWrittenFileWithNonUniformPoints) {
    const std::string path = mt::temp_path(".vts");
    std::ofstream os(path, std::ios::binary);
    os << "<?xml version=\"1.0\"?>\n"
          "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
          "<StructuredGrid WholeExtent=\"0 1 0 1 0 1\">\n"
          "<Piece Extent=\"0 1 0 1 0 1\">\n"
          "<Points>\n"
          "<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n"
          "0 0 0\n1 0 0\n0 1 0\n1 1 0\n0 0 1\n1 0 1\n0 1 1\n1.5 1.5 1.5\n"
          "</DataArray>\n"
          "</Points>\n"
          "</Piece>\n"
          "</StructuredGrid>\n"
          "</VTKFile>\n";
    os.close();

    const Mesh back = read_vts(path);
    EXPECT_EQ(back.NumPoints(), 8u);
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).NumCells(), 1u);
    // The last point is NOT on a regular grid -- exactly what .vti could not
    // express and .vts can.
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.Points(), 7 * 3 + 0), 1.5);
    std::remove(path.c_str());
}
