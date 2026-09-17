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
// VTK XML RectilinearGrid (roadmap §1 tier B4, v11.6.0).

// System includes
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtr.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::read_vtr;
using meshioplusplus::read_vtr_metadata;
using meshioplusplus::ReadError;
using meshioplusplus::write_vtr;
using meshioplusplus::WriteError;

namespace {

Mesh vtr_lattice(std::int64_t N, double Origin, double Spacing) {
    return meshioplusplus::grid({{N, N, N}}, {{Origin, Origin, Origin}},
                                {{Spacing, Spacing, Spacing}});
}

}  // namespace

TEST(Vtr, RoundTripsAUniformLatticeExactly) {
    for (bool binary : {false, true}) {
        Mesh m = vtr_lattice(3, -0.5, 0.25);
        const std::string path = mt::temp_path(".vtr");
        write_vtr(path, m, binary, /*zlib=*/false);
        const Mesh back = read_vtr(path);

        ASSERT_EQ(back.NumPoints(), m.NumPoints());
        for (std::size_t p = 0; p < m.NumPoints(); ++p)
            for (std::size_t c = 0; c < 3; ++c)
                EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.Points(), p * 3 + c),
                                 meshioplusplus::detail::read_double(m.Points(), p * 3 + c));
        ASSERT_EQ(back.NumCellBlocks(), 1u);
        EXPECT_EQ(back.Cells(0).Type(), "hexahedron");
        EXPECT_EQ(back.Cells(0).NumCells(), 27u);
        std::remove(path.c_str());
    }
}

TEST(Vtr, RefusesAMeshThatIsNotAUniformLattice) {
    Mesh m = mt::tet_mesh();
    const std::string path = mt::temp_path(".vtr");
    EXPECT_THROW(write_vtr(path, m, false, false), WriteError);
}

TEST(Vtr, MetadataAgreesWithARealRead) {
    Mesh m = vtr_lattice(2, 0.0, 1.0);
    const std::string path = mt::temp_path(".vtr");
    write_vtr(path, m, true, false);

    const meshioplusplus::MeshMetadata meta = read_vtr_metadata(path);
    const Mesh back = read_vtr(path);
    EXPECT_EQ(meta.mNumPoints, back.NumPoints());
    ASSERT_EQ(meta.mCellBlocks.size(), 1u);
    EXPECT_EQ(meta.mCellBlocks[0].mType, "hexahedron");
    EXPECT_EQ(meta.mCellBlocks[0].mNumCells, back.Cells(0).NumCells());
    std::remove(path.c_str());
}

// The identity .vtr exists for, over .vti: genuinely non-uniform per-axis
// spacing -- a graded grid, finer near x=0. The reader never checks
// uniformity at all, unlike the writer.
TEST(Vtr, ReadsAGenuinelyGradedGrid) {
    const std::string path = mt::temp_path(".vtr");
    std::ofstream os(path, std::ios::binary);
    os << "<?xml version=\"1.0\"?>\n"
          "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
          "<RectilinearGrid WholeExtent=\"0 2 0 1 0 1\">\n"
          "<Piece Extent=\"0 2 0 1 0 1\">\n"
          "<Coordinates>\n"
          "<DataArray type=\"Float64\" Name=\"x_coordinates\" format=\"ascii\">\n"
          "0 1 4\n</DataArray>\n"
          "<DataArray type=\"Float64\" Name=\"y_coordinates\" format=\"ascii\">\n"
          "0 10\n</DataArray>\n"
          "<DataArray type=\"Float64\" Name=\"z_coordinates\" format=\"ascii\">\n"
          "0 100\n</DataArray>\n"
          "</Coordinates>\n"
          "</Piece>\n"
          "</RectilinearGrid>\n"
          "</VTKFile>\n";
    os.close();

    const Mesh back = read_vtr(path);
    EXPECT_EQ(back.NumPoints(), 12u);  // 3 x 2 x 2
    ASSERT_EQ(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).NumCells(), 2u);  // 2 x 1 x 1
    // Point (i=2, j=0, k=0): x-fastest index 2.
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.Points(), 2 * 3 + 0), 4.0);
    // Point (i=0, j=1, k=0): x-fastest index 3.
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.Points(), 3 * 3 + 1), 10.0);
    // Point (i=0, j=0, k=1): x-fastest index 6.
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.Points(), 6 * 3 + 2), 100.0);
    std::remove(path.c_str());
}
