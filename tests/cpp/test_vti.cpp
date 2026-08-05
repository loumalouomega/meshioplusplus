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
// VTK XML ImageData.
//
// A round trip through our own writer is a WEAK oracle here and deliberately not
// the only one: the reader and the writer share a coordinate convention, so a
// consistently wrong `Origin` or a transposed extent survives one unchanged.
// The tests below therefore also (a) read hand-written files whose expected
// geometry is stated independently -- a non-zero Origin, an extent not starting
// at zero, an inverted extent, a foreign Piece -- and (b) assert on the raw
// written bytes.

// System includes
#include <array>
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
#include "meshioplusplus/detail/grid_lattice.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vti.hpp"
#include "meshioplusplus/operations/voxelize.hpp"

using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::read_vti;
using meshioplusplus::read_vti_metadata;
using meshioplusplus::ReadError;
using meshioplusplus::write_vti;
using meshioplusplus::WriteError;
namespace d = meshioplusplus::detail;

namespace {

std::string vti_write_temp(const std::string& rText) {
    const std::string path = mt::temp_path(".vti");
    std::ofstream os(path, std::ios::binary);
    os << rText;
    os.close();
    return path;
}

std::string vti_slurp(const std::string& rPath) {
    std::ifstream in(rPath, std::ios::binary);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

Mesh vti_lattice(std::int64_t N, double Origin, double Spacing) {
    return meshioplusplus::grid({{N, N, N}}, {{Origin, Origin, Origin}},
                                {{Spacing, Spacing, Spacing}});
}

}  // namespace

TEST(Vti, RoundTripsALatticeExactly) {
    for (bool binary : {false, true}) {
        Mesh m = vti_lattice(3, -0.5, 0.25);
        NDArray pd = NDArray::Uninit(meshioplusplus::DType::Float64, {m.NumPoints()});
        for (std::size_t p = 0; p < m.NumPoints(); ++p)
            pd.As<double>()[p] = static_cast<double>(p) * 0.5;
        m.AddPointData("f", std::move(pd));
        NDArray cd = NDArray::Uninit(meshioplusplus::DType::Int64, {27u});
        for (std::size_t c = 0; c < 27; ++c)
            cd.As<std::int64_t>()[c] = static_cast<std::int64_t>(c);
        std::vector<NDArray> blocks;
        blocks.push_back(std::move(cd));
        m.AddCellData("tag", std::move(blocks));

        const std::string path = mt::temp_path(".vti");
        write_vti(path, m, binary, /*zlib=*/false);
        const Mesh back = read_vti(path);

        ASSERT_EQ(back.NumPoints(), m.NumPoints());
        ASSERT_EQ(back.NumCellBlocks(), 1u);
        EXPECT_EQ(back.Cells(0).NumCells(), 27u);
        for (std::size_t p = 0; p < m.NumPoints(); ++p)
            for (std::size_t k = 0; k < 3; ++k)
                EXPECT_DOUBLE_EQ(d::read_double(back.Points(), p * 3 + k),
                                 d::read_double(m.Points(), p * 3 + k));
        ASSERT_TRUE(back.HasPointData("f"));
        for (std::size_t p = 0; p < m.NumPoints(); ++p)
            EXPECT_DOUBLE_EQ(d::read_double(back.PointData("f"), p), static_cast<double>(p) * 0.5);
        ASSERT_TRUE(back.HasCellData("tag"));
        for (std::size_t c = 0; c < 27; ++c)
            EXPECT_EQ(d::read_int(back.CellData("tag", 0), c), static_cast<std::int64_t>(c));
        std::remove(path.c_str());
    }
}

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
TEST(Vti, RoundTripsCompressed) {
    Mesh m = vti_lattice(4, 0.0, 1.0);
    NDArray pd = NDArray::Uninit(meshioplusplus::DType::Float64, {m.NumPoints()});
    for (std::size_t p = 0; p < m.NumPoints(); ++p)
        pd.As<double>()[p] = std::sin(static_cast<double>(p));
    m.AddPointData("f", std::move(pd));
    const std::string path = mt::temp_path(".vti");
    write_vti(path, m, /*binary=*/true, /*zlib=*/true);
    EXPECT_NE(vti_slurp(path).find("vtkZLibDataCompressor"), std::string::npos);
    const Mesh back = read_vti(path);
    for (std::size_t p = 0; p < m.NumPoints(); ++p)
        EXPECT_DOUBLE_EQ(d::read_double(back.PointData("f"), p),
                         d::read_double(m.PointData("f"), p));
    std::remove(path.c_str());
}
#endif

// The written bytes, not a read-back: a consistently wrong Origin/Spacing pair
// survives a round trip through this file's own reader unchanged.
TEST(Vti, TheWrittenAttributesAreTheGeometry) {
    const std::string path = mt::temp_path(".vti");
    write_vti(path, vti_lattice(3, -0.5, 0.25), /*binary=*/false, /*zlib=*/false);
    const std::string text = vti_slurp(path);
    std::remove(path.c_str());
    EXPECT_NE(text.find("type=\"ImageData\""), std::string::npos);
    EXPECT_NE(text.find("WholeExtent=\"0 3 0 3 0 3\""), std::string::npos);
    EXPECT_NE(text.find("Origin=\"-0.5 -0.5 -0.5\""), std::string::npos);
    EXPECT_NE(text.find("Spacing=\"0.25 0.25 0.25\""), std::string::npos);
    // No point array at all: the geometry IS the three attributes, and writing
    // points as well would make the file self-contradictory.
    EXPECT_EQ(text.find("<Points>"), std::string::npos);
}

// A hand-written file whose expected geometry is stated here rather than
// produced by our writer.
TEST(Vti, ReadsAForeignFileWithANonZeroOriginAndAShiftedExtent) {
    const std::string path = vti_write_temp(
        "<?xml version=\"1.0\"?>\n"
        "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        "<ImageData WholeExtent=\"2 4 0 1 0 1\" Origin=\"1 2 3\" Spacing=\"0.5 2 4\">\n"
        "<Piece Extent=\"2 4 0 1 0 1\">\n"
        "<PointData><DataArray type=\"Int32\" Name=\"i\" format=\"ascii\">\n"
        "0 1 2 3 4 5 6 7 8 9 10 11\n</DataArray></PointData>\n"
        "</Piece></ImageData></VTKFile>\n");
    const Mesh m = read_vti(path);
    std::remove(path.c_str());

    // 2 x 1 x 1 cells, so 3 x 2 x 2 = 12 points.
    ASSERT_EQ(m.NumPoints(), 12u);
    ASSERT_EQ(m.NumCellBlocks(), 1u);
    EXPECT_EQ(m.Cells(0).NumCells(), 2u);
    // The extent starts at (2, 0, 0), so the lo corner is Origin + start*Spacing
    // = (1 + 2*0.5, 2, 3) = (2, 2, 3). Dropping that offset would translate the
    // whole grid by one cell in x -- silently, and plausibly.
    EXPECT_DOUBLE_EQ(d::read_double(m.Points(), 0), 2.0);
    EXPECT_DOUBLE_EQ(d::read_double(m.Points(), 1), 2.0);
    EXPECT_DOUBLE_EQ(d::read_double(m.Points(), 2), 3.0);
    // x fastest: point 1 is one spacing along x.
    EXPECT_DOUBLE_EQ(d::read_double(m.Points(), 3), 2.5);
    // point 3 is the first of the next y row.
    EXPECT_DOUBLE_EQ(d::read_double(m.Points(), 3 * 3 + 1), 4.0);
    ASSERT_TRUE(m.HasPointData("i"));
    EXPECT_EQ(d::read_int(m.PointData("i"), 11), 11);
}

TEST(Vti, MetadataAgreesWithTheRealRead) {
    const std::string path = mt::temp_path(".vti");
    Mesh m = vti_lattice(3, -0.5, 0.25);
    NDArray pd = NDArray(meshioplusplus::DType::Float64, {m.NumPoints()});
    m.AddPointData("f", std::move(pd));
    write_vti(path, m, /*binary=*/false, /*zlib=*/false);
    const meshioplusplus::MeshMetadata meta = read_vti_metadata(path);
    const Mesh back = read_vti(path);
    std::remove(path.c_str());

    EXPECT_EQ(meta.mNumPoints, back.NumPoints());
    ASSERT_EQ(meta.mCellBlocks.size(), 1u);
    EXPECT_EQ(meta.mCellBlocks[0].mType, "hexahedron");
    EXPECT_EQ(meta.mCellBlocks[0].mNumCells, back.Cells(0).NumCells());
    ASSERT_EQ(meta.mPointDataNames.size(), 1u);
    EXPECT_EQ(meta.mPointDataNames[0], "f");
    // Unlike every other native metadata path this one CAN report the bounding
    // box, because the extent is the bounding box.
    EXPECT_TRUE(meta.mHasBBox);
    EXPECT_DOUBLE_EQ(meta.mBBoxMin[0], -0.5);
    EXPECT_DOUBLE_EQ(meta.mBBoxMax[0], 0.25);
}

TEST(Vti, RefusesAMeshThatIsNotALattice) {
    const std::string path = mt::temp_path(".vti");
    // A tetra mesh has no extent at all.
    EXPECT_THROW(write_vti(path, mt::tet_mesh(), false, false), WriteError);
    // A PARTIAL lattice is the interesting refusal: it is made of the right
    // cells in the right places, and ImageData still cannot express its holes.
    meshioplusplus::VoxelOptions vo;
    vo.mResolution = std::array<std::int64_t, 3>{{4, 4, 4}};
    vo.mFill = meshioplusplus::VoxelFill::Surface;
    const Mesh partial = meshioplusplus::voxelize(mt::tri_mesh(), vo).mMesh;
    EXPECT_THROW(write_vti(path, partial, false, false), WriteError);
    std::remove(path.c_str());
}

TEST(Vti, DeclinesTheConstructsItDoesNotImplement) {
    const char* cases[] = {
        // appended data
        "<VTKFile type=\"ImageData\"><AppendedData encoding=\"raw\">_</AppendedData>"
        "<ImageData WholeExtent=\"0 1 0 1 0 1\"><Piece Extent=\"0 1 0 1 0 1\"/>"
        "</ImageData></VTKFile>",
        // two pieces
        "<VTKFile type=\"ImageData\"><ImageData WholeExtent=\"0 1 0 1 0 1\">"
        "<Piece Extent=\"0 1 0 1 0 1\"/><Piece Extent=\"0 1 0 1 0 1\"/>"
        "</ImageData></VTKFile>",
        // a piece that is not the whole extent
        "<VTKFile type=\"ImageData\"><ImageData WholeExtent=\"0 2 0 1 0 1\">"
        "<Piece Extent=\"0 1 0 1 0 1\"/></ImageData></VTKFile>",
        // a rotated lattice
        "<VTKFile type=\"ImageData\"><ImageData WholeExtent=\"0 1 0 1 0 1\" "
        "Direction=\"0 1 0 -1 0 0 0 0 1\"><Piece Extent=\"0 1 0 1 0 1\"/>"
        "</ImageData></VTKFile>",
        // an inverted extent
        "<VTKFile type=\"ImageData\"><ImageData WholeExtent=\"2 0 0 1 0 1\">"
        "<Piece Extent=\"2 0 0 1 0 1\"/></ImageData></VTKFile>",
        // the wrong dataset type
        "<VTKFile type=\"UnstructuredGrid\"><UnstructuredGrid/></VTKFile>",
    };
    for (const char* text : cases) {
        const std::string path = vti_write_temp(text);
        EXPECT_THROW(read_vti(path), ReadError) << text;
        std::remove(path.c_str());
    }
}

TEST(Vti, AnArrayOfTheWrongLengthIsAnErrorRatherThanATruncation) {
    const std::string path = vti_write_temp(
        "<VTKFile type=\"ImageData\"><ImageData WholeExtent=\"0 1 0 1 0 1\">"
        "<Piece Extent=\"0 1 0 1 0 1\">"
        "<PointData><DataArray type=\"Int32\" Name=\"i\" format=\"ascii\">1 2 3</DataArray>"
        "</PointData></Piece></ImageData></VTKFile>");
    EXPECT_THROW(read_vti(path), ReadError);
    std::remove(path.c_str());
}

// The identity the format exists for: written and read back, a grid IS the same
// grid -- geometry included, which no other format in meshio++ can promise for a
// mesh whose only description was `field_data`.
TEST(Vti, TheRecoveredHeaderEqualsTheWrittenOne) {
    const Mesh m = vti_lattice(5, 1.25, 0.125);
    d::LatticeSpec before;
    ASSERT_TRUE(d::lattice_from_mesh(m, before));

    const std::string path = mt::temp_path(".vti");
    write_vti(path, m, /*binary=*/true, /*zlib=*/false);
    const Mesh back = read_vti(path);
    std::remove(path.c_str());

    d::LatticeSpec after;
    ASSERT_TRUE(d::lattice_from_mesh(back, after));
    for (std::size_t k = 0; k < 3; ++k) {
        EXPECT_EQ(after.mDims[k], before.mDims[k]);
        EXPECT_DOUBLE_EQ(after.mOrigin[k], before.mOrigin[k]);
        EXPECT_DOUBLE_EQ(after.mSpacing[k], before.mSpacing[k]);
    }
}

TEST(LatticeFromMesh, RejectsWhatIsNotADenseLattice) {
    d::LatticeSpec spec;
    // A tetra mesh.
    EXPECT_FALSE(d::lattice_from_mesh(mt::tet_mesh(), spec));
    // A lattice whose points have been permuted: identical plane sets, different
    // mesh. Checking only the plane positions would accept this.
    Mesh m = vti_lattice(2, 0.0, 1.0);
    NDArray pts = NDArray::Uninit(meshioplusplus::DType::Float64, {m.NumPoints(), 3u});
    for (std::size_t p = 0; p < m.NumPoints(); ++p) {
        const std::size_t src = m.NumPoints() - 1 - p;
        for (std::size_t k = 0; k < 3; ++k)
            pts.As<double>()[p * 3 + k] = d::read_double(m.Points(), src * 3 + k);
    }
    Mesh permuted;
    permuted.AssignPoints(std::move(pts));
    permuted.AddCellBlock("hexahedron", m.Cells(0).Conn());
    EXPECT_FALSE(d::lattice_from_mesh(permuted, spec));
    // A graded grid: uniform plane counts, non-uniform spacing.
    Mesh graded;
    std::vector<double> xs = {0.0, 1.0, 3.0};
    NDArray gp = NDArray::Uninit(meshioplusplus::DType::Float64, {27u, 3u});
    std::size_t idx = 0;
    for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t j = 0; j < 3; ++j)
            for (std::size_t i = 0; i < 3; ++i, ++idx) {
                gp.As<double>()[idx * 3 + 0] = xs[i];
                gp.As<double>()[idx * 3 + 1] = static_cast<double>(j);
                gp.As<double>()[idx * 3 + 2] = static_cast<double>(k);
            }
    graded.AssignPoints(std::move(gp));
    graded.AddCellBlock("hexahedron", vti_lattice(2, 0.0, 1.0).Cells(0).Conn());
    EXPECT_FALSE(d::lattice_from_mesh(graded, spec));
}
