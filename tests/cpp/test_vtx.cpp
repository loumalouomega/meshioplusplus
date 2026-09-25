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

// DOLFINx VTX (.bp) reader: synthetic files written through ADIOS2's C++ API in
// VTXWriter's layout (the real DOLFINx runs are the Python tests' fixtures).
#ifdef MESHIOPLUSPLUS_HAS_ADIOS2

// System includes
#include <array>
#include <atomic>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

// External includes
#include <adios2.h>
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/vtx.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::GhostPolicy;
using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;

class VtxTempDir {
public:
    VtxTempDir() {
        static std::atomic<unsigned> counter{0};
        mPath = std::filesystem::temp_directory_path() /
                ("meshio_vtx_" + std::to_string(counter++) + "_" +
                 std::to_string(::testing::UnitTest::GetInstance()->random_seed()));
        std::filesystem::create_directories(mPath);
    }
    ~VtxTempDir() {
        std::error_code ec;
        std::filesystem::remove_all(mPath, ec);
    }
    VtxTempDir(const VtxTempDir&) = delete;
    VtxTempDir& operator=(const VtxTempDir&) = delete;
    std::string operator/(const std::string& rName) const { return (mPath / rName).string(); }
    const std::filesystem::path& Path() const { return mPath; }

private:
    std::filesystem::path mPath;
};

const char* const kSchema = R"(<VTKFile type="UnstructuredGrid" version="0.1">
  <UnstructuredGrid>
    <Piece NumberOfPoints="NumberOfNodes" NumberOfCells="NumberOfCells">
      <Points><DataArray Name="geometry" /></Points>
      <Cells>
        <DataArray Name="connectivity" />
        <DataArray Name="types" />
      </Cells>
      <PointData>
        <DataArray Name="TIME">step</DataArray>
        <DataArray Name="vtkOriginalPointIds" />
        <DataArray Name="vtkGhostType" />
        <DataArray Name="u" />
      </PointData>
      <CellData><DataArray Name="rank" /></CellData>
    </Piece>
  </UnstructuredGrid>
</VTKFile>)";

// One "rank": a local block. Two unit tetrahedra of degree 2 split over two
// ranks that share the face (1, 2, 3): each rank holds its tetrahedron's ten
// nodes, the shared face's six as ghosts on rank 1.
struct Rank {
    std::vector<double> mX;            // n x 3
    std::vector<std::int64_t> mConn;   // 1 x 11, count-prefixed
    std::vector<std::int64_t> mIds;    // global node ids
    std::vector<std::uint8_t> mGhost;  // 1 = ghost
};

// The ten P2 nodes of the tetrahedron (a, b, c, d) in VTK order: corners, then
// the edge midpoints (ab, bc, ca, ad, bd, cd).
std::vector<std::array<double, 3>> vtx_p2_nodes(const std::array<std::array<double, 3>, 4>& rV) {
    auto mid = [&](int i, int j) {
        return std::array<double, 3>{(rV[i][0] + rV[j][0]) / 2, (rV[i][1] + rV[j][1]) / 2,
                                     (rV[i][2] + rV[j][2]) / 2};
    };
    return {rV[0],     rV[1],     rV[2],     rV[3],     mid(0, 1),
            mid(1, 2), mid(2, 0), mid(0, 3), mid(1, 3), mid(2, 3)};
}

std::vector<Rank> vtx_two_ranks() {
    const std::array<double, 3> p0{0, 0, 0}, p1{1, 0, 0}, p2{0, 1, 0}, p3{0, 0, 1}, p4{1, 1, 1};
    const auto a = vtx_p2_nodes({p0, p1, p2, p3});
    const auto b = vtx_p2_nodes({p4, p1, p2, p3});
    // Global ids: the first tetrahedron's nodes 0..9; the second's corner p4 is
    // 10, its edges (p4,p1) 11, (p1,p2) 5, (p2,p4) 12, (p4,p3) 13, (p1,p3) 8,
    // (p2,p3) 9 -- node k of the second tetrahedron in VTK order:
    const std::vector<std::int64_t> ids_b = {10, 1, 2, 3, 11, 5, 12, 13, 8, 9};
    const std::vector<std::uint8_t> ghost_b = {0, 1, 1, 1, 0, 1, 0, 0, 1, 1};
    Rank r0, r1;
    for (int k = 0; k < 10; ++k) {
        r0.mX.insert(r0.mX.end(), a[k].begin(), a[k].end());
        r1.mX.insert(r1.mX.end(), b[k].begin(), b[k].end());
        r0.mIds.push_back(k);
        r0.mGhost.push_back(0);
    }
    r1.mIds = ids_b;
    r1.mGhost = ghost_b;
    r0.mConn.push_back(10);
    r1.mConn.push_back(10);
    for (int k = 0; k < 10; ++k) {
        r0.mConn.push_back(k);
        r1.mConn.push_back(k);
    }
    return {r0, r1};
}

// Writes a VTX file: the mesh in step 0 only (VTXMeshPolicy.reuse), `u` (the
// global id plus 100 * step) and the ghost arrays every step. With
// `StrayBlock`, step 0's `vtkGhostType` gets a trailing one-value block.
void vtx_write(const std::string& rPath, int Steps, bool Schema = true, bool StrayBlock = false) {
    adios2::ADIOS adios;
    adios2::IO io = adios.DeclareIO("vtx_test_writer");
    io.SetEngine("BP5");
    adios2::Engine engine = io.Open(rPath, adios2::Mode::Write);
    if (Schema)
        io.DefineAttribute<std::string>("vtk.xml", kSchema);
    const std::vector<Rank> ranks = vtx_two_ranks();
    auto step_var = io.DefineVariable<double>("step");
    auto geometry = io.DefineVariable<double>("geometry", {}, {}, {10, 3});
    auto conn = io.DefineVariable<std::int64_t>("connectivity", {}, {}, {1, 11});
    auto types = io.DefineVariable<std::uint32_t>("types");
    auto ids = io.DefineVariable<std::int64_t>("vtkOriginalPointIds", {}, {}, {10});
    auto ghost = io.DefineVariable<std::uint8_t>("vtkGhostType", {}, {}, {10});
    auto u = io.DefineVariable<double>("u", {}, {}, {10});
    auto rank_var = io.DefineVariable<double>("rank", {}, {}, {1});
    for (int s = 0; s < Steps; ++s) {
        engine.BeginStep();
        engine.Put(step_var, 0.5 * s);
        if (s == 0)
            engine.Put(types, std::uint32_t{71});
        for (std::size_t r = 0; r < ranks.size(); ++r) {
            const Rank& rR = ranks[r];
            if (s == 0) {
                engine.Put(geometry, rR.mX.data(), adios2::Mode::Sync);
                engine.Put(conn, rR.mConn.data(), adios2::Mode::Sync);
                const double rv = static_cast<double>(r);
                engine.Put(rank_var, &rv, adios2::Mode::Sync);
            }
            engine.Put(ids, rR.mIds.data(), adios2::Mode::Sync);
            engine.Put(ghost, rR.mGhost.data(), adios2::Mode::Sync);
            std::vector<double> values;
            for (std::int64_t id : rR.mIds)
                values.push_back(static_cast<double>(id) + 100.0 * s);
            engine.Put(u, values.data(), adios2::Mode::Sync);
        }
        if (StrayBlock && s == 0) {
            ghost.SetSelection({{}, {1}});
            const std::uint8_t one = 1;
            engine.Put(ghost, &one, adios2::Mode::Sync);
            ghost.SetSelection({{}, {10}});
        }
        engine.EndStep();
    }
    engine.Close();
}

TEST(Vtx, ConcatenatesRankBlocksAndLowersLagrangeCells) {
    VtxTempDir dir;
    const std::string path = dir / "run.bp";
    vtx_write(path, 3);
    const Mesh mesh = meshioplusplus::read_vtx(path);
    ASSERT_EQ(mesh.NumPoints(), 20u);
    ASSERT_EQ(mesh.NumCellBlocks(), 1u);
    const auto cells = mesh.Cells(0);
    EXPECT_EQ(cells.Type(), "tetra10");
    ASSERT_EQ(cells.NumCells(), 2u);
    const auto* c = cells.Conn().As<std::int64_t>();
    EXPECT_EQ(c[0], 0);
    EXPECT_EQ(c[10], 10);  // the second rank's numbering is offset
    EXPECT_EQ(c[19], 19);
    const auto* rank = mesh.CellData("rank", 0).As<double>();
    EXPECT_EQ(rank[0], 0.0);
    EXPECT_EQ(rank[1], 1.0);
    EXPECT_EQ(mesh.FieldData("meshio:time").As<double>()[0], 0.0);
}

TEST(Vtx, ReusePolicyCarriesTheMeshForward) {
    VtxTempDir dir;
    const std::string path = dir / "run.bp";
    vtx_write(path, 3);
    ReadOptions opts;
    opts.mTimeStep = 2;
    const Mesh mesh = meshioplusplus::read_vtx(path, opts);
    EXPECT_EQ(mesh.NumPoints(), 20u);
    EXPECT_EQ(mesh.FieldData("meshio:time").As<double>()[0], 1.0);
    const auto* u = mesh.PointData("u").As<double>();
    EXPECT_EQ(u[0], 200.0);
    EXPECT_EQ(u[10], 210.0);
    EXPECT_EQ(meshioplusplus::vtx_time_values(path), (std::vector<double>{0.0, 0.5, 1.0}));
    opts.mTimeStep = 3;
    EXPECT_THROW(meshioplusplus::read_vtx(path, opts), ReadError);
}

TEST(Vtx, DropGhostsWeldsOntoOwners) {
    VtxTempDir dir;
    const std::string path = dir / "run.bp";
    vtx_write(path, 1);
    ReadOptions opts;
    opts.mGhosts = GhostPolicy::Drop;
    const Mesh mesh = meshioplusplus::read_vtx(path, opts);
    EXPECT_EQ(mesh.NumPoints(), 14u);  // 20 minus the shared face's six ghosts
    EXPECT_FALSE(mesh.HasPointData("vtkGhostType"));
    EXPECT_FALSE(mesh.HasPointData("vtkOriginalPointIds"));
    const auto* c = mesh.Cells(0).Conn().As<std::int64_t>();
    const auto* u = mesh.PointData("u").As<double>();
    // The second tetrahedron's corner 1 is the first one's corner 1 (id 1).
    EXPECT_EQ(c[11], c[1]);
    EXPECT_EQ(u[c[10]], 10.0);
}

TEST(Vtx, SkipsBlocksMatchingNoRank) {
    VtxTempDir dir;
    const std::string path = dir / "run.bp";
    vtx_write(path, 2, true, true);
    const Mesh mesh = meshioplusplus::read_vtx(path);
    ASSERT_TRUE(mesh.HasPointData("vtkGhostType"));
    EXPECT_EQ(mesh.PointData("vtkGhostType").Shape()[0], 20u);
}

TEST(Vtx, RefusesAFileWithoutSchema) {
    VtxTempDir dir;
    const std::string path = dir / "plain.bp";
    vtx_write(path, 1, false);
    try {
        meshioplusplus::read_vtx(path);
        FAIL() << "expected a ReadError";
    } catch (const ReadError& rE) {
        EXPECT_NE(std::string(rE.what()).find("adios4dolfinx"), std::string::npos);
    }
    EXPECT_THROW(meshioplusplus::read_vtx(dir / "missing.bp"), ReadError);
}

TEST(Vtx, MetadataFromTheBlockLists) {
    VtxTempDir dir;
    const std::string path = dir / "run.bp";
    vtx_write(path, 3);
    const auto meta = meshioplusplus::read_vtx_metadata(path, {});
    EXPECT_FALSE(meta.mFellBackToFullRead);
    EXPECT_EQ(meta.mNumPoints, 20u);
    ASSERT_EQ(meta.mCellBlocks.size(), 1u);
    EXPECT_EQ(meta.mCellBlocks[0].mType, "tetra10");
    EXPECT_EQ(meta.mCellBlocks[0].mNumCells, 2u);
    EXPECT_EQ(meta.mTimeValues, (std::vector<double>{0.0, 0.5, 1.0}));
    EXPECT_EQ(meta.mPointDataNames,
              (std::vector<std::string>{"u", "vtkGhostType", "vtkOriginalPointIds"}));
    EXPECT_EQ(meta.mCellDataNames, (std::vector<std::string>{"rank"}));
}

TEST(Vtx, RegistrySniffAndSequenceGlob) {
    VtxTempDir dir;
    vtx_write(dir / "run_0.bp", 1);
    vtx_write(dir / "run_1.bp", 2);
    EXPECT_EQ(meshioplusplus::resolve_format(dir / "run_0.bp", ""), "vtx");
    std::filesystem::rename(dir / "run_1.bp", dir / "noext");
    EXPECT_EQ(meshioplusplus::sniff_format(dir / "noext"), "vtx");
    std::filesystem::rename(dir / "noext", dir / "run_1.bp");
    EXPECT_TRUE(meshioplusplus::seq_format_may_have_steps("vtx"));
    meshioplusplus::SequenceInput in;
    in.mPattern = dir / "run_*.bp";
    const auto entries = meshioplusplus::sequence_expand(in);
    ASSERT_EQ(entries.size(), 3u);  // one step, then two
    EXPECT_EQ(std::filesystem::path(entries[0].mPath).filename(), "run_0.bp");
    EXPECT_EQ(std::filesystem::path(entries[2].mPath).filename(), "run_1.bp");
}

}  // namespace

#endif  // MESHIOPLUSPLUS_HAS_ADIOS2
