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
// VTK XML parallel indices `.pvtu` / `.pvtp` (v15.0.0): an index
// that declares the arrays and names one `.vtu`/`.vtp` piece per part.

// System includes
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
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
#include "meshioplusplus/formats/pvtp.hpp"
#include "meshioplusplus/formats/pvtu.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/registry.hpp"
#include "meshioplusplus/region.hpp"

using meshioplusplus::DType;
using meshioplusplus::GhostPolicy;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::read_pvtp;
using meshioplusplus::read_pvtp_metadata;
using meshioplusplus::read_pvtu;
using meshioplusplus::read_pvtu_metadata;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;
using meshioplusplus::RegionKind;
using meshioplusplus::write_pvtp;
using meshioplusplus::write_pvtp_pieces_codec;
using meshioplusplus::write_pvtu;
using meshioplusplus::write_pvtu_codec;
using meshioplusplus::write_pvtu_pieces_codec;
using meshioplusplus::WriteError;
using meshioplusplus::detail::VtkCodec;

namespace {

namespace fs = std::filesystem;

// Remove the index file and its sibling pieces directory (`<stem>/`).
void pvtu_cleanup(const std::string& rPath) {
    std::error_code ec;
    const fs::path p(rPath);
    fs::remove(p, ec);
    fs::remove_all(p.parent_path() / p.stem(), ec);
}

std::string pvtu_slurp(const std::string& rPath) {
    std::ifstream in(rPath);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

std::size_t pvtu_count(const std::string& rText, const std::string& rNeedle) {
    std::size_t n = 0;
    for (std::size_t at = rText.find(rNeedle); at != std::string::npos;
         at = rText.find(rNeedle, at + rNeedle.size()))
        ++n;
    return n;
}

std::size_t pvtu_num_cells(const Mesh& rMesh) {
    std::size_t n = 0;
    for (const auto cb : rMesh.CellRange())
        n += cb.NumCells();
    return n;
}

// An n x n square split into triangles, with a point array `u` and a cell array `c`.
Mesh pvtu_grid(int N) {
    const std::size_t np = static_cast<std::size_t>((N + 1) * (N + 1));
    const std::size_t nc = static_cast<std::size_t>(2 * N * N);
    NDArray pts(DType::Float64, {np, 3});
    NDArray u(DType::Float64, {np});
    double* p = pts.As<double>();
    for (int j = 0; j <= N; ++j)
        for (int i = 0; i <= N; ++i) {
            const std::size_t k = static_cast<std::size_t>(j * (N + 1) + i);
            p[3 * k] = i;
            p[3 * k + 1] = j;
            p[3 * k + 2] = 0.0;
            u.As<double>()[k] = i + 10.0 * j;
        }
    NDArray conn(DType::Int64, {nc, 3});
    std::int64_t* c = conn.As<std::int64_t>();
    std::size_t at = 0;  // two triangles per square: (a, b, e) and (a, e, d)
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
            const std::int64_t a = j * (N + 1) + i, b = a + 1, d = a + N + 1, e = d + 1;
            const std::int64_t tris[6] = {a, b, e, a, e, d};
            std::copy(tris, tris + 6, c + 3 * at);
            at += 2;
        }
    NDArray cd(DType::Float64, {nc});
    for (std::size_t i = 0; i < nc; ++i)
        cd.As<double>()[i] = static_cast<double>(i);

    Mesh m;
    m.AssignPoints(std::move(pts));
    m.AddCellBlock("triangle", std::move(conn));
    m.AddPointData("u", std::move(u));
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(cd));
    m.AddCellData("c", std::move(blocks));
    return m;
}

// `pvtu_grid` with `partition:part` attached, as `partition --labels-only` does.
Mesh pvtu_labelled(int N, int NParts) {
    Mesh m = pvtu_grid(N);
    meshioplusplus::PartitionOptions o;
    o.mNParts = NParts;
    o.mMethod = meshioplusplus::PartitionMethod::SFC;
    m.AddCellData("partition:part", meshioplusplus::partition_labels(m, o));
    return m;
}

std::vector<std::array<double, 3>> pvtu_sorted_points(const Mesh& rMesh) {
    const NDArray& pts = rMesh.Points();
    const double* p = pts.As<double>();
    std::vector<std::array<double, 3>> out;
    for (std::size_t i = 0; i < rMesh.NumPoints(); ++i)
        out.push_back({p[3 * i], p[3 * i + 1], p[3 * i + 2]});
    std::sort(out.begin(), out.end());
    return out;
}

Mesh pvtu_weld(const Mesh& rMesh) {
    meshioplusplus::CleanOptions co;
    co.weld = true;
    co.remove_orphans = false;
    co.drop_degenerate = false;
    co.drop_duplicate_cells = false;
    return meshioplusplus::clean(rMesh, co).mMesh;
}

// A one-triangle piece with point array `pName` of dtype `Dt`, for declaration tests.
Mesh pvtu_piece(DType Dt = DType::Float64, const char* pName = "u", std::size_t Comps = 1) {
    NDArray pts(DType::Float64, {3, 3});
    double* p = pts.As<double>();
    p[3] = 1.0;
    p[7] = 1.0;
    NDArray conn(DType::Int64, {1, 3});
    conn.As<std::int64_t>()[1] = 1;
    conn.As<std::int64_t>()[2] = 2;
    Mesh m;
    m.AssignPoints(std::move(pts));
    m.AddCellBlock("triangle", std::move(conn));
    m.AddPointData(pName, Comps == 1 ? NDArray(Dt, {3}) : NDArray(Dt, {3, Comps}));
    std::vector<NDArray> blocks;
    blocks.push_back(NDArray(DType::Float64, {1}));
    m.AddCellData("c", std::move(blocks));
    return m;
}

std::vector<const Mesh*> pvtu_ptrs(const std::vector<Mesh>& rMeshes) {
    std::vector<const Mesh*> out;
    for (const Mesh& m : rMeshes)
        out.push_back(&m);
    return out;
}

const VtkCodec kNone = VtkCodec::None;

}  // namespace

TEST(Pvtu, RoundTripsThreePartsAsThreeNamedCellRegions) {
    const Mesh m = pvtu_labelled(6, 3);
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu(path, m, /*binary=*/false, /*zlib=*/false);
    const Mesh back = read_pvtu(path);

    EXPECT_EQ(pvtu_num_cells(back), pvtu_num_cells(m));
    EXPECT_GT(back.NumPoints(), m.NumPoints());  // interface points once per piece
    ASSERT_EQ(back.NumRegions(), 3u);
    std::size_t entries = 0;
    for (std::size_t i = 0; i < 3; ++i) {
        EXPECT_EQ(back.Region(i).mKind, RegionKind::Cell);
        EXPECT_EQ(back.Region(i).mName, "piece_" + std::to_string(i));
        entries += back.Region(i).NumEntries();
    }
    EXPECT_EQ(entries, pvtu_num_cells(m));
    EXPECT_TRUE(back.HasPointData("u"));
    EXPECT_TRUE(back.HasCellData("c"));
    EXPECT_TRUE(back.HasCellData("partition:part"));
    EXPECT_EQ(back.NumFieldData(), 0u);

    pvtu_cleanup(path);
}

TEST(Pvtu, WritesTheDocumentedIndexAndPieces) {
    const Mesh m = pvtu_labelled(6, 3);
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu(path, m, /*binary=*/true, /*zlib=*/false);

    const std::string text = pvtu_slurp(path);
    EXPECT_NE(text.find("<VTKFile type=\"PUnstructuredGrid\""), std::string::npos);
    EXPECT_NE(text.find("<PUnstructuredGrid GhostLevel=\"0\">"), std::string::npos);
    EXPECT_NE(text.find("<PDataArray type=\"Float64\" Name=\"u\"/>"), std::string::npos);
    EXPECT_NE(text.find("<PDataArray type=\"Int64\" Name=\"partition:part\"/>"), std::string::npos);
    EXPECT_NE(text.find("<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\"/>"),
              std::string::npos);
    const fs::path idx(path);
    const std::string stem = idx.stem().string();
    for (int i = 0; i < 3; ++i) {
        char num[8];
        std::snprintf(num, sizeof(num), "%04d", i);
        const std::string name = stem + "_" + num + ".vtu";
        EXPECT_NE(text.find("<Piece Source=\"" + stem + "/" + name + "\"/>"), std::string::npos)
            << name;
        const std::string piece = (idx.parent_path() / stem / name).string();
        ASSERT_TRUE(fs::exists(piece)) << piece;
        EXPECT_GT(pvtu_num_cells(meshioplusplus::read_vtu(piece)), 0u);
    }
    pvtu_cleanup(path);
}

// The roadmap's "done when": partition -> .pvtu -> read-merge is the original,
// up to point ordering (welding fuses the interface points a partition duplicates).
TEST(Pvtu, PartitionThenPvtuThenReadMergeIsTheOriginalUpToPointOrdering) {
    const Mesh m = pvtu_grid(8);
    meshioplusplus::PartitionOptions o;
    o.mNParts = 4;
    o.mMethod = meshioplusplus::PartitionMethod::SFC;
    const meshioplusplus::PartitionResult result = meshioplusplus::partition(m, o);
    std::vector<const Mesh*> ptrs;
    for (const auto& piece : result.mPieces)
        ptrs.push_back(&piece.mMesh);

    const std::string path = mt::temp_path(".pvtu");
    write_pvtu_pieces_codec(path, ptrs, /*binary=*/true, VtkCodec::Zlib);
    const Mesh welded = pvtu_weld(read_pvtu(path));

    EXPECT_EQ(pvtu_num_cells(welded), pvtu_num_cells(m));
    EXPECT_EQ(pvtu_sorted_points(welded), pvtu_sorted_points(m));
    pvtu_cleanup(path);
}

TEST(Pvtu, NoPartArrayWritesOnePieceAndReadsWithoutARegion) {
    const Mesh m = pvtu_grid(3);
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu(path, m, false, false);
    EXPECT_EQ(pvtu_count(pvtu_slurp(path), "<Piece "), 1u);
    const Mesh back = read_pvtu(path);
    EXPECT_EQ(back.NumRegions(), 0u);
    EXPECT_EQ(back.NumPoints(), m.NumPoints());
    pvtu_cleanup(path);
}

TEST(Pvtu, AnEmptyPartKeyOrAbsentKeyWritesOnePiece) {
    Mesh m = pvtu_labelled(4, 3);
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu_codec(path, m, false, kNone, "");
    EXPECT_EQ(pvtu_count(pvtu_slurp(path), "<Piece "), 1u);
    write_pvtu_codec(path, m, false, kNone, "nope");
    EXPECT_EQ(pvtu_count(pvtu_slurp(path), "<Piece "), 1u);
    write_pvtu_codec(path, m, false, kNone, "partition:part");
    EXPECT_EQ(pvtu_count(pvtu_slurp(path), "<Piece "), 3u);
    pvtu_cleanup(path);
}

TEST(Pvtu, ANonIntegerPartArrayIsRefused) {
    Mesh m = pvtu_grid(3);
    std::vector<NDArray> blocks;
    blocks.push_back(NDArray(DType::Float64, {pvtu_num_cells(m)}));
    m.AddCellData("partition:part", std::move(blocks));
    const std::string path = mt::temp_path(".pvtu");
    try {
        write_pvtu(path, m, false, false);
        FAIL() << "expected a WriteError";
    } catch (const WriteError& e) {
        EXPECT_NE(std::string(e.what()).find("integer part ids"), std::string::npos) << e.what();
    }
    pvtu_cleanup(path);
}

TEST(Pvtu, PartsWithNoCellsAreKeptSoTheRestAreNotRenumbered) {
    Mesh m = pvtu_grid(3);
    const std::size_t nc = pvtu_num_cells(m);
    NDArray labels(DType::Int64, {nc});
    for (std::size_t i = 0; i < nc; i += 2)
        labels.As<std::int64_t>()[i] = 2;  // parts 0 and 2; part 1 owns nothing
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(labels));
    m.AddCellData("partition:part", std::move(blocks));

    const std::string path = mt::temp_path(".pvtu");
    write_pvtu(path, m, false, false);
    EXPECT_EQ(pvtu_count(pvtu_slurp(path), "<Piece "), 3u);
    const Mesh back = read_pvtu(path);
    ASSERT_EQ(back.NumRegions(), 3u);
    EXPECT_EQ(back.Region(1).NumEntries(), 0u);
    EXPECT_EQ(pvtu_num_cells(back), nc);

    ReadOptions o;
    o.mPiece = 1;
    o.mPieceSet = true;
    EXPECT_EQ(pvtu_num_cells(read_pvtu(path, o)), 0u);
    pvtu_cleanup(path);
}

TEST(Pvtu, PieceSelectionResolvesLikeEveryOtherPartitionedFormat) {
    const Mesh m = pvtu_labelled(6, 3);
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu(path, m, false, false);
    const Mesh all = read_pvtu(path);

    std::size_t sum = 0;
    for (int k = 0; k < 3; ++k) {
        ReadOptions o;
        o.mPiece = k;
        o.mPieceSet = true;
        const Mesh one = read_pvtu(path, o);
        EXPECT_EQ(one.NumRegions(), 0u);  // a selected piece carries no region
        EXPECT_EQ(pvtu_num_cells(one), all.Region(static_cast<std::size_t>(k)).NumEntries());
        sum += pvtu_num_cells(one);
    }
    EXPECT_EQ(sum, pvtu_num_cells(m));

    ReadOptions last;
    last.mPiece = -1;
    last.mPieceSet = true;
    EXPECT_EQ(pvtu_num_cells(read_pvtu(path, last)), all.Region(2).NumEntries());

    ReadOptions bad;
    bad.mPiece = 3;
    bad.mPieceSet = true;
    try {
        read_pvtu(path, bad);
        FAIL() << "expected a ReadError";
    } catch (const ReadError& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("piece 3 is out of range"), std::string::npos) << msg;
        EXPECT_NE(msg.find("3 pieces"), std::string::npos) << msg;
    }
    pvtu_cleanup(path);
}

// --- declarations -------------------------------------------------------------

TEST(Pvtu, AMismatchedDeclarationRefusesBeforeWritingAnything) {
    struct Case {
        Mesh mBad;
        const char* mNeedle;
    };
    std::vector<Case> cases;
    cases.push_back({pvtu_piece(DType::Float32),
                     "piece 1 declares point_data 'u' as Float32 with 1 component, but piece 0 "
                     "declares it as Float64 with 1 component"});
    cases.push_back({pvtu_piece(DType::Float64, "u", 3),
                     "point_data 'u' as Float64 with 3 components, but piece 0 declares it as "
                     "Float64 with 1 component"});
    cases.push_back({pvtu_piece(DType::Float64, "v"),
                     "piece 1 is missing point_data 'u', which piece 0 declares"});
    for (const Case& c : cases) {
        const std::string path = mt::temp_path(".pvtu");
        const Mesh good = pvtu_piece();
        try {
            write_pvtu_pieces_codec(path, {&good, &c.mBad}, false, kNone);
            FAIL() << "expected a WriteError for " << c.mNeedle;
        } catch (const WriteError& e) {
            const std::string msg = e.what();
            EXPECT_NE(msg.find(c.mNeedle), std::string::npos) << msg;
            EXPECT_NE(msg.find("every piece of a parallel index must declare identical arrays"),
                      std::string::npos)
                << msg;
        }
        EXPECT_FALSE(fs::exists(path));
        EXPECT_FALSE(fs::exists(fs::path(path).parent_path() / fs::path(path).stem()));
    }
}

TEST(Pvtu, AnEmptyPieceListAndAnExtraArrayAreRefused) {
    const std::string path = mt::temp_path(".pvtu");
    EXPECT_THROW(write_pvtu_pieces_codec(path, {}, false, kNone), WriteError);
    Mesh extra = pvtu_piece();
    extra.AddPointData("w", NDArray(DType::Float64, {3}));
    const Mesh good = pvtu_piece();
    try {
        write_pvtu_pieces_codec(path, {&good, &extra}, false, kNone);
        FAIL() << "expected a WriteError";
    } catch (const WriteError& e) {
        EXPECT_NE(std::string(e.what()).find("declares point_data 'w', which piece 0 does not"),
                  std::string::npos)
            << e.what();
    }
}

TEST(Pvtu, ACellFreePieceMayOmitCellArrays) {
    Mesh idle;
    idle.AssignPoints(NDArray::Uninit(DType::Float64, {0, 3}));
    idle.AddCellBlock("triangle", NDArray::Uninit(DType::Int64, {0, 3}));
    idle.AddPointData("u", NDArray::Uninit(DType::Float64, {0}));
    const Mesh good = pvtu_piece();
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu_pieces_codec(path, {&good, &idle, &good}, false, kNone);
    const Mesh back = read_pvtu(path);
    EXPECT_EQ(pvtu_num_cells(back), 2u);
    ASSERT_EQ(back.NumRegions(), 3u);
    EXPECT_EQ(back.Region(1).NumEntries(), 0u);
    pvtu_cleanup(path);
}

// --- ghosts -------------------------------------------------------------------

namespace {

meshioplusplus::PartitionResult pvtu_ghosted(const Mesh& rMesh, int NParts, int Layers) {
    meshioplusplus::PartitionOptions o;
    o.mNParts = NParts;
    o.mMethod = meshioplusplus::PartitionMethod::SFC;
    o.mGhostLayers = Layers;
    return meshioplusplus::partition(rMesh, o);
}

}  // namespace

TEST(Pvtu, GhostLayersBecomeVtkGhostTypeAndGhostLevel) {
    const Mesh m = pvtu_grid(6);
    const meshioplusplus::PartitionResult result = pvtu_ghosted(m, 3, 1);
    std::vector<const Mesh*> ptrs;
    for (const auto& piece : result.mPieces)
        ptrs.push_back(&piece.mMesh);
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu_pieces_codec(path, ptrs, false, kNone);

    const std::string text = pvtu_slurp(path);
    EXPECT_NE(text.find("<PUnstructuredGrid GhostLevel=\"1\">"), std::string::npos);
    EXPECT_EQ(pvtu_count(text, "<PDataArray type=\"UInt8\" Name=\"vtkGhostType\"/>"), 2u);

    const fs::path idx(path);
    const Mesh piece0 = meshioplusplus::read_vtu(
        (idx.parent_path() / idx.stem() / (idx.stem().string() + "_0000.vtu")).string());
    ASSERT_TRUE(piece0.HasCellData("vtkGhostType"));
    ASSERT_TRUE(piece0.HasPointData("vtkGhostType"));
    EXPECT_EQ(piece0.CellData("vtkGhostType", 0).Dtype(), DType::UInt8);
    EXPECT_EQ(piece0.PointData("vtkGhostType").Dtype(), DType::UInt8);
    EXPECT_TRUE(piece0.HasCellData("partition:ghost"));  // the layer number survives

    // cell flags are exactly (layer > 0)
    const NDArray& layers = result.mPieces[0].mMesh.CellData("partition:ghost", 0);
    const NDArray& flags = piece0.CellData("vtkGhostType", 0);
    ASSERT_EQ(layers.Size(), flags.Size());
    std::size_t ghost_cells = 0;
    for (std::size_t i = 0; i < flags.Size(); ++i) {
        const bool halo = layers.As<std::int64_t>()[i] > 0;
        EXPECT_EQ(flags.As<std::uint8_t>()[i], halo ? 1 : 0);
        ghost_cells += halo ? 1 : 0;
    }
    EXPECT_GT(ghost_cells, 0u);

    // a point is a duplicate exactly when no owned cell of the piece uses it
    std::vector<bool> owned(piece0.NumPoints(), false);
    const auto cb = result.mPieces[0].mMesh.Cells(0);
    for (std::size_t c = 0; c < cb.NumCells(); ++c)
        if (layers.As<std::int64_t>()[c] == 0)
            for (std::size_t v = 0; v < 3; ++v)
                owned[static_cast<std::size_t>(cb.Conn().As<std::int64_t>()[3 * c + v])] = true;
    const NDArray& pflags = piece0.PointData("vtkGhostType");
    for (std::size_t i = 0; i < owned.size(); ++i)
        EXPECT_EQ(pflags.As<std::uint8_t>()[i], owned[i] ? 0 : 1) << "point " << i;
    pvtu_cleanup(path);
}

TEST(Pvtu, GhostsAreKeptByDefaultAndDroppedOnRequest) {
    const Mesh m = pvtu_grid(6);
    const meshioplusplus::PartitionResult result = pvtu_ghosted(m, 3, 1);
    std::vector<const Mesh*> ptrs;
    std::size_t total = 0;
    for (const auto& piece : result.mPieces) {
        ptrs.push_back(&piece.mMesh);
        total += pvtu_num_cells(piece.mMesh);
    }
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu_pieces_codec(path, ptrs, true, VtkCodec::Zlib);

    const Mesh kept = read_pvtu(path);
    EXPECT_EQ(pvtu_num_cells(kept), total);
    EXPECT_GT(total, pvtu_num_cells(m));
    EXPECT_TRUE(kept.HasCellData("vtkGhostType"));
    EXPECT_TRUE(kept.HasPointData("vtkGhostType"));

    ReadOptions drop;
    drop.mGhosts = GhostPolicy::Drop;
    const Mesh dropped = read_pvtu(path, drop);
    EXPECT_EQ(pvtu_num_cells(dropped), pvtu_num_cells(m));
    EXPECT_FALSE(dropped.HasCellData("vtkGhostType"));
    EXPECT_FALSE(dropped.HasPointData("vtkGhostType"));
    EXPECT_FALSE(dropped.HasCellData("partition:ghost"));
    const Mesh welded = pvtu_weld(dropped);
    EXPECT_EQ(pvtu_sorted_points(welded), pvtu_sorted_points(m));

    ReadOptions one;
    one.mPiece = 0;
    one.mPieceSet = true;
    std::size_t owned = 0;
    const NDArray& layers = result.mPieces[0].mMesh.CellData("partition:ghost", 0);
    for (std::size_t i = 0; i < layers.Size(); ++i)
        owned += layers.As<std::int64_t>()[i] == 0 ? 1 : 0;
    one.mGhosts = GhostPolicy::Drop;
    EXPECT_EQ(pvtu_num_cells(read_pvtu(path, one)), owned);
    pvtu_cleanup(path);
}

TEST(Pvtu, ASuppliedVtkGhostTypeIsPassedThroughAndDropsRefinedCells) {
    Mesh piece = pvtu_piece();
    NDArray both(DType::UInt8, {1});
    both.As<std::uint8_t>()[0] = 8;  // REFINEDCELL
    std::vector<NDArray> cell_flags;
    cell_flags.push_back(std::move(both));
    piece.AddCellData("vtkGhostType", std::move(cell_flags));
    const Mesh plain = pvtu_piece();
    Mesh second = pvtu_piece();
    std::vector<NDArray> zero;
    zero.push_back(NDArray(DType::UInt8, {1}));
    second.AddCellData("vtkGhostType", std::move(zero));

    const std::string path = mt::temp_path(".pvtu");
    write_pvtu_pieces_codec(path, {&piece, &second}, false, kNone);
    const Mesh back = read_pvtu(path);
    EXPECT_EQ(pvtu_num_cells(back), 2u);
    ReadOptions drop;
    drop.mGhosts = GhostPolicy::Drop;
    EXPECT_EQ(pvtu_num_cells(read_pvtu(path, drop)), 1u);
    (void)plain;
    pvtu_cleanup(path);
}

TEST(Pvtu, WritingDoesNotMutateTheCallersPieces) {
    const Mesh m = pvtu_grid(6);
    const meshioplusplus::PartitionResult result = pvtu_ghosted(m, 3, 1);
    std::vector<const Mesh*> ptrs;
    for (const auto& piece : result.mPieces)
        ptrs.push_back(&piece.mMesh);
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu_pieces_codec(path, ptrs, false, kNone);
    for (const auto& piece : result.mPieces) {
        EXPECT_FALSE(piece.mMesh.HasCellData("vtkGhostType"));
        EXPECT_FALSE(piece.mMesh.HasPointData("vtkGhostType"));
    }
    pvtu_cleanup(path);
}

// --- reading hand-made and foreign files ----------------------------------------

namespace {

void pvtu_write_text(const fs::path& rPath, const std::string& rText) {
    std::ofstream out(rPath);
    out << rText;
}

std::string pvtu_index(const std::vector<std::string>& rSources,
                       const std::string& rRoot = "PUnstructuredGrid") {
    std::string t = "<?xml version=\"1.0\"?>\n<VTKFile type=\"" + rRoot +
                    "\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n<" +
                    rRoot +
                    " GhostLevel=\"0\">\n<PPoints><PDataArray type=\"Float64\" "
                    "NumberOfComponents=\"3\"/></PPoints>\n";
    for (const std::string& s : rSources)
        t += "<Piece Source=\"" + s + "\"/>\n";
    return t + "</" + rRoot + ">\n</VTKFile>\n";
}

// A directory for a hand-made index; removed on destruction.
struct PvtuDir {
    fs::path mDir;
    PvtuDir() {
        mDir = fs::path(mt::temp_path(".d"));
        fs::create_directories(mDir);
    }
    ~PvtuDir() {
        std::error_code ec;
        fs::remove_all(mDir, ec);
    }
    std::string operator/(const std::string& rName) const { return (mDir / rName).string(); }
};

}  // namespace

TEST(Pvtu, TheHeaderTypeOfTheIndexDoesNotGovernThePieces) {
    // The index says header_type="UInt64"; the piece was written with the
    // default UInt32. The index carries no arrays, so each piece's own header
    // is the authority.
    PvtuDir d;
    meshioplusplus::write_vtu(d / "a.vtu", pvtu_piece(), true, false);
    pvtu_write_text(d.mDir / "i.pvtu", pvtu_index({"a.vtu"}));
    EXPECT_EQ(pvtu_num_cells(read_pvtu(d / "i.pvtu")), 1u);
}

TEST(Pvtu, PiecePathsWithSpacesAndEntitiesResolve) {
    PvtuDir d;
    fs::create_directories(d.mDir / "my pieces & more");
    meshioplusplus::write_vtu(d / "my pieces & more/piece one.vtu", pvtu_piece(), false, false);
    pvtu_write_text(d.mDir / "i.pvtu", pvtu_index({"my pieces &amp; more/piece one.vtu"}));
    EXPECT_EQ(pvtu_num_cells(read_pvtu(d / "i.pvtu")), 1u);

    const Mesh m = pvtu_labelled(4, 2);
    write_pvtu(d / "my case.pvtu", m, false, false);
    EXPECT_TRUE(fs::exists(d.mDir / "my case" / "my case_0000.vtu"));
    EXPECT_EQ(pvtu_num_cells(read_pvtu(d / "my case.pvtu")), pvtu_num_cells(m));
}

TEST(Pvtu, AMissingPieceNamesTheAttributeAndTheFile) {
    PvtuDir d;
    pvtu_write_text(d.mDir / "i.pvtu", pvtu_index({"gone/piece.vtu"}));
    try {
        read_pvtu(d / "i.pvtu");
        FAIL() << "expected a ReadError";
    } catch (const ReadError& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("gone/piece.vtu"), std::string::npos) << msg;
        EXPECT_NE(msg.find("Source="), std::string::npos) << msg;
        EXPECT_NE(msg.find("does not exist"), std::string::npos) << msg;
    }
}

TEST(Pvtu, AnAbsolutePathFromAnotherMachineIsAnErrorNotAGuess) {
    PvtuDir d;
    meshioplusplus::write_vtu(d / "piece.vtu", pvtu_piece(), false, false);  // must NOT be read
    pvtu_write_text(d.mDir / "i.pvtu", pvtu_index({"/scratch/run7/piece.vtu"}));
    try {
        read_pvtu(d / "i.pvtu");
        FAIL() << "expected a ReadError";
    } catch (const ReadError& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("/scratch/run7/piece.vtu"), std::string::npos) << msg;
        EXPECT_NE(msg.find("not searched"), std::string::npos) << msg;
    }
}

TEST(Pvtu, APieceThatIsNotVtuOrVtpIsRefusedByNameAndAnIndexIsNotAPiece) {
    PvtuDir d;
    pvtu_write_text(d.mDir / "m.msh", "$MeshFormat\n");
    pvtu_write_text(d.mDir / "i.pvtu", pvtu_index({"m.msh"}));
    try {
        read_pvtu(d / "i.pvtu");
        FAIL() << "expected a ReadError";
    } catch (const ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("unsupported piece"), std::string::npos) << e.what();
    }
    // No recursion inside a parallel index, so no cycle is expressible.
    pvtu_write_text(d.mDir / "self.pvtu", pvtu_index({"self.pvtu"}));
    EXPECT_THROW(read_pvtu(d / "self.pvtu"), ReadError);
}

TEST(Pvtu, EmptyIndexReadsAsAnEmptyMeshAndWrongRootIsRefused) {
    PvtuDir d;
    pvtu_write_text(d.mDir / "e.pvtu", pvtu_index({}));
    const Mesh empty = read_pvtu(d / "e.pvtu");
    EXPECT_EQ(empty.NumPoints(), 0u);
    EXPECT_EQ(empty.NumCellBlocks(), 0u);

    pvtu_write_text(d.mDir / "w.pvtu", pvtu_index({}, "PPolyData"));
    try {
        read_pvtu(d / "w.pvtu");
        FAIL() << "expected a ReadError";
    } catch (const ReadError& e) {
        EXPECT_NE(std::string(e.what()).find("expected type PUnstructuredGrid"), std::string::npos)
            << e.what();
    }
    pvtu_write_text(d.mDir / "x.pvtu", "not xml");
    EXPECT_THROW(read_pvtu(d / "x.pvtu"), ReadError);
}

TEST(Pvtu, MetadataAgreesWithARealRead) {
    const Mesh m = pvtu_labelled(6, 3);
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu(path, m, true, true);
    const Mesh back = read_pvtu(path);
    const meshioplusplus::MeshMetadata meta = read_pvtu_metadata(path);

    EXPECT_EQ(meta.mFormat, "pvtu");
    EXPECT_EQ(meta.mNumPoints, back.NumPoints());
    ASSERT_EQ(meta.mCellBlocks.size(), back.NumCellBlocks());
    for (std::size_t b = 0; b < back.NumCellBlocks(); ++b) {
        EXPECT_EQ(meta.mCellBlocks[b].mType, back.Cells(b).Type());
        EXPECT_EQ(meta.mCellBlocks[b].mNumCells, back.Cells(b).NumCells());
    }
    EXPECT_EQ(meta.mPointDataNames, back.PointDataNames());
    EXPECT_EQ(meta.mCellDataNames, back.CellDataNames());
    pvtu_cleanup(path);
}

// --- field data and the ghost option through the registry ------------------------

TEST(Pvtu, FieldDataIsTheUnionAcrossPiecesNotNamespaced) {
    // Field data belongs to the dataset: every piece repeats it, and merge() would
    // rename the copies `0:TimeValue`, `1:TimeValue`, ... which nothing looks up.
    Mesh a = pvtu_piece(), b = pvtu_piece(), c = pvtu_piece();
    for (Mesh* m : {&a, &b, &c}) {
        NDArray t(DType::Float64, {1});
        *t.As<double>() = 0.25;
        m->AddFieldData("TimeValue", std::move(t));
    }
    NDArray extra(DType::Float64, {2});
    b.AddFieldData("only_in_b", std::move(extra));

    const std::string path = mt::temp_path(".pvtu");
    write_pvtu_pieces_codec(path, {&a, &b, &c}, false, kNone);
    const Mesh back = read_pvtu(path);
    EXPECT_EQ(back.FieldDataNames(), (std::vector<std::string>{"TimeValue", "only_in_b"}));
    EXPECT_EQ(*back.FieldData("TimeValue").As<double>(), 0.25);
    // ... and the summary names exactly what the read returns.
    EXPECT_EQ(read_pvtu_metadata(path).mFieldDataNames, back.FieldDataNames());
    pvtu_cleanup(path);
}

TEST(Pvtu, TheGhostOptionIsAReadOptionsMemberSoTheRegistryHonoursIt) {
    const Mesh m = pvtu_grid(6);
    const meshioplusplus::PartitionResult result = pvtu_ghosted(m, 3, 1);
    std::vector<const Mesh*> ptrs;
    for (const auto& piece : result.mPieces)
        ptrs.push_back(&piece.mMesh);
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu_pieces_codec(path, ptrs, false, kNone);

    ReadOptions keep;  // default: Keep -- and what a zero-initialized struct means
    const Mesh kept = meshioplusplus::registry_read(path, "pvtu", keep);
    ReadOptions drop;
    drop.mGhosts = GhostPolicy::Drop;
    const Mesh dropped = meshioplusplus::registry_read(path, "pvtu", drop);
    EXPECT_GT(pvtu_num_cells(kept), pvtu_num_cells(dropped));
    EXPECT_EQ(pvtu_num_cells(dropped), pvtu_num_cells(m));
    EXPECT_FALSE(dropped.HasCellData("vtkGhostType"));

    // Every other reader ignores it: a plain file with no halo is already the answer.
    const std::string vtu = mt::temp_path(".vtu");
    meshioplusplus::write_vtu(vtu, m, false, false);
    EXPECT_EQ(pvtu_num_cells(meshioplusplus::registry_read(vtu, "vtu", drop)), pvtu_num_cells(m));
    pvtu_cleanup(path);
}

// --- polyhedra ------------------------------------------------------------------

namespace {

// Two tetrahedra sharing a face, as polyhedra, with an Int64 cell array `pName`.
Mesh pvtu_two_polyhedra(const char* pName, std::int64_t First, std::int64_t Second) {
    NDArray pts(DType::Float64, {5, 3});
    const double xyz[5][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}};
    for (std::size_t i = 0; i < 5; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            pts.As<double>()[3 * i + j] = xyz[i][j];
    using Face = std::vector<std::int64_t>;
    const std::vector<std::vector<Face>> cells = {
        {{0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {0, 2, 3}},
        {{1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}},
    };
    Mesh m;
    m.AssignPoints(std::move(pts));
    m.AddPolyhedronBlock("polyhedron4", cells);
    NDArray ids(DType::Int64, {2});
    ids.As<std::int64_t>()[0] = First;
    ids.As<std::int64_t>()[1] = Second;
    std::vector<NDArray> blocks;
    blocks.push_back(std::move(ids));
    m.AddCellData(pName, std::move(blocks));
    return m;
}

}  // namespace

TEST(Pvtu, PolyhedronBlocksAreCarvedGhostedAndDropped) {
    const std::string path = mt::temp_path(".pvtu");
    write_pvtu(path, pvtu_two_polyhedra("partition:part", 0, 1), false, false);
    const Mesh back = read_pvtu(path);
    EXPECT_EQ(pvtu_num_cells(back), 2u);
    ASSERT_EQ(back.NumRegions(), 2u);
    EXPECT_TRUE(back.Cells(0).IsPolyhedron());
    pvtu_cleanup(path);

    const Mesh ghosted = pvtu_two_polyhedra("partition:ghost", 0, 1);
    write_pvtu_pieces_codec(path, {&ghosted}, false, kNone);
    const Mesh kept = read_pvtu(path);
    EXPECT_EQ(pvtu_num_cells(kept), 2u);
    EXPECT_TRUE(kept.HasCellData("vtkGhostType"));
    ReadOptions drop;
    drop.mGhosts = GhostPolicy::Drop;
    const Mesh dropped = read_pvtu(path, drop);
    EXPECT_EQ(pvtu_num_cells(dropped), 1u);
    EXPECT_FALSE(dropped.HasCellData("vtkGhostType"));
    pvtu_cleanup(path);
}

// --- .pvtp ----------------------------------------------------------------------

namespace {

Mesh pvtu_quads() {
    Mesh m = mt::quad_mesh();
    std::vector<NDArray> parts;
    NDArray p(DType::Int64, {m.Cells(0).NumCells()});
    for (std::size_t i = 0; i < p.Size(); ++i)
        p.As<std::int64_t>()[i] = static_cast<std::int64_t>(i % 2);
    parts.push_back(std::move(p));
    m.AddCellData("partition:part", std::move(parts));
    return m;
}

}  // namespace

TEST(Pvtp, RoundTripsOverVtpPieces) {
    const Mesh m = pvtu_quads();
    const std::string path = mt::temp_path(".pvtp");
    write_pvtp(path, m, true, false);
    const std::string text = pvtu_slurp(path);
    EXPECT_NE(text.find("<VTKFile type=\"PPolyData\""), std::string::npos);
    EXPECT_NE(text.find("<PPolyData GhostLevel=\"0\">"), std::string::npos);
    const fs::path idx(path);
    EXPECT_NE(
        text.find("Source=\"" + idx.stem().string() + "/" + idx.stem().string() + "_0000.vtp\""),
        std::string::npos);

    const Mesh back = read_pvtp(path);
    EXPECT_EQ(pvtu_num_cells(back), pvtu_num_cells(m));
    EXPECT_GE(back.NumRegions(), 1u);
    EXPECT_EQ(read_pvtp_metadata(path).mFormat, "pvtp");
    pvtu_cleanup(path);
}

TEST(Pvtp, SharesTheDeclarationCheck) {
    Mesh good = pvtu_quads();
    Mesh bad = pvtu_quads();
    bad.AddPointData("u", NDArray(DType::Float32, {bad.NumPoints()}));
    good.AddPointData("u", NDArray(DType::Float64, {good.NumPoints()}));
    const std::string path = mt::temp_path(".pvtp");
    try {
        write_pvtp_pieces_codec(path, {&good, &bad}, false, kNone);
        FAIL() << "expected a WriteError";
    } catch (const WriteError& e) {
        EXPECT_NE(std::string(e.what()).find("pvtp: piece 1 declares point_data 'u'"),
                  std::string::npos)
            << e.what();
    }
}
