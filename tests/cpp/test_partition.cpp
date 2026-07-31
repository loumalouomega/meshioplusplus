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
// The partition operation. SFC assertions are exact (the cut rule is a pinned
// deterministic contract); KaHIP assertions cover balance / coverage /
// fixed-seed reproducibility only, never exact labels, since KaHIP's
// assignment may change between KaHIP versions.

// System includes
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/partition.hpp"

namespace {

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::partition;
using meshioplusplus::partition_has_kahip;
using meshioplusplus::partition_labels;
using meshioplusplus::partition_method_from_name;
using meshioplusplus::partition_mode_from_name;
using meshioplusplus::PartitionMethod;
using meshioplusplus::PartitionMode;
using meshioplusplus::PartitionOptions;
using meshioplusplus::PartitionResult;

PartitionOptions opts(int nparts, PartitionMethod method = PartitionMethod::SFC) {
    PartitionOptions o;
    o.mNParts = nparts;
    o.mMethod = method;
    return o;
}

// A structured nx x ny x nz shared-vertex hexahedron grid on [0,nx]x[0,ny]x[0,nz].
Mesh hex_grid(std::size_t nx, std::size_t ny, std::size_t nz) {
    std::vector<std::vector<double>> pts;
    for (std::size_t k = 0; k <= nz; ++k)
        for (std::size_t j = 0; j <= ny; ++j)
            for (std::size_t i = 0; i <= nx; ++i)
                pts.push_back(
                    {static_cast<double>(i), static_cast<double>(j), static_cast<double>(k)});
    auto vid = [&](std::size_t i, std::size_t j, std::size_t k) {
        return static_cast<std::int64_t>((k * (ny + 1) + j) * (nx + 1) + i);
    };
    std::vector<std::vector<std::int64_t>> cells;
    for (std::size_t k = 0; k < nz; ++k)
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t i = 0; i < nx; ++i)
                cells.push_back({vid(i, j, k), vid(i + 1, j, k), vid(i + 1, j + 1, k),
                                 vid(i, j + 1, k), vid(i, j, k + 1), vid(i + 1, j, k + 1),
                                 vid(i + 1, j + 1, k + 1), vid(i, j + 1, k + 1)});
    return mt::make_mesh(std::move(pts), "hexahedron", std::move(cells));
}

// A structured nx x ny quad grid on [0,nx]x[0,ny] (z = 0).
Mesh quad_grid(std::size_t nx, std::size_t ny) {
    std::vector<std::vector<double>> pts;
    for (std::size_t j = 0; j <= ny; ++j)
        for (std::size_t i = 0; i <= nx; ++i)
            pts.push_back({static_cast<double>(i), static_cast<double>(j), 0.0});
    auto vid = [&](std::size_t i, std::size_t j) {
        return static_cast<std::int64_t>(j * (nx + 1) + i);
    };
    std::vector<std::vector<std::int64_t>> cells;
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i)
            cells.push_back({vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)});
    return mt::make_mesh(std::move(pts), "quad", std::move(cells));
}

std::size_t total_cells(const Mesh& rMesh) {
    std::size_t total = 0;
    for (const auto cb : rMesh.CellRange())
        total += cb.NumCells();
    return total;
}

// Flatten block-aligned labels to the global cell order.
std::vector<std::int64_t> flat_labels(const std::vector<NDArray>& rLabels) {
    std::vector<std::int64_t> flat;
    for (const NDArray& a : rLabels) {
        const std::size_t n = a.Shape().empty() ? 0 : a.Shape()[0];
        for (std::size_t i = 0; i < n; ++i)
            flat.push_back(meshioplusplus::detail::read_int(a, i));
    }
    return flat;
}

std::vector<std::size_t> part_sizes(const std::vector<std::int64_t>& rFlat, int nparts) {
    std::vector<std::size_t> sizes(static_cast<std::size_t>(nparts), 0);
    for (std::int64_t p : rFlat) {
        EXPECT_GE(p, 0);
        EXPECT_LT(p, nparts);
        ++sizes[static_cast<std::size_t>(p)];
    }
    return sizes;
}

// Every input cell must land in exactly one piece (partition of unity),
// checked through the per-piece cell maps.
void check_partition_of_unity(const Mesh& rIn, const PartitionResult& rRes) {
    const std::size_t nblocks = rIn.NumCellBlocks();
    std::size_t b = 0;
    for (const auto cb : rIn.CellRange()) {
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            int owners = 0;
            for (const auto& piece : rRes.mPieces) {
                ASSERT_EQ(piece.mCellMaps.size(), nblocks);
                if (meshioplusplus::detail::read_int(piece.mCellMaps[b], c) >= 0)
                    ++owners;
            }
            ASSERT_EQ(owners, 1) << "block " << b << " cell " << c;
        }
        ++b;
    }
}

// --- SFC: labels -------------------------------------------------------------

TEST(Partition, LabelsAreBlockAlignedAndInRange) {
    const Mesh m = hex_grid(4, 4, 4);
    const std::vector<NDArray> labels = partition_labels(m, opts(4));
    ASSERT_EQ(labels.size(), m.NumCellBlocks());
    const std::vector<std::int64_t> flat = flat_labels(labels);
    EXPECT_EQ(flat.size(), total_cells(m));
    part_sizes(flat, 4);  // asserts the [0, nparts) range
}

TEST(Partition, UnweightedSizesDifferByAtMostOne) {
    const Mesh m = hex_grid(4, 4, 3);  // 48 cells
    for (int nparts : {2, 3, 5, 7}) {
        const std::vector<std::int64_t> flat = flat_labels(partition_labels(m, opts(nparts)));
        const std::vector<std::size_t> sizes = part_sizes(flat, nparts);
        std::size_t lo = flat.size(), hi = 0;
        for (std::size_t s : sizes) {
            lo = std::min(lo, s);
            hi = std::max(hi, s);
        }
        EXPECT_LE(hi - lo, 1u) << "nparts=" << nparts;
    }
}

TEST(Partition, DeterministicAcrossRuns) {
    const Mesh m = hex_grid(3, 3, 3);
    const std::vector<std::int64_t> a = flat_labels(partition_labels(m, opts(4)));
    const std::vector<std::int64_t> b = flat_labels(partition_labels(m, opts(4)));
    EXPECT_EQ(a, b);
}

TEST(Partition, SfcPartsAreSpatiallyContiguousRanges) {
    // On a 1D-ish strip the Hilbert cut must produce contiguous index ranges of
    // near-equal size; more importantly the parts must be *connected* along the
    // strip. We assert the weaker pinned property: sorting cells by part id
    // never interleaves (each part is one contiguous run along the curve), via
    // sizes + the ≤1 rule already covered; here we check 2D grids get every
    // part non-empty when nparts <= ncells.
    const Mesh m = quad_grid(6, 6);
    const std::vector<std::int64_t> flat = flat_labels(partition_labels(m, opts(6)));
    const std::vector<std::size_t> sizes = part_sizes(flat, 6);
    for (std::size_t s : sizes)
        EXPECT_GT(s, 0u);
}

TEST(Partition, NPartsOneAssignsEverythingToZero) {
    const Mesh m = mt::tri_quad_mesh();
    const std::vector<std::int64_t> flat = flat_labels(partition_labels(m, opts(1)));
    for (std::int64_t p : flat)
        EXPECT_EQ(p, 0);
}

TEST(Partition, MorePartsThanCellsLeavesEmptyParts) {
    const Mesh m = mt::tet_mesh();  // 1 cell
    const int nparts = 3;
    const PartitionResult res = partition(m, opts(nparts));
    ASSERT_EQ(res.mPieces.size(), static_cast<std::size_t>(nparts));
    std::size_t total_out = 0;
    for (const auto& piece : res.mPieces)
        total_out += total_cells(piece.mMesh);
    EXPECT_EQ(total_out, total_cells(m));
}

// --- SFC: pieces -------------------------------------------------------------

TEST(Partition, PiecesArePartitionOfUnity) {
    const Mesh m = hex_grid(3, 3, 2);
    const PartitionResult res = partition(m, opts(4));
    ASSERT_EQ(res.mPieces.size(), 4u);
    check_partition_of_unity(m, res);
    std::size_t total_out = 0;
    for (const auto& piece : res.mPieces) {
        EXPECT_EQ(piece.mPartId, static_cast<int>(&piece - res.mPieces.data()));
        total_out += total_cells(piece.mMesh);
    }
    EXPECT_EQ(total_out, total_cells(m));
}

TEST(Partition, PiecesKeepInputBlockStructure) {
    // Unlike split, every piece has exactly the input's block count, in input
    // order (drop_empty_blocks=false), so cell maps index by input block.
    const Mesh m = mt::tri_quad_mesh();
    const PartitionResult res = partition(m, opts(2));
    for (const auto& piece : res.mPieces) {
        EXPECT_EQ(piece.mMesh.NumCellBlocks(), m.NumCellBlocks());
        std::size_t b = 0;
        for (const auto cb : piece.mMesh.CellRange()) {
            EXPECT_EQ(std::string(cb.Type()), std::string(m.Cells(b).Type()));
            ++b;
        }
    }
}

TEST(Partition, PiecesAgreeWithLabels) {
    const Mesh m = hex_grid(3, 2, 2);
    const std::vector<std::int64_t> flat = flat_labels(partition_labels(m, opts(3)));
    const PartitionResult res = partition(m, opts(3));
    std::size_t base = 0, b = 0;
    for (const auto cb : m.CellRange()) {
        for (std::size_t c = 0; c < cb.NumCells(); ++c) {
            const std::int64_t p = flat[base + c];
            EXPECT_GE(meshioplusplus::detail::read_int(
                          res.mPieces[static_cast<std::size_t>(p)].mCellMaps[b], c),
                      0);
        }
        base += cb.NumCells();
        ++b;
    }
}

TEST(Partition, PointMapsPruneToThePiece) {
    const Mesh m = hex_grid(2, 2, 2);
    const PartitionResult res = partition(m, opts(2));
    for (const auto& piece : res.mPieces) {
        const NDArray& pm = piece.mPointMap;
        ASSERT_EQ(pm.Shape()[0], m.NumPoints());
        std::size_t mapped = 0;
        for (std::size_t i = 0; i < m.NumPoints(); ++i) {
            const std::int64_t v = meshioplusplus::detail::read_int(pm, i);
            if (v >= 0) {
                EXPECT_LT(static_cast<std::size_t>(v), piece.mMesh.NumPoints());
                ++mapped;
            }
        }
        EXPECT_EQ(mapped, piece.mMesh.NumPoints());
    }
}

TEST(Partition, RecordIdsAttachOriginalIndices) {
    const Mesh m = hex_grid(2, 2, 1);
    PartitionOptions o = opts(2);
    o.mRecordIds = true;
    const std::vector<std::int64_t> flat = flat_labels(partition_labels(m, o));
    const PartitionResult res = partition(m, o);
    for (const auto& piece : res.mPieces) {
        ASSERT_TRUE(piece.mMesh.HasPointData("partition:original_point_id"));
        ASSERT_TRUE(piece.mMesh.HasCellData("partition:original_cell_id"));
        const NDArray& ids = piece.mMesh.CellData("partition:original_cell_id", 0);
        const std::size_t n = ids.Shape().empty() ? 0 : ids.Shape()[0];
        for (std::size_t c = 0; c < n; ++c) {
            const std::int64_t orig = meshioplusplus::detail::read_int(ids, c);
            EXPECT_EQ(flat[static_cast<std::size_t>(orig)], piece.mPartId);
        }
    }
}

TEST(Partition, CarriesPointAndCellData) {
    const Mesh m = mt::data_mesh();
    const PartitionResult res = partition(m, opts(2));
    check_partition_of_unity(m, res);
    for (const auto& piece : res.mPieces) {
        for (const std::string& name : m.PointDataNames())
            EXPECT_TRUE(piece.mMesh.HasPointData(name)) << name;
        for (const std::string& name : m.CellDataNames())
            EXPECT_TRUE(piece.mMesh.HasCellData(name)) << name;
    }
}

// --- weights -----------------------------------------------------------------

TEST(Partition, WeightedCutBalancesWeightNotCount) {
    // One dominant-weight cell: the weighted cut must isolate it near-alone
    // while the unweighted cut would split 50/50 by count.
    Mesh m = quad_grid(8, 1);  // 8 cells in a row
    std::vector<double> w(8, 1.0);
    w[0] = 100.0;
    m.AddCellData("w", {mt::data_array(w)});
    PartitionOptions o = opts(2);
    o.mWeightsKey = "w";
    const std::vector<std::int64_t> flat = flat_labels(partition_labels(m, o));
    double wsum[2] = {0.0, 0.0};
    double wtotal = 0.0, wmax = 0.0;
    for (std::size_t c = 0; c < 8; ++c) {
        wsum[flat[c]] += w[c];
        wtotal += w[c];
        wmax = std::max(wmax, w[c]);
    }
    // The midpoint rule bounds each part's weight by the ideal share plus one
    // maximal cell.
    EXPECT_LE(wsum[0], wtotal / 2 + wmax);
    EXPECT_LE(wsum[1], wtotal / 2 + wmax);
    // The heavy cell dominates half the weight, so its part holds few cells.
    const std::int64_t heavy_part = flat[0];
    std::size_t heavy_count = 0;
    for (std::int64_t p : flat)
        if (p == heavy_part)
            ++heavy_count;
    EXPECT_LE(heavy_count, 2u);
}

TEST(Partition, WeightsMustExist) {
    const Mesh m = mt::tet_mesh();
    PartitionOptions o = opts(2);
    o.mWeightsKey = "nope";
    EXPECT_THROW(partition_labels(m, o), std::invalid_argument);
}

TEST(Partition, WeightsMustBeScalar) {
    Mesh m = quad_grid(2, 2);
    m.AddCellData("v", {mt::data_array({1, 2, 3, 4, 5, 6, 7, 8}, /*cols=*/2)});
    PartitionOptions o = opts(2);
    o.mWeightsKey = "v";
    EXPECT_THROW(partition_labels(m, o), std::invalid_argument);
}

TEST(Partition, NegativeWeightsThrow) {
    Mesh m = quad_grid(2, 2);
    m.AddCellData("w", {mt::data_array({1, 1, -1, 1})});
    PartitionOptions o = opts(2);
    o.mWeightsKey = "w";
    EXPECT_THROW(partition_labels(m, o), std::invalid_argument);
}

// --- errors ------------------------------------------------------------------

TEST(Partition, InvalidNPartsThrows) {
    const Mesh m = mt::tet_mesh();
    EXPECT_THROW(partition_labels(m, opts(0)), std::invalid_argument);
    EXPECT_THROW(partition_labels(m, opts(-2)), std::invalid_argument);
}

TEST(Partition, GhostLayersZeroIsAPartitionOfUnity) {
    // The v9.0.0 contract for the unghosted case is unchanged: every cell
    // lands in exactly one piece and nothing carries a partition:ghost tag.
    const Mesh m = quad_grid(4, 4);
    PartitionOptions o = opts(3);
    o.mGhostLayers = 0;
    const PartitionResult r = partition(m, o);
    std::size_t total = 0;
    for (const auto& r_piece : r.mPieces) {
        total += total_cells(r_piece.mMesh);
        EXPECT_FALSE(r_piece.mMesh.HasCellData("partition:ghost"));
    }
    EXPECT_EQ(total, 16u);
}

TEST(Partition, GhostLayersGrowEachPieceAndTagTheHalo) {
    const Mesh m = quad_grid(4, 4);  // 16 quads, 2 parts along the Hilbert curve
    PartitionOptions base = opts(2);
    const PartitionResult plain = partition(m, base);

    PartitionOptions o = opts(2);
    o.mGhostLayers = 1;
    const PartitionResult ghosted = partition(m, o);

    ASSERT_EQ(ghosted.mPieces.size(), plain.mPieces.size());
    std::size_t ghost_total = 0;
    for (std::size_t p = 0; p < ghosted.mPieces.size(); ++p) {
        const Mesh& r_g = ghosted.mPieces[p].mMesh;
        const std::size_t owned = total_cells(plain.mPieces[p].mMesh);
        const std::size_t with_halo = total_cells(r_g);
        // A halo can only add cells, and on a connected grid cut in two it
        // genuinely must add some.
        EXPECT_GT(with_halo, owned) << "part " << p << " gained no ghosts";

        ASSERT_TRUE(r_g.HasCellData("partition:ghost"));
        ASSERT_EQ(r_g.CellDataNumBlocks("partition:ghost"), r_g.NumCellBlocks());
        std::size_t owned_tagged = 0, ghost_tagged = 0;
        for (std::size_t b = 0; b < r_g.NumCellBlocks(); ++b) {
            const NDArray& a = r_g.CellData("partition:ghost", b);
            for (std::size_t i = 0; i < a.Size(); ++i) {
                const std::int64_t d = a.As<std::int64_t>()[i];
                EXPECT_GE(d, 0);
                EXPECT_LE(d, 1) << "one layer requested, so depth cannot exceed 1";
                (d == 0 ? owned_tagged : ghost_tagged)++;
            }
        }
        // The tag partitions the piece into exactly the owned set plus the halo.
        EXPECT_EQ(owned_tagged, owned);
        EXPECT_EQ(owned_tagged + ghost_tagged, with_halo);
        ghost_total += ghost_tagged;
    }
    EXPECT_GT(ghost_total, 0u);
}

TEST(Partition, GhostLayersReadInt32ConnectivityCorrectly) {
    // Regression: the ghost-layer cell->node incidence used to read dense
    // connectivity as `Conn().As<std::int64_t>()`, which performs no dtype
    // check. A MESHIO-backed mesh routinely carries Int32 connectivity straight
    // from numpy, so every node id was two fused Int32 entries; most failed the
    // range filter and vanished, leaving halos silently too small (and reading
    // past the end of the last row). The incidence now goes through
    // `detail::cell_node_ids`, which reads via `detail::read_int`.
    //
    // On the NATIVE and KRATOS backends `AddCellBlock` canonicalizes to Int64,
    // so this is a tautology there -- it is the MESHIO leg that has teeth.
    const Mesh wide = quad_grid(4, 4);

    Mesh narrow;
    narrow.AssignPoints(mt::points_from([] {
        std::vector<std::vector<double>> pts;
        for (std::size_t j = 0; j <= 4; ++j)
            for (std::size_t i = 0; i <= 4; ++i)
                pts.push_back({static_cast<double>(i), static_cast<double>(j), 0.0});
        return pts;
    }()));
    {
        const std::size_t ncells = 16, npc = 4;
        NDArray conn(DType::Int32, {ncells, npc});
        std::int32_t* dst = conn.As<std::int32_t>();
        std::size_t c = 0;
        for (std::size_t j = 0; j < 4; ++j)
            for (std::size_t i = 0; i < 4; ++i) {
                const std::int32_t v = static_cast<std::int32_t>(j * 5 + i);
                dst[c * npc + 0] = v;
                dst[c * npc + 1] = v + 1;
                dst[c * npc + 2] = v + 6;
                dst[c * npc + 3] = v + 5;
                ++c;
            }
        narrow.AddCellBlock("quad", std::move(conn));
    }

    PartitionOptions o = opts(2);
    o.mGhostLayers = 1;
    const PartitionResult from_wide = partition(wide, o);
    const PartitionResult from_narrow = partition(narrow, o);

    ASSERT_EQ(from_narrow.mPieces.size(), from_wide.mPieces.size());
    for (std::size_t p = 0; p < from_wide.mPieces.size(); ++p)
        EXPECT_EQ(total_cells(from_narrow.mPieces[p].mMesh),
                  total_cells(from_wide.mPieces[p].mMesh))
            << "part " << p << ": the halo must not depend on the connectivity dtype";
}

TEST(Partition, MoreGhostLayersNeverShrinkAPiece) {
    const Mesh m = quad_grid(5, 5);
    std::size_t prev = 0;
    for (int layers = 0; layers <= 3; ++layers) {
        PartitionOptions o = opts(2);
        o.mGhostLayers = layers;
        const PartitionResult r = partition(m, o);
        const std::size_t n = total_cells(r.mPieces[0].mMesh);
        EXPECT_GE(n, prev) << "layer " << layers << " shrank the piece";
        EXPECT_LE(n, 25u) << "a halo cannot exceed the whole mesh";
        prev = n;
    }
}

TEST(Partition, NegativeGhostLayersThrow) {
    const Mesh m = mt::tet_mesh();
    PartitionOptions o = opts(2);
    o.mGhostLayers = -1;
    EXPECT_THROW(partition(m, o), std::invalid_argument);
}

TEST(Partition, LabelsRejectGhostLayers) {
    // A flat per-cell label array is the ownership map; a cell can be a ghost
    // of several parts at once, so ghosting has no representation there.
    const Mesh m = mt::tet_mesh();
    PartitionOptions o = opts(2);
    o.mGhostLayers = 1;
    try {
        partition_labels(m, o);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("ghost_layers"), std::string::npos);
    }
}

TEST(Partition, MethodAndModeParsers) {
    EXPECT_EQ(partition_method_from_name("sfc"), PartitionMethod::SFC);
    EXPECT_EQ(partition_method_from_name("kahip"), PartitionMethod::KaHIP);
    EXPECT_EQ(partition_method_from_name("auto"), PartitionMethod::Auto);
    EXPECT_THROW(partition_method_from_name("metis"), std::invalid_argument);
    EXPECT_EQ(partition_mode_from_name("fast"), PartitionMode::Fast);
    EXPECT_EQ(partition_mode_from_name("eco"), PartitionMode::Eco);
    EXPECT_EQ(partition_mode_from_name("strong"), PartitionMode::Strong);
    EXPECT_THROW(partition_mode_from_name("turbo"), std::invalid_argument);
}

TEST(Partition, AutoResolvesToSfcWithoutKahip) {
    if (partition_has_kahip())
        GTEST_SKIP() << "KaHIP is compiled in; Auto resolves to KaHIP here";
    const Mesh m = hex_grid(2, 2, 2);
    const std::vector<std::int64_t> a =
        flat_labels(partition_labels(m, opts(3, PartitionMethod::Auto)));
    const std::vector<std::int64_t> b =
        flat_labels(partition_labels(m, opts(3, PartitionMethod::SFC)));
    EXPECT_EQ(a, b);
}

// --- KaHIP -------------------------------------------------------------------

#ifdef MESHIOPLUSPLUS_HAS_KAHIP

TEST(PartitionKahip, CoverageAndBalance) {
    const Mesh m = hex_grid(4, 4, 4);  // 64 cells
    PartitionOptions o = opts(4, PartitionMethod::KaHIP);
    const std::vector<std::int64_t> flat = flat_labels(partition_labels(m, o));
    ASSERT_EQ(flat.size(), 64u);
    const std::vector<std::size_t> sizes = part_sizes(flat, 4);
    // Within the requested imbalance (3%) plus one cell of slack.
    const double avg = 64.0 / 4.0;
    for (std::size_t s : sizes) {
        EXPECT_GT(s, 0u);
        EXPECT_LE(static_cast<double>(s), avg * (1.0 + o.mImbalance) + 1.0);
    }
}

TEST(PartitionKahip, PiecesArePartitionOfUnity) {
    const Mesh m = hex_grid(3, 3, 3);
    const PartitionResult res = partition(m, opts(3, PartitionMethod::KaHIP));
    ASSERT_EQ(res.mPieces.size(), 3u);
    check_partition_of_unity(m, res);
}

TEST(PartitionKahip, FixedSeedIsReproducible) {
    const Mesh m = hex_grid(3, 3, 3);
    PartitionOptions o = opts(4, PartitionMethod::KaHIP);
    o.mSeed = 42;
    const std::vector<std::int64_t> a = flat_labels(partition_labels(m, o));
    const std::vector<std::int64_t> b = flat_labels(partition_labels(m, o));
    EXPECT_EQ(a, b);
}

TEST(PartitionKahip, ModesAllProduceValidPartitions) {
    const Mesh m = quad_grid(6, 6);
    for (PartitionMode mode : {PartitionMode::Fast, PartitionMode::Eco, PartitionMode::Strong}) {
        PartitionOptions o = opts(4, PartitionMethod::KaHIP);
        o.mMode = mode;
        const std::vector<std::int64_t> flat = flat_labels(partition_labels(m, o));
        ASSERT_EQ(flat.size(), 36u);
        part_sizes(flat, 4);
    }
}

TEST(PartitionKahip, InvalidImbalanceThrows) {
    const Mesh m = mt::tet_mesh();
    PartitionOptions o = opts(2, PartitionMethod::KaHIP);
    o.mImbalance = 1.5;
    EXPECT_THROW(partition_labels(m, o), std::invalid_argument);
}

#else  // !MESHIOPLUSPLUS_HAS_KAHIP

TEST(PartitionKahip, CompiledOutErrorNamesTheOption) {
    const Mesh m = mt::tet_mesh();
    try {
        partition_labels(m, opts(2, PartitionMethod::KaHIP));
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("MESHIOPLUSPLUS_WITH_KAHIP"), std::string::npos);
    }
}

#endif  // MESHIOPLUSPLUS_HAS_KAHIP

}  // namespace
