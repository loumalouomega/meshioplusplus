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
// Tests for the data array-management operation (rename / drop / keep).

// System includes
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/data_manage.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace {

using meshioplusplus::DataKey;
using meshioplusplus::DataLocation;
using meshioplusplus::DataManageOptions;
using meshioplusplus::DataManageResult;
using meshioplusplus::DataRename;
using meshioplusplus::Mesh;

TEST(DataManage, RenamePreservesValuesAndGeometry) {
    Mesh in = mt::data_mesh();
    Mesh out = meshioplusplus::data_rename(in, DataLocation::Point, "T", "temperature");

    EXPECT_FALSE(out.HasPointData("T"));
    ASSERT_TRUE(out.HasPointData("temperature"));
    ASSERT_EQ(out.PointData("temperature").Size(), in.PointData("T").Size());
    for (std::size_t i = 0; i < in.PointData("T").Size(); ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out.PointData("temperature"), i),
                         meshioplusplus::detail::read_double(in.PointData("T"), i));
    // Untouched arrays survive.
    EXPECT_TRUE(out.HasPointData("v"));
    EXPECT_TRUE(out.HasCellData("mat"));
    EXPECT_TRUE(out.HasFieldData("meta"));
    mt::expect_same_geometry(in, out);
}

TEST(DataManage, RenameMultiBlockCellData) {
    Mesh in = mt::data_mesh();
    Mesh out = meshioplusplus::data_rename(in, DataLocation::Cell, "mat", "material");
    ASSERT_TRUE(out.HasCellData("material"));
    EXPECT_FALSE(out.HasCellData("mat"));
    // Every block must survive the rename — the uniform API requires exactly
    // one cell_data array per cell block.
    ASSERT_EQ(out.CellDataNumBlocks("material"), out.NumCellBlocks());
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out.CellData("material", 0), 0), 1.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(out.CellData("material", 1), 0), 3.0);
}

TEST(DataManage, DropRemovesExactlyTheNamedArrays) {
    Mesh in = mt::data_mesh();
    Mesh out = meshioplusplus::data_drop(in, DataLocation::Point, {"T"});
    EXPECT_FALSE(out.HasPointData("T"));
    EXPECT_TRUE(out.HasPointData("v"));
    // Other locations untouched.
    EXPECT_TRUE(out.HasCellData("mat"));
    EXPECT_TRUE(out.HasCellData("tag"));
    EXPECT_TRUE(out.HasFieldData("meta"));
    mt::expect_same_geometry(in, out);
}

TEST(DataManage, KeepRetainsOnlyTheNamedSubset) {
    Mesh in = mt::data_mesh();
    Mesh out = meshioplusplus::data_keep(in, DataLocation::Cell, {"tag"});
    EXPECT_TRUE(out.HasCellData("tag"));
    EXPECT_FALSE(out.HasCellData("mat"));
    // A location the whitelist does not mention is left completely alone.
    EXPECT_TRUE(out.HasPointData("T"));
    EXPECT_TRUE(out.HasPointData("v"));
    EXPECT_TRUE(out.HasFieldData("meta"));
}

TEST(DataManage, KeepNothingDropsEverythingAtThatLocation) {
    Mesh in = mt::data_mesh();
    Mesh out = meshioplusplus::data_keep(in, DataLocation::Point, {});
    EXPECT_EQ(out.NumPointData(), 0u);
    EXPECT_TRUE(out.HasCellData("mat"));
}

TEST(DataManage, UnknownKeyListsAvailableKeys) {
    Mesh in = mt::data_mesh();
    try {
        meshioplusplus::data_drop(in, DataLocation::Point, {"nope"});
        FAIL() << "expected an exception";
    } catch (const std::invalid_argument& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("nope"), std::string::npos);
        EXPECT_NE(msg.find("point_data"), std::string::npos);
        // Every available key is named.
        EXPECT_NE(msg.find("T"), std::string::npos);
        EXPECT_NE(msg.find("v"), std::string::npos);
    }
}

TEST(DataManage, IgnoreMissingSilencesUnknownKeys) {
    Mesh in = mt::data_mesh();
    Mesh out = meshioplusplus::data_drop(in, DataLocation::Point, {"nope"}, /*IgnoreMissing=*/true);
    EXPECT_TRUE(out.HasPointData("T"));
}

TEST(DataManage, RenameOntoExistingNameThrows) {
    Mesh in = mt::data_mesh();
    EXPECT_THROW(meshioplusplus::data_rename(in, DataLocation::Point, "T", "v"),
                 std::invalid_argument);
}

TEST(DataManage, TwoRenamesToTheSameTargetThrow) {
    Mesh in = mt::data_mesh();
    DataManageOptions opts;
    opts.rename.push_back(DataRename{DataLocation::Point, "T", "x"});
    opts.rename.push_back(DataRename{DataLocation::Point, "v", "x"});
    EXPECT_THROW(meshioplusplus::data_manage(in, opts), std::invalid_argument);
}

TEST(DataManage, SwapNamesIsAllowed) {
    // T -> v and v -> T: neither target "already exists" in the result, because
    // both are renamed away in the same pass.
    Mesh in = mt::data_mesh();
    DataManageOptions opts;
    opts.rename.push_back(DataRename{DataLocation::Point, "T", "v"});
    opts.rename.push_back(DataRename{DataLocation::Point, "v", "T"});
    DataManageResult r = meshioplusplus::data_manage(in, opts);
    ASSERT_TRUE(r.mMesh.HasPointData("T"));
    ASSERT_TRUE(r.mMesh.HasPointData("v"));
    // The old "v" (a 3-vector) is now called "T".
    EXPECT_EQ(r.mMesh.PointData("T").Size(), in.PointData("v").Size());
    EXPECT_EQ(r.mMesh.PointData("v").Size(), in.PointData("T").Size());
}

TEST(DataManage, ReportsWhatWasDroppedAndRenamed) {
    Mesh in = mt::data_mesh();
    DataManageOptions opts;
    opts.drop.push_back(DataKey{DataLocation::Point, "T"});
    opts.rename.push_back(DataRename{DataLocation::Field, "meta", "metadata"});
    DataManageResult r = meshioplusplus::data_manage(in, opts);
    ASSERT_EQ(r.mDropped.size(), 1u);
    EXPECT_EQ(r.mDropped[0], "point_data:T");
    ASSERT_EQ(r.mRenamed.size(), 1u);
    EXPECT_EQ(r.mRenamed[0].first, "field_data:meta");
    EXPECT_EQ(r.mRenamed[0].second, "field_data:metadata");
}

TEST(DataManage, KeepThenDropThenRenameOrder) {
    Mesh in = mt::data_mesh();
    DataManageOptions opts;
    opts.keep.push_back(DataKey{DataLocation::Cell, "mat"});
    opts.keep.push_back(DataKey{DataLocation::Cell, "tag"});
    opts.drop.push_back(DataKey{DataLocation::Cell, "tag"});
    opts.rename.push_back(DataRename{DataLocation::Cell, "mat", "material"});
    DataManageResult r = meshioplusplus::data_manage(in, opts);
    EXPECT_TRUE(r.mMesh.HasCellData("material"));
    EXPECT_FALSE(r.mMesh.HasCellData("tag"));
    EXPECT_FALSE(r.mMesh.HasCellData("mat"));
}

TEST(DataManage, SetsToDataOverlapsOrderingAndUncovered) {
    using meshioplusplus::Region;
    using meshioplusplus::RegionKind;
    Mesh in = mt::data_mesh();
    in.AddRegion(Region("alpha", RegionKind::Point, mt::int_data_array({0, 1})));
    in.AddRegion(Region("zeta", RegionKind::Point, mt::int_data_array({1, 2})));
    in.AddRegion(Region("empty", RegionKind::Point, mt::int_data_array({})));
    in.AddRegion(Region("wall", RegionKind::Side, 2, 9,
                        meshioplusplus::NDArray(meshioplusplus::DType::Int64, {0, 2})));
    Mesh out = meshioplusplus::sets_to_data(in, DataLocation::Point, std::nullopt, ":",
                                            {"zeta", "empty", "alpha"});
    const auto& data = out.PointData("zeta:empty:alpha");
    EXPECT_EQ(data.Dtype(), meshioplusplus::DType::Int64);
    const auto* labels = data.As<std::int64_t>();
    EXPECT_EQ(labels[0], 2);
    EXPECT_EQ(labels[1], 2);
    EXPECT_EQ(labels[2], 0);
    EXPECT_EQ(labels[3], -1);
    ASSERT_EQ(out.NumRegions(), 1u);
    EXPECT_EQ(out.Region(0).mName, "wall");
    EXPECT_EQ(out.Region(0).mTag, 9);
    EXPECT_EQ(in.NumRegions(), 4u);
    EXPECT_TRUE(out.HasPointData("T"));
    EXPECT_TRUE(out.HasFieldData("meta"));
    mt::expect_same_geometry(in, out);
    EXPECT_THROW(meshioplusplus::sets_to_data(in, DataLocation::Field), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::sets_to_data(in, DataLocation::Point, "bad", "-",
                                              {"alpha", "alpha", "empty"}),
                 std::invalid_argument);
}

TEST(DataManage, CellSetsUseGlobalIndicesAndRetainBlockOrder) {
    using meshioplusplus::Region;
    using meshioplusplus::RegionKind;
    Mesh in = mt::data_mesh();
    in.AddRegion(Region("first", RegionKind::Cell, mt::int_data_array({0, 2})));
    in.AddRegion(Region("second", RegionKind::Cell, mt::int_data_array({1, 2})));
    Mesh out = meshioplusplus::sets_to_data(in, DataLocation::Cell);
    ASSERT_EQ(out.CellDataNumBlocks("first-second"), 2u);
    EXPECT_EQ(out.CellData("first-second", 0).As<std::int64_t>()[0], 0);
    EXPECT_EQ(out.CellData("first-second", 0).As<std::int64_t>()[1], 1);
    EXPECT_EQ(out.CellData("first-second", 1).As<std::int64_t>()[0], 1);
    Mesh back = meshioplusplus::data_to_sets(out, DataLocation::Cell, "first-second");
    ASSERT_EQ(back.NumRegions(), 2u);
    EXPECT_EQ(back.Region(0).mName, "first");
    ASSERT_EQ(back.Region(1).NumEntries(), 2u);
    EXPECT_EQ(back.Region(1).Entries()[0], 1);
    EXPECT_EQ(back.Region(1).Entries()[1], 2);
    EXPECT_FALSE(back.HasCellData("first-second"));
    EXPECT_TRUE(back.HasCellData("tag"));
    mt::expect_same_geometry(in, back);
}

TEST(DataManage, DataToSetsKeepsExactIntegerExtremaAndNames) {
    using meshioplusplus::DType;
    using meshioplusplus::NDArray;
    Mesh in = mt::data_mesh();
#if defined(MESHIOPLUSPLUS_MESH_BACKEND_MESHIO)
    // The other backends canonicalize all integer storage to signed Int64.
    NDArray values(DType::UInt64, {in.NumPoints()});
    std::fill_n(values.As<std::uint64_t>(), in.NumPoints(),
                std::numeric_limits<std::uint64_t>::max());
    values.As<std::uint64_t>()[1] = 0;
    values.As<std::uint64_t>()[2] = 9007199254740993ULL;
    in.AddPointData("ids", std::move(values));
    auto out = meshioplusplus::data_to_sets(in, DataLocation::Point, "ids");
    ASSERT_EQ(out.NumRegions(), 3u);
    bool found_max = false, found_precise = false;
    for (std::size_t i = 0; i < out.NumRegions(); ++i) {
        const auto& region = out.Region(i);
        found_max |= region.mName == "set-key-18446744073709551615";
        found_precise |= region.mName == "set-key-9007199254740993";
    }
    EXPECT_TRUE(found_max);
    EXPECT_TRUE(found_precise);
#endif
    NDArray signed_values(DType::Int64, {in.NumPoints()});
    std::fill_n(signed_values.As<std::int64_t>(), in.NumPoints(),
                std::numeric_limits<std::int64_t>::min());
    in.AddPointData("signed", std::move(signed_values));
    auto signed_out = meshioplusplus::data_to_sets(in, DataLocation::Point, "signed");
    // One tag uses the name obtained by splitting the key.
    ASSERT_EQ(signed_out.NumRegions(), 1u);
    EXPECT_EQ(signed_out.Region(0).mName, "signed");
    EXPECT_EQ(signed_out.Region(0).NumEntries(), in.NumPoints());
}

TEST(DataManage, DataToSetsRetainsExistingRegionIdentity) {
    using meshioplusplus::Region;
    using meshioplusplus::RegionKind;
    Mesh in = mt::data_mesh();
    in.AddRegion(Region("left", RegionKind::Point, 0, 17, mt::int_data_array({5})));
    in.AddPointData("left-left-right", mt::int_data_array({-1, 4, -1, 4, -1, 4}));
    auto out = meshioplusplus::data_to_sets(in, DataLocation::Point, "left-left-right");
    ASSERT_EQ(out.NumRegions(), 2u);
    EXPECT_EQ(out.Region(0).mName, "left");
    EXPECT_EQ(out.Region(0).mTag, 17);
    EXPECT_EQ(out.Region(0).mDim, 0);
    EXPECT_EQ(out.Region(0).NumEntries(), 3u);
    EXPECT_EQ(out.Region(1).mName, "right");
}

TEST(DataManage, SetsDataEmptyAndInvalidInputs) {
    Mesh in = mt::data_mesh();
    auto copy = meshioplusplus::sets_to_data(in, DataLocation::Point);
    EXPECT_EQ(copy.NumPointData(), in.NumPointData());
    EXPECT_NE(copy.Points().Data(), in.Points().Data());
    EXPECT_THROW(meshioplusplus::data_to_sets(in, DataLocation::Point, "T"), std::invalid_argument);
    EXPECT_THROW(meshioplusplus::data_to_sets(in, DataLocation::Point, "v"),
                 meshioplusplus::Unsupported);
    EXPECT_THROW(meshioplusplus::data_to_sets(in, DataLocation::Cell, "missing"),
                 std::invalid_argument);
    EXPECT_THROW(meshioplusplus::data_to_sets(in, DataLocation::Field, "meta"),
                 std::invalid_argument);
    Mesh empty;
    empty.AssignPoints(meshioplusplus::NDArray(meshioplusplus::DType::Float64, {0, 3}));
    empty.AddPointData("ids", mt::int_data_array({}));
    auto out = meshioplusplus::data_to_sets(empty, DataLocation::Point, "ids");
    EXPECT_FALSE(out.HasPointData("ids"));
    EXPECT_EQ(out.NumRegions(), 0u);
}

TEST(DataManage, SetsDataKeepsEmptyAndRaggedCellBlocks) {
    using meshioplusplus::Region;
    using meshioplusplus::RegionKind;
    Mesh in = mt::quad_mesh();
    in.AddCellBlock("triangle", meshioplusplus::NDArray(meshioplusplus::DType::Int64, {0, 3}));
    in.AddPolygonBlock("polygon", std::vector<std::vector<std::int64_t>>{{0, 1, 2}, {0, 1, 2, 3}});
    in.AddRegion(Region("set", RegionKind::Cell, mt::int_data_array({0, 2})));
    auto out = meshioplusplus::sets_to_data(in, DataLocation::Cell);
    ASSERT_EQ(out.CellDataNumBlocks("set"), 3u);
    EXPECT_EQ(out.CellData("set", 1).Size(), 0u);
    EXPECT_EQ(out.CellData("set", 2).As<std::int64_t>()[0], -1);
    EXPECT_EQ(out.CellData("set", 2).As<std::int64_t>()[1], 0);
    for (std::size_t row = 0; row < in.Cells(2).NumCells(); ++row) {
        const auto a = in.Cells(2).Row(row), b = out.Cells(2).Row(row);
        ASSERT_EQ(in.Cells(2).RowSize(row), out.Cells(2).RowSize(row));
        for (std::size_t i = 0; i < in.Cells(2).RowSize(row); ++i)
            EXPECT_EQ(a[i], b[i]);
    }
    mt::expect_same_geometry(in, out);
}

}  // namespace
