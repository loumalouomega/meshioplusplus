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
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/data_manage.hpp"

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

}  // namespace
