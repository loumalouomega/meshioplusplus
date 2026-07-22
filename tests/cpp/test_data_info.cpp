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
// Tests for the read-only per-array data summary.

// System includes
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/data_info.hpp"

namespace {

using meshioplusplus::data_array_info;
using meshioplusplus::data_info;
using meshioplusplus::DataArrayInfo;
using meshioplusplus::DataInfoReport;
using meshioplusplus::DataLocation;
using meshioplusplus::Mesh;

const DataArrayInfo* find(const DataInfoReport& rReport, DataLocation loc,
                          const std::string& rName) {
    for (const DataArrayInfo& a : rReport.mArrays)
        if (a.mLocation == loc && a.mName == rName)
            return &a;
    return nullptr;
}

TEST(DataInfo, CountsEveryArray) {
    DataInfoReport r = data_info(mt::data_mesh());
    EXPECT_EQ(r.mNumPointData, 2);  // T, v
    EXPECT_EQ(r.mNumCellData, 2);   // mat, tag
    EXPECT_EQ(r.mNumFieldData, 1);  // meta
    EXPECT_EQ(r.mArrays.size(), 5u);
}

TEST(DataInfo, ArraysAreGroupedByLocationThenSortedByName) {
    DataInfoReport r = data_info(mt::data_mesh());
    ASSERT_EQ(r.mArrays.size(), 5u);
    EXPECT_EQ(r.mArrays[0].mLocation, DataLocation::Point);
    EXPECT_EQ(r.mArrays[0].mName, "T");
    EXPECT_EQ(r.mArrays[1].mName, "v");
    EXPECT_EQ(r.mArrays[2].mLocation, DataLocation::Cell);
    EXPECT_EQ(r.mArrays[2].mName, "mat");
    EXPECT_EQ(r.mArrays[3].mName, "tag");
    EXPECT_EQ(r.mArrays[4].mLocation, DataLocation::Field);
    EXPECT_EQ(r.mArrays[4].mName, "meta");
}

TEST(DataInfo, ScalarPointArrayStatistics) {
    DataInfoReport r = data_info(mt::data_mesh());
    const DataArrayInfo* a = find(r, DataLocation::Point, "T");
    ASSERT_NE(a, nullptr);
    // T = {0, 1, 11, 10, 2, 12}
    EXPECT_EQ(a->mNumEntries, 6);
    EXPECT_EQ(a->mNumComponents, 1);
    EXPECT_EQ(a->mNumValues, 6);
    EXPECT_DOUBLE_EQ(a->mMin, 0.0);
    EXPECT_DOUBLE_EQ(a->mMax, 12.0);
    EXPECT_NEAR(a->mMean, 36.0 / 6.0, 1e-12);
    EXPECT_EQ(a->mNumNan, 0);
    EXPECT_EQ(a->mNumInf, 0);
    EXPECT_EQ(a->mNumFinite, 6);
}

TEST(DataInfo, VectorArrayReportsPerComponentStatistics) {
    DataInfoReport r = data_info(mt::data_mesh());
    const DataArrayInfo* a = find(r, DataLocation::Point, "v");
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->mNumComponents, 3);
    EXPECT_EQ(a->mNumEntries, 6);
    EXPECT_EQ(a->mNumValues, 18);
    ASSERT_EQ(a->mMinPerComponent.size(), 3u);
    // Component 0 of v is {1,0,0,1,2,0}
    EXPECT_DOUBLE_EQ(a->mMinPerComponent[0], 0.0);
    EXPECT_DOUBLE_EQ(a->mMaxPerComponent[0], 2.0);
    EXPECT_NEAR(a->mMeanPerComponent[0], 4.0 / 6.0, 1e-12);
    // Component 2 is {0,0,1,0,0,0}
    EXPECT_DOUBLE_EQ(a->mMaxPerComponent[2], 1.0);
}

TEST(DataInfo, CellDataSpansEveryBlock) {
    DataInfoReport r = data_info(mt::data_mesh());
    const DataArrayInfo* a = find(r, DataLocation::Cell, "mat");
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->mNumBlocks, 2);
    // 2 triangle cells + 1 quad cell
    EXPECT_EQ(a->mNumEntries, 3);
    EXPECT_EQ(a->mNumValues, 3);
    EXPECT_DOUBLE_EQ(a->mMin, 1.0);
    EXPECT_DOUBLE_EQ(a->mMax, 3.0);
    EXPECT_NEAR(a->mMean, 2.0, 1e-12);
    EXPECT_FALSE(a->mInconsistentBlocks);
}

TEST(DataInfo, ReportsDtypeAsStored) {
    DataInfoReport r = data_info(mt::data_mesh());
    const DataArrayInfo* t = find(r, DataLocation::Point, "T");
    ASSERT_NE(t, nullptr);
    EXPECT_TRUE(meshioplusplus::detail::is_float_dtype(t->mDtype));
    const DataArrayInfo* tag = find(r, DataLocation::Cell, "tag");
    ASSERT_NE(tag, nullptr);
    EXPECT_FALSE(meshioplusplus::detail::is_float_dtype(tag->mDtype));
}

TEST(DataInfo, CountsNaNAndInf) {
    Mesh m = mt::data_mesh();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    m.AddPointData("n", mt::data_array({1.0, nan, inf, -inf, 2.0, nan}));
    DataInfoReport r = data_info(m);
    const DataArrayInfo* a = find(r, DataLocation::Point, "n");
    ASSERT_NE(a, nullptr);
    EXPECT_EQ(a->mNumNan, 2);
    EXPECT_EQ(a->mNumInf, 2);
    EXPECT_EQ(a->mNumFinite, 2);
    EXPECT_EQ(a->mNumValues, 6);
    // Reductions ignore the non-finite values entirely.
    EXPECT_DOUBLE_EQ(a->mMin, 1.0);
    EXPECT_DOUBLE_EQ(a->mMax, 2.0);
    EXPECT_NEAR(a->mMean, 1.5, 1e-12);
}

TEST(DataInfo, AllNonFiniteArrayYieldsNaNStatistics) {
    Mesh m = mt::data_mesh();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    m.AddPointData("n", mt::data_array({nan, nan, nan, nan, nan, nan}));
    const DataArrayInfo a = data_array_info(m, DataLocation::Point, "n");
    EXPECT_EQ(a.mNumFinite, 0);
    EXPECT_TRUE(std::isnan(a.mMin));
    EXPECT_TRUE(std::isnan(a.mMax));
    EXPECT_TRUE(std::isnan(a.mMean));
}

TEST(DataInfo, EmptyMeshHasNoArrays) {
    Mesh m;
    DataInfoReport r = data_info(m);
    EXPECT_EQ(r.mArrays.size(), 0u);
    EXPECT_EQ(r.mNumPointData, 0);
    EXPECT_EQ(r.mNumCellData, 0);
    EXPECT_EQ(r.mNumFieldData, 0);
}

TEST(DataInfo, MeshWithGeometryButNoDataHasNoArrays) {
    DataInfoReport r = data_info(mt::tri_mesh());
    EXPECT_EQ(r.mArrays.size(), 0u);
}

TEST(DataInfo, SingleArrayLookup) {
    const DataArrayInfo a = data_array_info(mt::data_mesh(), DataLocation::Field, "meta");
    EXPECT_EQ(a.mName, "meta");
    EXPECT_EQ(a.mNumEntries, 3);
    EXPECT_DOUBLE_EQ(a.mMin, 1.0);
    EXPECT_DOUBLE_EQ(a.mMax, 3.0);
}

TEST(DataInfo, UnknownNameThrows) {
    EXPECT_THROW(data_array_info(mt::data_mesh(), DataLocation::Point, "nope"),
                 std::invalid_argument);
}

TEST(DataInfo, DoesNotModifyTheMesh) {
    Mesh in = mt::data_mesh();
    const std::size_t before = in.NumPointData();
    data_info(in);
    EXPECT_EQ(in.NumPointData(), before);
    EXPECT_TRUE(in.HasPointData("T"));
}

}  // namespace
