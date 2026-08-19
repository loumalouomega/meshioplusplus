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

// System includes
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/data_integrate.hpp"
#include "meshioplusplus/region.hpp"

namespace {

using meshioplusplus::data_integrate;
using meshioplusplus::DataIntegrateOptions;
using meshioplusplus::DataIntegrateReport;
using meshioplusplus::DType;
using meshioplusplus::FieldIntegralArray;
using meshioplusplus::FieldIntegralRegion;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;
using meshioplusplus::Region;
using meshioplusplus::RegionKind;

NDArray i64(const std::vector<std::int64_t>& rVals) {
    NDArray a = NDArray::Uninit(DType::Int64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<std::int64_t>()[i] = rVals[i];
    return a;
}

// A 1 x n row of unit quads in the z = 0 plane over [0, n] x [0, 1] -- each
// cell has area exactly 1, so a total is trivial to hand-verify.
Mesh quad_row(int n, double shift = 0.0) {
    std::vector<std::vector<double>> pts;
    for (int j = 0; j < 2; ++j)
        for (int i = 0; i <= n; ++i)
            pts.push_back({static_cast<double>(i) + shift, static_cast<double>(j) + shift, 0.0});
    auto pid = [n](int i, int j) { return static_cast<std::int64_t>(j * (n + 1) + i); };
    std::vector<std::vector<std::int64_t>> cells;
    for (int i = 0; i < n; ++i)
        cells.push_back({pid(i, 0), pid(i + 1, 0), pid(i + 1, 1), pid(i, 1)});
    return mt::make_mesh(std::move(pts), "quad", std::move(cells));
}

const FieldIntegralRegion* find_region(const FieldIntegralArray& rArr, const std::string& rName) {
    for (const auto& r : rArr.mRegions)
        if (r.mName == rName)
            return &r;
    return nullptr;
}

// --- basic reductions --------------------------------------------------------

TEST(DataIntegrate, UniformFieldTotalEqualsValueTimesDomainMeasure) {
    Mesh m = quad_row(4);  // 4 unit-area cells
    m.AddCellData("f", {mt::data_array({3.0, 3.0, 3.0, 3.0})});

    const DataIntegrateReport report = data_integrate(m);
    ASSERT_EQ(report.mArrays.size(), 1u);
    const FieldIntegralArray& arr = report.mArrays[0];
    EXPECT_EQ(arr.mName, "f");
    EXPECT_EQ(arr.mNumComponents, 1);
    EXPECT_DOUBLE_EQ(arr.mDomain.mDomainMeasurePerComponent[0], 4.0);
    EXPECT_DOUBLE_EQ(arr.mDomain.mTotalPerComponent[0], 12.0);
}

TEST(DataIntegrate, MeanRecoversTheConstant) {
    Mesh m = quad_row(3);
    m.AddCellData("f", {mt::data_array({5.0, 5.0, 5.0})});

    const DataIntegrateReport report = data_integrate(m);
    EXPECT_DOUBLE_EQ(report.mArrays[0].mDomain.mMeanPerComponent[0], 5.0);
}

TEST(DataIntegrate, TranslatingTheMeshDoesNotChangeTheTotal) {
    Mesh m = quad_row(5);
    Mesh shifted = quad_row(5, 1.0e8);
    std::vector<double> vals = {1.0, 2.0, 3.0, 4.0, 5.0};
    m.AddCellData("f", {mt::data_array(vals)});
    shifted.AddCellData("f", {mt::data_array(vals)});

    const double total = data_integrate(m).mArrays[0].mDomain.mTotalPerComponent[0];
    const double total_shifted = data_integrate(shifted).mArrays[0].mDomain.mTotalPerComponent[0];
    EXPECT_NEAR(total, total_shifted, 1e-6);
}

// --- exclusion rules ----------------------------------------------------------

TEST(DataIntegrate, RaggedOrUnsupportedCellsAreExcludedNotZeroed) {
    // A quad row plus a "vertex" block (dimension 0, no measure at all).
    Mesh m = quad_row(3);
    m.AddCellBlock("vertex", mt::conn_from({{0}, {1}}));
    m.AddCellData("f", {mt::data_array({1.0, 1.0, 1.0}), mt::data_array({99.0, 99.0})});

    const DataIntegrateReport report = data_integrate(m);
    const FieldIntegralRegion& dom = report.mArrays[0].mDomain;
    EXPECT_EQ(dom.mNumCells, 3);
    EXPECT_EQ(dom.mNumSkipped, 2);
    EXPECT_DOUBLE_EQ(dom.mDomainMeasurePerComponent[0], 3.0);
    EXPECT_DOUBLE_EQ(dom.mTotalPerComponent[0], 3.0);
}

TEST(DataIntegrate, NonFiniteValuesAreExcludedFromNumeratorAndDenominator) {
    Mesh m = quad_row(3);
    m.AddCellData("f", {mt::data_array({1.0, std::nan(""), 3.0})});

    const DataIntegrateReport report = data_integrate(m);
    const FieldIntegralRegion& dom = report.mArrays[0].mDomain;
    EXPECT_EQ(dom.mNumNanPerComponent[0], 1);
    // The NaN cell must not shrink the total, and must not appear in the
    // denominator either -- mean over the two remaining cells is 2.0.
    EXPECT_DOUBLE_EQ(dom.mTotalPerComponent[0], 4.0);
    EXPECT_DOUBLE_EQ(dom.mDomainMeasurePerComponent[0], 2.0);
    EXPECT_DOUBLE_EQ(dom.mMeanPerComponent[0], 2.0);
}

// --- regions --------------------------------------------------------------

TEST(DataIntegrate, PerRegionTotalsAreIndependentNotAPartition) {
    Mesh m = quad_row(4);  // cells 0,1,2,3
    m.AddCellData("f", {mt::data_array({1.0, 2.0, 3.0, 4.0})});
    // Overlapping regions: cell 1 is in both "a" and "b"; cell 3 in neither.
    m.AddRegion(Region("a", RegionKind::Cell, i64({0, 1})));
    m.AddRegion(Region("b", RegionKind::Cell, i64({1, 2})));

    const DataIntegrateReport report = data_integrate(m);
    const FieldIntegralArray& arr = report.mArrays[0];
    ASSERT_EQ(arr.mRegions.size(), 2u);

    const FieldIntegralRegion* a = find_region(arr, "a");
    const FieldIntegralRegion* b = find_region(arr, "b");
    ASSERT_NE(a, nullptr);
    ASSERT_NE(b, nullptr);
    EXPECT_DOUBLE_EQ(a->mTotalPerComponent[0], 3.0);  // cells 0 + 1
    EXPECT_DOUBLE_EQ(b->mTotalPerComponent[0], 5.0);  // cells 1 + 2
    // Cell 1 contributed to both -- not consumed by either.
    EXPECT_DOUBLE_EQ(arr.mDomain.mTotalPerComponent[0], 10.0);
}

TEST(DataIntegrate, PointAndSideRegionsAreSkipped) {
    Mesh m = quad_row(2);
    m.AddCellData("f", {mt::data_array({1.0, 1.0})});
    m.AddRegion(Region("pts", RegionKind::Point, i64({0, 1})));

    const DataIntegrateReport report = data_integrate(m);
    EXPECT_TRUE(report.mArrays[0].mRegions.empty());
}

// --- validation -------------------------------------------------------------

TEST(DataIntegrate, EmptyArrayNamesMeansEveryCellDataArray) {
    Mesh m = quad_row(2);
    m.AddCellData("f", {mt::data_array({1.0, 1.0})});
    m.AddCellData("g", {mt::data_array({2.0, 2.0})});

    const DataIntegrateReport report = data_integrate(m);
    EXPECT_EQ(report.mArrays.size(), 2u);
}

TEST(DataIntegrate, PointDataNameThrowsAndNamesTheFix) {
    Mesh m = quad_row(2);
    m.AddPointData("f", mt::data_array({1.0, 2.0, 3.0, 4.0, 5.0, 6.0}));

    DataIntegrateOptions opts;
    opts.mArrayNames = {"f"};
    try {
        data_integrate(m, opts);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("point_data_to_cell_data"), std::string::npos);
    }
}

TEST(DataIntegrate, UnknownArrayNameListsWhatExists) {
    Mesh m = quad_row(2);
    m.AddCellData("f", {mt::data_array({1.0, 1.0})});

    DataIntegrateOptions opts;
    opts.mArrayNames = {"nope"};
    try {
        data_integrate(m, opts);
        FAIL() << "expected std::invalid_argument";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("f"), std::string::npos);
    }
}

// --- determinism -------------------------------------------------------------

TEST(DataIntegrate, RepeatRunsByteIdentical) {
    Mesh m = quad_row(8);
    std::vector<double> vals;
    for (int i = 0; i < 8; ++i)
        vals.push_back(1.0 + 0.37 * i);
    m.AddCellData("f", {mt::data_array(vals)});

    const double first = data_integrate(m).mArrays[0].mDomain.mTotalPerComponent[0];
    for (int i = 0; i < 5; ++i)
        EXPECT_EQ(data_integrate(m).mArrays[0].mDomain.mTotalPerComponent[0], first);
}

}  // namespace
