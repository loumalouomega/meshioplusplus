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
//                   SPDX-License-Identifier: MIT
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
// Small versions of the thicken and hollow C API test and of the decimate deviation test
// (`test_c_api.cpp`, `test_decimate.cpp`). The full ones offset a sphere of 5,120 triangles
// and decimate one of 20,480 to a deviation of 0.005, which takes minutes each in the
// instrumented coverage build, so a pull request skips them there; these run the same calls
// on spheres of 320 and 1,280 triangles and check the same properties. Their suite name
// starts with `Quick`, which the job's exclusion does not match.

// System includes
#include <cmath>
#include <cstddef>
#include <cstdint>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/meshioplusplus.h"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/hausdorff.hpp"
#include "meshioplusplus/operations/primitives.hpp"

// The offset surface, the hollow wall and the thickened sphere through the C API.
TEST(QuickCApi, OffsetThickenAndHollowASphere) {
    mio_gen_opts g;
    mio_gen_opts_init(&g);
    g.counts[0] = 2;
    mio_mesh* s = mio_icosa_sphere(&g);
    ASSERT_NE(s, nullptr);
    mio_offset_report r;
    mio_mesh* grown = mio_offset_surface(s, 0.1, 0, &r);
    ASSERT_NE(grown, nullptr) << mio_last_error();
    EXPECT_EQ(r.num_folded, 0);
    EXPECT_GE(r.max_move, 0.1);
    EXPECT_EQ(mio_offset_warnings(nullptr, 0), 0);
    mio_mesh* shell = mio_hollow(s, 0.1, 0, &r);
    ASSERT_NE(shell, nullptr) << mio_last_error();
    EXPECT_EQ(r.num_rim_quads, 0);
    mio_mesh* thick = mio_thicken(s, 0.1, MIO_THICKEN_OUTWARD, 0, &r);
    ASSERT_NE(thick, nullptr) << mio_last_error();
    EXPECT_EQ(r.num_rim_quads, 0);
    EXPECT_EQ(mio_thicken(s, 0.1, 3, 0, &r), nullptr);
    EXPECT_EQ(mio_offset_surface(s, std::nan(""), 0, &r), nullptr);
    mio_mesh_free(thick);
    mio_mesh_free(shell);
    mio_mesh_free(grown);
    mio_mesh_free(s);
}

// Decimating to a deviation bound keeps it, halves the faces at least, is deterministic, and
// a tighter bound keeps more faces.
TEST(QuickDecimate, MaxDeviationHoldsTwoSidedOnTheCoverage) {
    using namespace meshioplusplus;
    const Mesh s = icosa_sphere(1.0, 3);
    DecimateOptions o;
    o.mMaxDeviation = 0.03;
    const DecimateResult r = decimate(s, o);
    EXPECT_TRUE(r.mDeviationMet);
    EXPECT_LE(r.mDeviation, 0.03);
    EXPECT_GT(r.mRounds, 1);
    EXPECT_LT(r.mMesh.Cells(0).NumCells(), s.Cells(0).NumCells() / 2);
    HausdorffOptions ho;
    ho.mFaceSamples = 2;
    EXPECT_LT(hausdorff_distance(r.mMesh, s, ho).mDistance, 0.07);
    const DecimateResult again = decimate(s, o);
    ASSERT_EQ(again.mMesh.Cells(0).NumCells(), r.mMesh.Cells(0).NumCells());
    const double* a = r.mMesh.Points().As<double>();
    const double* b = again.mMesh.Points().As<double>();
    for (std::size_t i = 0; i < 3 * r.mMesh.NumPoints(); ++i)
        ASSERT_EQ(a[i], b[i]);
    o.mMaxDeviation = 0.008;
    EXPECT_GT(decimate(s, o).mMesh.Cells(0).NumCells(), r.mMesh.Cells(0).NumCells());
}
