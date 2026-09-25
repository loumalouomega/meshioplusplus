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

// Tecplot SZL (.szplt) through TecIO: synthetic files written with TecIO's own
// writer must read as the same data in ASCII does.
#ifdef MESHIOPLUSPLUS_HAS_TECIO

// System includes
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

#include "TECIO.h"

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/szplt.hpp"
#include "meshioplusplus/formats/tecplot.hpp"
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/registry.hpp"

namespace {

using meshioplusplus::Mesh;
using meshioplusplus::ReadError;
using meshioplusplus::ReadOptions;

// Two transient quad zones (strand 1, t = 0 and 1): the second shares X, Y and
// the connectivity of the first; P nodal, Q cell-centred.
void szplt_write(const std::string& rPath) {
    void* h = nullptr;
    ASSERT_EQ(tecFileWriterOpen(rPath.c_str(), "t", "X,Y,P,Q", 1, 0, 2, nullptr, &h), 0);
    std::vector<int32_t> types(4, 2), loc = {1, 1, 1, 0}, passive(4, 0);
    for (int s = 0; s < 2; ++s) {
        std::vector<int32_t> share(4, 0);
        if (s == 1)
            share[0] = share[1] = 1;
        int32_t zone = 0;
        ASSERT_EQ(
            tecZoneCreateFE(h, s == 0 ? "first" : "second", 3, 6, 2, types.data(), share.data(),
                            loc.data(), passive.data(), s == 1 ? 1 : 0, 0, 0, &zone),
            0);
        ASSERT_EQ(tecZoneSetUnsteadyOptions(h, zone, static_cast<double>(s), 1), 0);
        if (s == 0) {
            const std::vector<double> x = {0, 1, 2, 0, 1, 2}, y = {0, 0, 0, 1, 1, 1};
            ASSERT_EQ(tecZoneVarWriteDoubleValues(h, zone, 1, 0, 6, x.data()), 0);
            ASSERT_EQ(tecZoneVarWriteDoubleValues(h, zone, 2, 0, 6, y.data()), 0);
        }
        std::vector<double> p;
        for (int n = 0; n < 6; ++n)
            p.push_back(n + 10.0 * s);
        const std::vector<double> q = {1.0 + s, 2.0 + s};
        ASSERT_EQ(tecZoneVarWriteDoubleValues(h, zone, 3, 0, 6, p.data()), 0);
        ASSERT_EQ(tecZoneVarWriteDoubleValues(h, zone, 4, 0, 2, q.data()), 0);
        if (s == 0) {
            const std::vector<int32_t> nodes = {1, 2, 5, 4, 2, 3, 6, 5};
            ASSERT_EQ(tecZoneNodeMapWrite32(h, zone, 0, 1, 8, nodes.data()), 0);
        }
    }
    ASSERT_EQ(tecFileWriterClose(&h), 0);
}

const char* const kAscii = R"(VARIABLES = "X" "Y" "P" "Q"
ZONE T="first", ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK, NODES=6, ELEMENTS=2, VARLOCATION=([4]=CELLCENTERED), SOLUTIONTIME=0, STRANDID=1
0 1 2 0 1 2
0 0 0 1 1 1
0 1 2 3 4 5
1 2
1 2 5 4
2 3 6 5
ZONE T="second", ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK, NODES=6, ELEMENTS=2, VARLOCATION=([4]=CELLCENTERED), VARSHARELIST=([1,2]=1), CONNECTIVITYSHAREZONE=1, SOLUTIONTIME=1, STRANDID=1
10 11 12 13 14 15
2 3
)";

TEST(Szplt, ReadsAsTheSameDataInAscii) {
    const std::string szl = mt::temp_path(".szplt");
    const std::string dat = mt::temp_path(".dat");
    szplt_write(szl);
    {
        std::ofstream out(dat);
        out << kAscii;
    }
    EXPECT_EQ(meshioplusplus::szplt_time_values(szl), (std::vector<double>{0.0, 1.0}));
    for (int step = 0; step < 2; ++step) {
        ReadOptions opts;
        opts.mTimeStep = step;
        const Mesh a = meshioplusplus::read_szplt(szl, opts);
        const Mesh b = meshioplusplus::read_tecplot(dat, opts);
        mt::expect_mesh_eq(a, b, 0.0);
        const auto* pa = a.PointData("P").As<double>();
        const auto* pb = b.PointData("P").As<double>();
        for (std::size_t i = 0; i < 6; ++i)
            EXPECT_EQ(pa[i], pb[i]);
        EXPECT_EQ(a.CellData("Q", 0).As<double>()[1], 2.0 + step);
    }
    ReadOptions opts;
    opts.mTimeStep = 2;
    EXPECT_THROW(meshioplusplus::read_szplt(szl, opts), ReadError);
    EXPECT_EQ(meshioplusplus::sniff_format(szl), "szplt");
    EXPECT_EQ(meshioplusplus::resolve_format(szl, ""), "szplt");
    std::filesystem::remove(szl);
    std::filesystem::remove(dat);
}

TEST(Szplt, MissingFileIsAReadError) {
    EXPECT_THROW(meshioplusplus::read_szplt("/nonexistent/x.szplt"), ReadError);
}

}  // namespace

#endif  // MESHIOPLUSPLUS_HAS_TECIO
