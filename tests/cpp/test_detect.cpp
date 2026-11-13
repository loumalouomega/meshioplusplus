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
// detect_format: the ranked format report (roadmap 2.4.9).

#include <cstdio>
#include <fstream>
#include <string>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/operations/detect.hpp"
#include "meshioplusplus/registry.hpp"


TEST(Detect, ASignatureIsCertain) {
    const std::string p = mt::temp_path("_detect.vtu");
    meshioplusplus::registry_writers().at("vtu")(p, mt::tri_mesh());
    const auto found = meshioplusplus::detect_format(p);
    ASSERT_FALSE(found.empty());
    EXPECT_EQ(found[0].mFormat, "vtu");
    EXPECT_EQ(found[0].mConfidence, "certain");
    EXPECT_EQ(meshioplusplus::detect_format_text(p).rfind("vtu\tcertain\tcontent signature\n", 0),
              0u);
    std::remove(p.c_str());
}

TEST(Detect, ASharedExtensionListsEveryCandidate) {
    const std::string p = mt::temp_path("_detect.msh");
    meshioplusplus::registry_writers().at("gmsh")(p, mt::tri_mesh());
    const auto found = meshioplusplus::detect_format(p);
    ASSERT_GE(found.size(), 2u);
    EXPECT_EQ(found[0].mFormat, "gmsh");
    EXPECT_EQ(found[0].mConfidence, "certain");
    for (std::size_t i = 1; i < found.size(); ++i)
        EXPECT_NE(found[i].mConfidence, "certain");
    std::remove(p.c_str());
}

TEST(Detect, AWrongExtensionIsAGuessAndNothingIsEmpty) {
    const std::string stl = mt::temp_path("_detect.stl");
    { std::ofstream(stl) << "not a mesh at all\n"; }
    const auto found = meshioplusplus::detect_format(stl);
    ASSERT_FALSE(found.empty());
    for (const auto& c : found)
        EXPECT_EQ(c.mConfidence, "guess");
    const std::string none = mt::temp_path("_detect.unknownext");
    { std::ofstream(none) << "?\n"; }
    EXPECT_TRUE(meshioplusplus::detect_format(none).empty());
    EXPECT_EQ(meshioplusplus::detect_format_text(none), "");
    std::remove(stl.c_str());
    std::remove(none.c_str());
}
