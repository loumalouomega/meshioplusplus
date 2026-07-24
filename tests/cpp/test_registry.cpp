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

// External includes
#include <gtest/gtest.h>

// System includes
#include <string>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/registry.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

// The registry (registry.hpp/registry.cpp) is the shared dispatch surface behind
// the flat WASM and C/Fortran bindings. It is only exercised indirectly by the
// per-format tests, so these cover its own logic directly.

TEST(Registry, ResolveFormatUsesExtensionDefault) {
    EXPECT_EQ(resolve_format("mesh.vtu", ""), "vtu");
    EXPECT_EQ(resolve_format("mesh.stl", ""), "stl");
    // Ambiguous extensions get the repo's import-order default.
    EXPECT_EQ(resolve_format("mesh.msh", ""), "gmsh");
    EXPECT_EQ(resolve_format("mesh.inp", ""), "abaqus");
}

TEST(Registry, ResolveFormatNodeElePairStayTetgen) {
    // .node/.ele are tetgen-first everywhere; only .poly defaults to triangle.
    EXPECT_EQ(resolve_format("mesh.node", ""), "tetgen");
    EXPECT_EQ(resolve_format("mesh.ele", ""), "tetgen");
    EXPECT_EQ(resolve_format("mesh.poly", ""), "triangle");
}

TEST(Registry, ResolveFormatExplicitOverridesExtension) {
    // An explicit format wins even when it disagrees with the extension.
    EXPECT_EQ(resolve_format("mesh.vtu", "gmsh"), "gmsh");
    EXPECT_EQ(resolve_format("mesh.node", "triangle"), "triangle");
    // ...and even when the extension is unknown.
    EXPECT_EQ(resolve_format("mesh.unknownext", "vtu"), "vtu");
}

TEST(Registry, ResolveFormatUnknownExtensionThrows) {
    EXPECT_THROW(resolve_format("mesh.nosuchext", ""), ReadError);
    EXPECT_THROW(resolve_format("noextension", ""), ReadError);
}

TEST(Registry, ReadersAndWritersContainCoreFormats) {
    const auto& readers = registry_readers();
    const auto& writers = registry_writers();
    for (const char* fmt : {"vtu", "vtk", "gmsh", "stl", "ply", "obj"}) {
        EXPECT_EQ(readers.count(fmt), 1u) << "missing reader: " << fmt;
        EXPECT_EQ(writers.count(fmt), 1u) << "missing writer: " << fmt;
    }
}

TEST(Registry, OpenfoamIsReadOnly) {
    // openfoam is read-only in the C++ core: reader present, writer absent.
    EXPECT_EQ(registry_readers().count("openfoam"), 1u);
    EXPECT_EQ(registry_writers().count("openfoam"), 0u);
}

TEST(Registry, Gmsh22IsWriteOnly) {
    // Reading auto-detects the version from the file itself, so there is no
    // separate "gmsh22" reader -- only "gmsh". Before this entry existed, no
    // flat binding (WASM/C API/Fortran) could select the 2.2 writer at all,
    // which is the only one that round-trips named Cell region MEMBERSHIP
    // (4.1 keeps only the name; see write_gmsh22's doc comment).
    EXPECT_EQ(registry_readers().count("gmsh22"), 0u);
    EXPECT_EQ(registry_writers().count("gmsh22"), 1u);
}

TEST(Registry, Gmsh22RoundTripsRegionMembershipThroughTheRegistryDispatch) {
    // Exercises the registry's own tables end to end -- the same tables the
    // flat bindings (WASM/C API/Fortran) dispatch through -- so this proves
    // the wiring is real, not just that write_gmsh22 itself works.
    Mesh mesh = mt::tet_mesh();
    NDArray entries = NDArray::Uninit(DType::Int64, {1});
    entries.As<std::int64_t>()[0] = 0;
    mesh.AddRegion(Region("solid", RegionKind::Cell, 3, 7, std::move(entries)));

    const std::string path = mt::temp_path(".msh");
    registry_writers().at("gmsh22")(path, mesh);
    Mesh back = registry_readers().at("gmsh")(path);

    ASSERT_EQ(back.NumRegions(), 1u);
    EXPECT_EQ(back.Region(0).mName, "solid");
    EXPECT_EQ(back.Region(0).mTag, 7);

    // The default "gmsh" (4.1) writer keeps only the NAME, not membership --
    // the documented gap that made a distinct, selectable "gmsh22" entry
    // necessary rather than merely a nicer default.
    const std::string path41 = mt::temp_path(".msh");
    registry_writers().at("gmsh")(path41, mesh);
    Mesh back41 = registry_readers().at("gmsh")(path41);
    EXPECT_EQ(back41.NumRegions(), 0u);
}

TEST(Registry, ExtensionDefaultsMapCommonExtensions) {
    const auto& ext = registry_extension_defaults();
    EXPECT_EQ(ext.at(".vtu"), "vtu");
    EXPECT_EQ(ext.at(".stl"), "stl");
    EXPECT_EQ(ext.at(".obj"), "obj");
    // Optional-dependency extensions are mapped even when compiled out.
    EXPECT_EQ(ext.count(".med"), 1u);
    EXPECT_EQ(ext.count(".e"), 1u);
}

TEST(Registry, CompiledOutReportsMissingDependency) {
    // Unknown formats are never "compiled out" -- only absent optional-dep ones.
    EXPECT_EQ(registry_compiled_out("definitely-not-a-format"), nullptr);
    EXPECT_EQ(registry_compiled_out("vtu"), nullptr);

    // HDF5-backed formats: present -> nullptr; absent -> names "HDF5". Keying off
    // the reader table keeps this correct whether or not HDF5 is in this build.
    const bool has_med = registry_readers().count("med") == 1u;
    if (has_med) {
        EXPECT_EQ(registry_compiled_out("med"), nullptr);
    } else {
        EXPECT_STREQ(registry_compiled_out("med"), "HDF5");
    }

    const bool has_exodus = registry_readers().count("exodus") == 1u;
    if (has_exodus) {
        EXPECT_EQ(registry_compiled_out("exodus"), nullptr);
    } else {
        EXPECT_STREQ(registry_compiled_out("exodus"), "netCDF");
    }
}
