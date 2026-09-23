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
#include "meshioplusplus/write_options.hpp"
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
    EXPECT_EQ(resolve_format("deck.k", ""), "lsdyna");
    EXPECT_EQ(resolve_format("deck.key", ""), "lsdyna");
    EXPECT_EQ(resolve_format("deck.dyn", ""), "lsdyna");
    EXPECT_EQ(resolve_format("mesh.mail", ""), "code_aster");
    EXPECT_EQ(resolve_format("model.pat", ""), "patran");
    EXPECT_EQ(resolve_format("model.out", ""), "patran");
    EXPECT_EQ(resolve_format("mesh.mphbin", ""), "mphbin");
    EXPECT_EQ(resolve_format("mesh.mphtxt", ""), "mphtxt");
    EXPECT_EQ(resolve_format("results.frd", ""), "frd");
    // `.h5` is MSC Nastran HDF5; the longer `.post.h5` stays GiD.
    EXPECT_EQ(resolve_format("job.h5", ""), "nastran_h5");
    EXPECT_EQ(resolve_format("job.post.h5", ""), "gid");
    // glTF is one format with two suffixes; the writer tells the containers apart.
    EXPECT_EQ(resolve_format("model.glb", ""), "gltf");
    EXPECT_EQ(resolve_format("model.gltf", ""), "gltf");
}

TEST(Registry, ElmerIsReadWriteByDirectoryNotExtension) {
    EXPECT_EQ(registry_readers().count("elmer"), 1u);
    EXPECT_EQ(registry_writers().count("elmer"), 1u);
    EXPECT_EQ(registry_readers_ex().count("elmer"), 1u);
    // A directory has no extension: resolve_format refuses, sniff_format finds it.
    EXPECT_THROW(resolve_format("some/mesh_dir", ""), ReadError);
}

TEST(Registry, GltfIsWriteOnly) {
    EXPECT_EQ(registry_writers().count("gltf"), 1u);
    EXPECT_EQ(registry_readers().count("gltf"), 0u);
    // The suffix decides the container, so an encoding request is refused by name.
    WriteOptions options;
    options.mEncoding = WriteEncoding::Binary;
    Mesh m;
    EXPECT_THROW(registry_write_ex("x.glb", m, "gltf", options), WriteError);
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

TEST(Registry, OpenfoamIsReadWriteAndResolvesByExtension) {
    // openfoam was read-only until v9.20.0 -- this test asserted that.
    EXPECT_EQ(registry_readers().count("openfoam"), 1u);
    EXPECT_EQ(registry_writers().count("openfoam"), 1u);
    // `.foam` also joined the extension table then; before that
    // `resolve_format` THREW on it, so no flat binding could reach a case by
    // extension in either direction.
    EXPECT_EQ(resolve_format("case.foam", ""), "openfoam");
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

    // The 4.1 writer records membership through $Entities, one entity per
    // CELL BLOCK (v11.5.0, roadmap §1 tier B3): a region covering only PART
    // of a block (as here -- "solid" tags cell 0 of the mesh's two
    // tetrahedra) cannot be represented that way and is dropped with a
    // warning, keeping only the $PhysicalNames row (see
    // Gmsh.RegionTagAllocation* in test_gmsh.cpp for the block-aligned case,
    // where 4.1 now does keep membership).
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
