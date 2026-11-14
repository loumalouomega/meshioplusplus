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
// The ZIP container (roadmap 2.4.5.1) and 3MF packages (2.4.5.2).

#include <cstdio>
#include <string>

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/threemf.hpp"
#include "../../src/cpp/src/detail/zip_archive.hpp"

TEST(ZipArchive, MembersRoundTripStoredAndDeflated) {
    meshioplusplus::detail::ZipWriter w;
    const std::string big(5000, 'a');
    w.Add("a.txt", "hello");
    w.Add("dir/b.txt", big);
    w.Add("c.bin", std::string("\0\1\2", 3), /*Deflate=*/false);
    const std::string bytes = w.Finish();
    const meshioplusplus::detail::ZipReader r(bytes);
    ASSERT_EQ(r.Members().size(), 3u);
    EXPECT_EQ(r.Read("a.txt"), "hello");
    EXPECT_EQ(r.Read("dir/b.txt"), big);
    EXPECT_EQ(r.Read("c.bin"), std::string("\0\1\2", 3));
    EXPECT_EQ(r.Find("nope"), nullptr);
    EXPECT_THROW(r.Read("nope"), meshioplusplus::ReadError);
    EXPECT_THROW(meshioplusplus::detail::ZipReader("not a zip archive at all, no"),
                 meshioplusplus::ReadError);
}

TEST(ThreeMf, ATriangleMeshRoundTrips) {
    const auto mesh = mt::tri_mesh();
    const std::string p = mt::temp_path("_part.3mf");
    meshioplusplus::write_3mf(p, mesh);
    const auto back = meshioplusplus::read_3mf(p);
    ASSERT_EQ(back.NumPoints(), mesh.NumPoints());
    EXPECT_EQ(back.Cells(0).NumCells(), mesh.Cells(0).NumCells());
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.FieldData("meshio:length_unit"), 0),
                     1e-3);
    EXPECT_THROW(meshioplusplus::write_3mf(p, mt::tet_mesh()), meshioplusplus::WriteError);
    std::remove(p.c_str());
}
