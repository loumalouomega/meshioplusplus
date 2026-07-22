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
#include <filesystem>
#include <fstream>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/xdmf.hpp"

namespace {
void rt(const mt::Mesh& mesh, const std::string& data_format) {
    mt::roundtrip([&](const std::string& p,
                      const mt::Mesh& m) { meshioplusplus::write_xdmf(p, m, data_format, -1); },
                  [](const std::string& p) { return meshioplusplus::read_xdmf(p); }, mesh, ".xdmf");
}
}  // namespace

TEST(Xdmf, Xml) {
    rt(mt::tri_mesh(), "XML");
    rt(mt::tet_mesh(), "XML");
    rt(mt::tri_quad_mesh(), "XML");  // Mixed topology
    rt(mt::tri_mesh_2d(), "XML");    // XY geometry
}
TEST(Xdmf, Binary) {
    rt(mt::tri_mesh(), "Binary");
    rt(mt::hex_mesh(), "Binary");
    rt(mt::tri_quad_mesh(), "Binary");
}
#ifdef MESHIOPLUSPLUS_HAS_HDF5
TEST(Xdmf, Hdf) {
    rt(mt::tri_mesh(), "HDF");
    rt(mt::tet_mesh(), "HDF");
    rt(mt::tri_quad_mesh(), "HDF");
}
#endif

TEST(Xdmf, ReadRejectsMissingRoot) {
    // Well-formed XML that lacks the <Xdmf> root must raise rather than be
    // treated as an empty mesh.
    std::string path = mt::temp_path(".xdmf");
    {
        std::ofstream f(path);
        f << "<?xml version=\"1.0\"?>\n<NotXdmf/>\n";
    }
    EXPECT_THROW(meshioplusplus::read_xdmf(path), meshioplusplus::ReadError);
    std::error_code ec;
    std::filesystem::remove(path, ec);
}
