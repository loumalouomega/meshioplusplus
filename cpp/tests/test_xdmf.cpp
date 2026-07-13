// XDMF round-trip tests: XML and Binary data paths always; HDF when built with
// MESHIO_HAS_HDF5.

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshio/formats/xdmf.hpp"

namespace {
void rt(const mt::Mesh& mesh, const std::string& data_format) {
    mt::roundtrip(
        [&](const std::string& p, const mt::Mesh& m) {
            meshio::write_xdmf(p, m, data_format, -1);
        },
        [](const std::string& p) { return meshio::read_xdmf(p); }, mesh, ".xdmf");
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
#ifdef MESHIO_HAS_HDF5
TEST(Xdmf, Hdf) {
    rt(mt::tri_mesh(), "HDF");
    rt(mt::tet_mesh(), "HDF");
    rt(mt::tri_quad_mesh(), "HDF");
}
#endif
