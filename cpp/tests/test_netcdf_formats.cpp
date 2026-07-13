// Round-trip tests for the netCDF-backed Exodus format. Compiled to nothing
// unless the extension is built with MESHIO_HAS_NETCDF.

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"

#ifdef MESHIO_HAS_NETCDF

#include "meshio/formats/exodus.hpp"

TEST(Exodus, Basic) {
    auto w = [](const std::string& p, const mt::Mesh& m) { meshio::write_exodus(p, m); };
    auto r = [](const std::string& p) { return meshio::read_exodus(p); };
    mt::roundtrip(w, r, mt::tri_mesh(), ".e");
    mt::roundtrip(w, r, mt::tet_mesh(), ".e");
    mt::roundtrip(w, r, mt::hex_mesh(), ".e");
    mt::roundtrip(w, r, mt::tri_quad_mesh(), ".e");
}

#endif  // MESHIO_HAS_NETCDF
