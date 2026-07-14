// Round-trip tests for the netCDF-backed Exodus format. Compiled to nothing
// unless the extension is built with MESHIOPLUSPLUS_HAS_NETCDF.

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"

#ifdef MESHIOPLUSPLUS_HAS_NETCDF

#include "meshioplusplus/formats/exodus.hpp"

TEST(Exodus, Basic) {
    auto w = [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_exodus(p, m); };
    auto r = [](const std::string& p) { return meshioplusplus::read_exodus(p); };
    mt::roundtrip(w, r, mt::tri_mesh(), ".e");
    mt::roundtrip(w, r, mt::tet_mesh(), ".e");
    mt::roundtrip(w, r, mt::hex_mesh(), ".e");
    mt::roundtrip(w, r, mt::tri_quad_mesh(), ".e");
}

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
