// Round-trip tests for the HDF5-backed formats. Compiled to nothing unless the
// extension is built with MESHIO_HAS_HDF5 (else the Python fallback handles
// these formats and there is no C++ path to test).

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"

#ifdef MESHIO_HAS_HDF5

#include "meshio/formats/cgns.hpp"
#include "meshio/formats/h5m.hpp"
#include "meshio/formats/hmf.hpp"
#include "meshio/formats/med.hpp"

TEST(Cgns, TetraCompressed) {
    for (int gzip : {-1, 4}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshio::write_cgns(p, m, gzip);
        };
        auto r = [](const std::string& p) { return meshio::read_cgns(p); };
        mt::roundtrip(w, r, mt::tet_mesh(), ".cgns");
    }
}

TEST(H5m, LineTriangleTetra) {
    for (int gzip : {-1, 4}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshio::write_h5m(p, m, true, gzip);
        };
        auto r = [](const std::string& p) { return meshio::read_h5m(p); };
        mt::roundtrip(w, r, mt::line_mesh(), ".h5m");
        mt::roundtrip(w, r, mt::tri_mesh(), ".h5m");
        mt::roundtrip(w, r, mt::tet_mesh(), ".h5m");
    }
}

TEST(Hmf, Basic) {
    for (int gzip : {-1, 4}) {
        auto w = [=](const std::string& p, const mt::Mesh& m) {
            meshio::write_hmf(p, m, gzip);
        };
        auto r = [](const std::string& p) { return meshio::read_hmf(p); };
        mt::roundtrip(w, r, mt::tri_mesh(), ".hmf");
        mt::roundtrip(w, r, mt::tet_mesh(), ".hmf");
        mt::roundtrip(w, r, mt::hex_mesh(), ".hmf");
    }
}

TEST(Med, Basic) {
    auto w = [](const std::string& p, const mt::Mesh& m) {
        meshio::write_med(p, m, meshio::MedInfo{});
    };
    auto r = [](const std::string& p) {
        meshio::MedInfo info;
        return meshio::read_med(p, info);
    };
    mt::roundtrip(w, r, mt::tri_mesh(), ".med");
    mt::roundtrip(w, r, mt::tet_mesh(), ".med");
    mt::roundtrip(w, r, mt::hex_mesh(), ".med");
}

#endif  // MESHIO_HAS_HDF5
