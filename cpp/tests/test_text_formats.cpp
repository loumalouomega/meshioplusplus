// Round-trip tests for the plain-text single-signature formats:
// medit, nastran, abaqus, avsucd, permas, tecplot, ugrid, netgen.

#include <gtest/gtest.h>

#include "mesh_fixtures.hpp"
#include "meshio/formats/abaqus.hpp"
#include "meshio/formats/avsucd.hpp"
#include "meshio/formats/medit.hpp"
#include "meshio/formats/nastran.hpp"
#include "meshio/formats/netgen.hpp"
#include "meshio/formats/permas.hpp"
#include "meshio/formats/tecplot.hpp"
#include "meshio/formats/ugrid.hpp"

// Generic helper for `void write_X(path, mesh)` / `Mesh read_X(path)`.
#define SIMPLE_RT(WRITER, READER, MESH, SUFFIX, ATOL)                             \
    mt::roundtrip([](const std::string& p, const mt::Mesh& m) { WRITER(p, m); },  \
                  [](const std::string& p) { return READER(p); }, MESH, SUFFIX, ATOL)

TEST(Medit, Basic) {
    SIMPLE_RT(meshio::write_medit_ascii, meshio::read_medit_ascii, mt::tri_mesh(), ".mesh", 1e-12);
    SIMPLE_RT(meshio::write_medit_ascii, meshio::read_medit_ascii, mt::tet_mesh(), ".mesh", 1e-12);
    SIMPLE_RT(meshio::write_medit_ascii, meshio::read_medit_ascii, mt::hex_mesh(), ".mesh", 1e-12);
    SIMPLE_RT(meshio::write_medit_ascii, meshio::read_medit_ascii, mt::quad_mesh(), ".mesh", 1e-12);
}

TEST(Nastran, Basic) {
    // The C++ reader is sentinel-gated to meshio-written files, so a self
    // round-trip is the supported path.
    SIMPLE_RT(meshio::write_nastran, meshio::read_nastran, mt::tri_mesh(), ".bdf", 1e-10);
    SIMPLE_RT(meshio::write_nastran, meshio::read_nastran, mt::tet_mesh(), ".bdf", 1e-10);
    SIMPLE_RT(meshio::write_nastran, meshio::read_nastran, mt::hex_mesh(), ".bdf", 1e-10);
}

TEST(Abaqus, Basic) {
    SIMPLE_RT(meshio::write_abaqus, meshio::read_abaqus, mt::tri_mesh(), ".inp", 1e-12);
    SIMPLE_RT(meshio::write_abaqus, meshio::read_abaqus, mt::tet_mesh(), ".inp", 1e-12);
    SIMPLE_RT(meshio::write_abaqus, meshio::read_abaqus, mt::hex_mesh(), ".inp", 1e-12);
}

TEST(Avsucd, Basic) {
    SIMPLE_RT(meshio::write_avsucd, meshio::read_avsucd, mt::tri_mesh(), ".avs", 1e-12);
    SIMPLE_RT(meshio::write_avsucd, meshio::read_avsucd, mt::tet_mesh(), ".avs", 1e-12);
    SIMPLE_RT(meshio::write_avsucd, meshio::read_avsucd, mt::hex_mesh(), ".avs", 1e-12);
}

TEST(Permas, Basic) {
    SIMPLE_RT(meshio::write_permas, meshio::read_permas, mt::tri_mesh(), ".post", 1e-12);
    SIMPLE_RT(meshio::write_permas, meshio::read_permas, mt::tri_quad_mesh(), ".post", 1e-12);
    SIMPLE_RT(meshio::write_permas, meshio::read_permas, mt::hex_mesh(), ".post", 1e-12);
}

TEST(Tecplot, SingleType) {
    SIMPLE_RT(meshio::write_tecplot, meshio::read_tecplot, mt::tri_mesh(), ".dat", 1e-12);
    SIMPLE_RT(meshio::write_tecplot, meshio::read_tecplot, mt::quad_mesh(), ".dat", 1e-12);
    SIMPLE_RT(meshio::write_tecplot, meshio::read_tecplot, mt::tet_mesh(), ".dat", 1e-12);
    SIMPLE_RT(meshio::write_tecplot, meshio::read_tecplot, mt::hex_mesh(), ".dat", 1e-12);
}

TEST(Ugrid, Basic) {
    SIMPLE_RT(meshio::write_ugrid, meshio::read_ugrid, mt::tri_mesh(), ".ugrid", 1e-12);
    SIMPLE_RT(meshio::write_ugrid, meshio::read_ugrid, mt::tet_mesh(), ".ugrid", 1e-12);
    SIMPLE_RT(meshio::write_ugrid, meshio::read_ugrid, mt::hex_mesh(), ".ugrid", 1e-12);
}

TEST(Netgen, Basic) {
    auto w = [](const std::string& p, const mt::Mesh& m) {
        meshio::write_netgen(p, m, ".16e");
    };
    auto r = [](const std::string& p) { return meshio::read_netgen(p); };
    mt::roundtrip(w, r, mt::tri_mesh(), ".vol");
    mt::roundtrip(w, r, mt::tet_mesh(), ".vol");
    mt::roundtrip(w, r, mt::hex_mesh(), ".vol");
}
