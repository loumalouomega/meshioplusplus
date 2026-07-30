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

// Project includes
#include "mesh_fixtures.hpp"

#ifdef MESHIOPLUSPLUS_HAS_NETCDF

#include <netcdf.h>

#include <cstddef>
#include <filesystem>
#include <set>
#include <string>

#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/formats/exodus.hpp"
#include "meshioplusplus/registry.hpp"

TEST(Exodus, Basic) {
    auto w = [](const std::string& p, const mt::Mesh& m) { meshioplusplus::write_exodus(p, m); };
    auto r = [](const std::string& p) { return meshioplusplus::read_exodus(p); };
    mt::roundtrip(w, r, mt::tri_mesh(), ".e");
    mt::roundtrip(w, r, mt::tet_mesh(), ".e");
    mt::roundtrip(w, r, mt::hex_mesh(), ".e");
    mt::roundtrip(w, r, mt::tri_quad_mesh(), ".e");
}

// The Exodus side->facet tables are a transcription from the Exodus II spec's
// side-node lists, and a wrong entry still yields a *valid-looking* facet index
// -- a Side region pointing confidently at the wrong face. So rather than trust
// the table, check the property it is supposed to have: the facet it names must
// have exactly the corner nodes the Exodus side did.
namespace {

/// The Exodus II spec's side->node lists (1-based node numbers within the cell).
struct ExodusSideNodes {
    const char* mCellType;
    int mNumSides;
    std::vector<std::vector<int>> mSides;
};

const std::vector<ExodusSideNodes>& exodus_side_node_lists() {
    static const std::vector<ExodusSideNodes> lists = {
        {"tetra", 4, {{1, 2, 4}, {2, 3, 4}, {1, 4, 3}, {1, 3, 2}}},
        {"hexahedron",
         6,
         {{1, 2, 6, 5}, {2, 3, 7, 6}, {3, 4, 8, 7}, {4, 1, 5, 8}, {1, 4, 3, 2}, {5, 6, 7, 8}}},
        {"wedge", 5, {{1, 2, 5, 4}, {2, 3, 6, 5}, {1, 4, 6, 3}, {1, 3, 2}, {4, 5, 6}}},
        {"pyramid", 5, {{1, 2, 5}, {2, 3, 5}, {3, 4, 5}, {4, 1, 5}, {1, 4, 3, 2}}},
    };
    return lists;
}

}  // namespace

TEST(Exodus, FaceIndexTablesMatchCellFaces) {
    for (const ExodusSideNodes& entry : exodus_side_node_lists()) {
        const meshioplusplus::CellType type = meshioplusplus::cell_type_from_name(entry.mCellType);
        const std::vector<meshioplusplus::detail::CellFaceDef>& faces =
            meshioplusplus::detail::cell_faces(type);
        ASSERT_EQ(static_cast<int>(faces.size()), entry.mNumSides) << entry.mCellType;

        std::set<int> mapped;
        for (int side = 1; side <= entry.mNumSides; ++side) {
            const int facet = meshioplusplus::exo_face_index(entry.mCellType, side);
            ASSERT_GE(facet, 0) << entry.mCellType << " side " << side << " is unmapped";
            ASSERT_LT(facet, static_cast<int>(faces.size()));

            // The corner-node SET must match; winding may legitimately differ,
            // since meshio++ orients every facet outward and Exodus does not.
            std::set<int> expected;
            for (int n : entry.mSides[static_cast<std::size_t>(side - 1)])
                expected.insert(n - 1);  // Exodus node numbers are 1-based
            std::set<int> got;
            const meshioplusplus::detail::CellFaceDef& face =
                faces[static_cast<std::size_t>(facet)];
            for (std::uint8_t i = 0; i < face.mNumCorners; ++i)
                got.insert(static_cast<int>(face.mNodes[i]));

            EXPECT_EQ(got, expected) << entry.mCellType << " Exodus side " << side
                                     << " maps to facet " << facet << ", which is a different face";
            mapped.insert(facet);
        }
        // Every facet must be hit exactly once: a table with a duplicate would
        // leave one real face unreachable from any side set.
        EXPECT_EQ(mapped.size(), faces.size())
            << entry.mCellType << ": the side->facet map is not a bijection";
    }
}

TEST(Exodus, FaceIndexRejectsUnknownPairs) {
    EXPECT_EQ(meshioplusplus::exo_face_index("hexahedron", 0), -1);
    EXPECT_EQ(meshioplusplus::exo_face_index("hexahedron", 7), -1);
    EXPECT_EQ(meshioplusplus::exo_face_index("line", 1), -1);
    // Higher-order variants share their linear base's facet ordering.
    EXPECT_EQ(meshioplusplus::exo_face_index("hexahedron20", 4),
              meshioplusplus::exo_face_index("hexahedron", 4));
    EXPECT_EQ(meshioplusplus::exo_face_index("tetra10", 3),
              meshioplusplus::exo_face_index("tetra", 3));
}

namespace {

/// Write a minimal Exodus file of one-node SPHERE elements, straight through the
/// netCDF C API so `elem_type`'s length can be chosen exactly.
///
/// `ElemTypeLen` of 7 for "SPHERE" is what NetCDF.jl -- and so PeriLab, and any
/// other Julia-written peridynamics output -- actually stores: the C string's
/// terminating NUL counted as part of the attribute. netCDF4-python strips it,
/// which is why only the C++ reader ever failed on these files and why this
/// fixture has to be built here rather than in the Python suite.
void write_sphere_file(const std::string& rPath, std::size_t ElemTypeLen, bool WithAttribute) {
    int ncid;
    ASSERT_EQ(nc_create(rPath.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid), NC_NOERR);
    int d_nodes, d_dim, d_el, d_npe, d_att, d_str;
    ASSERT_EQ(nc_def_dim(ncid, "num_nodes", 3, &d_nodes), NC_NOERR);
    ASSERT_EQ(nc_def_dim(ncid, "num_dim", 3, &d_dim), NC_NOERR);
    ASSERT_EQ(nc_def_dim(ncid, "num_el_in_blk1", 3, &d_el), NC_NOERR);
    ASSERT_EQ(nc_def_dim(ncid, "num_nod_per_el1", 1, &d_npe), NC_NOERR);
    ASSERT_EQ(nc_def_dim(ncid, "len_string", 33, &d_str), NC_NOERR);

    int coord_dims[2] = {d_dim, d_nodes};
    int v_coord;
    ASSERT_EQ(nc_def_var(ncid, "coord", NC_DOUBLE, 2, coord_dims, &v_coord), NC_NOERR);
    const double coords[9] = {0, 1, 2, 0, 0, 0, 0, 0, 0};
    ASSERT_EQ(nc_put_var_double(ncid, v_coord, coords), NC_NOERR);

    int conn_dims[2] = {d_el, d_npe};
    int v_conn;
    ASSERT_EQ(nc_def_var(ncid, "connect1", NC_INT, 2, conn_dims, &v_conn), NC_NOERR);
    // The NUL is *inside* the value, so the length is what decides what arrives.
    ASSERT_EQ(nc_put_att_text(ncid, v_conn, "elem_type", ElemTypeLen, "SPHERE\0"), NC_NOERR);
    const int conn[3] = {1, 2, 3};
    ASSERT_EQ(nc_put_var_int(ncid, v_conn, conn), NC_NOERR);

    if (WithAttribute) {
        ASSERT_EQ(nc_def_dim(ncid, "num_att_in_blk1", 1, &d_att), NC_NOERR);
        int att_dims[2] = {d_el, d_att};
        int v_att;
        ASSERT_EQ(nc_def_var(ncid, "attrib1", NC_DOUBLE, 2, att_dims, &v_att), NC_NOERR);
        const double radii[3] = {0.5, 0.25, 0.125};
        ASSERT_EQ(nc_put_var_double(ncid, v_att, radii), NC_NOERR);
        int name_dims[2] = {d_att, d_str};
        int v_name;
        ASSERT_EQ(nc_def_var(ncid, "attrib_name1", NC_CHAR, 2, name_dims, &v_name), NC_NOERR);
        std::size_t start[2] = {0, 0}, count[2] = {1, 6};
        ASSERT_EQ(nc_put_vara_text(ncid, v_name, start, count, "RADIUS"), NC_NOERR);
    }
    ASSERT_EQ(nc_close(ncid), NC_NOERR);
}

}  // namespace

// The regression behind
// https://github.com/loumalouomega/VSCode-MDPA-Preview/issues/63: a SPHERE block
// whose `elem_type` counts its terminating NUL used to fail the read with
// "unknown element type SPHERE" -- the NUL invisible in the message, since
// `what()` is a `const char*` that stops at it. Only the C++ path was affected,
// so the Python shim's fallback hid it everywhere except WASM.
TEST(Exodus, NulTerminatedElemTypeIsRead) {
    for (std::size_t len : {std::size_t(6), std::size_t(7)}) {
        std::string path = mt::temp_path(".e");
        write_sphere_file(path, len, false);
        meshioplusplus::Mesh mesh = meshioplusplus::read_exodus(path);
        ASSERT_EQ(mesh.NumCellBlocks(), 1u) << "elem_type length " << len;
        EXPECT_EQ(mesh.Cells(0).Type(), "vertex") << "elem_type length " << len;
        EXPECT_EQ(mesh.Cells(0).NumCells(), 3u);
        std::filesystem::remove(path);
    }
}

TEST(Exodus, ElementAttributesRoundTripAsCellData) {
    const std::string name = std::string(meshioplusplus::kExodusAttributePrefix) + "RADIUS";
    std::string path = mt::temp_path(".e");
    write_sphere_file(path, 7, true);
    meshioplusplus::Mesh mesh = meshioplusplus::read_exodus(path);
    std::filesystem::remove(path);

    ASSERT_TRUE(mesh.HasCellData(name));
    ASSERT_EQ(mesh.CellDataNumBlocks(name), 1u);
    const meshioplusplus::NDArray& radii = mesh.CellData(name, 0);
    ASSERT_EQ(radii.Size(), 3u);
    EXPECT_DOUBLE_EQ(radii.As<double>()[0], 0.5);
    EXPECT_DOUBLE_EQ(radii.As<double>()[2], 0.125);

    // ... and back out again: the attribute is the one thing a sphere viewer
    // needs, so losing it on write would make the read half pointless.
    std::string out = mt::temp_path(".e");
    meshioplusplus::write_exodus(out, mesh);
    meshioplusplus::Mesh back = meshioplusplus::read_exodus(out);
    std::filesystem::remove(out);
    ASSERT_TRUE(back.HasCellData(name));
    const meshioplusplus::NDArray& again = back.CellData(name, 0);
    ASSERT_EQ(again.Size(), 3u);
    EXPECT_DOUBLE_EQ(again.As<double>()[0], 0.5);
    EXPECT_DOUBLE_EQ(again.As<double>()[2], 0.125);
}

TEST(Exodus, CellDataRoundTripsAsElementVariables) {
    // v9.9.0: this writer emitted no `vals_elem_var` at all, so every ordinary
    // cell_data array was silently dropped while point_data round-tripped.
    std::string p = mt::temp_path(".e");
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}, {0, 2, 3}}));

    meshioplusplus::NDArray mat(meshioplusplus::DType::Float64, {2});
    mat.As<double>()[0] = 7.0;
    mat.As<double>()[1] = 8.0;
    m.AddCellData("material", {std::move(mat)});
    meshioplusplus::NDArray pd(meshioplusplus::DType::Float64, {4});
    for (std::size_t i = 0; i < 4; ++i)
        pd.As<double>()[i] = static_cast<double>(i);
    m.AddPointData("T", std::move(pd));

    meshioplusplus::write_exodus(p, m);
    meshioplusplus::Mesh back = meshioplusplus::read_exodus(p);
    std::filesystem::remove(p);

    ASSERT_TRUE(back.HasCellData("material"));
    ASSERT_EQ(back.CellDataNumBlocks("material"), 1u);
    const meshioplusplus::NDArray& got = back.CellData("material", 0);
    ASSERT_EQ(got.Size(), 2u);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(got, 0), 7.0);
    EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(got, 1), 8.0);
    ASSERT_TRUE(back.HasPointData("T"));  // unchanged by the new element vars
}

TEST(Exodus, MultiComponentCellDataRoundTrips) {
    // Also the regression test for a heap buffer overflow: the reader's
    // cell_data assembly allocated `{total}` (scalar) but memcpy'd each block's
    // full `Nbytes()`, so a k>1 element variable wrote k*n bytes into an
    // n-byte buffer. No real SEACAS file has one (standard Exodus element
    // variables are scalar per element), so nothing reached it until this
    // writer started emitting the trailing dims the nodal path already did.
    std::string p = mt::temp_path(".e");
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}, {0, 2, 3}}));
    meshioplusplus::NDArray s(meshioplusplus::DType::Float64, {2, 3});
    for (std::size_t i = 0; i < 6; ++i)
        s.As<double>()[i] = static_cast<double>(i);
    m.AddCellData("stress", {std::move(s)});

    meshioplusplus::write_exodus(p, m);
    meshioplusplus::Mesh back = meshioplusplus::read_exodus(p);
    std::filesystem::remove(p);

    ASSERT_TRUE(back.HasCellData("stress"));
    const meshioplusplus::NDArray& got = back.CellData("stress", 0);
    ASSERT_EQ(got.Shape().size(), 2u);
    EXPECT_EQ(got.Shape()[0], 2u);
    EXPECT_EQ(got.Shape()[1], 3u);
    for (std::size_t i = 0; i < 6; ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(got, i), static_cast<double>(i));
}

TEST(Exodus, BlockNamesRoundTripThroughCellRegions) {
    // The writer emitted no `eb_names`, so a block name became the reader's
    // synthetic "Block N" on the way back. It is now recovered from the Cell
    // region the reader itself derives from `eb_names`.
    std::string p = mt::temp_path(".e");
    meshioplusplus::Mesh m;
    m.AssignPoints(
        mt::points_from({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {2, 0, 0}, {2, 1, 0}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    m.AddCellBlock("quad", mt::conn_from({{1, 4, 5, 2}}));
    // Global block-major cell indices: block 0 is cell 0, block 1 is cell 1.
    meshioplusplus::NDArray a(meshioplusplus::DType::Int64, {1});
    a.As<std::int64_t>()[0] = 0;
    m.AddRegion(meshioplusplus::Region("skin", meshioplusplus::RegionKind::Cell, std::move(a)));
    meshioplusplus::NDArray b(meshioplusplus::DType::Int64, {1});
    b.As<std::int64_t>()[0] = 1;
    m.AddRegion(meshioplusplus::Region("body", meshioplusplus::RegionKind::Cell, std::move(b)));

    meshioplusplus::write_exodus(p, m);
    meshioplusplus::Mesh back = meshioplusplus::read_exodus(p);
    std::filesystem::remove(p);

    EXPECT_TRUE(back.HasRegion("skin", meshioplusplus::RegionKind::Cell));
    EXPECT_TRUE(back.HasRegion("body", meshioplusplus::RegionKind::Cell));
    EXPECT_FALSE(back.HasRegion("Block 1", meshioplusplus::RegionKind::Cell));
}

TEST(Exodus, TimeValueComesFromFieldData) {
    // The single step's recorded time was hard-coded 0; `exodus:time` supplies
    // it so one frame of a transient solve can be labelled correctly. A real
    // multi-step writer remains a follow-up (XdmfTimeSeriesWriter's shape).
    std::string p = mt::temp_path(".e");
    meshioplusplus::Mesh m = mt::tri_mesh();
    meshioplusplus::NDArray t(meshioplusplus::DType::Float64, {1});
    t.As<double>()[0] = 2.5;
    m.AddFieldData("exodus:time", std::move(t));
    meshioplusplus::write_exodus(p, m);

    meshioplusplus::MeshMetadata md = meshioplusplus::read_exodus_metadata(p);
    std::filesystem::remove(p);
    ASSERT_EQ(md.mTimeValues.size(), 1u);
    EXPECT_NEAR(md.mTimeValues[0], 2.5, 1e-6);  // stored as NC_FLOAT
}

TEST(Exodus, IsAnOptionsAwareReader) {
    // Registering exodus in registry_readers_ex()/registry_metadata_readers() is
    // what lets any flat binding pass a time step at all; before it, this was
    // false and `--time-step` had nowhere to go.
    EXPECT_TRUE(meshioplusplus::registry_reader_supports_options("exodus"));
}

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
