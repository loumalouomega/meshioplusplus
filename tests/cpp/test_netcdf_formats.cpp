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

TEST(Exodus, IsAnOptionsAwareReader) {
    // Registering exodus in registry_readers_ex()/registry_metadata_readers() is
    // what lets any flat binding pass a time step at all; before it, this was
    // false and `--time-step` had nowhere to go.
    EXPECT_TRUE(meshioplusplus::registry_reader_supports_options("exodus"));
}

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
