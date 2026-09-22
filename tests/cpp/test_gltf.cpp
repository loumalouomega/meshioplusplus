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

// System includes
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "meshioplusplus/detail/colormap.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/gltf.hpp"
#include "meshioplusplus/region.hpp"
#include "mesh_fixtures.hpp"

using namespace meshioplusplus;

namespace {

using mt::make_mesh;

std::string slurp(const std::string& rPath) {
    std::ifstream f(rPath, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::uint32_t u32(const std::string& rData, std::size_t Offset) {
    std::uint32_t v = 0;
    for (int k = 3; k >= 0; --k)
        v = (v << 8) | static_cast<unsigned char>(rData[Offset + static_cast<std::size_t>(k)]);
    return v;
}

// A decoded GLB: the JSON text (padding removed) and the BIN chunk.
struct Glb {
    std::string mJson;
    std::string mBin;
    std::size_t mFileSize = 0;
};

Glb read_glb(const std::string& rPath) {
    const std::string data = slurp(rPath);
    Glb glb;
    glb.mFileSize = data.size();
    EXPECT_GE(data.size(), 20u);
    EXPECT_EQ(u32(data, 0), 0x46546C67u);
    EXPECT_EQ(u32(data, 4), 2u);
    EXPECT_EQ(u32(data, 8), data.size()) << "the header length is the file length";
    std::size_t off = 12;
    int chunk = 0;
    while (off < data.size()) {
        const std::uint32_t len = u32(data, off);
        const std::uint32_t type = u32(data, off + 4);
        EXPECT_EQ(len % 4, 0u) << "chunks are 4-byte aligned";
        const std::string body = data.substr(off + 8, len);
        if (chunk == 0) {
            EXPECT_EQ(type, 0x4E4F534Au);
            std::size_t end = body.size();
            while (end > 0 && body[end - 1] == ' ')
                --end;
            glb.mJson = body.substr(0, end);
        } else {
            EXPECT_EQ(type, 0x004E4942u);
            glb.mBin = body;
        }
        off += 8 + len;
        ++chunk;
    }
    EXPECT_EQ(off, data.size());
    return glb;
}

std::size_t count_of(const std::string& rText, const std::string& rNeedle) {
    std::size_t n = 0;
    for (std::size_t pos = rText.find(rNeedle); pos != std::string::npos;
         pos = rText.find(rNeedle, pos + 1))
        ++n;
    return n;
}

bool has(const std::string& rText, const std::string& rNeedle) {
    return rText.find(rNeedle) != std::string::npos;
}

Mesh cube_quads() {
    return make_mesh(
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
        "quad",
        {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 4, 7, 3}, {1, 2, 6, 5}});
}

NDArray i64(const std::vector<std::int64_t>& rVals) {
    NDArray a = NDArray::Uninit(DType::Int64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<std::int64_t>()[i] = rVals[i];
    return a;
}

NDArray f64(const std::vector<double>& rVals) {
    NDArray a = NDArray::Uninit(DType::Float64, {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<double>()[i] = rVals[i];
    return a;
}

}  // namespace

TEST(Gltf, GlbContainerRules) {
    const std::string path = mt::temp_path(".glb");
    write_gltf(path, cube_quads());
    const Glb glb = read_glb(path);
    EXPECT_EQ(glb.mJson.rfind("{\"asset\":{\"generator\":\"Written by meshio++ v", 0), 0u);
    EXPECT_TRUE(has(glb.mJson, "\"version\":\"2.0\""));
    EXPECT_FALSE(has(glb.mJson, "\"uri\""));
    EXPECT_FALSE(glb.mBin.empty());
    std::remove(path.c_str());
}

TEST(Gltf, JsonContainerWritesASidecar) {
    const std::string path = mt::temp_path(".gltf");
    write_gltf(path, cube_quads());
    const std::string json = slurp(path);
    const std::string bin_path = path.substr(0, path.size() - 5) + ".bin";
    const std::string bin = slurp(bin_path);
    EXPECT_EQ(json.front(), '{');
    EXPECT_TRUE(has(json, "\"byteLength\":" + std::to_string(bin.size())));
    const std::size_t cut = bin_path.find_last_of('/');
    EXPECT_TRUE(has(json, "\"uri\":\"" + bin_path.substr(cut + 1) + "\""));

    // The suffix picks the container unless it is overridden.
    GltfWriteOptions forced;
    forced.mContainer = GltfContainer::Binary;
    write_gltf(path, cube_quads(), forced);
    EXPECT_EQ(slurp(path).substr(0, 4), "glTF");
    std::remove(path.c_str());
    std::remove(bin_path.c_str());
}

TEST(Gltf, ACubeOfQuadsHasACreasePerEdge) {
    const std::string path = mt::temp_path(".glb");
    write_gltf(path, cube_quads());
    const Glb glb = read_glb(path);
    EXPECT_TRUE(has(glb.mJson, "\"count\":36,\"type\":\"SCALAR\""));  // 12 triangles
    EXPECT_TRUE(has(glb.mJson, "\"count\":24,\"type\":\"VEC3\""));    // 6 faces x 4 corners
    EXPECT_TRUE(has(glb.mJson, "\"mode\":4"));
    EXPECT_TRUE(has(glb.mJson, "\"NORMAL\""));
    EXPECT_TRUE(has(glb.mJson, "\"max\":[0.5,0.5,0.5],\"min\":[-0.5,-0.5,-0.5]"))
        << "recentred on the bounding-box centre";

    // A smooth split angle shares the corners.
    GltfWriteOptions smooth;
    smooth.mSplitAngle = 180.0;
    write_gltf(path, cube_quads(), smooth);
    EXPECT_TRUE(has(read_glb(path).mJson, "\"count\":8,\"type\":\"VEC3\""));
    std::remove(path.c_str());
}

TEST(Gltf, TheRootTransformCarriesTheAxisAndTheOffset) {
    const std::string path = mt::temp_path(".glb");
    write_gltf(path, cube_quads());
    // z-up source: -90 degrees about x; the centre (0.5, 0.5, 0.5) -> (0.5, 0.5, -0.5)
    EXPECT_TRUE(has(read_glb(path).mJson,
                    "\"rotation\":[-0.70710678118654757,0,0,0.70710678118654757],"
                    "\"translation\":[0.5,0.5,-0.5]"));

    GltfWriteOptions y;
    y.mUpAxis = GltfUpAxis::Y;
    y.mScale = 0.001;
    write_gltf(path, cube_quads(), y);
    const std::string json = read_glb(path).mJson;
    EXPECT_FALSE(has(json, "\"rotation\""));
    EXPECT_TRUE(has(json, "\"scale\":[0.001,0.001,0.001]"));

    // A flat mesh is y-up by itself.
    Mesh flat = mt::tri_mesh_2d();
    write_gltf(path, flat);
    EXPECT_FALSE(has(read_glb(path).mJson, "\"rotation\""));
    std::remove(path.c_str());
}

TEST(Gltf, VolumesExportTheirSkin) {
    const std::string path = mt::temp_path(".glb");
    write_gltf(path,
               make_mesh({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, "tetra", {{0, 1, 2, 3}}));
    EXPECT_TRUE(has(read_glb(path).mJson, "\"count\":12,\"type\":\"VEC3\""));
    // two tets sharing a face: the shared face is interior, 6 boundary triangles
    write_gltf(path, mt::tet_mesh());
    EXPECT_TRUE(has(read_glb(path).mJson, "\"count\":18,\"type\":\"SCALAR\""));
    write_gltf(path, mt::hex_mesh());
    EXPECT_TRUE(has(read_glb(path).mJson, "\"count\":24,\"type\":\"VEC3\""));
    std::remove(path.c_str());
}

TEST(Gltf, ABoundaryPatchWinsOverTheSkinAndNodesFollowRegions) {
    Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}));
    m.AddCellBlock("tetra", mt::conn_from({{0, 1, 2, 3}}));
    m.AddCellBlock("triangle", mt::conn_from({{0, 1, 2}}));
    m.AddRegion(Region("fluid", RegionKind::Cell, i64({0})));
    m.AddRegion(Region("wall", RegionKind::Cell, i64({1})));
    const std::string path = mt::temp_path(".glb");
    write_gltf(path, m);
    const std::string json = read_glb(path).mJson;
    EXPECT_EQ(count_of(json, "\"name\":\"fluid\""), 2u);  // node + mesh
    EXPECT_EQ(count_of(json, "\"name\":\"wall\""), 2u);
    // the skin has 4 facets; the wall cell replaces one, so 3 + 1 triangles
    EXPECT_TRUE(has(json, "\"count\":9,\"type\":\"SCALAR\""));
    EXPECT_TRUE(has(json, "\"count\":3,\"type\":\"SCALAR\""));

    GltfWriteOptions one;
    one.mByRegion = false;
    write_gltf(path, m, one);
    EXPECT_TRUE(has(read_glb(path).mJson, "\"name\":\"mesh\""));
    std::remove(path.c_str());
}

TEST(Gltf, PointCloudsAndLinesUseTheirOwnPrimitives) {
    const std::string path = mt::temp_path(".glb");
    Mesh cloud;
    cloud.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 2, 0}}));
    write_gltf(path, cloud);
    std::string json = read_glb(path).mJson;
    EXPECT_TRUE(has(json, "\"mode\":0"));
    EXPECT_FALSE(has(json, "\"indices\""));

    write_gltf(path, mt::line_mesh());
    json = read_glb(path).mJson;
    EXPECT_TRUE(has(json, "\"mode\":1"));
    EXPECT_FALSE(has(json, "\"NORMAL\""));
    std::remove(path.c_str());
}

TEST(Gltf, PointDataBecomesCustomAttributesAndColour) {
    Mesh m = cube_quads();
    m.AddPointData("temperature", f64({0, 1, 2, 3, 4, 5, 6, 7}));
    NDArray bad = f64({0, 1, 2, 3, 4, 5, 6, 7});
    bad.As<double>()[2] = std::numeric_limits<double>::quiet_NaN();
    m.AddPointData("bad", std::move(bad));
    const std::string path = mt::temp_path(".glb");
    write_gltf(path, m);
    std::string json = read_glb(path).mJson;
    EXPECT_TRUE(has(json, "\"_TEMPERATURE\""));
    EXPECT_TRUE(has(json, "\"name\":\"temperature\""));
    EXPECT_FALSE(has(json, "\"_BAD\""));
    EXPECT_FALSE(has(json, "COLOR_0"));

    GltfWriteOptions color;
    color.mColorBy = "temperature";
    write_gltf(path, m, color);
    json = read_glb(path).mJson;
    EXPECT_TRUE(has(json, "\"COLOR_0\""));
    EXPECT_TRUE(has(json, "\"extensionsUsed\":[\"KHR_materials_unlit\"]"));
    EXPECT_TRUE(has(json, "\"baseColorFactor\":[1,1,1,1]"));

    color.mUnlit = false;
    write_gltf(path, m, color);
    EXPECT_FALSE(has(read_glb(path).mJson, "KHR_materials_unlit"));

    GltfWriteOptions off;
    off.mFields = false;
    write_gltf(path, m, off);
    EXPECT_FALSE(has(read_glb(path).mJson, "\"_TEMPERATURE\""));
    std::remove(path.c_str());
}

TEST(Gltf, BadOptionsAreRefusedByName) {
    const std::string path = mt::temp_path(".glb");
    const Mesh m = cube_quads();
    GltfWriteOptions o;
    o.mSplitAngle = 200.0;
    EXPECT_THROW(write_gltf(path, m, o), std::invalid_argument);
    o = GltfWriteOptions{};
    o.mScale = 0.0;
    EXPECT_THROW(write_gltf(path, m, o), std::invalid_argument);
    o = GltfWriteOptions{};
    o.mColorBy = "nope";
    EXPECT_THROW(write_gltf(path, m, o), std::invalid_argument);
    o.mColorBy = "x";
    o.mCmap = "nope";
    EXPECT_THROW(write_gltf(path, m, o), std::invalid_argument);
    EXPECT_THROW(gltf_container_from_name("zip"), std::invalid_argument);
    EXPECT_THROW(gltf_up_axis_from_name("w"), std::invalid_argument);

    const double inf = std::numeric_limits<double>::infinity();
    const Mesh bad =
        make_mesh({{0, 0, 0}, {1, 0, 0}, {1, inf, 0}, {0, 1, 0}}, "quad", {{0, 1, 2, 3}});
    EXPECT_THROW(write_gltf(path, bad), WriteError);
    std::remove(path.c_str());
}

TEST(Gltf, AnEmptyMeshIsAValidEmptyScene) {
    const std::string path = mt::temp_path(".glb");
    Mesh empty;
    write_gltf(path, empty);
    const Glb glb = read_glb(path);
    EXPECT_TRUE(has(glb.mJson, "\"scenes\":[{}]"));
    EXPECT_TRUE(glb.mBin.empty());
    std::remove(path.c_str());
}

TEST(Gltf, TheSrgbTableIsTheInverseTransferFunction) {
    const std::uint32_t* bits = detail::srgb_to_linear_bits();
    auto value = [&](int i) {
        float f;
        std::memcpy(&f, &bits[i], sizeof(f));
        return f;
    };
    EXPECT_EQ(value(0), 0.0f);
    EXPECT_EQ(value(255), 1.0f);
    EXPECT_NEAR(value(128), 0.2158605f, 1e-7);
    for (int i = 1; i < 256; ++i)
        EXPECT_GT(value(i), value(i - 1));
}
