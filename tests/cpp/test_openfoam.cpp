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
#include <atomic>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/openfoam.hpp"
#include "meshioplusplus/region.hpp"

namespace fs = std::filesystem;

namespace {

// Write a minimal single-hex ASCII polyMesh case; return the case dir.
fs::path make_hex_case() {
    static std::atomic<unsigned> counter{0};
    fs::path base = fs::temp_directory_path() / ("meshio_of_" + std::to_string(counter++));
    fs::path poly = base / "constant" / "polyMesh";
    fs::create_directories(poly);

    auto hdr = [](const std::string& cls, const std::string& obj) {
        return "FoamFile\n{\n format ascii;\n class " + cls + ";\n object " + obj + ";\n}\n";
    };

    std::ofstream(poly / "points")
        << hdr("vectorField", "points")
        << "8\n(\n(0 0 0)\n(1 0 0)\n(1 1 0)\n(0 1 0)\n(0 0 1)\n(1 0 1)\n(1 1 "
           "1)\n(0 1 1)\n)\n";
    std::ofstream(poly / "faces")
        << hdr("faceList", "faces")
        << "6\n(\n4(0 3 2 1)\n4(4 5 6 7)\n4(0 1 5 4)\n4(2 3 7 6)\n4(1 2 6 "
           "5)\n4(0 4 7 3)\n)\n";
    std::ofstream(poly / "owner") << hdr("labelList", "owner") << "6\n(\n0\n0\n0\n0\n0\n0\n)\n";
    std::ofstream(poly / "boundary")
        << hdr("polyBoundaryMesh", "boundary")
        << "3\n(\nbottom { type wall; nFaces 1; startFace 0; }\ntop { type "
           "wall; nFaces 1; startFace 1; }\nsides { type wall; nFaces 4; "
           "startFace 2; }\n)\n";

    std::ofstream(base / "case.foam") << "";
    return base;
}

fs::path temp_case_dir() {
    static std::atomic<unsigned> counter{0};
    return fs::temp_directory_path() / ("meshio_ofw_" + std::to_string(counter++));
}

// ==========================================================================
// An INDEPENDENT polyMesh parser.
//
// Every assertion below reads the written bytes with this, never with
// `read_openfoam`. That is the point: a writer and its own reader share
// conventions, so a shared misconception is invisible to a round trip. This
// parser knows only what the OpenFOAM file format says.
// ==========================================================================

struct FoamPatchRead {
    std::string mName, mType;
    std::int64_t mNFaces = 0, mStartFace = 0;
};

struct PolyMeshRead {
    std::vector<std::array<double, 3>> mPoints;
    std::vector<std::vector<std::int64_t>> mFaces;
    std::vector<std::int64_t> mOwner, mNeighbour;
    std::vector<FoamPatchRead> mPatches;
};

// Drop the FoamFile header block and any comments, then return the body.
std::string of_body(const fs::path& rPath) {
    std::ifstream f(rPath);
    EXPECT_TRUE(f.good()) << "missing file " << rPath.string();
    std::string all((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    // strip /* */ and //
    std::string s;
    for (std::size_t i = 0; i < all.size();) {
        if (i + 1 < all.size() && all[i] == '/' && all[i + 1] == '*') {
            const std::size_t e = all.find("*/", i + 2);
            i = e == std::string::npos ? all.size() : e + 2;
        } else if (i + 1 < all.size() && all[i] == '/' && all[i + 1] == '/') {
            const std::size_t e = all.find('\n', i + 2);
            i = e == std::string::npos ? all.size() : e;
        } else {
            s.push_back(all[i++]);
        }
    }
    const std::size_t h = s.find("FoamFile");
    if (h == std::string::npos)
        return s;
    const std::size_t open = s.find('{', h);
    int depth = 0;
    for (std::size_t p = open; p < s.size(); ++p) {
        if (s[p] == '{')
            ++depth;
        else if (s[p] == '}' && --depth == 0)
            return s.substr(p + 1);
    }
    return s;
}

// The list header the OpenFOAM ASCII format requires: a count on a line of its
// own, then '(' on a line of its own. Returns the count and advances past '('.
std::int64_t of_list_header(std::istringstream& rSs) {
    std::string line;
    std::int64_t n = -1;
    while (std::getline(rSs, line)) {
        std::string t = line;
        t.erase(0, t.find_first_not_of(" \t\r"));
        if (!t.empty())
            t.erase(t.find_last_not_of(" \t\r") + 1);
        if (t.empty())
            continue;
        if (n < 0) {
            EXPECT_EQ(t.find_first_not_of("0123456789"), std::string::npos)
                << "the count must be on a line of ITS OWN, got: " << t;
            n = std::atoll(t.c_str());
            continue;
        }
        EXPECT_EQ(t, "(") << "'(' must be on a line of its own, got: " << t;
        return n;
    }
    ADD_FAILURE() << "no list header found";
    return 0;
}

PolyMeshRead read_polymesh(const fs::path& rPoly) {
    PolyMeshRead m;
    {
        std::istringstream ss(of_body(rPoly / "points"));
        const std::int64_t n = of_list_header(ss);
        std::string line;
        while (static_cast<std::int64_t>(m.mPoints.size()) < n && std::getline(ss, line)) {
            for (char& c : line)
                if (c == '(' || c == ')')
                    c = ' ';
            std::istringstream ls(line);
            double x, y, z;
            if (ls >> x >> y >> z)
                m.mPoints.push_back({x, y, z});
        }
        EXPECT_EQ(static_cast<std::int64_t>(m.mPoints.size()), n);
    }
    {
        std::istringstream ss(of_body(rPoly / "faces"));
        const std::int64_t n = of_list_header(ss);
        std::string line;
        while (static_cast<std::int64_t>(m.mFaces.size()) < n && std::getline(ss, line)) {
            const std::size_t lp = line.find('('), rp = line.find(')');
            if (lp == std::string::npos || rp == std::string::npos)
                continue;
            // The <count>( prefix is part of the format; check it agrees.
            const std::int64_t declared = std::atoll(line.c_str());
            std::istringstream ls(line.substr(lp + 1, rp - lp - 1));
            std::vector<std::int64_t> row;
            std::int64_t v;
            while (ls >> v)
                row.push_back(v);
            EXPECT_EQ(declared, static_cast<std::int64_t>(row.size()))
                << "face row's count prefix disagrees with its node count";
            m.mFaces.push_back(std::move(row));
        }
        EXPECT_EQ(static_cast<std::int64_t>(m.mFaces.size()), n);
    }
    auto read_labels = [&](const char* name, std::vector<std::int64_t>& out) {
        std::istringstream ss(of_body(rPoly / name));
        const std::int64_t n = of_list_header(ss);
        std::string line;
        while (static_cast<std::int64_t>(out.size()) < n && std::getline(ss, line)) {
            std::istringstream ls(line);
            std::int64_t v;
            while (ls >> v)
                out.push_back(v);
        }
        EXPECT_EQ(static_cast<std::int64_t>(out.size()), n) << name << " is short";
    };
    read_labels("owner", m.mOwner);
    read_labels("neighbour", m.mNeighbour);
    {
        const std::string body = of_body(rPoly / "boundary");
        std::size_t scan = 0;
        while (true) {
            const std::size_t brace = body.find('{', scan);
            if (brace == std::string::npos)
                break;
            // The patch name is the last token before its '{'.
            std::istringstream hs(body.substr(scan, brace - scan));
            std::string tok, name;
            while (hs >> tok)
                if (tok != "(")
                    name = tok;
            int depth = 0;
            std::size_t close = std::string::npos;
            for (std::size_t q = brace; q < body.size(); ++q) {
                if (body[q] == '{') {
                    ++depth;
                } else if (body[q] == '}' && --depth == 0) {
                    close = q;
                    break;
                }
            }
            if (close == std::string::npos)
                break;
            const std::string blk = body.substr(brace + 1, close - brace - 1);
            auto word = [&](const char* key) {
                const std::size_t k = blk.find(key);
                if (k == std::string::npos)
                    return std::string();
                std::istringstream bs(blk.substr(k + std::strlen(key)));
                std::string w;
                bs >> w;
                if (!w.empty() && w.back() == ';')
                    w.pop_back();
                return w;
            };
            FoamPatchRead pr;
            pr.mName = name;
            pr.mType = word("type");
            pr.mNFaces = std::atoll(word("nFaces").c_str());
            pr.mStartFace = std::atoll(word("startFace").c_str());
            m.mPatches.push_back(pr);
            scan = close + 1;
        }
    }
    return m;
}

// Newell area vector of a face, straight from the parsed file.
std::array<double, 3> of_area(const PolyMeshRead& rM, std::size_t Face) {
    std::array<double, 3> n{0, 0, 0};
    const auto& r = rM.mFaces[Face];
    for (std::size_t i = 0; i < r.size(); ++i) {
        const auto& a = rM.mPoints[static_cast<std::size_t>(r[i])];
        const auto& b = rM.mPoints[static_cast<std::size_t>(r[(i + 1) % r.size()])];
        n[0] += (a[1] - b[1]) * (a[2] + b[2]);
        n[1] += (a[2] - b[2]) * (a[0] + b[0]);
        n[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    for (int i = 0; i < 3; ++i)
        n[i] *= 0.5;
    return n;
}

std::int64_t of_num_cells(const PolyMeshRead& rM) {
    std::int64_t n = -1;
    for (std::int64_t o : rM.mOwner)
        n = std::max(n, o);
    for (std::int64_t o : rM.mNeighbour)
        n = std::max(n, o);
    return n + 1;
}

}  // namespace

TEST(OpenFoam, SingleHexAscii) {
    fs::path base = make_hex_case();
    meshioplusplus::OpenFoamInfo info;
    meshioplusplus::Mesh mesh = meshioplusplus::read_openfoam((base / "case.foam").string(), info);

    EXPECT_EQ(mesh.NumPoints(), 8u);
    bool has_hex = false;
    std::size_t nquad = 0;
    for (const auto cb : mesh.CellRange()) {
        if (cb.Type() == "hexahedron")
            has_hex = true;
        if (cb.Type() == "quad")
            nquad += cb.NumCells();
    }
    EXPECT_TRUE(has_hex);
    EXPECT_EQ(nquad, 6u);
    // 3 boundary patches -> 3 negative family tags.
    EXPECT_EQ(info.mCellTags.size(), 3u);

    std::error_code ec;
    fs::remove_all(base, ec);
}

TEST(OpenFoam, ResolveViaCaseDir) {
    fs::path base = make_hex_case();
    meshioplusplus::OpenFoamInfo info;
    // Pass the case directory itself (not the .foam file).
    meshioplusplus::Mesh mesh = meshioplusplus::read_openfoam(base.string(), info);
    EXPECT_EQ(mesh.NumPoints(), 8u);
    std::error_code ec;
    fs::remove_all(base, ec);
}

// Roadmap §1 tier B2: multi-region case selection.
namespace {

// A single-hex polyMesh written directly under `<base>/constant/<region>/`,
// mirroring `make_hex_case` but for a named region rather than the plain
// single-region layout.
void write_hex_region(const fs::path& base, const std::string& region) {
    fs::path poly = base / "constant" / region / "polyMesh";
    fs::create_directories(poly);
    auto hdr = [](const std::string& cls, const std::string& obj) {
        return "FoamFile\n{\n format ascii;\n class " + cls + ";\n object " + obj + ";\n}\n";
    };
    std::ofstream(poly / "points")
        << hdr("vectorField", "points")
        << "8\n(\n(0 0 0)\n(1 0 0)\n(1 1 0)\n(0 1 0)\n(0 0 1)\n(1 0 1)\n(1 1 "
           "1)\n(0 1 1)\n)\n";
    std::ofstream(poly / "faces")
        << hdr("faceList", "faces")
        << "6\n(\n4(0 3 2 1)\n4(4 5 6 7)\n4(0 1 5 4)\n4(2 3 7 6)\n4(1 2 6 "
           "5)\n4(0 4 7 3)\n)\n";
    std::ofstream(poly / "owner") << hdr("labelList", "owner") << "6\n(\n0\n0\n0\n0\n0\n0\n)\n";
    std::ofstream(poly / "boundary")
        << hdr("polyBoundaryMesh", "boundary")
        << "1\n(\nallB { type wall; nFaces 6; startFace 0; }\n)\n";
}

fs::path make_multi_region_case() {
    static std::atomic<unsigned> counter{0};
    fs::path base = fs::temp_directory_path() / ("meshio_of_mr_" + std::to_string(counter++));
    write_hex_region(base, "fluid");
    write_hex_region(base, "solid");
    std::ofstream(base / "constant" / "regionProperties") << "FoamFile\n{\n}\nregions\n(\n);\n";
    return base;
}

}  // namespace

TEST(OpenFoam, MultiRegionWithNoRegionSetThrowsNamingTheRegions) {
    const fs::path base = make_multi_region_case();
    meshioplusplus::OpenFoamInfo info;
    bool threw = false;
    try {
        meshioplusplus::read_openfoam(base.string(), info);
    } catch (const meshioplusplus::ReadError& e) {
        threw = true;
        EXPECT_NE(std::string(e.what()).find("fluid"), std::string::npos);
        EXPECT_NE(std::string(e.what()).find("solid"), std::string::npos);
    }
    EXPECT_TRUE(threw);
    std::error_code ec;
    fs::remove_all(base, ec);
}

TEST(OpenFoam, MultiRegionSelectsTheNamedRegion) {
    const fs::path base = make_multi_region_case();
    meshioplusplus::OpenFoamInfo info;
    info.mRegion = "fluid";
    const meshioplusplus::Mesh mesh = meshioplusplus::read_openfoam(base.string(), info);
    EXPECT_EQ(mesh.NumPoints(), 8u);
    std::error_code ec;
    fs::remove_all(base, ec);
}

TEST(OpenFoam, MultiRegionAnUnknownRegionThrowsNamingIt) {
    const fs::path base = make_multi_region_case();
    meshioplusplus::OpenFoamInfo info;
    info.mRegion = "nope";
    EXPECT_THROW(meshioplusplus::read_openfoam(base.string(), info), meshioplusplus::ReadError);
    std::error_code ec;
    fs::remove_all(base, ec);
}

// The `polyMesh`-directory rule already covers this without `mRegion` at all.
TEST(OpenFoam, MultiRegionDirectPolyMeshPathNeedsNoRegionField) {
    const fs::path base = make_multi_region_case();
    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh mesh =
        meshioplusplus::read_openfoam((base / "constant" / "solid" / "polyMesh").string(), info);
    EXPECT_EQ(mesh.NumPoints(), 8u);
    std::error_code ec;
    fs::remove_all(base, ec);
}

// ==========================================================================
//                              WRITER
// ==========================================================================

namespace {

// N x N x N unit hexahedra.
meshioplusplus::Mesh hex_grid(int N) {
    meshioplusplus::Mesh m;
    std::vector<std::vector<double>> pts;
    const int P = N + 1;
    for (int i = 0; i < P; ++i)
        for (int j = 0; j < P; ++j)
            for (int k = 0; k < P; ++k)
                pts.push_back({double(i), double(j), double(k)});
    auto id = [&](int i, int j, int k) { return std::int64_t((i * P + j) * P + k); };
    std::vector<std::vector<std::int64_t>> rows;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            for (int k = 0; k < N; ++k)
                rows.push_back({id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k),
                                id(i, j + 1, k), id(i, j, k + 1), id(i + 1, j, k + 1),
                                id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)});
    m.AssignPoints(mt::points_from(pts));
    m.AddCellBlock("hexahedron", mt::conn_from(rows));
    return m;
}

fs::path write_case(const meshioplusplus::Mesh& rMesh,
                    const meshioplusplus::OpenFoamInfo& rInfo = {}) {
    const fs::path base = temp_case_dir();
    meshioplusplus::write_openfoam((base / "case.foam").string(), rMesh, rInfo);
    return base;
}

meshioplusplus::NDArray i64(const std::vector<std::int64_t>& rVals) {
    meshioplusplus::NDArray a = meshioplusplus::NDArray::Uninit(meshioplusplus::DType::Int64,
                                                                 {rVals.size()});
    for (std::size_t i = 0; i < rVals.size(); ++i)
        a.As<std::int64_t>()[i] = rVals[i];
    return a;
}

meshioplusplus::NDArray i64_pairs(const std::vector<std::int64_t>& rFlat) {
    meshioplusplus::NDArray a =
        meshioplusplus::NDArray::Uninit(meshioplusplus::DType::Int64, {rFlat.size() / 2, 2});
    for (std::size_t i = 0; i < rFlat.size(); ++i)
        a.As<std::int64_t>()[i] = rFlat[i];
    return a;
}

}  // namespace

// --- C0: the neighbour file --------------------------------------------
// OpenFOAM's `neighbour` holds ONLY internal faces. Our own reader also
// accepts a -1-padded full-length list, which is precisely why a round trip
// through it cannot be the oracle for this writer.
TEST(OpenFoamWrite, NeighbourFileHoldsOnlyInternalFaces) {
    const fs::path base = write_case(hex_grid(2));
    const PolyMeshRead m = read_polymesh(base / "constant" / "polyMesh");

    EXPECT_LT(m.mNeighbour.size(), m.mOwner.size()) << "neighbour must be shorter than owner";
    EXPECT_EQ(m.mOwner.size(), m.mFaces.size());
    for (std::int64_t v : m.mNeighbour)
        EXPECT_GE(v, 0) << "neighbour must not be -1-padded";
    std::error_code ec;
    fs::remove_all(base, ec);
}

// --- C1 / C3: the upper-triangular ordering ----------------------------
// `cell_faces(Hexahedron)`'s rows visit cell 0's neighbours in a genuinely
// non-monotone order, so a 2x2x2 grid fails this unless the writer sorts.
TEST(OpenFoamWrite, InternalFacesAreUpperTriangularAndSorted) {
    const fs::path base = write_case(hex_grid(2));
    const PolyMeshRead m = read_polymesh(base / "constant" / "polyMesh");

    ASSERT_FALSE(m.mNeighbour.empty());
    for (std::size_t i = 0; i < m.mNeighbour.size(); ++i) {
        EXPECT_LT(m.mOwner[i], m.mNeighbour[i]) << "face " << i << " has owner >= neighbour";
        if (i == 0)
            continue;
        const bool ok = m.mOwner[i - 1] < m.mOwner[i] ||
                        (m.mOwner[i - 1] == m.mOwner[i] && m.mNeighbour[i - 1] <= m.mNeighbour[i]);
        EXPECT_TRUE(ok) << "faces " << i - 1 << "," << i << " are out of order: ("
                        << m.mOwner[i - 1] << "," << m.mNeighbour[i - 1] << ") then ("
                        << m.mOwner[i] << "," << m.mNeighbour[i] << ")";
    }
    std::error_code ec;
    fs::remove_all(base, ec);
}

// --- C2: internal faces come first -------------------------------------
TEST(OpenFoamWrite, InternalFacesComeFirst) {
    const fs::path base = write_case(hex_grid(2));
    const PolyMeshRead m = read_polymesh(base / "constant" / "polyMesh");

    ASSERT_FALSE(m.mPatches.empty());
    std::int64_t first_boundary = m.mPatches.front().mStartFace;
    EXPECT_EQ(first_boundary, static_cast<std::int64_t>(m.mNeighbour.size()))
        << "the boundary range must start exactly where the internal one ends";
    std::error_code ec;
    fs::remove_all(base, ec);
}

// --- C4: normals point owner -> neighbour ------------------------------
TEST(OpenFoamWrite, InternalNormalsPointOwnerToNeighbour) {
    const fs::path base = write_case(hex_grid(2));
    const PolyMeshRead m = read_polymesh(base / "constant" / "polyMesh");

    // Cell centroids from the written topology alone.
    const std::int64_t nc = of_num_cells(m);
    std::vector<std::array<double, 3>> c(nc, {0, 0, 0});
    std::vector<double> w(nc, 0.0);
    for (std::size_t f = 0; f < m.mFaces.size(); ++f) {
        for (std::int64_t nid : m.mFaces[f]) {
            const auto& p = m.mPoints[static_cast<std::size_t>(nid)];
            for (int side = 0; side < 2; ++side) {
                const std::int64_t cell =
                    side == 0 ? m.mOwner[f]
                              : (f < m.mNeighbour.size() ? m.mNeighbour[f] : -1);
                if (cell < 0)
                    continue;
                for (int a = 0; a < 3; ++a)
                    c[static_cast<std::size_t>(cell)][a] += p[a];
                w[static_cast<std::size_t>(cell)] += 1.0;
            }
        }
    }
    for (std::int64_t i = 0; i < nc; ++i)
        for (int a = 0; a < 3; ++a)
            c[static_cast<std::size_t>(i)][a] /= w[static_cast<std::size_t>(i)];

    for (std::size_t f = 0; f < m.mNeighbour.size(); ++f) {
        const auto n = of_area(m, f);
        const auto& co = c[static_cast<std::size_t>(m.mOwner[f])];
        const auto& cn = c[static_cast<std::size_t>(m.mNeighbour[f])];
        const double d =
            n[0] * (cn[0] - co[0]) + n[1] * (cn[1] - co[1]) + n[2] * (cn[2] - co[2]);
        EXPECT_GT(d, 0.0) << "internal face " << f << " does not point owner->neighbour";
    }
    std::error_code ec;
    fs::remove_all(base, ec);
}

// --- C6: patches partition the boundary range exactly ------------------
// A patch spanning BOTH a triangle and a quad block is what breaks a writer
// that emits boundary faces in cell-block order.
TEST(OpenFoamWrite, BoundaryFacesAreContiguousPerPatch) {
    // One hex + one pyramid glued on top, so the boundary has both quads and
    // triangles, and we tag them into two interleaved patches.
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1},
                                    {0.5, 0.5, 2}}));
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}}));
    m.AddCellBlock("pyramid", mt::conn_from({{4, 5, 6, 7, 8}}));
    // boundary blocks: two triangles of the pyramid in patch A, one quad of the
    // hex in patch A too, and another quad in patch B.
    m.AddCellBlock("triangle", mt::conn_from({{4, 5, 8}, {6, 7, 8}}));
    m.AddCellBlock("quad", mt::conn_from({{0, 3, 2, 1}, {0, 1, 5, 4}}));
    m.AddCellData("cell_tags", {mt::int_data_array({0}), mt::int_data_array({0}),
                                mt::int_data_array({-1, -2}), mt::int_data_array({-1, -2})});

    meshioplusplus::OpenFoamInfo info;
    info.mCellTags[-1] = {"alpha"};
    info.mCellTags[-2] = {"beta"};
    info.mPatchTypes[-1] = "wall";

    const fs::path base = write_case(m, info);
    const PolyMeshRead r = read_polymesh(base / "constant" / "polyMesh");

    ASSERT_GE(r.mPatches.size(), 2u);
    EXPECT_EQ(r.mPatches[0].mName, "alpha");
    EXPECT_EQ(r.mPatches[0].mType, "wall");
    EXPECT_EQ(r.mPatches[1].mName, "beta");

    std::int64_t expect = static_cast<std::int64_t>(r.mNeighbour.size());
    for (const auto& p : r.mPatches) {
        EXPECT_EQ(p.mStartFace, expect) << "patch '" << p.mName << "' is not contiguous";
        EXPECT_GT(p.mNFaces, 0);
        expect += p.mNFaces;
    }
    EXPECT_EQ(expect, static_cast<std::int64_t>(r.mFaces.size()))
        << "patches do not cover the whole boundary range";
    std::error_code ec;
    fs::remove_all(base, ec);
}

// --- The strongest oracle available without OpenFOAM -------------------
// Computed from the written files alone, sharing no code with the writer's
// ordering logic: every cell's outward face areas must sum to zero, and its
// divergence-theorem volume must be the true one.
TEST(OpenFoamWrite, EveryWrittenCellIsClosedAndHasTheRightVolume) {
    meshioplusplus::Mesh m = hex_grid(2);
    const fs::path base = write_case(m);
    const PolyMeshRead r = read_polymesh(base / "constant" / "polyMesh");

    const std::int64_t nc = of_num_cells(r);
    ASSERT_EQ(nc, 8);
    std::vector<std::array<double, 3>> sum(nc, {0, 0, 0});
    std::vector<double> scale(nc, 0.0), vol(nc, 0.0);

    for (std::size_t f = 0; f < r.mFaces.size(); ++f) {
        const auto a = of_area(r, f);
        const double mag = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
        std::array<double, 3> fc{0, 0, 0};
        for (std::int64_t nid : r.mFaces[f])
            for (int k = 0; k < 3; ++k)
                fc[k] += r.mPoints[static_cast<std::size_t>(nid)][k] /
                         static_cast<double>(r.mFaces[f].size());
        // owner: outward (+). neighbour: the same face is inward for it (-).
        const std::int64_t own = r.mOwner[f];
        for (int k = 0; k < 3; ++k)
            sum[static_cast<std::size_t>(own)][k] += a[k];
        scale[static_cast<std::size_t>(own)] += mag;
        vol[static_cast<std::size_t>(own)] += (fc[0] * a[0] + fc[1] * a[1] + fc[2] * a[2]) / 3.0;
        if (f < r.mNeighbour.size()) {
            const std::int64_t nb = r.mNeighbour[f];
            for (int k = 0; k < 3; ++k)
                sum[static_cast<std::size_t>(nb)][k] -= a[k];
            scale[static_cast<std::size_t>(nb)] += mag;
            vol[static_cast<std::size_t>(nb)] -=
                (fc[0] * a[0] + fc[1] * a[1] + fc[2] * a[2]) / 3.0;
        }
    }
    for (std::int64_t i = 0; i < nc; ++i) {
        const auto& s = sum[static_cast<std::size_t>(i)];
        const double mag = std::sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
        EXPECT_LT(mag, 1e-9 * scale[static_cast<std::size_t>(i)])
            << "cell " << i << " is not closed";
        EXPECT_NEAR(vol[static_cast<std::size_t>(i)], 1.0, 1e-12)
            << "cell " << i << " has the wrong volume (or is inside out)";
    }
    std::error_code ec;
    fs::remove_all(base, ec);
}

// A polyhedron mesh is what this format exists for.
TEST(OpenFoamWrite, APolyhedronCellRoundTripsThroughTheWrittenFiles) {
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1}}));
    m.AddPolyhedronBlock("polyhedron8", {{{
                                             {0, 3, 2, 1},
                                             {4, 5, 6, 7},
                                             {0, 1, 5, 4},
                                             {2, 3, 7, 6},
                                             {1, 2, 6, 5},
                                             {0, 4, 7, 3},
                                         }}});
    const fs::path base = write_case(m);
    const PolyMeshRead r = read_polymesh(base / "constant" / "polyMesh");

    EXPECT_EQ(of_num_cells(r), 1);
    EXPECT_EQ(r.mFaces.size(), 6u);
    EXPECT_TRUE(r.mNeighbour.empty());

    double vol = 0.0;
    for (std::size_t f = 0; f < r.mFaces.size(); ++f) {
        const auto a = of_area(r, f);
        std::array<double, 3> fc{0, 0, 0};
        for (std::int64_t nid : r.mFaces[f])
            for (int k = 0; k < 3; ++k)
                fc[k] += r.mPoints[static_cast<std::size_t>(nid)][k] /
                         static_cast<double>(r.mFaces[f].size());
        vol += (fc[0] * a[0] + fc[1] * a[1] + fc[2] * a[2]) / 3.0;
    }
    EXPECT_NEAR(vol, 1.0, 1e-12);
    std::error_code ec;
    fs::remove_all(base, ec);
}

// A mesh from any other format carries no tags. It must still produce a
// loadable single-patch case -- blockMesh's own `defaultFaces` -- not an error
// and not patches invented from geometry.
TEST(OpenFoamWrite, AnUntaggedMeshGetsOneDefaultFacesPatch) {
    const fs::path base = write_case(hex_grid(1));
    const PolyMeshRead r = read_polymesh(base / "constant" / "polyMesh");

    ASSERT_EQ(r.mPatches.size(), 1u);
    EXPECT_EQ(r.mPatches[0].mName, "defaultFaces");
    EXPECT_EQ(r.mPatches[0].mType, "patch");
    EXPECT_EQ(r.mPatches[0].mNFaces, 6);
    EXPECT_EQ(r.mPatches[0].mStartFace, 0);
    std::error_code ec;
    fs::remove_all(base, ec);
}

// A type needing companion dictionary entries we cannot carry must be
// downgraded, not passed through: OpenFOAM refuses to LOAD such a case.
TEST(OpenFoamWrite, UnsafePatchTypeIsDowngradedToPatch) {
    meshioplusplus::Mesh m = hex_grid(1);
    m.AddCellBlock("quad", mt::conn_from({{0, 2, 3, 1}}));
    m.AddCellData("cell_tags", {mt::int_data_array({0}), mt::int_data_array({-1})});
    meshioplusplus::OpenFoamInfo info;
    info.mCellTags[-1] = {"inlet"};
    info.mPatchTypes[-1] = "cyclicAMI";

    const fs::path base = write_case(m, info);
    const PolyMeshRead r = read_polymesh(base / "constant" / "polyMesh");
    bool found = false;
    for (const auto& p : r.mPatches)
        if (p.mName == "inlet") {
            found = true;
            EXPECT_EQ(p.mType, "patch") << "cyclicAMI was written without its companion keys";
        }
    EXPECT_TRUE(found);
    std::error_code ec;
    fs::remove_all(base, ec);
}

// An inverted cell must be written correctly oriented, not rejected.
TEST(OpenFoamWrite, AnInvertedCellIsWrittenOutward) {
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1},
                                    {0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0}}));
    m.AddCellBlock("hexahedron", mt::conn_from({{0, 1, 2, 3, 4, 5, 6, 7}}));
    const fs::path base = write_case(m);
    const PolyMeshRead r = read_polymesh(base / "constant" / "polyMesh");

    double vol = 0.0;
    for (std::size_t f = 0; f < r.mFaces.size(); ++f) {
        const auto a = of_area(r, f);
        std::array<double, 3> fc{0, 0, 0};
        for (std::int64_t nid : r.mFaces[f])
            for (int k = 0; k < 3; ++k)
                fc[k] += r.mPoints[static_cast<std::size_t>(nid)][k] /
                         static_cast<double>(r.mFaces[f].size());
        vol += (fc[0] * a[0] + fc[1] * a[1] + fc[2] * a[2]) / 3.0;
    }
    EXPECT_NEAR(vol, 1.0, 1e-12) << "an inverted cell was written with inward normals";
    std::error_code ec;
    fs::remove_all(base, ec);
}

// A stale `neighbour` from a bigger previous case is one of the nastiest ways
// to corrupt a polyMesh, so the file is always written, even empty.
TEST(OpenFoamWrite, StalePolyMeshCompanionsAreReplacedOrRemoved) {
    const fs::path base = temp_case_dir();
    const fs::path poly = base / "constant" / "polyMesh";
    fs::create_directories(poly);
    std::ofstream(poly / "neighbour") << "GARBAGE FROM A PREVIOUS CASE\n";
    std::ofstream(poly / "cellZones") << "stale\n";

    // A single cell has NO internal faces at all -- the case where a writer is
    // most tempted to skip the file.
    meshioplusplus::write_openfoam((base / "case.foam").string(), hex_grid(1), {});

    const PolyMeshRead r = read_polymesh(poly);
    EXPECT_TRUE(r.mNeighbour.empty());
    EXPECT_FALSE(fs::exists(poly / "cellZones")) << "a stale companion file was left behind";
    EXPECT_TRUE(fs::exists(base / "case.foam")) << "the .foam marker was not written";
    std::error_code ec;
    fs::remove_all(base, ec);
}

// The weak oracle, kept deliberately small and labelled: the reader shares this
// writer's conventions, so a shared misconception is invisible here.
TEST(OpenFoamWrite, RoundTripsThroughOurOwnReader) {
    meshioplusplus::Mesh m = hex_grid(2);
    const fs::path base = write_case(m);

    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh back =
        meshioplusplus::read_openfoam((base / "case.foam").string(), info);
    EXPECT_EQ(back.NumPoints(), 27u);
    std::size_t nhex = 0;
    for (const auto cb : back.CellRange())
        if (cb.Type() == "hexahedron")
            nhex += cb.NumCells();
    EXPECT_EQ(nhex, 8u);
    std::error_code ec;
    fs::remove_all(base, ec);
}

// Roadmap §1 tier B2: cellZones/faceZones/pointZones as named Regions.
TEST(OpenFoamWrite, ZonesRoundTripAsNamedRegions) {
    meshioplusplus::Mesh m = hex_grid(2);  // 8 hex cells, 27 points
    m.AddRegion(meshioplusplus::Region("core", meshioplusplus::RegionKind::Cell, i64({0, 5})));
    m.AddRegion(meshioplusplus::Region("corners", meshioplusplus::RegionKind::Point, i64({0, 26})));
    // Cell 0's facet 0 (`cell_faces(Hexahedron)[0]` = local corners 0,4,7,3) is
    // its x=0 face -- a genuine boundary face for cell (0,0,0), so its owner is
    // unambiguous. Point ids 0,1,3,4 are that face's corners in hex_grid's own
    // numbering (`id(i,j,k) = (i*P+j)*P+k`, P=3): (0,0,0)=0, (0,0,1)=1,
    // (0,1,0)=3, (0,1,1)=4.
    m.AddRegion(
        meshioplusplus::Region("inlet", meshioplusplus::RegionKind::Side, i64_pairs({0, 0})));

    const fs::path base = write_case(m);
    EXPECT_TRUE(fs::exists(base / "constant" / "polyMesh" / "cellZones"));
    EXPECT_TRUE(fs::exists(base / "constant" / "polyMesh" / "pointZones"));
    EXPECT_TRUE(fs::exists(base / "constant" / "polyMesh" / "faceZones"));

    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh back =
        meshioplusplus::read_openfoam((base / "case.foam").string(), info);

    ASSERT_NE(back.FindRegion("core", meshioplusplus::RegionKind::Cell), meshioplusplus::Mesh::npos);
    const meshioplusplus::Region& core =
        back.Region(back.FindRegion("core", meshioplusplus::RegionKind::Cell));
    ASSERT_EQ(core.NumEntries(), 2u);
    EXPECT_EQ(core.Entries()[0], 0);
    EXPECT_EQ(core.Entries()[1], 5);

    ASSERT_NE(back.FindRegion("corners", meshioplusplus::RegionKind::Point),
             meshioplusplus::Mesh::npos);
    const meshioplusplus::Region& corners =
        back.Region(back.FindRegion("corners", meshioplusplus::RegionKind::Point));
    ASSERT_EQ(corners.NumEntries(), 2u);
    EXPECT_EQ(corners.Entries()[0], 0);
    EXPECT_EQ(corners.Entries()[1], 26);

    ASSERT_NE(back.FindRegion("inlet", meshioplusplus::RegionKind::Side), meshioplusplus::Mesh::npos);
    const meshioplusplus::Region& inlet =
        back.Region(back.FindRegion("inlet", meshioplusplus::RegionKind::Side));
    ASSERT_EQ(inlet.NumEntries(), 1u);
    EXPECT_EQ(inlet.Entries()[0], 0);  // global cell (single block, unambiguous)
    // The facet index itself need not survive numerically -- the reader's
    // hexahedron reconstruction is free to renumber local nodes -- but it must
    // still name the SAME geometric face: point ids {0, 1, 3, 4}.
    const std::int64_t facet = inlet.Entries()[1];
    ASSERT_GE(facet, 0);
    const auto back_hex = back.Cells(0);
    ASSERT_EQ(std::string(back_hex.Type()), "hexahedron");
    const auto& facedefs = meshioplusplus::detail::cell_faces(meshioplusplus::CellType::Hexahedron);
    ASSERT_LT(static_cast<std::size_t>(facet), facedefs.size());
    std::set<std::int64_t> face_pts;
    for (int c = 0; c < facedefs[static_cast<std::size_t>(facet)].mNumCorners; ++c)
        face_pts.insert(
            back_hex.Conn().As<std::int64_t>()[facedefs[static_cast<std::size_t>(facet)].mNodes[c]]);
    EXPECT_EQ(face_pts, (std::set<std::int64_t>{0, 1, 3, 4}));

    std::error_code ec;
    fs::remove_all(base, ec);
}

// Roadmap §1.1: binary write, all four label/scalar width combinations.
// Every combination must bit-exactly round-trip through the reader, which
// already had binary support for every width (only the writer was missing).
class OpenFoamBinaryWidths
    : public ::testing::TestWithParam<std::pair<int, int>> {};

TEST_P(OpenFoamBinaryWidths, RoundTripsBitExactly) {
    const auto [label_bits, scalar_bits] = GetParam();
    meshioplusplus::Mesh m = hex_grid(2);  // 8 hex cells, 27 points
    m.AddRegion(meshioplusplus::Region("core", meshioplusplus::RegionKind::Cell, i64({0, 5})));

    const fs::path base = temp_case_dir();
    meshioplusplus::OpenFoamWriteOptions wopts;
    wopts.mBinary = true;
    wopts.mLabelBits = label_bits;
    wopts.mScalarBits = scalar_bits;
    meshioplusplus::write_openfoam((base / "case.foam").string(), m, meshioplusplus::OpenFoamInfo{},
                                   wopts);

    // The header must actually say so, or the round trip below would silently
    // exercise the ASCII reader path instead of the binary one.
    {
        std::ifstream pf(base / "constant" / "polyMesh" / "points");
        std::string text((std::istreambuf_iterator<char>(pf)), std::istreambuf_iterator<char>());
        EXPECT_NE(text.find("format      binary;"), std::string::npos);
        EXPECT_NE(text.find("label=" + std::to_string(label_bits)), std::string::npos);
        EXPECT_NE(text.find("scalar=" + std::to_string(scalar_bits)), std::string::npos);
    }

    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh back =
        meshioplusplus::read_openfoam((base / "case.foam").string(), info);
    EXPECT_EQ(back.NumPoints(), 27u);
    std::size_t nhex = 0;
    for (const auto cb : back.CellRange())
        if (cb.Type() == "hexahedron")
            nhex += cb.NumCells();
    EXPECT_EQ(nhex, 8u);

    // Points must match to the written scalar width's own precision -- exact
    // for 64-bit, float-rounded for 32-bit (hex_grid's coordinates are small
    // integers, which float32 represents exactly, so this stays an exact
    // comparison either way rather than needing a tolerance).
    ASSERT_EQ(back.NumPoints(), m.NumPoints());
    for (std::size_t i = 0; i < back.NumPoints() * 3; ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(back.Points(), i),
                         meshioplusplus::detail::read_double(m.Points(), i));

    ASSERT_NE(back.FindRegion("core", meshioplusplus::RegionKind::Cell), meshioplusplus::Mesh::npos);
    const meshioplusplus::Region& core =
        back.Region(back.FindRegion("core", meshioplusplus::RegionKind::Cell));
    ASSERT_EQ(core.NumEntries(), 2u);
    EXPECT_EQ(core.Entries()[0], 0);
    EXPECT_EQ(core.Entries()[1], 5);

    std::error_code ec;
    fs::remove_all(base, ec);
}

INSTANTIATE_TEST_SUITE_P(WidthCombinations, OpenFoamBinaryWidths,
                        ::testing::Values(std::pair{32, 64}, std::pair{32, 32},
                                          std::pair{64, 64}, std::pair{64, 32}));

TEST(OpenFoamWrite, BinaryAndAsciiWriteOfTheSameMeshReadToIdenticalMeshes) {
    meshioplusplus::Mesh m = hex_grid(2);
    m.AddRegion(meshioplusplus::Region("core", meshioplusplus::RegionKind::Cell, i64({0, 5})));

    const fs::path ascii_base = write_case(m);
    const fs::path binary_base = temp_case_dir();
    meshioplusplus::OpenFoamWriteOptions wopts;
    wopts.mBinary = true;
    meshioplusplus::write_openfoam((binary_base / "case.foam").string(), m,
                                   meshioplusplus::OpenFoamInfo{}, wopts);

    meshioplusplus::OpenFoamInfo ai, bi;
    const meshioplusplus::Mesh a =
        meshioplusplus::read_openfoam((ascii_base / "case.foam").string(), ai);
    const meshioplusplus::Mesh b =
        meshioplusplus::read_openfoam((binary_base / "case.foam").string(), bi);

    EXPECT_EQ(a.NumPoints(), b.NumPoints());
    EXPECT_EQ(a.NumCellBlocks(), b.NumCellBlocks());
    for (std::size_t i = 0; i < a.NumPoints() * 3; ++i)
        EXPECT_DOUBLE_EQ(meshioplusplus::detail::read_double(a.Points(), i),
                         meshioplusplus::detail::read_double(b.Points(), i));
    for (std::size_t bi_ = 0; bi_ < a.NumCellBlocks(); ++bi_)
        for (std::size_t i = 0; i < a.Cells(bi_).Conn().Size(); ++i)
            EXPECT_EQ(a.Cells(bi_).Conn().As<std::int64_t>()[i],
                     b.Cells(bi_).Conn().As<std::int64_t>()[i]);

    std::error_code ec;
    fs::remove_all(ascii_base, ec);
    fs::remove_all(binary_base, ec);
}

TEST(OpenFoamWrite, BinaryPolyhedronRoundTrips) {
    meshioplusplus::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0},
                                    {1, 0, 0},
                                    {1, 1, 0},
                                    {0, 1, 0},
                                    {0, 0, 1},
                                    {1, 0, 1},
                                    {1, 1, 1},
                                    {0, 1, 1}}));
    // A cube with its top face split into two triangles: 7 faces over 8
    // points, which matches none of the reader's known (n_faces, n_points)
    // cell signatures (tetra 4/4, pyramid 5/5, wedge 5/6, hexahedron 6/8) --
    // unlike a plain 6-quad-face cube, which the reader would reclassify
    // back into a hexahedron. This keeps the ragged nfaced binary path
    // genuinely exercised end to end, not just decoded and then merged away.
    m.AddPolyhedronBlock("polyhedron7", {{{0, 3, 2, 1},
                                          {0, 1, 5, 4},
                                          {1, 2, 6, 5},
                                          {2, 3, 7, 6},
                                          {3, 0, 4, 7},
                                          {4, 5, 6},
                                          {4, 6, 7}}});

    const fs::path base = temp_case_dir();
    meshioplusplus::OpenFoamWriteOptions wopts;
    wopts.mBinary = true;
    meshioplusplus::write_openfoam((base / "case.foam").string(), m, meshioplusplus::OpenFoamInfo{},
                                   wopts);

    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh back =
        meshioplusplus::read_openfoam((base / "case.foam").string(), info);
    EXPECT_EQ(back.NumPoints(), 8u);
    // With no cell to be a neighbour, all 7 faces are also boundary faces,
    // so they additionally come back as five one-cell `quad` and two
    // one-cell `triangle` blocks (the single, untagged `defaultFaces` patch).
    std::size_t nvol = 0, nquad = 0, ntri = 0;
    for (const auto cb : back.CellRange()) {
        if (cb.IsPolyhedron())
            nvol += cb.NumCells();
        else if (cb.Type() == "quad")
            nquad += cb.NumCells();
        else if (cb.Type() == "triangle")
            ntri += cb.NumCells();
    }
    EXPECT_EQ(nvol, 1u);
    EXPECT_EQ(nquad, 5u);
    EXPECT_EQ(ntri, 2u);

    std::error_code ec;
    fs::remove_all(base, ec);
}

TEST(OpenFoamWrite, BinaryZonesRoundTripAsNamedRegions) {
    // Same fixture as ZonesRoundTripAsNamedRegions, but binary -- exercises
    // parse_zone_file_binary specifically (cellZones/pointZones/faceZones'
    // List<label> payload is raw bytes, everything else around it stays text).
    meshioplusplus::Mesh m = hex_grid(2);
    m.AddRegion(meshioplusplus::Region("core", meshioplusplus::RegionKind::Cell, i64({0, 5})));
    m.AddRegion(meshioplusplus::Region("corners", meshioplusplus::RegionKind::Point, i64({0, 26})));
    m.AddRegion(
        meshioplusplus::Region("inlet", meshioplusplus::RegionKind::Side, i64_pairs({0, 0})));

    const fs::path base = temp_case_dir();
    meshioplusplus::OpenFoamWriteOptions wopts;
    wopts.mBinary = true;
    meshioplusplus::write_openfoam((base / "case.foam").string(), m, meshioplusplus::OpenFoamInfo{},
                                   wopts);
    EXPECT_TRUE(fs::exists(base / "constant" / "polyMesh" / "cellZones"));
    EXPECT_TRUE(fs::exists(base / "constant" / "polyMesh" / "pointZones"));
    EXPECT_TRUE(fs::exists(base / "constant" / "polyMesh" / "faceZones"));

    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh back =
        meshioplusplus::read_openfoam((base / "case.foam").string(), info);

    ASSERT_NE(back.FindRegion("core", meshioplusplus::RegionKind::Cell), meshioplusplus::Mesh::npos);
    const meshioplusplus::Region& core =
        back.Region(back.FindRegion("core", meshioplusplus::RegionKind::Cell));
    ASSERT_EQ(core.NumEntries(), 2u);
    EXPECT_EQ(core.Entries()[0], 0);
    EXPECT_EQ(core.Entries()[1], 5);

    ASSERT_NE(back.FindRegion("corners", meshioplusplus::RegionKind::Point),
             meshioplusplus::Mesh::npos);
    const meshioplusplus::Region& corners =
        back.Region(back.FindRegion("corners", meshioplusplus::RegionKind::Point));
    ASSERT_EQ(corners.NumEntries(), 2u);
    EXPECT_EQ(corners.Entries()[0], 0);
    EXPECT_EQ(corners.Entries()[1], 26);

    ASSERT_NE(back.FindRegion("inlet", meshioplusplus::RegionKind::Side), meshioplusplus::Mesh::npos);
    const meshioplusplus::Region& inlet =
        back.Region(back.FindRegion("inlet", meshioplusplus::RegionKind::Side));
    ASSERT_EQ(inlet.NumEntries(), 1u);
    EXPECT_EQ(inlet.Entries()[0], 0);

    std::error_code ec;
    fs::remove_all(base, ec);
}

// Roadmap §1 tier B2: reconstructing a decomposed (`processor*/`) case.
//
// Two unit hexahedra sharing one face (cube A: x in [0,1], cube B: x in
// [1,2]), decomposed by hand into `processor0`/`processor1` -- one cell each
// -- with correctly-derived `*ProcAddressing` files. `write_openfoam` writes
// each processor's own isolated-cell polyMesh (6 boundary faces, no internal
// ones -- a single cell has no neighbour of its own); this test's job is
// purely the *ProcAddressing bookkeeping and the reconstruction that reads
// it back into ONE mesh with the shared face restored as a genuine internal
// face. The point/cell addressing is chosen so processor0's local ids ARE
// the global ids (an identity map is still a real exercise of the address
// lookup, just not a permutation of it); processor1's is a real offset/shift,
// so both addressing shapes are covered.
TEST(OpenFoamDecompose, ReconstructsTwoProcessorsSharingOneFace) {
    const fs::path base = temp_case_dir();

    // 12 global points: cube A's own 4 (x=0), the shared 4 (x=1), cube B's
    // own 4 (x=2).
    const std::vector<std::vector<double>> g = {
        {0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1},  // 0-3: A only
        {1, 0, 0}, {1, 1, 0}, {1, 1, 1}, {1, 0, 1},  // 4-7: shared (x=1)
        {2, 0, 0}, {2, 1, 0}, {2, 1, 1}, {2, 0, 1},  // 8-11: B only
    };

    auto write_processor = [&](int proc, const std::vector<std::vector<double>>& pts,
                               const std::vector<std::int64_t>& conn) {
        meshioplusplus::Mesh m;
        m.AssignPoints(mt::points_from(pts));
        m.AddCellBlock("hexahedron", mt::conn_from({conn}));
        const fs::path marker =
            base / ("processor" + std::to_string(proc)) / "case.foam";
        meshioplusplus::write_openfoam(marker.string(), m, {});
        return base / ("processor" + std::to_string(proc)) / "constant" / "polyMesh";
    };

    const fs::path poly0 = write_processor(0, {g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7]},
                                           {0, 1, 2, 3, 4, 5, 6, 7});
    const fs::path poly1 = write_processor(1, {g[4], g[5], g[6], g[7], g[8], g[9], g[10], g[11]},
                                           {0, 1, 2, 3, 4, 5, 6, 7});

    const PolyMeshRead p0 = read_polymesh(poly0), p1 = read_polymesh(poly1);
    ASSERT_EQ(p0.mFaces.size(), 6u);
    ASSERT_EQ(p1.mFaces.size(), 6u);
    ASSERT_TRUE(p0.mNeighbour.empty());  // a lone cell has no internal face of its own
    ASSERT_TRUE(p1.mNeighbour.empty());

    // The shared face is the one whose every corner sits at x=1, in EACH
    // processor's own written point coordinates.
    auto find_interface_face = [](const PolyMeshRead& p) -> std::int64_t {
        for (std::size_t f = 0; f < p.mFaces.size(); ++f) {
            bool all_x1 = true;
            for (std::int64_t nid : p.mFaces[f])
                if (std::abs(p.mPoints[static_cast<std::size_t>(nid)][0] - 1.0) > 1e-9)
                    all_x1 = false;
            if (all_x1)
                return static_cast<std::int64_t>(f);
        }
        return -1;
    };
    const std::int64_t iface0 = find_interface_face(p0);
    const std::int64_t iface1 = find_interface_face(p1);
    ASSERT_GE(iface0, 0);
    ASSERT_GE(iface1, 0);

    constexpr std::int64_t kInterfaceGlobalFace = 0;
    auto write_face_addr = [&](const fs::path& poly, std::size_t n, std::int64_t iface,
                               bool flip, std::int64_t& next_id) {
        std::vector<std::int64_t> addr(n);
        for (std::size_t f = 0; f < n; ++f)
            addr[f] = static_cast<std::int64_t>(f) == iface
                          ? (flip ? -(kInterfaceGlobalFace + 1) : (kInterfaceGlobalFace + 1))
                          : (next_id++) + 1;
        std::ofstream f(poly / "faceProcAddressing");
        f << "FoamFile\n{\n format ascii;\n class labelList;\n object "
             "faceProcAddressing;\n}\n"
          << addr.size() << "\n(\n";
        for (std::int64_t v : addr)
            f << v << "\n";
        f << ")\n";
    };
    std::int64_t next_id = 1;  // 0 is reserved for the interface face
    write_face_addr(poly0, p0.mFaces.size(), iface0, /*flip=*/false, next_id);
    write_face_addr(poly1, p1.mFaces.size(), iface1, /*flip=*/true, next_id);

    auto write_label_list = [](const fs::path& file, const std::vector<std::int64_t>& vals) {
        std::ofstream f(file);
        f << "FoamFile\n{\n format ascii;\n class labelList;\n object "
          << file.filename().string() << ";\n}\n"
          << vals.size() << "\n(\n";
        for (std::int64_t v : vals)
            f << v << "\n";
        f << ")\n";
    };
    write_label_list(poly0 / "pointProcAddressing", {0, 1, 2, 3, 4, 5, 6, 7});
    write_label_list(poly1 / "pointProcAddressing", {4, 5, 6, 7, 8, 9, 10, 11});
    write_label_list(poly0 / "cellProcAddressing", {0});
    write_label_list(poly1 / "cellProcAddressing", {1});
    write_label_list(poly0 / "boundaryProcAddressing", {0});  // defaultFaces -> global patch 0
    write_label_list(poly1 / "boundaryProcAddressing", {0});

    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh back = meshioplusplus::read_openfoam(base.string(), info);

    EXPECT_EQ(back.NumPoints(), 12u);
    std::size_t nhex = 0, nquad = 0;
    for (const auto cb : back.CellRange()) {
        if (cb.Type() == "hexahedron")
            nhex += cb.NumCells();
        else if (cb.Type() == "quad")
            nquad += cb.NumCells();
    }
    EXPECT_EQ(nhex, 2u) << "the shared face must have become one internal face, not two "
                           "boundary ones, leaving both cells intact";
    EXPECT_EQ(nquad, 10u) << "6 + 6 boundary faces minus the 2 that became internal";

    std::error_code ec;
    fs::remove_all(base, ec);
}

// Roadmap §1 tier B2: time-directory fields (internalField only, uniform and
// nonuniform, both cell and point homes).
TEST(OpenFoamFields, ReadsUniformAndNonuniformInternalFieldsFromTimeZero) {
    meshioplusplus::Mesh m = hex_grid(2);  // 8 hex cells, 27 points
    const fs::path base = write_case(m);
    const fs::path zero = base / "0";
    fs::create_directories(zero);

    std::ofstream(zero / "p") << "FoamFile\n{\n format ascii;\n class volScalarField;\n object "
                                 "p;\n}\ninternalField   uniform 42;\n";
    {
        std::ofstream f(zero / "U");
        f << "FoamFile\n{\n format ascii;\n class volVectorField;\n object U;\n}\n";
        f << "internalField   nonuniform List<vector>\n8\n(\n";
        for (int i = 0; i < 8; ++i)
            f << "(" << i << " " << (i * 2) << " " << (i * 3) << ")\n";
        f << ")\n;\n";
    }
    {
        std::ofstream f(zero / "pointDisp");
        f << "FoamFile\n{\n format ascii;\n class pointScalarField;\n object pointDisp;\n}\n";
        f << "internalField   nonuniform List<scalar>\n27\n(\n";
        for (int i = 0; i < 27; ++i)
            f << i << "\n";
        f << ")\n;\n";
    }

    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::ReadOptions opts;
    const meshioplusplus::Mesh mesh =
        meshioplusplus::read_openfoam((base / "case.foam").string(), opts, info);

    ASSERT_TRUE(mesh.HasCellData("p"));
    ASSERT_TRUE(mesh.HasCellData("U"));
    ASSERT_TRUE(mesh.HasPointData("pointDisp"));

    std::size_t hexblock = static_cast<std::size_t>(-1);
    for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b)
        if (mesh.Cells(b).Type() == "hexahedron") {
            hexblock = b;
            break;
        }
    ASSERT_NE(hexblock, static_cast<std::size_t>(-1));

    const meshioplusplus::NDArray& p = mesh.CellData("p", hexblock);
    for (std::size_t i = 0; i < 8; ++i)
        EXPECT_DOUBLE_EQ(p.As<double>()[i], 42.0);

    const meshioplusplus::NDArray& u = mesh.CellData("U", hexblock);
    for (int i = 0; i < 8; ++i) {
        EXPECT_DOUBLE_EQ(u.As<double>()[static_cast<std::size_t>(i) * 3 + 0], i);
        EXPECT_DOUBLE_EQ(u.As<double>()[static_cast<std::size_t>(i) * 3 + 1], i * 2);
        EXPECT_DOUBLE_EQ(u.As<double>()[static_cast<std::size_t>(i) * 3 + 2], i * 3);
    }

    const meshioplusplus::NDArray& pd = mesh.PointData("pointDisp");
    ASSERT_EQ(pd.Shape()[0], 27u);
    for (std::size_t i = 0; i < 27; ++i)
        EXPECT_DOUBLE_EQ(pd.As<double>()[i], static_cast<double>(i));

    // A boundary (non-volume) cell block carries NaN for a volume-only field.
    for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b) {
        if (b == hexblock)
            continue;
        const meshioplusplus::NDArray& bp = mesh.CellData("p", b);
        for (std::size_t i = 0; i < bp.Shape()[0]; ++i)
            EXPECT_TRUE(std::isnan(bp.As<double>()[i]));
    }

    std::error_code ec;
    fs::remove_all(base, ec);
}

TEST(OpenFoamFields, TimeStepSelectsTheDirectoryAndMetadataListsTheValues) {
    meshioplusplus::Mesh m = hex_grid(1);  // 1 hex cell, 8 points
    const fs::path base = write_case(m);
    for (const auto& [dir, val] : std::vector<std::pair<std::string, double>>{{"0", 0.0}, {"1", 10.0}}) {
        const fs::path td = base / dir;
        fs::create_directories(td);
        std::ofstream(td / "p") << "FoamFile\n{\n format ascii;\n class volScalarField;\n object "
                                   "p;\n}\ninternalField   uniform "
                                << val << ";\n";
    }

    meshioplusplus::OpenFoamInfo info;
    meshioplusplus::ReadOptions opts;
    opts.mTimeStep = 1;
    const meshioplusplus::Mesh mesh =
        meshioplusplus::read_openfoam((base / "case.foam").string(), opts, info);
    ASSERT_TRUE(mesh.HasCellData("p"));
    std::size_t hexblock = 0;
    for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b)
        if (mesh.Cells(b).Type() == "hexahedron") {
            hexblock = b;
            break;
        }
    EXPECT_DOUBLE_EQ(mesh.CellData("p", hexblock).As<double>()[0], 10.0);

    const meshioplusplus::MeshMetadata meta =
        meshioplusplus::read_openfoam_metadata((base / "case.foam").string(), meshioplusplus::ReadOptions{});
    ASSERT_EQ(meta.mTimeValues.size(), 2u);
    EXPECT_DOUBLE_EQ(meta.mTimeValues[0], 0.0);
    EXPECT_DOUBLE_EQ(meta.mTimeValues[1], 1.0);  // directory name "1", not the field's own value

    std::error_code ec;
    fs::remove_all(base, ec);
}

// checkMesh is the ONLY oracle that catches a convention error -- a globally
// inverted winding passes every internally-consistent check above. It is
// virtually never installed, so this skips loudly rather than silently.
TEST(OpenFoamWrite, CheckMeshAcceptsOurOutput) {
    const char* exe = std::getenv("MESHIO_OPENFOAM_CHECKMESH");
    if (!exe || !*exe)
        GTEST_SKIP() << "set MESHIO_OPENFOAM_CHECKMESH=/path/to/checkMesh to run the "
                        "authoritative OpenFOAM oracle";
    const fs::path base = write_case(hex_grid(3));
    const std::string cmd =
        std::string(exe) + " -case " + base.string() + " > " + (base / "check.log").string() +
        " 2>&1";
    const int rc = std::system(cmd.c_str());
    std::ifstream log(base / "check.log");
    const std::string text((std::istreambuf_iterator<char>(log)),
                           std::istreambuf_iterator<char>());
    EXPECT_EQ(rc, 0) << text;
    EXPECT_NE(text.find("Mesh OK"), std::string::npos) << text;
    std::error_code ec;
    fs::remove_all(base, ec);
}

// The reader fix that had to come with `type` parsing: a patch carrying a
// nested sub-dictionary was truncated at the FIRST '}'.
TEST(OpenFoam, ParseBoundaryHandlesNestedBracesAndReadsType) {
    const fs::path base = temp_case_dir();
    const fs::path poly = base / "constant" / "polyMesh";
    fs::create_directories(poly);
    auto hdr = [](const std::string& cls, const std::string& obj) {
        return "FoamFile\n{\n format ascii;\n class " + cls + ";\n object " + obj + ";\n}\n";
    };
    std::ofstream(poly / "points")
        << hdr("vectorField", "points")
        << "8\n(\n(0 0 0)\n(1 0 0)\n(1 1 0)\n(0 1 0)\n(0 0 1)\n(1 0 1)\n(1 1 1)\n(0 1 1)\n)\n";
    std::ofstream(poly / "faces")
        << hdr("faceList", "faces")
        << "6\n(\n4(0 3 2 1)\n4(4 5 6 7)\n4(0 1 5 4)\n4(2 3 7 6)\n4(1 2 6 5)\n4(0 4 7 3)\n)\n";
    std::ofstream(poly / "owner") << hdr("labelList", "owner") << "6\n(\n0\n0\n0\n0\n0\n0\n)\n";
    // The first patch nests `transform { ... }`, exactly as a real cyclicAMI
    // does. Taking the first '}' truncates it and invents garbage patches.
    std::ofstream(poly / "boundary")
        << hdr("polyBoundaryMesh", "boundary")
        << "2\n(\n"
           "inlet\n{\n    type            cyclicAMI;\n    transform\n    {\n        type    "
           "translational;\n    }\n    nFaces          2;\n    startFace       0;\n}\n"
           "walls\n{\n    type            wall;\n    nFaces          4;\n    startFace       "
           "2;\n}\n)\n";

    meshioplusplus::OpenFoamInfo info;
    meshioplusplus::Mesh mesh = meshioplusplus::read_openfoam(base.string(), info);
    EXPECT_EQ(info.mCellTags.size(), 2u) << "the nested sub-dictionary was mis-parsed";
    EXPECT_EQ(info.mCellTags.at(-1).front(), "inlet");
    EXPECT_EQ(info.mCellTags.at(-2).front(), "walls");
    // `type` must be read -- the Python reader has always done this.
    EXPECT_EQ(info.mPatchTypes.at(-1), "cyclicAMI");
    EXPECT_EQ(info.mPatchTypes.at(-2), "wall");
    std::error_code ec;
    fs::remove_all(base, ec);
}

// Roadmap 3.3.3: a decomposed case written by `write_openfoam_decomposed` reassembles to the mesh
// the serial writer produces (a 2x2x2 grid in two processors, so each owns internal faces of its
// own as well as the ones it shares), a repeated step is refused, and a region lands under
// `constant/<region>/polyMesh`.
namespace {

meshioplusplus::Mesh foam_test_grid(std::size_t N) {
    const std::size_t p = N + 1;
    meshioplusplus::NDArray pts(meshioplusplus::DType::Float64, {p * p * p, std::size_t{3}});
    for (std::size_t k = 0; k < p; ++k)
        for (std::size_t j = 0; j < p; ++j)
            for (std::size_t i = 0; i < p; ++i) {
                double* x = pts.As<double>() + 3 * (i + p * (j + p * k));
                x[0] = static_cast<double>(i);
                x[1] = static_cast<double>(j);
                x[2] = static_cast<double>(k);
            }
    meshioplusplus::NDArray conn(meshioplusplus::DType::Int64,
                                 {N * N * N, std::size_t{8}});
    std::size_t c = 0;
    const auto id = [&](std::size_t i, std::size_t j, std::size_t k) {
        return static_cast<std::int64_t>(i + p * (j + p * k));
    };
    for (std::size_t k = 0; k < N; ++k)
        for (std::size_t j = 0; j < N; ++j)
            for (std::size_t i = 0; i < N; ++i, ++c) {
                std::int64_t* row = conn.As<std::int64_t>() + 8 * c;
                const std::int64_t v[8] = {id(i, j, k),         id(i + 1, j, k),
                                           id(i + 1, j + 1, k), id(i, j + 1, k),
                                           id(i, j, k + 1),     id(i + 1, j, k + 1),
                                           id(i + 1, j + 1, k + 1), id(i, j + 1, k + 1)};
                std::copy(v, v + 8, row);
            }
    meshioplusplus::Mesh mesh;
    mesh.AssignPoints(std::move(pts));
    mesh.AddCellBlock("hexahedron", std::move(conn));
    return mesh;
}

}  // namespace

TEST(OpenFoamWriter, DecomposedCaseReassemblesAndRegionsAndSeriesWrite) {
    meshioplusplus::Mesh mesh = foam_test_grid(2);
    meshioplusplus::NDArray part(meshioplusplus::DType::Int64, {std::size_t{8}});
    for (std::size_t c = 0; c < 8; ++c)
        part.As<std::int64_t>()[c] = c < 4 ? 0 : 1;
    mesh.AddCellData("partition:part", {std::move(part)});

    const fs::path base = temp_case_dir();
    meshioplusplus::OpenFoamInfo info;
    EXPECT_THROW(meshioplusplus::write_openfoam_decomposed((base / "x.foam").string(),
                                                           foam_test_grid(2), info, {}),
                 meshioplusplus::WriteError);
    meshioplusplus::write_openfoam_decomposed((base / "case.foam").string(), mesh, info, {});
    ASSERT_TRUE(fs::exists(base / "processor1" / "constant" / "polyMesh" / "faceProcAddressing"));
    meshioplusplus::OpenFoamInfo read_info;
    const meshioplusplus::Mesh back =
        meshioplusplus::read_openfoam((base / "case.foam").string(), read_info);
    ASSERT_GE(back.NumCellBlocks(), 1u);
    EXPECT_EQ(back.Cells(0).Type(), "hexahedron");
    EXPECT_EQ(back.Cells(0).NumCells(), 8u);
    EXPECT_EQ(back.NumPoints(), 27u);

    meshioplusplus::write_openfoam_region((base / "mr.foam").string(), foam_test_grid(1), info, {},
                                          "fluid");
    EXPECT_TRUE(fs::exists(base / "constant" / "fluid" / "polyMesh" / "points"));
    EXPECT_THROW(meshioplusplus::write_openfoam_region((base / "mr.foam").string(),
                                                       foam_test_grid(1), info, {}, "a/b"),
                 meshioplusplus::WriteError);

    meshioplusplus::OpenFoamSeriesWriter series((base / "series" / "case.foam").string());
    series.Write(0.25, foam_test_grid(1));
    series.Write(0.5, foam_test_grid(1));
    series.Finalize();
    EXPECT_EQ(series.NumSteps(), 2u);
    EXPECT_THROW(series.Write(1.0, foam_test_grid(2)), meshioplusplus::WriteError);
    std::error_code ec;
    fs::remove_all(base, ec);
}

// Roadmap 3.3.3: zones in a decomposed case. Each processor writes every zone in its own ids; the
// reader merges them back through the `*ProcAddressing` lists, so the reconstruction names the same
// cells, points and faces as the serial case does.
namespace {

// A zone's entries as sorted rows (pairs for a side zone), so two meshes compare by content.
std::map<std::string, std::vector<std::int64_t>> zone_entries(const meshioplusplus::Mesh& rMesh) {
    std::map<std::string, std::vector<std::int64_t>> out;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        const std::size_t width = r.mKind == meshioplusplus::RegionKind::Side ? 2 : 1;
        std::vector<std::int64_t> rows;
        for (std::size_t k = 0; k < r.NumEntries(); ++k) {
            std::int64_t key = 0;
            for (std::size_t w = 0; w < width; ++w)
                key = key * 1000 + r.Entries()[k * width + w];
            rows.push_back(key);
        }
        std::sort(rows.begin(), rows.end());
        out[r.mName + ":" + std::to_string(static_cast<int>(r.mKind))] = std::move(rows);
    }
    return out;
}

}  // namespace

TEST(OpenFoamWriter, DecomposedZonesComeBackAsInTheSerialCase) {
    for (const bool binary : {false, true}) {
        meshioplusplus::Mesh mesh = foam_test_grid(3);
        meshioplusplus::NDArray part(meshioplusplus::DType::Int64, {std::size_t{27}});
        for (std::size_t c = 0; c < 27; ++c)
            part.As<std::int64_t>()[c] = static_cast<std::int64_t>(c / 9);
        mesh.AddCellData("partition:part", {std::move(part)});
        mesh.AddRegion(
            meshioplusplus::Region("core", meshioplusplus::RegionKind::Cell, i64({2, 11, 20})));
        mesh.AddRegion(meshioplusplus::Region("none", meshioplusplus::RegionKind::Cell, i64({})));
        mesh.AddRegion(
            meshioplusplus::Region("rim", meshioplusplus::RegionKind::Point, i64({0, 5, 17, 40})));
        // facets 0 and 5 of cells on either side of both processor boundaries
        mesh.AddRegion(
            meshioplusplus::Region("baffle", meshioplusplus::RegionKind::Side,
                                   i64_pairs({0, 0, 0, 5, 8, 5, 9, 0, 17, 5, 18, 0, 26, 5})));

        const fs::path base = temp_case_dir();
        meshioplusplus::OpenFoamInfo info;
        meshioplusplus::OpenFoamWriteOptions options;
        options.mBinary = binary;
        meshioplusplus::write_openfoam_decomposed((base / "dec" / "case.foam").string(), mesh, info,
                                                  options);
        meshioplusplus::write_openfoam((base / "ser" / "case.foam").string(), mesh, info, options);
        for (int p = 0; p < 3; ++p) {
            const fs::path poly =
                base / "dec" / ("processor" + std::to_string(p)) / "constant" / "polyMesh";
            EXPECT_TRUE(fs::exists(poly / "cellZones"));
            EXPECT_TRUE(fs::exists(poly / "faceZones"));
            EXPECT_TRUE(fs::exists(poly / "pointZones"));
        }
        meshioplusplus::OpenFoamInfo dec_info, ser_info;
        const meshioplusplus::Mesh dec =
            meshioplusplus::read_openfoam((base / "dec" / "case.foam").string(), dec_info);
        const meshioplusplus::Mesh ser =
            meshioplusplus::read_openfoam((base / "ser" / "case.foam").string(), ser_info);
        const auto got = zone_entries(dec);
        const auto want = zone_entries(ser);
        EXPECT_EQ(got.size(), 4u);
        EXPECT_EQ(got, want);
        EXPECT_EQ(
            got.at("core:" + std::to_string(static_cast<int>(meshioplusplus::RegionKind::Cell))),
            (std::vector<std::int64_t>{2, 11, 20}));
        EXPECT_EQ(
            got.at("baffle:" + std::to_string(static_cast<int>(meshioplusplus::RegionKind::Side)))
                .size(),
            7u);

        // Dropping the zones drops the files with them.
        meshioplusplus::Mesh plain = foam_test_grid(3);
        meshioplusplus::NDArray part2(meshioplusplus::DType::Int64, {std::size_t{27}});
        for (std::size_t c = 0; c < 27; ++c)
            part2.As<std::int64_t>()[c] = static_cast<std::int64_t>(c / 9);
        plain.AddCellData("partition:part", {std::move(part2)});
        meshioplusplus::write_openfoam_decomposed((base / "dec" / "case.foam").string(), plain,
                                                  info, options);
        EXPECT_FALSE(
            fs::exists(base / "dec" / "processor0" / "constant" / "polyMesh" / "cellZones"));
        std::error_code ec;
        fs::remove_all(base, ec);
    }
}

// Roadmap 3.3.3: boundary-patch values. The values of the boundary-face cell blocks are written
// as per-patch `boundaryField` entries and read back onto those blocks.
namespace {

// The rows of every 2D block of @p rMesh's cell-data array @p rName, keyed by the sorted corner
// coordinates of the face (so the comparison does not depend on the row or corner order).
std::map<std::string, std::vector<double>> boundary_face_values(const meshioplusplus::Mesh& rMesh,
                                                                const std::string& rName,
                                                                std::size_t Cols) {
    std::map<std::string, std::vector<double>> out;
    const double* pts = rMesh.Points().As<double>();
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const std::string type(rMesh.Cells(b).Type());
        if (type != "triangle" && type != "quad")
            continue;
        const meshioplusplus::NDArray& conn = rMesh.Cells(b).Conn();
        const std::size_t npc = rMesh.Cells(b).NodesPerCell();
        const meshioplusplus::NDArray& vals = rMesh.CellData(rName, b);
        for (std::size_t r = 0; r < rMesh.Cells(b).NumCells(); ++r) {
            std::vector<std::string> corners;
            for (std::size_t k = 0; k < npc; ++k) {
                const double* x =
                    pts + 3 * static_cast<std::size_t>(conn.As<std::int64_t>()[r * npc + k]);
                corners.push_back(std::to_string(x[0]) + "," + std::to_string(x[1]) + "," +
                                  std::to_string(x[2]));
            }
            std::sort(corners.begin(), corners.end());
            std::string key;
            for (const std::string& c : corners)
                key += c + ";";
            out[key].assign(vals.As<double>() + r * Cols, vals.As<double>() + (r + 1) * Cols);
        }
    }
    return out;
}

// A unit hex with its six faces as quads: tags -1 (y-min, x-min, z-min) and -2 (the rest), and
// scalar p / vector U values per quad that tell the faces apart (NaN without @p WithValues).
meshioplusplus::Mesh single_hex_with_boundary(bool WithValues, meshioplusplus::OpenFoamInfo& rInfo,
                                              double Offset = 0.0) {
    meshioplusplus::Mesh m = hex_grid(1);
    m.AddCellBlock("quad", mt::conn_from({{0, 4, 5, 1},  // y = 0
                                          {0, 2, 6, 4},  // z = 0
                                          {0, 1, 3, 2},
                                          {1, 5, 7, 3},
                                          {2, 3, 7, 6},
                                          {4, 6, 7, 5}}));
    m.AddCellData("cell_tags",
                  {mt::int_data_array({0}), mt::int_data_array({-1, -1, -1, -2, -2, -2})});
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> p_quad(6, nan), u_quad(18, nan);
    if (WithValues)
        for (std::size_t i = 0; i < 6; ++i) {
            p_quad[i] = 10.0 + Offset + static_cast<double>(i);
            u_quad[3 * i] = static_cast<double>(i);
            u_quad[3 * i + 1] = 2.0 * static_cast<double>(i);
            u_quad[3 * i + 2] = 3.0 * static_cast<double>(i) + 0.5;
        }
    m.AddCellData("p", {mt::data_array({1.5}), mt::data_array(p_quad)});
    m.AddCellData("U", {mt::data_array({1.0, 2.0, 3.0}, 3), mt::data_array(u_quad, 3)});
    rInfo.mCellTags[-1] = {"alpha"};
    rInfo.mCellTags[-2] = {"beta"};
    rInfo.mPatchTypes[-1] = "wall";
    return m;
}

std::string slurp(const fs::path& rPath) {
    std::ifstream f(rPath, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

meshioplusplus::Mesh read_case(const fs::path& rCase, meshioplusplus::OpenFoamInfo& rInfo,
                               const meshioplusplus::ReadOptions& rOptions = {}) {
    return meshioplusplus::read_openfoam(rCase.string(), rOptions, rInfo);
}

}  // namespace

class OpenFoamPatchValues : public ::testing::TestWithParam<bool> {};

TEST_P(OpenFoamPatchValues, WrittenPerPatchAndReadBackOntoTheBoundaryBlocks) {
    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh m = single_hex_with_boundary(true, info);
    meshioplusplus::OpenFoamWriteOptions opts;
    opts.mBinary = GetParam();
    const fs::path base = temp_case_dir();
    meshioplusplus::write_openfoam((base / "case.foam").string(), m, info, opts);

    // Both patches get explicit entries ahead of the catch-all.
    const std::string text = slurp(base / "0" / "p");
    EXPECT_NE(text.find("    alpha\n    {\n        type            calculated;\n"
                        "        value           nonuniform List<scalar>\n3\n("),
              std::string::npos);
    EXPECT_NE(text.find("    beta\n    {"), std::string::npos);
    EXPECT_LT(text.find("    beta\n"), text.find("\".*\""));

    meshioplusplus::OpenFoamInfo read_info;
    const meshioplusplus::Mesh back = read_case(base / "case.foam", read_info);
    EXPECT_EQ(boundary_face_values(back, "p", 1), boundary_face_values(m, "p", 1));
    EXPECT_EQ(boundary_face_values(back, "U", 3), boundary_face_values(m, "U", 3));
    ASSERT_EQ(boundary_face_values(back, "p", 1).size(), 6u);
    EXPECT_DOUBLE_EQ(back.CellData("p", 0).As<double>()[0], 1.5);
    std::error_code ec;
    fs::remove_all(base, ec);
}

INSTANTIATE_TEST_SUITE_P(Formats, OpenFoamPatchValues, ::testing::Values(false, true));

TEST(OpenFoamPatchValuesOutput, NoBoundaryValuesLeavesTheCatchAllAlone) {
    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh m = single_hex_with_boundary(false, info);  // NaN on the quads
    const fs::path base = temp_case_dir();
    meshioplusplus::write_openfoam((base / "case.foam").string(), m, info);
    const std::string text = slurp(base / "0" / "p");
    EXPECT_EQ(text.find("alpha"), std::string::npos);
    EXPECT_EQ(text.find("nonuniform", text.find("boundaryField")), std::string::npos);
    EXPECT_NE(text.find("boundaryField\n{\n    \".*\"\n"), std::string::npos);
    std::error_code ec;
    fs::remove_all(base, ec);
}

TEST(OpenFoamPatchValuesOutput, ConstraintPatchesGetNoValues) {
    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh m = single_hex_with_boundary(true, info);
    info.mPatchTypes[-2] = "symmetry";
    const fs::path base = temp_case_dir();
    meshioplusplus::write_openfoam((base / "case.foam").string(), m, info);
    const std::string text = slurp(base / "0" / "p");
    EXPECT_NE(text.find("    alpha\n"), std::string::npos);
    EXPECT_EQ(text.find("    beta\n"), std::string::npos);
    std::error_code ec;
    fs::remove_all(base, ec);
}

TEST(OpenFoamPatchValuesOutput, SeriesStepsKeepTheirOwnPatchValues) {
    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh step0 = single_hex_with_boundary(true, info);
    const meshioplusplus::Mesh step1 = single_hex_with_boundary(true, info, 100.0);
    const fs::path base = temp_case_dir();
    meshioplusplus::OpenFoamSeriesWriter series((base / "case.foam").string(), info);
    series.Write(0.0, step0);
    series.Write(1.0, step1);
    series.Finalize();
    meshioplusplus::ReadOptions opts;
    opts.mTimeStep = 1;
    meshioplusplus::OpenFoamInfo read_info;
    const meshioplusplus::Mesh back = read_case(base / "case.foam", read_info, opts);
    EXPECT_EQ(boundary_face_values(back, "p", 1), boundary_face_values(step1, "p", 1));
    EXPECT_NE(boundary_face_values(back, "p", 1), boundary_face_values(step0, "p", 1));
    std::error_code ec;
    fs::remove_all(base, ec);
}

TEST(OpenFoamPatchValuesOutput, RegionAndDecomposedCasesCarryPatchValues) {
    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh m = single_hex_with_boundary(true, info);
    const fs::path base = temp_case_dir();
    meshioplusplus::write_openfoam_region((base / "case.foam").string(), m, info, {}, "fluid");
    EXPECT_NE(slurp(base / "0" / "fluid" / "p").find("    alpha\n"), std::string::npos);
    meshioplusplus::OpenFoamInfo read_info;
    read_info.mRegion = "fluid";
    const meshioplusplus::Mesh back = read_case(base / "case.foam", read_info);
    EXPECT_EQ(boundary_face_values(back, "p", 1), boundary_face_values(m, "p", 1));

    // 2x2x2 grid in two processors: the faces between them become `processor` patches with no
    // values, every outer face keeps its own.
    meshioplusplus::OpenFoamInfo dinfo;
    meshioplusplus::Mesh grid = hex_grid(2);
    {
        const fs::path tmp = temp_case_dir();
        meshioplusplus::write_openfoam((tmp / "case.foam").string(), grid, dinfo);
        grid = read_case(tmp / "case.foam", dinfo);
        std::error_code rc;
        fs::remove_all(tmp, rc);
    }
    ASSERT_EQ(grid.NumCellBlocks(), 2u);  // 8 hexahedra, then the 24 boundary quads
    std::vector<double> pv(8), pq(24);
    for (std::size_t i = 0; i < 8; ++i)
        pv[i] = static_cast<double>(i);
    for (std::size_t i = 0; i < 24; ++i)
        pq[i] = 100.0 + static_cast<double>(i);
    grid.AddCellData("p", {mt::data_array(pv), mt::data_array(pq)});
    meshioplusplus::NDArray part(meshioplusplus::DType::Int64, {std::size_t{8}});
    for (std::size_t c = 0; c < 8; ++c)
        part.As<std::int64_t>()[c] = c < 4 ? 0 : 1;
    grid.AddCellData("partition:part",
                     {std::move(part), mt::int_data_array(std::vector<std::int32_t>(24, 0))});
    const fs::path dbase = temp_case_dir();
    meshioplusplus::write_openfoam_decomposed((dbase / "case.foam").string(), grid, dinfo, {});
    const std::string proc0 = slurp(dbase / "processor0" / "0" / "p");
    EXPECT_NE(proc0.find("defaultFaces"), std::string::npos);
    EXPECT_EQ(proc0.find("procBoundary"), std::string::npos);
    meshioplusplus::OpenFoamInfo dread;
    const meshioplusplus::Mesh dback = read_case(dbase / "case.foam", dread);
    EXPECT_EQ(boundary_face_values(dback, "p", 1), boundary_face_values(grid, "p", 1));
    ASSERT_EQ(boundary_face_values(dback, "p", 1).size(), 24u);
    std::error_code ec;
    fs::remove_all(base, ec);
    fs::remove_all(dbase, ec);
}

TEST(OpenFoamFields, ReadsBoundaryFieldValuesInEveryFormOpenFoamWritesThem) {
    // make_hex_case: bottom 1 face, top 1 face, sides 4 faces.
    const fs::path hand = make_hex_case();
    const fs::path zero = hand / "0";
    fs::create_directories(zero);
    std::ofstream(zero / "p")
        << "FoamFile\n{\n format ascii;\n class volScalarField;\n object p;\n}\n"
           "dimensions [0 0 0 0 0 0 0];\ninternalField uniform 42;\n"
           "boundaryField\n{\n"
           "    // a comment { with ( stray ) braces\n"
           "    bottom { type fixedValue; value uniform 7; }\n"
           "    top\n    {\n        type fixedValue;\n"
           "        value nonuniform List<scalar> 1(5);\n    }\n"
           "    sides\n    {\n        type mixed;\n"
           "        refGradient nonuniform List<scalar> 4(9 9 9 9);\n"
           "        valueFraction { type table; table ((0 1) (1 2)); }\n"
           "        value nonuniform List<scalar>\n4\n(\n1\n2\n3\n4\n)\n;\n    }\n"
           "    \".*\" { type calculated; value uniform 0; }\n}\n";
    std::ofstream(zero / "U")
        << "FoamFile\n{\n format ascii;\n class volVectorField;\n object U;\n}\n"
           "internalField uniform (1 0 0);\n"
           "boundaryField\n{\n"
           "    bottom { type zeroGradient; }\n"
           "    top { type fixedValue; value $internalField; }\n"
           "    sides { type fixedValue; value nonuniform List<vector> 4{(0 0 1)}; }\n}\n";
    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh mesh = read_case(hand / "case.foam", info);
    std::size_t quad = static_cast<std::size_t>(-1);
    for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b)
        if (mesh.Cells(b).Type() == "quad")
            quad = b;
    ASSERT_NE(quad, static_cast<std::size_t>(-1));
    const double* p = mesh.CellData("p", quad).As<double>();
    const double expect_p[6] = {7, 5, 1, 2, 3, 4};
    for (int i = 0; i < 6; ++i)
        EXPECT_DOUBLE_EQ(p[i], expect_p[i]) << "row " << i;
    const double* u = mesh.CellData("U", quad).As<double>();
    for (int i = 0; i < 2; ++i)
        for (int c = 0; c < 3; ++c)
            EXPECT_TRUE(std::isnan(u[3 * i + c])) << "row " << i;
    for (int i = 2; i < 6; ++i) {
        EXPECT_DOUBLE_EQ(u[3 * i], 0.0);
        EXPECT_DOUBLE_EQ(u[3 * i + 1], 0.0);
        EXPECT_DOUBLE_EQ(u[3 * i + 2], 1.0);
    }
    std::error_code ec;
    fs::remove_all(hand, ec);
}

// A field file whose `arch` names a width the decoders do not implement is refused, not divided by.
TEST(OpenFoamFields, AnUnsupportedArchWidthInAFieldFileIsRefused) {
    meshioplusplus::OpenFoamInfo info;
    const meshioplusplus::Mesh m = single_hex_with_boundary(true, info);
    meshioplusplus::OpenFoamWriteOptions opts;
    opts.mBinary = true;
    const fs::path base = temp_case_dir();
    meshioplusplus::write_openfoam((base / "case.foam").string(), m, info, opts);
    std::string text = slurp(base / "0" / "p");
    const std::size_t at = text.find("scalar=64");
    ASSERT_NE(at, std::string::npos);
    text.replace(at, 9, "scalar=6G");
    std::ofstream(base / "0" / "p", std::ios::binary) << text;
    meshioplusplus::OpenFoamInfo read_info;
    EXPECT_THROW(read_case(base / "case.foam", read_info), meshioplusplus::ReadError);
    std::error_code ec;
    fs::remove_all(base, ec);
}

// Roadmap 3.3.3: `constant/regionProperties`, written only on request and read back in order.
TEST(OpenFoamRegionProperties, WrittenOnRequestAndReadBackInOrder) {
    const fs::path base = temp_case_dir();
    meshioplusplus::OpenFoamInfo info;
    meshioplusplus::write_openfoam_region((base / "case.foam").string(), foam_test_grid(1), info,
                                          {}, "water");
    EXPECT_FALSE(fs::exists(base / "constant" / "regionProperties"));

    const std::vector<meshioplusplus::OpenFoamRegionGroup> groups = {{"fluid", {"water", "air"}},
                                                                     {"solid", {"heater"}}};
    meshioplusplus::write_openfoam_region_properties((base / "case.foam").string(), groups);
    std::ifstream in(base / "constant" / "regionProperties");
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(text.find("regions\n(\n    fluid (water air)\n    solid (heater)\n);"),
              std::string::npos);
    EXPECT_NE(text.find("object      regionProperties;"), std::string::npos);

    for (const std::string& path : {(base / "case.foam").string(), base.string()}) {
        const auto back = meshioplusplus::read_openfoam_region_properties(path);
        ASSERT_EQ(back.size(), 2u);
        EXPECT_EQ(back[0].mFamily, "fluid");
        EXPECT_EQ(back[0].mRegions, (std::vector<std::string>{"water", "air"}));
        EXPECT_EQ(back[1].mFamily, "solid");
        EXPECT_EQ(back[1].mRegions, (std::vector<std::string>{"heater"}));
    }
    std::error_code ec;
    fs::remove_all(base, ec);
}

TEST(OpenFoamRegionProperties, RefusesWhatASolverCouldNotReadAndMalformedFiles) {
    const fs::path base = temp_case_dir();
    const std::string path = (base / "case.foam").string();
    using Groups = std::vector<meshioplusplus::OpenFoamRegionGroup>;
    EXPECT_THROW(meshioplusplus::write_openfoam_region_properties(path, Groups{}),
                 meshioplusplus::WriteError);
    for (const Groups& bad :
         {Groups{{"fluid", {}}}, Groups{{"fluid", {"a b"}}}, Groups{{"fluid", {".hidden"}}},
          Groups{{"fluid", {"a/b"}}}, Groups{{"", {"a"}}},
          Groups{{"fluid", {"a"}}, {"fluid", {"b"}}}, Groups{{"fluid", {"a"}}, {"solid", {"a"}}}})
        EXPECT_THROW(meshioplusplus::write_openfoam_region_properties(path, bad),
                     meshioplusplus::WriteError);
    EXPECT_FALSE(fs::exists(base / "constant" / "regionProperties"));
    EXPECT_THROW(meshioplusplus::write_openfoam_region_properties(
                     (base / "constant" / "polyMesh").string(), Groups{{"fluid", {"a"}}}),
                 meshioplusplus::WriteError);

    EXPECT_THROW(meshioplusplus::read_openfoam_region_properties(path), meshioplusplus::ReadError);
    fs::create_directories(base / "constant");
    for (const char* body : {"regions ( fluid water );", "regions ( fluid (water );", "x 1;"}) {
        std::ofstream(base / "constant" / "regionProperties") << body;
        EXPECT_THROW(meshioplusplus::read_openfoam_region_properties(path),
                     meshioplusplus::ReadError)
            << body;
    }
    // the multi-region error names the regions the file lists when no directory exists yet
    std::ofstream(base / "constant" / "regionProperties") << "regions ( fluid (a b) );";
    meshioplusplus::OpenFoamInfo info;
    try {
        meshioplusplus::read_openfoam(path, info);
        FAIL() << "expected a ReadError";
    } catch (const meshioplusplus::ReadError& rException) {
        EXPECT_NE(std::string(rException.what()).find("a, b"), std::string::npos);
    }
    std::error_code ec;
    fs::remove_all(base, ec);
}
