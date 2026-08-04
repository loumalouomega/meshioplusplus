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
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "mesh_fixtures.hpp"
#include "meshioplusplus/formats/openfoam.hpp"

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
