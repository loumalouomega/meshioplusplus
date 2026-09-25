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
#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/parse_guard.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "face_cells_common.hpp"

// ANSYS Fluent `.msh` (also TGrid and GAMBIT meshes). A face row lists its
// nodes and the two cells it separates, `n0 .. nk c0 c1` (hexadecimal in
// ASCII), and a cell zone normally declares only a range of cell ids, so volume
// cells are rebuilt from their faces: the right-hand normal of `n0 .. nk`
// points into c0 (in 2-D, walking n0 -> n1 leaves c0 on the left), so a face is
// outward for c1 and reversed for c0. Twin of ansys/_ansys.py.

namespace meshioplusplus {

namespace {

// Cell element types with a fixed node count.
bool fluent_cell_type(std::int64_t ElementType, std::string& rType, int& rNodes) {
    switch (ElementType) {
        case 1:
            rType = "triangle";
            rNodes = 3;
            return true;
        case 2:
            rType = "tetra";
            rNodes = 4;
            return true;
        case 3:
            rType = "quad";
            rNodes = 4;
            return true;
        case 4:
            rType = "hexahedron";
            rNodes = 8;
            return true;
        case 5:
            rType = "pyramid";
            rNodes = 5;
            return true;
        case 6:
            rType = "wedge";
            rNodes = 6;
            return true;
        default:
            return false;
    }
}

// Face element types with a fixed node count; 0 (mixed) and 5 (polygonal)
// rows lead with their node count.
int fluent_face_nodes(std::int64_t FaceType) {
    return (FaceType >= 2 && FaceType <= 4) ? static_cast<int>(FaceType) : 0;
}

constexpr std::int64_t kFluentInterior = 2;

bool fluent_is_space(char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

struct FluentReader {
    const std::string& mD;
    std::size_t mP = 0;

    explicit FluentReader(const std::string& rD) : mD(rD) {}

    bool Eof() const { return mP >= mD.size(); }
    char At() const { return mD[mP]; }

    void SkipWs() {
        while (mP < mD.size() && fluent_is_space(mD[mP]))
            ++mP;
    }

    // Past `Depth` unmatched closing brackets; quoted text is skipped whole.
    void SkipBalanced(int Depth) {
        while (Depth > 0 && mP < mD.size()) {
            const char c = mD[mP++];
            if (c == '(') {
                ++Depth;
            } else if (c == ')') {
                --Depth;
            } else if (c == '"') {
                const std::size_t q = mD.find('"', mP);
                mP = q == std::string::npos ? mD.size() : q + 1;
            }
        }
    }

    // Past `End of Binary Section NNNN)` if it follows, else past the section's
    // own closing bracket.
    void SkipBinaryTrailer() {
        SkipWs();
        if (mD.compare(mP, 21, "End of Binary Section") == 0) {
            const std::size_t q = mD.find(')', mP);
            mP = q == std::string::npos ? mD.size() : q + 1;
        } else if (mP < mD.size() && mD[mP] == ')') {
            ++mP;
        }
    }

    // The `(a b c ...)` group after a section index, as hex integers; returns
    // whether a body follows.
    bool Header(std::vector<std::int64_t>& rValues) {
        SkipWs();
        if (Eof() || At() != '(')
            throw ReadError("Fluent: expected a section header");
        const std::size_t q = mD.find(')', mP);
        if (q == std::string::npos)
            throw ReadError("Fluent: unterminated section header");
        rValues.clear();
        for (std::size_t i = mP + 1; i < q;) {
            while (i < q && fluent_is_space(mD[i]))
                ++i;
            const std::size_t j = i;
            while (i < q && !fluent_is_space(mD[i]))
                ++i;
            if (i > j)
                rValues.push_back(std::strtoll(mD.substr(j, i - j).c_str(), nullptr, 16));
        }
        mP = q + 1;
        // A body opens with '(' after nothing but blanks; anything else (a
        // closing bracket, `End of Binary Section`) makes a declaration.
        const std::size_t start = mP;
        SkipWs();
        if (!Eof() && At() == '(') {
            const bool newline_before = mD.find('\n', start) < mP;
            ++mP;
            // meshio's writer puts the body's '(' last on the header line and
            // starts the data on the next line.
            if (!newline_before) {
                if (mD.compare(mP, 2, "\r\n") == 0)
                    mP += 2;
                else if (mP < mD.size() && mD[mP] == '\n')
                    ++mP;
            }
            return true;
        }
        mP = start;
        SkipBinaryTrailer();
        return false;
    }

    std::vector<std::string_view> AsciiBody() {
        const std::size_t q = mD.find(')', mP);
        if (q == std::string::npos)
            throw ReadError("Fluent: unterminated section body");
        std::vector<std::string_view> tokens;
        const std::string_view body(mD.data() + mP, q - mP);
        for (std::size_t i = 0; i < body.size();) {
            while (i < body.size() && fluent_is_space(body[i]))
                ++i;
            const std::size_t j = i;
            while (i < body.size() && !fluent_is_space(body[i]))
                ++i;
            if (i > j)
                tokens.push_back(body.substr(j, i - j));
        }
        mP = q + 1;
        SkipWs();
        if (!Eof() && At() == ')')
            ++mP;
        return tokens;
    }

    // `Count` little-endian values of `Size` bytes (floats when `Float`).
    template <class T>
    std::vector<T> Binary(std::size_t Count) {
        const std::size_t size = Count * sizeof(T);
        if (mP > mD.size() || Count > (mD.size() - mP) / sizeof(T))  // no wrap
            throw ReadError("Fluent: binary section runs past the end of the file");
        std::vector<T> out(Count);
        if (size)
            std::memcpy(out.data(), mD.data() + mP, size);
        mP += size;
        SkipWs();
        if (!Eof() && At() == ')')
            ++mP;
        SkipBinaryTrailer();
        return out;
    }

    // A mixed binary face body: each row leads with its node count.
    template <class T>
    std::vector<std::int64_t> BinaryMixedFaces(std::size_t Count) {
        std::vector<std::int64_t> values;
        for (std::size_t r = 0; r < Count; ++r) {
            if (mP + sizeof(T) > mD.size())
                throw ReadError("Fluent: binary section runs past the end of the file");
            T n;
            std::memcpy(&n, mD.data() + mP, sizeof(T));
            const std::size_t row = static_cast<std::size_t>(n) + 3;
            // Compared by division: row * sizeof(T) wraps for a corrupt n.
            if (n < 0 || row > (mD.size() - mP) / sizeof(T))
                throw ReadError("Fluent: binary section runs past the end of the file");
            for (std::size_t k = 0; k < row; ++k) {
                T v;
                std::memcpy(&v, mD.data() + mP + k * sizeof(T), sizeof(T));
                values.push_back(static_cast<std::int64_t>(v));
            }
            mP += row * sizeof(T);
        }
        SkipWs();
        if (!Eof() && At() == ')')
            ++mP;
        SkipBinaryTrailer();
        return values;
    }
};

std::int64_t fluent_hex(std::string_view Tok) {
    return std::strtoll(std::string(Tok).c_str(), nullptr, 16);
}

struct FluentNodeZone {
    std::int64_t mFirst, mLast;
    std::size_t mDim;
    std::vector<double> mPoints;
};

struct FluentFace {
    std::int64_t mZone, mBc, mC0, mC1;
    std::vector<std::int64_t> mNodes;
};

struct FluentCellZone {
    std::int64_t mZone, mFirst, mLast;
};

// (nodes, c0, c1) rows of a face body.
void fluent_face_rows(const std::vector<std::int64_t>& rValues, std::size_t Count,
                      std::int64_t FaceType, std::int64_t Zone, std::int64_t Bc,
                      std::vector<FluentFace>& rFaces) {
    const int k = fluent_face_nodes(FaceType);
    std::size_t i = 0;
    for (std::size_t r = 0; r < Count; ++r) {
        const std::size_t n = k ? static_cast<std::size_t>(k)
                                : (i < rValues.size() ? static_cast<std::size_t>(rValues[i++]) : 0);
        if (i + n + 2 > rValues.size())
            throw ReadError("Fluent: face section shorter than its header");
        FluentFace f{Zone, Bc, rValues[i + n], rValues[i + n + 1], {}};
        f.mNodes.assign(rValues.begin() + static_cast<std::ptrdiff_t>(i),
                        rValues.begin() + static_cast<std::ptrdiff_t>(i + n));
        rFaces.push_back(std::move(f));
        i += n + 2;
    }
}

}  // namespace

Mesh read_ansys(const std::string& rPath) {
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        throw ReadError("Could not open file: " + rPath);
    const std::string data((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    FluentReader rd(data);

    std::int64_t dim = 0;
    std::vector<FluentNodeZone> node_zones;
    std::vector<FluentFace> faces;
    std::vector<FluentCellZone> cell_zones;
    struct Legacy {
        std::string mType;
        std::size_t mRows, mCols;
        std::vector<std::int64_t> mConn;
    };
    std::vector<Legacy> legacy;
    std::map<std::int64_t, std::string> names;
    std::vector<std::int64_t> head;

    while (true) {
        rd.SkipWs();
        if (rd.Eof())
            break;
        const std::size_t start = rd.mP;
        if (rd.At() != '(')
            throw ReadError("Fluent: expected a section at byte " + std::to_string(start));
        ++rd.mP;
        rd.SkipWs();
        const std::size_t idx0 = rd.mP;
        while (!rd.Eof() && std::isdigit(static_cast<unsigned char>(rd.At())))
            ++rd.mP;
        const std::string index = data.substr(idx0, rd.mP - idx0);
        if (index.empty())
            throw ReadError("Fluent: expected a section at byte " + std::to_string(start));
        const std::string prefix = index.size() <= 2 ? "" : index.substr(0, index.size() - 2);
        const std::string core = index.size() <= 2 ? index : index.substr(index.size() - 2);

        if (index == "2") {
            rd.SkipWs();
            const std::size_t q = data.find(')', rd.mP);
            dim = std::strtoll(data.substr(rd.mP, q - rd.mP).c_str(), nullptr, 10);
            rd.mP = q == std::string::npos ? data.size() : q + 1;
        } else if (index == "39" || index == "45") {
            // (45 (id type name)(...)): the id is decimal here.
            rd.SkipWs();
            if (!rd.Eof() && rd.At() == '(') {
                const std::size_t q = data.find_first_of("()", rd.mP + 1);
                if (q != std::string::npos) {
                    auto iss =
                        detail::make_classic_istringstream(data.substr(rd.mP + 1, q - rd.mP - 1));
                    std::int64_t id = 0;
                    std::string type, name;
                    if (iss >> id >> type >> name)
                        names[id] = name;
                }
            }
            rd.SkipBalanced(1);
        } else if ((core == "10" || core == "12" || core == "13") &&
                   (prefix.empty() || prefix == "20" || prefix == "30")) {
            const bool has_body = rd.Header(head);
            if (head.size() < 4)
                continue;
            const std::int64_t zone = head[0], first = head[1], last = head[2];
            const std::size_t count =
                static_cast<std::size_t>(std::max<std::int64_t>(0, last - first + 1));
            auto ints = [&](std::size_t n) {
                std::vector<std::int64_t> v;
                if (prefix == "20") {
                    for (std::int32_t x : rd.Binary<std::int32_t>(n))
                        v.push_back(x);
                } else {
                    v = rd.Binary<std::int64_t>(n);
                }
                return v;
            };
            if (core == "10") {
                if (!has_body)
                    continue;
                const std::size_t nd =
                    head.size() > 4 ? static_cast<std::size_t>(head[4]) : (dim ? dim : 3);
                FluentNodeZone z{first, last, nd, {}};
                if (prefix.empty()) {
                    for (std::string_view t : rd.AsciiBody())
                        z.mPoints.push_back(detail::parse_double(std::string(t)));
                } else if (prefix == "20") {
                    for (float x : rd.Binary<float>(count * nd))
                        z.mPoints.push_back(x);
                } else {
                    z.mPoints = rd.Binary<double>(count * nd);
                }
                if (z.mPoints.size() != count * nd)
                    throw ReadError("Fluent: node section size mismatch");
                node_zones.push_back(std::move(z));
            } else if (core == "13") {
                if (!has_body || zone == 0)
                    continue;
                const std::int64_t bc = head[3];
                const std::int64_t face_type = head.size() > 4 ? head[4] : 0;
                std::vector<std::int64_t> values;
                const int k = fluent_face_nodes(face_type);
                if (prefix.empty()) {
                    for (std::string_view t : rd.AsciiBody())
                        values.push_back(fluent_hex(t));
                } else if (k) {
                    values = ints(count * static_cast<std::size_t>(k + 2));
                } else if (prefix == "20") {
                    values = rd.BinaryMixedFaces<std::int32_t>(count);
                } else {
                    values = rd.BinaryMixedFaces<std::int64_t>(count);
                }
                fluent_face_rows(values, count, face_type, zone, bc, faces);
            } else {
                const std::int64_t zone_type = head[3];
                const std::int64_t element_type = head.size() > 4 ? head[4] : 0;
                std::string type;
                int npc = 1;
                const bool fixed = fluent_cell_type(element_type, type, npc);
                std::vector<std::int64_t> body;
                if (has_body) {
                    if (prefix.empty()) {
                        for (std::string_view t : rd.AsciiBody())
                            body.push_back(fluent_hex(t));
                    } else {
                        body = ints(count * static_cast<std::size_t>(npc));
                    }
                }
                if (zone == 0 || zone_type == 0)
                    continue;  // the global declaration, or a dead zone
                if (has_body && fixed) {
                    if (body.size() != count * static_cast<std::size_t>(npc))
                        throw ReadError("Fluent: cell section size mismatch");
                    legacy.push_back({type, count, static_cast<std::size_t>(npc), std::move(body)});
                } else {
                    cell_zones.push_back({zone, first, last});
                }
            }
        } else if (prefix == "20" || prefix == "30") {
            // An unknown binary section: its body is not bracket-balanced.
            std::size_t q = data.find("End of Binary Section", rd.mP);
            q = q == std::string::npos ? q : data.find(')', q);
            rd.mP = q == std::string::npos ? data.size() : q + 1;
        } else {
            rd.SkipBalanced(1);
        }
    }

    if (node_zones.empty())
        throw ReadError("Fluent: no nodes");
    std::int64_t base = node_zones[0].mFirst, top = node_zones[0].mLast;
    std::size_t nd = 0;
    for (const auto& z : node_zones) {
        base = std::min(base, z.mFirst);
        top = std::max(top, z.mLast);
        nd = std::max(nd, z.mDim);
    }
    // Every node, cell and face id range is a count the file declares; none
    // can exceed the file's size in bytes, which bounds the allocations and
    // the loops over the ranges below.
    const std::size_t bytes = detail::file_bytes(rPath);
    const std::size_t npoints = detail::checked_count(top - base + 1, bytes, "Fluent", "node");
    for (const auto& z : node_zones)
        if (z.mLast < z.mFirst ||
            static_cast<std::size_t>(z.mLast - z.mFirst + 1) * z.mDim > z.mPoints.size())
            throw ReadError("Fluent: a node zone holds fewer coordinates than its id range");
    for (const auto& z : cell_zones)
        detail::checked_count(z.mLast - z.mFirst + 1, bytes, "Fluent", "cell");
    NDArray pts(DType::Float64, {npoints, nd});
    double* pp = pts.As<double>();
    std::fill(pp, pp + npoints * nd, 0.0);
    for (const auto& z : node_zones)
        for (std::size_t r = 0; r < static_cast<std::size_t>(z.mLast - z.mFirst + 1); ++r)
            for (std::size_t c = 0; c < z.mDim; ++c)
                pp[(static_cast<std::size_t>(z.mFirst - base) + r) * nd + c] =
                    z.mPoints[r * z.mDim + c];
    face_cells::P3 p3(npoints, {0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < npoints; ++i)
        for (std::size_t c = 0; c < std::min<std::size_t>(nd, 3); ++c)
            p3[i][c] = pp[i * nd + c];

    Mesh mesh;
    mesh.AssignPoints(std::move(pts));

    if (!legacy.empty()) {
        if (!faces.empty())
            log::warn("Fluent: cells with connectivity bodies; the face sections are ignored");
        for (auto& rc : legacy) {
            NDArray conn(DType::Int64, {rc.mRows, rc.mCols});
            std::int64_t* dp = conn.As<std::int64_t>();
            for (std::size_t k = 0; k < rc.mConn.size(); ++k)
                dp[k] = rc.mConn[k] - base;
            mesh.AddCellBlock(rc.mType, std::move(conn));
        }
        return mesh;
    }
    if (dim == 0)
        dim = static_cast<std::int64_t>(nd);

    std::map<std::int64_t, std::int64_t> zone_of;  // cell id -> zone, ascending
    for (const auto& z : cell_zones)
        for (std::int64_t c = z.mFirst; c <= z.mLast; ++c)
            zone_of[c] = z.mZone;
    if (cell_zones.empty())
        for (const auto& f : faces)
            for (std::int64_t c : {f.mC0, f.mC1})
                if (c)
                    zone_of.emplace(c, 0);
    auto live = [&](std::int64_t c) { return zone_of.count(c) > 0; };

    // Outward faces (3-D) or directed edges (2-D) per cell, in file order.
    std::unordered_map<std::int64_t, std::vector<face_cells::Face>> per_cell;
    for (const auto& f : faces) {
        face_cells::Face g(f.mNodes.size());
        for (std::size_t i = 0; i < g.size(); ++i) {
            g[i] = f.mNodes[i] - base;
            if (g[i] < 0 || static_cast<std::size_t>(g[i]) >= npoints)
                throw ReadError("Fluent: a face names a node outside the node zones");
        }
        if (live(f.mC1))
            per_cell[f.mC1].push_back(g);
        if (live(f.mC0))
            per_cell[f.mC0].emplace_back(g.rbegin(), g.rend());
    }

    // (zone, type) buckets in first-seen order.
    struct Bucket {
        std::int64_t mZone;
        std::string mType;
        std::vector<face_cells::Face> mConn;
        std::vector<std::vector<face_cells::Face>> mPoly;
    };
    std::vector<Bucket> volume, surface;
    auto bucket = [](std::vector<Bucket>& rList, std::int64_t Zone,
                     const std::string& rType) -> Bucket& {
        for (auto& b : rList)
            if (b.mZone == Zone && b.mType == rType)
                return b;
        rList.push_back({Zone, rType, {}, {}});
        return rList.back();
    };
    auto face_type_name = [](std::size_t N, bool Surface) -> std::string {
        if (N == 2 && Surface)
            return "line";
        if (N == 3)
            return "triangle";
        if (N == 4)
            return "quad";
        return "polygon" + std::to_string(N);
    };

    std::size_t skipped = 0;
    for (const auto& [cid, zone] : zone_of) {
        auto it = per_cell.find(cid);
        if (it == per_cell.end() || it->second.empty()) {
            ++skipped;
            continue;
        }
        const auto& cf = it->second;
        if (dim == 2) {
            std::vector<std::array<std::int64_t, 2>> edges;
            for (const auto& f : cf)
                if (!f.empty())
                    edges.push_back({f.front(), f.back()});
            face_cells::Face ring = face_cells::polygon_from_edges(edges, p3);
            if (ring.empty()) {
                ++skipped;
                continue;
            }
            bucket(volume, zone, face_type_name(ring.size(), false)).mConn.push_back(ring);
            continue;
        }
        auto [type, conn] = face_cells::reconstruct_cell(cf, p3);
        if (type == "polyhedron") {
            const std::string key =
                "polyhedron" + std::to_string(face_cells::unique_node_count(cf));
            bucket(volume, zone, key).mPoly.push_back(cf);
        } else if (conn.empty()) {
            ++skipped;
        } else {
            bucket(volume, zone, type).mConn.push_back(std::move(conn));
        }
    }
    if (skipped)
        log::warn("Fluent: {} cell(s) with no usable faces skipped", skipped);

    for (const auto& f : faces) {
        if (f.mBc == kFluentInterior && live(f.mC0) && live(f.mC1))
            continue;
        face_cells::Face g(f.mNodes.size());
        for (std::size_t i = 0; i < g.size(); ++i) {
            g[i] = f.mNodes[i] - base;
            if (g[i] < 0 || static_cast<std::size_t>(g[i]) >= npoints)
                throw ReadError("Fluent: a face names a node outside the node zones");
        }
        // Outward from the domain: the normal points into c0.
        if (live(f.mC0) && !live(f.mC1))
            std::reverse(g.begin(), g.end());
        bucket(surface, f.mZone, face_type_name(g.size(), true)).mConn.push_back(std::move(g));
    }

    std::vector<NDArray> zone_data;
    // (zone, dim) -> global cell ids, first-seen order
    std::vector<std::pair<std::pair<std::int64_t, int>, std::vector<std::int64_t>>> regions;
    std::int64_t offset = 0;
    for (auto* pGroup : {&volume, &surface}) {
        const int cdim = static_cast<int>(pGroup == &volume ? dim : dim - 1);
        for (auto& b : *pGroup) {
            std::size_t n;
            if (b.mType.rfind("polyhedron", 0) == 0) {
                n = b.mPoly.size();
                mesh.AddPolyhedronBlock(b.mType, std::move(b.mPoly));
            } else {
                n = b.mConn.size();
                const std::size_t k = n ? b.mConn[0].size() : 0;
                NDArray conn(DType::Int64, {n, k});
                std::int64_t* c = conn.As<std::int64_t>();
                for (std::size_t r = 0; r < n; ++r)
                    std::copy(b.mConn[r].begin(), b.mConn[r].end(), c + r * k);
                mesh.AddCellBlock(b.mType, std::move(conn));
            }
            NDArray z(DType::Int64, {n});
            std::fill(z.As<std::int64_t>(), z.As<std::int64_t>() + n, b.mZone);
            zone_data.push_back(std::move(z));
            auto key = std::make_pair(b.mZone, cdim);
            auto rit = std::find_if(regions.begin(), regions.end(),
                                    [&](const auto& r) { return r.first == key; });
            if (rit == regions.end()) {
                regions.push_back({key, {}});
                rit = regions.end() - 1;
            }
            for (std::size_t r = 0; r < n; ++r)
                rit->second.push_back(offset + static_cast<std::int64_t>(r));
            offset += static_cast<std::int64_t>(n);
        }
    }
    if (!zone_data.empty())
        mesh.AddCellData("ansys:zone", std::move(zone_data));
    for (auto& [key, ids] : regions) {
        auto nit = names.find(key.first);
        const std::string name =
            nit != names.end() ? nit->second : "zone_" + std::to_string(key.first);
        NDArray e(DType::Int64, {ids.size()});
        std::copy(ids.begin(), ids.end(), e.As<std::int64_t>());
        mesh.AddRegion(Region(name, RegionKind::Cell, key.second, key.first, std::move(e)));
    }
    return mesh;
}

namespace {

// One face of the written mesh: its nodes (0-based) wound so Fluent reads the
// right cells on either side, and the compact ids of those cells (-1: none).
struct FlFace {
    std::vector<std::int64_t> mNodes;
    std::int64_t mC0 = -1;
    std::int64_t mC1 = -1;
};

struct FlZone {
    std::int64_t mId = 0;
    std::string mType;                  // fluid, interior, wall
    std::vector<std::int64_t> mGlobal;  // the mesh cells it stands for (for its name)
    int mDim = 0;                       // their dimension
    std::vector<std::size_t> mMembers;  // compact cells or faces, in order
};

// Fluent element types of the cell zones (12) and the corners a cell keeps.
int fl_element_type(const std::string& rType) {
    if (rType.rfind("triangle", 0) == 0)
        return 1;
    if (rType.rfind("tetra", 0) == 0)
        return 2;
    if (rType.rfind("quad", 0) == 0)
        return 3;
    if (rType.rfind("hexahedron", 0) == 0)
        return 4;
    if (rType.rfind("pyramid", 0) == 0)
        return 5;
    if (rType.rfind("wedge", 0) == 0)
        return 6;
    return 7;  // polyhedral (3-D) or polygonal (2-D): defined by its faces
}

bool fl_is_linear(const std::string& rType) {
    static const char* const kLinear[] = {"triangle", "tetra", "quad",    "hexahedron",
                                          "pyramid",  "wedge", "polygon", "line"};
    for (const char* t : kLinear)
        if (rType == t)
            return true;
    return rType.rfind("polygon", 0) == 0 || rType.rfind("polyhedron", 0) == 0;
}

// The corner ring of a 2-D cell (a surface facet in 3-D, a cell in 2-D).
std::vector<std::int64_t> fl_ring(const Mesh::CellView& rBlock, std::size_t Cell) {
    std::vector<std::int64_t> ring;
    if (rBlock.IsRagged()) {
        const std::int64_t* p = rBlock.Row(Cell);
        ring.assign(p, p + rBlock.RowSize(Cell));
        return ring;
    }
    const std::string& t = rBlock.Type();
    const std::size_t npc = rBlock.NodesPerCell();
    std::size_t corners = npc;
    if (t.rfind("triangle", 0) == 0)
        corners = 3;
    else if (t.rfind("quad", 0) == 0)
        corners = 4;
    else if (t.rfind("line", 0) == 0)
        corners = 2;
    const NDArray& conn = rBlock.Conn();
    for (std::size_t k = 0; k < corners && k < npc; ++k)
        ring.push_back(detail::read_int(conn, Cell * npc + k));
    return ring;
}

int fl_dimension(const Mesh::CellView& rBlock) {
    if (rBlock.IsPolyhedron())
        return 3;
    if (rBlock.IsRagged() || rBlock.Type().rfind("polygon", 0) == 0)
        return 2;
    return cell_type_dimension(cell_type_from_name(rBlock.Type()));
}

std::string fl_hex(std::int64_t V) {
    char buf[24];
    std::snprintf(buf, sizeof(buf), "%llx", static_cast<unsigned long long>(V));
    return buf;
}

// A Fluent zone name: no blanks or parentheses.
std::string fl_name(std::string rName) {
    for (char& c : rName)
        if (c == ' ' || c == '\t' || c == '(' || c == ')' || c == '"')
            c = '_';
    return rName.empty() ? std::string("zone") : rName;
}

}  // namespace

void write_ansys(const std::string& rPath, const Mesh& rMesh, bool binary) {
    const std::size_t npoints = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = rMesh.PointDim();
    if (pdim != 2 && pdim != 3)
        throw WriteError("Fluent: can only write points of dimension 2 or 3");

    // The cells are the blocks of the mesh's highest dimension (2 or 3); the
    // blocks one dimension lower name boundary zones; anything else is dropped.
    int dim = 0;
    for (const auto cb : rMesh.CellRange())
        dim = std::max(dim, fl_dimension(cb));
    if (dim < 2)
        throw WriteError("Fluent: the mesh has no 2-D or 3-D cells");
    if (dim == 3 && pdim != 3)
        throw WriteError("Fluent: 3-D cells need 3-D points");
    // A 2-D Fluent mesh lies in the xy plane: z is dropped only when it is 0.
    if (dim == 2 && pdim == 3)
        for (std::size_t i = 0; i < npoints; ++i)
            if (detail::read_double(points, i * 3 + 2) != 0.0)
                throw WriteError("Fluent: a 2-D mesh must lie in the z = 0 plane (point " +
                                 std::to_string(i) +
                                 " has z != 0); Fluent has no 3-D surface meshes");
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const bool has_zone = rMesh.HasCellData("ansys:zone");
    auto zone_value = [&](std::size_t Block, std::size_t Cell, std::int64_t& rOut) {
        if (!has_zone)
            return false;
        const NDArray& z = rMesh.CellData("ansys:zone", Block);
        if (Cell >= z.Size())
            return false;
        rOut = detail::read_int(z, Cell);
        return rOut > 0;
    };

    // ---- faces and cells ----------------------------------------------------
    std::vector<FlFace> faces;
    std::vector<std::int64_t> cell_to_global;  // compact cell -> global cell
    std::vector<std::size_t> cell_block;       // compact cell -> block
    std::vector<std::size_t> surface_blocks;   // blocks of dimension dim - 1
    std::size_t dropped_blocks = 0, quadratic = 0;
    if (dim == 3) {
        const detail::GlobalFaces g = detail::build_global_faces(rMesh);
        if (g.mNumUnorientable != 0)
            log::warn(
                "Fluent: {} cell(s) are not closed, orientable solids; their faces are "
                "written as found",
                g.mNumUnorientable);
        if (g.mNumNonManifold != 0)
            log::warn("Fluent: {} face(s) are shared by more than two cells", g.mNumNonManifold);
        cell_to_global = g.mCellToGlobal;
        cell_block.resize(g.NumCells());
        for (std::size_t c = 0; c < g.NumCells(); ++c)
            cell_block[c] = static_cast<std::size_t>(
                std::upper_bound(bases.begin(), bases.end(), g.mCellToGlobal[c]) - bases.begin() -
                1);
        faces.resize(g.NumFaces());
        for (std::size_t f = 0; f < g.NumFaces(); ++f) {
            // Stored outward from the owner; Fluent's normal points into c0.
            const std::int64_t* ring = g.Face(f);
            faces[f].mNodes.assign(std::make_reverse_iterator(ring + g.FaceSize(f)),
                                   std::make_reverse_iterator(ring));
            faces[f].mC0 = g.mOwner[f];
            faces[f].mC1 = g.mNeighbour[f];
        }
        for (std::size_t b : g.mNonCellBlocks)
            (fl_dimension(rMesh.Cells(b)) == 2 ? surface_blocks.push_back(b)
                                               : void(++dropped_blocks));
    } else {
        // Each 2-D cell counter-clockwise (in x-y); an edge a -> b of its ring
        // has it on the left, where Fluent puts c0.
        std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> edge_of;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const auto cb = rMesh.Cells(b);
            const int d = fl_dimension(cb);
            if (d != 2) {
                (d == 1 ? surface_blocks.push_back(b) : void(++dropped_blocks));
                continue;
            }
            for (std::size_t i = 0; i < cb.NumCells(); ++i) {
                std::vector<std::int64_t> ring = fl_ring(cb, i);
                double area = 0.0;
                for (std::size_t k = 0; k < ring.size(); ++k) {
                    const std::size_t p = static_cast<std::size_t>(ring[k]);
                    const std::size_t q = static_cast<std::size_t>(ring[(k + 1) % ring.size()]);
                    area += detail::read_double(points, p * pdim) *
                                detail::read_double(points, q * pdim + 1) -
                            detail::read_double(points, q * pdim) *
                                detail::read_double(points, p * pdim + 1);
                }
                if (area < 0.0)
                    std::reverse(ring.begin(), ring.end());
                const std::int64_t cell = static_cast<std::int64_t>(cell_to_global.size());
                cell_to_global.push_back(bases[b] + static_cast<std::int64_t>(i));
                cell_block.push_back(b);
                for (std::size_t k = 0; k < ring.size(); ++k) {
                    const std::int64_t a = ring[k], e = ring[(k + 1) % ring.size()];
                    const auto key = std::minmax(a, e);
                    const auto it = edge_of.find(key);
                    if (it == edge_of.end()) {
                        edge_of.emplace(key, faces.size());
                        faces.push_back({{a, e}, cell, -1});
                    } else if (faces[it->second].mC1 < 0) {
                        faces[it->second].mC1 = cell;
                    }
                }
            }
        }
    }
    if (cell_to_global.empty())
        throw WriteError("Fluent: the mesh has no cells Fluent can hold");
    if (dropped_blocks != 0)
        log::warn("Fluent: {} cell block(s) of other dimensions are not written", dropped_blocks);
    for (std::size_t b : std::set<std::size_t>(cell_block.begin(), cell_block.end()))
        quadratic += fl_is_linear(rMesh.Cells(b).Type()) ? 0 : 1;
    if (quadratic != 0)
        log::warn("Fluent: cells are linear; the mid-side nodes of {} block(s) are not used",
                  quadratic);

    // ---- zones ---------------------------------------------------------------
    // ansys:zone values are kept as zone ids (the reader writes them); a block
    // without one, the interior faces, the leftover boundary faces and the
    // nodes get fresh ids.
    std::set<std::int64_t> used;
    std::vector<std::int64_t> cell_zone_id(cell_to_global.size(), 0);
    for (std::size_t c = 0; c < cell_to_global.size(); ++c) {
        std::int64_t z = 0;
        if (zone_value(cell_block[c],
                       static_cast<std::size_t>(cell_to_global[c] - bases[cell_block[c]]), z)) {
            cell_zone_id[c] = z;
            used.insert(z);
        }
    }
    const std::set<std::int64_t> cell_ids_used = used;
    std::map<std::pair<std::size_t, std::size_t>, std::int64_t> surface_explicit;  // (block, cell)
    for (std::size_t b : surface_blocks) {
        const auto cb = rMesh.Cells(b);
        for (std::size_t i = 0; i < cb.NumCells(); ++i) {
            std::int64_t z = 0;
            if (zone_value(b, i, z) && cell_ids_used.count(z) == 0) {
                surface_explicit[{b, i}] = z;
                used.insert(z);
            }
        }
    }
    std::int64_t next = used.empty() ? 1 : *used.rbegin() + 1;
    auto fresh = [&]() {
        while (used.count(next))
            ++next;
        used.insert(next);
        return next++;
    };
    std::map<std::size_t, std::int64_t> block_zone;  // a block without ansys:zone
    for (std::size_t c = 0; c < cell_to_global.size(); ++c)
        if (cell_zone_id[c] == 0) {
            auto it = block_zone.find(cell_block[c]);
            if (it == block_zone.end())
                it = block_zone.emplace(cell_block[c], fresh()).first;
            cell_zone_id[c] = it->second;
        }

    // Cell zones in first-seen order; Fluent numbers cells zone by zone.
    std::vector<FlZone> cell_zones;
    std::map<std::int64_t, std::size_t> cell_zone_pos;
    for (std::size_t c = 0; c < cell_to_global.size(); ++c) {
        auto it = cell_zone_pos.find(cell_zone_id[c]);
        if (it == cell_zone_pos.end()) {
            it = cell_zone_pos.emplace(cell_zone_id[c], cell_zones.size()).first;
            cell_zones.push_back({cell_zone_id[c], "fluid", {}, dim, {}});
        }
        cell_zones[it->second].mMembers.push_back(c);
        cell_zones[it->second].mGlobal.push_back(cell_to_global[c]);
    }
    std::vector<std::int64_t> fluent_cell(cell_to_global.size(), 0);
    {
        std::int64_t id = 1;
        for (const FlZone& z : cell_zones)
            for (std::size_t c : z.mMembers)
                fluent_cell[c] = id++;
    }

    // Boundary faces: the facet cell that matches one names its zone.
    std::vector<std::int64_t> face_zone(faces.size(), 0);
    std::vector<FlZone> face_zones;
    face_zones.push_back({fresh(), "interior", {}, dim - 1, {}});
    std::map<std::int64_t, std::size_t> face_zone_pos;
    std::unordered_map<std::string, std::size_t> face_by_key;
    auto key_of = [](std::vector<std::int64_t> rIds) {
        std::sort(rIds.begin(), rIds.end());
        std::string k;
        for (std::int64_t v : rIds)
            k += std::to_string(v) + ",";
        return k;
    };
    for (std::size_t f = 0; f < faces.size(); ++f)
        if (faces[f].mC1 < 0)
            face_by_key.emplace(key_of(faces[f].mNodes), f);
    std::size_t unmatched = 0;
    for (std::size_t b : surface_blocks) {
        const auto cb = rMesh.Cells(b);
        std::int64_t block_id = 0;
        for (std::size_t i = 0; i < cb.NumCells(); ++i) {
            const auto it = face_by_key.find(key_of(fl_ring(cb, i)));
            if (it == face_by_key.end()) {
                ++unmatched;
                continue;
            }
            const std::size_t f = it->second;
            if (face_zone[f] != 0)
                continue;
            std::int64_t z = 0;
            const auto e = surface_explicit.find({b, i});
            if (e != surface_explicit.end()) {
                z = e->second;
            } else {
                if (block_id == 0)
                    block_id = fresh();
                z = block_id;
            }
            face_zone[f] = z;
            auto pos = face_zone_pos.find(z);
            if (pos == face_zone_pos.end()) {
                pos = face_zone_pos.emplace(z, face_zones.size()).first;
                face_zones.push_back({z, "wall", {}, dim - 1, {}});
            }
            face_zones[pos->second].mGlobal.push_back(bases[b] + static_cast<std::int64_t>(i));
        }
    }
    if (unmatched != 0)
        log::warn(
            "Fluent: {} facet cell(s) are not on the boundary of the cells and are not "
            "written",
            unmatched);
    std::int64_t default_wall = 0;
    for (std::size_t f = 0; f < faces.size(); ++f) {
        if (faces[f].mC1 >= 0) {
            face_zones[0].mMembers.push_back(f);
            continue;
        }
        if (face_zone[f] == 0) {
            if (default_wall == 0) {
                default_wall = fresh();
                face_zone_pos.emplace(default_wall, face_zones.size());
                face_zones.push_back({default_wall, "wall", {}, dim - 1, {}});
            }
            face_zone[f] = default_wall;
        }
        face_zones[face_zone_pos.at(face_zone[f])].mMembers.push_back(f);
    }
    const std::int64_t node_zone = fresh();

    // Zone names: a region the reader would have made (same tag and
    // dimension), else a region with exactly these cells, else type_<id>.
    auto zone_name = [&](const FlZone& rZ) {
        std::vector<std::int64_t> want = rZ.mGlobal;
        std::sort(want.begin(), want.end());
        for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
            const meshioplusplus::Region& reg = rMesh.Region(r);
            if (reg.mKind == RegionKind::Cell && reg.mTag == rZ.mId && reg.mDim == rZ.mDim)
                return fl_name(reg.mName);
        }
        if (!want.empty())
            for (std::size_t r = 0; r < rMesh.NumRegions(); ++r) {
                const meshioplusplus::Region& reg = rMesh.Region(r);
                if (reg.mKind != RegionKind::Cell || reg.NumEntries() != want.size())
                    continue;
                if (std::equal(want.begin(), want.end(), reg.Entries()))
                    return fl_name(reg.mName);
            }
        return rZ.mType + "_" + std::to_string(rZ.mId);
    };

    if (rMesh.NumPointData() != 0 || rMesh.NumCellData() > (has_zone ? 1u : 0u) ||
        rMesh.NumFieldData() != 0)
        detail::provenance_note("data-dropped", "a Fluent mesh file holds no data arrays");

    // ---- write -------------------------------------------------------------------
    auto fh = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!fh)
        throw WriteError("Could not open file for writing: " + rPath);
    const std::size_t odim = static_cast<std::size_t>(dim);
    const std::size_t ncells = cell_to_global.size();
    fh << "(1 \"" << detail::provenance_lines(detail::SlotTier::SingleLine)[0] << "\")\n";
    fh << "(2 " << odim << ")\n";
    fh << "(10 (0 1 " << fl_hex(static_cast<std::int64_t>(npoints)) << " 0 " << odim << "))\n";
    fh << "(13 (0 1 " << fl_hex(static_cast<std::int64_t>(faces.size())) << " 0))\n";
    fh << "(12 (0 1 " << fl_hex(static_cast<std::int64_t>(ncells)) << " 0))\n";

    // Nodes.
    fh << "(" << (binary ? "3010" : "10") << " (" << fl_hex(node_zone) << " 1 "
       << fl_hex(static_cast<std::int64_t>(npoints)) << " 1 " << odim << ")"
       << (binary ? "\n(" : "(");
    if (binary) {
        for (std::size_t i = 0; i < npoints; ++i)
            for (std::size_t c = 0; c < odim; ++c) {
                const double v = detail::read_double(points, i * pdim + c);
                fh.write(reinterpret_cast<const char*>(&v), 8);
            }
        fh << ")\nEnd of Binary Section 3010)\n";
    } else {
        fh << "\n";
        char buf[32];
        for (std::size_t i = 0; i < npoints; ++i) {
            for (std::size_t c = 0; c < odim; ++c) {
                detail::snprintf_c(buf, sizeof(buf), "%.16e",
                                   detail::read_double(points, i * pdim + c));
                fh << buf << (c + 1 == odim ? "\n" : " ");
            }
        }
        fh << "))\n";
    }

    auto write_ints = [&](const char* pKey, const std::string& rHead,
                          const std::vector<std::vector<std::int64_t>>& rRows) {
        fh << "(" << (binary ? std::string("20") + pKey : std::string(pKey)) << " (" << rHead << ")"
           << (binary ? "\n(" : "(");
        if (binary) {
            for (const auto& row : rRows)
                for (std::int64_t v : row) {
                    if (v > std::numeric_limits<std::int32_t>::max())
                        throw WriteError("Fluent: an id does not fit a binary 32-bit section");
                    const std::int32_t w = static_cast<std::int32_t>(v);
                    fh.write(reinterpret_cast<const char*>(&w), 4);
                }
            fh << ")\nEnd of Binary Section 20" << pKey << ")\n";
            return;
        }
        fh << "\n";
        for (const auto& row : rRows) {
            for (std::size_t k = 0; k < row.size(); ++k)
                fh << fl_hex(row[k]) << (k + 1 == row.size() ? "\n" : " ");
        }
        fh << "))\n";
    };

    // Cell zones: one element type in the header, or 0 and a type per cell.
    std::int64_t first = 1;
    for (const FlZone& z : cell_zones) {
        std::vector<int> types;
        for (std::size_t c : z.mMembers)
            types.push_back(fl_element_type(rMesh.Cells(cell_block[c]).Type()));
        const bool mixed =
            std::any_of(types.begin(), types.end(), [&](int t) { return t != types.front(); });
        const std::int64_t last = first + static_cast<std::int64_t>(z.mMembers.size()) - 1;
        const std::string head = fl_hex(z.mId) + " " + fl_hex(first) + " " + fl_hex(last) + " 1 " +
                                 fl_hex(mixed ? 0 : types.front());
        if (!mixed) {
            fh << "(12 (" << head << "))\n";
        } else {
            std::vector<std::vector<std::int64_t>> rows;
            for (int t : types)
                rows.push_back({t});
            write_ints("12", head, rows);
        }
        first = last + 1;
    }

    // Face zones: interior (bc 2) first, then the walls (bc 3). A zone of one
    // face size lists bare nodes; a mixed one leads each face with its size.
    first = 1;
    for (const FlZone& z : face_zones) {
        if (z.mMembers.empty())
            continue;
        const std::size_t size0 = faces[z.mMembers.front()].mNodes.size();
        bool uniform = true;
        std::size_t largest = 0;
        for (std::size_t f : z.mMembers) {
            uniform = uniform && faces[f].mNodes.size() == size0;
            largest = std::max(largest, faces[f].mNodes.size());
        }
        const int ftype =
            uniform && size0 >= 2 && size0 <= 4 ? static_cast<int>(size0) : (largest > 4 ? 5 : 0);
        const std::int64_t last = first + static_cast<std::int64_t>(z.mMembers.size()) - 1;
        const std::string head = fl_hex(z.mId) + " " + fl_hex(first) + " " + fl_hex(last) + " " +
                                 (z.mType == "interior" ? "2" : "3") + " " + fl_hex(ftype);
        std::vector<std::vector<std::int64_t>> rows;
        rows.reserve(z.mMembers.size());
        for (std::size_t f : z.mMembers) {
            std::vector<std::int64_t> row;
            if (ftype == 0 || ftype == 5)
                row.push_back(static_cast<std::int64_t>(faces[f].mNodes.size()));
            for (std::int64_t n : faces[f].mNodes)
                row.push_back(n + 1);
            row.push_back(fluent_cell[static_cast<std::size_t>(faces[f].mC0)]);
            row.push_back(faces[f].mC1 < 0 ? 0
                                           : fluent_cell[static_cast<std::size_t>(faces[f].mC1)]);
            rows.push_back(std::move(row));
        }
        write_ints("13", head, rows);
        first = last + 1;
    }

    for (const FlZone& z : cell_zones)
        fh << "(45 (" << z.mId << " fluid " << zone_name(z) << ")())\n";
    for (const FlZone& z : face_zones)
        if (!z.mMembers.empty())
            fh << "(45 (" << z.mId << " " << z.mType << " " << zone_name(z) << ")())\n";
    if (!fh)
        throw WriteError("Fluent: failed writing " + rPath);
}

}  // namespace meshioplusplus
