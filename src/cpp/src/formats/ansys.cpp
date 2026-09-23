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
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/ansys.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
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
        if (mP > mD.size() || size > mD.size() - mP)
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
            if (n < 0 || row * sizeof(T) > mD.size() - mP)
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
    const std::size_t npoints = static_cast<std::size_t>(top - base + 1);
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
        for (std::size_t i = 0; i < g.size(); ++i)
            g[i] = f.mNodes[i] - base;
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
        for (std::size_t i = 0; i < g.size(); ++i)
            g[i] = f.mNodes[i] - base;
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

void write_ansys(const std::string& rPath, const Mesh& rMesh, bool binary) {
    auto fh = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!fh)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t npoints = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    const std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 0;
    if (dim != 2 && dim != 3)
        throw WriteError("ANSYS: can only write dimension 2 or 3");

    static const std::unordered_map<std::string, int> meshio_to_ansys = {
        {"triangle", 1},   {"tetra", 2},   {"quad", 3},
        {"hexahedron", 4}, {"pyramid", 5}, {"wedge", 6}};

    char hbuf[128];
    fh << "(1 \"" << detail::provenance_lines(detail::SlotTier::SingleLine)[0] << "\")\n";
    std::snprintf(hbuf, sizeof(hbuf), "(2 %zu)\n", dim);
    fh << hbuf;

    const std::size_t first_node_index = 1;
    std::snprintf(hbuf, sizeof(hbuf), "(10 (0 %zx %zx 0))\n", first_node_index, npoints);
    fh << hbuf;

    std::size_t total_cells = 0;
    for (const auto cb : rMesh.CellRange())
        total_cells += cb.NumCells();
    std::snprintf(hbuf, sizeof(hbuf), "(12 (0 1 %zx 0))\n", total_cells);
    fh << hbuf;

    // Nodes
    const char* nkey = binary ? "3010" : "10";
    std::snprintf(hbuf, sizeof(hbuf), "(%s (1 %zx %zx 1 %zx)(\n", nkey, first_node_index, npoints,
                  dim);
    fh << hbuf;
    if (binary) {
        for (std::size_t i = 0; i < npoints; ++i)
            for (std::size_t c = 0; c < dim; ++c) {
                double v = detail::read_double(points, i * dim + c);
                fh.write(reinterpret_cast<const char*>(&v), 8);
            }
        fh << "\n)";
        fh << "End of Binary Section 3010)\n";
    } else {
        char cbuf[32];
        for (std::size_t i = 0; i < npoints; ++i) {
            for (std::size_t c = 0; c < dim; ++c) {
                detail::snprintf_c(cbuf, sizeof(cbuf), "%.16e",
                                   detail::read_double(points, i * dim + c));
                fh << cbuf << (c + 1 == dim ? "" : " ");
            }
            fh << "\n";
        }
        fh << "))\n";
    }

    // Cells
    std::size_t first_index = 0;
    for (const auto cb : rMesh.CellRange()) {
        auto it = meshio_to_ansys.find(cb.Type());
        if (it == meshio_to_ansys.end())
            throw WriteError("ANSYS: illegal cell type '" + cb.Type() + "'");
        int ansys_type = it->second;
        std::size_t n = cb.NumCells();
        const NDArray& conn = cb.Conn();
        std::size_t ncols = detail::cols(conn);
        std::size_t last_index = first_index + n - 1;
        bool is_i32 = (conn.Dtype() == DType::Int32);
        const char* ckey = binary ? (is_i32 ? "2012" : "3012") : "12";
        std::snprintf(hbuf, sizeof(hbuf), "(%s (1 %zx %zx 1 %d)(\n", ckey, first_index, last_index,
                      ansys_type);
        fh << hbuf;
        if (binary) {
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t c = 0; c < ncols; ++c) {
                    std::int64_t v = detail::read_int(conn, r * ncols + c) + 1;
                    if (is_i32) {
                        std::int32_t v32 = static_cast<std::int32_t>(v);
                        fh.write(reinterpret_cast<const char*>(&v32), 4);
                    } else
                        fh.write(reinterpret_cast<const char*>(&v), 8);
                }
            fh << "\n)";
            std::snprintf(hbuf, sizeof(hbuf), "End of Binary Section %s)\n", ckey);
            fh << hbuf;
        } else {
            char cbuf[24];
            for (std::size_t r = 0; r < n; ++r) {
                for (std::size_t c = 0; c < ncols; ++c) {
                    std::snprintf(
                        cbuf, sizeof(cbuf), "%llx",
                        static_cast<unsigned long long>(detail::read_int(conn, r * ncols + c) + 1));
                    fh << cbuf << (c + 1 == ncols ? "" : " ");
                }
                fh << "\n";
            }
            fh << "))\n";
        }
        first_index = last_index + 1;
    }
}

}  // namespace meshioplusplus
