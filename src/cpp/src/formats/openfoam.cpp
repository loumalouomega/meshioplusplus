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
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/openfoam.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"

namespace fs = std::filesystem;

namespace meshioplusplus {

namespace {

using Face = std::vector<std::int64_t>;

struct FoamFormat {
    bool mBinary = false;
    int mLabelBytes = 8;
    int mScalarBytes = 8;
};

/**
 * @brief Whole-file access, mapped where that pays (detail/file_source.hpp).
 *
 * This replaces an `ostringstream` + `.str()` slurp, which paid for **two**
 * extra full-file copies on top of the read -- by far the worst of the
 * whole-file readers, and the reason this one benefits most from mapping.
 * Returns the source itself so the caller controls its lifetime; everything
 * below takes a view into it.
 */
detail::FileSource read_whole(const std::string& rPath) {
    try {
        return detail::FileSource(rPath);
    } catch (const ReadError&) {
        throw ReadError("Could not open OpenFOAM file: " + rPath);
    }
}

std::string openfoam_strip(const std::string& rS) {
    std::size_t a = rS.find_first_not_of(" \t\r\n");
    if (a == std::string::npos)
        return "";
    std::size_t b = rS.find_last_not_of(" \t\r\n");
    return rS.substr(a, b - a + 1);
}

// Parse the FoamFile header for format/arch (label/scalar byte widths).
FoamFormat detect_format(const std::string& rPath) {
    FoamFormat fmt;
    std::ifstream f(rPath, std::ios::binary);
    if (!f)
        return fmt;
    std::string line;
    while (std::getline(f, line)) {
        std::string s = openfoam_strip(line);
        // format <word>;
        std::size_t p = s.find("format");
        if (p == 0) {
            std::string rest = openfoam_strip(s.substr(6));
            if (!rest.empty() && rest.back() == ';')
                rest.pop_back();
            rest = openfoam_strip(rest);
            if (rest == "binary")
                fmt.mBinary = true;
            else if (rest == "ascii")
                fmt.mBinary = false;
        }
        if (s.rfind("arch", 0) == 0) {
            std::size_t lp = s.find("label=");
            if (lp != std::string::npos) {
                int bits = std::atoi(s.c_str() + lp + 6);
                if (bits)
                    fmt.mLabelBytes = bits / 8;
            }
            std::size_t sp = s.find("scalar=");
            if (sp != std::string::npos) {
                int bits = std::atoi(s.c_str() + sp + 7);
                if (bits)
                    fmt.mScalarBytes = bits / 8;
            }
        }
        if (s == "}")
            break;
    }
    return fmt;
}

// Strip C-style /* */ and // comments and drop the FoamFile { ... } block.
std::string strip_comments_and_header(std::string_view rText) {
    std::string out;
    out.reserve(rText.size());
    // remove /* */ and //
    for (std::size_t i = 0; i < rText.size();) {
        if (i + 1 < rText.size() && rText[i] == '/' && rText[i + 1] == '*') {
            std::size_t e = rText.find("*/", i + 2);
            i = (e == std::string::npos) ? rText.size() : e + 2;
        } else if (i + 1 < rText.size() && rText[i] == '/' && rText[i + 1] == '/') {
            std::size_t e = rText.find('\n', i + 2);
            i = (e == std::string::npos) ? rText.size() : e;
        } else {
            out.push_back(rText[i++]);
        }
    }
    // drop FoamFile { ... }
    std::istringstream ss(out);
    std::string line, result;
    bool in_header = false;
    int depth = 0;
    while (std::getline(ss, line)) {
        std::string s = openfoam_strip(line);
        if (s.find("FoamFile") != std::string::npos)
            in_header = true;
        if (in_header) {
            for (char c : s) {
                if (c == '{')
                    ++depth;
                else if (c == '}')
                    --depth;
            }
            if (depth <= 0)
                in_header = false;
            continue;
        }
        result += line;
        result.push_back('\n');
    }
    return result;
}

// ---- ASCII parsers ----

std::vector<std::array<double, 3>> parse_points_ascii(const std::string& rBody) {
    std::vector<std::array<double, 3>> pts;
    std::istringstream ss(rBody);
    std::string line;
    bool in_block = false;
    bool have_n = false;
    while (std::getline(ss, line)) {
        std::string s = openfoam_strip(line);
        if (s.empty())
            continue;
        if (!have_n && s.find_first_not_of("0123456789") == std::string::npos) {
            have_n = true;
            continue;
        }
        if (s == "(" && have_n) {
            in_block = true;
            continue;
        }
        if (s == ")" && in_block)
            break;
        if (in_block) {
            // extract up to 3 numbers from within parentheses
            std::string t = s;
            for (char& c : t)
                if (c == '(' || c == ')')
                    c = ' ';
            std::istringstream ns(t);
            double a, b, c;
            if (ns >> a >> b >> c)
                pts.push_back({a, b, c});
        }
    }
    return pts;
}

std::vector<Face> parse_faces_ascii(const std::string& rBody) {
    std::vector<Face> faces;
    std::istringstream ss(rBody);
    std::string line;
    bool in_block = false, have_n = false;
    while (std::getline(ss, line)) {
        std::string s = openfoam_strip(line);
        if (s.empty())
            continue;
        if (!have_n && s.find_first_not_of("0123456789") == std::string::npos) {
            have_n = true;
            continue;
        }
        if (s == "(" && have_n) {
            in_block = true;
            continue;
        }
        if (s == ")" && in_block)
            break;
        if (in_block) {
            // form: <count>(<ids...>)
            std::size_t lp = s.find('(');
            std::size_t rp = s.find(')', lp);
            if (lp == std::string::npos || rp == std::string::npos)
                continue;
            std::string inside = s.substr(lp + 1, rp - lp - 1);
            std::istringstream ns(inside);
            Face f;
            std::int64_t v;
            while (ns >> v)
                f.push_back(v);
            faces.push_back(std::move(f));
        }
    }
    return faces;
}

std::vector<std::int64_t> parse_int_list_ascii(const std::string& rBody) {
    std::vector<std::int64_t> out;
    std::istringstream ss(rBody);
    std::string line;
    bool in_block = false, have_n = false;
    while (std::getline(ss, line)) {
        std::string s = openfoam_strip(line);
        if (s.empty())
            continue;
        if (!have_n && s.find_first_not_of("0123456789") == std::string::npos) {
            have_n = true;
            continue;
        }
        if (s == "(") {
            in_block = true;
            continue;
        }
        if (s == ")")
            break;
        if (in_block) {
            std::istringstream ns(s);
            std::int64_t v;
            while (ns >> v)
                out.push_back(v);
        }
    }
    return out;
}

// Boundary patch descriptor. `mNFaces`/`mStartFace` deliberately mirror
// OpenFOAM's own on-disk `boundary` field names (`nFaces`/`startFace`).
struct Patch {
    std::string mName;
    std::string mType;  ///< the `type` entry; empty when the file omitted it
    std::int64_t mNFaces = 0;
    std::int64_t mStartFace = 0;
};

/**
 * @brief Read the value of key @p pKey from a `boundary` sub-dictionary body.
 *
 * `nFaces`/`startFace` are read with `atoll`, but `type` is a word, so it needs
 * real tokenising. The word-boundary guard matters: a bare `find("type")` also
 * matches `physicalType` and `patchType`, both of which are legal entries in the
 * same dictionary and neither of which is the patch's type.
 *
 * @return the value token, or "" when the key is absent.
 */
std::string openfoam_dict_word(const std::string& rBlock, const char* pKey) {
    const std::size_t klen = std::strlen(pKey);
    std::size_t p = 0;
    while ((p = rBlock.find(pKey, p)) != std::string::npos) {
        const bool left_ok = p == 0 || std::isspace(static_cast<unsigned char>(rBlock[p - 1])) ||
                             rBlock[p - 1] == ';';
        const std::size_t after = p + klen;
        const bool right_ok =
            after < rBlock.size() && std::isspace(static_cast<unsigned char>(rBlock[after]));
        if (!left_ok || !right_ok) {
            p = after;
            continue;
        }
        std::size_t a = after;
        while (a < rBlock.size() && std::isspace(static_cast<unsigned char>(rBlock[a])))
            ++a;
        std::size_t b = a;
        while (b < rBlock.size() && !std::isspace(static_cast<unsigned char>(rBlock[b])) &&
               rBlock[b] != ';')
            ++b;
        return rBlock.substr(a, b - a);
    }
    return "";
}

std::vector<Patch> parse_boundary(const std::string& rBody) {
    // Find `name { ... }` blocks with nFaces/startFace.
    std::vector<Patch> patches;
    std::size_t i = 0, n = rBody.size();
    auto skip_ws = [&](std::size_t& p) {
        while (p < n && std::isspace(static_cast<unsigned char>(rBody[p])))
            ++p;
    };
    while (i < n) {
        skip_ws(i);
        // read a token (patch name)
        std::size_t start = i;
        while (i < n && !std::isspace(static_cast<unsigned char>(rBody[i])) && rBody[i] != '{' &&
               rBody[i] != '(' && rBody[i] != ')')
            ++i;
        std::string name = rBody.substr(start, i - start);
        skip_ws(i);
        if (i < n && rBody[i] == '{') {
            // Match the brace by DEPTH, not by the first '}': real patches nest
            // (a `cyclicAMI` carries `transform { ... }`, a `mappedWall` carries
            // `sample { ... }`), and taking the first close truncates the block
            // and then resumes scanning from inside it, inventing patches.
            std::size_t close = std::string::npos;
            int depth = 0;
            for (std::size_t p = i; p < n; ++p) {
                if (rBody[p] == '{') {
                    ++depth;
                } else if (rBody[p] == '}') {
                    if (--depth == 0) {
                        close = p;
                        break;
                    }
                }
            }
            if (close == std::string::npos)
                break;
            std::string block = rBody.substr(i + 1, close - i - 1);
            Patch pt;
            pt.mName = name;
            pt.mType = openfoam_dict_word(block, "type");
            bool has_n = false, has_s = false;
            std::size_t np = block.find("nFaces");
            if (np != std::string::npos) {
                pt.mNFaces = std::atoll(block.c_str() + np + 6);
                has_n = true;
            }
            std::size_t sp = block.find("startFace");
            if (sp != std::string::npos) {
                pt.mStartFace = std::atoll(block.c_str() + sp + 9);
                has_s = true;
            }
            if (has_n && has_s && !name.empty())
                patches.push_back(pt);
            i = close + 1;
        } else if (i < n && (rBody[i] == '(' || rBody[i] == ')')) {
            ++i;  // skip list delimiters
        } else if (name.empty()) {
            ++i;
        }
    }
    return patches;
}

// ---- binary parsers ----

// Return (N, offset just after the outer '(').
std::pair<std::int64_t, std::size_t> data_start(std::string_view rRaw) {
    std::size_t end = rRaw.find('}');
    if (end == std::string::npos)
        throw ReadError("OpenFOAM: no FoamFile header");
    std::size_t lp = rRaw.find('(', end);
    if (lp == std::string::npos)
        throw ReadError("OpenFOAM: no data list '('");
    // last integer between end and lp
    std::int64_t n = 0;
    bool found = false;
    std::size_t i = end;
    while (i < lp) {
        if (std::isdigit(static_cast<unsigned char>(rRaw[i]))) {
            std::int64_t v = 0;
            while (i < lp && std::isdigit(static_cast<unsigned char>(rRaw[i])))
                v = v * 10 + (rRaw[i++] - '0');
            n = v;
            found = true;
        } else {
            ++i;
        }
    }
    if (!found)
        throw ReadError("OpenFOAM: no element count before '('");
    return {n, lp + 1};
}

template <typename T>
T read_le(const char* pP) {
    T v;
    std::memcpy(&v, pP, sizeof(T));
    return v;
}

std::vector<std::array<double, 3>> read_binary_points(std::string_view rRaw, int scalar_bytes) {
    auto [n, start] = data_start(rRaw);
    std::vector<std::array<double, 3>> pts(static_cast<std::size_t>(n));
    const char* base = rRaw.data() + start;
    for (std::int64_t i = 0; i < n; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::size_t off =
                (static_cast<std::size_t>(i) * 3 + j) * static_cast<std::size_t>(scalar_bytes);
            pts[i][j] = scalar_bytes == 4 ? static_cast<double>(read_le<float>(base + off))
                                          : read_le<double>(base + off);
        }
    }
    return pts;
}

std::vector<std::int64_t> read_binary_labels(std::string_view rRaw, int label_bytes) {
    auto [n, start] = data_start(rRaw);
    std::vector<std::int64_t> out(static_cast<std::size_t>(n));
    const char* base = rRaw.data() + start;
    for (std::int64_t i = 0; i < n; ++i) {
        std::size_t off = static_cast<std::size_t>(i) * static_cast<std::size_t>(label_bytes);
        out[i] = label_bytes == 4 ? static_cast<std::int64_t>(read_le<std::int32_t>(base + off))
                                  : read_le<std::int64_t>(base + off);
    }
    return out;
}

std::vector<Face> read_binary_faces(std::string_view rRaw, int label_bytes) {
    auto [nfaces, pos] = data_start(rRaw);
    std::vector<Face> faces(static_cast<std::size_t>(nfaces));
    std::size_t p = pos;
    for (std::int64_t i = 0; i < nfaces; ++i) {
        std::size_t lp = rRaw.find('(', p);
        if (lp == std::string::npos)
            throw ReadError("OpenFOAM: missing '(' in faces");
        std::int64_t count = std::atoll(std::string(rRaw.substr(p, lp - p)).c_str());
        std::size_t blob = lp + 1;
        Face f(static_cast<std::size_t>(count));
        for (std::int64_t j = 0; j < count; ++j) {
            std::size_t off =
                blob + static_cast<std::size_t>(j) * static_cast<std::size_t>(label_bytes);
            f[j] = label_bytes == 4
                       ? static_cast<std::int64_t>(read_le<std::int32_t>(rRaw.data() + off))
                       : read_le<std::int64_t>(rRaw.data() + off);
        }
        faces[i] = std::move(f);
        p = blob + static_cast<std::size_t>(count) * static_cast<std::size_t>(label_bytes) + 1;
    }
    return faces;
}

// ---- dispatch readers ----

std::vector<std::array<double, 3>> read_points(const fs::path& rPath) {
    FoamFormat fmt = detect_format(rPath.string());
    const detail::FileSource source = read_whole(rPath.string());
    const std::string_view raw = source.View();
    if (fmt.mBinary)
        return read_binary_points(raw, fmt.mScalarBytes);
    return parse_points_ascii(strip_comments_and_header(raw));
}

std::vector<Face> read_faces(const fs::path& rPath) {
    FoamFormat fmt = detect_format(rPath.string());
    const detail::FileSource source = read_whole(rPath.string());
    const std::string_view raw = source.View();
    if (fmt.mBinary)
        return read_binary_faces(raw, fmt.mLabelBytes);
    return parse_faces_ascii(strip_comments_and_header(raw));
}

std::vector<std::int64_t> read_int_list(const fs::path& rPath) {
    FoamFormat fmt = detect_format(rPath.string());
    const detail::FileSource source = read_whole(rPath.string());
    const std::string_view raw = source.View();
    if (fmt.mBinary)
        return read_binary_labels(raw, fmt.mLabelBytes);
    return parse_int_list_ascii(strip_comments_and_header(raw));
}

// ---- geometry ----

double triple(const std::array<double, 3>& rA, const std::array<double, 3>& rB,
              const std::array<double, 3>& rC) {
    // a . (b x c)
    double cx = rB[1] * rC[2] - rB[2] * rC[1];
    double cy = rB[2] * rC[0] - rB[0] * rC[2];
    double cz = rB[0] * rC[1] - rB[1] * rC[0];
    return rA[0] * cx + rA[1] * cy + rA[2] * cz;
}

std::array<double, 3> sub(const std::array<double, 3>& rA, const std::array<double, 3>& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}

std::size_t unique_node_count(const std::vector<Face>& rFaces) {
    std::unordered_set<std::int64_t> s;
    for (const auto& f : rFaces)
        for (std::int64_t v : f)
            s.insert(v);
    return s.size();
}

std::unordered_map<std::int64_t, std::unordered_set<std::int64_t>> node_adjacency(
    const std::vector<Face>& rFaces) {
    std::unordered_map<std::int64_t, std::unordered_set<std::int64_t>> adj;
    for (const auto& f : rFaces) {
        std::size_t m = f.size();
        for (std::size_t i = 0; i < m; ++i) {
            std::int64_t a = f[i], b = f[(i + 1) % m];
            adj[a].insert(b);
            adj[b].insert(a);
        }
    }
    return adj;
}

// Returns the ordered top ring, or empty if ambiguous.
std::vector<std::int64_t> match_top(const Face& rBottom, const std::vector<Face>& rOriented) {
    auto adj = node_adjacency(rOriented);
    std::unordered_set<std::int64_t> base(rBottom.begin(), rBottom.end());
    std::vector<std::int64_t> top;
    for (std::int64_t b : rBottom) {
        std::vector<std::int64_t> cand;
        for (std::int64_t x : adj[b])
            if (!base.count(x))
                cand.push_back(x);
        if (cand.size() != 1)
            return {};
        top.push_back(cand[0]);
    }
    return top;
}

using P3 = std::vector<std::array<double, 3>>;

Face build_tetra(const std::vector<Face>& rOriented, const P3& rP) {
    const Face& base = rOriented[0];
    std::unordered_set<std::int64_t> all;
    for (const auto& f : rOriented)
        for (std::int64_t v : f)
            all.insert(v);
    for (std::int64_t v : base)
        all.erase(v);
    std::int64_t apex = *all.begin();
    Face n = {base[0], base[1], base[2], apex};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[2]], rP[n[0]]), sub(rP[n[3]], rP[n[0]])) < 0)
        n = {base[0], base[2], base[1], apex};
    return n;
}

Face build_pyramid(const std::vector<Face>& rOriented, const P3& rP) {
    Face quad;
    for (const auto& f : rOriented)
        if (f.size() == 4) {
            quad = f;
            break;
        }
    std::unordered_set<std::int64_t> all;
    for (const auto& f : rOriented)
        for (std::int64_t v : f)
            all.insert(v);
    for (std::int64_t v : quad)
        all.erase(v);
    std::int64_t apex = *all.begin();
    Face n = {quad[0], quad[1], quad[2], quad[3], apex};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[3]], rP[n[0]]), sub(rP[n[4]], rP[n[0]])) < 0)
        n = {quad[0], quad[3], quad[2], quad[1], apex};
    return n;
}

Face build_wedge(const std::vector<Face>& rOriented, const P3& rP) {
    Face bottom;
    for (const auto& f : rOriented)
        if (f.size() == 3) {
            bottom = f;
            break;
        }
    std::vector<std::int64_t> top = match_top(bottom, rOriented);
    if (top.empty())
        return {};
    Face n = {bottom[0], bottom[1], bottom[2], top[0], top[1], top[2]};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[2]], rP[n[0]]), sub(rP[n[3]], rP[n[0]])) < 0)
        n = {bottom[0], bottom[2], bottom[1], top[0], top[2], top[1]};
    return n;
}

Face build_hexahedron(const std::vector<Face>& rOriented, const P3& rP) {
    Face bottom;
    for (const auto& f : rOriented)
        if (f.size() == 4) {
            bottom = f;
            break;
        }
    std::vector<std::int64_t> top = match_top(bottom, rOriented);
    if (top.empty())
        return {};
    Face n = {bottom[0], bottom[1], bottom[2], bottom[3], top[0], top[1], top[2], top[3]};
    if (triple(sub(rP[n[1]], rP[n[0]]), sub(rP[n[3]], rP[n[0]]), sub(rP[n[4]], rP[n[0]])) < 0)
        n = {bottom[0], bottom[3], bottom[2], bottom[1], top[0], top[3], top[2], top[1]};
    return n;
}

// Classify a cell. Returns {meshio type, connectivity}. For "polyhedron" the
// connectivity is empty (the caller keeps the oriented faces).
std::pair<std::string, Face> reconstruct_cell(const std::vector<Face>& rOriented, const P3& rP) {
    std::size_t nf = rOriented.size();
    std::size_t np = unique_node_count(rOriented);
    if (nf == 4 && np == 4)
        return {"tetra", build_tetra(rOriented, rP)};
    if (nf == 5 && np == 5)
        return {"pyramid", build_pyramid(rOriented, rP)};
    if (nf == 5 && np == 6)
        return {"wedge", build_wedge(rOriented, rP)};
    if (nf == 6 && np == 8)
        return {"hexahedron", build_hexahedron(rOriented, rP)};
    return {"polyhedron", {}};
}

}  // namespace

Mesh read_openfoam(const std::string& rPathIn, OpenFoamInfo& rInfo) {
    // resolve polyMesh directory
    fs::path path(rPathIn);
    fs::path poly;
    if (path.extension() == ".foam") {
        fs::path c = path.parent_path() / "constant" / "polyMesh";
        if (fs::exists(c))
            poly = c;
    }
    if (poly.empty() && path.filename() == "polyMesh" && fs::is_directory(path))
        poly = path;
    if (poly.empty()) {
        for (const fs::path& c : {path / "constant" / "polyMesh", path / "polyMesh"}) {
            if (fs::exists(c)) {
                poly = c;
                break;
            }
        }
    }
    if (poly.empty())
        throw ReadError(detail::format_compat(
            "Could not locate polyMesh from '{}'. Expected <case>/constant/polyMesh/.", rPathIn));
    log::info("Reading polyMesh from {}", poly.string());

    P3 points = read_points(poly / "points");
    std::vector<Face> faces = read_faces(poly / "faces");
    std::vector<std::int64_t> owner = read_int_list(poly / "owner");
    std::vector<std::int64_t> neighbour;
    if (fs::exists(poly / "neighbour"))
        neighbour = read_int_list(poly / "neighbour");
    std::vector<Patch> boundary;
    if (fs::exists(poly / "boundary"))
        boundary = parse_boundary(
            strip_comments_and_header(read_whole((poly / "boundary").string()).View()));

    std::int64_t owner_max = -1, neigh_max = -1;
    for (std::int64_t v : owner)
        owner_max = std::max(owner_max, v);
    for (std::int64_t v : neighbour)
        neigh_max = std::max(neigh_max, v);
    std::int64_t n_cells = owner.empty() ? 0 : std::max(owner_max, neigh_max) + 1;
    log::info("{} points, {} faces, {} cells, {} patches", points.size(), faces.size(), n_cells,
              boundary.size());

    // cell -> face ids
    std::vector<std::vector<std::int64_t>> cell_faces(static_cast<std::size_t>(n_cells));
    for (std::size_t fid = 0; fid < owner.size(); ++fid)
        cell_faces[static_cast<std::size_t>(owner[fid])].push_back(static_cast<std::int64_t>(fid));
    for (std::size_t fid = 0; fid < neighbour.size(); ++fid)
        if (neighbour[fid] >= 0)
            cell_faces[static_cast<std::size_t>(neighbour[fid])].push_back(
                static_cast<std::int64_t>(fid));

    // reconstruct volume cells
    std::vector<std::string> vol_order;
    std::map<std::string, std::vector<Face>> vol_buckets;
    // polyhedra grouped by unique node count -> "polyhedron<N>"
    std::vector<std::string> poly_order;
    std::map<std::string, std::vector<std::vector<Face>>> poly_buckets;

    // Per-cell geometric reconstruction is the expensive part and every cell
    // only reads faces/owner/points -> compute all cells in parallel into a
    // pre-sized result array, then do the (ordered) bucket grouping
    // sequentially.
    struct CellResult {
        std::string mType;         // "" = degenerate (skipped)
        Face mConn;                // named types
        std::vector<Face> mFaces;  // oriented faces, polyhedra only
    };
    std::vector<CellResult> results(static_cast<std::size_t>(n_cells));
    parallel_for(static_cast<std::size_t>(n_cells), [&](std::size_t cs) {
        const std::int64_t cid = static_cast<std::int64_t>(cs);
        std::vector<Face> oriented;
        for (std::int64_t fid : cell_faces[cs]) {
            Face f = faces[static_cast<std::size_t>(fid)];
            if (owner[static_cast<std::size_t>(fid)] != cid)
                std::reverse(f.begin(), f.end());
            oriented.push_back(std::move(f));
        }
        auto [mtype, conn] = reconstruct_cell(oriented, points);
        if (mtype == "polyhedron") {
            results[cs] = {"polyhedron", {}, std::move(oriented)};
        } else if (conn.empty()) {
            results[cs] = {};  // degenerate topology
        } else {
            results[cs] = {std::move(mtype), std::move(conn), {}};
        }
    });

    std::size_t n_skipped = 0;
    std::size_t n_polyhedra = 0;
    for (auto& res : results) {
        if (res.mType == "polyhedron") {
            std::size_t nn = unique_node_count(res.mFaces);
            std::string key = "polyhedron" + std::to_string(nn);
            if (!poly_buckets.count(key))
                poly_order.push_back(key);
            poly_buckets[key].push_back(std::move(res.mFaces));
            ++n_polyhedra;
        } else if (res.mType.empty()) {
            ++n_skipped;
        } else {
            if (!vol_buckets.count(res.mType))
                vol_order.push_back(res.mType);
            vol_buckets[res.mType].push_back(std::move(res.mConn));
        }
    }
    if (n_skipped > 0)
        log::warn("{} cell(s) skipped (degenerate topology).", n_skipped);
    if (n_polyhedra > 0)
        log::info("{} general polyhedron cell(s) found.", n_polyhedra);

    Mesh mesh;
    std::size_t npts = points.size();
    {
        NDArray pts(DType::Float64, {npts, 3});
        double* pdst = pts.As<double>();
        parallel_for(npts, [&](std::size_t i) {
            for (std::size_t j = 0; j < 3; ++j)
                pdst[i * 3 + j] = points[i][j];
        });
        mesh.AssignPoints(std::move(pts));
    }

    std::vector<NDArray> cell_tags;  // one per block, in final block order

    // rectangular volume blocks
    for (const std::string& t : vol_order) {
        const auto& rows = vol_buckets[t];
        std::size_t nc = rows.size();
        std::size_t k = nc ? rows[0].size() : 0;
        NDArray data(DType::Int64, {nc, k});
        std::int64_t* dp = data.As<std::int64_t>();
        parallel_for(nc, [&](std::size_t r) {
            for (std::size_t c = 0; c < k; ++c)
                dp[r * k + c] = rows[r][c];
        });
        mesh.AddCellBlock(t, std::move(data));
        cell_tags.emplace_back(DType::Int64, std::vector<std::size_t>{nc});  // zeros
    }
    // ragged polyhedron blocks
    for (const std::string& key : poly_order) {
        std::vector<std::vector<std::vector<std::int64_t>>> cells;
        for (const auto& cell : poly_buckets[key]) {
            std::vector<std::vector<std::int64_t>> ph;
            for (const auto& face : cell)
                ph.push_back(face);
            cells.push_back(std::move(ph));
        }
        std::size_t nc = cells.size();
        mesh.AddPolyhedronBlock(key, std::move(cells));
        cell_tags.emplace_back(DType::Int64, std::vector<std::size_t>{nc});  // zeros
    }

    // boundary cells grouped by size, with patch family tags
    std::map<int, std::vector<Face>> bysize;  // 3 -> triangles, 4 -> quads
    std::map<int, std::vector<std::int64_t>> tagsize;
    std::vector<Face> poly_faces;
    std::vector<std::int64_t> poly_tags;
    for (std::size_t pidx = 0; pidx < boundary.size(); ++pidx) {
        std::int64_t fam = -(static_cast<std::int64_t>(pidx) + 1);
        rInfo.mCellTags[fam] = {boundary[pidx].mName};
        if (!boundary[pidx].mType.empty())
            rInfo.mPatchTypes[fam] = boundary[pidx].mType;
        for (std::int64_t fid = boundary[pidx].mStartFace;
             fid < boundary[pidx].mStartFace + boundary[pidx].mNFaces; ++fid) {
            if (fid < 0 || static_cast<std::size_t>(fid) >= faces.size())
                continue;
            const Face& f = faces[static_cast<std::size_t>(fid)];
            if (f.size() == 3) {
                bysize[3].push_back(f);
                tagsize[3].push_back(fam);
            } else if (f.size() == 4) {
                bysize[4].push_back(f);
                tagsize[4].push_back(fam);
            } else {
                poly_faces.push_back(f);
                poly_tags.push_back(fam);
            }
        }
    }
    auto add_boundary_block = [&](const std::string& type, const std::vector<Face>& rows,
                                  const std::vector<std::int64_t>& tags) {
        std::size_t nc = rows.size();
        std::size_t k = nc ? rows[0].size() : 0;
        NDArray data(DType::Int64, {nc, k});
        NDArray tag(DType::Int64, {nc});
        std::int64_t* dp = data.As<std::int64_t>();
        std::int64_t* tp = tag.As<std::int64_t>();
        parallel_for(nc, [&](std::size_t r) {
            for (std::size_t c = 0; c < k; ++c)
                dp[r * k + c] = rows[r][c];
            tp[r] = tags[r];
        });
        mesh.AddCellBlock(type, std::move(data));
        cell_tags.push_back(std::move(tag));
    };
    if (!bysize[3].empty())
        add_boundary_block("triangle", bysize[3], tagsize[3]);
    if (!bysize[4].empty())
        add_boundary_block("quad", bysize[4], tagsize[4]);
    if (!poly_faces.empty()) {
        // group boundary polygons by vertex count -> polygon<N>
        std::map<std::size_t, std::vector<Face>> by_n;
        std::map<std::size_t, std::vector<std::int64_t>> tag_n;
        for (std::size_t i = 0; i < poly_faces.size(); ++i) {
            by_n[poly_faces[i].size()].push_back(poly_faces[i]);
            tag_n[poly_faces[i].size()].push_back(poly_tags[i]);
        }
        for (auto& kv : by_n)
            add_boundary_block("polygon" + std::to_string(kv.first), kv.second, tag_n[kv.first]);
    }

    if (!cell_tags.empty())
        mesh.AddCellData("cell_tags", std::move(cell_tags));
    return mesh;
}

// ==========================================================================
//                                  WRITER
// ==========================================================================

namespace {

/// One boundary patch as it will be written.
struct FoamPatchOut {
    std::string mName;
    std::string mType = "patch";
    std::int64_t mNFaces = 0;
    std::int64_t mStartFace = 0;
};

/// The written face order plus the patch table describing its tail.
struct FoamFaceOrder {
    std::vector<std::int64_t> mNewToOld;  ///< written face id -> GlobalFaces id
    std::int64_t mNumInternal = 0;
    std::vector<FoamPatchOut> mPatches;  ///< ascending mStartFace
};

/**
 * @brief Patch types that survive a round trip unchanged.
 *
 * Everything else needs companion dictionary entries `OpenFoamInfo` does not
 * carry -- `cyclic`/`cyclicAMI` need `neighbourPatch`, `processor` needs
 * `myProcNo`/`neighbProcNo`, `mapped*` needs `sample*` -- and OpenFOAM refuses
 * to *load* a case whose patch declares such a type without them. Downgrading
 * to `patch` yields a case that loads and solves with visibly wrong boundary
 * conditions, which is strictly better than one that does not open.
 */
bool foam_type_is_self_contained(const std::string& rType) {
    return rType == "patch" || rType == "wall" || rType == "symmetry" ||
           rType == "symmetryPlane" || rType == "empty" || rType == "wedge";
}

/// Resolve the polyMesh directory. One function for both directions, so the
/// reader's resolution and the writer's cannot drift apart.
fs::path foam_polymesh_dir(const fs::path& rPath, bool ForWrite) {
    if (rPath.extension() == ".foam") {
        const fs::path c = rPath.parent_path() / "constant" / "polyMesh";
        if (ForWrite || fs::exists(c))
            return c;
    }
    if (rPath.filename() == "polyMesh" && (ForWrite || fs::is_directory(rPath)))
        return rPath;
    if (!ForWrite) {
        for (const fs::path& c : {rPath / "constant" / "polyMesh", rPath / "polyMesh"}) {
            if (fs::exists(c))
                return c;
        }
        return {};
    }
    return rPath / "constant" / "polyMesh";
}

/// Standard FoamFile header. `detect_format` reads only `format` and `arch`,
/// but the rest is what makes the file legible to OpenFOAM itself.
void foam_write_header(std::ostream& rOs, const std::string& rClass, const std::string& rObject) {
    // The credit cell is fixed-width (48 chars before the closing box edge) so
    // the banner stays aligned regardless of how long the release string is.
    std::string credit = detail::kProvenanceTag;
    if (credit.size() < 48)
        credit.append(48 - credit.size(), ' ');
    rOs << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
           "| =========                 |                                                 |\n"
           "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
           "|  \\\\    /   O peration     |                                                 |\n"
           "|   \\\\  /    A nd           | "
        << credit
        << "|\n"
           "|    \\\\/     M anipulation  |                                                 |\n"
           "\\*---------------------------------------------------------------------------*/\n"
           "FoamFile\n"
           "{\n"
           "    version     2.0;\n"
           "    format      ascii;\n"
           "    class       "
        << rClass
        << ";\n"
           "    location    \"constant/polyMesh\";\n"
           "    object      "
        << rObject
        << ";\n"
           "}\n"
           "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
}

/**
 * @brief Assign every boundary face to a patch.
 *
 * @return per-`GlobalFaces`-face patch index, `-1` for internal or unassigned.
 */
struct FoamPatchAssignment {
    std::vector<std::int64_t> mFacePatch;
    std::vector<FoamPatchOut> mPatches;
    std::int64_t mNumOrphan = 0;         ///< 2D cell matching no face at all
    std::int64_t mNumInternalTagged = 0;  ///< 2D cell matching an INTERNAL face
};

FoamPatchAssignment foam_assign_patches(const Mesh& rMesh, const detail::GlobalFaces& rFaces,
                                        const OpenFoamInfo& rInfo) {
    FoamPatchAssignment out;
    out.mFacePatch.assign(rFaces.NumFaces(), -1);

    // Family ids in the reader's own order (ascending -fam == ascending patch
    // index), so an OpenFOAM round trip preserves the boundary file's order.
    std::vector<std::int64_t> fams;
    for (const auto& kv : rInfo.mCellTags)
        fams.push_back(kv.first);
    std::sort(fams.begin(), fams.end(), [](std::int64_t a, std::int64_t b) { return a > b; });

    std::unordered_map<std::int64_t, std::size_t> fam_to_patch;
    for (std::int64_t fam : fams) {
        FoamPatchOut p;
        const auto& names = rInfo.mCellTags.at(fam);
        p.mName = names.empty() ? ("patch" + std::to_string(-fam)) : names.front();
        const auto it = rInfo.mPatchTypes.find(fam);
        if (it != rInfo.mPatchTypes.end() && !it->second.empty()) {
            if (foam_type_is_self_contained(it->second)) {
                p.mType = it->second;
            } else {
                log::warn(
                    "OpenFOAM: patch '{}' has type '{}', which needs dictionary entries meshio++ "
                    "does not carry; writing 'patch' instead so the case still loads",
                    p.mName, it->second);
            }
        }
        fam_to_patch[fam] = out.mPatches.size();
        out.mPatches.push_back(std::move(p));
    }

    const detail::FaceLookup lookup(rFaces);
    const bool have_tags = rMesh.HasCellData("cell_tags");

    for (std::size_t block : rFaces.mNonCellBlocks) {
        const auto cb = rMesh.Cells(block);
        if (cell_type_dimension(cell_type_from_name(cb.Type())) != 2)
            continue;
        const NDArray* tags =
            have_tags && rMesh.CellDataNumBlocks("cell_tags") == rMesh.NumCellBlocks()
                ? &rMesh.CellData("cell_tags", block)
                : nullptr;

        std::vector<std::int64_t> ids;
        for (std::size_t i = 0; i < cb.NumCells(); ++i) {
            ids.clear();
            if (cb.IsRagged()) {
                const std::size_t n = cb.RowSize(i);
                const std::int64_t* row = cb.Row(i);
                ids.assign(row, row + n);
            } else {
                const std::size_t npc = cb.NodesPerCell();
                for (std::size_t k = 0; k < npc; ++k)
                    ids.push_back(detail::read_int(cb.Conn(), i * npc + k));
            }
            const std::int64_t fid = lookup.Find(ids.data(), ids.size());
            if (fid < 0) {
                ++out.mNumOrphan;
                continue;
            }
            if (rFaces.mNeighbour[static_cast<std::size_t>(fid)] >= 0) {
                // A face between two cells cannot also be a patch member --
                // OpenFOAM would refuse the case. Legitimate input (an interior
                // baffle), so drop it rather than throw.
                ++out.mNumInternalTagged;
                continue;
            }
            if (out.mFacePatch[static_cast<std::size_t>(fid)] >= 0)
                continue;  // first claim wins
            std::int64_t fam = 0;
            if (tags && i < tags->Shape()[0])
                fam = detail::read_int(*tags, i);
            const auto it = fam_to_patch.find(fam);
            if (fam < 0 && it != fam_to_patch.end())
                out.mFacePatch[static_cast<std::size_t>(fid)] =
                    static_cast<std::int64_t>(it->second);
        }
    }

    // Everything still unassigned joins `defaultFaces` -- blockMesh's own name
    // for exactly this, so the result is a valid single-patch case rather than
    // an error or an invented decomposition.
    std::int64_t n_unassigned = 0;
    for (std::size_t f = 0; f < rFaces.NumFaces(); ++f)
        if (rFaces.mNeighbour[f] < 0 && out.mFacePatch[f] < 0)
            ++n_unassigned;
    if (n_unassigned > 0) {
        const std::int64_t idx = static_cast<std::int64_t>(out.mPatches.size());
        FoamPatchOut p;
        p.mName = "defaultFaces";
        p.mType = "patch";
        out.mPatches.push_back(std::move(p));
        for (std::size_t f = 0; f < rFaces.NumFaces(); ++f)
            if (rFaces.mNeighbour[f] < 0 && out.mFacePatch[f] < 0)
                out.mFacePatch[f] = idx;
    }

    // An empty patch is legal but checkMesh flags it, and its emptiness always
    // means something went wrong upstream.
    std::vector<std::int64_t> counts(out.mPatches.size(), 0);
    for (std::size_t f = 0; f < rFaces.NumFaces(); ++f)
        if (out.mFacePatch[f] >= 0)
            ++counts[static_cast<std::size_t>(out.mFacePatch[f])];
    std::vector<FoamPatchOut> kept;
    std::vector<std::int64_t> remap(out.mPatches.size(), -1);
    for (std::size_t p = 0; p < out.mPatches.size(); ++p) {
        if (counts[p] == 0) {
            log::warn("OpenFOAM: patch '{}' has no faces and is not written",
                      out.mPatches[p].mName);
            continue;
        }
        remap[p] = static_cast<std::int64_t>(kept.size());
        kept.push_back(out.mPatches[p]);
    }
    for (std::size_t f = 0; f < rFaces.NumFaces(); ++f)
        if (out.mFacePatch[f] >= 0)
            out.mFacePatch[f] = remap[static_cast<std::size_t>(out.mFacePatch[f])];
    out.mPatches = std::move(kept);
    return out;
}

/// Order the faces the way OpenFOAM requires.
FoamFaceOrder foam_order_faces(const detail::GlobalFaces& rFaces,
                               const FoamPatchAssignment& rAssign) {
    FoamFaceOrder order;
    order.mPatches = rAssign.mPatches;

    std::vector<std::int64_t> internal, boundary;
    for (std::size_t f = 0; f < rFaces.NumFaces(); ++f)
        (rFaces.mNeighbour[f] >= 0 ? internal : boundary).push_back(static_cast<std::int64_t>(f));

    // Upper-triangular order. The face id is a THIRD key, not decoration: two
    // cells can share two distinct faces, so (owner, neighbour) alone is not a
    // strict weak ordering and std::sort would be undefined behaviour on it.
    std::sort(internal.begin(), internal.end(), [&](std::int64_t a, std::int64_t b) {
        const std::size_t ia = static_cast<std::size_t>(a), ib = static_cast<std::size_t>(b);
        if (rFaces.mOwner[ia] != rFaces.mOwner[ib])
            return rFaces.mOwner[ia] < rFaces.mOwner[ib];
        if (rFaces.mNeighbour[ia] != rFaces.mNeighbour[ib])
            return rFaces.mNeighbour[ia] < rFaces.mNeighbour[ib];
        return a < b;
    });
    // Boundary faces contiguous per patch, patches in table order.
    std::stable_sort(boundary.begin(), boundary.end(), [&](std::int64_t a, std::int64_t b) {
        return rAssign.mFacePatch[static_cast<std::size_t>(a)] <
               rAssign.mFacePatch[static_cast<std::size_t>(b)];
    });

    order.mNumInternal = static_cast<std::int64_t>(internal.size());
    order.mNewToOld = std::move(internal);
    order.mNewToOld.insert(order.mNewToOld.end(), boundary.begin(), boundary.end());

    std::int64_t start = order.mNumInternal;
    for (std::size_t p = 0; p < order.mPatches.size(); ++p) {
        std::int64_t n = 0;
        for (std::int64_t f : boundary)
            if (rAssign.mFacePatch[static_cast<std::size_t>(f)] == static_cast<std::int64_t>(p))
                ++n;
        order.mPatches[p].mStartFace = start;
        order.mPatches[p].mNFaces = n;
        start += n;
    }
    return order;
}

/**
 * @brief Check every clause of the polyMesh ordering contract.
 *
 * Runs in release builds too, deliberately: release is exactly where someone
 * writes a ten-million-cell case, the cost is a handful of flops per face
 * against ASCII formatting that costs far more, and a failure means we were
 * about to hand a solver a corrupt mesh.
 *
 * @return "" when valid, else the first violated clause, named.
 */
std::string foam_validate_order(const detail::GlobalFaces& rFaces, const FoamFaceOrder& rOrder,
                                const FoamPatchAssignment& rAssign, const NDArray& rPoints,
                                std::size_t PointDim, std::size_t NumPoints) {
    const std::size_t nf = rOrder.mNewToOld.size();
    if (nf != rFaces.NumFaces())
        return "C0: the written face list does not cover every face";

    // Cell centroids, for the two normal-direction clauses.
    std::vector<detail::Vec3> centroid(rFaces.NumCells(), detail::Vec3{0, 0, 0});
    for (std::size_t c = 0; c < rFaces.NumCells(); ++c) {
        detail::Vec3 acc{0, 0, 0};
        std::size_t n = 0;
        for (std::size_t k = 0; k < rFaces.NumCellFaces(c); ++k) {
            const std::size_t f =
                static_cast<std::size_t>(std::abs(rFaces.CellFaces(c)[k]) - 1);
            for (std::size_t j = 0; j < rFaces.FaceSize(f); ++j) {
                const detail::Vec3 p =
                    detail::read_point(rPoints, PointDim, rFaces.Face(f)[j]);
                acc[0] += p[0];
                acc[1] += p[1];
                acc[2] += p[2];
                ++n;
            }
        }
        if (n)
            for (int a = 0; a < 3; ++a)
                acc[a] /= static_cast<double>(n);
        centroid[c] = acc;
    }

    std::vector<detail::Vec3> ring;
    for (std::size_t i = 0; i < nf; ++i) {
        const std::size_t f = static_cast<std::size_t>(rOrder.mNewToOld[i]);
        const bool is_internal = i < static_cast<std::size_t>(rOrder.mNumInternal);

        // C2: internal faces first.
        if (is_internal != (rFaces.mNeighbour[f] >= 0))
            return detail::format_compat(
                "C2: face {} is {} but sits in the {} range", i,
                rFaces.mNeighbour[f] >= 0 ? "internal" : "boundary",
                is_internal ? "internal" : "boundary");

        // C7: node ids in range, ring big enough to bound an area.
        if (rFaces.FaceSize(f) < 3)
            return detail::format_compat("C7: face {} has fewer than three nodes", i);
        for (std::size_t k = 0; k < rFaces.FaceSize(f); ++k) {
            const std::int64_t id = rFaces.Face(f)[k];
            if (id < 0 || static_cast<std::size_t>(id) >= NumPoints)
                return detail::format_compat("C7: face {} references node {}, out of range", i,
                                             id);
        }

        ring.clear();
        for (std::size_t k = 0; k < rFaces.FaceSize(f); ++k)
            ring.push_back(detail::read_point(rPoints, PointDim, rFaces.Face(f)[k]));
        const detail::Vec3 nrm = detail::polygon_area_vector(ring.data(), ring.size());
        detail::Vec3 fc{0, 0, 0};
        for (const detail::Vec3& p : ring)
            for (int a = 0; a < 3; ++a)
                fc[a] += p[a] / static_cast<double>(ring.size());

        if (is_internal) {
            // C1: owner < neighbour.
            if (!(rFaces.mOwner[f] < rFaces.mNeighbour[f]))
                return detail::format_compat("C1: face {} has owner {} >= neighbour {}", i,
                                             rFaces.mOwner[f], rFaces.mNeighbour[f]);
            // C3: strictly increasing (owner, neighbour).
            if (i > 0) {
                const std::size_t g = static_cast<std::size_t>(rOrder.mNewToOld[i - 1]);
                const bool ok = rFaces.mOwner[g] < rFaces.mOwner[f] ||
                                (rFaces.mOwner[g] == rFaces.mOwner[f] &&
                                 rFaces.mNeighbour[g] <= rFaces.mNeighbour[f]);
                if (!ok)
                    return detail::format_compat(
                        "C3: face {} has (owner,neighbour)=({},{}) after ({},{})", i,
                        rFaces.mOwner[f], rFaces.mNeighbour[f], rFaces.mOwner[g],
                        rFaces.mNeighbour[g]);
            }
            // C4: normal points owner -> neighbour.
            const detail::Vec3& co = centroid[static_cast<std::size_t>(rFaces.mOwner[f])];
            const detail::Vec3& cn = centroid[static_cast<std::size_t>(rFaces.mNeighbour[f])];
            const double d = nrm[0] * (cn[0] - co[0]) + nrm[1] * (cn[1] - co[1]) +
                             nrm[2] * (cn[2] - co[2]);
            if (!(d > 0.0))
                return detail::format_compat(
                    "C4: internal face {} does not point from owner to neighbour", i);
        } else {
            // C5: boundary normal points out of the domain.
            const detail::Vec3& co = centroid[static_cast<std::size_t>(rFaces.mOwner[f])];
            const double d = nrm[0] * (fc[0] - co[0]) + nrm[1] * (fc[1] - co[1]) +
                             nrm[2] * (fc[2] - co[2]);
            if (!(d > 0.0))
                return detail::format_compat("C5: boundary face {} is wound inward", i);
        }
    }

    // C6: patches partition the boundary range exactly -- AND every face in a
    // patch's range really belongs to that patch. Checking only that the
    // start/count table tiles the range is not enough: the table is built by
    // counting, so it describes a contiguity the face order may simply not
    // have, and the resulting case is silently wrong.
    std::int64_t expect = rOrder.mNumInternal;
    for (std::size_t p = 0; p < rOrder.mPatches.size(); ++p) {
        const FoamPatchOut& patch = rOrder.mPatches[p];
        if (patch.mStartFace != expect)
            return detail::format_compat("C6: patch '{}' starts at {}, expected {}", patch.mName,
                                         patch.mStartFace, expect);
        for (std::int64_t i = patch.mStartFace; i < patch.mStartFace + patch.mNFaces; ++i) {
            const std::size_t f = static_cast<std::size_t>(rOrder.mNewToOld[
                static_cast<std::size_t>(i)]);
            if (rAssign.mFacePatch[f] != static_cast<std::int64_t>(p))
                return detail::format_compat(
                    "C6: face {} sits in patch '{}'s range but belongs to patch {}", i,
                    patch.mName, rAssign.mFacePatch[f]);
        }
        expect += patch.mNFaces;
    }
    if (expect != static_cast<std::int64_t>(nf))
        return detail::format_compat("C6: patches cover {} faces, expected {}",
                                     expect - rOrder.mNumInternal,
                                     static_cast<std::int64_t>(nf) - rOrder.mNumInternal);
    return "";
}

std::ofstream foam_open(const fs::path& rPath) {
    std::ofstream f(rPath, std::ios::binary);
    if (!f)
        throw WriteError("OpenFOAM: could not open for writing: " + rPath.string());
    return f;
}

}  // namespace

void write_openfoam(const std::string& rPath, const Mesh& rMesh, const OpenFoamInfo& rInfo) {
    const detail::GlobalFaces faces = detail::build_global_faces(rMesh);

    if (faces.NumCells() == 0)
        throw WriteError(
            "OpenFOAM: the mesh has no volume cells; a polyMesh needs at least one "
            "tetra/pyramid/wedge/hexahedron/polyhedron cell");
    // A 3D block we could not turn into faces would be a silently dropped solid.
    for (std::size_t block : faces.mNonCellBlocks) {
        const auto cb = rMesh.Cells(block);
        const std::string type(cb.Type());
        if (cell_type_dimension(cell_type_from_name(type)) == 3)
            throw WriteError(detail::format_compat(
                "OpenFOAM: cell type '{}' is 3D but has no face topology in meshio++, so writing "
                "it would silently drop those cells",
                type));
    }
    if (faces.mNumNonManifold > 0)
        throw WriteError(detail::format_compat(
            "OpenFOAM: {} face(s) are shared by three or more cells; a polyMesh face has at most "
            "an owner and one neighbour",
            faces.mNumNonManifold));
    if (faces.mNumUnorientable > 0)
        log::warn("OpenFOAM: {} cell(s) are not closed orientable surfaces; their faces are "
                  "written with the winding they arrived with",
                  faces.mNumUnorientable);
    if (faces.mNumFlipped > 0)
        log::info("OpenFOAM: rewound {} inverted cell(s) so their faces point outward",
                  faces.mNumFlipped);

    const FoamPatchAssignment assign = foam_assign_patches(rMesh, faces, rInfo);
    if (assign.mNumOrphan > 0)
        log::warn("OpenFOAM: {} boundary cell(s) match no cell face and are not written",
                  assign.mNumOrphan);
    if (assign.mNumInternalTagged > 0)
        log::warn("OpenFOAM: {} boundary cell(s) coincide with an INTERNAL face; OpenFOAM cannot "
                  "put such a face on a patch, so they are not written",
                  assign.mNumInternalTagged);

    const FoamFaceOrder order = foam_order_faces(faces, assign);
    const std::string bad = foam_validate_order(faces, order, assign, rMesh.Points(),
                                                rMesh.PointDim(), rMesh.NumPoints());
    if (!bad.empty())
        throw WriteError("OpenFOAM: internal error, the written face order violates " + bad);

    const fs::path poly = foam_polymesh_dir(fs::path(rPath), /*ForWrite=*/true);
    std::error_code ec;
    fs::create_directories(poly, ec);
    if (ec && !fs::is_directory(poly))
        throw WriteError("OpenFOAM: could not create directory " + poly.string() + ": " +
                         ec.message());
    if (fs::path(rPath).extension() == ".foam") {
        // The marker file is what makes the case openable by ParaView and by
        // this reader's own `.foam` branch.
        std::ofstream marker(rPath, std::ios::binary);
    }

    // Companion files this writer does not produce but OpenFOAM would read.
    // Leaving a stale one behind corrupts the case, so remove exactly these --
    // never the whole directory, which may hold a user's own files.
    for (const char* name : {"cellZones", "faceZones", "pointZones", "meshModifiers",
                             "boundaryProcAddressing", "cellProcAddressing",
                             "faceProcAddressing", "pointProcAddressing", "cellLevel",
                             "pointLevel", "level0Edge", "refinementHistory", "surfaceIndex"}) {
        std::error_code rc;
        if (fs::remove(poly / name, rc))
            log::info("OpenFOAM: removed stale {}", name);
    }

    const NDArray& pts = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t np = rMesh.NumPoints();

    {
        std::ofstream f = foam_open(poly / "points");
        foam_write_header(f, "vectorField", "points");
        // The count MUST be on a line of its own: every ASCII parser here takes
        // "the first line that is entirely digits" as the count, so `8(` would
        // be read as data and the list would come back EMPTY, not as an error.
        f << np << "\n(\n";
        f << std::setprecision(16);
        for (std::size_t i = 0; i < np; ++i) {
            const detail::Vec3 p = detail::read_point(pts, dim, static_cast<std::int64_t>(i));
            f << "(" << p[0] << " " << p[1] << " " << p[2] << ")\n";
        }
        f << ")\n";
    }
    {
        std::ofstream f = foam_open(poly / "faces");
        foam_write_header(f, "faceList", "faces");
        f << order.mNewToOld.size() << "\n(\n";
        for (std::int64_t old : order.mNewToOld) {
            const std::size_t fi = static_cast<std::size_t>(old);
            f << faces.FaceSize(fi) << "(";
            for (std::size_t k = 0; k < faces.FaceSize(fi); ++k)
                f << (k ? " " : "") << faces.Face(fi)[k];
            f << ")\n";
        }
        f << ")\n";
    }
    {
        std::ofstream f = foam_open(poly / "owner");
        foam_write_header(f, "labelList", "owner");
        f << order.mNewToOld.size() << "\n(\n";
        for (std::int64_t old : order.mNewToOld)
            f << faces.mOwner[static_cast<std::size_t>(old)] << "\n";
        f << ")\n";
    }
    {
        // Always written, even with zero entries: a stale `neighbour` left from
        // a previous, larger case is one of the nastiest ways to corrupt one.
        // OpenFOAM's `neighbour` holds ONLY internal faces -- our own reader
        // also accepts a -1-padded full-length list, which is exactly why a
        // round trip through it is a weak oracle for this writer.
        std::ofstream f = foam_open(poly / "neighbour");
        foam_write_header(f, "labelList", "neighbour");
        f << order.mNumInternal << "\n(\n";
        for (std::int64_t i = 0; i < order.mNumInternal; ++i)
            f << faces.mNeighbour[static_cast<std::size_t>(order.mNewToOld[
                  static_cast<std::size_t>(i)])]
              << "\n";
        f << ")\n";
    }
    {
        std::ofstream f = foam_open(poly / "boundary");
        foam_write_header(f, "polyBoundaryMesh", "boundary");
        f << order.mPatches.size() << "\n(\n";
        for (const FoamPatchOut& p : order.mPatches) {
            f << "    " << p.mName << "\n    {\n";
            f << "        type            " << p.mType << ";\n";
            f << "        nFaces          " << p.mNFaces << ";\n";
            f << "        startFace       " << p.mStartFace << ";\n";
            f << "    }\n";
        }
        f << ")\n";
    }

    log::info("Wrote polyMesh to {} ({} cells, {} faces, {} internal, {} patches)",
              poly.string(), faces.NumCells(), order.mNewToOld.size(), order.mNumInternal,
              order.mPatches.size());
}

}  // namespace meshioplusplus
