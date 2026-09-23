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
#include <bit>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
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
#include "meshioplusplus/detail/cell_faces.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "face_cells_common.hpp"

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
    auto f = detail::make_classic_ifstream(rPath, std::ios::binary);
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
            // OpenFOAM's own arch strings are "LSB;label=32;scalar=64" (the
            // "BSB" spelling is what a big-endian host's files would carry).
            // Binary bytes this reader decodes are always little-endian, so a
            // file naming anything else is refused by name rather than
            // silently misread.
            if (s.find("BSB") != std::string::npos)
                throw ReadError(
                    "OpenFOAM: big-endian ('BSB') binary files are not supported, only "
                    "little-endian ('LSB')");
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
    auto ss = detail::make_classic_istringstream(out);
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
    auto ss = detail::make_classic_istringstream(rBody);
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
            auto ns = detail::make_classic_istringstream(t);
            double a, b, c;
            if (ns >> a >> b >> c)
                pts.push_back({a, b, c});
        }
    }
    return pts;
}

std::vector<Face> parse_faces_ascii(const std::string& rBody) {
    std::vector<Face> faces;
    auto ss = detail::make_classic_istringstream(rBody);
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
            auto ns = detail::make_classic_istringstream(inside);
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
    auto ss = detail::make_classic_istringstream(rBody);
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
            auto ns = detail::make_classic_istringstream(s);
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

/**
 * @brief Scan top-level `name { ... }` blocks, matching braces by DEPTH.
 *
 * Shared by `parse_boundary` and the zone-file parsers below: both formats
 * are a flat `N ( name { ... } name { ... } ... )` list of named
 * sub-dictionaries. Depth matching (not "the first `}`") matters here too --
 * a patch's `transform`/`sample` sub-block would otherwise truncate it.
 *
 * @return `(name, block body)` pairs, in file order.
 */
std::vector<std::pair<std::string, std::string>> foam_named_blocks(const std::string& rBody) {
    std::vector<std::pair<std::string, std::string>> blocks;
    std::size_t i = 0, n = rBody.size();
    auto skip_ws = [&](std::size_t& p) {
        while (p < n && std::isspace(static_cast<unsigned char>(rBody[p])))
            ++p;
    };
    while (i < n) {
        skip_ws(i);
        std::size_t start = i;
        while (i < n && !std::isspace(static_cast<unsigned char>(rBody[i])) && rBody[i] != '{' &&
               rBody[i] != '(' && rBody[i] != ')')
            ++i;
        std::string name = rBody.substr(start, i - start);
        skip_ws(i);
        if (i < n && rBody[i] == '{') {
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
            if (!name.empty())
                blocks.emplace_back(name, rBody.substr(i + 1, close - i - 1));
            i = close + 1;
        } else if (i < n && (rBody[i] == '(' || rBody[i] == ')')) {
            ++i;  // skip list delimiters
        } else if (name.empty()) {
            ++i;
        }
    }
    return blocks;
}

std::vector<Patch> parse_boundary(const std::string& rBody) {
    std::vector<Patch> patches;
    for (const auto& [name, block] : foam_named_blocks(rBody)) {
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
        if (has_n && has_s)
            patches.push_back(pt);
    }
    return patches;
}

/// One named `cellZone`/`faceZone`/`pointZone` entry: a name plus its member
/// ids (cell/face/point ids, in the file's own numbering).
struct Zone {
    std::string mName;
    std::vector<std::int64_t> mIds;
};

/**
 * @brief Read the `List<label>` value of key @p pKey out of a zone block.
 *
 * Zone blocks look like `type cellZone; cellLabels List<label> 3(0 5 9);` --
 * `flipMap` (a `faceZone`-only `List<bool>`) is deliberately never read: a
 * flip only matters for a zone consumer that walks faces directionally
 * (cyclic AMI construction, e.g.), and `Region`'s `Side` entries carry no
 * orientation bit to hold it in. See doc/formats/openfoam.md.
 */
std::vector<std::int64_t> foam_zone_label_list(const std::string& rBlock, const char* pKey) {
    std::size_t p = rBlock.find(pKey);
    if (p == std::string::npos)
        return {};
    std::size_t lp = rBlock.find('(', p);
    if (lp == std::string::npos)
        return {};
    std::size_t rp = lp + 1;
    int depth = 1;
    while (rp < rBlock.size() && depth > 0) {
        if (rBlock[rp] == '(')
            ++depth;
        else if (rBlock[rp] == ')')
            --depth;
        ++rp;
    }
    const std::string inside = rBlock.substr(lp + 1, rp - lp - 2);
    auto ss = detail::make_classic_istringstream(inside);
    std::vector<std::int64_t> out;
    std::int64_t v;
    while (ss >> v)
        out.push_back(v);
    return out;
}

/// Parse an **ASCII** `cellZones`/`faceZones`/`pointZones` file body (the
/// comment-and-header-stripped text). `parse_zone_file_binary`, further
/// down, is the binary counterpart -- it cannot reuse this one, since
/// `strip_comments_and_header`'s comment-removal pass is unsafe to run over
/// a binary body.
std::vector<Zone> parse_zone_file(const std::string& rBody, const char* pLabelKey) {
    std::vector<Zone> zones;
    for (const auto& [name, block] : foam_named_blocks(rBody))
        zones.push_back({name, foam_zone_label_list(block, pLabelKey)});
    return zones;
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

/**
 * @brief Binary counterpart of `parse_zone_file`.
 *
 * A zone file's *structure* (zone count, names, `{`/`}`, the `type ...;`
 * line, the `List<label>` keyword and its own decimal count) stays plain
 * text even under `format binary;` -- only each zone's id payload is raw
 * bytes. That payload can legitimately contain a byte equal to `{`, `}` or
 * `/` (a small id's low byte routinely does, e.g. `123` as little-endian
 * `int32` starts with `0x7B` == `{`), so this walks the text directly on
 * the **raw, unstripped** file bytes -- never through `strip_comments_and_header`,
 * whose comment-removal pass scans the whole body and would desync on
 * exactly those bytes -- and explicitly skips `count * LabelBytes` bytes as
 * one opaque unit wherever it recognizes @p pLabelKey, instead of ever
 * scanning byte-by-byte across a blob for a delimiter.
 */
std::vector<Zone> parse_zone_file_binary(std::string_view rRaw, const char* pLabelKey,
                                         int LabelBytes) {
    std::vector<Zone> zones;
    auto [nzones, pos] = data_start(rRaw);
    std::size_t p = pos;
    for (std::int64_t z = 0; z < nzones; ++z) {
        while (p < rRaw.size() && std::isspace(static_cast<unsigned char>(rRaw[p])))
            ++p;
        const std::size_t name_start = p;
        while (p < rRaw.size() && !std::isspace(static_cast<unsigned char>(rRaw[p])))
            ++p;
        const std::string name(rRaw.substr(name_start, p - name_start));

        const std::size_t brace = rRaw.find('{', p);
        if (brace == std::string_view::npos)
            throw ReadError("OpenFOAM: zone '" + name + "' has no opening '{'");
        const std::size_t key_pos = rRaw.find(pLabelKey, brace + 1);
        if (key_pos == std::string_view::npos)
            throw ReadError("OpenFOAM: zone '" + name + "' has no '" + pLabelKey + "'");
        const std::size_t lparen = rRaw.find('(', key_pos);
        if (lparen == std::string_view::npos)
            throw ReadError("OpenFOAM: zone '" + name + "' has no '(' after '" + pLabelKey + "'");

        std::int64_t count = 0;
        bool found = false;
        for (std::size_t i = key_pos; i < lparen; ++i)
            if (std::isdigit(static_cast<unsigned char>(rRaw[i]))) {
                std::int64_t v = 0;
                while (i < lparen && std::isdigit(static_cast<unsigned char>(rRaw[i])))
                    v = v * 10 + (rRaw[i++] - '0');
                count = v;
                found = true;
            }
        if (!found)
            throw ReadError("OpenFOAM: zone '" + name + "' has no count before '('");

        std::vector<std::int64_t> ids(static_cast<std::size_t>(count));
        const char* base = rRaw.data() + lparen + 1;
        for (std::int64_t i = 0; i < count; ++i) {
            const std::size_t off = static_cast<std::size_t>(i) * static_cast<std::size_t>(LabelBytes);
            ids[static_cast<std::size_t>(i)] =
                LabelBytes == 4 ? static_cast<std::int64_t>(read_le<std::int32_t>(base + off))
                                : read_le<std::int64_t>(base + off);
        }
        zones.push_back({name, std::move(ids)});

        // Resume the text scan only *after* the raw blob -- everything from
        // here to this zone's own closing '}' (its `);` terminator, the
        // closing brace) is plain text again.
        const std::size_t blob_end =
            lparen + 1 + static_cast<std::size_t>(count) * static_cast<std::size_t>(LabelBytes);
        const std::size_t zone_close = rRaw.find('}', blob_end);
        if (zone_close == std::string_view::npos)
            throw ReadError("OpenFOAM: zone '" + name + "' has no closing '}'");
        p = zone_close + 1;
    }
    return zones;
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

// The face-to-cell kernel (triple, match_top, build_*, reconstruct_cell) is
// shared with the Fluent reader: formats/face_cells_common.hpp.
using face_cells::P3;
using face_cells::reconstruct_cell;
using face_cells::unique_node_count;

// ---- decomposed (processorN) cases (v11.4.0, roadmap §1 tier B2) ----

/// A single `polyMesh`'s raw, un-reconstructed file contents -- what
/// `read_openfoam` used to read inline before it could also come from
/// `reconstruct_decomposed`.
struct RawPolyMesh {
    P3 mPoints;
    std::vector<Face> mFaces;
    std::vector<std::int64_t> mOwner;
    std::vector<std::int64_t> mNeighbour;  ///< only the internal faces
    std::vector<Patch> mBoundary;
};

RawPolyMesh read_raw_polymesh(const fs::path& rPoly) {
    RawPolyMesh raw;
    raw.mPoints = read_points(rPoly / "points");
    raw.mFaces = read_faces(rPoly / "faces");
    raw.mOwner = read_int_list(rPoly / "owner");
    if (fs::exists(rPoly / "neighbour"))
        raw.mNeighbour = read_int_list(rPoly / "neighbour");
    if (fs::exists(rPoly / "boundary"))
        raw.mBoundary = parse_boundary(
            strip_comments_and_header(read_whole((rPoly / "boundary").string()).View()));
    return raw;
}

/// One processor's `polyMesh` plus the addressing lists that map its local
/// ids back onto the undecomposed case's global ones.
struct ProcMesh {
    RawPolyMesh mRaw;
    std::vector<std::int64_t> mPointAddr;     ///< local point -> global point id
    std::vector<std::int64_t> mCellAddr;      ///< local cell -> global cell id
    std::vector<std::int64_t> mFaceAddr;      ///< local face -> signed (global face id + 1)
    std::vector<std::int64_t> mBoundaryAddr;  ///< local patch -> global patch id, -1 if none
};

/// The processor-local face (and its orientation) claiming a given global
/// face id -- one entry for a face interior to one processor or an original
/// external boundary face, two for a face split by decomposition.
struct DecompFaceClaim {
    std::int64_t mProc = -1;
    std::int64_t mLocalFace = -1;
    bool mFlipped = false;
};

/// The local patch containing local face @p LocalFace, or `npos`.
std::size_t foam_patch_of_local_face(const std::vector<Patch>& rBoundary, std::int64_t LocalFace) {
    for (std::size_t p = 0; p < rBoundary.size(); ++p)
        if (LocalFace >= rBoundary[p].mStartFace &&
            LocalFace < rBoundary[p].mStartFace + rBoundary[p].mNFaces)
            return p;
    return static_cast<std::size_t>(-1);
}

/// The sorted processor indices of a decomposed case's `<root>/processorN/`
/// directories that carry a `constant/polyMesh` -- not assumed contiguous.
std::vector<std::size_t> foam_processor_ids(const fs::path& rCaseRoot) {
    std::vector<std::size_t> ids;
    std::error_code ec;
    for (const auto& entry : fs::directory_iterator(rCaseRoot, ec)) {
        const std::string name = entry.path().filename().string();
        if (name.rfind("processor", 0) != 0)
            continue;
        const std::string digits = name.substr(std::strlen("processor"));
        if (digits.empty() || digits.find_first_not_of("0123456789") != std::string::npos)
            continue;
        if (fs::exists(entry.path() / "constant" / "polyMesh"))
            ids.push_back(static_cast<std::size_t>(std::atoll(digits.c_str())));
    }
    std::sort(ids.begin(), ids.end());
    return ids;
}

/**
 * @brief Reassemble a decomposed case's `processorN` directories into one
 * global `RawPolyMesh`, mirroring what `reconstructParMesh` does on disk.
 *
 * Points/cells/faces are placed at the global ids their processor's own
 * `*ProcAddressing` files (plain `labelList`s, read like `owner`) name.
 * A global face claimed by exactly one processor is either an original
 * external boundary face or one interior to that processor alone; claimed by
 * two, it is an internal face `decomposePar` split at a processor boundary --
 * the positive-signed `faceProcAddressing` entry names the true owner side,
 * the negative-signed one the neighbour side (the sign meaning "this local
 * copy is stored reversed relative to the global orientation"). A boundary
 * face's global patch comes from its owning processor's `boundaryProcAddressing`
 * (missing/negative marks a `processor*` inter-rank patch, dropped -- it has no
 * counterpart in the original case).
 */
RawPolyMesh reconstruct_decomposed(const fs::path& rCaseRoot,
                                   const std::vector<std::size_t>& rProcIds) {
    std::vector<ProcMesh> procs(rProcIds.size());
    std::int64_t max_point = -1, max_cell = -1, max_face = -1;
    for (std::size_t k = 0; k < rProcIds.size(); ++k) {
        const fs::path poly =
            rCaseRoot / ("processor" + std::to_string(rProcIds[k])) / "constant" / "polyMesh";
        ProcMesh& pm = procs[k];
        pm.mRaw = read_raw_polymesh(poly);
        pm.mPointAddr = read_int_list(poly / "pointProcAddressing");
        pm.mCellAddr = read_int_list(poly / "cellProcAddressing");
        pm.mFaceAddr = read_int_list(poly / "faceProcAddressing");
        if (fs::exists(poly / "boundaryProcAddressing"))
            pm.mBoundaryAddr = read_int_list(poly / "boundaryProcAddressing");
        for (std::int64_t v : pm.mPointAddr)
            max_point = std::max(max_point, v);
        for (std::int64_t v : pm.mCellAddr)
            max_cell = std::max(max_cell, v);
        for (std::int64_t v : pm.mFaceAddr)
            max_face = std::max(max_face, std::abs(v) - 1);
    }
    const std::size_t n_points = static_cast<std::size_t>(max_point + 1);
    const std::size_t n_faces = static_cast<std::size_t>(max_face + 1);

    P3 points(n_points);
    for (const ProcMesh& pm : procs)
        for (std::size_t i = 0; i < pm.mPointAddr.size(); ++i)
            points[static_cast<std::size_t>(pm.mPointAddr[i])] = pm.mRaw.mPoints[i];

    // Every processor-local face that claims a given global id.
    std::vector<std::array<DecompFaceClaim, 2>> claims(n_faces);
    std::vector<std::uint8_t> n_claims(n_faces, 0);
    for (std::size_t k = 0; k < procs.size(); ++k) {
        const auto& fa = procs[k].mFaceAddr;
        for (std::size_t i = 0; i < fa.size(); ++i) {
            const std::size_t g = static_cast<std::size_t>(std::abs(fa[i]) - 1);
            const std::uint8_t slot = n_claims[g]++;
            if (slot < 2)
                claims[g][slot] = {static_cast<std::int64_t>(k), static_cast<std::int64_t>(i),
                                   fa[i] < 0};
        }
    }

    std::vector<Face> internal_faces;
    std::vector<std::int64_t> internal_owner, internal_neighbour;
    // Global patch id -> its (name, type) plus the member faces' node rings
    // and owner cells, in the order they are found.
    std::map<std::int64_t, Patch> patch_table;
    std::map<std::int64_t, std::vector<Face>> patch_faces;
    std::map<std::int64_t, std::vector<std::int64_t>> patch_owners;
    std::size_t n_dropped_processor_faces = 0, n_dropped_unclaimed = 0;

    for (std::size_t g = 0; g < n_faces; ++g) {
        if (n_claims[g] == 0) {
            ++n_dropped_unclaimed;  // an id `faceProcAddressing` never actually used
            continue;
        }
        const DecompFaceClaim* owner_claim = nullptr;
        const DecompFaceClaim* neigh_claim = nullptr;
        for (std::uint8_t k = 0; k < std::min<std::uint8_t>(n_claims[g], 2); ++k) {
            const DecompFaceClaim& c = claims[g][k];
            (c.mFlipped ? neigh_claim : owner_claim) = &c;
        }
        if (!owner_claim)
            owner_claim = &claims[g][0];  // defensive: both flipped should not happen

        const ProcMesh& op = procs[static_cast<std::size_t>(owner_claim->mProc)];
        const Face& lf = op.mRaw.mFaces[static_cast<std::size_t>(owner_claim->mLocalFace)];
        Face gf(lf.size());
        for (std::size_t k = 0; k < lf.size(); ++k)
            gf[k] = op.mPointAddr[static_cast<std::size_t>(lf[k])];
        const std::int64_t owner_cell = op.mCellAddr[static_cast<std::size_t>(
            op.mRaw.mOwner[static_cast<std::size_t>(owner_claim->mLocalFace)])];

        if (neigh_claim) {
            const ProcMesh& np = procs[static_cast<std::size_t>(neigh_claim->mProc)];
            const std::int64_t neigh_cell = np.mCellAddr[static_cast<std::size_t>(
                np.mRaw.mOwner[static_cast<std::size_t>(neigh_claim->mLocalFace)])];
            internal_faces.push_back(std::move(gf));
            internal_owner.push_back(owner_cell);
            internal_neighbour.push_back(neigh_cell);
            continue;
        }

        // A genuine boundary face: resolve its global patch via the owning
        // processor's own local patch + `boundaryProcAddressing`.
        const std::size_t local_patch =
            foam_patch_of_local_face(op.mRaw.mBoundary, owner_claim->mLocalFace);
        std::int64_t global_patch = -1;
        if (local_patch != static_cast<std::size_t>(-1) &&
            local_patch < op.mBoundaryAddr.size())
            global_patch = op.mBoundaryAddr[local_patch];
        const bool is_processor_patch =
            global_patch < 0 || (local_patch != static_cast<std::size_t>(-1) &&
                                 op.mRaw.mBoundary[local_patch].mType.rfind("processor", 0) == 0);
        if (is_processor_patch) {
            ++n_dropped_processor_faces;
            continue;
        }
        if (!patch_table.count(global_patch)) {
            Patch p = op.mRaw.mBoundary[local_patch];  // name/type only; counts recomputed below
            p.mNFaces = 0;
            p.mStartFace = 0;
            patch_table[global_patch] = p;
        }
        patch_faces[global_patch].push_back(std::move(gf));
        patch_owners[global_patch].push_back(owner_cell);
    }
    if (n_dropped_processor_faces > 0)
        log::info("OpenFOAM: reconstructed case drops {} inter-processor patch face(s)",
                  n_dropped_processor_faces);
    if (n_dropped_unclaimed > 0)
        log::warn("OpenFOAM: {} face id(s) in *ProcAddressing were never claimed",
                  n_dropped_unclaimed);

    RawPolyMesh out;
    out.mPoints = std::move(points);
    out.mFaces = std::move(internal_faces);
    out.mOwner = internal_owner;
    out.mNeighbour = std::move(internal_neighbour);
    for (auto& kv : patch_table) {
        kv.second.mStartFace = static_cast<std::int64_t>(out.mFaces.size());
        for (Face& f : patch_faces[kv.first])
            out.mFaces.push_back(std::move(f));
        for (std::int64_t c : patch_owners[kv.first])
            out.mOwner.push_back(c);
        kv.second.mNFaces = static_cast<std::int64_t>(out.mFaces.size()) - kv.second.mStartFace;
        out.mBoundary.push_back(kv.second);
    }
    return out;
}

// ---- time-directory fields (v11.4.0, roadmap §1 tier B2) ----

/// One field's shape: number of components, and whether it belongs on
/// points (`pointScalarField`/`pointVectorField`) rather than cells.
struct FoamFieldClass {
    int mComponents = 0;
    bool mIsPoint = false;
    bool mSupported = true;
};

FoamFieldClass foam_field_class(const std::string& rClass) {
    if (rClass == "volScalarField")
        return {1, false, true};
    if (rClass == "volVectorField")
        return {3, false, true};
    if (rClass == "volSymmTensorField")
        return {6, false, true};
    if (rClass == "volTensorField")
        return {9, false, true};
    if (rClass == "pointScalarField")
        return {1, true, true};
    if (rClass == "pointVectorField")
        return {3, true, true};
    return {0, false, false};  // e.g. surfaceScalarField: no cell/point home
}

/// A field value list: `mFlat` holds `mComponents` entries (uniform) or
/// `mCount * mComponents` (nonuniform), read by `foam_read_internal_field`.
struct FoamField {
    bool mUniform = false;
    std::int64_t mCount = 0;
    std::vector<double> mFlat;
};

/// Parse `uniform <value>` (`rText` starting AT the `uniform` keyword).
std::vector<double> foam_scan_uniform_value(std::string_view rText, int components) {
    std::size_t p = std::strlen("uniform");
    while (p < rText.size() && std::isspace(static_cast<unsigned char>(rText[p])))
        ++p;
    std::vector<double> out;
    if (components == 1) {
        out.push_back(detail::parse_double(std::string(rText.substr(p))));
        return out;
    }
    const std::size_t lp = rText.find('(', p);
    const std::size_t rp = rText.find(')', lp);
    if (lp == std::string::npos || rp == std::string::npos)
        return out;
    auto ss = detail::make_classic_istringstream(std::string(rText.substr(lp + 1, rp - lp - 1)));
    double v;
    while (ss >> v)
        out.push_back(v);
    return out;
}

/// Parse `nonuniform List<T>\n<N>\n(\n<entries>\n)` (`rText` starting AT the
/// `nonuniform` keyword) -- ASCII only; the binary variant is read directly
/// via `data_start` in `foam_read_internal_field`, which needs the whole raw
/// buffer rather than a text view.
FoamField foam_scan_nonuniform_list(std::string_view rText, int components) {
    FoamField out;
    const std::string text_owned(rText);
    auto ss = detail::make_classic_istringstream(text_owned);
    std::string line;
    bool have_n = false;
    std::int64_t n = 0;
    while (std::getline(ss, line)) {
        std::string s = openfoam_strip(line);
        if (s.empty())
            continue;
        if (!have_n) {
            if (s.find_first_not_of("0123456789") == std::string::npos) {
                n = std::atoll(s.c_str());
                have_n = true;
            }
            continue;
        }
        if (s == "(")
            break;
    }
    out.mCount = n;
    out.mFlat.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(components));
    for (std::int64_t i = 0; i < n && std::getline(ss, line);) {
        std::string s = openfoam_strip(line);
        if (s.empty())
            continue;
        if (components == 1) {
            out.mFlat.push_back(detail::parse_double(s));
        } else {
            for (char& c : s)
                if (c == '(' || c == ')')
                    c = ' ';
            auto ls = detail::make_classic_istringstream(s);
            double v;
            while (ls >> v)
                out.mFlat.push_back(v);
        }
        ++i;
    }
    return out;
}

/**
 * @brief Read a field file's `internalField`.
 *
 * `uniform` is always plain text, in both ASCII and binary field files (a
 * single small value is never worth binary-encoding), so it is scanned the
 * same way regardless of `FoamFormat`. `nonuniform` follows `points`/`faces`/
 * `owner`'s own dispatch: ASCII is line-scanned, binary reuses `data_start` --
 * safe here because nothing between the FoamFile header and `internalField`'s
 * own data list can introduce a stray `(` (`dimensions` uses `[...]`).
 */
FoamField foam_read_internal_field(const fs::path& rPath, int components) {
    const FoamFormat fmt = detect_format(rPath.string());
    const detail::FileSource source = read_whole(rPath.string());
    const std::string_view raw = source.View();
    const std::size_t kp = raw.find("internalField");
    if (kp == std::string::npos)
        throw ReadError("OpenFOAM: field file has no internalField: " + rPath.string());
    std::size_t p = kp + std::strlen("internalField");
    while (p < raw.size() && std::isspace(static_cast<unsigned char>(raw[p])))
        ++p;
    if (raw.compare(p, 7, "uniform") == 0) {
        FoamField out;
        out.mUniform = true;
        out.mFlat = foam_scan_uniform_value(raw.substr(p), components);
        return out;
    }
    if (raw.compare(p, 10, "nonuniform") != 0)
        throw ReadError("OpenFOAM: internalField is neither uniform nor nonuniform: " +
                        rPath.string());
    if (!fmt.mBinary)
        return foam_scan_nonuniform_list(raw.substr(p), components);

    auto [n, start] = data_start(raw);
    FoamField out;
    out.mCount = n;
    out.mFlat.resize(static_cast<std::size_t>(n) * static_cast<std::size_t>(components));
    const char* base = raw.data() + start;
    for (std::int64_t i = 0; i < n; ++i)
        for (int c = 0; c < components; ++c) {
            const std::size_t off = (static_cast<std::size_t>(i) * static_cast<std::size_t>(components) +
                                     static_cast<std::size_t>(c)) *
                                    static_cast<std::size_t>(fmt.mScalarBytes);
            out.mFlat[static_cast<std::size_t>(i) * static_cast<std::size_t>(components) +
                     static_cast<std::size_t>(c)] =
                fmt.mScalarBytes == 4 ? static_cast<double>(read_le<float>(base + off))
                                     : read_le<double>(base + off);
        }
    return out;
}

/// Parse a time directory's name as an OpenFOAM time value, requiring the
/// WHOLE name to be consumed (so `"0.1_backup"` is correctly not a time dir).
bool foam_parse_time_dir_name(const std::string& rName, double& rValue) {
    if (rName.empty())
        return false;
    const char* end = nullptr;
    const double v = detail::parse_double(rName.c_str(), end);
    if (end != rName.c_str() + rName.size())
        return false;
    rValue = v;
    return true;
}

/// One `<case_root>/<numeric>/` time directory: its parsed value and its
/// own (exact, on-disk) name -- kept together because a value re-formatted
/// back to text ("0.1" vs "0.100000") is not reliably the same string.
struct FoamTimeDir {
    double mValue = 0.0;
    std::string mName;
};

/// Time directories holding at least one regular file (a field), sorted
/// ascending by value. `0` is included only when it holds fields -- the
/// historical (field-free) `read_openfoam` never depended on a `0/`
/// directory existing at all.
std::vector<FoamTimeDir> foam_time_dirs(const fs::path& rCaseRoot) {
    std::vector<FoamTimeDir> dirs;
    std::error_code ec;
    for (const auto& entry : fs::directory_iterator(rCaseRoot, ec)) {
        if (!entry.is_directory())
            continue;
        const std::string name = entry.path().filename().string();
        double t = 0.0;
        if (!foam_parse_time_dir_name(name, t))
            continue;
        std::error_code ec2;
        bool has_field = false;
        for (const auto& f : fs::directory_iterator(entry.path(), ec2)) {
            if (f.is_regular_file()) {
                has_field = true;
                break;
            }
        }
        if (has_field)
            dirs.push_back({t, name});
    }
    std::sort(dirs.begin(), dirs.end(),
             [](const FoamTimeDir& a, const FoamTimeDir& b) { return a.mValue < b.mValue; });
    return dirs;
}

/// Field file names directly inside a time directory (regular files only,
/// non-recursive -- `uniform/`, `polyMesh/` and other sub-directories a
/// moving-mesh case may carry there are not field files).
std::vector<std::string> foam_field_files(const fs::path& rTimeDir) {
    std::vector<std::string> names;
    std::error_code ec;
    for (const auto& entry : fs::directory_iterator(rTimeDir, ec))
        if (entry.is_regular_file())
            names.push_back(entry.path().filename().string());
    std::sort(names.begin(), names.end());
    return names;
}

/// The `class` entry of a field file's `FoamFile` header (`volScalarField`,
/// …), read the same cheap line-scan way `detect_format` reads `format`/
/// `arch` -- no full parse needed just to classify the field.
std::string foam_field_file_class(const fs::path& rPath) {
    auto f = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!f)
        return {};
    std::string line;
    while (std::getline(f, line)) {
        const std::string s = openfoam_strip(line);
        if (s.rfind("class", 0) == 0) {
            std::string rest = openfoam_strip(s.substr(std::strlen("class")));
            if (!rest.empty() && rest.back() == ';')
                rest.pop_back();
            return openfoam_strip(rest);
        }
        if (s == "}")
            break;
    }
    return {};
}

}  // namespace

Mesh read_openfoam(const std::string& rPathIn, OpenFoamInfo& rInfo) {
    return read_openfoam(rPathIn, ReadOptions{}, rInfo);
}

Mesh read_openfoam(const std::string& rPathIn, const ReadOptions& rOptions, OpenFoamInfo& rInfo) {
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
    // Multi-region case (v11.4.0, roadmap §1 tier B2): no single
    // `constant/polyMesh`, but `<case>/constant/<region>/polyMesh` per
    // region. `case_root` is the case directory regardless of which of the
    // three `rPathIn` forms was given.
    const fs::path case_root = path.extension() == ".foam" ? path.parent_path() : path;
    if (poly.empty() && !rInfo.mRegion.empty()) {
        const fs::path c = case_root / "constant" / rInfo.mRegion / "polyMesh";
        if (fs::exists(c))
            poly = c;
        else
            throw ReadError(detail::format_compat(
                "OpenFOAM: region '{}' has no {}", rInfo.mRegion, c.string()));
    }
    if (poly.empty() && fs::exists(case_root / "constant" / "regionProperties")) {
        std::vector<std::string> regions;
        std::error_code ec;
        for (const auto& entry : fs::directory_iterator(case_root / "constant", ec)) {
            if (entry.is_directory() && fs::exists(entry.path() / "polyMesh"))
                regions.push_back(entry.path().filename().string());
        }
        std::sort(regions.begin(), regions.end());
        std::string joined;
        for (std::size_t i = 0; i < regions.size(); ++i)
            joined += (i ? ", " : "") + regions[i];
        throw ReadError(detail::format_compat(
            "'{}' is a multi-region case; set OpenFoamInfo::mRegion to one of: {}", rPathIn,
            joined));
    }

    // Decomposed case (v11.4.0, roadmap §1 tier B2): no single
    // `constant/polyMesh` and no `<region>` selected, but `processorN`
    // directories each carrying their own `constant/polyMesh`.
    bool decomposed = false;
    std::vector<std::size_t> proc_ids;
    if (poly.empty()) {
        proc_ids = foam_processor_ids(case_root);
        decomposed = !proc_ids.empty();
        if (decomposed)
            poly = case_root;  // for the "Reading polyMesh from" log line below only
    }
    if (poly.empty())
        throw ReadError(detail::format_compat(
            "Could not locate polyMesh from '{}'. Expected <case>/constant/polyMesh/.", rPathIn));
    log::info("Reading polyMesh from {}", poly.string());

    RawPolyMesh raw;
    if (decomposed) {
        log::info("OpenFOAM: reconstructing {} from {} processor director{}", case_root.string(),
                  proc_ids.size(), proc_ids.size() == 1 ? "y" : "ies");
        raw = reconstruct_decomposed(case_root, proc_ids);
    } else {
        raw = read_raw_polymesh(poly);
    }
    P3 points = std::move(raw.mPoints);
    std::vector<Face> faces = std::move(raw.mFaces);
    std::vector<std::int64_t> owner = std::move(raw.mOwner);
    std::vector<std::int64_t> neighbour = std::move(raw.mNeighbour);
    std::vector<Patch> boundary = std::move(raw.mBoundary);

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

    // Original OpenFOAM cell id -> (is_poly, bucket key, row within that
    // bucket), captured as cells are bucketed so zone regions can later
    // recover each cell's position in the final Mesh (see `orig_cell_to_global`
    // below). `row == npos` marks a skipped (degenerate) cell.
    constexpr std::size_t npos = static_cast<std::size_t>(-1);
    std::vector<std::tuple<bool, std::string, std::size_t>> placement(
        static_cast<std::size_t>(n_cells), std::tuple<bool, std::string, std::size_t>{false, "",
                                                                                       npos});

    std::size_t n_skipped = 0;
    std::size_t n_polyhedra = 0;
    for (std::size_t cid = 0; cid < results.size(); ++cid) {
        auto& res = results[cid];
        if (res.mType == "polyhedron") {
            std::size_t nn = unique_node_count(res.mFaces);
            std::string key = "polyhedron" + std::to_string(nn);
            if (!poly_buckets.count(key))
                poly_order.push_back(key);
            placement[cid] = {true, key, poly_buckets[key].size()};
            poly_buckets[key].push_back(std::move(res.mFaces));
            ++n_polyhedra;
        } else if (res.mType.empty()) {
            ++n_skipped;
        } else {
            if (!vol_buckets.count(res.mType))
                vol_order.push_back(res.mType);
            placement[cid] = {false, res.mType, vol_buckets[res.mType].size()};
            vol_buckets[res.mType].push_back(std::move(res.mConn));
        }
    }
    if (n_skipped > 0)
        log::warn("{} cell(s) skipped (degenerate topology).", n_skipped);
    if (n_polyhedra > 0)
        log::info("{} general polyhedron cell(s) found.", n_polyhedra);

    // Block index of each bucket in the order blocks are about to be added
    // (volume types first, in `vol_order`, then polyhedron buckets in
    // `poly_order`) -- the boundary 2D blocks added further down come after
    // both, so this table stays valid for `detail::block_bases` once those
    // volume/polyhedron blocks are on the mesh.
    std::unordered_map<std::string, std::size_t> vol_block_index, poly_block_index;
    for (std::size_t k = 0; k < vol_order.size(); ++k)
        vol_block_index[vol_order[k]] = k;
    for (std::size_t k = 0; k < poly_order.size(); ++k)
        poly_block_index[poly_order[k]] = vol_order.size() + k;

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

    // Original OpenFOAM cell id -> global (block-major) cell index, used by
    // both the zones-as-regions block below and time-directory field
    // attachment further down. `detail::block_bases(mesh)` only needs to be
    // right for the blocks added so far (volume + polyhedron); boundary (2D)
    // blocks are always appended after, so their prefix sums never change
    // what is computed here.
    std::vector<std::int64_t> orig_cell_to_global(static_cast<std::size_t>(n_cells), -1);
    {
        const std::vector<std::int64_t> bases = detail::block_bases(mesh);
        for (std::size_t cid = 0; cid < placement.size(); ++cid) {
            const auto& [is_poly, key, row] = placement[cid];
            if (row == npos)
                continue;
            const std::size_t block = is_poly ? poly_block_index.at(key) : vol_block_index.at(key);
            orig_cell_to_global[cid] =
                detail::block_row_to_global(bases, block, static_cast<std::int64_t>(row));
        }
    }

    // ---- zones as named regions (cellZones/faceZones/pointZones) ----
    {
        // Local facet index of face `fid` within its owner cell `cid`,
        // matching the order `AddPolyhedronBlock`/`AddCellBlock` end up
        // storing: for a polyhedron, that is simply `fid`'s position in
        // `cell_faces[cid]` (the same order the per-cell reconstruction loop
        // above walked to build `oriented`, and that `AddPolyhedronBlock`
        // preserves verbatim); a named type's reconstruction instead rewinds
        // into `cell_faces.hpp`'s canonical per-type order, so the position
        // is recovered by matching each canonical face's corner *set*
        // against `fid`'s -- orientation- and start-point-independent.
        auto local_facet = [&](std::size_t cid, std::int64_t fid) -> int {
            const auto& [is_poly, key, row] = placement[cid];
            if (row == npos)
                return -1;
            if (is_poly) {
                const auto& cf = cell_faces[cid];
                for (std::size_t k = 0; k < cf.size(); ++k)
                    if (cf[k] == fid)
                        return static_cast<int>(k);
                return -1;
            }
            const Face& conn = vol_buckets.at(key)[row];
            const auto& facedefs = detail::cell_faces(cell_type_from_name(key));
            const Face& fnodes = faces[static_cast<std::size_t>(fid)];
            const std::unordered_set<std::int64_t> target(fnodes.begin(), fnodes.end());
            for (std::size_t k = 0; k < facedefs.size(); ++k) {
                std::unordered_set<std::int64_t> cand;
                for (int c = 0; c < facedefs[k].mNumCorners; ++c)
                    cand.insert(conn[facedefs[k].mNodes[c]]);
                if (cand == target)
                    return static_cast<int>(k);
            }
            return -1;
        };

        auto read_zone_file = [&](const char* pFile, const char* pKey) {
            std::vector<Zone> zones;
            const fs::path zone_path = poly / pFile;
            if (!fs::exists(zone_path))
                return zones;
            const FoamFormat zone_fmt = detect_format(zone_path.string());
            if (zone_fmt.mBinary) {
                // strip_comments_and_header's comment-removal pass scans the
                // whole body and is unsafe over raw id bytes -- see
                // parse_zone_file_binary's own doc comment.
                const detail::FileSource source = read_whole(zone_path.string());
                zones = parse_zone_file_binary(source.View(), pKey, zone_fmt.mLabelBytes);
            } else {
                zones = parse_zone_file(
                    strip_comments_and_header(read_whole(zone_path.string()).View()), pKey);
            }
            return zones;
        };

        for (const Zone& z : read_zone_file("cellZones", "cellLabels")) {
            std::size_t n_dropped = 0;
            std::vector<std::int64_t> entries;
            entries.reserve(z.mIds.size());
            for (std::int64_t cid : z.mIds) {
                if (cid < 0 || cid >= n_cells ||
                    orig_cell_to_global[static_cast<std::size_t>(cid)] < 0) {
                    ++n_dropped;
                    continue;
                }
                entries.push_back(orig_cell_to_global[static_cast<std::size_t>(cid)]);
            }
            if (n_dropped > 0)
                log::warn("OpenFOAM: cellZone '{}' drops {} entr{} (degenerate cell)", z.mName,
                          n_dropped, n_dropped == 1 ? "y" : "ies");
            NDArray arr = NDArray::Uninit(DType::Int64, {entries.size()});
            std::copy(entries.begin(), entries.end(), arr.As<std::int64_t>());
            mesh.AddRegion(Region(z.mName, RegionKind::Cell, std::move(arr)));
        }

        for (const Zone& z : read_zone_file("pointZones", "pointLabels")) {
            std::vector<std::int64_t> entries;
            entries.reserve(z.mIds.size());
            for (std::int64_t pid : z.mIds)
                if (pid >= 0 && static_cast<std::size_t>(pid) < npts)
                    entries.push_back(pid);
            NDArray arr = NDArray::Uninit(DType::Int64, {entries.size()});
            std::copy(entries.begin(), entries.end(), arr.As<std::int64_t>());
            mesh.AddRegion(Region(z.mName, RegionKind::Point, std::move(arr)));
        }

        for (const Zone& z : read_zone_file("faceZones", "faceLabels")) {
            std::size_t n_dropped = 0;
            std::vector<std::int64_t> pairs;
            pairs.reserve(z.mIds.size() * 2);
            for (std::int64_t fid : z.mIds) {
                if (fid < 0 || static_cast<std::size_t>(fid) >= faces.size()) {
                    ++n_dropped;
                    continue;
                }
                const std::int64_t cid = owner[static_cast<std::size_t>(fid)];
                const std::int64_t global = orig_cell_to_global[static_cast<std::size_t>(cid)];
                const int facet = local_facet(static_cast<std::size_t>(cid), fid);
                if (global < 0 || facet < 0) {
                    ++n_dropped;
                    continue;
                }
                pairs.push_back(global);
                pairs.push_back(facet);
            }
            if (n_dropped > 0)
                log::warn("OpenFOAM: faceZone '{}' drops {} entr{} (degenerate owner cell)",
                          z.mName, n_dropped, n_dropped == 1 ? "y" : "ies");
            NDArray arr = NDArray::Uninit(DType::Int64, {pairs.size() / 2, 2});
            std::copy(pairs.begin(), pairs.end(), arr.As<std::int64_t>());
            mesh.AddRegion(Region(z.mName, RegionKind::Side, std::move(arr)));
        }
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

    // ---- time-directory fields (v11.4.0, roadmap §1 tier B2) ----
    if (rOptions.WantsAnyData()) {
        const std::vector<FoamTimeDir> time_dirs = foam_time_dirs(case_root);
        if (!time_dirs.empty()) {
            const std::size_t step = rOptions.ResolveTimeStep(time_dirs.size());
            const fs::path time_dir = case_root / time_dirs[step].mName;
            for (const std::string& field_name : foam_field_files(time_dir)) {
                if (!rOptions.WantsArray(field_name))
                    continue;
                const fs::path field_path = time_dir / field_name;
                const std::string cls = foam_field_file_class(field_path);
                const FoamFieldClass fc = foam_field_class(cls);
                if (!fc.mSupported) {
                    log::warn(
                        "OpenFOAM: field '{}' has class '{}', which has no point/cell home; skipped",
                        field_name, cls.empty() ? "?" : cls);
                    continue;
                }
                const FoamField values = foam_read_internal_field(field_path, fc.mComponents);

                if (fc.mIsPoint) {
                    if (!values.mUniform && static_cast<std::size_t>(values.mCount) != npts) {
                        log::warn(
                            "OpenFOAM: field '{}' has {} value(s), expected {} point(s); skipped",
                            field_name, values.mCount, npts);
                        continue;
                    }
                    NDArray arr = NDArray::Uninit(
                        DType::Float64, {npts, static_cast<std::size_t>(fc.mComponents)});
                    double* dst = arr.As<double>();
                    if (values.mUniform) {
                        for (std::size_t i = 0; i < npts; ++i)
                            for (int c = 0; c < fc.mComponents; ++c)
                                dst[i * static_cast<std::size_t>(fc.mComponents) +
                                    static_cast<std::size_t>(c)] = values.mFlat[static_cast<std::size_t>(c)];
                    } else {
                        std::copy(values.mFlat.begin(), values.mFlat.end(), dst);
                    }
                    mesh.AddPointData(field_name, std::move(arr));
                    continue;
                }

                if (!values.mUniform && values.mCount != n_cells) {
                    log::warn("OpenFOAM: field '{}' has {} value(s), expected {} cell(s); skipped",
                              field_name, values.mCount, n_cells);
                    continue;
                }
                std::vector<NDArray> blocks(mesh.NumCellBlocks());
                for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b) {
                    NDArray arr = NDArray::Uninit(
                        DType::Float64,
                        {mesh.Cells(b).NumCells(), static_cast<std::size_t>(fc.mComponents)});
                    double* dst = arr.As<double>();
                    std::fill(dst, dst + arr.Size(), std::numeric_limits<double>::quiet_NaN());
                    blocks[b] = std::move(arr);
                }
                for (std::size_t cid = 0; cid < placement.size(); ++cid) {
                    const auto& [is_poly, key, row] = placement[cid];
                    if (row == npos)
                        continue;
                    const std::size_t block =
                        is_poly ? poly_block_index.at(key) : vol_block_index.at(key);
                    double* dst = blocks[block].As<double>();
                    for (int c = 0; c < fc.mComponents; ++c)
                        dst[row * static_cast<std::size_t>(fc.mComponents) + static_cast<std::size_t>(c)] =
                            values.mUniform ? values.mFlat[static_cast<std::size_t>(c)]
                                            : values.mFlat[cid * static_cast<std::size_t>(fc.mComponents) +
                                                            static_cast<std::size_t>(c)];
                }
                mesh.AddCellData(field_name, std::move(blocks));
            }
        }
    }

    return mesh;
}

MeshMetadata read_openfoam_metadata(const std::string& rPathIn, const ReadOptions& rOptions) {
    const fs::path path(rPathIn);
    const fs::path case_root = path.extension() == ".foam" ? path.parent_path() : path;

    std::vector<double> times;
    for (const FoamTimeDir& t : foam_time_dirs(case_root))
        times.push_back(t.mValue);

    OpenFoamInfo info;
    MeshMetadata meta = metadata_from_mesh(read_openfoam(rPathIn, rOptions, info));
    meta.mTimeValues = std::move(times);
    meta.mFellBackToFullRead = true;
    return meta;
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
    return rType == "patch" || rType == "wall" || rType == "symmetry" || rType == "symmetryPlane" ||
           rType == "empty" || rType == "wedge";
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
void foam_write_header(std::ostream& rOs, const std::string& rClass, const std::string& rObject,
                       const OpenFoamWriteOptions& rOpts) {
    // The credit cell is fixed-width (48 chars before the closing box edge) so
    // the banner stays aligned regardless of how long the release string is.
    std::string credit = detail::provenance_lines(detail::SlotTier::Bounded)[0];
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
           "    format      "
        << (rOpts.mBinary ? "binary" : "ascii") << ";\n";
    if (rOpts.mBinary)
        rOs << "    arch        \"LSB;label=" << rOpts.mLabelBits
            << ";scalar=" << rOpts.mScalarBits << "\";\n";
    rOs << "    class       "
        << rClass
        << ";\n"
           "    location    \"constant/polyMesh\";\n"
           "    object      "
        << rObject
        << ";\n"
           "}\n"
           "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
}

// ---- binary body writers ----
// Raw little-endian bytes, matching read_binary_points/read_binary_labels/
// read_binary_faces's own expectations exactly (this file's read side, not
// invented independently) -- see foam_write_header's arch line for the
// widths a reader needs to parse these back.

void foam_write_binary_label(std::ostream& rOs, std::int64_t v, int LabelBits) {
    if (LabelBits == 32) {
        const std::int32_t v32 = static_cast<std::int32_t>(v);
        rOs.write(reinterpret_cast<const char*>(&v32), sizeof(v32));
    } else {
        rOs.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }
}

void foam_write_binary_scalar(std::ostream& rOs, double v, int ScalarBits) {
    if (ScalarBits == 32) {
        const float f = static_cast<float>(v);
        rOs.write(reinterpret_cast<const char*>(&f), sizeof(f));
    } else {
        rOs.write(reinterpret_cast<const char*>(&v), sizeof(v));
    }
}

/**
 * @brief Assign every boundary face to a patch.
 *
 * @return per-`GlobalFaces`-face patch index, `-1` for internal or unassigned.
 */
struct FoamPatchAssignment {
    std::vector<std::int64_t> mFacePatch;
    std::vector<FoamPatchOut> mPatches;
    std::int64_t mNumOrphan = 0;          ///< 2D cell matching no face at all
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
                detail::provenance_note(
                    "data-dropped", "patch '" + p.mName + "' downgraded from type '" + it->second +
                                        "' to 'patch' -- its companion dictionary entries "
                                        "are not carried");
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
            detail::provenance_note("regions-dropped", "patch '" + out.mPatches[p].mName +
                                                           "' dropped -- it has no faces");
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
            const std::size_t f = static_cast<std::size_t>(std::abs(rFaces.CellFaces(c)[k]) - 1);
            for (std::size_t j = 0; j < rFaces.FaceSize(f); ++j) {
                const detail::Vec3 p = detail::read_point(rPoints, PointDim, rFaces.Face(f)[j]);
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
            return detail::format_compat("C2: face {} is {} but sits in the {} range", i,
                                         rFaces.mNeighbour[f] >= 0 ? "internal" : "boundary",
                                         is_internal ? "internal" : "boundary");

        // C7: node ids in range, ring big enough to bound an area.
        if (rFaces.FaceSize(f) < 3)
            return detail::format_compat("C7: face {} has fewer than three nodes", i);
        for (std::size_t k = 0; k < rFaces.FaceSize(f); ++k) {
            const std::int64_t id = rFaces.Face(f)[k];
            if (id < 0 || static_cast<std::size_t>(id) >= NumPoints)
                return detail::format_compat("C7: face {} references node {}, out of range", i, id);
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
            const double d =
                nrm[0] * (cn[0] - co[0]) + nrm[1] * (cn[1] - co[1]) + nrm[2] * (cn[2] - co[2]);
            if (!(d > 0.0))
                return detail::format_compat(
                    "C4: internal face {} does not point from owner to neighbour", i);
        } else {
            // C5: boundary normal points out of the domain.
            const detail::Vec3& co = centroid[static_cast<std::size_t>(rFaces.mOwner[f])];
            const double d =
                nrm[0] * (fc[0] - co[0]) + nrm[1] * (fc[1] - co[1]) + nrm[2] * (fc[2] - co[2]);
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
            const std::size_t f =
                static_cast<std::size_t>(rOrder.mNewToOld[static_cast<std::size_t>(i)]);
            if (rAssign.mFacePatch[f] != static_cast<std::int64_t>(p))
                return detail::format_compat(
                    "C6: face {} sits in patch '{}'s range but belongs to patch {}", i, patch.mName,
                    rAssign.mFacePatch[f]);
        }
        expect += patch.mNFaces;
    }
    if (expect != static_cast<std::int64_t>(nf))
        return detail::format_compat("C6: patches cover {} faces, expected {}",
                                     expect - rOrder.mNumInternal,
                                     static_cast<std::int64_t>(nf) - rOrder.mNumInternal);
    return "";
}

auto foam_open(const fs::path& rPath) {
    auto f = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!f)
        throw WriteError("OpenFOAM: could not open for writing: " + rPath.string());
    return f;
}

/// One zone as it will be written: a name plus its member ids, already
/// converted to the OpenFOAM numbering (compact cell id / point id / written
/// face id -- see the three `foam_collect_*_zones` callers).
using FoamZoneOut = std::pair<std::string, std::vector<std::int64_t>>;

void foam_write_zone_file(const fs::path& rPath, const char* pClass, const char* pObject,
                          const char* pZoneType, const char* pLabelKey,
                          const std::vector<FoamZoneOut>& rZones,
                          const OpenFoamWriteOptions& rOpts) {
    auto f = foam_open(rPath);
    foam_write_header(f, pClass, pObject, rOpts);
    f << rZones.size() << "\n(\n";
    for (const auto& [name, ids] : rZones) {
        f << name << "\n{\n";
        f << "    type " << pZoneType << ";\n";
        f << "    " << pLabelKey << " List<label>\n    " << ids.size() << "\n    (";
        if (rOpts.mBinary) {
            for (std::int64_t id : ids)
                foam_write_binary_label(f, id, rOpts.mLabelBits);
            f << ");\n";
        } else {
            f << "\n";
            for (std::int64_t id : ids)
                f << "    " << id << "\n";
            f << "    );\n";
        }
        f << "}\n";
    }
    f << ")\n";
}

/// `Region`s of kind `Cell` -> `cellZones` entries: a global cell index maps
/// 1:1 onto a written OpenFOAM cell id via `rG2C` (`GlobalFaces::mCellToGlobal`
/// inverted) -- cells are never reordered on write, unlike faces.
std::vector<FoamZoneOut> foam_collect_cell_zones(
    const Mesh& rMesh, const std::unordered_map<std::int64_t, std::int64_t>& rG2C) {
    std::vector<FoamZoneOut> zones;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Cell)
            continue;
        std::vector<std::int64_t> ids;
        const std::int64_t* e = r.Entries();
        std::size_t n_dropped = 0;
        for (std::size_t k = 0; k < r.NumEntries(); ++k) {
            const auto it = rG2C.find(e[k]);
            if (it == rG2C.end()) {
                ++n_dropped;
                continue;
            }
            ids.push_back(it->second);
        }
        if (n_dropped > 0)
            log::warn("OpenFOAM: cellZone '{}' drops {} entr{} outside the volume mesh", r.mName,
                      n_dropped, n_dropped == 1 ? "y" : "ies");
        zones.emplace_back(r.mName, std::move(ids));
    }
    return zones;
}

/// `Region`s of kind `Point` -> `pointZones` entries: a point index needs no
/// conversion, since points are never reordered on write either.
std::vector<FoamZoneOut> foam_collect_point_zones(const Mesh& rMesh) {
    std::vector<FoamZoneOut> zones;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Point)
            continue;
        const std::int64_t* e = r.Entries();
        zones.emplace_back(r.mName, std::vector<std::int64_t>(e, e + r.NumEntries()));
    }
    return zones;
}

/// `Region`s of kind `Side` -> `faceZones` entries. A `(global cell, local
/// facet)` pair becomes a written face id via `rG2C` (global -> compact cell),
/// `GlobalFaces::CellFaces` (compact cell + local facet -> signed GlobalFaces
/// face id) and `rOldToNew` (`FoamFaceOrder::mNewToOld` inverted). `flipMap` is
/// never written -- see `foam_zone_label_list`'s doc comment on the read side.
std::vector<FoamZoneOut> foam_collect_face_zones(
    const Mesh& rMesh, const detail::GlobalFaces& rFaces,
    const std::unordered_map<std::int64_t, std::int64_t>& rG2C,
    const std::vector<std::int64_t>& rOldToNew) {
    std::vector<FoamZoneOut> zones;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Side)
            continue;
        std::vector<std::int64_t> ids;
        const std::int64_t* e = r.Entries();
        std::size_t n_dropped = 0;
        for (std::size_t k = 0; k < r.NumEntries(); ++k) {
            const std::int64_t global_cell = e[2 * k];
            const std::int64_t facet = e[2 * k + 1];
            const auto it = rG2C.find(global_cell);
            if (it == rG2C.end()) {
                ++n_dropped;
                continue;
            }
            const std::size_t compact = static_cast<std::size_t>(it->second);
            if (facet < 0 || static_cast<std::size_t>(facet) >= rFaces.NumCellFaces(compact)) {
                ++n_dropped;
                continue;
            }
            const std::int64_t signed_face = rFaces.CellFaces(compact)[facet];
            const std::size_t old_face = static_cast<std::size_t>(std::abs(signed_face) - 1);
            ids.push_back(rOldToNew[old_face]);
        }
        if (n_dropped > 0)
            log::warn("OpenFOAM: faceZone '{}' drops {} entr{} outside the volume mesh", r.mName,
                      n_dropped, n_dropped == 1 ? "y" : "ies");
        zones.emplace_back(r.mName, std::move(ids));
    }
    return zones;
}

}  // namespace

void write_openfoam(const std::string& rPath, const Mesh& rMesh, const OpenFoamInfo& rInfo,
                    const OpenFoamWriteOptions& rOpts) {
    if (rOpts.mBinary && std::endian::native == std::endian::big)
        throw WriteError(
            "OpenFOAM: binary write requested on a big-endian host; OpenFOAM binary files are "
            "little-endian only, and writing big-endian bytes under a 'LSB' arch header would "
            "silently corrupt every value read back");
    if (rOpts.mBinary && rOpts.mLabelBits != 32 && rOpts.mLabelBits != 64)
        throw WriteError("OpenFOAM: label_bits must be 32 or 64");
    if (rOpts.mBinary && rOpts.mScalarBits != 32 && rOpts.mScalarBits != 64)
        throw WriteError("OpenFOAM: scalar_bits must be 32 or 64");

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
        log::warn(
            "OpenFOAM: {} cell(s) are not closed orientable surfaces; their faces are "
            "written with the winding they arrived with",
            faces.mNumUnorientable);
    if (faces.mNumFlipped > 0)
        log::info("OpenFOAM: rewound {} inverted cell(s) so their faces point outward",
                  faces.mNumFlipped);

    const FoamPatchAssignment assign = foam_assign_patches(rMesh, faces, rInfo);
    if (assign.mNumOrphan > 0) {
        log::warn("OpenFOAM: {} boundary cell(s) match no cell face and are not written",
                  assign.mNumOrphan);
        detail::provenance_note("cells-dropped",
                                std::to_string(assign.mNumOrphan) +
                                    " boundary cell(s) dropped -- they match no cell face");
    }
    if (assign.mNumInternalTagged > 0) {
        log::warn(
            "OpenFOAM: {} boundary cell(s) coincide with an INTERNAL face; OpenFOAM cannot "
            "put such a face on a patch, so they are not written",
            assign.mNumInternalTagged);
        detail::provenance_note("cells-dropped",
                                std::to_string(assign.mNumInternalTagged) +
                                    " boundary cell(s) dropped -- they coincide with an internal "
                                    "face, which OpenFOAM cannot put on a patch");
    }

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
        auto marker = detail::make_classic_ofstream(rPath, std::ios::binary);
    }

    // Companion files this writer does not produce but OpenFOAM would read.
    // Leaving a stale one behind corrupts the case, so remove exactly these --
    // never the whole directory, which may hold a user's own files.
    // `cellZones`/`faceZones`/`pointZones` are handled separately below: this
    // writer produces them once the mesh carries the matching `Region` kind,
    // and only deletes a stale one when it no longer does (so an old zone
    // file is not left behind once its region is removed from the mesh).
    for (const char* name :
         {"meshModifiers", "boundaryProcAddressing", "cellProcAddressing", "faceProcAddressing",
          "pointProcAddressing", "cellLevel", "pointLevel", "level0Edge", "refinementHistory",
          "surfaceIndex"}) {
        std::error_code rc;
        if (fs::remove(poly / name, rc))
            log::info("OpenFOAM: removed stale {}", name);
    }

    const NDArray& pts = rMesh.Points();
    const std::size_t dim = rMesh.PointDim();
    const std::size_t np = rMesh.NumPoints();

    {
        auto f = foam_open(poly / "points");
        foam_write_header(f, "vectorField", "points", rOpts);
        // The count MUST be on a line of its own: every ASCII parser here takes
        // "the first line that is entirely digits" as the count, so `8(` would
        // be read as data and the list would come back EMPTY, not as an error.
        if (rOpts.mBinary) {
            f << np << "\n(";
            for (std::size_t i = 0; i < np; ++i) {
                const detail::Vec3 p = detail::read_point(pts, dim, static_cast<std::int64_t>(i));
                for (double c : p)
                    foam_write_binary_scalar(f, c, rOpts.mScalarBits);
            }
            f << ")\n";
        } else {
            f << np << "\n(\n";
            f << std::setprecision(16);
            for (std::size_t i = 0; i < np; ++i) {
                const detail::Vec3 p = detail::read_point(pts, dim, static_cast<std::int64_t>(i));
                f << "(" << p[0] << " " << p[1] << " " << p[2] << ")\n";
            }
            f << ")\n";
        }
    }
    {
        auto f = foam_open(poly / "faces");
        foam_write_header(f, "faceList", "faces", rOpts);
        // Non-contiguous, so each face is its own length-prefixed labelList,
        // in binary exactly as in ASCII -- the same shape read_binary_faces
        // (this file's own reader) expects, not CompactListList.
        if (rOpts.mBinary) {
            f << order.mNewToOld.size() << "\n(";
            for (std::int64_t old : order.mNewToOld) {
                const std::size_t fi = static_cast<std::size_t>(old);
                f << faces.FaceSize(fi) << "(";
                for (std::size_t k = 0; k < faces.FaceSize(fi); ++k)
                    foam_write_binary_label(f, faces.Face(fi)[k], rOpts.mLabelBits);
                f << ")";
            }
            f << ")\n";
        } else {
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
    }
    {
        auto f = foam_open(poly / "owner");
        foam_write_header(f, "labelList", "owner", rOpts);
        if (rOpts.mBinary) {
            f << order.mNewToOld.size() << "\n(";
            for (std::int64_t old : order.mNewToOld)
                foam_write_binary_label(f, faces.mOwner[static_cast<std::size_t>(old)],
                                        rOpts.mLabelBits);
            f << ")\n";
        } else {
            f << order.mNewToOld.size() << "\n(\n";
            for (std::int64_t old : order.mNewToOld)
                f << faces.mOwner[static_cast<std::size_t>(old)] << "\n";
            f << ")\n";
        }
    }
    {
        // Always written, even with zero entries: a stale `neighbour` left from
        // a previous, larger case is one of the nastiest ways to corrupt one.
        // OpenFOAM's `neighbour` holds ONLY internal faces -- our own reader
        // also accepts a -1-padded full-length list, which is exactly why a
        // round trip through it is a weak oracle for this writer.
        auto f = foam_open(poly / "neighbour");
        foam_write_header(f, "labelList", "neighbour", rOpts);
        if (rOpts.mBinary) {
            f << order.mNumInternal << "\n(";
            for (std::int64_t i = 0; i < order.mNumInternal; ++i)
                foam_write_binary_label(
                    f,
                    faces.mNeighbour[static_cast<std::size_t>(
                        order.mNewToOld[static_cast<std::size_t>(i)])],
                    rOpts.mLabelBits);
            f << ")\n";
        } else {
            f << order.mNumInternal << "\n(\n";
            for (std::int64_t i = 0; i < order.mNumInternal; ++i)
                f << faces.mNeighbour[static_cast<std::size_t>(
                         order.mNewToOld[static_cast<std::size_t>(i)])]
                  << "\n";
            f << ")\n";
        }
    }
    {
        auto f = foam_open(poly / "boundary");
        foam_write_header(f, "polyBoundaryMesh", "boundary", rOpts);
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

    // ---- zones from named regions (cellZones/faceZones/pointZones) ----
    {
        std::unordered_map<std::int64_t, std::int64_t> global_to_compact;
        for (std::size_t c = 0; c < faces.mCellToGlobal.size(); ++c)
            global_to_compact[faces.mCellToGlobal[c]] = static_cast<std::int64_t>(c);
        std::vector<std::int64_t> old_to_new(faces.NumFaces());
        for (std::size_t i = 0; i < order.mNewToOld.size(); ++i)
            old_to_new[static_cast<std::size_t>(order.mNewToOld[i])] = static_cast<std::int64_t>(i);

        auto write_or_remove = [&](const char* pFile, const std::vector<FoamZoneOut>& rZones,
                                   const char* pClass, const char* pZoneType,
                                   const char* pLabelKey) {
            if (rZones.empty()) {
                std::error_code rc;
                if (fs::remove(poly / pFile, rc))
                    log::info("OpenFOAM: removed stale {}", pFile);
                return;
            }
            foam_write_zone_file(poly / pFile, pClass, pFile, pZoneType, pLabelKey, rZones, rOpts);
        };

        write_or_remove("cellZones", foam_collect_cell_zones(rMesh, global_to_compact),
                        "cellZoneList", "cellZone", "cellLabels");
        write_or_remove("pointZones", foam_collect_point_zones(rMesh), "pointZoneList",
                        "pointZone", "pointLabels");
        write_or_remove("faceZones",
                        foam_collect_face_zones(rMesh, faces, global_to_compact, old_to_new),
                        "faceZoneList", "faceZone", "faceLabels");
    }

    log::info("Wrote polyMesh to {} ({} cells, {} faces, {} internal, {} patches)", poly.string(),
              faces.NumCells(), order.mNewToOld.size(), order.mNumInternal, order.mPatches.size());
}

void write_openfoam(const std::string& rPath, const Mesh& rMesh, const OpenFoamInfo& rInfo) {
    write_openfoam(rPath, rMesh, rInfo, OpenFoamWriteOptions{});
}

}  // namespace meshioplusplus
