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
// Conservative content-based format detection: match only unambiguous leading
// signatures, return "" otherwise. Signatures shared by several formats (the
// generic HDF5 magic, a headerless binary STL) are intentionally NOT claimed.
// A directory is sniffed by the files it holds (Elmer, OpenFOAM).

// System includes
#include <algorithm>
#include <bit>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/sniff.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/formats/marc.hpp"
#include "meshioplusplus/formats/lsdyna_d3plot.hpp"
#include "meshioplusplus/formats/z88.hpp"

namespace meshioplusplus {

namespace {

// Whether `hay` starts with `needle`.
bool sniff_starts_with(const std::string& rHay, const char* pNeedle) {
    const std::string needle(pNeedle);
    return rHay.size() >= needle.size() && rHay.compare(0, needle.size(), needle) == 0;
}

// Whether `hay` contains `needle`.
bool sniff_contains(const std::string& rHay, const char* pNeedle) {
    return rHay.find(pNeedle) != std::string::npos;
}

// The prefix with leading ASCII whitespace removed.
std::string sniff_lstrip(const std::string& rIn) {
    std::size_t i = 0;
    while (i < rIn.size() && std::isspace(static_cast<unsigned char>(rIn[i])) != 0)
        ++i;
    return rIn.substr(i);
}

// COMSOL files open with the version pair 0 1, a tag count and the first tag, a
// length-prefixed name. Binary: little-endian int32s, one per character too.
bool sniff_is_mphbin(const std::string& rHead) {
    if (rHead.size() < 20)
        return false;
    auto i32 = [&](std::size_t k) {
        const unsigned char* p = reinterpret_cast<const unsigned char*>(rHead.data() + k);
        return static_cast<std::int64_t>(static_cast<std::int32_t>(
            static_cast<std::uint32_t>(p[0]) | (static_cast<std::uint32_t>(p[1]) << 8) |
            (static_cast<std::uint32_t>(p[2]) << 16) | (static_cast<std::uint32_t>(p[3]) << 24)));
    };
    const std::int64_t first_char = i32(16);
    return i32(0) == 0 && i32(4) == 1 && i32(8) >= 1 && i32(8) <= 4096 && i32(12) >= 1 &&
           i32(12) <= 1024 && first_char > 32 && first_char < 127;
}

// Text: the same values as tokens, `#` comments skipped.
bool sniff_is_mphtxt(const std::string& rHead) {
    std::vector<std::string> tokens;
    std::size_t k = 0;
    while (k < rHead.size() && tokens.size() < 5) {
        const char c = rHead[k];
        if (c == '#') {
            while (k < rHead.size() && rHead[k] != '\n')
                ++k;
        } else if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            ++k;
        } else {
            const std::size_t b = k;
            while (k < rHead.size() && rHead[k] != ' ' && rHead[k] != '\t' && rHead[k] != '\n' &&
                   rHead[k] != '\r' && rHead[k] != '#')
                ++k;
            tokens.push_back(rHead.substr(b, k - b));
        }
    }
    auto count = [](const std::string& rT) {
        return !rT.empty() && rT.size() < 6 &&
               rT.find_first_not_of("0123456789") == std::string::npos && rT != "0";
    };
    return tokens.size() == 5 && tokens[0] == "0" && tokens[1] == "1" && count(tokens[2]) &&
           count(tokens[3]) &&
           ((tokens[4][0] >= 'a' && tokens[4][0] <= 'z') ||
            (tokens[4][0] >= 'A' && tokens[4][0] <= 'Z'));
}

// A directory-shaped mesh, recognised by the files the readers themselves look
// for: an Elmer mesh directory (`mesh.header`, a `partitioning.N` directory of
// `part.n.*` files, or a directory holding one) and an OpenFOAM case (the
// layouts the openfoam reader resolves). A directory that looks like both, or
// like neither, is "".
// Femap neutral file: a lone `-1` line (on the first or second line -- MYSTRAN
// writes a number before it), then the header block's id, 100.
bool sniff_is_femap(const std::string& rHead) {
    std::vector<std::string> lines;
    std::size_t pos = 0;
    while (pos < rHead.size() && lines.size() < 3) {
        std::size_t eol = rHead.find('\n', pos);
        if (eol == std::string::npos)
            break;
        std::string line = rHead.substr(pos, eol - pos);
        const std::size_t b = line.find_first_not_of(" \t\r");
        const std::size_t e = line.find_last_not_of(" \t\r");
        lines.push_back(b == std::string::npos ? std::string() : line.substr(b, e - b + 1));
        pos = eol + 1;
    }
    for (std::size_t k = 0; k + 1 < lines.size() && k < 2; ++k)
        if (lines[k] == "-1" && lines[k + 1] == "100")
            return true;
    return false;
}

// Patran 2 neutral file: the first card is a title (25) or summary (26) packet
// header in the fixed `(I2,8I8)` columns -- every field right-justified digits --
// announcing at least one data card. Matched on the unstripped head: the columns
// are the signature.
bool sniff_is_patran(const std::string& rHead) {
    if (rHead.size() < 26 || rHead[0] != '2' || (rHead[1] != '5' && rHead[1] != '6'))
        return false;
    std::size_t eol = rHead.find('\n');
    if (eol == std::string::npos)
        eol = rHead.size();
    std::string line = rHead.substr(0, eol);
    if (!line.empty() && line.back() == '\r')
        line.pop_back();
    if (line.size() < 26 || line.size() > 80)
        return false;
    std::int64_t kc = 0;
    for (std::size_t f = 0; 2 + 8 * f < line.size() && f < 8; ++f) {
        const std::string field = line.substr(2 + 8 * f, 8);
        const std::size_t first = field.find_first_not_of(' ');
        if (first == std::string::npos)
            return false;
        std::size_t k = first + (field[first] == '-' ? 1 : 0);
        if (k >= field.size())
            return false;
        for (; k < field.size(); ++k)
            if (field[k] < '0' || field[k] > '9')
                return false;
        if (f == 2)
            kc = std::stoll(field.substr(first));
    }
    return kc >= 1;
}

// Nastran OP2 written with PARAM,POST,-1: Fortran blocks (4-byte markers, either
// byte order) holding a one-word 3, a 3-word date, a one-word 7 and the 7-word
// tape code, in 4- or 8-byte words.
bool sniff_is_op2(const std::string& rHead) {
    for (bool big : {false, true}) {
        const auto u32 = [&](std::size_t At) -> std::int64_t {
            if (At + 4 > rHead.size())
                return -1;
            std::uint32_t v;
            std::memcpy(&v, rHead.data() + At, 4);
            if (big != (std::endian::native == std::endian::big))
                v = detail::bswap32(v);
            return static_cast<std::int32_t>(v);
        };
        // Fortran blocks with 4-byte markers: [n][payload][n].
        std::size_t pos = 0;
        std::vector<std::pair<std::size_t, std::int64_t>> blocks;  // (payload offset, size)
        while (blocks.size() < 4) {
            const std::int64_t n = u32(pos);
            if (n <= 0 || u32(pos + 4 + static_cast<std::size_t>(n)) != n)
                break;
            blocks.emplace_back(pos + 4, n);
            pos += 8 + static_cast<std::size_t>(n);
        }
        if (blocks.size() < 4)
            continue;
        const std::int64_t ws = blocks[0].second;
        if (ws != 4 && ws != 8)
            continue;
        const auto word = [&](std::size_t Block) {
            if (ws == 4)
                return u32(blocks[Block].first);
            std::uint64_t v;
            std::memcpy(&v, rHead.data() + blocks[Block].first, 8);
            if (big != (std::endian::native == std::endian::big))
                v = detail::bswap64(v);
            return static_cast<std::int64_t>(v);
        };
        if (word(0) == 3 && blocks[1].second == 3 * ws && blocks[2].second == ws && word(2) == 7 &&
            blocks[3].second == 7 * ws)
            return true;
    }
    return false;
}

// Abaqus results file, binary: a 4096-byte Fortran record (marker 4096 in
// either byte order) whose first 8-byte word is a record length and whose second
// is the key of a record Abaqus writes first (1921 release, 1922 heading, 1900
// element, 1901 node, 2000 increment start).
bool sniff_is_fil_key(std::int64_t Key) {
    return Key == 1921 || Key == 1922 || Key == 1900 || Key == 1901 || Key == 2000;
}

bool sniff_is_abaqus_fil_binary(const std::string& rHead) {
    if (rHead.size() < 20)
        return false;
    const auto u32 = [&](std::size_t k, bool big) {
        std::uint32_t v = 0;
        for (std::size_t b = 0; b < 4; ++b) {
            const std::uint32_t byte = static_cast<unsigned char>(rHead[k + b]);
            v |= big ? byte << (8 * (3 - b)) : byte << (8 * b);
        }
        return v;
    };
    for (const bool big : {false, true}) {
        if (u32(0, big) != 4096)
            continue;
        // An integer fills the first four bytes of its word (Fortran
        // EQUIVALENCE), or the low half of an 8-byte integer: the last four
        // bytes when big-endian.
        for (const std::size_t lo : {std::size_t{0}, std::size_t{4}}) {
            if (lo == 4 && !big)
                continue;
            const std::uint32_t length = u32(4 + lo, big);
            const std::uint32_t key = u32(12 + lo, big);
            if (length >= 2 && length <= 512 && sniff_is_fil_key(key))
                return true;
        }
    }
    return false;
}

// Abaqus results file, ASCII: a `*` record start, then the record length and key
// as `I` items (`I` + two-digit digit count + the digits).
bool sniff_is_abaqus_fil_ascii(const std::string& rStripped) {
    if (rStripped.size() < 12 || rStripped[0] != '*')
        return false;
    std::size_t pos = 1;
    std::int64_t items[2] = {0, 0};
    for (std::int64_t& item : items) {
        if (pos + 3 > rStripped.size() || rStripped[pos] != 'I')
            return false;
        const std::string width = rStripped.substr(pos + 1, 2);
        const std::size_t first = width.find_first_not_of(' ');
        if (first == std::string::npos || width[1] < '0' || width[1] > '9')
            return false;
        const int n = std::stoi(width.substr(first));
        if (n < 1 || n > 18 || pos + 3 + static_cast<std::size_t>(n) > rStripped.size())
            return false;
        const std::string digits = rStripped.substr(pos + 3, static_cast<std::size_t>(n));
        for (char c : digits)
            if (c < '0' || c > '9')
                return false;
        item = std::stoll(digits);
        pos += 3 + static_cast<std::size_t>(n);
    }
    return items[0] >= 2 && sniff_is_fil_key(items[1]);
}

// OpenRadioss starter deck: `#RADIOSS STARTER`, or `/BEGIN` as the first line
// that is not a `#`/`$` comment.
bool sniff_is_radioss(const std::string& rStripped) {
    std::size_t pos = 0;
    while (pos < rStripped.size()) {
        std::size_t eol = rStripped.find('\n', pos);
        if (eol == std::string::npos)
            eol = rStripped.size();
        const std::string line = rStripped.substr(pos, eol - pos);
        pos = eol + 1;
        if (sniff_starts_with(line, "#RADIOSS STARTER"))
            return true;
        if (line.empty() || line[0] == '#' || line[0] == '$' || line[0] == '\r')
            continue;
        return sniff_starts_with(line, "/BEGIN");
    }
    return false;
}

std::string sniff_directory(const std::filesystem::path& rDir) {
    namespace fs = std::filesystem;
    std::error_code ec;
    const auto is_file = [&](const fs::path& rPath) { return fs::is_regular_file(rPath, ec); };
    const auto is_dir = [&](const fs::path& rPath) { return fs::is_directory(rPath, ec); };
    const auto has_polymesh = [&](const fs::path& rPoly) {
        return is_file(rPoly / "owner") && is_file(rPoly / "faces");
    };

    const auto is_partitioning = [&](const fs::path& rPath) {
        return rPath.filename().string().rfind("partitioning.", 0) == 0 &&
               is_file(rPath / "part.1.header");
    };
    bool elmer = is_file(rDir / "mesh.header") || is_partitioning(rDir);
    for (auto it = fs::directory_iterator(rDir, ec);
         !elmer && !ec && it != fs::directory_iterator(); it.increment(ec))
        elmer = is_partitioning(it->path());
    const bool openfoam = (rDir.filename() == "polyMesh" && has_polymesh(rDir)) ||
                          has_polymesh(rDir / "constant" / "polyMesh") ||
                          has_polymesh(rDir / "polyMesh") ||
                          is_dir(rDir / "processor0" / "constant" / "polyMesh") ||
                          is_file(rDir / "constant" / "regionProperties");
    if (elmer == openfoam)
        return "";
    return elmer ? "elmer" : "openfoam";
}

}  // namespace

std::string sniff_format(const std::string& rPath) {
    {
        namespace fs = std::filesystem;
        std::error_code ec;
        const fs::path path(rPath);
        if (fs::is_directory(path, ec))
            return sniff_directory(path);
        // The header of an Elmer mesh directory stands for the directory.
        if (path.filename() == "mesh.header" && fs::is_regular_file(path, ec))
            return "elmer";
        // Z88's input and output files have fixed names.
        if (is_z88_filename(rPath) && fs::is_regular_file(path, ec))
            return "z88";
        if (is_d3plot_filename(rPath) && fs::is_regular_file(path, ec))
            return "lsdyna_d3plot";
    }
    auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
    if (!in)
        return "";
    char buf[512];
    in.read(buf, sizeof(buf));
    const std::string head(buf, static_cast<std::size_t>(in.gcount()));
    if (head.empty())
        return "";
    const std::string stripped = sniff_lstrip(head);

    // --- binary magics ---
    if (sniff_is_mphbin(head))
        return "mphbin";
    // Tecplot binary (.plt): "#!TDV" and a three-character version.
    if (head.size() >= 5 && head.compare(0, 5, "#!TDV") == 0)
        return "tecplot";
    // FEBio plot file: the magic 0x00464542, in either byte order.
    if (head.size() >= 4 && (head.compare(0, 4, std::string("BEF\0", 4)) == 0 ||
                             head.compare(0, 4, std::string("\0FEB", 4)) == 0))
        return "xplt";
    // Ansys MAPDL results: a 100-word integer record (length 100, flags
    // 0x80000000) whose first value is the file number, 12.
    if (head.size() >= 12 &&
        head.compare(0, 12, std::string("d\0\0\0\0\0\0\x80\x0c\0\0\0", 12)) == 0)
        return "ansys_rst";
    // libMesh XDR mesh: the version string, a big-endian length then "libMesh-".
    if (head.size() >= 12 && head.compare(0, 2, std::string("\0\0", 2)) == 0 &&
        head.compare(4, 8, "libMesh-") == 0)
        return "libmesh";
    if (sniff_is_abaqus_fil_binary(head))
        return "abaqus_fil";
    // Nastran OP2 (PARAM,POST,-1): a one-word record 3, a 3-word date, a
    // one-word record 7 and the 7-word tape code, any word size and order.
    if (sniff_is_op2(head))
        return "nastran_op2";
    // LS-DYNA d3plot: a plausible 64-word control block, any word size and order.
    if (is_d3plot_head(head.data(), head.size()))
        return "lsdyna_d3plot";
    // FEBio input: XML whose root is <febio_spec>.
    if (sniff_contains(head, "<febio_spec"))
        return "febio";
    // VTK XML formats begin (possibly after a BOM/whitespace) with "<?xml" or
    // directly a "<VTKFile" element carrying the grid type.
    if (sniff_contains(head, "VTKFile")) {
        // The parallel indices and the collection come first, and match the
        // quoted `type=` value: `PUnstructuredGrid` contains `UnstructuredGrid`
        // and `PPolyData` contains `PolyData`, so the loose substring checks
        // below would call an index a piece; and a bare `Collection` would also
        // match `vtkPartitionedDataSetCollection`. A parallel image, structured
        // or rectilinear index is refused outright rather than mistaken for its
        // serial twin. (v15.0.0; mirrors `_sniff.py`.)
        static const struct {
            const char* mValue;
            const char* mFormat;
        } kIndexTypes[] = {
            {"Collection", "pvd"}, {"PUnstructuredGrid", "pvtu"}, {"PPolyData", "pvtp"},
            {"PImageData", ""},    {"PStructuredGrid", ""},       {"PRectilinearGrid", ""},
        };
        for (const char quote : {'"', '\''})
            for (const auto& index_type : kIndexTypes) {
                const std::string needle = std::string("type=") + quote + index_type.mValue + quote;
                if (sniff_contains(head, needle.c_str()))
                    return index_type.mFormat;
            }
        if (sniff_contains(head, "UnstructuredGrid"))
            return "vtu";
        if (sniff_contains(head, "PolyData"))
            return "vtp";
        // Checked last of the four: the grid-type strings are disjoint, but a
        // future dataset type could contain another as a substring, and the
        // cheapest defence is to keep the most recently added one from
        // shadowing anything.
        if (sniff_contains(head, "ImageData"))
            return "vti";
        // v11.6.0, roadmap §1 tier B4.
        if (sniff_contains(head, "StructuredGrid"))
            return "vts";
        if (sniff_contains(head, "RectilinearGrid"))
            return "vtr";
        if (sniff_contains(head, "MultiBlockDataSet"))
            return "vtm";
    }
    if (sniff_starts_with(stripped, "<Xdmf") || sniff_contains(head, "<Xdmf"))
        return "xdmf";

    // --- line/text signatures ---
    if (sniff_starts_with(stripped, "# vtk DataFile"))
        return "vtk";
    if (sniff_starts_with(stripped, "$MeshFormat"))
        return "gmsh";
    // GiD postprocess. The results file is self-identifying. The geometry file
    // is matched on `MESH "` -- keyword, space AND opening quote -- because a
    // bare "MESH " prefix is exactly the generic English token this file's own
    // contract warns against claiming; the quote is what gidpost always writes
    // and what makes the match unambiguous. Neither `.post.bin` (a deflated
    // stream, no stable leading signature) nor `.post.h5` (the generic HDF5
    // magic, which this file deliberately never claims) is sniffable.
    if (sniff_starts_with(stripped, "GiD Post Results File"))
        return "gid";
    if (sniff_starts_with(stripped, "MESH \""))
        return "gid";
    if (sniff_is_mphtxt(head))
        return "mphtxt";
    // PLY: "ply" on its own first line.
    if (sniff_starts_with(stripped, "ply\n") || sniff_starts_with(stripped, "ply\r") ||
        stripped == "ply")
        return "ply";
    // OFF variants (OFF / COFF / NOFF / STOFF ...): a token ending in "OFF".
    if (sniff_starts_with(stripped, "OFF") || sniff_starts_with(stripped, "COFF") ||
        sniff_starts_with(stripped, "NOFF") || sniff_starts_with(stripped, "STOFF"))
        return "off";
    // PCL point clouds: the writer's first line is "# .PCD v0.7 - ...", and a header with
    // its comments stripped opens on VERSION.
    if (sniff_starts_with(stripped, "# .PCD") ||
        (sniff_starts_with(stripped, "VERSION") && sniff_contains(head, "\nFIELDS")))
        return "pcd";
    // CalculiX results: a lone "    1C" record, then the "1U" user header (or the "2C"
    // node block of a file without one). Matched on the unstripped head, so the fixed
    // columns are part of the signature.
    if (sniff_starts_with(head, "    1C")) {
        std::size_t eol = 6;
        while (eol < head.size() && (head[eol] == ' ' || head[eol] == '\r'))
            ++eol;
        if (eol < head.size() && head[eol] == '\n') {
            const std::string next = head.substr(eol + 1, 6);
            if (sniff_starts_with(next, "    1U") || sniff_starts_with(next, "    2C"))
                return "frd";
        }
    }
    // I-DEAS universal file: a lone "-1" line, then a known dataset number (a trailing
    // "b" marks a binary dataset such as 58b).
    if (sniff_starts_with(stripped, "-1")) {
        std::size_t eol = 2;
        while (eol < stripped.size() && (stripped[eol] == ' ' || stripped[eol] == '\r'))
            ++eol;
        if (eol < stripped.size() && stripped[eol] == '\n') {
            std::size_t a = eol + 1;
            while (a < stripped.size() && stripped[a] == ' ')
                ++a;
            std::size_t b = a;
            while (b < stripped.size() && stripped[b] >= '0' && stripped[b] <= '9')
                ++b;
            static const char* const kUnvIds[] = {"15",   "18",   "55",   "56",   "57",   "58",
                                                  "82",   "151",  "164",  "780",  "781",  "2400",
                                                  "2411", "2412", "2414", "2417", "2420", "2429",
                                                  "2430", "2432", "2435", "2452", "2467", "2477"};
            const std::string id = stripped.substr(a, b - a);
            const bool ends = b >= stripped.size() || stripped[b] == ' ' || stripped[b] == '\r' ||
                              stripped[b] == '\n' || stripped[b] == 'b' || stripped[b] == 'B';
            for (const char* known : kUnvIds)
                if (ends && id == known)
                    return "unv";
        }
    }
    // MFEM mesh: the first line names it (`MFEM mesh v1.x`; the non-conforming,
    // NURBS and INLINE kinds are MFEM's too, and read_mfem names why it refuses them).
    if (sniff_starts_with(stripped, "MFEM mesh v1.") ||
        sniff_starts_with(stripped, "MFEM NC mesh") ||
        sniff_starts_with(stripped, "MFEM NURBS mesh") ||
        sniff_starts_with(stripped, "MFEM INLINE mesh"))
        return "mfem";
    if (sniff_starts_with(stripped, "libMesh-"))
        return "libmesh";
    if (sniff_is_abaqus_fil_ascii(stripped))
        return "abaqus_fil";
    if (sniff_is_radioss(stripped))
        return "radioss";
    // MSC Marc formatted post file: the analysis title block opens it.
    if (sniff_starts_with(stripped, "=beg=50100"))
        return "marc_t19";
    // MSC Marc input deck: a Marc parameter opens it and an END, CONNECTIVITY or
    // COORDINATES line follows -- which can be past 512 bytes, so read on.
    if (!stripped.empty() && std::isalpha(static_cast<unsigned char>(stripped[0])) &&
        stripped.substr(0, stripped.find('\n')).find('=') == std::string::npos) {
        auto deck = detail::make_classic_ifstream(rPath, std::ios::binary);
        std::string text(65536, '\0');
        deck.read(text.data(), static_cast<std::streamsize>(text.size()));
        text.resize(static_cast<std::size_t>(deck.gcount()));
        if (is_marc_deck(text))
            return "marc";
    }
    if (sniff_is_femap(head))
        return "femap";
    if (sniff_is_patran(head))
        return "patran";
    // ASCII STL.
    if (sniff_starts_with(stripped, "solid "))
        return "stl";
    // Code_Aster .mail meshes open, after any `%` comment lines, with a TITRE or
    // COOR_1D/2D/3D block keyword.
    {
        std::size_t pos = 0;
        while (pos < stripped.size()) {
            std::size_t eol = stripped.find('\n', pos);
            if (eol == std::string::npos)
                eol = stripped.size();
            std::string line = stripped.substr(pos, eol - pos);
            pos = eol + 1;
            const std::size_t first = line.find_first_not_of(" \t\r");
            if (first == std::string::npos || line[first] == '%')
                continue;
            line = line.substr(first);
            std::transform(line.begin(), line.end(), line.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
            const std::size_t end = line.find_first_of(" \t\r,%");
            const std::string word = line.substr(0, end);
            if (word == "TITRE" || word == "COOR_1D" || word == "COOR_2D" || word == "COOR_3D")
                return "code_aster";
            break;
        }
    }
    // LS-DYNA keyword decks open with "*KEYWORD" after any `$` comment lines; this
    // runs before the Abaqus rule so a deck that opens with a keyword line other
    // than *NODE is never mistaken for one.
    {
        std::size_t pos = 0;
        while (pos < stripped.size()) {
            std::size_t eol = stripped.find('\n', pos);
            if (eol == std::string::npos)
                eol = stripped.size();
            const std::string line = stripped.substr(pos, eol - pos);
            pos = eol + 1;
            if (line.empty() || line[0] == '$' || line[0] == '\r')
                continue;
            std::string upper = line.substr(0, 8);
            std::transform(upper.begin(), upper.end(), upper.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
            if (sniff_starts_with(upper, "*KEYWORD"))
                return "lsdyna";
            break;
        }
    }
    // Abaqus input decks start with a keyword line "*Heading"/"*Node"/"*NODE".
    {
        std::string upper = stripped.substr(0, 8);
        std::transform(upper.begin(), upper.end(), upper.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        if (sniff_starts_with(upper, "*HEADING") || sniff_starts_with(upper, "*NODE"))
            return "abaqus";
    }

    return "";
}

}  // namespace meshioplusplus
