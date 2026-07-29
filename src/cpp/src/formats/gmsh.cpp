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
#include <functional>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/types.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

namespace {

// ---- type maps (subset; ported from gmsh/common.py) --------------------------
const std::unordered_map<int, std::string>& gmsh_to_meshio_type() {
    static const std::unordered_map<int, std::string> m = {
        {1, "line"},          {2, "triangle"},      {3, "quad"},           {4, "tetra"},
        {5, "hexahedron"},    {6, "wedge"},         {7, "pyramid"},        {8, "line3"},
        {9, "triangle6"},     {10, "quad9"},        {11, "tetra10"},       {12, "hexahedron27"},
        {13, "wedge18"},      {14, "pyramid14"},    {15, "vertex"},        {16, "quad8"},
        {17, "hexahedron20"}, {18, "wedge15"},      {19, "pyramid13"},     {21, "triangle10"},
        {23, "triangle15"},   {25, "triangle21"},   {26, "line4"},         {27, "line5"},
        {28, "line6"},        {29, "tetra20"},      {30, "tetra35"},       {31, "tetra56"},
        {36, "quad16"},       {37, "quad25"},       {38, "quad36"},        {62, "line7"},
        {63, "line8"},        {64, "line9"},        {65, "line10"},        {66, "line11"},
        {71, "tetra84"},      {72, "tetra120"},     {73, "tetra165"},      {74, "tetra220"},
        {75, "tetra286"},     {92, "hexahedron64"}, {93, "hexahedron125"},
    };
    return m;
}

const std::unordered_map<std::string, int>& meshio_to_gmsh_type() {
    static const std::unordered_map<std::string, int> m = [] {
        std::unordered_map<std::string, int> r;
        for (const auto& kv : gmsh_to_meshio_type())
            r[kv.second] = kv.first;
        return r;
    }();
    return m;
}

// Permutation P such that meshio_row[j] = gmsh_row[P[j]]; empty = identity.
const std::vector<int>& gmsh_to_meshio_perm(const std::string& rT) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"tetra10", {0, 1, 2, 3, 4, 5, 6, 7, 9, 8}},
        {"hexahedron20", {0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17, 10, 12, 14, 15}},
        {"hexahedron27", {0,  1,  2,  3,  4,  5,  6,  7,  8,  11, 13, 9,  16, 18,
                          19, 17, 10, 12, 14, 15, 22, 23, 21, 24, 20, 25, 26}},
        {"wedge15", {0, 1, 2, 3, 4, 5, 6, 9, 7, 12, 14, 13, 8, 10, 11}},
        {"pyramid13", {0, 1, 2, 3, 4, 5, 8, 10, 6, 7, 9, 11, 12}},
    };
    static const std::vector<int> empty;
    auto it = m.find(rT);
    return it == m.end() ? empty : it->second;
}

const std::vector<int>& meshio_to_gmsh_perm(const std::string& rT) {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"tetra10", {0, 1, 2, 3, 4, 5, 6, 7, 9, 8}},
        {"hexahedron20", {0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 16, 9, 17, 10, 18, 19, 12, 15, 13, 14}},
        {"hexahedron27", {0,  1,  2,  3,  4,  5,  6,  7,  8,  11, 16, 9,  17, 10,
                          18, 19, 12, 15, 13, 14, 24, 22, 20, 21, 23, 25, 26}},
        {"wedge15", {0, 1, 2, 3, 4, 5, 6, 8, 12, 7, 13, 14, 9, 11, 10}},
        {"pyramid13", {0, 1, 2, 3, 4, 5, 8, 9, 6, 10, 7, 11, 12}},
    };
    static const std::vector<int> empty;
    auto it = m.find(rT);
    return it == m.end() ? empty : it->second;
}

std::string gmsh_trim(const std::string& rS) {
    std::size_t b = 0, e = rS.size();
    while (b < e && std::isspace(static_cast<unsigned char>(rS[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(rS[e - 1])))
        --e;
    return rS.substr(b, e - b);
}

struct GmshCursor {
    // A view, not a reference to a std::string: the buffer may be a memory
    // mapping rather than an owned string (see detail/file_source.hpp).
    std::string_view mBuf;
    std::size_t mPos = 0;
    explicit GmshCursor(std::string_view b) : mBuf(b) {}
    bool eof() const { return mPos >= mBuf.size(); }

    std::string read_line() {
        std::size_t start = mPos;
        while (mPos < mBuf.size() && mBuf[mPos] != '\n')
            ++mPos;
        std::string line(mBuf.substr(start, mPos - start));
        if (mPos < mBuf.size())
            ++mPos;
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        return line;
    }
    std::string next_nonblank() {
        while (!eof()) {
            std::string l = read_line();
            if (!gmsh_trim(l).empty())
                return l;
        }
        return "";
    }
    void skip_to_end(const std::string& rEnv) {
        std::string target = "$End" + rEnv;
        while (!eof()) {
            if (gmsh_trim(read_line()) == target)
                return;
        }
    }
    double next_double() {
        // strtod scans for a terminator. A buffered source is a std::string
        // (NUL-terminated); a mapped one relies on the kernel zero-filling the
        // final partial page -- which is exactly why FileSource declines to map
        // files whose size is an exact page multiple.
        const char* base = mBuf.data();
        char* endp = nullptr;
        double v = std::strtod(base + mPos, &endp);
        if (endp == base + mPos)
            throw ReadError("Gmsh: expected a number");
        mPos = static_cast<std::size_t>(endp - base);
        return v;
    }
    std::int64_t next_int() { return static_cast<std::int64_t>(next_double()); }

    std::int32_t read_i32() {
        std::int32_t v;
        std::memcpy(&v, mBuf.data() + mPos, 4);
        mPos += 4;
        return v;
    }
    double read_f64() {
        double v;
        std::memcpy(&v, mBuf.data() + mPos, 8);
        mPos += 8;
        return v;
    }
    // Read an unsigned integer of `sz` bytes (little-endian host).
    std::uint64_t read_uint(int sz) {
        std::uint64_t v = 0;
        std::memcpy(&v, mBuf.data() + mPos, static_cast<std::size_t>(sz));
        mPos += static_cast<std::size_t>(sz);
        return v;
    }
};

struct EBlock {
    std::string mType;
    std::size_t mN = 0;
    std::size_t mCount = 0;
    std::size_t mNumTags = 0;
    std::vector<std::int64_t> mConn;  // count*n, 0-based gmsh ids
    std::vector<std::int64_t> mTags;  // count*num_tags
};

void read_physical_names(GmshCursor& rCur, std::unordered_map<std::string, NDArray>& rFieldData) {
    std::int64_t num = std::stoll(gmsh_trim(rCur.read_line()));
    for (std::int64_t i = 0; i < num; ++i) {
        std::string line = rCur.read_line();
        std::istringstream iss(line);
        long long dim, tag;
        iss >> dim >> tag;
        std::size_t q1 = line.find('"');
        std::size_t q2 = line.rfind('"');
        std::string name =
            (q1 != std::string::npos && q2 > q1) ? line.substr(q1 + 1, q2 - q1 - 1) : "";
        NDArray v(DType::Int64, {2});
        v.As<std::int64_t>()[0] = tag;  // physical number
        v.As<std::int64_t>()[1] = dim;
        rFieldData.emplace(name, std::move(v));
    }
    rCur.skip_to_end("PhysicalNames");
}

// --- physical groups <-> named regions ---------------------------------------
//
// Gmsh describes a physical group twice over: `$PhysicalNames` gives it a name,
// a dimension and an integer tag (which meshio++ stores as `field_data[name] =
// [tag, dim]`), and the per-element tag column says which cells belong to it
// (stored as the `gmsh:physical` cell_data). Both of those stay exactly as they
// were — every existing consumer is unaffected, and a mesh that carries them
// writes byte-identical bytes. A `Region` is *derived* from the pair on read,
// and consulted only to fill gaps on write. See doc/regions.md for the
// precedence rule this implements.

/// `$PhysicalNames` as `(tag, dim) -> name`, read back out of `field_data`.
std::map<std::pair<std::int64_t, int>, std::string> gmsh_physical_names(const Mesh& rMesh) {
    std::map<std::pair<std::int64_t, int>, std::string> out;
    for (const auto& name : rMesh.FieldDataNames()) {
        const NDArray& d = rMesh.FieldData(name);
        if (d.Size() < 2)
            continue;
        const std::int64_t tag = detail::read_int(d, 0);
        const int dim = static_cast<int>(detail::read_int(d, 1));
        out.emplace(std::make_pair(tag, dim), name);
    }
    return out;
}

/**
 * @brief Derive one `Cell` region per named physical group.
 *
 * The region's `mTag` is the gmsh physical tag and its `mDim` the group's
 * dimension, so a round-trip back to gmsh can restore both. Groups the file
 * never named in `$PhysicalNames` get no region — their tag still lives in the
 * `gmsh:physical` cell_data, which is untouched.
 *
 * Note that a dimension-0 physical group tags `vertex` *cells* in gmsh, not
 * points, so it too becomes a `Cell` region (with `mDim == 0`) rather than a
 * `Point` one.
 */
void gmsh_attach_regions(Mesh& rMesh) {
    if (!rMesh.HasCellData("gmsh:physical"))
        return;
    const std::map<std::pair<std::int64_t, int>, std::string> names = gmsh_physical_names(rMesh);
    if (names.empty())
        return;
    if (rMesh.CellDataNumBlocks("gmsh:physical") != rMesh.NumCellBlocks())
        return;  // partial tag data cannot be aligned with the blocks

    // Group cells by (tag, block dimension); the dimension disambiguates two
    // physical groups that legitimately share a tag across dimensions.
    std::map<std::pair<std::int64_t, int>, std::vector<std::int64_t>> members;
    std::size_t b = 0;
    std::int64_t base = 0;
    for (const auto cb : rMesh.CellRange()) {
        const NDArray& tags = rMesh.CellData("gmsh:physical", b);
        const std::size_t ncells = cb.NumCells();
        int dim = cell_type_dimension(cell_type_from_name(std::string(cb.Type())));
        if (dim < 0) {
            auto it = topological_dimension().find(std::string(cb.Type()));
            dim = it != topological_dimension().end() ? it->second : -1;
        }
        if (tags.Size() >= ncells)
            for (std::size_t c = 0; c < ncells; ++c)
                members[{detail::read_int(tags, c), dim}].push_back(base +
                                                                    static_cast<std::int64_t>(c));
        base += static_cast<std::int64_t>(ncells);
        ++b;
    }

    for (const auto& [key, name] : names) {
        auto it = members.find(key);
        if (it == members.end())
            continue;
        NDArray entries = NDArray::Uninit(DType::Int64, {it->second.size()});
        for (std::size_t i = 0; i < it->second.size(); ++i)
            entries.As<std::int64_t>()[i] = it->second[i];
        rMesh.AddRegion(Region(name, RegionKind::Cell, key.second, key.first, std::move(entries)));
    }
}

/**
 * @brief The `$PhysicalNames` rows to write: `field_data` first, then any
 * region that describes a group `field_data` does not.
 *
 * `field_data` winning is what keeps output byte-identical for every mesh that
 * already carried gmsh's own metadata; regions only add groups that came from
 * another format.
 * @return `(dim, tag, name)` rows, sorted — the order the writer emits.
 */
std::vector<std::tuple<long long, long long, std::string>> gmsh_physical_rows(const Mesh& rMesh) {
    std::vector<std::tuple<long long, long long, std::string>> rows;
    std::set<std::string> seen;
    for (const auto& name : rMesh.FieldDataNames()) {
        const NDArray& d = rMesh.FieldData(name);
        if (d.Size() < 2)
            continue;
        rows.emplace_back(detail::read_int(d, 1), detail::read_int(d, 0), name);
        seen.insert(name);
    }
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Cell || r.mTag < 0 || seen.count(r.mName))
            continue;
        rows.emplace_back(r.mDim, r.mTag, r.mName);
        seen.insert(r.mName);
    }
    std::sort(rows.begin(), rows.end());
    return rows;
}

/**
 * @brief Per-block `gmsh:physical` tag arrays synthesized from `Cell` regions.
 *
 * Only used when the mesh carries no `gmsh:physical` cell_data of its own — a
 * mesh read from another format. Cells in no tagged region get tag 0, which is
 * gmsh's "no physical group".
 * @return one Int64 array per cell block, or an empty vector when there is
 *         nothing to synthesize.
 */
std::vector<NDArray> gmsh_tags_from_regions(const Mesh& rMesh) {
    std::vector<NDArray> blocks;
    bool any = false;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind == RegionKind::Cell && r.mTag >= 0 && r.NumEntries() > 0)
            any = true;
    }
    if (!any)
        return blocks;

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    std::vector<std::int64_t> flat(static_cast<std::size_t>(detail::total_cells(bases)), 0);
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Cell || r.mTag < 0)
            continue;
        const std::int64_t* e = r.Entries();
        for (std::size_t k = 0; k < r.NumEntries(); ++k)
            if (e[k] >= 0 && e[k] < static_cast<std::int64_t>(flat.size()))
                flat[static_cast<std::size_t>(e[k])] = r.mTag;
    }
    blocks.reserve(rMesh.NumCellBlocks());
    for (std::size_t b = 0; b + 1 < bases.size(); ++b) {
        const std::size_t n = static_cast<std::size_t>(bases[b + 1] - bases[b]);
        NDArray a = NDArray::Uninit(DType::Int64, {n});
        for (std::size_t c = 0; c < n; ++c)
            a.As<std::int64_t>()[c] = flat[static_cast<std::size_t>(bases[b]) + c];
        blocks.push_back(std::move(a));
    }
    return blocks;
}

void read_nodes(GmshCursor& rCur, bool is_ascii, NDArray& rPoints,
                std::vector<std::int64_t>& rPointTags) {
    std::int64_t num = std::stoll(gmsh_trim(rCur.read_line()));
    rPoints = NDArray(DType::Float64, {static_cast<std::size_t>(num), 3});
    rPointTags.resize(num);
    double* pp = rPoints.As<double>();
    if (is_ascii) {
        for (std::int64_t i = 0; i < num; ++i) {
            rPointTags[i] = rCur.next_int();
            pp[i * 3 + 0] = rCur.next_double();
            pp[i * 3 + 1] = rCur.next_double();
            pp[i * 3 + 2] = rCur.next_double();
        }
    } else {
        for (std::int64_t i = 0; i < num; ++i) {
            rPointTags[i] = rCur.read_i32();
            pp[i * 3 + 0] = rCur.read_f64();
            pp[i * 3 + 1] = rCur.read_f64();
            pp[i * 3 + 2] = rCur.read_f64();
        }
    }
    rCur.skip_to_end("Nodes");
}

void append_element(std::vector<EBlock>& rBlocks, const std::string& rType, std::size_t n,
                    std::size_t num_tags, const std::int64_t* pTags, const std::int64_t* pNodes) {
    if (rBlocks.empty() || rBlocks.back().mType != rType || rBlocks.back().mNumTags != num_tags) {
        EBlock b;
        b.mType = rType;
        b.mN = n;
        b.mNumTags = num_tags;
        rBlocks.push_back(std::move(b));
    }
    EBlock& cur = rBlocks.back();
    for (std::size_t j = 0; j < num_tags; ++j)
        cur.mTags.push_back(pTags[j]);
    for (std::size_t j = 0; j < n; ++j)
        cur.mConn.push_back(pNodes[j] - 1);
    ++cur.mCount;
}

void read_elements(GmshCursor& rCur, bool is_ascii, std::vector<EBlock>& rBlocks) {
    std::int64_t total = std::stoll(gmsh_trim(rCur.read_line()));
    const auto& g2m = gmsh_to_meshio_type();
    const auto& nnpc = num_nodes_per_cell();

    if (is_ascii) {
        for (std::int64_t e = 0; e < total; ++e) {
            std::string line = rCur.read_line();
            std::istringstream iss(line);
            std::vector<std::int64_t> v;
            long long x;
            while (iss >> x)
                v.push_back(x);
            int gtype = static_cast<int>(v[1]);
            std::size_t num_tags = static_cast<std::size_t>(v[2]);
            auto it = g2m.find(gtype);
            if (it == g2m.end())
                throw ReadError("Gmsh element type " + std::to_string(gtype) +
                                " not supported by the C++ reader");
            std::size_t n = static_cast<std::size_t>(nnpc.at(it->second));
            append_element(rBlocks, it->second, n, num_tags, v.data() + 3, v.data() + 3 + num_tags);
        }
    } else {
        std::int64_t done = 0;
        while (done < total) {
            int gtype = rCur.read_i32();
            std::int32_t nelem = rCur.read_i32();
            std::int32_t num_tags = rCur.read_i32();
            auto it = g2m.find(gtype);
            if (it == g2m.end())
                throw ReadError("Gmsh element type " + std::to_string(gtype) +
                                " not supported by the C++ reader");
            std::size_t n = static_cast<std::size_t>(nnpc.at(it->second));
            std::vector<std::int64_t> tags(num_tags), nodes(n);
            for (std::int32_t k = 0; k < nelem; ++k) {
                rCur.read_i32();  // element id
                for (std::int32_t j = 0; j < num_tags; ++j)
                    tags[j] = rCur.read_i32();
                for (std::size_t j = 0; j < n; ++j)
                    nodes[j] = rCur.read_i32();
                append_element(rBlocks, it->second, n, num_tags, tags.data(), nodes.data());
            }
            done += nelem;
        }
    }
    rCur.skip_to_end("Elements");
}

// NodeData / ElementData
/**
 * @param rOpts when the section's name is not wanted, the value block is
 *        skipped wholesale via `skip_to_end` -- no allocation and no parsing of
 *        the `nitems * ncomp` values. The name sits at the top of the section,
 *        before the values, so this costs nothing to decide.
 */
void read_data(GmshCursor& rCur, const std::string& rTag, bool is_ascii,
               std::unordered_map<std::string, NDArray>& rOut, const ReadOptions& rOpts) {
    std::int64_t num_str = std::stoll(gmsh_trim(rCur.read_line()));
    std::string name;
    for (std::int64_t i = 0; i < num_str; ++i) {
        std::string s = gmsh_trim(rCur.read_line());
        if (i == 0) {
            // strip quotes
            std::size_t q1 = s.find('"'), q2 = s.rfind('"');
            name = (q1 != std::string::npos && q2 > q1) ? s.substr(q1 + 1, q2 - q1 - 1) : s;
        }
    }
    std::int64_t num_real = std::stoll(gmsh_trim(rCur.read_line()));
    for (std::int64_t i = 0; i < num_real; ++i)
        rCur.read_line();
    std::int64_t num_int = std::stoll(gmsh_trim(rCur.read_line()));
    std::vector<std::int64_t> itags(num_int);
    for (std::int64_t i = 0; i < num_int; ++i)
        itags[i] = std::stoll(gmsh_trim(rCur.read_line()));
    std::size_t ncomp = static_cast<std::size_t>(itags[1]);
    std::size_t nitems = static_cast<std::size_t>(itags[2]);

    if (!rOpts.WantsAnyData() || !rOpts.WantsArray(name)) {
        rCur.skip_to_end(rTag);  // never touch the nitems * ncomp values
        return;
    }

    NDArray data(DType::Float64, {nitems, ncomp});
    double* dp = data.As<double>();
    if (is_ascii) {
        for (std::size_t i = 0; i < nitems; ++i) {
            rCur.next_int();  // index
            for (std::size_t c = 0; c < ncomp; ++c)
                dp[i * ncomp + c] = rCur.next_double();
        }
    } else {
        for (std::size_t i = 0; i < nitems; ++i) {
            rCur.read_i32();  // index
            for (std::size_t c = 0; c < ncomp; ++c)
                dp[i * ncomp + c] = rCur.read_f64();
        }
    }
    rCur.skip_to_end(rTag);
    if (ncomp == 1)
        data.Reshape({nitems});
    rOut.emplace(name, std::move(data));
}

NDArray slice_rows(const NDArray& rA, std::size_t r0, std::size_t r1) {
    std::size_t nc = rA.Shape().size() >= 2 ? rA.Shape()[1] : 1;
    std::size_t isz = dtype_size(rA.Dtype());
    std::vector<std::size_t> shape = rA.Shape();
    shape[0] = r1 - r0;
    NDArray out(rA.Dtype(), shape);
    if (r1 > r0)
        std::memcpy(out.Data(), rA.Data() + r0 * nc * isz, (r1 - r0) * nc * isz);
    return out;
}

// ---- version 4.1 -------------------------------------------------------------

struct E41 {
    std::string mType;
    std::size_t mN = 0;
    std::size_t mCount = 0;
    int mEntityDim = 0;
    int mEntityTag = 0;
    std::int32_t mPhysicalTag = 0;  // 0 == gmsh's "no physical group"
    std::vector<std::int32_t> mBounding;
    NDArray mConn;  // (count, n) Int64, 0-based gmsh node ids; moved into the
                    // cell block directly when the tag remap is the identity.
};

/**
 * @brief `$Entities` physical tags and bounding entities.
 *
 * Both arrays are indexed `[0..3]` by entity **dimension** -- which the section
 * never states, being implicit in its ordering (points, curves, surfaces,
 * volumes) -- and then keyed by entity tag.
 *
 * `$Entities` is the only place format 4.1 records physical-group membership:
 * an `$Elements` block names an `(entityDim, entityTag)` pair and the physical
 * tag lives on the entity, so without this section a 4.1 file has no
 * `gmsh:physical` and therefore no named regions.
 */
struct GmshEntities41 {
    std::array<std::unordered_map<std::int32_t, std::vector<std::int32_t>>, 4> mPhysical;
    std::array<std::unordered_map<std::int32_t, std::vector<std::int32_t>>, 4> mBounding;
    /// Whether any entity in the file carries a physical tag at all.
    bool mAnyPhysical = false;
};

/**
 * @brief Parse `$Entities` (format 4.1).
 *
 * Layout, per dimension `d`, after the four `size_t` per-dimension counts:
 * `tag(int32)`, a bounding box of **3** doubles for `d == 0` and **6**
 * otherwise (discarded -- the reader has nowhere to put it and gmsh recomputes
 * it), `numPhysicalTags(size_t)` then that many `int32`, and for `d > 0`
 * `numBounding(size_t)` then that many `int32`. Bounding-entity tags are
 * **signed**: the sign carries the boundary's orientation.
 *
 * Deliberately serial -- the section is tiny next to `$Nodes`/`$Elements` and
 * the work is hash-map inserts.
 */
GmshEntities41 read_entities_41(GmshCursor& rCur, bool is_ascii, int data_size) {
    auto rd_size = [&]() -> std::int64_t {
        return is_ascii ? rCur.next_int() : static_cast<std::int64_t>(rCur.read_uint(data_size));
    };
    auto rd_int = [&]() -> std::int32_t {
        return is_ascii ? static_cast<std::int32_t>(rCur.next_int()) : rCur.read_i32();
    };
    auto skip_dbl = [&](int n) {
        if (is_ascii) {
            // The ascii body is a whitespace token stream, so the coordinates
            // still have to be consumed even though they are thrown away.
            for (int i = 0; i < n; ++i)
                rCur.next_double();
        } else {
            rCur.mPos += static_cast<std::size_t>(n) * 8;
        }
    };

    GmshEntities41 out;
    std::array<std::int64_t, 4> counts{};
    for (int d = 0; d < 4; ++d)
        counts[static_cast<std::size_t>(d)] = rd_size();

    for (int d = 0; d < 4; ++d) {
        const std::size_t dz = static_cast<std::size_t>(d);
        for (std::int64_t i = 0; i < counts[dz]; ++i) {
            const std::int32_t tag = rd_int();
            skip_dbl(d == 0 ? 3 : 6);  // bounding box
            const std::int64_t num_phys = rd_size();
            if (num_phys > 0) {
                std::vector<std::int32_t> phys(static_cast<std::size_t>(num_phys));
                for (std::int64_t k = 0; k < num_phys; ++k)
                    phys[static_cast<std::size_t>(k)] = rd_int();
                out.mAnyPhysical = true;
                out.mPhysical[dz].emplace(tag, std::move(phys));
            }
            if (d > 0) {
                const std::int64_t num_bnd = rd_size();
                std::vector<std::int32_t> bnd(static_cast<std::size_t>(num_bnd));
                for (std::int64_t k = 0; k < num_bnd; ++k)
                    bnd[static_cast<std::size_t>(k)] = rd_int();
                out.mBounding[dz].emplace(tag, std::move(bnd));
            }
        }
    }
    rCur.skip_to_end("Entities");
    return out;
}

void read_nodes_41(GmshCursor& rCur, bool is_ascii, int data_size, NDArray& rPoints,
                   std::vector<std::int64_t>& rTags,
                   std::vector<std::array<std::int64_t, 2>>& rDimTags) {
    auto rd_size = [&]() -> std::int64_t {
        return is_ascii ? rCur.next_int() : static_cast<std::int64_t>(rCur.read_uint(data_size));
    };
    auto rd_int = [&]() -> int {
        return is_ascii ? static_cast<int>(rCur.next_int()) : rCur.read_i32();
    };
    auto rd_dbl = [&]() -> double { return is_ascii ? rCur.next_double() : rCur.read_f64(); };

    std::int64_t num_blocks = rd_size();
    std::int64_t num_nodes = rd_size();
    rd_size();  // min tag
    rd_size();  // max tag
    rPoints = NDArray(DType::Float64, {static_cast<std::size_t>(num_nodes), 3});
    rTags.resize(num_nodes);
    rDimTags.resize(num_nodes);
    double* pp = rPoints.As<double>();

    std::size_t idx = 0;
    for (std::int64_t b = 0; b < num_blocks; ++b) {
        int dim = rd_int();
        int entity_tag = rd_int();
        int parametric = rd_int();
        if (parametric != 0)
            throw ReadError("parametric Gmsh nodes not supported");
        std::int64_t nb = rd_size();
        const std::size_t nbz = static_cast<std::size_t>(nb);
        if (!is_ascii && data_size == 8) {
            // Native-endian, contiguous: bulk-copy tags (u64) and coords (3*f64).
            std::memcpy(&rTags[idx], rCur.mBuf.data() + rCur.mPos, nbz * 8);
            rCur.mPos += nbz * 8;
            for (std::size_t i = 0; i < nbz; ++i)
                rTags[idx + i] -= 1;
            std::memcpy(pp + idx * 3, rCur.mBuf.data() + rCur.mPos, nbz * 3 * 8);
            rCur.mPos += nbz * 3 * 8;
        } else {
            for (std::int64_t i = 0; i < nb; ++i)
                rTags[idx + i] = rd_size() - 1;
            for (std::int64_t i = 0; i < nb; ++i) {
                pp[(idx + i) * 3 + 0] = rd_dbl();
                pp[(idx + i) * 3 + 1] = rd_dbl();
                pp[(idx + i) * 3 + 2] = rd_dbl();
            }
        }
        for (std::int64_t i = 0; i < nb; ++i)
            rDimTags[idx + i] = {dim, entity_tag};
        idx += static_cast<std::size_t>(nb);
    }
    rCur.skip_to_end("Nodes");
}

/**
 * @brief Parse `$Elements` (format 4.1).
 *
 * @param pEntities the file's `$Entities`, or `nullptr` when it had none. Real
 *        gmsh files always emit `$Entities` first, so resolving each block's
 *        physical tag inline is safe -- the same section ordering the Python
 *        reference assumes.
 */
void read_elements_41(GmshCursor& rCur, bool is_ascii, int data_size, std::vector<E41>& rBlocks,
                      const GmshEntities41* pEntities) {
    auto rd_size = [&]() -> std::int64_t {
        return is_ascii ? rCur.next_int() : static_cast<std::int64_t>(rCur.read_uint(data_size));
    };
    auto rd_int = [&]() -> int {
        return is_ascii ? static_cast<int>(rCur.next_int()) : rCur.read_i32();
    };

    std::int64_t num_blocks = rd_size();
    rd_size();  // num elements
    rd_size();  // min tag
    rd_size();  // max tag
    const auto& g2m = gmsh_to_meshio_type();
    const auto& nnpc = num_nodes_per_cell();

    for (std::int64_t b = 0; b < num_blocks; ++b) {
        int entity_dim = rd_int();
        int entity_tag = rd_int();
        int etype = rd_int();
        std::int64_t num_ele = rd_size();
        auto it = g2m.find(etype);
        if (it == g2m.end())
            throw ReadError("Gmsh element type " + std::to_string(etype) +
                            " not supported by the C++ reader");
        std::size_t n = static_cast<std::size_t>(nnpc.at(it->second));
        E41 blk;
        blk.mType = it->second;
        blk.mN = n;
        blk.mCount = static_cast<std::size_t>(num_ele);
        blk.mEntityDim = entity_dim;
        blk.mEntityTag = entity_tag;
        if (pEntities && entity_dim >= 0 && entity_dim < 4) {
            const std::size_t dz = static_cast<std::size_t>(entity_dim);
            auto pit = pEntities->mPhysical[dz].find(entity_tag);
            // Only the first physical tag: a gmsh entity may carry several, but
            // one cell_data column can hold one. The Python reference makes the
            // same choice.
            if (pit != pEntities->mPhysical[dz].end() && !pit->second.empty())
                blk.mPhysicalTag = pit->second[0];
            if (entity_dim > 0) {
                auto bit = pEntities->mBounding[dz].find(entity_tag);
                if (bit != pEntities->mBounding[dz].end())
                    blk.mBounding = bit->second;
            }
        }
        const std::size_t nez = static_cast<std::size_t>(num_ele);
        blk.mConn = NDArray(DType::Int64, {nez, n});
        std::int64_t* dst = blk.mConn.As<std::int64_t>();
        if (!is_ascii && data_size == 8) {
            // Each element is [tag, node0..node(n-1)] u64, native-endian and
            // contiguous. Decode the nodes straight from the slurped buffer into
            // the owning connectivity array (drop the tag), one parallel pass.
            const std::size_t stride = n + 1;
            const char* base = rCur.mBuf.data() + rCur.mPos;
            parallel_for_bw(nez, [&](std::size_t e) {
                const char* row = base + (e * stride + 1) * 8;  // skip element tag
                for (std::size_t j = 0; j < n; ++j) {
                    std::uint64_t v;
                    std::memcpy(&v, row + j * 8, 8);
                    dst[e * n + j] = static_cast<std::int64_t>(v) - 1;
                }
            });
            rCur.mPos += nez * stride * 8;
        } else {
            std::size_t p = 0;
            for (std::int64_t e = 0; e < num_ele; ++e) {
                rd_size();  // element tag
                for (std::size_t j = 0; j < n; ++j)
                    dst[p++] = rd_size() - 1;
            }
        }
        rBlocks.push_back(std::move(blk));
    }
    rCur.skip_to_end("Elements");
}

Mesh read_gmsh41_body(GmshCursor& rCur, bool is_ascii, int data_size, const ReadOptions& rOpts,
                      GmshInfo* pInfo) {
    NDArray points(DType::Float64, {0, 3});
    std::vector<std::int64_t> point_tags;
    std::vector<std::array<std::int64_t, 2>> dim_tags;
    std::vector<E41> eblocks;
    GmshEntities41 entities;
    bool have_entities = false;
    std::unordered_map<std::string, NDArray> field_data, point_data, cell_data_raw;

    while (!rCur.eof()) {
        std::string line = rCur.next_nonblank();
        if (line.empty())
            break;
        if (line[0] != '$')
            throw ReadError("Gmsh: unexpected line " + line);
        std::string env = gmsh_trim(line.substr(1));
        if (env == "PhysicalNames")
            read_physical_names(rCur, field_data);
        else if (env == "Entities") {
            entities = read_entities_41(rCur, is_ascii, data_size);
            have_entities = true;
        } else if (env == "Nodes")
            read_nodes_41(rCur, is_ascii, data_size, points, point_tags, dim_tags);
        else if (env == "Elements")
            read_elements_41(rCur, is_ascii, data_size, eblocks,
                             have_entities ? &entities : nullptr);
        else if (env == "Periodic")
            throw ReadError("Gmsh $Periodic not supported by the C++ reader");
        else if (env == "NodeData")
            read_data(rCur, "NodeData", is_ascii, point_data, rOpts);
        else if (env == "ElementData")
            read_data(rCur, "ElementData", is_ascii, cell_data_raw, rOpts);
        else
            rCur.skip_to_end(env);
    }

    // When node tags are contiguous 0..N-1 (the common case) the tag->row remap
    // is the identity, so we can skip building it *and* skip the random-access
    // gather below (the connectivity is already the final mesh indexing).
    bool remap_identity = true;
    for (std::size_t i = 0; i < point_tags.size(); ++i)
        if (point_tags[i] != static_cast<std::int64_t>(i)) {
            remap_identity = false;
            break;
        }
    std::vector<std::int64_t> remap;
    if (!remap_identity) {
        std::int64_t max_tag = 0;
        for (auto t : point_tags)
            max_tag = std::max(max_tag, t);
        remap.assign(static_cast<std::size_t>(max_tag) + 1, -1);
        // Scatter: node tags are unique, so writes never alias -> parallel.
        parallel_for_bw(point_tags.size(), [&](std::size_t i) {
            remap[static_cast<std::size_t>(point_tags[i])] = static_cast<std::int64_t>(i);
        });
    }

    Mesh mesh;
    mesh.AssignPoints(std::move(points));
    for (auto& kv : point_data)
        mesh.AddPointData(kv.first, std::move(kv.second));
    for (auto& kv : field_data)
        mesh.AddFieldData(kv.first, std::move(kv.second));

    // Node entity (dim, tag) -> gmsh:dim_tags point data. Structural rather
    // than user data, but it lands in mesh.point_data all the same, so it obeys
    // the same filter: "points_only" must mean no point_data, with no
    // exceptions the caller has to know about. (A points_only mesh is an
    // explicitly lossy request; preserving this for round-tripping would make
    // the contract inconsistent instead of useful.)
    if (rOpts.WantsAnyData() && rOpts.WantsArray("gmsh:dim_tags")) {
        NDArray dt(DType::Int64, {dim_tags.size(), 2});
        parallel_for_bw(dim_tags.size(), [&](std::size_t i) {
            dt.As<std::int64_t>()[i * 2 + 0] = dim_tags[i][0];
            dt.As<std::int64_t>()[i * 2 + 1] = dim_tags[i][1];
        });
        mesh.AddPointData("gmsh:dim_tags", std::move(dt));
    }

    std::vector<NDArray> geom_blocks, physical_blocks;
    for (auto& b : eblocks) {
        const std::vector<int>& perm = gmsh_to_meshio_perm(b.mType);
        const int* prm = perm.empty() ? nullptr : perm.data();
        if (remap_identity && !prm) {
            // Identity remap, no reorder -> the connectivity is already final:
            // move the owning (count, n) array straight into the cell block.
            mesh.AddCellBlock(b.mType, std::move(b.mConn));
        } else {
            NDArray data(DType::Int64, {b.mCount, b.mN});
            std::int64_t* dp = data.As<std::int64_t>();
            const std::int64_t* cn = b.mConn.As<std::int64_t>();
            if (remap_identity) {
                parallel_for_bw(b.mCount, [&](std::size_t r) {
                    for (std::size_t j = 0; j < b.mN; ++j)
                        dp[r * b.mN + j] = cn[r * b.mN + static_cast<std::size_t>(prm[j])];
                });
            } else {
                // Gather through the prebuilt read-only remap -> parallel by row.
                parallel_for_bw(b.mCount, [&](std::size_t r) {
                    for (std::size_t j = 0; j < b.mN; ++j) {
                        std::size_t src = prm ? static_cast<std::size_t>(prm[j]) : j;
                        dp[r * b.mN + j] = remap[static_cast<std::size_t>(cn[r * b.mN + src])];
                    }
                });
            }
            mesh.AddCellBlock(b.mType, std::move(data));
        }

        NDArray ge(DType::Int32, {b.mCount});
        std::int32_t* gep = ge.As<std::int32_t>();
        const std::int32_t etag = b.mEntityTag;
        parallel_for_bw(b.mCount, [&](std::size_t r) { gep[r] = etag; });
        geom_blocks.push_back(std::move(ge));

        // One array per block whenever the file tags anything at all, with 0
        // (gmsh's "no physical group") for the untagged blocks. Emitting only
        // the tagged ones -- what the Python reference does -- would leave
        // gmsh:physical shorter than mesh.cells, which the uniform mesh API
        // cannot express and gmsh_attach_regions rejects.
        if (entities.mAnyPhysical) {
            NDArray ph(DType::Int32, {b.mCount});
            std::int32_t* php = ph.As<std::int32_t>();
            const std::int32_t ptag = b.mPhysicalTag;
            parallel_for_bw(b.mCount, [&](std::size_t r) { php[r] = ptag; });
            physical_blocks.push_back(std::move(ph));
        }
    }

    for (auto& kv : cell_data_raw) {
        std::vector<NDArray> per_block;
        std::size_t offset = 0;
        for (const auto& b : eblocks) {
            per_block.push_back(slice_rows(kv.second, offset, offset + b.mCount));
            offset += b.mCount;
        }
        mesh.AddCellData(kv.first, std::move(per_block));
    }
    if (!physical_blocks.empty() && rOpts.WantsAnyData() && rOpts.WantsArray("gmsh:physical"))
        mesh.AddCellData("gmsh:physical", std::move(physical_blocks));
    if (!geom_blocks.empty() && rOpts.WantsAnyData() && rOpts.WantsArray("gmsh:geometrical"))
        mesh.AddCellData("gmsh:geometrical", std::move(geom_blocks));

    // Bounding entities are signed (the sign carries orientation), so they
    // cannot be a Region -- they ride the side channel instead, like MedInfo.
    if (pInfo && have_entities) {
        pInfo->mBoundingEntities.reserve(eblocks.size());
        for (const auto& b : eblocks)
            pInfo->mBoundingEntities.push_back(b.mBounding);
    }

    gmsh_attach_regions(mesh);

    return mesh;
}

/**
 * @brief Walk a 4.1 `$Nodes`/`$Elements` section far enough to describe it,
 *        skipping the bulk payload.
 *
 * Scoped to format 4.1 deliberately. 4.1 groups elements into typed blocks
 * whose headers carry the type and count, so a summary is nearly free. Format
 * 2.2's `$Elements` is a flat list where **every element carries its own type**,
 * so summarizing it would cost essentially a full parse -- there is no cheap
 * path to have, and pretending otherwise would just move the work around.
 * 2.2 therefore keeps the honest full-read fallback.
 */
struct GmshMeta41 {
    std::size_t mNumPoints = 0;
    std::vector<CellBlockInfo> mBlocks;
    /// Per block, its `$Elements` `(entityDim, entityTag)` -- what resolves the
    /// block's physical group against `$Entities`.
    std::vector<std::pair<int, std::int32_t>> mBlockEntities;
    std::vector<std::string> mPointDataNames;
    std::vector<std::string> mCellDataNames;
    bool mHasDimTags = false;
};

/// Consume the remainder of the current line.
void gmsh_finish_line(GmshCursor& rCur) {
    while (!rCur.eof() && rCur.mBuf[rCur.mPos] != '\n')
        ++rCur.mPos;
    if (!rCur.eof())
        ++rCur.mPos;
}

void gmsh_scan_nodes_41(GmshCursor& rCur, bool is_ascii, int data_size, GmshMeta41& rMeta) {
    auto rd_size = [&]() -> std::int64_t {
        return is_ascii ? rCur.next_int() : static_cast<std::int64_t>(rCur.read_uint(data_size));
    };
    const std::int64_t num_blocks = rd_size();
    const std::int64_t num_nodes = rd_size();
    rd_size();  // min tag
    rd_size();  // max tag
    rMeta.mNumPoints = static_cast<std::size_t>(num_nodes < 0 ? 0 : num_nodes);
    // Every 4.1 file yields gmsh:dim_tags, built from the per-node entity.
    rMeta.mHasDimTags = num_nodes > 0;

    if (is_ascii) {
        // Coordinates are the bulk here and are never needed for a summary, so
        // skip to the terminator rather than tokenizing them.
        gmsh_finish_line(rCur);
        rCur.skip_to_end("Nodes");
        return;
    }
    // Binary: advance by the exact payload size instead of scanning, since raw
    // bytes could otherwise contain anything that looks like a terminator.
    for (std::int64_t b = 0; b < num_blocks; ++b) {
        rCur.read_i32();  // dim
        rCur.read_i32();  // entity tag
        rCur.read_i32();  // parametric
        const std::int64_t in_block = static_cast<std::int64_t>(rCur.read_uint(data_size));
        rCur.mPos += static_cast<std::size_t>(in_block) * static_cast<std::size_t>(data_size);
        rCur.mPos += static_cast<std::size_t>(in_block) * 3u * 8u;
    }
    rCur.skip_to_end("Nodes");
}

void gmsh_scan_elements_41(GmshCursor& rCur, bool is_ascii, int data_size, GmshMeta41& rMeta) {
    auto rd_size = [&]() -> std::int64_t {
        return is_ascii ? rCur.next_int() : static_cast<std::int64_t>(rCur.read_uint(data_size));
    };
    const std::int64_t num_blocks = rd_size();
    rd_size();  // num elements
    rd_size();  // min tag
    rd_size();  // max tag
    const auto& g2m = gmsh_to_meshio_type();
    const auto& nnpc = num_nodes_per_cell();
    if (is_ascii)
        gmsh_finish_line(rCur);

    for (std::int64_t b = 0; b < num_blocks; ++b) {
        int etype = 0, entity_dim = 0;
        std::int32_t entity_tag = 0;
        std::int64_t num_ele = 0;
        if (is_ascii) {
            entity_dim = static_cast<int>(rCur.next_int());
            entity_tag = static_cast<std::int32_t>(rCur.next_int());
            etype = static_cast<int>(rCur.next_int());
            num_ele = rCur.next_int();
            gmsh_finish_line(rCur);
        } else {
            entity_dim = rCur.read_i32();
            entity_tag = rCur.read_i32();
            etype = rCur.read_i32();
            num_ele = static_cast<std::int64_t>(rCur.read_uint(data_size));
        }

        auto it = g2m.find(etype);
        if (it == g2m.end())
            throw ReadError("Gmsh element type " + std::to_string(etype) +
                            " not supported by the C++ reader");
        const std::size_t n = static_cast<std::size_t>(nnpc.at(it->second));

        CellBlockInfo info;
        info.mType = it->second;
        info.mNumCells = static_cast<std::size_t>(num_ele < 0 ? 0 : num_ele);
        info.mNodesPerCell = n;
        rMeta.mBlocks.push_back(std::move(info));
        rMeta.mBlockEntities.emplace_back(entity_dim, entity_tag);

        // Skip the block's connectivity: one line per element (ascii) or
        // num_ele * (tag + n nodes) fixed-width integers (binary).
        if (is_ascii) {
            for (std::int64_t e = 0; e < num_ele; ++e)
                gmsh_finish_line(rCur);
        } else {
            rCur.mPos +=
                static_cast<std::size_t>(num_ele) * (n + 1u) * static_cast<std::size_t>(data_size);
        }
    }
    rCur.skip_to_end("Elements");
}

/// Read a `$NodeData`/`$ElementData` section's name, then skip its values.
std::string gmsh_scan_data_name(GmshCursor& rCur, const std::string& rTag) {
    const std::int64_t num_str = std::stoll(gmsh_trim(rCur.read_line()));
    std::string name;
    for (std::int64_t i = 0; i < num_str; ++i) {
        const std::string line = gmsh_trim(rCur.read_line());
        if (i == 0) {
            const std::size_t q1 = line.find('"'), q2 = line.rfind('"');
            name = (q1 != std::string::npos && q2 > q1) ? line.substr(q1 + 1, q2 - q1 - 1) : line;
        }
    }
    rCur.skip_to_end(rTag);
    return name;
}
}  // namespace

Mesh read_gmsh(const std::string& rPath, const ReadOptions& rOpts) {
    GmshInfo unused;
    return read_gmsh(rPath, unused, rOpts);
}

Mesh read_gmsh(const std::string& rPath, GmshInfo& rInfo, const ReadOptions& rOpts) {
    // Memory-mapped where that pays (see detail/file_source.hpp), copied
    // otherwise. The source is function-local: every parsed value is copied
    // into owning mesh storage below, so nothing in the returned Mesh points
    // back into this buffer.
    const detail::FileSource source(rPath, rOpts.mMmap);
    const std::string_view buf = source.View();
    GmshCursor cur(buf);

    if (gmsh_trim(cur.read_line()) != "$MeshFormat")
        throw ReadError("Expected $MeshFormat");
    std::string fmt = cur.read_line();
    std::istringstream fss(fmt);
    std::string version;
    int file_type = 0, data_size = 8;
    fss >> version >> file_type >> data_size;
    bool is_ascii = (file_type == 0);
    if (!is_ascii) {
        cur.read_i32();  // endianness marker
        // consume trailing newline before $EndMeshFormat
        if (cur.mPos < buf.size() && buf[cur.mPos] == '\n')
            ++cur.mPos;
    }
    cur.skip_to_end("MeshFormat");

    if (version == "4.1" || version == "4")
        return read_gmsh41_body(cur, is_ascii, data_size, rOpts, &rInfo);
    if (version.rfind("2", 0) != 0)
        throw ReadError("C++ Gmsh reader handles versions 2.2 and 4.1 only");

    NDArray points(DType::Float64, {0, 3});
    std::vector<std::int64_t> point_tags;
    std::vector<EBlock> eblocks;
    std::unordered_map<std::string, NDArray> field_data, point_data, cell_data_raw;

    while (!cur.eof()) {
        std::string line = cur.next_nonblank();
        if (line.empty())
            break;
        if (line[0] != '$')
            throw ReadError("Gmsh: unexpected line " + line);
        std::string env = gmsh_trim(line.substr(1));
        if (env == "PhysicalNames")
            read_physical_names(cur, field_data);
        else if (env == "Nodes")
            read_nodes(cur, is_ascii, points, point_tags);
        else if (env == "Elements")
            read_elements(cur, is_ascii, eblocks);
        else if (env == "Periodic")
            throw ReadError("Gmsh $Periodic not supported by the C++ reader");
        else if (env == "NodeData")
            read_data(cur, "NodeData", is_ascii, point_data, rOpts);
        else if (env == "ElementData")
            read_data(cur, "ElementData", is_ascii, cell_data_raw, rOpts);
        else
            cur.skip_to_end(env);
    }

    // Build node-tag remap (gmsh ids are 1-based, possibly non-contiguous).
    std::int64_t max_tag = 0;
    for (auto t : point_tags)
        max_tag = std::max(max_tag, t - 1);
    std::vector<std::int64_t> remap(static_cast<std::size_t>(max_tag) + 1, -1);
    // Scatter: node tags are unique, so writes never alias -> parallel.
    parallel_for_bw(point_tags.size(), [&](std::size_t i) {
        remap[static_cast<std::size_t>(point_tags[i] - 1)] = static_cast<std::int64_t>(i);
    });

    Mesh mesh;
    mesh.AssignPoints(std::move(points));
    for (auto& kv : point_data)
        mesh.AddPointData(kv.first, std::move(kv.second));
    for (auto& kv : field_data)
        mesh.AddFieldData(kv.first, std::move(kv.second));

    // Determine which tag columns are present across all blocks.
    std::size_t min_tags = eblocks.empty() ? 0 : SIZE_MAX;
    for (const auto& b : eblocks)
        min_tags = std::min(min_tags, b.mNumTags);

    std::vector<NDArray> physical_blocks, geometrical_blocks;
    for (const auto& b : eblocks) {
        const std::vector<int>& perm = gmsh_to_meshio_perm(b.mType);
        NDArray data(DType::Int64, {b.mCount, b.mN});
        std::int64_t* dp = data.As<std::int64_t>();
        // Gather through the prebuilt read-only remap -> parallel over rows.
        parallel_for_bw(b.mCount, [&](std::size_t r) {
            for (std::size_t j = 0; j < b.mN; ++j) {
                std::size_t src = perm.empty() ? j : static_cast<std::size_t>(perm[j]);
                std::int64_t gid = b.mConn[r * b.mN + src];
                dp[r * b.mN + j] = remap[static_cast<std::size_t>(gid)];
            }
        });
        mesh.AddCellBlock(b.mType, std::move(data));

        if (min_tags >= 1) {
            NDArray ph(DType::Int32, {b.mCount});
            std::int32_t* php = ph.As<std::int32_t>();
            parallel_for_bw(b.mCount, [&](std::size_t r) {
                php[r] = static_cast<std::int32_t>(b.mTags[r * b.mNumTags + 0]);
            });
            physical_blocks.push_back(std::move(ph));
        }
        if (min_tags >= 2) {
            NDArray ge(DType::Int32, {b.mCount});
            std::int32_t* gep = ge.As<std::int32_t>();
            parallel_for_bw(b.mCount, [&](std::size_t r) {
                gep[r] = static_cast<std::int32_t>(b.mTags[r * b.mNumTags + 1]);
            });
            geometrical_blocks.push_back(std::move(ge));
        }
    }

    // Split ElementData (concatenated over blocks) back per block.
    for (auto& kv : cell_data_raw) {
        std::vector<NDArray> per_block;
        std::size_t offset = 0;
        for (const auto& b : eblocks) {
            per_block.push_back(slice_rows(kv.second, offset, offset + b.mCount));
            offset += b.mCount;
        }
        mesh.AddCellData(kv.first, std::move(per_block));
    }
    if (!physical_blocks.empty() && rOpts.WantsAnyData() && rOpts.WantsArray("gmsh:physical"))
        mesh.AddCellData("gmsh:physical", std::move(physical_blocks));
    if (!geometrical_blocks.empty() && rOpts.WantsAnyData() && rOpts.WantsArray("gmsh:geometrical"))
        mesh.AddCellData("gmsh:geometrical", std::move(geometrical_blocks));

    gmsh_attach_regions(mesh);

    return mesh;
}

// ---- writer ------------------------------------------------------------------

namespace {

void write_physical_names(std::ostream& rOs, const Mesh& rMesh) {
    // field_data first, then any named region describing a group field_data
    // does not — so a mesh carrying gmsh's own metadata writes byte-identical
    // bytes, and one whose groups came from another format still gets them.
    std::vector<std::tuple<long long, long long, std::string>> sortable =
        gmsh_physical_rows(rMesh);  // dim, num, name
    if (sortable.empty())
        return;
    rOs << "$PhysicalNames\n" << sortable.size() << "\n";
    for (auto& e : sortable)
        rOs << std::get<0>(e) << ' ' << std::get<1>(e) << " \"" << std::get<2>(e) << "\"\n";
    rOs << "$EndPhysicalNames\n";
}

// Writes the cell-data array named `rName` as one $ElementData-style section,
// concatenated across cell blocks.
void write_data(std::ostream& rOs, const char* pTag, const std::string& rName, const Mesh& rMesh,
                bool binary) {
    // Concatenate blocks.
    const std::size_t nblocks = rMesh.CellDataNumBlocks(rName);
    std::size_t total = 0, ncomp = 1;
    for (std::size_t k = 0; k < nblocks; ++k) {
        const NDArray& b = rMesh.CellData(rName, k);
        total += b.Shape().empty() ? 0 : b.Shape()[0];
        ncomp = b.Shape().size() >= 2 ? b.Shape()[1] : 1;
    }
    rOs << "$" << pTag << "\n1\n\"" << rName << "\"\n1\n0\n3\n0\n"
        << ncomp << "\n"
        << total << "\n";
    std::int64_t idx = 1;
    for (std::size_t k = 0; k < nblocks; ++k) {
        const NDArray& b = rMesh.CellData(rName, k);
        std::size_t rows = b.Shape().empty() ? 0 : b.Shape()[0];
        for (std::size_t r = 0; r < rows; ++r) {
            if (binary) {
                std::int32_t id = static_cast<std::int32_t>(idx);
                rOs.write(reinterpret_cast<const char*>(&id), 4);
                for (std::size_t c = 0; c < ncomp; ++c) {
                    double v = detail::read_double(b, r * ncomp + c);
                    rOs.write(reinterpret_cast<const char*>(&v), 8);
                }
            } else {
                rOs << idx;
                char buf[32];
                for (std::size_t c = 0; c < ncomp; ++c) {
                    std::snprintf(buf, sizeof(buf), " %.17g",
                                  detail::read_double(b, r * ncomp + c));
                    rOs << buf;
                }
                rOs << '\n';
            }
            ++idx;
        }
    }
    if (binary)
        rOs << '\n';
    rOs << "$End" << pTag << "\n";
}

/// One geometric entity on the 4.1 write path: the nodes it owns and, when one
/// exists, the cell block that supplies its physical tag and bounding entities.
struct GmshWriteEntity41 {
    int mDim = 0;
    std::int32_t mTag = 0;
    std::vector<std::int64_t> mNodes;  // ascending point indices
    std::size_t mCellBlock = static_cast<std::size_t>(-1);
};

/**
 * @brief Partition the points into 4.1 entity blocks and pair each with its
 *        cell block.
 *
 * The entity set is the **union** of the node entities (from
 * `point_data["gmsh:dim_tags"]`) and the cell-block entities (from
 * `cell_data["gmsh:geometrical"]`). Node entities are nearly a superset — gmsh
 * assigns one to every object needed to describe the geometry — but not quite:
 * a straight curve whose only nodes are its two endpoints has its nodes on the
 * *points*, so it owns none of its own while still carrying line elements
 * (`example/example.msh` has six such curves). Taking the union is what keeps
 * their physical tags and bounding entities from being dropped; gmsh is happy
 * with an entity that has no `$Nodes` block.
 *
 * Sorting the unique `(dim, tag)` pairs lexicographically is what makes the
 * per-dimension counts and the points/curves/surfaces/volumes emission order
 * agree by construction, so `$Entities` and `$Nodes` cannot disagree.
 *
 * @return one entry per entity, or an empty vector when the mesh carries no
 *         `gmsh:dim_tags` (the caller then writes its single legacy block and
 *         no `$Entities` at all).
 * @throws WriteError on an entity dimension above 3, or on two cell blocks
 *         claiming the same `(dim, entity tag)` — the entity would then have
 *         no single physical tag to write.
 */
std::vector<GmshWriteEntity41> gmsh_entity_blocks_41(
    const Mesh& rMesh, const std::function<int(const std::string&)>& rCellDim) {
    std::vector<GmshWriteEntity41> out;
    if (!rMesh.HasPointData("gmsh:dim_tags"))
        return out;
    const NDArray& dt = rMesh.PointData("gmsh:dim_tags");
    const std::size_t n = rMesh.NumPoints();
    if (n == 0)
        return out;
    // Normally an (N, 2) array, but the flat bindings (WASM, C API) carry data
    // arrays without a shape, so the same values arrive as a bare 2N-long row.
    // Both describe one (dim, tag) pair per point; refusing the flat one would
    // silently drop $Entities on exactly the surfaces this exists to serve.
    const std::size_t stride = dt.Shape().size() >= 2 ? dt.Shape()[1] : dt.Size() / n;
    if (stride < 2 || dt.Size() < n * stride)
        return out;

    // Per-block entity keys, needed both for the union below and for the
    // pairing pass afterwards.
    const bool has_geometrical = rMesh.HasCellData("gmsh:geometrical");
    std::vector<std::pair<int, std::int32_t>> cell_keys(rMesh.NumCellBlocks());
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci) {
        std::int32_t etag = 0;
        if (has_geometrical && ci < rMesh.CellDataNumBlocks("gmsh:geometrical") &&
            rMesh.CellData("gmsh:geometrical", ci).Size() > 0)
            etag = static_cast<std::int32_t>(
                detail::read_int(rMesh.CellData("gmsh:geometrical", ci), 0));
        cell_keys[ci] = {rCellDim(std::string(rMesh.Cells(ci).Type())), etag};
    }

    // Unique (dim, tag) pairs, sorted. Serial: the entity count is tiny next to
    // the node count, and a deterministic order is the whole point.
    std::vector<std::pair<int, std::int32_t>> keys(n);
    for (std::size_t i = 0; i < n; ++i)
        keys[i] = {static_cast<int>(detail::read_int(dt, i * stride + 0)),
                   static_cast<std::int32_t>(detail::read_int(dt, i * stride + 1))};
    std::vector<std::pair<int, std::int32_t>> uniq = keys;
    uniq.insert(uniq.end(), cell_keys.begin(), cell_keys.end());
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    for (const auto& k : uniq) {
        if (k.first < 0 || k.first > 3)
            throw WriteError("Gmsh 4.1 writer: entity dimension " + std::to_string(k.first) +
                             " is outside 0..3");
    }

    std::map<std::pair<int, std::int32_t>, std::size_t> index;
    out.reserve(uniq.size());
    for (const auto& k : uniq) {
        index.emplace(k, out.size());
        GmshWriteEntity41 e;
        e.mDim = k.first;
        e.mTag = k.second;
        out.push_back(std::move(e));
    }
    // Ascending node ids per entity, which is what the reader's tag->row remap
    // expects to see and what keeps the output stable.
    for (std::size_t i = 0; i < n; ++i)
        out[index.at(keys[i])].mNodes.push_back(static_cast<std::int64_t>(i));

    // Pair each entity with the cell block living on it.
    for (std::size_t ci = 0; ci < cell_keys.size(); ++ci) {
        const std::size_t ei = index.at(cell_keys[ci]);
        if (out[ei].mCellBlock != static_cast<std::size_t>(-1))
            throw WriteError("Gmsh 4.1 writer: two cell blocks share entity (dim " +
                             std::to_string(out[ei].mDim) + ", tag " +
                             std::to_string(out[ei].mTag) + ")");
        out[ei].mCellBlock = ci;
    }
    return out;
}

}  // namespace

void write_gmsh22(const std::string& rPath, const Mesh& rMesh, bool binary) {
    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t num_points = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    const std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 0;
    const std::size_t nblocks = rMesh.NumCellBlocks();

    // Tag cell data ("gmsh:physical"/"gmsh:geometrical") is written inline with
    // the elements; per-block zeros stand in when a tag column is absent.
    const bool has_physical = rMesh.HasCellData("gmsh:physical");
    const bool has_geometrical = rMesh.HasCellData("gmsh:geometrical");
    std::vector<NDArray> zeros_phys, zeros_geom;
    if (!has_physical) {
        // No gmsh:physical column of its own: synthesize one from any tagged
        // Cell regions, so a mesh whose groups came from another format still
        // writes real physical groups. With no such regions this yields the
        // per-block zeros it always did, and the output is byte-identical.
        zeros_phys = gmsh_tags_from_regions(rMesh);
        if (zeros_phys.empty())
            for (const auto cb : rMesh.CellRange())
                zeros_phys.emplace_back(DType::Int32, std::vector<std::size_t>{cb.NumCells()});
    }
    if (!has_geometrical)
        for (const auto cb : rMesh.CellRange())
            zeros_geom.emplace_back(DType::Int32, std::vector<std::size_t>{cb.NumCells()});

    os << "$MeshFormat\n2.2 " << (binary ? 1 : 0) << " 8\n";
    if (binary) {
        std::int32_t one = 1;
        os.write(reinterpret_cast<const char*>(&one), 4);
        os << '\n';
    }
    os << "$EndMeshFormat\n";

    write_physical_names(os, rMesh);

    // Nodes.
    os << "$Nodes\n" << num_points << "\n";
    if (binary) {
        for (std::size_t i = 0; i < num_points; ++i) {
            std::int32_t id = static_cast<std::int32_t>(i + 1);
            os.write(reinterpret_cast<const char*>(&id), 4);
            for (std::size_t c = 0; c < 3; ++c) {
                double v = (c < dim) ? detail::read_double(points, i * dim + c) : 0.0;
                os.write(reinterpret_cast<const char*>(&v), 8);
            }
        }
        os << '\n';
    } else {
        // %zu (up to 20 digits) + 3x %.16e (up to 24 chars each) + separators/'\n'/'\0'
        char buf[128];
        for (std::size_t i = 0; i < num_points; ++i) {
            double x = (0 < dim) ? detail::read_double(points, i * dim + 0) : 0.0;
            double y = (1 < dim) ? detail::read_double(points, i * dim + 1) : 0.0;
            double z = (2 < dim) ? detail::read_double(points, i * dim + 2) : 0.0;
            std::snprintf(buf, sizeof(buf), "%zu %.16e %.16e %.16e\n", i + 1, x, y, z);
            os << buf;
        }
    }
    os << "$EndNodes\n";

    // Elements.
    std::size_t total_cells = 0;
    for (const auto cb : rMesh.CellRange())
        total_cells += cb.NumCells();
    os << "$Elements\n" << total_cells << "\n";
    const auto& m2g = meshio_to_gmsh_type();
    std::size_t consecutive = 0;
    for (std::size_t k = 0; k < nblocks; ++k) {
        const auto cb = rMesh.Cells(k);
        auto it = m2g.find(cb.Type());
        if (it == m2g.end())
            throw WriteError("Gmsh writer: unsupported cell type " + cb.Type());
        int gtype = it->second;
        const NDArray& conn = cb.Conn();
        std::size_t n = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        const std::vector<int>& perm = meshio_to_gmsh_perm(cb.Type());
        std::size_t count = cb.NumCells();
        const NDArray& ph = has_physical ? rMesh.CellData("gmsh:physical", k) : zeros_phys[k];
        const NDArray& ge = has_geometrical ? rMesh.CellData("gmsh:geometrical", k) : zeros_geom[k];

        if (binary) {
            std::int32_t hdr[3] = {gtype, static_cast<std::int32_t>(count), 2};
            os.write(reinterpret_cast<const char*>(hdr), 12);
            for (std::size_t r = 0; r < count; ++r) {
                std::int32_t id = static_cast<std::int32_t>(consecutive + r + 1);
                std::int32_t t0 = static_cast<std::int32_t>(detail::read_int(ph, r));
                std::int32_t t1 = static_cast<std::int32_t>(detail::read_int(ge, r));
                os.write(reinterpret_cast<const char*>(&id), 4);
                os.write(reinterpret_cast<const char*>(&t0), 4);
                os.write(reinterpret_cast<const char*>(&t1), 4);
                for (std::size_t j = 0; j < n; ++j) {
                    std::size_t src = perm.empty() ? j : static_cast<std::size_t>(perm[j]);
                    std::int32_t node =
                        static_cast<std::int32_t>(detail::read_int(conn, r * n + src) + 1);
                    os.write(reinterpret_cast<const char*>(&node), 4);
                }
            }
        } else {
            for (std::size_t r = 0; r < count; ++r) {
                os << (consecutive + r + 1) << ' ' << gtype << " 2 " << detail::read_int(ph, r)
                   << ' ' << detail::read_int(ge, r);
                for (std::size_t j = 0; j < n; ++j) {
                    std::size_t src = perm.empty() ? j : static_cast<std::size_t>(perm[j]);
                    os << ' ' << (detail::read_int(conn, r * n + src) + 1);
                }
                os << '\n';
            }
        }
        consecutive += count;
    }
    if (binary)
        os << '\n';
    os << "$EndElements\n";

    for (const auto& name : rMesh.PointDataNames()) {
        if (name == "gmsh:dim_tags")
            continue;
        // Reusing write_data (cell-data-shaped) for point data is awkward; inline:
        const NDArray& d = rMesh.PointData(name);
        std::size_t ncomp = d.Shape().size() >= 2 ? d.Shape()[1] : 1;
        std::size_t rows = d.Shape().empty() ? 0 : d.Shape()[0];
        os << "$NodeData\n1\n\"" << name << "\"\n1\n0\n3\n0\n" << ncomp << "\n" << rows << "\n";
        char buf[32];
        for (std::size_t r = 0; r < rows; ++r) {
            if (binary) {
                std::int32_t id = static_cast<std::int32_t>(r + 1);
                os.write(reinterpret_cast<const char*>(&id), 4);
                for (std::size_t c = 0; c < ncomp; ++c) {
                    double v = detail::read_double(d, r * ncomp + c);
                    os.write(reinterpret_cast<const char*>(&v), 8);
                }
            } else {
                os << (r + 1);
                for (std::size_t c = 0; c < ncomp; ++c) {
                    std::snprintf(buf, sizeof(buf), " %.17g",
                                  detail::read_double(d, r * ncomp + c));
                    os << buf;
                }
                os << '\n';
            }
        }
        if (binary)
            os << '\n';
        os << "$EndNodeData\n";
    }

    for (const auto& name : rMesh.CellDataNames()) {
        if (name == "gmsh:physical" || name == "gmsh:geometrical" || name == "cell_tags")
            continue;
        write_data(os, "ElementData", name, rMesh, binary);
    }
}

void write_gmsh41(const std::string& rPath, const Mesh& rMesh, bool binary) {
    write_gmsh41(rPath, rMesh, binary, GmshInfo{});
}

void write_gmsh41(const std::string& rPath, const Mesh& rMesh, bool binary, const GmshInfo& rInfo) {
    std::ofstream os(rPath, std::ios::binary);
    if (!os)
        throw WriteError("Could not open file for writing: " + rPath);

    const std::size_t num_points = rMesh.NumPoints();
    const NDArray& points = rMesh.Points();
    const std::size_t dim = points.Shape().size() >= 2 ? points.Shape()[1] : 0;
    const int data_size = 8;

    auto put_u64 = [&](std::uint64_t v) { os.write(reinterpret_cast<const char*>(&v), 8); };
    auto put_i32 = [&](std::int32_t v) { os.write(reinterpret_cast<const char*>(&v), 4); };
    auto put_f64 = [&](double v) { os.write(reinterpret_cast<const char*>(&v), 8); };

    // "gmsh:geometrical" supplies the per-block entity tag below; the other
    // tag names are excluded from the $NodeData/$ElementData sections.
    const bool has_geometrical = rMesh.HasCellData("gmsh:geometrical");

    const auto& topo = topological_dimension();
    auto cell_dim = [&](const std::string& t) -> int {
        auto it = topo.find(t);
        return it == topo.end() ? 0 : it->second;
    };

    os << "$MeshFormat\n4.1 " << (binary ? 1 : 0) << " " << data_size << "\n";
    if (binary) {
        put_i32(1);
        os << '\n';
    }
    os << "$EndMeshFormat\n";

    write_physical_names(os, rMesh);

    // The 3-padded coordinates of one point, formatted the one way both the
    // single-block and per-entity paths use.
    auto put_coords_ascii = [&](std::size_t i) {
        char buf[128];  // 3x %.16e (up to 24 chars each) + separators/'\n'/'\0'
        double x = (0 < dim) ? detail::read_double(points, i * dim + 0) : 0.0;
        double y = (1 < dim) ? detail::read_double(points, i * dim + 1) : 0.0;
        double z = (2 < dim) ? detail::read_double(points, i * dim + 2) : 0.0;
        std::snprintf(buf, sizeof(buf), "%.16e %.16e %.16e\n", x, y, z);
        os << buf;
    };

    const std::vector<GmshWriteEntity41> entities = gmsh_entity_blocks_41(rMesh, cell_dim);

    // $Entities, and with it the per-entity $Nodes blocks, only when the mesh
    // says which entity each node belongs to. Without that the legacy
    // single-block output below is emitted byte for byte as before.
    if (!entities.empty()) {
        const bool has_physical = rMesh.HasCellData("gmsh:physical");
        std::array<std::size_t, 4> per_dim{};
        for (const auto& e : entities)
            ++per_dim[static_cast<std::size_t>(e.mDim)];

        os << "$Entities\n";
        if (binary) {
            for (int d = 0; d < 4; ++d)
                put_u64(per_dim[static_cast<std::size_t>(d)]);
        } else {
            os << per_dim[0] << " " << per_dim[1] << " " << per_dim[2] << " " << per_dim[3] << "\n";
        }
        for (const auto& e : entities) {
            const bool paired = e.mCellBlock != static_cast<std::size_t>(-1);
            // The bounding box is written as zeros, matching the Python
            // reference writer: recovering a real one would mean a min/max pass
            // over each entity's points, and gmsh recomputes it on load anyway.
            std::int32_t ptag = 0;
            bool has_ptag = false;
            if (paired && has_physical && e.mCellBlock < rMesh.CellDataNumBlocks("gmsh:physical") &&
                rMesh.CellData("gmsh:physical", e.mCellBlock).Size() > 0) {
                ptag = static_cast<std::int32_t>(
                    detail::read_int(rMesh.CellData("gmsh:physical", e.mCellBlock), 0));
                has_ptag = true;
            }
            const std::vector<std::int32_t>* bounds = nullptr;
            if (paired && e.mDim > 0 && e.mCellBlock < rInfo.mBoundingEntities.size() &&
                !rInfo.mBoundingEntities[e.mCellBlock].empty())
                bounds = &rInfo.mBoundingEntities[e.mCellBlock];

            if (binary) {
                put_i32(e.mTag);
                for (int k = 0; k < (e.mDim == 0 ? 3 : 6); ++k)
                    put_f64(0.0);
                put_u64(has_ptag ? 1u : 0u);
                if (has_ptag)
                    put_i32(ptag);
                if (e.mDim > 0) {
                    put_u64(bounds ? bounds->size() : 0u);
                    if (bounds)
                        for (std::int32_t b : *bounds)
                            put_i32(b);
                }
            } else {
                os << e.mTag << " ";
                os << (e.mDim == 0 ? "0 0 0 " : "0 0 0 0 0 0 ");
                if (has_ptag)
                    os << "1 " << ptag << " ";
                else
                    os << "0 ";
                if (e.mDim > 0) {
                    if (bounds) {
                        os << bounds->size() << " ";
                        for (std::int32_t b : *bounds)
                            os << b << " ";
                        os << "\n";
                    } else {
                        os << "0\n";
                    }
                } else {
                    os << "\n";
                }
            }
        }
        if (binary)
            os << '\n';
        os << "$EndEntities\n";
    }

    os << "$Nodes\n";
    if (entities.empty()) {
        // Legacy path: one entity block covering every node.
        int node_dim = rMesh.NumCellBlocks() == 0 ? 0 : cell_dim(rMesh.Cells(0).Type());
        if (binary) {
            put_u64(1);
            put_u64(num_points);
            put_u64(1);
            put_u64(num_points);
            put_i32(node_dim);
            put_i32(0);
            put_i32(0);
            put_u64(num_points);
            // Node tags 1..num_points and the (3-padded) coords, each as one write
            // instead of a stream call per scalar (native endianness).
            std::vector<std::uint64_t> ntags(num_points);
            for (std::size_t i = 0; i < num_points; ++i)
                ntags[i] = i + 1;
            os.write(reinterpret_cast<const char*>(ntags.data()),
                     static_cast<std::streamsize>(num_points * 8));
            std::vector<double> cbuf(num_points * 3, 0.0);
            detail::dispatch_dtype(points.Dtype(), [&]<class T>() {
                const T* src = points.As<T>();
                parallel_for_bw(num_points, [&](std::size_t i) {
                    for (std::size_t c = 0; c < dim && c < 3; ++c)
                        cbuf[i * 3 + c] = static_cast<double>(src[i * dim + c]);
                });
            });
            os.write(reinterpret_cast<const char*>(cbuf.data()),
                     static_cast<std::streamsize>(num_points * 3 * 8));
            os << '\n';
        } else {
            os << "1 " << num_points << " 1 " << num_points << "\n";
            os << node_dim << " 0 0 " << num_points << "\n";
            for (std::size_t i = 0; i < num_points; ++i)
                os << (i + 1) << "\n";
            for (std::size_t i = 0; i < num_points; ++i)
                put_coords_ascii(i);
        }
    } else {
        // One block per entity that owns nodes -- an entity may own none (see
        // gmsh_entity_blocks_41), and an empty block would be legal but is not
        // what the reference writer emits. Node tags stay the global 1-based
        // point ids, so each block's list is sparse: exactly the case the
        // reader's tag->row remap exists for.
        std::size_t num_node_blocks = 0;
        for (const auto& e : entities)
            if (!e.mNodes.empty())
                ++num_node_blocks;
        if (binary) {
            put_u64(num_node_blocks);
            put_u64(num_points);
            put_u64(1);
            put_u64(num_points);
        } else {
            os << num_node_blocks << " " << num_points << " 1 " << num_points << "\n";
        }
        for (const auto& e : entities) {
            if (e.mNodes.empty())
                continue;
            const std::size_t nb = e.mNodes.size();
            if (binary) {
                put_i32(e.mDim);
                put_i32(e.mTag);
                put_i32(0);  // parametric
                put_u64(nb);
                std::vector<std::uint64_t> ntags(nb);
                for (std::size_t i = 0; i < nb; ++i)
                    ntags[i] = static_cast<std::uint64_t>(e.mNodes[i]) + 1;
                os.write(reinterpret_cast<const char*>(ntags.data()),
                         static_cast<std::streamsize>(nb * 8));
                std::vector<double> cbuf(nb * 3, 0.0);
                detail::dispatch_dtype(points.Dtype(), [&]<class T>() {
                    const T* src = points.As<T>();
                    parallel_for_bw(nb, [&](std::size_t i) {
                        const std::size_t p = static_cast<std::size_t>(e.mNodes[i]);
                        for (std::size_t c = 0; c < dim && c < 3; ++c)
                            cbuf[i * 3 + c] = static_cast<double>(src[p * dim + c]);
                    });
                });
                os.write(reinterpret_cast<const char*>(cbuf.data()),
                         static_cast<std::streamsize>(nb * 3 * 8));
            } else {
                os << e.mDim << " " << e.mTag << " 0 " << nb << "\n";
                for (std::size_t i = 0; i < nb; ++i)
                    os << (e.mNodes[i] + 1) << "\n";
                for (std::size_t i = 0; i < nb; ++i)
                    put_coords_ascii(static_cast<std::size_t>(e.mNodes[i]));
            }
        }
        if (binary)
            os << '\n';
    }
    os << "$EndNodes\n";

    // Elements: one block per cell block.
    std::size_t total_cells = 0;
    for (const auto cb : rMesh.CellRange())
        total_cells += cb.NumCells();
    const auto& m2g = meshio_to_gmsh_type();
    os << "$Elements\n";
    if (binary) {
        put_u64(rMesh.NumCellBlocks());
        put_u64(total_cells);
        put_u64(1);
        put_u64(total_cells);
    } else {
        os << rMesh.NumCellBlocks() << " " << total_cells << " 1 " << total_cells << "\n";
    }
    std::size_t tag0 = 1;
    for (std::size_t ci = 0; ci < rMesh.NumCellBlocks(); ++ci) {
        const auto cb = rMesh.Cells(ci);
        auto it = m2g.find(cb.Type());
        if (it == m2g.end())
            throw WriteError("Gmsh writer: unsupported cell type " + cb.Type());
        int gtype = it->second;
        int bdim = cell_dim(cb.Type());
        int entity_tag =
            (has_geometrical && ci < rMesh.CellDataNumBlocks("gmsh:geometrical") &&
             rMesh.CellData("gmsh:geometrical", ci).Size() > 0)
                ? static_cast<int>(detail::read_int(rMesh.CellData("gmsh:geometrical", ci), 0))
                : 0;
        const NDArray& conn = cb.Conn();
        std::size_t n = conn.Shape().size() >= 2 ? conn.Shape()[1] : 1;
        const std::vector<int>& perm = meshio_to_gmsh_perm(cb.Type());
        std::size_t count = cb.NumCells();
        if (binary) {
            put_i32(bdim);
            put_i32(entity_tag);
            put_i32(gtype);
            put_u64(count);
            // One buffer per block: [tag, node0..node(n-1)] u64, native, one write.
            const int* prm = perm.empty() ? nullptr : perm.data();
            const std::size_t stride = n + 1;
            std::vector<std::uint64_t> ebuf(count * stride);
            const std::uint64_t base = tag0;
            detail::dispatch_dtype(conn.Dtype(), [&]<class T>() {
                const T* src = conn.As<T>();
                parallel_for_bw(count, [&](std::size_t r) {
                    std::uint64_t* o = ebuf.data() + r * stride;
                    o[0] = base + r;
                    for (std::size_t j = 0; j < n; ++j) {
                        std::size_t sc = prm ? static_cast<std::size_t>(prm[j]) : j;
                        o[j + 1] = static_cast<std::uint64_t>(src[r * n + sc]) + 1;
                    }
                });
            });
            os.write(reinterpret_cast<const char*>(ebuf.data()),
                     static_cast<std::streamsize>(ebuf.size() * 8));
        } else {
            os << bdim << " " << entity_tag << " " << gtype << " " << count << "\n";
            for (std::size_t r = 0; r < count; ++r) {
                os << (tag0 + r);
                for (std::size_t j = 0; j < n; ++j) {
                    std::size_t src = perm.empty() ? j : static_cast<std::size_t>(perm[j]);
                    os << " " << (detail::read_int(conn, r * n + src) + 1);
                }
                os << "\n";
            }
        }
        tag0 += count;
    }
    if (binary)
        os << '\n';
    os << "$EndElements\n";

    for (const auto& name : rMesh.PointDataNames()) {
        if (name == "gmsh:dim_tags")
            continue;
        const NDArray& d = rMesh.PointData(name);
        std::size_t ncomp = d.Shape().size() >= 2 ? d.Shape()[1] : 1;
        std::size_t rows = d.Shape().empty() ? 0 : d.Shape()[0];
        os << "$NodeData\n1\n\"" << name << "\"\n1\n0\n3\n0\n" << ncomp << "\n" << rows << "\n";
        char buf[32];
        for (std::size_t r = 0; r < rows; ++r) {
            if (binary) {
                std::int32_t id = static_cast<std::int32_t>(r + 1);
                os.write(reinterpret_cast<const char*>(&id), 4);
                for (std::size_t c = 0; c < ncomp; ++c)
                    put_f64(detail::read_double(d, r * ncomp + c));
            } else {
                os << (r + 1);
                for (std::size_t c = 0; c < ncomp; ++c) {
                    std::snprintf(buf, sizeof(buf), " %.17g",
                                  detail::read_double(d, r * ncomp + c));
                    os << buf;
                }
                os << '\n';
            }
        }
        if (binary)
            os << '\n';
        os << "$EndNodeData\n";
    }

    for (const auto& name : rMesh.CellDataNames()) {
        if (name == "gmsh:physical" || name == "gmsh:geometrical" || name == "cell_tags")
            continue;
        write_data(os, "ElementData", name, rMesh, binary);
    }
}

MeshMetadata read_gmsh_metadata(const std::string& rPath, const ReadOptions& rOpts) {
    const detail::FileSource source(rPath, rOpts.mMmap);
    GmshCursor cur(source.View());

    if (gmsh_trim(cur.read_line()) != "$MeshFormat")
        throw ReadError("Expected $MeshFormat");
    std::istringstream fss(cur.read_line());
    std::string version;
    int file_type = 0, data_size = 8;
    fss >> version >> file_type >> data_size;
    const bool is_ascii = (file_type == 0);
    if (!is_ascii) {
        cur.read_i32();  // endianness marker
        if (cur.mPos < source.Size() && source.View()[cur.mPos] == '\n')
            ++cur.mPos;
    }
    cur.skip_to_end("MeshFormat");

    // Only 4.1 has typed element blocks whose headers make a summary cheap;
    // 2.2 stores a per-element type, so declining here costs a full read but
    // never a wrong answer (registry_read_metadata falls back).
    if (version != "4.1" && version != "4")
        throw ReadError("Gmsh: only format 4.1 has a header-only metadata path");

    GmshMeta41 meta;
    GmshEntities41 entities;
    // $PhysicalNames and $Entities are both small and both precede $Elements,
    // so a summary can report gmsh:physical and the named regions without ever
    // touching a coordinate or a connectivity row.
    std::unordered_map<std::string, NDArray> field_data;
    while (!cur.eof()) {
        const std::string line = cur.next_nonblank();
        if (line.empty())
            break;
        if (line[0] != '$')
            throw ReadError("Gmsh: unexpected line " + line);
        const std::string env = gmsh_trim(line.substr(1));
        if (env == "Nodes")
            gmsh_scan_nodes_41(cur, is_ascii, data_size, meta);
        else if (env == "Elements")
            gmsh_scan_elements_41(cur, is_ascii, data_size, meta);
        else if (env == "PhysicalNames")
            read_physical_names(cur, field_data);
        else if (env == "Entities")
            entities = read_entities_41(cur, is_ascii, data_size);
        else if (env == "NodeData")
            meta.mPointDataNames.push_back(gmsh_scan_data_name(cur, "NodeData"));
        else if (env == "ElementData")
            meta.mCellDataNames.push_back(gmsh_scan_data_name(cur, "ElementData"));
        else if (env == "Periodic")
            throw ReadError("Gmsh $" + env + " not supported by the C++ reader");
        else
            cur.skip_to_end(env);
    }

    MeshMetadata out;
    out.mNumPoints = meta.mNumPoints;
    out.mPointDim = 3;  // gmsh coordinates are always 3-D
    out.mCellBlocks = std::move(meta.mBlocks);
    out.mPointDataNames = std::move(meta.mPointDataNames);
    out.mCellDataNames = std::move(meta.mCellDataNames);
    // The reader synthesizes these from the entity structure, so a summary must
    // list them too or it would disagree with a real read.
    if (meta.mHasDimTags)
        out.mPointDataNames.push_back("gmsh:dim_tags");
    if (!out.mCellBlocks.empty()) {
        out.mCellDataNames.push_back("gmsh:geometrical");
        if (entities.mAnyPhysical)
            out.mCellDataNames.push_back("gmsh:physical");
    }
    std::sort(out.mPointDataNames.begin(), out.mPointDataNames.end());
    std::sort(out.mCellDataNames.begin(), out.mCellDataNames.end());

    // Named regions, counted from the block headers alone. gmsh_attach_regions
    // groups by (physical tag, block topological dimension), so this does too --
    // a summary that disagreed with a real read would be worse than none.
    if (entities.mAnyPhysical && !field_data.empty()) {
        std::map<std::pair<std::int64_t, int>, std::string> names;
        for (const auto& kv : field_data) {
            if (kv.second.Size() < 2)
                continue;
            names.emplace(std::make_pair(detail::read_int(kv.second, 0),
                                         static_cast<int>(detail::read_int(kv.second, 1))),
                          kv.first);
        }
        std::map<std::pair<std::int64_t, int>, std::size_t> counts;
        for (std::size_t b = 0; b < out.mCellBlocks.size() && b < meta.mBlockEntities.size(); ++b) {
            const auto& [edim, etag] = meta.mBlockEntities[b];
            if (edim < 0 || edim >= 4)
                continue;
            const auto pit = entities.mPhysical[static_cast<std::size_t>(edim)].find(etag);
            if (pit == entities.mPhysical[static_cast<std::size_t>(edim)].end() ||
                pit->second.empty())
                continue;
            int dim = cell_type_dimension(cell_type_from_name(out.mCellBlocks[b].mType));
            if (dim < 0) {
                const auto tit = topological_dimension().find(out.mCellBlocks[b].mType);
                dim = tit != topological_dimension().end() ? tit->second : -1;
            }
            counts[{pit->second[0], dim}] += out.mCellBlocks[b].mNumCells;
        }
        for (const auto& [key, name] : names) {
            const auto cit = counts.find(key);
            if (cit == counts.end())
                continue;
            RegionSummary rs;
            rs.mName = name;
            rs.mKind = RegionKind::Cell;
            rs.mDim = key.second;
            rs.mTag = key.first;
            rs.mNumEntries = cit->second;
            out.mRegions.push_back(std::move(rs));
        }
    }

    out.mHasBBox = false;  // would require decoding the coordinates
    return out;
}

}  // namespace meshioplusplus
