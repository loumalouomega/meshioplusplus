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
// glTF 2.0 writer. See formats/gltf.hpp for the contract. The Python reference
// (src/python/meshioplusplus/gltf/_gltf.py) follows the same steps in the same
// order with the same arithmetic, so the two write identical bytes.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/gltf.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/colormap.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/face_color.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/surface_distance.hpp"
#include "meshioplusplus/detail/surface_normals.hpp"
#include "meshioplusplus/detail/ragged_csr.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/skin.hpp"

namespace meshioplusplus {
namespace {

using detail::Vec3;

constexpr const char* kGltfPrefix = "meshio++: gltf: ";
constexpr const char* kGltfSourcePoint = "gltf:source_point";
constexpr double kGltfHalfSqrt2 = 0.70710678118654757;

constexpr int kGltfComponentUInt32 = 5125;
constexpr int kGltfComponentFloat32 = 5126;
constexpr int kGltfTargetArray = 34962;
constexpr int kGltfTargetElementArray = 34963;

constexpr int kGltfModePoints = 0;
constexpr int kGltfModeLines = 1;
constexpr int kGltfModeTriangles = 4;

bool gltf_starts_with(const std::string& rText, const char* pPrefix) {
    return rText.rfind(pPrefix, 0) == 0;
}

// --------------------------------------------------------------------------- //
// The JSON writer: a fixed key order, no whitespace, %.17g numbers.
// --------------------------------------------------------------------------- //
class GltfJson {
public:
    void BeginObject() {
        Separate();
        mOut += '{';
        mStack.push_back(false);
    }
    void EndObject() {
        mStack.pop_back();
        mOut += '}';
    }
    void BeginArray() {
        Separate();
        mOut += '[';
        mStack.push_back(false);
    }
    void EndArray() {
        mStack.pop_back();
        mOut += ']';
    }
    void Key(const char* pKey) {
        Separate();
        AppendString(pKey);
        mOut += ':';
        mAfterKey = true;
    }
    void Str(const std::string& rText) {
        Separate();
        AppendString(rText);
    }
    void Int(std::int64_t Value) {
        Separate();
        mOut += std::to_string(Value);
    }
    void Num(double Value) {
        Separate();
        char buf[40];
        detail::snprintf_c(buf, sizeof(buf), "%.17g", Value);
        mOut += buf;
    }
    void Bool(bool Value) {
        Separate();
        mOut += Value ? "true" : "false";
    }
    const std::string& Text() const { return mOut; }

private:
    void Separate() {
        if (mAfterKey) {
            mAfterKey = false;
            return;
        }
        if (!mStack.empty()) {
            if (mStack.back())
                mOut += ',';
            mStack.back() = true;
        }
    }

    // Escapes quote, backslash and control characters, and replaces every
    // ill-formed UTF-8 byte with U+FFFD, so the document is always valid JSON
    // whatever bytes a name carries.
    void AppendString(const std::string& rText) {
        mOut += '"';
        const std::size_t n = rText.size();
        std::size_t i = 0;
        while (i < n) {
            const unsigned char c = static_cast<unsigned char>(rText[i]);
            if (c < 0x80) {
                switch (c) {
                    case '"':
                        mOut += "\\\"";
                        break;
                    case '\\':
                        mOut += "\\\\";
                        break;
                    case '\b':
                        mOut += "\\b";
                        break;
                    case '\f':
                        mOut += "\\f";
                        break;
                    case '\n':
                        mOut += "\\n";
                        break;
                    case '\r':
                        mOut += "\\r";
                        break;
                    case '\t':
                        mOut += "\\t";
                        break;
                    default:
                        if (c < 0x20) {
                            char buf[8];
                            detail::snprintf_c(buf, sizeof(buf), "\\u%04x",
                                               static_cast<unsigned>(c));
                            mOut += buf;
                        } else {
                            mOut += static_cast<char>(c);
                        }
                }
                ++i;
                continue;
            }
            std::size_t len = 0;
            unsigned lo = 0x80;
            unsigned hi = 0xBF;
            if (c >= 0xC2 && c <= 0xDF) {
                len = 2;
            } else if (c >= 0xE0 && c <= 0xEF) {
                len = 3;
                lo = c == 0xE0 ? 0xA0 : 0x80;
                hi = c == 0xED ? 0x9F : 0xBF;
            } else if (c >= 0xF0 && c <= 0xF4) {
                len = 4;
                lo = c == 0xF0 ? 0x90 : 0x80;
                hi = c == 0xF4 ? 0x8F : 0xBF;
            }
            bool ok = len != 0 && i + len <= n;
            for (std::size_t k = 1; ok && k < len; ++k) {
                const unsigned cc = static_cast<unsigned char>(rText[i + k]);
                ok = k == 1 ? (cc >= lo && cc <= hi) : (cc >= 0x80 && cc <= 0xBF);
            }
            if (ok) {
                mOut.append(rText, i, len);
                i += len;
            } else {
                mOut += "\xEF\xBF\xBD";
                ++i;
            }
        }
        mOut += '"';
    }

    std::string mOut;
    std::vector<bool> mStack;
    bool mAfterKey = false;
};

// The binary chunk: little-endian whatever the host is.
class GltfBin {
public:
    void U32(std::uint32_t Value) {
        for (int k = 0; k < 4; ++k)
            mBytes += static_cast<char>((Value >> (8 * k)) & 0xFFu);
    }
    void F32(float Value) {
        std::uint32_t bits;
        std::memcpy(&bits, &Value, sizeof(bits));
        U32(bits);
    }
    std::size_t Size() const { return mBytes.size(); }
    const std::string& Bytes() const { return mBytes; }

private:
    std::string mBytes;
};

// --------------------------------------------------------------------------- //
// Collecting the surface.
// --------------------------------------------------------------------------- //

// What the writer exports, before any grouping or numbering: triangles with the
// facet they were fanned from (the unit a smooth fan may not cross) and the
// input cell that owns them (regions and cell-data colours), line segments, and
// vertex cells. Every id is an input point id.
struct GltfSurface {
    std::vector<std::array<std::int64_t, 3>> mTriVertices;
    std::vector<std::int64_t> mTriFacet;
    std::vector<std::int64_t> mTriCell;
    std::vector<std::array<std::int64_t, 2>> mLineVertices;
    std::vector<std::int64_t> mLineCell;
    std::vector<std::int64_t> mVertexPoint;
    std::vector<std::int64_t> mVertexCell;
};

// The ids of one cell, ragged or rectangular.
void gltf_cell_ids(const Mesh::CellView& rBlock, std::size_t Cell,
                   std::vector<std::int64_t>& rIds) {
    if (rBlock.IsRagged()) {
        const std::int64_t* pRow = rBlock.Row(Cell);
        rIds.assign(pRow, pRow + rBlock.RowSize(Cell));
        return;
    }
    const NDArray& conn = rBlock.Conn();
    const std::size_t npc = rBlock.NodesPerCell();
    rIds.resize(npc);
    for (std::size_t k = 0; k < npc; ++k)
        rIds[k] = detail::read_int(conn, Cell * npc + k);
}

using GltfKey4 = std::array<std::int64_t, 4>;

GltfKey4 gltf_sorted_key(const std::int64_t* pIds, std::size_t Count) {
    GltfKey4 key{-1, -1, -1, -1};
    for (std::size_t k = 0; k < Count && k < 4; ++k)
        key[k] = pIds[k];
    std::sort(key.begin(),
              key.begin() + static_cast<std::ptrdiff_t>(std::min<std::size_t>(Count, 4)));
    return key;
}

// A copy of the geometry of @p rMesh -- points and every cell block -- without
// its regions, which `surface_extract` would only warn about dropping.
Mesh gltf_geometry_copy(const Mesh& rMesh) {
    Mesh out;
    out.AssignPoints(detail::data_owned_copy(rMesh.Points()));
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsRagged()) {
            detail::append_ragged_copy(cb, out);
        } else {
            out.AddCellBlock(std::string(cb.Type()), detail::data_owned_copy(cb.Conn()));
        }
    }
    return out;
}

// The surface of @p rMesh: the skin of its volume blocks, its 2-D cells, its
// lines and its vertex cells. A 2-D cell that coincides with a skin facet wins,
// so a boundary patch keeps its own cells instead of z-fighting with the skin.
GltfSurface gltf_collect(const Mesh& rMesh, const std::vector<std::int64_t>& rBases) {
    GltfSurface out;

    struct Flat2D {
        std::int64_t mCell;
        std::size_t mBegin;
        std::size_t mCount;
    };
    std::vector<Flat2D> flat;
    std::vector<std::int64_t> flat_ids;
    std::vector<GltfKey4> keys2d;
    std::vector<std::string> skipped;
    std::vector<std::int64_t> ids;
    bool has_volume = false;

    std::size_t bi = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t base = rBases[bi++];
        const std::string type(cb.Type());
        const std::size_t ncells = cb.NumCells();
        const CellType ct = cell_type_from_name(type);
        if (cb.IsPolyhedron() || cell_type_dimension(ct) == 3) {
            has_volume = true;
            continue;
        }
        int kind = 0;  // 0 surface, 1 line, 2 vertex
        std::size_t corners = 0;
        if (gltf_starts_with(type, "triangle")) {
            corners = 3;
        } else if (gltf_starts_with(type, "quad")) {
            corners = 4;
        } else if (gltf_starts_with(type, "polygon")) {
            corners = 0;
        } else if (gltf_starts_with(type, "line")) {
            kind = 1;
        } else if (type == "vertex") {
            kind = 2;
        } else {
            skipped.push_back(type);
            continue;
        }
        for (std::size_t c = 0; c < ncells; ++c) {
            gltf_cell_ids(cb, c, ids);
            const std::int64_t cell = base + static_cast<std::int64_t>(c);
            if (kind == 0) {
                std::size_t k = corners == 0 ? ids.size() : std::min(corners, ids.size());
                if (k < 3)
                    continue;
                flat.push_back({cell, flat_ids.size(), k});
                flat_ids.insert(flat_ids.end(), ids.begin(),
                                ids.begin() + static_cast<std::ptrdiff_t>(k));
                if (k == 3 || k == 4)
                    keys2d.push_back(gltf_sorted_key(ids.data(), k));
            } else if (kind == 1) {
                if (ids.size() < 2 || ids[0] == ids[1])
                    continue;
                out.mLineVertices.push_back({ids[0], ids[1]});
                out.mLineCell.push_back(cell);
            } else {
                if (ids.empty())
                    continue;
                out.mVertexPoint.push_back(ids[0]);
                out.mVertexCell.push_back(cell);
            }
        }
    }
    if (!skipped.empty()) {
        std::string joined;
        for (std::size_t i = 0; i < skipped.size(); ++i)
            joined += (i ? ", " : "") + skipped[i];
        log::warn("gltf: cell block(s) of type {} have no glTF equivalent; skipping.", joined);
        detail::provenance_note("cells-dropped",
                                "cell block(s) of type " + joined + " have no glTF equivalent");
    }
    std::sort(keys2d.begin(), keys2d.end());

    // The skin of the volume blocks.
    std::int64_t facet = 0;
    if (has_volume) {
        if (has_skinnable_cells(rMesh)) {
            Mesh vol = gltf_geometry_copy(rMesh);
            NDArray source = NDArray::Uninit(DType::Int64, {rMesh.NumPoints()});
            for (std::size_t i = 0; i < rMesh.NumPoints(); ++i)
                source.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
            vol.AddPointData(kGltfSourcePoint, std::move(source));
            const Mesh skin = detail::surface_extract(vol, true, true, true, "gltf");
            const NDArray& src = skin.PointData(kGltfSourcePoint);

            std::vector<std::int64_t> parents;
            const std::string parent_key = "surface:parent_cell";
            for (std::size_t b = 0; b < skin.CellDataNumBlocks(parent_key); ++b) {
                const NDArray& a = skin.CellData(parent_key, b);
                for (std::size_t r = 0; r < skin.Cells(b).NumCells(); ++r)
                    parents.push_back(detail::read_int(a, r));
            }

            std::size_t cell_index = 0;
            for (const auto cb : skin.CellRange()) {
                const std::string type(cb.Type());
                const std::size_t ncells = cb.NumCells();
                if (type != "triangle" && type != "quad") {
                    cell_index += ncells;
                    facet += static_cast<std::int64_t>(ncells);
                    continue;
                }
                const std::size_t npc = type == "quad" ? 4 : 3;
                const NDArray& conn = cb.Conn();
                std::array<std::int64_t, 4> corner{};
                for (std::size_t c = 0; c < ncells; ++c) {
                    for (std::size_t k = 0; k < npc; ++k)
                        corner[k] = detail::read_int(
                            src, static_cast<std::size_t>(detail::read_int(conn, c * npc + k)));
                    const std::int64_t parent = parents[cell_index++];
                    const std::int64_t this_facet = facet++;
                    const GltfKey4 key = gltf_sorted_key(corner.data(), npc);
                    if (std::binary_search(keys2d.begin(), keys2d.end(), key))
                        continue;
                    for (std::size_t k = 1; k + 1 < npc; ++k) {
                        const std::array<std::int64_t, 3> tri{corner[0], corner[k], corner[k + 1]};
                        if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2])
                            continue;
                        out.mTriVertices.push_back(tri);
                        out.mTriFacet.push_back(this_facet);
                        out.mTriCell.push_back(parent);
                    }
                }
            }
        } else {
            log::warn(
                "gltf: the volume cell blocks are not supported by the skin extractor; skipping.");
        }
    }

    // The 2-D cells, fanned like convert_cells(simplexify).
    for (std::size_t f = 0; f < flat.size(); ++f) {
        const std::int64_t* p = flat_ids.data() + flat[f].mBegin;
        const std::size_t k = flat[f].mCount;
        const std::int64_t this_facet = facet + static_cast<std::int64_t>(f);
        for (std::size_t j = 1; j + 1 < k; ++j) {
            const std::array<std::int64_t, 3> tri{p[0], p[j], p[j + 1]};
            if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2])
                continue;
            out.mTriVertices.push_back(tri);
            out.mTriFacet.push_back(this_facet);
            out.mTriCell.push_back(flat[f].mCell);
        }
    }
    return out;
}

// --------------------------------------------------------------------------- //
// One glTF node per cell region.
// --------------------------------------------------------------------------- //
struct GltfNodes {
    std::vector<std::string> mNames;
    std::vector<std::int32_t> mCellNode;  // per global cell
};

// Each cell goes to the smallest cell region that contains it (ties: the first in
// the canonical (kind, name, dim, tag) order), regions of one name share a node,
// and a cell in no region goes to "unassigned".
GltfNodes gltf_assign_nodes(const Mesh& rMesh, std::size_t TotalCells, bool ByRegion) {
    GltfNodes out;
    out.mCellNode.assign(TotalCells, 0);
    std::vector<std::size_t> regions;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        if (r.mKind == RegionKind::Cell && r.NumEntries() > 0)
            regions.push_back(i);
    }
    if (!ByRegion || regions.empty()) {
        out.mNames = {"mesh"};
        return out;
    }

    std::vector<std::int32_t> region_node(regions.size());
    for (std::size_t r = 0; r < regions.size(); ++r) {
        const std::string& name = rMesh.Region(regions[r]).mName;
        auto it = std::find(out.mNames.begin(), out.mNames.end(), name);
        if (it == out.mNames.end()) {
            out.mNames.push_back(name);
            it = out.mNames.end() - 1;
        }
        region_node[r] = static_cast<std::int32_t>(it - out.mNames.begin());
    }

    std::vector<std::int32_t> best(TotalCells, -1);
    std::vector<std::size_t> best_size(TotalCells, 0);
    std::vector<std::uint8_t> members(TotalCells, 0);
    for (std::size_t r = 0; r < regions.size(); ++r) {
        const meshioplusplus::Region& reg = rMesh.Region(regions[r]);
        const std::size_t size = reg.NumEntries();
        for (std::size_t e = 0; e < size; ++e) {
            const std::int64_t g = detail::read_int(reg.mEntries, e);
            if (g < 0 || static_cast<std::size_t>(g) >= TotalCells)
                continue;
            const std::size_t gi = static_cast<std::size_t>(g);
            if (members[gi] < 255)
                ++members[gi];
            if (best[gi] < 0 || size < best_size[gi]) {
                best[gi] = static_cast<std::int32_t>(r);
                best_size[gi] = size;
            }
        }
    }
    std::int64_t overlapping = 0;
    bool any_unassigned = false;
    for (std::size_t g = 0; g < TotalCells; ++g) {
        if (members[g] > 1)
            ++overlapping;
        if (best[g] < 0)
            any_unassigned = true;
    }
    if (any_unassigned)
        out.mNames.push_back("unassigned");
    for (std::size_t g = 0; g < TotalCells; ++g)
        out.mCellNode[g] = best[g] < 0 ? static_cast<std::int32_t>(out.mNames.size() - 1)
                                       : region_node[static_cast<std::size_t>(best[g])];
    if (overlapping > 0) {
        log::debug("gltf: {} cell(s) belong to more than one region; each goes to the smallest.",
                   overlapping);
        detail::provenance_note("regions-overlap",
                                std::to_string(overlapping) +
                                    " cell(s) are in more than one region and were exported "
                                    "under the smallest");
    }
    return out;
}

// --------------------------------------------------------------------------- //
// Primitives.
// --------------------------------------------------------------------------- //
struct GltfVertexKey {
    std::int64_t mPoint;
    std::int64_t mGroup;  // the smooth fan, or -1
    std::int64_t mCell;   // the colouring cell, or -1
    bool operator==(const GltfVertexKey& rOther) const {
        return mPoint == rOther.mPoint && mGroup == rOther.mGroup && mCell == rOther.mCell;
    }
};

struct GltfVertexKeyHash {
    std::size_t operator()(const GltfVertexKey& rKey) const {
        std::size_t h = 0;
        for (std::int64_t v : {rKey.mPoint, rKey.mGroup, rKey.mCell})
            h ^= std::hash<std::int64_t>{}(v) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

// Vertices are numbered by first use, so the output does not depend on hashing.
class GltfPrimBuilder {
public:
    explicit GltfPrimBuilder(int Mode) : mMode(Mode) {}
    void Add(const GltfVertexKey& rKey) {
        auto [it, inserted] = mLookup.emplace(rKey, static_cast<std::uint32_t>(mVertices.size()));
        if (inserted) {
            if (mVertices.size() >= std::numeric_limits<std::uint32_t>::max() - 1)
                throw WriteError(std::string(kGltfPrefix) +
                                 "a primitive has more vertices than a 32-bit index can address");
            mVertices.push_back(rKey);
        }
        if (mMode != kGltfModePoints)
            mIndices.push_back(it->second);
    }
    bool Empty() const { return mVertices.empty(); }

    int mMode;
    std::vector<std::uint32_t> mIndices;
    std::vector<GltfVertexKey> mVertices;

private:
    std::unordered_map<GltfVertexKey, std::uint32_t, GltfVertexKeyHash> mLookup;
};

struct GltfNodeBuilders {
    GltfPrimBuilder mTriangles{kGltfModeTriangles};
    GltfPrimBuilder mLines{kGltfModeLines};
    GltfPrimBuilder mPoints{kGltfModePoints};
};

// --------------------------------------------------------------------------- //
// Colour and fields.
// --------------------------------------------------------------------------- //

// Reduce one row to a scalar: the requested component, or the magnitude, summed
// left to right with one sqrt at the end -- what the Python reference does.
double gltf_scalarize(const NDArray& rArray, std::size_t Row, std::size_t NumComponents,
                      const std::optional<int>& rComponent) {
    if (NumComponents == 0)
        return std::nan("");
    const std::size_t base = Row * NumComponents;
    if (rComponent.has_value()) {
        const int c = *rComponent;
        if (c < 0 || static_cast<std::size_t>(c) >= NumComponents)
            throw std::invalid_argument(std::string(kGltfPrefix) + "component " +
                                        std::to_string(c) + " is out of range for an array with " +
                                        std::to_string(NumComponents) + " component(s)");
        return detail::read_double(rArray, base + static_cast<std::size_t>(c));
    }
    if (NumComponents == 1)
        return detail::read_double(rArray, base);
    double sum = 0.0;
    for (std::size_t k = 0; k < NumComponents; ++k) {
        const double v = detail::read_double(rArray, base + k);
        sum += v * v;
    }
    return std::sqrt(sum);
}

detail::Rgb gltf_parse_color(const std::string& rText) {
    auto hex = [](char c) -> int {
        if (c >= '0' && c <= '9')
            return c - '0';
        if (c >= 'a' && c <= 'f')
            return c - 'a' + 10;
        if (c >= 'A' && c <= 'F')
            return c - 'A' + 10;
        return -1;
    };
    bool ok = rText.size() == 7 && rText[0] == '#';
    for (std::size_t i = 1; ok && i < 7; ++i)
        ok = hex(rText[i]) >= 0;
    if (!ok)
        throw std::invalid_argument(std::string(kGltfPrefix) + "nan_color '" + rText +
                                    "' is not a '#rrggbb' colour");
    auto byte = [&](std::size_t i) {
        return static_cast<std::uint8_t>(hex(rText[i]) * 16 + hex(rText[i + 1]));
    };
    return detail::Rgb{byte(1), byte(3), byte(5)};
}

float gltf_linear(std::uint8_t Byte) {
    const std::uint32_t bits = detail::srgb_to_linear_bits()[Byte];
    float out;
    std::memcpy(&out, &bits, sizeof(out));
    return out;
}

struct GltfColorPlan {
    bool mActive = false;
    bool mCellMode = false;
    std::vector<double> mValues;  // per point, or per global cell
    double mLo = 0.0;
    double mHi = 1.0;
    const std::uint8_t* mpTable = nullptr;
    detail::Rgb mNan{128, 128, 128};

    // The three linear floats of vertex @p rKey.
    std::array<float, 3> Color(const GltfVertexKey& rKey) const {
        const double v = mValues[static_cast<std::size_t>(mCellMode ? rKey.mCell : rKey.mPoint)];
        detail::Rgb c = mNan;
        if (std::isfinite(v))
            c = detail::colormap_lookup(mpTable, detail::color_param(v, mLo, mHi));
        return {gltf_linear(c.mR), gltf_linear(c.mG), gltf_linear(c.mB)};
    }
};

std::string gltf_available_arrays(const Mesh& rMesh) {
    std::vector<std::string> names = rMesh.PointDataNames();
    for (const std::string& n : rMesh.CellDataNames())
        names.push_back(n);
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());
    std::string out;
    for (const std::string& n : names)
        out += (out.empty() ? "" : ", ") + n;
    return out.empty() ? "none" : out;
}

struct GltfField {
    std::string mName;
    std::string mAttr;
    std::size_t mNumComp = 1;
    std::vector<float> mData;  // one row per input point
};

// `temperature` -> `_TEMPERATURE`: upper-cased ASCII, anything else `_`, and a
// numeric suffix when two names collapse onto one.
std::string gltf_attr_name(const std::string& rName, std::vector<std::string>& rTaken) {
    std::string base = "_";
    for (unsigned char c : rName) {
        if (c >= 'a' && c <= 'z')
            c = static_cast<unsigned char>(c - 'a' + 'A');
        const bool ok = (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_';
        base += ok ? static_cast<char>(c) : '_';
    }
    std::string attr = base;
    for (int suffix = 2; std::find(rTaken.begin(), rTaken.end(), attr) != rTaken.end(); ++suffix)
        attr = base + "_" + std::to_string(suffix);
    rTaken.push_back(attr);
    return attr;
}

const char* gltf_type_name(std::size_t NumComp) {
    switch (NumComp) {
        case 1:
            return "SCALAR";
        case 2:
            return "VEC2";
        case 3:
            return "VEC3";
        default:
            return "VEC4";
    }
}

struct GltfAccessor {
    int mComponentType = kGltfComponentFloat32;
    std::size_t mCount = 0;
    std::size_t mNumComp = 1;
    bool mMinMax = false;
    std::array<double, 3> mMin{0, 0, 0};
    std::array<double, 3> mMax{0, 0, 0};
    std::string mName;
    int mTarget = kGltfTargetArray;
    std::size_t mOffset = 0;
    std::size_t mLength = 0;
};

struct GltfPrimOut {
    int mMode = kGltfModeTriangles;
    int mIndices = -1;
    std::vector<std::pair<std::string, int>> mAttrs;
};

struct GltfNodeOut {
    std::string mName;
    std::vector<GltfPrimOut> mPrims;
};

struct GltfRendered {
    std::string mJson;
    std::string mBin;
};

GltfRendered gltf_render(const Mesh& rMesh, const GltfWriteOptions& rOpt, bool Binary,
                         const std::string& rBinUri) {
    const std::size_t n = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    // Validate everything cheap before any real work.
    if (rOpt.mNormals && !(rOpt.mSplitAngle >= 0.0 && rOpt.mSplitAngle <= 180.0))
        throw std::invalid_argument(std::string(kGltfPrefix) +
                                    "split_angle must lie in [0, 180] degrees");
    if (!(rOpt.mScale > 0.0) || !std::isfinite(rOpt.mScale))
        throw std::invalid_argument(std::string(kGltfPrefix) + "scale must be a positive number");
    const bool color_active = !rOpt.mColorBy.empty();
    const std::uint8_t* table = nullptr;
    detail::Rgb nan_color{128, 128, 128};
    if (color_active) {
        table = detail::colormap_table(rOpt.mCmap);
        nan_color = gltf_parse_color(rOpt.mNanColor);
        if (rOpt.mVMin.has_value() && rOpt.mVMax.has_value() && *rOpt.mVMin > *rOpt.mVMax)
            throw std::invalid_argument(std::string(kGltfPrefix) + "vmin must not exceed vmax");
    }
    const bool color_point = color_active && rMesh.HasPointData(rOpt.mColorBy);
    const bool color_cell = color_active && !color_point && rMesh.HasCellData(rOpt.mColorBy);
    if (color_active && !color_point && !color_cell)
        throw std::invalid_argument(std::string(kGltfPrefix) + "color_by array '" + rOpt.mColorBy +
                                    "' is in neither point_data nor cell_data (available: " +
                                    gltf_available_arrays(rMesh) + ")");

    // The points, as doubles, 2-D padded with z = 0.
    detail::TriangleSoup soup;
    soup.mPoints.resize(n);
    for (std::size_t p = 0; p < n; ++p)
        soup.mPoints[p] = detail::read_point(rMesh.Points(), dim, static_cast<std::int64_t>(p));
    const std::vector<Vec3>& pts = soup.mPoints;

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::size_t total_cells = static_cast<std::size_t>(detail::total_cells(bases));
    GltfSurface surf = gltf_collect(rMesh, bases);
    const GltfNodes nodes = gltf_assign_nodes(rMesh, total_cells, rOpt.mByRegion);

    // Smooth fans over the whole soup, so smoothing crosses region borders.
    const std::size_t ntri = surf.mTriVertices.size();
    detail::VertexNormalGroups groups;
    if (rOpt.mNormals && ntri > 0) {
        soup.mVertices = surf.mTriVertices;
        soup.mSourceCell = surf.mTriFacet;
        soup.mCorners.reserve(ntri * 3);
        for (const auto& tri : surf.mTriVertices)
            for (std::int64_t v : tri)
                soup.mCorners.push_back(pts[static_cast<std::size_t>(v)]);
        groups = detail::vertex_normal_groups(soup, rOpt.mNormalWeight, rOpt.mSplitAngle);
    }

    // Primitives, one set per node, vertices numbered by first use.
    std::vector<GltfNodeBuilders> builders(nodes.mNames.size());
    auto node_of = [&](std::int64_t Cell) {
        return static_cast<std::size_t>(nodes.mCellNode[static_cast<std::size_t>(Cell)]);
    };
    for (std::size_t t = 0; t < ntri; ++t) {
        GltfPrimBuilder& b = builders[node_of(surf.mTriCell[t])].mTriangles;
        for (std::size_t i = 0; i < 3; ++i)
            b.Add({surf.mTriVertices[t][i],
                   rOpt.mNormals ? groups.mCornerGroup[t * 3 + i] : std::int64_t{-1},
                   color_cell ? surf.mTriCell[t] : std::int64_t{-1}});
    }
    for (std::size_t s = 0; s < surf.mLineVertices.size(); ++s) {
        GltfPrimBuilder& b = builders[node_of(surf.mLineCell[s])].mLines;
        for (std::int64_t v : surf.mLineVertices[s])
            b.Add({v, -1, color_cell ? surf.mLineCell[s] : std::int64_t{-1}});
    }
    for (std::size_t s = 0; s < surf.mVertexPoint.size(); ++s)
        builders[node_of(surf.mVertexCell[s])].mPoints.Add(
            {surf.mVertexPoint[s], -1, color_cell ? surf.mVertexCell[s] : std::int64_t{-1}});
    if (rMesh.NumCellBlocks() == 0)
        for (std::size_t p = 0; p < n; ++p)
            builders[0].mPoints.Add({static_cast<std::int64_t>(p), -1, -1});

    // The points that are exported: the bounding box, the up axis and the field
    // checks look at these and no others.
    std::vector<char> used(n, 0);
    for (const GltfNodeBuilders& nb : builders)
        for (const GltfPrimBuilder* pb : {&nb.mTriangles, &nb.mLines, &nb.mPoints})
            for (const GltfVertexKey& k : pb->mVertices)
                used[static_cast<std::size_t>(k.mPoint)] = 1;
    Vec3 lo{0, 0, 0};
    Vec3 hi{0, 0, 0};
    bool any_used = false;
    bool flat = dim < 3;
    for (std::size_t p = 0; p < n; ++p) {
        if (!used[p])
            continue;
        for (std::size_t k = 0; k < 3; ++k)
            if (!std::isfinite(pts[p][k]))
                throw WriteError(std::string(kGltfPrefix) + "point " + std::to_string(p) +
                                 " has a non-finite coordinate");
        if (!any_used) {
            lo = pts[p];
            hi = pts[p];
            any_used = true;
        } else {
            for (std::size_t k = 0; k < 3; ++k) {
                lo[k] = pts[p][k] < lo[k] ? pts[p][k] : lo[k];
                hi[k] = pts[p][k] > hi[k] ? pts[p][k] : hi[k];
            }
        }
        if (std::abs(pts[p][2]) > 1e-14)
            flat = false;
    }
    Vec3 center{0, 0, 0};
    if (rOpt.mRecenter && any_used)
        for (std::size_t k = 0; k < 3; ++k)
            center[k] = 0.5 * (lo[k] + hi[k]);

    // The colour values and their range over the exported vertices.
    GltfColorPlan color;
    if (color_active) {
        color.mActive = true;
        color.mCellMode = color_cell;
        color.mpTable = table;
        color.mNan = nan_color;
        if (color_point) {
            const NDArray& arr = rMesh.PointData(rOpt.mColorBy);
            const std::size_t ncomp = n == 0 ? 0 : arr.Size() / n;
            color.mValues.assign(n, std::nan(""));
            for (std::size_t p = 0; p < n; ++p)
                if (used[p])
                    color.mValues[p] = gltf_scalarize(arr, p, ncomp, rOpt.mComponent);
        } else {
            const std::size_t blocks = rMesh.NumCellBlocks();
            if (rMesh.CellDataNumBlocks(rOpt.mColorBy) != blocks)
                throw std::invalid_argument(std::string(kGltfPrefix) + "cell_data array '" +
                                            rOpt.mColorBy + "' has " +
                                            std::to_string(rMesh.CellDataNumBlocks(rOpt.mColorBy)) +
                                            " block(s) but the mesh has " + std::to_string(blocks));
            for (std::size_t b = 0; b < blocks; ++b) {
                const std::size_t nc = rMesh.Cells(b).NumCells();
                const NDArray& arr = rMesh.CellData(rOpt.mColorBy, b);
                const std::size_t ncomp = nc == 0 ? 0 : arr.Size() / nc;
                for (std::size_t r = 0; r < nc; ++r)
                    color.mValues.push_back(gltf_scalarize(arr, r, ncomp, rOpt.mComponent));
            }
        }
        bool seen = false;
        double vlo = 0.0;
        double vhi = 0.0;
        for (const GltfNodeBuilders& nb : builders)
            for (const GltfPrimBuilder* pb : {&nb.mTriangles, &nb.mLines, &nb.mPoints})
                for (const GltfVertexKey& k : pb->mVertices) {
                    const double v =
                        color.mValues[static_cast<std::size_t>(color_cell ? k.mCell : k.mPoint)];
                    if (!std::isfinite(v))
                        continue;
                    if (!seen) {
                        vlo = v;
                        vhi = v;
                        seen = true;
                    } else {
                        vlo = v < vlo ? v : vlo;
                        vhi = v > vhi ? v : vhi;
                    }
                }
        color.mLo = rOpt.mVMin.has_value() ? *rOpt.mVMin : vlo;
        color.mHi = rOpt.mVMax.has_value() ? *rOpt.mVMax : vhi;
        if (color.mLo > color.mHi)
            throw std::invalid_argument(std::string(kGltfPrefix) + "vmin must not exceed vmax");
    }

    // Raw fields: every one-to-four component point_data array, as float32.
    std::vector<GltfField> fields;
    if (rOpt.mFields) {
        std::vector<std::string> taken;
        for (const std::string& name : rMesh.PointDataNames()) {
            if (name == "normals")
                continue;
            const NDArray& arr = rMesh.PointData(name);
            const std::vector<std::size_t>& shape = arr.Shape();
            std::size_t ncomp = 0;
            if (n > 0 && !shape.empty() && shape[0] == n && shape.size() <= 2)
                ncomp = shape.size() == 1 ? 1 : shape[1];
            if (ncomp < 1 || ncomp > 4) {
                log::warn(
                    "gltf: point_data '{}' has no glTF attribute type (1 to 4 components "
                    "per point); skipping.",
                    name);
                detail::provenance_note(
                    "fields-dropped", "point_data '" + name + "' is not a 1 to 4 component array");
                continue;
            }
            GltfField f;
            f.mName = name;
            f.mNumComp = ncomp;
            f.mData.resize(n * ncomp);
            bool finite = true;
            for (std::size_t p = 0; p < n; ++p)
                for (std::size_t c = 0; c < ncomp; ++c) {
                    const float v = static_cast<float>(detail::read_double(arr, p * ncomp + c));
                    f.mData[p * ncomp + c] = v;
                    if (used[p] && !std::isfinite(v))
                        finite = false;
                }
            if (!finite) {
                log::warn(
                    "gltf: point_data '{}' has non-finite values, which glTF cannot hold; "
                    "skipping.",
                    name);
                detail::provenance_note("fields-dropped",
                                        "point_data '" + name + "' has non-finite values");
                continue;
            }
            f.mAttr = gltf_attr_name(name, taken);
            fields.push_back(std::move(f));
        }
    }

    // Point normals, for a POINTS primitive.
    std::vector<Vec3> point_normals;
    std::vector<char> point_normal_ok;
    if (rOpt.mNormals && rMesh.HasPointData("normals")) {
        const NDArray& arr = rMesh.PointData("normals");
        const std::vector<std::size_t>& shape = arr.Shape();
        if (n > 0 && shape.size() == 2 && shape[0] == n && shape[1] == 3) {
            point_normals.resize(n);
            point_normal_ok.assign(n, 0);
            for (std::size_t p = 0; p < n; ++p) {
                const Vec3 v{detail::read_double(arr, p * 3), detail::read_double(arr, p * 3 + 1),
                             detail::read_double(arr, p * 3 + 2)};
                const double len = detail::vec3_norm(v);
                if (std::isfinite(len) && len > 0.0) {
                    point_normals[p] = detail::vec3_scale(v, 1.0 / len);
                    point_normal_ok[p] = 1;
                }
            }
        }
    }

    // ---- the binary chunk and the accessors ----
    GltfBin bin;
    std::vector<GltfAccessor> accessors;
    auto add_floats = [&](const std::vector<float>& rData, std::size_t NumComp, bool MinMax,
                          const std::string& rName) {
        GltfAccessor a;
        a.mNumComp = NumComp;
        a.mCount = rData.size() / NumComp;
        a.mName = rName;
        a.mOffset = bin.Size();
        a.mLength = rData.size() * 4;
        if (MinMax && a.mCount > 0) {
            for (std::size_t k = 0; k < 3; ++k) {
                a.mMin[k] = rData[k];
                a.mMax[k] = rData[k];
            }
            for (std::size_t v = 1; v < a.mCount; ++v)
                for (std::size_t k = 0; k < 3; ++k) {
                    const double x = rData[v * 3 + k];
                    a.mMin[k] = x < a.mMin[k] ? x : a.mMin[k];
                    a.mMax[k] = x > a.mMax[k] ? x : a.mMax[k];
                }
            a.mMinMax = true;
        }
        for (float f : rData)
            bin.F32(f);
        accessors.push_back(std::move(a));
        return static_cast<int>(accessors.size() - 1);
    };
    auto add_indices = [&](const std::vector<std::uint32_t>& rIdx) {
        GltfAccessor a;
        a.mComponentType = kGltfComponentUInt32;
        a.mTarget = kGltfTargetElementArray;
        a.mCount = rIdx.size();
        a.mOffset = bin.Size();
        a.mLength = rIdx.size() * 4;
        for (std::uint32_t v : rIdx)
            bin.U32(v);
        accessors.push_back(std::move(a));
        return static_cast<int>(accessors.size() - 1);
    };

    std::vector<GltfNodeOut> out_nodes;
    for (std::size_t ni = 0; ni < builders.size(); ++ni) {
        GltfNodeOut node;
        node.mName = nodes.mNames[ni];
        const GltfPrimBuilder* prims[3] = {&builders[ni].mTriangles, &builders[ni].mLines,
                                           &builders[ni].mPoints};
        for (const GltfPrimBuilder* pb : prims) {
            if (pb->Empty())
                continue;
            const std::size_t nv = pb->mVertices.size();
            GltfPrimOut po;
            po.mMode = pb->mMode;
            if (pb->mMode != kGltfModePoints)
                po.mIndices = add_indices(pb->mIndices);

            std::vector<float> position(nv * 3);
            for (std::size_t v = 0; v < nv; ++v) {
                const Vec3& p = pts[static_cast<std::size_t>(pb->mVertices[v].mPoint)];
                for (std::size_t k = 0; k < 3; ++k)
                    position[v * 3 + k] = static_cast<float>(p[k] - center[k]);
            }
            po.mAttrs.emplace_back("POSITION", add_floats(position, 3, true, ""));

            if (pb->mMode == kGltfModeTriangles && rOpt.mNormals) {
                std::vector<float> normal(nv * 3);
                for (std::size_t v = 0; v < nv; ++v) {
                    const Vec3& g =
                        groups.mGroupNormal[static_cast<std::size_t>(pb->mVertices[v].mGroup)];
                    const bool defined = g[0] != 0.0 || g[1] != 0.0 || g[2] != 0.0;
                    const Vec3 u = defined ? g : Vec3{0.0, 0.0, 1.0};
                    for (std::size_t k = 0; k < 3; ++k)
                        normal[v * 3 + k] = static_cast<float>(u[k]);
                }
                po.mAttrs.emplace_back("NORMAL", add_floats(normal, 3, false, ""));
            } else if (pb->mMode == kGltfModePoints && !point_normals.empty()) {
                bool all_ok = true;
                for (const GltfVertexKey& k : pb->mVertices)
                    all_ok = all_ok && point_normal_ok[static_cast<std::size_t>(k.mPoint)];
                if (all_ok) {
                    std::vector<float> normal(nv * 3);
                    for (std::size_t v = 0; v < nv; ++v)
                        for (std::size_t k = 0; k < 3; ++k)
                            normal[v * 3 + k] = static_cast<float>(
                                point_normals[static_cast<std::size_t>(pb->mVertices[v].mPoint)]
                                             [k]);
                    po.mAttrs.emplace_back("NORMAL", add_floats(normal, 3, false, ""));
                } else {
                    log::warn(
                        "gltf: point_data 'normals' has a zero or non-finite row on the "
                        "exported points; NORMAL omitted.");
                }
            }

            if (color.mActive) {
                std::vector<float> rgb(nv * 3);
                for (std::size_t v = 0; v < nv; ++v) {
                    const std::array<float, 3> c = color.Color(pb->mVertices[v]);
                    for (std::size_t k = 0; k < 3; ++k)
                        rgb[v * 3 + k] = c[k];
                }
                po.mAttrs.emplace_back("COLOR_0", add_floats(rgb, 3, false, ""));
            }
            for (const GltfField& f : fields) {
                std::vector<float> data(nv * f.mNumComp);
                for (std::size_t v = 0; v < nv; ++v) {
                    const std::size_t p = static_cast<std::size_t>(pb->mVertices[v].mPoint);
                    for (std::size_t c = 0; c < f.mNumComp; ++c)
                        data[v * f.mNumComp + c] = f.mData[p * f.mNumComp + c];
                }
                po.mAttrs.emplace_back(f.mAttr, add_floats(data, f.mNumComp, false, f.mName));
            }
            node.mPrims.push_back(std::move(po));
        }
        if (!node.mPrims.empty())
            out_nodes.push_back(std::move(node));
    }
    if (out_nodes.empty())
        log::warn(
            "gltf: nothing to export (no surface, line or vertex cells); writing an empty "
            "scene.");

    // ---- the root transform ----
    GltfUpAxis up = rOpt.mUpAxis;
    if (up == GltfUpAxis::Auto)
        up = flat ? GltfUpAxis::Y : GltfUpAxis::Z;
    const double s = rOpt.mScale;
    const double sx = s * center[0];
    const double sy = s * center[1];
    const double sz = s * center[2];
    // world = R * (s * v) + R * (s * center), so the source coordinates come back.
    std::array<double, 3> translation{sx, sy, sz};
    if (up == GltfUpAxis::Z)
        translation = {sx, sz, -sy};
    else if (up == GltfUpAxis::X)
        translation = {-sy, sx, sz};
    for (double& t : translation)
        t += 0.0;  // -0 -> 0

    // ---- the JSON ----
    GltfJson js;
    js.BeginObject();
    js.Key("asset");
    js.BeginObject();
    const std::vector<std::string> provenance = detail::provenance_lines(detail::SlotTier::Block);
    js.Key("generator");
    js.Str(provenance.empty() ? std::string(detail::kProvenanceTag) : provenance[0]);
    js.Key("version");
    js.Str("2.0");
    if (provenance.size() > 1) {
        js.Key("extras");
        js.BeginObject();
        js.Key("meshioplusplus:provenance");
        js.BeginArray();
        for (const std::string& line : provenance)
            js.Str(line);
        js.EndArray();
        js.EndObject();
    }
    js.EndObject();
    const bool unlit = color.mActive && rOpt.mUnlit;
    if (unlit) {
        js.Key("extensionsUsed");
        js.BeginArray();
        js.Str("KHR_materials_unlit");
        js.EndArray();
    }
    js.Key("scene");
    js.Int(0);
    js.Key("scenes");
    js.BeginArray();
    js.BeginObject();
    if (!out_nodes.empty()) {
        js.Key("nodes");
        js.BeginArray();
        js.Int(0);
        js.EndArray();
    }
    js.EndObject();
    js.EndArray();

    if (!out_nodes.empty()) {
        js.Key("nodes");
        js.BeginArray();
        js.BeginObject();
        js.Key("name");
        js.Str("meshio++");
        js.Key("children");
        js.BeginArray();
        for (std::size_t i = 0; i < out_nodes.size(); ++i)
            js.Int(static_cast<std::int64_t>(i + 1));
        js.EndArray();
        if (up == GltfUpAxis::Z || up == GltfUpAxis::X) {
            js.Key("rotation");
            js.BeginArray();
            if (up == GltfUpAxis::Z) {
                js.Num(-kGltfHalfSqrt2);
                js.Num(0.0);
                js.Num(0.0);
                js.Num(kGltfHalfSqrt2);
            } else {
                js.Num(0.0);
                js.Num(0.0);
                js.Num(kGltfHalfSqrt2);
                js.Num(kGltfHalfSqrt2);
            }
            js.EndArray();
        }
        if (s != 1.0) {
            js.Key("scale");
            js.BeginArray();
            for (int k = 0; k < 3; ++k)
                js.Num(s);
            js.EndArray();
        }
        js.Key("translation");
        js.BeginArray();
        for (double t : translation)
            js.Num(t);
        js.EndArray();
        js.EndObject();
        for (std::size_t i = 0; i < out_nodes.size(); ++i) {
            js.BeginObject();
            js.Key("name");
            js.Str(out_nodes[i].mName);
            js.Key("mesh");
            js.Int(static_cast<std::int64_t>(i));
            js.EndObject();
        }
        js.EndArray();

        js.Key("meshes");
        js.BeginArray();
        for (const GltfNodeOut& node : out_nodes) {
            js.BeginObject();
            js.Key("name");
            js.Str(node.mName);
            js.Key("primitives");
            js.BeginArray();
            for (const GltfPrimOut& po : node.mPrims) {
                js.BeginObject();
                js.Key("attributes");
                js.BeginObject();
                for (const auto& [semantic, index] : po.mAttrs) {
                    js.Key(semantic.c_str());
                    js.Int(index);
                }
                js.EndObject();
                if (po.mIndices >= 0) {
                    js.Key("indices");
                    js.Int(po.mIndices);
                }
                js.Key("material");
                js.Int(0);
                js.Key("mode");
                js.Int(po.mMode);
                js.EndObject();
            }
            js.EndArray();
            js.EndObject();
        }
        js.EndArray();

        js.Key("materials");
        js.BeginArray();
        js.BeginObject();
        js.Key("name");
        js.Str("meshio++");
        js.Key("pbrMetallicRoughness");
        js.BeginObject();
        js.Key("baseColorFactor");
        js.BeginArray();
        const double base = color.mActive ? 1.0 : 0.8;
        js.Num(base);
        js.Num(base);
        js.Num(base);
        js.Num(1.0);
        js.EndArray();
        js.Key("metallicFactor");
        js.Num(0.0);
        js.Key("roughnessFactor");
        js.Num(0.5);
        js.EndObject();
        if (unlit) {
            js.Key("extensions");
            js.BeginObject();
            js.Key("KHR_materials_unlit");
            js.BeginObject();
            js.EndObject();
            js.EndObject();
        }
        js.Key("doubleSided");
        js.Bool(true);
        js.EndObject();
        js.EndArray();

        js.Key("accessors");
        js.BeginArray();
        for (std::size_t i = 0; i < accessors.size(); ++i) {
            const GltfAccessor& a = accessors[i];
            js.BeginObject();
            js.Key("bufferView");
            js.Int(static_cast<std::int64_t>(i));
            js.Key("componentType");
            js.Int(a.mComponentType);
            js.Key("count");
            js.Int(static_cast<std::int64_t>(a.mCount));
            js.Key("type");
            js.Str(gltf_type_name(a.mNumComp));
            if (a.mMinMax) {
                js.Key("max");
                js.BeginArray();
                for (double v : a.mMax)
                    js.Num(v);
                js.EndArray();
                js.Key("min");
                js.BeginArray();
                for (double v : a.mMin)
                    js.Num(v);
                js.EndArray();
            }
            if (!a.mName.empty()) {
                js.Key("name");
                js.Str(a.mName);
            }
            js.EndObject();
        }
        js.EndArray();

        js.Key("bufferViews");
        js.BeginArray();
        for (const GltfAccessor& a : accessors) {
            js.BeginObject();
            js.Key("buffer");
            js.Int(0);
            js.Key("byteOffset");
            js.Int(static_cast<std::int64_t>(a.mOffset));
            js.Key("byteLength");
            js.Int(static_cast<std::int64_t>(a.mLength));
            js.Key("target");
            js.Int(a.mTarget);
            js.EndObject();
        }
        js.EndArray();

        js.Key("buffers");
        js.BeginArray();
        js.BeginObject();
        js.Key("byteLength");
        js.Int(static_cast<std::int64_t>(bin.Size()));
        if (!Binary) {
            js.Key("uri");
            js.Str(rBinUri);
        }
        js.EndObject();
        js.EndArray();
    }
    js.EndObject();

    GltfRendered out;
    out.mJson = js.Text();
    out.mBin = bin.Bytes();
    return out;
}

std::string gltf_lower(std::string Text) {
    for (char& c : Text)
        if (c >= 'A' && c <= 'Z')
            c = static_cast<char>(c - 'A' + 'a');
    return Text;
}

bool gltf_has_suffix(const std::string& rText, const char* pSuffix) {
    const std::string s = pSuffix;
    return rText.size() >= s.size() && gltf_lower(rText.substr(rText.size() - s.size())) == s;
}

// The basename of the sidecar, percent-encoded (RFC 3986 unreserved characters
// are kept), as a glTF `uri`.
std::string gltf_bin_uri(const std::string& rBinPath) {
    const std::size_t cut = rBinPath.find_last_of("/\\");
    const std::string name = cut == std::string::npos ? rBinPath : rBinPath.substr(cut + 1);
    std::string out;
    for (unsigned char c : name) {
        const bool keep = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
                          (c >= '0' && c <= '9') || c == '-' || c == '.' || c == '_' || c == '~';
        if (keep) {
            out += static_cast<char>(c);
        } else {
            char buf[8];
            detail::snprintf_c(buf, sizeof(buf), "%%%02X", static_cast<unsigned>(c));
            out += buf;
        }
    }
    return out;
}

void gltf_put_u32(std::string& rOut, std::uint32_t Value) {
    for (int k = 0; k < 4; ++k)
        rOut += static_cast<char>((Value >> (8 * k)) & 0xFFu);
}

}  // namespace

GltfContainer gltf_container_from_name(const std::string& rName) {
    if (rName == "auto")
        return GltfContainer::Auto;
    if (rName == "glb" || rName == "binary")
        return GltfContainer::Binary;
    if (rName == "gltf" || rName == "json")
        return GltfContainer::Json;
    throw std::invalid_argument(std::string(kGltfPrefix) + "unknown container '" + rName +
                                "' (expected 'auto', 'glb' or 'gltf')");
}

GltfUpAxis gltf_up_axis_from_name(const std::string& rName) {
    if (rName == "auto")
        return GltfUpAxis::Auto;
    if (rName == "x")
        return GltfUpAxis::X;
    if (rName == "y")
        return GltfUpAxis::Y;
    if (rName == "z")
        return GltfUpAxis::Z;
    throw std::invalid_argument(std::string(kGltfPrefix) + "unknown up axis '" + rName +
                                "' (expected 'auto', 'x', 'y' or 'z')");
}

void write_gltf(const std::string& rPath, const Mesh& rMesh, const GltfWriteOptions& rOptions) {
    const bool binary =
        rOptions.mContainer == GltfContainer::Binary ||
        (rOptions.mContainer == GltfContainer::Auto && gltf_has_suffix(rPath, ".glb"));
    std::string bin_path;
    if (!binary)
        bin_path =
            (gltf_has_suffix(rPath, ".gltf") ? rPath.substr(0, rPath.size() - 5) : rPath) + ".bin";
    const GltfRendered r =
        gltf_render(rMesh, rOptions, binary, binary ? "" : gltf_bin_uri(bin_path));

    if (binary) {
        std::string json = r.mJson;
        while (json.size() % 4 != 0)
            json += ' ';
        std::string bin = r.mBin;
        while (bin.size() % 4 != 0)
            bin += '\0';
        const std::uint64_t total =
            12ull + 8ull + json.size() + (bin.empty() ? 0ull : 8ull + bin.size());
        if (total > std::numeric_limits<std::uint32_t>::max())
            throw WriteError(std::string(kGltfPrefix) +
                             "the file is larger than a GLB can hold (4 GiB)");
        std::string out;
        out.reserve(static_cast<std::size_t>(total));
        gltf_put_u32(out, 0x46546C67u);
        gltf_put_u32(out, 2u);
        gltf_put_u32(out, static_cast<std::uint32_t>(total));
        gltf_put_u32(out, static_cast<std::uint32_t>(json.size()));
        gltf_put_u32(out, 0x4E4F534Au);
        out += json;
        if (!bin.empty()) {
            gltf_put_u32(out, static_cast<std::uint32_t>(bin.size()));
            gltf_put_u32(out, 0x004E4942u);
            out += bin;
        }
        auto f = detail::make_classic_ofstream(rPath, std::ios::out | std::ios::binary);
        if (!f)
            throw WriteError(std::string(kGltfPrefix) + "cannot open '" + rPath + "' for writing");
        f.write(out.data(), static_cast<std::streamsize>(out.size()));
        return;
    }

    auto f = detail::make_classic_ofstream(rPath, std::ios::out | std::ios::binary);
    if (!f)
        throw WriteError(std::string(kGltfPrefix) + "cannot open '" + rPath + "' for writing");
    f.write(r.mJson.data(), static_cast<std::streamsize>(r.mJson.size()));
    if (!r.mBin.empty()) {
        auto b = detail::make_classic_ofstream(bin_path, std::ios::out | std::ios::binary);
        if (!b)
            throw WriteError(std::string(kGltfPrefix) + "cannot open '" + bin_path +
                             "' for writing");
        b.write(r.mBin.data(), static_cast<std::streamsize>(r.mBin.size()));
    }
}

}  // namespace meshioplusplus
