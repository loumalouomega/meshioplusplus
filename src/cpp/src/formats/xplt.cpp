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
// FEBio plot file (`.xplt`) reader. The chunk layout follows FEBio's writer
// (FEBioPlot/FEBioPlotFile.cpp, PltArchive.cpp) and FEBio Studio's reader
// (XPLTLib/xpltReader3.cpp); see doc/formats/xplt.md.
// Python twin: src/python/meshioplusplus/xplt/_xplt.py.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/xplt.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/facet_index.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/detail/node_order.hpp"
#include "meshioplusplus/detail/zlib_inflate.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

// Chunk ids (FEBioPlotFile.h).
constexpr std::uint32_t kXpltMagic = 0x00464542;
constexpr std::uint32_t kRoot = 0x01000000, kHeader = 0x01010000, kHdrVersion = 0x01010001,
                        kHdrCompression = 0x01010004, kDictionary = 0x01020000,
                        kDicItem = 0x01020001, kDicItemType = 0x01020002, kDicItemFmt = 0x01020003,
                        kDicItemName = 0x01020004, kDicItemArraySize = 0x01020005;
constexpr std::uint32_t kMesh = 0x01040000, kNodeSection = 0x01041000, kNodeHeader = 0x01041100,
                        kNodeSize = 0x01041101, kNodeCoords = 0x01041200,
                        kDomainSection = 0x01042000, kDomain = 0x01042100, kDomainHdr = 0x01042101,
                        kDomElemType = 0x01042102, kDomPartId = 0x01042103, kDomName = 0x01032105,
                        kDomElemList = 0x01042200, kElement = 0x01042201,
                        kSurfaceSection = 0x01043000, kSurface = 0x01043100,
                        kSurfaceHdr = 0x01043101, kSurfaceName = 0x01043104, kFaceList = 0x01043200,
                        kFace = 0x01043201, kNodesetSection = 0x01044000, kNodeset = 0x01044100,
                        kNodesetHdr = 0x01044101, kNodesetName = 0x01044103,
                        kNodesetList = 0x01044200, kPartsSection = 0x01045000, kPart = 0x01045100,
                        kPartId = 0x01045101, kPartName = 0x01045102,
                        kElementsetSection = 0x01046000, kElementset = 0x01046100,
                        kElementsetHdr = 0x01046101, kElementsetName = 0x01046103,
                        kElementsetList = 0x01046200, kFacetsetSection = 0x01047000,
                        kFacetset = 0x01047100, kFacetsetHdr = 0x01047101,
                        kFacetsetName = 0x01047103, kFacetsetList = 0x01047200, kFacet = 0x01047201,
                        kEdgeSection = 0x01048000, kEdge = 0x01048100, kEdgeHdr = 0x01048101,
                        kEdgeName = 0x01048104, kEdgeList = 0x01048200, kLine = 0x01048201;
constexpr std::uint32_t kState = 0x02000000, kStateHeader = 0x02010000, kStateTime = 0x02010002,
                        kStateStatus = 0x02010003, kStateData = 0x02020000,
                        kStateVariable = 0x02020001, kStateVarId = 0x02020002,
                        kStateVarData = 0x02020003;

enum class XpltGroup { Global, Node, Domain, Surface, Edge, None };

XpltGroup xplt_dictionary_group(std::uint32_t Id) {
    switch (Id) {
        case 0x01021000:
            return XpltGroup::Global;
        case 0x01023000:
            return XpltGroup::Node;
        case 0x01024000:
            return XpltGroup::Domain;
        case 0x01025000:
            return XpltGroup::Surface;
        case 0x01026000:
            return XpltGroup::Edge;
        default:
            return XpltGroup::None;
    }
}

XpltGroup xplt_data_group(std::uint32_t Id) {
    switch (Id) {
        case 0x02020100:
            return XpltGroup::Global;
        case 0x02020300:
            return XpltGroup::Node;
        case 0x02020400:
            return XpltGroup::Domain;
        case 0x02020500:
            return XpltGroup::Surface;
        case 0x02020600:
            return XpltGroup::Edge;
        default:
            return XpltGroup::None;
    }
}

[[noreturn]] void xplt_fail(const std::string& rWhat) {
    throw ReadError("FEBio .xplt: " + rWhat);
}

// PLT element type -> meshio++ type (or a FEBio-only name) and file node count.
struct XpltElem {
    const char* mType;
    std::size_t mNodes;
};

constexpr XpltElem kXpltElems[] = {
    {"hexahedron", 8}, {"wedge", 6},         {"tetra", 4},         {"quad", 4},
    {"triangle", 3},   {"line", 2},          {"hexahedron20", 20}, {"tetra10", 10},
    {"tet15", 15},     {"hexahedron27", 27}, {"triangle6", 6},     {"quad8", 8},
    {"quad9", 9},      {"wedge15", 15},      {"tet20", 20},        {"triangle10", 10},
    {"pyramid", 5},    {"tet5", 5},          {"pyramid13", 13},    {"line3", 3},
};

// FEBio-only types and what `mLenient` keeps of them (nullptr: nothing).
const char* xplt_downgrade(std::string_view Type, bool& rLossy) {
    rLossy = true;
    if (Type == "tet15")
        return "tetra10";
    if (Type == "tet5")
        return "tetra";
    if (Type == "tet20")
        return nullptr;
    rLossy = false;
    return nullptr;
}

std::size_t xplt_corners(std::size_t Nodes) {
    switch (Nodes) {
        case 6:
        case 7:
        case 10:
            return 3;
        case 8:
        case 9:
            return 4;
        default:
            return Nodes;
    }
}

const char* xplt_facet_type(std::size_t Nodes) {
    switch (Nodes) {
        case 3:
            return "triangle";
        case 4:
            return "quad";
        case 6:
            return "triangle6";
        case 7:
            return "triangle7";
        case 8:
            return "quad8";
        case 9:
            return "quad9";
        case 10:
            return "triangle10";
        default:
            return nullptr;
    }
}

int xplt_dim(std::string_view Type) {
    if (Type.substr(0, 4) == "line")
        return 1;
    if (Type.substr(0, 8) == "triangle" || Type.substr(0, 4) == "quad")
        return 2;
    return 3;
}

// A byte range in either the mapped file or an inflated chunk, with the
// file's byte order.
struct XpltView {
    std::string_view mData;
    bool mSwap = false;

    std::uint32_t U32(std::size_t Pos) const {
        std::uint32_t v;
        std::memcpy(&v, mData.data() + Pos, 4);
        return mSwap ? detail::bswap32(v) : v;
    }
    std::int32_t I32(std::size_t Pos) const { return static_cast<std::int32_t>(U32(Pos)); }
    float F32(std::size_t Pos) const {
        const std::uint32_t u = U32(Pos);
        float f;
        std::memcpy(&f, &u, 4);
        return f;
    }

    // Calls rFn(id, begin, end) for each chunk in [Begin, End).
    template <class F>
    void Chunks(std::size_t Begin, std::size_t End, F&& rFn) const {
        std::size_t pos = Begin;
        while (pos + 8 <= End) {
            const std::uint32_t id = U32(pos);
            const std::uint32_t size = U32(pos + 4);
            if (pos + 8 + size > End)
                xplt_fail("a chunk runs past its parent (truncated file?)");
            if (!rFn(id, pos + 8, pos + 8 + static_cast<std::size_t>(size)))
                return;
            pos += 8 + static_cast<std::size_t>(size);
        }
    }

    std::optional<std::pair<std::size_t, std::size_t>> Child(std::size_t Begin, std::size_t End,
                                                             std::uint32_t Id) const {
        std::optional<std::pair<std::size_t, std::size_t>> out;
        Chunks(Begin, End, [&](std::uint32_t CId, std::size_t A, std::size_t B) {
            if (CId == Id) {
                out = std::make_pair(A, B);
                return false;
            }
            return true;
        });
        return out;
    }

    std::optional<std::int32_t> Int(std::size_t Begin, std::size_t End, std::uint32_t Id) const {
        const auto c = Child(Begin, End, Id);
        if (!c || c->second - c->first < 4)
            return std::nullopt;
        return I32(c->first);
    }

    // A length-prefixed string.
    std::string String(std::size_t Begin, std::size_t End) const {
        if (End - Begin < 4)
            return {};
        const auto n = static_cast<std::size_t>(std::max<std::int32_t>(I32(Begin), 0));
        std::string s(mData.substr(Begin + 4, std::min(n, End - Begin - 4)));
        while (!s.empty() && s.back() == '\0')
            s.pop_back();
        return s;
    }

    // A 64-byte NUL-padded string.
    std::string Fixed(std::size_t Begin, std::size_t End) const {
        std::string s(mData.substr(Begin, End - Begin));
        const std::size_t nul = s.find('\0');
        return nul == std::string::npos ? s : s.substr(0, nul);
    }
};

// One top-level chunk: where it lives in the file, and whether it is a zlib
// stream there.
struct XpltTop {
    std::uint32_t mId = 0;
    std::size_t mOffset = 0;  // file offset of the chunk (or of its stream)
    bool mCompressed = false;
    std::size_t mBegin = 0, mEnd = 0;  // payload range, file offsets when raw
};

struct XpltItem {
    std::string mName;
    int mFmt = 1;  // 0 node, 1 item, 2 mult, 3 region, 4 matpoints
    std::size_t mWidth = 1;
};

struct XpltDomain {
    std::string mType;
    std::int32_t mPart = -1;
    std::string mName;
    std::vector<std::int64_t> mIds;
    std::vector<std::int64_t> mConn;  // file (FEBio) order, `mNodes` per element
    std::size_t mNodes = 0;
};

struct XpltFile {
    detail::FileSource mSource;
    XpltView mRaw;
    std::uint32_t mVersion = 0;
    std::map<XpltGroup, std::vector<XpltItem>> mItems;
    std::vector<XpltTop> mStates;
    std::vector<float> mTimes;
    std::vector<std::optional<std::int32_t>> mStatus;
    // Every mesh section in file order (a remeshed run writes one before the
    // states that use it), and the mesh each state uses.
    std::vector<XpltTop> mMeshes;
    std::vector<std::size_t> mStateMesh;

    explicit XpltFile(const std::string& rPath, const ReadOptions& rOptions)
        : mSource(rPath, rOptions.mMmap) {
        const std::string_view data = mSource.View();
        if (data.size() < 4)
            xplt_fail("file too short");
        std::uint32_t magic;
        std::memcpy(&magic, data.data(), 4);
        mRaw.mData = data;
        if (magic == kXpltMagic)
            mRaw.mSwap = false;
        else if (detail::bswap32(magic) == kXpltMagic)
            mRaw.mSwap = true;
        else
            xplt_fail("not an FEBio plot file (bad magic)");

        std::size_t pos = 4;
        bool compressed = false;
        std::size_t n_top = 0;
        std::size_t n_mesh = 0;
        std::string inflated;
        while (pos < data.size()) {
            if (n_top >= 2 && compressed) {
                std::size_t used = 0;
                try {
                    inflated = detail::zlib_inflate(data.substr(pos), 15, &used, "FEBio .xplt");
                } catch (const ReadError&) {
                    if (!detail::zlib_available())
                        throw;
                    log::warn("FEBio .xplt: the last state is truncated and was not read");
                    break;
                }
                if (inflated.size() < 8)
                    break;
                const XpltView view{inflated, mRaw.mSwap};
                XpltTop top{view.U32(0), pos, true, 8, 8 + view.U32(4)};
                if (top.mId == kMesh) {
                    ++n_mesh;
                    mMeshes.push_back(top);
                } else if (top.mId == kState) {
                    AddState(view, top);
                }
                pos += used;
                ++n_top;
                continue;
            }
            if (pos + 8 > data.size()) {
                log::warn("FEBio .xplt: trailing bytes after the last state were ignored");
                break;
            }
            const std::uint32_t id = mRaw.U32(pos);
            const std::size_t size = mRaw.U32(pos + 4);
            if (pos + 8 + size > data.size()) {
                if (id == kState) {
                    log::warn("FEBio .xplt: the last state is truncated and was not read");
                    break;
                }
                xplt_fail("the file is truncated");
            }
            XpltTop top{id, pos, false, pos + 8, pos + 8 + size};
            if (n_top == 0) {
                if (id != kRoot)
                    xplt_fail("the file does not start with its root section");
                ReadRoot(top.mBegin, top.mEnd, compressed);
            } else if (id == kMesh) {
                ++n_mesh;
                mMeshes.push_back(top);
            } else if (id == kState) {
                AddState(mRaw, top);
            }
            pos += 8 + size;
            ++n_top;
        }
        if (n_mesh == 0)
            xplt_fail("the file has no mesh");
    }

    // A mesh section's payload, inflated into `rStore` when it is compressed.
    XpltView MeshView(std::size_t Index, std::string& rStore, std::size_t& rBegin,
                      std::size_t& rEnd) const {
        const XpltTop& top = mMeshes[Index];
        if (!top.mCompressed) {
            rBegin = top.mBegin;
            rEnd = top.mEnd;
            return mRaw;
        }
        rStore = detail::zlib_inflate(mRaw.mData.substr(top.mOffset), 15, nullptr, "FEBio .xplt");
        const XpltView view{rStore, mRaw.mSwap};
        rBegin = 8;
        rEnd = std::min<std::size_t>(8 + view.U32(4), rStore.size());
        return view;
    }

    void ReadRoot(std::size_t Begin, std::size_t End, bool& rCompressed) {
        const auto header = mRaw.Child(Begin, End, kHeader);
        if (header) {
            mVersion = static_cast<std::uint32_t>(
                mRaw.Int(header->first, header->second, kHdrVersion).value_or(0));
            rCompressed = mRaw.Int(header->first, header->second, kHdrCompression).value_or(0) != 0;
        }
        if (mVersion < 0x0030 || mVersion > 0x00FF) {
            char buf[16];
            std::snprintf(buf, sizeof(buf), "0x%04x", mVersion);
            xplt_fail(std::string("plot file version ") + buf +
                      " is not supported (FEBio 3 and 4 write 0x0030 and later)");
        }
        const auto dict = mRaw.Child(Begin, End, kDictionary);
        if (!dict)
            return;
        mRaw.Chunks(
            dict->first, dict->second, [&](std::uint32_t GId, std::size_t A, std::size_t B) {
                const XpltGroup group = xplt_dictionary_group(GId);
                if (group == XpltGroup::None)
                    return true;
                auto& items = mItems[group];
                mRaw.Chunks(A, B, [&](std::uint32_t IId, std::size_t C, std::size_t D) {
                    if (IId != kDicItem)
                        return true;
                    const std::int32_t type = mRaw.Int(C, D, kDicItemType).value_or(0);
                    const std::int32_t fmt = mRaw.Int(C, D, kDicItemFmt).value_or(0);
                    const std::int32_t asize = mRaw.Int(C, D, kDicItemArraySize).value_or(0);
                    const auto named = mRaw.Child(C, D, kDicItemName);
                    std::string name =
                        named ? mRaw.Fixed(named->first, named->second) : std::string();
                    if (const std::size_t eq = name.find('='); eq != std::string::npos)
                        name = name.substr(eq + 1);
                    std::size_t width = 1;
                    switch (type) {
                        case 1:
                            width = 3;
                            break;
                        case 2:
                            width = 6;
                            break;
                        case 3:
                            width = 3;
                            break;
                        case 4:
                            width = 21;
                            break;
                        case 5:
                            width = 9;
                            break;
                        case 6:
                            width = static_cast<std::size_t>(std::max(asize, 1));
                            break;
                        case 7:
                            width = 3 * static_cast<std::size_t>(std::max(asize, 1));
                            break;
                        default:
                            width = 1;
                    }
                    items.push_back(XpltItem{name, fmt >= 0 && fmt <= 4 ? fmt : 1, width});
                    return true;
                });
                return true;
            });
    }

    void AddState(const XpltView& rView, const XpltTop& rTop) {
        float time = 0.0f;
        std::optional<std::int32_t> status;
        if (const auto hdr = rView.Child(rTop.mBegin, rTop.mEnd, kStateHeader)) {
            if (const auto t = rView.Child(hdr->first, hdr->second, kStateTime))
                time = rView.F32(t->first);
            status = rView.Int(hdr->first, hdr->second, kStateStatus);
        }
        mStates.push_back(rTop);
        mTimes.push_back(time);
        mStatus.push_back(status);
        mStateMesh.push_back(mMeshes.empty() ? 0 : mMeshes.size() - 1);
    }
};

std::vector<std::pair<std::size_t, std::vector<std::int64_t>>> xplt_facets(const XpltView& rView,
                                                                           std::size_t Begin,
                                                                           std::size_t End,
                                                                           std::uint32_t ListTag,
                                                                           std::uint32_t LeafTag) {
    std::vector<std::pair<std::size_t, std::vector<std::int64_t>>> out;
    const auto list = rView.Child(Begin, End, ListTag);
    if (!list)
        return out;
    rView.Chunks(list->first, list->second, [&](std::uint32_t Id, std::size_t A, std::size_t B) {
        if (Id != LeafTag || B - A < 8)
            return true;
        const auto nn = static_cast<std::size_t>(std::max<std::int32_t>(rView.I32(A + 4), 0));
        std::vector<std::int64_t> nodes;
        for (std::size_t k = 0; k < nn && A + 8 + 4 * k + 4 <= B; ++k)
            nodes.push_back(rView.I32(A + 8 + 4 * k));
        out.emplace_back(nn, std::move(nodes));
        return true;
    });
    return out;
}

struct XpltSurface {
    std::string mName;
    std::vector<std::pair<std::size_t, std::vector<std::int64_t>>> mFacets;  // (nodes, ids)
    bool mFacetSet = false;  // a facet set (FEBio 4) rather than a data surface
};

struct XpltRaw {
    std::vector<double> mCoords;
    std::vector<XpltDomain> mDomains;
    std::vector<XpltSurface> mSurfaces;
    std::vector<XpltSurface> mEdges;  // their "facets" are line segments
    std::vector<std::pair<std::string, std::vector<std::int64_t>>> mNodeSets, mElemSets;
    std::map<std::int32_t, std::string> mParts;
};

XpltRaw xplt_read_mesh(const XpltView& rView, std::size_t Begin, std::size_t End) {
    XpltRaw raw;
    bool have_nodes = false;
    rView.Chunks(Begin, End, [&](std::uint32_t Sid, std::size_t A, std::size_t B) {
        if (Sid == kNodeSection) {
            const auto head = rView.Child(A, B, kNodeHeader);
            const auto n = static_cast<std::size_t>(
                head ? std::max(rView.Int(head->first, head->second, kNodeSize).value_or(0), 0)
                     : 0);
            const auto coords = rView.Child(A, B, kNodeCoords);
            if (!coords)
                xplt_fail("the mesh has no node coordinates");
            const std::size_t words = (coords->second - coords->first) / 4;
            const std::size_t stride = words == 4 * n ? 4 : words == 3 * n ? 3 : 0;
            if (stride == 0)
                xplt_fail("the node coordinates do not match the node count");
            raw.mCoords.resize(3 * n);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t d = 0; d < 3; ++d)
                    raw.mCoords[3 * i + d] =
                        rView.F32(coords->first + 4 * (stride * i + (stride - 3) + d));
            have_nodes = true;
        } else if (Sid == kDomainSection) {
            rView.Chunks(A, B, [&](std::uint32_t Did, std::size_t C, std::size_t D) {
                if (Did != kDomain)
                    return true;
                const auto hdr = rView.Child(C, D, kDomainHdr);
                XpltDomain dom;
                const std::int32_t etype =
                    hdr ? rView.Int(hdr->first, hdr->second, kDomElemType).value_or(-1) : -1;
                if (etype < 0 || etype >= static_cast<std::int32_t>(std::size(kXpltElems)))
                    xplt_fail("unknown element type " + std::to_string(etype));
                dom.mType = kXpltElems[etype].mType;
                dom.mNodes = kXpltElems[etype].mNodes;
                dom.mPart = rView.Int(hdr->first, hdr->second, kDomPartId).value_or(-1);
                if (const auto named = rView.Child(hdr->first, hdr->second, kDomName))
                    dom.mName = rView.String(named->first, named->second);
                if (const auto list = rView.Child(C, D, kDomElemList)) {
                    rView.Chunks(list->first, list->second,
                                 [&](std::uint32_t Eid, std::size_t E, std::size_t F) {
                                     if (Eid != kElement || F - E < 4 * (dom.mNodes + 1))
                                         return true;
                                     dom.mIds.push_back(rView.I32(E));
                                     for (std::size_t k = 0; k < dom.mNodes; ++k)
                                         dom.mConn.push_back(rView.I32(E + 4 * (k + 1)));
                                     return true;
                                 });
                }
                raw.mDomains.push_back(std::move(dom));
                return true;
            });
        } else if (Sid == kSurfaceSection || Sid == kFacetsetSection) {
            const bool facetset = Sid == kFacetsetSection;
            rView.Chunks(A, B, [&](std::uint32_t Id, std::size_t C, std::size_t D) {
                if (Id != (facetset ? kFacetset : kSurface))
                    return true;
                std::string name;
                if (const auto hdr = rView.Child(C, D, facetset ? kFacetsetHdr : kSurfaceHdr))
                    if (const auto named = rView.Child(hdr->first, hdr->second,
                                                       facetset ? kFacetsetName : kSurfaceName))
                        name = rView.String(named->first, named->second);
                raw.mSurfaces.push_back(
                    {name,
                     xplt_facets(rView, C, D, facetset ? kFacetsetList : kFaceList,
                                 facetset ? kFacet : kFace),
                     facetset});
                return true;
            });
        } else if (Sid == kEdgeSection) {
            // An edge (FEBio 3.5+): named line segments, each [edge id, nodes,
            // n1 .. n3], the carrier of edge variables.
            rView.Chunks(A, B, [&](std::uint32_t Id, std::size_t C, std::size_t D) {
                if (Id != kEdge)
                    return true;
                std::string name;
                if (const auto hdr = rView.Child(C, D, kEdgeHdr))
                    if (const auto named = rView.Child(hdr->first, hdr->second, kEdgeName))
                        name = rView.String(named->first, named->second);
                raw.mEdges.push_back({name, xplt_facets(rView, C, D, kEdgeList, kLine), false});
                return true;
            });
        } else if (Sid == kNodesetSection || Sid == kElementsetSection) {
            const bool nodes = Sid == kNodesetSection;
            rView.Chunks(A, B, [&](std::uint32_t Id, std::size_t C, std::size_t D) {
                if (Id != (nodes ? kNodeset : kElementset))
                    return true;
                std::string name;
                if (const auto hdr = rView.Child(C, D, nodes ? kNodesetHdr : kElementsetHdr))
                    if (const auto named = rView.Child(hdr->first, hdr->second,
                                                       nodes ? kNodesetName : kElementsetName))
                        name = rView.String(named->first, named->second);
                std::vector<std::int64_t> members;
                if (const auto list = rView.Child(C, D, nodes ? kNodesetList : kElementsetList))
                    for (std::size_t p = list->first; p + 4 <= list->second; p += 4)
                        members.push_back(rView.I32(p));
                (nodes ? raw.mNodeSets : raw.mElemSets).emplace_back(name, std::move(members));
                return true;
            });
        } else if (Sid == kPartsSection) {
            rView.Chunks(A, B, [&](std::uint32_t Id, std::size_t C, std::size_t D) {
                if (Id != kPart)
                    return true;
                const std::int32_t pid = rView.Int(C, D, kPartId).value_or(-1);
                const auto named = rView.Child(C, D, kPartName);
                raw.mParts[pid] = named ? rView.Fixed(named->first, named->second) : std::string();
                return true;
            });
        }
        return true;
    });
    if (!have_nodes)
        xplt_fail("the file has no mesh");
    return raw;
}

// A region under construction: (tag, dim, entries), keyed by (kind, name).
struct XpltGroupEntry {
    std::int64_t mTag;
    int mDim;
    std::vector<std::int64_t> mEntries;
};

}  // namespace

Mesh read_xplt(const std::string& rPath, const ReadOptions& rOptions) {
    const XpltFile file(rPath, rOptions);
    // The step first: a remeshed run reads the mesh that step uses.
    const std::size_t n_states = file.mStates.size();
    if (n_states == 0 && rOptions.mTimeStep != 0 && rOptions.mTimeStep != -1)
        xplt_fail("time step " + std::to_string(rOptions.mTimeStep) +
                  " is out of range: the file has no states");
    const std::size_t index = n_states == 0 ? 0 : rOptions.ResolveTimeStep(n_states);
    std::string mesh_store;
    std::size_t mesh_begin = 0, mesh_end = 0;
    const XpltView mesh_view =
        file.MeshView(n_states == 0 ? 0 : file.mStateMesh[index], mesh_store, mesh_begin, mesh_end);
    const XpltRaw raw = xplt_read_mesh(mesh_view, mesh_begin, mesh_end);

    // --- mesh -----------------------------------------------------------------------
    Mesh mesh;
    NDArray points(DType::Float64, {raw.mCoords.size() / 3, 3});
    std::copy(raw.mCoords.begin(), raw.mCoords.end(), points.As<double>());
    mesh.AssignPoints(std::move(points));
    std::map<std::pair<int, std::string>, XpltGroupEntry> groups;
    const auto group = [&](RegionKind Kind, const std::string& rName, std::int64_t Tag,
                           int Dim) -> XpltGroupEntry& {
        return groups.try_emplace({static_cast<int>(Kind), rName}, XpltGroupEntry{Tag, Dim, {}})
            .first->second;
    };
    std::unordered_map<std::int64_t, std::int64_t> elem_index;
    std::int64_t base = 0;
    std::size_t lossy = 0;
    for (std::size_t k = 0; k < raw.mDomains.size(); ++k) {
        const XpltDomain& dom = raw.mDomains[k];
        std::string type = dom.mType;
        bool is_lossy = false;
        const char* keep = xplt_downgrade(type, is_lossy);
        if (is_lossy) {
            if (!rOptions.mLenient || !keep)
                xplt_fail("element type " + type + " has no meshio++ cell type" +
                          std::string(keep ? " (read with lenient to downgrade it)" : ""));
            ++lossy;
            type = keep;
        }
        const std::size_t kept =
            static_cast<std::size_t>(cell_type_num_nodes(cell_type_from_name(type)));
        const detail::NodeOrder* order = detail::node_order("febio", type);
        const std::size_t rows = dom.mIds.size();
        NDArray conn(DType::Int64, {rows, kept});
        std::int64_t* c = conn.As<std::int64_t>();
        for (std::size_t r = 0; r < rows; ++r)
            for (std::size_t j = 0; j < kept; ++j)
                c[r * kept + j] =
                    dom.mConn[r * dom.mNodes +
                              (order ? static_cast<std::size_t>(order->mToMeshio[j]) : j)];
        mesh.AddCellBlock(type, std::move(conn));
        for (std::size_t r = 0; r < rows; ++r)
            elem_index.emplace(dom.mIds[r], base + static_cast<std::int64_t>(r));
        // FEBio 4 names only solid domains; a nameless one takes the name of the
        // element set holding exactly its elements, else its part's.
        std::string name = dom.mName;
        if (name.empty()) {
            const std::set<std::int64_t> mine(dom.mIds.begin(), dom.mIds.end());
            for (const auto& [set_name, members] : raw.mElemSets)
                if (std::set<std::int64_t>(members.begin(), members.end()) == mine) {
                    name = set_name;
                    break;
                }
        }
        if (name.empty()) {
            const auto part = raw.mParts.find(dom.mPart);
            if (part != raw.mParts.end())
                name = part->second;
        }
        if (name.empty())
            name = "Part" + std::to_string(k + 1);
        XpltGroupEntry& g = group(RegionKind::Cell, name, dom.mPart, xplt_dim(type));
        g.mDim = std::max(g.mDim, xplt_dim(type));
        for (std::size_t r = 0; r < rows; ++r)
            g.mEntries.push_back(base + static_cast<std::int64_t>(r));
        base += static_cast<std::int64_t>(rows);
    }
    if (lossy)
        log::warn("FEBio .xplt: {} domain(s) downgraded to a meshio++ cell type", lossy);

    for (const auto& [name, members] : raw.mNodeSets) {
        XpltGroupEntry& g = group(RegionKind::Point, name, -1, -1);
        g.mEntries.insert(g.mEntries.end(), members.begin(), members.end());
    }
    for (const auto& [name, members] : raw.mElemSets) {
        XpltGroupEntry& g = group(RegionKind::Cell, name, -1, -1);
        for (std::int64_t id : members) {
            const auto it = elem_index.find(id);
            if (it != elem_index.end())
                g.mEntries.push_back(it->second);
        }
    }

    // Surfaces: a data surface is a block of facet cells per facet type (its
    // variables live there), plus a side region when its facets lie on the
    // cells; a facet set is only the side region, or facet blocks when some
    // facet lies on no cell.
    std::optional<detail::FacetIndex> faces;
    std::vector<std::vector<std::int64_t>> surface_cells;  // per data surface, per facet
    std::vector<std::pair<std::string, std::vector<std::int64_t>>> blocks_to_add;  // (type, rows)
    std::vector<std::int64_t> block_surface;  // data surface id (1-based) of each added block
    std::set<std::string> seen_sides;
    std::int64_t surface_id = 0;
    for (const XpltSurface& surf : raw.mSurfaces) {
        const std::string& name = surf.mName;
        if (!surf.mFacetSet)
            ++surface_id;
        if (surf.mFacets.empty()) {
            if (!surf.mFacetSet)
                surface_cells.emplace_back();
            continue;
        }
        if (!faces) {
            detail::FacetIndexOptions options;
            options.mSurfaceEdges = false;
            faces.emplace(mesh, options);
        }
        std::vector<std::int64_t> sides;
        bool all = true;
        for (const auto& [nn, nodes] : surf.mFacets) {
            const std::size_t corners = std::min(xplt_corners(nn), nodes.size());
            const detail::FacetHit* hit = faces->Find(nodes.data(), corners);
            if (!hit) {
                all = false;
                break;
            }
            sides.push_back(hit->mFirst.mCell);
            sides.push_back(hit->mFirst.mFacet);
        }
        if (all && seen_sides.insert(name).second) {
            XpltGroupEntry& g = group(RegionKind::Side, name, -1, 2);
            g.mEntries.insert(g.mEntries.end(), sides.begin(), sides.end());
        }
        if (all && surf.mFacetSet)
            continue;
        // Facet blocks, one per facet type in first-seen order; each facet's
        // global cell is recorded in facet order for the surface's data.
        std::vector<std::string> order;
        std::map<std::string, std::vector<std::size_t>> by_type;  // type -> facet indices
        for (std::size_t f = 0; f < surf.mFacets.size(); ++f) {
            const auto& [nn, nodes] = surf.mFacets[f];
            const char* type = xplt_facet_type(nn);
            if (!type || nodes.size() != nn)
                continue;
            if (!by_type.count(type))
                order.push_back(type);
            by_type[type].push_back(f);
        }
        std::vector<std::int64_t> cells(surf.mFacets.size(), -1);
        XpltGroupEntry& g = group(RegionKind::Cell, name, -1, 2);
        for (const std::string& type : order) {
            std::vector<std::int64_t> rows;
            for (std::size_t f : by_type[type]) {
                const auto& nodes = surf.mFacets[f].second;
                rows.insert(rows.end(), nodes.begin(), nodes.end());
                cells[f] = base;
                g.mEntries.push_back(base);
                ++base;
            }
            blocks_to_add.emplace_back(type, std::move(rows));
            block_surface.push_back(surf.mFacetSet ? 0 : surface_id);
        }
        if (!surf.mFacetSet)
            surface_cells.push_back(std::move(cells));
    }
    // Edges: a block of line cells per segment type, a cell region each, and
    // each segment's global cell for the edge's data.
    std::vector<std::vector<std::int64_t>> edge_cells;  // per edge, per segment
    std::vector<std::int64_t> block_edge(blocks_to_add.size(), 0);
    for (std::size_t e = 0; e < raw.mEdges.size(); ++e) {
        const XpltSurface& edge = raw.mEdges[e];
        std::vector<std::string> order;
        std::map<std::string, std::vector<std::size_t>> by_type;
        for (std::size_t f = 0; f < edge.mFacets.size(); ++f) {
            const auto& [nn, nodes] = edge.mFacets[f];
            const char* type = nn == 2 ? "line" : nn == 3 ? "line3" : nullptr;
            if (!type || nodes.size() != nn)
                continue;
            if (!by_type.count(type))
                order.push_back(type);
            by_type[type].push_back(f);
        }
        std::vector<std::int64_t> cells(edge.mFacets.size(), -1);
        XpltGroupEntry& g = group(RegionKind::Cell, edge.mName, -1, 1);
        g.mDim = std::max(g.mDim, 1);
        for (const std::string& type : order) {
            std::vector<std::int64_t> rows;
            for (std::size_t f : by_type[type]) {
                const auto& nodes = edge.mFacets[f].second;
                rows.insert(rows.end(), nodes.begin(), nodes.end());
                cells[f] = base;
                g.mEntries.push_back(base);
                ++base;
            }
            blocks_to_add.emplace_back(type, std::move(rows));
            block_surface.push_back(0);
            block_edge.push_back(static_cast<std::int64_t>(e + 1));
        }
        edge_cells.push_back(std::move(cells));
    }
    const std::size_t n_domain_blocks = mesh.NumCellBlocks();
    for (const auto& [type, flat] : blocks_to_add) {
        const std::size_t k =
            static_cast<std::size_t>(cell_type_num_nodes(cell_type_from_name(type)));
        NDArray conn(DType::Int64, {flat.size() / k, k});
        std::copy(flat.begin(), flat.end(), conn.As<std::int64_t>());
        mesh.AddCellBlock(type, std::move(conn));
    }
    if (rOptions.WantsAnyData() && rOptions.WantsArray("xplt:surface") &&
        std::any_of(block_surface.begin(), block_surface.end(),
                    [](std::int64_t v) { return v > 0; })) {
        std::vector<NDArray> ids;
        for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b) {
            const std::size_t n = mesh.Cells(b).NumCells();
            NDArray a(DType::Int64, {n});
            const std::int64_t v = b < n_domain_blocks ? 0 : block_surface[b - n_domain_blocks];
            std::fill(a.As<std::int64_t>(), a.As<std::int64_t>() + n, v);
            ids.push_back(std::move(a));
        }
        mesh.AddCellData("xplt:surface", std::move(ids));
    }
    if (rOptions.WantsAnyData() && rOptions.WantsArray("xplt:edge") &&
        std::any_of(block_edge.begin(), block_edge.end(), [](std::int64_t v) { return v > 0; })) {
        std::vector<NDArray> ids;
        for (std::size_t b = 0; b < mesh.NumCellBlocks(); ++b) {
            const std::size_t n = mesh.Cells(b).NumCells();
            NDArray a(DType::Int64, {n});
            const std::int64_t v = b < n_domain_blocks ? 0 : block_edge[b - n_domain_blocks];
            std::fill(a.As<std::int64_t>(), a.As<std::int64_t>() + n, v);
            ids.push_back(std::move(a));
        }
        mesh.AddCellData("xplt:edge", std::move(ids));
    }
    for (auto& [key, g] : groups) {
        const RegionKind kind = static_cast<RegionKind>(key.first);
        const bool side = kind == RegionKind::Side;
        NDArray entries = side ? NDArray(DType::Int64, {g.mEntries.size() / 2, 2})
                               : NDArray(DType::Int64, {g.mEntries.size()});
        std::copy(g.mEntries.begin(), g.mEntries.end(), entries.As<std::int64_t>());
        mesh.AddRegion(Region(key.second, kind, g.mDim, g.mTag, std::move(entries)));
    }

    // --- the chosen state --------------------------------------------------------------
    if (n_states == 0)
        return mesh;
    const auto scalar = [](double V, DType T) {
        NDArray a(T, {1});
        if (T == DType::Float64)
            a.As<double>()[0] = V;
        else
            a.As<std::int64_t>()[0] = static_cast<std::int64_t>(V);
        return a;
    };
    mesh.AddFieldData("meshio:time", scalar(file.mTimes[index], DType::Float64));
    mesh.AddFieldData("xplt:step", scalar(static_cast<double>(index), DType::Int64));
    if (file.mStatus[index])
        mesh.AddFieldData("xplt:status", scalar(*file.mStatus[index], DType::Int64));
    if (rOptions.mPointsOnly)
        return mesh;

    const XpltTop& top = file.mStates[index];
    std::string inflated;
    XpltView view = file.mRaw;
    std::size_t begin = top.mBegin, end = top.mEnd;
    if (top.mCompressed) {
        inflated =
            detail::zlib_inflate(file.mRaw.mData.substr(top.mOffset), 15, nullptr, "FEBio .xplt");
        view = XpltView{inflated, file.mRaw.mSwap};
        begin = 8;
        end = std::min<std::size_t>(8 + view.U32(4), inflated.size());
    }
    const auto data = view.Child(begin, end, kStateData);
    if (!data)
        return mesh;

    const std::size_t n_points = mesh.NumPoints();
    std::vector<std::size_t> sizes;
    for (const auto cb : mesh.CellRange())
        sizes.push_back(cb.NumCells());
    std::vector<std::string> skipped;
    const auto skip = [&](const std::string& rName) {
        if (std::find(skipped.begin(), skipped.end(), rName) == skipped.end())
            skipped.push_back(rName);
    };
    // Element-node variables: (sum, count) per point, in first-seen order.
    std::vector<std::string> sum_order;
    std::map<std::string, std::pair<std::vector<double>, std::vector<double>>> sums;
    std::map<std::string, std::size_t> sum_width;
    const auto f32s = [&](std::size_t A, std::size_t B) {
        std::vector<double> v((B - A) / 4);
        for (std::size_t i = 0; i < v.size(); ++i)
            v[i] = view.F32(A + 4 * i);
        return v;
    };
    view.Chunks(data->first, data->second, [&](std::uint32_t Gid, std::size_t C, std::size_t D) {
        const XpltGroup grp = xplt_data_group(Gid);
        if (grp == XpltGroup::None)
            return true;
        const auto items_it = file.mItems.find(grp);
        view.Chunks(C, D, [&](std::uint32_t Vid, std::size_t E, std::size_t F) {
            if (Vid != kStateVariable || items_it == file.mItems.end())
                return true;
            const std::int32_t var = view.Int(E, F, kStateVarId).value_or(0);
            if (var < 1 || static_cast<std::size_t>(var) > items_it->second.size())
                return true;
            const XpltItem& item = items_it->second[static_cast<std::size_t>(var - 1)];
            if (!rOptions.WantsArray(item.mName))
                return true;
            const auto vd = view.Child(E, F, kStateVarData);
            if (!vd)
                return true;
            std::vector<std::pair<std::uint32_t, std::vector<double>>> regions;
            view.Chunks(vd->first, vd->second,
                        [&](std::uint32_t Rid, std::size_t G, std::size_t H) {
                            regions.emplace_back(Rid, f32s(G, H));
                            return true;
                        });
            const std::size_t width = item.mWidth;
            const auto make = [&](std::size_t Rows, const double* pValues) {
                NDArray a = width == 1 ? NDArray(DType::Float64, {Rows})
                                       : NDArray(DType::Float64, {Rows, width});
                std::copy(pValues, pValues + Rows * width, a.As<double>());
                return a;
            };
            if (grp == XpltGroup::Surface || grp == XpltGroup::Edge) {
                // Region k is data surface (or edge) k: per facet or segment
                // (item), one value (region) -> its cells; per node (node, in
                // first-seen order over its facets or segments) or per facet
                // node (mult) -> averaged at the points, as element-node
                // values are.
                const bool on_edge = grp == XpltGroup::Edge;
                const auto& carrier_cells = on_edge ? edge_cells : surface_cells;
                if (item.mFmt == 1 || item.mFmt == 3) {
                    std::vector<std::vector<double>> blocks;
                    for (std::size_t n : sizes)
                        blocks.emplace_back(n * width, std::numeric_limits<double>::quiet_NaN());
                    const auto bases = detail::block_bases(mesh);
                    bool landed = false;
                    for (const auto& [rid, values] : regions) {
                        const std::size_t k = static_cast<std::size_t>(rid) - 1;
                        if (rid < 1 || k >= carrier_cells.size())
                            continue;
                        const auto& cells = carrier_cells[k];
                        if (values.size() < (item.mFmt == 3 ? width : cells.size() * width))
                            continue;
                        for (std::size_t f = 0; f < cells.size(); ++f) {
                            if (cells[f] < 0)
                                continue;
                            landed = true;
                            const std::size_t b = static_cast<std::size_t>(
                                std::upper_bound(bases.begin(), bases.end(), cells[f]) -
                                bases.begin() - 1);
                            const std::size_t r = static_cast<std::size_t>(cells[f] - bases[b]);
                            for (std::size_t w = 0; w < width; ++w)
                                blocks[b][r * width + w] =
                                    item.mFmt == 3 ? values[w] : values[f * width + w];
                        }
                    }
                    // A variable on no surface (edge) of this mesh gives no array.
                    if (!landed)
                        return true;
                    std::vector<NDArray> arrays;
                    for (std::size_t b = 0; b < sizes.size(); ++b)
                        arrays.push_back(make(sizes[b], blocks[b].data()));
                    mesh.AddCellData(item.mName, std::move(arrays));
                    return true;
                }
                if (item.mFmt != 0 && item.mFmt != 2) {
                    skip(item.mName);
                    return true;
                }
                const auto sums_of = [&]() -> auto& {
                    auto [it, fresh] =
                        sums.try_emplace(item.mName, std::vector<double>(n_points * width, 0.0),
                                         std::vector<double>(n_points, 0.0));
                    if (fresh) {
                        sum_order.push_back(item.mName);
                        sum_width[item.mName] = width;
                    }
                    return it->second;
                };
                std::size_t data_surface = 0;
                for (const XpltSurface& surf : on_edge ? raw.mEdges : raw.mSurfaces) {
                    if (surf.mFacetSet)
                        continue;
                    ++data_surface;
                    for (const auto& [rid, values] : regions) {
                        if (static_cast<std::size_t>(rid) != data_surface)
                            continue;
                        std::vector<std::int64_t> nodes;
                        std::unordered_set<std::int64_t> seen_node;
                        for (const auto& facet : surf.mFacets)
                            for (std::int64_t p : facet.second)
                                if (item.mFmt == 2 || seen_node.insert(p).second)
                                    nodes.push_back(p);
                        if (values.size() < nodes.size() * width || nodes.empty())
                            continue;
                        auto& [total, count] = sums_of();
                        for (std::size_t i = 0; i < nodes.size(); ++i) {
                            const auto p = static_cast<std::size_t>(nodes[i]);
                            if (p >= n_points)
                                continue;
                            for (std::size_t w = 0; w < width; ++w)
                                total[p * width + w] += values[i * width + w];
                            count[p] += 1.0;
                        }
                    }
                }
            } else if (grp == XpltGroup::Global) {
                if (!regions.empty()) {
                    const auto& v = regions.front().second;
                    NDArray a(DType::Float64, {v.size()});
                    std::copy(v.begin(), v.end(), a.As<double>());
                    mesh.AddFieldData(item.mName, std::move(a));
                }
            } else if (grp == XpltGroup::Node) {
                if (!regions.empty() && regions.front().second.size() >= n_points * width)
                    mesh.AddPointData(item.mName, make(n_points, regions.front().second.data()));
            } else if (item.mFmt == 1 || item.mFmt == 3) {
                std::vector<std::vector<double>> blocks;
                for (std::size_t n : sizes)
                    blocks.emplace_back(n * width, std::numeric_limits<double>::quiet_NaN());
                for (const auto& [rid, values] : regions) {
                    const std::size_t k = static_cast<std::size_t>(rid) - 1;
                    if (rid < 1 || k >= raw.mDomains.size())
                        continue;
                    const std::size_t ne = raw.mDomains[k].mIds.size();
                    if (values.size() < (item.mFmt == 3 ? width : ne * width))
                        continue;
                    for (std::size_t r = 0; r < ne; ++r)
                        for (std::size_t w = 0; w < width; ++w)
                            blocks[k][r * width + w] =
                                item.mFmt == 3 ? values[w] : values[r * width + w];
                }
                std::vector<NDArray> arrays;
                for (std::size_t b = 0; b < sizes.size(); ++b)
                    arrays.push_back(make(sizes[b], blocks[b].data()));
                mesh.AddCellData(item.mName, std::move(arrays));
            } else if (item.mFmt == 0 || item.mFmt == 2) {
                auto [it, fresh] =
                    sums.try_emplace(item.mName, std::vector<double>(n_points * width, 0.0),
                                     std::vector<double>(n_points, 0.0));
                if (fresh) {
                    sum_order.push_back(item.mName);
                    sum_width[item.mName] = width;
                }
                auto& [total, count] = it->second;
                for (const auto& [rid, values] : regions) {
                    const std::size_t k = static_cast<std::size_t>(rid) - 1;
                    if (rid < 1 || k >= raw.mDomains.size())
                        continue;
                    const auto& conn = raw.mDomains[k].mConn;
                    std::vector<std::int64_t> nodes;
                    if (item.mFmt == 0) {
                        std::unordered_map<std::int64_t, bool> seen_node;
                        for (std::int64_t p : conn)
                            if (seen_node.emplace(p, true).second)
                                nodes.push_back(p);
                    } else {
                        nodes = conn;
                    }
                    if (values.size() < nodes.size() * width)
                        continue;
                    for (std::size_t i = 0; i < nodes.size(); ++i) {
                        const auto p = static_cast<std::size_t>(nodes[i]);
                        if (p >= n_points)
                            continue;
                        for (std::size_t w = 0; w < width; ++w)
                            total[p * width + w] += values[i * width + w];
                        count[p] += 1.0;
                    }
                }
            } else {
                skip(item.mName);
            }
            return true;
        });
        return true;
    });
    for (const std::string& name : sum_order) {
        if (mesh.HasPointData(name)) {
            log::warn("FEBio .xplt: nodal '{}' kept; the element-node '{}' was dropped", name,
                      name);
            continue;
        }
        const std::size_t width = sum_width[name];
        auto& [total, count] = sums[name];
        for (std::size_t p = 0; p < n_points; ++p)
            for (std::size_t w = 0; w < width; ++w)
                total[p * width + w] = count[p] > 0 ? total[p * width + w] / count[p]
                                                    : std::numeric_limits<double>::quiet_NaN();
        NDArray a = width == 1 ? NDArray(DType::Float64, {n_points})
                               : NDArray(DType::Float64, {n_points, width});
        std::copy(total.begin(), total.end(), a.As<double>());
        mesh.AddPointData(name, std::move(a));
    }
    if (!skipped.empty()) {
        std::string list;
        for (const std::string& s : skipped)
            list += (list.empty() ? "" : ", ") + s;
        log::warn("FEBio .xplt: material-point variables are not read: {}", list);
    }
    return mesh;
}

MeshMetadata read_xplt_metadata(const std::string& rPath, const ReadOptions& rOptions) {
    ReadOptions options = rOptions;
    options.mPointsOnly = true;
    options.mTimeStep = 0;
    MeshMetadata meta = metadata_from_mesh(read_xplt(rPath, options));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "xplt";
    const XpltFile file(rPath, rOptions);
    meta.mTimeValues.assign(file.mTimes.begin(), file.mTimes.end());
    return meta;
}

}  // namespace meshioplusplus
