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
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes (private, not installed)
#include "xdmf_doc.hpp"
#include "xdmf_sets.hpp"

// Project includes
#include "meshioplusplus/formats/xdmf.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/xdmf_common.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"

#include "../detail/text_cursor.hpp"
#include "../detail/typed_view.hpp"

#ifdef MESHIOPLUSPLUS_HAS_HDF5
#include "meshioplusplus/detail/hdf5_util.hpp"
#endif

namespace fs = std::filesystem;

namespace meshioplusplus {

// The document-structure helpers live in the private `formats/xdmf_doc.hpp` so
// the transient writer's append path resolves a document exactly the way this
// reader does. They were local to this file through v9.1.0, which is how the
// writer came to carry a weaker transcription. Pulled in by name rather than
// qualified at each of the five call sites, so those stay byte-identical.
using xdmfdetail::xdmf_resolve;
using xdmfdetail::XdmfDoc;
// The reader has always spelled this one `parse_dims`; keep that at the call
// sites rather than churn them.
constexpr auto& parse_dims = xdmfdetail::xdmf_parse_dims;

namespace {

// ---- type maps (shared with HMF via detail/xdmf_common.hpp) ----

using xdmfcommon::attribute_type;
using xdmfcommon::concat_cell_data;
using xdmfcommon::dims_string;
using xdmfcommon::meshio_to_xdmf;
using xdmfcommon::numpy_to_xdmf_dtype;
using xdmfcommon::pack_mixed_topology;
using xdmfcommon::split_raw_cell_data;
using xdmfcommon::xdmf_to_meshio;

std::string xdmf_idx_to_meshio(int idx) {
    static const std::unordered_map<int, std::string> m = {
        {0x1, "vertex"},        {0x2, "line"},          {0x4, "triangle"},     {0x5, "quad"},
        {0x6, "tetra"},         {0x7, "pyramid"},       {0x8, "wedge"},        {0x9, "hexahedron"},
        {0x22, "line3"},        {0x23, "quad9"},        {0x24, "triangle6"},   {0x25, "quad8"},
        {0x26, "tetra10"},      {0x27, "pyramid13"},    {0x28, "wedge15"},     {0x29, "wedge18"},
        {0x30, "hexahedron20"}, {0x31, "hexahedron24"}, {0x32, "hexahedron27"}};
    auto it = m.find(idx);
    if (it == m.end())
        throw ReadError("XDMF: unknown mixed topology index");
    return it->second;
}

int xdmf_idx_num_nodes(int idx) {
    static const std::unordered_map<int, int> m = {
        {1, 1},     {2, 2},     {4, 3},     {5, 4},     {6, 4},     {7, 5},    {8, 6},
        {9, 8},     {11, 6},    {0x22, 3},  {0x23, 9},  {0x24, 6},  {0x25, 8}, {0x26, 10},
        {0x27, 13}, {0x28, 15}, {0x29, 18}, {0x30, 20}, {0x31, 24}, {0x32, 27}};
    auto it = m.find(idx);
    if (it == m.end())
        throw ReadError("XDMF: unknown mixed topology index");
    return it->second;
}

DType xdmf_to_dtype(const std::string& rDataType, const std::string& rPrecision) {
    int p = std::atoi(rPrecision.c_str());
    if (rDataType == "Int")
        return p == 1 ? DType::Int8 : p == 2 ? DType::Int16 : p == 4 ? DType::Int32 : DType::Int64;
    if (rDataType == "UInt")
        return p == 1   ? DType::UInt8
               : p == 2 ? DType::UInt16
               : p == 4 ? DType::UInt32
                        : DType::UInt64;
    return p == 4 ? DType::Float32 : DType::Float64;
}

// References are document-local. Resolve iteratively so cycles never recurse
// into the parser, and share this with metadata reads (which need the target's
// Dimensions rather than the reference node's usually absent attributes).
pugi::xml_node xdmf_resolve_data_item(pugi::xml_node node) {
    std::unordered_set<const void*> visited;
    while (node.attribute("Reference")) {
        if (!visited.insert(node.internal_object()).second)
            throw ReadError("XDMF: cyclic DataItem reference");
        const std::string ref = node.attribute("Reference").value();
        std::string xpath = ref == "XML" ? node.text().get() : ref;
        const std::size_t first = xpath.find_first_not_of(" \t\r\n");
        const std::size_t last = xpath.find_last_not_of(" \t\r\n");
        xpath = first == std::string::npos ? "" : xpath.substr(first, last - first + 1);
        if (xpath.empty() || xpath[0] != '/')
            throw ReadError("XDMF: DataItem reference must be an absolute XPath");
        try {
            const auto targets = node.root().select_nodes(xpath.c_str());
            if (targets.size() != 1 || std::string(targets[0].node().name()) != "DataItem")
                throw ReadError("XDMF: reference must select exactly one DataItem: " + xpath);
            node = targets[0].node();
        } catch (const pugi::xpath_exception& exc) {
            throw ReadError("XDMF: invalid reference XPath '" + xpath + "': " + exc.what());
        }
    }
    if (!node || std::string(node.name()) != "DataItem")
        throw ReadError("XDMF: missing DataItem");
    return node;
}

NDArray read_data_item(const pugi::xml_node& rItem, const fs::path& rBaseDir) {
    const pugi::xml_node rDi = xdmf_resolve_data_item(rItem);
    std::vector<std::size_t> dims = parse_dims(rDi.attribute("Dimensions").value());

    std::string data_type = "Float";
    if (rDi.attribute("DataType") && rDi.attribute("NumberType"))
        throw ReadError("XDMF: DataItem has both DataType and NumberType");
    if (rDi.attribute("DataType"))
        data_type = rDi.attribute("DataType").value();
    else if (rDi.attribute("NumberType"))
        data_type = rDi.attribute("NumberType").value();
    std::string precision = rDi.attribute("Precision") ? rDi.attribute("Precision").value() : "4";
    std::string fmt = rDi.attribute("Format").value();
    DType dt = xdmf_to_dtype(data_type, precision);

    std::size_t total = dims.empty() ? 0
                                     : std::accumulate(dims.begin(), dims.end(), std::size_t{1},
                                                       std::multiplies<>());

    if (fmt == "XML") {
        // The element's text read in place, token by token, and parsed straight
        // into the typed buffer with the dtype switch taken once (roadmap §3):
        // the lenient parse_double / strtoll / strtoull store_token used. A
        // short item leaves the rest zero, as the zero-filled array did.
        NDArray a = NDArray::Uninit(dt, dims);
        const std::string_view text = rDi.text().get();
        detail::dispatch_dtype(dt, [&]<class T>() {
            T* out = a.As<T>();
            std::size_t i = 0, pos = 0;
            while (i < total) {
                while (pos < text.size() && detail::text_is_blank(text[pos]))
                    ++pos;
                if (pos >= text.size())
                    break;
                const std::size_t b = pos;
                while (pos < text.size() && !detail::text_is_blank(text[pos]))
                    ++pos;
                const std::string_view tok = text.substr(b, pos - b);
                if constexpr (std::is_floating_point_v<T>)
                    out[i++] = static_cast<T>(detail::parse_double_prefix(tok));
                else if constexpr (std::is_signed_v<T>)
                    out[i++] = static_cast<T>(detail::strtoll_token(tok));
                else
                    out[i++] = static_cast<T>(detail::strtoull_token(tok));
            }
            if (i * sizeof(T) < a.Nbytes())
                std::memset(reinterpret_cast<unsigned char*>(out) + i * sizeof(T), 0,
                            a.Nbytes() - i * sizeof(T));
        });
        return a;
    }
    if (fmt == "Binary") {
        std::string rel = rDi.text().get();
        // trim whitespace
        std::size_t a0 = rel.find_first_not_of(" \t\r\n");
        std::size_t a1 = rel.find_last_not_of(" \t\r\n");
        std::string path = (a0 == std::string::npos) ? "" : rel.substr(a0, a1 - a0 + 1);
        auto bin = detail::make_classic_ifstream(path, std::ios::binary);
        if (!bin) {  // try relative to the xdmf file
            bin.open((rBaseDir / path).string(), std::ios::binary);
            if (!bin)
                throw ReadError("XDMF: could not open binary file " + path);
        }
        NDArray a(dt, dims);
        bin.read(reinterpret_cast<char*>(a.Data()), static_cast<std::streamsize>(a.Nbytes()));
        return a;
    }
    if (fmt != "HDF")
        throw ReadError("XDMF: unknown data format " + fmt);

#ifdef MESHIOPLUSPLUS_HAS_HDF5
    // "<file>.h5:/path/to/dataset", file path relative to the xdmf file.
    std::string info = rDi.text().get();
    std::size_t a0 = info.find_first_not_of(" \t\r\n");
    std::size_t a1 = info.find_last_not_of(" \t\r\n");
    info = (a0 == std::string::npos) ? "" : info.substr(a0, a1 - a0 + 1);
    std::size_t colon = info.find(':');
    if (colon == std::string::npos)
        throw ReadError("XDMF: malformed HDF reference '" + info + "'");
    std::string h5file = info.substr(0, colon);
    std::string h5path = info.substr(colon + 1);

    h5::SilenceErrors silence;
    fs::path full = rBaseDir / h5file;
    h5::Hid f = h5::open_file_read(full.string());
    NDArray a = h5::read_dataset(f, h5path);
    a.Reshape(dims);  // stored shape is authoritative in the XML
    return a;
#else
    throw ReadError("XDMF: HDF data format handled by Python fallback");
#endif
}

// Mixed-topology translation (ported from common.translate_mixed_cells).
// Appends one cell block per run of consecutive equal types onto `rMesh`.
void translate_mixed(const NDArray& rFlat, Mesh& rMesh) {
    std::size_t n = rFlat.Size();
    std::vector<int> types;
    std::vector<std::size_t> offsets;
    std::size_t r = 0;
    const detail::Int64View flat_values(rFlat);
    while (r < n) {
        int xt = static_cast<int>(flat_values[r]);
        types.push_back(xt);
        offsets.push_back(r);
        const auto nn = static_cast<std::size_t>(xdmf_idx_num_nodes(xt));
        // Polyvertex (1) and Polyline (2) carry their node count after the
        // type -- the writers emit it for both; reading it only for lines
        // misparsed every mixed topology holding a vertex.
        const std::size_t head = (xt == 1 || xt == 2) ? 2 : 1;
        if (r + head + nn > n)
            throw ReadError("XDMF: mixed topology ends inside a cell");
        if (head == 2 && static_cast<std::size_t>(flat_values[r + 1]) != nn)
            throw ReadError(xt == 1 ? "XDMF: only 1-point polyvertices supported"
                                    : "XDMF: only 2-point lines supported");
        r += head + nn;
    }
    // group consecutive equal types
    std::size_t start = 0;
    while (start < types.size()) {
        std::size_t end = start + 1;
        while (end < types.size() && types[end] == types[start])
            ++end;
        int xt = types[start];
        int nn = xdmf_idx_num_nodes(xt);
        std::size_t nrows = end - start;
        NDArray data = NDArray::Uninit(DType::Int64, {nrows, static_cast<std::size_t>(nn)});
        std::int64_t* dp = data.As<std::int64_t>();
        const std::size_t head = (xt == 1 || xt == 2) ? 2 : 1;
        // Rows are independent: copy them in parallel, one dtype switch per run.
        detail::dispatch_dtype(rFlat.Dtype(), [&]<class T>() {
            const T* src = rFlat.As<T>();
            parallel_for_bw(nrows, [&](std::size_t b) {
                const std::size_t base = offsets[start + b] + head;
                for (int j = 0; j < nn; ++j) {
                    if constexpr (std::is_floating_point_v<T>)
                        dp[b * nn + j] = detail::typed_view_int<T>(src[base + j]);
                    else
                        dp[b * nn + j] = static_cast<std::int64_t>(src[base + j]);
                }
            });
        });
        rMesh.AddCellBlock(xdmf_idx_to_meshio(xt), std::move(data));
        start = end;
    }
}

/** @brief Read a grid's `<Topology>`/`<Geometry>` into @p rMesh, ignoring the rest. */
void xdmf_read_geometry(const pugi::xml_node& rGrid, const fs::path& rBaseDir, Mesh& rMesh) {
    for (pugi::xml_node c : rGrid.children()) {
        const std::string tag = c.name();
        if (tag == "Topology") {
            std::string ctype = c.attribute("Type") ? c.attribute("Type").value()
                                                    : c.attribute("TopologyType").value();
            NDArray data = read_data_item(c.child("DataItem"), rBaseDir);
            if (ctype == "Mixed")
                translate_mixed(data, rMesh);
            else
                rMesh.AddCellBlock(xdmf_to_meshio(ctype), std::move(data));
        } else if (tag == "Geometry") {
            rMesh.AssignPoints(read_data_item(c.child("DataItem"), rBaseDir));
        }
    }
}

/** @brief Attach staged point data and split staged raw cell data onto @p rMesh. */
void xdmf_attach_data(Mesh& rMesh, std::vector<std::pair<std::string, NDArray>>& rPointData,
                      std::vector<std::pair<std::string, NDArray>>& rCellDataRaw) {
    for (auto& kv : rPointData)
        rMesh.AddPointData(kv.first, std::move(kv.second));

    std::vector<std::size_t> sizes;
    for (const auto cb : rMesh.CellRange())
        sizes.push_back(cb.NumCells());
    for (auto& kv : rCellDataRaw)
        rMesh.AddCellData(kv.first, split_raw_cell_data(kv.second, sizes));
}

/** @brief The `<Time Value=...>` of a step grid, or 0 when it carries none. */
double xdmf_step_time(const pugi::xml_node& rStep) {
    pugi::xml_node t = rStep.child("Time");
    return t ? t.attribute("Value").as_double(0.0) : 0.0;
}

// ---------------------------------------------------------------------------
// <Set> <-> regions (v16.27.0)
//
// SetType Node -> a Point region, Cell -> a Cell region, Face/Edge -> a Side
// region: the first DataItem holds the cell indices and the second the
// cell-local face or edge indices (the XDMF model's own layout), numbered as
// meshio++ numbers facets (`cell_faces`/`cell_edges`, doc/regions.md). The
// region's dim and tag ride in `<Information Name="meshio++:dim|tag">`.
// ---------------------------------------------------------------------------

std::vector<std::int64_t> xdmf_int_values(const NDArray& rArr) {
    std::vector<std::int64_t> out(rArr.Size());
    // One dtype switch for the array, and no converted copy beside `out`.
    detail::dispatch_dtype(rArr.Dtype(), [&]<class T>() {
        const T* src = rArr.As<T>();
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = detail::typed_view_int<T>(src[i]);
    });
    return out;
}

/// The `<Set>`'s region kind, or false for a SetType this reader skips.
bool xdmf_set_kind(const std::string& rSetType, RegionKind& rKind) {
    if (rSetType == "Node")
        rKind = RegionKind::Point;
    else if (rSetType == "Cell")
        rKind = RegionKind::Cell;
    else if (rSetType == "Face" || rSetType == "Edge")
        rKind = RegionKind::Side;
    else
        return false;
    return true;
}

void xdmf_set_dim_tag(const pugi::xml_node& rSet, int& rDim, std::int64_t& rTag) {
    for (pugi::xml_node info : rSet.children("Information")) {
        const std::string name = info.attribute("Name").value();
        if (name == "meshio++:dim")
            rDim = info.attribute("Value").as_int(-1);
        else if (name == "meshio++:tag")
            rTag = info.attribute("Value").as_llong(-1);
    }
}

/// An empty `<DataItem>` (`Dimensions="0"`) has no payload to read.
bool xdmf_empty_item(const pugi::xml_node& rDi) {
    return std::string(rDi.attribute("Dimensions").value()) == "0";
}

/// Read one `<Set>` into @p rRegions, merging a name already met (a Side
/// region written as a Face and an Edge set is one region).
void xdmf_read_set(const pugi::xml_node& rSet, const fs::path& rBaseDir,
                   std::vector<Region>& rRegions) {
    const std::string name = rSet.attribute("Name").value();
    const std::string set_type = rSet.attribute("SetType").value();
    RegionKind kind{};
    if (!xdmf_set_kind(set_type, kind)) {
        log::warn("xdmf: skipping set '{}' of SetType '{}'", name, set_type);
        return;
    }
    int dim = -1;
    std::int64_t tag = -1;
    xdmf_set_dim_tag(rSet, dim, tag);
    std::vector<pugi::xml_node> items;
    for (pugi::xml_node di : rSet.children("DataItem"))
        items.push_back(di);
    if (rSet.child("Attribute"))
        log::warn("xdmf: set '{}' carries attributes, which are not read", name);
    const std::size_t need = kind == RegionKind::Side ? 2 : 1;
    if (items.size() < need)
        throw ReadError("XDMF: set '" + name + "' of SetType '" + set_type + "' needs " +
                        std::to_string(need) + " DataItem(s)");
    std::vector<std::int64_t> flat;
    if (!xdmf_empty_item(items[0])) {
        const std::vector<std::int64_t> ids = xdmf_int_values(read_data_item(items[0], rBaseDir));
        if (kind != RegionKind::Side) {
            flat = ids;
        } else {
            const std::vector<std::int64_t> local =
                xdmf_int_values(read_data_item(items[1], rBaseDir));
            if (local.size() != ids.size())
                throw ReadError("XDMF: set '" + name + "' has " + std::to_string(ids.size()) +
                                " cells but " + std::to_string(local.size()) +
                                " local face/edge indices");
            for (std::size_t i = 0; i < ids.size(); ++i) {
                flat.push_back(ids[i]);
                flat.push_back(local[i]);
            }
        }
    }
    for (Region& r_region : rRegions)
        if (r_region.mName == name && r_region.mKind == kind) {
            std::vector<std::int64_t> merged = xdmf_int_values(r_region.mEntries);
            merged.insert(merged.end(), flat.begin(), flat.end());
            flat = std::move(merged);
            r_region = Region(name, kind, r_region.mDim, r_region.mTag, NDArray());
            break;
        }
    const std::size_t stride = kind == RegionKind::Side ? 2 : 1;
    NDArray entries = stride == 1 ? NDArray(DType::Int64, {flat.size()})
                                  : NDArray(DType::Int64, {flat.size() / 2, 2});
    std::copy(flat.begin(), flat.end(), entries.As<std::int64_t>());
    for (Region& r_region : rRegions)
        if (r_region.mName == name && r_region.mKind == kind) {
            r_region.mEntries = std::move(entries);
            return;
        }
    rRegions.emplace_back(name, kind, dim, tag, std::move(entries));
}

void xdmf_attach_regions(Mesh& rMesh, std::vector<Region>& rRegions) {
    for (Region& r_region : rRegions)
        rMesh.AddRegion(std::move(r_region));
}

void xdmf_read_information(const pugi::xml_node& rNode, Mesh& rMesh, const ReadOptions& rOpts) {
    pugi::xml_document info;
    if (!info.load_string(rNode.text().get()))
        throw ReadError("XDMF: malformed Information payload");
    for (pugi::xml_node entry : info.document_element().children()) {
        if (!entry.attribute("key") || !entry.attribute("dim"))
            throw ReadError("XDMF: Information entry needs key and dim");
        if (!rOpts.WantsArray(entry.attribute("key").value()))
            continue;
        std::int64_t tag = 0, dim = 0;
        detail::TextStream tag_stream(entry.text().get());
        detail::TextStream dim_stream(entry.attribute("dim").value());
        if (!(tag_stream >> tag) || !(dim_stream >> dim))
            throw ReadError("XDMF: invalid Information tag or dimension");
        std::string extra;
        if ((tag_stream >> extra) || (dim_stream >> extra))
            throw ReadError("XDMF: invalid Information tag or dimension");
        NDArray data(DType::Int64, {2});
        data.As<std::int64_t>()[0] = tag;
        data.As<std::int64_t>()[1] = dim;
        rMesh.AddFieldData(entry.attribute("key").value(), std::move(data));
    }
}
}  // namespace

Mesh read_xdmf(const std::string& rPath, const ReadOptions& rOpts) {
    pugi::xml_document doc;
    if (!doc.load_file(rPath.c_str()))
        throw ReadError("XDMF: could not parse " + rPath);
    XdmfDoc parsed = xdmf_resolve(doc);
    const bool want_data = rOpts.WantsAnyData();

    fs::path base_dir =
        fs::path(rPath).has_parent_path() ? fs::path(rPath).parent_path() : fs::path(".");

    Mesh mesh;
    std::vector<std::pair<std::string, NDArray>> point_data;  // preserve order
    std::vector<std::pair<std::string, NDArray>> cell_data_raw;

    if (!parsed.mSteps.empty()) {
        // Temporal collection: the geometry lives once in the mesh grid and the
        // requested step's <Grid> carries that step's attributes.
        xdmf_read_geometry(parsed.mMeshGrid, base_dir, mesh);
        std::vector<Region> regions;
        for (pugi::xml_node set : parsed.mMeshGrid.children("Set"))
            xdmf_read_set(set, base_dir, regions);
        xdmf_attach_regions(mesh, regions);
        const std::size_t k = rOpts.ResolveTimeStep(parsed.mSteps.size());
        for (pugi::xml_node c : parsed.mSteps[k].children()) {
            if (std::string(c.name()) != "Attribute")
                continue;  // <Time>, the xi:include placeholder, ...
            const std::string name = c.attribute("Name").value();
            const std::string center = c.attribute("Center").value();
            if (!want_data || !rOpts.WantsArray(name))
                continue;
            NDArray data = read_data_item(c.child("DataItem"), base_dir);
            if (center == "Node")
                point_data.emplace_back(name, std::move(data));
            else if (center == "Cell")
                cell_data_raw.emplace_back(name, std::move(data));
            else
                throw ReadError("XDMF: unknown attribute center " + center);
        }
        xdmf_attach_data(mesh, point_data, cell_data_raw);
        return mesh;
    }

    pugi::xml_node grid = parsed.mMeshGrid;
    std::vector<Region> regions;
    for (pugi::xml_node c : grid.children()) {
        std::string tag = c.name();
        if (tag == "Topology") {
            std::string ctype = c.attribute("Type") ? c.attribute("Type").value()
                                                    : c.attribute("TopologyType").value();
            pugi::xml_node di = c.child("DataItem");
            NDArray data = read_data_item(di, base_dir);
            if (ctype == "Mixed") {
                translate_mixed(data, mesh);
            } else {
                mesh.AddCellBlock(xdmf_to_meshio(ctype), std::move(data));
            }
        } else if (tag == "Geometry") {
            pugi::xml_node di = c.child("DataItem");
            mesh.AssignPoints(read_data_item(di, base_dir));
        } else if (tag == "Attribute") {
            std::string name = c.attribute("Name").value();
            std::string center = c.attribute("Center").value();
            // Name/Center are attributes, so an unwanted attribute is skipped
            // before its DataItem is touched -- for the HDF path that means the
            // dataset is never opened at all.
            if (!want_data || !rOpts.WantsArray(name))
                continue;
            pugi::xml_node di = c.child("DataItem");
            NDArray data = read_data_item(di, base_dir);
            if (center == "Node")
                point_data.emplace_back(name, std::move(data));
            else if (center == "Cell")
                cell_data_raw.emplace_back(name, std::move(data));
            else
                throw ReadError("XDMF: unknown attribute center " + center);
        } else if (tag == "Set") {
            xdmf_read_set(c, base_dir, regions);
        } else if (tag == "Information") {
            if (want_data)
                xdmf_read_information(c, mesh, rOpts);
        } else {
            throw ReadError("XDMF: unknown section " + tag);
        }
    }

    // Split raw cell data into per-block arrays (cell_data_from_raw).
    xdmf_attach_data(mesh, point_data, cell_data_raw);
    xdmf_attach_regions(mesh, regions);

    return mesh;
}

MeshMetadata read_xdmf_metadata(const std::string& rPath, const ReadOptions&) {
    pugi::xml_document doc;
    if (!doc.load_file(rPath.c_str()))
        throw ReadError("XDMF: could not parse " + rPath);
    XdmfDoc parsed = xdmf_resolve(doc);

    MeshMetadata meta;
    // XDMF is the best case for a summary: every <DataItem> declares its shape
    // in a `Dimensions` attribute, so counts are exact without reading any
    // payload -- and for the HDF path, without opening the .h5 file at all.
    if (!parsed.mSteps.empty()) {
        // A temporal collection's step count and time values are the whole point
        // of summarizing one, and both are XML attributes -- no payload at all.
        for (const pugi::xml_node& step : parsed.mSteps)
            meta.mTimeValues.push_back(xdmf_step_time(step));
        // Array names come from the first step: every step of a series carries
        // the same attributes, and reporting names that only exist at some other
        // step would be a worse answer than reporting the file's own shape.
        for (pugi::xml_node c : parsed.mSteps.front().children()) {
            if (std::string(c.name()) != "Attribute")
                continue;
            const std::string name = c.attribute("Name").value();
            const std::string center = c.attribute("Center").value();
            if (center == "Node")
                meta.mPointDataNames.push_back(name);
            else if (center == "Cell")
                meta.mCellDataNames.push_back(name);
            else
                throw ReadError("XDMF: unknown attribute center " + center);
        }
    }

    for (pugi::xml_node c : parsed.mMeshGrid.children()) {
        const std::string tag = c.name();
        if (tag == "Topology") {
            const std::string ctype = c.attribute("Type") ? c.attribute("Type").value()
                                                          : c.attribute("TopologyType").value();
            if (ctype == "Mixed")
                throw ReadError("XDMF: Mixed topology needs the full reader to be summarized");
            const std::vector<std::size_t> dims = parse_dims(
                xdmf_resolve_data_item(c.child("DataItem")).attribute("Dimensions").value());
            CellBlockInfo info;
            info.mType = xdmf_to_meshio(ctype);
            info.mNumCells = dims.empty() ? 0 : dims[0];
            info.mNodesPerCell = dims.size() >= 2 ? dims[1] : 0;
            meta.mCellBlocks.push_back(std::move(info));
        } else if (tag == "Geometry") {
            const std::vector<std::size_t> dims = parse_dims(
                xdmf_resolve_data_item(c.child("DataItem")).attribute("Dimensions").value());
            meta.mNumPoints = dims.empty() ? 0 : dims[0];
            meta.mPointDim = dims.size() >= 2 ? dims[1] : 3;
        } else if (!parsed.mSteps.empty() && tag != "Set") {
            // In a temporal file the mesh grid contributes geometry only; its
            // attributes (if it doubles as step 0) were already taken above, and
            // <Time>/xi:include are not sections this reader has to understand.
            continue;
        } else if (tag == "Attribute") {
            const std::string name = c.attribute("Name").value();
            const std::string center = c.attribute("Center").value();
            if (center == "Node")
                meta.mPointDataNames.push_back(name);
            else if (center == "Cell")
                meta.mCellDataNames.push_back(name);
            else
                throw ReadError("XDMF: unknown attribute center " + center);
        } else if (tag == "Set") {
            // Counted from the first DataItem's declared Dimensions; a Side
            // region written as a Face and an Edge set is one region.
            RegionKind kind{};
            if (!xdmf_set_kind(c.attribute("SetType").value(), kind))
                continue;
            RegionSummary rs;
            rs.mName = c.attribute("Name").value();
            rs.mKind = kind;
            xdmf_set_dim_tag(c, rs.mDim, rs.mTag);
            const std::vector<std::size_t> dims = parse_dims(
                xdmf_resolve_data_item(c.child("DataItem")).attribute("Dimensions").value());
            rs.mNumEntries = dims.empty() ? 0 : dims[0];
            bool merged = false;
            for (RegionSummary& r_prev : meta.mRegions)
                if (r_prev.mName == rs.mName && r_prev.mKind == rs.mKind) {
                    r_prev.mNumEntries += rs.mNumEntries;
                    merged = true;
                }
            if (!merged)
                meta.mRegions.push_back(std::move(rs));
        } else if (tag == "Information") {
            Mesh fields;
            xdmf_read_information(c, fields, {});
            meta.mFieldDataNames = fields.FieldDataNames();
        } else {
            throw ReadError("XDMF: unknown section " + tag);
        }
    }
    // Match the uniform API's sorted-name guarantee.
    std::sort(meta.mPointDataNames.begin(), meta.mPointDataNames.end());
    std::sort(meta.mCellDataNames.begin(), meta.mCellDataNames.end());
    std::sort(meta.mFieldDataNames.begin(), meta.mFieldDataNames.end());

    meta.mHasBBox = false;  // would require reading the Geometry payload
    return meta;
}

namespace {

// Append a <DataItem> under `parent` carrying `rArr`, storing the heavy data
// through `rStore`. The element side lives here (pugixml is build-only and no
// installed header may name it); the payload side is shared with the transient
// writer via xdmfcommon::DataItemStore.
void xdmf_add_data_item(pugi::xml_node parent, xdmfcommon::DataItemStore& rStore,
                        const NDArray& rArr) {
    auto [dtype_s, prec] = numpy_to_xdmf_dtype(rArr.Dtype());
    const std::string dims = dims_string(rArr);
    pugi::xml_node di = parent.append_child("DataItem");
    di.append_attribute("DataType") = dtype_s;
    di.append_attribute("Dimensions") = dims.c_str();
    di.append_attribute("Format") = rStore.DataFormat().c_str();
    di.append_attribute("Precision") = prec;
    di.text().set(rStore.Store(rArr).c_str());
}

/// One Int64 `<DataItem>`, or an empty one (`Dimensions="0"`) with no payload.
void xdmf_add_ids(pugi::xml_node parent, xdmfcommon::DataItemStore& rStore,
                  const std::vector<std::int64_t>& rIds) {
    if (rIds.empty()) {
        pugi::xml_node di = parent.append_child("DataItem");
        di.append_attribute("DataType") = "Int";
        di.append_attribute("Dimensions") = "0";
        di.append_attribute("Format") = "XML";
        di.append_attribute("Precision") = "8";
        return;
    }
    NDArray a(DType::Int64, {rIds.size()});
    std::copy(rIds.begin(), rIds.end(), a.As<std::int64_t>());
    xdmf_add_data_item(parent, rStore, a);
}

void xdmf_write_set(pugi::xml_node grid, xdmfcommon::DataItemStore& rStore, const Region& rRegion,
                    const char* pSetType, const std::vector<std::int64_t>& rIds,
                    const std::vector<std::int64_t>* pLocal) {
    pugi::xml_node set = grid.append_child("Set");
    set.append_attribute("Name") = rRegion.mName.c_str();
    set.append_attribute("SetType") = pSetType;
    if (rRegion.mDim != -1) {
        pugi::xml_node info = set.append_child("Information");
        info.append_attribute("Name") = "meshio++:dim";
        info.append_attribute("Value") = std::to_string(rRegion.mDim).c_str();
    }
    if (rRegion.mTag != -1) {
        pugi::xml_node info = set.append_child("Information");
        info.append_attribute("Name") = "meshio++:tag";
        info.append_attribute("Value") = std::to_string(rRegion.mTag).c_str();
    }
    xdmf_add_ids(set, rStore, rIds);
    if (pLocal)
        xdmf_add_ids(set, rStore, *pLocal);
}

}  // namespace

namespace xdmfdetail {

void xdmf_write_sets(pugi::xml_node grid, xdmfcommon::DataItemStore& rStore, const Mesh& rMesh) {
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r_region = rMesh.Region(i);
        const std::int64_t* e = r_region.Entries();
        const std::size_t n = r_region.NumEntries();
        if (r_region.mKind != RegionKind::Side) {
            xdmf_write_set(grid, rStore, r_region,
                           r_region.mKind == RegionKind::Point ? "Node" : "Cell",
                           std::vector<std::int64_t>(e, e + n), nullptr);
            continue;
        }
        // A facet of a 3-D cell is a face, of a 2-D one an edge: split by the
        // owning cell's dimension, one set each (the reader merges them).
        std::vector<std::int64_t> face_cells, face_local, edge_cells, edge_local;
        for (std::size_t k = 0; k < n; ++k) {
            const auto [b, row] = detail::global_to_block_row(bases, e[2 * k]);
            (void)row;
            int dim = 3;
            if (b != static_cast<std::size_t>(-1)) {
                const auto cb = rMesh.Cells(b);
                dim = cb.IsPolyhedron()
                          ? 3
                          : cell_type_dimension(cell_type_from_name(std::string(cb.Type())));
            }
            (dim == 2 ? edge_cells : face_cells).push_back(e[2 * k]);
            (dim == 2 ? edge_local : face_local).push_back(e[2 * k + 1]);
        }
        if (!face_cells.empty() || edge_cells.empty())
            xdmf_write_set(grid, rStore, r_region, "Face", face_cells, &face_local);
        if (!edge_cells.empty())
            xdmf_write_set(grid, rStore, r_region, "Edge", edge_cells, &edge_local);
    }
}

}  // namespace xdmfdetail

void write_xdmf(const std::string& rPath, const Mesh& rMesh, const std::string& rDataFormat,
                int gzip_level) {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    const bool hdf_ok = true;
#else
    const bool hdf_ok = false;
#endif
    if (rDataFormat != "XML" && rDataFormat != "Binary" && !(rDataFormat == "HDF" && hdf_ok))
        throw WriteError("XDMF C++ core cannot write data format " + rDataFormat);

    std::string base = rPath;
    std::size_t dot = base.find_last_of('.');
    if (dot != std::string::npos)
        base = base.substr(0, dot);

    xdmfcommon::DataItemStore store(rDataFormat, base, gzip_level);

    pugi::xml_document doc;
    pugi::xml_node xdmf = doc.append_child("Xdmf");
    xdmf.append_attribute("Version") = "3.0";
    pugi::xml_node domain = xdmf.append_child("Domain");
    pugi::xml_node grid = domain.append_child("Grid");
    grid.append_attribute("Name") = "Grid";

    // Geometry
    const NDArray& points = rMesh.Points();
    const std::size_t pdim = points.Shape().size() >= 2 ? points.Shape()[1] : 3;
    if (pdim > 3)
        throw WriteError("XDMF: can only write points up to dimension 3");
    const char* geo_type = (pdim == 1) ? "X" : (pdim == 2) ? "XY" : "XYZ";
    pugi::xml_node geo = grid.append_child("Geometry");
    geo.append_attribute("GeometryType") = geo_type;
    xdmf_add_data_item(geo, store, points);

    // Topology
    if (rMesh.NumCellBlocks() == 1) {
        const auto cb = rMesh.Cells(0);
        const NDArray& conn = cb.Conn();
        pugi::xml_node topo = grid.append_child("Topology");
        topo.append_attribute("TopologyType") = meshio_to_xdmf(cb.Type());
        topo.append_attribute("NumberOfElements") = std::to_string(cb.NumCells()).c_str();
        topo.append_attribute("NodesPerElement") = std::to_string(detail::cols(conn)).c_str();
        xdmf_add_data_item(topo, store, conn);
    } else if (rMesh.NumCellBlocks() > 1) {
        std::size_t total_cells = 0;
        NDArray cd = pack_mixed_topology(rMesh, total_cells);
        pugi::xml_node topo = grid.append_child("Topology");
        topo.append_attribute("TopologyType") = "Mixed";
        topo.append_attribute("NumberOfElements") = std::to_string(total_cells).c_str();
        xdmf_add_data_item(topo, store, cd);
    }

    // Point data (sorted key order for deterministic output)
    for (const auto& name : rMesh.PointDataNames()) {
        const NDArray& d = rMesh.PointData(name);
        pugi::xml_node att = grid.append_child("Attribute");
        att.append_attribute("Name") = name.c_str();
        att.append_attribute("AttributeType") = attribute_type(d.Shape()).c_str();
        att.append_attribute("Center") = "Node";
        xdmf_add_data_item(att, store, d);
    }

    // Cell data (concatenated across blocks: raw_from_cell_data)
    for (const auto& name : rMesh.CellDataNames()) {
        if (rMesh.CellDataNumBlocks(name) == 0)
            continue;
        NDArray raw = concat_cell_data(rMesh, name);
        pugi::xml_node att = grid.append_child("Attribute");
        att.append_attribute("Name") = name.c_str();
        att.append_attribute("AttributeType") = attribute_type(raw.Shape()).c_str();
        att.append_attribute("Center") = "Cell";
        xdmf_add_data_item(att, store, raw);
    }

    // Regions as <Set>s (see "<Set> <-> regions" above).
    xdmfdetail::xdmf_write_sets(grid, store, rMesh);

    if (!doc.save_file(rPath.c_str(), "  "))
        throw WriteError("XDMF: could not write " + rPath);
}

}  // namespace meshioplusplus
