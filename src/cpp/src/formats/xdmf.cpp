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
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes (private, not installed)
#include "xdmf_doc.hpp"

// Project includes
#include "meshioplusplus/formats/xdmf.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/xdmf_common.hpp"
#include "meshioplusplus/exceptions.hpp"

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

void store_token(NDArray& rA, std::size_t i, const std::string& rTok) {
    switch (rA.Dtype()) {
        case DType::Float32:
            rA.As<float>()[i] = std::strtof(rTok.c_str(), nullptr);
            break;
        case DType::Float64:
            rA.As<double>()[i] = std::strtod(rTok.c_str(), nullptr);
            break;
        case DType::Int8:
            rA.As<std::int8_t>()[i] =
                static_cast<std::int8_t>(std::strtoll(rTok.c_str(), nullptr, 10));
            break;
        case DType::Int16:
            rA.As<std::int16_t>()[i] =
                static_cast<std::int16_t>(std::strtoll(rTok.c_str(), nullptr, 10));
            break;
        case DType::Int32:
            rA.As<std::int32_t>()[i] =
                static_cast<std::int32_t>(std::strtoll(rTok.c_str(), nullptr, 10));
            break;
        case DType::Int64:
            rA.As<std::int64_t>()[i] = std::strtoll(rTok.c_str(), nullptr, 10);
            break;
        case DType::UInt8:
            rA.As<std::uint8_t>()[i] =
                static_cast<std::uint8_t>(std::strtoull(rTok.c_str(), nullptr, 10));
            break;
        case DType::UInt16:
            rA.As<std::uint16_t>()[i] =
                static_cast<std::uint16_t>(std::strtoull(rTok.c_str(), nullptr, 10));
            break;
        case DType::UInt32:
            rA.As<std::uint32_t>()[i] =
                static_cast<std::uint32_t>(std::strtoull(rTok.c_str(), nullptr, 10));
            break;
        case DType::UInt64:
            rA.As<std::uint64_t>()[i] = std::strtoull(rTok.c_str(), nullptr, 10);
            break;
    }
}

NDArray read_data_item(const pugi::xml_node& rDi, const fs::path& rBaseDir) {
    std::vector<std::size_t> dims = parse_dims(rDi.attribute("Dimensions").value());

    std::string data_type = "Float";
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
        NDArray a(dt, dims);
        std::istringstream iss(rDi.text().get());
        std::string tok;
        std::size_t i = 0;
        while (i < total && (iss >> tok))
            store_token(a, i++, tok);
        return a;
    }
    if (fmt == "Binary") {
        std::string rel = rDi.text().get();
        // trim whitespace
        std::size_t a0 = rel.find_first_not_of(" \t\r\n");
        std::size_t a1 = rel.find_last_not_of(" \t\r\n");
        std::string path = (a0 == std::string::npos) ? "" : rel.substr(a0, a1 - a0 + 1);
        std::ifstream bin(path, std::ios::binary);
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
    while (r < n) {
        int xt = static_cast<int>(detail::read_int(rFlat, r));
        types.push_back(xt);
        offsets.push_back(r);
        if (xt == 2) {  // polyline: next value is point count, must be 2
            if (detail::read_int(rFlat, r + 1) != 2)
                throw ReadError("XDMF: only 2-point lines supported");
            r += 1;
        }
        r += 1;
        r += static_cast<std::size_t>(xdmf_idx_num_nodes(xt));
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
        NDArray data(DType::Int64, {nrows, static_cast<std::size_t>(nn)});
        std::int64_t* dp = data.As<std::int64_t>();
        for (std::size_t b = 0; b < nrows; ++b) {
            std::size_t base = offsets[start + b] + (xt == 2 ? 2 : 1);
            for (int j = 0; j < nn; ++j)
                dp[b * nn + j] = detail::read_int(rFlat, base + j);
        }
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
        } else if (tag == "Information") {
            // field_data not handled by the C++ core
            throw ReadError("XDMF: Information section handled by Python fallback");
        } else {
            throw ReadError("XDMF: unknown section " + tag);
        }
    }

    // Split raw cell data into per-block arrays (cell_data_from_raw).
    xdmf_attach_data(mesh, point_data, cell_data_raw);

    return mesh;
}

MeshMetadata read_xdmf_metadata(const std::string& rPath, const ReadOptions&) {
    pugi::xml_document doc;
    if (!doc.load_file(rPath.c_str(), pugi::parse_minimal))
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
            const std::vector<std::size_t> dims =
                parse_dims(c.child("DataItem").attribute("Dimensions").value());
            CellBlockInfo info;
            info.mType = xdmf_to_meshio(ctype);
            info.mNumCells = dims.empty() ? 0 : dims[0];
            info.mNodesPerCell = dims.size() >= 2 ? dims[1] : 0;
            meta.mCellBlocks.push_back(std::move(info));
        } else if (tag == "Geometry") {
            const std::vector<std::size_t> dims =
                parse_dims(c.child("DataItem").attribute("Dimensions").value());
            meta.mNumPoints = dims.empty() ? 0 : dims[0];
            meta.mPointDim = dims.size() >= 2 ? dims[1] : 3;
        } else if (!parsed.mSteps.empty()) {
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
        } else if (tag == "Information") {
            throw ReadError("XDMF: Information section handled by Python fallback");
        } else {
            throw ReadError("XDMF: unknown section " + tag);
        }
    }
    // Match the uniform API's sorted-name guarantee.
    std::sort(meta.mPointDataNames.begin(), meta.mPointDataNames.end());
    std::sort(meta.mCellDataNames.begin(), meta.mCellDataNames.end());

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

}  // namespace

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

    if (!doc.save_file(rPath.c_str(), "  "))
        throw WriteError("XDMF: could not write " + rPath);
}

}  // namespace meshioplusplus
