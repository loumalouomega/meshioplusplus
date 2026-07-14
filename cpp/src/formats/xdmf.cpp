#include "meshioplusplus/formats/xdmf.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/xdmf_common.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "pugixml.hpp"

#ifdef MESHIOPLUSPLUS_HAS_HDF5
#include "meshioplusplus/detail/hdf5_util.hpp"
#endif

namespace fs = std::filesystem;

namespace meshioplusplus {

namespace {

// ---- type maps (shared with HMF via detail/xdmf_common.hpp) ----

using xdmfcommon::concat_cell_data;
using xdmfcommon::meshio_to_xdmf;
using xdmfcommon::split_raw_cell_data;
using xdmfcommon::xdmf_to_meshio;

int meshio_to_xdmf_index(const std::string& t) {
    static const std::unordered_map<std::string, int> m = {
        {"vertex", 0x1}, {"line", 0x2}, {"triangle", 0x4}, {"quad", 0x5},
        {"tetra", 0x6}, {"pyramid", 0x7}, {"wedge", 0x8}, {"hexahedron", 0x9},
        {"line3", 0x22}, {"quad9", 0x23}, {"triangle6", 0x24}, {"quad8", 0x25},
        {"tetra10", 0x26}, {"pyramid13", 0x27}, {"wedge15", 0x28}, {"wedge18", 0x29},
        {"hexahedron20", 0x30}, {"hexahedron24", 0x31}, {"hexahedron27", 0x32}};
    auto it = m.find(t);
    if (it == m.end()) throw WriteError("XDMF: cannot mix cell type " + t);
    return it->second;
}

std::string xdmf_idx_to_meshio(int idx) {
    static const std::unordered_map<int, std::string> m = {
        {0x1, "vertex"}, {0x2, "line"}, {0x4, "triangle"}, {0x5, "quad"},
        {0x6, "tetra"}, {0x7, "pyramid"}, {0x8, "wedge"}, {0x9, "hexahedron"},
        {0x22, "line3"}, {0x23, "quad9"}, {0x24, "triangle6"}, {0x25, "quad8"},
        {0x26, "tetra10"}, {0x27, "pyramid13"}, {0x28, "wedge15"}, {0x29, "wedge18"},
        {0x30, "hexahedron20"}, {0x31, "hexahedron24"}, {0x32, "hexahedron27"}};
    auto it = m.find(idx);
    if (it == m.end()) throw ReadError("XDMF: unknown mixed topology index");
    return it->second;
}

int xdmf_idx_num_nodes(int idx) {
    static const std::unordered_map<int, int> m = {
        {1, 1}, {2, 2}, {4, 3}, {5, 4}, {6, 4}, {7, 5}, {8, 6}, {9, 8}, {11, 6},
        {0x22, 3}, {0x23, 9}, {0x24, 6}, {0x25, 8}, {0x26, 10}, {0x27, 13},
        {0x28, 15}, {0x29, 18}, {0x30, 20}, {0x31, 24}, {0x32, 27}};
    auto it = m.find(idx);
    if (it == m.end()) throw ReadError("XDMF: unknown mixed topology index");
    return it->second;
}

std::pair<const char*, const char*> numpy_to_xdmf_dtype(DType dt) {
    switch (dt) {
        case DType::Int8: return {"Int", "1"};
        case DType::Int16: return {"Int", "2"};
        case DType::Int32: return {"Int", "4"};
        case DType::Int64: return {"Int", "8"};
        case DType::UInt8: return {"UInt", "1"};
        case DType::UInt16: return {"UInt", "2"};
        case DType::UInt32: return {"UInt", "4"};
        case DType::UInt64: return {"UInt", "8"};
        case DType::Float32: return {"Float", "4"};
        case DType::Float64: return {"Float", "8"};
    }
    return {"Float", "8"};
}

DType xdmf_to_dtype(const std::string& data_type, const std::string& precision) {
    int p = std::atoi(precision.c_str());
    if (data_type == "Int")
        return p == 1 ? DType::Int8 : p == 2 ? DType::Int16 : p == 4 ? DType::Int32 : DType::Int64;
    if (data_type == "UInt")
        return p == 1 ? DType::UInt8 : p == 2 ? DType::UInt16 : p == 4 ? DType::UInt32 : DType::UInt64;
    return p == 4 ? DType::Float32 : DType::Float64;
}

std::string attribute_type(const std::vector<std::size_t>& shape) {
    if (shape.size() == 1 || (shape.size() == 2 && shape[1] == 1)) return "Scalar";
    if (shape.size() == 2 && (shape[1] == 2 || shape[1] == 3)) return "Vector";
    if ((shape.size() == 2 && shape[1] == 9) ||
        (shape.size() == 3 && shape[1] == 3 && shape[2] == 3))
        return "Tensor";
    if (shape.size() == 2 && shape[1] == 6) return "Tensor6";
    return "Matrix";
}

std::vector<std::size_t> parse_dims(const std::string& s) {
    std::vector<std::size_t> dims;
    std::istringstream iss(s);
    std::int64_t v;
    while (iss >> v) dims.push_back(static_cast<std::size_t>(v));
    return dims;
}

void store_token(NDArray& a, std::size_t i, const std::string& tok) {
    switch (a.dtype()) {
        case DType::Float32: a.as<float>()[i] = std::strtof(tok.c_str(), nullptr); break;
        case DType::Float64: a.as<double>()[i] = std::strtod(tok.c_str(), nullptr); break;
        case DType::Int8: a.as<std::int8_t>()[i] = static_cast<std::int8_t>(std::strtoll(tok.c_str(), nullptr, 10)); break;
        case DType::Int16: a.as<std::int16_t>()[i] = static_cast<std::int16_t>(std::strtoll(tok.c_str(), nullptr, 10)); break;
        case DType::Int32: a.as<std::int32_t>()[i] = static_cast<std::int32_t>(std::strtoll(tok.c_str(), nullptr, 10)); break;
        case DType::Int64: a.as<std::int64_t>()[i] = std::strtoll(tok.c_str(), nullptr, 10); break;
        case DType::UInt8: a.as<std::uint8_t>()[i] = static_cast<std::uint8_t>(std::strtoull(tok.c_str(), nullptr, 10)); break;
        case DType::UInt16: a.as<std::uint16_t>()[i] = static_cast<std::uint16_t>(std::strtoull(tok.c_str(), nullptr, 10)); break;
        case DType::UInt32: a.as<std::uint32_t>()[i] = static_cast<std::uint32_t>(std::strtoull(tok.c_str(), nullptr, 10)); break;
        case DType::UInt64: a.as<std::uint64_t>()[i] = std::strtoull(tok.c_str(), nullptr, 10); break;
    }
}

NDArray read_data_item(const pugi::xml_node& di, const fs::path& base_dir) {
    std::vector<std::size_t> dims = parse_dims(di.attribute("Dimensions").value());

    std::string data_type = "Float";
    if (di.attribute("DataType")) data_type = di.attribute("DataType").value();
    else if (di.attribute("NumberType")) data_type = di.attribute("NumberType").value();
    std::string precision = di.attribute("Precision") ? di.attribute("Precision").value() : "4";
    std::string fmt = di.attribute("Format").value();
    DType dt = xdmf_to_dtype(data_type, precision);

    std::size_t total = dims.empty() ? 0 : std::accumulate(dims.begin(), dims.end(),
                                                            std::size_t{1}, std::multiplies<>());

    if (fmt == "XML") {
        NDArray a(dt, dims);
        std::istringstream iss(di.text().get());
        std::string tok;
        std::size_t i = 0;
        while (i < total && (iss >> tok)) store_token(a, i++, tok);
        return a;
    }
    if (fmt == "Binary") {
        std::string rel = di.text().get();
        // trim whitespace
        std::size_t a0 = rel.find_first_not_of(" \t\r\n");
        std::size_t a1 = rel.find_last_not_of(" \t\r\n");
        std::string path = (a0 == std::string::npos) ? "" : rel.substr(a0, a1 - a0 + 1);
        std::ifstream bin(path, std::ios::binary);
        if (!bin) {  // try relative to the xdmf file
            bin.open((base_dir / path).string(), std::ios::binary);
            if (!bin) throw ReadError("XDMF: could not open binary file " + path);
        }
        NDArray a(dt, dims);
        bin.read(reinterpret_cast<char*>(a.data()), static_cast<std::streamsize>(a.nbytes()));
        return a;
    }
    if (fmt != "HDF") throw ReadError("XDMF: unknown data format " + fmt);

#ifdef MESHIOPLUSPLUS_HAS_HDF5
    // "<file>.h5:/path/to/dataset", file path relative to the xdmf file.
    std::string info = di.text().get();
    std::size_t a0 = info.find_first_not_of(" \t\r\n");
    std::size_t a1 = info.find_last_not_of(" \t\r\n");
    info = (a0 == std::string::npos) ? "" : info.substr(a0, a1 - a0 + 1);
    std::size_t colon = info.find(':');
    if (colon == std::string::npos)
        throw ReadError("XDMF: malformed HDF reference '" + info + "'");
    std::string h5file = info.substr(0, colon);
    std::string h5path = info.substr(colon + 1);

    h5::SilenceErrors silence;
    fs::path full = base_dir / h5file;
    h5::Hid f = h5::open_file_read(full.string());
    NDArray a = h5::read_dataset(f, h5path);
    a.reshape(dims);  // stored shape is authoritative in the XML
    return a;
#else
    throw ReadError("XDMF: HDF data format handled by Python fallback");
#endif
}

// Mixed-topology translation (ported from common.translate_mixed_cells).
void translate_mixed(const NDArray& flat, std::vector<CellBlock>& out) {
    std::size_t n = flat.size();
    std::vector<int> types;
    std::vector<std::size_t> offsets;
    std::size_t r = 0;
    while (r < n) {
        int xt = static_cast<int>(detail::read_int(flat, r));
        types.push_back(xt);
        offsets.push_back(r);
        if (xt == 2) {  // polyline: next value is point count, must be 2
            if (detail::read_int(flat, r + 1) != 2)
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
        while (end < types.size() && types[end] == types[start]) ++end;
        int xt = types[start];
        int nn = xdmf_idx_num_nodes(xt);
        std::size_t nrows = end - start;
        NDArray data(DType::Int64, {nrows, static_cast<std::size_t>(nn)});
        std::int64_t* dp = data.as<std::int64_t>();
        for (std::size_t b = 0; b < nrows; ++b) {
            std::size_t base = offsets[start + b] + (xt == 2 ? 2 : 1);
            for (int j = 0; j < nn; ++j)
                dp[b * nn + j] = detail::read_int(flat, base + j);
        }
        out.emplace_back(xdmf_idx_to_meshio(xt), std::move(data));
        start = end;
    }
}

}  // namespace

Mesh read_xdmf(const std::string& path) {
    pugi::xml_document doc;
    if (!doc.load_file(path.c_str())) throw ReadError("XDMF: could not parse " + path);
    pugi::xml_node root = doc.child("Xdmf");
    if (!root) throw ReadError("XDMF: missing <Xdmf> root");
    std::string version = root.attribute("Version").value();
    if (!version.empty() && version[0] != '3')
        throw ReadError("XDMF: only version 3 handled by the C++ core");

    pugi::xml_node domain = root.child("Domain");
    pugi::xml_node grid = domain.child("Grid");
    if (!grid) throw ReadError("XDMF: missing <Grid>");

    fs::path base_dir = fs::path(path).has_parent_path() ? fs::path(path).parent_path()
                                                          : fs::path(".");

    Mesh mesh;
    std::vector<std::pair<std::string, NDArray>> point_data;  // preserve order
    std::vector<std::pair<std::string, NDArray>> cell_data_raw;

    for (pugi::xml_node c : grid.children()) {
        std::string tag = c.name();
        if (tag == "Topology") {
            std::string ctype = c.attribute("Type") ? c.attribute("Type").value()
                                                     : c.attribute("TopologyType").value();
            pugi::xml_node di = c.child("DataItem");
            NDArray data = read_data_item(di, base_dir);
            if (ctype == "Mixed") {
                translate_mixed(data, mesh.cells);
            } else {
                mesh.cells.emplace_back(xdmf_to_meshio(ctype), std::move(data));
            }
        } else if (tag == "Geometry") {
            pugi::xml_node di = c.child("DataItem");
            mesh.points = read_data_item(di, base_dir);
        } else if (tag == "Attribute") {
            std::string name = c.attribute("Name").value();
            std::string center = c.attribute("Center").value();
            pugi::xml_node di = c.child("DataItem");
            NDArray data = read_data_item(di, base_dir);
            if (center == "Node") point_data.emplace_back(name, std::move(data));
            else if (center == "Cell") cell_data_raw.emplace_back(name, std::move(data));
            else throw ReadError("XDMF: unknown attribute center " + center);
        } else if (tag == "Information") {
            // field_data not handled by the C++ core
            throw ReadError("XDMF: Information section handled by Python fallback");
        } else {
            throw ReadError("XDMF: unknown section " + tag);
        }
    }

    for (auto& kv : point_data) mesh.point_data.emplace(kv.first, std::move(kv.second));

    // Split raw cell data into per-block arrays (cell_data_from_raw).
    std::vector<std::size_t> sizes;
    for (const auto& cb : mesh.cells) sizes.push_back(cb.num_cells());
    for (auto& kv : cell_data_raw)
        mesh.cell_data.emplace(kv.first, split_raw_cell_data(kv.second, sizes));

    return mesh;
}

namespace {

struct XmlWriter {
    const std::string& data_format;
    std::string base;     // path without extension (for .bin / .h5 files)
    int counter = 0;
    int gzip_level = -1;  // HDF only
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    meshioplusplus::h5::Hid h5_file;   // lazily created sibling <base>.h5
    std::string h5_basename;
#endif

    // Append a <DataItem> under `parent` carrying `arr` in the chosen format.
    void add_data_item(pugi::xml_node parent, const NDArray& arr) {
        auto [dtype_s, prec] = numpy_to_xdmf_dtype(arr.dtype());
        std::string dims;
        for (std::size_t i = 0; i < arr.shape().size(); ++i) {
            if (i) dims += " ";
            dims += std::to_string(arr.shape()[i]);
        }
        pugi::xml_node di = parent.append_child("DataItem");
        di.append_attribute("DataType") = dtype_s;
        di.append_attribute("Dimensions") = dims.c_str();
        di.append_attribute("Format") = data_format.c_str();
        di.append_attribute("Precision") = prec;

        std::size_t rows = arr.shape().empty() ? 0 : arr.shape()[0];
        std::size_t cols = rows ? arr.size() / rows : 0;

        if (data_format == "Binary") {
            std::string fn = base + std::to_string(counter++) + ".bin";
            std::ofstream bf(fn, std::ios::binary);
            bf.write(reinterpret_cast<const char*>(arr.data()),
                     static_cast<std::streamsize>(arr.nbytes()));
            di.text().set(fn.c_str());
            return;
        }
        if (data_format == "HDF") {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
            if (!h5_file.valid()) {
                std::string h5_path = base + ".h5";
                h5_file = meshioplusplus::h5::create_file(h5_path);
                std::size_t slash = h5_path.find_last_of("/\\");
                h5_basename =
                    slash == std::string::npos ? h5_path : h5_path.substr(slash + 1);
            }
            std::string name = "data" + std::to_string(counter++);
            meshioplusplus::h5::write_dataset(h5_file, name, arr, gzip_level);
            di.text().set((h5_basename + ":/" + name).c_str());
            return;
#else
            throw WriteError("XDMF: HDF data format requires an HDF5-enabled build");
#endif
        }
        // XML inline
        std::string text = "\n";
        char buf[40];
        bool is_float = detail::is_float_dtype(arr.dtype());
        bool f32 = arr.dtype() == DType::Float32;
        for (std::size_t r = 0; r < rows; ++r) {
            for (std::size_t cc = 0; cc < cols; ++cc) {
                std::size_t i = r * cols + cc;
                if (is_float) {
                    std::snprintf(buf, sizeof(buf), f32 ? "%.7e" : "%.16e",
                                  detail::read_double(arr, i));
                } else {
                    std::snprintf(buf, sizeof(buf), "%lld",
                                  static_cast<long long>(detail::read_int(arr, i)));
                }
                if (cc) text += " ";
                text += buf;
            }
            text += "\n";
        }
        di.text().set(text.c_str());
    }
};

}  // namespace

void write_xdmf(const std::string& path, const Mesh& mesh,
                const std::string& data_format, int gzip_level) {
#ifdef MESHIOPLUSPLUS_HAS_HDF5
    const bool hdf_ok = true;
#else
    const bool hdf_ok = false;
#endif
    if (data_format != "XML" && data_format != "Binary" &&
        !(data_format == "HDF" && hdf_ok))
        throw WriteError("XDMF C++ core cannot write data format " + data_format);

    std::string base = path;
    std::size_t dot = base.find_last_of('.');
    if (dot != std::string::npos) base = base.substr(0, dot);

    XmlWriter w{data_format, base, 0, gzip_level};

    pugi::xml_document doc;
    pugi::xml_node xdmf = doc.append_child("Xdmf");
    xdmf.append_attribute("Version") = "3.0";
    pugi::xml_node domain = xdmf.append_child("Domain");
    pugi::xml_node grid = domain.append_child("Grid");
    grid.append_attribute("Name") = "Grid";

    // Geometry
    const std::size_t pdim = mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 3;
    if (pdim > 3) throw WriteError("XDMF: can only write points up to dimension 3");
    const char* geo_type = (pdim == 1) ? "X" : (pdim == 2) ? "XY" : "XYZ";
    pugi::xml_node geo = grid.append_child("Geometry");
    geo.append_attribute("GeometryType") = geo_type;
    w.add_data_item(geo, mesh.points);

    // Topology
    if (mesh.cells.size() == 1) {
        const CellBlock& cb = mesh.cells[0];
        pugi::xml_node topo = grid.append_child("Topology");
        topo.append_attribute("TopologyType") = meshio_to_xdmf(cb.type);
        topo.append_attribute("NumberOfElements") = std::to_string(cb.num_cells()).c_str();
        topo.append_attribute("NodesPerElement") =
            std::to_string(detail::cols(cb.data)).c_str();
        w.add_data_item(topo, cb.data);
    } else if (mesh.cells.size() > 1) {
        std::size_t total_cells = 0, total_len = 0;
        for (const auto& cb : mesh.cells) {
            std::size_t nc = cb.num_cells();
            std::size_t npc = detail::cols(cb.data);
            std::size_t prefix = (cb.type == "vertex" || cb.type == "line") ? 2 : 1;
            total_cells += nc;
            total_len += nc * (prefix + npc);
        }
        NDArray cd(DType::Int64, {total_len});
        std::int64_t* cp = cd.as<std::int64_t>();
        std::size_t pos = 0;
        for (const auto& cb : mesh.cells) {
            std::size_t nc = cb.num_cells();
            std::size_t npc = detail::cols(cb.data);
            int idx = meshio_to_xdmf_index(cb.type);
            std::size_t prefix = (cb.type == "vertex" || cb.type == "line") ? 2 : 1;
            for (std::size_t r = 0; r < nc; ++r) {
                for (std::size_t pq = 0; pq < prefix; ++pq) cp[pos++] = idx;
                for (std::size_t j = 0; j < npc; ++j)
                    cp[pos++] = detail::read_int(cb.data, r * npc + j);
            }
        }
        pugi::xml_node topo = grid.append_child("Topology");
        topo.append_attribute("TopologyType") = "Mixed";
        topo.append_attribute("NumberOfElements") = std::to_string(total_cells).c_str();
        w.add_data_item(topo, cd);
    }

    // Point data
    for (const auto& kv : mesh.point_data) {
        pugi::xml_node att = grid.append_child("Attribute");
        att.append_attribute("Name") = kv.first.c_str();
        att.append_attribute("AttributeType") = attribute_type(kv.second.shape()).c_str();
        att.append_attribute("Center") = "Node";
        w.add_data_item(att, kv.second);
    }

    // Cell data (concatenated across blocks: raw_from_cell_data)
    for (const auto& kv : mesh.cell_data) {
        if (kv.second.empty()) continue;
        NDArray raw = concat_cell_data(kv.second);
        pugi::xml_node att = grid.append_child("Attribute");
        att.append_attribute("Name") = kv.first.c_str();
        att.append_attribute("AttributeType") = attribute_type(raw.shape()).c_str();
        att.append_attribute("Center") = "Cell";
        w.add_data_item(att, raw);
    }

    if (!doc.save_file(path.c_str(), "  "))
        throw WriteError("XDMF: could not write " + path);
}

}  // namespace meshioplusplus
