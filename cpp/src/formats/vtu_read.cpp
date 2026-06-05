#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "pugixml.hpp"

#include "meshio/detail/value_io.hpp"
#include "meshio/detail/vtk_cells.hpp"
#include "meshio/detail/vtu_binary.hpp"
#include "meshio/exceptions.hpp"
#include "meshio/formats/vtu.hpp"
#include "meshio/types.hpp"
#include "meshio/vtk_common.hpp"

namespace meshio {

namespace {

DType dtype_from_vtu(const std::string& s) {
    if (s == "Float32") return DType::Float32;
    if (s == "Float64") return DType::Float64;
    if (s == "Int8") return DType::Int8;
    if (s == "Int16") return DType::Int16;
    if (s == "Int32") return DType::Int32;
    if (s == "Int64") return DType::Int64;
    if (s == "UInt8") return DType::UInt8;
    if (s == "UInt16") return DType::UInt16;
    if (s == "UInt32") return DType::UInt32;
    if (s == "UInt64") return DType::UInt64;
    throw ReadError("Illegal VTU data type '" + s + "'");
}

void store(NDArray& a, std::size_t i, double d, std::int64_t v, bool isflt) {
    switch (a.dtype()) {
        case DType::Float32: a.as<float>()[i] = static_cast<float>(d); break;
        case DType::Float64: a.as<double>()[i] = d; break;
        case DType::Int8: a.as<std::int8_t>()[i] = static_cast<std::int8_t>(v); break;
        case DType::Int16: a.as<std::int16_t>()[i] = static_cast<std::int16_t>(v); break;
        case DType::Int32: a.as<std::int32_t>()[i] = static_cast<std::int32_t>(v); break;
        case DType::Int64: a.as<std::int64_t>()[i] = v; break;
        case DType::UInt8: a.as<std::uint8_t>()[i] = static_cast<std::uint8_t>(v); break;
        case DType::UInt16: a.as<std::uint16_t>()[i] = static_cast<std::uint16_t>(v); break;
        case DType::UInt32: a.as<std::uint32_t>()[i] = static_cast<std::uint32_t>(v); break;
        case DType::UInt64: a.as<std::uint64_t>()[i] = static_cast<std::uint64_t>(v); break;
    }
    (void)isflt;
}

NDArray parse_ascii(const char* text, DType dt) {
    const bool isflt = detail::is_float_dtype(dt);
    std::vector<double> dv;
    std::vector<std::int64_t> iv;
    const char* p = text ? text : "";
    while (*p) {
        while (*p && std::isspace(static_cast<unsigned char>(*p))) ++p;
        if (!*p) break;
        char* endp = nullptr;
        if (isflt) {
            double x = std::strtod(p, &endp);
            if (endp == p) break;
            dv.push_back(x);
        } else {
            long long x = std::strtoll(p, &endp, 10);
            if (endp == p) break;
            iv.push_back(static_cast<std::int64_t>(x));
        }
        p = endp;
    }
    std::size_t n = isflt ? dv.size() : iv.size();
    NDArray a(dt, {n});
    for (std::size_t i = 0; i < n; ++i)
        store(a, i, isflt ? dv[i] : 0.0, isflt ? 0 : iv[i], isflt);
    return a;
}

std::string strip(const char* s) {
    std::string t = s ? s : "";
    std::size_t b = 0, e = t.size();
    while (b < e && std::isspace(static_cast<unsigned char>(t[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(t[e - 1]))) --e;
    return t.substr(b, e - b);
}

NDArray parse_binary(const std::string& text, DType dt, int compression,
                     std::size_t hsz) {
    std::vector<unsigned char> bytes;
    if (compression == 0)
        bytes = detail::vtu_decode_uncompressed(text.c_str(), text.size(), hsz);
    else
        bytes = detail::vtu_decode_zlib(text.c_str(), text.size(), hsz);
    std::size_t isz = dtype_size(dt);
    std::size_t n = isz ? bytes.size() / isz : 0;
    NDArray a(dt, {n});
    if (n) std::memcpy(a.data(), bytes.data(), n * isz);
    return a;
}

// compression: 0 = none, 1 = zlib. lzma/appended raise (handled by caller).
NDArray read_data_array(const pugi::xml_node& da, int compression, std::size_t hsz,
                        int& num_components) {
    std::string fmt = da.attribute("format").as_string("ascii");
    DType dt = dtype_from_vtu(da.attribute("type").as_string());
    num_components = da.attribute("NumberOfComponents").as_int(0);

    if (fmt == "ascii") return parse_ascii(da.text().get(), dt);
    if (fmt == "binary") return parse_binary(strip(da.text().get()), dt, compression, hsz);
    throw ReadError("VTU '" + fmt + "' data is not supported by the C++ reader");
}

std::vector<std::int64_t> to_int64(const NDArray& a) {
    std::vector<std::int64_t> v(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) v[i] = detail::read_int(a, i);
    return v;
}

}  // namespace

Mesh read_vtu(const std::string& path) {
    pugi::xml_document doc;
    pugi::xml_parse_result res = doc.load_file(path.c_str());
    if (!res)
        throw ReadError(std::string("VTU XML parse failed: ") + res.description());

    pugi::xml_node root = doc.child("VTKFile");
    if (!root) throw ReadError("Expected tag 'VTKFile'");
    if (std::string(root.attribute("type").as_string()) != "UnstructuredGrid")
        throw ReadError("Expected type UnstructuredGrid");

    int compression = 0;  // 0 none, 1 zlib
    std::string compressor = root.attribute("compressor").as_string("");
    if (compressor == "vtkZLibDataCompressor")
        compression = 1;
    else if (compressor == "vtkLZMADataCompressor")
        throw ReadError("lzma-compressed VTU not supported by the C++ reader");
    else if (!compressor.empty())
        throw ReadError("Unknown VTU compressor '" + compressor + "'");

    std::string header_type = root.attribute("header_type").as_string("UInt32");
    std::size_t hsz = (header_type == "UInt64") ? 8 : 4;

    pugi::xml_node grid = root.child("UnstructuredGrid");
    if (!grid) throw ReadError("No UnstructuredGrid found");

    // Appended data is not handled here -> let the Python reader take over.
    if (grid.parent().child("AppendedData") || root.child("AppendedData"))
        throw ReadError("appended VTU data not supported by the C++ reader");

    pugi::xml_node piece = grid.child("Piece");
    if (!piece) throw ReadError("No Piece found");
    // A single piece is supported; multiple pieces -> Python reader.
    if (piece.next_sibling("Piece"))
        throw ReadError("multi-piece VTU not supported by the C++ reader");

    std::size_t num_points =
        static_cast<std::size_t>(piece.attribute("NumberOfPoints").as_ullong());

    Mesh mesh;
    std::vector<std::int64_t> conn, offsets, types;
    std::map<std::string, NDArray> cell_data_raw;

    for (pugi::xml_node child : piece.children()) {
        std::string tag = child.name();
        if (tag == "Points") {
            pugi::xml_node da = child.child("DataArray");
            int nc = 0;
            NDArray pts = read_data_array(da, compression, hsz, nc);
            if (nc <= 0) nc = 3;
            pts.reshape({num_points, static_cast<std::size_t>(nc)});
            mesh.points = std::move(pts);
        } else if (tag == "Cells") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = read_data_array(da, compression, hsz, nc);
                if (name == "connectivity")
                    conn = to_int64(arr);
                else if (name == "offsets")
                    offsets = to_int64(arr);
                else if (name == "types")
                    types = to_int64(arr);
                else if (name == "faces" || name == "faceoffsets")
                    throw ReadError("polyhedron VTU not supported by the C++ reader");
            }
        } else if (tag == "PointData") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = read_data_array(da, compression, hsz, nc);
                if (nc > 1) arr.reshape({arr.size() / nc, static_cast<std::size_t>(nc)});
                mesh.point_data.emplace(name, std::move(arr));
            }
        } else if (tag == "CellData") {
            for (pugi::xml_node da : child.children("DataArray")) {
                int nc = 0;
                std::string name = da.attribute("Name").as_string();
                NDArray arr = read_data_array(da, compression, hsz, nc);
                if (nc > 1) arr.reshape({arr.size() / nc, static_cast<std::size_t>(nc)});
                cell_data_raw.emplace(name, std::move(arr));
            }
        }
    }

    detail::reconstruct_cells(conn, offsets, types, cell_data_raw, mesh.cells,
                              mesh.cell_data);
    return mesh;
}

}  // namespace meshio
