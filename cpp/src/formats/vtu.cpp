#include "meshio/formats/vtu.hpp"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "meshio/detail/value_io.hpp"
#include "meshio/detail/vtu_binary.hpp"
#include "meshio/exceptions.hpp"
#include "meshio/types.hpp"
#include "meshio/vtk_common.hpp"

namespace meshio {

namespace {

using detail::cols;
using detail::is_float_dtype;
using detail::read_double;
using detail::read_int;

const char* vtu_type_str(DType dt) {
    switch (dt) {
        case DType::Float32: return "Float32";
        case DType::Float64: return "Float64";
        case DType::Int8: return "Int8";
        case DType::Int16: return "Int16";
        case DType::Int32: return "Int32";
        case DType::Int64: return "Int64";
        case DType::UInt8: return "UInt8";
        case DType::UInt16: return "UInt16";
        case DType::UInt32: return "UInt32";
        case DType::UInt64: return "UInt64";
    }
    return "Float64";
}

void ascii_double(std::ostream& os, double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.11e", v);
    os << buf << '\n';
}

void ascii_ndarray(std::ostream& os, const NDArray& a) {
    const bool flt = is_float_dtype(a.dtype());
    const std::size_t n = a.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (flt)
            ascii_double(os, read_double(a, i));
        else
            os << read_int(a, i) << '\n';
    }
}

void append_elem(std::vector<unsigned char>& buf, const NDArray& a, std::size_t idx) {
    const std::size_t isz = dtype_size(a.dtype());
    const unsigned char* p = reinterpret_cast<const unsigned char*>(a.data()) + idx * isz;
    buf.insert(buf.end(), p, p + isz);
}

}  // namespace

void write_vtu(const std::string& path, const Mesh& mesh, bool binary, bool zlib) {
    for (const auto& cb : mesh.cells) {
        if (cb.type.rfind("polyhedron", 0) == 0)
            throw WriteError("C++ VTU writer does not support polyhedron cells");
    }

    std::ofstream os(path, std::ios::binary);
    if (!os) throw WriteError("Could not open file for writing: " + path);

    const std::size_t num_points = mesh.num_points();
    const std::size_t dim = mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;
    const std::size_t pt_isz = dtype_size(mesh.points.dtype());

    std::size_t total_cells = 0;
    for (const auto& cb : mesh.cells) total_cells += cb.num_cells();

    const char* fmt = binary ? "binary" : "ascii";

    auto da_header = [&](const char* type, const std::string& name, int ncomp) {
        os << "<DataArray type=\"" << type << "\" Name=\"" << name << "\"";
        if (ncomp > 0) os << " NumberOfComponents=\"" << ncomp << "\"";
        os << " format=\"" << fmt << "\">\n";
    };
    auto emit_bin = [&](const unsigned char* d, std::size_t n) {
        os << detail::vtu_encode_binary(d, n, zlib) << "\n";
    };

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
          "byte_order=\"LittleEndian\"";
    if (binary && zlib) os << " compressor=\"vtkZLibDataCompressor\"";
    os << ">\n";
    os << "<!--This file was created by meshio (C++ core)-->\n";
    os << "<UnstructuredGrid>\n";
    os << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\""
       << total_cells << "\">\n";

    // Points (3 components; pad 2D with zero z).
    os << "<Points>\n";
    da_header(vtu_type_str(mesh.points.dtype()), "Points", 3);
    if (binary) {
        std::vector<unsigned char> buf;
        buf.reserve(num_points * 3 * pt_isz);
        const std::vector<unsigned char> zero(pt_isz, 0);
        for (std::size_t r = 0; r < num_points; ++r)
            for (std::size_t c = 0; c < 3; ++c) {
                if (c < dim)
                    append_elem(buf, mesh.points, r * dim + c);
                else
                    buf.insert(buf.end(), zero.begin(), zero.end());
            }
        emit_bin(buf.data(), buf.size());
    } else {
        for (std::size_t r = 0; r < num_points; ++r)
            for (std::size_t c = 0; c < 3; ++c)
                ascii_double(os, (c < dim) ? read_double(mesh.points, r * dim + c) : 0.0);
    }
    os << "</DataArray>\n</Points>\n";

    if (!mesh.cells.empty()) {
        // Build connectivity / offsets / types (Int64).
        std::vector<std::int64_t> connectivity, offsets, types;
        const auto& tmap = meshio_to_vtk_type();
        std::int64_t running = 0;
        for (const auto& cb : mesh.cells) {
            const std::size_t nc = cb.num_cells();
            const std::size_t k = cols(cb.data);
            std::vector<int> order = meshio_to_vtk_order(cb.type);
            for (std::size_t r = 0; r < nc; ++r)
                for (std::size_t j = 0; j < k; ++j) {
                    std::size_t col = order.empty() ? j : static_cast<std::size_t>(order[j]);
                    connectivity.push_back(read_int(cb.data, r * k + col));
                }
            for (std::size_t r = 0; r < nc; ++r) {
                running += static_cast<std::int64_t>(k);
                offsets.push_back(running);
            }
            auto it = tmap.find(cb.type);
            if (it == tmap.end())
                throw WriteError("Unknown cell type for VTU: " + cb.type);
            for (std::size_t r = 0; r < nc; ++r) types.push_back(it->second);
        }

        auto emit_i64 = [&](const char* name, const std::vector<std::int64_t>& v) {
            da_header("Int64", name, 0);
            if (binary) {
                emit_bin(reinterpret_cast<const unsigned char*>(v.data()),
                         v.size() * sizeof(std::int64_t));
            } else {
                for (std::int64_t x : v) os << x << '\n';
            }
            os << "</DataArray>\n";
        };

        os << "<Cells>\n";
        emit_i64("connectivity", connectivity);
        emit_i64("offsets", offsets);
        emit_i64("types", types);
        os << "</Cells>\n";
    }

    if (!mesh.point_data.empty()) {
        os << "<PointData>\n";
        for (const auto& kv : mesh.point_data) {
            const NDArray& d = kv.second;
            int ncomp = (d.shape().size() == 2) ? static_cast<int>(cols(d)) : 0;
            da_header(vtu_type_str(d.dtype()), kv.first, ncomp);
            if (binary)
                emit_bin(reinterpret_cast<const unsigned char*>(d.data()), d.nbytes());
            else
                ascii_ndarray(os, d);
            os << "</DataArray>\n";
        }
        os << "</PointData>\n";
    }

    if (!mesh.cell_data.empty()) {
        os << "<CellData>\n";
        for (const auto& kv : mesh.cell_data) {
            const auto& blocks = kv.second;
            if (blocks.empty()) continue;
            const NDArray& first = blocks.front();
            int ncomp = (first.shape().size() == 2) ? static_cast<int>(cols(first)) : 0;
            da_header(vtu_type_str(first.dtype()), kv.first, ncomp);
            if (binary) {
                std::vector<unsigned char> buf;
                for (const auto& blk : blocks) {
                    const unsigned char* p =
                        reinterpret_cast<const unsigned char*>(blk.data());
                    buf.insert(buf.end(), p, p + blk.nbytes());
                }
                emit_bin(buf.data(), buf.size());
            } else {
                for (const auto& blk : blocks) ascii_ndarray(os, blk);
            }
            os << "</DataArray>\n";
        }
        os << "</CellData>\n";
    }

    os << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

}  // namespace meshio
