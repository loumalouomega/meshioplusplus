#include "meshio/formats/vtu.hpp"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "meshio/detail/value_io.hpp"
#include "meshio/exceptions.hpp"
#include "meshio/types.hpp"
#include "meshio/vtk_common.hpp"

namespace meshio {

namespace {

using detail::cols;
using detail::is_float_dtype;
using detail::read_double;
using detail::read_int;

// VTU "type" attribute string for an NDArray dtype.
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

void write_double(std::ostream& os, double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.11e", v);
    os << buf << '\n';
}

void write_int(std::ostream& os, std::int64_t v) {
    os << v << '\n';
}

// Append a whole NDArray's values (row-major) using meshio's ASCII formatting:
// "%.11e" for floats, "%d" for integers.
void write_ndarray_values(std::ostream& os, const NDArray& a) {
    const bool flt = is_float_dtype(a.dtype());
    const std::size_t n = a.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (flt)
            write_double(os, read_double(a, i));
        else
            write_int(os, read_int(a, i));
    }
}

}  // namespace

void write_vtu_ascii(const std::string& path, const Mesh& mesh) {
    for (const auto& cb : mesh.cells) {
        if (cb.type.rfind("polyhedron", 0) == 0)
            throw WriteError("C++ VTU writer does not support polyhedron cells");
    }

    std::ofstream os(path, std::ios::binary);
    if (!os) throw WriteError("Could not open file for writing: " + path);

    const std::size_t num_points = mesh.num_points();
    const std::size_t dim = mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;

    std::size_t total_cells = 0;
    for (const auto& cb : mesh.cells) total_cells += cb.num_cells();

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
          "byte_order=\"LittleEndian\">\n";
    os << "<!--This file was created by meshio (C++ core)-->\n";
    os << "<UnstructuredGrid>\n";
    os << "<Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\""
       << total_cells << "\">\n";

    // Points (always 3 components; pad 2D with a zero z).
    os << "<Points>\n";
    os << "<DataArray type=\"" << vtu_type_str(mesh.points.dtype())
       << "\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (std::size_t r = 0; r < num_points; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            double v = (c < dim) ? read_double(mesh.points, r * dim + c) : 0.0;
            write_double(os, v);
        }
    }
    os << "</DataArray>\n</Points>\n";

    // Cells: connectivity, offsets, types.
    if (!mesh.cells.empty()) {
        os << "<Cells>\n";

        // connectivity
        os << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        for (const auto& cb : mesh.cells) {
            const std::size_t nc = cb.num_cells();
            const std::size_t k = cols(cb.data);
            std::vector<int> order = meshio_to_vtk_order(cb.type);
            for (std::size_t r = 0; r < nc; ++r) {
                for (std::size_t j = 0; j < k; ++j) {
                    std::size_t col = order.empty() ? j : static_cast<std::size_t>(order[j]);
                    write_int(os, read_int(cb.data, r * k + col));
                }
            }
        }
        os << "</DataArray>\n";

        // offsets (cumulative end offsets)
        os << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        std::int64_t running = 0;
        for (const auto& cb : mesh.cells) {
            const std::size_t nc = cb.num_cells();
            const std::int64_t k = static_cast<std::int64_t>(cols(cb.data));
            for (std::size_t r = 0; r < nc; ++r) {
                running += k;
                write_int(os, running);
            }
        }
        os << "</DataArray>\n";

        // types
        os << "<DataArray type=\"Int64\" Name=\"types\" format=\"ascii\">\n";
        const auto& tmap = meshio_to_vtk_type();
        for (const auto& cb : mesh.cells) {
            auto it = tmap.find(cb.type);
            if (it == tmap.end())
                throw WriteError("Unknown cell type for VTU: " + cb.type);
            const std::int64_t vtk_id = it->second;
            for (std::size_t r = 0; r < cb.num_cells(); ++r) write_int(os, vtk_id);
        }
        os << "</DataArray>\n";

        os << "</Cells>\n";
    }

    // PointData
    if (!mesh.point_data.empty()) {
        os << "<PointData>\n";
        for (const auto& kv : mesh.point_data) {
            const NDArray& d = kv.second;
            os << "<DataArray type=\"" << vtu_type_str(d.dtype()) << "\" Name=\""
               << kv.first << "\"";
            if (d.shape().size() == 2) os << " NumberOfComponents=\"" << cols(d) << "\"";
            os << " format=\"ascii\">\n";
            write_ndarray_values(os, d);
            os << "</DataArray>\n";
        }
        os << "</PointData>\n";
    }

    // CellData (concatenate the per-block arrays for each name).
    if (!mesh.cell_data.empty()) {
        os << "<CellData>\n";
        for (const auto& kv : mesh.cell_data) {
            const auto& blocks = kv.second;
            if (blocks.empty()) continue;
            const NDArray& first = blocks.front();
            os << "<DataArray type=\"" << vtu_type_str(first.dtype()) << "\" Name=\""
               << kv.first << "\"";
            if (first.shape().size() == 2)
                os << " NumberOfComponents=\"" << cols(first) << "\"";
            os << " format=\"ascii\">\n";
            for (const auto& blk : blocks) write_ndarray_values(os, blk);
            os << "</DataArray>\n";
        }
        os << "</CellData>\n";
    }

    os << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
}

}  // namespace meshio
