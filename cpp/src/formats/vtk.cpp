#include "meshio/formats/vtk.hpp"

#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "meshio/detail/value_io.hpp"
#include "meshio/exceptions.hpp"
#include "meshio/vtk_common.hpp"

namespace meshio {

namespace {

using detail::cols;
using detail::is_float_dtype;
using detail::read_double;
using detail::read_int;

// VTK legacy dtype token for an NDArray dtype (numpy_to_vtk_dtype).
const char* vtk_dtype_str(DType dt) {
    switch (dt) {
        case DType::Float32: return "float";
        case DType::Float64: return "double";
        case DType::Int8: return "vtktypeint8";
        case DType::Int16: return "vtktypeint16";
        case DType::Int32: return "vtktypeint32";
        case DType::Int64: return "vtktypeint64";
        case DType::UInt8: return "vtktypeuint8";
        case DType::UInt16: return "vtktypeuint16";
        case DType::UInt32: return "vtktypeuint32";
        case DType::UInt64: return "vtktypeuint64";
    }
    return "double";
}

// Full-precision float so ASCII roundtrips exactly (test tol is 1e-15).
void write_double(std::ostream& os, double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    os << buf;
}

void write_field(std::ostream& os, const std::string& name, const NDArray& a,
                 std::size_t num_tuples, std::size_t num_components) {
    if (name.find(' ') != std::string::npos)
        throw WriteError("VTK doesn't support spaces in field names ('" + name + "').");
    os << name << ' ' << num_components << ' ' << num_tuples << ' '
       << vtk_dtype_str(a.dtype()) << '\n';
    const bool flt = is_float_dtype(a.dtype());
    const std::size_t n = a.size();
    for (std::size_t i = 0; i < n; ++i) {
        if (flt)
            write_double(os, read_double(a, i));
        else
            os << read_int(a, i);
        os << (i + 1 == n ? '\n' : ' ');
    }
    if (n == 0) os << '\n';
}

}  // namespace

void write_vtk_ascii_51(const std::string& path, const Mesh& mesh) {
    for (const auto& cb : mesh.cells) {
        if (cb.type.rfind("polyhedron", 0) == 0)
            throw WriteError("C++ VTK writer does not support polyhedron cells");
    }

    std::ofstream os(path, std::ios::binary);
    if (!os) throw WriteError("Could not open file for writing: " + path);

    const std::size_t num_points = mesh.num_points();
    const std::size_t dim = mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;

    std::size_t total_cells = 0;
    std::size_t total_idx = 0;
    for (const auto& cb : mesh.cells) {
        total_cells += cb.num_cells();
        total_idx += cb.data.size();
    }

    os << "# vtk DataFile Version 5.1\n";
    os << "written by meshio (C++ core)\n";
    os << "ASCII\n";
    os << "DATASET UNSTRUCTURED_GRID\n";

    // Points (always 3 components; pad 2D with zero z).
    os << "POINTS " << num_points << ' ' << vtk_dtype_str(mesh.points.dtype()) << '\n';
    for (std::size_t r = 0; r < num_points; ++r) {
        for (std::size_t c = 0; c < 3; ++c) {
            double v = (c < dim) ? read_double(mesh.points, r * dim + c) : 0.0;
            write_double(os, v);
            os << ((r + 1 == num_points && c == 2) ? '\n' : ' ');
        }
    }
    if (num_points == 0) os << '\n';

    // Cells: OFFSETS (num_cells + 1) and CONNECTIVITY (total_idx).
    os << "CELLS " << (total_cells + 1) << ' ' << total_idx << '\n';

    os << "OFFSETS vtktypeint64\n";
    std::int64_t running = 0;
    os << running << '\n';
    for (const auto& cb : mesh.cells) {
        const std::int64_t k = static_cast<std::int64_t>(cols(cb.data));
        for (std::size_t r = 0; r < cb.num_cells(); ++r) {
            running += k;
            os << running << '\n';
        }
    }

    os << "CONNECTIVITY vtktypeint64\n";
    for (const auto& cb : mesh.cells) {
        const std::size_t nc = cb.num_cells();
        const std::size_t k = cols(cb.data);
        std::vector<int> order = meshio_to_vtk_order(cb.type);
        for (std::size_t r = 0; r < nc; ++r) {
            for (std::size_t j = 0; j < k; ++j) {
                std::size_t col = order.empty() ? j : static_cast<std::size_t>(order[j]);
                os << read_int(cb.data, r * k + col) << '\n';
            }
        }
    }

    // Cell types.
    os << "CELL_TYPES " << total_cells << '\n';
    const auto& tmap = meshio_to_vtk_type();
    for (const auto& cb : mesh.cells) {
        auto it = tmap.find(cb.type);
        if (it == tmap.end()) throw WriteError("Unknown cell type for VTK: " + cb.type);
        for (std::size_t r = 0; r < cb.num_cells(); ++r) os << it->second << '\n';
    }

    // Point data.
    if (!mesh.point_data.empty()) {
        os << "POINT_DATA " << num_points << '\n';
        os << "FIELD FieldData " << mesh.point_data.size() << '\n';
        for (const auto& kv : mesh.point_data) {
            const NDArray& d = kv.second;
            write_field(os, kv.first, d, detail::rows(d), cols(d));
        }
    }

    // Cell data (concatenate per-block arrays for each name).
    if (!mesh.cell_data.empty()) {
        os << "CELL_DATA " << total_cells << '\n';
        os << "FIELD FieldData " << mesh.cell_data.size() << '\n';
        for (const auto& kv : mesh.cell_data) {
            const auto& blocks = kv.second;
            if (blocks.empty()) continue;
            std::size_t num_components = cols(blocks.front());
            if (num_components == 1 && blocks.front().shape().size() < 2)
                num_components = 1;
            // Header line.
            if (kv.first.find(' ') != std::string::npos)
                throw WriteError("VTK doesn't support spaces in field names ('" +
                                 kv.first + "').");
            os << kv.first << ' ' << num_components << ' ' << total_cells << ' '
               << vtk_dtype_str(blocks.front().dtype()) << '\n';
            const bool flt = is_float_dtype(blocks.front().dtype());
            std::size_t emitted = 0;
            for (const auto& blk : blocks) {
                const std::size_t n = blk.size();
                for (std::size_t i = 0; i < n; ++i) {
                    if (flt)
                        write_double(os, read_double(blk, i));
                    else
                        os << read_int(blk, i);
                    ++emitted;
                    os << ' ';
                }
            }
            os << '\n';
        }
    }
}

}  // namespace meshio
