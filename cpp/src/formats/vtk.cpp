#include "meshioplusplus/formats/vtk.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/vtk_common.hpp"

namespace meshioplusplus {

namespace {

using detail::cols;
using detail::is_float_dtype;
using detail::read_double;
using detail::read_int;

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

void ascii_double(std::ostream& os, double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.17g", v);
    os << buf;
}

// Write isz bytes of *p in reversed (big-endian) order.
void put_be(std::ostream& os, const unsigned char* p, std::size_t isz) {
    for (std::size_t i = isz; i-- > 0;) os.put(static_cast<char>(p[i]));
}
void put_be_i64(std::ostream& os, std::int64_t v) {
    put_be(os, reinterpret_cast<const unsigned char*>(&v), 8);
}
void put_be_i32(std::ostream& os, std::int32_t v) {
    put_be(os, reinterpret_cast<const unsigned char*>(&v), 4);
}

// Emit one NDArray element (index i) as big-endian raw bytes.
void put_be_elem(std::ostream& os, const NDArray& a, std::size_t i) {
    std::size_t isz = dtype_size(a.dtype());
    put_be(os, reinterpret_cast<const unsigned char*>(a.data()) + i * isz, isz);
}

// Byte-swap a whole array into a big-endian buffer (elements independent ->
// parallel), for a single os.write instead of per-element stream calls.
std::vector<unsigned char> be_buffer(const NDArray& a) {
    const std::size_t isz = dtype_size(a.dtype());
    const std::size_t n = a.size();
    const auto* src = reinterpret_cast<const unsigned char*>(a.data());
    std::vector<unsigned char> buf(n * isz);
    parallel_for(n, [&](std::size_t i) {
        for (std::size_t b = 0; b < isz; ++b)
            buf[i * isz + b] = src[i * isz + (isz - 1 - b)];
    });
    return buf;
}

void write_field_block(std::ostream& os, const std::string& name, DType dt,
                       std::size_t num_components, std::size_t num_tuples, bool binary,
                       const std::vector<const NDArray*>& blocks) {
    if (name.find(' ') != std::string::npos)
        throw WriteError("VTK doesn't support spaces in field names ('" + name + "').");
    os << name << ' ' << num_components << ' ' << num_tuples << ' ' << vtk_dtype_str(dt)
       << '\n';
    const bool flt = is_float_dtype(dt);
    for (const NDArray* blk : blocks) {
        if (binary) {
            std::vector<unsigned char> buf = be_buffer(*blk);
            os.write(reinterpret_cast<const char*>(buf.data()),
                     static_cast<std::streamsize>(buf.size()));
        } else {
            const std::size_t n = blk->size();
            for (std::size_t i = 0; i < n; ++i) {
                if (flt)
                    ascii_double(os, read_double(*blk, i));
                else
                    os << read_int(*blk, i);
                os << ' ';
            }
        }
    }
    os << '\n';
}

}  // namespace

void write_vtk(const std::string& path, const Mesh& mesh, bool binary, bool v51) {
    for (const auto& cb : mesh.cells)
        if (cb.type.rfind("polyhedron", 0) == 0)
            throw WriteError("C++ VTK writer does not support polyhedron cells");

    std::ofstream os(path, std::ios::binary);
    if (!os) throw WriteError("Could not open file for writing: " + path);

    const std::size_t num_points = mesh.num_points();
    const std::size_t dim = mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;
    const std::size_t pt_isz = dtype_size(mesh.points.dtype());

    std::size_t total_cells = 0, total_idx = 0;
    for (const auto& cb : mesh.cells) {
        total_cells += cb.num_cells();
        total_idx += cb.data.size();
    }

    os << (v51 ? "# vtk DataFile Version 5.1\n" : "# vtk DataFile Version 4.2\n");
    os << "written by meshio++ (C++ core)\n";
    os << (binary ? "BINARY\n" : "ASCII\n");
    os << "DATASET UNSTRUCTURED_GRID\n";

    // Points (3 components; pad 2D with zero z).
    os << "POINTS " << num_points << ' ' << vtk_dtype_str(mesh.points.dtype()) << '\n';
    if (binary) {
        // Pre-sized padded buffer, parallel byte-swap, then one write.
        const auto* src = reinterpret_cast<const unsigned char*>(mesh.points.data());
        std::vector<unsigned char> buf(num_points * 3 * pt_isz, 0);
        parallel_for(num_points, [&](std::size_t r) {
            for (std::size_t c = 0; c < dim && c < 3; ++c)
                for (std::size_t b = 0; b < pt_isz; ++b)
                    buf[(r * 3 + c) * pt_isz + b] =
                        src[(r * dim + c) * pt_isz + (pt_isz - 1 - b)];
        });
        os.write(reinterpret_cast<const char*>(buf.data()),
                 static_cast<std::streamsize>(buf.size()));
        os << '\n';
    } else {
        for (std::size_t r = 0; r < num_points; ++r)
            for (std::size_t c = 0; c < 3; ++c) {
                ascii_double(os, (c < dim) ? read_double(mesh.points, r * dim + c) : 0.0);
                os << ((r + 1 == num_points && c == 2) ? '\n' : ' ');
            }
        if (num_points == 0) os << '\n';
    }

    if (v51) {
        // Version 5.1: OFFSETS (num_cells + 1) and CONNECTIVITY (total_idx).
        os << "CELLS " << (total_cells + 1) << ' ' << total_idx << '\n';
        os << "OFFSETS vtktypeint64\n";
        std::int64_t running = 0;
        if (binary)
            put_be_i64(os, running);
        else
            os << running << '\n';
        for (const auto& cb : mesh.cells) {
            const std::int64_t k = static_cast<std::int64_t>(cols(cb.data));
            for (std::size_t r = 0; r < cb.num_cells(); ++r) {
                running += k;
                if (binary)
                    put_be_i64(os, running);
                else
                    os << running << '\n';
            }
        }
        if (binary) os << '\n';

        os << "CONNECTIVITY vtktypeint64\n";
        for (const auto& cb : mesh.cells) {
            const std::size_t nc = cb.num_cells();
            const std::size_t k = cols(cb.data);
            std::vector<int> order = meshio_to_vtk_order(cb.type);
            for (std::size_t r = 0; r < nc; ++r)
                for (std::size_t j = 0; j < k; ++j) {
                    std::size_t col = order.empty() ? j : static_cast<std::size_t>(order[j]);
                    std::int64_t v = read_int(cb.data, r * k + col);
                    if (binary)
                        put_be_i64(os, v);
                    else
                        os << v << '\n';
                }
        }
        if (binary) os << '\n';
    } else {
        // Version 4.2: interleaved [count, nodes...] per cell, as int32.
        os << "CELLS " << total_cells << ' ' << (total_idx + total_cells) << '\n';
        for (const auto& cb : mesh.cells) {
            const std::size_t nc = cb.num_cells();
            const std::int32_t k = static_cast<std::int32_t>(cols(cb.data));
            std::vector<int> order = meshio_to_vtk_order(cb.type);
            for (std::size_t r = 0; r < nc; ++r) {
                if (binary)
                    put_be_i32(os, k);
                else
                    os << k << '\n';
                for (int j = 0; j < k; ++j) {
                    std::size_t col = order.empty() ? j : static_cast<std::size_t>(order[j]);
                    std::int64_t v = read_int(cb.data, r * static_cast<std::size_t>(k) + col);
                    if (binary)
                        put_be_i32(os, static_cast<std::int32_t>(v));
                    else
                        os << v << '\n';
                }
            }
        }
        if (binary) os << '\n';
    }

    // Cell types.
    os << "CELL_TYPES " << total_cells << '\n';
    const auto& tmap = meshio_to_vtk_type();
    for (const auto& cb : mesh.cells) {
        auto it = tmap.find(cb.type);
        if (it == tmap.end()) throw WriteError("Unknown cell type for VTK: " + cb.type);
        for (std::size_t r = 0; r < cb.num_cells(); ++r) {
            if (binary)
                put_be_i32(os, it->second);
            else
                os << it->second << '\n';
        }
    }
    if (binary) os << '\n';

    // Point data.
    if (!mesh.point_data.empty()) {
        os << "POINT_DATA " << num_points << '\n';
        os << "FIELD FieldData " << mesh.point_data.size() << '\n';
        for (const auto& kv : mesh.point_data) {
            const NDArray& d = kv.second;
            std::size_t ncomp = cols(d);
            write_field_block(os, kv.first, d.dtype(), ncomp,
                              d.shape().empty() ? 0 : d.shape()[0], binary, {&d});
        }
    }

    // Cell data (concatenate per-block arrays for each name).
    if (!mesh.cell_data.empty()) {
        os << "CELL_DATA " << total_cells << '\n';
        os << "FIELD FieldData " << mesh.cell_data.size() << '\n';
        for (const auto& kv : mesh.cell_data) {
            const auto& blocks = kv.second;
            if (blocks.empty()) continue;
            std::vector<const NDArray*> ptrs;
            for (const auto& b : blocks) ptrs.push_back(&b);
            write_field_block(os, kv.first, blocks.front().dtype(),
                              cols(blocks.front()), total_cells, binary, ptrs);
        }
    }
}

}  // namespace meshioplusplus
