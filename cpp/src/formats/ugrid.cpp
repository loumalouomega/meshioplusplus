#include "meshioplusplus/formats/ugrid.hpp"

#include <cctype>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <iterator>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

// File flavour decoded from the penultimate filename suffix.
struct UgridType {
    bool ascii = true;
    bool fortran = false;   // Fortran record-length markers around each record
    bool big_endian = false;
    int float_size = 4;     // 4 or 8
    int int_size = 4;       // 4 or 8
};

UgridType resolve_type(const std::string& path) {
    // suffix table mirrors meshioplusplus.ugrid.file_types
    //   key -> {fortran, big_endian, float_size, int_size}
    struct Spec { bool fortran; bool big; int fs; int is; };
    static const std::map<std::string, Spec> table = {
        {"b8l", {false, true, 8, 8}}, {"b8", {false, true, 8, 4}},
        {"b4", {false, true, 4, 4}},  {"lb8l", {false, false, 8, 8}},
        {"lb8", {false, false, 8, 4}},{"lb4", {false, false, 4, 4}},
        {"r8", {true, true, 8, 4}},   {"r4", {true, true, 4, 4}},
        {"lr8", {true, false, 8, 4}}, {"lr4", {true, false, 4, 4}},
    };
    // penultimate dot-separated component, e.g. "test.lb8.ugrid" -> "lb8"
    std::vector<std::string> parts;
    std::size_t start = 0;
    for (std::size_t i = 0; i <= path.size(); ++i) {
        if (i == path.size() || path[i] == '.') {
            parts.push_back(path.substr(start, i - start));
            start = i + 1;
        }
    }
    UgridType ft;
    if (parts.size() > 1) {
        auto it = table.find(parts[parts.size() - 2]);
        if (it != table.end()) {
            ft.ascii = false;
            ft.fortran = it->second.fortran;
            ft.big_endian = it->second.big;
            ft.float_size = it->second.fs;
            ft.int_size = it->second.is;
        }
    }
    return ft;
}

// Host is assumed little-endian; swap when the file is big-endian.
inline void swap_bytes(char* p, int n) {
    for (int i = 0; i < n / 2; ++i) std::swap(p[i], p[n - 1 - i]);
}

inline std::int64_t bin_read_int(std::istream& in, int size, bool swap) {
    char buf[8];
    in.read(buf, size);
    if (in.gcount() != size) throw ReadError("UGRID: unexpected end of file");
    if (swap) swap_bytes(buf, size);
    if (size == 4) {
        std::int32_t v;
        std::memcpy(&v, buf, 4);
        return v;
    }
    std::int64_t v;
    std::memcpy(&v, buf, 8);
    return v;
}

inline double bin_read_float(std::istream& in, int size, bool swap) {
    char buf[8];
    in.read(buf, size);
    if (in.gcount() != size) throw ReadError("UGRID: unexpected end of file");
    if (swap) swap_bytes(buf, size);
    if (size == 4) {
        float v;
        std::memcpy(&v, buf, 4);
        return v;
    }
    double v;
    std::memcpy(&v, buf, 8);
    return v;
}

inline void bin_write_int(std::ostream& os, std::int64_t v, int size, bool swap) {
    char buf[8];
    if (size == 4) {
        std::int32_t t = static_cast<std::int32_t>(v);
        std::memcpy(buf, &t, 4);
    } else {
        std::memcpy(buf, &v, 8);
    }
    if (swap) swap_bytes(buf, size);
    os.write(buf, size);
}

inline void bin_write_float(std::ostream& os, double v, int size, bool swap) {
    char buf[8];
    if (size == 4) {
        float t = static_cast<float>(v);
        std::memcpy(buf, &t, 4);
    } else {
        std::memcpy(buf, &v, 8);
    }
    if (swap) swap_bytes(buf, size);
    os.write(buf, size);
}

// Volume element keywords, in UGRID write/read order, with node counts.
struct VolSpec { const char* type; int nverts; };
const VolSpec kVolume[] = {
    {"tetra", 4}, {"pyramid", 5}, {"wedge", 6}, {"hexahedron", 8},
};

}  // namespace

Mesh read_ugrid(const std::string& path) {
    UgridType ft = resolve_type(path);

    // Buffer for ascii tokenizing.
    std::ifstream in(path, std::ios::binary);
    if (!in) throw ReadError("Could not open file: " + path);

    std::string buf;
    std::size_t tok_pos = 0;
    if (ft.ascii) {
        buf.assign((std::istreambuf_iterator<char>(in)),
                   std::istreambuf_iterator<char>());
    }

    const bool swap = ft.big_endian;  // host little-endian

    auto next_token = [&]() -> std::string {
        while (tok_pos < buf.size() &&
               std::isspace(static_cast<unsigned char>(buf[tok_pos])))
            ++tok_pos;
        std::size_t s = tok_pos;
        while (tok_pos < buf.size() &&
               !std::isspace(static_cast<unsigned char>(buf[tok_pos])))
            ++tok_pos;
        if (s == tok_pos) throw ReadError("UGRID: unexpected end of file");
        return buf.substr(s, tok_pos - s);
    };
    auto next_int = [&]() -> std::int64_t {
        if (ft.ascii) return std::strtoll(next_token().c_str(), nullptr, 10);
        return bin_read_int(in, ft.int_size, swap);
    };
    auto next_float = [&]() -> double {
        if (ft.ascii) return std::strtod(next_token().c_str(), nullptr);
        return bin_read_float(in, ft.float_size, swap);
    };
    auto skip_marker = [&]() {
        if (ft.fortran) next_int();
    };

    skip_marker();
    std::int64_t counts[7];
    for (int i = 0; i < 7; ++i) counts[i] = next_int();
    skip_marker();

    const std::int64_t npoints = counts[0];
    const std::int64_t ntri = counts[1];
    const std::int64_t nquad = counts[2];

    DType fdt = (ft.float_size == 8) ? DType::Float64 : DType::Float32;
    DType idt = (ft.int_size == 8) ? DType::Int64 : DType::Int32;

    skip_marker();  // start of second Fortran record

    // Points (always 3 coordinates).
    Mesh mesh;
    mesh.points = NDArray(fdt, {static_cast<std::size_t>(npoints), 3});
    for (std::int64_t i = 0; i < npoints * 3; ++i) {
        double v = next_float();
        if (fdt == DType::Float64)
            mesh.points.as<double>()[i] = v;
        else
            mesh.points.as<float>()[i] = static_cast<float>(v);
    }

    auto store_int = [&](NDArray& a, std::int64_t i, std::int64_t v) {
        if (idt == DType::Int64)
            a.as<std::int64_t>()[i] = v;
        else
            a.as<std::int32_t>()[i] = static_cast<std::int32_t>(v);
    };

    std::vector<NDArray> refs;  // aligns with mesh.cells order

    // Surface connectivity: triangle then quad (1-based -> 0-based).
    const std::pair<const char*, int> surf[] = {{"triangle", 3}, {"quad", 4}};
    const std::int64_t surf_n[] = {ntri, nquad};
    for (int s = 0; s < 2; ++s) {
        std::int64_t n = surf_n[s];
        if (n == 0) continue;
        int k = surf[s].second;
        NDArray data(idt, {static_cast<std::size_t>(n), static_cast<std::size_t>(k)});
        for (std::int64_t i = 0; i < n * k; ++i) store_int(data, i, next_int() - 1);
        mesh.cells.emplace_back(surf[s].first, std::move(data));
    }

    // Surface boundary tags -> ugrid:ref.
    for (int s = 0; s < 2; ++s) {
        std::int64_t n = surf_n[s];
        if (n == 0) continue;
        NDArray ref(DType::Int64, {static_cast<std::size_t>(n)});
        for (std::int64_t i = 0; i < n; ++i) ref.as<std::int64_t>()[i] = next_int();
        refs.push_back(std::move(ref));
    }

    // Volume elements: tetra, pyramid (reorder), wedge, hexahedron.
    for (int vi = 0; vi < 4; ++vi) {
        std::int64_t n = counts[3 + vi];
        if (n == 0) continue;
        int k = kVolume[vi].nverts;
        NDArray data(idt, {static_cast<std::size_t>(n), static_cast<std::size_t>(k)});
        for (std::int64_t i = 0; i < n; ++i) {
            std::int64_t row[8];
            for (int j = 0; j < k; ++j) row[j] = next_int() - 1;
            if (std::string(kVolume[vi].type) == "pyramid") {
                // ugrid -> meshio: out[:, [1, 0, 3, 4, 2]]
                const int perm[5] = {1, 0, 3, 4, 2};
                for (int j = 0; j < 5; ++j)
                    store_int(data, i * 5 + j, row[perm[j]]);
            } else {
                for (int j = 0; j < k; ++j) store_int(data, i * k + j, row[j]);
            }
        }
        mesh.cells.emplace_back(kVolume[vi].type, std::move(data));
        // Volume elements carry zero ref tags.
        NDArray ref(DType::Int64, {static_cast<std::size_t>(n)});
        std::memset(ref.data(), 0, ref.nbytes());
        refs.push_back(std::move(ref));
    }

    skip_marker();  // end of second Fortran record

    mesh.cell_data.emplace("ugrid:ref", std::move(refs));
    return mesh;
}

void write_ugrid(const std::string& path, const Mesh& mesh) {
    UgridType ft = resolve_type(path);

    std::ofstream os(path, std::ios::binary);
    if (!os) throw WriteError("Could not open file for writing: " + path);

    const bool swap = ft.big_endian;

    // Resolve the single block index for each UGRID-known cell type.
    std::map<std::string, int> block_of;  // type -> index in mesh.cells
    for (std::size_t i = 0; i < mesh.cells.size(); ++i) {
        const std::string& t = mesh.cells[i].type;
        bool known = (t == "triangle" || t == "quad" || t == "tetra" ||
                      t == "pyramid" || t == "wedge" || t == "hexahedron");
        if (!known)
            throw WriteError("UGRID mesh format doesn't know " + t + " cells.");
        if (block_of.count(t))
            throw WriteError("Ugrid can only handle one cell block of a type.");
        block_of[t] = static_cast<int>(i);
    }

    auto count_of = [&](const char* t) -> std::int64_t {
        auto it = block_of.find(t);
        return it == block_of.end() ? 0 : detail::rows(mesh.cells[it->second].data);
    };

    const std::int64_t npoints = static_cast<std::int64_t>(mesh.num_points());
    std::int64_t counts[7] = {npoints,
                              count_of("triangle"),
                              count_of("quad"),
                              count_of("tetra"),
                              count_of("pyramid"),
                              count_of("wedge"),
                              count_of("hexahedron")};

    // First int cell-data array, used for surface boundary tags.
    const std::vector<NDArray>* labels = nullptr;
    for (const auto& kv : mesh.cell_data) {
        if (kv.second.empty()) continue;
        DType t = kv.second.front().dtype();
        if (t != DType::Float32 && t != DType::Float64) {
            labels = &kv.second;
            break;
        }
    }

    const std::size_t ncols =
        mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 3;

    // ---- ascii branch ----
    if (ft.ascii) {
        char fbuf[64];
        for (int i = 0; i < 7; ++i) os << counts[i] << (i == 6 ? '\n' : ' ');
        for (std::int64_t i = 0; i < npoints; ++i) {
            for (std::size_t c = 0; c < ncols; ++c) {
                std::snprintf(fbuf, sizeof(fbuf), "%.16g",
                              detail::read_double(mesh.points, i * ncols + c));
                os << fbuf << (c + 1 == ncols ? '\n' : ' ');
            }
        }
        const std::pair<const char*, int> surf[] = {{"triangle", 3}, {"quad", 4}};
        for (int s = 0; s < 2; ++s) {
            if (count_of(surf[s].first) == 0) continue;
            const CellBlock& cb = mesh.cells[block_of[surf[s].first]];
            int k = surf[s].second;
            std::int64_t n = detail::rows(cb.data);
            for (std::int64_t i = 0; i < n; ++i)
                for (int j = 0; j < k; ++j)
                    os << (detail::read_int(cb.data, i * k + j) + 1)
                       << (j + 1 == k ? '\n' : ' ');
        }
        for (int s = 0; s < 2; ++s) {
            const char* t = surf[s].first;
            std::int64_t n = count_of(t);
            if (n == 0) continue;
            int bi = block_of[t];
            const NDArray* lab =
                (labels && static_cast<std::size_t>(bi) < labels->size())
                    ? &(*labels)[bi]
                    : nullptr;
            for (std::int64_t i = 0; i < n; ++i)
                os << (lab ? detail::read_int(*lab, i) : 1) << '\n';
        }
        for (int vi = 0; vi < 4; ++vi) {
            const char* t = kVolume[vi].type;
            if (count_of(t) == 0) continue;
            const CellBlock& cb = mesh.cells[block_of[t]];
            int k = kVolume[vi].nverts;
            std::int64_t n = detail::rows(cb.data);
            for (std::int64_t i = 0; i < n; ++i) {
                if (std::string(t) == "pyramid") {
                    const int perm[5] = {1, 0, 4, 2, 3};  // meshio -> ugrid
                    for (int j = 0; j < 5; ++j)
                        os << (detail::read_int(cb.data, i * 5 + perm[j]) + 1)
                           << (j + 1 == 5 ? '\n' : ' ');
                } else {
                    for (int j = 0; j < k; ++j)
                        os << (detail::read_int(cb.data, i * k + j) + 1)
                           << (j + 1 == k ? '\n' : ' ');
                }
            }
        }
        return;
    }

    // ---- binary branch ----
    auto wint = [&](std::int64_t v) { bin_write_int(os, v, ft.int_size, swap); };
    auto wflt = [&](double v) { bin_write_float(os, v, ft.float_size, swap); };

    // Fortran record-length markers; values are not validated on read, so we
    // emit each record's nominal byte length in the file representation.
    const std::int64_t header_bytes = 7 * ft.int_size;
    std::int64_t body_bytes = npoints * 3 * ft.float_size;
    body_bytes += (counts[1] * 3 + counts[2] * 4 + counts[3] * 4 + counts[4] * 5 +
                   counts[5] * 6 + counts[6] * 8) *
                  ft.int_size;
    body_bytes += (counts[1] + counts[2]) * ft.int_size;  // surface tags

    if (ft.fortran) wint(header_bytes);
    for (int i = 0; i < 7; ++i) wint(counts[i]);
    if (ft.fortran) wint(header_bytes);

    if (ft.fortran) wint(body_bytes);
    for (std::int64_t i = 0; i < npoints; ++i)
        for (std::size_t c = 0; c < ncols; ++c)
            wflt(detail::read_double(mesh.points, i * ncols + c));

    const std::pair<const char*, int> surf[] = {{"triangle", 3}, {"quad", 4}};
    for (int s = 0; s < 2; ++s) {
        if (count_of(surf[s].first) == 0) continue;
        const CellBlock& cb = mesh.cells[block_of[surf[s].first]];
        int k = surf[s].second;
        std::int64_t n = detail::rows(cb.data);
        for (std::int64_t i = 0; i < n; ++i)
            for (int j = 0; j < k; ++j) wint(detail::read_int(cb.data, i * k + j) + 1);
    }
    for (int s = 0; s < 2; ++s) {
        const char* t = surf[s].first;
        std::int64_t n = count_of(t);
        if (n == 0) continue;
        int bi = block_of[t];
        const NDArray* lab =
            (labels && static_cast<std::size_t>(bi) < labels->size()) ? &(*labels)[bi]
                                                                      : nullptr;
        for (std::int64_t i = 0; i < n; ++i) wint(lab ? detail::read_int(*lab, i) : 1);
    }
    for (int vi = 0; vi < 4; ++vi) {
        const char* t = kVolume[vi].type;
        if (count_of(t) == 0) continue;
        const CellBlock& cb = mesh.cells[block_of[t]];
        int k = kVolume[vi].nverts;
        std::int64_t n = detail::rows(cb.data);
        for (std::int64_t i = 0; i < n; ++i) {
            if (std::string(t) == "pyramid") {
                const int perm[5] = {1, 0, 4, 2, 3};  // meshio -> ugrid
                for (int j = 0; j < 5; ++j)
                    wint(detail::read_int(cb.data, i * 5 + perm[j]) + 1);
            } else {
                for (int j = 0; j < k; ++j) wint(detail::read_int(cb.data, i * k + j) + 1);
            }
        }
    }
    if (ft.fortran) wint(body_bytes);
}

}  // namespace meshioplusplus
