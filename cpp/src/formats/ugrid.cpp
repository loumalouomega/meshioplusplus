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

// Project includes
#include "meshioplusplus/formats/ugrid.hpp"
#include "meshioplusplus/detail/byteswap.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"

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
// (bswap intrinsic — one instruction instead of a per-byte loop.)
inline void swap_bytes(char* p, int n) { detail::bswap_inplace(p, n); }

// Read exactly `n` bytes from the stream (throws on short read).
inline void read_exact(std::istream& in, char* dst, std::size_t n) {
    in.read(dst, static_cast<std::streamsize>(n));
    if (static_cast<std::size_t>(in.gcount()) != n)
        throw ReadError("UGRID: unexpected end of file");
}

// Scalar int read straight off the stream (header counts + Fortran markers
// only — every bulk section reads directly into its destination array).
inline std::int64_t stream_read_int(std::istream& in, int size, bool swap) {
    char tmp[8];
    read_exact(in, tmp, static_cast<std::size_t>(size));
    if (swap) swap_bytes(tmp, size);
    if (size == 4) {
        std::int32_t v;
        std::memcpy(&v, tmp, 4);
        return v;
    }
    std::int64_t v;
    std::memcpy(&v, tmp, 8);
    return v;
}

// Store one int/float of `size` bytes at `dst` (pre-sized output buffer).
inline void store_scalar_int(char* dst, std::int64_t v, int size, bool swap) {
    if (size == 4) {
        std::int32_t t = static_cast<std::int32_t>(v);
        std::memcpy(dst, &t, 4);
    } else {
        std::memcpy(dst, &v, 8);
    }
    if (swap) swap_bytes(dst, size);
}

inline void store_scalar_float(char* dst, double v, int size, bool swap) {
    if (size == 4) {
        float t = static_cast<float>(v);
        std::memcpy(dst, &t, 4);
    } else {
        std::memcpy(dst, &v, 8);
    }
    if (swap) swap_bytes(dst, size);
}

// Encode `count` floats from `data` into `dst` (float_size-wide, optional
// swap), one parallel pass. Verbatim memcpy when the widths match.
inline void bulk_write_floats(char* dst, const NDArray& data, std::size_t count,
                              int float_size, bool swap) {
    if (dtype_size(data.dtype()) == static_cast<std::size_t>(float_size) &&
        (data.dtype() == DType::Float64 || data.dtype() == DType::Float32)) {
        std::memcpy(dst, data.data(), count * static_cast<std::size_t>(float_size));
        if (swap)
            parallel_for_bw(count, [&](std::size_t i) {
                detail::bswap_inplace(dst + i * static_cast<std::size_t>(float_size),
                                      float_size);
            });
        return;
    }
    detail::dispatch_dtype(data.dtype(), [&]<class T>() {
        const T* s = data.as<T>();
        parallel_for_bw(count, [&](std::size_t i) {
            store_scalar_float(dst + i * static_cast<std::size_t>(float_size),
                               static_cast<double>(s[i]), float_size, swap);
        });
    });
}

// Encode an (nrows, k) integer block into `dst`: value = data[r*k + (perm ?
// perm[j] : j)] + shift, int_size-wide, optional swap. One parallel pass.
inline void bulk_write_ints(char* dst, const NDArray& data, std::size_t nrows,
                            std::size_t k, const int* perm, int int_size, bool swap,
                            std::int64_t shift) {
    detail::dispatch_dtype(data.dtype(), [&]<class T>() {
        const T* s = data.as<T>();
        parallel_for_bw(nrows, [&](std::size_t r) {
            char* row = dst + r * k * static_cast<std::size_t>(int_size);
            for (std::size_t j = 0; j < k; ++j) {
                std::size_t sc = perm ? static_cast<std::size_t>(perm[j]) : j;
                store_scalar_int(row + j * static_cast<std::size_t>(int_size),
                                 static_cast<std::int64_t>(s[r * k + sc]) + shift,
                                 int_size, swap);
            }
        });
    });
}

// Bulk-decode `count` floats (float_size bytes, little/big-endian) from buf at
// pos into dst (dtype fdt), one parallel pass. Replaces the per-value loop.
inline void bulk_read_floats(std::istream& in, std::size_t count, NDArray& dst,
                             int float_size, bool swap) {
    const std::size_t nbytes = count * static_cast<std::size_t>(float_size);
    // Fast path: the dst dtype matches float_size (true today — fdt is derived
    // from float_size), so the stream reads straight into the destination
    // array; big-endian files then get one in-place parallel bswap pass.
    if (dtype_size(dst.dtype()) == static_cast<std::size_t>(float_size)) {
        read_exact(in, reinterpret_cast<char*>(dst.data()), nbytes);
        if (swap) {
            char* d = reinterpret_cast<char*>(dst.data());
            parallel_for_bw(count, [&](std::size_t i) {
                detail::bswap_inplace(d + i * static_cast<std::size_t>(float_size),
                                      float_size);
            });
        }
    } else {
        std::vector<char> raw(nbytes);
        read_exact(in, raw.data(), nbytes);
        const char* base = raw.data();
        detail::dispatch_dtype(dst.dtype(), [&]<class T>() {
            T* d = dst.as<T>();
            parallel_for_bw(count, [&](std::size_t i) {
                char tmp[8];
                std::memcpy(tmp, base + i * float_size,
                            static_cast<std::size_t>(float_size));
                if (swap) swap_bytes(tmp, float_size);
                double v;
                if (float_size == 4) {
                    float t;
                    std::memcpy(&t, tmp, 4);
                    v = t;
                } else {
                    std::memcpy(&v, tmp, 8);
                }
                d[i] = static_cast<T>(v);
            });
        });
    }
}

// Bulk-decode a (nrows, k) integer block (int_size bytes, little/big-endian)
// from buf at pos into dst (dtype idt), applying `shift` (e.g. -1 for the
// 1-based->0-based conversion) and an optional per-row column permutation
// (dst column j <- source column perm[j]). One parallel pass over rows.
inline void bulk_read_ints(std::istream& in, std::size_t nrows, std::size_t k,
                           const int* perm, NDArray& dst, int int_size, bool swap,
                           std::int64_t shift) {
    const std::size_t total = nrows * k;
    const std::size_t nbytes = total * static_cast<std::size_t>(int_size);
    // Fast path: no column permutation AND the dst element width equals the
    // on-disk int width (true for connectivity, whose dtype is Int32/Int64 to
    // match int_size). The stream reads straight into the destination array; a
    // single parallel pass applies the byte-swap (big-endian only) and +shift.
    if (!perm && dtype_size(dst.dtype()) == static_cast<std::size_t>(int_size)) {
        read_exact(in, reinterpret_cast<char*>(dst.data()), nbytes);
        detail::dispatch_dtype(dst.dtype(), [&]<class T>() {
            T* d = dst.as<T>();
            if (swap || shift != 0)
                parallel_for_bw(total, [&](std::size_t i) {
                    if (swap) swap_bytes(reinterpret_cast<char*>(d + i), int_size);
                    d[i] = static_cast<T>(d[i] + shift);
                });
        });
    } else {
        // General strided path: dst column j <- source column perm[j] (or j),
        // with dtype conversion (e.g. surface tags are Int64 for a 4-byte file).
        std::vector<char> raw(nbytes);
        read_exact(in, raw.data(), nbytes);
        const char* base = raw.data();
        detail::dispatch_dtype(dst.dtype(), [&]<class T>() {
            T* d = dst.as<T>();
            parallel_for_bw(nrows, [&](std::size_t r) {
                for (std::size_t j = 0; j < k; ++j) {
                    std::size_t sc = perm ? static_cast<std::size_t>(perm[j]) : j;
                    char tmp[8];
                    std::memcpy(tmp, base + (r * k + sc) * int_size,
                                static_cast<std::size_t>(int_size));
                    if (swap) swap_bytes(tmp, int_size);
                    std::int64_t v;
                    if (int_size == 4) {
                        std::int32_t t;
                        std::memcpy(&t, tmp, 4);
                        v = t;
                    } else {
                        std::memcpy(&v, tmp, 8);
                    }
                    d[r * k + j] = static_cast<T>(v + shift);
                }
            });
        });
    }
}

// Volume element keywords, in UGRID write/read order, with node counts.
struct VolSpec { const char* type; int nverts; };
const VolSpec kVolume[] = {
    {"tetra", 4}, {"pyramid", 5}, {"wedge", 6}, {"hexahedron", 8},
};

}  // namespace

Mesh read_ugrid(const std::string& path) {
    UgridType ft = resolve_type(path);

    std::ifstream in(path, std::ios::binary);
    if (!in) throw ReadError("Could not open file: " + path);

    // ASCII slurps the file for tokenizing; binary streams each section
    // directly into its destination array (no whole-file intermediate).
    std::string buf;
    std::size_t tok_pos = 0;  // ascii tokenizer cursor
    if (ft.ascii) {
        in.seekg(0, std::ios::end);
        std::streamoff flen = in.tellg();
        in.seekg(0, std::ios::beg);
        if (flen > 0) {
            buf.resize(static_cast<std::size_t>(flen));
            in.read(buf.data(), flen);
        }
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
        return stream_read_int(in, ft.int_size, swap);
    };
    auto next_float = [&]() -> double {
        if (ft.ascii) return std::strtod(next_token().c_str(), nullptr);
        char tmp[8];
        read_exact(in, tmp, static_cast<std::size_t>(ft.float_size));
        if (swap) swap_bytes(tmp, ft.float_size);
        if (ft.float_size == 4) {
            float t;
            std::memcpy(&t, tmp, 4);
            return t;
        }
        double v;
        std::memcpy(&v, tmp, 8);
        return v;
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
    if (ft.ascii) {
        for (std::int64_t i = 0; i < npoints * 3; ++i) {
            double v = next_float();
            if (fdt == DType::Float64)
                mesh.points.as<double>()[i] = v;
            else
                mesh.points.as<float>()[i] = static_cast<float>(v);
        }
    } else {
        bulk_read_floats(in, static_cast<std::size_t>(npoints) * 3, mesh.points,
                         ft.float_size, swap);
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
        if (ft.ascii)
            for (std::int64_t i = 0; i < n * k; ++i) store_int(data, i, next_int() - 1);
        else
            bulk_read_ints(in, static_cast<std::size_t>(n),
                           static_cast<std::size_t>(k), nullptr, data, ft.int_size, swap, -1);
        mesh.cells.emplace_back(surf[s].first, std::move(data));
    }

    // Surface boundary tags -> ugrid:ref.
    for (int s = 0; s < 2; ++s) {
        std::int64_t n = surf_n[s];
        if (n == 0) continue;
        NDArray ref(DType::Int64, {static_cast<std::size_t>(n)});
        if (ft.ascii)
            for (std::int64_t i = 0; i < n; ++i) ref.as<std::int64_t>()[i] = next_int();
        else
            bulk_read_ints(in, static_cast<std::size_t>(n), 1, nullptr, ref,
                           ft.int_size, swap, 0);
        refs.push_back(std::move(ref));
    }

    // Volume elements: tetra, pyramid (reorder), wedge, hexahedron.
    for (int vi = 0; vi < 4; ++vi) {
        std::int64_t n = counts[3 + vi];
        if (n == 0) continue;
        int k = kVolume[vi].nverts;
        const bool is_pyramid = std::strcmp(kVolume[vi].type, "pyramid") == 0;
        // ugrid -> meshio pyramid node order: out[:, [1, 0, 3, 4, 2]].
        static const int pyramid_perm[5] = {1, 0, 3, 4, 2};
        const int* perm = is_pyramid ? pyramid_perm : nullptr;
        NDArray data(idt, {static_cast<std::size_t>(n), static_cast<std::size_t>(k)});
        if (ft.ascii) {
            for (std::int64_t i = 0; i < n; ++i) {
                std::int64_t row[8];
                for (int j = 0; j < k; ++j) row[j] = next_int() - 1;
                for (int j = 0; j < k; ++j)
                    store_int(data, i * k + j, row[perm ? perm[j] : j]);
            }
        } else {
            bulk_read_ints(in, static_cast<std::size_t>(n),
                           static_cast<std::size_t>(k), perm, data, ft.int_size, swap, -1);
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
    // Pre-size the whole file, encode each section with one parallel typed
    // pass at its computed offset, then a single os.write.
    const std::size_t is = static_cast<std::size_t>(ft.int_size);
    const std::size_t fs = static_cast<std::size_t>(ft.float_size);

    // Fortran record-length markers; values are not validated on read, so we
    // emit each record's nominal byte length in the file representation.
    const std::int64_t header_bytes = 7 * ft.int_size;
    std::int64_t body_bytes = npoints * 3 * ft.float_size;
    const std::int64_t conn_ints = counts[1] * 3 + counts[2] * 4 + counts[3] * 4 +
                                   counts[4] * 5 + counts[5] * 6 + counts[6] * 8;
    body_bytes += conn_ints * ft.int_size;
    body_bytes += (counts[1] + counts[2]) * ft.int_size;  // surface tags

    const std::size_t total_bytes =
        (ft.fortran ? 4 * is : 0) + 7 * is +
        static_cast<std::size_t>(npoints) * ncols * fs +
        static_cast<std::size_t>(conn_ints) * is +
        static_cast<std::size_t>(counts[1] + counts[2]) * is;
    std::vector<char> out(total_bytes);
    std::size_t off = 0;
    auto put_int = [&](std::int64_t v) {
        store_scalar_int(out.data() + off, v, ft.int_size, swap);
        off += is;
    };

    if (ft.fortran) put_int(header_bytes);
    for (int i = 0; i < 7; ++i) put_int(counts[i]);
    if (ft.fortran) put_int(header_bytes);

    if (ft.fortran) put_int(body_bytes);
    bulk_write_floats(out.data() + off, mesh.points,
                      static_cast<std::size_t>(npoints) * ncols, ft.float_size, swap);
    off += static_cast<std::size_t>(npoints) * ncols * fs;

    const std::pair<const char*, int> surf[] = {{"triangle", 3}, {"quad", 4}};
    for (int s = 0; s < 2; ++s) {
        std::int64_t n = count_of(surf[s].first);
        if (n == 0) continue;
        const CellBlock& cb = mesh.cells[block_of[surf[s].first]];
        const std::size_t k = static_cast<std::size_t>(surf[s].second);
        bulk_write_ints(out.data() + off, cb.data, static_cast<std::size_t>(n), k,
                        nullptr, ft.int_size, swap, +1);
        off += static_cast<std::size_t>(n) * k * is;
    }
    for (int s = 0; s < 2; ++s) {
        const char* t = surf[s].first;
        std::int64_t n = count_of(t);
        if (n == 0) continue;
        int bi = block_of[t];
        const NDArray* lab =
            (labels && static_cast<std::size_t>(bi) < labels->size()) ? &(*labels)[bi]
                                                                      : nullptr;
        char* base = out.data() + off;
        const std::size_t nz = static_cast<std::size_t>(n);
        if (lab) {
            detail::dispatch_dtype(lab->dtype(), [&]<class T>() {
                const T* sl = lab->as<T>();
                parallel_for_bw(nz, [&](std::size_t i) {
                    store_scalar_int(base + i * is, static_cast<std::int64_t>(sl[i]),
                                     ft.int_size, swap);
                });
            });
        } else {
            parallel_for_bw(nz, [&](std::size_t i) {
                store_scalar_int(base + i * is, 1, ft.int_size, swap);
            });
        }
        off += nz * is;
    }
    // meshio -> ugrid pyramid node order.
    static const int pyramid_perm_w[5] = {1, 0, 4, 2, 3};
    for (int vi = 0; vi < 4; ++vi) {
        const char* t = kVolume[vi].type;
        std::int64_t n = count_of(t);
        if (n == 0) continue;
        const CellBlock& cb = mesh.cells[block_of[t]];
        const std::size_t k = static_cast<std::size_t>(kVolume[vi].nverts);
        const int* perm = (std::strcmp(t, "pyramid") == 0) ? pyramid_perm_w : nullptr;
        bulk_write_ints(out.data() + off, cb.data, static_cast<std::size_t>(n), k, perm,
                        ft.int_size, swap, +1);
        off += static_cast<std::size_t>(n) * k * is;
    }
    if (ft.fortran) put_int(body_bytes);

    if (off != out.size())
        throw WriteError("UGRID: internal size mismatch while encoding");
    os.write(out.data(), static_cast<std::streamsize>(out.size()));
}

}  // namespace meshioplusplus
