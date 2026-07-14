#ifdef MESHIOPLUSPLUS_HAS_NETCDF

#include "meshioplusplus/formats/exodus.hpp"

#include <netcdf.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {

namespace {

void check(int status, const char* what, bool writing = false) {
    if (status != NC_NOERR) {
        std::string msg = std::string("Exodus/netCDF: ") + what + ": " +
                          nc_strerror(status);
        if (writing) throw WriteError(msg);
        throw ReadError(msg);
    }
}

const std::unordered_map<std::string, std::string>& exodus_to_meshio() {
    static const std::unordered_map<std::string, std::string> m = {
        {"SPHERE", "vertex"}, {"BEAM", "line"}, {"BEAM2", "line"},
        {"BEAM3", "line3"}, {"BAR2", "line"}, {"SHELL", "quad"},
        {"SHELL4", "quad"}, {"SHELL8", "quad8"}, {"SHELL9", "quad9"},
        {"QUAD", "quad"}, {"QUAD4", "quad"}, {"QUAD5", "quad5"},
        {"QUAD8", "quad8"}, {"QUAD9", "quad9"}, {"TRI", "triangle"},
        {"TRIANGLE", "triangle"}, {"TRI3", "triangle"}, {"TRI6", "triangle6"},
        {"TRI7", "triangle7"}, {"HEX", "hexahedron"}, {"HEXAHEDRON", "hexahedron"},
        {"HEX8", "hexahedron"}, {"HEX9", "hexahedron9"}, {"HEX20", "hexahedron20"},
        {"HEX27", "hexahedron27"}, {"TETRA", "tetra"}, {"TETRA4", "tetra4"},
        {"TET4", "tetra4"}, {"TETRA8", "tetra8"}, {"TETRA10", "tetra10"},
        {"TETRA14", "tetra14"}, {"PYRAMID", "pyramid"}, {"WEDGE", "wedge"}};
    return m;
}

// The Python reverse map is last-wins over dict order.
const std::unordered_map<std::string, std::string>& meshio_to_exodus() {
    static const std::unordered_map<std::string, std::string> m = {
        {"vertex", "SPHERE"}, {"line", "BAR2"}, {"line3", "BEAM3"},
        {"quad", "QUAD4"}, {"quad5", "QUAD5"}, {"quad8", "QUAD8"},
        {"quad9", "QUAD9"}, {"triangle", "TRI3"}, {"triangle6", "TRI6"},
        {"triangle7", "TRI7"}, {"hexahedron", "HEX8"}, {"hexahedron9", "HEX9"},
        {"hexahedron20", "HEX20"}, {"hexahedron27", "HEX27"}, {"tetra", "TETRA"},
        {"tetra4", "TET4"}, {"tetra8", "TETRA8"}, {"tetra10", "TETRA10"},
        {"tetra14", "TETRA14"}, {"pyramid", "PYRAMID"}, {"wedge", "WEDGE"}};
    return m;
}

nc_type nc_type_of(DType dt) {
    switch (dt) {
        case DType::Float32: return NC_FLOAT;
        case DType::Float64: return NC_DOUBLE;
        case DType::Int8: return NC_BYTE;
        case DType::Int16: return NC_SHORT;
        case DType::Int32: return NC_INT;
        case DType::Int64: return NC_INT64;
        case DType::UInt8: return NC_UBYTE;
        case DType::UInt16: return NC_USHORT;
        case DType::UInt32: return NC_UINT;
        case DType::UInt64: return NC_UINT64;
    }
    return NC_DOUBLE;
}

DType dtype_of(nc_type t) {
    switch (t) {
        case NC_FLOAT: return DType::Float32;
        case NC_DOUBLE: return DType::Float64;
        case NC_BYTE: return DType::Int8;
        case NC_SHORT: return DType::Int16;
        case NC_INT: return DType::Int32;
        case NC_INT64: return DType::Int64;
        case NC_UBYTE: return DType::UInt8;
        case NC_USHORT: return DType::UInt16;
        case NC_UINT: return DType::UInt32;
        case NC_UINT64: return DType::UInt64;
        default: throw ReadError("Exodus: unsupported netCDF variable type");
    }
}

// Read a whole variable (or a start/count hyperslab) into an NDArray.
NDArray read_var(int ncid, int varid, const std::vector<std::size_t>& start,
                 const std::vector<std::size_t>& count) {
    nc_type t;
    check(nc_inq_vartype(ncid, varid, &t), "inq_vartype");
    DType dt = dtype_of(t);
    std::vector<std::size_t> shape;
    for (std::size_t c : count) shape.push_back(c);
    NDArray out(dt, shape);
    if (out.size() > 0)
        check(nc_get_vara(ncid, varid, start.data(), count.data(), out.data()),
              "get_vara");
    return out;
}

std::vector<std::size_t> var_dims(int ncid, int varid) {
    int ndims;
    check(nc_inq_varndims(ncid, varid, &ndims), "inq_varndims");
    std::vector<int> dimids(ndims);
    check(nc_inq_vardimid(ncid, varid, dimids.data()), "inq_vardimid");
    std::vector<std::size_t> out;
    for (int d : dimids) {
        std::size_t len;
        check(nc_inq_dimlen(ncid, d, &len), "inq_dimlen");
        out.push_back(len);
    }
    return out;
}

// (n, len_string) char variable -> list of strings.
std::vector<std::string> read_names(int ncid, int varid) {
    std::vector<std::size_t> dims = var_dims(ncid, varid);
    std::size_t n = dims.size() >= 1 ? dims[0] : 0;
    std::size_t w = dims.size() >= 2 ? dims[1] : 0;
    std::vector<char> buf(n * w, '\0');
    if (n * w > 0) check(nc_get_var_text(ncid, varid, buf.data()), "get names");
    std::vector<std::string> out;
    for (std::size_t i = 0; i < n; ++i) {
        std::string s(buf.data() + i * w, strnlen(buf.data() + i * w, w));
        out.push_back(std::move(s));
    }
    return out;
}

// categorize() from _exodus.py: recombine <name>X/Y/Z triplets and
// <name>_R/_Z doubles.
struct Categorized {
    std::vector<std::pair<std::string, int>> single;
    std::vector<std::array<int, 2>> double_idx;
    std::vector<std::string> double_name;
    std::vector<std::array<int, 3>> triple_idx;
    std::vector<std::string> triple_name;
};

int index_of(const std::vector<std::string>& names, const std::string& s) {
    auto it = std::find(names.begin(), names.end(), s);
    return it == names.end() ? -1 : static_cast<int>(it - names.begin());
}

Categorized categorize(const std::vector<std::string>& names) {
    Categorized out;
    std::vector<bool> accounted(names.size(), false);
    for (std::size_t k = 0; k < names.size(); ++k) {
        if (accounted[k]) continue;
        const std::string& name = names[k];
        if (!name.empty() && name.back() == 'X') {
            int ix = static_cast<int>(k);
            int iy = index_of(names, name.substr(0, name.size() - 1) + "Y");
            int iz = index_of(names, name.substr(0, name.size() - 1) + "Z");
            // NB: Python checks truthiness, so index 0 counts as "not found".
            if (iy > 0 && iz > 0) {
                out.triple_idx.push_back({ix, iy, iz});
                out.triple_name.push_back(name.substr(0, name.size() - 1));
                accounted[ix] = accounted[iy] = accounted[iz] = true;
            } else {
                out.single.emplace_back(name, ix);
                accounted[ix] = true;
            }
        } else if (name.size() >= 2 && name.compare(name.size() - 2, 2, "_R") == 0) {
            int ir = static_cast<int>(k);
            int iz = index_of(names, name.substr(0, name.size() - 2) + "_Z");
            if (iz > 0) {
                out.double_idx.push_back({ir, iz});
                out.double_name.push_back(name.substr(0, name.size() - 2));
                accounted[ir] = accounted[iz] = true;
            } else {
                out.single.emplace_back(name, ir);
                accounted[ir] = true;
            }
        } else {
            out.single.emplace_back(name, static_cast<int>(k));
            accounted[k] = true;
        }
    }
    for (bool a : accounted)
        if (!a) throw ReadError("Exodus: inconsistent point data names");
    return out;
}

NDArray column_stack(const std::vector<const NDArray*>& cols) {
    std::size_t n = cols.empty() || cols[0]->shape().empty() ? 0 : cols[0]->shape()[0];
    NDArray out(cols[0]->dtype(), {n, cols.size()});
    for (std::size_t c = 0; c < cols.size(); ++c)
        for (std::size_t i = 0; i < n; ++i) {
            double v = detail::read_double(*cols[c], i);
            if (out.dtype() == DType::Float32)
                out.as<float>()[i * cols.size() + c] = static_cast<float>(v);
            else
                out.as<double>()[i * cols.size() + c] = v;
        }
    return out;
}

}  // namespace

Mesh read_exodus(const std::string& path) {
    int ncid;
    check(nc_open(path.c_str(), NC_NOWRITE, &ncid), "open");
    struct Closer {
        int id;
        ~Closer() { nc_close(id); }
    } closer{ncid};

    int nvars;
    check(nc_inq_nvars(ncid, &nvars), "inq_nvars");

    Mesh mesh;
    NDArray points_xyz;  // for coordx/y/z assembly
    std::size_t num_nodes = 0;
    {
        int dimid;
        if (nc_inq_dimid(ncid, "num_nodes", &dimid) == NC_NOERR)
            check(nc_inq_dimlen(ncid, dimid, &num_nodes), "num_nodes");
    }
    bool have_coord = false;
    points_xyz = NDArray(DType::Float64, {num_nodes, 3});

    std::vector<std::string> point_data_names, cell_data_names;
    std::map<int, NDArray> pd;                    // idx -> values (first step)
    std::map<int, std::map<int, NDArray>> cd;     // idx -> block -> values
    struct Block { std::string type; NDArray data; };
    std::vector<std::pair<int, Block>> blocks;    // connect{k} in numeric order

    for (int varid = 0; varid < nvars; ++varid) {
        char namebuf[NC_MAX_NAME + 1] = {0};
        check(nc_inq_varname(ncid, varid, namebuf), "inq_varname");
        std::string key(namebuf);
        std::vector<std::size_t> dims = var_dims(ncid, varid);

        if (key == "info_records" || key == "qa_records" || key == "ns_names" ||
            key.rfind("node_ns", 0) == 0) {
            // info + node sets live outside the conversion layer
            throw ReadError("Exodus: " + key + " handled by Python fallback");
        } else if (key.rfind("connect", 0) == 0) {
            char et[NC_MAX_NAME + 1] = {0};
            std::size_t attlen = 0;
            check(nc_inq_attlen(ncid, varid, "elem_type", &attlen), "elem_type len");
            check(nc_get_att_text(ncid, varid, "elem_type", et), "elem_type");
            std::string elem_type(et, attlen);
            std::transform(elem_type.begin(), elem_type.end(), elem_type.begin(),
                           [](unsigned char c) { return std::toupper(c); });
            auto it = exodus_to_meshio().find(elem_type);
            if (it == exodus_to_meshio().end())
                throw ReadError("Exodus: unknown element type " + elem_type);
            NDArray conn =
                read_var(ncid, varid, std::vector<std::size_t>(dims.size(), 0), dims);
            for (std::size_t i = 0; i < conn.size(); ++i) {
                switch (conn.dtype()) {
                    case DType::Int32: conn.as<std::int32_t>()[i] -= 1; break;
                    case DType::Int64: conn.as<std::int64_t>()[i] -= 1; break;
                    default: throw ReadError("Exodus: unexpected connectivity dtype");
                }
            }
            int blk = key.size() > 7 ? std::atoi(key.c_str() + 7) : 1;
            blocks.emplace_back(blk, Block{it->second, std::move(conn)});
        } else if (key == "coord") {
            NDArray coord =
                read_var(ncid, varid, std::vector<std::size_t>(dims.size(), 0), dims);
            std::size_t d = dims.size() >= 1 ? dims[0] : 0;
            std::size_t n = dims.size() >= 2 ? dims[1] : 0;
            mesh.points = NDArray(coord.dtype(), {n, d});
            for (std::size_t c = 0; c < d; ++c)
                for (std::size_t i = 0; i < n; ++i) {
                    if (coord.dtype() == DType::Float32)
                        mesh.points.as<float>()[i * d + c] = coord.as<float>()[c * n + i];
                    else
                        mesh.points.as<double>()[i * d + c] = coord.as<double>()[c * n + i];
                }
            have_coord = true;
        } else if (key == "coordx" || key == "coordy" || key == "coordz") {
            int c = key.back() - 'x';
            NDArray v =
                read_var(ncid, varid, std::vector<std::size_t>(dims.size(), 0), dims);
            for (std::size_t i = 0; i < num_nodes && i < v.size(); ++i)
                points_xyz.as<double>()[i * 3 + c] = detail::read_double(v, i);
        } else if (key == "name_nod_var") {
            point_data_names = read_names(ncid, varid);
        } else if (key.rfind("vals_nod_var", 0) == 0) {
            int idx = key.size() == 12 ? 0 : std::atoi(key.c_str() + 12) - 1;
            // dims: (time_step, ...) -> first step only
            std::vector<std::size_t> start(dims.size(), 0), count = dims;
            if (!count.empty()) count[0] = 1;
            NDArray v = read_var(ncid, varid, start, count);
            std::vector<std::size_t> shape(dims.begin() + 1, dims.end());
            v.reshape(shape);
            pd.emplace(idx, std::move(v));
        } else if (key == "name_elem_var") {
            cell_data_names = read_names(ncid, varid);
        } else if (key.rfind("vals_elem_var", 0) == 0) {
            // vals_elem_var(\d+)?(eb(\d+))?
            std::string rest = key.substr(13);
            int idx = 0, block = 0;
            std::size_t eb = rest.find("eb");
            std::string first = eb == std::string::npos ? rest : rest.substr(0, eb);
            if (!first.empty()) idx = std::atoi(first.c_str()) - 1;
            if (eb != std::string::npos) block = std::atoi(rest.c_str() + eb + 2) - 1;
            std::vector<std::size_t> start(dims.size(), 0), count = dims;
            if (!count.empty()) count[0] = 1;
            NDArray v = read_var(ncid, varid, start, count);
            std::vector<std::size_t> shape(dims.begin() + 1, dims.end());
            v.reshape(shape);
            cd[idx].emplace(block, std::move(v));
        }
        // all other variables (time_whole, coor_names, eb_prop1, ...) ignored
    }

    if (!have_coord) mesh.points = std::move(points_xyz);

    std::sort(blocks.begin(), blocks.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    for (auto& b : blocks) mesh.cells.emplace_back(b.second.type, std::move(b.second.data));

    // Point data with X/Y/Z + _R/_Z recombination.
    if (!point_data_names.empty()) {
        Categorized cat = categorize(point_data_names);
        for (const auto& kv : cat.single)
            mesh.point_data.emplace(kv.first, std::move(pd.at(kv.second)));
        for (std::size_t i = 0; i < cat.double_idx.size(); ++i)
            mesh.point_data.emplace(
                cat.double_name[i],
                column_stack({&pd.at(cat.double_idx[i][0]), &pd.at(cat.double_idx[i][1])}));
        for (std::size_t i = 0; i < cat.triple_idx.size(); ++i)
            mesh.point_data.emplace(
                cat.triple_name[i],
                column_stack({&pd.at(cat.triple_idx[i][0]), &pd.at(cat.triple_idx[i][1]),
                              &pd.at(cat.triple_idx[i][2])}));
    }

    // Cell data: concatenate blocks, then re-split by cell-block sizes.
    if (!cell_data_names.empty() && !cd.empty()) {
        std::vector<std::size_t> sizes;
        for (const auto& cb : mesh.cells) sizes.push_back(cb.num_cells());
        std::size_t name_i = 0;
        for (auto& kv : cd) {
            if (name_i >= cell_data_names.size()) break;
            const std::string& name = cell_data_names[name_i++];
            // concatenate in block order
            std::size_t total = 0;
            DType dt = kv.second.begin()->second.dtype();
            for (const auto& b : kv.second)
                total += b.second.shape().empty() ? 0 : b.second.shape()[0];
            NDArray all(dt, {total});
            std::size_t off = 0;
            for (const auto& b : kv.second) {
                std::memcpy(all.data() + off, b.second.data(), b.second.nbytes());
                off += b.second.nbytes();
            }
            // split
            std::vector<NDArray> out_blocks;
            std::size_t pos = 0;
            for (std::size_t s : sizes) {
                NDArray blk(dt, {s});
                std::memcpy(blk.data(), all.data() + pos * dtype_size(dt),
                            s * dtype_size(dt));
                pos += s;
                out_blocks.push_back(std::move(blk));
            }
            mesh.cell_data.emplace(name, std::move(out_blocks));
        }
    }

    return mesh;
}

void write_exodus(const std::string& path, const Mesh& mesh) {
    int ncid;
    check(nc_create(path.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid), "create", true);
    struct Closer {
        int id;
        ~Closer() { nc_close(id); }
    } closer{ncid};

    const std::size_t npts = mesh.num_points();
    const std::size_t pdim =
        mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;

    // global attributes
    {
        std::string title = "Created by meshio++ (C++ core)";
        check(nc_put_att_text(ncid, NC_GLOBAL, "title", title.size(), title.c_str()),
              "title", true);
        float v = 5.1f;
        check(nc_put_att_float(ncid, NC_GLOBAL, "version", NC_FLOAT, 1, &v), "version",
              true);
        check(nc_put_att_float(ncid, NC_GLOBAL, "api_version", NC_FLOAT, 1, &v),
              "api_version", true);
        long long w = 8;
        check(nc_put_att_longlong(ncid, NC_GLOBAL, "floating_point_word_size",
                                  NC_INT64, 1, &w),
              "fpws", true);
    }

    std::size_t total_elems = 0;
    for (const auto& cb : mesh.cells) total_elems += cb.num_cells();

    int d_nodes, d_dim, d_elem, d_blk, d_ns, d_str, d_line, d_four, d_time;
    check(nc_def_dim(ncid, "num_nodes", npts, &d_nodes), "def num_nodes", true);
    check(nc_def_dim(ncid, "num_dim", pdim, &d_dim), "def num_dim", true);
    check(nc_def_dim(ncid, "num_elem", total_elems, &d_elem), "def num_elem", true);
    check(nc_def_dim(ncid, "num_el_blk", mesh.cells.size(), &d_blk), "def num_el_blk",
          true);
    check(nc_def_dim(ncid, "num_node_sets", 0, &d_ns), "def num_node_sets", true);
    check(nc_def_dim(ncid, "len_string", 33, &d_str), "def len_string", true);
    check(nc_def_dim(ncid, "len_line", 81, &d_line), "def len_line", true);
    check(nc_def_dim(ncid, "four", 4, &d_four), "def four", true);
    check(nc_def_dim(ncid, "time_step", NC_UNLIMITED, &d_time), "def time_step", true);

    // dummy time step
    {
        int var;
        check(nc_def_var(ncid, "time_whole", NC_FLOAT, 1, &d_time, &var),
              "def time_whole", true);
        std::size_t start = 0, count = 1;
        float zero = 0.0f;
        check(nc_put_vara_float(ncid, var, &start, &count, &zero), "time_whole", true);
    }

    // coor_names
    {
        int dims[2] = {d_dim, d_str};
        int var;
        check(nc_def_var(ncid, "coor_names", NC_CHAR, 2, dims, &var), "coor_names",
              true);
        const char* names = "XYZ";
        for (std::size_t c = 0; c < pdim && c < 3; ++c) {
            std::size_t start[2] = {c, 0}, count[2] = {1, 1};
            check(nc_put_vara_text(ncid, var, start, count, &names[c]), "coor_names",
                  true);
        }
    }

    // coord (num_dim, num_nodes) = points^T
    {
        int dims[2] = {d_dim, d_nodes};
        int var;
        check(nc_def_var(ncid, "coord", nc_type_of(mesh.points.dtype()), 2, dims, &var),
              "def coord", true);
        NDArray t(mesh.points.dtype(), {pdim, npts});
        for (std::size_t c = 0; c < pdim; ++c)
            for (std::size_t i = 0; i < npts; ++i) {
                if (t.dtype() == DType::Float32)
                    t.as<float>()[c * npts + i] = mesh.points.as<float>()[i * pdim + c];
                else
                    t.as<double>()[c * npts + i] = mesh.points.as<double>()[i * pdim + c];
            }
        if (t.size() > 0) check(nc_put_var(ncid, var, t.data()), "coord", true);
    }

    // eb_prop1
    {
        int var;
        check(nc_def_var(ncid, "eb_prop1", NC_INT, 1, &d_blk, &var), "eb_prop1", true);
        std::vector<int> ids(mesh.cells.size());
        for (std::size_t k = 0; k < ids.size(); ++k) ids[k] = static_cast<int>(k);
        if (!ids.empty()) check(nc_put_var_int(ncid, var, ids.data()), "eb_prop1", true);
    }

    // connectivity blocks
    for (std::size_t k = 0; k < mesh.cells.size(); ++k) {
        const CellBlock& cb = mesh.cells[k];
        auto it = meshio_to_exodus().find(cb.type);
        if (it == meshio_to_exodus().end())
            throw WriteError("Exodus: unsupported cell type " + cb.type);
        std::string dim1 = "num_el_in_blk" + std::to_string(k + 1);
        std::string dim2 = "num_nod_per_el" + std::to_string(k + 1);
        int d1, d2;
        check(nc_def_dim(ncid, dim1.c_str(), cb.num_cells(), &d1), "blk dim", true);
        check(nc_def_dim(ncid, dim2.c_str(), detail::cols(cb.data), &d2), "blk dim",
              true);
        int dims[2] = {d1, d2};
        int var;
        std::string vname = "connect" + std::to_string(k + 1);
        check(nc_def_var(ncid, vname.c_str(), nc_type_of(cb.data.dtype()), 2, dims,
                         &var),
              "def connect", true);
        check(nc_put_att_text(ncid, var, "elem_type", it->second.size(),
                              it->second.c_str()),
              "elem_type", true);
        NDArray shifted(cb.data.dtype(), cb.data.shape());
        for (std::size_t i = 0; i < cb.data.size(); ++i) {
            std::int64_t v = detail::read_int(cb.data, i) + 1;
            switch (shifted.dtype()) {
                case DType::Int32: shifted.as<std::int32_t>()[i] = static_cast<std::int32_t>(v); break;
                case DType::Int64: shifted.as<std::int64_t>()[i] = v; break;
                default: throw WriteError("Exodus: unexpected connectivity dtype");
            }
        }
        if (shifted.size() > 0)
            check(nc_put_var(ncid, var, shifted.data()), "connect", true);
    }

    // point data
    if (!mesh.point_data.empty()) {
        int d_nnv;
        check(nc_def_dim(ncid, "num_nod_var", mesh.point_data.size(), &d_nnv),
              "num_nod_var", true);
        int name_var;
        {
            int dims[2] = {d_nnv, d_str};
            check(nc_def_var(ncid, "name_nod_var", NC_CHAR, 2, dims, &name_var),
                  "name_nod_var", true);
        }
        std::size_t k = 0;
        for (const auto& kv : mesh.point_data) {
            std::size_t start[2] = {k, 0};
            std::size_t count[2] = {1, std::min<std::size_t>(kv.first.size(), 33)};
            if (count[1] > 0)
                check(nc_put_vara_text(ncid, name_var, start, count, kv.first.c_str()),
                      "name_nod_var", true);

            const NDArray& data = kv.second;
            std::vector<int> dims = {d_time};
            for (std::size_t i = 0; i < data.shape().size(); ++i) {
                std::string dn = "dim_nod_var" + std::to_string(k) + std::to_string(i);
                int di;
                check(nc_def_dim(ncid, dn.c_str(), data.shape()[i], &di), "pd dim",
                      true);
                dims.push_back(di);
            }
            int var;
            std::string vname = "vals_nod_var" + std::to_string(k + 1);
            check(nc_def_var(ncid, vname.c_str(), nc_type_of(data.dtype()),
                             static_cast<int>(dims.size()), dims.data(), &var),
                  "def vals_nod_var", true);
            check(nc_def_var_fill(ncid, var, NC_NOFILL, nullptr), "nofill", true);
            std::vector<std::size_t> startv(dims.size(), 0), countv;
            countv.push_back(1);
            for (std::size_t s : data.shape()) countv.push_back(s);
            if (data.size() > 0)
                check(nc_put_vara(ncid, var, startv.data(), countv.data(), data.data()),
                      "vals_nod_var", true);
            ++k;
        }
    }

    // Node sets (point_sets) are not representable in the conversion layer;
    // the shim routes meshes with point_sets to the Python writer.
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
