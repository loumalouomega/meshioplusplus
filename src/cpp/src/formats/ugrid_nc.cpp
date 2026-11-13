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
// UGRID-conventions NetCDF meshes. See formats/ugrid_nc.hpp.

#ifdef MESHIOPLUSPLUS_HAS_NETCDF

// External includes
#include <netcdf.h>

// System includes
#include <algorithm>
#include <cctype>
#include <cstring>
#include <optional>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/formats/ugrid_nc.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "../detail/library_lock.hpp"
#include "../detail/netcdf_util.hpp"
#include "../detail/typed_view.hpp"
#include "library_preflight.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kUgWho = "UGRID netCDF";
constexpr const char* kUgNameAtt = "meshioplusplus_name";

void ug_check(int Status, const char* pWhat, bool Writing = false) {
    detail::nc::check(Status, pWhat, Writing, kUgWho);
}

bool ug_has_att(int Ncid, int Varid, const char* pName) {
    std::size_t len = 0;
    return nc_inq_attlen(Ncid, Varid, pName, &len) == NC_NOERR;
}

std::string ug_att(int Ncid, int Varid, const char* pName) {
    if (!ug_has_att(Ncid, Varid, pName))
        return {};
    nc_type t;
    ug_check(nc_inq_atttype(Ncid, Varid, pName, &t), pName);
    if (t != NC_CHAR)
        return {};
    return detail::nc::read_att_text(Ncid, Varid, pName, kUgWho);
}

std::optional<long long> ug_att_int(int Ncid, int Varid, const char* pName) {
    if (!ug_has_att(Ncid, Varid, pName))
        return std::nullopt;
    nc_type t;
    ug_check(nc_inq_atttype(Ncid, Varid, pName, &t), pName);
    if (t == NC_CHAR) {
        const std::string text = ug_att(Ncid, Varid, pName);
        try {
            return std::stoll(text);
        } catch (const std::exception&) {
            return std::nullopt;
        }
    }
    long long v = 0;
    ug_check(nc_get_att_longlong(Ncid, Varid, pName, &v), pName);
    return v;
}

std::vector<std::string> ug_split(const std::string& rText) {
    std::vector<std::string> out;
    std::string word;
    for (const char c : rText) {
        if (c == ' ' || c == '\t' || c == '\n') {
            if (!word.empty())
                out.push_back(word);
            word.clear();
        } else {
            word += c;
        }
    }
    if (!word.empty())
        out.push_back(word);
    return out;
}

int ug_varid(int Ncid, const std::string& rName) {
    int id = -1;
    if (nc_inq_varid(Ncid, rName.c_str(), &id) != NC_NOERR)
        throw ReadError("meshio++: ugrid_nc: the mesh names a variable '" + rName +
                        "' the file does not have");
    return id;
}

std::vector<int> ug_dimids(int Ncid, int Varid) {
    int nd = 0;
    ug_check(nc_inq_varndims(Ncid, Varid, &nd), "inq_varndims");
    std::vector<int> ids(static_cast<std::size_t>(nd));
    if (nd > 0)
        ug_check(nc_inq_vardimid(Ncid, Varid, ids.data()), "inq_vardimid");
    return ids;
}

int ug_dimid(int Ncid, const std::string& rName) {
    int id = -1;
    return nc_inq_dimid(Ncid, rName.c_str(), &id) == NC_NOERR ? id : -1;
}

/// A connectivity variable's rows (0-based, fill entries dropped) along @p Along.
std::vector<std::vector<std::int64_t>> ug_connectivity(int Ncid, const std::string& rName,
                                                       int Along, std::size_t NumNodes) {
    const int varid = ug_varid(Ncid, rName);
    const std::vector<int> dims = ug_dimids(Ncid, varid);
    if (dims.size() != 2)
        throw ReadError("meshio++: ugrid_nc: '" + rName + "' must be two-dimensional");
    const bool transposed = Along >= 0 && dims[1] == Along && dims[0] != Along;
    std::size_t a = 0, b = 0;
    ug_check(nc_inq_dimlen(Ncid, dims[0], &a), "inq_dimlen");
    ug_check(nc_inq_dimlen(Ncid, dims[1], &b), "inq_dimlen");
    std::vector<long long> raw(a * b, 0);
    if (!raw.empty())
        ug_check(nc_get_var_longlong(Ncid, varid, raw.data()), rName.c_str());
    const long long start = ug_att_int(Ncid, varid, "start_index").value_or(0);
    const std::optional<long long> fill = ug_att_int(Ncid, varid, "_FillValue");
    const std::size_t rows = transposed ? b : a, width = transposed ? a : b;
    std::vector<std::vector<std::int64_t>> out(rows);
    for (std::size_t r = 0; r < rows; ++r)
        for (std::size_t c = 0; c < width; ++c) {
            const long long v = transposed ? raw[c * b + r] : raw[r * b + c];
            if ((fill && v == *fill) || v == NC_FILL_INT || v == NC_FILL_INT64 || v < start)
                continue;
            const long long node = v - start;
            if (node >= static_cast<long long>(NumNodes))
                throw ReadError("meshio++: ugrid_nc: '" + rName + "' names node " +
                                std::to_string(v) + " of " + std::to_string(NumNodes));
            out[r].push_back(node);
        }
    return out;
}

/// An array with @p Rows rows of @p rLike's dtype and row shape, NaN (or 0) filled.
NDArray ug_padding(const NDArray& rLike, std::size_t Rows) {
    std::vector<std::size_t> shape = rLike.Shape();
    shape[0] = Rows;
    NDArray out(rLike.Dtype(), shape);
    detail::dispatch_dtype(rLike.Dtype(), [&]<class T>() {
        T* p = out.As<T>();
        const T value = std::is_floating_point_v<T> ? static_cast<T>(std::nan("")) : T{};
        std::fill(p, p + out.Size(), value);
    });
    return out;
}

/// Rows @p rIndex of @p rIn (rows along the first axis).
NDArray ug_take(const NDArray& rIn, const std::vector<std::size_t>& rIndex) {
    std::vector<std::size_t> shape = rIn.Shape();
    const std::size_t width = shape[0] ? rIn.Size() / shape[0] : 0;
    shape[0] = rIndex.size();
    NDArray out(rIn.Dtype(), shape);
    detail::dispatch_dtype(rIn.Dtype(), [&]<class T>() {
        const T* src = rIn.As<T>();
        T* dst = out.As<T>();
        for (std::size_t r = 0; r < rIndex.size(); ++r)
            std::copy_n(src + rIndex[r] * width, width, dst + r * width);
    });
    return out;
}

struct UgFile {
    int mNcid = -1;
    std::string mTopology;
    int mNodeDim = -1, mFaceDim = -1, mEdgeDim = -1, mTimeDim = -1;
    std::vector<double> mTimes;
};

/// The topology variable, its dimensions, and the time axis of its data.
UgFile ug_scan(int Ncid) {
    UgFile f;
    f.mNcid = Ncid;
    int nvars = 0;
    ug_check(nc_inq_nvars(Ncid, &nvars), "inq_nvars");
    int topo = -1;
    for (int v = 0; v < nvars && topo < 0; ++v)
        if (ug_att(Ncid, v, "cf_role") == "mesh_topology")
            topo = v;
    if (topo < 0) {
        const bool exodus = ug_dimid(Ncid, "num_nodes") >= 0 && ug_dimid(Ncid, "num_dim") >= 0;
        throw ReadError(exodus ? "meshio++: ugrid_nc: this netCDF file is an Exodus mesh; read "
                                 "it with file_format='exodus'"
                               : "meshio++: ugrid_nc: no variable has cf_role = "
                                 "\"mesh_topology\"; this netCDF file is not a UGRID mesh");
    }
    char name[NC_MAX_NAME + 1] = {0};
    ug_check(nc_inq_varname(Ncid, topo, name), "inq_varname");
    f.mTopology = name;
    if (ug_att_int(Ncid, topo, "topology_dimension").value_or(2) == 3)
        throw ReadError("meshio++: ugrid_nc: '" + f.mTopology +
                        "' is a volume (3-D) topology; meshio++ reads 1-D and 2-D UGRID meshes");
    const std::vector<std::string> coords = ug_split(ug_att(Ncid, topo, "node_coordinates"));
    if (coords.size() < 2)
        throw ReadError("meshio++: ugrid_nc: '" + f.mTopology +
                        "' needs at least two node_coordinates");
    f.mNodeDim = ug_dimids(Ncid, ug_varid(Ncid, coords[0])).at(0);
    const std::string face_dim = ug_att(Ncid, topo, "face_dimension");
    const std::string faces = ug_att(Ncid, topo, "face_node_connectivity");
    if (!faces.empty())
        f.mFaceDim = face_dim.empty() ? ug_dimids(Ncid, ug_varid(Ncid, faces)).at(0)
                                      : ug_dimid(Ncid, face_dim);
    const std::string edge_dim = ug_att(Ncid, topo, "edge_dimension");
    const std::string edges = ug_att(Ncid, topo, "edge_node_connectivity");
    if (!edges.empty())
        f.mEdgeDim = edge_dim.empty() ? ug_dimids(Ncid, ug_varid(Ncid, edges)).at(0)
                                      : ug_dimid(Ncid, edge_dim);
    // The time axis: the leading extra dimension of the first data variable that has one.
    for (int v = 0; v < nvars && f.mTimeDim < 0; ++v) {
        if (ug_att(Ncid, v, "mesh") != f.mTopology || ug_att(Ncid, v, "location").empty())
            continue;
        const std::vector<int> dims = ug_dimids(Ncid, v);
        if (dims.size() >= 2 && dims[0] != f.mNodeDim && dims[0] != f.mFaceDim &&
            dims[0] != f.mEdgeDim)
            f.mTimeDim = dims[0];
    }
    if (f.mTimeDim >= 0) {
        std::size_t n = 0;
        ug_check(nc_inq_dimlen(Ncid, f.mTimeDim, &n), "inq_dimlen");
        char dname[NC_MAX_NAME + 1] = {0};
        ug_check(nc_inq_dimname(Ncid, f.mTimeDim, dname), "inq_dimname");
        f.mTimes.assign(n, 0.0);
        int tv = -1;
        if (nc_inq_varid(Ncid, dname, &tv) == NC_NOERR && n > 0)
            ug_check(nc_get_var_double(Ncid, tv, f.mTimes.data()), "time");
        else
            for (std::size_t k = 0; k < n; ++k)
                f.mTimes[k] = static_cast<double>(k);
    }
    return f;
}

struct UgOpen {
    int mId = -1;
    explicit UgOpen(const std::string& rPath) {
        ug_check(nc_open(rPath.c_str(), NC_NOWRITE, &mId), "open");
    }
    ~UgOpen() {
        if (mId >= 0)
            nc_close(mId);
    }
    UgOpen(const UgOpen&) = delete;
    UgOpen& operator=(const UgOpen&) = delete;
};

}  // namespace

Mesh read_ugrid_nc(const std::string& rPath, const ReadOptions& rOptions) {
    detail::LibraryLock lock;
    detail::library_preflight_netcdf(rPath);
    const UgOpen file(rPath);
    const int ncid = file.mId;
    const UgFile f = ug_scan(ncid);
    const int topo = ug_varid(ncid, f.mTopology);
    const std::vector<std::string> coords = ug_split(ug_att(ncid, topo, "node_coordinates"));
    std::size_t nn = 0;
    ug_check(nc_inq_dimlen(ncid, f.mNodeDim, &nn), "inq_dimlen");
    NDArray pts(DType::Float64, {nn, 3});
    for (std::size_t c = 0; c < std::min<std::size_t>(coords.size(), 3); ++c) {
        std::vector<double> col(nn, 0.0);
        if (nn > 0)
            ug_check(nc_get_var_double(ncid, ug_varid(ncid, coords[c]), col.data()),
                     coords[c].c_str());
        for (std::size_t i = 0; i < nn; ++i)
            pts.As<double>()[i * 3 + c] = col[i];
    }
    Mesh mesh;
    mesh.AssignPoints(std::move(pts));

    // Faces grouped by node count: triangle, quad, then each larger count.
    std::vector<std::vector<std::int64_t>> faces;
    const std::string face_var = ug_att(ncid, topo, "face_node_connectivity");
    if (!face_var.empty())
        faces = ug_connectivity(ncid, face_var, f.mFaceDim, nn);
    std::map<std::size_t, std::vector<std::size_t>> by_count;
    for (std::size_t r = 0; r < faces.size(); ++r) {
        if (faces[r].size() < 3)
            throw ReadError("meshio++: ugrid_nc: face " + std::to_string(r) + " has " +
                            std::to_string(faces[r].size()) + " nodes");
        by_count[faces[r].size()].push_back(r);
    }
    std::vector<std::vector<std::size_t>> face_blocks;  // file face indices per block
    for (const auto& [count, rows] : by_count) {
        NDArray conn(DType::Int64, {rows.size(), count});
        for (std::size_t k = 0; k < rows.size(); ++k)
            std::copy(faces[rows[k]].begin(), faces[rows[k]].end(),
                      conn.As<std::int64_t>() + k * count);
        mesh.AddCellBlock(count == 3   ? "triangle"
                          : count == 4 ? "quad"
                                       : "polygon",
                          std::move(conn));
        face_blocks.push_back(rows);
    }
    std::size_t num_edges = 0;
    const std::string edge_var = ug_att(ncid, topo, "edge_node_connectivity");
    if (!edge_var.empty()) {
        const auto edges = ug_connectivity(ncid, edge_var, f.mEdgeDim, nn);
        num_edges = edges.size();
        if (num_edges > 0) {
            NDArray conn(DType::Int64, {num_edges, 2});
            for (std::size_t r = 0; r < num_edges; ++r) {
                if (edges[r].size() != 2)
                    throw ReadError("meshio++: ugrid_nc: edge " + std::to_string(r) +
                                    " does not have two nodes");
                conn.As<std::int64_t>()[r * 2] = edges[r][0];
                conn.As<std::int64_t>()[r * 2 + 1] = edges[r][1];
            }
            mesh.AddCellBlock("line", std::move(conn));
        }
    }

    const std::size_t step = rOptions.ResolveTimeStep(f.mTimes.size());
    if (!f.mTimes.empty()) {
        NDArray t(DType::Float64, {1});
        t.As<double>()[0] = f.mTimes[step];
        mesh.AddFieldData("meshio:time", std::move(t));
    }
    if (!rOptions.WantsAnyData())
        return mesh;

    // Data variables: node -> point data, face/edge -> cell data (merged by name).
    std::map<std::string, NDArray> face_data, edge_data;
    int nvars = 0;
    ug_check(nc_inq_nvars(ncid, &nvars), "inq_nvars");
    for (int v = 0; v < nvars; ++v) {
        if (ug_att(ncid, v, "mesh") != f.mTopology)
            continue;
        const std::string location = ug_att(ncid, v, "location");
        const int loc_dim = location == "node"   ? f.mNodeDim
                            : location == "face" ? f.mFaceDim
                            : location == "edge" ? f.mEdgeDim
                                                 : -1;
        char vname[NC_MAX_NAME + 1] = {0};
        ug_check(nc_inq_varname(ncid, v, vname), "inq_varname");
        std::string name = ug_att(ncid, v, kUgNameAtt);
        if (name.empty())
            name = vname;
        if (loc_dim < 0 || !rOptions.WantsArray(name))
            continue;
        const std::vector<int> dims = ug_dimids(ncid, v);
        std::vector<std::size_t> start(dims.size(), 0), count(dims.size(), 1);
        std::vector<std::size_t> shape;
        std::size_t k = 0;
        if (!dims.empty() && dims[0] == f.mTimeDim && f.mTimeDim >= 0) {
            start[0] = step;
            k = 1;
        }
        if (k >= dims.size() || dims[k] != loc_dim) {
            log::warn("ugrid_nc: '{}' is not laid out as ([time,] {}[, components]); skipped",
                      std::string(vname), location);
            continue;
        }
        for (std::size_t d = k; d < dims.size(); ++d) {
            ug_check(nc_inq_dimlen(ncid, dims[d], &count[d]), "inq_dimlen");
            shape.push_back(count[d]);
        }
        NDArray data = detail::nc::read_var(ncid, v, start, count, kUgWho, kUgWho);
        std::size_t width = 1;
        for (std::size_t d = 1; d < shape.size(); ++d)
            width *= shape[d];
        data.Reshape(width == 1 ? std::vector<std::size_t>{shape[0]}
                                : std::vector<std::size_t>{shape[0], width});
        if (location == "node")
            mesh.AddPointData(name, std::move(data));
        else if (location == "face")
            face_data.insert_or_assign(name, std::move(data));
        else
            edge_data.insert_or_assign(name, std::move(data));
    }
    std::vector<std::string> names;
    for (const auto& kv : face_data)
        names.push_back(kv.first);
    for (const auto& kv : edge_data)
        if (!face_data.count(kv.first))
            names.push_back(kv.first);
    for (const std::string& name : names) {
        const auto fit = face_data.find(name);
        const auto eit = edge_data.find(name);
        const NDArray& like = fit != face_data.end() ? fit->second : eit->second;
        std::vector<NDArray> blocks;
        for (const auto& rows : face_blocks)
            blocks.push_back(fit != face_data.end() ? ug_take(fit->second, rows)
                                                    : ug_padding(like, rows.size()));
        if (num_edges > 0)
            blocks.push_back(eit != edge_data.end() ? eit->second : ug_padding(like, num_edges));
        if (!blocks.empty())
            mesh.AddCellData(name, std::move(blocks));
    }
    return mesh;
}

MeshMetadata read_ugrid_nc_metadata(const std::string& rPath, const ReadOptions& rOptions) {
    detail::LibraryLock lock;
    MeshMetadata meta = metadata_from_mesh(read_ugrid_nc(rPath, rOptions));
    meta.mFellBackToFullRead = true;
    const UgOpen file(rPath);
    meta.mTimeValues = ug_scan(file.mId).mTimes;
    return meta;
}

namespace {

std::string ug_sanitize(const std::string& rName) {
    std::string out;
    for (const char c : rName)
        out += (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.')
                   ? c
                   : '_';
    if (out.empty() || !(std::isalpha(static_cast<unsigned char>(out[0])) || out[0] == '_'))
        out = "_" + out;
    return out;
}

/// Define one data variable over @p LocDim (and a component dimension when the array is
/// wider than one), queued in @p rPending for its values after `nc_enddef`.
void ug_define_data(int Ncid, const std::string& rName, const std::string& rVarName,
                    const std::string& rLocation, int LocDim, const NDArray& rArray,
                    std::vector<std::pair<int, NDArray>>& rPending) {
    std::vector<int> dims = {LocDim};
    const std::size_t rows = rArray.Shape().empty() ? 1 : rArray.Shape()[0];
    const std::size_t width = rows ? rArray.Size() / rows : 0;
    if (width != 1) {
        int comp = -1;
        ug_check(nc_def_dim(Ncid, ("n" + rVarName + "_comp").c_str(), width, &comp), "def_dim",
                 true);
        dims.push_back(comp);
    }
    int varid = -1;
    ug_check(nc_def_var(Ncid, rVarName.c_str(), detail::nc::nc_type_of(rArray.Dtype()),
                        static_cast<int>(dims.size()), dims.data(), &varid),
             rVarName.c_str(), true);
    ug_check(nc_put_att_text(Ncid, varid, "mesh", 4, "mesh"), "mesh", true);
    ug_check(nc_put_att_text(Ncid, varid, "location", rLocation.size(), rLocation.c_str()),
             "location", true);
    if (rVarName != rName)
        ug_check(nc_put_att_text(Ncid, varid, kUgNameAtt, rName.size(), rName.c_str()), kUgNameAtt,
                 true);
    rPending.emplace_back(varid, rArray);
}

NDArray ug_concat(const std::vector<NDArray>& rParts) {
    std::size_t rows = 0;
    for (const auto& p : rParts)
        rows += p.Shape().empty() ? 1 : p.Shape()[0];
    std::vector<std::size_t> shape = rParts.front().Shape();
    if (shape.empty())
        shape = {1};
    shape[0] = rows;
    NDArray out(rParts.front().Dtype(), shape);
    std::size_t at = 0;
    for (const auto& p : rParts) {
        if (p.Dtype() != out.Dtype())
            throw WriteError("meshio++: ugrid_nc: a cell data array changes dtype between blocks");
        std::memcpy(out.Data() + at, p.Data(), p.Nbytes());
        at += p.Nbytes();
    }
    return out;
}

}  // namespace

void write_ugrid_nc(const std::string& rPath, const Mesh& rMesh) {
    detail::LibraryLock lock;
    std::vector<std::size_t> face_blocks, line_blocks;
    std::size_t max_nodes = 0, num_faces = 0, num_edges = 0;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        const std::string type = cb.Type();
        if (type == "triangle" || type == "quad" || type == "polygon") {
            face_blocks.push_back(b);
            num_faces += cb.NumCells();
            for (std::size_t r = 0; r < cb.NumCells(); ++r)
                max_nodes = std::max(max_nodes, cb.IsRagged() ? cb.RowSize(r) : cb.NodesPerCell());
        } else if (type == "line") {
            line_blocks.push_back(b);
            num_edges += cb.NumCells();
        } else {
            throw WriteError("meshio++: ugrid_nc: a '" + type +
                             "' block has no UGRID 1-D or 2-D topology; write triangles, quads, "
                             "polygons and lines");
        }
    }
    int ncid = -1;
    ug_check(nc_create(rPath.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid), "create", true);
    struct Closer {
        int mId;
        ~Closer() { nc_close(mId); }
    } closer{ncid};
    const char* conv = "CF-1.8 UGRID-1.0";
    ug_check(nc_put_att_text(ncid, NC_GLOBAL, "Conventions", std::strlen(conv), conv),
             "Conventions", true);
    const std::size_t nn = rMesh.NumPoints();
    int node_dim = -1, face_dim = -1, max_dim = -1, edge_dim = -1, two_dim = -1;
    ug_check(nc_def_dim(ncid, "nmesh_node", nn, &node_dim), "def_dim", true);
    const detail::DoubleView pv(rMesh.Points());
    const std::size_t pd = rMesh.PointDim();
    bool with_z = false;
    for (std::size_t i = 0; pd > 2 && i < nn && !with_z; ++i)
        with_z = pv[i * pd + 2] != 0.0;
    int topo = -1;
    ug_check(nc_def_var(ncid, "mesh", NC_INT, 0, nullptr, &topo), "def_var", true);
    const auto text = [&](int Var, const char* pName, const std::string& rValue) {
        ug_check(nc_put_att_text(ncid, Var, pName, rValue.size(), rValue.c_str()), pName, true);
    };
    text(topo, "cf_role", "mesh_topology");
    text(topo, "long_name", "Topology data of the mesh");
    const int tdim = face_blocks.empty() ? 1 : 2;
    ug_check(nc_put_att_int(ncid, topo, "topology_dimension", NC_INT, 1, &tdim),
             "topology_dimension", true);
    text(topo, "node_coordinates",
         with_z ? "mesh_node_x mesh_node_y mesh_node_z" : "mesh_node_x mesh_node_y");
    std::vector<int> coord_vars;
    for (const char* axis : {"x", "y", "z"}) {
        if (std::string(axis) == "z" && !with_z)
            break;
        int v = -1;
        ug_check(nc_def_var(ncid, (std::string("mesh_node_") + axis).c_str(), NC_DOUBLE, 1,
                            &node_dim, &v),
                 "def_var", true);
        text(v, "standard_name", std::string("projection_") + axis + "_coordinate");
        coord_vars.push_back(v);
    }
    int face_var = -1, edge_var = -1;
    const int fill = -1, zero = 0;
    if (!face_blocks.empty()) {
        ug_check(nc_def_dim(ncid, "nmesh_face", num_faces, &face_dim), "def_dim", true);
        ug_check(nc_def_dim(ncid, "nmesh_face_max_nodes", max_nodes, &max_dim), "def_dim", true);
        const int fdims[2] = {face_dim, max_dim};
        ug_check(nc_def_var(ncid, "mesh_face_nodes", NC_INT, 2, fdims, &face_var), "def_var", true);
        ug_check(nc_def_var_fill(ncid, face_var, 0, &fill), "def_var_fill", true);
        text(face_var, "cf_role", "face_node_connectivity");
        ug_check(nc_put_att_int(ncid, face_var, "start_index", NC_INT, 1, &zero), "start_index",
                 true);
        text(topo, "face_node_connectivity", "mesh_face_nodes");
        text(topo, "face_dimension", "nmesh_face");
    }
    if (!line_blocks.empty()) {
        ug_check(nc_def_dim(ncid, "nmesh_edge", num_edges, &edge_dim), "def_dim", true);
        ug_check(nc_def_dim(ncid, "Two", 2, &two_dim), "def_dim", true);
        const int edims[2] = {edge_dim, two_dim};
        ug_check(nc_def_var(ncid, "mesh_edge_nodes", NC_INT, 2, edims, &edge_var), "def_var", true);
        text(edge_var, "cf_role", "edge_node_connectivity");
        ug_check(nc_put_att_int(ncid, edge_var, "start_index", NC_INT, 1, &zero), "start_index",
                 true);
        text(topo, "edge_node_connectivity", "mesh_edge_nodes");
        text(topo, "edge_dimension", "nmesh_edge");
    }

    std::vector<std::pair<int, NDArray>> pending;
    for (const std::string& name : rMesh.PointDataNames())
        ug_define_data(ncid, name, ug_sanitize(name), "node", node_dim, rMesh.PointData(name),
                       pending);
    for (const std::string& name : rMesh.CellDataNames()) {
        if (rMesh.CellDataNumBlocks(name) != rMesh.NumCellBlocks())
            continue;
        std::vector<NDArray> fparts, eparts;
        for (std::size_t b : face_blocks)
            fparts.push_back(rMesh.CellData(name, b));
        for (std::size_t b : line_blocks)
            eparts.push_back(rMesh.CellData(name, b));
        if (!fparts.empty())
            ug_define_data(ncid, name, ug_sanitize(name), "face", face_dim, ug_concat(fparts),
                           pending);
        if (!eparts.empty())
            ug_define_data(ncid, name, ug_sanitize(name) + (fparts.empty() ? "" : "_edge"), "edge",
                           edge_dim, ug_concat(eparts), pending);
    }
    ug_check(nc_enddef(ncid), "enddef", true);

    for (std::size_t c = 0; c < coord_vars.size(); ++c) {
        std::vector<double> col(nn);
        for (std::size_t i = 0; i < nn; ++i)
            col[i] = c < pd ? pv[i * pd + c] : 0.0;
        if (nn)
            ug_check(nc_put_var_double(ncid, coord_vars[c], col.data()), "put coordinates", true);
    }
    if (face_var >= 0) {
        std::vector<int> conn(num_faces * max_nodes, fill);
        std::size_t r0 = 0;
        for (std::size_t b : face_blocks) {
            const auto cb = rMesh.Cells(b);
            for (std::size_t r = 0; r < cb.NumCells(); ++r, ++r0) {
                if (cb.IsRagged()) {
                    const auto row = cb.Row(r);
                    for (std::size_t k = 0; k < cb.RowSize(r); ++k)
                        conn[r0 * max_nodes + k] = static_cast<int>(row[k]);
                } else {
                    const detail::Int64View cv(cb.Conn());
                    for (std::size_t k = 0; k < cb.NodesPerCell(); ++k)
                        conn[r0 * max_nodes + k] = static_cast<int>(cv[r * cb.NodesPerCell() + k]);
                }
            }
        }
        if (!conn.empty())
            ug_check(nc_put_var_int(ncid, face_var, conn.data()), "put faces", true);
    }
    if (edge_var >= 0) {
        std::vector<int> conn;
        for (std::size_t b : line_blocks) {
            const detail::Int64View cv(rMesh.Cells(b).Conn());
            for (std::size_t k = 0; k < rMesh.Cells(b).NumCells() * 2; ++k)
                conn.push_back(static_cast<int>(cv[k]));
        }
        if (!conn.empty())
            ug_check(nc_put_var_int(ncid, edge_var, conn.data()), "put edges", true);
    }
    for (const auto& [varid, data] : pending)
        if (data.Size() > 0)
            ug_check(nc_put_var(ncid, varid, data.Data()), "put data", true);
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
