#ifdef MESHIO_HAS_HDF5

#include "meshio/formats/med.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "meshio/detail/hdf5_util.hpp"
#include "meshio/detail/value_io.hpp"
#include "meshio/exceptions.hpp"
#include "meshio/types.hpp"

namespace meshio {

namespace {

const std::unordered_map<std::string, std::string>& meshio_to_med() {
    static const std::unordered_map<std::string, std::string> m = {
        {"vertex", "PO1"},   {"line", "SE2"},        {"line3", "SE3"},
        {"triangle", "TR3"}, {"triangle6", "TR6"},   {"quad", "QU4"},
        {"quad8", "QU8"},    {"tetra", "TE4"},       {"tetra10", "T10"},
        {"hexahedron", "HE8"}, {"hexahedron20", "H20"}, {"pyramid", "PY5"},
        {"pyramid13", "P13"}, {"wedge", "PE6"},      {"wedge15", "P15"}};
    return m;
}

const std::unordered_map<std::string, std::string>& med_to_meshio() {
    static const std::unordered_map<std::string, std::string> m = [] {
        std::unordered_map<std::string, std::string> out;
        for (const auto& kv : meshio_to_med()) out.emplace(kv.second, kv.first);
        return out;
    }();
    return m;
}

// Fixed-length (h5py np.bytes_-style) string attribute.
void write_attr_bytes(hid_t loc, const std::string& name, const std::string& value) {
    h5::Hid t(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(t, std::max<std::size_t>(1, value.size()));
    H5Tset_strpad(t, H5T_STR_NULLPAD);
    h5::Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    h5::Hid a(H5Acreate2(loc, name.c_str(), t, space, H5P_DEFAULT, H5P_DEFAULT),
              H5Aclose);
    if (!a.valid()) throw WriteError("MED: could not create attribute " + name);
    std::string buf = value.empty() ? std::string(1, '\0') : value;
    H5Awrite(a, t, buf.data());
}

void write_attr_double(hid_t loc, const std::string& name, double v) {
    h5::Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    h5::Hid a(H5Acreate2(loc, name.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT,
                         H5P_DEFAULT),
              H5Aclose);
    H5Awrite(a, H5T_NATIVE_DOUBLE, &v);
}

double read_attr_double(hid_t loc, const std::string& name) {
    h5::Hid a(H5Aopen(loc, name.c_str(), H5P_DEFAULT), H5Aclose);
    if (!a.valid()) throw ReadError("MED: missing attribute " + name);
    double v = 0;
    H5Aread(a, H5T_NATIVE_DOUBLE, &v);
    return v;
}

// Fortran-order (n, k) -> flat column-major buffer, applying `shift`.
NDArray flatten_f(const NDArray& a, std::int64_t shift) {
    const std::size_t n = detail::rows(a);
    const std::size_t k = detail::cols(a);
    NDArray out(a.dtype(), {n * k});
    for (std::size_t c = 0; c < k; ++c)
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t src = i * k + c, dst = c * n + i;
            if (detail::is_float_dtype(a.dtype())) {
                double v = detail::read_double(a, src);
                if (a.dtype() == DType::Float32)
                    out.as<float>()[dst] = static_cast<float>(v);
                else
                    out.as<double>()[dst] = v;
            } else {
                std::int64_t v = detail::read_int(a, src) + shift;
                switch (a.dtype()) {
                    case DType::Int32: out.as<std::int32_t>()[dst] = static_cast<std::int32_t>(v); break;
                    case DType::Int64: out.as<std::int64_t>()[dst] = v; break;
                    case DType::UInt32: out.as<std::uint32_t>()[dst] = static_cast<std::uint32_t>(v); break;
                    case DType::UInt64: out.as<std::uint64_t>()[dst] = static_cast<std::uint64_t>(v); break;
                    default: throw WriteError("MED: unexpected integer dtype");
                }
            }
        }
    return out;
}

// Flat column-major buffer -> (n, k) row-major, applying `shift` to ints.
NDArray unflatten_f(const NDArray& flat, std::size_t n, std::size_t k,
                    std::int64_t shift) {
    NDArray out(flat.dtype(), {n, k});
    for (std::size_t c = 0; c < k; ++c)
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t src = c * n + i, dst = i * k + c;
            if (detail::is_float_dtype(flat.dtype())) {
                double v = detail::read_double(flat, src);
                if (flat.dtype() == DType::Float32)
                    out.as<float>()[dst] = static_cast<float>(v);
                else
                    out.as<double>()[dst] = v;
            } else {
                std::int64_t v = detail::read_int(flat, src) + shift;
                switch (flat.dtype()) {
                    case DType::Int32: out.as<std::int32_t>()[dst] = static_cast<std::int32_t>(v); break;
                    case DType::Int64: out.as<std::int64_t>()[dst] = v; break;
                    case DType::UInt32: out.as<std::uint32_t>()[dst] = static_cast<std::uint32_t>(v); break;
                    case DType::UInt64: out.as<std::uint64_t>()[dst] = static_cast<std::uint64_t>(v); break;
                    default: throw ReadError("MED: unexpected integer dtype");
                }
            }
        }
    return out;
}

constexpr const char* kProfile = "MED_NO_PROFILE_INTERNAL";

// ---- families (point/cell tags) ----

std::map<std::int64_t, std::vector<std::string>> read_families(hid_t fas_group) {
    std::map<std::int64_t, std::vector<std::string>> families;
    for (const std::string& fam_name : h5::group_links(fas_group)) {
        h5::Hid fam = h5::open_group(fas_group, fam_name);
        std::int64_t set_id = h5::read_attr_int(fam, "NUM");
        h5::Hid gro = h5::open_group(fam, "GRO");
        std::int64_t n_subsets = h5::read_attr_int(gro, "NBR");
        NDArray nom = h5::read_dataset(gro, "NOM");  // (n_subsets, 80) int8
        std::vector<std::string> names;
        for (std::int64_t i = 0; i < n_subsets; ++i) {
            std::string s;
            for (int c = 0; c < 80; ++c) {
                char ch = static_cast<char>(detail::read_int(nom, i * 80 + c));
                if (ch == '\0') break;
                s += ch;
            }
            // strip
            std::size_t b = s.find_first_not_of(' ');
            std::size_t e = s.find_last_not_of(' ');
            names.push_back(b == std::string::npos ? std::string()
                                                   : s.substr(b, e - b + 1));
        }
        families.emplace(set_id, std::move(names));
    }
    return families;
}

void write_families(hid_t fm_group,
                    const std::map<std::int64_t, std::vector<std::string>>& tags) {
    for (const auto& kv : tags) {
        std::string fam_name = "FAM_" + std::to_string(kv.first);
        for (const auto& n : kv.second) fam_name += "_" + n;
        h5::Hid family = h5::create_group(fm_group, fam_name);
        h5::write_attr_int(family, "NUM", kv.first);
        h5::Hid gro = h5::create_group(family, "GRO");
        h5::write_attr_int(gro, "NBR", static_cast<std::int64_t>(kv.second.size()));
        // (n, 80) int8 dataset via an ARRAY datatype, like h5py's "80int8".
        hsize_t n = kv.second.size(), eighty = 80;
        h5::Hid at(H5Tarray_create2(H5T_STD_I8LE, 1, &eighty), H5Tclose);
        h5::Hid mt(H5Tarray_create2(H5T_NATIVE_INT8, 1, &eighty), H5Tclose);
        h5::Hid space(H5Screate_simple(1, &n, nullptr), H5Sclose);
        h5::Hid d(H5Dcreate2(gro, "NOM", at, space, H5P_DEFAULT, H5P_DEFAULT,
                             H5P_DEFAULT),
                  H5Dclose);
        std::vector<std::int8_t> buf(kv.second.size() * 80, 0);
        for (std::size_t i = 0; i < kv.second.size(); ++i)
            for (std::size_t c = 0; c < kv.second[i].size() && c < 80; ++c)
                buf[i * 80 + c] = static_cast<std::int8_t>(kv.second[i][c]);
        H5Dwrite(d, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    }
}

// ---- CHA fields ----

// Read one data support ("NOE" nodal or "NOE.<TYPE>"/"MAI.<TYPE>" per-cell).
NDArray read_field_support(hid_t supp_group, bool nodal) {
    std::string profile = h5::read_attr_string(supp_group, "PFL");
    if (profile != kProfile)
        throw ReadError("MED: non-default profiles handled by Python fallback");
    h5::Hid prof = h5::open_group(supp_group, profile);
    std::int64_t nbr = h5::read_attr_int(prof, "NBR");
    NDArray co = h5::read_dataset(prof, "CO");

    if (nodal) {
        std::size_t ncomp = nbr > 0 ? co.size() / static_cast<std::size_t>(nbr) : 1;
        NDArray v = unflatten_f(co, static_cast<std::size_t>(nbr), ncomp, 0);
        if (ncomp == 1) v.reshape({static_cast<std::size_t>(nbr)});
        return v;
    }
    std::int64_t nga = h5::read_attr_int(prof, "NGA");
    std::size_t n = static_cast<std::size_t>(nbr);
    std::size_t g = static_cast<std::size_t>(nga);
    std::size_t ncomp = (n * g) > 0 ? co.size() / (n * g) : 1;
    // reshape(n, nga, ncomp, order='F'): out[i,j,k] = co[i + n*j + n*g*k]
    NDArray out(co.dtype(), {n, g, ncomp});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < g; ++j)
            for (std::size_t k = 0; k < ncomp; ++k) {
                std::size_t src = i + n * j + n * g * k;
                std::size_t dst = (i * g + j) * ncomp + k;
                if (co.dtype() == DType::Float32)
                    out.as<float>()[dst] = co.as<float>()[src];
                else if (co.dtype() == DType::Float64)
                    out.as<double>()[dst] = co.as<double>()[src];
                else
                    out.as<std::int64_t>()[dst] = detail::read_int(co, src);
            }
    if (g == 1) {
        // (n, 1, ncomp) -> (n, ncomp); scalars -> (n,)
        out.reshape({n, ncomp});
        if (ncomp == 1) out.reshape({n});
    }
    return out;
}

std::string format_g(double t) {
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%g", t);
    return buf;
}

}  // namespace

Mesh read_med(const std::string& path, MedInfo& info) {
    h5::SilenceErrors silence;
    h5::Hid f = h5::open_file_read(path);

    h5::Hid ens = h5::open_group(f, "ENS_MAA");
    std::vector<std::string> meshes = h5::group_links(ens);
    if (meshes.size() != 1)
        throw ReadError("Must only contain exactly 1 mesh, found " +
                        std::to_string(meshes.size()) + ".");
    const std::string mesh_name = meshes[0];
    h5::Hid mesh_grp = h5::open_group(ens, mesh_name);

    std::int64_t dim = h5::read_attr_int(mesh_grp, "ESP");

    // Possible time-stepping indirection.
    h5::Hid data_grp;
    if (h5::exists(mesh_grp, "NOE")) {
        data_grp = std::move(mesh_grp);
    } else {
        std::vector<std::string> steps = h5::group_links(mesh_grp);
        if (steps.size() != 1)
            throw ReadError("Must only contain exactly 1 time-step, found " +
                            std::to_string(steps.size()) + ".");
        data_grp = h5::open_group(mesh_grp, steps[0]);
    }

    Mesh mesh;

    // Points
    h5::Hid noe = h5::open_group(data_grp, "NOE");
    {
        h5::Hid coo_ds(H5Dopen2(noe, "COO", H5P_DEFAULT), H5Dclose);
        if (!coo_ds.valid()) throw ReadError("MED: missing NOE/COO");
        std::int64_t n_points = h5::read_attr_int(coo_ds, "NBR");
        NDArray coo = h5::read_dataset(noe, "COO");
        mesh.points = unflatten_f(coo, static_cast<std::size_t>(n_points),
                                  static_cast<std::size_t>(dim), 0);
    }

    // Point tags
    if (h5::exists(noe, "FAM"))
        mesh.point_data.emplace("point_tags", h5::read_dataset(noe, "FAM"));

    // Families info
    h5::Hid fas = h5::exists(data_grp, "FAS") ? h5::open_group(data_grp, "FAS")
                                              : h5::Hid();
    if (!fas.valid()) {
        h5::Hid fas_root = h5::open_group(f, "FAS");
        fas = h5::open_group(fas_root, mesh_name);
    }
    if (h5::exists(fas, "NOEUD")) {
        h5::Hid noeud = h5::open_group(fas, "NOEUD");
        info.point_tags = read_families(noeud);
    }

    // Cells
    std::vector<std::string> cell_types;  // meshio names, in read order
    h5::Hid mai = h5::open_group(data_grp, "MAI");
    std::vector<NDArray> cell_tag_blocks;
    bool any_cell_tags = false;
    for (const std::string& med_type : h5::group_links(mai)) {
        auto it = med_to_meshio().find(med_type);
        if (it == med_to_meshio().end())
            throw ReadError("MED: unsupported cell type " + med_type);
        h5::Hid g = h5::open_group(mai, med_type);
        h5::Hid nod_ds(H5Dopen2(g, "NOD", H5P_DEFAULT), H5Dclose);
        if (!nod_ds.valid()) throw ReadError("MED: missing NOD for " + med_type);
        std::int64_t n_cells = h5::read_attr_int(nod_ds, "NBR");
        NDArray nod = h5::read_dataset(g, "NOD");
        std::size_t k = n_cells > 0 ? nod.size() / static_cast<std::size_t>(n_cells)
                                    : 0;
        mesh.cells.emplace_back(
            it->second, unflatten_f(nod, static_cast<std::size_t>(n_cells), k, -1));
        cell_types.push_back(it->second);

        if (h5::exists(g, "FAM")) {
            cell_tag_blocks.push_back(h5::read_dataset(g, "FAM"));
            any_cell_tags = true;
        }
    }
    if (any_cell_tags) {
        if (cell_tag_blocks.size() != mesh.cells.size())
            throw ReadError("MED: partial cell tags handled by Python fallback");
        mesh.cell_data.emplace("cell_tags", std::move(cell_tag_blocks));
    }

    if (h5::exists(fas, "ELEME")) {
        h5::Hid eleme = h5::open_group(fas, "ELEME");
        info.cell_tags = read_families(eleme);
    }

    // Fields (CHA)
    if (h5::exists(f, "CHA")) {
        if (h5::exists(f, "PROFILS"))
            throw ReadError("MED: profiles handled by Python fallback");
        h5::Hid cha = h5::open_group(f, "CHA");
        for (const std::string& fname : h5::group_links(cha)) {
            h5::Hid field = h5::open_group(cha, fname);
            if (h5::has_attr(field, "NOM")) {
                std::string nom = h5::read_attr_string(field, "NOM");
                std::vector<std::string> parts;
                std::size_t pos = 0;
                while (pos < nom.size()) {
                    while (pos < nom.size() && nom[pos] == ' ') ++pos;
                    std::size_t s = pos;
                    while (pos < nom.size() && nom[pos] != ' ') ++pos;
                    if (pos > s) parts.push_back(nom.substr(s, pos - s));
                }
                info.med_nom.push_back(std::move(parts));
            }

            std::vector<std::string> steps = h5::group_links(field);
            std::sort(steps.begin(), steps.end());
            for (std::size_t i = 0; i < steps.size(); ++i) {
                h5::Hid ts = h5::open_group(field, steps[i]);
                std::string name = fname;
                if (steps.size() > 1) {
                    double t = read_attr_double(ts, "PDT");
                    name = fname + "[" + std::to_string(i) + "] - " + format_g(t);
                }
                for (const std::string& supp : h5::group_links(ts)) {
                    h5::Hid sg = h5::open_group(ts, supp);
                    if (supp == "NOE") {
                        mesh.point_data.emplace(name, read_field_support(sg, true));
                    } else {
                        std::size_t dot = supp.find('.');
                        if (dot == std::string::npos)
                            throw ReadError("MED: unexpected field support " + supp);
                        auto it = med_to_meshio().find(supp.substr(dot + 1));
                        if (it == med_to_meshio().end())
                            throw ReadError("MED: unknown field cell type " + supp);
                        auto ct = std::find(cell_types.begin(), cell_types.end(),
                                            it->second);
                        if (ct == cell_types.end())
                            throw ReadError("MED: field on unknown cell block");
                        std::size_t idx =
                            static_cast<std::size_t>(ct - cell_types.begin());
                        auto& blocks = mesh.cell_data[name];
                        if (blocks.empty()) blocks.resize(cell_types.size());
                        blocks[idx] = read_field_support(sg, false);
                    }
                }
            }
        }
        // The conversion layer cannot represent missing per-block entries.
        for (const auto& kv : mesh.cell_data)
            for (const auto& b : kv.second)
                if (b.shape().empty() && kv.first != "cell_tags")
                    throw ReadError("MED: partial cell field handled by Python fallback");
    }

    return mesh;
}

namespace {

void write_field_data(hid_t cha, const std::string& mesh_name,
                      const std::vector<std::string>* field_name,
                      const std::string& name, const std::string& supp,
                      const NDArray& data, const std::string& med_type) {
    if (supp == "ELGA") return;  // unknown Gauss points: skipped (as in Python)

    constexpr const char* kStep = "0000000000000000000100000000000000000001";

    h5::Hid field, time_step;
    if (h5::exists(cha, name)) {
        field = h5::open_group(cha, name);
        std::vector<std::string> keys = h5::group_links(field);
        std::sort(keys.begin(), keys.end());
        time_step = h5::open_group(field, keys.back());
    } else {
        field = h5::create_group(cha, name);
        write_attr_bytes(field, "MAI", mesh_name);
        h5::write_attr_int(field, "TYP", 6);  // MED_FLOAT64
        write_attr_bytes(field, "UNI", "");
        write_attr_bytes(field, "UNT", "");
        std::size_t ncomp = data.ndim() <= 1 ? 1 : data.shape().back();
        h5::write_attr_int(field, "NCO", static_cast<std::int64_t>(ncomp));
        std::string nom;
        if (field_name) {
            for (const auto& n : *field_name) {
                char buf[20];
                std::snprintf(buf, sizeof(buf), "%-16s", n.c_str());
                nom += buf;
            }
        } else {
            nom = std::string(16, ' ');
        }
        write_attr_bytes(field, "NOM", nom);

        time_step = h5::create_group(field, kStep);
        h5::write_attr_int(time_step, "NDT", 1);
        h5::write_attr_int(time_step, "NOR", 1);
        write_attr_double(time_step, "PDT", 0.0);
        h5::write_attr_int(time_step, "RDT", -1);
        h5::write_attr_int(time_step, "ROR", -1);
    }

    std::string typ_name = supp == "NOEU" ? "NOE"
                           : supp == "ELNO" ? ("NOE." + med_type)
                                            : ("MAI." + med_type);
    h5::Hid typ = h5::create_group(time_step, typ_name);
    write_attr_bytes(typ, "GAU", "");
    write_attr_bytes(typ, "PFL", kProfile);
    h5::Hid prof = h5::create_group(typ, kProfile);
    std::size_t n = data.shape().empty() ? 0 : data.shape()[0];
    h5::write_attr_int(prof, "NBR", static_cast<std::int64_t>(n));
    std::int64_t nga = (supp == "ELNO" && data.ndim() >= 2)
                           ? static_cast<std::int64_t>(data.shape()[1])
                           : 1;
    h5::write_attr_int(prof, "NGA", nga);
    write_attr_bytes(prof, "GAU", "");

    // CO = data.flatten(order='F'): for shape (d0, d1, ..., dk) column-major.
    NDArray co(data.dtype(), {data.size()});
    const auto& shape = data.shape();
    std::size_t ndim = shape.size();
    // generic F-order flatten
    std::vector<std::size_t> idx(ndim, 0);
    for (std::size_t cnt = 0; cnt < data.size(); ++cnt) {
        // row-major source offset of current index tuple
        std::size_t src = 0;
        for (std::size_t d = 0; d < ndim; ++d) src = src * shape[d] + idx[d];
        if (detail::is_float_dtype(data.dtype())) {
            double v = detail::read_double(data, src);
            if (data.dtype() == DType::Float32)
                co.as<float>()[cnt] = static_cast<float>(v);
            else
                co.as<double>()[cnt] = v;
        } else {
            std::int64_t v = detail::read_int(data, src);
            switch (data.dtype()) {
                case DType::Int32: co.as<std::int32_t>()[cnt] = static_cast<std::int32_t>(v); break;
                case DType::Int64: co.as<std::int64_t>()[cnt] = v; break;
                case DType::UInt32: co.as<std::uint32_t>()[cnt] = static_cast<std::uint32_t>(v); break;
                case DType::UInt64: co.as<std::uint64_t>()[cnt] = static_cast<std::uint64_t>(v); break;
                default: throw WriteError("MED: unexpected field dtype");
            }
        }
        // increment index tuple in F-order (first axis fastest)
        for (std::size_t d = 0; d < ndim; ++d) {
            if (++idx[d] < shape[d]) break;
            idx[d] = 0;
        }
    }
    h5::write_dataset(prof, "CO", co);
}

}  // namespace

void write_med(const std::string& path, const Mesh& mesh, const MedInfo& info) {
    h5::SilenceErrors silence;

    // MED cannot have two blocks of the same type.
    for (std::size_t i = 0; i < mesh.cells.size(); ++i)
        for (std::size_t j = i + 1; j < mesh.cells.size(); ++j)
            if (mesh.cells[i].type == mesh.cells[j].type)
                throw WriteError(
                    "MED files cannot have two sections of the same cell type.");

    h5::Hid f = h5::create_file(path);

    h5::Hid infos = h5::create_group(f, "INFOS_GENERALES");
    h5::write_attr_int(infos, "MAJ", 3);
    h5::write_attr_int(infos, "MIN", 0);
    h5::write_attr_int(infos, "REL", 0);

    const std::string mesh_name = "mesh";
    const std::size_t dim =
        mesh.points.shape().size() >= 2 ? mesh.points.shape()[1] : 0;

    h5::Hid ens = h5::create_group(f, "ENS_MAA");
    h5::Hid med_mesh = h5::create_group(ens, mesh_name);
    h5::write_attr_int(med_mesh, "DIM", static_cast<std::int64_t>(dim));
    h5::write_attr_int(med_mesh, "ESP", static_cast<std::int64_t>(dim));
    h5::write_attr_int(med_mesh, "REP", 0);
    write_attr_bytes(med_mesh, "UNT", "");
    write_attr_bytes(med_mesh, "UNI", "");
    h5::write_attr_int(med_mesh, "SRT", 1);
    {
        const char* names[3] = {"X", "Y", "Z"};
        std::string nom;
        for (std::size_t c = 0; c < dim && c < 3; ++c) {
            char buf[20];
            std::snprintf(buf, sizeof(buf), "%-16s", names[c]);
            nom += buf;
        }
        write_attr_bytes(med_mesh, "NOM", nom);
    }
    write_attr_bytes(med_mesh, "DES", "Mesh created with meshio");
    h5::write_attr_int(med_mesh, "TYP", 0);

    h5::Hid time_step = h5::create_group(med_mesh, "-0000000000000000001-0000000000000000001");
    h5::write_attr_int(time_step, "CGT", 1);
    h5::write_attr_int(time_step, "NDT", -1);
    h5::write_attr_int(time_step, "NOR", -1);
    write_attr_double(time_step, "PDT", -1.0);

    // Points
    h5::Hid noe = h5::create_group(time_step, "NOE");
    h5::write_attr_int(noe, "CGT", 1);
    h5::write_attr_int(noe, "CGS", 1);
    write_attr_bytes(noe, "PFL", kProfile);
    {
        NDArray coo = flatten_f(mesh.points, 0);
        h5::write_dataset(noe, "COO", coo);
        h5::Hid d(H5Dopen2(noe, "COO", H5P_DEFAULT), H5Dclose);
        h5::write_attr_int(d, "CGT", 1);
        h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(mesh.num_points()));
    }
    auto pt = mesh.point_data.find("point_tags");
    if (pt != mesh.point_data.end()) {
        h5::write_dataset(noe, "FAM", pt->second);
        h5::Hid d(H5Dopen2(noe, "FAM", H5P_DEFAULT), H5Dclose);
        h5::write_attr_int(d, "CGT", 1);
        h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(mesh.num_points()));
    }

    // Cells
    h5::Hid mai = h5::create_group(time_step, "MAI");
    h5::write_attr_int(mai, "CGT", 1);
    auto ct = mesh.cell_data.find("cell_tags");
    for (std::size_t k = 0; k < mesh.cells.size(); ++k) {
        const CellBlock& cb = mesh.cells[k];
        auto it = meshio_to_med().find(cb.type);
        if (it == meshio_to_med().end())
            throw WriteError("MED: unsupported cell type " + cb.type);
        h5::Hid g = h5::create_group(mai, it->second);
        h5::write_attr_int(g, "CGT", 1);
        h5::write_attr_int(g, "CGS", 1);
        write_attr_bytes(g, "PFL", kProfile);
        NDArray nod = flatten_f(cb.data, +1);
        h5::write_dataset(g, "NOD", nod);
        {
            h5::Hid d(H5Dopen2(g, "NOD", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "CGT", 1);
            h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(cb.num_cells()));
        }
        if (ct != mesh.cell_data.end() && k < ct->second.size()) {
            h5::write_dataset(g, "FAM", ct->second[k]);
            h5::Hid d(H5Dopen2(g, "FAM", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "CGT", 1);
            h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(cb.num_cells()));
        }
    }

    // Families
    h5::Hid fas = h5::create_group(f, "FAS");
    h5::Hid families = h5::create_group(fas, mesh_name);
    h5::Hid family_zero = h5::create_group(families, "FAMILLE_ZERO");
    h5::write_attr_int(family_zero, "NUM", 0);
    if (!info.point_tags.empty()) {
        h5::Hid node = h5::create_group(families, "NOEUD");
        write_families(node, info.point_tags);
    }
    if (!info.cell_tags.empty()) {
        h5::Hid element = h5::create_group(families, "ELEME");
        write_families(element, info.cell_tags);
    }

    // Fields
    h5::Hid cha = h5::create_group(f, "CHA");
    std::size_t name_idx = 0;
    const auto& field_names = info.med_nom;

    for (const auto& kv : mesh.point_data) {
        if (kv.first == "point_tags") continue;
        const std::vector<std::string>* fn =
            (!field_names.empty() && name_idx < field_names.size())
                ? &field_names[name_idx]
                : nullptr;
        ++name_idx;
        write_field_data(cha, mesh_name, fn, kv.first, "NOEU", kv.second, "");
    }

    for (const auto& kv : mesh.cell_data) {
        if (kv.first == "cell_tags") continue;
        for (std::size_t k = 0; k < mesh.cells.size() && k < kv.second.size(); ++k) {
            const NDArray& data = kv.second[k];
            const std::string& med_type = meshio_to_med().at(mesh.cells[k].type);
            std::string supp;
            if (data.ndim() <= 2) {
                supp = "ELEM";
            } else if (data.shape()[1] ==
                       static_cast<std::size_t>(
                           num_nodes_per_cell().at(mesh.cells[k].type))) {
                supp = "ELNO";
            } else {
                supp = "ELGA";
            }
            const std::vector<std::string>* fn =
                (!field_names.empty() && name_idx < field_names.size())
                    ? &field_names[name_idx]
                    : nullptr;
            write_field_data(cha, mesh_name, fn, kv.first, supp, data, med_type);
        }
        ++name_idx;
    }
}

}  // namespace meshio

#endif  // MESHIO_HAS_HDF5
