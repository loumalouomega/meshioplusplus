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
#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// External includes
#ifdef MESHIOPLUSPLUS_HAS_EIGEN
#include <Eigen/Dense>
#endif

// Project includes
#include "meshioplusplus/formats/med.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

namespace {

const std::unordered_map<std::string, std::string>& meshio_to_med() {
    static const std::unordered_map<std::string, std::string> m = {
        {"vertex", "PO1"},     {"line", "SE2"},         {"line3", "SE3"},     {"line4", "SE4"},
        {"triangle", "TR3"},   {"triangle6", "TR6"},    {"triangle7", "TR7"}, {"quad", "QU4"},
        {"quad8", "QU8"},      {"quad9", "QU9"},        {"tetra", "TE4"},     {"tetra10", "T10"},
        {"hexahedron", "HE8"}, {"hexahedron20", "H20"}, {"pyramid", "PY5"},   {"pyramid13", "P13"},
        {"wedge", "PE6"},      {"wedge15", "P15"},      {"polygon", "POG"},   {"polygon2", "POG2"},
        {"polyhedron", "POE"}};
    return m;
}

// Quadratic 3D types share the meshio <-> MED orientation difference, but
// their permutations are not implemented; warn (like the Python reference)
// when reading or writing them unconverted.
void warn_unconverted_3d(const std::string& rCellType) {
    if (rCellType == "tetra10" || rCellType == "hexahedron20" || rCellType == "pyramid13" ||
        rCellType == "wedge15") {
        log::warn(
            "MED: orientation conversion for quadratic 3D cells '{}' is not yet "
            "implemented. These cells may be mis-oriented for MED tools (Salome, "
            "code_saturne, code_aster, etc.).",
            rCellType);
    }
}

// self-inverse meshio <-> MED node permutations (linear 3D types).
const std::unordered_map<std::string, std::vector<int>>& med_node_perm() {
    static const std::unordered_map<std::string, std::vector<int>> m = {
        {"tetra", {0, 1, 3, 2}},
        {"pyramid", {0, 3, 2, 1, 4}},
        {"wedge", {3, 4, 5, 0, 1, 2}},
        {"hexahedron", {4, 5, 6, 7, 0, 1, 2, 3}}};
    return m;
}

// (The former reorder_med_cells pass is fused into flatten_f/unflatten_f via
// their optional `perm` argument — one pass instead of two on both read+write.)

const std::unordered_map<std::string, std::string>& med_to_meshio() {
    static const std::unordered_map<std::string, std::string> m = [] {
        std::unordered_map<std::string, std::string> out;
        for (const auto& kv : meshio_to_med())
            out.emplace(kv.second, kv.first);
        return out;
    }();
    return m;
}

// Fixed-length (h5py np.bytes_-style) string attribute.
void write_attr_bytes(hid_t loc, const std::string& rName, const std::string& rValue) {
    h5::Hid t(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(t, std::max<std::size_t>(1, rValue.size()));
    H5Tset_strpad(t, H5T_STR_NULLPAD);
    h5::Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    h5::Hid a(H5Acreate2(loc, rName.c_str(), t, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw WriteError(detail::format_compat("MED: could not create attribute {}", rName));
    std::string buf = rValue.empty() ? std::string(1, '\0') : rValue;
    H5Awrite(a, t, buf.data());
}

void write_attr_double(hid_t loc, const std::string& rName, double v) {
    h5::Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    h5::Hid a(H5Acreate2(loc, rName.c_str(), H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT),
              H5Aclose);
    H5Awrite(a, H5T_NATIVE_DOUBLE, &v);
}

// Read a scalar double attribute, defaulting to 0.0 when absent.
double read_attr_double(hid_t loc, const std::string& rName) {
    if (!h5::has_attr(loc, rName))
        return 0.0;
    h5::Hid a(H5Aopen(loc, rName.c_str(), H5P_DEFAULT), H5Aclose);
    double v = 0.0;
    H5Aread(a, H5T_NATIVE_DOUBLE, &v);
    return v;
}

// Fortran-order (n, k) -> flat column-major buffer, applying `shift` to
// integer dtypes and (fused, same pass) an optional column permutation `perm`
// (the meshio->MED node reorder). Pure index transpose (memory-bandwidth bound).
NDArray flatten_f(const NDArray& rA, std::int64_t shift, const std::vector<int>* pPerm = nullptr) {
    const std::size_t n = detail::rows(rA);
    const std::size_t k = detail::cols(rA);
    const int* p = (pPerm && pPerm->size() == k) ? pPerm->data() : nullptr;
    NDArray out(rA.Dtype(), {n * k});
    detail::dispatch_dtype(rA.Dtype(), [&]<class T>() {
        const T* src = rA.As<T>();
        T* dst = out.As<T>();
#ifdef MESHIOPLUSPLUS_HAS_EIGEN
        if (!p && shift == 0) {
            // (n,k) row-major -> (n,k) col-major = Eigen storage-order convert.
            using RM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
            using CM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
            Eigen::Map<CM>(dst, n, k) = Eigen::Map<const RM>(src, n, k);
            return;
        }
#endif
        const T s = static_cast<T>(shift);
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t c = 0; c < k; ++c) {
                std::size_t sc = p ? static_cast<std::size_t>(p[c]) : c;
                if constexpr (std::is_floating_point_v<T>)
                    dst[c * n + i] = src[i * k + sc];
                else
                    dst[c * n + i] = static_cast<T>(src[i * k + sc] + s);
            }
        });
    });
    return out;
}

// Flat column-major buffer -> (n, k) row-major, applying `shift` to integer
// dtypes and (fused, in the same pass) an optional column permutation `perm`
// (the MED->meshio node reorder). Inverse transpose of flatten_f.
NDArray unflatten_f(const NDArray& rFlat, std::size_t n, std::size_t k, std::int64_t shift,
                    const std::vector<int>* pPerm = nullptr) {
    const int* p = (pPerm && pPerm->size() == k) ? pPerm->data() : nullptr;
    NDArray out(rFlat.Dtype(), {n, k});
    detail::dispatch_dtype(rFlat.Dtype(), [&]<class T>() {
        const T* src = rFlat.As<T>();
        T* dst = out.As<T>();
#ifdef MESHIOPLUSPLUS_HAS_EIGEN
        if (!p && shift == 0) {
            using RM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
            using CM = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
            Eigen::Map<RM>(dst, n, k) = Eigen::Map<const CM>(src, n, k);
            return;
        }
#endif
        const T s = static_cast<T>(shift);
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t c = 0; c < k; ++c) {
                std::size_t sc = p ? static_cast<std::size_t>(p[c]) : c;
                if constexpr (std::is_floating_point_v<T>)
                    dst[i * k + c] = src[sc * n + i];
                else
                    dst[i * k + c] = static_cast<T>(src[sc * n + i] + s);
            }
        });
    });
    return out;
}

// --- same-type cell block consolidation -------------------------------------
//
// MED cannot have two `MAI` sections of the same element type; MSH 4.1's
// canonical structure is one cell block per *entity*, so this bites every
// real 4.1 file. The Python reference does not reject -- it merges same-type
// blocks (`_med.py:811-869`), which is the compatibility baseline these
// mirror row for row: concatenate first, then do the one Fortran-order
// flatten (or FAM/NUM row concatenation) that formerly ran per block.

// Row-major (n_total, k) connectivity across `rBlockIdx`, promoting to Int64
// only when the contributing blocks disagree on dtype (mirrors numpy's
// concatenate upcasting; when they agree the merged array keeps that dtype,
// same as a single un-consolidated block would).
NDArray med_concat_conn_rows(const Mesh& rMesh, const std::vector<std::size_t>& rBlockIdx,
                             std::size_t k) {
    std::size_t n_total = 0;
    for (std::size_t bi : rBlockIdx)
        n_total += rMesh.Cells(bi).NumCells();
    DType dt = rMesh.Cells(rBlockIdx[0]).Conn().Dtype();
    for (std::size_t bi : rBlockIdx)
        if (rMesh.Cells(bi).Conn().Dtype() != dt) {
            dt = DType::Int64;
            break;
        }
    NDArray out(dt, {n_total, k});
    detail::dispatch_dtype(dt, [&]<class T>() {
        T* dst = out.As<T>();
        std::size_t row = 0;
        for (std::size_t bi : rBlockIdx) {
            const NDArray& conn = rMesh.Cells(bi).Conn();
            const std::size_t nb = rMesh.Cells(bi).NumCells();
            for (std::size_t i = 0; i < nb; ++i)
                for (std::size_t c = 0; c < k; ++c)
                    dst[(row + i) * k + c] = static_cast<T>(detail::read_int(conn, i * k + c));
            row += nb;
        }
    });
    return out;
}

// Flat 1-D concatenation of a per-block cell_data array (FAM / med:num)
// across `rBlockIdx`. Assumes a consistent dtype across the contributing
// blocks (the uniform mesh API's own invariant for a named cell_data array).
//
// Shape-preserving: `rName` may be a scalar array (cell_tags/med:num, always
// one component -- shape (n_total,)) OR an arbitrary multi-component
// cell_data array reached via write_cha_cell_field, which must keep its
// column count (shape (n_total, k)) or every component past the first is
// silently dropped.
NDArray med_concat_cell_data_rows(const Mesh& rMesh, const std::string& rName,
                                  const std::vector<std::size_t>& rBlockIdx) {
    std::size_t n_total = 0;
    for (std::size_t bi : rBlockIdx)
        n_total += detail::rows(rMesh.CellData(rName, bi));
    const NDArray& first = rMesh.CellData(rName, rBlockIdx[0]);
    const DType dt = first.Dtype();
    const std::size_t k = detail::cols(first);
    NDArray out(dt,
                k > 1 ? std::vector<std::size_t>{n_total, k} : std::vector<std::size_t>{n_total});
    detail::dispatch_dtype(dt, [&]<class T>() {
        T* dst = out.As<T>();
        std::size_t row = 0;
        for (std::size_t bi : rBlockIdx) {
            const NDArray& d = rMesh.CellData(rName, bi);
            const T* src = d.As<T>();
            const std::size_t nb = detail::rows(d);
            std::memcpy(dst + row * k, src, nb * k * sizeof(T));
            row += nb;
        }
    });
    return out;
}

// Same as med_concat_cell_data_rows, but over a standalone per-block
// NDArray vector rather than the mesh's own cell_data (used for the
// gmsh:physical/region-synthesized FAM blocks, which never reach the mesh).
NDArray med_concat_ndarray_rows(const std::vector<NDArray>& rBlocks,
                                const std::vector<std::size_t>& rBlockIdx) {
    std::size_t n_total = 0;
    for (std::size_t bi : rBlockIdx)
        n_total += detail::rows(rBlocks[bi]);
    const DType dt = rBlocks[rBlockIdx[0]].Dtype();
    NDArray out(dt, {n_total});
    detail::dispatch_dtype(dt, [&]<class T>() {
        T* dst = out.As<T>();
        std::size_t row = 0;
        for (std::size_t bi : rBlockIdx) {
            const NDArray& d = rBlocks[bi];
            const T* src = d.As<T>();
            const std::size_t nb = detail::rows(d);
            for (std::size_t i = 0; i < nb; ++i)
                dst[row + i] = src[i];
            row += nb;
        }
    });
    return out;
}

constexpr const char* kProfile = "MED_NO_PROFILE_INTERNAL";

// --- CHA (field) writing: the single-timestep common case only -------------
//
// Deferred to the Python writer (see the guard in write_med): multi-timestep
// metadata (med:step_meta), units (med:field_units), component names
// (med:nom), and the MED-4.1 optimization bitmask (LEN/LGC/LNA/... -- our own
// reader never reads these attributes, so their absence costs nothing for a
// meshio++ round-trip; it is a narrower interoperability gap with tools that
// use them, e.g. Salome/MEDCoupling). ELNO/ELGA (element-nodal/Gauss-point
// data) needs no guard at all: the uniform mesh API's cell_data is always
// (n,) or (n, k), so a 3-D "per-node-within-cell" shape cannot even be
// constructed here.

// MED's `NOM` field attribute: 16 characters per component, concatenated.
// A scalar field gets one blank 16-char slot (the historical output, kept
// byte-identical); a k>1 field gets MED's own default component spelling
// `V1..Vk`, each left-justified in 16 characters. Twin of `_med.py`'s
// `_create_component_names` + its `f"{n:<16}"` join.
std::string med_default_component_names(std::size_t ncomponents) {
    if (ncomponents <= 1)
        return std::string(16, ' ');
    std::string out;
    out.reserve(16 * ncomponents);
    for (std::size_t i = 0; i < ncomponents; ++i) {
        char buf[24];
        std::snprintf(buf, sizeof(buf), "%-16s", ("V" + std::to_string(i + 1)).c_str());
        out += buf;
    }
    return out;
}

int med_field_type_code(DType dt) {
    switch (dt) {
        case DType::Float32:
            return 4;  // MED_FLOAT32
        case DType::Float64:
            return 6;  // MED_FLOAT64
        case DType::Int32:
            return 24;  // MED_INT32
        default:
            // Every other dtype (Int8/16, UInt8/16/32/64) is widened to Int64
            // by med_field_widen() before the data is written, so 26
            // (MED_INT64) is the correct code for all of them too.
            return 26;  // MED_INT64
    }
}

// Widen any dtype MED's field TYP does not have a code for (Int8/16,
// UInt8/16/32/64) to Int64, element by element -- fields are not a hot path,
// so a fallback copy costs nothing worth optimizing away.
NDArray med_field_widen(const NDArray& rArr) {
    switch (rArr.Dtype()) {
        case DType::Float32:
        case DType::Float64:
        case DType::Int32:
        case DType::Int64:
            return rArr;
        default:
            break;
    }
    NDArray out(DType::Int64, rArr.Shape());
    for (std::size_t i = 0; i < rArr.Size(); ++i)
        out.As<std::int64_t>()[i] = detail::read_int(rArr, i);
    return out;
}

// The field group's shared attrs (name, type, component count, blank
// units/component-names -- the deferred metadata) and its one timestep
// subgroup (fixed ndt=1, nor=-1, pdt=0.0 -- the single-timestep case this
// path is scoped to). Returns the timestep group to write the actual support
// (NOE / MAI.<type>) subgroup into.
h5::Hid write_cha_field_header(hid_t cha, const std::string& rMeshName, const std::string& rName,
                               DType dt, std::size_t ncomponents) {
    h5::Hid field = h5::create_group(cha, rName);
    write_attr_bytes(field, "MAI", rMeshName);
    h5::write_attr_int(field, "TYP", med_field_type_code(dt));
    h5::write_attr_int(field, "NCO", static_cast<std::int64_t>(ncomponents));
    write_attr_bytes(field, "UNI", "");
    write_attr_bytes(field, "UNT", "");
    // MED's NOM is 16 characters PER COMPONENT, not a fixed 16. With no
    // component names to carry (the C++ Mesh has no `med:nom` channel -- see
    // write_med's comment on lenient_field_data), generate MED's own default
    // spelling `V1..Vk` for a multi-component field so a strict consumer
    // (Salome/MEDCoupling) reads k names rather than one blank. A scalar
    // field keeps the single 16-space string, so its bytes are unchanged.
    // Twinned in `_med.py`'s `_create_component_names`.
    write_attr_bytes(field, "NOM", med_default_component_names(ncomponents));

    char step_name[64];
    std::snprintf(step_name, sizeof(step_name), "%020lld%020lld", 1LL, -1LL);
    h5::Hid ts = h5::create_group(field, step_name);
    h5::write_attr_int(ts, "NDT", 1);
    h5::write_attr_int(ts, "NOR", -1);
    write_attr_double(ts, "PDT", 0.0);
    h5::write_attr_int(ts, "RDT", -1);
    h5::write_attr_int(ts, "ROR", -1);
    return ts;
}

// One "support" subgroup (NOE for nodal, MAI.<type> for a cell block):
// GAU/PFL attrs, the default-profile subgroup with NBR/NGA/GAU/CO.
void write_cha_support(hid_t ts, const std::string& rSupportName, const NDArray& rData) {
    h5::Hid typ = h5::create_group(ts, rSupportName);
    write_attr_bytes(typ, "GAU", "");
    write_attr_bytes(typ, "PFL", kProfile);
    h5::Hid profile = h5::create_group(typ, kProfile);
    h5::write_attr_int(profile, "NBR", static_cast<std::int64_t>(detail::rows(rData)));
    h5::write_attr_int(profile, "NGA", 1);
    write_attr_bytes(profile, "GAU", "");
    NDArray widened = med_field_widen(rData);
    NDArray flat = flatten_f(widened, 0);
    h5::write_dataset(profile, "CO", flat);
}

void write_cha_nodal_field(hid_t cha, const std::string& rMeshName, const std::string& rName,
                           const NDArray& rData) {
    h5::Hid ts = write_cha_field_header(cha, rMeshName, rName, rData.Dtype(), detail::cols(rData));
    write_cha_support(ts, "NOE", rData);
}

void write_cha_cell_field(hid_t cha, const std::string& rMeshName, const std::string& rName,
                          const Mesh& rMesh) {
    // The type/component count come from the first block that actually has
    // rows -- an empty mesh's block still has a dtype/shape, so this never
    // has nothing to report.
    DType dt = DType::Float64;
    std::size_t ncomponents = 1;
    bool found = false;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks() && !found; ++b) {
        const NDArray& d = rMesh.CellData(rName, b);
        if (detail::rows(d) > 0) {
            dt = d.Dtype();
            ncomponents = detail::cols(d);
            found = true;
        }
    }
    h5::Hid ts = write_cha_field_header(cha, rMeshName, rName, dt, ncomponents);

    // Group blocks by MED type first: two blocks of the same type share one
    // "MAI.<type>" support subgroup (mirrors the connectivity-writing loop's
    // own consolidation -- two calls to write_cha_support with the same
    // name would collide creating the group).
    std::vector<std::string> type_order;
    std::unordered_map<std::string, std::vector<std::size_t>> by_type;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        auto it = meshio_to_med().find(rMesh.Cells(b).Type());
        if (it == meshio_to_med().end())
            continue;  // unsupported cell type; already reported by the
                       // connectivity-writing loop, which runs first
        if (!by_type.count(it->second))
            type_order.push_back(it->second);
        by_type[it->second].push_back(b);
    }
    for (const std::string& med_type : type_order) {
        const std::vector<std::size_t>& idxs = by_type[med_type];
        std::size_t n_total = 0;
        for (std::size_t b : idxs)
            n_total += detail::rows(rMesh.CellData(rName, b));
        if (n_total == 0)
            continue;  // no contributing block has rows for this array
        write_cha_support(ts, "MAI." + med_type, med_concat_cell_data_rows(rMesh, rName, idxs));
    }
}

// ---- families (point/cell tags) ----

void read_families(hid_t fas_group, std::map<std::int64_t, std::vector<std::string>>& rFamilies,
                   std::map<std::int64_t, std::string>& rGroupNames) {
    for (const std::string& fam_name : h5::group_links(fas_group)) {
        h5::Hid fam = h5::open_group(fas_group, fam_name);
        std::int64_t set_id = h5::read_attr_int(fam, "NUM");
        rGroupNames[set_id] = fam_name;
        if (!h5::exists(fam, "GRO")) {
            rFamilies[set_id] = {};
            continue;
        }
        h5::Hid gro = h5::open_group(fam, "GRO");
        std::int64_t n_subsets = h5::read_attr_int(gro, "NBR");
        NDArray nom = h5::read_dataset(gro, "NOM");  // (n_subsets, 80) int8
        std::vector<std::string> names;
        for (std::int64_t i = 0; i < n_subsets; ++i) {
            std::string s;
            for (int c = 0; c < 80; ++c) {
                char ch = static_cast<char>(detail::read_int(nom, i * 80 + c));
                if (ch == '\0')
                    break;
                s += ch;
            }
            std::size_t b = s.find_first_not_of(' ');
            std::size_t e = s.find_last_not_of(' ');
            names.push_back(b == std::string::npos ? std::string() : s.substr(b, e - b + 1));
        }
        rFamilies.emplace(set_id, std::move(names));
    }
}

// Read a fixed-length string attribute (latin-1), stripped of spaces and NULs.
std::string read_attr_bytes(hid_t loc, const std::string& rName) {
    if (!h5::has_attr(loc, rName))
        return "";
    std::string s = h5::read_attr_string(loc, rName);
    // strip trailing NULs and surrounding spaces
    std::size_t z = s.find('\0');
    if (z != std::string::npos)
        s = s.substr(0, z);
    std::size_t b = s.find_first_not_of(' ');
    if (b == std::string::npos)
        return "";
    std::size_t e = s.find_last_not_of(' ');
    return s.substr(b, e - b + 1);
}

// Matches _write_families in _med.py: family link name from `group_names`
// (else "FAM_<id>_"), '/'->'_', capped at 64 bytes -> "FAM_<id>"; no GRO
// subgroup when the family has no named groups; GRO/NOM is an
// H5T_ARRAY{[80] char} dataset, one 80-char slot per name, space-padded.
void write_families(hid_t fm_group, const std::map<std::int64_t, std::vector<std::string>>& rTags,
                    const std::map<std::int64_t, std::string>& rGroupNames) {
    for (const auto& kv : rTags) {
        std::int64_t set_id = kv.first;
        const std::vector<std::string>& names = kv.second;
        auto git = rGroupNames.find(set_id);
        std::string gname =
            git != rGroupNames.end() ? git->second : ("FAM_" + std::to_string(set_id) + "_");
        for (char& c : gname)
            if (c == '/')
                c = '_';
        if (gname.size() > 64)
            gname = "FAM_" + std::to_string(set_id);

        h5::Hid family = h5::create_group(fm_group, gname);
        h5::write_attr_int(family, "NUM", set_id);
        if (names.empty())
            continue;

        h5::Hid gro = h5::create_group(family, "GRO");
        h5::write_attr_int(gro, "NBR", static_cast<std::int64_t>(names.size()));
        hsize_t n = names.size(), eighty = 80;
        h5::Hid at(H5Tarray_create2(H5T_STD_I8LE, 1, &eighty), H5Tclose);
        h5::Hid mt(H5Tarray_create2(H5T_NATIVE_INT8, 1, &eighty), H5Tclose);
        h5::Hid space(H5Screate_simple(1, &n, nullptr), H5Sclose);
        h5::Hid d(H5Dcreate2(gro, "NOM", at, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                  H5Dclose);
        std::vector<std::int8_t> buf(names.size() * 80, static_cast<std::int8_t>(' '));
        for (std::size_t i = 0; i < names.size(); ++i) {
            if (names[i].size() > 80)
                throw WriteError(detail::format_compat(
                    "Family name '{}' is too long for MED format (max 80 bytes).", names[i]));
            for (std::size_t c = 0; c < names[i].size(); ++c)
                buf[i * 80 + c] = static_cast<std::int8_t>(names[i][c]);
        }
        H5Dwrite(d, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    }
}

// ---- named regions <-> families (see doc/regions.md) ----
//
// Read direction: derive one Region per group *name* from the family tables
// already read into `rInfo` plus the per-point/per-cell tag id arrays
// already on `rMesh`. Mirrors `_families_to_point_sets`/`_families_to_cell_sets`
// in `_med.py` exactly, including their asymmetry: a point family with zero
// matching points is skipped entirely (no region at all), while a cell
// family always creates its named region -- even empty -- because that is
// what the Python readers have always done and C++/Python outputs must
// keep matching.
void med_attach_point_regions(Mesh& rMesh, const MedInfo& rInfo) {
    if (rInfo.mPointTags.empty() || !rMesh.HasPointData("point_tags"))
        return;
    const NDArray& fam = rMesh.PointData("point_tags");
    std::map<std::string, std::vector<std::int64_t>> by_name;
    for (const auto& kv : rInfo.mPointTags) {
        const std::int64_t fid = kv.first;
        const std::vector<std::string>& names = kv.second;
        std::vector<std::int64_t> matches;
        for (std::size_t i = 0; i < fam.Size(); ++i)
            if (detail::read_int(fam, i) == fid)
                matches.push_back(static_cast<std::int64_t>(i));
        if (matches.empty())
            continue;  // a family matching no point contributes no region --
                       // even one already seen under this name.
        for (const auto& name : names) {
            std::vector<std::int64_t>& dst = by_name[name];
            dst.insert(dst.end(), matches.begin(), matches.end());
        }
    }
    for (auto& kv : by_name) {
        NDArray arr = NDArray::Uninit(DType::Int64, {kv.second.size()});
        std::copy(kv.second.begin(), kv.second.end(), arr.As<std::int64_t>());
        rMesh.AddRegion(meshioplusplus::Region(kv.first, RegionKind::Point, std::move(arr)));
    }
}

void med_attach_cell_regions(Mesh& rMesh, const MedInfo& rInfo) {
    if (rInfo.mCellTags.empty() || !rMesh.HasCellData("cell_tags"))
        return;
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::int64_t total = detail::total_cells(bases);
    std::map<std::string, std::vector<std::int64_t>> by_name;
    // Every name named by any family gets a region, even an empty one --
    // this loop runs regardless of whether that family matches any cell.
    for (const auto& kv : rInfo.mCellTags)
        for (const auto& name : kv.second)
            by_name.try_emplace(name);
    for (const auto& kv : rInfo.mCellTags) {
        const std::int64_t fid = kv.first;
        const std::vector<std::string>& names = kv.second;
        if (names.empty())
            continue;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const NDArray& fam = rMesh.CellData("cell_tags", b);
            for (std::size_t i = 0; i < fam.Size(); ++i) {
                if (detail::read_int(fam, i) != fid)
                    continue;
                const std::int64_t g =
                    detail::block_row_to_global(bases, b, static_cast<std::int64_t>(i));
                if (g < 0 || g >= total)
                    continue;
                for (const auto& name : names)
                    by_name[name].push_back(g);
            }
        }
    }
    for (auto& kv : by_name) {
        NDArray arr = NDArray::Uninit(DType::Int64, {kv.second.size()});
        std::copy(kv.second.begin(), kv.second.end(), arr.As<std::int64_t>());
        rMesh.AddRegion(meshioplusplus::Region(kv.first, RegionKind::Cell, std::move(arr)));
    }
}

// Write direction: synthesize a per-point/per-cell family id array plus the
// family/group-name tables from Point/Cell regions -- a C++ port of
// `_ensure_med_families`'s combo logic in `_med.py`, matched step for step so
// both writers produce byte-identical FAS/FAM output for the same input
// mesh: one family per unique combination of region names a point/cell
// belongs to, ids assigned in first-encounter order scanning points (then
// cells) ascending, node families positive from +1, element families
// negative from -1. Only called when the mesh carries no native point_tags/
// cell_tags of its own -- native data always wins (see doc/regions.md), so a
// MED->MED round trip through this writer is unaffected.
bool med_point_regions_to_tags(const Mesh& rMesh, NDArray& rFamArray,
                               std::map<std::int64_t, std::vector<std::string>>& rTags,
                               std::map<std::int64_t, std::string>& rGroupNames) {
    const std::size_t n_points = rMesh.NumPoints();
    std::vector<std::set<std::string>> groups(n_points);
    bool any = false;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Point)
            continue;
        any = true;
        const std::int64_t* e = r.Entries();
        for (std::size_t k = 0; k < r.NumEntries(); ++k)
            if (e[k] >= 0 && static_cast<std::size_t>(e[k]) < n_points)
                groups[static_cast<std::size_t>(e[k])].insert(r.mName);
    }
    if (!any)
        return false;

    rFamArray = NDArray(DType::Int32, {n_points});
    std::int32_t* fam = rFamArray.As<std::int32_t>();
    std::fill(fam, fam + n_points, 0);
    std::map<std::set<std::string>, std::int64_t> combo_to_fam;
    std::int64_t next_fam = 1;  // node families: positive (MED spec)
    for (std::size_t i = 0; i < n_points; ++i) {
        if (groups[i].empty())
            continue;
        auto it = combo_to_fam.find(groups[i]);
        std::int64_t fid;
        if (it == combo_to_fam.end()) {
            fid = next_fam++;
            combo_to_fam.emplace(groups[i], fid);
            rTags[fid] = std::vector<std::string>(groups[i].begin(), groups[i].end());
            rGroupNames[fid] = "FAM_" + std::to_string(fid);
        } else {
            fid = it->second;
        }
        fam[i] = static_cast<std::int32_t>(fid);
    }
    return true;
}

bool med_cell_regions_to_tags(const Mesh& rMesh, std::vector<NDArray>& rFamBlocks,
                              std::map<std::int64_t, std::vector<std::string>>& rTags,
                              std::map<std::int64_t, std::string>& rGroupNames) {
    bool has_cell_regions = false;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        if (rMesh.Region(i).mKind == RegionKind::Cell)
            has_cell_regions = true;
    const bool has_gmsh_physical = rMesh.HasCellData("gmsh:physical");
    if (!has_cell_regions && !has_gmsh_physical)
        return false;

    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::int64_t total = detail::total_cells(bases);
    std::vector<std::set<std::string>> groups(static_cast<std::size_t>(total));

    // Named `Cell` regions -- this already covers every Gmsh 4.1 physical
    // group that was matched in $PhysicalNames (see gmsh_attach_regions in
    // gmsh.cpp, which derives exactly these regions from gmsh:physical +
    // field_data on read).
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const meshioplusplus::Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Cell)
            continue;
        const std::int64_t* e = r.Entries();
        for (std::size_t k = 0; k < r.NumEntries(); ++k)
            if (e[k] >= 0 && e[k] < total)
                groups[static_cast<std::size_t>(e[k])].insert(r.mName);
    }

    // Raw `gmsh:physical` ids, for ids NOT already exposed as a named `Cell`
    // region (avoids a duplicate group when Gmsh 4.1 carried both). This is
    // the C++ port of `_ensure_med_families`'s cell-side bridging in
    // `_med.py`, kept step for step: the readable name is `field_data`'s
    // name for that tag (first element only, no dim disambiguation -- unlike
    // `gmsh_physical_names`'s `(tag, dim)` key -- matching the Python
    // reference exactly), else `"group_<id>"`. Without this path a plain
    // Gmsh 2.2/4.0 mesh, or an un-named Gmsh 4.1 group, would silently lose
    // every physical tag on write.
    //
    // Known divergence from the Python reference: `field_data` there is
    // iterated in insertion order (last duplicate tag wins); here
    // `FieldDataNames()` is sorted, so a malformed mesh whose two group names
    // share one tag can pick a different winner. Only reachable with such a
    // mesh.
    if (has_gmsh_physical) {
        std::unordered_map<std::int64_t, std::string> id_to_name;
        for (const auto& name : rMesh.FieldDataNames()) {
            const NDArray& d = rMesh.FieldData(name);
            if (d.Size() >= 1)
                id_to_name[detail::read_int(d, 0)] = name;
        }
        const std::size_t nblocks =
            std::min(rMesh.NumCellBlocks(), rMesh.CellDataNumBlocks("gmsh:physical"));
        for (std::size_t b = 0; b < nblocks; ++b) {
            const NDArray& phys = rMesh.CellData("gmsh:physical", b);
            const std::int64_t base = bases[b];
            const std::size_t ncells = static_cast<std::size_t>(bases[b + 1] - bases[b]);
            for (std::size_t c = 0; c < ncells && c < phys.Size(); ++c) {
                const std::int64_t pid = detail::read_int(phys, c);
                if (pid == 0)
                    continue;  // family 0 -- no group
                auto nit = id_to_name.find(pid);
                std::string name;
                if (nit != id_to_name.end()) {
                    if (rMesh.HasRegion(nit->second, RegionKind::Cell))
                        continue;  // already captured as a named Cell region
                    name = nit->second;
                } else {
                    name = "group_" + std::to_string(pid);
                }
                groups[static_cast<std::size_t>(base) + c].insert(std::move(name));
            }
        }
    }

    std::map<std::set<std::string>, std::int64_t> combo_to_fam;
    std::int64_t next_fam = -1;  // element families: negative (MED spec)
    std::vector<std::int32_t> flat(static_cast<std::size_t>(total), 0);
    for (std::size_t g = 0; g < groups.size(); ++g) {
        if (groups[g].empty())
            continue;
        auto it = combo_to_fam.find(groups[g]);
        std::int64_t fid;
        if (it == combo_to_fam.end()) {
            fid = next_fam--;
            combo_to_fam.emplace(groups[g], fid);
            rTags[fid] = std::vector<std::string>(groups[g].begin(), groups[g].end());
            rGroupNames[fid] = "FAM_" + std::to_string(fid);
        } else {
            fid = it->second;
        }
        flat[g] = static_cast<std::int32_t>(fid);
    }

    rFamBlocks.reserve(rMesh.NumCellBlocks());
    for (std::size_t b = 0; b + 1 < bases.size(); ++b) {
        const std::size_t n = static_cast<std::size_t>(bases[b + 1] - bases[b]);
        NDArray block(DType::Int32, {n});
        std::int32_t* dst = block.As<std::int32_t>();
        for (std::size_t c = 0; c < n; ++c)
            dst[c] = flat[static_cast<std::size_t>(bases[b]) + c];
        rFamBlocks.push_back(std::move(block));
    }
    return true;
}

// A `Side` region has no MED equivalent (a facet is not a node or an
// element): warn and drop, like the KRATOS-backend precedent for region
// names MED-adjacent formats cannot represent structurally.
void med_warn_side_regions_dropped(const Mesh& rMesh) {
    std::size_t n = 0;
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i)
        if (rMesh.Region(i).mKind == RegionKind::Side)
            ++n;
    if (n > 0) {
        log::warn(
            "MED: {} side region(s) have no MED equivalent (a facet is not a node or an "
            "element) and were not written.",
            n);
        detail::provenance_note("regions-dropped",
                                std::to_string(n) +
                                    " side region(s) dropped -- a facet is neither a node nor "
                                    "an element in MED");
    }
}

// --- CHA (field) reading: the mirror image of write_cha_*, same scope -----
//
// Accepts only the exact shape write_med's CHA writer produces: one timestep
// group (the fixed ndt=1/nor=-1 key), the default profile, and either a
// single "NOE" (nodal) support or one-or-more "MAI.<type>" (cell) supports --
// never a mix of the two, since the writer never produces one. Anything else
// (multi-timestep, a named profile, an ELNO/ELGA support name) declines by
// throwing, exactly like the pre-existing unconditional CHA guard did, so a
// file the enhanced Python reader is needed for still gets it.

// Read one support subgroup's data ("CO" under its default-profile child),
// reshaped to `(rows,)` for a scalar field or `(rows, ncomponents)` otherwise
// -- the "1-D scalars stay 1-D" convention the rest of the core keeps.
/**
 * @brief Decline a `CHA` construct, or (under `mLenient`) skip it.
 *
 * Mirrors `mdpa_reject_or_skip`. Returns `false` when the caller should drop
 * the field and carry on, and only ever returns at all in lenient mode --
 * strict always throws, which is what keeps the Python shim falling back and
 * therefore the Python surface unchanged.
 */
void med_reject_or_skip(const std::string& rWhat, bool Lenient, MedInfo* pInfo) {
    if (!Lenient)
        throw ReadError("MED: " + rWhat +
                        " is handled by Python fallback (set ReadOptions::mLenient to skip it "
                        "instead)");
    log::warn("MED: skipping {} (ReadOptions::mLenient)", rWhat);
    if (pInfo)
        pInfo->mSkippedConstructs.push_back(rWhat);
}

/// Whether a support subgroup names a real (non-default) profile, i.e. its
/// `CO` covers a subset of the entities indexed by a separate list this reader
/// does not resolve. The caller decides whether that is fatal or a skip.
bool med_support_has_named_profile(hid_t support) {
    const std::string pfl = read_attr_bytes(support, "PFL");
    return !pfl.empty() && pfl != kProfile;
}

NDArray read_cha_support_data(hid_t support, std::int64_t ncomponents, std::size_t rows) {
    h5::Hid profile = h5::open_group(support, kProfile);
    NDArray flat = h5::read_dataset(profile, "CO");
    const std::size_t k = ncomponents > 0 ? static_cast<std::size_t>(ncomponents) : 1;
    if (flat.Size() != rows * k)
        throw ReadError("MED: field data size does not match its declared shape");
    NDArray out = unflatten_f(flat, rows, k, 0);
    if (k == 1)
        out.Reshape({rows});
    return out;
}

void read_cha_fields(hid_t cha, Mesh& rMesh, bool Lenient, const ReadOptions& rOptions,
                     MedInfo* pInfo) {
    std::unordered_map<std::string, std::size_t> med_to_block;
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        auto it = meshio_to_med().find(rMesh.Cells(b).Type());
        if (it != meshio_to_med().end())
            med_to_block.emplace(it->second, b);
    }

    for (const std::string& field_name : h5::group_links(cha)) {
        h5::Hid field = h5::open_group(cha, field_name);
        const std::int64_t ncomponents =
            h5::has_attr(field, "NCO") ? h5::read_attr_int(field, "NCO") : 1;

        // Units are a real piece of information the C++ Mesh cannot carry
        // (they are strings; `med:field_units` is a Python-only dict-valued
        // field_data convention). Strict declines rather than silently
        // reading the file without them; lenient reads them into MedInfo, so
        // nothing is lost even though nothing lands on the Mesh.
        const std::string uni = read_attr_bytes(field, "UNI");
        const std::string unt = read_attr_bytes(field, "UNT");
        if (!uni.empty() || !unt.empty()) {
            if (!Lenient)
                throw ReadError("MED: field '" + field_name +
                                "' declares units, handled by Python fallback (set "
                                "ReadOptions::mLenient to read it and report the units in "
                                "MedInfo::mFieldUnits instead)");
            log::warn(
                "MED: field '{}' declares units ('{}'/'{}'); reported in "
                "MedInfo::mFieldUnits (ReadOptions::mLenient)",
                field_name, uni, unt);
            if (pInfo)
                pInfo->mFieldUnits.emplace(field_name, std::make_pair(uni, unt));
        }

        // Every step this field carries, sorted by its group name -- which is
        // the zero-padded (NDT, NOR) pair, so name order IS step order. Their
        // PDTs go into MedInfo unconditionally, so a caller can see what
        // `mTimeStep` may select before asking for it.
        std::vector<std::string> steps = h5::group_links(field);
        std::sort(steps.begin(), steps.end());
        if (steps.empty())
            continue;  // a field group with no step carries nothing
        if (pInfo) {
            std::vector<double> times;
            times.reserve(steps.size());
            for (const std::string& s : steps) {
                h5::Hid g = h5::open_group(field, s);
                times.push_back(read_attr_double(g, "PDT"));
            }
            pInfo->mFieldTimeValues.emplace(field_name, std::move(times));
        }

        // Which step to read. An explicit non-default `mTimeStep` is honoured
        // whether or not `mLenient` is set: it is a request, not a fallback,
        // and no Python behaviour depends on it (the shim never passes one).
        // With the default step, a multi-step field is still a strict decline.
        std::size_t step_index = 0;
        if (steps.size() != 1) {
            if (rOptions.mTimeStep != 0) {
                step_index = rOptions.ResolveTimeStep(steps.size());
            } else if (!Lenient) {
                throw ReadError("MED: multi-timestep field '" + field_name + "' has " +
                                std::to_string(steps.size()) +
                                " steps and is handled by Python fallback (set "
                                "ReadOptions::mTimeStep to select one, or "
                                "ReadOptions::mLenient to take the first)");
            } else {
                log::warn(
                    "MED: field '{}' has {} timesteps; reading the first "
                    "(ReadOptions::mLenient -- set mTimeStep to choose)",
                    field_name, steps.size());
                if (pInfo)
                    pInfo->mSkippedConstructs.push_back("field '" + field_name + "' timesteps 2.." +
                                                        std::to_string(steps.size()));
            }
        } else if (rOptions.mTimeStep != 0) {
            step_index = rOptions.ResolveTimeStep(steps.size());
        }
        h5::Hid ts = h5::open_group(field, steps[step_index]);

        // A step whose NDT/NOR/PDT are not the write-side default (1/-1/0)
        // carries real information; silently reporting ndt=1 for a genuinely
        // tagged step 7 would be a wrong answer, not just an incomplete one.
        // Strict declines; lenient records it in MedInfo::mStepMeta.
        const std::int64_t ndt = h5::read_attr_int(ts, "NDT");
        const std::int64_t nor = h5::read_attr_int(ts, "NOR");
        const double pdt = read_attr_double(ts, "PDT");
        if (ndt != 1 || nor != -1 || pdt != 0.0) {
            if (!Lenient && rOptions.mTimeStep == 0)
                throw ReadError("MED: field '" + field_name +
                                "' has non-default timestep metadata, handled by Python fallback "
                                "(set ReadOptions::mLenient to read it and report the metadata in "
                                "MedInfo::mStepMeta instead)");
            if (pInfo)
                pInfo->mStepMeta.emplace(field_name, std::make_tuple(ndt, nor, pdt));
        }

        std::vector<std::string> supports = h5::group_links(ts);
        const bool is_nodal = std::find(supports.begin(), supports.end(), "NOE") != supports.end();

        if (is_nodal) {
            if (supports.size() != 1) {
                med_reject_or_skip("field '" + field_name + "' mixes nodal and cell support",
                                   Lenient, pInfo);
                continue;
            }
            h5::Hid noe = h5::open_group(ts, "NOE");
            if (med_support_has_named_profile(noe)) {
                med_reject_or_skip("field '" + field_name + "' on a named profile", Lenient, pInfo);
                continue;
            }
            rMesh.AddPointData(field_name,
                               read_cha_support_data(noe, ncomponents, rMesh.NumPoints()));
            continue;
        }

        std::vector<NDArray> per_block;
        per_block.reserve(rMesh.NumCellBlocks());
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b)
            per_block.emplace_back(DType::Float64, std::vector<std::size_t>{0});
        std::vector<bool> filled(rMesh.NumCellBlocks(), false);
        bool skip_field = false;
        for (const std::string& supp : supports) {
            // Anything not "NOE"/"MAI.<type>" is an ELNO/ELGA support: one
            // value per node-within-cell or per Gauss point, a 3-D shape the
            // uniform mesh API's (n,)/(n,k) cell_data cannot hold at all --
            // structurally unrepresentable, not merely unimplemented.
            if (supp.rfind("MAI.", 0) != 0) {
                med_reject_or_skip("field '" + field_name + "' support '" + supp +
                                       "' (ELNO/ELGA data has no (n,)/(n,k) representation)",
                                   Lenient, pInfo);
                skip_field = true;
                break;
            }
            const std::string med_type = supp.substr(4);
            auto bit = med_to_block.find(med_type);
            if (bit == med_to_block.end()) {
                med_reject_or_skip("field '" + field_name + "' names cell type '" + med_type +
                                       "' with no matching block",
                                   Lenient, pInfo);
                skip_field = true;
                break;
            }
            const std::size_t b = bit->second;
            h5::Hid grp = h5::open_group(ts, supp);
            if (med_support_has_named_profile(grp)) {
                med_reject_or_skip("field '" + field_name + "' on a named profile", Lenient, pInfo);
                skip_field = true;
                break;
            }
            per_block[b] = read_cha_support_data(grp, ncomponents, rMesh.Cells(b).NumCells());
            filled[b] = true;
        }
        if (skip_field)
            continue;
        // A block with no support subgroup carries no data for this field --
        // AddCellData still needs exactly one array per block, so it gets an
        // appropriately-shaped empty/zero one (a documented gap: the reader
        // cannot tell "genuinely zero cells' worth of data" apart from "this
        // block's rows were simply never written", so it reports zeros).
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            if (filled[b])
                continue;
            const std::size_t ncells = rMesh.Cells(b).NumCells();
            std::vector<std::size_t> shape =
                ncomponents > 1
                    ? std::vector<std::size_t>{ncells, static_cast<std::size_t>(ncomponents)}
                    : std::vector<std::size_t>{ncells};
            per_block[b] = NDArray(DType::Float64, shape);
        }
        rMesh.AddCellData(field_name, std::move(per_block));
    }
}

}  // namespace

// The single implementation behind both `read_med` overloads. `rOptions`
// default-constructed reproduces the historical strict behaviour exactly,
// which is what let the options overload land without touching any caller.
namespace {
Mesh med_read_impl(const std::string& rPath, MedInfo& rInfo, const ReadOptions& rOptions) {
    h5::SilenceErrors silence;
    h5::Hid f = h5::open_file_read(rPath);

    // MED data-model version: a file written by a MED major version newer
    // than the 4.x layout this reader implements gets a clear diagnosis
    // instead of an obscure "missing NOE/COO"-style structural error further
    // down. Older majors (the repo's own fixtures include MED 3.x files)
    // read exactly as before -- only *newer* is rejected.
    if (h5::exists(f, "INFOS_GENERALES")) {
        h5::Hid infos = h5::open_group(f, "INFOS_GENERALES");
        if (h5::has_attr(infos, "MAJ")) {
            std::int64_t maj = h5::read_attr_int(infos, "MAJ");
            if (maj > 4) {
                std::int64_t min = h5::has_attr(infos, "MIN") ? h5::read_attr_int(infos, "MIN") : 0;
                std::int64_t rel = h5::has_attr(infos, "REL") ? h5::read_attr_int(infos, "REL") : 0;
                throw ReadError(detail::format_compat(
                    "MED file '{}' was written by MED {}.{}.{}, newer than the MED 4.1 "
                    "data model this reader implements",
                    rPath, maj, min, rel));
            }
        }
    }

    h5::Hid ens = h5::open_group(f, "ENS_MAA");
    std::vector<std::string> meshes = h5::group_links(ens);
    if (meshes.size() != 1)
        throw ReadError(
            detail::format_compat("Must only contain exactly 1 mesh, found {}.", meshes.size()));
    const std::string mesh_name = meshes[0];
    h5::Hid mesh_grp = h5::open_group(ens, mesh_name);

    std::int64_t dim = h5::read_attr_int(mesh_grp, "ESP");

    // Mesh-level metadata attributes.
    rInfo.mMeshName = mesh_name;
    rInfo.mDescription = read_attr_bytes(mesh_grp, "DES");
    rInfo.mUnitTime = read_attr_bytes(mesh_grp, "UNT");
    rInfo.mUnitCoords = read_attr_bytes(mesh_grp, "UNI");

    // Possible time-stepping indirection.
    h5::Hid data_grp;
    if (h5::exists(mesh_grp, "NOE")) {
        data_grp = std::move(mesh_grp);
    } else {
        std::vector<std::string> steps = h5::group_links(mesh_grp);
        if (steps.size() != 1)
            throw ReadError(detail::format_compat(
                "Must only contain exactly 1 time-step, found {}.", steps.size()));
        data_grp = h5::open_group(mesh_grp, steps[0]);
    }

    Mesh mesh;

    // Points
    h5::Hid noe = h5::open_group(data_grp, "NOE");
    {
        h5::Hid coo_ds(H5Dopen2(noe, "COO", H5P_DEFAULT), H5Dclose);
        if (!coo_ds.Valid())
            throw ReadError("MED: missing NOE/COO");
        std::int64_t n_points = h5::read_attr_int(coo_ds, "NBR");
        NDArray coo = h5::read_dataset(noe, "COO");
        mesh.AssignPoints(
            unflatten_f(coo, static_cast<std::size_t>(n_points), static_cast<std::size_t>(dim), 0));
    }

    // Point tags
    if (h5::exists(noe, "FAM"))
        mesh.AddPointData("point_tags", h5::read_dataset(noe, "FAM"));

    // Global point numbering (NUM) -- optional; Salome/Code_Aster/Kratos
    // write it, this reader has ignored it entirely until now.
    if (h5::exists(noe, "NUM"))
        mesh.AddPointData("med:num", h5::read_dataset(noe, "NUM"));

    // Families info
    h5::Hid fas = h5::exists(data_grp, "FAS") ? h5::open_group(data_grp, "FAS") : h5::Hid();
    if (!fas.Valid()) {
        h5::Hid fas_root = h5::open_group(f, "FAS");
        fas = h5::open_group(fas_root, mesh_name);
    }
    if (h5::exists(fas, "NOEUD")) {
        h5::Hid noeud = h5::open_group(fas, "NOEUD");
        read_families(noeud, rInfo.mPointTags, rInfo.mPointTagGroups);
    }

    // Cells
    std::vector<std::string> cell_types;  // meshio names, in read order
    h5::Hid mai = h5::open_group(data_grp, "MAI");
    std::vector<NDArray> cell_tag_blocks;
    bool any_cell_tags = false;
    std::vector<NDArray> cell_num_blocks;
    std::size_t num_blocks_with_num = 0;
    // Cell-block order is significant (aligns cell_data / cell_sets); iterate in
    // HDF5 creation order to match the Python (h5py track_order) reader.
    for (const std::string& med_type : h5::group_links_crt(mai)) {
        auto it = med_to_meshio().find(med_type);
        if (it == med_to_meshio().end())
            throw ReadError(detail::format_compat("MED: unsupported cell type {}", med_type));
        h5::Hid g = h5::open_group(mai, med_type);

        if (med_type == "POE") {
            // MED_POLYHEDRON: NOD (flat nodes) + INN (face -> NOD) + IND
            // (cell -> face), all 1-based. Regrouped back into polyhedron<N>
            // by unique node count on the way out, the same bucketing the
            // OpenFOAM, EnSight and CGNS readers use.
            NDArray nod = h5::read_dataset(g, "NOD");
            NDArray inn = h5::read_dataset(g, "INN");
            NDArray ind = h5::read_dataset(g, "IND");
            const std::size_t ncells = ind.Size() > 0 ? ind.Size() - 1 : 0;
            std::vector<std::vector<std::vector<std::int64_t>>> cells(ncells);
            std::vector<std::size_t> node_counts(ncells, 0);
            for (std::size_t c = 0; c < ncells; ++c) {
                const std::int64_t f0 = detail::read_int(ind, c) - 1;
                const std::int64_t f1 = detail::read_int(ind, c + 1) - 1;
                std::vector<std::int64_t> uniq;
                for (std::int64_t f = f0; f < f1; ++f) {
                    const std::int64_t a = detail::read_int(inn, static_cast<std::size_t>(f)) - 1;
                    const std::int64_t b =
                        detail::read_int(inn, static_cast<std::size_t>(f) + 1) - 1;
                    std::vector<std::int64_t> face;
                    for (std::int64_t j = a; j < b; ++j)
                        face.push_back(detail::read_int(nod, static_cast<std::size_t>(j)) - 1);
                    uniq.insert(uniq.end(), face.begin(), face.end());
                    cells[c].push_back(std::move(face));
                }
                std::sort(uniq.begin(), uniq.end());
                uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
                node_counts[c] = uniq.size();
            }
            std::vector<std::size_t> order;
            std::map<std::size_t, std::vector<std::size_t>> groups;
            for (std::size_t c = 0; c < ncells; ++c) {
                if (groups.find(node_counts[c]) == groups.end())
                    order.push_back(node_counts[c]);
                groups[node_counts[c]].push_back(c);
            }
            for (std::size_t n : order) {
                std::vector<std::vector<std::vector<std::int64_t>>> group;
                group.reserve(groups[n].size());
                for (std::size_t c : groups[n])
                    group.push_back(std::move(cells[c]));
                const std::string tname = "polyhedron" + std::to_string(n);
                mesh.AddPolyhedronBlock(tname, std::move(group));
                cell_types.push_back(tname);
            }
        } else if (med_type == "POG" || med_type == "POG2") {
            // Ragged polygons: flat 1-based NOD + 1-based INN offsets.
            NDArray nod = h5::read_dataset(g, "NOD");
            NDArray inn = h5::read_dataset(g, "INN");
            std::size_t npoly = inn.Size() > 0 ? inn.Size() - 1 : 0;
            std::vector<std::vector<std::int64_t>> rows;
            for (std::size_t i = 0; i < npoly; ++i) {
                std::int64_t a = detail::read_int(inn, i) - 1;
                std::int64_t b = detail::read_int(inn, i + 1) - 1;
                std::vector<std::int64_t> row;
                for (std::int64_t j = a; j < b; ++j)
                    row.push_back(detail::read_int(nod, static_cast<std::size_t>(j)) - 1);
                rows.push_back(std::move(row));
            }
            mesh.AddPolygonBlock(it->second, std::move(rows));
            cell_types.push_back(it->second);
        } else {
            h5::Hid nod_ds(H5Dopen2(g, "NOD", H5P_DEFAULT), H5Dclose);
            if (!nod_ds.Valid())
                throw ReadError(detail::format_compat("MED: missing NOD for {}", med_type));
            std::int64_t n_cells = h5::read_attr_int(nod_ds, "NBR");
            NDArray nod = h5::read_dataset(g, "NOD");
            std::size_t k = n_cells > 0 ? nod.Size() / static_cast<std::size_t>(n_cells) : 0;
            warn_unconverted_3d(it->second);
            // Fuse the Fortran->C transpose (shift -1) with the MED->meshio
            // node reorder into a single pass over the connectivity.
            auto pit = med_node_perm().find(it->second);
            const std::vector<int>* perm =
                (pit != med_node_perm().end() && pit->second.size() == k) ? &pit->second : nullptr;
            NDArray data = unflatten_f(nod, static_cast<std::size_t>(n_cells), k, -1, perm);
            mesh.AddCellBlock(it->second, std::move(data));
            cell_types.push_back(it->second);
        }

        // One entry per block, always -- a block with no `FAM` gets a
        // placeholder filled in with family 0 below, so the array stays
        // aligned with `mesh.cells`.
        if (h5::exists(g, "FAM")) {
            cell_tag_blocks.push_back(h5::read_dataset(g, "FAM"));
            any_cell_tags = true;
        } else {
            cell_tag_blocks.emplace_back();
        }

        // Global cell numbering (NUM) -- optional, and only carried when
        // *every* block has it: a partial NUM array cannot be a mesh-wide
        // "global" numbering, and fabricating the missing entries (as the
        // Kratos MedApplication does with iota) would be a wrong answer, not
        // an honest gap.
        if (h5::exists(g, "NUM")) {
            cell_num_blocks.push_back(h5::read_dataset(g, "NUM"));
            ++num_blocks_with_num;
        } else {
            cell_num_blocks.emplace_back();
        }
    }
    if (any_cell_tags) {
        // A block the file left `FAM` off of belongs to no family, and MED
        // spells "no family" as id **0** -- so filling those blocks with zeros
        // is the file's own meaning, not a guess. This used to throw ("partial
        // cell tags handled by Python fallback"), which made a perfectly
        // ordinary Salome file unreadable wherever there is no Python; the
        // Python reference meanwhile appended only the blocks that had a
        // `FAM`, leaving `cell_data["cell_tags"]` *shorter* than `mesh.cells`
        // and so violating the one-array-per-block invariant. Both now do
        // this, so the two agree and both are well-formed.
        bool any_filled = false;
        for (std::size_t b = 0; b < cell_tag_blocks.size(); ++b) {
            if (cell_tag_blocks[b].Size() != 0)
                continue;
            const std::size_t n = mesh.Cells(b).NumCells();
            NDArray zeros(DType::Int32, {n});
            std::memset(zeros.Data(), 0, zeros.Nbytes());
            cell_tag_blocks[b] = std::move(zeros);
            any_filled = true;
        }
        if (any_filled)
            log::warn(
                "MED: some cell blocks carry no FAM dataset; those cells are reported as family 0 "
                "(MED's \"no family\").");
        mesh.AddCellData("cell_tags", std::move(cell_tag_blocks));
    }
    if (num_blocks_with_num > 0) {
        if (num_blocks_with_num == mesh.NumCellBlocks()) {
            mesh.AddCellData("med:num", std::move(cell_num_blocks));
        } else {
            log::warn(
                "MED: cell NUM is present on only {} of {} cell blocks; ignoring "
                "'med:num' for this mesh.",
                num_blocks_with_num, mesh.NumCellBlocks());
        }
    }

    if (h5::exists(fas, "ELEME")) {
        h5::Hid eleme = h5::open_group(fas, "ELEME");
        read_families(eleme, rInfo.mCellTags, rInfo.mCellTagGroups);
    }

    // Named regions derived from the family tables just read, one per group
    // name (see doc/regions.md). Kept independent of point_tags/cell_tags:
    // both representations are populated and neither is derived from the
    // other on this path.
    med_attach_point_regions(mesh, rInfo);
    med_attach_cell_regions(mesh, rInfo);

    // Fields (CHA): the single-timestep, default-profile common case is read
    // directly (see read_cha_fields). Anything past that scope -- units,
    // multi-timestep metadata, a named profile, ELNO/ELGA support -- throws
    // from inside it under the default (strict) options and defers the whole
    // file to Python; `ReadOptions::mLenient`/`mTimeStep` are what let a
    // caller with no Python fallback read it anyway.
    if (h5::exists(f, "CHA")) {
        h5::Hid cha = h5::open_group(f, "CHA");
        read_cha_fields(cha, mesh, rOptions.mLenient, rOptions, &rInfo);
    }

    return mesh;
}

}  // namespace

Mesh read_med(const std::string& rPath, MedInfo& rInfo) {
    return med_read_impl(rPath, rInfo, ReadOptions{});
}

Mesh read_med(const std::string& rPath, MedInfo& rInfo, const ReadOptions& rOptions) {
    return med_read_impl(rPath, rInfo, rOptions);
}

void write_med(const std::string& rPath, const Mesh& rMesh, const MedInfo& rInfo,
               const std::string& rMedVersion) {
    h5::SilenceErrors silence;

    // Fields (CHA): the single-timestep, no-profile, no-units common case is
    // written directly (see write_cha_nodal_field/write_cha_cell_field
    // below). The enhanced Python writer's multi-timestep metadata, units and
    // component names (med:step_meta/med:field_units/med:nom) are Python-only
    // conventions -- they are ordinary dicts/lists of strings, which cannot
    // become an NDArray, so they can never reach this Mesh through ANY
    // binding's field_data conversion (Python's own `med_write` pybind
    // wrapper uses `lenient_field_data=true` specifically because of this:
    // a non-numeric field_data entry is silently dropped before it gets
    // here, never thrown). The guard therefore has to live where those
    // conventions actually exist -- the Python shim (`med/__init__.py`),
    // which checks `mesh.field_data` itself before ever calling in here.
    //
    // `gmsh:physical` used to be an unconditional throw here ("handled by
    // Python fallback") -- but the information it carries is already
    // representable natively: med_cell_regions_to_tags (below) now folds it,
    // alongside any named `Cell` region, into the same family synthesis the
    // v9.6.0 regions work added. See `"gmsh:physical"` in the two CHA
    // skip-lists a few lines down -- it must never *also* be written as a
    // numeric field.

    // Named regions -> families (see doc/regions.md): only synthesized when
    // the mesh carries no native point_tags/cell_tags of its own -- a mesh
    // read from MED (or built with them directly) writes exactly as before.
    NDArray synth_point_fam;
    std::map<std::int64_t, std::vector<std::string>> synth_point_tags;
    std::map<std::int64_t, std::string> synth_point_group_names;
    const bool synthesized_point =
        !rMesh.HasPointData("point_tags") &&
        med_point_regions_to_tags(rMesh, synth_point_fam, synth_point_tags,
                                  synth_point_group_names);

    std::vector<NDArray> synth_cell_fam_blocks;
    std::map<std::int64_t, std::vector<std::string>> synth_cell_tags;
    std::map<std::int64_t, std::string> synth_cell_group_names;
    const bool synthesized_cell = !rMesh.HasCellData("cell_tags") &&
                                  med_cell_regions_to_tags(rMesh, synth_cell_fam_blocks,
                                                           synth_cell_tags, synth_cell_group_names);

    if (synthesized_point || synthesized_cell)
        med_warn_side_regions_dropped(rMesh);

    const std::map<std::int64_t, std::vector<std::string>>& point_tags =
        synthesized_point ? synth_point_tags : rInfo.mPointTags;
    const std::map<std::int64_t, std::string>& point_tag_groups =
        synthesized_point ? synth_point_group_names : rInfo.mPointTagGroups;
    const std::map<std::int64_t, std::vector<std::string>>& cell_tags =
        synthesized_cell ? synth_cell_tags : rInfo.mCellTags;
    const std::map<std::int64_t, std::string>& cell_tag_groups =
        synthesized_cell ? synth_cell_group_names : rInfo.mCellTagGroups;

    // Parse med_version -> MAJ.MIN.REL (default 4.1.0 on error).
    int maj = 4, min = 1, rel = 0;
    {
        int parts[3] = {4, 1, 0};
        std::size_t start = 0, idx = 0;
        bool ok = true;
        for (idx = 0; idx < 3; ++idx) {
            std::size_t dot = rMedVersion.find('.', start);
            std::string tok = rMedVersion.substr(
                start, dot == std::string::npos ? std::string::npos : dot - start);
            try {
                parts[idx] = std::stoi(tok);
            } catch (...) {
                ok = false;
                break;
            }
            if (dot == std::string::npos)
                break;
            start = dot + 1;
        }
        if (ok) {
            maj = parts[0];
            min = parts[1];
            rel = parts[2];
        }
    }

    h5::Hid f = h5::create_file(rPath);

    h5::Hid infos = h5::create_group(f, "INFOS_GENERALES");
    h5::write_attr_int(infos, "MAJ", maj);
    h5::write_attr_int(infos, "MIN", min);
    h5::write_attr_int(infos, "REL", rel);

    const std::string mesh_name = rInfo.mMeshName.empty() ? "mesh" : rInfo.mMeshName;
    const std::size_t dim = rMesh.PointDim();

    h5::Hid ens = h5::create_group(f, "ENS_MAA");
    h5::Hid med_mesh = h5::create_group(ens, mesh_name);
    h5::write_attr_int(med_mesh, "DIM", static_cast<std::int64_t>(dim));
    h5::write_attr_int(med_mesh, "ESP", static_cast<std::int64_t>(dim));
    h5::write_attr_int(med_mesh, "REP", 0);
    write_attr_bytes(med_mesh, "UNT", rInfo.mUnitTime);
    write_attr_bytes(med_mesh, "UNI", rInfo.mUnitCoords);
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
    write_attr_bytes(
        med_mesh, "DES",
        rInfo.mDescription.empty() ? "Mesh created with meshio++" : rInfo.mDescription);
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
        NDArray coo = flatten_f(rMesh.Points(), 0);
        h5::write_dataset(noe, "COO", coo);
        h5::Hid d(H5Dopen2(noe, "COO", H5P_DEFAULT), H5Dclose);
        h5::write_attr_int(d, "CGT", 1);
        h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(rMesh.NumPoints()));
    }
    if (rMesh.HasPointData("point_tags") || synthesized_point) {
        const NDArray& point_fam =
            rMesh.HasPointData("point_tags") ? rMesh.PointData("point_tags") : synth_point_fam;
        h5::write_dataset(noe, "FAM", point_fam);
        h5::Hid d(H5Dopen2(noe, "FAM", H5P_DEFAULT), H5Dclose);
        h5::write_attr_int(d, "CGT", 1);
        h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(rMesh.NumPoints()));
    }
    if (rMesh.HasPointData("med:num")) {
        h5::write_dataset(noe, "NUM", rMesh.PointData("med:num"));
        h5::Hid d(H5Dopen2(noe, "NUM", H5P_DEFAULT), H5Dclose);
        h5::write_attr_int(d, "CGT", 1);
        h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(rMesh.NumPoints()));
    }

    // Cells. MED cannot have two sections of the same type -- MSH 4.1's
    // canonical structure is one cell block per *entity*, so group blocks by
    // type (first-seen order) and write one section per type, consolidating
    // the contributing blocks' connectivity/FAM/NUM. Mirrors the Python
    // reference's write-time merge (_med.py:811-869), which is why this
    // never rejects a mesh the pre-v9.8.0 pairwise-type check used to.
    h5::Hid mai = h5::create_group(time_step, "MAI");
    h5::write_attr_int(mai, "CGT", 1);
    const bool has_cell_num = rMesh.HasCellData("med:num");
    const bool has_native_cell_tags = rMesh.HasCellData("cell_tags");
    const std::size_t native_tag_blocks =
        has_native_cell_tags ? rMesh.CellDataNumBlocks("cell_tags") : 0;
    const std::size_t num_blocks = has_cell_num ? rMesh.CellDataNumBlocks("med:num") : 0;

    std::vector<std::string> cell_type_order;
    std::unordered_map<std::string, std::vector<std::size_t>> blocks_by_type;
    for (std::size_t k = 0; k < rMesh.NumCellBlocks(); ++k) {
        std::string ctype(rMesh.Cells(k).Type());
        // Every `polyhedron<N>` collapses to one key. MED holds ONE section per
        // type inside a MAI group, so grouping on the exact type string would
        // try to create POE three times over for a mesh carrying polyhedron4,
        // polyhedron5 and polyhedron8 -- an invalid file, and a group-creation
        // failure rather than a clean error. The node count is a meshio++
        // bucketing convention, not part of the MED type.
        if (ctype.rfind("polyhedron", 0) == 0)
            ctype = "polyhedron";
        if (meshio_to_med().find(ctype) == meshio_to_med().end())
            throw WriteError(detail::format_compat("MED: unsupported cell type {}", ctype));
        if (!blocks_by_type.count(ctype))
            cell_type_order.push_back(ctype);
        blocks_by_type[ctype].push_back(k);
    }

    for (const std::string& ctype : cell_type_order) {
        const std::vector<std::size_t>& idxs = blocks_by_type[ctype];
        const bool is_ragged_type = (ctype == "polygon" || ctype == "polygon2");
        const bool is_polyhedron_type = (ctype == "polyhedron");

        // Blocks of the same type must agree on node count to be merged into
        // one section -- unlike merge.cpp's grouping (keyed on Type() alone,
        // trusting the first contributor), a disagreement here is a checked
        // error rather than a silently truncated/garbled NOD array.
        if (!is_ragged_type && !is_polyhedron_type) {
            const std::size_t npc = rMesh.Cells(idxs[0]).NodesPerCell();
            for (std::size_t bi : idxs)
                if (rMesh.Cells(bi).NodesPerCell() != npc)
                    throw WriteError(detail::format_compat(
                        "MED: cell type '{}' has blocks disagreeing on node count ({} vs {}); "
                        "cannot consolidate into one section",
                        ctype, npc, rMesh.Cells(bi).NodesPerCell()));
        }

        h5::Hid g = h5::create_group(mai, meshio_to_med().at(ctype));
        h5::write_attr_int(g, "CGT", 1);
        h5::write_attr_int(g, "CGS", 1);
        write_attr_bytes(g, "PFL", kProfile);

        std::size_t n_total = 0;
        for (std::size_t bi : idxs)
            n_total += rMesh.Cells(bi).NumCells();

        if (is_polyhedron_type) {
            // MED_POLYHEDRON needs THREE 1-based arrays where a polygon needs
            // two: NOD (every face's nodes, flat), INN (face -> start in NOD),
            // and IND (cell -> start in the face list). Concatenated across
            // every contributing block, in block order.
            std::vector<std::int64_t> nod;
            std::vector<std::int64_t> inn = {1};
            std::vector<std::int64_t> ind = {1};
            for (std::size_t bi : idxs) {
                const auto cb = rMesh.Cells(bi);
                for (std::size_t i = 0; i < cb.NumCells(); ++i) {
                    for (std::size_t f = 0; f < cb.NumFaces(i); ++f) {
                        const auto face = cb.Face(i, f);
                        for (std::size_t j = 0; j < face.second; ++j)
                            nod.push_back(face.first[j] + 1);
                        inn.push_back(inn.back() + static_cast<std::int64_t>(face.second));
                    }
                    ind.push_back(static_cast<std::int64_t>(inn.size()));
                }
            }
            auto to_arr = [](const std::vector<std::int64_t>& v) {
                NDArray a(DType::Int64, {v.size()});
                for (std::size_t i = 0; i < v.size(); ++i)
                    a.As<std::int64_t>()[i] = v[i];
                return a;
            };
            h5::write_dataset(g, "NOD", to_arr(nod));
            h5::write_dataset(g, "INN", to_arr(inn));
            h5::write_dataset(g, "IND", to_arr(ind));
            h5::Hid d(H5Dopen2(g, "NOD", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "CGT", 1);
            h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(n_total));
        } else if (is_ragged_type) {
            // Ragged: flat 1-based NOD + 1-based INN offsets, concatenated
            // across every contributing block in order.
            std::vector<std::int64_t> nod;
            std::vector<std::int64_t> inn = {1};
            for (std::size_t bi : idxs) {
                const auto cb = rMesh.Cells(bi);
                for (std::size_t i = 0; i < cb.NumCells(); ++i) {
                    const std::int64_t* row = cb.Row(i);
                    const std::size_t row_size = cb.RowSize(i);
                    for (std::size_t j = 0; j < row_size; ++j)
                        nod.push_back(row[j] + 1);
                    inn.push_back(inn.back() + static_cast<std::int64_t>(row_size));
                }
            }
            NDArray nod_a(DType::Int64, {nod.size()});
            for (std::size_t i = 0; i < nod.size(); ++i)
                nod_a.As<std::int64_t>()[i] = nod[i];
            NDArray inn_a(DType::Int64, {inn.size()});
            for (std::size_t i = 0; i < inn.size(); ++i)
                inn_a.As<std::int64_t>()[i] = inn[i];
            h5::write_dataset(g, "NOD", nod_a);
            h5::write_dataset(g, "INN", inn_a);
            h5::Hid d(H5Dopen2(g, "NOD", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "CGT", 1);
            h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(n_total));
        } else {
            warn_unconverted_3d(ctype);
            // Fuse the meshio->MED node reorder with the Fortran transpose
            // (shift +1) into a single pass (mirrors the read side), over the
            // concatenated connectivity of every contributing block.
            auto pit = med_node_perm().find(ctype);
            const std::vector<int>* perm = (pit != med_node_perm().end()) ? &pit->second : nullptr;
            NDArray conn = med_concat_conn_rows(rMesh, idxs, rMesh.Cells(idxs[0]).NodesPerCell());
            NDArray nod = flatten_f(conn, +1, perm);
            h5::write_dataset(g, "NOD", nod);
            h5::Hid d(H5Dopen2(g, "NOD", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "CGT", 1);
            h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(n_total));
        }

        // FAM: native cell_tags wins when it covers every contributing block
        // (the pre-existing per-block priority); else the synthesized
        // regions/gmsh:physical families, which always cover every block by
        // construction (synth_cell_fam_blocks has one entry per mesh block).
        bool native_covers_all = has_native_cell_tags;
        for (std::size_t bi : idxs)
            if (bi >= native_tag_blocks) {
                native_covers_all = false;
                break;
            }
        if (native_covers_all) {
            NDArray fam = med_concat_cell_data_rows(rMesh, "cell_tags", idxs);
            h5::write_dataset(g, "FAM", fam);
            h5::Hid d(H5Dopen2(g, "FAM", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "CGT", 1);
            h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(n_total));
        } else if (has_native_cell_tags) {
            log::warn(
                "MED: cell type '{}' has 'cell_tags' covering only some of its blocks; FAM not "
                "written for this section.",
                ctype);
            detail::provenance_note("data-dropped",
                                    "FAM not written for cell type '" + std::string(ctype) +
                                        "' -- cell_tags covers only some of its blocks");
        } else if (synthesized_cell) {
            NDArray fam = med_concat_ndarray_rows(synth_cell_fam_blocks, idxs);
            h5::write_dataset(g, "FAM", fam);
            h5::Hid d(H5Dopen2(g, "FAM", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "CGT", 1);
            h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(n_total));
        }

        // NUM: same all-or-nothing rule as global numbering elsewhere.
        bool num_covers_all = has_cell_num;
        for (std::size_t bi : idxs)
            if (bi >= num_blocks) {
                num_covers_all = false;
                break;
            }
        if (num_covers_all) {
            NDArray num = med_concat_cell_data_rows(rMesh, "med:num", idxs);
            h5::write_dataset(g, "NUM", num);
            h5::Hid d(H5Dopen2(g, "NUM", H5P_DEFAULT), H5Dclose);
            h5::write_attr_int(d, "CGT", 1);
            h5::write_attr_int(d, "NBR", static_cast<std::int64_t>(n_total));
        } else if (has_cell_num) {
            log::warn(
                "MED: cell type '{}' has 'med:num' covering only some of its blocks; NUM not "
                "written for this section.",
                ctype);
            detail::provenance_note("data-dropped",
                                    "NUM not written for cell type '" + std::string(ctype) +
                                        "' -- med:num covers only some of its blocks");
        }
    }

    // Families
    h5::Hid fas = h5::create_group(f, "FAS");
    h5::Hid families = h5::create_group(fas, mesh_name);
    h5::Hid family_zero = h5::create_group(families, "FAMILLE_ZERO");
    h5::write_attr_int(family_zero, "NUM", 0);
    if (!point_tags.empty()) {
        h5::Hid node = h5::create_group(families, "NOEUD");
        write_families(node, point_tags, point_tag_groups);
    }
    if (!cell_tags.empty()) {
        h5::Hid element = h5::create_group(families, "ELEME");
        write_families(element, cell_tags, cell_tag_groups);
    }

    // Fields (CHA) -- single-timestep common case only; see the guard above
    // and write_cha_nodal_field/write_cha_cell_field's own doc comments.
    bool has_point_fields = false;
    for (const auto& name : rMesh.PointDataNames())
        if (name != "point_tags" && name != "med:num") {
            has_point_fields = true;
            break;
        }
    bool has_cell_fields = false;
    for (const auto& name : rMesh.CellDataNames())
        if (name != "cell_tags" && name != "med:num") {
            has_cell_fields = true;
            break;
        }
    if (has_point_fields || has_cell_fields) {
        // A field's row count IS its entity count -- `write_cha_support`
        // writes `NBR = rows(data)` and the reader validates that against
        // `NumPoints()`/`NumCells()`, so a mis-shaped array (e.g. an (n,3)
        // vector flattened to (3n,) by a caller that lost its component
        // count) produces a file this very reader rejects. Catch it here,
        // where the array and both counts can be named, rather than emitting
        // an unreadable file and failing on the way back in.
        for (const auto& name : rMesh.PointDataNames()) {
            if (name == "point_tags" || name == "med:num")
                continue;
            const std::size_t rows = detail::rows(rMesh.PointData(name));
            if (rows != rMesh.NumPoints())
                throw WriteError(detail::format_compat(
                    "MED: point_data '{}' has {} rows but the mesh has {} points; a "
                    "multi-component field must be shaped (n_points, n_components), not "
                    "flattened",
                    name, rows, rMesh.NumPoints()));
        }
        for (const auto& name : rMesh.CellDataNames()) {
            if (name == "cell_tags" || name == "med:num")
                continue;
            for (std::size_t b = 0; b < rMesh.CellDataNumBlocks(name) && b < rMesh.NumCellBlocks();
                 ++b) {
                const std::size_t rows = detail::rows(rMesh.CellData(name, b));
                if (rows != rMesh.Cells(b).NumCells())
                    throw WriteError(detail::format_compat(
                        "MED: cell_data '{}' block {} has {} rows but that block has {} cells; a "
                        "multi-component field must be shaped (n_cells, n_components), not "
                        "flattened",
                        name, b, rows, rMesh.Cells(b).NumCells()));
            }
        }

        h5::Hid cha = h5::create_group(f, "CHA");
        for (const auto& name : rMesh.PointDataNames())
            if (name != "point_tags" && name != "med:num")
                write_cha_nodal_field(cha, mesh_name, name, rMesh.PointData(name));
        for (const auto& name : rMesh.CellDataNames())
            if (name != "cell_tags" && name != "med:num")
                write_cha_cell_field(cha, mesh_name, name, rMesh);
    }
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
