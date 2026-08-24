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
#ifdef MESHIOPLUSPLUS_HAS_NETCDF

// External includes
#include <netcdf.h>

// System includes
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/formats/exodus.hpp"
#include "meshioplusplus/detail/cell_index.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/types.hpp"

namespace meshioplusplus {

namespace {

void check(int status, const char* pWhat, bool writing = false) {
    if (status != NC_NOERR) {
        std::string msg = std::string("Exodus/netCDF: ") + pWhat + ": " + nc_strerror(status);
        if (writing)
            throw WriteError(msg);
        throw ReadError(msg);
    }
}

const std::unordered_map<std::string, std::string>& exodus_to_meshio() {
    static const std::unordered_map<std::string, std::string> m = {
        {"SPHERE", "vertex"},      {"BEAM", "line"},        {"BEAM2", "line"},
        {"BEAM3", "line3"},        {"BAR2", "line"},        {"SHELL", "quad"},
        {"SHELL4", "quad"},        {"SHELL8", "quad8"},     {"SHELL9", "quad9"},
        {"QUAD", "quad"},          {"QUAD4", "quad"},       {"QUAD5", "quad5"},
        {"QUAD8", "quad8"},        {"QUAD9", "quad9"},      {"TRI", "triangle"},
        {"TRIANGLE", "triangle"},  {"TRI3", "triangle"},    {"TRI6", "triangle6"},
        {"TRI7", "triangle7"},     {"HEX", "hexahedron"},   {"HEXAHEDRON", "hexahedron"},
        {"HEX8", "hexahedron"},    {"HEX9", "hexahedron9"}, {"HEX20", "hexahedron20"},
        {"HEX27", "hexahedron27"}, {"TETRA", "tetra"},      {"TETRA4", "tetra4"},
        {"TET4", "tetra4"},        {"TETRA8", "tetra8"},    {"TETRA10", "tetra10"},
        {"TETRA14", "tetra14"},    {"PYRAMID", "pyramid"},  {"WEDGE", "wedge"}};
    return m;
}

// The Python reverse map is last-wins over dict order.
const std::unordered_map<std::string, std::string>& meshio_to_exodus() {
    static const std::unordered_map<std::string, std::string> m = {
        {"vertex", "SPHERE"},      {"line", "BAR2"},          {"line3", "BEAM3"},
        {"quad", "QUAD4"},         {"quad5", "QUAD5"},        {"quad8", "QUAD8"},
        {"quad9", "QUAD9"},        {"triangle", "TRI3"},      {"triangle6", "TRI6"},
        {"triangle7", "TRI7"},     {"hexahedron", "HEX8"},    {"hexahedron9", "HEX9"},
        {"hexahedron20", "HEX20"}, {"hexahedron27", "HEX27"}, {"tetra", "TETRA"},
        {"tetra4", "TET4"},        {"tetra8", "TETRA8"},      {"tetra10", "TETRA10"},
        {"tetra14", "TETRA14"},    {"pyramid", "PYRAMID"},    {"wedge", "WEDGE"}};
    return m;
}

nc_type nc_type_of(DType dt) {
    switch (dt) {
        case DType::Float32:
            return NC_FLOAT;
        case DType::Float64:
            return NC_DOUBLE;
        case DType::Int8:
            return NC_BYTE;
        case DType::Int16:
            return NC_SHORT;
        case DType::Int32:
            return NC_INT;
        case DType::Int64:
            return NC_INT64;
        case DType::UInt8:
            return NC_UBYTE;
        case DType::UInt16:
            return NC_USHORT;
        case DType::UInt32:
            return NC_UINT;
        case DType::UInt64:
            return NC_UINT64;
    }
    return NC_DOUBLE;
}

DType dtype_of(nc_type t) {
    switch (t) {
        case NC_FLOAT:
            return DType::Float32;
        case NC_DOUBLE:
            return DType::Float64;
        case NC_BYTE:
            return DType::Int8;
        case NC_SHORT:
            return DType::Int16;
        case NC_INT:
            return DType::Int32;
        case NC_INT64:
            return DType::Int64;
        case NC_UBYTE:
            return DType::UInt8;
        case NC_USHORT:
            return DType::UInt16;
        case NC_UINT:
            return DType::UInt32;
        case NC_UINT64:
            return DType::UInt64;
        default:
            throw ReadError("Exodus: unsupported netCDF variable type");
    }
}

// Read a whole variable (or a start/count hyperslab) into an NDArray.
NDArray read_var(int ncid, int varid, const std::vector<std::size_t>& rStart,
                 const std::vector<std::size_t>& rCount) {
    nc_type t;
    check(nc_inq_vartype(ncid, varid, &t), "inq_vartype");
    DType dt = dtype_of(t);
    std::vector<std::size_t> shape;
    for (std::size_t c : rCount)
        shape.push_back(c);
    NDArray out(dt, shape);
    if (out.Size() > 0)
        check(nc_get_vara(ncid, varid, rStart.data(), rCount.data(), out.Data()), "get_vara");
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

// A char variable whose LAST dimension is the string width -> list of strings.
//
// Covers both shapes Exodus uses: the 2-D `(n, len_string)` name arrays
// (`eb_names`, `ns_names`, `ss_names`, `name_nod_var`, ...) and the 3-D
// `(num_qa_rec, four, len_string)` QA records, which flatten to 4 strings per
// record -- the same order the Python reference reader appends them in.
std::vector<std::string> read_names(int ncid, int varid) {
    std::vector<std::size_t> dims = var_dims(ncid, varid);
    if (dims.empty())
        return {};
    const std::size_t w = dims.back();
    std::size_t n = 1;
    for (std::size_t k = 0; k + 1 < dims.size(); ++k)
        n *= dims[k];
    std::vector<char> buf(n * w, '\0');
    if (n * w > 0)
        check(nc_get_var_text(ncid, varid, buf.data()), "get names");
    std::vector<std::string> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; ++i)
        out.emplace_back(buf.data() + i * w, strnlen(buf.data() + i * w, w));
    return out;
}

// A netCDF text attribute, with the padding some writers count as part of its
// own length trimmed off.
//
// NetCDF.jl -- what PeriLab and other Julia solvers write Exodus with -- stores
// the C string's terminating NUL *inside* the attribute, so a SPHERE block's
// `elem_type` arrives as the 7 characters "SPHERE\0" and matches no key in
// `exodus_to_meshio()`. The read then failed with "unknown element type SPHERE"
// -- the NUL invisible, since `std::runtime_error::what()` is a `const char*`
// that stops at it. netCDF4-python strips the NUL on the way in, so the Python
// reference never saw this and the shim's fallback hid it everywhere except
// WASM, which has no fallback. Fixed-width writers pad with spaces instead;
// both are trimmed here.
std::string read_att_text(int ncid, int varid, const char* pName) {
    std::size_t attlen = 0;
    check(nc_inq_attlen(ncid, varid, pName, &attlen), pName);
    std::string out(attlen, '\0');
    if (attlen > 0)
        check(nc_get_att_text(ncid, varid, pName, out.data()), pName);
    // The set is built with an explicit length: a `const char*` would end at
    // the NUL that is the whole point of this trim.
    const std::size_t last = out.find_last_not_of(std::string("\0 \t\r\n", 5));
    out.erase(last == std::string::npos ? 0 : last + 1);
    return out;
}

// Read a 1-D integer variable as int64 (ids, set members, side numbers).
std::vector<std::int64_t> read_ints(int ncid, int varid) {
    std::vector<std::size_t> dims = var_dims(ncid, varid);
    const std::size_t n = dims.empty() ? 0 : dims[0];
    std::vector<std::int64_t> out(n, 0);
    if (n > 0)
        check(nc_get_var_longlong(ncid, varid, reinterpret_cast<long long*>(out.data())),
              "get ints");
    return out;
}

// The numeric suffix of a variable name, e.g. "connect12" -> 12, "connect" -> 1.
int exo_suffix(const std::string& rKey, std::size_t Prefix) {
    return rKey.size() > Prefix ? std::atoi(rKey.c_str() + Prefix) : 1;
}

// The file's `time_whole` values, or empty when it records none.
//
// Read in its own pass rather than inside the main variable loop: the step must
// be resolved (and an out-of-range request rejected) *before* any data array is
// sliced, and netCDF gives no ordering guarantee that would put `time_whole`
// ahead of `vals_nod_var1`.
// The start offset for a data variable's time dimension.
//
// `time_whole` and a `vals_*` variable can legitimately disagree on length (a
// writer that crashed mid-step leaves a short data array), so the step resolved
// against `time_whole` is re-checked against the variable actually being
// sliced. Naming the variable matters: "step 4 of 3" is far easier to act on
// when it says which array came up short.
std::size_t exo_step_offset(std::size_t Step, std::size_t Available, const std::string& rKey) {
    if (Step >= Available)
        throw ReadError("Exodus: time step " + std::to_string(Step) + " is out of range for '" +
                        rKey + "', which has " + std::to_string(Available) +
                        (Available == 1 ? " step" : " steps"));
    return Step;
}

std::vector<double> exo_time_values(int ncid) {
    int varid;
    if (nc_inq_varid(ncid, "time_whole", &varid) != NC_NOERR)
        return {};
    std::vector<std::size_t> dims = var_dims(ncid, varid);
    const std::size_t n = dims.empty() ? 0 : dims[0];
    std::vector<double> out(n, 0.0);
    if (n > 0)
        check(nc_get_var_double(ncid, varid, out.data()), "time_whole");
    return out;
}

}  // namespace

int exo_face_index(const std::string& rCellType, int ExodusSide) {
    if (ExodusSide < 1)
        return -1;
    const int s = ExodusSide - 1;  // 0-based Exodus side number

    // Exodus side->node lists are from the Exodus II spec's side-set numbering
    // tables; the meshio++ rows are detail/cell_faces.cpp. Each entry is the
    // meshio++ facet index whose corner-node SET equals that Exodus side's.
    //
    // tetra: Exodus S1={1,2,4} S2={2,3,4} S3={1,4,3} S4={1,3,2}
    //        meshio++ 0={0,1,3} 1={1,2,3} 2={2,0,3} 3={0,2,1}
    static const int tetra[4] = {0, 1, 2, 3};
    // hexahedron: Exodus S1={1,2,6,5} S2={2,3,7,6} S3={3,4,8,7}
    //                    S4={4,1,5,8} S5={1,4,3,2} S6={5,6,7,8}
    //             meshio++ 0={0,4,7,3} 1={1,2,6,5} 2={0,1,5,4}
    //                      3={3,7,6,2} 4={0,3,2,1} 5={4,5,6,7}
    static const int hexa[6] = {2, 1, 3, 0, 4, 5};
    // wedge: Exodus S1={1,2,5,4} S2={2,3,6,5} S3={1,4,6,3} S4={1,3,2} S5={4,5,6}
    //        meshio++ 0={0,2,1} 1={3,4,5} 2={0,1,4,3} 3={1,2,5,4} 4={2,0,3,5}
    static const int wedge[5] = {2, 3, 4, 0, 1};
    // pyramid: Exodus S1={1,2,5} S2={2,3,5} S3={3,4,5} S4={4,1,5} S5={1,4,3,2}
    //          meshio++ 0={0,3,2,1} 1={0,1,4} 2={1,2,4} 3={2,3,4} 4={3,0,4}
    static const int pyramid[5] = {1, 2, 3, 4, 0};
    // 2-D elements: Exodus numbers the edges 1-2, 2-3, ... in the same order
    // detail/cell_edges.hpp walks them, so the mapping is the identity.
    static const int tri[3] = {0, 1, 2};
    static const int quad[4] = {0, 1, 2, 3};

    const int* table = nullptr;
    int count = 0;
    // Higher-order variants share their linear base's facet ordering.
    if (rCellType == "tetra" || rCellType == "tetra4" || rCellType == "tetra10" ||
        rCellType == "tetra14") {
        table = tetra;
        count = 4;
    } else if (rCellType == "hexahedron" || rCellType == "hexahedron20" ||
               rCellType == "hexahedron27") {
        table = hexa;
        count = 6;
    } else if (rCellType == "wedge" || rCellType == "wedge15" || rCellType == "wedge18") {
        table = wedge;
        count = 5;
    } else if (rCellType == "pyramid" || rCellType == "pyramid13" || rCellType == "pyramid14") {
        table = pyramid;
        count = 5;
    } else if (rCellType == "triangle" || rCellType == "triangle6" || rCellType == "triangle7") {
        table = tri;
        count = 3;
    } else if (rCellType == "quad" || rCellType == "quad8" || rCellType == "quad9") {
        table = quad;
        count = 4;
    }
    if (!table || s >= count)
        return -1;
    return table[s];
}

namespace {

// categorize() from _exodus.py: recombine <name>X/Y/Z triplets and
// <name>_R/_Z doubles.
struct Categorized {
    std::vector<std::pair<std::string, int>> mSingle;
    std::vector<std::array<int, 2>> mDoubleIdx;
    std::vector<std::string> mDoubleName;
    std::vector<std::array<int, 3>> mTripleIdx;
    std::vector<std::string> mTripleName;
};

int index_of(const std::vector<std::string>& rNames, const std::string& s) {
    auto it = std::find(rNames.begin(), rNames.end(), s);
    return it == rNames.end() ? -1 : static_cast<int>(it - rNames.begin());
}

Categorized categorize(const std::vector<std::string>& rNames) {
    Categorized out;
    std::vector<bool> accounted(rNames.size(), false);
    for (std::size_t k = 0; k < rNames.size(); ++k) {
        if (accounted[k])
            continue;
        const std::string& name = rNames[k];
        if (!name.empty() && name.back() == 'X') {
            int ix = static_cast<int>(k);
            int iy = index_of(rNames, name.substr(0, name.size() - 1) + "Y");
            int iz = index_of(rNames, name.substr(0, name.size() - 1) + "Z");
            // NB: Python checks truthiness, so index 0 counts as "not found".
            if (iy > 0 && iz > 0) {
                out.mTripleIdx.push_back({ix, iy, iz});
                out.mTripleName.push_back(name.substr(0, name.size() - 1));
                accounted[ix] = accounted[iy] = accounted[iz] = true;
            } else {
                out.mSingle.emplace_back(name, ix);
                accounted[ix] = true;
            }
        } else if (name.size() >= 2 && name.compare(name.size() - 2, 2, "_R") == 0) {
            int ir = static_cast<int>(k);
            int iz = index_of(rNames, name.substr(0, name.size() - 2) + "_Z");
            if (iz > 0) {
                out.mDoubleIdx.push_back({ir, iz});
                out.mDoubleName.push_back(name.substr(0, name.size() - 2));
                accounted[ir] = accounted[iz] = true;
            } else {
                out.mSingle.emplace_back(name, ir);
                accounted[ir] = true;
            }
        } else {
            out.mSingle.emplace_back(name, static_cast<int>(k));
            accounted[k] = true;
        }
    }
    for (bool a : accounted)
        if (!a)
            throw ReadError("Exodus: inconsistent point data names");
    return out;
}

// Build an Int64 NDArray of shape (n,) or (n, 2) from a flat vector.
NDArray exo_entries(const std::vector<std::int64_t>& rValues, std::size_t Stride) {
    const std::size_t n = Stride == 0 ? 0 : rValues.size() / Stride;
    NDArray out = Stride == 1 ? NDArray::Uninit(DType::Int64, {n})
                              : NDArray::Uninit(DType::Int64, {n, Stride});
    for (std::size_t k = 0; k < rValues.size(); ++k)
        out.As<std::int64_t>()[k] = rValues[k];
    return out;
}

// The name for set/block index `k`, or a stable synthetic one.
//
// `eb_names`/`ns_names`/`ss_names` are optional, and SEACAS writes them blank
// for unnamed groups -- but a group with no name is still a group, and dropping
// it would lose the only handle a consumer has on it. `"<Prefix> <id>"` keys off
// the id rather than the index so it stays meaningful.
std::string exo_group_name(const std::vector<std::string>& rNames,
                           const std::vector<std::int64_t>& rIds, std::size_t Index,
                           const char* pPrefix) {
    if (Index < rNames.size() && !rNames[Index].empty())
        return rNames[Index];
    const std::int64_t id =
        Index < rIds.size() ? rIds[Index] : static_cast<std::int64_t>(Index) + 1;
    return std::string(pPrefix) + " " + std::to_string(id);
}

// Turn Exodus's element blocks, node sets and side sets into named regions.
//
// Runs after the cell blocks are in the mesh, because Cell/Side entries are
// GLOBAL block-major cell indices and those only exist once the blocks are
// ordered. The bases come from detail::block_bases -- the single owner of that
// numbering -- rather than being re-derived here.
/**
 * @brief The inverse of `exo_add_regions`' element-block half: one name per
 *        cell block, taken from the `Cell` region that covers exactly it.
 *
 * The read side gives every `connect{k}` a `Cell` region whose entries are the
 * contiguous global range `[bases[k], bases[k+1])`, named from `eb_names`. So a
 * region matching that range exactly is that block's name, and writing it back
 * to `eb_names` is what makes a block name survive a round trip instead of
 * being replaced by the synthetic `"Block N"` the reader falls back to. A block
 * with no such region gets an empty name, which is exactly what SEACAS itself
 * writes when it has none.
 */
std::vector<std::string> exo_block_names_from_regions(const Mesh& rMesh) {
    std::vector<std::string> names(rMesh.NumCellBlocks());
    if (rMesh.NumRegions() == 0)
        return names;
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
        const Region& r = rMesh.Region(i);
        if (r.mKind != RegionKind::Cell)
            continue;
        for (std::size_t k = 0; k + 1 < bases.size(); ++k) {
            const std::size_t n = static_cast<std::size_t>(bases[k + 1] - bases[k]);
            if (!names[k].empty() || r.NumEntries() != n || n == 0)
                continue;
            // Entries are canonical (sorted, de-duplicated), so "covers exactly
            // this block" is a first/last check rather than a set comparison.
            const std::int64_t* e = r.Entries();
            if (e[0] == bases[k] && e[n - 1] == bases[k + 1] - 1)
                names[k] = r.mName;
        }
    }
    return names;
}

void exo_add_regions(Mesh& rMesh, const std::vector<std::string>& rBlockTypes,
                     const std::vector<std::string>& rEbNames,
                     const std::vector<std::int64_t>& rEbIds,
                     const std::vector<std::string>& rNsNames,
                     const std::vector<std::int64_t>& rNsIds,
                     const std::map<int, std::vector<std::int64_t>>& rNodeSets,
                     const std::vector<std::string>& rSsNames,
                     const std::vector<std::int64_t>& rSsIds,
                     const std::map<int, std::vector<std::int64_t>>& rSideElems,
                     const std::map<int, std::vector<std::int64_t>>& rSideSides) {
    const std::vector<std::int64_t> bases = detail::block_bases(rMesh);
    const std::int64_t total_cells = bases.empty() ? 0 : bases.back();

    // --- element blocks -> Cell regions -----------------------------------
    // One per connect{k}. The name and tag come from per-block arrays, which is
    // what keeps two blocks of the SAME element type distinguishable -- exactly
    // the case a materials assignment depends on.
    for (std::size_t k = 0; k < rBlockTypes.size() && k + 1 < bases.size(); ++k) {
        std::vector<std::int64_t> entries;
        entries.reserve(static_cast<std::size_t>(bases[k + 1] - bases[k]));
        for (std::int64_t g = bases[k]; g < bases[k + 1]; ++g)
            entries.push_back(g);
        const std::int64_t tag = k < rEbIds.size() ? rEbIds[k] : -1;
        const auto td = topological_dimension().find(rBlockTypes[k]);
        const int dim = td == topological_dimension().end() ? -1 : td->second;
        rMesh.AddRegion(Region(exo_group_name(rEbNames, rEbIds, k, "Block"), RegionKind::Cell, dim,
                               tag, exo_entries(entries, 1)));
    }

    // --- node sets -> Point regions ---------------------------------------
    const std::size_t npts = rMesh.NumPoints();
    std::size_t ns_index = 0;
    for (const auto& [k, nodes] : rNodeSets) {
        (void)k;
        std::vector<std::int64_t> entries;
        entries.reserve(nodes.size());
        for (std::int64_t id : nodes) {
            const std::int64_t p = id - 1;  // Exodus is 1-based
            if (p >= 0 && static_cast<std::size_t>(p) < npts)
                entries.push_back(p);
        }
        const std::int64_t tag = ns_index < rNsIds.size() ? rNsIds[ns_index] : -1;
        rMesh.AddRegion(Region(exo_group_name(rNsNames, rNsIds, ns_index, "Nodeset"),
                               RegionKind::Point, -1, tag, exo_entries(entries, 1)));
        ++ns_index;
    }

    // --- side sets -> Side regions ----------------------------------------
    // (global cell, local facet) pairs. Exodus numbers an element's sides its
    // own way, so the facet goes through exo_face_index; an unmappable pair is
    // skipped rather than stored pointing at the wrong face.
    std::size_t ss_index = 0;
    for (const auto& [k, elems] : rSideElems) {
        auto sit = rSideSides.find(k);
        const std::vector<std::int64_t> empty;
        const std::vector<std::int64_t>& sides = sit == rSideSides.end() ? empty : sit->second;
        std::vector<std::int64_t> pairs;
        pairs.reserve(elems.size() * 2);
        for (std::size_t i = 0; i < elems.size() && i < sides.size(); ++i) {
            const std::int64_t g = elems[i] - 1;  // Exodus is 1-based
            if (g < 0 || g >= total_cells)
                continue;
            const auto [b, row] = detail::global_to_block_row(bases, g);
            (void)row;
            if (b == static_cast<std::size_t>(-1) || b >= rBlockTypes.size())
                continue;
            const int facet = exo_face_index(rBlockTypes[b], static_cast<int>(sides[i]));
            if (facet < 0)
                continue;
            pairs.push_back(g);
            pairs.push_back(facet);
        }
        const std::int64_t tag = ss_index < rSsIds.size() ? rSsIds[ss_index] : -1;
        rMesh.AddRegion(Region(exo_group_name(rSsNames, rSsIds, ss_index, "Sideset"),
                               RegionKind::Side, -1, tag, exo_entries(pairs, 2)));
        ++ss_index;
    }
}

// Turn Exodus per-element attributes into `cell_data` under
// `kExodusAttributePrefix`.
//
// `attrib{k}` is `(num_el_in_blk{k}, num_att_in_blk{k})` -- one column per
// attribute, one row per element -- while cell_data is the other way round: one
// array per *name*, holding one sub-array per cell block. So the columns are
// transposed out here, and a block that does not carry a given attribute is
// filled with NaN. That fill is not just padding: it is the signal `write_exodus`
// reads back to leave the attribute out of that block again, which is what makes
// a mixed file round-trip.
//
// Values are always Float64 (Exodus attributes are floating point by definition,
// and NaN needs somewhere to live) regardless of the on-disk type.
void exo_add_attributes(Mesh& rMesh, const std::vector<int>& rBlockKeys,
                        const std::map<int, NDArray>& rAttribs,
                        const std::map<int, std::vector<std::string>>& rNames) {
    const std::size_t nblocks = rBlockKeys.size();
    if (nblocks == 0 || rAttribs.empty())
        return;

    // name -> block position -> that block's column. `std::map` keeps the
    // cell_data insertion order deterministic across runs and backends.
    std::map<std::string, std::map<std::size_t, std::vector<double>>> by_name;
    for (std::size_t b = 0; b < nblocks; ++b) {
        auto ait = rAttribs.find(rBlockKeys[b]);
        if (ait == rAttribs.end())
            continue;
        const NDArray& att = ait->second;
        const std::size_t rows = att.Shape().empty() ? 0 : att.Shape()[0];
        const std::size_t cols = att.Shape().size() >= 2 ? att.Shape()[1] : 1;
        auto nit = rNames.find(rBlockKeys[b]);
        for (std::size_t c = 0; c < cols; ++c) {
            std::string name;
            if (nit != rNames.end() && c < nit->second.size())
                name = nit->second[c];
            // An unnamed attribute is still an attribute -- SEACAS writes
            // `attrib_name{k}` blank often enough. Naming it by its 1-based
            // column is what lets two blocks agree on which one is "the first".
            if (name.empty())
                name = "attribute" + std::to_string(c + 1);
            std::vector<double> vals(rows, 0.0);
            for (std::size_t r = 0; r < rows; ++r)
                vals[r] = detail::read_double(att, r * cols + c);
            by_name[std::string(kExodusAttributePrefix) + name].emplace(b, std::move(vals));
        }
    }

    for (auto& kv : by_name) {
        std::vector<NDArray> per_block;
        per_block.reserve(nblocks);
        for (std::size_t b = 0; b < nblocks; ++b) {
            const std::size_t n = rMesh.Cells(b).NumCells();
            NDArray arr(DType::Float64, {n});
            auto it = kv.second.find(b);
            for (std::size_t r = 0; r < n; ++r)
                arr.As<double>()[r] = (it != kv.second.end() && r < it->second.size())
                                          ? it->second[r]
                                          : std::numeric_limits<double>::quiet_NaN();
            per_block.push_back(std::move(arr));
        }
        rMesh.AddCellData(kv.first, std::move(per_block));
    }
}

NDArray column_stack(const std::vector<const NDArray*>& rCols) {
    std::size_t n = rCols.empty() || rCols[0]->Shape().empty() ? 0 : rCols[0]->Shape()[0];
    NDArray out(rCols[0]->Dtype(), {n, rCols.size()});
    for (std::size_t c = 0; c < rCols.size(); ++c)
        for (std::size_t i = 0; i < n; ++i) {
            double v = detail::read_double(*rCols[c], i);
            if (out.Dtype() == DType::Float32)
                out.As<float>()[i * rCols.size() + c] = static_cast<float>(v);
            else
                out.As<double>()[i * rCols.size() + c] = v;
        }
    return out;
}

}  // namespace

Mesh read_exodus(const std::string& rPath, ExodusInfo& rInfo, const ReadOptions& rOptions) {
    int ncid;
    check(nc_open(rPath.c_str(), NC_NOWRITE, &ncid), "open");
    struct Closer {
        int mId;
        ~Closer() { nc_close(mId); }
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
    std::map<int, NDArray> pd;                 // idx -> values (selected step)
    std::map<int, std::map<int, NDArray>> cd;  // idx -> block -> values
    struct Block {
        std::string mType;
        NDArray mData;
    };
    std::vector<std::pair<int, Block>> blocks;  // connect{k} in numeric order

    // Region raw material, collected here and turned into regions after the
    // cell blocks are in the mesh -- a Cell/Side region needs global block-major
    // indices, which only exist once the blocks are ordered and added.
    std::vector<std::string> eb_names, ns_names, ss_names;
    std::vector<std::int64_t> eb_ids, ns_ids, ss_ids;
    std::map<int, std::vector<std::int64_t>> node_sets;   // set k -> node ids (1-based)
    std::map<int, std::vector<std::int64_t>> side_elems;  // set k -> element ids (1-based)
    std::map<int, std::vector<std::int64_t>> side_sides;  // set k -> Exodus side numbers

    // Per-element attributes, keyed by the same block number `connect{k}` uses.
    // Kept raw here and turned into cell_data once the blocks are ordered, since
    // a cell_data array is per *mesh block position*, not per file block id.
    std::map<int, NDArray> attribs;                        // block k -> (n_el, n_att)
    std::map<int, std::vector<std::string>> attrib_names;  // block k -> per-column names

    // Resolve the requested time step up front, so an out-of-range request
    // fails before any heavy array is decoded rather than midway through.
    const std::vector<double> time_values = exo_time_values(ncid);
    const std::size_t step = rOptions.ResolveTimeStep(time_values.size());

    for (int varid = 0; varid < nvars; ++varid) {
        char namebuf[NC_MAX_NAME + 1] = {0};
        check(nc_inq_varname(ncid, varid, namebuf), "inq_varname");
        std::string key(namebuf);
        std::vector<std::size_t> dims = var_dims(ncid, varid);

        if (key == "info_records" || key == "qa_records") {
            // Provenance strings. NDArray has no string dtype, so these ride the
            // ExodusInfo side channel (MedInfo's pattern) rather than the mesh.
            // Reading them used to throw, which made every file SEACAS/Cubit/
            // Sierra writes unreadable wherever no Python fallback exists.
            std::vector<std::string> recs = read_names(ncid, varid);
            rInfo.mInfoRecords.insert(rInfo.mInfoRecords.end(), recs.begin(), recs.end());
        } else if (key == "eb_names") {
            eb_names = read_names(ncid, varid);
        } else if (key == "ns_names") {
            ns_names = read_names(ncid, varid);
        } else if (key == "ss_names") {
            ss_names = read_names(ncid, varid);
        } else if (key == "eb_prop1") {
            eb_ids = read_ints(ncid, varid);
        } else if (key == "ns_prop1") {
            ns_ids = read_ints(ncid, varid);
        } else if (key == "ss_prop1") {
            ss_ids = read_ints(ncid, varid);
        } else if (key.rfind("node_ns", 0) == 0) {
            node_sets[exo_suffix(key, 7)] = read_ints(ncid, varid);
        } else if (key.rfind("elem_ss", 0) == 0) {
            side_elems[exo_suffix(key, 7)] = read_ints(ncid, varid);
        } else if (key.rfind("side_ss", 0) == 0) {
            side_sides[exo_suffix(key, 7)] = read_ints(ncid, varid);
        } else if (key.rfind("attrib_name", 0) == 0) {
            // Tested before "attrib", which is its prefix.
            attrib_names[exo_suffix(key, 11)] = read_names(ncid, varid);
        } else if (key.rfind("attrib", 0) == 0) {
            attribs[exo_suffix(key, 6)] =
                read_var(ncid, varid, std::vector<std::size_t>(dims.size(), 0), dims);
        } else if (key.rfind("connect", 0) == 0) {
            std::string elem_type = read_att_text(ncid, varid, "elem_type");
            std::transform(elem_type.begin(), elem_type.end(), elem_type.begin(),
                           [](unsigned char c) { return std::toupper(c); });
            auto it = exodus_to_meshio().find(elem_type);
            if (it == exodus_to_meshio().end())
                throw ReadError("Exodus: unknown element type " + elem_type);
            NDArray conn = read_var(ncid, varid, std::vector<std::size_t>(dims.size(), 0), dims);
            for (std::size_t i = 0; i < conn.Size(); ++i) {
                switch (conn.Dtype()) {
                    case DType::Int32:
                        conn.As<std::int32_t>()[i] -= 1;
                        break;
                    case DType::Int64:
                        conn.As<std::int64_t>()[i] -= 1;
                        break;
                    default:
                        throw ReadError("Exodus: unexpected connectivity dtype");
                }
            }
            int blk = key.size() > 7 ? std::atoi(key.c_str() + 7) : 1;
            blocks.emplace_back(blk, Block{it->second, std::move(conn)});
        } else if (key == "coord") {
            NDArray coord = read_var(ncid, varid, std::vector<std::size_t>(dims.size(), 0), dims);
            std::size_t d = dims.size() >= 1 ? dims[0] : 0;
            std::size_t n = dims.size() >= 2 ? dims[1] : 0;
            NDArray pts(coord.Dtype(), {n, d});
            for (std::size_t c = 0; c < d; ++c)
                for (std::size_t i = 0; i < n; ++i) {
                    if (coord.Dtype() == DType::Float32)
                        pts.As<float>()[i * d + c] = coord.As<float>()[c * n + i];
                    else
                        pts.As<double>()[i * d + c] = coord.As<double>()[c * n + i];
                }
            mesh.AssignPoints(std::move(pts));
            have_coord = true;
        } else if (key == "coordx" || key == "coordy" || key == "coordz") {
            int c = key.back() - 'x';
            NDArray v = read_var(ncid, varid, std::vector<std::size_t>(dims.size(), 0), dims);
            for (std::size_t i = 0; i < num_nodes && i < v.Size(); ++i)
                points_xyz.As<double>()[i * 3 + c] = detail::read_double(v, i);
        } else if (key == "name_nod_var") {
            point_data_names = read_names(ncid, varid);
        } else if (key.rfind("vals_nod_var", 0) == 0) {
            int idx = key.size() == 12 ? 0 : std::atoi(key.c_str() + 12) - 1;
            // dims: (time_step, ...) -> the one requested step
            std::vector<std::size_t> start(dims.size(), 0), count = dims;
            if (!count.empty()) {
                count[0] = 1;
                start[0] = exo_step_offset(step, dims[0], key);
            }
            NDArray v = read_var(ncid, varid, start, count);
            std::vector<std::size_t> shape(dims.begin() + 1, dims.end());
            v.Reshape(shape);
            pd.emplace(idx, std::move(v));
        } else if (key == "name_elem_var") {
            cell_data_names = read_names(ncid, varid);
        } else if (key.rfind("vals_elem_var", 0) == 0) {
            // vals_elem_var(\d+)?(eb(\d+))?
            std::string rest = key.substr(13);
            int idx = 0, block = 0;
            std::size_t eb = rest.find("eb");
            std::string first = eb == std::string::npos ? rest : rest.substr(0, eb);
            if (!first.empty())
                idx = std::atoi(first.c_str()) - 1;
            if (eb != std::string::npos)
                block = std::atoi(rest.c_str() + eb + 2) - 1;
            std::vector<std::size_t> start(dims.size(), 0), count = dims;
            if (!count.empty()) {
                count[0] = 1;
                start[0] = exo_step_offset(step, dims[0], key);
            }
            NDArray v = read_var(ncid, varid, start, count);
            std::vector<std::size_t> shape(dims.begin() + 1, dims.end());
            v.Reshape(shape);
            cd[idx].emplace(block, std::move(v));
        }
        // all other variables (time_whole, coor_names, eb_prop1, ...) ignored
    }

    if (!have_coord)
        mesh.AssignPoints(std::move(points_xyz));

    std::sort(blocks.begin(), blocks.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<std::string> block_types;
    std::vector<int> block_keys;  // the file's block number per mesh block position
    block_types.reserve(blocks.size());
    block_keys.reserve(blocks.size());
    for (auto& b : blocks) {
        block_types.push_back(b.second.mType);
        block_keys.push_back(b.first);
        mesh.AddCellBlock(b.second.mType, std::move(b.second.mData));
    }

    exo_add_regions(mesh, block_types, eb_names, eb_ids, ns_names, ns_ids, node_sets, ss_names,
                    ss_ids, side_elems, side_sides);
    exo_add_attributes(mesh, block_keys, attribs, attrib_names);

    // Point data with X/Y/Z + _R/_Z recombination.
    if (!point_data_names.empty()) {
        Categorized cat = categorize(point_data_names);
        for (const auto& kv : cat.mSingle)
            mesh.AddPointData(kv.first, std::move(pd.at(kv.second)));
        for (std::size_t i = 0; i < cat.mDoubleIdx.size(); ++i)
            mesh.AddPointData(cat.mDoubleName[i], column_stack({&pd.at(cat.mDoubleIdx[i][0]),
                                                                &pd.at(cat.mDoubleIdx[i][1])}));
        for (std::size_t i = 0; i < cat.mTripleIdx.size(); ++i)
            mesh.AddPointData(cat.mTripleName[i], column_stack({&pd.at(cat.mTripleIdx[i][0]),
                                                                &pd.at(cat.mTripleIdx[i][1]),
                                                                &pd.at(cat.mTripleIdx[i][2])}));
    }

    // Cell data: concatenate blocks, then re-split by cell-block sizes.
    if (!cell_data_names.empty() && !cd.empty()) {
        std::vector<std::size_t> sizes;
        for (const auto cb : mesh.CellRange())
            sizes.push_back(cb.NumCells());
        std::size_t name_i = 0;
        for (auto& kv : cd) {
            if (name_i >= cell_data_names.size())
                break;
            const std::string& name = cell_data_names[name_i++];
            // concatenate in block order. The component count comes from the
            // variable's own trailing dimensions -- this used to be hard-coded
            // scalar (`{total}` then `{s}`), which silently truncated a
            // multi-component element variable to its first component. Standard
            // Exodus element variables are scalar per element, so no real
            // SEACAS file exercised that; meshio++'s own writer emits the
            // trailing dims (as the nodal path already did), so it does.
            DType dt = kv.second.begin()->second.Dtype();
            std::size_t ncomp = 1;
            {
                const std::vector<std::size_t>& s0 = kv.second.begin()->second.Shape();
                for (std::size_t d = 1; d < s0.size(); ++d)
                    ncomp *= s0[d];
            }
            std::size_t total = 0;
            for (const auto& b : kv.second)
                total += b.second.Shape().empty() ? 0 : b.second.Shape()[0];
            NDArray all(dt, ncomp > 1 ? std::vector<std::size_t>{total, ncomp}
                                      : std::vector<std::size_t>{total});
            std::size_t off = 0;
            for (const auto& b : kv.second) {
                std::memcpy(all.Data() + off, b.second.Data(), b.second.Nbytes());
                off += b.second.Nbytes();
            }
            // split
            const std::size_t row_bytes = ncomp * dtype_size(dt);
            std::vector<NDArray> out_blocks;
            std::size_t pos = 0;
            for (std::size_t s : sizes) {
                NDArray blk(dt, ncomp > 1 ? std::vector<std::size_t>{s, ncomp}
                                          : std::vector<std::size_t>{s});
                std::memcpy(blk.Data(), all.Data() + pos * row_bytes, s * row_bytes);
                pos += s;
                out_blocks.push_back(std::move(blk));
            }
            mesh.AddCellData(name, std::move(out_blocks));
        }
    }

    // The time of the step actually returned, so the writer's
    // `field_data["exodus:time"]` closes into a round trip rather than being
    // write-only. `read_metadata` still owns the *whole* list of steps -- this
    // is the one value this mesh is a snapshot at.
    if (step < time_values.size()) {
        NDArray tv(DType::Float64, {1});
        *reinterpret_cast<double*>(tv.Data()) = time_values[step];
        mesh.AddFieldData("exodus:time", std::move(tv));
    }

    return mesh;
}

Mesh read_exodus(const std::string& rPath, const ReadOptions& rOptions) {
    // The provenance strings have nowhere to go on this path -- the flat
    // bindings have no `info` slot. Dropping them is what `registry.cpp` already
    // does for MedInfo, and is not a reason to fail the read.
    ExodusInfo info;
    return read_exodus(rPath, info, rOptions);
}

MeshMetadata read_exodus_metadata(const std::string& rPath, const ReadOptions& rOptions) {
    // No native cheap path yet: Exodus's counts live in dimensions that would be
    // cheap to walk, but the cell-block *types* come from per-variable
    // `elem_type` attributes, so a summary still has to visit every connect{k}.
    // What this override does buy is `mTimeValues`, which a full read discards
    // and which is the only way to discover how many steps `mTimeStep` may name.
    ExodusInfo info;
    MeshMetadata meta = metadata_from_mesh(read_exodus(rPath, info, rOptions));
    meta.mFellBackToFullRead = true;

    int ncid;
    check(nc_open(rPath.c_str(), NC_NOWRITE, &ncid), "open");
    struct Closer {
        int mId;
        ~Closer() { nc_close(mId); }
    } closer{ncid};
    meta.mTimeValues = exo_time_values(ncid);
    return meta;
}

void write_exodus(const std::string& rPath, const Mesh& rMesh) {
    int ncid;
    check(nc_create(rPath.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid), "create", true);
    struct Closer {
        int mId;
        ~Closer() { nc_close(mId); }
    } closer{ncid};

    const NDArray& points = rMesh.Points();
    const std::size_t npts = rMesh.NumPoints();
    const std::size_t pdim = rMesh.PointDim();

    // global attributes
    {
        std::string title = detail::kProvenanceTag;
        check(nc_put_att_text(ncid, NC_GLOBAL, "title", title.size(), title.c_str()), "title",
              true);
        float v = 5.1f;
        check(nc_put_att_float(ncid, NC_GLOBAL, "version", NC_FLOAT, 1, &v), "version", true);
        check(nc_put_att_float(ncid, NC_GLOBAL, "api_version", NC_FLOAT, 1, &v), "api_version",
              true);
        long long w = 8;
        check(nc_put_att_longlong(ncid, NC_GLOBAL, "floating_point_word_size", NC_INT64, 1, &w),
              "fpws", true);
    }

    std::size_t total_elems = 0;
    for (const auto cb : rMesh.CellRange())
        total_elems += cb.NumCells();

    int d_nodes, d_dim, d_elem, d_blk, d_ns, d_str, d_line, d_four, d_time;
    check(nc_def_dim(ncid, "num_nodes", npts, &d_nodes), "def num_nodes", true);
    check(nc_def_dim(ncid, "num_dim", pdim, &d_dim), "def num_dim", true);
    check(nc_def_dim(ncid, "num_elem", total_elems, &d_elem), "def num_elem", true);
    check(nc_def_dim(ncid, "num_el_blk", rMesh.NumCellBlocks(), &d_blk), "def num_el_blk", true);
    check(nc_def_dim(ncid, "num_node_sets", 0, &d_ns), "def num_node_sets", true);
    check(nc_def_dim(ncid, "len_string", 33, &d_str), "def len_string", true);
    check(nc_def_dim(ncid, "len_line", 81, &d_line), "def len_line", true);
    check(nc_def_dim(ncid, "four", 4, &d_four), "def four", true);
    check(nc_def_dim(ncid, "time_step", NC_UNLIMITED, &d_time), "def time_step", true);

    // The single time step this writer emits. A `Mesh` is one state, so one
    // step is all there is to write; what changed in v9.9.0 is that its
    // recorded time is no longer hard-coded to 0 -- `field_data["exodus:time"]`
    // (the `<format>:<thing>` convention `med:num` established) supplies it, so
    // a caller writing one frame of a transient solve can label it correctly.
    // A genuine multi-step writer is a separate object with its own lifecycle,
    // the shape `XdmfTimeSeriesWriter` already has; that remains a follow-up.
    {
        int var;
        check(nc_def_var(ncid, "time_whole", NC_FLOAT, 1, &d_time, &var), "def time_whole", true);
        std::size_t start = 0, count = 1;
        float t = 0.0f;
        if (rMesh.HasFieldData("exodus:time")) {
            const NDArray& tv = rMesh.FieldData("exodus:time");
            if (tv.Size() >= 1)
                t = static_cast<float>(detail::read_double(tv, 0));
            if (tv.Size() > 1)
                log::warn(
                    "Exodus: field_data[\"exodus:time\"] has {} values but this writer emits a "
                    "single time step; using the first.",
                    tv.Size());
        }
        check(nc_put_vara_float(ncid, var, &start, &count, &t), "time_whole", true);
    }

    // coor_names
    {
        int dims[2] = {d_dim, d_str};
        int var;
        check(nc_def_var(ncid, "coor_names", NC_CHAR, 2, dims, &var), "coor_names", true);
        const char* names = "XYZ";
        for (std::size_t c = 0; c < pdim && c < 3; ++c) {
            std::size_t start[2] = {c, 0}, count[2] = {1, 1};
            check(nc_put_vara_text(ncid, var, start, count, &names[c]), "coor_names", true);
        }
    }

    // coord (num_dim, num_nodes) = points^T
    {
        int dims[2] = {d_dim, d_nodes};
        int var;
        check(nc_def_var(ncid, "coord", nc_type_of(points.Dtype()), 2, dims, &var), "def coord",
              true);
        NDArray t(points.Dtype(), {pdim, npts});
        for (std::size_t c = 0; c < pdim; ++c)
            for (std::size_t i = 0; i < npts; ++i) {
                if (t.Dtype() == DType::Float32)
                    t.As<float>()[c * npts + i] = points.As<float>()[i * pdim + c];
                else
                    t.As<double>()[c * npts + i] = points.As<double>()[i * pdim + c];
            }
        if (t.Size() > 0)
            check(nc_put_var(ncid, var, t.Data()), "coord", true);
    }

    // eb_prop1
    {
        int var;
        check(nc_def_var(ncid, "eb_prop1", NC_INT, 1, &d_blk, &var), "eb_prop1", true);
        std::vector<int> ids(rMesh.NumCellBlocks());
        for (std::size_t k = 0; k < ids.size(); ++k)
            ids[k] = static_cast<int>(k);
        if (!ids.empty())
            check(nc_put_var_int(ncid, var, ids.data()), "eb_prop1", true);
    }

    // eb_names, recovered from the Cell regions the reader derives from them.
    // Written only when at least one block has a name; an all-blank array is
    // what SEACAS emits when it has none, and omitting it keeps output for a
    // region-less mesh byte-identical to pre-v9.9.0.
    {
        const std::vector<std::string> block_names = exo_block_names_from_regions(rMesh);
        const bool any = std::any_of(block_names.begin(), block_names.end(),
                                     [](const std::string& s) { return !s.empty(); });
        if (any) {
            int dims[2] = {d_blk, d_str};
            int var;
            check(nc_def_var(ncid, "eb_names", NC_CHAR, 2, dims, &var), "def eb_names", true);
            for (std::size_t k = 0; k < block_names.size(); ++k) {
                if (block_names[k].empty())
                    continue;
                std::size_t start[2] = {k, 0};
                std::size_t count[2] = {1, std::min<std::size_t>(block_names[k].size(), 33)};
                check(nc_put_vara_text(ncid, var, start, count, block_names[k].c_str()), "eb_names",
                      true);
            }
        }
    }

    // connectivity blocks
    std::vector<int> block_elem_dims(rMesh.NumCellBlocks(), -1);
    for (std::size_t k = 0; k < rMesh.NumCellBlocks(); ++k) {
        const auto cb = rMesh.Cells(k);
        auto it = meshio_to_exodus().find(cb.Type());
        if (it == meshio_to_exodus().end())
            throw WriteError("Exodus: unsupported cell type " + cb.Type());
        const NDArray& conn = cb.Conn();
        std::string dim1 = "num_el_in_blk" + std::to_string(k + 1);
        std::string dim2 = "num_nod_per_el" + std::to_string(k + 1);
        int d1, d2;
        check(nc_def_dim(ncid, dim1.c_str(), cb.NumCells(), &d1), "blk dim", true);
        check(nc_def_dim(ncid, dim2.c_str(), detail::cols(conn), &d2), "blk dim", true);
        block_elem_dims[k] = d1;
        int dims[2] = {d1, d2};
        int var;
        std::string vname = "connect" + std::to_string(k + 1);
        check(nc_def_var(ncid, vname.c_str(), nc_type_of(conn.Dtype()), 2, dims, &var),
              "def connect", true);
        check(nc_put_att_text(ncid, var, "elem_type", it->second.size(), it->second.c_str()),
              "elem_type", true);
        NDArray shifted(conn.Dtype(), conn.Shape());
        for (std::size_t i = 0; i < conn.Size(); ++i) {
            std::int64_t v = detail::read_int(conn, i) + 1;
            switch (shifted.Dtype()) {
                case DType::Int32:
                    shifted.As<std::int32_t>()[i] = static_cast<std::int32_t>(v);
                    break;
                case DType::Int64:
                    shifted.As<std::int64_t>()[i] = v;
                    break;
                default:
                    throw WriteError("Exodus: unexpected connectivity dtype");
            }
        }
        if (shifted.Size() > 0)
            check(nc_put_var(ncid, var, shifted.Data()), "connect", true);
    }

    // Per-element attributes: every cell_data array whose name starts with
    // `kExodusAttributePrefix` goes back out as a column of `attrib{k}`, named in
    // `attrib_name{k}`. Everything else in cell_data is left alone (this writer
    // emits no `vals_elem_var` at all -- a separate, pre-existing gap).
    {
        const std::string prefix(kExodusAttributePrefix);
        std::vector<std::string> att_names;
        for (const auto& name : rMesh.CellDataNames())
            if (name.rfind(prefix, 0) == 0)
                att_names.push_back(name);

        for (std::size_t k = 0; k < rMesh.NumCellBlocks() && !att_names.empty(); ++k) {
            const std::size_t n = rMesh.Cells(k).NumCells();
            std::vector<std::string> names;
            std::vector<const NDArray*> cols;
            for (const auto& full : att_names) {
                const NDArray& arr = rMesh.CellData(full, k);
                // Product of ALL trailing dims, not just `detail::cols()`:
                // an `(n,1,3)` array has cols == 1 and would otherwise slip
                // past and be silently truncated to its first component. This
                // matches `_exodus.py`'s `prod(shape[1:]) != 1` twin.
                std::size_t trailing = 1;
                for (std::size_t d = 1; d < arr.Shape().size(); ++d)
                    trailing *= arr.Shape()[d];
                if (trailing != 1)
                    throw WriteError("Exodus: element attribute '" + full +
                                     "' must be scalar (one value per element)");
                // A block whose values are all non-finite never carried this
                // attribute -- that NaN is exactly what the reader fills in for a
                // block the file left it out of. Skipping it here is what makes a
                // mixed file round-trip instead of gaining NaN attributes.
                bool any_finite = false;
                for (std::size_t r = 0; r < n && !any_finite; ++r)
                    any_finite = std::isfinite(detail::read_double(arr, r));
                if (!any_finite)
                    continue;
                names.push_back(full.substr(prefix.size()));
                cols.push_back(&arr);
            }
            if (names.empty())
                continue;

            std::string dim_el = "num_el_in_blk" + std::to_string(k + 1);
            std::string dim_att = "num_att_in_blk" + std::to_string(k + 1);
            int d_el, d_att;
            check(nc_inq_dimid(ncid, dim_el.c_str(), &d_el), "num_el_in_blk", true);
            check(nc_def_dim(ncid, dim_att.c_str(), names.size(), &d_att), "def num_att_in_blk",
                  true);

            int dims[2] = {d_el, d_att};
            int var;
            std::string vname = "attrib" + std::to_string(k + 1);
            check(nc_def_var(ncid, vname.c_str(), NC_DOUBLE, 2, dims, &var), "def attrib", true);
            NDArray flat(DType::Float64, {n, names.size()});
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t c = 0; c < names.size(); ++c)
                    flat.As<double>()[r * names.size() + c] = detail::read_double(*cols[c], r);
            if (flat.Size() > 0)
                check(nc_put_var(ncid, var, flat.Data()), "attrib", true);

            int name_dims[2] = {d_att, d_str};
            int name_var;
            std::string nname = "attrib_name" + std::to_string(k + 1);
            check(nc_def_var(ncid, nname.c_str(), NC_CHAR, 2, name_dims, &name_var),
                  "def attrib_name", true);
            for (std::size_t c = 0; c < names.size(); ++c) {
                std::size_t start[2] = {c, 0};
                std::size_t count[2] = {1, std::min<std::size_t>(names[c].size(), 33)};
                if (count[1] > 0)
                    check(nc_put_vara_text(ncid, name_var, start, count, names[c].c_str()),
                          "attrib_name", true);
            }
        }
    }

    // Cell data -> element variables (`name_elem_var` + one
    // `vals_elem_var{j}eb{k}` per (variable, block)), new in v9.9.0: this
    // writer previously emitted none at all, so every ordinary `cell_data`
    // array was silently dropped while `point_data` round-tripped. The
    // `kExodusAttributePrefix` arrays are excluded -- they are constant-in-time
    // per-element *attributes* and already went out as `attrib{k}` above, which
    // is a different Exodus concept and the reason that prefix has to be
    // explicit.
    {
        const std::string prefix(kExodusAttributePrefix);
        std::vector<std::string> var_names;
        for (const auto& name : rMesh.CellDataNames())
            if (name.rfind(prefix, 0) != 0)
                var_names.push_back(name);

        if (!var_names.empty()) {
            int d_nev;
            check(nc_def_dim(ncid, "num_elem_var", var_names.size(), &d_nev), "num_elem_var", true);
            int name_var;
            {
                int dims[2] = {d_nev, d_str};
                check(nc_def_var(ncid, "name_elem_var", NC_CHAR, 2, dims, &name_var),
                      "def name_elem_var", true);
            }
            // The truth table: which variable exists on which block. Every
            // cell_data array covers every block by the uniform API's own
            // invariant, so it is all ones -- but real Exodus readers expect
            // the variable to be present, so it is written rather than assumed.
            {
                int dims[2] = {d_blk, d_nev};
                int tab;
                check(nc_def_var(ncid, "elem_var_tab", NC_INT, 2, dims, &tab), "def elem_var_tab",
                      true);
                std::vector<int> ones(rMesh.NumCellBlocks() * var_names.size(), 1);
                if (!ones.empty())
                    check(nc_put_var_int(ncid, tab, ones.data()), "elem_var_tab", true);
            }

            for (std::size_t j = 0; j < var_names.size(); ++j) {
                const std::string& name = var_names[j];
                {
                    std::size_t start[2] = {j, 0};
                    std::size_t count[2] = {1, std::min<std::size_t>(name.size(), 33)};
                    if (count[1] > 0)
                        check(nc_put_vara_text(ncid, name_var, start, count, name.c_str()),
                              "name_elem_var", true);
                }
                if (rMesh.CellDataNumBlocks(name) != rMesh.NumCellBlocks()) {
                    log::warn(
                        "Exodus: cell_data '{}' covers {} of {} blocks; not written as an element "
                        "variable.",
                        name, rMesh.CellDataNumBlocks(name), rMesh.NumCellBlocks());
                    continue;
                }
                for (std::size_t k = 0; k < rMesh.NumCellBlocks(); ++k) {
                    const NDArray& data = rMesh.CellData(name, k);
                    // Trailing dims become extra netCDF dimensions, exactly as
                    // the nodal-variable path above already does for a
                    // multi-component point field -- so a vector cell field
                    // round-trips through this reader (and the Python twin's,
                    // which reshapes from the same dims) rather than being
                    // dropped. Standard Exodus element variables are scalar
                    // per element, so a k>1 array here is a meshio++ extension
                    // of the same kind the nodal path already is.
                    std::vector<int> dims = {d_time, block_elem_dims[k]};
                    for (std::size_t i = 1; i < data.Shape().size(); ++i) {
                        std::string dn = "dim_elem_var" + std::to_string(j) + "_" +
                                         std::to_string(k) + "_" + std::to_string(i);
                        int di;
                        check(nc_def_dim(ncid, dn.c_str(), data.Shape()[i], &di), "cd dim", true);
                        dims.push_back(di);
                    }
                    int var;
                    std::string vname =
                        "vals_elem_var" + std::to_string(j + 1) + "eb" + std::to_string(k + 1);
                    check(nc_def_var(ncid, vname.c_str(), nc_type_of(data.Dtype()),
                                     static_cast<int>(dims.size()), dims.data(), &var),
                          "def vals_elem_var", true);
                    check(nc_def_var_fill(ncid, var, NC_NOFILL, nullptr), "nofill", true);
                    std::vector<std::size_t> startv(dims.size(), 0), countv;
                    countv.push_back(1);
                    for (std::size_t s : data.Shape())
                        countv.push_back(s);
                    if (data.Size() > 0)
                        check(nc_put_vara(ncid, var, startv.data(), countv.data(), data.Data()),
                              "vals_elem_var", true);
                }
            }
        }
    }

    // point data
    if (rMesh.NumPointData() > 0) {
        int d_nnv;
        check(nc_def_dim(ncid, "num_nod_var", rMesh.NumPointData(), &d_nnv), "num_nod_var", true);
        int name_var;
        {
            int dims[2] = {d_nnv, d_str};
            check(nc_def_var(ncid, "name_nod_var", NC_CHAR, 2, dims, &name_var), "name_nod_var",
                  true);
        }
        std::size_t k = 0;
        // Sorted key order: assigns the on-disk variable index (slot k) and
        // name deterministically, independent of the map's storage order.
        for (const auto& name : rMesh.PointDataNames()) {
            std::size_t start[2] = {k, 0};
            std::size_t count[2] = {1, std::min<std::size_t>(name.size(), 33)};
            if (count[1] > 0)
                check(nc_put_vara_text(ncid, name_var, start, count, name.c_str()), "name_nod_var",
                      true);

            const NDArray& data = rMesh.PointData(name);
            std::vector<int> dims = {d_time};
            for (std::size_t i = 0; i < data.Shape().size(); ++i) {
                std::string dn = "dim_nod_var" + std::to_string(k) + std::to_string(i);
                int di;
                check(nc_def_dim(ncid, dn.c_str(), data.Shape()[i], &di), "pd dim", true);
                dims.push_back(di);
            }
            int var;
            std::string vname = "vals_nod_var" + std::to_string(k + 1);
            check(nc_def_var(ncid, vname.c_str(), nc_type_of(data.Dtype()),
                             static_cast<int>(dims.size()), dims.data(), &var),
                  "def vals_nod_var", true);
            check(nc_def_var_fill(ncid, var, NC_NOFILL, nullptr), "nofill", true);
            std::vector<std::size_t> startv(dims.size(), 0), countv;
            countv.push_back(1);
            for (std::size_t s : data.Shape())
                countv.push_back(s);
            if (data.Size() > 0)
                check(nc_put_vara(ncid, var, startv.data(), countv.data(), data.Data()),
                      "vals_nod_var", true);
            ++k;
        }
    }

    // Node sets (point_sets) are not representable in the conversion layer;
    // the shim routes meshes with point_sets to the Python writer.
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
