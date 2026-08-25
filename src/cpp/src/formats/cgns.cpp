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
#include <cstring>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/cgns.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/face_mesh.hpp"
#include "meshioplusplus/detail/format_compat.hpp"
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {

namespace {

// --- node encoding: CGNS/ADF-over-HDF5's own convention, not an ad hoc one -
//
// Every CGNS node is an HDF5 group carrying four fixed attributes (name/
// label/type/flags) plus, when it has a payload, a dataset literally named
// " data" (leading space -- cgnslib's own ADFH.c: `#define D_DATA " data"`).
// h5::write_attr_string (variable-length UTF-8) and MED's write_attr_bytes
// (NULLPAD, size = strlen) both write the wrong encoding for this: CGNS uses
// fixed-size (33/33/3 bytes) NULLTERM strings.

void cgns_write_attr_str(hid_t loc, const std::string& rName, const std::string& rValue,
                         std::size_t size) {
    h5::Hid t(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(t, size);
    H5Tset_strpad(t, H5T_STR_NULLTERM);
    h5::Hid space(H5Screate(H5S_SCALAR), H5Sclose);
    h5::Hid a(H5Acreate2(loc, rName.c_str(), t, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw WriteError(detail::format_compat("CGNS: could not create attribute '{}'", rName));
    std::string buf(size, '\0');
    std::memcpy(buf.data(), rValue.data(), std::min(rValue.size(), size));
    H5Awrite(a, t, buf.data());
}

void cgns_write_attr_flags(hid_t loc) {
    hsize_t dim = 1;
    h5::Hid space(H5Screate_simple(1, &dim, nullptr), H5Sclose);
    h5::Hid a(H5Acreate2(loc, "flags", H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose);
    if (!a.Valid())
        throw WriteError("CGNS: could not create 'flags' attribute");
    std::int32_t v = 1;
    H5Awrite(a, H5T_NATIVE_INT32, &v);
}

// The ADF two-character type code for a node's " data" payload dtype.
std::string cgns_type_code(DType dt) {
    switch (dt) {
        case DType::Int32:
            return "I4";
        case DType::Int64:
            return "I8";
        case DType::Float32:
            return "R4";
        case DType::Float64:
            return "R8";
        default:
            throw WriteError("CGNS: unsupported dtype for a CGNS node payload");
    }
}

void cgns_write_node_attrs(hid_t loc, const std::string& rName, const std::string& rLabel,
                           const std::string& rType) {
    cgns_write_attr_str(loc, "name", rName, 33);
    cgns_write_attr_str(loc, "label", rLabel, 33);
    cgns_write_attr_str(loc, "type", rType, 3);
    cgns_write_attr_flags(loc);
}

// The file and every group must track HDF5 link (and attribute) *creation*
// order -- cgnslib's has_child/has_data (ADFH.c) iterate H5_INDEX_CRT_ORDER
// with no name-order fallback, so a file created with the library defaults
// (H5P_DEFAULT) is invisible to a real CGNS reader even though every node
// and attribute is individually well-formed.

h5::Hid cgns_create_file(const std::string& rPath) {
    h5::Hid fcpl(H5Pcreate(H5P_FILE_CREATE), H5Pclose);
    H5Pset_link_creation_order(fcpl, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
    H5Pset_attr_creation_order(fcpl, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
    hid_t f = H5Fcreate(rPath.c_str(), H5F_ACC_TRUNC, fcpl, H5P_DEFAULT);
    if (f < 0)
        throw WriteError(detail::format_compat("CGNS: could not create file '{}'", rPath));
    return h5::Hid(f, H5Fclose);
}

h5::Hid cgns_create_group(hid_t loc, const std::string& rName) {
    h5::Hid gcpl(H5Pcreate(H5P_GROUP_CREATE), H5Pclose);
    H5Pset_link_creation_order(gcpl, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
    H5Pset_attr_creation_order(gcpl, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
    hid_t g = H5Gcreate2(loc, rName.c_str(), H5P_DEFAULT, gcpl, H5P_DEFAULT);
    if (g < 0)
        throw WriteError(detail::format_compat("CGNS: could not create group '{}'", rName));
    return h5::Hid(g, H5Gclose);
}

// A CHAR (ADF "C1") payload: `rS`'s bytes, zero-padded/truncated to exactly
// `total` bytes (NOT necessarily NUL-terminated -- CGNS text data like
// ZoneType's "Unstructured" carries no terminator, only its declared length).
NDArray cgns_padded_int8(const std::string& rS, std::size_t total) {
    NDArray out(DType::Int8, {total});
    std::int8_t* dst = out.As<std::int8_t>();
    std::memset(dst, 0, total);
    std::memcpy(dst, rS.data(), std::min(rS.size(), total));
    return out;
}

std::string cgns_hdf5_version_string() {
    unsigned maj = 0, min = 0, rel = 0;
    H5get_libversion(&maj, &min, &rel);
    std::ostringstream os;
    os << "HDF5 Version " << maj << "." << min << "." << rel;
    return os.str();
}

// --- meshio++ <-> CGNS ElementType_t table ----------------------------------
//
// Derived from the SIDS element-numbering-conventions edge/face tables and
// cross-checked against VTK's own CGNS translator tables (doc/formats/
// cgns.md records the sources and the per-type derivation). Permutations are
// involutions (a block swap), applied identically on read and write:
// `dst[c] = src[p[c]]`, exactly `med_node_perm()`'s convention.
//
// Deliberately NOT here (a `WriteError`/`ReadError` names them instead of
// guessing): the cubic/quartic Lagrange families (line4, triangle10, quad16,
// tetra20, hexahedron64, ...), whose SIDS interior-node order is shown only
// in image-only figures with no text source to verify against, and anything
// ragged (polygon/polyhedron), which CGNS has no fixed-size representation
// for at all (MIXED/NGON_n/NFACE_n sections are not written by meshio++).
struct CgnsTypeInfo {
    std::string mCgnsName;
    int mCode;
    std::vector<int> mPerm;  // empty = identity
};

const std::unordered_map<std::string, CgnsTypeInfo>& cgns_type_table() {
    static const std::unordered_map<std::string, CgnsTypeInfo> m = {
        {"vertex", {"NODE", 2, {}}},
        {"line", {"BAR_2", 3, {}}},
        {"line3", {"BAR_3", 4, {}}},
        {"triangle", {"TRI_3", 5, {}}},
        {"triangle6", {"TRI_6", 6, {}}},
        {"quad", {"QUAD_4", 7, {}}},
        {"quad8", {"QUAD_8", 8, {}}},
        {"quad9", {"QUAD_9", 9, {}}},
        {"tetra", {"TETRA_4", 10, {}}},
        {"tetra10", {"TETRA_10", 11, {}}},
        {"pyramid", {"PYRA_5", 12, {}}},
        {"pyramid14", {"PYRA_14", 13, {}}},
        {"wedge", {"PENTA_6", 14, {}}},
        {"wedge15", {"PENTA_15", 15, {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11}}},
        {"wedge18",
         {"PENTA_18", 16, {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11, 15, 16, 17}}},
        {"hexahedron", {"HEXA_8", 17, {}}},
        {"hexahedron20",
         {"HEXA_20", 18, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15}}},
        {"hexahedron27", {"HEXA_27", 19, {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 16, 17,
                                          18, 19, 12, 13, 14, 15, 24, 22, 21, 23, 20, 25, 26}}},
        // PYRA_13's code (21) is non-monotonic -- appended to ElementType_t
        // after MIXED in the CGNS enum -- not a typo.
        {"pyramid13", {"PYRA_13", 21, {}}},
    };
    return m;
}

const std::unordered_map<int, std::string>& cgns_code_to_meshio() {
    static const std::unordered_map<int, std::string> m = [] {
        std::unordered_map<int, std::string> out;
        for (const auto& kv : cgns_type_table())
            out.emplace(kv.second.mCode, kv.first);
        return out;
    }();
    return m;
}

const std::unordered_map<int, std::string>& cgns_code_to_name() {
    static const std::unordered_map<int, std::string> m = {
        {2, "NODE"},      {3, "BAR_2"},     {4, "BAR_3"},    {5, "TRI_3"},     {6, "TRI_6"},
        {7, "QUAD_4"},    {8, "QUAD_8"},    {9, "QUAD_9"},   {10, "TETRA_4"},  {11, "TETRA_10"},
        {12, "PYRA_5"},   {13, "PYRA_14"},  {14, "PENTA_6"}, {15, "PENTA_15"}, {16, "PENTA_18"},
        {17, "HEXA_8"},   {18, "HEXA_20"},  {19, "HEXA_27"}, {20, "MIXED"},    {21, "PYRA_13"},
        {22, "NGON_n"},   {23, "NFACE_n"},  {24, "BAR_4"},   {26, "TRI_10"},   {28, "QUAD_16"},
        {30, "TETRA_20"}, {36, "PENTA_40"}, {39, "HEXA_64"},
    };
    return m;
}

// Apply an optional column permutation and additive shift to an (n, k)
// row-major connectivity array (CGNS is element-major/row-major, unlike
// MED's Fortran storage, so no transpose is needed -- only the permutation).
// Self-inverse (`p == p^-1` for every table entry above), so the same `pPerm`
// serves both write (shift=+1) and read (shift=-1).
NDArray cgns_permute_conn(const NDArray& rConn, std::size_t n, std::size_t k, std::int64_t shift,
                          const std::vector<int>* pPerm) {
    const int* p = (pPerm && pPerm->size() == k) ? pPerm->data() : nullptr;
    NDArray out(rConn.Dtype(), {n, k});
    detail::dispatch_dtype(rConn.Dtype(), [&]<class T>() {
        const T* src = rConn.As<T>();
        T* dst = out.As<T>();
        const T s = static_cast<T>(shift);
        parallel_for_bw(n, [&](std::size_t i) {
            for (std::size_t c = 0; c < k; ++c) {
                const std::size_t sc = p ? static_cast<std::size_t>(p[c]) : c;
                dst[i * k + c] = static_cast<T>(src[i * k + sc] + s);
            }
        });
    });
    return out;
}

// --- legacy (pre-v9.8.0) layout: read-only, byte-for-byte the old reader ---
Mesh cgns_read_legacy(hid_t f) {
    if (!h5::exists(f, "Base"))
        throw ReadError("Expected \"Base\" in file. Malformed CGNS?");
    h5::Hid base = h5::open_group(f, "Base");
    if (!h5::exists(base, "Zone1"))
        throw ReadError("Expected \"Zone1\" in \"Base\". Malformed CGNS?");
    h5::Hid zone = h5::open_group(base, "Zone1");

    h5::Hid coords = h5::open_group(zone, "GridCoordinates");
    h5::Hid gx = h5::open_group(coords, "CoordinateX");
    h5::Hid gy = h5::open_group(coords, "CoordinateY");
    h5::Hid gz = h5::open_group(coords, "CoordinateZ");
    NDArray x = h5::read_dataset(gx, " data");
    NDArray y = h5::read_dataset(gy, " data");
    NDArray z = h5::read_dataset(gz, " data");

    const std::size_t n = x.Shape().empty() ? 0 : x.Shape()[0];
    Mesh mesh;
    NDArray pts(DType::Float64, {n, 3});
    double* pp = pts.As<double>();
    for (std::size_t i = 0; i < n; ++i) {
        pp[i * 3 + 0] = detail::read_double(x, i);
        pp[i * 3 + 1] = detail::read_double(y, i);
        pp[i * 3 + 2] = detail::read_double(z, i);
    }
    mesh.AssignPoints(std::move(pts));

    h5::Hid elems = h5::open_group(zone, "GridElements");
    h5::Hid rng = h5::open_group(elems, "ElementRange");
    h5::Hid conn = h5::open_group(elems, "ElementConnectivity");
    NDArray range = h5::read_dataset(rng, " data");
    NDArray flat = h5::read_dataset(conn, " data");

    if (range.Size() < 2)
        throw ReadError("CGNS: malformed ElementRange");
    std::int64_t idx_max = detail::read_int(range, 1);
    if (idx_max <= 0 || flat.Size() % static_cast<std::size_t>(idx_max) != 0)
        throw ReadError("CGNS: malformed ElementConnectivity");
    std::size_t k = flat.Size() / static_cast<std::size_t>(idx_max);
    if (k != 4)
        throw ReadError("Can only read tetrahedra.");

    NDArray cells(flat.Dtype(), {static_cast<std::size_t>(idx_max), k});
    for (std::size_t i = 0; i < flat.Size(); ++i) {
        std::int64_t v = detail::read_int(flat, i) - 1;
        switch (cells.Dtype()) {
            case DType::Int32:
                cells.As<std::int32_t>()[i] = static_cast<std::int32_t>(v);
                break;
            case DType::Int64:
                cells.As<std::int64_t>()[i] = v;
                break;
            case DType::UInt32:
                cells.As<std::uint32_t>()[i] = static_cast<std::uint32_t>(v);
                break;
            case DType::UInt64:
                cells.As<std::uint64_t>()[i] = static_cast<std::uint64_t>(v);
                break;
            default:
                throw ReadError("CGNS: unexpected connectivity dtype");
        }
    }
    mesh.AddCellBlock("tetra", std::move(cells));
    return mesh;
}

// Every reserved payload dataset this format writes (" data", " format",
// " hdf5version") starts with a leading space by the ADF-over-HDF5
// convention; no legitimate CGNS *node* (a group) is named that way, so this
// is what lets group-enumeration loops below skip a location's own payload
// dataset without mis-opening it as a group.
bool cgns_is_reserved_name(const std::string& rName) {
    return !rName.empty() && rName[0] == ' ';
}

// The two CGNS ElementType_t codes that describe a cell by its FACES rather
// than by a fixed node count.
constexpr int kCgnsNgon = 22;   // NGON_n:  a face, given as a node list
constexpr int kCgnsNface = 23;  // NFACE_n: a cell, given as SIGNED face ids

// meshio++ spells a jagged 2D block `polygon`/`polygon2`/`polygon<N>` and a 3D
// one `polyhedron<N>`. Both families map to the face-based sections, so the
// test is on the NAME prefix -- see the note at the writer's validation loop
// for why `IsRagged()` is the wrong question.
bool cgns_is_polygon_type(const std::string& rType) {
    return rType.rfind("polygon", 0) == 0;
}

bool cgns_is_polyhedron_type(const std::string& rType) {
    return rType.rfind("polyhedron", 0) == 0;
}

/// One NGON_n / NFACE_n section, read into flat CSR.
struct CgnsPolySection {
    std::string mName;
    int mCode = 0;
    std::int64_t mFirst = 0, mLast = 0;
    std::vector<std::int64_t> mOffsets;  ///< `count + 1` entries into mData
    std::vector<std::int64_t> mData;     ///< 1-based node ids, or SIGNED face element ids
    std::size_t Count() const { return mOffsets.empty() ? 0 : mOffsets.size() - 1; }
};

/**
 * @brief Read one face-based element section.
 *
 * meshio++ writes the CGNS >= 4.0 layout (`ElementStartOffset` beside
 * `ElementConnectivity`) and reads only that. A CGNS 3.x file instead prefixes
 * each row with its own length inline, and normalising the two is exactly what
 * `cg_poly_elements_read` exists for -- so a 3.x file is refused **by name**,
 * pointing at the optional cgnslib backend, rather than misread.
 */
CgnsPolySection cgns_read_poly_section(h5::Hid& rSect, const std::string& rName, int Code) {
    CgnsPolySection out;
    out.mName = rName;
    out.mCode = Code;

    h5::Hid rng = h5::open_group(rSect, "ElementRange");
    NDArray range = h5::read_dataset(rng, " data");
    out.mFirst = detail::read_int(range, 0);
    out.mLast = detail::read_int(range, 1);
    const std::size_t n = out.mLast >= out.mFirst
                              ? static_cast<std::size_t>(out.mLast - out.mFirst + 1)
                              : 0;

    bool has_offset = false;
    for (const std::string& l : h5::group_links(rSect))
        if (l == "ElementStartOffset")
            has_offset = true;
    if (!has_offset)
        throw ReadError(detail::format_compat(
            "CGNS: section '{}' is {} but has no ElementStartOffset, so it uses the CGNS 3.x "
            "inline-length layout; meshio++'s own reader handles only the 4.0 layout -- rebuild "
            "with -DMESHIOPLUSPLUS_WITH_CGNSLIB=ON to read it",
            rName, Code == kCgnsNgon ? "NGON_n" : "NFACE_n"));

    h5::Hid og = h5::open_group(rSect, "ElementStartOffset");
    NDArray off = h5::read_dataset(og, " data");
    if (off.Size() != n + 1)
        throw ReadError(detail::format_compat(
            "CGNS: section '{}' declares {} elements but ElementStartOffset has {} entries "
            "(expected {})",
            rName, n, off.Size(), n + 1));
    out.mOffsets.resize(off.Size());
    for (std::size_t i = 0; i < off.Size(); ++i)
        out.mOffsets[i] = detail::read_int(off, i);

    h5::Hid cg = h5::open_group(rSect, "ElementConnectivity");
    NDArray conn = h5::read_dataset(cg, " data");
    if (!out.mOffsets.empty() &&
        out.mOffsets.back() != static_cast<std::int64_t>(conn.Size()))
        throw ReadError(detail::format_compat(
            "CGNS: section '{}' has ElementStartOffset ending at {} but ElementConnectivity has "
            "{} entries",
            rName, out.mOffsets.back(), conn.Size()));
    out.mData.resize(conn.Size());
    for (std::size_t i = 0; i < conn.Size(); ++i)
        out.mData[i] = detail::read_int(conn, i);
    return out;
}

/**
 * @brief Write one face-based element section.
 *
 * CGNS >= 4.0 stores a variable-length section as `ElementStartOffset`
 * (`nElems + 1` cumulative offsets) beside `ElementConnectivity`; CGNS 3.x
 * instead prefixed each row with its own length inline. meshio++ writes the
 * 4.0 layout **only** -- a writer picks its layout, so the version split that
 * makes `cg_poly_elements_read` worth having on the READ side simply does not
 * arise here.
 */
void cgns_write_poly_section(h5::Hid& rZone, const std::string& rName, int Code,
                             std::int64_t First, std::int64_t Last,
                             const std::vector<std::int64_t>& rOffsets,
                             const std::vector<std::int64_t>& rData, int GzipLevel) {
    h5::Hid sect = cgns_create_group(rZone, rName);
    cgns_write_node_attrs(sect, rName, "Elements_t", "I4");
    {
        NDArray sdata(DType::Int32, {2});
        sdata.As<std::int32_t>()[0] = Code;
        sdata.As<std::int32_t>()[1] = 0;  // ElementSizeBoundary: unsorted
        h5::write_dataset(sect, " data", sdata);
    }
    h5::Hid rng = cgns_create_group(sect, "ElementRange");
    cgns_write_node_attrs(rng, "ElementRange", "IndexRange_t", cgns_type_code(DType::Int64));
    {
        NDArray rdata(DType::Int64, {2});
        rdata.As<std::int64_t>()[0] = First;
        rdata.As<std::int64_t>()[1] = Last;
        h5::write_dataset(rng, " data", rdata);
    }
    {
        NDArray off(DType::Int64, {rOffsets.size()});
        std::copy(rOffsets.begin(), rOffsets.end(), off.As<std::int64_t>());
        h5::Hid g = cgns_create_group(sect, "ElementStartOffset");
        cgns_write_node_attrs(g, "ElementStartOffset", "DataArray_t",
                              cgns_type_code(DType::Int64));
        h5::write_dataset(g, " data", off, GzipLevel);
    }
    {
        NDArray conn(DType::Int64, {rData.size()});
        std::copy(rData.begin(), rData.end(), conn.As<std::int64_t>());
        h5::Hid g = cgns_create_group(sect, "ElementConnectivity");
        cgns_write_node_attrs(g, "ElementConnectivity", "DataArray_t",
                              cgns_type_code(DType::Int64));
        h5::write_dataset(g, " data", conn, GzipLevel);
    }
}

// Whether the root has a child group whose "label" attribute is
// "CGNSBase_t" -- the structural (not extension/version-based) discriminator
// between the spec layout and the pre-v9.8.0 legacy one.
bool cgns_is_spec_layout(hid_t f) {
    for (const std::string& name : h5::group_links(f)) {
        if (cgns_is_reserved_name(name))
            continue;
        h5::Hid g = h5::open_group(f, name);
        if (h5::has_attr(g, "label") && h5::read_attr_string(g, "label") == "CGNSBase_t")
            return true;
    }
    return false;
}

// --- FlowSolution_t: point/cell data ----------------------------------------
//
// CGNS has NO component concept: one `DataArray_t` is one scalar, and there is
// no `NumberOfComponents` anywhere in the SIDS. A k>1 meshio++ array is
// therefore split into k sibling `DataArray_t` nodes suffixed `_0.._k-1` and
// re-joined on read from a contiguous run -- a documented meshio++ convention
// (see doc/formats/cgns.md), not something SIDS specifies.
constexpr const char* kCgnsVertexSolution = "FlowSolution";
constexpr const char* kCgnsCellSolution = "FlowSolutionCells";

/// `<base>_<i>` for a component of a multi-component array, or `<base>` when
/// the array is scalar. The single owner of the suffix convention.
std::string cgns_component_name(const std::string& rBase, std::size_t k, std::size_t i) {
    return k <= 1 ? rBase : rBase + "_" + std::to_string(i);
}

/// Write one `DataArray_t` holding component `i` of `rArr` (stride `k`).
void cgns_write_solution_array(hid_t sol, const std::string& rName, const NDArray& rArr,
                               std::size_t rows, std::size_t k, std::size_t i, int gzip_level) {
    NDArray col(DType::Float64, {rows});
    double* dst = col.As<double>();
    for (std::size_t r = 0; r < rows; ++r)
        dst[r] = detail::read_double(rArr, r * k + i);
    h5::Hid g = cgns_create_group(sol, rName);
    cgns_write_node_attrs(g, rName, "DataArray_t", "R8");
    h5::write_dataset(g, " data", col, gzip_level);
}

/// The `FlowSolution_t` node plus its `GridLocation_t` child, or an invalid
/// handle when the caller has nothing to write.
h5::Hid cgns_create_solution(hid_t zone, const std::string& rName, const std::string& rLocation) {
    h5::Hid sol = cgns_create_group(zone, rName);
    cgns_write_node_attrs(sol, rName, "FlowSolution_t", "MT");
    h5::Hid gl = cgns_create_group(sol, "GridLocation");
    cgns_write_node_attrs(gl, "GridLocation", "GridLocation_t", "C1");
    // Like ZoneType's payload: the declared length carries the string, with no
    // trailing NUL.
    h5::write_dataset(gl, " data", cgns_padded_int8(rLocation, rLocation.size()));
    return sol;
}

/// Re-join a zone's `DataArray_t` children into meshio++ arrays, undoing the
/// `_0.._k-1` split. Returns `name -> (rows, k)` interleaved Float64 arrays.
/// `rExpectedRows` is `NVertex`/`NCell`; a mismatch is a `ReadError`.
std::vector<std::pair<std::string, NDArray>> cgns_read_solution(hid_t sol,
                                                                const std::string& rSolName,
                                                                std::size_t expected_rows) {
    // Collect the raw arrays in name order first; `group_links` is name-sorted,
    // so a `_0.._k-1` run arrives contiguous and in component order for k < 10.
    // For k >= 10 the lexicographic order would interleave, so components are
    // placed by their parsed index rather than by arrival.
    std::vector<std::string> raw_names;
    for (const std::string& child : h5::group_links(sol)) {
        if (cgns_is_reserved_name(child) || child == "GridLocation")
            continue;
        h5::Hid g = h5::open_group(sol, child);
        if (!(h5::has_attr(g, "label") && h5::read_attr_string(g, "label") == "DataArray_t"))
            continue;
        raw_names.push_back(child);
    }

    // Group by base name: a trailing `_<digits>` marks a component of a
    // multi-component array, anything else is a scalar in its own right.
    std::map<std::string, std::map<std::size_t, std::string>> components;
    std::vector<std::string> base_order;
    for (const std::string& n : raw_names) {
        std::string base = n;
        std::size_t idx = 0;
        const std::size_t us = n.rfind('_');
        if (us != std::string::npos && us + 1 < n.size() &&
            n.find_first_not_of("0123456789", us + 1) == std::string::npos) {
            base = n.substr(0, us);
            idx = static_cast<std::size_t>(std::stoull(n.substr(us + 1)));
        }
        if (!components.count(base))
            base_order.push_back(base);
        components[base][idx] = n;
    }

    std::vector<std::pair<std::string, NDArray>> out;
    for (const std::string& base : base_order) {
        const std::map<std::size_t, std::string>& parts = components[base];
        // A genuine multi-component array is a contiguous 0..k-1 run. Anything
        // else (a lone `foo_7`, a gap) is treated as a scalar under its own
        // literal name -- guessing would invent components.
        bool contiguous = true;
        std::size_t expect = 0;
        for (const auto& kv : parts) {
            if (kv.first != expect++)
                contiguous = false;
        }
        const std::size_t k = parts.size();
        if (!contiguous || (k == 1 && parts.begin()->first != 0)) {
            for (const auto& kv : parts) {
                h5::Hid g = h5::open_group(sol, kv.second);
                NDArray a = h5::read_dataset(g, " data");
                if (a.Size() != expected_rows)
                    throw ReadError(detail::format_compat(
                        "CGNS: '{}/{}' has {} values but the zone has {} entities at this "
                        "GridLocation",
                        rSolName, kv.second, a.Size(), expected_rows));
                out.emplace_back(kv.second, std::move(a));
            }
            continue;
        }

        std::vector<NDArray> cols;
        cols.reserve(k);
        for (const auto& kv : parts) {
            h5::Hid g = h5::open_group(sol, kv.second);
            NDArray a = h5::read_dataset(g, " data");
            if (a.Size() != expected_rows)
                throw ReadError(detail::format_compat(
                    "CGNS: '{}/{}' has {} values but the zone has {} entities at this "
                    "GridLocation",
                    rSolName, kv.second, a.Size(), expected_rows));
            cols.push_back(std::move(a));
        }
        NDArray joined(DType::Float64, k > 1 ? std::vector<std::size_t>{expected_rows, k}
                                             : std::vector<std::size_t>{expected_rows});
        double* dst = joined.As<double>();
        for (std::size_t i = 0; i < k; ++i)
            for (std::size_t r = 0; r < expected_rows; ++r)
                dst[r * k + i] = detail::read_double(cols[i], r);
        out.emplace_back(base, std::move(joined));
    }
    return out;
}

}  // namespace

void write_cgns(const std::string& rPath, const Mesh& rMesh, int gzip_level) {
    h5::SilenceErrors silence;

    const std::size_t point_dim = rMesh.PointDim();
    if (point_dim > 3)
        throw WriteError(
            detail::format_compat("CGNS: PhysicalDimension must be 1..3, got {}", point_dim));

    // Validate every cell type up front (before any file is created) and
    // compute CellDim = the max topological dimension over all blocks.
    int cell_dim = 0;
    // Blocks that go out as NGON_n (2D: their cells ARE faces) and as
    // NGON_n + NFACE_n (3D: their cells are lists of signed face ids).
    // Classified by TYPE NAME, not `IsRagged()`: a *uniform* polygon block
    // (every cell the same node count) stores rectangularly and so is not
    // structurally ragged, yet still has no fixed-size ElementType_t.
    std::vector<std::size_t> ngon_blocks, nface_blocks;
    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi) {
        const auto cb = rMesh.Cells(bi);
        const std::string ctype(cb.Type());
        if (cgns_is_polygon_type(ctype)) {
            ngon_blocks.push_back(bi);
            if (cell_dim < 2)
                cell_dim = 2;
            continue;
        }
        if (cgns_is_polyhedron_type(ctype)) {
            nface_blocks.push_back(bi);
            cell_dim = 3;
            continue;
        }
        if (cb.IsRagged())
            throw WriteError(detail::format_compat(
                "CGNS: cell type '{}' is a ragged block with no CGNS representation (only "
                "polygon*/polyhedron* map to NGON_n/NFACE_n)",
                ctype));
        auto it = cgns_type_table().find(ctype);
        if (it == cgns_type_table().end()) {
            const CellType ct = cell_type_from_name(ctype);
            if (ct != CellType::Custom)
                throw WriteError(detail::format_compat(
                    "CGNS: cell type '{}' maps to a CGNS ElementType_t but its CGNS node "
                    "ordering is not yet verified in meshio++; refusing to write a guessed "
                    "ordering",
                    ctype));
            throw WriteError(detail::format_compat(
                "CGNS: cell type '{}' has no fixed-size CGNS ElementType_t equivalent", ctype));
        }
        const int d = cell_type_dimension(cell_type_from_name(ctype));
        if (d > cell_dim)
            cell_dim = d;
    }

    int phys_dim = static_cast<int>(std::max(point_dim, static_cast<std::size_t>(cell_dim)));
    phys_dim = std::clamp(phys_dim, 1, 3);
    const int cgns_cell_dim = cell_dim > 0 ? cell_dim : phys_dim;

    // NCell counts the zone's cells -- the blocks at CellDim. The face-based
    // families must be counted explicitly: `cell_type_from_name` does not know
    // `polygon<N>`/`polyhedron<N>`, so they would report dimension 0 and NCell
    // would come out too small. That is not cosmetic: cgnscheck sizes its cell
    // arrays from NCell and then reads NFACE_n, so an undercount corrupts its
    // heap rather than producing a diagnostic (found exactly that way).
    std::size_t n_cells_at_dim = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::string ctype(cb.Type());
        int d;
        if (cgns_is_polygon_type(ctype))
            d = 2;
        else if (cgns_is_polyhedron_type(ctype))
            d = 3;
        else
            d = cell_type_dimension(cell_type_from_name(ctype));
        if (d == cgns_cell_dim)
            n_cells_at_dim += cb.NumCells();
    }

    h5::Hid f = cgns_create_file(rPath);
    cgns_write_attr_str(f, "name", "HDF5 MotherNode", 33);
    cgns_write_attr_str(f, "label", "Root Node of HDF5 File", 33);
    cgns_write_attr_str(f, "type", "MT", 3);
    h5::write_dataset(f, " format", cgns_padded_int8("IEEE_LITTLE_32", 15));
    h5::write_dataset(f, " hdf5version", cgns_padded_int8(cgns_hdf5_version_string(), 33));

    h5::Hid cgver = cgns_create_group(f, "CGNSLibraryVersion");
    cgns_write_node_attrs(cgver, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4");
    NDArray ver(DType::Float32, {1});
    // 4.0 is not cosmetic when a face-based section is present: BELOW it,
    // NGON_n/NFACE_n use the 3.x layout, where each row is prefixed inline by
    // its own length instead of being described by ElementStartOffset. cgnslib
    // switches on this number, so declaring 3.1 while writing 4.0 arrays makes
    // its reader splice offsets into the connectivity and then corrupt its own
    // heap -- found with cgnscheck, which aborted rather than diagnosing it.
    // Files with no such section read identically under either number, so they
    // keep 3.1 and their bytes are unchanged.
    ver.As<float>()[0] = (ngon_blocks.empty() && nface_blocks.empty()) ? 3.1f : 4.0f;
    h5::write_dataset(cgver, " data", ver);

    h5::Hid base = cgns_create_group(f, "Base");
    cgns_write_node_attrs(base, "Base", "CGNSBase_t", "I4");
    {
        NDArray basedata(DType::Int32, {2});
        basedata.As<std::int32_t>()[0] = cgns_cell_dim;
        basedata.As<std::int32_t>()[1] = phys_dim;
        h5::write_dataset(base, " data", basedata);
    }

    const std::size_t n_points = rMesh.NumPoints();
    const bool wide = n_points > static_cast<std::size_t>(INT32_MAX) ||
                      n_cells_at_dim > static_cast<std::size_t>(INT32_MAX) ||
                      !ngon_blocks.empty() || !nface_blocks.empty();
    const DType zone_dt = wide ? DType::Int64 : DType::Int32;

    h5::Hid zone = cgns_create_group(base, "Zone1");
    cgns_write_node_attrs(zone, "Zone1", "Zone_t", cgns_type_code(zone_dt));
    {
        NDArray zdata(zone_dt, {3, 1});
        detail::dispatch_dtype(zone_dt, [&]<class T>() {
            T* d = zdata.As<T>();
            d[0] = static_cast<T>(n_points);
            d[1] = static_cast<T>(n_cells_at_dim);
            d[2] = static_cast<T>(0);
        });
        h5::write_dataset(zone, " data", zdata);
    }

    h5::Hid zt = cgns_create_group(zone, "ZoneType");
    cgns_write_node_attrs(zt, "ZoneType", "ZoneType_t", "C1");
    h5::write_dataset(zt, " data", cgns_padded_int8("Unstructured", 12));

    h5::Hid coords = cgns_create_group(zone, "GridCoordinates");
    cgns_write_node_attrs(coords, "GridCoordinates", "GridCoordinates_t", "MT");
    {
        const NDArray& points = rMesh.Points();
        const DType coord_dt = points.Dtype() == DType::Float32 ? DType::Float32 : DType::Float64;
        const char* names[3] = {"CoordinateX", "CoordinateY", "CoordinateZ"};
        const std::size_t n_coords = phys_dim >= 3 ? 3 : 2;
        for (std::size_t c = 0; c < n_coords; ++c) {
            h5::Hid g = cgns_create_group(coords, names[c]);
            NDArray col(coord_dt, {n_points});
            detail::dispatch_dtype(coord_dt, [&]<class T>() {
                T* dst = col.As<T>();
                parallel_for_bw(n_points, [&](std::size_t i) {
                    dst[i] = static_cast<T>(
                        (c < point_dim) ? detail::read_double(points, i * point_dim + c) : 0.0);
                });
            });
            cgns_write_node_attrs(g, names[c], "DataArray_t", cgns_type_code(coord_dt));
            h5::write_dataset(g, " data", col, gzip_level);
        }
    }

    // Elements: one Elements_t section per meshio++ cell block, in mesh
    // order, with contiguous non-overlapping 1-based ElementRanges. Unlike
    // MED, sections of the same type are NOT consolidated -- CGNS has no
    // "two sections, one type" restriction.
    std::int64_t next = 1;
    std::unordered_map<std::string, int> type_counts;

    // ---- the face-based sections come first ---------------------------------
    // NFACE_n references faces by ELEMENT id, so every face must already have
    // one: the NGON_n section is emitted ahead of the NFACE_n sections that
    // point into it. Faces are deduplicated across the polyhedral blocks ONLY
    // (see build_global_faces' block-filtered overload) -- a hexahedron in the
    // same mesh keeps its own HEXA_8 section, and putting its faces here would
    // leave NGON_n elements no cell ever references.
    std::vector<std::int64_t> face_element_id;  // GlobalFaces id -> 1-based element id
    if (!nface_blocks.empty()) {
        const detail::GlobalFaces gf = detail::build_global_faces(rMesh, nface_blocks);
        if (gf.mNumNonManifold > 0)
            log::warn("CGNS: {} face(s) are shared by three or more cells; NFACE_n still "
                      "references them correctly, but the mesh is non-manifold",
                      gf.mNumNonManifold);

        std::vector<std::int64_t> off{0}, data;
        off.reserve(gf.NumFaces() + 1);
        for (std::size_t f = 0; f < gf.NumFaces(); ++f) {
            for (std::size_t j = 0; j < gf.FaceSize(f); ++j)
                data.push_back(gf.Face(f)[j] + 1);  // CGNS node ids are 1-based
            off.push_back(static_cast<std::int64_t>(data.size()));
        }
        const std::int64_t first = next;
        const std::int64_t last = next + static_cast<std::int64_t>(gf.NumFaces()) - 1;
        next = last + 1;
        face_element_id.resize(gf.NumFaces());
        for (std::size_t f = 0; f < gf.NumFaces(); ++f)
            face_element_id[f] = first + static_cast<std::int64_t>(f);
        cgns_write_poly_section(zone, "NGON_n_1", kCgnsNgon, first, last, off, data, gzip_level);

        // One NFACE_n per polyhedral block, so the reader can rebuild the
        // blocks rather than one merged one.
        std::size_t compact = 0;
        int nface_idx = 0;
        for (std::size_t bi : nface_blocks) {
            const std::size_t nc = rMesh.Cells(bi).NumCells();
            if (nc == 0) {
                compact += nc;
                continue;
            }
            std::vector<std::int64_t> coff{0}, cdata;
            for (std::size_t c = compact; c < compact + nc; ++c) {
                for (std::size_t j = 0; j < gf.NumCellFaces(c); ++j) {
                    const std::int64_t sid = gf.CellFaces(c)[j];
                    const std::int64_t fid = face_element_id[static_cast<std::size_t>(
                        (sid > 0 ? sid : -sid) - 1)];
                    cdata.push_back(sid > 0 ? fid : -fid);  // sign = "traverse reversed"
                }
                coff.push_back(static_cast<std::int64_t>(cdata.size()));
            }
            compact += nc;
            const std::int64_t cf = next;
            const std::int64_t cl = next + static_cast<std::int64_t>(nc) - 1;
            next = cl + 1;
            cgns_write_poly_section(zone, detail::format_compat("NFACE_n_{}", ++nface_idx),
                                    kCgnsNface, cf, cl, coff, cdata, gzip_level);
        }
    }
    // A 2D jagged block is itself a face list, so it needs no dedup at all.
    {
        int idx = nface_blocks.empty() ? 0 : 1;
        for (std::size_t bi : ngon_blocks) {
            const auto cb = rMesh.Cells(bi);
            const std::size_t nc = cb.NumCells();
            if (nc == 0)
                continue;
            std::vector<std::int64_t> off{0}, data;
            for (std::size_t i = 0; i < nc; ++i) {
                if (cb.IsRagged()) {
                    const std::int64_t* row = cb.Row(i);
                    for (std::size_t j = 0; j < cb.RowSize(i); ++j)
                        data.push_back(row[j] + 1);
                } else {
                    const std::size_t npc = cb.NodesPerCell();
                    for (std::size_t j = 0; j < npc; ++j)
                        data.push_back(detail::read_int(cb.Conn(), i * npc + j) + 1);
                }
                off.push_back(static_cast<std::int64_t>(data.size()));
            }
            const std::int64_t first = next;
            const std::int64_t last = next + static_cast<std::int64_t>(nc) - 1;
            next = last + 1;
            cgns_write_poly_section(zone, detail::format_compat("NGON_n_{}", ++idx), kCgnsNgon,
                                    first, last, off, data, gzip_level);
        }
    }

    for (std::size_t k = 0; k < rMesh.NumCellBlocks(); ++k) {
        const auto cb = rMesh.Cells(k);
        const std::string ctype(cb.Type());
        if (cgns_is_polygon_type(ctype) || cgns_is_polyhedron_type(ctype))
            continue;  // already emitted above as NGON_n / NFACE_n
        const CgnsTypeInfo& info = cgns_type_table().at(ctype);
        const std::size_t npc = cb.NodesPerCell();
        const std::size_t nc = cb.NumCells();
        if (nc == 0)
            continue;  // a zero-length ElementRange is not representable

        const std::int64_t first = next;
        const std::int64_t last = next + static_cast<std::int64_t>(nc) - 1;
        next = last + 1;

        const int idx = ++type_counts[info.mCgnsName];
        const std::string section_name = detail::format_compat("{}_{}", info.mCgnsName, idx);

        h5::Hid sect = cgns_create_group(zone, section_name);
        cgns_write_node_attrs(sect, section_name, "Elements_t", "I4");
        {
            NDArray sdata(DType::Int32, {2});
            sdata.As<std::int32_t>()[0] = info.mCode;
            sdata.As<std::int32_t>()[1] = 0;  // ElementSizeBoundary: unsorted
            h5::write_dataset(sect, " data", sdata);
        }

        const bool section_wide =
            (cb.Conn().Dtype() != DType::Int32) || (n_points > static_cast<std::size_t>(INT32_MAX));
        const DType out_dt = section_wide ? DType::Int64 : DType::Int32;

        h5::Hid rng = cgns_create_group(sect, "ElementRange");
        cgns_write_node_attrs(rng, "ElementRange", "IndexRange_t", cgns_type_code(out_dt));
        {
            NDArray rdata(out_dt, {2});
            detail::dispatch_dtype(out_dt, [&]<class T>() {
                T* d = rdata.As<T>();
                d[0] = static_cast<T>(first);
                d[1] = static_cast<T>(last);
            });
            h5::write_dataset(rng, " data", rdata);
        }

        const std::vector<int>* perm = info.mPerm.empty() ? nullptr : &info.mPerm;
        NDArray conn = cgns_permute_conn(cb.Conn(), nc, npc, +1, perm);
        // Widen to the section's chosen dtype if it disagrees with Conn()'s.
        if (conn.Dtype() != out_dt) {
            NDArray widened(out_dt, {nc, npc});
            detail::dispatch_dtype(out_dt, [&]<class T>() {
                T* dst = widened.As<T>();
                for (std::size_t i = 0; i < nc * npc; ++i)
                    dst[i] = static_cast<T>(detail::read_int(conn, i));
            });
            conn = std::move(widened);
        }
        conn.Reshape({nc * npc});

        h5::Hid ec = cgns_create_group(sect, "ElementConnectivity");
        cgns_write_node_attrs(ec, "ElementConnectivity", "DataArray_t", cgns_type_code(out_dt));
        h5::write_dataset(ec, " data", conn, gzip_level);
    }

    // FlowSolution_t: point_data at "Vertex", cell_data at "CellCenter".
    // Names come from the sorted *DataNames() accessors, so output order is
    // deterministic. field_data has no CGNS home (it is neither per-vertex nor
    // per-cell) and is not written.
    {
        std::vector<std::string> pnames = rMesh.PointDataNames();
        if (!pnames.empty()) {
            h5::Hid sol = cgns_create_solution(zone, kCgnsVertexSolution, "Vertex");
            for (const std::string& name : pnames) {
                const NDArray& a = rMesh.PointData(name);
                const std::size_t rows = detail::rows(a);
                const std::size_t k = detail::cols(a);
                if (rows != n_points) {
                    log::warn(
                        "CGNS: point_data '{}' has {} rows but the mesh has {} points; not "
                        "written.",
                        name, rows, n_points);
                    continue;
                }
                for (std::size_t i = 0; i < k; ++i)
                    cgns_write_solution_array(sol, cgns_component_name(name, k, i), a, rows, k, i,
                                              gzip_level);
            }
        }
    }
    {
        // A zone-wide CellCenter array has one value per zone cell, but
        // meshio++'s cell_data is per BLOCK -- and only blocks at CellDim are
        // zone cells. Concatenating block-major is therefore only well defined
        // when every block is at CellDim; a mixed-dimension mesh (e.g. tets
        // plus boundary triangles) has no way to distribute the array back
        // across blocks on read without inventing values, so it is skipped
        // with a warning rather than written wrongly.
        bool all_at_cell_dim = rMesh.NumCellBlocks() > 0;
        for (const auto cb : rMesh.CellRange())
            if (cell_type_dimension(cell_type_from_name(std::string(cb.Type()))) != cgns_cell_dim)
                all_at_cell_dim = false;

        std::vector<std::string> cnames = rMesh.CellDataNames();
        if (!cnames.empty() && !all_at_cell_dim) {
            log::warn(
                "CGNS: this mesh mixes cell blocks of different topological dimensions, so a "
                "zone-wide CellCenter FlowSolution cannot be distributed back across them; {} "
                "cell_data array(s) not written.",
                cnames.size());
            detail::provenance_note(
                "data-dropped",
                std::to_string(cnames.size()) +
                    " cell_data array(s) not written -- CGNS's CellCenter FlowSolution is "
                    "per-zone and this mesh mixes topological dimensions");
        } else if (!cnames.empty()) {
            h5::Hid sol = cgns_create_solution(zone, kCgnsCellSolution, "CellCenter");
            for (const std::string& name : cnames) {
                if (rMesh.CellDataNumBlocks(name) != rMesh.NumCellBlocks()) {
                    log::warn("CGNS: cell_data '{}' covers {} of {} blocks; not written.", name,
                              rMesh.CellDataNumBlocks(name), rMesh.NumCellBlocks());
                    continue;
                }
                const std::size_t k = detail::cols(rMesh.CellData(name, 0));
                std::size_t rows = 0;
                bool consistent = true;
                for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
                    const NDArray& d = rMesh.CellData(name, b);
                    if (detail::cols(d) != k || detail::rows(d) != rMesh.Cells(b).NumCells())
                        consistent = false;
                    rows += detail::rows(d);
                }
                if (!consistent || rows != n_cells_at_dim) {
                    log::warn(
                        "CGNS: cell_data '{}' does not line up with the zone's cells; not "
                        "written.",
                        name);
                    continue;
                }
                // Concatenate the blocks block-major into one interleaved
                // (n_cells, k) buffer, then split it per component.
                NDArray merged(DType::Float64, k > 1 ? std::vector<std::size_t>{rows, k}
                                                     : std::vector<std::size_t>{rows});
                double* dst = merged.As<double>();
                std::size_t row = 0;
                for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
                    const NDArray& d = rMesh.CellData(name, b);
                    const std::size_t nb = detail::rows(d);
                    for (std::size_t r = 0; r < nb; ++r)
                        for (std::size_t i = 0; i < k; ++i)
                            dst[(row + r) * k + i] = detail::read_double(d, r * k + i);
                    row += nb;
                }
                for (std::size_t i = 0; i < k; ++i)
                    cgns_write_solution_array(sol, cgns_component_name(name, k, i), merged, rows, k,
                                              i, gzip_level);
            }
        }
    }
}

Mesh read_cgns(const std::string& rPath) {
#ifdef MESHIOPLUSPLUS_HAS_CGNSLIB
    // With cgnslib built, IT is the reader: the input is not ours, and the MLL
    // reaches things this raw-HDF5 path fundamentally cannot -- the ADF
    // container, links, multiple bases, and NGON_n/NFACE_n polyhedral sections.
    //
    // The pre-v9.8.0 legacy layout has no ADF node attributes at all, so the
    // MLL rejects it; that ONE case falls through to the hand-rolled path
    // below. This is a narrow, specific fallback and not a blanket catch: a
    // genuine MLL error still surfaces, or a corrupt file would be silently
    // re-read by a reader that cannot diagnose it either.
    try {
        return read_cgns_mll(rPath);
    } catch (const ReadError&) {
        // fall through to the structural probe, which either reads the legacy
        // layout or reports its own (more specific) error
    }
#endif
    h5::SilenceErrors silence;
    h5::Hid f = h5::open_file_read(rPath);

    if (!cgns_is_spec_layout(f))
        return cgns_read_legacy(f);

    // Locate the first CGNSBase_t child of the root (structural, not by
    // name); warn (never silently drop data) if there is more than one.
    std::string base_name;
    int n_bases = 0;
    for (const std::string& name : h5::group_links(f)) {
        if (cgns_is_reserved_name(name))
            continue;
        h5::Hid g = h5::open_group(f, name);
        if (h5::has_attr(g, "label") && h5::read_attr_string(g, "label") == "CGNSBase_t") {
            if (n_bases == 0)
                base_name = name;
            ++n_bases;
        }
    }
    if (n_bases > 1)
        log::warn("CGNS: file has {} CGNSBase_t nodes; only the first ('{}') is read.", n_bases,
                  base_name);

    h5::Hid base = h5::open_group(f, base_name);

    Mesh mesh;
    NDArray points;  // accumulated across zones
    std::vector<NDArray> point_chunks;
    std::int64_t point_offset = 0;
    std::size_t point_dim_out = 3;

    // FlowSolution_t is only read for a single-zone file: across several zones
    // the point/cell arrays would have to be concatenated in exactly the order
    // the zones happen to be listed in, and a solution present on only some
    // zones has no defensible filler. Our own writer emits one zone.
    std::size_t n_zones = 0;
    for (const std::string& zname : h5::group_links(base)) {
        if (cgns_is_reserved_name(zname))
            continue;
        h5::Hid z = h5::open_group(base, zname);
        if (h5::has_attr(z, "label") && h5::read_attr_string(z, "label") == "Zone_t")
            ++n_zones;
    }

    for (const std::string& zname : h5::group_links(base)) {
        if (cgns_is_reserved_name(zname))
            continue;  // Base's own " data" payload ([CellDim, PhysDim])
        h5::Hid zone = h5::open_group(base, zname);
        if (!(h5::has_attr(zone, "label") && h5::read_attr_string(zone, "label") == "Zone_t"))
            continue;

        if (h5::exists(zone, "ZoneType")) {
            h5::Hid zt = h5::open_group(zone, "ZoneType");
            if (h5::exists(zt, " data")) {
                NDArray raw = h5::read_dataset(zt, " data");
                std::string zt_name(raw.Size(), '\0');
                for (std::size_t i = 0; i < raw.Size(); ++i)
                    zt_name[i] = static_cast<char>(detail::read_int(raw, i));
                if (zt_name != "Unstructured")
                    throw ReadError(detail::format_compat(
                        "CGNS: zone '{}' has ZoneType '{}'; only Unstructured zones are "
                        "supported.",
                        zname, zt_name));
            }
        }

        h5::Hid coords = h5::open_group(zone, "GridCoordinates");
        std::vector<std::string> axes;
        for (const char* ax : {"CoordinateX", "CoordinateY", "CoordinateZ"})
            if (h5::exists(coords, ax))
                axes.push_back(ax);
        if (axes.empty())
            throw ReadError(detail::format_compat("CGNS: zone '{}' has no GridCoordinates", zname));
        point_dim_out = std::max<std::size_t>(2, axes.size());

        std::vector<NDArray> cols;
        std::size_t n_zone_points = 0;
        for (const std::string& ax : axes) {
            h5::Hid g = h5::open_group(coords, ax);
            NDArray c = h5::read_dataset(g, " data");
            n_zone_points = c.Shape().empty() ? 0 : c.Shape()[0];
            cols.push_back(std::move(c));
        }
        NDArray zpts(DType::Float64, {n_zone_points, point_dim_out});
        double* pp = zpts.As<double>();
        for (std::size_t i = 0; i < n_zone_points; ++i)
            for (std::size_t d = 0; d < point_dim_out; ++d)
                pp[i * point_dim_out + d] = d < cols.size() ? detail::read_double(cols[d], i) : 0.0;
        point_chunks.push_back(std::move(zpts));

        // Elements_t sections found by label (never by name), ordered
        // ascending by ElementRange[0] -- reproduces write_cgns's own block
        // order exactly and is immune to HDF5 iteration order.
        struct Sect {
            std::string mName;
            std::int64_t mFirst;
        };
        // Cell counts of the blocks this zone contributes, in emission order --
        // what a zone-wide CellCenter FlowSolution is split back across.
        std::vector<std::size_t> zone_block_cells;

        std::vector<Sect> sects;
        for (const std::string& sname : h5::group_links(zone)) {
            if (cgns_is_reserved_name(sname))
                continue;  // Zone1's own " data" payload ([NVertex, NCell, 0])
            h5::Hid s = h5::open_group(zone, sname);
            if (!(h5::has_attr(s, "label") && h5::read_attr_string(s, "label") == "Elements_t"))
                continue;
            h5::Hid rng = h5::open_group(s, "ElementRange");
            NDArray range = h5::read_dataset(rng, " data");
            if (range.Size() < 2)
                throw ReadError(detail::format_compat(
                    "CGNS: section '{}' has a malformed ElementRange", sname));
            sects.push_back({sname, detail::read_int(range, 0)});
        }
        std::sort(sects.begin(), sects.end(),
                  [](const Sect& a, const Sect& b) { return a.mFirst < b.mFirst; });

        // --- pre-pass: the face-based (NGON_n / NFACE_n) sections ------------
        // Read ahead of the main loop for two reasons. An NFACE_n may legally
        // precede the NGON_n it points into, and -- the load-bearing one -- an
        // NGON_n becomes a `polygon` cell block ONLY when no NFACE_n references
        // it: otherwise it is the shared face pool, not a set of cells, and
        // emitting it as one would double every polyhedron's geometry.
        std::map<std::int64_t, CgnsPolySection> poly;  // keyed by ElementRange start
        bool any_nface = false;
        for (const Sect& sec : sects) {
            h5::Hid s = h5::open_group(zone, sec.mName);
            NDArray sdata = h5::read_dataset(s, " data");
            if (sdata.Size() < 1)
                continue;
            const int code = static_cast<int>(detail::read_int(sdata, 0));
            if (code != kCgnsNgon && code != kCgnsNface)
                continue;
            poly.emplace(sec.mFirst, cgns_read_poly_section(s, sec.mName, code));
            any_nface = any_nface || code == kCgnsNface;
        }
        // Which NGON element ids are referenced by some NFACE_n cell.
        std::set<std::int64_t> ngon_used;
        if (any_nface) {
            for (const auto& kv : poly) {
                if (kv.second.mCode != kCgnsNface)
                    continue;
                for (std::int64_t v : kv.second.mData)
                    ngon_used.insert(v < 0 ? -v : v);
            }
        }

        for (const Sect& sec : sects) {
            h5::Hid s = h5::open_group(zone, sec.mName);
            NDArray sdata = h5::read_dataset(s, " data");
            if (sdata.Size() < 1)
                throw ReadError(detail::format_compat(
                    "CGNS: section '{}' has a malformed Elements_t descriptor", sec.mName));
            const int code = static_cast<int>(detail::read_int(sdata, 0));

            if (code == kCgnsNgon || code == kCgnsNface) {
                const CgnsPolySection& ps = poly.at(sec.mFirst);
                const std::size_t nc = ps.Count();
                if (nc == 0)
                    continue;

                if (code == kCgnsNgon) {
                    // A face pool referenced by some NFACE_n is not a set of
                    // cells; emitting it as one would duplicate every
                    // polyhedron's geometry.
                    bool referenced = false;
                    for (std::int64_t e = ps.mFirst; e <= ps.mLast && !referenced; ++e)
                        referenced = ngon_used.count(e) != 0;
                    if (referenced)
                        continue;

                    std::vector<std::vector<std::int64_t>> rows(nc);
                    for (std::size_t i = 0; i < nc; ++i) {
                        const std::size_t lo = static_cast<std::size_t>(ps.mOffsets[i]);
                        const std::size_t hi = static_cast<std::size_t>(ps.mOffsets[i + 1]);
                        rows[i].reserve(hi - lo);
                        for (std::size_t j = lo; j < hi; ++j)
                            rows[i].push_back(ps.mData[j] - 1 + point_offset);
                    }
                    // Grouped by node count, the convention the OpenFOAM and
                    // EnSight readers already use.
                    std::map<std::size_t, std::vector<std::vector<std::int64_t>>> by_n;
                    for (auto& r : rows)
                        by_n[r.size()].push_back(std::move(r));
                    for (auto& kv : by_n) {
                        const std::size_t cnt = kv.second.size();
                        mesh.AddPolygonBlock("polygon" + std::to_string(kv.first),
                                             std::move(kv.second));
                        zone_block_cells.push_back(cnt);
                    }
                    continue;
                }

                // NFACE_n: dereference each signed face id into its node ring,
                // reversing where the sign says so.
                std::vector<std::vector<std::vector<std::int64_t>>> cells(nc);
                for (std::size_t i = 0; i < nc; ++i) {
                    const std::size_t lo = static_cast<std::size_t>(ps.mOffsets[i]);
                    const std::size_t hi = static_cast<std::size_t>(ps.mOffsets[i + 1]);
                    for (std::size_t j = lo; j < hi; ++j) {
                        const std::int64_t sid = ps.mData[j];
                        const std::int64_t fid = sid < 0 ? -sid : sid;
                        // Locate the NGON section holding element id `fid`.
                        const CgnsPolySection* src = nullptr;
                        for (const auto& kv : poly) {
                            if (kv.second.mCode == kCgnsNgon && fid >= kv.second.mFirst &&
                                fid <= kv.second.mLast) {
                                src = &kv.second;
                                break;
                            }
                        }
                        if (!src)
                            throw ReadError(detail::format_compat(
                                "CGNS: section '{}' references face element {}, which no NGON_n "
                                "section defines",
                                sec.mName, fid));
                        const std::size_t fi = static_cast<std::size_t>(fid - src->mFirst);
                        const std::size_t flo = static_cast<std::size_t>(src->mOffsets[fi]);
                        const std::size_t fhi = static_cast<std::size_t>(src->mOffsets[fi + 1]);
                        std::vector<std::int64_t> ring;
                        ring.reserve(fhi - flo);
                        for (std::size_t k = flo; k < fhi; ++k)
                            ring.push_back(src->mData[k] - 1 + point_offset);
                        if (sid < 0)
                            std::reverse(ring.begin(), ring.end());
                        cells[i].push_back(std::move(ring));
                    }
                }
                // Grouped by DISTINCT node count -> polyhedron<N>, matching the
                // OpenFOAM/EnSight/cgnslib readers.
                std::map<std::size_t, std::vector<std::vector<std::vector<std::int64_t>>>> by_n;
                for (auto& c : cells) {
                    std::set<std::int64_t> uniq;
                    for (const auto& f : c)
                        uniq.insert(f.begin(), f.end());
                    by_n[uniq.size()].push_back(std::move(c));
                }
                for (auto& kv : by_n) {
                    const std::size_t cnt = kv.second.size();
                    mesh.AddPolyhedronBlock("polyhedron" + std::to_string(kv.first),
                                            std::move(kv.second));
                    zone_block_cells.push_back(cnt);
                }
                continue;
            }
            if (code == 20) {
                const auto& names = cgns_code_to_name();
                auto nit = names.find(code);
                throw ReadError(detail::format_compat(
                    "CGNS: element section '{}' has ElementType {} ({}); MIXED sections are not "
                    "supported.",
                    sec.mName, nit != names.end() ? nit->second : "?", code));
            }
            const auto& code_map = cgns_code_to_meshio();
            auto tit = code_map.find(code);
            if (tit == code_map.end()) {
                const auto& names = cgns_code_to_name();
                auto nit = names.find(code);
                if (nit != names.end())
                    throw ReadError(detail::format_compat(
                        "CGNS: element section '{}' has type {} ({}), whose node ordering is "
                        "not yet verified in meshio++.",
                        sec.mName, nit->second, code));
                throw ReadError(detail::format_compat(
                    "CGNS: element section '{}' has unknown ElementType code {}.", sec.mName,
                    code));
            }
            const std::string& meshio_type = tit->second;
            const CgnsTypeInfo& info = cgns_type_table().at(meshio_type);

            h5::Hid rng = h5::open_group(s, "ElementRange");
            NDArray range = h5::read_dataset(rng, " data");
            const std::int64_t first = detail::read_int(range, 0);
            const std::int64_t last = detail::read_int(range, 1);
            if (last < first)
                throw ReadError(detail::format_compat(
                    "CGNS: section '{}' has an empty/overlapping ElementRange [{}, {}]", sec.mName,
                    first, last));
            const std::size_t nc = static_cast<std::size_t>(last - first + 1);
            const std::size_t npc =
                static_cast<std::size_t>(cell_type_num_nodes(cell_type_from_name(meshio_type)));

            h5::Hid conn = h5::open_group(s, "ElementConnectivity");
            NDArray flat = h5::read_dataset(conn, " data");
            if (flat.Size() != nc * npc)
                throw ReadError(detail::format_compat(
                    "CGNS: section '{}' declares {} elements of {} ({} nodes) but "
                    "ElementConnectivity has {} entries",
                    sec.mName, nc, info.mCgnsName, nc * npc, flat.Size()));

            const std::vector<int>* perm = info.mPerm.empty() ? nullptr : &info.mPerm;
            NDArray out = cgns_permute_conn(flat, nc, npc, -1, perm);
            if (point_offset != 0) {
                detail::dispatch_dtype(out.Dtype(), [&]<class T>() {
                    T* d = out.As<T>();
                    const T off = static_cast<T>(point_offset);
                    for (std::size_t i = 0; i < nc * npc; ++i)
                        d[i] = static_cast<T>(d[i] + off);
                });
            }
            mesh.AddCellBlock(meshio_type, std::move(out));
            zone_block_cells.push_back(nc);
        }

        // FlowSolution_t (see the writer): "Vertex" arrays become point_data,
        // "CellCenter" arrays become cell_data split back across this zone's
        // blocks in ElementRange order. GridLocation absent => "Vertex", the
        // SIDS default.
        if (n_zones == 1) {
            std::size_t zone_total_cells = 0;
            for (std::size_t c : zone_block_cells)
                zone_total_cells += c;

            for (const std::string& child : h5::group_links(zone)) {
                if (cgns_is_reserved_name(child))
                    continue;
                h5::Hid sol = h5::open_group(zone, child);
                if (!(h5::has_attr(sol, "label") &&
                      h5::read_attr_string(sol, "label") == "FlowSolution_t"))
                    continue;

                std::string location = "Vertex";
                if (h5::exists(sol, "GridLocation")) {
                    h5::Hid gl = h5::open_group(sol, "GridLocation");
                    if (h5::exists(gl, " data")) {
                        NDArray raw = h5::read_dataset(gl, " data");
                        std::string s(raw.Size(), '\0');
                        for (std::size_t i = 0; i < raw.Size(); ++i)
                            s[i] = static_cast<char>(detail::read_int(raw, i));
                        while (!s.empty() && (s.back() == '\0' || s.back() == ' '))
                            s.pop_back();
                        location = s;
                    }
                }

                if (location == "Vertex") {
                    for (auto& [name, arr] : cgns_read_solution(sol, child, n_zone_points))
                        mesh.AddPointData(name, std::move(arr));
                } else if (location == "CellCenter") {
                    for (auto& [name, arr] : cgns_read_solution(sol, child, zone_total_cells)) {
                        // Split the zone-wide array back across the blocks.
                        const std::size_t k = detail::cols(arr);
                        std::vector<NDArray> blocks;
                        blocks.reserve(zone_block_cells.size());
                        std::size_t row = 0;
                        for (std::size_t nb : zone_block_cells) {
                            NDArray b(DType::Float64, k > 1 ? std::vector<std::size_t>{nb, k}
                                                            : std::vector<std::size_t>{nb});
                            double* dst = b.As<double>();
                            for (std::size_t r = 0; r < nb; ++r)
                                for (std::size_t i = 0; i < k; ++i)
                                    dst[r * k + i] = detail::read_double(arr, (row + r) * k + i);
                            row += nb;
                            blocks.push_back(std::move(b));
                        }
                        mesh.AddCellData(name, std::move(blocks));
                    }
                } else {
                    log::warn(
                        "CGNS: FlowSolution '{}' has GridLocation '{}'; only Vertex and "
                        "CellCenter are supported, so it was not read.",
                        child, location);
                }
            }
        } else if (n_zones > 1) {
            for (const std::string& child : h5::group_links(zone)) {
                if (cgns_is_reserved_name(child))
                    continue;
                h5::Hid sol = h5::open_group(zone, child);
                if (h5::has_attr(sol, "label") &&
                    h5::read_attr_string(sol, "label") == "FlowSolution_t") {
                    log::warn(
                        "CGNS: file has {} zones; FlowSolution '{}' on zone '{}' was not read "
                        "(cross-zone field concatenation is not supported).",
                        n_zones, child, zname);
                    break;
                }
            }
        }

        point_offset += static_cast<std::int64_t>(n_zone_points);
    }

    if (point_chunks.empty())
        throw ReadError(
            detail::format_compat("CGNS: base '{}' has no Unstructured zones", base_name));
    if (point_chunks.size() == 1) {
        mesh.AssignPoints(std::move(point_chunks.front()));
    } else {
        std::size_t total = 0;
        for (const NDArray& c : point_chunks)
            total += c.Shape()[0];
        NDArray all(DType::Float64, {total, point_dim_out});
        double* dst = all.As<double>();
        std::size_t row = 0;
        for (const NDArray& c : point_chunks) {
            const double* src = c.As<double>();
            const std::size_t n = c.Shape()[0];
            std::memcpy(dst + row * point_dim_out, src, n * point_dim_out * sizeof(double));
            row += n;
        }
        mesh.AssignPoints(std::move(all));
    }
    return mesh;
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
