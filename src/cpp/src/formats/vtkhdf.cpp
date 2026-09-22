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
#include <filesystem>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/vtkhdf.hpp"
#include "meshioplusplus/formats/vtkhdf_time_series.hpp"
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/mesh_carve.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtk_cells.hpp"
#include "meshioplusplus/detail/vtk_xml.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/region.hpp"
#include "meshioplusplus/types.hpp"
#include "meshioplusplus/vtk_common.hpp"

namespace meshioplusplus {

namespace {

using detail::cols;
using detail::read_double;
using detail::read_int;
using h5::Hid;

constexpr int kVtkhdfKnownMinor = 8;
constexpr std::int64_t kVtkhdfPolyhedron = 42;
constexpr std::size_t kVtkhdfGzipMinBytes = 4096;
const char* const kVtkhdfRoot = "VTKHDF";
const char* const kVtkhdfProvenanceAttr = "meshioplusplus:provenance";
const char* const kVtkhdfCategories[4] = {"Vertices", "Lines", "Polygons", "Strips"};

using I64 = std::int64_t;
using I64Vec = std::vector<std::int64_t>;

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------
std::string vtkhdf_group_name(hid_t loc) {
    const ssize_t len = H5Iget_name(loc, nullptr, 0);
    if (len <= 0)
        return "/";
    std::string name(static_cast<std::size_t>(len), '\0');
    H5Iget_name(loc, name.data(), static_cast<std::size_t>(len) + 1);
    return name;
}

I64Vec vtkhdf_ints(hid_t g, const std::string& rName) {
    if (!h5::exists(g, rName))
        throw ReadError("meshio++: vtkhdf: '" + vtkhdf_group_name(g) + "' has no '" + rName +
                        "' dataset");
    return detail::vtu_to_int64(h5::read_dataset(g, rName));
}

I64Vec vtkhdf_int_rows(hid_t g, const std::string& rName, std::size_t Row0, std::size_t Count) {
    return detail::vtu_to_int64(h5::read_dataset_rows(g, rName, Row0, Count));
}

I64 vtkhdf_sum(const I64Vec& rV, std::size_t Lo, std::size_t Hi) {
    I64 s = 0;
    for (std::size_t i = Lo; i < Hi && i < rV.size(); ++i)
        s += rV[i];
    return s;
}

std::size_t vtkhdf_to_size(I64 V, const char* pWhat) {
    if (V < 0)
        throw ReadError(std::string("meshio++: vtkhdf: negative ") + pWhat);
    return static_cast<std::size_t>(V);
}

/// Exclusive prefix sums: where each piece starts inside a concatenation.
I64Vec vtkhdf_starts(const I64Vec& rCounts) {
    I64Vec out(rCounts.size(), 0);
    I64 acc = 0;
    for (std::size_t i = 0; i < rCounts.size(); ++i) {
        out[i] = acc;
        acc += rCounts[i];
    }
    return out;
}

bool vtkhdf_is_dataset(hid_t loc, const std::string& rName) {
    Hid d(H5Dopen2(loc, rName.c_str(), H5P_DEFAULT), H5Dclose);
    return d.Valid();
}

/// Names of the datasets directly under group `rGroup` of `loc`, in sorted order.
std::vector<std::string> vtkhdf_dataset_names(hid_t loc, const std::string& rGroup) {
    std::vector<std::string> out;
    if (!h5::exists(loc, rGroup))
        return out;
    Hid g = h5::open_group(loc, rGroup);
    for (const auto& name : h5::group_links(g))
        if (vtkhdf_is_dataset(g, name))
            out.push_back(name);
    return out;
}

/// The oldest VTKHDF version covering the content: {major, minor}.
std::pair<int, int> vtkhdf_needed_version(VtkhdfType Type, bool HasPoly) {
    if (HasPoly)
        return {2, 5};
    if (Type == VtkhdfType::PartitionedDataSetCollection || Type == VtkhdfType::MultiBlockDataSet)
        return {2, 1};
    if (Type == VtkhdfType::PolyData)
        return {2, 0};
    return {1, 0};
}

const char* vtkhdf_type_name(VtkhdfType Type) {
    switch (Type) {
        case VtkhdfType::UnstructuredGrid:
            return "UnstructuredGrid";
        case VtkhdfType::PolyData:
            return "PolyData";
        case VtkhdfType::PartitionedDataSetCollection:
            return "PartitionedDataSetCollection";
        case VtkhdfType::MultiBlockDataSet:
            return "MultiBlockDataSet";
    }
    return "UnstructuredGrid";
}

bool vtkhdf_is_leaf(const std::string& rType) {
    return rType == "UnstructuredGrid" || rType == "PolyData";
}

bool vtkhdf_is_composite(const std::string& rType) {
    return rType == "PartitionedDataSetCollection" || rType == "MultiBlockDataSet";
}

std::string vtkhdf_friendly_type(I64 VtkType) {
    switch (VtkType) {
        case 0:
            return "empty cell (VTK type 0)";
        case 2:
            return "poly-vertex (VTK type 2)";
        case 4:
            return "poly-line (VTK type 4)";
        case 6:
            return "triangle-strip (VTK type 6)";
        case 11:
            return "voxel (VTK type 11)";
        default:
            return "VTK type " + std::to_string(VtkType);
    }
}

NDArray vtkhdf_i64_array(const I64Vec& rV) {
    NDArray a = NDArray::Uninit(DType::Int64, {rV.size()});
    if (!rV.empty())
        std::memcpy(a.Data(), rV.data(), rV.size() * sizeof(I64));
    return a;
}

/// Rows `rIdx` of `rA` (rank >= 1), as an owning array.
NDArray vtkhdf_gather_rows(const NDArray& rA, const std::vector<std::size_t>& rIdx) {
    std::vector<std::size_t> shape = rA.Shape();
    if (shape.empty())
        shape.push_back(1);
    const std::size_t row_elems = shape[0] == 0 ? 0 : rA.Size() / shape[0];
    shape[0] = rIdx.size();
    NDArray out = NDArray::Uninit(rA.Dtype(), shape);
    const std::size_t row_bytes = row_elems * dtype_size(rA.Dtype());
    const std::byte* src = rA.Data();
    std::byte* dst = out.Data();
    for (std::size_t i = 0; i < rIdx.size(); ++i)
        if (row_bytes)
            std::memcpy(dst + i * row_bytes, src + rIdx[i] * row_bytes, row_bytes);
    return out;
}

/// Concatenate arrays along axis 0. Equal dtypes copy; mixed dtypes widen to
/// Float64 (any float involved) or Int64, as numpy's promotion would.
NDArray vtkhdf_concat_rows(const std::vector<const NDArray*>& rParts, const std::string& rWhat) {
    if (rParts.empty())
        return NDArray();
    bool same = true;
    bool any_float = false;
    for (const NDArray* p : rParts) {
        same = same && p->Dtype() == rParts[0]->Dtype();
        any_float = any_float || p->Dtype() == DType::Float32 || p->Dtype() == DType::Float64;
    }
    std::vector<std::size_t> tail = rParts[0]->Shape();
    if (tail.empty())
        tail.push_back(1);
    std::size_t rows = 0;
    for (const NDArray* p : rParts) {
        std::vector<std::size_t> s = p->Shape();
        if (s.empty())
            s.push_back(1);
        if (!std::equal(s.begin() + 1, s.end(), tail.begin() + 1, tail.end()))
            throw ReadError("meshio++: vtkhdf: cannot merge '" + rWhat +
                            "': its shape differs between pieces");
        rows += s[0];
    }
    const DType dt = same ? rParts[0]->Dtype() : (any_float ? DType::Float64 : DType::Int64);
    tail[0] = rows;
    NDArray out = NDArray::Uninit(dt, tail);
    std::size_t at = 0;
    for (const NDArray* p : rParts) {
        if (p->Size() == 0)
            continue;
        if (same) {
            std::memcpy(out.Data() + at, p->Data(), p->Nbytes());
            at += p->Nbytes();
        } else {
            for (std::size_t i = 0; i < p->Size(); ++i) {
                if (dt == DType::Float64)
                    reinterpret_cast<double*>(out.Data() + at)[i] = read_double(*p, i);
                else
                    reinterpret_cast<I64*>(out.Data() + at)[i] = read_int(*p, i);
            }
            at += p->Size() * dtype_size(dt);
        }
    }
    return out;
}

/// Per-piece `count+1` boundary runs -> global per-cell END offsets.
I64Vec vtkhdf_join_offsets(const I64Vec& rRaw, const I64Vec& rCounts, const I64Vec& rTotals) {
    if (rCounts.empty())
        return {};
    I64 expected = 0;
    for (I64 c : rCounts)
        expected += c + 1;
    if (static_cast<I64>(rRaw.size()) != expected)
        throw ReadError("meshio++: vtkhdf: an Offsets dataset holds " +
                        std::to_string(rRaw.size()) + " entries; the piece counts require " +
                        std::to_string(expected));
    I64Vec out;
    out.reserve(static_cast<std::size_t>(expected) - rCounts.size());
    std::size_t at = 0;
    I64 base = 0;
    for (std::size_t j = 0; j < rCounts.size(); ++j) {
        ++at;  // drop this piece's leading boundary
        for (I64 i = 0; i < rCounts[j]; ++i)
            out.push_back(rRaw[at++] + base);
        base += rTotals[j];
    }
    return out;
}

// ---------------------------------------------------------------------------
// Steps
// ---------------------------------------------------------------------------
class VtkhdfSteps {
public:
    explicit VtkhdfSteps(hid_t G) {
        if (!h5::exists(G, "Steps"))
            return;
        mGroup = h5::open_group(G, "Steps");
        mTransient = true;
        if (!h5::has_attr(mGroup, "NSteps"))
            throw ReadError("meshio++: vtkhdf: '" + vtkhdf_group_name(G) + "/Steps' has no NSteps");
        const I64 n = h5::read_attr_int(mGroup, "NSteps");
        if (n < 1)
            throw ReadError("meshio++: vtkhdf: NSteps is zero; the file has no data");
        mCount = static_cast<std::size_t>(n);
        if (h5::exists(mGroup, "Values")) {
            const NDArray v = h5::read_dataset(mGroup, "Values");
            if (v.Size() != mCount)
                throw ReadError("meshio++: vtkhdf: Steps/Values has " + std::to_string(v.Size()) +
                                " entries but NSteps is " + std::to_string(mCount));
            for (std::size_t i = 0; i < mCount; ++i)
                mValues.push_back(read_double(v, i));
        }
    }

    bool Transient() const { return mTransient; }
    /// Steps a caller may pick from: a static file has exactly one.
    std::size_t Count() const { return mTransient ? mCount : 1; }
    std::optional<double> Time(std::size_t K) const {
        if (K < mValues.size())
            return mValues[K];
        return std::nullopt;
    }
    const std::vector<double>& Values() const { return mValues; }

    /// `Steps/<name>[k]` (or `[k, col]`); nullopt if absent or too short.
    std::optional<I64> Scalar(const std::string& rName, std::size_t K, std::size_t Col = 0) const {
        if (!mTransient || !h5::exists(mGroup, rName))
            return std::nullopt;
        const auto shape = h5::dataset_shape(mGroup, rName);
        if (K >= shape[0])
            return std::nullopt;
        const I64Vec row = vtkhdf_int_rows(mGroup, rName, K, 1);
        if (row.empty())
            return std::nullopt;
        return row[std::min(Col, row.size() - 1)];
    }

    /// Start row of array `rName` for step `K`; nullopt marks a static array.
    std::optional<I64> DataStart(const std::string& rKind, const std::string& rName,
                                 std::size_t K) const {
        if (!mTransient)
            return std::nullopt;
        const std::string path = rKind + "DataOffsets/" + rName;
        if (!h5::exists(mGroup, rKind + "DataOffsets") || !h5::exists(mGroup, path))
            return std::nullopt;
        const auto shape = h5::dataset_shape(mGroup, path);
        if (K >= shape[0])
            return std::nullopt;
        return vtkhdf_int_rows(mGroup, path, K, 1)[0];
    }

    /// `(ncomp, ntuples)` of a field-data array at step `K`, if recorded.
    std::optional<std::pair<I64, I64>> FieldSize(const std::string& rName, std::size_t K) const {
        if (!mTransient)
            return std::nullopt;
        const std::string path = "FieldDataSizes/" + rName;
        if (!h5::exists(mGroup, "FieldDataSizes") || !h5::exists(mGroup, path))
            return std::nullopt;
        const auto shape = h5::dataset_shape(mGroup, path);
        if (K >= shape[0])
            return std::nullopt;
        const I64Vec row = vtkhdf_int_rows(mGroup, path, K, 1);
        if (row.size() < 2)
            return std::nullopt;
        return std::make_pair(row[0], row[1]);
    }

    /// `(first piece, piece count)` of step `K` among `Total` pieces.
    std::pair<std::size_t, std::size_t> Parts(std::size_t K, std::size_t Total) const {
        if (!mTransient)
            return {0, Total};
        std::optional<I64> nparts = Scalar("NumberOfParts", K);
        std::size_t np;
        if (nparts)
            np = vtkhdf_to_size(*nparts, "NumberOfParts");
        else
            np = (Total >= mCount && Total % mCount == 0) ? Total / mCount : Total;
        std::optional<I64> part0 = Scalar("PartOffsets", K);
        std::size_t p0;
        if (part0)
            p0 = vtkhdf_to_size(*part0, "PartOffsets");
        else
            p0 = (Total >= mCount * np) ? K * np : 0;
        if (p0 + np > Total)
            throw ReadError("meshio++: vtkhdf: step " + std::to_string(K) + " names pieces [" +
                            std::to_string(p0) + ", " + std::to_string(p0 + np) +
                            ") but the file holds " + std::to_string(Total));
        return {p0, np};
    }

private:
    Hid mGroup;
    bool mTransient = false;
    std::size_t mCount = 1;
    std::vector<double> mValues;
};

// ---------------------------------------------------------------------------
// One UnstructuredGrid / PolyData group, one step, a window of its pieces
// ---------------------------------------------------------------------------
struct VtkhdfLeaf {
    NDArray mPoints;
    I64Vec mConn;   ///< window-global point ids
    I64Vec mEnds;   ///< one END offset per cell into mConn
    I64Vec mTypes;  ///< VTK type id per cell
    bool mHasPoly = false;
    I64Vec mFaceConn, mFaceEnd, mToFaces, mPolyEnd;
    I64Vec mPieceOfCell;  ///< window piece index of each cell
    std::size_t mNumPieces = 0;
    std::map<std::string, NDArray> mPointData;
    std::map<std::string, NDArray> mCellRaw;
    std::map<std::string, NDArray> mFieldData;
    std::optional<double> mTime;
};

struct VtkhdfTopology {
    I64Vec mConn, mEnds, mCounts;
    I64 mC0 = 0;
};

VtkhdfTopology vtkhdf_gather_topology(hid_t Grp, const VtkhdfSteps& rSteps, std::size_t K,
                                      std::size_t Col, std::size_t Part0, std::size_t First,
                                      std::size_t Cnt, const I64Vec& rPtStarts) {
    const I64Vec ncells_all = vtkhdf_ints(Grp, "NumberOfCells");
    const I64Vec nconn_all = vtkhdf_ints(Grp, "NumberOfConnectivityIds");
    const std::size_t lo = Part0 + First, hi = Part0 + First + Cnt;
    if (hi > ncells_all.size() || hi > nconn_all.size())
        throw ReadError("meshio++: vtkhdf: '" + vtkhdf_group_name(Grp) +
                        "' has fewer pieces than the points describe");
    const I64 c_step =
        rSteps.Scalar("CellOffsets", K, Col).value_or(vtkhdf_sum(ncells_all, 0, Part0));
    const I64 n_step =
        rSteps.Scalar("ConnectivityIdOffsets", K, Col).value_or(vtkhdf_sum(nconn_all, 0, Part0));
    VtkhdfTopology t;
    t.mCounts.assign(ncells_all.begin() + static_cast<std::ptrdiff_t>(lo),
                     ncells_all.begin() + static_cast<std::ptrdiff_t>(hi));
    const I64Vec nconn(nconn_all.begin() + static_cast<std::ptrdiff_t>(lo),
                       nconn_all.begin() + static_cast<std::ptrdiff_t>(hi));
    t.mC0 = c_step + vtkhdf_sum(ncells_all, Part0, lo);
    const I64 n0 = n_step + vtkhdf_sum(nconn_all, Part0, lo);
    const I64 c_total = vtkhdf_sum(t.mCounts, 0, t.mCounts.size());
    const I64 n_total = vtkhdf_sum(nconn, 0, nconn.size());
    if (t.mC0 < 0 || n0 < 0)
        throw ReadError("meshio++: vtkhdf: negative offset in the Steps tables");
    t.mConn = vtkhdf_int_rows(Grp, "Connectivity", static_cast<std::size_t>(n0),
                              static_cast<std::size_t>(n_total));
    const I64Vec offs = vtkhdf_int_rows(Grp, "Offsets", static_cast<std::size_t>(t.mC0) + lo,
                                        static_cast<std::size_t>(c_total) + Cnt);
    // piece-local -> window-global point ids
    std::size_t at = 0;
    for (std::size_t j = 0; j < nconn.size(); ++j)
        for (I64 i = 0; i < nconn[j]; ++i)
            t.mConn[at++] += rPtStarts[j];
    t.mEnds = vtkhdf_join_offsets(offs, t.mCounts, nconn);
    return t;
}

void vtkhdf_gather_faces(hid_t G, std::size_t Part0, std::size_t First, std::size_t Cnt,
                         const I64Vec& rPtStarts, const VtkhdfTopology& rTopo, VtkhdfLeaf& rLeaf) {
    if (!h5::exists(G, "FaceConnectivity"))
        return;
    const I64Vec nfaces_all = vtkhdf_ints(G, "NumberOfFaces");
    const I64Vec nfconn_all = vtkhdf_ints(G, "NumberOfFaceConnectivityIds");
    const I64Vec npf_all = vtkhdf_ints(G, "NumberOfPolyhedronToFaceIds");
    const std::size_t lo = Part0 + First, hi = Part0 + First + Cnt;
    // Derived from the counts, not from Steps/FaceConnectivityOffsets and friends:
    // VTK's writer leaves the last piece out of those tables for partitioned steps.
    const I64 f0 = vtkhdf_sum(nfaces_all, 0, lo);
    const I64 fc0 = vtkhdf_sum(nfconn_all, 0, lo);
    const I64 pf0 = vtkhdf_sum(npf_all, 0, lo);
    const I64Vec nfaces(nfaces_all.begin() + static_cast<std::ptrdiff_t>(lo),
                        nfaces_all.begin() + static_cast<std::ptrdiff_t>(hi));
    const I64Vec nfconn(nfconn_all.begin() + static_cast<std::ptrdiff_t>(lo),
                        nfconn_all.begin() + static_cast<std::ptrdiff_t>(hi));
    const I64Vec npf(npf_all.begin() + static_cast<std::ptrdiff_t>(lo),
                     npf_all.begin() + static_cast<std::ptrdiff_t>(hi));
    rLeaf.mHasPoly = true;

    rLeaf.mFaceConn = vtkhdf_int_rows(G, "FaceConnectivity", static_cast<std::size_t>(fc0),
                                      static_cast<std::size_t>(vtkhdf_sum(nfconn, 0, Cnt)));
    std::size_t at = 0;
    for (std::size_t j = 0; j < nfconn.size(); ++j)
        for (I64 i = 0; i < nfconn[j]; ++i)
            rLeaf.mFaceConn[at++] += rPtStarts[j];

    const I64Vec foffs =
        vtkhdf_int_rows(G, "FaceOffsets", static_cast<std::size_t>(f0) + lo,
                        static_cast<std::size_t>(vtkhdf_sum(nfaces, 0, Cnt)) + Cnt);
    rLeaf.mFaceEnd = vtkhdf_join_offsets(foffs, nfaces, nfconn);

    rLeaf.mToFaces = vtkhdf_int_rows(G, "PolyhedronToFaces", static_cast<std::size_t>(pf0),
                                     static_cast<std::size_t>(vtkhdf_sum(npf, 0, Cnt)));
    const I64Vec fstarts = vtkhdf_starts(nfaces);
    at = 0;
    for (std::size_t j = 0; j < npf.size(); ++j)
        for (I64 i = 0; i < npf[j]; ++i)
            rLeaf.mToFaces[at++] += fstarts[j];

    const I64 c_total = vtkhdf_sum(rTopo.mCounts, 0, rTopo.mCounts.size());
    const I64Vec poffs =
        vtkhdf_int_rows(G, "PolyhedronOffsets", static_cast<std::size_t>(rTopo.mC0) + lo,
                        static_cast<std::size_t>(c_total) + Cnt);
    rLeaf.mPolyEnd = vtkhdf_join_offsets(poffs, rTopo.mCounts, npf);
}

NDArray vtkhdf_read_span(hid_t G, const std::string& rKind, const std::string& rName,
                         const VtkhdfSteps& rSteps, std::size_t K, I64 InStep, std::size_t Count) {
    const std::string path = rKind + "Data/" + rName;
    const I64 start = rSteps.DataStart(rKind, rName, K).value_or(0) + InStep;
    if (start < 0)
        throw ReadError("meshio++: vtkhdf: negative start row for " + path);
    const auto shape = h5::dataset_shape(G, path);
    if (static_cast<std::size_t>(start) + Count > shape[0])
        throw ReadError("meshio++: vtkhdf: " + path + " has " + std::to_string(shape[0]) +
                        " rows; step " + std::to_string(K) + " needs rows [" +
                        std::to_string(start) + ", " + std::to_string(start + Count) + ")");
    return h5::read_dataset_rows(G, path, static_cast<std::size_t>(start), Count);
}

void vtkhdf_read_field_data(hid_t G, const VtkhdfSteps& rSteps, std::size_t K,
                            const ReadOptions& rOpts, VtkhdfLeaf& rLeaf) {
    if (rOpts.mPointsOnly || !h5::exists(G, "FieldData"))
        return;
    for (const auto& name : vtkhdf_dataset_names(G, "FieldData")) {
        if (!rOpts.WantsArray(name))
            continue;
        const std::string path = "FieldData/" + name;
        const auto start = rSteps.DataStart("Field", name, K);
        const auto size = rSteps.FieldSize(name, K);
        if (start && size) {
            const I64 ncomp = size->first, ntuples = size->second;
            NDArray arr = h5::read_dataset_rows(G, path, static_cast<std::size_t>(*start),
                                                static_cast<std::size_t>(ntuples));
            if (ncomp > 1 && arr.Shape().size() == 1)
                arr.Reshape({static_cast<std::size_t>(ntuples), static_cast<std::size_t>(ncomp)});
            rLeaf.mFieldData.emplace(name, std::move(arr));
        } else {
            rLeaf.mFieldData.emplace(name, h5::read_dataset(G, path));
        }
    }
}

/// Read step `K` of group `G`; `Only` < 0 keeps every piece of the step.
VtkhdfLeaf vtkhdf_read_leaf(hid_t G, const std::string& rKind, const VtkhdfSteps& rSteps,
                            std::size_t K, std::int64_t Only, const ReadOptions& rOpts) {
    const I64Vec npts_all = vtkhdf_ints(G, "NumberOfPoints");
    const auto [part0, nparts] = rSteps.Parts(K, npts_all.size());
    std::size_t first = 0, cnt = nparts;
    if (Only >= 0) {
        if (static_cast<std::size_t>(Only) >= nparts)
            throw ReadError("meshio++: vtkhdf: piece " + std::to_string(Only) +
                            " is out of range (" + std::to_string(nparts) + ")");
        first = static_cast<std::size_t>(Only);
        cnt = 1;
    }
    const std::size_t lo = part0 + first, hi = part0 + first + cnt;
    const I64Vec npts(npts_all.begin() + static_cast<std::ptrdiff_t>(lo),
                      npts_all.begin() + static_cast<std::ptrdiff_t>(hi));
    const I64 p_step = rSteps.Scalar("PointOffsets", K).value_or(vtkhdf_sum(npts_all, 0, part0));
    const I64 in_step_pts = vtkhdf_sum(npts_all, part0, lo);
    const I64 n_pts = vtkhdf_sum(npts, 0, npts.size());
    if (p_step < 0)
        throw ReadError("meshio++: vtkhdf: negative PointOffsets");

    VtkhdfLeaf leaf;
    leaf.mNumPieces = cnt;
    leaf.mPoints =
        h5::read_dataset_rows(G, "Points", static_cast<std::size_t>(p_step + in_step_pts),
                              static_cast<std::size_t>(n_pts));
    if (leaf.mPoints.Shape().size() != 2)
        throw ReadError("meshio++: vtkhdf: Points must be a two-dimensional dataset");
    const I64Vec pt_starts = vtkhdf_starts(npts);

    I64 c_in_step = 0;
    if (rKind == "UnstructuredGrid") {
        VtkhdfTopology t = vtkhdf_gather_topology(G, rSteps, K, 0, part0, first, cnt, pt_starts);
        leaf.mConn = std::move(t.mConn);
        leaf.mEnds = std::move(t.mEnds);
        leaf.mTypes =
            vtkhdf_int_rows(G, "Types", static_cast<std::size_t>(t.mC0),
                            static_cast<std::size_t>(vtkhdf_sum(t.mCounts, 0, t.mCounts.size())));
        vtkhdf_gather_faces(G, part0, first, cnt, pt_starts, t, leaf);
        for (std::size_t j = 0; j < t.mCounts.size(); ++j)
            leaf.mPieceOfCell.insert(leaf.mPieceOfCell.end(),
                                     static_cast<std::size_t>(t.mCounts[j]), static_cast<I64>(j));
        c_in_step = vtkhdf_sum(vtkhdf_ints(G, "NumberOfCells"), part0, lo);
    } else {
        // Topology is stored per category (Vertices, Lines, Polygons, Strips), but
        // VTK's cell order -- and so the order of CellData -- is piece-major:
        // piece 0's cells in category order, then piece 1's (measured with
        // vtkHDFWriter). Gather each category, then interleave.
        struct Row {
            I64 mPiece, mCat, mStart, mSize;
        };
        std::vector<Row> rows;
        I64Vec conn_all;
        I64Vec piece_cells(npts_all.size(), 0);
        for (std::size_t col = 0; col < 4; ++col) {
            if (!h5::exists(G, kVtkhdfCategories[col]))
                continue;
            Hid sub = h5::open_group(G, kVtkhdfCategories[col]);
            const I64Vec nc_all = vtkhdf_ints(sub, "NumberOfCells");
            for (std::size_t i = 0; i < nc_all.size() && i < piece_cells.size(); ++i)
                piece_cells[i] += nc_all[i];
            VtkhdfTopology t =
                vtkhdf_gather_topology(sub, rSteps, K, col, part0, first, cnt, pt_starts);
            const I64 base = static_cast<I64>(conn_all.size());
            conn_all.insert(conn_all.end(), t.mConn.begin(), t.mConn.end());
            std::size_t cell = 0;
            for (std::size_t j = 0; j < t.mCounts.size(); ++j) {
                for (I64 i = 0; i < t.mCounts[j]; ++i, ++cell) {
                    const I64 end = t.mEnds[cell];
                    const I64 begin = cell == 0 ? 0 : t.mEnds[cell - 1];
                    rows.push_back(
                        {static_cast<I64>(j), static_cast<I64>(col), base + begin, end - begin});
                }
            }
        }
        std::stable_sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
            return a.mPiece != b.mPiece ? a.mPiece < b.mPiece : a.mCat < b.mCat;
        });
        I64 acc = 0;
        for (const Row& r : rows) {
            for (I64 i = 0; i < r.mSize; ++i)
                leaf.mConn.push_back(conn_all[static_cast<std::size_t>(r.mStart + i)]);
            acc += r.mSize;
            leaf.mEnds.push_back(acc);
            I64 t;
            if (r.mCat == 0)
                t = r.mSize == 1 ? 1 : 2;
            else if (r.mCat == 1)
                t = r.mSize == 2 ? 3 : 4;
            else if (r.mCat == 2)
                t = r.mSize == 3 ? 5 : (r.mSize == 4 ? 9 : 7);
            else
                t = 6;
            leaf.mTypes.push_back(t);
            leaf.mPieceOfCell.push_back(r.mPiece);
        }
        c_in_step = vtkhdf_sum(piece_cells, part0, lo);
    }

    const std::size_t n_cells = leaf.mTypes.size();
    if (!rOpts.mPointsOnly) {
        for (const auto& name : vtkhdf_dataset_names(G, "CellData"))
            if (rOpts.WantsArray(name))
                leaf.mCellRaw.emplace(
                    name, vtkhdf_read_span(G, "Cell", name, rSteps, K, c_in_step, n_cells));
        for (const auto& name : vtkhdf_dataset_names(G, "PointData"))
            if (rOpts.WantsArray(name))
                leaf.mPointData.emplace(name,
                                        vtkhdf_read_span(G, "Point", name, rSteps, K, in_step_pts,
                                                         static_cast<std::size_t>(n_pts)));
    }
    vtkhdf_read_field_data(G, rSteps, K, rOpts, leaf);
    if (rSteps.Transient())
        leaf.mTime = rSteps.Time(K);
    return leaf;
}

// ---------------------------------------------------------------------------
// cells
// ---------------------------------------------------------------------------
void vtkhdf_append_cell_data(Mesh& rMesh, const VtkhdfLeaf& rLeaf, std::size_t R0, std::size_t R1) {
    for (const auto& kv : rLeaf.mCellRaw)
        rMesh.AppendCellData(kv.first, detail::slice_rows(kv.second, R0, R1));
}

void vtkhdf_append_cell_data_rows(Mesh& rMesh, const VtkhdfLeaf& rLeaf,
                                  const std::vector<std::size_t>& rRows) {
    for (const auto& kv : rLeaf.mCellRaw)
        rMesh.AppendCellData(kv.first, vtkhdf_gather_rows(kv.second, rRows));
}

/// Decode polyhedra `[Lo, Hi)` into `polyhedron<N>` blocks (N = unique node count),
/// bucketed in first-seen order -- the convention the VTU, OpenFOAM, MED and CGNS
/// readers all use -- so the run may be reordered; `rPerm` records where each lands.
std::size_t vtkhdf_polyhedron_run(Mesh& rMesh, const VtkhdfLeaf& rLeaf, std::size_t Lo,
                                  std::size_t Hi, I64Vec& rPerm, std::size_t At) {
    if (!rLeaf.mHasPoly)
        throw ReadError(
            "meshio++: vtkhdf: a cell has VTK type 42 (polyhedron) but the file carries no "
            "FaceConnectivity/FaceOffsets/PolyhedronToFaces/PolyhedronOffsets datasets");
    std::vector<std::vector<I64Vec>> cells;
    std::vector<std::size_t> counts;
    for (std::size_t c = Lo; c < Hi; ++c) {
        if (c >= rLeaf.mPolyEnd.size())
            throw ReadError("meshio++: vtkhdf: PolyhedronOffsets is shorter than the cell count");
        const I64 ps = c > 0 ? rLeaf.mPolyEnd[c - 1] : 0;
        const I64 pe = rLeaf.mPolyEnd[c];
        if (ps < 0 || pe < ps || static_cast<std::size_t>(pe) > rLeaf.mToFaces.size())
            throw ReadError("meshio++: vtkhdf: PolyhedronOffsets entry out of range");
        std::vector<I64Vec> faces;
        I64Vec uniq;
        for (I64 i = ps; i < pe; ++i) {
            const I64 f = rLeaf.mToFaces[static_cast<std::size_t>(i)];
            if (f < 0 || static_cast<std::size_t>(f) >= rLeaf.mFaceEnd.size())
                throw ReadError("meshio++: vtkhdf: PolyhedronToFaces names a missing face");
            const I64 fs = f > 0 ? rLeaf.mFaceEnd[static_cast<std::size_t>(f) - 1] : 0;
            const I64 fe = rLeaf.mFaceEnd[static_cast<std::size_t>(f)];
            if (fs < 0 || fe < fs || static_cast<std::size_t>(fe) > rLeaf.mFaceConn.size())
                throw ReadError("meshio++: vtkhdf: FaceOffsets entry out of range");
            faces.emplace_back(rLeaf.mFaceConn.begin() + fs, rLeaf.mFaceConn.begin() + fe);
            uniq.insert(uniq.end(), faces.back().begin(), faces.back().end());
        }
        std::sort(uniq.begin(), uniq.end());
        uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
        counts.push_back(uniq.size());
        cells.push_back(std::move(faces));
    }
    std::vector<std::size_t> order;
    std::map<std::size_t, std::vector<std::size_t>> groups;
    for (std::size_t i = 0; i < cells.size(); ++i) {
        if (groups.find(counts[i]) == groups.end())
            order.push_back(counts[i]);
        groups[counts[i]].push_back(i);
    }
    for (std::size_t n : order) {
        const auto& idx = groups[n];
        std::vector<std::vector<I64Vec>> group;
        std::vector<std::size_t> rows;
        group.reserve(idx.size());
        for (std::size_t i : idx) {
            group.push_back(std::move(cells[i]));
            rows.push_back(Lo + i);
            rPerm[Lo + i] = static_cast<I64>(At++);
        }
        rMesh.AddPolyhedronBlock("polyhedron" + std::to_string(n), std::move(group));
        vtkhdf_append_cell_data_rows(rMesh, rLeaf, rows);
    }
    return At;
}

/// Cells from the flat VTK arrays; returns `perm` (file cell -> block-major
/// index, or -1 when `Lenient` dropped it). Polyhedra bucket per piece.
I64Vec vtkhdf_build_cells(Mesh& rMesh, const VtkhdfLeaf& rLeaf, bool Lenient) {
    const auto& vmap = vtk_to_meshio_type();
    const auto& nmap = num_nodes_per_cell();
    const std::size_t n = rLeaf.mTypes.size();
    I64Vec perm(n, 0);
    std::size_t at = 0;
    std::size_t a = 0;
    while (a < n) {
        std::size_t b = a + 1;
        while (b < n && rLeaf.mTypes[b] == rLeaf.mTypes[a])
            ++b;
        const I64 vtk_type = rLeaf.mTypes[a];
        if (vtk_type == kVtkhdfPolyhedron) {
            std::size_t lo = a;
            for (std::size_t i = a + 1; i <= b; ++i) {
                if (i == b || rLeaf.mPieceOfCell[i] != rLeaf.mPieceOfCell[i - 1]) {
                    at = vtkhdf_polyhedron_run(rMesh, rLeaf, lo, i, perm, at);
                    lo = i;
                }
            }
            a = b;
            continue;
        }
        auto it = vmap.find(static_cast<int>(vtk_type));
        if (it == vmap.end() || vtk_type == 0) {
            if (!Lenient)
                throw ReadError("meshio++: vtkhdf: " + vtkhdf_friendly_type(vtk_type) +
                                " has no meshio++ cell type; set lenient to skip such cells");
            log::warn("meshio++: vtkhdf: skipping {} cell(s) of {}.", b - a,
                      vtkhdf_friendly_type(vtk_type));
            for (std::size_t i = a; i < b; ++i)
                perm[i] = -1;
            a = b;
            continue;
        }
        const std::string& meshio_type = it->second;
        const I64 first = a > 0 ? rLeaf.mEnds[a - 1] : 0;
        if (is_special_cell(meshio_type)) {
            // Split further wherever the per-cell node count changes.
            std::size_t i = a;
            while (i < b) {
                const I64 prev = i == 0 ? 0 : rLeaf.mEnds[i - 1];
                const I64 sz = rLeaf.mEnds[i] - prev;
                std::size_t j = i + 1;
                while (j < b && rLeaf.mEnds[j] - rLeaf.mEnds[j - 1] == sz)
                    ++j;
                const std::size_t m = j - i;
                NDArray data = NDArray::Uninit(DType::Int64, {m, static_cast<std::size_t>(sz)});
                I64* out = data.As<I64>();
                for (std::size_t r = 0; r < m; ++r)
                    for (I64 c = 0; c < sz; ++c)
                        out[r * static_cast<std::size_t>(sz) + static_cast<std::size_t>(c)] =
                            rLeaf.mConn[static_cast<std::size_t>(rLeaf.mEnds[i + r] - sz + c)];
                rMesh.AddCellBlock(meshio_type, std::move(data));
                vtkhdf_append_cell_data(rMesh, rLeaf, i, j);
                for (std::size_t r = i; r < j; ++r)
                    perm[r] = static_cast<I64>(at++);
                i = j;
            }
        } else {
            auto nit = nmap.find(meshio_type);
            if (nit == nmap.end())
                throw ReadError("meshio++: vtkhdf: '" + meshio_type +
                                "' has no fixed node count and is not a polygon or Lagrange cell");
            const std::size_t nn = static_cast<std::size_t>(nit->second);
            const std::size_t m = b - a;
            for (std::size_t r = 0; r < m; ++r) {
                const I64 prev = r == 0 ? first : rLeaf.mEnds[a + r - 1];
                if (rLeaf.mEnds[a + r] - prev != static_cast<I64>(nn))
                    throw ReadError("meshio++: vtkhdf: '" + meshio_type + "' cells must have " +
                                    std::to_string(nn) + " nodes, found " +
                                    std::to_string(rLeaf.mEnds[a + r] - prev));
            }
            const std::vector<int> order = vtk_to_meshio_order(static_cast<int>(vtk_type));
            NDArray data = NDArray::Uninit(DType::Int64, {m, nn});
            I64* out = data.As<I64>();
            for (std::size_t r = 0; r < m; ++r) {
                const std::size_t base = static_cast<std::size_t>(rLeaf.mEnds[a + r]) - nn;
                for (std::size_t j = 0; j < nn; ++j)
                    out[r * nn + j] =
                        rLeaf
                            .mConn[base + (order.empty() ? j : static_cast<std::size_t>(order[j]))];
            }
            rMesh.AddCellBlock(meshio_type, std::move(data));
            vtkhdf_append_cell_data(rMesh, rLeaf, a, b);
            for (std::size_t r = a; r < b; ++r)
                perm[r] = static_cast<I64>(at++);
        }
        a = b;
    }
    return perm;
}

/// Build the mesh from a leaf. `rPieceNames[p]` names window piece `p` (empty
/// vector = attach no regions).
Mesh vtkhdf_assemble(VtkhdfLeaf& rLeaf, const ReadOptions& rOpts,
                     const std::vector<std::string>& rPieceNames) {
    Mesh mesh;
    mesh.AssignPoints(std::move(rLeaf.mPoints));
    const I64Vec perm = vtkhdf_build_cells(mesh, rLeaf, rOpts.mLenient);
    for (auto& kv : rLeaf.mPointData)
        mesh.AddPointData(kv.first, std::move(kv.second));
    for (auto& kv : rLeaf.mFieldData)
        mesh.AddFieldData(kv.first, std::move(kv.second));
    if (rLeaf.mTime) {
        NDArray t = NDArray::Uninit(DType::Float64, {1});
        t.As<double>()[0] = *rLeaf.mTime;
        mesh.AddFieldData(kSequenceTimeKey, std::move(t));
    }
    if (!rPieceNames.empty()) {
        std::vector<I64Vec> entries(rPieceNames.size());
        for (std::size_t i = 0; i < perm.size(); ++i)
            if (perm[i] >= 0)
                entries[static_cast<std::size_t>(rLeaf.mPieceOfCell[i])].push_back(perm[i]);
        for (std::size_t p = 0; p < rPieceNames.size(); ++p)
            // tag = the piece's position in the file, which a writer sorts by to restore
            // the order (the mesh stores regions sorted by name, so it cannot carry it).
            mesh.AddRegion(Region(rPieceNames[p], RegionKind::Cell, -1,
                                  static_cast<std::int64_t>(p), vtkhdf_i64_array(entries[p])));
    }
    return mesh;
}

// ---------------------------------------------------------------------------
// composites
// ---------------------------------------------------------------------------
struct VtkhdfBlock {
    std::string mName;  ///< the link name (the block's name)
    std::string mPath;  ///< the assembly path, used when names collide
    Hid mGroup;
};

void vtkhdf_walk_assembly(hid_t Node, const std::string& rPrefix, std::vector<VtkhdfBlock>& rOut) {
    for (const auto& name : h5::group_links_crt(Node)) {
        if (h5::is_soft_link(Node, name)) {
            if (H5Oexists_by_name(Node, name.c_str(), H5P_DEFAULT) <= 0)
                throw ReadError("meshio++: vtkhdf: Assembly link '" + rPrefix + name +
                                "' points at '" + h5::soft_link_target(Node, name) +
                                "', which does not exist");
            VtkhdfBlock b;
            b.mName = name;
            b.mPath = rPrefix + name;
            b.mGroup = Hid(H5Gopen2(Node, name.c_str(), H5P_DEFAULT), H5Gclose);
            rOut.push_back(std::move(b));
        } else {
            Hid child(H5Gopen2(Node, name.c_str(), H5P_DEFAULT), H5Gclose);
            if (child.Valid())
                vtkhdf_walk_assembly(child, rPrefix + name + "/", rOut);
        }
    }
}

std::vector<VtkhdfBlock> vtkhdf_composite_blocks(hid_t G, const std::string& rKind) {
    std::vector<VtkhdfBlock> blocks;
    if (h5::exists(G, "Assembly")) {
        Hid asm_ = h5::open_group(G, "Assembly");
        vtkhdf_walk_assembly(asm_, "", blocks);
    } else {
        for (const auto& name : h5::group_links_crt(G)) {
            if (name == "Assembly" || name == "Steps" || name == "FieldData")
                continue;
            Hid child(H5Gopen2(G, name.c_str(), H5P_DEFAULT), H5Gclose);
            if (child.Valid() && h5::has_attr(child, "Type")) {
                VtkhdfBlock b;
                b.mName = b.mPath = name;
                b.mGroup = std::move(child);
                blocks.push_back(std::move(b));
            }
        }
    }
    std::set<std::string> seen;
    bool collide = false;
    for (const auto& b : blocks)
        collide = !seen.insert(b.mName).second || collide;
    if (collide)
        for (auto& b : blocks)
            b.mName = b.mPath;
    if (rKind == "PartitionedDataSetCollection" && !blocks.empty()) {
        bool all_indexed = true;
        for (const auto& b : blocks)
            all_indexed = all_indexed && h5::has_attr(b.mGroup, "Index");
        if (all_indexed) {
            std::vector<std::pair<I64, std::size_t>> keyed;
            for (std::size_t i = 0; i < blocks.size(); ++i)
                keyed.emplace_back(h5::read_attr_int(blocks[i].mGroup, "Index"), i);
            std::stable_sort(keyed.begin(), keyed.end(),
                             [](const auto& a, const auto& b) { return a.first < b.first; });
            std::vector<VtkhdfBlock> sorted;
            for (const auto& kv : keyed)
                sorted.push_back(std::move(blocks[kv.second]));
            blocks = std::move(sorted);
        }
    }
    return blocks;
}

VtkhdfSteps vtkhdf_check_no_time(hid_t G, const std::string& rWhere) {
    VtkhdfSteps steps(G);
    if (steps.Count() > 1)
        throw ReadError(
            "meshio++: vtkhdf: " + rWhere + " carries " + std::to_string(steps.Count()) +
            " time steps; transient composite (PartitionedDataSetCollection / "
            "MultiBlockDataSet) datasets are not supported -- only UnstructuredGrid and "
            "PolyData files are transient");
    return steps;
}

/// Concatenate leaves without reordering or welding; arrays not present in every
/// leaf are dropped with a warning. The first leaf's field data is kept.
VtkhdfLeaf vtkhdf_concat_leaves(std::vector<VtkhdfLeaf>& rLeaves) {
    VtkhdfLeaf out;
    if (rLeaves.empty())
        return out;
    std::vector<const NDArray*> pts;
    I64 pt_off = 0, conn_off = 0, piece_off = 0;
    for (auto& l : rLeaves) {
        pts.push_back(&l.mPoints);
        out.mHasPoly = out.mHasPoly || l.mHasPoly;
        for (I64 c : l.mConn)
            out.mConn.push_back(c + pt_off);
        for (I64 e : l.mEnds)
            out.mEnds.push_back(e + conn_off);
        out.mTypes.insert(out.mTypes.end(), l.mTypes.begin(), l.mTypes.end());
        for (I64 p : l.mPieceOfCell)
            out.mPieceOfCell.push_back(p + piece_off);
        pt_off += static_cast<I64>(l.mPoints.Shape().empty() ? 0 : l.mPoints.Shape()[0]);
        conn_off += static_cast<I64>(l.mConn.size());
        piece_off += static_cast<I64>(l.mNumPieces);
        out.mNumPieces += l.mNumPieces;
    }
    if (out.mHasPoly) {
        // Every leaf's polyhedron arrays are window-global, so each shifts by the
        // totals of everything before it; a leaf without polyhedra still owes its
        // cells an empty span in mPolyEnd.
        I64 pto = 0, fco = 0, fo = 0, tfo = 0;
        for (auto& l : rLeaves) {
            for (I64 c : l.mFaceConn)
                out.mFaceConn.push_back(c + pto);
            for (I64 e : l.mFaceEnd)
                out.mFaceEnd.push_back(e + fco);
            for (I64 f : l.mToFaces)
                out.mToFaces.push_back(f + fo);
            if (l.mHasPoly) {
                for (I64 e : l.mPolyEnd)
                    out.mPolyEnd.push_back(e + tfo);
            } else {
                out.mPolyEnd.insert(out.mPolyEnd.end(), l.mTypes.size(), tfo);
            }
            pto += static_cast<I64>(l.mPoints.Shape().empty() ? 0 : l.mPoints.Shape()[0]);
            fco += static_cast<I64>(l.mFaceConn.size());
            fo += static_cast<I64>(l.mFaceEnd.size());
            tfo += static_cast<I64>(l.mToFaces.size());
        }
    }
    out.mPoints = vtkhdf_concat_rows(pts, "Points");

    auto common = [&](auto member, const char* pKind) {
        std::map<std::string, std::vector<const NDArray*>> merged;
        std::set<std::string> everything;
        for (auto& l : rLeaves)
            for (auto& kv : l.*member)
                everything.insert(kv.first);
        std::map<std::string, NDArray> result;
        std::vector<std::string> dropped;
        for (const auto& name : everything) {
            std::vector<const NDArray*> parts;
            for (auto& l : rLeaves) {
                auto it = (l.*member).find(name);
                if (it == (l.*member).end())
                    break;
                parts.push_back(&it->second);
            }
            if (parts.size() == rLeaves.size())
                result.emplace(name, vtkhdf_concat_rows(parts, name));
            else
                dropped.push_back(name);
        }
        if (!dropped.empty()) {
            std::string list;
            for (const auto& d : dropped)
                list += (list.empty() ? "" : ", ") + d;
            log::warn("meshio++: vtkhdf: dropping {} array(s) not present in every piece: {}",
                      pKind, list);
        }
        return result;
    };
    out.mPointData = common(&VtkhdfLeaf::mPointData, "point");
    out.mCellRaw = common(&VtkhdfLeaf::mCellRaw, "cell");
    out.mFieldData = std::move(rLeaves.front().mFieldData);
    return out;
}

Mesh vtkhdf_read_composite(hid_t G, const std::string& rKind, const ReadOptions& rOpts) {
    std::vector<VtkhdfBlock> blocks = vtkhdf_composite_blocks(G, rKind);
    const VtkhdfSteps root_steps = vtkhdf_check_no_time(G, kVtkhdfRoot);
    (void)rOpts.ResolveTimeStep(root_steps.Count());  // 0 / -1 are the only steps of a static file

    struct Planned {
        std::string mType;
        std::size_t mNumParts;
    };
    std::vector<Planned> plan;
    for (const auto& b : blocks) {
        const std::string btype =
            h5::has_attr(b.mGroup, "Type") ? h5::read_attr_string(b.mGroup, "Type") : "";
        if (!vtkhdf_is_leaf(btype))
            throw ReadError("meshio++: vtkhdf: composite block '" + b.mName + "' has Type '" +
                            btype + "'; only UnstructuredGrid and PolyData blocks are supported");
        (void)vtkhdf_check_no_time(b.mGroup, "block '" + b.mName + "'");
        plan.push_back({btype, vtkhdf_ints(b.mGroup, "NumberOfPoints").size()});
    }

    VtkhdfLeaf root_fd;
    vtkhdf_read_field_data(G, root_steps, 0, rOpts, root_fd);

    if (rOpts.mPieceSet) {
        std::size_t total = 0;
        for (const auto& p : plan)
            total += p.mNumParts;
        std::size_t flat = rOpts.ResolvePiece(total);
        for (std::size_t i = 0; i < plan.size(); ++i) {
            if (flat < plan[i].mNumParts) {
                const VtkhdfSteps steps(blocks[i].mGroup);
                VtkhdfLeaf leaf = vtkhdf_read_leaf(blocks[i].mGroup, plan[i].mType, steps, 0,
                                                   static_cast<std::int64_t>(flat), rOpts);
                for (auto& kv : root_fd.mFieldData)
                    leaf.mFieldData.emplace(kv.first, std::move(kv.second));
                return vtkhdf_assemble(leaf, rOpts, {});
            }
            flat -= plan[i].mNumParts;
        }
    }

    std::vector<VtkhdfLeaf> leaves;
    std::vector<std::string> piece_names;
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        const VtkhdfSteps steps(blocks[i].mGroup);
        leaves.push_back(vtkhdf_read_leaf(blocks[i].mGroup, plan[i].mType, steps, 0, -1, rOpts));
        for (std::size_t j = 0; j < plan[i].mNumParts; ++j)
            piece_names.push_back(plan[i].mNumParts == 1
                                      ? blocks[i].mName
                                      : blocks[i].mName + "/piece_" + std::to_string(j));
    }
    VtkhdfLeaf merged = vtkhdf_concat_leaves(leaves);
    for (auto& kv : root_fd.mFieldData)
        merged.mFieldData.emplace(kv.first, std::move(kv.second));
    return vtkhdf_assemble(merged, rOpts, piece_names);
}

// ---------------------------------------------------------------------------
// entry: open + validate
// ---------------------------------------------------------------------------
void vtkhdf_check_version(hid_t Root) {
    if (!h5::has_attr(Root, "Version"))
        throw ReadError("meshio++: vtkhdf: '" + vtkhdf_group_name(Root) +
                        "' has no Version attribute");
    const I64Vec v = h5::read_attr_int_array(Root, "Version");
    if (v.size() != 2)
        throw ReadError("meshio++: vtkhdf: Version must be two integers");
    if (v[0] != 1 && v[0] != 2)
        throw ReadError("meshio++: vtkhdf: unsupported VTKHDF version " + std::to_string(v[0]) +
                        "." + std::to_string(v[1]) + "; this build reads 1.x and 2.x");
    if (v[0] == 2 && v[1] > kVtkhdfKnownMinor)
        log::warn(
            "meshio++: vtkhdf: file declares VTKHDF {}.{}, newer than the 2.{} this build "
            "knows; reading what it recognizes.",
            v[0], v[1], kVtkhdfKnownMinor);
}

Hid vtkhdf_open(const std::string& rPath, Hid& rRoot, std::string& rKind) {
    Hid f(H5Fopen(rPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    if (!f.Valid())
        throw ReadError("meshio++: vtkhdf: cannot open '" + rPath + "' as an HDF5 file");
    if (!h5::exists(f, kVtkhdfRoot))
        throw ReadError("meshio++: vtkhdf: '" + rPath + "' is an HDF5 file with no /" +
                        kVtkhdfRoot +
                        " group; if it holds another HDF5 format, pass that format explicitly");
    rRoot = h5::open_group(f, kVtkhdfRoot);
    rKind = h5::has_attr(rRoot, "Type") ? h5::read_attr_string(rRoot, "Type") : "";
    vtkhdf_check_version(rRoot);
    return f;
}

// ---------------------------------------------------------------------------
// writing
// ---------------------------------------------------------------------------
void vtkhdf_put(hid_t Loc, const std::string& rName, const NDArray& rArr, int Gzip) {
    const int level = (Gzip >= 0 && rArr.Nbytes() >= kVtkhdfGzipMinBytes) ? Gzip : -1;
    h5::write_dataset(Loc, rName, rArr, level);
}

void vtkhdf_put_i64(hid_t Loc, const std::string& rName, const I64Vec& rV, int Gzip) {
    vtkhdf_put(Loc, rName, vtkhdf_i64_array(rV), Gzip);
}

void vtkhdf_check_name(const std::string& rName, const char* pWhat) {
    if (rName.empty() || rName.find('/') != std::string::npos)
        throw WriteError(std::string("meshio++: vtkhdf: ") + pWhat + " name '" + rName +
                         "' cannot be an HDF5 dataset name (empty or containing '/')");
}

/// A rank <= 2 view of `rArr` (trailing dimensions folded into columns). No copy.
NDArray vtkhdf_flat_view(const NDArray& rArr) {
    if (rArr.Shape().size() <= 2)
        return NDArray::MakeView(rArr.Dtype(),
                                 rArr.Shape().empty() ? std::vector<std::size_t>{1} : rArr.Shape(),
                                 const_cast<std::byte*>(rArr.Data()));
    std::size_t cols_total = 1;
    for (std::size_t i = 1; i < rArr.Shape().size(); ++i)
        cols_total *= rArr.Shape()[i];
    return NDArray::MakeView(rArr.Dtype(), {rArr.Shape()[0], cols_total},
                             const_cast<std::byte*>(rArr.Data()));
}

/// The points as `(n, 3)` in their own dtype (floats stay floats), rows `pUsed` only.
NDArray vtkhdf_points3(const NDArray& rPts, std::size_t NumPoints, std::size_t Dim,
                       const std::vector<std::size_t>* pUsed) {
    if (Dim > 3)
        throw WriteError("meshio++: vtkhdf: points must be (n, 1..3), got dimension " +
                         std::to_string(Dim));
    const DType dt = (rPts.Dtype() == DType::Float32 || rPts.Dtype() == DType::Float64)
                         ? rPts.Dtype()
                         : DType::Float64;
    const std::size_t n = pUsed ? pUsed->size() : NumPoints;
    NDArray out(dt, {n, 3});  // zero-filled: the padded columns stay 0
    for (std::size_t r = 0; r < n; ++r) {
        const std::size_t src = pUsed ? (*pUsed)[r] : r;
        for (std::size_t c = 0; c < Dim; ++c) {
            if (dt == DType::Float32)
                out.As<float>()[r * 3 + c] = static_cast<float>(read_double(rPts, src * Dim + c));
            else
                out.As<double>()[r * 3 + c] = read_double(rPts, src * Dim + c);
        }
    }
    return out;
}

struct VtkhdfCells {
    I64Vec mConn, mSizes, mTypes;
    I64Vec mFaceConn, mFaceSizes, mFacesPerCell;
    bool mHasPoly = false;
    std::vector<std::size_t> mRows;  ///< global (block-major) index of each collected cell
};

/// Flat VTK arrays for the mesh's cells, blocks visited in `rBlockOrder`; only
/// cells with `pMask[g]` set when a mask is given. Node ids are the mesh's own.
VtkhdfCells vtkhdf_collect(const Mesh& rMesh, const std::vector<std::size_t>& rBlockOrder,
                           const std::vector<std::size_t>& rBases, const std::vector<char>* pMask) {
    const auto& tmap = meshio_to_vtk_type();
    VtkhdfCells out;
    for (std::size_t bi : rBlockOrder) {
        const auto cb = rMesh.Cells(bi);
        const std::size_t nc = cb.NumCells();
        const std::string type = cb.Type();
        auto wanted = [&](std::size_t r) { return !pMask || (*pMask)[rBases[bi] + r] != 0; };
        if (cb.IsPolyhedron()) {
            out.mHasPoly = true;
            for (std::size_t r = 0; r < nc; ++r) {
                if (!wanted(r))
                    continue;
                I64Vec uniq;
                std::size_t nfaces = cb.NumFaces(r);
                for (std::size_t f = 0; f < nfaces; ++f) {
                    const auto face = cb.Face(r, f);
                    out.mFaceSizes.push_back(static_cast<I64>(face.second));
                    for (std::size_t j = 0; j < face.second; ++j) {
                        out.mFaceConn.push_back(face.first[j]);
                        uniq.push_back(face.first[j]);
                    }
                }
                std::sort(uniq.begin(), uniq.end());
                uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
                out.mConn.insert(out.mConn.end(), uniq.begin(), uniq.end());
                out.mSizes.push_back(static_cast<I64>(uniq.size()));
                out.mTypes.push_back(kVtkhdfPolyhedron);
                out.mFacesPerCell.push_back(static_cast<I64>(nfaces));
                out.mRows.push_back(rBases[bi] + r);
            }
            continue;
        }
        const bool polygon = type.rfind("polygon", 0) == 0;
        if (cb.IsRagged()) {
            for (std::size_t r = 0; r < nc; ++r) {
                if (!wanted(r))
                    continue;
                for (std::size_t j = 0; j < cb.RowSize(r); ++j)
                    out.mConn.push_back(cb.Row(r)[j]);
                out.mSizes.push_back(static_cast<I64>(cb.RowSize(r)));
                out.mTypes.push_back(tmap.at("polygon"));
                out.mFacesPerCell.push_back(0);
                out.mRows.push_back(rBases[bi] + r);
            }
            continue;
        }
        auto it = polygon ? tmap.find("polygon") : tmap.find(type);
        if (it == tmap.end())
            throw WriteError("meshio++: vtkhdf: cell type '" + type + "' has no VTK cell type id");
        const NDArray& conn = cb.Conn();
        const std::size_t k = cols(conn);
        const std::vector<int> order = polygon ? std::vector<int>{} : meshio_to_vtk_order(type);
        for (std::size_t r = 0; r < nc; ++r) {
            if (!wanted(r))
                continue;
            for (std::size_t j = 0; j < k; ++j) {
                const std::size_t col = order.empty() ? j : static_cast<std::size_t>(order[j]);
                out.mConn.push_back(read_int(conn, r * k + col));
            }
            out.mSizes.push_back(static_cast<I64>(k));
            out.mTypes.push_back(it->second);
            out.mFacesPerCell.push_back(0);
            out.mRows.push_back(rBases[bi] + r);
        }
    }
    return out;
}

/// Cell data over the collected rows, in collected order (`rRows` are global indices).
NDArray vtkhdf_cell_data_rows(const Mesh& rMesh, const std::string& rName,
                              const std::vector<std::size_t>& rRows,
                              const std::vector<std::size_t>& rBases) {
    const std::size_t nblocks = rMesh.NumCellBlocks();
    if (rMesh.CellDataNumBlocks(rName) != nblocks)
        throw WriteError("meshio++: vtkhdf: cell_data '" + rName + "' has data for " +
                         std::to_string(rMesh.CellDataNumBlocks(rName)) + " of " +
                         std::to_string(nblocks) + " cell blocks");
    const NDArray& first = rMesh.CellData(rName, 0);
    std::vector<std::size_t> shape = first.Shape();
    if (shape.empty())
        shape.push_back(1);
    std::size_t row_elems = 1;
    for (std::size_t i = 1; i < shape.size(); ++i)
        row_elems *= shape[i];
    shape[0] = rRows.size();
    NDArray out = NDArray::Uninit(first.Dtype(), shape);
    const std::size_t row_bytes = row_elems * dtype_size(first.Dtype());
    for (std::size_t i = 0; i < rRows.size(); ++i) {
        const std::size_t g = rRows[i];
        const std::size_t bi =
            static_cast<std::size_t>(
                std::upper_bound(rBases.begin(),
                                 rBases.begin() + static_cast<std::ptrdiff_t>(nblocks), g) -
                rBases.begin()) -
            1;
        const NDArray& blk = rMesh.CellData(rName, bi);
        if (blk.Dtype() != first.Dtype())
            throw WriteError("meshio++: vtkhdf: cell_data '" + rName +
                             "' changes dtype between cell blocks");
        if (row_bytes)
            std::memcpy(out.Data() + i * row_bytes, blk.Data() + (g - rBases[bi]) * row_bytes,
                        row_bytes);
    }
    return out;
}

void vtkhdf_write_data_groups(hid_t Grp, const Mesh& rMesh, const VtkhdfCells& rCells,
                              const std::vector<std::size_t>& rBases,
                              const std::vector<std::size_t>* pUsed, int Gzip) {
    Hid pd = h5::create_group(Grp, "PointData");
    for (const auto& name : rMesh.PointDataNames()) {
        vtkhdf_check_name(name, "point_data");
        const NDArray& arr = rMesh.PointData(name);
        if (pUsed) {
            const NDArray g = vtkhdf_gather_rows(arr, *pUsed);
            vtkhdf_put(pd, name, vtkhdf_flat_view(g), Gzip);
        } else {
            vtkhdf_put(pd, name, vtkhdf_flat_view(arr), Gzip);
        }
    }
    Hid cd = h5::create_group(Grp, "CellData");
    for (const auto& name : rMesh.CellDataNames()) {
        vtkhdf_check_name(name, "cell_data");
        if (rMesh.CellDataNumBlocks(name) == 0)
            continue;
        const NDArray rows = vtkhdf_cell_data_rows(rMesh, name, rCells.mRows, rBases);
        vtkhdf_put(cd, name, vtkhdf_flat_view(rows), Gzip);
    }
    Hid fd = h5::create_group(Grp, "FieldData");
    for (const auto& name : rMesh.FieldDataNames()) {
        vtkhdf_check_name(name, "field_data");
        vtkhdf_put(fd, name, vtkhdf_flat_view(rMesh.FieldData(name)), Gzip);
    }
}

I64Vec vtkhdf_boundaries(const I64Vec& rSizes) {
    I64Vec out(rSizes.size() + 1, 0);
    for (std::size_t i = 0; i < rSizes.size(); ++i)
        out[i + 1] = out[i] + rSizes[i];
    return out;
}

void vtkhdf_write_type_and_version(hid_t Grp, const char* pType, std::pair<int, int> Version) {
    h5::write_attr_int_array(Grp, "Version", {Version.first, Version.second});
    h5::write_attr_string_fixed(Grp, "Type", pType);
}

std::vector<std::size_t> vtkhdf_block_bases(const Mesh& rMesh) {
    std::vector<std::size_t> bases(rMesh.NumCellBlocks() + 1, 0);
    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi)
        bases[bi + 1] = bases[bi] + rMesh.Cells(bi).NumCells();
    return bases;
}

/// The geometry datasets of an UnstructuredGrid group: `NumberOf*`, Points, Connectivity,
/// Offsets, Types and (with polyhedra) the four face arrays. `pUsed` prunes the points.
void vtkhdf_write_ug_geometry(hid_t Grp, const Mesh& rMesh, const VtkhdfCells& rCells,
                              const std::vector<std::size_t>* pUsed, int Gzip) {
    const VtkhdfCells& cells = rCells;
    const std::size_t n_points = pUsed ? pUsed->size() : rMesh.NumPoints();
    vtkhdf_put_i64(Grp, "NumberOfPoints", {static_cast<I64>(n_points)}, -1);
    vtkhdf_put_i64(Grp, "NumberOfCells", {static_cast<I64>(cells.mTypes.size())}, -1);
    vtkhdf_put_i64(Grp, "NumberOfConnectivityIds", {static_cast<I64>(cells.mConn.size())}, -1);
    vtkhdf_put(Grp, "Points",
               vtkhdf_points3(rMesh.Points(), rMesh.NumPoints(), rMesh.PointDim(), pUsed), Gzip);
    vtkhdf_put_i64(Grp, "Connectivity", cells.mConn, Gzip);
    vtkhdf_put_i64(Grp, "Offsets", vtkhdf_boundaries(cells.mSizes), Gzip);
    NDArray types = NDArray::Uninit(DType::UInt8, {cells.mTypes.size()});
    for (std::size_t i = 0; i < cells.mTypes.size(); ++i)
        types.As<std::uint8_t>()[i] = static_cast<std::uint8_t>(cells.mTypes[i]);
    vtkhdf_put(Grp, "Types", types, Gzip);
    if (cells.mHasPoly) {
        I64Vec to_faces(cells.mFaceSizes.size());
        for (std::size_t i = 0; i < to_faces.size(); ++i)
            to_faces[i] = static_cast<I64>(i);
        vtkhdf_put_i64(Grp, "NumberOfFaces", {static_cast<I64>(cells.mFaceSizes.size())}, -1);
        vtkhdf_put_i64(Grp, "NumberOfFaceConnectivityIds",
                       {static_cast<I64>(cells.mFaceConn.size())}, -1);
        vtkhdf_put_i64(Grp, "NumberOfPolyhedronToFaceIds", {static_cast<I64>(to_faces.size())}, -1);
        vtkhdf_put_i64(Grp, "FaceConnectivity", cells.mFaceConn, Gzip);
        vtkhdf_put_i64(Grp, "FaceOffsets", vtkhdf_boundaries(cells.mFaceSizes), Gzip);
        vtkhdf_put_i64(Grp, "PolyhedronToFaces", to_faces, Gzip);
        vtkhdf_put_i64(Grp, "PolyhedronOffsets", vtkhdf_boundaries(cells.mFacesPerCell), Gzip);
    }
}

/// An UnstructuredGrid group for `pCells` (sorted global cell indices) or, when
/// null, for the whole mesh. A subset is pruned to the points it uses.
void vtkhdf_write_ug_group(hid_t Grp, const Mesh& rMesh, const std::vector<std::size_t>* pCells,
                           int Gzip, std::pair<int, int> Version) {
    vtkhdf_write_type_and_version(Grp, "UnstructuredGrid", Version);
    const std::vector<std::size_t> bases = vtkhdf_block_bases(rMesh);
    std::vector<std::size_t> order(rMesh.NumCellBlocks());
    for (std::size_t i = 0; i < order.size(); ++i)
        order[i] = i;
    std::vector<char> mask;
    if (pCells) {
        mask.assign(bases.back(), 0);
        for (std::size_t g : *pCells)
            mask[g] = 1;
    }
    VtkhdfCells cells = vtkhdf_collect(rMesh, order, bases, pCells ? &mask : nullptr);

    std::vector<std::size_t> used;
    if (pCells) {
        // Prune to the referenced points and renumber.
        std::vector<char> seen(rMesh.NumPoints(), 0);
        for (I64 c : cells.mConn)
            seen[static_cast<std::size_t>(c)] = 1;
        for (I64 c : cells.mFaceConn)
            seen[static_cast<std::size_t>(c)] = 1;
        std::vector<I64> remap(rMesh.NumPoints(), -1);
        for (std::size_t i = 0; i < seen.size(); ++i)
            if (seen[i]) {
                remap[i] = static_cast<I64>(used.size());
                used.push_back(i);
            }
        for (I64& c : cells.mConn)
            c = remap[static_cast<std::size_t>(c)];
        for (I64& c : cells.mFaceConn)
            c = remap[static_cast<std::size_t>(c)];
    }
    vtkhdf_write_ug_geometry(Grp, rMesh, cells, pCells ? &used : nullptr, Gzip);
    vtkhdf_write_data_groups(Grp, rMesh, cells, bases, pCells ? &used : nullptr, Gzip);
}

int vtkhdf_pd_category(const std::string& rType) {
    if (rType == "vertex")
        return 0;
    if (rType == "line")
        return 1;
    if (rType == "triangle" || rType == "quad" || rType.rfind("polygon", 0) == 0)
        return 2;
    throw WriteError("meshio++: vtkhdf: PolyData cannot hold '" + rType +
                     "' cells (only vertex, line, triangle, quad and polygon); write "
                     "UnstructuredGrid to keep them");
}

void vtkhdf_write_polydata_group(hid_t Grp, const Mesh& rMesh, int Gzip,
                                 std::pair<int, int> Version) {
    std::vector<int> kinds;
    for (std::size_t bi = 0; bi < rMesh.NumCellBlocks(); ++bi)
        kinds.push_back(vtkhdf_pd_category(rMesh.Cells(bi).Type()));
    // VTK's canonical PolyData cell order is Vertices, Lines, Polygons, Strips: regroup
    // the blocks (stably); cell_data follows the same order.
    std::vector<std::size_t> order;
    for (int want = 0; want < 3; ++want)
        for (std::size_t bi = 0; bi < kinds.size(); ++bi)
            if (kinds[bi] == want)
                order.push_back(bi);
    const std::vector<std::size_t> bases = vtkhdf_block_bases(rMesh);
    vtkhdf_write_type_and_version(Grp, "PolyData", Version);
    vtkhdf_put_i64(Grp, "NumberOfPoints", {static_cast<I64>(rMesh.NumPoints())}, -1);
    vtkhdf_put(Grp, "Points",
               vtkhdf_points3(rMesh.Points(), rMesh.NumPoints(), rMesh.PointDim(), nullptr), Gzip);

    I64Vec conn[3], sizes[3];
    for (std::size_t bi : order) {
        const auto cb = rMesh.Cells(bi);
        const int cat = kinds[bi];
        for (std::size_t r = 0; r < cb.NumCells(); ++r) {
            std::size_t k;
            if (cb.IsRagged()) {
                k = cb.RowSize(r);
                for (std::size_t j = 0; j < k; ++j)
                    conn[cat].push_back(cb.Row(r)[j]);
            } else {
                const NDArray& c = cb.Conn();
                k = cols(c);
                for (std::size_t j = 0; j < k; ++j)
                    conn[cat].push_back(read_int(c, r * k + j));
            }
            if (cat == 0 && k != 1)
                throw WriteError("meshio++: vtkhdf: vertex cells must have one node");
            if (cat == 1 && k != 2)
                throw WriteError("meshio++: vtkhdf: line cells must have two nodes");
            sizes[cat].push_back(static_cast<I64>(k));
        }
    }
    for (std::size_t col = 0; col < 4; ++col) {
        Hid sub = h5::create_group(Grp, kVtkhdfCategories[col]);
        const I64Vec empty;
        const I64Vec& c = col < 3 ? conn[col] : empty;
        const I64Vec& s = col < 3 ? sizes[col] : empty;
        vtkhdf_put_i64(sub, "NumberOfCells", {static_cast<I64>(s.size())}, -1);
        vtkhdf_put_i64(sub, "NumberOfConnectivityIds", {static_cast<I64>(c.size())}, -1);
        vtkhdf_put_i64(sub, "Connectivity", c, Gzip);
        vtkhdf_put_i64(sub, "Offsets", vtkhdf_boundaries(s), Gzip);
    }
    VtkhdfCells cells;
    for (std::size_t bi : order)
        for (std::size_t r = 0; r < rMesh.Cells(bi).NumCells(); ++r)
            cells.mRows.push_back(bases[bi] + r);
    vtkhdf_write_data_groups(Grp, rMesh, cells, bases, nullptr, Gzip);
}

// VTKHDF's composite carving is now the shared detail::carve_by_region (roadmap
// §1.1, adopted so EnSight Gold's multi-part writer does not reimplement the same
// region-vs-block fallback). `detail::MeshPart` is VTKHDF's own former `VtkhdfPiece`
// under a shared name; `StrictNames=true, RejectSlash=true` reproduce this format's
// original behaviour exactly -- a name collision or a `/` in a region name is a
// `WriteError`, never a silent rename of an Assembly link.
using VtkhdfPiece = detail::MeshPart;

std::vector<VtkhdfPiece> vtkhdf_carve(const Mesh& rMesh) {
    return detail::carve_by_region(rMesh, "vtkhdf", /*StrictNames=*/true, /*RejectSlash=*/true);
}

std::pair<int, int> vtkhdf_resolve_version(VtkhdfVersion Requested, VtkhdfType Type, bool HasPoly) {
    const auto needed = vtkhdf_needed_version(Type, HasPoly);
    if (Requested.mMajor == 0 && Requested.mMinor == 0)
        return std::max(needed, std::make_pair(2, 0));
    const std::pair<int, int> v{Requested.mMajor, Requested.mMinor};
    if ((v.first != 1 && v.first != 2) || v.second < 0)
        throw WriteError(
            "meshio++: vtkhdf: version must be (major, minor) with major 1 or 2, got (" +
            std::to_string(v.first) + ", " + std::to_string(v.second) + ")");
    if (v.first == 2 && v.second > kVtkhdfKnownMinor)
        throw WriteError("meshio++: vtkhdf: version " + std::to_string(v.first) + "." +
                         std::to_string(v.second) + " is newer than the 2." +
                         std::to_string(kVtkhdfKnownMinor) + " this build knows");
    if (v < needed) {
        const std::string feature = HasPoly ? std::string("polyhedral cells")
                                            : std::string("Type '") + vtkhdf_type_name(Type) + "'";
        throw WriteError("meshio++: vtkhdf: " + feature + " needs VTKHDF " +
                         std::to_string(needed.first) + "." + std::to_string(needed.second) +
                         " or newer; version " + std::to_string(v.first) + "." +
                         std::to_string(v.second) + " cannot express it");
    }
    return v;
}

}  // namespace

// ---------------------------------------------------------------------------
// public API
// ---------------------------------------------------------------------------
void write_vtkhdf(const std::string& rPath, const Mesh& rMesh, int GzipLevel, VtkhdfType Type,
                  VtkhdfVersion Version) {
    if (GzipLevel > 9)
        throw WriteError("meshio++: vtkhdf: gzip level must be 0-9");
    bool has_poly = false;
    for (const auto cb : rMesh.CellRange())
        has_poly = has_poly || cb.IsPolyhedron();
    const auto ver = vtkhdf_resolve_version(Version, Type, has_poly);

    h5::SilenceErrors silence;
    Hid f = h5::create_file(rPath);
    Hid root = h5::create_group_crt(f, kVtkhdfRoot);
    if (Type == VtkhdfType::UnstructuredGrid) {
        vtkhdf_write_ug_group(root, rMesh, nullptr, GzipLevel, ver);
    } else if (Type == VtkhdfType::PolyData) {
        vtkhdf_write_polydata_group(root, rMesh, GzipLevel, ver);
    } else {
        const std::vector<VtkhdfPiece> pieces = vtkhdf_carve(rMesh);
        vtkhdf_write_type_and_version(root, vtkhdf_type_name(Type), ver);
        const bool pdc = Type == VtkhdfType::PartitionedDataSetCollection;
        for (std::size_t i = 0; i < pieces.size(); ++i) {
            // A composite's root may hold only blocks and the Assembly (vtkHDFReader reads
            // every other root group as a block and then fails), so mesh-wide field_data
            // rides on each block's own FieldData group, which the block writer does.
            Hid bg = h5::create_group_crt(root, pieces[i].mName);
            if (pdc)
                h5::write_attr_int(bg, "Index", static_cast<I64>(i));
            vtkhdf_write_ug_group(bg, rMesh, &pieces[i].mCells, GzipLevel, ver);
        }
        Hid asm_ = h5::create_group_crt(root, "Assembly");
        for (const auto& p : pieces) {
            const std::string target = std::string("/") + kVtkhdfRoot + "/" + p.mName;
            if (pdc) {
                Hid node = h5::create_group_crt(asm_, p.mName);
                h5::create_soft_link(node, p.mName, target);
            } else {
                h5::create_soft_link(asm_, p.mName, target);
            }
        }
    }
    std::string provenance;
    for (const auto& line : detail::provenance_lines(detail::SlotTier::Block))
        provenance += (provenance.empty() ? "" : "\n") + line;
    if (!provenance.empty())
        h5::write_attr_string(root, kVtkhdfProvenanceAttr, provenance);
}

Mesh read_vtkhdf(const std::string& rPath) {
    return read_vtkhdf(rPath, ReadOptions{});
}

Mesh read_vtkhdf(const std::string& rPath, const ReadOptions& rOpts) {
    h5::SilenceErrors silence;
    Hid root;
    std::string kind;
    Hid f = vtkhdf_open(rPath, root, kind);
    if (vtkhdf_is_composite(kind))
        return vtkhdf_read_composite(root, kind, rOpts);
    if (!vtkhdf_is_leaf(kind))
        throw ReadError("meshio++: vtkhdf: Type '" + kind +
                        "' is not supported; this build reads UnstructuredGrid, PolyData, "
                        "PartitionedDataSetCollection, MultiBlockDataSet");
    const VtkhdfSteps steps(root);
    const std::size_t k = rOpts.ResolveTimeStep(steps.Count());
    const std::size_t nparts = steps.Parts(k, vtkhdf_ints(root, "NumberOfPoints").size()).second;
    std::int64_t only = -1;
    if (rOpts.mPieceSet) {
        only = static_cast<std::int64_t>(rOpts.ResolvePiece(nparts));
    }
    VtkhdfLeaf leaf = vtkhdf_read_leaf(root, kind, steps, k, only, rOpts);
    std::vector<std::string> names;
    if (only < 0 && nparts > 1)
        for (std::size_t j = 0; j < nparts; ++j)
            names.push_back("piece_" + std::to_string(j));
    return vtkhdf_assemble(leaf, rOpts, names);
}

MeshMetadata read_vtkhdf_metadata(const std::string& rPath, const ReadOptions& rOpts) {
    h5::SilenceErrors silence;
    Hid root;
    std::string kind;
    Hid f = vtkhdf_open(rPath, root, kind);
    std::vector<double> times;
    if (vtkhdf_is_leaf(kind))
        times = VtkhdfSteps(root).Values();

    MeshMetadata meta;
    bool native = false;
    if (kind == "UnstructuredGrid") {
        const VtkhdfSteps steps(root);
        const std::size_t k = rOpts.ResolveTimeStep(steps.Count());
        const I64Vec npts_all = vtkhdf_ints(root, "NumberOfPoints");
        const auto [part0, nparts] = steps.Parts(k, npts_all.size());
        const I64Vec ncells_all = vtkhdf_ints(root, "NumberOfCells");
        const I64 c0 = steps.Scalar("CellOffsets", k).value_or(vtkhdf_sum(ncells_all, 0, part0));
        const I64 n_cells = vtkhdf_sum(ncells_all, part0, part0 + nparts);
        const I64Vec types = vtkhdf_int_rows(root, "Types", static_cast<std::size_t>(c0),
                                             static_cast<std::size_t>(n_cells));
        // Polyhedra bucket by node count (needs the faces) and variable-size types need
        // the offsets: both are answered by a full read instead.
        const bool has_poly =
            std::find(types.begin(), types.end(), kVtkhdfPolyhedron) != types.end();
        if (!has_poly && !detail::cells_need_offsets(types)) {
            native = true;
            meta.mNumPoints = static_cast<std::size_t>(vtkhdf_sum(npts_all, part0, part0 + nparts));
            const auto shape = h5::dataset_shape(root, "Points");
            meta.mPointDim = shape.size() > 1 ? shape[1] : 3;
            meta.mCellBlocks = detail::summarize_cells({}, types);
            meta.mPointDataNames = vtkhdf_dataset_names(root, "PointData");
            meta.mCellDataNames = vtkhdf_dataset_names(root, "CellData");
            meta.mFieldDataNames = vtkhdf_dataset_names(root, "FieldData");
            meta.mHasBBox = false;
            meta.mFellBackToFullRead = false;
        }
    }
    if (!native) {
        meta = metadata_from_mesh(read_vtkhdf(rPath, rOpts));
        meta.mFellBackToFullRead = true;
    }
    meta.mTimeValues = times;
    meta.mFormat = "vtkhdf";
    return meta;
}

// ---------------------------------------------------------------------------
// transient writer
// ---------------------------------------------------------------------------
namespace {

void vtkhdf_set_int_attr(hid_t Loc, const char* pName, I64 Value) {
    if (h5::has_attr(Loc, pName))
        H5Adelete(Loc, pName);
    h5::write_attr_int(Loc, pName, Value);
}

using VtkhdfSeriesInfo = std::pair<DType, std::vector<std::size_t>>;  // dtype, trailing dims
using VtkhdfSeriesArray = std::pair<std::string, NDArray>;

/// Fold trailing dimensions into columns so the array is rank <= 2 (owning arrays only).
void vtkhdf_fold_to_rank2(NDArray& rArr) {
    if (rArr.Shape().size() <= 2)
        return;
    std::size_t rest = 1;
    for (std::size_t i = 1; i < rArr.Shape().size(); ++i)
        rest *= rArr.Shape()[i];
    rArr.Reshape({rArr.Shape()[0], rest});
}

VtkhdfSeriesInfo vtkhdf_series_info(const NDArray& rArr) {
    return {rArr.Dtype(), rArr.Shape().size() == 2 ? std::vector<std::size_t>{rArr.Shape()[1]}
                                                   : std::vector<std::size_t>{}};
}

/// One step's rows per chunk: the natural unit of "random access per step", capped so a
/// huge mesh does not get a multi-hundred-megabyte chunk.
std::size_t vtkhdf_series_chunk_rows(std::size_t Rows, const VtkhdfSeriesInfo& rInfo) {
    std::size_t row_elems = 1;
    for (std::size_t e : rInfo.second)
        row_elems *= e;
    const std::size_t row_bytes = std::max<std::size_t>(1, row_elems * dtype_size(rInfo.first));
    const std::size_t cap = std::max<std::size_t>(1, (std::size_t{4} << 20) / row_bytes);
    return std::max<std::size_t>(1, std::min(Rows, cap));
}

NDArray vtkhdf_series_scalar(I64 V) {
    NDArray a = NDArray::Uninit(DType::Int64, {1});
    a.As<I64>()[0] = V;
    return a;
}

NDArray vtkhdf_series_named(const VtkhdfTimeSeriesWriter::NamedArray& rArr, std::size_t Rows,
                            const char* pWhat) {
    if (rArr.mNumComponents == 0 || rArr.mValues.size() != Rows * rArr.mNumComponents)
        throw WriteError(std::string("meshio++: vtkhdf: ") + pWhat + " '" + rArr.mName +
                         "' holds " + std::to_string(rArr.mValues.size()) + " values; " +
                         std::to_string(Rows) + " x " + std::to_string(rArr.mNumComponents) +
                         " are required");
    NDArray a =
        NDArray::Uninit(DType::Float64, rArr.mNumComponents == 1
                                            ? std::vector<std::size_t>{Rows}
                                            : std::vector<std::size_t>{Rows, rArr.mNumComponents});
    if (!rArr.mValues.empty())
        std::memcpy(a.Data(), rArr.mValues.data(), rArr.mValues.size() * sizeof(double));
    return a;
}

}  // namespace

struct VtkhdfTimeSeriesWriter::Impl {
    Hid mFile;
    std::string mPath;
    int mGzip = -1;
    bool mFinalized = false;
    bool mAutoFlush = true;
    bool mGeometry = false;
    bool mAdopted = false;  // the grid came from an existing file (Append)
    bool mNamesFixed = false;
    bool mHasPoly = false;
    std::size_t mNumPoints = 0;
    std::size_t mNumCells = 0;
    std::size_t mNumSteps = 0;
    std::map<std::string, VtkhdfSeriesInfo> mPointInfo, mCellInfo, mFieldInfo;
    std::map<std::string, I64> mFieldTuples;  // tuples written so far, per field array

    Hid Root() const { return h5::open_group(mFile, kVtkhdfRoot); }

    void Flush() {
        if (mFile.Valid() && H5Fflush(mFile, H5F_SCOPE_GLOBAL) < 0)
            throw WriteError("meshio++: vtkhdf: could not flush '" + mPath + "'");
    }

    void OpenExisting() {
        mFile = Hid(H5Fopen(mPath.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose);
        if (!mFile.Valid())
            throw WriteError("meshio++: vtkhdf: cannot open '" + mPath + "' to append to it");
        const std::string why = "meshio++: vtkhdf: '" + mPath +
                                "' is not a transient UnstructuredGrid this writer can continue";
        if (!h5::exists(mFile, kVtkhdfRoot))
            throw WriteError(why);
        Hid root = Root();
        if (!h5::has_attr(root, "Type") ||
            h5::read_attr_string(root, "Type") != "UnstructuredGrid" || !h5::exists(root, "Steps"))
            throw WriteError(why);
        try {
            const I64Vec npts = vtkhdf_ints(root, "NumberOfPoints");
            const I64Vec ncells = vtkhdf_ints(root, "NumberOfCells");
            if (npts.size() != 1 || ncells.size() != 1)
                throw WriteError(why + " (it is partitioned)");
            mNumPoints = static_cast<std::size_t>(npts[0]);
            mNumCells = static_cast<std::size_t>(ncells[0]);
            mHasPoly = h5::exists(root, "FaceConnectivity");
            Hid steps = h5::open_group(root, "Steps");
            mNumSteps = static_cast<std::size_t>(h5::read_attr_int(steps, "NSteps"));
            const auto adopt = [&](const char* pGroup,
                                   std::map<std::string, VtkhdfSeriesInfo>& rInfo) {
                for (const auto& name : vtkhdf_dataset_names(root, pGroup)) {
                    Hid d(H5Dopen2(root, (std::string(pGroup) + "/" + name).c_str(), H5P_DEFAULT),
                          H5Dclose);
                    Hid dt(H5Dget_type(d), H5Tclose);
                    auto shape = h5::dataset_shape(root, std::string(pGroup) + "/" + name);
                    shape.erase(shape.begin());
                    rInfo[name] = {h5::dtype_from_h5(dt), shape};
                }
            };
            adopt("PointData", mPointInfo);
            adopt("CellData", mCellInfo);
            adopt("FieldData", mFieldInfo);
            for (const auto& kv : mFieldInfo)
                mFieldTuples[kv.first] =
                    static_cast<I64>(h5::dataset_num_rows(root, "FieldData/" + kv.first));
        } catch (const ReadError& e) {
            throw WriteError(why + " (" + e.what() + ")");
        }
        mNamesFixed = mNumSteps > 0;
        mGeometry = true;
        mAdopted = true;
    }

    void WriteGeometry(const Mesh& rMesh) {
        if (mFinalized)
            throw WriteError("meshio++: vtkhdf: the series is already finalized");
        if (mAdopted) {
            mAdopted = false;  // the one allowed re-statement of the grid
            if (rMesh.NumPoints() != mNumPoints)
                throw WriteError("meshio++: vtkhdf: the appended series has " +
                                 std::to_string(mNumPoints) + " points but the mesh has " +
                                 std::to_string(rMesh.NumPoints()));
            std::size_t cells = 0;
            for (const auto cb : rMesh.CellRange())
                cells += cb.NumCells();
            if (cells != mNumCells)
                throw WriteError("meshio++: vtkhdf: the appended series has " +
                                 std::to_string(mNumCells) + " cells but the mesh has " +
                                 std::to_string(cells));
            return;
        }
        if (mGeometry)
            throw WriteError("meshio++: vtkhdf: WritePointsCells was already called");
        const std::vector<std::size_t> bases = vtkhdf_block_bases(rMesh);
        std::vector<std::size_t> order(rMesh.NumCellBlocks());
        for (std::size_t i = 0; i < order.size(); ++i)
            order[i] = i;
        const VtkhdfCells cells = vtkhdf_collect(rMesh, order, bases, nullptr);
        mHasPoly = cells.mHasPoly;
        Hid root = h5::create_group_crt(mFile, kVtkhdfRoot);
        vtkhdf_write_type_and_version(
            root, "UnstructuredGrid",
            vtkhdf_resolve_version({}, VtkhdfType::UnstructuredGrid, mHasPoly));
        vtkhdf_write_ug_geometry(root, rMesh, cells, nullptr, mGzip);
        h5::create_group(root, "PointData");
        h5::create_group(root, "CellData");
        h5::create_group(root, "FieldData");
        std::string provenance;
        for (const auto& line : detail::provenance_lines(detail::SlotTier::Block))
            provenance += (provenance.empty() ? "" : "\n") + line;
        if (!provenance.empty())
            h5::write_attr_string(root, kVtkhdfProvenanceAttr, provenance);
        mNumPoints = rMesh.NumPoints();
        mNumCells = cells.mTypes.size();
        mGeometry = true;
        if (mAutoFlush)
            Flush();
    }

    static void CheckSet(const char* pWhat,
                         const std::map<std::string, VtkhdfSeriesInfo>& rExpected,
                         const std::vector<VtkhdfSeriesArray>& rGot, std::size_t Step) {
        for (const auto& [name, arr] : rGot) {
            auto it = rExpected.find(name);
            if (it == rExpected.end())
                throw WriteError(
                    "meshio++: vtkhdf: step " + std::to_string(Step) + " introduces " + pWhat +
                    " '" + name +
                    "', which the first step did not have; the set of array names is fixed "
                    "at the first WriteData");
            const VtkhdfSeriesInfo info = vtkhdf_series_info(arr);
            if (info != it->second)
                throw WriteError("meshio++: vtkhdf: step " + std::to_string(Step) +
                                 " changes the dtype or "
                                 "component count of " +
                                 pWhat + " '" + name + "'");
        }
        for (const auto& kv : rExpected) {
            const bool present =
                std::any_of(rGot.begin(), rGot.end(),
                            [&](const VtkhdfSeriesArray& a) { return a.first == kv.first; });
            if (!present)
                throw WriteError("meshio++: vtkhdf: step " + std::to_string(Step) + " is missing " +
                                 pWhat + " '" + kv.first + "', which the first step had");
        }
    }

    void CreateStepTables(hid_t Root, const std::vector<VtkhdfSeriesArray>& rPoint,
                          const std::vector<VtkhdfSeriesArray>& rCell,
                          const std::vector<VtkhdfSeriesArray>& rField) {
        Hid steps = h5::create_group(Root, "Steps");
        vtkhdf_set_int_attr(steps, "NSteps", 0);
        h5::create_appendable_dataset(steps, "Values", DType::Float64, {}, 256, -1);
        std::vector<const char*> zero_tables = {"PartOffsets", "PointOffsets", "CellOffsets",
                                                "ConnectivityIdOffsets"};
        if (mHasPoly)
            for (const char* t :
                 {"FaceConnectivityOffsets", "FaceOffsetsOffsets", "PolyhedronToFaceIdOffsets"})
                zero_tables.push_back(t);
        zero_tables.push_back("NumberOfParts");
        for (const char* t : zero_tables)
            h5::create_appendable_dataset(steps, t, DType::Int64, {}, 256, -1);
        h5::create_group(steps, "PointDataOffsets");
        h5::create_group(steps, "CellDataOffsets");
        h5::create_group(steps, "FieldDataOffsets");
        h5::create_group(steps, "FieldDataSizes");
        Hid pd = h5::open_group(Root, "PointData");
        for (const auto& [name, arr] : rPoint) {
            const VtkhdfSeriesInfo info = vtkhdf_series_info(arr);
            h5::create_appendable_dataset(pd, name, info.first, info.second,
                                          vtkhdf_series_chunk_rows(mNumPoints, info), mGzip);
            h5::create_appendable_dataset(steps, "PointDataOffsets/" + name, DType::Int64, {}, 256,
                                          -1);
            mPointInfo[name] = info;
        }
        Hid cd = h5::open_group(Root, "CellData");
        for (const auto& [name, arr] : rCell) {
            const VtkhdfSeriesInfo info = vtkhdf_series_info(arr);
            h5::create_appendable_dataset(cd, name, info.first, info.second,
                                          vtkhdf_series_chunk_rows(mNumCells, info), mGzip);
            h5::create_appendable_dataset(steps, "CellDataOffsets/" + name, DType::Int64, {}, 256,
                                          -1);
            mCellInfo[name] = info;
        }
        Hid fd = h5::open_group(Root, "FieldData");
        for (const auto& [name, arr] : rField) {
            const VtkhdfSeriesInfo info = vtkhdf_series_info(arr);
            h5::create_appendable_dataset(
                fd, name, info.first, info.second,
                vtkhdf_series_chunk_rows(arr.Shape().empty() ? 1 : arr.Shape()[0], info), mGzip);
            h5::create_appendable_dataset(steps, "FieldDataOffsets/" + name, DType::Int64, {}, 256,
                                          -1);
            h5::create_appendable_dataset(steps, "FieldDataSizes/" + name, DType::Int64, {2}, 256,
                                          -1);
            mFieldInfo[name] = info;
            mFieldTuples[name] = 0;
        }
        mNamesFixed = true;
    }

    void AppendStep(double Time, const std::vector<VtkhdfSeriesArray>& rPoint,
                    const std::vector<VtkhdfSeriesArray>& rCell,
                    const std::vector<VtkhdfSeriesArray>& rField) {
        if (mFinalized)
            throw WriteError("meshio++: vtkhdf: the series is already finalized");
        if (!mGeometry)
            throw WriteError("meshio++: vtkhdf: WritePointsCells must be called before WriteData");
        for (const auto& [name, arr] : rPoint)
            if (arr.Shape().empty() || arr.Shape()[0] != mNumPoints)
                throw WriteError("meshio++: vtkhdf: point_data '" + name + "' has " +
                                 std::to_string(arr.Shape().empty() ? 0 : arr.Shape()[0]) +
                                 " rows; the series has " + std::to_string(mNumPoints) + " points");
        for (const auto& [name, arr] : rCell)
            if (arr.Shape().empty() || arr.Shape()[0] != mNumCells)
                throw WriteError("meshio++: vtkhdf: cell_data '" + name + "' has " +
                                 std::to_string(arr.Shape().empty() ? 0 : arr.Shape()[0]) +
                                 " rows; the series has " + std::to_string(mNumCells) + " cells");
        Hid root = Root();
        if (mNamesFixed) {
            CheckSet("point_data", mPointInfo, rPoint, mNumSteps);
            CheckSet("cell_data", mCellInfo, rCell, mNumSteps);
            CheckSet("field_data", mFieldInfo, rField, mNumSteps);
        } else {
            CreateStepTables(root, rPoint, rCell, rField);
        }
        const I64 step = static_cast<I64>(mNumSteps);
        NDArray time = NDArray::Uninit(DType::Float64, {1});
        time.As<double>()[0] = Time;
        h5::append_rows(root, "Steps/Values", time);
        for (const char* t :
             {"PartOffsets", "PointOffsets", "CellOffsets", "ConnectivityIdOffsets"})
            h5::append_rows(root, std::string("Steps/") + t, vtkhdf_series_scalar(0));
        h5::append_rows(root, "Steps/NumberOfParts", vtkhdf_series_scalar(1));
        if (mHasPoly)
            for (const char* t :
                 {"FaceConnectivityOffsets", "FaceOffsetsOffsets", "PolyhedronToFaceIdOffsets"})
                h5::append_rows(root, std::string("Steps/") + t, vtkhdf_series_scalar(0));
        for (const auto& [name, arr] : rPoint) {
            h5::append_rows(root, "PointData/" + name, arr);
            h5::append_rows(root, "Steps/PointDataOffsets/" + name,
                            vtkhdf_series_scalar(step * static_cast<I64>(mNumPoints)));
        }
        for (const auto& [name, arr] : rCell) {
            h5::append_rows(root, "CellData/" + name, arr);
            h5::append_rows(root, "Steps/CellDataOffsets/" + name,
                            vtkhdf_series_scalar(step * static_cast<I64>(mNumCells)));
        }
        for (const auto& [name, arr] : rField) {
            const I64 tuples = static_cast<I64>(arr.Shape().empty() ? 1 : arr.Shape()[0]);
            const I64 ncomp = arr.Shape().size() == 2 ? static_cast<I64>(arr.Shape()[1]) : 1;
            h5::append_rows(root, "FieldData/" + name, arr);
            h5::append_rows(root, "Steps/FieldDataOffsets/" + name,
                            vtkhdf_series_scalar(mFieldTuples[name]));
            NDArray sizes = NDArray::Uninit(DType::Int64, {1, 2});
            sizes.As<I64>()[0] = ncomp;
            sizes.As<I64>()[1] = tuples;
            h5::append_rows(root, "Steps/FieldDataSizes/" + name, sizes);
            mFieldTuples[name] += tuples;
        }
        ++mNumSteps;
        Hid steps = h5::open_group(root, "Steps");
        vtkhdf_set_int_attr(steps, "NSteps", static_cast<I64>(mNumSteps));
        if (mAutoFlush)
            Flush();
    }

    void Finalize() {
        if (mFinalized)
            return;
        mFinalized = true;
        if (mFile.Valid()) {
            const int rc = H5Fflush(mFile, H5F_SCOPE_GLOBAL);
            mFile.Reset();
            if (rc < 0)
                throw WriteError("meshio++: vtkhdf: could not flush '" + mPath + "'");
        }
    }
};

VtkhdfTimeSeriesWriter::VtkhdfTimeSeriesWriter(const std::string& rPath, int GzipLevel,
                                               VtkhdfSeriesMode Mode)
    : mImpl(std::make_unique<Impl>()) {
    if (GzipLevel > 9)
        throw WriteError("meshio++: vtkhdf: gzip level must be 0-9");
    h5::SilenceErrors silence;
    mImpl->mPath = rPath;
    mImpl->mGzip = GzipLevel;
    if (Mode == VtkhdfSeriesMode::Append && std::filesystem::exists(rPath))
        mImpl->OpenExisting();
    else
        mImpl->mFile = h5::create_file(rPath);
}

VtkhdfTimeSeriesWriter::~VtkhdfTimeSeriesWriter() {
    if (!mImpl)
        return;
    try {
        h5::SilenceErrors silence;
        mImpl->Finalize();
    } catch (...) {
        // an exception must not leave a destructor; Finalize() is how a caller sees it
    }
}

VtkhdfTimeSeriesWriter::VtkhdfTimeSeriesWriter(VtkhdfTimeSeriesWriter&&) noexcept = default;
VtkhdfTimeSeriesWriter& VtkhdfTimeSeriesWriter::operator=(
    VtkhdfTimeSeriesWriter&& rOther) noexcept {
    if (this != &rOther) {
        try {
            if (mImpl)
                mImpl->Finalize();
        } catch (...) {
        }
        mImpl = std::move(rOther.mImpl);
    }
    return *this;
}

void VtkhdfTimeSeriesWriter::WritePointsCells(const Mesh& rMesh) {
    if (!mImpl)
        throw WriteError("meshio++: vtkhdf: this writer was moved from");
    h5::SilenceErrors silence;
    mImpl->WriteGeometry(rMesh);
}

void VtkhdfTimeSeriesWriter::WriteData(double Time, const Mesh& rMesh) {
    if (!mImpl)
        throw WriteError("meshio++: vtkhdf: this writer was moved from");
    h5::SilenceErrors silence;
    std::vector<VtkhdfSeriesArray> point, cell, field;
    for (const auto& name : rMesh.PointDataNames()) {
        vtkhdf_check_name(name, "point_data");
        point.emplace_back(name, vtkhdf_flat_view(rMesh.PointData(name)));
    }
    const std::vector<std::size_t> bases = vtkhdf_block_bases(rMesh);
    std::vector<std::size_t> rows(bases.back());
    for (std::size_t i = 0; i < rows.size(); ++i)
        rows[i] = i;
    for (const auto& name : rMesh.CellDataNames()) {
        vtkhdf_check_name(name, "cell_data");
        if (rMesh.CellDataNumBlocks(name) == 0)
            continue;
        NDArray a = vtkhdf_cell_data_rows(rMesh, name, rows, bases);
        vtkhdf_fold_to_rank2(a);
        cell.emplace_back(name, std::move(a));
    }
    for (const auto& name : rMesh.FieldDataNames()) {
        if (name == kSequenceTimeKey)
            continue;  // the step's time is `Time`, recorded in Steps/Values
        vtkhdf_check_name(name, "field_data");
        field.emplace_back(name, vtkhdf_flat_view(rMesh.FieldData(name)));
    }
    mImpl->AppendStep(Time, point, cell, field);
}

void VtkhdfTimeSeriesWriter::WriteData(double Time, const std::vector<NamedArray>& rPointData,
                                       const std::vector<NamedArray>& rCellData) {
    if (!mImpl)
        throw WriteError("meshio++: vtkhdf: this writer was moved from");
    h5::SilenceErrors silence;
    std::vector<VtkhdfSeriesArray> point, cell, field;
    for (const auto& a : rPointData) {
        vtkhdf_check_name(a.mName, "point_data");
        point.emplace_back(a.mName, vtkhdf_series_named(a, mImpl->mNumPoints, "point_data"));
    }
    for (const auto& a : rCellData) {
        vtkhdf_check_name(a.mName, "cell_data");
        cell.emplace_back(a.mName, vtkhdf_series_named(a, mImpl->mNumCells, "cell_data"));
    }
    mImpl->AppendStep(Time, point, cell, field);
}

void VtkhdfTimeSeriesWriter::Flush() {
    if (mImpl && !mImpl->mFinalized) {
        h5::SilenceErrors silence;
        mImpl->Flush();
    }
}

void VtkhdfTimeSeriesWriter::SetAutoFlush(bool Enable) {
    if (mImpl)
        mImpl->mAutoFlush = Enable;
}

bool VtkhdfTimeSeriesWriter::AutoFlush() const {
    return mImpl && mImpl->mAutoFlush;
}

void VtkhdfTimeSeriesWriter::Finalize() {
    if (mImpl) {
        h5::SilenceErrors silence;
        mImpl->Finalize();
    }
}

std::size_t VtkhdfTimeSeriesWriter::NumSteps() const {
    return mImpl ? mImpl->mNumSteps : 0;
}

bool VtkhdfTimeSeriesWriter::Finalized() const {
    return !mImpl || mImpl->mFinalized;
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
