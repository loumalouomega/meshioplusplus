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
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "../detail/library_lock.hpp"
#include "../detail/typed_view.hpp"
#include "meshioplusplus/formats/nastran_h5.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/nastran_model.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/ndarray.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "meshioplusplus/read_options.hpp"
#include "meshioplusplus/region.hpp"

namespace meshioplusplus {

namespace {

constexpr const char* kNh5Grid = "/NASTRAN/INPUT/NODE/GRID";
constexpr const char* kNh5Elements = "/NASTRAN/INPUT/ELEMENT";
constexpr const char* kNh5Properties = "/NASTRAN/INPUT/PROPERTY";
constexpr const char* kNh5Domains = "/NASTRAN/RESULT/DOMAINS";
constexpr const char* kNh5Nodal = "/NASTRAN/RESULT/NODAL";
constexpr const char* kNh5Elemental = "/NASTRAN/RESULT/ELEMENTAL";

[[noreturn]] void nh5_fail(const std::string& rMessage) {
    throw ReadError("MSC Nastran HDF5: " + rMessage);
}

bool nh5_has_member(const std::vector<h5::CompoundMember>& rMembers, const std::string& rName) {
    for (const h5::CompoundMember& m : rMembers)
        if (m.mName == rName)
            return true;
    return false;
}

bool nh5_is_dataset(hid_t loc, const std::string& rName) {
    h5::SilenceErrors quiet;
    h5::Hid obj(H5Oopen(loc, rName.c_str(), H5P_DEFAULT), H5Oclose);
    return obj.Valid() && H5Iget_type(obj) == H5I_DATASET;
}

bool nh5_is_group(hid_t loc, const std::string& rName) {
    h5::SilenceErrors quiet;
    h5::Hid obj(H5Oopen(loc, rName.c_str(), H5P_DEFAULT), H5Oclose);
    return obj.Valid() && H5Iget_type(obj) == H5I_GROUP;
}

/** Sorted names of the datasets (or groups) directly under a path; none if it is absent. */
std::vector<std::string> nh5_children(hid_t file, const std::string& rPath, bool Groups) {
    std::vector<std::string> out;
    if (!h5::exists(file, rPath) || !nh5_is_group(file, rPath))
        return out;
    const h5::Hid g = h5::open_group(file, rPath);
    for (const std::string& name : h5::link_names(g)) {
        const std::string full = rPath + "/" + name;
        if (Groups ? nh5_is_group(file, full) : nh5_is_dataset(file, full))
            out.push_back(name);
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::size_t nh5_rows(hid_t file, const std::string& rPath) {
    return h5::dataset_num_rows(file, rPath);
}

/** An integer member, all rows, as int64 (shape kept). */
std::vector<std::int64_t> nh5_int_member(hid_t file, const std::string& rPath,
                                         const std::string& rMember, std::size_t Row0,
                                         std::size_t Count, std::size_t* pWidth = nullptr) {
    const DType as = DType::Int64;
    const NDArray a = h5::read_compound_member(file, rPath, rMember, Row0, Count, &as);
    const std::int64_t* p = a.As<std::int64_t>();
    if (pWidth != nullptr) {
        std::size_t width = 1;
        for (std::size_t k = 1; k < a.Shape().size(); ++k)
            width *= a.Shape()[k];
        *pWidth = width;
    }
    return std::vector<std::int64_t>(p, p + a.Size());
}

/** A float member as double, reduced to its first entry per row when it is an array. */
std::vector<double> nh5_float_member(hid_t file, const std::string& rPath,
                                     const std::string& rMember, std::size_t Row0,
                                     std::size_t Count) {
    const DType as = DType::Float64;
    const NDArray a = h5::read_compound_member(file, rPath, rMember, Row0, Count, &as);
    const double* p = a.As<double>();
    const std::size_t width = Count == 0 ? 1 : a.Size() / Count;
    std::vector<double> out(Count);
    for (std::size_t i = 0; i < Count; ++i)
        out[i] = p[i * width];
    return out;
}

/** A float member as double, every entry: (Count * width) values and the width. */
std::vector<double> nh5_float_member_all(hid_t file, const std::string& rPath,
                                         const std::string& rMember, std::size_t Row0,
                                         std::size_t Count, std::size_t& rWidth) {
    const DType as = DType::Float64;
    const NDArray a = h5::read_compound_member(file, rPath, rMember, Row0, Count, &as);
    rWidth = Count == 0 ? 1 : a.Size() / Count;
    const double* p = a.As<double>();
    return std::vector<double>(p, p + a.Size());
}

/** A fixed-length string member, trailing blanks and NULs removed. */
std::vector<std::string> nh5_string_member(hid_t file, const std::string& rPath,
                                           const std::string& rMember, std::size_t Row0,
                                           std::size_t Count) {
    h5::Hid d(H5Dopen2(file, rPath.c_str(), H5P_DEFAULT), H5Dclose);
    // This read is a raw H5Dread, so the chunks it reaches are verified here, as hdf5_util does.
    h5::verify_chunk_rows(d, rPath, Row0, Count);
    h5::Hid ftype(H5Dget_type(d), H5Tclose);
    const int index = H5Tget_member_index(ftype, rMember.c_str());
    if (index < 0)
        nh5_fail(rPath + " has no " + rMember + " column");
    h5::Hid mtype(H5Tget_member_type(ftype, static_cast<unsigned>(index)), H5Tclose);
    if (H5Tget_class(mtype) != H5T_STRING || H5Tis_variable_str(mtype) > 0)
        nh5_fail(rPath + " " + rMember + " is not a fixed-length string");
    const std::size_t len = H5Tget_size(mtype);
    std::vector<std::string> out(Count);
    if (Count == 0)
        return out;
    h5::Hid str(H5Tcopy(H5T_C_S1), H5Tclose);
    H5Tset_size(str, len);
    H5Tset_strpad(str, H5T_STR_NULLPAD);
    h5::Hid mem_type(H5Tcreate(H5T_COMPOUND, len), H5Tclose);
    if (!mem_type.Valid() || H5Tinsert(mem_type, rMember.c_str(), 0, str) < 0)
        nh5_fail("could not build a memory type for " + rPath + " " + rMember);
    h5::Hid space(H5Dget_space(d), H5Sclose);
    const hsize_t start = Row0;
    const hsize_t count = Count;
    H5Sselect_hyperslab(space, H5S_SELECT_SET, &start, nullptr, &count, nullptr);
    h5::Hid mem(H5Screate_simple(1, &count, nullptr), H5Sclose);
    std::vector<char> buf(len * Count);
    if (H5Dread(d, mem_type, mem, space, H5P_DEFAULT, buf.data()) < 0)
        nh5_fail("failed reading " + rPath + " " + rMember);
    for (std::size_t i = 0; i < Count; ++i) {
        std::string v(buf.data() + i * len, len);
        while (!v.empty() && (v.back() == ' ' || v.back() == '\0'))
            v.pop_back();
        out[i] = std::move(v);
    }
    return out;
}

bool nh5_is_float(const h5::CompoundMember& rM) {
    return rM.mNumeric && (rM.mDtype == DType::Float32 || rM.mDtype == DType::Float64);
}

/** One row of /NASTRAN/RESULT/DOMAINS. */
struct Nh5Domain {
    std::int64_t mId = 0;
    std::int64_t mSubcase = 0;
    std::int64_t mStep = 0;
    std::int64_t mAnalysis = 0;
    double mTimeFreqEigr = 0.0;
    double mEigi = 0.0;
    std::int64_t mMode = 0;
};

/** One result table and its INDEX: domain id -> (position, length). */
struct Nh5ResultTable {
    std::string mPath;   // /NASTRAN/RESULT/NODAL/<T> or .../ELEMENTAL/<G>/<T>
    std::string mGroup;  // "" for a nodal table, <G> for an elemental one
    std::string mName;   // <T>
    bool mNodal = true;
    std::string mKey;  // "ID" or "EID"
    std::vector<h5::CompoundMember> mMembers;
    std::map<std::int64_t, std::pair<std::size_t, std::size_t>> mIndex;
};

/** The file's structure: everything but the result payloads. */
struct Nh5File {
    h5::Hid mFile;
    std::vector<std::int64_t> mGridIds;
    std::unordered_map<std::int64_t, std::size_t> mGridIndex;
    std::unordered_set<std::int64_t> mScalarPoints;  // SPOINT/EPOINT ids
    std::vector<Nh5Domain> mSteps;                   // result domains, DOMAINS-row order
    std::vector<Nh5ResultTable> mTables;

    explicit Nh5File(const std::string& rPath) {
        {
            h5::SilenceErrors quiet;
            try {
                mFile = h5::open_file_read(rPath);
            } catch (const ReadError&) {
                nh5_fail("'" + rPath + "' is not an HDF5 file that can be opened");
            }
        }
        if (!h5::exists(mFile, "/NASTRAN") || !h5::exists(mFile, "/NASTRAN/INPUT") ||
            !h5::exists(mFile, "/NASTRAN/INPUT/NODE") || !nh5_is_dataset(mFile, kNh5Grid))
            nh5_fail("'" + rPath +
                     "' has no /NASTRAN/INPUT/NODE/GRID table; it is not an MSC Nastran HDF5 "
                     "result file");
        {
            const h5::Hid nastran = h5::open_group(mFile, "/NASTRAN");
            if (h5::has_attr(nastran, "VERSION")) {
                std::string version = h5::read_attr_string(nastran, "VERSION");
                std::string lower = version;
                std::transform(lower.begin(), lower.end(), lower.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                const auto first = lower.find_first_not_of(" \t");
                if (first == std::string::npos || lower.compare(first, 3, "msc") != 0)
                    nh5_fail("'" + rPath + "' was written by '" + version +
                             "', not MSC Nastran; other vendors' HDF5 schemas are not supported");
            }
        }
        ReadGrids();
        ReadResultTables();
    }

    void ReadGrids() {
        const std::size_t n = nh5_rows(mFile, kNh5Grid);
        mGridIds = nh5_int_member(mFile, kNh5Grid, "ID", 0, n);
        mGridIndex.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
            if (!mGridIndex.emplace(mGridIds[i], i).second)
                nh5_fail("GRID " + std::to_string(mGridIds[i]) + " is defined twice");
        for (const char* table : {"/NASTRAN/INPUT/NODE/SPOINT", "/NASTRAN/INPUT/NODE/EPOINT"}) {
            if (!nh5_is_dataset(mFile, table) ||
                !nh5_has_member(h5::compound_members(mFile, table), "ID"))
                continue;
            for (std::int64_t id : nh5_int_member(mFile, table, "ID", 0, nh5_rows(mFile, table)))
                mScalarPoints.insert(id);
        }
    }

    void AddTable(const std::string& rPath, const std::string& rGroup, const std::string& rName,
                  bool Nodal) {
        Nh5ResultTable t;
        t.mPath = rPath;
        t.mGroup = rGroup;
        t.mName = rName;
        t.mNodal = Nodal;
        t.mMembers = h5::compound_members(mFile, rPath);
        // A table without DOMAIN_ID is a lookup table (ENERGY/IDENT), not a result.
        if (!nh5_has_member(t.mMembers, "DOMAIN_ID"))
            return;
        if (!Nodal && nh5_has_member(t.mMembers, "EID"))
            t.mKey = "EID";
        else if (nh5_has_member(t.mMembers, "ID"))
            t.mKey = "ID";
        else {
            log::warn("MSC Nastran HDF5: result table {} has no ID/EID column; skipped", rPath);
            return;
        }
        const std::string index = "/INDEX" + rPath;
        if (!nh5_is_dataset(mFile, index)) {
            log::warn("MSC Nastran HDF5: result table {} has no {} table; skipped", rPath, index);
            return;
        }
        const std::size_t n = nh5_rows(mFile, index);
        const auto dom = nh5_int_member(mFile, index, "DOMAIN_ID", 0, n);
        const auto pos = nh5_int_member(mFile, index, "POSITION", 0, n);
        const auto len = nh5_int_member(mFile, index, "LENGTH", 0, n);
        const std::size_t rows = nh5_rows(mFile, rPath);
        for (std::size_t i = 0; i < n; ++i) {
            if (pos[i] < 0 || len[i] < 0 ||
                static_cast<std::size_t>(pos[i]) + static_cast<std::size_t>(len[i]) > rows)
                nh5_fail(index + " gives rows outside " + rPath);
            t.mIndex[dom[i]] = {static_cast<std::size_t>(pos[i]), static_cast<std::size_t>(len[i])};
        }
        mTables.push_back(std::move(t));
    }

    void ReadResultTables() {
        for (const std::string& name : nh5_children(mFile, kNh5Nodal, false))
            AddTable(std::string(kNh5Nodal) + "/" + name, "", name, true);
        for (const std::string& group : nh5_children(mFile, kNh5Elemental, true)) {
            const std::string gpath = std::string(kNh5Elemental) + "/" + group;
            for (const std::string& name : nh5_children(mFile, gpath, false))
                AddTable(gpath + "/" + name, group, name, false);
        }
        if (mTables.empty())
            return;
        if (!nh5_is_dataset(mFile, kNh5Domains))
            nh5_fail("the file has result tables but no /NASTRAN/RESULT/DOMAINS table");
        const std::size_t n = nh5_rows(mFile, kNh5Domains);
        const auto members = h5::compound_members(mFile, kNh5Domains);
        auto ints = [&](const char* pName) {
            return nh5_has_member(members, pName) ? nh5_int_member(mFile, kNh5Domains, pName, 0, n)
                                                  : std::vector<std::int64_t>(n, 0);
        };
        auto reals = [&](const char* pName) {
            return nh5_has_member(members, pName)
                       ? nh5_float_member(mFile, kNh5Domains, pName, 0, n)
                       : std::vector<double>(n, 0.0);
        };
        if (!nh5_has_member(members, "ID"))
            nh5_fail("/NASTRAN/RESULT/DOMAINS has no ID column");
        const auto id = ints("ID");
        const auto subcase = ints("SUBCASE");
        const auto step = ints("STEP");
        const auto analysis = ints("ANALYSIS");
        const auto mode = ints("MODE");
        const auto tfe = reals("TIME_FREQ_EIGR");
        const auto eigi = reals("EIGI");
        std::set<std::int64_t> used;
        for (const Nh5ResultTable& t : mTables)
            for (const auto& kv : t.mIndex)
                used.insert(kv.first);
        std::set<std::int64_t> known;
        for (std::size_t i = 0; i < n; ++i) {
            known.insert(id[i]);
            if (used.count(id[i]) == 0)
                continue;
            Nh5Domain d;
            d.mId = id[i];
            d.mSubcase = subcase[i];
            d.mStep = step[i];
            d.mAnalysis = analysis[i];
            d.mTimeFreqEigr = tfe[i];
            d.mEigi = eigi[i];
            d.mMode = mode[i];
            mSteps.push_back(d);
            used.erase(id[i]);  // a repeated DOMAINS row is one step
        }
        for (std::int64_t missing : used)
            if (known.count(missing) == 0)
                nh5_fail("a result table names domain " + std::to_string(missing) +
                         ", which /NASTRAN/RESULT/DOMAINS does not define");
    }
};

NDArray nh5_scalar(DType Type, double Value) {
    NDArray out(Type, {std::size_t{1}});
    if (Type == DType::Float64)
        *out.As<double>() = Value;
    else
        *out.As<std::int64_t>() = static_cast<std::int64_t>(Value);
    return out;
}

NDArray nh5_int_array(const std::vector<std::int64_t>& rV) {
    NDArray out(DType::Int64, {rV.size()});
    std::copy(rV.begin(), rV.end(), out.As<std::int64_t>());
    return out;
}

/** The nodal arrays a table yields: output name -> member names (1 or 3). */
std::vector<std::pair<std::string, std::vector<std::string>>> nh5_nodal_outputs(
    const Nh5ResultTable& rT) {
    std::vector<std::pair<std::string, std::vector<std::string>>> out;
    std::set<std::string> taken{rT.mKey, "DOMAIN_ID"};
    std::set<std::string> floats;
    for (const h5::CompoundMember& m : rT.mMembers)
        if (nh5_is_float(m) && taken.count(m.mName) == 0)
            floats.insert(m.mName);
    const std::string suffix = "_CPLX";
    const bool cplx = rT.mName.size() > suffix.size() &&
                      rT.mName.compare(rT.mName.size() - suffix.size(), suffix.size(), suffix) == 0;
    const std::string base = cplx ? rT.mName.substr(0, rT.mName.size() - suffix.size()) : rT.mName;
    auto group = [&](const std::string& rName, std::vector<std::string> members) {
        for (const std::string& m : members)
            if (floats.count(m) == 0)
                return;
        for (const std::string& m : members)
            floats.erase(m);
        out.emplace_back(rName, std::move(members));
    };
    if (cplx) {
        group(base + "_real", {"XR", "YR", "ZR"});
        group(base + "_imag", {"XI", "YI", "ZI"});
        group(base + "_ROT_real", {"RXR", "RYR", "RZR"});
        group(base + "_ROT_imag", {"RXI", "RYI", "RZI"});
    } else {
        group(base, {"X", "Y", "Z"});
        group(base + "_ROT", {"RX", "RY", "RZ"});
    }
    if (floats.size() == 1 && *floats.begin() == "VALUE") {
        out.emplace_back(rT.mName, std::vector<std::string>{"VALUE"});
        floats.clear();
    }
    // Remaining members in declaration order.
    for (const h5::CompoundMember& m : rT.mMembers)
        if (floats.count(m.mName) != 0)
            out.emplace_back(rT.mName + ":" + m.mName, std::vector<std::string>{m.mName});
    return out;
}

/** The CORD1R/C/S and CORD2R/C/S tables of /NASTRAN/INPUT/COORDINATE_SYSTEM. */
std::vector<detail::NastranCoordCard> nh5_coord_cards(hid_t File) {
    constexpr const char* kDir = "/NASTRAN/INPUT/COORDINATE_SYSTEM";
    std::vector<detail::NastranCoordCard> out;
    static const std::pair<const char*, int> kTables[] = {
        {"CORD2R", 1}, {"CORD2C", 2}, {"CORD2S", 3}, {"CORD1R", 1}, {"CORD1C", 2}, {"CORD1S", 3}};
    for (const auto& [name, type] : kTables) {
        const std::string path = std::string(kDir) + "/" + name;
        if (!nh5_is_dataset(File, path))
            continue;
        const auto members = h5::compound_members(File, path);
        const std::size_t n = nh5_rows(File, path);
        const bool by_grids = name[4] == '1';
        const std::vector<const char*> needed =
            by_grids ? std::vector<const char*>{"CID", "G1", "G2", "G3"}
                     : std::vector<const char*>{"CID", "RID", "A1", "A2", "A3", "B1",
                                                "B2",  "B3",  "C1", "C2", "C3"};
        bool complete = true;
        for (const char* m : needed)
            complete = complete && nh5_has_member(members, m);
        if (!complete) {
            log::warn("MSC Nastran HDF5: {} lacks the expected columns; its systems are not read",
                      path);
            continue;
        }
        const auto cid = nh5_int_member(File, path, "CID", 0, n);
        std::vector<std::vector<std::int64_t>> ints;
        std::vector<std::vector<double>> reals;
        for (std::size_t k = 1; k < needed.size(); ++k) {
            if (by_grids || k == 1)
                ints.push_back(nh5_int_member(File, path, needed[k], 0, n));
            else
                reals.push_back(nh5_float_member(File, path, needed[k], 0, n));
        }
        for (std::size_t i = 0; i < n; ++i) {
            detail::NastranCoordCard c;
            c.mCid = cid[i];
            c.mType = type;
            c.mByGrids = by_grids;
            if (by_grids) {
                for (int k = 0; k < 3; ++k)
                    c.mGrids[k] = ints[static_cast<std::size_t>(k)][i];
            } else {
                c.mRid = ints[0][i];
                for (std::size_t k = 0; k < 9; ++k)
                    c.mAbc[k] = reals[k][i];
            }
            out.push_back(c);
        }
    }
    return out;
}

std::string nh5_join(const std::vector<std::string>& rV) {
    std::string out;
    for (const std::string& s : rV)
        out += (out.empty() ? "" : ", ") + s;
    return out;
}

}  // namespace

Mesh read_nastran_h5(const std::string& rPath, const ReadOptions& rOpts) {
    detail::LibraryLock lock;
    const Nh5File file(rPath);
    const hid_t f = file.mFile;
    Mesh mesh;

    // --- points ---------------------------------------------------------------
    const std::size_t npts = file.mGridIds.size();
    detail::NastranCoordSystems systems;
    std::vector<std::int64_t> cd;
    {
        const DType as = DType::Float64;
        NDArray x = h5::read_compound_member(f, kNh5Grid, "X", 0, npts, &as);
        if (x.Shape().size() != 2 || x.Shape()[1] != 3)
            nh5_fail("GRID X is not a 3-vector");
        mesh.AssignPoints(std::move(x));
        const auto members = h5::compound_members(f, kNh5Grid);
        auto frame = [&](const char* pName) {
            return nh5_has_member(members, pName) ? nh5_int_member(f, kNh5Grid, pName, 0, npts)
                                                  : std::vector<std::int64_t>(npts, 0);
        };
        cd = frame("CD");
        systems = detail::nastran_apply_frames(mesh, nh5_coord_cards(f), file.mGridIds, frame("CP"),
                                               cd, "MSC Nastran HDF5");
    }
    const NDArray basic_points = mesh.Points();
    const double* basic = basic_points.As<double>();

    // --- cells and property regions ------------------------------------------------
    std::vector<detail::NastranCardRows> cards;
    std::vector<std::string> skipped_cards;
    for (const std::string& card : nh5_children(f, kNh5Elements, false)) {
        const detail::NastranCardSpec* spec = detail::nastran_card_spec(card);
        if (spec == nullptr) {
            skipped_cards.push_back(card);
            continue;
        }
        const std::string path = std::string(kNh5Elements) + "/" + card;
        const auto members = h5::compound_members(f, path);
        const std::size_t n = nh5_rows(f, path);
        if (!nh5_has_member(members, "EID"))
            nh5_fail(path + " has no EID column");
        detail::NastranCardRows rows;
        rows.mCard = card;
        rows.mEid = nh5_int_member(f, path, "EID", 0, n);
        rows.mPid = nh5_has_member(members, "PID") ? nh5_int_member(f, path, "PID", 0, n)
                                                   : std::vector<std::int64_t>(n, -1);
        if (nh5_has_member(members, "G")) {
            rows.mNodes = nh5_int_member(f, path, "G", 0, n, &rows.mWidth);
        } else {
            const char* a = nh5_has_member(members, "GA") ? "GA" : "G1";
            const char* b = nh5_has_member(members, "GA") ? "GB" : "G2";
            if (!nh5_has_member(members, a) || !nh5_has_member(members, b))
                nh5_fail(path + " has no G, GA/GB or G1/G2 columns");
            const auto ga = nh5_int_member(f, path, a, 0, n);
            const auto gb = nh5_int_member(f, path, b, 0, n);
            rows.mWidth = 2;
            rows.mNodes.resize(2 * n);
            for (std::size_t i = 0; i < n; ++i) {
                rows.mNodes[2 * i] = ga[i];
                rows.mNodes[2 * i + 1] = gb[i];
            }
        }
        if (rows.mWidth < spec->mLinearNodes)
            nh5_fail(path + " has " + std::to_string(rows.mWidth) + " node columns, " + card +
                     " needs " + std::to_string(spec->mLinearNodes));
        cards.push_back(std::move(rows));
    }
    if (!skipped_cards.empty())
        log::warn("MSC Nastran HDF5: skipped element tables with no cell type: {}",
                  nh5_join(skipped_cards));
    std::map<std::int64_t, std::string> ptype;
    for (const std::string& prop : nh5_children(f, kNh5Properties, false)) {
        const std::string path = std::string(kNh5Properties) + "/" + prop;
        if (!nh5_has_member(h5::compound_members(f, path), "PID"))
            continue;
        for (std::int64_t p : nh5_int_member(f, path, "PID", 0, nh5_rows(f, path)))
            ptype.emplace(p, prop);
    }
    // Grouped properties (PCOMP/IDENTITY) live one level deeper.
    for (const std::string& prop : nh5_children(f, kNh5Properties, true)) {
        const std::string path = std::string(kNh5Properties) + "/" + prop + "/IDENTITY";
        if (!nh5_is_dataset(f, path) || !nh5_has_member(h5::compound_members(f, path), "PID"))
            continue;
        for (std::int64_t p : nh5_int_member(f, path, "PID", 0, nh5_rows(f, path)))
            ptype.emplace(p, prop);
    }
    const detail::NastranCells model = detail::nastran_add_cells(
        mesh, cards, file.mGridIndex, file.mScalarPoints, ptype, "MSC Nastran HDF5");
    const auto& cell_index = model.mCellIndex;
    const std::size_t ncells = model.mNumCells;

    // --- the step -------------------------------------------------------------
    if (file.mSteps.empty()) {
        rOpts.ResolveTimeStep(0);  // refuses anything but the first/last "step" of no steps
        return mesh;
    }
    const Nh5Domain& dom = file.mSteps[rOpts.ResolveTimeStep(file.mSteps.size())];
    mesh.AddFieldData(kSequenceTimeKey, nh5_scalar(DType::Float64, dom.mTimeFreqEigr));
    mesh.AddFieldData("nastran:domain", nh5_scalar(DType::Int64, static_cast<double>(dom.mId)));
    mesh.AddFieldData("nastran:subcase",
                      nh5_scalar(DType::Int64, static_cast<double>(dom.mSubcase)));
    mesh.AddFieldData("nastran:step", nh5_scalar(DType::Int64, static_cast<double>(dom.mStep)));
    mesh.AddFieldData("nastran:analysis",
                      nh5_scalar(DType::Int64, static_cast<double>(dom.mAnalysis)));
    mesh.AddFieldData("nastran:mode", nh5_scalar(DType::Int64, static_cast<double>(dom.mMode)));
    mesh.AddFieldData("nastran:eigi", nh5_scalar(DType::Float64, dom.mEigi));
    if (!rOpts.WantsAnyData())
        return mesh;

    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::string> used_point_names;
    // Elemental arrays accumulate across tables: <G>:<M> -> per-global-cell values.
    std::map<std::string, std::vector<double>> cell_arrays;
    std::vector<std::string> cell_order;
    // Multi-valued cell arrays (per ply, per station, per corner or element
    // node) gather (cell, column, value) and are laid out once every table is
    // read, as (cells, columns) with NaN where a cell has no value.
    struct Wide {
        std::vector<std::size_t> mCell, mCol;
        std::vector<double> mValue;
    };
    std::map<std::string, Wide> wide;
    auto push = [&](const std::string& rName, std::size_t Cell, std::size_t Col, double V) {
        Wide& w = wide[rName];
        w.mCell.push_back(Cell);
        w.mCol.push_back(Col);
        w.mValue.push_back(V);
    };
    // The position of GRID `Grid` in the connectivity of global cell `Cell`.
    constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
    auto node_position = [&](std::size_t Cell, std::int64_t Grid) -> std::size_t {
        const auto g = file.mGridIndex.find(Grid);
        if (g == file.mGridIndex.end())
            return npos;
        const std::size_t b = static_cast<std::size_t>(
            std::upper_bound(model.mOffsets.begin(), model.mOffsets.end(), Cell) -
            model.mOffsets.begin() - 1);
        const NDArray& conn = mesh.Cells(b).Conn();
        const std::size_t width = conn.Shape()[1];
        const std::size_t row = Cell - model.mOffsets[b];
        for (std::size_t k = 0; k < width; ++k)
            if (static_cast<std::size_t>(detail::read_int(conn, row * width + k)) == g->second)
                return k;
        return npos;
    };
    auto has_int_member = [](const Nh5ResultTable& rT, const char* pName, std::size_t Dims) {
        for (const h5::CompoundMember& m : rT.mMembers)
            if (m.mName == pName && m.mNumeric && m.mDtype != DType::Float32 &&
                m.mDtype != DType::Float64 && m.mDims.size() == Dims)
                return true;
        return false;
    };
    auto starts = [](const std::string& rName, const char* pPrefix) {
        return rName.rfind(pPrefix, 0) == 0;
    };

    // The tables with several values per element or node: false when `rT` is not one.
    auto read_multi = [&](const Nh5ResultTable& rT, std::size_t Row0, std::size_t Count) -> bool {
        const bool ply = !rT.mNodal && has_int_member(rT, "PLY", 0);
        const bool bars = !rT.mNodal && (rT.mName == "BARS" || rT.mName == "BARS_CPLX");
        const bool arrays = !rT.mNodal && has_int_member(rT, "GRID", 1);
        const bool grid_force = rT.mNodal && rT.mName == "GRID_FORCE" &&
                                has_int_member(rT, "EID", 0) &&
                                nh5_has_member(rT.mMembers, "ELNAME");
        if (!ply && !bars && !arrays && !grid_force)
            return false;
        std::vector<std::string> floats;
        for (const h5::CompoundMember& m : rT.mMembers)
            if (nh5_is_float(m) && m.mName != rT.mKey && m.mName != "DOMAIN_ID")
                floats.push_back(m.mName);
        const auto keys = nh5_int_member(f, rT.mPath, rT.mKey, Row0, Count);

        if (grid_force) {
            // Element rows (EID > 0) are the forces on each element node, as
            // (cells, nodes) in the cell's node order; the others (*TOTALS*,
            // APP-LOAD, F-OF-SPC, ...) are point data named after ELNAME.
            // Both are in the GRID's output system (CD), rotated to basic.
            const auto eids = nh5_int_member(f, rT.mPath, "EID", Row0, Count);
            const auto elname = nh5_string_member(f, rT.mPath, "ELNAME", Row0, Count);
            std::vector<std::vector<double>> v;
            v.reserve(floats.size());
            for (const std::string& m : floats)
                v.push_back(nh5_float_member(f, rT.mPath, m, Row0, Count));
            const bool triplets = floats.size() == 6;
            std::map<std::string, std::vector<double>> totals;
            std::vector<std::string> total_order;
            for (std::size_t r = 0; r < Count; ++r) {
                const auto g = file.mGridIndex.find(keys[r]);
                if (g == file.mGridIndex.end())
                    continue;
                const std::size_t p = g->second;
                std::vector<double> row(floats.size());
                for (std::size_t c = 0; c < floats.size(); ++c)
                    row[c] = v[c][r];
                if (triplets)
                    detail::nastran_rotate_to_basic(systems, cd[p], basic + 3 * p, row.data(), 2);
                if (eids[r] > 0) {
                    const auto c = cell_index.find(eids[r]);
                    if (c == cell_index.end())
                        continue;
                    const std::size_t pos = node_position(c->second, keys[r]);
                    if (pos == npos)
                        continue;
                    for (std::size_t k = 0; k < floats.size(); ++k) {
                        const std::string name = rT.mName + ":" + floats[k];
                        if (rOpts.WantsArray(name))
                            push(name, c->second, pos, row[k]);
                    }
                    continue;
                }
                std::string label;
                for (char ch : elname[r])
                    if (ch != ' ' && ch != '*')
                        label += ch;
                for (std::size_t k = 0; k < floats.size(); ++k) {
                    const std::string name = rT.mName + ":" + label + ":" + floats[k];
                    if (!rOpts.WantsArray(name))
                        continue;
                    auto& values = totals[name];
                    if (values.empty()) {
                        values.assign(npts, std::numeric_limits<double>::quiet_NaN());
                        total_order.push_back(name);
                    }
                    values[p] = row[k];
                }
            }
            for (const std::string& name : total_order) {
                NDArray a(DType::Float64, {npts});
                std::copy(totals[name].begin(), totals[name].end(), a.As<double>());
                mesh.AddPointData(name, std::move(a));
            }
            return true;
        }

        std::vector<std::size_t> target(Count, npos);
        for (std::size_t r = 0; r < Count; ++r) {
            const auto c = cell_index.find(keys[r]);
            if (c != cell_index.end())
                target[r] = c->second;
        }
        if (ply || bars) {
            // One row per ply (column PLY - 1) or per station along a bar
            // (columns in row order).
            std::vector<std::size_t> column(Count, 0);
            if (ply) {
                const auto plies = nh5_int_member(f, rT.mPath, "PLY", Row0, Count);
                for (std::size_t r = 0; r < Count; ++r)
                    column[r] = plies[r] >= 1 ? static_cast<std::size_t>(plies[r] - 1) : npos;
            } else {
                std::unordered_map<std::int64_t, std::size_t> seen;
                for (std::size_t r = 0; r < Count; ++r)
                    column[r] = seen[keys[r]]++;
            }
            const std::string suffix = ply ? "@ply" : "@station";
            for (const std::string& m : floats) {
                const std::string name = rT.mGroup + ":" + m + suffix;
                if (!rOpts.WantsArray(name))
                    continue;
                const auto v = nh5_float_member(f, rT.mPath, m, Row0, Count);
                for (std::size_t r = 0; r < Count; ++r)
                    if (target[r] != npos && column[r] != npos)
                        push(name, target[r], column[r], v[r]);
            }
            return true;
        }

        // One row per element with arrays: entry 0 is the centre (or end A of
        // a beam), which stays the plain `<G>:<M>` value. A BEAM's entries are
        // its 11 stations; the others' entries 1.. are corners, placed at their
        // GRID's position in the cell.
        {
            std::unordered_set<std::int64_t> once;
            for (std::int64_t k : keys)
                if (!once.insert(k).second) {
                    log::warn("MSC Nastran HDF5: {} has several rows per element; skipped",
                              rT.mPath);
                    return true;
                }
        }
        const bool beam = starts(rT.mName, "BEAM");
        std::size_t gw = 1;
        const auto grids = nh5_int_member(f, rT.mPath, "GRID", Row0, Count, &gw);
        // A beam station with no GRID and no distance (SD) was not output.
        std::size_t sw = 0;
        const auto sd = beam && nh5_has_member(rT.mMembers, "SD")
                            ? nh5_float_member_all(f, rT.mPath, "SD", Row0, Count, sw)
                            : std::vector<double>();
        for (const std::string& m : floats) {
            const std::string centre = rT.mGroup + ":" + m;
            const std::string name = centre + (beam ? "@station" : "@corner");
            const bool want_centre = rOpts.WantsArray(centre);
            const bool want_wide = rOpts.WantsArray(name);
            if (!want_centre && !want_wide)
                continue;
            std::size_t w = 1;
            const auto v = nh5_float_member_all(f, rT.mPath, m, Row0, Count, w);
            if (want_centre) {
                auto& values = cell_arrays[centre];
                if (values.empty()) {
                    values.assign(ncells, std::numeric_limits<double>::quiet_NaN());
                    cell_order.push_back(centre);
                }
                for (std::size_t r = 0; r < Count; ++r)
                    if (target[r] != npos)
                        values[target[r]] = v[r * w];
            }
            if (!want_wide || w < 2)
                continue;
            for (std::size_t r = 0; r < Count; ++r) {
                if (target[r] == npos)
                    continue;
                for (std::size_t k = beam ? 0 : 1; k < w; ++k) {
                    std::size_t col = k;
                    if (beam && k > 0 && w == gw && sw == w && grids[r * gw + k] == 0 &&
                        sd[r * sw + k] == 0.0)
                        continue;
                    if (!beam) {
                        if (w != gw || grids[r * gw + k] <= 0)
                            continue;
                        col = node_position(target[r], grids[r * gw + k]);
                        if (col == npos)
                            continue;
                    }
                    push(name, target[r], col, v[r * w + k]);
                }
            }
        }
        return true;
    };

    for (const Nh5ResultTable& t : file.mTables) {
        const auto it = t.mIndex.find(dom.mId);
        if (it == t.mIndex.end())
            continue;
        const auto [row0, count] = it->second;
        if (read_multi(t, row0, count))
            continue;

        // Which outputs of this table are wanted, before any payload is read.
        std::vector<std::pair<std::string, std::vector<std::string>>> outputs;
        if (t.mNodal) {
            for (auto& o : nh5_nodal_outputs(t)) {
                std::string name = o.first;
                for (int k = 2; std::find(used_point_names.begin(), used_point_names.end(), name) !=
                                used_point_names.end();
                     ++k)
                    name = o.first + "_" + std::to_string(k);
                if (rOpts.WantsArray(name))
                    outputs.emplace_back(name, std::move(o.second));
            }
        } else {
            for (const h5::CompoundMember& m : t.mMembers)
                if (nh5_is_float(m) && m.mName != t.mKey && m.mName != "DOMAIN_ID" &&
                    rOpts.WantsArray(t.mGroup + ":" + m.mName))
                    outputs.emplace_back(t.mGroup + ":" + m.mName,
                                         std::vector<std::string>{m.mName});
        }
        if (outputs.empty())
            continue;

        const auto keys = nh5_int_member(f, t.mPath, t.mKey, row0, count);
        {
            std::unordered_set<std::int64_t> seen;
            bool repeated = false;
            for (std::int64_t k : keys)
                repeated = repeated || !seen.insert(k).second;
            if (repeated) {
                log::warn("MSC Nastran HDF5: {} has several rows per {}; skipped", t.mPath,
                          t.mNodal ? "node" : "element");
                continue;
            }
        }
        // Row -> target index (point or global cell), or npos when not in the mesh.
        std::vector<std::size_t> target(count, npos);
        std::size_t mapped = 0;
        for (std::size_t r = 0; r < count; ++r) {
            if (t.mNodal) {
                const auto g = file.mGridIndex.find(keys[r]);
                if (g != file.mGridIndex.end())
                    target[r] = g->second;
            } else {
                const auto c = cell_index.find(keys[r]);
                if (c != cell_index.end())
                    target[r] = c->second;
            }
            mapped += target[r] != npos ? 1 : 0;
        }
        if (mapped == 0)
            continue;

        for (auto& [name, members] : outputs) {
            const std::size_t nc = members.size();
            if (t.mNodal) {
                NDArray data =
                    nc == 1 ? NDArray(DType::Float64, {npts}) : NDArray(DType::Float64, {npts, nc});
                double* out = data.As<double>();
                std::fill(out, out + data.Size(), nan);
                for (std::size_t c = 0; c < nc; ++c) {
                    const auto v = nh5_float_member(f, t.mPath, members[c], row0, count);
                    for (std::size_t r = 0; r < count; ++r)
                        if (target[r] != npos)
                            out[target[r] * nc + c] = v[r];
                }
                // Vector results are in each GRID's output system (CD).
                if (nc == 3)
                    for (std::size_t r = 0; r < count; ++r)
                        if (target[r] != npos)
                            detail::nastran_rotate_to_basic(systems, cd[target[r]],
                                                            basic + 3 * target[r],
                                                            out + 3 * target[r], 1);
                used_point_names.push_back(name);
                mesh.AddPointData(name, std::move(data));
            } else {
                auto& values = cell_arrays[name];
                if (values.empty()) {
                    values.assign(ncells, nan);
                    cell_order.push_back(name);
                }
                const auto v = nh5_float_member(f, t.mPath, members[0], row0, count);
                for (std::size_t r = 0; r < count; ++r)
                    if (target[r] != npos)
                        values[target[r]] = v[r];
            }
        }
    }
    for (const std::string& name : cell_order) {
        const std::vector<double>& values = cell_arrays[name];
        std::vector<NDArray> per_block;
        for (std::size_t b = 0; b < model.mSizes.size(); ++b) {
            NDArray a(DType::Float64, {model.mSizes[b]});
            std::copy(
                values.begin() + static_cast<std::ptrdiff_t>(model.mOffsets[b]),
                values.begin() + static_cast<std::ptrdiff_t>(model.mOffsets[b] + model.mSizes[b]),
                a.As<double>());
            per_block.push_back(std::move(a));
        }
        mesh.AddCellData(name, std::move(per_block));
    }
    for (const auto& [name, w] : wide) {
        std::size_t width = 0;
        for (std::size_t col : w.mCol)
            width = std::max(width, col + 1);
        std::vector<NDArray> per_block;
        for (std::size_t b = 0; b < model.mSizes.size(); ++b) {
            NDArray a(DType::Float64, {model.mSizes[b], width});
            std::fill(a.As<double>(), a.As<double>() + a.Size(), nan);
            per_block.push_back(std::move(a));
        }
        for (std::size_t i = 0; i < w.mCell.size(); ++i) {
            const std::size_t b = static_cast<std::size_t>(
                std::upper_bound(model.mOffsets.begin(), model.mOffsets.end(), w.mCell[i]) -
                model.mOffsets.begin() - 1);
            per_block[b].As<double>()[(w.mCell[i] - model.mOffsets[b]) * width + w.mCol[i]] =
                w.mValue[i];
        }
        mesh.AddCellData(name, std::move(per_block));
        NDArray layout(DType::Int64, {std::size_t{2}});
        layout.As<std::int64_t>()[0] = static_cast<std::int64_t>(width);
        layout.As<std::int64_t>()[1] = 1;
        mesh.AddFieldData("nastran:layout:" + name, std::move(layout));
    }
    return mesh;
}

MeshMetadata read_nastran_h5_metadata(const std::string& rPath, const ReadOptions& /*rOpts*/) {
    detail::LibraryLock lock;
    const Nh5File file(rPath);
    // No header-only path: the model is read in full, with step 0's results.
    MeshMetadata meta = metadata_from_mesh(read_nastran_h5(rPath, ReadOptions{}));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "nastran_h5";
    meta.mTimeValues.clear();
    for (const Nh5Domain& d : file.mSteps)
        meta.mTimeValues.push_back(d.mTimeFreqEigr);
    return meta;
}


// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------
//
// The twin of `nastran_h5/_write.py`: the same tables from the same rules.

namespace {

[[noreturn]] void nh5w_fail(const std::string& rMessage) {
    throw WriteError("Nastran HDF5: " + rMessage);
}

enum Nh5wProp { kNh5wNoProp = -1, kNh5wProd = 0, kNh5wPshell, kNh5wPsolid, kNh5wPtube };

constexpr int kNh5wHexa20[20] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 12, 13, 14, 15};
constexpr int kNh5wPenta15[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11};

struct Nh5wCard {
    const char* mName;
    std::size_t mLinear;  // nodes of the linear card
    std::size_t mWidth;   // columns of G
    const int* mPerm;     // quadratic cell -> Nastran's order
    Nh5wProp mProp;
    const char* mFamily;  // result tables' name; null: no results
    bool mHasPid;
};

// In card-name order.
constexpr Nh5wCard kNh5wCards[] = {
    {"CHEXA", 8, 20, kNh5wHexa20, kNh5wPsolid, "HEXA", true},
    {"CONM2", 1, 1, nullptr, kNh5wNoProp, nullptr, false},
    {"CONROD", 2, 2, nullptr, kNh5wNoProp, "CONROD", false},
    {"CPENTA", 6, 15, kNh5wPenta15, kNh5wPsolid, "PENTA", true},
    {"CPYRAM", 5, 13, nullptr, kNh5wPsolid, "PYRAM", true},
    {"CQUAD", 4, 9, nullptr, kNh5wPshell, "QUADR", true},
    {"CQUAD4", 4, 4, nullptr, kNh5wPshell, "QUAD4", true},
    {"CQUAD8", 4, 8, nullptr, kNh5wPshell, "QUAD8", true},
    {"CROD", 2, 2, nullptr, kNh5wProd, "ROD", true},
    {"CTETRA", 4, 10, nullptr, kNh5wPsolid, "TETRA", true},
    {"CTRIA3", 3, 3, nullptr, kNh5wPshell, "TRIA3", true},
    {"CTRIA6", 3, 6, nullptr, kNh5wPshell, "TRIA6", true},
    {"CTUBE", 2, 2, nullptr, kNh5wPtube, "TUBE", true},
};
constexpr std::size_t kNh5wNumCards = sizeof(kNh5wCards) / sizeof(kNh5wCards[0]);

std::size_t nh5w_card_index(std::string_view Name) {
    for (std::size_t i = 0; i < kNh5wNumCards; ++i)
        if (Name == kNh5wCards[i].mName)
            return i;
    return kNh5wNumCards;
}

struct Nh5wTypeSpec {
    const char* mType;
    const char* mCard;  // null: a line takes CROD, CTUBE or CONROD by its property
    std::size_t mNodes;
    std::int64_t mPid;  // the property id of a block with no 'nastran:pid'
};
constexpr Nh5wTypeSpec kNh5wTypes[] = {
    {"vertex", "CONM2", 1, -1},    {"line", nullptr, 2, 3},        {"triangle", "CTRIA3", 3, 1},
    {"triangle6", "CTRIA6", 6, 1}, {"quad", "CQUAD4", 4, 1},       {"quad8", "CQUAD8", 8, 1},
    {"quad9", "CQUAD", 9, 1},      {"tetra", "CTETRA", 4, 2},      {"tetra10", "CTETRA", 10, 2},
    {"pyramid", "CPYRAM", 5, 2},   {"pyramid13", "CPYRAM", 13, 2}, {"wedge", "CPENTA", 6, 2},
    {"wedge15", "CPENTA", 15, 2},  {"hexahedron", "CHEXA", 8, 2},  {"hexahedron20", "CHEXA", 20, 2},
};

const Nh5wTypeSpec* nh5w_type_spec(const std::string& rType) {
    for (const Nh5wTypeSpec& t : kNh5wTypes)
        if (rType == t.mType)
            return &t;
    return nullptr;
}

struct Nh5wBlock {
    std::size_t mIndex = 0;
    std::size_t mCount = 0;
    std::size_t mNodes = 0;
    std::vector<std::int64_t> mConn;  // 0-based, mCount x mNodes
    std::vector<std::int64_t> mEids;
    std::vector<std::int64_t> mPids;
    std::vector<std::size_t> mCards;
};

struct Nh5wRow {
    std::int64_t mEid;
    std::size_t mBlock;
    std::size_t mRow;
};

using Nh5wSignature = std::vector<std::tuple<std::string, std::size_t, std::size_t>>;

Nh5wSignature nh5w_signature(const Mesh& rMesh) {
    Nh5wSignature sig;
    sig.emplace_back("points", rMesh.NumPoints(), 0);
    for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
        const auto cb = rMesh.Cells(b);
        sig.emplace_back(cb.Type(), cb.NumCells(), cb.NodesPerCell());
    }
    return sig;
}

std::int64_t nh5w_ptube_pid(const std::string& rName) {
    static const std::string prefix = "PTUBE_";
    if (rName.size() <= prefix.size() || rName.compare(0, prefix.size(), prefix) != 0)
        return -1;
    std::int64_t pid = 0;
    for (std::size_t i = prefix.size(); i < rName.size(); ++i) {
        if (rName[i] < '0' || rName[i] > '9' || pid > 100000000000000000)
            return -1;
        pid = pid * 10 + (rName[i] - '0');
    }
    return pid;
}

/// The mesh as the file holds it: GRIDs, element cards and properties.
class Nh5wModel {
public:
    explicit Nh5wModel(const Mesh& rMesh) {
        mSignature = nh5w_signature(rMesh);
        mNpoints = rMesh.NumPoints();
        {
            const NDArray& pts = rMesh.Points();
            const std::size_t dim = rMesh.PointDim();
            const detail::DoubleView view(pts);
            mPoints.assign(mNpoints * 3, 0.0);
            for (std::size_t i = 0; i < mNpoints; ++i)
                for (std::size_t k = 0; k < 3 && k < dim; ++k)
                    mPoints[3 * i + k] = view.Data()[i * dim + k];
        }
        std::set<std::int64_t> ptube;
        for (std::size_t i = 0; i < rMesh.NumRegions(); ++i) {
            const std::int64_t pid = nh5w_ptube_pid(rMesh.Region(i).mName);
            if (pid >= 0)
                ptube.insert(pid);
        }
        const bool have_eid = rMesh.HasCellData("nastran:eid");
        const bool have_pid = rMesh.HasCellData("nastran:pid");
        std::int64_t sequential = 1;
        for (std::size_t b = 0; b < rMesh.NumCellBlocks(); ++b) {
            const auto cb = rMesh.Cells(b);
            const Nh5wTypeSpec* spec = nh5w_type_spec(cb.Type());
            if (spec == nullptr) {
                if (cb.NumCells() > 0)
                    mSkipped.insert(cb.Type());
                continue;
            }
            Nh5wBlock blk;
            blk.mIndex = b;
            blk.mCount = cb.NumCells();
            blk.mNodes = cb.NodesPerCell();
            if (blk.mCount > 0 && blk.mNodes != spec->mNodes)
                nh5w_fail("a '" + cb.Type() + "' block needs " + std::to_string(spec->mNodes) +
                          " nodes per cell");
            const detail::Int64View conn(cb.Conn());
            blk.mConn.assign(conn.Data(), conn.Data() + blk.mCount * blk.mNodes);
            for (const std::int64_t node : blk.mConn)
                if (node < 0 || static_cast<std::size_t>(node) >= mNpoints)
                    nh5w_fail("a '" + cb.Type() +
                              "' cell refers to a point the mesh does not have");
            blk.mEids.resize(blk.mCount);
            blk.mPids.resize(blk.mCount);
            if (have_eid) {
                const NDArray& a = rMesh.CellData("nastran:eid", b);
                if (a.Size() != blk.mCount)
                    nh5w_fail("'nastran:eid' must have one value per cell");
                const detail::Int64View view(a);
                blk.mEids.assign(view.Data(), view.Data() + blk.mCount);
            } else {
                for (std::size_t i = 0; i < blk.mCount; ++i)
                    blk.mEids[i] = sequential + static_cast<std::int64_t>(i);
            }
            sequential += static_cast<std::int64_t>(blk.mCount);
            if (have_pid) {
                const NDArray& a = rMesh.CellData("nastran:pid", b);
                if (a.Size() != blk.mCount)
                    nh5w_fail("'nastran:pid' must have one value per cell");
                const detail::Int64View view(a);
                blk.mPids.assign(view.Data(), view.Data() + blk.mCount);
            } else {
                blk.mPids.assign(blk.mCount, spec->mPid);
            }
            blk.mCards.resize(blk.mCount);
            for (std::size_t i = 0; i < blk.mCount; ++i) {
                if (spec->mCard != nullptr) {
                    blk.mCards[i] = nh5w_card_index(spec->mCard);
                } else {
                    const std::int64_t pid = blk.mPids[i];
                    blk.mCards[i] = nh5w_card_index(pid <= 0                ? "CONROD"
                                                    : ptube.count(pid) != 0 ? "CTUBE"
                                                                            : "CROD");
                }
            }
            mBlocks.push_back(std::move(blk));
        }
        for (const Nh5wBlock& blk : mBlocks)
            for (const std::int64_t e : blk.mEids)
                if (e <= 0)
                    nh5w_fail("element ids must be positive");
        {
            // A card's table is keyed by element id; two cards may share one.
            std::vector<std::pair<std::size_t, std::int64_t>> keyed;
            for (const Nh5wBlock& blk : mBlocks)
                for (std::size_t r = 0; r < blk.mCount; ++r)
                    keyed.emplace_back(blk.mCards[r], blk.mEids[r]);
            std::sort(keyed.begin(), keyed.end());
            if (std::adjacent_find(keyed.begin(), keyed.end()) != keyed.end())
                nh5w_fail("element ids are not unique");
        }
        mRows.resize(kNh5wNumCards);
        for (std::size_t bi = 0; bi < mBlocks.size(); ++bi)
            for (std::size_t r = 0; r < mBlocks[bi].mCount; ++r)
                mRows[mBlocks[bi].mCards[r]].push_back({mBlocks[bi].mEids[r], bi, r});
        for (auto& rows : mRows)
            std::stable_sort(rows.begin(), rows.end(), [](const Nh5wRow& rA, const Nh5wRow& rB) {
                return rA.mEid < rB.mEid;
            });
        for (const Nh5wBlock& blk : mBlocks)
            for (std::size_t r = 0; r < blk.mCount; ++r) {
                const Nh5wProp prop = kNh5wCards[blk.mCards[r]].mProp;
                if (prop != kNh5wNoProp)
                    mProps[prop].insert(blk.mPids[r]);
            }
        static const char* const names[] = {"PROD", "PSHELL", "PSOLID", "PTUBE"};
        std::map<std::int64_t, int> used;
        for (const int prop : {kNh5wPshell, kNh5wPsolid, kNh5wProd, kNh5wPtube})
            for (const std::int64_t p : mProps[prop]) {
                const auto it = used.find(p);
                if (it != used.end())
                    nh5w_fail("property id " + std::to_string(p) + " is used by " +
                              names[it->second] + " and " + names[prop] +
                              " elements: shells, solids and rods need ids of their own");
                used[p] = prop;
            }
        for (const auto& [pid, prop] : used)
            if (pid <= 0)
                nh5w_fail("property ids must be positive");
    }

    const Nh5wSignature& Signature() const { return mSignature; }
    std::size_t NumPoints() const { return mNpoints; }
    const std::vector<double>& Points() const { return mPoints; }
    const std::vector<Nh5wBlock>& Blocks() const { return mBlocks; }
    const std::vector<Nh5wRow>& Rows(std::size_t Card) const { return mRows[Card]; }
    const std::set<std::int64_t>& Props(int Prop) const { return mProps[Prop]; }
    const std::set<std::string>& Skipped() const { return mSkipped; }
    bool HasProperties() const {
        for (const auto& s : mProps)
            if (!s.empty())
                return true;
        return false;
    }

private:
    Nh5wSignature mSignature;
    std::size_t mNpoints = 0;
    std::vector<double> mPoints;
    std::vector<Nh5wBlock> mBlocks;
    std::vector<std::vector<Nh5wRow>> mRows;
    std::set<std::int64_t> mProps[4];
    std::set<std::string> mSkipped;
};

// -- steps ---------------------------------------------------------------------------

using Nh5wMembers = std::vector<std::pair<std::string, std::vector<double>>>;

struct Nh5wStep {
    std::int64_t mDomain = 1, mSubcase = 1, mStep = 0, mAnalysis = 1, mMode = 1;
    double mTime = 0.0, mEigi = 0.0;
    std::map<std::string, Nh5wMembers> mNodal;                           // table -> members
    std::map<std::pair<std::string, std::size_t>, Nh5wMembers> mElemental;  // (group, card)
    std::string mSignature;  // the tables and their members, in a canonical text
};

std::int64_t nh5w_scalar_int(const Mesh& rMesh, const std::string& rName, std::int64_t Default) {
    if (!rMesh.HasFieldData(rName))
        return Default;
    const NDArray& a = rMesh.FieldData(rName);
    if (a.Size() == 0)
        return Default;
    return detail::Int64View(a).Data()[0];
}

double nh5w_scalar_double(const Mesh& rMesh, const std::string& rName, double Default) {
    if (!rMesh.HasFieldData(rName))
        return Default;
    const NDArray& a = rMesh.FieldData(rName);
    if (a.Size() == 0)
        return Default;
    return detail::DoubleView(a).Data()[0];
}

double nh5w_time_of(const Mesh& rMesh) { return nh5w_scalar_double(rMesh, kSequenceTimeKey, 0.0); }

bool nh5w_starts_with(const std::string& rS, const char* pPrefix) {
    return rS.rfind(pPrefix, 0) == 0;
}

bool nh5w_ends_with(const std::string& rS, const std::string& rSuffix) {
    return rS.size() >= rSuffix.size() &&
           rS.compare(rS.size() - rSuffix.size(), rSuffix.size(), rSuffix) == 0;
}

/// Insert or replace a member in a table's ordered members.
void nh5w_put(Nh5wMembers& rMembers, const std::string& rName, std::vector<double> Values) {
    for (auto& [name, values] : rMembers)
        if (name == rName) {
            values = std::move(Values);
            return;
        }
    rMembers.emplace_back(rName, std::move(Values));
}

bool nh5w_quiet_point(const std::string& rName) {
    return rName == "nastran:cp" || rName == "nastran:cd" || nh5w_starts_with(rName, "nastran:");
}

/// `{table: {member: (npoints,)}}` from the point data, members in order.
void nh5w_nodal_members(const Mesh& rMesh, std::size_t NumPoints, std::set<std::string>& rTaken,
                        std::map<std::string, Nh5wMembers>& rTables) {
    const std::vector<std::string> names = rMesh.PointDataNames();  // sorted
    auto values_of = [&](const std::string& rName) {
        const NDArray& a = rMesh.PointData(rName);
        const detail::DoubleView view(a);
        return std::vector<double>(view.Data(), view.Data() + a.Size());
    };
    // The three columns of a (n, 3) array.
    auto vector_of = [&](const std::string& rName) {
        const NDArray& a = rMesh.PointData(rName);
        if (a.Size() != NumPoints * 3)
            nh5w_fail("'" + rName + "' must have 3 components per point");
        rTaken.insert(rName);
        const std::vector<double> flat = values_of(rName);
        std::array<std::vector<double>, 3> cols;
        for (std::size_t k = 0; k < 3; ++k) {
            cols[k].resize(NumPoints);
            for (std::size_t i = 0; i < NumPoints; ++i)
                cols[k][i] = flat[i * 3 + k];
        }
        return cols;
    };
    for (const std::string& name : names) {
        if (nh5w_quiet_point(name))
            continue;
        if (nh5w_ends_with(name, "_real") &&
            rMesh.HasPointData(name.substr(0, name.size() - 5) + "_imag")) {
            const std::string base = name.substr(0, name.size() - 5);
            const bool rot = nh5w_ends_with(base, "_ROT");
            const std::string table = (rot ? base.substr(0, base.size() - 4) : base) + "_CPLX";
            const std::string prefix = rot ? "R" : "";
            const auto re = vector_of(base + "_real");
            const auto im = vector_of(base + "_imag");
            const char axes[] = {'X', 'Y', 'Z'};
            for (std::size_t k = 0; k < 3; ++k) {
                nh5w_put(rTables[table], prefix + axes[k] + "R", re[k]);
                nh5w_put(rTables[table], prefix + axes[k] + "I", im[k]);
            }
        }
    }
    for (const std::string& name : names) {
        if (rTaken.count(name) || nh5w_quiet_point(name))
            continue;
        const NDArray& a = rMesh.PointData(name);
        const std::size_t colon = name.find(':');
        if (colon != std::string::npos) {
            if (a.Size() != NumPoints)
                continue;
            nh5w_put(rTables[name.substr(0, colon)], name.substr(colon + 1), values_of(name));
            rTaken.insert(name);
        } else if (nh5w_ends_with(name, "_ROT") && a.Size() == NumPoints * 3) {
            const std::string table = name.substr(0, name.size() - 4);
            const auto cols = vector_of(name);
            const char* axes[] = {"RX", "RY", "RZ"};
            for (std::size_t k = 0; k < 3; ++k)
                nh5w_put(rTables[table], axes[k], cols[k]);
        } else if (a.Size() == NumPoints * 3) {
            const auto cols = vector_of(name);
            const char* axes[] = {"X", "Y", "Z"};
            for (std::size_t k = 0; k < 3; ++k)
                nh5w_put(rTables[name], axes[k], cols[k]);
        } else if (a.Size() == NumPoints) {
            nh5w_put(rTables[name], "VALUE", values_of(name));
            rTaken.insert(name);
        }
    }
}

bool nh5w_is_group_name(const std::string& rS) {
    if (rS.empty() || !(rS[0] >= 'A' && rS[0] <= 'Z'))
        return false;
    for (const char c : rS)
        if (!((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_'))
            return false;
    return true;
}

/// A cell-data array of the cells of a card, in the card's row order; false when
/// no block of the card has it.
bool nh5w_gather(const Mesh& rMesh, const std::string& rName, const Nh5wModel& rModel,
                 std::size_t Card, std::vector<double>& rOut) {
    if (!rMesh.HasCellData(rName))
        return false;
    const std::size_t nblocks = rMesh.CellDataNumBlocks(rName);
    std::vector<std::vector<double>> firsts(rModel.Blocks().size());
    const auto& rows = rModel.Rows(Card);
    rOut.assign(rows.size(), std::numeric_limits<double>::quiet_NaN());
    std::map<std::size_t, std::vector<double>> by_block;
    for (const Nh5wRow& row : rows) {
        const Nh5wBlock& blk = rModel.Blocks()[row.mBlock];
        if (by_block.find(row.mBlock) == by_block.end()) {
            if (blk.mIndex >= nblocks)
                nh5w_fail("'" + rName + "' has no array for every cell block");
            const NDArray& a = rMesh.CellData(rName, blk.mIndex);
            const std::size_t n = blk.mCount;
            if (n && (a.Size() < n || a.Size() % n))
                nh5w_fail("'" + rName + "' must have one value per cell");
            const detail::DoubleView view(a);
            const std::size_t width = n ? a.Size() / n : 1;
            std::vector<double> first(n);
            for (std::size_t i = 0; i < n; ++i)
                first[i] = view.Data()[i * width];
            by_block.emplace(row.mBlock, std::move(first));
        }
    }
    for (std::size_t i = 0; i < rows.size(); ++i)
        rOut[i] = by_block[rows[i].mBlock][rows[i].mRow];
    return true;
}

void nh5w_element_members(const Mesh& rMesh, const Nh5wModel& rModel, std::set<std::string>& rTaken,
                          std::map<std::pair<std::string, std::size_t>, Nh5wMembers>& rOut) {
    std::map<std::string, std::vector<std::pair<std::string, std::string>>> groups;
    for (const std::string& name : rMesh.CellDataNames()) {
        if (name == "nastran:eid" || name == "nastran:pid")
            continue;
        const std::size_t colon = name.find(':');
        if (colon == std::string::npos)
            continue;
        const std::string group = name.substr(0, colon);
        const std::string member = name.substr(colon + 1);
        if (!nh5w_is_group_name(group) || member.find('@') != std::string::npos)
            continue;
        groups[group].emplace_back(member, name);
    }
    for (const auto& [group, members] : groups) {
        for (std::size_t card = 0; card < kNh5wNumCards; ++card) {
            if (kNh5wCards[card].mFamily == nullptr || rModel.Rows(card).empty())
                continue;
            Nh5wMembers columns;
            for (const auto& [member, name] : members) {
                std::vector<double> column;
                if (nh5w_gather(rMesh, name, rModel, card, column)) {
                    rTaken.insert(name);
                    columns.emplace_back(member, std::move(column));
                }
            }
            if (!columns.empty())
                rOut[{group, card}] = std::move(columns);
        }
    }
}

Nh5wStep nh5w_make_step(const Mesh& rMesh, const Nh5wModel& rModel, double Time, std::size_t Index,
                        std::set<std::string>& rTaken) {
    Nh5wStep step;
    step.mDomain = static_cast<std::int64_t>(Index) + 1;
    step.mSubcase = nh5w_scalar_int(rMesh, "nastran:subcase", 1);
    step.mStep = nh5w_scalar_int(rMesh, "nastran:step", static_cast<std::int64_t>(Index));
    step.mAnalysis = nh5w_scalar_int(rMesh, "nastran:analysis", Time == 0.0 ? 1 : 6);
    step.mTime = Time;
    step.mEigi = nh5w_scalar_double(rMesh, "nastran:eigi", 0.0);
    step.mMode = nh5w_scalar_int(rMesh, "nastran:mode", 1);
    nh5w_nodal_members(rMesh, rModel.NumPoints(), rTaken, step.mNodal);
    nh5w_element_members(rMesh, rModel, rTaken, step.mElemental);
    for (const auto& [table, members] : step.mNodal) {
        step.mSignature += "N:" + table;
        for (const auto& [m, v] : members)
            step.mSignature += "," + m;
        step.mSignature += ";";
    }
    for (const auto& [key, members] : step.mElemental) {
        step.mSignature += "E:" + key.first + "/" + std::to_string(key.second);
        for (const auto& [m, v] : members)
            step.mSignature += "," + m;
        step.mSignature += ";";
    }
    return step;
}

void nh5w_warn_dropped(const Mesh& rMesh, const Nh5wModel& rModel,
                       const std::set<std::string>& rTaken) {
    std::set<std::string> names;
    for (const std::string& n : rMesh.PointDataNames())
        if (!rTaken.count(n) && n != "nastran:cp" && n != "nastran:cd")
            names.insert(n);
    for (const std::string& n : rMesh.CellDataNames())
        if (!rTaken.count(n) && n != "nastran:eid" && n != "nastran:pid")
            names.insert(n);
    for (const std::string& n : rMesh.FieldDataNames())
        if (n != kSequenceTimeKey && !nh5w_starts_with(n, "nastran:") && !rTaken.count(n))
            names.insert(n);
    if (!names.empty()) {
        std::string joined;
        for (const std::string& n : names)
            joined += (joined.empty() ? "" : ", ") + n;
        log::warn("Nastran HDF5: arrays without a place in the file are not written: {}", joined);
    }
    if (!rModel.Skipped().empty()) {
        std::string joined;
        for (const std::string& n : rModel.Skipped())
            joined += (joined.empty() ? "" : ", ") + n;
        log::warn("Nastran HDF5: cell blocks with no element card are not written: {}", joined);
    }
}

// -- HDF5 tables ---------------------------------------------------------------------

struct Nh5wColumn {
    std::string mName;
    bool mFloat = false;
    bool mString = false;  // a fixed 4-byte string per row (mS holds the bytes)
    std::size_t mArity = 1;
    std::vector<std::int64_t> mI;
    std::vector<double> mF;
    std::string mS;
    std::size_t Rows() const {
        return (mString ? mS.size() / 4 : mFloat ? mF.size() : mI.size()) / mArity;
    }
    std::size_t Bytes() const { return mString ? 4 * mArity : 8 * mArity; }
};

Nh5wColumn nh5w_ints(std::string Name, std::vector<std::int64_t> Values, std::size_t Arity = 1) {
    Nh5wColumn c;
    c.mName = std::move(Name);
    c.mArity = Arity;
    c.mI = std::move(Values);
    return c;
}

Nh5wColumn nh5w_string4(std::string Name, const std::string& rValue, std::size_t Rows) {
    Nh5wColumn c;
    c.mName = std::move(Name);
    c.mString = true;
    std::string cell = rValue.substr(0, 4);
    cell.resize(4, '\0');
    for (std::size_t r = 0; r < Rows; ++r)
        c.mS += cell;
    return c;
}

Nh5wColumn nh5w_floats(std::string Name, std::vector<double> Values, std::size_t Arity = 1) {
    Nh5wColumn c;
    c.mName = std::move(Name);
    c.mFloat = true;
    c.mArity = Arity;
    c.mF = std::move(Values);
    return c;
}

/// Create the group at an absolute path, with the groups above it, unless it exists.
void nh5w_require_group(hid_t File, const std::string& rPath) {
    if (nh5_is_group(File, rPath))
        return;
    h5::Hid lcpl(H5Pcreate(H5P_LINK_CREATE), H5Pclose);
    H5Pset_create_intermediate_group(lcpl, 1);
    h5::Hid g(H5Gcreate2(File, rPath.c_str(), lcpl, H5P_DEFAULT, H5P_DEFAULT), H5Gclose);
    if (!g.Valid())
        nh5w_fail("could not create the group " + rPath);
}

/// A compound dataset at an absolute path, from its columns (equal row counts).
void nh5w_write_table(hid_t File, const std::string& rPath, const std::vector<Nh5wColumn>& rColumns) {
    const std::size_t rows = rColumns.empty() ? 0 : rColumns.front().Rows();
    std::size_t record = 0;
    std::vector<std::size_t> offsets;
    for (const Nh5wColumn& c : rColumns) {
        if (c.Rows() != rows)
            nh5w_fail("internal error: column " + c.mName + " of " + rPath + " has " +
                      std::to_string(c.Rows()) + " rows, not " + std::to_string(rows));
        offsets.push_back(record);
        record += c.Bytes();
    }
    h5::Hid file_type(H5Tcreate(H5T_COMPOUND, record), H5Tclose);
    h5::Hid mem_type(H5Tcreate(H5T_COMPOUND, record), H5Tclose);
    std::vector<h5::Hid> keep;
    for (std::size_t k = 0; k < rColumns.size(); ++k) {
        const Nh5wColumn& c = rColumns[k];
        hid_t file_member = c.mFloat ? H5T_IEEE_F64LE : H5T_STD_I64LE;
        hid_t mem_member = c.mFloat ? H5T_NATIVE_DOUBLE : H5T_NATIVE_INT64;
        if (c.mString) {
            keep.emplace_back(H5Tcopy(H5T_C_S1), H5Tclose);
            H5Tset_size(keep.back(), 4);
            H5Tset_strpad(keep.back(), H5T_STR_NULLPAD);
            file_member = mem_member = keep.back();
        }
        if (c.mArity > 1) {
            const hsize_t dims[1] = {static_cast<hsize_t>(c.mArity)};
            keep.emplace_back(H5Tarray_create2(file_member, 1, dims), H5Tclose);
            file_member = keep.back();
            keep.emplace_back(H5Tarray_create2(mem_member, 1, dims), H5Tclose);
            mem_member = keep.back();
        }
        H5Tinsert(file_type, c.mName.c_str(), offsets[k], file_member);
        H5Tinsert(mem_type, c.mName.c_str(), offsets[k], mem_member);
    }
    std::vector<char> buffer(record * rows);
    for (std::size_t k = 0; k < rColumns.size(); ++k) {
        const Nh5wColumn& c = rColumns[k];
        for (std::size_t r = 0; r < rows; ++r)
            std::memcpy(buffer.data() + r * record + offsets[k],
                        c.mString  ? static_cast<const void*>(c.mS.data() + r * 4 * c.mArity)
                        : c.mFloat ? static_cast<const void*>(c.mF.data() + r * c.mArity)
                                   : static_cast<const void*>(c.mI.data() + r * c.mArity),
                        c.Bytes());
    }
    const hsize_t dim = static_cast<hsize_t>(rows);
    h5::Hid space(H5Screate_simple(1, &dim, nullptr), H5Sclose);
    h5::Hid dset(H5Dcreate2(File, rPath.c_str(), file_type, space, H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT),
                 H5Dclose);
    if (!dset.Valid())
        nh5w_fail("could not create the table " + rPath);
    if (rows > 0 && H5Dwrite(dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0)
        nh5w_fail("could not write the table " + rPath);
}

std::vector<std::int64_t> nh5w_iota(std::size_t N, std::int64_t First = 1) {
    std::vector<std::int64_t> v(N);
    for (std::size_t i = 0; i < N; ++i)
        v[i] = First + static_cast<std::int64_t>(i);
    return v;
}

// Columns the MSC tables carry beyond the ones the reader uses, and the value this writer
// gives them: kind 'i' (int64), 'f' (float64) or 's' (a 4-byte string, `mText`), `mArity`
// values per row.
struct Nh5wSpecCol {
    const char* mName;
    char mKind;
    std::size_t mArity;
    double mValue;
    const char* mText;
};

std::vector<Nh5wColumn> nh5w_filled(std::size_t Rows, const std::vector<Nh5wSpecCol>& rSpec) {
    std::vector<Nh5wColumn> out;
    for (const Nh5wSpecCol& c : rSpec) {
        if (c.mKind == 's')
            out.push_back(nh5w_string4(c.mName, c.mText, Rows));
        else if (c.mKind == 'f')
            out.push_back(nh5w_floats(c.mName, std::vector<double>(Rows * c.mArity, c.mValue),
                                      c.mArity));
        else
            out.push_back(nh5w_ints(
                c.mName, std::vector<std::int64_t>(Rows * c.mArity, static_cast<std::int64_t>(c.mValue)),
                c.mArity));
    }
    return out;
}

const std::vector<Nh5wSpecCol>& nh5w_element_extra(std::string_view Card) {
    static const std::vector<Nh5wSpecCol> none;
    static const std::vector<Nh5wSpecCol> quad4 = {
        {"THETA", 'f', 1, 0.0, ""}, {"ZOFFS", 'f', 1, 0.0, ""}, {"TFLAG", 'i', 1, 0, ""},
        {"T", 'f', 4, 0.0, ""},     {"MCID", 'i', 1, -1, ""}};
    static const std::vector<Nh5wSpecCol> quad8 = {
        {"T", 'f', 4, 0.0, ""}, {"THETA", 'f', 1, 0.0, ""}, {"ZOFFS", 'f', 1, 0.0, ""},
        {"TFLAG", 'i', 1, 0, ""}, {"MCID", 'i', 1, -1, ""}};
    static const std::vector<Nh5wSpecCol> tria3 = {
        {"THETA", 'f', 1, 0.0, ""}, {"ZOFFS", 'f', 1, 0.0, ""}, {"TFLAG", 'i', 1, 0, ""},
        {"T", 'f', 3, 0.0, ""},     {"MCID", 'i', 1, -1, ""}};
    static const std::vector<Nh5wSpecCol> tria6 = {
        {"THETA", 'f', 1, 0.0, ""}, {"ZOFFS", 'f', 1, 0.0, ""}, {"T", 'f', 3, 0.0, ""},
        {"TFLAG", 'i', 1, 0, ""},   {"MCID", 'i', 1, -1, ""}};
    if (Card == "CQUAD" || Card == "CQUAD4")
        return quad4;
    if (Card == "CQUAD8")
        return quad8;
    if (Card == "CTRIA3")
        return tria3;
    if (Card == "CTRIA6")
        return tria6;
    return none;
}

const std::vector<Nh5wSpecCol>& nh5w_property_columns(int Prop) {
    static const std::vector<Nh5wSpecCol> prod = {{"MID", 'i', 1, 1, ""}, {"A", 'f', 1, 1.0, ""},
                                                  {"J", 'f', 1, 0.0, ""}, {"C", 'f', 1, 0.0, ""},
                                                  {"NSM", 'f', 1, 0.0, ""}};
    static const std::vector<Nh5wSpecCol> pshell = {
        {"MID1", 'i', 1, 1, ""},  {"T", 'f', 1, 1.0, ""},   {"MID2", 'i', 1, -1, ""},
        {"BK", 'f', 1, 1.0, ""},  {"MID3", 'i', 1, -1, ""}, {"TS", 'f', 1, 5.0 / 6.0, ""},
        {"NSM", 'f', 1, 0.0, ""}, {"Z1", 'f', 1, -0.5, ""}, {"Z2", 'f', 1, 0.5, ""},
        {"MID4", 'i', 1, -1, ""}};
    static const std::vector<Nh5wSpecCol> psolid = {
        {"MID", 'i', 1, 1, ""},   {"CORDM", 'i', 1, 0, ""}, {"IN", 'i', 1, 0, ""},
        {"STRESS", 'i', 1, 0, ""}, {"ISOP", 'i', 1, 0, ""},  {"FCTN", 's', 1, 0, "SMEC"}};
    static const std::vector<Nh5wSpecCol> ptube = {{"MID", 'i', 1, 1, ""}, {"OD", 'f', 1, 1.0, ""},
                                                   {"T", 'f', 1, 0.0, ""}, {"NSM", 'f', 1, 0.0, ""}};
    switch (Prop) {
        case kNh5wProd:
            return prod;
        case kNh5wPshell:
            return pshell;
        case kNh5wPsolid:
            return psolid;
        default:
            return ptube;
    }
}

void nh5w_write_model(hid_t File, const Nh5wModel& rModel) {
    const std::size_t n = rModel.NumPoints();
    auto domain = [](std::size_t Rows) {
        return nh5w_ints("DOMAIN_ID", std::vector<std::int64_t>(Rows, 1));
    };
    nh5w_require_group(File, "/NASTRAN/INPUT");
    nh5w_write_table(File, "/NASTRAN/INPUT/DOMAINS",
                     {nh5w_ints("ID", {1}), nh5w_ints("SE", {0}), nh5w_ints("AFPM", {0}),
                      nh5w_ints("TRMC", {0}), nh5w_ints("MODULE", {1})});
    nh5w_require_group(File, "/NASTRAN/INPUT/NODE");
    nh5w_write_table(File, "/NASTRAN/INPUT/NODE/GRID",
                     {nh5w_ints("ID", nh5w_iota(n)), nh5w_ints("CP", std::vector<std::int64_t>(n, 0)),
                      nh5w_floats("X", rModel.Points(), 3),
                      nh5w_ints("CD", std::vector<std::int64_t>(n, 0)),
                      nh5w_ints("PS", std::vector<std::int64_t>(n, 0)),
                      nh5w_ints("SEID", std::vector<std::int64_t>(n, 0)), domain(n)});
    nh5w_require_group(File, "/NASTRAN/INPUT/ELEMENT");
    for (std::size_t c = 0; c < kNh5wNumCards; ++c) {
        const auto& rows = rModel.Rows(c);
        if (rows.empty())
            continue;
        const Nh5wCard& card = kNh5wCards[c];
        std::vector<std::int64_t> eid, pid, nodes(rows.size() * card.mWidth, 0);
        for (std::size_t r = 0; r < rows.size(); ++r) {
            const Nh5wBlock& blk = rModel.Blocks()[rows[r].mBlock];
            eid.push_back(rows[r].mEid);
            pid.push_back(blk.mPids[rows[r].mRow]);
            const std::int64_t* conn = blk.mConn.data() + rows[r].mRow * blk.mNodes;
            for (std::size_t k = 0; k < blk.mNodes; ++k) {
                const std::size_t at =
                    (blk.mNodes > card.mLinear && card.mPerm != nullptr) ? card.mPerm[k] : k;
                nodes[r * card.mWidth + at] = conn[k] + 1;
            }
        }
        std::vector<Nh5wColumn> columns;
        const std::string_view name = card.mName;
        if (name == "CONROD") {
            std::vector<std::int64_t> g1, g2;
            for (std::size_t r = 0; r < rows.size(); ++r) {
                g1.push_back(nodes[r * 2]);
                g2.push_back(nodes[r * 2 + 1]);
            }
            columns.push_back(nh5w_ints("EID", eid));
            columns.push_back(nh5w_ints("G1", g1));
            columns.push_back(nh5w_ints("G2", g2));
            for (auto& col : nh5w_filled(rows.size(), {{"MID", 'i', 1, 1, ""},
                                                       {"A", 'f', 1, 1.0, ""},
                                                       {"J", 'f', 1, 0.0, ""},
                                                       {"C", 'f', 1, 0.0, ""},
                                                       {"NSM", 'f', 1, 0.0, ""}}))
                columns.push_back(std::move(col));
        } else if (name == "CONM2") {
            columns.push_back(nh5w_ints("EID", eid));
            columns.push_back(nh5w_ints("G", nodes));
            for (auto& col : nh5w_filled(rows.size(), {{"CID", 'i', 1, 0, ""},
                                                       {"M", 'f', 1, 1.0, ""},
                                                       {"X1", 'f', 1, 0.0, ""},
                                                       {"X2", 'f', 1, 0.0, ""},
                                                       {"X3", 'f', 1, 0.0, ""},
                                                       {"I1", 'f', 1, 0.0, ""},
                                                       {"I2", 'f', 2, 0.0, ""},
                                                       {"I3", 'f', 3, 0.0, ""}}))
                columns.push_back(std::move(col));
        } else {
            columns.push_back(nh5w_ints("EID", eid));
            columns.push_back(nh5w_ints("PID", pid));
            columns.push_back(nh5w_ints("G", nodes, card.mWidth));
            for (auto& col : nh5w_filled(rows.size(), nh5w_element_extra(name)))
                columns.push_back(std::move(col));
        }
        columns.push_back(domain(rows.size()));
        nh5w_write_table(File, std::string("/NASTRAN/INPUT/ELEMENT/") + card.mName, columns);
    }
    static const char* const names[] = {"PROD", "PSHELL", "PSOLID", "PTUBE"};
    nh5w_require_group(File, "/NASTRAN/INPUT/PROPERTY");
    for (const int prop : {kNh5wProd, kNh5wPshell, kNh5wPsolid, kNh5wPtube}) {
        const auto& pids = rModel.Props(prop);
        if (pids.empty())
            continue;
        std::vector<Nh5wColumn> columns;
        columns.push_back(nh5w_ints("PID", std::vector<std::int64_t>(pids.begin(), pids.end())));
        for (auto& col : nh5w_filled(pids.size(), nh5w_property_columns(prop)))
            columns.push_back(std::move(col));
        columns.push_back(domain(pids.size()));
        nh5w_write_table(File, std::string("/NASTRAN/INPUT/PROPERTY/") + names[prop], columns);
    }
    if (rModel.HasProperties()) {
        nh5w_require_group(File, "/NASTRAN/INPUT/MATERIAL");
        std::vector<Nh5wColumn> columns = nh5w_filled(
            1, {{"MID", 'i', 1, 1, ""},   {"E", 'f', 1, 1.0, ""},   {"G", 'f', 1, 0.0, ""},
                {"NU", 'f', 1, 0.3, ""},  {"RHO", 'f', 1, 0.0, ""}, {"A", 'f', 1, 0.0, ""},
                {"TREF", 'f', 1, 0.0, ""}, {"GE", 'f', 1, 0.0, ""},  {"ST", 'f', 1, 0.0, ""},
                {"SC", 'f', 1, 0.0, ""},  {"SS", 'f', 1, 0.0, ""},  {"MCSID", 'i', 1, 0, ""}});
        columns.push_back(domain(1));
        nh5w_write_table(File, "/NASTRAN/INPUT/MATERIAL/MAT1", columns);
    }
}

/// One result table over the steps that have rows for it, and its INDEX twin.
/// `Columns(step)` gives the step's `(member, values)` pairs (ids first).
struct Nh5wTableRows {
    std::vector<std::int64_t> mDomain;
    std::vector<Nh5wColumn> mColumns;
    std::vector<std::int64_t> mIndexDomain, mIndexPosition, mIndexLength;
};

void nh5w_emit(hid_t File, const std::string& rParent, const std::string& rName,
               const std::vector<std::pair<std::int64_t, std::vector<Nh5wColumn>>>& rPerStep) {
    Nh5wTableRows out;
    std::size_t position = 0;
    for (const auto& [domain, columns] : rPerStep) {
        const std::size_t count = columns.front().Rows();
        if (count == 0)
            continue;
        if (out.mColumns.empty())
            for (const Nh5wColumn& c : columns) {
                Nh5wColumn empty;
                empty.mName = c.mName;
                empty.mFloat = c.mFloat;
                empty.mString = c.mString;
                empty.mArity = c.mArity;
                out.mColumns.push_back(std::move(empty));
            }
        for (std::size_t k = 0; k < columns.size(); ++k) {
            auto& dst = out.mColumns[k];
            dst.mI.insert(dst.mI.end(), columns[k].mI.begin(), columns[k].mI.end());
            dst.mF.insert(dst.mF.end(), columns[k].mF.begin(), columns[k].mF.end());
            dst.mS += columns[k].mS;
        }
        out.mDomain.insert(out.mDomain.end(), count, domain);
        out.mIndexDomain.push_back(domain);
        out.mIndexPosition.push_back(static_cast<std::int64_t>(position));
        out.mIndexLength.push_back(static_cast<std::int64_t>(count));
        position += count;
    }
    if (out.mColumns.empty())
        return;
    out.mColumns.push_back(nh5w_ints("DOMAIN_ID", out.mDomain));
    nh5w_require_group(File, "/NASTRAN/RESULT/" + rParent);
    nh5w_write_table(File, "/NASTRAN/RESULT/" + rParent + "/" + rName, out.mColumns);
    nh5w_require_group(File, "/INDEX/NASTRAN/RESULT/" + rParent);
    nh5w_write_table(File, "/INDEX/NASTRAN/RESULT/" + rParent + "/" + rName,
                     {nh5w_ints("DOMAIN_ID", out.mIndexDomain),
                      nh5w_ints("POSITION", out.mIndexPosition),
                      nh5w_ints("LENGTH", out.mIndexLength)});
}

void nh5w_write_results(hid_t File, const Nh5wModel& rModel, const std::vector<Nh5wStep>& rSteps) {
    bool any = false;
    for (const Nh5wStep& s : rSteps)
        any = any || !s.mNodal.empty() || !s.mElemental.empty();
    if (!any)
        return;
    nh5w_require_group(File, "/NASTRAN/RESULT");
    {
        std::vector<std::int64_t> id, subcase, step, analysis, mode;
        std::vector<double> time, eigi;
        for (const Nh5wStep& s : rSteps) {
            id.push_back(s.mDomain);
            subcase.push_back(s.mSubcase);
            step.push_back(s.mStep);
            analysis.push_back(s.mAnalysis);
            time.push_back(s.mTime);
            eigi.push_back(s.mEigi);
            mode.push_back(s.mMode);
        }
        const std::vector<std::int64_t> zeros(rSteps.size(), 0);
        std::vector<Nh5wColumn> columns = {
            nh5w_ints("ID", id),         nh5w_ints("SUBCASE", subcase),
            nh5w_ints("STEP", step),     nh5w_ints("ANALYSIS", analysis),
            nh5w_floats("TIME_FREQ_EIGR", time), nh5w_floats("EIGI", eigi),
            nh5w_ints("MODE", mode)};
        for (const char* extra : {"DESIGN_CYCLE", "RANDOM", "SE", "AFPM", "TRMC", "INSTANCE",
                                  "MODULE", "SUBSTEP", "IMPFID"})
            columns.push_back(nh5w_ints(extra, zeros));
        nh5w_write_table(File, "/NASTRAN/RESULT/DOMAINS", columns);
    }
    std::set<std::string> tables;
    for (const Nh5wStep& s : rSteps)
        for (const auto& [t, m] : s.mNodal)
            tables.insert(t);
    for (const std::string& table : tables) {
        std::vector<std::pair<std::int64_t, std::vector<Nh5wColumn>>> per_step;
        for (const Nh5wStep& s : rSteps) {
            const auto it = s.mNodal.find(table);
            if (it == s.mNodal.end())
                continue;
            const Nh5wMembers& members = it->second;
            std::vector<std::size_t> keep;
            for (std::size_t i = 0; i < rModel.NumPoints(); ++i) {
                bool all_nan = true;
                for (const auto& [m, v] : members)
                    all_nan = all_nan && std::isnan(v[i]);
                if (!all_nan)
                    keep.push_back(i);
            }
            std::vector<Nh5wColumn> columns;
            std::vector<std::int64_t> ids;
            for (const std::size_t i : keep)
                ids.push_back(static_cast<std::int64_t>(i) + 1);
            columns.push_back(nh5w_ints("ID", ids));
            for (const auto& [m, v] : members) {
                std::vector<double> picked;
                for (const std::size_t i : keep)
                    picked.push_back(v[i]);
                columns.push_back(nh5w_floats(m, std::move(picked)));
            }
            per_step.emplace_back(s.mDomain, std::move(columns));
        }
        nh5w_emit(File, "NODAL", table, per_step);
    }
    std::set<std::pair<std::string, std::size_t>> keys;
    for (const Nh5wStep& s : rSteps)
        for (const auto& [k, m] : s.mElemental)
            keys.insert(k);
    for (const auto& key : keys) {
        const std::size_t card = key.second;
        std::vector<std::pair<std::int64_t, std::vector<Nh5wColumn>>> per_step;
        const auto& rows = rModel.Rows(card);
        for (const Nh5wStep& s : rSteps) {
            const auto it = s.mElemental.find(key);
            if (it == s.mElemental.end())
                continue;
            const Nh5wMembers& members = it->second;
            std::vector<std::size_t> keep;
            for (std::size_t i = 0; i < rows.size(); ++i) {
                bool all_nan = true;
                for (const auto& [m, v] : members)
                    all_nan = all_nan && std::isnan(v[i]);
                if (!all_nan)
                    keep.push_back(i);
            }
            std::vector<Nh5wColumn> columns;
            std::vector<std::int64_t> eids;
            for (const std::size_t i : keep)
                eids.push_back(rows[i].mEid);
            columns.push_back(nh5w_ints(key.first == "ENERGY" ? "ID" : "EID", eids));
            for (const auto& [m, v] : members) {
                std::vector<double> picked;
                for (const std::size_t i : keep)
                    picked.push_back(v[i]);
                columns.push_back(nh5w_floats(m, std::move(picked)));
            }
            per_step.emplace_back(s.mDomain, std::move(columns));
        }
        nh5w_emit(File, "ELEMENTAL/" + key.first, kNh5wCards[card].mFamily, per_step);
    }
}

void nh5w_assemble(const std::string& rPath, const Nh5wModel& rModel,
                   const std::vector<Nh5wStep>& rSteps) {
    detail::LibraryLock lock;
    h5::SilenceErrors quiet;
    h5::Hid file = h5::create_file(rPath);
    if (!file.Valid())
        nh5w_fail("could not create " + rPath);
    h5::write_attr_int_array(file, "SCHEMA", {20200});
    nh5w_write_model(file, rModel);
    nh5w_write_results(file, rModel, rSteps);
}

}  // namespace

void write_nastran_h5(const std::string& rPath, const Mesh& rMesh) {
    const detail::ProvenanceSlotlessWrite slotless;
    const Nh5wModel model(rMesh);
    std::set<std::string> taken;
    std::vector<Nh5wStep> steps;
    steps.push_back(nh5w_make_step(rMesh, model, nh5w_time_of(rMesh), 0, taken));
    nh5w_warn_dropped(rMesh, model, taken);
    nh5w_assemble(rPath, model, steps);
}

struct NastranH5SeriesWriter::Impl {
    std::string mPath;
    std::optional<Nh5wModel> mModel;
    std::vector<Nh5wStep> mSteps;
    bool mDone = false;
};

NastranH5SeriesWriter::NastranH5SeriesWriter(const std::string& rPath)
    : mpImpl(std::make_unique<Impl>()) {
    mpImpl->mPath = rPath;
}

NastranH5SeriesWriter::~NastranH5SeriesWriter() = default;
NastranH5SeriesWriter::NastranH5SeriesWriter(NastranH5SeriesWriter&&) noexcept = default;
NastranH5SeriesWriter& NastranH5SeriesWriter::operator=(NastranH5SeriesWriter&&) noexcept = default;

void NastranH5SeriesWriter::Write(double Time, const Mesh& rMesh) {
    if (!mpImpl)
        nh5w_fail("series: the writer was moved from");
    Impl& s = *mpImpl;
    const std::size_t index = s.mSteps.size();
    if (!s.mModel) {
        s.mModel.emplace(rMesh);
    } else if (s.mModel->Signature() != nh5w_signature(rMesh)) {
        nh5w_fail("series: step " + std::to_string(index) +
                  "'s points or cells differ from the first step's; a results file holds one "
                  "model -- write one file per step with '{step}'");
    }
    std::set<std::string> taken;
    Nh5wStep step = nh5w_make_step(rMesh, *s.mModel, Time, index, taken);
    if (s.mSteps.empty()) {
        nh5w_warn_dropped(rMesh, *s.mModel, taken);
    } else if (step.mSignature != s.mSteps.front().mSignature) {
        nh5w_fail("series: step " + std::to_string(index) +
                  "'s result arrays differ from the first step's; a result table holds the same "
                  "members in every domain");
    }
    s.mSteps.push_back(std::move(step));
}

std::size_t NastranH5SeriesWriter::NumSteps() const noexcept {
    return mpImpl ? mpImpl->mSteps.size() : 0;
}

void NastranH5SeriesWriter::Finalize() {
    if (!mpImpl || mpImpl->mDone)
        return;
    if (mpImpl->mSteps.empty())
        nh5w_fail("series: no step was written to " + mpImpl->mPath);
    mpImpl->mDone = true;
    nh5w_assemble(mpImpl->mPath, *mpImpl->mModel, mpImpl->mSteps);
}

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
