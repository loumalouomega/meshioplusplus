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
constexpr const char* kNh5Materials = "/NASTRAN/INPUT/MATERIAL";
constexpr const char* kNh5Domains = "/NASTRAN/RESULT/DOMAINS";
constexpr const char* kNh5Nodal = "/NASTRAN/RESULT/NODAL";
constexpr const char* kNh5Elemental = "/NASTRAN/RESULT/ELEMENTAL";

[[noreturn]] void nh5_fail(const std::string& rMessage) {
    throw ReadError("MSC Nastran HDF5: " + rMessage);
}

// The property and material members `nastran:<CARD>_<id>:<MEMBER>` field data sets and
// reads: the field data name, the table's column when it is called otherwise, whether it is an
// integer (the MID ones), and the value the writer gives a member the mesh has none for.
struct Nh5Member {
    const char* mName;
    const char* mColumn;  // null: the same
    bool mInt;
    double mDefault;
};
struct Nh5ValueCard {
    const char* mCard;
    const char* mKey;  // the id column
    std::vector<Nh5Member> mMembers;
};
const std::vector<Nh5ValueCard>& nh5_value_cards() {
    static const std::vector<Nh5ValueCard> cards = {
        {"PSHELL",
         "PID",
         {{"MID1", nullptr, true, 1}, {"T", nullptr, false, 1.0}, {"NSM", nullptr, false, 0.0}}},
        {"PROD",
         "PID",
         {{"MID", nullptr, true, 1},
          {"A", nullptr, false, 1.0},
          {"J", nullptr, false, 0.0},
          {"C", nullptr, false, 0.0},
          {"NSM", nullptr, false, 0.0}}},
        {"PTUBE",
         "PID",
         {{"MID", nullptr, true, 1},
          {"OD", nullptr, false, 1.0},
          {"T", nullptr, false, 0.0},
          {"NSM", nullptr, false, 0.0}}},
        {"PSOLID", "PID", {{"MID", nullptr, true, 1}}},
        {"PBAR",
         "PID",
         {{"MID", nullptr, true, 1},
          {"A", nullptr, false, 1.0},
          {"I1", nullptr, false, 0.0},
          {"I2", nullptr, false, 0.0},
          {"J", nullptr, false, 0.0},
          {"NSM", nullptr, false, 0.0},
          {"C1", nullptr, false, 0.0},
          {"C2", nullptr, false, 0.0},
          {"D1", nullptr, false, 0.0},
          {"D2", nullptr, false, 0.0},
          {"E1", nullptr, false, 0.0},
          {"E2", nullptr, false, 0.0},
          {"F1", nullptr, false, 0.0},
          {"F2", nullptr, false, 0.0},
          {"K1", nullptr, false, 0.0},
          {"K2", nullptr, false, 0.0},
          {"I12", nullptr, false, 0.0}}},
        {"PELAS",
         "PID",
         {{"K", nullptr, false, 1.0}, {"GE", nullptr, false, 0.0}, {"S", nullptr, false, 0.0}}},
        {"PDAMP", "PID", {{"B", nullptr, false, 1.0}}},
        {"PVISC", "PID", {{"CE", nullptr, false, 1.0}, {"CR", nullptr, false, 0.0}}},
        {"PSHEAR",
         "PID",
         {{"MID", nullptr, true, 1},
          {"T", nullptr, false, 1.0},
          {"NSM", nullptr, false, 0.0},
          {"F1", nullptr, false, 0.0},
          {"F2", nullptr, false, 0.0}}},
        {"MAT1",
         "MID",
         {{"E", nullptr, false, 1.0},
          {"G", nullptr, false, 0.0},
          {"NU", nullptr, false, 0.3},
          {"RHO", nullptr, false, 0.0},
          {"A", nullptr, false, 0.0},
          {"TREF", nullptr, false, 0.0},
          {"GE", nullptr, false, 0.0}}},
        {"MAT2",
         "MID",
         {{"G11", nullptr, false, 1.0},
          {"G12", nullptr, false, 0.0},
          {"G13", nullptr, false, 0.0},
          {"G22", nullptr, false, 1.0},
          {"G23", nullptr, false, 0.0},
          {"G33", nullptr, false, 1.0},
          {"RHO", nullptr, false, 0.0},
          {"A1", nullptr, false, 0.0},
          {"A2", nullptr, false, 0.0},
          {"A3", "A12", false, 0.0},
          {"TREF", nullptr, false, 0.0},
          {"GE", nullptr, false, 0.0},
          {"ST", nullptr, false, 0.0},
          {"SC", nullptr, false, 0.0},
          {"SS", nullptr, false, 0.0},
          {"MCSID", nullptr, true, 0}}},
        {"MAT4",
         "MID",
         {{"K", nullptr, false, 1.0},
          {"CP", nullptr, false, 0.0},
          {"RHO", nullptr, false, 1.0},
          {"H", nullptr, false, 0.0},
          {"MU", nullptr, false, 0.0},
          {"HGEN", nullptr, false, 1.0},
          {"REFENTH", nullptr, false, 0.0},
          {"TCH", nullptr, false, 0.0},
          {"TDELTA", nullptr, false, 0.0},
          {"QLAT", nullptr, false, 0.0}}},
    };
    return cards;
}

/// Whether @p rCard is a material card (its table is below /MATERIAL).
bool nh5_is_material(const std::string& rCard) {
    return rCard.compare(0, 3, "MAT") == 0;
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
    std::vector<std::vector<std::int64_t>> comps;  // per card: two components per element
    std::vector<std::vector<double>> bars;         // per card: 21 doubles per element, or none
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
        {
            // A spring's or damper's components; the values the cell data keep: a CBAR's and
            // CBEAM's orientation, G0, pins and offsets (the CBAR layout of the shared model),
            // a CELAS2's K, GE and S, a CDAMP2's B. Each element takes 21 columns.
            const auto ints = [&](const char* pName) {
                return nh5_has_member(members, pName) ? nh5_int_member(f, path, pName, 0, n)
                                                      : std::vector<std::int64_t>(n, 0);
            };
            const auto floats = [&](const char* pName, std::size_t Width) {
                std::size_t w = 1;
                if (!nh5_has_member(members, pName))
                    return std::vector<double>(n * Width, 0.0);
                std::vector<double> v = nh5_float_member_all(f, path, pName, 0, n, w);
                v.resize(n * Width, 0.0);
                return v;
            };
            std::vector<std::int64_t> c(2 * n, 0);
            std::vector<double> bar;
            if (card == "CELAS1" || card == "CDAMP1" || card == "CELAS2" || card == "CDAMP2") {
                const auto c1 = ints("C1");
                const auto c2 = ints("C2");
                for (std::size_t i = 0; i < n; ++i) {
                    c[2 * i] = c1[i];
                    c[2 * i + 1] = c2[i];
                }
            }
            const auto put = [&](std::size_t I, const std::vector<double>& rRow) {
                bar.resize((I + 1) * 21, 0.0);
                std::copy(rRow.begin(), rRow.end(), bar.begin() + static_cast<long>(I * 21));
            };
            if (card == "CELAS2") {
                const auto k = floats("K", 1), ge = floats("GE", 1), sv = floats("S", 1);
                for (std::size_t i = 0; i < n; ++i)
                    put(i, {k[i], ge[i], sv[i]});
            } else if (card == "CDAMP2") {
                const auto b = floats("B", 1);
                for (std::size_t i = 0; i < n; ++i)
                    put(i, {b[i]});
            } else if (card == "CBAR") {
                const auto flag = ints("FLAG"), go = ints("GO"), pa = ints("PA"), pb = ints("PB");
                const auto x1 = floats("X1", 1), x2 = floats("X2", 1), x3 = floats("X3", 1);
                const char* const w_names[6] = {"W1A", "W2A", "W3A", "W1B", "W2B", "W3B"};
                std::vector<std::vector<double>> w;
                for (const char* name : w_names)
                    w.push_back(floats(name, 1));
                for (std::size_t i = 0; i < n; ++i) {
                    const bool g0 = flag[i] == 2;
                    std::vector<double> row = {
                        g0 ? 0.0 : x1[i],           g0 ? 0.0 : x2[i],
                        g0 ? 0.0 : x3[i],           g0 ? static_cast<double>(go[i]) : 0.0,
                        static_cast<double>(pa[i]), static_cast<double>(pb[i])};
                    for (std::size_t k = 0; k < 6; ++k)
                        row.push_back(w[k][i]);
                    put(i, row);
                }
            } else if (card == "CBEAM") {
                const auto flag = ints("F"), g0v = ints("G0"), pa = ints("PA"), pb = ints("PB");
                const auto x = floats("X", 3), wa = floats("WA", 3), wb = floats("WB", 3);
                for (std::size_t i = 0; i < n; ++i) {
                    const std::int64_t low = flag[i] & 3;
                    std::vector<double> row;
                    for (std::size_t k = 0; k < 3; ++k)
                        row.push_back(low <= 1 ? x[3 * i + k] : 0.0);
                    row.push_back(low == 2 ? static_cast<double>(g0v[i]) : 0.0);
                    row.push_back(static_cast<double>(pa[i]));
                    row.push_back(static_cast<double>(pb[i]));
                    for (std::size_t k = 0; k < 3; ++k)
                        row.push_back(wa[3 * i + k]);
                    for (std::size_t k = 0; k < 3; ++k)
                        row.push_back(wb[3 * i + k]);
                    put(i, row);
                }
            }
            comps.push_back(std::move(c));
            bars.push_back(std::move(bar));
        }
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
    detail::NastranCellIds ids;
    const detail::NastranCells model =
        detail::nastran_add_cells(mesh, cards, file.mGridIndex, file.mScalarPoints, ptype,
                                  "MSC Nastran HDF5", comps, bars, &ids);
    if (!ids.mCard.empty()) {
        mesh.AddCellData("nastran:card", std::move(ids.mCard));
        mesh.AddCellData("nastran:c1", std::move(ids.mC1));
        mesh.AddCellData("nastran:c2", std::move(ids.mC2));
    }
    if (ids.mHasBar) {
        mesh.AddCellData("nastran:x", std::move(ids.mX));
        mesh.AddCellData("nastran:g0", std::move(ids.mG0));
        mesh.AddCellData("nastran:pin", std::move(ids.mPin));
        mesh.AddCellData("nastran:w", std::move(ids.mW));
    }
    if (ids.mHasSpring) {
        mesh.AddCellData("nastran:k", std::move(ids.mK));
        mesh.AddCellData("nastran:ge", std::move(ids.mGe));
        mesh.AddCellData("nastran:s", std::move(ids.mS));
    }
    if (ids.mHasDamper)
        mesh.AddCellData("nastran:b", std::move(ids.mB));
    const auto& cell_index = model.mCellIndex;
    const std::size_t ncells = model.mNumCells;
    for (const Nh5ValueCard& card : nh5_value_cards()) {  // property and material values
        const std::string path = (nh5_is_material(card.mCard) ? std::string(kNh5Materials)
                                                              : std::string(kNh5Properties)) +
                                 "/" + card.mCard;
        if (!nh5_is_dataset(f, path))
            continue;
        const auto members = h5::compound_members(f, path);
        if (!nh5_has_member(members, card.mKey))
            continue;
        const std::size_t rows = nh5_rows(f, path);
        const auto ids = nh5_int_member(f, path, card.mKey, 0, rows);
        for (const Nh5Member& m : card.mMembers) {
            const std::string member = m.mName;
            const std::string column = m.mColumn != nullptr ? m.mColumn : m.mName;
            const bool is_int = m.mInt;
            if (!nh5_has_member(members, column))
                continue;
            std::vector<double> floats;
            std::vector<std::int64_t> ints;
            if (is_int)
                ints = nh5_int_member(f, path, column, 0, rows);
            else
                floats = nh5_float_member(f, path, column, 0, rows);
            for (std::size_t r = 0; r < ids.size(); ++r) {
                const std::string name = std::string("nastran:") + card.mCard + "_" +
                                         std::to_string(ids[r]) + ":" + member;
                if (mesh.HasFieldData(name))
                    continue;
                NDArray a(is_int ? DType::Int64 : DType::Float64, {1});
                if (is_int)
                    a.As<std::int64_t>()[0] = ints[r];
                else
                    a.As<double>()[0] = floats[r];
                mesh.AddFieldData(name, std::move(a));
            }
        }
    }
    // PBEAM: MID, NSEG (NSEGS) and CCF, then SO, XXB, A, I1, I2, I12, J and NSM per station and
    // the outline points C1, C2, D1 and D2 as (n, 2) rows (the columns C1/C2, D1/D2, E1/E2,
    // F1/F2). The stations are NSEGS + 1: the first in slot 0, the last in slot 10, any
    // between them in the slots after the first.
    if (const std::string path = std::string(kNh5Properties) + "/PBEAM"; nh5_is_dataset(f, path)) {
        const auto members = h5::compound_members(f, path);
        const std::size_t rows = nh5_rows(f, path);
        if (nh5_has_member(members, "PID") && nh5_has_member(members, "MID") &&
            nh5_has_member(members, "NSEGS")) {
            const auto pid = nh5_int_member(f, path, "PID", 0, rows);
            const auto mid = nh5_int_member(f, path, "MID", 0, rows);
            const auto nsegs = nh5_int_member(f, path, "NSEGS", 0, rows);
            const auto ccf = nh5_has_member(members, "CCF")
                                 ? nh5_int_member(f, path, "CCF", 0, rows)
                                 : std::vector<std::int64_t>(rows, 0);
            const char* const stations[] = {"SO", "XXB", "A", "I1", "I2", "I12", "J", "NSM"};
            const char* const points[4][2] = {
                {"C1", "C2"}, {"D1", "D2"}, {"E1", "E2"}, {"F1", "F2"}};
            const char* const point_names[] = {"C1", "C2", "D1", "D2"};
            std::map<std::string, std::vector<double>> columns;  // column -> rows x 11
            for (const char* name : stations)
                if (nh5_has_member(members, name)) {
                    std::size_t w = 1;
                    columns[name] = nh5_float_member_all(f, path, name, 0, rows, w);
                    columns[name].resize(rows * 11, 0.0);
                }
            for (const auto& pair : points)
                for (const char* name : pair)
                    if (nh5_has_member(members, name)) {
                        std::size_t w = 1;
                        columns[name] = nh5_float_member_all(f, path, name, 0, rows, w);
                        columns[name].resize(rows * 11, 0.0);
                    }
            for (std::size_t r = 0; r < rows; ++r) {
                const std::string prefix = "nastran:PBEAM_" + std::to_string(pid[r]) + ":";
                if (nsegs[r] < 1 || nsegs[r] > 10 || mesh.HasFieldData(prefix + "MID"))
                    continue;
                const auto n = static_cast<std::size_t>(nsegs[r]) + 1;
                const auto slot = [&](std::size_t I) {
                    return I == 0 ? std::size_t{0} : I == n - 1 ? std::size_t{10} : I;
                };
                const auto put_int = [&](const char* pName, std::int64_t V) {
                    NDArray a(DType::Int64, {1});
                    a.As<std::int64_t>()[0] = V;
                    mesh.AddFieldData(prefix + pName, std::move(a));
                };
                put_int("MID", mid[r]);
                put_int("NSEG", nsegs[r]);
                put_int("CCF", ccf[r]);
                for (const char* name : stations) {
                    const auto it = columns.find(name);
                    if (it == columns.end())
                        continue;
                    NDArray a(DType::Float64, {n});
                    for (std::size_t i = 0; i < n; ++i)
                        a.As<double>()[i] = it->second[r * 11 + slot(i)];
                    mesh.AddFieldData(prefix + name, std::move(a));
                }
                for (std::size_t p = 0; p < 4; ++p) {
                    const auto u = columns.find(points[p][0]);
                    const auto v = columns.find(points[p][1]);
                    if (u == columns.end() || v == columns.end())
                        continue;
                    NDArray a(DType::Float64, {n, std::size_t{2}});
                    for (std::size_t i = 0; i < n; ++i) {
                        a.As<double>()[2 * i] = u->second[r * 11 + slot(i)];
                        a.As<double>()[2 * i + 1] = v->second[r * 11 + slot(i)];
                    }
                    mesh.AddFieldData(prefix + point_names[p], std::move(a));
                }
            }
        }
    }
    // PCOMP (the grouped IDENTITY and PLY tables): Z0, NSM, SB, TREF, GE, FT and SYM (0), and
    // PLY, the plies as (n, 4) rows of MID, T, THETA and SOUT.
    {
        const std::string ident = std::string(kNh5Properties) + "/PCOMP/IDENTITY";
        const std::string plies = std::string(kNh5Properties) + "/PCOMP/PLY";
        if (nh5_is_dataset(f, ident) && nh5_is_dataset(f, plies)) {
            const auto members = h5::compound_members(f, ident);
            const auto ply_members = h5::compound_members(f, plies);
            if (nh5_has_member(members, "PID") && nh5_has_member(members, "PLY_POS") &&
                nh5_has_member(members, "PLY_LEN") && nh5_has_member(ply_members, "MID") &&
                nh5_has_member(ply_members, "T") && nh5_has_member(ply_members, "THETA") &&
                nh5_has_member(ply_members, "SOUT")) {
                const std::size_t rows = nh5_rows(f, ident);
                const std::size_t total = nh5_rows(f, plies);
                const auto pid = nh5_int_member(f, ident, "PID", 0, rows);
                const auto pos = nh5_int_member(f, ident, "PLY_POS", 0, rows);
                const auto len = nh5_int_member(f, ident, "PLY_LEN", 0, rows);
                const auto mid = nh5_int_member(f, plies, "MID", 0, total);
                const auto t = nh5_float_member(f, plies, "T", 0, total);
                const auto theta = nh5_float_member(f, plies, "THETA", 0, total);
                const auto sout = nh5_int_member(f, plies, "SOUT", 0, total);
                const char* const scalars[] = {"Z0", "NSM", "SB", "TREF", "GE"};
                std::map<std::string, std::vector<double>> columns;
                for (const char* name : scalars)
                    if (nh5_has_member(members, name))
                        columns[name] = nh5_float_member(f, ident, name, 0, rows);
                const auto ft = nh5_has_member(members, "FT")
                                    ? nh5_int_member(f, ident, "FT", 0, rows)
                                    : std::vector<std::int64_t>(rows, 0);
                for (std::size_t r = 0; r < rows; ++r) {
                    const std::string prefix = "nastran:PCOMP_" + std::to_string(pid[r]) + ":";
                    if (mesh.HasFieldData(prefix + "PLY") || len[r] <= 0 || pos[r] < 0 ||
                        static_cast<std::size_t>(pos[r] + len[r]) > total)
                        continue;
                    const auto count = static_cast<std::size_t>(len[r]);
                    NDArray a(DType::Float64, {count, std::size_t{4}});
                    for (std::size_t j = 0; j < count; ++j) {
                        const auto at = static_cast<std::size_t>(pos[r]) + j;
                        a.As<double>()[4 * j] = static_cast<double>(mid[at]);
                        a.As<double>()[4 * j + 1] = t[at];
                        a.As<double>()[4 * j + 2] = theta[at];
                        a.As<double>()[4 * j + 3] = static_cast<double>(sout[at]);
                    }
                    mesh.AddFieldData(prefix + "PLY", std::move(a));
                    for (const char* name : scalars) {
                        const auto it = columns.find(name);
                        if (it == columns.end())
                            continue;
                        NDArray v(DType::Float64, {1});
                        v.As<double>()[0] = it->second[r];
                        mesh.AddFieldData(prefix + name, std::move(v));
                    }
                    NDArray f_t(DType::Int64, {1});
                    f_t.As<std::int64_t>()[0] = ft[r];
                    mesh.AddFieldData(prefix + "FT", std::move(f_t));
                    NDArray sym(DType::Int64, {1});
                    sym.As<std::int64_t>()[0] = 0;
                    mesh.AddFieldData(prefix + "SYM", std::move(sym));
                }
            }
        }
    }
    // GRID ids that are not 1 to n, as point data: only when point data is wanted
    // (`points_only` drops it, and `arrays` names what to keep).
    if (rOpts.WantsAnyData() && rOpts.WantsArray("nastran:nid"))
        for (std::size_t i = 0; i < file.mGridIds.size(); ++i)
            if (file.mGridIds[i] != static_cast<std::int64_t>(i + 1)) {  // not 1 to n: keep them
                NDArray nid(DType::Int64, {file.mGridIds.size()});
                std::copy(file.mGridIds.begin(), file.mGridIds.end(), nid.As<std::int64_t>());
                mesh.AddPointData("nastran:nid", std::move(nid));
                break;
            }

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

/// MSC's table column sets, frozen from the fixtures (`tools/gen_nastran_h5_schema.py`): the
/// path below `/NASTRAN` and the columns as `NAME:i8`, `f8`, `S4` or `S8`, with `[n]` for an
/// array of `n`.
struct Nh5wSchemaEntry {
    const char* mPath;
    const char* mColumns;
};
// BEGIN GENERATED SCHEMA (tools/gen_nastran_h5_schema.py)
constexpr Nh5wSchemaEntry kNh5wSchema[] = {
    {"INPUT/ELEMENT/CBAR",
     "EID:i8 PID:i8 GA:i8 GB:i8 FLAG:i8 X1:f8 X2:f8 X3:f8 GO:i8 PA:i8 PB:i8 W1A:f8 W2A:f8 W3A:f8 "
     "W1B:f8 W2B:f8 W3B:f8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CBEAM",
     "EID:i8 PID:i8 GA:i8 GB:i8 SA:i8 SB:i8 X:f8[3] G0:i8 F:i8 PA:i8 PB:i8 WA:f8[3] WB:f8[3] "
     "DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CDAMP1", "EID:i8 PID:i8 G1:i8 G2:i8 C1:i8 C2:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CDAMP2", "EID:i8 B:f8 G1:i8 G2:i8 C1:i8 C2:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CDAMP3", "EID:i8 PID:i8 S1:i8 S2:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CDAMP4", "EID:i8 B:f8 S1:i8 S2:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CELAS1", "EID:i8 PID:i8 G1:i8 G2:i8 C1:i8 C2:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CELAS2", "EID:i8 K:f8 G1:i8 G2:i8 C1:i8 C2:i8 GE:f8 S:f8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CELAS3", "EID:i8 PID:i8 S1:i8 S2:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CELAS4", "EID:i8 K:f8 S1:i8 S2:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CHBDYE",
     "EID:i8 EID2:i8 SIDE:i8 IVIEWF:i8 IVIEWB:i8 RADMIDF:i8 RADMIDB:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CHEXA", "EID:i8 PID:i8 G:i8[20] DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CONM2",
     "EID:i8 G:i8 CID:i8 M:f8 X1:f8 X2:f8 X3:f8 I1:f8 I2:f8[2] I3:f8[3] DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CONROD", "EID:i8 G1:i8 G2:i8 MID:i8 A:f8 J:f8 C:f8 NSM:f8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CPENTA", "EID:i8 PID:i8 G:i8[15] DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CQUAD4",
     "EID:i8 PID:i8 G:i8[4] THETA:f8 ZOFFS:f8 TFLAG:i8 T:f8[4] MCID:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CQUAD8",
     "EID:i8 PID:i8 G:i8[8] T:f8[4] THETA:f8 ZOFFS:f8 TFLAG:i8 MCID:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CQUADR",
     "EID:i8 PID:i8 G:i8[4] THETA:f8 ZOFFS:f8 TFLAG:i8 T:f8[4] MCID:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CROD", "EID:i8 PID:i8 G:i8[2] DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CSHEAR", "EID:i8 PID:i8 G:i8[4] DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CTETRA", "EID:i8 PID:i8 G:i8[10] DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CTRIA3",
     "EID:i8 PID:i8 G:i8[3] THETA:f8 ZOFFS:f8 TFLAG:i8 T:f8[3] MCID:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CTRIA6",
     "EID:i8 PID:i8 G:i8[6] THETA:f8 ZOFFS:f8 T:f8[3] TFLAG:i8 MCID:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CTRIAR",
     "EID:i8 PID:i8 G:i8[3] THETA:f8 ZOFFS:f8 TFLAG:i8 T:f8[3] MCID:i8 DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CTUBE", "EID:i8 PID:i8 G:i8[2] DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/CVISC", "EID:i8 PID:i8 G:i8[2] DOMAIN_ID:i8"},
    {"INPUT/ELEMENT/PLOTEL", "EID:i8 G:i8[2] DOMAIN_ID:i8"},
    {"INPUT/MATERIAL/MAT1",
     "MID:i8 E:f8 G:f8 NU:f8 RHO:f8 A:f8 TREF:f8 GE:f8 ST:f8 SC:f8 SS:f8 MCSID:i8 DOMAIN_ID:i8"},
    {"INPUT/MATERIAL/MAT2",
     "MID:i8 G11:f8 G12:f8 G13:f8 G22:f8 G23:f8 G33:f8 RHO:f8 A1:f8 A2:f8 A12:f8 TREF:f8 GE:f8 "
     "ST:f8 SC:f8 SS:f8 MCSID:i8 DOMAIN_ID:i8"},
    {"INPUT/MATERIAL/MAT4",
     "MID:i8 K:f8 CP:f8 RHO:f8 H:f8 MU:f8 HGEN:f8 REFENTH:f8 TCH:f8 TDELTA:f8 QLAT:f8 "
     "DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PBAR",
     "PID:i8 MID:i8 A:f8 I1:f8 I2:f8 J:f8 NSM:f8 FE:f8 C1:f8 C2:f8 D1:f8 D2:f8 E1:f8 E2:f8 F1:f8 "
     "F2:f8 K1:f8 K2:f8 I12:f8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PBEAM",
     "PID:i8 MID:i8 NSEGS:i8 CCF:i8 CWELD:i8 SO:f8[11] XXB:f8[11] A:f8[11] I1:f8[11] I2:f8[11] "
     "I12:f8[11] J:f8[11] NSM:f8[11] C1:f8[11] C2:f8[11] D1:f8[11] D2:f8[11] E1:f8[11] E2:f8[11] "
     "F1:f8[11] F2:f8[11] K1:f8 K2:f8 S1:f8 S2:f8 NSIA:f8 NSIB:f8 CWA:f8 CWB:f8 M1A:f8 M2A:f8 "
     "M1B:f8 M2B:f8 N1A:f8 N2A:f8 N1B:f8 N2B:f8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PBEND",
     "PID:i8 MID:i8 A:f8 I1:f8 I2:f8 J:f8 FSI:i8 RM:f8 T:f8 P:f8 RB:f8 THETAB:f8 C1:f8 C2:f8 D1:f8 "
     "D2:f8 E1:f8 E2:f8 F1:f8 F2:f8 K1:f8 K2:f8 NSM:f8 RC:f8 ZC:f8 DELTAN:i8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PBUSH1D",
     "PID:i8 K:f8 C:f8 M:f8 ALPHA:f8 SA:f8 EA:f8 TYPEA:i8 CVT:f8 CVC:f8 EXPVT:f8 EXPVC:f8 IDTSU:i8 "
     "IDTCU:i8 IDTSUD:i8 IDCSUD:i8 TYPES:i8 IDTS:i8 IDCS:i8 IDTDU1:i8 IDCDU1:i8 TYPED:i8 IDTD1:i8 "
     "IDTD2:i8 IDTDV1:i8 IDCDV1:i8 TYPEG:i8 IDTG:i8 IDCG:i8 IDTDU2:i8 IDCDU2:i8 IDTDV2:i8 "
     "IDCDV2:i8 TYPEF:i8 IDTF:i8 IDCF:i8 UT:f8 UC:f8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PCOMP/IDENTITY",
     "PID:i8 NPLIES:i8 Z0:f8 NSM:f8 SB:f8 FT:i8 TREF:f8 GE:f8 PLY_POS:i8 PLY_LEN:i8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PCOMP/PLY", "MID:i8 T:f8 THETA:f8 SOUT:i8"},
    {"INPUT/PROPERTY/PCONV",
     "PID:i8 MID:i8 FORM:i8 EXPF:f8 FTYPE:i8 HCF1:f8 HCF2:f8 HCF3:f8 HCF4:f8 HCF5:f8 HCF6:f8 "
     "HCF7:f8 HCF8:f8 TID:i8 CHLEN:f8 GIDIN:i8 CE:i8 E:f8[3] DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PDAMP", "PID:i8 B:f8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PELAS", "PID:i8 K:f8 GE:f8 S:f8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PELAST", "PID:i8 TKID:i8 TGEID:i8 TKNID:i8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PROD", "PID:i8 MID:i8 A:f8 J:f8 C:f8 NSM:f8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PSHEAR", "PID:i8 MID:i8 T:f8 NSM:f8 F1:f8 F2:f8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PSHELL",
     "PID:i8 MID1:i8 T:f8 MID2:i8 BK:f8 MID3:i8 TS:f8 NSM:f8 Z1:f8 Z2:f8 MID4:i8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PSOLID",
     "PID:i8 MID:i8 CORDM:i8 IN:i8 STRESS:i8 ISOP:i8 FCTN:S4 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PTUBE", "PID:i8 MID:i8 OD:f8 T:f8 NSM:f8 DOMAIN_ID:i8"},
    {"INPUT/PROPERTY/PVISC", "PID:i8 CE:f8 CR:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/BAR",
     "EID:i8 BM1A:f8 BM2A:f8 BM1B:f8 BM2B:f8 TS1:f8 TS2:f8 AF:f8 TRQ:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/BARS",
     "EID:i8 SD:f8 BM1:f8 BM2:f8 TS1:f8 TS2:f8 AF:f8 TRQ:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/BAR_CPLX",
     "EID:i8 BM1AR:f8 BM2AR:f8 BM1BR:f8 BM2BR:f8 TS1R:f8 TS2R:f8 AFR:f8 TRQR:f8 BM1AI:f8 BM2AI:f8 "
     "BM1BI:f8 BM2BI:f8 TS1I:f8 TS2I:f8 AFI:f8 TRQI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/BEAM",
     "EID:i8 GRID:i8[11] SD:f8[11] BM1:f8[11] BM2:f8[11] TS1:f8[11] TS2:f8[11] AF:f8[11] "
     "TTRQ:f8[11] WTRQ:f8[11] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/BEAM_CPLX",
     "EID:i8 GRID:i8[11] SD:f8[11] BM1R:f8[11] BM2R:f8[11] TS1R:f8[11] TS2R:f8[11] AFR:f8[11] "
     "TTRQR:f8[11] WTRQR:f8[11] BM1I:f8[11] BM2I:f8[11] TS1I:f8[11] TS2I:f8[11] AFI:f8[11] "
     "TTRQI:f8[11] WTRQI:f8[11] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/CONROD", "EID:i8 AF:f8 TRQ:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/CONROD_CPLX",
     "EID:i8 AFR:f8 AFI:f8 TRQR:f8 TRQI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP1_CPLX", "EID:i8 FR:f8 FI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP2_CPLX", "EID:i8 FR:f8 FI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP3_CPLX", "EID:i8 FR:f8 FI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/DAMP4_CPLX", "EID:i8 FR:f8 FI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1", "EID:i8 F:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS1_CPLX", "EID:i8 FR:f8 FI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS2", "EID:i8 F:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS2_CPLX", "EID:i8 FR:f8 FI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS3", "EID:i8 F:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS3_CPLX", "EID:i8 FR:f8 FI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS4", "EID:i8 F:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/ELAS4_CPLX", "EID:i8 FR:f8 FI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/GRAD_FLUX",
     "EID:i8 XGRAD:f8 YGRAD:f8 ZGRAD:f8 XFLUX:f8 YFLUX:f8 ZFLUX:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/HBDYE",
     "EID:i8 FAPPLIED:f8 FREECONV:f8 FORCECON:f8 FRAD:f8 FTOTAL:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_CN",
     "EID:i8 TERM:S4 GRID:i8[5] MX:f8[5] MY:f8[5] MXY:f8[5] BMX:f8[5] BMY:f8[5] BMXY:f8[5] "
     "TX:f8[5] TY:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD4_CN_CPLX",
     "EID:i8 TERM:S4 GRID:i8[5] MXR:f8[5] MYR:f8[5] MXYR:f8[5] BMXR:f8[5] BMYR:f8[5] BMXYR:f8[5] "
     "TXR:f8[5] TYR:f8[5] MXI:f8[5] MYI:f8[5] MXYI:f8[5] BMXI:f8[5] BMYI:f8[5] BMXYI:f8[5] "
     "TXI:f8[5] TYI:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD8",
     "EID:i8 TERM:S4 GRID:i8[5] MX:f8[5] MY:f8[5] MXY:f8[5] BMX:f8[5] BMY:f8[5] BMXY:f8[5] "
     "TX:f8[5] TY:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/QUAD8_CPLX",
     "EID:i8 TERM:S4 GRID:i8[5] MXR:f8[5] MYR:f8[5] MXYR:f8[5] BMXR:f8[5] BMYR:f8[5] BMXYR:f8[5] "
     "TXR:f8[5] TYR:f8[5] MXI:f8[5] MYI:f8[5] MXYI:f8[5] BMXI:f8[5] BMYI:f8[5] BMXYI:f8[5] "
     "TXI:f8[5] TYI:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/QUADR",
     "EID:i8 TERM:S4 GRID:i8[5] MX:f8[5] MY:f8[5] MXY:f8[5] BMX:f8[5] BMY:f8[5] BMXY:f8[5] "
     "TX:f8[5] TY:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/QUADR_CPLX",
     "EID:i8 TERM:S4 GRID:i8[5] MXR:f8[5] MYR:f8[5] MXYR:f8[5] BMXR:f8[5] BMYR:f8[5] BMXYR:f8[5] "
     "TXR:f8[5] TYR:f8[5] MXI:f8[5] MYI:f8[5] MXYI:f8[5] BMXI:f8[5] BMYI:f8[5] BMXYI:f8[5] "
     "TXI:f8[5] TYI:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/ROD", "EID:i8 AF:f8 TRQ:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/ROD_CPLX",
     "EID:i8 AFR:f8 AFI:f8 TRQR:f8 TRQI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/SHEAR",
     "EID:i8 F41:f8 F21:f8 F12:f8 F32:f8 F23:f8 F43:f8 F34:f8 F14:f8 KF1:f8 S12:f8 KF2:f8 S23:f8 "
     "KF3:f8 S34:f8 KF4:f8 S41:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/SHEAR_CPLX",
     "EID:i8 F41R:f8 F21R:f8 F12R:f8 F32R:f8 F23R:f8 F43R:f8 F34R:f8 F14R:f8 KF1R:f8 S12R:f8 "
     "KF2R:f8 S23R:f8 KF3R:f8 S34R:f8 KF4R:f8 S41R:f8 F41I:f8 F21I:f8 F12I:f8 F32I:f8 F23I:f8 "
     "F43I:f8 F34I:f8 F14I:f8 KF1I:f8 S12I:f8 KF2I:f8 S23I:f8 KF3I:f8 S34I:f8 KF4I:f8 S41I:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA3",
     "EID:i8 MX:f8 MY:f8 MXY:f8 BMX:f8 BMY:f8 BMXY:f8 TX:f8 TY:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA3_CPLX",
     "EID:i8 MXR:f8 MYR:f8 MXYR:f8 BMXR:f8 BMYR:f8 BMXYR:f8 TXR:f8 TYR:f8 MXI:f8 MYI:f8 MXYI:f8 "
     "BMXI:f8 BMYI:f8 BMXYI:f8 TXI:f8 TYI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA6",
     "EID:i8 TERM:S4 GRID:i8[4] MX:f8[4] MY:f8[4] MXY:f8[4] BMX:f8[4] BMY:f8[4] BMXY:f8[4] "
     "TX:f8[4] TY:f8[4] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/TRIA6_CPLX",
     "EID:i8 TERM:S4 GRID:i8[4] MXR:f8[4] MYR:f8[4] MXYR:f8[4] BMXR:f8[4] BMYR:f8[4] BMXYR:f8[4] "
     "TXR:f8[4] TYR:f8[4] MXI:f8[4] MYI:f8[4] MXYI:f8[4] BMXI:f8[4] BMYI:f8[4] BMXYI:f8[4] "
     "TXI:f8[4] TYI:f8[4] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/TRIAR",
     "EID:i8 TERM:S4 GRID:i8[4] MX:f8[4] MY:f8[4] MXY:f8[4] BMX:f8[4] BMY:f8[4] BMXY:f8[4] "
     "TX:f8[4] TY:f8[4] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/TRIAR_CPLX",
     "EID:i8 TERM:S4 GRID:i8[4] MXR:f8[4] MYR:f8[4] MXYR:f8[4] BMXR:f8[4] BMYR:f8[4] BMXYR:f8[4] "
     "TXR:f8[4] TYR:f8[4] MXI:f8[4] MYI:f8[4] MXYI:f8[4] BMXI:f8[4] BMYI:f8[4] BMXYI:f8[4] "
     "TXI:f8[4] TYI:f8[4] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/TUBE", "EID:i8 AF:f8 TRQ:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/TUBE_CPLX",
     "EID:i8 AFR:f8 AFI:f8 TRQR:f8 TRQI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ELEMENT_FORCE/VISC_CPLX",
     "EID:i8 AFR:f8 AFI:f8 TRQR:f8 TRQI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/ENERGY/IDENT",
     "IDENT:i8 ELNAME:S8 ETOTAL:f8 CVALRES:i8 ESUBT:f8 ETOTPOS:f8 ETOTNEG:f8"},
    {"RESULT/ELEMENTAL/ENERGY/STRAIN_ELEM", "ID:i8 ENERGY:f8 PCT:f8 DEN:f8 IDENT:i8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/BAR",
     "EID:i8 X1A:f8 X2A:f8 X3A:f8 X4A:f8 AX:f8 MAXA:f8 MINA:f8 MST:f8 X1B:f8 X2B:f8 X3B:f8 X4B:f8 "
     "MAXB:f8 MINB:f8 MSC:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/BARS",
     "EID:i8 SD:f8 XC:f8 XD:f8 XE:f8 XF:f8 AX:f8 MAX:f8 MIN:f8 MS:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/BAR_CPLX",
     "EID:i8 X1AR:f8 X2AR:f8 X3AR:f8 X4AR:f8 AXR:f8 X1AI:f8 X2AI:f8 X3AI:f8 X4AI:f8 AXI:f8 X1BR:f8 "
     "X2BR:f8 X3BR:f8 X4BR:f8 X1BI:f8 X2BI:f8 X3BI:f8 X4BI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/BEAM",
     "EID:i8 GRID:i8[11] SD:f8[11] XC:f8[11] XD:f8[11] XE:f8[11] XF:f8[11] MAX:f8[11] MIN:f8[11] "
     "MST:f8[11] MSC:f8[11] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/BEAM_CPLX",
     "EID:i8 GRID:i8[11] SD:f8[11] XCR:f8[11] XDR:f8[11] XER:f8[11] XFR:f8[11] XCI:f8[11] "
     "XDI:f8[11] XEI:f8[11] XFI:f8[11] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/CONROD", "EID:i8 A:f8 MSA:f8 T:f8 MST:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/CONROD_CPLX", "EID:i8 AR:f8 AI:f8 TR:f8 TI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/ELAS1", "EID:i8 S:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/ELAS1_CPLX", "EID:i8 SR:f8 SI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/ELAS2", "EID:i8 S:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/ELAS2_CPLX", "EID:i8 SR:f8 SI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/ELAS3", "EID:i8 S:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/ELAS3_CPLX", "EID:i8 SR:f8 SI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/HEXA",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[9] X:f8[9] Y:f8[9] Z:f8[9] TXY:f8[9] TYZ:f8[9] "
     "TZX:f8[9] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/HEXA_CPLX",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[9] XR:f8[9] YR:f8[9] ZR:f8[9] TXYR:f8[9] TYZR:f8[9] "
     "TZXR:f8[9] XI:f8[9] YI:f8[9] ZI:f8[9] TXYI:f8[9] TYZI:f8[9] TZXI:f8[9] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/PENTA",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[7] X:f8[7] Y:f8[7] Z:f8[7] TXY:f8[7] TYZ:f8[7] "
     "TZX:f8[7] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/PENTA_CPLX",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[7] XR:f8[7] YR:f8[7] ZR:f8[7] TXYR:f8[7] TYZR:f8[7] "
     "TZXR:f8[7] XI:f8[7] YI:f8[7] ZI:f8[7] TXYI:f8[7] TYZI:f8[7] TZXI:f8[7] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/QUAD4_COMP",
     "EID:i8 PLY:i8 X1:f8 Y1:f8 T1:f8 L1:f8 L2:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/QUAD4_COMP_CPLX",
     "EID:i8 PLY:i8 X1R:f8 Y1R:f8 T1R:f8 L1R:f8 L2R:f8 X1I:f8 Y1I:f8 T1I:f8 L1I:f8 L2I:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/QUAD8",
     "EID:i8 TERM:S4 GRID:i8[5] FD1:f8[5] X1:f8[5] Y1:f8[5] TXY1:f8[5] FD2:f8[5] X2:f8[5] Y2:f8[5] "
     "TXY2:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/QUAD8_CPLX",
     "EID:i8 TERM:S4 GRID:i8[5] FD1:f8[5] X1R:f8[5] X1I:f8[5] Y1R:f8[5] Y1I:f8[5] TXY1R:f8[5] "
     "TXY1I:f8[5] FD2:f8[5] X2R:f8[5] X2I:f8[5] Y2R:f8[5] Y2I:f8[5] TXY2R:f8[5] TXY2I:f8[5] "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/QUADR_COMP",
     "EID:i8 PLY:i8 X1:f8 Y1:f8 T1:f8 L1:f8 L2:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/QUADR_COMP_CPLX",
     "EID:i8 PLY:i8 X1R:f8 Y1R:f8 T1R:f8 L1R:f8 L2R:f8 X1I:f8 Y1I:f8 T1I:f8 L1I:f8 L2I:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/QUADR_CPLX",
     "EID:i8 TERM:S4 GRID:i8[5] FD1:f8[5] X1R:f8[5] X1I:f8[5] Y1R:f8[5] Y1I:f8[5] TXY1R:f8[5] "
     "TXY1I:f8[5] FD2:f8[5] X2R:f8[5] X2I:f8[5] Y2R:f8[5] Y2I:f8[5] TXY2R:f8[5] TXY2I:f8[5] "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/QUAD_CN",
     "EID:i8 TERM:S4 GRID:i8[5] FD1:f8[5] X1:f8[5] Y1:f8[5] TXY1:f8[5] FD2:f8[5] X2:f8[5] Y2:f8[5] "
     "TXY2:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/QUAD_CN_CPLX",
     "EID:i8 TERM:S4 GRID:i8[5] FD1:f8[5] X1R:f8[5] X1I:f8[5] Y1R:f8[5] Y1I:f8[5] TXY1R:f8[5] "
     "TXY1I:f8[5] FD2:f8[5] X2R:f8[5] X2I:f8[5] Y2R:f8[5] Y2I:f8[5] TXY2R:f8[5] TXY2I:f8[5] "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/ROD", "EID:i8 A:f8 MSA:f8 T:f8 MST:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/ROD_CPLX", "EID:i8 AR:f8 AI:f8 TR:f8 TI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TETRA",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[5] X:f8[5] Y:f8[5] Z:f8[5] TXY:f8[5] TYZ:f8[5] "
     "TZX:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TETRA_CPLX",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[5] XR:f8[5] YR:f8[5] ZR:f8[5] TXYR:f8[5] TYZR:f8[5] "
     "TZXR:f8[5] XI:f8[5] YI:f8[5] ZI:f8[5] TXYI:f8[5] TYZI:f8[5] TZXI:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TRIA3",
     "EID:i8 FD1:f8 X1:f8 Y1:f8 TXY1:f8 FD2:f8 X2:f8 Y2:f8 TXY2:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TRIA3_COMP",
     "EID:i8 PLY:i8 X1:f8 Y1:f8 T1:f8 L1:f8 L2:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TRIA3_COMP_CPLX",
     "EID:i8 PLY:i8 X1R:f8 Y1R:f8 T1R:f8 L1R:f8 L2R:f8 X1I:f8 Y1I:f8 T1I:f8 L1I:f8 L2I:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TRIA3_CPLX",
     "EID:i8 FD1:f8 X1R:f8 X1I:f8 Y1R:f8 Y1I:f8 TXY1R:f8 TXY1I:f8 FD2:f8 X2R:f8 X2I:f8 Y2R:f8 "
     "Y2I:f8 TXY2R:f8 TXY2I:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TRIA6",
     "EID:i8 TERM:S4 GRID:i8[4] FD1:f8[4] X1:f8[4] Y1:f8[4] TXY1:f8[4] FD2:f8[4] X2:f8[4] Y2:f8[4] "
     "TXY2:f8[4] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TRIA6_CPLX",
     "EID:i8 TERM:S4 GRID:i8[4] FD1:f8[4] X1R:f8[4] X1I:f8[4] Y1R:f8[4] Y1I:f8[4] TXY1R:f8[4] "
     "TXY1I:f8[4] FD2:f8[4] X2R:f8[4] X2I:f8[4] Y2R:f8[4] Y2I:f8[4] TXY2R:f8[4] TXY2I:f8[4] "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TRIAR_COMP",
     "EID:i8 PLY:i8 X1:f8 Y1:f8 T1:f8 L1:f8 L2:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TRIAR_COMP_CPLX",
     "EID:i8 PLY:i8 X1R:f8 Y1R:f8 T1R:f8 L1R:f8 L2R:f8 X1I:f8 Y1I:f8 T1I:f8 L1I:f8 L2I:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TRIAR_CPLX",
     "EID:i8 TERM:S4 GRID:i8[4] FD1:f8[4] X1R:f8[4] X1I:f8[4] Y1R:f8[4] Y1I:f8[4] TXY1R:f8[4] "
     "TXY1I:f8[4] FD2:f8[4] X2R:f8[4] X2I:f8[4] Y2R:f8[4] Y2I:f8[4] TXY2R:f8[4] TXY2I:f8[4] "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TUBE", "EID:i8 AS:f8 MSA:f8 TS:f8 MST:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRAIN/TUBE_CPLX", "EID:i8 ASR:f8 ASI:f8 TSR:f8 TSI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/BAR",
     "EID:i8 X1A:f8 X2A:f8 X3A:f8 X4A:f8 AX:f8 MAXA:f8 MINA:f8 MST:f8 X1B:f8 X2B:f8 X3B:f8 X4B:f8 "
     "MAXB:f8 MINB:f8 MSC:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/BARS",
     "EID:i8 SD:f8 XC:f8 XD:f8 XE:f8 XF:f8 AX:f8 MAX:f8 MIN:f8 MS:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/BAR_CPLX",
     "EID:i8 X1AR:f8 X2AR:f8 X3AR:f8 X4AR:f8 AXR:f8 X1AI:f8 X2AI:f8 X3AI:f8 X4AI:f8 AXI:f8 X1BR:f8 "
     "X2BR:f8 X3BR:f8 X4BR:f8 X1BI:f8 X2BI:f8 X3BI:f8 X4BI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/BEAM",
     "EID:i8 GRID:i8[11] SD:f8[11] XC:f8[11] XD:f8[11] XE:f8[11] XF:f8[11] MAX:f8[11] MIN:f8[11] "
     "MST:f8[11] MSC:f8[11] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/BEAM_CPLX",
     "EID:i8 GRID:i8[11] SD:f8[11] XCR:f8[11] XDR:f8[11] XER:f8[11] XFR:f8[11] XCI:f8[11] "
     "XDI:f8[11] XEI:f8[11] XFI:f8[11] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/CONROD", "EID:i8 A:f8 MSA:f8 T:f8 MST:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/CONROD_CPLX", "EID:i8 AR:f8 AI:f8 TR:f8 TI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/ELAS1", "EID:i8 S:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/ELAS1_CPLX", "EID:i8 SR:f8 SI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/ELAS2", "EID:i8 S:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/ELAS2_CPLX", "EID:i8 SR:f8 SI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/ELAS3", "EID:i8 S:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/ELAS3_CPLX", "EID:i8 SR:f8 SI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/HEXA",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[9] X:f8[9] Y:f8[9] Z:f8[9] TXY:f8[9] TYZ:f8[9] "
     "TZX:f8[9] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/HEXA_CPLX",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[9] XR:f8[9] YR:f8[9] ZR:f8[9] TXYR:f8[9] TYZR:f8[9] "
     "TZXR:f8[9] XI:f8[9] YI:f8[9] ZI:f8[9] TXYI:f8[9] TYZI:f8[9] TZXI:f8[9] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/PENTA",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[7] X:f8[7] Y:f8[7] Z:f8[7] TXY:f8[7] TYZ:f8[7] "
     "TZX:f8[7] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/PENTA_CPLX",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[7] XR:f8[7] YR:f8[7] ZR:f8[7] TXYR:f8[7] TYZR:f8[7] "
     "TZXR:f8[7] XI:f8[7] YI:f8[7] ZI:f8[7] TXYI:f8[7] TYZI:f8[7] TZXI:f8[7] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/QUAD4_COMP",
     "EID:i8 PLY:i8 X1:f8 Y1:f8 T1:f8 L1:f8 L2:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/QUAD4_COMP_CPLX",
     "EID:i8 PLY:i8 X1R:f8 Y1R:f8 T1R:f8 L1R:f8 L2R:f8 X1I:f8 Y1I:f8 T1I:f8 L1I:f8 L2I:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/QUAD8",
     "EID:i8 TERM:S4 GRID:i8[5] FD1:f8[5] X1:f8[5] Y1:f8[5] TXY1:f8[5] FD2:f8[5] X2:f8[5] Y2:f8[5] "
     "TXY2:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/QUAD8_CPLX",
     "EID:i8 TERM:S4 GRID:i8[5] FD1:f8[5] X1R:f8[5] X1I:f8[5] Y1R:f8[5] Y1I:f8[5] TXY1R:f8[5] "
     "TXY1I:f8[5] FD2:f8[5] X2R:f8[5] X2I:f8[5] Y2R:f8[5] Y2I:f8[5] TXY2R:f8[5] TXY2I:f8[5] "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/QUADR_COMP",
     "EID:i8 PLY:i8 X1:f8 Y1:f8 T1:f8 L1:f8 L2:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/QUADR_COMP_CPLX",
     "EID:i8 PLY:i8 X1R:f8 Y1R:f8 T1R:f8 L1R:f8 L2R:f8 X1I:f8 Y1I:f8 T1I:f8 L1I:f8 L2I:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/QUADR_CPLX",
     "EID:i8 TERM:S4 GRID:i8[5] FD1:f8[5] X1R:f8[5] X1I:f8[5] Y1R:f8[5] Y1I:f8[5] TXY1R:f8[5] "
     "TXY1I:f8[5] FD2:f8[5] X2R:f8[5] X2I:f8[5] Y2R:f8[5] Y2I:f8[5] TXY2R:f8[5] TXY2I:f8[5] "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/QUAD_CN",
     "EID:i8 TERM:S4 GRID:i8[5] FD1:f8[5] X1:f8[5] Y1:f8[5] TXY1:f8[5] FD2:f8[5] X2:f8[5] Y2:f8[5] "
     "TXY2:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/QUAD_CN_CPLX",
     "EID:i8 TERM:S4 GRID:i8[5] FD1:f8[5] X1R:f8[5] X1I:f8[5] Y1R:f8[5] Y1I:f8[5] TXY1R:f8[5] "
     "TXY1I:f8[5] FD2:f8[5] X2R:f8[5] X2I:f8[5] Y2R:f8[5] Y2I:f8[5] TXY2R:f8[5] TXY2I:f8[5] "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/ROD", "EID:i8 A:f8 MSA:f8 T:f8 MST:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/ROD_CPLX", "EID:i8 AR:f8 AI:f8 TR:f8 TI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/SHEAR", "EID:i8 TMAX:f8 TAVG:f8 MS:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/SHEAR_CPLX",
     "EID:i8 TMAXR:f8 TMAXI:f8 TAVGR:f8 TAVGI:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TETRA",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[5] X:f8[5] Y:f8[5] Z:f8[5] TXY:f8[5] TYZ:f8[5] "
     "TZX:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TETRA_CPLX",
     "EID:i8 CID:i8 CTYPE:S4 NODEF:i8 GRID:i8[5] XR:f8[5] YR:f8[5] ZR:f8[5] TXYR:f8[5] TYZR:f8[5] "
     "TZXR:f8[5] XI:f8[5] YI:f8[5] ZI:f8[5] TXYI:f8[5] TYZI:f8[5] TZXI:f8[5] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TRIA3",
     "EID:i8 FD1:f8 X1:f8 Y1:f8 TXY1:f8 FD2:f8 X2:f8 Y2:f8 TXY2:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TRIA3_COMP",
     "EID:i8 PLY:i8 X1:f8 Y1:f8 T1:f8 L1:f8 L2:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TRIA3_COMP_CPLX",
     "EID:i8 PLY:i8 X1R:f8 Y1R:f8 T1R:f8 L1R:f8 L2R:f8 X1I:f8 Y1I:f8 T1I:f8 L1I:f8 L2I:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TRIA3_CPLX",
     "EID:i8 FD1:f8 X1R:f8 X1I:f8 Y1R:f8 Y1I:f8 TXY1R:f8 TXY1I:f8 FD2:f8 X2R:f8 X2I:f8 Y2R:f8 "
     "Y2I:f8 TXY2R:f8 TXY2I:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TRIA6",
     "EID:i8 TERM:S4 GRID:i8[4] FD1:f8[4] X1:f8[4] Y1:f8[4] TXY1:f8[4] FD2:f8[4] X2:f8[4] Y2:f8[4] "
     "TXY2:f8[4] DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TRIA6_CPLX",
     "EID:i8 TERM:S4 GRID:i8[4] FD1:f8[4] X1R:f8[4] X1I:f8[4] Y1R:f8[4] Y1I:f8[4] TXY1R:f8[4] "
     "TXY1I:f8[4] FD2:f8[4] X2R:f8[4] X2I:f8[4] Y2R:f8[4] Y2I:f8[4] TXY2R:f8[4] TXY2I:f8[4] "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TRIAR_COMP",
     "EID:i8 PLY:i8 X1:f8 Y1:f8 T1:f8 L1:f8 L2:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TRIAR_COMP_CPLX",
     "EID:i8 PLY:i8 X1R:f8 Y1R:f8 T1R:f8 L1R:f8 L2R:f8 X1I:f8 Y1I:f8 T1I:f8 L1I:f8 L2I:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TRIAR_CPLX",
     "EID:i8 TERM:S4 GRID:i8[4] FD1:f8[4] X1R:f8[4] X1I:f8[4] Y1R:f8[4] Y1I:f8[4] TXY1R:f8[4] "
     "TXY1I:f8[4] FD2:f8[4] X2R:f8[4] X2I:f8[4] Y2R:f8[4] Y2I:f8[4] TXY2R:f8[4] TXY2I:f8[4] "
     "DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TUBE", "EID:i8 AS:f8 MSA:f8 TS:f8 MST:f8 DOMAIN_ID:i8"},
    {"RESULT/ELEMENTAL/STRESS/TUBE_CPLX", "EID:i8 ASR:f8 ASI:f8 TSR:f8 TSI:f8 DOMAIN_ID:i8"},
    {"RESULT/NODAL/APPLIED_LOAD", "ID:i8 X:f8 Y:f8 Z:f8 RX:f8 RY:f8 RZ:f8 DOMAIN_ID:i8"},
    {"RESULT/NODAL/APPLIED_LOAD_CPLX",
     "ID:i8 XR:f8 YR:f8 ZR:f8 RXR:f8 RYR:f8 RZR:f8 XI:f8 YI:f8 ZI:f8 RXI:f8 RYI:f8 RZI:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/NODAL/DISPLACEMENT", "ID:i8 X:f8 Y:f8 Z:f8 RX:f8 RY:f8 RZ:f8 DOMAIN_ID:i8"},
    {"RESULT/NODAL/DISPLACEMENT_CPLX",
     "ID:i8 XR:f8 YR:f8 ZR:f8 RXR:f8 RYR:f8 RZR:f8 XI:f8 YI:f8 ZI:f8 RXI:f8 RYI:f8 RZI:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/NODAL/EIGENVECTOR", "ID:i8 X:f8 Y:f8 Z:f8 RX:f8 RY:f8 RZ:f8 DOMAIN_ID:i8"},
    {"RESULT/NODAL/EIGENVECTOR_CPLX",
     "ID:i8 XR:f8 YR:f8 ZR:f8 RXR:f8 RYR:f8 RZR:f8 XI:f8 YI:f8 ZI:f8 RXI:f8 RYI:f8 RZI:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/NODAL/GRID_FORCE",
     "ID:i8 EID:i8 ELNAME:S8 F1:f8 F2:f8 F3:f8 M1:f8 M2:f8 M3:f8 DOMAIN_ID:i8"},
    {"RESULT/NODAL/KINETIC_ENERGY",
     "ID:i8 TYPE:i8 KET1:f8 KET2:f8 KET3:f8 KER1:f8 KER2:f8 KER3:f8 DOMAIN_ID:i8"},
    {"RESULT/NODAL/MPC_FORCE", "ID:i8 X:f8 Y:f8 Z:f8 RX:f8 RY:f8 RZ:f8 DOMAIN_ID:i8"},
    {"RESULT/NODAL/SPC_FORCE", "ID:i8 X:f8 Y:f8 Z:f8 RX:f8 RY:f8 RZ:f8 DOMAIN_ID:i8"},
    {"RESULT/NODAL/SPC_FORCE_CPLX",
     "ID:i8 XR:f8 YR:f8 ZR:f8 RXR:f8 RYR:f8 RZR:f8 XI:f8 YI:f8 ZI:f8 RXI:f8 RYI:f8 RZI:f8 "
     "DOMAIN_ID:i8"},
    {"RESULT/NODAL/TEMPERATURE", "ID:i8 VALUE:f8 DOMAIN_ID:i8"},
    {"RESULT/NODAL/VELOCITY", "ID:i8 X:f8 Y:f8 Z:f8 RX:f8 RY:f8 RZ:f8 DOMAIN_ID:i8"},
};
// END GENERATED SCHEMA

enum Nh5wProp {
    kNh5wNoProp = -1,
    kNh5wProd = 0,
    kNh5wPshell,
    kNh5wPsolid,
    kNh5wPtube,
    kNh5wPbar,
    kNh5wPbeam,
    kNh5wPelas,
    kNh5wPdamp,
    kNh5wPvisc,
    kNh5wPshear,
    kNh5wPcomp,
    kNh5wNumProps
};
// The property cards' names, in the order of `Nh5wProp`.
constexpr const char* kNh5wPropNames[] = {"PROD",  "PSHELL", "PSOLID", "PTUBE",  "PBAR", "PBEAM",
                                          "PELAS", "PDAMP",  "PVISC",  "PSHEAR", "PCOMP"};
// The order a clash between two families' property ids is reported in.
constexpr int kNh5wUsedOrder[] = {kNh5wPshell, kNh5wPsolid, kNh5wProd,  kNh5wPtube,
                                  kNh5wPbar,   kNh5wPbeam,  kNh5wPelas, kNh5wPdamp,
                                  kNh5wPvisc,  kNh5wPshear, kNh5wPcomp};

constexpr int kNh5wHexa20[20] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                 10, 11, 16, 17, 18, 19, 12, 13, 14, 15};
constexpr int kNh5wPenta15[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 9, 10, 11};

struct Nh5wCard {
    const char* mName;
    std::int64_t mCode;   // the card's type number: `nastran:card`
    std::size_t mLinear;  // nodes of the linear card
    std::size_t mWidth;   // columns of G
    const int* mPerm;     // quadratic cell -> Nastran's order
    Nh5wProp mProp;
    const char* mFamily;  // result tables' name; null: no results
    bool mHasPid;
};

// In card-name order.
constexpr Nh5wCard kNh5wCards[] = {
    {"CBAR", 24, 2, 2, nullptr, kNh5wPbar, "BARS", true},
    {"CBEAM", 54, 2, 2, nullptr, kNh5wPbeam, "BEAM", true},
    {"CDAMP1", 2, 2, 2, nullptr, kNh5wPdamp, nullptr, true},
    {"CDAMP2", 3, 2, 2, nullptr, kNh5wNoProp, nullptr, false},
    {"CELAS1", 6, 2, 2, nullptr, kNh5wPelas, "ELAS1", true},
    {"CELAS2", 7, 2, 2, nullptr, kNh5wNoProp, "ELAS2", false},
    {"CHEXA", 73, 8, 20, kNh5wHexa20, kNh5wPsolid, "HEXA", true},
    {"CONM2", 15, 1, 1, nullptr, kNh5wNoProp, nullptr, false},
    {"CONROD", 16, 2, 2, nullptr, kNh5wNoProp, "CONROD", false},
    {"CPENTA", 41, 6, 15, kNh5wPenta15, kNh5wPsolid, "PENTA", true},
    {"CPYRAM", 172, 5, 13, nullptr, kNh5wPsolid, "PYRAM", true},
    {"CQUAD", 91, 4, 9, nullptr, kNh5wPshell, "QUADR", true},
    {"CQUAD4", 51, 4, 4, nullptr, kNh5wPshell, "QUAD4", true},
    {"CQUAD8", 47, 4, 8, nullptr, kNh5wPshell, "QUAD8", true},
    {"CROD", 30, 2, 2, nullptr, kNh5wProd, "ROD", true},
    {"CSHEAR", 31, 4, 4, nullptr, kNh5wPshear, "SHEAR", true},
    {"CTETRA", 55, 4, 10, nullptr, kNh5wPsolid, "TETRA", true},
    {"CTRIA3", 59, 3, 3, nullptr, kNh5wPshell, "TRIA3", true},
    {"CTRIA6", 48, 3, 6, nullptr, kNh5wPshell, "TRIA6", true},
    {"CTUBE", 37, 2, 2, nullptr, kNh5wPtube, "TUBE", true},
    {"CVISC", 39, 2, 2, nullptr, kNh5wPvisc, nullptr, true},
    {"PLOTEL", 52, 2, 2, nullptr, kNh5wNoProp, nullptr, false},
};
constexpr std::size_t kNh5wNumCards = sizeof(kNh5wCards) / sizeof(kNh5wCards[0]);

std::size_t nh5w_card_index(std::string_view Name) {
    for (std::size_t i = 0; i < kNh5wNumCards; ++i)
        if (Name == kNh5wCards[i].mName)
            return i;
    return kNh5wNumCards;
}

/// The card a line or vertex takes from its `nastran:card` number, or nullptr when the OP2
/// writer holds no cell type for it. Cards the HDF5 writer has no table for come back too
/// (CBUSH, CGAP, CMASS1, CMASS2, CONM1); `nh5w_card_index` tells them apart.
const char* nh5w_line_card(std::int64_t Code) {
    switch (Code) {
        case 30:
            return "CROD";
        case 16:
            return "CONROD";
        case 37:
            return "CTUBE";
        case 6:
            return "CELAS1";
        case 2:
            return "CDAMP1";
        case 15:
            return "CONM2";
        case 24:
            return "CBAR";
        case 54:
            return "CBEAM";
        case 26:
            return "CBUSH";
        case 3:
            return "CDAMP2";
        case 7:
            return "CELAS2";
        case 39:
            return "CVISC";
        case 52:
            return "PLOTEL";
        case 19:
            return "CGAP";
        case 10:
            return "CMASS1";
        case 11:
            return "CMASS2";
        case 14:
            return "CONM1";
        default:
            return nullptr;
    }
}

/// Whether a cell block of @p rType may hold an element of card @p Card.
bool nh5w_block_holds(const std::string& rType, std::string_view Card) {
    if (rType == "vertex")
        return Card == "CONM2" || Card == "CELAS1" || Card == "CDAMP1" || Card == "CELAS2" ||
               Card == "CDAMP2" || Card == "CMASS1" || Card == "CMASS2" || Card == "CONM1";
    if (rType == "line")
        return Card == "CROD" || Card == "CONROD" || Card == "CTUBE" || Card == "CELAS1" ||
               Card == "CDAMP1" || Card == "CBAR" || Card == "CBEAM" || Card == "CBUSH" ||
               Card == "CELAS2" || Card == "CDAMP2" || Card == "CVISC" || Card == "PLOTEL" ||
               Card == "CGAP" || Card == "CMASS1" || Card == "CMASS2";
    return false;
}

/// The card a `nastran:card` number names (its GEOM2 type), for the warnings.
std::string nh5w_card_name(std::int64_t Code) {
    static const std::pair<const char*, std::int64_t> kTypes[] = {
        {"CBAR", 24},    {"CBEAM", 54},  {"CBUSH", 26},  {"CDAMP1", 2},  {"CDAMP2", 3},
        {"CELAS1", 6},   {"CELAS2", 7},  {"CGAP", 19},   {"CHEXA", 73},  {"CMASS1", 10},
        {"CMASS2", 11},  {"CONM1", 14},  {"CONM2", 15},  {"CONROD", 16}, {"CPENTA", 41},
        {"CPYRAM", 172}, {"CQUAD", 91},  {"CQUAD4", 51}, {"CQUAD8", 47}, {"CQUADR", 80},
        {"CROD", 30},    {"CSHEAR", 31}, {"CTETRA", 55}, {"CTRIA3", 59}, {"CTRIA6", 48},
        {"CTRIAR", 92},  {"CTUBE", 37},  {"CVISC", 39},  {"PLOTEL", 52},
    };
    for (const auto& [name, number] : kTypes)
        if (number == Code)
            return name;
    return "card " + std::to_string(Code);
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

/// The columns of a row's card fields: a CBAR's X1..X3, G0, PA, PB and W1A..W3B (12), a CELAS2's
/// K, GE and S, a CDAMP2's B (the card's first columns).
constexpr std::size_t kNh5wBarWidth = 12;

struct Nh5wBlock {
    std::size_t mIndex = 0;
    std::size_t mCount = 0;
    std::size_t mNodes = 0;
    std::vector<std::int64_t> mConn;  // 0-based, mCount x mNodes
    std::vector<std::int64_t> mEids;
    std::vector<std::int64_t> mPids;
    std::vector<std::size_t> mCards;
    std::vector<std::int64_t> mComps;  // two per cell: a spring's or damper's components
    std::vector<double> mBar;          // kNh5wBarWidth per cell
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

/// The GRID id of each point: `point_data["nastran:nid"]`, else 1 to n.
std::vector<std::int64_t> nh5w_node_ids(const Mesh& rMesh, std::size_t N) {
    std::vector<std::int64_t> ids(N);
    if (!rMesh.HasPointData("nastran:nid")) {
        for (std::size_t i = 0; i < N; ++i)
            ids[i] = static_cast<std::int64_t>(i + 1);
        return ids;
    }
    const NDArray& a = rMesh.PointData("nastran:nid");
    if (a.Size() != N || detail::is_float_dtype(a.Dtype()))
        nh5w_fail("'nastran:nid' must hold one integer per point");
    const detail::Int64View v(a);
    ids.assign(v.Data(), v.Data() + N);
    std::vector<std::int64_t> sorted = ids;
    std::sort(sorted.begin(), sorted.end());
    if (!sorted.empty() && sorted.front() < 1)
        nh5w_fail("GRID ids must be positive");
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
        nh5w_fail("GRID ids are not unique");
    return ids;
}

using Nh5wValues = std::map<std::pair<std::string, std::int64_t>, std::map<std::string, double>>;

/// `{(card, id): {member: value}}` from the `nastran:<CARD>_<id>:<MEMBER>` field data; keys
/// that name no card or member here are left alone.
Nh5wValues nh5w_property_values(const Mesh& rMesh) {
    Nh5wValues out;
    for (const std::string& name : rMesh.FieldDataNames()) {
        if (name.rfind("nastran:", 0) != 0)
            continue;
        const std::size_t colon = name.rfind(':');
        const std::size_t under = name.rfind('_', colon);
        if (colon <= 8 || under == std::string::npos || under < 8 || under + 1 >= colon)
            continue;
        const std::string card = name.substr(8, under - 8);
        const std::string member = name.substr(colon + 1);
        const std::string digits = name.substr(under + 1, colon - under - 1);
        if (digits.find_first_not_of("0123456789") != std::string::npos || digits.size() > 15)
            continue;
        for (const Nh5ValueCard& c : nh5_value_cards()) {
            if (card != c.mCard)
                continue;
            for (const Nh5Member& m : c.mMembers)
                if (member == m.mName) {
                    const NDArray& a = rMesh.FieldData(name);
                    if (a.Size() != 1)
                        nh5w_fail("'" + name + "' must hold one value");
                    out[{card, std::stoll(digits)}][member] =
                        m.mInt ? static_cast<double>(detail::Int64View(a).Data()[0])
                               : detail::DoubleView(a).Data()[0];
                }
        }
    }
    return out;
}

/// The PCOMP scalars the field data gives (`Z0`, `NSM`, `SB`, `FT`, `TREF`, `GE`), by id.
Nh5wValues nh5w_pcomp_values(const Mesh& rMesh) {
    Nh5wValues out;
    for (const std::string& name : rMesh.FieldDataNames()) {
        static const std::string prefix = "nastran:PCOMP_";
        if (name.rfind(prefix, 0) != 0)
            continue;
        const std::size_t colon = name.find(':', prefix.size());
        if (colon == std::string::npos)
            continue;
        const std::string digits = name.substr(prefix.size(), colon - prefix.size());
        const std::string member = name.substr(colon + 1);
        if (digits.empty() || digits.size() > 15 ||
            digits.find_first_not_of("0123456789") != std::string::npos)
            continue;
        if (member != "Z0" && member != "NSM" && member != "SB" && member != "FT" &&
            member != "TREF" && member != "GE")
            continue;
        const NDArray& a = rMesh.FieldData(name);
        if (a.Size() != 1)
            nh5w_fail("'" + name + "' must hold one value");
        out[{"PCOMP", std::stoll(digits)}][member] =
            member == "FT" ? static_cast<double>(detail::Int64View(a).Data()[0])
                           : detail::DoubleView(a).Data()[0];
    }
    return out;
}

using Nh5wPlies = std::map<std::int64_t, std::vector<std::array<double, 4>>>;

/// The plies of each PCOMP id: `nastran:PCOMP_<id>:PLY`, `(n, 4)` rows of MID, T, THETA, SOUT.
Nh5wPlies nh5w_plies(const Mesh& rMesh) {
    Nh5wPlies out;
    const std::string prefix = "nastran:PCOMP_";
    const std::string suffix = ":PLY";
    for (const std::string& name : rMesh.FieldDataNames()) {
        if (name.size() <= prefix.size() + suffix.size() || name.rfind(prefix, 0) != 0 ||
            name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
            continue;
        const std::string digits =
            name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
        if (digits.empty() || digits.size() > 15 ||
            digits.find_first_not_of("0123456789") != std::string::npos)
            continue;
        const NDArray& a = rMesh.FieldData(name);
        const std::vector<std::size_t>& shape = a.Shape();
        if (shape.size() != 2 || shape[1] != 4 || shape[0] == 0)
            nh5w_fail("'" + name + "' must hold (plies, 4) rows of MID, T, THETA and SOUT");
        const detail::DoubleView view(a);
        auto& plies = out[std::stoll(digits)];
        for (std::size_t r = 0; r < shape[0]; ++r) {
            const double* row = view.Data() + 4 * r;
            plies.push_back({row[0], row[1], row[2], row[3]});
        }
    }
    return out;
}

using Nh5wArrays = std::map<std::int64_t, std::map<std::string, std::vector<double>>>;

/// `nastran:PBEAM_<id>:<MEMBER>` field data by id, each member flattened.
Nh5wArrays nh5w_beam_arrays(const Mesh& rMesh) {
    const std::string prefix = "nastran:PBEAM_";
    Nh5wArrays out;
    for (const std::string& name : rMesh.FieldDataNames()) {
        if (name.rfind(prefix, 0) != 0)
            continue;
        const std::size_t colon = name.find(':', prefix.size());
        if (colon == std::string::npos || colon == prefix.size())
            continue;
        const std::string digits = name.substr(prefix.size(), colon - prefix.size());
        const std::string member = name.substr(colon + 1);
        if (digits.find_first_not_of("0123456789") != std::string::npos || digits.size() > 15)
            continue;
        if (member.empty() ||
            member.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") != std::string::npos)
            continue;
        const NDArray& a = rMesh.FieldData(name);
        std::vector<double> values(a.Size());
        if (a.Dtype() == DType::Int64) {
            const detail::Int64View view(a);
            for (std::size_t i = 0; i < values.size(); ++i)
                values[i] = static_cast<double>(view.Data()[i]);
        } else {
            const detail::DoubleView view(a);
            std::copy(view.Data(), view.Data() + values.size(), values.begin());
        }
        out[std::stoll(digits)][member] = std::move(values);
    }
    return out;
}

/// The mesh as the file holds it: GRIDs, element cards and properties.
class Nh5wModel {
public:
    explicit Nh5wModel(const Mesh& rMesh) {
        mSignature = nh5w_signature(rMesh);
        mValues = nh5w_property_values(rMesh);
        for (auto& [key, members] : nh5w_pcomp_values(rMesh))
            mValues[key].insert(members.begin(), members.end());
        mPlies = nh5w_plies(rMesh);
        mBeams = nh5w_beam_arrays(rMesh);
        // MAT5, MAT8 and MAT9 values have no table here: their ids are MAT1s.
        for (const std::string& name : rMesh.FieldDataNames())
            for (const char* card : {"MAT5", "MAT8", "MAT9"}) {
                const std::string prefix = std::string("nastran:") + card + "_";
                if (name.rfind(prefix, 0) != 0)
                    continue;
                const std::size_t colon = name.find(':', prefix.size());
                if (colon == std::string::npos)
                    continue;
                const std::string digits = name.substr(prefix.size(), colon - prefix.size());
                if (!digits.empty() && digits.size() <= 15 &&
                    digits.find_first_not_of("0123456789") == std::string::npos)
                    mOtherMaterials.emplace(std::make_pair(std::string(card), std::stoll(digits)));
            }
        mNpoints = rMesh.NumPoints();
        mNodeIds = nh5w_node_ids(rMesh, mNpoints);
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
        const bool have_card = rMesh.HasCellData("nastran:card");
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
            std::vector<std::int64_t> codes;
            if (have_card) {
                const NDArray& a = rMesh.CellData("nastran:card", b);
                if (a.Size() != blk.mCount)
                    nh5w_fail("'nastran:card' must have one value per cell");
                const detail::Int64View view(a);
                codes.assign(view.Data(), view.Data() + blk.mCount);
            }
            const bool line_or_vertex = cb.Type() == "line" || cb.Type() == "vertex";
            std::vector<bool> keep(blk.mCount, true);
            blk.mCards.assign(blk.mCount, kNh5wNumCards);
            for (std::size_t i = 0; i < blk.mCount; ++i) {
                if (have_card && line_or_vertex) {
                    const char* name = nh5w_line_card(codes[i]);
                    if (name == nullptr) {
                        keep[i] = false;
                        mUnwritten.insert(nh5w_card_name(codes[i]));
                    } else if (!nh5w_block_holds(cb.Type(), name)) {
                        nh5w_fail("a '" + cb.Type() + "' block cannot hold " + name + " cards");
                    } else if (nh5w_card_index(name) == kNh5wNumCards) {
                        keep[i] = false;  // a card with no table here
                        mUnwritten.insert(name);
                    } else {
                        blk.mCards[i] = nh5w_card_index(name);
                    }
                } else if (have_card && codes[i] == 31) {
                    // A quad is a CSHEAR; no other cell type holds one.
                    if (cb.Type() == "quad") {
                        blk.mCards[i] = nh5w_card_index("CSHEAR");
                    } else {
                        keep[i] = false;
                        mUnwritten.insert("CSHEAR");
                    }
                } else if (spec->mCard != nullptr) {
                    blk.mCards[i] = nh5w_card_index(spec->mCard);
                } else {
                    const std::int64_t pid = blk.mPids[i];
                    blk.mCards[i] = nh5w_card_index(pid <= 0                ? "CONROD"
                                                    : ptube.count(pid) != 0 ? "CTUBE"
                                                                            : "CROD");
                }
            }
            const std::size_t cbar = nh5w_card_index("CBAR");
            const std::size_t cbeam = nh5w_card_index("CBEAM");
            const std::size_t celas2 = nh5w_card_index("CELAS2");
            const std::size_t cdamp2 = nh5w_card_index("CDAMP2");
            for (std::size_t i = 0; i < blk.mCount; ++i)
                if (keep[i] && blk.mCards[i] == cbar && !rMesh.HasCellData("nastran:x") &&
                    !rMesh.HasCellData("nastran:g0"))
                    nh5w_fail("a CBAR needs its orientation: nastran:x, or nastran:g0 for a G0");
            for (std::size_t i = 0; i < blk.mCount; ++i) {
                if (keep[i] && blk.mCards[i] == celas2 && !rMesh.HasCellData("nastran:k"))
                    nh5w_fail("a CELAS2 needs its value: nastran:k");
                if (keep[i] && blk.mCards[i] == cdamp2 && !rMesh.HasCellData("nastran:b"))
                    nh5w_fail("a CDAMP2 needs its value: nastran:b");
            }
            for (std::size_t i = 0; i < blk.mCount; ++i)
                if (keep[i] && blk.mCards[i] == cbeam && !rMesh.HasCellData("nastran:x") &&
                    !rMesh.HasCellData("nastran:g0"))
                    nh5w_fail("a CBEAM needs its orientation: nastran:x, or nastran:g0 for a G0");
            std::size_t kept = 0;
            for (const bool k : keep)
                kept += k ? 1 : 0;
            if (kept == 0)
                continue;
            // A spring's or damper's components: `nastran:c1` and `nastran:c2`, else 1 at each
            // end. A vertex's second end is ground, so 0.
            std::vector<std::int64_t> c1(blk.mCount, 1);
            std::vector<std::int64_t> c2(blk.mCount, 1);
            const auto read_component = [&](const char* pName, std::vector<std::int64_t>& rValues) {
                if (!rMesh.HasCellData(pName))
                    return;
                const NDArray& a = rMesh.CellData(pName, b);
                if (a.Size() != blk.mCount)
                    nh5w_fail(std::string("'") + pName + "' must have one value per cell");
                const detail::Int64View view(a);
                rValues.assign(view.Data(), view.Data() + blk.mCount);
            };
            read_component("nastran:c1", c1);
            read_component("nastran:c2", c2);
            // The card fields of each row, from the cell data of its own card only.
            blk.mBar.assign(blk.mCount * kNh5wBarWidth, 0.0);
            const auto read_field = [&](const char* pName, std::initializer_list<std::size_t> Cards,
                                        std::size_t Lo, std::size_t Width) {
                if (!rMesh.HasCellData(pName))
                    return;
                const NDArray& a = rMesh.CellData(pName, b);
                if (a.Size() != blk.mCount * Width)
                    nh5w_fail(std::string("'") + pName + "' must have " + std::to_string(Width) +
                              " value(s) per cell");
                const detail::DoubleView view(a);
                for (std::size_t i = 0; i < blk.mCount; ++i)
                    if (std::find(Cards.begin(), Cards.end(), blk.mCards[i]) != Cards.end())
                        for (std::size_t k = 0; k < Width; ++k)
                            blk.mBar[i * kNh5wBarWidth + Lo + k] = view.Data()[i * Width + k];
            };
            read_field("nastran:x", {cbar, cbeam}, 0, 3);
            read_field("nastran:g0", {cbar, cbeam}, 3, 1);
            read_field("nastran:pin", {cbar, cbeam}, 4, 2);
            read_field("nastran:w", {cbar, cbeam}, 6, 6);
            read_field("nastran:k", {celas2}, 0, 1);
            read_field("nastran:ge", {celas2}, 1, 1);
            read_field("nastran:s", {celas2}, 2, 1);
            read_field("nastran:b", {cdamp2}, 0, 1);
            blk.mComps.resize(2 * blk.mCount);
            for (std::size_t i = 0; i < blk.mCount; ++i) {
                blk.mComps[2 * i] = c1[i];
                blk.mComps[2 * i + 1] = blk.mNodes == 1 ? 0 : c2[i];
            }
            if (kept < blk.mCount) {
                const auto compact = [&keep](auto& rValues, std::size_t Stride) {
                    std::size_t to = 0;
                    for (std::size_t i = 0; i < keep.size(); ++i)
                        if (keep[i])
                            for (std::size_t s = 0; s < Stride; ++s)
                                rValues[to++] = rValues[i * Stride + s];
                    rValues.resize(to);
                };
                compact(blk.mConn, blk.mNodes);
                compact(blk.mEids, 1);
                compact(blk.mPids, 1);
                compact(blk.mCards, 1);
                compact(blk.mComps, 2);
                compact(blk.mBar, kNh5wBarWidth);
                blk.mCount = kept;
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
        // A shell whose property has plies is a PCOMP, not a PSHELL.
        for (auto it = mProps[kNh5wPshell].begin(); it != mProps[kNh5wPshell].end();) {
            if (mPlies.count(*it) != 0) {
                mProps[kNh5wPcomp].insert(*it);
                it = mProps[kNh5wPshell].erase(it);
            } else {
                ++it;
            }
        }
        std::map<std::int64_t, int> used;
        for (const int prop : kNh5wUsedOrder)
            for (const std::int64_t p : mProps[prop]) {
                const auto it = used.find(p);
                if (it != used.end())
                    nh5w_fail("property id " + std::to_string(p) + " is used by " +
                              kNh5wPropNames[it->second] + " and " + kNh5wPropNames[prop] +
                              " elements: shells, solids and rods need ids of their own");
                used[p] = prop;
            }
        for (const auto& [pid, prop] : used)
            if (pid <= 0)
                nh5w_fail("property ids must be positive");
        MaterialCards();  // refuses a material id given as two cards
    }

    const Nh5wSignature& Signature() const { return mSignature; }
    std::size_t NumPoints() const { return mNpoints; }
    const std::vector<double>& Points() const { return mPoints; }
    const std::vector<std::int64_t>& NodeIds() const { return mNodeIds; }
    const Nh5wValues& Values() const { return mValues; }
    const Nh5wPlies& Plies() const { return mPlies; }
    const Nh5wArrays& Beams() const { return mBeams; }
    const std::vector<Nh5wBlock>& Blocks() const { return mBlocks; }
    const std::vector<Nh5wRow>& Rows(std::size_t Card) const { return mRows[Card]; }
    const std::set<std::int64_t>& Props(int Prop) const { return mProps[Prop]; }
    const std::set<std::string>& Skipped() const { return mSkipped; }
    const std::set<std::string>& Unwritten() const { return mUnwritten; }

    /// The card of each material id the field data gives (`nastran:MAT1_<id>:*`, `MAT2` or
    /// `MAT4`); an id given as two cards is an error. An id with none is a MAT1.
    std::map<std::int64_t, std::string> MaterialCards() const {
        std::set<std::pair<std::string, std::int64_t>> given = mOtherMaterials;
        for (const auto& [key, members] : mValues)
            if (key.first == "MAT1" || key.first == "MAT2" || key.first == "MAT4")
                given.insert(key);
        std::map<std::int64_t, std::string> cards;
        for (const auto& [card, ident] : given) {
            const auto [it, added] = cards.emplace(ident, card);
            if (!added)
                nh5w_fail("material id " + std::to_string(ident) + " is given as both " +
                          it->second + " and " + card + " field data");
        }
        return cards;
    }

private:
    Nh5wSignature mSignature;
    std::size_t mNpoints = 0;
    std::vector<double> mPoints;
    std::vector<std::int64_t> mNodeIds;
    Nh5wValues mValues;
    Nh5wPlies mPlies;
    Nh5wArrays mBeams;
    std::vector<Nh5wBlock> mBlocks;
    std::vector<std::vector<Nh5wRow>> mRows;
    std::set<std::int64_t> mProps[kNh5wNumProps];
    std::set<std::string> mSkipped;
    std::set<std::string> mUnwritten;  // element cards the file has no table for
    std::set<std::pair<std::string, std::int64_t>> mOtherMaterials;
};

// -- steps ---------------------------------------------------------------------------

using Nh5wMembers = std::vector<std::pair<std::string, std::vector<double>>>;

/// A column of an HDF5 table: int64 or float64 values, or a 4-byte string per row, with
/// `mArity` entries per row.
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

/// One element result table of a step: its rows' element ids and its columns (`GRID` and
/// `PLY` as integers, the members as floats; arity 1 or the array's length).
struct Nh5wTable {
    std::vector<std::int64_t> mEid;
    std::vector<Nh5wColumn> mColumns;
};

using Nh5wKey = std::tuple<std::string, std::size_t, std::string>;  // group, card, variant

struct Nh5wStep {
    std::int64_t mDomain = 1, mSubcase = 1, mStep = 0, mAnalysis = 1, mMode = 1;
    double mTime = 0.0, mEigi = 0.0;
    std::map<std::string, Nh5wMembers> mNodal;  // table -> members
    std::map<Nh5wKey, Nh5wTable> mElemental;
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

double nh5w_time_of(const Mesh& rMesh) {
    return nh5w_scalar_double(rMesh, kSequenceTimeKey, 0.0);
}

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
            // MSC's column order: every real component, then every imaginary one.
            for (std::size_t k = 0; k < 3; ++k)
                nh5w_put(rTables[table], prefix + axes[k] + "R", re[k]);
            for (std::size_t k = 0; k < 3; ++k)
                nh5w_put(rTables[table], prefix + axes[k] + "I", im[k]);
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

/// A cell-data array of the cells of a card, in the card's row order, as its first column;
/// false when no block of the card has it.
bool nh5w_gather(const Mesh& rMesh, const std::string& rName, const Nh5wModel& rModel,
                 std::size_t Card, std::vector<double>& rOut) {
    if (!rMesh.HasCellData(rName))
        return false;
    const std::size_t nblocks = rMesh.CellDataNumBlocks(rName);
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

/// A `(cells, columns)` cell-data array of the cells of a card, in the card's row order, NaN
/// where a block lacks columns: `rWidth` is the widest block's. False when no block of the card
/// has it.
bool nh5w_gather_columns(const Mesh& rMesh, const std::string& rName, const Nh5wModel& rModel,
                         std::size_t Card, std::vector<double>& rOut, std::size_t& rWidth) {
    if (!rMesh.HasCellData(rName))
        return false;
    const std::size_t nblocks = rMesh.CellDataNumBlocks(rName);
    const auto& rows = rModel.Rows(Card);
    std::map<std::size_t, std::pair<std::size_t, std::vector<double>>> by_block;
    rWidth = 0;
    for (const Nh5wRow& row : rows) {
        const Nh5wBlock& blk = rModel.Blocks()[row.mBlock];
        if (by_block.find(row.mBlock) != by_block.end())
            continue;
        if (blk.mIndex >= nblocks)
            nh5w_fail("'" + rName + "' has no array for every cell block");
        const NDArray& a = rMesh.CellData(rName, blk.mIndex);
        const std::size_t n = blk.mCount;
        if (n && (a.Size() < n || a.Size() % n))
            nh5w_fail("'" + rName + "' must have one value per cell and column");
        const detail::DoubleView view(a);
        const std::size_t w = n ? a.Size() / n : 1;
        rWidth = std::max(rWidth, w);
        by_block.emplace(row.mBlock, std::make_pair(w, std::vector<double>(
                                                           view.Data(), view.Data() + a.Size())));
    }
    rWidth = std::max<std::size_t>(rWidth, 1);
    rOut.assign(rows.size() * rWidth, std::numeric_limits<double>::quiet_NaN());
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& [w, v] = by_block[rows[i].mBlock];
        for (std::size_t j = 0; j < w; ++j)
            rOut[i * rWidth + j] = v[rows[i].mRow * w + j];
    }
    return true;
}

/// A `@corner` cell-data array of the cells of a card, in the card's row order, as `Width`
/// entries per row: entry 0 is left NaN for the centre and entry j + 1 holds column j (the
/// cell's j-th node).
bool nh5w_gather_wide(const Mesh& rMesh, const std::string& rName, const Nh5wModel& rModel,
                      std::size_t Card, std::size_t Width, std::vector<double>& rOut) {
    if (!rMesh.HasCellData(rName))
        return false;
    const std::size_t nblocks = rMesh.CellDataNumBlocks(rName);
    const auto& rows = rModel.Rows(Card);
    rOut.assign(rows.size() * Width, std::numeric_limits<double>::quiet_NaN());
    std::map<std::size_t, std::pair<std::size_t, std::vector<double>>> by_block;
    for (const Nh5wRow& row : rows) {
        if (by_block.find(row.mBlock) != by_block.end())
            continue;
        const Nh5wBlock& blk = rModel.Blocks()[row.mBlock];
        if (blk.mIndex >= nblocks)
            nh5w_fail("'" + rName + "' has no array for every cell block");
        const NDArray& a = rMesh.CellData(rName, blk.mIndex);
        const std::size_t n = blk.mCount;
        if (n && (a.Size() < n || a.Size() % n))
            nh5w_fail("'" + rName + "' must have one value per cell and node");
        const detail::DoubleView view(a);
        by_block.emplace(row.mBlock,
                         std::make_pair(n ? a.Size() / n : 1,
                                        std::vector<double>(view.Data(), view.Data() + a.Size())));
    }
    for (std::size_t i = 0; i < rows.size(); ++i) {
        const auto& [w, v] = by_block[rows[i].mBlock];
        for (std::size_t j = 0; j < w && j + 1 < Width; ++j)
            rOut[i * Width + 1 + j] = v[rows[i].mRow * w + j];
    }
    return true;
}

bool nh5w_is_station_card(std::size_t Card) {
    return std::string_view(kNh5wCards[Card].mName) == "CBAR" ||
           std::string_view(kNh5wCards[Card].mName) == "CBEAM";
}

/// The shells a composite's plies are results of, and the family their ply tables take.
const char* nh5w_ply_family(std::size_t Card) {
    const std::string_view name = kNh5wCards[Card].mName;
    if (name == "CQUAD4")
        return "QUAD4";
    if (name == "CTRIA3")
        return "TRIA3";
    if (name == "CQUAD")
        return "QUADR";
    if (name == "CQUAD8")
        return "QUAD8";
    if (name == "CTRIA6")
        return "TRIA6";
    return nullptr;
}

std::vector<std::int64_t> nh5w_eids(const Nh5wModel& rModel, std::size_t Card) {
    std::vector<std::int64_t> out;
    for (const Nh5wRow& row : rModel.Rows(Card))
        out.push_back(row.mEid);
    return out;
}

/// The names of the cell-data arrays `<G>:<M><suffix>` ('' for the plain ones), by group and
/// member, both sorted.
std::map<std::string, std::map<std::string, std::string>> nh5w_split_names(
    const Mesh& rMesh, const std::string& rSuffix) {
    std::map<std::string, std::map<std::string, std::string>> found;
    for (const std::string& name : rMesh.CellDataNames()) {
        const std::size_t colon = name.find(':');
        if (colon == std::string::npos)
            continue;
        std::string base = name;
        if (!rSuffix.empty()) {
            if (!nh5w_ends_with(name, rSuffix))
                continue;
            base = name.substr(0, name.size() - rSuffix.size());
        }
        const std::string group = base.substr(0, colon);
        const std::string member = base.substr(colon + 1);
        if (!nh5w_is_group_name(group) || member.empty() || member.find('@') != std::string::npos)
            continue;
        found[group][member] = name;
    }
    return found;
}

Nh5wColumn* nh5w_find(std::vector<Nh5wColumn>& rColumns, const std::string& rName) {
    for (Nh5wColumn& c : rColumns)
        if (c.mName == rName)
            return &c;
    return nullptr;
}

void nh5w_element_members(const Mesh& rMesh, const Nh5wModel& rModel, std::set<std::string>& rTaken,
                          std::map<Nh5wKey, Nh5wTable>& rOut) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    // The centre values.
    for (const auto& [group, members] : nh5w_split_names(rMesh, "")) {
        for (std::size_t card = 0; card < kNh5wNumCards; ++card) {
            if (kNh5wCards[card].mFamily == nullptr || rModel.Rows(card).empty() ||
                nh5w_is_station_card(card))
                continue;
            Nh5wTable table;
            for (const auto& [member, name] : members) {
                std::vector<double> column;
                if (nh5w_gather(rMesh, name, rModel, card, column)) {
                    rTaken.insert(name);
                    table.mColumns.push_back(nh5w_floats(member, std::move(column)));
                }
            }
            if (!table.mColumns.empty()) {
                table.mEid = nh5w_eids(rModel, card);
                rOut[{group, card, ""}] = std::move(table);
            }
        }
    }
    // `<G>:<M>@corner` (cells, nodes): the corner entries of the array member `M`, whose entry
    // 0 is the centre value `<G>:<M>`.
    for (const std::string& name : rMesh.CellDataNames()) {
        if (!nh5w_ends_with(name, "@corner"))
            continue;
        const std::string base = name.substr(0, name.size() - 7);
        const std::size_t colon = base.find(':');
        if (colon == std::string::npos)
            continue;
        const std::string group = base.substr(0, colon);
        const std::string member = base.substr(colon + 1);
        if (!nh5w_is_group_name(group) || member.empty())
            continue;
        for (std::size_t card = 0; card < kNh5wNumCards; ++card) {
            const auto& rows = rModel.Rows(card);
            if (kNh5wCards[card].mFamily == nullptr || rows.empty() || nh5w_is_station_card(card))
                continue;
            std::size_t width = 1;
            for (const Nh5wRow& row : rows)
                width = std::max(width, rModel.Blocks()[row.mBlock].mNodes + 1);
            std::vector<double> wide;
            if (!nh5w_gather_wide(rMesh, name, rModel, card, width, wide))
                continue;
            Nh5wTable& table = rOut[{group, card, ""}];
            if (table.mEid.empty())
                table.mEid = nh5w_eids(rModel, card);
            if (nh5w_find(table.mColumns, "GRID") == nullptr) {
                std::vector<std::int64_t> grid(rows.size() * width, 0);
                for (std::size_t i = 0; i < rows.size(); ++i) {
                    const Nh5wBlock& blk = rModel.Blocks()[rows[i].mBlock];
                    for (std::size_t j = 0; j < blk.mNodes; ++j)
                        grid[i * width + 1 + j] = rModel.NodeIds()[static_cast<std::size_t>(
                            blk.mConn[rows[i].mRow * blk.mNodes + j])];
                }
                table.mColumns.push_back(nh5w_ints("GRID", std::move(grid), width));
            }
            for (auto m = table.mColumns.begin(); m != table.mColumns.end(); ++m)
                if (m->mName == member && m->mArity == 1 && m->mFloat) {
                    for (std::size_t i = 0; i < rows.size(); ++i)
                        wide[i * width] = m->mF[i];
                    table.mColumns.erase(m);
                    break;
                }
            table.mColumns.push_back(nh5w_floats(member, std::move(wide), width));
            rTaken.insert(name);
        }
    }
    // `<G>:<M>@ply` of the shells as rows of a `<FAMILY>_COMP` table: one per ply with any
    // value, `PLY` (from 1) and the members.
    for (const auto& [group, members] : nh5w_split_names(rMesh, "@ply")) {
        for (std::size_t card = 0; card < kNh5wNumCards; ++card) {
            if (nh5w_ply_family(card) == nullptr || rModel.Rows(card).empty())
                continue;
            std::vector<std::pair<std::string, std::vector<double>>> arrays;
            std::vector<std::size_t> widths;
            for (const auto& [member, name] : members) {
                std::vector<double> wide;
                std::size_t w = 0;
                if (nh5w_gather_columns(rMesh, name, rModel, card, wide, w)) {
                    arrays.emplace_back(member, std::move(wide));
                    widths.push_back(w);
                    rTaken.insert(name);
                }
            }
            if (arrays.empty())
                continue;
            const std::size_t width = *std::max_element(widths.begin(), widths.end());
            const std::size_t n = rModel.Rows(card).size();
            // gather_columns pads to the block's own width: re-pad to the widest
            for (std::size_t a = 0; a < arrays.size(); ++a)
                if (widths[a] != width) {
                    std::vector<double> padded(n * width, nan);
                    for (std::size_t i = 0; i < n; ++i)
                        std::copy(arrays[a].second.begin() + static_cast<long>(i * widths[a]),
                                  arrays[a].second.begin() + static_cast<long>((i + 1) * widths[a]),
                                  padded.begin() + static_cast<long>(i * width));
                    arrays[a].second = std::move(padded);
                }
            const std::vector<std::int64_t> all_eids = nh5w_eids(rModel, card);
            Nh5wTable table;
            std::vector<std::int64_t> ply;
            std::vector<std::pair<std::size_t, std::size_t>> rows;
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t c = 0; c < width; ++c) {
                    bool shown = false;
                    for (const auto& [m, v] : arrays)
                        shown = shown || !std::isnan(v[i * width + c]);
                    if (shown)
                        rows.emplace_back(i, c);
                }
            for (const auto& [i, c] : rows) {
                table.mEid.push_back(all_eids[i]);
                ply.push_back(static_cast<std::int64_t>(c) + 1);
            }
            table.mColumns.push_back(nh5w_ints("PLY", std::move(ply)));
            for (const auto& [m, v] : arrays) {
                std::vector<double> column;
                for (const auto& [i, c] : rows)
                    column.push_back(v[i * width + c]);
                table.mColumns.push_back(nh5w_floats(m, std::move(column)));
            }
            rOut[{group, card, "ply"}] = std::move(table);
        }
    }
    // The bars and beams: a `BEAM` row per beam with 11-slot arrays (`GRID` names its two ends,
    // `SD` the station along it) and a `BARS` row per station of a bar, from
    // `<G>:<M>@station` (entry 0 is end A) and the plain `<G>:<M>` (end A).
    const auto plain = nh5w_split_names(rMesh, "");
    const auto stations = nh5w_split_names(rMesh, "@station");
    std::set<std::string> station_groups;
    for (const auto& [group, members] : plain)
        station_groups.insert(group);
    for (const auto& [group, members] : stations)
        station_groups.insert(group);
    for (const std::string& group : station_groups) {
        for (std::size_t card = 0; card < kNh5wNumCards; ++card) {
            if (!nh5w_is_station_card(card) || rModel.Rows(card).empty())
                continue;
            const auto& rows = rModel.Rows(card);
            const std::size_t n = rows.size();
            const auto plain_group = plain.find(group);
            const auto station_group = stations.find(group);
            std::set<std::string> member_names;
            if (plain_group != plain.end())
                for (const auto& [m, name] : plain_group->second)
                    member_names.insert(m);
            if (station_group != stations.end())
                for (const auto& [m, name] : station_group->second)
                    member_names.insert(m);
            std::vector<std::pair<std::string, std::pair<std::size_t, std::vector<double>>>> arrays;
            for (const std::string& member : member_names) {
                std::vector<double> wide;
                std::size_t w = 0;
                bool have = false;
                if (station_group != stations.end() && station_group->second.count(member)) {
                    const std::string& name = station_group->second.at(member);
                    have = nh5w_gather_columns(rMesh, name, rModel, card, wide, w);
                    rTaken.insert(name);
                }
                if (plain_group != plain.end() && plain_group->second.count(member)) {
                    const std::string& name = plain_group->second.at(member);
                    std::vector<double> single;
                    if (nh5w_gather(rMesh, name, rModel, card, single)) {
                        rTaken.insert(name);
                        if (!have) {
                            wide = std::move(single);
                            w = 1;
                            have = true;
                        } else {
                            for (std::size_t i = 0; i < n; ++i)
                                if (std::isnan(wide[i * w]))
                                    wide[i * w] = single[i];
                        }
                    }
                }
                if (have)
                    arrays.push_back({member, {w, std::move(wide)}});
            }
            if (arrays.empty())
                continue;
            const bool has_sd = std::any_of(arrays.begin(), arrays.end(),
                                            [](const auto& A) { return A.first == "SD"; });
            Nh5wTable table;
            if (std::string_view(kNh5wCards[card].mName) == "CBEAM") {
                table.mEid = nh5w_eids(rModel, card);
                std::vector<bool> shown(n * 11, false);
                std::vector<Nh5wColumn> members;
                for (const auto& [member, wv] : arrays) {
                    const auto& [w, v] = wv;
                    std::vector<double> full(n * 11, nan);
                    for (std::size_t i = 0; i < n; ++i)
                        for (std::size_t k = 0; k < std::min<std::size_t>(11, w); ++k) {
                            full[i * 11 + k] = v[i * w + k];
                            if (!std::isnan(full[i * 11 + k]))
                                shown[i * 11 + k] = true;
                        }
                    members.push_back(nh5w_floats(member, std::move(full), 11));
                }
                std::vector<std::int64_t> grid(n * 11, 0);
                for (std::size_t i = 0; i < n; ++i) {
                    const Nh5wBlock& blk = rModel.Blocks()[rows[i].mBlock];
                    const std::int64_t* conn = blk.mConn.data() + rows[i].mRow * blk.mNodes;
                    grid[i * 11] = rModel.NodeIds()[static_cast<std::size_t>(conn[0])];
                    grid[i * 11 + 10] = rModel.NodeIds()[static_cast<std::size_t>(conn[1])];
                }
                table.mColumns.push_back(nh5w_ints("GRID", std::move(grid), 11));
                if (!has_sd) {  // the station along the beam: 0 to 1 over the 11 slots
                    std::vector<double> sd(n * 11, 0.0);
                    for (std::size_t i = 0; i < n; ++i)
                        for (std::size_t k = 0; k < 11; ++k)
                            sd[i * 11 + k] =
                                shown[i * 11 + k] ? static_cast<double>(k) / 10.0 : 0.0;
                    table.mColumns.push_back(nh5w_floats("SD", std::move(sd), 11));
                }
                for (Nh5wColumn& c : members)
                    table.mColumns.push_back(std::move(c));
            } else {
                std::size_t width = 1;
                for (const auto& [member, wv] : arrays)
                    width = std::max(width, wv.first);
                std::vector<std::pair<std::size_t, std::size_t>> rowcols;
                for (std::size_t i = 0; i < n; ++i)
                    for (std::size_t c = 0; c < width; ++c) {
                        bool shown = false;
                        for (const auto& [member, wv] : arrays)
                            shown =
                                shown || (c < wv.first && !std::isnan(wv.second[i * wv.first + c]));
                        if (shown)
                            rowcols.emplace_back(i, c);
                    }
                const std::vector<std::int64_t> all_eids = nh5w_eids(rModel, card);
                for (const auto& [i, c] : rowcols)
                    table.mEid.push_back(all_eids[i]);
                if (!has_sd) {
                    std::vector<double> sd;
                    for (const auto& [i, c] : rowcols)
                        sd.push_back(static_cast<double>(c) /
                                     static_cast<double>(std::max<std::size_t>(width - 1, 1)));
                    table.mColumns.push_back(nh5w_floats("SD", std::move(sd)));
                }
                for (const auto& [member, wv] : arrays) {
                    std::vector<double> column;
                    for (const auto& [i, c] : rowcols)
                        column.push_back(c < wv.first ? wv.second[i * wv.first + c] : nan);
                    table.mColumns.push_back(nh5w_floats(member, std::move(column)));
                }
            }
            rOut[{group, card, ""}] = std::move(table);
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
    for (const auto& [key, table] : step.mElemental) {
        step.mSignature += "E:" + std::get<0>(key) + "/" + std::to_string(std::get<1>(key)) + "/" +
                           std::get<2>(key);
        for (const Nh5wColumn& c : table.mColumns)
            step.mSignature += "," + c.mName;
        step.mSignature += ";";
    }
    return step;
}

void nh5w_warn_dropped(const Mesh& rMesh, const Nh5wModel& rModel,
                       const std::set<std::string>& rTaken) {
    static const char* const kQuiet[] = {
        "nastran:eid", "nastran:pid", "nastran:card", "nastran:c1", "nastran:c2",
        "nastran:x",   "nastran:g0",  "nastran:pin",  "nastran:w",  "nastran:k",
        "nastran:ge",  "nastran:s",   "nastran:b"};
    std::set<std::string> names;
    for (const std::string& n : rMesh.PointDataNames())
        if (!rTaken.count(n) && n != "nastran:cp" && n != "nastran:cd" && n != "nastran:nid")
            names.insert(n);
    for (const std::string& n : rMesh.CellDataNames())
        if (!rTaken.count(n) &&
            std::find(std::begin(kQuiet), std::end(kQuiet), n) == std::end(kQuiet))
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
    if (!rModel.Unwritten().empty()) {
        std::string joined;
        for (const std::string& n : rModel.Unwritten())
            joined += (joined.empty() ? "" : ", ") + n;
        log::warn("Nastran HDF5: element cards with no table here are left out: {}", joined);
    }
}

// -- HDF5 tables ---------------------------------------------------------------------

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
void nh5w_write_table(hid_t File, const std::string& rPath,
                      const std::vector<Nh5wColumn>& rColumns) {
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
    h5::Hid dset(
        H5Dcreate2(File, rPath.c_str(), file_type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
        H5Dclose);
    if (!dset.Valid())
        nh5w_fail("could not create the table " + rPath);
    if (rows > 0 && H5Dwrite(dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()) < 0)
        nh5w_fail("could not write the table " + rPath);
}

// MSC's columns of a table: the schema the fixtures froze.
struct Nh5wSchemaColumn {
    std::string mName;
    char mKind;  // 'i' int64, 'f' float64, 's' 4-byte string
    std::size_t mArity;
};

/// The columns of the table at @p rPath (below `/NASTRAN`), or null for a table the fixtures
/// do not hold.
const std::vector<Nh5wSchemaColumn>* nh5w_schema(const std::string& rPath) {
    static std::map<std::string, std::vector<Nh5wSchemaColumn>> cache = [] {
        std::map<std::string, std::vector<Nh5wSchemaColumn>> out;
        for (const Nh5wSchemaEntry& entry : kNh5wSchema) {
            std::vector<Nh5wSchemaColumn> columns;
            const std::string text = entry.mColumns;
            std::size_t at = 0;
            while (at < text.size()) {
                std::size_t end = text.find(' ', at);
                if (end == std::string::npos)
                    end = text.size();
                const std::string item = text.substr(at, end - at);
                at = end + 1;
                const std::size_t colon = item.find(':');
                const std::string name = item.substr(0, colon);
                std::string kind = item.substr(colon + 1);
                std::size_t arity = 1;
                const std::size_t bracket = kind.find('[');
                if (bracket != std::string::npos) {
                    arity = static_cast<std::size_t>(std::stoul(kind.substr(bracket + 1)));
                    kind = kind.substr(0, bracket);
                }
                columns.push_back({name, kind[0] == 'S' ? 's' : kind == "f8" ? 'f' : 'i', arity});
            }
            out[entry.mPath] = std::move(columns);
        }
        return out;
    }();
    const auto it = cache.find(rPath);
    return it == cache.end() ? nullptr : &it->second;
}

/// Column @p rSrc as MSC's column of kind @p Kind and @p Arity for @p Rows rows: the first
/// entries of a wider column, zeros after a narrower one.
Nh5wColumn nh5w_fit(const Nh5wColumn& rSrc, const std::string& rName, char Kind, std::size_t Arity,
                    std::size_t Rows) {
    if (Kind == 's') {
        Nh5wColumn c = rSrc;
        c.mName = rName;
        return c;
    }
    Nh5wColumn out;
    out.mName = rName;
    out.mFloat = Kind == 'f';
    out.mArity = Arity;
    const std::size_t wide = Rows == 0 ? 1 : std::max<std::size_t>(rSrc.mArity, 1);
    const auto value = [&](std::size_t R, std::size_t K) -> double {
        const std::size_t at = R * wide + K;
        return rSrc.mFloat ? rSrc.mF[at] : static_cast<double>(rSrc.mI[at]);
    };
    const auto ivalue = [&](std::size_t R, std::size_t K) -> std::int64_t {
        const std::size_t at = R * wide + K;
        return rSrc.mFloat ? static_cast<std::int64_t>(rSrc.mF[at]) : rSrc.mI[at];
    };
    for (std::size_t r = 0; r < Rows; ++r)
        for (std::size_t k = 0; k < Arity; ++k) {
            const bool have = k < wide;
            if (out.mFloat)
                out.mF.push_back(have ? value(r, k) : 0.0);
            else
                out.mI.push_back(have ? ivalue(r, k) : 0);
        }
    return out;
}

/// The columns of the table at @p rPath: MSC's columns in MSC's order, each from @p Provided
/// or its default (zeros, 'CEN/' for a shell's TERM, 'GRID' for a solid's CTYPE, its corner
/// count for NODEF), then any provided column MSC's table lacks, before `DOMAIN_ID`. A table the
/// fixtures do not hold is its provided columns.
std::vector<Nh5wColumn> nh5w_fill(const std::string& rPath, std::vector<Nh5wColumn> Provided,
                                  std::size_t Rows) {
    const std::vector<Nh5wSchemaColumn>* spec = nh5w_schema(rPath);
    if (spec == nullptr)
        return Provided;
    std::size_t grid_arity = 0;
    for (const Nh5wSchemaColumn& c : *spec)
        if (c.mName == "GRID")
            grid_arity = c.mArity;
    std::vector<Nh5wColumn> out;
    std::set<std::string> names;
    for (const Nh5wSchemaColumn& c : *spec) {
        names.insert(c.mName);
        Nh5wColumn* given = nh5w_find(Provided, c.mName);
        if (given != nullptr) {
            out.push_back(nh5w_fit(*given, c.mName, c.mKind, c.mArity, Rows));
            continue;
        }
        if (c.mName == "TERM") {
            out.push_back(nh5w_string4("TERM", "CEN/", Rows));
        } else if (c.mName == "CTYPE") {
            out.push_back(nh5w_string4("CTYPE", "GRID", Rows));
        } else if (c.mKind == 's') {
            out.push_back(nh5w_string4(c.mName, "", Rows));
        } else if (c.mName == "NODEF") {
            out.push_back(nh5w_ints(
                "NODEF", std::vector<std::int64_t>(Rows, grid_arity > 0 ? grid_arity - 1 : 0)));
        } else if (c.mName == "DOMAIN_ID") {
            out.push_back(nh5w_ints("DOMAIN_ID", std::vector<std::int64_t>(Rows, 1)));
        } else if (c.mKind == 'f') {
            out.push_back(
                nh5w_floats(c.mName, std::vector<double>(Rows * c.mArity, 0.0), c.mArity));
        } else {
            out.push_back(
                nh5w_ints(c.mName, std::vector<std::int64_t>(Rows * c.mArity, 0), c.mArity));
        }
    }
    std::vector<Nh5wColumn> extras;
    for (Nh5wColumn& c : Provided)
        if (names.count(c.mName) == 0)
            extras.push_back(std::move(c));
    if (!extras.empty()) {
        Nh5wColumn last;
        const bool has_last = !out.empty() && out.back().mName == "DOMAIN_ID";
        if (has_last) {
            last = std::move(out.back());
            out.pop_back();
        }
        for (Nh5wColumn& c : extras)
            out.push_back(std::move(c));
        if (has_last)
            out.push_back(std::move(last));
    }
    return out;
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

/// Set row `Row` of the columns named in `rValues` (a float or an integer column).
void nh5w_set_values(std::vector<Nh5wColumn>& rColumns, std::size_t Row,
                     const std::map<std::string, double>& rValues) {
    for (Nh5wColumn& c : rColumns) {
        const auto it = rValues.find(c.mName);
        if (it == rValues.end() || c.mString || c.mArity != 1)
            continue;
        if (c.mFloat)
            c.mF[Row] = it->second;
        else
            c.mI[Row] = static_cast<std::int64_t>(it->second);
    }
}

std::vector<Nh5wColumn> nh5w_filled(std::size_t Rows, const std::vector<Nh5wSpecCol>& rSpec) {
    std::vector<Nh5wColumn> out;
    for (const Nh5wSpecCol& c : rSpec) {
        if (c.mKind == 's')
            out.push_back(nh5w_string4(c.mName, c.mText, Rows));
        else if (c.mKind == 'f')
            out.push_back(
                nh5w_floats(c.mName, std::vector<double>(Rows * c.mArity, c.mValue), c.mArity));
        else
            out.push_back(nh5w_ints(
                c.mName,
                std::vector<std::int64_t>(Rows * c.mArity, static_cast<std::int64_t>(c.mValue)),
                c.mArity));
    }
    return out;
}

const std::vector<Nh5wSpecCol>& nh5w_element_extra(std::string_view Card) {
    static const std::vector<Nh5wSpecCol> none;
    static const std::vector<Nh5wSpecCol> quad4 = {{"THETA", 'f', 1, 0.0, ""},
                                                   {"ZOFFS", 'f', 1, 0.0, ""},
                                                   {"TFLAG", 'i', 1, 0, ""},
                                                   {"T", 'f', 4, 0.0, ""},
                                                   {"MCID", 'i', 1, -1, ""}};
    static const std::vector<Nh5wSpecCol> quad8 = {{"T", 'f', 4, 0.0, ""},
                                                   {"THETA", 'f', 1, 0.0, ""},
                                                   {"ZOFFS", 'f', 1, 0.0, ""},
                                                   {"TFLAG", 'i', 1, 0, ""},
                                                   {"MCID", 'i', 1, -1, ""}};
    static const std::vector<Nh5wSpecCol> tria3 = {{"THETA", 'f', 1, 0.0, ""},
                                                   {"ZOFFS", 'f', 1, 0.0, ""},
                                                   {"TFLAG", 'i', 1, 0, ""},
                                                   {"T", 'f', 3, 0.0, ""},
                                                   {"MCID", 'i', 1, -1, ""}};
    static const std::vector<Nh5wSpecCol> tria6 = {{"THETA", 'f', 1, 0.0, ""},
                                                   {"ZOFFS", 'f', 1, 0.0, ""},
                                                   {"T", 'f', 3, 0.0, ""},
                                                   {"TFLAG", 'i', 1, 0, ""},
                                                   {"MCID", 'i', 1, -1, ""}};
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
    static const std::vector<Nh5wSpecCol> prod = {{"MID", 'i', 1, 1, ""},
                                                  {"A", 'f', 1, 1.0, ""},
                                                  {"J", 'f', 1, 0.0, ""},
                                                  {"C", 'f', 1, 0.0, ""},
                                                  {"NSM", 'f', 1, 0.0, ""}};
    static const std::vector<Nh5wSpecCol> pshell = {
        {"MID1", 'i', 1, 1, ""},  {"T", 'f', 1, 1.0, ""},   {"MID2", 'i', 1, -1, ""},
        {"BK", 'f', 1, 1.0, ""},  {"MID3", 'i', 1, -1, ""}, {"TS", 'f', 1, 5.0 / 6.0, ""},
        {"NSM", 'f', 1, 0.0, ""}, {"Z1", 'f', 1, -0.5, ""}, {"Z2", 'f', 1, 0.5, ""},
        {"MID4", 'i', 1, -1, ""}};
    static const std::vector<Nh5wSpecCol> psolid = {
        {"MID", 'i', 1, 1, ""},    {"CORDM", 'i', 1, 0, ""}, {"IN", 'i', 1, 0, ""},
        {"STRESS", 'i', 1, 0, ""}, {"ISOP", 'i', 1, 0, ""},  {"FCTN", 's', 1, 0, "SMEC"}};
    static const std::vector<Nh5wSpecCol> ptube = {{"MID", 'i', 1, 1, ""},
                                                   {"OD", 'f', 1, 1.0, ""},
                                                   {"T", 'f', 1, 0.0, ""},
                                                   {"NSM", 'f', 1, 0.0, ""}};
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

/// The 11 slots of a PBEAM station.
std::size_t nh5w_pbeam_slot(std::size_t I, std::size_t N) {
    return I == 0 ? std::size_t{0} : I == N - 1 ? std::size_t{10} : I;
}

/// A PBEAM's station arrays as `(11, 16)` rows of SO, XXB, A, I1, I2, I12, J, NSM and the outline
/// points C, D, E, F (the columns C1...F2), from the id's field data; an id with none gets one
/// segment of constant section, with A, I1, I2 and J of 1.
struct Nh5wPbeam {
    std::int64_t mMid = 1, mNseg = 1, mCcf = 1;
    std::vector<double> mSlots;  // 11 x 16
};

Nh5wPbeam nh5w_pbeam(std::int64_t Pid, const Nh5wArrays& rArrays) {
    static const std::map<std::string, std::vector<double>> kNone;
    const auto beam = rArrays.find(Pid);
    const auto& given = beam == rArrays.end() ? kNone : beam->second;
    const auto scalar = [&](const char* pName, double Default) {
        const auto it = given.find(pName);
        return it == given.end() || it->second.empty() ? Default : it->second[0];
    };
    Nh5wPbeam out;
    out.mMid = static_cast<std::int64_t>(scalar("MID", 1));
    out.mNseg = static_cast<std::int64_t>(scalar("NSEG", 1));
    out.mCcf = static_cast<std::int64_t>(scalar("CCF", 1));
    if (out.mNseg < 1 || out.mNseg > 10)
        nh5w_fail("PBEAM " + std::to_string(Pid) + " has " + std::to_string(out.mNseg) +
                  " segments; a PBEAM holds 1 to 10");
    const auto n = static_cast<std::size_t>(out.mNseg) + 1;
    static const char* const kStation[] = {"SO", "XXB", "A", "I1", "I2", "I12", "J", "NSM"};
    static const char* const kPoint[] = {"C1", "C2", "D1", "D2"};
    const auto values = [&](const char* pName, std::size_t Count, std::vector<double> Default) {
        const auto it = given.find(pName);
        if (it == given.end())
            return Default;
        if (it->second.size() != Count)
            nh5w_fail("'nastran:PBEAM_" + std::to_string(Pid) + ":" + pName + "' must have " +
                      std::to_string(Count) + " value(s)");
        return it->second;
    };
    out.mSlots.assign(11 * 16, 0.0);
    for (std::size_t m = 0; m < 8; ++m) {
        std::vector<double> def(
            n, std::string_view(kStation[m]) == "NSM" || std::string_view(kStation[m]) == "I12"
                   ? 0.0
                   : 1.0);
        if (std::string_view(kStation[m]) == "XXB") {
            const double step = n > 1 ? 1.0 / static_cast<double>(n - 1) : 0.0;
            for (std::size_t i = 0; i < n; ++i)
                def[i] = static_cast<double>(i) * step;
            def.back() = n > 1 ? 1.0 : 0.0;
        }
        const std::vector<double> v = values(kStation[m], n, def);
        for (std::size_t i = 0; i < n; ++i)
            out.mSlots[nh5w_pbeam_slot(i, n) * 16 + m] = v[i];
    }
    for (std::size_t p = 0; p < 4; ++p) {
        const std::vector<double> v = values(kPoint[p], 2 * n, std::vector<double>(2 * n, 0.0));
        for (std::size_t i = 0; i < n; ++i) {
            out.mSlots[nh5w_pbeam_slot(i, n) * 16 + 8 + 2 * p] = v[2 * i];
            out.mSlots[nh5w_pbeam_slot(i, n) * 16 + 9 + 2 * p] = v[2 * i + 1];
        }
    }
    return out;
}

/// The columns of a card the writer fills from the schema (see `nh5w_fill`).
std::vector<Nh5wColumn> nh5w_bar_card(const Nh5wCard& rCard, const std::vector<std::int64_t>& rEid,
                                      const std::vector<std::int64_t>& rPid,
                                      const std::vector<std::int64_t>& rNodes,
                                      const std::vector<std::int64_t>& rComps,
                                      const std::vector<double>& rBar) {
    const std::size_t n = rEid.size();
    const std::size_t w = rCard.mWidth;
    const auto node = [&](std::size_t K) {
        std::vector<std::int64_t> out(n);
        for (std::size_t i = 0; i < n; ++i)
            out[i] = rNodes[i * w + K];
        return out;
    };
    const auto comp = [&](std::size_t K) {
        std::vector<std::int64_t> out(n);
        for (std::size_t i = 0; i < n; ++i)
            out[i] = rComps[2 * i + K];
        return out;
    };
    const auto bar = [&](std::size_t Lo, std::size_t Width) {
        std::vector<double> out(n * Width);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = 0; k < Width; ++k)
                out[i * Width + k] = rBar[i * kNh5wBarWidth + Lo + k];
        return out;
    };
    const auto bar_int = [&](std::size_t Lo) {
        std::vector<std::int64_t> out(n);
        for (std::size_t i = 0; i < n; ++i)
            out[i] = static_cast<std::int64_t>(rBar[i * kNh5wBarWidth + Lo]);
        return out;
    };
    const std::string_view name = rCard.mName;
    std::vector<Nh5wColumn> out;
    out.push_back(nh5w_ints("EID", rEid));
    if (name == "CBAR" || name == "CBEAM") {
        std::vector<std::int64_t> flag(n);
        const std::vector<std::int64_t> g0 = bar_int(3);
        for (std::size_t i = 0; i < n; ++i)
            flag[i] = g0[i] != 0 ? 2 : 1;
        out.push_back(nh5w_ints("PID", rPid));
        out.push_back(nh5w_ints("GA", node(0)));
        out.push_back(nh5w_ints("GB", node(1)));
        if (name == "CBAR") {
            out.push_back(nh5w_ints("FLAG", flag));
            for (std::size_t k = 0; k < 3; ++k)
                out.push_back(nh5w_floats("X" + std::to_string(k + 1), bar(k, 1)));
            out.push_back(nh5w_ints("GO", g0));
            out.push_back(nh5w_ints("PA", bar_int(4)));
            out.push_back(nh5w_ints("PB", bar_int(5)));
            const char* const w_names[6] = {"W1A", "W2A", "W3A", "W1B", "W2B", "W3B"};
            for (std::size_t k = 0; k < 6; ++k)
                out.push_back(nh5w_floats(w_names[k], bar(6 + k, 1)));
        } else {
            out.push_back(nh5w_floats("X", bar(0, 3), 3));
            out.push_back(nh5w_ints("G0", g0));
            out.push_back(nh5w_ints("F", flag));
            out.push_back(nh5w_ints("PA", bar_int(4)));
            out.push_back(nh5w_ints("PB", bar_int(5)));
            out.push_back(nh5w_floats("WA", bar(6, 3), 3));
            out.push_back(nh5w_floats("WB", bar(9, 3), 3));
        }
        return out;
    }
    if (name == "CELAS2" || name == "CDAMP2") {
        out.push_back(nh5w_floats(name == "CELAS2" ? "K" : "B", bar(0, 1)));
    } else if (name == "CELAS1" || name == "CDAMP1") {
        out.push_back(nh5w_ints("PID", rPid));
    } else if (name == "CSHEAR" || name == "CVISC") {
        out.push_back(nh5w_ints("PID", rPid));
        out.push_back(nh5w_ints("G", rNodes, w));
        return out;
    } else {  // PLOTEL
        out.push_back(nh5w_ints("G", rNodes, w));
        return out;
    }
    out.push_back(nh5w_ints("G1", node(0)));
    out.push_back(nh5w_ints("G2", node(1)));
    out.push_back(nh5w_ints("C1", comp(0)));
    out.push_back(nh5w_ints("C2", comp(1)));
    if (name == "CELAS2") {
        out.push_back(nh5w_floats("GE", bar(1, 1)));
        out.push_back(nh5w_floats("S", bar(2, 1)));
    }
    return out;
}

void nh5w_write_materials(hid_t File, const Nh5wModel& rModel) {
    // The materials the properties name (MID, default 1) and those given.
    std::set<std::int64_t> mids;
    const auto mid_of = [&](const char* pCard, std::int64_t Pid, const char* pMember) {
        const auto given = rModel.Values().find({pCard, Pid});
        if (given != rModel.Values().end()) {
            const auto it = given->second.find(pMember);
            if (it != given->second.end())
                return static_cast<std::int64_t>(it->second);
        }
        return std::int64_t{1};
    };
    for (const std::int64_t pid : rModel.Props(kNh5wPshell))
        mids.insert(mid_of("PSHELL", pid, "MID1"));
    for (const int prop : {kNh5wProd, kNh5wPsolid, kNh5wPtube, kNh5wPbar, kNh5wPshear})
        for (const std::int64_t pid : rModel.Props(prop))
            mids.insert(mid_of(kNh5wPropNames[prop], pid, "MID"));
    for (const std::int64_t pid : rModel.Props(kNh5wPbeam))
        mids.insert(nh5w_pbeam(pid, rModel.Beams()).mMid);
    for (const std::int64_t pid : rModel.Props(kNh5wPcomp))
        for (const auto& ply : rModel.Plies().at(pid))
            mids.insert(static_cast<std::int64_t>(ply[0]));
    const std::map<std::int64_t, std::string> cards = rModel.MaterialCards();
    // any card given an id (the other cards' values are not written: those ids are MAT1s)
    for (const auto& [mid, card] : cards)
        mids.insert(mid);
    {
        std::set<std::int64_t> other;
        for (const auto& [key, members] : rModel.Values())
            if (key.first == "MAT1")
                other.insert(key.second);
    }
    // MAT5, MAT8 and MAT9 ids were folded into `cards` by the model; list them with the rest.
    if (mids.empty())
        return;
    nh5w_require_group(File, "/NASTRAN/INPUT/MATERIAL");
    std::map<std::string, std::vector<std::int64_t>> by_card = {
        {"MAT1", {}}, {"MAT2", {}}, {"MAT4", {}}};
    for (const std::int64_t mid : mids) {
        const auto at = cards.find(mid);
        const std::string card = at == cards.end() ? "MAT1" : at->second;
        by_card[by_card.count(card) ? card : "MAT1"].push_back(mid);
    }
    for (const auto& [card, ids] : by_card) {
        if (ids.empty())
            continue;
        std::vector<Nh5wColumn> columns;
        if (card == "MAT1") {
            columns = nh5w_filled(ids.size(), {{"MID", 'i', 1, 1, ""},
                                               {"E", 'f', 1, 1.0, ""},
                                               {"G", 'f', 1, 0.0, ""},
                                               {"NU", 'f', 1, 0.3, ""},
                                               {"RHO", 'f', 1, 0.0, ""},
                                               {"A", 'f', 1, 0.0, ""},
                                               {"TREF", 'f', 1, 0.0, ""},
                                               {"GE", 'f', 1, 0.0, ""},
                                               {"ST", 'f', 1, 0.0, ""},
                                               {"SC", 'f', 1, 0.0, ""},
                                               {"SS", 'f', 1, 0.0, ""},
                                               {"MCSID", 'i', 1, 0, ""}});
            std::size_t k = 0;
            for (const std::int64_t mid : ids) {
                std::map<std::string, double> values;
                const auto given = rModel.Values().find({"MAT1", mid});
                if (given != rModel.Values().end())
                    values = given->second;
                values["MID"] = static_cast<double>(mid);
                nh5w_set_values(columns, k, values);
                ++k;
            }
            columns.push_back(nh5w_ints("DOMAIN_ID", std::vector<std::int64_t>(ids.size(), 1)));
        } else {
            std::vector<Nh5wColumn> provided = {nh5w_ints("MID", ids)};
            for (const Nh5ValueCard& c : nh5_value_cards()) {
                if (card != c.mCard)
                    continue;
                for (const Nh5Member& m : c.mMembers) {
                    std::vector<double> floats;
                    std::vector<std::int64_t> ints;
                    for (const std::int64_t mid : ids) {
                        double v = m.mDefault;
                        const auto given = rModel.Values().find({card, mid});
                        if (given != rModel.Values().end()) {
                            const auto it = given->second.find(m.mName);
                            if (it != given->second.end())
                                v = it->second;
                        }
                        floats.push_back(v);
                        ints.push_back(static_cast<std::int64_t>(v));
                    }
                    const std::string column = m.mColumn != nullptr ? m.mColumn : m.mName;
                    provided.push_back(m.mInt ? nh5w_ints(column, std::move(ints))
                                              : nh5w_floats(column, std::move(floats)));
                }
            }
            columns = nh5w_fill("INPUT/MATERIAL/" + card, std::move(provided), ids.size());
        }
        nh5w_write_table(File, "/NASTRAN/INPUT/MATERIAL/" + card, columns);
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
    nh5w_write_table(
        File, "/NASTRAN/INPUT/NODE/GRID",
        {nh5w_ints("ID", rModel.NodeIds()), nh5w_ints("CP", std::vector<std::int64_t>(n, 0)),
         nh5w_floats("X", rModel.Points(), 3), nh5w_ints("CD", std::vector<std::int64_t>(n, 0)),
         nh5w_ints("PS", std::vector<std::int64_t>(n, 0)),
         nh5w_ints("SEID", std::vector<std::int64_t>(n, 0)), domain(n)});
    nh5w_require_group(File, "/NASTRAN/INPUT/ELEMENT");
    for (std::size_t c = 0; c < kNh5wNumCards; ++c) {
        const auto& rows = rModel.Rows(c);
        if (rows.empty())
            continue;
        const Nh5wCard& card = kNh5wCards[c];
        std::vector<std::int64_t> eid, pid, nodes(rows.size() * card.mWidth, 0);
        std::vector<std::int64_t> comps;
        std::vector<double> bar;
        for (std::size_t r = 0; r < rows.size(); ++r) {
            const Nh5wBlock& blk = rModel.Blocks()[rows[r].mBlock];
            eid.push_back(rows[r].mEid);
            pid.push_back(blk.mPids[rows[r].mRow]);
            comps.push_back(blk.mComps[2 * rows[r].mRow]);
            comps.push_back(blk.mComps[2 * rows[r].mRow + 1]);
            bar.insert(bar.end(),
                       blk.mBar.begin() + static_cast<long>(rows[r].mRow * kNh5wBarWidth),
                       blk.mBar.begin() + static_cast<long>((rows[r].mRow + 1) * kNh5wBarWidth));
            const std::int64_t* conn = blk.mConn.data() + rows[r].mRow * blk.mNodes;
            for (std::size_t k = 0; k < blk.mNodes; ++k) {
                const std::size_t at =
                    (blk.mNodes > card.mLinear && card.mPerm != nullptr) ? card.mPerm[k] : k;
                if (at < card.mWidth)
                    nodes[r * card.mWidth + at] =
                        rModel.NodeIds()[static_cast<std::size_t>(conn[k])];
            }
        }
        const std::string_view name = card.mName;
        if (name == "CBAR" || name == "CBEAM" || name == "CDAMP1" || name == "CDAMP2" ||
            name == "CELAS1" || name == "CELAS2" || name == "CSHEAR" || name == "CVISC" ||
            name == "PLOTEL") {
            nh5w_write_table(
                File, std::string("/NASTRAN/INPUT/ELEMENT/") + card.mName,
                nh5w_fill(std::string("INPUT/ELEMENT/") + card.mName,
                          nh5w_bar_card(card, eid, pid, nodes, comps, bar), rows.size()));
            continue;
        }
        std::vector<Nh5wColumn> columns;
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
    nh5w_require_group(File, "/NASTRAN/INPUT/PROPERTY");
    for (const int prop : {kNh5wProd, kNh5wPshell, kNh5wPsolid, kNh5wPtube}) {
        const auto& pids = rModel.Props(prop);
        if (pids.empty())
            continue;
        std::vector<Nh5wColumn> columns;
        columns.push_back(nh5w_ints("PID", std::vector<std::int64_t>(pids.begin(), pids.end())));
        for (auto& col : nh5w_filled(pids.size(), nh5w_property_columns(prop)))
            columns.push_back(std::move(col));
        std::size_t k = 0;
        for (const std::int64_t pid : pids) {
            const auto given = rModel.Values().find({kNh5wPropNames[prop], pid});
            if (given != rModel.Values().end())
                nh5w_set_values(columns, k, given->second);
            ++k;
        }
        columns.push_back(domain(pids.size()));
        nh5w_write_table(File, std::string("/NASTRAN/INPUT/PROPERTY/") + kNh5wPropNames[prop],
                         columns);
    }
    for (const int prop : {kNh5wPbar, kNh5wPelas, kNh5wPdamp, kNh5wPvisc, kNh5wPshear}) {
        const auto& pids = rModel.Props(prop);
        if (pids.empty())
            continue;
        std::vector<Nh5wColumn> provided = {
            nh5w_ints("PID", std::vector<std::int64_t>(pids.begin(), pids.end()))};
        for (const Nh5ValueCard& c : nh5_value_cards()) {
            if (std::string_view(c.mCard) != kNh5wPropNames[prop])
                continue;
            for (const Nh5Member& m : c.mMembers) {
                std::vector<double> floats;
                std::vector<std::int64_t> ints;
                for (const std::int64_t pid : pids) {
                    double v = m.mDefault;
                    const auto given = rModel.Values().find({kNh5wPropNames[prop], pid});
                    if (given != rModel.Values().end()) {
                        const auto it = given->second.find(m.mName);
                        if (it != given->second.end())
                            v = it->second;
                    }
                    floats.push_back(v);
                    ints.push_back(static_cast<std::int64_t>(v));
                }
                provided.push_back(m.mInt ? nh5w_ints(m.mName, std::move(ints))
                                          : nh5w_floats(m.mName, std::move(floats)));
            }
        }
        nh5w_write_table(File, std::string("/NASTRAN/INPUT/PROPERTY/") + kNh5wPropNames[prop],
                         nh5w_fill(std::string("INPUT/PROPERTY/") + kNh5wPropNames[prop],
                                   std::move(provided), pids.size()));
    }
    if (const auto& pids = rModel.Props(kNh5wPbeam); !pids.empty()) {
        std::vector<std::int64_t> pid_list(pids.begin(), pids.end()), mid, nsegs, ccf;
        static const char* const kNames[] = {"SO", "XXB", "A",  "I1", "I2", "I12", "J",  "NSM",
                                             "C1", "C2",  "D1", "D2", "E1", "E2",  "F1", "F2"};
        std::vector<std::vector<double>> slots(16);
        for (const std::int64_t pid : pids) {
            const Nh5wPbeam beam = nh5w_pbeam(pid, rModel.Beams());
            mid.push_back(beam.mMid);
            nsegs.push_back(beam.mNseg);
            ccf.push_back(beam.mCcf);
            for (std::size_t m = 0; m < 16; ++m)
                for (std::size_t s = 0; s < 11; ++s)
                    slots[m].push_back(beam.mSlots[s * 16 + m]);
        }
        std::vector<Nh5wColumn> provided = {nh5w_ints("PID", pid_list), nh5w_ints("MID", mid),
                                            nh5w_ints("NSEGS", nsegs), nh5w_ints("CCF", ccf)};
        for (std::size_t m = 0; m < 16; ++m)
            provided.push_back(nh5w_floats(kNames[m], std::move(slots[m]), 11));
        nh5w_write_table(File, "/NASTRAN/INPUT/PROPERTY/PBEAM",
                         nh5w_fill("INPUT/PROPERTY/PBEAM", std::move(provided), pids.size()));
    }
    if (const auto& pids = rModel.Props(kNh5wPcomp); !pids.empty()) {
        nh5w_require_group(File, "/NASTRAN/INPUT/PROPERTY/PCOMP");
        std::vector<std::int64_t> pid_list(pids.begin(), pids.end()), nplies, pos, len, ft;
        std::map<std::string, std::vector<double>> scalars;
        std::vector<double> mid, t, theta, sout;
        std::int64_t position = 0;
        for (const std::int64_t pid : pids) {
            const auto& plies = rModel.Plies().at(pid);
            nplies.push_back(static_cast<std::int64_t>(plies.size()));
            pos.push_back(position);
            len.push_back(static_cast<std::int64_t>(plies.size()));
            position += static_cast<std::int64_t>(plies.size());
            const auto given = rModel.Values().find({"PCOMP", pid});
            for (const char* member : {"Z0", "NSM", "SB", "TREF", "GE"}) {
                double v = 0.0;
                if (given != rModel.Values().end()) {
                    const auto it = given->second.find(member);
                    if (it != given->second.end())
                        v = it->second;
                }
                scalars[member].push_back(v);
            }
            std::int64_t f = 0;
            if (given != rModel.Values().end()) {
                const auto it = given->second.find("FT");
                if (it != given->second.end())
                    f = static_cast<std::int64_t>(it->second);
            }
            ft.push_back(f);
            for (const auto& ply : plies) {
                mid.push_back(ply[0]);
                t.push_back(ply[1]);
                theta.push_back(ply[2]);
                sout.push_back(ply[3]);
            }
        }
        std::vector<Nh5wColumn> head = {nh5w_ints("PID", pid_list),
                                        nh5w_ints("NPLIES", nplies),
                                        nh5w_floats("Z0", scalars["Z0"]),
                                        nh5w_floats("NSM", scalars["NSM"]),
                                        nh5w_floats("SB", scalars["SB"]),
                                        nh5w_ints("FT", ft),
                                        nh5w_floats("TREF", scalars["TREF"]),
                                        nh5w_floats("GE", scalars["GE"]),
                                        nh5w_ints("PLY_POS", pos),
                                        nh5w_ints("PLY_LEN", len)};
        nh5w_write_table(File, "/NASTRAN/INPUT/PROPERTY/PCOMP/IDENTITY",
                         nh5w_fill("INPUT/PROPERTY/PCOMP/IDENTITY", std::move(head), pids.size()));
        std::vector<std::int64_t> mid_i(mid.begin(), mid.end()), sout_i(sout.begin(), sout.end());
        std::vector<Nh5wColumn> ply = {
            nh5w_ints("MID", std::move(mid_i)), nh5w_floats("T", std::move(t)),
            nh5w_floats("THETA", std::move(theta)), nh5w_ints("SOUT", std::move(sout_i))};
        nh5w_write_table(File, "/NASTRAN/INPUT/PROPERTY/PCOMP/PLY",
                         nh5w_fill("INPUT/PROPERTY/PCOMP/PLY", std::move(ply), mid.size()));
    }
    nh5w_write_materials(File, rModel);
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
                if (c.mName == "DOMAIN_ID")
                    continue;  // the table's last column: every row's step
                Nh5wColumn empty;
                empty.mName = c.mName;
                empty.mFloat = c.mFloat;
                empty.mString = c.mString;
                empty.mArity = c.mArity;
                out.mColumns.push_back(std::move(empty));
            }
        std::size_t dst_at = 0;
        for (const Nh5wColumn& src : columns) {
            if (src.mName == "DOMAIN_ID")
                continue;
            auto& dst = out.mColumns[dst_at++];
            dst.mI.insert(dst.mI.end(), src.mI.begin(), src.mI.end());
            dst.mF.insert(dst.mF.end(), src.mF.begin(), src.mF.end());
            dst.mS += src.mS;
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
    nh5w_write_table(
        File, "/INDEX/NASTRAN/RESULT/" + rParent + "/" + rName,
        {nh5w_ints("DOMAIN_ID", out.mIndexDomain), nh5w_ints("POSITION", out.mIndexPosition),
         nh5w_ints("LENGTH", out.mIndexLength)});
}

/// The name of a result table: the card's family (MSC calls a CQUAD4's stresses and strains
/// `QUAD_CN` and its forces `QUAD4_CN`), `_COMP` for plies.
std::string nh5w_table_name(const std::string& rGroup, std::size_t Card,
                            const std::string& rVariant) {
    if (rVariant == "ply")
        return std::string(nh5w_ply_family(Card)) + "_COMP";
    if (std::string_view(kNh5wCards[Card].mName) == "CQUAD4")
        return rGroup == "ELEMENT_FORCE" ? "QUAD4_CN" : "QUAD_CN";
    return kNh5wCards[Card].mFamily;
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
        std::vector<Nh5wColumn> columns = {nh5w_ints("ID", id),
                                           nh5w_ints("SUBCASE", subcase),
                                           nh5w_ints("STEP", step),
                                           nh5w_ints("ANALYSIS", analysis),
                                           nh5w_floats("TIME_FREQ_EIGR", time),
                                           nh5w_floats("EIGI", eigi),
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
                ids.push_back(rModel.NodeIds()[i]);
            columns.push_back(nh5w_ints("ID", ids));
            for (const auto& [m, v] : members) {
                std::vector<double> picked;
                for (const std::size_t i : keep)
                    picked.push_back(v[i]);
                columns.push_back(nh5w_floats(m, std::move(picked)));
            }
            per_step.emplace_back(
                s.mDomain, nh5w_fill("RESULT/NODAL/" + table, std::move(columns), keep.size()));
        }
        nh5w_emit(File, "NODAL", table, per_step);
    }
    std::set<Nh5wKey> keys;
    for (const Nh5wStep& s : rSteps)
        for (const auto& [k, t] : s.mElemental)
            keys.insert(k);
    for (const Nh5wKey& key : keys) {
        const std::string& group = std::get<0>(key);
        const std::string name = nh5w_table_name(group, std::get<1>(key), std::get<2>(key));
        std::vector<std::pair<std::int64_t, std::vector<Nh5wColumn>>> per_step;
        for (const Nh5wStep& s : rSteps) {
            const auto it = s.mElemental.find(key);
            if (it == s.mElemental.end())
                continue;
            const Nh5wTable& table = it->second;
            const std::size_t rows = table.mEid.size();
            std::vector<std::size_t> keep;
            for (std::size_t i = 0; i < rows; ++i) {
                bool all_nan = true;
                for (const Nh5wColumn& c : table.mColumns)
                    if (c.mFloat)
                        for (std::size_t k = 0; k < c.mArity; ++k)
                            all_nan = all_nan && std::isnan(c.mF[i * c.mArity + k]);
                if (!all_nan)
                    keep.push_back(i);
            }
            std::vector<Nh5wColumn> provided;
            std::vector<std::int64_t> eids;
            for (const std::size_t i : keep)
                eids.push_back(table.mEid[i]);
            provided.push_back(nh5w_ints(group == "ENERGY" ? "ID" : "EID", eids));
            for (const Nh5wColumn& c : table.mColumns) {
                Nh5wColumn picked;
                picked.mName = c.mName;
                picked.mFloat = c.mFloat;
                picked.mArity = c.mArity;
                for (const std::size_t i : keep) {
                    if (c.mFloat)
                        picked.mF.insert(picked.mF.end(),
                                         c.mF.begin() + static_cast<long>(i * c.mArity),
                                         c.mF.begin() + static_cast<long>((i + 1) * c.mArity));
                    else
                        picked.mI.insert(picked.mI.end(),
                                         c.mI.begin() + static_cast<long>(i * c.mArity),
                                         c.mI.begin() + static_cast<long>((i + 1) * c.mArity));
                }
                provided.push_back(std::move(picked));
            }
            per_step.emplace_back(s.mDomain, nh5w_fill("RESULT/ELEMENTAL/" + group + "/" + name,
                                                       std::move(provided), keep.size()));
        }
        nh5w_emit(File, "ELEMENTAL/" + group, name, per_step);
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
