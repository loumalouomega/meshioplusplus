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
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/nastran_h5.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/hdf5_util.hpp"
#include "meshioplusplus/detail/nastran_model.hpp"
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

std::string nh5_join(const std::vector<std::string>& rV) {
    std::string out;
    for (const std::string& s : rV)
        out += (out.empty() ? "" : ", ") + s;
    return out;
}

}  // namespace

Mesh read_nastran_h5(const std::string& rPath, const ReadOptions& rOpts) {
    const Nh5File file(rPath);
    const hid_t f = file.mFile;
    Mesh mesh;

    // --- points ---------------------------------------------------------------
    const std::size_t npts = file.mGridIds.size();
    {
        const DType as = DType::Float64;
        NDArray x = h5::read_compound_member(f, kNh5Grid, "X", 0, npts, &as);
        if (x.Shape().size() != 2 || x.Shape()[1] != 3)
            nh5_fail("GRID X is not a 3-vector");
        mesh.AssignPoints(std::move(x));
        const auto members = h5::compound_members(f, kNh5Grid);
        for (const char* frame : {"CP", "CD"})
            if (nh5_has_member(members, frame))
                detail::nastran_add_frame(mesh, frame, nh5_int_member(f, kNh5Grid, frame, 0, npts),
                                          "MSC Nastran HDF5");
    }

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
    for (const Nh5ResultTable& t : file.mTables) {
        const auto it = t.mIndex.find(dom.mId);
        if (it == t.mIndex.end())
            continue;
        const auto [row0, count] = it->second;

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
        constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
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
    return mesh;
}

MeshMetadata read_nastran_h5_metadata(const std::string& rPath, const ReadOptions& /*rOpts*/) {
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

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
