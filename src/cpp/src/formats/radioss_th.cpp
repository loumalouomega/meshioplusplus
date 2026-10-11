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
// System includes
#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/radioss_th.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fortran_records.hpp"
#include "meshioplusplus/detail/provenance.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/operations/sequence.hpp"
#include "../detail/open_source.hpp"
#include "../detail/typed_view.hpp"

namespace meshioplusplus {

namespace {

// The layout follows OpenRadioss's th_to_csv converter (MIT); no code is copied.

// Big-endian Fortran records with 4-byte markers.
constexpr detail::FortranRecordLayout kThLayout{4, true};

[[noreturn]] void th_fail(const std::string& rMessage) {
    throw ReadError("Radioss time history: " + rMessage);
}

std::uint32_t th_be32(const char* p) {
    const auto* u = reinterpret_cast<const unsigned char*>(p);
    return (static_cast<std::uint32_t>(u[0]) << 24) | (static_cast<std::uint32_t>(u[1]) << 16) |
           (static_cast<std::uint32_t>(u[2]) << 8) | static_cast<std::uint32_t>(u[3]);
}

std::int64_t th_int(const char* p) {
    return static_cast<std::int32_t>(th_be32(p));
}

double th_float(const char* p) {
    const std::uint32_t bits = th_be32(p);
    float f = 0.0F;
    std::memcpy(&f, &bits, sizeof f);
    return f;
}

// The name of a global variable by its code (th_to_csv's column titles).
std::string th_global_name(std::int64_t Code) {
    static const char* kNames[] = {"internal_energy",
                                   "kinetic_energy",
                                   "x_momentum",
                                   "y_momentum",
                                   "z_momentum",
                                   "mass",
                                   "time_step",
                                   "rotation_energy",
                                   "external_work",
                                   "spring_energy",
                                   "contact_energy",
                                   "hourglass_energy",
                                   "elastic_contact_energy",
                                   "frictional_contact_energy",
                                   "damping_contact_energy",
                                   "plastic_work",
                                   "added_mass",
                                   "percentage_added_mass",
                                   "inlet_mass",
                                   "outlet_mass",
                                   "inlet_energy",
                                   "outlet_energy"};
    if (Code >= 1 && Code <= static_cast<std::int64_t>(std::size(kNames)))
        return kNames[Code - 1];
    return "var" + std::to_string(Code);
}

// The name of a part or subset variable by its code (Radioss's keywords).
std::string th_part_name(std::int64_t Code) {
    static const char* kNames[] = {
        "IE",    "KE",    "XMOM",   "YMOM", "ZMOM", "MASS", "HE",  "TURBKE", "XCG", "YCG", "ZCG",
        "XXMOM", "YYMOM", "ZZMOM",  "IXX",  "IYY",  "IZZ",  "IXY", "IYZ",    "IZX", "RIE", "KERB",
        "RKERB", "RKE",   "ERODED", "",     "",     "HEAT", "VX",  "VY",     "VZ",  "PW"};
    if (Code >= 1 && Code <= static_cast<std::int64_t>(std::size(kNames)) && kNames[Code - 1][0])
        return kNames[Code - 1];
    return "var" + std::to_string(Code);
}

// A TH group's kind by its type (the /TH/<kind> keyword it came from).
std::string th_group_kind(std::int64_t Type) {
    switch (Type) {
        case 0:
            return "node";
        case 1:
            return "brick";
        case 2:
            return "quad";
        case 3:
            return "shell";
        case 4:
            return "truss";
        case 5:
            return "beam";
        case 6:
            return "spring";
        case 7:
            return "sh3n";
        case 51:
            return "sphcel";
        case 101:
            return "inter";
        case 102:
            return "rwall";
        case 103:
            return "rbody";
        case 104:
            return "sectio";
        case 107:
            return "monvol";
        case 108:
            return "accel";
        default:
            return "type" + std::to_string(Type);
    }
}

struct ThEntity {
    std::int64_t mId = 0;
    std::vector<std::int64_t> mCodes;
};

struct ThGroup {
    std::int64_t mId = 0, mType = 0;
    std::vector<std::int64_t> mIds, mCodes;
};

// One time-history file: its descriptions and the records of every complete
// output.
struct ThFile {
    std::optional<detail::FileSource> mSource;
    std::string_view mData;
    std::vector<detail::FortranRecord> mRecords;  // payload (offset, size)
    std::size_t mNext = 0;                        // next record to read
    std::vector<std::int64_t> mGlobals;
    std::vector<ThEntity> mParts, mSubsets;
    std::vector<ThGroup> mGroups;
    std::vector<double> mFactors;       // mass, length, time (file version > 3050)
    std::vector<std::size_t> mOutputs;  // each output's first record
    std::vector<double> mTimes;

    // mData views mSource's bytes: pinned in place.
    ThFile(const ThFile&) = delete;
    ThFile& operator=(const ThFile&) = delete;
    explicit ThFile(const std::string& rPath) {
        try {
            mSource.emplace(rPath);
        } catch (const ReadError&) {
            th_fail("cannot open '" + rPath + "'");
        }
        mData = mSource->View();
        if (!is_radioss_th_head(mData.data(), mData.size()))
            th_fail("'" + rPath + "' is not a time-history file");
        // Split the records; a truncated last one ends the file.
        mRecords = detail::fortran_records(mData.data(), mData.size(), kThLayout,
                                           "Radioss time history", /*Lenient=*/true);
        Parse();
    }

    // The next record, which must hold Size bytes (any when Size is npos).
    const char* Take(std::size_t Size, const char* pWhat) {
        if (mNext >= mRecords.size())
            th_fail(std::string("the file ends before ") + pWhat);
        const auto [offset, length] = mRecords[mNext];
        if (Size != static_cast<std::size_t>(-1) && length != Size)
            th_fail(std::string("unexpected record size for ") + pWhat);
        ++mNext;
        return mData.data() + offset;
    }

    std::vector<std::int64_t> Ints(std::size_t N, const char* pWhat) {
        if (N == 0)
            return {};
        // Take first: it checks the record holds exactly 4 * N bytes, so a
        // corrupt N never sizes the allocation (4 * N wrapping included).
        if (N > mData.size() / 4)
            th_fail(std::string("implausible size for ") + pWhat);
        const char* p = Take(4 * N, pWhat);
        std::vector<std::int64_t> out(N);
        for (std::size_t i = 0; i < N; ++i)
            out[i] = th_int(p + 4 * i);
        return out;
    }

    static std::size_t Count(std::int64_t N, const char* pWhat) {
        if (N < 0 || N > (1 << 28))
            th_fail(std::string("implausible ") + pWhat + " count");
        return static_cast<std::size_t>(N);
    }

    void Parse() {
        const std::int64_t version = th_int(Take(84, "the title"));
        const std::size_t title = version >= 4021 ? 100 : version >= 3041 ? 80 : 40;
        Take(80, "the date");
        if (version > 3050) {
            Take(4, "the additional records");
            Take(4, "the title length");
            const char* f = Take(12, "the unit factors");
            mFactors = {th_float(f), th_float(f + 4), th_float(f + 8)};
        }
        const std::vector<std::int64_t> counts = Ints(6, "the hierarchy counts");
        const std::size_t nparts = Count(counts[0], "part"), nmats = Count(counts[1], "material"),
                          ngeos = Count(counts[2], "property"), nsubs = Count(counts[3], "subset"),
                          ngroups = Count(counts[4], "TH group"),
                          nglob = Count(counts[5], "global variable");
        mGlobals = Ints(nglob, "the global variables");
        for (std::size_t i = 0; i < nparts; ++i) {
            const char* p = Take(4 + title + 16, "a part");
            ThEntity part;
            part.mId = th_int(p);
            part.mCodes = Ints(Count(th_int(p + 4 + title + 12), "part variable"), "a part");
            mParts.push_back(std::move(part));
        }
        for (std::size_t i = 0; i < nmats + ngeos; ++i)
            Take(4 + title, "a material or property");
        for (std::size_t i = 0; i < nsubs; ++i) {
            const char* p = Take(20 + title, "a subset");
            ThEntity subset;
            subset.mId = th_int(p);
            Ints(Count(th_int(p + 8), "child subset"), "a subset");
            Ints(Count(th_int(p + 12), "subset part"), "a subset");
            subset.mCodes = Ints(Count(th_int(p + 16), "subset variable"), "a subset");
            mSubsets.push_back(std::move(subset));
        }
        for (std::size_t i = 0; i < ngroups; ++i) {
            const char* p = Take(20 + title, "a TH group");
            ThGroup group;
            group.mId = th_int(p);
            group.mType = th_int(p + 4);
            const std::size_t n = Count(th_int(p + 12), "TH group entity");
            for (std::size_t k = 0; k < n; ++k)
                group.mIds.push_back(th_int(Take(4 + title, "a TH group entity")));
            group.mCodes = Ints(Count(th_int(p + 16), "TH group variable"), "a TH group");
            mGroups.push_back(std::move(group));
        }
        // The outputs: time, globals, parts, subsets, one record per group. A
        // run stopped mid-write leaves an incomplete last output, not read.
        std::size_t nparts_vars = 0, nsubs_vars = 0;
        for (const ThEntity& e : mParts)
            nparts_vars += e.mCodes.size();
        for (const ThEntity& e : mSubsets)
            nsubs_vars += e.mCodes.size();
        std::vector<std::size_t> sizes;
        if (nglob > 0)
            sizes.push_back(4 * nglob);
        if (nparts_vars > 0)
            sizes.push_back(4 * nparts_vars);
        if (nsubs_vars > 0)
            sizes.push_back(4 * nsubs_vars);
        for (const ThGroup& g : mGroups)
            sizes.push_back(4 * g.mIds.size() * g.mCodes.size());
        while (mNext + 1 + sizes.size() <= mRecords.size()) {
            if (mRecords[mNext].mSize != 4)
                th_fail("unexpected record size for an output's time");
            for (std::size_t k = 0; k < sizes.size(); ++k)
                if (mRecords[mNext + 1 + k].mSize != sizes[k])
                    th_fail("unexpected record size in the output at record " +
                            std::to_string(mNext + 1 + k));
            mOutputs.push_back(mNext);
            mTimes.push_back(th_float(mData.data() + mRecords[mNext].mOffset));
            mNext += 1 + sizes.size();
        }
        if (mOutputs.empty())
            th_fail("the file has no complete output");
    }

    const char* Record(std::size_t Index) const { return mData.data() + mRecords[Index].mOffset; }
};

void th_scalar(Mesh& rMesh, const ReadOptions& rOpts, const std::string& rName, double Value) {
    if (!rOpts.WantsArray(rName))
        return;
    NDArray a(DType::Float64, {1});
    a.As<double>()[0] = Value;
    rMesh.AddFieldData(rName, std::move(a));
}

void th_ints(Mesh& rMesh, const ReadOptions& rOpts, const std::string& rName,
             const std::vector<std::int64_t>& rValues) {
    if (!rOpts.WantsArray(rName))
        return;
    NDArray a(DType::Int64, {rValues.size()});
    for (std::size_t i = 0; i < rValues.size(); ++i)
        a.As<std::int64_t>()[i] = rValues[i];
    rMesh.AddFieldData(rName, std::move(a));
}

// Part or subset variables: one scalar per (entity, variable), in file order.
void th_entities(Mesh& rMesh, const ReadOptions& rOpts, const std::string& rKind,
                 const std::vector<ThEntity>& rEntities, const char* pValues) {
    std::size_t k = 0;
    for (const ThEntity& e : rEntities)
        for (const std::int64_t code : e.mCodes)
            th_scalar(
                rMesh, rOpts,
                "radioss_th:" + rKind + ":" + std::to_string(e.mId) + ":" + th_part_name(code),
                th_float(pValues + 4 * k++));
}

}  // namespace

bool is_radioss_th_filename(const std::string& rPath) {
    const std::string name = std::filesystem::path(rPath).filename().string();
    const std::size_t n = name.size();
    return n >= 4 && name.find('.') == std::string::npos && name[n - 3] == 'T' &&
           name[n - 2] >= '0' && name[n - 2] <= '9' && name[n - 1] >= '0' && name[n - 1] <= '9';
}

bool is_radioss_th_head(const char* pHead, std::size_t Size) {
    if (Size < 8 || th_be32(pHead) != 84)
        return false;
    const std::int64_t version = th_int(pHead + 4);
    if (version < 1000 || version > 99999)
        return false;
    // The title record's end marker, then the 80-byte date record.
    return Size < 96 || (th_be32(pHead + 88) == 84 && th_be32(pHead + 92) == 80);
}

Mesh read_radioss_th(const std::string& rPath, const ReadOptions& rOpts) {
    const ThFile f(rPath);
    const std::size_t index = rOpts.ResolveTimeStep(f.mOutputs.size());
    std::size_t r = f.mOutputs[index];

    Mesh mesh;
    mesh.AssignPoints(NDArray(DType::Float64, {0, 3}));
    mesh.AddCellBlock("vertex", NDArray(DType::Int64, {0, 1}));
    {
        NDArray t(DType::Float64, {1});
        t.As<double>()[0] = f.mTimes[index];
        mesh.AddFieldData(kSequenceTimeKey, std::move(t));
    }
    if (rOpts.mPointsOnly)
        return mesh;

    ++r;  // the time
    if (!f.mGlobals.empty()) {
        const char* p = f.Record(r++);
        for (std::size_t i = 0; i < f.mGlobals.size(); ++i)
            th_scalar(mesh, rOpts, "radioss_th:global:" + th_global_name(f.mGlobals[i]),
                      th_float(p + 4 * i));
    }
    bool any = false;
    for (const ThEntity& e : f.mParts)
        any = any || !e.mCodes.empty();
    if (any)
        th_entities(mesh, rOpts, "part", f.mParts, f.Record(r++));
    any = false;
    for (const ThEntity& e : f.mSubsets)
        any = any || !e.mCodes.empty();
    if (any)
        th_entities(mesh, rOpts, "subset", f.mSubsets, f.Record(r++));
    for (const ThGroup& g : f.mGroups) {
        const char* p = f.Record(r++);
        const std::string name =
            "radioss_th:" + th_group_kind(g.mType) + ":" + std::to_string(g.mId);
        if (mesh.HasFieldData(name))
            continue;
        if (rOpts.WantsArray(name)) {
            const std::size_t n = g.mIds.size(), m = g.mCodes.size();
            NDArray a(DType::Float64, {n, m});
            for (std::size_t i = 0; i < n * m; ++i)
                a.As<double>()[i] = th_float(p + 4 * i);
            mesh.AddFieldData(name, std::move(a));
        }
        th_ints(mesh, rOpts, name + ":ids", g.mIds);
        th_ints(mesh, rOpts, name + ":variables", g.mCodes);
    }
    if (!f.mFactors.empty() && rOpts.WantsArray("radioss_th:unit_factors")) {
        NDArray a(DType::Float64, {3});
        for (std::size_t i = 0; i < 3; ++i)
            a.As<double>()[i] = f.mFactors[i];
        mesh.AddFieldData("radioss_th:unit_factors", std::move(a));
    }
    return mesh;
}

std::vector<double> radioss_th_time_values(const std::string& rPath) {
    return ThFile(rPath).mTimes;
}

MeshMetadata read_radioss_th_metadata(const std::string& rPath, const ReadOptions& rOpts) {
    ReadOptions options = rOpts;
    options.mPointsOnly = true;
    options.mTimeStep = 0;
    MeshMetadata meta = metadata_from_mesh(read_radioss_th(rPath, options));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "radioss_th";
    meta.mTimeValues = radioss_th_time_values(rPath);
    return meta;
}

// ---------------------------------------------------------------------------
// The writer: the reader's inverse (see doc/formats/radioss_th.md).
// ---------------------------------------------------------------------------

namespace {

// Version 3040 (40-character titles, no unit factors) is what OpenRadioss itself
// writes and `th_to_csv` reads; with unit factors the file is 4021 (100-character
// titles): the reader's two branches that a writer can reach.
constexpr std::int32_t kThwVersionPlain = 3040;
constexpr std::int32_t kThwVersionFactors = 4021;

struct ThwVar {
    std::int64_t mCode = 0;
    std::string mName;
    float mValue = 0.0F;
};

struct ThwEntity {
    std::int64_t mId = 0;
    std::vector<ThwVar> mVars;
};

struct ThwGroup {
    std::int64_t mType = 0, mId = 0;
    std::string mKind;
    std::vector<std::int64_t> mIds, mCodes;
    std::vector<float> mValues;  // entities x variables, row-major
};

// One output's worth of a mesh: everything the file stores.
struct ThwStep {
    float mTime = 0.0F;
    std::optional<std::array<float, 3>> mFactors;
    std::vector<ThwVar> mGlobals;
    std::vector<ThwEntity> mParts, mSubsets;
    std::vector<ThwGroup> mGroups;
};

[[noreturn]] void thw_fail(const std::string& rMessage) {
    throw WriteError("Radioss time history: " + rMessage);
}

// A whole-string decimal integer that fits an int32 (a leading '-' only).
bool thw_parse_int(std::string_view Text, std::int64_t& rOut) {
    if (Text.empty())
        return false;
    for (std::size_t i = 0; i < Text.size(); ++i) {
        const char c = Text[i];
        if (!((c >= '0' && c <= '9') || (i == 0 && c == '-' && Text.size() > 1)))
            return false;
    }
    std::int64_t v = 0;
    const auto [end, ec] = std::from_chars(Text.data(), Text.data() + Text.size(), v);
    if (ec != std::errc() || end != Text.data() + Text.size())
        return false;
    if (v < INT32_MIN || v > INT32_MAX)
        return false;
    rOut = v;
    return true;
}

// The code a variable name stands for: found in the reader's own table (so the
// two directions cannot drift), else `var<code>`.
template <class NameFn>
bool thw_code(std::string_view Name, NameFn Fn, std::int64_t& rCode) {
    for (std::int64_t code = 1; code <= 64; ++code) {
        if (Fn(code) == Name) {
            rCode = code;
            return true;
        }
    }
    return Name.rfind("var", 0) == 0 && thw_parse_int(Name.substr(3), rCode) && rCode >= 1;
}

// The TH group type a kind name stands for (`type<N>` for one the reader has no
// name for).
bool thw_kind(std::string_view Name, std::int64_t& rType) {
    for (std::int64_t type = 0; type <= 200; ++type) {
        if (th_group_kind(type) == Name) {
            rType = type;
            return true;
        }
    }
    return Name.rfind("type", 0) == 0 && thw_parse_int(Name.substr(4), rType) && rType >= 0;
}

std::vector<std::string_view> thw_split(std::string_view Name) {
    std::vector<std::string_view> parts;
    std::size_t start = 0;
    for (;;) {
        const std::size_t colon = Name.find(':', start);
        if (colon == std::string_view::npos) {
            parts.push_back(Name.substr(start));
            return parts;
        }
        parts.push_back(Name.substr(start, colon - start));
        start = colon + 1;
    }
}

bool thw_is_scalar(const NDArray& rArray, float& rValue) {
    if (rArray.Size() != 1)
        return false;
    const detail::DoubleView values(rArray);
    rValue = static_cast<float>(values[0]);
    return true;
}

// An int32 (an id or a variable code) from a value that must be integral.
std::int64_t thw_int_value(double Value, const std::string& rWhat) {
    if (!(Value >= static_cast<double>(INT32_MIN) && Value <= static_cast<double>(INT32_MAX)) ||
        Value != std::floor(Value))
        thw_fail(rWhat + " holds a value that is not an integer in 32-bit range");
    return static_cast<std::int64_t>(Value);
}

std::vector<std::int64_t> thw_int_array(const NDArray& rArray, std::size_t Expected,
                                        const std::string& rName) {
    if (rArray.Size() != Expected)
        thw_fail("'" + rName + "' holds " + std::to_string(rArray.Size()) + " values, expected " +
                 std::to_string(Expected));
    const detail::DoubleView values(rArray);
    std::vector<std::int64_t> out(Expected);
    for (std::size_t i = 0; i < Expected; ++i)
        out[i] = thw_int_value(values[i], "'" + rName + "'");
    return out;
}

bool thw_var_less(const ThwVar& rA, const ThwVar& rB) {
    return rA.mCode != rB.mCode ? rA.mCode < rB.mCode : rA.mName < rB.mName;
}

// The mesh's `radioss_th:*` field data as one output. @p rDropped receives the
// names of the field data that has no place in the file.
ThwStep thw_collect(const Mesh& rMesh, const std::optional<double>& rTime,
                    std::vector<std::string>& rDropped) {
    struct Partial {
        const NDArray* mpData = nullptr;
        const NDArray* mpIds = nullptr;
        const NDArray* mpCodes = nullptr;
        std::vector<std::string> mNames;
    };
    ThwStep step;
    std::map<std::int64_t, std::vector<ThwVar>> parts, subsets;
    std::map<std::tuple<std::int64_t, std::int64_t, std::string>, Partial> groups;
    for (const std::string& name : rMesh.FieldDataNames()) {
        const NDArray& array = rMesh.FieldData(name);
        if (name == kSequenceTimeKey) {
            if (!rTime && array.Size() > 0) {
                const detail::DoubleView values(array);
                step.mTime = static_cast<float>(values[0]);
            }
            continue;
        }
        const std::vector<std::string_view> p = thw_split(name);
        bool taken = false;
        std::int64_t id = 0, code = 0, type = 0;
        float value = 0.0F;
        if (p.size() >= 2 && p[0] == "radioss_th") {
            if (p.size() == 2 && p[1] == "unit_factors") {
                if (array.Size() != 3)
                    thw_fail("'radioss_th:unit_factors' must hold three values");
                const detail::DoubleView values(array);
                step.mFactors = std::array<float, 3>{static_cast<float>(values[0]),
                                                     static_cast<float>(values[1]),
                                                     static_cast<float>(values[2])};
                taken = true;
            } else if (p.size() == 3 && p[1] == "global") {
                if (thw_code(p[2], th_global_name, code) && thw_is_scalar(array, value)) {
                    step.mGlobals.push_back({code, std::string(p[2]), value});
                    taken = true;
                }
            } else if (p.size() == 4 && (p[1] == "part" || p[1] == "subset")) {
                if (thw_parse_int(p[2], id) && thw_code(p[3], th_part_name, code) &&
                    thw_is_scalar(array, value)) {
                    (p[1] == "part" ? parts : subsets)[id].push_back(
                        {code, std::string(p[3]), value});
                    taken = true;
                }
            } else if ((p.size() == 3 || p.size() == 4) && thw_kind(p[1], type) &&
                       thw_parse_int(p[2], id) &&
                       (p.size() == 3 || p[3] == "ids" || p[3] == "variables")) {
                Partial& g = groups[{type, id, std::string(p[1])}];
                if (p.size() == 3)
                    g.mpData = &array;
                else if (p[3] == "ids")
                    g.mpIds = &array;
                else
                    g.mpCodes = &array;
                g.mNames.push_back(name);
                taken = true;
            }
        }
        if (!taken)
            rDropped.push_back(name);
    }
    if (rTime)
        step.mTime = static_cast<float>(*rTime);

    std::sort(step.mGlobals.begin(), step.mGlobals.end(), thw_var_less);
    for (auto* source : {&parts, &subsets}) {
        for (auto& [id, vars] : *source) {
            std::sort(vars.begin(), vars.end(), thw_var_less);
            (source == &parts ? step.mParts : step.mSubsets).push_back({id, std::move(vars)});
        }
    }
    for (const auto& [key, g] : groups) {
        const std::string base =
            "radioss_th:" + std::get<2>(key) + ":" + std::to_string(std::get<1>(key));
        if (!g.mpData) {
            // companions of an array that is not there
            rDropped.insert(rDropped.end(), g.mNames.begin(), g.mNames.end());
            continue;
        }
        const NDArray& data = *g.mpData;
        if (data.Ndim() != 2)
            thw_fail("'" + base + "' must be a 2-D array (entities x variables)");
        const std::size_t n = data.Shape()[0], m = data.Shape()[1];
        ThwGroup group;
        group.mType = std::get<0>(key);
        group.mId = std::get<1>(key);
        group.mKind = std::get<2>(key);
        if (g.mpIds) {
            group.mIds = thw_int_array(*g.mpIds, n, base + ":ids");
        } else {
            group.mIds.resize(n);
            for (std::size_t i = 0; i < n; ++i)
                group.mIds[i] = static_cast<std::int64_t>(i + 1);
        }
        if (g.mpCodes) {
            group.mCodes = thw_int_array(*g.mpCodes, m, base + ":variables");
        } else {
            group.mCodes.resize(m);
            for (std::size_t i = 0; i < m; ++i)
                group.mCodes[i] = static_cast<std::int64_t>(i + 1);
        }
        const detail::DoubleView values(data);
        group.mValues.resize(n * m);
        for (std::size_t i = 0; i < n * m; ++i)
            group.mValues[i] = static_cast<float>(values[i]);
        step.mGroups.push_back(std::move(group));
    }
    std::sort(rDropped.begin(), rDropped.end());
    return step;
}

template <class Entities>
bool thw_same_entities(const Entities& rA, const Entities& rB) {
    if (rA.size() != rB.size())
        return false;
    for (std::size_t i = 0; i < rA.size(); ++i) {
        if (rA[i].mId != rB[i].mId || rA[i].mVars.size() != rB[i].mVars.size())
            return false;
        for (std::size_t k = 0; k < rA[i].mVars.size(); ++k)
            if (rA[i].mVars[k].mCode != rB[i].mVars[k].mCode)
                return false;
    }
    return true;
}

// Whether two outputs have the same variables (the values and time may differ).
bool thw_same_layout(const ThwStep& rA, const ThwStep& rB) {
    if (rA.mFactors.has_value() != rB.mFactors.has_value() ||
        rA.mGlobals.size() != rB.mGlobals.size() || rA.mGroups.size() != rB.mGroups.size())
        return false;
    for (std::size_t i = 0; i < rA.mGlobals.size(); ++i)
        if (rA.mGlobals[i].mCode != rB.mGlobals[i].mCode)
            return false;
    if (!thw_same_entities(rA.mParts, rB.mParts) || !thw_same_entities(rA.mSubsets, rB.mSubsets))
        return false;
    for (std::size_t i = 0; i < rA.mGroups.size(); ++i) {
        const ThwGroup& a = rA.mGroups[i];
        const ThwGroup& b = rB.mGroups[i];
        if (a.mType != b.mType || a.mId != b.mId || a.mKind != b.mKind || a.mIds != b.mIds ||
            a.mCodes != b.mCodes)
            return false;
    }
    return true;
}

// Big-endian Fortran records, built field by field.
class ThwBytes {
public:
    void Int(std::int64_t Value) {
        detail::fortran_put_int32(mPayload, static_cast<std::int32_t>(Value), true);
    }
    void Float(float Value) { detail::fortran_put_float32(mPayload, Value, true); }
    // Text padded with blanks (or cut) to Width characters.
    void Text(const std::string& rText, std::size_t Width) {
        std::string t = rText.substr(0, Width);
        t.resize(Width, ' ');
        mPayload += t;
    }
    // Close the record being built.
    void End() {
        detail::fortran_append_record(mFile, mPayload.data(), mPayload.size(), kThLayout);
        mPayload.clear();
    }
    void Ints(const std::vector<std::int64_t>& rValues) {
        for (const std::int64_t v : rValues)
            Int(v);
        End();
    }
    std::string mFile;

private:
    std::string mPayload;
};

template <class Entities>
std::size_t thw_count_vars(const Entities& rEntities) {
    std::size_t n = 0;
    for (const auto& e : rEntities)
        n += e.mVars.size();
    return n;
}

std::vector<std::int64_t> thw_codes(const std::vector<ThwVar>& rVars) {
    std::vector<std::int64_t> codes;
    for (const ThwVar& v : rVars)
        codes.push_back(v.mCode);
    return codes;
}

// The descriptions: titles, the unit factors, the hierarchy and the variable
// codes of every part, subset and TH group.
std::string thw_header(const ThwStep& rStep) {
    const bool factors = rStep.mFactors.has_value();
    const std::size_t title = factors ? 100 : 40;
    ThwBytes b;
    b.Int(factors ? kThwVersionFactors : kThwVersionPlain);
    b.Text(detail::provenance_lines(detail::SlotTier::Bounded)[0], 80);
    b.End();
    b.Text("meshio++", 80);
    b.End();
    if (factors) {
        b.Int(1);
        b.End();
        b.Int(static_cast<std::int64_t>(title));
        b.End();
        for (const float f : *rStep.mFactors)
            b.Float(f);
        b.End();
    }
    b.Ints({static_cast<std::int64_t>(rStep.mParts.size()), 0, 0,
            static_cast<std::int64_t>(rStep.mSubsets.size()),
            static_cast<std::int64_t>(rStep.mGroups.size()),
            static_cast<std::int64_t>(rStep.mGlobals.size())});
    if (!rStep.mGlobals.empty())
        b.Ints(thw_codes(rStep.mGlobals));
    for (const ThwEntity& part : rStep.mParts) {
        b.Int(part.mId);
        b.Text("PART " + std::to_string(part.mId), title);
        b.Int(0);
        b.Int(0);
        b.Int(0);
        b.Int(static_cast<std::int64_t>(part.mVars.size()));
        b.End();
        if (!part.mVars.empty())
            b.Ints(thw_codes(part.mVars));
    }
    for (const ThwEntity& subset : rStep.mSubsets) {
        b.Int(subset.mId);
        b.Int(0);
        b.Int(0);  // child subsets
        b.Int(0);  // parts
        b.Int(static_cast<std::int64_t>(subset.mVars.size()));
        b.Text("SUBSET " + std::to_string(subset.mId), title);
        b.End();
        if (!subset.mVars.empty())
            b.Ints(thw_codes(subset.mVars));
    }
    for (const ThwGroup& group : rStep.mGroups) {
        b.Int(group.mId);
        b.Int(group.mType);
        b.Int(0);
        b.Int(static_cast<std::int64_t>(group.mIds.size()));
        b.Int(static_cast<std::int64_t>(group.mCodes.size()));
        b.Text("TH " + group.mKind + " " + std::to_string(group.mId), title);
        b.End();
        for (const std::int64_t id : group.mIds) {
            b.Int(id);
            b.Text(group.mKind + " " + std::to_string(id), title);
            b.End();
        }
        if (!group.mCodes.empty())
            b.Ints(group.mCodes);
    }
    return std::move(b.mFile);
}

// One output: the time, the globals, the part and subset values, one record per
// TH group.
std::string thw_output(const ThwStep& rStep) {
    ThwBytes b;
    b.Float(rStep.mTime);
    b.End();
    if (!rStep.mGlobals.empty()) {
        for (const ThwVar& v : rStep.mGlobals)
            b.Float(v.mValue);
        b.End();
    }
    for (const auto* entities : {&rStep.mParts, &rStep.mSubsets}) {
        if (thw_count_vars(*entities) == 0)
            continue;
        for (const ThwEntity& e : *entities)
            for (const ThwVar& v : e.mVars)
                b.Float(v.mValue);
        b.End();
    }
    for (const ThwGroup& group : rStep.mGroups) {
        for (const float v : group.mValues)
            b.Float(v);
        b.End();
    }
    return std::move(b.mFile);
}

void thw_warn_dropped(const std::vector<std::string>& rDropped) {
    if (rDropped.empty())
        return;
    std::string names;
    for (std::size_t i = 0; i < rDropped.size(); ++i)
        names += (i ? ", " : "") + rDropped[i];
    log::warn("Radioss time history: field data without a place in the file is not written: {}",
              names);
}

}  // namespace

void write_radioss_th(const std::string& rPath, const Mesh& rMesh) {
    std::vector<std::string> dropped;
    const ThwStep step = thw_collect(rMesh, std::nullopt, dropped);
    thw_warn_dropped(dropped);
    auto out = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!out)
        throw WriteError("Could not open file for writing: " + rPath);
    const std::string bytes = thw_header(step) + thw_output(step);
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    out.close();
    if (!out)
        throw WriteError("Radioss time history: failed writing " + rPath);
}

struct RadiossThSeriesWriter::Impl {
    std::string mPath;
    decltype(detail::make_classic_ofstream("")) mOut;  // opened by make_classic_ofstream
    std::optional<ThwStep> mFirst;
    std::size_t mSteps = 0;
};

RadiossThSeriesWriter::RadiossThSeriesWriter(const std::string& rPath)
    : mpImpl(std::make_unique<Impl>()) {
    mpImpl->mPath = rPath;
    mpImpl->mOut = detail::make_classic_ofstream(rPath, std::ios::binary);
    if (!mpImpl->mOut)
        throw WriteError("Could not open file for writing: " + rPath);
}

RadiossThSeriesWriter::~RadiossThSeriesWriter() = default;
RadiossThSeriesWriter::RadiossThSeriesWriter(RadiossThSeriesWriter&&) noexcept = default;
RadiossThSeriesWriter& RadiossThSeriesWriter::operator=(RadiossThSeriesWriter&&) noexcept = default;

void RadiossThSeriesWriter::Write(double Time, const Mesh& rMesh) {
    if (!mpImpl)
        thw_fail("series: the writer was moved from");
    Impl& s = *mpImpl;
    std::vector<std::string> dropped;
    const ThwStep step = thw_collect(rMesh, Time, dropped);
    std::string bytes;
    if (!s.mFirst) {
        thw_warn_dropped(dropped);
        bytes = thw_header(step);
        s.mFirst = step;
    } else if (!thw_same_layout(*s.mFirst, step)) {
        thw_fail("series: step " + std::to_string(s.mSteps) +
                 "'s variables differ from the first step's; a time-history file holds one "
                 "set -- write one file per step with '{step}'");
    }
    bytes += thw_output(step);
    s.mOut.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    if (!s.mOut)
        thw_fail("series: failed writing " + s.mPath);
    ++s.mSteps;
}

std::size_t RadiossThSeriesWriter::NumSteps() const noexcept {
    return mpImpl ? mpImpl->mSteps : 0;
}

void RadiossThSeriesWriter::Finalize() {
    if (!mpImpl || !mpImpl->mOut.is_open())
        return;
    if (mpImpl->mSteps == 0)
        thw_fail("series: no step was written to " + mpImpl->mPath);
    mpImpl->mOut.close();
    if (!mpImpl->mOut)
        thw_fail("series: failed writing " + mpImpl->mPath);
}

}  // namespace meshioplusplus
