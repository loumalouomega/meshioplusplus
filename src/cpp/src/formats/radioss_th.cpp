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
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/radioss_th.hpp"
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/sequence.hpp"

namespace meshioplusplus {

namespace {

// The layout follows OpenRadioss's th_to_csv converter (MIT); no code is copied.

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
    std::string mData;
    std::vector<std::pair<std::size_t, std::size_t>> mRecords;  // (offset, length)
    std::size_t mNext = 0;                                      // next record to read
    std::vector<std::int64_t> mGlobals;
    std::vector<ThEntity> mParts, mSubsets;
    std::vector<ThGroup> mGroups;
    std::vector<double> mFactors;       // mass, length, time (file version > 3050)
    std::vector<std::size_t> mOutputs;  // each output's first record
    std::vector<double> mTimes;

    explicit ThFile(const std::string& rPath) {
        auto in = detail::make_classic_ifstream(rPath, std::ios::binary);
        if (!in)
            th_fail("cannot open '" + rPath + "'");
        mData.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
        if (!is_radioss_th_head(mData.data(), mData.size()))
            th_fail("'" + rPath + "' is not a time-history file");
        // Split the records; a truncated last one ends the file.
        for (std::size_t p = 0; p + 4 <= mData.size();) {
            const std::size_t n = th_be32(mData.data() + p);
            if (p + 8 + n > mData.size() || th_be32(mData.data() + p + 4 + n) != n)
                break;
            mRecords.emplace_back(p + 4, n);
            p += 8 + n;
        }
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
        std::vector<std::int64_t> out(N);
        if (N == 0)
            return out;
        const char* p = Take(4 * N, pWhat);
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
            if (mRecords[mNext].second != 4)
                th_fail("unexpected record size for an output's time");
            for (std::size_t k = 0; k < sizes.size(); ++k)
                if (mRecords[mNext + 1 + k].second != sizes[k])
                    th_fail("unexpected record size in the output at record " +
                            std::to_string(mNext + 1 + k));
            mOutputs.push_back(mNext);
            mTimes.push_back(th_float(mData.data() + mRecords[mNext].first));
            mNext += 1 + sizes.size();
        }
        if (mOutputs.empty())
            th_fail("the file has no complete output");
    }

    const char* Record(std::size_t Index) const { return mData.data() + mRecords[Index].first; }
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

}  // namespace meshioplusplus
