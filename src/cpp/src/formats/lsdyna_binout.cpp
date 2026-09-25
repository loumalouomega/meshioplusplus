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
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <system_error>
#include <array>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/formats/lsdyna_binout.hpp"
#include "meshioplusplus/detail/file_source.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/operations/sequence.hpp"

namespace meshioplusplus {

namespace {

namespace fs = std::filesystem;

// The LSDA layout follows lasso-python's lsda_py3 (BSD-3), itself LSTC's.
constexpr std::int64_t kLbCd = 2, kLbVariable = 4, kLbBegin = 5, kLbOffset = 7, kLbLink = 11;

[[noreturn]] void lb_fail(const std::string& rMessage) {
    throw ReadError("LS-DYNA binout: " + rMessage);
}

// One LSDA file: its field sizes, byte order and bytes.
struct LbFile {
    detail::FileSource mSource;
    std::size_t mLength = 0, mOffset = 0, mCommand = 0, mType = 0, mStart = 0;
    bool mLittle = true;

    explicit LbFile(const std::string& rPath) : mSource(rPath) {
        const char* d = mSource.Data();
        if (!is_binout_head(d, mSource.Size()))
            lb_fail("'" + rPath + "' is not an LSDA (binout) file");
        const auto byte = [&](std::size_t I) { return static_cast<unsigned char>(d[I]); };
        mStart = byte(0);
        mLength = byte(1);
        mOffset = byte(2);
        mCommand = byte(3);
        mType = byte(4);
        mLittle = byte(5) == 1;
    }

    std::uint64_t UInt(std::size_t Pos, std::size_t N) const {
        if (Pos + N > mSource.Size())
            lb_fail("a record runs past the end of the file");
        const auto* p = reinterpret_cast<const unsigned char*>(mSource.Data() + Pos);
        std::uint64_t v = 0;
        for (std::size_t i = 0; i < N; ++i)
            v |= static_cast<std::uint64_t>(p[mLittle ? i : N - 1 - i]) << (8 * i);
        return v;
    }
};

struct LbVar {
    std::int64_t mType = 0;
    std::uint64_t mOffset = 0, mLength = 0;
    const LbFile* mpFile = nullptr;
};

struct LbDir {
    std::map<std::string, std::unique_ptr<LbDir>> mDirs;
    std::map<std::string, LbVar> mVars;
};

// The symbol table of a binout and its continuation files (`%001`...).
class LbArchive {
public:
    explicit LbArchive(const std::string& rPath) {
        std::vector<std::string> paths{rPath};
        std::error_code ec;
        const fs::path path(rPath);
        const std::string stem = path.filename().string() + "%";
        std::vector<std::string> more;
        for (const auto& entry : fs::directory_iterator(
                 path.parent_path().empty() ? fs::path(".") : path.parent_path(), ec)) {
            const std::string name = entry.path().filename().string();
            if (name.size() > stem.size() && name.compare(0, stem.size(), stem) == 0 &&
                std::isdigit(static_cast<unsigned char>(name[stem.size()])))
                more.push_back(entry.path().string());
        }
        std::sort(more.begin(), more.end());
        paths.insert(paths.end(), more.begin(), more.end());
        for (const std::string& p : paths)
            mFiles.push_back(std::make_unique<LbFile>(p));
        for (const auto& f : mFiles)
            ReadSymbols(*f);
    }

    const LbDir* Dir(const std::vector<std::string>& rPath) const {
        const LbDir* d = &mRoot;
        for (const std::string& part : rPath) {
            const auto it = d->mDirs.find(part);
            if (it == d->mDirs.end())
                return nullptr;
            d = it->second.get();
        }
        return d;
    }

    const LbVar* Var(const std::vector<std::string>& rPath) const {
        if (rPath.empty())
            return nullptr;
        const LbDir* d = Dir(std::vector<std::string>(rPath.begin(), rPath.end() - 1));
        if (d == nullptr)
            return nullptr;
        const auto it = d->mVars.find(rPath.back());
        return it == d->mVars.end() ? nullptr : &it->second;
    }

    // A variable's values (a link's target's); rInteger says whether they
    // are integers, rNumeric whether they are numbers (1-byte types are text).
    std::vector<double> Values(const LbVar& rVar, bool& rInteger, bool& rNumeric) const {
        const LbVar* v = &rVar;
        for (int seen = 0; v->mType == kLbLink && seen < 16; ++seen) {
            const std::string target = Raw(*v);
            std::vector<std::string> parts;
            std::size_t a = 0;
            while (a <= target.size()) {
                const std::size_t b = target.find('/', a);
                const std::string part =
                    target.substr(a, b == std::string::npos ? std::string::npos : b - a);
                if (!part.empty())
                    parts.push_back(part);
                if (b == std::string::npos)
                    break;
                a = b + 1;
            }
            v = Var(parts);
            if (v == nullptr)
                lb_fail("a link points to '" + target + "', which is not a variable");
        }
        static const std::size_t kSize[] = {0, 1, 2, 4, 8, 1, 2, 4, 8, 4, 8, 1};
        const std::int64_t t = v->mType;
        if (t < 1 || t > 11)
            lb_fail("a variable has unknown type " + std::to_string(t));
        rInteger = t <= 8;
        rNumeric = kSize[t] > 1;
        const std::string raw = Raw(*v);
        const LbFile& f = *v->mpFile;
        const std::size_t n = kSize[t], count = raw.size() / n;
        std::vector<double> out(count);
        for (std::size_t i = 0; i < count; ++i) {
            std::uint64_t u = 0;
            for (std::size_t b = 0; b < n; ++b)
                u |= static_cast<std::uint64_t>(
                         static_cast<unsigned char>(raw[i * n + (f.mLittle ? b : n - 1 - b)]))
                     << (8 * b);
            switch (t) {
                case 1:
                    out[i] = static_cast<double>(static_cast<std::int8_t>(u));
                    break;
                case 2:
                    out[i] = static_cast<double>(static_cast<std::int16_t>(u));
                    break;
                case 3:
                    out[i] = static_cast<double>(static_cast<std::int32_t>(u));
                    break;
                case 4:
                    out[i] = static_cast<double>(static_cast<std::int64_t>(u));
                    break;
                case 9: {
                    float x;
                    const auto w = static_cast<std::uint32_t>(u);
                    std::memcpy(&x, &w, 4);
                    out[i] = static_cast<double>(x);
                    break;
                }
                case 10: {
                    double x;
                    std::memcpy(&x, &u, 8);
                    out[i] = x;
                    break;
                }
                default:  // unsigned
                    out[i] = static_cast<double>(u);
            }
        }
        return out;
    }

    LbDir mRoot;

private:
    std::string Raw(const LbVar& rVar) const {
        const LbFile& f = *rVar.mpFile;
        static const std::size_t kSize[] = {0, 1, 2, 4, 8, 1, 2, 4, 8, 4, 8, 1};
        const std::size_t head = f.mLength + f.mCommand + f.mType;
        const auto name_len = static_cast<std::size_t>(f.UInt(rVar.mOffset + head, 1));
        const std::size_t pos = rVar.mOffset + head + 1 + name_len;
        const std::size_t size = kSize[std::clamp<std::int64_t>(rVar.mType, 0, 11)] * rVar.mLength;
        if (pos + size > f.mSource.Size())
            lb_fail("a variable's data runs past the end of the file");
        return std::string(f.mSource.Data() + pos, size);
    }

    void ReadSymbols(const LbFile& rF) {
        std::size_t pos = rF.mStart;
        const std::size_t lc = rF.mLength + rF.mCommand;
        if (rF.UInt(pos + rF.mLength, rF.mCommand) != static_cast<std::uint64_t>(kLbOffset))
            return;
        pos += lc;
        LbDir* cwd = &mRoot;
        std::vector<std::string> cwd_path;
        for (;;) {
            const std::uint64_t offset = rF.UInt(pos, rF.mOffset);
            if (offset == 0)
                return;
            pos = offset;
            if (rF.UInt(pos + rF.mLength, rF.mCommand) != static_cast<std::uint64_t>(kLbBegin))
                return;
            pos += lc;
            for (;;) {
                const std::uint64_t length = rF.UInt(pos, rF.mLength);
                const std::uint64_t cmd = rF.UInt(pos + rF.mLength, rF.mCommand);
                if (length < lc)
                    lb_fail("a symbol table record is shorter than its header");
                const std::size_t body = length - lc;
                pos += lc;
                if (cmd == static_cast<std::uint64_t>(kLbCd)) {
                    if (pos + body > rF.mSource.Size())
                        lb_fail("a record runs past the end of the file");
                    Cd(cwd, cwd_path, std::string(rF.mSource.Data() + pos, body));
                    pos += body;
                } else if (cmd == static_cast<std::uint64_t>(kLbVariable)) {
                    const std::size_t tail = rF.mType + rF.mOffset + rF.mLength;
                    if (body < tail || pos + body > rF.mSource.Size())
                        lb_fail("a variable record is malformed");
                    const std::string name(rF.mSource.Data() + pos, body - tail);
                    const std::size_t p = pos + body - tail;
                    LbVar v;
                    v.mType = static_cast<std::int64_t>(rF.UInt(p, rF.mType));
                    v.mOffset = rF.UInt(p + rF.mType, rF.mOffset);
                    v.mLength = rF.UInt(p + rF.mType + rF.mOffset, rF.mLength);
                    v.mpFile = &rF;
                    cwd->mVars[name] = v;
                    pos += body;
                } else {  // the end of this table: the next table's offset follows
                    break;
                }
            }
        }
    }

    void Cd(LbDir*& rpCwd, std::vector<std::string>& rCwdPath, std::string Path) {
        if (!Path.empty() && Path[0] == '/') {
            rpCwd = &mRoot;
            rCwdPath.clear();
            Path.erase(0, 1);
        }
        std::size_t a = 0;
        while (a <= Path.size()) {
            const std::size_t b = Path.find('/', a);
            const std::string part =
                Path.substr(a, b == std::string::npos ? std::string::npos : b - a);
            if (part == "..") {
                if (!rCwdPath.empty())
                    rCwdPath.pop_back();
                rpCwd = &mRoot;
                for (const std::string& p : rCwdPath)
                    rpCwd = rpCwd->mDirs[p].get();
            } else if (!part.empty() && part != ".") {
                auto& child = rpCwd->mDirs[part];
                if (!child)
                    child = std::make_unique<LbDir>();
                rpCwd = child.get();
                rCwdPath.push_back(part);
            }
            if (b == std::string::npos)
                break;
            a = b + 1;
        }
    }

    std::vector<std::unique_ptr<LbFile>> mFiles;
};

// Whether a folder name is an output's: `d` and digits (`d000001`); `deforc`
// is a database.
bool lb_is_step(const std::string& rName) {
    if (rName.size() < 2 || rName[0] != 'd')
        return false;
    return std::all_of(rName.begin() + 1, rName.end(),
                       [](char C) { return std::isdigit(static_cast<unsigned char>(C)) != 0; });
}

std::vector<std::string> lb_steps(const LbArchive& rA, const std::vector<std::string>& rDb) {
    std::vector<std::string> out;
    if (const LbDir* d = rA.Dir(rDb))
        for (const auto& [name, child] : d->mDirs)
            if (lb_is_step(name))
                out.push_back(name);  // std::map: already in order
    return out;
}

// Every database with step folders, as its path: nodout, elout/beam...
void lb_walk(const LbDir& rDir, std::vector<std::string>& rPath,
             std::vector<std::vector<std::string>>& rOut) {
    bool steps = false;
    for (const auto& [name, child] : rDir.mDirs)
        steps = steps || lb_is_step(name);
    if (steps && !rPath.empty())
        rOut.push_back(rPath);
    for (const auto& [name, child] : rDir.mDirs) {
        if (lb_is_step(name) || name == "metadata")
            continue;
        rPath.push_back(name);
        lb_walk(*child, rPath, rOut);
        rPath.pop_back();
    }
}

std::optional<double> lb_time(const LbArchive& rA, std::vector<std::string> Path) {
    Path.push_back("time");
    const LbVar* v = rA.Var(Path);
    if (v == nullptr)
        return std::nullopt;
    bool integer = false, numeric = false;
    const auto values = rA.Values(*v, integer, numeric);
    if (values.empty())
        return std::nullopt;
    return values[0];
}

struct LbBinout {
    LbArchive mArchive;
    std::vector<std::vector<std::string>> mDatabases;
    std::vector<std::string> mMain;
    std::vector<std::string> mSteps;
    std::vector<double> mTimes;
    // the other databases' outputs: (times, folders) in time order
    std::vector<std::pair<std::vector<std::string>,
                          std::pair<std::vector<double>, std::vector<std::string>>>>
        mOthers;

    explicit LbBinout(const std::string& rPath) : mArchive((lb_exists(rPath), rPath)) {
        std::vector<std::string> path;
        lb_walk(mArchive.mRoot, path, mDatabases);
        if (mDatabases.empty())
            lb_fail("the file holds no database with outputs");
        const std::vector<std::string> nodout{"nodout"};
        mMain = std::find(mDatabases.begin(), mDatabases.end(), nodout) != mDatabases.end()
                    ? nodout
                    : mDatabases.front();
        mSteps = lb_steps(mArchive, mMain);
        for (const std::string& s : mSteps) {
            std::vector<std::string> p = mMain;
            p.push_back(s);
            mTimes.push_back(lb_time(mArchive, p).value_or(0.0));
        }
        for (const auto& db : mDatabases) {
            if (db == mMain)
                continue;
            std::vector<std::pair<double, std::string>> pairs;
            for (const std::string& s : lb_steps(mArchive, db)) {
                std::vector<std::string> p = db;
                p.push_back(s);
                if (const auto t = lb_time(mArchive, p))
                    pairs.emplace_back(*t, s);
            }
            std::sort(pairs.begin(), pairs.end());
            std::vector<double> times;
            std::vector<std::string> folders;
            for (const auto& [t, s] : pairs) {
                times.push_back(t);
                folders.push_back(s);
            }
            mOthers.push_back({db, {std::move(times), std::move(folders)}});
        }
    }

    static bool lb_exists(const std::string& rPath) {
        std::error_code ec;
        if (!fs::is_regular_file(rPath, ec))
            lb_fail("'" + rPath + "' does not exist");
        return true;
    }
};

NDArray lb_array(const std::vector<double>& rV, bool Integer) {
    if (Integer) {
        NDArray a(DType::Int64, {rV.size()});
        for (std::size_t i = 0; i < rV.size(); ++i)
            a.As<std::int64_t>()[i] = static_cast<std::int64_t>(rV[i]);
        return a;
    }
    NDArray a(DType::Float64, {rV.size()});
    std::copy(rV.begin(), rV.end(), a.As<double>());
    return a;
}

std::string lb_join(const std::vector<std::string>& rPath) {
    std::string out;
    for (const std::string& p : rPath)
        out += (out.empty() ? "" : "/") + p;
    return out;
}

// The index of the last of Times (sorted) not after Time, their float32
// roundings compared (a binout writes times in single precision); npos if none.
std::size_t lb_latest(const std::vector<double>& rTimes, double Time) {
    const float t = static_cast<float>(Time);
    std::size_t k = 0;
    while (k < rTimes.size() && static_cast<float>(rTimes[k]) <= t)
        ++k;
    return k == 0 ? static_cast<std::size_t>(-1) : k - 1;
}

void lb_field_data(Mesh& rMesh, const LbArchive& rA, const std::vector<std::string>& rDb,
                   const std::string& rStep, const ReadOptions& rOpts, const double* pTime) {
    const std::string prefix = "binout:" + lb_join(rDb) + ":";
    if (pTime != nullptr && rOpts.WantsArray(prefix + "time")) {
        NDArray a(DType::Float64, {1});
        a.As<double>()[0] = *pTime;
        rMesh.AddFieldData(prefix + "time", std::move(a));
    }
    std::vector<std::string> meta = rDb;
    meta.push_back("metadata");
    meta.push_back("ids");
    if (const LbVar* ids = rA.Var(meta)) {
        bool integer = false, numeric = false;
        const auto v = rA.Values(*ids, integer, numeric);
        if (numeric && rOpts.WantsArray(prefix + "ids"))
            rMesh.AddFieldData(prefix + "ids", lb_array(v, true));
    }
    std::vector<std::string> folder = rDb;
    folder.push_back(rStep);
    const LbDir* d = rA.Dir(folder);
    if (d == nullptr)
        return;
    for (const auto& [name, var] : d->mVars) {
        if (name == "time" || !rOpts.WantsArray(prefix + name))
            continue;
        bool integer = false, numeric = false;
        const auto v = rA.Values(var, integer, numeric);
        if (numeric)
            rMesh.AddFieldData(prefix + name, lb_array(v, integer));
    }
}

}  // namespace

bool is_binout_filename(const std::string& rPath) {
    std::string name = fs::path(rPath).filename().string();
    for (char& c : name)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (name.compare(0, 6, "binout") != 0)
        return false;
    return std::all_of(name.begin() + 6, name.end(),
                       [](char C) { return std::isdigit(static_cast<unsigned char>(C)) != 0; });
}

bool is_binout_head(const char* pHead, std::size_t Size) {
    if (Size < 8)
        return false;
    const auto b = [&](std::size_t I) { return static_cast<unsigned char>(pHead[I]); };
    const std::size_t size = b(0), lsize = b(1), osize = b(2);
    if (size < 8 || (lsize != 4 && lsize != 8) || (osize != 4 && osize != 8))
        return false;
    if (b(3) != 1 || b(4) != 1 || b(5) > 1)
        return false;
    if (Size < size + lsize + 1)
        return true;
    return b(size + lsize) == kLbOffset;
}

Mesh read_lsdyna_binout(const std::string& rPath, const ReadOptions& rOpts) {
    const LbBinout b(rPath);
    const LbArchive& a = b.mArchive;
    const std::size_t index = rOpts.ResolveTimeStep(b.mSteps.size());
    const std::string& step = b.mSteps[index];
    const double time = b.mTimes[index];
    std::vector<std::string> folder = b.mMain;
    folder.push_back(step);
    const LbDir* d = a.Dir(folder);
    const bool nodout = b.mMain == std::vector<std::string>{"nodout"};

    Mesh mesh;
    std::vector<double> ids;
    std::size_t n = 0;
    if (nodout) {
        if (const LbVar* v = a.Var({"nodout", "metadata", "ids"})) {
            bool integer = false, numeric = false;
            ids = a.Values(*v, integer, numeric);
        }
        n = ids.size();
    }
    NDArray points(DType::Float64, {n, 3});
    std::fill(points.As<double>(), points.As<double>() + points.Size(),
              std::numeric_limits<double>::quiet_NaN());
    static const char* kCoords[] = {"x_coordinate", "y_coordinate", "z_coordinate"};
    if (nodout && d != nullptr) {
        bool all = n > 0;
        for (const char* c : kCoords)
            all = all && d->mVars.count(c) != 0;
        if (all)
            for (std::size_t c = 0; c < 3; ++c) {
                bool integer = false, numeric = false;
                const auto v = a.Values(d->mVars.at(kCoords[c]), integer, numeric);
                for (std::size_t i = 0; i < n && i < v.size(); ++i)
                    points.As<double>()[3 * i + c] = v[i];
            }
    }
    mesh.AssignPoints(std::move(points));
    NDArray conn(DType::Int64, {n, 1});
    for (std::size_t i = 0; i < n; ++i)
        conn.As<std::int64_t>()[i] = static_cast<std::int64_t>(i);
    mesh.AddCellBlock("vertex", std::move(conn));
    {
        NDArray t(DType::Float64, {1});
        t.As<double>()[0] = time;
        mesh.AddFieldData(kSequenceTimeKey, std::move(t));
    }
    if (nodout && !ids.empty())
        mesh.AddPointData("lsdyna:nid", lb_array(ids, true));
    if (rOpts.mPointsOnly)
        return mesh;

    if (nodout && d != nullptr) {
        // nodout's per-node variables: (point data name, its x, y, z variables)
        static const std::pair<const char*, std::array<const char*, 3>> kVectors[] = {
            {"displacement", {"x_displacement", "y_displacement", "z_displacement"}},
            {"rotation", {"rx_displacement", "ry_displacement", "rz_displacement"}},
            {"velocity", {"x_velocity", "y_velocity", "z_velocity"}},
            {"rotational_velocity", {"rx_velocity", "ry_velocity", "rz_velocity"}},
            {"acceleration", {"x_acceleration", "y_acceleration", "z_acceleration"}},
            {"rotational_acceleration", {"rx_acceleration", "ry_acceleration", "rz_acceleration"}},
        };
        std::vector<std::string> used{"x_coordinate", "y_coordinate", "z_coordinate", "time"};
        for (const auto& [name, parts] : kVectors) {
            bool all = true;
            for (const char* p : parts)
                all = all && d->mVars.count(p) != 0;
            if (!all)
                continue;
            used.insert(used.end(), parts.begin(), parts.end());
            if (!rOpts.WantsArray(name))
                continue;
            NDArray arr(DType::Float64, {n, 3});
            std::fill(arr.As<double>(), arr.As<double>() + arr.Size(),
                      std::numeric_limits<double>::quiet_NaN());
            for (std::size_t c = 0; c < 3; ++c) {
                bool integer = false, numeric = false;
                const auto v = a.Values(d->mVars.at(parts[c]), integer, numeric);
                for (std::size_t i = 0; i < n && i < v.size(); ++i)
                    arr.As<double>()[3 * i + c] = v[i];
            }
            mesh.AddPointData(name, std::move(arr));
        }
        for (const auto& [name, var] : d->mVars) {
            if (std::find(used.begin(), used.end(), name) != used.end())
                continue;
            bool integer = false, numeric = false;
            const auto v = a.Values(var, integer, numeric);
            if (!numeric)
                continue;
            if (v.size() == n && n > 1) {
                if (rOpts.WantsArray(name))
                    mesh.AddPointData(name, lb_array(v, false));
            } else if (rOpts.WantsArray("binout:nodout:" + name)) {
                mesh.AddFieldData("binout:nodout:" + name, lb_array(v, integer));
            }
        }
    } else {
        lb_field_data(mesh, a, b.mMain, step, rOpts, nullptr);
    }
    // every other database's latest output at or before this time
    for (const auto& [db, entry] : b.mOthers) {
        const auto& [times, folders] = entry;
        const std::size_t k = lb_latest(times, time);
        if (k != static_cast<std::size_t>(-1))
            lb_field_data(mesh, a, db, folders[k], rOpts, &times[k]);
    }
    return mesh;
}

std::vector<double> lsdyna_binout_time_values(const std::string& rPath) {
    return LbBinout(rPath).mTimes;
}

MeshMetadata read_lsdyna_binout_metadata(const std::string& rPath, const ReadOptions& rOpts) {
    ReadOptions options = rOpts;
    options.mPointsOnly = true;
    options.mTimeStep = 0;
    MeshMetadata meta = metadata_from_mesh(read_lsdyna_binout(rPath, options));
    meta.mFellBackToFullRead = true;
    meta.mFormat = "lsdyna_binout";
    meta.mTimeValues = lsdyna_binout_time_values(rPath);
    return meta;
}

}  // namespace meshioplusplus
