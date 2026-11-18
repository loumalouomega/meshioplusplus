// Container benchmark (roadmap 5.22.1): the workloads of the core's own maps
// and sets, lifted from the code, each on the standard container and on the
// candidates that could replace it. No external framework: std::chrono,
// a warmup and the median of a few runs, like bench_ops.cpp.
//
//   meshioplusplus_bench_containers [--workload W] [--container C] [--n N]...
//                                   [--runs R]
//
// One CSV row per (workload, container, n): `workload,container,n,median_s,
// checksum,maxrss_kb`. The checksum is the same for every container of a
// workload (a candidate that gives another is wrong, not fast). `maxrss_kb` is
// the process's peak so far: run one (workload, container) per process, as
// tools/bench_containers.sh does, for a number that belongs to it.

// System includes
#include <sys/resource.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes (private, not installed)
#include "../src/detail/edge_use_table.hpp"

namespace {

using std::int64_t;
using std::size_t;
using std::uint64_t;

// ---------------------------------------------------------------- utilities

struct Rng {
    uint64_t mState;
    explicit Rng(uint64_t Seed) : mState(Seed) {}
    uint64_t Next() {
        mState += 0x9E3779B97F4A7C15ull;
        uint64_t z = mState;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    size_t Below(size_t N) { return static_cast<size_t>(Next() % N); }
};

inline uint64_t mix64(uint64_t H) {
    H ^= H >> 33;
    H *= 0xFF51AFD7ED558CCDull;
    H ^= H >> 33;
    H *= 0xC4CEB9FE1A85EC53ull;
    H ^= H >> 33;
    return H;
}

// ------------------------------------------------------ candidate containers

/// Open addressing with linear probing in one array: the weld's table
/// (detail/weld.hpp) made generic. Insert-or-find only; no erase.
template <class K, class V, class Hash>
class FlatMap {
public:
    explicit FlatMap(size_t Expected = 0) { Rehash(Expected < 8 ? 16 : 2 * Expected); }
    /// A reference to the value of `Key`, inserted as `V{}` when absent.
    V& operator[](const K& rKey) {
        if (2 * (mCount + 1) > mSlots.size())
            Rehash(2 * mSlots.size());
        Slot& s = Probe(rKey);
        if (!s.mUsed) {
            s.mUsed = true;
            s.mKey = rKey;
            s.mValue = V{};
            ++mCount;
        }
        return s.mValue;
    }
    const V* Find(const K& rKey) const {
        const Slot& s = Probe(rKey);
        return s.mUsed ? &s.mValue : nullptr;
    }
    size_t size() const { return mCount; }
    /// Call `Fn(key, value)` for every entry (in slot order, not key order).
    template <class Fn>
    void ForEach(Fn&& rFn) const {
        for (const Slot& s : mSlots)
            if (s.mUsed)
                rFn(s.mKey, s.mValue);
    }

private:
    struct Slot {
        K mKey{};
        V mValue{};
        bool mUsed = false;
    };
    size_t mCount = 0;
    std::vector<Slot> mSlots;
    Slot& Probe(const K& rKey) {
        const size_t mask = mSlots.size() - 1;
        size_t i = Hash{}(rKey)&mask;
        while (mSlots[i].mUsed && !(mSlots[i].mKey == rKey))
            i = (i + 1) & mask;
        return mSlots[i];
    }
    const Slot& Probe(const K& rKey) const { return const_cast<FlatMap*>(this)->Probe(rKey); }
    void Rehash(size_t Cap) {
        size_t cap = 16;
        while (cap < Cap)
            cap <<= 1;
        std::vector<Slot> old(cap);
        old.swap(mSlots);
        mCount = 0;
        for (Slot& s : old)
            if (s.mUsed)
                (*this)[s.mKey] = s.mValue;
    }
};

struct IntHash {
    size_t operator()(int64_t K) const {
        return static_cast<size_t>(mix64(static_cast<uint64_t>(K)));
    }
};

/// Id to index for a reader: a vector indexed by `id - min` when the ids are
/// dense (their range within four times their count), the flat map otherwise.
class IdIndexMap {
public:
    explicit IdIndexMap(const std::vector<int64_t>& rIds) : mSparse(rIds.size()) {
        if (rIds.empty())
            return;
        const auto [lo, hi] = std::minmax_element(rIds.begin(), rIds.end());
        const uint64_t range = static_cast<uint64_t>(*hi - *lo) + 1;
        if (range <= 4 * rIds.size()) {
            mDense = true;
            mMin = *lo;
            mTable.assign(range, -1);
            for (size_t i = 0; i < rIds.size(); ++i)
                mTable[static_cast<size_t>(rIds[i] - mMin)] = static_cast<int64_t>(i);
        } else {
            for (size_t i = 0; i < rIds.size(); ++i)
                mSparse[rIds[i]] = static_cast<int64_t>(i);
        }
    }
    int64_t At(int64_t Id) const {
        if (mDense) {
            const int64_t k = Id - mMin;
            return k < 0 || static_cast<size_t>(k) >= mTable.size()
                       ? -1
                       : mTable[static_cast<size_t>(k)];
        }
        const int64_t* v = mSparse.Find(Id);
        return v ? *v : -1;
    }

private:
    bool mDense = false;
    int64_t mMin = 0;
    std::vector<int64_t> mTable;
    FlatMap<int64_t, int64_t, IntHash> mSparse;
};

// -------------------------------------------------------------- the harness

struct Options {
    std::string mWorkload, mContainer;
    std::vector<size_t> mSizes{10000, 160000, 750000};
    int mRuns = 5;
};

using Body = std::function<uint64_t()>;  // returns the checksum

double rss_kb() {
    rusage u{};
    getrusage(RUSAGE_SELF, &u);
    return static_cast<double>(u.ru_maxrss);
}

void run(const Options& rO, const char* pWorkload, const char* pContainer, size_t N,
         const Body& rBody) {
    if (!rO.mWorkload.empty() && rO.mWorkload != pWorkload)
        return;
    if (!rO.mContainer.empty() && rO.mContainer != pContainer)
        return;
    uint64_t sum = rBody();  // warmup
    std::vector<double> t;
    for (int r = 0; r < rO.mRuns; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        const uint64_t s = rBody();
        t.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        if (s != sum) {
            std::fprintf(stderr, "%s/%s: checksum changed between runs\n", pWorkload, pContainer);
            std::exit(2);
        }
    }
    std::sort(t.begin(), t.end());
    std::printf("%s,%s,%zu,%.6f,%016llx,%.0f\n", pWorkload, pContainer, N, t[t.size() / 2],
                static_cast<unsigned long long>(sum), rss_kb());
    std::fflush(stdout);
}

// ---------------------------------------------------------------- workloads

// A reader's id-to-index map: build it from `n` ids in file order, then look
// up the ids of a connectivity, eight per entry (a hexahedron's corners).
struct IdWorkload {
    std::vector<int64_t> mIds;
    std::vector<int64_t> mLookups;
};

IdWorkload make_ids(size_t N, bool Sparse) {
    IdWorkload w;
    w.mIds.resize(N);
    Rng rng(7);
    for (size_t i = 0; i < N; ++i)
        w.mIds[i] = Sparse ? static_cast<int64_t>(7 * i + 3) : static_cast<int64_t>(i + 1);
    if (Sparse)
        for (size_t i = N; i > 1; --i)
            std::swap(w.mIds[i - 1], w.mIds[rng.Below(i)]);
    w.mLookups.resize(8 * N);
    for (int64_t& id : w.mLookups)
        id = w.mIds[rng.Below(N)];
    return w;
}

void id_workloads(const Options& rO, size_t N, bool Sparse) {
    const IdWorkload w = make_ids(N, Sparse);
    const char* name = Sparse ? "idmap_sparse" : "idmap_dense";
    run(rO, name, "std_unordered_map", N, [&] {
        std::unordered_map<int64_t, int64_t> m;
        for (size_t i = 0; i < w.mIds.size(); ++i)
            m[w.mIds[i]] = static_cast<int64_t>(i);
        uint64_t s = 0;
        for (int64_t id : w.mLookups)
            s += static_cast<uint64_t>(m.find(id)->second);
        return s;
    });
    run(rO, name, "std_unordered_map_reserved", N, [&] {
        std::unordered_map<int64_t, int64_t> m;
        m.reserve(w.mIds.size());
        for (size_t i = 0; i < w.mIds.size(); ++i)
            m[w.mIds[i]] = static_cast<int64_t>(i);
        uint64_t s = 0;
        for (int64_t id : w.mLookups)
            s += static_cast<uint64_t>(m.find(id)->second);
        return s;
    });
    run(rO, name, "std_map", N, [&] {
        std::map<int64_t, int64_t> m;
        for (size_t i = 0; i < w.mIds.size(); ++i)
            m[w.mIds[i]] = static_cast<int64_t>(i);
        uint64_t s = 0;
        for (int64_t id : w.mLookups)
            s += static_cast<uint64_t>(m.find(id)->second);
        return s;
    });
    run(rO, name, "flat_map", N, [&] {
        FlatMap<int64_t, int64_t, IntHash> m(w.mIds.size());
        for (size_t i = 0; i < w.mIds.size(); ++i)
            m[w.mIds[i]] = static_cast<int64_t>(i);
        uint64_t s = 0;
        for (int64_t id : w.mLookups)
            s += static_cast<uint64_t>(*m.Find(id));
        return s;
    });
    run(rO, name, "id_index_map", N, [&] {
        const IdIndexMap m(w.mIds);
        uint64_t s = 0;
        for (int64_t id : w.mLookups)
            s += static_cast<uint64_t>(m.At(id));
        return s;
    });
}

// The edges of a triangulated grid with about N points, counted per
// undirected edge (`EdgeUseTable`'s job in the operations).
struct EdgeWorkload {
    std::vector<std::array<int64_t, 3>> mTris;
};

EdgeWorkload make_tris(size_t N) {
    EdgeWorkload w;
    size_t side = 2;
    while (side * side < N)
        ++side;
    for (size_t j = 0; j + 1 < side; ++j)
        for (size_t i = 0; i + 1 < side; ++i) {
            const auto v = [&](size_t a, size_t b) { return static_cast<int64_t>(b * side + a); };
            w.mTris.push_back({v(i, j), v(i + 1, j), v(i + 1, j + 1)});
            w.mTris.push_back({v(i, j), v(i + 1, j + 1), v(i, j + 1)});
        }
    return w;
}

uint64_t edge_summary(uint64_t Edges, uint64_t Squares) {
    return Edges * 1000003ull + Squares;
}

struct PairHash {
    size_t operator()(uint64_t K) const { return static_cast<size_t>(mix64(K)); }
};

void edge_workloads(const Options& rO, size_t N) {
    const EdgeWorkload w = make_tris(N);
    const auto pack = [](int64_t A, int64_t B) {
        const auto lo = static_cast<uint64_t>(std::min(A, B)),
                   hi = static_cast<uint64_t>(std::max(A, B));
        return (lo << 32) | hi;
    };
    run(rO, "edge_count", "std_map_pair", N, [&] {
        std::map<std::pair<int64_t, int64_t>, int> m;
        for (const auto& t : w.mTris)
            for (int k = 0; k < 3; ++k)
                ++m[{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])}];
        uint64_t sq = 0;
        for (const auto& e : m)
            sq += static_cast<uint64_t>(e.second) * static_cast<uint64_t>(e.second);
        return edge_summary(m.size(), sq);
    });
    run(rO, "edge_count", "std_unordered_map_packed", N, [&] {
        std::unordered_map<uint64_t, int> m;
        for (const auto& t : w.mTris)
            for (int k = 0; k < 3; ++k)
                ++m[pack(t[k], t[(k + 1) % 3])];
        uint64_t sq = 0;
        for (const auto& e : m)
            sq += static_cast<uint64_t>(e.second) * static_cast<uint64_t>(e.second);
        return edge_summary(m.size(), sq);
    });
    run(rO, "edge_count", "flat_map_packed", N, [&] {
        FlatMap<uint64_t, int, PairHash> m(2 * w.mTris.size());
        for (const auto& t : w.mTris)
            for (int k = 0; k < 3; ++k)
                ++m[pack(t[k], t[(k + 1) % 3])];
        uint64_t sq = 0;
        m.ForEach(
            [&](uint64_t, int c) { sq += static_cast<uint64_t>(c) * static_cast<uint64_t>(c); });
        return edge_summary(m.size(), sq);
    });
    run(rO, "edge_count", "edge_use_table", N, [&] {
        meshioplusplus::detail::EdgeUseTable<int64_t> t;
        t.Reserve(3 * w.mTris.size());
        for (size_t f = 0; f < w.mTris.size(); ++f)
            for (int k = 0; k < 3; ++k)
                t.Add(w.mTris[f][k], w.mTris[f][(k + 1) % 3], static_cast<int64_t>(f));
        t.Build();
        uint64_t sq = 0;
        for (size_t e = 0; e < t.NumEdges(); ++e)
            sq += static_cast<uint64_t>(t.Count(e)) * static_cast<uint64_t>(t.Count(e));
        return edge_summary(t.NumEdges(), sq);
    });
}

// The weld's cell table: quantised points keyed by a 3-integer grid cell,
// numbered in first-seen order.
struct GridKey {
    int64_t x, y, z;
    bool operator==(const GridKey& r) const { return x == r.x && y == r.y && z == r.z; }
};
struct GridHash {
    size_t operator()(const GridKey& k) const {
        uint64_t h = static_cast<uint64_t>(k.x) * 0x9E3779B97F4A7C15ull;
        h ^= static_cast<uint64_t>(k.y) * 0xC2B2AE3D27D4EB4Full + (h >> 29);
        h ^= static_cast<uint64_t>(k.z) * 0x165667B19E3779F9ull + (h >> 31);
        return static_cast<size_t>(mix64(h));
    }
};

void point_workloads(const Options& rO, size_t N) {
    std::vector<GridKey> keys(N);
    Rng rng(11);
    const size_t side = static_cast<size_t>(std::cbrt(static_cast<double>(N))) + 1;
    for (GridKey& k : keys)
        k = {static_cast<int64_t>(rng.Below(side)), static_cast<int64_t>(rng.Below(side)),
             static_cast<int64_t>(rng.Below(side))};
    run(rO, "point_key", "std_unordered_map", N, [&] {
        std::unordered_map<GridKey, int64_t, GridHash> m;
        uint64_t s = 0;
        for (const GridKey& k : keys) {
            const auto it = m.try_emplace(k, static_cast<int64_t>(m.size())).first;
            s += static_cast<uint64_t>(it->second);
        }
        return s;
    });
    run(rO, "point_key", "flat_map", N, [&] {
        FlatMap<GridKey, int64_t, GridHash> m(N);
        uint64_t s = 0;
        for (const GridKey& k : keys) {
            // First-seen numbering: a new key takes the current size.
            const size_t before = m.size();
            int64_t& v = m[k];
            if (m.size() != before)
                v = static_cast<int64_t>(before);
            s += static_cast<uint64_t>(v);
        }
        return s;
    });
}

// A twelve-entry string table (the cell-type names), looked up N times.
void string_workloads(const Options& rO, size_t N) {
    static const char* const kNames[] = {"vertex", "line",      "triangle", "quad",
                                         "tetra",  "pyramid",   "wedge",    "hexahedron",
                                         "line3",  "triangle6", "tetra10",  "hexahedron20"};
    std::vector<std::string> picks(N);
    Rng rng(5);
    for (std::string& s : picks)
        s = kNames[rng.Below(12)];
    run(rO, "string_key", "std_unordered_map", N, [&] {
        std::unordered_map<std::string, int> m;
        for (int i = 0; i < 12; ++i)
            m[kNames[i]] = i;
        uint64_t s = 0;
        for (const std::string& k : picks)
            s += static_cast<uint64_t>(m.find(k)->second);
        return s;
    });
    run(rO, "string_key", "std_map", N, [&] {
        std::map<std::string, int> m;
        for (int i = 0; i < 12; ++i)
            m[kNames[i]] = i;
        uint64_t s = 0;
        for (const std::string& k : picks)
            s += static_cast<uint64_t>(m.find(k)->second);
        return s;
    });
    run(rO, "string_key", "sorted_array", N, [&] {
        std::vector<std::pair<std::string_view, int>> a;
        for (int i = 0; i < 12; ++i)
            a.push_back({kNames[i], i});
        std::sort(a.begin(), a.end());
        uint64_t s = 0;
        for (const std::string& k : picks) {
            const auto it =
                std::lower_bound(a.begin(), a.end(), std::string_view(k),
                                 [](const auto& p, std::string_view v) { return p.first < v; });
            s += static_cast<uint64_t>(it->second);
        }
        return s;
    });
}

// Build once from unique keys in random order, then walk in key order.
void ordered_workloads(const Options& rO, size_t N) {
    std::vector<int64_t> keys(N);
    for (size_t i = 0; i < N; ++i)
        keys[i] = static_cast<int64_t>(i) * 3;
    Rng rng(3);
    for (size_t i = N; i > 1; --i)
        std::swap(keys[i - 1], keys[rng.Below(i)]);
    run(rO, "ordered_iteration", "std_map", N, [&] {
        std::map<int64_t, int64_t> m;
        for (size_t i = 0; i < keys.size(); ++i)
            m[keys[i]] = static_cast<int64_t>(i);
        uint64_t s = 0, pos = 1;
        for (const auto& e : m)
            s += static_cast<uint64_t>(e.second) * pos++;
        return s;
    });
    run(rO, "ordered_iteration", "sorted_vector", N, [&] {
        std::vector<std::pair<int64_t, int64_t>> v;
        v.reserve(keys.size());
        for (size_t i = 0; i < keys.size(); ++i)
            v.push_back({keys[i], static_cast<int64_t>(i)});
        std::sort(v.begin(), v.end());
        uint64_t s = 0, pos = 1;
        for (const auto& e : v)
            s += static_cast<uint64_t>(e.second) * pos++;
        return s;
    });
}

// A flag per entry: set a random half, then test all.
void flag_workloads(const Options& rO, size_t N) {
    std::vector<size_t> picks(N / 2);
    Rng rng(9);
    for (size_t& p : picks)
        p = rng.Below(N);
    run(rO, "flags", "vector_bool", N, [&] {
        std::vector<bool> f(N, false);
        for (size_t p : picks)
            f[p] = true;
        uint64_t s = 0;
        for (size_t i = 0; i < N; ++i)
            s += f[i] ? i : 0;
        return s;
    });
    run(rO, "flags", "vector_uint8", N, [&] {
        std::vector<std::uint8_t> f(N, 0);
        for (size_t p : picks)
            f[p] = 1;
        uint64_t s = 0;
        for (size_t i = 0; i < N; ++i)
            s += f[i] ? i : 0;
        return s;
    });
}

}  // namespace

int main(int Argc, char** pArgv) {
    Options o;
    std::vector<size_t> sizes;
    for (int i = 1; i < Argc; ++i) {
        const std::string a = pArgv[i];
        if (a == "--workload" && i + 1 < Argc)
            o.mWorkload = pArgv[++i];
        else if (a == "--container" && i + 1 < Argc)
            o.mContainer = pArgv[++i];
        else if (a == "--n" && i + 1 < Argc)
            sizes.push_back(static_cast<size_t>(std::strtoull(pArgv[++i], nullptr, 10)));
        else if (a == "--runs" && i + 1 < Argc)
            o.mRuns = std::atoi(pArgv[++i]);
        else {
            std::fprintf(stderr,
                         "usage: %s [--workload W] [--container C] [--n N]... [--runs R]\n"
                         "workloads: idmap_dense idmap_sparse edge_count point_key string_key "
                         "ordered_iteration flags\n",
                         pArgv[0]);
            return 1;
        }
    }
    if (!sizes.empty())
        o.mSizes = sizes;
    std::printf("workload,container,n,median_s,checksum,maxrss_kb\n");
    for (const size_t n : o.mSizes) {
        id_workloads(o, n, false);
        id_workloads(o, n, true);
        edge_workloads(o, n);
        point_workloads(o, n);
        string_workloads(o, n);
        ordered_workloads(o, n);
        flag_workloads(o, n);
    }
    return 0;
}
