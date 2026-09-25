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

/**
 * @file bench_ops.cpp
 * @brief Timings of the mesh operations over a size sweep (roadmap §3/§4).
 *
 * The companion of `bench_backends.cpp` for operations rather than I/O:
 * `extract_surface`, `smooth`, `refine`, `merge`, `clean`, `compute_sdf`,
 * `decimate`, `partition` and `reorder` on the same structured tetrahedral
 * cube, at several sizes. The parallel backend is a compile-time choice, so a
 * binary measures one backend; `tools/bench_ops.sh` builds SEQ, OpenMP and TBB
 * trees and sweeps `OMP_NUM_THREADS`. Warmup plus the median of N runs,
 * `std::chrono::steady_clock`, no framework.
 *
 * Output, one CSV row per (tier, op) on stdout:
 * `backend,threads,op,cells,median_s,runs`.
 *
 * Usage: `meshioplusplus_bench_ops [--tier S|M|L|XL]... [--ops a,b,...] [--runs N]`.
 * The default tiers are S, M and L (about 10k, 160k and 750k tetrahedra);
 * XL (about 10M) is opt-in, since it needs several GB.
 */

// System includes
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <thread>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/mesh_api.hpp"
#include "meshioplusplus/operations/clean.hpp"
#include "meshioplusplus/operations/decimate.hpp"
#include "meshioplusplus/operations/merge.hpp"
#include "meshioplusplus/operations/partition.hpp"
#include "meshioplusplus/operations/refine.hpp"
#include "meshioplusplus/operations/reorder.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/operations/smooth.hpp"
#include "meshioplusplus/operations/surface.hpp"
#include "meshioplusplus/parallel.hpp"

using meshioplusplus::DType;
using meshioplusplus::Mesh;
using meshioplusplus::NDArray;

namespace {

/** @brief Structured tet cube: (n+1)^3 shared vertices, 6 tets per hex (Kuhn). */
Mesh bench_ops_tet_cube(std::size_t n) {
    const std::size_t np = n + 1;
    NDArray pts = NDArray::Uninit(DType::Float64, {np * np * np, 3});
    double* p = pts.As<double>();
    for (std::size_t k = 0; k < np; ++k)
        for (std::size_t j = 0; j < np; ++j)
            for (std::size_t i = 0; i < np; ++i) {
                const std::size_t idx = (k * np + j) * np + i;
                p[idx * 3 + 0] = static_cast<double>(i) / static_cast<double>(n);
                p[idx * 3 + 1] = static_cast<double>(j) / static_cast<double>(n);
                p[idx * 3 + 2] = static_cast<double>(k) / static_cast<double>(n);
            }
    static const int tets[6][4][3] = {
        {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {1, 1, 1}},
        {{0, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {1, 1, 1}},
        {{0, 0, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}}, {{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {1, 1, 1}},
    };
    NDArray conn = NDArray::Uninit(DType::Int64, {6 * n * n * n, 4});
    std::int64_t* c = conn.As<std::int64_t>();
    std::size_t t = 0;
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i)
                for (int s = 0; s < 6; ++s, ++t)
                    for (int v = 0; v < 4; ++v)
                        c[t * 4 + v] = static_cast<std::int64_t>(
                            ((k + tets[s][v][2]) * np + (j + tets[s][v][1])) * np +
                            (i + tets[s][v][0]));
    Mesh m;
    m.AssignPoints(std::move(pts));
    m.AddCellBlock("tetra", std::move(conn));
    return m;
}

double bench_ops_median(const std::function<void()>& rFn, int Runs) {
    rFn();  // warmup
    std::vector<double> ts;
    for (int r = 0; r < Runs; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        rFn();
        ts.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
    std::sort(ts.begin(), ts.end());
    return ts[ts.size() / 2];
}

}  // namespace

int main(int argc, char** argv) {
    namespace mio = meshioplusplus;
    const std::map<std::string, std::size_t> tiers = {{"S", 12}, {"M", 30}, {"L", 50}, {"XL", 120}};
    std::vector<std::string> chosen;
    std::vector<std::string> ops;
    int runs = 3;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--tier" && i + 1 < argc) {
            chosen.push_back(argv[++i]);
        } else if (a == "--ops" && i + 1 < argc) {
            const std::string list = argv[++i];
            for (std::size_t at = 0; at <= list.size();) {
                const std::size_t comma = std::min(list.find(',', at), list.size());
                if (comma > at)
                    ops.push_back(list.substr(at, comma - at));
                at = comma + 1;
            }
        } else if (a == "--runs" && i + 1 < argc) {
            runs = std::max(1, std::atoi(argv[++i]));
        } else {
            std::fprintf(stderr, "usage: %s [--tier S|M|L|XL]... [--ops a,b,...] [--runs N]\n",
                         argv[0]);
            return 2;
        }
    }
    if (chosen.empty())
        chosen = {"S", "M", "L"};
    // OMP_NUM_THREADS is honoured by OpenMP itself; TBB has no such variable,
    // so the same one caps it here, keeping tools/bench_ops.sh's sweep uniform.
    const char* env_threads = std::getenv("OMP_NUM_THREADS");
    const std::string threads =
        env_threads ? env_threads : std::to_string(std::thread::hardware_concurrency());
#if defined(MESHIOPLUSPLUS_PARALLEL_TBB)
    std::optional<tbb::global_control> cap;
    if (env_threads && std::atoi(env_threads) > 0)
        cap.emplace(tbb::global_control::max_allowed_parallelism,
                    static_cast<std::size_t>(std::atoi(env_threads)));
#endif
    const auto wanted = [&](const char* pOp) {
        return ops.empty() || std::find(ops.begin(), ops.end(), pOp) != ops.end();
    };

    std::printf("backend,threads,op,cells,median_s,runs\n");
    for (const std::string& tier : chosen) {
        auto it = tiers.find(tier);
        if (it == tiers.end()) {
            std::fprintf(stderr, "unknown tier '%s'\n", tier.c_str());
            return 2;
        }
        const Mesh volume = bench_ops_tet_cube(it->second);
        const Mesh surface = mio::extract_surface(volume);
        const std::size_t ncells = volume.Cells(0).NumCells();
        const auto row = [&](const char* pOp, const std::function<void()>& rFn) {
            if (!wanted(pOp))
                return;
            const double s = bench_ops_median(rFn, runs);
            std::printf("%s,%s,%s,%zu,%.6f,%d\n", mio::parallel_backend_name(), threads.c_str(),
                        pOp, ncells, s, runs);
            std::fflush(stdout);
        };
        row("extract_surface", [&] { (void)mio::extract_surface(volume); });
        row("smooth", [&] {
            mio::SmoothOptions o;
            o.mIterations = 10;
            (void)mio::smooth(surface, o);
        });
        row("refine", [&] { (void)mio::refine(volume); });
        row("merge", [&] {
            mio::MergeOptions o;
            o.weld = true;
            (void)mio::merge({&volume, &volume}, o);
        });
        row("clean", [&] { (void)mio::clean(volume); });
        row("compute_sdf", [&] {
            mio::SdfOptions o;
            o.mResolution = std::array<std::int64_t, 3>{48, 48, 48};
            (void)mio::compute_sdf(surface, o);
        });
        row("decimate", [&] {
            mio::DecimateOptions o;
            o.mTargetRatio = 0.5;
            (void)mio::decimate(surface, o);
        });
        row("partition", [&] {
            mio::PartitionOptions o;
            o.mNParts = 8;
            (void)mio::partition(volume, o);
        });
        row("reorder", [&] { (void)mio::reorder(volume); });
    }
    return 0;
}
