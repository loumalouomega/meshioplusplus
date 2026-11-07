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
 * @file test_library_lock.cpp
 * @brief The process-wide lock that serialises HDF5, netCDF and the other
 * non-thread-safe libraries (roadmap 3.4.1.1, `detail/library_lock.hpp`).
 *
 * Several threads write and read their own files at once: without the lock the
 * HDF5 library's unsynchronised global state corrupts a file or crashes (run
 * under the sanitizers to see it), with it every round trip matches a serial one.
 */

// System includes
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <string>
#include <thread>
#include <vector>

// External includes
#include <gtest/gtest.h>

// Project includes
#include "../../src/cpp/src/detail/library_lock.hpp"
#include "mesh_fixtures.hpp"
#include "meshioplusplus/formats/vtkhdf.hpp"
#include "meshioplusplus/formats/vtkhdf_time_series.hpp"
#include "meshioplusplus/formats/xdmf.hpp"
#include "meshioplusplus/formats/xdmf_time_series.hpp"
#ifdef MESHIOPLUSPLUS_HAS_NETCDF
#include "meshioplusplus/formats/exodus.hpp"
#endif

namespace {

constexpr int kThreads = 8;

mt::Mesh lock_mesh(int Seed) {
    mt::Mesh m;
    m.AssignPoints(mt::points_from({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}}));
    mt::NDArray conn = mt::NDArray::Uninit(meshioplusplus::DType::Int64, {2, 4});
    const std::int64_t rows[8] = {0, 1, 2, 3, 1, 2, 3, 4};
    std::copy(rows, rows + 8, conn.As<std::int64_t>());
    m.AddCellBlock("tetra", std::move(conn));
    const double v = Seed;
    m.AddPointData("u", mt::data_array({v, v + 1, v + 2, v + 3, v + 4}));
    return m;
}

double first_u(const mt::Mesh& rMesh) {
    return meshioplusplus::detail::read_double(rMesh.PointData("u"), 0);
}

/// Runs `Body(i, path)` on `kThreads` threads, each with its own file, and
/// counts the threads whose body threw or returned false.
template <class TBody>
int run_threads(const std::string& rSuffix, TBody Body) {
    std::vector<std::string> paths;
    for (int i = 0; i < kThreads; ++i)
        paths.push_back(mt::temp_path(rSuffix));
    std::atomic<int> failures{0};
    std::vector<std::thread> threads;
    for (int i = 0; i < kThreads; ++i)
        threads.emplace_back([&, i] {
            try {
                if (!Body(i, paths[static_cast<std::size_t>(i)]))
                    ++failures;
            } catch (...) {
                ++failures;
            }
        });
    for (std::thread& t : threads)
        t.join();
    std::error_code ec;
    for (const std::string& p : paths) {
        std::filesystem::remove(p, ec);
        std::filesystem::remove(std::filesystem::path(p).replace_extension(".h5"), ec);
    }
    return failures.load();
}

}  // namespace

TEST(LibraryLock, IsRecursive) {
    meshioplusplus::detail::LibraryLock outer;
    meshioplusplus::detail::LibraryLock inner;  // the same thread nests: must not deadlock
    SUCCEED();
}

TEST(LibraryLock, ExcludesOtherThreads) {
    std::atomic<bool> got{false};
    std::thread t;
    {
        meshioplusplus::detail::LibraryLock held;
        t = std::thread([&] {
            meshioplusplus::detail::LibraryLock lock;
            got = true;
        });
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        EXPECT_FALSE(got.load());
    }  // `held` drops here; the thread can proceed
    t.join();
    EXPECT_TRUE(got.load());
}

// Roadmap 5.8.2: the parallel deflate and inflate give the lock up. Every level
// the thread holds is released and taken back, and another thread gets in.
TEST(LibraryLock, UnlockReleasesEveryLevelAndRestoresThem) {
    using meshioplusplus::detail::library_lock_depth;
    EXPECT_EQ(library_lock_depth(), 0);
    {
        const meshioplusplus::detail::LibraryUnlock none;  // holding nothing: a no-op
        EXPECT_EQ(library_lock_depth(), 0);
    }
    meshioplusplus::detail::LibraryLock outer;
    meshioplusplus::detail::LibraryLock inner;
    EXPECT_EQ(library_lock_depth(), 2);
    std::atomic<bool> got{false};
    {
        const meshioplusplus::detail::LibraryUnlock unlock;
        EXPECT_EQ(library_lock_depth(), 0);
        std::thread t([&] {
            meshioplusplus::detail::LibraryLock lock;
            got = true;
        });
        t.join();  // would deadlock if a level were still held
    }
    EXPECT_TRUE(got.load());
    EXPECT_EQ(library_lock_depth(), 2);
    // Re-taken: another thread waits again.
    std::atomic<bool> got2{false};
    std::thread t2([&] {
        meshioplusplus::detail::LibraryLock lock;
        got2 = true;
    });
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    EXPECT_FALSE(got2.load());
    {
        // Leave the scope's two levels by hand so t2 can finish before the
        // guards above unwind.
        const meshioplusplus::detail::LibraryUnlock unlock;
        t2.join();
    }
    EXPECT_TRUE(got2.load());
}

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// Large enough that the gzip datasets span several chunks, so the writes and
// reads compress and inflate with the lock given up while other threads use
// HDF5.
TEST(LibraryLock, ConcurrentGzipVtkhdfRoundTripsWithTheLockGivenUp) {
    EXPECT_EQ(0, run_threads(".vtkhdf", [](int I, const std::string& rPath) {
                  const std::size_t n = 300000;  // point rows: 7 MB of float64
                  mt::NDArray pts = mt::NDArray::Uninit(meshioplusplus::DType::Float64, {n, 3});
                  mt::NDArray u = mt::NDArray::Uninit(meshioplusplus::DType::Float64, {n});
                  for (std::size_t k = 0; k < n; ++k) {
                      for (int d = 0; d < 3; ++d)
                          pts.As<double>()[k * 3 + d] = static_cast<double>((k * (d + 3)) % 1000);
                      u.As<double>()[k] = static_cast<double>(I) + static_cast<double>(k % 7);
                  }
                  mt::Mesh big;
                  big.AssignPoints(std::move(pts));
                  mt::NDArray conn = mt::NDArray::Uninit(meshioplusplus::DType::Int64, {1, 4});
                  const std::int64_t row[4] = {0, 1, 2, 3};
                  std::copy(row, row + 4, conn.As<std::int64_t>());
                  big.AddCellBlock("tetra", std::move(conn));
                  big.AddPointData("u", std::move(u));
                  meshioplusplus::write_vtkhdf(rPath, big, 6);
                  const mt::Mesh back = meshioplusplus::read_vtkhdf(rPath);
                  if (back.NumPoints() != n)
                      return false;
                  const mt::NDArray& bu = back.PointData("u");
                  for (std::size_t k = 0; k < n; k += 997)
                      if (meshioplusplus::detail::read_double(bu, k) !=
                          static_cast<double>(I) + static_cast<double>(k % 7))
                          return false;
                  return true;
              }));
}

TEST(LibraryLock, ConcurrentVtkhdfRoundTrips) {
    EXPECT_EQ(0, run_threads(".vtkhdf", [](int I, const std::string& rPath) {
                  meshioplusplus::write_vtkhdf(rPath, lock_mesh(I));
                  const mt::Mesh back = meshioplusplus::read_vtkhdf(rPath);
                  return first_u(back) == static_cast<double>(I);
              }));
}

TEST(LibraryLock, ConcurrentXdmfHdfRoundTrips) {
    EXPECT_EQ(0, run_threads(".xdmf", [](int I, const std::string& rPath) {
                  meshioplusplus::write_xdmf(rPath, lock_mesh(I), "HDF", -1);
                  const mt::Mesh back = meshioplusplus::read_xdmf(rPath);
                  return first_u(back) == static_cast<double>(I);
              }));
}

TEST(LibraryLock, ConcurrentSeriesWritersAndDestructors) {
    EXPECT_EQ(0, run_threads(".vtkhdf", [](int I, const std::string& rPath) {
                  {
                      // Never finalized by hand: the destructor does it, under the lock.
                      meshioplusplus::VtkhdfTimeSeriesWriter w(rPath);
                      w.WritePointsCells(lock_mesh(I));
                      for (int k = 0; k < 3; ++k)
                          w.WriteData(0.5 * k, lock_mesh(I + k));
                  }
                  meshioplusplus::ReadOptions o;
                  o.mTimeStep = 2;
                  return first_u(meshioplusplus::read_vtkhdf(rPath, o)) ==
                         static_cast<double>(I + 2);
              }));
}

TEST(LibraryLock, ConcurrentXdmfSeriesWritersAndDestructors) {
    EXPECT_EQ(0, run_threads(".xdmf", [](int I, const std::string& rPath) {
                  {
                      meshioplusplus::XdmfTimeSeriesWriter w(rPath, "HDF");
                      w.WritePointsCells(lock_mesh(I));
                      for (int k = 0; k < 3; ++k)
                          w.WriteData(0.5 * k, lock_mesh(I + k));
                  }
                  meshioplusplus::ReadOptions o;
                  o.mTimeStep = 2;
                  return first_u(meshioplusplus::read_xdmf(rPath, o)) == static_cast<double>(I + 2);
              }));
}

TEST(LibraryLock, SeriesWriterNestsInsideAHeldLock) {
    // A caller (a pipeline, a binding) may already hold the lock: the writer's own
    // entry points take it again on the same thread.
    const std::string path = mt::temp_path(".vtkhdf");
    {
        meshioplusplus::detail::LibraryLock held;
        meshioplusplus::VtkhdfTimeSeriesWriter w(path);
        w.WritePointsCells(lock_mesh(0));
        w.WriteData(0.0, lock_mesh(0));
        w.Finalize();
    }
    EXPECT_EQ(first_u(meshioplusplus::read_vtkhdf(path)), 0.0);
    std::filesystem::remove(path);
}

#endif  // MESHIOPLUSPLUS_HAS_HDF5

#ifdef MESHIOPLUSPLUS_HAS_NETCDF

TEST(LibraryLock, ConcurrentExodusRoundTripsAndSeriesWriters) {
    EXPECT_EQ(0, run_threads(".exo", [](int I, const std::string& rPath) {
                  meshioplusplus::write_exodus(rPath, lock_mesh(I));
                  const mt::Mesh back = meshioplusplus::read_exodus(rPath);
                  if (first_u(back) != static_cast<double>(I))
                      return false;
                  const std::string series = rPath + ".series.exo";
                  {
                      meshioplusplus::ExodusTimeSeriesWriter w(series);
                      w.WritePointsCells(lock_mesh(I));
                      w.WriteData(0.0, lock_mesh(I));
                  }
                  const bool ok = meshioplusplus::read_exodus(series).NumPoints() == 5;
                  std::error_code ec;
                  std::filesystem::remove(series, ec);
                  return ok;
              }));
}

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
