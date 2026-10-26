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
 * @file fuzz_mesh_cad.cpp
 * @brief libFuzzer harness over the mesh generator rather than a reader
 * (roadmap 6.1.5, doc/fuzzing.md): two targets in one source.
 *
 *  - **`mesh_cad_surface`.** The input is a Kratos CAD JSON file, read with
 *    `read_cad_json_info` (a `ReadError` is a refusal, as for the reader
 *    target) and meshed at a size derived from the input's bytes and the
 *    model's box, so a run stays bounded (`max_points` 200 000). The
 *    contract: a result, with every face meshed or reported by status and
 *    reason, or `std::invalid_argument` (a model with no faces) or
 *    `std::length_error` (the budget). A crash, a sanitizer report, a hang,
 *    any other exception, or a face reported with an internal error of the
 *    triangulation is a finding.
 *  - **`triangulate_polygon_2d`.** The input decodes into one to four loops
 *    of points on a small integer grid, perturbed in the ways inexact
 *    triangulators fail on (Shewchuk 1997; Kettner et al. 2008): points a
 *    few ulps apart, a few ulps off a line, on one line, or raw doubles. The
 *    contract: a mesh whose every triangle is counter-clockwise by the exact
 *    `orient2d`, or `std::invalid_argument` (a degenerate, open or crossing
 *    boundary). Anything else is a finding.
 *
 * The target is taken from, in order: `MIO_FUZZ_TARGET` or `MIO_FUZZ_FORMAT`
 * (what tools/fuzz/run_campaign.sh sets), a `-format=<name>` or
 * `-target=<name>` argument (`op_` prefix allowed, as the regression
 * directories are named), or the binary's name (`meshioplusplus_fuzz_<name>`).
 * `meshioplusplus_fuzz_ops_replay` is the same harness with replay_main.cpp's `main()`, so ctest
 * replays the committed inputs under `tests/fuzz/regressions/op_<name>/`.
 */

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

#include <unistd.h>

// Project includes
#include "../../src/cpp/src/detail/predicates.hpp"
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/cad_json.hpp"
#include "meshioplusplus/operations/mesh_cad.hpp"

namespace {

enum class FuzzOpTarget { None, MeshCadSurface, TriangulatePolygon2d };

FuzzOpTarget gTarget = FuzzOpTarget::None;
std::filesystem::path gDir;
std::filesystem::path gPath;

FuzzOpTarget fuzz_op_parse(std::string Name) {
    if (Name.rfind("op_", 0) == 0)
        Name.erase(0, 3);
    if (Name == "mesh_cad_surface")
        return FuzzOpTarget::MeshCadSurface;
    if (Name == "triangulate_polygon_2d")
        return FuzzOpTarget::TriangulatePolygon2d;
    return FuzzOpTarget::None;
}

FuzzOpTarget fuzz_op_from_args(int argc, char** argv) {
    for (const char* var : {"MIO_FUZZ_TARGET", "MIO_FUZZ_FORMAT"})
        if (const char* env = std::getenv(var); env && *env)
            return fuzz_op_parse(env);
    for (int i = 1; i < argc; ++i)
        for (const char* flag : {"-format=", "-target="})
            if (std::strncmp(argv[i], flag, std::strlen(flag)) == 0)
                return fuzz_op_parse(argv[i] + std::strlen(flag));
    if (argc > 0) {
        const std::string name = std::filesystem::path(argv[0]).filename().string();
        const std::string prefix = "meshioplusplus_fuzz_";
        if (name.rfind(prefix, 0) == 0)
            return fuzz_op_parse(name.substr(prefix.size()));
    }
    return FuzzOpTarget::None;
}

void fuzz_op_cleanup() {
    std::error_code ec;
    std::filesystem::remove_all(gDir, ec);
}

[[noreturn]] void fuzz_op_finding(const char* pWhat, const std::string& rDetail) {
    std::fprintf(stderr, "meshioplusplus_fuzz_ops: %s: %s\n", pWhat, rDetail.c_str());
    std::abort();
}

/// FNV-1a over the input: the size and options are a function of the bytes,
/// so a finding replays exactly.
std::uint64_t fuzz_op_hash(const std::uint8_t* pData, std::size_t Size) {
    std::uint64_t h = 1469598103934665603ull;
    for (std::size_t i = 0; i < Size; ++i) {
        h ^= pData[i];
        h *= 1099511628211ull;
    }
    return h;
}

void fuzz_mesh_cad_surface(const std::uint8_t* pData, std::size_t Size) {
    {
        std::ofstream out(gPath, std::ios::binary | std::ios::trunc);
        out.write(reinterpret_cast<const char*>(pData), static_cast<std::streamsize>(Size));
    }
    meshioplusplus::BrepInfo info;
    try {
        info = meshioplusplus::read_cad_json_info(gPath.string());
    } catch (const meshioplusplus::ReadError&) {
        return;  // the reader's own refusal; fuzz_read covers the parser
    }
    // The size: the control nets' box over 4 to 64 pieces, chosen by the
    // input's hash, with the chordal cap on for half the inputs.
    double lo[3] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[3] = {-lo[0], -lo[1], -lo[2]};
    for (const meshioplusplus::NurbsPatch& r_p : info.mPatches) {
        const meshioplusplus::NDArray& c = r_p.mControlPoints;
        if (c.Shape().size() != 2 || c.Shape()[1] != 3)
            continue;
        for (std::size_t i = 0; i < 3 * c.Shape()[0]; ++i) {
            const double v = meshioplusplus::detail::read_double(c, i);
            if (std::isfinite(v)) {
                lo[i % 3] = std::min(lo[i % 3], v);
                hi[i % 3] = std::max(hi[i % 3], v);
            }
        }
    }
    double diag = 0.0;
    if (lo[0] <= hi[0])
        diag = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
                         (hi[2] - lo[2]) * (hi[2] - lo[2]));
    const std::uint64_t hash = fuzz_op_hash(pData, Size);
    meshioplusplus::MeshCadSurfaceOptions o;
    o.mTargetSize = std::isfinite(diag) && diag > 0.0 ? diag / double(4 + hash % 61) : 1.0;
    if (!(o.mTargetSize > 0.0) || !std::isfinite(o.mTargetSize))
        o.mTargetSize = 1.0;
    if ((hash >> 8) & 1)
        o.mMaxChordalDeviation = 0.01 * o.mTargetSize;
    // Refinement (v17.15.0) for three inputs in four, smoothing for two in
    // three, so the refinement loop and the guarded moves see mutated models.
    static constexpr double kAngles[4] = {0.0, 20.0, 25.0, 30.0};
    o.mMinAngleDeg = kAngles[(hash >> 16) % 4];
    o.mSmoothIterations = static_cast<std::int64_t>((hash >> 20) % 3);
    o.mSmoothMethod = ((hash >> 24) & 1) ? meshioplusplus::MeshCadSmoothing::Odt
                                         : meshioplusplus::MeshCadSmoothing::Laplacian;
    // Curved elements (v17.16.0) of order 1, 2 or 3, a third of the inputs
    // each, so the node placement and the validity check see them too.
    o.mOrder = 1 + static_cast<std::int32_t>((hash >> 28) % 3);
    // Sizing and pole rings (v17.18.0), each for one input in four, so the
    // first pass, the field and the ring placement see mutated models.
    if (((hash >> 32) & 3) == 0)
        o.mCurvatureDeviation = 0.02 * o.mTargetSize;
    if (((hash >> 34) & 3) == 0)
        o.mCellsAcrossGaps = 2.0;
    if (((hash >> 36) & 3) == 0)
        o.mPoleRings = 2;
    // A mesh must fit libFuzzer's 10 s timeout under ASan, UBSan and
    // coverage, about ten times slower than a release build. A face sized to
    // a quarter of h and refined to 30 degrees (the reweighted-edge seed)
    // took 49 s instrumented at 200000 points; at 10000 it takes 5 s, and
    // below about 8500 its lattice alone would be refused.
    o.mMaxPoints = 10000;
    // Given sizes (v17.19.0) for one input in four: the first face and the
    // first edge at a quarter of the target size.
    meshioplusplus::MeshCadSizeSources sources;
    if (((hash >> 38) & 3) == 0) {
        if (!info.mFaces.empty())
            sources.mFaceSizes = {{0, 0.25 * o.mTargetSize}};
        if (!info.mEdges.empty())
            sources.mEdgeSizes = {{0, 0.25 * o.mTargetSize}};
    }
    // Virtual edges and a varying metric (v17.19.0), each for one input in
    // four (virtual edges only at order 1, where they are allowed).
    sources.mVirtualEdges = ((hash >> 40) & 3) == 0 && o.mOrder == 1;
    if (((hash >> 42) & 3) == 0)
        o.mMetric = meshioplusplus::MeshCadMetric::Varying;
    // The discrete reparametrisation for one input in four, virtual faces and
    // model tolerances for one in eight each (v17.20.0; the first two order 1
    // only).
    meshioplusplus::MeshCadSurfaceExtras extras;
    if (((hash >> 44) & 3) == 0 && o.mOrder == 1)
        o.mReparamConditionNumber = 4.0;
    extras.mVirtualFaces = ((hash >> 46) & 7) == 0 && o.mOrder == 1;
    extras.mModelTolerances = ((hash >> 49) & 7) == 0;
    try {
        const meshioplusplus::MeshCadSurfaceResult r =
            meshioplusplus::mesh_cad_surface(info, o, sources, extras);
        for (const meshioplusplus::MeshCadFaceRecord& r_f : r.mFaces)
            if (r_f.mStatus == meshioplusplus::MeshCadFaceStatus::Other &&
                r_f.mReason.find("internal error") != std::string::npos)
                fuzz_op_finding("a face reported an internal error", r_f.mReason);
    } catch (const std::invalid_argument&) {
    } catch (const std::length_error&) {
    } catch (const std::exception& e) {
        fuzz_op_finding("mesh_cad_surface leaked an exception",
                        std::string(typeid(e).name()) + ": " + e.what());
    }
}

/// Bytes read in order; zero past the end.
class FuzzOpBytes {
public:
    FuzzOpBytes(const std::uint8_t* pData, std::size_t Size) : pData(pData), mSize(Size) {}
    std::uint8_t Byte() { return mAt < mSize ? pData[mAt++] : 0; }
    double Raw() {
        std::uint64_t bits = 0;
        for (int k = 0; k < 8; ++k)
            bits |= static_cast<std::uint64_t>(Byte()) << (8 * k);
        double d = 0.0;
        std::memcpy(&d, &bits, sizeof d);
        return d;
    }

private:
    const std::uint8_t* pData;
    std::size_t mSize;
    std::size_t mAt = 0;
};

/// `X` moved by `K` ulps (towards +inf for K > 0).
double fuzz_op_ulps(double X, int K) {
    for (; K > 0; --K)
        X = std::nextafter(X, std::numeric_limits<double>::infinity());
    for (; K < 0; ++K)
        X = std::nextafter(X, -std::numeric_limits<double>::infinity());
    return X;
}

void fuzz_triangulate_polygon_2d(const std::uint8_t* pData, std::size_t Size) {
    using Vec2 = std::array<double, 2>;
    FuzzOpBytes in(pData, Size);
    const int loops = 1 + in.Byte() % 4;
    const int mode = in.Byte() % 4;
    const int size_sel = in.Byte() % 8;
    std::vector<std::vector<Vec2>> poly;
    double lo[2] = {std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()};
    double hi[2] = {-lo[0], -lo[1]};
    for (int l = 0; l < loops; ++l) {
        const int n = 3 + in.Byte() % 30;
        std::vector<Vec2> loop;
        for (int i = 0; i < n; ++i) {
            Vec2 p{};
            if (mode == 3) {
                p = {in.Raw(), in.Raw()};
            } else {
                p = {static_cast<std::int8_t>(in.Byte()) / 16.0,
                     static_cast<std::int8_t>(in.Byte()) / 16.0};
                if (mode == 1) {  // a few ulps off the grid: near-duplicates, near-collinear
                    p[0] = fuzz_op_ulps(p[0], static_cast<std::int8_t>(in.Byte()) % 5);
                    p[1] = fuzz_op_ulps(p[1], static_cast<std::int8_t>(in.Byte()) % 5);
                } else if (mode == 2) {  // on the line y = x / 3, which no double represents
                    p[1] = p[0] / 3.0;
                    p[1] = fuzz_op_ulps(p[1], static_cast<std::int8_t>(in.Byte()) % 3);
                }
            }
            loop.push_back(p);
            for (int k = 0; k < 2; ++k)
                if (std::isfinite(p[static_cast<std::size_t>(k)])) {
                    lo[k] = std::min(lo[k], p[static_cast<std::size_t>(k)]);
                    hi[k] = std::max(hi[k], p[static_cast<std::size_t>(k)]);
                }
        }
        poly.push_back(std::move(loop));
    }
    const double diag = lo[0] <= hi[0] ? std::hypot(hi[0] - lo[0], hi[1] - lo[1]) : 0.0;
    // 0 (the boundary alone), or the box over 4 to 128 pieces: bounded work.
    double h = 0.0;
    if (size_sel > 0 && std::isfinite(diag) && diag > 0.0)
        h = diag / double(4 << (size_sel - 1));
    if (!std::isfinite(h))
        h = 0.0;
    // A byte past the loops (0 past the end, so older inputs keep their
    // meaning) asks for refinement (v17.15.0), bounded by max_points.
    static constexpr double kAngles[4] = {0.0, 20.7, 25.0, 30.0};
    meshioplusplus::TriangulatePolygon2dOptions o;
    o.mTargetSize = h;
    o.mMinAngleDeg = kAngles[in.Byte() % 4];
    o.mMaxPoints = 20000;
    try {
        const meshioplusplus::Mesh m = meshioplusplus::triangulate_polygon_2d(poly, o);
        const double* x = m.Points().As<double>();
        const auto cells = m.Cells(0);
        const meshioplusplus::NDArray& c = cells.Conn();
        for (std::size_t t = 0; t < cells.NumCells(); ++t) {
            const double* a = x + 3 * meshioplusplus::detail::read_int(c, 3 * t);
            const double* b = x + 3 * meshioplusplus::detail::read_int(c, 3 * t + 1);
            const double* d = x + 3 * meshioplusplus::detail::read_int(c, 3 * t + 2);
            if (!(meshioplusplus::detail::orient2d(a, b, d) > 0.0))
                fuzz_op_finding("triangulate_polygon_2d", "a triangle is not counter-clockwise");
        }
    } catch (const std::invalid_argument&) {
    } catch (const std::exception& e) {
        fuzz_op_finding("triangulate_polygon_2d leaked an exception",
                        std::string(typeid(e).name()) + ": " + e.what());
    }
}

}  // namespace

/**
 * @brief Picks the target and creates the per-process scratch directory.
 * @return 0; exits with status 2 on an unknown target.
 */
extern "C" int LLVMFuzzerInitialize(int* pArgc, char*** pArgv) {
    gTarget = fuzz_op_from_args(*pArgc, *pArgv);
    if (gTarget == FuzzOpTarget::None) {
        std::fprintf(stderr,
                     "meshioplusplus_fuzz_ops: unknown target; known: "
                     "mesh_cad_surface triangulate_polygon_2d\n");
        std::exit(2);
    }
    gDir = std::filesystem::temp_directory_path() /
           ("mio-fuzz-ops-" + std::to_string(static_cast<long>(::getpid())));
    std::filesystem::create_directories(gDir);
    gPath = gDir / "input.cad.json";
    std::atexit(fuzz_op_cleanup);
    return 0;
}

/** @brief Runs one input through the chosen operation; aborts on a finding. */
extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* pData, std::size_t size) {
    if (gTarget == FuzzOpTarget::MeshCadSurface)
        fuzz_mesh_cad_surface(pData, size);
    else
        fuzz_triangulate_polygon_2d(pData, size);
    return 0;
}
