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
// The surface-boolean spike's measurement harness (roadmap 7.8,
// doc/spikes.md): runs union, intersection and difference on a pair of
// closed surfaces and prints one JSON line per operation with what the spike
// measures. Not an installed program; built with the tests.
//
//   meshioplusplus_boolean_spike A B [--write PREFIX] [--case NAME]
//   meshioplusplus_boolean_spike --synthetic SEED N [--only NAME] [--write PREFIX]
//
// The synthetic corpus is generated here from a seed (spheres, boxes, tori
// under seeded rigid motions, and the degenerate contacts) so a run repeats
// on any machine; tools/boolean_corpus.py drives the Thingi10K pairs.

// System includes
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/operations/find_intersections.hpp"
#include "meshioplusplus/operations/primitives.hpp"
#include "meshioplusplus/operations/sdf.hpp"
#include "meshioplusplus/formats/vtu.hpp"
#include "meshioplusplus/registry.hpp"

// Project includes (private, not installed)
#include "../src/detail/mesh_boolean.hpp"

namespace mio = meshioplusplus;

namespace {

using BsVec = std::array<double, 3>;

/// A portable generator: splitmix64, so a seed gives the same corpus everywhere.
struct BsRng {
    std::uint64_t mState;
    std::uint64_t Next() {
        std::uint64_t z = (mState += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    double Uniform(double Lo, double Hi) {
        return Lo + (Hi - Lo) * static_cast<double>(Next() >> 11) * (1.0 / 9007199254740992.0);
    }
};

mio::Mesh bs_triangles(const std::vector<BsVec>& rP,
                       const std::vector<std::array<std::int64_t, 3>>& rT) {
    mio::Mesh m;
    mio::NDArray p = mio::NDArray::Uninit(mio::DType::Float64, {rP.size(), std::size_t{3}});
    for (std::size_t i = 0; i < rP.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            p.As<double>()[3 * i + k] = rP[i][k];
    m.AssignPoints(std::move(p));
    mio::NDArray c = mio::NDArray::Uninit(mio::DType::Int64, {rT.size(), std::size_t{3}});
    for (std::size_t t = 0; t < rT.size(); ++t)
        for (std::size_t k = 0; k < 3; ++k)
            c.As<std::int64_t>()[3 * t + k] = rT[t][k];
    m.AddCellBlock("triangle", std::move(c));
    return m;
}

/// The box [Lo, Hi] with each face an N x N grid of triangle pairs, outward.
mio::Mesh bs_box(const BsVec& rLo, const BsVec& rHi, int N) {
    std::vector<BsVec> pts;
    std::vector<std::array<std::int64_t, 3>> tris;
    std::vector<std::int64_t> index((N + 1) * (N + 1) * (N + 1), -1);
    const auto id = [&](int I, int J, int K) {
        std::int64_t& r = index[static_cast<std::size_t>((I * (N + 1) + J) * (N + 1) + K)];
        if (r < 0) {
            r = static_cast<std::int64_t>(pts.size());
            pts.push_back({rLo[0] + (rHi[0] - rLo[0]) * I / N, rLo[1] + (rHi[1] - rLo[1]) * J / N,
                           rLo[2] + (rHi[2] - rLo[2]) * K / N});
        }
        return r;
    };
    // Each face: fixed axis `a` at side `s`, the grid over the other two.
    for (int a = 0; a < 3; ++a)
        for (int s = 0; s < 2; ++s)
            for (int u = 0; u < N; ++u)
                for (int v = 0; v < N; ++v) {
                    const auto corner = [&](int Du, int Dv) {
                        int c[3];
                        c[a] = s * N;
                        c[(a + 1) % 3] = u + Du;
                        c[(a + 2) % 3] = v + Dv;
                        return id(c[0], c[1], c[2]);
                    };
                    const std::int64_t p00 = corner(0, 0), p10 = corner(1, 0), p11 = corner(1, 1),
                                       p01 = corner(0, 1);
                    // (a+1, a+2, a) is right-handed: the face at s = 1 turns
                    // counter-clockwise seen from outside.
                    if (s == 1) {
                        tris.push_back({p00, p10, p11});
                        tris.push_back({p00, p11, p01});
                    } else {
                        tris.push_back({p00, p11, p10});
                        tris.push_back({p00, p01, p11});
                    }
                }
    return bs_triangles(pts, tris);
}

/// A torus about z: major radius R, minor r, Nu x Nv quads split, outward.
mio::Mesh bs_torus(double R, double r, int Nu, int Nv) {
    const double pi = 3.14159265358979323846;
    std::vector<BsVec> pts;
    std::vector<std::array<std::int64_t, 3>> tris;
    for (int i = 0; i < Nu; ++i)
        for (int j = 0; j < Nv; ++j) {
            const double u = 2.0 * pi * i / Nu, v = 2.0 * pi * j / Nv;
            pts.push_back({(R + r * std::cos(v)) * std::cos(u), (R + r * std::cos(v)) * std::sin(u),
                           r * std::sin(v)});
        }
    const auto at = [&](int I, int J) {
        return static_cast<std::int64_t>(((I + Nu) % Nu) * Nv + (J + Nv) % Nv);
    };
    for (int i = 0; i < Nu; ++i)
        for (int j = 0; j < Nv; ++j) {
            tris.push_back({at(i, j), at(i + 1, j), at(i + 1, j + 1)});
            tris.push_back({at(i, j), at(i + 1, j + 1), at(i, j + 1)});
        }
    return bs_triangles(pts, tris);
}

/// The mesh moved: rotated by the unit quaternion `rQ`, then translated.
mio::Mesh bs_move(const mio::Mesh& rM, const std::array<double, 4>& rQ, const BsVec& rT) {
    const double w = rQ[0], x = rQ[1], y = rQ[2], z = rQ[3];
    const double r[3][3] = {{1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)},
                            {2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)},
                            {2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)}};
    std::vector<BsVec> pts(rM.NumPoints());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const double px = mio::detail::read_double(rM.Points(), 3 * i);
        const double py = mio::detail::read_double(rM.Points(), 3 * i + 1);
        const double pz = mio::detail::read_double(rM.Points(), 3 * i + 2);
        for (std::size_t k = 0; k < 3; ++k)
            pts[i][k] = r[k][0] * px + r[k][1] * py + r[k][2] * pz + rT[k];
    }
    std::vector<std::array<std::int64_t, 3>> tris;
    const auto cells = rM.Cells(0);
    for (std::size_t c = 0; c < cells.NumCells(); ++c)
        tris.push_back({mio::detail::read_int(cells.Conn(), 3 * c),
                        mio::detail::read_int(cells.Conn(), 3 * c + 1),
                        mio::detail::read_int(cells.Conn(), 3 * c + 2)});
    return bs_triangles(pts, tris);
}

std::array<double, 4> bs_quat(BsRng& rRng) {
    double q[4], n = 0.0;
    do {
        n = 0.0;
        for (double& r_c : q) {
            r_c = rRng.Uniform(-1.0, 1.0);
            n += r_c * r_c;
        }
    } while (n > 1.0 || n < 1e-6);
    n = std::sqrt(n);
    return {q[0] / n, q[1] / n, q[2] / n, q[3] / n};
}

BsVec bs_dir(BsRng& rRng) {
    const auto q = bs_quat(rRng);
    return {q[1], q[2], q[3]};
}

std::string bs_num(double V) {
    char buf[40];
    mio::detail::snprintf_c(buf, sizeof(buf), "%.9g", V);
    return buf;
}

std::string bs_escape(const std::string& rS) {
    std::string out;
    for (const char c : rS)
        out += (c == '"' || c == '\\') ? std::string("\\") + c : std::string(1, c);
    return out;
}

/// `--only` and `--write` of the synthetic mode: run only the cases whose
/// name holds `g_only`, and write their inputs and results under `g_write`.
std::string g_only, g_write;

/// Run the three operations on one pair and print a JSON line each.
void bs_run(const std::string& rCase, const mio::Mesh& rA, const mio::Mesh& rB,
            std::string rWrite) {
    if (!g_only.empty() && rCase.find(g_only) == std::string::npos)
        return;
    if (rWrite.empty() && !g_write.empty()) {
        rWrite = g_write + rCase;
        mio::write_vtu(rWrite + "_a.vtu", rA, false, false);
        mio::write_vtu(rWrite + "_b.vtu", rB, false, false);
    }
    const double va = mio::detail::closed_surface_volume(rA);
    const double vb = mio::detail::closed_surface_volume(rB);
    double vol[3] = {0.0, 0.0, 0.0};
    bool ok[3] = {false, false, false};
    struct Row {
        std::string mJson;
    };
    std::vector<std::string> rows;
    const char* names[3] = {"union", "intersection", "difference"};
    const mio::detail::BooleanOp ops[3] = {mio::detail::BooleanOp::Union,
                                           mio::detail::BooleanOp::Intersection,
                                           mio::detail::BooleanOp::Difference};
    for (int o = 0; o < 3; ++o) {
        const auto t0 = std::chrono::steady_clock::now();
        mio::detail::BooleanResult r = mio::detail::mesh_boolean_spike(rA, rB, ops[o]);
        const double secs =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        bool watertight = false, clean = false;
        if (r.mOk) {
            ok[o] = true;
            vol[o] = mio::detail::closed_surface_volume(r.mMesh);
            if (r.mReport.mNumOutputTriangles > 0) {
                watertight = mio::surface_watertight_check(r.mMesh).mWatertight;
                mio::FindIntersectionsOptions fo;
                fo.mAttachFlags = false;
                fo.mSegments = false;
                clean = mio::find_intersections(r.mMesh, fo).mClean;
            } else {
                watertight = clean = true;  // an empty result is closed
            }
            if (!rWrite.empty() && r.mReport.mNumOutputTriangles > 0)
                mio::write_vtu(rWrite + "_" + names[o] + ".vtu", r.mMesh, false, false);
        }
        const auto& rep = r.mReport;
        std::string j = "{\"case\":\"" + bs_escape(rCase) + "\",\"op\":\"" + names[o] +
                        "\",\"ok\":" + (r.mOk ? "true" : "false") + ",\"failure\":\"" +
                        bs_escape(r.mFailure) +
                        "\",\"watertight\":" + (watertight ? "true" : "false") +
                        ",\"self_intersection_free\":" + (clean ? "true" : "false") +
                        ",\"seconds\":" + bs_num(secs) +
                        ",\"triangles_a\":" + std::to_string(rA.Cells(0).NumCells()) +
                        ",\"triangles_b\":" + std::to_string(rB.Cells(0).NumCells()) +
                        ",\"intersecting_pairs\":" + std::to_string(rep.mNumIntersectingPairs) +
                        ",\"coplanar_pairs\":" + std::to_string(rep.mNumCoplanarPairs) +
                        ",\"implicit_points\":" + std::to_string(rep.mNumImplicitPoints) +
                        ",\"segments\":" + std::to_string(rep.mNumSegments) +
                        ",\"cut_triangles\":" + std::to_string(rep.mNumCutTriangles) +
                        ",\"patches\":" + std::to_string(rep.mNumPatchesA + rep.mNumPatchesB) +
                        ",\"ray_retries\":" + std::to_string(rep.mNumRayRetries) +
                        ",\"label_conflicts\":" + std::to_string(rep.mNumLabelConflicts) +
                        ",\"output_triangles\":" + std::to_string(rep.mNumOutputTriangles) +
                        ",\"rounding_inverted\":" + std::to_string(rep.mNumRoundingInverted) +
                        ",\"seconds_pairs\":" + bs_num(rep.mSecondsPairs) +
                        ",\"seconds_arrangement\":" + bs_num(rep.mSecondsArrangement) +
                        ",\"seconds_triangulate\":" + bs_num(rep.mSecondsTriangulate) +
                        ",\"seconds_classify\":" + bs_num(rep.mSecondsClassify) +
                        ",\"volume\":" + bs_num(vol[o]);
        rows.push_back(j);
    }
    // The volume identities, when all three ran.
    const double scale = std::fabs(va) + std::fabs(vb);
    const double e_sum = ok[0] && ok[1] ? std::fabs(vol[0] + vol[1] - va - vb) / scale : -1.0;
    const double e_diff = ok[1] && ok[2] ? std::fabs(vol[2] - (va - vol[1])) / scale : -1.0;
    for (std::string& r_j : rows)
        std::cout << r_j << ",\"volume_a\":" << bs_num(va) << ",\"volume_b\":" << bs_num(vb)
                  << ",\"identity_error_union\":" << bs_num(e_sum)
                  << ",\"identity_error_difference\":" << bs_num(e_diff) << "}\n";
}

void bs_synthetic(std::uint64_t Seed, int N) {
    BsRng rng{Seed};
    // The degenerate contacts first, fixed.
    const mio::Mesh unit = bs_box({0, 0, 0}, {1, 1, 1}, 2);
    bs_run("axis_boxes_overlap", unit, bs_box({0.5, 0.5, 0.5}, {1.5, 1.5, 1.5}, 2), "");
    bs_run("boxes_shared_face", unit, bs_box({1, 0, 0}, {2, 1, 1}, 2), "");
    bs_run("boxes_overlapping_coplanar_face", unit, bs_box({0.5, 0, 0.5}, {1.5, 1, 1.5}, 2), "");
    const mio::Mesh sphere = mio::icosa_sphere(1.0, 2);
    bs_run("identical_spheres", sphere, sphere, "");
    bs_run("nested_spheres", sphere, mio::icosa_sphere(0.4, 2), "");
    bs_run("disjoint_spheres", sphere, mio::icosa_sphere(1.0, 2, {3.0, 0.0, 0.0}), "");
    bs_run("touching_spheres_at_a_vertex", sphere, bs_move(sphere, {1, 0, 0, 0}, {0.0, 0.0, 2.0}),
           "");
    bs_run("box_vertex_on_box_face", unit,
           bs_move(unit, {0.9238795325112867, 0.0, 0.0, 0.3826834323650898}, {0.5, -0.2, 0.25}),
           "");
    // The seeded classes.
    for (int i = 0; i < N; ++i) {
        const int kind = i % 5;
        const std::string name = "seed" + std::to_string(Seed) + "_" + std::to_string(i);
        if (kind == 0) {
            const BsVec d = bs_dir(rng);
            const double s = rng.Uniform(0.2, 1.8);
            bs_run(name + "_spheres_offset", mio::icosa_sphere(1.0, 3),
                   mio::icosa_sphere(rng.Uniform(0.5, 1.2), 2, {s * d[0], s * d[1], s * d[2]}), "");
        } else if (kind == 1) {
            const BsVec d = bs_dir(rng);
            bs_run(name + "_sphere_box", mio::icosa_sphere(1.0, 3),
                   bs_move(bs_box({-0.6, -0.6, -0.6}, {0.6, 0.6, 0.6}, 4), bs_quat(rng),
                           {0.8 * d[0], 0.8 * d[1], 0.8 * d[2]}),
                   "");
        } else if (kind == 2) {
            const BsVec d = bs_dir(rng);
            bs_run(name + "_boxes_rotated", bs_box({-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5}, 3),
                   bs_move(bs_box({-0.5, -0.5, -0.5}, {0.5, 0.5, 0.5}, 3), bs_quat(rng),
                           {0.5 * d[0], 0.5 * d[1], 0.5 * d[2]}),
                   "");
        } else if (kind == 3) {
            // Two tori cutting each other: the second turned a quarter about x,
            // its tube crossing the first's twice.
            const double tilt = rng.Uniform(-0.2, 0.2);
            const double h = std::sqrt(0.5);
            bs_run(name + "_tori_crossing", bs_torus(1.0, 0.3, 48, 16),
                   bs_move(
                       bs_torus(1.0, 0.3, 48, 16), {h * std::cos(tilt), h, 0.0, h * std::sin(tilt)},
                       {rng.Uniform(0.45, 0.6), rng.Uniform(-0.1, 0.1), rng.Uniform(-0.1, 0.1)}),
                   "");
        } else {
            const BsVec d = bs_dir(rng);
            const double s = rng.Uniform(0.3, 1.2);
            bs_run(name + "_torus_sphere", bs_torus(1.0, 0.4, 40, 16),
                   mio::icosa_sphere(rng.Uniform(0.4, 0.9), 3, {s * d[0], s * d[1], s * d[2]}), "");
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        std::vector<std::string> args(argv + 1, argv + argc);
        if (args.size() >= 3 && args[0] == "--synthetic") {
            for (std::size_t i = 3; i + 1 < args.size(); i += 2) {
                if (args[i] == "--only")
                    g_only = args[i + 1];
                else if (args[i] == "--write")
                    g_write = args[i + 1];
                else
                    throw std::invalid_argument("unknown option " + args[i]);
            }
            bs_synthetic(std::stoull(args[1]), std::stoi(args[2]));
            return 0;
        }
        if (args.size() >= 2) {
            std::string write, name = args[0] + "|" + args[1];
            for (std::size_t i = 2; i + 1 < args.size(); i += 2) {
                if (args[i] == "--write")
                    write = args[i + 1];
                else if (args[i] == "--case")
                    name = args[i + 1];
                else
                    throw std::invalid_argument("unknown option " + args[i]);
            }
            const mio::Mesh a =
                mio::registry_read(args[0], mio::resolve_format(args[0], ""), mio::ReadOptions{});
            const mio::Mesh b =
                mio::registry_read(args[1], mio::resolve_format(args[1], ""), mio::ReadOptions{});
            bs_run(name, a, b, write);
            return 0;
        }
        std::cerr << "usage: meshioplusplus_boolean_spike A B [--write PREFIX] [--case NAME]\n"
                     "       meshioplusplus_boolean_spike --synthetic SEED N [--only NAME] "
                     "[--write PREFIX]\n";
        return 2;
    } catch (const std::exception& rE) {
        std::cerr << "boolean_spike: " << rE.what() << "\n";
        return 1;
    }
}
