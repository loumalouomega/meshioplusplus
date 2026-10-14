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
#pragma once

// Analytic surfaces for the segment / fit_primitives / step tests: quads on a
// unit cube skin and on surfaces of revolution, built from indices so every
// point lies exactly on its surface (up to rounding).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "mesh_fixtures.hpp"

namespace analytic {

constexpr double kPi = 3.14159265358979323846;

// The closed, outward-wound skin of the unit cube, n x n quads per face.
inline meshioplusplus::Mesh box_skin(int n) {
    std::map<std::tuple<long, long, long>, std::int64_t> ids;
    std::vector<std::vector<double>> pts;
    auto pid = [&](long i, long j, long k) {
        const auto key = std::make_tuple(i, j, k);
        auto it = ids.find(key);
        if (it != ids.end())
            return it->second;
        const std::int64_t id = static_cast<std::int64_t>(pts.size());
        ids[key] = id;
        pts.push_back({double(i) / n, double(j) / n, double(k) / n});
        return id;
    };
    std::vector<std::vector<std::int64_t>> quads;
    for (int axis = 0; axis < 3; ++axis)
        for (int side = 0; side < 2; ++side) {
            const int u = (axis + 1) % 3, v = (axis + 2) % 3;
            for (int a = 0; a < n; ++a)
                for (int b = 0; b < n; ++b) {
                    std::vector<std::int64_t> ring;
                    const int corners[4][2] = {{a, b}, {a + 1, b}, {a + 1, b + 1}, {a, b + 1}};
                    for (const auto& c : corners) {
                        long x[3];
                        x[axis] = side * n;
                        x[u] = c[0];
                        x[v] = c[1];
                        ring.push_back(pid(x[0], x[1], x[2]));
                    }
                    // (u, v, axis) is right-handed, so this ring faces +axis.
                    if (side == 0)
                        ring = {ring[0], ring[3], ring[2], ring[1]};
                    quads.push_back(ring);
                }
        }
    return mt::make_mesh(pts, "quad", quads);
}

// Quads on the surface of revolution about z through the (r, z) profile.
inline meshioplusplus::Mesh revolve(const std::vector<std::pair<double, double>>& rProfile,
                                    int N = 32) {
    std::vector<std::vector<double>> pts;
    for (const auto& rz : rProfile)
        for (int k = 0; k < N; ++k) {
            const double a = 2.0 * kPi * k / N;
            pts.push_back({rz.first * std::cos(a), rz.first * std::sin(a), rz.second});
        }
    std::vector<std::vector<std::int64_t>> quads;
    for (std::size_t i = 0; i + 1 < rProfile.size(); ++i)
        for (int k = 0; k < N; ++k) {
            const std::int64_t k1 = (k + 1) % N;
            const std::int64_t r0 = static_cast<std::int64_t>(i) * N,
                               r1 = static_cast<std::int64_t>(i + 1) * N;
            quads.push_back({r0 + k, r0 + k1, r1 + k1, r1 + k});
        }
    return mt::make_mesh(pts, "quad", quads);
}

inline meshioplusplus::Mesh cylinder(double R = 1.0, double H = 2.0, int M = 6) {
    std::vector<std::pair<double, double>> prof;
    for (int i = 0; i <= M; ++i)
        prof.push_back({R, H * i / M});
    return revolve(prof);
}

inline meshioplusplus::Mesh cone(double HalfAngleDeg = 20.0, int M = 6) {
    const double t = std::tan(HalfAngleDeg * kPi / 180.0);
    std::vector<std::pair<double, double>> prof;
    for (int i = 0; i <= M; ++i) {
        const double z = 1.0 + 2.0 * i / M;
        prof.push_back({t * z, z});
    }
    return revolve(prof);
}

inline meshioplusplus::Mesh sphere_band(double R = 1.5, int M = 12) {
    std::vector<std::pair<double, double>> prof;
    for (int i = 0; i <= M; ++i) {
        const double th = 0.3 + (kPi - 0.6) * i / M;
        prof.push_back({R * std::sin(th), -R * std::cos(th)});
    }
    return revolve(prof);
}

inline meshioplusplus::Mesh torus(double Major = 3.0, double Minor = 1.0, int N = 32, int M = 24) {
    std::vector<std::vector<double>> pts;
    for (int i = 0; i < M; ++i) {
        const double b = 2.0 * kPi * i / M;
        for (int k = 0; k < N; ++k) {
            const double a = 2.0 * kPi * k / N;
            const double rr = Major + Minor * std::cos(b);
            pts.push_back({rr * std::cos(a), rr * std::sin(a), Minor * std::sin(b)});
        }
    }
    std::vector<std::vector<std::int64_t>> quads;
    for (int i = 0; i < M; ++i)
        for (int k = 0; k < N; ++k) {
            const std::int64_t i1 = (i + 1) % M, k1 = (k + 1) % N;
            quads.push_back(
                {std::int64_t(i) * N + k, std::int64_t(i) * N + k1, i1 * N + k1, i1 * N + k});
        }
    return mt::make_mesh(pts, "quad", quads);
}

// A closed surface of revolution about z: the (r, z) profile revolved into
// quads split into triangles, with triangle fans closing each end whose
// radius is not zero at the axis points (0, 0, z_first) and (0, 0, z_last).
// Wound outwards when the profile runs from the bottom pole to the top one.
inline meshioplusplus::Mesh closed_revolve(const std::vector<std::pair<double, double>>& rProfile,
                                           double ZBottom, double ZTop, int N = 32) {
    std::vector<std::vector<double>> pts;
    pts.push_back({0.0, 0.0, ZBottom});
    for (const auto& rz : rProfile)
        for (int k = 0; k < N; ++k) {
            const double a = 2.0 * kPi * k / N;
            pts.push_back({rz.first * std::cos(a), rz.first * std::sin(a), rz.second});
        }
    const std::int64_t top = static_cast<std::int64_t>(pts.size());
    pts.push_back({0.0, 0.0, ZTop});
    auto ring = [&](std::size_t I, int K) {
        return 1 + static_cast<std::int64_t>(I) * N + (K % N);
    };
    const std::size_t last = rProfile.size() - 1;
    std::vector<std::vector<std::int64_t>> tris;
    for (int k = 0; k < N; ++k)
        tris.push_back({0, ring(0, k + 1), ring(0, k)});
    for (int k = 0; k < N; ++k)
        tris.push_back({top, ring(last, k), ring(last, k + 1)});
    for (std::size_t i = 0; i < last; ++i)
        for (int k = 0; k < N; ++k) {
            tris.push_back({ring(i, k), ring(i, k + 1), ring(i + 1, k + 1)});
            tris.push_back({ring(i, k), ring(i + 1, k + 1), ring(i + 1, k)});
        }
    return mt::make_mesh(pts, "triangle", tris);
}

// A capped cylinder of radius R from z = 0 to H: two flat caps and a wall.
inline meshioplusplus::Mesh capped_cylinder(double R = 1.0, double H = 2.0, int M = 6, int N = 32) {
    std::vector<std::pair<double, double>> prof;
    for (int i = 1; i <= 3; ++i)
        prof.push_back({R * i / 3.0, 0.0});
    for (int i = 1; i <= M; ++i)
        prof.push_back({R, H * i / M});
    for (int i = 2; i >= 1; --i)
        prof.push_back({R * i / 3.0, H});
    return closed_revolve(prof, 0.0, H, N);
}

// A capped band of revolution whose radius waves along z, r = 1 + 0.15
// sin(2 pi z) for z in [0, 2], closed by flat caps: no primitive explains the
// wall, and being closed around the axis, neither its best-fit plane nor an
// unrolled cylinder parametrises it.
inline meshioplusplus::Mesh wavy_capped(int N = 96, int M = 64) {
    std::vector<std::pair<double, double>> prof;
    for (int i = 1; i <= 3; ++i)
        prof.push_back({i / 3.0, 0.0});
    for (int i = 1; i <= M; ++i) {
        const double z = 2.0 * i / M;
        prof.push_back({1.0 + 0.15 * std::sin(2.0 * kPi * z), z});
    }
    for (int i = 2; i >= 1; --i)
        prof.push_back({i / 3.0, 2.0});
    return closed_revolve(prof, 0.0, 2.0, N);
}

// A triangle skin under construction: points, and quads split into two
// triangles each wound so its normal points along the quad's expected
// outward direction.
struct SkinBuilder {
    std::vector<std::vector<double>> mPts;
    std::vector<std::vector<std::int64_t>> mTris;

    std::int64_t Add(double X, double Y, double Z) {
        mPts.push_back({X, Y, Z});
        return static_cast<std::int64_t>(mPts.size() - 1);
    }
    // Triangle (a, b, c), flipped when its normal points against `rOut`.
    void Tri(std::int64_t A, std::int64_t B, std::int64_t C, const std::vector<double>& rOut) {
        const auto& p = mPts[static_cast<std::size_t>(A)];
        const auto& q = mPts[static_cast<std::size_t>(B)];
        const auto& r = mPts[static_cast<std::size_t>(C)];
        const double u[3] = {q[0] - p[0], q[1] - p[1], q[2] - p[2]};
        const double v[3] = {r[0] - p[0], r[1] - p[1], r[2] - p[2]};
        const double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                             u[0] * v[1] - u[1] * v[0]};
        if (n[0] * rOut[0] + n[1] * rOut[1] + n[2] * rOut[2] >= 0.0)
            mTris.push_back({A, B, C});
        else
            mTris.push_back({A, C, B});
    }
    // Quad (a, b, c, d) in cyclic order, split along a-c.
    void Quad(std::int64_t A, std::int64_t B, std::int64_t C, std::int64_t D,
              const std::vector<double>& rOut) {
        Tri(A, B, C, rOut);
        Tri(A, C, D, rOut);
    }
    meshioplusplus::Mesh Build() const { return mt::make_mesh(mPts, "triangle", mTris); }
};

// A cylindrical boss of radius R leaning Alpha degrees from the normal of a
// round plate (radius L, thickness T, top at z = 0), capped at axial height
// H, its foot rounded by a concave rolling-ball fillet of radius Rf. The
// fillet's spine -- the plane z = Rf meeting the cylinder of radius R + Rf --
// is an ellipse, so the fillet is a canal surface no primitive explains,
// tangent to a plane and a cylinder that fit_primitives does explain.
inline meshioplusplus::Mesh inclined_boss(double Alpha = 20.0, double R = 0.6, double Rf = 0.2,
                                          double H = 1.5, double L = 1.6, double T = 0.4,
                                          int N = 96, int K = 8) {
    const double al = Alpha * kPi / 180.0, sa = std::sin(al), ca = std::cos(al);
    const double ax[3] = {sa, 0.0, ca}, e1[3] = {ca, 0.0, -sa};
    SkinBuilder b;
    std::vector<std::vector<std::int64_t>> fillet(static_cast<std::size_t>(K + 1));
    std::vector<std::vector<double>> centre(static_cast<std::size_t>(N)), unit(centre.size());
    std::vector<double> foot(static_cast<std::size_t>(N));  // axial height of the cylinder contact
    for (int i = 0; i < N; ++i) {
        const double th = 2.0 * kPi * i / N;
        const double u[3] = {std::cos(th) * e1[0], std::sin(th), std::cos(th) * e1[2]};
        const double h = (Rf + (R + Rf) * sa * std::cos(th)) / ca;
        const double c[3] = {h * ax[0] + (R + Rf) * u[0], h * ax[1] + (R + Rf) * u[1],
                             h * ax[2] + (R + Rf) * u[2]};
        centre[static_cast<std::size_t>(i)] = {c[0], c[1], c[2]};
        unit[static_cast<std::size_t>(i)] = {u[0], u[1], u[2]};
        foot[static_cast<std::size_t>(i)] = h;
        // From the plate contact (-z) to the cylinder contact (-u), the short way.
        const double w = std::acos(std::clamp(u[2], -1.0, 1.0));
        for (int k = 0; k <= K; ++k) {
            const double s = static_cast<double>(k) / K;
            const double f0 = std::sin((1.0 - s) * w) / std::sin(w), f1 = std::sin(s * w) / std::sin(w);
            const double d[3] = {-f1 * u[0], -f1 * u[1], -f0 - f1 * u[2]};
            fillet[static_cast<std::size_t>(k)].push_back(
                b.Add(c[0] + Rf * d[0], c[1] + Rf * d[1], c[2] + Rf * d[2]));
        }
    }
    auto ring_quads = [&](const std::vector<std::int64_t>& rA, const std::vector<std::int64_t>& rB,
                          const auto& rOut) {
        const std::size_t n = rA.size();
        for (std::size_t i = 0; i < n; ++i) {
            const std::size_t j = (i + 1) % n;
            b.Quad(rA[i], rA[j], rB[j], rB[i], rOut(i));
        }
    };
    // The fillet: outward is towards the ball centre.
    for (int k = 0; k < K; ++k)
        ring_quads(fillet[static_cast<std::size_t>(k)], fillet[static_cast<std::size_t>(k + 1)],
                   [&](std::size_t I) {
                       const auto& p = b.mPts[static_cast<std::size_t>(
                           fillet[static_cast<std::size_t>(k)][I])];
                       const auto& c = centre[I];
                       return std::vector<double>{c[0] - p[0], c[1] - p[1], c[2] - p[2]};
                   });
    // The boss wall, from the fillet up to the cap.
    const int mw = 8;
    std::vector<std::int64_t> prev = fillet[static_cast<std::size_t>(K)];
    for (int j = 1; j <= mw; ++j) {
        std::vector<std::int64_t> row;
        for (int i = 0; i < N; ++i) {
            const double t = foot[static_cast<std::size_t>(i)] +
                             (H - foot[static_cast<std::size_t>(i)]) * j / mw;
            const auto& u = unit[static_cast<std::size_t>(i)];
            row.push_back(b.Add(t * ax[0] + R * u[0], t * ax[1] + R * u[1], t * ax[2] + R * u[2]));
        }
        ring_quads(prev, row, [&](std::size_t I) { return unit[I]; });
        prev = row;
    }
    // The cap.
    const std::int64_t top = b.Add(H * ax[0], H * ax[1], H * ax[2]);
    for (int i = 0; i < N; ++i)
        b.Tri(top, prev[static_cast<std::size_t>(i)], prev[static_cast<std::size_t>((i + 1) % N)],
              {ax[0], ax[1], ax[2]});
    // The plate top: rings from the fillet's foot out to the rim, at the
    // foot's own polar angles.
    const int ma = 6;
    std::vector<std::int64_t> inner = fillet[0];
    std::vector<double> psi(static_cast<std::size_t>(N)), rho0(psi.size());
    for (int i = 0; i < N; ++i) {
        const auto& p = b.mPts[static_cast<std::size_t>(inner[static_cast<std::size_t>(i)])];
        psi[static_cast<std::size_t>(i)] = std::atan2(p[1], p[0]);
        rho0[static_cast<std::size_t>(i)] = std::hypot(p[0], p[1]);
    }
    for (int j = 1; j <= ma; ++j) {
        std::vector<std::int64_t> row;
        for (int i = 0; i < N; ++i) {
            const double rr = rho0[static_cast<std::size_t>(i)] +
                              (L - rho0[static_cast<std::size_t>(i)]) * j / ma;
            row.push_back(b.Add(rr * std::cos(psi[static_cast<std::size_t>(i)]),
                                rr * std::sin(psi[static_cast<std::size_t>(i)]), 0.0));
        }
        ring_quads(inner, row, [](std::size_t) { return std::vector<double>{0.0, 0.0, 1.0}; });
        inner = row;
    }
    // The plate's rim and bottom.
    std::vector<std::int64_t> lower;
    for (int i = 0; i < N; ++i)
        lower.push_back(b.Add(L * std::cos(psi[static_cast<std::size_t>(i)]),
                              L * std::sin(psi[static_cast<std::size_t>(i)]), -T));
    ring_quads(inner, lower, [&](std::size_t I) {
        return std::vector<double>{std::cos(psi[I]), std::sin(psi[I]), 0.0};
    });
    const std::int64_t bottom = b.Add(0.0, 0.0, -T);
    for (int i = 0; i < N; ++i)
        b.Tri(bottom, lower[static_cast<std::size_t>(i)],
              lower[static_cast<std::size_t>((i + 1) % N)], {0.0, 0.0, -1.0});
    return b.Build();
}

// A rod about z whose surface is r = f(z - K theta), K = Pitch / (2 pi),
// one period of f given as nodes (rW[i], rR[i]): grid lines follow the
// helices through every node, and each end is closed by a fan to a point on
// the axis Depth pitches inside the rod -- a shallow dish, so a cap never
// meets the flank facing the same way at a dihedral the crease test would
// cross -- wound against the loop's edges in the thread's own triangles
// (part of each loop runs along the axis, where no direction test can tell).
inline meshioplusplus::Mesh helical_rod(const std::vector<double>& rW,
                                        const std::vector<double>& rR, double Pitch, int Turns,
                                        int N, double Depth = 1.0) {
    const double k = Pitch / (2.0 * kPi);
    const int per = static_cast<int>(rW.size());  // samples per period
    const int rows = per * Turns;                 // grid rows j = 0..rows
    SkinBuilder b;
    std::vector<std::vector<std::int64_t>> id(static_cast<std::size_t>(N));
    for (int i = 0; i < N; ++i) {
        const double th = 2.0 * kPi * i / N;
        for (int j = 0; j <= rows; ++j) {
            const double wj = rW[static_cast<std::size_t>(j % per)] + Pitch * (j / per);
            const double rj = rR[static_cast<std::size_t>(j % per)];
            id[static_cast<std::size_t>(i)].push_back(
                b.Add(rj * std::cos(th), rj * std::sin(th), wj + k * th));
        }
    }
    // Node (i, j); column N is column 0 a period up.
    auto node = [&](int I, int J) -> std::int64_t {
        if (I < N)
            return id[static_cast<std::size_t>(I)][static_cast<std::size_t>(J)];
        return J + per <= rows ? id[0][static_cast<std::size_t>(J + per)] : -1;
    };
    auto outward = [&](std::int64_t A) {
        const auto& p = b.mPts[static_cast<std::size_t>(A)];
        return std::vector<double>{p[0], p[1], 0.0};
    };
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < rows; ++j) {
            const std::int64_t a = node(i, j), bb = node(i + 1, j), c = node(i + 1, j + 1),
                               d = node(i, j + 1);
            if (bb < 0 || c < 0)
                continue;
            // The quad's outward direction: radial, and along the axis where
            // the quad is a flank (its rows differ in radius).
            const auto& pa = b.mPts[static_cast<std::size_t>(a)];
            const auto& pd = b.mPts[static_cast<std::size_t>(d)];
            const double ra = std::hypot(pa[0], pa[1]), rd = std::hypot(pd[0], pd[1]);
            std::vector<double> out = outward(a);
            const double rl = std::hypot(out[0], out[1]);
            out = {out[0] / rl, out[1] / rl, rd > ra + 1e-12 ? -1.0 : rd < ra - 1e-12 ? 1.0 : 0.0};
            b.Quad(a, bb, c, d, out);
        }
    auto fan = [&](const std::vector<std::int64_t>& rLoop, bool Reverse) {
        double lo = 1e300, hi = -1e300;
        for (std::int64_t v : rLoop) {
            lo = std::min(lo, b.mPts[static_cast<std::size_t>(v)][2]);
            hi = std::max(hi, b.mPts[static_cast<std::size_t>(v)][2]);
        }
        const std::int64_t c = b.Add(0.0, 0.0, Reverse ? hi + Depth * Pitch : lo - Depth * Pitch);
        for (std::size_t i = 0; i < rLoop.size(); ++i) {
            const std::int64_t p = rLoop[i], q = rLoop[(i + 1) % rLoop.size()];
            if (Reverse)
                b.mTris.push_back({c, q, p});
            else
                b.mTris.push_back({c, p, q});
        }
    };
    std::vector<std::int64_t> low, high;
    for (int i = 0; i < N; ++i)
        low.push_back(node(i, 0));
    for (int j = per; j >= 1; --j)
        low.push_back(node(0, j));
    for (int i = 0; i < N; ++i)
        high.push_back(node(i, rows));
    for (int j = rows - 1; j >= rows - per; --j)
        high.push_back(node(N - 1, j));
    fan(low, true);
    fan(high, false);
    return b.Build();
}

// A rod with a trapezoidal thread about z: root (radius R0), a rising flank,
// crest (radius R1), a falling flank, over Turns turns. The flanks are ruled
// helicoids, the root and crest helical strips of cylinders.
inline meshioplusplus::Mesh threaded_rod(double Pitch = 0.5, double R0 = 0.8, double R1 = 1.0,
                                         int Turns = 4, int N = 64) {
    const double seg[4] = {0.26 * Pitch, 0.24 * Pitch, 0.26 * Pitch, 0.24 * Pitch};
    const int sub[4] = {2, 3, 2, 3};
    std::vector<double> w, r;
    double at = 0.0;
    for (int s = 0; s < 4; ++s) {
        for (int i = 0; i < sub[s]; ++i) {
            const double t = static_cast<double>(i) / sub[s];
            w.push_back(at + seg[s] * t);
            r.push_back(s == 0   ? R0
                        : s == 1 ? R0 + (R1 - R0) * t
                        : s == 2 ? R1
                                 : R1 - (R1 - R0) * t);
        }
        at += seg[s];
    }
    return helical_rod(w, r, Pitch, Turns, N);
}

// A rod with threaded_rod's thread but no cylinder anywhere: the root and
// the crest are circular arcs of radius Rr across the trapezoid's flats, the
// crest bulging out and the root in, meeting the straight flanks at a crease
// (tests/python/analytic_surfaces.py's rounded_thread_rod).
inline meshioplusplus::Mesh rounded_thread_rod(double Pitch = 0.5, double R0 = 0.8, double R1 = 1.0,
                                               double Rr = 0.2, int Turns = 4, int N = 64) {
    const double seg[4] = {0.26 * Pitch, 0.24 * Pitch, 0.26 * Pitch, 0.24 * Pitch};
    const int sub[4] = {4, 3, 4, 3};
    const double half = 0.13 * Pitch, d = std::sqrt(Rr * Rr - half * half);
    std::vector<double> w, r;
    double at = 0.0;
    for (int s = 0; s < 4; ++s) {
        for (int i = 0; i < sub[s]; ++i) {
            const double t = static_cast<double>(i) / sub[s];
            const double x = seg[s] * t - half;  // from the flat's middle
            w.push_back(at + seg[s] * t);
            r.push_back(s == 0   ? R0 + d - std::sqrt(Rr * Rr - x * x)
                        : s == 1 ? R0 + (R1 - R0) * t
                        : s == 2 ? R1 - d + std::sqrt(Rr * Rr - x * x)
                                 : R1 - (R1 - R0) * t);
        }
        at += seg[s];
    }
    return helical_rod(w, r, Pitch, Turns, N);
}

namespace detail_fixture {

using V3 = std::array<double, 3>;
inline V3 v_add(const V3& rA, const V3& rB) {
    return {rA[0] + rB[0], rA[1] + rB[1], rA[2] + rB[2]};
}
inline V3 v_sub(const V3& rA, const V3& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1], rA[2] - rB[2]};
}
inline V3 v_mul(double S, const V3& rA) {
    return {S * rA[0], S * rA[1], S * rA[2]};
}
inline double v_dot(const V3& rA, const V3& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}
inline V3 v_cross(const V3& rA, const V3& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}
inline V3 v_unit(const V3& rA) {
    return v_mul(1.0 / std::sqrt(v_dot(rA, rA)), rA);
}

// Rows of cross-section loops (the same count per row, in order) as quads,
// the ends closed by fans to each end loop's centroid; a loop point with an
// empty outward direction faces away from its row's centre.
inline meshioplusplus::Mesh loft_closed(
    const std::vector<std::vector<std::pair<V3, std::vector<double>>>>& rLoops,
    const std::vector<V3>& rCentres) {
    SkinBuilder b;
    std::vector<std::vector<std::int64_t>> id(rLoops.size());
    for (std::size_t i = 0; i < rLoops.size(); ++i)
        for (const auto& pv : rLoops[i])
            id[i].push_back(b.Add(pv.first[0], pv.first[1], pv.first[2]));
    const std::size_t count = rLoops[0].size();
    for (std::size_t i = 0; i + 1 < rLoops.size(); ++i)
        for (std::size_t j = 0; j < count; ++j) {
            const std::size_t jn = (j + 1) % count;
            std::vector<double> out = rLoops[i][j].second;
            if (out.empty()) {
                const V3 d = v_sub(rLoops[i][j].first, rCentres[i]);
                out = {d[0], d[1], d[2]};
            }
            b.Quad(id[i][j], id[i + 1][j], id[i + 1][jn], id[i][jn], out);
        }
    for (std::size_t e = 0; e < 2; ++e) {
        const std::size_t i = e == 0 ? 0 : rLoops.size() - 1;
        V3 mid{0.0, 0.0, 0.0};
        for (const auto& pv : rLoops[i])
            mid = v_add(mid, pv.first);
        mid = v_mul(1.0 / static_cast<double>(count), mid);
        const std::int64_t c = b.Add(mid[0], mid[1], mid[2]);
        for (std::size_t j = 0; j < count; ++j)
            b.Tri(c, id[i][j], id[i][(j + 1) % count], {e == 0 ? -1.0 : 1.0, 0.0, 0.0});
    }
    return b.Build();
}

}  // namespace detail_fixture

// The block [0, Length] x [0, Width] x [0, Height] whose edge at y = Width,
// z = Height is rounded by a rolling ball of radius r(x) = Rm + Ra cos(pi x /
// Length), varying along it. Each row of the fillet is the exact
// characteristic circle of the ball family, from its foot on the side to its
// foot on the top; r' vanishes at both ends, so the end caps are planes.
inline meshioplusplus::Mesh variable_fillet_block(double Length = 2.0, double Width = 1.0,
                                                  double Height = 1.0, double Rm = 0.25,
                                                  double Ra = 0.05, int N = 48, int K = 16,
                                                  int M = 6) {
    using namespace detail_fixture;
    const double W = Width, H = Height, L = Length;
    std::vector<std::vector<std::pair<V3, std::vector<double>>>> loops;
    std::vector<V3> centres;
    for (int i = 0; i <= N; ++i) {
        const double x = L * i / N;
        const double r = Rm + Ra * std::cos(kPi * x / L);
        const double dr = -Ra * kPi / L * std::sin(kPi * x / L);
        const V3 c{x, W - r, H - r};
        const V3 cd{1.0, -dr, -dr};
        const double cdl = std::sqrt(v_dot(cd, cd));
        const V3 t = v_mul(1.0 / cdl, cd);
        const double delta = -r * dr / cdl;
        const V3 cc = v_add(c, v_mul(delta, t));
        const double rho = std::sqrt(r * r - delta * delta);
        const V3 f_side{x, W, H - r}, f_top{x, W - r, H};
        const V3 a = v_mul(1.0 / rho, v_sub(f_side, cc));
        const V3 b2 = v_sub(f_top, cc);
        const V3 bv = v_unit(v_sub(b2, v_mul(v_dot(b2, a), a)));
        const double om = std::atan2(v_dot(b2, bv), v_dot(b2, a));
        std::vector<std::pair<V3, std::vector<double>>> loop;
        for (int j = 0; j < M; ++j)  // the bottom, y from 0 to W
            loop.push_back({{x, W * j / M, 0.0}, {0.0, 0.0, -1.0}});
        for (int j = 0; j < M; ++j)  // the side, z from 0 to H - r
            loop.push_back({{x, W, (H - r) * j / M}, {0.0, 1.0, 0.0}});
        for (int j = 0; j < K; ++j)  // the fillet, side foot to top foot
            loop.push_back({v_add(cc, v_mul(rho, v_add(v_mul(std::cos(om * j / K), a),
                                                       v_mul(std::sin(om * j / K), bv)))),
                            {}});
        for (int j = 0; j < M; ++j)  // the top, y from W - r to 0
            loop.push_back({{x, (W - r) * (1.0 - static_cast<double>(j) / M), H}, {0.0, 0.0, 1.0}});
        for (int j = 0; j < M; ++j)  // the back, z from H to 0
            loop.push_back({{x, 0.0, H * (1.0 - static_cast<double>(j) / M)}, {0.0, -1.0, 0.0}});
        loops.push_back(std::move(loop));
        centres.push_back(c);
    }
    return loft_closed(loops, centres);
}

// A bar along x whose top and front faces are both freeform, meeting at an
// edge rounded by a ball of radius R rolling along the spine (x, 0, Wave cos(4
// pi x / L)): the top waves with the spine, and the front leans in and out by
// up to Twist degrees, so the angle between the two faces changes along the
// edge. Both faces are ruled surfaces through the ball's feet, their rulings
// normal to its radius there, so tangent to every ball by construction; the
// back, the bottom and the ends are planes.
inline meshioplusplus::Mesh wavy_bar(double Length = 3.0, double Wave = 0.06, double Twist = 3.0,
                                     double R = 0.2, double Back = 1.0, double Depth = 1.0,
                                     int N = 96, int K = 8, int M = 8) {
    using namespace detail_fixture;
    const double L = Length;
    std::vector<std::vector<std::pair<V3, std::vector<double>>>> loops;
    std::vector<V3> centres;
    for (int i = 0; i <= N; ++i) {
        const double x = L * i / N;
        const V3 c{x, 0.0, Wave * std::cos(4.0 * kPi * x / L)};
        const V3 t = v_unit({1.0, 0.0, -Wave * 4.0 * kPi / L * std::sin(4.0 * kPi * x / L)});
        const V3 side{0.0, 1.0, 0.0};
        const V3 up = v_cross(t, side);
        const double psi = Twist * kPi / 180.0 * std::sin(2.0 * kPi * x / L);
        const V3 m1 = up;
        const V3 m2 = v_add(v_mul(std::cos(psi), side), v_mul(std::sin(psi), up));
        const V3 q1 = v_add(c, v_mul(R, m1)), q2 = v_add(c, v_mul(R, m2));
        // The rulings, away from the edge, normal to the ball's radius at the
        // foot and tilted along the spine so they keep their x.
        V3 w1 = v_cross(m1, t);
        w1 = w1[1] < 0.0 ? w1 : v_mul(-1.0, w1);
        w1 = v_sub(w1, v_mul(w1[0] / t[0], t));
        V3 w2 = v_cross(m2, t);
        w2 = w2[2] < 0.0 ? w2 : v_mul(-1.0, w2);
        w2 = v_sub(w2, v_mul(w2[0] / t[0], t));
        const V3 top_far = v_add(q1, v_mul((q1[1] + Back) / -w1[1], w1));
        const V3 front_low = v_add(q2, v_mul((q2[2] + Depth) / -w2[2], w2));
        const V3 corner{top_far[0], -Back, -Depth};
        const double om = std::acos(std::clamp(v_dot(m1, m2), -1.0, 1.0));
        std::vector<std::pair<V3, std::vector<double>>> loop;
        auto lerp = [](const V3& rA, const V3& rB, double F) {
            return v_add(rA, v_mul(F, v_sub(rB, rA)));
        };
        for (int j = 0; j < M; ++j)  // the back, up to the top's far edge
            loop.push_back({lerp(corner, top_far, static_cast<double>(j) / M), {0.0, -1.0, 0.0}});
        for (int j = 0; j < M; ++j)  // the top, towards the edge
            loop.push_back({lerp(top_far, q1, static_cast<double>(j) / M), {m1[0], m1[1], m1[2]}});
        for (int j = 0; j < K; ++j) {  // the blend, top foot to front foot
            const double f0 = std::sin((1.0 - static_cast<double>(j) / K) * om) / std::sin(om);
            const double f1 = std::sin(static_cast<double>(j) / K * om) / std::sin(om);
            loop.push_back({v_add(c, v_mul(R, v_add(v_mul(f0, m1), v_mul(f1, m2)))), {}});
        }
        for (int j = 0; j < M; ++j)  // the front, down to the bottom
            loop.push_back(
                {lerp(q2, front_low, static_cast<double>(j) / M), {m2[0], m2[1], m2[2]}});
        for (int j = 0; j < M; ++j)  // the bottom, back to the corner
            loop.push_back({lerp(front_low, corner, static_cast<double>(j) / M), {0.0, 0.0, -1.0}});
        loops.push_back(std::move(loop));
        centres.push_back(c);
    }
    return loft_closed(loops, centres);
}

// The turned part of tests/python/analytic_surfaces.py: a flat base, a
// cylinder wall, a toroidal fillet tangent to the wall and to a flat
// shoulder, a cone and a spherical dome, closed by fans at the poles.
inline meshioplusplus::Mesh turned_part(int N = 64) {
    std::vector<std::pair<double, double>> prof;
    for (int i = 1; i <= 4; ++i)
        prof.push_back({2.0 * i / 4, 0.0});
    for (int i = 1; i <= 8; ++i)
        prof.push_back({2.0, 1.5 * i / 8});
    for (int i = 1; i <= 8; ++i) {
        const double a = (kPi / 2) * i / 8;
        prof.push_back({1.5 + 0.5 * std::cos(a), 1.5 + 0.5 * std::sin(a)});
    }
    for (int i = 1; i <= 3; ++i)
        prof.push_back({1.5 - 0.5 * i / 3, 2.0});
    for (int i = 1; i <= 6; ++i)
        prof.push_back({1.0 - 0.4 * i / 6, 2.0 + 1.0 * i / 6});
    const double zc = 1.625, big_r = 1.5, a0 = std::asin(0.6 / big_r);
    for (int i = 1; i <= 5; ++i) {
        const double a = a0 * (1.0 - i / 6.0);
        prof.push_back({big_r * std::sin(a), zc + big_r * std::cos(a)});
    }
    return closed_revolve(prof, 0.0, zc + big_r, N);
}

}  // namespace analytic
