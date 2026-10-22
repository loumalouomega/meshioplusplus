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

#include "meshioplusplus/detail/data_ops.hpp"
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
//
// Taper and PitchRate generalise the screw (roadmap 6.3.3): a node at the
// unwrapped angle T (its column's angle plus a turn per period up) sits at
// z = w + g(T), g(T) = K T + PitchRate T^2 / (8 pi^2) -- the pitch grows by
// PitchRate per turn -- and r = f + Taper (w + g(T)), the profile cut on a
// cone whose radius grows by Taper per unit length. Both 0 is the rod above.
inline meshioplusplus::Mesh helical_rod(const std::vector<double>& rW,
                                        const std::vector<double>& rR, double Pitch, int Turns,
                                        int N, double Depth = 1.0, double Taper = 0.0,
                                        double PitchRate = 0.0) {
    const double k = Pitch / (2.0 * kPi);
    const int per = static_cast<int>(rW.size());  // samples per period
    const int rows = per * Turns;                 // grid rows j = 0..rows
    SkinBuilder b;
    std::vector<std::vector<std::int64_t>> id(static_cast<std::size_t>(N));
    for (int i = 0; i < N; ++i) {
        const double th = 2.0 * kPi * i / N;
        for (int j = 0; j <= rows; ++j) {
            if (Taper == 0.0 && PitchRate == 0.0) {
                const double wj = rW[static_cast<std::size_t>(j % per)] + Pitch * (j / per);
                const double rj = rR[static_cast<std::size_t>(j % per)];
                id[static_cast<std::size_t>(i)].push_back(
                    b.Add(rj * std::cos(th), rj * std::sin(th), wj + k * th));
                continue;
            }
            const double big = th + 2.0 * kPi * (j / per);
            const double g = k * big + PitchRate * big * big / (8.0 * kPi * kPi);
            const double w = rW[static_cast<std::size_t>(j % per)];
            const double rj = rR[static_cast<std::size_t>(j % per)] + Taper * (w + g);
            id[static_cast<std::size_t>(i)].push_back(
                b.Add(rj * std::cos(th), rj * std::sin(th), w + g));
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

// threaded_rod's trapezoidal thread cut on a cone: its radius grows by Taper
// per unit length (a pipe thread), so the root and the crest are helical
// strips of cones and the flanks ruled helicoids (roadmap 6.3.3).
inline meshioplusplus::Mesh tapered_thread_rod(double Pitch = 0.5, double R0 = 0.8, double R1 = 1.0,
                                               double Taper = 0.0625, int Turns = 4, int N = 64) {
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
    return helical_rod(w, r, Pitch, Turns, N, 1.0, Taper);
}

// rounded_thread_rod's thread cut on a cone, as tapered_thread_rod's: no
// cylinder or cone anywhere, so the screw, its taper and the profile are all
// found from the flanks themselves.
inline meshioplusplus::Mesh tapered_rounded_thread_rod(double Pitch = 0.5, double R0 = 0.8,
                                                       double R1 = 1.0, double Rr = 0.2,
                                                       double Taper = 0.0625, int Turns = 4,
                                                       int N = 64) {
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
    return helical_rod(w, r, Pitch, Turns, N, 1.0, Taper);
}

// threaded_rod's trapezoidal thread with a pitch growing by PitchRate per
// turn (roadmap 6.3.3). A period runs from the middle of one root to the
// middle of the next, so the root -- a helical strip of the cylinder R0 --
// takes up the growing gap and the flanks stay ruled helicoids, whose
// advance is quadratic in the angle.
inline meshioplusplus::Mesh variable_pitch_thread_rod(double Pitch = 0.5, double R0 = 0.8,
                                                      double R1 = 1.0, double PitchRate = 0.05,
                                                      int Turns = 4, int N = 64) {
    const double seg[5] = {0.13 * Pitch, 0.24 * Pitch, 0.26 * Pitch, 0.24 * Pitch, 0.13 * Pitch};
    const int sub[5] = {1, 3, 2, 3, 1};
    std::vector<double> w, r;
    double at = 0.0;
    for (int s = 0; s < 5; ++s) {
        for (int i = 0; i < sub[s]; ++i) {
            const double t = static_cast<double>(i) / sub[s];
            w.push_back(at + seg[s] * t);
            r.push_back(s == 0 || s == 4 ? R0
                        : s == 1         ? R0 + (R1 - R0) * t
                        : s == 2         ? R1
                                         : R1 - (R1 - R0) * t);
        }
        at += seg[s];
    }
    return helical_rod(w, r, Pitch, Turns, N, 1.0, 0.0, PitchRate);
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

// The block [0, Length] x [0, Width] x [0, Height] whose edge at y = Width,
// z = Height is cut by a blend from the top, D1 in from the edge, to the
// side, d2(x) down from it, d2 growing linearly from D2a to D2b along the
// edge (roadmap 6.3.2). Each cross-section, in the plane x = const, is the
// rational quadratic from the top foot through the corner (x, Width,
// Height) to the side foot with weight Weight: 0 a chord -- a chamfer of
// varying width, a hyperbolic paraboloid -- and more a conic arc tangent to
// both faces, a circle only where the legs are equal and Weight is
// sin(pi/4). The triangles carry their faces as `segment:patch` (0 bottom, 1
// side, 2 top, 3 back, 4 and 5 the ends, 6 the blend): a conic blend is
// tangent to both faces and its curvature varies across it, which no
// segmentation by curvature isolates.
inline meshioplusplus::Mesh section_blend_block(double Weight = 0.0, double Length = 2.0,
                                                double Width = 1.0, double Height = 1.0,
                                                double D1 = 0.3, double D2a = 0.2, double D2b = 0.4,
                                                int N = 40, int K = 12, int M = 6) {
    using namespace detail_fixture;
    const double W = Width, H = Height, L = Length;
    std::vector<std::vector<std::pair<V3, std::vector<double>>>> loops;
    std::vector<V3> centres;
    for (int i = 0; i <= N; ++i) {
        const double x = L * i / N;
        const double d2 = D2a + (D2b - D2a) * x / L;
        const V3 top{x, W - D1, H}, side{x, W, H - d2}, apex{x, W, H};
        const V3 c{x, W - D1, H - d2};  // inside the blend: its points face away
        std::vector<std::pair<V3, std::vector<double>>> loop;
        for (int j = 0; j < M; ++j)  // the bottom, y from 0 to W
            loop.push_back({{x, W * j / M, 0.0}, {0.0, 0.0, -1.0}});
        for (int j = 0; j < M; ++j)  // the side, z from 0 to H - d2
            loop.push_back({{x, W, (H - d2) * j / M}, {0.0, 1.0, 0.0}});
        for (int j = 0; j < K; ++j) {  // the blend, side foot to top foot
            // Evenly along a chamfer's chord; by the conic's own parameter else.
            const double u = static_cast<double>(j) / K;
            const double b0 = Weight > 0.0 ? (1.0 - u) * (1.0 - u) : 1.0 - u,
                         b1 = 2.0 * Weight * u * (1.0 - u), b2 = Weight > 0.0 ? u * u : u,
                         den = b0 + b1 + b2;
            loop.push_back(
                {v_mul(1.0 / den, v_add(v_add(v_mul(b0, side), v_mul(b1, apex)), v_mul(b2, top))),
                 {}});
        }
        for (int j = 0; j < M; ++j)  // the top, y from W - D1 to 0
            loop.push_back(
                {{x, (W - D1) * (1.0 - static_cast<double>(j) / M), H}, {0.0, 0.0, 1.0}});
        for (int j = 0; j < M; ++j)  // the back, z from H to 0
            loop.push_back({{x, 0.0, H * (1.0 - static_cast<double>(j) / M)}, {0.0, -1.0, 0.0}});
        loops.push_back(std::move(loop));
        centres.push_back(c);
    }
    meshioplusplus::Mesh mesh = loft_closed(loops, centres);
    const double* p = mesh.Points().As<double>();
    const std::int64_t* conn = mesh.Cells(0).Conn().As<std::int64_t>();
    const std::size_t nc = mesh.Cells(0).NumCells();
    meshioplusplus::NDArray lab =
        meshioplusplus::NDArray::Uninit(meshioplusplus::DType::Int64, {nc});
    for (std::size_t t = 0; t < nc; ++t) {
        double m[3] = {0.0, 0.0, 0.0};
        for (std::size_t k = 0; k < 3; ++k)
            for (std::size_t d = 0; d < 3; ++d)
                m[d] += p[static_cast<std::size_t>(conn[t * 3 + k]) * 3 + d] / 3.0;
        const double e = 1e-9;
        lab.As<std::int64_t>()[t] = std::fabs(m[2]) < e       ? 0
                                    : std::fabs(m[1] - W) < e ? 1
                                    : std::fabs(m[2] - H) < e ? 2
                                    : std::fabs(m[1]) < e     ? 3
                                    : std::fabs(m[0]) < e     ? 4
                                    : std::fabs(m[0] - L) < e ? 5
                                                              : 6;
    }
    std::vector<meshioplusplus::NDArray> blocks;
    blocks.push_back(std::move(lab));
    mesh.AddCellData("segment:patch", std::move(blocks));
    return mesh;
}

// The box [-Sx, Sx] x [-Sy, Sy] x [-Sz, Sz] with every edge and corner
// rounded: the cube [-1, 1]^3 rounded by a ball of radius R (each point the
// nearest point of the inner cube [R - 1, 1 - R]^3 plus R along the way out),
// scaled by (Sx, Sy, Sz). Scaled unequally, its 12 edge blends are conic
// blends (elliptic cross-sections tangent to both faces) and its 8 corners
// ellipsoid patches, each a vertex blend tangent to its three edge blends
// (roadmap 6.3.1); the faces stay planes. Every face of the cube is gridded
// with lines on the blends' boundaries (Nf cells across a flat, Nb across
// half a blend, evenly in angle), so no triangle straddles two patches, and
// the triangles carry the patches as `segment:patch`: no segmentation by
// curvature can isolate a tangent blend whose curvature varies across it.
inline meshioplusplus::Mesh rounded_box(double R = 0.3, double Sx = 1.0, double Sy = 1.3,
                                        double Sz = 0.8, int Nf = 4, int Nb = 6) {
    const double inner = 1.0 - R;
    std::vector<double> g;  // the grid's coordinates along one edge of a face
    for (int k = Nb; k >= 1; --k)
        g.push_back(-(inner + R * std::tan(0.25 * kPi * k / Nb)));
    for (int k = 0; k <= Nf; ++k)
        g.push_back(-inner + 2.0 * inner * k / Nf);
    for (int k = 1; k <= Nb; ++k)
        g.push_back(inner + R * std::tan(0.25 * kPi * k / Nb));
    const double scale[3] = {Sx, Sy, Sz};
    SkinBuilder b;
    std::map<std::array<long long, 3>, std::int64_t> ids;
    auto vertex = [&](const std::array<double, 3>& rX) {
        const std::array<long long, 3> key{std::llround(rX[0] * 1e9), std::llround(rX[1] * 1e9),
                                           std::llround(rX[2] * 1e9)};
        auto it = ids.find(key);
        if (it != ids.end())
            return it->second;
        double c[3], d[3], len = 0.0;
        for (int k = 0; k < 3; ++k) {
            c[k] = std::clamp(rX[static_cast<std::size_t>(k)], -inner, inner);
            d[k] = rX[static_cast<std::size_t>(k)] - c[k];
            len += d[k] * d[k];
        }
        len = std::sqrt(len);
        const std::int64_t id =
            b.Add(scale[0] * (c[0] + R * d[0] / len), scale[1] * (c[1] + R * d[1] / len),
                  scale[2] * (c[2] + R * d[2] / len));
        ids.emplace(key, id);
        return id;
    };
    // A triangle's patch: the signs of the inner cube's faces its centroid
    // lies beyond (one: a face, two: an edge blend, three: a corner).
    std::vector<std::int64_t> label;
    auto patch_of = [&](const std::array<double, 3>& rC) {
        std::int64_t code = 0;
        for (int k = 0; k < 3; ++k) {
            const double x = rC[static_cast<std::size_t>(k)];
            code = code * 3 + (x > inner ? 2 : x < -inner ? 0 : 1);
        }
        return code;
    };
    const int n = static_cast<int>(g.size());
    for (int a = 0; a < 3; ++a)
        for (int sign : {-1, 1}) {
            const int u = (a + 1) % 3, v = (a + 2) % 3;
            auto at = [&](int I, int J) {
                std::array<double, 3> x{};
                x[static_cast<std::size_t>(a)] = sign;
                x[static_cast<std::size_t>(u)] = g[static_cast<std::size_t>(I)];
                x[static_cast<std::size_t>(v)] = g[static_cast<std::size_t>(J)];
                return x;
            };
            std::vector<double> out(3, 0.0);
            out[static_cast<std::size_t>(a)] = sign;
            for (int i = 0; i + 1 < n; ++i)
                for (int j = 0; j + 1 < n; ++j) {
                    std::array<double, 3> c{};
                    for (const auto& x : {at(i, j), at(i + 1, j + 1)})
                        for (int k = 0; k < 3; ++k)
                            c[static_cast<std::size_t>(k)] += 0.5 * x[static_cast<std::size_t>(k)];
                    b.Quad(vertex(at(i, j)), vertex(at(i + 1, j)), vertex(at(i + 1, j + 1)),
                           vertex(at(i, j + 1)), out);
                    label.push_back(patch_of(c));
                    label.push_back(patch_of(c));
                }
        }
    meshioplusplus::Mesh mesh = b.Build();
    meshioplusplus::NDArray lab =
        meshioplusplus::NDArray::Uninit(meshioplusplus::DType::Int64, {label.size()});
    std::copy(label.begin(), label.end(), lab.As<std::int64_t>());
    std::vector<meshioplusplus::NDArray> blocks;
    blocks.push_back(std::move(lab));
    mesh.AddCellData("segment:patch", std::move(blocks));
    return mesh;
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

// The plate [-a, a]^2 x [0, h] of tests/python/analytic_surfaces.py, with a
// through-hole of radius r at its centre and its corner (a, a) rounded by a
// quarter circle of radius rc tangent to both sides: a CAD part tessellated as
// CAD exporters do, rings of vertices on the hole and the outline (k segments
// along the fillet) joined by strips, every vertex on its surfaces. Six planes
// and two cylinders; the fillet meets its sides at smooth edges.
inline meshioplusplus::Mesh filleted_plate(double A = 2.0, double H = 0.5, double R = 0.6,
                                           double Rc = 0.8, int N = 48, int K = 8) {
    const double c = A - Rc;
    std::vector<double> angles;
    for (int i = 0; i < N; ++i)
        angles.push_back(2.0 * kPi * i / N);
    for (int j = 0; j < 3; ++j)
        angles.push_back(kPi * (0.75 + 0.5 * j));  // the sharp corners
    for (int j = 0; j <= K; ++j) {  // along the fillet, from tangent point to tangent point
        const double phi = 0.5 * kPi * j / K;
        angles.push_back(std::atan2(c + Rc * std::sin(phi), c + Rc * std::cos(phi)));
    }
    std::sort(angles.begin(), angles.end());
    angles.erase(std::unique(angles.begin(), angles.end()), angles.end());
    const double lo = std::atan2(c, A), hi = std::atan2(A, c);
    SkinBuilder b;
    std::vector<std::array<std::int64_t, 4>> rings;
    for (double th : angles) {
        const double d[2] = {std::cos(th), std::sin(th)};
        double t = A / std::max(std::fabs(d[0]), std::fabs(d[1]));
        if (lo <= th && th <= hi) {
            const double dc = d[0] * c + d[1] * c;
            t = dc + std::sqrt(dc * dc - 2.0 * c * c + Rc * Rc);
        }
        rings.push_back({b.Add(R * d[0], R * d[1], 0.0), b.Add(R * d[0], R * d[1], H),
                         b.Add(t * d[0], t * d[1], 0.0), b.Add(t * d[0], t * d[1], H)});
    }
    const std::size_t m = rings.size();
    for (std::size_t i = 0; i < m; ++i) {
        const auto& r0 = rings[i];
        const auto& r1 = rings[(i + 1) % m];
        const double mid =
            0.5 * (angles[i] + angles[(i + 1) % m] + (i + 1 == m ? 2.0 * kPi : 0.0));
        const std::vector<double> out{std::cos(mid), std::sin(mid), 0.0};
        const std::vector<double> in{-out[0], -out[1], 0.0};
        b.Quad(r0[1], r0[3], r1[3], r1[1], {0.0, 0.0, 1.0});
        b.Quad(r0[0], r0[2], r1[2], r1[0], {0.0, 0.0, -1.0});
        b.Quad(r0[2], r1[2], r1[3], r0[3], out);
        b.Quad(r0[0], r1[0], r1[1], r0[1], in);
    }
    return b.Build();
}

// ---------------------------------------------------------------------------
// Freeform fixtures (v17.12.0): skins whose faces no primitive explains, for
// fit_cad's freeform modes. Built from indices and closed-form functions only,
// so `tests/python/analytic_surfaces.py` builds the same shapes; the
// trig-free ones (the monkey saddle, the S band, the hood, the Catmull-Clark
// cages and the genus-2 implicit) bit for bit, since NumPy may dispatch
// `sin`/`cos`/`exp` to SIMD kernels whose last bits differ.
// ---------------------------------------------------------------------------

// The unit cube's skin as triangles, n x n quads per face split along a-c,
// with `Warp` a monotone reparametrisation of x and y that leaves every face
// planar but grades the triangles (anisotropic near the middle).
inline meshioplusplus::Mesh triangle_box(int n, bool Warp = false) {
    const meshioplusplus::Mesh quads = box_skin(n);
    const meshioplusplus::NDArray& p = quads.Points();
    SkinBuilder b;
    for (std::size_t i = 0; i < quads.NumPoints(); ++i) {
        double x = p.As<double>()[i * 3], y = p.As<double>()[i * 3 + 1];
        if (Warp) {
            x += 0.35 * x * (1.0 - x) * (1.0 - 2.0 * x);
            y += 0.35 * y * (1.0 - y) * (1.0 - 2.0 * y);
        }
        b.Add(x, y, p.As<double>()[i * 3 + 2]);
    }
    const meshioplusplus::NDArray& conn = quads.Cells(0).Conn();
    for (std::size_t q = 0; q < quads.Cells(0).NumCells(); ++q) {
        const std::int64_t* r = conn.As<std::int64_t>() + q * 4;
        b.mTris.push_back({r[0], r[1], r[2]});
        b.mTris.push_back({r[0], r[2], r[3]});
    }
    return b.Build();
}

// A closed box whose top (z = 1) is lifted to z = 1 + F(x, y). The sides'
// points follow in proportion (z (1 + F)), only in z, so the four sides and
// the bottom stay planar and a top lowered by more than a row's height does
// not fold its sides' top row over (v17.12.0; a point of the top and of a
// side where F = 0 is where it was).
template <class Height>
inline meshioplusplus::Mesh heightfield_box(int n, Height F, bool Warp = false) {
    meshioplusplus::Mesh m = triangle_box(n, Warp);
    meshioplusplus::NDArray q =
        meshioplusplus::NDArray::Uninit(meshioplusplus::DType::Float64, {m.NumPoints(), 3});
    double* a = q.As<double>();
    std::copy(m.Points().As<double>(), m.Points().As<double>() + m.NumPoints() * 3, a);
    for (std::size_t i = 0; i < m.NumPoints(); ++i)
        if (a[i * 3 + 2] > 0.0)
            a[i * 3 + 2] = a[i * 3 + 2] * (1.0 + F(a[i * 3], a[i * 3 + 1]));
    m.AssignPoints(std::move(q));
    return m;
}

// The monkey saddle z = x^3 - 3 x y^2 over the top (trig-free).
inline meshioplusplus::Mesh monkey_saddle_box(int n = 24, bool Warp = false) {
    return heightfield_box(
        n,
        [](double X, double Y) {
            const double u = 2.0 * X - 1.0, v = 2.0 * Y - 1.0;
            return 0.15 * (u * u * u - 3.0 * u * v * v);
        },
        Warp);
}

// A Gaussian bump on a plate.
inline meshioplusplus::Mesh gaussian_bump_box(int n = 24, bool Warp = false) {
    return heightfield_box(
        n,
        [](double X, double Y) {
            const double r2 = (X - 0.5) * (X - 0.5) + (Y - 0.5) * (Y - 0.5);
            return 0.3 * std::exp(-r2 / 0.03);
        },
        Warp);
}

// An S band: a smootherstep rise across x, sheared along y (trig-free).
inline meshioplusplus::Mesh s_band_box(int n = 24, bool Warp = false) {
    return heightfield_box(
        n,
        [](double X, double Y) {
            const double s = X * X * X * (X * (6.0 * X - 15.0) + 10.0);
            return 0.25 * (s - 0.5) * (1.0 + 0.5 * (Y - 0.5));
        },
        Warp);
}

// A car-hood-like freeform top closed by planes (trig-free).
inline meshioplusplus::Mesh hood_box(int n = 24, bool Warp = false) {
    return heightfield_box(
        n,
        [](double X, double Y) {
            return 0.35 * 16.0 * X * (1.0 - X) * Y * (1.0 - Y) * (1.0 + 0.4 * X);
        },
        Warp);
}

// A genus-1 organic blob: a torus whose tube radius ripples around it.
inline meshioplusplus::Mesh organic_torus(int N = 48, int M = 24) {
    SkinBuilder b;
    for (int i = 0; i < M; ++i) {
        const double t = 2.0 * kPi * i / M;
        for (int k = 0; k < N; ++k) {
            const double a = 2.0 * kPi * k / N;
            const double r = 0.35 * (1.0 + 0.2 * std::cos(2.0 * t) * std::sin(3.0 * a));
            const double rr = 1.0 + r * std::cos(t);
            b.Add(rr * std::cos(a), rr * std::sin(a), 0.7 * r * std::sin(t));
        }
    }
    for (int i = 0; i < M; ++i)
        for (int k = 0; k < N; ++k) {
            const std::int64_t i1 = (i + 1) % M, k1 = (k + 1) % N;
            const std::int64_t a = std::int64_t(i) * N + k, bb = std::int64_t(i) * N + k1;
            const std::int64_t c = i1 * N + k1, d = i1 * N + k;
            // (a, b, c, d) winds outward, as `torus` does.
            b.mTris.push_back({a, bb, c});
            b.mTris.push_back({a, c, d});
        }
    return b.Build();
}

// The zero level of `F` (negative inside) over an n^3 grid of the box
// [Lo, Hi] (chosen so no node lands on the level set: a crossing at a node
// would make a triangle of no area, and no orientation), each cube cut into the six Kuhn tetrahedra
// around its main diagonal, by marching tetrahedra: a crossing on each grid edge the level set cuts
// (shared by every tetrahedron on that edge, so the skin is closed), a node at exactly 0 counted
// outside, no triangle dropped, each wound along the gradient of the field's linear interpolant
// over its tetrahedron.
template <class Field>
inline meshioplusplus::Mesh implicit_skin(Field F, std::array<double, 3> Lo,
                                          std::array<double, 3> Hi, std::int64_t n) {
    const std::int64_t m = n + 1;
    const auto node = [&](std::int64_t i, std::int64_t j, std::int64_t k) {
        return (k * m + j) * m + i;
    };
    std::vector<std::array<double, 3>> xyz(static_cast<std::size_t>(m * m * m));
    std::vector<double> f(xyz.size());
    for (std::int64_t k = 0; k < m; ++k)
        for (std::int64_t j = 0; j < m; ++j)
            for (std::int64_t i = 0; i < m; ++i) {
                const std::size_t id = static_cast<std::size_t>(node(i, j, k));
                xyz[id] = {
                    Lo[0] + (Hi[0] - Lo[0]) * static_cast<double>(i) / static_cast<double>(n),
                    Lo[1] + (Hi[1] - Lo[1]) * static_cast<double>(j) / static_cast<double>(n),
                    Lo[2] + (Hi[2] - Lo[2]) * static_cast<double>(k) / static_cast<double>(n)};
                f[id] = F(xyz[id][0], xyz[id][1], xyz[id][2]);
            }
    SkinBuilder b;
    std::map<std::pair<std::int64_t, std::int64_t>, std::int64_t> crossing;
    // The level set never needs the quad's own diagonal to be consistent:
    // both halves take the same outward direction.
    const auto cross = [&](std::int64_t A, std::int64_t B) {
        const auto key = std::make_pair(std::min(A, B), std::max(A, B));
        auto it = crossing.find(key);
        if (it != crossing.end())
            return it->second;
        // From the inside node (negative) to the outside one.
        const std::int64_t in =
            f[static_cast<std::size_t>(key.first)] < 0.0 ? key.first : key.second;
        const std::int64_t out = in == key.first ? key.second : key.first;
        const double fi = f[static_cast<std::size_t>(in)], fo = f[static_cast<std::size_t>(out)];
        const double t = fi / (fi - fo);
        const auto& p = xyz[static_cast<std::size_t>(in)];
        const auto& q = xyz[static_cast<std::size_t>(out)];
        const std::int64_t id =
            b.Add(p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]), p[2] + t * (q[2] - p[2]));
        crossing.emplace(key, id);
        return id;
    };
    // The six Kuhn tetrahedra: (0, e_a, e_a + e_b, 1) for each axis order.
    const int orders[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
    for (std::int64_t k = 0; k < n; ++k)
        for (std::int64_t j = 0; j < n; ++j)
            for (std::int64_t i = 0; i < n; ++i)
                for (const auto& o : orders) {
                    std::int64_t c[3] = {i, j, k};
                    std::int64_t v[4];
                    v[0] = node(c[0], c[1], c[2]);
                    ++c[o[0]];
                    v[1] = node(c[0], c[1], c[2]);
                    ++c[o[1]];
                    v[2] = node(c[0], c[1], c[2]);
                    ++c[o[2]];
                    v[3] = node(c[0], c[1], c[2]);
                    std::vector<std::int64_t> ins, outs;
                    for (std::int64_t x : v)
                        (f[static_cast<std::size_t>(x)] < 0.0 ? ins : outs).push_back(x);
                    if (ins.empty() || outs.empty())
                        continue;
                    // Outward: the gradient of the field's linear interpolant
                    // over the tetrahedron, whose level set the cut is.
                    std::vector<double> out(3, 0.0);
                    for (int a = 0; a < 3; ++a)
                        out[static_cast<std::size_t>(o[a])] =
                            (f[static_cast<std::size_t>(v[a + 1])] -
                             f[static_cast<std::size_t>(v[a])]) /
                            (Hi[static_cast<std::size_t>(o[a])] -
                             Lo[static_cast<std::size_t>(o[a])]);
                    if (ins.size() == 1 || outs.size() == 1) {
                        const bool one_in = ins.size() == 1;
                        const std::int64_t apex = one_in ? ins[0] : outs[0];
                        const std::vector<std::int64_t>& rest = one_in ? outs : ins;
                        b.Tri(cross(apex, rest[0]), cross(apex, rest[1]), cross(apex, rest[2]),
                              out);
                    } else {
                        // A quad: in0-out0, in0-out1, in1-out1, in1-out0.
                        b.Quad(cross(ins[0], outs[0]), cross(ins[0], outs[1]),
                               cross(ins[1], outs[1]), cross(ins[1], outs[0]), out);
                    }
                }
    return b.Build();
}

// A sphere with a handle (genus 1): a sphere and a torus through it, joined by
// a polynomial smooth minimum.
inline meshioplusplus::Mesh sphere_with_handle(std::int64_t n = 32) {
    return implicit_skin(
        [](double X, double Y, double Z) {
            const double sphere = std::sqrt(X * X + Y * Y + Z * Z) - 0.6;
            const double q = std::sqrt((X - 0.6) * (X - 0.6) + Z * Z) - 0.42;
            const double ring = std::sqrt(q * q + Y * Y) - 0.11;
            const double k = 0.08, h = std::max(k - std::fabs(sphere - ring), 0.0) / k;
            return std::min(sphere, ring) - h * h * k * 0.25;
        },
        {-0.8123, -0.8077, -0.7951}, {1.2111, 0.8173, 0.8029}, n);
}

// A double torus (genus 2): a tube around the lemniscate (trig-free).
inline meshioplusplus::Mesh double_torus(std::int64_t n = 40) {
    return implicit_skin(
        [](double X, double Y, double Z) {
            const double r2 = X * X + Y * Y;
            const double l = r2 * r2 - X * X + Y * Y;
            return l * l + 4.0 * Z * Z - 0.02;
        },
        {-1.3117, -0.6531, -0.4073}, {1.3089, 0.6563, 0.4049}, n);
}

// Catmull-Clark subdivision of a closed, outward-wound polygon cage, `Levels`
// steps, every point then moved to its limit position (the limit stencil of
// a quad mesh), and the quads split along a-c: an exact sampling of the limit
// surface, with an extraordinary vertex for every cage vertex of valence
// other than 4 and every cage face of other than 4 sides.
inline meshioplusplus::Mesh catmull_clark_limit(std::vector<std::array<double, 3>> P,
                                                std::vector<std::vector<std::int64_t>> F,
                                                int Levels) {
    using Edge = std::pair<std::int64_t, std::int64_t>;
    for (int level = 0; level < Levels; ++level) {
        const std::size_t nv = P.size();
        std::vector<std::array<double, 3>> fp(F.size(), {0.0, 0.0, 0.0});
        for (std::size_t f = 0; f < F.size(); ++f) {
            for (std::int64_t v : F[f])
                for (int j = 0; j < 3; ++j)
                    fp[f][j] += P[static_cast<std::size_t>(v)][j];
            for (int j = 0; j < 3; ++j)
                fp[f][j] /= static_cast<double>(F[f].size());
        }
        std::map<Edge, std::vector<std::size_t>> edges;  // sorted: a fixed numbering
        for (std::size_t f = 0; f < F.size(); ++f)
            for (std::size_t k = 0; k < F[f].size(); ++k) {
                const std::int64_t a = F[f][k], b = F[f][(k + 1) % F[f].size()];
                edges[{std::min(a, b), std::max(a, b)}].push_back(f);
            }
        std::map<Edge, std::int64_t> edge_id;
        std::vector<std::array<double, 3>> next(P.size() + F.size() + edges.size());
        for (std::size_t f = 0; f < F.size(); ++f)
            next[nv + f] = fp[f];
        std::int64_t e = static_cast<std::int64_t>(nv + F.size());
        std::vector<std::array<double, 3>> fsum(nv, {0.0, 0.0, 0.0}), rsum(nv, {0.0, 0.0, 0.0});
        std::vector<int> valence(nv, 0), nfaces(nv, 0);
        for (const auto& kv : edges) {
            const auto& a = P[static_cast<std::size_t>(kv.first.first)];
            const auto& b = P[static_cast<std::size_t>(kv.first.second)];
            std::array<double, 3> ep{};
            for (int j = 0; j < 3; ++j)
                ep[j] = 0.25 * (a[j] + b[j] + fp[kv.second[0]][j] + fp[kv.second[1]][j]);
            next[static_cast<std::size_t>(e)] = ep;
            edge_id[kv.first] = e++;
            for (std::int64_t v : {kv.first.first, kv.first.second}) {
                ++valence[static_cast<std::size_t>(v)];
                for (int j = 0; j < 3; ++j)
                    rsum[static_cast<std::size_t>(v)][j] += 0.5 * (a[j] + b[j]);
            }
        }
        for (std::size_t f = 0; f < F.size(); ++f)
            for (std::int64_t v : F[f]) {
                ++nfaces[static_cast<std::size_t>(v)];
                for (int j = 0; j < 3; ++j)
                    fsum[static_cast<std::size_t>(v)][j] += fp[f][j];
            }
        for (std::size_t v = 0; v < nv; ++v) {
            const double n = valence[v];
            for (int j = 0; j < 3; ++j)
                next[v][j] =
                    (fsum[v][j] / nfaces[v] + 2.0 * rsum[v][j] / n + (n - 3.0) * P[v][j]) / n;
        }
        std::vector<std::vector<std::int64_t>> quads;
        for (std::size_t f = 0; f < F.size(); ++f) {
            const std::size_t k = F[f].size();
            for (std::size_t i = 0; i < k; ++i) {
                const std::int64_t prev = F[f][(i + k - 1) % k], cur = F[f][i],
                                   nxt = F[f][(i + 1) % k];
                quads.push_back({cur, edge_id[{std::min(cur, nxt), std::max(cur, nxt)}],
                                 static_cast<std::int64_t>(nv + f),
                                 edge_id[{std::min(prev, cur), std::max(prev, cur)}]});
            }
        }
        P = std::move(next);
        F = std::move(quads);
    }
    // The limit stencil (n^2 v + 4 sum e + sum f) / (n (n + 5)) of a quad mesh.
    std::vector<std::array<double, 3>> esum(P.size(), {0.0, 0.0, 0.0}),
        dsum(P.size(), {0.0, 0.0, 0.0});
    std::vector<int> val(P.size(), 0);
    for (const auto& q : F)
        for (std::size_t i = 0; i < 4; ++i) {
            const auto v = static_cast<std::size_t>(q[i]);
            const auto& nb = P[static_cast<std::size_t>(q[(i + 1) % 4])];
            const auto& dg = P[static_cast<std::size_t>(q[(i + 2) % 4])];
            ++val[v];
            for (int j = 0; j < 3; ++j) {
                esum[v][j] += nb[j];
                dsum[v][j] += dg[j];
            }
        }
    SkinBuilder b;
    for (std::size_t v = 0; v < P.size(); ++v) {
        const double n = val[v];
        double x[3];
        for (int j = 0; j < 3; ++j)
            x[j] = (n * n * P[v][j] + 4.0 * esum[v][j] + dsum[v][j]) / (n * (n + 5.0));
        b.Add(x[0], x[1], x[2]);
    }
    for (const auto& q : F) {
        b.mTris.push_back({q[0], q[1], q[2]});
        b.mTris.push_back({q[0], q[2], q[3]});
    }
    return b.Build();
}

// A Catmull-Clark limit surface with extraordinary vertices of `Valence`: a
// cube cage (valence 3 at its corners) or an n-gonal prism (valence n at its
// cap centres, 3 at its corners).
inline meshioplusplus::Mesh cc_cage(int Valence, int Levels = 4) {
    std::vector<std::array<double, 3>> P;
    std::vector<std::vector<std::int64_t>> F;
    if (Valence == 3) {
        P = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
             {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
        F = {{0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}};
    } else {
        // Cap polygons via integer-free rational points would avoid trig; the
        // n-gon's corners are on the unit circle, so the cage itself uses cos/sin.
        const std::int64_t n = Valence;
        for (int z = 0; z < 2; ++z)
            for (std::int64_t k = 0; k < n; ++k) {
                const double a = 2.0 * kPi * static_cast<double>(k) / static_cast<double>(n);
                P.push_back({std::cos(a), std::sin(a), z ? 0.8 : -0.8});
            }
        std::vector<std::int64_t> bottom, top;
        for (std::int64_t k = 0; k < n; ++k) {
            bottom.push_back(n - 1 - k);
            top.push_back(n + k);
        }
        F.push_back(bottom);
        F.push_back(top);
        for (std::int64_t k = 0; k < n; ++k)
            F.push_back({k, (k + 1) % n, n + (k + 1) % n, n + k});
    }
    return catmull_clark_limit(std::move(P), std::move(F), Levels);
}

// A loaf: a superellipse profile |x|^4 + |z|^4 = 1 over a flat bottom,
// extruded along y and swelling towards the middle, so its top and its two
// curved sides run into each other tangentially (G1) and no primitive
// explains them; the ends and the bottom are planar.
inline meshioplusplus::Mesh g1_loaf(int N = 32, int L = 24) {
    SkinBuilder b;
    std::vector<std::vector<std::int64_t>> rings;
    for (int j = 0; j <= L; ++j) {
        const double y = static_cast<double>(j) / L;
        const double swell = 1.0 + 0.25 * 4.0 * y * (1.0 - y);
        std::vector<std::int64_t> ring;
        for (int i = 0; i <= N; ++i) {
            // t in [0, pi]: x = cos^(1/2), z = sin^(1/2) keeps |x|^4 + |z|^4 = 1.
            const double t = kPi * i / N, c = std::cos(t), s = std::sin(t);
            const double x = (c < 0.0 ? -1.0 : 1.0) * std::sqrt(std::fabs(c));
            const double z = std::sqrt(std::fabs(s));
            ring.push_back(b.Add(0.5 * x, 2.0 * y, 0.6 * swell * z));
        }
        rings.push_back(ring);
    }
    for (int j = 0; j < L; ++j)
        for (int i = 0; i < N; ++i) {
            const std::int64_t a = rings[j][i], bb = rings[j][i + 1], c = rings[j + 1][i + 1],
                               d = rings[j + 1][i];
            const double t = kPi * (i + 0.5) / N;
            b.Quad(a, bb, c, d, {std::cos(t), 0.0, std::sin(t)});
        }
    // The bottom: a strip between the two rims (i = 0 and i = N at z = 0).
    for (int j = 0; j < L; ++j)
        b.Quad(rings[j][0], rings[j][N], rings[j + 1][N], rings[j + 1][0], {0.0, 0.0, -1.0});
    // The ends: fans from the profile's centroid.
    for (int e = 0; e < 2; ++e) {
        const int j = e ? L : 0;
        const double y = 2.0 * j / L;
        const double swell =
            1.0 + 0.25 * 4.0 * (static_cast<double>(j) / L) * (1.0 - static_cast<double>(j) / L);
        const std::int64_t c = b.Add(0.0, y, 0.3 * swell);
        for (int i = 0; i < N; ++i)
            b.Tri(c, rings[j][i], rings[j][i + 1], {0.0, e ? 1.0 : -1.0, 0.0});
        b.Tri(c, rings[j][N], rings[j][0], {0.0, e ? 1.0 : -1.0, 0.0});
    }
    return b.Build();
}

// splitmix64 (Steele, Lea and Flood 2014): a seeded, platform-independent
// stream, the generator `refine_patches` uses for its RANSAC split.
inline std::uint64_t splitmix64(std::uint64_t& rState) {
    std::uint64_t z = (rState += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}

// Gaussian noise of deviation `Sigma` on every coordinate: the Irwin-Hall sum
// of twelve uniforms less six per draw, no libm call, so the C++ and NumPy
// fixtures carry the same noise. A clone of `rSkin` with the moved points (the
// KRATOS backend's meshes are not copy-constructible).
inline meshioplusplus::Mesh gaussian_noise(const meshioplusplus::Mesh& rSkin, double Sigma,
                                           std::uint64_t Seed = 7) {
    meshioplusplus::NDArray q =
        meshioplusplus::NDArray::Uninit(meshioplusplus::DType::Float64, {rSkin.NumPoints(), 3});
    double* a = q.As<double>();
    std::copy(rSkin.Points().As<double>(), rSkin.Points().As<double>() + rSkin.NumPoints() * 3, a);
    std::uint64_t state = Seed;
    for (std::size_t i = 0; i < rSkin.NumPoints() * 3; ++i) {
        double g = -6.0;
        for (int k = 0; k < 12; ++k)
            g += static_cast<double>(splitmix64(state) >> 11) * 0x1.0p-53;
        a[i] += Sigma * g;
    }
    meshioplusplus::Mesh out = meshioplusplus::detail::clone_mesh(rSkin);
    out.AssignPoints(std::move(q));
    return out;
}

}  // namespace analytic
