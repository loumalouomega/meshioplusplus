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
