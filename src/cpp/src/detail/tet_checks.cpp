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
// The volume tolerance and the face-overlap count of tetrahedra.

// System includes
#include <algorithm>
#include <cmath>
#include <cstddef>

// Project includes
#include "tet_checks.hpp"

namespace meshioplusplus::detail {

std::pair<double, double> tet_volume_and_tolerance(const std::vector<Vec3>& rPoints,
                                                   const TetNodes& rTet,
                                                   const double ZeroVolumeFactor) {
    // Relative precision of the written coordinates (six significant digits): a volume below its
    // effect on the longest edge squared is noise.
    constexpr double kCoordinatePrecision = 1.0e-6;
    const auto point = [&](const std::int64_t I) -> const Vec3& {
        return rPoints[static_cast<std::size_t>(I)];
    };
    double longest = 0.0;
    double largest_coordinate = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (const double c : point(rTet[i]))
            largest_coordinate = std::max(largest_coordinate, std::abs(c));
        for (int j = i + 1; j < 4; ++j)
            longest = std::max(longest, vec3_norm(vec3_sub(point(rTet[j]), point(rTet[i]))));
    }
    // The arithmetic of the octree volume fix this was lifted from (the operations'
    // outputs are pinned bit for bit): (b - a) x (c - a) . (d - a).
    const double volume = vec3_dot(vec3_cross(vec3_sub(point(rTet[1]), point(rTet[0])),
                                              vec3_sub(point(rTet[2]), point(rTet[0]))),
                                   vec3_sub(point(rTet[3]), point(rTet[0]))) /
                          6.0;
    return {volume, ZeroVolumeFactor * longest * longest * longest +
                        kCoordinatePrecision * largest_coordinate * longest * longest};
}

std::int64_t count_tet_overlaps(const std::vector<Vec3>& rPoints,
                                const std::vector<TetNodes>& rTets,
                                const std::vector<std::uint8_t>* pSkip) {
    // The four faces of a tetrahedron, each opposite the node of the same index.
    static constexpr int kFaces[4][3] = {{0, 1, 3}, {1, 2, 3}, {2, 0, 3}, {0, 2, 1}};
    using Index = std::int64_t;
    const auto point = [&](const Index I) -> const Vec3& {
        return rPoints[static_cast<std::size_t>(I)];
    };
    std::vector<std::pair<std::array<Index, 3>, std::pair<std::size_t, Index>>> faces;
    faces.reserve(rTets.size() * 4);
    for (std::size_t t = 0; t < rTets.size(); ++t) {
        if (pSkip && (*pSkip)[t])
            continue;
        for (const auto& facet : kFaces) {
            std::array<Index, 3> key{rTets[t][facet[0]], rTets[t][facet[1]], rTets[t][facet[2]]};
            std::sort(key.begin(), key.end());
            // The node of the tetrahedron which is not on the face.
            Index opposite = -1;
            for (const Index v : rTets[t])
                if (v != key[0] && v != key[1] && v != key[2])
                    opposite = v;
            faces.push_back({key, {t, opposite}});
        }
    }
    std::sort(faces.begin(), faces.end());
    std::int64_t overlaps = 0;
    // Only a face used by exactly two tetrahedra is a pair: a run of equal keys of any other
    // length (a face three tetrahedra share) is skipped whole.
    for (std::size_t i = 0, run_end = 0; i + 1 < faces.size(); i = run_end) {
        run_end = i + 1;
        while (run_end < faces.size() && faces[run_end].first == faces[i].first)
            ++run_end;
        if (run_end - i != 2)
            continue;
        const Vec3& a = point(faces[i].first[0]);
        const Vec3 normal = vec3_cross(vec3_sub(point(faces[i].first[1]), a),
                                       vec3_sub(point(faces[i].first[2]), a));
        const double side_0 = vec3_dot(normal, vec3_sub(point(faces[i].second.second), a));
        const double side_1 = vec3_dot(normal, vec3_sub(point(faces[i + 1].second.second), a));
        overlaps += side_0 * side_1 > 0.0;
    }
    return overlaps;
}

}  // namespace meshioplusplus::detail
