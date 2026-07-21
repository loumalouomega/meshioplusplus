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
#include <cmath>
#include <string>

// Project includes
#include "meshioplusplus/detail/projection.hpp"
#include "meshioplusplus/detail/value_io.hpp"

namespace meshioplusplus {
namespace detail {

CameraBasis camera_basis(double azimuth, double elevation, double roll) {
    constexpr double deg = 3.141592653589793 / 180.0;
    const double az = azimuth * deg;
    const double el = elevation * deg;
    const double caz = std::cos(az);
    const double saz = std::sin(az);
    const double cel = std::cos(el);
    const double sel = std::sin(el);
    const double w0 = cel * caz;
    const double w1 = cel * saz;
    const double w2 = sel;

    double u0, u1, u2;
    if (std::fabs(w2) > 1.0 - 1.0e-12) {
        // Looking straight along z: use y as the up reference.
        // u = normalize(cross((0,1,0), w))
        const double n = std::sqrt(w2 * w2 + w0 * w0);
        u0 = w2 / n;
        u1 = 0.0;
        u2 = -w0 / n;
    } else {
        // u = normalize(cross((0,0,1), w))
        const double n = std::sqrt(w1 * w1 + w0 * w0);
        u0 = -w1 / n;
        u1 = w0 / n;
        u2 = 0.0;
    }
    // v = cross(w, u)
    double v0 = w1 * u2 - w2 * u1;
    double v1 = w2 * u0 - w0 * u2;
    double v2 = w0 * u1 - w1 * u0;

    if (roll != 0.0) {
        const double cr = std::cos(roll * deg);
        const double sr = std::sin(roll * deg);
        const double ru0 = cr * u0 + sr * v0;
        const double ru1 = cr * u1 + sr * v1;
        const double ru2 = cr * u2 + sr * v2;
        const double rv0 = -sr * u0 + cr * v0;
        const double rv1 = -sr * u1 + cr * v1;
        const double rv2 = -sr * u2 + cr * v2;
        u0 = ru0;
        u1 = ru1;
        u2 = ru2;
        v0 = rv0;
        v1 = rv1;
        v2 = rv2;
    }
    return CameraBasis{{u0, u1, u2}, {v0, v1, v2}, {w0, w1, w2}};
}

ProjectedSurface project_surface(const Mesh& rMesh, double azimuth, double elevation, double roll) {
    const CameraBasis cam = camera_basis(azimuth, elevation, roll);
    const NDArray& points = rMesh.Points();
    const std::size_t num_points = rMesh.NumPoints();
    const std::size_t dim = rMesh.PointDim();

    ProjectedSurface out;
    out.mX.resize(num_points);
    out.mY.resize(num_points);
    std::vector<double> depth(num_points);
    for (std::size_t i = 0; i < num_points; ++i) {
        const double px = (0 < dim) ? read_double(points, i * dim + 0) : 0.0;
        const double py = (1 < dim) ? read_double(points, i * dim + 1) : 0.0;
        const double pz = (2 < dim) ? read_double(points, i * dim + 2) : 0.0;
        out.mX[i] = px * cam.mU[0] + py * cam.mU[1] + pz * cam.mU[2];
        out.mY[i] = px * cam.mV[0] + py * cam.mV[1] + pz * cam.mV[2];
        depth[i] = px * cam.mW[0] + py * cam.mW[1] + pz * cam.mW[2];
    }

    // Counts over EVERY block, including the ones skipped below, so mSourceCell
    // indexes the projected mesh's cells globally -- the same convention
    // "surface:parent_cell" uses, which is what lets the two indices compose.
    std::int64_t global_cell_base = 0;
    for (const auto cb : rMesh.CellRange()) {
        const std::int64_t block_base = global_cell_base;
        global_cell_base += static_cast<std::int64_t>(cb.NumCells());
        const std::string& type = cb.Type();
        std::uint8_t n_corner = 0;
        bool is_line = false;
        if (type == "line") {
            n_corner = 2;
            is_line = true;
        } else if (type == "triangle" || type == "triangle6") {
            n_corner = 3;
        } else if (type == "quad" || type == "quad8" || type == "quad9") {
            n_corner = 4;
        } else {
            continue;
        }
        const NDArray& conn = cb.Conn();
        const std::size_t ncols = cols(conn);
        const std::size_t n = cb.NumCells();
        for (std::size_t r = 0; r < n; ++r) {
            ProjectedFace face;
            face.mNodes = {-1, -1, -1, -1};
            face.mNumNodes = n_corner;
            face.mIsLine = is_line;
            face.mSourceCell = block_base + static_cast<std::int64_t>(r);
            double d = 0.0;
            for (std::uint8_t k = 0; k < n_corner; ++k) {
                const std::int64_t p = read_int(conn, r * ncols + k);
                face.mNodes[k] = p;
                d += depth[static_cast<std::size_t>(p)];
            }
            face.mDepth = d / static_cast<double>(n_corner);
            out.mFaces.push_back(face);
        }
    }

    std::stable_sort(
        out.mFaces.begin(), out.mFaces.end(),
        [](const ProjectedFace& rA, const ProjectedFace& rB) { return rA.mDepth < rB.mDepth; });
    return out;
}

}  // namespace detail
}  // namespace meshioplusplus
