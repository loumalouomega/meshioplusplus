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
// The external tetrahedral meshers (roadmap 7.7.2) and quad remeshers (7.7.3): TetGen and Gmsh,
// when installed, run as separate processes on a region's surface and read back, registered by name
// (`tetrahedralize(..., backend="tetgen")`). Running them as programs, never linking them, is the
// route the licence decision on TetGen (AGPL) and Gmsh (GPL) leaves open. Never in WebAssembly.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iterator>
#include <map>
#include <stdexcept>
#include <string>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/detail/classic_stream.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/formats/gmsh.hpp"
#include "meshioplusplus/formats/tetgen.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/log.hpp"
#include "meshioplusplus/providers.hpp"

// Project includes (private, not installed)
#include "bounding_box_utils.hpp"
#include "geometry_utils.hpp"
#include "subprocess.hpp"

namespace meshioplusplus {

namespace {

using EmPoint = std::array<double, 3>;

// The run's directory, removed afterwards unless MESHIOPLUSPLUS_KEEP_TEMP is set.
struct EmWorkDir {
    std::string mPath;
    explicit EmWorkDir(const std::string& rTag) : mPath(detail::make_temp_dir(rTag)) {}
    EmWorkDir(const EmWorkDir&) = delete;
    EmWorkDir& operator=(const EmWorkDir&) = delete;
    ~EmWorkDir() {
        const char* p_keep = std::getenv("MESHIOPLUSPLUS_KEEP_TEMP");
        if (p_keep && *p_keep && std::string(p_keep) != "0") {
            log::info("external mesher: kept '{}'", mPath);
            return;
        }
        std::error_code ec;
        std::filesystem::remove_all(mPath, ec);
    }
    std::string File(const std::string& rName) const {
        return (std::filesystem::path(mPath) / rName).string();
    }
};

std::string em_num(double Value) {
    char buf[40];
    detail::snprintf_c(buf, sizeof(buf), "%.17g", Value);
    return buf;
}

// The last lines of the program's log, for an error message.
std::string em_log_tail(const std::string& rLog) {
    auto in = detail::make_classic_ifstream(rLog);
    std::vector<std::string> lines;
    for (std::string line; std::getline(in, line);)
        if (!line.empty())
            lines.push_back(line);
    std::string out;
    for (std::size_t i = lines.size() > 6 ? lines.size() - 6 : 0; i < lines.size(); ++i)
        out += (out.empty() ? "" : " | ") + lines[i];
    return out;
}

void em_run(const std::string& rName, const std::vector<std::string>& rArgv,
            const EmWorkDir& rDir) {
    const std::string log = rDir.File("mesher.log");
    const int code = detail::run_process(rArgv, rDir.mPath, log);
    if (code != 0)
        throw std::runtime_error(rName + " failed (exit code " + std::to_string(code) +
                                 "): " + em_log_tail(log));
}

/**
 * The region's result from the mesher's points and tetrahedra: the request's
 * points first, exactly as given, each matched to the mesher's copy of it
 * (which a text round trip may have rounded in the last digit), then the
 * points the mesher added. `mBoundaryPreserved` holds when every request
 * point came back and the tetrahedra's boundary is the request's triangles.
 */
TetMeshRegionResult em_assemble(const TetMeshRegionRequest& rRequest,
                                const std::vector<EmPoint>& rOut,
                                const std::vector<std::array<std::int64_t, 4>>& rTets) {
    const std::size_t n_in = rRequest.mPoints.size();
    const double diag =
        detail::BoundingBoxUtils::Diagonal(detail::BoundingBoxUtils::FromPoints(rRequest.mPoints));
    const double tol = std::max(1e-9 * diag, 1e-300);
    const double cell = 4.0 * tol;
    struct KeyHash {
        std::size_t operator()(const std::array<std::int64_t, 3>& rK) const noexcept {
            return static_cast<std::size_t>(rK[0] * 73856093LL ^ rK[1] * 19349663LL ^
                                            rK[2] * 83492791LL);
        }
    };
    const auto key = [&](const EmPoint& rP) {
        return std::array<std::int64_t, 3>{static_cast<std::int64_t>(std::floor(rP[0] / cell)),
                                           static_cast<std::int64_t>(std::floor(rP[1] / cell)),
                                           static_cast<std::int64_t>(std::floor(rP[2] / cell))};
    };
    std::unordered_map<std::array<std::int64_t, 3>, std::vector<std::size_t>, KeyHash> grid;
    for (std::size_t i = 0; i < n_in; ++i)
        grid[key(rRequest.mPoints[i])].push_back(i);

    TetMeshRegionResult result;
    result.mPoints = rRequest.mPoints;
    std::vector<std::int64_t> map(rOut.size(), -1);
    std::vector<bool> matched(n_in, false);
    for (std::size_t o = 0; o < rOut.size(); ++o) {
        const auto k = key(rOut[o]);
        std::int64_t best = -1;
        double best_d = tol * tol;
        for (std::int64_t dx = -1; dx <= 1; ++dx)
            for (std::int64_t dy = -1; dy <= 1; ++dy)
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    const auto it = grid.find({k[0] + dx, k[1] + dy, k[2] + dz});
                    if (it == grid.end())
                        continue;
                    for (const std::size_t i : it->second) {
                        double d = 0.0;
                        for (std::size_t c = 0; c < 3; ++c)
                            d += (rOut[o][c] - rRequest.mPoints[i][c]) *
                                 (rOut[o][c] - rRequest.mPoints[i][c]);
                        if (d <= best_d && !matched[i]) {
                            best_d = d;
                            best = static_cast<std::int64_t>(i);
                        }
                    }
                }
        if (best >= 0) {
            matched[static_cast<std::size_t>(best)] = true;
            map[o] = best;
        } else {
            map[o] = static_cast<std::int64_t>(result.mPoints.size());
            result.mPoints.push_back(rOut[o]);
        }
    }
    result.mBoundaryPreserved =
        std::all_of(matched.begin(), matched.end(), [](bool B) { return B; });
    result.mTets.reserve(rTets.size());
    for (const auto& r_t : rTets) {
        std::array<std::int64_t, 4> t{};
        for (std::size_t j = 0; j < 4; ++j) {
            if (r_t[j] < 0 || static_cast<std::size_t>(r_t[j]) >= map.size())
                throw std::runtime_error("the mesher returned a tetrahedron naming no point");
            t[j] = map[static_cast<std::size_t>(r_t[j])];
        }
        result.mTets.push_back(t);
    }
    // The boundary: faces of one tetrahedron only, against the request's.
    if (result.mBoundaryPreserved) {
        std::map<std::array<std::int64_t, 3>, int> faces;
        static constexpr int kFace[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
        for (const auto& r_t : result.mTets)
            for (const auto& r_f : kFace) {
                std::array<std::int64_t, 3> f{r_t[r_f[0]], r_t[r_f[1]], r_t[r_f[2]]};
                std::sort(f.begin(), f.end());
                ++faces[f];
            }
        std::vector<std::array<std::int64_t, 3>> got, want;
        for (const auto& r_e : faces)
            if (r_e.second == 1)
                got.push_back(r_e.first);
        for (const auto& r_t : rRequest.mTriangles) {
            std::array<std::int64_t, 3> f = r_t;
            std::sort(f.begin(), f.end());
            want.push_back(f);
        }
        std::sort(want.begin(), want.end());
        result.mBoundaryPreserved = got == want;
    }
    return result;
}

/**
 * The region's shells (triangles connected across edges), the outer one first
 * (the largest bounding box), each as triangle indices; and for each inner
 * shell (a hole of the region) a point inside it: a triangle's centre moved a
 * little along its normal, which faces out of the region, so into the hole.
 */
struct EmShells {
    std::vector<std::vector<std::size_t>> mShells;
    std::vector<EmPoint> mHolePoints;
};

EmShells em_shells(const TetMeshRegionRequest& rRequest) {
    const std::size_t nt = rRequest.mTriangles.size();
    std::vector<std::size_t> parent(nt);
    for (std::size_t t = 0; t < nt; ++t)
        parent[t] = t;
    const auto find = [&](std::size_t A) {
        while (parent[A] != A)
            A = parent[A] = parent[parent[A]];
        return A;
    };
    std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> edge_owner;
    for (std::size_t t = 0; t < nt; ++t)
        for (std::size_t k = 0; k < 3; ++k) {
            std::int64_t a = rRequest.mTriangles[t][k], b = rRequest.mTriangles[t][(k + 1) % 3];
            if (a > b)
                std::swap(a, b);
            const auto [it, fresh] = edge_owner.emplace(std::make_pair(a, b), t);
            if (!fresh)
                parent[find(t)] = find(it->second);
        }
    std::map<std::size_t, std::size_t> root_index;
    EmShells out;
    for (std::size_t t = 0; t < nt; ++t) {
        const std::size_t r = find(t);
        const auto [it, fresh] = root_index.emplace(r, out.mShells.size());
        if (fresh)
            out.mShells.emplace_back();
        out.mShells[it->second].push_back(t);
    }
    const auto extent = [&](const std::vector<std::size_t>& rShell) {
        EmPoint lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
        for (const std::size_t t : rShell)
            for (const std::int64_t v : rRequest.mTriangles[t])
                for (std::size_t c = 0; c < 3; ++c) {
                    lo[c] = std::min(lo[c], rRequest.mPoints[static_cast<std::size_t>(v)][c]);
                    hi[c] = std::max(hi[c], rRequest.mPoints[static_cast<std::size_t>(v)][c]);
                }
        return (hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) +
               (hi[2] - lo[2]) * (hi[2] - lo[2]);
    };
    std::size_t outer = 0;
    for (std::size_t s = 1; s < out.mShells.size(); ++s)
        if (extent(out.mShells[s]) > extent(out.mShells[outer]))
            outer = s;
    std::swap(out.mShells[0], out.mShells[outer]);
    for (std::size_t s = 1; s < out.mShells.size(); ++s) {
        // The largest triangle of the shell, for a step well inside its size.
        std::size_t best = out.mShells[s][0];
        double best_area = -1.0;
        EmPoint best_n{};
        for (const std::size_t t : out.mShells[s]) {
            const auto& r_t = rRequest.mTriangles[t];
            const EmPoint& a = rRequest.mPoints[static_cast<std::size_t>(r_t[0])];
            const EmPoint& b = rRequest.mPoints[static_cast<std::size_t>(r_t[1])];
            const EmPoint& c = rRequest.mPoints[static_cast<std::size_t>(r_t[2])];
            const EmPoint n = detail::GeometryUtils::TriangleNormal(a, b, c);
            const double area = detail::MathUtils::Norm(n);
            if (area > best_area) {
                best_area = area;
                best = t;
                best_n = n;
            }
        }
        const auto& r_t = rRequest.mTriangles[best];
        EmPoint centre{};
        for (const std::int64_t v : r_t)
            for (std::size_t c = 0; c < 3; ++c)
                centre[c] += rRequest.mPoints[static_cast<std::size_t>(v)][c] / 3.0;
        const double len = best_area > 0.0 ? best_area : 1.0;
        // A step of a tenth of the triangle's size, along the unit normal.
        const double step = 0.1 * std::sqrt(0.5 * best_area);
        out.mHolePoints.push_back({centre[0] + step * best_n[0] / len,
                                   centre[1] + step * best_n[1] / len,
                                   centre[2] + step * best_n[2] / len});
    }
    return out;
}

// TetGen 1.5/1.6 meshes a piecewise-linear complex (`-p`) without adding
// points on its boundary (`-Y`), quietly (`-Q`), with quality (`-q`) when the
// optimisation passes are asked for and a volume bound (`-a`) from the size.
TetMeshRegionResult em_tetgen(const std::string& rExe, const TetMeshRegionRequest& rRequest) {
    const EmWorkDir dir("tetgen");
    Mesh surface;
    {
        const std::size_t n = rRequest.mPoints.size();
        NDArray p(DType::Float64, {n, std::size_t{3}});
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t c = 0; c < 3; ++c)
                p.As<double>()[3 * i + c] = rRequest.mPoints[i][c];
        surface.AssignPoints(std::move(p));
        NDArray t(DType::Int64, {rRequest.mTriangles.size(), std::size_t{3}});
        for (std::size_t f = 0; f < rRequest.mTriangles.size(); ++f)
            for (std::size_t c = 0; c < 3; ++c)
                t.As<std::int64_t>()[3 * f + c] = rRequest.mTriangles[f][c];
        surface.AddCellBlock("triangle", std::move(t));
    }
    const std::string in = dir.File("in.smesh");
    write_tetgen(in, surface);
    {
        // The holes: write_tetgen's `.smesh` ends with empty hole and region
        // lists, which a region with inner shells replaces.
        const EmShells shells = em_shells(rRequest);
        if (!shells.mHolePoints.empty()) {
            auto text_in = detail::make_classic_ifstream(in);
            std::string text((std::istreambuf_iterator<char>(text_in)),
                             std::istreambuf_iterator<char>());
            text_in.close();
            const std::string tail = "# part 3: holes\n0\n";
            const std::size_t at = text.rfind(tail);
            if (at == std::string::npos)
                throw std::runtime_error("tetgen: unexpected .smesh layout");
            std::string holes =
                "# part 3: holes\n" + std::to_string(shells.mHolePoints.size()) + "\n";
            for (std::size_t h = 0; h < shells.mHolePoints.size(); ++h)
                holes += std::to_string(h) + " " + em_num(shells.mHolePoints[h][0]) + " " +
                         em_num(shells.mHolePoints[h][1]) + " " + em_num(shells.mHolePoints[h][2]) +
                         "\n";
            text.replace(at, tail.size(), holes);
            auto text_out = detail::make_classic_ofstream(in, std::ios::binary);
            text_out << text;
        }
    }
    std::string switches = "-pYQ";
    if (rRequest.mOptimizeSteps > 0)
        switches += "q";
    if (rRequest.mMaxSize > 0.0) {
        // The volume of a regular tetrahedron of that edge.
        const double h = rRequest.mMaxSize;
        switches += "a" + em_num(h * h * h / (6.0 * std::sqrt(2.0)));
    }
    em_run("tetgen", {rExe, switches, in}, dir);
    const Mesh out = read_tetgen(dir.File("in.1.node"));
    std::vector<EmPoint> pts(out.NumPoints());
    const NDArray& r_p = out.Points();
    for (std::size_t i = 0; i < pts.size(); ++i)
        for (std::size_t c = 0; c < 3; ++c)
            pts[i][c] = detail::read_double(r_p, 3 * i + c);
    std::vector<std::array<std::int64_t, 4>> tets;
    for (std::size_t b = 0; b < out.NumCellBlocks(); ++b) {
        const auto cells = out.Cells(b);
        if (cells.Type() != "tetra")
            continue;
        for (std::size_t c = 0; c < cells.NumCells(); ++c)
            tets.push_back({detail::read_int(cells.Conn(), 4 * c),
                            detail::read_int(cells.Conn(), 4 * c + 1),
                            detail::read_int(cells.Conn(), 4 * c + 2),
                            detail::read_int(cells.Conn(), 4 * c + 3)});
    }
    return em_assemble(rRequest, pts, tets);
}

// Gmsh meshes the volume a discrete surface (an ASCII STL, full precision)
// bounds, keeping its triangles: `Surface Loop` and `Volume` over it, `-3`.
TetMeshRegionResult em_gmsh(const std::string& rExe, const TetMeshRegionRequest& rRequest) {
    const EmWorkDir dir("gmsh");
    const EmShells shells = em_shells(rRequest);
    {
        // One `solid` per shell: Gmsh makes each a discrete surface, numbered
        // from 1 in file order.
        auto stl = detail::make_classic_ofstream(dir.File("in.stl"));
        for (const auto& r_shell : shells.mShells) {
            stl << "solid region\n";
            for (const std::size_t t_index : r_shell) {
                const auto& r_t = rRequest.mTriangles[t_index];
                const EmPoint& a = rRequest.mPoints[static_cast<std::size_t>(r_t[0])];
                const EmPoint& b = rRequest.mPoints[static_cast<std::size_t>(r_t[1])];
                const EmPoint& c = rRequest.mPoints[static_cast<std::size_t>(r_t[2])];
                const EmPoint n = detail::MathUtils::UnitByDivision(
                    detail::GeometryUtils::TriangleNormal(a, b, c));
                stl << "facet normal " << em_num(n[0]) << " " << em_num(n[1]) << " " << em_num(n[2])
                    << "\n outer loop\n";
                for (const EmPoint* p_v : {&a, &b, &c})
                    stl << "  vertex " << em_num((*p_v)[0]) << " " << em_num((*p_v)[1]) << " "
                        << em_num((*p_v)[2]) << "\n";
                stl << " endloop\nendfacet\n";
            }
            stl << "endsolid region\n";
        }
        auto geo = detail::make_classic_ofstream(dir.File("in.geo"));
        geo << "Merge \"in.stl\";\n";
        std::string loops;
        for (std::size_t k = 1; k <= shells.mShells.size(); ++k) {
            geo << "Surface Loop(" << k << ") = {" << k << "};\n";
            loops += (k > 1 ? "," : "") + std::to_string(k);
        }
        geo << "Volume(1) = {" << loops << "};\n";
        if (rRequest.mOptimizeSteps <= 0)
            geo << "Mesh.Optimize = 0;\n";
    }
    std::vector<std::string> argv = {rExe,      "-3",    "-v", "2",
                                     "-format", "msh41", "-o", dir.File("out.msh")};
    if (rRequest.mMaxSize > 0.0) {
        argv.push_back("-clmax");
        argv.push_back(em_num(rRequest.mMaxSize));
    }
    argv.push_back(dir.File("in.geo"));
    em_run("gmsh", argv, dir);
    const Mesh out = read_gmsh(dir.File("out.msh"));
    std::vector<EmPoint> pts(out.NumPoints());
    const NDArray& r_p = out.Points();
    for (std::size_t i = 0; i < pts.size(); ++i)
        for (std::size_t c = 0; c < 3; ++c)
            pts[i][c] = detail::read_double(r_p, 3 * i + c);
    std::vector<std::array<std::int64_t, 4>> tets;
    for (std::size_t b = 0; b < out.NumCellBlocks(); ++b) {
        const auto cells = out.Cells(b);
        if (cells.Type() != "tetra")
            continue;
        for (std::size_t c = 0; c < cells.NumCells(); ++c)
            tets.push_back({detail::read_int(cells.Conn(), 4 * c),
                            detail::read_int(cells.Conn(), 4 * c + 1),
                            detail::read_int(cells.Conn(), 4 * c + 2),
                            detail::read_int(cells.Conn(), 4 * c + 3)});
    }
    if (tets.empty())
        throw std::runtime_error("gmsh produced no tetrahedra: " +
                                 em_log_tail(dir.File("mesher.log")));
    return em_assemble(rRequest, pts, tets);
}

// --- quad remeshers (roadmap 7.7.3) ---------------------------------------------

// The surface as a Wavefront OBJ at full precision.
void em_write_obj(const std::string& rPath, const QuadRemeshRequest& rRequest) {
    auto obj = detail::make_classic_ofstream(rPath);
    for (const EmPoint& r_p : rRequest.mPoints)
        obj << "v " << em_num(r_p[0]) << " " << em_num(r_p[1]) << " " << em_num(r_p[2]) << "\n";
    for (const auto& r_t : rRequest.mTriangles)
        obj << "f " << r_t[0] + 1 << " " << r_t[1] + 1 << " " << r_t[2] + 1 << "\n";
    if (!obj)
        throw std::runtime_error("cannot write '" + rPath + "'");
}

// The `v` and `f` records of the program's OBJ (`f` corners as `i`, `i/t` or
// `i/t/n`, negative indices relative): quads, and triangles where the program
// left them. A larger face is refused.
QuadRemeshResult em_read_obj(const std::string& rPath, const char* pProgram) {
    auto in = detail::make_classic_ifstream(rPath);
    if (!in)
        throw std::runtime_error(std::string(pProgram) + " wrote no '" + rPath + "'");
    QuadRemeshResult out;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 2 || line[1] != ' ' || (line[0] != 'v' && line[0] != 'f'))
            continue;
        std::vector<std::string> tokens;
        for (std::size_t at = 2; at < line.size();) {
            while (at < line.size() && (line[at] == ' ' || line[at] == '\t' || line[at] == '\r'))
                ++at;
            std::size_t end = at;
            while (end < line.size() && line[end] != ' ' && line[end] != '\t' && line[end] != '\r')
                ++end;
            if (end > at)
                tokens.push_back(line.substr(at, end - at));
            at = end;
        }
        if (line[0] == 'v') {
            if (tokens.size() < 3)
                throw std::runtime_error(std::string(pProgram) +
                                         " wrote a vertex without 3 coordinates");
            EmPoint p{};
            for (std::size_t c = 0; c < 3; ++c) {
                p[c] = detail::parse_double(tokens[c]);
                if (!std::isfinite(p[c]))
                    throw std::runtime_error(std::string(pProgram) + " wrote a non-finite vertex");
            }
            out.mPoints.push_back(p);
            continue;
        }
        std::vector<std::int64_t> face;
        for (const std::string& r_tok : tokens) {
            const std::string head = r_tok.substr(0, r_tok.find('/'));
            char* p_end = nullptr;
            const long long v = std::strtoll(head.c_str(), &p_end, 10);
            if (head.empty() || p_end != head.c_str() + head.size() || v == 0)
                throw std::runtime_error(std::string(pProgram) + " wrote a bad face corner '" +
                                         r_tok + "'");
            const auto n = static_cast<long long>(out.mPoints.size());
            face.push_back(v > 0 ? v - 1 : n + v);
        }
        for (const std::int64_t v : face)
            if (v < 0 || v >= static_cast<std::int64_t>(out.mPoints.size()))
                throw std::runtime_error(std::string(pProgram) + " wrote a face naming no vertex");
        if (face.size() == 4)
            out.mQuads.push_back({face[0], face[1], face[2], face[3]});
        else if (face.size() == 3)
            out.mTriangles.push_back({face[0], face[1], face[2]});
        else
            throw std::runtime_error(std::string(pProgram) + " wrote a face of " +
                                     std::to_string(face.size()) + " corners");
    }
    if (out.mQuads.empty() && out.mTriangles.empty())
        throw std::runtime_error(std::string(pProgram) + " wrote no face");
    return out;
}

// QuadriFlow (Huang et al. 2018): `quadriflow -i in.obj -o out.obj -f N`, an
// all-quad mesh.
QuadRemeshResult em_quadriflow(const std::string& rExe, const QuadRemeshRequest& rRequest) {
    const EmWorkDir dir("quadriflow");
    em_write_obj(dir.File("in.obj"), rRequest);
    std::vector<std::string> argv = {rExe};
    if (rRequest.mSharp)
        argv.push_back("-sharp");
    if (rRequest.mBoundary)
        argv.push_back("-boundary");
    if (rRequest.mAdaptive)
        argv.push_back("-adaptive");
    argv.insert(argv.end(),
                {"-seed", std::to_string(rRequest.mSeed), "-i", dir.File("in.obj"), "-o",
                 dir.File("out.obj"), "-f", std::to_string(rRequest.mTargetFaces)});
    em_run("quadriflow", argv, dir);
    return em_read_obj(dir.File("out.obj"), "quadriflow");
}

// Instant Meshes (Jakob et al. 2015) in batch mode: `-o out.obj -f N -d`,
// deterministic; `-D` quad-dominant, `-b` boundaries kept, `-c` creases.
QuadRemeshResult em_instant_meshes(const std::string& rExe, const QuadRemeshRequest& rRequest) {
    const EmWorkDir dir("instant_meshes");
    em_write_obj(dir.File("in.obj"), rRequest);
    std::vector<std::string> argv = {
        rExe, "-d", "-o", dir.File("out.obj"), "-f", std::to_string(rRequest.mTargetFaces)};
    if (rRequest.mDominant)
        argv.push_back("-D");
    if (rRequest.mBoundary)
        argv.push_back("-b");
    if (rRequest.mSharp) {
        argv.push_back("-c");
        argv.push_back("30");
    }
    argv.push_back(dir.File("in.obj"));
    em_run("instant_meshes", argv, dir);
    return em_read_obj(dir.File("out.obj"), "instant_meshes");
}

}  // namespace

namespace providers {

std::vector<TetMesherProvider> make_external_tet_meshers() {
    std::vector<TetMesherProvider> out;
    if (!detail::subprocess_available())
        return out;
    const char* p_off = std::getenv("MESHIOPLUSPLUS_NO_EXTERNAL_MESHERS");
    if (p_off && *p_off && std::string(p_off) != "0")
        return out;
    const std::string tetgen = detail::find_executable("tetgen", "MESHIOPLUSPLUS_TETGEN");
    if (!tetgen.empty()) {
        TetMesherProvider p;
        p.mName = "tetgen";
        p.mVersion = "external: " + tetgen;
        p.mMeshRegion = [tetgen](const TetMeshRegionRequest& rRequest) {
            return em_tetgen(tetgen, rRequest);
        };
        out.push_back(std::move(p));
    }
    const std::string gmsh = detail::find_executable("gmsh", "MESHIOPLUSPLUS_GMSH");
    if (!gmsh.empty()) {
        TetMesherProvider p;
        p.mName = "gmsh";
        p.mVersion = "external: " + gmsh;
        p.mMeshRegion = [gmsh](const TetMeshRegionRequest& rRequest) {
            return em_gmsh(gmsh, rRequest);
        };
        out.push_back(std::move(p));
    }
    return out;
}

std::vector<QuadRemesherProvider> make_external_quad_remeshers() {
    std::vector<QuadRemesherProvider> out;
    if (!detail::subprocess_available())
        return out;
    const char* p_off = std::getenv("MESHIOPLUSPLUS_NO_EXTERNAL_MESHERS");
    if (p_off && *p_off && std::string(p_off) != "0")
        return out;
    const std::string quadriflow =
        detail::find_executable("quadriflow", "MESHIOPLUSPLUS_QUADRIFLOW");
    if (!quadriflow.empty()) {
        QuadRemesherProvider p;
        p.mName = "quadriflow";
        p.mVersion = "external: " + quadriflow;
        p.mRemesh = [quadriflow](const QuadRemeshRequest& rRequest) {
            return em_quadriflow(quadriflow, rRequest);
        };
        out.push_back(std::move(p));
    }
    std::string instant =
        detail::find_executable("Instant Meshes", "MESHIOPLUSPLUS_INSTANT_MESHES");
    for (const char* p_alt : {"instant-meshes", "InstantMeshes"})
        if (instant.empty())
            instant = detail::find_executable(p_alt);
    if (!instant.empty()) {
        QuadRemesherProvider p;
        p.mName = "instant_meshes";
        p.mVersion = "external: " + instant;
        p.mRemesh = [instant](const QuadRemeshRequest& rRequest) {
            return em_instant_meshes(instant, rRequest);
        };
        out.push_back(std::move(p));
    }
    return out;
}

}  // namespace providers

}  // namespace meshioplusplus
