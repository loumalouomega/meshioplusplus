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
// Curves to STEP. See detail/step_wireframe.hpp.

// System includes
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

// Project includes
#include "step_wireframe.hpp"
#include "polyline_chain.hpp"
#include "sketch_pieces.hpp"
#include "step_stream.hpp"
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/detail/length_unit.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace detail {

namespace {

constexpr const char* kSwPrefix = "meshio++: step: ";
constexpr double kSwPi = 3.14159265358979323846;

std::array<double, 3> sw_cross(const std::array<double, 3>& rA, const std::array<double, 3>& rB) {
    return {rA[1] * rB[2] - rA[2] * rB[1], rA[2] * rB[0] - rA[0] * rB[2],
            rA[0] * rB[1] - rA[1] * rB[0]};
}

double sw_dot(const std::array<double, 3>& rA, const std::array<double, 3>& rB) {
    return rA[0] * rB[0] + rA[1] * rB[1] + rA[2] * rB[2];
}

}  // namespace

bool step_is_wireframe(const Mesh& rMesh) {
    bool lines = false;
    for (const auto cb : rMesh.CellRange()) {
        if (cb.IsPolyhedron())
            return false;
        const std::string type(cb.Type());
        if (type.rfind("polygon", 0) == 0)
            return false;
        const int dim = cell_type_dimension(cell_type_from_name(type));
        if (dim >= 2)
            return false;
        lines = lines || (dim == 1 && cb.NumCells() > 0);
    }
    return lines;
}

void step_write_wireframe(const std::string& rPath, const Mesh& rMesh) {
    StepStream st(rPath, "meshio++ wireframe", mesh_length_unit(rMesh));
    auto point = [&](const std::array<double, 3>& rP) {
        return st.Ent("CARTESIAN_POINT(''," + StepStream::Triple(rP.data()) + ")");
    };
    auto direction = [&](const std::array<double, 3>& rD) {
        return st.Ent("DIRECTION(''," + StepStream::Triple(rD.data()) + ")");
    };
    std::vector<std::int64_t> curves;
    std::string name;
    if (has_sketch(rMesh)) {
        name = "sketch";
        for (const SketchWorldPiece& r_p : sketch_world_pieces(rMesh, kSwPrefix)) {
            if (r_p.mType == 1) {
                std::array<double, 3> d{r_p.mEnd[0] - r_p.mStart[0], r_p.mEnd[1] - r_p.mStart[1],
                                        r_p.mEnd[2] - r_p.mStart[2]};
                const double len = std::sqrt(sw_dot(d, d));
                if (!(len > 0.0))
                    continue;
                for (double& r_x : d)
                    r_x /= len;
                const std::int64_t vec =
                    st.Ent("VECTOR(''," + StepStream::Ref(direction(d)) + ",1.)");
                const std::int64_t line = st.Ent("LINE(''," + StepStream::Ref(point(r_p.mStart)) +
                                                 "," + StepStream::Ref(vec) + ")");
                curves.push_back(st.Ent("TRIMMED_CURVE(''," + StepStream::Ref(line) +
                                        ",(PARAMETER_VALUE(0.)),(PARAMETER_VALUE(" +
                                        step_real(len) + ")),.T.,.PARAMETER.)"));
                continue;
            }
            if (r_p.mType != 2 && r_p.mType != 3)
                throw WriteError(std::string(kSwPrefix) + "a sketch piece of unknown type " +
                                 std::to_string(r_p.mType));
            // An arc runs counter-clockwise about its axis: the normal, or its
            // reverse for a clockwise sweep.
            std::array<double, 3> axis = r_p.mNormal;
            if (r_p.mType == 2 && r_p.mSweepDeg < 0.0)
                axis = {-axis[0], -axis[1], -axis[2]};
            const std::int64_t place = st.Ent(
                "AXIS2_PLACEMENT_3D(''," + StepStream::Ref(point(r_p.mCentre)) + "," +
                StepStream::Ref(direction(axis)) + "," + StepStream::Ref(direction(r_p.mE1)) + ")");
            const std::int64_t circle =
                st.Ent("CIRCLE(''," + StepStream::Ref(place) + "," + step_real(r_p.mRadius) + ")");
            if (r_p.mType == 3) {
                curves.push_back(circle);
                continue;
            }
            const std::array<double, 3> y = sw_cross(axis, r_p.mE1);
            const std::array<double, 3> s{r_p.mStart[0] - r_p.mCentre[0],
                                          r_p.mStart[1] - r_p.mCentre[1],
                                          r_p.mStart[2] - r_p.mCentre[2]};
            double t0 = std::atan2(sw_dot(s, y), sw_dot(s, r_p.mE1));
            if (t0 < 0.0)
                t0 += 2.0 * kSwPi;
            const double t1 = t0 + std::fabs(r_p.mSweepDeg) * kSwPi / 180.0;
            curves.push_back(st.Ent("TRIMMED_CURVE(''," + StepStream::Ref(circle) +
                                    ",(PARAMETER_VALUE(" + step_real(t0) + ")),(PARAMETER_VALUE(" +
                                    step_real(t1) + ")),.T.,.PARAMETER.)"));
        }
    } else {
        name = "sections";
        std::vector<std::array<std::int64_t, 2>> segments;
        for (const auto cb : rMesh.CellRange()) {
            if (cell_type_from_name(std::string(cb.Type())) != CellType::Line)
                continue;
            for (std::size_t c = 0; c < cb.NumCells(); ++c)
                segments.push_back({read_int(cb.Conn(), 2 * c), read_int(cb.Conn(), 2 * c + 1)});
        }
        const std::size_t np = rMesh.NumPoints();
        const std::size_t dim = rMesh.PointDim();
        for (const auto& r_s : segments)
            for (const std::int64_t v : r_s)
                if (v < 0 || static_cast<std::size_t>(v) >= np)
                    throw WriteError(std::string(kSwPrefix) + "a line cell names point " +
                                     std::to_string(v) + " of " + std::to_string(np));
        const NDArray& pts = rMesh.Points();
        for (const PolylineChain& r_c : polyline_chains(segments, np)) {
            std::string list;
            for (const std::int64_t v : r_c.mPoints) {
                std::array<double, 3> p{};
                for (std::size_t k = 0; k < dim && k < 3; ++k)
                    p[k] = read_double(pts, static_cast<std::size_t>(v) * dim + k);
                list += (list.empty() ? "" : ",") + StepStream::Ref(point(p));
            }
            curves.push_back(st.Ent("POLYLINE('',(" + list + "))"));
        }
    }
    if (curves.empty())
        throw WriteError(std::string(kSwPrefix) + "no curves to write");
    std::string list;
    for (const std::int64_t c : curves)
        list += (list.empty() ? "" : ",") + StepStream::Ref(c);
    const std::int64_t set = st.Ent("GEOMETRIC_CURVE_SET('" + name + "',(" + list + "))");
    st.FinishMixed({}, {}, {set});
}

}  // namespace detail
}  // namespace meshioplusplus
