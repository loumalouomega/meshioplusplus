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
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/operations/fit_curves.hpp"
#include "../detail/lsq_fit.hpp"
#include "meshioplusplus/detail/data_ops.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/operations/slice.hpp"

namespace meshioplusplus {

namespace {

using FcV2 = std::array<double, 2>;
constexpr double kFcPi = 3.14159265358979323846;
// An arc or a circle steps at most this far between its points (30 degrees,
// a twelve-sided tessellation); a coarser one is a polygon, its sides lines.
constexpr double kFcMaxArcStep = kFcPi / 6.0 + 1e-9;
constexpr const char* kFcPrefix = "meshio++: fit_curves: ";

// A fitted run of chain points [mI, mJ]: a line (point and unit direction) or
// a circle (centre and radius).
struct FcFit {
    SketchPieceType mType = SketchPieceType::Line;
    std::size_t mI = 0, mJ = 0;
    FcV2 mP{{0.0, 0.0}}, mD{{1.0, 0.0}};
    FcV2 mC{{0.0, 0.0}};
    double mR = 0.0;
    double mResidual = 0.0;
    double mSweep = 0.0;  // radians, signed, along the run
};

double fc_cross(const FcV2& rA, const FcV2& rB) {
    return rA[0] * rB[1] - rA[1] * rB[0];
}
FcV2 fc_sub(const FcV2& rA, const FcV2& rB) {
    return {rA[0] - rB[0], rA[1] - rB[1]};
}
double fc_norm(const FcV2& rA) {
    return std::sqrt(rA[0] * rA[0] + rA[1] * rA[1]);
}

bool fc_line(const std::vector<FcV2>& rP, std::size_t I, std::size_t J, FcFit& rOut) {
    const double n = static_cast<double>(J - I + 1);
    FcV2 c{0.0, 0.0};
    for (std::size_t k = I; k <= J; ++k) {
        c[0] += rP[k][0];
        c[1] += rP[k][1];
    }
    c = {c[0] / n, c[1] / n};
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t k = I; k <= J; ++k) {
        const FcV2 d = fc_sub(rP[k], c);
        sxx += d[0] * d[0];
        sxy += d[0] * d[1];
        syy += d[1] * d[1];
    }
    // The direction of largest spread, oriented from the first point to the last.
    const double theta = 0.5 * std::atan2(2.0 * sxy, sxx - syy);
    FcV2 d{std::cos(theta), std::sin(theta)};
    const FcV2 run = fc_sub(rP[J], rP[I]);
    if (d[0] * run[0] + d[1] * run[1] < 0.0)
        d = {-d[0], -d[1]};
    if (!(sxx + syy > 0.0)) {
        const double len = fc_norm(run);
        if (!(len > 0.0))
            return false;
        d = {run[0] / len, run[1] / len};
    }
    double res = 0.0;
    for (std::size_t k = I; k <= J; ++k)
        res = std::max(res, std::fabs(fc_cross(d, fc_sub(rP[k], c))));
    rOut = FcFit{};
    rOut.mType = SketchPieceType::Line;
    rOut.mI = I;
    rOut.mJ = J;
    rOut.mP = c;
    rOut.mD = d;
    rOut.mResidual = res;
    return true;
}

// The signed angle swept going through the points I..J about rC.
double fc_sweep(const std::vector<FcV2>& rP, std::size_t I, std::size_t J, const FcV2& rC) {
    double total = 0.0;
    for (std::size_t k = I; k < J; ++k) {
        const FcV2 a = fc_sub(rP[k], rC), b = fc_sub(rP[k + 1], rC);
        total += std::atan2(fc_cross(a, b), a[0] * b[0] + a[1] * b[1]);
    }
    return total;
}

// The largest angle, about the centre, between consecutive points I..J: a
// sampled arc steps a little at a time, a polygon's corners a lot (a
// rectangle's four corners lie on one circle, ninety degrees apart).
double fc_max_step(const std::vector<FcV2>& rP, std::size_t I, std::size_t J, const FcV2& rC) {
    double most = 0.0;
    for (std::size_t k = I; k < J; ++k) {
        const FcV2 a = fc_sub(rP[k], rC), b = fc_sub(rP[k + 1], rC);
        most = std::max(most, std::fabs(std::atan2(fc_cross(a, b), a[0] * b[0] + a[1] * b[1])));
    }
    return most;
}

bool fc_circle(const std::vector<FcV2>& rP, std::size_t I, std::size_t J, FcFit& rOut) {
    std::vector<double> xs, ys;
    for (std::size_t k = I; k <= J; ++k) {
        xs.push_back(rP[k][0]);
        ys.push_back(rP[k][1]);
    }
    double cx = 0.0, cy = 0.0, r = 0.0;
    if (!detail::lsq_circle_kasa(xs, ys, cx, cy, r))
        return false;
    detail::lsq_circle_geometric(xs, ys, cx, cy, r);
    double res = 0.0;
    for (std::size_t k = 0; k < xs.size(); ++k)
        res = std::max(res, std::fabs(std::hypot(xs[k] - cx, ys[k] - cy) - r));
    rOut = FcFit{};
    rOut.mType = SketchPieceType::Arc;
    rOut.mI = I;
    rOut.mJ = J;
    rOut.mC = {cx, cy};
    rOut.mR = r;
    rOut.mResidual = res;
    rOut.mSweep = fc_sweep(rP, I, J, rOut.mC);
    return true;
}

struct FcContext {
    const std::vector<FcV2>* mpP = nullptr;
    double mTol = 0.0;
    double mExtent = 1.0;
    std::size_t mMinArcPoints = 4;
    double mMinSweep = 0.0;  // radians
};

// One run as a line or an arc within the tolerance, the line preferred.
bool fc_fit_run(const FcContext& rC, std::size_t I, std::size_t J, FcFit& rOut) {
    const std::vector<FcV2>& p = *rC.mpP;
    FcFit line;
    if (fc_line(p, I, J, line) && line.mResidual <= rC.mTol) {
        rOut = line;
        return true;
    }
    if (J - I + 1 < rC.mMinArcPoints)
        return false;
    // An arc is sampled evenly enough: no chord longer than four times the
    // run's median (a long chord is a straight side the arc would bridge).
    std::vector<double> chords;
    for (std::size_t k = I; k < J; ++k)
        chords.push_back(fc_norm(fc_sub(p[k + 1], p[k])));
    std::vector<double> sorted = chords;
    std::nth_element(sorted.begin(),
                     sorted.begin() + static_cast<std::ptrdiff_t>(sorted.size() / 2), sorted.end());
    const double median = sorted[sorted.size() / 2];
    if (*std::max_element(chords.begin(), chords.end()) > 4.0 * median)
        return false;
    FcFit arc;
    if (fc_circle(p, I, J, arc) && arc.mResidual <= rC.mTol &&
        std::fabs(arc.mSweep) >= rC.mMinSweep && arc.mR < 1e6 * rC.mExtent &&
        fc_max_step(p, I, J, arc.mC) <= kFcMaxArcStep) {
        rOut = arc;
        return true;
    }
    return false;
}

void fc_split(const FcContext& rC, std::size_t I, std::size_t J, std::vector<FcFit>& rOut) {
    FcFit f;
    if (J == I + 1 || fc_fit_run(rC, I, J, f)) {
        if (J == I + 1)
            fc_line(*rC.mpP, I, J, f);
        rOut.push_back(f);
        return;
    }
    const std::vector<FcV2>& p = *rC.mpP;
    const FcV2 chord = fc_sub(p[J], p[I]);
    const double len = fc_norm(chord);
    std::size_t best = I + 1;
    double far = -1.0;
    for (std::size_t k = I + 1; k < J; ++k) {
        const FcV2 d = fc_sub(p[k], p[I]);
        const double dist = len > 0.0 ? std::fabs(fc_cross(chord, d)) / len : fc_norm(d);
        if (dist > far) {
            far = dist;
            best = k;
        }
    }
    fc_split(rC, I, best, rOut);
    fc_split(rC, best, J, rOut);
}

void fc_merge(const FcContext& rC, std::vector<FcFit>& rPieces) {
    bool changed = true;
    while (changed) {
        changed = false;
        for (std::size_t k = 0; k + 1 < rPieces.size(); ++k) {
            FcFit u;
            if (fc_fit_run(rC, rPieces[k].mI, rPieces[k + 1].mJ, u)) {
                rPieces[k] = u;
                rPieces.erase(rPieces.begin() + static_cast<std::ptrdiff_t>(k + 1));
                changed = true;
                break;
            }
        }
    }
}

FcV2 fc_project(const FcFit& rF, const FcV2& rQ) {
    if (rF.mType == SketchPieceType::Line) {
        const FcV2 d = fc_sub(rQ, rF.mP);
        const double t = d[0] * rF.mD[0] + d[1] * rF.mD[1];
        return {rF.mP[0] + t * rF.mD[0], rF.mP[1] + t * rF.mD[1]};
    }
    const FcV2 d = fc_sub(rQ, rF.mC);
    const double len = fc_norm(d);
    if (!(len > 0.0))
        return {rF.mC[0] + rF.mR, rF.mC[1]};
    return {rF.mC[0] + rF.mR * d[0] / len, rF.mC[1] + rF.mR * d[1] / len};
}

// The intersections of two pieces' entities.
std::vector<FcV2> fc_intersections(const FcFit& rA, const FcFit& rB) {
    std::vector<FcV2> out;
    const bool la = rA.mType == SketchPieceType::Line, lb = rB.mType == SketchPieceType::Line;
    if (la && lb) {
        const double den = fc_cross(rA.mD, rB.mD);
        if (std::fabs(den) < 1e-12)
            return out;
        const double t = fc_cross(fc_sub(rB.mP, rA.mP), rB.mD) / den;
        out.push_back({rA.mP[0] + t * rA.mD[0], rA.mP[1] + t * rA.mD[1]});
        return out;
    }
    if (la != lb) {
        const FcFit& l = la ? rA : rB;
        const FcFit& c = la ? rB : rA;
        const FcV2 f = fc_sub(l.mP, c.mC);
        const double b = f[0] * l.mD[0] + f[1] * l.mD[1];
        const double disc = b * b - (f[0] * f[0] + f[1] * f[1] - c.mR * c.mR);
        if (disc < 0.0)
            return out;
        for (const double s : {-1.0, 1.0}) {
            const double t = -b + s * std::sqrt(disc);
            out.push_back({l.mP[0] + t * l.mD[0], l.mP[1] + t * l.mD[1]});
        }
        return out;
    }
    const FcV2 d = fc_sub(rB.mC, rA.mC);
    const double dist = fc_norm(d);
    if (!(dist > 0.0) || dist > rA.mR + rB.mR || dist < std::fabs(rA.mR - rB.mR))
        return out;
    const double a = (rA.mR * rA.mR - rB.mR * rB.mR + dist * dist) / (2.0 * dist);
    const double h = std::sqrt(std::max(0.0, rA.mR * rA.mR - a * a));
    const FcV2 m{rA.mC[0] + a * d[0] / dist, rA.mC[1] + a * d[1] / dist};
    out.push_back({m[0] - h * d[1] / dist, m[1] + h * d[0] / dist});
    out.push_back({m[0] + h * d[1] / dist, m[1] - h * d[0] / dist});
    return out;
}

// Turn a line within `Snap` radians of a multiple of 45 degrees onto it.
void fc_snap(FcFit& rF, const std::vector<FcV2>& rP, double Snap) {
    if (rF.mType != SketchPieceType::Line || !(Snap > 0.0))
        return;
    const double theta = std::atan2(rF.mD[1], rF.mD[0]);
    const double step = kFcPi / 4.0;
    const double k = std::round(theta / step);
    if (std::fabs(theta - k * step) > Snap)
        return;
    const FcV2 a = fc_project(rF, rP[rF.mI]), b = fc_project(rF, rP[rF.mJ]);
    rF.mP = {0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])};
    rF.mD = {std::cos(k * step), std::sin(k * step)};
    double res = 0.0;
    for (std::size_t q = rF.mI; q <= rF.mJ; ++q)
        res = std::max(res, std::fabs(fc_cross(rF.mD, fc_sub(rP[q], rF.mP))));
    rF.mResidual = res;
}

double fc_angle(const FcV2& rC, const FcV2& rQ) {
    return std::atan2(rQ[1] - rC[1], rQ[0] - rC[0]);
}

// The arc's sweep from its start to its end, in the direction it runs.
double fc_end_sweep(const FcFit& rF, const FcV2& rStart, const FcV2& rEnd) {
    double s = fc_angle(rF.mC, rEnd) - fc_angle(rF.mC, rStart);
    if (rF.mSweep >= 0.0) {
        while (s <= 0.0)
            s += 2.0 * kFcPi;
        while (s > 2.0 * kFcPi)
            s -= 2.0 * kFcPi;
    } else {
        while (s >= 0.0)
            s -= 2.0 * kFcPi;
        while (s < -2.0 * kFcPi)
            s += 2.0 * kFcPi;
    }
    return s;
}

// The pieces of one chain of 2-D points (closed: the last equals the first).
std::vector<FcFit> fc_chain(const FcContext& rC, bool Closed) {
    const std::vector<FcV2>& p = *rC.mpP;
    const std::size_t m = p.size() - 1;
    std::vector<FcFit> pieces;
    if (Closed && p.size() >= rC.mMinArcPoints + 1) {
        FcFit circle;
        if (fc_circle(p, 0, m, circle) && circle.mResidual <= rC.mTol &&
            circle.mR < 1e6 * rC.mExtent && fc_max_step(p, 0, m, circle.mC) <= kFcMaxArcStep) {
            circle.mType = SketchPieceType::Circle;
            pieces.push_back(circle);
            return pieces;
        }
    }
    if (Closed && m >= 2) {
        std::size_t far = 1;
        double best = -1.0;
        for (std::size_t k = 1; k < m; ++k) {
            const double d = fc_norm(fc_sub(p[k], p[0]));
            if (d > best) {
                best = d;
                far = k;
            }
        }
        fc_split(rC, 0, far, pieces);
        fc_split(rC, far, m, pieces);
    } else {
        fc_split(rC, 0, m, pieces);
    }
    fc_merge(rC, pieces);
    return pieces;
}

NDArray fc_int64(const std::vector<std::int64_t>& rV) {
    NDArray out = NDArray::Uninit(DType::Int64, {rV.size()});
    std::copy(rV.begin(), rV.end(), out.As<std::int64_t>());
    return out;
}

}  // namespace

const char* sketch_piece_type_name(SketchPieceType Type) {
    switch (Type) {
        case SketchPieceType::Line:
            return "line";
        case SketchPieceType::Arc:
            return "arc";
        case SketchPieceType::Circle:
            return "circle";
    }
    return "line";
}

FitCurvesResult fit_curves(const Mesh& rSections, const FitCurvesOptions& rOptions) {
    if (!(rOptions.mTolerance >= 0.0) || !std::isfinite(rOptions.mTolerance))
        throw std::invalid_argument(std::string(kFcPrefix) + "tolerance must be finite and >= 0");
    if (!(rOptions.mSnapAngleDeg >= 0.0 && rOptions.mSnapAngleDeg < 22.5))
        throw std::invalid_argument(std::string(kFcPrefix) + "snap_angle must lie in [0, 22.5)");
    if (rOptions.mMinArcPoints < 3 || !(rOptions.mMinArcSweepDeg >= 0.0))
        throw std::invalid_argument(std::string(kFcPrefix) +
                                    "min_arc_points must be at least 3 and min_arc_sweep >= 0");
    for (std::size_t b = 0; b < rSections.NumCellBlocks(); ++b)
        if (rSections.Cells(b).Type() != "line")
            throw std::invalid_argument(std::string(kFcPrefix) +
                                        "the input must be a line mesh (slice(chain=True) or "
                                        "sections), not " +
                                        rSections.Cells(b).Type());
    if (rSections.NumCellBlocks() > 0 && !rSections.HasCellData(kSectionChainName))
        throw std::invalid_argument(std::string(kFcPrefix) +
                                    "the input carries no section:chain; chain it first "
                                    "(slice(chain=True) or sections)");
    const std::size_t npts = rSections.NumPoints();
    const std::size_t dim = rSections.PointDim();
    std::vector<detail::Vec3> pts(npts, detail::Vec3{0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < npts; ++i)
        for (std::size_t a = 0; a < dim && a < 3; ++a)
            pts[i][a] = detail::read_double(rSections.Points(), i * dim + a);
    FitCurvesResult result;
    detail::Vec3 lo{0, 0, 0}, hi{0, 0, 0};
    for (std::size_t i = 0; i < npts; ++i)
        for (std::size_t a = 0; a < 3; ++a) {
            lo[a] = i == 0 ? pts[i][a] : std::min(lo[a], pts[i][a]);
            hi[a] = i == 0 ? pts[i][a] : std::max(hi[a], pts[i][a]);
        }
    const double diag = detail::vec3_norm(detail::vec3_sub(hi, lo));
    result.mTolerance = rOptions.mTolerance > 0.0 ? rOptions.mTolerance : 1e-3 * diag;

    // The chains: their segments in cell order, and the offset group of each.
    struct FcChain {
        std::vector<std::int64_t> mPoints;
        std::vector<std::size_t> mSegCell;   // global cell of each segment
        std::vector<std::size_t> mSegStart;  // index in mPoints of its first point
        bool mClosed = false;
        double mOffset = 0.0;
        bool mHasOffset = false;
    };
    std::map<std::int64_t, FcChain> chains;
    std::size_t base = 0;
    const bool offsets = rSections.HasCellData(kSectionOffsetName);
    std::vector<std::int64_t> piece_of_cell;
    for (std::size_t b = 0; b < rSections.NumCellBlocks(); ++b) {
        const auto cells = rSections.Cells(b);
        const NDArray& conn = cells.Conn();
        const NDArray& ids = rSections.CellData(kSectionChainName, b);
        for (std::size_t s = 0; s < cells.NumCells(); ++s) {
            FcChain& r_c = chains[detail::read_int(ids, s)];
            const std::int64_t a = detail::read_int(conn, 2 * s),
                               q = detail::read_int(conn, 2 * s + 1);
            if (r_c.mPoints.empty() || r_c.mPoints.back() != a)
                r_c.mPoints.push_back(a);
            r_c.mSegStart.push_back(r_c.mPoints.size() - 1);
            r_c.mPoints.push_back(q);
            r_c.mSegCell.push_back(base + s);
            if (offsets) {
                r_c.mOffset = detail::read_double(rSections.CellData(kSectionOffsetName, b), s);
                r_c.mHasOffset = true;
            }
        }
        base += cells.NumCells();
    }
    piece_of_cell.assign(base, -1);
    for (auto& [id, r_c] : chains)
        r_c.mClosed = r_c.mPoints.size() > 3 && r_c.mPoints.front() == r_c.mPoints.back();

    // The frames: one per offset group (one for all without offsets).
    std::map<double, std::int64_t> frame_of_offset;
    std::vector<std::vector<std::int64_t>> group_points;
    for (const auto& [id, r_c] : chains) {
        const double key = r_c.mHasOffset ? r_c.mOffset : 0.0;
        auto it = frame_of_offset.find(key);
        if (it == frame_of_offset.end()) {
            it = frame_of_offset.emplace(key, static_cast<std::int64_t>(group_points.size())).first;
            group_points.emplace_back();
        }
        auto& r_g = group_points[static_cast<std::size_t>(it->second)];
        r_g.insert(r_g.end(), r_c.mPoints.begin(), r_c.mPoints.end());
    }
    struct FcFrame {
        detail::Vec3 mO, mE1, mE2, mN;
    };
    std::vector<FcFrame> frames;
    const double given = detail::vec3_norm(rOptions.mNormal);
    for (const auto& r_g : group_points) {
        detail::Vec3 c{0, 0, 0};
        for (const std::int64_t p : r_g)
            c = detail::vec3_add(c, pts[static_cast<std::size_t>(p)]);
        c = detail::vec3_scale(c, 1.0 / static_cast<double>(std::max<std::size_t>(1, r_g.size())));
        detail::Vec3 n{0, 0, 1};
        if (given > 0.0) {
            n = detail::vec3_scale(rOptions.mNormal, 1.0 / given);
        } else {
            std::array<double, 9> cov{};
            for (const std::int64_t p : r_g) {
                const detail::Vec3 d = detail::vec3_sub(pts[static_cast<std::size_t>(p)], c);
                for (std::size_t i = 0; i < 3; ++i)
                    for (std::size_t j = 0; j < 3; ++j)
                        cov[3 * i + j] += d[i] * d[j];
            }
            std::array<double, 3> w{};
            std::array<double, 9> v{};
            detail::lsq_sym3_eigen(cov.data(), w.data(), v.data());
            n = {v[0], v[1], v[2]};
        }
        FcFrame f;
        f.mN = n;
        f.mO = detail::vec3_scale(n, detail::vec3_dot(n, c));
        detail::lsq_frame(n.data(), f.mE1.data(), f.mE2.data());
        frames.push_back(f);
        std::array<double, 16> m{};
        for (std::size_t r = 0; r < 3; ++r) {
            m[4 * r + 0] = f.mE1[r];
            m[4 * r + 1] = f.mE2[r];
            m[4 * r + 2] = f.mN[r];
            m[4 * r + 3] = f.mO[r];
        }
        m[15] = 1.0;
        result.mFrames.push_back(m);
    }

    FcContext ctx;
    ctx.mTol = result.mTolerance;
    ctx.mExtent = std::max(diag, 1e-300);
    ctx.mMinArcPoints = static_cast<std::size_t>(rOptions.mMinArcPoints);
    ctx.mMinSweep = rOptions.mMinArcSweepDeg * kFcPi / 180.0;
    const double snap = rOptions.mSnapAngleDeg * kFcPi / 180.0;
    double planar = 0.0;
    for (const auto& [id, r_c] : chains) {
        const std::int64_t fi = frame_of_offset[r_c.mHasOffset ? r_c.mOffset : 0.0];
        const FcFrame& f = frames[static_cast<std::size_t>(fi)];
        std::vector<FcV2> p2(r_c.mPoints.size());
        for (std::size_t k = 0; k < p2.size(); ++k) {
            const detail::Vec3 d =
                detail::vec3_sub(pts[static_cast<std::size_t>(r_c.mPoints[k])], f.mO);
            p2[k] = {detail::vec3_dot(d, f.mE1), detail::vec3_dot(d, f.mE2)};
            planar = std::max(planar, std::fabs(detail::vec3_dot(d, f.mN)));
        }
        if (p2.size() < 2)
            continue;
        // Points on the straight segment between their neighbours (a chord's
        // midpoint where a triangle's diagonal was cut) carry no shape: fit
        // the others, and map their runs back.
        std::vector<std::size_t> keep;
        const double flat = 1e-9 * ctx.mExtent;
        for (std::size_t k = 0; k < p2.size(); ++k) {
            const bool end = k == 0 || k + 1 == p2.size();
            if (!end) {
                const FcV2& a = p2[keep.back()];
                const FcV2& c = p2[k + 1];
                const FcV2 ac = fc_sub(c, a), ab = fc_sub(p2[k], a);
                const double len = fc_norm(ac);
                const double t = len > 0.0 ? (ab[0] * ac[0] + ab[1] * ac[1]) / (len * len) : -1.0;
                if (len > 0.0 && t > 0.0 && t < 1.0 && std::fabs(fc_cross(ac, ab)) / len <= flat)
                    continue;
            }
            keep.push_back(k);
        }
        std::vector<FcV2> q2(keep.size());
        for (std::size_t k = 0; k < keep.size(); ++k)
            q2[k] = p2[keep[k]];
        const bool closed = r_c.mClosed && q2.size() > 3;
        ctx.mpP = &q2;
        std::vector<FcFit> pieces = fc_chain(ctx, closed);
        if (closed && pieces.size() > 1) {
            // The seam: start the loop at the last piece's start and fit again;
            // keep that when it needs fewer pieces.
            const std::size_t m = q2.size() - 1, r = pieces.back().mI;
            std::vector<FcV2> rot(m + 1);
            std::vector<std::size_t> rot_keep(m + 1);
            for (std::size_t t = 0; t <= m; ++t) {
                rot[t] = q2[(r + t) % m];
                rot_keep[t] = keep[(r + t) % m];
            }
            ctx.mpP = &rot;
            std::vector<FcFit> again = fc_chain(ctx, true);
            if (again.size() < pieces.size()) {
                pieces = std::move(again);
                q2 = std::move(rot);
                keep = std::move(rot_keep);
            }
            ctx.mpP = &q2;
        }
        // The chain points a run spans, in order (a loop's run may wrap).
        const std::size_t morig = p2.size() - 1;
        std::vector<std::int64_t> seg_at(p2.size(), -1);
        for (std::size_t s = 0; s < r_c.mSegStart.size(); ++s)
            seg_at[r_c.mSegStart[s]] = static_cast<std::int64_t>(s);
        const auto span = [&](const FcFit& rF) {
            std::vector<std::size_t> out{keep[rF.mI]};
            for (std::size_t k = rF.mI; k < rF.mJ; ++k) {
                std::size_t a = keep[k], b = keep[k + 1];
                if (b <= a && r_c.mClosed)
                    b += morig;
                for (std::size_t q = a + 1; q <= b; ++q)
                    out.push_back(r_c.mClosed ? q % morig : q);
            }
            return out;
        };
        for (FcFit& r_f : pieces) {
            // Residuals over every chain point the run spans.
            for (const std::size_t q : span(r_f)) {
                const double d = r_f.mType == SketchPieceType::Line
                                     ? std::fabs(fc_cross(r_f.mD, fc_sub(p2[q], r_f.mP)))
                                     : std::fabs(fc_norm(fc_sub(p2[q], r_f.mC)) - r_f.mR);
                r_f.mResidual = std::max(r_f.mResidual, d);
            }
        }
        for (FcFit& r_f : pieces)
            fc_snap(r_f, q2, snap);
        // The ends: each run's end points projected, then the shared corners.
        std::vector<FcV2> starts(pieces.size()), ends(pieces.size());
        for (std::size_t k = 0; k < pieces.size(); ++k) {
            starts[k] = fc_project(pieces[k], q2[pieces[k].mI]);
            ends[k] = fc_project(pieces[k], q2[pieces[k].mJ]);
        }
        const std::size_t joints = closed ? pieces.size() : pieces.size() - 1;
        if (pieces.size() > 1)
            for (std::size_t k = 0; k < joints; ++k) {
                const std::size_t q = (k + 1) % pieces.size();
                const FcV2& joint = q2[pieces[k].mJ];
                const std::vector<FcV2> cand = fc_intersections(pieces[k], pieces[q]);
                double best = 10.0 * ctx.mTol;
                bool found = false;
                FcV2 at{0.0, 0.0};
                for (const FcV2& r_x : cand) {
                    const double d = fc_norm(fc_sub(r_x, joint));
                    if (d <= best) {
                        best = d;
                        at = r_x;
                        found = true;
                    }
                }
                if (found) {
                    ends[k] = at;
                    starts[q] = at;
                }
            }
        for (std::size_t k = 0; k < pieces.size(); ++k) {
            const FcFit& r_f = pieces[k];
            const std::vector<std::size_t> pts_of = span(r_f);
            SketchPiece sp;
            sp.mType = r_f.mType;
            sp.mStart = starts[k];
            sp.mEnd = ends[k];
            sp.mResidual = r_f.mResidual;
            sp.mNumPoints = static_cast<std::int64_t>(pts_of.size());
            sp.mChain = id;
            sp.mFrame = fi;
            if (r_f.mType != SketchPieceType::Line) {
                sp.mCentre = r_f.mC;
                sp.mRadius = r_f.mR;
                sp.mSweepDeg = r_f.mType == SketchPieceType::Circle
                                   ? 360.0
                                   : fc_end_sweep(r_f, starts[k], ends[k]) * 180.0 / kFcPi;
            }
            const auto index = static_cast<std::int64_t>(result.mPieces.size());
            for (std::size_t t = 0; t + 1 < pts_of.size(); ++t) {
                const std::int64_t seg = seg_at[pts_of[t]];
                if (seg >= 0)
                    piece_of_cell[r_c.mSegCell[static_cast<std::size_t>(seg)]] = index;
            }
            result.mPieces.push_back(sp);
            result.mMaxResidual = std::max(result.mMaxResidual, sp.mResidual);
            result.mNumLines += sp.mType == SketchPieceType::Line ? 1 : 0;
            result.mNumArcs += sp.mType == SketchPieceType::Arc ? 1 : 0;
            result.mNumCircles += sp.mType == SketchPieceType::Circle ? 1 : 0;
        }
    }
    if (planar > result.mTolerance)
        result.mWarnings.push_back("a chain lies up to " + std::to_string(planar) +
                                   " off its sketch plane; its pieces fit its projection");

    result.mMesh = detail::clone_mesh(rSections);
    std::vector<NDArray> blocks;
    std::size_t at = 0;
    for (std::size_t b = 0; b < rSections.NumCellBlocks(); ++b) {
        const std::size_t n = rSections.Cells(b).NumCells();
        blocks.push_back(fc_int64(std::vector<std::int64_t>(
            piece_of_cell.begin() + static_cast<std::ptrdiff_t>(at),
            piece_of_cell.begin() + static_cast<std::ptrdiff_t>(at + n))));
        at += n;
    }
    if (!blocks.empty())
        result.mMesh.AddCellData(kSketchPieceName, std::move(blocks));
    const std::size_t np = result.mPieces.size();
    std::vector<std::int64_t> types, chain_ids, frame_ids;
    NDArray params = NDArray::Uninit(DType::Float64, {np, std::size_t{10}});
    for (std::size_t k = 0; k < np; ++k) {
        const SketchPiece& r_p = result.mPieces[k];
        types.push_back(static_cast<std::int64_t>(r_p.mType));
        chain_ids.push_back(r_p.mChain);
        frame_ids.push_back(r_p.mFrame);
        const double row[10] = {r_p.mStart[0],  r_p.mStart[1],
                                r_p.mEnd[0],    r_p.mEnd[1],
                                r_p.mCentre[0], r_p.mCentre[1],
                                r_p.mRadius,    r_p.mSweepDeg,
                                r_p.mResidual,  static_cast<double>(r_p.mNumPoints)};
        std::copy(row, row + 10, params.As<double>() + 10 * k);
    }
    NDArray frm = NDArray::Uninit(DType::Float64, {result.mFrames.size(), std::size_t{16}});
    for (std::size_t k = 0; k < result.mFrames.size(); ++k)
        std::copy(result.mFrames[k].begin(), result.mFrames[k].end(), frm.As<double>() + 16 * k);
    result.mMesh.AddFieldData(kSketchTypeName, fc_int64(types));
    result.mMesh.AddFieldData(kSketchParamsName, std::move(params));
    result.mMesh.AddFieldData(kSketchChainName, fc_int64(chain_ids));
    result.mMesh.AddFieldData(kSketchFrameName, fc_int64(frame_ids));
    result.mMesh.AddFieldData(kSketchFramesName, std::move(frm));
    return result;
}

}  // namespace meshioplusplus
