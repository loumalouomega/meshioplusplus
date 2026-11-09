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
// Profile solids through a provider (roadmap 7.7.4). See
// operations/profile_solid.hpp.

// System includes
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

// Project includes
#include "meshioplusplus/operations/profile_solid.hpp"
#include "../detail/bounding_box_utils.hpp"
#include "../detail/sketch_pieces.hpp"
#include "meshioplusplus/detail/fast_number.hpp"
#include "meshioplusplus/detail/geometry.hpp"
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/providers.hpp"

namespace meshioplusplus {

namespace {

using detail::Vec3;

[[noreturn]] void ps_fail(const std::string& rWhy) {
    throw std::invalid_argument("meshio++: profile_solid: " + rWhy);
}

double ps_diagonal(const Mesh& rMesh) {
    const NDArray& p = rMesh.Points();
    std::vector<Vec3> pts(rMesh.NumPoints());
    for (std::size_t i = 0; i < pts.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            pts[i][k] = detail::read_double(p, 3 * i + k);
    return pts.empty()
               ? 0.0
               : detail::BoundingBoxUtils::Diagonal(detail::BoundingBoxUtils::FromPoints(pts));
}

// The area a loop encloses, in its plane (signed by the normal).
double ps_loop_area(const std::vector<ProfilePiece>& rLoop, const Vec3& rNormal) {
    Vec3 sum{0.0, 0.0, 0.0};
    const auto add = [&](const Vec3& rA, const Vec3& rB) {
        sum = detail::vec3_add(sum, detail::vec3_cross(rA, rB));
    };
    for (const ProfilePiece& r_p : rLoop) {
        if (r_p.mType == 3)
            return std::acos(-1.0) * r_p.mRadius * r_p.mRadius;
        if (r_p.mType == 1) {
            add(r_p.mStart, r_p.mEnd);
            continue;
        }
        // An arc, by 32 chords.
        const Vec3 u = detail::vec3_sub(r_p.mStart, r_p.mCentre);
        const Vec3 w = detail::vec3_cross(rNormal, u);
        Vec3 prev = r_p.mStart;
        for (int k = 1; k <= 32; ++k) {
            const double t = r_p.mSweepDeg * std::acos(-1.0) / 180.0 * k / 32.0;
            const Vec3 q =
                detail::vec3_add(r_p.mCentre, detail::vec3_add(detail::vec3_scale(u, std::cos(t)),
                                                               detail::vec3_scale(w, std::sin(t))));
            add(prev, q);
            prev = q;
        }
    }
    return 0.5 * detail::vec3_dot(sum, rNormal);
}

}  // namespace

ProfileFeatureKind profile_feature_from_name(const std::string& rName) {
    if (rName == "extrude")
        return ProfileFeatureKind::Extrude;
    if (rName == "revolve")
        return ProfileFeatureKind::Revolve;
    if (rName == "fillet")
        return ProfileFeatureKind::Fillet;
    if (rName == "chamfer")
        return ProfileFeatureKind::Chamfer;
    if (rName == "shell")
        return ProfileFeatureKind::Shell;
    ps_fail("a feature is extrude, revolve, fillet, chamfer or shell, not '" + rName + "'");
}

std::vector<ProfileFeature> profile_features_from_text(const std::string& rText) {
    const auto split = [](const std::string& rS, char Sep) {
        std::vector<std::string> out;
        std::size_t at = 0;
        while (at <= rS.size()) {
            const std::size_t end = std::min(rS.find(Sep, at), rS.size());
            out.push_back(rS.substr(at, end - at));
            at = end + 1;
        }
        return out;
    };
    std::vector<ProfileFeature> out;
    for (const std::string& r_item : split(rText, ';')) {
        if (r_item.empty())
            continue;
        const std::vector<std::string> parts = split(r_item, ':');
        const auto number = [&](const std::string& rS) {
            const char* end = nullptr;
            const double v = detail::parse_double(rS.c_str(), end);
            if (rS.empty() || end != rS.c_str() + rS.size())
                ps_fail("cannot read '" + rS + "' in the feature '" + r_item + "'");
            return v;
        };
        const auto vec3 = [&](const std::string& rS) {
            const std::vector<std::string> c = split(rS, ',');
            if (c.size() != 3)
                ps_fail("'" + rS + "' in the feature '" + r_item + "' is not x,y,z");
            return std::array<double, 3>{number(c[0]), number(c[1]), number(c[2])};
        };
        ProfileFeature f;
        f.mKind = profile_feature_from_name(parts[0]);
        switch (f.mKind) {
            case ProfileFeatureKind::Extrude:
                if (parts.size() < 2 || parts.size() > 3 ||
                    (parts.size() == 3 && parts[2] != "symmetric"))
                    ps_fail("an extrusion is extrude:D or extrude:D:symmetric, not '" + r_item +
                            "'");
                f.mDistance = number(parts[1]);
                f.mSymmetric = parts.size() == 3;
                break;
            case ProfileFeatureKind::Revolve:
                if (parts.size() < 3 || parts.size() > 4)
                    ps_fail("a revolution is revolve:OX,OY,OZ:DX,DY,DZ[:ANGLE], not '" + r_item +
                            "'");
                f.mAxisOrigin = vec3(parts[1]);
                f.mAxisDirection = vec3(parts[2]);
                if (parts.size() == 4)
                    f.mAngleDeg = number(parts[3]);
                break;
            case ProfileFeatureKind::Fillet:
                if (parts.size() != 2)
                    ps_fail("a fillet is fillet:R, not '" + r_item + "'");
                f.mRadius = number(parts[1]);
                break;
            case ProfileFeatureKind::Chamfer:
            case ProfileFeatureKind::Shell:
                if (parts.size() != 2)
                    ps_fail("'" + r_item + "' takes one distance");
                f.mDistance = number(parts[1]);
                break;
        }
        out.push_back(f);
    }
    return out;
}

ProfileSolidResult profile_solid(const Mesh& rSketch, const std::string& rStepPath,
                                 const ProfileSolidOptions& rOptions, const Mesh* pScan) {
    if (!detail::has_sketch(rSketch))
        ps_fail("the mesh carries no sketch (run fit_curves on a chained section first)");
    const std::vector<ProfileFeature>& feats = rOptions.mFeatures;
    if (feats.empty() || (feats[0].mKind != ProfileFeatureKind::Extrude &&
                          feats[0].mKind != ProfileFeatureKind::Revolve))
        ps_fail("the first feature is extrude or revolve");
    for (std::size_t k = 0; k < feats.size(); ++k) {
        const ProfileFeature& f = feats[k];
        if (k > 0 &&
            (f.mKind == ProfileFeatureKind::Extrude || f.mKind == ProfileFeatureKind::Revolve))
            ps_fail("only the first feature makes the solid");
        if (!std::isfinite(f.mDistance) || !std::isfinite(f.mRadius) || !std::isfinite(f.mAngleDeg))
            ps_fail("a feature's value is not finite");
        if (f.mKind == ProfileFeatureKind::Revolve &&
            (!(f.mAngleDeg > 0.0) || f.mAngleDeg > 360.0 ||
             detail::vec3_norm(
                 Vec3{f.mAxisDirection[0], f.mAxisDirection[1], f.mAxisDirection[2]}) == 0.0))
            ps_fail("a revolve needs an angle in (0, 360] and an axis direction");
        if (f.mKind == ProfileFeatureKind::Fillet && !(f.mRadius > 0.0))
            ps_fail("a fillet needs a positive radius");
        if ((f.mKind == ProfileFeatureKind::Chamfer || f.mKind == ProfileFeatureKind::Shell) &&
            !(f.mDistance > 0.0))
            ps_fail("a chamfer or a shell needs a positive distance");
    }

    // The pieces, chained into closed loops.
    const std::vector<detail::SketchWorldPiece> pieces =
        detail::sketch_world_pieces(rSketch, "meshio++: profile_solid: ");
    if (pieces.empty())
        ps_fail("the sketch has no piece");
    ProfileSolidRequest req;
    req.mNormal = pieces[0].mNormal;
    req.mOrigin = pieces[0].mOrigin;
    const double diag = ps_diagonal(rSketch);
    const double tol = 1e-6 * std::max(diag, 1.0);
    std::vector<char> used(pieces.size(), 0);
    const auto as_piece = [](const detail::SketchWorldPiece& rW, bool Reverse) {
        ProfilePiece p;
        p.mType = static_cast<std::int32_t>(rW.mType);
        p.mStart = Reverse ? rW.mEnd : rW.mStart;
        p.mEnd = Reverse ? rW.mStart : rW.mEnd;
        p.mCentre = rW.mCentre;
        p.mRadius = rW.mRadius;
        p.mSweepDeg = Reverse ? -rW.mSweepDeg : rW.mSweepDeg;
        return p;
    };
    const Vec3 n{req.mNormal[0], req.mNormal[1], req.mNormal[2]};
    for (std::size_t s = 0; s < pieces.size(); ++s) {
        if (used[s])
            continue;
        used[s] = 1;
        std::vector<ProfilePiece> loop{as_piece(pieces[s], false)};
        if (pieces[s].mType == 3) {
            req.mLoops.push_back(std::move(loop));
            continue;
        }
        const auto close = [&](const std::array<double, 3>& rA, const std::array<double, 3>& rB) {
            return detail::vec3_norm(detail::vec3_sub(Vec3{rA[0], rA[1], rA[2]},
                                                      Vec3{rB[0], rB[1], rB[2]})) <= tol;
        };
        while (!close(loop.back().mEnd, loop.front().mStart)) {
            bool found = false;
            for (std::size_t q = 0; q < pieces.size() && !found; ++q) {
                if (used[q] || pieces[q].mType == 3)
                    continue;
                if (close(pieces[q].mStart, loop.back().mEnd)) {
                    loop.push_back(as_piece(pieces[q], false));
                    used[q] = found = true;
                } else if (close(pieces[q].mEnd, loop.back().mEnd)) {
                    loop.push_back(as_piece(pieces[q], true));
                    used[q] = found = true;
                }
            }
            if (!found)
                ps_fail("the sketch's pieces do not close into loops");
        }
        req.mLoops.push_back(std::move(loop));
    }
    // The outer loop (the largest) first, wound counter-clockwise about the
    // normal; the holes clockwise.
    std::vector<double> area(req.mLoops.size());
    for (std::size_t l = 0; l < req.mLoops.size(); ++l)
        area[l] = ps_loop_area(req.mLoops[l], n);
    std::size_t outer = 0;
    for (std::size_t l = 1; l < area.size(); ++l)
        if (std::fabs(area[l]) > std::fabs(area[outer]))
            outer = l;
    std::swap(req.mLoops[0], req.mLoops[outer]);
    std::swap(area[0], area[outer]);
    for (std::size_t l = 0; l < req.mLoops.size(); ++l) {
        const bool want_positive = l == 0;
        if ((area[l] > 0.0) != want_positive && req.mLoops[l][0].mType != 3) {
            std::vector<ProfilePiece>& r_loop = req.mLoops[l];
            std::reverse(r_loop.begin(), r_loop.end());
            for (ProfilePiece& r_p : r_loop) {
                std::swap(r_p.mStart, r_p.mEnd);
                r_p.mSweepDeg = -r_p.mSweepDeg;
            }
        }
    }

    // An extrusion of 0 runs to the scan's extent along the normal.
    req.mFeatures = feats;
    ProfileFeature& first = req.mFeatures[0];
    if (first.mKind == ProfileFeatureKind::Extrude && first.mDistance == 0.0) {
        if (!pScan)
            ps_fail("an extrusion of 0 runs to the scan's extent: pass the scan");
        const NDArray& p = pScan->Points();
        double lo = std::numeric_limits<double>::infinity(), hi = -lo;
        for (std::size_t i = 0; i < pScan->NumPoints(); ++i) {
            const Vec3 q{detail::read_double(p, 3 * i) - req.mOrigin[0],
                         detail::read_double(p, 3 * i + 1) - req.mOrigin[1],
                         detail::read_double(p, 3 * i + 2) - req.mOrigin[2]};
            const double h = detail::vec3_dot(q, n);
            lo = std::min(lo, h);
            hi = std::max(hi, h);
        }
        first.mDistance = std::fabs(hi) >= std::fabs(lo) ? hi : lo;
        if (!(std::fabs(first.mDistance) > tol))
            ps_fail("the scan does not extend off the sketch's plane");
    }
    req.mStepPath = rStepPath;
    req.mDeflection = rOptions.mDeflection > 0.0 ? rOptions.mDeflection
                                                 : 1e-3 * (pScan ? ps_diagonal(*pScan) : diag);
    const std::shared_ptr<const ProfileSolidProvider> provider =
        profile_solid_provider(rOptions.mProvider);
    if (!provider || !provider->mBuild)
        throw Unsupported(
            "meshio++: profile_solid: no profile-solid provider" +
            (rOptions.mProvider.empty() ? std::string() : " named '" + rOptions.mProvider + "'") +
            " is registered; it needs OpenCASCADE (a build with MESHIOPLUSPLUS_WITH_OCCT, or the "
            "meshioplusplus-occt package)");
    ProfileSolidResult result = provider->mBuild(req);
    result.mProvider = provider->mName;
    return result;
}

}  // namespace meshioplusplus
