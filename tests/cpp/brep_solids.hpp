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

// Closed B-rep solids built in code (the C++ suite has no test-data path):
// a box, a capped cylinder and a sphere of two hemispheres, each face an
// exact NURBS from detail/brep_nurbs (planes, a cylinder band closed through
// a seam, hemispheres closed through a seam and a pole), assembled with
// shared edges the way the native STEP writer assembles its solid: one
// `BrepRef` per use of an edge, the seams and poles edges of their face alone.

// System includes
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "../../src/cpp/src/detail/brep_nurbs.hpp"

namespace brep_solids {

namespace mio = meshioplusplus;
using mio::detail::BrepNurbsSurface;
using mio::detail::BrepNurbsTrim;
using mio::detail::PatchVec;

constexpr double kPi = 3.14159265358979323846;

struct SolidFace {
    BrepNurbsSurface mSurface;
    std::vector<std::vector<BrepNurbsTrim>> mLoops;
    bool mSense = true;
};

// Build every face and share the edges `rEdgeCurves` (index = BrepNurbsTrim::mEdge).
inline mio::BrepInfo assemble(const std::vector<SolidFace>& rFaces,
                              const std::vector<mio::BrepCurve>& rEdgeCurves, double Tol = 1e-9) {
    mio::BrepInfo info;
    const std::size_t ns = rEdgeCurves.size();
    std::vector<std::vector<mio::BrepRef>> refs(ns);
    std::vector<mio::BrepCurve> seam_curves;
    std::vector<bool> seam_has;
    for (std::size_t i = 0; i < rFaces.size(); ++i) {
        mio::detail::BrepNurbsFace bf;
        std::string why;
        if (!mio::detail::brep_nurbs_face(rFaces[i].mSurface, rFaces[i].mLoops, rFaces[i].mSense,
                                          Tol, bf, why))
            throw std::runtime_error("brep_solids: face " + std::to_string(i) + ": " + why);
        mio::BrepFace f = bf.mFace;
        f.mId = static_cast<std::int64_t>(i) + 1;
        f.mPatch = static_cast<std::int64_t>(i);
        const std::size_t seam_base = ns + seam_curves.size();
        std::size_t k = 0;
        for (const mio::BrepLoop& l : f.mLoops)
            for (const mio::BrepTrim& t : l.mTrims) {
                const auto& ref = bf.mRefs[k++];
                const std::size_t e = ref.mEdge >= 0
                                          ? static_cast<std::size_t>(ref.mEdge)
                                          : seam_base + static_cast<std::size_t>(ref.mSeam);
                if (e >= refs.size())
                    refs.resize(e + 1);
                mio::BrepRef r;
                r.mId = f.mId;
                r.mTrimIndex = t.mTrimIndex;
                r.mRelativeDirection = ref.mAlong ? 1 : 0;
                refs[e].push_back(r);
            }
        for (std::size_t j = 0; j < bf.mSeams.size(); ++j) {
            seam_curves.push_back(bf.mSeams[j]);
            seam_has.push_back(bf.mSeamHasCurve[j]);
        }
        info.mPatches.push_back(bf.mPatch);
        info.mFaces.push_back(f);
    }
    refs.resize(ns + seam_curves.size());
    for (std::size_t e = 0; e < refs.size(); ++e) {
        mio::BrepEdge be;
        be.mId = static_cast<std::int64_t>(rFaces.size() + e) + 1;
        be.mHasCurve = e < ns || seam_has[e - ns];
        if (be.mHasCurve)
            be.mCurve = e < ns ? rEdgeCurves[e] : seam_curves[e - ns];
        be.mTopology = refs[e];
        info.mEdges.push_back(be);
    }
    return info;
}

inline BrepNurbsTrim line_trim(const PatchVec& rA, const PatchVec& rB, std::int64_t Edge,
                               bool Along) {
    BrepNurbsTrim t;
    t.mPoints = {rA, rB};
    t.mEdge = Edge;
    t.mAlong = Along;
    return t;
}

// The circle of radius R at height Z about the z axis, from angle A0 by Sweep.
inline BrepNurbsTrim circle_trim(double R, double Z, double Sweep, std::int64_t Edge, bool Along,
                                 bool OnPlane, double A0 = 0.0) {
    BrepNurbsTrim t;
    t.mAt = [=](double T) {
        const double a = A0 + Sweep * T;
        return PatchVec{R * std::cos(a), R * std::sin(a), Z};
    };
    t.mEdge = Edge;
    t.mAlong = Along;
    if (OnPlane) {
        t.mCircle = true;
        t.mCentre = {0.0, 0.0, Z};
        t.mRadius = R;
    }
    return t;
}

inline BrepNurbsSurface plane(const PatchVec& rOrigin, const PatchVec& rX, const PatchVec& rY,
                              const PatchVec& rZ) {
    BrepNurbsSurface s;
    s.mKind = 1;
    s.mOrigin = rOrigin;
    s.mX = rX;
    s.mY = rY;
    s.mZ = rZ;
    return s;
}

// The box [0, Sx] x [0, Sy] x [0, Sz]: six planes, twelve shared edges.
inline mio::BrepInfo box(double Sx = 1.0, double Sy = 1.0, double Sz = 1.0) {
    std::array<PatchVec, 8> v;
    for (int i = 0; i < 8; ++i)
        v[static_cast<std::size_t>(i)] = {(i & 1) ? Sx : 0.0, (i & 2) ? Sy : 0.0,
                                          (i & 4) ? Sz : 0.0};
    // Each face's corners counter-clockwise seen from outside, its frame.
    struct Quad {
        std::array<int, 4> mC;
        PatchVec mX, mY, mZ;
    };
    const std::array<Quad, 6> quads = {{{{0, 2, 3, 1}, {1, 0, 0}, {0, -1, 0}, {0, 0, -1}},
                                        {{4, 5, 7, 6}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                                        {{0, 1, 5, 4}, {1, 0, 0}, {0, 0, 1}, {0, -1, 0}},
                                        {{2, 6, 7, 3}, {1, 0, 0}, {0, 0, -1}, {0, 1, 0}},
                                        {{0, 4, 6, 2}, {0, 0, 1}, {0, 1, 0}, {-1, 0, 0}},
                                        {{1, 3, 7, 5}, {0, 1, 0}, {0, 0, 1}, {1, 0, 0}}}};
    std::vector<std::array<int, 2>> edge_ends;
    std::vector<mio::BrepCurve> curves;
    const auto edge = [&](int A, int B) {
        const int lo = std::min(A, B), hi = std::max(A, B);
        for (std::size_t e = 0; e < edge_ends.size(); ++e)
            if (edge_ends[e][0] == lo && edge_ends[e][1] == hi)
                return static_cast<std::int64_t>(e);
        edge_ends.push_back({lo, hi});
        curves.push_back(mio::detail::brep_curve_line(v[static_cast<std::size_t>(lo)],
                                                      v[static_cast<std::size_t>(hi)]));
        return static_cast<std::int64_t>(edge_ends.size() - 1);
    };
    std::vector<SolidFace> faces;
    for (const Quad& q : quads) {
        SolidFace f;
        f.mSurface = plane(v[static_cast<std::size_t>(q.mC[0])], q.mX, q.mY, q.mZ);
        f.mLoops.emplace_back();
        for (int k = 0; k < 4; ++k) {
            const int a = q.mC[static_cast<std::size_t>(k)];
            const int b = q.mC[static_cast<std::size_t>((k + 1) % 4)];
            f.mLoops[0].push_back(line_trim(v[static_cast<std::size_t>(a)],
                                            v[static_cast<std::size_t>(b)], edge(a, b), a < b));
        }
        faces.push_back(f);
    }
    return assemble(faces, curves);
}

// The cylinder of radius R from z = 0 to z = H, capped: a band closed through
// a seam and two discs, sharing the two rims (closed edges from angle 0).
// `TopCapStart` starts the top disc's rim elsewhere, which a model may do: the
// disc's circle is a loop of its own. `CapsFirst` lists the discs before the
// band, so the first use of each rim is a disc's.
inline mio::BrepInfo capped_cylinder(double R = 1.0, double H = 2.0, double TopCapStart = 0.0,
                                     bool CapsFirst = false) {
    std::vector<mio::BrepCurve> curves = {
        mio::detail::brep_curve_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, R, 0.0, 2 * kPi),
        mio::detail::brep_curve_arc({0, 0, H}, {1, 0, 0}, {0, 1, 0}, R, 0.0, 2 * kPi)};
    std::vector<SolidFace> faces(3);
    faces[0].mSurface.mKind = 3;
    faces[0].mSurface.mRadius = R;
    faces[0].mLoops = {{circle_trim(R, 0.0, 2 * kPi, 0, true, false)},
                       {circle_trim(R, H, -2 * kPi, 1, false, false)}};
    faces[1].mSurface = plane({0, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 0, -1});
    faces[1].mLoops = {{circle_trim(R, 0.0, -2 * kPi, 0, false, true)}};
    faces[2].mSurface = plane({0, 0, H}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    faces[2].mLoops = {{circle_trim(R, H, 2 * kPi, 1, true, true, TopCapStart)}};
    if (CapsFirst)
        std::rotate(faces.begin(), faces.begin() + 1, faces.end());
    return assemble(faces, curves);
}

// Two bands of radius R stacked from z = 0 to H and H to 2H, sharing the rim
// at z = H, capped. The upper band's seam is a quarter turn round from the
// lower's, so the shared rim's two uses start at different points and
// neither is a loop of its own (each runs between its band's seams): an
// inconsistency no rotation resolves.
inline mio::BrepInfo stacked_bands(double R = 1.0, double H = 1.0) {
    std::vector<mio::BrepCurve> curves = {
        mio::detail::brep_curve_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, R, 0.0, 2 * kPi),
        mio::detail::brep_curve_arc({0, 0, H}, {1, 0, 0}, {0, 1, 0}, R, 0.0, 2 * kPi),
        mio::detail::brep_curve_arc({0, 0, 2 * H}, {1, 0, 0}, {0, 1, 0}, R, 0.0, 2 * kPi)};
    const double q = 0.5 * kPi;
    std::vector<SolidFace> faces(4);
    faces[0].mSurface.mKind = 3;
    faces[0].mSurface.mRadius = R;
    faces[0].mLoops = {{circle_trim(R, 0.0, 2 * kPi, 0, true, false)},
                       {circle_trim(R, H, -2 * kPi, 1, false, false)}};
    faces[1].mSurface.mKind = 3;
    faces[1].mSurface.mRadius = R;
    faces[1].mSurface.mX = {0, 1, 0};
    faces[1].mSurface.mY = {-1, 0, 0};
    faces[1].mLoops = {{circle_trim(R, H, 2 * kPi, 1, true, false, q)},
                       {circle_trim(R, 2 * H, -2 * kPi, 2, false, false, q)}};
    faces[2].mSurface = plane({0, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 0, -1});
    faces[2].mLoops = {{circle_trim(R, 0.0, -2 * kPi, 0, false, true)}};
    faces[3].mSurface = plane({0, 0, 2 * H}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    faces[3].mLoops = {{circle_trim(R, 2 * H, 2 * kPi, 2, true, true)}};
    return assemble(faces, curves);
}

// The sphere of radius R about the origin: two hemispheres sharing the
// equator, each closed through a seam and a pole.
inline mio::BrepInfo sphere(double R = 1.0) {
    std::vector<mio::BrepCurve> curves = {
        mio::detail::brep_curve_arc({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, R, 0.0, 2 * kPi)};
    std::vector<SolidFace> faces(2);
    for (auto& f : faces) {
        f.mSurface.mKind = 2;
        f.mSurface.mRadius = R;
    }
    faces[0].mLoops = {{circle_trim(R, 0.0, 2 * kPi, 0, true, false)}};
    faces[1].mLoops = {{circle_trim(R, 0.0, -2 * kPi, 0, false, false)}};
    return assemble(faces, curves);
}

// One planar square face [0, S]^2 at z = 0 (normal +z), its four edges.
inline mio::BrepInfo square(double S = 1.0) {
    const std::array<PatchVec, 4> v = {{{0, 0, 0}, {S, 0, 0}, {S, S, 0}, {0, S, 0}}};
    std::vector<mio::BrepCurve> curves;
    SolidFace f;
    f.mSurface = plane(v[0], {1, 0, 0}, {0, 1, 0}, {0, 0, 1});
    f.mLoops.emplace_back();
    for (std::size_t k = 0; k < 4; ++k) {
        curves.push_back(mio::detail::brep_curve_line(v[k], v[(k + 1) % 4]));
        f.mLoops[0].push_back(line_trim(v[k], v[(k + 1) % 4], static_cast<std::int64_t>(k), true));
    }
    return assemble({f}, curves);
}

// The unit box as six untrimmed bilinear patches and no edges or loops: each
// face is its parameter rectangle, oriented outward (a x b along the
// outward normal), and only the trims' images tie the faces together.
inline mio::BrepInfo patch_box() {
    struct Side {
        PatchVec mO, mA, mB;
    };
    const std::array<Side, 6> sides = {{{{0, 0, 0}, {0, 1, 0}, {1, 0, 0}},
                                        {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}},
                                        {{0, 0, 0}, {1, 0, 0}, {0, 0, 1}},
                                        {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},
                                        {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}},
                                        {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}};
    mio::BrepInfo info;
    for (std::size_t f = 0; f < sides.size(); ++f) {
        mio::NurbsPatch p;
        p.mDegrees = {1, 1};
        p.mKnots = {{0, 0, 1, 1}, {0, 0, 1, 1}};
        p.mControlPoints = mio::NDArray(mio::DType::Float64, {4, 3});
        p.mWeights = mio::NDArray(mio::DType::Float64, {4});
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 2; ++i)
                for (int k = 0; k < 3; ++k)
                    p.mControlPoints.As<double>()[3 * (2 * j + i) + k] =
                        sides[f].mO[static_cast<std::size_t>(k)] +
                        i * sides[f].mA[static_cast<std::size_t>(k)] +
                        j * sides[f].mB[static_cast<std::size_t>(k)];
        for (int i = 0; i < 4; ++i)
            p.mWeights.As<double>()[i] = 1.0;
        info.mPatches.push_back(p);
        mio::BrepFace face;
        face.mId = static_cast<std::int64_t>(f) + 1;
        face.mPatch = static_cast<std::int64_t>(f);
        info.mFaces.push_back(face);
    }
    return info;
}

}  // namespace brep_solids
