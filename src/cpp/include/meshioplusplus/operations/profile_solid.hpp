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

/**
 * @file operations/profile_solid.hpp
 * @brief Profile solids through a CAD provider (roadmap 7.7.4): a closed
 * `fit_curves` sketch extruded along its plane's normal or revolved about an
 * axis, then filleted, chamfered or shelled, written as STEP and tessellated.
 *
 * The design (`doc/profile_solid.md`): the smallest scripted API, a list of
 * features applied in order to one sketch, and no more. The sketch is a line
 * mesh with `fit_curves`' `sketch:*` field data; its pieces (lines, arcs,
 * circles) are chained into closed loops in world coordinates, the outer loop
 * the one of largest area, the others holes. The first feature makes the
 * solid (`Extrude` or `Revolve`); the rest modify it (`Fillet` and `Chamfer`
 * round or bevel every edge, `Shell` hollows it to a wall, removing the face
 * the extrusion ended at). The core chains the loops, and a named provider
 * (`add_profile_solid_provider`; OpenCASCADE's is `"occt"`) builds the solid,
 * writes it and tessellates it. Without a provider the call throws
 * `Unsupported`, naming the add-on. The interactive parametric modeller
 * (sketch constraints, a feature history with rebuilds) is out (roadmap
 * non-goals).
 */

// System includes
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

enum class ProfileFeatureKind : std::int32_t {
    Extrude = 0,  ///< along the sketch's normal by `mDistance` (0: to the scan's extent)
    Revolve = 1,  ///< about `mAxisOrigin`, `mAxisDirection` by `mAngleDeg`
    Fillet = 2,   ///< every edge rounded by `mRadius`
    Chamfer = 3,  ///< every edge bevelled by `mDistance`
    Shell = 4     ///< hollowed to the wall `mDistance`, the end face removed
};

/// `"extrude"`, `"revolve"`, `"fillet"`, `"chamfer"` or `"shell"`.
MESHIOPLUSPLUS_API ProfileFeatureKind profile_feature_from_name(const std::string& rName);

struct ProfileFeature {
    ProfileFeatureKind mKind = ProfileFeatureKind::Extrude;
    double mDistance = 0.0;   ///< extrude length (negative: against the normal); chamfer; wall
    bool mSymmetric = false;  ///< extrude half each way
    double mAngleDeg = 360.0;
    std::array<double, 3> mAxisOrigin{0.0, 0.0, 0.0};
    std::array<double, 3> mAxisDirection{0.0, 0.0, 1.0};
    double mRadius = 0.0;
};

/**
 * @brief Features from their text form, `;`-separated, as the CLIs and the
 * pipeline take them: `extrude:D` (`extrude:D:symmetric`), `revolve:OX,OY,OZ:
 * DX,DY,DZ[:ANGLE]`, `fillet:R`, `chamfer:D`, `shell:W`.
 * @throws std::invalid_argument naming the item it cannot read.
 */
MESHIOPLUSPLUS_API std::vector<ProfileFeature> profile_features_from_text(const std::string& rText);

struct ProfileSolidOptions {
    std::vector<ProfileFeature> mFeatures;
    /// The provider by name; empty takes the first registered.
    std::string mProvider;
    /// The tessellation's chordal deflection; at most 0 takes `1e-3` of the
    /// sketch's (or the scan's) diagonal.
    double mDeflection = 0.0;
};

/// One piece of a chained loop, in world coordinates.
struct ProfilePiece {
    std::int32_t mType = 1;  ///< 1 line, 2 arc (from start to end about the normal), 3 circle
    std::array<double, 3> mStart{}, mEnd{}, mCentre{};
    double mRadius = 0.0;
    double mSweepDeg = 0.0;  ///< an arc's sweep, counter-clockwise about the normal positive
};

/// What a provider is asked to build.
struct ProfileSolidRequest {
    /// The closed loops, the outer one first, each piece's end the next's start.
    std::vector<std::vector<ProfilePiece>> mLoops;
    std::array<double, 3> mOrigin{}, mNormal{};  ///< the sketch plane
    std::vector<ProfileFeature> mFeatures;       ///< extrusions to the extent resolved
    std::string mStepPath;                       ///< empty: no file
    double mDeflection = 0.0;
};

struct ProfileSolidResult {
    Mesh mMesh;  ///< the solid's tessellation
    double mVolume = 0.0;
    double mArea = 0.0;
    std::int64_t mNumFaces = 0;
    bool mValid = false;  ///< the provider's own check of the solid
    std::string mProvider;
    std::vector<std::string> mWarnings;
};

struct ProfileSolidProvider {
    std::string mName;
    std::function<ProfileSolidResult(const ProfileSolidRequest&)> mBuild;
};

/// Register a provider under its name (a second with the name replaces it).
MESHIOPLUSPLUS_API void add_profile_solid_provider(ProfileSolidProvider Provider);
MESHIOPLUSPLUS_API void remove_profile_solid_provider(const std::string& rName) noexcept;
/// The provider by name, or the first when empty; null when none.
MESHIOPLUSPLUS_API std::shared_ptr<const ProfileSolidProvider> profile_solid_provider(
    const std::string& rName = {}) noexcept;
MESHIOPLUSPLUS_API std::vector<std::string> profile_solid_provider_names();

/**
 * @brief Build a profile solid from a sketch (see the file comment); @p pScan,
 * when given, sets an `Extrude` of distance 0 to the scan's extent along the
 * sketch's normal.
 * @throws Unsupported without a provider; std::invalid_argument on a mesh
 * without a sketch, loops that do not close, no feature or a first feature
 * that is not `Extrude` or `Revolve`, or a value out of range.
 */
MESHIOPLUSPLUS_API ProfileSolidResult profile_solid(const Mesh& rSketch,
                                                    const std::string& rStepPath,
                                                    const ProfileSolidOptions& rOptions,
                                                    const Mesh* pScan = nullptr);

}  // namespace meshioplusplus
