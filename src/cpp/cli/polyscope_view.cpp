// SPDX-License-Identifier: MIT
/**
 * @file polyscope_view.cpp
 * @brief Registering a @ref ViewPayload with Polyscope.
 *
 * The only translation unit in the project that includes a Polyscope header,
 * and it is compiled only when `MESHIOPLUSPLUS_WITH_POLYSCOPE=ON`.
 * `view_payload.cpp` — where all the mapping logic lives — is deliberately
 * Polyscope-free so it builds and is tested everywhere. The Python side splits
 * at exactly the same seam (`_to_polyscope_payload` vs `_view_polyscope`).
 *
 * Each structure type has its own quantity vocabulary, and they are **not**
 * uniform: surface meshes say vertex/face, volume meshes vertex/cell, curve
 * networks node/edge, and a point cloud takes no qualifier at all because
 * there is nowhere else its data could live. That is why this is four small
 * explicit functions rather than one template — the method *names* differ, so
 * a template would fight the API rather than use it. (On the Python side the
 * same asymmetry is a runtime `TypeError`, which is why that code carries an
 * explicit per-kind vocabulary instead.)
 */

#include "polyscope_view.hpp"

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/screenshot.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"
#include "polyscope/volume_mesh.h"

#include "view_payload.hpp"

namespace meshioplusplus::cli {

namespace {

/// De-interleave a vector/colour array; the payload keeps them flat so it can
/// stay free of Polyscope's types.
std::vector<std::array<double, 3>> psview_triples(const std::vector<double>& rFlat) {
    std::vector<std::array<double, 3>> out(rFlat.size() / 3);
    for (std::size_t i = 0; i < out.size(); ++i)
        out[i] = {rFlat[i * 3], rFlat[i * 3 + 1], rFlat[i * 3 + 2]};
    return out;
}

/**
 * @brief Headless backends, most preferred first.
 *
 * Asked for **by name** rather than through Polyscope's
 * allow-headless-backends option: not every build has EGL, and when it is
 * absent that route falls through to GLFW and fails without a display. Naming
 * them lets the error say which were tried. Same reasoning, and same list, as
 * `_viewer._HEADLESS_BACKENDS`.
 */
constexpr const char* kHeadlessBackends[] = {"openGL3_egl", "openGL_mock"};

void psview_init(bool Headless) {
    // Polyscope refuses to re-initialize with a *different* backend, so the
    // first initialization in a process wins.
    if (polyscope::isInitialized())
        return;
    if (!Headless) {
        polyscope::init();
        return;
    }
    std::string tried;
    for (const char* backend : kHeadlessBackends) {
        try {
            polyscope::init(backend);
            return;
        } catch (const std::exception& e) {
            if (!tried.empty())
                tried += "; ";
            tried += std::string(backend) + ": " + e.what();
        }
    }
    throw std::runtime_error(
        "screenshot: no headless rendering backend is available (tried " + tried +
        "). Install EGL, or run under xvfb.");
}

/// Pin an explicit colormap range when the payload asked for one — it only
/// does so when the array holds non-finite values.
template <class Q>
void psview_apply_range(Q* pQuantity, const Quantity& rSpec) {
    if (rSpec.mHasRange)
        pQuantity->setMapRange({rSpec.mMin, rSpec.mMax});
}

void psview_add_surface(polyscope::SurfaceMesh* pMesh, const ViewPayload& rPayload) {
    for (const auto& q : rPayload.mQuantities) {
        const bool onVertices = q.mOn == QuantityOn::Vertices;
        switch (q.mKind) {
            case QuantityKind::Scalar:
                if (onVertices)
                    psview_apply_range(pMesh->addVertexScalarQuantity(q.mName, q.mValues), q);
                else
                    psview_apply_range(pMesh->addFaceScalarQuantity(q.mName, q.mValues), q);
                break;
            case QuantityKind::Vector:
                if (onVertices)
                    pMesh->addVertexVectorQuantity(q.mName, psview_triples(q.mValues));
                else
                    pMesh->addFaceVectorQuantity(q.mName, psview_triples(q.mValues));
                break;
            case QuantityKind::Color:
                if (onVertices)
                    pMesh->addVertexColorQuantity(q.mName, psview_triples(q.mValues));
                else
                    pMesh->addFaceColorQuantity(q.mName, psview_triples(q.mValues));
                break;
        }
    }
}

void psview_add_volume(polyscope::VolumeMesh* pMesh, const ViewPayload& rPayload) {
    for (const auto& q : rPayload.mQuantities) {
        const bool onVertices = q.mOn == QuantityOn::Vertices;
        switch (q.mKind) {
            case QuantityKind::Scalar:
                if (onVertices)
                    psview_apply_range(pMesh->addVertexScalarQuantity(q.mName, q.mValues), q);
                else
                    psview_apply_range(pMesh->addCellScalarQuantity(q.mName, q.mValues), q);
                break;
            case QuantityKind::Vector:
                if (onVertices)
                    pMesh->addVertexVectorQuantity(q.mName, psview_triples(q.mValues));
                else
                    pMesh->addCellVectorQuantity(q.mName, psview_triples(q.mValues));
                break;
            case QuantityKind::Color:
                if (onVertices)
                    pMesh->addVertexColorQuantity(q.mName, psview_triples(q.mValues));
                else
                    pMesh->addCellColorQuantity(q.mName, psview_triples(q.mValues));
                break;
        }
    }
}

void psview_add_curve(polyscope::CurveNetwork* pCurve, const ViewPayload& rPayload) {
    for (const auto& q : rPayload.mQuantities) {
        const bool onNodes = q.mOn == QuantityOn::Vertices;
        switch (q.mKind) {
            case QuantityKind::Scalar:
                if (onNodes)
                    psview_apply_range(pCurve->addNodeScalarQuantity(q.mName, q.mValues), q);
                else
                    psview_apply_range(pCurve->addEdgeScalarQuantity(q.mName, q.mValues), q);
                break;
            case QuantityKind::Vector:
                if (onNodes)
                    pCurve->addNodeVectorQuantity(q.mName, psview_triples(q.mValues));
                else
                    pCurve->addEdgeVectorQuantity(q.mName, psview_triples(q.mValues));
                break;
            case QuantityKind::Color:
                if (onNodes)
                    pCurve->addNodeColorQuantity(q.mName, psview_triples(q.mValues));
                else
                    pCurve->addEdgeColorQuantity(q.mName, psview_triples(q.mValues));
                break;
        }
    }
}

void psview_add_points(polyscope::PointCloud* pCloud, const ViewPayload& rPayload) {
    // No qualifier: a point cloud has nowhere else data could live.
    for (const auto& q : rPayload.mQuantities) {
        switch (q.mKind) {
            case QuantityKind::Scalar:
                psview_apply_range(pCloud->addScalarQuantity(q.mName, q.mValues), q);
                break;
            case QuantityKind::Vector:
                pCloud->addVectorQuantity(q.mName, psview_triples(q.mValues));
                break;
            case QuantityKind::Color:
                pCloud->addColorQuantity(q.mName, psview_triples(q.mValues));
                break;
        }
    }
}

/// Register the payload and attach its quantities.
void psview_register(const ViewPayload& rPayload, const std::string& rName) {
    switch (rPayload.mKind) {
        case ViewKind::Points:
            psview_add_points(polyscope::registerPointCloud(rName, rPayload.mVertices),
                              rPayload);
            return;
        case ViewKind::Curve:
            psview_add_curve(
                polyscope::registerCurveNetwork(rName, rPayload.mVertices, rPayload.mEdges),
                rPayload);
            return;
        case ViewKind::Volume: {
            polyscope::VolumeMesh* pMesh = nullptr;
            if (!rPayload.mTets.empty())
                pMesh = polyscope::registerTetMesh(rName, rPayload.mVertices, rPayload.mTets);
            else if (!rPayload.mHexes.empty())
                pMesh = polyscope::registerHexMesh(rName, rPayload.mVertices, rPayload.mHexes);
            else
                // Already padded, in *our* block order -- see ViewPayload.
                pMesh = polyscope::registerVolumeMesh(rName, rPayload.mVertices,
                                                      rPayload.mMixedCells);
            psview_add_volume(pMesh, rPayload);
            return;
        }
        case ViewKind::Surface:
        case ViewKind::Auto:
            psview_add_surface(
                polyscope::registerSurfaceMesh(rName, rPayload.mVertices, rPayload.mFaces),
                rPayload);
            return;
    }
}

/// Fail naming the available quantities, rather than silently ignoring the flag.
void psview_check_color_by(const ViewPayload& rPayload, const std::string& rColorBy) {
    if (rColorBy.empty())
        return;
    std::string available;
    for (const auto& q : rPayload.mQuantities) {
        if (q.mName == rColorBy)
            return;
        if (!available.empty())
            available += ", ";
        available += q.mName;
    }
    throw std::invalid_argument("view: no data array named '" + rColorBy +
                                "' (available: " + (available.empty() ? "none" : available) +
                                ")");
}

/// Emit the payload's audit trail, so nothing lossy happens silently.
void psview_report(const ViewPayload& rPayload) {
    for (const auto& note : rPayload.mNotes)
        std::fprintf(stderr, "warning: view: %s\n", note.c_str());
}

}  // namespace

bool has_polyscope() noexcept { return true; }

void view_mesh(const Mesh& rMesh, ViewKind Kind, const std::string& rColorBy,
               const std::string& rName) {
    const ViewPayload payload = build_view_payload(rMesh, Kind);
    psview_check_color_by(payload, rColorBy);
    psview_report(payload);

    psview_init(/*Headless=*/false);
    psview_register(payload, rName);
    polyscope::show();
}

void screenshot_mesh(const Mesh& rMesh, ViewKind Kind, const std::string& rColorBy,
                     const std::string& rName, const std::string& rOutPath, int Width,
                     int Height, bool Transparent) {
    const ViewPayload payload = build_view_payload(rMesh, Kind);
    psview_check_color_by(payload, rColorBy);
    psview_report(payload);

    psview_init(/*Headless=*/true);
    polyscope::view::setWindowSize(Width, Height);
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    psview_register(payload, rName);
    polyscope::view::resetCameraToHomeView();
    polyscope::screenshot(rOutPath, Transparent);
    polyscope::removeAllStructures();
}

}  // namespace meshioplusplus::cli
