// SPDX-License-Identifier: MIT
/**
 * @file view_payload.hpp
 * @brief Mapping a mesh onto the arrays a Polyscope structure needs.
 *
 * **This is the C++ twin of `src/python/meshioplusplus/_viewer.py`'s
 * `_to_polyscope_payload` -- KEEP THEM IN SYNC.** Two implementations exist
 * because Python needs one that works with no compiled core (the pure-Python
 * fallback contract), so moving the logic here would relocate the duplication
 * rather than remove it. `tests/cpp/test_view_payload.cpp` and
 * `tests/python/test_viewer.py` assert the same invariants, which is what makes them
 * a parity check rather than two unrelated suites.
 *
 * Deliberately free of any Polyscope header: this half is pure mapping, so it
 * compiles and is tested in the **default** build with no OpenGL, no GLFW and
 * no window. `polyscope_view.cpp` is the half that registers the result and is
 * the only thing guarded by `MESHIOPLUSPLUS_HAS_POLYSCOPE`. The Python side
 * splits at exactly the same seam.
 *
 * Two constraints drive the whole design, and getting either wrong yields a
 * mesh that renders but is silently coloured wrong.
 *
 * **1. A Polyscope volume mesh holds only tetrahedra and hexahedra.**
 * `registerTetMesh`/`registerHexMesh`/`registerVolumeMesh` take nothing else,
 * so wedges, pyramids and the higher-order family go through
 * `convert_cells(Simplexify)` first -- but hexahedra deliberately do not, since
 * decomposing them destroys exactly the element structure a user opens a
 * viewer to look at.
 *
 * **2. `convert_cells` preserves cell data block-for-block;
 * `extract_surface` destroys it.** So the volume path re-reads the converted
 * mesh's own cell data, while the surface-of-a-solid path uses
 * `extract_surface(recordParentIds=true)` and `gather_cell_data_onto_surface`.
 */
#ifndef MESHIOPLUSPLUS_CLI_VIEW_PAYLOAD_HPP
#define MESHIOPLUSPLUS_CLI_VIEW_PAYLOAD_HPP

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus::cli {

/// What kind of Polyscope structure a mesh becomes.
enum class ViewKind {
    Auto,     ///< Pick by the mesh's highest cell dimension.
    Surface,  ///< Boundary of a solid, or the 2D cells as they are.
    Volume,   ///< Requires 3D cells.
    Curve,    ///< 1D cells.
    Points,   ///< The point cloud alone.
};

/// Parse a `--kind` value; throws `std::invalid_argument` on anything else.
ViewKind view_kind_from_name(const std::string& rName);

/// How one data array is rendered.
enum class QuantityKind {
    Scalar,  ///< One value per element.
    Vector,  ///< Three components, drawn as arrows.
    Color,   ///< Three components in 0..1, drawn as RGB.
};

/// Where a quantity lives on its structure.
enum class QuantityOn {
    Vertices,
    Faces,  ///< Surface meshes.
    Cells,  ///< Volume meshes.
    Edges,  ///< Curve networks.
};

/// One renderable data array, already reduced to what Polyscope accepts.
struct Quantity {
    std::string mName;  ///< After collision resolution.
    QuantityKind mKind = QuantityKind::Scalar;
    QuantityOn mOn = QuantityOn::Vertices;
    /// `N` values for a scalar, `3N` interleaved for a vector/colour.
    std::vector<double> mValues;
    /// Explicit colormap range, set only when the array holds non-finite
    /// entries -- Polyscope's autoscale would otherwise be driven by them.
    bool mHasRange = false;
    double mMin = 0.0;
    double mMax = 0.0;
};

/**
 * @brief Everything needed to register one Polyscope structure.
 *
 * Exactly one of `mFaces`, `mEdges`, `{mTets, mHexes}` and `mMixedCells` is
 * populated, selected by `mKind`.
 */
struct ViewPayload {
    ViewKind mKind = ViewKind::Surface;
    /// `(P, 3)`, always three columns whatever the input dimension.
    std::vector<std::array<double, 3>> mVertices;

    /// Ragged, for `Surface` -- Polyscope takes nested vectors directly, so
    /// mixed triangle/quad and polygon blocks need no triangulation.
    std::vector<std::vector<std::uint32_t>> mFaces;
    /// For `Curve`.
    std::vector<std::array<std::uint32_t, 2>> mEdges;
    /// For `Volume` when every 3D block is `tetra`.
    std::vector<std::array<std::uint32_t, 4>> mTets;
    /// For `Volume` when every 3D block is `hexahedron`.
    std::vector<std::array<std::uint32_t, 8>> mHexes;
    /**
     * @brief For `Volume` with both types, padded with `kPadIndex`.
     *
     * Built here rather than by handing Polyscope separate tet/hex lists:
     * `registerTetHexMesh` concatenates tets-then-hexes regardless of input
     * order, so the row order -- and therefore the meaning of every per-cell
     * quantity -- would be Polyscope's rather than ours.
     *
     * Note the sentinel differs from the Python twin's `-1`: the container
     * there is a signed numpy array, here it is `uint32_t`, and Polyscope's own
     * marker is `INVALID_IND_32`. Same meaning, different spelling.
     */
    std::vector<std::array<std::uint32_t, 8>> mMixedCells;

    std::vector<Quantity> mQuantities;
    /// Human-readable log of every lossy or surprising step. Empty when the
    /// mesh was rendered exactly as given.
    std::vector<std::string> mNotes;

    /// How many faces/cells/edges this payload draws.
    std::size_t NumPrimitives() const;
};

/// Padding for unused slots of a mixed-cell row; matches Polyscope's
/// `INVALID_IND_32` without needing its headers here.
inline constexpr std::uint32_t kPadIndex = 0xFFFFFFFFu;

/**
 * @brief Map a mesh onto the arrays a Polyscope structure needs.
 *
 * Pure: registers nothing, touches no window, and does not modify `rMesh`.
 *
 * @param rMesh the mesh to map
 * @param Kind which representation to build
 * @return the payload, including a `mNotes` audit trail of every lossy step
 * @throws std::invalid_argument on `Volume` without 3D cells, or a polyhedron
 *         block under `Volume`
 */
ViewPayload build_view_payload(const Mesh& rMesh, ViewKind Kind = ViewKind::Auto);

}  // namespace meshioplusplus::cli

#endif  // MESHIOPLUSPLUS_CLI_VIEW_PAYLOAD_HPP
