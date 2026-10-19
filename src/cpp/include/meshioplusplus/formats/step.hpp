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
 * @file formats/step.hpp
 * @brief STEP (ISO 10303-21, AP214) writer: a closed surface as a faceted
 * B-rep solid a CAD tool opens -- the CAD hand-off of the skin-to-CAD work;
 * and, through the optional OpenCASCADE backend, the STEP reader.
 *
 * ### Reading
 *
 * `read_step` needs the OpenCASCADE backend (`MESHIOPLUSPLUS_WITH_OCCT`, LGPL-2.1 with the OCCT exception, bring your
 * own) and throws `Unsupported` naming the flag without it. The file's shape
 * is triangulated by `BRepMesh` at `CadReadOptions::mDeflection` into one
 * `triangle` block with the `cad:*` cell data of `BrepInfo`; the overload
 * taking a `BrepInfo` also returns the shape as NURBS -- every face's surface,
 * its trimming loops and the edges' curves and topology -- which
 * `write_cad_json` writes as Kratos CAD JSON. See doc/formats/step.md.
 *
 * ### Writing
 *
 * Each surface facet becomes an `ADVANCED_FACE` on a `PLANE`, bounded by an
 * `EDGE_LOOP` of `ORIENTED_EDGE`s over `EDGE_CURVE`s on `LINE`s; an edge two
 * facets share is one `EDGE_CURVE`, so the shell is topologically closed with
 * no sewing. Every edge-connected piece of the surface is one
 * `CLOSED_SHELL` inside one `MANIFOLD_SOLID_BREP`, collected in an
 * `ADVANCED_BREP_SHAPE_REPRESENTATION` under the usual AP214 product
 * structure, in millimetres.
 *
 *  - **Input**: `triangle`, `quad` and `polygon` blocks (quadratic surface
 *    cells contribute their corners), or -- when the mesh has volume cells --
 *    the skin of those (`extract_skin`), as the `stl` writer takes it. Lines
 *    and vertices are dropped with a warning.
 *  - **Planarity**: a facet of four or more corners further than `1e-9` of
 *    its size from its own plane is fanned into triangles from its first
 *    corner, so every face lies exactly on its `PLANE`.
 *  - **Closed and oriented**: an open, non-manifold or inconsistently wound
 *    edge throws `WriteError`, naming the first few as point-id pairs (run
 *    `repair` first). A piece wound inwards (negative signed volume) is
 *    written reversed, with a warning, so every solid has positive volume.
 *  - **Bytes**: coordinates are written with `%.17g` (with the decimal point
 *    STEP requires) and the header time stamp is fixed, so a file is
 *    reproducible and the numpy twin writes the same bytes.
 *
 * ### Faces (since v17.9.0)
 *
 * `StepWriteOptions::mFaces` chooses how the facets become faces:
 *
 *  - `Facet`: one face per facet, as described above (the only behaviour
 *    before v17.9.0, and the only one of the numpy twin);
 *  - `Planar`: edge-connected coplanar facets merge into one `PLANE` face,
 *    holes as inner bounds, and the chains between two planes are cut at
 *    their corners into `LINE`s -- a CAD-exported STL's flat faces come back
 *    as the faces they were;
 *  - `Fitted`: the patches `segment` and `fit_primitives` annotated become
 *    `PLANE`, `SPHERICAL_SURFACE`, `CYLINDRICAL_SURFACE`, `CONICAL_SURFACE`
 *    and `TOROIDAL_SURFACE` faces, bounded by `LINE`, `CIRCLE` and polyline
 *    B-spline edges, and the rest coplanar groups as with `Planar`;
 *  - `Auto` (the default): `Fitted` when the mesh carries `segment:patch` and
 *    the `fit:patch_*` table, else `Planar`.
 *
 * No CAD kernel is involved: see `detail/step_brep.hpp` for how the faces,
 * edges and loops are built, and `fit_cad` for the OpenCASCADE route.
 */

// System includes
#include <cstdint>
#include <string>

// Project includes
#include "meshioplusplus/brep.hpp"
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/// How the STEP writer turns facets into faces (since v17.9.0).
enum class StepFaces : std::int32_t {
    /// `Fitted` when the mesh carries fitted patches, else `Planar`.
    Auto = 0,
    /// Edge-connected coplanar facets merge into one planar face.
    Planar = 1,
    /// Fitted patches become analytic faces; the rest as with `Planar`.
    Fitted = 2,
    /// One planar face per facet.
    Facet = 3,
};

/// Options of `write_step` (since v17.9.0).
struct StepWriteOptions {
    /// How facets become faces.
    StepFaces mFaces = StepFaces::Auto;
    /// `Fitted` only: how far, as a length, a patch's points may lie from its
    /// surface for it to be written on it; at most 0 takes three times the
    /// largest accepted `fit:patch_tolerance`, else `1e-3` of the diagonal.
    double mTolerance = 0.0;
};

/**
 * @brief Write a closed surface (or a volume mesh's skin) as a STEP solid,
 * one face per surface (`StepFaces::Auto`).
 * @param rPath the `.step` / `.stp` file to create.
 * @param rMesh the mesh to write.
 * @throws WriteError on an open, non-manifold or inconsistently wound surface,
 *         a degenerate facet, a non-finite coordinate, or no facets at all.
 */
MESHIOPLUSPLUS_API void write_step(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Write a closed surface as a STEP solid with @p rOptions (since v17.9.0).
 * @throws WriteError as `write_step`, and on an option out of range.
 */
MESHIOPLUSPLUS_API void write_step(const std::string& rPath, const Mesh& rMesh,
                                   const StepWriteOptions& rOptions);

/// Read a STEP file as its triangulation (see "Reading" above).
/// @throws Unsupported in a build without OpenCASCADE; ReadError on a file
///         OpenCASCADE cannot read or transfer.
MESHIOPLUSPLUS_API Mesh read_step(const std::string& rPath);

/// Read a STEP file as its triangulation and its B-rep.
MESHIOPLUSPLUS_API Mesh read_step(const std::string& rPath, BrepInfo& rInfo,
                                  const CadReadOptions& rOptions = {});

}  // namespace meshioplusplus
