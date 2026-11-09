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
 * @file formats/dxf.hpp
 * @brief DXF R2000 (AC1015) drawings of line meshes and sketches (roadmap
 * 7.6.4).
 *
 * The writer emits an ASCII DXF R2000 file: the HEADER (`$ACADVER`,
 * `$HANDSEED`, `$INSUNITS` from `meshio:length_unit`), the TABLES and BLOCKS a
 * reader expects (the model and paper space block records and layouts), the
 * ENTITIES and an OBJECTS root dictionary. A mesh carrying a sketch
 * (`fit_curves`' `sketch:*` field data) is written as its pieces -- `LINE`,
 * `ARC` and `CIRCLE` in world coordinates, the arcs and circles in the object
 * coordinate system of their sketch plane's normal (the DXF arbitrary-axis
 * algorithm) -- and any other line mesh as one `LINE` per segment; vertex
 * cells are skipped and other cells refused by name. The numbers are written
 * with `%.17g`, so the Python reference (`dxf/_dxf.py`) writes the same bytes.
 *
 * The reader takes `LINE`, `ARC`, `CIRCLE` (64 segments a full turn),
 * `LWPOLYLINE` (with bulges) and `SPLINE` (de Boor, eight samples a control
 * point) from the ENTITIES section and refuses every other entity by name;
 * coincident points are merged, and each segment carries the Int64 cell data
 * `dxf:entity` (its entity's index) and `dxf:type` (1 line, 2 arc, 3 circle,
 * 4 polyline, 5 spline).
 * @see doc/formats/dxf.md
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/** @brief Read a DXF file's lines, arcs, circles, polylines and splines. */
MESHIOPLUSPLUS_API Mesh read_dxf(const std::string& rPath);

/** @brief Write a line mesh or a sketch as ASCII DXF R2000. */
MESHIOPLUSPLUS_API void write_dxf(const std::string& rPath, const Mesh& rMesh);

}  // namespace meshioplusplus
