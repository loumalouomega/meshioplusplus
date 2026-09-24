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
 * @file libmesh.hpp
 * @brief libMesh `.xda` (ASCII) / `.xdr` (XDR binary) mesh reader.
 *
 * libMesh's native mesh file (`XdrIO`), also what MOOSE's `--mesh-only` and
 * checkpoint meshes use. Both encodings carry one value stream: a version
 * string (`libMesh-0.7.0+` ... `libMesh-1.8.0`), the element and node counts,
 * four "inline or not" flags (boundary conditions, subdomain, processor and
 * p-level ids), per-field integer sizes (0.9.2+), subdomain names, one
 * connectivity block per refinement level, the coordinates, then the side
 * sets, node sets (0.9.2+) and edge and shell-face sets (1.1.0+). `.xdr` is
 * big-endian XDR: 4-byte integers and lengths, 8-byte header integers from
 * 1.3.0 on, strings padded to 4 bytes.
 *
 * The mesh:
 *  - Only *active* elements (leaves of the refinement tree) become cells; with
 *    refinement levels present, `libmesh:level` is kept as cell data.
 *  - The subdomain id is the `libmesh:subdomain` cell data and a `Cell` region
 *    per subdomain, named by the file's subdomain map (else
 *    `subdomain_<id>`, tag = id). An inline p-level is `libmesh:p_level`.
 *  - Side sets become `Side` regions (named by the sideset map, else
 *    `boundary_<id>`). A side of a refined element is carried down to the
 *    active descendants whose sides lie on it. Node sets become `Point`
 *    regions (`nodeset_<id>` when unnamed). Edge and shell-face sets are
 *    skipped with a warning.
 *  - HEX20/HEX27/PRISM15/PRISM18 are reordered through the `"libmesh"` tables
 *    of `detail/node_order.hpp`; the other shapes are already in meshio++'s
 *    order. TET14, PRISM20/21 and PYRAMID18 keep the nodes meshio++ has a type
 *    for (tetra10, wedge18, pyramid14), with a warning; shell and subdivision
 *    variants read as their base shape, NODEELEM as a vertex. Infinite
 *    elements are skipped with a warning.
 *  - Node ids that no element uses (libMesh writes their coordinates as NaN)
 *    are dropped; the original ids are then kept as `libmesh:id` point data.
 *
 * Legacy pre-`libMesh` files (`DEAL 003`, `LIBM 0`) are refused, as libMesh
 * itself does. Compressed `.xda.gz`/`.xdr.bz2` files must be decompressed
 * first. See doc/formats/libmesh.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Read a libMesh `.xda` or `.xdr` mesh.
 *
 * The encoding comes from the content, not the extension: a file whose first
 * four bytes are an XDR string length is read as XDR.
 * @param rPath filesystem path to read
 * @return the mesh, with subdomains, side sets and node sets as regions
 * @throws ReadError if the file can't be read, is truncated, is a legacy
 *         pre-`libMesh` file, names an unknown element type or node, or holds
 *         polygon/polyhedron elements
 */
MESHIOPLUSPLUS_API Mesh read_libmesh(const std::string& rPath);

}  // namespace meshioplusplus
