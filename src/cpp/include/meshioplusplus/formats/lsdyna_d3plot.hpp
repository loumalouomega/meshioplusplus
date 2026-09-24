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
 * @file lsdyna_d3plot.hpp
 * @brief LS-DYNA binary state database (the `d3plot` family) reader.
 *
 * The file LS-DYNA writes for post-processing (LS-PrePost, lasso-python and
 * most public crash datasets start from it). It is word-addressed -- 4-byte
 * words, or 8 in a double-precision run, in either byte order, sniffed from
 * the control block -- and split across `d3plot`, `d3plot01`, `d3plot02`...:
 * the family is one stream of states, none of which straddles a file. Open
 * the base file; a numbered member is refused with a message naming it.
 *
 * The mesh comes from the geometry section: the initial coordinates (the user
 * node ids as `point_data["lsdyna:nid"]`), then one cell block per element
 * family and cell type, in the order solids, thick shells, beams, shells.
 * Degenerate 8-node solids and thick shells collapse through
 * `detail::collapse_brick` (the keyword reader's rules), a shell with
 * `n4 == n3` is a triangle, ten-node solids (`NEL8 < 0`) are `tetra10` and
 * eight-node shells (`NEL48`) `quad8`. Element user ids are
 * `cell_data["lsdyna:eid"]` and user part ids `cell_data["lsdyna:part"]`; every
 * part with elements is a cell region tagged with its id and named by its
 * title, or `Part <id>` without one.
 *
 * Every state is a step: `time_step` selects one (0 = first, negative from the
 * end), its time is `field_data["meshio:time"]` with `lsdyna:state`. It
 * becomes
 *  - point data: `displacement` (current minus initial coordinates; the points
 *    stay where the geometry put them), `velocity`, `acceleration`,
 *    `temperature`, `heat_flux`, `mass_scaling`, `temperature_gradient`,
 *    `residual_forces`, `residual_moments`, whichever the header flags;
 *  - cell data: `stress`, `effective_plastic_strain`, `history_variables`,
 *    `strain`, `plastic_strain_tensor` and `thermal_strain_tensor` per
 *    integration point (solids: 1 or 8) or through-thickness layer (shells,
 *    thick shells: `MAXINT`; beams: their integration points) as
 *    `(cells, points * components)` with `field_data["lsdyna_d3plot:layout:<name>"]
 *    = [points, components]`, NaN where a family has fewer points or no such
 *    variable; `strain_inner`/`strain_outer`, `thickness`, `internal_energy`,
 *    `shell_bending_moment`, `shell_shear_force`, `shell_normal_force`,
 *    `shell_element_variables` for shells, `beam_axial_force`,
 *    `beam_shear_force`, `beam_bending_moment`, `beam_torsion_moment`,
 *    `beam_axial_stress`, `beam_shear_stress`, `beam_axial_strain` for beams,
 *    `thermal_variables` for solids;
 *  - the deletion flags as an int8 mask `lsdyna:alive` (1 = active, 0 =
 *    deleted), per cell (`MDLOPT = 2`) or per point (`MDLOPT = 1`); deleted
 *    cells stay in the mesh;
 *  - field data: `global_kinetic_energy`, `global_internal_energy`,
 *    `global_total_energy`, `global_velocity` and the per-part
 *    `part_internal_energy`, `part_kinetic_energy`, `part_velocity`,
 *    `part_mass`, `part_hourglass_energy`.
 *
 * Since v16.12.0 also: `d3part` files, 20/27-node hexahedra, SPH particles
 * (vertices with their variables), airbag particles and rigid road segments
 * (on points after the nodes), rigid-body motion (field data). Refused by
 * name: `intfor` files, femzip-compressed files, two-dimensional databases,
 * CFD/multi-solver data, adaptive remeshing and 21/15/20/40/64-node solids.
 *
 * See doc/formats/lsdyna_d3plot.md.
 */

// System includes
#include <cstddef>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Read one state of a d3plot family.
 * @param rPath the family's base file (`d3plot`)
 * @param rOpts `mTimeStep` selects the state; `mArrays`/`mPointsOnly` narrow
 *        the data
 * @return the mesh with the selected state's results
 * @throws ReadError if the file is not a d3plot, is truncated, uses a refused
 *         feature, is a numbered member of a family, or `mTimeStep` is out of
 *         range
 */
MESHIOPLUSPLUS_API Mesh read_lsdyna_d3plot(const std::string& rPath, const ReadOptions& rOpts = {});

/**
 * @brief The time of every state of the family, in family order.
 */
MESHIOPLUSPLUS_API std::vector<double> lsdyna_d3plot_time_values(const std::string& rPath);

/**
 * @brief Metadata (counts, data names, `mTimeValues`) of a d3plot family.
 */
MESHIOPLUSPLUS_API MeshMetadata read_lsdyna_d3plot_metadata(const std::string& rPath,
                                                            const ReadOptions& rOpts = {});

/**
 * @brief Whether a path names a d3plot base file: its file name is `d3plot`
 *        (any case). Format inference uses it, as `z88i1.txt` names Z88.
 */
MESHIOPLUSPLUS_API bool is_d3plot_filename(const std::string& rPath);

/**
 * @brief Whether the first bytes of a file are a d3plot control block.
 * @param pHead the file's first bytes
 * @param Size how many there are (a control block needs 256, or 512 in
 *        double precision)
 */
MESHIOPLUSPLUS_API bool is_d3plot_head(const char* pHead, std::size_t Size);

}  // namespace meshioplusplus
