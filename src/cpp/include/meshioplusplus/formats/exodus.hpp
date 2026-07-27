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
 * @file exodus.hpp
 * @brief Exodus II (.e/.exo/.ex2) C++ reader/writer, stored in netCDF using
 *        its classic variable/dimension conventions.
 *
 * Key variables: `coord(num_dim, num_nodes)` (transposed relative to
 * meshio++'s `(n, dim)` layout) or separate `coordx`/`coordy`/`coordz`
 * (both accepted on read); `eb_prop1(num_el_blk)` arbitrary distinct block
 * ids; `connect{k}(num_el_in_blk{k}, num_nod_per_el{k})` per element block
 * with a text `elem_type` attribute and 1-based node indices;
 * `name_nod_var`/`vals_nod_var{k}` point-data (first timestep only);
 * `name_elem_var`/`vals_elem_var{idx}[eb{block}]` cell data, concatenated
 * across blocks then re-split by target cell-block size. Compiled in only
 * when `MESHIOPLUSPLUS_HAS_NETCDF` is defined; otherwise the Python
 * `netCDF4` fallback handles this format. See doc/formats/exodus.md.
 */

#ifdef MESHIOPLUSPLUS_HAS_NETCDF

// System includes
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Exodus provenance strings, which the conversion layer cannot carry.
 *
 * `qa_records` (who wrote the file, with what version, when) and `info_records`
 * (free-form notes) are text, and `NDArray` has no string dtype -- so unlike
 * every other thing a reader produces, these cannot ride on the `Mesh`. They
 * travel in this side-channel struct instead, the same shape `MedInfo` and
 * `OpenFoamInfo` already use: the pybind binding `setattr`s the result onto the
 * Python `Mesh` as `info`, and the flat bindings (C, Fortran, WASM) construct a
 * local and drop it, exactly as `registry.cpp` already does for `MedInfo`.
 *
 * This is why reading no longer throws on these variables. Every file SEACAS,
 * Cubit or Sierra writes carries `qa_records`, so treating them as an
 * unsupported construct made the format unreadable wherever there is no Python
 * fallback to defer to -- which is to say, all of WASM.
 */
struct ExodusInfo {
    /**
     * @brief `info_records` followed by `qa_records`, flattened.
     *
     * The order reproduces the Python reference reader's exactly (it appends to
     * one `info` list in netCDF variable order, four strings per QA record), so
     * `mesh.info` is identical whichever path produced it.
     */
    std::vector<std::string> mInfoRecords;
};

/**
 * @brief `cell_data` name prefix that marks an Exodus per-element attribute.
 *
 * Exodus stores a fixed number of floating-point *attributes* per element of a
 * block (`attrib{k}`, named by `attrib_name{k}`) alongside its connectivity --
 * the standard home for a `SPHERE`/`CIRCLE` element's radius, a beam's cross
 * section, a shell's thickness. They are per-cell values, so meshio++ carries
 * them as ordinary `cell_data`, under this prefix.
 *
 * The prefix is what makes the round trip unambiguous in both directions: on
 * read it keeps an attribute from colliding with a same-named element *variable*
 * (`name_elem_var`, a genuinely different concept -- attributes are constant in
 * time, element variables are per-time-step), and on write it is the only signal
 * telling the writer which `cell_data` arrays belong in `attrib{k}`.
 *
 * So a `RADIUS` attribute reads back as `cell_data["exodus:attr:RADIUS"]`.
 */
inline constexpr const char* kExodusAttributePrefix = "exodus:attr:";

/**
 * @brief Write `mesh` as an Exodus II (netCDF classic) file.
 *
 * Writes global attrs (`title`, `version=5.1f`, `api_version=5.1f`,
 * `floating_point_word_size=8`), a dummy single `0.0` `time_whole` step, one
 * `connect{k}` variable per cell block (element type mapped through the
 * canonical meshio++ -> Exodus reverse table, e.g. `hexahedron -> HEX8`,
 * `tetra -> TETRA`, `tetra4 -> TET4` as a distinct entry from plain
 * `tetra`), and point_data/cell_data as `vals_nod_var`/`vals_elem_var`
 * variables.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError if a cell block's type has no entry in the meshio++ ->
 *         Exodus type table, or if the connectivity dtype is unsupported
 * @note the shim only attempts this C++ path when `mesh.point_sets` is
 *       empty — the C++ writer has no support for Exodus node sets at all
 */
MESHIOPLUSPLUS_API void write_exodus(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read an Exodus II (netCDF classic) file.
 *
 * Reads coordinates (either `coord` or `coordx`/`coordy`/`coordz`), one cell
 * block per `connect{k}` variable (via its `elem_type` attribute and the
 * Exodus -> meshio++ type table), and point_data/cell_data — with the
 * point-data name recombination quirk `categorize()` reproduces on purpose
 * from the reference implementation: names ending `X`/`Y`/`Z` (or `_R`/`_Z`)
 * are stacked into a 3- (or 2-) component vector when a sibling exists, but
 * the "sibling found" check uses Python truthiness on the found variable
 * index, so index `0` is treated the same as "not found" — a latent
 * reference-implementation edge case deliberately preserved rather than
 * fixed, so the two implementations agree.
 *
 * ## Named regions
 *
 * Element blocks, node sets and side sets all become #Region s -- the three
 * things Exodus spells separately that meshio++ spells one way:
 *  - `connect{k}` -> `RegionKind::Cell`, named from `eb_names` (falling back to
 *    `"Block <id>"`), tagged with its `eb_prop1` id. Two blocks of the *same*
 *    element type stay distinguishable, which is the whole point of the id.
 *  - `node_ns{k}` -> `RegionKind::Point`, named from `ns_names`, tagged from
 *    `ns_prop1`.
 *  - `elem_ss{k}`/`side_ss{k}` -> `RegionKind::Side`, named from `ss_names`, as
 *    `(global cell, local facet)` pairs. Exodus numbers an element's sides its
 *    own way, so the facet column is remapped through `exo_face_index` rather
 *    than stored raw -- see that function for the per-type tables.
 *
 * ## Time steps
 *
 * `ReadOptions::mTimeStep` selects which step of `time_whole` the data arrays
 * come from (0 = first, the historical behaviour; negative counts from the
 * end). An out-of-range request throws naming the available count.
 *
 * @param rPath filesystem path to read
 * @param rInfo receives the `qa_records`/`info_records` provenance strings
 * @param rOptions per-call reader options; defaulted, so the historical
 *        behaviour is what a plain `read_exodus(path)` still does
 * @return the read Mesh, with regions as above
 * @throws ReadError if a variable has an unsupported netCDF type, point-data
 *         names are inconsistent, a `connect{k}` names an unknown Exodus
 *         element type, the connectivity dtype is unsupported, or the
 *         requested time step is out of range
 * @note point_data keys ending X/Y/Z or _R/_Z may be recombined into vector
 *       arrays; cell_data is split per cell block by node count
 */
MESHIOPLUSPLUS_API Mesh read_exodus(const std::string& rPath, ExodusInfo& rInfo, const ReadOptions& rOptions = {});

/**
 * @brief Read an Exodus II file, discarding the provenance side channel.
 *
 * The `ReadFn`-shaped overload the registry and the flat bindings use.
 */
MESHIOPLUSPLUS_API Mesh read_exodus(const std::string& rPath, const ReadOptions& rOptions = {});

/**
 * @brief Summarize an Exodus II file without materializing its data arrays.
 *
 * Fills `MeshMetadata::mTimeValues` from `time_whole`, which is what makes the
 * available step count discoverable before committing to a read.
 *
 * @param rPath filesystem path to read
 * @param rOptions per-call reader options
 * @return the summary
 */
MESHIOPLUSPLUS_API MeshMetadata read_exodus_metadata(const std::string& rPath, const ReadOptions& rOptions = {});

/**
 * @brief Map an Exodus 1-based side number to a meshio++ local facet index.
 *
 * Exodus numbers an element's sides in its own order, which is *not*
 * `detail::cell_faces`' order -- storing the raw number would silently point a
 * `Side` region at the wrong face. The per-type tables live in the definition,
 * with the node lists they were derived from, mirroring `abq_face_index` in
 * `abaqus.cpp`; a gtest pins each entry against `cell_faces` so a transcription
 * slip cannot survive.
 *
 * @param rCellType meshio++ cell type name of the owning element
 * @param ExodusSide 1-based Exodus side number
 * @return the local facet index, or -1 when the pair has no mapping
 */
MESHIOPLUSPLUS_API int exo_face_index(const std::string& rCellType, int ExodusSide);

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
