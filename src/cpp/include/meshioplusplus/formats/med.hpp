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
 * @file med.hpp
 * @brief MED/Salome (.med) HDF5-backed C++ reader/writer.
 *
 * MED (Salome/Code-Aster) is the most structurally involved format
 * meshio++ supports. On disk it is HDF5 with groups written in order:
 * `INFOS_GENERALES` (MAJ/MIN/REL version triple), `ENS_MAA/<mesh_name>`
 * (mesh-level attrs `DIM`/`ESP`/`UNT`/`UNI`/`NOM`/`DES`) holding a single
 * time-step group with `NOE` (nodes: `COO` Fortran-order-flattened
 * coordinates, optional `FAM` per-point family id) and `MAI` (one group per
 * cell block keyed by MED type, e.g. `HE8`, each with a Fortran-order
 * 1-based `NOD` connectivity and optional `FAM`), plus `FAS/<mesh_name>`
 * (family/group definitions: `FAMILLE_ZERO`, `NOEUD/FAM_<id>_.../GRO/NOM`,
 * `ELEME/...` with the same layout) and `CHA/<field_name>` (fields).
 * Coordinate/connectivity arrays are stored **Fortran-ordered**
 * (column-major); this C++ implementation flattens/unflattens explicitly
 * since C++ has no native Fortran-order array type (`med.cpp`'s transpose
 * uses Eigen when `MESHIOPLUSPLUS_HAS_EIGEN`, else a hand-written
 * transpose).
 *
 * **What the C++ path handles** (matching the Python output byte-for-byte):
 * points, point/cell tags, families with `GRO` group names, **named regions**
 * derived from those families (one `Region` per group name — `Point` from
 * `NOEUD`, `Cell` from `ELEME`, `dim`/`tag` left unspecified since a name may
 * span several family ids; see `med_attach_point_regions`/
 * `med_attach_cell_regions` in `med.cpp` and doc/regions.md) and, on write,
 * synthesized back into families when the mesh carries no native
 * `point_tags`/`cell_tags` of its own (`med_point_regions_to_tags`/
 * `med_cell_regions_to_tags` — native data always wins, so a MED→MED round
 * trip is unaffected), an optional **`INFOS_GENERALES` version check**
 * (a file written by MED major version > 4 is rejected with a named error
 * instead of an obscure structural one further down), optional **`NUM`**
 * global point/cell numbering (`point_data`/`cell_data["med:num"]`; cell
 * `NUM` is only carried when *every* block has it), mesh-level metadata
 * (`mesh_name`/`description`/`unit_time`/`unit_coords`/
 * `point_tag_groups`/`cell_tag_groups`, all carried via #MedInfo), the
 * fixed node-orientation permutations for linear 3D types (`tetra`,
 * `pyramid`, `wedge`, `hexahedron` — see `_med_node_perm` in
 * doc/formats/med.md), and ragged `POG`/`POG2` polygon cell blocks (CSR
 * `NOD`+`INN` offset arrays, copied — not zero-copy — across the C++/Python
 * boundary since ragged data cannot be viewed in place). The `MAI` cell
 * blocks are iterated in HDF5 **creation order** (matching h5py's
 * `track_order`) since block order must align with `cell_data`/`cell_sets`.
 *
 * **What always falls back to Python** (the C++ functions `throw` and the
 * `meshioplusplus.med` shim catches and retries with the pure-Python/h5py
 * implementation): a `CHA` **field** past the single-timestep, no-profile,
 * no-units common case (MED-4.1 bitmask attributes,
 * `field_data["med:field_units"]`/`["med:step_meta"]`, and multi-timestep
 * field-name grouping are Python-only — see `read_cha_fields`/
 * `write_cha_nodal_field`/`write_cha_cell_field`), the `gmsh:physical`→family
 * **bridging** performed on write, non-default **profiles** / `ELGA`
 * support, and **multi-mesh** files (`read_med_multi`/`write_med_multi`,
 * which have no C++ equivalent at all). Quadratic 3D types (`tetra10`,
 * `hexahedron20`, `pyramid13`, `wedge15`) share the linear types' orientation
 * convention but have no implemented corners+midpoints permutation yet —
 * they round-trip unconverted (a warning is logged the first time one is
 * seen); see doc/formats/med.md for the planned fix.
 */

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Side-channel struct carrying MED mesh-level metadata that the
 *        zero-copy Mesh conversion layer cannot carry (custom attributes
 *        on the Python `meshio.Mesh`, not `point_data`/`cell_data`/
 *        `field_data`). The binding layer `setattr`s these fields onto the
 *        Python `Mesh` on read, and reads them back off it on write.
 */
struct MedInfo {
    /**
     * Per-point family/tag membership: `set_id -> [subset_name, ...]`,
     * corresponding to Python `mesh.point_tags`. Populated on read from
     * `FAS/NOEUD` family group names; consumed on write to build the
     * `NOEUD` family definitions (and each point's `FAM` id, which is
     * itself carried in `point_data["point_tags"]`, not here).
     */
    std::map<std::int64_t, std::vector<std::string>> mPointTags;
    /**
     * Per-cell-block family/tag membership: `set_id -> [subset_name, ...]`,
     * corresponding to Python `mesh.cell_tags`. Populated on read from
     * `FAS/ELEME` family group names (including the synthetic families
     * created for OpenFOAM boundary patches or Gmsh physical groups when
     * bridged elsewhere); consumed on write analogously to `point_tags`.
     */
    std::map<std::int64_t, std::vector<std::string>> mCellTags;
    /** `field_data["med:nom"]` — one component-name list per field, in
     * field-iteration order (point_data fields first, then cell_data
     * fields); each field's `NOM` attribute is the 16-char-padded
     * concatenation of its entry here. */
    std::vector<std::vector<std::string>> mMedNom;  // field_data["med:nom"]

    // Mesh-level metadata attributes (custom attributes on the Python Mesh).
    /** Python `mesh.mesh_name` — the `ENS_MAA` group name and the mesh's
     * `NOM` attribute value; defaults to `"mesh"` when absent. */
    std::string mMeshName = "mesh";
    /** Python `mesh.description` — the `DES` attribute; defaults to
     * `"Mesh created with meshio++"` on write when unset. */
    std::string mDescription;
    /** Python `mesh.unit_time` — the `UNT` attribute (physical unit of the
     * time axis, e.g. `"s"`). */
    std::string mUnitTime;
    /** Python `mesh.unit_coords` — the `UNI` attribute (physical unit of
     * the coordinate axes, e.g. `"m"`). Values round-trip through
     * `latin-1` and are stripped of surrounding whitespace/NUL padding on
     * read. */
    std::string mUnitCoords;
    // set_id -> family link name (e.g. "FAM_2_Side").
    /** `set_id -> "FAM_<id>..."` short family link name, mirroring Python
     * `mesh.point_tag_groups`; always present (possibly empty) after any
     * Python `read()` regardless of whether the source file had a `FAS`
     * section. */
    std::map<std::int64_t, std::string> mPointTagGroups;
    /** `set_id -> "FAM_<id>..."` short family link name for cell-block
     * families, mirroring Python `mesh.cell_tag_groups`; same defaulting
     * behavior as `point_tag_groups`. */
    std::map<std::int64_t, std::string> mCellTagGroups;

    // --- populated only under ReadOptions::mLenient (see read_med) ----------
    /**
     * @brief What `ReadOptions::mLenient` skipped, in the order it was met.
     *
     * Mirrors `MdpaInfo::mSkippedConstructs`. Empty after a strict read, since
     * a strict read either represents everything or throws. Each entry names
     * the field and the construct, e.g.
     * `"field 'v' on a named profile"`.
     */
    std::vector<std::string> mSkippedConstructs;
    /**
     * @brief `field name -> (UNI, UNT)` — a field's physical and time unit.
     *
     * The C++ `Mesh` has nowhere to put these (they are strings, and
     * `NDArray` has no string dtype; Python carries them as the dict-valued
     * `field_data["med:field_units"]`, which cannot cross any binding). Under
     * `mLenient` they are read here rather than thrown on, so a caller with no
     * Python fallback gets the field values *and* can still see its units.
     */
    std::map<std::string, std::pair<std::string, std::string>> mFieldUnits;
    /**
     * @brief `field name -> (NDT, NOR, PDT)` of the step actually read.
     *
     * Recorded whenever it is not the write-side default `(1, -1, 0.0)`, i.e.
     * exactly when a strict read would have declined. Python's counterpart is
     * the dict-valued `field_data["med:step_meta"]`.
     */
    std::map<std::string, std::tuple<std::int64_t, std::int64_t, double>> mStepMeta;
    /**
     * @brief `field name -> the PDT of every timestep the field carries`.
     *
     * Always filled (one entry per `CHA` field), so a caller can discover
     * what `ReadOptions::mTimeStep` may select before issuing the request —
     * the same reason `MeshMetadata::mTimeValues` exists for exodus.
     */
    std::map<std::string, std::vector<double>> mFieldTimeValues;
};

/**
 * @brief Read a MED (.med) HDF5 file into a Mesh, handling the
 *        mesh-representation subset described in the file-level docs.
 *
 * Reads points (un-transposing the Fortran-ordered `COO` dataset), the
 * `MAI` cell blocks in HDF5 creation order, per-point/per-cell `FAM` tag
 * arrays (exposed as `point_data["point_tags"]`/`cell_data["cell_tags"]`),
 * optional `NUM` global numbering (`point_data`/`cell_data["med:num"]`, cell
 * `NUM` only when every block has it), family/group names from `FAS`
 * (searched first under the mesh's own time-step group, then at the top
 * level) — attached additionally as named `Point`/`Cell` regions, one per
 * group name (see `med_attach_point_regions`/`med_attach_cell_regions`) —
 * mesh-level metadata, the fixed node-orientation permutation for linear 3D
 * types, and ragged `POG`/`POG2` polygon blocks (materialized as a copied
 * `list`-like ragged `CellBlock` since they cannot be represented as a
 * rectangular NDArray without loss).
 *
 * @param rPath filesystem path to the .med file to read
 * @param rInfo output side-channel struct populated with tags, families,
 *        and mesh-level metadata (see #MedInfo)
 * @return the read Mesh (points, cells, point_data["point_tags"],
 *         cell_data["cell_tags"], named regions, arbitrary named point/cell
 *         data from `CHA` fields except those excluded below)
 * @throws ReadError — on a file written by MED major version > 4; on a `CHA`
 *         field past the single-timestep/no-profile/no-units common case
 *         (units, multi-timestep metadata, a named profile, or ELNO/ELGA
 *         support); on multi-mesh files; on malformed/unsupported HDF5
 *         layout. Callers (the Python shim) catch this and retry with the
 *         pure-Python/h5py reader.
 */
MESHIOPLUSPLUS_API Mesh read_med(const std::string& rPath, MedInfo& rInfo);

/**
 * @brief `read_med` with read options — the overload that makes a MED file
 *        with enhanced `CHA` fields readable where there is no Python.
 *
 * Two options change anything here; the narrowing ones (`mPointsOnly`,
 * `mDataArrays`, `mMmap`) are applied by the Python layer after any read, as
 * for every other format.
 *
 * **`mLenient`** downgrades every `CHA` decline listed under the plain
 * overload to a `log::warn` plus a recorded entry in
 * `MedInfo::mSkippedConstructs`, so the mesh, its tags/families/regions and
 * every *representable* field still come back. Units and non-default step
 * metadata are not skipped but **read into** `MedInfo::mFieldUnits` /
 * `mStepMeta`; only a construct with no representation at all (a named
 * profile, an ELNO/ELGA support, a field mixing nodal and cell support) causes
 * that one field to be dropped. This is the same mechanism, and the same
 * rationale, as `ReadOptions::mLenient` for MDPA: strict is what the Python
 * shim uses (so it still falls back and the Python surface is unchanged),
 * lenient is the only way a real Salome/Code_Aster file reads at all from the
 * C API, Fortran, Julia, R, WASM or the native CLI.
 *
 * **`mTimeStep`** selects which timestep of a multi-step field to read
 * (0 = first, negative counts from the end), resolved per field via
 * `ResolveTimeStep`. A *non-default* `mTimeStep` is honoured whether or not
 * `mLenient` is set — it is an explicit request, and no Python behaviour
 * depends on it. With the default step and no `mLenient`, a multi-step field
 * still throws, exactly as before. Every field's available step times are
 * reported in `MedInfo::mFieldTimeValues` regardless.
 *
 * @param rPath filesystem path to the .med file to read
 * @param rInfo output side-channel struct (see #MedInfo)
 * @param rOptions read options; only `mLenient` and `mTimeStep` are consulted
 * @return the read Mesh
 * @throws ReadError as the plain overload, minus whatever `mLenient` and
 *         `mTimeStep` cover; a `mTimeStep` out of range throws naming the
 *         field and its step count rather than clamping.
 */
MESHIOPLUSPLUS_API Mesh read_med(const std::string& rPath, MedInfo& rInfo,
                                 const ReadOptions& rOptions);

/**
 * @brief Write a Mesh to a MED (.med) HDF5 file, handling the
 *        mesh-representation subset described in the file-level docs.
 *
 * Writes `INFOS_GENERALES` (parsed from `med_version`, falling back to
 * `4, 1, 0` if unparsable), `ENS_MAA/<rInfo.mMeshName>` with points
 * (Fortran-order-flattened) and one `MAI/<MED type>` group per cell block
 * (rejecting up front with `WriteError` if two blocks share a MED type,
 * since MED cannot represent that), `FAS` family definitions built from
 * `rInfo.mPointTags`/`rInfo.mCellTags` — or, when the mesh carries no native
 * `point_tags`/`cell_tags` of its own, synthesized from any `Point`/`Cell`
 * regions the mesh carries (`med_point_regions_to_tags`/
 * `med_cell_regions_to_tags`; native data always wins, so a MED→MED round
 * trip through this writer is unaffected; `Side` regions have no MED
 * equivalent and are dropped with a warning) — a family with no groups omits
 * `GRO` entirely, and the fixed node-orientation permutation applied to
 * linear 3D cell types before writing `NOD`. Optional `NUM` global numbering
 * is written when `point_data`/`cell_data["med:num"]` is present. Ragged
 * `polygon`/`polygon2` blocks are written as `POG`/`POG2` CSR data. Family
 * names longer than 80 bytes after `latin-1` encoding raise `WriteError`
 * rather than truncating.
 *
 * @param rPath filesystem path to the .med file to create/overwrite
 * @param rMesh the mesh to write
 * @param rInfo side-channel struct supplying tags, families, and mesh-level
 *        metadata (see #MedInfo); read from the Python Mesh's custom
 *        attributes by the binding layer
 * @param rMedVersion the `MAJ.MIN.REL` triple written to
 *        `INFOS_GENERALES` (default `"4.1.0"`)
 * @throws WriteError — if the mesh carries `CHA`-worthy fields (any
 *         point_data/cell_data beyond `point_tags`/`cell_tags`/`med:num`
 *         that this path doesn't handle), `gmsh:physical` bridging is
 *         needed, two cell blocks share one MED type, or a family name
 *         exceeds 80 bytes. Callers (the Python shim) catch this and retry
 *         with the pure-Python/h5py writer.
 * @note point_data/cell_data keys produced/consumed: `"point_tags"`,
 *       `"cell_tags"`, `"med:num"`.
 */
MESHIOPLUSPLUS_API void write_med(const std::string& rPath, const Mesh& rMesh, const MedInfo& rInfo,
                                  const std::string& rMedVersion = "4.1.0");

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
