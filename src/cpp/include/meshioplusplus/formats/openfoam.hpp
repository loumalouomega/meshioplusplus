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
 * @file openfoam.hpp
 * @brief OpenFOAM polyMesh (read-only) C++ reader.
 *
 * A polyMesh is a directory of sibling "FoamFile"-headered files (`points`,
 * `faces`, `owner`, `neighbour`, `boundary`), each ASCII or binary
 * (little-endian only; `label=32/64`, `scalar=32/64` per the file's `arch`
 * header string). `points` (`vectorField`) and `owner`/`neighbour`
 * (`labelList`) are flat contiguous buffers read directly; `faces`
 * (`faceList`) is non-contiguous (each face is its own length-prefixed
 * `labelList`) and is read via a two-pass CSR gather bounded in peak
 * memory. `boundary` is a `patch_name -> {type, nFaces, startFace}` table
 * parsed with a brace-matching regex.
 *
 * Cells are reconstructed from the owner/neighbour/face topology: each
 * face is oriented outward from its owning cell (reversed if the cell is
 * that face's neighbour), then classified by `(n_faces, n_points)` into
 * `tetra` (4,4), `pyramid` (5,5), `wedge` (5,6), `hexahedron` (6,8) — each
 * with a dedicated orientation-fixing builder that flips node order if a
 * scalar triple product comes out negative — or, for any other signature,
 * a general `polyhedron<N>` **ragged** cell block (ragged data crosses the
 * C++/Python boundary as a copied list of face-node arrays, never
 * zero-copy). Boundary faces become `triangle`/`quad`/`polygon<N>` blocks,
 * one per patch/size combination.
 *
 * `write_openfoam` (v9.20.0) is the inverse, and is the **only writer in
 * meshio++ that takes a directory path** — it creates
 * `<case>/constant/polyMesh/` and writes all five files. It goes through
 * `detail/face_mesh.hpp`'s global face table, which is also what CGNS's
 * `NFACE_n` writer uses. ASCII only; a binary polyMesh is a follow-up.
 *
 * Only mesh topology is read or written; OpenFOAM field files (`U`, `p`,
 * `T`, …) under a case's time directories are never touched by this module,
 * so no `point_data`/`field_data` is ever produced or consumed.
 */

// System includes
#include <cstdint>
#include <map>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Side-channel struct carrying OpenFOAM boundary-patch tag data
 *        that the zero-copy Mesh conversion layer cannot carry (`Python
 *        mesh.cell_tags` is a custom Mesh attribute, not `cell_data`). The
 *        binding layer `setattr`s this onto the returned Python `Mesh`.
 */
struct OpenFoamInfo {
    // MED-style negative family id -> {patch name}.
    /**
     * `family_id -> [patch_name]`, mirroring Python `mesh.cell_tags`. Each
     * boundary patch gets a distinct negative "MED-style family id"
     * `-(patch_index+1)` (assigned once per patch and reused across
     * whichever face-size cell blocks that patch's faces fall into); the
     * matching `cell_data["cell_tags"]` array on the returned Mesh holds
     * `0` for every volume-cell block and the patch's family id for its
     * boundary-face blocks. This lets a subsequent MED write bridge patch
     * names through the same family mechanism used for Gmsh physical
     * groups (see doc/formats/med.md).
     */
    std::map<std::int64_t, std::vector<std::string>> mCellTags;

    /**
     * `family_id -> patch type` — the `type` entry of the on-disk `boundary`
     * file (`patch`, `wall`, `symmetry`, `symmetryPlane`, `empty`, `wedge`,
     * `cyclic`, …), keyed by the **same** negative family id as #mCellTags
     * rather than by patch name: that key is this format's primary key
     * everywhere else, and two maps keyed differently invite a bad join.
     *
     * A patch whose `type` the file omitted simply has no entry, and the
     * writer then emits `patch` — OpenFOAM's base type, always safe.
     * Deliberately **not** `wall`: `wall` selects wall functions and
     * `nut*WallFunction` boundary behaviour, so guessing it would silently
     * change a solve's physics.
     *
     * Types needing companion keys this struct cannot carry (`cyclic`'s
     * `neighbourPatch`, `processor`'s `myProcNo`, `mappedWall`'s `sample`,
     * …) are **downgraded to `patch` with a warning** on write: emitting
     * them bare produces a case OpenFOAM refuses to load, whereas a
     * downgraded case loads and solves with boundary conditions the user
     * can see and fix.
     */
    std::map<std::int64_t, std::string> mPatchTypes;
};

// `path` may be a `.foam` marker file, a case directory, or a polyMesh
// directory (resolved like the Python reader).
/**
 * @brief Read an OpenFOAM polyMesh into a Mesh.
 *
 * `path` may be a `.foam` marker file (looks for
 * `<parent>/constant/polyMesh`), a directory literally named `polyMesh`
 * (used as-is), or any other directory (checked for `constant/polyMesh`
 * then `polyMesh` as subdirectories) — resolved identically to the Python
 * reader's `_resolve_polymesh`. Reconstructs volume cells
 * (tetra/pyramid/wedge/hexahedron/general polyhedron) and boundary faces
 * (triangle/quad/polygon) from the `points`/`faces`/`owner`/`neighbour`/
 * `boundary` files, auto-detecting ASCII vs binary and label/scalar width
 * per file. Degenerate volume cells that match a named type's
 * `(n_faces, n_points)` signature but whose topology doesn't cleanly
 * resolve are silently skipped (logged as a warning count) rather than
 * demoted to a general polyhedron.
 *
 * @param rPath a `.foam` file, case directory, or polyMesh directory
 * @param rInfo output side-channel struct populated with boundary-patch
 *        family ids and names (see #OpenFoamInfo)
 * @return the read Mesh: points, volume + boundary cell blocks,
 *         `cell_data["cell_tags"]` (0 for volume blocks, a per-patch
 *         negative id for boundary blocks), `mesh.point_tags` always set
 *         to `{}` (OpenFOAM has no point-tag concept; present only for
 *         interface symmetry with the MED-derived tag convention) — no
 *         point_data or field_data
 * @throws ReadError / std::filesystem-related errors if no polyMesh
 *         directory can be resolved, or on a malformed/unsupported file;
 *         callers (the Python shim) catch this and retry with the
 *         pure-Python reader
 */
MESHIOPLUSPLUS_API Mesh read_openfoam(const std::string& rPath, OpenFoamInfo& rInfo);

/**
 * @brief Write a Mesh as an OpenFOAM polyMesh case.
 *
 * The **only meshio++ writer that creates a directory**: @p rPath is
 * resolved exactly as `read_openfoam` resolves it (a `.foam` marker file, a
 * directory literally named `polyMesh`, or any other directory taken as the
 * case root), and `<case>/constant/polyMesh/` is created if absent. A
 * `.foam` target additionally gets its empty marker file written, which is
 * what makes the case openable by ParaView and by this reader.
 *
 * Volume cells become the `faces`/`owner`/`neighbour` triple via
 * `detail::build_global_faces`, which repairs each cell's winding — so an
 * inverted cell is written correctly oriented rather than rejected by
 * `checkMesh`. All four ordering rules OpenFOAM requires (internal faces
 * first, `owner < neighbour`, faces sorted by owner then neighbour, normals
 * pointing owner→neighbour) are enforced and then re-validated before
 * anything is written; a violation is a `WriteError` naming the rule,
 * because it means an internal bug rather than bad input.
 *
 * Boundary patches are recovered from `cell_data["cell_tags"]`'s negative
 * values together with @p rInfo. A mesh carrying none — anything converted
 * from another format — gets a single `defaultFaces` patch of type `patch`,
 * which is what `blockMesh` itself produces and yields a loadable case.
 * Patches are **not** synthesized from geometry.
 *
 * ASCII only: a binary polyMesh is a documented follow-up, so an explicit
 * binary request fails by name rather than silently writing ASCII.
 *
 * @param rPath a `.foam` file, case directory, or polyMesh directory
 * @param rMesh the mesh to write; ragged polyhedron blocks are supported
 *        and are in fact this format's native cell shape
 * @param rInfo patch names and types, keyed by the same negative family ids
 *        `read_openfoam` produces (see #OpenFoamInfo). An empty struct is
 *        valid and yields the `defaultFaces` case above.
 * @throws WriteError if the mesh has no volume cells, if a 3D block cannot
 *         contribute faces (which would silently drop solids), if a face is
 *         shared by three or more cells, or if the output directory cannot
 *         be created
 */
MESHIOPLUSPLUS_API void write_openfoam(const std::string& rPath, const Mesh& rMesh,
                                       const OpenFoamInfo& rInfo);

}  // namespace meshioplusplus
