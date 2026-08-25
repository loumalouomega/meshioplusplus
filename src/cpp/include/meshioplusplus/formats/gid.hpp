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
 * @file gid.hpp
 * @brief GiD postprocess format reader and writer.
 *
 * GiD (https://www.gidsimulation.com/) is a pre/postprocessor widely used in
 * the same FE community meshio++ serves. Its postprocess format is **written**
 * through CIMNE's gidpost C library, vendored as a hardcopy at
 * `src/cpp/third_party/gidpost/` (BSD-2-Clause-Views) — see that directory's
 * `README.meshioplusplus.md` for what was vendored and why.
 *
 * gidpost's public API has **zero read functions**, so the **reader** is
 * meshio++'s own code against the on-disk grammar, entirely independent of the
 * vendored library. That independence has a visible consequence worth stating
 * up front: the reader's real dependencies are **per flavour** — ASCII needs
 * nothing at all, the compressed-binary and gzipped-ASCII flavours need zlib,
 * and the HDF5 flavour needs HDF5 — whereas *writing* any flavour needs
 * gidpost, which is itself hard-gated on zlib. So `gid` is **readable in
 * strictly more build configurations than it is writable**: the
 * statically-linked release CLI binaries and the Windows wheels build with
 * zlib off and cannot write GiD at all, yet still read the ASCII flavour.
 * `gid_available` reports the write side, `gid_readable` the read side; they
 * genuinely differ, which is why there are two.
 *
 * Three on-disk flavours, chosen by `GidMode` or inferred from the path:
 *  - **Ascii**: two sibling files, `<stem>.post.msh` (geometry) +
 *    `<stem>.post.res` (results), human-readable.
 *  - **Binary**: one deflated file, `<stem>.post.bin`.
 *  - **Hdf5**: one HDF5 file, `<stem>.post.h5` (needs
 *    `MESHIOPLUSPLUS_WITH_HDF5=ON` in addition to gidpost itself).
 *
 * Hard-gated on zlib: gidpost's own `gidpostInt.h` includes `<zlib.h>`
 * unconditionally and the binary flavour is always deflated, so zlib is a
 * build prerequisite for the whole library, not for one of its three modes.
 * `MESHIOPLUSPLUS_WITH_GIDPOST` (default ON) auto-disables — with an
 * actionable `STATUS` message, never a hard configure error — whenever zlib
 * is unavailable; compiled out, `write_gid` still exists and throws naming
 * both CMake flags (the `partition_kahip_parts` contract), so `.post.msh`
 * can never silently resolve to another format.
 *
 * ## Mesh mapping
 *
 * GiD has exactly ten element types (`GiD_ElementType`); higher-order
 * variants share a type with a larger node count. One meshio cell block
 * becomes one named GiD mesh (`"<celltype>_<blockindex>"`), since GiD
 * results reference meshes by name. Only the four linear types plus their
 * one-level-quadratic siblings whose node ordering has been independently
 * verified against GiD's own geometry (pinned by `tests/cpp/test_gid.cpp`'s
 * `GidOrdering` suite, never by a round trip through a reader that does not
 * exist) are supported: `vertex`, `line`/`line3`, `triangle`/`triangle6`,
 * `quad`/`quad8`/`quad9`, `tetra`/`tetra10`, `hexahedron`/`hexahedron20`,
 * `wedge`, `pyramid`. Everything else — `hexahedron27`, `wedge15`,
 * `pyramid13` (orderings not yet verified), `polygon`/`polyhedron` (GiD has
 * no such type), every `VTK_LAGRANGE_*` and higher-degree Lagrange type —
 * throws a `WriteError` naming the offending type, never a silent drop or a
 * guessed permutation.
 *
 * Points are written once (the first mesh's coordinate block; every
 * subsequent mesh gets an empty `Coordinates`/`End Coordinates` pair, which
 * gidpost's own state machine requires), always as 3-D (z padded with 0 for
 * a 2-D mesh, since gidpost's coordinate writer always emits three columns
 * regardless of the declared dimension). Element ids are globally unique
 * 1-based integers across every block (gidpost numbers `1..n` *per mesh*
 * otherwise, which collides the moment a mesh has more than one block).
 *
 * An integral `cell_data` array named `"gmsh:physical"` is written as each
 * element's material id (`GiD_fWriteElementsIdMatBlock`) and is excluded
 * from the result output, since it is already in the geometry file. No
 * other key is consulted for material ids in v1.
 *
 * ## Data mapping
 *
 * `point_data` is written `GiD_OnNodes`. `GiD_ResultLocation` has no
 * "on cells" concept, so `cell_data` is written `GiD_OnGaussPoints` against
 * a synthetic one-point Gauss-point set per block (`"gp_<mesh_name>"`) — the
 * standard GiD idiom for a per-element field, and how Kratos's own GiD
 * writer represents one. An array spanning several blocks becomes several
 * result blocks sharing one result name but different Gauss-point sets.
 *
 * `GiD_ResultType`'s valid component counts are irregular (Scalar 1; Vector
 * 2/3/4; Matrix 3/6; MainMatrix 12; ...), and gidpost does not validate an
 * invalid count itself — it silently emits a malformed file. This writer
 * validates instead: 1 component → `GiD_Scalar`; 2 or 3 → `GiD_Vector`;
 * anything else splits into that many named `GiD_Scalar` results
 * (`"<name>_1"` .. `"<name>_k"`), recorded with a `provenance_note`. A
 * 6-component array is deliberately **not** mapped to `GiD_Matrix` even
 * though stress tensors are GiD's canonical use case: meshio's `(n,6)`
 * carries no declaration that it *is* a symmetric tensor, and GiD's own
 * component order (`xx,yy,xy,zz,xz,yz`) differs from meshio/VTK's
 * (`xx,yy,zz,xy,yz,xz`) — mapping on shape alone would silently permute six
 * possibly-unrelated scalars. Splitting is lossless and unambiguous; a
 * tensor-aware side channel is a documented follow-up.
 *
 * Named regions, `field_data`, and multi-step results are not carried
 * (each dropped with a `provenance_note` or `log::warn`); this writes
 * exactly one step (`rStepValue`).
 *
 * ## Provenance
 *
 * `SlotTier::Block`: `provenance_lines(SlotTier::Block)` is rendered one
 * line per `GiD_fWriteMeshUserAttribute`/`GiD_fWriteResultUserAttribute`
 * call (`"meshio++"` for line 0, `"meshio++_<n>"` for any further line),
 * which gidpost itself renders as the `# Name: value` comment its own
 * source documents. There is **no pure-Python reference writer** for this
 * format — gidpost cannot be cheaply reimplemented in pure Python, and a
 * second implementation of the discrete node-ordering permutations risks
 * disagreeing near the exact cases the bytes tests exist to pin — so
 * `meshioplusplus.gid` is a C++-core-only surface (the `openfoam` writer's
 * precedent, generalized: here there is no fallback engine at all, not
 * merely one that is unconditionally used).
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/// Which of gidpost's three (of four) post-file flavours to emit.
enum class GidMode : int {
    /// Infer from `rPath`'s extension: `.post.bin` -> Binary, `.post.h5` ->
    /// Hdf5, anything else (including `.post.msh`/`.post.res`) -> Ascii.
    /// The default, and what every binding but a caller passing `mode=`
    /// explicitly uses.
    Auto = 0,
    /// Two sibling files, `<stem>.post.msh` + `<stem>.post.res`.
    Ascii = 1,
    /// One deflated file, `<stem>.post.bin`.
    Binary = 2,
    /// One HDF5 file, `<stem>.post.h5`. Needs `MESHIOPLUSPLUS_WITH_HDF5=ON`
    /// in addition to `MESHIOPLUSPLUS_WITH_GIDPOST=ON`.
    Hdf5 = 3,
};

/**
 * @brief Parses a `GidMode` from its lower-case name (`"auto"`, `"ascii"`,
 * `"binary"`, `"hdf5"`), the flat-binding/CLI/Python spelling.
 * @throws std::invalid_argument on an unrecognized name.
 */
MESHIOPLUSPLUS_API GidMode gid_mode_from_name(const std::string& rName);

/**
 * @brief Whether this build can emit @p mode.
 *
 * `Ascii`/`Binary` follow `MESHIOPLUSPLUS_HAS_GIDPOST`; `Hdf5` additionally
 * needs `MESHIOPLUSPLUS_HAS_GIDPOST_HDF5` (gidpost's HDF5 flavour, gated on
 * `MESHIOPLUSPLUS_WITH_HDF5` at configure time). `GidMode::Auto` always
 * returns the same as `Ascii`, since that is what an extension-less/`.post`
 * path resolves to.
 */
MESHIOPLUSPLUS_API bool gid_available(GidMode mode);

/**
 * @brief Whether this build can **read** @p mode.
 *
 * Deliberately distinct from @ref gid_available, which reports the *write*
 * side. Reading needs no gidpost: `Ascii` is always readable, `Binary` needs
 * `MESHIOPLUSPLUS_HAS_ZLIB` (the flavour is a gzip stream), and `Hdf5` needs
 * `MESHIOPLUSPLUS_HAS_HDF5`. `GidMode::Auto` answers as `Ascii`, since that is
 * what an unrecognized extension resolves to.
 *
 * @note `gid_available`'s meaning is the one it shipped with in v10.18.0 and
 * is deliberately unchanged — silently redefining it to mean "readable" would
 * be a behaviour change for an existing consumer.
 */
MESHIOPLUSPLUS_API bool gid_readable(GidMode mode);

/**
 * @brief The CMake flags that would enable @p mode, for an actionable
 * error message.
 */
MESHIOPLUSPLUS_API std::string gid_build_option(GidMode mode);

/**
 * @brief Write a mesh as a GiD postprocess file (set).
 *
 * @param rPath filesystem path — any of `.post.msh`, `.post.res`,
 *        `.post.bin`, `.post.h5`, or an arbitrary path when @p mode is not
 *        `Auto`. For `GidMode::Ascii` this names (or derives) the sibling
 *        `.post.msh`/`.post.res` pair; for `Binary`/`Hdf5` it names the
 *        single output file directly.
 * @param rMesh the mesh to write.
 * @param mode which flavour to emit; `Auto` (default) infers it from
 *        `rPath`'s extension.
 * @param rAnalysisName the GiD "analysis name" every result is grouped
 *        under (gidpost has no sane default of its own).
 * @param stepValue the single time/load step every result is written at.
 * @throws WriteError when this build cannot emit @p mode (naming
 *         `-DMESHIOPLUSPLUS_WITH_GIDPOST=ON` and
 *         `-DMESHIOPLUSPLUS_WITH_ZLIB=ON`, plus `-DMESHIOPLUSPLUS_WITH_HDF5=ON`
 *         for `Hdf5`), on an unmappable cell type (see the cell-mapping
 *         section above), when the mesh exceeds `INT_MAX` points, cells, or
 *         node indices (gidpost's connectivity is 32-bit), or when gidpost
 *         itself reports a failure.
 */
MESHIOPLUSPLUS_API void write_gid(const std::string& rPath, const Mesh& rMesh,
                                  GidMode mode = GidMode::Auto,
                                  const std::string& rAnalysisName = "meshio++",
                                  double stepValue = 1.0);

/**
 * @brief Read a GiD postprocess file.
 *
 * The flavour is resolved from @p rPath's extension and then **confirmed
 * against the leading bytes**, because the extension alone cannot tell: a
 * `.post.msh` may legitimately be gzipped (gidpost's `GiD_PostAsciiZipped`
 * writes the same ASCII text through `gzprintf`), and a `.post.bin` is a gzip
 * stream wrapping a binary record layout.
 *
 * **Sibling policy** (the `triangle` `.node`/`.ele` precedent): for the ASCII
 * flavour the geometry file `<stem>.post.msh` is **mandatory** and the results
 * file `<stem>.post.res` is **optional** — a mesh with no results file reads
 * back as geometry only. Passing the `.post.res` path directly derives and
 * reads the `.post.msh`; *its* absence is an error, since results alone carry
 * no geometry.
 *
 * Real-world files differ from meshio++'s own output in two ways this reader
 * handles deliberately, neither of which our writer can produce: the full node
 * table may be **repeated in every `MESH` block** (rather than written once
 * with empty blocks thereafter), and element ids may **restart at 1 per
 * block** (rather than being globally unique). Node ids are therefore
 * accumulated into one global table de-duplicated by id, and element ids are
 * tracked per block.
 *
 * @param rPath the file to read (any of the four `.post.*` spellings).
 * @param rOptions selective-read options. `mTimeStep` selects one step of a
 *        multi-step results file and is honoured natively (it cannot be
 *        emulated after the fact); `mMmap` is honoured for the ASCII flavour;
 *        the narrowing options are left to the shared caller-side filter.
 * @return the read Mesh.
 * @throws ReadError on malformed input, on a cell type with no verified GiD
 *         ordering (`hexahedron27`/`wedge15`/`pyramid13`) or no meshio++
 *         mapping (`Sphere`/`Circle`), on an out-of-range time step, or when
 *         this build cannot read the resolved flavour (naming the missing
 *         CMake flag).
 */
MESHIOPLUSPLUS_API Mesh read_gid(const std::string& rPath, const ReadOptions& rOptions = {});

/**
 * @brief Summarize a GiD postprocess file without materializing its arrays.
 *
 * Reports point/cell counts and block shapes from the `MESH` headers, and the
 * available step values from the results file. **Declines** — by throwing
 * `ReadError` — for anything it cannot summarize cheaply, which costs the
 * caller a slower full read rather than a failure (`registry_read_metadata`
 * catches it and falls back).
 */
MESHIOPLUSPLUS_API MeshMetadata read_gid_metadata(const std::string& rPath,
                                                  const ReadOptions& rOptions = {});

}  // namespace meshioplusplus
