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
 * validates instead.
 *
 * **Which type an array is written as** is chosen in one of two ways:
 *
 * 1. **Declared**, via a `field_data` entry `"gid:result_type:<name>"`
 *    (`kGidResultTypePrefix`) holding a `GidResultType` value. The count is
 *    checked against that type's legal counts and an illegal one is a
 *    `WriteError` naming the array — never a silent fallback.
 * 2. **Inferred**, when no declaration is present: 1 component →
 *    `GiD_Scalar`; 2 or 3 → `GiD_Vector`; anything else splits into that
 *    many named `GiD_Scalar` results (`"<name>_1"` .. `"<name>_k"`),
 *    recorded with a `provenance_note`.
 *
 * An undeclared 6-component array is deliberately **not** inferred as
 * `GiD_Matrix`, even though stress tensors are GiD's canonical use case:
 * `(n,6)` alone is genuinely ambiguous (it could equally be `Matrix:6`,
 * `ComplexMatrix:3` or `ComplexVector:6`), so inferring would silently pick
 * one meaning. Declaring is the way to say which. Note the component order
 * itself needs no conversion: GiD's `Matrix:6` is `Sxx Syy Szz Sxy Syz Sxz`,
 * which is exactly meshio/VTK's symmetric-tensor order.
 *
 * Named regions and multi-step results are not carried (each dropped with a
 * `provenance_note` or `log::warn`); this writes exactly one step
 * (`rStepValue`). `field_data` is not written either, but *is* read, for the
 * `"gid:result_type:*"` keys above.
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
#include <cstddef>
#include <cstdint>
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
    /// Two sibling files like `Ascii`, but gzipped -- gidpost's
    /// `GiD_PostAsciiZipped`, which is the SAME ASCII text through `gzprintf`.
    ///
    /// Reading has always handled it (the reader sniffs the gzip magic and
    /// inflates, since the flavour is textually identical), so this is the
    /// write counterpart alone. It needs no dependency `Ascii` does not:
    /// gidpost hard-requires zlib regardless.
    ///
    /// **`Auto` never resolves to this.** No extension can express "zipped" --
    /// a gzipped file still ends `.post.msh` -- so the mode is explicit-only.
    /// Inferring it would also change what every existing `.post.msh` write
    /// produces, which is exactly what `Auto` must not do.
    AsciiZipped = 4,
};

/**
 * @brief What kind of quantity a result array holds, in GiD's own vocabulary.
 *
 * meshio++'s `Mesh` has no way to say "this array is a symmetric tensor" or
 * "this array is complex" -- `NDArray` has neither a complex dtype nor a
 * string dtype -- so a caller declares it out of band, through a `field_data`
 * entry named `kGidResultTypePrefix + "<array name>"` holding one of these
 * values. That mechanism was chosen over a typed side-channel struct
 * (`MedInfo`'s shape) precisely because the registry's `(path, mesh)` writers
 * cannot carry a side channel, so it would have been invisible from the CLI,
 * WASM, the C API and Fortran; a `field_data` key reaches all of them with no
 * per-binding code, exactly as `gmsh:physical` and `med:num` already do.
 *
 * The enumerator values deliberately equal gidpost's own `GiD_ResultType`, so
 * the writer is a `static_cast`; `gid.cpp` `static_assert`s every one of them
 * against it (the `mio_cell_type` precedent), making drift a compile error.
 * The enum is declared here rather than reusing `GiD_ResultType` because this
 * header must stay free of gidpost, which the *reader* is deliberately built
 * without -- reading ASCII GiD needs neither gidpost nor zlib.
 *
 * ## Component counts and order
 *
 * Counts are gidpost's (`_ResultTypeInfo`); the orders are quoted from
 * CIMNE's GiD Customization Manual. Values are stored **verbatim in GiD's
 * order** -- meshio++ does not reinterpret them, having no canonical complex
 * or tensor layout of its own to convert to.
 *
 * | Type | Legal counts | Order |
 * |---|---|---|
 * | `Scalar` | 1 | the value |
 * | `Vector` | 2, 3, 4 | X, Y, Z, \|V\| (4th = signed modulus) |
 * | `Matrix` | 3, 6 | 3: Sxx Syy Sxy / 6: Sxx Syy Szz Sxy Syz Sxz |
 * | `PlainDeformationMatrix` | 4 | Sxx Syy Sxy Szz |
 * | `MainMatrix` | 12 | Si Sii Siii, then the three eigenvectors |
 * | `LocalAxes` | 3 | euler_ang_1..3 |
 * | `ComplexScalar` | 2 | real, imag |
 * | `ComplexVector` | 4, 6 | **interleaved**: x_re x_im y_re y_im [z_re z_im] |
 * | `ComplexMatrix` | 6, 12 | **blocked**: every real, then every imag |
 *
 * Two of those are worth not rediscovering. `Matrix:6` is **already**
 * meshio/VTK's symmetric-tensor order, so a stress tensor needs no
 * permutation in either direction. And `ComplexVector` interleaves real and
 * imaginary parts per component while `ComplexMatrix` blocks them -- the same
 * family, opposite conventions -- so neither may be inferred from the other.
 */
enum class GidResultType : std::int64_t {
    Scalar = 0,
    Vector = 1,
    Matrix = 2,
    PlainDeformationMatrix = 3,
    MainMatrix = 4,
    LocalAxes = 5,
    ComplexScalar = 6,
    ComplexVector = 7,
    ComplexMatrix = 8,
};

/**
 * @brief `field_data` key prefix declaring an array's `GidResultType`.
 *
 * The full key is this prefix plus the array's own name, and its value is a
 * single integer `GidResultType`. `field_data` is global rather than
 * per-location, so one key covers an array of that name wherever it appears:
 * if the same name exists in both `point_data` and `cell_data`, the single
 * declaration applies to both and must be legal for both counts.
 */
inline constexpr const char* kGidResultTypePrefix = "gid:result_type:";

/**
 * @brief `field_data` key prefix declaring how many Gauss points per element a
 * `cell_data` array carries.
 *
 * GiD writes a per-element result against a **Gauss-point set** of G points,
 * emitting G rows per element. meshio++'s `cell_data` has no
 * per-point-within-cell axis (the same limit MED's ELNO/ELGA documents), so a
 * G-point, k-component array is stored **flat** as `(ncells, G*k)`, laid out
 * Gauss-point-major -- `[gp0_c0..gp0_ck-1, gp1_c0..gp1_ck-1, ...]`, which is
 * the order GiD's own `Values` rows arrive in, so neither direction re-packs.
 *
 * The full key is this prefix plus the array's own name; its value is a single
 * integer G. It is written **only when it carries information**, i.e. never
 * for `G == 1` -- so an ordinary per-element array is a plain `(ncells, k)`
 * with no declaration at all, and its bytes are unchanged by this mechanism's
 * existence. Without the declaration a `(ncells, 3)` array is genuinely
 * ambiguous (a 3-component vector at one Gauss point, or a scalar at three),
 * which is exactly why the declaration is required rather than inferred.
 *
 * Interacts with `kGidResultTypePrefix`: the result type's legal component
 * counts are checked against **k**, not `G*k`. A `Matrix` (k=6) at G=3 is a
 * `(n, 18)` array whose declared type is still validated as 6 components.
 */
inline constexpr const char* kGidGaussPointsPrefix = "gid:gauss_points:";

/**
 * @brief `field_data` key prefix supplying a Gauss-point set's **natural
 * coordinates**, for counts GiD cannot compute itself.
 *
 * GiD accepts `Natural Coordinates: Internal` -- letting it place the points
 * -- only for specific counts per element type (triangle 1/3/6, quadrilateral
 * 1/4/9, tetrahedra 1/4/10, hexahedra 1/8/27, prism 1/6, pyramid 1/5; line
 * elements accept any count, equally spaced). Any **other** G must declare
 * `Natural Coordinates: Given` and list the points explicitly.
 *
 * The full key is this prefix plus `"<meshio++ cell type>:<G>"` -- keyed by
 * `(cell type, G)` rather than by array name, because that is what a GiD
 * Gauss-point set actually depends on: two arrays sharing a cell block and a
 * G share one set, and would otherwise have to repeat identical coordinates.
 * The value is `G * dim` doubles, point-major, where `dim` is 2 for
 * triangle/quadrilateral and 3 for tetrahedra/hexahedra/prism/pyramid.
 *
 * Coordinate ranges are GiD's own: 0..1 for triangle, tetrahedra and prism;
 * -1..1 for quadrilateral, hexahedra and pyramid. **Line elements cannot use
 * `Given` at all** -- GiD forbids it -- so they are Internal-only, which is no
 * restriction since they already accept any count.
 *
 * Supplying coordinates for a G that GiD *could* have computed is honoured
 * (the set is written `Given`); omitting them for a G it cannot is a
 * `WriteError` naming the type and its legal Internal counts.
 */
inline constexpr const char* kGidGaussCoordsPrefix = "gid:gauss_coords:";

/// The GiD spelling of @p type (`"Matrix"`, `"ComplexVector"`, ...), as it
/// appears in a `.post.res` `Result` header.
MESHIOPLUSPLUS_API const char* gid_result_type_name(GidResultType type);

/**
 * @brief Parses a `GidResultType` from its GiD spelling.
 * @throws std::invalid_argument on an unrecognized name.
 */
MESHIOPLUSPLUS_API GidResultType gid_result_type_from_name(const std::string& rName);

/// Whether @p k is a component count @p type accepts (see the table above).
MESHIOPLUSPLUS_API bool gid_result_dim_is_legal(GidResultType type, std::size_t k);

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
