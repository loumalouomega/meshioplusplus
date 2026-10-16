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
 * @file marc.hpp
 * @brief MSC Marc input decks (`.dat`) and formatted post files (`.t19`), read.
 *
 *  - An input deck's model definition makes the mesh: `COORDINATES` (the last
 *    definition of a node wins) and `CONNECTIVITY`, in fixed fields (5/10
 *    columns, twice as wide after the `EXTENDED` parameter, reals with or
 *    without an exponent letter) or free (comma) format. Every Marc element
 *    type lists its nodes in meshio++'s order; an 8-node brick with repeated
 *    nodes is a wedge, pyramid or tetrahedron. `"marc:element"` and
 *    `"marc:type"` cell data keep the element number and Marc type; elements of
 *    a type meshio++ has no cell for are skipped with a warning.
 *  - `DEFINE ELEMENT SET` / `DEFINE NODE SET` become cell and point regions:
 *    numbers, `a TO b [BY c]` ranges and earlier sets' names, combined by `AND`,
 *    `EXCEPT` and `INTERSECT`, continued by a trailing `C`. Other set kinds
 *    (edges, faces ...) are skipped with a warning. History definition is not
 *    read.
 *  - A `.t19` post file (Volume D, PLDUMP2000) gives the same mesh and sets,
 *    and one step per increment: its time (a frequency or buckling factor in a
 *    modal, harmonic or buckling analysis) is `field_data["meshio:time"]` with
 *    `"marc:increment"` and `"marc:subincrement"`; nodal vectors are point data
 *    named as the file names them; element post codes are cell data per
 *    integration point, `(cells, points * components)` flattened point-major
 *    with `field_data["marc:layout:<name>"]` = `[points, components]`
 *    (`(cells[, components])` and no layout when an element has one), a
 *    tensor's six codes combined as `xx yy zz xy yz zx`.
 *    Increments that remesh the model are refused.
 * See doc/formats/marc.md.
 */

// System includes
#include <cstddef>
#include <memory>
#include <string>
#include <string_view>

// Project includes
#include "meshioplusplus/export.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/read_options.hpp"

namespace meshioplusplus {

/**
 * @brief Whether `Head` (the start of a `.dat` file) is a Marc input deck: its
 *        first keyword line is a Marc parameter with no `=` (unlike a Tecplot
 *        header), and a `CONNECTIVITY`, `COORDINATES` or `END` line follows.
 */
MESHIOPLUSPLUS_API bool is_marc_deck(std::string_view Head);

/**
 * @brief Read an MSC Marc input deck.
 * @throws ReadError for a file that is not a Marc deck, a malformed field, an
 *         element naming an undefined node, or a set naming an undefined member
 */
MESHIOPLUSPLUS_API Mesh read_marc(const std::string& rPath);

/**
 * @brief Write an MSC Marc input deck (`.dat`) in `EXTENDED` format.
 *
 * The parameter section (`TITLE`, `SIZING`, one `ELEMENTS` line per type,
 * `END`), then `CONNECTIVITY`, `COORDINATES`, a `DEFINE NODE SET` per point
 * region and a `DEFINE ELEMENT SET` per cell region, and `END OPTION`. A
 * cell's Marc type is its `marc:type` when that fits it, else a default (the
 * solids and shells in 3-D, the plane-strain elements in a planar mesh; a
 * pyramid as a degenerate brick); element numbers are `marc:element` when
 * positive and unique. Face and edge sets come back from their `marc:` field
 * data. `.dat` is Tecplot's for a write by extension: name the format.
 * Other cells, data and side regions are dropped with a warning.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError for a mesh with no points or points of dimension above 3
 * @note Since v16.17.0. The Python twin writes the same bytes.
 */
MESHIOPLUSPLUS_API void write_marc(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read one increment of an MSC Marc formatted post file.
 * @param rOptions `mTimeStep` picks the increment (negative counts from the
 *        end); `mPointsOnly`/`mDataArrays` narrow the data read
 * @throws ReadError for a file that is not a post file, a remeshing increment,
 *         or an out-of-range time step
 */
/**
 * @brief Write one increment's MSC Marc formatted post file (`.t19`): the mesh,
 *        its sets and, when present, element and nodal results, at
 *        `field_data["meshio:time"]`. The model is the input deck writer's
 *        (cells with a Marc element type, `marc:element` numbers, point and
 *        cell regions as node and element sets). Point arrays of 1, 2, 3, 6 or 9
 *        components are nodal results; cell arrays are element post codes (per
 *        integration point where `marc:layout:<name>` says so). The Python twin
 *        writes the same bytes.
 * @throws WriteError when the mesh has no cell with a Marc element type
 */
MESHIOPLUSPLUS_API void write_marc_t19(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief An increment per mesh in one post file, written by `Finalize`; every
 *        mesh must have the first one's points, cells and arrays. The Python
 *        twin is `meshioplusplus.marc.SeriesWriter`.
 */
class MESHIOPLUSPLUS_API MarcT19SeriesWriter {
public:
    explicit MarcT19SeriesWriter(const std::string& rPath);
    ~MarcT19SeriesWriter();

    MarcT19SeriesWriter(const MarcT19SeriesWriter&) = delete;
    MarcT19SeriesWriter& operator=(const MarcT19SeriesWriter&) = delete;
    MarcT19SeriesWriter(MarcT19SeriesWriter&&) noexcept;
    MarcT19SeriesWriter& operator=(MarcT19SeriesWriter&&) noexcept;

    /// Add one increment. @throws WriteError on a mismatch with the first one.
    void Write(double Time, const Mesh& rMesh);

    /// The number of increments added so far.
    std::size_t NumSteps() const noexcept;

    /// Write the file. Idempotent. @throws WriteError when no increment was added.
    void Finalize();

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
};

MESHIOPLUSPLUS_API Mesh read_marc_t19(const std::string& rPath, const ReadOptions& rOptions = {});

/**
 * @brief Summarise an MSC Marc formatted post file, its increment times included.
 */
MESHIOPLUSPLUS_API MeshMetadata read_marc_t19_metadata(const std::string& rPath,
                                                       const ReadOptions& rOptions = {});

}  // namespace meshioplusplus
