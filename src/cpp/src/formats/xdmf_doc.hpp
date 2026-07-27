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
 * @file formats/xdmf_doc.hpp
 * @brief The XDMF document-structure helpers `read_xdmf` and
 * `XdmfTimeSeriesWriter` both need: locating a temporal collection and its
 * static mesh grid, and recovering the point/cell counts from a grid already
 * on disk.
 *
 * @note This is a **private** header, deliberately under `src/cpp/src/` rather
 * than `src/cpp/include/`, and the first of its kind in the tree (the nearest
 * precedent is `src/cpp/cli/view_payload.hpp`). It cannot live in
 * `detail/xdmf_common.hpp`, which is *installed*: these declarations name
 * `pugi::xml_node` by value, pugixml is a build-only vendored dependency that
 * no installed header may include, and a forward declaration does not help
 * because `XdmfDoc` stores a `std::vector<pugi::xml_node>` (a complete type is
 * required at class definition). The amalgamator picks it up automatically --
 * it resolves a quoted include relative to the including file's directory and
 * emits any header a `.cpp` needs into the implementation prelude.
 *
 * Sharing matters here, and not only on principle: through v9.1.0 the transient
 * writer's append path re-transcribed a *weaker* copy of `xdmf_resolve` that
 * skipped the `Version` check and recognised a static grid only if it was
 * literally named `"mesh"`, so appending to any other producer's series
 * silently appended a second static grid.
 */

// System includes
#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

// External includes
#include "pugixml.hpp"

// Project includes
#include "meshioplusplus/exceptions.hpp"

namespace meshioplusplus {
namespace xdmfdetail {

// Helpers here are `xdmf_`-prefixed because the single-header amalgamation
// concatenates every src/cpp/src/**.cpp into one translation unit.

/**
 * @brief Parse a `Dimensions="..."` attribute into its extents.
 * @param rS The attribute value, whitespace-separated integers.
 * @return The extents, outermost first; empty when @p rS holds no integer.
 */
inline std::vector<std::size_t> xdmf_parse_dims(const std::string& rS) {
    std::vector<std::size_t> dims;
    std::istringstream iss(rS);
    std::int64_t v;
    while (iss >> v)
        dims.push_back(static_cast<std::size_t>(v));
    return dims;
}

/**
 * @brief A resolved XDMF document: its static mesh grid, and its time steps.
 *
 * `mSteps` empty means "a plain single-grid file", which is the historical case
 * and takes the historical code path unchanged.
 */
struct XdmfDoc {
    pugi::xml_node mMeshGrid;            ///< Grid holding `<Topology>`/`<Geometry>`.
    pugi::xml_node mCollection;          ///< The temporal collection, if the file has one.
    std::vector<pugi::xml_node> mSteps;  ///< Temporal-collection children, in file order.
};

/**
 * @brief Validate the document and locate its mesh grid (and time steps).
 *
 * Shared by the mesh reader, the metadata reader and the transient writer's
 * append path, so a summary can never accept a file `read_xdmf` would reject
 * and an append can never disagree with a read about which grid is the mesh. A
 * temporal collection (`GridType="Collection" CollectionType="Temporal"`, what
 * `XdmfTimeSeriesWriter` emits) stores the static geometry once in a sibling
 * `Uniform` grid and one `<Grid>` per step inside the collection; the step
 * grids reference the geometry through an `xi:include` this reader ignores,
 * resolving it structurally instead -- the same thing the Python
 * `TimeSeriesReader` does, and for the same reason: a generic XInclude pass
 * would have to implement XPointer.
 *
 * @param rDoc The parsed document.
 * @return The resolved grids.
 * @throws ReadError if the document is not an XDMF 3 document with a `<Grid>`,
 *         or if a temporal collection carries no mesh grid at all.
 */
inline XdmfDoc xdmf_resolve(const pugi::xml_document& rDoc) {
    pugi::xml_node root = rDoc.child("Xdmf");
    if (!root)
        throw ReadError("XDMF: missing <Xdmf> root");
    std::string version = root.attribute("Version").value();
    if (!version.empty() && version[0] != '3')
        throw ReadError("XDMF: only version 3 handled by the C++ core");

    pugi::xml_node domain = root.child("Domain");
    pugi::xml_node first, uniform, collection;
    for (pugi::xml_node g : domain.children("Grid")) {
        if (!first)
            first = g;
        const std::string gtype = g.attribute("GridType").value();
        if (gtype == "Collection" &&
            std::string(g.attribute("CollectionType").value()) == "Temporal") {
            if (!collection)
                collection = g;
        } else if (gtype == "Uniform" && !uniform) {
            uniform = g;
        }
    }
    if (!first)
        throw ReadError("XDMF: missing <Grid>");

    XdmfDoc out;
    out.mCollection = collection;
    if (!collection) {
        out.mMeshGrid = first;
        return out;
    }
    for (pugi::xml_node s : collection.children("Grid"))
        out.mSteps.push_back(s);
    out.mMeshGrid = uniform;
    if (!out.mMeshGrid) {
        // No sibling mesh grid: take the first uniform grid inside the
        // collection, which is where a writer that repeats the geometry per
        // step puts it.
        for (pugi::xml_node s : out.mSteps) {
            if (std::string(s.attribute("GridType").value()) == "Uniform") {
                out.mMeshGrid = s;
                break;
            }
        }
    }
    if (!out.mMeshGrid)
        throw ReadError("XDMF: temporal collection carries no mesh grid");
    return out;
}

/** @brief Point/cell counts recovered from a mesh grid that is already on disk. */
struct XdmfGridCounts {
    std::size_t mNumPoints = 0;
    std::size_t mNumCells = 0;
    bool mPointsKnown = false;
    bool mCellsKnown = false;
};

/**
 * @brief Recover a grid's point and cell counts from its attributes alone.
 *
 * No payload is read: every `<DataItem>` declares its shape, so this costs one
 * attribute lookup even for an HDF series (the `.h5` is never opened).
 *
 * Cells come from `<Topology>`'s `NumberOfElements` **first**, and only then
 * from the child `<DataItem>`'s `Dimensions`. The order is load-bearing: for a
 * `TopologyType="Mixed"` grid the `DataItem` holds one flat packed array whose
 * length is not the cell count, and `NumberOfElements` is the only correct
 * source. That is also why this cannot simply call `read_xdmf_metadata`, which
 * declines Mixed outright.
 *
 * @param rMeshGrid The static mesh grid, e.g. `xdmf_resolve(doc).mMeshGrid`.
 * @return The counts, each flagged with whether anything readable was found.
 */
inline XdmfGridCounts xdmf_grid_counts(const pugi::xml_node& rMeshGrid) {
    XdmfGridCounts out;
    if (!rMeshGrid)
        return out;

    if (pugi::xml_node geo = rMeshGrid.child("Geometry")) {
        const std::vector<std::size_t> dims =
            xdmf_parse_dims(geo.child("DataItem").attribute("Dimensions").value());
        if (!dims.empty()) {
            out.mNumPoints = dims[0];
            out.mPointsKnown = true;
        }
    }

    if (pugi::xml_node topo = rMeshGrid.child("Topology")) {
        if (pugi::xml_attribute n = topo.attribute("NumberOfElements")) {
            out.mNumCells = static_cast<std::size_t>(n.as_ullong());
            out.mCellsKnown = true;
        } else {
            const std::vector<std::size_t> dims =
                xdmf_parse_dims(topo.child("DataItem").attribute("Dimensions").value());
            if (!dims.empty()) {
                out.mNumCells = dims[0];
                out.mCellsKnown = true;
            }
        }
    } else {
        // A grid with no <Topology> is points-only, which WritePointsCells
        // emits for a mesh with zero cell blocks: zero cells, known exactly.
        out.mCellsKnown = true;
    }
    return out;
}

}  // namespace xdmfdetail
}  // namespace meshioplusplus
