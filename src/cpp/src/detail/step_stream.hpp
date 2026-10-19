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
 * @file detail/step_stream.hpp
 * @brief The AP214 file around a STEP B-rep: the header, the product
 * structure, the units and the geometric context, then the entities a writer
 * adds one by one (sequential ids), and the shape representation that
 * collects its solids. Shared by the faceted `step` writer and the
 * kernel-free B-rep writer (`step_brep`), so both write the same frame.
 *
 * A **core-private** header (the `crease_edges.hpp` precedent): no installed
 * header names it, and it adds nothing to the API or the ABI.
 */

// System includes
#include <cstdint>
#include <string>
#include <vector>

namespace meshioplusplus {
namespace detail {

/// A STEP REAL: `%.17g`, with the decimal point the grammar requires and an
/// upper-case exponent ("3" -> "3.", "1e-05" -> "1.E-05").
std::string step_real(double V);

/// A STEP string body: printable ASCII, apostrophes doubled.
std::string step_string(const std::string& rText);

/// The entities of one STEP file, written to memory and then to disk.
class StepStream {
public:
    /**
     * @brief Start a file: the header (with the provenance comment), the
     * product named after @p rPath's basename, millimetres and radians, and
     * the world placement.
     * @param rDescription the `FILE_DESCRIPTION` text.
     */
    StepStream(const std::string& rPath, const std::string& rDescription);

    /// Add an entity; returns its id.
    std::int64_t Ent(const std::string& rBody);
    /// `#id`.
    static std::string Ref(std::int64_t Id);
    /// `(x,y,z)` of three REALs.
    static std::string Triple(const double* pV);
    /// Close the representation over @p rSolids (MANIFOLD_SOLID_BREP ids) and
    /// write the file (nothing, for an empty path).
    /// @throws WriteError when the file cannot be opened.
    void Finish(const std::vector<std::int64_t>& rSolids);

private:
    std::string mPath;
    std::string mOut;
    std::int64_t mNext = 1;
    std::int64_t mShape = 0;    // PRODUCT_DEFINITION_SHAPE
    std::int64_t mContext = 0;  // the geometric representation context
    std::int64_t mOrigin = 0;   // the world AXIS2_PLACEMENT_3D
};

}  // namespace detail
}  // namespace meshioplusplus
