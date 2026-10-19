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

// STEP and IGES reading through OpenCASCADE. See formats/step.hpp and
// formats/iges.hpp.

// System includes
#include <cmath>
#include <stdexcept>
#include <string>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/iges.hpp"
#include "meshioplusplus/formats/step.hpp"

#ifdef MESHIOPLUSPLUS_HAS_OCCT
// Project includes (private, not installed)
#include "../detail/library_lock.hpp"
#include "../detail/occt_brep.hpp"
#endif

namespace meshioplusplus {

namespace {

#ifdef MESHIOPLUSPLUS_HAS_OCCT
Mesh cad_read_shape(const std::string& rPath, bool Iges, BrepInfo& rInfo,
                    const CadReadOptions& rOptions) {
    const char* fmt = Iges ? "iges" : "step";
    if (!std::isfinite(rOptions.mDeflection) || !(rOptions.mAngularDeflection > 0.0) ||
        !std::isfinite(rOptions.mAngularDeflection))
        throw std::invalid_argument(std::string("meshio++: ") + fmt +
                                    ": the deflection must be finite and the angular "
                                    "deflection positive");
    detail::LibraryLock lock;
    detail::OcctQuiet quiet;
    detail::OcctNurbs nurbs = detail::occt_nurbs(detail::occt_read_cad_file(rPath, Iges, fmt));
    Mesh mesh = detail::occt_tessellate(nurbs.mShape, nurbs.mInfo, rOptions.mDeflection,
                                        rOptions.mAngularDeflection);
    rInfo = std::move(nurbs.mInfo);
    return mesh;
}
#else
[[noreturn]] void cad_read_unsupported(const char* pFormat) {
    throw Unsupported(std::string("meshio++: ") + pFormat +
                      ": reading needs the OpenCASCADE backend, which the PyPI wheels and "
                      "the release binaries carry; rebuild with -DMESHIOPLUSPLUS_WITH_OCCT=ON "
                      "(OpenCASCADE 7.6 or later, build/build-occt.sh builds one)");
}
#endif

}  // namespace

bool cad_read_has_occt() noexcept {
#ifdef MESHIOPLUSPLUS_HAS_OCCT
    return true;
#else
    return false;
#endif
}

Mesh read_step(const std::string& rPath) {
    BrepInfo info;
    return read_step(rPath, info);
}

Mesh read_step(const std::string& rPath, BrepInfo& rInfo, const CadReadOptions& rOptions) {
#ifdef MESHIOPLUSPLUS_HAS_OCCT
    return cad_read_shape(rPath, false, rInfo, rOptions);
#else
    (void)rPath;
    (void)rInfo;
    (void)rOptions;
    cad_read_unsupported("step");
#endif
}

Mesh read_iges(const std::string& rPath) {
    BrepInfo info;
    return read_iges(rPath, info);
}

Mesh read_iges(const std::string& rPath, BrepInfo& rInfo, const CadReadOptions& rOptions) {
#ifdef MESHIOPLUSPLUS_HAS_OCCT
    return cad_read_shape(rPath, true, rInfo, rOptions);
#else
    (void)rPath;
    (void)rInfo;
    (void)rOptions;
    cad_read_unsupported("iges");
#endif
}

}  // namespace meshioplusplus
