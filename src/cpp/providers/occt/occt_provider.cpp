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

// The OpenCASCADE CAD kernel as a provider: STEP and IGES reading and the
// OpenCASCADE backend of fit_cad. Compiled into the core in a monolithic
// build, and into the meshioplusplus-occt add-on in a split one. See
// providers.hpp.

#ifdef MESHIOPLUSPLUS_HAS_OCCT

// System includes
#include <cmath>
#include <stdexcept>
#include <string>

// External includes (OpenCASCADE)
#include <Standard_Failure.hxx>
#include <Standard_Version.hxx>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/providers.hpp"

// Project includes (private, not installed)
#include "detail/fit_cad_native.hpp"
#include "detail/library_lock.hpp"
#include "occt_brep.hpp"
#include "occt_provider.hpp"

namespace meshioplusplus {

namespace {

Mesh occt_read_shape(const std::string& rPath, bool Iges, BrepInfo& rInfo,
                     const CadReadOptions& rOptions) {
    const char* fmt = Iges ? "iges" : "step";
    if (!std::isfinite(rOptions.mDeflection) || !(rOptions.mAngularDeflection > 0.0) ||
        !std::isfinite(rOptions.mAngularDeflection))
        throw std::invalid_argument(std::string("meshio++: ") + fmt +
                                    ": the deflection must be finite and the angular "
                                    "deflection positive");
    detail::LibraryLock lock;
    detail::OcctQuiet quiet;
    try {
        detail::OcctNurbs nurbs = detail::occt_nurbs(detail::occt_read_cad_file(rPath, Iges, fmt));
        Mesh mesh = detail::occt_tessellate(nurbs.mShape, nurbs.mInfo, rOptions.mDeflection,
                                            rOptions.mAngularDeflection);
        rInfo = std::move(nurbs.mInfo);
        return mesh;
    } catch (const Standard_Failure& rErr) {
        // Never OpenCASCADE's own exception, which is no std::exception.
        throw ReadError(std::string("meshio++: ") + fmt + ": OpenCASCADE failed converting '" +
                        rPath + "': " + rErr.GetMessageString());
    }
}

}  // namespace

namespace providers {

CadKernelProvider make_occt_provider() {
    CadKernelProvider provider;
    provider.mName = "occt";
    provider.mVersion = OCC_VERSION_COMPLETE;
    provider.mReadStep = [](const std::string& rPath, BrepInfo& rInfo,
                            const CadReadOptions& rOptions) {
        return occt_read_shape(rPath, false, rInfo, rOptions);
    };
    provider.mReadIges = [](const std::string& rPath, BrepInfo& rInfo,
                            const CadReadOptions& rOptions) {
        return occt_read_shape(rPath, true, rInfo, rOptions);
    };
    provider.mFitCad = detail::occt_fit_cad;
    return provider;
}

void register_occt() {
    register_cad_kernel(make_occt_provider());
}

}  // namespace providers

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_OCCT
