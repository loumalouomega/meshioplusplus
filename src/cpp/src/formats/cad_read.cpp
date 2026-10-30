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

// STEP and IGES reading through the registered CAD kernel (OpenCASCADE). See
// formats/step.hpp, formats/iges.hpp and providers.hpp; the OpenCASCADE kernel
// is providers/occt/.

// System includes
#include <functional>
#include <string>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/formats/iges.hpp"
#include "meshioplusplus/formats/step.hpp"
#include "meshioplusplus/providers.hpp"

namespace meshioplusplus {

namespace {

[[noreturn]] void cad_read_unsupported(const char* pFormat) {
    throw Unsupported(std::string("meshio++: ") + pFormat +
                      ": reading needs a CAD kernel (OpenCASCADE), which the release binaries "
                      "carry and the Python wheels provide through the meshioplusplus-occt "
                      "package (pip install meshioplusplus[occt]); a source build needs "
                      "-DMESHIOPLUSPLUS_WITH_OCCT=ON (OpenCASCADE 7.6 or later, "
                      "build/build-occt.sh builds one)");
}

Mesh cad_read(const char* pFormat, bool Iges, const std::string& rPath, BrepInfo& rInfo,
              const CadReadOptions& rOptions) {
    const auto p_kernel = cad_kernel();
    if (!p_kernel || !(Iges ? p_kernel->mReadIges : p_kernel->mReadStep))
        cad_read_unsupported(pFormat);
    return (Iges ? p_kernel->mReadIges : p_kernel->mReadStep)(rPath, rInfo, rOptions);
}

}  // namespace

bool cad_read_has_occt() noexcept {
    return cad_kernel() != nullptr;
}

Mesh read_step(const std::string& rPath) {
    BrepInfo info;
    return read_step(rPath, info);
}

Mesh read_step(const std::string& rPath, BrepInfo& rInfo, const CadReadOptions& rOptions) {
    return cad_read("step", false, rPath, rInfo, rOptions);
}

Mesh read_iges(const std::string& rPath) {
    BrepInfo info;
    return read_iges(rPath, info);
}

Mesh read_iges(const std::string& rPath, BrepInfo& rInfo, const CadReadOptions& rOptions) {
    return cad_read("iges", true, rPath, rInfo, rOptions);
}

}  // namespace meshioplusplus
