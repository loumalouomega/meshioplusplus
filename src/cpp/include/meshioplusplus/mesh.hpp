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
 * @file mesh.hpp
 * @brief Compile-time mesh-backend dispatch: selects which in-memory mesh
 * structure `meshioplusplus::Mesh` is.
 *
 * meshio++ has three interchangeable mesh backends, selected at build time
 * by the `MESHIOPLUSPLUS_MESH_BACKEND` CMake option (exactly one of the
 * `MESHIOPLUSPLUS_MESH_BACKEND_*` macros is defined — mirroring the
 * `MESHIOPLUSPLUS_PARALLEL_*` parallel-backend pattern in `parallel.hpp`):
 *
 *  - **MESHIO** (`backends/meshio_mesh.hpp`, the default): the
 *    meshio-mirroring `Mesh`/`CellBlock` over dtype-erased `NDArray`s.
 *    Required when the pybind11 extension is built — the zero-copy numpy
 *    boundary (`bindings/python/np_conversions.hpp`) is written against it.
 *  - **NATIVE** (`backends/native_mesh.hpp`): canonical statically-typed
 *    storage — Float64 points, Int64 connectivity, `CellType` enum,
 *    CSR-shaped ragged blocks. The fastest pure-C++ consumer surface; used
 *    by the WebAssembly build.
 *  - **KRATOS** (`backends/kratos_mesh.hpp`): a Kratos-Multiphysics-style
 *    `ModelPart` (Nodes/Elements/Conditions/SubModelParts) behind the same
 *    API, for near-costless exchange with Kratos (see `kratos_bridge.hpp`).
 *
 * All three implement the uniform format-facing API documented in
 * `mesh_api.hpp`; format code compiles unchanged under any of them. To add
 * a backend: add one CMake branch defining a new
 * `MESHIOPLUSPLUS_MESH_BACKEND_<NAME>` macro, one `#elif` below, and a
 * `backends/<name>_mesh.hpp` implementing the API.
 */

// Exactly one backend macro may be defined. Two would silently take the first
// branch below and ignore the rest, giving this TU a different `Mesh` from the
// library it links against -- see detail/mesh_backend_check.hpp, which also
// catches the opposite mistake (none defined, so the MESHIO default applies to
// a TU linking a NATIVE/KRATOS build).
#if (defined(MESHIOPLUSPLUS_MESH_BACKEND_MESHIO) + defined(MESHIOPLUSPLUS_MESH_BACKEND_NATIVE) + \
     defined(MESHIOPLUSPLUS_MESH_BACKEND_KRATOS)) > 1
#error \
    "meshio++: more than one MESHIOPLUSPLUS_MESH_BACKEND_* macro is defined. Define exactly the one the library you link was built with (CMake consumers inherit it automatically from meshioplusplus::core*)."
#endif

// Project includes
#include "mesh_api.hpp"
// Referenced by nobody here on purpose: it plants the link-time sentinel in
// every TU that includes mesh.hpp, which is exactly the set that could get the
// backend wrong. IWYU pragma: keep
#include "detail/mesh_backend_check.hpp"
// Likewise, for the second axis: same TU set, same technique, different
// question -- "were these headers and this library the same meshio++?" rather
// than "the same backend?". See detail/abi_version_check.hpp. IWYU pragma: keep
#include "detail/abi_version_check.hpp"

#if defined(MESHIOPLUSPLUS_MESH_BACKEND_NATIVE)
#include "backends/native_mesh.hpp"
namespace meshioplusplus {
using Mesh = NativeMesh;
}
#elif defined(MESHIOPLUSPLUS_MESH_BACKEND_KRATOS)
#include "backends/kratos_mesh.hpp"
namespace meshioplusplus {
using Mesh = KratosMesh;
}
#else  // MESHIOPLUSPLUS_MESH_BACKEND_MESHIO (and the no-macro default)
#include "backends/meshio_mesh.hpp"
// backends/meshio_mesh.hpp defines `struct Mesh` directly (no alias) so the
// pybind11 binding layer sees literally the same type as before the
// backends existed.
#endif
