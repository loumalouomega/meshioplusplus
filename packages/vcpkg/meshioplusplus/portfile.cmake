# vcpkg overlay port for the meshio++ C API (libmeshioplusplus).
#
# Packages the installable C API only (-DMESHIOPLUSPLUS_BUILD_C_API=ON,
# -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF) -- the same config-package
# (meshioplusplus::meshioplusplus) + pkg-config the standalone `cmake --install`
# produces. Eigen (a git submodule, absent from the release tarball) stays off,
# so the MED transpose uses the hand-written fallback; pugixml and gidpost
# (the GiD postprocess writer's backend) are both vendored hardcopies and so,
# unlike Eigen, ARE present in the release tarball this port pulls.
#
# SHA512 must be refreshed on every release tag: it is the hash of the
# https://github.com/loumalouomega/meshioplusplus/archive/v${VERSION}.tar.gz
# tarball. `vcpkg install meshioplusplus --overlay-ports=packages/vcpkg` prints the
# expected value on mismatch; the packages.yml CI workflow computes it on tag.
vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO loumalouomega/meshioplusplus
    REF "v${VERSION}"
    SHA512 0
    HEAD_REF main
)

vcpkg_check_features(
    OUT_FEATURE_OPTIONS FEATURE_OPTIONS
    FEATURES
        hdf5    MESHIOPLUSPLUS_WITH_HDF5
        netcdf  MESHIOPLUSPLUS_WITH_NETCDF
        zlib    MESHIOPLUSPLUS_WITH_ZLIB
        zstd    MESHIOPLUSPLUS_WITH_ZSTD
        lz4     MESHIOPLUSPLUS_WITH_LZ4
        kahip   MESHIOPLUSPLUS_WITH_KAHIP
        cgnslib MESHIOPLUSPLUS_WITH_CGNSLIB
        gidpost MESHIOPLUSPLUS_WITH_GIDPOST
        cxx-api MESHIOPLUSPLUS_INSTALL_CPP
)

# Which mesh backends the cxx-api feature builds. vcpkg features are booleans,
# not an enum, so the backend set is expressed as one feature per non-default
# backend rather than as a single MESHIOPLUSPLUS_MESH_BACKEND value; MESHIO is
# always included because it is the default meshioplusplus::core points at.
set(MESHIOPLUSPLUS_CXX_BACKENDS "MESHIO")
if("cxx-api-native" IN_LIST FEATURES)
    list(APPEND MESHIOPLUSPLUS_CXX_BACKENDS "NATIVE")
endif()
if("cxx-api-kratos" IN_LIST FEATURES)
    list(APPEND MESHIOPLUSPLUS_CXX_BACKENDS "KRATOS")
endif()
string(REPLACE ";" "\\;" MESHIOPLUSPLUS_CXX_BACKENDS_ARG "${MESHIOPLUSPLUS_CXX_BACKENDS}")

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DMESHIOPLUSPLUS_BUILD_C_API=ON
        -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF
        -DMESHIOPLUSPLUS_WITH_EIGEN=OFF
        # Same reason as Eigen: nlohmann/json is a git submodule absent from the
        # release tarball, so the pipeline JSON entry points raise by name.
        -DMESHIOPLUSPLUS_WITH_JSON=OFF
        "-DMESHIOPLUSPLUS_INSTALL_CPP_BACKENDS=${MESHIOPLUSPLUS_CXX_BACKENDS_ARG}"
        ${FEATURE_OPTIONS}
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME meshioplusplus CONFIG_PATH lib/cmake/meshioplusplus)
vcpkg_fixup_pkgconfig()
vcpkg_copy_pdbs()

# No headers belong in the debug tree in either configuration -- vcpkg installs
# one copy under the release prefix and both builds share it.
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

file(INSTALL "${CMAKE_CURRENT_LIST_DIR}/usage"
     DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
