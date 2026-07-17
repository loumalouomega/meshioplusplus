# vcpkg overlay port for the meshio++ C API (libmeshioplusplus).
#
# Packages the installable C API only (-DMESHIOPLUSPLUS_BUILD_C_API=ON,
# -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF) -- the same config-package
# (meshioplusplus::meshioplusplus) + pkg-config the standalone `cmake --install`
# produces. Eigen (a git submodule, absent from the release tarball) stays off,
# so the MED transpose uses the hand-written fallback; pugixml is vendored.
#
# SHA512 must be refreshed on every release tag: it is the hash of the
# https://github.com/loumalouomega/meshioplusplus/archive/v${VERSION}.tar.gz
# tarball. `vcpkg install meshioplusplus --overlay-ports=ports` prints the
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
)

vcpkg_cmake_configure(
    SOURCE_PATH "${SOURCE_PATH}"
    OPTIONS
        -DMESHIOPLUSPLUS_BUILD_C_API=ON
        -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF
        -DMESHIOPLUSPLUS_WITH_EIGEN=OFF
        ${FEATURE_OPTIONS}
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME meshioplusplus CONFIG_PATH lib/cmake/meshioplusplus)
vcpkg_fixup_pkgconfig()
vcpkg_copy_pdbs()

# The library is shared-only; no headers belong in the debug tree.
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")

file(INSTALL "${CMAKE_CURRENT_LIST_DIR}/usage"
     DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}")
vcpkg_install_copyright(FILE_LIST "${SOURCE_PATH}/LICENSE")
