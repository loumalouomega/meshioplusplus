# Conan 2.x recipe for the meshio++ C API (libmeshioplusplus).
#
# This packages the installable C API only (-DMESHIOPLUSPLUS_BUILD_C_API=ON,
# -DMESHIOPLUSPLUS_BUILD_PYTHON=OFF): the same relocatable CMake config-package
# (`meshioplusplus::meshioplusplus`) + pkg-config file the standalone
# `cmake --install` produces. The Python wheel is published to PyPI separately
# (see pyproject.toml / .github/workflows/wheels.yml) and is unrelated to this.
#
# `version` must track pyproject.toml's `version` and CMakeLists.txt's
# `project(... VERSION ...)` -- bump all three together on a release.
#
# Notes on optional dependencies:
#   * pugixml is vendored in-tree (cpp/third_party/pugixml) -- no requirement.
#   * Eigen is a git submodule (cpp/third_party/eigen), absent from a source
#     export, so `with_eigen` defaults False and the MED transpose uses the
#     hand-written fallback loop. Enabling it would need a small CMake change to
#     accept an external Eigen3::Eigen (tracked as a follow-up).
#   * HDF5/netCDF/zlib are PRIVATE to the shared library (the installed surface
#     is the C header alone), but they must be present at build and run time.

import os

from conan import ConanFile
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain, cmake_layout
from conan.tools.files import copy


class MeshioplusplusConan(ConanFile):
    name = "meshioplusplus"
    version = "6.6.2"
    license = "MIT"
    description = "C++ core for the meshio++ mesh I/O library (installable C API)"
    homepage = "https://github.com/loumalouomega/meshioplusplus"
    url = "https://github.com/loumalouomega/meshioplusplus"
    topics = ("mesh", "fem", "file-formats", "scientific-computing", "hpc")

    settings = "os", "compiler", "build_type", "arch"
    # The C API shared library is hardcoded SHARED in CMake today, so there is
    # no `shared` option (v1). A static build is a documented follow-up.
    options = {
        "fPIC": [True, False],
        "with_hdf5": [True, False],
        "with_netcdf": [True, False],
        "with_zlib": [True, False],
        "with_eigen": [True, False],
        "fortran": [True, False],
    }
    default_options = {
        "fPIC": True,
        "with_hdf5": True,
        "with_netcdf": True,
        "with_zlib": True,
        "with_eigen": False,  # submodule not in a source export -> fallback transpose
        "fortran": False,
    }

    # Everything the C API build needs. pugixml is inside cpp/third_party; the
    # Eigen submodule is deliberately not exported (with_eigen defaults off).
    exports_sources = (
        "CMakeLists.txt",
        "LICENSE",
        "cpp/include/*",
        "cpp/src/*",
        "cpp/third_party/pugixml/*",
        "bindings_c/*",
        "bindings_fortran/*",
        "cmake/*",
    )

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

    def configure(self):
        # Header-only / C consumer surface: no C++ settings leak from deps.
        self.settings.rm_safe("compiler.cppstd")

    def layout(self):
        cmake_layout(self)

    def requirements(self):
        if self.options.with_zlib:
            self.requires("zlib/[>=1.2.11 <2]")
        if self.options.with_netcdf:
            self.requires("netcdf/[>=4.8 <5]")
        if self.options.with_hdf5:
            # netcdf transitively pins an exact, older hdf5 (netcdf/4.8.1 ->
            # hdf5/1.14.3) than the top of our floating range, which Conan
            # reports as a hard version conflict. When both are on, pin our
            # direct hdf5 to netcdf's version so the graph resolves and we still
            # link hdf5 ourselves; otherwise float within the compatible range.
            # (The `packages` CI now runs on every recipe change, so a future
            # netcdf bumping its hdf5 pin is caught here rather than downstream.)
            if self.options.with_netcdf:
                self.requires("hdf5/1.14.3")
            else:
                self.requires("hdf5/[>=1.14 <2]")

    def generate(self):
        deps = CMakeDeps(self)
        deps.generate()
        tc = CMakeToolchain(self)
        tc.cache_variables["MESHIOPLUSPLUS_BUILD_C_API"] = True
        tc.cache_variables["MESHIOPLUSPLUS_BUILD_PYTHON"] = False
        tc.cache_variables["MESHIOPLUSPLUS_BUILD_FORTRAN"] = bool(self.options.fortran)
        tc.cache_variables["MESHIOPLUSPLUS_WITH_HDF5"] = bool(self.options.with_hdf5)
        tc.cache_variables["MESHIOPLUSPLUS_WITH_NETCDF"] = bool(
            self.options.with_netcdf
        )
        tc.cache_variables["MESHIOPLUSPLUS_WITH_ZLIB"] = bool(self.options.with_zlib)
        tc.cache_variables["MESHIOPLUSPLUS_WITH_EIGEN"] = bool(self.options.with_eigen)
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()
        copy(
            self,
            "LICENSE",
            src=self.source_folder,
            dst=os.path.join(self.package_folder, "licenses"),
        )

    def package_info(self):
        self.cpp_info.libs = ["meshioplusplus"]
        # Match the installed find_package config so downstream
        # find_package(meshioplusplus) + meshioplusplus::meshioplusplus resolve.
        self.cpp_info.set_property("cmake_file_name", "meshioplusplus")
        self.cpp_info.set_property(
            "cmake_target_name", "meshioplusplus::meshioplusplus"
        )
        self.cpp_info.set_property("pkg_config_name", "meshioplusplus")
        if self.options.fortran:
            self.cpp_info.libs.insert(0, "meshioplusplus_fortran")
