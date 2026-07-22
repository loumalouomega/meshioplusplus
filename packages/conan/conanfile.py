# Conan test package (packages/conan/): builds a tiny C consumer against the packaged C API and
# runs it, proving the config-package + target name (meshioplusplus::meshioplusplus)
# resolve for a downstream find_package.
import os

from conan import ConanFile
from conan.tools.build import can_run
from conan.tools.cmake import CMake, cmake_layout


class MeshioplusplusTestConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeToolchain", "CMakeDeps", "VirtualRunEnv"
    test_type = "explicit"

    def requirements(self):
        self.requires(self.tested_reference_str)

    def layout(self):
        cmake_layout(self)

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def test(self):
        if can_run(self):
            self.run(
                os.path.join(self.cpp.build.bindir, "test_consumer"), env="conanrun"
            )
