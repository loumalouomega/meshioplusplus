cmake_minimum_required(VERSION 3.15...3.30)

project(
  meshio_core
  LANGUAGES CXX
  DESCRIPTION "C++ core for the meshio mesh I/O library")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# scikit-build-core provides the right Python; locate the module-only component.
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)
find_package(ZLIB REQUIRED)

# C++ core sources (grows as formats are ported).
file(GLOB_RECURSE MESHIO_CORE_SOURCES
    bindings/*.cpp
    cpp/src/*.cpp
    cpp/third_party/pugixml/pugixml.cpp)

pybind11_add_module(_core ${MESHIO_CORE_SOURCES})
target_include_directories(
  _core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include
                ${CMAKE_CURRENT_SOURCE_DIR}/bindings
                ${CMAKE_CURRENT_SOURCE_DIR}/cpp/third_party/pugixml)
target_compile_features(_core PRIVATE cxx_std_17)
target_link_libraries(_core PRIVATE ZLIB::ZLIB)

# Land the extension next to the pure-Python package inside the wheel.
install(TARGETS _core DESTINATION meshio)
