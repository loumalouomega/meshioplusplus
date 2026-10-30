# The meshioplusplus-occt add-on's extension (MESHIOPLUSPLUS_ADDON=occt).
#
# Built by packages/python/meshioplusplus-occt against the installed base
# wheel: the core is found through the CMake package the wheel ships
# (meshioplusplus.get_cmake_dir()), never compiled again. The OpenCASCADE
# provider sources are the same ones a monolithic build compiles into the core.
# See doc/python_packages_plan.md.
cmake_minimum_required(VERSION 3.15)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

if(NOT meshioplusplus_DIR)
  execute_process(COMMAND "${Python_EXECUTABLE}" -c
                  "import meshioplusplus; print(meshioplusplus.get_cmake_dir())"
                  OUTPUT_VARIABLE meshioplusplus_DIR OUTPUT_STRIP_TRAILING_WHITESPACE
                  RESULT_VARIABLE _mio_rc)
  if(NOT _mio_rc EQUAL 0)
    message(FATAL_ERROR "meshio++: the base meshioplusplus wheel is not installed in the "
                        "build environment (pip install meshioplusplus==${PROJECT_VERSION})")
  endif()
endif()
find_package(meshioplusplus ${PROJECT_VERSION} EXACT CONFIG REQUIRED COMPONENTS CXX)
find_package(OpenCASCADE CONFIG REQUIRED)
if(OpenCASCADE_VERSION VERSION_LESS 7.6)
  message(FATAL_ERROR "meshio++: OpenCASCADE ${OpenCASCADE_VERSION} is older than 7.6")
endif()

# The toolkits, as the core's own probe picks them: 7.8 renamed the STEP and
# IGES ones into TKDESTEP and TKDEIGES.
set(_mio_occt_libs TKernel TKMath TKG2d TKG3d TKBRep TKGeomBase TKGeomAlgo TKTopAlgo
                   TKShHealing TKBO TKBool TKMesh TKPrim)
if(TARGET TKDESTEP OR OpenCASCADE_VERSION VERSION_GREATER_EQUAL 7.8)
  list(APPEND _mio_occt_libs TKDESTEP TKXSBase TKDEIGES)
else()
  list(APPEND _mio_occt_libs TKSTEP TKSTEPBase TKSTEPAttr TKSTEP209 TKXSBase TKIGES)
endif()

file(GLOB _mio_occt_sources src/cpp/providers/occt/*.cpp src/cpp/providers/occt/addon/*.cpp)
pybind11_add_module(_occt ${_mio_occt_sources})
target_compile_definitions(_occt PRIVATE MESHIOPLUSPLUS_HAS_OCCT)
target_include_directories(_occt PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/src)
target_include_directories(_occt SYSTEM PRIVATE ${OpenCASCADE_INCLUDE_DIR})
target_link_libraries(_occt PRIVATE meshioplusplus::core_meshio ${_mio_occt_libs})
set_target_properties(_occt PROPERTIES
  CXX_STANDARD 20
  INSTALL_RPATH "$<IF:$<BOOL:${APPLE}>,@loader_path/../meshioplusplus/lib;@loader_path,$ORIGIN/../meshioplusplus/lib:$ORIGIN>")

# manylinux's libstdc++ is older than the compiler's: linked statically, the
# runtime left exported, as in the core (CMakeLists.txt, MESHIOPLUSPLUS_PYTHON_SHARED_CORE).
option(MESHIOPLUSPLUS_STATIC_RUNTIME "Link libgcc and libstdc++ statically (GNU)" OFF)
if(MESHIOPLUSPLUS_STATIC_RUNTIME AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_options(_occt PRIVATE -static-libgcc -static-libstdc++)
endif()

install(TARGETS _occt DESTINATION meshioplusplus_occt)
