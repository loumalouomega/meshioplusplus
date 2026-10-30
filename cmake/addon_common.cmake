# The preamble every add-on wheel's build shares (cmake/addon_<name>.cmake):
# Python and pybind11, and the installed base wheel's CMake package, which the
# extension links instead of compiling the core again. See
# doc/python_packages_plan.md.
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

# manylinux's libstdc++ is older than the compiler's: linked statically, the
# runtime left exported, as in the core (CMakeLists.txt, MESHIOPLUSPLUS_PYTHON_SHARED_CORE).
option(MESHIOPLUSPLUS_STATIC_RUNTIME "Link libgcc and libstdc++ statically (GNU)" OFF)

# Finishes an add-on's extension target: C++20, an RPATH to the base wheel's
# core beside the add-on's own libraries, the static runtime, the install rule.
function(meshioplusplus_finish_addon target package)
  set_target_properties(${target} PROPERTIES
    CXX_STANDARD 20
    INSTALL_RPATH "$<IF:$<BOOL:${APPLE}>,@loader_path/../meshioplusplus/lib;@loader_path,$ORIGIN/../meshioplusplus/lib:$ORIGIN>")
  if(MESHIOPLUSPLUS_STATIC_RUNTIME AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_options(${target} PRIVATE -static-libgcc -static-libstdc++)
  endif()
  install(TARGETS ${target} DESTINATION ${package})
endfunction()
