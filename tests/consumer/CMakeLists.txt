# Standalone consumer of an INSTALLED meshio++. Deliberately a separate project:
# it is never add_subdirectory()'d from the main tree, and it must build with the
# meshio++ source tree unreachable -- everything it needs comes from the install
# prefix via find_package().
#
#   cmake -S tests/consumer -B build-consumer \
#         -Dmeshioplusplus_DIR=<prefix>/lib/cmake/meshioplusplus
#   cmake --build build-consumer && ./build-consumer/consumer_smoketest
#
# MESHIOPLUSPLUS_CONSUMER_BACKEND picks which backend variant to link
# (meshio|native|kratos); unset means meshioplusplus::core, i.e. whichever
# backend the install treats as its default.
cmake_minimum_required(VERSION 3.15)
project(meshioplusplus_consumer_smoketest LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# doc/cpp_api.md pins COMPONENTS CXX with an EXACT version, and CMake's
# SameMajorVersion EXACT mode is a full-STRING compare -- "9.2" never matches
# "9.2.0", so the documented line has to spell all three components. Through
# v9.1.0 the docs printed the two-component form, which could never succeed and
# which nothing here exercised. CI now passes the real version (from
# `pkg-config --modversion` against the install itself, so nothing is hard-coded
# and nothing needs bumping at release time) to make that line a gate.
# Empty -- the default -- keeps the historical unversioned call, so every
# existing invocation is unchanged.
set(MESHIOPLUSPLUS_CONSUMER_EXACT_VERSION "" CACHE STRING
    "Request this exact version from find_package (empty: no version request)")

if(MESHIOPLUSPLUS_CONSUMER_EXACT_VERSION)
  find_package(meshioplusplus ${MESHIOPLUSPLUS_CONSUMER_EXACT_VERSION} EXACT
               CONFIG REQUIRED COMPONENTS CXX C)
else()
  find_package(meshioplusplus CONFIG REQUIRED COMPONENTS CXX C)
endif()

set(MESHIOPLUSPLUS_CONSUMER_BACKEND "" CACHE STRING
    "Backend variant to link: meshio, native, kratos, or empty for meshioplusplus::core")

if(MESHIOPLUSPLUS_CONSUMER_BACKEND)
  set(_core_target meshioplusplus::core_${MESHIOPLUSPLUS_CONSUMER_BACKEND})
  if(NOT TARGET ${_core_target})
    message(FATAL_ERROR
      "meshio++ consumer: ${_core_target} is not in this install "
      "(it shipped: ${MESHIOPLUSPLUS_MESH_BACKENDS})")
  endif()
else()
  set(_core_target meshioplusplus::core)
endif()

message(STATUS "meshio++ consumer: linking ${_core_target}")
message(STATUS "meshio++ consumer: default backend  ${MESHIOPLUSPLUS_MESH_BACKEND}")
message(STATUS "meshio++ consumer: parallel backend ${MESHIOPLUSPLUS_PARALLEL_BACKEND}")
message(STATUS "meshio++ consumer: backends shipped ${MESHIOPLUSPLUS_MESH_BACKENDS}")

# MESHIOPLUSPLUS_NO_STD_SPAN guards only the <span> include and the inline
# NativeMesh::ConnSpan() accessor -- no member variable, so it is ABI-neutral --
# and the guard is #ifndef. A consumer whose own Boost uBLAS collides with
# <span> under MSVC must therefore be able to define it themselves even against
# a prefix that did not. This option is how CI proves that.
option(MESHIOPLUSPLUS_CONSUMER_NO_STD_SPAN
       "Define MESHIOPLUSPLUS_NO_STD_SPAN on the consumer side only" OFF)

add_executable(consumer_smoketest main.cpp)

if(MESHIOPLUSPLUS_CONSUMER_NO_STD_SPAN)
  target_compile_definitions(consumer_smoketest PRIVATE MESHIOPLUSPLUS_NO_STD_SPAN)
endif()

# Both components in one binary, on purpose: proving they coexist is part of what
# this test is for.
target_link_libraries(consumer_smoketest
  PRIVATE ${_core_target} meshioplusplus::meshioplusplus)

enable_testing()
add_test(NAME consumer_smoketest COMMAND consumer_smoketest)
