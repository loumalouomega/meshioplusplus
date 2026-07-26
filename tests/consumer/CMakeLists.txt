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

find_package(meshioplusplus CONFIG REQUIRED COMPONENTS CXX C)

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

add_executable(consumer_smoketest main.cpp)

# Both components in one binary, on purpose: proving they coexist is part of what
# this test is for.
target_link_libraries(consumer_smoketest
  PRIVATE ${_core_target} meshioplusplus::meshioplusplus)

enable_testing()
add_test(NAME consumer_smoketest COMMAND consumer_smoketest)
