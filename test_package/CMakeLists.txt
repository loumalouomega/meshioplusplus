cmake_minimum_required(VERSION 3.15)
project(test_package LANGUAGES C)

find_package(meshioplusplus CONFIG REQUIRED)

add_executable(test_consumer test_consumer.c)
target_link_libraries(test_consumer PRIVATE meshioplusplus::meshioplusplus)
