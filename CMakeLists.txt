cmake_minimum_required(VERSION 3.20.0)
project(cas_fem_jim CXX)

find_package(Catch2 3 REQUIRED CONFIG)
find_package(Kokkos 4.6 REQUIRED CONFIG)

add_subdirectory(math_utils)
add_subdirectory(solver)
add_subdirectory(tests)