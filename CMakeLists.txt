cmake_minimum_required(VERSION 3.22.0)
project(cas_fem_jim CXX)

set(Kokkos_DIR "$ENV{HOME}/kokkos-install/lib/cmake/Kokkos")
find_package(Kokkos REQUIRED)

find_package(Catch2 3 REQUIRED)



add_subdirectory(math_utils)
add_subdirectory(solver)
add_subdirectory(tests)