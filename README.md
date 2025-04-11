# Assignment 2 - FEM Solver

by Jay Gaiardelli, Ickbum Kim, Mikel Gica

## Build Instructions

To build on an appropriate CCI node, first run the following.
```
module load cmake
module load gcc
module load spectrum-mpi
module load cuda/11.2
```
If not on CCI, follow the appropriate steps to setup CMake and parallelism.

To configure and compile, run `autoconfig.sh`.

`autoconfig.sh` has parameters `autoconfig.sh <CMAKE_BUILD_TYPE> <CMAKE_CXX_COMPILER> <HOST_BACKEND> <DEVICE_BACKEND>`.

where `<CMAKE_BUILD_TYPE>` is the CMake build type (e.g. `Debug`); `<CMAKE_CXX_COMPILER>` is the C++ compiler (e.g. `g++`); `<HOST_BACKEND>` and `<DEVICE_BACKEND>` are [valid Kokkos backends](https://kokkos.org/kokkos-core-wiki/get-started/configuration-guide.html#keywords-backends) (e.g. to use `Kokkos_ENABLE_CUDA` to `ON` set `<DEVICE_BACKEND>` to `CUDA`).

for systems with no working CUDA support, run (for example):

`./autoconfig.sh DEBUG g++ THREADS NOCUDA`

To compile without configuring, run `cmake --build build -j$(nproc)`.

## Testing

To run tests manually, run the following from the root directory.

`ctest --test-dir ./build/tests`

**Note:** To test new features, please prefer to implement a Catch2 test instead of a main function.