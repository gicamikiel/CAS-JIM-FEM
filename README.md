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

To compile without configuring, run `cmake --build build -j$(nproc)`.