#ifndef SOLVER_ELEMENT
#define SOLVER_ELEMENT

#include <Kokkos_Macros.hpp>

enum CoordinateParametric {
    NONE = 0,
    XI = 1,
    ETA = 2
};

// Use templates instead of virtual functions as vtables are ill-behaved in the GPU
// see https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Kokkos-and-Virtual-Functions.html

template <typename ElmType, typename IntRule>
class FiniteElement : public ElmType, public IntRule {
    public:
        // Definitions must be here because of how templates work

        KOKKOS_FUNCTION int getGlobalRow(int node) const {
            return 0;
        }

        KOKKOS_FUNCTION int getGlobalCol(int node) const {
            return 0;
        }

        KOKKOS_FUNCTION double computeMatrixEntry(int node) const {
            return 0;
        }
};

#endif