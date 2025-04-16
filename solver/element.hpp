#ifndef SOLVER_ELEMENT
#define SOLVER_ELEMENT

#include <Kokkos_Macros.hpp>

enum CoordinateParametric {
    NONE = 0,
    XI = 1,
    ETA = 2
};

// Cannot use virtual functions for GPU code
// https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Kokkos-and-Virtual-Functions.html
// Try using CRTP
// https://en.cppreference.com/w/cpp/language/crtp

template <class T>
class FiniteElement {
    public:
        // Definitions must be here because of how templates work

        // INTERFACE FUNCTIONS
        
        int getNumNodes() const {
            return static_cast<T*>(this)->getNumNodes(); 
        };

        double getXi(int idx) const { 
            return static_cast<T*>(this)->getXi(idx); 
        };

        double getEta(int idx) const {
            return static_cast<T*>(this)->getEta(idx);
        };

        double shape(int idx, double xi, double eta) const {
            return static_cast<T*>(this)->shape(idx, xi, eta);
        };

        double shapeD(int idx, double xi, double eta, CoordinateParametric dim) const {
            return static_cast<T*>(this)->shape(idx, xi, eta, dim);
        };


        // MATH AND BOOK-KEEPING

        int getGlobalRow() const;

        int getGlobalCol() const;

        double computeMatrixEntry() const;

    protected:
        KOKKOS_INLINE_FUNCTION FiniteElement() = default;
        KOKKOS_INLINE_FUNCTION ~FiniteElement() = default;
};

#endif