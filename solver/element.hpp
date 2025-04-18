#ifndef SOLVER_ELEMENT
#define SOLVER_ELEMENT

#include <Kokkos_Macros.hpp>
#include <Kokkos_Core.hpp>
#include <iostream>

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

        KOKKOS_FUNCTION
        double computeJacobianAt(double xi, double eta, Kokkos::View<double*> xCoords, Kokkos::View<double*> yCoords) const {
            double dx_dxi = 0;
            double dx_deta = 0;
            double dy_dxi = 0;
            double dy_deta = 0;
            for (std::size_t k=0; k<this->numNodes; k++) {
                // create del x/del xi
                dx_dxi += this->shapeD(k, xi, eta, XI)*xCoords(k);
                // create del x/del eta
                dx_deta += this->shapeD(k, xi, eta, ETA)*xCoords(k);
                // create del y/del xi
                dy_dxi += this->shapeD(k, xi, eta, XI)*yCoords(k);
                // create del y/del eta
                dy_deta += this->shapeD(k, xi, eta, ETA)*yCoords(k);
                printf("%.2f, %.2f \n", xCoords(k), yCoords(k));
            }
            return dx_dxi*dy_deta-dx_deta*dy_dxi;
        }

        double stiffnessIntegrand(std::size_t subA, std::size_t subB, double xi, double eta, Kokkos::View<double*> xCoords, Kokkos::View<double*> yCoords) const {
            double out = this->shapeD(subA, xi, eta, XI)*this->shapeD(subB, xi, eta, XI)
                            +this->shape(subA, xi, eta, ETA)*this->shapeD(subB, xi, eta, ETA);
            return out/computeJacobianAt(xi, eta, xCoords, yCoords);
        }

        std::size_t stiffnessEvalsN() const {
            return uniqueStiffnessEntriesN()*integrationGridN();
        }

        std::size_t uniqueStiffnessEntriesN() const {
            // local stiffness matrix is symmetric so we only need to compute some of it
            // numPoints + numPoints-1 + numPoints-2 + ... + 1
            return this->numNodes*(this->numNodes+1)/2;
        }

        std::size_t integrationGridN() const {
            return this->numPoints*this->numPoints;
        }

        void unwrapStiffnessN(std::size_t i, double parEvalPt[2], std::size_t localIdx[2]) const {
            std::size_t integrationIdx = i % integrationGridN();
            std::size_t uniqueEntryIdx = i / integrationGridN();

            parEvalPt[0] = this->points[integrationIdx % this->numPoints];
            parEvalPt[1] = this->points[integrationIdx / this->numPoints];
            
            localStiffnessIdx(uniqueEntryIdx, localIdx);
        }

        // based on https://stackoverflow.com/questions/242711/algorithm-for-index-numbers-of-triangular-matrix-coefficients
        void localStiffnessIdx(std::size_t i, std::size_t idx[2]) const {
            std::size_t row = 0;
            std::size_t delta = this->numNodes - 1;
            std::size_t x;
            for( x = delta; x < i; x += delta-- ){
                row++;
            }
            idx[0] = row;
            idx[1] = this->numNodes + i - x - 1;
        }
};

#endif