#ifndef SOLVER_ELEMENT
#define SOLVER_ELEMENT

#include<Kokkos_Macros.hpp>

enum CoordinateParametric {
    NONE = 0,
    XI = 1,
    ETA = 2
};

// Use templates instead of virtual functions as vtables are ill-behaved in the GPU
// see https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Kokkos-and-Virtual-Functions.html

template <typename ElmType, typename IntRule>
class FiniteElementDef : public ElmType, public IntRule {
    public:
        // Definitions must be here because of how templates work

        // Functions for evaluating jacobian

        KOKKOS_FUNCTION
        double jacAddendDelxDelXi(std::size_t node, double xi, double eta, double xa) const {
            return this->shapeD(node, xi, eta, XI)*xa;
        }

        KOKKOS_FUNCTION
        double jacAddendDelyDelEta(std::size_t node, double xi, double eta, double ya) const {
            return this->shapeD(node, xi, eta, ETA)*ya;
        }

        KOKKOS_FUNCTION
        double jacAddendDelxDelEta(std::size_t node, double xi, double eta, double xa) const {
            return this->shapeD(node, xi, eta, ETA)*xa;
        }

        KOKKOS_FUNCTION
        double jacAddendDelyDelXi(std::size_t node, double xi, double eta, double ya) const {
            return this->shapeD(node, xi, eta, XI)*ya;
        }

        KOKKOS_FUNCTION
        double jacAddendEvalsN() const {
            return 4 * // 4 derivatives required to compute Jacobian
                this->numNodes * // each derivative requires an addend at each node
                this->numQuad; // need to evaluate Jacobian at quadrature points
        }

        KOKKOS_FUNCTION
        void unwrapAddendEvalsN(std::size_t i, std::size_t &node, std::size_t &quadPt, std::size_t &derivatve) const {
            node = i % this->numNodes;
            quadPt = (i / this->numNodes) % this->numQuad;
            derivatve = ((i / this->numNodes) / this->numQuad) % 4;
        }

        // Functions for evaluating stiffness matrix
        KOKKOS_FUNCTION
        double stiffnessIntegrand(std::size_t subA, std::size_t subB, double xi, double eta, double jac) const {
            double out = this->shapeD(subA, xi, eta, XI)*this->shapeD(subB, xi, eta, XI)
                            +this->shape(subA, xi, eta, ETA)*this->shapeD(subB, xi, eta, ETA);
            return out/jac;
        }

        KOKKOS_FUNCTION
        double stiffnessQuadEval(std::size_t subA, std::size_t subB, std::size_t quadPt, double jac) {
            return this->quadWt[quadPt] * stiffnessIntegrand(subA, subB, this->quadXi[quadPt], this->quadEta[quadPt], jac);
        }

        KOKKOS_FUNCTION
        std::size_t uniqueStiffnessEntriesN() const {
            // local stiffness matrix is symmetric so we only need to compute some of it
            // numPoints + numPoints-1 + numPoints-2 + ... + 1
            return this->numNodes*(this->numNodes+1)/2;
        }

        KOKKOS_FUNCTION
        std::size_t stiffnessEvalsN() const {
            return uniqueStiffnessEntriesN()*this->numQuad;
        }

        KOKKOS_FUNCTION
        void unwrapStiffnessN(std::size_t i, std::size_t &quadPt, std::size_t &localRow, std::size_t &localCol) const {
            std::size_t uniqueEntryIdx = i % uniqueStiffnessEntriesN();
            quadPt = i / uniqueStiffnessEntriesN();
            
            localStiffnessIdx(uniqueEntryIdx, localRow, localCol);
        }

        KOKKOS_FUNCTION
        // based on https://stackoverflow.com/a/242924
        void localStiffnessIdx(std::size_t i, std::size_t &localRow, std::size_t &localCol) const {
            localRow = 0;
            std::size_t delta = this->numNodes - 1;
            std::size_t x;
            for( x = delta; x < i; x += delta-- ){
                localRow++;
            }
            localCol = this->numNodes + i - x - 1;
        }
};

#endif