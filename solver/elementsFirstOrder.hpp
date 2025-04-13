#ifndef SOLVER_ELEMENTS_FIRST_ORDER
#define SOLVER_ELEMENTS_FIRST_ORDER

#include <element.hpp>

enum CoordinateParametric {
    NONE = 0,
    XI = 1,
    ETA = 2
};

class FirstOrderQuad : public FiniteElement {
    public:
        int getNumNodes() const override;

        double getXi(int idx) const;
        double getEta(int idx) const;
        double shape(int idx, double xi, double eta) const;
        double shapeD(int idx, double xi, double eta, CoordinateParametric dim) const;

    private:
        double parXi[4] = {-1, 1, 1, -1};
        double parEta[4] = {-1, -1, 1, 1};
};

#endif