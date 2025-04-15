#ifndef SOLVER_ELEMENTS_FIRST_ORDER
#define SOLVER_ELEMENTS_FIRST_ORDER

#include <element.hpp>

class FirstOrderQuad : public FiniteElement {
    public:
        int getNumNodes() const override;
        double getXi(int idx) const override;
        double getEta(int idx) const override;
        double shape(int idx, double xi, double eta) const override;
        double shapeD(int idx, double xi, double eta, CoordinateParametric dim) const override;

    private:
        double parXi[4] = {-1, 1, 1, -1};
        double parEta[4] = {-1, -1, 1, 1};
};

class FirstOrderTri : public FiniteElement {
    public:
        int getNumNodes() const override;
        double getXi(int idx) const override;
        double getEta(int idx) const override;
        double shape(int idx, double xi, double eta) const override;
        double shapeD(int idx, double xi, double eta, CoordinateParametric dim) const override;

    private:
        double parXi[3] = {0, 1, 0};
        double parEta[3] = {0, 0, 1};
};

#endif