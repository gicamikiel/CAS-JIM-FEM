#ifndef SOLVER_ELEMENTS_FIRST_ORDER
#define SOLVER_ELEMENTS_FIRST_ORDER

#include <element.hpp>

class FirstOrderQuad : public FiniteElement<FirstOrderQuad> {
    public:
        int getNumNodes() const;
        double getXi(int idx) const;
        double getEta(int idx) const;
        double shape(int idx, double xi, double eta) const;
        double shapeD(int idx, double xi, double eta, CoordinateParametric dim) const;

    private:
        double const parXi[4] = {-1, 1, 1, -1};
        double const parEta[4] = {-1, -1, 1, 1};
};

class FirstOrderTri : public FiniteElement<FirstOrderTri> {
    public:
        int getNumNodes() const;
        double getXi(int idx) const;
        double getEta(int idx) const;
        double shape(int idx, double xi, double eta) const;
        double shapeD(int idx, double xi, double eta, CoordinateParametric dim) const;

    private:
        double const parXi[3] = {0, 1, 0};
        double const parEta[3] = {0, 0, 1};
};

#endif