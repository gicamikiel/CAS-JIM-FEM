#include <element.hpp>

int FiniteElement::getNumNodes() const { return 0; }
double FiniteElement::getXi(int idx) const { return 0; };
double FiniteElement::getEta(int idx) const { return 0; };
double FiniteElement::shape(int idx, double xi, double eta) const { return 0; };
double FiniteElement::shapeD(int idx, double xi, double eta, CoordinateParametric dim) const { return 0; };