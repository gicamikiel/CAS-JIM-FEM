#ifndef SOLVER_ELEMENTS_FIRST_ORDER
#define SOLVER_ELEMENTS_FIRST_ORDER

#include <element.hpp>
#include <quadratureRules.hpp>

struct FirstOrderQuad {
    std::size_t const numNodes = 4;
    double const parXi[4] = {-1, 1, 1, -1};
    double const parEta[4] = {-1, -1, 1, 1};

    KOKKOS_FUNCTION double shape(std::size_t node, double xi, double eta) const;
    KOKKOS_FUNCTION double shapeD(std::size_t node, double xi, double eta, CoordinateParametric dim) const;
};

struct FirstOrderTri {
    std::size_t const numNodes = 4;
    double const parXi[3] = {0, 1, 0};
    double const parEta[3] = {0, 0, 1};

    KOKKOS_FUNCTION double shape(std::size_t node, double xi, double eta) const;
    KOKKOS_FUNCTION double shapeD(std::size_t node, double xi, double eta, CoordinateParametric dim) const;
};

#endif