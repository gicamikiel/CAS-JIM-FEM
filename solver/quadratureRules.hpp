#ifndef SOLVER_QUADRATURE_RULES
#define SOLVER_QUADRATURE_RULES

#include <cmath>

struct TwoByTwoGaussLegendre {
    std::size_t const numQuad = 4;
    double const quadXi[4] = {1.0/std::sqrt(3), -1.0/std::sqrt(3), 1.0/std::sqrt(3), 1.0/std::sqrt(3)};
    double const quadEta[4] = {1.0/std::sqrt(3), 1.0/std::sqrt(3), -1.0/std::sqrt(3), 1.0/std::sqrt(3)};
    double const quadWt[4] = {1.0, 1.0, 1.0, 1.0};
};

#endif