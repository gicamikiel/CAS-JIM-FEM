#ifndef SOLVER_QUADRATURE_RULES
#define SOLVER_QUADRATURE_RULES

#include <cmath>

struct TwoByTwoGaussLegendre {
    std::size_t const numQuad = 4;
    //double const quadXi[4] = {1.0/std::sqrt(3), -1.0/std::sqrt(3), 1.0/std::sqrt(3), -1.0/std::sqrt(3)};
    //double const quadEta[4] = {1.0/std::sqrt(3), 1.0/std::sqrt(3), -1.0/std::sqrt(3), -1.0/std::sqrt(3)};
    #define XI1 0.577350269189626
    double const quadXi[4] = {XI1, -XI1, XI1, -XI1};
    double const quadEta[4] = {XI1, XI1, -XI1, -XI1};
    double const quadWt[4] = {1.0, 1.0, 1.0, 1.0};
};

struct ThreePointTriangle {
    std::size_t const numQuad = 3;
    //double const quadXi[4] = {1.0/std::sqrt(3), -1.0/std::sqrt(3), 1.0/std::sqrt(3), -1.0/std::sqrt(3)};
    //double const quadEta[4] = {1.0/std::sqrt(3), 1.0/std::sqrt(3), -1.0/std::sqrt(3), -1.0/std::sqrt(3)};
    #define L1 .666666666666666
    #define L2 .166666666666666
    #define WT 0.333333333333333
    double const quadXi[4] = {L1, L2, L2};
    double const quadEta[4] = {L2, L1, L2};
    double const quadWt[4] = {WT, WT, WT};
};

#endif