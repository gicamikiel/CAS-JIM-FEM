struct TwoPointGaussLegendre {
    std::size_t const numPoints = 2;
    double const points[2] = {-1.0/std::sqrt(3), 1.0/std::sqrt(3)};
    double const weight[2] = {1.0, 1.0};
};