#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <functional>
#include <cmath>

#include <elementsFirstOrder.hpp>

using namespace Catch::Matchers;

// FirstOrderQuad

TEST_CASE("N_n is 1 at node n and all other N_i are 0 (quads)") {
    FirstOrderQuad quad;

    int node = GENERATE(0, 1, 2, 3);
    int function = GENERATE(0, 1, 2, 3);
    double xi = quad.getXi(node);
    double eta = quad.getEta(node);

    INFO("Node: " << node << " (xi, eta: = " << xi << "," << eta << ")");
    INFO("Function Subscript: " << function);
    REQUIRE_THAT(node == function ? 1.0 : 0.0, WithinRel(quad.shape(function, xi, eta)));
}

TEST_CASE("check N_n derivatives consistent with N_n using finite difference approximations (quads)") {
    FirstOrderQuad quad;

    int function = GENERATE(0, 1, 2, 3);
    double xi0 = GENERATE(-1/std::sqrt(3), 1/std::sqrt(3));
    double eta0 = GENERATE(-1/std::sqrt(3), 1/std::sqrt(3));
    double dx = 0.1;

    double analyticXi = quad.shapeD(function, xi0, eta0, XI);
    double numericXi = quad.shape(function, xi0+dx, eta0)-quad.shape(function, xi0-dx, eta0);
    numericXi /= 2*dx;

    double analyticEta = quad.shapeD(function, xi0, eta0, ETA);
    double numericEta = quad.shape(function, xi0, eta0+dx)-quad.shape(function, xi0, eta0-dx);
    numericEta /=2*dx;

    INFO("Shape Function: " << function);
    REQUIRE_THAT(analyticXi, WithinRel(numericXi));
    REQUIRE_THAT(numericEta, WithinRel(analyticEta));
}

// FirstOrderTri

TEST_CASE("N_n is 1 at node n and all other N_i are 0 (tri)") {
    FirstOrderTri tri;

    int node = GENERATE(0, 1, 2);
    int function = GENERATE(0, 1, 2);
    double xi = tri.getXi(node);
    double eta = tri.getEta(node);

    INFO("Node: " << node << " (xi, eta: = " << xi << "," << eta << ")");
    INFO("Function Subscript: " << function);
    REQUIRE_THAT(node == function ? 1.0 : 0.0, WithinRel(tri.shape(function, xi, eta)));
}

TEST_CASE("check N_n derivatives consistent with N_n using finite difference approximations (tri)") {
    FirstOrderTri tri;

    int function = GENERATE(0, 1, 2);
    double xi0 = GENERATE(-1/std::sqrt(3), 1/std::sqrt(3));
    double eta0 = GENERATE(-1/std::sqrt(3), 1/std::sqrt(3));
    double dx = 0.1;

    double analyticXi = tri.shapeD(function, xi0, eta0, XI);
    double numericXi = tri.shape(function, xi0+dx, eta0)-tri.shape(function, xi0-dx, eta0);
    numericXi /= 2*dx;

    double analyticEta = tri.shapeD(function, xi0, eta0, ETA);
    double numericEta = tri.shape(function, xi0, eta0+dx)-tri.shape(function, xi0, eta0-dx);
    numericEta /=2*dx;

    INFO("Shape Function: " << function);
    REQUIRE_THAT(analyticXi, WithinRel(numericXi));
    REQUIRE_THAT(numericEta, WithinRel(analyticEta));
}