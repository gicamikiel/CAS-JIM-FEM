#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <Kokkos_Core.hpp>
#include <cmath>

#include <elementsFirstOrder.hpp>

using namespace Catch::Matchers;

typedef FiniteElementDef<FirstOrderQuad, TwoByTwoGaussLegendre> Quad;
typedef FiniteElementDef<FirstOrderTri, ThreePointTriangle> Tri;

#ifdef KOKKOS_ENABLE_CUDA
TEST_CASE("Test if class members are callable from the GPU") {
    Quad quad;
    Tri tri;

    int dummy = 0;
    Kokkos::parallel_reduce("test-quad-parallelizable", 1, KOKKOS_LAMBDA (const int i, int& count) {
        quad.shape(0, 0, 0);
        tri.shape(0, 0, 0);
    }, dummy);
}
#endif

// FirstOrderQuad

TEST_CASE("N_n is 1 at node n and all other N_i are 0 (quads)") {
    Quad quad;

    int node = GENERATE(0, 1, 2, 3);
    int function = GENERATE(0, 1, 2, 3);
    double xi = quad.parXi[node];
    double eta = quad.parEta[node];

    INFO("Node: " << node << " (xi, eta: = " << xi << "," << eta << ")");
    INFO("Function Subscript: " << function);
    REQUIRE_THAT(node == function ? 1.0 : 0.0, WithinRel(quad.shape(function, xi, eta)));
}

TEST_CASE("check N_n derivatives consistent with N_n using finite difference approximations (quads)") {
    Quad quad;

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

TEST_CASE("check all N_n adds up to 1 (quads)") {
    Quad quad;

    double xi0 = GENERATE(-1/std::sqrt(3), 1/std::sqrt(3));
    double eta0 = GENERATE(-1/std::sqrt(3), 1/std::sqrt(3));

    double sum = 0;
    for(int i=0; i<quad.numNodes; i++) {
        sum += quad.shape(i, xi0, eta0);
    }

    INFO("Point: " << xi0 << ", " << eta0);
    REQUIRE_THAT(1.0, WithinRel(sum));
}

// FirstOrderTri

TEST_CASE("N_n is 1 at node n and all other N_i are 0 (tri)") {
    Tri tri;

    int node = GENERATE(0, 1, 2);
    int function = GENERATE(0, 1, 2);
    double xi = tri.parXi[node];
    double eta = tri.parEta[node];

    INFO("Node: " << node << " (xi, eta: = " << xi << "," << eta << ")");
    INFO("Function Subscript: " << function);
    REQUIRE_THAT(node == function ? 1.0 : 0.0, WithinRel(tri.shape(function, xi, eta)));
}

TEST_CASE("check N_n derivatives consistent with N_n using finite difference approximations (tri)") {
    Tri tri;

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

TEST_CASE("check all N_n adds up to 1 (tris)") {
    Tri tris;

    double xi0 = GENERATE(-1/std::sqrt(3), 1/std::sqrt(3));
    double eta0 = GENERATE(-1/std::sqrt(3), 1/std::sqrt(3));

    double sum = 0;
    for(int i=0; i<tris.numNodes; i++) {
        sum += tris.shape(i, xi0, eta0);
    }

    INFO("Point: " << xi0 << ", " << eta0);
    REQUIRE_THAT(1.0, WithinRel(sum));
}