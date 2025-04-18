#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Kokkos_Core.hpp>
#include <cmath>

#include <elementsFirstOrder.hpp>
#include <readGmsh.hpp>

using namespace Catch::Matchers;

typedef FiniteElementDef<FirstOrderQuad, TwoByTwoGaussLegendre> Quad;
typedef FiniteElementDef<FirstOrderTri, TwoByTwoGaussLegendre> Tri;

TEST_CASE("Test on handwritten mesh") {
    /*
    DOF indices
    5--4--3
    |  |  |
    0--1--2
    */
    int numNodes = 4;
    Kokkos::View<int*> globalElm1("global_elm1_test", numNodes);
    Kokkos::View<double*> xCoordsElm1("xcoord_elm1_test", numNodes);
    Kokkos::View<double*> yCoordsElm1("ycoord_elm1_test", numNodes);

    auto globalElm1_host = Kokkos::create_mirror_view(globalElm1);
    auto xCoordsElm1_host = Kokkos::create_mirror_view(xCoordsElm1);
    auto yCoordsElm1_host =Kokkos::create_mirror_view(yCoordsElm1);

    globalElm1_host(3) = 5; globalElm1_host(2) = 4;
    globalElm1_host(0) = 0; globalElm1_host(1) = 1;

    xCoordsElm1_host(3) = 0; yCoordsElm1_host(3) = 2; xCoordsElm1_host(2) = 2; yCoordsElm1_host(2) = 2; 
    xCoordsElm1_host(0) = 0; yCoordsElm1_host(0) = 0; xCoordsElm1_host(1) = 2; yCoordsElm1_host(1) = 0; 

    Kokkos::deep_copy(globalElm1, globalElm1_host);
    Kokkos::deep_copy(xCoordsElm1, xCoordsElm1_host);
    Kokkos::deep_copy(yCoordsElm1, yCoordsElm1_host);

    Kokkos::View<int*> globalElm2("global_elm2_test", numNodes);
    Kokkos::View<double*> xCoordsElm2("xcoord_elm2_test", numNodes);
    Kokkos::View<double*> yCoordsElm2("ycoord_elm2_test", numNodes);

    auto globalElm2_host = Kokkos::create_mirror_view(globalElm2);
    auto xCoordsElm2_host = Kokkos::create_mirror_view(xCoordsElm2);
    auto yCoordsElm2_host =Kokkos::create_mirror_view(yCoordsElm2);

    globalElm2_host(3) = 4; globalElm2_host(2) = 3;
    globalElm2_host(0) = 1; globalElm2_host(1) = 2;

    xCoordsElm2_host(3) = 2; yCoordsElm2_host(3) = 2; xCoordsElm2_host(2) = 4; yCoordsElm2_host(2) = 2; 
    xCoordsElm2_host(0) = 2; yCoordsElm2_host(0) = 0; xCoordsElm2_host(1) = 4; yCoordsElm2_host(1) = 0; 

    Kokkos::deep_copy(globalElm2, globalElm2_host);
    Kokkos::deep_copy(xCoordsElm2, xCoordsElm2_host);
    Kokkos::deep_copy(yCoordsElm2, yCoordsElm2_host);

    Quad quad;

    /*
    SECTION("Evaluate Jacobian") {
        double result = 0;
        double testxi = 0;
        double testeta = 0;
        Kokkos::parallel_reduce("test-jacobian", 1, KOKKOS_LAMBDA (const int i, double& out) {
            out = quad.computeJacobianAt(testxi, testeta, xCoordsElm1, yCoordsElm1);
        }, result);
        CAPTURE(quad.shapeD(0, testxi, testeta, XI), quad.shapeD(0, testxi, testeta, ETA));
        CAPTURE(quad.shapeD(1, testxi, testeta, XI), quad.shapeD(1, testxi, testeta, ETA));
        CAPTURE(quad.shapeD(2, testxi, testeta, XI), quad.shapeD(2, testxi, testeta, ETA));
        CAPTURE(quad.shapeD(3, testxi, testeta, XI), quad.shapeD(3, testxi, testeta, ETA));
        REQUIRE_THAT(1, WithinRel(result));
    }
    */
}

TEST_CASE("Check local matrix flattening") {
    Quad quad;

    SECTION("Check upper triangular index") {
        std::size_t row;
        std::size_t col;

        REQUIRE(quad.uniqueStiffnessEntriesN() == 10);
        /*
        0 1 2 3
        - 4 5 6
        - - 7 8
        - - - 9
        */
        quad.localStiffnessIdx(0, row, col); REQUIRE(row == 0); REQUIRE(col == 0);
        quad.localStiffnessIdx(1, row, col); REQUIRE(row == 0); REQUIRE(col == 1);
        quad.localStiffnessIdx(2, row, col); REQUIRE(row == 0); REQUIRE(col == 2);
        quad.localStiffnessIdx(3, row, col); REQUIRE(row == 0); REQUIRE(col == 3);
        quad.localStiffnessIdx(4, row, col); REQUIRE(row == 1); REQUIRE(col == 1);
        quad.localStiffnessIdx(5, row, col); REQUIRE(row == 1); REQUIRE(col == 2);
        quad.localStiffnessIdx(6, row, col); REQUIRE(row == 1); REQUIRE(col == 3);
        quad.localStiffnessIdx(7, row, col); REQUIRE(row == 2); REQUIRE(col == 2);
        quad.localStiffnessIdx(8, row, col); REQUIRE(row == 2); REQUIRE(col == 3);
        quad.localStiffnessIdx(9, row, col); REQUIRE(row == 3); REQUIRE(col == 3);
    }
}

// test that needs manual checking of printf output
// run with "./build/tests/unit_tests [manual]"
TEST_CASE("Check unwrapping", "[.manual]") {
    Quad quad;

    printf("Stiffness Unwrap\n");
    for(std::size_t i = 0; i < quad.stiffnessEvalsN(); i++) {
        std::size_t quadPoint;
        std::size_t row;
        std::size_t col;
        quad.unwrapStiffnessN(i, quadPoint, row, col);
        printf("%d: qp=%d r,c=(%d, %d) \n", (unsigned int)i, (unsigned int)quadPoint, (unsigned int)row, (unsigned int)col);
    }

    printf("Derivative Addend Unwrap\n");
    for(std::size_t i = 0; i < quad.stiffnessEvalsN(); i++) {
        std::size_t node;
        std::size_t quadPt;
        std::size_t derivatve;
        quad.unwrapAddendEvalsN(i, node, quadPt, derivatve);
        printf("%d: node=%d qp=%d derivative=%d \n", (unsigned int)i, (unsigned int)node, (unsigned int)quadPt, (unsigned int)derivatve);
    }
}