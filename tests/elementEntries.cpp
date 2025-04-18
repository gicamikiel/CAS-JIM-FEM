#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Kokkos_Core.hpp>
#include <cmath>

#include <elementsFirstOrder.hpp>
#include <quadratureRules.hpp>
#include <readGmsh.hpp>

using namespace Catch::Matchers;

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

    FiniteElement<FirstOrderQuad, TwoPointGaussLegendre> quad;

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
}

TEST_CASE("Test local matrix flattening") {
    FiniteElement<FirstOrderQuad, TwoPointGaussLegendre> quad;

    SECTION("Check upper triangular index") {
        std::size_t test[2];
        quad.localStiffnessIdx(0, test);

        REQUIRE(quad.uniqueStiffnessEntriesN() == 10);
        REQUIRE(quad.integrationGridN() == 4);
        /*
        0 1 2 3
        - 4 5 6
        - - 7 8
        - - - 9
        */
        REQUIRE(test[0] == 0); REQUIRE(test[1] == 0);
        quad.localStiffnessIdx(1, test);
        REQUIRE(test[0] == 0); REQUIRE(test[1] == 1);
        quad.localStiffnessIdx(2, test);
        REQUIRE(test[0] == 0); REQUIRE(test[1] == 2);
        quad.localStiffnessIdx(3, test);
        REQUIRE(test[0] == 0); REQUIRE(test[1] == 3);
        quad.localStiffnessIdx(4, test);
        REQUIRE(test[0] == 1); REQUIRE(test[1] == 1);
        quad.localStiffnessIdx(5, test);
        REQUIRE(test[0] == 1); REQUIRE(test[1] == 2);
        quad.localStiffnessIdx(6, test);
        REQUIRE(test[0] == 1); REQUIRE(test[1] == 3);
        quad.localStiffnessIdx(7, test);
        REQUIRE(test[0] == 2); REQUIRE(test[1] == 2);
        quad.localStiffnessIdx(8, test);
        REQUIRE(test[0] == 2); REQUIRE(test[1] == 3);
        quad.localStiffnessIdx(9, test);
        REQUIRE(test[0] == 3); REQUIRE(test[1] == 3);
    }

    SECTION("Check unwrapping") {
        Kokkos::View<std::size_t*> global_tags("globalTags_test", 4);
        auto host_globalTags = Kokkos::create_mirror_view(global_tags);
        host_globalTags(0) = 1; 
        host_globalTags(1) = 2;
        host_globalTags(2) = 3;
        host_globalTags(3) = 4; 
        for(std::size_t i = 0; i < quad.uniqueStiffnessEntriesN(); i++) {
            double evalPt[2];
            std::size_t localIdx[2]; 
            quad.unwrapStiffnessN(i, evalPt, localIdx);
            CAPTURE(i, evalPt[0], evalPt[1], localIdx[0], localIdx[1]);
        }
    }
}