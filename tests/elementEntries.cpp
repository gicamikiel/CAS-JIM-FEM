#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <Kokkos_Core.hpp>
#include <cmath>

#include <elementsFirstOrder.hpp>
#include <quadratureRules.hpp>
#include <readGmsh.hpp>

TEST_CASE("Test Jacobian computation") {
    int numNodes = 4;
    Kokkos::View<std::size_t*, Kokkos::HostSpace> host_nodeTags("nodeTags_test", numNodes);
    Kokkos::View<double**, Kokkos::HostSpace> host_nodeCoords("nodeCoords_test", numNodes, 3);

    host_nodeTags(0) = 1; host_nodeTags(1) = 2; host_nodeTags(2) = 3;
    host_nodeCoords(0, 0) = -2.0; host_nodeCoords(1, 0) =  2.0; host_nodeCoords(2, 0) = 2.0; host_nodeCoords(3, 0) = -2.0;
    host_nodeCoords(0, 1) = -2.0; host_nodeCoords(1, 1) = -2.0; host_nodeCoords(2, 1) = 2.0; host_nodeCoords(3, 1) =  2.0;

    FiniteElement<FirstOrderQuad, TwoPointGaussLegendre> quad;
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