#include <catch2/catch_test_macros.hpp>
#include <Kokkos_Core.hpp>

// file for random Kokkos experiments

// run with "./build/tests/unit_tests [experiments]"
TEST_CASE("Random view experiments", "[.experiments]") {
    int numNodes = 4;
    Kokkos::View<int*> globalElm1("global_elm1_test", numNodes);
    Kokkos::View<double*> xCoordsElm1("xcoord_elm1_test", numNodes);
    Kokkos::View<double*> yCoordsElm1("ycoord_elm1_test", numNodes);

    auto globalElm1_host = Kokkos::create_mirror_view(globalElm1);
    auto xCoordsElm1_host = Kokkos::create_mirror_view(xCoordsElm1);
    auto yCoordsElm1_host = Kokkos::create_mirror_view(yCoordsElm1);

    globalElm1_host(3) = 5; globalElm1_host(2) = 4;
    globalElm1_host(0) = 0; globalElm1_host(1) = 1;

    xCoordsElm1_host(3) = 0; yCoordsElm1_host(3) = 2; xCoordsElm1_host(2) = 2; yCoordsElm1_host(2) = 2; 
    xCoordsElm1_host(0) = 0; yCoordsElm1_host(0) = 0; xCoordsElm1_host(1) = 2; yCoordsElm1_host(1) = 0; 

    Kokkos::deep_copy(globalElm1, globalElm1_host);
    Kokkos::deep_copy(xCoordsElm1, xCoordsElm1_host);
    Kokkos::deep_copy(yCoordsElm1, yCoordsElm1_host);

    /* no
    SECTION("can I put a subview into a view") {
        Kokkos::View<double[2]> subviewView = Kokkos::subview(globalElm1_host, std::pair<int,int>(0,2));
    }
    */

    SECTION("const subview?") {
        const auto subviewView = Kokkos::subview(globalElm1_host, std::pair<int,int>(0,2));
    }
}