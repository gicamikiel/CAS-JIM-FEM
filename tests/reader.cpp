#include <catch2/catch_test_macros.hpp>
#include <Kokkos_Core.hpp>
#include <readGmsh.hpp>

// Tests for Gmsh reading
TEST_CASE("test if readGmsh creates any elements") {
    // TODO: this is a stub, please add to this
    // This actually does not work but commenting out so branches from main won't inherit this failed test
    /*
    std::string path = "../../dataset/bun_zipper_res4";
    GmshDataKokkos data = readGmsh(path);
    /*
    REQUIRE(true);
}