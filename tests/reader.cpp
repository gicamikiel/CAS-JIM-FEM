#include <catch2/catch_test_macros.hpp>
#include <Kokkos_Core.hpp>
#include <readGmsh.hpp>

// Tests for Gmsh reading
TEST_CASE("test if readGmsh creates any elements") {
    // TODO: this is a stub, please add to this
    int argc = 0;
    Kokkos::initialize(argc, nullptr);
    {
        std::string path = "../../dataset/bun_zipper_res4";
        GmshDataKokkos data = readGmsh(path);
    }
    Kokkos::finalize();
    REQUIRE(true);
}