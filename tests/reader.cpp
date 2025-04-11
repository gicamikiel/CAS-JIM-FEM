#include <catch2/catch_test_macros.hpp>
#include <Kokkos_Core.hpp>
#include <readGmsh.hpp>

// Tests for Gmsh reading
TEST_CASE("test if readGmsh creates any elements") {
    std::string path = "../../dataset/bun_zipper_res4";
    GmshData data = readGmsh(path);
    REQUIRE(data.nodes.size() > 0);
    REQUIRE(data.triangles.size());
    REQUIRE(data.quadrangles.size() > 0);
    // note: will need much more than this
}