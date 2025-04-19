#include <catch2/catch_test_macros.hpp>
#include "readGmsh.hpp"

TEST_CASE("Just checking") {
    std::string path = "../../dataset/bun_zipper_res4";
    GmshDataKokkos data = readGmsh(path);

    // Check the number of nodes
    REQUIRE(data.nodeTags.extent(0) == 453);
    REQUIRE(data.nodeCoords.extent(0) == 453); // Should match nodeTags
    REQUIRE(data.nodeCoords.extent(1) == 3);   // Each node has x, y, z

    REQUIRE(true);
}
