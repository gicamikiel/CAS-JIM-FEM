#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "readKappa.hpp"
#include <Kokkos_Core.hpp>

TEST_CASE("Kappa values should all be 0.1") {
    Kokkos::initialize();
    {
        std::string path = "../../dataset/kappa.txt";
        KappaDataKokkos data = readKappa(path);

        auto h_kappas = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), data.kappas);

        for (std::size_t i = 0; i < h_kappas.extent(0); ++i) {
            REQUIRE(h_kappas(i) == Catch::Approx(0.1).epsilon(1e-10));
        }
    }
    Kokkos::finalize();
}
